#define ESMF_ERR_ABORT(rc) \
if (rc /= ESMF_SUCCESS) write(0,*) 'rc=',rc,__FILE__,__LINE__; if(ESMF_LogFoundError(rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!-----------------------------------------------------------------------
!
  module module_fcst_grid_comp
!
!-----------------------------------------------------------------------
!***  Forecast gridded component.
!-----------------------------------------------------------------------
!***
!***  HISTORY
!***
!       Apr 2017:  J. Wang  - initial code for forecast grid component
!
!---------------------------------------------------------------------------------
!
  use mpi
  use esmf

  use time_manager_mod,   only: time_type, set_calendar_type, set_time,    &
                                set_date, month_name,                      &
                                operator(+), operator(-), operator (<),    &
                                operator (>), operator (/=), operator (/), &
                                operator (==), operator (*),               &
                                THIRTY_DAY_MONTHS, JULIAN, GREGORIAN,      &
                                NOLEAP, NO_CALENDAR,                       &
                                date_to_string, get_date, get_time

  use  atmos_model_mod,   only: atmos_model_init, atmos_model_end,         &
                                get_atmos_model_ungridded_dim,             &
                                update_atmos_model_dynamics,               &
                                update_atmos_radiation_physics,            &
                                update_atmos_model_state,                  &
                                atmos_data_type, atmos_model_restart,      &
                                atmos_model_exchange_phase_1,              &
                                atmos_model_exchange_phase_2,              &
                                addLsmask2grid

  use constants_mod,      only: constants_init
  use fms_mod,            only: error_mesg, fms_init, fms_end,             &
                                write_version_number, uppercase

  use mpp_mod,            only: mpp_init, mpp_pe, mpp_npes, mpp_root_pe,   &
                                mpp_error, FATAL, WARNING, NOTE
  use mpp_mod,            only: mpp_clock_id, mpp_clock_begin

  use mpp_io_mod,         only: mpp_open, mpp_close, MPP_DELETE

  use mpp_domains_mod,    only: mpp_get_compute_domains, domain2D
  use sat_vapor_pres_mod, only: sat_vapor_pres_init

  use diag_manager_mod,   only: diag_manager_init, diag_manager_end,       &
                                diag_manager_set_time_end

  use data_override_mod,  only: data_override_init
  use fv_nggps_diags_mod, only: fv_dyn_bundle_setup
  use fv3gfs_io_mod,      only: fv_phys_bundle_setup

  use fms_io_mod,         only: field_exist, read_data

  use atmosphere_mod,     only: atmosphere_control_data

  use module_fv3_io_def,  only: num_pes_fcst, num_files, filename_base,    &
                                nbdlphys, iau_offset
  use module_fv3_config,  only: dt_atmos, fcst_mpi_comm, fcst_ntasks,      &
                                quilting, calendar,                        &
                                cplprint_flag, restart_endfcst

  use get_stochy_pattern_mod, only: write_stoch_restart_atm
!
!-----------------------------------------------------------------------
!
  implicit none
!
!-----------------------------------------------------------------------
!
  private
!
!---- model defined-types ----

  type(atmos_data_type), save :: Atmos

  type(ESMF_Grid)             :: fcstGrid
  integer                     :: num_atmos_calls, intrm_rst

!----- coupled model data -----

  integer :: calendar_type = -99
  integer :: date_init(6)
  integer :: numLevels     = 0
  integer :: numSoilLayers = 0
  integer :: numTracers    = 0

  integer :: frestart(999)
!
!-----------------------------------------------------------------------
!
  public SetServices, fcstGrid
  public numLevels, numSoilLayers, numTracers
!
  contains
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
  subroutine SetServices(fcst_comp, rc)
!
    type(ESMF_GridComp)  :: fcst_comp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    call ESMF_GridCompSetEntryPoint(fcst_comp, ESMF_METHOD_INITIALIZE, &
                                    userRoutine=fcst_initialize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
    call ESMF_GridCompSetEntryPoint(fcst_comp, ESMF_METHOD_RUN, &
                                    userRoutine=fcst_run_phase_1, phase=1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompSetEntryPoint(fcst_comp, ESMF_METHOD_RUN, &
                                    userRoutine=fcst_run_phase_2, phase=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
    call ESMF_GridCompSetEntryPoint(fcst_comp, ESMF_METHOD_FINALIZE, &
                                    userRoutine=fcst_finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine SetServices
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
  subroutine fcst_initialize(fcst_comp, importState, exportState, clock, rc)
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE FORECAST GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
    type(esmf_GridComp)                    :: fcst_comp
    type(ESMF_State)                       :: importState, exportState
    type(esmf_Clock)                       :: clock
    integer,intent(out)                    :: rc
!
!***  local variables
!
    type(ESMF_VM)                          :: vm
    integer                                :: tl, i, j
    integer,dimension(2,6)                 :: decomptile                  !define delayout for the 6 cubed-sphere tiles
    integer,dimension(2)                   :: regdecomp                   !define delayout for the nest grid
    type(ESMF_FieldBundle)                 :: fieldbundle
!
    type(ESMF_Time)                        :: CurrTime, StartTime, StopTime
    type(ESMF_TimeInterval)                :: RunDuration
    type(ESMF_Config)                      :: cf

    integer                                :: Run_length
    integer,dimension(6)                   :: date, date_end
!
    character(len=9) :: month
    integer :: initClock, unit, total_inttime
    integer :: mype
    character(4) dateSY
    character(2) dateSM,dateSD,dateSH,dateSN,dateSS
    character(len=esmf_maxstr) name_FB, name_FB1
    character(len=80) :: dateS

    character(256)                         :: gridfile
    type(ESMF_FieldBundle),dimension(:), allocatable  :: fieldbundlephys

    real(kind=8) :: mpi_wtime, timeis

    type(ESMF_DELayout) :: delayout
    type(ESMF_DistGrid) :: distgrid
    real(ESMF_KIND_R8),dimension(:,:), pointer :: glatPtr, glonPtr
    real(ESMF_KIND_R8),parameter :: dtor = 180.0_ESMF_KIND_R8 / 3.1415926535897931_ESMF_KIND_R8
    integer :: jsc, jec, isc, iec, nlev
    type(domain2D)  :: domain
    integer :: n, fcstNpes, tmpvar
    logical :: freq_restart, fexist
    integer, allocatable, dimension(:) :: isl, iel, jsl, jel
    integer, allocatable, dimension(:,:,:) :: deBlockList
    integer :: tlb(2), tub(2)

    type(ESMF_Decomp_Flag)  :: decompflagPTile(2,6)

    integer               :: TileLayout(2)
    integer               :: nestRootPet, npes(1), peListSize(1)
    integer, allocatable  :: petMap(:)

    integer                       :: num_restart_interval, restart_starttime
    real,dimension(:),allocatable :: restart_interval
    type(time_type)               :: Time_init, Time, Time_step, Time_end, &
                                     Time_restart, Time_step_restart
    type(time_type)               :: iautime
    integer                       :: io_unit, calendar_type_res, date_res(6), date_init_res(6)

!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
    timeis = mpi_wtime()
    rc     = ESMF_SUCCESS
!
    call ESMF_VMGetCurrent(vm=vm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_VMGet(vm=vm, localPet=mype, mpiCommunicator=fcst_mpi_comm, &
                    petCount=fcst_ntasks, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (mype == 0) write(*,*)'in fcst comp init, fcst_ntasks=',fcst_ntasks

    CF = ESMF_ConfigCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_ConfigLoadFile(config=CF ,filename='model_configure' ,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    num_restart_interval = ESMF_ConfigGetLen(config=CF, label ='restart_interval:',rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (mype == 0) print *,'af nems config,num_restart_interval=',num_restart_interval
    if (num_restart_interval<=0) num_restart_interval = 1
    allocate(restart_interval(num_restart_interval))
    restart_interval = 0
    call ESMF_ConfigGetAttribute(CF,valueList=restart_interval,label='restart_interval:', &
                                 count=num_restart_interval, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (mype == 0) print *,'af nems config,restart_interval=',restart_interval
!
    call fms_init(fcst_mpi_comm)
    call mpp_init()
    initClock = mpp_clock_id( 'Initialization' )
    call mpp_clock_begin (initClock) !nesting problem

    call constants_init
    call sat_vapor_pres_init

    select case( uppercase(trim(calendar)) )
    case( 'JULIAN' )
        calendar_type = JULIAN
    case( 'GREGORIAN' )
        calendar_type = GREGORIAN
    case( 'NOLEAP' )
        calendar_type = NOLEAP
    case( 'THIRTY_DAY' )
        calendar_type = THIRTY_DAY_MONTHS
    case( 'NO_CALENDAR' )
        calendar_type = NO_CALENDAR
    case default
        call mpp_error ( FATAL, 'fcst_initialize: calendar must be one of '// &
                                'JULIAN|GREGORIAN|NOLEAP|THIRTY_DAY|NO_CALENDAR.' )
    end select

    call set_calendar_type (calendar_type)
!
!-----------------------------------------------------------------------
!***  set atmos time
!-----------------------------------------------------------------------
!
    call ESMF_ClockGet(clock, CurrTime=CurrTime, StartTime=StartTime, &
                       StopTime=StopTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    date_init = 0
    call ESMF_TimeGet (StartTime,                      &
                       YY=date_init(1), MM=date_init(2), DD=date_init(3), &
                       H=date_init(4),  M =date_init(5), S =date_init(6), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    Time_init  = set_date (date_init(1), date_init(2), date_init(3), &
                           date_init(4), date_init(5), date_init(6))
    if (mype == 0) write(*,'(A,6I5)') 'StartTime=',date_init

    date=0
    call ESMF_TimeGet (CurrTime,                           &
                       YY=date(1), MM=date(2), DD=date(3), &
                       H=date(4),  M =date(5), S =date(6), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    Time = set_date (date(1), date(2), date(3),  &
                     date(4), date(5), date(6))
    if (mype == 0) write(*,'(A,6I5)') 'CurrTime =',date

    date_end=0
    call ESMF_TimeGet (StopTime,                                       &
                       YY=date_end(1), MM=date_end(2), DD=date_end(3), &
                       H=date_end(4),  M =date_end(5), S =date_end(6), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    Time_end   = set_date (date_end(1), date_end(2), date_end(3),  &
                           date_end(4), date_end(5), date_end(6))
    if (mype == 0) write(*,'(A,6I5)') 'StopTime =',date_end

!------------------------------------------------------------------------
!   If this is a restarted run ('INPUT/coupler.res' file exists),
!   compare date and date_init to the values in 'coupler.res'

    if (mype == 0) then
      inquire(FILE='INPUT/coupler.res', EXIST=fexist)
      if (fexist) then  ! file exists, this is a restart run

        call ESMF_UtilIOUnitGet(unit=io_unit, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        open(unit=io_unit, file='INPUT/coupler.res', status='old', action='read', err=998)
        read (io_unit,*,err=999) calendar_type_res
        read (io_unit,*) date_init_res
        read (io_unit,*) date_res
        close(io_unit)

        if(date_res(1) == 0 .and. date_init_res(1) /= 0) date_res = date_init_res

        if(mype == 0) write(*,'(A,6(I4))') 'INPUT/coupler.res: date_init=',date_init_res
        if(mype == 0) write(*,'(A,6(I4))') 'INPUT/coupler.res: date     =',date_res

        if (calendar_type /= calendar_type_res) then
          write(0,'(A)')      'fcst_initialize ERROR: calendar_type /= calendar_type_res'
          write(0,'(A,6(I4))')'                       calendar_type     = ', calendar_type
          write(0,'(A,6(I4))')'                       calendar_type_res = ', calendar_type_res
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        endif

        if (.not. ALL(date_init.EQ.date_init_res)) then
          write(0,'(A)')      'fcst_initialize ERROR: date_init /= date_init_res'
          write(0,'(A,6(I4))')'                       date_init     = ', date_init
          write(0,'(A,6(I4))')'                       date_init_res = ', date_init_res
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        endif

        if (.not. ALL(date.EQ.date_res)) then
          write(0,'(A)')      'fcst_initialize ERROR: date /= date_res'
          write(0,'(A,6(I4))')'                       date     = ', date
          write(0,'(A,6(I4))')'                       date_res = ', date_res
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        endif

  999 continue
  998 continue

      endif ! fexist
    endif ! mype == 0

    RunDuration = StopTime - CurrTime

    CALL ESMF_TimeIntervalGet(RunDuration, S=Run_length, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
    call diag_manager_init (TIME_INIT=date)
    call diag_manager_set_time_end(Time_end)
!
    Time_step = set_time (dt_atmos,0)
    num_atmos_calls = Run_length / dt_atmos
    if (mype == 0) write(*,*)'num_atmos_calls=',num_atmos_calls,'time_init=', &
                    date_init,'time=',date,'time_end=',date_end,'dt_atmos=',dt_atmos, &
                    'Run_length=',Run_length

! set up forecast time array that controls when to write out restart files
    frestart = 0
    call get_time(Time_end - Time_init, total_inttime)
! set iau offset time
    Atmos%iau_offset    = iau_offset
    if(iau_offset > 0 ) then
      iautime =  set_time(iau_offset * 3600, 0)
    endif
! if the second item is -1, the first number is frequency
    freq_restart = .false.
    if(num_restart_interval == 2) then
      if(restart_interval(2)== -1) freq_restart = .true.
    endif
    if(freq_restart) then
      if(restart_interval(1) >= 0) then
        tmpvar = restart_interval(1) * 3600
        Time_step_restart = set_time (tmpvar, 0)
        if(iau_offset > 0 ) then
          Time_restart = Time_init + iautime + Time_step_restart
          frestart(1) = tmpvar + iau_offset *3600
        else
          Time_restart = Time_init + Time_step_restart
          frestart(1) = tmpvar
        endif
        if(restart_interval(1) > 0) then
          i = 2
          do while ( Time_restart < Time_end )
            frestart(i) = frestart(i-1) + tmpvar
            Time_restart = Time_restart + Time_step_restart
             i = i + 1
          enddo
        endif
      endif
! otherwise it is an array with forecast time at which the restart files will be written out
    else if(num_restart_interval >= 1) then
      if(num_restart_interval == 1 .and. restart_interval(1) == 0 ) then
        frestart(1) = total_inttime
      else
        if(iau_offset > 0 ) then
          restart_starttime = iau_offset *3600
        else
          restart_starttime = 0
        endif
        do i=1,num_restart_interval
          frestart(i) = restart_interval(i) * 3600. + restart_starttime
        enddo
      endif
    endif
! if to write out restart at the end of forecast
    restart_endfcst = .false.
    if ( ANY(frestart(:) == total_inttime) ) restart_endfcst = .true.
    if (mype == 0) print *,'frestart=',frestart(1:10)/3600, 'restart_endfcst=',restart_endfcst, &
      'total_inttime=',total_inttime
! if there is restart writing during integration
    intrm_rst         = 0
    if (frestart(1)>0) intrm_rst = 1
!
!----- write time stamps (for start time and end time) ------

     call mpp_open( unit, 'time_stamp.out', nohdrs=.TRUE. )
     month = month_name(date(2))
     if ( mpp_pe() == mpp_root_pe() ) write (unit,20) date, month(1:3)
     month = month_name(date_end(2))
     if ( mpp_pe() == mpp_root_pe() ) write (unit,20) date_end, month(1:3)
     call mpp_close (unit)
 20  format (6i4,2x,a3)
!
!------ initialize component models ------

     call  atmos_model_init (Atmos, Time_init, Time, Time_step)
!
     inquire(FILE='data_table', EXIST=fexist)
     if (fexist) then
       call data_override_init()
     endif
!-----------------------------------------------------------------------
!---- open and close dummy file in restart dir to check if dir exists --

      if (mpp_pe() == 0 ) then
         call mpp_open( unit, 'RESTART/file' )
         call mpp_close(unit, MPP_DELETE)
      endif
!
!
!-----------------------------------------------------------------------
!*** create grid for output fields
!*** first try: Create cubed sphere grid from file
!-----------------------------------------------------------------------
!
      call mpp_error(NOTE, 'before create fcst grid')

      gridfile = "grid_spec.nc" ! default

      if (field_exist("INPUT/grid_spec.nc", "atm_mosaic_file")) then
        call read_data("INPUT/grid_spec.nc", "atm_mosaic_file", gridfile)
      endif

      if (mpp_pe() == mpp_root_pe()) &
      write(*, *) 'create fcst grid: mype,regional,nested=',mype,Atmos%regional,Atmos%nested

      ! regional-only without nests
      if( Atmos%regional .and. .not. Atmos%nested ) then

        call atmosphere_control_data (isc, iec, jsc, jec, nlev)

        domain   = Atmos%domain
        fcstNpes = Atmos%layout(1)*Atmos%layout(2)
        allocate(isl(fcstNpes), iel(fcstNpes), jsl(fcstNpes), jel(fcstNpes))
        allocate(deBlockList(2,2,fcstNpes))
        call mpp_get_compute_domains(domain,xbegin=isl,xend=iel,ybegin=jsl,yend=jel)
        do n=1,fcstNpes
           deBlockList(:,1,n) = (/ isl(n),iel(n) /)
           deBlockList(:,2,n) = (/ jsl(n),jel(n) /)
        end do
        delayout = ESMF_DELayoutCreate(petMap=(/(i,i=0,fcstNpes-1)/), rc=rc); ESMF_ERR_ABORT(rc)
        distgrid = ESMF_DistGridCreate(minIndex=(/1,1/), &
                                         maxIndex=(/Atmos%mlon,Atmos%mlat/), &
                                         delayout=delayout, &
                                         deBlockList=deBlockList, rc=rc); ESMF_ERR_ABORT(rc)

        fcstGrid = ESMF_GridCreateNoPeriDim(regDecomp=(/Atmos%layout(1),Atmos%layout(2)/), &
                                              minIndex=(/1,1/), &
                                              maxIndex=(/Atmos%mlon,Atmos%mlat/), &
                                              gridAlign=(/-1,-1/), &
                                              decompflag=(/ESMF_DECOMP_SYMMEDGEMAX,ESMF_DECOMP_SYMMEDGEMAX/), &
                                              name="fcst_grid", &
                                              indexflag=ESMF_INDEX_DELOCAL, &
                                              rc=rc); ESMF_ERR_ABORT(rc)

        ! add and define "center" coordinate values
        call ESMF_GridAddCoord(fcstGrid, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc); ESMF_ERR_ABORT(rc)
        call ESMF_GridGetCoord(fcstGrid, coordDim=1, staggerLoc=ESMF_STAGGERLOC_CENTER, &
                                 farrayPtr=glonPtr, rc=rc); ESMF_ERR_ABORT(rc)
        call ESMF_GridGetCoord(fcstGrid, coordDim=2, staggerLoc=ESMF_STAGGERLOC_CENTER, &
                                 farrayPtr=glatPtr, rc=rc); ESMF_ERR_ABORT(rc)

        do j = jsc, jec
          do i = isc, iec
            glonPtr(i-isc+1,j-jsc+1) = Atmos%lon(i-isc+1,j-jsc+1) * dtor
            glatPtr(i-isc+1,j-jsc+1) = Atmos%lat(i-isc+1,j-jsc+1) * dtor
          enddo
        enddo

        ! add and define "corner" coordinate values
        call ESMF_GridAddCoord(fcstGrid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
                               rc=rc); ESMF_ERR_ABORT(rc)
        call ESMF_GridGetCoord(fcstGrid, coordDim=1, staggerLoc=ESMF_STAGGERLOC_CORNER, &
                               totalLBound=tlb, totalUBound=tub, &
                               farrayPtr=glonPtr, rc=rc); ESMF_ERR_ABORT(rc)
        glonPtr(tlb(1):tub(1),tlb(2):tub(2)) = &
           Atmos%lon_bnd(tlb(1):tub(1),tlb(2):tub(2)) * dtor
        call ESMF_GridGetCoord(fcstGrid, coordDim=2, staggerLoc=ESMF_STAGGERLOC_CORNER, &
                               totalLBound=tlb, totalUBound=tub, &
                               farrayPtr=glatPtr, rc=rc); ESMF_ERR_ABORT(rc)
        glatPtr(tlb(1):tub(1),tlb(2):tub(2)) = &
          Atmos%lat_bnd(tlb(1):tub(1),tlb(2):tub(2)) * dtor

        call mpp_error(NOTE, 'after create fcst grid for regional-only')

      else ! not regional only

        if (.not. Atmos%regional .and. .not. Atmos%nested ) then  !! global only

          do tl=1,6
              decomptile(1,tl) = Atmos%layout(1)
              decomptile(2,tl) = Atmos%layout(2)
              decompflagPTile(:,tl) = (/ESMF_DECOMP_SYMMEDGEMAX,ESMF_DECOMP_SYMMEDGEMAX/)
          enddo
          fcstGrid = ESMF_GridCreateMosaic(filename="INPUT/"//trim(gridfile),                                   &
                                             regDecompPTile=decomptile,tileFilePath="INPUT/",                   &
                                             decompflagPTile=decompflagPTile,                                   &
                                             staggerlocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
                                             name='fcst_grid', rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call mpp_error(NOTE, 'after create fcst grid for global-only with INPUT/'//trim(gridfile))

        else !! global-nesting or regional-nesting

          if (mype == 0) TileLayout = Atmos%layout
          call ESMF_VMBroadcast(vm, bcstData=TileLayout, count=2, &
                                  rootPet=0, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if (mype == 0) npes(1) = mpp_npes()
          call ESMF_VMBroadcast(vm, bcstData=npes, count=1, &
                                  rootPet=0, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if ( npes(1) == TileLayout(1) * TileLayout(2) * 6 ) then
            ! global-nesting
            nestRootPet = npes(1)
            gridfile="grid.nest02.tile7.nc"
          else if ( npes(1) == TileLayout(1) * TileLayout(2) ) then
            ! regional-nesting
            nestRootPet = npes(1)
            gridfile="grid.nest02.tile2.nc"
          else
            call mpp_error(FATAL, 'Inconsistent nestRootPet and Atmos%layout')
          endif

          if (mype == nestRootPet) then
            if (nestRootPet /= Atmos%pelist(1)) then
              write(0,*)'error in fcst_initialize: nestRootPet /= Atmos%pelist(1)'
              write(0,*)'error in fcst_initialize: nestRootPet = ',nestRootPet
              write(0,*)'error in fcst_initialize: Atmos%pelist(1) = ',Atmos%pelist(1)
              ESMF_ERR_ABORT(100)
            endif
          endif

          ! nest rootPet shares peList with others
          if (mype == nestRootPet) peListSize(1) = size(Atmos%pelist)
          call ESMF_VMBroadcast(vm, bcstData=peListSize, count=1, rootPet=nestRootPet, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          ! nest rootPet shares layout with others
          if (mype == nestRootPet) regDecomp = Atmos%layout
          call ESMF_VMBroadcast(vm, bcstData=regDecomp, count=2, rootPet=nestRootPet, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          ! prepare petMap variable
          allocate(petMap(peListSize(1)))
          if (mype == nestRootPet) petMap = Atmos%pelist
          ! do the actual broadcast of the petMap
          call ESMF_VMBroadcast(vm, bcstData=petMap, count=peListSize(1), rootPet=nestRootPet, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          ! create the DELayout that maps DEs to the PETs in the petMap
          delayout = ESMF_DELayoutCreate(petMap=petMap, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          ! create the nest Grid by reading it from file but use DELayout
          fcstGrid = ESMF_GridCreate(filename="INPUT/"//trim(gridfile),                                   &
                                       fileformat=ESMF_FILEFORMAT_GRIDSPEC, regDecomp=regDecomp,          &
                                       decompflag=(/ESMF_DECOMP_SYMMEDGEMAX,ESMF_DECOMP_SYMMEDGEMAX/),    &
                                       delayout=delayout, isSphere=.false., indexflag=ESMF_INDEX_DELOCAL, &
                                       rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call mpp_error(NOTE, 'after create fcst grid with INPUT/'//trim(gridfile))

        endif

      endif
!
      !! FIXME
      if ( .not. Atmos%nested ) then  !! global only
        call addLsmask2grid(fcstGrid, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!       print *,'call addLsmask2grid after fcstGrid, rc=',rc
      endif

!test to write out vtk file:
!     if( cplprint_flag ) then
!       call ESMF_GridWriteVTK(fcstGrid, staggerloc=ESMF_STAGGERLOC_CENTER,  &
!                              filename='fv3cap_fv3Grid', rc=rc)
!       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!     endif
!
! Write grid to netcdf file
      if( cplprint_flag ) then
        call wrt_fcst_grid(fcstGrid, "diagnostic_FV3_fcstGrid.nc", &
                                 regridArea=.TRUE., rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      endif

! Add gridfile Attribute to the exportState
      call ESMF_AttributeAdd(exportState, convention="NetCDF", purpose="FV3", &
                               attrList=(/"gridfile"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="gridfile", value=trim(gridfile), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! Add dimension Attributes to Grid
      call ESMF_AttributeAdd(fcstGrid, convention="NetCDF", purpose="FV3",  &
                               attrList=(/"ESMF:gridded_dim_labels"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(fcstGrid, convention="NetCDF", purpose="FV3", &
                               name="ESMF:gridded_dim_labels", valueList=(/"grid_xt", "grid_yt"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
! Add time Attribute to the exportState
      call ESMF_AttributeAdd(exportState, convention="NetCDF", purpose="FV3", &
        attrList=(/ "time               ", &
                      "time:long_name     ", &
                      "time:units         ", &
                      "time:cartesian_axis", &
                      "time:calendar_type ", &
                      "time:calendar      " /), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="time", value=real(0,ESMF_KIND_R8), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      write(dateSY,'(I4.4)')date_init(1)
      write(dateSM,'(I2.2)')date_init(2)
      write(dateSD,'(I2.2)')date_init(3)
      write(dateSH,'(I2.2)')date_init(4)
      write(dateSN,'(I2.2)')date_init(5)
      write(dateSS,'(I2.2)')date_init(6)

      dateS="hours since "//dateSY//'-'//dateSM//'-'//dateSD//' '//dateSH//':'//    &
            dateSN//":"//dateSS
      if (mype == 0) write(*,*)'dateS=',trim(dateS),'date_init=',date_init

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="time:units", value=trim(dateS), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="time:long_name", value="time", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="time:cartesian_axis", value="T", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="time:calendar_type", value=uppercase(trim(calendar)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="time:calendar", value=uppercase(trim(calendar)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
! Create FieldBundle for Fields that need to be regridded bilinear
      if( quilting ) then

        do i=1,num_files
!
         name_FB = filename_base(i)
!
         if( i==1 ) then
! for dyn
           name_FB1 = trim(name_FB)//'_bilinear'
           fieldbundle = ESMF_FieldBundleCreate(name=trim(name_FB1),rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
           if (mype == 0) write(*,*)'af create fcst fieldbundle, name=',trim(name_FB),'rc=',rc

           call fv_dyn_bundle_setup(Atmos%axes,          &
                                    fieldbundle, fcstGrid, quilting, rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

           ! Add the field to the importState so parent can connect to it
           call ESMF_StateAdd(exportState, (/fieldbundle/), rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

         else if( i==2 ) then
! for phys
           nbdlphys = 2
           allocate(fieldbundlephys(nbdlphys))
           do j=1, nbdlphys
             if( j==1 ) then
               name_FB1 = trim(name_FB)//'_nearest_stod'
             else
               name_FB1 = trim(name_FB)//'_bilinear'
             endif
             fieldbundlephys(j) = ESMF_FieldBundleCreate(name=trim(name_FB1),rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
             if (mype == 0) write(*,*)'af create fcst fieldbundle, name=',trim(name_FB1),'rc=',rc
           enddo
!
           call fv_phys_bundle_setup(Atmos%diag, Atmos%axes, &
                                     fieldbundlephys, fcstGrid, quilting, nbdlphys)
!
           ! Add the field to the importState so parent can connect to it
           do j=1,nbdlphys
             call ESMF_StateAdd(exportState, (/fieldbundlephys(j)/), rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
           enddo

         else

           write(0,*)' unknown name_FB ', trim(name_FB)
           ESMF_ERR_ABORT(101)

         endif
!
        enddo

!end qulting
      endif

      call get_atmos_model_ungridded_dim(nlev=numLevels,         &
                                         nsoillev=numSoilLayers, &
                                         ntracers=numTracers)

      if (mype == 0) write(*,*)'fcst_initialize total time: ', mpi_wtime() - timeis
!
!-----------------------------------------------------------------------
!
   end subroutine fcst_initialize
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
   subroutine fcst_run_phase_1(fcst_comp, importState, exportState,clock,rc)
!
!-----------------------------------------------------------------------
!***  the run step for the fcst gridded component.
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp)        :: fcst_comp
      type(ESMF_State)           :: importState, exportState
      type(ESMF_Clock)           :: clock
      integer,intent(out)        :: rc
!
!***  local variables
!
      integer                    :: mype, na
      integer(kind=ESMF_KIND_I8) :: ntimestep_esmf
      real(kind=8)               :: mpi_wtime, tbeg1
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      tbeg1 = mpi_wtime()
      rc    = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
      call ESMF_GridCompGet(fcst_comp, localpet=mype, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
      call ESMF_ClockGet(clock, advanceCount=NTIMESTEP_ESMF, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      na = NTIMESTEP_ESMF
!
!-----------------------------------------------------------------------
! *** call fcst integration subroutines

      call update_atmos_model_dynamics (Atmos)

      call update_atmos_radiation_physics (Atmos)

      call atmos_model_exchange_phase_1 (Atmos, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if (mype == 0) write(*,*)"PASS: fcstRUN phase 1, na = ",na, ' time is ', mpi_wtime()-tbeg1
!
!-----------------------------------------------------------------------
!
   end subroutine fcst_run_phase_1
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
   subroutine fcst_run_phase_2(fcst_comp, importState, exportState,clock,rc)
!
!-----------------------------------------------------------------------
!***  the run step for the fcst gridded component.
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp)        :: fcst_comp
      type(ESMF_State)           :: importState, exportState
      type(ESMF_Clock)           :: clock
      integer,intent(out)        :: rc
!
!***  local variables
!
      integer                    :: mype, na, date(6), seconds
      integer(kind=ESMF_KIND_I8) :: ntimestep_esmf
      character(len=64)          :: timestamp
      integer                    :: unit
      real(kind=8)               :: mpi_wtime, tbeg1
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      tbeg1 = mpi_wtime()
      rc    = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
      call ESMF_GridCompGet(fcst_comp, localpet=mype, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_ClockGet(clock, advanceCount=NTIMESTEP_ESMF, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      na = NTIMESTEP_ESMF
!
!-----------------------------------------------------------------------
! *** call fcst integration subroutines

      call atmos_model_exchange_phase_2 (Atmos, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call update_atmos_model_state (Atmos, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      !--- intermediate restart
      if (intrm_rst>0) then
        if (na /= num_atmos_calls-1) then
          call get_time(Atmos%Time - Atmos%Time_init, seconds)
          if (ANY(frestart(:) == seconds)) then
            if (mype == 0) write(*,*)'write out restart at na=',na,' seconds=',seconds,  &
                                     'integration lenght=',na*dt_atmos/3600.

            timestamp = date_to_string (Atmos%Time)
            call atmos_model_restart(Atmos, timestamp)
            call write_stoch_restart_atm('RESTART/'//trim(timestamp)//'.atm_stoch.res.nc')

            !----- write restart file ------
            if (mpp_pe() == mpp_root_pe())then
                call get_date (Atmos%Time, date(1), date(2), date(3),  &
                                                       date(4), date(5), date(6))
                call mpp_open( unit, 'RESTART/'//trim(timestamp)//'.coupler.res', nohdrs=.TRUE. )
                write( unit, '(i6,8x,a)' )calendar_type, &
                     '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'

                write( unit, '(6i6,8x,a)' )date_init, &
                     'Model start time:   year, month, day, hour, minute, second'
                write( unit, '(6i6,8x,a)' )date, &
                     'Current model time: year, month, day, hour, minute, second'
                call mpp_close(unit)
            endif
          endif
        endif
      endif

      if (mype == 0) write(*,*)"PASS: fcstRUN phase 2, na = ",na, ' time is ', mpi_wtime()-tbeg1
!
!-----------------------------------------------------------------------
!
   end subroutine fcst_run_phase_2
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
   subroutine fcst_finalize(fcst_comp, importState, exportState,clock,rc)
!
!-----------------------------------------------------------------------
!***  finalize the forecast grid component.
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp)        :: fcst_comp
      type(ESMF_State)           :: importState, exportState
      type(ESMF_Clock)           :: clock
      integer,intent(out)        :: rc
!
!***  local variables
!
      integer                    :: mype
      integer                    :: unit
      integer,dimension(6)       :: date
      real(kind=8)               :: mpi_wtime, tbeg1
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      tbeg1 = mpi_wtime()
      rc    = ESMF_SUCCESS

      call ESMF_GridCompGet(fcst_comp, localpet=mype, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call atmos_model_end (Atmos)

!*** write restart file
      if( restart_endfcst ) then
        call get_date (Atmos%Time, date(1), date(2), date(3),  &
                               date(4), date(5), date(6))
        call mpp_open( unit, 'RESTART/coupler.res', nohdrs=.TRUE. )
        if (mpp_pe() == mpp_root_pe())then
          write( unit, '(i6,8x,a)' )calendar_type, &
              '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'

          write( unit, '(6i6,8x,a)' )date_init, &
              'Model start time:   year, month, day, hour, minute, second'
          write( unit, '(6i6,8x,a)' )date, &
              'Current model time: year, month, day, hour, minute, second'
        endif
        call mpp_close(unit)
      endif

      call diag_manager_end (Atmos%Time)

      call fms_end

      if (mype == 0) write(*,*)'fcst_finalize total time: ', mpi_wtime() - tbeg1
!
!-----------------------------------------------------------------------
!
  end subroutine fcst_finalize
!
!#######################################################################
!-- write forecast grid to NetCDF file for diagnostics
!
  subroutine wrt_fcst_grid(grid, fileName, relaxedflag, regridArea, rc)
    type(ESMF_Grid), intent(in)                      :: grid
    character(len=*), intent(in), optional           :: fileName
    logical, intent(in), optional                    :: relaxedflag
    logical, intent(in), optional                    :: regridArea
    integer, intent(out)                             :: rc
!
!***  local variables
!
    logical                     :: ioCapable
    logical                     :: doItFlag
    character(len=64)           :: lfileName
    character(len=64)           :: gridName
    type(ESMF_Array)            :: array
    type(ESMF_ArrayBundle)      :: arraybundle
    logical                     :: isPresent
    logical                     :: hasCorners
    logical                     :: lRegridArea
    type(ESMF_Field)            :: areaField
    type(ESMF_FieldStatus_Flag) :: areaFieldStatus

    ioCapable = (ESMF_IO_PIO_PRESENT .and. &
      (ESMF_IO_NETCDF_PRESENT .or. ESMF_IO_PNETCDF_PRESENT))
    doItFlag = .true.
    if (present(relaxedFlag)) then
      doItFlag = .not.relaxedflag .or. (relaxedflag.and.ioCapable)
    endif

    if (doItFlag) then
      ! Process optional arguments
      if (present(fileName)) then
        lfileName = trim(fileName)
      else
        call ESMF_GridGet(grid, name=gridName, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        lfileName = trim(gridName)//".nc"
      endif
      if (present(regridArea)) then
        lRegridArea = regridArea
      else
        lRegridArea = .FALSE.
      endif

      ! Create bundle for storing output
      arraybundle = ESMF_ArrayBundleCreate(rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      ! -- Centers --
      call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
        isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (isPresent) then
        call ESMF_GridGetCoord(grid, coordDim=1, &
          staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArraySet(array, name="lon_center", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_GridGetCoord(grid, coordDim=2, &
          staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArraySet(array, name="lat_center", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      endif

      ! -- Corners --
      call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
        isPresent=hasCorners, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (hasCorners) then
        call ESMF_GridGetCoord(grid, coordDim=1, &
          staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
        if (.not. ESMF_LogFoundError(rc, line=__LINE__, file=__FILE__)) then
          call ESMF_ArraySet(array, name="lon_corner", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        endif
        call ESMF_GridGetCoord(grid, coordDim=2, &
          staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
        if (.not. ESMF_LogFoundError(rc, line=__LINE__, file=__FILE__)) then
          call ESMF_ArraySet(array, name="lat_corner", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        endif
        if (lRegridArea) then
          areaField = ESMF_FieldCreate(grid=grid, &
            typekind=ESMF_TYPEKIND_R8, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
          call ESMF_FieldRegridGetArea(areaField, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          call ESMF_FieldGet(areaField, array=array, rc=rc)
          if (.not. ESMF_LogFoundError(rc, line=__LINE__, file=__FILE__)) then
            call ESMF_ArraySet(array, name="regrid_area", rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=__FILE__)) return
            call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=__FILE__)) return
          endif
        endif
      endif

      ! -- Mask --
      call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_MASK, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (isPresent) then
        call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
          itemflag=ESMF_GRIDITEM_MASK, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArraySet(array, name="mask", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      endif

      ! -- Area --
      call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_AREA, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (isPresent) then
        call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
          itemflag=ESMF_GRIDITEM_AREA, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArraySet(array, name="area", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      endif

      ! Write array bundle to grid file
      call ESMF_ArrayBundleWrite(arraybundle, fileName=trim(lfileName), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      ! Clean-up
      if (lRegridArea) then
        call ESMF_FieldGet(areaField, status=areaFieldStatus, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        if (areaFieldStatus.eq.ESMF_FIELDSTATUS_COMPLETE) then
          call ESMF_FieldDestroy(areaField, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        endif
      endif
      call ESMF_ArrayBundleDestroy(arraybundle,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif
  end subroutine wrt_fcst_grid
!
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
end module  module_fcst_grid_comp
!
!----------------------------------------------------------------------------
