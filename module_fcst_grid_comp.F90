#ifdef __PGI
#define ESMF_ERR_ABORT(rc) \
if (rc /= ESMF_SUCCESS) write(0,*) 'rc=',rc,__FILE__,__LINE__; call ESMF_Finalize(endflag=ESMF_END_ABORT)
#else
#define ESMF_ERR_ABORT(rc) \
if (rc /= ESMF_SUCCESS) write(0,*) 'rc=',rc,__FILE__,__LINE__; if(ESMF_LogFoundError(rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
#endif

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

  use diag_manager_mod,   only: diag_manager_init, diag_manager_end, &
                                diag_manager_set_time_end

  use data_override_mod,  only: data_override_init
  use fv_nggps_diags_mod, only: fv_dyn_bundle_setup
  use fv3gfs_io_mod,      only: fv_phys_bundle_setup

  use fms_io_mod,         only: field_exist, read_data

  use atmosphere_mod,     only: atmosphere_control_data
!
  use module_fv3_io_def, only:  num_pes_fcst, num_files, filename_base, nbdlphys, &
                                iau_offset
  use module_fv3_config, only:  dt_atmos, calendar, fcst_mpi_comm, fcst_ntasks,   &
                                quilting, calendar_type,                          &
                                cplprint_flag, force_date_from_configure,         &
                                restart_endfcst
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

  type atmos_internalstate_type
    type(atmos_data_type)  :: Atm
    type(time_type)        :: Time_atmos, Time_init, Time_end,  &
                              Time_step_atmos, Time_step_ocean, &
                              Time_restart, Time_step_restart,  &
                              Time_atstart
    integer :: num_atmos_calls, ret, intrm_rst
  end type

  type atmos_internalstate_wrapper
    type(atmos_internalstate_type), pointer :: ptr
  end type

  type(atmos_internalstate_type),pointer,save :: atm_int_state
  type(atmos_internalstate_wrapper),save      :: wrap
  type(ESMF_VM),save                          :: VM
  type(ESMF_Grid)                             :: fcstGrid

!----- coupled model data -----

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
!***  INITIALIZE THE WRITE GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
    type(esmf_GridComp)                    :: fcst_comp
    type(ESMF_State)                       :: importState, exportState
    type(esmf_Clock)                       :: clock
    integer,intent(out)                    :: rc
!
!***  LOCAL VARIABLES
!
    integer                                :: tl, i, j
    integer,dimension(2,6)                 :: decomptile                  !define delayout for the 6 cubed-sphere tiles
    integer,dimension(2)                   :: regdecomp                   !define delayout for the nest grid
    type(ESMF_FieldBundle)                 :: fieldbundle
!
    type(ESMF_Time)                        :: CurrTime, StartTime, StopTime
    type(ESMF_TimeInterval)                :: RunDuration, TimeElapsed
    type(ESMF_Config)                      :: cf

    integer                                :: Run_length
    integer,dimension(6)                   :: date, date_end
!
    character(len=9) :: month
    integer :: initClock, unit, nfhour, total_inttime
    integer :: mype
    character(3) cfhour
    character(4) dateSY
    character(2) dateSM,dateSD,dateSH,dateSN,dateSS
    character(len=esmf_maxstr) name_FB, name_FB1
    character(len=80) :: dateS
    real,    allocatable, dimension(:,:) :: glon_bnd, glat_bnd

    character(256)                         :: gridfile
    type(ESMF_FieldBundle),dimension(:), allocatable  :: fieldbundlephys

    real(8) :: mpi_wtime, timeis

    type(ESMF_DELayout) :: delayout
    type(ESMF_DistGrid) :: distgrid
    real(ESMF_KIND_R8),dimension(:,:), pointer :: glatPtr, glonPtr
    real(ESMF_KIND_R8),parameter :: dtor = 180.0_ESMF_KIND_R8 / 3.1415926535897931_ESMF_KIND_R8
    integer :: jsc, jec, isc, iec, nlev
    type(domain2D)  :: domain
    type(time_type) :: iautime
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
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
    timeis = mpi_wtime()
    rc     = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE WRITE COMPONENT'S INTERNAL STATE.
!-----------------------------------------------------------------------
!
    allocate(atm_int_state,stat=rc)
!
!-----------------------------------------------------------------------
!***  ATTACH THE INTERNAL STATE TO THE WRITE COMPONENT.
!-----------------------------------------------------------------------
!
    wrap%ptr => atm_int_state
    call ESMF_GridCompSetInternalState(fcst_comp, wrap, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
    call ESMF_VMGetCurrent(vm=VM,rc=RC)
    call ESMF_VMGet(vm=VM, localPet=mype, mpiCommunicator=fcst_mpi_comm, &
                    petCount=fcst_ntasks, rc=rc)
    if (mype == 0) write(0,*)'in fcst comp init, fcst_ntasks=',fcst_ntasks
!
    CF = ESMF_ConfigCreate(rc=rc)
    call ESMF_ConfigLoadFile(config=CF ,filename='model_configure' ,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
    num_restart_interval = ESMF_ConfigGetLen(config=CF, label ='restart_interval:',rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if(mype == 0) print *,'af nems config,num_restart_interval=',num_restart_interval
    if (num_restart_interval<=0) num_restart_interval = 1
    allocate(restart_interval(num_restart_interval))
    restart_interval = 0
    call ESMF_ConfigGetAttribute(CF,valueList=restart_interval,label='restart_interval:', &
      count=num_restart_interval, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if(mype == 0) print *,'af nems config,restart_interval=',restart_interval

!
    call fms_init(fcst_mpi_comm)
    call mpp_init()
    initClock = mpp_clock_id( 'Initialization' )
    call mpp_clock_begin (initClock) !nesting problem

    call fms_init
    call constants_init
    call sat_vapor_pres_init
!
    if ( force_date_from_configure ) then

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

    endif
!
    call set_calendar_type (calendar_type         )
!
!-----------------------------------------------------------------------
!***  set atmos time
!-----------------------------------------------------------------------
!
    call ESMF_ClockGet(clock, CurrTime=CurrTime, StartTime=StartTime, &
                       StopTime=StopTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    RunDuration = StopTime - CurrTime

    date_init = 0
    call ESMF_TimeGet (StartTime,                      &
                       YY=date_init(1), MM=date_init(2), DD=date_init(3), &
                       H=date_init(4),  M =date_init(5), S =date_init(6), RC=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if ( date_init(1) == 0 ) date_init = date
    atm_int_state%Time_init  = set_date (date_init(1), date_init(2), date_init(3), &
                                         date_init(4), date_init(5), date_init(6))
    if(mype==0) write(*,'(A,6I5)') 'StartTime=',date_init

    date=0
    call ESMF_TimeGet (CurrTime,                           &
                       YY=date(1), MM=date(2), DD=date(3), &
                       H=date(4),  M =date(5), S =date(6), RC=rc )
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if(mype==0) write(*,'(A,6I5)') 'CurrTime =',date

    atm_int_state%Time_atmos = set_date (date(1), date(2), date(3),  &
                                         date(4), date(5), date(6))

    date_end=0
    call ESMF_TimeGet (StopTime,                                       &
                       YY=date_end(1), MM=date_end(2), DD=date_end(3), &
                       H=date_end(4),  M =date_end(5), S =date_end(6), RC=rc )
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if ( date_end(1) == 0 ) date_end = date
    atm_int_state%Time_end   = set_date (date_end(1), date_end(2), date_end(3),  &
                                         date_end(4), date_end(5), date_end(6))
    if(mype==0) write(*,'(A,6I5)') 'StopTime =',date_end
!
    call diag_manager_set_time_end(atm_int_state%Time_end)
!
    CALL ESMF_TimeIntervalGet(RunDuration, S=Run_length, RC=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
    call diag_manager_init (TIME_INIT=date)
    call diag_manager_set_time_end(atm_int_state%Time_end)
!
    atm_int_state%Time_step_atmos = set_time (dt_atmos,0)
    atm_int_state%num_atmos_calls = Run_length / dt_atmos
    atm_int_state%Time_atstart = atm_int_state%Time_atmos
    if (mype == 0) write(0,*)'num_atmos_calls=',atm_int_state%num_atmos_calls,'time_init=', &
                    date_init,'time_atmos=',date,'time_end=',date_end,'dt_atmos=',dt_atmos, &
                    'Run_length=',Run_length

! set up forecast time array that controls when to write out restart files
    frestart = 0
    call get_time(atm_int_state%Time_end - atm_int_state%Time_init,total_inttime)
! set iau offset time
    atm_int_state%Atm%iau_offset    = iau_offset
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
        atm_int_state%Time_step_restart = set_time (tmpvar, 0)
        if(iau_offset > 0 ) then
          atm_int_state%Time_restart = atm_int_state%Time_init + iautime + atm_int_state%Time_step_restart
          frestart(1) = tmpvar + iau_offset *3600
        else
          atm_int_state%Time_restart = atm_int_state%Time_init + atm_int_state%Time_step_restart
          frestart(1) = tmpvar
        endif
        if(restart_interval(1) > 0) then
          i = 2
          do while ( atm_int_state%Time_restart < atm_int_state%Time_end )
            frestart(i) = frestart(i-1) + tmpvar
            atm_int_state%Time_restart = atm_int_state%Time_restart + atm_int_state%Time_step_restart
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
    atm_int_state%intrm_rst         = 0
    if (frestart(1)>0) atm_int_state%intrm_rst = 1
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

     call  atmos_model_init (atm_int_state%Atm,  atm_int_state%Time_init, &
                             atm_int_state%Time_atmos, atm_int_state%Time_step_atmos)
!
     inquire(FILE='data_table', EXIST=fexist)
     if (fexist) then
       call data_override_init ( ) ! Atm_domain_in  = Atm%domain, &
                                   ! Ice_domain_in  = Ice%domain, &
                                   ! Land_domain_in = Land%domain )
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
      write(*, *) 'create fcst grid: mype,regional,nested=',mype,atm_int_state%Atm%regional,atm_int_state%Atm%nested

      ! regional-only without nests
      if( atm_int_state%Atm%regional .and. .not. atm_int_state%Atm%nested ) then

        call atmosphere_control_data (isc, iec, jsc, jec, nlev)

        domain   = atm_int_state%Atm%domain
        fcstNpes = atm_int_state%Atm%layout(1)*atm_int_state%Atm%layout(2)
        allocate(isl(fcstNpes), iel(fcstNpes), jsl(fcstNpes), jel(fcstNpes))
        allocate(deBlockList(2,2,fcstNpes))
        call mpp_get_compute_domains(domain,xbegin=isl,xend=iel,ybegin=jsl,yend=jel)
        do n=1,fcstNpes
           deBlockList(:,1,n) = (/ isl(n),iel(n) /)
           deBlockList(:,2,n) = (/ jsl(n),jel(n) /)
        end do
        delayout = ESMF_DELayoutCreate(petMap=(/(i,i=0,fcstNpes-1)/), rc=rc); ESMF_ERR_ABORT(rc)
        distgrid = ESMF_DistGridCreate(minIndex=(/1,1/), &
                                         maxIndex=(/atm_int_state%Atm%mlon,atm_int_state%Atm%mlat/), &
                                         delayout=delayout, &
                                         deBlockList=deBlockList, rc=rc); ESMF_ERR_ABORT(rc)

        fcstGrid = ESMF_GridCreateNoPeriDim(regDecomp=(/atm_int_state%Atm%layout(1),atm_int_state%Atm%layout(2)/), &
                                              minIndex=(/1,1/), &
                                              maxIndex=(/atm_int_state%Atm%mlon,atm_int_state%Atm%mlat/), &
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
            glonPtr(i-isc+1,j-jsc+1) = atm_int_state%Atm%lon(i-isc+1,j-jsc+1) * dtor
            glatPtr(i-isc+1,j-jsc+1) = atm_int_state%Atm%lat(i-isc+1,j-jsc+1) * dtor
          enddo
        enddo

        ! add and define "corner" coordinate values
        call ESMF_GridAddCoord(fcstGrid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
                               rc=rc); ESMF_ERR_ABORT(rc)
        call ESMF_GridGetCoord(fcstGrid, coordDim=1, staggerLoc=ESMF_STAGGERLOC_CORNER, &
                               totalLBound=tlb, totalUBound=tub, &
                               farrayPtr=glonPtr, rc=rc); ESMF_ERR_ABORT(rc)
        glonPtr(tlb(1):tub(1),tlb(2):tub(2)) = &
           atm_int_state%Atm%lon_bnd(tlb(1):tub(1),tlb(2):tub(2)) * dtor
        call ESMF_GridGetCoord(fcstGrid, coordDim=2, staggerLoc=ESMF_STAGGERLOC_CORNER, &
                               totalLBound=tlb, totalUBound=tub, &
                               farrayPtr=glatPtr, rc=rc); ESMF_ERR_ABORT(rc)
        glatPtr(tlb(1):tub(1),tlb(2):tub(2)) = &
          atm_int_state%Atm%lat_bnd(tlb(1):tub(1),tlb(2):tub(2)) * dtor

        call mpp_error(NOTE, 'after create fcst grid for regional-only')

      else ! not regional only

        if (.not. atm_int_state%Atm%regional .and. .not. atm_int_state%Atm%nested ) then  !! global only

          do tl=1,6
              decomptile(1,tl) = atm_int_state%Atm%layout(1)
              decomptile(2,tl) = atm_int_state%Atm%layout(2)
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

          if (mype==0) TileLayout = atm_int_state%Atm%layout
          call ESMF_VMBroadcast(vm, bcstData=TileLayout, count=2, &
                                  rootPet=0, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if (mype==0) npes(1) = mpp_npes()
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
            call mpp_error(FATAL, 'Inconsistent nestRootPet and Atm%layout')
          endif

          if (mype == nestRootPet) then
            if (nestRootPet /= atm_int_state%Atm%pelist(1)) then
              write(0,*)'error in fcst_initialize: nestRootPet /= atm_int_state%Atm%pelist(1)'
              write(0,*)'error in fcst_initialize: nestRootPet = ',nestRootPet
              write(0,*)'error in fcst_initialize: atm_int_state%Atm%pelist(1) = ',atm_int_state%Atm%pelist(1)
              ESMF_ERR_ABORT(100)
            endif
          endif

          ! nest rootPet shares peList with others
          if (mype == nestRootPet) peListSize(1) = size(atm_int_state%Atm%pelist)
          call ESMF_VMBroadcast(vm, bcstData=peListSize, count=1, rootPet=nestRootPet, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          ! nest rootPet shares layout with others
          if (mype == nestRootPet) regDecomp = atm_int_state%Atm%layout
          call ESMF_VMBroadcast(vm, bcstData=regDecomp, count=2, rootPet=nestRootPet, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          ! prepare petMap variable
          allocate(petMap(peListSize(1)))
          if (mype == nestRootPet) petMap = atm_int_state%Atm%pelist
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
      if ( .not. atm_int_state%Atm%nested ) then  !! global only
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
      if (mype == 0) write(0,*)'dateS=',trim(dateS),'date_init=',date_init

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="time:units", value=trim(dateS), rc=rc)
!                              name="time:units", value="hours since 2016-10-03 00:00:00", rc=rc)
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
           if (mype == 0) write(0,*)'af create fcst fieldbundle, name=',trim(name_FB),'rc=',rc
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

           call fv_dyn_bundle_setup(atm_int_state%Atm%axes,          &
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
             if (mype == 0) write(0,*)'af create fcst fieldbundle, name=',trim(name_FB1),'rc=',rc
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
           enddo
!
           call fv_phys_bundle_setup(atm_int_state%Atm%diag, atm_int_state%Atm%axes, &
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
!
!-----------------------------------------------------------------------
!
      IF(rc /= ESMF_SUCCESS) THEN
        WRITE(0,*)"FAIL: Fcst_Initialize."
!      ELSE
!        WRITE(0,*)"PASS: Fcst_Initialize."
      ENDIF
!
      if (mype == 0) write(0,*)'in fcst,init total time: ', mpi_wtime() - timeis
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
!-----------------------------------------------------------------------
!***  local variables
!
      integer                    :: i,j, mype, na, date(6)
      character(20)              :: compname

      type(ESMF_Time)            :: currtime
      integer(kind=ESMF_KIND_I8) :: ntimestep_esmf
      character(len=64)          :: timestamp
!
!-----------------------------------------------------------------------
!
      real(kind=8)   :: mpi_wtime, tbeg1
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      tbeg1 = mpi_wtime()
      rc    = esmf_success
!
!-----------------------------------------------------------------------
!
      call ESMF_GridCompGet(fcst_comp, name=compname, localpet=mype, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
      call ESMF_ClockGet(clock, advanceCount=NTIMESTEP_ESMF, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      na = NTIMESTEP_ESMF
!
!-----------------------------------------------------------------------
! *** call fcst integration subroutines

      call get_date (atm_int_state%Time_atmos, date(1), date(2), date(3),  &
                     date(4), date(5), date(6))
      atm_int_state%Time_atmos = atm_int_state%Time_atmos + atm_int_state%Time_step_atmos

      call update_atmos_model_dynamics (atm_int_state%Atm)

      call update_atmos_radiation_physics (atm_int_state%Atm)

      call atmos_model_exchange_phase_1 (atm_int_state%Atm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!-----------------------------------------------------------------------
!
!     IF(RC /= ESMF_SUCCESS) THEN
!       if(mype==0) WRITE(0,*)"FAIL: fcst_RUN"
!      ELSE
        if(mype==0) WRITE(*,*)"PASS: fcstRUN phase 1, na = ",na, ' time is ', mpi_wtime()-tbeg1
!     ENDIF
!
!-----------------------------------------------------------------------
!
   end subroutine fcst_run_phase_1
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
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
!-----------------------------------------------------------------------
!***  local variables
!
      integer                    :: i,j, mype, na, date(6), seconds
      character(20)              :: compname

      type(time_type)            :: restart_inctime
      type(ESMF_Time)            :: currtime
      integer(kind=ESMF_KIND_I8) :: ntimestep_esmf
      character(len=64)          :: timestamp
!
!-----------------------------------------------------------------------
!
      real(kind=8)   :: mpi_wtime, tbeg1
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      tbeg1 = mpi_wtime()
      rc    = esmf_success
!
!-----------------------------------------------------------------------
!
      call ESMF_GridCompGet(fcst_comp, name=compname, localpet=mype, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
      call ESMF_ClockGet(clock, advanceCount=NTIMESTEP_ESMF, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      na = NTIMESTEP_ESMF
      if (mype == 0) write(0,*)'in fcst run phase 2, na=',na
!
!-----------------------------------------------------------------------
! *** call fcst integration subroutines

      call atmos_model_exchange_phase_2 (atm_int_state%Atm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call update_atmos_model_state (atm_int_state%Atm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!--- intermediate restart
      if (atm_int_state%intrm_rst>0) then
        if (na /= atm_int_state%num_atmos_calls-1) then
          call get_time(atm_int_state%Time_atmos - atm_int_state%Time_init, seconds)
          if (ANY(frestart(:) == seconds)) then
            if (mype == 0) write(0,*)'write out restart at na=',na,' seconds=',seconds,  &
               'integration lenght=',na*dt_atmos/3600.
            timestamp = date_to_string (atm_int_state%Time_atmos)
            call atmos_model_restart(atm_int_state%Atm, timestamp)
            call write_stoch_restart_atm('RESTART/'//trim(timestamp)//'.atm_stoch.res.nc')

            call wrt_atmres_timestamp(atm_int_state,timestamp)
          endif
        endif
      endif
!
!-----------------------------------------------------------------------
!
!     IF(RC /= ESMF_SUCCESS) THEN
!      if(mype==0) WRITE(0,*)"FAIL: fcst_RUN"
!      ELSE
       if(mype==0) WRITE(*,*)"PASS: fcstRUN phase 2, na = ",na, ' time is ', mpi_wtime()-tbeg1
!     ENDIF
!
!-----------------------------------------------------------------------
!
   end subroutine fcst_run_phase_2
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
   subroutine fcst_finalize(fcst_comp, importState, exportState,clock,rc)
!
!-----------------------------------------------------------------------
!***  finalize the forecast grid component.
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp)            :: fcst_comp
      type(ESMF_State)               :: importState, exportState
      type(ESMF_Clock)               :: clock
      integer,intent(out)            :: rc
!
!***  local variables
!
      integer :: unit
      integer,dimension(6)           :: date

      real(8) mpi_wtime, tfs, tfe
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      tfs = mpi_wtime()
      rc  = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  retrieve the fcst component's esmf internal state
!-----------------------------------------------------------------------
!
      call ESMF_GridCompGetInternalState(fcst_comp, wrap, rc)
      atm_int_state => wrap%ptr
!
!-----------------------------------------------------------------------
!
      call atmos_model_end (atm_int_state%atm)
!
!*** check time versus expected ending time

      if (atm_int_state%Time_atmos /= atm_int_state%Time_end)  &
        call error_mesg ('program coupler',  &
                         'final time does not match expected ending time', WARNING)

!*** write restart file
      if( restart_endfcst ) then
        call get_date (atm_int_state%Time_atmos, date(1), date(2), date(3),  &
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
!
      call diag_manager_end(atm_int_state%Time_atmos )

      call fms_end
!
!-----------------------------------------------------------------------
!
      IF(RC /= ESMF_SUCCESS)THEN
        WRITE(0,*)'FAIL: Write_Finalize.'
!      ELSE
!        WRITE(0,*)'PASS: Write_Finalize.'
      ENDIF
!
      tfe = mpi_wtime()
!      print *,'fms end time: ', tfe-tfs
!-----------------------------------------------------------------------
!
  end subroutine fcst_finalize
!
!#######################################################################
!-- change name from coupler_res to wrt_res_stamp to avoid confusion,
!-- here we only write out atmos restart time stamp
!
  subroutine wrt_atmres_timestamp(atm_int_state,timestamp)
    type(atmos_internalstate_type), intent(in) :: atm_int_state
    character(len=32), intent(in) :: timestamp

    integer :: unit, date(6)

!----- compute current date ------

    call get_date (atm_int_state%Time_atmos, date(1), date(2), date(3),  &
                                 date(4), date(5), date(6))

!----- write restart file ------

    if (mpp_pe() == mpp_root_pe())then
        call mpp_open( unit, 'RESTART/'//trim(timestamp)//'.coupler.res', nohdrs=.TRUE. )
        write( unit, '(i6,8x,a)' )calendar_type, &
             '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'

        write( unit, '(6i6,8x,a)' )date_init, &
             'Model start time:   year, month, day, hour, minute, second'
        write( unit, '(6i6,8x,a)' )date, &
             'Current model time: year, month, day, hour, minute, second'
        call mpp_close(unit)
    endif
  end subroutine wrt_atmres_timestamp
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
!-----------------------------------------------------------------------
!***  local variables
!
    logical                     :: ioCapable
    logical                     :: doItFlag
    character(len=64)           :: lfileName
    character(len=64)           :: gridName
    type(ESMF_Array)            :: array
    type(ESMF_ArrayBundle)      :: arraybundle
    logical                     :: isPresent
    integer                     :: stat
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
