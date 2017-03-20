!--------------- FV3GFS solo model -----------------
!
!*** The FV3 atmosphere grid component nuopc cap
!
! Author:  Jun Wang@noaa.gov
!
! revision history
! 11 Oct 2016: J. Wang          Initial code
!

module fv3gfs_cap_mod

  use time_manager_mod,  only: time_type, set_calendar_type, set_time,    &
                               set_date, days_in_month, month_name,       &
                               operator(+), operator (<), operator (>),   &
                               operator (/=), operator (/), operator (==),&
                               operator (*), THIRTY_DAY_MONTHS, JULIAN,   &
                               NOLEAP, NO_CALENDAR, date_to_string,       &
                               get_date

  use  atmos_model_mod,  only: atmos_model_init, atmos_model_end,  &
                               update_atmos_model_dynamics,        &
                               update_atmos_radiation_physics,     &
                               update_atmos_model_state,           &
                               atmos_data_type, atmos_model_restart

  use constants_mod,     only: constants_init
  use       fms_mod,     only: open_namelist_file, file_exist, check_nml_error,  &
                               error_mesg, fms_init, fms_end, close_file,        &
                               write_version_number, uppercase

  use mpp_mod,           only: mpp_init, mpp_pe, mpp_root_pe, mpp_npes, mpp_get_current_pelist, &
                               mpp_set_current_pelist, stdlog, mpp_error, NOTE, FATAL, WARNING
  use mpp_mod,           only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, mpp_sync

  use mpp_io_mod,        only: mpp_open, mpp_close, &
                               MPP_NATIVE, MPP_RDONLY, MPP_DELETE

  use mpp_domains_mod,   only: mpp_get_global_domain, mpp_global_field, CORNER
  use memutils_mod,      only: print_memuse_stats
  use sat_vapor_pres_mod,only: sat_vapor_pres_init

  use  diag_manager_mod, only: diag_manager_init, diag_manager_end, &
                               get_base_date, diag_manager_set_time_end

  use data_override_mod, only: data_override_init


  use ESMF
  use NUOPC
  use NUOPC_Model,       only: &
    model_routine_SS      => SetServices, &
    model_label_Advance   => label_Advance, &
    model_label_Finalize  => label_Finalize

  use time_utils_mod

  implicit none
  private
  public SetServices

!-----------------------------------------------------------------------

  character(len=128) :: version = '$Id: coupler_main.F90,v 19.0.4.1.2.3 2014/09/09 23:51:59 Rusty.Benson Exp $'
  character(len=128) :: tag = '$Name: ulm_201505 $'

!---- model defined-types ----

  type atmos_internalstate_type
    type(atmos_data_type)  :: Atm
    type (time_type)       :: Time_atmos, Time_init, Time_end,  &
                              Time_step_atmos, Time_step_ocean, &
                              Time_restart, Time_step_restart
    integer :: num_atmos_calls, ret, intrm_rst
  end type

  type atmos_internalstate_wrapper
    type(atmos_internalstate_type), pointer :: ptr
  end type

  type(atmos_internalstate_type),pointer,save :: atm_int_state
  type(atmos_internalstate_wrapper),save      :: wrap
  type(ESMF_Clock),save                       :: clock_fv3

!----- coupled model date -----

  integer :: date_init(6)
  integer :: calendar_type = -99
  integer :: restart_interval
  logical :: profile_memory = .true.

!-----------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------
  !------------------- Solo fv3gfs code starts here ----------------------
  !-----------------------------------------------------------------------

  subroutine SetServices(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname='(fv3gfs_cap:SetServices)'

    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! initialization, switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeRealize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! model advance method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! model finalize method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
      specRoutine=atmos_model_finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine SetServices

  !-----------------------------------------------------------------------------

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    
    character(len=10)                         :: value

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv00p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_AttributeGet(gcomp, name="ProfileMemory", value=value, defaultValue="true", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    profile_memory=(trim(value)/="false")


!    ! A restart_interval value of 0 means no restart will be written.
!    call ESMF_AttributeGet(gcomp, name="restart_interval", value=value, defaultValue="0", &
!      convention="NUOPC", purpose="Instance", rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
!
!    restart_interval = ESMF_UtilString2Int(value, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

!    if(restart_interval < 0) then
!      call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
!        msg="FV3_CAP: ATM attribute: restart_interval cannot be negative.", &
!        line=__LINE__, &
!        file=__FILE__, rcToReturn=rc)
!      return
!    endif
!    call ESMF_LogWrite('FV3_CAP:restart_interval = '//trim(value), ESMF_LOGMSG_INFO, rc=dbrc)  
    
  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    type(ESMF_GridComp)                    :: gcomp
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Clock)                       :: clock
    integer, intent(out)                   :: rc

    type(ESMF_VM)                          :: vm
    type(ESMF_Time)                        :: CurrTime, TINI, StopTime
    type(ESMF_TimeInterval)                :: TINT, RunDuration
    type(ESMF_Config)                      :: cf


    type(time_type)                        :: Run_len      ! length of experiment
    type(time_type)                        :: Time
    type(time_type)                        :: Time_restart
    type(time_type)                        :: DT
    integer                                :: isc,iec,jsc,jec
    integer                                :: dt_atmos,Run_length
    integer,dimension(6)                   :: date, date_end
    integer                                :: mpi_comm_atm
    integer                                :: res_intvl

    type(ESMF_Grid)                        :: gridIn
    type(ESMF_Grid)                        :: gridOut

    integer                                :: npet, npet_x, npet_y
    logical                                ::  force_date_from_configure 
    character(len=*),parameter  :: subname='(mom_cap:InitializeAdvertise)'
    character(len=9) :: month
    character(len=17) :: calendar = '                 '
    integer :: initClock, mainClock, termClock, unit, gnlon, gnlat
    real,    allocatable, dimension(:,:) :: glon_bnd, glat_bnd

!jw debug
    character(ESMF_MAXSTR)          :: name
    integer ::petcount,tasks,mype,tstep,date1(6)
    type (time_type)       :: Time_init1


    rc = ESMF_SUCCESS

    call ESMF_GridCompGet(gcomp,name=name,vm=vm,rc=rc)
!    print *,'infv3_cap,initAdvertize,name=',trim(name),'rc=',rc
!    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm, mpiCommunicator=mpi_comm_atm,petCount=petcount, &
             localpet = mype,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    clock_fv3=clock
!
    call ESMF_ClockGet(clock_fv3, currTIME=CurrTime, TimeStep=TINT,      &
                       StartTime=TINI, RunDuration=RunDuration, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    CALL ESMF_TimeIntervalGet(TINT,s=tstep,rc=RC)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    date=0
    call ESMF_TimeGet (CurrTime,                    &
                       YY=date(1), MM=date(2), DD=date(3), &
                       H=date(4),  M =date(5), S =date(6), RC=rc )
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if(mype==0) print *,'CurrTime=',date

    date_init=0
    call ESMF_TimeGet (TINI,                      &
                       YY=date_init(1), MM=date_init(2), DD=date_init(3), &
                       H=date_init(4),  M =date_init(5), S =date_init(6), RC=rc )
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if ( date_init(1) == 0 ) date_init = date
    if(mype==0) print *,'InitTime=',date_init

    StopTime = CurrTime + RunDuration
    date_end=0
    call ESMF_TimeGet (StopTime,                      &
                       YY=date_end(1), MM=date_end(2), DD=date_end(3), &
                       H=date_end(4),  M =date_end(5), S =date_end(6), RC=rc )
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if ( date_end(1) == 0 ) date_end = date
    if(mype==0) print *,'StopTime=',date_end,'date=',date

!
!    CALL ESMF_GridCompGet(gcomp, config=CF, rc=RC)
    CF=ESMF_ConfigCreate(rc=RC)
    CALL ESMF_ConfigLoadFile(config=CF ,filename='model_configure' ,rc=RC)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!
    CALL ESMF_ConfigGetAttribute(config=CF,value=restart_interval, &
                                 label ='restart_interval:',rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!
    CALL ESMF_ConfigGetAttribute(config=CF,value=calendar, &
                                 label ='calendar:',rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!
    CALL ESMF_ConfigGetAttribute(config=CF,value=dt_atmos, &
                                 label ='dt_atmos:',rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

! set internal state
    allocate(atm_int_state, stat = rc)
    wrap%ptr => atm_int_state
!
    call ESMF_GridCompSetInternalState(gcomp, wrap, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!#######################################################################
! set fvs variables
!
     call fms_init(mpi_comm_atm)
     call mpp_init()
     initClock = mpp_clock_id( 'Initialization' )
     call mpp_clock_begin (initClock) !nesting problem

     call fms_init
     call constants_init
     call sat_vapor_pres_init
!
     if (file_exist('INPUT/coupler.res')) then
       call mpp_open( unit, 'INPUT/coupler.res', action=MPP_RDONLY )
       read (unit,*,err=999) calendar_type
       read (unit,*) date_init
       read (unit,*) date
       goto 998 !back to fortran-4
     ! read old-style coupler.res
   999 call mpp_close (unit)
       call mpp_open (unit, 'INPUT/coupler.res', action=MPP_RDONLY, form=MPP_NATIVE)
       read (unit) calendar_type
       read (unit) date
   998 call mpp_close(unit)
       force_date_from_configure = .false.
       if(mype==0) print *,'restart, date=',date
      else
!-- set restart date from configure file
!-- need to add current date in model_configure
!-- otherwise, assume it is not a restart run
       force_date_from_configure = .true.
      endif

      if ( force_date_from_configure ) then

       select case( uppercase(trim(calendar)) )
       case( 'JULIAN' )
           calendar_type = JULIAN
       case( 'NOLEAP' )
           calendar_type = NOLEAP
       case( 'THIRTY_DAY' )
           calendar_type = THIRTY_DAY_MONTHS
       case( 'NO_CALENDAR' )
           calendar_type = NO_CALENDAR
       case default
           call mpp_error ( FATAL, 'COUPLER_MAIN: coupler_nml entry calendar must '// &
                                   'be one of JULIAN|NOLEAP|THIRTY_DAY|NO_CALENDAR.' )
       end select

      endif

      call set_calendar_type (calendar_type         )

      call diag_manager_init (TIME_INIT=date)
!
      Time_init1  = set_date (date_init(1), date_init(2), date_init(3), &
                            date_init(4), date_init(5), date_init(6))
      atm_int_state%Time_init  = set_date (date_init(1), date_init(2), date_init(3), &
                            date_init(4), date_init(5), date_init(6))

      atm_int_state%Time_atmos = set_date (date(1), date(2), date(3),  &
                            date(4), date(5), date(6))
      call ESMF_TimeSet(time=CurrTime,yy=date(1),mm=date(2),dd=date(3),h=date(4), &
                       m=date(5),s=date(6),rc=rc)
!reset CurrTime in clock
      call ESMF_ClockSet(clock, currTIME=CurrTime, rc=rc)
!jw test
      date1=0
      call ESMF_TimeGet (CurrTime,                    &
                        YY=date1(1), MM=date1(2), DD=date1(3), &
                        H=date1(4),  M =date1(5), S =date1(6), RC=rc )
      if(mype==0) print *,'aft rsstm date1=',date1
      RunDuration = StopTime - CurrTime
!
      CALL ESMF_TimeIntervalGet(RunDuration, S=Run_length, RC=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
!
      atm_int_state%Time_end   = set_date (date_end(1), date_end(2), date_end(3),  &
                            date_end(4), date_end(5), date_end(6))
!
      call diag_manager_set_time_end(atm_int_state%Time_end)
!-----------------------------------------------------------------------
!----- write time stamps (for start time and end time) ------

      call mpp_open( unit, 'time_stamp.out', nohdrs=.TRUE. )
      month = month_name(date(2))
      if ( mpp_pe() == mpp_root_pe() ) write (unit,20) date, month(1:3)
      month = month_name(date_end(2))
      if ( mpp_pe() == mpp_root_pe() ) write (unit,20) date_end, month(1:3)
      call mpp_close (unit)
  20  format (6i4,2x,a3)

!
     atm_int_state%Time_step_atmos = set_time (dt_atmos,0)
     atm_int_state%num_atmos_calls = Run_length / dt_atmos
     if(mype==0) print *,'num_atmos_calls=',atm_int_state%num_atmos_calls,'time_init=', &
       date_init,'time_atmos=',date,'time_end=',date_end,'dt_atmos=',dt_atmos, &
       'Run_length=',Run_length
     res_intvl=restart_interval*3600
     atm_int_state%Time_step_restart = set_time (res_intvl, 0)
     atm_int_state%Time_restart = atm_int_state%Time_atmos + atm_int_state%Time_step_restart
     atm_int_state%intrm_rst = .false.
     if (res_intvl>0) atm_int_state%intrm_rst = .true.
!
!------ initialize component models ------

     call  atmos_model_init (atm_int_state%Atm,  atm_int_state%Time_init, &
                             atm_int_state%Time_atmos, atm_int_state%Time_step_atmos)
!
     call mpp_get_global_domain(atm_int_state%Atm%Domain, xsize=gnlon, ysize=gnlat)
     allocate ( glon_bnd(gnlon+1,gnlat+1), glat_bnd(gnlon+1,gnlat+1) )
     call mpp_global_field(atm_int_state%Atm%Domain, atm_int_state%Atm%lon_bnd, glon_bnd, position=CORNER)
     call mpp_global_field(atm_int_state%Atm%Domain, atm_int_state%Atm%lat_bnd, glat_bnd, position=CORNER)

     call data_override_init ( ) ! Atm_domain_in  = Atm%domain, &
                                 ! Ice_domain_in  = Ice%domain, &
                                 ! Land_domain_in = Land%domain )

!-----------------------------------------------------------------------
!---- open and close dummy file in restart dir to check if dir exists --

      if (mpp_pe() == 0 ) then
         call mpp_open( unit, 'RESTART/file' )
         call mpp_close(unit, MPP_DELETE)
      endif
!
!-----------------------------------------------------------------------
! atmos clock may need to update for coupled system as the clock might be
! changed from restart. For the time bring, it is OK for standalone fv3
!

  end subroutine InitializeAdvertise
!
  !-----------------------------------------------------------------------------

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Grid)         :: grid

    rc = ESMF_SUCCESS

    ! nothing is realized in the import/export States

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)                    :: gcomp
    integer, intent(out)                   :: rc
    
    ! local variables
    type(ESMF_Clock)                       :: clock
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Time)                        :: currTime
    type(ESMF_TimeInterval)                :: timeStep
    type(ESMF_Time)                        :: startTime
    type(ESMF_TimeInterval)                :: time_elapsed
    integer(ESMF_KIND_I8)                  :: n_interval, time_elapsed_sec
    character(len=64)                      :: timestamp

    ! define some time types 
    type(time_type)                        :: Time        

    integer :: na,i,j,i1,j1,dth,dtm,dts
    character(len=*),parameter  :: subname='(fv3_cap:ModelAdvance)'
    character(240)              :: msgString
!jw debug
    character(ESMF_MAXSTR)          :: name
    type(ESMF_VM)                :: vm
    integer :: mype,date(6)


    rc = ESMF_SUCCESS
    if(profile_memory) call ESMF_VMLogMemInfo("Entering FV3 Model_ADVANCE: ")
    
    call ESMF_GridCompGet(gcomp,name=name,vm=vm,rc=rc)
!    print *,'infv3_cap,initAdvertize,name=',trim(name),'rc=',rc
!    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm, localpet = mype,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
!    if(mype==0) print *,'n fv3_cap,in model run,compget,rc=',rc
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!    call ESMF_GridCompGetInternalState(gcomp, wrap, rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
!
!    Atm_int_state => wrap%ptr

    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in 
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.
    
    call ESMF_ClockPrint(clock_fv3, options="currTime", &
      preString="------>Advancing FV3 from: ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_ClockGet(clock_fv3, startTime=startTime, currTime=currTime, &
      timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_TimePrint(currTime + timeStep, &
      preString="--------------------------------> to: ", &
      unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_TimeIntervalGet(timeStep, h=dth, m=dtm, s=dts, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!    if(mype==0) print *,'n fv3_cap,in model run,timestep,rc=',rc,'dth=',dth, &
!      'dtm=',dtm, 'dts=',dts,'num_atmos_calls=',atm_int_state%num_atmos_calls

!---
    atm_int_state%Time_atmos = esmf2fms_time(currTime)
    call get_date (atm_int_state%Time_atmos, date(1), date(2), date(3),  &
                                 date(4), date(5), date(6))
    call get_date (atm_int_state%Time_step_atmos, date(1), date(2), date(3),  &
                                 date(4), date(5), date(6))
!
    do na = 1, atm_int_state%num_atmos_calls

      atm_int_state%Time_atmos = atm_int_state%Time_atmos + atm_int_state%Time_step_atmos

      call update_atmos_model_dynamics (atm_int_state%Atm)

      call update_atmos_radiation_physics (atm_int_state%Atm)

      call update_atmos_model_state (atm_int_state%Atm)

!--- intermediate restart
      if (atm_int_state%intrm_rst) then
        if ((na /= atm_int_state%num_atmos_calls) .and.   &
           (atm_int_state%Time_atmos == atm_int_state%Time_restart)) then
          timestamp = date_to_string (atm_int_state%Time_restart)
          call atmos_model_restart(atm_int_state%Atm, timestamp)

          call wrt_atmres_timestamp(atm_int_state,timestamp)
          atm_int_state%Time_restart = atm_int_state%Time_restart + atm_int_state%Time_step_restart
        endif
      endif

      call print_memuse_stats('after full step')
!    if(mype==0) print *,'n fv3_cap,in model run,end one time step,na=',na

    enddo

    if(profile_memory) call ESMF_VMLogMemInfo("Leaving FV3 Model_ADVANCE: ")
  end subroutine ModelAdvance

  !-----------------------------------------------------------------------------

  subroutine atmos_model_finalize(gcomp, rc)

    ! input arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
!
    integer :: unit, date(6)
    
!    print *,'FV3: --- finalize called ---'
    rc = ESMF_SUCCESS

    call atmos_model_end (atm_int_state%atm)
!
!----- check time versus expected ending time ----

    if (atm_int_state%Time_atmos /= atm_int_state%Time_end)  &
      call error_mesg ('program coupler',  &
           'final time does not match expected ending time', WARNING)

!----- write restart file ------

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
!
    call diag_manager_end(atm_int_state%Time_atmos )

    call fms_end

!    print *,'FV3: --- finalize completed ---'

  end subroutine atmos_model_finalize

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

  !-----------------------------------------------------------------------------

end module fv3gfs_cap_mod
