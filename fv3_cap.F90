!--------------- FV3GFS solo model -----------------
!
!*** The FV3 atmosphere grid component nuopc cap
!
! Author:  Jun Wang@noaa.gov
!
! revision history
! 11 Oct 2016: J. Wang          Initial code
! 18 Apr 2017: J. Wang          set up fcst grid component and write grid components
! 24 Jul 2017: J. Wang          initialization and time stepping changes for coupling
! 02 Nov 2017: J. Wang          Use Gerhard's transferable RouteHandle
!

module fv3gfs_cap_mod

  use ESMF
  use NUOPC
  use NUOPC_Model,            only: model_routine_SS => SetServices,         &
                                    SetVM,                                   &
                                    routine_Run,                             &
                                    label_Advertise,                         &
                                    label_RealizeProvided,                   &
                                    label_Advance,                           &
                                    label_CheckImport,                       &
                                    label_SetRunClock,                       &
                                    label_TimestampExport,                   &
                                    label_Finalize,                          &
                                    NUOPC_ModelGet
!
  use module_fv3_config,      only: quilting, output_fh,                     &
                                    nfhout, nfhout_hf, nsout, dt_atmos,      &
                                    calendar,                                &
                                    cplprint_flag,output_1st_tstep_rst,      &
                                    first_kdt

  use module_fv3_io_def,      only: num_pes_fcst,write_groups,               &
                                    num_files, filename_base,                &
                                    wrttasks_per_group, n_group,             &
                                    lead_wrttask, last_wrttask,              &
                                    nsout_io, iau_offset, lflname_fulltime
!
  use module_fcst_grid_comp,  only: fcstSS => SetServices,                   &
                                    fcstGrid, numLevels, numSoilLayers,      &
                                    numTracers, mygrid, grid_number_on_all_pets

  use module_wrt_grid_comp,   only: wrtSS => SetServices
!
  use module_cplfields,       only: nExportFields, exportFields, exportFieldsInfo, &
                                    nImportFields, importFields, importFieldsInfo, &
                                    importFieldsValid, queryImportFields

  use module_cplfields,       only: realizeConnectedCplFields
  use module_cap_cpl,         only: diagnose_cplFields

  use atmos_model_mod,        only: setup_exportdata

  implicit none
  private
  public SetServices
!
!-----------------------------------------------------------------------
!

  type(ESMF_GridComp)                         :: fcstComp
  type(ESMF_State)                            :: fcstState
  type(ESMF_FieldBundle), allocatable         :: fcstFB(:)
  integer, save                               :: FBCount

  type(ESMF_GridComp),    allocatable         :: wrtComp(:)
  type(ESMF_State),       allocatable         :: wrtState(:)
  type(ESMF_FieldBundle), allocatable         :: wrtFB(:,:)

  type(ESMF_RouteHandle), allocatable         :: routehandle(:,:)

  logical                                     :: profile_memory = .true.

  integer                                     :: mype = -1
  integer                                     :: dbug = 0

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
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSpecialize(gcomp, specLabel=label_Advertise, specRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=label_RealizeProvided, specRoutine=InitializeRealize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! model advance method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=label_Advance, &
                              specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! checking the import fields is a bit more complex because of coldstart option
    call ESMF_MethodRemove(gcomp, label_CheckImport, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=label_CheckImport, &
                              specRoutine=fv3_checkimport, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! setup Run/Advance phase: phase1
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
                                 phaseLabelList=(/"phase1"/), userRoutine=routine_Run, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=label_Advance, &
                              specPhaseLabel="phase1", specRoutine=ModelAdvance_phase1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! setup Run/Advance phase: phase2
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
                                 phaseLabelList=(/"phase2"/), userRoutine=routine_Run, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=label_Advance, &
                              specPhaseLabel="phase2", specRoutine=ModelAdvance_phase2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! specializations to set fv3 cap run clock (model clock)
    call ESMF_MethodRemove(gcomp, label=label_SetRunClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=label_SetRunClock, &
                                     specRoutine=ModelSetRunClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! specializations required to support 'inline' run sequences
    call NUOPC_CompSpecialize(gcomp, specLabel=label_CheckImport, &
                              specPhaseLabel="phase1", specRoutine=fv3_checkimport, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=label_TimestampExport, &
                              specPhaseLabel="phase1", specRoutine=TimestampExport_phase1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=label_CheckImport, &
                              specPhaseLabel="phase2", specRoutine=NUOPC_NoOp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! model finalize method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=label_Finalize, &
                              specRoutine=ModelFinalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine SetServices

!-----------------------------------------------------------------------------

  subroutine InitializeAdvertise(gcomp, rc)

    type(ESMF_GridComp)                    :: gcomp
    integer, intent(out)                   :: rc

! local variables
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Clock)                       :: clock

    character(len=10)                      :: value
    character(240)                         :: msgString
    logical                                :: isPresent, isSet
    type(ESMF_VM)                          :: vm, fcstVM
    type(ESMF_Time)                        :: currTime, startTime
    type(ESMF_TimeInterval)                :: timeStep, rsthour
    type(ESMF_Config)                      :: cf
    type(ESMF_RegridMethod_Flag)           :: regridmethod

    integer                                :: i, j, k, urc, ist
    integer                                :: noutput_fh, nfh, nfh2
    integer                                :: petcount
    integer                                :: nfhmax_hf
    real                                   :: nfhmax
    real                                   :: output_startfh, outputfh, outputfh2(2)
    logical                                :: loutput_fh, lfreq
    character(ESMF_MAXSTR)                 :: name
    integer,dimension(:), allocatable      :: petList, fcstPetList, originPetList, targetPetList
    character(len=esmf_maxstr),allocatable :: fcstItemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: fcstItemTypeList(:)
    character(20)                          :: cwrtcomp
    integer                                :: isrcTermProcessing
    type(ESMF_Info)                        :: parentInfo, childInfo

    character(len=*),parameter             :: subname='(fv3_cap:InitializeAdvertise)'
    real(kind=8)                           :: MPI_Wtime, timeis, timerhs
!
!------------------------------------------------------------------------
!
    rc = ESMF_SUCCESS
    timeis = MPI_Wtime()

    call ESMF_GridCompGet(gcomp,name=name,vm=vm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_VMGet(vm, petCount=petcount, localpet=mype, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! query for importState and exportState
    call NUOPC_ModelGet(gcomp, driverClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeGet(gcomp, name="ProfileMemory", value=value, defaultValue="false", &
                           convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    profile_memory = (trim(value)/="false")

    call ESMF_AttributeGet(gcomp, name="DumpFields", value=value, defaultValue="false", &
                           convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    cplprint_flag = (trim(value)=="true")
    write(msgString,'(A,l6)') trim(subname)//' cplprint_flag = ',cplprint_flag
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)

    ! Read in cap debug flag
    call NUOPC_CompAttributeGet(gcomp, name='dbug_flag', value=value, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (isPresent .and. isSet) then
     read(value,*) dbug
    end if
    write(msgString,'(A,i6)') trim(subname)//' dbug = ',dbug
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)

!------------------------------------------------------------------------
! get config variables
!
    CF = ESMF_ConfigCreate(rc=rc)
    call ESMF_ConfigLoadFile(config=CF ,filename='model_configure' ,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
    call ESMF_ConfigGetAttribute(config=CF,value=calendar, &
                                 label ='calendar:', &
                                 default='gregorian',rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
    call ESMF_ConfigGetAttribute(config=CF,value=quilting, &
                                 label ='quilting:',rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_ConfigGetAttribute(config=CF,value=output_1st_tstep_rst, &
                                 default=.false., label ='output_1st_tstep_rst:',rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_ConfigGetAttribute(config=CF,value=iau_offset,default=0,label ='iau_offset:',rc=rc)
    if (iau_offset < 0) iau_offset=0

    noutput_fh = ESMF_ConfigGetLen(config=CF, label ='output_fh:',rc=rc)

    if(mype == 0) print *,'af nems config,quilting=',quilting,' calendar=', trim(calendar),' iau_offset=',iau_offset, &
      ' noutput_fh=',noutput_fh
!
    nfhout = 0 ; nfhmax_hf = 0 ; nfhout_hf = 0 ; nsout = 0
    if ( quilting ) then
      call ESMF_ConfigGetAttribute(config=CF,value=write_groups, &
                                   label ='write_groups:',rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
      call ESMF_ConfigGetAttribute(config=CF,value=wrttasks_per_group, &
                                   label ='write_tasks_per_group:',rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_ConfigGetAttribute(config=CF,value=isrcTermProcessing, default=0, &
                                   label ='isrcTermProcessing:',rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if(mype == 0) print *,'af nems config,quilting=',quilting,' write_groups=', &
        write_groups,wrttasks_per_group,' isrcTermProcessing=', isrcTermProcessing
!
      call ESMF_ConfigGetAttribute(config=CF,value=num_files, &
                                   label ='num_files:',rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
      allocate(filename_base(num_files))
      call ESMF_ConfigFindLabel(CF,'filename_base:',rc=rc)
      do i=1,num_files
        call ESMF_ConfigGetAttribute(config=CF,value=filename_base(i), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      enddo

! variables for output
      call ESMF_ConfigGetAttribute(config=CF, value=nfhout,   label ='nfhout:',   default=-1,rc=rc)
      call ESMF_ConfigGetAttribute(config=CF, value=nfhmax_hf,label ='nfhmax_hf:',default=-1,rc=rc)
      call ESMF_ConfigGetAttribute(config=CF, value=nfhout_hf,label ='nfhout_hf:',default=-1,rc=rc)
      call ESMF_ConfigGetAttribute(config=CF, value=nsout,    label ='nsout:',    default=-1,rc=rc)
      nsout_io = nsout
!
      if(mype==0) print *,'af nems config,nfhout,nsout=',nfhout,nfhmax_hf,nfhout_hf, nsout,noutput_fh

    endif ! quilting
!
    call ESMF_ConfigGetAttribute(config=CF, value=dt_atmos, label ='dt_atmos:',   rc=rc)
    call ESMF_ConfigGetAttribute(config=CF, value=nfhmax,   label ='nhours_fcst:',rc=rc)
    if(mype == 0) print *,'af nems config,dt_atmos=',dt_atmos,'nfhmax=',nfhmax

    call ESMF_TimeIntervalSet(timeStep, s=dt_atmos, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    first_kdt = 1
    if( output_1st_tstep_rst) then
      rsthour   = currTime - StartTime
      first_kdt = nint(rsthour/timeStep) + 1
    endif

!
!#######################################################################
! set up fcst grid component
!
!----------------------------------------------------------------------
!*** create fv3 atm tasks and quilt servers
!-----------------------------------------------------------------------
!
! create fcst grid component

    if( quilting ) then
      num_pes_fcst = petcount - write_groups * wrttasks_per_group
    else
      num_pes_fcst = petcount
    endif
    allocate(fcstPetList(num_pes_fcst))
    do j=1, num_pes_fcst
      fcstPetList(j) = j - 1
    enddo
    fcstComp = ESMF_GridCompCreate(petList=fcstPetList, name='fv3_fcst', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
    ! copy attributes from fv3cap component to fcstComp
    call ESMF_InfoGetFromHost(gcomp, info=parentInfo, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_InfoGetFromHost(fcstComp, info=childInfo, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_InfoUpdate(lhs=childInfo, rhs=parentInfo, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! use the generic SetVM method to do resource and threading control
    call ESMF_GridCompSetVM(fcstComp, SetVM, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
    call ESMF_GridCompSetServices(fcstComp, fcstSS, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

! obtain fcst VM
    call ESMF_GridCompGet(fcstComp, vm=fcstVM, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
! create fcst state
    fcstState = ESMF_StateCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! call fcst Initialize (including creating fcstgrid and fcst fieldbundle)
    call ESMF_GridCompInitialize(fcstComp, exportState=fcstState,    &
                                 clock=clock, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
!
! reconcile the fcstComp's export state
    call ESMF_StateReconcile(fcstState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
! determine number elements in fcstState
    call ESMF_StateGet(fcstState, itemCount=FBCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if(mype == 0) print *,'af fcstCom FBCount= ',FBcount
!
! set start time for output
      output_startfh = 0.
!
!-----------------------------------------------------------------------
!***  create and initialize Write component(s).
!-----------------------------------------------------------------------
!
    if( quilting ) then

      allocate(fcstFB(FBCount), fcstItemNameList(FBCount), fcstItemTypeList(FBCount))
      allocate(wrtComp(write_groups), wrtState(write_groups) )
      allocate(wrtFB(FBCount,write_groups), routehandle(FBCount,write_groups))
      allocate(lead_wrttask(write_groups), last_wrttask(write_groups))
      allocate(petList(wrttasks_per_group))
      allocate(originPetList(num_pes_fcst+wrttasks_per_group))
      allocate(targetPetList(num_pes_fcst+wrttasks_per_group))
      if(mype == 0) print *,'af allco wrtComp,write_groups=',write_groups

! pull out the item names and item types from fcstState
      call ESMF_StateGet(fcstState, itemNameList=fcstItemNameList, &
                         itemTypeList=fcstItemTypeList, &
                        !itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                         rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! loop over all items in the fcstState and collect all FieldBundles
      do i=1, FBcount
        if (fcstItemTypeList(i) == ESMF_STATEITEM_FIELDBUNDLE) then
          ! access the FieldBundle
          call ESMF_StateGet(fcstState, itemName=fcstItemNameList(i), &
                             fieldbundle=fcstFB(i), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!          if(mype==0.or.mype==144) print *,'af fcstFB,i=',i,'name=',trim(fcstItemNameList(i))
        else
        !***### anything but a FieldBundle in the state is unexpected here
          call ESMF_LogSetError(ESMF_RC_ARG_BAD,                                 &
                                msg="Only FieldBundles supported in fcstState.", &
                                line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
        endif
      enddo
!
      k = num_pes_fcst
      timerhs = MPI_Wtime()
      do i=1, write_groups

! prepare petList for wrtComp(i)
        lead_wrttask(i) = k
        do j=1, wrttasks_per_group
          petList(j) = k + j-1
        enddo
        k = k + wrttasks_per_group
        last_wrttask(i) = k - 1
!        if(mype==0)print *,'af wrtComp(i)=',i,'k=',k

! prepare name of the wrtComp(i)
        write(cwrtcomp,"(A,I2.2)") "wrtComp_", i
! create wrtComp(i)
        wrtComp(i) = ESMF_GridCompCreate(petList=petList, name=trim(cwrtcomp), rc=rc)
!      print *,'af wrtComp(i)=',i,'name=',trim(cwrtcomp),'rc=',rc
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! copy attributes from fv3cap component to wrtComp
        call ESMF_InfoGetFromHost(wrtComp(i), info=childInfo, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        call ESMF_InfoUpdate(lhs=childInfo, rhs=parentInfo, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! use the generic SetVM method to do resource and threading control
        call ESMF_GridCompSetVM(wrtComp(i), SetVM, userRc=urc, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

! call into wrtComp(i) SetServices
        call ESMF_GridCompSetServices(wrtComp(i), wrtSS, userRc=urc, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

! add configuration file
        call ESMF_GridCompSet(gridcomp=wrtComp(i),config=CF,rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! create wrtstate(i)
        wrtstate(i) = ESMF_StateCreate(rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! add the fcst FieldBundles to the wrtState(i) so write component can
! use this info to create mirror objects
        call ESMF_AttributeCopy(fcstState, wrtState(i), &
                                attcopy=ESMF_ATTCOPY_REFERENCE, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_StateAdd(wrtState(i), fcstFB, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! call into wrtComp(i) Initialize
        call ESMF_GridCompInitialize(wrtComp(i), importState=wrtstate(i), &
                                     clock=clock, phase=1, userRc=urc, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

! remove fcst FieldBundles from the wrtState(i) because done with it
        call ESMF_StateRemove(wrtState(i), fcstItemNameList, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! reconcile the wrtComp(i)'s export state
        call ESMF_StateReconcile(wrtState(i), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if(mype==0) print *,'af wrtState reconcile, FBcount=',FBcount

        call ESMF_AttributeCopy(fcstState, wrtState(i), &
                                attcopy=ESMF_ATTCOPY_REFERENCE, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! loop over all FieldBundle in the states and precompute Regrid operation
        do j=1, FBcount

          ! access the mirrored FieldBundle in the wrtState(i)
          call ESMF_StateGet(wrtState(i),                                   &
                             itemName="mirror_"//trim(fcstItemNameList(j)), &
                             fieldbundle=wrtFB(j,i), rc=rc)
          if(mype == 0) print *,'af get wrtfb=',"mirror_"//trim(fcstItemNameList(j)),' rc=',rc
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! determine regridmethod
          if (index(fcstItemNameList(j),"_bilinear") >0 )  then
            regridmethod = ESMF_REGRIDMETHOD_BILINEAR
          else if (index(fcstItemNameList(j),"_patch") >0)  then
            regridmethod = ESMF_REGRIDMETHOD_PATCH
          else if (index(fcstItemNameList(j),"_nearest_stod") >0) then
            regridmethod = ESMF_REGRIDMETHOD_NEAREST_STOD
          else if (index(fcstItemNameList(j),"_nearest_dtos") >0) then
            regridmethod = ESMF_REGRIDMETHOD_NEAREST_DTOS
          else if (index(fcstItemNameList(j),"_conserve") >0) then
            regridmethod = ESMF_REGRIDMETHOD_CONSERVE
          else
            call ESMF_LogSetError(ESMF_RC_ARG_BAD,                          &
                                  msg="Unable to determine regrid method.", &
                                  line=__LINE__, file=__FILE__, rcToReturn=rc)
            return
          endif

          call ESMF_LogWrite('bf FieldBundleRegridStore', ESMF_LOGMSG_INFO, rc=rc)
          write(msgString,"(A,I2.2,',',I2.2,A)") "calling into wrtFB(",j,i, ") FieldBundleRegridStore()...."
          call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)

          if (i==1) then
! this is a Store() for the first wrtComp -> must do the Store()
            call ESMF_FieldBundleRegridStore(fcstFB(j), wrtFB(j,1),                                    &
                                             regridMethod=regridmethod, routehandle=routehandle(j,1),  &
                                             unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,                &
                                             srcTermProcessing=isrcTermProcessing, rc=rc)

!           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            if (rc /= ESMF_SUCCESS) then
              write(0,*)'fv3_cap.F90:InitializeAdvertise error in ESMF_FieldBundleRegridStore'
              call ESMF_LogWrite('fv3_cap.F90:InitializeAdvertise error in ESMF_FieldBundleRegridStore', ESMF_LOGMSG_ERROR, rc=rc)
              call ESMF_Finalize(endflag=ESMF_END_ABORT)
            endif
            call ESMF_LogWrite('af FieldBundleRegridStore', ESMF_LOGMSG_INFO, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            originPetList(1:num_pes_fcst)  = fcstPetList(:)
            originPetList(num_pes_fcst+1:) = petList(:)

          else
            targetPetList(1:num_pes_fcst)  = fcstPetList(:)
            targetPetList(num_pes_fcst+1:) = petList(:)
            routehandle(j,i) = ESMF_RouteHandleCreate(routehandle(j,1),            &
                                                      originPetList=originPetList, &
                                                      targetPetList=targetPetList, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          endif
          write(msgString,"(A,I2.2,',',I2.2,A)") "... returned from wrtFB(",j,i, ") FieldBundleRegridStore()."
          call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
        enddo  ! j=1, FBcount

! end write_groups
      enddo   ! i=1, write_groups
      if(mype==0) print *,'in fv3cap init, time wrtcrt/regrdst',MPI_Wtime()-timerhs
      deallocate(petList)
      deallocate(originPetList)
      deallocate(targetPetList)
!
!---------------------------------------------------------------------------------
!---  set up output forecast time array
!
!--- get current forecast length
      if(iau_offset > 0) then
        output_startfh = iau_offset
      endif
      if(mype==0) print *,'in fv3 cap init, output_startfh=',output_startfh,'nsout=',nsout, &
        'iau_offset=',iau_offset,'nfhmax_hf=',nfhmax_hf,'nfhout_hf=',nfhout_hf, &
        'nfhout=',nfhout
!
!--- set up output_fh with output forecast hours
! if the run does not have iau, it will have output after first step integration as fh00
! if the run has iau, it will start output at fh=00 at the cycle time (usually StartTime+IAU_offsetTI)
      if(nsout > 0) then
!--- use nsout for output frequency nsout*dt_atmos
        nfh = 0
        if( nfhmax > output_startfh ) nfh = nint((nfhmax-output_startfh)/(nsout*dt_atmos/3600.))+1
        if(nfh >0) then
          allocate(output_fh(nfh))
          if( output_startfh == 0) then
            output_fh(1) = dt_atmos/3600.
          else
            output_fh(1) = output_startfh
          endif
          do i=2,nfh
            output_fh(i) = (i-1)*nsout*dt_atmos/3600. + output_startfh
          enddo
        endif
      elseif (nfhmax_hf > 0 ) then
!--- use high frequency output and low frequency for output forecast time
        nfh = 0
        if( nfhout_hf>0 .and. nfhmax_hf>output_startfh) nfh = nint((nfhmax_hf-output_startfh)/nfhout_hf)+1
        nfh2 = 0
        if( nfhout>0 .and. nfhmax>nfhmax_hf) nfh2 = nint((nfhmax-nfhmax_hf)/nfhout)
        if( nfh+nfh2 > 0) then
          allocate(output_fh(nfh+nfh2))
          if( output_startfh == 0) then
            output_fh(1) = dt_atmos/3600.
          else
            output_fh(1) = output_startfh
          endif
          do i=2,nfh
            output_fh(i) = (i-1)*nfhout_hf + output_startfh
          enddo
          do i=1,nfh2
            output_fh(nfh+i) = nfhmax_hf + i*nfhout
          enddo
        endif
      elseif (nfhout > 0 ) then
!--- use one output freqency
        nfh = 0
        if( nfhout > 0 .and. nfhmax>output_startfh) nfh = nint((nfhmax-output_startfh)/nfhout) + 1
        if( nfh > 0 ) then
          allocate(output_fh(nfh))
          if( output_startfh == 0) then
            output_fh(1) = dt_atmos/3600.
          else
            output_fh(1) = output_startfh
          endif
          do i=2,nfh
            output_fh(i) = (i-1)*nfhout + output_startfh
          enddo
        endif
      endif
!
!-----------------------------------------------------------------------
!***  SET THE FIRST WRITE GROUP AS THE FIRST ONE TO ACT.
!-----------------------------------------------------------------------
!
      n_group = 1
!
!end quilting
    endif
!
!-- set up output forecast time if output_fh is specified
    if (noutput_fh > 0 ) then
!--- use output_fh to sepcify output forecast time
      loutput_fh = .true.
      lflname_fulltime = .false.
      if(noutput_fh == 1) then
        call ESMF_ConfigGetAttribute(CF,value=outputfh,label='output_fh:', rc=rc)
        if(outputfh == -1) loutput_fh = .false.
      endif
      if( loutput_fh ) then
        lfreq = .false.
        if( allocated(output_fh)) deallocate(output_fh)
        if(noutput_fh == 2) then
          call ESMF_ConfigGetAttribute(CF,valueList=outputfh2,label='output_fh:', &
             count=noutput_fh, rc=rc)
          if(outputfh2(2) == -1) then
            !--- output_hf is output frequency, the second item is -1
            lfreq = .true.
            nfh = 0
            if( nfhmax>output_startfh) nfh = nint((nfhmax-output_startfh)/outputfh2(1)) + 1
            if( nfh > 0) then
              allocate(output_fh(nfh))
              if( output_startfh == 0) then
                output_fh(1) = dt_atmos/3600.
              else
                output_fh(1) = output_startfh
              endif
              do i=2,nfh
                output_fh(i) = (i-1)*outputfh2(1) + output_startfh
                ! Except fh000, which is the first time output, if any other of the
                ! output time is not integer hour, set lflname_fulltime to be true, so the
                ! history file names will contain the full time stamp (HHH-MM-SS).
                if(.not.lflname_fulltime) then
                  if(mod(nint(output_fh(i)*3600.),3600) /= 0) lflname_fulltime = .true.
                endif
              enddo
            endif
          endif
        endif
        if( noutput_fh /= 2 .or. .not. lfreq ) then
          allocate(output_fh(noutput_fh))
          output_fh = 0
          call ESMF_ConfigGetAttribute(CF,valueList=output_fh,label='output_fh:', &
             count=noutput_fh, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          if( output_startfh == 0) then
            ! If the output time in output_fh array contains first time stamp output,
            ! check the rest of output time, otherwise, check all the output time.
            ! If any of them is not integer hour, the history file names will
            ! contain the full time stamp (HHH-MM-SS)
            ist = 1
            if(output_fh(1)==0) then
              output_fh(1) = dt_atmos/3600.
              ist= 2
            endif
            do i=ist,noutput_fh
              if(.not.lflname_fulltime) then
                if(mod(nint(output_fh(i)*3600.),3600) /= 0) lflname_fulltime = .true.
              endif
            enddo
          else
            do i=1,noutput_fh
              output_fh(i) = output_startfh + output_fh(i)
              ! When output_startfh >0, check all the output time, if any of
              ! them is not integer hour, set lflname_fulltime to be true. The
              ! history file names will contain the full time stamp (HHH-MM-SS).
              if(.not.lflname_fulltime) then
                if(mod(nint(output_fh(i)*3600.),3600) /= 0) lflname_fulltime = .true.
              endif
            enddo
          endif
        endif
      endif ! end loutput_fh
    endif
    if(mype==0) print *,'output_fh=',output_fh(1:size(output_fh)),'lflname_fulltime=',lflname_fulltime
!
    ! --- advertise Fields in importState and exportState -------------------

    ! importable fields:
    do i = 1, size(importFieldsInfo)
      call NUOPC_Advertise(importState, &
                           StandardName=trim(importFieldsInfo(i)%name), &
                           SharePolicyField='share', vm=fcstVM, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    end do

    ! exportable fields:
    do i = 1, size(exportFieldsInfo)
      call NUOPC_Advertise(exportState, &
                           StandardName=trim(exportFieldsInfo(i)%name), &
                           SharePolicyField='share', vm=fcstVM, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    end do

    if(mype==0) print *,'in fv3_cap, aft import, export fields in atmos'
    if(mype==0) print *,'in fv3_cap, init time=',MPI_Wtime()-timeis
!-----------------------------------------------------------------------
!
  end subroutine InitializeAdvertise

!-----------------------------------------------------------------------------

  subroutine InitializeRealize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*),parameter :: subname='(fv3gfs_cap:InitializeRealize)'
    type(ESMF_State)           :: importState, exportState
    logical                    :: isPetLocal

    rc = ESMF_SUCCESS

    ! query for importState and exportState
    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! --- conditionally realize or remove Fields in importState and exportState -------------------

    isPetLocal = ESMF_GridCompIsPetLocal(fcstComp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,  file=__FILE__)) return

    if (isPetLocal) then

      ! -- realize connected fields in exportState
      call realizeConnectedCplFields(exportState, fcstGrid(mygrid),                  &
                                     numLevels, numSoilLayers, numTracers,           &
                                     exportFieldsInfo, 'FV3 Export', exportFields, 0.0_ESMF_KIND_R8, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,  file=__FILE__)) return

      ! -- initialize export fields if applicable
      call setup_exportdata(rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,  file=__FILE__)) return

      ! -- realize connected fields in importState
      call realizeConnectedCplFields(importState, fcstGrid(mygrid),                  &
                                     numLevels, numSoilLayers, numTracers,           &
                                     importFieldsInfo, 'FV3 Import', importFields, 9.99e20_ESMF_KIND_R8, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,  file=__FILE__)) return

    end if

  end subroutine InitializeRealize

!-----------------------------------------------------------------------------

  subroutine ModelAdvance(gcomp, rc)

    type(ESMF_GridComp)         :: gcomp
    integer, intent(out)        :: rc

!-----------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (profile_memory) call ESMF_VMLogMemInfo("Entering FV3 ModelAdvance: ")

    call ModelAdvance_phase1(gcomp, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ModelAdvance_phase2(gcomp, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (profile_memory) call ESMF_VMLogMemInfo("Leaving FV3 ModelAdvance: ")

  end subroutine ModelAdvance

!-----------------------------------------------------------------------------

  subroutine ModelAdvance_phase1(gcomp, rc)
    type(ESMF_GridComp)         :: gcomp
    integer, intent(out)        :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    integer                     :: urc
    logical                     :: fcstpe
    character(len=*),parameter  :: subname='(fv3_cap:ModelAdvance_phase1)'
    character(240)              :: msgString

!-----------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if(profile_memory) call ESMF_VMLogMemInfo("Entering FV3 ModelAdvance_phase1: ")

    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_ClockPrint(clock, options="currTime", &
                         preString="entering FV3_ADVANCE phase1 with clock current: ", &
                         unit=msgString)
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
    call ESMF_ClockPrint(clock, options="startTime", &
                         preString="entering FV3_ADVANCE phase1 with clock start:   ", &
                         unit=msgString)
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
    call ESMF_ClockPrint(clock, options="stopTime", &
                         preString="entering FV3_ADVANCE phase1 with clock stop:    ", &
                         unit=msgString)
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

    call ESMF_GridCompRun(fcstComp, exportState=fcstState, clock=clock, phase=1, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    if( dbug > 0 .or. cplprint_flag ) then
         fcstpe = .false.
         if( mype < num_pes_fcst ) fcstpe = .true.
         call diagnose_cplFields(gcomp, clock, fcstpe, cplprint_flag, dbug, 'import')
    endif

    if (profile_memory) call ESMF_VMLogMemInfo("Leaving FV3 ModelAdvance_phase1: ")

  end subroutine ModelAdvance_phase1

!-----------------------------------------------------------------------------

  subroutine ModelAdvance_phase2(gcomp, rc)
    type(ESMF_GridComp)         :: gcomp
    integer, intent(out)        :: rc

    ! local variables
    type(ESMF_Time)             :: currTime
    type(ESMF_TimeInterval)     :: timeStep
    type(ESMF_Time)             :: startTime
    type(ESMF_TimeInterval)     :: time_elapsed

    integer                     :: na, j, urc
    integer                     :: nfseconds
    logical                     :: fcstpe
    character(len=*),parameter  :: subname='(fv3_cap:ModelAdvance_phase2)'

    character(240)              :: msgString

    type(ESMF_Clock)            :: clock, clock_out

!-----------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if(profile_memory) call ESMF_VMLogMemInfo("Entering FV3 ModelAdvance_phase2: ")

    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompRun(fcstComp, exportState=fcstState, clock=clock, phase=2, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    clock_out = ESMF_ClockCreate(clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_ClockAdvance(clock_out, rc = RC)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!-------------------------------------------------------------------------------
!*** if it is output time, call data transfer and write grid comp run
    if( quilting ) then

      call ESMF_ClockGet(clock_out, startTime=startTime, currTime=currTime, &
                         timeStep=timeStep, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      time_elapsed  = currTime - startTime
      na = nint(time_elapsed/timeStep)
      call ESMF_TimeIntervalGet(time_elapsed, s=nfseconds, rc=rc)

      output: if (ANY(nint(output_fh(:)*3600.0) == nfseconds)) then
!
       if (mype == 0 .or. mype == lead_wrttask(1)) print *,' aft fcst run output time=',nfseconds, &
       'FBcount=',FBcount,'na=',na

        call ESMF_VMEpochEnter(epoch=ESMF_VMEpoch_Buffer, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        do j=1, FBCount

          call ESMF_FieldBundleRegrid(fcstFB(j), wrtFB(j,n_group),         &
                                      routehandle=routehandle(j, n_group), &
                                      termorderflag=(/ESMF_TERMORDER_SRCSEQ/), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
        enddo

        call ESMF_VMEpochExit(rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_LogWrite('Model Advance: before wrtcomp run ', ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_GridCompRun(wrtComp(n_group), importState=wrtState(n_group), clock=clock_out, userRc=urc, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call ESMF_LogWrite('Model Advance: after wrtcomp run ', ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if (n_group == write_groups) then
          n_group = 1
        else
          n_group = n_group + 1
        endif

      endif output

    endif ! quilting

    call ESMF_ClockPrint(clock, options="currTime", &
                         preString="leaving FV3_ADVANCE phase2 with clock current: ", &
                         unit=msgString)
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
    call ESMF_ClockPrint(clock, options="startTime", &
                         preString="leaving FV3_ADVANCE phase2 with clock start:   ", &
                         unit=msgString)
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
    call ESMF_ClockPrint(clock, options="stopTime", &
                         preString="leaving FV3_ADVANCE phase2 with clock stop:    ", &
                         unit=msgString)
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

    if( dbug > 0 .or. cplprint_flag ) then
      fcstpe = .false.
      if( mype < num_pes_fcst ) fcstpe = .true.
      call diagnose_cplFields(gcomp, clock_out, fcstpe, cplprint_flag, dbug, 'export')
    end if

    if(profile_memory) call ESMF_VMLogMemInfo("Leaving FV3 ModelAdvance_phase2: ")

  end subroutine ModelAdvance_phase2

!-----------------------------------------------------------------------------

  subroutine ModelSetRunClock(gcomp, rc)

    type(ESMF_GridComp)         :: gcomp
    integer, intent(out)        :: rc

    ! local variables
    type(ESMF_Clock)            :: dclock, mclock
    type(ESMF_TimeInterval)     :: dtimestep, mtimestep
    type(ESMF_Time)             :: mcurrtime, mstoptime

!-----------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_ClockGet(dclock, timeStep=dtimestep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_ClockGet(mclock, currTime=mcurrtime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_TimeIntervalSet(mtimestep,s=dt_atmos,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    mstoptime = mcurrtime + dtimestep

    call ESMF_ClockSet(mclock, timeStep=mtimestep, stopTime=mstoptime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine ModelSetRunClock

!-----------------------------------------------------------------------------

  subroutine fv3_checkimport(gcomp, rc)

!***  Check the import state fields

    ! input arguments
    type(ESMF_GridComp)        :: gcomp
    integer, intent(out)       :: rc

    ! local variables
    character(len=*),parameter :: subname='(fv3gfs_cap:fv3_checkimport)'
    integer                    :: n, nf
    type(ESMF_Clock)           :: clock
    type(ESMF_Time)            :: currTime, invalidTime
    type(ESMF_State)           :: importState
    logical                    :: isValid
    type(ESMF_Field),pointer   :: fieldList(:)
    character(len=128)         :: fldname
    character(esmf_maxstr)     :: msgString
    integer                    :: date(6)

    rc = ESMF_SUCCESS

    ! query the Component for its clock
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! get the current time out of the clock
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    date(1:6) = 0
    call ESMF_TimeGet(time=currTime,yy=date(1),mm=date(2),dd=date(3),h=date(4), &
                      m=date(5),s=date(6),rc=rc)
!   if(mype==0) print *,'in fv3_checkimport, currtime=',date(1:6)

    ! set up invalid time (by convention)
    call ESMF_TimeSet(invalidTime, yy=99999999, mm=01, dd=01, &
                      h=00, m=00, s=00, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    nullify(fieldList)
    call NUOPC_GetStateMemberLists(importState, fieldList=fieldList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! set the importFieldsValid flag
    ! associated(fieldList) will be false if there are no fields

    importFieldsValid(:) = .true.
    if (associated(fieldList)) then
!     if(mype==0) print *,'in fv3_checkimport, inside associated(fieldList)'
      do n = 1,size(fieldList)
        call ESMF_FieldGet(fieldList(n), name=fldname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        ! check if import field carries a valid timestamp
        call NUOPC_GetTimestamp(fieldList(n), isValid=isValid, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if (isValid) then
          ! if timestamp is set, check if it is valid
          isValid = .not.NUOPC_IsAtTime(fieldList(n), invalidTime, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        end if

        ! store field status in internal array
        nf = queryImportFields(fldname)
        importFieldsValid(nf) = isValid

        if (isValid) then
          ! check if field is current
          isValid = NUOPC_IsAtTime(fieldList(n), currTime, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          if (.not.isValid) then
            call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
                                  msg="NUOPC INCOMPATIBILITY DETECTED: Import Field " &
                                      // trim(fldname) // " not at current time", &
                                  line=__LINE__, file=__FILE__, rcToReturn=rc)
            return
          end if
        end if
        write(msgString,'(A,2i4,l3)') "fv3_checkimport "//trim(fldname),n,nf,importFieldsValid(nf)
        call ESMF_LogWrite(msgString,ESMF_LOGMSG_INFO,rc=rc)
      enddo
    endif

  end subroutine fv3_checkimport

!-----------------------------------------------------------------------------

  subroutine TimestampExport_phase1(gcomp, rc)

    ! input arguments
    type(ESMF_GridComp)        :: gcomp
    integer, intent(out)       :: rc

    ! local variables
    character(len=*),parameter :: subname='(fv3gfs_cap:TimestampExport_phase1)'
    type(ESMF_Clock)           :: driverClock, modelClock
    type(ESMF_State)           :: exportState

    rc = ESMF_SUCCESS

    ! get driver and model clock
    call NUOPC_ModelGet(gcomp, driverClock=driverClock, &
                        modelClock=modelClock, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! reset model clock to initial time
    call NUOPC_CheckSetClock(modelClock, driverClock, forceCurrTime=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! update timestamp on export Fields
    call NUOPC_SetTimestamp(exportState, modelClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine TimestampExport_phase1

!-----------------------------------------------------------------------------

  subroutine ModelFinalize(gcomp, rc)

    ! input arguments
    type(ESMF_GridComp)        :: gcomp
    integer, intent(out)       :: rc

    ! local variables
    character(len=*),parameter :: subname='(fv3gfs_cap:ModelFinalize)'
    integer                    :: i, urc
    type(ESMF_VM)              :: vm
    real(kind=8)               :: MPI_Wtime, timeffs
!
!-----------------------------------------------------------------------------
!*** finialize forecast

    timeffs = MPI_Wtime()
    rc = ESMF_SUCCESS
!
    call ESMF_GridCompGet(gcomp,vm=vm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
!*** finalize grid comps
    if( quilting ) then
      do i = 1, write_groups
        call ESMF_GridCompFinalize(wrtComp(i), importState=wrtstate(i),userRc=urc, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      enddo
    endif

    call ESMF_GridCompFinalize(fcstComp, exportState=fcststate,userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
!
!*** destroy grid comps
    if( quilting ) then
      do i = 1, write_groups
        call ESMF_StateDestroy(wrtState(i), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        call ESMF_GridCompDestroy(wrtComp(i), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      enddo
    endif

    call ESMF_StateDestroy(fcstState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_GridCompDestroy(fcstComp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
    if(mype==0)print *,' wrt grid comp destroy time=',MPI_Wtime()-timeffs

  end subroutine ModelFinalize
!
!-----------------------------------------------------------------------------

end module fv3gfs_cap_mod
