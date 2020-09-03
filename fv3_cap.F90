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
  use NUOPC_Model,            only: model_routine_SS        => SetServices,       &
                                    model_routine_Run       => routine_Run,       &
                                    model_label_Advance     => label_Advance,     &
                                    model_label_CheckImport => label_CheckImport, &
                                    model_label_Finalize    => label_Finalize,    &
                                    NUOPC_ModelGet
!
  use module_fv3_config,      only: quilting, restart_interval,              &
                                    nfhout, nfhout_hf, nsout, dt_atmos,      &
                                    nfhmax, nfhmax_hf,output_hfmax,          &
                                    output_interval,output_interval_hf,      &
                                    alarm_output_hf, alarm_output,           &
                                    calendar, calendar_type, cpl,            &
                                    force_date_from_configure,               &
                                    cplprint_flag,output_1st_tstep_rst,      &
                                    first_kdt,num_restart_interval       

  use module_fv3_io_def,      only: num_pes_fcst,write_groups,app_domain,    &
                                    num_files, filename_base,                &
                                    wrttasks_per_group, n_group,             &
                                    lead_wrttask, last_wrttask,              &
                                    output_grid, output_file,                &
                                    imo,jmo,ichunk2d,jchunk2d,write_nemsioflip,&
                                    ichunk3d,jchunk3d,kchunk3d,              &
                                    write_fsyncflag, nsout_io,               &
                                    cen_lon, cen_lat, ideflate,              &
                                    lon1, lat1, lon2, lat2, dlon, dlat,      &
                                    stdlat1, stdlat2, dx, dy, iau_offset,    &
                                    nbits
!
  use module_fcst_grid_comp,  only: fcstSS => SetServices,                   &
                                    fcstGrid, numLevels, numSoilLayers,      &
                                    numTracers, num_diag_sfc_emis_flux,      &
                                    num_diag_down_flux,                      &
                                    num_diag_type_down_flux,                 &
                                    num_diag_burn_emis_flux, num_diag_cmass

  use module_wrt_grid_comp,   only: wrtSS => SetServices
!
  use module_cplfields,       only: nExportFields,    exportFields,          &
                                    exportFieldsList, exportFieldTypes,      &
                                    exportFieldShare,                        &
                                    nImportFields,    importFields,          &
                                    importFieldsList, importFieldTypes,      &
                                    importFieldShare, importFieldsValid,     &
                                    queryFieldList,   fillExportFields,      &
                                    exportData
  use module_cap_cpl,         only: realizeConnectedCplFields,               &
                                    clock_cplIntval, diagnose_cplFields      


  implicit none
  private
  public SetServices
!
!-----------------------------------------------------------------------
!
  type(ESMF_Clock),save                       :: clock_fv3

  type(ESMF_GridComp)                         :: fcstComp
  type(ESMF_State)                            :: fcstState
  character(len=esmf_maxstr),allocatable      :: fcstItemNameList(:)
  type(ESMF_StateItem_Flag), allocatable      :: fcstItemTypeList(:)
  type(ESMF_FieldBundle),    allocatable      :: fcstFB(:)
  integer, save                               :: FBCount

  type(ESMF_GridComp),    allocatable         :: wrtComp(:)
  type(ESMF_State),       allocatable         :: wrtState(:)
  type(ESMF_FieldBundle), allocatable         :: wrtFB(:,:)

  type(ESMF_RouteHandle), allocatable         :: routehandle(:,:)
  integer, allocatable                        :: fcstPetList(:)

  logical :: profile_memory = .true.

  character(len=160) :: nuopcMsg
  integer            :: timeslice = 0
  integer            :: fcstmype
  integer            :: dbug = 0

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

    ! initialization, switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
                                    userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
                                 phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
                                 phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! model advance method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
                              specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! checking the import fields is a bit more complex because of coldstart option
    call ESMF_MethodRemove(gcomp, model_label_CheckImport, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_CheckImport, &
                              specRoutine=fv3_checkimport, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! setup Run/Advance phase: phase1
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
                                 phaseLabelList=(/"phase1"/), userRoutine=model_routine_Run, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
                              specPhaseLabel="phase1", specRoutine=ModelAdvance_phase1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! setup Run/Advance phase: phase2
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
                                 phaseLabelList=(/"phase2"/), userRoutine=model_routine_Run, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
                              specPhaseLabel="phase2", specRoutine=ModelAdvance_phase2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! specializations required to support 'inline' run sequences
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_CheckImport, &
                              specPhaseLabel="phase1", specRoutine=NUOPC_NoOp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(gcomp, specLabel="ModelBase_TimestampExport", &
                              specPhaseLabel="phase1", specRoutine=TimestampExport_phase1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_CheckImport, &
                              specPhaseLabel="phase2", specRoutine=NUOPC_NoOp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! model finalize method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
                              specRoutine=atmos_model_finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine SetServices

  !-----------------------------------------------------------------------------

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc

    character(len=10)     :: value
    character(240)        :: msgString
    logical               :: isPresent, isSet
    character(len=*),parameter  :: subname='(fv3gfs_cap:InitializeP0)'

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
                                  acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeGet(gcomp, name="ProfileMemory", value=value, defaultValue="true", &
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

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    type(ESMF_GridComp)                    :: gcomp
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Clock)                       :: clock
    integer, intent(out)                   :: rc

!
! local variables
    type(ESMF_VM)                          :: vm
    type(ESMF_Time)                        :: CurrTime, starttime, StopTime
    type(ESMF_Time)                        :: alarm_output_hf_ring, alarm_output_ring
    type(ESMF_Time)                        :: alarm_output_hf_stop, alarm_output_stop
    type(ESMF_TimeInterval)                :: RunDuration, timeStep, rsthour, IAU_offsetTI
    type(ESMF_Config)                      :: cf
    type(ESMF_RegridMethod_Flag)           :: regridmethod
    type(ESMF_TimeInterval)                :: earthStep
    integer(ESMF_KIND_I4)                  :: nhf, nrg

    integer,dimension(6)                   :: date, date_init
    integer                                :: mpi_comm_atm
    integer                                :: i, j, k, io_unit, urc, ierr
    integer                                :: petcount, mype
    integer                                :: num_output_file
    logical                                :: isPetLocal
    logical                                :: OPENED
    character(ESMF_MAXSTR)                 :: name
    logical                                :: fcstpe
    logical,dimension(:), allocatable      :: wrtpe
    integer,dimension(:), allocatable      :: petList, originPetList, targetPetList
    character(20)                          :: cwrtcomp
    character(160)                         :: msg
    integer                                :: isrcTermProcessing

    character(len=*),parameter  :: subname='(fv3_cap:InitializeAdvertise)'
    integer nfmout, nfsout , nfmout_hf, nfsout_hf
    real(kind=8) :: MPI_Wtime, timewri, timeis,timeie,timerhs, timerhe
!
!------------------------------------------------------------------------
!
    rc = ESMF_SUCCESS
    timeis = mpi_wtime()

    call ESMF_GridCompGet(gcomp,name=name,vm=vm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_VMGet(vm, mpiCommunicator=mpi_comm_atm,petCount=petcount, &
                    localpet = mype,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    fcstmype = mype
!    print *,'in fv3_cap,initAdvertize,name=',trim(name),'mpi_comm=',mpi_comm_atm, &
!       'petcount=',petcount,'mype=',mype
!
! create an instance clock for fv3
    clock_fv3 = ESMF_ClockCreate(clock, rc=RC)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
!------------------------------------------------------------------------
! get config variables
!
    CF = ESMF_ConfigCreate(rc=RC)
    CALL ESMF_ConfigLoadFile(config=CF ,filename='model_configure' ,rc=RC)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
    num_restart_interval = ESMF_ConfigGetLen(config=CF, label ='restart_interval:',rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if(mype == 0) print *,'af nems config,num_restart_interval=',num_restart_interval
    if (num_restart_interval<=0) num_restart_interval = 1
    allocate(restart_interval(num_restart_interval))
    restart_interval = 0
    CALL  ESMF_ConfigGetAttribute(CF,valueList=restart_interval,label='restart_interval:', &
      count=num_restart_interval, rc=RC)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if(mype == 0) print *,'af nems config,restart_interval=',restart_interval
!
    CALL ESMF_ConfigGetAttribute(config=CF,value=calendar, &
                                 label ='calendar:', &
                                 default='gregorian',rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
    CALL ESMF_ConfigGetAttribute(config=CF,value=cpl,default=.false.,label ='cpl:',rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
    CALL ESMF_ConfigGetAttribute(config=CF,value=quilting, &
                                 label ='quilting:',rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    CALL ESMF_ConfigGetAttribute(config=CF,value=output_1st_tstep_rst, &
                                 default=.false., label ='output_1st_tstep_rst:',rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_ConfigGetAttribute(config=CF,value=iau_offset,default=0,label ='iau_offset:',rc=rc)
    if (iau_offset < 0) iau_offset=0

    ! chunksizes for netcdf_parallel
    call ESMF_ConfigGetAttribute(config=CF,value=ichunk2d,default=0,label ='ichunk2d:',rc=rc)
    call ESMF_ConfigGetAttribute(config=CF,value=jchunk2d,default=0,label ='jchunk2d:',rc=rc)
    call ESMF_ConfigGetAttribute(config=CF,value=ichunk3d,default=0,label ='ichunk3d:',rc=rc)
    call ESMF_ConfigGetAttribute(config=CF,value=jchunk3d,default=0,label ='jchunk3d:',rc=rc)
    call ESMF_ConfigGetAttribute(config=CF,value=kchunk3d,default=0,label ='kchunk3d:',rc=rc)
    
    ! zlib compression flag
    call ESMF_ConfigGetAttribute(config=CF,value=ideflate,default=0,label ='ideflate:',rc=rc)
    if (ideflate < 0) ideflate=0

    call ESMF_ConfigGetAttribute(config=CF,value=nbits,default=0,label ='nbits:',rc=rc)
    ! nbits quantization level for lossy compression (must be between 1 and 31)
    ! 1 is most compression, 31 is least. If outside this range, set to zero
    ! which means use lossless compression.
    if (nbits < 1 .or. nbits > 31)  nbits=0  ! lossless compression (no quantization)

    if(mype == 0) print *,'af nems config,quilting=',quilting,'calendar=', trim(calendar),' iau_offset=',iau_offset
    if(mype == 0) print *,'af nems config,ideflate=',ideflate,'nbits=',nbits
!
    nfhout = 0 ; nfhmax_hf = 0 ; nfhout_hf = 0 ; nsout = 0
    if ( quilting ) then
      CALL ESMF_ConfigGetAttribute(config=CF,value=write_groups, &
                                   label ='write_groups:',rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
      CALL ESMF_ConfigGetAttribute(config=CF,value=wrttasks_per_group, &
                                   label ='write_tasks_per_group:',rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      CALL ESMF_ConfigGetAttribute(config=CF,value=app_domain, default="global", &
                                   label ='app_domain:',rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      CALL ESMF_ConfigGetAttribute(config=CF,value=isrcTermProcessing, default=0, &
                                   label ='isrcTermProcessing:',rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if(mype == 0) print *,'af nems config,quilting=',quilting,'write_groups=', &
        write_groups,wrttasks_per_group,'calendar=',trim(calendar),'calendar_type=',calendar_type, &
        'isrcTermProcessing=', isrcTermProcessing
!
      CALL ESMF_ConfigGetAttribute(config=CF,value=num_files, &
                                   label ='num_files:',rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
      allocate(filename_base(num_files))
      CALL ESMF_ConfigFindLabel(CF,'filename_base:',rc=RC)
      do i=1,num_files
        CALL ESMF_ConfigGetAttribute(config=CF,value=filename_base(i), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      enddo

      allocate(output_file(num_files))
      num_output_file = ESMF_ConfigGetLen(config=CF, label ='output_file:',rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      if (num_files == num_output_file) then
        CALL ESMF_ConfigGetAttribute(CF,valueList=output_file,label='output_file:', &
             count=num_files, rc=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        do i = 1, num_files
          if(output_file(i) /= "netcdf" .and. output_file(i) /= "netcdf_parallel") then
            write(0,*)"fv3_cap.F90: only netcdf and netcdf_parallel are allowed for multiple values of output_file"
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
          endif
        enddo
      else if ( num_output_file == 1) then
        CALL ESMF_ConfigGetAttribute(CF,valuelist=output_file,label='output_file:', count=1, rc=RC)
        output_file(1:num_files) = output_file(1)
      else
        output_file(1:num_files) = 'netcdf'
      endif
      if(mype == 0) then
        print *,'af nems config,num_files=',num_files
        do i=1,num_files
           print *,'num_file=',i,'filename_base= ',trim(filename_base(i)),&
           ' output_file= ',trim(output_file(i))
        enddo
      endif
!
! variables for alarms
      call ESMF_ConfigGetAttribute(config=CF, value=nfhout,   label ='nfhout:',   rc=rc)
      call ESMF_ConfigGetAttribute(config=CF, value=nfhmax_hf,label ='nfhmax_hf:',rc=rc)
      call ESMF_ConfigGetAttribute(config=CF, value=nfhout_hf,label ='nfhout_hf:',rc=rc)
      call ESMF_ConfigGetAttribute(config=CF, value=nsout,    label ='nsout:',rc=rc)
      nsout_io = nsout
      if(mype==0) print *,'af nems config,nfhout,nsout=',nfhout,nfhmax_hf,nfhout_hf, nsout

! variables for I/O options
      call ESMF_ConfigGetAttribute(config=CF, value=output_grid, label ='output_grid:',rc=rc)
      if (mype == 0) then
        print *,'output_grid=',trim(output_grid)
      end if
      write_nemsioflip =.false.
      write_fsyncflag  =.false.

      if(trim(output_grid) == 'gaussian_grid' .or. trim(output_grid) == 'global_latlon') then
        call ESMF_ConfigGetAttribute(config=CF, value=imo, label ='imo:',rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=jmo, label ='jmo:',rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=write_nemsioflip, label ='write_nemsioflip:', default=.true., rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=write_fsyncflag,  label ='write_fsyncflag:', default=.true., rc=rc)
        if (mype == 0) then
          print *,'imo=',imo,'jmo=',jmo
          print *,'write_nemsioflip=',write_nemsioflip,'write_fsyncflag=',write_fsyncflag
        end if
      else if(trim(output_grid) == 'regional_latlon') then
        call ESMF_ConfigGetAttribute(config=CF, value=lon1, label ='lon1:',rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=lat1, label ='lat1:',rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=lon2, label ='lon2:',rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=lat2, label ='lat2:',rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=dlon, label ='dlon:',rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=dlat, label ='dlat:',rc=rc)
        imo = (lon2-lon1)/dlon + 1
        jmo = (lat2-lat1)/dlat + 1
        if (mype == 0) then
          print *,'lon1=',lon1,' lat1=',lat1
          print *,'lon2=',lon2,' lat2=',lat2
          print *,'dlon=',dlon,' dlat=',dlat
          print *,'imo =',imo, ' jmo=',jmo
        end if
      else if (trim(output_grid) == 'rotated_latlon') then
        call ESMF_ConfigGetAttribute(config=CF, value=cen_lon, label ='cen_lon:',rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=cen_lat, label ='cen_lat:',rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=lon1,    label ='lon1:',   rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=lat1,    label ='lat1:',   rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=lon2,    label ='lon2:',   rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=lat2,    label ='lat2:',   rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=dlon,    label ='dlon:',   rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=dlat,    label ='dlat:',   rc=rc)
        imo = (lon2-lon1)/dlon + 1
        jmo = (lat2-lat1)/dlat + 1
        if (mype == 0) then
          print *,'lon1=',lon1,' lat1=',lat1
          print *,'lon2=',lon2,' lat2=',lat2
          print *,'dlon=',dlon,' dlat=',dlat
          print *,'imo =',imo, ' jmo=',jmo
        end if
      else if (trim(output_grid) == 'lambert_conformal') then
        call ESMF_ConfigGetAttribute(config=CF, value=cen_lon, label ='cen_lon:',rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=cen_lat, label ='cen_lat:',rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=stdlat1, label ='stdlat1:',rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=stdlat2, label ='stdlat2:',rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=imo,     label ='nx:',     rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=jmo,     label ='ny:',     rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=lon1,    label ='lon1:',   rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=lat1,    label ='lat1:',   rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=dx,      label ='dx:',     rc=rc)
        call ESMF_ConfigGetAttribute(config=CF, value=dy,      label ='dy:',     rc=rc)
        if (mype == 0) then
          print *,'cen_lon=',cen_lon,' cen_lat=',cen_lat
          print *,'stdlat1=',stdlat1,' stdlat2=',stdlat2
          print *,'lon1=',lon1,' lat1=',lat1
          print *,'nx=',imo, ' ny=',jmo
          print *,'dx=',dx,' dy=',dy
        endif
      endif ! output_grid

    endif ! quilting
!
    call ESMF_ConfigGetAttribute(config=CF, value=dt_atmos, label ='dt_atmos:',   rc=rc)
    call ESMF_ConfigGetAttribute(config=CF, value=nfhmax,   label ='nhours_fcst:',rc=rc)
    if(mype == 0) print *,'af nems config,dt_atmos=',dt_atmos,'nfhmax=',nfhmax
    call ESMF_TimeIntervalSet(timeStep,s=dt_atmos,rc=rc)
    call ESMF_ClockSet(clock_fv3,timeStep=timeStep, rc=rc)
!
!------------------------------------------------------------------------
! may need to set currTime for restart
!
    call ESMF_ClockGet(clock_fv3, currTIME=CurrTime,  StartTime=startTime,    &
                       RunDuration=RunDuration, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    StopTime = startTime + RunDuration

! *** read restart time from restart file
    do i=751,899
       inquire(i, opened=OPENED)
       if(.not. OPENED)then
         io_unit = i
         exit
       endif
    enddo
!
    date = 0 ; date_init = 0
    force_date_from_configure = .true.
!
    open(unit=io_unit, file=trim('INPUT/coupler.res'),status="old",err=998 )
    read (io_unit,*,err=999) calendar_type
    read (io_unit,*) date_init
    read (io_unit,*) date
    close(io_unit)
    force_date_from_configure = .false.
!
    if(date(1) == 0 .and. date_init(1) /= 0) date = date_init
    if(mype == 0) print *,'bf clock_fv3,date=',date,'date_init=',date_init

    call ESMF_VMbroadcast(vm, date, 6, 0)
    call ESMF_TimeSet(time=CurrTime,yy=date(1),mm=date(2),dd=date(3),h=date(4), &
                      m=date(5),s=date(6),rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
999 continue
998 continue
!    if(mype==0) print *,'final date =',date,'date_init=',date_init

!reset CurrTime in clock
    call ESMF_ClockSet(clock_fv3, currTIME=CurrTime, startTime=startTime,  &
                       stopTime=stopTime, timeStep=timeStep, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
 
    !
    !Under NUOPC, the EARTH driver clock is a separate instance from the
    ! - fv3 clock. However, the fv3 clock may have been reset from restart
    ! - therefore the EARTH driver clock must also be adjusted.
    ! - Affected: currTime, timeStep
    call ESMF_ClockGet(clock, timeStep=earthStep, rc=RC)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    if (earthStep > (stopTime-currTime)) earthStep = stopTime - currTime
    call ESMF_ClockSet(clock, currTime=currTime, timeStep=earthStep, rc=RC)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Set fv3 component clock as copy of EARTH clock.
    call NUOPC_CompSetClock(gcomp, clock, rc=RC)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Read in the FV3 coupling interval
    if ( cpl ) then
      call clock_cplIntval(gcomp, CF)
    endif
!
    first_kdt = 1
    if( output_1st_tstep_rst) then
      rsthour   = CurrTime - StartTime
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

    fcstpe = .false.
    num_pes_fcst = petcount - write_groups * wrttasks_per_group
    allocate(fcstPetList(num_pes_fcst))
    do j=1, num_pes_fcst
      fcstPetList(j) = j - 1
      if(mype == fcstPetlist(j)) fcstpe = .true.
    enddo
    fcstComp = ESMF_GridCompCreate(petList=fcstPetList, name='fv3_fcst', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
    call ESMF_GridCompSetServices(fcstComp, fcstSS, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

! create fcst state
    fcstState = ESMF_StateCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! call fcst Initialize (including creating fcstgrid and fcst fieldbundle)
    call ESMF_GridCompInitialize(fcstComp, exportState=fcstState,    &
                                 clock=clock_fv3, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
!
! reconcile the fcstComp's import state
    call ESMF_StateReconcile(fcstState, attreconflag= ESMF_ATTRECONCILE_ON, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
! determine number elements in fcstState
    call ESMF_StateGet(fcstState, itemCount=FBCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if(mype == 0) print *,'af fcstCom FBCount= ',FBcount
!
! allocate arrays
    allocate(fcstFB(FBCount), fcstItemNameList(FBCount), fcstItemTypeList(FBCount))
!
! pull out the item names and item types from fcstState
    call ESMF_StateGet(fcstState, itemNameList=fcstItemNameList, &
                       itemTypeList=fcstItemTypeList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
! loop over all items in the fcstState and collect all FieldBundles
    do i=1, FBcount
      if (fcstItemTypeList(i) == ESMF_STATEITEM_FIELDBUNDLE) then
        ! access the FieldBundle
        call ESMF_StateGet(fcstState, itemName=fcstItemNameList(i), &
                           fieldbundle=fcstFB(i), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!    if(mype==0.or.mype==144) print *,'af fcstFB,i=',i,'name=',trim(fcstItemNameList(i))
      else

      !***### anything but a FieldBundle in the state is unexpected here
          call ESMF_LogSetError(ESMF_RC_ARG_BAD,                                 &
                                msg="Only FieldBundles supported in fcstState.", &
                                line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
      endif
    enddo
!
!-----------------------------------------------------------------------
!***  create and initialize Write component(s).
!-----------------------------------------------------------------------
!
    if( quilting ) then

!
      allocate(wrtComp(write_groups), wrtState(write_groups) )
      allocate(wrtFB(5,write_groups), routehandle(5,write_groups))
      allocate(lead_wrttask(write_groups), last_wrttask(write_groups))
!
      allocate(petList(wrttasks_per_group))
      if(mype == 0) print *,'af allco wrtComp,write_groups=',write_groups

! set up ESMF time interval at center of iau window 
      call ESMF_TimeIntervalSet(IAU_offsetTI, h=iau_offset, rc=rc)
!
      allocate(originPetList(num_pes_fcst+wrttasks_per_group))
      allocate(targetPetList(num_pes_fcst+wrttasks_per_group))
!
      k = num_pes_fcst
      timerhs = mpi_wtime()
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
                                     clock=clock_fv3, phase=1, userRc=urc, rc=rc)
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
          write(msg,"(A,I2.2,',',I2.2,A)") "calling into wrtFB(",j,i, ") FieldBundleRegridStore()...."
          call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)

          if (i==1) then
! this is a Store() for the first wrtComp -> must do the Store()
            timewri = mpi_wtime()

            call ESMF_FieldBundleRegridStore(fcstFB(j), wrtFB(j,i),                                    &
                                             regridMethod=regridmethod, routehandle=routehandle(j,i),  &
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
          write(msg,"(A,I2.2,',',I2.2,A)") "... returned from wrtFB(",j,i, ") FieldBundleRegridStore()."
          call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)
        enddo

! end write_groups
      enddo
      if(mype==0) print *,'in fv3cap init, time wrtcrt/regrdst',mpi_wtime()-timerhs
      deallocate(petList)
      deallocate(originPetList)
      deallocate(targetPetList)
!
!---------------------------------------------------------------------------------
!---  SET UP ALARM
!
!--- for every time step output, overwrite nfhout

      if(nsout > 0) then
        nfhout = int(nsout*dt_atmos/3600.)
        nfmout = int((nsout*dt_atmos-nfhout*3600.)/60.)
        nfsout = int(nsout*dt_atmos-nfhout*3600.-nfmout*60)
      else
        nfmout = 0
        nfsout = 0
      endif
      call ESMF_TimeIntervalSet(output_interval, h=nfhout, m=nfmout,  s=nfsout, rc=rc)
      if(mype==0) print *,'af set up output_interval,rc=',rc,'nfhout=',nfhout,nfmout,nfsout

      if (nfhmax_hf > 0 .and. nsout <= 0) then

        nfmout_hf = 0; nfsout_hf = 0
        call ESMF_TimeIntervalSet(output_interval_hf, h=nfhout_hf, m=nfmout_hf, &
                                  s=nfsout_hf, rc=rc)
        call ESMF_TimeIntervalSet(output_hfmax, h=nfhmax_hf, m=0, s=0, rc=rc)
        alarm_output_hf_stop = starttime + output_hfmax + output_interval_hf 
        if (currtime <= starttime+output_hfmax) then
          nhf = (currtime-starttime)/output_interval_hf
          alarm_output_hf_ring = startTime + (nhf+1_ESMF_KIND_I4)*output_interval_hf
          if(iau_offset > 0) then
            alarm_output_hf_ring = startTime + IAU_offsetTI
            if( currtime > alarm_output_hf_ring ) then
              alarm_output_hf_ring = startTime + (nhf+1_ESMF_KIND_I4)*output_interval_hf
            endif
          endif
          alarm_output_hf = ESMF_AlarmCreate(clock_fv3,name='ALARM_OUTPUT_HF',  &
                                             ringTime =alarm_output_hf_ring,    &
                                             ringInterval =output_interval_hf,  &  !<-- Time interval between
                                             stoptime =alarm_output_hf_stop,    &  !<-- Time interval between
                                             ringTimeStepCount=1,               &  !<-- The Alarm rings for this many timesteps
                                             sticky           =.false.,         &  !<-- Alarm does not ring until turned off
                                             rc               =RC)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          alarm_output_ring = startTime + output_hfmax + output_interval
        else
          nrg = (currtime-starttime-output_hfmax)/output_interval
          alarm_output_ring = startTime + output_hfmax + (nrg+1_ESMF_KIND_I4) * output_interval
        endif
      else
        nrg = (currtime-starttime)/output_interval
        alarm_output_ring = startTime + (nrg+1_ESMF_KIND_I4) * output_interval
        if(iau_offset > 0) then
          alarm_output_ring = startTime + IAU_offsetTI
          if( currtime > alarm_output_ring ) then
            alarm_output_ring = startTime + (nrg+1_ESMF_KIND_I4) * output_interval
          endif
        endif
      endif

      call ESMF_TimeIntervalSet(output_interval, h=nfhout, m=nfmout, &
                                s=nfsout, rc=rc)
      alarm_output = ESMF_AlarmCreate(clock_fv3, name  ='ALARM_OUTPUT',    &
                                      ringTime         =alarm_output_ring, & !<-- Forecast/Restart start time (ESMF)
                                      ringInterval     =output_interval,   & !<-- Time interval between
                                      ringTimeStepCount=1,                 & !<-- The Alarm rings for this many timesteps
                                      sticky           =.false.,           & !<-- Alarm does not ring until turned off
                                      rc               =RC)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
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
    ! --- advertise Fields in importState and exportState -------------------

    if( cpl ) then

      isPetLocal = ESMF_GridCompIsPetLocal(fcstComp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if (isPetLocal) then
    
        ! importable fields:
        do i = 1, size(ImportFieldsList)
          if (importFieldShare(i)) then
            call NUOPC_Advertise(importState,         &
                                 StandardName=trim(ImportFieldsList(i)), &
                                 SharePolicyField="share", rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          else
            call NUOPC_Advertise(importState,         &
                                 StandardName=trim(ImportFieldsList(i)), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          end if
        end do
      
        ! exportable fields:
        do i = 1, size(exportFieldsList)
          if (exportFieldShare(i)) then
            call NUOPC_Advertise(exportState,                            &
                                 StandardName=trim(exportFieldsList(i)), &
                                 SharePolicyField="share", rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          else
            call NUOPC_Advertise(exportState,         &
                                 StandardName=trim(exportFieldsList(i)), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          end if
        end do
      
      endif
      if(mype==0) print *,'in fv3_cap, aft import, export fields in atmos'
    endif

    if(mype==0) print *,'in fv3_cap, init time=',mpi_wtime()-timeis
!-----------------------------------------------------------------------
!
  end subroutine InitializeAdvertise

  !-----------------------------------------------------------------------------

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    logical :: isPetLocal
    integer :: n

    rc = ESMF_SUCCESS

    ! --- conditionally realize or remove Fields in importState and exportState -------------------

    if ( cpl ) then

      isPetLocal = ESMF_GridCompIsPetLocal(fcstComp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,  file=__FILE__)) return

      if (isPetLocal) then

        ! -- realize connected fields in exportState
        call realizeConnectedCplFields(exportState, fcstGrid,                                                &
                                       numLevels, numSoilLayers, numTracers, num_diag_sfc_emis_flux,         &
                                       num_diag_down_flux, num_diag_type_down_flux, num_diag_burn_emis_flux, &
                                       num_diag_cmass, exportFieldsList, exportFieldTypes, 'FV3 Export',     &
                                       exportFields, rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,  file=__FILE__)) return

        ! -- realize connected fields in importState
        call realizeConnectedCplFields(importState, fcstGrid,                                                &
                                       numLevels, numSoilLayers, numTracers, num_diag_sfc_emis_flux,         &
                                       num_diag_down_flux, num_diag_type_down_flux, num_diag_burn_emis_flux, &
                                       num_diag_cmass, importFieldsList, importFieldTypes, 'FV3 Import',     &
                                       importFields, rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,  file=__FILE__)) return
      end if
!jw
      call fillExportFields(exportData)
    endif

  end subroutine InitializeRealize

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)                    :: gcomp
    integer, intent(out)                   :: rc
    
    ! local variables
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Clock)                       :: clock
    type(ESMF_Time)                        :: currTime
    type(ESMF_TimeInterval)                :: timeStep
    type(ESMF_Time)                        :: startTime, stopTime
    type(ESMF_TimeInterval)                :: time_elapsed
    integer(ESMF_KIND_I8)                  :: n_interval, time_elapsed_sec
!
    integer :: na, i, urc
    logical :: isAlarmEnabled, isAlarmRinging, lalarm, reconcileFlag
    character(len=*),parameter  :: subname='(fv3_cap:ModelAdvance)'
    character(240)              :: msgString
    character(240)              :: import_timestr, export_timestr
!jw debug
    character(ESMF_MAXSTR)      :: name
    integer :: mype,date(6), fieldcount, fcst_nfld
    real(kind=ESMF_KIND_R4), pointer  :: dataPtr(:,:,:), dataPtr2d(:,:)
    character(64)  :: fcstbdl_name
    real(kind=8)   :: MPI_Wtime
    real(kind=8)   :: timeri, timewri, timewr, timerhi, timerh


!-----------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (profile_memory) call ESMF_VMLogMemInfo("Entering FV3 Model_ADVANCE: ")

    timeri = mpi_wtime()
!    
    call ESMF_GridCompGet(gcomp, name=name, localpet=mype, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in 
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.
    
    call ESMF_ClockPrint(clock_fv3, options="currTime", &
                         preString="------>Advancing FV3 from: ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!-----------------------------------------------------------------------
!***  Use the internal Clock set by NUOPC layer for FV3 but update stopTime
!-----------------------------------------------------------------------

    ! Component internal Clock gets updated per NUOPC rules
    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! The stopTime will be updated to be the next coupling time
    call ESMF_ClockGet(clock, currTime=currTime, stopTime=stopTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Set the coupling time to be stopTime in Clock that FV3 core uses
    call ESMF_ClockSet(clock_fv3, currTime=currTime, stopTime=stopTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_ClockPrint(clock_fv3, options="currTime", &
                         preString="entering FV3_ADVANCE with clock_fv3 current: ", &
                         unit=nuopcMsg)
    call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
    call ESMF_ClockPrint(clock_fv3, options="startTime", &
                         preString="entering FV3_ADVANCE with clock_fv3 start:   ", &
                         unit=nuopcMsg)
    call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
    call ESMF_ClockPrint(clock_fv3, options="stopTime", &
                         preString="entering FV3_ADVANCE with clock_fv3 stop:    ", &
                         unit=nuopcMsg)
    call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)

    call ESMF_ClockGet(clock_fv3, startTime=startTime, currTime=currTime, &
                       timeStep=timeStep, stopTime=stopTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!    if(mype==0)  print *,'total steps=', nint((stopTime-startTime)/timeStep)
!    if(mype==lead_wrttask(1))  print *,'on wrt lead,total steps=', nint((stopTime-startTime)/timeStep)
    call ESMF_TimeGet(time=stopTime,yy=date(1),mm=date(2),dd=date(3),h=date(4), &
                      m=date(5),s=date(6),rc=rc)
!     if(mype==0) print *,'af clock,stop date=',date
!     if(mype==lead_wrttask(1)) print *,'on wrt lead,af clock,stop date=',date
    call ESMF_TimeIntervalGet(timeStep,yy=date(1),mm=date(2),d=date(3),h=date(4), &
                      m=date(5),s=date(6),rc=rc)
!     if(mype==0) print *,'af clock,timestep date=',date
!     if(mype==lead_wrttask(1)) print *,'on wrt lead,af clock,timestep date=',date
!
      call ESMF_ClockGet(clock_fv3, currTime=currTime, timeStep=timeStep, rc=rc)
      call ESMF_TimeGet(currTime,          timestring=import_timestr, rc=rc)
      call ESMF_TimeGet(currTime+timestep, timestring=export_timestr, rc=rc)

!
!-----------------------------------------------------------------------------
!*** integration loop

    reconcileFlag = .true.
    call esmf_clockget(clock_fv3, timestep=timestep, starttime=starttime, &
                       currtime=currtime, rc=rc)

    integrate: do while(.NOT.ESMF_ClockIsStopTime(clock_fv3, rc = RC))
!
!*** for forecast tasks

      timewri = mpi_wtime()
      call ESMF_LogWrite('Model Advance: before fcstcomp run ', ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_GridCompRun(fcstComp, exportState=fcstState, clock=clock_fv3, &
                            phase=1, userRc=urc, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      if ( cpl ) then
       ! assign import_data called during phase=1
       if( dbug > 0 .or. cplprint_flag ) then
        call diagnose_cplFields(gcomp, importState, exportstate, clock_fv3,    &
                              cplprint_flag, dbug, 'import', import_timestr)
       endif
      endif

      call ESMF_GridCompRun(fcstComp, exportState=fcstState, clock=clock_fv3, &
                            phase=2, userRc=urc, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call ESMF_LogWrite('Model Advance: after fcstcomp run ', ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_ClockAdvance(clock = clock_fv3, rc = RC)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call esmf_clockget(clock_fv3, currtime=currtime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      time_elapsed  = currtime - starttime
      na = nint(time_elapsed/timeStep)
!
!    if(mype==0) print *,'in fv3_cap,in model run, advance,na=',na

!-------------------------------------------------------------------------------
!*** if alarms ring, call data transfer and write grid comp run
     if( quilting ) then

       lalarm = .false.
       if (nfhmax_hf > 0) then

         if(currtime <= starttime+output_hfmax) then
           isAlarmEnabled = ESMF_AlarmIsEnabled(alarm = ALARM_OUTPUT_HF, rc = RC)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
           if(isAlarmEnabled) then
             isAlarmRinging = ESMF_AlarmIsRinging(alarm = ALARM_OUTPUT_HF,rc = Rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
             if (isAlarmRinging) LALARM = .true.
           endif
         else
           isAlarmEnabled = ESMF_AlarmIsEnabled(alarm = ALARM_OUTPUT, rc = RC)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
           if(isAlarmEnabled) then
             isAlarmRinging = ESMF_AlarmIsRinging(alarm = ALARM_OUTPUT,rc = Rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
             if (isAlarmRinging) LALARM = .true.
           endif
         endif
       endif
!
       isAlarmEnabled = ESMF_AlarmIsEnabled(alarm = ALARM_OUTPUT, rc = RC)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       if(isAlarmEnabled) then
         isAlarmRinging = ESMF_AlarmIsRinging(alarm = ALARM_OUTPUT,rc = Rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
         if (isAlarmRinging) LALARM = .true.
       endif
!      if (mype == 0 .or. mype == lead_wrttask(1)) print *,' aft fcst run lalarm=',lalarm, &
!      'FBcount=',FBcount,'na=',na

       output: IF(lalarm .or. na==first_kdt ) then

         timerhi = mpi_wtime()
!         if (mype == 0 .or. mype == lead_wrttask(1)) print *,' aft fcst run alarm is on, na=',na,'mype=',mype
         
         call ESMF_VMEpochEnter(epoch=ESMF_VMEpoch_Buffer, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

         do i=1, FBCount
!
! get fcst fieldbundle
!
           call ESMF_FieldBundleRegrid(fcstFB(i), wrtFB(i,n_group),         &
                                       routehandle=routehandle(i, n_group), &
                                       termorderflag=(/ESMF_TERMORDER_SRCSEQ/), rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
!end FBcount
         enddo
         call ESMF_VMEpochExit(rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
         timerh = mpi_wtime()
         if (mype == 0 .or. mype == lead_wrttask(n_group)) print *,'aft fieldbundleregrid,na=',na,  &
           ' time=', timerh- timerhi

!      if(mype==0 .or. mype==lead_wrttask(1))  print *,'on wrt bf wrt run, na=',na
          call ESMF_LogWrite('Model Advance: before wrtcomp run ', ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          timerhi = mpi_wtime()
          call ESMF_GridCompRun(wrtComp(n_group), importState=wrtState(n_group), clock=clock_fv3,userRc=urc,rc=rc)
          timerh = mpi_wtime()
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
!       if (mype == 0 .or. mype == lead_wrttask(n_group)) print *,'aft wrtgridcomp run,na=',na,  &
!        ' time=', timerh- timerhi

          call ESMF_LogWrite('Model Advance: after wrtcomp run ', ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!       if (mype == 0 .or. mype == lead_wrttask(n_group)) print *,'fv3_cap,aft model advance,na=', &
!       na,' time=', mpi_wtime()- timewri


          if(n_group == write_groups) then
            n_group = 1
          else
            n_group = n_group + 1
          endif

        endif output

! end quilting
      endif

!      if (mype == 0 .or. mype == 1536 .or. mype==2160) then
!        print *,'fv3_cap,end integrate,na=',na,' time=',mpi_wtime()- timewri
!      endif

!*** end integreate loop
    enddo integrate
!
!jw for coupled, check clock and dump import and export state
    if ( cpl ) then
      if( dbug > 0 .or. cplprint_flag ) then
       call diagnose_cplFields(gcomp, importState, exportstate, clock_fv3,    &
                              cplprint_flag, dbug, 'export', export_timestr) 
     end if
    endif

    if (mype==0) print *,'fv3_cap,end integrate,na=',na,' time=',mpi_wtime()- timeri

    if (profile_memory) call ESMF_VMLogMemInfo("Leaving FV3 Model_ADVANCE: ")

  end subroutine ModelAdvance

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance_phase1(gcomp, rc)
    type(ESMF_GridComp)                    :: gcomp
    integer, intent(out)                   :: rc
    
    ! local variables
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Clock)                       :: clock
    type(ESMF_Time)                        :: currTime
    type(ESMF_TimeInterval)                :: timeStep
    type(ESMF_Time)                        :: startTime, stopTime
    type(ESMF_TimeInterval)                :: time_elapsed
    integer(ESMF_KIND_I8)                  :: n_interval, time_elapsed_sec
!
    integer :: na, i, urc
    logical :: lalarm, reconcileFlag
    character(len=*),parameter  :: subname='(fv3_cap:ModelAdvance_phase1)'
    character(240)              :: msgString
!jw debug
    character(ESMF_MAXSTR)            :: name
    integer                           :: mype,date(6), fieldcount, fcst_nfld
    real(kind=ESMF_KIND_R4), pointer  :: dataPtr(:,:,:), dataPtr2d(:,:)
    character(64) :: fcstbdl_name
    real(kind=8)  :: MPI_Wtime
    real(kind=8)  :: timewri, timewr, timerhi, timerh

!-----------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if(profile_memory) call ESMF_VMLogMemInfo("Entering FV3 Model_ADVANCE phase1: ")
!    
    call ESMF_GridCompGet(gcomp, name=name, localpet=mype, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Expecting to be called by NUOPC run method exactly once for every coupling
    ! step.
    ! Also expecting the coupling step to be identical to the timeStep for 
    ! clock_fv3.
    
    call ESMF_ClockPrint(clock_fv3, options="currTime", &
                         preString="------>Advancing FV3 from: ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    
!-----------------------------------------------------------------------
!***  Use the internal Clock set by NUOPC layer for FV3 but update stopTime
!-----------------------------------------------------------------------

    ! Component internal Clock gets updated per NUOPC rules
    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! The stopTime will be updated to be the next external coupling time
    call ESMF_ClockGet(clock, currTime=currTime, stopTime=stopTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Set the FV3-OCN coupling time to be stopTime in Clock that FV3 core uses
    call ESMF_ClockSet(clock_fv3, currTime=currTime, stopTime=stopTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_ClockPrint(clock_fv3, options="currTime", &
                         preString="entering FV3_ADVANCE phase1 with clock_fv3 current: ", &
                         unit=nuopcMsg)
    call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
    call ESMF_ClockPrint(clock_fv3, options="startTime", &
                         preString="entering FV3_ADVANCE phase1 with clock_fv3 start:   ", &
                         unit=nuopcMsg)
    call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
    call ESMF_ClockPrint(clock_fv3, options="stopTime", &
                         preString="entering FV3_ADVANCE phase1 with clock_fv3 stop:    ", &
                         unit=nuopcMsg)
    call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)

    call ESMF_ClockGet(clock_fv3, startTime=startTime, currTime=currTime, &
                       timeStep=timeStep, stopTime=stopTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!    if(mype==0)  print *,'total steps=', nint((stopTime-startTime)/timeStep)
!    if(mype==lead_wrttask(1))  print *,'on wrt lead,total steps=', nint((stopTime-startTime)/timeStep)
    call ESMF_TimeGet(time=stopTime,yy=date(1),mm=date(2),dd=date(3),h=date(4), &
                      m=date(5),s=date(6),rc=rc)
!     if(mype==0) print *,'af clock,stop date=',date
!     if(mype==lead_wrttask(1)) print *,'on wrt lead,af clock,stop date=',date
    call ESMF_TimeIntervalGet(timeStep,yy=date(1),mm=date(2),d=date(3),h=date(4), &
                              m=date(5),s=date(6),rc=rc)
!     if(mype==0) print *,'af clock,timestep date=',date
!     if(mype==lead_wrttask(1)) print *,'on wrt lead,af clock,timestep date=',date
!

!-----------------------------------------------------------------------------
!*** no integration loop here!

    reconcileFlag = .true.

!*** for forecast tasks

    call ESMF_LogWrite('Model Advance phase1: before fcstcomp run ', ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompRun(fcstComp, exportState=fcstState, clock=clock_fv3, &
                          phase=1, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    call ESMF_LogWrite('Model Advance phase1: after fcstcomp run ', ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (profile_memory) call ESMF_VMLogMemInfo("Leaving FV3 Model_ADVANCE phase1: ")

  end subroutine ModelAdvance_phase1

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance_phase2(gcomp, rc)
    type(ESMF_GridComp)                    :: gcomp
    integer, intent(out)                   :: rc

    ! local variables
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Clock)                       :: clock
    type(ESMF_Time)                        :: currTime
    type(ESMF_TimeInterval)                :: timeStep
    type(ESMF_Time)                        :: startTime, stopTime
    type(ESMF_TimeInterval)                :: time_elapsed
    integer(ESMF_KIND_I8)                  :: n_interval, time_elapsed_sec
!
    integer :: na, i, urc
    logical :: isAlarmEnabled, isAlarmRinging, lalarm, reconcileFlag
    character(len=*),parameter  :: subname='(fv3_cap:ModelAdvance_phase2)'
    character(240)              :: msgString
!jw debug
    character(ESMF_MAXSTR)            :: name
    integer                           :: mype,date(6), fieldcount, fcst_nfld
    real(kind=ESMF_KIND_R4), pointer  :: dataPtr(:,:,:), dataPtr2d(:,:)
    character(64)  :: fcstbdl_name
    real(kind=8)   :: MPI_Wtime
    real(kind=8)   :: timewri, timewr, timerhi, timerh

!-----------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if(profile_memory) &
      call ESMF_VMLogMemInfo("Entering FV3 Model_ADVANCE phase2: ")
!
    call ESMF_GridCompGet(gcomp, name=name, localpet=mype, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!-----------------------------------------------------------------------------
!*** no integration loop

    reconcileFlag = .true.

!
!*** for forecast tasks

      timewri = mpi_wtime()
      call ESMF_LogWrite('Model Advance phase2: before fcstcomp run ', ESMF_LOGMSG_INFO, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_GridCompRun(fcstComp, exportState=fcstState, clock=clock_fv3, &
                            phase=2, userRc=urc, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call ESMF_LogWrite('Model Advance phase2: after fcstcomp run ', ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_ClockAdvance(clock = clock_fv3, rc = RC)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_ClockGet(clock_fv3, startTime=startTime, currTime=currTime, &
                         timeStep=timeStep, stopTime=stopTime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      time_elapsed  = currtime - starttime
      na = nint(time_elapsed/timeStep)
!
     if(mype==0) print *,'n fv3_cap,in model run, advance,na=',na

!-------------------------------------------------------------------------------
!*** if alarms ring, call data transfer and write grid comp run
     if( quilting ) then

       lalarm = .false.
       if (nfhmax_hf > 0) then

         if(currtime <= starttime+output_hfmax) then
           isAlarmEnabled = ESMF_AlarmIsEnabled(alarm = ALARM_OUTPUT_HF, rc = RC)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
           if(isAlarmEnabled) then
             isAlarmRinging = ESMF_AlarmIsRinging(alarm = ALARM_OUTPUT_HF,rc = Rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
             if (isAlarmRinging) LALARM = .true.
           endif
         else
           isAlarmEnabled = ESMF_AlarmIsEnabled(alarm = ALARM_OUTPUT, rc = RC)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
           if(isAlarmEnabled) then
             isAlarmRinging = ESMF_AlarmIsRinging(alarm = ALARM_OUTPUT,rc = Rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
             if (isAlarmRinging) LALARM = .true.
           endif
         endif

       endif
!
       isAlarmEnabled = ESMF_AlarmIsEnabled(alarm = ALARM_OUTPUT, rc = RC)
       if(isAlarmEnabled) then
         isAlarmRinging = ESMF_AlarmIsRinging(alarm = ALARM_OUTPUT,rc = Rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
         if (isAlarmRinging) LALARM = .true.
       endif
       if (mype == 0 .or. mype == lead_wrttask(1)) print *,' aft fcst run lalarm=',lalarm, &
                                                           'FBcount=',FBcount,'na=',na

       output: IF(lalarm .or. na==first_kdt ) then

         timerhi = mpi_wtime()
         call ESMF_VMEpochEnter(epoch=ESMF_VMEpoch_Buffer, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
         do i=1, FBCount
!
! get fcst fieldbundle
!
           call ESMF_FieldBundleRegrid(fcstFB(i), wrtFB(i,n_group),    &
                                       routehandle=routehandle(i, n_group), rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
!end FBcount
         enddo
         call ESMF_VMEpochExit(rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
         if (mype == 0 .or. mype == lead_wrttask(n_group)) print *,'aft fieldbundleregrid,na=',na,  &
                                                                   ' time=', timerh- timerhi

!        if(mype==0 .or. mype==lead_wrttask(1))  print *,'on wrt bf wrt run, na=',na
         call ESMF_LogWrite('Model Advance: before wrtcomp run ', ESMF_LOGMSG_INFO, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

         timerhi = mpi_wtime()
         call ESMF_GridCompRun(wrtComp(n_group), importState=wrtState(n_group), clock=clock_fv3,userRc=urc,rc=rc)

         timerh = mpi_wtime()

         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
         if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

         if (mype == 0 .or. mype == lead_wrttask(n_group)) print *,'aft wrtgridcomp run,na=',na,  &
                                                                   ' time=', timerh- timerhi

         call ESMF_LogWrite('Model Advance: after wrtcomp run ', ESMF_LOGMSG_INFO, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

         if (mype == 0 .or. mype == lead_wrttask(n_group)) print *,'fv3_cap,aft model advance phase2,na=', &
                                                                   na,' time=', mpi_wtime()- timewri

         if(n_group == write_groups) then
           n_group = 1
         else
           n_group = n_group + 1
         endif

       endif output

! end quilting
     endif

!
!jw check clock
     call ESMF_ClockPrint(clock_fv3, options="currTime", &
                          preString="leaving FV3_ADVANCE phase2 with clock_fv3 current: ", &
                          unit=nuopcMsg)
     call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
     call ESMF_ClockPrint(clock_fv3, options="startTime", &
                          preString="leaving FV3_ADVANCE phase2 with clock_fv3 start:   ", &
                          unit=nuopcMsg)
     call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
     call ESMF_ClockPrint(clock_fv3, options="stopTime", &
                          preString="leaving FV3_ADVANCE phase2 with clock_fv3 stop:    ", &
                          unit=nuopcMsg)
     call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)

     if(profile_memory) &
       call ESMF_VMLogMemInfo("Leaving FV3 Model_ADVANCE phase2: ")

  end subroutine ModelAdvance_phase2

  !-----------------------------------------------------------------------------

  subroutine fv3_checkimport(gcomp, rc)
!
!***  Check the import state fields
!
    type(ESMF_GridComp)   :: gcomp
    integer, intent(out)  :: rc
!
    integer :: n, nf
    type(ESMF_Clock)   :: clock
    type(ESMF_Time)    :: currTime, invalidTime
    type(ESMF_State)   :: importState
    logical            :: timeCheck1,timeCheck2
    type(ESMF_Field),pointer  :: fieldList(:)
    character(len=128) :: fldname
    character(esmf_maxstr) :: MESSAGE_CHECK
!jwtest
    integer date(6)

    ! query the Component for its clock
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! get the current time out of the clock
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!jwtest:
    date(1:6) = 0
    call ESMF_TimeGet(time=CurrTime,yy=date(1),mm=date(2),dd=date(3),h=date(4), &
                      m=date(5),s=date(6),rc=rc)
!   if(fcstmype==0) print *,'in fv3_checkimport, currtime=',date(1:6)

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
!     if(fcstmype==0) print *,'in fv3_checkimport, inside associated(fieldList)'
      do n = 1,size(fieldList)
        call ESMF_FieldGet(fieldList(n), name=fldname, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        nf = queryFieldList(ImportFieldsList,fldname)
        timeCheck1 = NUOPC_IsAtTime(fieldList(n), invalidTime, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if (timeCheck1) then
          importFieldsValid(nf) = .false.
!         if(fcstmype==0) print *,'in fv3_checkimport,',trim(fldname),' is set unvalid, nf=',nf,' at time',date(1:6)
        else
          timeCheck2 = NUOPC_IsAtTime(fieldList(n), currTime, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if (.not.timeCheck2) then
            !TODO: introduce and use INCOMPATIBILITY return codes!!!!
            call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
                                  msg="NUOPC INCOMPATIBILITY DETECTED: Import Field not at current time", &
                                  line=__LINE__, file=__FILE__, rcToReturn=rc)
              return
          endif
        endif
        write(MESSAGE_CHECK,'(A,2i4,l3)') &
          "FV3_CHECKIMPORT "//trim(fldname),n,nf,importFieldsValid(nf)
        CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
      enddo
    endif

  end subroutine fv3_checkimport

  !-----------------------------------------------------------------------------

  subroutine TimestampExport_phase1(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)     :: driverClock, modelClock
    type(ESMF_State)     :: exportState

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

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine atmos_model_finalize(gcomp, rc)

    ! input arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
!
    integer              :: i, unit, date(6), mype, urc
    type(ESMF_VM)        :: vm
    real(kind=8)         :: MPI_Wtime, timeffs
!
!-----------------------------------------------------------------------------
!*** finialize forecast

    timeffs = mpi_wtime()
    rc = ESMF_SUCCESS
!
    call ESMF_GridCompGet(gcomp,vm=vm,rc=rc)
    call ESMF_VMGet(vm, localpet = mype,rc=rc)
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
    if(mype==0)print *,' wrt grid comp destroy time=',mpi_wtime()-timeffs


  end subroutine atmos_model_finalize
!#######################################################################
!
!
!-----------------------------------------------------------------------------

end module fv3gfs_cap_mod
