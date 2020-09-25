module module_cap_cpl
!
!*** this module contains the debug subroutines for fv3 coupled run
!
! revision history
!  12 Mar 2018: J. Wang       Pull coupled subroutines from fv3_cap.F90 to this module
!
  use esmf
  use NUOPC
!
  implicit none
  private
  public clock_cplIntval
  public realizeConnectedInternCplField
  public realizeConnectedCplFields
  public diagnose_cplFields
!
  contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

    subroutine clock_cplIntval(gcomp, CF)

      type(ESMF_GridComp)      :: gcomp
      type(ESMF_Config)        :: CF
!
      real(ESMF_KIND_R8)       :: medAtmCouplingIntervalSec
      type(ESMF_Clock)         :: fv3Clock
      type(ESMF_TimeInterval)  :: fv3Step
      integer                  :: rc
!
      call ESMF_ConfigGetAttribute(config=CF, value=medAtmCouplingIntervalSec, &
                                   label="atm_coupling_interval_sec:", default=-1.0_ESMF_KIND_R8, rc=RC)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)

      if (medAtmCouplingIntervalSec > 0._ESMF_KIND_R8) then ! The coupling time step is provided
        call ESMF_TimeIntervalSet(fv3Step, s_r8=medAtmCouplingIntervalSec, rc=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        call ESMF_GridCompGet(gcomp, clock=fv3Clock, rc=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        call ESMF_ClockSet(fv3Clock, timestep=fv3Step, rc=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif

    end subroutine clock_cplIntval

  !-----------------------------------------------------------------------------

    subroutine realizeConnectedInternCplField(state, field, standardName, grid, rc)

      type(ESMF_State)                :: state
      type(ESMF_Field), optional      :: field
      character(len=*), optional      :: standardName
      type(ESMF_Grid), optional       :: grid
      integer, intent(out), optional  :: rc

      ! local variables
      character(len=80)               :: fieldName
      type(ESMF_ArraySpec)            :: arrayspec
      integer                         :: i, localrc
      logical                         :: isConnected
      real(ESMF_KIND_R8), pointer     :: fptr(:,:)

      if (present(rc)) rc = ESMF_SUCCESS

      fieldName = standardName  ! use standard name as field name

      !! Create fields using wam2dmesh if they are WAM fields
      isConnected = NUOPC_IsConnected(state, fieldName=fieldName, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      if (isConnected) then

        field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, name=fieldName, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
        call NUOPC_Realize(state, field=field, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call ESMF_FieldGet(field, farrayPtr=fptr, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        fptr=0._ESMF_KIND_R8 ! zero out the entire field
        call NUOPC_SetAttribute(field, name="Updated", value="true", rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      else
        ! remove a not connected Field from State
        call ESMF_StateRemove(state, (/fieldName/), rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      endif

    end subroutine realizeConnectedInternCplField

  !-----------------------------------------------------------------------------

    subroutine realizeConnectedCplFields(state, grid,                                      &
                                         numLevels, numSoilLayers, numTracers,             &
                                         num_diag_sfc_emis_flux, num_diag_down_flux,       &
                                         num_diag_type_down_flux, num_diag_burn_emis_flux, &
                                         num_diag_cmass, fieldNames, fieldTypes, state_tag,&
                                         fieldList, rc)

      type(ESMF_State),            intent(inout)  :: state
      type(ESMF_Grid),                intent(in)  :: grid
      integer,                        intent(in)  :: numLevels
      integer,                        intent(in)  :: numSoilLayers
      integer,                        intent(in)  :: numTracers
      integer,                        intent(in)  :: num_diag_sfc_emis_flux
      integer,                        intent(in)  :: num_diag_down_flux
      integer,                        intent(in)  :: num_diag_type_down_flux
      integer,                        intent(in)  :: num_diag_burn_emis_flux
      integer,                        intent(in)  :: num_diag_cmass
      character(len=*), dimension(:), intent(in)  :: fieldNames
      character(len=*), dimension(:), intent(in)  :: fieldTypes
      character(len=*),               intent(in)  :: state_tag                              !< Import or export.
      type(ESMF_Field), dimension(:), intent(out) :: fieldList
      integer,                        intent(out) :: rc

      ! local variables
      integer          :: item
      logical          :: isConnected
      type(ESMF_Field) :: field

      ! begin
      rc = ESMF_SUCCESS

      if (size(fieldNames) /= size(fieldTypes)) then
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SIZE, &
          msg="fieldNames and fieldTypes must have same size.", line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      end if

      do item = 1, size(fieldNames)
        isConnected = NUOPC_IsConnected(state, fieldName=trim(fieldNames(item)), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        if (isConnected) then
          call ESMF_StateGet(state, field=field, itemName=trim(fieldNames(item)), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          call ESMF_FieldEmptySet(field, grid=grid, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          select case (fieldTypes(item))
            case ('l','layer')
              call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, &
                                       ungriddedLBound=(/1/), ungriddedUBound=(/numLevels/), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            case ('i','interface')
              call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, &
                                       ungriddedLBound=(/1/), ungriddedUBound=(/numLevels+1/), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            case ('t','tracer')
              call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, &
                                       ungriddedLBound=(/1, 1/), ungriddedUBound=(/numLevels, numTracers/), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            case ('u','tracer_up_flux')
              call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, &
                                       ungriddedLBound=(/1/), ungriddedUBound=(/num_diag_sfc_emis_flux/), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            case ('d','tracer_down_flx')
              call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, &
                                       ungriddedLBound=(/1, 1/),        &
                                       ungriddedUBound=(/num_diag_down_flux, num_diag_type_down_flux/), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            case ('b','tracer_anth_biom_emission')
              call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, &
                                       ungriddedLBound=(/1/), ungriddedUBound=(/num_diag_burn_emis_flux/), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            case ('c','tracer_column_mass_density')
              call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, &
                                       ungriddedLBound=(/1/), ungriddedUBound=(/num_diag_cmass/), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            case ('s','surface')
              call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            case ('g','soil')
              call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, &
                                       ungriddedLBound=(/1/), ungriddedUBound=(/numSoilLayers/), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            case default
              call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
                msg="exportFieldType = '"//trim(fieldTypes(item))//"' not recognized", &
                line=__LINE__, file=__FILE__, rcToReturn=rc)
              return
          end select
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          ! -- zero out field 
          call ESMF_FieldFill(field, dataFillScheme="const", const1=0._ESMF_KIND_R8, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          ! -- save field
          fieldList(item) = field
          call ESMF_LogWrite('realizeConnectedCplFields '//trim(state_tag)//' Field '//trim(fieldNames(item))  &
                             // ' is connected ', ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=rc)
        else
          ! remove a not connected Field from State
          call ESMF_StateRemove(state, (/trim(fieldNames(item))/), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          call ESMF_LogWrite('realizeConnectedCplFields '//trim(state_tag)//' Field '//trim(fieldNames(item))  &
                             // ' is not connected ', ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=rc)
        end if
      end do

    end subroutine realizeConnectedCplFields

  !-----------------------------------------------------------------------------

    subroutine diagnose_cplFields(gcomp, importState, exportstate, clock_fv3,    &
         statewrite_flag, stdiagnose_flag, state_tag, timestr)

      type(ESMF_GridComp), intent(in)       :: gcomp
      type(ESMF_State)                      :: importState, exportstate
      type(ESMF_Clock),intent(in)           :: clock_fv3
      logical, intent(in)                   :: statewrite_flag
      integer, intent(in)                   :: stdiagnose_flag
      character(len=*),         intent(in)  :: state_tag                        !< Import or export.
      character(len=*),         intent(in)  :: timestr                          !< Import or export.
      integer                               :: timeslice = 1 
!
      character(len=160) :: nuopcMsg
      character(len=160) :: filename
      integer :: rc
!
      call ESMF_ClockPrint(clock_fv3, options="currTime",                            &
                           preString="leaving FV3_ADVANCE with clock_fv3 current: ", &
                           unit=nuopcMsg)
!      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      call ESMF_ClockPrint(clock_fv3, options="startTime",                           &
                           preString="leaving FV3_ADVANCE with clock_fv3 start:   ", &
                           unit=nuopcMsg)
!      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      call ESMF_ClockPrint(clock_fv3, options="stopTime",                            &
                           preString="leaving FV3_ADVANCE with clock_fv3 stop:    ", &
                           unit=nuopcMsg)
!      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)

      if(trim(state_tag) .eq. 'import')then 
        call ESMF_GridCompGet(gcomp, importState=importState, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if(stdiagnose_flag > 0)then
         call state_diagnose(importState, ':IS', rc=rc)
        end if

        ! Dump Fields out
        if (statewrite_flag) then
          write(filename,'(A)') 'fv3_cap_import_'//trim(timestr)//'_'
          call State_RWFields_tiles(importState,trim(filename), timeslice, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        end if
      end if

      if(trim(state_tag) .eq. 'export')then
        call ESMF_GridCompGet(gcomp, exportState=exportState, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if(stdiagnose_flag > 0)then
         call state_diagnose(exportState, ':ES', rc=rc)
        end if

        ! Dump Fields out
        if (statewrite_flag) then
          write(filename,'(A)') 'fv3_cap_export_'//trim(timestr)//'_'
          call State_RWFields_tiles(exportState,trim(filename), timeslice, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        end if
      end if

    end subroutine diagnose_cplFields

  !-----------------------------------------------------------------------------

    subroutine ESMFPP_RegridWriteState(state, fileName, timeslice, rc)

      type(ESMF_State), intent(in)          :: state
      character(len=*), intent(in)          :: fileName
      integer, intent(in)                   :: timeslice
      integer, intent(out)                  :: rc

      ! local
      type(ESMF_Field)                       :: field
      type(ESMF_Grid)                        :: outGrid
      integer                                :: i, icount
      character(64), allocatable             :: itemNameList(:)
      type(ESMF_StateItem_Flag), allocatable :: typeList(:)

      rc = ESMF_SUCCESS

      ! 1degx1deg
      outGrid = ESMF_GridCreate1PeriDimUfrm(maxIndex=(/360,180/), &
                                            minCornerCoord=(/0.0_ESMF_KIND_R8,-90.0_ESMF_KIND_R8/), &
                                            maxCornerCoord=(/360.0_ESMF_KIND_R8,90.0_ESMF_KIND_R8/), &
                                            staggerLocList=(/ESMF_STAGGERLOC_CORNER, ESMF_STAGGERLOC_CENTER/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_StateGet(state, itemCount=icount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      allocate(typeList(icount), itemNameList(icount))
      call ESMF_StateGet(state, itemTypeList=typeList, itemNameList=itemNameList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      do i = 1, icount
        if(typeList(i) == ESMF_STATEITEM_FIELD) then
          call ESMF_LogWrite("RegridWrite Field Name Initiated: "//trim(itemNameList(i)), ESMF_LOGMSG_INFO)
          call ESMF_StateGet(state, itemName=itemNameList(i), field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          call ESMFPP_RegridWrite(field, outGrid, ESMF_REGRIDMETHOD_BILINEAR, &
                                  fileName//trim(itemNameList(i))//'.nc', trim(itemNameList(i)), timeslice, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          call ESMF_LogWrite("RegridWrite Field Name done: "//trim(itemNameList(i)), ESMF_LOGMSG_INFO)
        endif
      enddo

      deallocate(typeList, itemNameList)

      call ESMF_GridDestroy(outGrid,noGarbage=.true., rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    end subroutine ESMFPP_RegridWriteState

    subroutine ESMFPP_RegridWrite(inField, outGrid, regridMethod, fileName, fieldName, timeslice, rc)

      ! input arguments
      type(ESMF_Field), intent(in)             :: inField
      type(ESMF_Grid), intent(in)              :: outGrid
      type(ESMF_RegridMethod_Flag), intent(in) :: regridMethod
      character(len=*), intent(in)             :: filename
      character(len=*), intent(in)             :: fieldName
      integer,          intent(in)             :: timeslice
      integer,          intent(inout)          :: rc

      ! local variables
      integer                                  :: srcTermProcessing
      type(ESMF_Routehandle)                   :: rh
      type(ESMF_Field)                         :: outField

      outField = ESMF_FieldCreate(outGrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      ! Perform entire regridding arithmetic on the destination PET
      srcTermProcessing = 0
      ! For other options for the regrid operation, please refer to:
      ! http://www.earthsystemmodeling.org/esmf_releases/last_built/ESMF_refdoc/node5.html#SECTION050366000000000000000
      call ESMF_FieldRegridStore(inField, outField, regridMethod=regridMethod, &
                                 unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,    &
                                 srcTermProcessing=srcTermProcessing, Routehandle=rh, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      ! Use fixed ascending order for the sum terms based on their source
      ! sequence index to ensure bit-for-bit reproducibility
      call ESMF_FieldRegrid(inField, outField, Routehandle=rh, &
                            termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_FieldWrite(outField, fileName, variableName=fieldName, timeslice=timeslice, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_FieldRegridRelease(routehandle=rh, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_FieldDestroy(outField,noGarbage=.true., rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      rc = ESMF_SUCCESS

    end subroutine ESMFPP_RegridWrite

  !-----------------------------------------------------------------------------

  ! This subroutine requires ESMFv8 - for coupled FV3
    subroutine State_RWFields_tiles(state,filename,timeslice,rc)

      type(ESMF_State), intent(in)          :: state
      character(len=*), intent(in)          :: fileName
      integer, intent(in)                   :: timeslice
      integer, intent(out)                  :: rc

      ! local
      type(ESMF_Field)                       :: firstESMFFLD
      type(ESMF_Field),allocatable           :: flds(:)
      type(ESMF_GridComp) :: IOComp
      type(ESMF_Grid) :: gridFv3

      character(len=256) :: msgString
      integer                                :: i, icount, ifld
      integer                                :: fieldcount, firstfld
      character(64), allocatable             :: itemNameList(:), fldNameList(:)
      type(ESMF_StateItem_Flag), allocatable :: typeList(:)

      character(len=*),parameter :: subname='(module_cap_cpl:State_RWFields_tiles)'

      ! local variables

      rc = ESMF_SUCCESS
      !call ESMF_LogWrite(trim(subname)//trim(filename)//": called",
      !ESMF_LOGMSG_INFO, rc=rc)

      call ESMF_StateGet(state, itemCount=icount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      allocate(typeList(icount), itemNameList(icount))
      call ESMF_StateGet(state, itemTypeList=typeList, itemNameList=itemNameList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      ! find first stateitem that is a field and the count of fields
      firstfld = 0; fieldcount = 0
      do i = icount,1,-1
        if(typeList(i) == ESMF_STATEITEM_FIELD) firstfld = i
        if(typeList(i) == ESMF_STATEITEM_FIELD) fieldcount = fieldcount + 1
      enddo
      !write(msgString,*) trim(subname)//' icount = ',icount," fieldcount =
      !",fieldcount," firstfld = ",firstfld
      !call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)

      allocate(flds(fieldCount),fldNameList(fieldCount))
      ifld = 1
      do i = 1, icount
        if(typeList(i) == ESMF_STATEITEM_FIELD) then
          fldNameList(ifld) = itemNameList(i)
          ifld = ifld + 1
        endif
      enddo

      call ESMF_LogWrite(trim(subname)//": write "//trim(filename)//"tile1-tile6", ESMF_LOGMSG_INFO, rc=rc)
      ! get first field
      call ESMF_StateGet(state, itemName=itemNameList(firstfld), field=firstESMFFLD, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out

      call ESMF_FieldGet(firstESMFFLD, grid=gridFv3, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out

      IOComp = ESMFIO_Create(gridFv3, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
      call ESMF_LogWrite(trim(subname)//": write "//trim(filename), ESMF_LOGMSG_INFO, rc=rc)

      do ifld=1, fieldCount
        call ESMF_StateGet(state, itemName=fldNameList(ifld), field=flds(ifld), rc=rc)
      enddo

      call ESMFIO_Write(IOComp, filename, flds, filePath='./', rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out

! -- Finalize ESMFIO
      deallocate(flds)
      deallocate(fldNameList)
      call ESMFIO_Destroy(IOComp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) call ESMF_Finalize()

      !call ESMF_LogWrite(trim(subname)//trim(filename)//": finished",
      !ESMF_LOGMSG_INFO, rc=rc)

    end subroutine State_RWFields_tiles

  !-----------------------------------------------------------------------------

  subroutine state_diagnose(State,string, rc)
    ! ----------------------------------------------
    ! Diagnose status of state
    ! ----------------------------------------------
    type(ESMF_State), intent(inout) :: State
    character(len=*), intent(in), optional :: string
    integer, intent(out), optional  :: rc

    ! local variables
    integer                       :: i,j,n
    integer                       :: itemCount
    character(len=64) ,pointer    :: itemNameList(:)
    character(len=64)             :: lstring
    character(len=256)            :: tmpstr

    type(ESMF_Field)            :: lfield
    type(ESMF_StateItem_Flag)   :: itemType
    real(ESMF_KIND_R8), pointer :: dataPtr2d(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr3d(:,:,:)
    integer                     :: lrc, dimCount
    character(len=*),parameter  :: subname='(FV3: state_diagnose)'

    lstring = ''
    if (present(string)) then
       lstring = trim(string)
    endif

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call ESMF_StateGet(State, itemCount=itemCount, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    allocate(itemNameList(itemCount))
   
    call ESMF_StateGet(State, itemNameList=itemNameList, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
   
    do n = 1, itemCount
       call ESMF_StateGet(State, itemName=trim(itemNameList(n)), itemType=itemType, rc=lrc)
       if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
   
     if(itemType == ESMF_STATEITEM_FIELD)then
       call ESMF_StateGet(State, itemName=trim(itemNameList(n)), field=lfield, rc=lrc)
       if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
  
       call ESMF_FieldGet(lfield, dimCount=dimcount, rc=lrc)
       if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       
       if(dimcount == 2)then
         call ESMF_FieldGet(lfield, farrayPtr=dataPtr2d, rc=lrc)
         if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
   
         write(tmpstr,'(A,3g14.7)') trim(subname)//' '//trim(lstring)//':'//trim(itemNameList(n))//'  ', &
           minval(dataPtr2d),maxval(dataPtr2d),sum(dataPtr2d)
         call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=lrc)
       else
         call ESMF_FieldGet(lfield, farrayPtr=dataPtr3d, rc=lrc)
         if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
  
         write(tmpstr,'(A,3g14.7)') trim(subname)//' '//trim(lstring)//':'//trim(itemNameList(n))//'  ', &
           minval(dataPtr3d),maxval(dataPtr3d),sum(dataPtr3d)
         call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=lrc)
       end if
     end if
    enddo
     deallocate(itemNameList)

    if (present(rc)) rc = lrc
    call ESMF_LogWrite(subname//' exit', ESMF_LOGMSG_INFO)

  end subroutine state_diagnose

  !-----------------------------------------------------------------------------

end module module_cap_cpl
