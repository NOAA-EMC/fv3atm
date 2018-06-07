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
  public Dump_cplFields
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
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)

      if (medAtmCouplingIntervalSec>0._ESMF_KIND_R8) then
        ! The coupling time step was provided
        call ESMF_TimeIntervalSet(fv3Step, s_r8=medAtmCouplingIntervalSec, rc=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        call ESMF_GridCompGet(gcomp, clock=fv3Clock, rc=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        call ESMF_ClockSet(fv3Clock, timestep=fv3Step, rc=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
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
      integer                         :: i
      real(ESMF_KIND_R8), pointer     :: fptr(:,:)

      if (present(rc)) rc = ESMF_SUCCESS

      fieldName = standardName  ! use standard name as field name

      !! Create fields using wam2dmesh if they are WAM fields
      if (NUOPC_IsConnected(state, fieldName=fieldName)) then

        field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, name=fieldName, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call NUOPC_Realize(state, field=field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return

        call ESMF_FieldGet(field, farrayPtr=fptr, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        fptr=0.d0 ! zero out the entire field
        call NUOPC_SetAttribute(field, name="Updated", value="true", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

      else
        ! remove a not connected Field from State
        call ESMF_StateRemove(state, (/fieldName/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

    end subroutine realizeConnectedInternCplField

  !-----------------------------------------------------------------------------

    subroutine realizeConnectedCplFields(state, grid, &
      numLevels, numSoilLayers, numTracers, fieldNames, fieldTypes, fieldList, rc)

      type(ESMF_State),            intent(inout)  :: state
      type(ESMF_Grid),                intent(in)  :: grid
      integer,                        intent(in)  :: numLevels
      integer,                        intent(in)  :: numSoilLayers
      integer,                        intent(in)  :: numTracers
      character(len=*), dimension(:), intent(in)  :: fieldNames
      character(len=*), dimension(:), intent(in)  :: fieldTypes
      type(ESMF_Field), dimension(:), intent(out) :: fieldList
      integer,                        intent(out) :: rc

      ! local variables
      integer          :: item
      type(ESMF_Field) :: field

      ! begin
      rc = ESMF_SUCCESS

      if (size(fieldNames) /= size(fieldTypes)) then
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_SIZE, &
          msg="fieldNames and fieldTypes must have same size.", &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      end if

      do item = 1, size(fieldNames)
        if (NUOPC_IsConnected(state, fieldName=trim(fieldNames(item)))) then
          select case (fieldTypes(item))
            case ('l','layer')
              field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, &
                name=trim(fieldNames(item)), &
                ungriddedLBound=(/1/), ungriddedUBound=(/numLevels/), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out
            case ('i','interface')
              field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, &
                name=trim(fieldNames(item)), &
                ungriddedLBound=(/1/), ungriddedUBound=(/numLevels+1/), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out
            case ('t','tracer')
              field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, &
                name=trim(fieldNames(item)), &
                ungriddedLBound=(/1,1/), ungriddedUBound=(/numLevels, numTracers/), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out
            case ('s','surface')
              field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, &
                name=trim(fieldNames(item)), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out
            case ('g','soil')
              field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, &
                name=trim(fieldNames(item)), &
                ungriddedLBound=(/1/), ungriddedUBound=(/numSoilLayers/), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out
            case default
              call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
                msg="exportFieldType = '"//trim(fieldTypes(item))//"' not recognized", &
                line=__LINE__, file=__FILE__, rcToReturn=rc)
              return
          end select
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return
          ! -- save field
          fieldList(item) = field
        else
          ! remove a not connected Field from State
          call ESMF_StateRemove(state, (/trim(fieldNames(item))/), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        end if
      end do

    end subroutine realizeConnectedCplFields

  !-----------------------------------------------------------------------------

    subroutine Dump_cplFields(gcomp, importState, exportstate, clock_fv3,    &
         statewrite_flag, timeslice)

      type(ESMF_GridComp), intent(in)       :: gcomp
      type(ESMF_State)                      :: importState, exportstate
      type(ESMF_Clock),intent(in)           :: clock_fv3
      logical, intent(in)                   :: statewrite_flag
      integer                               :: timeslice
!
      character(len=160) :: nuopcMsg
      integer :: rc
!
      call ESMF_ClockPrint(clock_fv3, options="currTime", &
        preString="leaving FV3_ADVANCE with clock_fv3 current: ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      call ESMF_ClockPrint(clock_fv3, options="startTime", &
        preString="leaving FV3_ADVANCE with clock_fv3 start:   ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      call ESMF_ClockPrint(clock_fv3, options="stopTime", &
        preString="leaving FV3_ADVANCE with clock_fv3 stop:    ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)

      ! Dumping Fields out
      if (statewrite_flag) then
        timeslice = timeslice + 1
        call ESMF_GridCompGet(gcomp, importState=importState, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMFPP_RegridWriteState(importState, "fv3_cap_import_", timeslice, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_GridCompGet(gcomp, exportState=exportState, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        call ESMFPP_RegridWriteState(exportState, "fv3_cap_export_", timeslice, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
!
    end subroutine Dump_cplFields

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

      outGrid = ESMF_GridCreate1PeriDimUfrm( maxIndex=(/180,360/), &
        minCornerCoord=(/0.0_ESMF_KIND_R8,-90.0_ESMF_KIND_R8/), &
        maxCornerCoord=(/360.0_ESMF_KIND_R8,90.0_ESMF_KIND_R8/), &
        staggerLocList=(/ESMF_STAGGERLOC_CORNER, ESMF_STAGGERLOC_CENTER/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call ESMF_StateGet(state, itemCount=icount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      allocate(typeList(icount), itemNameList(icount))
      call ESMF_StateGet(state, itemTypeList=typeList, itemNameList=itemNameList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      do i = 1, icount
        if(typeList(i) == ESMF_STATEITEM_FIELD) then
          call ESMF_LogWrite("RegridWrite Field Name Initiated: "//trim(itemNameList(i)), ESMF_LOGMSG_INFO)
          call ESMF_StateGet(state, itemName=itemNameList(i), field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          call ESMFPP_RegridWrite(field, outGrid, ESMF_REGRIDMETHOD_BILINEAR, &
            fileName//trim(itemNameList(i))//'.nc', trim(itemNameList(i)), timeslice, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          call ESMF_LogWrite("RegridWrite Field Name done: "//trim(itemNameList(i)), ESMF_LOGMSG_INFO)
        endif
      enddo

      deallocate(typeList, itemNameList)

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

      ! local arguments
      type(ESMF_Routehandle)                   :: rh
      type(ESMF_Field)                         :: outField

      outField = ESMF_FieldCreate(outGrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! For other options for the regrid operation, please refer to:
      ! http://www.earthsystemmodeling.org/esmf_releases/last_built/ESMF_refdoc/node5.html#SECTION050366000000000000000
      call ESMF_FieldRegridStore(inField, outField, regridMethod=regridMethod, &
        unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
        Routehandle=rh, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call ESMF_FieldRegrid(inField, outField, Routehandle=rh, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call ESMF_FieldWrite(outField, fileName, variableName=fieldName, timeslice=timeslice, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      rc = ESMF_SUCCESS

    end subroutine ESMFPP_RegridWrite


  !-----------------------------------------------------------------------------

end module module_cap_cpl
