module module_cplfields

  !-----------------------------------------------------------------------------
  ! This module contains the fv3 Coupling Fields: export and import
  !
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC

  implicit none

  private

! Export Fields ----------------------------------------
  integer,          public, parameter :: NexportFields = 2
  type(ESMF_Field), target, public    :: exportFields(NexportFields)
  character(len=*), public, parameter :: exportFieldsList(NexportFields) = (/ &
       "inst_zonal_wind_height10m                ", &
       "inst_merid_wind_height10m                "  &
  /)
  ! Field types should be provided for proper handling
  ! according to the table below:
  !  g : soil levels (3D)
  !  i : interface (3D)
  !  l : model levels (3D)
  !  s : surface (2D)
  !  t : tracers (4D)
  character(len=*), public, parameter :: exportFieldTypes(NexportFields) = (/ &
       "s","s"                                      &
  /)
  ! Set exportFieldShare to .true. if field is provided as memory reference
  ! to coupled components
  logical, public, parameter :: exportFieldShare(NexportFields) = (/ &
       .false.,.false.                              &
  /)
  real(kind=8), allocatable, public :: exportData(:,:,:)

! Import Fields ----------------------------------------
  integer,          public, parameter :: NimportFields = 1
  logical,          public            :: importFieldsValid(NimportFields)
  type(ESMF_Field), target, public    :: importFields(NimportFields)
  character(len=*), public, parameter :: importFieldsList(NimportFields) = (/ &
       "land_mask                              "    &
  /)
  character(len=*), public, parameter :: importFieldTypes(NimportFields) = (/ &
       "s"                                  &
  /)
  ! Set importFieldShare to .true. if field is provided as memory reference
  ! from coupled components
  logical, public, parameter :: importFieldShare(NimportFields) = (/ &
       .false.                                  &
  /)

  ! Methods
  public fillExportFields
  public queryFieldList
  public cplFieldGet

!-----------------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------------

  subroutine fillExportFields(data_a2oi, rc)
    ! Fill updated data into the export Fields.
    real(kind=8), target, intent(in)            :: data_a2oi(:,:,:)
    integer, intent(out), optional              :: rc

    integer                                     :: localrc
    integer                                     :: n,dimCount
    logical                                     :: isCreated
    type(ESMF_TypeKind_Flag)                    :: datatype
    real(kind=ESMF_KIND_R4), dimension(:,:), pointer   :: datar42d
    real(kind=ESMF_KIND_R8), dimension(:,:), pointer   :: datar82d
    
!
    if (present(rc)) rc=ESMF_SUCCESS

    do n=1, size(exportFields)
      isCreated = ESMF_FieldIsCreated(exportFields(n), rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      if (isCreated) then
! set data 
        call ESMF_FieldGet(exportFields(n), dimCount=dimCount, typekind=datatype, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
        if ( datatype == ESMF_TYPEKIND_R8) then
           if ( dimCount == 2) then
             call ESMF_FieldGet(exportFields(n),farrayPtr=datar82d,localDE=0, rc=localrc)
             if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
             datar82d = data_a2oi(:,:,n)
           endif
        else if ( datatype == ESMF_TYPEKIND_R4) then
           if ( dimCount == 2) then
             call ESMF_FieldGet(exportFields(n),farrayPtr=datar82d,localDE=0, rc=localrc)
             if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
             datar42d = data_a2oi(:,:,n)
           endif
        endif
      endif
    enddo
  end subroutine fillExportFields
!
!------------------------------------------------------------------------------
!
  integer function queryFieldList(fieldlist, fieldname, abortflag, rc)
    ! returns integer index of first found fieldname in fieldlist
    ! by default, will abort if field not found, set abortflag to false
    ! to turn off the abort.
    ! return value of < 1 means the field was not found

    character(len=*),intent(in) :: fieldlist(:)
    character(len=*),intent(in) :: fieldname
    logical, optional           :: abortflag
    integer, optional           :: rc

    integer :: n
    logical :: labort

    labort = .true.
    if (present(abortflag)) then
      labort = abortflag
    endif

    queryFieldList = 0
    n = 1
    do while (queryFieldList < 1 .and. n <= size(fieldlist))
      if (trim(fieldlist(n)) == trim(fieldname)) then
        queryFieldList = n
      else
        n = n + 1
      endif
    enddo

    if (labort .and. queryFieldList < 1) then
      call ESMF_LogWrite('queryFieldList ABORT on fieldname '//trim(fieldname), &
                          ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=rc)
      CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    endif
  end function queryFieldList
!
!------------------------------------------------------------------------------
!
  subroutine cplStateGet(state, fieldList, fieldCount, rc)

    character(len=*), intent(in)            :: state
    type(ESMF_Field), pointer,     optional :: fieldList(:)
    integer,          intent(out), optional :: fieldCount
    integer,          intent(out), optional :: rc

    !--- begin
    if (present(rc)) rc = ESMF_SUCCESS

    select case (trim(state))
      case ('import','i')
        if (present(fieldList )) fieldList  => importFields
        if (present(fieldCount)) fieldCount =  size(importFields)
      case ('export','o')
        if (present(fieldList )) fieldList  => exportFields
        if (present(fieldCount)) fieldCount =  size(exportFields)
      case default
        call ESMF_LogSetError(ESMF_RC_ARG_OUTOFRANGE, &
          msg="state argument can only be import(i)/export(o).", &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
    end select

  end subroutine cplStateGet


  subroutine cplFieldGet(state, name, localDe, &
                         farrayPtr2d, farrayPtr3d, farrayPtr4d, rc)

    character(len=*),   intent(in)            :: state
    character(len=*),   intent(in)            :: name
    integer,            intent(in),  optional :: localDe
    real(ESMF_KIND_R8), pointer,     optional :: farrayPtr2d(:,:)
    real(ESMF_KIND_R8), pointer,     optional :: farrayPtr3d(:,:,:)
    real(ESMF_KIND_R8), pointer,     optional :: farrayPtr4d(:,:,:,:)
    integer,            intent(out), optional :: rc

    !--- local variables
    integer                    :: localrc
    integer                    :: de, item, fieldCount, rank
    logical                    :: isCreated
    type(ESMF_Field), pointer  :: fieldList(:)
    character(len=ESMF_MAXSTR) :: fieldName

    !--- begin
    if (present(rc)) rc = ESMF_SUCCESS

    if (present(farrayPtr2d)) nullify(farrayPtr2d)
    if (present(farrayPtr3d)) nullify(farrayPtr3d)
    if (present(farrayPtr4d)) nullify(farrayPtr4d)

    de = 0
    if (present(localDe)) de = localDe

    call cplStateGet(state, fieldList=fieldList, fieldCount=fieldCount, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    do item = 1, fieldCount
      isCreated = ESMF_FieldIsCreated(fieldList(item), rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      if (isCreated) then
        call ESMF_FieldGet(fieldList(item), name=fieldName, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
        if (trim(fieldName) == trim(name)) then
          call ESMF_FieldGet(fieldList(item), rank=rank, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
          select case (rank)
            case (2)
              if (present(farrayPtr2d)) then
                call ESMF_FieldGet(fieldList(item), localDe=de, farrayPtr=farrayPtr2d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
              end if
            case (3)
              if (present(farrayPtr3d)) then
                call ESMF_FieldGet(fieldList(item), localDe=de, farrayPtr=farrayPtr3d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
              end if
            case (4)
              if (present(farrayPtr4d)) then
                call ESMF_FieldGet(fieldList(item), localDe=de, farrayPtr=farrayPtr4d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
              end if
            case default
              call ESMF_LogSetError(ESMF_RC_NOT_IMPL, msg="field rank should be 2, 3, or 4.", &
                                    line=__LINE__, file=__FILE__, rcToReturn=rc)
              return
          end select
          exit
        end if
      end if
    end do

  end subroutine cplFieldGet
!
!------------------------------------------------------------------------------
!
end module module_cplfields
