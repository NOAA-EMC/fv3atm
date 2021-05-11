module module_cplfields

  !-----------------------------------------------------------------------------
  ! This module contains the fv3 Coupling Fields: export and import
  !
  !-----------------------------------------------------------------------------

  use ESMF

  implicit none

  private

  type, public :: FieldInfo
    character(len=41) :: name
    character(len=1) :: type
    character(len=9) :: shared
  end type

! Export Fields ----------------------------------------
  integer,          public, parameter :: NexportFields = 72
  type(ESMF_Field), target, public    :: exportFields(NexportFields)

  type(FieldInfo), dimension(NexportFields), public, parameter :: exportFieldsInfo = [ &
    FieldInfo("inst_pres_interface                      ", "i", "share    " ), &
    FieldInfo("inst_pres_levels                         ", "l", "share    " ), &
    FieldInfo("inst_geop_interface                      ", "i", "share    " ), &
    FieldInfo("inst_geop_levels                         ", "l", "share    " ), &
    FieldInfo("inst_temp_levels                         ", "l", "share    " ), &
    FieldInfo("inst_zonal_wind_levels                   ", "l", "share    " ), &
    FieldInfo("inst_merid_wind_levels                   ", "l", "share    " ), &
    FieldInfo("inst_omega_levels                        ", "l", "share    " ), &
    FieldInfo("inst_tracer_mass_frac                    ", "t", "share    " ), &
    FieldInfo("soil_type                                ", "s", "share    " ), &
    FieldInfo("inst_pbl_height                          ", "s", "share    " ), &
    FieldInfo("surface_cell_area                        ", "s", "share    " ), &
    FieldInfo("inst_convective_rainfall_amount          ", "s", "share    " ), &
    FieldInfo("inst_exchange_coefficient_heat_levels    ", "l", "share    " ), &
    FieldInfo("inst_spec_humid_conv_tendency_levels     ", "l", "share    " ), &
    FieldInfo("inst_friction_velocity                   ", "s", "share    " ), &
    FieldInfo("inst_rainfall_amount                     ", "s", "share    " ), &
    FieldInfo("inst_soil_moisture_content               ", "g", "share    " ), &
    FieldInfo("inst_up_sensi_heat_flx                   ", "s", "share    " ), &
    FieldInfo("inst_lwe_snow_thickness                  ", "s", "share    " ), &
    FieldInfo("vegetation_type                          ", "s", "share    " ), &
    FieldInfo("inst_vegetation_area_frac                ", "s", "share    " ), &
    FieldInfo("inst_surface_roughness                   ", "s", "share    " ), &
    FieldInfo("mean_zonal_moment_flx_atm                ", "s", "not share" ), &
    FieldInfo("mean_merid_moment_flx_atm                ", "s", "not share" ), &
    FieldInfo("mean_sensi_heat_flx                      ", "s", "not share" ), &
    FieldInfo("mean_laten_heat_flx                      ", "s", "not share" ), &
    FieldInfo("mean_down_lw_flx                         ", "s", "not share" ), &
    FieldInfo("mean_down_sw_flx                         ", "s", "not share" ), &
    FieldInfo("mean_prec_rate                           ", "s", "not share" ), &
    FieldInfo("inst_zonal_moment_flx                    ", "s", "not share" ), &
    FieldInfo("inst_merid_moment_flx                    ", "s", "not share" ), &
    FieldInfo("inst_sensi_heat_flx                      ", "s", "not share" ), &
    FieldInfo("inst_laten_heat_flx                      ", "s", "not share" ), &
    FieldInfo("inst_down_lw_flx                         ", "s", "not share" ), &
    FieldInfo("inst_down_sw_flx                         ", "s", "share    " ), &
    FieldInfo("inst_temp_height2m                       ", "s", "not share" ), &
    FieldInfo("inst_spec_humid_height2m                 ", "s", "not share" ), &
    FieldInfo("inst_zonal_wind_height10m                ", "s", "not share" ), &
    FieldInfo("inst_merid_wind_height10m                ", "s", "not share" ), &
    FieldInfo("inst_temp_height_surface                 ", "s", "share    " ), &
    FieldInfo("inst_pres_height_surface                 ", "s", "not share" ), &
    FieldInfo("inst_surface_height                      ", "s", "not share" ), &
    FieldInfo("mean_net_lw_flx                          ", "s", "not share" ), &
    FieldInfo("mean_net_sw_flx                          ", "s", "not share" ), &
    FieldInfo("inst_net_lw_flx                          ", "s", "not share" ), &
    FieldInfo("inst_net_sw_flx                          ", "s", "not share" ), &
    FieldInfo("mean_down_sw_ir_dir_flx                  ", "s", "not share" ), &
    FieldInfo("mean_down_sw_ir_dif_flx                  ", "s", "not share" ), &
    FieldInfo("mean_down_sw_vis_dir_flx                 ", "s", "not share" ), &
    FieldInfo("mean_down_sw_vis_dif_flx                 ", "s", "not share" ), &
    FieldInfo("inst_down_sw_ir_dir_flx                  ", "s", "not share" ), &
    FieldInfo("inst_down_sw_ir_dif_flx                  ", "s", "not share" ), &
    FieldInfo("inst_down_sw_vis_dir_flx                 ", "s", "not share" ), &
    FieldInfo("inst_down_sw_vis_dif_flx                 ", "s", "not share" ), &
    FieldInfo("mean_net_sw_ir_dir_flx                   ", "s", "not share" ), &
    FieldInfo("mean_net_sw_ir_dif_flx                   ", "s", "not share" ), &
    FieldInfo("mean_net_sw_vis_dir_flx                  ", "s", "not share" ), &
    FieldInfo("mean_net_sw_vis_dif_flx                  ", "s", "not share" ), &
    FieldInfo("inst_net_sw_ir_dir_flx                   ", "s", "not share" ), &
    FieldInfo("inst_net_sw_ir_dif_flx                   ", "s", "not share" ), &
    FieldInfo("inst_net_sw_vis_dir_flx                  ", "s", "not share" ), &
    FieldInfo("inst_net_sw_vis_dif_flx                  ", "s", "not share" ), &
    FieldInfo("inst_land_sea_mask                       ", "s", "share    " ), &
    FieldInfo("inst_temp_height_lowest                  ", "s", "not share" ), &
    FieldInfo("inst_spec_humid_height_lowest            ", "s", "not share" ), &
    FieldInfo("inst_zonal_wind_height_lowest            ", "s", "not share" ), &
    FieldInfo("inst_merid_wind_height_lowest            ", "s", "not share" ), &
    FieldInfo("inst_pres_height_lowest                  ", "s", "not share" ), &
    FieldInfo("inst_height_lowest                       ", "s", "not share" ), &
    FieldInfo("mean_fprec_rate                          ", "s", "not share" ), &
    FieldInfo("openwater_frac_in_atm                    ", "s", "not share" ) ]

  real(kind=8), allocatable, public :: exportData(:,:,:)

! Import Fields ----------------------------------------
  integer,          public, parameter :: NimportFields = 17
  logical,          public            :: importFieldsValid(NimportFields)
  type(ESMF_Field), target, public    :: importFields(NimportFields)

  type(FieldInfo), dimension(NimportFields), public, parameter :: importFieldsInfo = [ &
    FieldInfo("inst_tracer_mass_frac                    ", "t", "share    " ), &
    FieldInfo("land_mask                                ", "s", "not share" ), &
    FieldInfo("sea_ice_surface_temperature              ", "s", "not share" ), &
    FieldInfo("sea_surface_temperature                  ", "s", "not share" ), &
    FieldInfo("ice_fraction                             ", "s", "not share" ), &
    FieldInfo("mean_up_lw_flx_ice                       ", "s", "not share" ), &
    FieldInfo("mean_laten_heat_flx_atm_into_ice         ", "s", "not share" ), &
    FieldInfo("mean_sensi_heat_flx_atm_into_ice         ", "s", "not share" ), &
    FieldInfo("stress_on_air_ice_zonal                  ", "s", "not share" ), &
    FieldInfo("stress_on_air_ice_merid                  ", "s", "not share" ), &
    FieldInfo("mean_ice_volume                          ", "s", "not share" ), &
    FieldInfo("mean_snow_volume                         ", "s", "not share" ), &
    FieldInfo("inst_tracer_up_surface_flx               ", "u", "share    " ), &
    FieldInfo("inst_tracer_down_surface_flx             ", "d", "share    " ), &
    FieldInfo("inst_tracer_clmn_mass_dens               ", "c", "share    " ), &
    FieldInfo("inst_tracer_anth_biom_flx                ", "b", "share    " ), &
    FieldInfo("wave_z0_roughness_length                 ", "s", "not share" ) ]

  ! Methods
  public fillExportFields
  public queryImportFields, queryExportFields
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
    character(len=ESMF_MAXSTR)                  :: fieldName
    real(kind=ESMF_KIND_R4), dimension(:,:), pointer   :: datar42d
    real(kind=ESMF_KIND_R8), dimension(:,:), pointer   :: datar82d

!
    if (present(rc)) rc=ESMF_SUCCESS

    do n=1, size(exportFields)
      isCreated = ESMF_FieldIsCreated(exportFields(n), rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      if (isCreated) then
! set data
        call ESMF_FieldGet(exportFields(n), name=fieldname, dimCount=dimCount, typekind=datatype, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
        !print *,'in fillExportFields, field created n=',n,size(exportFields),'name=', trim(fieldname)
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
  integer function queryExportFields(fieldname, abortflag)

    character(len=*),intent(in) :: fieldname
    logical, optional           :: abortflag

    queryExportFields = queryFieldList(exportFieldsInfo, fieldname, abortflag)

  end function queryExportFields

  integer function queryImportFields(fieldname, abortflag)

    character(len=*),intent(in) :: fieldname
    logical, optional           :: abortflag

    queryImportFields = queryFieldList(importFieldsInfo, fieldname, abortflag)

  end function queryImportFields


  integer function queryFieldList(fieldsInfo, fieldname, abortflag)
    ! returns integer index of first found fieldname in fieldlist
    ! by default, will abort if field not found, set abortflag to false
    ! to turn off the abort.
    ! return value of < 1 means the field was not found

    type(FieldInfo) ,intent(in) :: fieldsInfo(:)
    character(len=*),intent(in) :: fieldname
    logical, optional           :: abortflag

    integer :: n
    logical :: labort
    integer :: rc

    labort = .true.
    if (present(abortflag)) then
      labort = abortflag
    endif

    queryFieldList = 0
    n = 1
    do while (queryFieldList < 1 .and. n <= size(fieldsInfo))
      if (trim(fieldsInfo(n)%name) == trim(fieldname)) then
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
