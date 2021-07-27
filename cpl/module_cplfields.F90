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
  end type

! Export Fields ----------------------------------------

  ! Please specify fields as: FieldInfo("standard_name", "type")
  ! Field types should be provided according to the table below:
  !  g : soil levels (3D)
  !  i : interface (3D)
  !  l : model levels (3D)
  !  s : surface (2D)
  !  t : tracers (4D)
  integer,          public, parameter :: NexportFields = 105
  type(ESMF_Field), target, public    :: exportFields(NexportFields)

  type(FieldInfo), dimension(NexportFields), public, parameter :: exportFieldsInfo = [ &
    FieldInfo("inst_pres_interface                      ", "i"), &
    FieldInfo("inst_pres_levels                         ", "l"), &
    FieldInfo("inst_geop_interface                      ", "i"), &
    FieldInfo("inst_geop_levels                         ", "l"), &
    FieldInfo("inst_temp_levels                         ", "l"), &
    FieldInfo("inst_zonal_wind_levels                   ", "l"), &
    FieldInfo("inst_merid_wind_levels                   ", "l"), &
    FieldInfo("inst_omega_levels                        ", "l"), &
    FieldInfo("inst_tracer_mass_frac                    ", "t"), &
    FieldInfo("soil_type                                ", "s"), &
    FieldInfo("inst_pbl_height                          ", "s"), &
    FieldInfo("surface_cell_area                        ", "s"), &
    FieldInfo("inst_convective_rainfall_amount          ", "s"), &
    FieldInfo("inst_exchange_coefficient_heat_levels    ", "l"), &
    FieldInfo("inst_spec_humid_conv_tendency_levels     ", "l"), &
    FieldInfo("inst_ice_nonconv_tendency_levels         ", "l"), &
    FieldInfo("inst_liq_nonconv_tendency_levels         ", "l"), &
    FieldInfo("inst_cloud_frac_levels                   ", "l"), &
    FieldInfo("inst_friction_velocity                   ", "s"), &
    FieldInfo("inst_rainfall_amount                     ", "s"), &
    FieldInfo("inst_soil_moisture_content               ", "g"), &
    FieldInfo("inst_surface_soil_wetness                ", "s"), &
    FieldInfo("inst_up_sensi_heat_flx                   ", "s"), &
    FieldInfo("inst_lwe_snow_thickness                  ", "s"), &
    FieldInfo("vegetation_type                          ", "s"), &
    FieldInfo("inst_vegetation_area_frac                ", "s"), &
    FieldInfo("inst_surface_roughness                   ", "s"), &
    FieldInfo("mean_zonal_moment_flx_atm                ", "s"), &
    FieldInfo("mean_merid_moment_flx_atm                ", "s"), &
    FieldInfo("mean_sensi_heat_flx                      ", "s"), &
    FieldInfo("mean_laten_heat_flx                      ", "s"), &
    FieldInfo("mean_down_lw_flx                         ", "s"), &
    FieldInfo("mean_down_sw_flx                         ", "s"), &
    FieldInfo("mean_prec_rate                           ", "s"), &
    FieldInfo("inst_zonal_moment_flx                    ", "s"), &
    FieldInfo("inst_merid_moment_flx                    ", "s"), &
    FieldInfo("inst_sensi_heat_flx                      ", "s"), &
    FieldInfo("inst_laten_heat_flx                      ", "s"), &
    FieldInfo("inst_down_lw_flx                         ", "s"), &
    FieldInfo("inst_down_sw_flx                         ", "s"), &
    FieldInfo("inst_temp_height2m                       ", "s"), &
    FieldInfo("inst_spec_humid_height2m                 ", "s"), &
    FieldInfo("inst_zonal_wind_height10m                ", "s"), &
    FieldInfo("inst_merid_wind_height10m                ", "s"), &
    FieldInfo("inst_temp_height_surface                 ", "s"), &
    FieldInfo("inst_pres_height_surface                 ", "s"), &
    FieldInfo("inst_surface_height                      ", "s"), &
    FieldInfo("mean_net_lw_flx                          ", "s"), &
    FieldInfo("mean_net_sw_flx                          ", "s"), &
    FieldInfo("inst_net_lw_flx                          ", "s"), &
    FieldInfo("inst_net_sw_flx                          ", "s"), &
    FieldInfo("mean_down_sw_ir_dir_flx                  ", "s"), &
    FieldInfo("mean_down_sw_ir_dif_flx                  ", "s"), &
    FieldInfo("mean_down_sw_vis_dir_flx                 ", "s"), &
    FieldInfo("mean_down_sw_vis_dif_flx                 ", "s"), &
    FieldInfo("inst_down_sw_ir_dir_flx                  ", "s"), &
    FieldInfo("inst_down_sw_ir_dif_flx                  ", "s"), &
    FieldInfo("inst_down_sw_vis_dir_flx                 ", "s"), &
    FieldInfo("inst_down_sw_vis_dif_flx                 ", "s"), &
    FieldInfo("mean_net_sw_ir_dir_flx                   ", "s"), &
    FieldInfo("mean_net_sw_ir_dif_flx                   ", "s"), &
    FieldInfo("mean_net_sw_vis_dir_flx                  ", "s"), &
    FieldInfo("mean_net_sw_vis_dif_flx                  ", "s"), &
    FieldInfo("inst_net_sw_ir_dir_flx                   ", "s"), &
    FieldInfo("inst_net_sw_ir_dif_flx                   ", "s"), &
    FieldInfo("inst_net_sw_vis_dir_flx                  ", "s"), &
    FieldInfo("inst_net_sw_vis_dif_flx                  ", "s"), &
    FieldInfo("inst_land_sea_mask                       ", "s"), &
    FieldInfo("inst_temp_height_lowest                  ", "s"), &
    FieldInfo("inst_spec_humid_height_lowest            ", "s"), &
    FieldInfo("inst_zonal_wind_height_lowest            ", "s"), &
    FieldInfo("inst_merid_wind_height_lowest            ", "s"), &
    FieldInfo("inst_pres_height_lowest                  ", "s"), &
    FieldInfo("inst_height_lowest                       ", "s"), &
    FieldInfo("mean_fprec_rate                          ", "s"), &
    FieldInfo("openwater_frac_in_atm                    ", "s"), &
    FieldInfo("ice_fraction_in_atm                      ", "s"), &
    FieldInfo("lake_fraction                            ", "s"), &
    FieldInfo("ocean_fraction                           ", "s"), &
    FieldInfo("surface_snow_area_fraction               ", "s"), &


    !  For JEDI
    ! dynamics
    FieldInfo("u                                        ", "l"), &
    FieldInfo("v                                        ", "l"), &
    FieldInfo("ua                                       ", "l"), &
    FieldInfo("va                                       ", "l"), &
    FieldInfo("t                                        ", "l"), &
    FieldInfo("delp                                     ", "l"), &
    FieldInfo("sphum                                    ", "l"), &
    FieldInfo("ice_wat                                  ", "l"), &
    FieldInfo("liq_wat                                  ", "l"), &
    FieldInfo("o3mr                                     ", "l"), &
    FieldInfo("phis                                     ", "s"), &
    FieldInfo("u_srf                                    ", "s"), &
    FieldInfo("v_srf                                    ", "s"), &
    ! physics
    FieldInfo("slmsk                                    ", "s"), &
    FieldInfo("weasd                                    ", "s"), &
    FieldInfo("tsea                                     ", "s"), &
    FieldInfo("vtype                                    ", "s"), &
    FieldInfo("stype                                    ", "s"), &
    FieldInfo("vfrac                                    ", "s"), &
    FieldInfo("stc                                      ", "g"), &
    FieldInfo("smc                                      ", "g"), &
    FieldInfo("snwdph                                   ", "s"), &
    FieldInfo("f10m                                     ", "s"), &
    FieldInfo("zorl                                     ", "s"), &
    FieldInfo("t2m                                      ", "s") ]

! Import Fields ----------------------------------------
  integer,          public, parameter :: NimportFields = 42
  logical,          public            :: importFieldsValid(NimportFields)
  type(ESMF_Field), target, public    :: importFields(NimportFields)

  type(FieldInfo), dimension(NimportFields), public, parameter :: importFieldsInfo = [ &
    FieldInfo("inst_tracer_mass_frac                    ", "t"), &
    FieldInfo("land_mask                                ", "s"), &
    FieldInfo("sea_ice_surface_temperature              ", "s"), &
    FieldInfo("sea_surface_temperature                  ", "s"), &
    FieldInfo("ice_fraction                             ", "s"), &
    FieldInfo("mean_up_lw_flx_ice                       ", "s"), &
    FieldInfo("mean_laten_heat_flx_atm_into_ice         ", "s"), &
    FieldInfo("mean_sensi_heat_flx_atm_into_ice         ", "s"), &
    FieldInfo("stress_on_air_ice_zonal                  ", "s"), &
    FieldInfo("stress_on_air_ice_merid                  ", "s"), &
    FieldInfo("mean_ice_volume                          ", "s"), &
    FieldInfo("mean_snow_volume                         ", "s"), &
    FieldInfo("inst_ice_ir_dif_albedo                   ", "s"), &
    FieldInfo("inst_ice_ir_dir_albedo                   ", "s"), &
    FieldInfo("inst_ice_vis_dif_albedo                  ", "s"), &
    FieldInfo("inst_ice_vis_dir_albedo                  ", "s"), &
    FieldInfo("wave_z0_roughness_length                 ", "s"), &

    !  For JEDI
    ! dynamics
    FieldInfo("u                                        ", "l"), &
    FieldInfo("v                                        ", "l"), &
    FieldInfo("ua                                       ", "l"), &
    FieldInfo("va                                       ", "l"), &
    FieldInfo("t                                        ", "l"), &
    FieldInfo("delp                                     ", "l"), &
    FieldInfo("sphum                                    ", "l"), &
    FieldInfo("ice_wat                                  ", "l"), &
    FieldInfo("liq_wat                                  ", "l"), &
    FieldInfo("o3mr                                     ", "l"), &
    FieldInfo("phis                                     ", "s"), &
    FieldInfo("u_srf                                    ", "s"), &
    FieldInfo("v_srf                                    ", "s"), &
    ! physics
    FieldInfo("slmsk                                    ", "s"), &
    FieldInfo("weasd                                    ", "s"), &
    FieldInfo("tsea                                     ", "s"), &
    FieldInfo("vtype                                    ", "s"), &
    FieldInfo("stype                                    ", "s"), &
    FieldInfo("vfrac                                    ", "s"), &
    FieldInfo("stc                                      ", "g"), &
    FieldInfo("smc                                      ", "g"), &
    FieldInfo("snwdph                                   ", "s"), &
    FieldInfo("f10m                                     ", "s"), &
    FieldInfo("zorl                                     ", "s"), &
    FieldInfo("t2m                                      ", "s") ]

  ! Methods
  public queryImportFields, queryExportFields
  public cplFieldGet

!-----------------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------------
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
