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
  integer,          public, parameter :: NexportFields = 72
  type(ESMF_Field), target, public    :: exportFields(NexportFields)
  character(len=*), public, parameter :: exportFieldsList(NexportFields) = (/ &
       "inst_pres_interface                      ", &
       "inst_pres_levels                         ", &
       "inst_geop_interface                      ", &
       "inst_geop_levels                         ", &
       "inst_temp_levels                         ", &
       "inst_zonal_wind_levels                   ", &
       "inst_merid_wind_levels                   ", &
       "inst_omega_levels                        ", &
       "inst_tracer_mass_frac                    ", &
       "soil_type                                ", &
       "inst_pbl_height                          ", &
       "surface_cell_area                        ", &
       "inst_convective_rainfall_amount          ", &
       "inst_exchange_coefficient_heat_levels    ", &
       "inst_spec_humid_conv_tendency_levels     ", &
       "inst_friction_velocity                   ", &
       "inst_rainfall_amount                     ", &
       "inst_soil_moisture_content               ", &
       "inst_up_sensi_heat_flx                   ", &
       "inst_lwe_snow_thickness                  ", &
       "vegetation_type                          ", &
       "inst_vegetation_area_frac                ", &
       "inst_surface_roughness                   ", &
       "mean_zonal_moment_flx_atm                ", &
       "mean_merid_moment_flx_atm                ", &
       "mean_sensi_heat_flx                      ", &
       "mean_laten_heat_flx                      ", &
       "mean_down_lw_flx                         ", &
       "mean_down_sw_flx                         ", &
       "mean_prec_rate                           ", &
       "inst_zonal_moment_flx                    ", &
       "inst_merid_moment_flx                    ", &
       "inst_sensi_heat_flx                      ", &
       "inst_laten_heat_flx                      ", &
       "inst_down_lw_flx                         ", &
       "inst_down_sw_flx                         ", &
       "inst_temp_height2m                       ", &
       "inst_spec_humid_height2m                 ", &
       "inst_zonal_wind_height10m                ", &
       "inst_merid_wind_height10m                ", &
       "inst_temp_height_surface                 ", &
       "inst_pres_height_surface                 ", &
       "inst_surface_height                      ", &
       "mean_net_lw_flx                          ", &
       "mean_net_sw_flx                          ", &
       "inst_net_lw_flx                          ", &
       "inst_net_sw_flx                          ", &
       "mean_down_sw_ir_dir_flx                  ", &
       "mean_down_sw_ir_dif_flx                  ", &
       "mean_down_sw_vis_dir_flx                 ", &
       "mean_down_sw_vis_dif_flx                 ", &
       "inst_down_sw_ir_dir_flx                  ", &
       "inst_down_sw_ir_dif_flx                  ", &
       "inst_down_sw_vis_dir_flx                 ", &
       "inst_down_sw_vis_dif_flx                 ", &
       "mean_net_sw_ir_dir_flx                   ", &
       "mean_net_sw_ir_dif_flx                   ", &
       "mean_net_sw_vis_dir_flx                  ", &
       "mean_net_sw_vis_dif_flx                  ", &
       "inst_net_sw_ir_dir_flx                   ", &
       "inst_net_sw_ir_dif_flx                   ", &
       "inst_net_sw_vis_dir_flx                  ", &
       "inst_net_sw_vis_dif_flx                  ", &
       "inst_land_sea_mask                       ", &
       "inst_temp_height_lowest                  ", &
       "inst_spec_humid_height_lowest            ", &
       "inst_zonal_wind_height_lowest            ", &
       "inst_merid_wind_height_lowest            ", &
       "inst_pres_height_lowest                  ", &
       "inst_height_lowest                       ", &
       "mean_fprec_rate                          ", &
       "openwater_frac_in_atm                    "  &
!      "northward_wind_neutral                   ", &
!      "eastward_wind_neutral                    ", &
!      "upward_wind_neutral                      ", &
!      "temp_neutral                             ", &
!      "O_Density                                ", &
!      "O2_Density                               ", &
!      "N2_Density                               ", &
!      "height                                   "  &
  /)
  ! Field types should be provided for proper handling
  ! according to the table below:
  !  g : soil levels (3D)
  !  i : interface (3D)
  !  l : model levels (3D)
  !  s : surface (2D)
  !  t : tracers (4D)
  character(len=*), public, parameter :: exportFieldTypes(NexportFields) = (/ &
       "i","l","i","l","l","l","l","l","t", &
       "s","s","s","s","l","l","s","s","g", &
       "s","s","s","s","s","s","s","s",     &
       "s","s","s","s","s","s","s","s",     &
       "s","s","s","s","s","s","s","s",     &
       "s","s","s","s","s","s","s","s",     &
       "s","s","s","s","s","s","s","s",     &
       "s","s","s","s","s","s","s","s",     &
       "s","s","s","s","s","s"              &
!      "l","l","l","l","l","l","l","s",     &
  /)
  ! Set exportFieldShare to .true. if field is provided as memory reference
  ! to coupled components
  logical, public, parameter :: exportFieldShare(NexportFields) = (/ &
       .true. ,.true. ,.true. ,.true. ,.true. , &
       .true. ,.true. ,.true. ,.true. ,.true. , &
       .true. ,.true. ,.true. ,.true. ,.true. , &
       .true. ,.true. ,.true. ,.true. ,.true. , &
       .true. ,.true. ,.true. ,.false.,.false., &
       .false.,.false.,.false.,.false.,.false., &
       .false.,.false.,.false.,.false.,.false. , &
       .true. ,.false.,.false.,.false.,.false. , &
       .true. ,.false.,.false.,.false.,.false., &
       .false.,.false.,.false.,.false.,.false., &
       .false.,.false.,.false.,.false.,.false., &
       .false.,.false.,.false.,.false.,.false., &
       .false.,.false.,.false.,.true. ,.false., &
       .false.,.false.,.false.,.false.,.false., &
       .false.,.false.                          &
!      .false.,.false.,.false.,.false.,.false., &
!      .false.,.false.,.false.                  &
  /)
  real(kind=8), allocatable, public :: exportData(:,:,:)

! Import Fields ----------------------------------------
  integer,          public, parameter :: NimportFields = 17
  logical,          public            :: importFieldsValid(NimportFields)
  type(ESMF_Field), target, public    :: importFields(NimportFields)
  character(len=*), public, parameter :: importFieldsList(NimportFields) = (/ &
       "inst_tracer_mass_frac                  ", &
       "land_mask                              ", &
       "sea_ice_surface_temperature            ", &
       "sea_surface_temperature                ", &
       "ice_fraction                           ", &
!      "inst_ice_ir_dif_albedo                 ", &
!      "inst_ice_ir_dir_albedo                 ", &
!      "inst_ice_vis_dif_albedo                ", &
!      "inst_ice_vis_dir_albedo                ", &
       "mean_up_lw_flx_ice                     ", &
       "mean_laten_heat_flx_atm_into_ice       ", &
       "mean_sensi_heat_flx_atm_into_ice       ", &
!      "mean_evap_rate                         ", &
       "stress_on_air_ice_zonal                ", &
       "stress_on_air_ice_merid                ", &
       "mean_ice_volume                        ", &
       "mean_snow_volume                       ", &
       "inst_tracer_up_surface_flx             ", &
       "inst_tracer_down_surface_flx           ", &
       "inst_tracer_clmn_mass_dens             ", &
       "inst_tracer_anth_biom_flx              ", &
       "wave_z0_roughness_length               "  &
  /)
  character(len=*), public, parameter :: importFieldTypes(NimportFields) = (/ &
       "t",                                 &
       "s","s","s","s","s",                 &
       "s","s","s","s","s",                 &
       "s","u","d","c","b",                 &
       "s"                                  &
  /)
  ! Set importFieldShare to .true. if field is provided as memory reference
  ! from coupled components
  logical, public, parameter :: importFieldShare(NimportFields) = (/ &
       .true. ,                                 &
       .false.,.false.,.false.,.false.,.false., &
       .false.,.false.,.false.,.false.,.false., &
       .false.,.true. ,.true. ,.true. ,.true. , &
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
