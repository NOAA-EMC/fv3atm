module module_cplfields

  !-----------------------------------------------------------------------------
  ! This module contains the fv3 Coupling Fields: export and import
  !
  !-----------------------------------------------------------------------------

  use ESMF
  implicit none

  private

  integer, parameter :: MAXNAMELEN = 128
  
! Export Fields ----------------------------------------
!  integer, public, parameter :: NexportFields = 56
  integer, public, parameter  :: NexportFields = 48
  type(ESMF_Field), public    :: exportFields(NexportFields)
  real(kind=8), allocatable,public     :: exportData(:,:,:)
  character(len=40), public, parameter :: exportFieldsList(NexportFields) = (/ &
      "mean_zonal_moment_flx                  ", &
      "mean_merid_moment_flx                  ", &
      "mean_sensi_heat_flx                    ", &
      "mean_laten_heat_flx                    ", &
      "mean_down_lw_flx                       ", &
      "mean_down_sw_flx                       ", &
      "mean_prec_rate                         ", &
      "inst_zonal_moment_flx                  ", &
      "inst_merid_moment_flx                  ", &
      "inst_sensi_heat_flx                    ", &
      "inst_laten_heat_flx                    ", &
      "inst_down_lw_flx                       ", &
      "inst_down_sw_flx                       ", &
      "inst_temp_height2m                     ", &
      "inst_spec_humid_height2m               ", &
      "inst_zonal_wind_height10m              ", &
      "inst_merid_wind_height10m              ", &
      "inst_temp_height_surface               ", &
      "inst_pres_height_surface               ", &
      "inst_surface_height                    ", &
      "mean_net_lw_flx                        ", &
      "mean_net_sw_flx                        ", &
      "inst_net_lw_flx                        ", &
      "inst_net_sw_flx                        ", &
      "mean_down_sw_ir_dir_flx                ", &
      "mean_down_sw_ir_dif_flx                ", &
      "mean_down_sw_vis_dir_flx               ", &
      "mean_down_sw_vis_dif_flx               ", &
      "inst_down_sw_ir_dir_flx                ", &
      "inst_down_sw_ir_dif_flx                ", &
      "inst_down_sw_vis_dir_flx               ", &
      "inst_down_sw_vis_dif_flx               ", &
      "mean_net_sw_ir_dir_flx                 ", &
      "mean_net_sw_ir_dif_flx                 ", &
      "mean_net_sw_vis_dir_flx                ", &
      "mean_net_sw_vis_dif_flx                ", &
      "inst_net_sw_ir_dir_flx                 ", &
      "inst_net_sw_ir_dif_flx                 ", &
      "inst_net_sw_vis_dir_flx                ", &
      "inst_net_sw_vis_dif_flx                ", &
      "inst_land_sea_mask                     ", &
      "inst_temp_height_lowest                ", &
      "inst_spec_humid_height_lowest          ", &
      "inst_zonal_wind_height_lowest          ", &
      "inst_merid_wind_height_lowest          ", &
      "inst_pres_height_lowest                ", &
      "inst_height_lowest                     ", &
      "mean_fprec_rate                        " &
!      "northward_wind_neutral                 ", &
!      "eastward_wind_neutral                  ", &
!      "upward_wind_neutral                    ", &
!      "temp_neutral                           ", &
!      "O_Density                              ", &
!      "O2_Density                             ", &
!      "N2_Density                             ", &
!      "height                                 "  &
  /)

! Import Fields ----------------------------------------
!  integer, public, parameter :: NimportFields = 16
  integer, public, parameter :: NimportFields = 11
  type(ESMF_Field), public   :: importFields(NimportFields)
  character(len=40), public, parameter :: importFieldsList(NimportFields) = (/ &
       "land_mask                              ", &
       "surface_temperature                    ", &
       "sea_surface_temperature                ", &
       "ice_fraction                           ", &
!      "inst_ice_ir_dif_albedo                 ", &
!      "inst_ice_ir_dir_albedo                 ", &
!      "inst_ice_vis_dif_albedo                ", &
!      "inst_ice_vis_dir_albedo                ", &
       "mean_up_lw_flx                         ", &
       "mean_laten_heat_flx                    ", &
       "mean_sensi_heat_flx                    ", &
!      "mean_evap_rate                         ", &
       "mean_zonal_moment_flx                  ", &
       "mean_merid_moment_flx                  ", &
       "mean_ice_volume                        ", &
       "mean_snow_volume                       "  &
  /)

  ! Methods
  public fillExportFields
  public queryFieldList

!-----------------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------------

  subroutine fillExportFields(data_a2oi, rc)
    ! Fill updated data into the export Fields.
    real(kind=8), target, intent(in)            :: data_a2oi(:,:,:)
    integer, intent(out), optional              :: rc

    integer                               :: n,dimCount
    type(ESMF_TypeKind_Flag)              :: datatype
    real(kind=ESMF_KIND_R4), dimension(:,:), pointer   :: datar42d
    real(kind=ESMF_KIND_R8), dimension(:,:), pointer   :: datar82d
    
!
    if (present(rc)) rc=ESMF_SUCCESS

    do n=1, size(exportFields)
      if (ESMF_FieldIsCreated(exportFields(n))) then
! set data 
        call ESMF_FieldGet(exportFields(n), dimCount=dimCount ,typekind=datatype, rc=rc)
        if ( datatype == ESMF_TYPEKIND_R8) then
           if ( dimCount == 2) then
             call ESMF_FieldGet(exportFields(n),farrayPtr=datar82d,localDE=0, rc=rc)
             datar82d=data_a2oi(:,:,n)
           endif
        else if ( datatype == ESMF_TYPEKIND_R4) then
           if ( dimCount == 2) then
             call ESMF_FieldGet(exportFields(n),farrayPtr=datar82d,localDE=0, rc=rc)
             datar42d=data_a2oi(:,:,n)
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
    !   to turn off the abort.
    ! return value of < 1 means the field was not found

    character(len=*),intent(in) :: fieldlist(:)
    character(len=*),intent(in) :: fieldname
    logical, optional :: abortflag
    integer, optional :: rc

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
!------------------------------------------------------------------------------
  end module module_cplfields
