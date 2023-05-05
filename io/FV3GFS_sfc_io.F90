module FV3GFS_sfc_io

  use block_control_mod,  only: block_control_type
  use fms2_io_mod,        only: FmsNetcdfDomainFile_t, unlimited, &
       register_axis, register_restart_field,       &
       register_variable_attribute, register_field, &
       read_restart, write_restart, write_data,     &
       get_global_io_domain_indices, variable_exists
  use FV3GFS_common_io,   only: copy_from_GFS_Data, copy_to_GFS_Data
  use GFS_typedefs,       only: GFS_sfcprop_type, GFS_control_type, kind_phys
  use GFS_restart,        only: GFS_restart_type
  use mpp_mod,            only: mpp_error,  mpp_pe, mpp_root_pe, &
       mpp_chksum, NOTE,   FATAL
  use physcons,           only: con_tice          !saltwater freezing temp (K)

  implicit none
  private

  public :: Sfc_io_data_type
  public :: Sfc_io_fill_2d_names, Sfc_io_fill_3d_names, Sfc_io_calculate_indices, Sfc_io_final

  type Sfc_io_data_type
    integer, public :: nvar2o = 0
    integer, public :: nvar3 = 0
    integer, public :: nvar2r = 0
    integer, public :: nvar2mp = 0
    integer, public :: nvar3mp = 0
    integer, public :: nvar2l = 0
    integer, public :: nvar2m = 0
    integer, public :: nvar_before_lake = 0

    ! The lsoil flag is only meaningful when reading:;
    logical, public :: is_lsoil = .false.

    ! SYNONYMS: Some nvar variables had two names in FV3GFS_io.F90. They have
    ! only one name here. The "_s" is redundant because this file only has
    ! surface restart variables.
    !
    !  - nvar2m = nvar_s2m
    !  - nvar2o = nvar_s2o
    !  - nvar2r = nvar_s2r
    !  - nvar2mp = nvar_s2mp
    !  - nvar3mp = nvar_s3mp

    real(kind=kind_phys), pointer, dimension(:,:,:), public :: var2 => null()
    real(kind=kind_phys), pointer, dimension(:,:,:), public :: var3ice => null()
    real(kind=kind_phys), pointer, dimension(:,:,:,:), public :: var3 => null()
    real(kind=kind_phys), pointer, dimension(:,:,:,:), public :: var3sn => null()
    real(kind=kind_phys), pointer, dimension(:,:,:,:), public :: var3eq => null()
    real(kind=kind_phys), pointer, dimension(:,:,:,:), public :: var3zn => null()

    character(len=32), pointer, dimension(:), public :: name2 => null()

    character(len=32), pointer, dimension(:), public :: name3 => null()
  contains
    procedure, public :: fill_2d_names => Sfc_io_fill_2d_names
    procedure, public :: fill_3d_names => Sfc_io_fill_3d_names
    procedure, public :: calculate_indices => Sfc_io_calculate_indices
    final :: Sfc_io_final
  end type Sfc_io_data_type

contains

  function Sfc_io_calculate_indices(sfc, Model, for_write, warm_start)
    implicit none
    class(Sfc_io_data_type)             :: sfc
    type(GFS_control_type),      intent(in) :: Model
    logical :: Sfc_io_calculate_indices
    logical, intent(in) :: for_write, warm_start

    integer :: nvar2m, nvar2o, nvar3, nvar2r, nvar2mp, nvar3mp, nvar2l
    integer :: nvar_before_lake

    nvar2m = 48
    if (Model%use_cice_alb .or. Model%lsm == Model%lsm_ruc) then
      nvar2m = nvar2m + 4
      !nvar2m = nvar2m + 5
    endif
    if (Model%cplwav) then
      nvar2m = nvar2m + 1
    endif
    if (Model%nstf_name(1) > 0) then
      nvar2o = 18
    else
      nvar2o = 0
    endif
    if (Model%lsm == Model%lsm_ruc .and. warm_start) then
      if (Model%rdlai) then
        nvar2r = 13
      else
        nvar2r = 12
      endif
      nvar3  = 5
    else
      if(for_write .and. Model%rdlai) then
        nvar2r = 1
      else
        nvar2r = 0
      endif
      nvar3  = 3
    endif
    if (Model%lsm == Model%lsm_noahmp) then
      nvar2mp = 29
      nvar3mp = 5
    else
      nvar2mp = 0
      nvar3mp = 0
    endif
    !CLM Lake and Flake
    if (Model%lkm > 0 .and. Model%iopt_lake==Model%iopt_lake_flake) then
      nvar2l = 10
    else
      nvar2l = 0
    endif

    nvar_before_lake=nvar2m+nvar2o+nvar2r+nvar2mp

    Sfc_io_calculate_indices = &
         nvar2m /= sfc%nvar2m .or. &
         nvar2o /= sfc%nvar2o .or. &
         nvar3 /= sfc%nvar3 .or. &
         nvar2r /= sfc%nvar2r .or. &
         nvar2mp /= sfc%nvar2mp .or. &
         nvar3mp /= sfc%nvar3mp .or. &
         nvar2l /= sfc%nvar2l .or. &
         nvar2m /= sfc%nvar2m .or. &
         nvar_before_lake /= sfc%nvar_before_lake

    sfc%nvar2m = nvar2m
    sfc%nvar2o = nvar2o
    sfc%nvar3 = nvar3
    sfc%nvar2r = nvar2r
    sfc%nvar2mp = nvar2mp
    sfc%nvar3mp = nvar3mp
    sfc%nvar2l = nvar2l
    sfc%nvar2m = nvar2m
    sfc%nvar_before_lake = nvar_before_lake

  end function Sfc_io_calculate_indices

  subroutine Sfc_io_fill_3d_names(sfc,Model,warm_start)
    implicit none
    class(Sfc_io_data_type)           :: sfc
    type(GFS_control_type),    intent(in) :: Model
    logical, intent(in) :: warm_start
    integer :: nt

    !--- names of the 3d variables to save
    if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp .or. (.not.warm_start)) then
      !--- names of the 3D variables to save
      sfc%name3(1) = 'stc'
      sfc%name3(2) = 'smc'
      sfc%name3(3) = 'slc'
      if (Model%lsm == Model%lsm_noahmp) then
        sfc%name3(4) = 'snicexy'
        sfc%name3(5) = 'snliqxy'
        sfc%name3(6) = 'tsnoxy'
        sfc%name3(7) = 'smoiseq'
        sfc%name3(8) = 'zsnsoxy'
      endif
    else if (Model%lsm == Model%lsm_ruc) then
      !--- names of the 3D variables to save
      sfc%name3(1) = 'tslb'
      sfc%name3(2) = 'smois'
      sfc%name3(3) = 'sh2o'
      sfc%name3(4) = 'smfr'
      sfc%name3(5) = 'flfr'
    end if
    sfc%name3(0) = 'tiice'
  end subroutine Sfc_io_fill_3d_names

  subroutine Sfc_io_fill_2d_names(sfc,Model,warm_start)
    implicit none
    class(Sfc_io_data_type)           :: sfc
    type(GFS_control_type),    intent(in) :: Model
    logical, intent(in) :: warm_start
    integer :: nt

    !--- names of the 2D variables to save
    nt=0
    nt=nt+1 ; sfc%name2(nt) = 'slmsk'
    nt=nt+1 ; sfc%name2(nt) = 'tsea'    !tsfc
    nt=nt+1 ; sfc%name2(nt) = 'sheleg'  !weasd
    nt=nt+1 ; sfc%name2(nt) = 'tg3'
    nt=nt+1 ; sfc%name2(nt) = 'zorl'
    nt=nt+1 ; sfc%name2(nt) = 'alvsf'
    nt=nt+1 ; sfc%name2(nt) = 'alvwf'
    nt=nt+1 ; sfc%name2(nt) = 'alnsf'
    nt=nt+1 ; sfc%name2(nt) = 'alnwf'
    nt=nt+1 ; sfc%name2(nt) = 'facsf'
    nt=nt+1 ; sfc%name2(nt) = 'facwf'
    nt=nt+1 ; sfc%name2(nt) = 'vfrac'
    nt=nt+1 ; sfc%name2(nt) = 'canopy'
    nt=nt+1 ; sfc%name2(nt) = 'f10m'
    nt=nt+1 ; sfc%name2(nt) = 't2m'
    nt=nt+1 ; sfc%name2(nt) = 'q2m'
    nt=nt+1 ; sfc%name2(nt) = 'vtype'
    nt=nt+1 ; sfc%name2(nt) = 'stype'
    nt=nt+1 ; sfc%name2(nt) = 'uustar'
    nt=nt+1 ; sfc%name2(nt) = 'ffmm'
    nt=nt+1 ; sfc%name2(nt) = 'ffhh'
    nt=nt+1 ; sfc%name2(nt) = 'hice'
    nt=nt+1 ; sfc%name2(nt) = 'fice'
    nt=nt+1 ; sfc%name2(nt) = 'tisfc'
    nt=nt+1 ; sfc%name2(nt) = 'tprcp'
    nt=nt+1 ; sfc%name2(nt) = 'srflag'
    nt=nt+1 ; sfc%name2(nt) = 'snwdph'  !snowd
    nt=nt+1 ; sfc%name2(nt) = 'shdmin'
    nt=nt+1 ; sfc%name2(nt) = 'shdmax'
    nt=nt+1 ; sfc%name2(nt) = 'slope'
    nt=nt+1 ; sfc%name2(nt) = 'snoalb'
    !--- variables below here are optional
    nt=nt+1 ; sfc%name2(nt) = 'sncovr'
    nt=nt+1 ; sfc%name2(nt) = 'snodl' !snowd on land portion of a cell
    nt=nt+1 ; sfc%name2(nt) = 'weasdl'!weasd on land portion of a cell
    nt=nt+1 ; sfc%name2(nt) = 'tsfc'  !tsfc composite
    nt=nt+1 ; sfc%name2(nt) = 'tsfcl' !temp on land portion of a cell
    nt=nt+1 ; sfc%name2(nt) = 'zorlw' !zorl on water portion of a cell
    nt=nt+1 ; sfc%name2(nt) = 'zorll' !zorl on land portion of a cell
    nt=nt+1 ; sfc%name2(nt) = 'zorli' !zorl on ice portion of a cell
    nt=nt+1 ; sfc%name2(nt) = 'albdirvis_lnd'
    nt=nt+1 ; sfc%name2(nt) = 'albdirnir_lnd'
    nt=nt+1 ; sfc%name2(nt) = 'albdifvis_lnd'
    nt=nt+1 ; sfc%name2(nt) = 'albdifnir_lnd'
    nt=nt+1 ; sfc%name2(nt) = 'emis_lnd'
    nt=nt+1 ; sfc%name2(nt) = 'emis_ice'
    nt=nt+1 ; sfc%name2(nt) = 'sncovr_ice'
    nt=nt+1 ; sfc%name2(nt) = 'snodi' ! snowd on ice portion of a cell
    nt=nt+1 ; sfc%name2(nt) = 'weasdi'! weasd on ice portion of a cell

    if (Model%use_cice_alb .or. Model%lsm == Model%lsm_ruc) then
      nt=nt+1 ; sfc%name2(nt) = 'albdirvis_ice'
      nt=nt+1 ; sfc%name2(nt) = 'albdifvis_ice'
      nt=nt+1 ; sfc%name2(nt) = 'albdirnir_ice'
      nt=nt+1 ; sfc%name2(nt) = 'albdifnir_ice'
    endif

    if(Model%cplwav) then
      nt=nt+1 ; sfc%name2(nt) = 'zorlwav' !zorl from wave component
      sfc%nvar2m = nt
    endif

    if (Model%nstf_name(1) > 0) then
      !--- NSSTM inputs only needed when (nstf_name(1) > 0) .and. (nstf_name(2)) == 0)
      nt=nt+1 ; sfc%name2(nt) = 'tref'
      nt=nt+1 ; sfc%name2(nt) = 'z_c'
      nt=nt+1 ; sfc%name2(nt) = 'c_0'
      nt=nt+1 ; sfc%name2(nt) = 'c_d'
      nt=nt+1 ; sfc%name2(nt) = 'w_0'
      nt=nt+1 ; sfc%name2(nt) = 'w_d'
      nt=nt+1 ; sfc%name2(nt) = 'xt'
      nt=nt+1 ; sfc%name2(nt) = 'xs'
      nt=nt+1 ; sfc%name2(nt) = 'xu'
      nt=nt+1 ; sfc%name2(nt) = 'xv'
      nt=nt+1 ; sfc%name2(nt) = 'xz'
      nt=nt+1 ; sfc%name2(nt) = 'zm'
      nt=nt+1 ; sfc%name2(nt) = 'xtts'
      nt=nt+1 ; sfc%name2(nt) = 'xzts'
      nt=nt+1 ; sfc%name2(nt) = 'd_conv'
      nt=nt+1 ; sfc%name2(nt) = 'ifd'
      nt=nt+1 ; sfc%name2(nt) = 'dt_cool'
      nt=nt+1 ; sfc%name2(nt) = 'qrain'
    endif
    !
    ! Only needed when Noah MP LSM is used - 29 2D
    !
    if (Model%lsm == Model%lsm_noahmp) then
      nt=nt+1 ; sfc%name2(nt) = 'snowxy'
      nt=nt+1 ; sfc%name2(nt) = 'tvxy'
      nt=nt+1 ; sfc%name2(nt) = 'tgxy'
      nt=nt+1 ; sfc%name2(nt) = 'canicexy'
      nt=nt+1 ; sfc%name2(nt) = 'canliqxy'
      nt=nt+1 ; sfc%name2(nt) = 'eahxy'
      nt=nt+1 ; sfc%name2(nt) = 'tahxy'
      nt=nt+1 ; sfc%name2(nt) = 'cmxy'
      nt=nt+1 ; sfc%name2(nt) = 'chxy'
      nt=nt+1 ; sfc%name2(nt) = 'fwetxy'
      nt=nt+1 ; sfc%name2(nt) = 'sneqvoxy'
      nt=nt+1 ; sfc%name2(nt) = 'alboldxy'
      nt=nt+1 ; sfc%name2(nt) = 'qsnowxy'
      nt=nt+1 ; sfc%name2(nt) = 'wslakexy'
      nt=nt+1 ; sfc%name2(nt) = 'zwtxy'
      nt=nt+1 ; sfc%name2(nt) = 'waxy'
      nt=nt+1 ; sfc%name2(nt) = 'wtxy'
      nt=nt+1 ; sfc%name2(nt) = 'lfmassxy'
      nt=nt+1 ; sfc%name2(nt) = 'rtmassxy'
      nt=nt+1 ; sfc%name2(nt) = 'stmassxy'
      nt=nt+1 ; sfc%name2(nt) = 'woodxy'
      nt=nt+1 ; sfc%name2(nt) = 'stblcpxy'
      nt=nt+1 ; sfc%name2(nt) = 'fastcpxy'
      nt=nt+1 ; sfc%name2(nt) = 'xsaixy'
      nt=nt+1 ; sfc%name2(nt) = 'xlaixy'
      nt=nt+1 ; sfc%name2(nt) = 'taussxy'
      nt=nt+1 ; sfc%name2(nt) = 'smcwtdxy'
      nt=nt+1 ; sfc%name2(nt) = 'deeprechxy'
      nt=nt+1 ; sfc%name2(nt) = 'rechxy'
    else if (Model%lsm == Model%lsm_ruc .and. warm_start) then
      nt=nt+1 ; sfc%name2(nt) = 'wetness'
      nt=nt+1 ; sfc%name2(nt) = 'clw_surf_land'
      nt=nt+1 ; sfc%name2(nt) = 'clw_surf_ice'
      nt=nt+1 ; sfc%name2(nt) = 'qwv_surf_land'
      nt=nt+1 ; sfc%name2(nt) = 'qwv_surf_ice'
      nt=nt+1 ; sfc%name2(nt) = 'tsnow_land'
      nt=nt+1 ; sfc%name2(nt) = 'tsnow_ice'
      nt=nt+1 ; sfc%name2(nt) = 'snowfall_acc_land'
      nt=nt+1 ; sfc%name2(nt) = 'snowfall_acc_ice'
      nt=nt+1 ; sfc%name2(nt) = 'sfalb_lnd'
      nt=nt+1 ; sfc%name2(nt) = 'sfalb_lnd_bck'
      nt=nt+1 ; sfc%name2(nt) = 'sfalb_ice'
      if (Model%rdlai) then
        nt=nt+1 ; sfc%name2(nt) = 'lai'
      endif
    else if (Model%lsm == Model%lsm_ruc .and. Model%rdlai) then
      nt=nt+1 ; sfc%name2(nt) = 'lai'
    endif

    if (Model%lkm > 0 .and. Model%iopt_lake==Model%iopt_lake_flake) then
      nt=nt+1 ; sfc%name2(nt) = 'T_snow'
      nt=nt+1 ; sfc%name2(nt) = 'T_ice'
      nt=nt+1 ; sfc%name2(nt) = 'h_ML'
      nt=nt+1 ; sfc%name2(nt) = 't_ML'
      nt=nt+1 ; sfc%name2(nt) = 't_mnw'
      nt=nt+1 ; sfc%name2(nt) = 'h_talb'
      nt=nt+1 ; sfc%name2(nt) = 't_talb'
      nt=nt+1 ; sfc%name2(nt) = 't_bot1'
      nt=nt+1 ; sfc%name2(nt) = 't_bot2'
      nt=nt+1 ; sfc%name2(nt) = 'c_t'
    endif
  end subroutine Sfc_io_fill_2d_names

  subroutine Sfc_io_final(sfc)
    implicit none
    type(Sfc_io_data_type)             :: sfc

    sfc%nvar2m=0
    sfc%nvar2o=0
    sfc%nvar2l=0
    sfc%nvar3=0
    sfc%nvar2r=0
    sfc%nvar2mp=0
    sfc%nvar3mp=0
    sfc%nvar2m=0
    sfc%nvar_before_lake=0
    sfc%is_lsoil=.false.

    ! This #define reduces code length by a lot
#define IF_ASSOC_DEALLOC_NULL(var) \
    if(associated(sfc%var)) then ; \
      deallocate(sfc%var) ; \
      nullify(sfc%var) ; \
    endif

    IF_ASSOC_DEALLOC_NULL(var2)
    IF_ASSOC_DEALLOC_NULL(var3ice)
    IF_ASSOC_DEALLOC_NULL(var3)
    IF_ASSOC_DEALLOC_NULL(var3sn)
    IF_ASSOC_DEALLOC_NULL(var3eq)
    IF_ASSOC_DEALLOC_NULL(var3zn)
    IF_ASSOC_DEALLOC_NULL(name2)
    IF_ASSOC_DEALLOC_NULL(name3)

#undef IF_ASSOC_DEALLOC_NULL

  end subroutine Sfc_io_final
end module FV3GFS_sfc_io
