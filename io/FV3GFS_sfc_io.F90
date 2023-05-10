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
  public :: Sfc_io_fill_2d_names, Sfc_io_fill_3d_names, Sfc_io_allocate_arrays, &
       Sfc_io_register_axes, Sfc_io_write_axes, Sfc_io_register_2d_fields, &
       Sfc_io_register_3d_fields, Sfc_io_copy_to_grid, Sfc_io_copy_from_grid, &
       Sfc_io_apply_safeguards, Sfc_io_final

  real(kind=kind_phys), parameter :: timin = 173.0_kind_phys  ! minimum temperature allowed for snow/ice
  real(kind_phys), parameter:: min_lake_orog = 200.0_kind_phys
  real(kind_phys), parameter:: zero = 0, one = 1

  type Sfc_io_data_type
    integer, private :: nvar2o = 0
    integer, private :: nvar3 = 0
    integer, private :: nvar2r = 0
    integer, private :: nvar2mp = 0
    integer, private :: nvar3mp = 0
    integer, private :: nvar2l = 0
    integer, private :: nvar2m = 0
    integer, private :: nvar_before_lake = 0

    ! The lsoil flag is only meaningful when reading:;
    logical, private :: is_lsoil = .false.

    ! SYNONYMS: Some nvar variables had two names in FV3GFS_io.F90. They have
    ! only one name here. The "_s" is redundant because this file only has
    ! surface restart variables.
    !
    !  - nvar2m = nvar_s2m
    !  - nvar2o = nvar_s2o
    !  - nvar2r = nvar_s2r
    !  - nvar2mp = nvar_s2mp
    !  - nvar3mp = nvar_s3mp

    real(kind=kind_phys), pointer, dimension(:,:,:), private :: var2 => null()
    real(kind=kind_phys), pointer, dimension(:,:,:), private :: var3ice => null()
    real(kind=kind_phys), pointer, dimension(:,:,:,:), private :: var3 => null()
    real(kind=kind_phys), pointer, dimension(:,:,:,:), private :: var3sn => null()
    real(kind=kind_phys), pointer, dimension(:,:,:,:), private :: var3eq => null()
    real(kind=kind_phys), pointer, dimension(:,:,:,:), private :: var3zn => null()

    character(len=32), pointer, dimension(:), private :: name2 => null()
    character(len=32), pointer, dimension(:), private :: name3 => null()

  contains

    procedure, public :: allocate_arrays => Sfc_io_allocate_arrays
    procedure, public :: register_axes => Sfc_io_register_axes
    procedure, public :: write_axes => Sfc_io_write_axes
    procedure, public :: register_2d_fields => Sfc_io_register_2d_fields
    procedure, public :: register_3d_fields => Sfc_io_register_3d_fields
    procedure, public :: fill_2d_names => Sfc_io_fill_2d_names
    procedure, public :: fill_3d_names => Sfc_io_fill_3d_names
    procedure, public :: init_fields => Sfc_io_init_fields
    procedure, public :: copy_to_grid => Sfc_io_copy_to_grid
    procedure, public :: copy_from_grid => Sfc_io_copy_from_grid
    procedure, public :: apply_safeguards => Sfc_io_apply_safeguards

    procedure, private :: calculate_indices => Sfc_io_calculate_indices

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

  function Sfc_io_allocate_arrays(sfc, Model, Atm_block, for_write, warm_start)
    implicit none
    class(Sfc_io_data_type)             :: sfc
    type(block_control_type),    intent(in) :: Atm_block
    type(GFS_control_type),      intent(in) :: Model
    logical :: Sfc_io_allocate_arrays
    logical, intent(in) :: for_write, warm_start

    integer :: isc, iec, jsc, jec, npz, nx, ny

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)

    Sfc_io_allocate_arrays = sfc%calculate_indices(Model, for_write, warm_start)
    Sfc_io_allocate_arrays = Sfc_io_allocate_arrays .or. .not. associated(sfc%name2)

    if(Sfc_io_allocate_arrays) then
      !--- allocate the various containers needed for restarts
      allocate(sfc%name2(sfc%nvar2m+sfc%nvar2o+sfc%nvar2mp+sfc%nvar2r+sfc%nvar2l))
      allocate(sfc%name3(0:sfc%nvar3+sfc%nvar3mp))
      allocate(sfc%var2(nx,ny,sfc%nvar2m+sfc%nvar2o+sfc%nvar2mp+sfc%nvar2r+sfc%nvar2l))

      ! Note that this may cause problems with RUC LSM for coldstart runs from GFS data
      ! if the initial conditions do contain this variable, because Model%kice is 9 for
      ! RUC LSM, but tiice in the initial conditions will only have two vertical layers
      allocate(sfc%var3ice(nx,ny,Model%kice))

      if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp .or. (.not.warm_start)) then
        allocate(sfc%var3(nx,ny,Model%lsoil,sfc%nvar3))
      elseif (Model%lsm == Model%lsm_ruc) then
        allocate(sfc%var3(nx,ny,Model%lsoil_lsm,sfc%nvar3))
      endif

      sfc%var2   = -9999.0_kind_phys
      sfc%var3   = -9999.0_kind_phys
      sfc%var3ice= -9999.0_kind_phys

      if (Model%lsm == Model%lsm_noahmp) then
        allocate(sfc%var3sn(nx,ny,-2:0,4:6))
        allocate(sfc%var3eq(nx,ny,1:4,7:7))
        allocate(sfc%var3zn(nx,ny,-2:4,8:8))

        sfc%var3sn = -9999.0_kind_phys
        sfc%var3eq = -9999.0_kind_phys
        sfc%var3zn = -9999.0_kind_phys
      endif
    endif

    if(Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_flake .and. Model%me==0) then
      if(size(sfc%name2)/=sfc%nvar_before_lake+10) then
3814    format("ERROR: size mismatch size(sfc%name2)=",I0," /= nvar_before_lake+10=",I0)
        write(0,3814) size(sfc%name2),sfc%nvar_before_lake+10
      endif
    endif
  end function Sfc_io_allocate_arrays

  subroutine Sfc_io_register_axes(sfc, Model, Sfc_restart, for_write, warm_start)
    implicit none
    class(Sfc_io_data_type)             :: sfc
    type(GFS_control_type),      intent(in) :: Model
    type(FmsNetcdfDomainFile_t) :: Sfc_restart
    logical, intent(in) :: for_write, warm_start

    if(.not.for_write) then
      sfc%is_lsoil = .false.
    endif

    if(.not.warm_start .and. .not.for_write) then
      if( variable_exists(Sfc_restart,"lsoil") ) then
        if(.not.for_write) then
          sfc%is_lsoil=.true.
        endif
        call register_axis(Sfc_restart, 'lon', 'X')
        call register_axis(Sfc_restart, 'lat', 'Y')
        call register_axis(Sfc_restart, 'lsoil', dimension_length=Model%lsoil)
      else
        call register_axis(Sfc_restart, 'xaxis_1', 'X')
        call register_axis(Sfc_restart, 'yaxis_1', 'Y')
        call register_axis(Sfc_restart, 'zaxis_1', dimension_length=4)
        call register_axis(Sfc_restart, 'Time', 1)
      end if
    else
      call register_axis(Sfc_restart, 'xaxis_1', 'X')
      call register_axis(Sfc_restart, 'yaxis_1', 'Y')
      call register_axis(Sfc_restart, 'zaxis_1', dimension_length=Model%kice)
      if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp) then
        call register_axis(Sfc_restart, 'zaxis_2', dimension_length=Model%lsoil)
      else if(Model%lsm == Model%lsm_ruc .and. .not.for_write) then
        ! Possible bug. This is only defined on read, not write.
        call register_axis(Sfc_restart, 'zaxis_2', dimension_length=Model%lsoil_lsm)
      endif
      if(Model%lsm == Model%lsm_noahmp) then
        call register_axis(Sfc_restart, 'zaxis_3', dimension_length=3)
        call register_axis(Sfc_restart, 'zaxis_4', dimension_length=7)
      end if
      call register_axis(Sfc_restart, 'Time', unlimited)
    endif
  end subroutine Sfc_io_register_axes

  subroutine Sfc_io_write_axes(sfc, Model, Sfc_restart)
    implicit none
    class(Sfc_io_data_type)             :: sfc
    type(GFS_control_type),      intent(in) :: Model
    type(FmsNetcdfDomainFile_t) :: Sfc_restart

    integer, allocatable :: buffer(:)
    integer :: i, is, ie
    logical :: mand

    call register_field(Sfc_restart, 'xaxis_1', 'double', (/'xaxis_1'/))
    call register_variable_attribute(Sfc_restart, 'xaxis_1', 'cartesian_axis', 'X', str_len=1)
    call get_global_io_domain_indices(Sfc_restart, 'xaxis_1', is, ie, indices=buffer)
    call write_data(Sfc_restart, "xaxis_1", buffer)
    deallocate(buffer)

    call register_field(Sfc_restart, 'yaxis_1', 'double', (/'yaxis_1'/))
    call register_variable_attribute(Sfc_restart, 'yaxis_1', 'cartesian_axis', 'Y', str_len=1)
    call get_global_io_domain_indices(Sfc_restart, 'yaxis_1', is, ie, indices=buffer)
    call write_data(Sfc_restart, "yaxis_1", buffer)
    deallocate(buffer)

    call register_field(Sfc_restart, 'zaxis_1', 'double', (/'zaxis_1'/))
    call register_variable_attribute(Sfc_restart, 'zaxis_1', 'cartesian_axis', 'Z', str_len=1)
    allocate( buffer(Model%kice) )
    do i=1, Model%kice
      buffer(i) = i
    end do
    call write_data(Sfc_restart, 'zaxis_1', buffer)
    deallocate(buffer)

    if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp) then
      call register_field(Sfc_restart, 'zaxis_2', 'double', (/'zaxis_2'/))
      call register_variable_attribute(Sfc_restart, 'zaxis_2', 'cartesian_axis', 'Z', str_len=1)
      allocate( buffer(Model%lsoil) )
      do i=1, Model%lsoil
        buffer(i)=i
      end do
      call write_data(Sfc_restart, 'zaxis_2', buffer)
      deallocate(buffer)
    endif

    if(Model%lsm == Model%lsm_noahmp) then
      call register_field(Sfc_restart, 'zaxis_3', 'double', (/'zaxis_3'/))
      call register_variable_attribute(Sfc_restart, 'zaxis_3', 'cartesian_axis', 'Z', str_len=1)
      allocate(buffer(3))
      do i=1, 3
        buffer(i) = i
      end do
      call write_data(Sfc_restart, 'zaxis_3', buffer)
      deallocate(buffer)

      call register_field(Sfc_restart, 'zaxis_4', 'double', (/'zaxis_4'/))
      call register_variable_attribute(Sfc_restart, 'zaxis_4', 'cartesian_axis' ,'Z', str_len=1)
      allocate(buffer(7))
      do i=1, 7
        buffer(i)=i
      end do
      call write_data(Sfc_restart, 'zaxis_4', buffer)
      deallocate(buffer)
    end if
    call register_field(Sfc_restart, 'Time', 'double', (/'Time'/))
    call register_variable_attribute(Sfc_restart, 'Time', 'cartesian_axis', 'T', str_len=1)
    call write_data( Sfc_restart, 'Time', 1)
  end subroutine Sfc_io_write_axes

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

  subroutine Sfc_io_register_2d_fields(sfc,Model,Sfc_restart,for_write,warm_start)
    implicit none
    class(Sfc_io_data_type)             :: sfc
    type(GFS_control_type),      intent(in) :: Model
    type(FmsNetcdfDomainFile_t) :: Sfc_restart
    logical, intent(in) :: for_write, warm_start

    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p  => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p  => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p1 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p2 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p3 => NULL()
    integer :: num
    logical :: mand

    character(len=7) :: time2d(3)

    if(for_write) then
      time2d=(/'xaxis_1','yaxis_1','Time   '/)
    else
      time2d=(/'Time   ','yaxis_1','xaxis_1'/)
    endif

    !--- register the 2D fields
    do num = 1,sfc%nvar2m
      var2_p => sfc%var2(:,:,num)
      if (trim(sfc%name2(num)) == 'sncovr' .or. trim(sfc%name2(num)) == 'tsfcl' .or. trim(sfc%name2(num)) == 'zorll' &
           .or. trim(sfc%name2(num)) == 'zorli' .or. trim(sfc%name2(num)) == 'zorlwav' &
           .or. trim(sfc%name2(num)) == 'snodl' .or. trim(sfc%name2(num)) == 'weasdl'  &
           .or. trim(sfc%name2(num)) == 'snodi' .or. trim(sfc%name2(num)) == 'weasdi'  &
           .or. trim(sfc%name2(num)) == 'tsfc'  .or. trim(sfc%name2(num)) ==  'zorlw'  &
           .or. trim(sfc%name2(num)) == 'albdirvis_lnd' .or. trim(sfc%name2(num)) == 'albdirnir_lnd' &
           .or. trim(sfc%name2(num)) == 'albdifvis_lnd' .or. trim(sfc%name2(num)) == 'albdifnir_lnd' &
           .or. trim(sfc%name2(num)) == 'albdirvis_ice' .or. trim(sfc%name2(num)) == 'albdirnir_ice' &
           .or. trim(sfc%name2(num)) == 'albdifvis_ice' .or. trim(sfc%name2(num)) == 'albdifnir_ice' &
           .or. trim(sfc%name2(num)) == 'emis_lnd'      .or. trim(sfc%name2(num)) == 'emis_ice'      &
           .or. trim(sfc%name2(num)) == 'sncovr_ice') then
        if(.not.for_write .and. sfc%is_lsoil) then
          call register_restart_field(Sfc_restart, sfc%name2(num), var2_p, dimensions=(/'lat','lon'/), is_optional=.true.)
        else
          call register_restart_field(Sfc_restart, sfc%name2(num), var2_p, dimensions=time2d,&
               &is_optional=.true.)
        end if
      else
        if(.not.for_write .and. sfc%is_lsoil) then
          call register_restart_field(Sfc_restart,sfc%name2(num),var2_p, dimensions=(/'lat','lon'/))
        else
          call register_restart_field(Sfc_restart,sfc%name2(num),var2_p, dimensions=time2d)
        end if
      endif
    enddo

    if (Model%nstf_name(1) > 0) then
      mand = .false.
      if (Model%nstf_name(2) == 0) mand = .true.
      do num = sfc%nvar2m+1,sfc%nvar2m+sfc%nvar2o
        var2_p => sfc%var2(:,:,num)
        if(sfc%is_lsoil) then
          call register_restart_field(Sfc_restart, sfc%name2(num), var2_p, dimensions=(/'lat','lon'/), is_optional=.not.mand)
        else
          call register_restart_field(Sfc_restart, sfc%name2(num), var2_p, dimensions=time2d, is_optional=.not.mand)
        endif
      enddo
    endif

    if (Model%lsm == Model%lsm_ruc) then ! sfc%nvar2mp = 0
      do num = sfc%nvar2m+sfc%nvar2o+1, sfc%nvar2m+sfc%nvar2o+sfc%nvar2r
        var2_p => sfc%var2(:,:,num)
        if(sfc%is_lsoil) then
          call register_restart_field(Sfc_restart, sfc%name2(num), var2_p, dimensions=(/'lat','lon'/) )
        else
          call register_restart_field(Sfc_restart, sfc%name2(num), var2_p, dimensions=time2d)
        end if
      enddo
    endif ! mp/ruc

    ! Noah MP register only necessary only lsm = 2, not necessary has values
    if ( (for_write .and. Model%lsm == Model%lsm_noahmp) &
         .or. (.not. for_write .and. sfc%nvar2mp > 0) ) then
      mand = for_write
      do num = sfc%nvar2m+sfc%nvar2o+1,sfc%nvar2m+sfc%nvar2o+sfc%nvar2mp
        var2_p => sfc%var2(:,:,num)
        if(sfc%is_lsoil) then
          call register_restart_field(Sfc_restart, sfc%name2(num), var2_p, dimensions=(/'lat','lon'/), is_optional=.not.mand)
        else
          call register_restart_field(Sfc_restart, sfc%name2(num), var2_p, dimensions=time2d, is_optional=.not.mand)
        end if
      enddo
    endif ! noahmp

    ! Flake
    if (Model%lkm > 0 .and. Model%iopt_lake==Model%iopt_lake_flake) then
      mand = for_write
      do num = sfc%nvar_before_lake+1,sfc%nvar_before_lake+sfc%nvar2l
        var2_p => sfc%var2(:,:,num)
        if(sfc%is_lsoil) then
          call register_restart_field(Sfc_restart, sfc%name2(num),var2_p,dimensions=(/'lat','lon'/), is_optional=.not.mand) 
        else
          call register_restart_field(Sfc_restart, sfc%name2(num),var2_p,dimensions=time2d, is_optional=.not.mand)
        endif
      enddo
    endif

  end subroutine Sfc_io_register_2d_fields

  subroutine Sfc_io_register_3d_fields(sfc,Model,Sfc_restart,for_write,warm_start)
    implicit none
    class(Sfc_io_data_type)             :: sfc
    type(GFS_control_type),      intent(in) :: Model
    type(FmsNetcdfDomainFile_t) :: Sfc_restart
    logical, intent(in) :: for_write, warm_start

    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p  => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p1 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p2 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p3 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_fr => NULL()
    integer :: num
    logical :: mand

    character(len=7), parameter :: xyz1_time(4) = (/'xaxis_1', 'yaxis_1', 'zaxis_1', 'Time   '/)
    character(len=7), parameter :: xyz2_time(4) = (/'xaxis_1', 'yaxis_1', 'zaxis_2', 'Time   '/)
    character(len=7), parameter :: xyz3_time(4) = (/'xaxis_1', 'yaxis_1', 'zaxis_3', 'Time   '/)
    character(len=7), parameter :: xyz4_time(4) = (/'xaxis_1', 'yaxis_1', 'zaxis_4', 'Time   '/)

    !--- register the 3D fields
    var3_p => sfc%var3ice(:,:,:)
    call register_restart_field(Sfc_restart, sfc%name3(0), var3_p, dimensions=xyz1_time, is_optional=.true.)

    if(.not. for_write) then
      do num = 1,sfc%nvar3
        var3_p => sfc%var3(:,:,:,num)
        if ( warm_start ) then
          call register_restart_field(Sfc_restart, sfc%name3(num), var3_p, dimensions=(/'xaxis_1', 'yaxis_1', 'lsoil  ', 'Time   '/),&
               &is_optional=.true.)
        else
          if(sfc%is_lsoil) then
            call register_restart_field(Sfc_restart, sfc%name3(num), var3_p, dimensions=(/'lat  ', 'lon  ', 'lsoil'/), is_optional=.true.)
          else
            call register_restart_field(Sfc_restart, sfc%name3(num), var3_p, dimensions=xyz2_time,&
                 &is_optional=.true.)
          end if
        end if
      enddo
    elseif(Model%lsm == Model%lsm_ruc) then
      do num = 1,sfc%nvar3
        var3_p => sfc%var3(:,:,:,num)
        call register_restart_field(Sfc_restart, sfc%name3(num), var3_p, dimensions=xyz1_time)
      enddo
      nullify(var3_p)
    else ! writing something other than ruc
      do num = 1,sfc%nvar3
        var3_p => sfc%var3(:,:,:,num)
        call register_restart_field(Sfc_restart, sfc%name3(num), var3_p, dimensions=xyz2_time)
      enddo
      nullify(var3_p)
    endif

    if (Model%lsm == Model%lsm_noahmp) then
      mand = for_write
      do num = sfc%nvar3+1,sfc%nvar3+3
        var3_p1 => sfc%var3sn(:,:,:,num)
        call register_restart_field(Sfc_restart, sfc%name3(num), var3_p1, dimensions=xyz3_time, is_optional=.not.mand)
      enddo

      var3_p2 => sfc%var3eq(:,:,:,7)
      call register_restart_field(Sfc_restart, sfc%name3(7), var3_p2, dimensions=xyz2_time, is_optional=.not.mand)

      var3_p3 => sfc%var3zn(:,:,:,8)
      call register_restart_field(Sfc_restart, sfc%name3(8), var3_p3, dimensions=xyz4_time, is_optional=.not.mand)
    endif   !mp

  end subroutine Sfc_io_register_3d_fields

  subroutine Sfc_io_init_fields(sfc,Model)
    implicit none
    class(Sfc_io_data_type)           :: sfc
    type(GFS_control_type),    intent(in) :: Model

    !--- Noah MP define arbitrary value (number layers of snow) to indicate
    !coldstart(sfcfile doesn't include noah mp fields) or not

    if (Model%lsm == Model%lsm_noahmp) then
      sfc%var2(1,1,sfc%nvar2m+19) = -66666.0_kind_phys
    endif
  end subroutine Sfc_io_init_fields

  subroutine Sfc_io_copy_to_grid(sfc, Model, Atm_block, Sfcprop, warm_start, override_frac_grid)
    !--- interface variable definitions
    implicit none

    class(Sfc_io_data_type)             :: sfc
    type(GFS_sfcprop_type),      intent(in) :: Sfcprop(:)
    type(block_control_type),    intent(in) :: Atm_block
    type(GFS_control_type),      intent(in) :: Model
    logical, intent(in) :: warm_start
    logical, intent(out), optional :: override_frac_grid

    integer :: i, j, k, nb, ix, lsoil, num, nt
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer, allocatable :: ii1(:), jj1(:)
    real(kind_phys) :: ice

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)

    !   write(0,*)' stype read in min,max=',minval(sfc_var2(:,:,35)),maxval(sfc_var2(:,:,35)),' sfc_name2=',sfc_name2(35)
    !   write(0,*)' stype read in min,max=',minval(sfc_var2(:,:,18)),maxval(sfc_var2(:,:,18))
    !   write(0,*)' sfc_var2=',sfc_var2(:,:,12)

    !$omp parallel do default(shared) private(i, j, nb, ix, nt, ii1, jj1, lsoil)
    block_loop: do nb = 1, Atm_block%nblks
      allocate(ii1(Atm_block%blksz(nb)))
      allocate(jj1(Atm_block%blksz(nb)))
      ii1=Atm_block%index(nb)%ii - isc + 1
      jj1=Atm_block%index(nb)%jj - jsc + 1

      nt=0

      !--- 2D variables
      !    ------------
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%slmsk)   !--- slmsk
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tsfco)   !--- tsfc (tsea in sfc file)
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%weasd)   !--- weasd (sheleg in sfc file)
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tg3)     !--- tg3
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%zorl)    !--- zorl composite
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%alvsf)   !--- alvsf
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%alvwf)   !--- alvwf
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%alnsf)   !--- alnsf
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%alnwf)   !--- alnwf
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%facsf)   !--- facsf
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%facwf)   !--- facwf
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%vfrac)   !--- vfrac
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%canopy)  !--- canopy
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%f10m)    !--- f10m
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%t2m)     !--- t2m
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%q2m)     !--- q2m
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%vtype)   !--- vtype
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%stype)   !--- stype
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%uustar)  !--- uustar
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%ffmm)    !--- ffmm
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%ffhh)    !--- ffhh
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%hice)    !--- hice
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%fice)    !--- fice
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tisfc)   !--- tisfc
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tprcp)   !--- tprcp
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%srflag)  !--- srflag
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%snowd)   !--- snowd (snwdph in the file)
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%shdmin)  !--- shdmin
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%shdmax)  !--- shdmax
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%slope)   !--- slope
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%snoalb)  !--- snoalb
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%sncovr)  !--- sncovr
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%snodl)   !--- snodl (snowd on land  portion of a cell)
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%weasdl)  !--- weasdl (weasd on land  portion of a cell)
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tsfc)    !--- tsfc composite
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tsfcl)   !--- tsfcl  (temp on land portion of a cell)
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%zorlw)   !--- zorlw (zorl on water portion of a cell)
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%zorll)   !--- zorll (zorl on land portion of a cell)
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%zorli)   !--- zorli (zorl on ice  portion of a cell)
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%albdirvis_lnd)
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%albdirnir_lnd)
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%albdifvis_lnd)
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%albdifnir_lnd)
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%emis_lnd)
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%emis_ice)
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%sncovr_ice)
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%snodi)   !--- snodi (snowd on ice  portion of a cell)
      call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%weasdi)  !--- weasdi (weasd on ice  portion of a cell)
      if (Model%use_cice_alb .or. Model%lsm == Model%lsm_ruc) then
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%albdirvis_ice)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%albdifvis_ice)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%albdirnir_ice)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%albdifnir_ice)
        !         call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%sfalb_ice)
      endif
      if(Model%cplwav) then
        !tgs - the following line is a bug. It should be nt = nt
        !nt = sfc%nvar2m-1 ! Next item will be at sfc%nvar2m
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%zorlwav) !--- (zorl from wave model)
      else
        Sfcprop(nb)%zorlwav  = Sfcprop(nb)%zorlw
      endif

      if(present(override_frac_grid)) then
        override_frac_grid=Model%frac_grid
      endif

      do_lsi_fractions: do ix = 1, Atm_block%blksz(nb)
        if (Sfcprop(nb)%stype(ix) == 14 .or. Sfcprop(nb)%stype(ix) <= 0) then
          Sfcprop(nb)%landfrac(ix) = zero
          Sfcprop(nb)%stype(ix) = 0
          if (Sfcprop(nb)%lakefrac(ix) > zero) then
            Sfcprop(nb)%lakefrac(ix) = one
          endif
        endif

        if_frac_grid: if (Model%frac_grid) then
          if (Sfcprop(nb)%landfrac(ix) > -999.0_kind_phys) then
            Sfcprop(nb)%slmsk(ix) = ceiling(Sfcprop(nb)%landfrac(ix)-1.0e-6)
            if (Sfcprop(nb)%slmsk(ix) == 1 .and. Sfcprop(nb)%stype(ix) == 14) &
                 Sfcprop(nb)%slmsk(ix) = 0
            if (Sfcprop(nb)%lakefrac(ix) > zero) then
              Sfcprop(nb)%oceanfrac(ix) = zero ! lake & ocean don't coexist in a cell
              if (nint(Sfcprop(nb)%slmsk(ix)) /= 1) then
                if(Sfcprop(nb)%fice(ix) >= Model%min_lakeice) then
                  Sfcprop(nb)%slmsk(ix) = 2
                else
                  Sfcprop(nb)%slmsk(ix) = 0
                endif
              endif
            else
              Sfcprop(nb)%lakefrac(ix)  = zero
              Sfcprop(nb)%oceanfrac(ix) = one - Sfcprop(nb)%landfrac(ix)
              if (nint(Sfcprop(nb)%slmsk(ix)) /= 1) then
                if (Sfcprop(nb)%fice(ix) >= Model%min_seaice) then
                  Sfcprop(nb)%slmsk(ix) = 2
                else
                  Sfcprop(nb)%slmsk(ix) = 0
                endif
              endif
            endif
          else
            if(present(override_frac_grid)) then
              override_frac_grid = .false.
            endif
            if (nint(Sfcprop(nb)%slmsk(ix)) == 1) then
              Sfcprop(nb)%landfrac(ix)  = one
              Sfcprop(nb)%lakefrac(ix)  = zero
              Sfcprop(nb)%oceanfrac(ix) = zero
            else
              if (Sfcprop(nb)%slmsk(ix) < 0.1_kind_phys .or. Sfcprop(nb)%slmsk(ix) > 1.9_kind_phys) then
                Sfcprop(nb)%landfrac(ix) = zero
                if (Sfcprop(nb)%oro_uf(ix) > min_lake_orog) then   ! lakes
                  Sfcprop(nb)%lakefrac(ix)  = one
                  Sfcprop(nb)%oceanfrac(ix) = zero
                else                                               ! ocean
                  Sfcprop(nb)%lakefrac(ix)  = zero
                  Sfcprop(nb)%oceanfrac(ix) = one
                endif
              endif
            endif
          endif
        else                                             ! not a fractional grid
          if (Sfcprop(nb)%landfrac(ix) > -999.0_kind_phys) then
            if (Sfcprop(nb)%lakefrac(ix) > zero) then
              Sfcprop(nb)%oceanfrac(ix) = zero
              Sfcprop(nb)%landfrac(ix)  = zero
              Sfcprop(nb)%lakefrac(ix)  = one
              Sfcprop(nb)%slmsk(ix)     = zero
              if (Sfcprop(nb)%fice(ix) >= Model%min_lakeice) Sfcprop(nb)%slmsk(ix) = 2.0
            else
              Sfcprop(nb)%slmsk(ix) = nint(Sfcprop(nb)%landfrac(ix))
              if (Sfcprop(nb)%stype(ix) <= 0 .or. Sfcprop(nb)%stype(ix) == 14) &
                   Sfcprop(nb)%slmsk(ix) = zero
              if (nint(Sfcprop(nb)%slmsk(ix)) == 0) then
                Sfcprop(nb)%oceanfrac(ix) = one
                Sfcprop(nb)%landfrac(ix)  = zero
                Sfcprop(nb)%lakefrac(ix)  = zero
                if (Sfcprop(nb)%fice(ix) >= Model%min_seaice) Sfcprop(nb)%slmsk(ix) = 2.0
              else
                Sfcprop(nb)%landfrac(ix)  = one
                Sfcprop(nb)%lakefrac(ix)  = zero
                Sfcprop(nb)%oceanfrac(ix) = zero
              endif
            endif
          else
            if (nint(Sfcprop(nb)%slmsk(ix)) == 1 .and. Sfcprop(nb)%stype(ix) > 0      &
                 .and. Sfcprop(nb)%stype(ix) /= 14) then
              Sfcprop(nb)%landfrac(ix)  = one
              Sfcprop(nb)%lakefrac(ix)  = zero
              Sfcprop(nb)%oceanfrac(ix) = zero
            else
              Sfcprop(nb)%slmsk(ix)    = zero
              Sfcprop(nb)%landfrac(ix) = zero
              if (Sfcprop(nb)%oro_uf(ix) > min_lake_orog) then   ! lakes
                Sfcprop(nb)%lakefrac(ix) = one
                Sfcprop(nb)%oceanfrac(ix) = zero
                if (Sfcprop(nb)%fice(ix) > Model%min_lakeice) Sfcprop(nb)%slmsk(ix) = 2.0
              else                                       ! ocean
                Sfcprop(nb)%lakefrac(ix)  = zero
                Sfcprop(nb)%oceanfrac(ix) = one
                if (Sfcprop(nb)%fice(ix) > Model%min_seaice) Sfcprop(nb)%slmsk(ix) = 2.0
              endif
            endif
          endif
        endif if_frac_grid
      enddo do_lsi_fractions

      if (warm_start .and. Model%kdt > 1) then
        do ix = 1, Atm_block%blksz(nb)
          Sfcprop(nb)%slmsk(ix)  = sfc%var2(ii1(ix),jj1(ix),1)    !--- slmsk
        enddo
      endif

      !
      !--- NSSTM variables
      !tgs - the following line is a bug that will show if(Model%cplwav) = true
      !nt = sfc%nvar2m 
      if (Model%nstf_name(1) > 0) then
        if (Model%nstf_name(2) == 1) then             ! nsst spinup
          !--- nsstm tref
          nt = nt + 18
          Sfcprop(nb)%tref    = Sfcprop(nb)%tsfco
          Sfcprop(nb)%z_c     = zero
          Sfcprop(nb)%c_0     = zero
          Sfcprop(nb)%c_d     = zero
          Sfcprop(nb)%w_0     = zero
          Sfcprop(nb)%w_d     = zero
          Sfcprop(nb)%xt      = zero
          Sfcprop(nb)%xs      = zero
          Sfcprop(nb)%xu      = zero
          Sfcprop(nb)%xv      = zero
          Sfcprop(nb)%xz      = 20.0_kind_phys
          Sfcprop(nb)%zm      = zero
          Sfcprop(nb)%xtts    = zero
          Sfcprop(nb)%xzts    = zero
          Sfcprop(nb)%d_conv  = zero
          Sfcprop(nb)%ifd     = zero
          Sfcprop(nb)%dt_cool = zero
          Sfcprop(nb)%qrain   = zero
        elseif (Model%nstf_name(2) == 0) then         ! nsst restart
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tref)  !--- nsstm tref
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%z_c)  !--- nsstm z_c
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%c_0)  !--- nsstm c_0
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%c_d)  !--- nsstm c_d
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%w_0)  !--- nsstm w_0
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%w_d)  !--- nsstm w_d
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xt)  !--- nsstm xt
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xs)  !--- nsstm xs
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xu)  !--- nsstm xu
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xv) !--- nsstm xv
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xz) !--- nsstm xz
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%zm) !--- nsstm zm
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xtts) !--- nsstm xtts
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xzts) !--- nsstm xzts
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%d_conv) !--- nsstm d_conv
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%ifd) !--- nsstm ifd
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%dt_cool) !--- nsstm dt_cool
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%qrain) !--- nsstm qrain
        endif
      endif

      if (Model%lsm == Model%lsm_ruc .and. warm_start) then
        !--- Extra RUC variables
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%wetness)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%clw_surf_land)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%clw_surf_ice)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%qwv_surf_land)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%qwv_surf_ice)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tsnow_land)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tsnow_ice)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%snowfallac_land)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%snowfallac_ice)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%sfalb_lnd)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%sfalb_lnd_bck)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%sfalb_ice)
        if (Model%rdlai) then
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xlaixy)
        endif
      else if (Model%lsm == Model%lsm_ruc) then
        ! Initialize RUC snow cover on ice from snow cover
        Sfcprop(nb)%sncovr_ice = Sfcprop(nb)%sncovr
        if (Model%rdlai) then
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xlaixy)
        end if
      elseif (Model%lsm == Model%lsm_noahmp) then
        !--- Extra Noah MP variables
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%snowxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tvxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tgxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%canicexy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%canliqxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%eahxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tahxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%cmxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%chxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%fwetxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%sneqvoxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%alboldxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%qsnowxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%wslakexy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%zwtxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%waxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%wtxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%lfmassxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%rtmassxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%stmassxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%woodxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%stblcpxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%fastcpxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xsaixy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xlaixy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%taussxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%smcwtdxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%deeprechxy)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%rechxy)
      endif
      if (Model%lkm > 0 .and. Model%iopt_lake==Model%iopt_lake_flake) then
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%T_snow)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%T_ice)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%h_ML)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%t_ML)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%t_mnw)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%h_talb)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%t_talb)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%t_bot1)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%t_bot2)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%c_t)
      endif
      if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp .or. (.not.warm_start)) then
        !--- 3D variables
        nt=0
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,1,Model%lsoil,sfc%var3,Sfcprop(nb)%stc)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,1,Model%lsoil,sfc%var3,Sfcprop(nb)%smc)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,1,Model%lsoil,sfc%var3,Sfcprop(nb)%slc)

        if (Model%lsm == Model%lsm_noahmp) then
          ! These use weird indexing which is lost during a Fortran subroutine call, so we use loops instead:
          nt=nt+1
          do lsoil = -2, 0
            do ix = 1, Atm_block%blksz(nb)
              Sfcprop(nb)%snicexy(ix,lsoil) = sfc%var3sn(ii1(ix),jj1(ix),lsoil,nt)
            enddo
          enddo

          nt=nt+1
          do lsoil = -2, 0
            do ix = 1, Atm_block%blksz(nb)
              Sfcprop(nb)%snliqxy(ix,lsoil) = sfc%var3sn(ii1(ix),jj1(ix),lsoil,nt)
            enddo
          enddo

          nt=nt+1
          do lsoil = -2, 0
            do ix = 1, Atm_block%blksz(nb)
              Sfcprop(nb)%tsnoxy(ix,lsoil)  = sfc%var3sn(ii1(ix),jj1(ix),lsoil,nt)
            enddo
          enddo

          nt=nt+1
          do lsoil = 1, 4
            do ix = 1, Atm_block%blksz(nb)
              Sfcprop(nb)%smoiseq(ix,lsoil)  = sfc%var3eq(ii1(ix),jj1(ix),lsoil,nt)
            enddo
          enddo

          nt=nt+1
          do lsoil = -2, 4
            do ix = 1, Atm_block%blksz(nb)
              Sfcprop(nb)%zsnsoxy(ix,lsoil)  = sfc%var3zn(ii1(ix),jj1(ix),lsoil,nt)
            enddo
          enddo
        endif

      else if (Model%lsm == Model%lsm_ruc) then
        !--- 3D variables
        nt=0
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,1,Model%lsoil_lsm,sfc%var3,Sfcprop(nb)%tslb)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,1,Model%lsoil_lsm,sfc%var3,Sfcprop(nb)%smois)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,1,Model%lsoil_lsm,sfc%var3,Sfcprop(nb)%sh2o)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,1,Model%lsoil_lsm,sfc%var3,Sfcprop(nb)%keepsmfr)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,1,Model%lsoil_lsm,sfc%var3,Sfcprop(nb)%flag_frsoil)
      endif

      do k = 1,Model%kice
        do ix = 1, Atm_block%blksz(nb)
          Sfcprop(nb)%tiice(ix,k) = sfc%var3ice(ii1(ix),jj1(ix),k)   !--- internal ice temp
        enddo
      enddo

      deallocate(ii1,jj1)

    end do block_loop
  end subroutine Sfc_io_copy_to_grid

  subroutine Sfc_io_copy_from_grid(sfc, Model, Atm_block, Sfcprop)
    !--- interface variable definitions
    implicit none

    class(Sfc_io_data_type)             :: sfc
    type(GFS_sfcprop_type),      intent(in) :: Sfcprop(:)
    type(block_control_type),    intent(in) :: Atm_block
    type(GFS_control_type),      intent(in) :: Model

    integer :: i, j, k, nb, ix, lsoil, num, nt
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer, allocatable :: ii1(:), jj1(:)
    real(kind_phys) :: ice

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)

    !$omp parallel do default(shared) private(i, j, nb, ix, nt, ii1, jj1, lsoil, k, ice)
    block_loop: do nb = 1, Atm_block%nblks
      allocate(ii1(Atm_block%blksz(nb)))
      allocate(jj1(Atm_block%blksz(nb)))
      ii1=Atm_block%index(nb)%ii - isc + 1
      jj1=Atm_block%index(nb)%jj - jsc + 1

      nt=0

      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%slmsk) !--- slmsk
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tsfco) !--- tsfc (tsea in sfc file)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%weasd) !--- weasd (sheleg in sfc file)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tg3)   !--- tg3
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%zorl)  !--- zorl
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%alvsf) !--- alvsf
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%alvwf) !--- alvwf
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%alnsf) !--- alnsf
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%alnwf) !--- alnwf
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%facsf) !--- facsf
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%facwf) !--- facwf
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%vfrac) !--- vfrac
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%canopy)!--- canopy
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%f10m)  !--- f10m
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%t2m)   !--- t2m
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%q2m)   !--- q2m

      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%vtype) !--- vtype
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%stype) !--- stype
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%uustar)!--- uustar
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%ffmm)  !--- ffmm
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%ffhh)  !--- ffhh
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%hice)  !--- hice
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%fice)  !--- fice
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tisfc) !--- tisfc
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tprcp) !--- tprcp
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%srflag)!--- srflag
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%snowd) !--- snowd (snwdph in the file)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%shdmin)!--- shdmin
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%shdmax)!--- shdmax
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%slope) !--- slope
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%snoalb)!--- snoalb
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%sncovr) !--- sncovr
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%snodl)  !--- snodl (snowd on land)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%weasdl) !--- weasdl (weasd on land)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tsfc)   !--- tsfc composite
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tsfcl)  !--- tsfcl (temp on land)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%zorlw)  !--- zorl (zorl on water)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%zorll)  !--- zorll (zorl on land)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%zorli)  !--- zorli (zorl on ice)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%albdirvis_lnd)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%albdirnir_lnd)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%albdifvis_lnd)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%albdifnir_lnd)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%emis_lnd)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%emis_ice)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%sncovr_ice)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%snodi)  !--- snodi (snowd on ice)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%weasdi) !--- weasdi (weasd on ice)
      if (Model%use_cice_alb .or. Model%lsm == Model%lsm_ruc) then
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%albdirvis_ice)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%albdifvis_ice)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%albdirnir_ice)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%albdifnir_ice)
        !        sfc%var2(i,j,53) = Sfcprop(nb)%sfalb_ice(ix)
      endif
      if (Model%cplwav) then
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%zorlwav) !--- zorlwav (zorl from wav)
      endif
      !--- NSSTM variables
      if (Model%nstf_name(1) > 0) then
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tref)   !--- nsstm tref
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%z_c)    !--- nsstm z_c
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%c_0)    !--- nsstm c_0
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%c_d)    !--- nsstm c_d
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%w_0)    !--- nsstm w_0
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%w_d)    !--- nsstm w_d
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xt)     !--- nsstm xt
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xs)     !--- nsstm xs
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xu)     !--- nsstm xu
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xv)     !--- nsstm xv
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xz)     !--- nsstm xz
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%zm)     !--- nsstm zm
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xtts)   !--- nsstm xtts
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xzts)   !--- nsstm xzts
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%d_conv) !--- nsstm d_conv
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%ifd)    !--- nsstm ifd
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%dt_cool)!--- nsstm dt_cool
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%qrain)  !--- nsstm qrain

        ! FIXME convert negative zero (-0.0) to zero (0.0)
        do j=1,ny
          do i=1,nx
            if(sfc%var2(i,j,nt) == 0.0) sfc%var2(i,j,nt) = 0.0
          end do
        end do
      endif

      if (Model%lsm == Model%lsm_ruc) then
        !--- Extra RUC variables
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%wetness)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%clw_surf_land)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%clw_surf_ice)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%qwv_surf_land)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%qwv_surf_ice)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tsnow_land)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tsnow_ice)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%snowfallac_land)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%snowfallac_ice)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%sfalb_lnd)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%sfalb_lnd_bck)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%sfalb_ice)
        if (Model%rdlai) then
          call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xlaixy)
        endif
      else if (Model%lsm == Model%lsm_noahmp) then
        !--- Extra Noah MP variables
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%snowxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tvxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tgxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%canicexy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%canliqxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%eahxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%tahxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%cmxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%chxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%fwetxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%sneqvoxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%alboldxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%qsnowxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%wslakexy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%zwtxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%waxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%wtxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%lfmassxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%rtmassxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%stmassxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%woodxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%stblcpxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%fastcpxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xsaixy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%xlaixy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%taussxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%smcwtdxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%deeprechxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop(nb)%rechxy)
      endif

      do k = 1,Model%kice
        do ix = 1, Atm_block%blksz(nb)
          ice=Sfcprop(nb)%tiice(ix,k)
          if(ice<one) then
            sfc%var3ice(ii1(ix),jj1(ix),k) = zero
          else
            sfc%var3ice(ii1(ix),jj1(ix),k) = ice
          endif
        enddo
      enddo

      if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp) then
        !--- 3D variables
        nt=0
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var3,Sfcprop(nb)%stc)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var3,Sfcprop(nb)%smc)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var3,Sfcprop(nb)%slc)
        ! 5 Noah MP 3D
        if (Model%lsm == Model%lsm_noahmp) then

          ! These arrays use bizarre indexing, which does not pass
          ! through function calls in Fortran, so we use loops here:

          nt=nt+1
          do lsoil = -2,0
            do ix = 1, Atm_block%blksz(nb)
              sfc%var3sn(ii1(ix),jj1(ix),lsoil,nt) = Sfcprop(nb)%snicexy(ix,lsoil)
            enddo
          enddo

          nt=nt+1
          do lsoil = -2,0
            do ix = 1, Atm_block%blksz(nb)
              sfc%var3sn(ii1(ix),jj1(ix),lsoil,nt) = Sfcprop(nb)%snliqxy(ix,lsoil)
            enddo
          enddo

          nt=nt+1
          do lsoil = -2,0
            do ix = 1, Atm_block%blksz(nb)
              sfc%var3sn(ii1(ix),jj1(ix),lsoil,nt) = Sfcprop(nb)%tsnoxy(ix,lsoil)
            enddo
          enddo

          nt=nt+1
          do lsoil = 1,Model%lsoil
            do ix = 1, Atm_block%blksz(nb)
              sfc%var3eq(ii1(ix),jj1(ix),lsoil,nt)  = Sfcprop(nb)%smoiseq(ix,lsoil)
            enddo
          enddo

          nt=nt+1
          do lsoil = -2,4
            do ix = 1, Atm_block%blksz(nb)
              sfc%var3zn(ii1(ix),jj1(ix),lsoil,nt)  = Sfcprop(nb)%zsnsoxy(ix,lsoil)
            enddo
          enddo

        endif  ! Noah MP
      else if (Model%lsm == Model%lsm_ruc) then
        !--- 3D variables
        nt=0
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var3,Sfcprop(nb)%tslb)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var3,Sfcprop(nb)%smois)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var3,Sfcprop(nb)%sh2o)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var3,Sfcprop(nb)%keepsmfr)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc%var3,Sfcprop(nb)%flag_frsoil)
      end if

      do k = 1,Model%kice
        do ix = 1, Atm_block%blksz(nb)
          Sfcprop(nb)%tiice(ix,k) = sfc%var3ice(ii1(ix),jj1(ix),k)   !--- internal ice temp
        enddo
      enddo
      
      deallocate(ii1,jj1)
    enddo block_loop

  end subroutine Sfc_io_copy_from_grid

  subroutine Sfc_io_apply_safeguards(sfc, Model, Atm_block, Sfcprop)
    !--- interface variable definitions
    implicit none

    class(Sfc_io_data_type)             :: sfc
    type(GFS_sfcprop_type),      intent(in) :: Sfcprop(:)
    type(block_control_type),    intent(in) :: Atm_block
    type(GFS_control_type),      intent(in) :: Model

    integer :: i, j, k, nb, ix, lsoil, num, nt
    integer :: isc, iec, jsc, jec, npz, nx, ny
    real(kind_phys) :: ice, tem, tem1

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)

    ! so far: At cold start everything is 9999.0, warm start snowxy has values
    !         but the 3D of snow fields are not available because not allocated yet.
    !         ix,nb loops may be consolidate with the Noah MP isnowxy init
    !         restore traditional vars first,we need some of them to init snow fields
    !         snow depth to actual snow layers; so we can allocate and register
    !         note zsnsoxy is from -2:4 - isnowxy is from 0:-2, but we need
    !         exact snow layers to pass 3D fields correctly, snow layers are
    !         different fro grid to grid, we have to init point by point/grid.
    !         It has to be done after the weasd is available
    !         sfc%var2(1,1,32) is the first; we need this to allocate snow related fields

    i = Atm_block%index(1)%ii(1) - isc + 1
    j = Atm_block%index(1)%jj(1) - jsc + 1

    if (sfc%var2(i,j,33) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing snodl')
      !$omp parallel do default(shared) private(nb, ix, tem)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          if (Sfcprop(nb)%landfrac(ix) > zero) then
            tem = one / (Sfcprop(nb)%fice(ix)*(one-Sfcprop(nb)%landfrac(ix))+Sfcprop(nb)%landfrac(ix))
            Sfcprop(nb)%snodl(ix)  = Sfcprop(nb)%snowd(ix) * tem
          else
            Sfcprop(nb)%snodl(ix)  = zero
          endif
        enddo
      enddo
    endif

    if (sfc%var2(i,j,34) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing weasdl')
      !$omp parallel do default(shared) private(nb, ix, tem)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          if (Sfcprop(nb)%landfrac(ix) > zero) then
            tem = one / (Sfcprop(nb)%fice(ix)*(one-Sfcprop(nb)%landfrac(ix))+Sfcprop(nb)%landfrac(ix))
            Sfcprop(nb)%weasdl(ix) = Sfcprop(nb)%weasd(ix) * tem
          else
            Sfcprop(nb)%weasdl(ix) = zero
          endif
        enddo
      enddo
    endif

    if (sfc%var2(i,j,36) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing tsfcl')
      !$omp parallel do default(shared) private(nb, ix)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          Sfcprop(nb)%tsfcl(ix) = Sfcprop(nb)%tsfco(ix) !--- compute tsfcl from existing variables
        enddo
      enddo
    endif

    if (sfc%var2(i,j,37) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing zorlw')
      !$omp parallel do default(shared) private(nb, ix)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          if (Sfcprop(nb)%landfrac(ix) < one .and. Sfcprop(nb)%fice(ix) < one) then
            Sfcprop(nb)%zorlw(ix) = min(Sfcprop(nb)%zorl(ix), 0.317)
          endif
        enddo
      enddo
    endif

    if (sfc%var2(i,j,38) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing zorll')
      !$omp parallel do default(shared) private(nb, ix)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          Sfcprop(nb)%zorll(ix) = Sfcprop(nb)%zorl(ix) !--- compute zorll from existing variables
        enddo
      enddo
    endif

    if (sfc%var2(i,j,39) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing zorli')
      !$omp parallel do default(shared) private(nb, ix)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          if (Sfcprop(nb)%fice(ix)*(one-Sfcprop(nb)%landfrac(ix)) > zero) then
            Sfcprop(nb)%zorli(ix) = one
          endif
        enddo
      enddo
    endif

    if (sfc%var2(i,j,45) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing emis_ice')
      !$omp parallel do default(shared) private(nb, ix)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          Sfcprop(nb)%emis_ice(ix) = 0.96
        enddo
      enddo
    endif

    if (sfc%var2(i,j,46) < -9990.0_kind_phys .and. Model%lsm /= Model%lsm_ruc) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing sncovr_ice')
      !$omp parallel do default(shared) private(nb, ix)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          !         Sfcprop(nb)%sncovr_ice(ix) = Sfcprop(nb)%sncovr(ix)
          Sfcprop(nb)%sncovr_ice(ix) = zero
        enddo
      enddo
    endif

    if (sfc%var2(i,j,47) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing snodi')
      !$omp parallel do default(shared) private(nb, ix, tem)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          if (Sfcprop(nb)%fice(ix) > zero) then
            tem = one / (Sfcprop(nb)%fice(ix)*(one-Sfcprop(nb)%landfrac(ix))+Sfcprop(nb)%landfrac(ix))
            Sfcprop(nb)%snodi(ix)  = min(Sfcprop(nb)%snowd(ix) * tem, 3.0)
          else
            Sfcprop(nb)%snodi(ix)  = zero
          endif
        enddo
      enddo
    endif

    if (sfc%var2(i,j,48) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing weasdi')
      !$omp parallel do default(shared) private(nb, ix, tem)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          if (Sfcprop(nb)%fice(ix) > zero) then
            tem = one / (Sfcprop(nb)%fice(ix)*(one-Sfcprop(nb)%landfrac(ix))+Sfcprop(nb)%landfrac(ix))
            Sfcprop(nb)%weasdi(ix)  = Sfcprop(nb)%weasd(ix)*tem
          else
            Sfcprop(nb)%weasdi(ix)  = zero
          endif
        enddo
      enddo
    endif

    if (Model%use_cice_alb) then
      if (sfc%var2(i,j,49) < -9990.0_kind_phys) then
        !$omp parallel do default(shared) private(nb, ix)
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            if (Sfcprop(nb)%oceanfrac(ix) > zero .and. &
                 Sfcprop(nb)%fice(ix) >= Model%min_seaice) then
              Sfcprop(nb)%albdirvis_ice(ix) = 0.6_kind_phys
              Sfcprop(nb)%albdifvis_ice(ix) = 0.6_kind_phys
              Sfcprop(nb)%albdirnir_ice(ix) = 0.6_kind_phys
              Sfcprop(nb)%albdifnir_ice(ix) = 0.6_kind_phys
            endif
          enddo
        enddo
      endif

    endif

    ! Fill in composite tsfc for coldstart runs - must happen after tsfcl is computed
    compute_tsfc_for_colstart: if (sfc%var2(i,j,35) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing composite tsfc')
      if(Model%frac_grid) then ! 3-way composite
        !$omp parallel do default(shared) private(nb, ix, tem, tem1)
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            Sfcprop(nb)%tsfco(ix) = max(con_tice, Sfcprop(nb)%tsfco(ix)) ! this may break restart reproducibility
            tem1 = one - Sfcprop(nb)%landfrac(ix)
            tem  = tem1 * Sfcprop(nb)%fice(ix) ! tem = ice fraction wrt whole cell
            Sfcprop(nb)%tsfc(ix) = Sfcprop(nb)%tsfcl(ix) * Sfcprop(nb)%landfrac(ix) &
                 + Sfcprop(nb)%tisfc(ix) * tem                      &
                 + Sfcprop(nb)%tsfco(ix) * (tem1-tem)
          enddo
        enddo
      else
        !$omp parallel do default(shared) private(nb, ix, tem)
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            if (Sfcprop(nb)%slmsk(ix) == 1) then
              Sfcprop(nb)%tsfc(ix) = Sfcprop(nb)%tsfcl(ix)
            else
              tem = one - Sfcprop(nb)%fice(ix)
              Sfcprop(nb)%tsfc(ix) = Sfcprop(nb)%tisfc(ix) * Sfcprop(nb)%fice(ix) &
                   + Sfcprop(nb)%tsfco(ix) * tem
            endif
          enddo
        enddo
      endif
    endif compute_tsfc_for_colstart

    if (sfc%var2(i,j,sfc%nvar2m) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing zorlwav')
      !$omp parallel do default(shared) private(nb, ix)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          Sfcprop(nb)%zorlwav(ix) = Sfcprop(nb)%zorl(ix) !--- compute zorlwav from existing variables
        enddo
      enddo
    endif

    if (nint(sfc%var3ice(1,1,1)) == -9999) then    !--- initialize internal ice temp from layer 1 and 2 soil temp
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing tiice')
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          Sfcprop(nb)%tiice(ix,1) = max(timin, min(con_tice, Sfcprop(nb)%stc(ix,1)))
          Sfcprop(nb)%tiice(ix,2) = max(timin, min(con_tice, Sfcprop(nb)%stc(ix,2)))
        enddo
      enddo
    endif

  end subroutine Sfc_io_apply_safeguards

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
