module FV3GFS_clm_lake_io
  use GFS_typedefs,       only: GFS_sfcprop_type, GFS_control_type, &
                                GFS_data_type, kind_phys
  use GFS_restart,        only: GFS_restart_type
  use GFS_diagnostics,    only: GFS_externaldiag_type
  use block_control_mod,  only: block_control_type
  use fms2_io_mod,        only: FmsNetcdfDomainFile_t, unlimited,      &
                                open_file, close_file,                 &
                                register_axis, register_restart_field, &
                                register_variable_attribute, register_field, &
                                read_restart, write_restart, write_data,     &
                                get_global_io_domain_indices, variable_exists
  use FV3GFS_common_io

  implicit none

  type clm_lake_data_type
    ! The clm_lake_data_type derived type is a class that stores
    ! temporary arrays used to read or write CLM Lake model restart
    ! and axis variables. It can safely be declared and unused, but
    ! you should only call these routines if the CLM Lake Model was
    ! (or will be) used by this execution of the FV3. It is the
    ! responsibility of the caller to ensure the necessary data is in
    ! Sfc_restart, Sfcprop, and Model.

    ! All 2D variables needed for a restart
    real(kind_phys), pointer, private, dimension(:,:) :: &
         T_snow=>null(), T_ice=>null(), &
         lake_snl2d=>null(), lake_h2osno2d=>null(), lake_tsfc=>null(), clm_lakedepth=>null(), &
         lake_savedtke12d=>null(), lake_sndpth2d=>null(), clm_lake_initialized=>null()

    ! All 3D variables needed for a restart
    real(kind_phys), pointer, private, dimension(:,:,:) :: &
         lake_z3d=>null(), lake_dz3d=>null(), lake_soil_watsat3d=>null(), &
         lake_csol3d=>null(), lake_soil_tkmg3d=>null(), lake_soil_tkdry3d=>null(), &
         lake_soil_tksatu3d=>null(), lake_snow_z3d=>null(), lake_snow_dz3d=>null(), &
         lake_snow_zi3d=>null(), lake_h2osoi_vol3d=>null(), lake_h2osoi_liq3d=>null(), &
         lake_h2osoi_ice3d=>null(), lake_t_soisno3d=>null(), lake_t_lake3d=>null(), &
         lake_icefrac3d=>null(),  lake_clay3d=>null(), lake_sand3d=>null()

  contains

    ! register_axes calls registers_axis on Sfc_restart for all required axes
    procedure, public :: register_axes => clm_lake_register_axes

    ! allocate_data allocates all of the pointers in this object
    procedure, public :: allocate_data => clm_lake_allocate_data

    ! register_fields calls register_field on Sfc_restart for all CLM Lake model restart variables
    procedure, public :: register_fields => clm_lake_register_fields

    ! deallocate_data deallocates all pointers, allowing this object to be used repeatedly.
    ! It is safe to call deallocate_data if no data has been allocated.
    procedure, public :: deallocate_data => clm_lake_deallocate_data

    ! write_axes writes variables to Sfc_restart, with the name of
    ! each axis, containing the appropriate information
    procedure, public :: write_axes => clm_lake_write_axes

    ! copy_from_grid copies from Sfcprop to internal pointers (declared above)
    procedure, public :: copy_from_grid => clm_lake_copy_from_grid

    ! copy_from_grid copies from internal pointers (declared above) to Sfcprop
    procedure, public :: copy_to_grid => clm_lake_copy_to_grid

    ! send field bundles in restart quilt server
    procedure, public :: bundle_fields => clm_lake_bundle_fields

    ! A fortran 2003 compliant compiler will call clm_lake_final
    ! automatically when an object of this type goes out of
    ! scope. This will deallocate any arrays via a call to
    ! deallocate_data. It is safe to call this routine if no data has
    ! been allocated.
    final :: clm_lake_final
  end type clm_lake_data_type

   CONTAINS
  subroutine clm_lake_allocate_data(data,Model)
    ! Deallocate all data, and reallocate to the size specified in Model
    implicit none
    class(clm_lake_data_type) :: data
    type(GFS_control_type),   intent(in) :: Model

    integer :: nx, ny
    call data%deallocate_data

    nx=Model%nx
    ny=Model%ny

    allocate(data%T_snow(nx,ny))
    allocate(data%T_ice(nx,ny))
    allocate(data%lake_snl2d(nx,ny))
    allocate(data%lake_h2osno2d(nx,ny))
    allocate(data%lake_tsfc(nx,ny))
    allocate(data%lake_savedtke12d(nx,ny))
    allocate(data%lake_sndpth2d(nx,ny))
    allocate(data%clm_lakedepth(nx,ny))
    allocate(data%clm_lake_initialized(nx,ny))
    
    allocate(data%lake_z3d(nx,ny,Model%nlevlake_clm_lake))
    allocate(data%lake_dz3d(nx,ny,Model%nlevlake_clm_lake))
    allocate(data%lake_soil_watsat3d(nx,ny,Model%nlevlake_clm_lake))
    allocate(data%lake_csol3d(nx,ny,Model%nlevlake_clm_lake))
    allocate(data%lake_soil_tkmg3d(nx,ny,Model%nlevlake_clm_lake))
    allocate(data%lake_soil_tkdry3d(nx,ny,Model%nlevlake_clm_lake))
    allocate(data%lake_soil_tksatu3d(nx,ny,Model%nlevlake_clm_lake))
    allocate(data%lake_snow_z3d(nx,ny,Model%nlevsnowsoil1_clm_lake))
    allocate(data%lake_snow_dz3d(nx,ny,Model%nlevsnowsoil1_clm_lake))
    allocate(data%lake_snow_zi3d(nx,ny,Model%nlevsnowsoil_clm_lake))
    allocate(data%lake_h2osoi_vol3d(nx,ny,Model%nlevsnowsoil1_clm_lake))
    allocate(data%lake_h2osoi_liq3d(nx,ny,Model%nlevsnowsoil1_clm_lake))
    allocate(data%lake_h2osoi_ice3d(nx,ny,Model%nlevsnowsoil1_clm_lake))
    allocate(data%lake_t_soisno3d(nx,ny,Model%nlevsnowsoil1_clm_lake))
    allocate(data%lake_t_lake3d(nx,ny,Model%nlevlake_clm_lake))
    allocate(data%lake_icefrac3d(nx,ny,Model%nlevlake_clm_lake))
    allocate(data%lake_clay3d(nx,ny,Model%nlevsoil_clm_lake))
    allocate(data%lake_sand3d(nx,ny,Model%nlevsoil_clm_lake))
  end subroutine clm_lake_allocate_data

  subroutine clm_lake_register_axes(data,Model,Sfc_restart)
    ! Register all five axes needed by CLM Lake restart data
    implicit none
    class(clm_lake_data_type) :: data
    type(GFS_control_type),      intent(in) :: Model
    type(FmsNetcdfDomainFile_t) :: Sfc_restart
    call register_axis(Sfc_restart, 'levlake_clm_lake', dimension_length=Model%nlevlake_clm_lake)

    call register_axis(Sfc_restart, 'levsoil_clm_lake', dimension_length=Model%nlevsoil_clm_lake)

    call register_axis(Sfc_restart, 'levsnow_clm_lake', dimension_length=Model%nlevsnow_clm_lake)

    call register_axis(Sfc_restart, 'levsnowsoil_clm_lake', dimension_length=Model%nlevsnowsoil_clm_lake)

    call register_axis(Sfc_restart, 'levsnowsoil1_clm_lake', dimension_length=Model%nlevsnowsoil1_clm_lake)
  end subroutine clm_lake_register_axes

  subroutine clm_lake_write_axes(data, Model, Sfc_restart)
    ! Create variables with the name name as each clm_lake axis, and
    ! fill the variable with the appropriate indices
    implicit none
    class(clm_lake_data_type) :: data
    type(GFS_control_type),      intent(in) :: Model
    type(FmsNetcdfDomainFile_t) :: Sfc_restart
    real(kind_phys) :: levlake_clm_lake(Model%nlevlake_clm_lake)
    real(kind_phys) :: levsoil_clm_lake(Model%nlevsoil_clm_lake)
    real(kind_phys) :: levsnow_clm_lake(Model%nlevsnow_clm_lake)
    real(kind_phys) :: levsnowsoil_clm_lake(Model%nlevsnowsoil_clm_lake)
    real(kind_phys) :: levsnowsoil1_clm_lake(Model%nlevsnowsoil1_clm_lake)
    integer :: i
    call register_field(Sfc_restart, 'levlake_clm_lake', 'double', (/'levlake_clm_lake'/))
    call register_variable_attribute(Sfc_restart, 'levlake_clm_lake', 'cartesian_axis' ,'Z', str_len=1)

    call register_field(Sfc_restart, 'levsoil_clm_lake', 'double', (/'levsoil_clm_lake'/))
    call register_variable_attribute(Sfc_restart, 'levsoil_clm_lake', 'cartesian_axis' ,'Z', str_len=1)

    call register_field(Sfc_restart, 'levsnow_clm_lake', 'double', (/'levsnow_clm_lake'/))
    call register_variable_attribute(Sfc_restart, 'levsnow_clm_lake', 'cartesian_axis' ,'Z', str_len=1)

    call register_field(Sfc_restart, 'levsnowsoil_clm_lake', 'double', (/'levsnowsoil_clm_lake'/))
    call register_variable_attribute(Sfc_restart, 'levsnowsoil_clm_lake', 'cartesian_axis' ,'Z', str_len=1)

    call register_field(Sfc_restart, 'levsnowsoil1_clm_lake', 'double', (/'levsnowsoil1_clm_lake'/))
    call register_variable_attribute(Sfc_restart, 'levsnowsoil1_clm_lake', 'cartesian_axis' ,'Z', str_len=1)

    do i=1,Model%nlevlake_clm_lake
      levlake_clm_lake(i) = i
    enddo
    do i=1,Model%nlevsoil_clm_lake
      levsoil_clm_lake(i) = i
    enddo
    do i=1,Model%nlevsnow_clm_lake
      levsnow_clm_lake(i) = i
    enddo
    do i=-Model%nlevsnow_clm_lake,Model%nlevsoil_clm_lake
      levsnowsoil_clm_lake(i+Model%nlevsnow_clm_lake+1) = i
    enddo
    do i=-Model%nlevsnow_clm_lake+1,Model%nlevsoil_clm_lake
      levsnowsoil1_clm_lake(i+Model%nlevsnow_clm_lake) = i
    enddo

    call write_data(Sfc_restart, 'levlake_clm_lake', levlake_clm_lake)
    call write_data(Sfc_restart, 'levsoil_clm_lake', levsoil_clm_lake)
    call write_data(Sfc_restart, 'levsnow_clm_lake', levsnow_clm_lake)
    call write_data(Sfc_restart, 'levsnowsoil_clm_lake', levsnowsoil_clm_lake)
    call write_data(Sfc_restart, 'levsnowsoil1_clm_lake', levsnowsoil1_clm_lake)
  end subroutine clm_lake_write_axes

  subroutine clm_lake_copy_from_grid(data, Model, Atm_block, Sfcprop)
    ! Copies from Sfcprop variables to the corresponding data temporary variables.
    ! Terrible things will happen if you don't call data%allocate_data first.
    implicit none
    class(clm_lake_data_type) :: data
    type(GFS_sfcprop_type),   intent(in) :: Sfcprop(:)
    type(GFS_control_type),   intent(in) :: Model
    type(block_control_type), intent(in) :: Atm_block

    integer :: nb, ix, isc, jsc, i, j
    isc = Model%isc
    jsc = Model%jsc

    ! Copy data to temporary arrays

!$omp parallel do default(shared) private(i, j, nb, ix)
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1

        data%T_snow(i,j) = Sfcprop(nb)%T_snow(ix)
        data%T_ice(i,j) = Sfcprop(nb)%T_ice(ix)
        data%lake_snl2d(i,j) = Sfcprop(nb)%lake_snl2d(ix)
        data%lake_h2osno2d(i,j) = Sfcprop(nb)%lake_h2osno2d(ix)
        data%lake_tsfc(i,j) = Sfcprop(nb)%lake_tsfc(ix)
        data%lake_savedtke12d(i,j) = Sfcprop(nb)%lake_savedtke12d(ix)
        data%lake_sndpth2d(i,j) = Sfcprop(nb)%lake_sndpth2d(ix)
        data%clm_lakedepth(i,j) = Sfcprop(nb)%clm_lakedepth(ix)
        data%clm_lake_initialized(i,j) = Sfcprop(nb)%clm_lake_initialized(ix)

        data%lake_z3d(i,j,:) = Sfcprop(nb)%lake_z3d(ix,:)
        data%lake_dz3d(i,j,:) = Sfcprop(nb)%lake_dz3d(ix,:)
        data%lake_soil_watsat3d(i,j,:) = Sfcprop(nb)%lake_soil_watsat3d(ix,:)
        data%lake_csol3d(i,j,:) = Sfcprop(nb)%lake_csol3d(ix,:)
        data%lake_soil_tkmg3d(i,j,:) = Sfcprop(nb)%lake_soil_tkmg3d(ix,:)
        data%lake_soil_tkdry3d(i,j,:) = Sfcprop(nb)%lake_soil_tkdry3d(ix,:)
        data%lake_soil_tksatu3d(i,j,:) = Sfcprop(nb)%lake_soil_tksatu3d(ix,:)
        data%lake_snow_z3d(i,j,:) = Sfcprop(nb)%lake_snow_z3d(ix,:)
        data%lake_snow_dz3d(i,j,:) = Sfcprop(nb)%lake_snow_dz3d(ix,:)
        data%lake_snow_zi3d(i,j,:) = Sfcprop(nb)%lake_snow_zi3d(ix,:)
        data%lake_h2osoi_vol3d(i,j,:) = Sfcprop(nb)%lake_h2osoi_vol3d(ix,:)
        data%lake_h2osoi_liq3d(i,j,:) = Sfcprop(nb)%lake_h2osoi_liq3d(ix,:)
        data%lake_h2osoi_ice3d(i,j,:) = Sfcprop(nb)%lake_h2osoi_ice3d(ix,:)
        data%lake_t_soisno3d(i,j,:) = Sfcprop(nb)%lake_t_soisno3d(ix,:)
        data%lake_t_lake3d(i,j,:) = Sfcprop(nb)%lake_t_lake3d(ix,:)
        data%lake_icefrac3d(i,j,:) = Sfcprop(nb)%lake_icefrac3d(ix,:)
        data%lake_clay3d(i,j,:) = Sfcprop(nb)%lake_clay3d(ix,:)
        data%lake_sand3d(i,j,:) = Sfcprop(nb)%lake_sand3d(ix,:)
      enddo
    enddo
  end subroutine clm_lake_copy_from_grid

  subroutine clm_lake_copy_to_grid(data, Model, Atm_block, Sfcprop)
    ! Copies from data temporary variables to the corresponding Sfcprop variables.
    ! Terrible things will happen if you don't call data%allocate_data first.
    implicit none
    class(clm_lake_data_type) :: data
    type(GFS_sfcprop_type),   intent(in) :: Sfcprop(:)
    type(GFS_control_type),   intent(in) :: Model
    type(block_control_type), intent(in) :: Atm_block

    integer :: nb, ix, isc, jsc, i, j
    isc = Model%isc
    jsc = Model%jsc

    ! Copy data to temporary arrays

!$omp parallel do default(shared) private(i, j, nb, ix)
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1

        Sfcprop(nb)%T_snow(ix) = data%T_snow(i,j)
        Sfcprop(nb)%T_ice(ix) = data%T_ice(i,j)
        Sfcprop(nb)%lake_snl2d(ix) = data%lake_snl2d(i,j)
        Sfcprop(nb)%lake_h2osno2d(ix) = data%lake_h2osno2d(i,j)
        Sfcprop(nb)%lake_tsfc(ix) = data%lake_tsfc(i,j)
        Sfcprop(nb)%lake_savedtke12d(ix) = data%lake_savedtke12d(i,j)
        Sfcprop(nb)%lake_sndpth2d(ix) = data%lake_sndpth2d(i,j)
        Sfcprop(nb)%clm_lakedepth(ix) = data%clm_lakedepth(i,j)
        Sfcprop(nb)%clm_lake_initialized(ix) = data%clm_lake_initialized(i,j)

        Sfcprop(nb)%lake_z3d(ix,:) = data%lake_z3d(i,j,:)
        Sfcprop(nb)%lake_dz3d(ix,:) = data%lake_dz3d(i,j,:)
        Sfcprop(nb)%lake_soil_watsat3d(ix,:) = data%lake_soil_watsat3d(i,j,:)
        Sfcprop(nb)%lake_csol3d(ix,:) = data%lake_csol3d(i,j,:)
        Sfcprop(nb)%lake_soil_tkmg3d(ix,:) = data%lake_soil_tkmg3d(i,j,:)
        Sfcprop(nb)%lake_soil_tkdry3d(ix,:) = data%lake_soil_tkdry3d(i,j,:)
        Sfcprop(nb)%lake_soil_tksatu3d(ix,:) = data%lake_soil_tksatu3d(i,j,:)
        Sfcprop(nb)%lake_snow_z3d(ix,:) = data%lake_snow_z3d(i,j,:)
        Sfcprop(nb)%lake_snow_dz3d(ix,:) = data%lake_snow_dz3d(i,j,:)
        Sfcprop(nb)%lake_snow_zi3d(ix,:) = data%lake_snow_zi3d(i,j,:)
        Sfcprop(nb)%lake_h2osoi_vol3d(ix,:) = data%lake_h2osoi_vol3d(i,j,:)
        Sfcprop(nb)%lake_h2osoi_liq3d(ix,:) = data%lake_h2osoi_liq3d(i,j,:)
        Sfcprop(nb)%lake_h2osoi_ice3d(ix,:) = data%lake_h2osoi_ice3d(i,j,:)
        Sfcprop(nb)%lake_t_soisno3d(ix,:) = data%lake_t_soisno3d(i,j,:)
        Sfcprop(nb)%lake_t_lake3d(ix,:) = data%lake_t_lake3d(i,j,:)
        Sfcprop(nb)%lake_icefrac3d(ix,:) = data%lake_icefrac3d(i,j,:)
        Sfcprop(nb)%lake_clay3d(ix,:) = data%lake_clay3d(i,j,:)
        Sfcprop(nb)%lake_sand3d(ix,:) = data%lake_sand3d(i,j,:)
      enddo
    enddo
  end subroutine clm_lake_copy_to_grid

  subroutine clm_lake_register_fields(data, Sfc_restart)
    ! Registers all restart fields needed by the CLM Lake Model.
    ! Terrible things will happen if you don't call data%allocate_data
    ! and data%register_axes first.
    implicit none
    class(clm_lake_data_type) :: data
    type(FmsNetcdfDomainFile_t) :: Sfc_restart

    ! Register 2D fields
    call register_restart_field(Sfc_restart, 'T_snow', data%T_snow, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)
    call register_restart_field(Sfc_restart, 'T_ice', data%T_ice, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)
    call register_restart_field(Sfc_restart, 'lake_snl2d', data%lake_snl2d, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)
    call register_restart_field(Sfc_restart, 'lake_h2osno2d', data%lake_h2osno2d, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)
    call register_restart_field(Sfc_restart, 'lake_tsfc', data%lake_tsfc, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)
    call register_restart_field(Sfc_restart, 'lake_savedtke12d', data%lake_savedtke12d, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)
    call register_restart_field(Sfc_restart, 'lake_sndpth2d', data%lake_sndpth2d, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)
    call register_restart_field(Sfc_restart, 'clm_lakedepth', data%clm_lakedepth, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)
    call register_restart_field(Sfc_restart, 'clm_lake_initialized', data%clm_lake_initialized, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)

    ! Register 3D fields
    call register_restart_field(Sfc_restart, 'lake_z3d', data%lake_z3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
                      'levlake_clm_lake     ', 'Time                 '/), is_optional=.true.)
    call register_restart_field(Sfc_restart, 'lake_dz3d', data%lake_dz3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
                      'levlake_clm_lake     ', 'Time                 '/), is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_soil_watsat3d', data%lake_soil_watsat3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
                      'levlake_clm_lake     ', 'Time                 '/), is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_csol3d', data%lake_csol3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
                      'levlake_clm_lake     ', 'Time                 '/), is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_soil_tkmg3d', data%lake_soil_tkmg3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
                      'levlake_clm_lake     ', 'Time                 '/), is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_soil_tkdry3d', data%lake_soil_tkdry3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
                      'levlake_clm_lake     ', 'Time                 '/), is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_soil_tksatu3d', data%lake_soil_tksatu3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
                      'levlake_clm_lake     ', 'Time                 '/), is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_snow_z3d', data%lake_snow_z3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
                      'levsnowsoil1_clm_lake', 'Time                 '/), is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_snow_dz3d', data%lake_snow_dz3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
                      'levsnowsoil1_clm_lake', 'Time                 '/), is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_snow_zi3d', data%lake_snow_zi3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
                      'levsnowsoil_clm_lake ', 'Time                 '/), is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_h2osoi_vol3d', data%lake_h2osoi_vol3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
                      'levsnowsoil1_clm_lake', 'Time                 '/), is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_h2osoi_liq3d', data%lake_h2osoi_liq3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
                      'levsnowsoil1_clm_lake', 'Time                 '/), is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_h2osoi_ice3d', data%lake_h2osoi_ice3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
                      'levsnowsoil1_clm_lake', 'Time                 '/), is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_t_soisno3d', data%lake_t_soisno3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
                      'levsnowsoil1_clm_lake', 'Time                 '/), is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_t_lake3d', data%lake_t_lake3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
                      'levlake_clm_lake     ', 'Time                 '/), is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_icefrac3d', data%lake_icefrac3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
                      'levlake_clm_lake     ', 'Time                 '/), is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_clay3d', data%lake_clay3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
                      'levsoil_clm_lake     ', 'Time                 '/), is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_sand3d', data%lake_sand3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
                      'levsoil_clm_lake     ', 'Time                 '/), is_optional=.true.)
  end subroutine clm_lake_register_fields


  subroutine clm_lake_bundle_fields(clm_lake, bundle, grid, Model, outputfile)
    use esmf
    use GFS_typedefs, only: GFS_control_type
    implicit none
    class(Clm_lake_data_type)             :: clm_lake
    type(ESMF_FieldBundle),intent(inout)        :: bundle
    type(ESMF_Grid),intent(inout)               :: grid
    type(GFS_control_type),          intent(in) :: Model
    character(*), intent(in)                    :: outputfile

    real(kind_phys),dimension(:,:),pointer :: temp_r2d
    real(kind_phys),dimension(:,:,:),pointer :: temp_r3d
    integer :: num

    ! Register 2D fields
    call create_2d_field_and_add_to_bundle(clm_lake%T_snow, "T_snow", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(clm_lake%T_ice, 'T_ice', trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(clm_lake%lake_snl2d, "lake_snl2d", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(clm_lake%lake_h2osno2d, "lake_h2osno2d", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(clm_lake%lake_tsfc, "lake_tsfc", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(clm_lake%lake_savedtke12d, "lake_savedtke12d", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(clm_lake%lake_sndpth2d, "lake_sndpth2d", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(clm_lake%clm_lakedepth, "clm_lakedepth", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(clm_lake%clm_lake_initialized, "clm_lake_initialized", trim(outputfile), grid, bundle)

    ! Register 3D fields
    call create_3d_field_and_add_to_bundle(clm_lake%lake_z3d, 'lake_z3d', 'levlake_clm_lake', Model%nlevlake_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_dz3d, 'lake_dz3d', 'levlake_clm_lake', Model%nlevlake_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_soil_watsat3d, 'lake_soil_watsat3d', 'levlake_clm_lake', Model%nlevlake_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_csol3d, 'lake_csol3d', 'levlake_clm_lake', Model%nlevlake_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_soil_tkmg3d, 'lake_soil_tkmg3d', 'levlake_clm_lake', Model%nlevlake_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_soil_tkdry3d, 'lake_soil_tkdry3d', 'levlake_clm_lake', Model%nlevlake_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_soil_tksatu3d, 'lake_soil_tksatu3d', 'levlake_clm_lake', Model%nlevlake_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_snow_z3d, 'lake_snow_z3d', 'levsnowsoil1_clm_lake', Model%nlevsnowsoil1_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_snow_dz3d, 'lake_snow_dz3d', 'levsnowsoil1_clm_lake', Model%nlevsnowsoil1_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_snow_zi3d, 'lake_snow_zi3d', 'levsnowsoil_clm_lake', Model%nlevsnowsoil_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_h2osoi_vol3d, 'lake_h2osoi_vol3d', 'levsnowsoil1_clm_lake', Model%nlevsnowsoil1_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_h2osoi_liq3d, 'lake_h2osoi_liq3d', 'levsnowsoil1_clm_lake', Model%nlevsnowsoil1_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_h2osoi_ice3d, 'lake_h2osoi_ice3d', 'levsnowsoil1_clm_lake', Model%nlevsnowsoil1_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_t_soisno3d, 'lake_t_soisno3d', 'levsnowsoil1_clm_lake', Model%nlevsnowsoil1_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_t_lake3d, 'lake_t_lake3d', 'levlake_clm_lake', Model%nlevlake_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_icefrac3d, 'lake_icefrac3d', 'levlake_clm_lake', Model%nlevlake_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_clay3d, 'lake_clay3d', 'levsoil_clm_lake', Model%nlevsoil_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_sand3d, 'lake_sand3d', 'levsoil_clm_lake', Model%nlevsoil_clm_lake, trim(outputfile), grid, bundle)

  end subroutine Clm_lake_bundle_fields

  subroutine clm_lake_final(data)
    ! Final routine for clm_lake_data_type, called automatically when
    ! an object of that type goes out of scope.  This is simply a
    ! wrapper around data%deallocate_data().
    implicit none
    type(clm_lake_data_type) :: data
    call clm_lake_deallocate_data(data)
  end subroutine clm_lake_final

  subroutine clm_lake_deallocate_data(data)
    ! Deallocates all data used, and nullifies the pointers. The data
    ! object can safely be used again after this call. This is also
    ! the implementation of the clm_lake_data_type final routine.
    implicit none
    class(clm_lake_data_type) :: data

    ! Deallocate and nullify any associated pointers

    ! This #define reduces code length by a lot
#define IF_ASSOC_DEALLOC_NULL(var) \
    if(associated(data%var)) then ; \
      deallocate(data%var) ; \
      nullify(data%var) ; \
    endif
    
    IF_ASSOC_DEALLOC_NULL(T_snow)
    IF_ASSOC_DEALLOC_NULL(T_ice)
    IF_ASSOC_DEALLOC_NULL(lake_snl2d)
    IF_ASSOC_DEALLOC_NULL(lake_h2osno2d)
    IF_ASSOC_DEALLOC_NULL(lake_tsfc)
    IF_ASSOC_DEALLOC_NULL(lake_savedtke12d)
    IF_ASSOC_DEALLOC_NULL(lake_sndpth2d)
    IF_ASSOC_DEALLOC_NULL(clm_lakedepth)
    IF_ASSOC_DEALLOC_NULL(clm_lake_initialized)

    IF_ASSOC_DEALLOC_NULL(lake_z3d)
    IF_ASSOC_DEALLOC_NULL(lake_dz3d)
    IF_ASSOC_DEALLOC_NULL(lake_soil_watsat3d)
    IF_ASSOC_DEALLOC_NULL(lake_csol3d)
    IF_ASSOC_DEALLOC_NULL(lake_soil_tkmg3d)
    IF_ASSOC_DEALLOC_NULL(lake_soil_tkdry3d)
    IF_ASSOC_DEALLOC_NULL(lake_soil_tksatu3d)
    IF_ASSOC_DEALLOC_NULL(lake_snow_z3d)
    IF_ASSOC_DEALLOC_NULL(lake_snow_dz3d)
    IF_ASSOC_DEALLOC_NULL(lake_snow_zi3d)
    IF_ASSOC_DEALLOC_NULL(lake_h2osoi_vol3d)
    IF_ASSOC_DEALLOC_NULL(lake_h2osoi_liq3d)
    IF_ASSOC_DEALLOC_NULL(lake_h2osoi_ice3d)
    IF_ASSOC_DEALLOC_NULL(lake_t_soisno3d)
    IF_ASSOC_DEALLOC_NULL(lake_t_lake3d)
    IF_ASSOC_DEALLOC_NULL(lake_icefrac3d)
    IF_ASSOC_DEALLOC_NULL(lake_clay3d)
    IF_ASSOC_DEALLOC_NULL(lake_sand3d)
    
#undef IF_ASSOC_DEALLOC_NULL
  end subroutine clm_lake_deallocate_data

end module FV3GFS_clm_lake_io
