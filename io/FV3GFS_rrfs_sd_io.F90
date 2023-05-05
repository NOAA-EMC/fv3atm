module rrfs_sd_io
  use block_control_mod,  only: block_control_type
  use fms2_io_mod,        only: FmsNetcdfDomainFile_t, unlimited,      &
                                open_file, close_file,                 &
                                register_axis, register_restart_field, &
                                register_variable_attribute, register_field, &
                                read_restart, write_restart, write_data,     &
                                get_global_io_domain_indices, variable_exists
  use GFS_typedefs,       only: GFS_sfcprop_type, GFS_control_type, &
                                GFS_data_type, kind_phys
  use GFS_restart,        only: GFS_restart_type
  use GFS_diagnostics,    only: GFS_externaldiag_type

  implicit none

  private

  public :: rrfs_sd_state_type, rrfs_sd_state_register_axis, rrfs_sd_state_write_axis, &
       rrfs_sd_state_fill_data, rrfs_sd_state_register_fields, rrfs_sd_state_deallocate_data, &
       rrfs_sd_state_copy_to_temporaries, rrfs_sd_state_copy_from_temporaries, &
       rrfs_sd_state_final

  type rrfs_sd_state_type
    ! The rrfs_sd_state_type stores temporary arrays used to read or
    ! write RRFS-SD restart and axis variables.

    real(kind_phys), pointer, private, dimension(:,:) :: & ! i,j variables
         emdust=>null(), emseas=>null(), emanoc=>null(), fhist=>null(), coef_bb_dc=>null()

    real(kind_phys), pointer, private, dimension(:,:,:) :: &
         fire_in=>null() ! i, j, fire_aux_data_levels

  contains
    procedure, public :: register_axis => rrfs_sd_state_register_axis ! register fire_aux_data_levels axis
    procedure, public :: write_axis => rrfs_sd_state_write_axis ! write fire_aux_data_levels variable
    procedure, public :: allocate_data => rrfs_sd_state_allocate_data ! allocate all pointers
    procedure, public :: fill_data => rrfs_sd_state_fill_data ! fill data with default values
    procedure, public :: register_fields => rrfs_sd_state_register_fields ! register rrfs_sd fields
    procedure, public :: deallocate_data => rrfs_sd_state_deallocate_data ! deallocate pointers
    procedure, public :: copy_to_temporaries => rrfs_sd_state_copy_to_temporaries ! Copy Sfcprop to arrays
    procedure, public :: copy_from_temporaries => rrfs_sd_state_copy_from_temporaries ! Copy arrays to Sfcprop
    final :: rrfs_sd_state_final ! Destructor; calls deallocate_data
  end type rrfs_sd_state_type

contains

  subroutine rrfs_sd_state_register_axis(data,Model,Sfc_restart)
    implicit none
    class(rrfs_sd_state_type) :: data
    type(FmsNetcdfDomainFile_t) :: Sfc_restart
    type(GFS_control_type),      intent(in) :: Model
    call register_axis(Sfc_restart, 'fire_aux_data_levels', &
         dimension_length=Model%fire_aux_data_levels)
  end subroutine rrfs_sd_state_register_axis

  subroutine rrfs_sd_state_write_axis(data,Model,Sfc_restart)
    implicit none
    class(rrfs_sd_state_type) :: data
    type(FmsNetcdfDomainFile_t) :: Sfc_restart
    type(GFS_control_type),      intent(in) :: Model
    real(kind_phys) :: fire_aux_data_levels(Model%fire_aux_data_levels)
    integer :: i

    call register_field(Sfc_restart, 'fire_aux_data_levels', 'double', (/'fire_aux_data_levels'/))
    call register_variable_attribute(Sfc_restart, 'fire_aux_data_levels', 'cartesian_axis' ,'Z', str_len=1)

    do i=1,Model%fire_aux_data_levels
      fire_aux_data_levels(i) = i
    enddo

    call write_data(Sfc_restart, 'fire_aux_data_levels', fire_aux_data_levels)
  end subroutine rrfs_sd_state_write_axis

  subroutine rrfs_sd_state_allocate_data(data,Model)
    implicit none
    class(rrfs_sd_state_type) :: data
    type(GFS_control_type),   intent(in) :: Model
    integer :: nx, ny

    call data%deallocate_data

    nx=Model%nx
    ny=Model%ny

    allocate(data%emdust(nx,ny))
    allocate(data%emseas(nx,ny))
    allocate(data%emanoc(nx,ny))
    allocate(data%fhist(nx,ny))
    allocate(data%coef_bb_dc(nx,ny))

    allocate(data%fire_in(nx,ny,Model%fire_aux_data_levels))

  end subroutine rrfs_sd_state_allocate_data

  subroutine rrfs_sd_state_fill_data(data, Model, Sfcprop, Atm_block)
    ! Fills all temporary variables with default values.
    ! Terrible things will happen if you don't call data%allocate_data first.
    ! IMPORTANT: This must match the corresponding code in sfcprop_create in
    ! GFS_typedefs.F90
    implicit none
    class(rrfs_sd_state_type) :: data
    type(GFS_sfcprop_type),   intent(in) :: Sfcprop(:)
    type(GFS_control_type),   intent(in) :: Model
    type(block_control_type), intent(in) :: Atm_block

    integer :: nb, ix, isc, jsc, i, j

    isc = Model%isc
    jsc = Model%jsc

!$omp parallel do default(shared) private(i, j, nb, ix)
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1

        data%emdust(i,j) = 0
        data%emseas(i,j) = 0
        data%emanoc(i,j) = 0
        data%fhist(i,j) = 1.
        data%coef_bb_dc(i,j) = 0

        data%fire_in(i,j,:) = 0
      end do
    end do
  end subroutine rrfs_sd_state_fill_data

  subroutine rrfs_sd_state_register_fields(data,Sfc_restart)
    ! Registers all restart fields needed by the RRFS-SD
    ! Terrible things will happen if you don't call data%allocate_data
    ! and data%register_axes first.
    implicit none
    class(rrfs_sd_state_type) :: data
    type(FmsNetcdfDomainFile_t) :: Sfc_restart

    ! Register 2D fields
    call register_restart_field(Sfc_restart, 'emdust', data%emdust, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)
    call register_restart_field(Sfc_restart, 'emseas', data%emseas, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)
    call register_restart_field(Sfc_restart, 'emanoc', data%emanoc, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)
    call register_restart_field(Sfc_restart, 'fhist', data%fhist, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)
    call register_restart_field(Sfc_restart, 'coef_bb_dc', data%coef_bb_dc, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)

    ! Register 3D field
    call register_restart_field(Sfc_restart, 'fire_in', data%fire_in, &
         dimensions=(/'xaxis_1             ', 'yaxis_1             ', &
                      'fire_aux_data_levels', 'Time                '/), &
         is_optional=.true.)
  end subroutine rrfs_sd_state_register_fields

  subroutine rrfs_sd_state_final(data)
    ! Final routine for rrfs_sd_state_type, called automatically when
    ! an object of that type goes out of scope.  This is a wrapper
    ! around data%deallocate_data() with necessary syntactic
    ! differences.
    implicit none
    type(rrfs_sd_state_type) :: data
    call rrfs_sd_state_deallocate_data(data)
  end subroutine rrfs_sd_state_final

  subroutine rrfs_sd_state_deallocate_data(data)
    ! Deallocates all data used, and nullifies the pointers. The data
    ! object can safely be used again after this call. This is also
    ! the implementation of the rrfs_sd_state_deallocate_data final routine.
    implicit none
    class(rrfs_sd_state_type) :: data

    ! This #define reduces code length by a lot
#define IF_ASSOC_DEALLOC_NULL(var) \
    if(associated(data%var)) then ; \
      deallocate(data%var) ; \
      nullify(data%var) ; \
    endif

    IF_ASSOC_DEALLOC_NULL(emdust)
    IF_ASSOC_DEALLOC_NULL(emseas)
    IF_ASSOC_DEALLOC_NULL(emanoc)
    IF_ASSOC_DEALLOC_NULL(fhist)
    IF_ASSOC_DEALLOC_NULL(coef_bb_dc)

    IF_ASSOC_DEALLOC_NULL(fire_in)

    ! Undefine this to avoid cluttering the cpp scope:
#undef IF_ASSOC_DEALLOC_NULL
  end subroutine rrfs_sd_state_deallocate_data

  subroutine rrfs_sd_state_copy_from_temporaries(data, Model, Sfcprop, Atm_block)
    implicit none
    class(rrfs_sd_state_type) :: data
    type(GFS_sfcprop_type),   intent(in) :: Sfcprop(:)
    type(GFS_control_type),   intent(in) :: Model
    type(block_control_type), intent(in) :: Atm_block

    integer :: nb, ix, isc, jsc, i, j

    isc = Model%isc
    jsc = Model%jsc

!$omp parallel do default(shared) private(i, j, nb, ix)
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1

        Sfcprop(nb)%emdust(ix) = data%emdust(i,j)
        Sfcprop(nb)%emseas(ix) = data%emseas(i,j)
        Sfcprop(nb)%emanoc(ix) = data%emanoc(i,j)
        Sfcprop(nb)%fhist(ix) = data%fhist(i,j)
        Sfcprop(nb)%coef_bb_dc(ix) = data%coef_bb_dc(i,j)

        Sfcprop(nb)%fire_in(ix,:) = data%fire_in(i,j,:)
      enddo
    enddo
  end subroutine rrfs_sd_state_copy_from_temporaries

  subroutine rrfs_sd_state_copy_to_temporaries(data, Model, Sfcprop, Atm_block)
    implicit none
    class(rrfs_sd_state_type) :: data
    type(GFS_sfcprop_type),   intent(in) :: Sfcprop(:)
    type(GFS_control_type),   intent(in) :: Model
    type(block_control_type), intent(in) :: Atm_block

    integer :: nb, ix, isc, jsc, i, j

    isc = Model%isc
    jsc = Model%jsc

!$omp parallel do default(shared) private(i, j, nb, ix)
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1

        data%emdust(i,j) = Sfcprop(nb)%emdust(ix)
        data%emseas(i,j) = Sfcprop(nb)%emseas(ix)
        data%emanoc(i,j) = Sfcprop(nb)%emanoc(ix)
        data%fhist(i,j) = Sfcprop(nb)%fhist(ix)
        data%coef_bb_dc(i,j) = Sfcprop(nb)%coef_bb_dc(ix)

        data%fire_in(i,j,:) = Sfcprop(nb)%fire_in(ix,:)
      enddo
    enddo
  end subroutine rrfs_sd_state_copy_to_temporaries

end module rrfs_sd_io
