module module_block_data

  use ESMF,              only: ESMF_KIND_R8, ESMF_SUCCESS, &
                               ESMF_RC_PTR_NOTALLOC, ESMF_RC_VAL_OUTOFRANGE
  use GFS_typedefs,      only: kind_phys
  use block_control_mod, only: block_control_type

  implicit none

  interface block_data_copy
    module procedure block_copy_1d_to_2d_r8
    module procedure block_copy_2d_to_2d_r8
    module procedure block_copy_2d_to_3d_r8
    module procedure block_copy_3d_to_3d_r8
    module procedure block_copy_1dslice_to_2d_r8
    module procedure block_copy_3dslice_to_3d_r8
  end interface block_data_copy

  interface block_data_fill
    module procedure block_fill_2d_r8
    module procedure block_fill_3d_r8
  end interface block_data_fill

  interface block_data_copy_or_fill
    module procedure block_copy_or_fill_1d_to_2d_r8
    module procedure block_copy_or_fill_2d_to_3d_r8
    module procedure block_copy_or_fill_1dslice_to_2d_r8
  end interface block_data_copy_or_fill

  interface block_data_combine_fractions
    module procedure block_combine_frac_1d_to_2d_r8
  end interface block_data_combine_fractions

  interface block_atmos_copy
    module procedure block_array_copy_2d_to_2d_r8
    module procedure block_array_copy_3d_to_3d_r8
    module procedure block_array_copy_3dslice_to_3d_r8
  end interface block_atmos_copy

  private

  public :: block_atmos_copy

  public :: block_data_copy
  public :: block_data_fill
  public :: block_data_copy_or_fill
  public :: block_data_combine_fractions

contains

  ! -- copy: 1D to 2D

  subroutine block_copy_1d_to_2d_r8(destin_ptr, source_ptr, block, block_index, scale_factor, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind_phys),           pointer     :: source_ptr(:)
    type(block_control_type),  intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind_phys), optional, intent(in)  :: scale_factor
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer         :: localrc
    integer         :: i, ib, ix, j, jb
    real(kind_phys) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
      factor = 1._kind_phys
      if (present(scale_factor)) factor = scale_factor
      do ix = 1, block%blksz(block_index)
        ib = block%index(block_index)%ii(ix)
        jb = block%index(block_index)%jj(ix)
        i = ib - block%isc + 1
        j = jb - block%jsc + 1
        destin_ptr(i,j) = factor * source_ptr(ix)
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_1d_to_2d_r8

  ! -- copy: 1D slice to 2D

  subroutine block_copy_1dslice_to_2d_r8(destin_ptr, source_ptr, slice, block, block_index, scale_factor, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind_phys),           pointer     :: source_ptr(:,:)
    integer,                   intent(in)  :: slice
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind_phys), optional, intent(in)  :: scale_factor
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer         :: localrc
    integer         :: i, ib, ix, j, jb
    real(kind_phys) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
      localrc = ESMF_RC_VAL_OUTOFRANGE
      if (slice > 0 .and. slice <= size(source_ptr, dim=2)) then
        factor = 1._kind_phys
        if (present(scale_factor)) factor = scale_factor
        do ix = 1, block%blksz(block_index)
          ib = block%index(block_index)%ii(ix)
          jb = block%index(block_index)%jj(ix)
          i = ib - block%isc + 1
          j = jb - block%jsc + 1
          destin_ptr(i,j) = factor * source_ptr(ix,slice)
        enddo
        localrc = ESMF_SUCCESS
      end if
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_1dslice_to_2d_r8

  ! -- copy: 2D to 3D

  subroutine block_copy_2d_to_3d_r8(destin_ptr, source_ptr, block, block_index, scale_factor, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real(kind_phys),           pointer     :: source_ptr(:,:)
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind_phys), optional, intent(in)  :: scale_factor
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer         :: localrc
    integer         :: i, ib, ix, j, jb, k
    real(kind_phys) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
      factor = 1._kind_phys
      if (present(scale_factor)) factor = scale_factor
      do k = 1, size(source_ptr, dim=2)
        do ix = 1, block%blksz(block_index)
          ib = block%index(block_index)%ii(ix)
          jb = block%index(block_index)%jj(ix)
          i = ib - block%isc + 1
          j = jb - block%jsc + 1
          destin_ptr(i,j,k) = factor * source_ptr(ix,k)
        enddo
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_2d_to_3d_r8

  ! -- copy: 2D to 2D

  subroutine block_copy_2d_to_2d_r8(destin_ptr, source_ptr, block, block_index, scale_factor, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind_phys),           pointer     :: source_ptr(:,:)
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind_phys), optional, intent(in)  :: scale_factor
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer         :: localrc
    integer         :: i, ib, ix, j, jb
    real(kind_phys) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
      factor = 1._kind_phys
      if (present(scale_factor)) factor = scale_factor
      do ix = 1, block%blksz(block_index)
        ib = block%index(block_index)%ii(ix)
        jb = block%index(block_index)%jj(ix)
        i = ib - block%isc + 1
        j = jb - block%jsc + 1
        destin_ptr(i,j) = factor * source_ptr(ib,jb)
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_2d_to_2d_r8

  subroutine block_array_copy_2d_to_2d_r8(destin_ptr, source_arr, block, block_index, scale_factor, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real,                      intent(in)  :: source_arr(:,:)
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real,    optional,         intent(in)  :: scale_factor
    integer, optional,         intent(out) :: rc

    ! -- local variables
    integer         :: localrc
    integer         :: i, ib, ix, j, jb
    real            :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr)) then
      factor = 1._kind_phys
      if (present(scale_factor)) factor = scale_factor
      do ix = 1, block%blksz(block_index)
        ib = block%index(block_index)%ii(ix)
        jb = block%index(block_index)%jj(ix)
        i = ib - block%isc + 1
        j = jb - block%jsc + 1
        destin_ptr(i,j) = factor * source_arr(ib,jb)
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_array_copy_2d_to_2d_r8

  ! -- copy: 3D to 3D

  subroutine block_copy_3d_to_3d_r8(destin_ptr, source_ptr, block, block_index, scale_factor, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real(kind_phys),           pointer     :: source_ptr(:,:,:)
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind_phys), optional, intent(in)  :: scale_factor
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer         :: localrc
    integer         :: i, ib, ix, j, jb, k
    real(kind_phys) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
      factor = 1._kind_phys
      if (present(scale_factor)) factor = scale_factor
      do k = 1, size(source_ptr, dim=3)
        do ix = 1, block%blksz(block_index)
          ib = block%index(block_index)%ii(ix)
          jb = block%index(block_index)%jj(ix)
          i = ib - block%isc + 1
          j = jb - block%jsc + 1
          destin_ptr(i,j,k) = factor * source_ptr(ib,jb,k)
        enddo
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_3d_to_3d_r8

  subroutine block_array_copy_3d_to_3d_r8(destin_ptr, source_arr, block, block_index, scale_factor, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real,                      intent(in)  :: source_arr(:,:,:)
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real,    optional,         intent(in)  :: scale_factor
    integer, optional,         intent(out) :: rc

    ! -- local variables
    integer         :: localrc
    integer         :: i, ib, ix, j, jb, k
    real            :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr)) then
      factor = 1._kind_phys
      if (present(scale_factor)) factor = scale_factor
      do k = 1, size(source_arr, dim=3)
        do ix = 1, block%blksz(block_index)
          ib = block%index(block_index)%ii(ix)
          jb = block%index(block_index)%jj(ix)
          i = ib - block%isc + 1
          j = jb - block%jsc + 1
          destin_ptr(i,j,k) = factor * source_arr(ib,jb,k)
        enddo
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_array_copy_3d_to_3d_r8

  ! -- copy: 3D slice to 3D

  subroutine block_copy_3dslice_to_3d_r8(destin_ptr, source_ptr, slice, block, block_index, scale_factor, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real(kind_phys),           pointer     :: source_ptr(:,:,:,:)
    integer,                   intent(in)  :: slice
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real(kind_phys), optional, intent(in)  :: scale_factor
    integer,         optional, intent(out) :: rc

    ! -- local variables
    integer         :: localrc
    integer         :: i, ib, ix, j, jb, k
    real(kind_phys) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. associated(source_ptr)) then
      localrc = ESMF_RC_VAL_OUTOFRANGE
      if (slice > 0 .and. slice <= size(source_ptr, dim=4)) then
        factor = 1._kind_phys
        if (present(scale_factor)) factor = scale_factor
        do k = 1, size(source_ptr, dim=3)
          do ix = 1, block%blksz(block_index)
            ib = block%index(block_index)%ii(ix)
            jb = block%index(block_index)%jj(ix)
            i = ib - block%isc + 1
            j = jb - block%jsc + 1
            destin_ptr(i,j,k) = factor * source_ptr(ib,jb,k,slice)
          enddo
        enddo
        localrc = ESMF_SUCCESS
      end if
    end if

    if (present(rc)) rc = localrc

  end subroutine block_copy_3dslice_to_3d_r8

  subroutine block_array_copy_3dslice_to_3d_r8(destin_ptr, source_arr, slice, block, block_index, scale_factor, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real,                      intent(in)  :: source_arr(:,:,:,:)
    integer,                   intent(in)  :: slice
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    real,    optional,         intent(in)  :: scale_factor
    integer, optional,         intent(out) :: rc

    ! -- local variables
    integer         :: localrc
    integer         :: i, ib, ix, j, jb, k
    real            :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr)) then
      localrc = ESMF_RC_VAL_OUTOFRANGE
      if (slice > 0 .and. slice <= size(source_arr, dim=4)) then
        factor = 1._kind_phys
        if (present(scale_factor)) factor = scale_factor
        do k = 1, size(source_arr, dim=3)
          do ix = 1, block%blksz(block_index)
            ib = block%index(block_index)%ii(ix)
            jb = block%index(block_index)%jj(ix)
            i = ib - block%isc + 1
            j = jb - block%jsc + 1
            destin_ptr(i,j,k) = factor * source_arr(ib,jb,k,slice)
          enddo
        enddo
        localrc = ESMF_SUCCESS
      end if
    end if

    if (present(rc)) rc = localrc

  end subroutine block_array_copy_3dslice_to_3d_r8

  ! -- fill: 2D

  subroutine block_fill_2d_r8(destin_ptr, fill_value, block, block_index, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(ESMF_KIND_R8),        intent(in)  :: fill_value
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    integer, optional,         intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: i, ib, ix, j, jb

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr)) then
      do ix = 1, block%blksz(block_index)
        ib = block%index(block_index)%ii(ix)
        jb = block%index(block_index)%jj(ix)
        i = ib - block%isc + 1
        j = jb - block%jsc + 1
        destin_ptr(i,j) = fill_value
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_fill_2d_r8

  ! -- fill: 3D

  subroutine block_fill_3d_r8(destin_ptr, fill_value, block, block_index, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real(ESMF_KIND_R8),        intent(in)  :: fill_value
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    integer, optional,         intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: i, ib, ix, j, jb, k

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr)) then
      do k = 1, size(destin_ptr, dim=3)
        do ix = 1, block%blksz(block_index)
          ib = block%index(block_index)%ii(ix)
          jb = block%index(block_index)%jj(ix)
          i = ib - block%isc + 1
          j = jb - block%jsc + 1
          destin_ptr(i,j,k) = fill_value
        enddo
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_fill_3d_r8

  ! -- copy/fill: 1D to 2D

  subroutine block_copy_or_fill_1d_to_2d_r8(destin_ptr, source_ptr, fill_value, block, block_index, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind_phys),           pointer     :: source_ptr(:)
    real(ESMF_KIND_R8),        intent(in)  :: fill_value
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    integer, optional,         intent(out) :: rc

    ! -- begin
    if (present(rc)) rc = ESMF_RC_PTR_NOTALLOC

    if (associated(destin_ptr)) then
      if (associated(source_ptr)) then
        call block_copy_1d_to_2d_r8(destin_ptr, source_ptr, block, block_index, rc=rc)
      else
        call block_fill_2d_r8(destin_ptr, fill_value, block, block_index, rc=rc)
      end if
    end if

  end subroutine block_copy_or_fill_1d_to_2d_r8

  ! -- copy/fill: 1D slice to 2D

  subroutine block_copy_or_fill_1dslice_to_2d_r8(destin_ptr, source_ptr, slice, fill_value, block, block_index, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind_phys),           pointer     :: source_ptr(:,:)
    integer,                   intent(in)  :: slice
    real(ESMF_KIND_R8),        intent(in)  :: fill_value
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    integer, optional,         intent(out) :: rc

    ! -- begin
    if (present(rc)) rc = ESMF_RC_PTR_NOTALLOC

    if (associated(destin_ptr)) then
      if (associated(source_ptr)) then
        call block_copy_1dslice_to_2d_r8(destin_ptr, source_ptr, slice, block, block_index, rc=rc)
      else
        call block_fill_2d_r8(destin_ptr, fill_value, block, block_index, rc=rc)
      end if
    end if

  end subroutine block_copy_or_fill_1dslice_to_2d_r8

  ! -- copy/fill: 2D to 3D

  subroutine block_copy_or_fill_2d_to_3d_r8(destin_ptr, source_ptr, fill_value, block, block_index, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:,:)
    real(kind_phys),           pointer     :: source_ptr(:,:)
    real(ESMF_KIND_R8),        intent(in)  :: fill_value
    type (block_control_type), intent(in)  :: block
    integer,                   intent(in)  :: block_index
    integer, optional,         intent(out) :: rc

    ! -- begin
    if (present(rc)) rc = ESMF_RC_PTR_NOTALLOC

    if (associated(destin_ptr)) then
      if (associated(source_ptr)) then
        call block_copy_2d_to_3d_r8(destin_ptr, source_ptr, block, block_index, rc=rc)
      else
        call block_fill_3d_r8(destin_ptr, fill_value, block, block_index, rc=rc)
      end if
    end if

  end subroutine block_copy_or_fill_2d_to_3d_r8

  ! -- combine: 1D to 2D

  subroutine block_combine_frac_1d_to_2d_r8(destin_ptr, fract1_ptr, fract2_ptr, block, block_index, rc)

    ! -- arguments
    real(ESMF_KIND_R8),        pointer     :: destin_ptr(:,:)
    real(kind_phys),           pointer     :: fract1_ptr(:)
    real(kind_phys),           pointer     :: fract2_ptr(:)
    type(block_control_type),  intent(in)  :: block
    integer,                   intent(in)  :: block_index
    integer, optional,         intent(out) :: rc

    ! -- local variables
    integer         :: localrc
    integer         :: i, ib, ix, j, jb
    real(kind_phys) :: factor

    ! -- begin
    localrc = ESMF_RC_PTR_NOTALLOC
    if (associated(destin_ptr) .and. &
        associated(fract1_ptr) .and. associated(fract2_ptr)) then
      do ix = 1, block%blksz(block_index)
        ib = block%index(block_index)%ii(ix)
        jb = block%index(block_index)%jj(ix)
        i = ib - block%isc + 1
        j = jb - block%jsc + 1
        destin_ptr(i,j) = fract1_ptr(ix) * (1._kind_phys - fract2_ptr(ix))
      enddo
      localrc = ESMF_SUCCESS
    end if

    if (present(rc)) rc = localrc

  end subroutine block_combine_frac_1d_to_2d_r8

end module module_block_data
