module GFS_init

  use machine,                  only: kind_phys
  use GFS_typedefs,             only: GFS_init_type,                       &
                                      GFS_statein_type, GFS_stateout_type, &
                                      GFS_sfcprop_type, GFS_coupling_type, &
                                      GFS_control_type, GFS_grid_type,     &
                                      GFS_tbd_type,     GFS_cldprop_type,  &
                                      GFS_radtend_type, GFS_diag_type,     &
                                      GFS_interstitial_type

  implicit none

  private

!----------------
! Public entities
!----------------
  public  GFS_initialize              !< GFS initialization routine

  CONTAINS
!*******************************************************************************************


!--------------
! GFS initialze
!--------------
  subroutine GFS_initialize (Model, Statein, Stateout, Sfcprop,     &
                             Coupling, Grid, Tbd, Cldprop, Radtend, & 
                             Diag, Interstitial, communicator,      &
                             ntasks, Init_parm)

#ifdef _OPENMP
    use omp_lib
#endif

    !--- interface variables
    type(GFS_control_type),      intent(inout) :: Model
    type(GFS_statein_type),      intent(inout) :: Statein(:)
    type(GFS_stateout_type),     intent(inout) :: Stateout(:)
    type(GFS_sfcprop_type),      intent(inout) :: Sfcprop(:)
    type(GFS_coupling_type),     intent(inout) :: Coupling(:)
    type(GFS_grid_type),         intent(inout) :: Grid(:)
    type(GFS_tbd_type),          intent(inout) :: Tbd(:)
    type(GFS_cldprop_type),      intent(inout) :: Cldprop(:)
    type(GFS_radtend_type),      intent(inout) :: Radtend(:)
    type(GFS_diag_type),         intent(inout) :: Diag(:)
    type(GFS_interstitial_type), intent(inout) :: Interstitial(:)
    integer,                     intent(in)    :: communicator
    integer,                     intent(in)    :: ntasks
    type(GFS_init_type),         intent(in)    :: Init_parm

    !--- local variables
    integer :: nb
    integer :: nblks
    integer :: nt
    integer :: nthrds
    logical :: non_uniform_blocks
    integer :: ix

    nblks = size(Init_parm%blksz)

#ifdef _OPENMP
    nthrds = omp_get_max_threads()
#else
    nthrds = 1
#endif

    !--- set control properties (including namelist read)
    call Model%init (Init_parm%nlunit, Init_parm%fn_nml,           &
                     Init_parm%me, Init_parm%master,               &
                     Init_parm%logunit, Init_parm%isc,             &
                     Init_parm%jsc, Init_parm%nx, Init_parm%ny,    &
                     Init_parm%levs, Init_parm%cnx, Init_parm%cny, &
                     Init_parm%gnx, Init_parm%gny,                 &
                     Init_parm%dt_dycore, Init_parm%dt_phys,       &
                     Init_parm%iau_offset, Init_parm%bdat,         &
                     Init_parm%cdat, Init_parm%nwat,               &
                     Init_parm%tracer_names,                       &
                     Init_parm%tracer_types,                       &
                     Init_parm%input_nml_file, Init_parm%tile_num, &
                     Init_parm%blksz, Init_parm%ak, Init_parm%bk,  &
                     Init_parm%restart, Init_parm%hydrostatic,     &
                     communicator, ntasks, nthrds)

    do nb = 1,nblks
      ix = Init_parm%blksz(nb)
      call Statein  (nb)%create (ix, Model)
      call Stateout (nb)%create (ix, Model)
      call Sfcprop  (nb)%create (ix, Model)
      call Coupling (nb)%create (ix, Model)
      call Grid     (nb)%create (ix, Model)
      call Tbd      (nb)%create (ix, Model)
      call Cldprop  (nb)%create (ix, Model)
      call Radtend  (nb)%create (ix, Model)
!--- internal representation of diagnostics
      call Diag     (nb)%create (ix, Model)
    enddo

! This logic deals with non-uniform block sizes for CCPP. When non-uniform block sizes
! are used, it is required that only the last block has a different (smaller) size than
! all other blocks. This is the standard in FV3. If this is the case, set non_uniform_blocks
! to .true. and initialize nthreads+1 elements of the interstitial array. The extra element
! will be used by the thread that runs over the last, smaller block.
    if (minval(Init_parm%blksz)==maxval(Init_parm%blksz)) then
       non_uniform_blocks = .false.
    elseif (all(minloc(Init_parm%blksz)==(/size(Init_parm%blksz)/))) then
       non_uniform_blocks = .true.
    else
       write(0,'(2a)') 'For non-uniform blocksizes, only the last element ', &
                       'in Init_parm%blksz can be different from the others'
       stop
    endif

! Initialize the Interstitial data type in parallel so that
! each thread creates (touches) its Interstitial(nt) first.
!$OMP parallel do default (shared) &
!$OMP            schedule (static,1) &
!$OMP            private  (nt)
    do nt=1,nthrds
      call Interstitial (nt)%create (maxval(Init_parm%blksz), Model)
    enddo
!$OMP end parallel do

    if (non_uniform_blocks) then
      call Interstitial (nthrds+1)%create (Init_parm%blksz(nblks), Model)
    end if

    !--- populate the grid components
    call GFS_grid_populate (Grid, Init_parm%xlon, Init_parm%xlat, Init_parm%area)

  end subroutine GFS_initialize

!------------------
! GFS_grid_populate
!------------------
  subroutine GFS_grid_populate (Grid, xlon, xlat, area)
    use physcons, only: pi => con_pi

    implicit none

    type(GFS_grid_type)              :: Grid(:)
    real(kind=kind_phys), intent(in) :: xlon(:,:)
    real(kind=kind_phys), intent(in) :: xlat(:,:)
    real(kind=kind_phys), intent(in) :: area(:,:)
    real(kind=kind_phys), parameter  :: rad2deg = 180.0_kind_phys/pi

    !--- local variables
    integer :: nb, ix, blksz, i, j

    blksz = size(Grid(1)%xlon)

    nb = 1
    ix = 0
    do j = 1,size(xlon,2)
      do i = 1,size(xlon,1)
        ix=ix+1
        if (ix > blksz) then
          nb = nb + 1
          ix = 1
        endif
        Grid(nb)%xlon(ix)   = xlon(i,j)
        Grid(nb)%xlat(ix)   = xlat(i,j)
        Grid(nb)%xlat_d(ix) = xlat(i,j) * rad2deg
        Grid(nb)%xlon_d(ix) = xlon(i,j) * rad2deg
        Grid(nb)%sinlat(ix) = sin(Grid(nb)%xlat(ix))
        Grid(nb)%coslat(ix) = sqrt(1.0_kind_phys - Grid(nb)%sinlat(ix)*Grid(nb)%sinlat(ix))
        Grid(nb)%area(ix)   = area(i,j)
        Grid(nb)%dx(ix)     = sqrt(area(i,j))
      enddo
    enddo

  end subroutine GFS_grid_populate

end module GFS_init
