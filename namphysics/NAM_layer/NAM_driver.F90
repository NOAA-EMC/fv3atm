module GFS_driver

  use machine,                  only: kind_phys
  use GFS_typedefs,             only: GFS_init_type,                       &
                                      GFS_statein_type, GFS_stateout_type, &
                                      GFS_sfcprop_type, GFS_coupling_type, &
                                      GFS_control_type, GFS_grid_type,     &
                                      GFS_tbd_type,     GFS_cldprop_type,  &
                                      GFS_radtend_type, GFS_diag_type
  use module_radiation_driver,  only: GFS_radiation_driver
  use module_solver_grid_comp,  only: solver_run

  implicit none

  private

  integer, allocatable :: blksz(:)

!----------------
! Public entities
!----------------
  public  GFS_initialize              !< GFS initialization routine
  public  GFS_time_vary_step          !< perform operations needed prior radiation or physics
  public  GFS_radiation_driver        !< radiation_driver (was grrad)
  public  solver_run
  public  GFS_stochastic_driver       !< stochastic physics


  CONTAINS
!*******************************************************************************************


!--------------
! GFS initialze
!--------------
  subroutine GFS_initialize (Model, Statein, Stateout, Sfcprop,    &
                             Coupling, Grid, Tbd, Cldprop, Radtend, & 
                             Diag, Init_parm)

       USE MODULE_RA_RRTM,      ONLY: RRTM_INIT
       USE MODULE_RA_GFDL,      ONLY: GFDL_INIT

    !--- interface variables
    type(GFS_control_type),   intent(inout) :: Model
    type(GFS_statein_type),   intent(inout) :: Statein(:)
    type(GFS_stateout_type),  intent(inout) :: Stateout(:)
    type(GFS_sfcprop_type),   intent(inout) :: Sfcprop(:)
    type(GFS_coupling_type),  intent(inout) :: Coupling(:)
    type(GFS_grid_type),      intent(inout) :: Grid(:)
    type(GFS_tbd_type),       intent(inout) :: Tbd(:)
    type(GFS_cldprop_type),   intent(inout) :: Cldprop(:)
    type(GFS_radtend_type),   intent(inout) :: Radtend(:)
    type(GFS_diag_type),      intent(inout) :: Diag(:)
    type(GFS_init_type),      intent(in)    :: Init_parm

    !--- local variables
    integer :: nb,nblks,l
    real(kind=kind_phys) :: gmt,pt
    real(kind=kind_phys), allocatable :: smid(:)
    real(kind=kind_phys), allocatable :: sful(:)
    real(kind=kind_phys), parameter   :: p_ref = 101325.0d0
    logical :: LSASHAL=.false.


    nblks = size(Init_parm%blksz)
    allocate (blksz(nblks))
    blksz(:) = Init_parm%blksz(:)

    !--- set control properties (including namelist read)
    call Model%init (Init_parm%nlunit, Init_parm%fn_nml,           &
                     Init_parm%me, Init_parm%master,               &
                     Init_parm%logunit, Init_parm%isc,             &
                     Init_parm%jsc, Init_parm%nx, Init_parm%ny,    &
                     Init_parm%levs, Init_parm%cnx, Init_parm%cny, &
                     Init_parm%gnx, Init_parm%gny,                 &
                     Init_parm%dt_dycore, Init_parm%dt_phys,       &
                     Init_parm%iau_offset, Init_parm%bdat,         &
                     Init_parm%cdat, Init_parm%tracer_names,       &
                     Init_parm%input_nml_file, Init_parm%tile_num, &
                     Init_parm%blksz)

    do nb = 1,nblks
      call Statein  (nb)%create (Init_parm%blksz(nb), Model)
      call Stateout (nb)%create (Init_parm%blksz(nb), Model)
      call Sfcprop  (nb)%create (Init_parm%blksz(nb), Model)
      call Coupling (nb)%create (Init_parm%blksz(nb), Model)
      call Grid     (nb)%create (Init_parm%blksz(nb), Model)
      call Tbd      (nb)%create (Init_parm%blksz(nb), Model)
      call Cldprop  (nb)%create (Init_parm%blksz(nb), Model)
      call Radtend  (nb)%create (Init_parm%blksz(nb), Model)
      !--- internal representation of diagnostics
      call Diag     (nb)%create (Init_parm%blksz(nb), Model)
    enddo

    !--- populate the grid components
    call GFS_grid_populate (Grid, Init_parm%xlon, Init_parm%xlat, Init_parm%area)

    allocate(smid(Model%levr  ))
    allocate(sful(Model%levr+1))
    do l=1,Model%levr+1
      sful(l) = (Init_parm%ak(l) + Init_parm%bk(l) * p_ref - Init_parm%ak(Model%levr+1)) &
               / (p_ref - Init_parm%ak(Model%levr+1))
    enddo
    do l=1,Model%levr
      smid(l)=(sful(l)+sful(L+1))*0.5
    enddo
    pt=Init_parm%ak(Model%levr+1) + Init_parm%bk(Model%levr+1)

    !--- initialize NAM radiation
    IF(TRIM(Model%LONGWAVE)=='gfdl' .AND. TRIM(Model%SHORTWAVE)=='gfdl')THEN

      GMT=REAL(Model%jdat(5))
      CALL GFDL_INIT(SFUL,SMID,PT*1.0E-3                              &
                    ,Model%jdat(1),Model%jdat(2),Model%jdat(3),GMT    &
                    ,Model%CO2TF,Model%levr)

    ELSEIF(TRIM(Model%LONGWAVE)=='rrtm' .AND. TRIM(Model%SHORTWAVE)=='rrtm')THEN

      call rad_initialize_nmmb                                        &
        (SFUL,Model%levr,Model%ICTM,Model%ISOL,Model%ICO2,Model%IAER  &
        ,Model%IAER_MDL,Model%IALB,Model%IEMS,Model%NTCW              &
        ,Model%NP3D,Model%NTOZ,Model%IOVR_SW,Model%IOVR_LW            &
        ,Model%ISUBC_SW,Model%ISUBC_LW,Model%ICLIQ_SW,Model%ICICE_SW  &
        ,Model%ICLIQ_LW,Model%ICICE_LW,LSASHAL,Model%CRICK_PROOF      &
        ,Model%CCNORM,Model%NORAD_PRECIP,Model%IFLIP,Model%me)
      CALL RRTM_INIT(SMID,PT*1.0E-3,Model%levr)

    ELSE

      write(0,*)' Unknown radiation scheme',TRIM(Model%LONGWAVE),TRIM(Model%SHORTWAVE)
      STOP
    ENDIF

    deallocate (smid)
    deallocate (sful)
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  end subroutine GFS_initialize


!-------------------------------------------------------------------------
  subroutine GFS_time_vary_step (Model, Statein, Stateout, Sfcprop, Coupling, & 
                                 Grid, Tbd, Cldprop, Radtend, Diag)

    implicit none

    !--- interface variables
    type(GFS_control_type),   intent(inout) :: Model
    type(GFS_statein_type),   intent(inout) :: Statein(:)
    type(GFS_stateout_type),  intent(inout) :: Stateout(:)
    type(GFS_sfcprop_type),   intent(inout) :: Sfcprop(:)
    type(GFS_coupling_type),  intent(inout) :: Coupling(:)
    type(GFS_grid_type),      intent(inout) :: Grid(:)
    type(GFS_tbd_type),       intent(inout) :: Tbd(:)
    type(GFS_cldprop_type),   intent(inout) :: Cldprop(:)
    type(GFS_radtend_type),   intent(inout) :: Radtend(:)
    type(GFS_diag_type),      intent(inout) :: Diag(:)
    !--- local variables
    real(kind=kind_phys) :: rinc(5)
    real(kind=kind_phys) :: sec

    rinc(1:5)   = 0
    call w3difdat(Model%jdat,Model%idat,4,rinc)
    sec = rinc(4)
    Model%kdt   = nint((sec + Model%dtp)/Model%dtp)

  end subroutine GFS_time_vary_step

!-------------------------------------------------------------------------
  subroutine GFS_stochastic_driver (Model, Statein, Stateout, Sfcprop, Coupling, &
                                    Grid, Tbd, Cldprop, Radtend, Diag)

    implicit none

    !--- interface variables
    type(GFS_control_type),   intent(in   ) :: Model
    type(GFS_statein_type),   intent(in   ) :: Statein
    type(GFS_stateout_type),  intent(in   ) :: Stateout
    type(GFS_sfcprop_type),   intent(in   ) :: Sfcprop
    type(GFS_coupling_type),  intent(inout) :: Coupling
    type(GFS_grid_type),      intent(in   ) :: Grid
    type(GFS_tbd_type),       intent(in   ) :: Tbd
    type(GFS_cldprop_type),   intent(in   ) :: Cldprop
    type(GFS_radtend_type),   intent(in   ) :: Radtend
    type(GFS_diag_type),      intent(inout) :: Diag
    !--- local variables
  end subroutine GFS_stochastic_driver



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!------------------
! GFS_grid_populate
!------------------
  subroutine GFS_grid_populate (Grid, xlon, xlat, area)

    implicit none

    type(GFS_grid_type)              :: Grid(:)
    real(kind=kind_phys), intent(in) :: xlon(:,:)
    real(kind=kind_phys), intent(in) :: xlat(:,:)
    real(kind=kind_phys), intent(in) :: area(:,:)

    !--- local variables
    integer :: nb, ix, blksz, i, j

    blksz = size(Grid(1)%xlon)

    nb = 1
    ix = 0
    do j = 1,size(xlon,2)
      do i = 1,size(xlon,1)
        ix=ix+1
        if (ix .gt. blksz) then
          nb = nb + 1
          ix = 1
        endif
        Grid(nb)%xlon(ix)   = xlon(i,j)
        Grid(nb)%xlat(ix)   = xlat(i,j)
        Grid(nb)%area(ix)   = area(i,j)
        Grid(nb)%dx(ix)     = sqrt(area(i,j))
      enddo
    enddo

  end subroutine GFS_grid_populate

end module GFS_driver

