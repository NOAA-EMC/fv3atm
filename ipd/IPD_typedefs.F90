module IPD_typedefs

!---------------------------------------------------------
! Physics/Radiation types used to create various IPD types
!---------------------------------------------------------
  use physics_abstraction_layer, only: IPD_control_type => control_type,     &
                                       IPD_init_type    => init_type,        &
                                       IPD_restart_type => restart_type,     &
                                       IPD_diag_type    => diagnostic_type,  &
                                       IPD_kind_phys    => kind_phys
#ifdef CCPP
  use physics_abstraction_layer, only: IPD_interstitial_type => interstitial_type
#endif

!---------------------------------------------------------
! Physics/Radiation types used to create the IPD_data_type
!---------------------------------------------------------
  use physics_abstraction_layer, only: statein_type,  stateout_type,         &
                                       sfcprop_type,  coupling_type,         &
                                       grid_type,     tbd_type,              &
                                       cldprop_type,  radtend_type,          &
                                       intdiag_type
#ifdef CCPP
  use physics_abstraction_layer, only: IPD_data_type    => data_type
#endif

!-------------------------------------------------
! Physics/Radiation routines to pass to IPD_driver
!-------------------------------------------------
  use physics_abstraction_layer,  only: initialize,          &
                                        diagnostic_populate, &
                                        restart_populate

#ifndef CCPP
!-------------------------------------------------------
!  IPD_data_type 
!    container of physics data types that can be blocked
!-------------------------------------------------------
  type IPD_data_type
    type(statein_type)  :: Statein
    type(stateout_type) :: Stateout
    type(sfcprop_type)  :: Sfcprop
    type(coupling_type) :: Coupling
    type(grid_type)     :: Grid
    type(tbd_type)      :: Tbd
    type(cldprop_type)  :: Cldprop
    type(radtend_type)  :: Radtend
    type(intdiag_type)  :: Intdiag
  end type IPD_data_type
#endif


!------------------------------------------------------
!  IPD function procedures 
!    definitions for scalar(0d) and vector(1d) versions
!------------------------------------------------------
  abstract interface
    subroutine IPD_func0d_proc (Control, Statein, Stateout,  &
                                Sfcprop, Coupling, Grid,     &
                                Tbd, Cldprop, Radtend,       &
                                Intdiag)
      import :: IPD_control_type, statein_type, stateout_type,    &
                sfcprop_type, coupling_type, grid_type, tbd_type, & 
                cldprop_type, radtend_type, intdiag_type
      type(IPD_control_type), intent(inout) :: Control
      type(statein_type),     intent(inout) :: Statein
      type(stateout_type),    intent(inout) :: Stateout
      type(sfcprop_type),     intent(inout) :: Sfcprop
      type(coupling_type),    intent(inout) :: Coupling
      type(grid_type),        intent(inout) :: Grid
      type(tbd_type),         intent(inout) :: Tbd
      type(cldprop_type),     intent(inout) :: Cldprop
      type(radtend_type),     intent(inout) :: Radtend
      type(intdiag_type),     intent(inout) :: Intdiag
    end subroutine IPD_func0d_proc

    subroutine IPD_func1d_proc (Control, Statein, Stateout,  &
                                Sfcprop, Coupling, Grid,     &
                                Tbd, Cldprop, Radtend,       &
                                Intdiag)
      import :: IPD_control_type, statein_type, stateout_type,    &
                sfcprop_type, coupling_type, grid_type, tbd_type, & 
                cldprop_type, radtend_type, intdiag_type
      type(IPD_control_type), intent(inout) :: Control
      type(statein_type),     intent(inout) :: Statein(:)
      type(stateout_type),    intent(inout) :: Stateout(:)
      type(sfcprop_type),     intent(inout) :: Sfcprop(:)
      type(coupling_type),    intent(inout) :: Coupling(:)
      type(grid_type),        intent(inout) :: Grid(:)
      type(tbd_type),         intent(inout) :: Tbd(:)
      type(cldprop_type),     intent(inout) :: Cldprop(:)
      type(radtend_type),     intent(inout) :: Radtend(:)
      type(intdiag_type),     intent(inout) :: Intdiag(:)
    end subroutine IPD_func1d_proc
  end interface


!------------------------------------------------
! SAMPLE var_subtype
!   pointers to two and three dimensional objects
!------------------------------------------------
!  type var_subtype
!    real(kind=kind_phys), pointer :: var2p(:)   => null()  !< 2D data saved in packed format [dim(ix)]
!    real(kind=kind_phys), pointer :: var3p(:,:) => null()  !< 3D data saved in packed format [dim(ix,levs)]
!  end type var_subtype
!   
!--------------------------------------------------
! SAMPLE restart_type to import as IPD_restart_type
!   data necessary for reproducible restarts
!--------------------------------------------------
!  type restart_type
!    integer           :: num2d                    !< current number of registered 2D restart variables
!    integer           :: num3d                    !< current number of registered 3D restart variables
!    character(len=32), allocatable :: name2d(:)   !< variable name as it will appear in the restart file
!    character(len=32), allocatable :: name3d(:)   !< variable name as it will appear in the restart file
!    type(var_subtype), allocatable :: data(:,:)   !< holds pointers to data in packed format (allocated to (nblks,max(2d/3dfields))
!  end type restart_type
!
!--------------------------------------------------
! SAMPLE diagnostic_type to import as IPD_diag_type
!   fields targetted as diagnostic output
!--------------------------------------------------
!  type diag_type
!    character(len=32)     :: name           !< variable name in source
!    character(len=32)     :: output_name    !< output name for variable
!    character(len=32)     :: mod_name       !< module name (e.g. physics, radiation, etc)
!    character(len=32)     :: file_name      !< output file name for variable
!    character(len=128)    :: desc           !< long description of field
!    character(len=32)     :: unit           !< units associated with fields
!    character(len=32)     :: type_stat_proc !< type of statistic processing:
!                                            !< average, accumulation, maximal, minimal, etc.
!    character(len=32)     :: level_type     !< vertical level of the field
!    integer               :: level          !< vertical level(s)
!    real(kind=kind_phys)  :: cnvfac         !< conversion factors to output in specified units
!    real(kind=kind_phys)  :: zhour          !< forecast hour when bucket was last emptied for statistical processing
!    real(kind=kind_phys)  :: fcst_hour      !< current forecast hour (same as fhour)
!    type(var_subtype), allocatable :: data(:) !< holds pointers to data in packed format (allocated to nblks)
!  end type diag_type


!------------------------
! IPD public declarations
!------------------------
  public IPD_kind_phys
  public IPD_control_type
  public IPD_data_type
  public IPD_restart_type
  public IPD_diag_type
  public IPD_init_type
#ifdef CCPP
  public IPD_interstitial_type
#endif


!-----------------------------------
! public declarations for IPD_driver
!-----------------------------------
  public initialize
  public diagnostic_populate
  public restart_populate

  CONTAINS
!*******************************************************************************************

end module IPD_typedefs
