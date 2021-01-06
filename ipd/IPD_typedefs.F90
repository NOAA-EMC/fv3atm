module IPD_typedefs

!---------------------------------------------------------
! Physics/Radiation types used to create various IPD types
!---------------------------------------------------------
  use physics_abstraction_layer, only: IPD_control_type      => control_type,      &
                                       IPD_init_type         => init_type,         &
                                       IPD_restart_type      => restart_type,      &
                                       IPD_diag_type         => diagnostic_type,   &
                                       IPD_interstitial_type => interstitial_type, &
                                       IPD_data_type         => data_type,         &
                                       IPD_kind_phys         => kind_phys

!---------------------------------------------------------
! Physics/Radiation types used to create the IPD_data_type
!---------------------------------------------------------
  use physics_abstraction_layer, only: statein_type,  stateout_type,         &
                                       sfcprop_type,  coupling_type,         &
                                       grid_type,     tbd_type,              &
                                       cldprop_type,  radtend_type,          &
                                       intdiag_type

!-------------------------------------------------
! Physics/Radiation routines to pass to IPD_driver
!-------------------------------------------------
  use physics_abstraction_layer,  only: initialize,          &
                                        diagnostic_populate, &
                                        restart_populate

!------------------------
! IPD public declarations
!------------------------
  public IPD_kind_phys
  public IPD_control_type
  public IPD_data_type
  public IPD_restart_type
  public IPD_diag_type
  public IPD_init_type
  public IPD_interstitial_type

!-----------------------------------
! public declarations for IPD_driver
!-----------------------------------
  public initialize
  public diagnostic_populate
  public restart_populate

  CONTAINS
!*******************************************************************************************

end module IPD_typedefs
