module IPD_driver

  use IPD_typedefs,               only: IPD_kind_phys,     IPD_init_type,    &
                                        IPD_control_type,  IPD_data_type,    &
                                        IPD_diag_type,     IPD_restart_type, &
                                        IPD_interstitial_type,               &
                                        initialize,                          &
                                        diagnostic_populate,                 &
                                        restart_populate

  implicit none

!------------------------------------------------------!
!  IPD containers                                      !
!------------------------------------------------------!
!  type(GFS_control_type)              :: IPD_Control  !
!  type(IPD_data_type)     allocatable :: IPD_Data(:)  !
!  type(IPD_diag_type),                :: IPD_Diag(:)  !
!  type(IPD_restart_type),             :: IPD_Restart  !
!------------------------------------------------------!

!----------------
! Public Entities
!----------------
! functions
  public IPD_initialize, IPD_initialize_rst

  CONTAINS
!*******************************************************************************************


!----------------
!  IPD Initialize 
!----------------
  subroutine IPD_initialize (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, &
                             IPD_Interstitial, communicator, ntasks, IPD_init_parm)
    type(IPD_control_type), intent(inout) :: IPD_Control
    type(IPD_data_type),    intent(inout) :: IPD_Data(:)
    type(IPD_diag_type),    intent(inout) :: IPD_Diag(:)
    type(IPD_restart_type), intent(inout) :: IPD_Restart
    type(IPD_interstitial_type), intent(inout) :: IPD_Interstitial(:)
    integer, intent(in)                   :: communicator
    integer, intent(in)                   :: ntasks
    type(IPD_init_type),    intent(in)    :: IPD_init_parm

    !--- initialize the physics suite
    call initialize (IPD_Control, IPD_Data(:)%Statein, IPD_Data(:)%Stateout,      &
                     IPD_Data(:)%Sfcprop, IPD_Data(:)%Coupling, IPD_Data(:)%Grid, &
                     IPD_Data(:)%Tbd, IPD_Data(:)%Cldprop, IPD_Data(:)%Radtend,   &
                     IPD_Data(:)%Intdiag, IPD_Interstitial(:), communicator,      &
                     ntasks, IPD_init_parm)

    !--- populate/associate the Diag container elements
    call diagnostic_populate (IPD_Diag, IPD_Control, IPD_Data%Statein, IPD_Data%Stateout,   &
                                        IPD_Data%Sfcprop, IPD_Data%Coupling, IPD_Data%Grid, &
                                        IPD_Data%Tbd, IPD_Data%Cldprop, IPD_Data%Radtend,   &
                                        IPD_Data%Intdiag, IPD_init_parm)

  end subroutine IPD_initialize

!----------------
!  IPD Initialize phase_rst
!----------------
  subroutine IPD_initialize_rst (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, IPD_init_parm)
    type(IPD_control_type), intent(inout) :: IPD_Control
    type(IPD_data_type),    intent(inout) :: IPD_Data(:)
    type(IPD_diag_type),    intent(inout) :: IPD_Diag(:)
    type(IPD_restart_type), intent(inout) :: IPD_Restart
    type(IPD_init_type),    intent(in)    :: IPD_init_parm

    !--- allocate and populate/associate the Restart container elements
    call restart_populate (IPD_Restart, IPD_Control, IPD_Data%Statein, IPD_Data%Stateout,   &
                                        IPD_Data%Sfcprop, IPD_Data%Coupling, IPD_Data%Grid, &
                                        IPD_Data%Tbd, IPD_Data%Cldprop, IPD_Data%Radtend,   &
                                        IPD_Data%Intdiag, IPD_init_parm, IPD_Diag)

  end subroutine IPD_initialize_rst

end module IPD_driver
