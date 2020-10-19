module IPD_driver

  use IPD_typedefs,               only: IPD_kind_phys,     IPD_init_type,    &
                                        IPD_control_type,  IPD_data_type,    &
                                        IPD_diag_type,     IPD_restart_type, &
                                        IPD_func0d_proc,   IPD_func1d_proc,  &
                                        initialize,                          &
                                        diagnostic_populate,                 &
                                        restart_populate
#ifdef CCPP
  use IPD_typedefs,               only: IPD_interstitial_type
#endif

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
  public IPD_step 

  CONTAINS
!*******************************************************************************************


!----------------
!  IPD Initialize 
!----------------
#ifdef CCPP
  subroutine IPD_initialize (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, &
                             IPD_Interstitial, communicator, ntasks, IPD_init_parm)
#else
  subroutine IPD_initialize (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, IPD_init_parm)
#endif
    type(IPD_control_type), intent(inout) :: IPD_Control
    type(IPD_data_type),    intent(inout) :: IPD_Data(:)
    type(IPD_diag_type),    intent(inout) :: IPD_Diag(:)
    type(IPD_restart_type), intent(inout) :: IPD_Restart
#ifdef CCPP
    type(IPD_interstitial_type), intent(inout) :: IPD_Interstitial(:)
    integer, intent(in)                   :: communicator
    integer, intent(in)                   :: ntasks
#endif
    type(IPD_init_type),    intent(in)    :: IPD_init_parm

    !--- initialize the physics suite
    call initialize (IPD_Control, IPD_Data(:)%Statein, IPD_Data(:)%Stateout,      &
                     IPD_Data(:)%Sfcprop, IPD_Data(:)%Coupling, IPD_Data(:)%Grid, &
                     IPD_Data(:)%Tbd, IPD_Data(:)%Cldprop, IPD_Data(:)%Radtend,   &
#ifdef CCPP
                     IPD_Data(:)%Intdiag, IPD_Interstitial(:), communicator,      &
                     ntasks, IPD_init_parm)
#else
                     IPD_Data(:)%Intdiag, IPD_init_parm)
#endif


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

!----------------------------------------------------------
!  IPD step
!    runs the given routine/function pointed to by IPD_func 
!----------------------------------------------------------
  subroutine IPD_step (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, IPD_func0d, IPD_func1d)
    type(IPD_control_type),     intent(inout) :: IPD_Control
    type(IPD_data_type),        intent(inout) :: IPD_Data(:)
    type(IPD_diag_type),        intent(inout) :: IPD_Diag(:)
    type(IPD_restart_type),     intent(inout) :: IPD_Restart
    procedure(IPD_func0d_proc), intent(in), optional, pointer :: IPD_func0d
    procedure(IPD_func1d_proc), intent(in), optional, pointer :: IPD_func1d

    if (size(IPD_Data,1) == 1 .and. PRESENT(IPD_func0d)) then
      call IPD_func0d (IPD_Control, IPD_Data(1)%Statein, IPD_Data(1)%Stateout,      &
                       IPD_Data(1)%Sfcprop, IPD_Data(1)%Coupling, IPD_Data(1)%Grid, &
                       IPD_Data(1)%Tbd, IPD_Data(1)%Cldprop, IPD_Data(1)%Radtend,   &
                       IPD_Data(1)%Intdiag)
    else
      call IPD_func1d (IPD_Control, IPD_Data(:)%Statein, IPD_Data(:)%Stateout,      &
                       IPD_Data(:)%Sfcprop, IPD_Data(:)%Coupling, IPD_Data(:)%Grid, &
                       IPD_Data(:)%Tbd, IPD_Data(:)%Cldprop, IPD_Data(:)%Radtend,   &
                       IPD_Data(:)%Intdiag)
    endif

  end subroutine IPD_step

end module IPD_driver
