!========================================!
      module module_radiation_driver        !
!........................................!
      use GFS_typedefs,              only: GFS_statein_type,             &
                                           GFS_stateout_type,            &
                                           GFS_sfcprop_type,             &
                                           GFS_coupling_type,            &
                                           GFS_control_type,             &
                                           GFS_grid_type,                &
                                           GFS_tbd_type,                 &
                                           GFS_cldprop_type,             &
                                           GFS_radtend_type,             &
                                           GFS_diag_type
!
      implicit   none
      private
      public radinit, radupdate, GFS_radiation_driver
! =================
      contains
! =================
!-----------------------------------
      subroutine radinit
      end subroutine radinit
!-----------------------------------
      subroutine radupdate
      end subroutine radupdate
!-----------------------------------
      subroutine GFS_radiation_driver                             &
         (Model, Statein, Stateout, Sfcprop, Coupling, Grid, Tbd, &
          Cldprop, Radtend, Diag)

      implicit none

      type(GFS_control_type),         intent(in)    :: Model
      type(GFS_statein_type),         intent(in)    :: Statein
      type(GFS_stateout_type),        intent(inout) :: Stateout
      type(GFS_sfcprop_type),         intent(in)    :: Sfcprop
      type(GFS_coupling_type),        intent(inout) :: Coupling
      type(GFS_grid_type),            intent(in)    :: Grid
      type(GFS_tbd_type),             intent(in)    :: Tbd
      type(GFS_cldprop_type),         intent(in)    :: Cldprop
      type(GFS_radtend_type),         intent(inout) :: Radtend
      type(GFS_diag_type),            intent(inout) :: Diag
        
      end subroutine GFS_radiation_driver
!----------------------------------------
      end module module_radiation_driver !
!========================================!
