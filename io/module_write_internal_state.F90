!-----------------------------------------------------------------------
      module write_internal_state
!
!-----------------------------------------------------------------------
!***  the internal state of the write component.
!-----------------------------------------------------------------------
!***
!***  revision history
!***
!       Feb 2017:  J. Wang - Initial code
!
!-----------------------------------------------------------------------
!
      use esmf
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type output_grid_info
        integer :: im, jm, lm
        integer :: i_start,i_end, j_start,j_end
        real,dimension(:,:),allocatable  :: lonPtr, latPtr
        integer,dimension(:),allocatable :: i_start_wrtgrp, i_end_wrtgrp, j_start_wrtgrp, j_end_wrtgrp
        real    :: latse, latnw, lonse, lonnw
        real    :: latstart, latlast, lonstart, lonlast
      end type output_grid_info

      type wrt_internal_state

!--------------------------------
! pe information and task layout
!--------------------------------
!
      integer :: mype
      integer :: petcount
!
!--------------------
!*** grid information
!--------------------
      type(esmf_grid) :: wrtgrid

      type(output_grid_info) ,dimension(:), allocatable :: out_grid_info
!
!--------------------------
!*** file bundle for output
!--------------------------
      integer :: FBCount
      character(128),dimension(:),allocatable     :: wrtFB_names
!
!-----------------------------------------------------------------------
!***  THE OUTPUT FILE
!-----------------------------------------------------------------------
!
      integer                 :: num_files
!
!-----------------------------------------------------------------------
!***  THE OUTPUT FILE
!-----------------------------------------------------------------------
!
      type(ESMF_FieldBundle),dimension(:),allocatable  :: wrtFB
!
!-------------------------------------
!***  Times used in history filenames
!-------------------------------------
!
      type(ESMF_Time)         :: io_basetime
      integer                 :: idate(7)
      integer                 :: fdate(7)
!
!-----------------------------------------
!***  I/O direction flags (Read or Write)
!-----------------------------------------
!
      logical :: output_history
!
!-----------------------------------------
!***  POST flags and required variables
!-----------------------------------------
!
      logical                  :: write_dopost
      character(80)            :: post_namelist
!
      integer                  :: fhzero
      integer                  :: ntrac
      integer                  :: ncld
      integer                  :: nsoil
      integer                  :: imp_physics
      integer                  :: dtp
      real,dimension(:),allocatable :: ak,bk
!-----------------------------------------------------------------------
!
      end type wrt_internal_state
!
!-----------------------------------------------------------------------
!***  THIS STATE IS SUPPORTED BY C POINTERS BUT NOT F90 POINTERS
!***  THEREFORE WE NEED THIS WRAP.
!-----------------------------------------------------------
!
      type write_wrap
        type(wrt_internal_state),pointer :: write_int_state
      end type write_wrap

!-----------------------------------------------------------
!
      end module write_internal_state
!
!-----------------------------------------------------------
