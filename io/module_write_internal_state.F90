!> @file 
!> @brief Internal state of the write component.
!> @author J. Wang/G. Theurich @date Feb, 2017

!> @brief Internal state of the write component.
!>
!> @author Jun Wang @date Feb, 2017
      module write_internal_state
      use esmf
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type output_grid_info
        integer :: im !< ???
        integer :: jm !< ???
        integer :: lm !< ???
        integer :: i_start !< ???
        integer :: i_end !< ???
        integer :: j_start !< ???
        integer :: j_end !< ???
        real,dimension(:,:),allocatable  :: lonPtr !< ???
        real,dimension(:,:),allocatable  :: latPtr !< ???
        integer,dimension(:),allocatable :: i_start_wrtgrp !< ???
        integer,dimension(:),allocatable :: i_end_wrtgrp !< ???
        integer,dimension(:),allocatable :: j_start_wrtgrp !< ???
        integer,dimension(:),allocatable :: j_end_wrtgrp !< ???
        real    :: latse !< ???
        real    :: latnw !< ???
        real    :: lonse !< ???
        real    :: lonnw !< ???
        real    :: latstart !< ???
        real    :: latlast !< ???
        real    :: lonstart !< ???
        real    :: lonlast !< ???
      end type output_grid_info

      type wrt_internal_state

!--------------------------------
! pe information and task layout
!--------------------------------
!
      integer :: mype !< ???
      integer :: petcount !< ???
!
!--------------------
!*** grid information
!--------------------
      type(esmf_grid) :: wrtgrid !< ???

      type(output_grid_info) ,dimension(:), allocatable :: out_grid_info !< ???
!
!--------------------------
!*** file bundle for output
!--------------------------
      integer :: FBCount !< ???
!
!-----------------------------------------------------------------------
!***  THE OUTPUT FILE
!-----------------------------------------------------------------------
!
      integer                 :: num_files !< ???
!
!-----------------------------------------------------------------------
!***  THE OUTPUT FILE
!-----------------------------------------------------------------------
!
      type(ESMF_FieldBundle),dimension(:),allocatable  :: wrtFB !< ???
!
!-------------------------------------
!***  Times used in history filenames
!-------------------------------------
!
      type(ESMF_Time)         :: io_basetime !< ???
      integer                 :: idate(7) !< ???
      integer                 :: fdate(7) !< ???
!
!-----------------------------------------
!***  I/O direction flags (Read or Write)
!-----------------------------------------
!
      logical :: output_history !< ???
!
!-----------------------------------------
!***  POST flags and required variables
!-----------------------------------------
!
      logical                  :: write_dopost !< ???
      character(80)            :: post_namelist !< ???
!
      integer                  :: fhzero !< ???
      integer                  :: ntrac !< ???
      integer                  :: ncld !< ???
      integer                  :: nsoil !< ???
      integer                  :: imp_physics !< ???
      integer                  :: dtp !< ???
      real,dimension(:),allocatable :: ak !< ???
      real,dimension(:),allocatable :: bk !< ???
!-----------------------------------------------------------------------
!
      end type wrt_internal_state
!

      !> This state is supported by c pointers but not f90 pointers
      !> therefore we need this wrap.
      type write_wrap
        type(wrt_internal_state),pointer :: write_int_state  !< ???
      end type write_wrap

!-----------------------------------------------------------
!
      end module write_internal_state
!
!-----------------------------------------------------------
