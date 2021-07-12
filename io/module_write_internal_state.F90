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
      character(64)   :: output_grid
      type(esmf_grid) :: wrtgrid
!
!-----------------------------
!***  full domain information
!-----------------------------
!
      integer :: im
      integer :: jm
      integer :: lm
!
!-----------------------------
!***  subdomain domain information
!-----------------------------
!
      integer :: lat_start, lon_start
      integer :: lat_end, lon_end
      real    :: latstart, latlast, lonstart, lonlast
      integer,dimension(:),allocatable :: lat_start_wrtgrp
      integer,dimension(:),allocatable :: lat_end_wrtgrp
      real,dimension(:,:),allocatable  :: lonPtr, latPtr
!
!--------------------------
!*** file bundle for output
!--------------------------
      integer :: FBCount
      integer,dimension(:), allocatable           :: ncount_attribs
      integer,dimension(:), allocatable           :: ncount_fields
      character(128),dimension(:),allocatable     :: wrtFB_names
      character(128),dimension(:,:),allocatable   :: field_names
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
      type(ESMF_TimeInterval) :: io_currtimediff
      real                    :: nfhour
      integer                 :: idate(7)
      integer                 :: fdate(7)
!
!-----------------------------------------
!***  I/O direction flags (Read or Write)
!-----------------------------------------
!
      logical :: output_history
      logical :: write_netcdfflag
!
!-----------------------------------------
!***  POST flags and required variables
!-----------------------------------------
!
      logical                  :: write_dopost
      integer                  :: post_nlunit
      character(80)            :: post_namelist
!
      integer                  :: post_maptype
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
