!> @file
!> @brief The internal state of the write component.
!> @author Jun Wang @date Feb, 2017

!> @brief The internal state of the write component.
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
        integer :: im !< output grid global I dimension size
        integer :: jm !< output grid global J dimension size
        integer :: lm !< output grid global L dimension size
        integer :: i_start !< output grid lower bound of I dimension on current PE
        integer :: i_end !< output grid upper bound of I dimension on current PE
        integer :: j_start !< output grid lower bound of J dimension on current PE
        integer :: j_end !< output grid upper bound of J dimension on current PE
        real,dimension(:,:),allocatable  :: lonPtr !< output grid longitudes
        real,dimension(:,:),allocatable  :: latPtr !< output grid latitudes
        integer,dimension(:),allocatable :: i_start_wrtgrp !< I dimension lower bound of all wrire groups
        integer,dimension(:),allocatable :: i_end_wrtgrp !< I dimension upper bound of all wrire groups
        integer,dimension(:),allocatable :: j_start_wrtgrp !< J dimension lower bound of all wrire groups
        integer,dimension(:),allocatable :: j_end_wrtgrp !< J dimension upper bound of all wrire groups
        real    :: latse !< output grid South East corner latitude
        real    :: latnw !< output grid North West corner latitude
        real    :: lonse !< output grid South East corner longitude
        real    :: lonnw !< output grid North West corner longitude
        real    :: latstart !< output grid start latitude
        real    :: latlast !< output grid last latitude
        real    :: lonstart !< output grid start logitude
        real    :: lonlast !< output grid last longitude
      end type output_grid_info

      type wrt_internal_state

!--------------------------------
! pe information and task layout
!--------------------------------
!
      integer :: mype !< MPI rank
      integer :: petcount !< Number of PEs.
!
!--------------------
!*** grid information
!--------------------
      type(esmf_grid) :: wrtgrid !< ESMF output grid

      type(output_grid_info) ,dimension(:), allocatable :: out_grid_info !< Array of output_grid_info for all domains
!
!--------------------------
!*** file bundle for output
!--------------------------
      integer :: FBCount !< Numebr of output ESMF field bundles
!
!-----------------------------------------------------------------------
!***  THE OUTPUT FILE
!-----------------------------------------------------------------------
!
      integer                 :: num_files !< number of output files
!
!-----------------------------------------------------------------------
!***  THE OUTPUT FILE
!-----------------------------------------------------------------------
!
      type(ESMF_FieldBundle),dimension(:),allocatable  :: wrtFB !< ESMF write field bundles
!
!-------------------------------------
!***  Times used in history filenames
!-------------------------------------
!
      type(ESMF_Time)         :: io_basetime !< ESMF clock's starting time
      integer                 :: idate(7) !< Forecast initial time
      integer                 :: fdate(7) !< Forecast current time
!
!-----------------------------------------
!***  I/O direction flags (Read or Write)
!-----------------------------------------
!
      logical :: output_history !< True if history output is requested
!
!-----------------------------------------
!***  POST flags and required variables
!-----------------------------------------
!
      logical                  :: write_dopost !< True if inline post is requested
      character(80)            :: post_namelist !< filename of the inline post namelist
!
      integer                  :: fhzero !< hours between clearing of diagnostic buckets
      integer                  :: ntrac !< number of tracers
      integer                  :: ncld !< number of hydrometeors
      integer                  :: nsoil !< number of soil layers
      integer                  :: imp_physics !< choice of microphysics scheme
      integer                  :: dtp !< physics timestep
      real,dimension(:),allocatable :: ak !< a parameter for sigma pressure level calculations
      real,dimension(:),allocatable :: bk !< b parameter for sigma pressure level calculations
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
        type(wrt_internal_state),pointer :: write_int_state  !< Write grid component internal state
      end type write_wrap

!-----------------------------------------------------------
!
      end module write_internal_state
!
!-----------------------------------------------------------
