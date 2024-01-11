! This program provides a unit test to verify that the variables in
! module_fv3_io_def.F90 are accessible.

program test_fv3_io_def

  use module_fv3_io_def

  implicit none

  integer myint

  num_pes_fcst = 2
  if (num_pes_fcst .ne. 2) stop 1
  wrttasks_per_group = 4
  if (wrttasks_per_group .ne. 4) stop 1
  write_groups = 8
  if (write_groups .ne. 8) stop 1
  n_group = 16
  if (n_group .ne. 16) stop 1
  num_files = 32
  if (num_files .ne. 32) stop 1
  nbdlphys = 64
  if (nbdlphys .ne. 64) stop 1
  iau_offset = 128
  if (iau_offset .ne. 128) stop 1

  lflname_fulltime = .true.
  if (.not. lflname_fulltime) stop 2
  time_unlimited = .true.
  if (.not. time_unlimited) stop 2

  allocate(filename_base(1))
  filename_base(1)="abc"
!  output_file="abc"

!  integer,dimension(:),allocatable     :: lead_wrttask, last_wrttask
!
!  character(len=esmf_maxstr),dimension(:),allocatable :: output_grid
!  integer,dimension(:),allocatable  :: imo,jmo
!  real,dimension(:),allocatable     :: cen_lon, cen_lat
!  real,dimension(:),allocatable     :: lon1, lat1, lon2, lat2, dlon, dlat
!  real,dimension(:),allocatable     :: stdlat1, stdlat2, dx, dy
!  integer,dimension(:),allocatable  :: ideflate, quantize_nsd, zstandard_level
!  character(len=esmf_maxstr),dimension(:),allocatable :: quantize_mode
!  integer,dimension(:),allocatable  :: ichunk2d, jchunk2d, ichunk3d, jchunk3d, kchunk3d


end program test_fv3_io_def
