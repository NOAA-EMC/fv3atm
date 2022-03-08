module module_fv3_io_def
!
!*** fv3 io related configration variables
!
! revision history
! 01/2017   Jun Wang    Initial code
!
!------------------------------------------------------------------------
!
  use esmf, only     : esmf_maxstr
  implicit none

  integer           :: num_pes_fcst
  integer           :: wrttasks_per_group, write_groups
  integer           :: n_group
  integer           :: num_files
  integer           :: nbdlphys
  integer           :: nsout_io, iau_offset
  logical           :: lflname_fulltime

  character(len=esmf_maxstr),dimension(:),allocatable :: filename_base
  character(len=esmf_maxstr),dimension(:),allocatable :: output_file

  integer,dimension(:),allocatable     :: lead_wrttask, last_wrttask

  character(len=esmf_maxstr),dimension(:),allocatable :: output_grid
  integer,dimension(:),allocatable  :: imo,jmo
  real,dimension(:),allocatable     :: cen_lon, cen_lat
  real,dimension(:),allocatable     :: lon1, lat1, lon2, lat2, dlon, dlat
  real,dimension(:),allocatable     :: stdlat1, stdlat2, dx, dy
  integer,dimension(:),allocatable  :: ideflate, nbits
  integer,dimension(:),allocatable  :: ichunk2d, jchunk2d, ichunk3d, jchunk3d, kchunk3d

end module module_fv3_io_def
