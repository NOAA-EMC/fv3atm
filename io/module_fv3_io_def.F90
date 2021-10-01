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
!
  integer           :: num_pes_fcst
  integer           :: wrttasks_per_group, write_groups
  integer           :: n_group
  integer           :: num_files
  character(len=esmf_maxstr)    :: app_domain
  character(len=esmf_maxstr)    :: output_grid
  integer           :: imo,jmo
  integer           :: ichunk2d,jchunk2d,ichunk3d,jchunk3d,kchunk3d
  integer           :: nbdlphys
  integer           :: nsout_io, iau_offset, ideflate, nbits
  logical           :: lflname_fulltime
  real              :: cen_lon, cen_lat, lon1, lat1, lon2, lat2, dlon, dlat
  real              :: stdlat1, stdlat2, dx, dy
  character(len=esmf_maxstr),dimension(:),allocatable :: filename_base
  character(len=esmf_maxstr),dimension(:),allocatable :: output_file
!
  integer,dimension(:),allocatable     :: lead_wrttask, last_wrttask
!
  end module module_fv3_io_def

