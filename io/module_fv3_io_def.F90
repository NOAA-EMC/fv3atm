  module module_fv3_io_def
!
!*** fv3 io related configration variables
!
! revision history
! 01/2017   Jun Wang    Initial code
!
!------------------------------------------------------------------------
!
  implicit none
!
  integer           :: num_pes_fcst
  integer           :: wrttasks_per_group, write_groups
  integer           :: n_group
  logical           :: write_nemsioflip
  logical           :: write_fsyncflag
  integer           :: num_files
  character(255)    :: app_domain
  character(255)    :: output_grid
  character(255)    :: output_file
  integer           :: imo,jmo
  integer           :: nbdlphys
  integer           :: nsout_io, iau_offset, ideflate, nbits
  real              :: cen_lon, cen_lat, lon1, lat1, lon2, lat2, dlon, dlat
  real              :: stdlat1, stdlat2, dx, dy
  character(255),dimension(:),allocatable :: filename_base
!
  integer,dimension(:),allocatable     :: lead_wrttask, last_wrttask
!
  end module module_fv3_io_def

