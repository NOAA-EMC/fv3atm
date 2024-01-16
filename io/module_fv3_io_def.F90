!> @file
!> @brief fv3 I/O related configration variables.
!> @author Jun Wang @date 01/2017

!> fv3 I/O related configration variables.
!>
!> @author Jun Wang @date 01/2017
module module_fv3_io_def
  use esmf, only     : esmf_maxstr
  implicit none

  !> Number of processors used in forecast run.
  integer           :: num_pes_fcst

  !> Number of write tasks per group.  
  integer           :: wrttasks_per_group

  !> ???  
  integer           :: write_groups 

  !> ???  
  integer           :: n_group 

  !> ???  
  integer           :: num_files 

  !> ???  
  integer           :: nbdlphys 

  !> ???  
  integer           :: iau_offset 

  !> ???  
  logical           :: lflname_fulltime 

  !> ???  
  logical           :: time_unlimited 


  !> ???  
  character(len=esmf_maxstr),dimension(:),allocatable :: filename_base  

  !> ???  
  character(len=esmf_maxstr),dimension(:),allocatable :: output_file 


  !> ???  
  integer,dimension(:),allocatable     :: lead_wrttask 

  !> ???  
  integer,dimension(:),allocatable     :: last_wrttask 


  !> ???  
  character(len=esmf_maxstr),dimension(:),allocatable :: output_grid 

  !> ???  
  integer,dimension(:),allocatable  :: imo 

  !> ???  
  integer,dimension(:),allocatable  :: jmo 

  !> ???  
  real,dimension(:),allocatable     :: cen_lon 

  !> ???  
  real,dimension(:),allocatable     :: cen_lat 

  !> ???  
  real,dimension(:),allocatable     :: lon1 

  !> ???  
  real,dimension(:),allocatable     :: lat1 

  !> ???  
  real,dimension(:),allocatable     :: lon2 

  !> ???  
  real,dimension(:),allocatable     :: lat2 

  !> ???  
  real,dimension(:),allocatable     :: dlon 

  !> ???  
  real,dimension(:),allocatable     :: dlat 

  !> ???  
  real,dimension(:),allocatable     :: stdlat1 

  !> ???  
  real,dimension(:),allocatable     :: stdlat2 

  !> ???  
  real,dimension(:),allocatable     :: dx 

  !> ???  
  real,dimension(:),allocatable     :: dy 

  !> Deflate level to use, 0 means no deflate.
  integer,dimension(:),allocatable  :: ideflate

  !> Number of significant digits for lossy compression.
  integer,dimension(:),allocatable  :: quantize_nsd

  !> Zstandard compression level, 0 means no zstandard compression.
  integer,dimension(:),allocatable  :: zstandard_level

  !> Quantize mode to use for lossy compression.
  character(len=esmf_maxstr),dimension(:),allocatable :: quantize_mode

  !> Chunk size in i dimension for 2D data.
  integer,dimension(:),allocatable  :: ichunk2d

  !> Chunk size in j dimension for 2D data.
  integer,dimension(:),allocatable  :: jchunk2d

  !> Chunk size in i dimension for 3D data.
  integer,dimension(:),allocatable  :: ichunk3d

  !> Chunk size in j dimension for 3D data.
  integer,dimension(:),allocatable  :: jchunk3d

  !> Chunk size in k dimension for 3D data.
  integer,dimension(:),allocatable  :: kchunk3d

end module module_fv3_io_def
