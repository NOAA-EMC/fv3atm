!> @file
!> @brief fv3 configure variables from model_configure.
!> @author Jun Wang @date 01/2017

!> fv3 configure variables from model_configure.
!>
!> @author Jun Wang @date 01/2017
  module module_fv3_config
  use esmf

  implicit none


  !> ???  
  integer                  :: dt_atmos

  !> ???  
  integer                  :: first_kdt

  !> ???  
  integer                  :: fcst_mpi_comm

  !> ???  
  integer                  :: fcst_ntasks


  !> ???  
  integer                  :: cpl_grid_id

  !> ???  
  logical                  :: cplprint_flag

  !> ???  
  logical                  :: quilting

  !> ???  
  logical                  :: quilting_restart


  !> ???  
  real,dimension(:),allocatable                   :: output_fh

  !> ???  
  character(esmf_maxstr),dimension(:),allocatable :: filename_base

  !> ???  
  character(17)            :: calendar='                 '

  end module module_fv3_config
