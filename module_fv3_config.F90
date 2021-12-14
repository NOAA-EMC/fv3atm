
  module module_fv3_config
!------------------------------------------------------------------------
!
!*** fv3 configure variables from model_configure
!
! revision history
! 01/2017   Jun Wang    Initial code
!
!------------------------------------------------------------------------
!
  use esmf

  implicit none
!
  integer                  :: nfhout, nfhout_hf, nsout, dt_atmos
  integer                  :: first_kdt
  integer                  :: fcst_mpi_comm, fcst_ntasks
!
  logical                  :: cplprint_flag
  logical                  :: quilting, output_1st_tstep_rst
  logical                  :: restart_endfcst
!
  real,dimension(:),allocatable                   :: output_fh
  character(esmf_maxstr),dimension(:),allocatable :: filename_base
  character(17)            :: calendar='                 '
!
  end module module_fv3_config
