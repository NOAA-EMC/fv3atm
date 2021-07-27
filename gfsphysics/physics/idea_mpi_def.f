! VAY-2015
!   collection of the mpi-related variables and subs
!   for the WAM-mpi debugging      .... 
!    mpi_WAM_quit(iret, message)   with iulog-unit output 
!                                  can be => NEMS-out
      MODULE IDEA_MPI_def
!
! or use mpi_def.f from  /src/atmos/gsm/phys
!/scratch1/portfolios/NCEPDEV/swpc/save/Valery.Yudin/NEMS/wam_april_2015/src/atmos/gsm/libutil/module_gfs_mpi_def.f90
      implicit none
      include 'mpif.h' 
!      use machine,   ONLY: KIND_io4, KIND_ior
!      use module_gfs_mpi_def             
      integer, parameter :: masterproc =0
      integer :: mpi_err
      integer :: mpi_id
      integer :: mpi_num_procs  
      INTEGER :: info
      INTEGER :: MPI_COMM_ALL
!
      integer stat(MPI_STATUS_SIZE)
      INTEGER ::  MC_COMP
      INTEGER ::  MC_IO
      INTEGER ::  MPI_COMM_ALL_DUP
!
!
      INTEGER :: num_pes
      INTEGER :: nnum_pes_fcst
      INTEGER :: first_fcst_pe
      INTEGER :: last_fcst_pe 

      INTEGER :: mpi_comm_inter
      INTEGER :: mpi_comm_comp
      INTEGER :: mpi_inter_b
!=======================================================      
!   mpi_def.f in ../gsm/phys  -lines
!
!      integer MPI_R_IO, MPI_R_MPI, MPI_R_DEF, MPI_A_DEF
!     &,       MPI_R_IO_R,MPI_R_MPI_R
!      PARAMETER (MPI_R_IO  =MPI_REAL4)
!      PARAMETER (MPI_R_IO_R=MPI_REAL8)
!
!      PARAMETER (MPI_R_DEF=MPI_REAL8)
!      PARAMETER (MPI_A_DEF=MPI_REAL8)
!========================================================
      CONTAINS
!
      SUBROUTINE mpi_WAM_quit(iret, message)
      use IDEA_IO_UNITS, only : iulog
      implicit none
!
      integer :: iret

      character(len=*) message
      write(iulog,*) trim(message),' in CALL  mpi_WAM_quit ',iret
!
!      CALL MPI_ABORT(MPI_COMM_WORLD,rc,ierr)
!
      CALL MPI_ABORT(MC_COMP, iret, info)
      END SUBROUTINE mpi_WAM_quit
!
      END MODULE IDEA_MPI_def
