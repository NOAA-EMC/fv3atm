!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
module post_gfs

  use module_fv3_io_def,    only : wrttasks_per_group,filename_base
  use write_internal_state, only : wrt_internal_state

  implicit none

  public  post_run_gfs, post_getattr_gfs

  contains

  subroutine post_run_gfs(wrt_int_state,mypei,mpicomp,lead_write,      &
             mynfhr,mynfmin,mynfsec)
!
!  revision history:
!     Jul 2019    J. Wang      create interface to run inline post for FV3
!
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(wrt_internal_state),intent(in)       :: wrt_int_state
      integer,intent(in)                        :: mypei
      integer,intent(in)                        :: mpicomp
      integer,intent(in)                        :: lead_write
      integer,intent(in)                        :: mynfhr
      integer,intent(in)                        :: mynfmin
      integer,intent(in)                        :: mynfsec
!
      print *,'in stub post_run_gfs - not supported on this machine, return'
!
    end subroutine post_run_gfs
!
!-----------------------------------------------------------------------
!
    subroutine post_getattr_gfs(wrt_int_state, fldbundle)
!
      use esmf
!
      implicit none
!
      type(wrt_internal_state),intent(inout)    :: wrt_int_state
      type(ESMF_FieldBundle), intent(in)        :: fldbundle
!
!
      print *,'in stub post_getattr_gfs - not supported on this machine, return'
!
    end subroutine post_getattr_gfs


    end module post_gfs
