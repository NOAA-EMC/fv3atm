!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
module inline_post

  use module_fv3_io_def,    only : wrttasks_per_group,filename_base
  use write_internal_state, only : wrt_internal_state

  implicit none

  public  inline_post_run, inline_post_getattr

  contains

  subroutine inline_post_run(wrt_int_state,mypei,mpicomp,lead_write,      &
             mynfhr,mynfmin,mynfsec)
!
!  revision history:
!     Oct 2020    J. Wang      create interface to run inline post
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
      print *,'in stub inline_post_run - not supported on this machine, return'
!
    end subroutine inline_post_run
!
!-----------------------------------------------------------------------
!
    subroutine inline_post_getattr(wrt_int_state)
!
      implicit none
!
      type(wrt_internal_state),intent(inout)    :: wrt_int_state
!
!
      print *,'in stub inline_post_getattr - not supported on this machine, return'
!
    end subroutine inline_post_getattr


    end module inline_post
