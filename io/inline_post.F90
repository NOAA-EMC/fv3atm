!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
module inline_post

  use module_fv3_io_def,    only : wrttasks_per_group,filename_base,      &
                                   output_grid
  use write_internal_state, only : wrt_internal_state
  use post_gfs,             only : post_getattr_gfs, post_run_gfs
  use post_regional,        only : post_getattr_regional, post_run_regional

  implicit none

  public  inline_post_run, inline_post_getattr

  contains

  subroutine inline_post_run(wrt_int_state,mypei,mpicomp,lead_write,      &
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
      if(mypei == 0) print *,'inline_post_run, output_grid=',trim(output_grid)
      if(trim(output_grid) == 'gaussian_grid'                &
        .or. trim(output_grid) == 'global_latlon') then
          call post_run_gfs(wrt_int_state, mypei, mpicomp, lead_write, &
                            mynfhr, mynfmin,mynfsec)
      else if( trim(output_grid) == 'regional_latlon'          &
        .or.  trim(output_grid) == 'rotated_latlon'            &
        .or.  trim(output_grid) == 'lambert_conformal') then
      if(mypei == 0) print *,'inline_post_run, call post_run_regional'
          call post_run_regional(wrt_int_state, mypei, mpicomp, lead_write, &
                            mynfhr, mynfmin,mynfsec)
        endif

!
    end subroutine inline_post_run
!
!-----------------------------------------------------------------------
!
    subroutine inline_post_getattr(wrt_int_state)
!
      use esmf
!
      implicit none
!
      type(wrt_internal_state),intent(inout)    :: wrt_int_state
!
        if(trim(output_grid) == 'gaussian_grid'                &
          .or. trim(output_grid) == 'global_latlon') then
            call post_getattr_gfs(wrt_int_state)
        else if( trim(output_grid) == 'regional_latlon'          &
          .or.  trim(output_grid) == 'rotated_latlon'            &
          .or.  trim(output_grid) == 'lambert_conformal') then
            call post_getattr_regional(wrt_int_state)
        endif
!
    end subroutine inline_post_getattr


    end module inline_post
