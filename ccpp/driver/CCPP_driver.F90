module CCPP_driver

#ifdef STATIC
! For static builds, the ccpp_physics_{init,run,finalize} calls
! are not pointing to code in the CCPP framework, but to auto-generated
! ccpp_suite_cap and ccpp_group_*_cap modules behind a ccpp_static_api
  use ccpp_api,           only: ccpp_t,                              &
                                ccpp_init,                           &
                                ccpp_finalize,                       &
                                ccpp_initialized
  use ccpp_static_api,    only: ccpp_physics_init,                   &
                                ccpp_physics_run,                    &
                                ccpp_physics_finalize
#else
  use ccpp_api,           only: ccpp_t,                              &
                                ccpp_init,                           &
                                ccpp_finalize,                       &
                                ccpp_physics_init,                   &
                                ccpp_physics_run,                    &
                                ccpp_physics_finalize,               &
                                ccpp_field_add,                      &
                                ccpp_initialized,                    &
                                ccpp_error
#endif

  use CCPP_data,          only: cdata_tile,                          &
                                cdata_domain,                        &
                                cdata_block,                         &
                                ccpp_suite,                          &
                                GFS_control

#ifndef STATIC
! Begin include auto-generated list of modules for ccpp
#include "ccpp_modules_slow_physics.inc"
! End include auto-generated list of modules for ccpp
  use iso_c_binding,      only: c_loc
#endif

  implicit none

!--------------------------------------------------------!
!  Pointer to CCPP containers defined in CCPP_data       !
!--------------------------------------------------------!
  type(ccpp_t), pointer :: cdata => null()

!--------------------------------------------------------!
!  Flag for non-uniform block sizes (last block smaller) !
!  and number of OpenMP threads (with special thread     !
!  number nthrdsX in case of non-uniform block sizes)    !
!--------------------------------------------------------!
  logical :: non_uniform_blocks
  integer :: nthrds, nthrdsX

!----------------
! Public Entities
!----------------
! functions
  public CCPP_step
! module variables
  public non_uniform_blocks

  CONTAINS
!*******************************************************************************************

!-------------------------------
!  CCPP step
!-------------------------------
  subroutine CCPP_step (step, nblks, ierr)

#ifdef OPENMP
    use omp_lib
#endif

    implicit none

    character(len=*),         intent(in)  :: step
    integer,          target, intent(in)  :: nblks
    integer,                  intent(out) :: ierr
    ! Local variables
    integer :: nb, nt, ntX
    integer :: ierr2

    ierr = 0

    if (trim(step)=="init") then

      ! Get and set number of OpenMP threads (module
      ! variable) that are available to run physics
#ifdef OPENMP
      nthrds = omp_get_max_threads()
#else
      nthrds = 1
#endif

      ! For non-uniform blocksizes, we use index nthrds+1
      ! for the interstitial data type with different length
      if (non_uniform_blocks) then
        nthrdsX = nthrds+1
      else
        nthrdsX = nthrds
      end if

      !--- Initialize CCPP framework, if cdata_tile for fast physics
      !    is already initialized, use its suite to avoid reading the
      !    SDF multiple times. Initialize cdata for the entire domain:
      if (ccpp_initialized(cdata_tile)) then
        call ccpp_init(trim(ccpp_suite), cdata_domain, ierr=ierr, cdata_target=cdata_tile)
      else
        call ccpp_init(trim(ccpp_suite), cdata_domain, ierr=ierr)
      end if
      if (ierr/=0) then
        write(0,*) 'An error occurred in ccpp_init'
        return
      end if

      ! For physics running over the entire domain, block and thread
      ! number are not used; set to safe values
      cdata_domain%blk_no = 1
      cdata_domain%thrd_no = 1

#ifndef STATIC
      ! Associate cdata with the global/domain cdata structure;
      ! this is needed because ccpp_fields.inc uses 'cdata' in
      ! its ccpp_field_add statements.
      associate_domain: associate (cdata => cdata_domain)
! Begin include auto-generated list of calls to ccpp_field_add
#include "ccpp_fields_slow_physics.inc"
! End include auto-generated list of calls to ccpp_field_add
      end associate associate_domain
#endif

      ! Allocate cdata structures for blocks and threads
      allocate(cdata_block(1:nblks,1:nthrdsX))

      ! Loop over all blocks and threads
      do nt=1,nthrdsX
        do nb=1,nblks
          !--- Initialize CCPP framework for blocks/threads, use suite from scalar cdata
          !    to avoid reading the SDF multiple times. If cdata_tile is initialized, use
          !    this version (since cdata_domain is just a copy), otherwise use cdata_domain
          if (ccpp_initialized(cdata_tile)) then
            call ccpp_init(trim(ccpp_suite), cdata_block(nb,nt), ierr=ierr, cdata_target=cdata_tile)
          else
            call ccpp_init(trim(ccpp_suite), cdata_block(nb,nt), ierr=ierr, cdata_target=cdata_domain)
          end if
          if (ierr/=0) then
            write(0,'(2(a,i4))') "An error occurred in ccpp_init for block ", nb, " and thread ", nt
            return
          end if

          ! Assign the correct block and thread numbers
          cdata_block(nb,nt)%blk_no = nb
          cdata_block(nb,nt)%thrd_no = nt

#ifndef STATIC
          ! Associate cdata with the cdata structure for this block
          ! and thread; this is needed because ccpp_fields.inc uses
          ! 'cdata' in its ccpp_field_add statements.
          associate_block: associate (cdata => cdata_block(nb,nt))
! Begin include auto-generated list of calls to ccpp_field_add
#include "ccpp_fields_slow_physics.inc"
! End include auto-generated list of calls to ccpp_field_add
          end associate associate_block
#endif

        end do
      end do

   else if (trim(step)=="physics_init") then

      ! Since the physics init steps are independent of the blocking structure,
      ! we can use cdata_domain here. Since we don't use threading on the outside,
      ! we can allow threading inside the time_vary routines.
      GFS_control%nthreads = nthrds

#ifdef STATIC
      call ccpp_physics_init(cdata_domain, suite_name=trim(ccpp_suite), ierr=ierr)
#else
      call ccpp_physics_init(cdata_domain, ierr=ierr)
#endif
      if (ierr/=0) then
        write(0,'(a)') "An error occurred in ccpp_physics_init"
        write(0,'(a)') trim(cdata_domain%errmsg)
        return
      end if

   else if (trim(step)=="time_vary") then

      ! Since the time_vary steps only use data structures for all blocks (except the
      ! CCPP-internal variables ccpp_error_flag and ccpp_error_message, which are defined
      ! for all cdata structures independently), we can use cdata_domain here.
      ! Since we don't use threading on the outside, we can allow threading
      ! inside the time_vary routines.
      GFS_control%nthreads = nthrds

#ifdef STATIC
      call ccpp_physics_run(cdata_domain, suite_name=trim(ccpp_suite), group_name="time_vary", ierr=ierr)
#else
      call ccpp_physics_run(cdata_domain, group_name="time_vary", ierr=ierr)
#endif
      if (ierr/=0) then
        write(0,'(a)') "An error occurred in ccpp_physics_run for group time_vary"
        write(0,'(a)') trim(cdata_domain%errmsg)
        return
      end if

   ! Radiation and stochastic physics
   else if (trim(step)=="radiation" .or. trim(step)=="physics" .or. trim(step)=="stochastics") then

      ! Set number of threads available to physics schemes to one,
      ! because threads are used on the outside for blocking
      GFS_control%nthreads = 1

!$OMP parallel num_threads (nthrds)      &
!$OMP          default (shared)          &
!$OMP          private (nb,nt,ntX,ierr2) &
!$OMP          reduction (+:ierr)
#ifdef OPENMP
      nt = omp_get_thread_num()+1
#else
      nt = 1
#endif
!$OMP do schedule (dynamic,1)
      do nb = 1,nblks
        ! For non-uniform blocks, the last block has a different (shorter)
        ! length than the other blocks; use special CCPP_Interstitial(nthrdsX)
        if (non_uniform_blocks .and. nb==nblks) then
            ntX = nthrdsX
        else
            ntX = nt
        end if
        !--- Call CCPP radiation/physics/stochastics group
#ifdef STATIC
        call ccpp_physics_run(cdata_block(nb,ntX), suite_name=trim(ccpp_suite), group_name=trim(step), ierr=ierr2)
#else
        call ccpp_physics_run(cdata_block(nb,ntX), group_name=trim(step), ierr=ierr2)
#endif
        if (ierr2/=0) then
           write(0,'(2a,3(a,i4),a)') "An error occurred in ccpp_physics_run for group ", trim(step), &
                                     ", block ", nb, " and thread ", nt, " (ntX=", ntX, ")"
           write(0,'(a)') trim(cdata_block(nb,nt)%errmsg)
           ierr = ierr + ierr2
        end if
      end do
!$OMP end do

!$OMP end parallel
      if (ierr/=0) return

   ! Finalize
   else if (trim(step)=="finalize") then

      ! Loop over blocks, don't use threading on the outside but allowing threading
      ! inside the finalization, similar to what is done for the initialization
      GFS_control%nthreads = nthrds

      ! Fast physics are finalized in atmosphere_end, loop over
      ! all blocks and threads to finalize all other physics
      do nt=1,nthrdsX
        do nb=1,nblks
          !--- Finalize CCPP physics
#ifdef STATIC
          call ccpp_physics_finalize(cdata_block(nb,nt), suite_name=trim(ccpp_suite), ierr=ierr)
#else
          call ccpp_physics_finalize(cdata_block(nb,nt), ierr=ierr)
#endif
          if (ierr/=0) then
            write(0,'(a,i4,a,i4)') "An error occurred in ccpp_physics_finalize for block ", nb, " and thread ", nt
            write(0,'(a)') trim(cdata_block(nb,nt)%errmsg)
            return
          end if
          !--- Finalize CCPP framework for blocks/threads
          call ccpp_finalize(cdata_block(nb,nt), ierr=ierr)
          if (ierr/=0) then
            write(0,'(a,i4,a,i4)') "An error occurred in ccpp_finalize for block ", nb, " and thread ", nt
            return
          end if
        end do
      end do

      ! Deallocate cdata structure for blocks and threads
      deallocate(cdata_block)

      !--- Finalize CCPP framework for domain
      call ccpp_finalize(cdata_domain, ierr=ierr)
      if (ierr/=0) then
         write(0,'(a)') "An error occurred in ccpp_finalize"
         return
      end if

      !--- Finalize CCPP framework for fast physics last (all other frameworks point to cdata_tile's suite)
      if (ccpp_initialized(cdata_tile)) then
         call ccpp_finalize(cdata_tile, ierr)
         if (ierr/=0) then
            write(0,'(a)') "An error occurred in ccpp_finalize"
            return
         end if
      end if

    else

      write(0,'(2a)') 'Error, undefined CCPP step ', trim(step)
      ierr = 1
      return

    end if

  end subroutine CCPP_step

end module CCPP_driver
