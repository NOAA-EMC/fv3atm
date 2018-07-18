module IPD_CCPP_driver

  use IPD_typedefs,       only: IPD_control_type, IPD_fastphys_type

  use ccpp_api,           only: ccpp_t,                              &
                                ccpp_init,                           &
                                ccpp_finalize,                       &
                                ccpp_physics_init,                   &
                                ccpp_physics_run,                    &
                                ccpp_physics_finalize,               &
                                ccpp_field_add
  
! Begin include auto-generated list of modules for ccpp
#include "ccpp_modules.inc"
! End include auto-generated list of modules for ccpp

  use iso_c_binding,      only: c_loc

  implicit none

!------------------------------------------------------!
!  CCPP container                                      !
!------------------------------------------------------!
  type(ccpp_t), save, target :: cdata

!----------------
! Public Entities
!----------------
! functions
  public IPD_CCPP_step

  CONTAINS
!*******************************************************************************************

!-------------------------------
!  IPD step generalized for CCPP
!-------------------------------
  subroutine IPD_CCPP_step (step, IPD_Control, IPD_Fastphys, ccpp_suite, ierr)

#ifdef OPENMP
    use omp_lib
#endif

    implicit none

    character(len=*),                intent(in)              :: step
    type(IPD_control_type),  target, intent(inout), optional :: IPD_Control
    type(IPD_fastphys_type), target, intent(inout), optional :: IPD_Fastphys
    character(len=*),                intent(in),    optional :: ccpp_suite
    integer,                         intent(out)             :: ierr

    ierr = 0

    if (trim(step)=="init") then

      if (.not.present(IPD_Control)) then
          write(0,*) 'Optional argument IPD_Control required for IPD-CCPP init step'
          ierr = 1
          return
        end if

      if (.not.present(IPD_Fastphys)) then
          write(0,*) 'Optional argument IPD_Fastphys required for IPD-CCPP init step'
          ierr = 1
          return
        end if

      if (.not.present(ccpp_suite)) then
        write(0,*) 'Optional argument ccpp_suite required for IPD-CCPP init step'
        ierr = 1
        return
      end if

      !--- Initialize CCPP framework
      call ccpp_init(trim(ccpp_suite), cdata, ierr)
      if (ierr/=0) then
        write(0,*) 'An error occurred in ccpp_init'
        return
      end if

! Begin include auto-generated list of calls to ccpp_field_add
#include "ccpp_fields.inc"
! End include auto-generated list of calls to ccpp_field_add

      !--- Initialize CCPP physics
      call ccpp_physics_init(cdata, ierr)
      if (ierr/=0) then
        write(0,*) 'An error occurred in ccpp_physics_init'
        return
      end if

    else if (trim(step)=="fast_physics") then

      call ccpp_physics_run(cdata, group_name='fast_physics', ierr=ierr)
      if (ierr/=0) then
         write(0,'(a)') "An error occurred in ccpp_physics_run for group fast_physics"
         return
      end if

    ! Finalize
    else if (trim(step)=="finalize") then

      !--- Finalize CCPP physics
      call ccpp_physics_finalize(cdata, ierr)
      if (ierr/=0) then
         write(0,'(a)') "An error occurred in ccpp_physics_finalize"
         return
      end if
      !--- Finalize CCPP framework
      call ccpp_finalize(cdata, ierr)
      if (ierr/=0) then
         write(0,'(a)') "An error occurred in ccpp_finalize"
         return
      end if

    else

      write(0,'(2a)') 'Error, undefined IPD step ', trim(step)
      ierr = 1
      return

    end if

  end subroutine IPD_CCPP_step

end module IPD_CCPP_driver
