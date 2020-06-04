!>  \file sfc_cice.f
!!  This file contains the sfc_sice for coupling to CICE
!> \defgroup sfc_sice for coupling to CICE
!! @{
!!  \section diagram Calling Hierarchy Diagram
!!  \section intraphysics Intraphysics Communication
!!
!> \brief Brief description of the subroutine
!!
!! \section arg_table_cice_run Arguments
!! | local var name | longname                                           | description                        | units   | rank | type    |    kind   | intent | optional |
!! |----------------|----------------------------------------------------|------------------------------------|---------|------|---------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                             | horizontal loop extent, start at 1 | index   |    0 | integer |           | in     | F        |
!!
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!!  @{
!
      module module_sfc_cice
      use machine , only : kind_phys
      use physcons, only : hvap   => con_hvap,  cp => con_cp,           &
     &                     rvrdm1 => con_fvirt, rd => con_rd
      implicit none
      contains
!
!-----------------------------------
      subroutine sfc_cice                                               &
!...................................
!  ---  inputs:
     &     ( im, t1, q1, cm, ch, prsl1,                                 &
     &       wind, flag_cice, flag_iter, dqsfc, dtsfc,                  &
     &       dusfc, dvsfc,                                              &
!  ---  outputs:
     &       qsurf, cmm, chh, evap, hflx, stress )

! ===================================================================== !
!  description:                                                         !
!  Sep 2015  --  Xingren Wu created from sfc_sice for coupling to CICE  !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call sfc_cice                                                      !
!       inputs:                                                         !
!          ( im, t1, q1, cm, ch, prsl1,                                 !
!            wind, flag_cice, flag_iter, dqsfc, dtsfc,                  !
!       outputs:                                                        !
!            qsurf, cmm, chh, evap, hflx)                               !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:
!     im, - integer, horiz dimension
!!    u1, v1   - real, u/v component of surface layer wind
!     t1       - real, surface layer mean temperature ( k )
!     q1       - real, surface layer mean specific humidity
!     cm       - real, surface exchange coeff for momentum (m/s)
!     ch       - real, surface exchange coeff heat & moisture(m/s)
!     prsl1    - real, surface layer mean pressure
!     islimsk  - integer, sea/land/ice mask
!     wind     - real, wind speed (m/s)
!     flag_iter- logical
!     dqsfc    - real, latent heat flux
!     dtsfc    - real, sensible heat flux
!     dusfc    - real, zonal momentum stress
!     dvsfc    - real, meridional momentum stress
!     dvsfc    - real, sensible heat flux
!  outputs:
!     qsurf    - real, specific humidity at sfc
!     cmm      - real, ?
!     chh      - real, ?
!     evap     - real, evaperation from latent heat
!     hflx     - real, sensible heat
!     stress   - real, surface stress
!  ====================    end of description    =====================  !
!
!
!  ---  constant parameters:
      real(kind=kind_phys), parameter :: cpinv = 1.0/cp
      real(kind=kind_phys), parameter :: hvapi = 1.0/hvap

!  ---  inputs:
      integer, intent(in) :: im

!     real (kind=kind_phys), dimension(im), intent(in) :: u1, v1,       &
      real (kind=kind_phys), dimension(im), intent(in) ::               &
     &       t1, q1, cm, ch, prsl1, wind, dqsfc, dtsfc, dusfc, dvsfc

      logical,                intent(in) :: flag_cice(im), flag_iter(im)

!  ---  outputs:
      real (kind=kind_phys), dimension(im), intent(out) :: qsurf,       &
     &                                  cmm, chh, evap, hflx, stress

!  ---  locals:

      real (kind=kind_phys) :: rho, tem

      integer :: i
 
      logical :: flag(im)
!
      do i = 1, im
        flag(i) = flag_cice(i) .and. flag_iter(i)
      enddo
!
      do i = 1, im
        if (flag(i)) then

          rho    = prsl1(i)                                             &
     &           / (rd * t1(i) * (1.0 + rvrdm1*max(q1(i), 1.0e-8)))

          cmm(i) = wind(i) * cm(i)
          chh(i) = wind(i) * ch(i) * rho

          qsurf(i)  = q1(i) + dqsfc(i) / (hvap*chh(i))
          tem       = 1.0 / rho
          hflx(i)   = dtsfc(i) * tem * cpinv
          evap(i)   = dqsfc(i) * tem * hvapi
          stress(i) = sqrt(dusfc(i)*dusfc(i) + dvsfc(i)*dvsfc(i)) * tem
        endif
      enddo

      return
!-----------------------------------
      end subroutine sfc_cice
!-----------------------------------

!-----------------------------------
      end module module_sfc_cice
!-----------------------------------
!> @}
