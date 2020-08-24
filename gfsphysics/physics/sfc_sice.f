!>  \file sfc_sice.f
!!  This file contains the GFS thermodynamics surface ice model.

!> \defgroup GFS_Ice GFS Thermodynamics Surface Ice
!! @{
!!  \brief Brief description of the parameterization
!!  \section diagram Calling Hierarchy Diagram
!!  \section intraphysics Intraphysics Communication

!> \brief Brief description of the subroutine
!!
!! \section arg_table_sice_run Arguments
!! | local var name | longname                                              | description                        | units   | rank | type    |    kind   | intent | optional |
!! |----------------|-------------------------------------------------------|------------------------------------|---------|------|---------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                | horizontal loop extent, start at 1 | index   |    0 | integer |           | in     | F        |
!!
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!!  @{
      module module_sfc_sice
      use machine , only : kind_phys
      use funcphys, only : fpvs
      use physcons, only :    sbc => con_sbc,  hvap => con_hvap,        &
     &                      tgice => con_tice,   cp => con_cp,          &
     &                        eps => con_eps, epsm1 => con_epsm1,       &
     &                     rvrdm1 => con_fvirt, t0c => con_t0c,         &
     &                         rd => con_rd
      implicit none
      contains
!-----------------------------------
      subroutine sfc_sice                                               &
!...................................
!  ---  inputs:
     &     ( im, km, ps, t1, q1, delt,                                  &
!    &     ( im, km, ps, u1, v1, t1, q1, delt,                          &
     &       sfcemis, dlwflx, sfcnsw, sfcdsw, srflag,                   &
     &       cm, ch, prsl1, prslki, islimsk, wind,                      &
     &       flag_iter, lprnt, ipr, cimin,                              &
!  ---  input/outputs:
     &       hice, fice, tice, weasd, tskin, tprcp, stc, ep,            &
!  ---  outputs:
     &       snwdph, qsurf, snowmt, gflux, cmm, chh, evap, hflx         &
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call sfc_sice                                                      !
!       inputs:                                                         !
!          ( im, km, ps, t1, q1, delt,                                  !
!!         ( im, km, ps, u1, v1, t1, q1, delt,                          !
!            sfcemis, dlwflx, sfcnsw, sfcdsw, srflag,                   !
!            cm, ch, prsl1, prslki, islimsk, wind,                      !
!            flag_iter,                                                 !
!       input/outputs:                                                  !
!            hice, fice, tice, weasd, tskin, tprcp, stc, ep,            !
!       outputs:                                                        !
!            snwdph, qsurf, snowmt, gflux, cmm, chh, evap, hflx )       !
!                                                                       !
!  subprogram called:  ice3lay.                                         !
!                                                                       !
!  program history log:                                                 !
!         2005  --  xingren wu created  from original progtm and added  !
!                     two-layer ice model                               !
!         200x  -- sarah lu    added flag_iter                          !
!    oct  2006  -- h. wei      added cmm and chh to output              !
!         2007  -- x. wu modified for mom4 coupling (i.e. cpldice)      !
!                                    (not used anymore)                 !
!         2007  -- s. moorthi micellaneous changes                      !
!    may  2009  -- y.-t. hou   modified to include surface emissivity   !
!                     effect on lw radiation. replaced the confusing    !
!                     slrad with sfc net sw sfcnsw (dn-up). reformatted !
!                     the code and add program documentation block.     !
!    sep  2009 -- s. moorthi removed rcl, changed pressure units and    !
!                     further optimized                                 !
!    jan  2015 -- x. wu change "cimin = 0.15" for both                  !
!                     uncoupled and coupled case                        !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im, km   - integer, horiz dimension and num of soil layers   1    !
!     ps       - real, surface pressure                            im   !
!     t1       - real, surface layer mean temperature ( k )        im   !
!     q1       - real, surface layer mean specific humidity        im   !
!     delt     - real, time interval (second)                      1    !
!     sfcemis  - real, sfc lw emissivity ( fraction )              im   !
!     dlwflx   - real, total sky sfc downward lw flux ( w/m**2 )   im   !
!     sfcnsw   - real, total sky sfc netsw flx into ground(w/m**2) im   !
!     sfcdsw   - real, total sky sfc downward sw flux ( w/m**2 )   im   !
!     srflag   - real, snow/rain fraction for precipitation        im   !
!     cm       - real, surface exchange coeff for momentum (m/s)   im   !
!     ch       - real, surface exchange coeff heat & moisture(m/s) im   !
!     prsl1    - real, surface layer mean pressure                 im   !
!     prslki   - real,                                             im   !
!     islimsk  - integer, sea/land/ice mask (=0/1/2)               im   !
!     wind     - real,                                             im   !
!     flag_iter- logical,                                          im   !
!                                                                       !
!  input/outputs:                                                       !
!     hice     - real, sea-ice thickness                           im   !
!     fice     - real, sea-ice concentration                       im   !
!     tice     - real, sea-ice surface temperature                 im   !
!     weasd    - real, water equivalent accumulated snow depth (mm)im   !
!     tskin    - real, ground surface skin temperature ( k )       im   !
!     tprcp    - real, total precipitation                         im   !
!     stc      - real, soil temp (k)                              im,km !
!     ep       - real, potential evaporation                       im   !
!                                                                       !
!  outputs:                                                             !
!     snwdph   - real, water equivalent snow depth (mm)            im   !
!     qsurf    - real, specific humidity at sfc                    im   !
!     snowmt   - real, snow melt (m)                               im   !
!     gflux    - real, soil heat flux (w/m**2)                     im   !
!     cmm      - real, surface exchange coeff for momentum(m/s)    im   !
!     chh      - real, surface exchange coeff heat&moisture (m/s)  im   !
!     evap     - real, evaperation from latent heat flux           im   !
!     hflx     - real, sensible heat flux                          im   !
!                                                                       !
! ===================================================================== !
!
!
!  ---  constant parameters:
      integer,              parameter :: kmi   = 2                  ! 2-layer of ice
      real(kind=kind_phys), parameter :: zero  = 0.0_kind_phys
      real(kind=kind_phys), parameter :: one   = 1.0_kind_phys
      real(kind=kind_phys), parameter :: cpinv = one/cp
      real(kind=kind_phys), parameter :: hvapi = one/hvap
      real(kind=kind_phys), parameter :: elocp = hvap/cp
      real(kind=kind_phys), parameter :: himax = 8.0_kind_phys    ! maximum ice thickness allowed
      real(kind=kind_phys), parameter :: himin = 0.1_kind_phys    ! minimum ice thickness required
      real(kind=kind_phys), parameter :: hsmax = 2.0_kind_phys    ! maximum snow depth allowed
      real(kind=kind_phys), parameter :: timin = 173.0_kind_phys  ! minimum temperature allowed for snow/ice
      real(kind=kind_phys), parameter :: albfw = 0.06_kind_phys   ! albedo for lead
      real(kind=kind_phys), parameter :: dsi   = one/0.33_kind_phys
      real(kind=kind_phys), parameter :: qmin  = 1.0e-8_kind_phys

!  ---  inputs:
      integer, intent(in) :: im, km, ipr
      logical, intent(in) :: lprnt

      real (kind=kind_phys), dimension(im), intent(in) :: ps,           &
     &       t1, q1, sfcemis, dlwflx, sfcnsw, sfcdsw, srflag, cm, ch,   &
     &       prsl1, prslki, wind

      integer, dimension(im), intent(in) :: islimsk
      real (kind=kind_phys),  intent(in) :: delt, cimin

      logical, intent(in) :: flag_iter(im)

!  ---  input/outputs:
      real (kind=kind_phys), dimension(im), intent(inout) :: hice,      &
     &       fice, tice, weasd, tskin, tprcp, ep

      real (kind=kind_phys), dimension(im,km), intent(inout) :: stc

!  ---  outputs:
      real (kind=kind_phys), dimension(im), intent(inout) :: snwdph,      &
     &       qsurf, snowmt, gflux, cmm, chh, evap, hflx

!  ---  locals:
      real (kind=kind_phys), dimension(im) :: ffw, evapi, evapw,        &
     &       sneti, snetw, hfd, hfi,                                    &
!    &       hflxi, hflxw, sneti, snetw, qssi, qssw, hfd, hfi, hfw,     &
     &       focn, snof,                                   rch, rho,    &
     &       snowd, theta1

      real (kind=kind_phys) :: t12, t14, tem, stsice(im,kmi)
     &,                        hflxi, hflxw, q0, qs1, qssi, qssw


      integer :: i, k

      logical :: flag(im)
!
!===> ...  begin here
!
!  --- ...  set flag for sea-ice

      do i = 1, im
        flag(i) = (islimsk(i) == 2) .and. flag_iter(i)
        if (flag_iter(i) .and. islimsk(i) < 2) then
          hice(i) = zero
          fice(i) = zero
        endif
      enddo
!
      do i = 1, im
        if (flag(i)) then
          if (srflag(i) > zero) then
            ep(i)    = ep(i)*(one-srflag(i))
            weasd(i) = weasd(i) + 1.0d3*tprcp(i)*srflag(i)
            tprcp(i) = tprcp(i)*(one-srflag(i))
          endif
        endif
      enddo
!  --- ...  update sea ice temperature

      do k = 1, kmi
        do i = 1, im
          if (flag(i)) then
            stsice(i,k) = stc(i,k)
          endif
        enddo
      enddo

!  --- ...  initialize variables. all units are supposedly m.k.s. unless specifie
!           psurf is in pascals, wind is wind speed, theta1 is adiabatic surface
!           temp from level 1, rho is density, qs1 is sat. hum. at level1 and qss
!           is sat. hum. at surface
!           convert slrad to the civilized unit from langley minute-1 k-4

      do i = 1, im
        if (flag(i)) then
!         psurf(i) = 1000.0 * ps(i)
!         ps1(i)   = 1000.0 * prsl1(i)

!         dlwflx has been given a negative sign for downward longwave
!         sfcnsw is the net shortwave flux (direction: dn-up)

          q0        = max(q1(i), qmin)
!         tsurf(i)  = tskin(i)
          theta1(i) = t1(i) * prslki(i)
          rho(i)    = prsl1(i) / (rd*t1(i)*(one+rvrdm1*q0))
          qs1       = fpvs(t1(i))
          qs1       = max(eps*qs1 / (prsl1(i) + epsm1*qs1), qmin)
          q0        = min(qs1, q0)

          if (fice(i) < cimin) then
            print *,'warning: ice fraction is low:', fice(i)
            fice(i) = cimin
            tice(i) = tgice
            tskin(i)= tgice
            print *,'fix ice fraction: reset it to:', fice(i)
          endif
          ffw(i)    = one - fice(i)

          qssi = fpvs(tice(i))
          qssi = eps*qssi / (ps(i) + epsm1*qssi)
          qssw = fpvs(tgice)
          qssw = eps*qssw / (ps(i) + epsm1*qssw)

!  --- ...  snow depth in water equivalent is converted from mm to m unit

          snowd(i) = weasd(i) * 0.001_kind_phys
!         flagsnw(i) = .false.

!  --- ...  when snow depth is less than 1 mm, a patchy snow is assumed and
!           soil is allowed to interact with the atmosphere.
!           we should eventually move to a linear combination of soil and
!           snow under the condition of patchy snow.

!  --- ...  rcp = rho cp ch v

          cmm(i) = cm(i)  * wind(i)
          chh(i) = rho(i) * ch(i) * wind(i)
          rch(i) = chh(i) * cp

!  --- ...  sensible and latent heat flux over open water & sea ice

          evapi(i) = elocp * rch(i) * (qssi - q0)
          evapw(i) = elocp * rch(i) * (qssw - q0)
!         evap(i)  = fice(i)*evapi(i) + ffw(i)*evapw(i)

          snetw(i) = sfcdsw(i) * (one - albfw)
          snetw(i) = min(3.0_kind_phys*sfcnsw(i)                        &
     &             / (one+2.0_kind_phys*ffw(i)), snetw(i))
          sneti(i) = (sfcnsw(i) - ffw(i)*snetw(i)) / fice(i)

          t12 = tice(i) * tice(i)
          t14 = t12 * t12

!  --- ...  hfi = net non-solar and upir heat flux @ ice surface

          hfi(i) = -dlwflx(i) + sfcemis(i)*sbc*t14 + evapi(i)           &
     &           + rch(i)*(tice(i) - theta1(i))
          hfd(i) = 4.0_kind_phys*sfcemis(i)*sbc*tice(i)*t12             &
     &           + (one + elocp*eps*hvap*qs1/(rd*t12)) * rch(i)


          t12 = tgice * tgice
          t14 = t12 * t12

!  --- ...  hfw = net heat flux @ water surface (within ice)

!         hfw(i) = -dlwflx(i) + sfcemis(i)*sbc*t14 + evapw(i)           &
!    &           + rch(i)*(tgice - theta1(i)) - snetw(i)

          focn(i) = 2.0_kind_phys   ! heat flux from ocean - should be from ocn model
          snof(i) = zero    ! snowfall rate - snow accumulates in gbphys

          hice(i) = max( min( hice(i), himax ), himin )
          snowd(i) = min( snowd(i), hsmax )

          if (snowd(i) > (2.0_kind_phys*hice(i))) then
            print *, 'warning: too much snow :',snowd(i)
            snowd(i) = hice(i) + hice(i)
            print *,'fix: decrease snow depth to:',snowd(i)
          endif
        endif
      enddo

      call ice3lay
!  ---  inputs:                                                         !
     &     ( im, kmi, fice, flag, hfi, hfd, sneti, focn, delt,          !
     &       lprnt, ipr,
!  ---  outputs:                                                        !
     &       snowd, hice, stsice, tice, snof, snowmt, gflux )           !

      do i = 1, im
        if (flag(i)) then
          if (tice(i) < timin) then
            print *,'warning: snow/ice temperature is too low:',tice(i)
     &,             ' i=',i
            tice(i) = timin
            print *,'fix snow/ice temperature: reset it to:',tice(i)
          endif

          if (stsice(i,1) < timin) then
            print *,'warning: layer 1 ice temp is too low:',stsice(i,1)
     &,             ' i=',i
            stsice(i,1) = timin
            print *,'fix layer 1 ice temp: reset it to:',stsice(i,1)
          endif

          if (stsice(i,2) < timin) then
            print *,'warning: layer 2 ice temp is too low:',stsice(i,2)
            stsice(i,2) = timin
            print *,'fix layer 2 ice temp: reset it to:',stsice(i,2)
          endif

          tskin(i) = tice(i)*fice(i) + tgice*ffw(i)
        endif
      enddo

      do k = 1, kmi
        do i = 1, im
          if (flag(i)) then
            stc(i,k) = min(stsice(i,k), t0c)
          endif
        enddo
      enddo

      do i = 1, im
        if (flag(i)) then
!  --- ...  calculate sensible heat flux (& evap over sea ice)

          hflxi    = rch(i) * (tice(i) - theta1(i))
          hflxw    = rch(i) * (tgice - theta1(i))
          hflx(i)  = fice(i)*hflxi    + ffw(i)*hflxw
          evap(i)  = fice(i)*evapi(i) + ffw(i)*evapw(i)
!
!  --- ...  the rest of the output

          qsurf(i) = q1(i) + evap(i) / (elocp*rch(i))

!  --- ...  convert snow depth back to mm of water equivalent

          weasd(i)  = snowd(i) * 1000.0_kind_phys
          snwdph(i) = weasd(i) * dsi             ! snow depth in mm

          tem     = one / rho(i)
          hflx(i) = hflx(i) * tem * cpinv
          evap(i) = evap(i) * tem * hvapi
        endif
      enddo
!
      return
      end subroutine sfc_sice


!-----------------------------------
!> \brief Brief description of the subroutine
!!
!-----------------------------------
      subroutine ice3lay
!...................................
!  ---  inputs:
     &     ( im, kmi, fice, flag, hfi, hfd, sneti, focn, delt,          &
     &       lprnt, ipr,
!  ---  input/outputs:
     &       snowd, hice, stsice, tice, snof,                           &
!  ---  outputs:
     &       snowmt, gflux                                              &
     &     )

!**************************************************************************
!                                                                         *
!            three-layer sea ice vertical thermodynamics                  *
!                                                                         *
! based on:  m. winton, "a reformulated three-layer sea ice model",       *
! journal of atmospheric and oceanic technology, 2000                     *
!                                                                         *
!                                                                         *
!        -> +---------+ <- tice - diagnostic surface temperature ( <= 0c )*
!       /   |         |                                                   *
!   snowd   |  snow   | <- 0-heat capacity snow layer                     *
!       \   |         |                                                   *
!        => +---------+                                                   *
!       /   |         |                                                   *
!      /    |         | <- t1 - upper 1/2 ice temperature; this layer has *
!     /     |         |         a variable (t/s dependent) heat capacity  *
!   hice    |...ice...|                                                   *
!     \     |         |                                                   *
!      \    |         | <- t2 - lower 1/2 ice temp. (fixed heat capacity) *
!       \   |         |                                                   *
!        -> +---------+ <- base of ice fixed at seawater freezing temp.   *
!                                                                         *
!  =====================  defination of variables  =====================  !
!                                                                         !
!  inputs:                                                         size   !
!     im, kmi  - integer, horiz dimension and num of ice layers      1    !
!     fice     - real, sea-ice concentration                         im   !
!     flag     - logical, ice mask flag                              1    !
!     hfi      - real, net non-solar and heat flux @ surface(w/m^2)  im   !
!     hfd      - real, heat flux derivatice @ sfc (w/m^2/deg-c)      im   !
!     sneti    - real, net solar incoming at top  (w/m^2)            im   !
!     focn     - real, heat flux from ocean    (w/m^2)               im   !
!     delt     - real, timestep                (sec)                 1    !
!                                                                         !
!  input/outputs:                                                         !
!     snowd    - real, surface pressure                              im   !
!     hice     - real, sea-ice thickness                             im   !
!     stsice   - real, temp @ midpt of ice levels  (deg c)          im,kmi!     
!     tice     - real, surface temperature     (deg c)               im   !
!     snof     - real, snowfall rate           (m/sec)               im   !
!                                                                         !
!  outputs:                                                               !
!     snowmt   - real, snow melt during delt   (m)                   im   !
!     gflux    - real, conductive heat flux    (w/m^2)               im   !
!                                                                         !
!  locals:                                                                !
!     hdi      - real, ice-water interface     (m)                        !
!     hsni     - real, snow-ice                (m)                        !
!                                                                         !
! ======================================================================= !
!

!  ---  constant parameters: (properties of ice, snow, and seawater)
      real (kind=kind_phys), parameter :: ds   = 330.0_kind_phys  ! snow (ov sea ice) density (kg/m^3)
      real (kind=kind_phys), parameter :: dw   =1000.0_kind_phys  ! fresh water density  (kg/m^3)
      real (kind=kind_phys), parameter :: dsdw = ds/dw
      real (kind=kind_phys), parameter :: dwds = dw/ds
      real (kind=kind_phys), parameter :: ks   = 0.31_kind_phys   ! conductivity of snow   (w/mk)
      real (kind=kind_phys), parameter :: i0   = 0.3_kind_phys    ! ice surface penetrating solar fraction
      real (kind=kind_phys), parameter :: ki   = 2.03_kind_phys   ! conductivity of ice  (w/mk)
      real (kind=kind_phys), parameter :: di   = 917.0_kind_phys  ! density of ice   (kg/m^3)
      real (kind=kind_phys), parameter :: didw = di/dw
      real (kind=kind_phys), parameter :: dsdi = ds/di
      real (kind=kind_phys), parameter :: ci   = 2054.0_kind_phys ! heat capacity of fresh ice (j/kg/k)
      real (kind=kind_phys), parameter :: li   = 3.34e5_kind_phys ! latent heat of fusion (j/kg-ice)
      real (kind=kind_phys), parameter :: si   = 1.0_kind_phys    ! salinity of sea ice
      real (kind=kind_phys), parameter :: mu   = 0.054_kind_phys  ! relates freezing temp to salinity
      real (kind=kind_phys), parameter :: tfi  = -mu*si           ! sea ice freezing temp = -mu*salinity
      real (kind=kind_phys), parameter :: tfw  = -1.8_kind_phys   ! tfw - seawater freezing temp (c)
      real (kind=kind_phys), parameter :: tfi0 = tfi-0.0001_kind_phys
      real (kind=kind_phys), parameter :: dici = di*ci
      real (kind=kind_phys), parameter :: dili = di*li
      real (kind=kind_phys), parameter :: dsli = ds*li
      real (kind=kind_phys), parameter :: ki4  = ki*4.0_kind_phys

      real (kind=kind_phys), parameter :: zero = 0.0_kind_phys
      real (kind=kind_phys), parameter :: half = 0.5_kind_phys
      real (kind=kind_phys), parameter :: one  = 1.0_kind_phys
      real (kind=kind_phys), parameter :: four = 4.0_kind_phys

!  ---  inputs:
      integer, intent(in) :: im, kmi, ipr
      logical             :: lprnt

      real (kind=kind_phys), dimension(im), intent(in) :: fice, hfi,    &
     &       hfd, sneti, focn

      real (kind=kind_phys),  intent(in) :: delt

      logical, dimension(im), intent(in) :: flag

!  ---  input/outputs:
      real (kind=kind_phys), dimension(im), intent(inout) :: snowd,     &
     &       hice, tice, snof

      real (kind=kind_phys), dimension(im,kmi), intent(inout) :: stsice

!  ---  outputs:
      real (kind=kind_phys), dimension(im), intent(out) :: snowmt,      &
     &       gflux

!  ---  locals:

      real (kind=kind_phys) :: dt2, dt4, dt6, h1, h2, dh, wrk, wrk1,    &
     &                         dt2i, hdi, hsni, ai, bi, a1, b1, a10, b10&
     &,                        c1, ip, k12, k32, tsf, f1, tmelt, bmelt

      integer :: i
!
!===> ...  begin here
!
      dt2  =  delt + delt
      dt4  =  dt2  + dt2
      dt6  =  dt2  + dt4
      dt2i = one / dt2

      do i = 1, im
        if (flag(i)) then
          snowd(i) = snowd(i) * dwds
          hdi      = (dsdw*snowd(i) + didw*hice(i))

          if (hice(i) < hdi) then
            snowd(i) = snowd(i) + hice(i) - hdi
            hsni     = (hdi - hice(i)) * dsdi
            hice (i) = hice(i) + hsni
          endif

          snof(i)     = snof(i) * dwds
          tice(i)     = tice(i) - t0c
          stsice(i,1) = min(stsice(i,1)-t0c, tfi0)     ! degc
          stsice(i,2) = min(stsice(i,2)-t0c, tfi0)     ! degc

          ip = i0 * sneti(i)         ! ip +v (in winton ip=-i0*sneti as sol -v)
          if (snowd(i) > zero) then
            tsf = zero
            ip  = zero
          else
            tsf = tfi
            ip  = i0 * sneti(i)      ! ip +v here (in winton ip=-i0*sneti)
          endif
          tice(i) = min(tice(i), tsf)

!  --- ...  compute ice temperature

          bi   = hfd(i)
          ai   = hfi(i) - sneti(i) + ip - tice(i)*bi  ! +v sol input here
          k12  = ki4*ks / (ks*hice(i) + ki4*snowd(i))
          k32  = (ki+ki) / hice(i)

          wrk    = one / (dt6*k32 + dici*hice(i))
          a10    = dici*hice(i)*dt2i + k32*(dt4*k32 + dici*hice(i))*wrk
          b10    = -di*hice(i) * (ci*stsice(i,1) + li*tfi/stsice(i,1))  &
     &           * dt2i - ip                                            &
     &           - k32*(dt4*k32*tfw + dici*hice(i)*stsice(i,2)) * wrk

          wrk1  = k12 / (k12 + bi)
          a1    = a10 + bi * wrk1
          b1    = b10 + ai * wrk1
          c1    = dili * tfi * dt2i * hice(i)

          stsice(i,1) = -(sqrt(b1*b1 - four*a1*c1) + b1)/(a1+a1)
          tice(i) = (k12*stsice(i,1) - ai) / (k12 + bi)
  
          if (tice(i) > tsf) then
            a1 = a10 + k12
            b1 = b10 - k12*tsf
            stsice(i,1) = -(sqrt(b1*b1 - four*a1*c1) + b1)/(a1+a1)
            tice(i) = tsf
            tmelt   = (k12*(stsice(i,1)-tsf) - (ai+bi*tsf)) * delt
          else
            tmelt    = zero
            snowd(i) = snowd(i) + snof(i)*delt
          endif

          stsice(i,2) = (dt2*k32*(stsice(i,1) + tfw + tfw)              &
     &                +  dici*hice(i)*stsice(i,2)) * wrk

          bmelt = (focn(i) + ki4*(stsice(i,2) - tfw)/hice(i)) * delt

!  --- ...  resize the ice ...

          h1 = half * hice(i)
          h2 = half * hice(i)

!  --- ...  top ...

          if (tmelt <= snowd(i)*dsli) then
            snowmt(i) = tmelt / dsli
            snowd (i) = snowd(i) - snowmt(i)
          else
            snowmt(i) = snowd(i)
            h1 = h1 - (tmelt - snowd(i)*dsli)                           &
     &         / (di * (ci - li/stsice(i,1)) * (tfi - stsice(i,1)))
            snowd(i) = zero
          endif

!  --- ...  and bottom

          if (bmelt < zero) then
            dh = -bmelt / (dili + dici*(tfi - tfw))
            stsice(i,2) = (h2*stsice(i,2) + dh*tfw) / (h2 + dh)
            h2 = h2 + dh
          else
            h2 = h2 - bmelt / (dili + dici*(tfi - stsice(i,2)))
          endif

!  --- ...  if ice remains, even up 2 layers, else, pass negative energy back in snow

          hice(i) = h1 + h2

          if (hice(i) > zero) then
            if (h1 > half*hice(i)) then
              f1 = one - (h2+h2) / hice(i)
              stsice(i,2) = f1 * (stsice(i,1) + li*tfi/(ci*stsice(i,1)))&
     &                    + (one - f1)*stsice(i,2)

              if (stsice(i,2) > tfi) then
                hice(i) = hice(i) - h2*ci*(stsice(i,2) - tfi)/ (li*delt)
                stsice(i,2) = tfi
              endif
            else
              f1 = (h1+h1) / hice(i)
              stsice(i,1) = f1 * (stsice(i,1) + li*tfi/(ci*stsice(i,1)))&
     &                    + (one - f1)*stsice(i,2)
              stsice(i,1) = (stsice(i,1) - sqrt(stsice(i,1)*stsice(i,1) &
     &                    - four*tfi*li/ci)) * half
            endif

            k12      = ki4*ks / (ks*hice(i) + ki4*snowd(i))
            gflux(i) = k12 * (stsice(i,1) - tice(i))
          else
            snowd(i) = snowd(i) + (h1*(ci*(stsice(i,1) - tfi)           &
     &               - li*(one - tfi/stsice(i,1)))                      &
     &               + h2*(ci*(stsice(i,2) - tfi) - li)) / li

            hice(i)     = max(zero, snowd(i)*dsdi)
            snowd(i)    = zero
            stsice(i,1) = tfw
            stsice(i,2) = tfw
            gflux(i)    = zero
          endif   ! end if_hice_block

          gflux(i)    = fice(i) * gflux(i)
          snowmt(i)   = snowmt(i) * dsdw
          snowd(i)    = snowd(i) * dsdw
          tice(i)     = tice(i)     + t0c
          stsice(i,1) = stsice(i,1) + t0c
          stsice(i,2) = stsice(i,2) + t0c
        endif   ! end if_flag_block
      enddo   ! end do_i_loop

      return
!...................................
      end subroutine ice3lay
!-----------------------------------

!-----------------------------------
      end module module_sfc_sice
!> @}
