      module module_sfc_diff

      use machine , only : kind_phys
      use physcons, grav => con_g
      real (kind=kind_phys), parameter :: ca=0.4d0  ! ca - von karman constant

      contains
      subroutine sfc_diff(im,ps,t1,q1,z1, wind,                 !intent(in)
     &                    prsl1,prslki,                         !intent(in)
     &                    sigmaf,vegtype,shdmax,ivegsrc,        !intent(in)
     &                    z0pert,ztpert,        ! mg, sfc-perts !intent(in)
     &                    flag_iter,redrag,                     !intent(in)
     &                    u10m,v10m,sfc_z0_type,  !hafs,z0 type !intent(in)
     &                    wet,dry,icy,                          !intent(in)
     &                    tskin, tsurf, snwdph, z0rl, z0rlw, ustar
!
     &,                   cm, ch, rb, stress, fm, fh, fm10, fh2)
!
      use physcons, rvrdm1 => con_fvirt
      implicit none
!
! 1 - land, 2 - ice, 3 - water
! --------  -------- ---------
      integer, intent(in) :: im, ivegsrc
      integer, intent(in) :: sfc_z0_type ! option for calculating surface roughness length over ocean
      integer, dimension(im), intent(in) :: vegtype

      logical, intent(in) :: redrag ! reduced drag coeff. flag for high wind over sea (j.han)
      logical, dimension(im), intent(in) :: flag_iter, wet, dry, icy

      real(kind=kind_phys), dimension(im), intent(in)    :: u10m,v10m
      real(kind=kind_phys), dimension(im), intent(in)    ::
     &                    ps,t1,q1,z1,prsl1,prslki,wind,
     &                    sigmaf,shdmax,
     &                    z0pert,ztpert ! mg, sfc-perts
      real(kind=kind_phys), dimension(im,3), intent(in)    ::
     &                    tskin, tsurf, snwdph

      real(kind=kind_phys), dimension(im), intent(in) :: z0rlw
      real(kind=kind_phys), dimension(im,3), intent(inout) ::
     &                       z0rl, ustar

! 1 - land, 2 - ice, 3 - water
! --------  -------- ---------
      real(kind=kind_phys), dimension(im,3), intent(out)   ::
     &                       cm, ch, rb, stress, fm, fh, fm10, fh2
!
!     locals
!
      integer   i
!
      real(kind=kind_phys) :: rat,   thv1, restar, wind10m,
     &                        czilc, tem1, tem2, virtfac

      real(kind=kind_phys) :: tvs, z0, z0max, ztmax
!
      real(kind=kind_phys), parameter ::
     &              one=1.0d0, zero=0.0d0, half=0.5d0, qmin=1.0d-8
     &,             charnock=.014d0, z0s_max=.317d-2       ! a limiting value at high winds over sea
     &,             zmin=1.0d-6
     &,             vis=1.4d-5, rnu=1.51d-5, visi=one/vis
     &,             log01=log(0.01), log05=log(0.05), log07=log(0.07)

!     parameter (charnock=.014,ca=.4)!c ca is the von karman constant
!     parameter (alpha=5.,a0=-3.975,a1=12.32,b1=-7.755,b2=6.041)
!     parameter (a0p=-7.941,a1p=24.75,b1p=-8.705,b2p=7.899,vis=1.4e-5)

!     real(kind=kind_phys) aa1,bb1,bb2,cc,cc1,cc2,arnu
!     parameter (aa1=-1.076,bb1=.7045,cc1=-.05808)
!     parameter (bb2=-.1954,cc2=.009999)
!     parameter (arnu=.135*rnu)
!
!    z0s_max=.196e-2 for u10_crit=25 m/s
!    z0s_max=.317e-2 for u10_crit=30 m/s
!    z0s_max=.479e-2 for u10_crit=35 m/s
!
! mbek -- toga-coare flux algorithm
!     parameter (rnu=1.51e-5,arnu=0.11*rnu)
!
!  initialize variables. all units are supposedly m.k.s. unless specified
!  ps is in pascals, wind is wind speed,
!  surface roughness length is converted to m from cm
!

!       write(0,*)'in sfc_diff, sfc_z0_type=',sfc_z0_type

      do i=1,im
        if(flag_iter(i)) then
          virtfac = one + rvrdm1 * max(q1(i),qmin)
          thv1    = t1(i) * prslki(i) * virtfac

!  compute stability dependent exchange coefficients
!  this portion of the code is presently suppressed
!
          if (dry(i)) then ! Some land
            tvs   = half * (tsurf(i,1)+tskin(i,1)) * virtfac
            z0max = max(zmin, min(0.01d0 * z0rl(i,1), z1(i)))
!** xubin's new z0  over land
            tem1 = one  - shdmax(i)
            tem2 = tem1 * tem1
            tem1 = one  - tem2

            if( ivegsrc == 1 ) then

              if (vegtype(i) == 10) then
                z0max = exp( tem2*log01 + tem1*log07 )
              elseif (vegtype(i) == 6) then
                z0max = exp( tem2*log01 + tem1*log05 )
              elseif (vegtype(i) == 7) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max = 0.01d0
              elseif (vegtype(i) == 16) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max = 0.01d0
              else
                z0max = exp( tem2*log01 + tem1*log(z0max) )
              endif

            elseif (ivegsrc == 2 ) then

              if (vegtype(i) == 7) then
                z0max = exp( tem2*log01 + tem1*log07 )
              elseif (vegtype(i) == 8) then
                z0max = exp( tem2*log01 + tem1*log05 )
              elseif (vegtype(i) == 9) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max = 0.01d0
              elseif (vegtype(i) == 11) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max = 0.01d0
              else
                z0max = exp( tem2*log01 + tem1*log(z0max) )
              endif

            endif
! mg, sfc-perts: add surface perturbations to z0max over land
            if (z0pert(i) /= zero ) then
              z0max = z0max * (10.0d0**z0pert(i))
            endif

            z0max = max(z0max, zmin)

!           czilc = 10.0 ** (- (0.40/0.07) * z0) ! fei's canopy height dependance of czil
            czilc = 0.8d0

            tem1  = 1.0d0 - sigmaf(i)
            ztmax = z0max*exp( - tem1*tem1
     &              * czilc*ca*sqrt(ustar(i,1)*(0.01/1.5e-05)))
 
! mg, sfc-perts: add surface perturbations to ztmax/z0max ratio over land
            if (ztpert(i) /= zero) then
              ztmax = ztmax * (10.0d0**ztpert(i))
            endif
            ztmax = max(ztmax, zmin)
!
            call stability
!  ---  inputs:
     &       (z1(i), snwdph(i,1), thv1, wind(i), z0max, ztmax, tvs,
!  ---  outputs:
     &        rb(i,1), fm(i,1), fh(i,1), fm10(i,1), fh2(i,1),
     &        cm(i,1), ch(i,1), stress(i,1), ustar(i,1))
          endif ! Dry points

          if (icy(i)) then ! Some ice
            tvs   = half * (tsurf(i,2)+tskin(i,2)) * virtfac
            z0max = max(zmin, min(0.01d0 * z0rl(i,2), z1(i)))
!** xubin's new z0  over land and sea ice
            tem1 = one  - shdmax(i)
            tem2 = tem1 * tem1
            tem1 = one  - tem2

            if( ivegsrc == 1 ) then

              z0max = exp( tem2*log01 + tem1*log(z0max) )
            elseif (ivegsrc == 2 ) then
              z0max = exp( tem2*log01 + tem1*log(z0max) )
            endif

            z0max = max(z0max, zmin)

!           czilc = 10.0 ** (- (0.40/0.07) * z0) ! fei's canopy height
!           dependance of czil

            czilc = 0.8d0

            tem1 = 1.0d0 - sigmaf(i)
            ztmax = z0max*exp( - tem1*tem1
     &              * czilc*ca*sqrt(ustar(i,2)*(0.01/1.5e-05)))
            ztmax = max(ztmax, 1.0e-6)
!
            call stability
!  ---  inputs:
     &       (z1(i), snwdph(i,2), thv1, wind(i), z0max, ztmax, tvs,
!  ---  outputs:
     &        rb(i,2), fm(i,2), fh(i,2), fm10(i,2), fh2(i,2),
     &        cm(i,2), ch(i,2), stress(i,2), ustar(i,2))
          endif ! Icy points

! BWG: Everything from here to end of subroutine was after
!      the stuff now put into "stability"

          if (wet(i)) then ! Some open ocean
            tvs        = half * (tsurf(i,3)+tskin(i,3)) * virtfac
            z0         = 0.01d0 * z0rl(i,3)
            z0max      = max(zmin, min(z0,z1(i)))
            ustar(i,3) = sqrt(grav * z0 / charnock)
            wind10m    = sqrt(u10m(i)*u10m(i)+v10m(i)*v10m(i))

!**  test xubin's new z0

!           ztmax  = z0max

            restar = max(ustar(i,3)*z0max*visi, 0.000001d0)

!           restar = log(restar)
!           restar = min(restar,5.)
!           restar = max(restar,-5.)
!           rat    = aa1 + (bb1 + cc1*restar) * restar
!           rat    = rat    / (1. + (bb2 + cc2*restar) * restar))
!  rat taken from zeng, zhao and dickinson 1997

            rat   = min(7.0d0, 2.67d0 * sqrt(sqrt(restar)) - 2.57d0)
            ztmax = max(z0max * exp(-rat), zmin)
!
            if (sfc_z0_type == 6) then
              call znot_t_v6(wind10m, ztmax)   ! 10-m wind,m/s, ztmax(m)
            else if (sfc_z0_type == 7) then
              call znot_t_v7(wind10m, ztmax)   ! 10-m wind,m/s, ztmax(m)
            else if (sfc_z0_type > 0) then
              write(0,*)'no option for sfc_z0_type=',sfc_z0_type
              stop
            endif
!
            call stability
!  ---  inputs:
     &       (z1(i), snwdph(i,3), thv1, wind(i), z0max, ztmax, tvs,
!  ---  outputs:
     &        rb(i,3), fm(i,3), fh(i,3), fm10(i,3), fh2(i,3),
     &        cm(i,3), ch(i,3), stress(i,3), ustar(i,3))
!
!  update z0 over ocean
!
            if (sfc_z0_type >= 0) then
              if (sfc_z0_type == 0) then
                z0 = (charnock / grav) * ustar(i,3) * ustar(i,3)

! mbek -- toga-coare flux algorithm
!               z0 = (charnock / grav) * ustar(i)*ustar(i) +  arnu/ustar(i)
!  new implementation of z0
!               cc = ustar(i) * z0 / rnu
!               pp = cc / (1. + cc)
!               ff = grav * arnu / (charnock * ustar(i) ** 3)
!               z0 = arnu / (ustar(i) * ff ** pp)

                if (redrag) then
                  z0rl(i,3) = 100.0d0 * max(min(z0, z0s_max), 1.0d-7)
                else
                  z0rl(i,3) = 100.0d0 * max(min(z0, 0.1d0), 1.0d-7)
                endif

              elseif (sfc_z0_type == 6) then   ! wang
                 call znot_m_v6(wind10m, z0)   ! wind, m/s, z0, m
                 z0rl(i,3) = 100.0d0 * z0      ! cm
              elseif (sfc_z0_type == 7) then   ! wang
                 call znot_m_v7(wind10m, z0)   ! wind, m/s, z0, m
                 z0rl(i,3) = 100.0d0 * z0      ! cm
              else
                 z0rl(i,3) = 1.0d-4
              endif

            elseif (z0rlw(i) < 1.0d-7) then
              z0 = (charnock / grav) * ustar(i,3) * ustar(i,3)

              if (redrag) then
                z0rl(i,3) = 100.0d0 * max(min(z0, z0s_max), 1.0d-7)
              else
                z0rl(i,3) = 100.0d0 * max(min(z0, 0.1d0), 1.0d-7)
              endif

            endif

          endif              ! end of if(open ocean)
!
        endif                ! end of if(flagiter) loop
      enddo

      return
      end subroutine sfc_diff


!----------------------------------------
      subroutine stability
!........................................
!  ---  inputs:
     &     ( z1, snwdph, thv1, wind, z0max, ztmax, tvs,                 &
!  ---  outputs:
     &       rb, fm, fh, fm10, fh2, cm, ch, stress, ustar)
!-----

!  ---  inputs:
      real(kind=kind_phys), intent(in) ::                               &
     &       z1, snwdph, thv1, wind, z0max, ztmax, tvs

!  ---  outputs:
      real(kind=kind_phys), intent(out) ::                              &
     &       rb, fm, fh, fm10, fh2, cm, ch, stress, ustar

!  ---  locals:
      real(kind=kind_phys), parameter :: alpha=5.0d0, a0=-3.975d0       &
     &,             a1=12.32d0, alpha4=4.0d0*alpha                      &
     &,             b1=-7.755d0,  b2=6.041d0,  alpha2=alpha+alpha       &
     &,             beta=1.0d0                                          &
     &,             a0p=-7.941d0, a1p=24.75d0, b1p=-8.705d0, b2p=7.899d0&
     &,             ztmin1=-999.0d0, zero=0.0d0, one=1.0d0

      real(kind=kind_phys) aa,     aa0,    bb,     bb0, dtv,   adtv,
     &                     hl1,    hl12,   pm,     ph,  pm10,  ph2,
     &                     z1i,
     &                     fms,    fhs,    hl0,    hl0inf, hlinf,
     &                     hl110,  hlt,    hltinf, olinf,
     &                     tem1,   tem2, ztmax1

          z1i = one / z1

          tem1   = z0max/z1
          if (abs(one-tem1) > 1.0d-6) then
            ztmax1 = - beta*log(tem1)/(alpha2*(one-tem1))
          else
            ztmax1 = 99.0d0
          endif
          if( z0max < 0.05d0 .and. snwdph < 10.0d0 ) ztmax1 = 99.0d0

!  compute stability indices (rb and hlinf)

          dtv     = thv1 - tvs
          adtv    = max(abs(dtv),0.001d0)
          dtv     = sign(1.,dtv) * adtv
          rb      = max(-5000.0d0, (grav+grav) * dtv * z1
     &            / ((thv1 + tvs) * wind * wind))
          tem1    = one / z0max
          tem2    = one / ztmax
          fm      = log((z0max+z1)  * tem1)
          fh      = log((ztmax+z1)  * tem2)
          fm10    = log((z0max+10.0d0) * tem1)
          fh2     = log((ztmax+2.0d0)  * tem2)
          hlinf   = rb * fm * fm / fh
          hlinf   = min(max(hlinf,ztmin1),ztmax1)
!
!  stable case
!
          if (dtv >= zero) then
            hl1 = hlinf
            if(hlinf > 0.25d0) then
              tem1   = hlinf * z1i
              hl0inf = z0max * tem1
              hltinf = ztmax * tem1
              aa     = sqrt(one + alpha4 * hlinf)
              aa0    = sqrt(one + alpha4 * hl0inf)
              bb     = aa
              bb0    = sqrt(one + alpha4 * hltinf)
              pm     = aa0 - aa + log( (aa + one)/(aa0 + one) )
              ph     = bb0 - bb + log( (bb + one)/(bb0 + one) )
              fms    = fm - pm
              fhs    = fh - ph
              hl1    = fms * fms * rb / fhs
              hl1    = min(max(hl1, ztmin1), ztmax1)
            endif
!
!  second iteration
!
            tem1  = hl1 * z1i
            hl0   = z0max * tem1
            hlt   = ztmax * tem1
            aa    = sqrt(one + alpha4 * hl1)
            aa0   = sqrt(one + alpha4 * hl0)
            bb    = aa
            bb0   = sqrt(one + alpha4 * hlt)
            pm    = aa0 - aa + log( (one+aa)/(one+aa0) )
            ph    = bb0 - bb + log( (one+bb)/(one+bb0) )
            hl110 = hl1 * 10.0d0 * z1i
            hl110 = min(max(hl110, ztmin1), ztmax1)
            aa    = sqrt(one + alpha4 * hl110)
            pm10  = aa0 - aa + log( (one+aa)/(one+aa0) )
            hl12  = (hl1+hl1) * z1i
            hl12  = min(max(hl12,ztmin1),ztmax1)
!           aa    = sqrt(one + alpha4 * hl12)
            bb    = sqrt(one + alpha4 * hl12)
            ph2   = bb0 - bb + log( (one+bb)/(one+bb0) )
!
!  unstable case - check for unphysical obukhov length
!
          else                          ! dtv < 0 case
            olinf = z1 / hlinf
            tem1  = 50.0d0 * z0max
            if(abs(olinf) <= tem1) then
              hlinf = -z1 / tem1
              hlinf = min(max(hlinf,ztmin1),ztmax1)
            endif
!
!  get pm and ph
!
            if (hlinf >= -0.5d0) then
              hl1   = hlinf
              pm    = (a0  + a1*hl1)  * hl1  / (one+ (b1+b2*hl1)  *hl1)
              ph    = (a0p + a1p*hl1) * hl1  / (one+ (b1p+b2p*hl1)*hl1)
              hl110 = hl1 * 10.0d0 * z1i
              hl110 = min(max(hl110, ztmin1), ztmax1)
              pm10  = (a0 + a1*hl110) * hl110/(one+(b1+b2*hl110)*hl110)
              hl12  = (hl1+hl1) * z1i
              hl12  = min(max(hl12, ztmin1), ztmax1)
              ph2   = (a0p + a1p*hl12) * hl12/(one+(b1p+b2p*hl12)*hl12)
            else                       ! hlinf < 0.05
              hl1   = -hlinf
              tem1  = one / sqrt(hl1)
              pm    = log(hl1) + 2.0d0 * sqrt(tem1) - .8776d0
              ph    = log(hl1) + 0.5d0 * tem1 + 1.386d0
!             pm    = log(hl1) + 2.0 * hl1 ** (-.25) - .8776
!             ph    = log(hl1) + 0.5 * hl1 ** (-.5) + 1.386
              hl110 = hl1 * 10.0d0 * z1i
              hl110 = min(max(hl110, ztmin1), ztmax1)
              pm10  = log(hl110) + 2.0d0 / sqrt(sqrt(hl110)) - 0.8776d0
!             pm10  = log(hl110) + 2. * hl110 ** (-.25) - .8776
              hl12  = (hl1+hl1) * z1i
              hl12  = min(max(hl12, ztmin1), ztmax1)
              ph2   = log(hl12) + 0.5d0 / sqrt(hl12) + 1.386d0
!             ph2   = log(hl12) + .5 * hl12 ** (-.5) + 1.386
            endif

          endif          ! end of if (dtv >= 0 ) then loop
!
!  finish the exchange coefficient computation to provide fm and fh
!
          fm        = fm - pm
          fh        = fh - ph
          fm10      = fm10 - pm10
          fh2       = fh2 - ph2
          cm        = ca * ca / (fm * fm)
          ch        = ca * ca / (fm * fh)
          tem1      = 0.00001d0 / z1
          cm        = max(cm, tem1)
          ch        = max(ch, tem1)
          stress    = cm * wind * wind
          ustar     = sqrt(stress)

      return
!.................................
      end subroutine stability
!---------------------------------


!! add fitted z0,zt curves for hurricane application (used in HWRF/HMON)
!! Weiguo Wang, 2019-0425

      SUBROUTINE znot_m_v6(uref, znotm)
      use machine , only : kind_phys
      IMPLICIT NONE
! Calculate areodynamical roughness over water with input 10-m wind
! For low-to-moderate winds, try to match the Cd-U10 relationship from COARE V3.5 (Edson et al. 2013)
! For high winds, try to fit available observational data
!
! Bin Liu, NOAA/NCEP/EMC 2017
! 
! uref(m/s)   :   wind speed at 10-m height
! znotm(meter):   areodynamical roughness scale over water
!

      REAL(kind=kind_phys), INTENT(IN) :: uref
      REAL(kind=kind_phys), INTENT(OUT):: znotm
      real(kind=kind_phys), parameter  :: p13 = -1.296521881682694e-02,
     &      p12 =  2.855780863283819e-01, p11 = -1.597898515251717e+00,
     &      p10 = -8.396975715683501e+00,

     &      p25 =  3.790846746036765e-10, p24 =  3.281964357650687e-09,
     &      p23 =  1.962282433562894e-07, p22 = -1.240239171056262e-06,
     &      p21 =  1.739759082358234e-07, p20 =  2.147264020369413e-05,

     &      p35 =  1.840430200185075e-07, p34 = -2.793849676757154e-05,
     &      p33 =  1.735308193700643e-03, p32 = -6.139315534216305e-02,
     &      p31 =  1.255457892775006e+00, p30 = -1.663993561652530e+01,

     &      p40 =  4.579369142033410e-04
  

       if (uref >= 0.0 .and.  uref <= 6.5 ) then
        znotm = exp(p10 + uref * (p11 + uref * (p12 + uref*p13))) 
       elseif (uref > 6.5 .and. uref <= 15.7) then
        znotm = p20 + uref * (p21 + uref * (p22 + uref * (p23
     &              + uref * (p24 + uref * p25))))
       elseif (uref > 15.7 .and. uref <= 53.0) then
        znotm = exp( p30 + uref * (p31 + uref * (p32 + uref * (p33
     &                   + uref * (p34 + uref * p35)))))
       elseif ( uref > 53.0) then
         znotm = p40
       else
        print*, 'Wrong input uref value:',uref
       endif

       END SUBROUTINE znot_m_v6

      SUBROUTINE znot_t_v6(uref, znott)
      use machine , only : kind_phys
      IMPLICIT NONE
! Calculate scalar roughness over water with input 10-m wind
! For low-to-moderate winds, try to match the Ck-U10 relationship from COARE algorithm
! For high winds, try to retain the Ck-U10 relationship of FY2015 HWRF
!
! Bin Liu, NOAA/NCEP/EMC 2017
!
! uref(m/s)   :   wind speed at 10-m height
! znott(meter):   scalar roughness scale over water
!

      REAL(kind=kind_phys), INTENT(IN) :: uref
      REAL(kind=kind_phys), INTENT(OUT):: znott
      real(kind=kind_phys), parameter  :: p00 =  1.100000000000000e-04,
     &      p15 = -9.144581627678278e-10, p14 =  7.020346616456421e-08,
     &      p13 = -2.155602086883837e-06, p12 =  3.333848806567684e-05,
     &      p11 = -2.628501274963990e-04, p10 =  8.634221567969181e-04,

     &      p25 = -8.654513012535990e-12, p24 =  1.232380050058077e-09,
     &      p23 = -6.837922749505057e-08, p22 =  1.871407733439947e-06,
     &      p21 = -2.552246987137160e-05, p20 =  1.428968311457630e-04,

     &      p35 =  3.207515102100162e-12, p34 = -2.945761895342535e-10,
     &      p33 =  8.788972147364181e-09, p32 = -3.814457439412957e-08,
     &      p31 = -2.448983648874671e-06, p30 =  3.436721779020359e-05,

     &      p45 = -3.530687797132211e-11, p44 =  3.939867958963747e-09,
     &      p43 = -1.227668406985956e-08, p42 = -1.367469811838390e-05,
     &      p41 =  5.988240863928883e-04, p40 = -7.746288511324971e-03,

     &      p56 = -1.187982453329086e-13, p55 =  4.801984186231693e-11,
     &      p54 = -8.049200462388188e-09, p53 =  7.169872601310186e-07,
     &      p52 = -3.581694433758150e-05, p51 =  9.503919224192534e-04,
     &      p50 = -1.036679430885215e-02,

     &      p60 =  4.751256171799112e-05

      if (uref >= 0.0 .and. uref < 5.9 ) then
         znott = p00
      elseif (uref >= 5.9 .and. uref <= 15.4) then
         znott = p10 + uref * (p11 + uref * (p12 + uref * (p13
     &               + uref * (p14 + uref * p15))))
      elseif (uref > 15.4 .and. uref <= 21.6) then
         znott = p20 + uref * (p21 + uref * (p22 + uref * (p23 
     &               + uref * (p24 + uref * p25))))
      elseif (uref > 21.6 .and. uref <= 42.2) then
         znott = p30 + uref * (p31 + uref * (p32 + uref * (p33 
     &               + uref * (p34 + uref * p35))))
      elseif ( uref > 42.2 .and. uref <= 53.3) then
         znott = p40 + uref * (p41 + uref * (p42 + uref * (p43 
     &               + uref * (p44 + uref * p45))))
      elseif ( uref > 53.3 .and. uref <= 80.0) then
         znott = p50 + uref * (p51 + uref * (p52 + uref * (p53 
     &               + uref * (p54 + uref * (p55 + uref * p56)))))
      elseif ( uref > 80.0) then
         znott = p60
      else
         print*, 'Wrong input uref value:',uref
      endif

      END SUBROUTINE znot_t_v6


      SUBROUTINE znot_m_v7(uref, znotm)
      use machine , only : kind_phys
      IMPLICIT NONE
! Calculate areodynamical roughness over water with input 10-m wind
! For low-to-moderate winds, try to match the Cd-U10 relationship from COARE V3.5 (Edson et al. 2013)
! For high winds, try to fit available observational data
! Comparing to znot_t_v6, slightly decrease Cd for higher wind speed
!
! Bin Liu, NOAA/NCEP/EMC 2018
!
! uref(m/s)   :   wind speed at 10-m height
! znotm(meter):   areodynamical roughness scale over water
!

      REAL(kind=kind_phys), INTENT(IN) :: uref
      REAL(kind=kind_phys), INTENT(OUT):: znotm

      real(kind=kind_phys), parameter  :: p13 = -1.296521881682694e-02,
     &      p12 =  2.855780863283819e-01, p11 = -1.597898515251717e+00,
     &      p10 = -8.396975715683501e+00,

     &      p25 =  3.790846746036765e-10, p24 =  3.281964357650687e-09,
     &      p23 =  1.962282433562894e-07, p22 = -1.240239171056262e-06,
     &      p21 =  1.739759082358234e-07, p20 =  2.147264020369413e-05,

     &      p35 =  1.897534489606422e-07, p34 = -3.019495980684978e-05,
     &      p33 =  1.931392924987349e-03, p32 = -6.797293095862357e-02,
     &      p31 =  1.346757797103756e+00, p30 = -1.707846930193362e+01,

     &      p40 =  3.371427455376717e-04

      if (uref >= 0.0 .and.  uref <= 6.5 ) then
        znotm = exp( p10 + uref * (p11 + uref * (p12 + uref * p13)))
      elseif (uref > 6.5 .and. uref <= 15.7) then
        znotm = p20 + uref * (p21 + uref * (p22 + uref * (p23
     &              + uref * (p24 + uref * p25))))
      elseif (uref > 15.7 .and. uref <= 53.0) then
        znotm = exp( p30 + uref * (p31 + uref * (p32 + uref * (p33
     &                   + uref * (p34 + uref * p35)))))
      elseif ( uref > 53.0) then
        znotm = p40
      else
        print*, 'Wrong input uref value:',uref
      endif

      END SUBROUTINE znot_m_v7
      SUBROUTINE znot_t_v7(uref, znott)
      use machine , only : kind_phys
      IMPLICIT NONE
! Calculate scalar roughness over water with input 10-m wind
! For low-to-moderate winds, try to match the Ck-U10 relationship from COARE algorithm
! For high winds, try to retain the Ck-U10 relationship of FY2015 HWRF
! To be compatible with the slightly decreased Cd for higher wind speed
!
! Bin Liu, NOAA/NCEP/EMC 2018
!
! uref(m/s)   :   wind speed at 10-m height
! znott(meter):   scalar roughness scale over water
!

      REAL(kind=kind_phys), INTENT(IN) :: uref
      REAL(kind=kind_phys), INTENT(OUT):: znott

      real(kind=kind_phys), parameter  :: p00 =  1.100000000000000e-04,

     &      p15 = -9.193764479895316e-10, p14 =  7.052217518653943e-08,
     &      p13 = -2.163419217747114e-06, p12 =  3.342963077911962e-05,
     &      p11 = -2.633566691328004e-04, p10 =  8.644979973037803e-04,

     &      p25 = -9.402722450219142e-12, p24 =  1.325396583616614e-09,
     &      p23 = -7.299148051141852e-08, p22 =  1.982901461144764e-06,
     &      p21 = -2.680293455916390e-05, p20 =  1.484341646128200e-04,

     &      p35 =  7.921446674311864e-12, p34 = -1.019028029546602e-09,
     &      p33 =  5.251986927351103e-08, p32 = -1.337841892062716e-06,
     &      p31 =  1.659454106237737e-05, p30 = -7.558911792344770e-05,

     &      p45 = -2.694370426850801e-10, p44 =  5.817362913967911e-08,
     &      p43 = -5.000813324746342e-06, p42 =  2.143803523428029e-04,
     &      p41 = -4.588070983722060e-03, p40 =  3.924356617245624e-02,

     &      p56 = -1.663918773476178e-13, p55 =  6.724854483077447e-11,
     &      p54 = -1.127030176632823e-08, p53 =  1.003683177025925e-06,
     &      p52 = -5.012618091180904e-05, p51 =  1.329762020689302e-03,
     &      p50 = -1.450062148367566e-02, p60 =  6.840803042788488e-05

        if (uref >= 0.0 .and. uref < 5.9 ) then
           znott = p00
        elseif (uref >= 5.9 .and. uref <= 15.4) then
           znott = p10 + uref * (p11 + uref * (p12 + uref * (p13
     &                 + uref * (p14 + uref * p15))))
        elseif (uref > 15.4 .and. uref <= 21.6) then
           znott = p20 + uref * (p21 + uref * (p22 + uref * (p23
     &                 + uref * (p24 + uref * p25))))
        elseif (uref > 21.6 .and. uref <= 42.6) then
           znott = p30 + uref * (p31 + uref * (p32 + uref * (p33
     &                 + uref * (p34 + uref * p35))))
        elseif ( uref > 42.6 .and. uref <= 53.0) then
           znott = p40 + uref * (p41 + uref * (p42 + uref * (p43
     &                 + uref * (p44 + uref * p45))))
        elseif ( uref > 53.0 .and. uref <= 80.0) then
           znott = p50 + uref * (p51 + uref * (p52 + uref * (p53
     &                 + uref * (p54 + uref * (p55 + uref * p56)))))
        elseif ( uref > 80.0) then
           znott = p60
        else
           print*, 'Wrong input uref value:',uref
        endif

        END SUBROUTINE znot_t_v7


!---------------------------------
      end module module_sfc_diff
