      subroutine sfc_diff(im,ps,u1,v1,t1,q1,z1,
     &                    snwdph,tskin,z0rl,cm,ch,rb,
     &                    prsl1,prslki,islimsk,
     &                    stress,fm,fh,
     &                    ustar,wind,ddvel,fm10,fh2,
     &                    sigmaf,vegtype,shdmax,ivegsrc,
     &                    z0pert,ztpert,                        ! mg, sfc-perts
     &                    u10m,v10m,     !wang
     &                    tsurf,flag_iter,redrag)
!
      use machine , only : kind_phys
      use funcphys, only : fpvs
      use physcons, grav => con_g,       cp => con_cp
     &,             rvrdm1 => con_fvirt, rd => con_rd
     &,             eps => con_eps, epsm1 => con_epsm1

      implicit none
!
      integer              im, ivegsrc
      real(kind=kind_phys), dimension(im) :: ps,  u1, v1, t1, q1, z1
     &,                                      tskin, z0rl, cm,  ch, rb
     &,                                      prsl1, prslki, stress
     &,                                      fm, fh, ustar, wind, ddvel
     &,                                      fm10, fh2, sigmaf, shdmax
     &,                                      tsurf, snwdph
     &,                                      z0pert,ztpert               ! mg, sfc-perts
      real(kind=kind_phys), dimension(im) :: u10m,v10m, wind10m

      integer, dimension(im)              ::  vegtype, islimsk

      logical   flag_iter(im) ! added by s.lu
      logical   redrag        ! reduced drag coeff. flag for high wind over sea (j.han)
!
!     locals
!

       logical run_tc    ! wang, use hurricane-obs-based z0
      integer   i
!
      real(kind=kind_phys) aa,     aa0,    bb,     bb0, dtv,   adtv,qs1,
     &                     hl1,    hl12,   pm,     ph,  pm10,  ph2, rat,
     &                     thv1,   tvs,    z1i,    z0,  z0max, ztmax,
     &                     fms,    fhs,    hl0,    hl0inf, hlinf,
     &                     hl110,  hlt,    hltinf, olinf,
     &                     restar, czilc,  tem1,   tem2, ztmax1
!
      real(kind=kind_phys), parameter ::
     &              charnock=.014, ca=.4  ! ca - von karman constant
     &,             z0s_max=.317e-2       ! a limiting value at high winds over sea
! Jili Dong modify z0s_max
!     &,             z0s_max=.196e-2       ! a limiting value at high winds over sea                     


     &,             alpha=5.,   a0=-3.975, a1=12.32, alpha4=4.0*alpha
     &,             b1=-7.755,  b2=6.041,  alpha2=alpha+alpha, beta=1.0
     &,             a0p=-7.941, a1p=24.75, b1p=-8.705, b2p=7.899
     &,             vis=1.4e-5, rnu=1.51e-5, visi=1.0/vis

     &,             log01=log(0.01), log05=log(0.05), log07=log(0.07)
     &,             ztmin1=-999.0

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

! Weiguo Wang added 20190425
!      run_tc=.true.       ! use obs-based roughness length ~ 10m wind
      run_tc=.false.     ! not

      do i=1,im
        if(flag_iter(i)) then 
          wind(i) = max(sqrt(u1(i)*u1(i) + v1(i)*v1(i))
     &                + max(0.0, min(ddvel(i), 30.0)), 1.0)
          tem1    = 1.0 + rvrdm1 * max(q1(i),1.e-8)
          thv1    = t1(i) * prslki(i) * tem1
          tvs     = 0.5 * (tsurf(i)+tskin(i)) * tem1
          qs1     = fpvs(t1(i))
          qs1     = max(1.0e-8, eps * qs1 / (prsl1(i) + epsm1 * qs1))

          z0      = 0.01 * z0rl(i)
          z0max   = max(1.0e-6, min(z0,z1(i)))
          z1i     = 1.0 / z1(i)

          wind10m(i) = max(sqrt( u10m(i)*u10m(i) + v10m(i)*v10m(i) ),
     &                 1.0)
!  compute stability dependent exchange coefficients
!  this portion of the code is presently suppressed
!

          if(islimsk(i) == 0) then            ! over ocean
            ustar(i) = sqrt(grav * z0 / charnock)

!**  test xubin's new z0

!           ztmax  = z0max

            restar = max(ustar(i)*z0max*visi, 0.000001)

!           restar = log(restar)
!           restar = min(restar,5.)
!           restar = max(restar,-5.)
!           rat    = aa1 + (bb1 + cc1*restar) * restar
!           rat    = rat    / (1. + (bb2 + cc2*restar) * restar))
!  rat taken from zeng, zhao and dickinson 1997

            rat    = min(7.0, 2.67 * sqrt(sqrt(restar)) - 2.57)
            ztmax  = z0max * exp(-rat)

! Weiguo Wang, 2019-0425, use fitted zt from obs
            if (run_tc) call znot_t_v7(wind10m(i),ztmax)   ! 10-m wind, m/s, ztmax, m
!

          else                                ! over land and sea ice
!** xubin's new z0  over land and sea ice
            tem1 = 1.0 - shdmax(i)
            tem2 = tem1 * tem1
            tem1 = 1.0  - tem2
          
            if( ivegsrc == 1 ) then

              if (vegtype(i) == 10) then
                z0max = exp( tem2*log01 + tem1*log07 )
              elseif (vegtype(i) == 6) then
                z0max = exp( tem2*log01 + tem1*log05 )
              elseif (vegtype(i) == 7) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max = 0.01
              elseif (vegtype(i) == 16) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max = 0.01
              else
                z0max = exp( tem2*log01 + tem1*log(z0max) )
              endif

            elseif (ivegsrc == 2 ) then

                if (vegtype(i) == 7) then
                  z0max = exp( tem2*log01 + tem1*log07 )
                elseif (vegtype(i) == 8) then
                  z0max = exp( tem2*log01 + tem1*log05 )
                elseif (vegtype(i) == 9) then
!                 z0max = exp( tem2*log01 + tem1*log01 )
                  z0max = 0.01
                elseif (vegtype(i) == 11) then
!                 z0max = exp( tem2*log01 + tem1*log01 )
                  z0max = 0.01
                else
                  z0max = exp( tem2*log01 + tem1*log(z0max) )
                endif

            endif


! mg, sfc-perts: add surface perturbations to z0max over land
            if ( islimsk(i) == 1 .and. z0pert(i) /= 0.0 ) then
              z0max = z0max * (10.**z0pert(i))
            endif
 
            z0max = max(z0max,1.0e-6)

!           czilc = 10.0 ** (- (0.40/0.07) * z0) ! fei's canopy height dependance of czil
            czilc = 0.8

            tem1 = 1.0 - sigmaf(i)
            ztmax = z0max*exp( - tem1*tem1
     &                         * czilc*ca*sqrt(ustar(i)*(0.01/1.5e-05)))

! mg, sfc-perts: add surface perturbations to ztmax/z0max ratio over land
            if ( islimsk(i) == 1  .and. ztpert(i) /= 0.0) then
              ztmax = ztmax * (10.**ztpert(i))
            endif


          endif       ! end of if(islimsk(i) == 0) then

          ztmax  = max(ztmax,1.0e-6)
          tem1   = z0max/z1(i)
          if (abs(1.0-tem1) > 1.0e-6) then
            ztmax1 = - beta*log(tem1)/(alpha2*(1.-tem1))
          else
            ztmax1 = 99.0
          endif
          if( z0max < 0.05 .and. snwdph(i) < 10.0 ) ztmax1 = 99.0


!  compute stability indices (rb and hlinf)

          dtv     = thv1 - tvs
          adtv    = max(abs(dtv),0.001)
          dtv     = sign(1.,dtv) * adtv
          rb(i)   = max(-5000.0, (grav+grav) * dtv * z1(i)
     &            / ((thv1 + tvs) * wind(i) * wind(i)))
          tem1    = 1.0 / z0max
          tem2    = 1.0 / ztmax
          fm(i)   = log((z0max+z1(i)) * tem1)
          fh(i)   = log((ztmax+z1(i)) * tem2)
          fm10(i) = log((z0max+10.)   * tem1)
          fh2(i)  = log((ztmax+2.)    * tem2)
          hlinf   = rb(i) * fm(i) * fm(i) / fh(i)
          hlinf   = min(max(hlinf,ztmin1),ztmax1)
!
!  stable case
!
          if (dtv >= 0.0) then
            hl1 = hlinf
            if(hlinf > .25) then
              tem1   = hlinf * z1i
              hl0inf = z0max * tem1
              hltinf = ztmax * tem1
              aa     = sqrt(1. + alpha4 * hlinf)
              aa0    = sqrt(1. + alpha4 * hl0inf)
              bb     = aa
              bb0    = sqrt(1. + alpha4 * hltinf)
              pm     = aa0 - aa + log( (aa + 1.)/(aa0 + 1.) )
              ph     = bb0 - bb + log( (bb + 1.)/(bb0 + 1.) )
              fms    = fm(i) - pm
              fhs    = fh(i) - ph
              hl1    = fms * fms * rb(i) / fhs
              hl1    = min(max(hl1, ztmin1), ztmax1)
            endif
!
!  second iteration
!
            tem1  = hl1 * z1i
            hl0   = z0max * tem1
            hlt   = ztmax * tem1
            aa    = sqrt(1. + alpha4 * hl1)
            aa0   = sqrt(1. + alpha4 * hl0)
            bb    = aa
            bb0   = sqrt(1. + alpha4 * hlt)
            pm    = aa0 - aa + log( (1.0+aa)/(1.0+aa0) )
            ph    = bb0 - bb + log( (1.0+bb)/(1.0+bb0) )
            hl110 = hl1 * 10. * z1i
            hl110 = min(max(hl110, ztmin1), ztmax1)
            aa    = sqrt(1. + alpha4 * hl110)
            pm10  = aa0 - aa + log( (1.0+aa)/(1.0+aa0) )
            hl12  = (hl1+hl1) * z1i
            hl12  = min(max(hl12,ztmin1),ztmax1)
!           aa    = sqrt(1. + alpha4 * hl12)
            bb    = sqrt(1. + alpha4 * hl12)
            ph2   = bb0 - bb + log( (1.0+bb)/(1.0+bb0) )
!
!  unstable case - check for unphysical obukhov length
!
          else                          ! dtv < 0 case
            olinf = z1(i) / hlinf
            tem1  = 50.0 * z0max
            if(abs(olinf) <= tem1) then
              hlinf = -z1(i) / tem1
              hlinf = min(max(hlinf,ztmin1),ztmax1)
            endif
!
!  get pm and ph
!
            if (hlinf >= -0.5) then
              hl1   = hlinf
              pm    = (a0  + a1*hl1)  * hl1   / (1.+ (b1+b2*hl1)  *hl1)
              ph    = (a0p + a1p*hl1) * hl1   / (1.+ (b1p+b2p*hl1)*hl1)
              hl110 = hl1 * 10. * z1i
              hl110 = min(max(hl110, ztmin1), ztmax1)
              pm10  = (a0 + a1*hl110) * hl110 / (1.+(b1+b2*hl110)*hl110)
              hl12  = (hl1+hl1) * z1i
              hl12  = min(max(hl12, ztmin1), ztmax1)
              ph2   = (a0p + a1p*hl12) * hl12 / (1.+(b1p+b2p*hl12)*hl12)
            else                       ! hlinf < 0.05
              hl1   = -hlinf
              tem1  = 1.0 / sqrt(hl1)
              pm    = log(hl1) + 2. * sqrt(tem1) - .8776
              ph    = log(hl1) + .5 * tem1 + 1.386
!             pm    = log(hl1) + 2.0 * hl1 ** (-.25) - .8776
!             ph    = log(hl1) + 0.5 * hl1 ** (-.5) + 1.386
              hl110 = hl1 * 10. * z1i
              hl110 = min(max(hl110, ztmin1), ztmax1)
              pm10  = log(hl110) + 2.0 / sqrt(sqrt(hl110)) - .8776
!             pm10  = log(hl110) + 2. * hl110 ** (-.25) - .8776
              hl12  = (hl1+hl1) * z1i
              hl12  = min(max(hl12, ztmin1), ztmax1)
              ph2   = log(hl12) + 0.5 / sqrt(hl12) + 1.386
!             ph2   = log(hl12) + .5 * hl12 ** (-.5) + 1.386
            endif

          endif          ! end of if (dtv >= 0 ) then loop
!
!  finish the exchange coefficient computation to provide fm and fh
!
          fm(i)     = fm(i) - pm
          fh(i)     = fh(i) - ph
          fm10(i)   = fm10(i) - pm10
          fh2(i)    = fh2(i) - ph2
          cm(i)     = ca * ca / (fm(i) * fm(i))
          ch(i)     = ca * ca / (fm(i) * fh(i))
          tem1      = 0.00001/z1(i)
          cm(i)     = max(cm(i), tem1)
          ch(i)     = max(ch(i), tem1)
          stress(i) = cm(i) * wind(i) * wind(i)
          ustar(i)  = sqrt(stress(i))
!
!  update z0 over ocean
!
          if(islimsk(i) == 0) then
            z0 = (charnock / grav) * ustar(i) * ustar(i)

! mbek -- toga-coare flux algorithm
!           z0 = (charnock / grav) * ustar(i)*ustar(i) +  arnu/ustar(i)
!  new implementation of z0
!           cc = ustar(i) * z0 / rnu
!           pp = cc / (1. + cc)
!           ff = grav * arnu / (charnock * ustar(i) ** 3)
!           z0 = arnu / (ustar(i) * ff ** pp)

            if (redrag) then
              z0rl(i) = 100.0 * max(min(z0, z0s_max), 1.e-7)
            else
              z0rl(i) = 100.0 * max(min(z0,.1), 1.e-7)
            endif
!! Weiguo Wang
          if (run_tc) then 
           call znot_m_v7(wind10m(i),z0)   ! wind, m/s, z0, m 
           z0rl(i) = 100.0 * z0        ! in cm 
          endif
!!
          endif
        endif                ! end of if(flagiter) loop
      enddo

      return
      end

!! add fitted z0,zt curves for hurricane application (used in HWRF/HMON)
!! Weiguo Wang, 2019-0425
       SUBROUTINE znot_m_v7(uref,znotm)
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

      REAL, INTENT(IN) :: uref
      REAL, INTENT(OUT):: znotm
      REAL             :: p13, p12, p11, p10
      REAL             :: p25, p24, p23, p22, p21, p20
      REAL             :: p35, p34, p33, p32, p31, p30
      REAL             :: p40

       p13 = -1.296521881682694e-02
       p12 =  2.855780863283819e-01
       p11 = -1.597898515251717e+00
       p10 = -8.396975715683501e+00

       p25 =  3.790846746036765e-10
       p24 =  3.281964357650687e-09
       p23 =  1.962282433562894e-07
       p22 = -1.240239171056262e-06
       p21 =  1.739759082358234e-07
       p20 =  2.147264020369413e-05


       p35 =  1.897534489606422e-07
       p34 = -3.019495980684978e-05
       p33 =  1.931392924987349e-03
       p32 = -6.797293095862357e-02
       p31 =  1.346757797103756e+00
       p30 = -1.707846930193362e+01

       p40 =  3.371427455376717e-04

       if (uref >= 0.0 .and.  uref <= 6.5 ) then
        znotm = exp( p10 + p11*uref + p12*uref**2 + p13*uref**3)
       elseif (uref > 6.5 .and. uref <= 15.7) then
        znotm = p25*uref**5 + p24*uref**4 + p23*uref**3 + 
     &          p22*uref**2 + p21*uref + p20
       elseif (uref > 15.7 .and. uref <= 53.0) then
        znotm = exp( p35*uref**5 + p34*uref**4 + p33*uref**3 
     &          + p32*uref**2 + p31*uref + p30 )
       elseif ( uref > 53.0) then
        znotm = p40
       else
        print*, 'Wrong input uref value:',uref
       endif

      END SUBROUTINE znot_m_v7
      SUBROUTINE znot_t_v7(uref,znott)
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

        REAL, INTENT(IN) :: uref
        REAL, INTENT(OUT):: znott

        REAL             :: p00
        REAL             :: p15, p14, p13, p12, p11, p10
        REAL             :: p25, p24, p23, p22, p21, p20
        REAL             :: p35, p34, p33, p32, p31, p30
        REAL             :: p45, p44, p43, p42, p41, p40
        REAL             :: p56, p55, p54, p53, p52, p51, p50
        REAL             :: p60

         p00 =  1.100000000000000e-04

         p15 = -9.193764479895316e-10
          p14 =  7.052217518653943e-08
         p13 = -2.163419217747114e-06
         p12 =  3.342963077911962e-05
         p11 = -2.633566691328004e-04
         p10 =  8.644979973037803e-04

         p25 = -9.402722450219142e-12
         p24 =  1.325396583616614e-09
         p23 = -7.299148051141852e-08
         p22 =  1.982901461144764e-06
         p21 = -2.680293455916390e-05
         p20 =  1.484341646128200e-04

         p35 =  7.921446674311864e-12
         p34 = -1.019028029546602e-09
         p33 =  5.251986927351103e-08
         p32 = -1.337841892062716e-06
         p31 =  1.659454106237737e-05
         p30 = -7.558911792344770e-05

         p45 = -2.694370426850801e-10
          p44 =  5.817362913967911e-08
         p43 = -5.000813324746342e-06
         p42 =  2.143803523428029e-04
         p41 = -4.588070983722060e-03
         p40 =  3.924356617245624e-02

        p56 = -1.663918773476178e-13
        p55 =  6.724854483077447e-11
        p54 = -1.127030176632823e-08
         p53 =  1.003683177025925e-06
        p52 = -5.012618091180904e-05
        p51 =  1.329762020689302e-03
        p50 = -1.450062148367566e-02

        p60 =  6.840803042788488e-05

        if (uref >= 0.0 .and. uref < 5.9 ) then
            znott = p00
         elseif (uref >= 5.9 .and. uref <= 15.4) then
           znott = p15*uref**5 + p14*uref**4 + p13*uref**3 + 
     &             p12*uref**2 + p11*uref + p10
         elseif (uref > 15.4 .and. uref <= 21.6) then
           znott = p25*uref**5 + p24*uref**4 + p23*uref**3 + 
     &             p22*uref**2 + p21*uref + p20
         elseif (uref > 21.6 .and. uref <= 42.6) then
           znott = p35*uref**5 + p34*uref**4 + p33*uref**3 + 
     &             p32*uref**2 + p31*uref + p30
         elseif ( uref > 42.6 .and. uref <= 53.0) then
           znott = p45*uref**5 + p44*uref**4 + p43*uref**3 + 
     &             p42*uref**2 + p41*uref + p40
         elseif ( uref > 53.0 .and. uref <= 80.0) then
           znott = p56*uref**6 + p55*uref**5 + p54*uref**4 + 
     &             p53*uref**3 + p52*uref**2 + p51*uref + p50
         elseif ( uref > 80.0) then
           znott = p60
        else
           print*, 'Wrong input uref value:',uref
         endif

        END SUBROUTINE znot_t_v7


 
  
