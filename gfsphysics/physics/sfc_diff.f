      module module_sfc_diff

      use machine , only : kind_phys
      use physcons, grav => con_g
      real (kind=kind_phys), parameter :: ca=.4  ! ca - von karman constant

      contains
      subroutine sfc_diff(im,ps,u1,v1,t1,q1,z1,                 !intent(in)
     &                    prsl1,prslki,ddvel,                   !intent(in)
     &                    sigmaf,vegtype,shdmax,ivegsrc,        !intent(in)
     &                    z0pert,ztpert,        ! mg, sfc-perts !intent(in)
     &                    flag_iter,redrag,                     !intent(in)
     &                    wet,dry,icy,fice,                     !intent(in)
     &                    tskin_ocn, tskin_lnd, tskin_ice,      !intent(in)
     &                    tsurf_ocn, tsurf_lnd, tsurf_ice,      !intent(in)
     &                   snwdph_ocn,snwdph_lnd,snwdph_ice,      !intent(in)
!
     &                     z0rl_ocn,  z0rl_lnd,  z0rl_ice,      !intent(inout)
     &                    ustar_ocn, ustar_lnd, ustar_ice,      !intent(inout)
!
     &                       cm_ocn,    cm_lnd,    cm_ice,      !intent(out)
     &                       ch_ocn,    ch_lnd,    ch_ice,      !intent(out)
     &                       rb_ocn,    rb_lnd,    rb_ice,      !intent(out)
     &                   stress_ocn,stress_lnd,stress_ice,      !intent(out)
     &                       fm_ocn,    fm_lnd,    fm_ice,      !intent(out)
     &                       fh_ocn,    fh_lnd,    fh_ice,      !intent(out)
     &                     fm10_ocn,  fm10_lnd,  fm10_ice,      !intent(out)
     &                      fh2_ocn,   fh2_lnd,   fh2_ice,      !intent(out)
     &                      wind)                               !intent(out)
!
      use funcphys, only : fpvs
      use physcons, rvrdm1 => con_fvirt
     &,             eps => con_eps, epsm1 => con_epsm1
      implicit none
!
      integer, intent(in) :: im, ivegsrc
      integer, dimension(im), intent(in) :: vegtype

      logical, intent(in) :: redrag ! reduced drag coeff. flag for high wind over sea (j.han)
      logical, dimension(im), intent(in) :: flag_iter, dry, wet, icy

      real(kind=kind_phys), dimension(im), intent(in)    ::
     &                    ps,u1,v1,t1,q1,z1,prsl1,prslki,ddvel,
     &                    sigmaf,shdmax,
     &                    z0pert,ztpert ! mg, sfc-perts
      real(kind=kind_phys), dimension(im), intent(in)    ::
     &                    tskin_ocn, tskin_lnd, tskin_ice,
     &                    tsurf_ocn, tsurf_lnd, tsurf_ice,
     &                   snwdph_ocn,snwdph_lnd,snwdph_ice,
     &                    fice

      real(kind=kind_phys), dimension(im), intent(inout) ::
     &                     z0rl_ocn,  z0rl_lnd,  z0rl_ice,
     &                    ustar_ocn, ustar_lnd, ustar_ice

      real(kind=kind_phys), dimension(im), intent(out)   ::
     &                       cm_ocn,    cm_lnd,    cm_ice,
     &                       ch_ocn,    ch_lnd,    ch_ice,
     &                       rb_ocn,    rb_lnd,    rb_ice,
     &                   stress_ocn,stress_lnd,stress_ice,
     &                       fm_ocn,    fm_lnd,    fm_ice,
     &                       fh_ocn,    fh_lnd,    fh_ice,
     &                     fm10_ocn,  fm10_lnd,  fm10_ice,
     &                      fh2_ocn,   fh2_lnd,   fh2_ice,
     &                      wind
!
!     locals
!
      integer   i
!
      real(kind=kind_phys) :: qs1,  rat, thv1, restar,
     &                      czilc, tem1, tem2

      real(kind=kind_phys) :: tvs_ocn,  tvs_lnd,  tvs_ice,
     &                         z0_ocn,   z0_lnd,   z0_ice,
     &                      z0max_ocn,z0max_lnd,z0max_ice,
     &                      ztmax_ocn,ztmax_lnd,ztmax_ice
!
      real(kind=kind_phys), parameter ::
     &              charnock=.014, z0s_max=.317e-2       ! a limiting value at high winds over sea
     &,             vis=1.4e-5, rnu=1.51e-5, visi=1.0/vis
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
      do i=1,im
        ztmax_ocn = 0.; ztmax_lnd = 0.; ztmax_ice = 0.
        if(flag_iter(i)) then 
          wind(i) = max(sqrt(u1(i)*u1(i) + v1(i)*v1(i))
     &                + max(0.0, min(ddvel(i), 30.0)), 1.0)
          tem1    = 1.0 + rvrdm1 * max(q1(i),1.e-8)
          thv1    = t1(i) * prslki(i) * tem1
          tvs_ocn = 0.5 * (tsurf_ocn(i)+tskin_ocn(i)) * tem1
          tvs_lnd = 0.5 * (tsurf_lnd(i)+tskin_lnd(i)) * tem1
          tvs_ice = 0.5 * (tsurf_ice(i)+tskin_ice(i)) * tem1
          qs1     = fpvs(t1(i))
          qs1     = max(1.0e-8, eps * qs1 / (prsl1(i) + epsm1 * qs1))

          z0_ocn      = 0.01 * z0rl_ocn(i)
          z0max_ocn   = max(1.0e-6, min(z0_ocn,z1(i)))
          z0_lnd      = 0.01 * z0rl_lnd(i)
          z0max_lnd   = max(1.0e-6, min(z0_lnd,z1(i)))
          z0_ice      = 0.01 * z0rl_ice(i)
          z0max_ice   = max(1.0e-6, min(z0_ice,z1(i)))

!  compute stability dependent exchange coefficients
!  this portion of the code is presently suppressed
!

          if (wet(i) .and. fice(i) < 1.) then ! some open ocean
            ustar_ocn(i) = sqrt(grav * z0_ocn / charnock)

!**  test xubin's new z0

!           ztmax  = z0max

            restar = max(ustar_ocn(i)*z0max_ocn*visi, 0.000001)

!           restar = log(restar)
!           restar = min(restar,5.)
!           restar = max(restar,-5.)
!           rat    = aa1 + (bb1 + cc1*restar) * restar
!           rat    = rat    / (1. + (bb2 + cc2*restar) * restar))
!  rat taken from zeng, zhao and dickinson 1997

            rat    = min(7.0, 2.67 * sqrt(sqrt(restar)) - 2.57)
            ztmax_ocn  = z0max_ocn * exp(-rat)
          endif ! Open ocean
          if (dry(i) .or. icy(i)) then ! over land or sea ice
!** xubin's new z0  over land and sea ice
            tem1 = 1.0 - shdmax(i)
            tem2 = tem1 * tem1
            tem1 = 1.0  - tem2
          
            if( ivegsrc == 1 ) then

              if (vegtype(i) == 10) then
                z0max_lnd = exp( tem2*log01 + tem1*log07 )
              elseif (vegtype(i) == 6) then
                z0max_lnd = exp( tem2*log01 + tem1*log05 )
              elseif (vegtype(i) == 7) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max_lnd = 0.01
              elseif (vegtype(i) == 16) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max_lnd = 0.01
              else
                z0max_lnd = exp( tem2*log01 + tem1*log(z0max_lnd) )
              endif

            elseif (ivegsrc == 2 ) then

                if (vegtype(i) == 7) then
                  z0max_lnd = exp( tem2*log01 + tem1*log07 )
                elseif (vegtype(i) == 8) then
                  z0max_lnd = exp( tem2*log01 + tem1*log05 )
                elseif (vegtype(i) == 9) then
!                 z0max = exp( tem2*log01 + tem1*log01 )
                  z0max_lnd = 0.01
                elseif (vegtype(i) == 11) then
!                 z0max = exp( tem2*log01 + tem1*log01 )
                  z0max_lnd = 0.01
                else
                  z0max_lnd = exp( tem2*log01 + tem1*log(z0max_lnd) )
                endif

            endif ! over land or sea ice

            z0max_ice = z0max_lnd

! mg, sfc-perts: add surface perturbations to z0max over land
            if (dry(i) .and. z0pert(i) /= 0.0 ) then
              z0max_lnd = z0max_lnd * (10.**z0pert(i))
            endif
 
            z0max_lnd = max(z0max_lnd,1.0e-6)
            z0max_ice = max(z0max_ice,1.0e-6)

!           czilc = 10.0 ** (- (0.40/0.07) * z0) ! fei's canopy height dependance of czil
            czilc = 0.8

            tem1 = 1.0 - sigmaf(i)
            ztmax_lnd = z0max_lnd*exp( - tem1*tem1
     &                     * czilc*ca*sqrt(ustar_lnd(i)*(0.01/1.5e-05)))
            ztmax_ice = z0max_ice*exp( - tem1*tem1
     &                     * czilc*ca*sqrt(ustar_ice(i)*(0.01/1.5e-05)))


! mg, sfc-perts: add surface perturbations to ztmax/z0max ratio over land
            if (dry(i) .and. ztpert(i) /= 0.0) then
              ztmax_lnd = ztmax_lnd * (10.**ztpert(i))
            endif


          endif       ! end of if(sfctype flags) then

          ztmax_ocn  = max(ztmax_ocn,1.0e-6)
          ztmax_lnd  = max(ztmax_lnd,1.0e-6)
          ztmax_ice  = max(ztmax_ice,1.0e-6)

! BWG begin "stability" block, 2019-03-23
      if (wet(i) .and. fice(i) < 1.) then ! Some open ocean
          call stability
!  ---  inputs:                                                  
     &     (z1(i),snwdph_ocn(i),thv1,wind(i),
     &      z0max_ocn,ztmax_ocn,tvs_ocn,
!  ---  outputs:
     &      rb_ocn(i),fm_ocn(i),fh_ocn(i),fm10_ocn(i),fh2_ocn(i),
     &      cm_ocn(i),ch_ocn(i),stress_ocn(i),ustar_ocn(i))
      endif ! Open ocean points

      if (dry(i)) then ! Some land
          call stability
!  ---  inputs:                                                  
     &     (z1(i),snwdph_lnd(i),thv1,wind(i),
     &      z0max_lnd,ztmax_lnd,tvs_lnd,
!  ---  outputs:
     &      rb_lnd(i),fm_lnd(i),fh_lnd(i),fm10_lnd(i),fh2_lnd(i),
     &      cm_lnd(i),ch_lnd(i),stress_lnd(i),ustar_lnd(i))
      endif ! Dry points

      if (icy(i)) then ! Some ice
          call stability
!  ---  inputs:                                                  
     &     (z1(i),snwdph_ice(i),thv1,wind(i),
     &      z0max_ice,ztmax_ice,tvs_ice,
!  ---  outputs:
     &      rb_ice(i),fm_ice(i),fh_ice(i),fm10_ice(i),fh2_ice(i),
     &      cm_ice(i),ch_ice(i),stress_ice(i),ustar_ice(i))
      endif ! Icy points

! BWG: Everything from here to end of subroutine was after
!      the stuff now put into "stability"

!
!  update z0 over ocean
!
          if (wet(i) .and. fice(i) < 1.) then
            z0_ocn = (charnock / grav) * ustar_ocn(i) * ustar_ocn(i)

! mbek -- toga-coare flux algorithm
!           z0 = (charnock / grav) * ustar(i)*ustar(i) +  arnu/ustar(i)
!  new implementation of z0
!           cc = ustar(i) * z0 / rnu
!           pp = cc / (1. + cc)
!           ff = grav * arnu / (charnock * ustar(i) ** 3)
!           z0 = arnu / (ustar(i) * ff ** pp)

            if (redrag) then
              z0rl_ocn(i) = 100.0 * max(min(z0_ocn, z0s_max), 1.e-7)
            else
              z0rl_ocn(i) = 100.0 * max(min(z0_ocn,.1), 1.e-7)
            endif
          endif              ! end of if(open ocean)
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
      real(kind=kind_phys), parameter :: alpha=5., a0=-3.975            &
     &,             a1=12.32, alpha4=4.0*alpha
     &,             b1=-7.755,  b2=6.041,  alpha2=alpha+alpha, beta=1.0
     &,             a0p=-7.941, a1p=24.75, b1p=-8.705, b2p=7.899
     &,             ztmin1=-999.0

      real(kind=kind_phys) aa,     aa0,    bb,     bb0, dtv,   adtv,
     &                     hl1,    hl12,   pm,     ph,  pm10,  ph2,
     &                     z1i,
     &                     fms,    fhs,    hl0,    hl0inf, hlinf,
     &                     hl110,  hlt,    hltinf, olinf,
     &                     tem1,   tem2, ztmax1

          z1i = 1.0 / z1

          tem1   = z0max/z1
          if (abs(1.0-tem1) > 1.0e-6) then
            ztmax1 = - beta*log(tem1)/(alpha2*(1.-tem1))
          else
            ztmax1 = 99.0
          endif
          if( z0max < 0.05 .and. snwdph < 10.0 ) ztmax1 = 99.0

!  compute stability indices (rb and hlinf)

          dtv     = thv1 - tvs
          adtv    = max(abs(dtv),0.001)
          dtv     = sign(1.,dtv) * adtv
          rb      = max(-5000.0, (grav+grav) * dtv * z1
     &            / ((thv1 + tvs) * wind * wind))
          tem1    = 1.0 / z0max
          tem2    = 1.0 / ztmax
          fm      = log((z0max+z1) * tem1)
          fh      = log((ztmax+z1) * tem2)
          fm10    = log((z0max+10.)   * tem1)
          fh2     = log((ztmax+2.)    * tem2)
          hlinf   = rb * fm * fm / fh
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
            olinf = z1 / hlinf
            tem1  = 50.0 * z0max
            if(abs(olinf) <= tem1) then
              hlinf = -z1 / tem1
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
          fm        = fm - pm
          fh        = fh - ph
          fm10      = fm10 - pm10
          fh2       = fh2 - ph2
          cm        = ca * ca / (fm * fm)
          ch        = ca * ca / (fm * fh)
          tem1      = 0.00001/z1
          cm        = max(cm, tem1)
          ch        = max(ch, tem1)
          stress    = cm * wind * wind
          ustar     = sqrt(stress)

      return
!.................................
      end subroutine stability
!---------------------------------

!---------------------------------
      end module module_sfc_diff
