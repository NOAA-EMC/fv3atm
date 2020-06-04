       module module_sfc_drv
       contains
!>  \file sfc_drv.f
!!  This file contains the NOAH land surface scheme.
!> \defgroup NOAH NOAH Land Surface
!! @{
!!
!!  The Noah LSM (Chen et al., 1996; Koren et al., 1999; Ek et al., 2003) is targeted for moderate complexity and good computational efficiency for numerical weather prediction and climate models. Thus, it omits subgrid surface tiling and uses a single-layer snowpack. The surface energy balance is solved via a Penman-based approximation for latent heat flux. The Noah model includes packages to simulate soil moisture, soil ice, soil temperature, skin temperature, snow depth, snow water equivalent, energy fluxes such as latent heat, sensible heat and ground heat, and water fluxes such as evaporation and total runoff. The Noah surface infiltration scheme follows that of Schaake et al. (1996) for its treatment of the subgrid variability of precipitation and soil moisture.
!!
!!  On 31 May and 14 June 2005, NCEP extensively upgraded the land-surface component of its Global Forecast System (GFS), including its Global Data Assimilation System (GDAS). The Noah LSM upgrade includes an increase from two (10, 190 cm thick) to four soil layers (10, 30, 60, 100 cm thick), addition of frozen soil physics, new formulations for infiltration and runoff (giving more runoff for unsaturated soils), revised physics of the snowpack and its influence on surface heat fluxes and albedo, tuning and adding canopy resistance parameters, allowing spatially varying root depth, revised treatment of ground heat flux and soil thermal conductivity, reformulation for dependence of direct surface evaporation on first layer soil moisture, and improved seasonality of green vegetation cover. The frozen soil physics includes soil heat sinks/sources from freezing/thawing and influences vertical transport of soil moisture, soil thermal conductivity and heat capacity, and surface infiltration. The prognostic states of snowpack depth and liquid soil moisture were added to the already present prognostic states of snowpack water-equivalent (SWE), total soil moisture (liquid plus frozen), soil temperature, canopy water, and skin temperature. SWE divided by the snowpack depth gives the snowpack density. Total soil moisture minus liquid soil moisture gives the frozen soil moisture (Mitchell et al. 2005)
!!
!!  The addition of Noah LSM greatly reduced the two prominent biases in land-surface processes: 1) an early depletion of snowpack; and 2) a high bias in both surface evaporation and precipitation in the warm season in non-arid mid-latitudes. However, a lower tropospheric warm bias as well as increased surface sensible heat flux emerged, particularly over the arid areas during the daytime. Extensive tests attributed this bias mainly to improper treatment of the thermal roughness length. In May 2011, a new thermal roughness length formulation, which assigned a smaller value for the thermal roughness length compared to the momentum roughness length, was implemented. This greatly reduced the warm surface air temperature bias and the cold skin temperature bias over the arid areas during the daytime (Wei et al. 2009; Zheng et al. 2012).
!!
!!  In January 2015, CFS/GLDAS soil moisture climatology at T574 was used for soil moisture nudge to replace the out-of-date coarse resolution bucket soil moisture climatology; a dependence of the ratio of the thermal and momentum roughness on vegetation type was added to address the land-atmosphere coupling strength; a look-up table based on vegetation type was used to replace 1.0 degree momentum roughness length climatology. After this implementation summer warm/dry biases were found over cropland/grassland areas. Some evaporation-related parameters were refined to increase the evaporation to address this issue. The refinement was implemented in May 2016.
!!
!!  In July 2017, new high-resolution MODIS-based snow-free albedo, maximum snow albedo, soil type and vegetation type were used to address the cold biases over the snow area and the blockiness of surface fields due to the coarse resolution data of soil type and vegetation type. The surface layer parameterization scheme was upgraded to modify the roughness-length formulation and introduce a stability parameter constraint in the Monin-Obukhov similarity theory to prevent the land-atmosphere system from decoupling which causes the rapid temperature drop during the sunset (Zheng et al. 2017).
!!
!!  \section diagram Calling Hierarchy Diagram
!!  \section intraphysics Intraphysics Communication
!!
!> \brief Brief description of the subroutine
!!
!!
!! \section arg_table_Noah_run Arguments
!! | local var name | longname                                           | description                        | units   | rank | type    |    kind   | intent | optional |
!! |----------------|----------------------------------------------------|------------------------------------|---------|------|---------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                             | horizontal loop extent, start at 1 | index   |    0 | integer |           | in     | F        |
!!
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!!  @{
! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!      call sfc_drv                                                     !
!  ---  inputs:                                                         !
!          ( im, km, ps, t1, q1, soiltyp, vegtype, sigmaf,              !
!            sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm, ch,          !
!            prsl1, prslki, zf, land, wind,  slopetyp,                  !
!            shdmin, shdmax, snoalb, sfalb, flag_iter, flag_guess,      !
!            lheatstrg, isot, ivegsrc,                                  !
!  ---  in/outs:                                                        !
!            weasd, snwdph, tskin, tprcp, srflag, smc, stc, slc,        !
!            canopy, trans, tsurf, zorl,                                !
!  ---  outputs:                                                        !
!            sncovr1, qsurf, gflux, drain, evap, hflx, ep, runoff,      !
!            cmm, chh, evbs, evcw, sbsno, snowc, stm, snohf,            !
!            smcwlt2, smcref2, wet1 )                                   !
!                                                                       !
!                                                                       !
!  subprogram called:  sflx                                             !
!                                                                       !
!  program history log:                                                 !
!         xxxx  --             created                                  !
!         200x  -- sarah lu    modified                                 !
!    oct  2006  -- h. wei      modified                                 !
!    apr  2009  -- y.-t. hou   modified to include surface emissivity   !
!                     effect on lw radiation. replaced the comfussing   !
!                     slrad (net sw + dlw) with sfc net sw snet=dsw-usw !
!    sep  2009  -- s. moorthi modification to remove rcl and unit change!
!    nov  2011  -- sarah lu    corrected wet1 calculation
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im       - integer, horiz dimention and num of used pts      1    !
!     km       - integer, vertical soil layer dimension            1    !
!     ps       - real, surface pressure (pa)                       im   !
!     t1       - real, surface layer mean temperature (k)          im   !
!     q1       - real, surface layer mean specific humidity        im   !
!     soiltyp  - integer, soil type (integer index)                im   !
!     vegtype  - integer, vegetation type (integer index)          im   !
!     sigmaf   - real, areal fractional cover of green vegetation  im   !
!     sfcemis  - real, sfc lw emissivity ( fraction )              im   !
!     dlwflx   - real, total sky sfc downward lw flux ( w/m**2 )   im   !
!     dswflx   - real, total sky sfc downward sw flux ( w/m**2 )   im   !
!     snet     - real, total sky sfc netsw flx into ground(w/m**2) im   !
!     delt     - real, time interval (second)                      1    !
!     tg3      - real, deep soil temperature (k)                   im   !
!     cm       - real, surface exchange coeff for momentum (m/s)   im   !
!     ch       - real, surface exchange coeff heat & moisture(m/s) im   !
!     prsl1    - real, sfc layer 1 mean pressure (pa)              im   !
!     prslki   - real,                                             im   !
!     zf       - real, height of bottom layer (m)                  im   !
!     land     - logical, = T if a point with any land             im   !
!     wind     - real, wind speed (m/s)                            im   !
!     slopetyp - integer, class of sfc slope (integer index)       im   !
!     shdmin   - real, min fractional coverage of green veg        im   !
!     shdmax   - real, max fractnl cover of green veg (not used)   im   !
!     snoalb   - real, upper bound on max albedo over deep snow    im   !
!     sfalb    - real, mean sfc diffused sw albedo (fractional)    im   !
!     flag_iter- logical,                                          im   !
!     flag_guess-logical,                                          im   !
!     lheatstrg- logical, flag for canopy heat storage             1    !
!                         parameterization                              !
!     isot     - integer, sfc soil type data source zobler or statsgo   !
!     ivegsrc  - integer, sfc veg type data source umd or igbp          !
!                                                                       !
!  input/outputs:                                                       !
!     weasd    - real, water equivalent accumulated snow depth (mm) im  !
!     snwdph   - real, snow depth (water equiv) over land          im   !
!     tskin    - real, ground surface skin temperature ( k )       im   !
!     tprcp    - real, total precipitation                         im   !
!     srflag   - real, snow/rain flag for precipitation            im   !
!     smc      - real, total soil moisture content (fractional)   im,km !
!     stc      - real, soil temp (k)                              im,km !
!     slc      - real, liquid soil moisture                       im,km !
!     canopy   - real, canopy moisture content (m)                 im   !
!     trans    - real, total plant transpiration (m/s)             im   !
!     tsurf    - real, surface skin temperature (after iteration)  im   !
!     zorl     - real, surface roughness                           im   !
!                                                                       !
!  outputs:                                                             !
!     sncovr1  - real, snow cover over land (fractional)           im   !
!     qsurf    - real, specific humidity at sfc                    im   !
!     gflux    - real, soil heat flux (w/m**2)                     im   !
!     drain    - real, subsurface runoff (mm/s)                    im   !
!     evap     - real, evaperation from latent heat flux           im   !
!     hflx     - real, sensible heat flux                          im   !
!     ep       - real, potential evaporation                       im   !
!     runoff   - real, surface runoff (m/s)                        im   !
!     cmm      - real,                                             im   !
!     chh      - real,                                             im   !
!     evbs     - real, direct soil evaporation (m/s)               im   !
!     evcw     - real, canopy water evaporation (m/s)              im   !
!     sbsno    - real, sublimation/deposit from snopack (m/s)      im   !
!     snowc    - real, fractional snow cover                       im   !
!     stm      - real, total soil column moisture content (m)      im   !
!     snohf    - real, snow/freezing-rain latent heat flux (w/m**2)im   !
!     smcwlt2  - real, dry soil moisture threshold                 im   !
!     smcref2  - real, soil moisture threshold                     im   !
!     wet1     - real, normalized soil wetness                     im   !
!                                                                       !
!  ====================    end of description    =====================  !

!-----------------------------------
      subroutine sfc_drv                                                &
!...................................
!  ---  inputs:
     &     ( im, km, ps, t1, q1, soiltyp, vegtype, sigmaf,              &
     &       sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm, ch,          &
     &       prsl1, prslki, zf, land, wind, slopetyp,                   &
     &       shdmin, shdmax, snoalb, sfalb, flag_iter, flag_guess,      &
     &       lheatstrg, isot, ivegsrc,                                  &
     &       bexppert, xlaipert, vegfpert,pertvegf,                     &  ! sfc perts, mgehne
!  ---  in/outs:
     &       weasd, snwdph, tskin, tprcp, srflag, smc, stc, slc,        &
     &       canopy, trans, tsurf, zorl,                                &
!  ---  outputs:
     &       sncovr1, qsurf, gflux, drain, evap, hflx, ep, runoff,      &
     &       cmm, chh, evbs, evcw, sbsno, snowc, stm, snohf,            &
     &       smcwlt2, smcref2, wet1                                     &
     &     )
!
      use machine , only : kind_phys
      use funcphys, only : fpvs
      use physcons, only : grav   => con_g,    cp   => con_cp,          &
     &                     hvap   => con_hvap, rd   => con_rd,          &
     &                     eps    => con_eps, epsm1 => con_epsm1,       &
     &                     rvrdm1 => con_fvirt

      use surface_perturbation, only : ppfbet

      implicit none

!  ---  constant parameters:
      real(kind=kind_phys), parameter :: cpinv   = 1.0/cp
      real(kind=kind_phys), parameter :: hvapi   = 1.0/hvap
      real(kind=kind_phys), parameter :: elocp   = hvap/cp
      real(kind=kind_phys), parameter :: rhoh2o  = 1000.0
      real(kind=kind_phys), parameter :: a2      = 17.2693882
      real(kind=kind_phys), parameter :: a3      = 273.16
      real(kind=kind_phys), parameter :: a4      = 35.86
      real(kind=kind_phys), parameter :: a23m4   = a2*(a3-a4)

      real(kind=kind_phys), save         :: zsoil_noah(4)
      data zsoil_noah / -0.1, -0.4, -1.0, -2.0 /

!  ---  input:
      integer, intent(in) :: im, km, isot, ivegsrc
      real (kind=kind_phys), dimension(5), intent(in) :: pertvegf

      integer, dimension(im), intent(in) :: soiltyp, vegtype, slopetyp

      real (kind=kind_phys), dimension(im), intent(in) :: ps,           &
     &       t1, q1, sigmaf, sfcemis, dlwflx, dswsfc, snet, tg3, cm,    &
     &       ch, prsl1, prslki, wind, shdmin, shdmax,                   &
     &       snoalb, sfalb, zf,
     &       bexppert, xlaipert, vegfpert

      real (kind=kind_phys),  intent(in) :: delt

      logical, dimension(im), intent(in) :: flag_iter, flag_guess, land

      logical, intent(in) :: lheatstrg

!  ---  in/out:
      real (kind=kind_phys), dimension(im), intent(inout) :: weasd,     &
     &       snwdph, tskin, tprcp, srflag, canopy, trans, tsurf, zorl

      real (kind=kind_phys), dimension(im,km), intent(inout) ::         &
     &       smc, stc, slc

!  ---  output:
      real (kind=kind_phys), dimension(im), intent(out) :: sncovr1,     &
     &       qsurf, gflux, drain, evap, hflx, ep, runoff, cmm, chh,     &
     &       evbs, evcw, sbsno, snowc, stm, snohf, smcwlt2, smcref2,    &
     &       wet1

!  ---  locals:
      real (kind=kind_phys), dimension(im) :: rch, rho,                 &
     &       q0, qs1, theta1,       weasd_old, snwdph_old,              &
     &       tprcp_old, srflag_old, tskin_old, canopy_old

      real (kind=kind_phys), dimension(km) :: et, sldpth, stsoil,       &
     &       smsoil, slsoil

      real (kind=kind_phys), dimension(im,km) :: zsoil, smc_old,        &
     &       stc_old, slc_old

      real (kind=kind_phys) :: alb, albedo, beta, chx, cmx, cmc,        &
     &       dew, drip, dqsdt2, ec, edir, ett, eta, esnow, etp,         &
     &       flx1, flx2, flx3, ffrozp, lwdn, pc, prcp, ptu, q2,         &
     &       q2sat, solnet, rc, rcs, rct, rcq, rcsoil, rsmin,           &
     &       runoff1, runoff2, runoff3, sfcspd, sfcprs, sfctmp,         &
     &       sfcems, sheat, shdfac, shdmin1d, shdmax1d, smcwlt,         &
     &       smcdry, smcref, smcmax, sneqv, snoalb1d, snowh,            &
     &       snomlt, sncovr, soilw, soilm, ssoil, tsea, th2, tbot,      &
     &       xlai, zlvl, swdn, tem, z0, bexpp, xlaip, vegfp,            &
     &       mv,sv,alphav,betav,vegftmp

      integer :: couple, ice, nsoil, nroot, slope, stype, vtype
      integer :: i, k, iflag

!
!===> ...  begin here
!
!  --- ...  save land-related prognostic fields for guess run

      do i = 1, im
        if (land(i) .and. flag_guess(i)) then
          weasd_old(i)  = weasd(i)
          snwdph_old(i) = snwdph(i)
          tskin_old(i)  = tskin(i)
          canopy_old(i) = canopy(i)
          tprcp_old(i)  = tprcp(i)
          srflag_old(i) = srflag(i)

          do k = 1, km
            smc_old(i,k) = smc(i,k)
            stc_old(i,k) = stc(i,k)
            slc_old(i,k) = slc(i,k)
          enddo
        endif   ! land & flag_guess
      enddo

!  --- ...  initialization block

      do i = 1, im
        if (flag_iter(i) .and. land(i)) then
          ep(i)     = 0.0
          evap (i)  = 0.0
          hflx (i)  = 0.0
          gflux(i)  = 0.0
          drain(i)  = 0.0
          canopy(i) = max(canopy(i), 0.0)

          evbs (i)  = 0.0
          evcw (i)  = 0.0
          trans(i)  = 0.0
          sbsno(i)  = 0.0
          snowc(i)  = 0.0
          snohf(i)  = 0.0
        endif   ! flag_iter & land
      enddo

!  --- ...  initialize variables

      do i = 1, im
        if (flag_iter(i) .and. land(i)) then
          q0(i)   = max(q1(i), 1.e-8)   !* q1=specific humidity at level 1 (kg/kg)
          theta1(i) = t1(i) * prslki(i) !* adiabatic temp at level 1 (k)

          rho(i) = prsl1(i) / (rd*t1(i)*(1.0+rvrdm1*q0(i)))
          qs1(i) = fpvs( t1(i) )        !* qs1=sat. humidity at level 1 (kg/kg)
          qs1(i) = max(eps*qs1(i) / (prsl1(i)+epsm1*qs1(i)), 1.e-8)
          q0 (i) = min(qs1(i), q0(i))
        endif   ! flag_iter & land
      enddo

      do i = 1, im
        if (flag_iter(i) .and. land(i)) then
          do k = 1, km
            zsoil(i,k) = zsoil_noah(k)
          enddo
        endif   ! flag_iter & land
      enddo

      do i = 1, im
        if (flag_iter(i) .and. land(i)) then

!  --- ...  noah: prepare variables to run noah lsm
!   1. configuration information (c):
!      ------------------------------
!    couple  - couple-uncouple flag (=1: coupled, =0: uncoupled)
!    ffrozp  - flag for snow-rain detection (1.=all snow, 0.=all rain, 0-1 mixed)
!    ice     - sea-ice flag (=1: sea-ice, =0: land)
!    dt      - timestep (sec) (dt should not exceed 3600 secs) = delt
!    zlvl    - height (m) above ground of atmospheric forcing variables
!    nsoil   - number of soil layers (at least 2)
!    sldpth  - the thickness of each soil layer (m)

          couple = 1                      ! run noah lsm in 'couple' mode
! use srflag directly to allow fractional rain/snow
!          if     (srflag(i) == 1.0) then  ! snow phase
!            ffrozp = 1.0
!          elseif (srflag(i) == 0.0) then  ! rain phase
!            ffrozp = 0.0
!          endif
          ffrozp = srflag(i)
          ice = 0

          zlvl = zf(i)

          nsoil = km
          sldpth(1) = - zsoil(i,1)
          do k = 2, km
            sldpth(k) = zsoil(i,k-1) - zsoil(i,k)
          enddo

!   2. forcing data (f):
!      -----------------
!    lwdn    - lw dw radiation flux (w/m2)
!    solnet  - net sw radiation flux (dn-up) (w/m2)
!    sfcprs  - pressure at height zlvl above ground (pascals)
!    prcp    - precip rate (kg m-2 s-1)
!    sfctmp  - air temperature (k) at height zlvl above ground
!    th2     - air potential temperature (k) at height zlvl above ground
!    q2      - mixing ratio at height zlvl above ground (kg kg-1)

          lwdn   = dlwflx(i)         !..downward lw flux at sfc in w/m2
          swdn   = dswsfc(i)         !..downward sw flux at sfc in w/m2
          solnet = snet(i)           !..net sw rad flx (dn-up) at sfc in w/m2
          sfcems = sfcemis(i)

          sfcprs = prsl1(i) 
          prcp   = rhoh2o * tprcp(i) / delt
          sfctmp = t1(i)  
          th2    = theta1(i)
          q2     = q0(i)

!   3. other forcing (input) data (i):
!      ------------------------------
!    sfcspd  - wind speed (m s-1) at height zlvl above ground
!    q2sat   - sat mixing ratio at height zlvl above ground (kg kg-1)
!    dqsdt2  - slope of sat specific humidity curve at t=sfctmp (kg kg-1 k-1)

          sfcspd = wind(i)
          q2sat  =  qs1(i)
          dqsdt2 = q2sat * a23m4/(sfctmp-a4)**2

!   4. canopy/soil characteristics (s):
!      --------------------------------
!    vegtyp  - vegetation type (integer index)                       -> vtype
!    soiltyp - soil type (integer index)                             -> stype
!    slopetyp- class of sfc slope (integer index)                    -> slope
!    shdfac  - areal fractional coverage of green vegetation (0.0-1.0)
!    shdmin  - minimum areal fractional coverage of green vegetation -> shdmin1d
!    ptu     - photo thermal unit (plant phenology for annuals/crops)
!    alb     - backround snow-free surface albedo (fraction)
!    snoalb  - upper bound on maximum albedo over deep snow          -> snoalb1d
!    tbot    - bottom soil temperature (local yearly-mean sfc air temp)

          vtype  = vegtype(i)
          stype  = soiltyp(i)
          slope  = slopetyp(i)
          shdfac = sigmaf(i)

!  perturb vegetation fraction that goes into sflx, use the same
!  perturbation strategy as for albedo (percentile matching)
          vegfp  = vegfpert(i)                    ! sfc-perts, mgehne
          if (pertvegf(1) > 0.0) then
                ! compute beta distribution parameters for vegetation fraction
                mv = shdfac
                sv = pertvegf(1)*mv*(1.-mv)
                alphav = mv*mv*(1.0-mv)/(sv*sv)-mv
                betav  = alphav*(1.0-mv)/mv
! compute beta distribution value corresponding
! to the given percentile albPpert to use as new albedo
                call ppfbet(vegfp,alphav,betav,iflag,vegftmp)
                shdfac = vegftmp
          endif
! *** sfc-perts, mgehne

          shdmin1d = shdmin(i)
          shdmax1d = shdmax(i)
          snoalb1d = snoalb(i)

          ptu  = 0.0
          alb  = sfalb(i)
          tbot = tg3(i)

!   5. history (state) variables (h):
!      ------------------------------
!    cmc     - canopy moisture content (m)
!    t1      - ground/canopy/snowpack) effective skin temperature (k)   -> tsea
!    stc(nsoil) - soil temp (k)                                         -> stsoil
!    smc(nsoil) - total soil moisture content (volumetric fraction)     -> smsoil
!    sh2o(nsoil)- unfrozen soil moisture content (volumetric fraction)  -> slsoil
!    snowh   - actual snow depth (m)
!    sneqv   - liquid water-equivalent snow depth (m)
!    albedo  - surface albedo including snow effect (unitless fraction)
!    ch      - surface exchange coefficient for heat and moisture (m s-1) -> chx
!    cm      - surface exchange coefficient for momentum (m s-1)          -> cmx

          cmc  = canopy(i) * 0.001           ! convert from mm to m
          tsea = tsurf(i)                    ! clu_q2m_iter

          do k = 1, km
            stsoil(k) = stc(i,k)
            smsoil(k) = smc(i,k)
            slsoil(k) = slc(i,k)
          enddo

          snowh = snwdph(i) * 0.001         ! convert from mm to m
          sneqv = weasd(i)  * 0.001         ! convert from mm to m
          if (sneqv /= 0.0 .and. snowh == 0.0) then
            snowh = 10.0 * sneqv
          endif

          chx    = ch(i)  * wind(i)              ! compute conductance
          cmx    = cm(i)  * wind(i)
          chh(i) = chx * rho(i)
          cmm(i) = cmx

!  ---- ... outside sflx, roughness uses cm as unit
          z0 = zorl(i)/100.
!  ---- mgehne, sfc-perts
          bexpp  = bexppert(i)                   ! sfc perts, mgehne
          xlaip  = xlaipert(i)                   ! sfc perts, mgehne

!  --- ...  call noah lsm

          call sflx                                                     &
!  ---  inputs:
     &     ( nsoil, couple, ice, ffrozp, delt, zlvl, sldpth,            &
     &       swdn, solnet, lwdn, sfcems, sfcprs, sfctmp,                &
     &       sfcspd, prcp, q2, q2sat, dqsdt2, th2, ivegsrc,             &
     &       vtype, stype, slope, shdmin1d, alb, snoalb1d,              &
     &       bexpp, xlaip,                                              & ! sfc-perts, mgehne
     &       lheatstrg,                                                 &
!  ---  input/outputs:
     &       tbot, cmc, tsea, stsoil, smsoil, slsoil, sneqv, chx, cmx,  &
     &       z0,                                                        &
!  ---  outputs:
     &       nroot, shdfac, snowh, albedo, eta, sheat, ec,              &
     &       edir, et, ett, esnow, drip, dew, beta, etp, ssoil,         &
     &       flx1, flx2, flx3, runoff1, runoff2, runoff3,               &
     &       snomlt, sncovr, rc, pc, rsmin, xlai, rcs, rct, rcq,        &
     &       rcsoil, soilw, soilm, smcwlt, smcdry, smcref, smcmax)

!  --- ...  noah: prepare variables for return to parent mode
!   6. output (o):
!      -----------
!    eta     - actual latent heat flux (w m-2: positive, if upward from sfc)
!    sheat   - sensible heat flux (w m-2: positive, if upward from sfc)
!    beta    - ratio of actual/potential evap (dimensionless)
!    etp     - potential evaporation (w m-2)
!    ssoil   - soil heat flux (w m-2: negative if downward from surface)
!    runoff1 - surface runoff (m s-1), not infiltrating the surface
!    runoff2 - subsurface runoff (m s-1), drainage out bottom

          evap(i)  = eta
          hflx(i)  = sheat
          gflux(i) = ssoil

          evbs(i)  = edir
          evcw(i)  = ec
          trans(i) = ett
          sbsno(i) = esnow
          snowc(i) = sncovr
          stm(i)   = soilm * 1000.0 ! unit conversion (from m to kg m-2)
          snohf(i) = flx1 + flx2 + flx3

          smcwlt2(i) = smcwlt
          smcref2(i) = smcref

          ep(i)      = etp
          tsurf(i)   = tsea

          do k = 1, km
            stc(i,k) = stsoil(k) 
            smc(i,k) = smsoil(k)
            slc(i,k) = slsoil(k)
          enddo
          wet1(i) = smsoil(1) / smcmax !Sarah Lu added 09/09/2010 (for GOCART)

!  --- ...  unit conversion (from m s-1 to mm s-1 and kg m-2 s-1)
          runoff(i)  = runoff1 * 1000.0
          drain (i)  = runoff2 * 1000.0

!  --- ...  unit conversion (from m to mm)
          canopy(i)  = cmc   * 1000.0
          snwdph(i)  = snowh * 1000.0
          weasd(i)   = sneqv * 1000.0
          sncovr1(i) = sncovr
!  ---- ... outside sflx, roughness uses cm as unit (update after snow's
!  effect)
          zorl(i) = z0*100.

!  --- ...  do not return the following output fields to parent model
!    ec      - canopy water evaporation (m s-1)
!    edir    - direct soil evaporation (m s-1)
!    et(nsoil)-plant transpiration from a particular root layer (m s-1)
!    ett     - total plant transpiration (m s-1)
!    esnow   - sublimation from (or deposition to if <0) snowpack (m s-1)
!    drip    - through-fall of precip and/or dew in excess of canopy
!              water-holding capacity (m)
!    dew     - dewfall (or frostfall for t<273.15) (m)
!    beta    - ratio of actual/potential evap (dimensionless)
!    flx1    - precip-snow sfc (w m-2)
!    flx2    - freezing rain latent heat flux (w m-2)
!    flx3    - phase-change heat flux from snowmelt (w m-2)
!    snomlt  - snow melt (m) (water equivalent)
!    sncovr  - fractional snow cover (unitless fraction, 0-1)
!    runoff3 - numerical trunctation in excess of porosity (smcmax)
!              for a given soil layer at the end of a time step
!    rc      - canopy resistance (s m-1)
!    pc      - plant coefficient (unitless fraction, 0-1) where pc*etp
!              = actual transp
!    xlai    - leaf area index (dimensionless)
!    rsmin   - minimum canopy resistance (s m-1)
!    rcs     - incoming solar rc factor (dimensionless)
!    rct     - air temperature rc factor (dimensionless)
!    rcq     - atmos vapor pressure deficit rc factor (dimensionless)
!    rcsoil  - soil moisture rc factor (dimensionless)
!    soilw   - available soil moisture in root zone (unitless fraction
!              between smcwlt and smcmax)
!    soilm   - total soil column moisture content (frozen+unfrozen) (m)
!    smcwlt  - wilting point (volumetric)
!    smcdry  - dry soil moisture threshold where direct evap frm top
!              layer ends (volumetric)
!    smcref  - soil moisture threshold where transpiration begins to
!              stress (volumetric)
!    smcmax  - porosity, i.e. saturated value of soil moisture
!              (volumetric)
!    nroot   - number of root layers, a function of veg type, determined
!              in subroutine redprm.

        endif   ! flag_iter and flag
      enddo   ! end do_i_loop

!   --- ...  compute qsurf (specific humidity at sfc)

      do i = 1, im
        if (flag_iter(i) .and. land(i)) then
          rch(i)   = rho(i) * cp * ch(i) * wind(i)
          qsurf(i) = q1(i)  + evap(i) / (elocp * rch(i))
        endif   ! flag_iter & flag
      enddo

      do i = 1, im
        if (flag_iter(i) .and. land(i)) then
          tem     = 1.0 / rho(i)
          hflx(i) = hflx(i) * tem * cpinv
          evap(i) = evap(i) * tem * hvapi
        endif   ! flag_iter & flag
      enddo

!  --- ...  restore land-related prognostic fields for guess run

      do i = 1, im
        if (land(i)) then
          if (flag_guess(i)) then
            weasd(i)  = weasd_old(i)
            snwdph(i) = snwdph_old(i)
            tskin(i)  = tskin_old(i)
            canopy(i) = canopy_old(i)
            tprcp(i)  = tprcp_old(i)
            srflag(i) = srflag_old(i)

            do k = 1, km
              smc(i,k) = smc_old(i,k)
              stc(i,k) = stc_old(i,k)
              slc(i,k) = slc_old(i,k)
            enddo
          else    ! flag_guess = F
            tskin(i) = tsurf(i)
          endif   ! flag_guess
        endif     ! flag

      enddo
!
      return
!...................................
      end subroutine sfc_drv
!-----------------------------------
!> @}
!> @}
       end module module_sfc_drv
