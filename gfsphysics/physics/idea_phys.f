
! WAM-IPE Physics Driver
!
      subroutine idea_phys(im,ix,levs,prsi,prsl,                        
     &                     adu,adv,adt,adr,ntrac,dtp,lat,               
     &                     solhr,slag,sdec,cdec,sinlat,coslat,          
     &                     xlon,xlat,oro,cozen,swh,hlw,dt6dt,           
     &                 thermodyn_id,sfcpress_id,gen_coord_hybrid,mpi_id,
     &                     fhour, kstep)
!    &                     mpi_ior,mpi_comm, fhour, kstep,
!    &                     gzmt, gmmt, gjhr, gshr, go2dr)
!-----------------------------------------------------------------------
!SK 2019Mar08
! -- Unused arguments mpi_ior, mpi_comm, gzmt, gmmt, gjhr, gshr, go2dr
! -- are released from the argument list. Arrays gzmt, gmmt, gjhr, gshr,
! -- and go2dr continue to be local arrays (Original code).
!-----------------------------------------------------------------------
! add temp, wind changes due to viscosity and thermal conductivity
! also solar heating
! Apr 06 2012   Henry Juang, initial implement for NEMS
! Jul 26 2012   Jun Wang, add mpi info
! Sep 06 2012   Jun Wang, add changing pressure to cb
! Dec    2012   Jun Wang, change to new rad_merge (from Rashid and Fei)
! May    2013   Jun Wang, tmp updated after rad_merge
! Jun    2013   S. Moorthi Some optimization and cosmetic changes
! Oct    2013   Henry Juang, correct the sequence to get prsi from model top
! Sep    2017   Weiyu Yang, add IPE back coupling to WAM code.
!=========================================================================================
!
! 2016/17     Upgrades of WAM physics
!
! Feb    2017   VAY, Manual updates of NEMS_WAM201606 =>NEMS_WAM201612 with
!               WAM physics, major updates
!
!               0) restructure of idea_phys.f (~600 lines) => (~200 lines)
!                  as  single sub-ne,  WAM-physics interface with NEMS-physics
!               1) cb => Pa (SI units across all WAM physics)
!               2) idea_tracer as suggested by TFR with Jo2 of TIME-GCM
!               3) idea_ion  with TIROS and correction of RAA 1/dip *dip
!               4) Four options to manipulate with SOLAR-GEO drivers
!               5) Fixing SNOE-NO and NO cooling
!               6) New  SHEAT with O2-photolysis and SRB inside => OLD in idea_solar_2014.f
!               7) Updating input sola flux paramterization for SMIN at ~F107 = 70
!               8) Updates in idea_o2_o3.f w/o SRB
!               9) Radiation merge boundaries for SHEAT with SRB+SRC+LYA below ~50 km
!              10) WAM-calendar to advance "time-dependent" input drivers.
!              11) Various minor errors that can "crash" WAM physics, like (mismatches of IM/IX)
!              12) First 120 steps (6-hr) run with MSIS composition 
!                  if  integer, parameter  :: irlx_msis=1, ; 0 - IC-WAM
!
! Not included from
! Research versions: a) Unified GW Physics; b) IAU in WAM-physics;
!                    c) Diurnal O3 cycles ; d) Nudging/Assimilating the "MSIS"-composition
!=========================================================================================
!
!
!
      use idea_wam_control, only : SPW_DRIVERS
!
      use idea_composition, only : nlev_h2o,nlevc_h2o, nlev_co2
      use idea_composition, only : mmr_min, amo, amo2, amo3, amn2
      use wam_ion,          only : idea_ion
!     use wam_f107_kp_mod,  only : f107_wy, kp_wy, kdt_interval, interpolate_weight, kpa_wy, f107d_wy, hp_wy, hpi_wy 
      use wam_f107_kp_mod,  only : f107_wy, kp_wy, kdt_interval,
     &                     interpolate_weight, kpa_wy, f107d_wy, nhp_wy,
     &                             shpi_wy, shp_wy, nhpi_wy, 
     &                   swbz_wy, swvel_wy, swbt_wy, swang_wy, swden_wy 

!     Changed by Zhuxiao.Li(05/2017) for back to the path to read in F10.7 and Kp
!     from solar_in namelist instead of from wam_f107_kp_mod.f
!     use  wam_f107_kp_mod, only : f107_fix, f107d_fix, kp_fix
      use  idea_solar_input,only : f107_fix, f107a_fix, kp_fix 
!
      use  wam_f107_kp_mod, only : fix_spweather_data
!     use  wam_f107_kp_mod, only : swpcf107_fix, swpcf107d_fix, swpckp_fix
!
!
!SK   use IDEA_MPI_def,      only : mpi_WAM_quit, mpi_id
!SK   use IDEA_MPI_def,      only : mpi_err, MPI_COMM_ALL, info      ! or use "me mpi_or mpi_comm

      use wam_date_calendar, only : idat_wam, irhour_wam
      use wam_date_calendar, only : CURRENT_NCEP_JDAT, ndwam
!  
!vay 0/2015  itheia = 1 no mpi-broadcast from mpi_id =0
!
      use idea_solar_input,  only : itheia

      use idea_solar_input,  only : solar_wamstep_advance
      use idea_solar_input,  only : solar_waccmx_advance
      use idea_solar_input,  only : idea_solar_fix
      use idea_solar_input,  only : wf107_s,  wkp_s
      use idea_solar_input,  only : wf107a_s, wap_s
      use idea_solar_input,  only : weuv_s, nwafix  
!      use idea_solar_input,  only : f107_fix, kp_fix, f107a_fix
! 
      use idea_imf_input,    only : idea_imf_fix, imf_wamstep_advance
      use idea_imf_input,    only : Bz_s, By_s, Swden_s, Swvel_s    ! only for debugging & comment later
      use idea_imf_input,    only : w_ndx1_imf, w_ndx2_imf          ! inly for debugging
      use module_IPE_to_WAM, only : lowst_ipe_level, ipe_to_wam_coupling
!
!      use GW_unified,  only: IMPL_UNIF_GW    ! 'GFS_91L' or 'WAM150L'
!
      implicit none
! Argument
      integer, intent(in) :: im              ! number of data points in adt (first dim)
      integer, intent(in) :: ix              ! max data points in adt (first dim)
      integer, intent(in) :: levs            ! number of pressure levels
      integer, intent(in) :: lat             ! latitude index
      integer, intent(in) :: ntrac           ! number of tracer
      integer, intent(in) :: mpi_id          ! my pe
!sk2019Mar08_B2
!     integer, intent(in) :: mpi_ior         ! mpi real for io
!     integer, intent(in) :: mpi_comm        ! mpi communicator
!sk2019Mar08_E2
!
      integer, intent(in) :: kstep           ! # of model steps from gloopb.f 
      real   , intent(in) :: fhour           ! # model time in hours since fhini in...gloopb.f 
      real,    intent(in) :: dtp             ! time step in second     
!
      real, intent(in)    :: prsi(ix,levs+1) ! pressure
      real, intent(in)    :: prsl(ix,levs)   ! pressure
!SK2020Oct10
!     real, intent(in)    :: hlw(ix,levs)    ! long wave rad (K/s)
!     real, intent(in)    :: swh(ix,levs)    ! short wave rad (K/s)
      real, intent(inout) :: hlw(ix,levs)    ! long wave rad (K/s)
      real, intent(inout) :: swh(ix,levs)    ! short wave rad (K/s)
!
      real, intent(in)    :: solhr,slag,sdec,cdec ! for solar zenith angle
      real, intent(in)    :: cozen(im)       ! daytime avg cos zenith angle see radiation_astronomy.f
      real, intent(in)    :: oro(im)         ! surface height (m)
      real, intent(in)    :: xlon(im)        ! (0-360)/180.*PI, radians
      real, intent(in)    :: xlat(im)        ! (90 => -90)/180.*PI, radians
      real, intent(in)    :: coslat(im)
      real, intent(in)    :: sinlat(im)
      real, intent(inout) :: adr(ix,levs,ntrac) ! tracer
      real, intent(inout) :: adt(ix,levs)       ! temperature
      real, intent(inout) :: adu(ix,levs)       ! W-E u
      real, intent(inout) :: adv(ix,levs)       ! S-N v
!
      real, intent(inout) :: dt6dt(ix,levs,6)   ! diagnostic 3D-array ....never used 
!                                               !
! (2)-wtot-merged (SH-LW, H2O, CO2, O2,O3), (4-5-6) for Strobel +Cooling
! (1)-MT_SHEAT(EUV+?),  (3) - Joule heating
!
      integer,intent(in)  :: thermodyn_id, sfcpress_id
      logical,intent(in)  :: gen_coord_hybrid
!
! Local variables
!      real,parameter      :: pa2cb=0.001,cb2pa=1000.
!
      integer, parameter  :: ntrac_i=2                  ! number of 2 WAM chem. tracers (O-O2)
!
!     real    :: f107_curdt, f107d_curdt, kp_curdt, f107a_fix    
      real    :: f107_curdt, f107d_curdt, kp_curdt, kpa_curdt, 
     &           nhp_curdt, nhpi_curdt
      real    :: swbz_curdt, swvel_curdt, swbt_curdt, swang_curdt, 
     &           shp_curdt, shpi_curdt
      real    :: swden_curdt
      integer :: Mjdat(ndwam)                           ! IDAT_WAM + FHOUR
      real    :: Hcur                                   !  current hour+min+sec real 
!
! locals 2D (ix, levs) -rule  1D can be only (im)   
  

      real  :: cospass(im), xmu(im)
!================================================
! Composition related locals
!
      real  :: nair(ix,levs)
      real ::  rho(ix,levs) 
      real  :: am(ix,levs)
!
! VAY-201701, options for the MSIS-constrained thermosphere
!             and O3-diurnal variations 
!msis
!       integer, parameter  :: irlx_msis=1
!       real, dimension(ix, levs)  :: T_MSIS, O3P_MSIS
!       real, dimension(ix, levs)  :: O2_MSIS, N2_MSIS
!       real                       :: Nair_msis
!
      real  :: am29(ix,levs), amu_msis
      real  :: o3_n(ix,levs)    ! transported ozone with updates of 
                                ! oxygen chemistry with diurnal variations of O2, O3, O
      real  :: o3_ng(ix,levs)   ! ozone_gsm(1:k71) + o3glob(k71+1:levs)  

!
      real  :: o_n(ix,levs),o2_n(ix,levs),n2_n(ix,levs)    
!================================================
! Tendency-related locals
      real  :: dudt(ix,levs),dvdt(ix,levs),dtdt(ix,levs)
      real, dimension(ix,lowst_ipe_level:levs)  ::
     &                            gzmt, gmmt, gjhr, gshr, go2dr
!SK2020Sep6
!     real  :: adu_(ix,levs),adv_(ix,levs),adt_(ix,levs),
!    &         adr_(ix,levs,ntrac)
!
!================================================
! Radiance-related locals
!
      real  :: dtRAD(ix,levs), dtco2c(ix,levs),dtco2h(ix,levs), 
     &         dth2oh(ix,levs),dth2oc(ix,levs),
     &         dto3(ix,levs),  wtot(ix, levs)  
!================================================
! Thermodynamics-related locals
       real ::  cp(ix,levs)                
       real ::  grav(ix,levs) ,zg(ix,levs)     
       real ::  prslk(ix,levs), prsik(ix,levs+1), exner(ix,levs) !
       real ::  phil(ix,levs),   phii(ix,levs+1)
       real ::  rdelp(ix, levs)                   ! k=1, plevs ..rdelp(i,k) = 1./(PRSI(i,k) - PRSI(i,k+1))
                                                  !
       real ::  delPa(ix, levs)                   ! thickness in Pa
!================================================
! geo-solar-related vars
       real  :: utsec,   sda
!
       real  :: maglat(im),maglon(im),btot(im),                    
     &          dipang(im),essa(im)
!
       integer  :: dayno                           ! ddd of year
       integer i,k,l, j1,j2

       real, dimension(ix)    :: xpk_low, xpk_high
       integer, dimension(ix) :: plow, phigh
!================================================
!sk2019Mar08_B1
!      logical, parameter :: return_from_here = .true.
!      logical, parameter :: skprnt = .true.
       logical, save :: skprnt
       real :: sk1,sk2,sk3
!
       skprnt = mpi_id.eq.0
!     if (return_from_here) goto 1111
!sk2019Mar08_E1
!================================================
! SK2020Sep2 Zero out dudt, dvdt, and dtdt to zero at start
! SK2020Sep6 Save adu, adv, adt to debug IDEA_PHYS
      dudt(:,:) = 0.0
      dvdt(:,:) = 0.0
      dtdt(:,:) = 0.0
!     dtRAD(:,:) = 0.0
!     dtco2c(:,:) = 0.0
!     dtco2h(:,:) = 0.0
!     dth2oh(:,:) = 0.0
!     dth2oc(:,:) = 0.0
!     dto3(:,:) = 0.0
!     wtot(:,:) = 0.0
!     do k=1,levs
!       do i=1,im
!         adu_(i,k) = adu(i,k)
!         adv_(i,k) = adv(i,k)
!         adt_(i,k) = adt(i,k)
!       enddo
!     enddo
!     do l=1,ntrac
!     do k=1,levs
!       do i=1,im
!         adr_(i,k,l) = adr(i,k,l)
!       enddo
!     enddo
!     enddo
!
!SK2020Sep21
! Reset incoming adr negative values to mmr_min
!
      do l = 1,ntrac
        do k = 1,levs
          do i = 1,im
            adr(i,k,l) = max(adr(i,k,l),mmr_min)
          enddo
        enddo
      enddo
!sk2020sep20
      sk1 = minval(adr)
      if (sk1 .lt. 0)
     &print*,'sep21:idea_phys:input/min_adr='
     &  ,sk1,' me= ',mpi_id
!sk
!================================================
!
! Start with Space Wea /Climate Drivers
!
!===============================================
! option 1:  SPW_DRIVERS ='swpc_fst'
!
       if (trim(SPW_DRIVERS)=='swpc_fst') then
          f107_curdt  = f107_wy (kdt_interval) * interpolate_weight  + 
     &              f107_wy (kdt_interval+1) * (1-interpolate_weight)
          kp_curdt    = kp_wy   (kdt_interval) * interpolate_weight  + 
     &              kp_wy   (kdt_interval+1) * (1-interpolate_weight)
          f107d_curdt = f107d_wy(kdt_interval) * interpolate_weight  + 
     &              f107d_wy(kdt_interval+1) * (1-interpolate_weight)
          kpa_curdt   = kpa_wy  (kdt_interval) * interpolate_weight  + 
     &              kpa_wy  (kdt_interval+1) * (1-interpolate_weight)
          nhp_curdt   = nhp_wy  (kdt_interval) * interpolate_weight  + 
     &              nhp_wy  (kdt_interval+1) * (1-interpolate_weight)
          nhpi_curdt  = nhpi_wy (kdt_interval) * interpolate_weight  + 
     &              nhpi_wy (kdt_interval+1) * (1-interpolate_weight)
          shp_curdt   = shp_wy  (kdt_interval) * interpolate_weight  + 
     &              shp_wy  (kdt_interval+1) * (1-interpolate_weight)
          shpi_curdt  = shpi_wy (kdt_interval) * interpolate_weight  + 
     &              shpi_wy (kdt_interval+1) * (1-interpolate_weight)
          swbt_curdt  = swbt_wy (kdt_interval) * interpolate_weight  + 
     &              swbt_wy (kdt_interval+1) * (1-interpolate_weight)
          swang_curdt = swang_wy(kdt_interval) * interpolate_weight  + 
     &              swang_wy(kdt_interval+1) * (1-interpolate_weight)
          swvel_curdt = swvel_wy(kdt_interval) * interpolate_weight  + 
     &              swvel_wy(kdt_interval+1) * (1-interpolate_weight)
          swbz_curdt  = swbz_wy (kdt_interval) * interpolate_weight  + 
     &              swbz_wy (kdt_interval+1) * (1-interpolate_weight)
          swden_curdt = swden_wy(kdt_interval) * interpolate_weight  +
     &                 swden_wy(kdt_interval+1) * (1-interpolate_weight)
       else
          kpa_curdt   = 0.0
          shp_curdt   = 0.0
          shpi_curdt  = 0.0
          nhp_curdt   = 0.0
          nhpi_curdt  = 0.0
          swbt_curdt  = 0.0
          swang_curdt = 0.0
          swvel_curdt = 0.0
          swbz_curdt  = 0.0
          swden_curdt = 0.0
       endif
       if (skprnt) print*,' idea_phys: ok001'
!================================================
! option 2:  SPW_DRIVERS ='sair_wam', only year of 2012
!
       if (trim(SPW_DRIVERS)=='sair_wam') then
!--------------------------------------------------------------------
! Optional Call of WAM solar-geo inputs
!--------------------------------------------------------------------
!
       if (idat_wam(1) == 2012 .and. trim(SPW_DRIVERS)=='sair_wam') then
!========================================================================
!vay-08/2015 Advance in time Solar-geo inputs (1-hr or 3-hr or daily)
!vay-01/2016 Advance in time IMF-inputs (5-min)
!----------------------------------------------------------------------  
!       if (idea_imf_fix.eq.0) then
!           CALL IMF_wamstep_advance(mpi_id, Mjdat, Hcur)
!         if (me == 0 .and. kstep .le. 1) print *, Bz_s, By_s,' Bz_s .... By_s, me= ',me
!       endif
!----------------------------------------------------------------------  
          if (idea_solar_fix.le.1) then
          CALL solar_wamstep_advance(mpi_id, Mjdat, Hcur)
      if (skprnt) print*,
     &           ' in idea_phys: <-- solar_wamstep_advance, me= ',mpi_id
!      
!          if (idea_solar_fix.le.1.and.itheia.eq.0) then 
!           call mpi_bcast(wap_s  , 1,    MPI_REAL8,0, MPI_COMM_ALL,info)
!           call mpi_bcast(wkp_s  , 1,    MPI_REAL8,0, MPI_COMM_ALL,info)
!           call mpi_bcast(wf107_s, 1,    MPI_REAL8,0, MPI_COMM_ALL,info)
!           call mpi_bcast(wf107a_s, 1,   MPI_REAL8,0, MPI_COMM_ALL,info)
!          endif
!          if (idea_solar_fix.eq.0) call mpi_bcast(wEUV_s, nwafix,MPI_REAL8,0, MPI_COMM_ALL,info) 
       ENDIF 
!
       ENDIF    ! idat_wam(1) == 2012 . +  "sair_wam":"f107_s, f107a_s, ap_s, kp_s, euv_s(37)"
!for now fixed
            f107_curdt  = f107_fix
            kp_curdt    = kp_fix
            kpa_curdt   = kp_fix
            f107d_curdt = f107a_fix
!
      ENDIF     !'sair_wam'
      if (skprnt) print*,' idea_phys: ok002'
!===============================================
! option 3:  SPW_DRIVERS ='wam_climate': fixed values
!
       if (trim(SPW_DRIVERS)=='wam_climate') then
!          call fix_spweather_data
          f107_curdt  = f107_fix
          kp_curdt    = kp_fix
          kpa_curdt    = kp_fix
          f107d_curdt = f107a_fix
       endif
!
!===============================================
       if (trim(SPW_DRIVERS)=='cires_wam') then
!
! option 4:  SPW_DRIVERS ='cires_wam_climate'
!
! Real data for F107/Kp [2000-2016]: option from VAY available under request
!         
!         CALL solar_waccmx_advance(mpi_id, Mjdat, Hcur, Kstep)   !wf107_s,  wkp_s
!
!         f107_curdt  = wf107_s
!         kp_curdt    = wkp_s
!         f107d_curdt = wf107a_s
!
! For standard WAM-IPE version, fixed values from "solar_in" namelist
!
!         call fix_spweather_data
          f107_curdt = f107_fix
          kp_curdt   = kp_fix
          kpa_curdt   = kp_fix
          f107d_curdt = f107a_fix
        endif
      if (skprnt) print*,' idea_phys: ok003'
! 
        if (mpi_id == 0 .and. kstep <= 1) then
            print *
            print *, kdt_interval, interpolate_weight
            print *, f107_curdt, f107d_curdt, kp_curdt, kpa_curdt, 
     &               nhp_curdt, nhpi_curdt
            print *
            print *, swbt_curdt, swang_curdt, swvel_curdt, swbz_curdt,
     &               shp_curdt, shpi_curdt, 'f107-kp data'
            print *, 'idea_phys'
!           print *, 'VAY-GW:',trim(IMPL_UNIF_GW)
            print *
            print *, 'ID-phys SPW-drivers option: ', trim(SPW_DRIVERS)
            print *
         endif
      if (skprnt) print*,' idea_phys: ok004'
!SK2020Aug31
!     if (return_from_here) goto 1111
!==================================================
! in all WAM-subs pass: f107_curdt, f107d_curdt, kp_curdt
!    defined above  from 4-cases (CLIM-SWPC-SAIR-CIRES)
!
! change prsi and prsl to centibar from pascal
!  WHO INITIATED this pa2cb-MESS, now everything in SI and PA, 1/m3, kg/m3
!      do k=1,levs
!        do i=1,im
!          prsi(i,k) = prsi(i,k)*pa2cb
!          prsl(i,k) = prsl(i,k)*pa2cb
!        enddo
!      enddo
!      do i=1,im
!        prsi(i,levs+1) = prsi(i,levs+1)*pa2cb
!      enddo
!      call get_phi(im,ix,levs,ntrac,adt,adr,                            
!     &             thermodyn_id, sfcpress_id,                           
!     &             gen_coord_hybrid,                                    
!     &             prsi,prsik,prsl,prslk,phii,phil)
!      call phi2z(im,ix,levs,phil,oro,zg,grav)
! prsi and prsl - mass of layers are constant in the pressure-based system
!        no needs to do second loop for 1-st & back-transfer 
!       ....prsi(i,k)    = prsi(i,k)*cb2pa
!============================================================ delete
!      print *, 'vay-prsl PA ', maxval(prsl), minval(prsl), me, ' PE '
!      print *, 'vay-prsi PA ', maxval(prsi), minval(prsi), me, ' PE '
!      print *, 'VAY-PRSI-PRSL', maxval(prsi), maxval(Prsl)
!
!VAY, Nov 2016: single call for: phii,phil, zg, grav, exner, delpa
!
      if (skprnt) print*,
     &           ' in idea_phys: --> get_exner_phi_zgrav, me= ',mpi_id
      CALL get_exner_phi_zgrav(ix, im, levs, ntrac,
     & adr, adt,
     & prsl,  prsi, rdelp, prsik, prslk,phii,phil, 
     & thermodyn_id, sfcpress_id, gen_coord_hybrid, 
     & oro, zg, grav, exner, delpa) 
      if (skprnt) print*,
     &           ' in idea_phys: <-- get_exner_phi_zgrav, me= ',mpi_id

      if (skprnt) print*,' idea_phys: ok005'
!SK2020Aug31
!     if (return_from_here) goto 1111
!
!      print *,'VAY:Grav-g81',maxval(prsi),maxval(zg),maxval(nint(grav))   
!      call mpi_WAM_quit(me, 'after get_exner_phi_zgrav in idea_phys.f ')
!
!    get presolar-related: cos solar zenith angle (instant)
!
!
      if (skprnt) print*,
     &           ' in idea_phys: --> presolar, me= ',mpi_id
      call presolar(im,ix,solhr,slag,                                   
     &              sinlat,coslat,sdec,cdec,xlon,xlat                   
     &              ,cospass,dayno,utsec,sda                            
     &              ,maglat,maglon,btot,dipang,essa)
      if (skprnt) print*,
     &           ' in idea_phys: <-- presolar, me= ',mpi_id

      if (skprnt) print*,' idea_phys: ok006'
!SK2020Aug31
!     if (return_from_here) goto 1111
! ===============================================================
!     Initialize WAM major tracers by MSIS-00
! Disable for WAM-IPE as rquested by NM & TFR Fen 22/2017 
! Good ICs for O-O2 Composition will be delivered by New Chemistry
!      as expected by WAM & IPE managers
!
!      IF (irlx_msis == 1 .and. kstep.lt.120) then
!
!      CALL WAM_MSIS(ix, im, levs, Zg, Xlat, Xlon,  
!     & dayno, utsec,f107_curdt, f107d_curdt, KP_curdt,
!     & T_MSIS, O2_MSIS, N2_MSIS, O3P_MSIS)
!      
!       do k=1, levs
!        do i=1, im
!        amu_msis =(adr(i,k,4)/amo+adr(i,k,5)/amo2+
!     & (1.-adr(i,k,4)-adr(i,k,5))/amn2)*amo     !amo/amair
!       nair_msis = O3P_MSIS(i,k)+N2_MSIS(i,k)+O2_MSIS(i,k)
!        adr(i,k,4) = .75*O3P_MSIS(i,k)/nair_msis*amu_msis +
!     &      .25*adr(i,k,4) 
!        adr(i,k,5) = .75*O2_MSIS(i,k)/nair_msis*amu_msis*2. +
!     &      .25*adr(i,k,5) 
!          if (adr(i,k,4) .lt. 0.0) adr(i,k,4) = mmr_min
!        enddo
!       enddo
!      ENDIF
!===============================================================
!get composition at layers (1/m3) and rho (kg/m3) or/and
!
!     call idea_tracer(im,ix,levs,ntrac,ntrac_i,
!     &                 grav,prsi,prsl,adt,adr, dtp,      
!     &                 o_n,o2_n,n2_n,nair,rho,am,
!     &                 cospass, dayno, zg, me)
!===============================================================
! Get all needed vertical related parameters for merging the 
! coupling back IPE arrays into the WAM model arrays.
!===============================================================
      IF(ipe_to_wam_coupling) THEN
        call get_vertical_parameters_for_merge(gzmt,  
     &    im, ix, lowst_ipe_level, levs, plow, phigh, 
     &    xpk_low, xpk_high, prsl)
      END IF
!=================================================================
!=> compute rho
!=> mid-layers: [o_n, o2_n, o3_n, n2_n, nair]
!=================================================================
      if (skprnt) print*,
     &           ' in idea_phys: --> idea_tracer, me= ',mpi_id
!sk2020sep20
      sk1 = minval(adr)
      if (sk1 .lt. 0)
     &print*,'sep21:idea_phys>idea_tracer:input/min_adr='
     &  ,sk1,' me= ',mpi_id
!sk
      call idea_tracer(im,ix,levs,ntrac,2,grav,prsi,prsl,adt,adr,       
     &                 dtp,o_n,o2_n,o3_n, n2_n,nair,rho,am, am29,
     & cospass, dayno, zg, f107_curdt, f107d_curdt, mpi_id, go2dr, plow,
     & phigh, xpk_low, xpk_high)
!sk2020sep20
      sk1 = minval(adr)
      if (sk1 .lt. 0)
     &print*,'sep21:idea_phys<idea_tracer:input/min_adr='
     &  ,sk1,' me= ',mpi_id
!sk
      if (skprnt) print*,
     &           ' in idea_phys: <-- idea_tracer, me= ',mpi_id
      if (skprnt) print*,' idea_phys: ok007'
!SK2020Aug31
!     if (return_from_here) goto 1111
!     if ( me == 0) print *, maxval(am29), minval(am29), 'VAY-am29C
!
!
! calculate cp and precompute [1/cp/rho =array] for dT/dt = Q/cp/rho
!=================================================================
      if (skprnt) print*,
     &           ' in idea_phys: --> getcp, me= ',mpi_id
!SK  &           ' in idea_phys: --> getcp_idea, me= ',mpi_id
!SK   call getcp_idea(im,ix,levs,ntrac,adr,cp,                          
!SK  &                thermodyn_id,gen_coord_hybrid)
      call get_cp(im,ix,levs,ntrac,adr,cp)
      if (skprnt) print*,
     &           ' in idea_phys: <-- getcp, me= ',mpi_id
!SK  &           ' in idea_phys: <-- getcp_idea, me= ',mpi_id
      if (skprnt) print*,' idea_phys: ok008'
!SK2020Aug31
!     if (return_from_here) goto 1111
!=================================================================
! dissipation +GW physics/turbulent eddies
! VAY-2016/11  version will be available
!   after additional validation and submission
!   of the paper, Unified GW physics in WAM
! "restart bit-by-bit" moved its develop-nt
!  to NEMS/WAM-201601
!============================================

!      if  (IMPL_UNIF_GW.eq.'WAM150L') then
!        call  idea_gwp_dissip(im,ix,levs,ntrac,grav,prsi,prsl,       
!     &   adu,adv,adt,adr, o_n,o2_n,n2_n,dtp,cp,dt6dt,
!     &   delpa, prslk, exner, phil, phii, xlat, me, kstep) 
            
!      else
!--------------------------------------------------------------------------
!  dissipation of 2013-idea: 
!  only molecular diss. (viscosity + conductivity) no eddy diff. for tracers
!  dissipation of 2013-idea:
      if (skprnt) print*,
     &           ' in idea_phys: --> idea_phys_dissipation, me= ',mpi_id
!sk2020sep10
      sk1 = minval(adt)
      if (sk1 .lt. 0) 
     &print*,'sep21:idea_phys>idea_phys_dissip-on:input/min_te='
     &  ,minval(adt),' me= ',mpi_id
!sk
      call idea_phys_dissipation(im,ix,levs,grav,prsi,prsl,             
     &      adu,adv,adt,o_n,o2_n,n2_n,dtp,cp,dt6dt)
!sk2020sep10
      sk1 = minval(adt)
      if (sk1 .lt. 0)
     &print*,'sep21:idea_phys<idea_phys_dissip-on:output/min_te='
     &  ,minval(adt),' me= ',mpi_id
!sk
      if (skprnt) print*,
     &           ' in idea_phys: <-- idea_phys_dissipation, me= ',mpi_id

      if (skprnt) print*,' idea_phys: ok009'
!SK2020Sep2
!     if (return_from_here) goto 1111
!
!--------------------------------------------------------------------------
!      endif
!============================================
!  end of dissipation +GW physics/turbulence
!============================================
!! VAY-11/2016  In WAM-2017 "validated" version:
! Do not UPDATE "adt-Temp-re" just collect Q-tendencies
!   sum them and MERGE with the  LA/GSM heating rates
!     and then UPDATE .............. adt !!!    
! to avoid "abrupt" heating changes      
!============================================
!
!
!  Start WAM radiation...o3_n ...compilation of global O3 + O3_gsm[1:k71]
!
      if (skprnt) print*,
     &           ' in idea_phys: --> o3pro, me= ',mpi_id
      call o3pro(im,ix,levs,ntrac,adr,am,nair,o3_ng)
      if (skprnt) print*,
     &           ' in idea_phys: <-- o3pro, me= ',mpi_id

      if (skprnt) print*,' idea_phys: ok010'
!SK2020Sep8
!     if (return_from_here) goto 1111
!
! get solar heating (EUV, UV-SRC-SRV-Lya) and NO cooling
! 
      if (skprnt) print*,
     &           ' in idea_phys: --> idea_sheat, me= ',mpi_id
!sk2020sep10
!     sk1 = minval(adt)
!     if (sk1 .lt. 0) 
!    &print*,'sk2020sep10:idea_phys>idea_sheat:input/min_te='
!    &  ,minval(adt),' me= ',mpi_id
!sk
      call idea_sheat(im,ix,levs,adt,dtRad,cospass,o_n,o2_n,o3_n,n2_n, 
     &                rho, cp,lat,dayno,prsl,zg,grav,am,maglat,dt6dt,
     &                f107_curdt, f107d_curdt, kpa_curdt)
!sk2020sep10
!     sk1 = minval(adt)
!     if (sk1 .lt. 0) 
!    &print*,'sk2020sep10:idea_phys<idea_sheat:output/min_te='
!    &  ,minval(adt),' me= ',mpi_id
!sk
      if (skprnt) print*,
     &           ' in idea_phys: <-- idea_sheat, me= ',mpi_id

      if (skprnt) print*,' idea_phys: ok011'
!
! Merge the  IPE back coupling WAM dtrad array into WAM.
!-------------------------------------------------------
      IF(ipe_to_wam_coupling) THEN
        DO k = lowst_ipe_level, levs
          DO i = 1, im
            GJHR(i, k) = GJHR(i, k) / cp(i, k)
            GSHR(i, k) = GSHR(i, k) / cp(i, k)
          END DO
        END DO

        call idea_merge_ipe_to_wam(GSHR, dtrad,
     &    im, ix, levs, lowst_ipe_level, prsl, plow, phigh,
     &    xpk_low, xpk_high)
      END IF
!-------------------------------------------------------
!     if (return_from_here) goto 1110
! 
! radiation
! co2 cooling, heating
      if (skprnt) print*,
     &           ' in idea_phys: --> idea_co2, me= ',mpi_id
      call idea_co2(im,ix,levs,nlev_co2,ntrac,cp,adr,adt,          
     &              dtco2c,cospass,dtco2h)
!SK   call idea_co2(im,ix,levs,nlev_co2,ntrac,grav,cp,adr,adt,          
!SK  &              dtco2c,cospass,dtco2h)
      if (skprnt) print*,
     &           ' in idea_phys: <-- idea_co2, me= ',mpi_id
      if (skprnt) print*,' idea_phys: ok012'

!SK2020Sep16
!     if (return_from_here) goto 1110
! h2o cooling heating 110-41 down ward
      if (skprnt) print*,
     &           ' in idea_phys: --> idea_h2o, me= ',mpi_id
      call idea_h2o(im,ix,levs,nlev_h2o,nlevc_h2o,ntrac,grav,cp,        
     &              adr,adt,dth2oh,cospass,dth2oc)
      if (skprnt) print*,
     &           ' in idea_phys: <-- idea_h2o, me= ',mpi_id
      if (skprnt) print*,' idea_phys: ok013'

!SK2020Sep16
!     if (return_from_here) goto 1110
! 
! o2 o3 heating
!
      if (skprnt) print*,
     &           ' in idea_phys: --> o3pro, me= ',mpi_id
      call o3pro(im,ix,levs,ntrac,adr,am,nair,o3_n)
      if (skprnt) print*,
     &           ' in idea_phys: <-- o3pro, me= ',mpi_id
      if (skprnt) print*,' idea_phys: ok014'
!SK2020Sep16
!     if (return_from_here) goto 1110
!======================================================================
! here 2 options: 1) o3_n (24-hr MLT O3) or 2) o3_ng (global prof > k71)
!
! for O3RAD =[o3_n, o3_ng] for Q in Hartley-Huggins-Herzberg-Chappius
!
!======================================================================
      if (skprnt) print*,
     &           ' in idea_phys: --> idea_o2_o3, me= ',mpi_id
      call idea_o2_o3(im,ix,levs,cospass,adt,o2_n,o3_ng,rho,cp,
     &                zg,grav,dto3)
      if (skprnt) print*,
     &           ' in idea_phys: <-- idea_o2_o3, me= ',mpi_id
      if (skprnt) print*,' idea_phys: ok015'
!SK2020Sep16
!
1110  continue
! get xmu as in "dcyc2.f"
!
      do i=1,im
!coszen (im)  - real, avg of cosz over daytime sw call interval
! daytime, XCOSZ=cosine of solar zenith angle at current time
!        if ( xcosz(i) > f_eps .and. coszen(i) > f_eps ) then
!          xmu(i) = xcosz(i) / coszen(i) f_eps =0.0001
!day
        if(cospass(i) > 0.0001 .and. cozen(i) > 0.0001) then
          xmu(i) = cospass(i)/cozen(i)
        else
!night
          xmu(i) = 0.
        endif
      enddo

      if (skprnt) print*,' idea_phys: ok016'
!SK2020Sep1
!     if (return_from_here) goto 1111
!====================================================================== 
! merge "ALL" heating rates dtdt-JouleHR ... dtRad (SRB+EUV-CNO)
!  here someone can put 1/cp/rho - factor don't divide inside heat-subs
!
! VAY-2017: new add thermospheric/MLT dtrad to merge take-out dt6dt
!
!====================================================================== 
      if (skprnt) print*,
     &           ' in idea_phys: --> rad_merge, me= ',mpi_id
      call rad_merge(im,ix,levs,hlw,swh,prsi,prsl,wtot,
     &               xmu,dtrad,dtco2c,dtco2h,dth2oh,dth2oc,dto3)
      if (skprnt) print*,
     &           ' in idea_phys: <-- rad_merge, me= ',mpi_id

      if (skprnt) print*,' idea_phys: ok017'
!SK2020Sep1
!     if (return_from_here) goto 1111

      do k=1,levs
        do i=1,im
          adt(i,k)     = adt(i,k) + dtp*wtot(i,k)
!sk06072018  Rashid's approach
          hlw(i,k) = wtot(i,k)
          swh(i,k) = 0.0
! dt6dt
          dt6dt(i,k,2) = wtot(i,k)
! 
        enddo
      enddo

      if (skprnt) print*,' idea_phys: b4 ok018'
!     if (return_from_here) goto 1111
!=========================================================================
! the last piece "simple" ionospheric block with empirical ION-RE moddels
!=========================================================================
      if (skprnt) print*,
     &           ' in idea_phys: --> idea_ion, me= ',mpi_id
      call idea_ion(mpi_id,prsl,solhr,cospass,zg, grav,o_n,o2_n,n2_n,cp,
     &              adu,adv,adt,dudt,dvdt,dtdt,rho,xlat,xlon,ix,im,levs,
     &              dayno,utsec,sda,maglon,maglat,btot,dipang,essa,
     &              f107_curdt, f107d_curdt, kp_curdt, nhp_curdt, 
     &              nhpi_curdt, shp_curdt, shpi_curdt, SPW_DRIVERS,
     &              swbz_curdt, swvel_curdt)
      if (skprnt) print*,
     &           ' in idea_phys: <-- idea_ion, me= ',mpi_id
! Merge the  IPE back coupling WAM dudt, dvdt and dtdt arrays into WAM.
!----------------------------------------------------------------------
      IF(ipe_to_wam_coupling) THEN
        call idea_merge_ipe_to_wam(GZMT, dudt,
     &    im, ix, levs, lowst_ipe_level, prsl, plow, phigh,
     &    xpk_low, xpk_high)
        call idea_merge_ipe_to_wam(GMMT, dvdt,
     &    im, ix, levs, lowst_ipe_level, prsl, plow, phigh,
     &    xpk_low, xpk_high)
        call idea_merge_ipe_to_wam(GJHR, dtdt,
     &    im, ix, levs, lowst_ipe_level, prsl, plow, phigh,
     &    xpk_low, xpk_high)
      END IF

      if (skprnt) print*,' idea_phys: ok018'
!SK2020Sep1
!     if (return_from_here) goto 1111
!
1111  continue
      do k=1,levs
        do i=1,im
          adu(i,k) = adu(i,k) + dtp*dudt(i,k)
          adv(i,k) = adv(i,k) + dtp*dvdt(i,k)
          adt(i,k) = adt(i,k) + dtp*dtdt(i,k)
        enddo
      enddo
!======================= WAM-IPE physics is completed ========
!     return
!=============================================================
!
!   debug print-outs
!
!================================================================
!1111  continue
!SK2020Sep6 to debug IDEA_PHYS
!     do k=1,levs
!       do i=1,im
!         adu(i,k) = adu_(i,k)
!         adv(i,k) = adv_(i,k)
!         adt(i,k) = adt_(i,k)
!       enddo
!     enddo
!     do l=1,ntrac
!     do k=1,levs
!       do i=1,im
!         adr(i,k,l) = adr_(i,k,l)
!       enddo
!     enddo
!     enddo

!1111  continue
      if (skprnt) then
       print *, 'idea_phys-T: ', maxval(adT), minval(adT),' me= ',mpi_id
       print *, 'idea_phys-U: ', maxval(adU), minval(adu),' me= ',mpi_id
       print *, 'idea_phys-V: ', maxval(adV), minval(adV),' me= ',mpi_id
      endif
      return
      end subroutine idea_phys
