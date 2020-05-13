!=============================================================================
! Rashid Akmaev, Dec 2017: Sun-Earth distance (SED) factor, some cleaning
!
! VAY-DEC 2016: Revision of idea_solar_heating.f with separation on
!               a) idea_solar_2014.f    - old non-tiegcm "idea"
!               b) idea_solar_heating.f - new tiegcm, SRB +Lya +slant columns
!               c) idea_solar_init.f    - time-inv fixed data and updated eff.
!                                         with 3 options to treat solar drivers
! Feb 2018     Zhuxiao Li and Tzu-Wei Fang, reading the 24hr ave. kp
!              (from driving parameter file) from idea_phys.f instead of kp. 
! Contains : idea_sheat
!            presolar  
!            solar_heat_dissociation_TIEGCM
!            idea_dissociation
!            get_slantcolumns
!=============================================================================
      subroutine idea_sheat(im,ix,levs,te,dt,cospass,
     & o_n,o2_n,o3_n,n2_n,     
     &ro,cp,lat,dayno,prsl,zg,grav,am,maglat,dt6dt, f107, f107d, kpa)
!----------------------------------------------------------------------------
! calculete solar heating, NO coooling from 2Pa up
!----------------------------------------------------------------------------
!
      use idea_solar, only : rgas, avgd, amo2, amo, amn2
      use idea_solar, only : nps
!                        
  
      use machine, only : kind_phys
      implicit none
      logical, parameter  :: sheat_hao=.true.
!
      integer, intent(in) :: im       !number of data piont in te
      integer, intent(in) :: ix       !maxmum data point in te
      integer, intent(in) :: levs     !number of press level
      integer, intent(in) :: lat      ! latitude index
      integer, intent(in) :: dayno    ! calender day
!
!VAY-2016
      real, intent(in)    :: f107, f107d, kpa ! solar-geo drivers
!
      real, intent(in)    :: te(ix,levs)  !temperature
      real, intent(in)    :: cospass(im)  ! cos zenith angle
      real, intent(in)    :: maglat(im)   ! in radians for NO-snoe
      real, intent(in)    :: cp(ix,levs)  ! 
      real, intent(in)    :: o_n(ix,levs) !number density of O(/cm3) 
      real, intent(in)    :: o2_n(ix,levs)!number density of O2
      real, intent(in)    :: o3_n(ix,levs)!number density of O2
      real, intent(in)    :: n2_n(ix,levs)!number density of N2
      real, intent(in)    :: am(ix,levs)  !mass of mix (kg)
      real, intent(in)    :: prsl(ix,levs)!layer press (Pa)
      real, intent(in)    :: zg(ix,levs)!layer height (m)
      real, intent(in)    :: grav(ix,levs)! (m/s2)
      real, intent(in)    :: ro(ix,levs)  ! density (kg/m3) 
      real, intent(inout) :: dt6dt(ix,levs,6)  ! 
      real, intent(out)    :: dt(ix,levs) ! (K/s) solar heating rate
! Locals
      integer  i,k
      real t(levs),n2(levs), o(levs),o2(levs),o3(levs)
      real :: ho(levs), ho2(levs),hn2(levs)
      real :: sheat(levs),qno(levs),no_new(levs),
     &amm(levs),prr(levs),alt(levs),nn(levs),sh1(levs),sh2(levs)
      real  :: ht(levs)
      real  :: dissociation_rate(levs), Jo3(levs)

      real :: vay_rgas_o, vay_rgas_o2,vay_rgas_n2 

      real :: vay_avgd, tg_vay, vay_cprho, rcpro, ron2

      vay_rgas_o =  1.e3*rgas/amo
      vay_rgas_o2 = 1.e3*rgas/amo2
      vay_rgas_n2 = 1.e3*rgas/amn2
      vay_avgd    = 1.e3*avgd 
      ron2 = amo/amn2          
      do i=1,im
        do k=1,levs
          o(k)=o_n(i,k)                            !/m3 in idea_phys
          o2(k)=o2_n(i,k)   
          o3(k)=o3_n(i,k)        
          n2(k)=n2_n(i,k)     
          t(k) = te(i,k)                 
          tg_vay = te(i,k)/grav(i,k)
          ho(k) =vay_rgas_o*tg_vay                  !m
          ho2(k)=.5*ho(k) 
          hn2(k)=ron2*ho(k)  
          ht(k) =  zg(i,k)
          alt(k) = 1.e-3*ht(k)
          prr(k)=prsl(i,k)
          nn(k)=o(k)+o2(k)+o3(k)+n2(k)
          amm(k)=am(i,k)*vay_avgd 
        enddo
!                  
! get heating

!
      IF (sheat_hao) then
       call solar_heat_dissociation_TIEGCM(levs,nps,o,o2,o3,n2,           
     &     ho,ho2,hn2, f107,f107d,cospass(i),dayno,         
     &     ht,sheat,sh1,sh2,dissociation_rate, Jo3)
       ELSE
         print *, ' OLD-SOLAR_EUV_SRC-2014'
!        call solar_heat(levs,nps,o,o2,n2,ho,ho2,hn2,effeuv,effuv,       
!     &   f107, cospass(i),sheat,sh1,sh2)
       ENDIF

!        print *, 'dissociation_rate=', dissociation_rate
!        print *, 'Jo3=', Jo3
!
!         print *, 'VAY idea_dissociation index ', i
!!        call getno(1,1,levs,maglat(i),dayno,alt,prr,nn,amm,no_new,
!!     &      f107, kpa)

       call getno1d(levs,f107,kpa,maglat(i),dayno,
     &      alt,prr,nn,amm,no_new)
!!     
!        print *, 'after call getno1d'
!        print *, 'f107', f107
!        print *, 'kpa=', kpa
              
!!        call COOLNO1(levs,nps,t,o,no_new,qno)   

       call WAM_COOLNO1(levs, nps, t, o, no_new, qno)                  
! 
        do k=nps,levs
!    
          rcpro = 1./(cp(i,k)*ro(i,k))
          dt6dt(i,k,1)=qno(k)*rcpro
          dt6dt(i,k,3)=sh1(k)*rcpro
          dt6dt(i,k,4)=sh2(k)*rcpro
!net Q
          dt(i,k)=(sheat(k)-qno(k))*rcpro
! 
        enddo ! k
        do k=1,nps-1
          dt(i,k)=0.
          dt6dt(i,k,1)=0.
          dt6dt(i,k,3)=0.
          dt6dt(i,k,4)=0.
          dt6dt(i,k,5)=0.
        enddo ! k

      enddo !i
! 
      return
      end
!---------------
      SUBROUTINE presolar(IM,IX,SOLHR,SLAG,                             
     &                   SINLAT,COSLAT,SDEC,CDEC,xlon,xlat              
     &                   ,XMU,dayno,utsec,sda                           
     &                   ,maglat,maglon,btot,dipang,essa)
!------------------------------------------------------------------------
! calculate solar zenith angle
!------------------------------------------------------------------------
      USE MACHINE     , ONLY : kind_phys
      USE PHYSCONS, PI => con_PI
      use date_def
      use wam_ion, only : getmag
      use wam_date_calendar, only : curddd_wam
      implicit none
!Argument
! input
      integer              IM,IX
      real(kind=kind_phys) sdec,slag,solhr,cdec
!vay-bug: IM not IX
      real(kind=kind_phys) SINLAT(IM),COSLAT(IM),XLON(IM),xlat(IM)
! output
      real     XMU(IM)       !cos solar zenith angle
! Output Magnetic and electric parameters 
!     REAL, INTENT(OUT) :: elx(im)
!     REAL, INTENT(OUT) :: ely(im)     !electric field
      REAL, INTENT(OUT) :: maglon(im)  !magnetic longitude (rad)
      REAL, INTENT(OUT) :: maglat(im)  !magnetic latitude (rad)
      REAL, INTENT(OUT) :: btot(im)    !mapgnetic field strength
      REAL, INTENT(OUT) :: dipang(im)  !Dip angle (degree)
      REAL, INTENT(OUT) :: essa(im)    !magnetic local time
      REAL, INTENT(OUT) :: sda         ! solar declination angle (rad)
! Output time parameters 
      REAL, INTENT(OUT) :: utsec       !universal time
      INTEGER, INTENT(OUT) :: dayno    !calendar day
! local vareable
      INTEGER   i,idat(8),jdat(8),jdow,jday 
      real(kind=kind_phys) cns,ss,cc,ch,rinc(5),ty
!  COMPUTE COSINE OF SOLAR ZENITH ANGLE FOR BOTH HEMISPHERES.
      CNS = PI*(SOLHR-12.)/12.+SLAG
      DO I=1,IM
        SS     = SINLAT(I) * SDEC
        CC     = COSLAT(I) * CDEC
        CH     = CC * COS(XLON(I)+CNS)
        XMU(I) = CH + SS
      ENDDO
! get day number year number UTsec
      idat=0
      idat(1)=idate(4)
      idat(2)=idate(2)
      idat(3)=idate(3)
      idat(5)=idate(1)
      rinc=0.
      rinc(2)=fhour
      call w3movdat(rinc,idat,jdat)
      call w3doxdat(jdat,jdow,dayno,jday)
!
!     dayno = curddd_wam
!     print*,'VAY',dayno,fhour
      utsec=solhr*3600.
! get solar declination angle
      ty = (dayno+15.5)*12./365.
      IF ( ty > 12.0 ) ty = ty - 12.0
      sda = ATAN(0.434*SIN(PI/6.0*(ty-3.17)))
!
!
      call getmag(im,utsec,xlat,xlon,sda,                            
     &btot,dipang,maglon,maglat,essa)
      btot=btot*1.e-9
      RETURN
      END
!
!
      SUBROUTINE solar_heat_dissociation_TIEGCM(np,nps,O,O2,O3,N2,
     &  HO,HO2,HN2,  F107,F107d,COSPASS,dayno,height,
     &  sheat,sh1,sh2, O2dissociation_rate, O3dissociation_rate)

! RAA Dec 2017: SED correction factors, some cleanup
!vay-2015   solar_heat_dissociation_TIEGCM   1D_height subroutine

       use idea_composition, only  : R0 => REARTH, pi
       use idea_solar,       only  : euv37, nsp_euv ! EUV-flux(nsp_euv=37)
       use idea_solar,       only  : o2_scale_factor,  SRBEFF  ! (levs)-arrays
       use idea_solar,       only  : effuv, effeuv             ! CTIP-orig
       use idea_solar,       only  : eff_src, eff_srb, eff_lya ! (levs)-arrays
       use idea_solar,       only  : rwpcc                     ! pcc/wavelength
       use idea_solar,       only  : csao, csao2, csan2, csao3
       use idea_solar,       only  : csio, csio2, csin2
       use idea_solar,       only  : csdo2, csdeo2
!       use idea_solar,       only  : eccentric
       use idea_solar,       only  : sfmin, afac
       use idea_solar,       only  : nwaves, nwaves_euv
       use idea_solar,       only  : lyman_a_num,nwaves_src
!-------------------------------------------------------------------------
! calculate solar heating and dissociation rates
! Fluxes and cross sections from QRJ.F in TIEGCM (TWFang, Jan 2015) => 
!        moved to idea_sola_init and passed through "idea_solar"
!-------------------------------------------------------------------------
!  **
!  calculates solar heating and dissociation rates 
!     from EUV and UV(SRC+LYA+SRB)
!  
!  Input:
!  O atomic oxygen number density profile m-3
!  O2 molecular oxygen number density profile m-3
!  N2 molecular nitrogen number density profile m-3
!  HO atomic oxygen scale height profile m
!  HO2 molecular oxygen scale height profile m
!  HN2 molecular nitrogen scale height profile m
!  F107 solar flux 
!  F107d 81 day mean F10.7
!  COSPASS cosine of solar zenith angle
!  dayno day of year ( 1 to 366)
!
!  Output:
!  SHEAT, sh1(euv), sh2(uv) heating rates J/m3
!  O2dissociation_rate dissociation rates s-1
!---------------------------------------------------------------------------
      implicit none
      integer, parameter  :: idea_solar_fix=1 
      logical, parameter  :: para_schemes=.false.
!
! input
!
      integer, intent(in) :: np      ! number of pressure levels
      integer, intent(in) :: nps     ! pressure index to start
      real, intent(in)    :: o(np),o2(np),n2(np)       ! number density/m3
      real, intent(in)    :: o3(np)                    ! ozone density 1/m3
      real, intent(in)    :: ho(np),ho2(np),hn2(np)    ! scale height(m)

      real, intent(in)    :: f107       ! f10.7cm
      real, intent(in)    :: f107d      ! 81 day mean f10.7
      real, intent(in)    :: cospass    ! cos zenith angle
      real, intent(in)    :: dayno      ! day of year
      real, intent(in)    :: height(np) !layer height (m)

! output

      real, intent(out)   :: sheat(np),sh1(np),sh2(np)   ! W/m3 heating rate
      real, intent(out)   :: O2dissociation_rate(np)
      real, intent(out)   :: O3dissociation_rate(np)

! locals 

      real sPAEUV              ! VAY scalar instead of 2D-array PAEUV
      real Z, coschi
      real seco,seco2, seco3, secn2
      real wo,wo2, wo3, wn2
      real tau,tauo,tauo2,tauo3,taun2 
      real nightfac
      real  attenuation, local_flux
      real rnight_o,rnight_o2,rnight_o3,rnight_n2

! Solar variability factor for SRB dissociation calculation
      REAL FLUX(NWAVES)
      REAL fmxfmn
      REAL pind
      real, dimension(np) :: sco2, sco, scn2,sco3 ! slant columns o2/o/n2
! local  O2 dissociation rates
      REAL pre_loc, PAEUV_loc,  JSRC_loc,                       
     &     JLYA_loc, JSRB_loc,  JEUV_loc
! O3- local
      REAL :: Jo3_euvloc,  Jo3_srcloc, Jo3_lyaloc, pdnolya
      REAL :: Jo3_harloc, Jo3_chhloc
! Lyman-a and SRC, SRB and EUV O2 dissociation rate coefficients
      REAL JLYA(np), JSRC(np), JSRB(np),JEUV(np)
! SRB
      REAL srband, srband_fac,srbheat(np)
      real, dimension(np) :: shsrc, shsrb, shlya

! vay-2015

      real :: vay1, vay2, vay3, vay4, rnight
!   
      real :: vay_fmxfmn, vay_srband_fac, vay_srb3
      real Zrnight_o2(np)
      real :: sfeps
      real :: temp1, sco3t
!
      integer i,j,jinv

! RAA: solar 10-cm flux normalized to 1 AU
      real f107_1au, f107d_1au

        nightfac=1.e-9         ! Ratio of nightime ionisation to sec=1.

! RAA: reintroduce the SED factor, normalize solar flux to 1 AU
!        sfeps   = 1.0
        sfeps   = (1.0+0.0167*COS(2.0*pi*(dayno-3.)/366.0))**2
        f107_1au = f107/sfeps
        f107d_1au = f107d/sfeps

!vay-2017 fmxfmn=(F107-71.)/(220-71.) does not make sense for F107~70.
        fmxfmn=(f107_1au-65.)/165.   
        srband_fac=0.784591675
        vay_fmxfmn=(1.+0.11*fmxfmn)
        vay_srband_fac=(1.+srband_fac)*.1

! RAA: replace eccentric with sfeps
!        vay_srb3=vay_srband_fac * vay_fmxfmn * eccentric
        vay_srb3=vay_srband_fac * vay_fmxfmn * sfeps

        pind=0.5*(f107_1au+f107d_1au) -80.
!

      do J = 1, NWAVES  
        jinv = nwaves+1-J
!
         if (idea_solar_fix.ge.1) then 
!parameterized EUV 
! RAA: normalize flux by SED
!           flux(J)=sfmin(jinv)*max(0.8, (1.0+afac(jinv)*pind) )
           flux(J)=sfmin(jinv)*max(0.8, (1.0+afac(jinv)*pind) )*sfeps
         else              
! solar_fix == 0 use observed EUV               
            flux(J)=euv37(j)                                          
         endif
         IF (flux(J) .LT. 0.0) flux(J) = 0.0
      enddo
!      
      goto 7771
      COSCHI=COSPASS
      rnight=1.0
      if(coschi.lt.0.07)then
        rnight_o=1.e-6
        SECO=1./0.07      ! was 1.0 in 2013
      else
        rnight_o=1.0
        SECO=1./COSCHI         
      end if
        SECO2=SECO
        SECN2=SECO
        rnight_o2 =rnight_o
        rnight_n2 =rnight_o
 7771  continue
 
          sheat(1:np)=0.
          sh1(1:np)=0.
          sh2(1:np)=0.
          shsrc(:) =0.
          shsrb(:) =0.
          shlya(:) =0.
          O2dissociation_rate(1:np)=0. 
          O3dissociation_rate(1:np)=0. 
!
! compute slant columns in 1/cm2
!
      call get_slantcolumns(np, nps, R0, height, cospass,
     &     HO, HO2, HN2, O, O2, O3,N2, sco, sco2, sco3, scn2)
!
! add extra: (seco, rnight_o)-arrays to avoid 2-nd calls of SUB_CHAPMAN
!    
      do i=nps,np
!                                   Set total height in m Rad+Zmeters
       z = R0 + height(i)

       jlya(i)   = 0.
       jsrc(i)   = 0.
       jsrb(i)   = 0.
       jeuv(i)   = 0.
       srband    = 0.
       srbheat(i)= 0.

       pre_loc     = 0.
       JSRC_loc    = 0.
       JSRB_loc    = 0.
       JLYA_loc    = 0.
       JEUV_loc    = 0.
       Jo3_euvloc    = 0.
       Jo3_lyaloc    = 0.
       Jo3_srcloc    = 0.
       TAU=0.
!  **
! calculate sec(ZA), incorporating Chapmann grazing incidence function
      CALL SUB_CHAPMAN(cospass, HO(i), z, nightfac, seco, rnight_o) 
      CALL SUB_CHAPMAN(cospass, HO2(i),z, nightfac, seco2,rnight_o2)
      CALL SUB_CHAPMAN(cospass, HO(i)*.333,z,nightfac,seco3,rnight_o3)
      CALL SUB_CHAPMAN(cospass, HN2(i),z, nightfac, secn2,rnight_n2) 
       Zrnight_o2(i) = rnight_o2

! For tau convert column abundance from m^-2 to cm^-2 (10^-4)
        WO= O(i) *HO(i) *SECO*1.e-4
        WO2=O2(i)*HO2(i)*SECO2*1.e-4
        WN2=N2(i)*HN2(i)*SECN2*1.e-4
        WO3 =sco3(i)

!  loop over all 37 bands
       sh1(i)=0.
       sh2(i)=0.
       O2dissociation_rate(i)=0.
       O3dissociation_rate(i)=0.
      do J = 1, NWAVES
         TAUO =CSAO(J)*WO
         TAUO2=CSAO2(J)*WO2
         TAUN2=CSAN2(J)*WN2
         TAUO3=CSAO3(J)*WO3
         TAU=TAUO+TAUO2+TAUN2

! use parameters: solar_tuny = 1.e-36, min_flux =1.e-20, tau_max=200.
      IF (tau .ge.  1.e-36 .AND. tau <= 200.) then 
           attenuation = exp(-tau)
      else  if (tau < 1.e-36 ) then  
            attenuation = 1.0
      else  
            attenuation=0.
      ENDIF

! RAA: Remove eccentric, flux is already normalized to SED
!          local_flux=flux(J)*attenuation*eccentric
          local_flux=flux(J)*attenuation

          IF(local_flux < 1.e-20 ) local_flux=0.

! EUV heating calculation(W/m3)
          IF(J <= NWAVES_EUV) THEN
 
! For SI units convert flux to m^-2 (x10^4)
!  convert cross sections to m^2 (x10^-4)
         sPAEUV=local_flux*(CSAO(J)*O(i)*rnight_o+                  
     &    CSAO2(J)*O2(i)*rnight_o2+CSAN2(J)*N2(i)*rnight_n2)*RWPCC(J)
         sh1(i)    =sh1(i)  + Spaeuv

! EUV JO2
          Jeuv_loc = Jeuv_loc+local_flux*(CSDO2(J)+CSDeO2(J))

! EUV JO3
          Jo3_euvloc = Jo3_euvloc +local_flux*CSAO3(J)
         ELSE

!  UV SRC/LYa channels 
         sPAEUV=local_flux*CSAO2(J)*O2(i)*RWPCC(J)*rnight_o2

! calculate JSRC, excluding JLYA
           IF(J /= lyman_a_num) then 
              JSRC_loc=JSRC_loc+local_flux*CSAO2(J)
              shsrc(i)  =shsrc(i)  +spaeuv
              Jo3_srcloc = Jo3_srcloc +local_flux*CSAO3(J)
           ELSE                                                

! calculate JLYA 
              JLYA_loc=local_flux*CSAO2(J)
              shlya(i)  =  spaeuv                              !* effuv(i) 
              pdnolya = (0.68431  *exp(-8.22114e-21*sco2(i))+
     &             0.229841 *exp(-1.77556e-20*sco2(i))+
     &             0.0865412*exp(-8.22112e-21*sco2(i)))*
     &             flux(lyman_a_num)
              Jo3_lyaloc =2.27e-17*pdnolya
           ENDIF                                               ! src/lya
         ENDIF          ! UV/EUV
        enddo           ! end of wavelength loop - J-index

!==========================================================
!  qtotal(k,i,lat) = qtotal(k,i,lat)+ho2src(k,i)+ho2srb(k,i) 
!  Calculate O2 Schumunn Runge band heating  see "o2srbc.F"
!==========================================================
      IF (WO2 < 1.E18) THEN
         srband=O2(i)*2.43E-25  !!1.e-19*1.e-6  Strobel-78 exp(21) 20% accuracy
      ELSE
         srband=O2(i)*1.e-6/(0.66*WO2+3.44E9*SQRT(WO2))
      ENDIF
      
      srband =srband * vay_srb3 ! (vay_srband_fac*vay_fmxfmn*eccentric)
      if(srband.lt.0.)   srband  = 0.

!      sum srb+src+lya heating into the total UV heat
!=====================================================
!      apply Heat-efficiency factors & rnight_o2
      shsrb(i) = srband*rnight_o2
      sh1(i)   = sh1(i)  *effeuv(i)                
      shsrc(i) = shsrc(i)*effuv(i)                 
      shlya(i)  = shlya(i)*effuv(i)                
      sh2(i)   = shsrc(i) + shlya(i) +shsrb(i)     

! total
      sheat(i) = sh1(i) +sh2(i) ! EUV + UV

!===============================================
! calculate JO2 due to SRB
      IF (WO2.LT.1.E19) THEN

! RAA: add SED factor
!     JSRB_loc = vay_fmxfmn*1.1E-7*EXP(-1.97E-10*(WO2**0.522))  
         JSRB_loc = vay_fmxfmn*1.1E-7*EXP(-1.97E-10*(WO2**0.522))*sfeps
         
      ELSE
!     JSRB_loc =vay_fmxfmn*1.45E8*(WO2**(-0.83))
         JSRB_loc =vay_fmxfmn*1.45E8*(WO2**(-0.83))*sfeps
      ENDIF

! sum dissociation rates, eccentric is inclued in the local_flux and SRB
 
      JSRC(i)   = JSRC_loc
      JSRB(i)   = JSRB_loc          
      JLYA(i)   = JLYA_loc
      JEUV(i)   = JEUV_loc
!     
      O2dissociation_rate(i) =O2_scale_factor(i)*rnight_o2*
     &     (JSRC_loc+JLYA_loc+JSRB_loc+JEUV_loc)
!
! Jo3_harloc & Jo3_chhloc
!         
      sco3t = sco3(i)                        
      if (sco3t < 1.e+5) sco3t = 1.e+5
      temp1 = 1.0e-3*exp(-1.5577e-13*sco3t**0.6932)
      if (sco3t < 1.6e+20) 
     &     temp1 = 1.04e-2*exp(-1.0217e-6*sco3t**0.3587)
! Hartley bands of O3:
!          Jo3_harloc = sfeps*0.68*
!     &      (temp1 + 1.085*exp(-1.4912e-14*sco2(i)**0.5298)*
!     &      ( 4.053e-4  *exp(-8.1381e-16*sco3t**0.8856) +
!     &        4.700e-6  *exp(-1.5871e-14*sco3t**0.7665)   )*
!     &       *exp(-1.4655e-25*sco2(i)**1.0743) )
      Jo3_harloc = 
     &     (temp1*0.68+exp(-1.4912e-14*sco2(i)**0.5298)*
     &     (4.053e-4  *exp(-8.1381e-16*sco3t**0.8856)+
     &     4.7e-6    *exp(-1.5871e-14*sco3t**0.7665))*
     &     1.085     *exp(-1.4655e-25*sco2(i)**1.0743)*0.68)*sfeps   
!
! Chappius and Huggins bands:
!
      Jo3_chhloc = sfeps*
     &     (    4.5e-4*exp(-3.4786e-11*sco2(i)**0.3366
     &     -1.0061e-20*sco3t  **0.9719)
     &     + ( 7.5e-4 *exp(-2.7663e-19*sco3t**1.0801 )
     &     +2.5e-4/(1.+1.5772e-18  *sco3t**0.9516))
     &     *exp(-1.0719e-10*sco2(i)**0.3172 )  )
!
      O3dissociation_rate(i) = rnight_o3*(Jo3_EUVloc+
     &     Jo3_SRCloc+Jo3_LYAloc+Jo3_harloc+Jo3_chhloc)
      enddo    ! height loop
!
!
      do i=1,nps-1
         sco3t =  exp(-abs(float(i-nps))*0.2) 
!
! Jo3 (O3P+O1D) below ~50 km is ~ constant with small decrease 
!

         O2dissociation_rate(i) = O2dissociation_rate(nps)*sco3t 
         O3dissociation_rate(i) = O3dissociation_rate(nps)*sco3t  
      enddo
!
      RETURN
      END SUBROUTINE solar_heat_dissociation_TIEGCM
!
!=====================================================
!
!tzw- 10/2015 new subroutine idea_dissociation, initial implement => 1D-arr
!vay  oct/2016  move to dissociation_rate2d, needed by ion_tracer and
!     correct errors of 2015
!
      subroutine idea_dissociation_jo3(im,ix,levs,te,cospass,
     & o_n,o2_n,o3_n, n2_n,
     & dayno,zg,grav, f107, f107d, Jo2_2d, Jo3_2d)
!----------------------------------------------------------------------------
! calculete solar dissociation of O2 (UV+EUV)
!----------------------------------------------------------------------------
!   
      use idea_solar, only : nps

!      use idea_solar, only : f107, f107a
      use idea_solar, only : avgd, rgas, amo, amn2, amno, amo2

      use machine, only    : kind_phys
      implicit none
!input
      integer, intent(in) :: im           !number of data piont in te
      integer, intent(in) :: ix           !maxmum data points reserved for 2D-arrays
      integer, intent(in) :: levs         !number of press level
      integer, intent(in) :: dayno        ! calender day
      real, intent(in)    :: te(ix,levs)     !temperature
      real, intent(in)    :: cospass(im)  ! cos zenith angle
      real, intent(in)    :: o_n(ix,levs) !number density of O(1/m3)
      real, intent(in)    :: o2_n(ix,levs)!number density of O2
      real, intent(in)    :: o3_n(ix,levs)!number density of O3
      real, intent(in)    :: n2_n(ix,levs)!number density of N2
      real, intent(in)    :: zg(ix,levs)  !layer height (m)
      real, intent(in)    :: grav(ix,levs)! (m/s2)
      real, intent(in)    :: f107, f107d
!
! VAY out dissociation_rate2d
!
      real, intent(out)   :: Jo2_2d(ix, levs), Jo3_2d(ix, levs)
!
!locals
      integer :: i,k
      real    :: n2(levs), o(levs),o2(levs), o3(levs)
      real    :: ho(levs), ho2(levs),hn2(levs)
      real    :: sheat(levs),sh1(levs),sh2(levs)
      real    :: ht(levs)
      real    ::  jo2_1d(levs), jo3_1d(levs)
! TWFANG
       real :: vay_rgas_o, vay_rgas_o2,vay_rgas_n2 
!
      real :: tg_vay, rodn2
!nullify
!     
       Jo2_2d (:,:) =0.0 
       Jo3_2d (:,:) =0.0    
!
      vay_rgas_o =  1.e3*rgas/amo
!      vay_rgas_o2 = 1.e3*rgas/amo2
!      vay_rgas_n2 = 1.e3*rgas/amn2
      rodn2 = amo/amn2
!
      do i=1,im
        do k=1,levs
          o(k)=o_n(i,k)                             !/m3
          o2(k)=o2_n(i,k)
          o3(k)=o3_n(i,k)
          n2(k)=n2_n(i,k)
      
          ht(k)  = zg(i,k)                          !layer height (m)
          tg_vay = te(i,k)/grav(i,k)
          ho(k) =vay_rgas_o*tg_vay                  !m
          ho2(k)=.5*ho(k) 
          hn2(k)= rodn2*ho(k)  
!          ho(k)=1.e3*rgas*t(k)/(amo*grav(i,k))     !m
!
!         
        enddo
! 
        call solar_heat_dissociation_TIEGCM(levs,nps,o,o2,o3,n2,
     &     ho,ho2,hn2,f107,f107d,cospass(i),dayno,
     &     ht,sheat,sh1,sh2, jo2_1d, jo3_1d)
! 
!VAY-oct 2016 .....dissociation_rate2d
!VAY-jan 2017    Jo2_2d & Jo3_2d based on TIME-GCM/2015
!
          Jo2_2d(i,nps:levs) = jo2_1d(nps:levs)
          Jo3_2d(i,nps:levs) = jo3_1d(nps:levs)
!
      enddo        !i-hor index
!      print *, 'VAY idea_dissociation '
      return
!
      end subroutine idea_dissociation_jo3
!
      subroutine idea_dissociation_jo2(im,ix,levs,te,cospass,
     & o_n,o2_n,o3_n, n2_n,
     & dayno,zg,grav, f107, f107d, Jo2_2d)
!----------------------------------------------------------------------------
! calculete solar dissociation of O2 (UV+EUV)
!----------------------------------------------------------------------------
!   
      use idea_solar, only : nps

!      use idea_solar, only : f107, f107a
      use idea_solar, only : avgd, rgas, amo, amn2, amno, amo2

      use machine, only    : kind_phys
      implicit none
!input
      integer, intent(in) :: im           !number of data piont in te
      integer, intent(in) :: ix           !maxmum data points reserved for 2D-arrays
      integer, intent(in) :: levs         !number of press level
      integer, intent(in) :: dayno        ! calender day
      real, intent(in)    :: te(ix,levs)     !temperature
      real, intent(in)    :: cospass(im)  ! cos zenith angle
      real, intent(in)    :: o_n(ix,levs) !number density of O(1/m3)
      real, intent(in)    :: o2_n(ix,levs)!number density of O2
      real, intent(in)    :: o3_n(ix,levs)!number density of O3
      real, intent(in)    :: n2_n(ix,levs)!number density of N2
      real, intent(in)    :: zg(ix,levs)  !layer height (m)
      real, intent(in)    :: grav(ix,levs)! (m/s2)
      real, intent(in)    :: f107, f107d
!
! VAY out dissociation_rate2d
!
      real, intent(out)   :: Jo2_2d(ix, levs)
!
!locals
      integer :: i,k
      real    :: n2(levs), o(levs),o2(levs), o3(levs)
      real    :: ho(levs), ho2(levs),hn2(levs)
      real    :: sheat(levs),sh1(levs),sh2(levs)
      real    :: ht(levs)
      real    ::  jo2_1d(levs), jo3_1d(levs)
! TWFANG
       real :: vay_rgas_o, vay_rgas_o2,vay_rgas_n2 
!
      real :: tg_vay, rodn2
!nullify
!     
       Jo2_2d (:,:) =0.0 
!
      vay_rgas_o =  1.e3*rgas/amo
!      vay_rgas_o2 = 1.e3*rgas/amo2
!      vay_rgas_n2 = 1.e3*rgas/amn2
      rodn2 = amo/amn2
!
      do i=1,im
        do k=1,levs
          o(k)=o_n(i,k)                             !/m3
          o2(k)=o2_n(i,k)
          o3(k)=o3_n(i,k)
          n2(k)=n2_n(i,k)
      
          ht(k)  = zg(i,k)                          !layer height (m)
          tg_vay = te(i,k)/grav(i,k)
          ho(k) =vay_rgas_o*tg_vay                  !m
          ho2(k)=.5*ho(k) 
          hn2(k)= rodn2*ho(k)  
!          ho(k)=1.e3*rgas*t(k)/(amo*grav(i,k))     !m
!
!         
        enddo
! 
        call solar_heat_dissociation_TIEGCM(levs,nps,o,o2,o3,n2,
     &     ho,ho2,hn2,f107,f107d,cospass(i),dayno,
     &     ht,sheat,sh1,sh2, jo2_1d, jo3_1d)
! 
!VAY-oct 2016 .....dissociation_rate2d
!VAY-jan 2017    Jo2_2d & Jo3_2d based on TIME-GCM/2015
!
          Jo2_2d(i,nps:levs) = jo2_1d(nps:levs)
!
      enddo        !i-hor index
!      print *, 'VAY idea_dissociation '
      return
!
      end subroutine idea_dissociation_jo2
!
      subroutine get_slantcolumns(levs, lev1, R0, Zgi, cospass,
     & HO, HO2, HN2, xO, xO2, xo3, xN2, sco, sco2, sco3, scn2)
       IMPLICIT NONE
!
! compute slant columns for 3-major species
!    
      integer,intent(in) :: levs,lev1
      real, intent(in)   :: R0, Zgi(levs)
      real, intent(in)   :: cospass
      real, intent(in), dimension(levs)   :: HO, HO2, HN2           ! m
      real, intent(in), dimension(levs)   :: xO, xO2, xN2,xo3       ! 1/m3
      real, intent(out), dimension(levs)  :: scO, scO2,sco3,scN2    ! 1/cm2
      real :: z ! meters
      real, parameter :: nightfac=1.e-6
      real, parameter :: smfac=1.e-4
      real :: seco,  rnight_o
      real :: seco2, rnight_o2
      real :: seco3, rnight_o3
      real :: secn2, rnight_n2

      real, dimension(levs)  :: vcO, vcO2,vco3,vcN2               ! 1/cm2
      real :: dzm, Ho3
      real :: rodfac
      integer :: i,k,j
!
! from top 2 bottom 
!rodfac=35./sqrt(1224.*cosz(i)**2+1.
      k=levs
!vert
      ho3 =HO(levs)/3.
      vco(levs)  = xo(levs)*smfac*HO(levs)
      vco2(levs) = xo2(levs)*smfac*HO2(levs)
      vco3(levs) = xo3(levs)*smfac*Ho3
      vcn2(levs) = xn2(levs)*smfac*HN2(levs)

         z = R0 +Zgi(k)
        CALL SUB_CHAPMAN(cospass, HO(k),  z,  nightfac, seco, rnight_o) 
        CALL SUB_CHAPMAN(cospass, HO2(k), z, nightfac, seco2, rnight_o2)
        CALL SUB_CHAPMAN(cospass, HN2(k), z, nightfac, secN2, rnight_n2) 
        CALL SUB_CHAPMAN(cospass, Ho3,    z, nightfac, seco3, rnight_o3) 
! slant
        sco(k)=  vco(k)*SECO*smfac
        sco2(k)= vco2(k)*SECO2*smfac
        sco3(k)= vco3(k)*SECO3*smfac
        scn2(k)= vcn2(k)*SECN2*smfac

      do k=levs-1, lev1, -1
         dzm = .5*(zgi(k+1)-zgi(k))
         vco(k)  = vco(k+1)  + (xo(k+1)+xo(k))*dzm 
         vco2(k) = vco2(k+1) + (xo2(k+1)+xo2(k))*dzm
         vcn2(k) = vcn2(k+1) + (xn2(k+1)+xn2(k))*dzm 
         vco3(k) = vco3(k+1) + (xo3(k+1)+xo3(k))*dzm 
         ho3 =HO(k)*0.3333
         z = R0 +Zgi(k)
      CALL SUB_CHAPMAN(cospass, HO(k), z,  nightfac,  seco, rnight_o ) 
      CALL SUB_CHAPMAN(cospass, HO2(k), z, nightfac, seco2, rnight_o2)
      CALL SUB_CHAPMAN(cospass, HO3,    z, nightfac, seco3, rnight_o3)
      CALL SUB_CHAPMAN(cospass, HN2(k), z, nightfac, secN2, rnight_n2) 
!
!  transform Vcol  => Scol
!       rodfac=35./sqrt(1224.*cospass(i)*cospass(i)+1.)
!
        sco(k)=  vco(k) *SECO *smfac
        sco2(k)= vco2(k)*SECO2*smfac
        sco3(k)= vco3(k)*SECO3*smfac
        scn2(k)= vcn2(k)*SECN2*smfac
!
!        WN2=N2(i)*HN2(i)*SECN2*1.e-4
!
      enddo
      end subroutine get_slantcolumns
!
!
      SUBROUTINE SUB_CHAPMAN(COSCHI,SCALE_HT,HT, NFAC, seco, rnight)
!=============================================================
! see Smith & Smith JGR 1972 for details,  v77, N19, p. 3592
!
! VAY oct 26/2016 review and corrections for undefined "PI"
!     ... etc..FUNCTION CHAPMANN(COSCHI,SCALE_HT,HT)
!! function is correponds to sec(zenith angle) 1./coschi
! bug was non-defined pi2 = pi/2 see original CTIPE-module
!=============================================================

      use idea_composition, only : PI, PI2, Pid2, R0 => REARTH, R_2_D
      IMPLICIT NONE     ! adding this .... major bugs are on compilation
!
      REAL, INTENT(IN)  :: SCALE_HT           ! meters
      REAL, INTENT(IN)  :: HT                 ! meters
      REAL, INTENT(IN)  :: COSCHI             !dim-ls
      REAL, INTENT(IN)  :: NFAC               ! 1.e-6 ????
      REAL, INTENT(OUT) :: SECO               ! 1/cosCHI
      REAL, INTENT(OUT) :: RNIGHT             ! 1.e-6 or less

!      REAL, PARAMETER :: R0    = 6.370E06    ! radius of earth (m)

      REAL ::  Y_ERR, RAD_TO_Z, CHI, CHID
      REAL ::  ERFCL                          ! don't use ERFC=> ERFCL
      REAL ::  SCHI
      CHI  = ACOS(COSCHI)
      CHID = CHI*R_2_D
!
! calculate chapman bit for angles over 75 degrees
!
      SCHI = sin(chi)

      IF(chid.GT.75. .and. chid .LT. 105.) THEN

          RAD_TO_Z = HT/SCALE_HT
!
! Step1: calc error function ERFC(X =SQRT(0.5*RAD_TO_Z) * ABS(COSCHI))
!
          Y_ERR = SQRT(0.5*RAD_TO_Z) * ABS(COSCHI)
          IF (Y_ERR.LE.8.0) THEN
            ERFCL = (1.0606963 + 0.55643831* Y_ERR) /            
     &      (1.0619898 + 1.7245609* Y_ERR + Y_ERR *Y_ERR)
          ELSE
            ERFCL = 0.56498823 / (0.06651874 + Y_ERR)
          ENDIF

! step2. Calculate chapmann             for solar zenith angles <= 90
!
       IF(CHID.LE.90.0)THEN
         SECO = SQRT(PId2 * RAD_TO_Z) * ERFCL
        ELSE
!                                      for solar  zenith angles > 90 (equation 15)
        SECO = SQRT(PI2 * RAD_TO_Z)*                      
     &   ( SQRT(SCHI)*EXP(RAD_TO_Z*(1-SCHI)) -.5*ERFCL )
!         write(6,*)'chapman over 90', CHID,CHAPMANN
        ENDIF

      ELSE                             ! out [75-105 degrees]range
!                                      CHAPMANN-SECO ~  sec(zenith angle)
          SECO = 1./COSCHI
      ENDIF
!
! for nighttime SECO < 0
!
          rnight=1.0
       IF(SECO < 0) THEN
          rnight=NFAC        !nightfac=1.e-6
          SECO=1.
       ENDIF

      END SUBROUTINE SUB_Chapman
