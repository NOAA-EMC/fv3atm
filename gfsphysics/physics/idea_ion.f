! Apr 06 2012 Henry Juang, initial implement for nems
! Dec    2012    Jun Wang, move init out of column physics
!
!====================================================================
! Nov    2015   VAY  SUB-NE IDEA_ION_INIT(levs) & idea_ion_input.f 
! Nov    2016   VAY  Upgrades related to  idea_ion_input.f
!                                          idea_ion_empirmodels.f   
! Jul    2017   Zhuxiao Li, Rashid and Tim add the Joule heating factor(JH_fac) to 
!               include the seasonal variation and semiannual variation of the
!               Joule Heating.
! Mar    2018   Zhuxiao Li and Tzu-Wei Fang add the new features to put
! the driving parameters into the idea_geteb before call get_efield. 
!====================================================================
!=                    GetIonParams                      =
!
! it "employs" Maute efield.f with read of the "external"
! vay-2015 data on unit=10 !!!!!!! danger, lunit should be replaced > 100  
!     see .... subroutine read_acoef
!
!    not clear is it active in WAM
!    the best is to " use efield only :  "name of variables employed in WAM ed1/ed2 etc... '
!      curdate_wam = jdat
!      curhour_wam = hr_jdat 
!      curyear_wam = jdat(1)
!      curmonth_wam = jdat(2)
!      curddd_wam =  JDAY - JDAY0101 +1
!      curday_wam =  jdat(3)
!      curjulday_wam  = JDAY
!      curutsec_wam = nint(hr_jdat*3600.) 
!========================================================
      SUBROUTINE IDEA_ION_INIT(me, master, levs)
!SK   SUBROUTINE IDEA_ION_INIT(levs)
!
      use IDEA_IO_UNITS, only : nml_ion, nlun_ion, ch100
!SK   use IDEA_MPI_def,  only : mpi_id, mpi_err, MPI_COMM_ALL,info
!
! subroutines
      use IDEA_ION_INPUT, only : ion_read_wam_init, tiros_read_wam_init
      use IDEA_ION_INPUT, only : ion_read_namelist,precomp_iondata_fixed
      use IDEA_ION_INPUT, only : tiros_activity_fixnam
!
!data-ion
      use IDEA_ION_INPUT, only : nxmag, nymag, glon, glat
      use IDEA_ION_INPUT, only : cormag, btot, dipang
!data-tiros
      use IDEA_ION_INPUT, only : NT_21, NT_20, NT_7, N_FLX, N_BND
      use IDEA_ION_INPUT, only : EMAPS, CMAPS, DJSPECTRA
      use IDEA_ION_INPUT, only : EMAPS1, CMAPS1, DJSPECTRA1
!data-imf
      use IDEA_IMF_INPUT, only :  IMF_read_wam_init, idea_imf_fix   
!
      implicit none
!SK   include 'mpif.h'
!
      integer          :: levs, me, master
      character(ch100) :: nc_file, ncfile_fpath, nctiros_fpath
      character(ch100) :: ncimf_fpath 
      integer          :: n3d, n2d
      logical          :: me_eq_master
!
      me_eq_master = me.eq.master
!
!SK    call MPI_COMM_RANK (MPI_COMM_WORLD, mpi_id, mpi_err)
!
! read_in NML_ION ( ION, TIROS & IMF data)
!
      if (me_eq_master)
     & print*,' 01:idea_ion_init-->ion_read_namelist;me= ',me
      CALL ion_read_namelist
!SK  &(nml_ion,nlun_ion,ncfile_fpath,nctiros_fpath,ncimf_fpath,mpi_id)
     &(nml_ion,nlun_ion,ncfile_fpath,nctiros_fpath,ncimf_fpath,me)
      if (me_eq_master)
     & print*,' 01:idea_ion_init<--ion_read_namelist;me= ',me
!
! READ-IMF data only for year of 2012 
! in /scratch3/NCEPDEV/swpc/save/Valery.Yudin/BASE_SVN/BASE_WAM_DATA/Solar_2012
! filename:  wam_nems_imf_2012_dp105410.nc
!
      if(me_eq_master)
     & print*,' 02:idea_ion_init->imf_read_wam_namelist;me=',me
       if (idea_imf_fix <=1 )  
     &     CALL IMF_read_wam_init(ncimf_fpath, me)
      if(me_eq_master)
     & print*,' 02:idea_ion_init<-imf_read_wam_namelist;me=',me

      if(me_eq_master)
     & print*,' 03:idea_ion_init->precomp_iondata_fixed;me=',me
           CALL precomp_iondata_fixed
      if(me_eq_master)
     & print*,' 03:idea_ion_init<-precomp_iondata_fixed;me=',me
!
!
!ion-2d-data replaces old  .....subroutine interp_field(ix,im,rlat,rlon,cormago,btoto,dipango)
!                               real cormag(20,91),btot(20,91),dipang(20,91),glat(91),glon(20)
!
      if (me_eq_master)
     & print*,' 04:idea_ion_init-->ion_read_wam_init;me= ',me
       CALL ion_read_wam_init(ncfile_fpath, me)       ! Tirosd/iondata_tjr.nc 
      if (me_eq_master)
     & print*,' 04:idea_ion_init<--ion_read_wam_init;me= ',me
!
       if (tiros_activity_fixnam > 0) then
!
!read nc or ascii formatted files
!
      if(me_eq_master)
     & print*,' 05:idea_ion_init->tiros_read_wam_init;me= ',me
          CALL tiros_read_wam_init(nctiros_fpath, me) ! Tiros/tiros_tjr.nc
      if(me_eq_master)
     & print*,' 05:idea_ion_init<-tiros_read_wam_init;me= ',me
!tiros-txt     call tiros_init(emaps,cmaps,djspectra)
      if(me_eq_master)
     & print*,' 06:idea_ion_init->tiros_init;me= ',me
          CALL tiros_init(emaps1,cmaps1,djspectra1)       ! Tiros ascci files
      if(me_eq_master)
     & print*,' 06:idea_ion_init<-tiros_init;me= ',me
       endif
!
       RETURN
!
! more fancy read on "selected" PE for input-files and do "mpi_bcast"
!
!SK     n2d  = nxmag*nymag
!SK     n3d  = NT_21*NT_20*NT_7   
!SK     call mpi_bcast(emaps  , n3d,    MPI_REAL8,0, MPI_COMM_ALL,info)
!       call mpi_bcast(cmaps  , n3d,    MPI_REAL8,0, MPI_COMM_ALL,info)
!       call mpi_bcast
!    &      (djspectra, N_FLX*N_BND, MPI_REAL8,0, MPI_COMM_ALL,info)
!       call mpi_bcast(cormag  , n2d,    MPI_REAL8,0, MPI_COMM_ALL,info)
!       call mpi_bcast(btot    , n2d,    MPI_REAL8,0, MPI_COMM_ALL,info)   
!       call mpi_bcast(dipang  , n2d,    MPI_REAL8,0, MPI_COMM_ALL,info)  
!       call mpi_bcast(glon, nxmag,    MPI_REAL8,0, MPI_COMM_ALL,info) 
!SK     call mpi_bcast(glat, nymag,    MPI_REAL8,0, MPI_COMM_ALL,info)    
!
 
      END SUBROUTINE IDEA_ION_INIT
!
!
      MODULE WAM_ION
!
!old         use idea_composition, f107 =>f107_idea, kp =>kp_idea
!old-wy      use wam_f107_kp_mod, only: f107, kp, kdt_3h
!      use idea_solar_input,  only : f107 => wf107_s, kp => wkp_s ....
!      the last "use" is deassembled to keep SWPC-F107/Kp forecasts activs

!SK   use module_physics_driver, only : mpi_id => mpi_me
!SK   use IDEA_MPI_def,      only : mpi_id
!
!      use wam_date_calendar, only : curday_wam, curmonth_wam, curddd_wam 
!      use wam_date_calendar, only : curyear_wam, curutsec_wam 
      use date_def, only          : idate                      ! substitute of WAM-calendar
      use idea_composition,  only : k91, ELCH
      use idea_composition,  only : PI,  PI2, DTR, R_2_D, fac_lst, PID2
      use idea_composition,  only : pi_24hr
!
!      use idea_ion_input, only : cormag(20,91),btot(20,91),dipang(20,91),glat(91),glon(20)
!
      CONTAINS
!
      subroutine idea_ion(mpi_id,pres,solhr,cospass,zg,grav,o_n,o2_n,
     &  n2_n,cp,adu,adv,adt,dudt,dvdt,dtdt,rho,rlat,rlon,ix,im,levs,
     &  dayno,utsec,sda,maglon,maglat,btot,dipang,essa,
     &  f107, f107d, kp, nhp, nhpi, shp, shpi, SPW_DRIVERS,
     &  swbz, swvel)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! driver      dtdt(i,k)=jh(i,k)/cp(i,k), dudt dvdt
!              ion darge and Joule heating
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      use namelist_wamphysics_def, only : JH0, JH_semiann,
     &                                    JH_ann, JH_tanh
      
      implicit none
      REAL  , INTENT(IN)   :: f107, f107d, kp  ! solar-geo inputs from different WAM applications
                                               ! CLIMATE-SWPC-RDATA, controlled in idea_phys
      REAL  , INTENT(IN)   :: nhp, nhpi, shp, shpi ! solar-geo inputs from wam_f107_kp.txt
      REAL  , INTENT(IN)   :: swbz, swvel ! solar-geo inputs from wam_f107_kp.txt

      Character,INTENT(IN) :: SPW_DRIVERS      ! SPACE weather/climate driver
!
!      REAL, PARAMETER :: DTR=3.141592653/180.0
      REAL, PARAMETER :: pi = 3.141592653
      INTEGER, INTENT(IN)     :: mpi_id
      INTEGER, INTENT(IN)     :: ix !longitude dim size
      INTEGER, INTENT(IN)     :: im !number of logitude
      INTEGER, INTENT(IN)     :: levs ! number of pres grid
      INTEGER, INTENT(IN)     :: dayno !calender day
 
      REAL, INTENT(IN)     :: pres(ix,levs) ! pressure, Pa

      REAL, INTENT(IN)     :: o_n(ix,levs)  ! number density O (/m3)
      REAL, INTENT(IN)     :: o2_n(ix,levs)
      REAL, INTENT(IN)     :: n2_n(ix,levs)
      REAL, INTENT(IN)     :: rho(ix,levs)  ! mass density (kg/m3)
      REAL, INTENT(IN)     :: zg(ix,levs)   !  height (m)
      REAL, INTENT(IN)     :: grav(ix,levs) !  gravity (m/s**2)
      REAL, INTENT(IN)     :: cp(ix,levs)   !  (J/kg/k)
      REAL, INTENT(IN)     :: cospass(im)! cos solar zenith angle (rad) 
      REAL, INTENT(IN)     :: rlat(im) ! latitude (rad)
      REAL, INTENT(IN)     :: rlon(im) ! longitude (rad)
      REAL, INTENT(IN)     :: solhr     ! universal time (h)
!  
! state
!
      REAL, INTENT(IN)     :: adt(ix,levs)  ! temperature (k)
      REAL, INTENT(IN)     :: adu(ix,levs)  ! zonal wind (m/s)
      REAL, INTENT(IN)     :: adv(ix,levs)  ! meridional wind (m/s)

! input Magnetic and electric parameters 
! 
      REAL, INTENT(in) :: maglon(im)  !magnetic longitude (rad)
      REAL, INTENT(in) :: maglat(im)  !magnetic latitude (rad)
      REAL, INTENT(in) :: btot(im)    !mapgnetic field strength
      REAL, INTENT(in) :: dipang(im)  !Dip angle (degree)
      REAL, INTENT(in) :: essa(im)    !magnetic local time
      REAL, INTENT(in) :: sda         ! solar declination angle (rad)
      REAL, INTENT(in) :: utsec       !universal time
! output
      REAL, INTENT(out)     :: dtdt(ix,levs)  ! temperature change (k/s)
      REAL, INTENT(out)     :: dudt(ix,levs)  ! zonal wind change (m/s2)
      REAL, INTENT(out)     :: dvdt(ix,levs)  ! meridional change wind (m/s2)
!
! local
!
! define jh_fac, storm_fac Zhuxiao Li 
      real rlt(im),sza(im),jh(ix,levs),rinc(5),jh_fac, st_fac, VBz

      INTEGER   i,k
      logical, save :: skprnt
!SK2020Oct5
!     skprnt = mpi_id.eq.0
      skprnt = .false.

! get VBz swvel*swbz

      VBz = swvel*swbz

! get sza in rad
      sza=acos(cospass)
! get local solar time in rad:  rlt=(rlon/(15.*pi/180.)+solhr)/24.*2.*pi
!                               rlt=( rlon*R_2_D/15. + solhr )*(2.*pi/24.)

       rlt =(solhr +rlon*fac_lst)*pi_24hr
!
! get ion_drag & joule heating
!
!
!      if (mpi_id == 0) then
!      print *, 'vay-pres-ion', maxval(pres)
!      print *, 'vay-rho-ion',  maxval(rho)
!      endif
      if (skprnt) then
       print*,' in idea_ion:max(pres)=',maxval(pres),' me=',mpi_id
       print*,' in idea_ion:max(rho)=',maxval(rho),' me=',mpi_id
       print*,' in idea_ion: --> GetIonParams, me=',mpi_id
      endif
      call GetIonParams(mpi_id,pres,
!    &   dayno,utsec,F107,KP,sda,sza,rlat,zg,grav,      
     &   dayno,utsec,F107,f107d,KP,NHP,NHPI,spw_drivers,sda,sza,rlat,zg,
     &   grav,o_n, o2_n, n2_n,adu,adv,adt,rho,rlt,rlon,ix,im,levs,k91,
     &   btot,dipang,maglon,maglat,essa,                                
     &   dudt,dvdt,jh) 
      if (skprnt) then
       print*,' in idea_ion: <-- GetIonParams, me=',mpi_id
       print*,' in idea_ion:F107= ',F107,' F107d= ',f107d,' me=',mpi_id
       print*,' in idea_ion:NHP= ',NHP,' NHPI= ',NHPI,' me=',mpi_id
      endif
!      if (mpi_id == 0) then
!       print *, 'F107=  ', F107, 'F107d=  ', f107d
!       print *, 'NHP=  ', NHP, 'NHPI=  ', NHPI
!      endif

! update to ........ K/sec
      do i=1,im

!   Joule heating factor to consider the seasonal variation and
!   semiannual variation, Zhuxiao.Li

!  JH0_5
!           jh_fac = 1.75+0.5*tanh(2.*rlat(i))*cos((dayno+9.)*2.*pi/365.)
!     &                  +0.5*(cos(4.*pi*(dayno-80.)/365.))

        jh_fac = JH0+JH_tanh*tanh(2.*rlat(i))*cos((dayno+9.)*2.*pi/365.)
     &              +JH_semiann*(cos(4.*pi*(dayno-80.)/365.))

! VBz adjustment

          if (abs(VBz).le.5000.) then
             st_fac = 1.
          else
             st_fac = (25000.+5000.)/(25000.+ abs(VBz))
          endif

      do k=1,levs
!           dtdt(i,k)=jh(i,k)*jh_fac/cp(i,k)
           dtdt(i,k)=jh(i,k)*jh_fac*st_fac/cp(i,k)
      enddo
      enddo
      return
!
      END SUBROUTINE idea_ion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE GetIonParams(me,pres, 
!     &   dayno,utsec,f107,kp,sda,sza,rlat,ht,grav, 
     &   dayno,utsec,f107,f107d,kp,hp,hpi,spw_drivers,sda,sza,rlat,ht,
     &   grav,o_n, o2_n, n2_n,adu,adv,adt,rho,rlt,rlon,ix,im,levs,lev1,
     &   btot,dipang,maglon,maglat,essa,                                
     &   dudt,dvdt,jh) 
!      use physcons,  pi => con_pi
!    
       use IDEA_ION_INPUT, only : EMAPS1, CMAPS1, DJSPECTRA1
       use IDEA_ION_INPUT, only : tiros_switch
       use IDEA_ION_INPUT, only : GW_fixnam, tiros_activity_fixnam 
       use idea_composition, only : pid12  
! 
      implicit none
      INTEGER,    INTENT(IN)     :: me
      INTEGER,    INTENT(IN)     :: dayno  !day 
      INTEGER,    INTENT(IN)     :: ix !longitude dim size
      INTEGER,    INTENT(IN)     :: im !number of logitude
      INTEGER,    INTENT(IN)     :: levs ! number of pres grid
      INTEGER,    INTENT(IN)     :: lev1 ! lowest pres level to start k91 from idea_comp...
!
!      REAL, PARAMETER :: DTR=3.141592653/180.0
!      REAL, PARAMETER :: ELCH=1.602e-19
      REAL, INTENT(IN)     :: pres(ix,levs) ! pressure, Pa
      REAL, INTENT(IN)     :: o_n(ix,levs)  ! number density O (/m3)
      REAL, INTENT(IN)     :: o2_n(ix,levs)
      REAL, INTENT(IN)     :: n2_n(ix,levs)
      REAL, INTENT(IN)     :: adt(ix,levs)  ! temperature (k)
      REAL, INTENT(IN)     :: adu(ix,levs)  ! zonal wind (m/s)
      REAL, INTENT(IN)     :: adv(ix,levs)  ! meridional wind (m/s)
      REAL, INTENT(IN)     :: rho(ix,levs)  ! mass density (kg/m3)
      REAL, INTENT(IN)     :: ht(ix,levs)   ! geopotential height (m)
      REAL, INTENT(IN)     :: grav(ix,levs) !  gravity (m/s**2)
      REAL, INTENT(IN)     :: f107, f107d, kp 
      REAL, INTENT(IN)     :: hp, hpi
      CHARACTER,INTENT(IN) :: spw_drivers
      REAL, INTENT(IN)     :: sda      ! solar declination angle (rad)
      REAL, INTENT(IN)     :: sza(im)  ! solar zenith angle (rad) 
      REAL, INTENT(IN)     :: rlt(im)  ! local time (rad) 
      REAL, INTENT(IN)     :: rlat(im) ! latitude (rad)
      REAL, INTENT(IN)     :: rlon(im) ! longitude (rad)
      REAL, INTENT(IN)     :: utsec    ! universal time (s)
!     REAL, INTENT(in) :: elx(im)
!     REAL, INTENT(in) :: ely(im)     !electric field
      REAL, INTENT(in) :: maglon(im)  !magnetic longitude (rad)
      REAL, INTENT(in) :: maglat(im)  !magnetic latitude (rad)
      REAL, INTENT(in) :: btot(im)    !mapgnetic field strength
      REAL, INTENT(in) :: dipang(im)  !Dip angle (degree)
      REAL, INTENT(in) :: essa(im)    !magnetic local time


      REAL, INTENT(OUT) :: dvdt(ix,levs) !(m/s2)
      REAL, INTENT(OUT) :: dudt(ix,levs) !(m/s2)
      REAL, INTENT(OUT) :: jh(ix,levs)   ! (J/kg/s)
! local
      real ht1(levs),v1(levs),nden(levs),o2n(levs),on(levs),            
     &    n2n(levs),elx(im),ely(im),ssa,elz(im),ee1(im),                
     &    ee2(im),cosdif,sindif,sdip,cdip,btheta,bphi,elecx,            
     &    elecy,dif,dlat,dlon
      INTEGER   k,i,n,tiros_activity_level  
!     Ion drag variables :
! 
!     teff(levs)      1d local array of temperature
!     pion1(levs)     number density O+
!     pion2(levs)     number desntiy NO+
!     pion3(levs)     number density O2+
!     r                
!     sigped           pedersen conductivity
!     sighall          hall conductivity
!     jphi(levs)     eastward curreil
!     jth(levs)      southward
!     rvin(levs)     Ion/Neutral collision  frequency param    
!     ramin(levs)    mean ion mass
!     a5               Meridional ion drag term
!     b5               Zonal ion drag term
!     c7               Joule heating term
!     
      REAL      :: teff(levs), pion1(levs), pion2(levs)
      REAL      :: sigped, sighal, pion3(levs)
      REAL      :: rvin(levs), ramin(levs)
      REAL      :: r, brad, bth, dip
      REAL      :: a5, b5, c7
      REAL      :: jth,jrad  
      REAL      :: jphi,   GW 
      REAL      :: eden(ix,levs)  !electron density
      REAL      :: eden_chiu(ix,levs)  !electron density from CHIU
      REAL      :: eden_aurora(ix,levs)  !electron density from TIROS
!
      real,dimension(levs) :: pr1,rho1,grav1, adt1
      real,dimension(levs) :: o3p_1,o2_1,n2_1,aur_1
!      REAL                 :: eflux, ch
!      real      :: emaps1(21,20,7),cmaps1(21,20,7),djspectra1(15,21)
!SK
      logical :: skprnt
      skprnt = me.eq.0
!  
!===================================================================
!         Calculate Electric Field and magnetic field         
!===================================================================
! VAY-2016: (im, ix) => (im)
!
      if(skprnt) print*,' in GetIonParams-->idea_geteb'
      call idea_geteb(im, dayno,utsec,f107,f107d,kp,maglat,maglon,
     &     essa,ee1,ee2)
      if(skprnt) print*,' in GetIonParams<--idea_geteb'
!     ee1=0.
!     ee2=0.
! ===================================================================
! =                   Calculate Electron Density                    =
! ===================================================================
! CHIU ionosphere for electron density (Earth_chiu_model.f90). 
! vay-out-NOV1        CALL tiros_init(emaps,cmaps,djspectra)
!
!       print *, 'tiros_switch  ', tiros_switch
!       if(tiros_switch ==1)  CALL tiros_init(emaps1,cmaps1,djspectra1)
!
! set the power index 1 to 10 eventually read in
! set the number of gigawatts of auroral power (used if activity level = 10)
!
      IF (trim(SPW_DRIVERS)=='swpc_fst') then
         tiros_activity_level = hpi
         GW = hp
      ELSE          
         tiros_activity_level = tiros_activity_fixnam
         GW = GW_fixnam
      ENDIF
!
      DO i = 1,im 
!
         do k=1,levs
           ht1(k)=ht(i,k)
         enddo
!
      if(skprnt.and.i.eq.1)print*,' in GetIonParams-->EARTH_CHIU_MODEL'
         CALL EARTH_CHIU_MODEL(sda,sza(i),maglat(i),                    
     &         maglon(i),rlt(i), rlat(i), f107,                         
     &         dipang(i)*DTR, dayno, ht1, eden_chiu,i,lev1,             
     &         levs,ix)
      if(skprnt.and.i.eq.1)print*,' in GetIonParams<--EARTH_CHIU_MODEL'
!
! tiros_ionize returns ionization rates for O, O2, and N2 for a given
! geomagnetic latitude GL and magnetic local time MLT based on
! TIROS/NOAA statistical maps of energy influx and characteristic energy
! the ionization rates assumes an observed spectrum based on the
! characteristic energy CH at the location, the TIROS spectrum is
! between 300eV - 100keV, below 300eV assumes a Maxwellian where
! the average energy of the distribution is CH
! the model atmosphere should be provided by calling program (e.g., IPE)
! a sample model profile is read in from ipe_nh1_data and ipe_nh2_data
!
! input
! geomagnetic latitude mlat in radians
! magnetic hour angle from noon essa in degrees
!      essa = essa_in
! output eden_aurora units number/m3   ion density from aurora
!
!
      if(skprnt.and.i.eq.1)print*,
     &  ' in GetIonParams:tiros_switch=',tiros_switch
   
      if(tiros_switch ==0) THEN
        do k=1, levs
        pr1(k) =pres(i,k)
        rho1(k)=rho(i,k)
        grav1(k) = grav(i,k)
        adt1(k) = adt(i,k)
        o3p_1(k) =o_n(i,k)
        o2_1(k) =o2_n(i,k)
        n2_1(k)=n2_n(i,k)
        enddo
      if(skprnt.and.i.eq.1) print*,' in GetIonParams-->tiros_ionize'
      call tiros_ionize(me,lev1,levs, pr1, rho1, ht1,
     &     grav1,o3p_1,o2_1, n2_1,adt1,maglat(i),essa(i),
     &     tiros_activity_level, GW, aur_1)
      if(skprnt.and.i.eq.1) print*,' in GetIonParams<--tiros_ionize'
        do k=1, levs
         eden_aurora(i,k) = aur_1(k)
        enddo
!
      ELSE IF (tiros_switch ==1) THEN
      if(skprnt.and.i.eq.1) print*,
     & ' in GetIonParams-->tiros_ionize_data'
      call tiros_ionize_data(me,pres(i,:), lev1,levs,ht1,emaps1,cmaps1,
     &      djspectra1, grav(i,:),o_n(i,:),o2_n(i,:),
     &      n2_n(i,:),adt(i,:),maglat(i),essa(i),tiros_activity_level,
     &      GW, eden_aurora(i,:) )    !don't use, eflux, ch)
      if(skprnt.and.i.eq.1) print*,
     &   ' in GetIonParams<--tiros_ionize_data'
      ENDIF
!
!
      eden(i,:)=sqrt(eden_chiu(i,:)**2   + eden_aurora(i,:)**2)
      ENDDO
!      if (mpi_id == 0) then
!      print*,'chiuok'
!      print*,'eden',eden(1,lev1:levs)
!      ENDIF

!=================================================================
!         Calculate Ion Drag, Joule Heating and Particle        =
!                     Precipitation terms                       =
!=================================================================
      DO i = 1, im 
!======================================Horizontal dep-nt loops
!  
          dip     = dipang(i)*DTR
          sdip    =sin(dip)                              !new
          cdip    =cos(dip)                              !new
!
          elecx   =ee2(i)*sdip
          elecy   =ee1(i)
!
          ssa=rlon(i)+(utsec/3600.-12.)*pid12              ! radians  ????? for South Pole
          dif     =essa(i)*DTR -   ssa                     !check unit
          cosdif  =cos(dif)
          sindif  =sin(dif)
!
          elz(i)  =-1.*ee2(i)*cdip
          if(sdip.ge.0.) then
          elx(i)  =elecx*cosdif-elecy*sindif
          ely(i)  =elecx*sindif+elecy*cosdif
          else
          elx(i)  =elecx*cosdif+elecy*sindif
          ely(i)  =-1.*elecx*sindif+elecy*cosdif
          endif
!
!
          btheta  =-1.*btot(i)*cdip*cosdif               !new
          bphi    =-1.*btot(i)*cdip*sindif               !new
          bth     = btot(i)                              ! In teslas, so no *1.e-9
          brad    = -1.*bth*sdip
!====================================== Vertical dep-nt loops
        do k = lev1,levs
          nden(k)=o_n(i,k)+o2_n(i,k)+n2_n(i,k)
          teff(k) = adt(i,k)
          v1(k)=           -adv(i,k)      ! v1 positive south
          on(k)=o_n(i,k)
          o2n(k)=o2_n(i,k)
          n2n(k)=n2_n(i,k)
        enddo
        do k=lev1,levs
!          dudt(i,k)  = 0.
!          dvdt(i,k)  = 0.
!          jh(i,k) = 0.
          pion1(k) = on(k)
          pion2(k) = 0.5*(o2n(k)+n2n(k)) 
          pion3(k) = 0.5*(o2n(k)+n2n(k))
         enddo
! Get ion neutral collision frequency
      if(skprnt.and.i.eq.1) print*,' in GetIonParams-->IONNEUT'
      CALL IONNEUT(on,o2n,n2n,pion1,pion2,pion3,
     &             teff,rvin,ramin,levs,lev1)
      if(skprnt.and.i.eq.1) print*,' in GetIonParams<--IONNEUT'
!     print*,'ionneut ok'
! Calculate ion drag and electron deposition
! jth                     - N/S electrical conductivity
! jphi                    - E/W electrical conductivity

        do k=lev1,levs 
          r       = (ramin(k)*rvin(k))/(ELCH*bth)
          sigped  = (eden(i,k)*ELCH*r)/(bth*(1.0+r**2))
          sighal  = sigped*r
          jphi = sigped*(ely(i)-v1(k)*brad)                             
     &         - sighal*(elx(i) + adu(i,k)*brad)*sdip         !new fix of 201612=> /sdip => *sdip
          jth  = sigped*(elx(i) + adu(i,k)*brad) +                      
     &         sighal*(ely(i)-v1(k)*brad)*sdip                !new
          jrad = sigped*(elz(i) - adu(i,k)*btheta+v1(k)*bphi)           
     &        -sighal*(ely(i)-v1(k)*brad)*cdip                !new
          a5 =(jphi*brad-jrad*bphi)/rho(i,k)                  !new
             b5 = -1.*(jth*brad-jrad*btheta)/rho(i,k)         !new
          c7 =(jth*(elx(i)+adu(i,k)*brad)+                              
     &         jphi*(ely(i)-v1(k)*brad)+jrad*elz(i))/rho(i,k) !new
! Calculation of ion drag terms END
          dvdt(i,k)   = -1.* a5
          dudt(i,k)   =  b5
          jh(i,k)   = c7
        enddo
        do k=1,lev1-1 
          dvdt(i,k) = 0.
          dudt(i,k) = 0.
          jh(i,k) = 0.
        enddo
      ENDDO                          ! i-hor.-pixel loop
!
!        if (mpi_id == 0) then
!           print *, 'vay-ion GetIonParams'
!        endif
!
      if(skprnt) print*,' in GetIonParams: RETURN'
      RETURN
      END SUBROUTINE GetIonParams
!
      SUBROUTINE IONNEUT(P1,P2,P3,PI1,PI2,PI3, T,VIN,A_MIn,NMAx,n0)
      INTEGER :: NMAx, n0
      REAL    ::  P1(NMAx) , P2(NMAx) , P3(NMAx) 
      REAL    ::   T(NMAx) , VIN(NMAx) ,  A_MIn(NMAx)
      REAL    :: PI1(NMAx) , PI2(NMAx), PI3(NMAx)
!
!vay-2015 : sum => sum_vay
!           Amin => A_min , ftiny ~ 1.e-90                                                          
      REAL    :: sum_vay ,  v1 , v2 
      INTEGER :: n
      REAL    ::     a(3) , b(3) 
      REAL, parameter :: amu=1.66E-27  
      REAL, parameter :: factor=1.0
      REAL, parameter :: ftiny=1.e-90
      REAL    :: mi1 , mi2, mi3, summol
      REAL    ::  logtn    
      DATA mi1 , mi2, mi3/16. , 30., 32./
!********************************************************************
! The following a,b, are cooeficients used to caculate ion-neutral
! collision frequency. Tim's Thesis  3.5a,3.5b.  mjh 1.9.97
!********************************************************************
      
      DATA a/3.42E-11 , 6.66E-10 , 6.82E-10/
      DATA b/2.44E-10 , 4.28E-10 , 4.34E-10/
      
      DO 100 n = n0 , NMAx
         summol = PI2(n) + PI3(n)
         if(summol.lt.ftiny) summol=0.0
         sum_vay = PI1(n) + PI2(n) + PI3(n)
         v2 = b(1)*P1(n) + b(2)*P2(n) + b(3)*P3(n)
         logtn =LOG10(T(n))
         v1 = a(3)*P3(n) + a(2)*P2(n) + a(1)*P1(n)*factor*SQRT(T(n))    
     &       *(1.08 - 0.139*LOGTN+4.51E-03*LOGTN*LOGTN)
         if(v1.lt.ftiny) v1=0.0    ! v1 = max(v1, ftiny)
         if(v2.lt.ftiny) v2=0.0    ! v1 = max(v2, ftiny)
!        if(pi1(n).lt.1.e-90) pi1(n)=0.0
!     if(iout.eq.1) write(6,*) 'here 5',n
         VIN(n) = (v1*PI1(n)+v2*summol)*1.E-06/sum_vay
         A_MIn(n) = (PI1(n)*mi1+PI2(n)*mi2+PI3(n)*mi3)*amu/sum_vay
 100  CONTINUE

      RETURN
      END SUBROUTINE IONNEUT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TIROS empirical Model


! idea
!
      subroutine idea_geteb(im, dayno,utsec,f107,f107a,kp,maglat,maglon,
     &   essa,ee1,ee2)
      use efield_wam      !  iday,iyear,iday_m,imo,f107d,by,bz,ut,v_sw     
!
! vay-2016: another interpolation for the Electric fields
!                           to model-based maglon-maglat  
!     
      use idea_imf_input, only : idea_imf_fix
      use idea_imf_input, only : bz_s, by_s, swvel_s
!
!     use idea_imf_input, only : bz_fix, by_fix, swvel_fix (Fixme: undefined)
! 
      use date_def   
!
      implicit none
      integer, intent(in) :: im  ! number of data points in efield 
      integer, intent(in) :: dayno  ! calender day
      real, intent(in)    :: utsec  ! second
      real, intent(in)    :: f107, f107a  ! 
      real, intent(in)    :: kp  ! 
      real, intent(in)    :: maglat(im)  ! magnetic latitude (rad)
      real, intent(in)    :: maglon(im)  ! magnetic longitude (rad)
      real, intent(in)    :: essa(im)    ! degree
      real, intent(out)   :: ee1(im)     ! electric field x direction mV/m
      real, intent(out)   :: ee2(im)     ! electric field y direction mV/m
!
!     character*(*), intent(in) ::   dir    ! directory located coef files
! local
      integer ::   i,k,iref,jref

      real  :: utsec_last
      real  :: dx,dy,aa,bb,maglond,maglatd,                      
     &    ed11(0:nmlon,0:nmlat),ed22(0:nmlon,0:nmlat)
      real  :: v_sw
!
      data utsec_last/-1./
      save utsec_last,ed11,ed22
!
! initiate
! calculate efield only if diff time step
!
      print*,' in idea_geteb: utsec,utsec_last', utsec, utsec_last
      if(utsec.ne.utsec_last) then
        utsec_last=utsec
!     use f107d directly from observation sheet, revised by Tzu-Wei and
!     Zhuxiao 
        f107d=f107a
        ut=utsec/3600.
        iday = dayno                   ! day of year
        imo=idate(2)
        iday_m=idate(3) 
      print*,' in idea_geteb: 001'
! wam-calendar
!     curday_wam, curmonth_wam, curddd_wam 
! 
!        iday   =  curddd_wam
!        imo    =  curmonth_wam
!        iday_m =  curday_wam
!        iyear  =  curyear_wam 
!       print *, 'VAY-geteb ', iday,iyear,iday_m,imo, ut
!
! use efield         !  iday,iyear,iday_m,imo, ut, f107d, by, bz, v_sw 
! Para-n for bz       as default setting: by, bz, v_sw 
!
        bz = .433726 - kp*(.0849999*kp + .0810363)                      
     &        + f107d*(.00793738 - .00219316*kp)
        by=0.
        v_sw = 450.00
!
!Time-dependent IMF only for year of 2012
! 
!        if (idea_imf_fix == 0) then
!          bz = bz_s
!          by = by_s
!          v_sw = swvel_s      
!        endif

! Tzu-Wei: print values to make sure they are right
!       bz= -8.
!       print *,'check IMF Bz=',bz,'Kp=',kp,'f107=',f107d
!================================================   call from "module efield"
      print*,' in idea_geteb: --> get_efield'
        CALL  get_efield
      print*,' in idea_geteb: <-- get_efield'

!       print*,'www'
!       print'(8f10.4)',potent(0:180,68)
        ed11=ed1
        ed22=ed2
!     print*,'ed2',ed2(149,65)
      endif
!
      do k=1,im
        maglatd=maglat(k)*R_2_D
!
        jref=0
        dy=0.0
!efield.f:     &         ylonm, ylatm   ! magnetic longitudes/latitudes (degc .....ylatm(0:nmlat)ylonm(0:nmlon)
      do i=0,nmlat-1
!         if(maglatd.ge.ylatm1(i)-90..and.maglatd.le.ylatm1(i+1)-90.)   
      IF(maglatd.ge.ylatm (i)-90..and.maglatd.le.ylatm (i+1)-90.) then
            jref=i
!hmhj       dy=(maglatd-ylatm1(i)+90.)/(ylatm1(i+1)-ylatm1(i))
            dy=(maglatd-ylatm (i)+90.)/(ylatm (i+1)-ylatm (i))
      endif
      enddo 
!       print*,'wwwlat',k,maglatd,jref,dy
!       maglond=maglon(k)/pi*180.
        maglond=essa(k)+180.
        if(maglond.lt.0.)   maglond=maglond+360.
        if(maglond.gt.360.) maglond=maglond-360.
! ylonm = [0,360]
        iref=0
        dx=0.0
        do i=0,nmlon-1
          if(maglond.ge.ylonm (i).and.maglond.le.ylonm (i+1)) then
            iref=i
            dx=(maglond-ylonm (i))/(ylonm (i+1)-ylonm (i))
          endif
        enddo 
!       print*,'wwwlon',k,maglond,iref,dx
        aa=(1.-dx)*ed11(iref,jref)+dx*ed11(iref+1,jref)
        bb=(1.-dx)*ed11(iref,jref+1)+dx*ed11(iref+1,jref+1)
        ee1(k)=(1.-dy)*aa+dy*bb
        aa=(1.-dx)*ed22(iref,jref)+dx*ed22(iref+1,jref)
        bb=(1.-dx)*ed22(iref,jref+1)+dx*ed22(iref+1,jref+1)
        ee2(k)=(1.-dy)*aa+dy*bb
!       if(ely(k).gt.100.) print*,'ely',utsec,ed22(iref,jref),          
!    &ed22(iref+1,jref),ed22(iref+1,jref),ed22(iref+1,jref+1),          
!    &maglond,maglatd,iref,jref
!       if(ely(k).gt.100.) ely(k)=0.
      enddo
!     ee1=1000.*ee1
!     ee2=1000.*ee2
! correct direction? 365.25day?
!     print*,' in idea_geteb: RETURN'
      return
      end subroutine idea_geteb
!
!r===========================================================
!r=           Earth Electric and Magnetic Field
!r===========================================================
!r
      SUBROUTINE getmag(im,utsec,rlat,rlon,sda,                      
     &btot,dipang,maglon,maglat,essa)
!
! Oct 2016 take out IX from all 1D-subs
!
!      use physcons, pi => con_pi
      IMPLICIT NONE

!     REAL(prcn) high_lat_limit  Limit in degrees above 
!     which foster used. Below this limit, Richmond 
!     field used
!      real, PARAMETER ::high_lat_limit=60.
!      real, PARAMETER:: DTR = pi/180.
!      real, PARAMETER:: ELCH = 1.062e-19
! Input parameters

      INTEGER,    INTENT(IN)     :: im  !number of longitude
      REAL,       INTENT(IN)     :: utsec !UT second
      REAL,       INTENT(IN)     :: rlat(im)  ! geo latitude (rad)
      REAL,       INTENT(IN)     :: rlon(im)  ! geo longitude (rad)
      REAL,       INTENT(IN)     :: sda  ! solar diclination angle (rad)
! Output Magnetic and electric parameters 
!     REAL, INTENT(OUT) :: elx(im)
!     REAL, INTENT(OUT) :: ely(im)     !electric field
!
      REAL, INTENT(OUT) :: maglon(im)  !magnetic longitude (rad)
      REAL, INTENT(OUT) :: maglat(im)  !magnetic latitude (rad)
      REAL, INTENT(OUT) :: btot(im)    !mapgnetic field strength
      REAL, INTENT(OUT) :: dipang(im)  !Dip angle (degree)
      REAL, INTENT(OUT) :: essa(im)    !magnetic local time
! Local 
      real cormag(im),cmorg(im)
      integer i
! set elx ely zero first
!     elx=0.
!     ely=0.
! get cormag btot dipang in grid
      call interp_field(im,rlat,rlon,cormag,btot,dipang)
! get maglon,maglat
      call SPOLE(im,RLAT,rlon,utsec,SDA,maglon,ESSA,CMORG)

      do i=1,im
         maglat(i)=pid2-cormag(i)*DTR
!         print *, i, maglat(i)*R_2_D, 'VAY_MAG_deg ' 
      enddo
      return
      end SUBROUTINE getmag
!      
      SUBROUTINE SPOLE(im,RLAT,PHIR,utsec,SDA,PHIMR,ESSA,CMORG)
      implicit none
      real, PARAMETER    ::PI=3.141592653,DTR=PI/180.
      integer,intent(in) :: im          ! number of longitude 
      real,   intent(in) :: rlat(im)        !geo latitude (rad)
      real,   intent(in) :: phir(im)    !geo longitude (rad)
      real,   intent(in) :: utsec       !UT second
      real,   intent(in) :: sda         !solar declination angle (rad)
      real,   intent(out):: phimr(im)   !maglongitude (rad)        
      real,   intent(out):: essa(im)    !magnetic local time    
      real,   intent(out):: cmorg(im)          
! local variables
      real th,th1,phi1,sinth,sinth1,costh1,sinph1,cosph1,ac1,bc1,cc1,   
     & ac2,bc2,cc2,phim,ssp,sspr,csda,as1,bs1,cs1,as2,bs2,cs2,gml,      
     & cmag
      integer i
!
      do i=1,im
      th=pi/2.-rlat(i)
!
! SET POLE COORD. FOR EACH HEMIS.
!
      IF (RLAT(i).GE.0.0) THEN
       TH1=9.25*DTR
       PHI1=-78.0*DTR
      ELSE
       TH1=16.32*DTR
       PHI1=-54.0*DTR
      END IF
!
      SINTH=SIN(TH)
      SINTH1=SIN(TH1)
      COSTH1=COS(TH1)
      SINPH1=SIN(PHI1)
      COSPH1=COS(PHI1)
!
!     do i=1,im
      AC1=SINTH*COS(PHIR(i))
      BC1=SINTH*SIN(PHIR(i))
      CC1=COS(TH)
      AC2=AC1*COSTH1*COSPH1+BC1*COSTH1*SINPH1-CC1*SINTH1
      IF((ABS(AC2)).LT.0.001)AC2=0.001
      BC2=-AC1*SINPH1+BC1*COSPH1
      CC2=AC1*SINTH1*COSPH1+BC1*SINTH1*SINPH1+CC1*COSTH1
      CMORG(i)=ACOS(CC2)
      PHIMR(i)=ATAN2(BC2,AC2)
      PHIM=PHIMR(i)/DTR
!     SSP=360.-utsec/240.
      SSP=180.-utsec/240.
      SSPR=SSP*DTR
      CSDA=PI/2.-SDA
      AS1=COS(SSPR)*SIN(CSDA)
      BS1=SIN(SSPR)*SIN(CSDA)
      CS1=COS(CSDA)
      AS2=AS1*COSTH1*COSPH1+BS1*COSTH1*SINPH1-CS1*SINTH1
      IF((ABS(AS2)).LT.0.001)AS2=0.001
      BS2=-AS1*SINPH1+BS1*COSPH1
      CS2=AS1*SINTH1*COSPH1+BS1*SINTH1*SINPH1+CS1*COSTH1
      GML=ATAN2(BS2,AS2)/DTR
      ESSA(i)=PHIM-GML
      enddo
      RETURN
      END SUBROUTINE SPOLE
!
      subroutine interp_field(im,rlat,rlon,cormago,btoto,dipango)
!
!  VAY DANGER !!!!!!
! interp works only for [20,91]  46 center + fixed ddlat=180/90? and ddlon=360/20?
!
!
      USE IDEA_ION_INPUT, only :
     & cormag, btot, dipang, glat, glon, nxmag,nymag
!
      implicit none
      
      integer,intent(in)  :: im            ! number of longitude
      real,   intent(in)  :: rlat(im)      ! latitude (rad)
      real,   intent(in)  :: rlon(im)      ! longitude (rad)
!
      real,   intent(out) :: cormago(im),btoto(im),dipango(im) ! field value 
!
! local variable
!      real cormag(20,91),btot(20,91),dipang(20,91),glat(91),glon(20)
      real dll,dl,ddlat,ddlon,a1,a2,b1,b2,aa,bb
      integer i,iref,jref,jref1
      integer ::  jcen, ixdim 
!
! lat lon interval
      ddlat= 3.4906585033333331E-002
      ddlon= 0.3141592653000000
      jcen =46 
       ixdim =20    
! 
      do i=1,im
! get latitude index
        iref=int(rlat(i)/ddlat)+  jcen
        dl=(rlat(i)-glat(iref))/ddlat
! print*,iref,dl
! get longitude index
        jref=int(rlon(i)/ddlon)+1
        jref1=jref+1
        if(jref1.gt.ixdim) jref1=jref1 - ixdim
        dll=(rlon(i)-glon(jref))/ddlon
! print*,i,jref,jref1,dll
!
        a1=cormag(jref,iref)
        a2=cormag(jref1,iref)
        b1=cormag(jref,iref+1)
        b2=cormag(jref1,iref+1)
        aa=(1.-dll)*a1+dll*a2
        bb=(1.-dll)*b1+dll*b2
        cormago(i)=(1.-dl)*aa+dl*bb
!
        a1=btot(jref,iref)
        a2=btot(jref1,iref)
        b1=btot(jref,iref+1)
        b2=btot(jref1,iref+1)
        aa=(1.-dll)*a1+dll*a2
        bb=(1.-dll)*b1+dll*b2
        btoto(i)=(1.-dl)*aa+dl*bb
!
        a1=dipang(jref,iref)
        a2=dipang(jref1,iref)
        b1=dipang(jref,iref+1)
        b2=dipang(jref1,iref+1)
        aa=(1.-dll)*a1+dll*a2
        bb=(1.-dll)*b1+dll*b2
        dipango(i)=(1.-dl)*aa+dl*bb
!
      enddo
      return
      end subroutine interp_field

      END MODULE WAM_ION
!
