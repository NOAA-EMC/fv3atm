!**  $Id: chiu_model.f90,v 1.1.1.1 2006/06/04 18:19:13 cwplot Exp $
!r
!r   chui_model.f       Chiu ionosphere, to return electron density.
!r                      Converted to run F90, but not changed. mjh
!r   tiros_ionize_data  added the scaling of qiont with large GW and high tiros_activity_level
!r                      Zhuxiao Li, 8/31/2017 
!r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE Earth_CHIU_MODEL(sda,sza,thmag,phimr,rlt,rlat,         
     &f107, dip, nday,ht1d,eden3d,ilon,lev1,ht_dim,lon_dim)
!
!VAY-2016, cosmetic corrections for IMPLICIT none (ty and ilon)
! +        programming style with sin & cos
!
      use idea_composition, only : pi, dtr, pid12, pid6, RYDAYS,
     &         pid2, pid3, pid4, pid9, pid18 
      IMPLICIT NONE
!
      REAL,     INTENT(IN)   :: sda, sza, thmag, phimr, rlt
      REAL,     INTENT(IN)   :: rlat, f107, dip
      INTEGER,  INTENT(IN)   :: nday,ht_dim,lon_dim
      INTEGER,  INTENT(IN)   :: lev1                    ! first level of ion-re k91, see idea_composition
      REAL,     INTENT(IN)   :: ht1d(ht_dim)
! out  
      REAL,     INTENT(OUT)  :: eden3d(lon_dim,ht_dim)  ! should 1D-array no neeeds for ilon-index
!locals
      REAL                   :: ty
      INTEGER                :: ilon, i, n
!
!*** Start of declarations inserted by SPAG
      REAL :: abstmg , beta , cbp , cosrlt , costmg , cosza , dipf ,    
     &    e , f , flong , g , g5 , g6 , g7 , g8 , gel , gel1 ,    
     &    gsm,  a, rh
      REAL :: P , pb , qel , rd , rgamma , RHO , rk , rl ,         
     &    rt , s , sap ,sintmg 
      REAL ::  ty1 , ty2 , u , V , w , wr , x , y, alp,rr,fz, fn
 
!*** End of declarations inserted by SPAG
  
!- define parameters

  
      DIMENSION f(3) , pb(3) , s(3) , rd(3) , rl(3) , rt(3) , e(3) ,    
     & u(3) , V(3) , P(3) , flong(3) , dipf(3), alp(3), a(3)            
     & ,rh(3), rr(3), fz(3), fn(3)
  
      REAL :: z(ht_dim)
  
!  
! absolutely no idea what these are. Imported from tucan.f ????
!
      DATA alp/.5 , .5 , 1./
      DATA p/110. , 180. , 0./
      DATA a/1.36 , 2.44 , 0.66/
      DATA rh/10. , 34. , 0./
  
      rho = (f107-50.)*0.01
      ty = (nday+15.5)*12.*RYDAYS       !/365.
      IF ( ty > 12.0 ) ty = ty - 12.0
  
      abstmg = ABS(THMag)
      cosza =  COS(SZA)
      sintmg = SIN(THMag)
      costmg = COS(THMag)
      cosrlt = COS(RLT)


      ty1 = SIN(PId12*TY)
      ty2 = COS(PId6 *TY)
      P(1) = 110.
      P(2) = 180.
      f(1) = 0.0
      f(2) = 0.0
      pb(1) = 1.0
      pb(2) = 1.0
      s(1) = SQRT(1.0+1.15*RHO)
      s(2) = SQRT(1.0+1.24*RHO+0.25*RHO*RHO)
      rl(1) = 1.0
      rl(2) = 1.0
      e(1) = 1.0
      e(2) = 1.0
      flong(1) = 1.0
      flong(2) = 1.0
      dipf(1) = 1.0
      dipf(2) = 1.0
      g5 = SIN(PHImr)
      g6 = SIN(PHImr*.5)
      g7 = SQRT(ABS(g5))
      g8 = COS(.5*(PHImr-PI))
      sap = SIN(SDA)*sintmg
      f(3) = EXP(-(2.92*SIN(PId2-abstmg))**6)
      IF ( THMag <= 0.0 ) THEN
         cbp = 0.0
         IF ( g7 /= 0.0 ) cbp = ty1*(0.5*g6-0.5*g5-g6**8)-(1.0+ty1)     
     &    *ty2*g5/g7*EXP(-4.0*g6*g6)
         pb(3) = (2.5+2.0*RHO+ty2*(0.5+(1.3+0.5*RHO)*g8**4)             
     &    +(1.3+0.5*RHO)*COS(RLT-PI*(1.0+cbp)))                         
     &    *(1.0+0.4*ty1*ty1*EXP(-ty1*g8**4))
      ELSE
        wr = EXP(-1.2*(COS(THMag-DTR*23.5*cosrlt)-costmg))
        pb(3) = (2.0+1.0*RHO)*wr*(1.0+0.3*ty1)
      ENDIF
      s(3) = (1.0+RHO+0.204*RHO**2+0.05*RHO**3)
      IF ( RHO > 1.1 ) s(3) = 2.41 + 1.53*(s(3)-2.41)*(sintmg)**2
      P(3) = 240 + 75.0*RHO + 83.0*RHO*sap*costmg +                     
     & 30.0*COS(RLT-4.5*ABS(THMag)-PI)                                  
     & + 10.0*costmg*COS(PId3*(TY-4.5))
      rd(1) = EXP(2.0*(cosza/ABS(cosza)*SQRT(ABS(cosza))-1.0))
      rd(2) = EXP((1.0+0.5*LOG(1.0+30.0*RHO))*(cosza/ABS(cosza)*SQRT(ABS 
     & (cosza))-1.0))
      rd(3) = (0.9+0.32*sap)*(1.0+sap*(COS(RLT+PId4))**2)             
     & *EXP(-1.1*(1.0+COS(RLT-0.873)))
      qel = 1.0 - 0.15*EXP(-SQRT((12.0*THMag+1.05)**2+(TY/2.0-3.0)**2))
      rl(3) = (1.2-0.5*(costmg)**2)                                     
     & *(1.0+0.05*RHO*(sintmg)**3*COS(PId6*TY))                       
     & *(EXP(3.0*COS(0.5*THMag*(SIN(RLT)-1.0))))*qel
      w = COS(RLAt+SDA*cosrlt) - COS(RLAt)
      rt(1) = EXP(-0.4*w)
      rt(2) = EXP(-0.25*w)
      beta = 1.3 + 0.139*RHO**2 + 0.009*RHO**3
      rk = 1.0 + 0.085*(COS(THMag-PId6)*(COS(PId12*(TY-2.0)))       
     & **3+COS(THMag+PId4)*(COS(PId12*(TY-8.0)))**2)
      x = 0.7*(rk+0.178*RHO**2/s(3)*COS(PId3*(TY-4.3)))               
     & *EXP(-beta*(COS(THMag+SDA*cosrlt)-costmg))
      y = 0.2*(1.0-SIN(abstmg-0.524))*(1.0+0.6*COS(PId3*(TY-3.94)))   
     & *COS(PId6*(TY-1.0)) + (0.13-0.06*SIN(ABS(abstmg-PId9)))      
     & *COS(PId3*(TY-4.5)) - (0.15+0.3*SIN(abstmg))*(1.-cosrlt)        
     & **0.25*(COS(THMag+SDA))**3
      rt(3) = x + y/s(3)
      g = (1.0+0.6*SQRT(RHO)-0.2*RHO)*EXP(0.25*(1.0+COS(RLT-4.01)))
      gel = (costmg)**8*(COS(abstmg-0.262))**12
      gel1 = 1.0 + 0.05*(0.5-COS(PId3*TY)+COS(PId6*TY))
      e(3) = (1.0-0.4*(costmg)**10)                                     
     & *(1.0+0.6*(costmg)**10*(COS(RLT+PId4))**2)*(1.0+g*gel)         
     & *gel1
      rgamma = 1.0 + 0.03*(0.5-COS(PId3*TY)+COS(PId6*TY))
      gsm=0.15-(1.0+RHO)*(SIN(THMag/2.0))**2*EXP(-0.33*(TY-6.0)**2)
      flong(3) = 1.0 + 0.1*(costmg)**3*COS(2.0*(PHImr-7.0*PId18))
      dipf(3)=rgamma*(1.0+gsm*EXP(-18.0*(ABS(DIP)-2.0*PId9)**2))
      DO i = 1 , 3
         u(i) = s(i)*rd(i)*rl(i)*rt(i)*e(i)*flong(i)*dipf(i)
         V(i) = f(i)*pb(i) + (1.0-f(i))*u(i)
      enddo
!
! vertica loop
!
      DO n = lev1 , ht_dim
         z(n) = ht1d(n)/1000.
         IF ( z(n) <= p(3) ) rh(3) = 2.0*(20.0+0.1*z(n))
         IF ( z(n) > p(3)  ) rh(3) = 2.0*(20.0+0.1*p(3))
      DO i = 1 , 3
        rr(i) = (z(n)-p(i))/rh(i)
        fz(i) = EXP(alp(i)*(1.0-rr(i)-EXP(-rr(i))))
        fn(i) = a(i)*fz(i)*v(i)
      enddo
       eden3d(ilon,n) = (fn(1)+fn(2)+fn(3))*1.E11
      ENDDO
      do n=1,lev1-1
            eden3d(ilon,n)=0.
      enddo

      RETURN
      END SUBROUTINE EARTH_CHIU_MODEL
!
      subroutine interp2_ionfield(im,rlat,rlon,cormago,btoto,dipango)
!
!  VAY DANGER !!!!!!
! interp works only for given [20,91]  with 
!         j=46 center + fixed ddlat=180/90? and ddlon=360/20?
!
!
!      USE IDEA_ION_INPUT, only :
!     & cormag, btot, dipang, glat, glon, nxmag,nymag
!
      implicit none
      
      integer,intent(in)  :: im            ! number of longitude
      real,   intent(in)  :: rlat(im)      ! latitude (rad)
      real,   intent(in)  :: rlon(im)      ! longitude (rad)
!
      real,   intent(out) :: cormago(im),btoto(im),dipango(im) ! field value 
!
! local variable
      real cormag(20,91),btot(20,91),dipang(20,91),glat(91),glon(20)
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
      end subroutine interp2_ionfield

!      call tiros_ionize_data
!     &(pres(i,:), lev1,levs,ht1,emaps1,cmaps1,
!     &      djspectra1, grav(i,:),o_n(i,:),o2_n(i,:),
!     &      n2_n(i,:),adt(i,:),maglat(i),essa(i),tiros_activity_level,
!     &      GW, eden_aurora(i,:) )  
!
      subroutine tiros_ionize_data
     & (mpi_id, pres, lev1,levs,z,emaps,cmaps,djspectra,
     & grav,on,o2n, n2n,tn,gm_lat,essa1,tiros_activity_level,GW,
     & eden_aurora1D)
!    &   ,eflux,ch)     
!
!vay-2015: pass den = rho from the TOP-level program take-out comput/arrays ntot & meanmass
!          version with data/xxx/-statements.......eden_aurora1D   3-density due to aurora
!
!         output: of tiros_ionize_data
!
!SK   use idea_mpi_def,      only : mpi_id
      use idea_composition,  only :  DTR, ELCH, R_2_d, PI
!     use tirosdata
      implicit none
!
      INTEGER :: j, i, m, l, tiros_activity_level, iband
      INTEGER :: levs, lev1
      INTEGER, parameter :: jmaxwell = 6

      INTEGER, intent(in) :: mpi_id
      real, intent(in) :: pres(levs)
      real :: pres1(levs)

      real ::  Z(levs),GRAV(levs),ON(levs),
     &        o2n(levs),n2n(levs),gl,mlt,GW,gm_lat
      real :: bz, gscon, amu, e0
!
      real  :: emaps(21,20,7),cmaps(21,20,7),djspectra(15,21)
      real ::  NTOT(levs),meanmass(levs)
!
      real ::  QIONT(levs),RATIO(21),RLAM(21)
     &,den(levs),dl_lower,dl_upper,qiont_lower,qiont_upper
     &,tn(levs),mo,mo2,mn2,alpha
     &,rno,RANGE_en,pr,ratioz,rlamz,mh,mhe,q
     &,qiont_O(levs),qiont_O2(levs),qiont_N2(levs)
      real :: width(15),en(15),TE11(21),TE15(21),width_maxwell
      real :: ionchr(21),ratio_ch,en_maxwell(jmaxwell),dl(jmaxwell)
     &,qion_maxwell(levs),eden_aurora1D(levs),lognpres(8),
     &ion_recomb(8),logpres,rr
c
      real :: ch , chi , dfac , diff , dprof , ed ,
     &     eflux , essa1 , qdmsp , QT(levs) ,
     &     ri , rj , th , swbz , offset, THMagd
      INTEGER i1 , i2 , j1 , j2 , k , kk , ld , n , nn , jj , jjj
      data en/.37,.6,.92,1.37,2.01,2.91,4.19,6.,8.56,12.18,
     &17.3,24.49,36.66,54.77,81.82/
      data width/.158,.315,.315,.63,.631,1.261,1.26,2.522,
     &2.522,5.043,5.043,10.,14.81,22.13,33.06/
      data RLAM/1.49,1.52,1.51,1.48,1.43,1.37,1.30,1.22,
     11.12,1.01,0.895,0.785,0.650,0.540,0.415,0.320,0.225,
     20.14,0.08,0.04,0.0/
c
      DATA ionchr/.378 , .458 , .616 , .773 , .913 , 1.088 , 1.403 ,
     &     1.718 , 2.033 , 2.349 , 2.979 , 3.610 , 4.250 , 4.780 ,
     &     6.130 , 7.392 , 8.653 , 9.914 , 12.436 , 14.957 , 17.479/
      data lognpres/-3.425,-4.835,-5.918,-7.066,-7.784,-8.366,-9.314,
     &-10.507/
      data ion_recomb/3.20e-13,3.20e-13,2.75e-13,1.45e-13,1.13e-13,
     %8.30e-14,3.70e-14,2.00e-14/
!      DATA dprof/4.23E19 , 5.62E19 , 5.77E19 , 5.70E19 , 1.04E19 ,
!     &     1.03E20 , 1.22E20 , 1.23E20 , 0.00E19 , 8.36E19 , 2.37E20 ,
!     &     2.61E20 , 0.00E19 , 0.00E18 , 3.07E20 , 5.26E20 , 0.00E19 ,
!     &     0.00E18 , 0.00E18 , 8.57E20/
!      DATA qdmsp/15*0.0/
    
      do 16 m=1,21
   16 ratio(m) = (m-1)*0.05
      do iband=1,21
      te15(iband)=0.0
      te11(iband)=0.0
! the ionization rates will need to be normalized to TE11, which is
! the energy flux between 300eV and 20keV, which is provided by the
! TIROS energy influx maps emaps, rather than the energy from
! 300eV to 100keV, which is what the spectra were normalized to
!
! check the energy influx is normalized to 1 erg/cm2/s
      do 17 m=1,15
   17 TE15(iband)=TE15(iband)+djspectra(m,iband)*en(m)*width(m)*1.6E-06
! normalize with the energy influx 300eV to 20keV
      do 18 m=1,11
   18 TE11(iband)=TE11(iband)+djspectra(m,iband)*en(m)*width(m)*1.6E-06
!      print *, iband, TE11(iband), TE15(iband)
      enddo
      bz = 1.38e-23
      gscon = 8.314e3
      mo = 16.
      mo2 = 32.
      mn2 = 28.
      mh = 1.
      mhe = 4.
      amu = 1.661e-27
      E0=0.035
      WIDTH_maxwell=0.050
      do j = 1,jmaxwell
      en_maxwell(j) = j*0.05 - 0.025
      enddo
! initialize qiont, etc
      do i=1,levs
      qiont(i) = 0.0
      qion_maxwell(i) = 0.0
      qiont_O(i) = 0.0
      qiont_O2(i) = 0.0
      qiont_N2(i) = 0.0
      eden_aurora1D(i) = 0.0
      enddo
! convert magnetic latitude from radians to degrees
      thmagd = gm_lat * R_2_D
!     print *, 'gm_lat   thmagd  essa1', gm_lat, thmagd, essa1
      th = abs(thmagd) - 50.
      if(abs(thmagd).le.50.) goto 200
! calculate magnetic hour angle from noon in gregrees
!      essa1 = (mlt + 12.)*15.
! now passed essa1 directly
      IF ( essa1.GE.360.0 ) THEN
          essa1 = essa1 - 360.0
      ELSEIF ( essa1.LT.0.0 ) THEN
          essa1 = essa1 + 360.
      ENDIF
cc  **
      l = tiros_activity_level - 2
      IF ( l.LT.1 ) l = 1
      IF ( l.GT.7 ) l = 7
! define dfac to scale qiont later with large GW and tiros_activity_level
! Added by Zhuxiao.Li
      dfac = 1.0
      IF (tiros_activity_level.gt.9 .and. GW.gt.96.0)
     &   dfac = GW/96.0
 
cc  **
      ri = essa1/18.0 + 11.
      i1 = ri                   ! i1 =int(ri) ?
      ri = ri - i1
      IF ( i1.GT.20 ) i1 = i1 - 20
      i2 = i1 + 1
      IF ( i2.GT.20 ) i2 = i2 - 20
      rj = th/2. + 1.
      j1 = rj                   ! j1 =int(rj) ?
      rj = rj - j1
      j2 = j1 + 1
!
      eflux = rj*ri*EMAps(j2,i2,l) + (1.-rj)*ri*EMAps(j1,i2,l)
     &        + rj*(1.-ri)*EMAps(j2,i1,l) + (1.-rj)*(1.-ri)
     &        *EMAps(j1,i1,l)
      eflux = 10.**(eflux)/1000.
!      print *, 'eflux   ', eflux
!
      ch = rj*ri*CMAps(j2,i2,l) + (1.-rj)*ri*CMAps(j1,i2,l) + rj*(1.-ri)
     &     *CMAps(j2,i1,l) + (1.-rj)*(1.-ri)*CMAps(j1,i1,l)
!      print *, 'ch   ', ch
! validation tests:
! to compare with figure 4 or 5 F-R and Evans 1987
! set ch to 5 different mean energies and
! set eflux to 1.0 mW/m2
!       ch = 2.98
!       eflux=1.0
! a useful thing to compare is the ionization rate profile for ch=2.98
! with the equivalent output assuming a Maxwellian spectrum using the
! other code ionize_ipe_3 with the same mean energy of 2.98
! in this case the profiles are similar, at other other values of ch
! they can be quite different.
!
!      print *, 'set for test ch eflux', ch, eflux
!
      IF ( ch.LT.0.378 ) ch = 0.379
!      IF ( ch.GT.17.479 ) WRITE (6,99001) ch
!
      DO 300 kk = 2 , 21
         IF ( ch.LE.ionchr(kk) ) THEN
            k = kk - 1
            GOTO 400
         ENDIF
 300  CONTINUE
 400  chi = ch - ionchr(k)
      diff = ionchr(kk) - ionchr(k)
      ratio_ch = chi/diff
!      if(ratio_ch.gt.1.) print *, 'ratio_ch out of bounds', ratio_ch
c
c
99001 FORMAT ('  ch value outof bound in tiros',f10.6)
!!! loop through ipe height levels
      do 2000 i=lev1,levs
! stop at 1000km altitude
      if(z(i)*1.e-3 .gt. 1000.) goto 200                ! zwam < 1000, vay       
! set up neutral parameters
! ntot  total number density m-3
! pres  pressure Pa
! meanmass amu
! den neutral mass density kg/m3
!      grav(i)=-gr(i)/100.
!
! rewrite
!............ pres =constant......den = rho !!! all computed before
!             should be passed to "tiros", neutrals are fixed here

      ntot(i) = on(i)+o2n(i)+n2n(i)
!      pres1(i) = ntot(i)*bz*tn(i)
      meanmass(i) = (on(i)*mo+o2n(i)*mo2+n2n(i)*mn2)/ntot(i)
      den(i) = pres(i)*meanmass(i)/(gscon*tn(i))
!
! calculate ion recombination rate
!      data lognpres/-3.425,-4.835,-5.918,-7.066,-7.784,-8.366,-9.314,-10.507/
      logpres = log(pres(i))

      if (logpres.ge.-4.835)then
      rr=3.20e-13
      goto 450
      endif
      if (logpres.le.-10.507)then
      rr=2.00e-14
      goto 450
      endif
      do jjj=3,8
      if(lognpres(jjj).le.logpres)then
      jj=jjj-1
      goto 451
      endif
      enddo
  451 continue

!err      rr = ion_recomb(jj)-logpres*(ion_recomb(jjj)-ion_recomb(jj))/
!err     &(lognpres(jjj)-lognpres(jj))
!update Vay-2016/10
!NEMS.x             0000000000ADF0A4  tiros_ionize_data         449/450  idea_ion_empirmodels.f
!ion_recomb(jjj =>  ion_recomb(jj)
!
      rr = ion_recomb(jjj)-(lognpres(jjj)-logpres)*
     & (ion_recomb(jjj)-ion_recomb(jj))/(lognpres(jjj)-lognpres(jj))
  450 continue
!
!

 3000 format(1x,i5,2f10.2,1p3e9.2)
! loop through all energy bands
      DO 10 L=1,15
!      DL(L)=RNO*en(l)*EXP(-EN(L)/ALPHA)
      DL_lower = djspectra(L,K)
      DL_upper = djspectra(L,KK)
      RANGE_en=4.57E-05*EN(L)**1.75
      PR=RANGE_en*grav(i)
      RATIOZ=PRES(i)/PR
      IF(RATIOZ.GT.1.0) GOTO 20
      DO 12 M=1,21
      IF(RATIOZ.GT.RATIO(M)) GOTO 12
      RLAMZ=RLAM(M-1)+(RATIOZ-RATIO(M-1))*(RLAM(M)-RLAM(M-1))/
     1(RATIO(M)-RATIO(M-1))
      GOTO 13
   12 CONTINUE
   13 CONTINUE
      GOTO 21
   20 RLAMZ=0.0
   21 CONTINUE
      QIONt_lower=den(i)*EN(L)*RLAMZ*DL_lower*WIDTH(l)*1.E7/RANGE_en/E0
      QIONt_upper=den(i)*EN(L)*RLAMZ*DL_upper*WIDTH(l)*1.E7/RANGE_en/E0
      QIONt(i)=qiont(i)+(ratio_ch*qiont_upper+(1.-ratio_ch)*qiont_lower)
     &*eflux/(ratio_ch*te11(kk)+(1.-ratio_ch)*te11(k))
   11 CONTINUE
   10 continue
! add 0 - 300eV as Maxwellian with ch, and eflux
      alpha = ch/2.
      RNO=eflux*6.24E12/2./ALPHA**3
      DO 110 L=1,jmaxwell
      DL(L)=RNO*en_maxwell(l)*EXP(-EN_maxwell(L)/ALPHA)
      RANGE_en=4.57E-05*EN_maxwell(L)**1.75
      PR=RANGE_en*grav(i)
      RATIOZ=PRES(i)/PR
      IF(RATIOZ.GT.1.0) GOTO 120
      DO 112 M=1,21
      IF(RATIOZ.GT.RATIO(M)) GOTO 112
      RLAMZ=RLAM(M-1)+(RATIOZ-RATIO(M-1))*(RLAM(M)-RLAM(M-1))/
     1(RATIO(M)-RATIO(M-1))
      GOTO 113
  112 CONTINUE
  113 CONTINUE
      GOTO 121
  120 RLAMZ=0.0
  121 CONTINUE
      qion_maxwell(i)=qion_maxwell(i)+den(i)*en_maxwell(L)*RLAMZ*DL(L)
     &*WIDTH_maxwell/RANGE_en/E0
  110 CONTINUE
      qiont(i) = qiont(i) + qion_maxwell(i)

!  the qiont scaled by dfac with large GW and tiros_activity level
!  added by Zhuxiao
       qiont(i) = qiont(i)*dfac
!
!vay-2016, extra security >0 and non-zero
!
      if (rr.gt.0.) then 
         eden_aurora1D(i)=sqrt(qiont(i)/rr)
      else
         eden_aurora1D(i)=0.
      endif
!
      q=qiont(i)/(0.92*n2n(i)+1.5*o2n(i)+0.56*on(i))
      qiont_O(i)=(0.5*o2n(i)+0.56*on(i))*q
      qiont_O2(i)=o2n(i)*q
      qiont_N2(i)=0.92*n2n(i)*q
 2000 continue
  200 continue
!
      RETURN
      END subroutine tiros_ionize_data
!
      subroutine tiros_ionize(mpi_id, lev1,levs, pres, den,
     &  z, grav,on,o2n,n2n,tn,
     &  gm_lat,essa1,tiros_activity_level,GW, eden_aurora1D)
!
! Version with IDEA_ION_INPUT
!
!Oct 2016 VAY: a) Err-Interp logpress
!              b) RANGE (special FORT function) => RANGE_en
!              c) Clarify IF ( ch.LT.0.378 ) ch = 0.379             ! ???????? 378-379
!              d)        IF ( ch.GT.17.479 ) WRITE (6,99001) ch    ! out of bounds ???
!
!              e) P-WAM = constabt, suggestion RLAMZ(15, levs) -----constant with time
!
      USE IDEA_ION_INPUT, only : NT_21, NT_20, NT_7, N_FLX, N_BND 
      USE IDEA_ION_INPUT, only : emaps, cmaps, djspectra
!
      USE IDEA_ION_INPUT, only : jmaxwell, E0     ! 6-21-15
      USE IDEA_ION_INPUT, only : RATIO, RLAM, ion_recomb,lognpres, 
     &                           width, en, TE11,TE15, ionchr, 
     &                           width_maxwell, en_maxwell
!
!
      use idea_composition, only : PI, PI2, Pid2, DTR, R_2_D
!SK   use idea_mpi_def,     only : mpi_id
      IMPLICIT NONE
!
      INTEGER, intent(in) ::  mpi_id
      INTEGER, intent(in) ::  levs, lev1
      INTEGER, intent(in) ::  tiros_activity_level
      real ::  essa1 
      real ::  GW  
      real ::  Z(levs), grav(levs), PRES(levs), den(levs)
      real ::  TN(levs)
      real ::  ON(levs),o2n(levs),n2n(levs)
      real ::  gl, mlt,  gm_lat
!out
      real, intent(out) ::   eden_aurora1D(levs) 
!
!
!local
!
!
      real :: DL (jmaxwell)
      real ::  alpha
      real ::  QIONT(levs)
      real ::  dl_lower,dl_upper,qiont_lower,qiont_upper
     &     rno,RANGE_en,pr,ratioz,rlamz,q, ratio_ch
      real :: QMID
      real :: qiont_O(levs),qiont_O2(levs),qiont_N2(levs), QT(levs)
      
      real :: qion_maxwell(levs)
      real :: QIONt_upper
      real :: RNO

      real :: ch , chi , dfac , diff , dprof , ed ,logpres,rr,
     &     eflux, qdmsp,
     &     ri , rj , th , swbz , offset, THMagd

      INTEGER i1 , i2 , j1 , j2 , k , kk , ld , n , nn , jj , jjj
      INTEGER :: j, i, m, l, iband 
! initialize qiont, etc
      do i=1,levs
      qiont(i) = 0.0
      qion_maxwell(i) = 0.0
      qiont_O(i) = 0.0
      qiont_O2(i) = 0.0
      qiont_N2(i) = 0.0
!...................................... Zero Output
      eden_aurora1D(i) = 0.0
      enddo
    
! convert magnetic latitude from radians to degrees
      thmagd = gm_lat * R_2_D
!!!!      print *, 'gm_lat   thmagd  essa1', gm_lat, thmagd, essa1
      th = abs(thmagd) - 50.
      if(abs(thmagd).le.50.) goto 200
! calculate magnetic hour angle from noon in gregrees
!      essa1 = (mlt + 12.)*15.
! now passed essa1 directly
      IF ( essa1.GE.360.0 ) THEN
          essa1 = essa1 - 360.0
      ELSEIF ( essa1.LT.0.0 ) THEN
          essa1 = essa1 + 360.
      ENDIF
!
      l = tiros_activity_level - 2
!
! WAM-limits from 1 to 7
!
      IF ( l.LT.1 ) l = 1
      IF ( l.GT.7 ) l = 7
! define dfac to scale qiont later with large GW and tiros_activity_level
! Added by Zhuxiao.Li
      dfac = 1.0
      IF (tiros_activity_level.gt.9 .and. GW.gt.96.0)
     &   dfac = GW/96.0
!
! What's this VAY
!
      ri = essa1/18.0 + 11.
      i1 = int(ri)                   ! i1 = ri
      ri = ri - i1                   ! ri = ri -int(ri) ???? NOW
!
      IF ( i1.GT.20 ) i1 = i1 - 20
      i2 = i1 + 1
      IF ( i2.GT.20 ) i2 = i2 - 20
      rj = th/2. + 1.
      j1 = int(rj)                   ! j1 =rj
      rj = rj - j1
      j2 = j1 + 1
! both rj <1 & ri < 1
      eflux = rj*ri*EMAps(j2,i2,l) + (1.-rj)*ri*EMAps(j1,i2,l)
     &        + rj*(1.-ri)*EMAps(j2,i1,l) + (1.-rj)*(1.-ri)
     &        *EMAps(j1,i1,l)
!
      eflux = 10.**(eflux) * 1.e-3
      if (mpi_id == 0)  print *, 'ion-eflux   ', eflux
!
!
      ch = rj*ri*CMAps(j2,i2,l) + (1.-rj)*ri*CMAps(j1,i2,l) + rj*(1.-ri)
     &     *CMAps(j2,i1,l) + (1.-rj)*(1.-ri)*CMAps(j1,i1,l)
!
      if (mpi_id == 0) print *, 'ion-ch   ', ch
! validation tests:
! to compare with figure 4 or 5 F-R and Evans 1987
! set ch to 5 different mean energies and
! set eflux to 1.0 mW/m2
!       ch = 2.98
!       eflux=1.0
! a useful thing to compare is the ionization rate profile for ch=2.98
! with the equivalent output assuming a Maxwellian spectrum using the
! other code ionize_ipe_3 with the same mean energy of 2.98
! in this case the profiles are similar, at other other values of ch
! they can be quite different.
!       if (mpi_id == 0) then
!       print *, 'set for test ch eflux', ch, eflux
!       endif
      IF ( ch.LT.0.378 ) ch = 0.379             ! ???????? 378-379
      IF ( ch.GT.17.479 ) WRITE (6,99001) ch
!
      DO 300 kk = 2 , 21
         IF ( ch.LE.ionchr(kk) ) THEN
            k = kk - 1
            GOTO 400
         ENDIF
 300  CONTINUE
 400  chi = ch - ionchr(k)
      diff = ionchr(kk) - ionchr(k)
      ratio_ch = chi/diff
!      if(ratio_ch.gt.1.) print *, 'ratio_ch out of bounds', ratio_ch
!
99001 FORMAT ('  ch value outof bound in tiros',f10.6)

!!! loop through ipe height levels OR WAM-levels

      do 2000 i=lev1,levs
! 
! set up neutral parameters
! ntot  total number density m-3
! pres  pressure Pa
! meanmass amu
! den neutral mass density kg/m3
!      grav(i)=-gr(i)/100.
!
!!!!VAY OCT-2015:  WHAT's this  ??????????? ALL DONE before in WAM-physics
!
!      ntot = (on(i)+o2n(i)+n2n(i))
!      meanmass(i) = (on(i)*mo+o2n(i)*mo2+n2n(i)*mn2)/ntot
!      den(i) = pres(i)*meanmass(i)/(gscon*tn(i))
!      pres(i) = ntot*bz*tn(i)        ! this is fixed in pressure-based WAM system !!!!!!
!
! calculate ion recombination rate ....CONSTANT with  TIME ????? pres = constant in WAM !!!!!
!
!                                      SUGGESTION.. USE CASE STATEMENTS or IF-else
! stop at 1000km altitude
      if(z(i)*1.e-3 .gt. 1000.) goto 200   ! do we need this Z-WAM < 650 km
!
! new Z-interpolation   for rr > 0 (apparently)
!
      logpres = log(pres(i))
      if (logpres.ge.-4.835)then
      rr=3.20e-13
      goto 450
      endif
      if (logpres.le.-10.507)then
      rr=2.00e-14
      goto 450
      endif
      do jjj=3,8
         if(lognpres(jjj).le.logpres)then
         jj=jjj-1
         goto 451
         endif
      enddo
  451 continue
!
! if it does not change with time introduce RR(levs) in the ION_INIT
!      all goto's ...... take cpu-time   
!      looks like interpol in pressure
!Fix VAY OCT-2016
!
      rr = ion_recomb(jjj)-(lognpres(jjj)-logpres)*
     & (ion_recomb(jjj)-ion_recomb(jj))/(lognpres(jjj)-lognpres(jj))
!
!err      rr = ion_recomb(jj)-logpres*(ion_recomb(jjj)-ion_recomb(jj))/
!err     &(lognpres(jjj)-lognpres(jj))
!
  450 continue
!      print *, i, z(i), pres(i), ntot(i), meanmass(i), den(i)
!      write(6,3000) i, z(i), meanmass(i), pres(i), ntot(i), den(i)
 3000 format(1x,i5,2f10.2,1p3e9.2)
!
! loop through all energy bands
      DO 10 L=1,15
!      DL(L)=RNO*en(l)*EXP(-EN(L)/ALPHA)
      DL_lower = djspectra(L,K)
      DL_upper = djspectra(L,KK)

      RANGE_en=4.57E-05*EN(L)**1.75
      PR=RANGE_en*grav(i)
      RATIOZ=PRES(i)/PR
      IF(RATIOZ.GT.1.0) GOTO 20
      DO 12 M=1,21
      IF(RATIOZ.GT.RATIO(M)) GOTO 12
      RLAMZ=RLAM(M-1)+(RATIOZ-RATIO(M-1))*(RLAM(M)-RLAM(M-1))/
     &      (RATIO(M)-RATIO(M-1))
      GOTO 13
   12 CONTINUE
   13 CONTINUE
      GOTO 21
   20 RLAMZ=0.0
   21 CONTINUE
! ...................... again interpolation in RATIO,,,RATIOZ
! similar suggestion RLAMZ(15, levs) -----constant with time
!
      QMID = den(i)*EN(L)*RLAMZ*WIDTH(l)*1.E7/RANGE_en/E0
      QIONt_lower= Qmid*DL_lower
      QIONt_upper= Qmid*DL_upper
      QIONt(i)=qiont(i)+(ratio_ch*qiont_upper+(1.-ratio_ch)*qiont_lower)
     &*eflux/(ratio_ch*te11(kk)+(1.-ratio_ch)*te11(k))
   11 CONTINUE

   10 continue     ! loop through all energy bands

! add 0 - 300eV as Maxwellian with ch, and eflux

      if (ch.le.0.) then
         print *, ' TIROS: CH=0 ', CH
         ch =2.98
      endif
      alpha = .5*ch

      RNO=eflux*6.24E12/2./ALPHA**3
!
      DO 110 L=1,jmaxwell
      DL(L)=RNO*en_maxwell(l)*EXP(-EN_maxwell(L)/ALPHA)
      RANGE_en=4.57E-05*EN_maxwell(L)**1.75
      PR=RANGE_en*grav(i)
      RATIOZ=PRES(i)/PR
      IF(RATIOZ.GT.1.0) GOTO 120
      DO 112 M=1,21
      IF(RATIOZ.GT.RATIO(M)) GOTO 112
      RLAMZ=RLAM(M-1)+(RATIOZ-RATIO(M-1))*(RLAM(M)-RLAM(M-1))/
     &      (RATIO(M)-RATIO(M-1))
      GOTO 113
  112 CONTINUE
  113 CONTINUE
      GOTO 121
  120 RLAMZ=0.0
  121 CONTINUE
!......................similar RLAMZ_JMX(levs, 6)
!                       
      qion_maxwell(i)=qion_maxwell(i)
     &  +den(i)*en_maxwell(L)*RLAMZ*DL(L)*WIDTH_maxwell/RANGE_en/E0
!
  110 CONTINUE   !  loop for300eV as Maxwellian with ch, and eflux
!
      qiont(i) = qiont(i) + qion_maxwell(i)
!  the qiont scaled by dfac with large GW and tiros_activity level
!  added by Zhuxiao
       qiont(i) = qiont(i)*dfac

!vay-2016
      if (rr.gt.0.) then 
         eden_aurora1D(i)=sqrt(qiont(i)/rr)
      else
         eden_aurora1D(i)=0.
      endif
!
! proxies for charge from NEUTRALS ????
!
      q=qiont(i)/(0.92*n2n(i)+1.5*o2n(i)+0.56*on(i))
      qiont_O(i)=(0.5*o2n(i)+0.56*on(i))*q
      qiont_O2(i)=o2n(i)*q
      qiont_N2(i)=0.92*n2n(i)*q
 2000 continue
  200 continue 
!
      RETURN
      END subroutine tiros_ionize
