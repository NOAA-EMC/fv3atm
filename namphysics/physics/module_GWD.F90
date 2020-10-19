!-----------------------------------------------------------------------
!
      MODULE MODULE_GWD
!
!-----------------------------------------------------------------------
!
!***  Module for Gravity Wave Drag (GWD) and Mountain Blocking (MB)
!
!***  Initially incorporated into the WRF NMM from the GFS by B. Ferrier
!***  in April/May 2007.  
!***  Ratko added in NMM-B (July '08)
!***  Mountain blocking made more robust in July 2016 by BSF
!
!***  Search for "ORIGINAL DOCUMENTATION BLOCK" for further description.
!
!-----------------------------------------------------------------------
!
      INTEGER, PARAMETER :: KIND_PHYS=SELECTED_REAL_KIND(13,60) ! the '60' maps to 64-bit real
      INTEGER,PRIVATE,SAVE :: NMTVR, IDBG, JDBG
      REAL (KIND=KIND_PHYS),PRIVATE :: DELTIM,RDELTIM
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
!-- Initialize variables used in GWD + MB
!
      SUBROUTINE GWD_init (DTPHS,RESTRT                                 &
                           ,CLEFFAMP,DPHD                               &
                           ,CLEFF                                       &
                           ,CEN_LAT,CEN_LON                             &
                           ,GLAT,GLON                                   &
                           ,CROT,SROT,HANGL                             &
                           ,IMS,IME,JMS,JME,LM)

!
      IMPLICIT NONE
!
!== INPUT:
!-- CEN_LAT, CEN_LON - central latitude, longitude (degrees)
!-- DPHD - DY in degrees latitude
!-- CLEFFAMP - amplification factor for variable CLEFF
!-- RESTRT - logical flag for restart file (true) or WRF input file (false)
!-- GLAT, GLON - central latitude, longitude at mass points (radians)
!-- CROT, SROT - cosine and sine of the angle between Earth and model coordinates
!-- HANGL  - angle of the mountain range w/r/t east (convert to degrees)
!
!-- Saved variables within module:
!-- NMTVR - number of input 2D orographic fields
!-- GRAV = gravitational acceleration
!-- DELTIM - physics time step (s)
!-- RDELTIM - reciprocal of physics time step (s)
!
!
      REAL, INTENT(IN) :: DTPHS,CEN_LAT,CEN_LON,DPHD,CLEFFAMP
      REAL, INTENT(OUT) :: CLEFF
      LOGICAL, INTENT(IN) :: RESTRT
      REAL, INTENT(IN), DIMENSION (ims:ime,jms:jme) :: GLON,GLAT
      REAL, INTENT(OUT), DIMENSION (ims:ime,jms:jme) :: CROT,SROT
      REAL, INTENT(INOUT), DIMENSION (ims:ime,jms:jme) :: HANGL
      INTEGER, INTENT(IN) :: IMS,IME,JMS,JME,LM
!
!-- Local variables:
!
      REAL, PARAMETER :: DPHD0=0.108      !-- Nominal 12-km grid resolution
      REAL :: DTR,RTD,TPH0,TLM0,CTPH0,STPH0,RELM,SRLM,CRLM,SPH,CPH,TPH,RCTPH
      INTEGER :: I,J
!
!---------------------------------------------------------------------
!
      NMTVR=14            !-- 14 input fields for orography
      DELTIM=DTPHS
      RDELTIM=1./DTPHS
!
!-- Calculate angle of rotation between Earth and model coordinates,
!   but pass back out cosine (CROT) and sine (SROT) of this angle
!
      DTR=ACOS(-1.)/180. !-- convert from degrees to radians
      TPH0=CEN_LAT*DTR
      TLM0=CEN_LON*DTR
      CTPH0=COS(TPH0)
      STPH0=SIN(TPH0)
      DO J=JMS,JME
        DO I=IMS,IME
          RELM=GLON(I,J)-TLM0
          SRLM=SIN(RELM)
          CRLM=COS(RELM)
          SPH=SIN(GLAT(I,J))
          CPH=COS(GLAT(I,J))
          TPH=ASIN(CTPH0*SPH-STPH0*CPH*CRLM)
          RCTPH=1.0/COS(TPH)
          SROT(I,J)=STPH0*SRLM*RCTPH
          CROT(I,J)=(CTPH0*CPH+STPH0*SPH*CRLM)*RCTPH
        ENDDO    !-- I
      ENDDO      !-- J
!-- Convert from radians to degrees
      RTD=1./DTR         !-- convert from radians to degrees
      DO J=JMS,JME
        DO I=IMS,IME
          HANGL(I,J)=RTD*HANGL(I,J)  !-- convert to degrees (+/-90 deg)
        ENDDO    !-- I
      ENDDO      !-- J
!
!-- Scale cleff to be w/r/t a nominal value of 1e-5 for 12-km grid Launcher runs
!     * Launcher 12-km runs used cleff=0.5e-5*SQRT(IMX/192), where IMX=IDE-1=705
!       for the 12-km air quality (na12aq) domain
!
      CLEFF=CLEFFAMP*1.E-5*SQRT(DPHD/DPHD0)
!     write(0,*) 'dphd,cleff=',dphd,cleff
!
      END SUBROUTINE GWD_init
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE GWD_driver(DTPHS,U,V,T,Q,Z,DP,PINT,PMID,EXNR,KPBL      &
                           ,HSTDV,HCNVX,HASYW,HASYS,HASYSW,HASYNW       &
                           ,HLENW,HLENS,HLENSW,HLENNW                   &
                           ,HANGL,HANIS,HSLOP,HZMAX,CROT,SROT           &
                           ,CDMB,CLEFF,SIGFAC,FACTOP,RLOLEV,DPMIN       &
                           ,DUDT,DVDT                                   &
                           ,IMS,IME,JMS,JME,LM)
!
!== INPUT:
!-- U, V - zonal (U), meridional (V) winds at mass points (m/s)
!-- T, Q - temperature (C), specific humidity (kg/kg)
!-- DP - pressure thickness (Pa)
!-- Z - geopotential height (m)
!-- PINT, PMID - interface and midlayer pressures, respectively (Pa)
!-- EXNR - (p/p0)**(Rd/Cp)
!-- KPBL - vertical index at PBL top
!-- HSTDV - orographic standard deviation
!-- HCNVX - normalized 4th moment of the orographic convexity
!-- Template for the next two sets of 4 arrays:
!             NWD  1   2   3   4   5   6   7   8
!              WD  W   S  SW  NW   E   N  NE  SE
!-- Orographic asymmetry (HASYx, x=1-4) for upstream & downstream flow (4 planes)
!-- * HASYW - orographic asymmetry for upstream & downstream flow in W-E plane
!-- * HASYS - orographic asymmetry for upstream & downstream flow in S-N plane
!-- * HASYSW - orographic asymmetry for upstream & downstream flow in SW-NE plane
!-- * HASYNW - orographic asymmetry for upstream & downstream flow in NW-SE plane
!-- Orographic length scale or mountain width (4 planes)
!-- * HLENW - orographic length scale for upstream & downstream flow in W-E plane
!-- * HLENS - orographic length scale for upstream & downstream flow in S-N plane
!-- * HLENSW - orographic length scale for upstream & downstream flow in SW-NE plane
!-- * HLENNW - orographic length scale for upstream & downstream flow in NW-SE plane
!-- HANGL  - angle (degrees) of the mountain range w/r/t east
!-- HANIS - anisotropy/aspect ratio of orography
!-- HSLOP - slope of orography
!-- HZMAX - max height above mean orography
!-- CROT, SROT - cosine & sine of the angle between Earth & model coordinates
!
!== OUTPUT:
!-- DUDT, DVDT - zonal, meridional wind tendencies
!-- UGWDsfc, VGWDsfc - zonal, meridional surface wind stresses (N/m**2)
!
!== INPUT indices:
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!
!-- INPUT variables:
!
      INTEGER, INTENT(IN) :: ims,ime,jms,jme,lm

      REAL, INTENT(IN):: DTPHS,CDMB,CLEFF,SIGFAC,FACTOP,RLOLEV,DPMIN
      REAL, INTENT(IN), DIMENSION (ims:ime, jms:jme, 1:lm+1) ::        &
     &                                   Z, PINT
      REAL, INTENT(IN), DIMENSION (ims:ime, jms:jme, 1:lm) ::        &
     &                                   U,V,T,Q,DP,PMID,EXNR
      REAL, INTENT(IN), DIMENSION (ims:ime, jms:jme) :: HSTDV,HCNVX     &
     &      ,HASYW,HASYS,HASYSW,HASYNW,HLENW,HLENS,HLENSW,HLENNW,HANGL  &
     &      ,HANIS,HSLOP,HZMAX,CROT,SROT
      INTEGER, INTENT(IN), DIMENSION (ims:ime, jms:jme) :: KPBL
!
!-- OUTPUT variables:
!
      REAL, INTENT(OUT), DIMENSION (ims:ime, jms:jme, 1:lm) ::       &
     &                                                        DUDT,DVDT
!--- when NPS is done with GWD, add wind stresses in output
      REAL,              DIMENSION (ims:ime, jms:jme) :: UGWDsfc,VGWDsfc
!
!-- Local variables
!-- DUsfc, DVsfc - zonal, meridional wind stresses (diagnostics)
!
      INTEGER, PARAMETER :: IM=1    !-- Reduces changes in subroutine GWPDS
      REAL(KIND=KIND_PHYS), PARAMETER :: G=9.806, GHALF=.5*G            &
     &,                                  THRESH=1.E-6, dtlarge=1.
      INTEGER, DIMENSION (IM) :: LPBL
      REAL(KIND=KIND_PHYS), DIMENSION (IM,4) :: OA4,CLX4
      REAL(KIND=KIND_PHYS), DIMENSION (IM) :: DUsfc,DVsfc               &
     &,                              HPRIME,OC,THETA,SIGMA,GAMMA,ELVMAX
      REAL(KIND=KIND_PHYS), DIMENSION (IM,1:LM) :: DUDTcol,DVDTcol   &
     &,                    Ucol,Vcol,Tcol,Qcol,DPcol,Pcol,EXNcol,PHIcol
      REAL(KIND=KIND_PHYS), DIMENSION (IM,1:LM+1) :: PINTcol,PHILIcol
      INTEGER :: I,J,IJ,K,KFLIP,IMX,Imid,Jmid
      REAL :: Ugeo,Vgeo,Umod,Vmod, TERRtest,TERRmin
      REAL(KIND=KIND_PHYS) :: TEST
!
!--------------------------  Executable below  -------------------------
!
!-- Initialize variables
!
      DELTIM=DTPHS
      RDELTIM=1./DTPHS
!
      DO K=1,LM
      DO J=JMS,JME
      DO I=IMS,IME
        DUDT(I,J,K)=0.
        DVDT(I,J,K)=0.
      ENDDO
      ENDDO
      ENDDO
!
      DO J=JMS,JME
      DO I=IMS,IME
        UGWDsfc(I,J)=0.
        VGWDsfc(I,J)=0.
      ENDDO
      ENDDO
!
      IMX=IDE-1
!
      DO J=JMS,JME
        DO I=IMS,IME
          if (kpbl(i,j)<=1 .or. kpbl(i,j)>LM) CYCLE  !-- Go to next grid column
!
!-- Initial test to see if GWD calculations should be made, otherwise skip
!
          TERRtest=HZMAX(I,J)+SIGFAC*HSTDV(I,J)
          TERRmin=Z(I,J,LM)-Z(I,J,LM+1)
          IF (TERRtest < TERRmin) CYCLE             !-- Go to next grid column
!
!-- For debugging:
!
          DO K=1,LM
            KFLIP = LM-K+1
            DUDTcol(IM,K)=0.
            DVDTcol(IM,K)=0.
!
!-- Transform/rotate winds from model to geodetic (Earth) coordinates
!
            Ucol(IM,K)=U(I,J,KFLIP)*CROT(I,J)+V(I,J,KFLIP)*SROT(I,J)
            Vcol(IM,K)=V(I,J,KFLIP)*CROT(I,J)-U(I,J,KFLIP)*SROT(I,J)
!
            Tcol(IM,K)=T(I,J,KFLIP)
            Qcol(IM,K)=Q(I,J,KFLIP)
!
!-- Convert from Pa to centibars, which is what's used in subroutine GWD_col
!
            DPcol(IM,K)=.001*DP(I,J,KFLIP)
            PINTcol(IM,K)=.001*PINT(I,J,KFLIP+1)
            Pcol(IM,K)=.001*PMID(I,J,KFLIP)
            EXNcol(IM,K)=EXNR(I,J,KFLIP)
!
!-- Next 2 fields are geopotential above the surface at the lower interface 
!   and at midlayer
!
            PHILIcol(IM,K)=G*(Z(I,J,KFLIP+1)-Z(I,J,LM+1))
            PHIcol(IM,K)=GHALF*(Z(I,J,KFLIP+1)+Z(I,J,KFLIP))-G*Z(I,J,LM+1)
          ENDDO   !- K
!
          PINTcol(IM,LM+1)=.001*PINT(I,J,1)
          PHILIcol(IM,LM+1)=G*(Z(I,J,1)-Z(I,J,LM+1))
!
!-- Terrain-specific inputs:
!
          HPRIME(IM)=HSTDV(I,J)   !-- standard deviation of orography
          OC(IM)=HCNVX(I,J)       !-- Normalized convexity
          OA4(IM,1)=HASYW(I,J)    !-- orographic asymmetry in W-E plane
          OA4(IM,2)=HASYS(I,J)    !-- orographic asymmetry in S-N plane
          OA4(IM,3)=HASYSW(I,J)   !-- orographic asymmetry in SW-NE plane
          OA4(IM,4)=HASYNW(I,J)   !-- orographic asymmetry in NW-SE plane
          CLX4(IM,1)=HLENW(I,J)   !-- orographic length scale in W-E plane
          CLX4(IM,2)=HLENS(I,J)   !-- orographic length scale in S-N plane
          CLX4(IM,3)=HLENSW(I,J)  !-- orographic length scale in SW-NE plane
          CLX4(IM,4)=HLENNW(I,J)  !-- orographic length scale in NW-SE plane
          THETA(IM)=HANGL(I,J)       !
          SIGMA(IM)=HSLOP(I,J)       !
          GAMMA(IM)=HANIS(I,J)       !
          ELVMAX(IM)=HZMAX(I,J)      !
          LPBL(IM)=LM+1-KPBL(I,J)      !
!
!-- Output (diagnostics)
!
          DUsfc(IM)=0.             !-- U wind stress
          DVsfc(IM)=0.             !-- V wind stress
!
!=======================================================================
!
          CALL GWD_col(DVDTcol,DUDTcol, DUsfc,DVsfc                     & ! Output
     &,              Ucol,Vcol,Tcol,Qcol,PINTcol,DPcol,Pcol,EXNcol      & ! Met input
     &,              PHILIcol,PHIcol                                    & ! Met input
     &,              HPRIME,OC,OA4,CLX4,THETA,SIGMA,GAMMA,ELVMAX        & ! Topo input
     &,              CDMB,CLEFF,SIGFAC,FACTOP,RLOLEV,DPMIN              & ! Tunable coefficients
     &,              LPBL,IMX,IM,LM)                                      ! Indices
!
!=======================================================================
!
          DO K=1,LM
            KFLIP=LM-K+1
            TEST=ABS(DUDTcol(IM,K))+ABS(DVDTcol(IM,K))
            IF (TEST > THRESH) THEN
!
!-- First update winds in geodetic coordinates
!
              Ugeo=Ucol(IM,K)+DUDTcol(IM,K)*DELTIM
              Vgeo=Vcol(IM,K)+DVDTcol(IM,K)*DELTIM
!
!-- Transform/rotate winds from geodetic back to model coordinates
!
              Umod=Ugeo*CROT(I,J)-Vgeo*SROT(I,J)
              Vmod=Ugeo*SROT(I,J)+Vgeo*CROT(I,J)
!
!-- Calculate wind tendencies from the updated model winds
!
              DUDT(I,J,KFLIP)=RDELTIM*(Umod-U(I,J,KFLIP))
              DVDT(I,J,KFLIP)=RDELTIM*(Vmod-V(I,J,KFLIP))
!
!dtest=abs(dudt(i,k,j))+abs(dvdt(i,k,j))
            ENDIF     !- IF (TEST > THRESH) THEN
!
          ENDDO   !- K
!
!-- Transform/rotate surface wind stresses from geodetic to model coordinates
!
          UGWDsfc(I,J)=DUsfc(IM)*CROT(I,J)-DVsfc(IM)*SROT(I,J)
          VGWDsfc(I,J)=DUsfc(IM)*SROT(I,J)+DVsfc(IM)*CROT(I,J)
!
        ENDDO     !- I
      ENDDO       !- J
!
      END SUBROUTINE GWD_driver
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE GWD_col (A,B, DUsfc,DVsfc                              &  !-- Output
     &, U1,V1,T1,Q1, PRSI,DEL,PRSL,PRSLK, PHII,PHIL                     &  !-- Met inputs
     &, HPRIME,OC,OA4,CLX4,THETA,SIGMA,GAMMA,ELVMAX                     &  !-- Topo inputs
     &, CDMB,CLEFF,SIGFAC,FACTOP,RLOLEV,DPMIN                           &  !-- Tunable coefficients
     &, KPBL,IMX,IM,LM)                                                    !-- Input
!
!-- "A", "B" (from GFS) in GWD_col are DVDTcol, DUDTcol, respectively in GWD_driver
!
!=== Output fields
!
!-- A (DVDT), B (DUDT) - output zonal & meridional wind tendencies in Earth coordinates (m s^-2)
!-- DUsfc, DVsfc - surface zonal meridional wind stresses in Earth coordinates (m s^-1?)
!
!=== Input fields
!
!-- U1, V1 - zonal, meridional wind (m/s)
!-- T1 - temperature (deg K)
!-- Q1 - specific humidity (kg/kg)
!-- PRSI - lower interface pressure in centibars (1000 Pa)
!-- DEL - pressure thickness of layer in centibars (1000 Pa)
!-- PRSL - midlayer pressure in centibars (1000 Pa)
!-- PRSLK - Exner function, (P/P0)**(Rd/CP)
!-- PHII - lower interface geopotential in mks units
!-- PHIL - midlayer geopotential in mks units
!-- KDT - number of time steps into integration for diagnostics
!-- HPRIME - orographic standard deviation
!-- OC - normalized 4th moment of the orographic convexity
!-- OA4 - orographic asymmetry for upstream & downstream flow measured 
!         along 4 vertical planes (W-E, S-N, SW-NE, NW-SE)
!-- CLX4 - orographic length scale or mountain width measured along
!          4 vertical planes (W-E, S-N, SW-NE, NW-SE)
!-- THETA - angle of the mountain range w/r/t east
!-- SIGMA - slope of orography
!-- GAMMA - anisotropy/aspect ratio
!-- ELVMAX - max height above mean orography
!-- KPBL(IM) - vertical index at the top of the PBL
!-- IMX - points in a grid row
!-- KM - number of vertical levels
!
!#######################################################################
!##################  ORIGINAL DOCUMENTATION BLOCK  #####################
!######  The following comments are from the original GFS code  ########
!#######################################################################
!   ********************************************************************
! ----->  I M P L E M E N T A T I O N    V E R S I O N   <----------
!
!          --- Not in this code --  History of GWDP at NCEP----
!              ----------------     -----------------------
!  VERSION 3  MODIFIED FOR GRAVITY WAVES, LOCATION: .FR30(V3GWD)  *J*
!---       3.1 INCLUDES VARIABLE SATURATION FLUX PROFILE CF ISIGST
!---       3.G INCLUDES PS COMBINED W/ PH (GLAS AND GFDL)
!-----         ALSO INCLUDED IS RI  SMOOTH OVER A THICK LOWER LAYER
!-----         ALSO INCLUDED IS DECREASE IN DE-ACC AT TOP BY 1/2
!-----     THE NMC GWD INCORPORATING BOTH GLAS(P&S) AND GFDL(MIGWD)
!-----        MOUNTAIN INDUCED GRAVITY WAVE DRAG 
!-----    CODE FROM .FR30(V3MONNX) FOR MONIN3
!-----        THIS VERSION (06 MAR 1987)
!-----        THIS VERSION (26 APR 1987)    3.G
!-----        THIS VERSION (01 MAY 1987)    3.9
!-----    CHANGE TO FORTRAN 77 (FEB 1989)     --- HANN-MING HENRY JUANG
!----- 
!
!   VERSION 4
!                ----- This code -----
!
!-----   MODIFIED TO IMPLEMENT THE ENHANCED LOW TROPOSPHERIC GRAVITY
!-----   WAVE DRAG DEVELOPED BY KIM AND ARAKAWA(JAS, 1995).
!        Orographic Std Dev (hprime), Convexity (OC), Asymmetry (OA4)
!        and Lx (CLX4) are input topographic statistics needed.
!
!-----   PROGRAMMED AND DEBUGGED BY HONG, ALPERT AND KIM --- JAN 1996.
!-----   debugged again - moorthi and iredell --- may 1998.
!-----
!       Further Cleanup, optimization and modification
!                                       - S. Moorthi May 98, March 99.
!-----   modified for usgs orography data (ncep office note 424)
!        and with several bugs fixed  - moorthi and hong --- july 1999.
!
!-----   Modified & implemented into NRL NOGAPS
!                                       - Young-Joon Kim, July 2000
!-----
!   VERSION lm MB  (6): oz fix 8/2003
!                ----- This code -----
!
!------   Changed to include the Lott and Miller Mtn Blocking
!         with some modifications by (*j*)  4/02
!        From a Principal Coordinate calculation using the
!        Hi Res 8 minute orography, the Angle of the
!        mtn with that to the East (x) axis is THETA, the slope
!        parameter SIGMA. The anisotropy is in GAMMA - all  are input
!        topographic statistics needed.  These are calculated off-line
!        as a function of model resolution in the fortran code ml01rg2.f,
!        with script mlb2.sh.   (*j*)
!-----   gwdps_mb.f version (following lmi) elvmax < hncrit (*j*)
!        MB3a expt to enhance elvmax mtn hgt see sigfac & hncrit
!-----
!----------------------------------------------------------------------C
!
      IMPLICIT NONE
!
!-- INPUT:
!
      INTEGER, INTENT(IN) :: IM,IMX,LM
      REAL,    INTENT(IN) :: CDMB,CLEFF,SIGFAC,FACTOP,RLOLEV,DPMIN

      REAL(kind=kind_phys), INTENT(IN), DIMENSION(IM,1:LM) ::        &
     &                                 U1,V1,T1,Q1,DEL,PRSL,PRSLK,PHIL
      REAL(kind=kind_phys), INTENT(IN), DIMENSION(IM,1:LM+1) ::      &
     &                                                       PRSI,PHII
      REAL(kind=kind_phys), INTENT(IN), DIMENSION(IM,4) :: OA4,CLX4
      REAL(kind=kind_phys), INTENT(IN), DIMENSION(IM) ::                &
     &                              HPRIME,OC,THETA,SIGMA,GAMMA,ELVMAX
      INTEGER, INTENT(IN), DIMENSION(IM) :: KPBL
!
!-- OUTPUT:
!
      REAL(kind=kind_phys), INTENT(INOUT), DIMENSION(IM,1:LM) :: A,B
      REAL(kind=kind_phys), INTENT(INOUT), DIMENSION(IM) :: DUsfc,DVsfc
!
!-----------------------------------------------------------------------
!-- LOCAL variables:
!-----------------------------------------------------------------------
!
!     Some constants
!
!
      REAL(kind=kind_phys), PARAMETER :: PI=3.1415926535897931        &
     &,        G=9.806, CP=1004.6, RD=287.04, RV=461.6                &
     &,        FV=RV/RD-1., RDI=1./RD, GOR=G/RD, GR2=G*GOR, GOCP=G/CP &
     &,        ROG=1./G, ROG2=ROG*ROG, ROGN=-1000.*ROG                &
     &,        DW2MIN=1., RIMIN=-100., RIC=0.25, BNV2MIN=1.0E-5       &
     &,        EFMIN=0.0, EFMAX=10.0, hpmax=2400.0                    & !-- hpmax was 200.0
     &,        FRC=1.0, CE=0.8, CEOFRC=CE/FRC, frmax=100.             &
     &,        CG=0.5, GMAX=1.0, CRITAC=5.0E-4, VELEPS=1.0            &
     &,        HZERO=0., HONE=1., HE_2=.01, HE_1=0.1                  & 
     &,        PHY180=180.,RAD_TO_DEG=PHY180/PI,DEG_TO_RAD=PI/PHY180  & 
!
!-- Lott & Miller mountain blocking => aka "lm mtn blocking"
!
!  hncrit set to 8000m and sigfac added to enhance elvmax mtn hgt
     &,  hncrit=8000.      &    ! Max value in meters for ELVMAX (*j*)
     &,  hminmt=50.        &    ! min mtn height (*j*)
     &,  hpmin=1.          &    ! minimum hprime (std dev of terrain)
     &,  hstdmin=25.       &    ! min orographic std dev in height
     &,  minwnd=0.1             ! min wind component (*j*)
!rv  &,  dpmin=5.0              ! Minimum thickness of the reference layer (centibars)
!rv                             ! values of dpmin=0, 20 have also been used
!
      integer, parameter :: mdir=8
      real(kind=kind_phys), parameter :: FDIR=mdir/(PI+PI)
!
!-- Template:
!             NWD  1   2   3   4   5   6   7   8
!              WD  W   S  SW  NW   E   N  NE  SE
!
      integer,save :: nwdir(mdir)
      data nwdir /6,7,5,8,2,3,1,4/
!
      LOGICAL ICRILV(IM)
!
!----   MOUNTAIN INDUCED GRAVITY WAVE DRAG
!
!
! for lm mtn blocking
      real(kind=kind_phys), DIMENSION(IM) :: WK,PE,EK,ZBK,UP,TAUB,XN    &
     & ,YN,UBAR,VBAR,ULOW,OA,CLX,ROLL,ULOI,DTFAC,XLINV,DELKS,DELKS1     &
     & ,SCOR,BNV2bar, ELEVMX   ! ,PSTAR
!
      real(kind=kind_phys), DIMENSION(IM,1:LM) ::                    &
     &                      BNV2LM,DB,ANG,UDS,BNV2,RI_N,TAUD,RO,VTK,VTJ
      real(kind=kind_phys), DIMENSION(IM,1:LM-1) :: VELCO
      real(kind=kind_phys), DIMENSION(IM,1:LM+1) :: TAUP
      real(kind=kind_phys), DIMENSION(LM-1) :: VELKO
!
      integer, DIMENSION(IM) ::                                         &
     &                 kref,kint,iwk,iwk2,ipt,kreflm,iwklm,iptlm,idxzb
!
! for lm mtn blocking
!
      real(kind=kind_phys) :: ZLEN, DBTMP, R, PHIANG,       DBIM        &
     &,                   xl,     rcsks, bnv,   fr                      &
     &,                   brvf,   tem,   tem1,  tem2, temc, temv        &
     &,                   wdir,   ti,    rdz,   dw2,   shr2, bvf2       &
     &,                   rdelks, wtkbj, efact, coefm, gfobnv           &
     &,                   scork,  rscor, hd,    fro,   rim,  sira       &
     &,                   dtaux,  dtauy, pkp1log, pklog, cosang, sinang
!
      integer :: ncnt, kmm1, kmm2, lcap, lcapp1, kbps, kbpsp1,kbpsm1    &
     &, kmps, kmpsp1, idir, nwd, i, j, k, klcap, kp1, kmpbl, npt, npr   &
     &, idxm1, ktrial, klevm1, kmll,kmds, KM                            &
!     &, ihit,jhit                                                       &
     &, ME              !-- processor element for debugging

!
!-----------------------------------------------------------------------
!
      KM = LM
      npr = 0
      DO I = 1, IM
         DUsfc(I) = 0.
         DVsfc(I) = 0.
!
!-- ELEVMX is a local array that could be changed below
!
         ELEVMX(I) = ELVMAX(I)
      ENDDO
!
!-- Note that A, B already set to zero as DUDTcol, DVDTcol in subroutine GWD_driver
!
      ipt = 0
      npt = 0
      IF (NMTVR >= 14) then 
        DO I = 1,IM
          IF (elvmax(i) > HMINMT .AND. hprime(i) > hpmin) then
             npt = npt + 1
             ipt(npt) = i
          ENDIF
        ENDDO
      ELSE
        DO I = 1,IM
          IF (hprime(i) > hpmin) then
            npt = npt + 1
            ipt(npt) = i
          ENDIF
        ENDDO
      ENDIF    !-- IF (NMTVR >= 14) then 
!

!
!-- Note important criterion for immediately exiting routine!
!
      IF (npt <= 0) RETURN     ! No gwd/mb calculation done!
!
      do i=1,npt
        IDXZB(i) = 0
      enddo
!
      DO K = 1, KM
      DO I = 1, IM
      DB(I,K) = 0.
      ANG(I,K) = 0.
      UDS(I,K) = 0.
      ENDDO
      ENDDO
!
      KMM1   = KM - 1
      KMM2   = KM - 2
      LCAP   = KM
      LCAPP1 = LCAP + 1
!
!
      IF (NMTVR .eq. 14) then 
! ----  for lm and gwd calculation points
!
! --- iwklm is the level above the height of the mountain.
! --- idxzb is the level of the dividing streamline.
! INITIALIZE DIVIDING STREAMLINE (DS) CONTROL VECTOR
!
        do i=1,npt
          iwklm(i) = 2
          kreflm(i) = 0
        enddo
!
!
! start lm mtn blocking (mb) section
!
!..............................
!..............................
!
!  (*j*)  11/03:  test upper limit on KMLL=km - 1
!      then do not need hncrit -- test with large hncrit first.
!       KMLL  = km / 2 ! maximum mtnlm height : # of vertical levels / 2
        KMLL = kmm1
! --- No mtn should be as high as KMLL (so we do not have to start at 
! --- the top of the model but could do calc for all levels).
!

        DO I = 1, npt
          j = ipt(i)
          ELEVMX(J) = min (ELEVMX(J) + sigfac * hprime(j), hncrit)
        ENDDO

        DO K = 1,KMLL
          DO I = 1, npt
            j = ipt(i)
! --- interpolate to max mtn height for index, iwklm(I) wk[gz]
! --- ELEVMX is limited to hncrit because to hi res topo30 orog.
            pkp1log =  phil(j,k+1) * ROG
            pklog =  phil(j,k) * ROG
            if ( ( ELEVMX(j) .le.  pkp1log ) .and.                      &
     &           ( ELEVMX(j) .ge.   pklog  ) ) THEN
! ---        wk for diags but can be saved and reused.  
               wk(i)  = G * ELEVMX(j) / ( phil(j,k+1) - phil(j,k) )
               iwklm(I)  =  MAX(iwklm(I), k+1 ) 

            endif
!
! ---        find at prsl levels large scale environment variables
! ---        these cover all possible mtn max heights
            VTJ(I,K)  = T1(J,K)  * (1.+FV*Q1(J,K))  ! virtual temperature
            VTK(I,K)  = VTJ(I,K) / PRSLK(J,K)       ! potential temperature
            RO(I,K)   = RDI * PRSL(J,K) / VTJ(I,K)  ! DENSITY (1.e-3 kg m^-3)

          ENDDO    !-- DO I = 1, npt
!
        ENDDO      !-- DO K = 1,KMLL
!
        klevm1 = KMLL - 1
        DO K = 1, klevm1  
          DO I = 1, npt
           j   = ipt(i)
            RDZ  = g   / ( phil(j,k+1) - phil(j,k) )
! ---                               Brunt-Vaisala Frequency
            BNV2LM(I,K) = (G+G) * RDZ * ( VTK(I,K+1)-VTK(I,K) )         &
     &                     / ( VTK(I,K+1)+VTK(I,K) )
            bnv2lm(i,k) = max( bnv2lm(i,k), bnv2min )

          ENDDO
        ENDDO
!
        DO I = 1, npt
          J   = ipt(i)
          DELKS(I)  = 1.0 / (PRSI(J,1) - PRSI(J,iwklm(i)))
          DELKS1(I) = 1.0 / (PRSL(J,1) - PRSL(J,iwklm(i)))
          UBAR (I)  = 0.0
          VBAR (I)  = 0.0
          ROLL (I)  = 0.0
          PE   (I)  = 0.0
          EK   (I)  = 0.0
          BNV2bar(I) = (PRSL(J,1)-PRSL(J,2)) * DELKS1(I) * BNV2LM(I,1)
        ENDDO
!
! --- find the dividing stream line height 
! --- starting from the level above the max mtn downward
! --- iwklm(i) is the k-index of mtn elevmx elevation
!
        DO Ktrial = KMLL, 1, -1
          DO I = 1, npt
             IF ( Ktrial .LT. iwklm(I) .and. kreflm(I) .eq. 0 ) then
                kreflm(I) = Ktrial
             ENDIF
          ENDDO
        ENDDO
!
! --- in the layer kreflm(I) to 1 find PE (which needs N, ELEVMX)
! ---  make averages, guess dividing stream (DS) line layer.
! ---  This is not used in the first cut except for testing and
! --- is the vert ave of quantities from the surface to mtn top.
!   

        DO I = 1, npt
          DO K = 1, Kreflm(I)
            J        = ipt(i)
            RDELKS     = DEL(J,K) * DELKS(I)
            UBAR(I)    = UBAR(I)  + RDELKS * U1(J,K) ! trial Mean U below 
            VBAR(I)    = VBAR(I)  + RDELKS * V1(J,K) ! trial Mean V below 
            ROLL(I)    = ROLL(I)  + RDELKS * RO(I,K) ! trial Mean RO below 
            RDELKS     = (PRSL(J,K)-PRSL(J,K+1)) * DELKS1(I)
            BNV2bar(I) = BNV2bar(I) + BNV2lm(I,K) * RDELKS
! --- these vert ave are for diags, testing and GWD to follow (*j*).

          ENDDO
        ENDDO

!
! --- integrate to get PE in the trial layer.
! --- Need the first layer where PE>EK - as soon as 
! --- IDXZB is not 0 we have a hit and Zb is found.
!
        DO I = 1, npt
          J = ipt(i)

          DO K = iwklm(I), 1, -1
            PHIANG   =  atan2(V1(J,K),U1(J,K))*RAD_TO_DEG
! --- THETA ranges from -+90deg |_ to the mtn "largest topo variations"
            ANG(I,K) = ( THETA(J) - PHIANG )
            if ( ANG(I,K) .gt.  90. ) ANG(I,K) = ANG(I,K) - 180.
            if ( ANG(I,K) .lt. -90. ) ANG(I,K) = ANG(I,K) + 180.
            ANG(I,K) = ANG(I,K)*DEG_TO_RAD
            UDS(I,K) = MAX(SQRT(U1(J,K)*U1(J,K)+V1(J,K)*V1(J,K)), minwnd)
! --- Test to see if we found Zb previously
            IF (IDXZB(I) .eq. 0 ) then
              PE(I) = PE(I) + BNV2lm(I,K) *                             &
     &           ( G * ELEVMX(J) - phil(J,K) ) *                        &
     &           ( PHII(J,K+1) - PHII(J,K) ) * ROG2

! --- Kinetic energy (KE): Wind projected on the line perpendicular to 
!     mtn range, U(Zb(K)).  KE is at the layer Zb
              UP(I) = UDS(I,K)*cos(ANG(I,K))
              EK(I) = 0.5*UP(I)*UP(I) 
! --- Dividing Stream lime  is found when PE =exceeds EK.
              IF ( PE(I) .ge.  EK(I) ) IDXZB(I) = K
! --- Then mtn blocked flow is between Zb=k(IDXZB(I)) and surface
!
            ENDIF     !-- IF (IDXZB(I) .eq. 0 ) then
          ENDDO       !-- DO K = iwklm(I), 1, -1
        ENDDO         !-- DO I = 1, npt
!
        DO I = 1, npt
          J    = ipt(i)
! --- Calc if N constant in layers (Zb guess) - a diagnostic only.
          ZBK(I) =  ELEVMX(J) - SQRT(UBAR(I)**2 + VBAR(I)**2)/BNV2bar(I)
        ENDDO
!
! --- The drag for mtn blocked flow
! 
        DO I = 1, npt
          J = ipt(i)
          ZLEN = 0.
          IF ( IDXZB(I) .gt. 0 ) then 
            DO K = IDXZB(I), 1, -1
              IF (PHIL(J,IDXZB(I)) > PHIL(J,K)) THEN
                ZLEN = SQRT( ( PHIL(J,IDXZB(I))-PHIL(J,K) ) /           &
     &                       ( PHIL(J,K ) + G * HPRIME(J) ) )
! --- lm eq 14:
                cosang=cos(ANG(I,K))
                sinang=sin(ANG(I,K))
                R = (cosang*cosang + GAMMA(J)*sinang*sinang) / &
     &              (GAMMA(J)*cosang*cosang + sinang*sinang)
! --- (negative of DB -- see sign at tendency); CDMB adjusts mountain blocking
! --- modified lm eq 15:
                DBTMP = 0.25 *  CDMB *                                  &
     &                  MAX( 2. - 1. / R, HZERO ) * SIGMA(J) *          &
     &                  MAX(HZERO, cosang, GAMMA(J)*sinang) *           &
     &                  ZLEN / hprime(J) 
                DB(I,K) = MAX(HZERO, DBTMP*UDS(I,K) )
              ENDIF        !-- IF (PHIL(J,IDXZB(I)) > PHIL(J,K) .AND. DEN > 0.) THEN
            ENDDO          !-- DO K = IDXZB(I), 1, -1
          endif
        ENDDO              !-- DO I = 1, npt
!
!.............................
!.............................
! end  mtn blocking section
!
      ENDIF      !-- IF ( NMTVR .eq. 14) then 
!
!.............................
!.............................
!
      KMPBL  = km / 2 ! maximum pbl height : # of vertical levels / 2

      DO K = 1,KM
        DO I =1,npt
          J         = ipt(i)
          VTJ(I,K)  = T1(J,K)  * (1.+FV*Q1(J,K))
          VTK(I,K)  = VTJ(I,K) / PRSLK(J,K)
          RO(I,K)   = RDI * PRSL(J,K) / VTJ(I,K)  ! DENSITY 
          TAUP(I,K) = 0.0
        ENDDO
      ENDDO

      DO K = 1,KMM1
        DO I =1,npt
          J         = ipt(i)
          TI        = 2.0 / (T1(J,K)+T1(J,K+1))
          TEM       = TI  / (PRSL(J,K)-PRSL(J,K+1))
          RDZ       = g   / (phil(j,k+1) - phil(j,k))
          TEM1      = U1(J,K) - U1(J,K+1)
          TEM2      = V1(J,K) - V1(J,K+1)
          DW2       = (TEM1*TEM1 + TEM2*TEM2)
          SHR2      = MAX(DW2,DW2MIN) * RDZ * RDZ
          BVF2      = G*(GOCP+RDZ*(VTJ(I,K+1)-VTJ(I,K))) * TI
!-- ri_n is Richardson Number, BNV2 is Brunt-Vaisala Frequency
          ri_n(I,K) = MAX(BVF2/SHR2,RIMIN)   ! Richardson number
          BNV2(I,K) = (G+G) * RDZ * (VTK(I,K+1)-VTK(I,K))               &
     &                            / (VTK(I,K+1)+VTK(I,K))
          bnv2(i,k) = max( bnv2(i,k), bnv2min )
        ENDDO     !-- DO K = 1,KMM1
      ENDDO       !-- DO I =1,npt
!
      do i=1,npt
        iwk(i) = 2
      enddo

      DO K=3,KMPBL
        DO I=1,npt
          j   = ipt(i)
          tem = (prsi(j,1) - prsi(j,k))
          if (tem .lt. dpmin) iwk(i) = k
        enddo
      enddo
!
      KBPS = 1
      KMPS = KM
      DO I=1,npt
        J         = ipt(i)
        kref(I)   = MAX(IWK(I), KPBL(J)+1 ) ! reference level 
        DELKS(I)  = 1.0 / (PRSI(J,1) - PRSI(J,kref(I)))
        DELKS1(I) = 1.0 / (PRSL(J,1) - PRSL(J,kref(I)))
        UBAR (I)  = 0.0
        VBAR (I)  = 0.0
        ROLL (I)  = 0.0
        KBPS      = MAX(KBPS,  kref(I))
        KMPS      = MIN(KMPS,  kref(I))
        BNV2bar(I) = (PRSL(J,1)-PRSL(J,2)) * DELKS1(I) * BNV2(I,1)
      ENDDO
!
      KBPSP1 = KBPS + 1
      KBPSM1 = KBPS - 1
      DO K = 1,KBPS
        DO I = 1,npt
          IF (K .LT. kref(I)) THEN
            J          = ipt(i)
            RDELKS     = DEL(J,K) * DELKS(I)
            UBAR(I)    = UBAR(I)  + RDELKS * U1(J,K)   ! Mean U below kref
            VBAR(I)    = VBAR(I)  + RDELKS * V1(J,K)   ! Mean V below kref
            ROLL(I)    = ROLL(I)  + RDELKS * RO(I,K)   ! Mean RO below kref
            RDELKS     = (PRSL(J,K)-PRSL(J,K+1)) * DELKS1(I)
            BNV2bar(I) = BNV2bar(I) + BNV2(I,K) * RDELKS
          ENDIF
        ENDDO
      ENDDO
!
!     FIGURE OUT LOW-LEVEL HORIZONTAL WIND DIRECTION AND FIND 'OA'
!
!             NWD  1   2   3   4   5   6   7   8
!              WD  W   S  SW  NW   E   N  NE  SE
!
      DO I = 1,npt
        J      = ipt(i)
        wdir   = atan2(UBAR(I),VBAR(I)) + pi
        idir   = mod(nint(fdir*wdir),mdir) + 1
        nwd    = nwdir(idir)
        OA(I)  = (1-2*INT( (NWD-1)/4 )) * OA4(J,MOD(NWD-1,4)+1)
        CLX(I) = CLX4(J,MOD(NWD-1,4)+1)
      ENDDO
!
!-----XN,YN            "LOW-LEVEL" WIND PROJECTIONS IN ZONAL
!                                    & MERIDIONAL DIRECTIONS
!-----ULOW             "LOW-LEVEL" WIND MAGNITUDE -        (= U)
!-----BNV2             BNV2 = N**2
!-----TAUB             BASE MOMENTUM FLUX
!-----= -(RO * U**3/(N*XL)*GF(FR) FOR N**2 > 0
!-----= 0.                        FOR N**2 < 0
!-----FR               FROUDE    =   N*HPRIME / U
!-----G                GMAX*FR**2/(FR**2+CG/OC)
!
!-----INITIALIZE SOME ARRAYS
!
      DO I = 1,npt
        XN(I)     = 0.0
        YN(I)     = 0.0
        TAUB (I)  = 0.0
        ULOW (I)  = 0.0
        DTFAC(I)  = 1.0
        ICRILV(I) = .FALSE. ! INITIALIZE CRITICAL LEVEL CONTROL VECTOR
!
!----COMPUTE THE "LOW LEVEL" WIND MAGNITUDE (M/S)
!
        ULOW(I) = MAX(SQRT(UBAR(I)*UBAR(I) + VBAR(I)*VBAR(I)), HONE)
        ULOI(I) = 1.0 / ULOW(I)
      ENDDO
!
      DO  K = 1,KMM1
        DO  I = 1,npt
          J            = ipt(i)
          VELCO(I,K)   = 0.5*ULOI(I)*( (U1(J,K)+U1(J,K+1))*UBAR(I)       &
     &                                +(V1(J,K)+V1(J,K+1))*VBAR(I) )
        ENDDO
      ENDDO
!
!   find the interface level of the projected wind where
!   low levels & upper levels meet above pbl
!
      do i=1,npt
        kint(i) = km
      enddo
      do k = 1,kmm1
        do i = 1,npt
          IF (K .GT. kref(I)) THEN
            if(velco(i,k) .lt. veleps .and. kint(i) .eq. km) then
              kint(i) = k+1
            endif
          endif
        enddo
      enddo
!  WARNING  KINT = KREF !!!!!!!!!
      do i=1,npt
        kint(i) = kref(i)
      enddo
!
!
      DO I = 1,npt
        J      = ipt(i)
        BNV    = SQRT( BNV2bar(I) )
        FR     = BNV     * ULOI(I) * min(HPRIME(J),hpmax)
        FR     = MIN(FR, FRMAX)
        XN(I)  = UBAR(I) * ULOI(I)
        YN(I)  = VBAR(I) * ULOI(I)
!
!     Compute the base level stress and store it in TAUB
!     CALCULATE ENHANCEMENT FACTOR, NUMBER OF MOUNTAINS & ASPECT
!     RATIO CONST. USE SIMPLIFIED RELATIONSHIP BETWEEN STANDARD
!     DEVIATION & CRITICAL HGT
!
        EFACT    = (OA(I) + 2.) ** (CEOFRC*FR)
        EFACT    = MIN( MAX(EFACT,EFMIN), EFMAX )
!
        COEFM    = (1. + CLX(I)) ** (OA(I)+1.)
!
        XLINV(I) = COEFM * CLEFF
!
        TEM      = FR    * FR * OC(J)
        GFOBNV   = GMAX  * TEM / ((TEM + CG)*BNV)  ! G/N0
!
        TAUB(I)  = XLINV(I) * ROLL(I) * ULOW(I) * ULOW(I)               &
     &           * ULOW(I)  * GFOBNV  * EFACT         ! BASE FLUX Tau0
!
        K        = MAX(1, kref(I)-1)
!-- Change in the GFS version of gwdps.f
        TEM      = MAX(VELCO(I,K)*VELCO(I,K), HE_1)
        SCOR(I)  = BNV2(I,K) / TEM  ! Scorer parameter below ref level
      ENDDO    !-- DO I = 1,npt
!                                                                       
!----SET UP BOTTOM VALUES OF STRESS
!
      DO K = 1, KBPS
        DO I = 1,npt
          IF (K .LE. kref(I)) TAUP(I,K) = TAUB(I)
        ENDDO
      ENDDO

!
!   Now compute vertical structure of the stress.
!
      DO K = KMPS, KMM1                   ! Vertical Level K Loop!
        KP1 = K + 1
        DO I = 1, npt
!
!-----UNSTABLE LAYER IF RI < RIC
!-----UNSTABLE LAYER IF UPPER AIR VEL COMP ALONG SURF VEL <=0 (CRIT LAY)
!---- AT (U-C)=0. CRIT LAYER EXISTS AND BIT VECTOR SHOULD BE SET (.LE.)
!
          IF (K .GE. kref(I)) THEN
            ICRILV(I) = ICRILV(I) .OR. ( ri_n(I,K) .LT. RIC)            &
     &                            .OR. (VELCO(I,K) .LE. 0.0)
          ENDIF
        ENDDO
!
        DO I = 1,npt
          IF (K .GE. kref(I))   THEN
!
            IF (.NOT.ICRILV(I) .AND. TAUP(I,K) .GT. 0.0 ) THEN
              TEMV = 1.0 / max(VELCO(I,K), HE_2)
!             IF (OA(I) .GT. 0. .AND.  PRSI(ipt(i),KP1).GT.RLOLEV) THEN
              IF (OA(I).GT.0. .AND. kp1 .lt. kint(i)) THEN
                SCORK   = BNV2(I,K) * TEMV * TEMV
                RSCOR   = MIN(HONE, SCORK / SCOR(I))
                SCOR(I) = SCORK
              ELSE 
                RSCOR   = 1.
              ENDIF
!
              BRVF = SQRT(BNV2(I,K))        ! Brunt-Vaisala Frequency
              TEM1 = XLINV(I)*(RO(I,KP1)+RO(I,K))*BRVF*0.5              &
     &                       * max(VELCO(I,K),HE_2)
              HD   = SQRT(TAUP(I,K) / TEM1)
              FRO  = BRVF * HD * TEMV
!
!    RIM is the  MINIMUM-RICHARDSON NUMBER BY SHUTTS (1985)
!
              TEM2   = SQRT(ri_n(I,K))
              TEM    = 1. + TEM2 * FRO
              RIM    = ri_n(I,K) * (1.-FRO) / (TEM * TEM)
!
!    CHECK STABILITY TO EMPLOY THE SATURATION HYPOTHESIS
!    OF LINDZEN (1981) EXCEPT AT TROPOSPHERIC DOWNSTREAM REGIONS
!
!                                       ----------------------
              IF (RIM .LE. RIC .AND.                                    &
     &           (OA(I) .LE. 0. .OR.  kp1 .ge. kint(i) )) THEN
                 TEMC = 2.0 + 1.0 / TEM2
                 HD   = VELCO(I,K) * (2.*SQRT(TEMC)-TEMC) / BRVF
                 TAUP(I,KP1) = TEM1 * HD * HD
              ELSE 
                 TAUP(I,KP1) = TAUP(I,K) * RSCOR
              ENDIF
              taup(i,kp1) = min(taup(i,kp1), taup(i,k))

            ENDIF    !-- IF (.NOT.ICRILV(I) .AND. TAUP(I,K) .GT. 0.0 ) THEN
          ENDIF      !-- IF (K .GE. kref(I))   THEN
        ENDDO        !-- DO I = 1,npt
      ENDDO          !-- DO K = KMPS, KMM1
!
      IF(LCAP .LE. KM) THEN

         DO KLCAP = LCAPP1, KM+1
            DO I = 1,npt
              SIRA          = PRSI(ipt(I),KLCAP) / PRSI(ipt(I),LCAP)
              TAUP(I,KLCAP) = SIRA * TAUP(I,LCAP)
            ENDDO
         ENDDO
      ENDIF
!
!     Calculate - (g/p*)*d(tau)/d(sigma) and Decel terms DTAUX, DTAUY
!
      DO K = 1,KM
        DO I = 1,npt
          TAUD(I,K) = G*(TAUP(I,K+1)-TAUP(I,K))/DEL(ipt(I),K)
        ENDDO
      ENDDO

!
!------LIMIT DE-ACCELERATION (MOMENTUM DEPOSITION ) AT TOP TO 1/2 VALUE
!------THE IDEA IS SOME STUFF MUST GO OUT THE TOP
!
      DO KLCAP = LCAP, KM
         DO I = 1,npt
            TAUD(I,KLCAP) = TAUD(I,KLCAP) * FACTOP
         ENDDO
      ENDDO
!
!------IF THE GRAVITY WAVE DRAG WOULD FORCE A CRITICAL LINE IN THE
!------LAYERS BELOW SIGMA=RLOLEV DURING THE NEXT DELTIM TIMESTEP,
!------THEN ONLY APPLY DRAG UNTIL THAT CRITICAL LINE IS REACHED.
!
      DO K = 1,KMM1
        DO I = 1,npt
           IF (K .GT. kref(I) .and. PRSI(ipt(i),K) .GE. RLOLEV) THEN
             IF(TAUD(I,K).NE.0.) THEN
               TEM = DELTIM * TAUD(I,K)
               DTFAC(I) = MIN(DTFAC(I),ABS(VELCO(I,K)/TEM))

             ENDIF
           ENDIF
        ENDDO
      ENDDO

      DO K = 1,KM
        DO I = 1,npt
          J          = ipt(i)
          TAUD(I,K)  = TAUD(I,K) * DTFAC(I)
          DTAUX      = TAUD(I,K) * XN(I)
          DTAUY      = TAUD(I,K) * YN(I)
! ---  lm mb (*j*)  changes overwrite GWD
          if ( K .lt. IDXZB(I) .AND. IDXZB(I) .ne. 0 ) then
            DBIM = DB(I,K) / (1.+DB(I,K)*DELTIM)

!dbg 20160708
if(DBIM<0.) write(0,"(a,3g12.4)") 'DBIM,DB,DELTIM=',DBIM,DB(I,K),DELTIM

            DBIM = MAX(DBIM, HZERO)
            A(J,K)  = - DBIM * V1(J,K) + A(J,K)
            B(J,K)  = - DBIM * U1(J,K) + B(J,K)
            DUsfc(J)   = DUsfc(J) - DBIM * V1(J,K) * DEL(J,K)
            DVsfc(J)   = DVsfc(J) - DBIM * U1(J,K) * DEL(J,K)
          else
            A(J,K)     = DTAUY     + A(J,K)
            B(J,K)     = DTAUX     + B(J,K)
            DUsfc(J)   = DUsfc(J)  + DTAUX * DEL(J,K)
            DVsfc(J)   = DVsfc(J)  + DTAUY * DEL(J,K)
          endif
        ENDDO      !-- DO I = 1,npt
      ENDDO        !-- DO K = 1,KM
!
      DO I = 1,npt
        J          = ipt(i)
        DUsfc(J) = ROGN * DUsfc(J)
        DVsfc(J) = ROGN * DVsfc(J)
      ENDDO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GWD_col
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_GWD
!
!-----------------------------------------------------------------------
