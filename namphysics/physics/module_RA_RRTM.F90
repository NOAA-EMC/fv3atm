!
      MODULE MODULE_RA_RRTM
!
!-----------------------------------------------------------------------
!
!***  THE RADIATION DRIVERS AND PACKAGES
!
!-----------------------------------------------------------------------
!
      USE MODULE_KINDS
!
      use physparam, only : icldflg, ioznflg, kind_phys, icmphys        &
     &                    , ivflip, lcrick, lcnorm, lnoprec, lsashal    &
     &                    , IALBflg
      use physcons,  only : con_eps, con_epsm1, con_fvirt               &
     &                    , con_t0c, con_ttp, con_g, con_rd

      USE MODULE_CONSTANTS, ONLY : R,CP,PI,EPSQ,STBOLT,EP_2
      USE MODULE_MP_FER_HIRES, ONLY : FPVS

      use module_radiation_driver_nmmb,  only : grrad_nmmb

      use module_radsw_parameters,  only : topfsw_type, sfcfsw_type
      use module_radlw_parameters,  only : topflw_type, sfcflw_type

      use module_radiation_clouds_nmmb,  only : NF_CLDS                 &
     &                                   , progcld1, progcld2, diagcld1 &
     &                                   , progcld8

!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: RRTM, RRTM_INIT
!
!-----------------------------------------------------------------------
!
!--- Used for Gaussian look up tables
!

!Moved to here from out of the RRTM routine to fix races found by Intel Inspector, jm 20131222,20131226
      LOGICAL       :: CNCLD=.TRUE.
      LOGICAL       :: OPER=.TRUE.
!$OMP THREADPRIVATE(CNCLD,OPER)

      REAL, PRIVATE,PARAMETER :: XSDmax=3.1, DXSD=.01
      INTEGER, PRIVATE,PARAMETER :: NXSD=XSDmax/DXSD
      REAL, DIMENSION(NXSD),PRIVATE,SAVE :: AXSD
      REAL, PRIVATE :: RSQR
      LOGICAL, PRIVATE, SAVE :: SDprint=.FALSE.

!-------------------------------
      INTEGER, SAVE, DIMENSION(3)     :: LTOP
      REAL,SAVE,DIMENSION(4) :: PTOPC
!--------------------------------
!
      REAL, PARAMETER ::         &
     &   RHgrd=1.00              & !--- RH (unitless) for onset of condensation
     &,  TRAD_ice=273.15-30.     & !--- Very tunable parameter
     &,  ABSCOEF_W=800.          & !--- Very tunable parameter
     &,  ABSCOEF_I=500.          & !--- Very tunable parameter
     &,  Qconv=0.1e-3            & !--- Very tunable parameter

     &,  CTauCW=ABSCOEF_W*Qconv  &
     &,  CTauCI=ABSCOEF_I*Qconv

!-- Set to TRUE to bogus in small amounts of convective clouds into the
!   input cloud calculations, but only if all of the following conditions 
!   are met:
!     (1) The maximum condensate mixing ratio is < QWmax.
!     (2) Only shallow convection is present, do not apply to deep convection.
!     (3) Only apply if the depth of shallow convection is between 
!         CU_DEEP_MIN (50 hPa) and CU_DEEP_MAX (200 hPa).
!     (4) Convective precipitation rate must be <0.01 mm/h.  
!
      LOGICAL, SAVE :: CUCLD=.FALSE.     &   ! was .TRUE.
     &                ,SUBGRID=.TRUE.
!
!-- After several tuning experiments, a value for QW_CU=0.003 g/kg should 
!   produce a cloud fraction of O(25%) and a SW reduction of O(100 W/m**2) 
!   for shallow convection with a maximum depth/thickness of O(200 hPa).
!-- QW_Cu=0.003 g/kg, which translates to a 3% cloud fraction in 
!   subroutine progcld2 at each model layer near line 960 in 
!   radiation_clouds.f, which translates to a O(25%) total cloud fraction
!   in the lower atmosphere (i.e., for "low-level" cloud fractions).
!
      REAL, PARAMETER :: QW_Cu=0.003E-3,QWmax=1.E-7,CUPPT_min=1.e-5   &
                        ,CU_DEEP_MIN=50.E2,CU_DEEP_MAX=200.E2

!--- set default quantities for new version of progcld2

      real (kind=kind_phys), parameter :: cclimit = 0.001, cclimit2=0.05
      real (kind=kind_phys), parameter :: recwat_def = 5.0    ! default liq radius to 5 microns at <0C
      real (kind=kind_phys), parameter :: recice_def = 10.0   ! default ice radius to 10 microns
      real (kind=kind_phys), parameter :: rerain_def = 100.0  ! default rain radius to 100 microns
      real (kind=kind_phys), parameter :: resnow_def = 50.0   ! default snow radius to 50 microns


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE RADIATION PACKAGE OPTIONS
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE RRTM (NTIMESTEP,DT,JDAT                                &
     &                    ,NPHS,GLAT,GLON                               &
     &                    ,NRADS,NRADL                                  &
     &                    ,P8W,P_PHY                                    & !rv prsi, prsl in Pa
     &                    ,T,Q,CW,O3                                    &
     &                    ,ALBEDO                                       &
     &                    ,ALBVB,ALBNB,ALBVD,ALBND                      & ! MODIS albedos
     &                    ,F_ICE,F_RAIN                                 &
     &                    ,QC,QI,QS,QR,QG,NI                            &
     &                    ,F_QC,F_QI,F_QS,F_QR,F_QG,F_NI                &
     &                    ,CLD_FRACTION                                 &
     &                    ,SM,CLDFRA                                    &
     &                    ,RLWTT,RSWTT                                  &
     &                    ,RLWIN,RSWIN                                  &
     &                    ,RSWINC,RLWINC,RSWOUC,RLWOUC                  &
     &                    ,RLWOUCtoa,RSWOUCtoa                          &
     &                    ,RSWOUT                                       &
     &                    ,RLWTOA,RSWTOA                                &
     &                    ,CZMEAN,SIGT4                                 &
     &                    ,CFRACL,CFRACM,CFRACH                         &
     &                    ,ACFRST                                       &
     &                    ,ACFRCV                                       &
     &                    ,CUPPT,SNOWC,SI                               & !was SNOW
     &                    ,HTOP,HBOT                                    &
     &                    ,TSKIN,Z0,SICE,F_RIMEF,MXSNAL,STDH,VVL        &
     &                    ,IMS,IME,JMS,JME                              &
     &                    ,ITS,ITE,JTS,JTE                              &
     &                    ,LM                                           &
     &                    ,SOLCON                                       &
     &                    ,MYPE )
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: IME,IMS,ITE,ITS                             &
     &                     ,JME,JMS,JTE,JTS                             &
     &                     ,LM,MYPE                                     &
     &                     ,NTIMESTEP                                   &
     &                     ,NPHS,NRADL,NRADS                            &
     &                     ,CLD_FRACTION
!
      INTEGER,INTENT(IN) :: JDAT(8)
!
      REAL,INTENT(IN) :: DT

      real (kind=kind_phys), INTENT(IN) :: SOLCON
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: CUPPT               &
                                                   ,GLAT,GLON           &
                                                   ,SM,SNOWC,SI

!-- MODIS albedos
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: ALBVB,ALBNB         & ! vis+uv & near IR beam
                                                   ,ALBVD,ALBND           ! vis+uv & near IR diffuse
!
      REAL,DIMENSION(IMS:IME,JMS:JME,LM),INTENT(IN) :: CW,O3,Q,T
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: F_ICE,F_RAIN
!
      REAL,DIMENSION(IMS:IME,1:LM+1),INTENT(IN) :: P8W
      REAL,DIMENSION(IMS:IME,1:LM)  ,INTENT(IN) :: P_PHY
!
!-- Update ALBEDO array if MODIS albedos are used
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: ACFRCV,ACFRST    &
                                                      ,RLWIN,RLWTOA     &
                                                      ,RSWIN,RSWOUT     &
                                                      ,HBOT,HTOP        &
                                                      ,ALBEDO           & ! allow to modify
                                                      ,RSWINC,RLWINC    &
                                                      ,RSWOUC,RLWOUC    &
                                                      ,RLWOUCtoa        &
                                                      ,RSWOUCtoa        &
                                                      ,RSWTOA
!
      REAL,DIMENSION(IMS:IME,JMS:JME,LM),INTENT(INOUT) :: RLWTT,RSWTT
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: CFRACH,CFRACL    &
                                                      ,CFRACM,CZMEAN    &
                                                      ,SIGT4
!
      LOGICAL,INTENT(IN) :: F_QC,F_QS,F_QI,F_QR,F_QG,F_NI

      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: QC,QS          &
     &                     ,QI,QR,QG,NI 
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: CLDFRA
!
       REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: TSKIN,Z0,SICE      &
                                                    ,MXSNAL,STDH        
!
       REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: F_RIMEF
       REAL,DIMENSION(IMS:IME,        1:LM),INTENT(IN) :: VVL
!
       real(kind=kind_phys),DIMENSION(IMS:IME,JMS:JME) :: COSZDG    ! future output
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!


!
      LOGICAL ::  LSSAV=.TRUE.
                ! logical flag for store 3-d cloud field
                ! ** need to be .TRUE. for non-zero FLUXR_V off GRRAD

      LOGICAL :: LPRNT=.FALSE.
!
      LOGICAL :: LSLWR, LSSWR

      INTEGER,PARAMETER :: NFLUXR=39

!==========================================================================
!  ---  constants
!  -- From eq. (5) on p. 2434 in McFarquhar & Heymsfield (1996)
!==========================================================================

      real (kind=kind_phys), parameter ::                               &
     &                 re_50C=1250.0/9.917, re_40C=1250.0/9.337,        &
     &                 re_30C=1250.0/9.208, re_20C=1250.0/9.387

!==========================================================================
!  Special for the lwrad to enhence the emissivity
!  it is similar to *CPATHFAC4LW to odcld in radlw  (Hsin-Mu Lin, 20140520)
!==========================================================================

      real(kind=kind_phys), PARAMETER :: CPATHFAC4LW=1.5
!
      real(kind=kind_phys), PARAMETER :: QMIN=1.0e-10, QME5=1.0e-7,       &
                                         QME6=1.0e-7
!
!-- WARNING: NTRAC must be large enough to account for 
!   different hydrometeor species +2, for ozone & aerosols
!
      INTEGER,PARAMETER :: NTRAC=9   ! GR1 dimension for ozone, aerosol, & clouds
      INTEGER,PARAMETER :: NTCW =3   ! ARRAY INDEX LOCATION FOR CLOUD CONDENSATE
      INTEGER,PARAMETER :: NCLDX=1   ! only used when ntcw .gt. 0

!
      INTEGER :: NUMX, NUMY, NFXR, KFLIP, I, L, J, K

      INTEGER :: ICWP, NTOZ
!
      real(kind=kind_phys),DIMENSION((ITE-ITS+1)) :: RTvR
      real(kind=kind_phys) ::  DTSW, DTLW, FHSWR, FHLWR, ARG_CW, QSTMP
      real(kind=kind_phys) ::  DUMMY
!
      real(kind=kind_phys),DIMENSION((ITE-ITS+1)) ::                          &
                             FLGMIN_L, CV, CVB, CVT, HPRIME_V, TSEA,      &
                             TISFC, FICE, ZORL, SLMSK, SNWDPH, SNCOVR,    &
                             SNOALB, ALVSF1, ALNSF1, ALVWF1, ALNWF1,      &
                             FACSF1, FACWF1, SFCNSW, SFCDSW, SFALB,       &
                             SFCDLW, TSFLW, TOAUSW, TOADSW, SFCCDSW,      &
                             TOAULW, SFCUSW, COSZEN_V, COSZDG_V,          &
                             SEMIS, XLAT, XLON, SINLAT, COSLAT,           &
                             CVB1, CVT1, tem1d, ES

                           !===================================
                           ! SEMIS: surface lw emissivity
                           !        is intended output in GLOOPR
                           !        ** not NMMB in RRTM driver
                           !===================================

      INTEGER, DIMENSION((ITE-ITS+1)) :: ICSDSW, ICSDLW

!---  variables of instantaneous calculated toa/sfc radiation fluxes
!
      type (topfsw_type), dimension((ITE-ITS+1)) :: TOPFSW
      type (sfcfsw_type), dimension((ITE-ITS+1)) :: SFCFSW

      type (topflw_type), dimension((ITE-ITS+1)) :: TOPFLW
      type (sfcflw_type), dimension((ITE-ITS+1)) :: SFCFLW

      real(kind=kind_phys),DIMENSION((ITE-ITS+1),LM) ::                     &
     &                      CLDCOV_V,PRSL,PRSLK,GT,GQ, F_ICEC,     &
                            F_RAINC,R_RIME,TAUCLOUDS,CLDF

      real(kind=kind_phys),DIMENSION((ITE-ITS+1),LM+1) :: PRSI, plvl
!
      real(kind=kind_phys),DIMENSION((ITE-ITS+1),5)  :: CLDSA_V
!
      real(kind=kind_phys),DIMENSION((ITE-ITS+1),NFLUXR) :: FLUXR_V
!
      real(kind=kind_phys),DIMENSION((ITE-ITS+1),LM,NTRAC) :: GR1   
!
      real(kind=kind_phys),DIMENSION((ITE-ITS+1),LM) :: SWH, HLW
!
      real(kind=kind_phys),DIMENSION((ITE-ITS+1),LM,NF_CLDS) :: clouds
!
      real(kind=kind_phys),DIMENSION((ITE-ITS+1),LM) ::                     &
                            plyr, qlyr, rhly, qstl, vvel, clw, tvly,     &
                            qc2d, qi2d, qs2d, ni2d,                     &
                            clwf, clw2, cldtot, tem2d,                  &
                            qcwat, qcice, qrain, rsden,                 &
                            cwp, cip, crp, csp, rew, rei, rer, res
!
      INTEGER, DIMENSION((ITE-ITS+1),3) :: mbota, mtopa
!
      INTEGER :: JDOY, JDAY, JDOW, MMM, MMP, MM, IRET, MONEND,          &
                 MON1, IS2, ISX, KPD9, IS1, NN, MON2, MON, IS,          &  
                 LUGB, LEN, JMSK, IMSK       
!
      REAL :: WV,QICE,QCLD,CLFR,ESAT,QSAT,RHUM,RHtot,ARG,SDM,           &
               PMOD,CONVPRATE,CLSTP,P1,P2,CC1,CC2,CLDMAX,CL1,CL2,       &
               CR1,DPCL,PRS1,PRS2,DELP,TCLD,CTau,CFSmax,CFCmax,         &
               CFRAVG,TDUM,CU_DEPTH
!
      INTEGER :: IXSD,NTSPH,NRADPP,NC,NMOD,LCNVT,LCNVB,NLVL,MALVL,      &
                 LLTOP,LLBOT,KBT2,KTH1,KBT1,KTH2,KTOP1,LM1,LL,LMP1,     &
                 JCX,LV,NP3D

!
      REAL, PARAMETER :: EPSQ1=1.E-5,EPSQ2=1.E-8,EPSO3=1.E-10,H0=0.,    &
                         H1=1.,HALF=.5,CUPRATE=24.*1000.,               &
                         HPINC=HALF*1.E1, CLFRmin=0.01, TAUCmax=4.161,  &
                         XSDmin=-XSDmax, DXSD1=-DXSD, STSDM=0.01,       & 
                         CVSDM=.04,DXSD2=HALF*DXSD,DXSD2N=-DXSD2,PCLDY=0.25
!
      REAL,DIMENSION(10),SAVE :: CC,PPT

! moved out of routine into module and made threadprivate (jm 20131226, see above)
!jm      LOGICAL, SAVE :: CNCLD=.TRUE.
!jm      LOGICAL, SAVE :: OPER=.TRUE.

!
      DATA CC/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/
      DATA PPT/0.,.14,.31,.70,1.6,3.4,7.7,17.,38.,85./
!
      REAL,DIMENSION(0:LM)  :: CLDAMT
!
      LOGICAL :: BITX,BITY,BITZ,BITW,BIT1,BIT2,NEW_CLOUD,CU_cloud((ITE-ITS+1))
!
      REAL :: CTHK(3)
      DATA CTHK/20000.0,20000.0,20000.0/
! 
      REAL,DIMENSION(ITS:ITE,JTS:JTE,3):: CLDCFR
      INTEGER,DIMENSION(ITS:ITE,JTS:JTE,3):: MBOT,MTOP

      REAL,DIMENSION(ITS:ITE,JTS:JTE):: CUTOP,CUBOT
      
      REAL,DIMENSION(ITS:ITE,JTS:JTE,LM) :: TauCI,CSMID,CCMID
!
      INTEGER,DIMENSION(ITS:ITE,JTS:JTE,LM+1) :: KTOP, KBTM

      REAL,DIMENSION(ITS:ITE,JTS:JTE,LM+1) :: CAMT
!
      INTEGER,DIMENSION(ITS:ITE,JTS:JTE) :: NCLDS, KCLD
!
      REAL,DIMENSION(ITS:ITE,JTS:JTE,LM) :: TAUTOTAL
!
      INTEGER :: NKTP, NBTM, NCLD, LML, ihr, imin, INDX_CW
!
      REAL :: CLFR1, TauC

      real (kind=kind_phys), parameter :: f24 = 24.0     ! hours/day
      real (kind=kind_phys) :: SFCALBEDO((ITE-ITS+1)), SMX((ITE-ITS+1)),        &
     &                         fhr, solhr, tem1, tem2, tem3,            &
     &                         clwmin, clwm, clwt, onemrh, value

      integer ii, nf
      integer,external :: omp_get_thread_num

!--------------------------------------------------------------------------------------------------
!
!***THIS SUBROUTINE SELECTS AND PREPARES THE NECESSARY INPUTS FOR GRRAD (GFS RRTM DRIVER)
!
!   GRRAD IS CALLED COLUMN BY COLUMN
!
!INPUTS/OUPUTS OF GRRAD: 
!    INPUT VARIABLES:                                                   !
!      PRSI  (LM+1)    : MODEL LEVEL PRESSURE IN CB (KPA)               !
!      PRSL  (LM)      : MODEL LAYER MEAN PRESSURE IN CB (KPA)          !
!      PRSLK (LM)      : Exner function (dimensionless)                 !
!      GT    (LM)      : MODEL LAYER MEAN TEMPERATURE IN K              !
!      GQ    (LM)      : LAYER SPECIFIC HUMIDITY IN GM/GM               !
!      GR1   (LM,NTRAC): TRACER ARRAY (OZONE, AEROSOL, Various Hydrometeors) !
!      VVL   (LM)      : LAYER MEAN VERTICAL VELOCITY IN CB/SEC         ! !rv now in Pa/s
!      SLMSK (1)       : SEA/LAND MASK ARRAY (SEA:0,LAND:1,SEA-ICE:2)   !
!      XLON,XLAT       : GRID LONGITUDE/LATITUDE IN RADIANS             !
!      TSEA  (1)       : SURFACE TEMPERATURE IN K                       !
!      SNWDPH (1)       : SNOW DEPTH WATER EQUIVALENT IN MM              !
!      SNCOVR(1)       : SNOW COVER IN FRACTION                         !
!      SNOALB(1)       : MAXIMUM SNOW ALBEDO IN FRACTION                !
!      ZORL  (1)       : SURFACE ROUGHNESS IN CM                        !
!      HPRIM_V (1)       : TOPOGRAPHIC STANDARD DEVIATION IN M            !
!      ALVSF1 (1)       : MEAN VIS ALBEDO WITH STRONG COSZ DEPENDENCY    !
!      ALNSF1 (1)       : MEAN NIR ALBEDO WITH STRONG COSZ DEPENDENCY    !
!      ALVWF1 (1)       : MEAN VIS ALBEDO WITH WEAK COSZ DEPENDENCY      !
!      ALNWF1 (1)       : MEAN NIR ALBEDO WITH WEAK COSZ DEPENDENCY      !
!      FACSF1 (1)       : FRACTIONAL COVERAGE WITH STRONG COSZ DEPENDEN  !
!      FACWF1 (1)       : FRACTIONAL COVERAGE WITH WEAK COSZ DEPENDENCY  !
!      FICE  (1)       : ICE FRACTION OVER OPEN WATER GRID              !
!      TISFC (1)       : SURFACE TEMPERATURE OVER ICE FRACTION          !
!      SOLCON          : SOLAR CONSTANT (SUN-EARTH DISTANT ADJUSTED)    !
!
!-- Following 5 quantities are defined within a local 'tile':
!      SINLAT_t        : SINE OF LATITUDE                               !
!      COSLAT_t        : COSINE OF LATITUDE                             !
!      XLON_t          : LONGITUDE                                      !
!      COSZEN_t        : MEAN COS OF ZENITH ANGLE OVER RAD CALL PERIOD  !
!      COSZDG_t        : MEAN COS OF ZENITH ANGLE OVER RAD CALL PERIOD  !
!
!      CV    (1)       : FRACTION OF CONVECTIVE CLOUD                   ! !not used
!      CVT, CVB (1)    : CONVECTIVE CLOUD TOP/BOTTOM PRESSURE IN CB     ! !not used
!      IOVRSW/IOVRLW   : CONTROL FLAG FOR CLOUD OVERLAP (SW/LW RAD)     !
!                        =0 RANDOM OVERLAPPING CLOUDS                   !
!                        =1 MAX/RAN OVERLAPPING CLOUDS                  !
!      F_ICEC (LM)     : FRACTION OF CLOUD ICE  (IN FERRIER SCHEME)     !
!      F_RAINC(LM)     : FRACTION OF RAIN WATER (IN FERRIER SCHEME)     !
!      RRIME  (LM)     : MASS RATIO OF TOTAL TO UNRIMED ICE ( >= 1 )    !
!      FLGMIN_L(1)     : MINIMIM LARGE ICE FRACTION                     !
!                        =8 THOMPSON MICROPHYSICS SCHEME                ! G. Thompson 23Feb2013
!      NTCW            : =0 NO CLOUD CONDENSATE CALCULATED              !
!                        >0 ARRAY INDEX LOCATION FOR CLOUD CONDENSATE   !
!      NCLDX           : ONLY USED WHEN NTCW .GT. 0                     !
!      NTOZ            : =0 CLIMATOLOGICAL OZONE PROFILE                !
!                        >0 INTERACTIVE OZONE PROFILE                   ! !does not work currently
!      NTRAC           : DIMENSION VERIABLE FOR ARRAY GR1               !
!      NFXR            : SECOND DIMENSION OF INPUT/OUTPUT ARRAY FLUXR   !
!      DTLW, DTSW      : TIME DURATION FOR LW/SW RADIATION CALL IN SEC  !
!      LSSAV           : LOGICAL FLAG FOR STORE 3-D CLOUD FIELD         !
!      LM              : VERTICAL LAYER DIMENSION                       !
!      MYPE            : CONTROL FLAG FOR PARALLEL PROCESS              !
!      LPRNT           : CONTROL FLAG FOR DIAGNOSTIC PRINT OUT          !
!      TAUCLOUDS(LM)   : CLOUD OPTICAL DEPTH FROM NMMB (ferrier+bmj)    ! !new
!      CLDF(LM)        : CLOUD FRACTION FROM NMMB (ferrier+bmj)         ! !new
!                                                                       !
!    OUTPUT VARIABLES:                                                  !
!      SWH (LM)       : TOTAL SKY SW HEATING RATE IN K/SEC              !
!      SFCNSW(1)      : TOTAL SKY SURFACE NET SW FLUX IN W/M**2         !
!      SFCDSW(1)      : TOTAL SKY SURFACE DOWNWARD SW FLUX IN W/M**2    !
!      SFALB (1)      : MEAN SURFACE DIFFUSED ALBEDO                    !
!      HLW (LM)       : TOTAL SKY LW HEATING RATE IN K/SEC              !
!      SFCDLW(1)      : TOTAL SKY SURFACE DOWNWARD LW FLUX IN W/M**2    !
!      TSFLW (1)      : SURFACE AIR TEMP DURING LW CALCULATION IN K     !
!
!      TOAUSW (IM)    : TOTAL SKY TOA UPWARD SW FLUX IN W/M**2         ! !new
!      TOADSW (IM)    : TOTAL SKY TOA DOWNWARD SW FLUX IN W/M**2       ! !new
!      SFCCDSW(IM)    : CLEAR SKY SURFACE SW DOWNWARD FLUX IN W/M**2   ! !new
!      TOAULW (IM)    : TOTAL SKY TOA LW FLUX W/M**2                   ! !new
!      SFCUSW (IM)    : TOTAL SKY SURFACE SW UPWARD FLUX IN W/M**2     ! !new
!                                                                       !
!    INPUT AND OUTPUT VARIABLES:                                        !
!      FLUXR_V (IX,NFXR) : TO SAVE 2-D FIELDS                           !
!                          (bucket)                                     !
!      CLDSA_V(IX,5)     : TO SAVE 2-D CLOUD FRACTION. L/M/H/TOT/BL     !
!                          (instantaneous)                              !
!      CLDCOV_V(IX,LM)   : TO SAVE 3-D CLOUD FRACTION                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!SELECT OPTIONS IN GRRAD
!
      NFXR=NFLUXR    ! second dimension of input/output array fluxr (FLUXR_V)

      ICSDSW(:)=0    ! auxiliary special cloud related array for SW
                     ! *** not used in this version of code ***
                     ! can be any value at this moment
      ICSDLW(:)=0    ! auxiliary special cloud related array for LW
                     ! *** not used in this version of code ***
                     ! can be any value at this moment

      ICWP = icldflg
      NTOZ = ioznflg

      np3d = ICMPHYS
      LMP1 = LM+1

!------------------------------
! for np3d=5 (Lin, 20150601)
!------------------------------

      IF (ICMPHYS == 5 ) THEN
         ICWP = -1
      ENDIF
!
!=========================================================================
!
      IF (ICWP/=-1 .AND. CNCLD) THEN
         CNCLD=.FALSE.        !-- used when ICWP=1, 0
      ENDIF
!
!--- Cloud water index for the GR1 array
!
      IF (F_NI) THEN
        INDX_CW=4  !-- Thompson 
      ELSE
        INDX_CW=3  !-- All others (as of Oct 2014)
      ENDIF
!
!CLOUDS
!
!----------------------CONVECTION--------------------------------------
!  NRADPP IS THE NUMBER OF TIME STEPS TO ACCUMULATE CONVECTIVE PRECIP
!     FOR RADIATION
!   NOTE: THIS WILL NOT WORK IF NRADS AND NRADL ARE DIFFERENT UNLESS
!         THEY ARE INTEGER MULTIPLES OF EACH OTHER
!  CLSTP IS THE NUMBER OF HOURS OF THE ACCUMULATION PERIOD
!
      NTSPH=NINT(3600./DT)
      NRADPP=MIN(NRADS,NRADL)
      CLSTP=1.0*NRADPP/NTSPH
      CONVPRATE=CUPRATE/CLSTP

      IF (ICWP>0 .AND. CUCLD) CONVPRATE=1000./CLSTP    !-- convert to mm/h
!
      LM1=LM-1
!
      DO J=JTS,JTE
      DO I=ITS,ITE
        DO K=1,LM
          CCMID(I,J,K)=0.
          CSMID(I,J,K)=0.
        ENDDO
      ENDDO
      ENDDO
!
! --- initialize for non Thompson cloud fraction (used only in gfdl type)
!     for thompson cloud fraction, "CLDFRC" is direct INPUT
!     
      IF (CLD_FRACTION==0) THEN
         DO K=1,LM
         DO J=JTS,JTE
         DO I=ITS,ITE
            CLDFRA(I,J,K)=0.
         ENDDO
         ENDDO
         ENDDO
      ENDIF

! ---- 

      DO K=1,LM
      DO J=JTS,JTE
      DO I=ITS,ITE
         TAUTOTAL(I,J,K)=0.
      ENDDO
      ENDDO
      ENDDO
!
      DO J=JTS,JTE
      DO I=ITS,ITE
          CFRACH(I,J)=0.
          CFRACL(I,J)=0.
          CFRACM(I,J)=0.
          CZMEAN(I,J)=0.
          SIGT4(I,J)=0.
      ENDDO
      ENDDO
!
      DO K=1,3
      DO J=JTS,JTE
      DO I=ITS,ITE
        CLDCFR(I,J,K)=0.
        MTOP(I,J,K)=0
        MBOT(I,J,K)=0
      ENDDO
      ENDDO
      ENDDO
!
      DO J=JTS,JTE
      DO I=ITS,ITE
       CUTOP(I,J)=LM+1-HTOP(I,J)
       CUBOT(I,J)=LM+1-HBOT(I,J)
      ENDDO
      ENDDO
!      
!-----------------------------------------------------------------------
!---  COMPUTE GRID-SCALE CLOUD COVER FOR RADIATION  (Ferrier, Nov '04)
!
!--- Assumes Gaussian-distributed probability density functions (PDFs) for
!    total relative humidity (RHtot) within the grid for convective and
!    grid-scale cloud processes.  The standard deviation of RHtot is assumed
!    to be larger for convective clouds than grid-scale (stratiform) clouds.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      ICWP_Test: IF (ICWP==-1) THEN   !-- *** Start of old NAM/GFDL cloud inputs ***
!-----------------------------------------------------------------------
       DO J=JTS,JTE
       DO I=ITS,ITE 
!
        DO 255 L=1,LM
!
          WV=MAX(EPSQ,Q(I,J,L))/(1.-MAX(EPSQ,Q(I,J,L)))   !-- Water vapor mixing ratio
          QICE=MAX(QS(I,J,L),0.)                          !-- Ice mixing ratio
          QCLD=QICE+MAX(QS(I,J,L),0.)                     !-- Total cloud water + ice mixing ratio
!rv------------------------------------
!rv   This should be temporary fix!!!!!
!rv   New (currently operational) calculation of cloud fraction is
!rv   causing different results with different decomposition
!rv   We should find cause of this!!!!!
!rv------------------------------------
          OPER_flag: IF (OPER) THEN
!rv------------------------------------
!-- From model tuning experiments vs CLAVR grid-to-grid verification:
!-- 100% cloud fractions at 0.01 g/kg (1.e-5 kg/kg) cloud mixing ratios
!-- 10% cloud fractions at 1.e-4 g/kg (1.e-7 kg/kg) cloud mixing ratios
!-- 1% cloud fractions at 1.e-6 g/kg (1.e-9 kg/kg) cloud mixing ratios
!
            CLFR=MIN(H1, MAX(H0,1.e5*QCLD))
            CLFR=SQRT(CLFR)
            IF (CLFR>=CLFRmin) CSMID(I,J,L)=CLFR
!rv------------------------------------
          else OPER_flag
!rv------------------------------------

!
            IF (QCLD .LE. EPSQ) GO TO 255                               !--- Skip if no condensate is present
            CLFR=H0
!
            WV=MAX(EPSQ,Q(I,J,L))/(1.-MAX(EPSQ,Q(I,J,L)))
!
!--- Saturation vapor pressure w/r/t water ( >=0C ) or ice ( <0C )
!
            ESAT=1000.*FPVS(T(I,J,L))                                   !--- Saturation vapor pressure (Pa)
            ESAT=MIN(ESAT, 0.99*P_PHY(I,L) )                            !--- Put limits on ESAT
            QSAT=EP_2*ESAT/(P_PHY(I,L)-ESAT)                            !--- Saturation mixing ratio
!
            RHUM=WV/QSAT                                                !--- Relative humidity
!
!--- Revised cloud cover parameterization (temporarily ignore rain)
!
            RHtot=(WV+QCLD)/QSAT                                        !--- Total relative humidity
!
            LCNVT=NINT(CUTOP(I,J))
            LCNVT=MIN(LM,LCNVT)
            LCNVB=NINT(CUBOT(I,J))
            LCNVB=MIN(LM,LCNVB)
            IF (L.GE.LCNVT .AND. L.LE.LCNVB) THEN
               SDM=CVSDM
            ELSE
               SDM=STSDM
            ENDIF
            ARG=(RHtot-RHgrd)/SDM
            IF (ARG.LE.DXSD2 .AND. ARG.GE.DXSD2N) THEN
               CLFR=HALF
            ELSE IF (ARG .GT. DXSD2) THEN
               IF (ARG .GE. XSDmax) THEN
                  CLFR=H1
               ELSE
                  IXSD=INT(ARG/DXSD+HALF)
                  IXSD=MIN(NXSD, MAX(IXSD,1))
                  CLFR=HALF+AXSD(IXSD)
               ENDIF              !--- End IF (ARG .GE. XSDmax)
            ELSE
               IF (ARG .LE. XSDmin) THEN
                  CLFR=H0
               ELSE
                  IXSD=INT(ARG/DXSD1+HALF)
                  IXSD=MIN(NXSD, MAX(IXSD,1))
                  CLFR=HALF-AXSD(IXSD)
                  IF (CLFR .LT. CLFRmin) CLFR=H0
               ENDIF        !--- End IF (ARG .LE. XSDmin)
            ENDIF           !--- IF (ARG.LE.DXSD2 .AND. ARG.GE.DXSD2N)
            CSMID(I,J,L)=CLFR
!rv------------------------------------
          endif  OPER_flag
!rv------------------------------------
!
255     CONTINUE         !--- End DO L=1,LM

       ENDDO ! End DO I=ITS,ITE
       ENDDO ! End DO J=JTS,JTE

!***********************************************************************
!******************  END OF GRID-SCALE CLOUD FRACTIONS  ****************


!***********************************************************************
!---  COMPUTE CONVECTIVE CLOUD COVER FOR RADIATION
!
!--- The parameterization of Slingo (1987, QJRMS, Table 1, p. 904) is
!    used for convective cloud fraction as a function of precipitation
!    rate.  Cloud fractions have been increased by 20% for each rainrate
!    interval so that shallow, nonprecipitating convection is ascribed a
!    constant cloud fraction of 0.1  (Ferrier, Feb '02).
!***********************************************************************
!
       GFDL_Conv: IF (CNCLD) THEN

        DO J=JTS,JTE
         DO I=ITS,ITE
!
!***  CLOUD TOPS AND BOTTOMS COME FROM CUCNVC
!     Convective clouds need to be at least 2 model layers thick
!
          IF (CUBOT(I,J)-CUTOP(I,J) .GT. 1.0) THEN
!--- Compute convective cloud fractions if appropriate  (Ferrier, Feb '02)
            CLFR=CC(1)
            PMOD=CUPPT(I,J)*CONVPRATE
            IF (PMOD .GT. PPT(1)) THEN
              DO NC=1,10
                IF(PMOD.GT.PPT(NC)) NMOD=NC
              ENDDO
              IF (NMOD .GE. 10) THEN
                CLFR=CC(10)
              ELSE
                CC1=CC(NMOD)
                CC2=CC(NMOD+1)
                P1=PPT(NMOD)
                P2=PPT(NMOD+1)
                CLFR=CC1+(CC2-CC1)*(PMOD-P1)/(P2-P1)
              ENDIF      !--- End IF (NMOD .GE. 10) ...
              CLFR=MIN(H1, CLFR)
            ENDIF        !--- End IF (PMOD .GT. PPT(1)) ...
!
!***  ADD LVL TO BE CONSISTENT WITH OTHER WORKING ARRAYS
!
            LCNVT=NINT(CUTOP(I,J))
            LCNVT=MIN(LM,LCNVT)
            LCNVB=NINT(CUBOT(I,J))
            LCNVB=MIN(LM,LCNVB)
!
!--- Build in small amounts of subgrid-scale convective condensate
!    (simple assumptions), but only if the convective cloud fraction
!    exceeds that of the grid-scale cloud fraction
!
            DO L=LCNVT,LCNVB
              ARG=MAX(H0, H1-CSMID(I,J,L))
              CCMID(I,J,L)=MIN(ARG,CLFR)
            ENDDO           !--- End DO LL=LCNVT,LCNVB
          ENDIF             !--- IF (CUBOT(I,J)-CUTOP(I,J) .GT. 1.0) ...
         ENDDO               ! End DO I=ITS,ITE
        ENDDO                ! End DO J=JTS,JTE
       ENDIF  GFDL_Conv      !--- End IF (CNCLD) ...
!
!*********************************************************************
!***************  END OF CONVECTIVE CLOUD FRACTIONS  *****************
!*********************************************************************
!***
!*** INITIALIZE ARRAYS FOR USES LATER
!***

       DO I=ITS,ITE
       DO J=JTS,JTE
!
         LML=LM
!***
!*** NOTE: LAYER=1 IS THE SURFACE, AND LAYER=2 IS THE FIRST CLOUD
!***       LAYER ABOVE THE SURFACE AND SO ON.
!***
         KTOP(I,J,1)=LM+1
         KBTM(I,J,1)=LM+1
         CAMT(I,J,1)=1.0
         KCLD(I,J)=2
!
         DO 510 L=2,LM+1
           CAMT(I,J,L)=0.0
           KTOP(I,J,L)=1
           KBTM(I,J,L)=1
  510    CONTINUE
!### End changes so far
!***
!*** NOW CALCULATE THE AMOUNT, TOP, BOTTOM AND TYPE OF EACH CLOUD LAYER
!*** CLOUD TYPE=1: STRATIFORM CLOUD
!***       TYPE=2: CONVECTIVE CLOUD
!*** WHEN BOTH CONVECTIVE AND STRATIFORM CLOUDS EXIST AT THE SAME POINT,
!*** SELECT CONVECTIVE CLOUD WITH THE HIGHER CLOUD FRACTION.
!*** CLOUD LAYERS ARE SEPARATED BY TOTAL ABSENCE OF CLOUDINESS.
!*** NOTE: THERE IS ONLY ONE CONVECTIVE CLOUD LAYER IN ONE COLUMN.
!*** KTOP AND KBTM ARE THE TOP AND BOTTOM OF EACH CLOUD LAYER IN TERMS
!*** OF MODEL LEVEL.
!***
         NEW_CLOUD=.TRUE.
!
      DO L=2,LML
        LL=LML-L+1                                  !-- Model layer
        CLFR=MAX(CCMID(I,J,LL),CSMID(I,J,LL))       !-- Cloud fraction in layer
        CLFR1=MAX(CCMID(I,J,LL+1),CSMID(I,J,LL+1))  !-- Cloud fraction in lower layer
!-------------------
        IF (CLFR .GE. CLFRMIN) THEN
!--- Cloud present at level
          IF (NEW_CLOUD) THEN
!--- New cloud layer
            IF(L==2.AND.CLFR1>=CLFRmin)THEN
              KBTM(I,J,KCLD(I,J))=LL+1
              CAMT(I,J,KCLD(I,J))=CLFR1
            ELSE
              KBTM(I,J,KCLD(I,J))=LL
              CAMT(I,J,KCLD(I,J))=CLFR
            ENDIF
            NEW_CLOUD=.FALSE.
          ELSE
!--- Existing cloud layer
            CAMT(I,J,KCLD(I,J))=AMAX1(CAMT(I,J,KCLD(I,J)), CLFR)
          ENDIF        ! End IF (NEW_CLOUD .EQ. 0) ...
        ELSE IF (CLFR1 .GE. CLFRMIN) THEN
!--- Cloud is not present at level but did exist at lower level, then ...
          IF (L .EQ. 2) THEN
!--- For the case of ground fog
           KBTM(I,J,KCLD(I,J))=LL+1
           CAMT(I,J,KCLD(I,J))=CLFR1
          ENDIF
          KTOP(I,J,KCLD(I,J))=LL+1
          NEW_CLOUD=.TRUE.
          KCLD(I,J)=KCLD(I,J)+1
          CAMT(I,J,KCLD(I,J))=0.0
        ENDIF
!-------------------
      ENDDO      !--- End DO L loop
!***
!*** THE REAL NUMBER OF CLOUD LAYERS IS (THE FIRST IS THE GROUND;
!*** THE LAST IS THE SKY):
!***
      NCLDS(I,J)=KCLD(I,J)-2
      NCLD=NCLDS(I,J)
!***
!***  NOW CALCULATE CLOUD RADIATIVE PROPERTIES
!***
      IF(NCLD.GE.1)THEN
!***
!*** NOTE: THE FOLLOWING CALCULATIONS, THE UNIT FOR PRESSURE IS MB!!!
!***
        DO NC=2,NCLD+1
!
        TauC=0.    !--- Total optical depth for each cloud layer (solar & longwave)
        NKTP=LM+1
        NBTM=0
        BITX=CAMT(I,J,NC).GE.CLFRMIN
        NKTP=MIN(NKTP,KTOP(I,J,NC))
        NBTM=MAX(NBTM,KBTM(I,J,NC))
!
        DO LL=NKTP,NBTM
          L=NBTM-LL+NKTP 
          IF(LL.GE.KTOP(I,J,NC).AND.LL.LE.KBTM(I,J,NC).AND.BITX)THEN
            PRS1=P8W(I,L)*0.01 
            PRS2=P8W(I,L+1)*0.01
            DELP=PRS2-PRS1
!
            CTau=0.
!-- For crude estimation of convective cloud optical depths
            IF (CCMID(I,J,L) .GE. CLFRmin) THEN
              IF (T(I,J,L) .GE. TRAD_ice) THEN
                CTau=CTauCW            !--- Convective cloud water
              ELSE
                CTau=CTauCI            !--- Convective ice
              ENDIF
            ENDIF
!
!-- For crude estimation of grid-scale cloud optical depths
!
!--   => The following 2 lines were intended to reduce cloud optical depths further
!        than what's parameterized in the NAM and what's theoretically justified
            CTau=CTau+ABSCOEF_W*QC(I,J,L)+ABSCOEF_I*QS(I,J,L)

            TAUTOTAL(I,J,L)=CTau*DELP                          !Total model level cloud optical depth
            CLDFRA(I,J,L)=MAX(CCMID(I,J,LL),CSMID(I,J,LL))     !Cloud fraction at model level           
            TauC=TauC+DELP*CTau                                !Total cloud optical depth as in GFDL
!
          ENDIF      !--- End IF(LL.GE.KTOP(I,NC) ....
        ENDDO        !--- End DO LL
!
      ENDDO
!
      ENDIF       ! NCLD.GE.1
!
      ENDDO  !  DO I=ITS,ITE
      ENDDO  !  DO J=JTS,JTE
!-----------------------------------------------------------------------
      ENDIF  ICWP_Test   !*** End of Old NAM/GFDL cloud inputs ***
!-----------------------------------------------------------------------

      FHSWR=(NRADS*DT)/3600.        ! [h]
      FHLWR=(NRADL*DT)/3600.        ! [h]
      DTLW =(NRADL*DT)              ! [s]
      DTSW =(NRADS*DT)              ! [s]
      LSSWR=MOD(NTIMESTEP,NRADS)==0
      LSLWR=MOD(NTIMESTEP,NRADL)==0

!==========================================================================
!  Similar to GFS "gloopr.f" line #370,  #413
!  The following block is from old "radiation_astronomy_nmmb.f"
!==========================================================================

      ihr   = JDAT(5)
      imin  = JDAT(6)

!  --- ...  hour of forecast time

      !  solhr = mod( float(ihr), f24 )    ! previous version

      !=== the new calculatuion will eliminate the time lag due to
      !    "jdate(5)" handled by ESMF  (201208)

      fhr = float(ihr)+float(imin)/60.
      solhr = mod( fhr, f24 )

!..........................................................................
!
!==========================================================================
! Main domain loop: calling grrad
!==========================================================================
!     

      DO J=JTS,JTE  !start grrad loop column by column



       ! if ( GLON(I,J) >= 0.0 ) then
           XLON(1:(ITE-ITS+1)) = GLON(ITS:ITE,J)
       ! else
       !    XLON(1) = GLON(I,J) + PI        ! if in -pi->+pi convert to 0->2pi
       ! endif

       XLAT(1:(ITE-ITS+1)) = GLAT(ITS:ITE,J)

       SINLAT(1:(ITE-ITS+1)) = SIN ( XLAT(1:(ITE-ITS+1)) )
       COSLAT(1:(ITE-ITS+1)) = COS ( XLAT(1:(ITE-ITS+1)) )

       TSEA(1:(ITE-ITS+1)) = TSKIN(ITS:ITE,J)
       TISFC(1:(ITE-ITS+1))= TSKIN(ITS:ITE,J)                  ! change later if necessary
       ZORL(1:(ITE-ITS+1)) = Z0(ITS:ITE,J)*100.d0
       SNWDPH(1:(ITE-ITS+1))=SI(ITS:ITE,J)                    ! snwdph[mm]
       SNCOVR(1:(ITE-ITS+1))=SNOWC(ITS:ITE,J)
       SNOALB(1:(ITE-ITS+1))=MXSNAL(ITS:ITE,J)
       HPRIME_V(1:(ITE-ITS+1))=STDH(ITS:ITE,J)

       WHERE (SICE(ITS:ITE,J).GT.0.5)              ! slmsk - ocean  - 0
         SLMSK(1:(ITE-ITS+1))= 2.0d0                   !         land   - 1
         FICE(1:(ITE-ITS+1))=SICE(ITS:ITE,J)           ! change this later
       ELSEWHERE                                   !         seaice - 2
         SLMSK(1:(ITE-ITS+1))= 1.0d0-SM(ITS:ITE,J)     !
         FICE(1:(ITE-ITS+1))= 0.0d0                    ! change this later
       ENDWHERE
!
!---
!
      FLGMIN_L(1:(ITE-ITS+1))= 0.20d0 ! --- for ferrier

      CV (1:(ITE-ITS+1))=0.d0         ! not in use
      CVB(1:(ITE-ITS+1))=0.d0         ! not in use
      CVT(1:(ITE-ITS+1))=0.d0         ! not in use

      PRSI(1:(ITE-ITS+1),1)=P8W(ITS:ITE,1)/1000.                                ! [kPa]
!
      DO L=1,LM
        PRSI(1:(ITE-ITS+1),L+1)=P8W(ITS:ITE,L+1)/1000.                          ! (pressure on interface) [kPa]
        PRSL(1:(ITE-ITS+1),L)=P_PHY(ITS:ITE,L)/1000.                            ! (pressure on mid-layer) [kPa] 
        PRSLK(1:(ITE-ITS+1),L)=(PRSL(1:(ITE-ITS+1),L)*0.01d0)**(R/CP)
        RTvR(1:(ITE-ITS+1))=1./(R*(Q(ITS:ITE,J,L)*0.608+1.- &
                   CW(ITS:ITE,J,L))*T(ITS:ITE,J,L))
        GT(1:(ITE-ITS+1),L)=T(ITS:ITE,J,L)
        GQ(1:(ITE-ITS+1),L)=Q(ITS:ITE,J,L)
!
!--- GR1(:,:,1) - ozone
!    GR1(:,:,2) - reserved for prognostic aerosols in the future
!    GR1(:,:,3) - total condensate
!    GR1(:,:,4-9) - hydrometeor species from Thompson scheme
!
        DO ii=1,NTRAC
          GR1(1:(ITE-ITS+1),L,ii)=0.d0
        ENDDO
!
        IF (NTOZ>0) GR1(1:(ITE-ITS+1),l,1)=MAX(O3(ITS:ITE,J,L),EPSO3)
!
        CLDCOV_V(1:(ITE-ITS+1),L)=0.d0                     ! used for prognostic cloud
        TAUCLOUDS(1:(ITE-ITS+1),L)=TAUTOTAL(ITS:ITE,J,L)    ! CLOUD OPTICAL DEPTH (ICWP==-1)
        CLDF(1:(ITE-ITS+1),L)=CLDFRA(ITS:ITE,J,L)           ! CLOUD FRACTION

        GR1(1:(ITE-ITS+1),L,3)=CW(ITS:ITE,J,L)              ! total condensate
!
!----------
!
        thompson_test: IF (F_NI) THEN
!
!----------
!-- IF F_NI=true, then this means the microphysics is Thompson only.
!   Other progcld"X" drivers must be introduced into grrad_nmmb.f
!   in order to use GR1(:,:,N) where N>=4.  (BSF, Oct 2014)
!
!-- Warnings from Thompson:
!.. This is an awful way to deal with different physics having different
!.. number of species.  Something must eventually be done to resolve this
!.. section to be more flexible.  For now, we are directly passing each
!.. species in the water array into the GR1 array for use in the RRTM
!.. radiation scheme, but that requires a priori knowledge of which species
!.. is which index number at some later time.  This is far from optimal,
!.. but we proceed anyway.  Future developers be careful.
!.. If the WATER species include separate hydrometeor species, then
!.. fill in other elements even if unused.  Thompson microphysics
!.. will utilize certain elements when computing cloud optical depth.
!----------
!
           GR1(1:(ITE-ITS+1),L,4)=QC(ITS:ITE,J,L)
           GR1(1:(ITE-ITS+1),L,5)=QI(ITS:ITE,J,L)
           GR1(1:(ITE-ITS+1),L,6)=QS(ITS:ITE,J,L)
           GR1(1:(ITE-ITS+1),L,7)=QR(ITS:ITE,J,L)
           GR1(1:(ITE-ITS+1),L,8)=QG(ITS:ITE,J,L)
           GR1(1:(ITE-ITS+1),L,9)=NI(ITS:ITE,J,L)
!----------
        ELSE  thompson_test        ! for non-Thompson microphysics
!----------
           F_ICEC (1:(ITE-ITS+1),L)=max(0.0, min(1.0, F_ICE (ITS:ITE,J,L) ) )
           F_RAINC(1:(ITE-ITS+1),L)=max(0.0, min(1.0, F_RAIN(ITS:ITE,J,L) ) )
           R_RIME (1:(ITE-ITS+1),L)=max(1.0, F_RIMEF(ITS:ITE,J,L) )
!----------
        ENDIF  thompson_test
!----------

      ENDDO
!
!
      subgrid_cloud: IF (SUBGRID) THEN
!
!-- Build in tiny amounts of subgrid-scale cloud when no cloud is
!   present and RH > 95%.
!-- Note GR1(ii,L,3) is total condensate for all microphysics schemes
!
        thompson_testx: IF (F_NI) THEN
!
!-- Build in tiny amounts of subgrid-scale cloud for the Thompson scheme
!
          DO L=1,LM
          DO I=1,(ITE-ITS+1)
            WV=GQ(I,L)/(1.-GQ(I,L))                !- Water vapor mixing ratio
            TCLD=REAL(GT(I,L))                     !- Temperature (deg K)
            ESAT=FPVS(TCLD)                        !- Saturation vapor pressure (kPa)
            P1=REAL(PRSL(I,L))                     !- Pressure (kPa)
            ESAT=MIN(ESAT, 0.99*P1)                !- Limit saturation vapor pressure

IF(P1<1.E-2) WRITE(6,"(a,3i4,2g11.4)") 'I,J,L,PRSL,E_sat=',I,J,L,P1,ESAT   !dbg

            QSAT=EP_2*ESAT/(P1-ESAT)               !- Saturation mixing ratio
            RHUM=WV/QSAT                           !- Relative humidity
            IF (GR1(I,L,3)<EPSQ .AND. RHUM>0.95) THEN
              ARG=MIN(0.01, RHUM-0.95)*QSAT
              GR1(I,L,3)=MIN(0.01E-3, ARG)
              IF (TCLD>TRAD_ICE) THEN
                GR1(I,L,4)=GR1(I,L,3)
              ELSE
                GR1(I,L,5)=GR1(I,L,3)
              ENDIF
            ENDIF      !- IF (GR1(ii,L,3)<EPSQ ...
          ENDDO        !- DO I
          ENDDO        !- DO L
!
        ELSE  thompson_testx
!
!-- Build in tiny amounts of subgrid-scale cloud for other microphysics schemes
!
          DO L=1,LM 
          DO I=1,(ITE-ITS+1)
            WV=GQ(I,L)/(1.-GQ(I,L))                !- Water vapor mixing ratio
            TCLD=REAL(GT(I,L))                     !- Temperature (deg K)
            ESAT=FPVS(TCLD)                        !- Saturation vapor pressure (kPa)
            P1=REAL(PRSL(I,L))                     !- Pressure (kPa)
            ESAT=MIN(ESAT, 0.99*P1)                !- Limit saturation vapor pressure

IF(P1<1.E-2) WRITE(6,"(a,3i4,2g11.4)") 'I,J,L,PRSL,E_sat=',I,J,L,P1,ESAT   !dbg

            QSAT=EP_2*ESAT/(P1-ESAT)               !- Saturation mixing ratio
            RHUM=WV/QSAT                           !- Relative humidity
            IF (GR1(I,L,3)<EPSQ .AND. RHUM>0.95) THEN
              ARG=MIN(0.01, RHUM-0.95)*QSAT
              GR1(I,L,3)=MIN(0.01E-3, ARG)
              IF (TCLD>TRAD_ICE) THEN
                F_ICEC(I,L)=0.
              ELSE
                F_ICEC(I,L)=1.
              ENDIF
              F_RAINC(I,L)=0.
              R_RIME(I,L)=1.
            ENDIF
          ENDDO
          ENDDO

        ENDIF  thompson_testx

      ENDIF  subgrid_cloud
!
!-- Bogus in tiny amounts of shallow convection, but only if there are no
!   grid-scale clouds nor convective precipitation present.  Arrays CUTOP,
!   CUBOT are flipped with 1 at the top & LM at the surface (BSF, 7/18/2012)
!-- There are extra, nested IF statements to filter conditions as an extra
!   layer of caution.  
!
      CU_cloud=.FALSE.
      CU_Bogus1: IF (CUCLD) THEN
       DO I=ITS,ITE
         ii = i-ITS+1
         LCNVT=MIN(LM, NINT(CUTOP(I,J)) )   !-- Convective cloud top
         LCNVB=MIN(LM, NINT(CUBOT(I,J)) )   !-- Convective cloud base
         CU_DEPTH=0.
         CU_Index: IF (LCNVB-LCNVT>1) THEN
            CU_DEPTH=1000.*(PRSL(ii,LCNVB)-PRSL(ii,LCNVT))   !- Pa
            CU_Deep: IF (CU_DEPTH>=CU_DEEP_MIN .AND. CU_DEPTH<=CU_DEEP_MAX) THEN
               QCLD=MAXVAL( GR1(ii,1:LM,3) )  !- Maximum *total condensate*
               PMOD=CUPPT(I,J)*CONVPRATE
               CU_Clds: IF (QCLD<QWmax .AND. PMOD<=CUPPT_min) THEN
                  CU_cloud(ii)=.TRUE.
                  DO L=LCNVT,LCNVB
                    GR1(ii,L,3)=GR1(ii,L,3)+QW_Cu  !- for *cloud water*
                    IF (F_NI) GR1(ii,L,INDX_CW)=GR1(ii,L,INDX_CW)+QW_Cu    !- *total condensate* in Thompson
                  ENDDO
               ENDIF CU_Clds
            ENDIF CU_Deep
         ENDIF CU_Index
       ENDDO
      ENDIF CU_Bogus1
!
      DO NC=1,5
        CLDSA_V(1:(ITE-ITS+1),NC)=0.d0                 !used for prognostic cloud
      ENDDO
      DO NC=1,NFLUXR
        FLUXR_V(1:(ITE-ITS+1),NC)=0.d0                 !used for prognostic cloud
      ENDDO
!
!---

!#######################################################################
!=======================================================================
!
!=== Proposed cloud properties for radiation (moved from grrad)
!

!  --- ...  prepare atmospheric profiles for radiation input
!           convert pressure unit from Pa to mb
!=============================================================
! ** NOTE: For NMMB, the conversion is "from cb to mb"
!=============================================================

      do k = 1, LM
        do i = 1, (ITE-ITS+1)
          plvl(i,k) = 10.0 * prsi(i,k)   ! cb to mb
          plyr(i,k) = 10.0 * prsl(i,k)   ! cb to mb
        enddo

!  --- ...  compute relative humidity

          ! Note for qs (**Saturation specific humidity):
          !   ** prsl unit in NMMB is kPa

! note this loop won't vectorize because of the call to fpvs

        do i = 1, (ITE-ITS+1)
           TCLD=REAL(GT(I,K))
           ES(i) = min( prsl(i,k), fpvs(TCLD) )    ! fpvs in KPa
        enddo

!DEC$ SIMD

        do i = 1, (ITE-ITS+1)
          QSTMP = max( QMIN, con_eps*ES(i)/(PRSL(I,K) + con_epsm1*ES(i)) )
          rhly(i,k) = max( 0.0, min( 1.0, max(QMIN, GQ(i,k))/QSTMP ) )
          qstl(i,k) = QSTMP
        enddo
      enddo

      do i = 1, (ITE-ITS+1)
        plvl(i,LMP1) = 10.0 * prsi(i,LMP1)  ! cb to mb
      enddo

      do i = 1, (ITE-ITS+1)
         tem1d(i) = QME6
      enddo

!
! --- virtual temp, tvly
!
      if (ivflip == 0) then               ! input data from toa to sfc
        do k = 1, LM
          do i = 1, (ITE-ITS+1)
            qlyr(i,k) = max( tem1d(i), GQ(i,k) )
            tem1d(i)  = min( QME5, qlyr(i,k) )
            tvly(i,k) = GT(i,k) * (1.0 + con_fvirt*qlyr(i,k))    ! virtual temp in K
          enddo
        enddo
      else                               ! input data from sfc to toa
        do k = LM, 1, -1
          do i = 1, (ITE-ITS+1)
            qlyr(i,k) = max( tem1d(i), GQ(i,k) )
            tem1d(i)  = min( QME5, qlyr(i,k) )
            tvly(i,k) = GT(i,k) * (1.0 + con_fvirt*qlyr(i,k))     ! virtual temp in K
          enddo
        enddo
      endif                              ! end_if_ivflip

!############################################################
!  --- Start Prognostic and Diagnostic cloud scheme choice
!############################################################

      if (ntcw > 0) then                   ! prognostic cloud scheme

        do k = 1, LM
          do i = 1, (ITE-ITS+1)
            clw(i,k) = 0.0
          enddo

          do jcx = 1, ncldx
            lv = ntcw + jcx - 1
            do i = 1, (ITE-ITS+1)
               clw(i,k) = clw(i,k) + GR1(i,k,lv)    ! cloud condensate amount
            enddo
          enddo
        enddo

        do k = 1, LM
          do i = 1, (ITE-ITS+1)
            if ( clw(i,k) < EPSQ ) clw(i,k) = 0.0
          enddo
        enddo

        if (np3d == 4) then              ! zhao/moorthi's prognostic cloud scheme

          call progcld1                                                 &
!  ---  inputs:
     &     ( plyr,plvl,gt,tvly,qlyr,qstl,rhly,clw,                      &
     &       xlat,xlon,slmsk,                                           &
     &       (ITE-ITS+1), LM, LMP1,                                         &
!  ---  outputs:
     &       clouds,cldsa_v,mtopa,mbota                                 &
     &      )


        elseif (np3d == 3) then          ! ferrier's microphysics

      !****************************************************************
      !========================  20161012 =============================
      ! --- Proposed cloud properties (moved from radiation_clouds)
      !

          do nf=1,NF_CLDS
          do k = 1, LM
          do i = 1, (ITE-ITS+1)
            clouds(i,k,nf) = 0.0
          enddo
          enddo
          enddo

          do k = 1, LM
          do i = 1, (ITE-ITS+1)
            cldtot(i,k) = 0.0
          enddo
          enddo

       !-- lcrick

          if ( lcrick ) then
            do i = 1, (ITE-ITS+1)
              clwf(i,1)  = 0.75*clw(i,1)  + 0.25*clw(i,2)
              clwf(i,LM) = 0.75*clw(i,LM) + 0.25*clw(i,LM-1)
            enddo

            do k = 2, LM-1
            do i = 1, (ITE-ITS+1)
              clwf(i,K)=0.25*clw(i,k-1) + 0.5*clw(i,k) + 0.25*clw(i,k+1)
            enddo
            enddo

          else

            do k = 1, LM
            do i = 1, (ITE-ITS+1)
              clwf(i,k) = clw(i,k)
            enddo
            enddo
          endif

!  ---  separate cloud condensate into liquid, ice, and rain types, and
!       save the liquid+ice condensate in array clw2 for later
!       calculation of cloud fraction

          do k = 1, LM
          do i = 1, (ITE-ITS+1)
            tem2d (i,k) = GT(i,k) - con_t0c
          enddo
          enddo

          do k = 1, LM
          do i = 1, (ITE-ITS+1)
            if (tem2d(i,k) > -40.0) then
              qcice(i,k) = clwf(i,k) * F_ICEC(i,k)
              tem1       = clwf(i,k) - qcice(i,k)
              qrain(i,k) = tem1 * F_RAINC(i,k)
              qcwat(i,k) = tem1 - qrain(i,k)
              clw2 (i,k) = qcwat(i,k) + qcice(i,k)
            else
              qcice(i,k) = clwf(i,k)
              qrain(i,k) = 0.0
              qcwat(i,k) = 0.0
              clw2 (i,k) = clwf(i,k)
            endif
          enddo
          enddo

!=======================================================
!--- Get cloud fraction at each grod point : cldtot
!=======================================================

!-----------------------------------------------------------------------
          cldfracs:  if (cld_fraction==0) then   !** default non Thompson cloud fraction
!-----------------------------------------------------------------------

          if ( ivflip == 0 ) then       ! input data from toa to sfc

            clwmin = 0.0

            if (.not. lsashal) then
              do k = LM, 1, -1
              do i = 1, (ITE-ITS+1)
                clwt = 1.0e-8 * (plyr(i,k)*0.001)

                if (clw2(i,k) > clwt) then
!
!-- The following are valid at 1000 hPa, actual values are normalized by pressure
!-- 100% cloud fractions at 0.01 g/kg cloud mixing ratios
!-- 10% cloud fractions at 0.001 g/kg cloud mixing ratios
!-- 1% cloud fractions at 0.0001 g/kg cloud mixing ratios
!
                  tem1 = 1.0e5*clw2(i,k)*plyr(i,k)*0.001

                  cldtot(i,k) = min(1.0, tem1)
                endif
              enddo
              enddo
            else
              do k = LM, 1, -1
              do i = 1, (ITE-ITS+1)
                clwt = 2.0e-6 * (plyr(i,k)*0.001)

                if (clw2(i,k) > clwt) then
                  onemrh= max( 1.e-10, 1.0-rhly(i,k) )
                  clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

                  tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)    !jhan
                  tem1  = 100.0 / tem1

                  value = max( min( tem1*(clw2(i,k)-clwm), 50.0 ), 0.0 )
                  tem2  = sqrt( sqrt(rhly(i,k)) )

                  cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
                endif
              enddo
              enddo
            endif

          else                          ! input data from sfc to toa

            clwmin = 0.0e-6

            if (.not. lsashal) then
              do k = 1, LM
              do i = 1, (ITE-ITS+1)
                clwt = 2.0e-6 * (plyr(i,k)*0.001)

                if (clw2(i,k) > clwt) then
                  onemrh= max( 1.e-10, 1.0-rhly(i,k) )
                  clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

                  tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
                  tem1  = 2000.0 / tem1

                  value = max( min( tem1*(clw2(i,k)-clwm), 50.0 ), 0.0 )
                  tem2  = sqrt( sqrt(rhly(i,k)) )

                  cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
                endif
              enddo
              enddo
            else
              do k = 1, LM
              do i = 1, (ITE-ITS+1)
                clwt = 2.0e-6 * (plyr(i,k)*0.001)

                if (clw2(i,k) > clwt) then
                  onemrh= max( 1.e-10, 1.0-rhly(i,k) )
                  clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

                  tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)   !jhan
                  tem1  = 100.0 / tem1

                  value = max( min( tem1*(clw2(i,k)-clwm), 50.0 ), 0.0 )
                  tem2  = sqrt( sqrt(rhly(i,k)) )

                  cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
                endif
              enddo
              enddo
            endif

          endif                                ! end_if_ivflip

       !----------------------------------------
          else         !  Thompson cloud fraction

            do k = 1, LM
            do i = 1, (ITE-ITS+1)
              cldtot(i,k) = cldf(i,k)
            enddo
            enddo
       !----------------------------------------

!-----------------------------------------------------------------------
          endif cldfracs   !*** End of choice of cloud fraction ***
!-----------------------------------------------------------------------

!=========================================================
!-- get optical path and effective radius of hydrometeors
!=========================================================

!  --- initilization

          do k = 1, LM
          do i = 1, (ITE-ITS+1)
            cwp   (i,k) = 0.0
            cip   (i,k) = 0.0
            crp   (i,k) = 0.0
            csp   (i,k) = 0.0
            rew   (i,k) = recwat_def    ! default cloud water radius
            rei   (i,k) = recice_def    ! default cloud ice radius
            rer   (i,k) = rerain_def    ! default rain radius
            res   (i,k) = resnow_def    ! default snow radius
          enddo
          enddo

!====================================
!-- BSF 20120319      call  rsipath2

          call  rsipath2_tmp                                            &
!  ---  inputs:
     &     ( plyr, plvl, GT, qlyr, qcwat, qcice, qrain, R_RIME,         &
     &       (ITE-ITS+1), LM, ivflip, FLGMIN_L,                             &
!  ---  outputs:
     &       cwp, cip, crp, csp, rew, rer, res, rsden                   &
     &     )

!  =====
          do k = 1, LM
          do i = 1, (ITE-ITS+1)
            if (cldtot(i,k) < cclimit) then
              cldtot(i,k) = 0.0
              cwp(i,k)    = 0.0
              cip(i,k)    = 0.0
              crp(i,k)    = 0.0
              csp(i,k)    = 0.0
            endif
          enddo
          enddo

       !-- When lnoprec = .true. snow/rain has no impact on radiation

          if ( lnoprec ) then
            do k = 1, LM
            do i = 1, (ITE-ITS+1)
              crp(i,k) = 0.0
              csp(i,k) = 0.0
            enddo
            enddo
          endif

       !-- lcnorm

          if ( lcnorm ) then
            do k = 1, LM
            do i = 1, (ITE-ITS+1)
              if (cldtot(i,k) >= cclimit) then
                tem1 = 1.0 / max(cclimit2, cldtot(i,k))
                cwp(i,k) = cwp(i,k) * tem1
                cip(i,k) = cip(i,k) * tem1
                crp(i,k) = crp(i,k) * tem1
                csp(i,k) = csp(i,k) * tem1
              endif
            enddo
            enddo
          endif

       !-- get tem2d

          if ( ivflip == 0 ) then          ! input data from toa to sfc
            do k = 1, LM
            do i = 1, (ITE-ITS+1)
              tem2d(i,k) = (con_g * plyr(i,k))                          &
     &                   / (con_rd* (plvl(i,k+1) - plvl(i,k)))
            enddo
            enddo
          else                             ! input data from sfc to toa
            do k = 1, LM
            do i = 1, (ITE-ITS+1)
              tem2d(i,k) = (con_g * plyr(i,k))                          &
     &                   / (con_rd* (plvl(i,k) - plvl(i,k+1)))
            enddo
            enddo
          endif                            ! end_if_ivflip


!  --- Get effective ice cloud droplet radius

          do k = 1, LM
          do i = 1, (ITE-ITS+1)
            tem2 = cip(i,k)

            if (tem2 > 0.0) then

!-- tem3 is ice water content (IWC) in g m^-3

              tem3 = tem2d(i,k) * tem2 / tvly(i,k)

!-- From eq. (5) on p. 2434 in McFarquhar & Heymsfield (1996).
!  1) Based on ice clouds from tropical convective outflows
!  2) Formulation follows Fu (1996) with calculations consistent
!     with the generalized effective size for hexagonal ice
!     crystals, Dge.  See comment below regarding 0.6495 factor.

              tem1 = GT(i,k) - con_ttp
              if (tem1 < -50.0) then
                tem2 = re_50C*tem3**0.109
              elseif (tem1 < -40.0) then
                tem2 = re_40C*tem3**0.08
              elseif (tem1 < -30.0) then
                tem2 = re_30C*tem3**0.055
              else
                tem2 = re_20C*tem3**0.031
              endif

!-- 0.6495 is used to convert from Dge to re following eq. (3.12) in Fu (1996)
!   since Re values are multiplied by 1.5396 in subroutine cldprop in
!   radsw_main.f

              rei(i,k)   = max(recice_def, tem2*0.6495)

            endif
          enddo
          enddo

          do k = 1, LM
          do i = 1, (ITE-ITS+1)
            clouds(i,k,1) = cldtot(i,k)
            clouds(i,k,2) = cwp(i,k)
            clouds(i,k,3) = rew(i,k)
!
!--- Combine the radiative impacts of snow and ice to be treated as
!    cloud ice, but with corrections made to rei above (Oct 2013).
!
            clouds(i,k,4) = cip(i,k)
            clouds(i,k,5) = rei(i,k)
            clouds(i,k,6) = crp(i,k)
            clouds(i,k,7) = rer(i,k)
          enddo
          enddo

      !
      ! --- end of cloud properties (moved from radiation_clouds)
      !================================================================
      !****************************************************************

          call progcld2                                                 &
!  ---  inputs:
     &     ( plyr, xlat, (ITE-ITS+1), LM, cldtot,                           &
!  ---  outputs:
     &       cldsa_v,mtopa,mbota                                        &
     &      )

        elseif (np3d == 5) then         ! nmmb ferrier+bmj  (GFDL type diagnostic)
!
! ======================================================================
!  The original GFS RRTM will use "call progcld3".
! ** After merging the nmmb & GFS RRTM, "np3d=5" is temporarily used
!    for the GFDL diagnostic cloud type which is in the original meso
!    nmmb                                    (Hsin-mu Lin 12/28/2011)
! ======================================================================
!
!         call progcld3                                                 &
!  ---  inputs:
!    &     ( plyr,plvl,gt,qlyr,qstl,rhly,clw,                           &
!    &       xlat,xlon,slmsk, F_ICEC,F_RAINC,R_RIME,FLGMIN_L,           &
!    &       (ITE-ITS+1), LM, LMP1, iflip, iovrsw,                          &
!  ---  outputs:
!    &       clouds,cldsa_v,mtopa,mbota                                 &
!    &      )

          do k = 1, LM                                          ! GFDL type
            do i = 1, (ITE-ITS+1)                                   ! GFDL type
              clouds(i,k,1) = cldf(i,k)                         ! GFDL type
              clouds(i,k,2) = tauclouds(i,k)                    ! GFDL type
              clouds(i,k,3) = 0.99                              ! GFDL type
              clouds(i,k,4) = 0.84                              ! GFDL type
            enddo                                               ! GFDL type
          enddo                                                 ! GFDL type

!..This is far from optimal!  The specific array element numbers are
!.. potentially flexible for cloud water, cloud ice, snow mixing ratios
!.. and cloud ice number concentration.  Therefore, numerical values of
!.. 4, 5, 6, 9 are hard-wired when they should be made flexible.

        elseif (np3d == 8) then            ! Thompson microphysics

          do k = 1, LM
            do i = 1, (ITE-ITS+1)
              qc2d(i,k) = gr1(i,k,4)
              qi2d(i,k) = gr1(i,k,5)
              qs2d(i,k) = gr1(i,k,6)

              if (qc2d(i,k) .LE. EPSQ) qc2d(i,k)=0.0
              if (qi2d(i,k) .LE. EPSQ) qi2d(i,k)=0.0
              if (qs2d(i,k) .LE. EPSQ) qs2d(i,k)=0.0

              ni2d(i,k) = MAX(1.E-8, gr1(i,k,9))
            enddo
          enddo

          call progcld8 (plyr,plvl, gt, qlyr, qc2d, qi2d, qs2d, ni2d,   &
     &                   xlat, (ITE-ITS+1), LM, LMP1,                       &
     &                   cldf, cld_fraction,                            &
     &                   clouds, cldsa_v, mtopa, mbota)

        endif                              ! end if_np3d
!
      else                                 ! diagnostic cloud scheme

        do i = 1, (ITE-ITS+1)
          cvt1(i) = 10.0 * cvt(i)
          cvb1(i) = 10.0 * cvb(i)
        enddo

        do k = 1, LM
          do i = 1, (ITE-ITS+1)
            vvel(i,k) = 0.01 * vvl (i,k) !rv conv Pa to mb
          enddo
        enddo

!  ---  compute diagnostic cloud related quantities

        call diagcld1                                                   &
!  ---  inputs:
     &     ( plyr,plvl,gt,rhly,vvel,cv,cvt1,cvb1,                       &
     &       xlat,xlon,slmsk,                                           &
     &       (ITE-ITS+1), LM, LMP1,                                     &
!  ---  outputs:
     &       clouds,cldsa_v,mtopa,mbota                                 &
     &      )

      endif                                ! end_if_ntcw

      do k = 1, LM
        do i = 1, (ITE-ITS+1)
          CLDCOV_V(i,k) = clouds(i,k,1)
        enddo
      enddo
!
!=== end of Proposed cloud properties for radiation (moved from grrad)
!
!=======================================================================
!#######################################################################


      SFCALBEDO(1:(ITE-ITS+1)) = ALBEDO(ITS:ITE,J)
      SMX(1:(ITE-ITS+1)) = SM(ITS:ITE,J)

!--- MODIS albedos

      IF (IALBflg == 1) THEN
         ALVSF1(1:(ITE-ITS+1))=ALBVB(ITS:ITE,J)   !- vis+uv beam (direct) albedo at 60 deg
         ALNSF1(1:(ITE-ITS+1))=ALBNB(ITS:ITE,J)   !- near IR beam (direct) albedo at 60 deg
         ALVWF1(1:(ITE-ITS+1))=ALBVD(ITS:ITE,J)   !- vis+uv diffuse albedo
         ALNWF1(1:(ITE-ITS+1))=ALBND(ITS:ITE,J)   !- near IR diffuse albedo
      ELSE
         ALVSF1(1:(ITE-ITS+1))=SFCALBEDO(1:(ITE-ITS+1))    !- Matthews old broadband albedo
         ALNSF1(1:(ITE-ITS+1))=SFCALBEDO(1:(ITE-ITS+1))
         ALVWF1(1:(ITE-ITS+1))=SFCALBEDO(1:(ITE-ITS+1))
         ALNWF1(1:(ITE-ITS+1))=SFCALBEDO(1:(ITE-ITS+1))
      ENDIF

      FACSF1(1:(ITE-ITS+1))=1.-SMX(1:(ITE-ITS+1))     !- 1 if land (SMX=0), 0 if sea (SMX=1)
      FACWF1(1:(ITE-ITS+1))=0.                        !- new code sums FACSF1+FACWF1 over land

!  --- ...  calling radiation driver

      call grrad_nmmb                                                   &
!!  ---  inputs:
           ( PRSI,PLYR,PLVL,PRSLK,GT,GQ,GR1,SLMSK,RHLY,                 &
             QLYR,TVLY,                                                 &
             XLON,XLAT,TSEA,SNWDPH,SNCOVR,SNOALB,ZORL,HPRIME_V,         &
             ALVSF1,ALNSF1,ALVWF1,ALNWF1,FACSF1,FACWF1,                 &  ! processed inside grrad
             SFCALBEDO,SMX,                                             &  ! input for albedo cal
             FICE,TISFC,                                                &
             SINLAT,COSLAT,SOLHR, JDAT, SOLCON,                         &
             FHSWR ,NRADS,                                              &  ! extra input
             ICSDSW,ICSDLW,NTOZ,NTRAC,NFXR,                             &
             CPATHFAC4LW,                                               &  ! enhance factor of cloud depth for LW
             DTLW,DTSW,LSSWR,LSLWR,LSSAV,                               &
             (ITE-ITS+1), (ITE-ITS+1), LM, MYPE, LPRNT, 0, 0,                   &
             CLOUDS,CLDSA_V,mtopa,mbota,                                &  !! clouds properties for radiation
!!  ---  outputs:
             SWH,TOPFSW,SFCFSW,SFALB,COSZEN_V,COSZDG_V,                 &
             HLW,TOPFLW,SFCFLW,TSFLW,SEMIS,                             &
!!  ---  input/output:
             FLUXR_V                                                    &
!! ---  optional outputs:
 !           ,HTRSWB,HTRLWB                                              &
           )

      COSZDG(ITS:ITE,J) = COSZDG_V (1:(ITE-ITS+1))
      CZMEAN(ITS:ITE,J) = COSZEN_V (1:(ITE-ITS+1))

      DO L=1,LM
        RLWTT(ITS:ITE,J,L)=HLW(1:(ITE-ITS+1),L)
        RSWTT(ITS:ITE,J,L)=SWH(1:(ITE-ITS+1),L)
      ENDDO

!=========================================================
! modify this section by using TOPFSW,SFCFSW,TOPFLW,SFCFLW
! instead of TOAUSW,TOADSW,SFCCDSW,TOAULW,SFCUSW
!=========================================================

      ! RLWIN(I,J)=SFCDLW(1)
      ! RSWIN(I,J)=SFCDSW(1)
      ! RSWINC(I,J)=SFCCDSW(1)
      ! RSWOUT(I,J)=RSWIN(I,J)*SFALB(1)
      ! RLWTOA(I,J)=TOAULW(1)
      ! RSWTOA(I,J)=TOAUSW(1)

      ! RSWOUT(I,J)=RSWIN(I,J)*SFALB(1)

      RLWIN(ITS:ITE,J) =SFCFLW(1:(ITE-ITS+1))%dnfxc
      RSWIN(ITS:ITE,J) =SFCFSW(1:(ITE-ITS+1))%dnfxc
      RSWOUT(ITS:ITE,J)=SFCFSW(1:(ITE-ITS+1))%upfxc

      RSWINC(ITS:ITE,J)=SFCFSW(1:(ITE-ITS+1))%dnfx0
      RLWINC(ITS:ITE,J)=SFCFLW(1:(ITE-ITS+1))%dnfx0
      RSWOUC(ITS:ITE,J)=SFCFSW(1:(ITE-ITS+1))%upfx0
      RLWOUC(ITS:ITE,J)=SFCFLW(1:(ITE-ITS+1))%upfx0

      RLWOUCtoa(ITS:ITE,J)=TOPFLW(1:(ITE-ITS+1))%upfx0
      RSWOUCtoa(ITS:ITE,J)=TOPFSW(1:(ITE-ITS+1))%upfx0

      RLWTOA(ITS:ITE,J)=TOPFLW(1:(ITE-ITS+1))%upfxc*DT
      RSWTOA(ITS:ITE,J)=TOPFSW(1:(ITE-ITS+1))%upfxc*DT

!-- Modify ALBEDO array if using MODIS albedos

      IF (IALBflg == 1) THEN
         DO I=ITS,ITE          !-- RANGE
            IF (RSWIN(I,J) > 0.) THEN
               ALBEDO(I,J)=RSWOUT(I,J)/RSWIN(I,J)
            ELSE
               DUMMY=I-ITS+1   !-- from 1 to LENIVEC
               ALBEDO(I,J)=SFALB(DUMMY)
            ENDIF
         ENDDO
      ENDIF

!=================================================================
! For non GFDL type cloud (use cloud fields from outputs of GRRAD)
!=================================================================

      IF (ICWP /= -1) THEN
         IF ( LSSAV ) THEN
         !===========================================================
         ! Eliminate cloud fraction form GR1 & RH<95% (20140334, Lin)
         ! EPSQ2=1.e-8 (and not EPSQ=1.e-12) based on multiple tests
         !===========================================================



            DO I=1,(ITE-ITS+1)
              ARG_CW = MAXVAL( CW((its+I-1),J,1:LM) )  !- for *total condensate*
              ! ARG_CW = MAXVAL( GR1(I,1:LM,3) )  !- for *total condensate*

              IF (ARG_CW<EPSQ2) THEN  
                 DO L=1,LM
                    CLDCOV_V(I,L) = 0.d0
                 ENDDO
                 DO NC=1,5
                    CLDSA_V(I,NC) = 0.d0
                 ENDDO
              ENDIF
            ENDDO
         !===== end of eliminating extra cloud fraction =====
         !=========================================================

            DO L=1,LM
               CLDFRA(ITS:ITE,J,L)=CLDCOV_V(1:(ITE-ITS+1),L)
               CSMID(ITS:ITE,J,L)=CLDCOV_V(1:(ITE-ITS+1),L)
            ENDDO

            CFRACL(ITS:ITE,J)=CLDSA_V(1:(ITE-ITS+1),1)
            CFRACM(ITS:ITE,J)=CLDSA_V(1:(ITE-ITS+1),2)
            CFRACH(ITS:ITE,J)=CLDSA_V(1:(ITE-ITS+1),3)
!
!@@@ To Do:  @@@
!@@@ Add CFRACT array to calculate the total cloud fraction, replace
!    the instantaneous cloud fraction in the post
!
            ACFRST(ITS:ITE,J)=ACFRST(ITS:ITE,J) + CLDSA_V(1:(ITE-ITS+1),4)
!-- Added a time-averaged convective cloud fraction calculation
            WHERE (CU_cloud(1:(ITE-ITS+1))) ACFRCV(ITS:ITE,J)=ACFRCV(ITS:ITE,J)+CLDSA_V(1:(ITE-ITS+1),4)
         ELSE
            PRINT *, '*** CLDFRA=0, need to set LSSAV=TRUE'
            STOP
         ENDIF
      ENDIF

!-----------------------------------------------------------------------
!
      ENDDO     ! --- END J LOOP for grrad


!
!-----------------------------------------------------------------------
!***  LONGWAVE
!-----------------------------------------------------------------------
!
      IF(MOD(NTIMESTEP,NRADL)==0)THEN
        DO J=JTS,JTE
          DO I=ITS,ITE
!
            TDUM=T(I,J,LM)
            SIGT4(I,J)=STBOLT*TDUM*TDUM*TDUM*TDUM
!
          ENDDO
        ENDDO
      ENDIF 
!-----------------------------------------------------------------------


!
!*** --------------------------------------------------------------------------
!***  DETERMINE THE FRACTIONAL CLOUD COVERAGE FOR HIGH, MID
!***  AND LOW OF CLOUDS FROM THE CLOUD COVERAGE AT EACH LEVEL
!***
!***  NOTE: THIS IS FOR DIAGNOSTICS ONLY!!!
!***
!***
!
!----------------------------------------------------------------------------
      ICWP_Test2: IF (ICWP==-1) THEN   !-- *** Start of old NAM/GFDL cloud ***
!----------------------------------------------------------------------------

       DO J=JTS,JTE
       DO I=ITS,ITE
!!
       DO L=0,LM
         CLDAMT(L)=0.
       ENDDO
!!
!!***  NOW GOES LOW, MIDDLE, HIGH
!!
       DO 480 NLVL=1,3
       CLDMAX=0.
       MALVL=LM
       LLTOP=LM+1-LTOP(NLVL)   !!!!COMES FROM GFDL INIT
!!***
!!***  GO TO THE NEXT CLOUD LAYER IF THE TOP OF THE CLOUD-TYPE IN
!!***  QUESTION IS BELOW GROUND OR IS IN THE LOWEST LAYER ABOVE GROUND.
!!***
       IF(LLTOP.GE.LM)GO TO 480
!!
       IF(NLVL.GT.1)THEN
         LLBOT=LM+1-LTOP(NLVL-1)-1
         LLBOT=MIN(LLBOT,LM1)
       ELSE
         LLBOT=LM1
       ENDIF
!!
       DO 435 L=LLTOP,LLBOT
       CLDAMT(L)=AMAX1(CSMID(I,J,L),CCMID(I,J,L))
       IF(CLDAMT(L).GT.CLDMAX)THEN
         MALVL=L
         CLDMAX=CLDAMT(L)
       ENDIF
   435 CONTINUE
!!*********************************************************************
!! NOW, CALCULATE THE TOTAL CLOUD FRACTION IN THIS PRESSURE DOMAIN
!! USING THE METHOD DEVELOPED BY Y.H., K.A.C. AND A.K. (NOV., 1992).
!! IN THIS METHOD, IT IS ASSUMED THAT SEPERATED CLOUD LAYERS ARE
!! RADOMLY OVERLAPPED AND ADJACENT CLOUD LAYERS ARE MAXIMUM OVERLAPPED.
!! VERTICAL LOCATION OF EACH TYPE OF CLOUD IS DETERMINED BY THE THICKEST
!! CONTINUING CLOUD LAYERS IN THE DOMAIN.
!!*********************************************************************
       CL1=0.0
       CL2=0.0
       KBT1=LLBOT
       KBT2=LLBOT
       KTH1=0
       KTH2=0
!!
       DO 450 LL=LLTOP,LLBOT
       L=LLBOT-LL+LLTOP
       BIT1=.FALSE.
       CR1=CLDAMT(L)
       BITX=(P8W(I,L).GE.PTOPC(NLVL+1)).AND.                           &
      &     (P8W(I,L).LT.PTOPC(NLVL)).AND.                             &
      &     (CLDAMT(L).GT.0.0)
       BIT1=BIT1.OR.BITX
       IF(.NOT.BIT1)GO TO 450
!!***
!!***  BITY=T: FIRST CLOUD LAYER; BITZ=T:CONSECUTIVE CLOUD LAYER
!!***  NOTE:  WE ASSUME THAT THE THICKNESS OF EACH CLOUD LAYER IN THE
!!***         DOMAIN IS LESS THAN 200 MB TO AVOID TOO MUCH COOLING OR
!!***         HEATING. SO WE SET CTHK(NLVL)=200*E2. BUT THIS LIMIT MAY
!!***         WORK WELL FOR CONVECTIVE CLOUDS. MODIFICATION MAY BE
!!***         NEEDED IN THE FUTURE.
!!***
       BITY=BITX.AND.(KTH2.LE.0)
       BITZ=BITX.AND.(KTH2.GT.0)
!!
       IF(BITY)THEN
         KBT2=L
         KTH2=1
       ENDIF
!!
       IF(BITZ)THEN
         KTOP1=KBT2-KTH2+1
         DPCL=P_PHY(I,KBT2)-P_PHY(I,KTOP1)
         IF(DPCL.LT.CTHK(NLVL))THEN
           KTH2=KTH2+1
         ELSE
           KBT2=KBT2-1
         ENDIF
       ENDIF
       IF(BITX)CL2=AMAX1(CL2,CR1)
!!***
!!*** AT THE DOMAIN BOUNDARY OR SEPARATED CLD LAYERS, RANDOM OVERLAP.
!!*** CHOOSE THE THICKEST OR THE LARGEST FRACTION AMT AS THE CLD
!!*** LAYER IN THAT DOMAIN.
!!***
       BIT2=.FALSE.
       BITY=BITX.AND.(CLDAMT(L-1).LE.0.0.OR. &
            P8W(I,L-1).LT.PTOPC(NLVL+1))
       BITZ=BITY.AND.CL1.GT.0.0
       BITW=BITY.AND.CL1.LE.0.0
       BIT2=BIT2.OR.BITY
       IF(.NOT.BIT2)GO TO 450
!!
!!
       IF(BITZ)THEN
         KBT1=INT((CL1*KBT1+CL2*KBT2)/(CL1+CL2))
         KTH1=INT((CL1*KTH1+CL2*KTH2)/(CL1+CL2))+1
         CL1=CL1+CL2-CL1*CL2
       ENDIF
!!
       IF(BITW)THEN
         KBT1=KBT2
         KTH1=KTH2
         CL1=CL2
       ENDIF
!!
       IF(BITY)THEN
         KBT2=LLBOT
         KTH2=0
         CL2=0.0
       ENDIF
  450 CONTINUE
!
        CLDCFR(I,J,NLVL)=AMIN1(1.0,CL1)
        MTOP(I,J,NLVL)=MIN(KBT1,KBT1-KTH1+1)
        MBOT(I,J,NLVL)=KBT1

  480 CONTINUE

      ENDDO ! End DO I=ITS,ITE
      ENDDO ! End DO J=ITS,JTE

!!
      DO J=JTS,JTE
      DO I=ITS,ITE

        CFRACL(I,J)=CLDCFR(I,J,1)
        CFRACM(I,J)=CLDCFR(I,J,2)
        CFRACH(I,J)=CLDCFR(I,J,3)

        IF(CNCLD)THEN
          CFSmax=0.   !-- Maximum cloud fraction (stratiform component)
          CFCmax=0.   !-- Maximum cloud fraction (convective component)
          DO L=1,LM
            CFSmax=MAX(CFSmax, CSMID(I,J,L) )
            CFCmax=MAX(CFCmax, CCMID(I,J,L) )
          ENDDO
          ACFRST(I,J)=ACFRST(I,J)+CFSmax
          ACFRCV(I,J)=ACFRCV(I,J)+CFCmax
        ELSE
  !--- Count only locations with grid-scale cloudiness, ignore convective clouds
  !    (option not used, but if so set to the total cloud fraction)
          CFRAVG=1.-(1.-CFRACL(I,J))*(1.-CFRACM(I,J))*(1.-CFRACH(I,J))
          ACFRST(I,J)=ACFRST(I,J)+CFRAVG
        ENDIF

      ENDDO  !  DO I=ITS,ITE
      ENDDO  !  DO J=JTS,JTE

!-----------------------------------------------------------------------
      ENDIF  ICWP_Test2   !*** End of Old NAM/GFDL cloud ***
!-----------------------------------------------------------------------


      END SUBROUTINE RRTM


!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE RRTM_INIT(SHALF,PPTOP,LM)
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: LM
      REAL,INTENT(IN) :: PPTOP
      REAL,DIMENSION(LM),INTENT(IN) :: SHALF
!
      INTEGER :: I,J,N
      REAL :: PCLD,XSD,SQR2PI
      REAL :: SSLP=1013.25
      REAL, PARAMETER :: PTOP_HI=150.,PTOP_MID=350.,PTOP_LO=642.,       &
     &                   PLBTM=105000.
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!***  INITIALIZE DIAGNOSTIC LOW,MIDDLE,HIGH CLOUD LAYER PRESSURE LIMITS.
!
      LTOP(1)=0
      LTOP(2)=0
      LTOP(3)=0
!
      DO N=1,LM
        PCLD=(SSLP-PPTOP*10.)*SHALF(N)+PPTOP*10.
        IF(PCLD>=PTOP_LO)LTOP(1)=N
        IF(PCLD>=PTOP_MID)LTOP(2)=N
        IF(PCLD>=PTOP_HI)LTOP(3)=N
      ENDDO
!***
!***  ASSIGN THE PRESSURES FOR CLOUD DOMAIN BOUNDARIES
!***
      PTOPC(1)=PLBTM
      PTOPC(2)=PTOP_LO*100.
      PTOPC(3)=PTOP_MID*100.
      PTOPC(4)=PTOP_HI*100.
!
!---  Calculate the area under the Gaussian curve at the start of the
!---  model run and build the look up table AXSD
!
      SQR2PI=SQRT(2.*PI)
      RSQR=1./SQR2PI
      DO I=1,NXSD
        XSD=REAL(I)*DXSD
        AXSD(I)=GAUSIN(XSD)
      ENDDO
!
!-----------------------------------------------------------------------
      END SUBROUTINE RRTM_INIT
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!-- BSF 20120319: Change w/r/t original rsipath2 is that this code uses
!      simpler approximations for determining water paths and effective
!      radius for cloud ice and snow.  The original rsipath2 code
!      requires 'subroutine GSMCONST' within module_bfmicrophysics.f,
!      which is not called within the NMMB.
!

      subroutine rsipath2_tmp                                           &

!  ---  inputs:
     &     ( plyr, plvl, tlyr, qlyr, qcwat, qcice, qrain, rrime,        &
     &       IM, LEVS, iflip, flgmin,                                   &
!  ---  outputs:
     &       cwatp, cicep, rainp, snowp, recwat, rerain, resnow, snden  &
     &     )

! =================   subprogram documentation block   ================ !
!                                                                       !
! abstract:  this program is a modified version of ferrier's original   !
!   "rsipath" subprogram.  it computes layer's cloud liquid, ice, rain, !
!   and snow water condensate path and the partical effective radius    !
!   for liquid droplet, rain drop, and snow flake.                      !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IM,LEVS) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IM,LEVS+1):model level pressure in mb (100Pa)                !
!   tlyr  (IM,LEVS) : model layer mean temperature in k                 !
!   qlyr  (IM,LEVS) : layer specific humidity in gm/gm                  !
!   qcwat (IM,LEVS) : layer cloud liquid water condensate amount        !
!   qcice (IM,LEVS) : layer cloud ice water condensate amount           !
!   qrain (IM,LEVS) : layer rain drop water amount                      !
!   rrime (IM,LEVS) : mass ratio of total to unrimed ice ( >= 1 )       !
!   IM              : horizontal dimention                              !
!   LEVS            : vertical layer dimensions                         !
!   iflip           : control flag for in/out vertical indexing         !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   flgmin          : Minimum large ice fraction                        !
!   lprnt           : logical check print control flag                  !
!                                                                       !
! output variables:                                                     !
!   cwatp (IM,LEVS) : layer cloud liquid water path          (g/m**2)   !
!   cicep (IM,LEVS) : layer cloud ice water path             (g/m**2)   !
!   rainp (IM,LEVS) : layer rain water path                  (g/m**2)   !
!   snowp (IM,LEVS) : layer snow water path                  (g/m**2)   !
!   recwat(IM,LEVS) : layer cloud eff radius for liqid water (micron)   !
!   rerain(IM,LEVS) : layer rain water effective radius      (micron)   !
!   resnow(IM,LEVS) : layer snow flake effective radius      (micron)   !
!   snden (IM,LEVS) : 1/snow density                                    !
!                                                                       !
!                                                                       !
! usage:     call rsipath2                                              !
!                                                                       !
! subroutines called:  none                                             !
!                                                                       !
! program history log:                                                  !
!      xx-xx-2001   b. ferrier     - original program                   !
!      xx-xx-2004   s. moorthi     - modified for use in gfs model      !
!      05-20-2004   y. hou         - modified, added vertical index flag!
!                     to reduce data flipping, and rearrange code to    !
!                     be comformable with radiation part programs.      !
!      02-24-2012   b. ferrier     - simple, temporary fix              !
!                     to separate cloud ice & snow                      !
!                                                                       !
!  ====================    end of description    =====================  !
!

      implicit none

!  ---  constant parameter:
      real, parameter :: CEXP= 1./3., EPSQ=1.E-12

! CN0r0=2511.54=1.E6/(3.1415*1000.*8.e6)**.25
      real, parameter :: CN0r0=2511.54    !-- N0r=8.e6 m^-4, RHOL=1000 kg m^-3

!
!-- Cloud droplet distribution:
!     Reff=0.5*Deff, (1)
!       Deff=Dmean (monodisperse) or 3*Dmean (exponential), (2)
!       Reff=0.5*Dmean (monodisperse) or 1.5*Dmean (exponential), (3)
!
!     Assume monodisperse below ....
!       Dmean=mean diameter=1.e6*(6*rho*Qcw/(pi*TNW*RHOL))**(1/3) in microns (4)
!     Combining (1)-(4),
!       Reff=0.5e6*(6*rho*Qcw/(pi*Ncw*RHOL))**(1/3) in microns, (5)
!       where rho*Qcw in kg/m**3, Ncw in m^-3
!     RHOL=1000 kg m^-3, Ncw=1.e6*TNW, TNW in cm^-3 (6)
!     Substitute (6) into (5),
!       Reff=620.35*(Qcw/TNW)**(1/3) for rho*Qcw in kg/m**3, TNW in cm^-3 (monodisperse) (7)

      real, parameter :: TNW=200.                    !--  Droplet # cm^-3

!  ---  inputs:
      real(kind_phys), dimension(:,:), intent(in) ::                    &
     &       plyr, plvl, tlyr, qlyr, qcwat, qcice, qrain, rrime

      integer, intent(in) :: IM, LEVS, iflip
      real(kind_phys), dimension(:),   intent(in) :: flgmin

!  ---  output:
      real(kind_phys), dimension(:,:), intent(out) ::                   &
     &       cwatp, cicep, rainp, snowp, recwat, rerain, resnow, snden

!  ---  locals:

      real(kind_phys) :: dsnow, qsnow, qclice, fsmall, xsimass, pfac,   &
     &                   nlice, xli, nlimax, dum, tem,                  &
     &                   rho, cpath, rc, totcnd, tc, recw1

      integer :: i, k, indexs, ksfc, k1
!
!===>  ...  begin here
!
      recw1=620.35/TNW**CEXP    !-- Monodisperse

      !--- hydrometeor's optical path

      do k = 1, LEVS
        do i = 1, IM
           cwatp(i,k) = 0.0
           cicep(i,k) = 0.0
           rainp(i,k) = 0.0
           snowp(i,k) = 0.0
           snden(i,k) = 0.0
           resnow(i,k) = resnow_def
        enddo
      enddo

!  ---  set up pressure related arrays, convert unit from mb to cb (10Pa)
!       cause the rest part uses cb in computation

      if (iflip == 0) then        ! data from toa to sfc
        ksfc = levs + 1
        k1   = 0
      else                        ! data from sfc to top
        ksfc = 1
        k1   = 1
      endif                       ! end_if_iflip
!
      do k = 1, LEVS
        do i = 1, IM
          totcnd = qcwat(i,k) + qcice(i,k) + qrain(i,k)
          qsnow = 0.0
          if(totcnd > EPSQ) then

!  ---  air density (rho, kg/m**3), temperature (tc, deg C)
!       model mass thickness (cpath, g/m**2)

            rho   = 100. * plyr(i,k)                                     &
     &            / (con_rd* tlyr(i,k) * (1.0 + con_fvirt*qlyr(i,k)))

           ! ---- convert presure unit from mb to Pa ==> x100
           ! ---- convert NT unit from kg/sec^2 to g/sec^2 ==> x1000
           ! ---- combine the conversion ==> "100000.0"

            cpath = abs(plvl(i,k+1) - plvl(i,k)) * (100000.0 / con_g)
            tc    = tlyr(i,k) - con_t0c

!! cloud water
!
!  ---  effective radius (recwat) & total water path (cwatp):
!       assume monodisperse distribution of droplets (see above)

            if (qcwat(i,k) > 0.0) then
              tem         = recw1*(rho*qcwat(i,k))**CEXP
              recwat(i,k) = MAX(recwat_def, tem)
              cwatp (i,k) = cpath * qcwat(i,k)           ! cloud water path
            endif

!! rain
!
!  ---  effective radius (rerain) & total water path (rainp):
!       factor of 1.5 accounts for r**3/r**2 moments for exponentially
!       distributed drops in effective radius calculations
!       (from m.d. chou's code provided to y.-t. hou)

            if (qrain(i,k) > 0.0) then
              tem         = CN0r0 * sqrt(sqrt(rho*qrain(i,k)))
              rerain(i,k) = 1.5*max(50., min(1000.,tem))
              rainp (i,k) = cpath * qrain(i,k)           ! rain water path
            endif

!-- Treat all ice as "cloud ice" and no snow (Oct 2013)

            cicep (i,k) = cpath * qcice(i,k)      ! cloud ice path

          endif                                   ! end if_totcnd block

        enddo
      enddo
!
!...................................
      end subroutine rsipath2_tmp
!-----------------------------------

!----------------------------------------------------------------------
!
      REAL FUNCTION GAUSIN(xsd)
      REAL, PARAMETER :: crit=1.e-3
      REAL A1,A2,RN,B1,B2,B3,SUM,xsd
!
!  This function calculate area under the Gaussian curve between mean
!  and xsd # of standard deviation (03/22/2004  Hsin-mu Lin)
!
      a1=xsd*RSQR
      a2=exp(-0.5*xsd**2)
      rn=1.
      b1=1.
      b2=1.
      b3=1.
      sum=1.
      do while (b2 .gt. crit)
         rn=rn+1.
         b2=xsd**2/(2.*rn-1.)
         b3=b1*b2
         sum=sum+b3
         b1=b3
      enddo
      GAUSIN=a1*a2*sum
      RETURN
      END FUNCTION GAUSIN
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      END MODULE MODULE_RA_RRTM
!
!-----------------------------------------------------------------------
