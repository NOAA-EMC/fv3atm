!-----------------------------------------------------------------------
!
      MODULE MODULE_TURBULENCE
!
!-----------------------------------------------------------------------
!
!***  THE OUTER DRIVER FOR THE SFC LAYER, LSM, AND FULL 3-D TURBULENCE
!***  PLUS THE WRF TURBULENCE DRIVER AND THE VARIOUS TURBULENCE SCHEMES.
!
!-----------------------------------------------------------------------
! HISTORY LOG:
!
!   2008-07-28  Vasic - Turned off counters (now computed in
!                       SET_INTERNAL_STATE_PHY).
!   2009-10-26  Jovic - Remove WRF driver from TURBL
!   2010-09-10  Weiguo Wang - add GFS PBL option
!   2010-10-06  Weiguo Wang - add RSWTT, RLWTT, used by GFS PBL
!   2014-06-24  Weiguo Wang - add GFDL surface layer
!   2014-06-25  Weiguo Wang - add GFSPBL for Hurricane option
!   2016-08-29  Weiguo Wang - add GFSPBLEDMF version for Hurrican applications
!-----------------------------------------------------------------------
!
      USE MODULE_KINDS
!
      USE MODULE_LS_NOAHLSM
      USE MODULE_LS_LISS
      USE MODULE_GWD
!
      USE MODULE_CONSTANTS,ONLY : CP,ELIV,ELWV,EPSQ                     &
                                 ,G,P608,PI,PQ0,R_D,R_V,RHOWATER        &
                                 ,STBOLT,CAPPA   &
                                 ,EP_1,EP_2
!
      USE MODULE_SF_JSFC,ONLY : JSFC
      USE MODULE_SF_GFDL,ONLY : SF_GFDL
      USE MODULE_BL_MYJPBL,ONLY : MYJPBL
      USE MODULE_BL_GFSPBL,ONLY : GFSPBL
      USE MODULE_BL_GFSPBLHUR,ONLY : GFSPBLHUR
      USE MODULE_BL_GFSPBLEDMFHUR,ONLY : GFSPBLEDMFHUR
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: TURBL
      PUBLIC :: LSMSCHEME, LISSSCHEME, GFDLSLABSCHEME
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE TURBULENCE OPTIONS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT),PARAMETER :: MYJPBLSCHEME=2
      INTEGER(kind=KINT),PARAMETER :: GFSPBLSCHEME=9
      INTEGER(kind=KINT),PARAMETER :: GFSPBLHURSCHEME=93
      INTEGER(kind=KINT),PARAMETER :: GFSPBLEDMFHURSCHEME=94
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE SURFACE LAYER OPTIONS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT),PARAMETER :: JSFCSCHEME=2   &
                                     ,GFDLSFCSCHEME=88
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE LANDSURFACE OPTIONS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT),PARAMETER :: LSMSCHEME   =2                    &
                                     ,LISSSCHEME  =101                  &
                                     ,GFDLSLABSCHEME=0
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE TURBL(NTSD,DT,NPHS                                     &
                      ,U_PHY,V_PHY                                      &
                      ,PHINT,PHMID                                      &
                      ,SM,CZEN,CZMEAN,SIGT4,RLWIN,RSWIN,RADOT           &
                      ,RLWTT,RSWTT                                      &
                      ,T,Q,CWM,F_ICE,F_RAIN,F_RIMEF,SR                  &
                      ,Q2,DTDT,DUDT,DVDT                                &
                      ,THS,TSFC,SST,PREC,SNO,SNOWC                      &
                      ,QC,QR,QI,QS,QG                                   &
                      ,F_QC,F_QR,F_QI,F_QS,F_QG                         &
                      ,FIS,Z0,Z0BASE,USTAR,PBLH,XLEN_MIX,RMOL           &
                      ,AKHS,AKMS,AKHS_OUT,AKMS_OUT                      &
                      ,THZ0,QZ0,UZ0,VZ0,QSH                             &
                      ,STC,SMC,CMC,SMSTAV,SMSTOT,SSROFF,BGROFF          &
                      ,IVGTYP,ISLTYP,VEGFRC,GRNFLX                      &
                      ,SFCEXC,ACSNOW,ACSNOM,SNOPCX,SICE,TG,SOILTB       &
                      ,ALBASE,MXSNAL,ALBEDO,SH2O,SI,EPSR                &
                      ,U10,V10,TH10,Q10,TSHLTR,QSHLTR,PSHLTR,PSFC_OUT   &
                      ,T2,TWBS,QWBS,SFCSHX,SFCLHX,SFCEVP,EPI            &
                      ,POTEVP,POTFLX,SUBSHX,EVBSA,EVCWA                 &
                      ,APHTIM,ARDSW,ARDLW,ASRFC                         &
                      ,CROT,SROT,MIXHT                                  &
                      ,HSTDV,HCNVX,HASYW,HASYS,HASYSW,HASYNW,HLENW      &
                      ,HLENS,HLENSW,HLENNW,HANGL,HANIS,HSLOP,HZMAX      &
                      ,CDMB,CLEFF,SIGFAC,FACTOP,RLOLEV,DPMIN            &
                      ,RSWOUT,RSWTOA,RLWTOA                             &
                      ,ASWIN,ASWOUT,ASWTOA                              &
                      ,ALWIN,ALWOUT                                     &
                      ,GWDFLG                                           &
                      ,PCPFLG,DDATA                                     & ! PRECIP ASSIM
                      ,UCMCALL,IVEGSRC                                  &
                      ,TURBULENCE,SFC_LAYER                             &
                      ,LAND_SURFACE                                     &
                      ,MICROPHYSICS                                     &
                      ,LISS_RESTART                                     &
                      ,VAR_RIC,COEF_RIC_L,COEF_RIC_S,DISHEAT            &
                      ,ALPHA,SFENTH                                     &
                      ,IMS,IME,LM)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    TURBL       TURBULENCE OUTER DRIVER
!   PRGRMMR: BLACK           ORG: W/NP22     DATE: 02-04-19
!
! ABSTRACT:
!     TURBL DRIVES THE TURBULENCE SCHEMES
!
! PROGRAM HISTORY LOG (with changes to called routines) :
!   95-03-15  JANJIC     - ORIGINATOR OF THE SUBROUTINES CALLED
!   BLACK & JANJIC       - ORIGINATORS OF THE DRIVER
!   95-03-28  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!   96-03-29  BLACK      - ADDED EXTERNAL EDGE; REMOVED SCRCH COMMON
!   96-07-19  MESINGER   - ADDED Z0 EFFECTIVE
!   98-??-??  TUCCILLO   - MODIFIED FOR CLASS VIII PARALLELISM
!   98-10-27  BLACK      - PARALLEL CHANGES INTO MOST RECENT CODE
!   02-01-10  JANJIC     - MOIST TURBULENCE (DRIVER, MIXLEN, VDIFH)
!   02-01-10  JANJIC     - VERT. DIF OF Q2 INCREASED (Grenier & Bretherton)
!   02-02-02  JANJIC     - NEW SFCDIF
!   02-04-19  BLACK      - ORIGINATOR OF THIS OUTER DRIVER FOR WRF
!   02-05-03  JANJIC     - REMOVAL OF SUPERSATURATION AT 2m AND 10m
!   04-11-18  BLACK      - THREADED
!   06-10-25  BLACK      - BUILT INTO NMMB PHYSICS COMPONENT
!   08-07-28  VASIC      - Turned off counters (now computed in
!                            SET_INTERNAL_STATE_PHY).
!   08-08     JANJIC     - Synchronize WATER array and Q.
!
!
! USAGE: CALL TURBL FROM PHY_RUN
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM
!$$$
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT),PARAMETER  :: JMS=1,JME=1
!
      INTEGER(kind=KINT),INTENT(IN) :: IMS,IME,LM                       &
                                      ,NPHS,NTSD,UCMCALL,IVEGSRC
!
      INTEGER(kind=KINT),DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: ISLTYP &
                                                                 ,IVGTYP
!
      REAL(kind=KFPT),INTENT(IN) :: DT,CDMB,CLEFF,SIGFAC,FACTOP,RLOLEV,DPMIN
!
      REAL(kind=KFPT),INTENT(INOUT) :: APHTIM,ARDSW,ARDLW,ASRFC
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: ALBASE &
                                                                 ,MXSNAL
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: CROT,SROT &
                           ,HSTDV,HCNVX,HASYW,HASYS,HASYSW,HASYNW,HLENW  &
                           ,HLENS,HLENSW,HLENNW,HANGL,HANIS,HSLOP,HZMAX
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: CZEN     &
                                                              ,CZMEAN   &
                                                              ,FIS      &
                                                              ,RLWIN    &
                                                              ,RLWTOA   &
                                                              ,RSWIN    &
                                                              ,RSWOUT   &
                                                              ,RSWTOA   &
                                                              ,SICE     &
                                                              ,SIGT4    &
                                                              ,SST,TG   &
                                                              ,VEGFRC
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: EPSR  &
                                                                 ,SM,SR
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: GRNFLX,QWBS,RADOT  &
                                                               ,SFCEXC,SMSTAV      &
                                                               ,SOILTB,TWBS
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: ACSNOM,ACSNOW    &
                                                                 ,AKHS,AKMS        &
                                                                 ,ALBEDO           &
                                                                 ,ALWIN,ALWOUT     &
                                                                 ,ASWIN            &
                                                                 ,ASWOUT,ASWTOA    &
                                                                 ,BGROFF,CMC       &
                                                                 ,EVBSA,EVCWA      &
                                                                 ,PBLH,POTEVP      &
                                                                 ,POTFLX,PREC      &
                                                                 ,QSH,QZ0,RMOL     &
                                                                 ,SFCEVP,EPI       &
                                                                 ,SFCLHX,SFCSHX    &
                                                                 ,SI,SMSTOT        &
                                                                 ,SNO,SNOPCX       &
                                                                 ,SNOWC            &
                                                                 ,SSROFF,SUBSHX    &
                                                                 ,T2,THS,THZ0      &
                                                                 ,TSFC             &
                                                                 ,USTAR,UZ0,VZ0    &
                                                                 ,Z0,Z0BASE
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: AKHS_OUT,AKMS_OUT  &
                                                               ,MIXHT,PSHLTR       &
                                                               ,Q10,QSHLTR         &
                                                               ,TH10,TSHLTR        &
                                                               ,U10,V10
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,1:LM)  ,INTENT(IN) :: PHMID
      REAL(kind=KFPT),DIMENSION(IMS:IME,1:LM+1),INTENT(IN) :: PHINT
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: CWM    &
                                                                      ,U_PHY,V_PHY &
                                                                      ,Q,Q2   &
                                                                      ,T,QC,QI,QR,QS,QG
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) ::  F_ICE  &
                                                                       ,F_RAIN
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(OUT) :: DTDT,DUDT,DVDT &
                                                                    ,XLEN_MIX
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,NUM_SOIL_LAYERS),INTENT(INOUT) :: SH2O,SMC,STC
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) ::  RSWTT, RLWTT &
                                                                    ,F_RIMEF
!
      LOGICAL(kind=KLOG),INTENT(INOUT) :: LISS_RESTART
!
      CHARACTER(99),INTENT(IN) :: LAND_SURFACE,MICROPHYSICS             &
                                 ,SFC_LAYER,TURBULENCE
!
!  For precip assimilation:
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: DDATA
      LOGICAL(kind=KLOG),INTENT(IN) :: F_QC,F_QR,F_QI,F_QS,F_QG
      LOGICAL(kind=KLOG),INTENT(IN) :: GWDFLG,PCPFLG

! FOR Hurricane PBL/SFCLAY
      REAL(kind=KFPT),INTENT(IN) :: SFENTH, ALPHA,VAR_RIC, COEF_RIC_L  &
                                   ,COEF_RIC_S ! from namelist
      LOGICAL(kind=KLOG),INTENT(IN)  :: DISHEAT    ! should be from namelist


!
!---------------------
!***  Local Variables
!---------------------
!
      REAL(kind=KFPT),PARAMETER :: XLV=ELWV
!
      INTEGER(kind=KINT) :: I,IEND,IJ,ISTR,IW,J,K,KOUNT_ALL         &
                           ,LENGTH_ROW,N,NRDL,NRL,NWL,SST_UPDATE
!
      INTEGER(kind=KINT) :: PBL_PHYSICS,SFCLAY_PHYSICS,SURFACE_PHYSICS
!
      INTEGER(kind=KINT) :: NUM_ROAD_LAYERS                             &
                           ,NUM_ROOF_LAYERS                             &
                           ,NUM_WALL_LAYERS
!
      INTEGER(kind=KINT),DIMENSION(IMS:IME,JMS:JME) :: UTYPE_URB2D,LPBL
!
      REAL(kind=KFPT) :: ALTITUDE,DECLIN_URB,DTBL,DTMIN,DTPHS,DZHALF           &
                        ,FACTOR,FACTRL,G_INV                                   &
                        ,QL,QLOWX,QCW,QRain,QCI,QSnow,QGraup                   &
                        ,RDTPHS,ROG,RXNER,SNO_FACTR                            &
                        ,TL,TLMH,TSFC2,XLVRW
!
      REAL(kind=KFPT),DIMENSION(1:LM-1) :: EPSL
      REAL(kind=KFPT),DIMENSION(1:LM)   :: EPSQ2
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME) ::CHKLOWQ,CHS,CHS2,CPM,CQS2      &
                                                   ,CT,CWMLOW,EVBS,EVCW           &
                                                   ,ELFLX,EXNSFC,FACTRS,FLHC,FLQC &
                                                   ,GZ1OZ0,PSFC_OUT,PSIH,PSIM     &
                                                   ,Q2X,QGH,QLOW,RAIN,RAINBL      &
                                                   ,RLW_DN_SFC,RSW_NET_SFC        &
                                                   ,RSW_DN_SFC,RIMEF              &
                                                   ,SFCEVPX,SFCZ,SNOW,SNOWH       &
                                                   ,TH2X,THLOW,TLOW               &
                                                   ,VGFRCK,XLAND

!! added for hwrf, tentatively
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME) :: MZNT,WSPD,  &
                                                    TAUX,TAUY   &
                                                   ,RC2D
      INTEGER(kind=KINT),DIMENSION(IMS:IME,JMS:JME) ::  KPBL2D
!     REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM) :: EXCH_M
      INTEGER(kind=KINT) :: NTSFLG

!! hwrf

!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM) :: DELP                    &
                                                        ,DZ,EXNER,RR,TH
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM+1) :: Z
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM) :: RTHBLTEN       &
                                                        ,RQBLTEN        &
                                                        ,RQCBLTEN       &
                                                        ,RQIBLTEN
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM) :: DUDT_GWD,DVDT_GWD

      REAL(kind=KFPT),DIMENSION(1:NUM_SOIL_LAYERS) :: DZB,DZR,DZG
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,1:NUM_SOIL_LAYERS,JMS:JME) :: TBL_URB3D   &
                                                           ,TGL_URB3D             &
                                                           ,TRL_URB3D
!
      REAL(kind=KFPT), DIMENSION( IMS:IME, JMS:JME ) :: AKMS_URB2D      &
                                                       ,COSZ_URB2D      &
                                                       ,FRC_URB2D       &
                                                       ,G_URB2D         &
                                                       ,GZ1OZ0_URB2D    &
                                                       ,LH_URB2D        &
                                                       ,OMG_URB2D       &
                                                       ,PSIH_URB2D      &
                                                       ,PSIM_URB2D      &
                                                       ,Q2_URB2D        &
                                                       ,QC_URB2D        &
                                                       ,RN_URB2D        &
                                                       ,SH_URB2D        &
                                                       ,TB_URB2D        &
                                                       ,TC_URB2D        &
                                                       ,TG_URB2D        &
                                                       ,TH2_URB2D       &
                                                       ,TR_URB2D        &
                                                       ,TS_URB2D        &
                                                       ,U10_URB2D       &
                                                       ,UC_URB2D        &
                                                       ,UST_URB2D       &
                                                       ,V10_URB2D       &
                                                       ,XLAT_URB2D      &
                                                       ,XXXB_URB2D      &
                                                       ,XXXC_URB2D      &
                                                       ,XXXG_URB2D      &
                                                       ,XXXR_URB2D

      REAL(kind=KFPT), DIMENSION( IMS:IME, JMS:JME )  :: RIB   ! Bulk Richardson Number
!
      LOGICAL(kind=KLOG) :: FRPCPN,MYJ,WARM_RAIN,FER_MIC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  Translate the package options in the config file needed by
!***  the Turbulence to their analogs in the WRF Registry so that
!***  the WRF surface and PBL drivers remain untouched.
!-----------------------------------------------------------------------
!

!     write(0,*)'select PBL=',TURBULENCE
      SELECT CASE (TRIM(TURBULENCE))
        CASE ('myj')
          PBL_PHYSICS=MYJPBLSCHEME
        CASE ('gfs')
          PBL_PHYSICS=GFSPBLSCHEME
        CASE ('gfshur')
          PBL_PHYSICS=GFSPBLHURSCHEME
        CASE ('gfsedmfhur')
          PBL_PHYSICS=GFSPBLEDMFHURSCHEME
        CASE DEFAULT
          WRITE(0,*)' User selected TURBULENCE=',TRIM(TURBULENCE)
          WRITE(0,*)' Improper selection of Turbulence scheme in TURBL'
          STOP !rv or call ABORT
      END SELECT
!
      SELECT CASE (TRIM(SFC_LAYER))
        CASE ('myj')
          SFCLAY_PHYSICS=JSFCSCHEME
        CASE ('gfs')
          SFCLAY_PHYSICS=JSFCSCHEME
        CASE ('gfdl')
          SFCLAY_PHYSICS=GFDLSFCSCHEME
        CASE DEFAULT
          WRITE(0,*)' User selected SFC_LAYER=',TRIM(SFC_LAYER)
          WRITE(0,*)' Improper selection of Surface Layer scheme in TURBL'
          STOP !rv or call ABORT
      END SELECT
!
      SELECT CASE (TRIM(LAND_SURFACE))
        CASE ('noah')
          SURFACE_PHYSICS=LSMSCHEME
        case ('liss')
          surface_physics=LISSSCHEME
        case ('gfdlslab')
          surface_physics=GFDLSLABSCHEME
        CASE DEFAULT
          WRITE(0,*)' User selected LAND_SURFACE=',TRIM(LAND_SURFACE)
          WRITE(0,*)' Improper selection of Land Surface scheme in TURBL'
          STOP !rv or call ABORT
      END SELECT
!
      IF(TRIM(MICROPHYSICS)=='fer' .OR. TRIM(MICROPHYSICS)=='fer_hires')THEN
        FER_MIC=.TRUE.
      ELSE
        FER_MIC=.FALSE.
      ENDIF
!
!rv
      EPSQ2(:)=0.02
      EPSL (:)=0.1
!rv   EPSL=SQRT(EPSQ2(K-1)*0.5)
!
!.......................................................................
!$omp parallel do private(j,k,i)
!.......................................................................
      DO K=1,LM
      DO J=JMS,JME
      DO I=IMS,IME
        DTDT(I,J,K)=0.
        DUDT(I,J,K)=0.
        DVDT(I,J,K)=0.
      ENDDO
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
      DTPHS=NPHS*DT
      RDTPHS=1./DTPHS
      G_INV=1./G
      ROG=R_D*G_INV
      FACTOR=-XLV*RHOWATER/DTPHS
!
      SST_UPDATE=0
!
!$omp parallel do private(j,i)
      DO J=JMS,JME
      DO I=IMS,IME
        RMOL(I,J)=0.     !Reciprocal of Monin-Obukhov length
        SFCEVPX(I,J)=0.  !Dummy for accumulated latent energy, not flux
      ENDDO
      ENDDO
!$omp end parallel do
!
      IF(SURFACE_PHYSICS==99.OR.SURFACE_PHYSICS==LISSSCHEME)THEN
        SNO_FACTR=1.
!$omp parallel do private(j,i)
        DO J=JMS,JME
        DO I=IMS,IME
          SNOWC(I,J)=0.
          IF(SNO(I,J)>0.) SNOWC(I,J)=1.
        ENDDO
        ENDDO
!$omp end parallel do
      ELSE
        SNO_FACTR=0.001
      ENDIF
!
!$omp parallel do private(j,i)
      DO J=JMS,JME
      DO I=IMS,IME
        VGFRCK(I,J)=100.*VEGFRC(I,J)
        SNOW(I,J)=SNO(I,J)
        SNOWH(I,J)=SI(I,J)*SNO_FACTR
        XLAND(I,J)=SM(I,J)+1.
        T2(I,J)=TSFC(I,J)
      ENDDO
      ENDDO
!$omp end parallel do
!
      IF(NTSD==0)THEN
!$omp parallel do private(j,i)
        DO J=JMS,JME
        DO I=IMS,IME
          Z0BASE(I,J)=Z0(I,J)
          IF(SM(I,J)>0.5.AND.SICE(I,J)>0.5)THEN  !Bandaid
            SM(I,J)=0.
          ENDIF
        ENDDO
        ENDDO
!$omp end parallel do
      ENDIF
!
!$omp parallel do private(k,j,i)
      DO K=1,LM
      DO J=JMS,JME
      DO I=IMS,IME
        Q2(i,j,k)=MAX(Q2(i,j,k),EPSQ2(K))
        DZ(I,J,K)=0.
      ENDDO
      ENDDO
      ENDDO
!$omp end parallel do
!
!$omp parallel do private(k,j,i)
      DO K=1,LM+1
      DO J=JMS,JME
      DO I=IMS,IME
        Z(I,J,K)=0.
      ENDDO
      ENDDO
      ENDDO
!$omp end parallel do

!-----------------------------------------------------------------------
!***  Prepare needed arrays
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do                                                       &
!$omp private(j,i,tlmh,factrl,k,ql,tl                                   &
!$omp        ,rxner),SCHEDULE(dynamic)
!.......................................................................
      DO J=JMS,JME
      DO I=IMS,IME
!
        EXNSFC(I,J)=(1.E5/PHINT(I,LM+1))**CAPPA
        THS(I,J)=(SST(I,J)*EXNSFC(I,J))*SM(I,J)+THS(I,J)*(1.-SM(I,J))
        TSFC(I,J)=THS(I,J)/EXNSFC(I,J)
        SFCZ(I,J)=FIS(I,J)*G_INV
        IF(PCPFLG.AND.DDATA(I,J).LT.100.)THEN
          RAIN(I,J)=DDATA(I,J)*RHOWATER
        ELSE
          RAIN(I,J)=PREC(I,J)*RHOWATER
        ENDIF
        RAINBL(I,J)=0.
        TH2X(I,J)=T(I,J,LM)*(1.E5/PHMID(I,LM))**CAPPA
        Q2X(I,J)=Q(I,J,LM)
        RIMEF(I,J)=MAX(1., F_RIMEF(I,J,LM))
!
!-----------------------------------------------------------------------
!*** Long and shortwave flux at ground surface
!-----------------------------------------------------------------------
!
        IF(CZMEAN(I,J)>0.)THEN
          FACTRS(I,J)=CZEN(I,J)/CZMEAN(I,J)
        ELSE
          FACTRS(I,J)=0.
        ENDIF
!
        IF(SIGT4(I,J)>0.)THEN
          TLMH=T(I,J,LM)
          FACTRL=STBOLT*TLMH*TLMH*TLMH*TLMH/SIGT4(I,J)
        ELSE
          FACTRL=0.
        ENDIF
!
!- RLWIN/RSWIN - downward longwave/shortwave at the surface
!
        RLW_DN_SFC(I,J)=RLWIN(I,J)*FACTRL
        RSW_NET_SFC(I,J)=(RSWIN(I,J)-RSWOUT(I,J))*FACTRS(I,J)
!
!- Instantaneous downward solar for NMM_LSM
!
        RSW_DN_SFC(I,J)=RSWIN(I,J)*FACTRS(I,J)
!
!-----------------------------------------------------------------------
!***  Fill the arrays for calling the inner driver.
!-----------------------------------------------------------------------
!
        Z(I,J,LM+1)=SFCZ(I,J)
!
!-----------------------------------------------------------------------
!***  Fill vertical working arrays.
!-----------------------------------------------------------------------
!
        DO K=LM,1,-1
!
          QL=MAX(Q(I,J,K),EPSQ)
          TL=T(I,J,K)
!
          RR(I,J,K)=PHMID(I,K)/(R_D*TL)
          RXNER=(1.E5/PHMID(I,K))**CAPPA
          EXNER(I,J,K)=1./RXNER
          TH(I,J,K)=TL*RXNER
!
          RQCBLTEN(I,J,K)=0.
          RQIBLTEN(I,J,K)=0.
          RTHBLTEN(I,J,K)=0.
          RQBLTEN(I,J,K)=0.
!
          DZ(I,J,K)=T(I,J,K)*(P608*QL+1.)*R_D                           &
                    *(PHINT(I,K+1)-PHINT(I,K))                          &
                    /(PHMID(I,K)*G)
          Z(I,J,K)=Z(I,J,K+1)+DZ(I,J,K)
!
          DELP(I,J,K)=PHINT(I,K+1)-PHINT(I,K)
!
        ENDDO
      ENDDO
      ENDDO

!.......................................................................
!$omp end parallel do
!.......................................................................
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,qlowx)
!.......................................................................
      DO J=JMS,JME
      DO I=IMS,IME
        TWBS(I,J)=0.
        QWBS(I,J)=0.
        THLOW(I,J)=TH(I,J,LM)
        TLOW(I,J)=T(I,J,LM)
        QLOW(I,J)=MAX(Q(I,J,LM),EPSQ)
        QLOWX=QLOW(I,J)/(1.-QLOW(I,J))
        QLOW(I,J)=QLOWX/(1.+QLOWX)
        CWMLOW(I,J)=CWM(I,J,LM)
        PBLH(I,J)=MAX(PBLH(I,J),0.)
        PBLH(I,J)=MIN(PBLH(I,J),Z(I,J,LM))
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!
!***  Call surface layer and land surface physics.
!
!-----------------------------------------------------------------------
!
      DTMIN = 0.
      DTBL  = 0.
!
      DO J=JMS,JME
      DO I=IMS,IME
        QGH(I,J)  = 0.
        CHS(I,J)  = 0.
        CPM(I,J)  = 0.
        CHS2(I,J) = 0.
        RAINBL(I,J) = RAINBL(I,J) + RAIN(I,J)
        RAINBL(I,J) = MAX (RAINBL(I,J), 0.0)
      ENDDO
      ENDDO
!
!------------
! Update SST
!------------
!
      IF (SST_UPDATE == 1) THEN
        DO J=JMS,JME
        DO I=IMS,IME
          IF(XLAND(I,J)>1.5)TSFC(I,J)=SST(I,J)
        ENDDO
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
      sfc_and_sfclyr: IF (NTSD==1 .OR. MOD(NTSD,NPHS)==0) THEN
!-----------------------------------------------------------------------
!
        MYJ = .FALSE.
        FRPCPN = .FALSE.
        DTMIN=DT/60.
!-----------------------------------------------------------------------
! Surface schemes need PBL time step for updates and accumulations.
! Assume these schemes provide no tendencies.
!-----------------------------------------------------------------------
!
        DTBL=DT*NPHS
!
!--------------------
!*** Save old values
!--------------------
!
        DO J=JMS,JME
        DO I=IMS,IME
          PSFC_OUT(I,J)=PHINT(I,LM+1)
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  SFCLAY_PHYSICS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

          sfclay_select: SELECT CASE(SFCLAY_PHYSICS)

            CASE (JSFCSCHEME)

              MYJ =.TRUE.
              CALL JSFC(NTSD,EPSQ2,SFCZ,DZ,                             &
                        PHMID,PHINT,TH,T,                               &
                        Q,QC,                                           &
                        U_PHY,V_PHY,Q2,                                 &
                        TSFC,QSH,THZ0,QZ0,UZ0,VZ0,                      &
                        XLAND,                                          &
                        VEGFRC,SNOWC,                                   & !added 5/17/2013
                        USTAR,Z0,Z0BASE,PBLH,RMOL,                      &
                        AKHS,AKMS,                                      &
                        CHS,CHS2,CQS2,TWBS,QWBS,ELFLX,FLHC,FLQC,        &
                        QGH,CPM,CT,                                     &
                        U10,V10,T2,TH2X,TSHLTR,TH10,Q2X,QSHLTR,Q10,     &
                        PSHLTR,RIB,                                     &
                        IMS,IME,JMS,JME,LM)

            CASE (GFDLSFCSCHEME)

              MYJ =.FALSE.
! not sure why this if-block is needed, just copied from hwrf
              NTSFLG=0
!20140813         if (SFCLAY_PHYSICS == 88 ) NTSFLG=1
!20140813 IF we use GFDL sfc, we have to set NTSFLG=1, meaning update TSK
!20140813 IF NOT, then some other surface model may take care of updating TSK
              if( surface_physics == GFDLSLABSCHEME ) NTSFLG=1
!
              CALL SF_GFDL(U3D=U_PHY,V3D=V_PHY,T3D=T,QV3D=Q,              &
                      P3D=PHMID,CP=CP, ROVCP=R_D/CP,R=R_D,                &
                      XLV=XLV,PSFC=PSFC_OUT,CHS=CHS,CHS2=CHS2,            &
                      CQS2=CQS2, CPM=CPM,                                 &
                      DT=DTBL, SMOIS=SMC,num_soil_layers=NUM_SOIL_LAYERS, &
                      ISLTYP=ISLTYP,ZNT=z0,                               &
                      MZNT=MZNT,                                          &
                      UST=USTAR,PSIM=PSIM,PSIH=PSIH,                      &
                      XLAND=XLAND,HFX=TWBS,QFX=QWBS,                      &
                      TAUX=TAUX,TAUY=TAUY,LH=ELFLX,                       &
                      GSW=RSW_DN_SFC,GLW=RLW_DN_SFC,TSK=TSFC,             &
                      FLHC=FLHC,FLQC=FLQC,                                &
                      QGH=QGH,QSFC=QSH,U10=U10,V10=V10,                   &
                      GZ1OZ0=GZ1OZ0,WSPD=WSPD,BR=RIB,ISFFLX=1,            &
                      EP1=EP_1,EP2=EP_2,KARMAN=0.4,NTSFLG=NTSFLG,         &
                      SFENTH=SFENTH,                                      &
                      ims=ims,ime=ime, jms=jms,jme=jme, lm=lm)
              DO I=ims,ime
               DO J=jms,jme
              CHKLOWQ(I,J)=1.0
               enddo
              enddo
            CASE DEFAULT

              WRITE(0,*)'The sfclay option does not exist: SFCLAY_PHYSICS = ', SFCLAY_PHYSICS

          END SELECT sfclay_select

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  SURFACE_PHYSICS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

          sfc_select: SELECT CASE(SURFACE_PHYSICS)

            CASE (LISSSCHEME)

              CALL LISS(DZ,Q,PHINT,RR,                                  &
                        T,TH,TSFC,CHS,                                  &
                        TWBS,QWBS,QGH,RSW_DN_SFC,RLW_DN_SFC,ELFLX,RMOL, &
                        SMSTAV,SMSTOT,SSROFF,                           &
                        BGROFF,IVGTYP,ISLTYP,VGFRCK,SFCEVPX,POTEVP,     &
                        GRNFLX,SFCEXC,ACSNOW,ACSNOM,SNOPCX,             &
                        ALBASE,TG,XLAND,SICE,QZ0,                       &
                        TH2X,Q2X,SNOWC,CQS2,QSH,SOILTB,CHKLOWQ,RAINBL,  &
                        NUM_SOIL_LAYERS,DTBL,DZSOIL,NTSD,               &
                        SMC,STC,SNOW,CMC,CPM,CAPPA,SR,                  &
                        ALBEDO,MXSNAL,SH2O,SNOWH,                       &
                        LISS_RESTART,                                   &
                        IMS,IME,JMS,JME,LM)

            CASE (LSMSCHEME)

              FRPCPN=.TRUE.
              NRL=1
              NWL=1
              NRDL=1
              FRC_URB2D=0.0

              CALL NOAHLSM(DZ,Q,PHINT,T,TSFC,                           &
                           TWBS,QWBS,ELFLX,GRNFLX,QGH,                  &
                           RSW_NET_SFC,RSW_DN_SFC,RLW_DN_SFC,           &
                           EVBS,EVCW,SMSTAV,SMSTOT,                     &
                           SSROFF,BGROFF,IVGTYP,ISLTYP,VGFRCK,          &
                           ALBEDO,ALBASE,Z0,Z0BASE,TG,XLAND,SICE,EPSR,  &
                           SNOWC,QSH,RAINBL,                            &
                           NUM_SOIL_LAYERS,DTBL,DZSOIL,NTSD,            &
                           SMC,STC,SNOW,CMC,                            &
                           CHS, CHS2, CQS2, CPM,CAPPA,                  &
                           SR,RIMEF, CHKLOWQ,QZ0,                       &
                           MYJ,FRPCPN,                                  &
                           SH2O,SNOWH,                                  & !H
                           U_PHY,V_PHY,                                 & !I
                           MXSNAL,                                      & !I
                           ACSNOM,ACSNOW,                               & !O
                           SNOPCX,                                      & !O
                           POTEVP, RIB,                                 & !O Added Bulk Richardson No.
                           IMS,IME,JMS,JME,LM,                          &
                           UCMCALL,IVEGSRC,                             &
! Optional urban
                           TR_URB2D,TB_URB2D,TG_URB2D,TC_URB2D,         & !H urban
                           QC_URB2D,UC_URB2D,                           & !H urban
                           XXXR_URB2D,XXXB_URB2D,XXXG_URB2D,XXXC_URB2D, & !H urban
                           TRL_URB3D,TBL_URB3D,TGL_URB3D,               & !H urban
                           SH_URB2D,LH_URB2D,G_URB2D,RN_URB2D,TS_URB2D, & !H urban
                           PSIM_URB2D,PSIH_URB2D,U10_URB2D,V10_URB2D,   & !O urban
                           GZ1OZ0_URB2D, AKMS_URB2D,                    & !O urban
                           TH2_URB2D,Q2_URB2D,USTAR,                    & !O urban
                           DECLIN_URB,COSZ_URB2D,OMG_URB2D,             & !I urban
                           XLAT_URB2D,                                  & !I urban
                           NRL, NWL,                                    & !I urban
                           NRDL, DZR, DZB, DZG,                         & !I urban
                           FRC_URB2D, UTYPE_URB2D                       & ! urban
                            )

              DO J=JMS,JME
              DO I=IMS,IME
                SFCEVPX(I,J)= SFCEVPX(I,J) + QWBS(I,J)*DTBL
                SFCEXC(I,J)= CHS(I,J)
                SOILTB(I,J)= STC(I,J,NUM_SOIL_LAYERS) !  nmmlsm vrbl., here only for output
              ENDDO
              ENDDO

              CALL SFCDIAGS(TWBS,QWBS,TSFC,QSH,CHS2,CQS2,T2,TH2X,Q2X,   &
                            PSFC_OUT,CP,R_d,CAPPA,                      &
                            IMS,IME,JMS,JME,LM)

              urban: IF(UCMCALL==1) THEN
!
                IW = 1
                IF( IVEGSRC == 1 ) IW = 13
!
                DO J=JMS,JME
                DO I=IMS,IME
                  IF( IVGTYP(I,J) == IW .OR. IVGTYP(I,J) == 31 .OR.   &
                      IVGTYP(I,J) == 32 .OR. IVGTYP(I,J) == 33 ) THEN
!
                    T2(I,J)   = FRC_URB2D(I,J)*TH2_URB2D(I,J)         &
                               +(1-FRC_URB2D(I,J))*T2(I,J)
                    TH2X(I,J) = T2(I,J)*(1.E5/PSFC_OUT(I,J))**CAPPA
                    Q2X(I,J)  = FRC_URB2D(i,j)*Q2_URB2D(I,J)          &
                               +(1-FRC_URB2D(I,J))* Q2X(I,J)
                    U10(I,J)  = U10_URB2D(I,J)
                    V10(I,J)  = V10_URB2D(I,J)
                    PSIM(I,J) = PSIM_URB2D(I,J)
                    PSIH(I,J) = PSIH_URB2D(I,J)
                    GZ1OZ0(I,J) = GZ1OZ0_URB2D(I,J)
                    AKHS(I,J) = CHS(I,J)
                    AKMS(I,J) = AKMS_URB2D(I,J)
!
                  END IF
                ENDDO
                ENDDO
              ENDIF  urban

            CASE(GFDLSLABSCHEME)
!            write(0,*)'GFDL SLAB LSM model is selected, whcih is included in GFDLSURFACE module'

            CASE DEFAULT

              WRITE(0,*) 'The surface option not exist: SURFACE_PHYSICS = ', SURFACE_PHYSICS

          END SELECT sfc_select

!-----------------------------------------------------
! Reset RAINBL in mm (Accumulation between PBL calls)
!-----------------------------------------------------
          DO J=JMS,JME
          DO I=IMS,IME
            RAINBL(I,J) = 0.
          ENDDO
          ENDDO

!
!-----------------------------------------------------------------------
      ENDIF  sfc_and_sfclyr
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!***  CALL FREE ATMOSPHERE TURBULENCE
!
!-----------------------------------------------------------------------
!
!***  The surface exchange coefficients AKHS and AKMS are actually
!***  multiplied by half the depth of the lowest layer.  We must retain
!***  those values for the next timestep so use auxilliary arrays for
!***  the output.
!
      DO J=JMS,JME
      DO I=IMS,IME
        DZHALF=0.5*DZ(I,J,LM)
        AKHS_OUT(I,J)=AKHS(I,J)*DZHALF
        AKMS_OUT(I,J)=AKMS(I,J)*DZHALF
      ENDDO
      ENDDO

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  PBL_PHYSICS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      pbl_select: SELECT CASE(PBL_PHYSICS)

        CASE (MYJPBLSCHEME)

          IF (NTSD == 1 .OR. MOD(NTSD,NPHS) == 0) THEN

              CALL MYJPBL(DT=DT,NPHS=NPHS,EPSL=EPSL,EPSQ2=EPSQ2         &
                         ,HT=SFCZ,STDH=HSTDV,DZ=DZ                      &
                         ,PMID=PHMID,PINH=PHINT,TH=TH,T=T,EXNER=EXNER &
                         ,Q=Q                                           &
                         ,CWM=QC                                        &
                         ,U=U_PHY,V=V_PHY                               &
                         ,TSK=TSFC,QSFC=QSH,CHKLOWQ=CHKLOWQ,THZ0=THZ0   &
                         ,QZ0=QZ0,UZ0=UZ0,VZ0=VZ0                       &
                         ,XLAND=XLAND,SICE=SICE,SNOW=SNOW               &
                         ,Q2=Q2,USTAR=USTAR,Z0=Z0                       &
!                        ,Q2=Q2,EXCH_H=EXCH_H,USTAR=USTAR,Z0=Z0         &
                         ,EL_MYJ=XLEN_MIX,PBLH=PBLH,KPBL=LPBL,CT=CT     &
                         ,AKHS=AKHS,AKMS=AKMS,ELFLX=ELFLX,MIXHT=MIXHT   &
                         ,RUBLTEN=DUDT                                  &
                         ,RVBLTEN=DVDT                                  &
                         ,RTHBLTEN=RTHBLTEN                             &
                         ,RQBLTEN=RQBLTEN                               &
                         ,RQCBLTEN=RQCBLTEN                             &
                         ,IMS=IMS,IME=IME,JMS=JMS,JME=JME,LM=LM)

          END IF

        CASE (GFSPBLSCHEME)    !! Wang 09-10-2010 added GFS PBL driver

          IF (NTSD == 1 .OR. MOD(NTSD,NPHS) == 0) THEN

              CALL GFSPBL(DT=DT,NPHS=NPHS,DP=DELP,AIRDEN=RR              &
                         ,RIB=RIB                                        &
                         ,PHMID=PHMID,PHINT=PHINT,T=T,ZINT=Z             &
                         ,Q=Q,QC=QC,QI=QI                                &
                         ,F_QC=F_QC,F_QI=F_QI                            &
                         ,U=U_PHY,V=V_PHY                                &
                         ,USTAR=USTAR                                    &
                         ,SHEAT=TWBS, LHEAT=QWBS*XLV*CHKLOWQ             &
                    !     ,SHEAT=TWBS, LHEAT=QWBS*XLV*CHKLOWQ            &    !! After testing, TWBS is regular
                                                                              !surface heat flux (i.e., up is +)
                         ,XLAND=XLAND                                    &
                         ,AKHS=AKHS,AKMS=AKMS                            &
                         ,THZ0=THZ0,QZ0=QZ0                              &
                         ,QSFC=QSH                                       &
                         ,TSK=TSFC,SNOW=SNOW,SICE=SICE,CHKLOWQ=CHKLOWQ   &
                         ,FACTRS=FACTRS,RSWTT=RSWTT,RLWTT=RLWTT          &    !! radiative heating
                         ,PBLH=PBLH,PBLK=LPBL                            &
                         ,MIXHT=MIXHT                                    &
                         ,RUBLTEN=DUDT                                   &
                         ,RVBLTEN=DVDT                                   &
                         ,RTHBLTEN=RTHBLTEN                              &
                         ,RQBLTEN=RQBLTEN                                &
                         ,RQCBLTEN=RQCBLTEN                              &
                         ,RQIBLTEN=RQIBLTEN                              &
                         ,IMS=IMS,IME=IME,JMS=JMS,JME=JME,LM=LM)

          END IF

        CASE(GFSPBLHURSCHEME)

          IF (NTSD == 1 .OR. MOD(NTSD,NPHS) == 0) THEN
           CALL GFSPBLHUR(U3D=U_PHY,V3D=V_PHY,TH3D=TH,T3D=T              &
                         ,QV3D=Q,QC3D=QC,QI3D=QS*0.0,P3D=PHMID,PI3D=EXNER    &
                         ,RUBLTEN=DUDT                                   &
                         ,RVBLTEN=DVDT                                &
                         ,RTHBLTEN=RTHBLTEN  &
                         ,RQVBLTEN=RQBLTEN   &
                         ,RQCBLTEN=RQCBLTEN  &
                         ,RQIBLTEN=RQIBLTEN  &
                         ,CP=CP,G=G,ROVCP=CAPPA,R=R_D,ROVG=R_D/G   &
                         ,F_QC=F_QC,F_QI=F_QI &
                   ,dz8w=DZ,z=Z,PSFC=PSFC_OUT &
                   ,UST=USTAR,PBL=PBLH,PSIM=PSIM,PSIH=PSIH &
                   ,HFX=TWBS,QFX=QWBS,TSK=TSFC&
                   ,GZ1OZ0=GZ1OZ0,WSPD=WSPD,BR=RIB  &
                   ,DT=DT*NPHS,KPBL2D=KPBL2D,EP1=EP_1,KARMAN=0.4  &
                   ,DISHEAT=DISHEAT  &
                   ,ALPHA=ALPHA     &
                   ,VAR_RIC=VAR_RIC &
                   ,U10=U10,V10=V10,ZNT=z0,MZNT=MZNT,rc2d=rc2d &
!                  ,DKU3D=EXCH_M,DKT3D=EXCH_H &
                   ,coef_ric_l=coef_ric_l &
                   ,coef_ric_s=coef_ric_s &
                   ,xland=xland           &
                   ,ims=ims,ime=ime,jms=jms,jme=jme,lm=lm)
            ENDIF

        CASE(GFSPBLEDMFHURSCHEME)
          IF (NTSD == 1 .OR. MOD(NTSD,NPHS) == 0) THEN
           CALL GFSPBLEDMFHUR(U3D=U_PHY,V3D=V_PHY,TH3D=TH,T3D=T      &
                   ,QV3D=Q,QC3D=QC,QI3D=QS*0.0,P3D=PHMID,PI3D=EXNER  &
                   ,RUBLTEN=DUDT                                     &
                   ,RVBLTEN=DVDT                                     &
                   ,RTHBLTEN=RTHBLTEN                                &
                   ,RQVBLTEN=RQBLTEN                                 &
                   ,RQCBLTEN=RQCBLTEN                                &
                   ,RQIBLTEN=RQIBLTEN                                &
                   ,CP=CP,G=G,ROVCP=CAPPA,R=R_D,ROVG=R_D/G           &
                   ,F_QC=F_QC,F_QI=F_QI                              &
                   ,dz8w=DZ,z=Z,PSFC=PSFC_OUT                        &
                   ,UST=USTAR,PBL=PBLH,PSIM=PSIM,PSIH=PSIH           &
                   ,HFX=TWBS,QFX=QWBS,TSK=TSFC                       &
                   ,GZ1OZ0=GZ1OZ0,WSPD=WSPD,BR=RIB                   &
                   ,DT=DT*NPHS,KPBL2D=KPBL2D,EP1=EP_1,KARMAN=0.4     &
                   ,DISHEAT=DISHEAT                                  &
                   ,ALPHA=ALPHA                                      &
                   ,VAR_RIC=VAR_RIC                                  &
                   ,U10=U10,V10=V10,ZNT=z0,MZNT=MZNT & !,rc2d=rc2d   &
!                  ,DKU3D=EXCH_M,DKT3D=EXCH_H                        &
                   ,coef_ric_l=coef_ric_l                            &
                   ,coef_ric_s=coef_ric_s                            &
                   ,xland=xland                                      &
                   ,rswtt=rswtt,rlwtt=rlwtt                          &
                   ,ims=ims,ime=ime,jms=jms,jme=jme,lm=lm)
            ENDIF

        CASE DEFAULT

          WRITE(0,*)'The pbl option does not exist: pbl_physics = ', pbl_physics

      END SELECT pbl_select

!-----------------------------------------------------------------------
!***  Note that the exchange coefficients for heat EXCH_H coming out of
!***  PBL_DRIVER are defined at the tops of the layers KTS to KTE-1
!***  if MODULE_BL_MYJPBL was invoked.
!-----------------------------------------------------------------------
!
      IF(TRIM(SFC_LAYER)/='myj')THEN
        DO J=JMS,JME
        DO I=IMS,IME
          IF(TSHLTR(I,J)<200..OR.TSHLTR(I,J)>350.)THEN
            WRITE(0,*)'Troublesome TSHLTR...I,J,TSHLTR,PSHLTR: ',       &
               I,J,TSHLTR(I,J),PSHLTR(I,J)
          ENDIF
        ENDDO
        ENDDO
      ENDIF
!
      IF(SURFACE_PHYSICS==99 .OR. SURFACE_PHYSICS==LISSSCHEME)THEN
        SNO_FACTR=1.
      ELSE
        SNO_FACTR=1000.
      ENDIF

!$omp parallel do private(j,i)
      DO J=JMS,JME
      DO I=IMS,IME
        SNO(I,J)=SNOW(I,J)
        SI(I,J)=SNOWH(I,J)*SNO_FACTR
      ENDDO
      ENDDO
!$omp end parallel do

!
!-----------------------------------------------------------------------
!***  Diagnostic radiation accumulation.
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do private(j,i,tsfc2)
!.......................................................................
      DO J=JMS,JME
      DO I=IMS,IME
!-- Remove the next 2 lines and uncomment "!was" lines below if not correct
        TSFC2=TSFC(I,J)*TSFC(I,J)
        RADOT(I,J)=EPSR(I,J)*STBOLT*TSFC2*TSFC2
!       RADOT(I,J)=-1.*EPSR(I,J)*STBOLT*TSFC2*TSFC2
!
        ASWIN (I,J)=ASWIN (I,J)+RSWIN(I,J)*FACTRS(I,J)*DT
        ASWOUT(I,J)=ASWOUT(I,J)+RSWOUT(I,J)*FACTRS(I,J)*DT
        ASWTOA(I,J)=ASWTOA(I,J)+RSWTOA(I,J)*FACTRS(I,J)
        ALWIN (I,J)=ALWIN (I,J)+RLW_DN_SFC(I,J)*DT
        ALWOUT(I,J)=ALWOUT(I,J)+RADOT (I,J)*DT
        THS(I,J)=TSFC(I,J)*EXNSFC(I,J)
        PREC(I,J)=0.
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!=======================================================================
!===  Begin gravity wave drag (GWD) and mountain blocking (MB)  ========
!=======================================================================
!
      IF(GWDFLG) THEN

        CALL GWD_DRIVER(DTPHS,U_PHY,V_PHY,T,Q                           &
                       ,Z,DELP                                          &
                       ,PHINT,PHMID,EXNER                               &
                       ,LPBL                                            &
                       ,HSTDV,HCNVX,HASYW,HASYS                         &
                       ,HASYSW,HASYNW,HLENW                             &
                       ,HLENS,HLENSW,HLENNW                             &
                       ,HANGL,HANIS,HSLOP,HZMAX                         &
                       ,CROT,SROT                                       &
                       ,CDMB,CLEFF,SIGFAC,FACTOP,RLOLEV,DPMIN           &
                       ,DUDT_GWD,DVDT_GWD                               &
                       ,IMS,IME,JMS,JME,LM)

!$omp parallel do private(k,j,i)
        DO K=1,LM
        DO J=JMS,JME
        DO I=IMS,IME
          DUDT(I,J,K)=DUDT(I,J,K)+DUDT_GWD(I,J,K)
          DVDT(I,J,K)=DVDT(I,J,K)+DVDT_GWD(I,J,K)
        ENDDO
        ENDDO
        ENDDO
!$omp end parallel do

      ENDIF
!
!=======================================================================
!=====  End gravity wave drag (GWD) and mountain blocking (MB)  ========
!=======================================================================
!
!-----------------------------------------------------------------------
!***  Update temperature, specific humidity, cloud, and q2.
!-----------------------------------------------------------------------
!

!.......................................................................
!$omp parallel do                                                       &
!$omp& private(j,k,i,QCW,QRain,QCI,QSnow,QGraup)
!.......................................................................
      DO K=1,LM
        DO J=JMS,JME
          DO I=IMS,IME
!rv  ------------  update is now in the SOLVER ------------
!rv         U_PHY(I,J,K)=U_PHY(I,J,K)+DUDT(I,J,K)*DTPHS
!rv         V_PHY(I,J,K)=V_PHY(I,J,K)+DVDT(I,J,K)*DTPHS
!rv         T(I,J,K)=T(I,J,K)+DTPHS*RTHBLTEN(I,J,K)*EXNER(I,J,K)
            DTDT(I,J,K)=RTHBLTEN(I,J,K)*EXNER(I,J,K)
            Q(I,J,K)=Q(I,J,K)+RQBLTEN(I,J,K)*DTPHS
            QC(I,J,K)=MAX(0.,QC(I,J,K)+RQCBLTEN(I,J,K)*DTPHS )
            QCW=QC(I,J,K)
            QRain=0.
            QCI=0.
            QSnow=0.
            QGraup=0.
            IF(F_QR) QRain=QR(I,J,K)
            IF(F_QI) THEN
              QCI=MAX(0.,QI(I,J,K)+RQIBLTEN(I,J,K)*DTPHS )
              QI(I,J,K)=QCI
            ENDIF
            IF(F_QS) QSnow=QS(I,J,K)
            IF(F_QG .AND. .NOT.FER_MIC) QGraup=QG(I,J,K)
!-- Couple CWM, F_ice, & F_rain arrays
            CWM(I,J,K)=QCW+QRain+QCI+QSnow+QGraup
            F_ICE(I,J,K)=0.
            F_RAIN(I,J,K)=0.
            IF(CWM(I,J,K)>EPSQ) F_ICE(I,J,K)=(QCI+QSnow+QGraup)/CWM(I,J,K)
            IF(QRain>EPSQ) F_RAIN(I,J,K)=QRain/(QCW+QRain)
          ENDDO
        ENDDO

      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***
!***  Save surface-related fields.
!***
!-----------------------------------------------------------------------
!
      XLVRW=DTPHS/(XLV*RHOWATER)
!.......................................................................
!$omp parallel do private(j,i)
!.......................................................................
      DO J=JMS,JME
      DO I=IMS,IME
!
!-----------------------------------------------------------------------
!***  Instantaneous sensible and latent heat fluX
!-----------------------------------------------------------------------
!
        QWBS(I,J)=  QWBS(I,J)*XLV*CHKLOWQ(I,J)
!
!-----------------------------------------------------------------------
!***  Accumulated quantities.
!***  In opnl LSM, SFCEVP appears to be in units of
!***  meters of liquid water.  It is coming from
!***  WRF module as kg/m**2.
!-----------------------------------------------------------------------
!
        SFCSHX(I,J)=SFCSHX(I,J)+TWBS(I,J)*DT
        SFCLHX(I,J)=SFCLHX(I,J)+QWBS(I,J)*DT
        EPI(I,J)   =-QWBS(I,J)*XLVRW
        SFCEVP(I,J)=SFCEVP(I,J)+EPI(I,J)
        POTEVP(I,J)=POTEVP(I,J)-QWBS(I,J)*SM(I,J)*XLVRW
        POTFLX(I,J)=POTEVP(I,J)*FACTOR
        SUBSHX(I,J)=SUBSHX(I,J)+GRNFLX(I,J)*DT
        EVBSA(I,J) =EVBSA(I,J)+EVBS(I,J)*DT
        EVCWA(I,J) =EVCWA(I,J)+EVCW(I,J)*DT
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  COUNTERS
!-----------------------------------------------------------------------
!
        APHTIM = APHTIM + 1.
        ARDSW  = ARDSW  + 1.
        ARDLW  = ARDLW  + 1.
        ASRFC  = ASRFC  + 1.

      END SUBROUTINE TURBL
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_TURBULENCE
!
!-----------------------------------------------------------------------
