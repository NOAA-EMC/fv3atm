!-----------------------------------------------------------------------
!
      MODULE MODULE_RADIATION
!
!-----------------------------------------------------------------------
!
!***  THE RADIATION DRIVERS AND PACKAGES
!
!---------------------
!--- Modifications ---
!---------------------
! 2010-04-02 Vasic - Removed WFR driver
!-----------------------------------------------------------------------
!
      USE MODULE_KINDS
      USE MODULE_RA_GFDL,ONLY   : GFDL,CAL_MON_DAY,ZENITH
      USE MODULE_RA_RRTM,ONLY   : RRTM
      USE MODULE_CONSTANTS,ONLY : CAPPA,CP,EPSQ,G,P608,PI,R_D,STBOLT
!
      USE module_mp_thompson, ONLY : cal_cldfra3
      use module_radiation_driver_nmmb,  only : radupdate_nmmb
      use machine, only : kind_phys
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: RADIATION
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE RADIATION PACKAGE OPTIONS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!***  Shortwave
!
      INTEGER,PARAMETER  :: GFDLSWSCHEME=99                             &  !<--- (GFDL)
                           ,SWRADSCHEME=1                               &  !<--- (Dudhia, WRF)
                           ,GSFCSWSCHEME=2                              &  !<--- (Goddard, WRF)
                           ,RRTMSWSCHEME=3                                 !<--- (RRTM)
!
!***  Longwave
!
      INTEGER,PARAMETER  :: GFDLLWSCHEME=99                             &  !<--- (GFDL)
                           ,RRTMLWSCHEME=3                                 !<--- (RRTM)
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE RADIATION(ITIMESTEP,DT,JULDAY,JULYR,XTIME,JULIAN       &
     &                    ,IHRST,NPHS,GLAT,GLON                         &
     &                    ,NRADS,NRADL                                  &
     &                    ,PRSI,PRSL                                    &
     &                    ,T,Q                                          &
     &                    ,THS,ALBEDO                                   &
     &                    ,QC,QR,QI,QS,QG,NI                            &
     &                    ,F_QC,F_QR,F_QI,F_QS,F_QG,F_NI                &
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
     &                    ,CUPPT,SNOW                                   &
     &                    ,HTOP,HBOT                                    &
     &                    ,SHORTWAVE,LONGWAVE                           &
     &                    ,CLDFRACTION                                  &
     &                    ,DYH                                          &
     &                    ,JDAT                                         &
     &                    ,CW,O3                                        &
     &                    ,F_ICE,F_RAIN                                 &
     &                    ,F_RIMEF                                      &
     &                    ,SI,TSKIN                                     &
     &                    ,Z0,SICE                                      &
     &                    ,MXSNAL                                       &
     &                    ,STDH,VVL                                     &
     &                    ,ALBVB,ALBNB                                  &  ! vis+uv & near IR beam
     &                    ,ALBVD,ALBND                                  &  ! vis+uv & near IR diffuse
     &                    ,SNOWC                                        &
     &                    ,LM,IMS,IME,MYPE)
!-----------------------------------------------------------------------
!***  NOTE ***
! RLWIN  - downward longwave at the surface (=GLW)
! RSWIN  - downward shortwave at the surface (=XXX)
! RSWINC - CLEAR-SKY downward shortwave at the surface (=SWDOWN, new for AQ)
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    RADIATION   RADIATION OUTER DRIVER
!   PRGRMMR: BLACK           ORG: W/NP22     DATE: 2002-06-04
!
! ABSTRACT:
!     RADIATION SERVES AS THE INTERFACE BETWEEN THE NMMB PHYSICS COMPONENT
!     AND THE WRF RADIATION DRIVER.
!
! PROGRAM HISTORY LOG:
!   02-06-04  BLACK      - ORIGINATOR
!   02-09-09  WOLFE      - CONVERTING TO GLOBAL INDEXING
!   04-11-18  BLACK      - THREADED
!   06-07-20  BLACK      - INCORPORATED INTO NMMB PHYSICS COMPONENT
!   08-08     JANJIC     - Synchronize WATER array and Q.
!   08-11-23  janjic     - general hybrid coordinate
!
! USAGE: CALL RADIATION FROM PHY_RUN
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM
!$$$
!-----------------------------------------------------------------------
!
      INTEGER,PARAMETER :: JMS=1,JME=1
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER,INTENT(IN) :: LM,IMS,IME,MYPE                             &
                           ,IHRST,ITIMESTEP,JULDAY,JULYR                &
                           ,NPHS,NRADL,NRADS
!
      INTEGER,INTENT(IN) :: JDAT(8)
!
      REAL,INTENT(IN) :: DT,JULIAN,XTIME
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: CUPPT               &
                                                   ,ALBVB,ALBNB         & ! MODIS direct beam albedos
                                                   ,ALBVD,ALBND         & ! MODIS diffuse albedos
                                                   ,DYH                 &
                                                   ,GLAT,GLON           &
                                                   ,SM                  &
                                                   ,SNOW,SNOWC,THS,SI   &
                                                   ,TSKIN,Z0,SICE       &
                                                   ,MXSNAL,STDH
!
      REAL,DIMENSION(IMS:IME,JMS:JME,LM),INTENT(IN) :: Q,T,CW,O3        &
                                                      ,F_ICE,F_RAIN     &
                                                      ,F_RIMEF
!
      REAL,DIMENSION(IMS:IME,LM)  ,INTENT(IN) :: VVL,PRSL
      REAL,DIMENSION(IMS:IME,LM+1),INTENT(IN) :: PRSI
!
!-- Modify ALBEDO array if MODIS albedos are used in the RRTM
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
      REAL,DIMENSION(IMS:IME,JMS:JME,LM),INTENT(INOUT) :: RLWTT,RSWTT,QC,QI,QS,QR,QG,NI
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: CFRACH,CFRACL      &
                                                    ,CFRACM,CZMEAN      &
                                                    ,SIGT4
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(OUT) :: CLDFRA
!
      LOGICAL,INTENT(IN) :: F_QC,F_QR,F_QI,F_QS,F_QG,F_NI
!
      CHARACTER(99),INTENT(IN) :: LONGWAVE,SHORTWAVE,CLDFRACTION
!
!---------------------
!***  Local Variables
!---------------------
!
!.......................................................................
      INTEGER :: I,II,J,JDAY,JMONTH,K,NRAD
!
      INTEGER :: LW_PHYSICS=0,SW_PHYSICS=0,CLD_FRACTION
!
      INTEGER,DIMENSION(3) :: IDAT
!
      REAL :: DAYI,GMT,HOUR,PLYR,RADT,TIMES,TDUM,QIdum,QLdum
!
      real (kind=kind_phys) :: SLAG, SDEC, CDEC, SOLCON, DTSW, DTX
!
      REAL,DIMENSION(1:LM) :: QL
!
      REAL,DIMENSION(IMS:IME,JMS:JME) :: GSW                            &
     &                                  ,TOT,TSFC,XLAND                 &
     &                                  ,GLW,SWDOWN,SWDOWNC,CZEN        &
     &                                  ,CUPPTR,gridkm
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM+1) :: PHINT
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM) :: PI3D                      &
                                             ,THRATEN,THRATENLW         &
                                             ,THRATENSW                 &
                                             ,PRL,RHO,QV                &
                                             ,QCW,QCI,QSNOW,NCI         &
                                             ,QTdum,FIdum,FRdum
!
      LOGICAL :: GFDL_LW, GFDL_SW, LSSWR
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!*****
!***** NOTE: THIS IS HARDWIRED FOR CALLS TO LONGWAVE AND SHORTWAVE
!*****       AT EQUAL INTERVALS
!*****
!-----------------------------------------------------------------------
!
      NRAD=NRADS
      RADT=DT*NRADS/60.
!
!-----------------------------------------------------------------------
!***  NOTE:  THE NMMB HAS IJK STORAGE WITH LAYER 1 AT THE TOP.
!***         THE WRF PHYSICS DRIVERS HAVE IKJ STORAGE WITH LAYER 1
!***         AT THE BOTTOM.
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do private (j,i,k,plyr,ql)
!.......................................................................
      DO J=JMS,JME
      DO I=IMS,IME
!
        XLAND(I,J)=SM(I,J)+1.
!
!-----------------------------------------------------------------------
!***  FILL THE SINGLE-COLUMN INPUT
!-----------------------------------------------------------------------
!
        DO K=1,LM
!
          PLYR=PRSL(I,K)
!
          QL(K)=AMAX1(Q(I,J,K),EPSQ)
!
          PHINT(I,J,K)=PRSI(I,K)
          PI3D(I,J,K)=(PLYR*1.E-5)**CAPPA          ! Exner funtion
!
          THRATEN(I,J,K)=0.
          THRATENLW(I,J,K)=0.
          THRATENSW(I,J,K)=0.

          PRL(I,J,K)=PLYR                                     ! model layer pressure
          RHO(I,J,K)=PLYR/(R_D*T(I,J,K)*(1.+P608*Q(I,J,K)))   ! Air density (kg/m**3)
        ENDDO
!
        PHINT(I,J,LM+1)=PRSI(I,LM+1)
!
!-----------------------------------------------------------------------
!
        TSFC(I,J)=THS(I,J)*(PHINT(I,J,LM+1)*1.E-5)**CAPPA
!
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
      GMT=REAL(IHRST)
!
!.......................................................................
!$omp parallel do private (k,j,i)
!.......................................................................
        DO K=1,LM
        DO J=JMS,JME
        DO I=IMS,IME
          CLDFRA(I,J,K)=0.
        ENDDO
        ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!.......................................................................
!$omp parallel do private (j,i)
!.......................................................................
      DO J=JMS,JME
        DO I=IMS,IME
          CFRACH(I,J)=0.
          CFRACL(I,J)=0.
          CFRACM(I,J)=0.
          CZMEAN(I,J)=0.
          SIGT4(I,J)=0.
          SWDOWN(I,J)=0.    ! TOTAL (clear+cloudy sky) shortwave down at the surface
          SWDOWNC(I,J)=0.   ! CLEAR SKY shortwave down at the surface
          GSW(I,J)=0.       ! Net (down - up) total (clear+cloudy sky) shortwave at the surface
          GLW(I,J)=0.       ! Total longwave down at the surface
          CUPPTR(I,J)=CUPPT(I,J)   ! Temporary array set to zero in radiation
        ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  SYNCHRONIZE MIXING RATIO IN WATER ARRAY WITH SPECIFIC HUMIDITY.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!***  CALL THE INNER DRIVER.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!***  A PRIMARY MODIFICATION TO THE WRF DRIVER IS THE SPECIFICATION
!***  OF THE PACKAGES IN THE SELECT_CASE BLOCKS BEING CHANGED FROM
!***  INTEGERS (LISTED IN THE PHYSICS SECTION OF THE WRF REGISTRY)
!***  TO CHARACTERS (AS DEFINED IN THE ESMF CONFIG FILE).
!
!-----------------------------------------------------------------------
!***  TRANSLATE THE RADIATION OPTIONS IN THE CONFIG FILE TO THEIR
!***  ANALOGS IN THE WRF REGISTRY SO THAT THE WRF RADIATION DRIVER
!***  REMAINS UNTOUCHED.
!-----------------------------------------------------------------------
!
      SELECT CASE (TRIM(SHORTWAVE))
        CASE ('gfdl')
          SW_PHYSICS=99
        CASE ('dudh')
          SW_PHYSICS=1
        CASE ('gsfc')
          SW_PHYSICS=2
        CASE ('rrtm')
          SW_PHYSICS=3
        CASE DEFAULT
          WRITE(0,*)' User selected SHORTWAVE=',TRIM(SHORTWAVE)
          WRITE(0,*)' Improper selection of SW scheme in RADIATION'
          STOP ! or call abort
      END SELECT

      SELECT CASE (TRIM(LONGWAVE))
        CASE ('gfdl')
          LW_PHYSICS=99
        CASE ('rrtm')
          LW_PHYSICS=3
        CASE DEFAULT
          WRITE(0,*)' User selected LONGWAVE=',TRIM(LONGWAVE)
          WRITE(0,*)' Improper selection of LW scheme in RADIATION'
          STOP ! or call abort
      END SELECT

!==========================================================================
! Put "call radupdate_nmmb" here for threading safe
!==========================================================================
!---- for forcast purpose IDAT=JDAT

       DTSW =NRADS*DT                  ! [s]
       LSSWR=MOD(ITIMESTEP,NRADS)==0

      IF (SW_PHYSICS .EQ. 3 .or. LW_PHYSICS .EQ. 3) THEN
         DTX =DT
         call radupdate_nmmb                                          &
!  ---   inputs:
     &      ( JDAT, JDAT, DTSW, DTX, LSSWR, MYPE,                     &
!  ---   outputs:
     &        SLAG, SDEC, CDEC, SOLCON                                &
     &      )
      ENDIF

!==========================================================================
!==========================================================================


!
!-----------------------------------------------------------------------
!     CALL RADIATION_DRIVER
!-----------------------------------------------------------------------

   IF (LW_PHYSICS .EQ. 0 .AND. SW_PHYSICS .EQ. 0)         RETURN

   IF (ITIMESTEP .EQ. 1 .OR. MOD(ITIMESTEP,NRAD) .EQ. 0) THEN
     GFDL_LW = .FALSE.
     GFDL_SW = .FALSE.
!-----------------------------------------------------------------------
!***  Initialize Data
!-----------------------------------------------------------------------
!
     DO J=JMS,JME
     DO I=IMS,IME
        GSW(I,J)=0.
        GLW(I,J)=0.
        SWDOWN(I,J)=0.
     ENDDO
     ENDDO
!
     DO K=1,LM
     DO J=JMS,JME
     DO I=IMS,IME
         THRATEN(I,J,K)=0.
     ENDDO
     ENDDO
     ENDDO

!-----------------------------------------------------------------------
!
     lwrad_gfdl_select: SELECT CASE(lw_physics)
!
!-----------------------------------------------------------------------

        CASE (GFDLLWSCHEME)

!-- Do nothing, since cloud fractions (with partial cloudiness effects)
!-- are defined in GFDL LW/SW schemes and do not need to be initialized.

        CASE (RRTMLWSCHEME)

!-- Do nothing, since cloud fractions is calculated in RRTM

        CASE DEFAULT

          CALL CAL_CLDFRA(CLDFRA,QC,QI,F_QC,F_QI,IMS,IME,JMS,JME,LM)

     END SELECT lwrad_gfdl_select


!-----------------------------------------------------------------------
!
     lwrad_select: SELECT CASE(lw_physics)
!
!-----------------------------------------------------------------------
        CASE (RRTMLWSCHEME)


          !==== cloud fraction modification (HM Lin, 201503) ===========
          !
          !--- use Thompson cloud fraction

           cfr3_select: SELECT CASE (TRIM(CLDFRACTION))

              CASE ('thompson')        ! -- THOMPSON CLOUD FRACTION
                IF(MYPE==0)THEN
                  write(6,*) 'DEBUG-GT: using thompson cloud fraction scheme'
                ENDIF
                CLD_FRACTION=88
!
!--- Use dummy arrays QCW, QCI, QSNOW, NCI for Thompson cloud fraction scheme
!    These arrays are updated in cal_cldfra3, and the adjust cloud fields are
!    provided as input to RRTM and used by the radiation, but they are **not
!    used** to change (update) the model arrays QC, QI, QS, and NI (those
!    remain unchanged; BSF 4/13/2015).
!
                DO K=1,LM
                  DO J=JMS,JME
                    DO I=IMS,IME
                      QV(I,J,K)=Q(I,J,K)/(1.-Q(I,J,K))
                      QCW(I,J,K)=QC(I,J,K)
                      QCI(I,J,K)=0.
                      QSNOW(I,J,K)=0.
                      NCI(I,J,K)=0.
                      IF (F_QI) QCI(I,J,K)=QI(I,J,K)
                      IF (F_QS) QSNOW(I,J,K)=QS(I,J,K)
                      IF (F_NI) NCI(I,J,K)=NI(I,J,K)
                    ENDDO
                  ENDDO
                ENDDO
!
                DO J=JMS,JME
                DO I=IMS,IME
                  gridkm(i,j) = DYH(i,j)/1000.          ! convert m to km
                ENDDO
                ENDDO

                CALL cal_cldfra3(CLDFRA,                           &
                                 QV,QCW,QCI,QSNOW,F_NI,NCI,        &
                                 PRL,T,RHO,XLAND, gridkm,          & ! note:12.=gridkm is only for testing
                                 IMS,IME,JMS,JME,1,LM)
!
                DO K=1,LM
                  DO J=JMS,JME
                    DO I=IMS,IME
                      FIdum(I,J,K)=F_ICE(I,J,K)
                      FRdum(I,J,K)=F_RAIN(I,J,K)
                      QLdum=QCW(I,J,K)
                      QIdum=QCI(I,J,K)+QSNOW(I,J,K)
                      IF (F_QR) QLdum=QLdum+QR(I,J,K)
                      IF (F_QG) QIdum=QIdum+QG(I,J,K)
                      QTdum(I,J,K)=QLdum+QIdum
                      IF (QTdum(I,J,K)>0.) FIdum(I,J,K)=QIdum/QTdum(I,J,K)
                      IF (QLdum>0.) FRdum(I,J,K)=QR(I,J,K)/QLdum
                    ENDDO
                  ENDDO
                ENDDO

              CASE DEFAULT
                CLD_FRACTION=0

                DO K=1,LM
                  DO J=JMS,JME
                    DO I=IMS,IME
                      QTdum(I,J,K)=CW(I,J,K)
                      FIdum(I,J,K)=F_ICE(I,J,K)
                      FRdum(I,J,K)=F_RAIN(I,J,K)
                      QCW(I,J,K)=QC(I,J,K)
                      IF (F_QI) QCI(I,J,K)=QI(I,J,K)
                      IF (F_QS) QSNOW(I,J,K)=QS(I,J,K)
                      IF (F_NI) NCI(I,J,K)=NI(I,J,K)
                    ENDDO
                  ENDDO
                ENDDO

           END SELECT cfr3_select

!-----------------------------------------------------------------------
!
                CALL RRTM(ITIMESTEP,DT,JDAT                       &
                 ,NPHS,GLAT,GLON                                  &
                 ,NRADS,NRADL                                     &
                 ,PRSI,PRSL                                       &
                 ,T,Q,QTdum,O3                                    &  ! QTdum was CW
                 ,ALBEDO                                          &
                 ,ALBVB,ALBNB,ALBVD,ALBND                         &  ! MODIS albedos
                 ,FIdum,FRdum                                     &  ! FIdum,FRdum were F_ICE,F_RAIN
                 ,QCW,QCI,QSNOW,QR,QG,NCI                         &  ! QCW,QCI,QSNOW,NCI were QC,QI,QS,NI
                 ,F_QC,F_QI,F_QS,F_QR,F_QG,F_NI                   &
                 ,CLD_FRACTION                                    &
                 ,SM,CLDFRA                                       &
                 ,RLWTT,RSWTT                                     &
                 ,RLWIN,RSWIN                                     &
                 ,RSWINC,RLWINC,RSWOUC,RLWOUC,RLWOUCtoa,RSWOUCtoa &
                 ,RSWOUT                                          &
                 ,RLWTOA,RSWTOA                                   &
                 ,CZMEAN,SIGT4                                    &
                 ,CFRACL,CFRACM,CFRACH                            &
                 ,ACFRST                                          &
                 ,ACFRCV                                          &
                 ,CUPPT,SNOWC,SI                                  & ! was SNOW
                 ,HTOP,HBOT                                       &
                 ,TSKIN,Z0,SICE,F_RIMEF,MXSNAL,STDH,VVL           &
                 ,IMS,IME,JMS,JME                                 &
                 ,IMS,IME,JMS,JME                                 &
                 ,LM                                              &
                 ,SOLCON                                          &
                 ,MYPE )

        CASE (GFDLLWSCHEME)

                 gfdl_lw  = .true.
                 CALL GFDL(                                         &
                  DT=dt,XLAND=xland                                 &
                 ,PHINT=phint,T=t                                   &
                 ,Q=Q                                               &
                 ,QW=QC                                             &
                 ,QI=QI                                             &
                 ,QS=QS                                             &
                 ,F_QC=F_QC,F_QI=F_QI,F_QS=F_QS                     &
                 ,TSK2D=tsfc,GLW=GLW,RSWIN=SWDOWN,GSW=GSW           &
                 ,RSWINC=SWDOWNC,CLDFRA=CLDFRA,PI3D=PI3D            &
                 ,GLAT=glat,GLON=glon,HTOP=htop,HBOT=hbot           &
                 ,ALBEDO=albedo,CUPPT=cupptr                        &
                 ,SNOW=snow,G=g,GMT=gmt                             &
                 ,NSTEPRA=nrad,NPHS=nphs,ITIMESTEP=itimestep        &
                 ,XTIME=xtime,JULIAN=julian                         &
                 ,JULYR=julyr,JULDAY=julday                         &
                 ,GFDL_LW=gfdl_lw,GFDL_SW=gfdl_sw                   &
                 ,CFRACL=cfracl,CFRACM=cfracm,CFRACH=cfrach         &
                 ,ACFRST=acfrst                                     &
                 ,ACFRCV=acfrcv                                     &
                 ,RSWTOA=rswtoa,RLWTOA=rlwtoa,CZMEAN=czmean         &
                 ,THRATEN=thraten,THRATENLW=thratenlw               &
                 ,THRATENSW=thratensw                               &
                 ,IDS=ims,IDE=ime, JDS=jms,JDE=jme, KDS=1,KDE=lm+1  &
                 ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=1,KME=lm+1  &
                 ,ITS=IMS,ITE=IME, JTS=JMS,JTE=JME, KTS=1,KTE=lm    &
                                                                    )

        CASE DEFAULT

             WRITE(0,*)'The longwave option does not exist: lw_physics = ', lw_physics
             STOP ! or call abort

!-----------------------------------------------------------------------

     END SELECT lwrad_select

!-----------------------------------------------------------------------
!
     swrad_select: SELECT CASE(sw_physics)
!
!-----------------------------------------------------------------------

        CASE (SWRADSCHEME)
!!!          CALL SWRAD()

        CASE (GSFCSWSCHEME)
!!!          CALL GSFCSWRAD()

        CASE (RRTMSWSCHEME)

!-- Already called complete RRTM SW/LW scheme in LW part of driver
!!!          CALL RRTM()

        CASE (GFDLSWSCHEME)

                 gfdl_sw = .true.
                 CALL GFDL(                                         &
                  DT=dt,XLAND=xland                                 &
                 ,PHINT=phint,T=t                                   &
                 ,Q=Q                                               &
                 ,QW=QC                                             &
                 ,QI=QI                                             &
                 ,QS=QS                                             &
                 ,F_QC=F_QC,F_QI=F_QI,F_QS=F_QS                     &
                 ,TSK2D=tsfc,GLW=GLW,RSWIN=SWDOWN,GSW=GSW           &
                 ,RSWINC=SWDOWNC,CLDFRA=CLDFRA,PI3D=PI3D            &
                 ,GLAT=glat,GLON=glon,HTOP=htop,HBOT=hbot           &
                 ,ALBEDO=albedo,CUPPT=cupptr                        &
                 ,SNOW=snow,G=g,GMT=gmt                             &
                 ,NSTEPRA=nrad,NPHS=nphs,ITIMESTEP=itimestep        &
                 ,XTIME=xtime,JULIAN=julian                         &
                 ,JULYR=julyr,JULDAY=julday                         &
                 ,GFDL_LW=gfdl_lw,GFDL_SW=gfdl_sw                   &
                 ,CFRACL=cfracl,CFRACM=cfracm,CFRACH=cfrach         &
                 ,ACFRST=acfrst                                     &
                 ,ACFRCV=acfrcv                                     &
                 ,RSWTOA=rswtoa,RLWTOA=rlwtoa,CZMEAN=czmean         &
                 ,THRATEN=thraten,THRATENLW=thratenlw               &
                 ,THRATENSW=thratensw                               &
                 ,IDS=ims,IDE=ime, JDS=jms,JDE=jme, KDS=1,KDE=lm+1  &
                 ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=1,KME=lm+1  &
                 ,ITS=IMS,ITE=IME, JTS=JMS,JTE=JME, KTS=1,KTE=lm    &
                                                                    )

        CASE DEFAULT

             WRITE(0,*)'The shortwave option does not exist: sw_physics = ', sw_physics
             STOP ! or call abort

!-----------------------------------------------------------------------

     END SELECT swrad_select

!-----------------------------------------------------------------------
!
   ENDIF
!
!-----------------------------------------------------------------------
!
        IF(TRIM(SHORTWAVE)=='rrtm')THEN
          RETURN
        ENDIF
!
!-----------------------------------------------------------------------
!
!***  UPDATE FLUXES AND TEMPERATURE TENDENCIES.
!
!-----------------------------------------------------------------------
!***  SHORTWAVE
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      IF(MOD(ITIMESTEP,NRADS)==0)THEN
!-----------------------------------------------------------------------
!
        IF(TRIM(SHORTWAVE)/='gfdl')THEN
!
!-----------------------------------------------------------------------
!***  COMPUTE CZMEAN FOR NON-GFDL SHORTWAVE
!-----------------------------------------------------------------------
!
          DO J=JMS,JME
          DO I=IMS,IME
            CZMEAN(I,J)=0.
            TOT(I,J)=0.
          ENDDO
          ENDDO
!
          CALL CAL_MON_DAY(JULDAY,JULYR,JMONTH,JDAY)
          IDAT(1)=JMONTH
          IDAT(2)=JDAY
          IDAT(3)=JULYR
!
          DO II=0,NRADS,NPHS
            TIMES=ITIMESTEP*DT+II*DT
            CALL ZENITH(TIMES,DAYI,HOUR,IDAT,IHRST,GLON,GLAT,CZEN       &
     &                 ,IMS,IME,JMS,JME)
            DO J=JMS,JME
            DO I=IMS,IME
              IF(CZEN(I,J)>0.)THEN
                CZMEAN(I,J)=CZMEAN(I,J)+CZEN(I,J)
                TOT(I,J)=TOT(I,J)+1.
              ENDIF
            ENDDO
            ENDDO
!
          ENDDO
!
          DO J=JMS,JME
          DO I=IMS,IME
            IF(TOT(I,J)>0.)CZMEAN(I,J)=CZMEAN(I,J)/TOT(I,J)
          ENDDO
          ENDDO
!
!-----------------------------------------------------------------------
!***  COMPUTE TOTAL SFC SHORTWAVE DOWN FOR NON-GFDL SCHEMES
!-----------------------------------------------------------------------
!
          DO J=JMS,JME
          DO I=IMS,IME
!
            SWDOWN(I,J)=GSW(I,J)/(1.-ALBEDO(I,J))
!--- No value currently available for clear-sky solar fluxes from
!    non GFDL schemes, though it's needed for air quality forecasts.
!    For the time being, set to the total downward solar fluxes.
            SWDOWNC(I,J)=SWDOWN(I,J)
!
          ENDDO
          ENDDO
!
        ENDIF   !End non-GFDL/non-RRTM block
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do private(i,j,k)
!.......................................................................
        DO J=JMS,JME
          DO I=IMS,IME
!
            RSWIN(I,J)=SWDOWN(I,J)
            RSWINC(I,J)=SWDOWNC(I,J)
            RSWOUT(I,J)=SWDOWN(I,J)-GSW(I,J)
!
            DO K=1,LM
              RSWTT(I,J,K)=THRATENSW(I,J,K)*PI3D(I,J,K)
            ENDDO
!
          ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  LONGWAVE
!-----------------------------------------------------------------------
!
      IF(MOD(ITIMESTEP,NRADL)==0)THEN
!
!.......................................................................
!$omp parallel do private(i,j,k,tdum)
!.......................................................................
        DO J=JMS,JME
          DO I=IMS,IME
!
            TDUM=T(I,J,LM)
            SIGT4(I,J)=STBOLT*TDUM*TDUM*TDUM*TDUM
!
            DO K=1,LM
              RLWTT(I,J,K)=THRATENLW(I,J,K)*PI3D(I,J,K)
            ENDDO
!
            RLWIN(I,J)=GLW(I,J)
!
          ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE RADIATION
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
   SUBROUTINE cal_cldfra(CLDFRA,QC,QI,F_QC,F_QI,ims,ime,jms,jme,lm)
!---------------------------------------------------------------------
   IMPLICIT NONE
!---------------------------------------------------------------------
   INTEGER,  INTENT(IN   )   :: ims,ime, jms,jme, lm

!
   REAL, DIMENSION( ims:ime, jms:jme, 1:lm ), INTENT(OUT  ) ::       &
                                                             CLDFRA

   REAL, DIMENSION( ims:ime, jms:jme, 1:lm ), INTENT(IN   ) ::       &
                                                                 QI, &
                                                                 QC

   LOGICAL,INTENT(IN) :: F_QC,F_QI

   REAL thresh
   INTEGER:: i,j,k
! !DESCRIPTION:
! Compute cloud fraction from input ice and cloud water fields
! if provided.
!
! Whether QI or QC is active or not is determined from the logical
! switches f_qi and f_qc. They are passed in to the routine
! to enable testing to see if QI and QC represent active fields.
!
!---------------------------------------------------------------------
     thresh=1.0e-6

     IF ( f_qi .AND. f_qc ) THEN
!.......................................................................
!$omp parallel do private(k,j,i)
!.......................................................................
        DO k = 1,lm
        DO j = jms,jme
        DO i = ims,ime
           IF ( QC(i,j,k)+QI(I,j,k) .gt. thresh) THEN
              CLDFRA(i,j,k)=1.
           ELSE
              CLDFRA(i,j,k)=0.
           ENDIF
        ENDDO
        ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
     ELSE IF ( f_qc ) THEN
!.......................................................................
!$omp parallel do private(k,j,i)
!.......................................................................
        DO k = 1,lm
        DO j = jms,jme
        DO i = ims,ime
           IF ( QC(i,j,k) .gt. thresh) THEN
              CLDFRA(i,j,k)=1.
           ELSE
              CLDFRA(i,j,k)=0.
           ENDIF
        ENDDO
        ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
     ELSE
!
!.......................................................................
!$omp parallel do private(k,j,i)
!.......................................................................
        DO k = 1,lm
        DO j = jms,jme
        DO i = ims,ime
           CLDFRA(i,j,k)=0.
        ENDDO
        ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
     ENDIF

   END SUBROUTINE cal_cldfra
!---------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!---------------------------------------------------------------------
!
      END MODULE MODULE_RADIATION
!
!-----------------------------------------------------------------------
