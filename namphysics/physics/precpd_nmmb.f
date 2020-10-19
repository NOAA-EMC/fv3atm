!-- Changes:
!   1) SR array, "snow ratio" (fraction of frozen precipitation at the surface),
!      was added as output from the routine.
!   2) Replaces SR "index" (=-1 Snow, =0 Rain/Snow, =1 Rain), which wasn't used.
!
       SUBROUTINE PRECPD_nmmb (IM,IX,KM,DT,DEL,PRSL,PS,Q,CWM,T,RN,SR    &
     &,                   rainp,u00k,psautco,prautco,evpco,wminco       &
     &,                   lprnt,jpr)
!
!
!     ******************************************************************
!     *                                                                *
!     *           SUBROUTINE FOR PRECIPITATION PROCESSES               *
!     *           FROM SUSPENDED CLOUD WATER/ICE                       *
!     *                                                                *
!     ******************************************************************
!     *                                                                *
!     *  Originally CREATED BY  Q. ZHAO                JAN. 1995       *
!     *                         -------                                *    
!     *  Modified and rewritten by Shrinivas Moorthi   Oct. 1998       *
!     *                            -----------------                   *
!     *  and                       Hua-Lu Pan                          *
!     *                            ----------                          *
!     *                                                                *
!     *  References:                                                   *
!     *                                                                *
!     *  Zhao and Carr (1997), Monthly Weather Review (August)         *
!     *  Sundqvist et al., (1989) Monthly Weather review. (August)     *
!     *                                                                *
!     ******************************************************************
!
!     In this code vertical indexing runs from surface to top of the
!     model
!
!     Argument List:
!     --------------
!       IM         : Inner dimension over which calculation is made
!       IX         : Maximum inner dimension
!       KM         : Number of vertical levels
!       DT         : Time step in seconds
!       DEL(KM)    : Pressure layer thickness (Bottom to top)
!       PRSL(KM)   : Pressure values for model layers (bottom to top)
!       PS(IM)     : Surface pressure (centibars)
!       Q(IX,KM)   : Specific Humidity (Updated in the code)
!       CWM(IX,KM) : Condensate mixing ratio (Updated in the code)
!       T(IX,KM)   : Temperature       (Updated in the code)
!       RN(IM)     : Precipitation over one time-step DT (m/DT)
!old    SR(IM)     : Index (=-1 Snow, =0 Rain/Snow, =1 Rain)
!new    SR(IM)     : "Snow ratio", ratio of snow to total precipitation
!       TCW(IM)    : Vertically integrated liquid water (Kg/m**2)
!       CLL(IX,KM) : Cloud cover
!hchuang RN(IM) unit in m per time step
!        precipitation rate conversion 1 mm/s = 1 kg/m2/s
!
      USE MACHINE , ONLY : kind_phys
      USE FUNCPHYS , ONLY : fpvs
      USE PHYSCONS, grav => con_g, HVAP => con_HVAP, HFUS => con_HFUS   &
     &,             TTP => con_TTP, CP => con_CP                        &
     &,             EPS => con_eps, EPSM1 => con_epsm1
      implicit none
!     include 'constant.h'
!
      real (kind=kind_phys) G,      H1,    H2,   H1000                  &
     &,                     H1000G, D00,   D125, D5                     &
     &,                     ELWV,   ELIV,  ROW                          &
     &,                     EPSQ,   DLDT,  TM10, ELIW                   &
     &,                     RCP,    RROW
       PARAMETER (G=grav,         H1=1.E0,     H2=2.E0,     H1000=1000.0  &
     &,           H1000G=H1000/G, D00=0.E0,    D125=.125E0, D5=0.5E0      &
     &,           ELWV=HVAP,      ELIV=HVAP+HFUS,   ROW=1.E3              &
     &,           EPSQ=2.E-12,    DLDT=2274.E0,TM10=TTP-10.0              &
     &,           ELIW=ELIV-ELWV, RCP=H1/CP,   RROW=H1/ROW)

      real(kind=kind_phys), parameter :: cons_0=0.0,     cons_p01=0.01  &
     &,                                  cons_20=20.0                   &
     &,                                  cons_m30=-30.0, cons_50=50.0
!
      integer IM, IX, KM, LAT, jpr
      real (kind=kind_phys) Q(IX,KM),   T(IX,KM),    CWM(IX,KM)         &
     &,                                 DEL(IX,KM),  PRSL(IX,KM)        &
!    &,                     CLL(IM,KM), DEL(IX,KM),  PRSL(IX,KM)        &
     &,                     PS(IM),     RN(IM),      SR(IM)             &
     &,                     TCW(IM),    DT                              &
!hchuang code change [+1L] : add record to record information in vertical in
!                       addition to total column PRECRL
     &,                     RAINP(IM,KM), RNP(IM),                      &
     &                      psautco, prautco, evpco, wminco(2)
!
!
      real (kind=kind_phys) ERR(IM),      ERS(IM),     PRECRL(IM)       &
     &,                     PRECSL(IM),   PRECRL1(IM), PRECSL1(IM)      &
     &,                     RQ(IM),       CONDT(IM)                     &
     &,                     CONDE(IM),    RCONDE(IM),  TMT0(IM)         &
     &,                     WMIN(IM,KM),  WMINK(IM),   PRES(IM)         &
     &,                     WMINI(IM,KM), CCR(IM),     CCLIM(KM)        &
     &,                     TT(IM),       QQ(IM),      WW(IM)           &
     &,                     WFIX(KM),     U00K(IM,KM), ES(IM)           &
     &,                     Zaodt
!
      integer IW(IM,KM), IPR(IM), IWL(IM),     IWL1(IM)
!
       LOGICAL COMPUT(IM)
       logical lprnt
!
      real (kind=kind_phys) ke,   rdt,  us, cclimit, climit, cws, csm1  &
     &,                     crs1, crs2, cr, aa2,     dtcp,   c00, cmr   &
     &,                     tem,  c1,   c2, wwn                         &
!    &,                     tem,  c1,   c2, u00b,    u00t,   wwn        &
     &,                     precrk, precsk, pres1,   qk,     qw,  qi    &
     &,                     ai,     bi, qint, fiw, wws, cwmk, expf      &
     &,                     psaut, psaci, amaxcm, tem1, tem2            &
     &,                     tmt0k, tmt15, psm1, psm2, ppr               &
     &,                     rprs,  erk,   pps, sid, rid, amaxps         &
     &,                     praut, pracw, fi, qc, amaxrq, rqkll
      integer i, k, ihpr, n
!
!-----------------------Preliminaries ---------------------------------
!
!     DO K=1,KM
!       DO I=1,IM
!         CLL(I,K) = 0.0
!       ENDDO
!     ENDDO
!
      RDT     = H1 / DT
!     KE      = 2.0E-5  ! commented on 09/10/99  -- OPR value
!     KE      = 2.0E-6
!     KE      = 1.0E-5
!!!   KE      = 5.0E-5
!!    KE      = 7.0E-5
      KE      = evpco
!     KE      = 7.0E-5
      US      = H1
      CCLIMIT = 1.0E-3
      CLIMIT  = 1.0E-20
      CWS     = 0.025
!
      zaodt   = 800.0 * RDT
!
      CSM1    = 5.0000E-8   * zaodt
      CRS1    = 5.00000E-6  * zaodt
      CRS2    = 6.66600E-10 * zaodt
      CR      = 5.0E-4      * zaodt
      AA2     = 1.25E-3     * zaodt
!
      ke      = ke * sqrt(rdt)
!     ke      = ke * sqrt(zaodt)
!
      DTCP    = DT * RCP
!
!     C00 = 1.5E-1 * DT
!     C00 = 10.0E-1 * DT
!     C00 = 3.0E-1 * DT          !05/09/2000
!     C00 = 1.0E-4 * DT          !05/09/2000
      C00 = prautco * DT         !05/09/2000
      CMR = 1.0 / 3.0E-4
!     CMR = 1.0 / 5.0E-4
!     C1  = 100.0
      C1  = 300.0
      C2  = 0.5
!
!
!--------CALCULATE C0 AND CMR USING LC AT PREVIOUS STEP-----------------
!
      DO K=1,KM
        DO I=1,IM
          tem   = (prsl(i,k)*0.00001)
!         tem   = sqrt(tem)
          IW(I,K)    = 0.0
!         wmin(i,k)  = 1.0e-5 * tem
!         wmini(i,k) = 1.0e-5 * tem       ! Testing for RAS
!

          wmin(i,k)  = wminco(1) * tem
          wmini(i,k) = wminco(2) * tem


          rainp(i,k) = 0.0

        ENDDO
      ENDDO
      DO I=1,IM
!       C0(I)  = 1.5E-1
!       CMR(I) = 3.0E-4
!
        IWL1(I)    = 0
        PRECRL1(I) = D00
        PRECSL1(I) = D00
        COMPUT(I)  = .FALSE.
        RN(I)      = D00
        SR(I)      = D00
        ccr(i)     = D00
!
        RNP(I)     = D00
      ENDDO
!------------SELECT COLUMNS WHERE RAIN CAN BE PRODUCED--------------
      DO K=1, KM-1
        DO I=1,IM
          tem = min(wmin(i,k), wmini(i,k))
          IF (CWM(I,K) .GT. tem) COMPUT(I) = .TRUE.
        ENDDO
      ENDDO
      IHPR = 0
      DO I=1,IM
        IF (COMPUT(I)) THEN
           IHPR      = IHPR + 1
           IPR(IHPR) = I
        ENDIF
      ENDDO
!***********************************************************************
!-----------------BEGINING OF PRECIPITATION CALCULATION-----------------
!***********************************************************************
!     DO K=KM-1,2,-1
      DO K=KM,1,-1
        DO N=1,IHPR
          PRECRL(N) = PRECRL1(N)
          PRECSL(N) = PRECSL1(N)
          ERR  (N)  = D00
          ERS  (N)  = D00
          IWL  (N)  = 0
!
          I         = IPR(N)
          TT(N)     = T(I,K)
          QQ(N)     = Q(I,K)
          WW(N)     = CWM(I,K)
          WMINK(N)  = WMIN(I,K)
          PRES(N)   = prSL(I,K)
!
          PRECRK = MAX(cons_0,    PRECRL1(N))
          PRECSK = MAX(cons_0,    PRECSL1(N))
          WWN    = MAX(WW(N), CLIMIT)
!         IF (WWN .GT. WMINK(N) .OR. (PRECRK+PRECSK) .GT. D00) THEN
          IF (WWN .GT. CLIMIT .OR. (PRECRK+PRECSK) .GT. D00) THEN
            COMPUT(N) = .TRUE.
          ELSE
            COMPUT(N) = .FALSE.
          ENDIF
        ENDDO
!
!       es(1:IHPR) = fpvs(TT(1:IHPR))
        DO N=1,IHPR
          IF (COMPUT(N)) THEN
            I = IPR(N)
            CONDE(N)  = (DT/G) * DEL(I,K)
            CONDT(N)  = CONDE(N) * RDT
            RCONDE(N) = H1 / CONDE(N)
            QK        = MAX(EPSQ,  QQ(N))
            TMT0(N)   = TT(N) - 273.16
            WWN       = MAX(WW(N), CLIMIT)
!
!           PL = PRES(N) * 0.01
!           CALL QSATD(TT(N), PL, QC)
!           RQ(N) = MAX(QQ(N), EPSQ) / MAX(QC, 1.0E-10)
!           RQ(N) = MAX(1.0E-10, RQ(N))           ! -- RELATIVE HUMIDITY---
!
!  the global qsat computation is done in Pa
            pres1   = pres(n) 
!           QW      = es(N)
            QW      = min(pres1, fpvs(TT(N)))
            QW      = EPS * QW / (PRES1 + EPSM1 * QW)
            QW      = MAX(QW,EPSQ)
!
!           TMT15 = MIN(TMT0(N), cons_m15)
!           AI    = 0.008855
!           BI    = 1.0
!           IF (TMT0(N) .LT. -20.0) THEN
!             AI = 0.007225
!             BI = 0.9674
!           ENDIF
!           QI   = QW * (BI + AI*MIN(TMT0(N),cons_0))
!           QINT = QW * (1.-0.00032*TMT15*(TMT15+15.))
!
            qi   = qw
            qint = qw
!           IF (TMT0(N).LE.-40.) QINT = QI
!
!-------------------ICE-WATER ID NUMBER IW------------------------------
            IF(TMT0(N).LT.-15.) THEN
               FI = QK - U00K(I,K)*QI
               IF(FI.GT.D00.OR.WWN.GT.CLIMIT) THEN
                  IWL(N) = 1
               ELSE
                  IWL(N) = 0
               ENDIF
!           ENDIF
            ELSEIF (TMT0(N).GE.0.) THEN
               IWL(N) = 0
!
!           IF(TMT0(N).LT.0.0.AND.TMT0(N).GE.-15.0) THEN
            ELSE
              IWL(N) = 0
              IF(IWL1(N).EQ.1.AND.WWN.GT.CLIMIT) IWL(N)=1
            ENDIF
!
!           IF(TMT0(N).GE.0.) THEN
!              IWL(N) = 0
!           ENDIF
!----------------THE SATUATION SPECIFIC HUMIDITY------------------------
            FIW   = FLOAT(IWL(N))
            QC    = (H1-FIW)*QINT + FIW*QI
!----------------THE RELATIVE HUMIDITY----------------------------------
            IF(QC .LE. 1.0E-10) THEN
               RQ(N) = D00
            ELSE
               RQ(N) = QK / QC
            ENDIF
!----------------CLOUD COVER RATIO CCR----------------------------------
            IF(RQ(N).LT.U00K(I,K)) THEN
                   CCR(N)=D00
            ELSEIF(RQ(N).GE.US) THEN
                   CCR(N)=US
            ELSE
                 RQKLL=MIN(US,RQ(N))
                 CCR(N)= H1-SQRT((US-RQKLL)/(US-U00K(I,K)))
            ENDIF
!
          ENDIF
        ENDDO
!-------------------ICE-WATER ID NUMBER IWL------------------------------
!       DO N=1,IHPR
!         IF (COMPUT(N) .AND.  (WW(N) .GT. CLIMIT)) THEN
!           IF (TMT0(N) .LT. -15.0
!    *         .OR. (TMT0(N) .LT. 0.0 .AND. IWL1(N) .EQ. 1))
!    *                                      IWL(N) = 1
!             CLL(IPR(N),K) = 1.0                           ! Cloud Cover!
!             CLL(IPR(N),K) = MIN(1.0, WW(N)*CCLIM(K))      ! Cloud Cover!
!         ENDIF
!       ENDDO
!
!---   PRECIPITATION PRODUCTION --  Auto Conversion and Accretion
!
        DO N=1,IHPR
          IF (COMPUT(N) .AND. CCR(N) .GT. 0.0) THEN
            WWS    = WW(N)
            CWMK   = MAX(cons_0, WWS)
!           AMAXCM = MAX(cons_0, CWMK - WMINK(N))
            IF (IWL(N) .EQ. 1) THEN                 !  Ice Phase
               AMAXCM = MAX(cons_0, CWMK - WMINI(IPR(N),K))
               EXPF      = DT * EXP(0.025*TMT0(N))
               PSAUT     = MIN(CWMK, psautco*EXPF*AMAXCM)

!              PSAUT     = MIN(CWMK, 2.0E-3*EXPF*AMAXCM)
!              PSAUT     = MIN(CWMK, 1.0E-3*EXPF*AMAXCM)
!              PSAUT     = MIN(CWMK, 7.5E-4*EXPF*AMAXCM)
!!!!!!!        PSAUT     = MIN(CWMK, 7.0E-4*EXPF*AMAXCM)
!b             PSAUT     = MIN(CWMK, 6.5E-4*EXPF*AMAXCM)
!!!!           PSAUT     = MIN(CWMK, 6.0E-4*EXPF*AMAXCM)
!              PSAUT     = MIN(CWMK, 5.0E-4*EXPF*AMAXCM)
!              PSAUT     = MIN(CWMK, 4.0E-4*EXPF*AMAXCM)

               WW(N)     = WW(N) - PSAUT
               CWMK      = MAX(cons_0, WW(N))
!              CWMK      = MAX(cons_0, WW(N)-wmini(ipr(n),k))
               PSACI     = MIN(CWMK, AA2*EXPF*PRECSL1(N)*CWMK)

               WW(N)     = WW(N) - PSACI
 
               PRECSL(N) = PRECSL(N) + (WWS - WW(N)) * CONDT(N)
            ELSE                                    !  Liquid Water
!
!          For using Sundqvist precip formulation of rain
!
               AMAXCM    = MAX(cons_0, CWMK - WMINK(N))
!!             AMAXCM    = CWMK
               TEM1      = PRECSL1(N) + PRECRL1(N)
               TEM2      = MIN(MAX(cons_0, 268.0-TT(N)), cons_20)
               TEM       = (1.0+C1*SQRT(TEM1*RDT)) * (1+C2*SQRT(TEM2))
!
               TEM2      = AMAXCM * CMR * TEM / max(CCR(N),cons_p01)
               TEM2      = MIN(cons_50, TEM2*TEM2)
               PRAUT     = C00  * TEM * AMAXCM * (1.0-EXP(-TEM2))
               PRAUT     = MIN(PRAUT, CWMK)
               WW(N)     = WW(N) - PRAUT
!
!          Below is for Zhao's precip formulation (water)
!
!              AMAXCM    = MAX(cons_0, CWMK - WMINK(N))
!              PRAUT     = MIN(CWMK, C00*AMAXCM*AMAXCM)
!              WW(N)     = WW(N) - PRAUT
!
!              CWMK      = MAX(cons_0, WW(N))
!              TEM1      = PRECSL1(N) + PRECRL1(N)
!              PRACW     = MIN(CWMK, CR*DT*TEM1*CWMK)
!              WW(N)     = WW(N) - PRACW
!
               PRECRL(N) = PRECRL(N) + (WWS - WW(N)) * CONDT(N)
!
!hchuang code change [+1L] : add record to record information in vertical
! TURN RNP in unit of WW (CWM and Q, kg/kg ???)
               RNP(N) = RNP(N) + (WWS - WW(N))
            ENDIF
          ENDIF
        ENDDO
!
!-----EVAPORATION OF PRECIPITATION-------------------------
!**** ERR & ERS POSITIVE--->EVAPORATION-- NEGTIVE--->CONDENSATION
!
        DO N=1,IHPR
          IF (COMPUT(N)) THEN
            I      = IPR(N)
            QK     = MAX(EPSQ,  QQ(N))
            TMT0K  = MAX(cons_m30, TMT0(N))
            PRECRK = MAX(cons_0,    PRECRL(N))
            PRECSK = MAX(cons_0,    PRECSL(N))
            AMAXRQ = MAX(cons_0,    U00K(I,K)-RQ(N)) * CONDE(N)
!----------------------------------------------------------------------
! INCREASE THE EVAPORATION FOR STRONG/LIGHT PREC
!----------------------------------------------------------------------
            PPR    = KE * AMAXRQ * SQRT(PRECRK)
!           PPR    = KE * AMAXRQ * SQRT(PRECRK*RDT)
            IF (TMT0(N) .GE. 0.) THEN
              PPS = 0.
            ELSE
              PPS = (CRS1+CRS2*TMT0K) * AMAXRQ * PRECSK / U00K(I,K)
            END IF
!---------------CORRECT IF OVER-EVAPO./COND. OCCURS--------------------
            ERK=PRECRK+PRECSK
            IF(RQ(N).GE.1.0E-10)  ERK = AMAXRQ * QK * RDT / RQ(N)
            IF (PPR+PPS .GT. ABS(ERK)) THEN
               RPRS   = ERK / (PRECRK+PRECSK)
               PPR    = PRECRK * RPRS
               PPS    = PRECSK * RPRS
            ENDIF
            PPR       = MIN(PPR, PRECRK)
            PPS       = MIN(PPS, PRECSK)
            ERR(N)    = PPR * RCONDE(N)
            ERS(N)    = PPS * RCONDE(N)
            PRECRL(N) = PRECRL(N) - PPR
!hchuang code change [+1L] : add record to record information in vertical
! Use ERR for kg/kg/DT not the PPR (mm/DT=kg/m2/DT)
!
            RNP(N) = RNP(N) - ERR(N)
!
            PRECSL(N) = PRECSL(N) - PPS
          ENDIF
        ENDDO
!--------------------MELTING OF THE SNOW--------------------------------
        DO N=1,IHPR
          IF (COMPUT(N)) THEN
            IF (TMT0(N) .GT. 0.) THEN
               AMAXPS = MAX(cons_0,    PRECSL(N))
               PSM1   = CSM1 * TMT0(N) * TMT0(N) * AMAXPS
               PSM2   = CWS * CR * MAX(cons_0, WW(N)) * AMAXPS
               PPR    = (PSM1 + PSM2) * CONDE(N)
               IF (PPR .GT. AMAXPS) THEN
                 PPR  = AMAXPS
                 PSM1 = AMAXPS * RCONDE(N)
               ENDIF
               PRECRL(N) = PRECRL(N) + PPR
!
!hchuang code change [+1L] : add record to record information in vertical
! TURN PPR (mm/DT=kg/m2/DT) to kg/kg/DT -> PPR/air density (kg/m3)
               RNP(N) = RNP(N) + PPR * RCONDE(N)
!
               PRECSL(N) = PRECSL(N) - PPR
            ELSE
               PSM1 = D00
            ENDIF
!
!---------------UPDATE T AND Q------------------------------------------
            TT(N) = TT(N) - DTCP * (ELWV*ERR(N)+ELIV*ERS(N)+ELIW*PSM1)
            QQ(N) = QQ(N) + DT * (ERR(N)+ERS(N))
          ENDIF
        ENDDO
!
        DO N=1,IHPR
          IWL1(N)    = IWL(N)
          PRECRL1(N) = MAX(cons_0, PRECRL(N))
          PRECSL1(N) = MAX(cons_0, PRECSL(N))
          I          = IPR(N)
          T(I,K)     = TT(N)
          Q(I,K)     = QQ(N)
          CWM(I,K)   = WW(N)
          IW(I,K)    = IWL(N)
!hchuang code change [+1L] : add record to record information in vertical
! RNP = PRECRL1*RCONDE(N) unit in kg/kg/DT
!
          RAINP(I,K) = RNP(N)
        ENDDO
!
!  move water from vapor to liquid should the liquid amount be negative
!
        do i = 1, im
          if (cwm(i,k) < 0.) then
            tem      = q(i,k) + cwm(i,k)
            if (tem >= 0.0) then
              q(i,k)   = tem
              t(i,k)   = t(i,k) - elwv * rcp * cwm(i,k)
              cwm(i,k) = 0.
            elseif (q(i,k) > 0.0) then
              cwm(i,k) = tem
              t(i,k)   = t(i,k) + elwv * rcp * q(i,k)
              q(i,k)   = 0.0
            endif
          endif
        enddo
!
      ENDDO                               ! K loop ends here!
!**********************************************************************
!-----------------------END OF PRECIPITATION PROCESSES-----------------
!**********************************************************************
!
      DO N=1,IHPR
        I = IPR(N)
        RN(I) = (PRECRL1(N)  + PRECSL1(N)) * RROW  ! Precip at surface
!old!
!old!----SR=1 IF SFC PREC IS RAIN ; ----SR=-1 IF SFC PREC IS SNOW
!old!----SR=0 FOR BOTH OF THEM OR NO SFC PREC
!old!
!old        RID = 0.
!old        SID = 0.
!old        IF (PRECRL1(N) .GE. 1.E-13) RID = 1.
!old        IF (PRECSL1(N) .GE. 1.E-13) SID = -1.
!old        SR(I) = RID + SID  ! SR=1 --> Rain, SR=-1 -->Snow, SR=0 -->Both
!
!new: SR = 'snow ratio', fraction of frozen precipitation (1 time step)
!
        RID=PRECRL1(N)+PRECSL1(N)
        IF (RID<1.E-13) THEN
           SR(I)=0.
        ELSE
           SR(I)=PRECSL1(N)/RID
        ENDIF

      ENDDO
!
      RETURN
      END
