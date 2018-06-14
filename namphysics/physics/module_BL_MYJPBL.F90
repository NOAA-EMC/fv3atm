!-----------------------------------------------------------------------
!
      MODULE MODULE_BL_MYJPBL
!
!-----------------------------------------------------------------------
!
!***  THE MYJ PBL SCHEME
!
!-----------------------------------------------------------------------
!
      USE MODULE_KINDS
!
      USE MODULE_CONSTANTS,ONLY : A2,A3,A4,CP,ELIV,ELWV,ELIWV &
                                 ,EP_1,EPSQ  &
                                 ,G,P608,PI,PQ0,R_D,R_V,RHOWATER        &
                                 ,STBOLT,CAPPA
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC:: MYJPBL_INIT, MYJPBL
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  FOR MYJ TURBULENCE
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      REAL(KIND=KFPT),PARAMETER:: &
       ELEVFC=0.6
!
      REAL(KIND=KFPT),PARAMETER:: &
       VKARMAN=0.4 &
!
      ,XLS=ELIV,XLV=ELWV &
      ,RLIVWV=XLS/XLV,ELOCP=2.72E6/CP &
!
      ,EPS1=1.E-12,EPS2=0. &
      ,EPSRU=1.E-7,EPSRS=1.E-7 &
      ,EPSTRB=1.E-24 &
      ,FH=1.10 &
!
      ,ALPH=0.30,BETA=1./273.,EL0MAX=1000.,EL0MIN=1. &
!      ,ELFC=0.5,GAM1=0.2222222222222222222 &
!      ,ELFC=0.23*0.25,GAM1=0.2222222222222222222 &
      ,ELFC=1.,GAM1=0.2222222222222222222 &
!
      ,A1=0.659888514560862645 &
      ,A2X=0.6574209922667784586 &
      ,B1=11.87799326209552761 &
      ,B2=7.226971804046074028 &
      ,C1=0.000830955950095854396 &
      ,ELZ0=0.,ESQ=5.0 &
!
      ,SEAFC=0.98,PQ0SEA=PQ0*SEAFC &
!
      ,BTG=BETA*G &
      ,ESQHF=0.5*5.0 &
      ,RB1=1./B1
!
      REAL(KIND=KFPT),PARAMETER:: &
       ADNH= 9.*A1*A2X*A2X*(12.*A1+3.*B2)*BTG*BTG &
      ,ADNM=18.*A1*A1*A2X*(B2-3.*A2X)*BTG &
      ,ANMH=-9.*A1*A2X*A2X*BTG*BTG &
      ,ANMM=-3.*A1*A2X*(3.*A2X+3.*B2*C1+18.*A1*C1-B2)*BTG &
      ,BDNH= 3.*A2X*(7.*A1+B2)*BTG &
      ,BDNM= 6.*A1*A1 &
      ,BEQH= A2X*B1*BTG+3.*A2X*(7.*A1+B2)*BTG &
      ,BEQM=-A1*B1*(1.-3.*C1)+6.*A1*A1 &
      ,BNMH=-A2X*BTG &
      ,BNMM=A1*(1.-3.*C1) &
      ,BSHH=9.*A1*A2X*A2X*BTG &
      ,BSHM=18.*A1*A1*A2X*C1 &
      ,BSMH=-3.*A1*A2X*(3.*A2X+3.*B2*C1+12.*A1*C1-B2)*BTG &
      ,CESH=A2X &
      ,CESM=A1*(1.-3.*C1) &
      ,CNV=EP_1*G/BTG
!
!-----------------------------------------------------------------------
!***  FREE TERM IN THE EQUILIBRIUM EQUATION FOR (L/Q)**2
!-----------------------------------------------------------------------
!
      REAL(KIND=KFPT),PARAMETER:: &
       AEQH=9.*A1*A2X*A2X*B1*BTG*BTG &
           +9.*A1*A2X*A2X*(12.*A1+3.*B2)*BTG*BTG &
      ,AEQM=3.*A1*A2X*B1*(3.*A2X+3.*B2*C1+18.*A1*C1-B2) &
           *BTG+18.*A1*A1*A2X*(B2-3.*A2X)*BTG
!
!-----------------------------------------------------------------------
!***  FORBIDDEN TURBULENCE AREA
!-----------------------------------------------------------------------
!
      REAL(KIND=KFPT),PARAMETER:: &
       REQU=-AEQH/AEQM &
      ,EPSGH=1.E-9,EPSGM=REQU*EPSGH
!
!-----------------------------------------------------------------------
!***  NEAR ISOTROPY FOR SHEAR TURBULENCE, WW/Q2 LOWER LIMIT
!-----------------------------------------------------------------------
!
      REAL(KIND=KFPT),PARAMETER:: &
       UBRYL=(18.*REQU*A1*A1*A2X*B2*C1*BTG &
             +9.*A1*A2X*A2X*B2*BTG*BTG) &
            /(REQU*ADNM+ADNH) &
      ,UBRY=(1.+EPSRS)*UBRYL,UBRY3=3.*UBRY
!
      REAL(KIND=KFPT),PARAMETER:: &
       AUBH=27.*A1*A2X*A2X*B2*BTG*BTG-ADNH*UBRY3 &
      ,AUBM=54.*A1*A1*A2X*B2*C1*BTG  -ADNM*UBRY3 &
      ,BUBH=(9.*A1*A2X+3.*A2X*B2)*BTG-BDNH*UBRY3 &
      ,BUBM=18.*A1*A1*C1             -BDNM*UBRY3 &
      ,CUBR=1.                       -     UBRY3 &
      ,RCUBR=1./CUBR
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---LOOK-UP TABLES------------------------------------------------------
INTEGER(KIND=KINT),PARAMETER:: &
 ITBL=401 &                  ! CONVECTION TABLES, DIMENSION 1
,JTBL=1201 &                 ! CONVECTION TABLES, DIMENSION 2
,KERFM=301 &                 ! SIZE OF ERF HALF TABLE
,KERFM2=KERFM-2              ! INTERNAL POINTS OF ERF HALF TABLE

REAL(KIND=KFPT),PARAMETER:: &
 PL=2500. &                  ! LOWER BOUND OF PRESSURE RANGE
,PH=105000. &                ! UPPER BOUND OF PRESSURE RANGE
,THL=210. &                  ! LOWER BOUND OF POTENTIAL TEMPERATURE RANGE
,THH=365. &                  ! UPPER BOUND OF POTENTIAL TEMPERATURE RANGE
,XEMIN=0. &                  ! LOWER BOUND OF ERF HALF TABLE
,XEMAX=3.                    ! UPPER BOUND OF ERF HALF TABLE

REAL(KIND=KFPT),PRIVATE,SAVE:: &
 RDP &                       ! SCALING FACTOR FOR PRESSURE
,RDQ &                       ! SCALING FACTOR FOR HUMIDITY
,RDTH &                      ! SCALING FACTOR FOR POTENTIAL TEMPERATURE
,RDTHE &                     ! SCALING FACTOR FOR EQUIVALENT POT. TEMPERATURE
,RDXE                        ! ERF HALF TABLE SCALING FACTOR

REAL(KIND=KFPT),DIMENSION(1:ITBL),PRIVATE,SAVE:: &
 STHE &                      ! RANGE FOR EQUIVALENT POTENTIAL TEMPERATURE
,THE0                        ! BASE FOR EQUIVALENT POTENTIAL TEMPERATURE           

REAL(KIND=KFPT),DIMENSION(1:JTBL),PRIVATE,SAVE:: &
 QS0 &                       ! BASE FOR SATURATION SPECIFIC HUMIDITY
,SQS                         ! RANGE FOR SATURATION SPECIFIC HUMIDITY

REAL(KIND=KFPT),DIMENSION(1:KERFM),PRIVATE,SAVE:: &
 HERFF                       ! HALF ERF TABLE

REAL(KIND=KFPT),DIMENSION(1:ITBL,1:JTBL),PRIVATE,SAVE:: &
 PTBL                        ! SATURATION PRESSURE TABLE

REAL(KIND=KFPT),DIMENSION(1:JTBL,1:ITBL),PRIVATE,SAVE:: &
 TTBL                        ! TEMPERATURE TABLE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
! REFERENCES:  JANJIC (2001), NCEP OFFICE NOTE 437
!
! ABSTRACT:
!     MYJ UPDATES THE TURBULENT KINETIC ENERGY WITH THE PRODUCTION/
!     DISSIPATION TERM AND THE VERTICAL DIFFUSION TERM
!     (USING AN IMPLICIT FORMULATION) FROM MELLOR-YAMADA
!     LEVEL 2.5 AS EXTENDED BY JANJIC.  EXCHANGE COEFFICIENTS FOR
!     THE SURFACE LAYER ARE COMPUTED FROM THE MONIN-OBUKHOV THEORY.
!     THE TURBULENT VERTICAL EXCHANGE IS THEN EXECUTED.
!
!-----------------------------------------------------------------------
      SUBROUTINE MYJPBL(DT,NPHS,EPSL,EPSQ2,HT,STDH,DZ &
                       ,PMID,PINH,TH,T,EXNER,Q,CWM,U,V &
                       ,TSK,QSFC,CHKLOWQ,THZ0,QZ0,UZ0,VZ0 &
                       ,XLAND,SICE,SNOW &
                       ,Q2,       USTAR,Z0,EL_MYJ,PBLH,KPBL,CT &
!                      ,Q2,EXCH_H,USTAR,Z0,EL_MYJ,PBLH,KPBL,CT &
                       ,AKHS,AKMS,ELFLX,MIXHT &
                       ,RUBLTEN,RVBLTEN,RTHBLTEN,RQBLTEN,RQCBLTEN &
                       ,IMS,IME,JMS,JME,LM)
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER(KIND=KINT),INTENT(IN):: &
       IMS,IME,JMS,JME,LM
!
      INTEGER(KIND=KINT),INTENT(IN):: &
       NPHS
!
      INTEGER(KIND=KINT),DIMENSION(IMS:IME,JMS:JME),INTENT(OUT):: &
       KPBL
!
      REAL(KIND=KFPT),INTENT(IN):: &
       DT
!
      real(kind=kfpt),dimension(1:lm-1),intent(in):: EPSL
      real(kind=kfpt),dimension(1:lm),intent(in):: EPSQ2
!
      REAL(KIND=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(IN):: &
       HT,SICE,SNOW,STDH  &
      ,TSK,XLAND &
      ,CHKLOWQ,ELFLX
!
      REAL(KIND=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN):: &
       DZ,EXNER,PMID,Q,CWM,U,V,T,TH
!
      REAL(KIND=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM+1),INTENT(IN):: &
       PINH
!
      REAL(KIND=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(OUT):: &
       MIXHT &
      ,PBLH
!
      REAL(KIND=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(OUT):: &
       EL_MYJ
!
      REAL(KIND=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(OUT):: &
       RQCBLTEN &
      ,RUBLTEN,RVBLTEN &
      ,RTHBLTEN,RQBLTEN
!
      REAL(KIND=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT):: &
       AKHS,AKMS
!
      REAL(KIND=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT):: &
       CT,QSFC,QZ0     &
      ,THZ0,USTAR      &
      ,UZ0,VZ0,Z0
!
      REAL(KIND=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT):: &
       Q2
!      EXCH_H &
!     ,Q2
!
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER(KIND=KINT):: &
       I,IQTB,ITTB,J,K,LLOW,LMH,LMXL
!
      INTEGER(KIND=KINT),DIMENSION(IMS:IME,JMS:JME):: &
       LPBL
!
      REAL(KIND=KFPT):: &
       AKHS_DENS,AKMS_DENS,BQ,BQS00K,BQS10K &
      ,DCDT,DELTAZ,DQDT,DTDIF,DTDT,DTTURBL &
      ,P00K,P01K,P10K,P11K,PELEVFC,PP1,PSFC,PSP,PTOP &
      ,QBT,QFC1,QLOW,QQ1,QX &
      ,RDTTURBL,RG,RSQDT,RXNERS,RXNSFC &
      ,SEAMASK,SQ,SQS00K,SQS10K &
      ,THBT,THNEW,THOLD,TQ,TTH &
      ,ULOW,VLOW,RSTDH,STDFAC,ZSF,ZSX,ZSY,ZUV
!
      REAL(KIND=KFPT),DIMENSION(1:LM):: &
       CWMK,PK,PSK,Q2K,QK,RHOK,RXNERK,THEK,THK,THVK,TK,UK,VK
!
      REAL(KIND=KFPT),DIMENSION(1:LM-1):: &
       AKHK,AKMK,DCOL,EL,GH,GM
!
      REAL(KIND=KFPT),DIMENSION(1:LM+1):: &
       ZHK
!
      REAL(KIND=KFPT),DIMENSION(IMS:IME,JMS:JME):: &
       THSK
!
      REAL(KIND=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM):: &
       RXNER,THV
!
      REAL(KIND=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM-1):: &
       AKH,AKM
!
      REAL(KIND=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM+1):: &
       ZINT
!
!***  Begin debugging
      REAL(KIND=KFPT):: ZSL_DIAG
      INTEGER(KIND=KINT):: IMD,JMD,PRINT_DIAG
!***  End debugging
!
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!
!***  Begin debugging
      IMD=(IMS+IME)/2
      JMD=(JMS+JME)/2
!***  End debugging
!
!***  MAKE PREPARATIONS
!
!----------------------------------------------------------------------
      STDFAC=1.
!----------------------------------------------------------------------
      DTTURBL=DT*NPHS
      RDTTURBL=1./DTTURBL
      RSQDT=SQRT(RDTTURBL)
      DTDIF=DTTURBL
      RG=1./G
!
      DO K=1,LM-1
      DO J=JMS,JME
      DO I=IMS,IME
        AKM(I,J,K)=0.
      ENDDO
      ENDDO
      ENDDO
!
      DO K=1,LM+1
      DO J=JMS,JME
      DO I=IMS,IME
        ZINT(I,J,K)=0.
      ENDDO
      ENDDO
      ENDDO
!
      DO J=JMS,JME
      DO I=IMS,IME
        ZINT(I,J,LM+1)=HT(I,J)     ! Z AT BOTTOM OF LOWEST SIGMA LAYER
      ENDDO
      ENDDO
!
      DO K=LM,1,-1
        DO J=JMS,JME
          DO I=IMS,IME
            ZINT(I,J,K)=ZINT(I,J,K+1)+DZ(I,J,K)
            RXNER(I,J,K)=1./EXNER(I,J,K)
            THV(I,J,K)=(Q(I,J,K)*0.608+(1.-CWM(I,J,K)))*TH(I,J,K)
          ENDDO
        ENDDO
      ENDDO
!
      DO J=JMS,JME
        DO I=IMS,IME
          EL_MYJ(I,J,LM)=0.
        ENDDO
      ENDDO
!
!----------------------------------------------------------------------
!.......................................................................
!ZJ$OMP PARALLEL DO &
!ZJ$OMP PRIVATE(J,I,LMH,PTOP,PSFC,SEAMASK,K,TK,THVK,QK,Q2K,RXNERK,     &
!ZJ$OMP         PK,UK,VK,Q2K,ZHK,LMXL,GM,GH,EL,AKMK,AKHK,DELTAZ),   &
!ZJ$OMP         SCHEDULE(DYNAMIC)
!.......................................................................
!----------------------------------------------------------------------
      setup_integration:  DO J=JMS,JME
!----------------------------------------------------------------------
!
        DO I=IMS,IME
!
          LMH=LM
!
          PTOP=PINH(I,J,1)
          PSFC=PINH(I,J,LMH+1)
!
!***  CONVERT LAND MASK (1 FOR SEA; 0 FOR LAND)
!
          SEAMASK=XLAND(I,J)-1.
!
!***  FILL 1-D VERTICAL ARRAYS
!
          DO K=LM,1,-1
            PK(K)=PMID(I,J,K)
            TK(K)=T(I,J,K)
            QK(K)=Q(I,J,K)
            THVK(K)=THV(I,J,K)
            RXNERK(K)=RXNER(I,J,K)
            UK(K)=U(I,J,K)
            VK(K)=V(I,J,K)
            Q2K(K)=Q2(I,J,K)
!
!***  COMPUTE THE HEIGHTS OF THE LAYER INTERFACES
!
            ZHK(K)=ZINT(I,J,K)
!
          ENDDO
          ZHK(LM+1)=HT(I,J)          ! Z AT BOTTOM OF LOWEST SIGMA LAYER
!
!***  POTENTIAL INSTABILITY
!
          PELEVFC=PMID(I,J,LMH)*ELEVFC
!
          DO K=LMH,1,-1
!-----------------------------------------------------------------------
            IF(K==LMH .OR. PMID(I,J,K)>PELEVFC) THEN
!---PREPARATION FOR SEARCH FOR MAX CAPE---------------------------------
              QBT=QK(K)
              THBT=TH(I,J,K)
              TTH=(THBT-THL)*RDTH
              QQ1=TTH-AINT(TTH)
              ITTB=INT(TTH)+1
!---KEEPING INDICES WITHIN THE TABLE------------------------------------
              IF(ITTB.LT.1)THEN
                ITTB=1
                QQ1=0.
              ELSE IF(ITTB.GE.JTBL)THEN
                ITTB=JTBL-1
                QQ1=0.
              ENDIF
!---BASE AND SCALING FACTOR FOR SPEC. HUMIDITY--------------------------
              BQS00K=QS0(ITTB)
              SQS00K=SQS(ITTB)
              BQS10K=QS0(ITTB+1)
              SQS10K=SQS(ITTB+1)
!--------------SCALING SPEC. HUMIDITY & TABLE INDEX---------------------
              BQ=(BQS10K-BQS00K)*QQ1+BQS00K
              SQ=(SQS10K-SQS00K)*QQ1+SQS00K
              TQ=(QBT-BQ)/SQ*RDQ
              PP1=TQ-AINT(TQ)
              IQTB=INT(TQ)+1
!----------------KEEPING INDICES WITHIN THE TABLE-----------------------
              IF(IQTB.LT.1)THEN
                IQTB=1
                PP1=0.
              ELSEIF(IQTB.GE.ITBL)THEN
                IQTB=ITBL-1
                PP1=0.
              ENDIF
!--------------SATURATION PRESSURE AT FOUR SURROUNDING TABLE PTS.-------
              P00K=PTBL(IQTB  ,ITTB  )
              P10K=PTBL(IQTB+1,ITTB  )
              P01K=PTBL(IQTB  ,ITTB+1)
              P11K=PTBL(IQTB+1,ITTB+1)
!--------------SATURATION POINT VARIABLES AT THE BOTTOM-----------------
              PSP=P00K+(P10K-P00K)*PP1+(P01K-P00K)*QQ1 &
                 +(P00K-P10K-P01K+P11K)*PP1*QQ1
              RXNERS=(1.E5/PSP)**CAPPA
              THEK(K)=THBT*EXP(ELOCP*QBT*RXNERS/THBT)
              PSK (K)=PSP
!-----------------------------------------------------------------------
            ELSE
!-----------------------------------------------------------------------
              THEK(K)=THEK(K+1)
              PSK (K)=PINH(I,J,1)
!-----------------------------------------------------------------------
            ENDIF
!-----------------------------------------------------------------------
          ENDDO
!
!***  Begin debugging
!         IF(I==IMD.AND.J==JMD)THEN
!           PRINT_DIAG=1
!         ELSE
!           PRINT_DIAG=0
!         ENDIF
!         IF(I==227.AND.J==363)PRINT_DIAG=2
!***  End debugging
!
!----------------------------------------------------------------------
!***
!***  FIND THE MIXING LENGTH
!***
          CALL MIXLEN(LMH,RSQDT,UK,VK,THVK,THEK &
                     ,Q2K,EPSL,EPSQ2,ZHK,PK,PSK,RXNERK,GM,GH,EL &
                     ,PBLH(I,J),LPBL(I,J),LMXL,CT(I,J),MIXHT(I,J) &
                     ,I,J,LM)
!
!----------------------------------------------------------------------
!***
!***  SOLVE FOR THE PRODUCTION/DISSIPATION OF
!***  THE TURBULENT KINETIC ENERGY
!***
!
          CALL PRODQ2(LMH,DTTURBL,USTAR(I,J),GM,GH,EL,Q2K  &
                     ,EPSL,EPSQ2,I,J,LM)
!
!----------------------------------------------------------------------
!*** THE MODEL LAYER (COUNTING UPWARD) CONTAINING THE TOP OF THE PBL
!----------------------------------------------------------------------
!
          KPBL(I,J)=LPBL(I,J)
!
!----------------------------------------------------------------------
!***
!***  FIND THE EXCHANGE COEFFICIENTS IN THE FREE ATMOSPHERE
!***
          CALL DIFCOF(LMH,LMXL,GM,GH,EL,TK,Q2K,ZHK,AKMK,AKHK,I,J,LM &
                     ,PRINT_DIAG)
!
!***  COUNTING DOWNWARD FROM THE TOP, THE EXCHANGE COEFFICIENTS AKH
!***  ARE DEFINED ON THE BOTTOMS OF THE LAYERS 1 TO LM-1.  COUNTING
!***  COUNTING UPWARD FROM THE BOTTOM, THOSE SAME COEFFICIENTS EXCH_H
!***  ARE DEFINED ON THE TOPS OF THE LAYERS 1 TO LM-1.
!
          DO K=1,LM-1
            AKH(I,J,K)=AKHK(K)
            AKM(I,J,K)=AKMK(K)
            DELTAZ=0.5*(ZHK(K)-ZHK(K+2))
!           EXCH_H(I,J,K)=AKHK(K)*DELTAZ
          ENDDO
!
!----------------------------------------------------------------------
!***
!***  CARRY OUT THE VERTICAL DIFFUSION OF
!***  TURBULENT KINETIC ENERGY
!***
!
          CALL VDIFQ(LMH,DTDIF,Q2K,EL,ZHK,I,J,LM)
!
!***  SAVE THE NEW Q2 AND MIXING LENGTH.
!
          DO K=1,LM
            Q2(I,J,K)=MAX(Q2K(K),EPSQ2(K))
            IF(K<LM)EL_MYJ(I,J,K)=EL(K)   ! EL IS NOT DEFINED AT LM
          ENDDO
!
        ENDDO
!
!----------------------------------------------------------------------
!
      ENDDO setup_integration
!
!.......................................................................
!ZJ$OMP END PARALLEL DO
!.......................................................................
!----------------------------------------------------------------------
!
!***  CONVERT SURFACE SENSIBLE TEMPERATURE TO POTENTIAL TEMPERATURE.
!
      DO J=JMS,JME
      DO I=IMS,IME
        PSFC=PINH(I,J,LM+1)
        THSK(I,J)=TSK(I,J)*(1.E5/PSFC)**CAPPA
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!.......................................................................
!ZJ$OMP PARALLEL DO  PRIVATE(I,J,K &
!ZJ$OMP  & ,THK,QK,CWMK,ZHK,RHOK &
!ZJ$OMP  & ,AKHK,SEAMASK,LLOW,AKHS_DENS,QFC1,QLOW,PSFC,RXNSFC,LMH &
!ZJ$OMP  & ,THOLD,THNEW,DTDT,DQDT,DCDT,ZSL_DIAG,AKMK,AKMS_DENS &
!ZJ$OMP  & ,UK,VK &
!ZJ$OMP  & ),SCHEDULE(DYNAMIC)
!.......................................................................
!----------------------------------------------------------------------
      main_integration:  DO J=JMS,JME
!----------------------------------------------------------------------
!
        DO I=IMS,IME
!
!***  FILL 1-D VERTICAL ARRAYS
!
          DO K=LM,1,-1
            THK(K)=TH(I,J,K)
            QK(K)=Q(I,J,K)
            CWMK(K)=CWM(I,J,K)
            ZHK(K)=ZINT(I,J,K)
            RHOK(K)=PMID(I,J,K)/(R_D*T(I,J,K)*(1.+P608*QK(K)-CWMK(K)))
          ENDDO
!
!***  COUNTING DOWNWARD FROM THE TOP, THE EXCHANGE COEFFICIENTS AKH
!***  ARE DEFINED ON THE BOTTOMS OF THE LAYERS 1 TO LM-1.  THESE COEFFICIENTS
!***  ARE ALSO MULTIPLIED BY THE DENSITY AT THE BOTTOM INTERFACE LEVEL.
!
          DO K=1,LM-1
            AKHK(K)=AKH(I,J,K)*0.5*(RHOK(K)+RHOK(K+1))
          ENDDO
!
          ZHK(LM+1)=ZINT(I,J,LM+1)
!
          SEAMASK=XLAND(I,J)-1.
          THZ0(I,J)=(1.-SEAMASK)*THSK(I,J)+SEAMASK*THZ0(I,J)
!
          LLOW=LM
          AKHS_DENS=AKHS(I,J)*RHOK(LM)
!
          IF(SEAMASK<0.5)THEN
            QFC1=XLV*CHKLOWQ(I,J)*AKHS_DENS
!
            IF(SNOW(I,J)>0..OR.SICE(I,J)>0.5)THEN
              QFC1=QFC1*RLIVWV
            ENDIF
!
            IF(QFC1>0.)THEN
              QLOW=QK(LM)
              QSFC(I,J)=QLOW+ELFLX(I,J)/QFC1
            ENDIF
!
          ELSE
            PSFC=PINH(I,J,LM+1)
            RXNSFC=(1.E5/PSFC)**CAPPA

            QSFC(I,J)=PQ0SEA/PSFC                                      &
     &         *EXP(A2*(THSK(I,J)-A3*RXNSFC)/(THSK(I,J)-A4*RXNSFC))
          ENDIF
!
          QZ0 (I,J)=(1.-SEAMASK)*QSFC(I,J)+SEAMASK*QZ0 (I,J)
!
          LMH=LM
!
!----------------------------------------------------------------------
!***  CARRY OUT THE VERTICAL DIFFUSION OF
!***  TEMPERATURE AND WATER VAPOR
!----------------------------------------------------------------------
!
          CALL VDIFH(DTDIF,LMH,THZ0(I,J),QZ0(I,J) &
                    ,AKHS_DENS,CHKLOWQ(I,J),CT(I,J) &
                    ,THK,QK,CWMK,AKHK,ZHK,RHOK,I,J,LM)
!----------------------------------------------------------------------
!***
!***  COMPUTE PRIMARY VARIABLE TENDENCIES
!***
          DO K=1,LM
            RTHBLTEN(I,J,K)=(THK(K)-TH(I,J,K))*RDTTURBL
            RQBLTEN(I,J,K)=(QK(K)-Q(I,J,K))*RDTTURBL
            RQCBLTEN(I,J,K)=(CWMK(K)-CWM(I,J,K))*RDTTURBL
          ENDDO
!
!*** Begin debugging
!         IF(I==IMD.AND.J==JMD)THEN
!           PRINT_DIAG=0
!         ELSE
!           PRINT_DIAG=0
!         ENDIF
!         IF(I==227.AND.J==363)PRINT_DIAG=0
!*** End debugging
!
        PSFC=.01*PINH(I,J,LM+1)
        ZSL_DIAG=0.5*DZ(I,J,LM)
!
!*** Begin debugging
!         IF(PRINT_DIAG==1)THEN
!
!           WRITE(6,"(A, 2I5, 2I3, 2F8.2, F6.2, 2F8.2)") &
!           '{TURB4 I,J, KPBL, KMXL, PSFC, ZSFC, ZSL, ZPBL, ZMXL = ' &
!           , I, J, KPBL(I,J), LM-LMXL+1, PSFC, ZHK(LMH+1), ZSL_DIAG  &
!           , PBLH(I,J), ZHK(LMXL)-ZHK(LMH+1)
!           WRITE(6,"(A, 2F7.2, F7.3, 3E11.4)") &
!           '{TURB4 TSK, THSK, QZ0, Q**2_0, AKHS, EXCH_0 = ' &
!           , TSK(I,J)-273.15, THSK(I,J), 1000.*QZ0(I,J) &
!           , Q2(I,1,J), AKHS(I,J), AKHS(I,J)*ZSL_DIAG
!           WRITE(6,"(A)") &
!           '{TURB5 K, PMID, PINH_1, TC, TH, DTH, GH, GM, EL, Q**2, AKH, EXCH_H, DZ, DP'
!           DO K=1,LM/2
!             WRITE(6,"(A,I3, 2F8.2, 2F8.3, 3E12.4, 4E11.4, F7.2, F6.2)") &
!            '{TURB5 ', K, .01*PMID(I,K,J),.01*PINH(I,K,J), T(I,K,J)-273.15 &
!            , TH(I,K,J), DTTURBL*RTHBLTEN(I,K,J), GH(K), GM(K) &
!            , EL_MYJ(I,K,J), Q2(I,K+1,J), AKH(I,K,J) &
!            , EXCH_H(I,K,J), DZ(I,K,J), .01*(PINH(I,K,J)-PINH(I,K+1,J))
!           ENDDO
!
!         ELSEIF(PRINT_DIAG==2)THEN
!
!           WRITE(6,"(A, 2I5, 2I3, 2F8.2, F6.2, 2F8.2)") &
!           '}TURB4 I,J, KPBL, KMXL, PSFC, ZSFC, ZSL, ZPBL, ZMXL = ' &
!           , I, J, KPBL(I,J), LM-LMXL+1, PSFC, ZHK(LMH+1), ZSL_DIAG  &
!           , PBLH(I,J), ZHK(LMXL)-ZHK(LMH+1)
!           WRITE(6,"(A, 2F7.2, F7.3, 3E11.4)") &
!           '}TURB4 TSK, THSK, QZ0, Q**2_0, AKHS, EXCH_0 = ' &
!           , TSK(I,J)-273.15, THSK(I,J), 1000.*QZ0(I,J) &
!           , Q2(I,1,J), AKHS(I,J), AKHS(I,J)*ZSL_DIAG
!           WRITE(6,"(A)") &
!           '}TURB5 K, PMID, PINH_1, TC, TH, DTH, GH, GM, EL, Q**2, AKH, EXCH_H, DZ, DP'
!           DO K=1,LM/2
!             WRITE(6,"(A,I3, 2F8.2, 2F8.3, 3E12.4, 4E11.4, F7.2, F6.2)") &
!            '}TURB5 ', K, .01*PMID(I,K,J),.01*PINH(I,K,J), T(I,K,J)-273.15 &
!            , TH(I,K,J), DTTURBL*RTHBLTEN(I,K,J), GH(K), GM(K) &
!            , EL_MYJ(I,K,J), Q2(I,K+1,J), AKH(I,K,J) &
!            , EXCH_H(I,K,J), DZ(I,K,J), .01*(PINH(I,K,J)-PINH(I,K+1,J))
!           ENDDO
!         ENDIF
!*** End debugging
!
!----------------------------------------------------------------------
!
          SEAMASK=XLAND(I,J)-1.
!
          IF(SEAMASK.LT.0.5.AND.STDH(I,J).GT.1.) THEN
            RSTDH=1./STDH(I,J)
          ELSE
            RSTDH=0.
          ENDIF
          ZHK(LM+1)=ZINT(I,J,LM+1)
          ZSF=STDH(I,J)*STDFAC+ZHK(LM+1)
!
!----------------------------------------------------------------------
!
!***  FILL 1-D VERTICAL ARRAYS
!
          DO K=1,LM-1
            AKMK(K)=AKM(I,J,K)
            AKMK(K)=AKMK(K)*(RHOK(K)+RHOK(K+1))*0.5
          ENDDO
!
          AKMS_DENS=AKMS(I,J)*RHOK(LM)
!
          DO K=LM,1,-1
            UK(K)=U(I,J,K)
            VK(K)=V(I,J,K)
            ZHK(K)=ZINT(I,J,K)
          ENDDO
          ZHK(LM+1)=ZINT(I,J,LM+1)
!
!----------------------------------------------------------------------
!
          DO K=1,LM-1
!jun23            IF(SEAMASK.GT.0.5) THEN
!jun23              DCOL(K)=0.
!jun23            ELSE
!jun23              ZUV=(ZHK(K)+ZHK(K+1))*0.5
!jun23              IF(ZUV.GT.ZSF) THEN
!jun23                DCOL(K)=0.
!jun23              ELSE
!jun23                DCOL(K)=HERF((((ZUV-ZHK(LM+1))*RSTDH)**2)*0.5)
!jun23              ENDIF
!jun23            ENDIF
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
            DCOL(K)=0. !ZJ
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
          ENDDO
!
!----------------------------------------------------------------------
!***  CARRY OUT THE VERTICAL DIFFUSION OF
!***  VELOCITY COMPONENTS
!----------------------------------------------------------------------
!
          CALL VDIFV(LMH,DTDIF,UZ0(I,J),VZ0(I,J) &
     &              ,AKMS_DENS,DCOL,UK,VK,AKMK,ZHK,RHOK,I,J,LM)
!
!----------------------------------------------------------------------
!***
!***  COMPUTE PRIMARY VARIABLE TENDENCIES
!***
          DO K=1,LM
            RUBLTEN(I,J,K)=(UK(K)-U(I,J,K))*RDTTURBL
            RVBLTEN(I,J,K)=(VK(K)-V(I,J,K))*RDTTURBL
          ENDDO
!
        ENDDO
!----------------------------------------------------------------------
!
      ENDDO main_integration
!JAA!ZJ$OMP END PARALLEL DO
!
!----------------------------------------------------------------------
!
      END SUBROUTINE MYJPBL
!
!----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------
                          SUBROUTINE MIXLEN                            &
!----------------------------------------------------------------------
!   ******************************************************************
!   *                                                                *
!   *                   LEVEL 2.5 MIXING LENGTH                      *
!   *                                                                *
!   ******************************************************************
!
      (LMH,RSQDT,U,V,THV,THE,Q2,EPSL,EPSQ2,Z,P,PS,RXNER &
      ,GM,GH,EL,PBLH,LPBL,LMXL,CT,MIXHT,I,J,LM)
!
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER(KIND=KINT),INTENT(IN):: &
       LMH,I,J,LM
!
      REAL(KIND=KFPT),INTENT(IN):: &
       RSQDT
!
      INTEGER(KIND=KINT),INTENT(OUT):: &
       LMXL,LPBL
!
      real(kind=kfpt),dimension(1:lm-1),intent(in):: EPSL
      REAL(KIND=KFPT),DIMENSION(1:LM),INTENT(IN):: &
       P,PS,Q2,EPSQ2,RXNER,THE,THV,U,V
!
      REAL(KIND=KFPT),DIMENSION(1:LM+1),INTENT(IN):: &
       Z
!
      REAL(KIND=KFPT),INTENT(OUT):: &
       MIXHT &
      ,PBLH
!
      REAL(KIND=KFPT),DIMENSION(1:LM-1),INTENT(OUT):: &
       EL,GH,GM
!
      REAL(KIND=KFPT),INTENT(INOUT):: CT
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER(KIND=KINT):: &
       K,LPBLM
!
      REAL(KIND=KFPT):: &
       ADEN,BDEN,AUBR,BUBR,BLMX,CUBRY,DTHV,DZ &
      ,EL0,ELOQ2X,GHL,GML &
      ,QOL2ST,QOL2UN,QDZL &
      ,RDZ,SQ,SREL,SZQ,VKRMZ,WCON
!
      REAL(KIND=KFPT),DIMENSION(1:LM):: &
       Q1
!
      REAL(KIND=KFPT),DIMENSION(1:LM-1):: &
       ELM,REL
!
!----------------------------------------------------------------------
!***********************************************************************
!--------1---------2---------3---------4---------5---------6---------7--
      CUBRY=UBRY*1.5  !*2.
!--------------FIND THE HEIGHT OF THE PBL-------------------------------
      LPBL=LMH
!
      DO K=LMH-1,1,-1
        if(q2(k)-epsq2(k)+epsq2(lm).le.epsq2(lm)*fh) then
          LPBL=K
          GO TO 110
        ENDIF
      ENDDO
!
      LPBL=1
!
!--------------THE HEIGHT OF THE PBL------------------------------------
!
 110  PBLH=Z(LPBL+1)-Z(LMH+1)
!
!-----------------------------------------------------------------------
      DO K=1,LMH
        Q1(K)=0.
      ENDDO
!-----------------------------------------------------------------------
      DO K=1,LMH-1
        DZ=(Z(K)-Z(K+2))*0.5
        RDZ=1./DZ
        GML=((U(K)-U(K+1))**2+(V(K)-V(K+1))**2)*RDZ*RDZ
        GM(K)=MAX(GML,EPSGM)
!
        DTHV=THV(K)-THV(K+1)
!----------------------------------------------------------------------
        IF(DTHV.GT.0.) THEN
          IF(THE(K+1).GT.THE(K)) THEN
            IF(PS(K+1).GT.P(K)) THEN                                       !>12KM
!
              WCON=(P(K+1)-PS(K+1))/(P(K+1)-P(K))
!
              if( &
                 (q2(k).gt.epsq2(k)) .and. &
                  (q2(k)*cubry.gt.(dz*wcon*rsqdt)**2) &
                ) then 
!
                 DTHV=(THE(K)-THE(K+1))+DTHV
!
              ENDIF
            ENDIF
          ENDIF
        ENDIF
!--------------------------------------------------------------------------
!
        GHL=DTHV*RDZ
        IF(ABS(GHL)<=EPSGH)GHL=EPSGH
        GH(K)=GHL
      ENDDO
!
      CT=0.
!
!----------------------------------------------------------------------
!***  FIND MAXIMUM MIXING LENGTHS AND THE LEVEL OF THE PBL TOP
!----------------------------------------------------------------------
!
      LMXL=LMH
!
      DO K=1,LMH-1
        GML=GM(K)
        GHL=GH(K)
!
        IF(GHL>=EPSGH)THEN
          IF(GML/GHL<=REQU)THEN
            ELM(K)=EPSL(K)
            LMXL=K+1
          ELSE
            AUBR=(AUBM*GML+AUBH*GHL)*GHL
            BUBR= BUBM*GML+BUBH*GHL
            QOL2ST=(-0.5*BUBR+SQRT(BUBR*BUBR*0.25-AUBR*CUBR))*RCUBR
            ELOQ2X=1./MAX(EPSGH, QOL2ST)
            ELM(K)=MAX(SQRT(ELOQ2X*Q2(K)),EPSL(K))
          ENDIF
        ELSE
          ADEN=(ADNM*GML+ADNH*GHL)*GHL
          BDEN= BDNM*GML+BDNH*GHL
          QOL2UN=-0.5*BDEN+SQRT(BDEN*BDEN*0.25-ADEN)
          ELOQ2X=1./(QOL2UN+EPSRU)       ! REPSR1/QOL2UN
          ELM(K)=MAX(SQRT(ELOQ2X*Q2(K)),EPSL(K))
        ENDIF
      ENDDO
!
      IF(ELM(LMH-1)==EPSL(LMH-1))LMXL=LMH
!
!----------------------------------------------------------------------
!***  THE HEIGHT OF THE MIXED LAYER
!----------------------------------------------------------------------
!
      BLMX=Z(LMXL)-Z(LMH+1)
      MIXHT=BLMX
!
!----------------------------------------------------------------------
      DO K=LPBL,LMH
        Q1(K)=SQRT(Q2(K))
      ENDDO
!----------------------------------------------------------------------
      SZQ=0.
      SQ =0.
!
      DO K=1,LMH-1
        QDZL=(Q1(K)+Q1(K+1))*(Z(K+1)-Z(K+2))
        SZQ=(Z(K+1)+Z(K+2)-Z(LMH+1)-Z(LMH+1))*QDZL+SZQ
        SQ=QDZL+SQ
      ENDDO
!
!----------------------------------------------------------------------
!***  COMPUTATION OF ASYMPTOTIC L IN BLACKADAR FORMULA
!----------------------------------------------------------------------
!
      EL0=MIN(ALPH*SZQ*0.5/SQ,EL0MAX)
      EL0=MAX(EL0            ,EL0MIN)
!
!----------------------------------------------------------------------
!***  ABOVE THE PBL TOP
!----------------------------------------------------------------------
!
      LPBLM=MAX(LPBL-1,1)
!
      DO K=1,LPBLM
        EL(K)=MIN((Z(K)-Z(K+2))*ELFC,ELM(K))
        REL(K)=EL(K)/ELM(K)
      ENDDO
!
!----------------------------------------------------------------------
!***  INSIDE THE PBL
!----------------------------------------------------------------------
!
      IF(LPBL<LMH)THEN
        DO K=LPBL,LMH-1
          VKRMZ=(Z(K+1)-Z(LMH+1))*VKARMAN
          EL(K)=MIN(VKRMZ/(VKRMZ/EL0+1.),ELM(K))
          REL(K)=EL(K)/ELM(K)
        ENDDO
      ENDIF
!
      DO K=LPBL+1,LMH-2
        SREL=MIN(((REL(K-1)+REL(K+1))*0.5+REL(K))*0.5,REL(K))
        EL(K)=MAX(SREL*ELM(K),EPSL(K))
      ENDDO
!
!----------------------------------------------------------------------
      END SUBROUTINE MIXLEN
!----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------
                          SUBROUTINE PRODQ2                            &
!----------------------------------------------------------------------
!   ******************************************************************
!   *                                                                *
!   *            LEVEL 2.5 Q2 PRODUCTION/DISSIPATION                 *
!   *                                                                *
!   ******************************************************************
!
      (LMH,DTTURBL,USTAR,GM,GH,EL,Q2,EPSL,EPSQ2,I,J,LM)
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER(KIND=KINT),INTENT(IN):: &
       LMH,I,J,LM
!
      REAL(KIND=KFPT),INTENT(IN):: &
       DTTURBL,USTAR
!
      REAL(KIND=KFPT),DIMENSION(1:LM-1),INTENT(IN):: &
       GH,GM
!
      real(kind=kfpt),dimension(1:lm-1),intent(in):: EPSL
      real(kind=kfpt),dimension(1:lm),intent(in):: EPSQ2
!
      REAL(KIND=KFPT),DIMENSION(1:LM-1),INTENT(INOUT):: &
       EL
!
      REAL(KIND=KFPT),DIMENSION(1:LM),INTENT(INOUT):: &
       Q2
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER(KIND=KINT):: &
       K
!
      REAL(KIND=KFPT):: &
       ADEN,AEQU,ANUM,ARHS,BDEN,BEQU,BNUM,BRHS,CDEN,CRHS &
      ,DLOQ1,ELOQ11,ELOQ12,ELOQ13,ELOQ21,ELOQ22,ELOQ31,ELOQ32 &
      ,ELOQ41,ELOQ42,ELOQ51,ELOQ52,ELOQN,EQOL2,GHL,GML &
      ,RDEN1,RDEN2,RHS2,RHSP1,RHSP2,RHST2
!
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!
      main_integration: DO K=1,LMH-1
        GML=GM(K)
        GHL=GH(K)
!
!----------------------------------------------------------------------
!***  COEFFICIENTS OF THE EQUILIBRIUM EQUATION
!----------------------------------------------------------------------
!
        AEQU=(AEQM*GML+AEQH*GHL)*GHL
        BEQU= BEQM*GML+BEQH*GHL
!
!----------------------------------------------------------------------
!***  EQUILIBRIUM SOLUTION FOR L/Q
!----------------------------------------------------------------------
!
        EQOL2=-0.5*BEQU+SQRT(BEQU*BEQU*0.25-AEQU)
!
!----------------------------------------------------------------------
!***  IS THERE PRODUCTION/DISSIPATION ?
!----------------------------------------------------------------------
!
        IF((GML+GHL*GHL<=EPSTRB)                                       &
     &   .OR.(GHL>=EPSGH.AND.GML/GHL<=REQU)                            &
     &   .OR.(EQOL2<=EPS2))THEN
!
!----------------------------------------------------------------------
!***  NO TURBULENCE
!----------------------------------------------------------------------
!
          Q2(K)=EPSQ2(K)
          EL(K)=EPSL(K)
!----------------------------------------------------------------------
!
        ELSE
!
!----------------------------------------------------------------------
!***  TURBULENCE
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!***  COEFFICIENTS OF THE TERMS IN THE NUMERATOR
!----------------------------------------------------------------------
!
          ANUM=(ANMM*GML+ANMH*GHL)*GHL
          BNUM= BNMM*GML+BNMH*GHL
!
!----------------------------------------------------------------------
!***  COEFFICIENTS OF THE TERMS IN THE DENOMINATOR
!----------------------------------------------------------------------
!
          ADEN=(ADNM*GML+ADNH*GHL)*GHL
          BDEN= BDNM*GML+BDNH*GHL
          CDEN= 1.
!
!----------------------------------------------------------------------
!***  COEFFICIENTS OF THE NUMERATOR OF THE LINEARIZED EQ.
!----------------------------------------------------------------------
!
          ARHS=-(ANUM*BDEN-BNUM*ADEN)*2.
          BRHS=- ANUM*4.
          CRHS=- BNUM*2.
!
!----------------------------------------------------------------------
!***  INITIAL VALUE OF L/Q
!----------------------------------------------------------------------
!
          DLOQ1=EL(K)/SQRT(Q2(K))
!
!----------------------------------------------------------------------
!***  FIRST ITERATION FOR L/Q, RHS=0
!----------------------------------------------------------------------
!
          ELOQ21=1./EQOL2
          ELOQ11=SQRT(ELOQ21)
          ELOQ31=ELOQ21*ELOQ11
          ELOQ41=ELOQ21*ELOQ21
          ELOQ51=ELOQ21*ELOQ31
!
!----------------------------------------------------------------------
!***  1./DENOMINATOR
!----------------------------------------------------------------------
!
          RDEN1=1./(ADEN*ELOQ41+BDEN*ELOQ21+CDEN)
!
!----------------------------------------------------------------------
!***  D(RHS)/D(L/Q)
!----------------------------------------------------------------------
!
          RHSP1=(ARHS*ELOQ51+BRHS*ELOQ31+CRHS*ELOQ11)*RDEN1*RDEN1
!
!----------------------------------------------------------------------
!***  FIRST-GUESS SOLUTION
!----------------------------------------------------------------------
!
          ELOQ12=ELOQ11+(DLOQ1-ELOQ11)*EXP(RHSP1*DTTURBL)
          ELOQ12=MAX(ELOQ12,EPS1)
!
!----------------------------------------------------------------------
!***  SECOND ITERATION FOR L/Q
!----------------------------------------------------------------------
!
          ELOQ22=ELOQ12*ELOQ12
          ELOQ32=ELOQ22*ELOQ12
          ELOQ42=ELOQ22*ELOQ22
          ELOQ52=ELOQ22*ELOQ32
!
!----------------------------------------------------------------------
!***  1./DENOMINATOR
!----------------------------------------------------------------------
!
          RDEN2=1./(ADEN*ELOQ42+BDEN*ELOQ22+CDEN)
          RHS2 =-(ANUM*ELOQ42+BNUM*ELOQ22)*RDEN2+RB1
          RHSP2= (ARHS*ELOQ52+BRHS*ELOQ32+CRHS*ELOQ12)*RDEN2*RDEN2
          RHST2=RHS2/RHSP2
!
!----------------------------------------------------------------------
!***  CORRECTED SOLUTION
!----------------------------------------------------------------------
!
          ELOQ13=ELOQ12-RHST2+(RHST2+DLOQ1-ELOQ12)*EXP(RHSP2*DTTURBL)
          ELOQ13=AMAX1(ELOQ13,EPS1)
!
!----------------------------------------------------------------------
!***  TWO ITERATIONS IS ENOUGH IN MOST CASES ...
!----------------------------------------------------------------------
!
          ELOQN=ELOQ13
!
          IF(ELOQN>EPS1)THEN
            Q2(K)=EL(K)*EL(K)/(ELOQN*ELOQN)
            Q2(K)=AMAX1(Q2(K),EPSQ2(K))
            IF(Q2(K)==EPSQ2(K))THEN
              EL(K)=EPSL(K)
            ENDIF
          ELSE
            Q2(K)=EPSQ2(K)
            EL(K)=EPSL(K)
          ENDIF
!
!----------------------------------------------------------------------
!***  END OF TURBULENT BRANCH
!----------------------------------------------------------------------
!
        ENDIF
!----------------------------------------------------------------------
!***  END OF PRODUCTION/DISSIPATION LOOP
!----------------------------------------------------------------------
!
      ENDDO main_integration
!
!----------------------------------------------------------------------
!***  LOWER BOUNDARY CONDITION FOR Q2
!----------------------------------------------------------------------
!
      Q2(LMH)=AMAX1(B1**(2./3.)*USTAR*USTAR,EPSQ2(LMH))
!----------------------------------------------------------------------
!
      END SUBROUTINE PRODQ2
!
!----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------
                           SUBROUTINE DIFCOF                           &
!   ******************************************************************
!   *                                                                *
!   *                LEVEL 2.5 DIFFUSION COEFFICIENTS                *
!   *                                                                *
!   ******************************************************************
      (LMH,LMXL,GM,GH,EL,T,Q2,Z,AKM,AKH,I,J,LM,PRINT_DIAG)
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER(KIND=KINT),INTENT(IN):: &
       LMH,LMXL,I,J,LM
!
      REAL(KIND=KFPT),DIMENSION(1:LM),INTENT(IN):: &
       Q2,T
!
      REAL(KIND=KFPT),DIMENSION(1:LM-1),INTENT(IN):: &
       EL,GH,GM
!
      REAL(KIND=KFPT),DIMENSION(1:LM+1),INTENT(IN):: &
       Z
!
      REAL(KIND=KFPT),DIMENSION(1:LM-1),INTENT(OUT):: &
       AKH,AKM
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER(KIND=KINT):: &
       K,KINV
!
      REAL(KIND=KFPT):: &
       ADEN,AKMIN,BDEN,BESH,BESM,CDEN,D2T,ELL,ELOQ2,ELOQ4,ELQDZ &
      ,ESH,ESM,GHL,GML,Q1L,RDEN,RDZ
!
!*** Begin debugging
      INTEGER(KIND=KINT),INTENT(IN):: PRINT_DIAG
!     REAL(KIND=KFPT):: D2TMIN
!*** End debugging
!
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!
      DO K=1,LMH-1
        ELL=EL(K)
!
        ELOQ2=ELL*ELL/Q2(K)
        ELOQ4=ELOQ2*ELOQ2
!
        GML=GM(K)
        GHL=GH(K)
!
!----------------------------------------------------------------------
!***  COEFFICIENTS OF THE TERMS IN THE DENOMINATOR
!----------------------------------------------------------------------
!
        ADEN=(ADNM*GML+ADNH*GHL)*GHL
        BDEN= BDNM*GML+BDNH*GHL
        CDEN= 1.
!
!----------------------------------------------------------------------
!***  COEFFICIENTS FOR THE SM DETERMINANT
!----------------------------------------------------------------------
!
        BESM=BSMH*GHL
!
!----------------------------------------------------------------------
!***  COEFFICIENTS FOR THE SH DETERMINANT
!----------------------------------------------------------------------
!
        BESH=BSHM*GML+BSHH*GHL
!
!----------------------------------------------------------------------
!***  1./DENOMINATOR
!----------------------------------------------------------------------
!
        RDEN=1./(ADEN*ELOQ4+BDEN*ELOQ2+CDEN)
!
!----------------------------------------------------------------------
!***  SM AND SH
!----------------------------------------------------------------------
!
        ESM=(BESM*ELOQ2+CESM)*RDEN
        ESH=(BESH*ELOQ2+CESH)*RDEN
!
!----------------------------------------------------------------------
!***  DIFFUSION COEFFICIENTS
!----------------------------------------------------------------------
!
        RDZ=2./(Z(K)-Z(K+2))
        Q1L=SQRT(Q2(K))
        ELQDZ=ELL*Q1L*RDZ
        AKM(K)=ELQDZ*ESM
        AKH(K)=ELQDZ*ESH
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!        AKM(K)=MAX(AKM(K),RDZ*3.)
!        AKH(K)=MAX(AKH(K),RDZ*3.)
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!----------------------------------------------------------------------
      ENDDO
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!***  INVERSIONS
!----------------------------------------------------------------------
!
!     IF(LMXL==LMH)THEN
!       KINV=LMH
!       D2TMIN=0.
!
!       DO K=LMH/2,LMH-1
!         D2T=T(K-1)-2.*T(K)+T(K+1)
!         IF(D2T<D2TMIN)THEN
!           D2TMIN=D2T
!           IF(D2T<0)KINV=K
!         ENDIF
!       ENDDO
!
!       IF(KINV<LMH)THEN
!         DO K=KINV-1,LMH-1
!           RDZ=2./(Z(K)-Z(K+2))
!           AKMIN=0.5*RDZ
!           AKM(K)=MAX(AKM(K),AKMIN)
!           AKH(K)=MAX(AKH(K),AKMIN)
!         ENDDO
!
!*** Begin debugging
!         IF(PRINT_DIAG>0)THEN
!           WRITE(6,"(A,3I3)") '{TURB1 LMXL,LMH,KINV=',LMXL,LMH,KINV
!           WRITE(6,"(A,3I3)") '}TURB1 LMXL,LMH,KINV=',LMXL,LMH,KINV
!           IF(PRINT_DIAG==1)THEN
!             WRITE(6,"(A)") &
!               '{TURB3 K, T, D2T, RDZ, Z(K), Z(K+2), AKMIN, AKH '
!           ELSE
!             WRITE(6,"(A)") &
!               '}TURB3 K, T, D2T, RDZ, Z(K), Z(K+2), AKMIN, AKH '
!           ENDIF
!           DO K=LMH-1,KINV-1,-1
!             D2T=T(K-1)-2.*T(K)+T(K+1)
!             RDZ=2./(Z(K)-Z(K+2))
!             AKMIN=0.5*RDZ
!             IF(PRINT_DIAG==1)THEN
!               WRITE(6,"(A,I3,F8.3,2E12.5,2F9.2,2E12.5)") '{TURB3 ' &
!               ,K,T(K)-273.15,D2T,RDZ,Z(K),Z(K+2),AKMIN,AKH(K)
!             ELSE
!               WRITE(6,"(A,I3,F8.3,2E12.5,2F9.2,2E12.5)") '}TURB3 ' &
!               ,K,T(K)-273.15,D2T,RDZ,Z(K),Z(K+2),AKMIN,AKH(K)
!             ENDIF
!           ENDDO
!         ENDIF     !- IF (PRINT_DIAG > 0) THEN
!       ENDIF       !- IF(KINV<LMH)THEN
!*** End debugging
!
!     ENDIF
!----------------------------------------------------------------------
!
      END SUBROUTINE DIFCOF
!
!----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------
                           SUBROUTINE VDIFQ                            &
!   ******************************************************************
!   *                                                                *
!   *               VERTICAL DIFFUSION OF Q2 (TKE)                   *
!   *                                                                *
!   ******************************************************************
      (LMH,DTDIF,Q2,EL,Z,I,J,LM)
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER(KIND=KINT),INTENT(IN):: &
       LMH,I,J,LM
!
      REAL(KIND=KFPT),INTENT(IN):: &
       DTDIF
!
      REAL(KIND=KFPT),DIMENSION(1:LM-1),INTENT(IN):: &
       EL
!
      REAL(KIND=KFPT),DIMENSION(1:LM+1),INTENT(IN):: &
       Z
!
      REAL(KIND=KFPT),DIMENSION(1:LM),INTENT(INOUT):: &
       Q2
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER(KIND=KINT):: &
       K
!
      REAL(KIND=KFPT):: &
       ADEN,AKQS,BDEN,BESH,BESM,CDEN,CF,DTOZS,ELL,ELOQ2,ELOQ4 &
      ,ELQDZ,ESH,ESM,ESQHF,GHL,GML,Q1L,RDEN,RDZ
!
      REAL(KIND=KFPT),DIMENSION(1:LM-2):: &
       AKQ,CM,CR,DTOZ,RSQ2
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!***
!***  VERTICAL TURBULENT DIFFUSION
!***
!----------------------------------------------------------------------
      ESQHF=0.5*ESQ
!
      DO K=1,LMH-2
        DTOZ(K)=(DTDIF+DTDIF)/(Z(K)-Z(K+2))
        AKQ(K)=SQRT((Q2(K)+Q2(K+1))*0.5)*(EL(K)+EL(K+1))*ESQHF         &
     &        /(Z(K+1)-Z(K+2))
        CR(K)=-DTOZ(K)*AKQ(K)
      ENDDO
!
      CM(1)=DTOZ(1)*AKQ(1)+1.
      RSQ2(1)=Q2(1)
!
      DO K=1+1,LMH-2
        CF=-DTOZ(K)*AKQ(K-1)/CM(K-1)
        CM(K)=-CR(K-1)*CF+(AKQ(K-1)+AKQ(K))*DTOZ(K)+1.
        RSQ2(K)=-RSQ2(K-1)*CF+Q2(K)
      ENDDO
!
      DTOZS=(DTDIF+DTDIF)/(Z(LMH-1)-Z(LMH+1))
      AKQS=SQRT((Q2(LMH-1)+Q2(LMH))*0.5)*(EL(LMH-1)+ELZ0)*ESQHF        &
     &    /(Z(LMH)-Z(LMH+1))
!
      CF=-DTOZS*AKQ(LMH-2)/CM(LMH-2)
!
      Q2(LMH-1)=(DTOZS*AKQS*Q2(LMH)-RSQ2(LMH-2)*CF+Q2(LMH-1))          &
     &        /((AKQ(LMH-2)+AKQS)*DTOZS-CR(LMH-2)*CF+1.)
!
      DO K=LMH-2,1,-1
        Q2(K)=(-CR(K)*Q2(K+1)+RSQ2(K))/CM(K)
      ENDDO
!----------------------------------------------------------------------
!
      END SUBROUTINE VDIFQ
!
!----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!---------------------------------------------------------------------
      SUBROUTINE VDIFH(DTDIF,LMH,THZ0,QZ0,RKHS,CHKLOWQ,CT &
                      ,TH,Q,CWM,RKH,Z,RHO,I,J,LM)
!     ***************************************************************
!     *                                                             *
!     *         VERTICAL DIFFUSION OF MASS VARIABLES                *
!     *                                                             *
!     ***************************************************************
!---------------------------------------------------------------------
!
      IMPLICIT NONE
!
!---------------------------------------------------------------------
      INTEGER(KIND=KINT),INTENT(IN):: &
       LMH,I,J,LM
!
      REAL(KIND=KFPT),INTENT(IN):: &
       CHKLOWQ,CT,DTDIF,QZ0,RKHS,THZ0
!
      REAL(KIND=KFPT),DIMENSION(1:LM-1),INTENT(IN):: &
       RKH
!
      REAL(KIND=KFPT),DIMENSION(1:LM),INTENT(IN):: &
       RHO
!
      REAL(KIND=KFPT),DIMENSION(1:LM+1),INTENT(IN):: &
       Z
!
      REAL(KIND=KFPT),DIMENSION(1:LM),INTENT(INOUT):: &
       CWM,Q,TH
!
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER(KIND=KINT):: &
       K
!
      REAL(KIND=KFPT):: &
       CF,CMB,CMCB,CMQB,CMTB,CTHF,DTOZL,DTOZS &
      ,RCML,RKHH,RKQS,RSCB,RSQB,RSTB
!
      REAL(KIND=KFPT),DIMENSION(1:LM-1):: &
       CM,CR,DTOZ,RKCT,RSC,RSQ,RST
!
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
      CTHF=0.5*CT
!
      DO K=1,LMH-1
        DTOZ(K)=DTDIF/(Z(K)-Z(K+1))
        CR(K)=-DTOZ(K)*RKH(K)
        RKCT(K)=RKH(K)*(Z(K)-Z(K+2))*CTHF
      ENDDO
!
      CM(1)=DTOZ(1)*RKH(1)+RHO(1)
!----------------------------------------------------------------------
      RST(1)=TH (1)*RHO(1)-RKCT(1)*DTOZ(1)
      RSQ(1)=Q  (1)*RHO(1)
      RSC(1)=CWM(1)*RHO(1)
!----------------------------------------------------------------------
      DO K=1+1,LMH-1
        DTOZL=DTOZ(K)
        CF=-DTOZL*RKH(K-1)/CM(K-1)
        CM(K)=-CR(K-1)*CF+(RKH(K-1)+RKH(K))*DTOZL+RHO(K)
        RST(K)=-RST(K-1)*CF+(RKCT(K-1)-RKCT(K))*DTOZL+TH(K)*RHO(K)
        RSQ(K)=-RSQ(K-1)*CF+Q(K)  *RHO(K)
        RSC(K)=-RSC(K-1)*CF+CWM(K)*RHO(K)
      ENDDO
!
      DTOZS=DTDIF/(Z(LMH)-Z(LMH+1))
      RKHH=RKH(LMH-1)
!
      CF=-DTOZS*RKHH/CM(LMH-1)
      RKQS=RKHS*CHKLOWQ
!
      CMB=CR(LMH-1)*CF
      CMTB=-CMB+(RKHH+RKHS)*DTOZS+RHO(LMH)
      CMQB=-CMB+(RKHH+RKQS)*DTOZS+RHO(LMH)
      CMCB=-CMB+(RKHH     )*DTOZS+RHO(LMH)
!
      RSTB=-RST(LMH-1)*CF+RKCT(LMH-1)*DTOZS+TH(LMH)*RHO(LMH)
      RSQB=-RSQ(LMH-1)*CF+Q(LMH)  *RHO(LMH)
      RSCB=-RSC(LMH-1)*CF+CWM(LMH)*RHO(LMH)
!----------------------------------------------------------------------
      TH(LMH) =(DTOZS*RKHS*THZ0+RSTB)/CMTB
      Q(LMH)  =(DTOZS*RKQS*QZ0 +RSQB)/CMQB
      CWM(LMH)=(                RSCB)/CMCB
!----------------------------------------------------------------------
      DO K=LMH-1,1,-1
        RCML=1./CM(K)
        TH(K) =(-CR(K)* TH(K+1)+RST(K))*RCML
        Q(K)  =(-CR(K)*  Q(K+1)+RSQ(K))*RCML
        CWM(K)=(-CR(K)*CWM(K+1)+RSC(K))*RCML
      ENDDO
!----------------------------------------------------------------------
!
      END SUBROUTINE VDIFH
!
!---------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!---------------------------------------------------------------------
      SUBROUTINE VDIFV(LMH,DTDIF,UZ0,VZ0,RKMS,D,U,V,RKM,Z,RHO,I,J,LM)
!     ***************************************************************
!     *                                                             *
!     *        VERTICAL DIFFUSION OF VELOCITY COMPONENTS            *
!     *                                                             *
!     ***************************************************************
!---------------------------------------------------------------------
!
      IMPLICIT NONE
!
!---------------------------------------------------------------------
      INTEGER(KIND=KINT),INTENT(IN):: &
       LMH,I,J,LM
!
      REAL(KIND=KFPT),INTENT(IN):: &
       RKMS,DTDIF,UZ0,VZ0
!
      REAL(KIND=KFPT),DIMENSION(1:LM-1),INTENT(IN):: &
       D,RKM
!
      REAL(KIND=KFPT),DIMENSION(1:LM),INTENT(IN):: &
       RHO
!
      REAL(KIND=KFPT),DIMENSION(1:LM+1),INTENT(IN):: &
       Z
!
      REAL(KIND=KFPT),DIMENSION(1:LM),INTENT(INOUT):: &
       U,V
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER(KIND=KINT):: &
       K
!
      REAL(KIND=KFPT):: &
       CF,DTOZAK,DTOZL,DTOZS,RCML,RCMVB,RHOK,RKMH
!
      REAL(KIND=KFPT),DIMENSION(1:LM-1):: &
       CM,CR,DTOZ,RSU,RSV
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
      DO K=1,LMH-1
        DTOZ(K)=DTDIF/(Z(K)-Z(K+1))
        CR(K)=-DTOZ(K)*RKM(K)
      ENDDO
!
      RHOK=RHO(1)
      CM(1)=DTOZ(1)*RKM(1)+RHOK
      RSU(1)=U(1)*RHOK
      RSV(1)=V(1)*RHOK
!----------------------------------------------------------------------
      DO K=2,LMH-1
        DTOZL=DTOZ(K)
        CF=-DTOZL*RKM(K-1)/CM(K-1)
        RHOK=RHO(K)
        CM(K)=-CR(K-1)*CF+(RKM(K-1)+RKM(K)+D(K)*RKMS)*DTOZL+RHOK
        RSU(K)=-RSU(K-1)*CF+U(K)*RHOK
        RSV(K)=-RSV(K-1)*CF+V(K)*RHOK
      ENDDO
!----------------------------------------------------------------------
      DTOZS=DTDIF/(Z(LMH)-Z(LMH+1))
      RKMH=RKM(LMH-1)
!
      CF=-DTOZS*RKMH/CM(LMH-1)
      RHOK=RHO(LMH)
      RCMVB=1./((RKMH+RKMS)*DTOZS-CR(LMH-1)*CF+RHOK)
      DTOZAK=DTOZS*RKMS
!----------------------------------------------------------------------
      U(LMH)=(DTOZAK*UZ0-RSU(LMH-1)*CF+U(LMH)*RHOK)*RCMVB
      V(LMH)=(DTOZAK*VZ0-RSV(LMH-1)*CF+V(LMH)*RHOK)*RCMVB
!----------------------------------------------------------------------
      DO K=LMH-1,1,-1
        RCML=1./CM(K)
        U(K)=(-CR(K)*U(K+1)+RSU(K))*RCML
        V(K)=(-CR(K)*V(K+1)+RSV(K))*RCML
      ENDDO
!----------------------------------------------------------------------
!
      END SUBROUTINE VDIFV
!
!-----------------------------------------------------------------------
!
!=======================================================================
!
      SUBROUTINE MYJPBL_INIT
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
      CALL HERFTBL
      CALL TABLEPT
      CALL TABLETT
!-----------------------------------------------------------------------
!
      END SUBROUTINE MYJPBL_INIT
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        SUBROUTINE HERFTBL
!     ******************************************************************
!     *                                                                *
!     *               POSITIVE HALF OF ERF FUNCTION TABLE              *
!     *               RESPONSIBLE PERSON: Z.JANJIC                     *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IMPLICIT NONE

!--LOCAL VARIABLES------------------------------------------------------
INTEGER(KIND=KINT):: &
 K                           ! INDEX

REAL(KIND=KFPT):: &
 DXE &                       ! ARGUMENT INCREMENT
,DXH &                       !
,RDEN &                      !
,X &                         ! ARGUMENT
,XGMIN                       !

REAL(KIND=KFPT),DIMENSION(1:KERFM+1):: &
 GAUSS                       !

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      DXE=(XEMAX-XEMIN)/FLOAT(KERFM-1)
      DXH=DXE*0.5
      XGMIN=XEMIN-DXH
!
      RDXE=1./DXE
!
      DO K=1,KERFM+1
        X=(K-1)*DXE+XGMIN
        GAUSS(K)=EXP(-X*X*0.5)
      ENDDO
!
      RDEN=1./SQRT(PI+PI)
      HERFF(1)=0.
      DO K=2,KERFM
        HERFF(K)=((GAUSS(K)+GAUSS(K+1))*DXH)*RDEN+HERFF(K-1)
      ENDDO
!
      DO K=1,KERFM
        HERFF(K)=0.5-HERFF(K)
      ENDDO
!-----------------------------------------------------------------------
                        ENDSUBROUTINE HERFTBL
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        FUNCTION HERF(X)
!     ******************************************************************
!     *                                                                *
!     *               ERF FUNCTION HALF-TABLE                          *
!     *               RESPONSIBLE PERSON: Z.JANJIC                     *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IMPLICIT NONE

!-----------------------------------------------------------------------
REAL(KIND=KFPT):: &
 HERF

!--LOCAL VARIABLES------------------------------------------------------
INTEGER(KIND=KINT):: &
 K                           ! INDEX

REAL(KIND=KFPT):: &
 AK &                        ! POSITION IN TABLE
,X                           ! ARGUMENT
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      AK=(X-XEMIN)*RDXE
      K=INT(AK)
      K=MAX(K,0)
      K=MIN(K,KERFM2)
!
      HERF=(HERFF(K+2)-HERFF(K+1))*(AK-REAL(K))+HERFF(K+1)
!-----------------------------------------------------------------------
                        ENDFUNCTION HERF
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        SUBROUTINE TABLEPT
!     ******************************************************************
!     *                                                                *
!     *    GENERATES THE TABLE FOR FINDING PRESSURE FROM               *
!     *    SATURATION SPECIFIC HUMIDITY AND POTENTIAL TEMPERATURE      *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IMPLICIT NONE

!-----------------------------------------------------------------------
REAL(KIND=KFPT),PARAMETER:: &
 EPS=1.E-10
!--LOCAL VARIABLES------------------------------------------------------
INTEGER(KIND=KINT):: &
 KTH &                       ! INDEX
,KP                          ! INDEX

REAL(KIND=KFPT):: &
 RXNER &                     ! 1./EXNER FUNCTION
,DTH &                       ! POTENTIAL TEMPERATURE STEP
,DP &                        ! PRESSURE STEP
,DQS &                       ! SATURATION SPECIFIC HUMIDITY STEP
,P &                         ! PRESSURE
,QS0K &                      ! BASE VALUE FOR SATURATION HUMIDITY
,SQSK &                      ! SATURATION SPEC. HUMIDITY RANGE
,TH                          ! POTENTIAL TEMPERATURE

REAL(KIND=KFPT),DIMENSION(1:ITBL):: &
 APP &                       ! TEMPORARY
,AQP &                       ! TEMPORARY
,PNEW &                      ! NEW PRESSURES
,POLD &                      ! OLD PRESSURE
,QSNEW &                     ! NEW SATURATION SPEC. HUMIDITY
,QSOLD &                     ! OLD SATURATION SPEC. HUMIDITY
,Y2P                         ! TEMPORARY
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      DTH=(THH-THL)/REAL(JTBL-1)
      DP=(PH-PL)/REAL(ITBL-1)
      RDTH=1./DTH
!-----------------------------------------------------------------------
      TH=THL-DTH
      DO KTH=1,JTBL
        TH=TH+DTH
        P=PL-DP
        DO KP=1,ITBL
          P=P+DP
          RXNER=(100000./P)**CAPPA
          QSOLD(KP)=PQ0/P*EXP(A2*(TH-A3*RXNER)/(TH-A4*RXNER))
          POLD(KP)=P
        ENDDO
!
        QS0K=QSOLD(1)
        SQSK=QSOLD(ITBL)-QSOLD(1)
        QSOLD(1)=0.
        QSOLD(ITBL)=1.
!
        DO KP=2,ITBL-1
          QSOLD(KP)=(QSOLD(KP)-QS0K)/SQSK
!WWWWWWWWWWWWWW FIX DUE TO 32 BIT PRECISION LIMITATION WWWWWWWWWWWWWWWWW
          IF((QSOLD(KP)-QSOLD(KP-1)).LT.EPS) THEN
            QSOLD(KP)=QSOLD(KP-1)+EPS
          ENDIF
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
        ENDDO
!
        QS0(KTH)=QS0K
        SQS(KTH)=SQSK
!
        QSNEW(1)=0.
        QSNEW(ITBL)=1.
        DQS=1./REAL(ITBL-1)
        RDQ=1./DQS
!
        DO KP=2,ITBL-1
          QSNEW(KP)=QSNEW(KP-1)+DQS
        ENDDO
!
        Y2P(1)=0.
        Y2P(ITBL)=0.
!
        CALL SPLINE(JTBL,ITBL,QSOLD,POLD,Y2P,ITBL,QSNEW,PNEW,APP,AQP)
!
        DO KP=1,ITBL
          PTBL(KP,KTH)=PNEW(KP)
        ENDDO
!-----------------------------------------------------------------------
      ENDDO
!-----------------------------------------------------------------------
                        ENDSUBROUTINE TABLEPT
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        SUBROUTINE TABLETT
!     ******************************************************************
!     *                                                                *
!     *    GENERATES THE TABLE FOR FINDING TEMPERATURE FROM            *
!     *    PRESSURE AND EQUIVALENT POTENTIAL TEMPERATURE               *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IMPLICIT NONE

!-----------------------------------------------------------------------
REAL(KIND=KFPT),PARAMETER:: &
 EPS=1.E-10
!--LOCAL VARIABLES------------------------------------------------------
INTEGER(KIND=KINT):: &
 KTH &                       ! INDEX
,KP                          ! INDEX

REAL(KIND=KFPT):: &
 RXNER &                     ! 1./EXNER FUNCTION
,DTH &                       ! POTENTIAL TEMPERATURE STEP
,DP &                        ! PRESSURE STEP
,DTHE &                      ! EQUIVALENT POT. TEMPERATURE STEP
,P &                         ! PRESSURE
,QS &                        ! SATURATION SPECIFIC HUMIDITY
,THE0K &                     ! BASE VALUE FOR EQUIVALENT POT. TEMPERATURE
,STHEK &                     ! EQUIVALENT POT. TEMPERATURE RANGE
,TH                          ! POTENTIAL TEMPERATURE

REAL(KIND=KFPT),DIMENSION(1:JTBL):: &
 APT &                       ! TEMPORARY
,AQT &                       ! TEMPORARY
,TNEW &                      ! NEW TEMPERATURE
,TOLD &                      ! OLD TEMPERATURE
,THENEW &                    ! NEW EQUIVALENT POTENTIAL TEMPERATURE
,THEOLD &                    ! OLD EQUIVALENT POTENTIAL TEMPERATURE
,Y2T                         ! TEMPORARY
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      DTH=(THH-THL)/REAL(JTBL-1)
      DP=(PH-PL)/REAL(ITBL-1)
      RDP=1./DP
!-----------------------------------------------------------------------
      P=PL-DP
      DO KP=1,ITBL
        P=P+DP
        TH=THL-DTH
        DO KTH=1,JTBL
          TH=TH+DTH
          RXNER=(100000./P)**CAPPA
          QS=PQ0/P*EXP(A2*(TH-A3*RXNER)/(TH-A4*RXNER))
          TOLD(KTH)=TH/RXNER
          THEOLD(KTH)=TH*EXP(ELIWV*QS/(CP*TOLD(KTH)))
        ENDDO
!
        THE0K=THEOLD(1)
        STHEK=THEOLD(JTBL)-THEOLD(1)
        THEOLD(1)=0.
        THEOLD(JTBL)=1.
!
        DO KTH=2,JTBL-1
          THEOLD(KTH)=(THEOLD(KTH)-THE0K)/STHEK
!WWWWWWWWWWWWWW FIX DUE TO 32 BIT PRECISION LIMITATION WWWWWWWWWWWWWWWWW
          IF((THEOLD(KTH)-THEOLD(KTH-1)).LT.EPS) THEN
            THEOLD(KTH)=THEOLD(KTH-1)+EPS
          ENDIF
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
        ENDDO
!
        THE0(KP)=THE0K
        STHE(KP)=STHEK
!
        THENEW(1)=0.
        THENEW(JTBL)=1.
        DTHE=1./REAL(JTBL-1)
        RDTHE=1./DTHE
!
        DO KTH=2,JTBL-1
          THENEW(KTH)=THENEW(KTH-1)+DTHE
        ENDDO
!
        Y2T(1)=0.
        Y2T(JTBL)=0.
!
        CALL SPLINE(JTBL,JTBL,THEOLD,TOLD,Y2T,JTBL,THENEW,TNEW,APT,AQT)
!
        DO KTH=1,JTBL
          TTBL(KTH,KP)=TNEW(KTH)
        ENDDO
!-----------------------------------------------------------------------
      ENDDO
!-----------------------------------------------------------------------
                        ENDSUBROUTINE TABLETT
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      subroutine spline(jtbx,nold,xold,yold,y2,nnew,xnew,ynew,p,q)
!   ********************************************************************
!   *                                                                  *
!   *  this is a one-dimensional cubic spline fitting routine          *
!   *  programed for a small scalar machine.                           *
!   *                                                                  *
!   *  programer z. janjic                                             *
!   *                                                                  *
!   *  nold - number of given values of the function.  must be ge 3.   *
!   *  xold - locations of the points at which the values of the       *
!   *         function are given.  must be in ascending order.         *
!   *  yold - the given values of the function at the points xold.     *
!   *  y2   - the second derivatives at the points xold.  if natural   *
!   *         spline is fitted y2(1)=0. and y2(nold)=0. must be        *
!   *         specified.                                               *
!   *  nnew - number of values of the function to be calculated.       *
!   *  xnew - locations of the points at which the values of the       *
!   *         function are calculated.  xnew(k) must be ge xold(1)     *
!   *         and le xold(nold).                                       *
!   *  ynew - the values of the function to be calculated.             *
!   *  p, q - auxiliary vectors of the length nold-2.                  *
!   *                                                                  *
!   ********************************************************************
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer(kind=kint),intent(in):: &
       jtbx,nnew,nold

      real(kind=kfpt),dimension(jtbx),intent(in):: &
       xnew,xold,yold

      real(kind=kfpt),dimension(jtbx),intent(inout):: &
       p,q,y2

      real(kind=kfpt),dimension(jtbx),intent(out):: &
       ynew
!
      integer(kind=kint):: &
       k,k1,k2,kold,noldm1

      real(kind=kfpt):: &
       ak,bk,ck,den,dx,dxc,dxl,dxr,dydxl,dydxr &
      ,rdx,rtdxc,x,xk,xsq,y2k,y2kp1
!-----------------------------------------------------------------------
      noldm1=nold-1
!
      dxl=xold(2)-xold(1)
      dxr=xold(3)-xold(2)
      dydxl=(yold(2)-yold(1))/dxl
      dydxr=(yold(3)-yold(2))/dxr
      rtdxc=0.5/(dxl+dxr)
!
      p(1)= rtdxc*(6.*(dydxr-dydxl)-dxl*y2(1))
      q(1)=-rtdxc*dxr
!
      if(nold==3)go to 150
!-----------------------------------------------------------------------
      k=3
!
  100 dxl=dxr
      dydxl=dydxr
      dxr=xold(k+1)-xold(k)
      dydxr=(yold(k+1)-yold(k))/dxr
      dxc=dxl+dxr
      den=1./(dxl*q(k-2)+dxc+dxc)
!
      p(k-1)= den*(6.*(dydxr-dydxl)-dxl*p(k-2))
      q(k-1)=-den*dxr
!
      k=k+1
      if(k<nold)go to 100
!-----------------------------------------------------------------------
  150 k=noldm1
!
  200 y2(k)=p(k-1)+q(k-1)*y2(k+1)
!
      k=k-1
      if(k>1)go to 200
!-----------------------------------------------------------------------
      k1=1
!
  300 xk=xnew(k1)
!
      do 400 k2=2,nold
!
      if(xold(k2)>xk)then
        kold=k2-1
        go to 450
      endif
!
  400 continue
!
      ynew(k1)=yold(nold)
      go to 600
!
  450 if(k1==1)go to 500
      if(k==kold)go to 550
!
  500 k=kold
!
      y2k=y2(k)
      y2kp1=y2(k+1)
      dx=xold(k+1)-xold(k)
      rdx=1./dx
!
      ak=.1666667*rdx*(y2kp1-y2k)
      bk=0.5*y2k
      ck=rdx*(yold(k+1)-yold(k))-.1666667*dx*(y2kp1+y2k+y2k)
!
  550 x=xk-xold(k)
      xsq=x*x
!
      ynew(k1)=ak*xsq*x+bk*xsq+ck*x+yold(k)
!
  600 k1=k1+1
      if(k1<=nnew)go to 300
!-----------------------------------------------------------------------
      end subroutine spline
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      END MODULE MODULE_BL_MYJPBL
!
!-----------------------------------------------------------------------
