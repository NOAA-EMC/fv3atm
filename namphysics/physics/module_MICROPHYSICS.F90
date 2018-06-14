!-----------------------------------------------------------------------
!
      MODULE MODULE_MICROPHYSICS_NMM
!
!-----------------------------------------------------------------------
!
!***  THE MICROPHYSICS DRIVERS AND PACKAGES

!      11-06-2009 W. Wang put NAM micorphysics into a single module
!      02-10-2010 W. Wang added wsm6
!-----------------------------------------------------------------------
!
! HISTORY LOG:
!
!    11-06-2009 W. Wang - Put NAM/Ferrier microphysics into
!                         a single module.
!
!-----------------------------------------------------------------------
!
      USE MODULE_KINDS
!
      USE MODULE_CONSTANTS,ONLY : CICE,CLIQ,CPV,EP_1,EP_2,EPSILON,G     &
                                 ,P608,PSAT,R_D,R_V,RHOAIR0,RHOWATER    &
                                 ,SVPT0,XLF,XLV                         &
                                 ,CAPPA,CP,EPSQ
!
! MP options
      USE MODULE_MP_ETANEW
      USE MODULE_MP_FER_HIRES
      USE MODULE_MP_WSM6
      USE MODULE_MP_THOMPSON
      USE MODULE_MP_GFS
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: GSMDRIVE,UPDATE_WATER,TQADJUST
!
!-----------------------------------------------------------------------
!
      INTEGER :: MYPE
      REAL, PRIVATE,PARAMETER ::                                     &
!--- Physical constants follow:
           XLS=2.834E6,R_G=1./G
!
      CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE GSMDRIVE(ITIMESTEP,DT,NPHS                             &
                         ,SM,FIS                                        &
                         ,PRSI,P_PHY                                    &
                         ,T,Q,CWM                                       &
                         ,TRAIN,SR                                      &
                         ,F_ICE,F_RAIN,F_RIMEF                          &
                         ,QC,QR,QI,QS,QG,NI,NR                          &
                         ,F_QC,F_QR,F_QI,F_QS,F_QG,F_NI,F_NR            &
                         ,has_reqc, has_reqi, has_reqs                  &
                         ,PREC,ACPREC,AVRAIN                            &
                         ,refl_10cm                                     &
                         ,re_cloud,re_ice,re_snow                       &
                         ,MICROPHYSICS                                  &
                         ,RHGRD                                         &
                         ,TP1,QP1,PSP1                                  &
                         ,IMS,IME,LM)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    GSMDRIVE    MICROPHYSICS OUTER DRIVER
!   PRGRMMR: BLACK           ORG: W/NP22     DATE: 02-03-26
!
! ABSTRACT:
!     RADIATION SERVES AS THE INTERFACE BETWEEN THE NMMB PHYSICS COMPONENT
!     AND THE WRF MICROPHYSICS DRIVER.
!
! PROGRAM HISTORY LOG:
!   02-03-26  BLACK      - ORIGINATOR
!   04-11-18  BLACK      - THREADED
!   06-07-31  BLACK      - BUILT INTO NMMB PHYSICS COMPONENT
!   08-08     JANJIC     - Synchronize WATER array and Q.
!
! USAGE: CALL GSMDRIVE FROM PHY_RUN
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
      INTEGER,PARAMETER :: JMS=1,JME=1,D_SS=1
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER,INTENT(IN) :: ITIMESTEP,NPHS                              &
                           ,IMS,IME,LM                                  &
                           ,has_reqc,has_reqi,has_reqs
!
!     LOGICAL,INTENT(IN) :: USE_RADAR
!
      REAL,INTENT(IN) :: DT,RHGRD
!
      REAL,INTENT(INOUT) :: AVRAIN
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM,D_SS)  :: MPRATES
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: FIS,SM
!
!     REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: DFI_TTEN
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: ACPREC,PREC
!
      REAL,DIMENSION(IMS:IME,1:LM)  ,INTENT(IN) :: P_PHY
      REAL,DIMENSION(IMS:IME,1:LM+1),INTENT(IN) :: PRSI
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: refl_10cm
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: re_cloud, re_ice, re_snow
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: CWM,Q,T     &
                                                           ,TRAIN       &
                                            ,QC,QI,QR,QS,QG,NI,NR
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: F_ICE       &
                                                           ,F_RAIN      &
                                                           ,F_RIMEF
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: SR
!
      CHARACTER(99),INTENT(IN) :: MICROPHYSICS
!
      LOGICAL,INTENT(IN) :: F_QC,F_QR,F_QI,F_QS,F_QG,F_NI,F_NR
!
!*** GFS microphysics
      REAL, DIMENSION(IMS:IME,JMS:JME,1:LM), INTENT(INOUT) :: TP1,QP1
      REAL, DIMENSION(IMS:IME,JMS:JME), INTENT(INOUT)      :: PSP1
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: I,IJ,J,K,MP_PHYSICS,N,NTSD
!
      INTEGER,DIMENSION(IMS:IME,JMS:JME) :: LOWLYR
!
      REAL :: DTPHS,PCPCOL,QW,RDTPHS,TNEW
      REAL :: MP_TTEN,mytten
!
      REAL,DIMENSION(1:LM) :: QL,TL
!
      REAL,DIMENSION(IMS:IME,JMS:JME) :: CUBOT,CUTOP,RAINNC,RAINNCV     &
                                        ,SNOWNC,SNOWNCV,XLAND           &
                                        ,graupelnc,graupelncv
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM) :: DZ,PI_PHY                 &
                                             ,RR,TH_PHY,QV
!
      LOGICAL :: WARM_RAIN,F_QT,USE_QV
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      NTSD=ITIMESTEP
      DTPHS=NPHS*DT
      RDTPHS=1./DTPHS
      AVRAIN=AVRAIN+1.
!
!-----------------------------------------------------------------------
!***  NOTE:  THE NMMB HAS IJK STORAGE WITH LAYER 1 AT THE TOP.
!***         THE WRF PHYSICS DRIVERS HAVE IKJ STORAGE WITH LAYER 1
!***         AT THE BOTTOM.
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(j,i,k,ql,tl)
!.......................................................................
      DO J=JMS,JME
      DO I=IMS,IME
!
        LOWLYR(I,J)=1
        XLAND(I,J)=SM(I,J)+1.
!
!-----------------------------------------------------------------------
!***   FILL RAINNC WITH ZERO (NORMALLY CONTAINS THE NONCONVECTIVE
!***                          ACCUMULATED RAIN BUT NOT YET USED BY NMM)
!***   COULD BE OBTAINED FROM ACPREC AND CUPREC (ACPREC-CUPREC)
!-----------------------------------------------------------------------
!..The NC variables were designed to hold simulation total accumulations
!.. whereas the NCV variables hold timestep only values, so change below
!.. to zero out only the timestep amount preparing to go into each
!.. micro routine while allowing NC vars to accumulate continually.
!.. But, the fact is, the total accum variables are local, never saved
!.. nor written so they go nowhere at the moment.
!
        RAINNC (I,J)=0. ! NOT YET USED BY NMM
        RAINNCv(I,J)=0.
        SNOWNCv(I,J)=0.
        graupelncv(i,j) = 0.0
!
!-----------------------------------------------------------------------
!***  FILL THE SINGLE-COLUMN INPUT
!-----------------------------------------------------------------------
!
        DO K=LM,1,-1   ! We are moving down from the top in the flipped arrays
!
          TL(K)=T(I,J,K)
          QL(K)=AMAX1(Q(I,J,K),EPSQ)
!
          RR(I,J,K)=P_PHY(I,K)/(R_D*TL(K)*(P608*QL(K)+1.))
          PI_PHY(I,J,K)=(P_PHY(I,K)*1.E-5)**CAPPA
          TH_PHY(I,J,K)=TL(K)/PI_PHY(I,J,K)
          DZ(I,J,K)=(PRSI(I,K+1)-PRSI(I,K))*R_G/RR(I,J,K)
!
        ENDDO    !- DO K=LM,1,-1
!
      ENDDO    !- DO I=IMS,IME
      ENDDO    !- DO J=JMS,JME
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  IF NEEDED, UPDATE WATER VAPOR RATIO FROM SPECIFIC HUMIDITY.
!-----------------------------------------------------------------------
!
      IF(TRIM(MICROPHYSICS)=='wsm6' .OR. TRIM(MICROPHYSICS)=='thompson')THEN
        USE_QV=.TRUE.    !-- Initialize QV, update Q & CWM at the end
      ELSE
        USE_QV=.FALSE.
      ENDIF
!
      IF(USE_QV) THEN
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,k)
!.......................................................................
        DO K=1,LM
          DO J=JMS,JME
          DO I=IMS,IME
            QV(I,J,K)=Q(I,J,K)/(1.-Q(I,J,K))
          ENDDO
          ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
      ENDIF
!
!-----------------------------------------------------------------------
!
!***  CALL MICROPHYSICS
!
!-----------------------------------------------------------------------
!
!---------------------------------------------------------------------
!  Check for microphysics type.  We need a clean way to
!  specify these things!
!---------------------------------------------------------------------
!
        micro_select: SELECT CASE (TRIM(MICROPHYSICS))
!
          CASE ('fer')
            CALL ETAMP_NEW(                                             &
                   ITIMESTEP=ntsd,DT=dtphs                              &
                  ,DZ8W=dz,RHO_PHY=rr,P_PHY=p_phy,PI_PHY=pi_phy         &
                  ,TH_PHY=th_phy                                        &
                  ,Q=Q,QC=QC,QS=QS,QR=QR,QT=cwm                         &
                  ,LOWLYR=LOWLYR,SR=SR                                  &
                  ,F_ICE_PHY=F_ICE,F_RAIN_PHY=F_RAIN                    &
                  ,F_RIMEF_PHY=F_RIMEF                                  &
                  ,RAINNC=rainnc,RAINNCV=rainncv                        &
                  ,IMS=IMS,IME=IME, JMS=JMS,JME=JME, LM=LM              &
                  ,D_SS=d_ss,MPRATES=mprates)
          CASE ('fer_hires')
!---------------------------------------------------------------------
!*** Update the rime factor array after 3d advection
!---------------------------------------------------------------------
              DO K=1,LM
              DO J=JMS,JME
              DO I=IMS,IME
                IF (QG(I,J,K)>EPSQ .AND. QS(I,J,K)>EPSQ) THEN
                  F_RIMEF(I,J,K)=MIN(50.,MAX(1.,QG(I,J,K)/QS(I,J,K)))
                ELSE
                  F_RIMEF(I,J,K)=1.
                ENDIF
              ENDDO
              ENDDO
              ENDDO
!---------------------------------------------------------------------

            CALL FER_HIRES(                                             &
                   ITIMESTEP=ntsd,DT=dtphs,RHgrd=RHGRD                  &
                  ,DZ8W=dz,RHO_PHY=rr,P_PHY=p_phy,PI_PHY=pi_phy         &
                  ,TH_PHY=th_phy                                        &
                  ,Q=Q,QC=QC,QS=QS,QR=QR,QT=cwm                         &
                  ,LOWLYR=LOWLYR,SR=SR                                  &
                  ,F_ICE_PHY=F_ICE,F_RAIN_PHY=F_RAIN                    &
                  ,F_RIMEF_PHY=F_RIMEF                                  &
                  ,RAINNC=rainnc,RAINNCV=rainncv                        &
                  ,IMS=IMS,IME=IME,JMS=JMS,JME=JME,LM=LM                &
                  ,D_SS=d_ss,MPRATES=mprates                            &
                  ,refl_10cm=refl_10cm)

!---------------------------------------------------------------------
!*** Calculate graupel from snow array and rime factor
!---------------------------------------------------------------------
              DO K=1,LM
              DO J=JMS,JME
              DO I=IMS,IME
                QG(I,J,K)=QS(I,J,K)*F_RIMEF(I,J,K)
              ENDDO
              ENDDO
              ENDDO
!---------------------------------------------------------------------
          CASE ('gfs')
            CALL GFSMP(DT=dtphs,                                        &
                   dz8w=dz,rho_phy=rr,p_phy=p_phy,pi_phy=pi_phy,        &
                   th_phy=th_phy,                                       &
                   SR=SR,QT=CWM, F_ICE_PHY=F_ICE,                       &
                   RAINNC=RAINNC,RAINNCV=RAINNCV,                       &
                   Q=Q,QC=QC,QI=QI,                                     &
                   F_QC=F_QC,F_QI=F_QI,                                 &
                   TP1=TP1,QP1=QP1,PSP1=PSP1,                           &
                   IMS=IMS,IME=IME, JMS=JMS,JME=JME, LM=LM )
          CASE ('wsm6')
             CALL wsm6(                                                 &
                  TH=th_phy                                             &
                 ,Q=QV                                                  &
                 ,QC=QC                                                 &
                 ,QR=QR                                                 &
                 ,QI=QI                                                 &
                 ,QS=QS                                                 &
                 ,QG=QG                                                 &
                 ,DEN=rr,PII=pi_phy,P=p_phy,DELZ=dz                     &
                 ,DELT=dtphs,G=g,CPD=cp,CPV=cpv                         &
                 ,RD=r_d,RV=r_v,T0C=svpt0                               &
                 ,EP1=ep_1, EP2=ep_2, QMIN=epsilon                      &
                 ,XLS=xls, XLV0=xlv, XLF0=xlf                           &
                 ,DEN0=rhoair0, DENR=rhowater                           &
                 ,CLIQ=cliq,CICE=cice,PSAT=psat                         &
                 ,RAIN=rainnc ,RAINNCV=rainncv                          &
                 ,SNOW=snownc ,SNOWNCV=snowncv                          &
                 ,SR=sr                                                 &
                 ,GRAUPEL=graupelnc ,GRAUPELNCV=graupelncv              &
                 ,IMS=IMS,IME=IME, JMS=JMS,JME=JME, LM=LM               &
                 ,D_SS=d_ss,MPRATES=mprates)
          CASE ('thompson')
!+---+-----------------------------------------------------------------+
!            write(6,*)'DEBUG-GT, calling mp_gt_driver'
             CALL mp_gt_driver(                                         &
                  qv=qv                                                 &
                 ,qc=qc                                                 &
                 ,qr=qr                                                 &
                 ,qi=qi                                                 &
                 ,qs=qs                                                 &
                 ,qg=qg                                                 &
                 ,ni=ni                                                 &
                 ,nr=nr                                                 &
                 ,TH=th_phy,PII=pi_phy,P=p_phy,dz=dz,dt_in=dtphs        &
                 ,itimestep=ntsd                                        &
                 ,RAINNC=rainnc ,RAINNCV=rainncv                        &
                 ,SNOWNC=snownc ,SNOWNCV=snowncv                        &
                 ,GRAUPELNC=graupelnc ,GRAUPELNCV=graupelncv            &
                 ,SR=sr                                                 &
                 ,refl_10cm=refl_10cm(ims,jms,1)                        &
                 ,diagflag=.true.                                       &
                 ,do_radar_ref=1                                        &
                 ,re_cloud=re_cloud(ims,jms,1)                          &
                 ,re_ice=re_ice(ims,jms,1)                              &
                 ,re_snow=re_snow(ims,jms,1)                            &
                 ,has_reqc=has_reqc                                     &
                 ,has_reqi=has_reqi                                     &
                 ,has_reqs=has_reqs                                     &
                 ,IMS=IMS,IME=IME, JMS=JMS,JME=JME, KTS=1,KTE=LM        &
                 ,D_SS=d_ss,MPRATES=mprates                         )
!
!+---+-----------------------------------------------------------------+

          CASE DEFAULT
            WRITE(0,*)' The microphysics option does not exist: MICROPHYSICS = ',TRIM(MICROPHYSICS)

        END SELECT micro_select

!
!-----------------------------------------------------------------------
!
      IF(USE_QV) THEN    !-- Update Q & CWM for WSM6 & Thompson microphysics
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,k)
!.......................................................................
        DO K=1,LM
          DO J=JMS,JME
          DO I=IMS,IME
            Q(I,J,K)=QV(I,J,K)/(1.+QV(I,J,K))
            CWM(I,J,K)=QC(i,j,k)+QR(i,j,k)+QI(i,j,k)+QS(i,j,k)+QG(i,j,k)
          ENDDO
          ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
      ENDIF
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,k,TNEW,MP_TTEN)
!.......................................................................
      DO K=1,LM
        DO J=JMS,JME
        DO I=IMS,IME
!
!-----------------------------------------------------------------------
!***  UPDATE TEMPERATURE, SPECIFIC HUMIDITY, CLOUD WATER, AND HEATING.
!-----------------------------------------------------------------------
!
          TNEW=TH_PHY(I,J,K)*PI_PHY(I,J,K)
          TRAIN(I,J,K)=TRAIN(I,J,K)+(TNEW-T(I,J,K))*RDTPHS
!         IF (USE_RADAR) THEN
!           MP_TTEN=(TNEW-T(I,J,K))*RDTPHS
!           IF(DFI_TTEN(I,J,K)>MP_TTEN.AND.DFI_TTEN(I,J,K)<0.01        &
!                                     .AND.MP_TTEN<0.0018)THEN
!             MP_TTEN=DFI_TTEN(I,J,K)
!           END IF
!           T(I,J,K)=T(I,J,K)+MP_TTEN/RDTPHS
!         ELSE
            T(I,J,K)=TNEW
!         ENDIF
        ENDDO
        ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  UPDATE PRECIPITATION
!-----------------------------------------------------------------------
!
!jaa!$omp parallel do                                                       &
!jaa!$omp& private(i,j,pcpcol)
      DO J=JMS,JME
      DO I=IMS,IME
        PCPCOL=RAINNCV(I,J)*1.E-3
        PREC(I,J)=PREC(I,J)+PCPCOL
        ACPREC(I,J)=ACPREC(I,J)+PCPCOL
!
! NOTE: RAINNC IS ACCUMULATED INSIDE MICROPHYSICS BUT NMM ZEROES IT OUT ABOVE
!       SINCE IT IS ONLY A LOCAL ARRAY FOR NOW
!
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GSMDRIVE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE UPDATE_WATER(CWM,F_ICE,F_RAIN,F_RIMEF                  &
                             ,T,QC,QR,QS,QI,QG                          &
                             ,MICROPHYSICS,SPEC_ADV,NTIMESTEP           &
                             ,LM,IMS,IME)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    UPDATE_WATER          UPDATE WATER ARRAY
!   PRGRMMR: FERRIER         ORG: NP22     DATE: 3 AUG 2009
!
! ABSTRACT:
!     UPDATE WATER ARRAY FOR FERRIER MICROPHYSICS
!
! PROGRAM HISTORY LOG (with changes to called routines) :
!   2009-08     FERRIER     - Synchronize WATER array with CWM, F_rain, F_ice arrays
!
! USAGE: CALL UPDATE_WATER FROM PHY_RUN
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!-----------------------------------------------------------------------
      USE MODULE_CONSTANTS,ONLY : EPSQ,TIW
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
!
      INTEGER,PARAMETER :: JMS=1,JME=1
!
!----------------------
!-- Argument Variables
!----------------------
!
      INTEGER,INTENT(IN) :: NTIMESTEP,LM,IMS,IME
!
      CHARACTER(99),INTENT(IN) :: MICROPHYSICS
!
      LOGICAL,INTENT(IN) :: SPEC_ADV
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: CWM         &
                                                           ,F_ICE       &
                                                           ,F_RAIN      &
                                                           ,F_RIMEF     &
                                                           ,T,QC,QR     &
                                                           ,QS,QI,QG
!
!--------------------
!--  Local Variables
!--------------------
!
      INTEGER :: I,J,K
      REAL :: FRACTION, LIQW, OLDCWM
      LOGICAL :: CLD_INIT
      LOGICAL :: deep_ice
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      IF(NTIMESTEP<=1)THEN
        CLD_INIT=.TRUE.
      ELSE
        CLD_INIT=.FALSE.
      ENDIF
!
!----------------------------------------------------------------------
!-- Couple 2 sets of condensed water arrays for different microphysics:
!   QC,QR,QS, etc. arrays <=> CWM,F_ice,F_rain,F_RimeF 3D arrays
!----------------------------------------------------------------------
!
      SELECT CASE ( TRIM(MICROPHYSICS) )
!
!----------------------------------------------------------------------
        CASE ('fer','fer_hires')  !-- Update fields for Ferrier microphysics
!----------------------------------------------------------------------
!
          spec_adv_fer: IF (.NOT.SPEC_ADV .OR. CLD_INIT) THEN
!-- Update WATER arrays when advecting only total condensate (spec_adv=F)
!   or at the initial time step
            DO K=1,LM
             DO J=JMS,JME
              DO I=IMS,IME
                IF (CWM(I,J,K)>EPSQ) THEN
                  LIQW=(1.-F_ice(I,J,K))*CWM(I,J,K)
                  QC(I,J,K)=(1.-F_rain(I,J,K))*LIQW
                  QR(I,J,K)=F_rain(I,J,K)*LIQW
                  QS(I,J,K)=F_ice(I,J,K)*CWM(I,J,K)
                ELSE
                  QC(I,J,K)=0.
                  QR(I,J,K)=0.
                  QS(I,J,K)=0.
                ENDIF
              ENDDO
             ENDDO
            ENDDO
!
          ELSE spec_adv_fer
!-- Update CWM,F_ICE,F_RAIN arrays from separate species advection (spec_adv=T)
            DO K=1,LM
             DO J=JMS,JME
              DO I=IMS,IME
                CWM(I,J,K)=QC(I,J,K)+QR(I,J,K)+QS(I,J,K)
                IF (QS(I,J,K)>EPSQ) THEN
                  F_ICE(I,J,K)=QS(I,J,K)/CWM(I,J,K)
                ELSE
                  F_ICE(I,J,K)=0.0
                ENDIF
                IF (QR(I,J,K)>EPSQ) THEN
                  F_RAIN(I,J,K)=QR(I,J,K)/(QC(I,J,K)+QR(I,J,K))
                ELSE
                  F_RAIN(I,J,K)=0.
                ENDIF
              ENDDO
             ENDDO
            ENDDO
          ENDIF spec_adv_fer
!
!----------------------------------------------------------------------
        CASE ('gfs')       !-- Update fields for GFS microphysics
!----------------------------------------------------------------------
!
          spec_adv_gfs: IF (.NOT.SPEC_ADV .OR. CLD_INIT) THEN
            cld_init_gfs: IF (CLD_INIT) THEN
!-- Initialize F_ICE, F_RAIN, & F_RIMEF arrays
              IF (SPEC_ADV) THEN
                WRITE(0,*) 'Never ran GFS microphysics with SPEC_ADV=T.'   &
                          ,'  Use at your own risk.'
              ENDIF
              DO K=1,LM
               DO J=JMS,JME
                DO I=IMS,IME
                  F_RAIN(I,J,K)=0.
                  F_RIMEF(I,J,K)=1.
                  IF (CWM(I,J,K)>EPSQ .AND. T(I,J,K)<233.15) THEN
                    F_ICE(I,J,K)=1.
                  ELSE
                    F_ICE(I,J,K)=0.
                  ENDIF
                ENDDO
               ENDDO
              ENDDO
            ENDIF  cld_init_gfs
!-- Update WATER arrays (QC,QI) when advecting only total condensate (spec_adv=F)
!   or initialize them at the start of the forecast (CLD_INIT=T).
            DO K=1,LM
             DO J=JMS,JME
              DO I=IMS,IME
                IF (CWM(I,J,K)>EPSQ) THEN
                  QC(I,J,K)=(1.-F_ice(I,J,K))*CWM(I,J,K)
                  QI(I,J,K)=F_ice(I,J,K)*CWM(I,J,K)
                ELSE
                  QC(I,J,K)=0.
                  QI(I,J,K)=0.
                ENDIF
              ENDDO
             ENDDO
            ENDDO
          ELSE spec_adv_gfs
!-- Update CWM, F_ICE arrays from separate species advection (spec_adv=T)
            DO K=1,LM
             DO J=JMS,JME
              DO I=IMS,IME
                CWM(I,J,K)=QC(I,J,K)+QI(I,J,K)
                IF (CWM(I,J,K)>EPSQ) THEN
                  F_ICE(I,J,K)=QI(I,J,K)/CWM(I,J,K)
                ELSE
                  F_ICE(I,J,K)=0.
                ENDIF
              ENDDO
             ENDDO
            ENDDO
          ENDIF  spec_adv_gfs
!
!----------------------------------------------------------------------
        CASE ('wsm6')      !-- Update fields for WSM6 microphysics
!----------------------------------------------------------------------
!
          init_adv_wsm6: IF (CLD_INIT) THEN
!-- Assume only cloud ice is present at initial time
            DO K=1,LM
             DO J=JMS,JME
              DO I=IMS,IME
                QS(I,J,K)=0.0
                QG(I,J,K)=0.0
                IF (CWM(I,J,K)>EPSQ) THEN
                  LIQW=(1.-F_ice(I,J,K))*CWM(I,J,K)
                  QC(I,J,K)=(1.-F_rain(I,J,K))*LIQW
                  QR(I,J,K)=F_rain(I,J,K)*LIQW
                  QI(I,J,K)=F_ice(I,J,K)*CWM(I,J,K)
                ELSE
                  QC(I,J,K)=0.
                  QR(I,J,K)=0.
                  QI(I,J,K)=0.
                ENDIF
              ENDDO
             ENDDO
            ENDDO
          ELSE init_adv_wsm6
            notspec_adv_wsm6: IF (.NOT.SPEC_ADV) THEN
!-- Update WATER arrays (QC,QR,...) when advecting only total condensate (spec_adv=F).
!-- Assume fraction of each water category is unchanged by advection.
              DO K=1,LM
               DO J=JMS,JME
                DO I=IMS,IME
                  OLDCWM=QC(I,J,K)+QR(I,J,K)   &
                        +QI(I,J,K)+QS(I,J,K)   &
                        +QG(I,J,K)
                  IF (OLDCWM>EPSQ) THEN
                    FRACTION=CWM(I,J,K)/OLDCWM
                    QC(I,J,K)=FRACTION*QC(I,J,K)
                    QR(I,J,K)=FRACTION*QR(I,J,K)
                    QI(I,J,K)=FRACTION*QI(I,J,K)
                    QS(I,J,K)=FRACTION*QS(I,J,K)
                    QG(I,J,K)=FRACTION*QG(I,J,K)
                  ELSE
                    QC(I,J,K)=0.0
                    QR(I,J,K)=0.0
                    QI(I,J,K)=0.0
                    QS(I,J,K)=0.0
                    QG(I,J,K)=0.0
                    IF (T(I,J,K)<233.15) THEN
                      QI(I,J,K)=CWM(I,J,K)
                    ELSE
                      QC(I,J,K)=CWM(I,J,K)
                    ENDIF
                  ENDIF
                ENDDO
               ENDDO
              ENDDO
            ENDIF  notspec_adv_wsm6
!
!-- Couple QC,QR,... <=> CWM,F_ice,F_rain,F_RimeF arrays
!-- Update CWM,F_XXX arrays from separate species advection (spec_adv=T)
!
            DO K=1,LM
             DO J=JMS,JME
              DO I=IMS,IME
                CWM(I,J,K)=QC(I,J,K)+QR(I,J,K)      &
                          +QI(I,J,K)+QS(I,J,K)      &
                          +QG(I,J,K)
                IF (CWM(I,J,K)>EPSQ) THEN
                  LIQW=QI(I,J,K)+QS(I,J,K)+QG(I,J,K)
                  F_ICE(I,J,K)=LIQW/CWM(I,J,K)
                ELSE
                  F_ICE(I,J,K)=0.
                ENDIF
                IF (QR(I,J,K)>EPSQ) THEN
                  F_RAIN(I,J,K)=QR(I,J,K)/(QC(I,J,K)+QR(I,J,K))
                ELSE
                  F_RAIN(I,J,K)=0.
                ENDIF
                IF (QG(I,J,K)>EPSQ) THEN
!-- Update F_RIMEF: assume 5x higher graupel density (500 kg/m**3) vs snow (100 kg/m**3)
                  LIQW=5.*QG(I,J,K)+QS(I,J,K)
                  F_RIMEF(I,J,K)=LIQW/(QS(I,J,K)+QG(I,J,K))
                ELSE
                  F_RIMEF(I,J,K)=1.
                ENDIF
              ENDDO
             ENDDO
            ENDDO
!
          ENDIF init_adv_wsm6
!
!----------------------------------------------------------------------
        CASE ('thompson')   !-- Update fields for Thompson microphysics
!----------------------------------------------------------------------
!
!+---+-----------------------------------------------------------------+
!..The CLD_INIT test provides a way to translate initial values of CWM
!.. into coomponent species of cloud water, rain, and ice, but not snow
!.. or graupel. Thompson MP will pretty rapidly make snow from the
!.. cloud ice field.  Next IF-test is whether individual species
!.. advection is enabled, which almost certainly should be the case when
!.. picking this scheme.  In this case, the separate species are summed
!.. into the CWM and ice, rain, and rime variables are computed only for
!.. consistency with other schemes.  But, if single species advection is
!.. not enabled, then each t-step the CWM array needs to be split into
!.. component species to prepare MP routine to have some semblance of
!.. proper individual species.  Again, this is strongly discouraged.
!+---+-----------------------------------------------------------------+
          spec_adv_thompson: IF (CLD_INIT) THEN
             DO K=1,LM
                DO J=JMS,JME
                DO I=IMS,IME
                   QS(I,J,K)=0.0
                   QG(I,J,K)=0.0
                   IF (CWM(I,J,K) .gt. EPSQ) THEN
                      LIQW=(1.-F_ice(I,J,K))*CWM(I,J,K)
                      QC(I,J,K)=(1.-F_rain(I,J,K))*LIQW
                      QR(I,J,K)=F_rain(I,J,K)*LIQW
                      QI(I,J,K)=F_ice(I,J,K)*CWM(I,J,K)
                   ELSE
                      QC(I,J,K)=0.
                      QR(I,J,K)=0.
                      QI(I,J,K)=0.
                   ENDIF
                ENDDO
                ENDDO
             ENDDO
          ELSE IF(SPEC_ADV) THEN  spec_adv_thompson
             DO K=1,LM
                DO J=JMS,JME
                DO I=IMS,IME
                   CWM(I,J,K) = QC(I,J,K)+QR(I,J,K)     &
                              + QI(I,J,K)                       &
                              + QS(I,J,K)+QG(I,J,K)
                   IF (CWM(I,J,K) .gt. EPSQ) THEN
                      LIQW = MAX(0., CWM(I,J,K) - QI(I,J,K)     &
                                                - QS(I,J,K)     &
                                                - QG(I,J,K))
                      F_ICE(I,J,K) = MAX(0., 1.0 - LIQW/CWM(I,J,K))
                      IF (QR(I,J,K) .gt. EPSQ) THEN
                         F_RAIN(I,J,K) = QR(I,J,K)              &
                                 / (QC(I,J,K)+QR(I,J,K))
                      ELSE
                         F_RAIN(I,J,K)=0.
                      ENDIF
                      IF (QG(I,J,K) .gt. EPSQ) THEN
                         F_RIMEF(I,J,K) = (5.*QG(I,J,K)         &
                                        +     QS(I,J,K))        &
                                        / (QS(I,J,K)            &
                                        +  QG(I,J,K))
                      ELSE
                         F_RIMEF(I,J,K)=1.
                      ENDIF
                   ELSE
                      F_ICE(I,J,K) = 0.
                      F_RAIN(I,J,K)=0.
                      F_RIMEF(I,J,K)=1.
                      CWM(I,J,K) = 0.
                   ENDIF
                ENDDO
                ENDDO
             ENDDO
          ELSE  spec_adv_thompson
            ! write(0,*) 'WARNING: This option is STRONGLY DISCOURAGED'
            ! write(0,*) '  please consider using full advection of all'
            ! write(0,*) '  species when picking Thompson microphysics.'
             DO J=JMS,JME
             DO I=IMS,IME
                DO K=LM,1,-1
                   deep_ice = .false.
                   IF (CWM(I,J,K) .gt. EPSQ) THEN
                      OLDCWM  = QC(I,J,K)+QR(I,J,K)     &
                              + QI(I,J,K)                       &
                              + QS(I,J,K)+QG(I,J,K)
                      IF (OLDCWM .gt. EPSQ) THEN
                         LIQW = MAX(0., OLDCWM - QI(I,J,K)      &
                                               - QS(I,J,K)      &
                                               - QG(I,J,K))
                         F_ICE(I,J,K) = MAX(0., 1.0 - LIQW/OLDCWM)
                         IF (QR(I,J,K) .gt. EPSQ) THEN
                            F_RAIN(I,J,K) = QR(I,J,K)           &
                                 / (QC(I,J,K)+QR(I,J,K))
                         ELSE
                            F_RAIN(I,J,K)=0.
                         ENDIF
                         IF (QG(I,J,K) .gt. EPSQ) THEN
                            F_RIMEF(I,J,K) = (5.*QG(I,J,K)      &
                                           +     QS(I,J,K))     &
                                           / (QS(I,J,K)         &
                                           +  QG(I,J,K))
                         ELSE
                            F_RIMEF(I,J,K)=1.
                         ENDIF
                         LIQW = MAX(0., (1.-F_ICE(I,J,K))*CWM(I,J,K))
                         QR(I,J,K) = LIQW*F_RAIN(I,J,K)*CWM(I,J,K)
                         QC(I,J,K) = LIQW*(1.-F_RAIN(I,J,K))*CWM(I,J,K)
                         IF (QG(I,J,K) .gt. EPSQ) THEN
                            FRACTION = MAX(0., MIN(QG(I,J,K)            &
                                       / (QG(I,J,K)+QS(I,J,K)), 1.) )
                         ELSE
                            FRACTION = 0.
                         ENDIF
                         QG(I,J,K) = FRACTION*F_ICE(I,J,K)*CWM(I,J,K)
                         QI(I,J,K) = 0.1*(1.-FRACTION)*F_ICE(I,J,K)*CWM(I,J,K)
                         QS(I,J,K) = 0.9*(1.-FRACTION)*F_ICE(I,J,K)*CWM(I,J,K)

                      ELSE       ! Below, the condensate is all new here
                         QC(I,J,K) = 0.0
                         QI(I,J,K) = 0.0
                         QR(I,J,K) = 0.0
                         QS(I,J,K) = 0.0
                         QG(I,J,K) = 0.0
                         IF (T(I,J,K) .le. 235.15) THEN
                            QI(I,J,K) = 0.5*CWM(I,J,K)
                            QS(I,J,K) = 0.5*CWM(I,J,K)
                         ELSEIF (T(I,J,K) .le. 258.15) THEN
                            QI(I,J,K) = 0.1*CWM(I,J,K)
                            QS(I,J,K) = 0.9*CWM(I,J,K)
                            deep_ice = .true.
                         ELSEIF (T(I,J,K) .le. 275.15) THEN
                            if (deep_ice .and. T(I,J,K).lt.273.15) then
                               QS(I,J,K) = CWM(I,J,K)
                            elseif (deep_ice .and. T(I,J,K).lt.274.15) then
                               QS(I,J,K) = 0.333*CWM(I,J,K)
                               QR(I,J,K) = 0.667*CWM(I,J,K)
                            elseif (deep_ice) then
                               QS(I,J,K) = 0.1*CWM(I,J,K)
                               QR(I,J,K) = 0.9*CWM(I,J,K)
                            else
                               QC(I,J,K) = CWM(I,J,K)
                            endif
                         ELSE
                            QC(I,J,K) = CWM(I,J,K)
                         ENDIF
                         LIQW = MAX(0., CWM(I,J,K) - QI(I,J,K)  &
                                                   - QS(I,J,K)  &
                                                   - QG(I,J,K))
                         IF (CWM(I,J,K) .gt. EPSQ) THEN
                            F_ICE(I,J,K) = (1.0-LIQW)/CWM(I,J,K)
                         ELSE
                            F_ICE(I,J,K) = 0.
                         ENDIF
                         IF (QR(I,J,K) .gt. EPSQ) THEN
                            F_RAIN(I,J,K) = QR(I,J,K)           &
                                    / (QC(I,J,K)+QR(I,J,K))
                         ELSE
                            F_RAIN(I,J,K)=0.
                         ENDIF
                         IF (QG(I,J,K) .gt. EPSQ) THEN
                            F_RIMEF(I,J,K) = (5.*QG(I,J,K)      &
                                           +     QS(I,J,K))     &
                                           / (QS(I,J,K)         &
                                           +  QG(I,J,K))
                         ELSE
                            F_RIMEF(I,J,K)=1.
                         ENDIF
                      ENDIF
                   ELSE
                      QC(I,J,K) = 0.0
                      QR(I,J,K) = 0.0
                      QI(I,J,K) = 0.0
                      QS(I,J,K) = 0.0
                      QG(I,J,K) = 0.0
                      F_ICE(I,J,K) = 0.0
                      F_RAIN(I,J,K) = 0.0
                      F_RIMEF(I,J,K) = 1.0
                   ENDIF
                ENDDO
             ENDDO
             ENDDO
          ENDIF  spec_adv_thompson

!
!----------------------------------------------------------------------
        CASE DEFAULT
!----------------------------------------------------------------------
!
          IF (CLD_INIT) THEN
            WRITE(0,*) 'Do nothing for default option'
          ENDIF
!
      END SELECT   ! MICROPHYSICS
!
!----------------------------------------------------------------------
!
      END SUBROUTINE UPDATE_WATER
!
!----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------

      SUBROUTINE TQADJUST(T,Q,QC,CWM,F_ICE,F_RAIN                       &
                         ,PRSI,PRSL                                     &
                         ,SPEC_ADV,RHgrd                                &
                         ,LM,IMS,IME)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    TQADJUST             TQADJUST
!   PRGRMMR: FERRIER         ORG: NP22     DATE: 5 APR 2016
!
! ABSTRACT:
!     Smooth temperature profiles when lapse rates exceed dry adiabatic
!     above PBL, prevent supersaturation with respect to water.
!
! PROGRAM HISTORY LOG (with changes to called routines) :
!   2016-04     FERRIER, JANJIC  - Smooth T profiles, prevent supersaturation
!
! USAGE: CALL TQADJUST FROM PHY_RUN
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!-----------------------------------------------------------------------
!
      USE MODULE_CONSTANTS,ONLY : CAPPA,CP,EP_2,EPSQ,R_d,R_v,CPV,CLIQ,  &
                                  A2,A4,PSAT,XLV,TIW
      USE MODULE_MP_ETANEW, ONLY : FPVS0
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
!
      INTEGER,PARAMETER :: JMS=1,JME=1
!
!----------------------
!-- Input argument variables
!----------------------
!
      INTEGER,INTENT(IN) :: LM,IMS,IME

      REAL,DIMENSION(IMS:IME,1:LM)  ,INTENT(IN) :: PRSL
      REAL,DIMENSION(IMS:IME,1:LM+1),INTENT(IN) :: PRSI
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) ::             &
                                                T,Q,QC,CWM,F_ICE,F_RAIN
      REAL,INTENT(IN) :: RHgrd
      LOGICAL,INTENT(IN) :: SPEC_ADV
!
!--  Local Variables
!
      INTEGER :: I,J,K,LM2,KMIX,KBOT,KTOP,ITmax,ITER,ITRmax
      REAL :: TK,PP,QV,QCW,ESW,QSW,DQsat,SSat,COND, &
              Qrain,Qice,Qliq
      REAL,DIMENSION(1:LM) :: Tcol,Pcol,QVcol,QCcol,EXNcol,THcol,DPcol, &
                              DTHcol,Fcol
      LOGICAL :: LRFilt,SSFilt
!
      REAL,PARAMETER :: SupSat=0.001, SubSat=-SupSat, DTHthresh=-0.01,  &
        TTP=TIW+0.01, XA=(CLIQ-CPV)/R_V, XB=XA+XLV/(R_V*TTP),           &
        XLV1=XLV/CP, XLV2=XLV1*XLV/R_V
!
!-----------------------------------------------------------------------
!
      ITmax=LM/5
      LM2=LM-2
!
!-----------------------------------------------------------------------
!--  Main loop through I, J,  ------------------------------------------
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
        DO I=IMS,IME
!
          LRFilt=.FALSE.       ! Lapse rate flag (full column)
          SSFilt=.FALSE.       ! Supersaturation flag
          IF(SPEC_ADV) THEN
            DO K=1,LM
              QCcol(K)=QC(I,J,K)
            ENDDO
          ELSE
            DO K=1,LM
              QCcol(K)=CWM(I,J,K)*(1.-F_ICE(I,J,K))*(1.-F_RAIN(I,J,K))
            ENDDO
          ENDIF
          DO K=1,LM
            Tcol(K)=T(I,J,K)
            Pcol(K)=PRSL(I,K)
            QVcol(K)=Q(I,J,K)/(1.-Q(I,J,K))        ! Water vapor mixing ratio
            EXNcol(K)=(1.E5/Pcol(K))**CAPPA
          ENDDO
!
!-----------------------------------------------------------------------
!-- Ferrier-Aligo condensation/evaporation algorithm - 1st of 2 times
!-----------------------------------------------------------------------
!
SSadj1:   DO K=1,LM
            TK=Tcol(K)                                     ! Temperature (deg K)
            PP=Pcol(K)                                     ! Pressure (Pa)
            QV=QVcol(K)                                    ! Water vapor mixing ratio
            QCW=QCcol(K)                                   ! Cloud water mixing ratio
            ESW=PSAT*EXP(A2*(TK-TTP)/(TK-A4))              ! Magnus Tetens
            ESW=MIN(ESW,0.99*PP)                           ! Saturation vapor pressure (water)
            QSW=RHgrd*EP_2*ESW/(PP-ESW)                    ! Saturation mixing ratio (water)
            DQsat=QV-QSW                                   ! Excess QV above saturation
            SSat=DQsat/QSW                                 ! Grid-scale supersaturation ratio
SSrem1:     IF(SSat>SupSat .OR.                     &      ! Remove supersaturation if SSat>0.1%
               (QCW>EPSQ .AND. SSat<SubSat) ) THEN         ! Adjust to saturation if SSat<0.1% w/ cloud water
              SSFilt=.TRUE.                                ! Supersaturation flag
cond_iter1:   DO ITER=1,10                                 ! Usually converges in <=3 iterations
                COND=DQsat/(1.+XLV2*QSW/(TK*TK))           ! Asai (1965, J. Japan)
                COND=MAX(COND, -QCW)                       ! Limit cloud water evaporation
                TK=TK+XLV1*COND                            ! Update temperature
                QV=QV-COND                                 ! Update water vapor mixing ratio
                QCW=QCW+COND                               ! Update cloud water mixing ratio
                ESW=PSAT*EXP(A2*(TK-TTP)/(TK-A4))          ! Magnus Tetens
                ESW=MIN(ESW,0.99*PP)                       ! Saturation vapor pressure (water)
                QSW=RHgrd*EP_2*ESW/(PP-ESW)                ! Water saturation mixing ratio
                DQsat=QV-QSW                               ! Excess QV above saturation
                SSat=DQsat/QSW                             ! Grid-scale supersaturation ratio
                IF (SSat>=SubSat .AND. SSat<=SupSat) EXIT  ! Exit if -0.1%<SSat<0.1%
                IF (SSat<SubSat .AND. QCW<=EPSQ)     EXIT  ! Exit if SSat<-0.1% & no cloud water
              ENDDO  cond_iter1                            ! 1st *cond*ensation *iter*ation
              IF (ITER<=10) THEN
                Tcol(K)=TK
                QVcol(K)=QV
                QCcol(K)=QCW
              ENDIF
            ENDIF  SSrem1
          ENDDO  SSadj1
!
          DO K=1,LM
            THcol(K)=Tcol(K)*EXNcol(K)
          ENDDO
!
          DTHcol(1)=1.
          DO K=2,LM
            DTHcol(K)=THcol(K-1)-THcol(K)
          ENDDO
!
          KMIX=0
          DO K=LM2,2,-1
            IF(DTHcol(K)>0.) THEN
!-- Start above the well-mixed layer immediately above the
!   surface where theta may decrease with height
              KMIX=K
              EXIT
            ENDIF
          ENDDO
!
!*************************
LRadjust: IF (KMIX>2) THEN
!*************************
!
            KTOP=0
            DO K=3,KMIX
              IF(DTHcol(K)<DTHthresh) THEN
                KTOP=K-1            !- Level at top of highest unstable layer
                EXIT
              ENDIF
            ENDDO
!
!-------------------------
Maybe_mix:  IF(KTOP>0) THEN
!-------------------------
!
              KBOT=0
              DO K=KMIX,2,-1
                IF(DTHcol(K)<DTHthresh) THEN
                  KBOT=K            !- Lowest unstable layer
                  EXIT
                ENDIF
              ENDDO
              IF(KBOT>0) THEN
                LRFilt=.TRUE.       !- For the full column (any layer)
                ITRmax=ITmax
              ELSE
                ITRmax=0            !- Do not mix
              ENDIF
!
              DO K=1,LM
                DPcol(K)=(PRSI(I,K+1)-PRSI(I,K)) ! Hydrostatic pressure thickness
                Fcol(K)=THcol(K)             !- Fcol, modified theta
              ENDDO
!
!- - - - - - - - - - - - -
Mix_lyrs:     DO ITER=1,ITRmax
!- - - - - - - - - - - - -
                DO K=KTOP,KBOT
                  IF(DTHcol(K)<DTHthresh) THEN
                    IF(DTHcol(K+1)<DTHthresh) THEN
!-- Mix 3 layers if current layer and layers above and below are unstable
                      Fcol(K)=(THcol(K-1)*DPcol(K-1)+THcol(K)*DPcol(K)    &
                               +THcol(K+1)*DPcol(K+1))/                   &
                              (DPcol(K-1)+DPcol(K)+DPcol(K+1))
                    ELSE
!-- Mix with higher layer if current layer is unstable
                      Fcol(K)=(THcol(K-1)*DPcol(K-1)+THcol(K)*DPcol(K))/  &
                              (DPcol(K-1)+DPcol(K))
                    ENDIF
                  ELSE IF(DTHcol(K+1)<DTHthresh) THEN
!-- Mix with lower layer if it is unstable
                    Fcol(K)=(THcol(K)*DPcol(K)+THcol(K+1)*DPcol(K+1))/    &
                            (DPcol(K)+DPcol(K+1))
                  ENDIF
!-- Do nothing if the current layer or the layer below is not unstable
                ENDDO
!
                DO K=KTOP,KBOT
                  THcol(K)=Fcol(K)
                ENDDO
                DO K=KTOP,KBOT
                  DTHcol(K)=THcol(K-1)-THcol(K)
                ENDDO
!
                KTOP=0
                DO K=3,KMIX
                  IF(DTHcol(K)<DTHthresh) THEN
                    KTOP=K-1         !- Level at top of highest unstable layer
                    EXIT
                  ENDIF
                ENDDO
!
                IF(KTOP<=0) EXIT Mix_lyrs   !- Exit with no unstable layer in column
!
                KBOT=0
                DO K=KMIX,2,-1
                  IF(DTHcol(K)<DTHthresh) THEN
                    KBOT=K           !- Lowest unstable layer
                    EXIT
                  ENDIF
                ENDDO
!- - - - - - - - - - - - -
              ENDDO  Mix_lyrs   !DO ITER
!- - - - - - - - - - - - -
              DO K=1,LM
                Tcol(K)=THcol(K)/EXNcol(K)     !-- Update T
              ENDDO
!-------------------------
            ENDIF  Maybe_mix    !IF(KTOP>0)
!-------------------------
!*************************
          ENDIF  LRadjust
!*************************
!
!-----------------------------------------------------------------------
!-- Ferrier-Aligo condensation/evaporation algorithm - 2nd of 2 times
!-----------------------------------------------------------------------
!
SSadj2:   DO K=1,LM
            TK=Tcol(K)                                     ! Temperature (deg K)
            PP=Pcol(K)                                     ! Pressure (Pa)
            QV=QVcol(K)                                    ! Water vapor mixing ratio
            QCW=QCcol(K)                                   ! Cloud water mixing ratio
!            TREF=TTP/TK                                    ! WSM6
!            ESW=PSAT*EXP(LOG(TREF)*(XA))*EXP(XB*(1.-TREF)) ! WSM6
!            ESW=1000.*FPVS0(TK)                            ! Old global tables
            ESW=PSAT*EXP(A2*(TK-TTP)/(TK-A4))              ! Magnus Tetens
!            TREF=TK-TIW                                    ! Bolton (1980)
!           ESW=611.2*EXP(17.67*TREF/(TREF+243.5))          ! Bolton (1980)
            ESW=MIN(ESW,0.99*PP)                           ! Saturation vapor pressure (water)
            QSW=RHgrd*EP_2*ESW/(PP-ESW)                    ! Saturation mixing ratio (water)
            DQsat=QV-QSW                                   ! Excess QV above saturation
            SSat=DQsat/QSW                                 ! Grid-scale supersaturation ratio
SSrem2:     IF(SSat>SupSat .OR.                     &      ! Remove supersaturation if SSat>0.1%
               (QCW>EPSQ .AND. SSat<SubSat) ) THEN         ! Adjust to saturation if SSat<0.1% w/ cloud water
              SSFilt=.TRUE.                                ! Supersaturation flag
cond_iter2:   DO ITER=1,10                                 ! Usually converges in <=3 iterations
                COND=DQsat/(1.+XLV2*QSW/(TK*TK))           ! Asai (1965, J. Japan)
                COND=MAX(COND, -QCW)                       ! Limit cloud water evaporation
                TK=TK+XLV1*COND                            ! Update temperature
                QV=QV-COND                                 ! Update water vapor mixing ratio
                QCW=QCW+COND                               ! Update cloud water mixing ratio
!                TREF=TTP/TK                                ! WSM6
!                ESW=PSAT*EXP(LOG(TREF)*(XA))*EXP(XB*(1.-TREF)) ! WSM6
!                ESW=1000.*FPVS0(TK)                        ! Old global tables
                ESW=PSAT*EXP(A2*(TK-TTP)/(TK-A4))          ! Magnus Tetens
!                TREF=TK-TIW                                ! Bolton (1980)
!                ESW=611.2*EXP(17.67*TREF/(TREF+243.5))     ! Bolton (1980)
                ESW=MIN(ESW,0.99*PP)                       ! Saturation vapor pressure (water)
                QSW=RHgrd*EP_2*ESW/(PP-ESW)                ! Water saturation mixing ratio
                DQsat=QV-QSW                               ! Excess QV above saturation
                SSat=DQsat/QSW                             ! Grid-scale supersaturation ratio
                IF (SSat>=SubSat .AND. SSat<=SupSat) EXIT  ! Exit if -0.1%<SSat<0.1%
                IF (SSat<SubSat .AND. QCW<=EPSQ)     EXIT  ! Exit if SSat<-0.1% & no cloud water
              ENDDO  cond_iter2                            ! 2nd *cond*ensation *iter*ation
              IF (ITER<=10) THEN
                Tcol(K)=TK
                QVcol(K)=QV
                QCcol(K)=QCW
              ENDIF
            ENDIF  SSrem2
          ENDDO  SSadj2
!
!#######################################################################
!-- Update 3D arrays
!#######################################################################
!
adjust1:  IF (LRFilt .OR. SSFilt) THEN
            DO K=1,LM
              T(I,J,K)=Tcol(K)    !- Update T
            ENDDO
          ENDIF  adjust1
!
adjust2:  IF (SSFilt) THEN
            DO K=1,LM             !- Update Q
              Q(I,J,K)=QVcol(K)/(1.+QVcol(K))
            ENDDO
            IF(SPEC_ADV) THEN
              DO K=1,LM           !- Update QC
                QC(I,J,K)=QCcol(K)
              ENDDO
            ELSE
              DO K=1,LM           !- Update CWM, F_ICE, F_RAIN
                Qrain=CWM(I,J,K)*(1.-F_ICE(I,J,K))*F_RAIN(I,J,K)
                Qliq=QCcol(K)+Qrain
                Qice=CWM(I,J,K)*F_ICE(I,J,K)
                CWM(I,J,K)=Qliq+Qice
                IF(CWM(I,J,K)>EPSQ) F_ICE(I,J,K)=Qice/CWM(I,J,K)
                IF(Qliq>EPSQ) F_RAIN(I,J,K)=Qrain/Qliq
              ENDDO
            ENDIF
          ENDIF  adjust2
!
        ENDDO   !- I
      ENDDO     !- J
!-----------------------------------------------------------------------
!
      END SUBROUTINE TQADJUST
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!

      END MODULE MODULE_MICROPHYSICS_NMM

!-----------------------------------------------------------------------
