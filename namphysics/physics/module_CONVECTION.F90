!-----------------------------------------------------------------------
!
      MODULE MODULE_CONVECTION
!
!-----------------------------------------------------------------------
!
!***  THE CONVECTION DRIVERS AND PACKAGES
!
!-----------------------------------------------------------------------
!
      USE MODULE_KINDS
!
      USE MODULE_CONSTANTS ,ONLY : CAPPA,EPSQ,G,R_D
      USE MODULE_CU_BMJ
      USE MODULE_CU_SAS
      USE MODULE_CU_SASHUR
      USE MODULE_CU_SCALE
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: CUCNVC
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE CUCNVC(NTSD,DT,NCNVC,NRADS,NRADL,MINUTES_HISTORY       &
                       ,ENTRAIN,NEWALL,NEWSWAP,NEWUPUP,NODEEP           &
                       ,FRES,FR,FSL,FSS                                 &
                       ,CLDEFI                                          &
                       ,U_PHY,V_PHY                                     &
                       ,F_ICE,F_RAIN                                    &
                       ,QC,QR,QI,QS,QG                                  &
                       ,F_QC,F_QR,F_QI,F_QS,F_QG                        &
                       ,PHINT,PHMID                                     &
                       ,area                                            &
                       ,T,Q,CWM,TCUCN                                   &
                       ,VVL                                             &
                       ,FIS                                             &
                       ,PREC,ACPREC,CUPREC,CUPPT,CPRATE                 &
                       ,CNVBOT,CNVTOP,SM                                &
                       ,HTOP,HTOPD,HTOPS                                &
                       ,HBOT,HBOTD,HBOTS                                &
                       ,AVCNVC,ACUTIM                                   &
                       ,RSWIN,RSWOUT                                    &
                       ,CONVECTION,MICROPHYSICS                         &
                       ,SICE,QWBS,TWBS,PBLH,DTDT_PHY,DUDT_PHY,DVDT_PHY  &
                       ,MOMMIX,PGCON,SAS_MASS_FLUX                      &
                       ,SHALCONV,SHAL_PGCON                             &
                       ,IMS,IME,LM                                      &
                                                    )
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    CUCNVC      CONVECTIVE PRECIPITATION OUTER DRIVER
!   PRGRMMR: BLACK           ORG: W/NP22     DATE: 02-03-21
!
! ABSTRACT:
!     CUCVNC DRIVES THE WRF CONVECTION SCHEMES
!
! PROGRAM HISTORY LOG:
!   02-03-21  BLACK      - ORIGINATOR
!   04-11-18  BLACK      - THREADED
!   06-10-11  BLACK      - BUILT INTO UMO PHYSICS COMPONENT
!   08-08     JANJIC     - Synchronize WATER array and Q.
!   10-10-26  WEIGUO WANG - add GFS SAS convection
!   14-06-19  WEIGUO WANG - add hurricane SAS (moved from hwrf)
!   16-08-29  WEIGUO WANG - add scale-aware convection schemes
! USAGE: CALL CUCNVC FROM PHY_RUN
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
      integer(kind=kint),parameter:: jms=1,jme=1
!-----------------------------------------------------------------------

      character(99),intent(in):: &
       convection,microphysics
!
      logical(kind=klog),intent(in):: &
       entrain,newall,newswap,newupup,nodeep &
      ,f_qc,f_qr,f_qi,f_qs,f_qg
!
      integer(kind=kint),intent(in):: &
       ims,ime,lm &
      ,ncnvc,minutes_history &
      ,nrads,nradl,ntsd
!
      real(kind=kfpt),intent(in):: &
       dt,fres,fr,fsl,fss
!
      real(kind=kfpt),intent(inout):: &
       acutim,avcnvc
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
       fis &
      ,rswin,rswout,sm,area
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout):: &
       acprec,cldefi &
      ,cnvbot,cnvtop &
      ,cuppt,cuprec &
      ,hbot,htop &
      ,hbotd,htopd &
      ,hbots,htops &
      ,prec,cprate
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
       sice,qwbs,twbs,pblh  !fOR SAS

      REAL(kind=kfpt), OPTIONAL, INTENT(IN) :: &
              PGCON,sas_mass_flux,shal_pgcon,mommix,shalconv         !sashur
!
      real(kind=kfpt),dimension(ims:ime,1:lm)  ,intent(in):: &
       vvl,phmid
!
      real(kind=kfpt),dimension(ims:ime,1:lm+1),intent(in):: &
       phint
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out):: &
       dudt_phy,dvdt_phy,dtdt_phy
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
       q,t &
      ,f_ice &
      ,f_rain &
      ,cwm &
      ,u_phy,v_phy &
      ,tcucn,QC,QI,QR,QS,QG
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      logical(kind=klog):: &
       warm_rain,F_QGr
!
      logical(kind=klog),dimension(ims:ime,jms:jme):: &
       cu_act_flag
!
      integer(kind=kint):: &
       i,j &
      ,k &
      ,mnto &
      ,n,ncubot,ncutop,n_timstps_output
!
      integer(kind=kint),dimension(ims:ime,jms:jme):: &
       LBOT,LTOP
!
      real(kind=kfpt):: &
       cf_hi,dtcnvc,fice,frain,g_inv &
      ,pcpcol,ql,ql_k,rdtcnvc &
      ,QCW,QCI,QRain,QSnow,QGraup &
      ,tl
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME):: &
       CUBOT,CUTOP,NCA &
      ,RAINC,RAINCV,SFCZ,XLAND
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:LM):: &
       DZ,exner &
      ,th,rr &
      ,RQCCUTEN,RQRCUTEN &
      ,RQICUTEN,RQSCUTEN &
      ,RQCUTEN,RTHCUTEN &
      ,RQGCUTEN
!-----------------------------------------------------------------------
      REAL(kind=kfpt) :: TCHANGE
!-----------------------------------------------------------------------
!***********************************************************************
!
!-----------------------------------------------------------------------
!***  RESET THE HBOT/HTOP CONVECTIVE CLOUD BOTTOM (BASE) AND TOP ARRAYS
!***  USED IN RADIATION.  THEY STORE THE MAXIMUM VERTICAL LIMITS OF
!***  CONVECTIVE CLOUD BETWEEN RADIATION CALLS.  THESE ARRAYS ARE OUT
!***  OF THE WRF PHYSICS AND THUS THEIR VALUES INCREASE UPWARD.
!***  CUPPT IS THE ACCUMULATED CONVECTIVE PRECIPITATION BETWEEN
!***  RADIATION CALLS.
!-----------------------------------------------------------------------
!
      IF(MOD(NTSD,NRADS)==0.OR.MOD(NTSD,NRADL)==0)THEN
         DO J=JMS,JME
         DO I=IMS,IME
           HTOP(I,J)=0.
           HBOT(I,J)=REAL(LM+1)
           CUPPT(I,J)=0.
         ENDDO
         ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
      IF(MOD(NTSD,NCNVC)/=0.AND.CONVECTION=='bmj')RETURN
      IF(MOD(NTSD,NCNVC)/=0.AND.CONVECTION=='sas')RETURN
      IF(MOD(NTSD,NCNVC)/=0.AND.CONVECTION=='sashur')RETURN
      IF(MOD(NTSD,NCNVC)/=0.AND.CONVECTION=='scalecu')RETURN
!-----------------------------------------------------------------------
!
      IF(MICROPHYSICS=='fer' .OR. MICROPHYSICS=='fer_hires') THEN
         F_QGr=.FALSE.
      ELSE
         F_QGr=F_QG
      ENDIF
!
!-----------------------------------------------------------------------
!***  GENERAL PREPARATION
!-----------------------------------------------------------------------
!
      AVCNVC=AVCNVC+1.
      ACUTIM=ACUTIM+1.
!
      DTCNVC=NCNVC*DT
      RDTCNVC=1./DTCNVC
      G_INV=1./G
!
!.......................................................................
!zj$omp parallel do &
!zj$omp& private(j,i,k,ql,tl)
!.......................................................................
      DO J=JMS,JME
      DO I=IMS,IME
!
        RAINCV(I,J)=0.
        RAINC(I,J)=0.
        XLAND(I,J)=SM(I,J)+1.
        NCA(I,J)=0.
        SFCZ(I,J)=FIS(I,J)*G_INV
!
        CUTOP(I,J)=999.
        CUBOT(I,J)=999.
!
!-----------------------------------------------------------------------
!***  FILL VERTICAL WORKING ARRAYS.
!-----------------------------------------------------------------------
!
        DO K=1,LM
!
          QL=MAX(Q(I,J,K),EPSQ)
          TL=T(I,J,K)
          RR(I,J,K)=PHMID(I,K)/(R_D*TL*(.608*ql+1.))
          T(I,J,K)=TL
!
          EXNER(I,J,K)=(PHMID(I,K)*1.E-5)**CAPPA
          TH(I,J,K)=TL/EXNER(I,J,K)
!
        ENDDO
      ENDDO
      ENDDO
!.......................................................................
!zj$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  Compute velocity components at mass points.
!-----------------------------------------------------------------------
!
!.......................................................................
!zj$omp parallel do &
!zj$omp& private(j,i,k)
!.......................................................................
      do k=1,lm
        do j=jms,jme
          do i=ims,ime
            RTHCUTEN(I,J,K)=0.
            RQCUTEN(I,J,K)=0.
            RQCCUTEN(I,J,K)=0.
            RQRCUTEN(I,J,K)=0.
            RQICUTEN(I,J,K)=0.
            RQSCUTEN(I,J,K)=0.
            RQGCUTEN(I,J,K)=0.
            dtdt_phy(i,j,k)=0.
            dudt_phy(i,j,k)=0.
            dvdt_phy(i,j,k)=0.
          enddo
        enddo
      enddo
!.......................................................................
!zj$omp end parallel do
!.......................................................................
!-----------------------------------------------------------------------
!.......................................................................
!zj$omp parallel do                                                       &
!zj$omp private(i,j,k,ql_k)
!.......................................................................
      DO J=JMS,JME
        DO I=IMS,IME
          DZ(I,J,LM)=T(I,J,LM)*(.608*Q(I,J,LM)+1.)*R_D &
                    *(PHINT(I,LM+1)-PHINT(I,LM)) &
                    /(PHMID(I,LM)*G)
        ENDDO
!
        DO K=LM-1,1,-1
        DO I=IMS,IME
          QL_K=MAX(Q(I,J,K),EPSQ)
          DZ(I,J,K)=T(I,J,K)*(.608*QL_K+1.)*R_D &
                    *(PHINT(I,K+1)-PHINT(I,K)) &
                    /(PHMID(I,K)*G)
        ENDDO
        ENDDO
!
      ENDDO
!
!-----------------------------------------------------------------------
!
!***  SINGLE-COLUMN CONVECTION
!
!-----------------------------------------------------------------------
!
      IF(trim(CONVECTION) == 'bmj') then
!
        call bmjdrv( &
                     ims,ime,jms,jme,lm &
                    ,entrain,newall,newswap,newupup,nodeep &
                    ,fres,fr,fsl,fss &
                    ,dt,ntsd,ncnvc &
                    ,raincv,cutop,cubot &
                    ,th,t,q,u_phy,v_phy,dudt_phy,dvdt_phy &
                    ,phint,phmid,exner &
                    ,cldefi,xland,cu_act_flag &
                  ! optional
                    ,rthcuten,rqcuten &
                    )
!
      ELSEIF(trim(CONVECTION) == 'sas') then
!
       call sasdrv( &
                   ims,ime,jms,jme,lm &
                  ,dt,ntsd,ncnvc &
                  ,th,t,sice,twbs,qwbs,pblh,u_phy,v_phy &
                  ,q,qc,qr,qi,qs,qg &
                  ,f_qc,f_qr,f_qi,f_qs,F_QGr &
                  ,phint,phmid,exner,rr,dz &
                  ,xland,cu_act_flag &
                  ,vvl &
                  ,raincv,cutop,cubot &
                  ,dudt_phy,dvdt_phy &
                  ! optional
                  ,rthcuten, rqcuten &
                  ,rqccuten, rqrcuten &
                  ,rqicuten, rqscuten &
                  ,rqgcuten  &
                  )
!
      ELSEIF(trim(CONVECTION) == 'sashur') then
!
       call sasdrv_hur( &
                   ims,ime,jms,jme,lm &
                  ,dt,ntsd,ncnvc &
                  ,th,t,sice,vvl,twbs,qwbs,pblh,u_phy,v_phy &
                  ,q,qc,qr,qi,qs,qg &
                  ,f_qc,f_qr,f_qi,f_qs,F_QGr &
                  ,phint,phmid,exner,rr,dz &
                  ,xland,cu_act_flag &
                  ,MOMMIX,PGCON,SAS_MASS_FLUX   &
                  ,SHALCONV,SHAL_PGCON          &
                  ,raincv,cutop,cubot &
                  ,dudt_phy,dvdt_phy &
                  ! optional
                  ,rthcuten, rqcuten &
                  ,rqccuten, rqrcuten &
                  ,rqicuten, rqscuten &
                  ,rqgcuten  &
                  )
!
      ELSEIF(trim(CONVECTION) == 'scalecu') then
!
       call scalecudrv( &
                   ims,ime,jms,jme,lm &
                  ,dt,ntsd,ncnvc &
                  ,th,t,sice,vvl,twbs,qwbs,pblh,u_phy,v_phy &
                  ,q,qc,qr,qi,qs,qg &
                  ,f_qc,f_qr,f_qi,f_qs,F_QGr &
                  ,phint,phmid,exner,rr,dz &
                  ,xland,cu_act_flag &
                  ,area               &
                  ,MOMMIX,PGCON,SAS_MASS_FLUX   &
                  ,SHALCONV,SHAL_PGCON          &
                  ,raincv,cutop,cubot &
                  ,dudt_phy,dvdt_phy &
                  ! optional
                  ,rthcuten, rqcuten &
                  ,rqccuten, rqrcuten &
                  ,rqicuten, rqscuten &
                  ,rqgcuten  &
                  )
!
      END IF
!
!-----------------------------------------------------------------------
!
!***  CNVTOP/CNVBOT HOLD THE MAXIMUM VERTICAL LIMITS OF CONVECTIVE CLOUD
!***  BETWEEN HISTORY OUTPUT TIMES.  HBOTS/HTOPS STORE SIMILIAR INFORMATION
!***  FOR SHALLOW (NONPRECIPITATING) CONVECTION, AND HBOTD/HTOPD ARE FOR
!***  DEEP (PRECIPITATING) CONVECTION.
!
      CF_HI=REAL(MINUTES_HISTORY)/60.
      N_TIMSTPS_OUTPUT=NINT(3600.*CF_HI/DT)
      MNTO=MOD(NTSD,N_TIMSTPS_OUTPUT)
!
      IF(MNTO>0.AND.MNTO<=NCNVC)THEN
        DO J=JMS,JME
        DO I=IMS,IME
          CNVBOT(I,J)=REAL(LM+1.)
          CNVTOP(I,J)=0.
          HBOTD(I,J)=REAL(LM+1.)
          HTOPD(I,J)=0.
          HBOTS(I,J)=REAL(LM+1.)
          HTOPS(I,J)=0.
        ENDDO
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!.......................................................................
!zj$omp parallel do                                                       &
!zj$omp& private(j,k,i,tchange,pcpcol,ncubot,ncutop,QCW,QRain,QCI,QSnow,QGraup)
!.......................................................................
!-----------------------------------------------------------------------
      do j=jms,jme
      do i=ims,ime
!-----------------------------------------------------------------------
!
!***  UPDATE TEMPERATURE, SPECIFIC HUMIDITY, AND HEATING.
!
        DO K=1,LM
!
!rv  ------------  update is now in the SOLVER ------------
          DTDT_PHY(I,J,K)=RTHCUTEN(I,J,K)*exner(I,J,K)
!rv       T(I,J,K)=T(I,J,K)+DTDT_PHY(I,J,K)*DTCNVC
          Q(I,J,K)=Q(I,J,K)+RQCUTEN(I,J,K)*DTCNVC
!rv       U_PHY(I,J,K)=U_PHY(I,J,K)+DUDT_PHY(I,J,K)*DTCNVC
!rv       V_PHY(I,J,K)=V_PHY(I,J,K)+DVDT_PHY(I,J,K)*DTCNVC
          TCUCN(I,J,K)=TCUCN(I,J,K)+DTDT_PHY(I,J,K)

sas_test: IF(CONVECTION=='sas') THEN
            QC(I,J,K)=QC(I,J,K)+DTCNVC*RQCCUTEN(I,J,K)
            QCW=QC(I,J,K)
            QRain=0.
            QCI=0.
            QSnow=0.
            QGraup=0.
            IF(F_QR) THEN
              QR(I,J,K)=QR(I,J,K)+DTCNVC*RQRCUTEN(I,J,K)
              QRain=QR(I,J,K)
            ENDIF
            IF(F_QI) THEN
              QI(I,J,K)=QI(I,J,K)+DTCNVC*RQICUTEN(I,J,K)
              QCI=QI(I,J,K)
            ENDIF
            IF(F_QS) THEN
              QS(I,J,K)=QS(I,J,K)+DTCNVC*RQSCUTEN(I,J,K)
              QSnow=QS(I,J,K)
            ENDIF
            IF(F_QGr) THEN
              QG(I,J,K)=QG(I,J,K)+DTCNVC*RQGCUTEN(I,J,K)
              QGraup=QG(I,J,K)
            ENDIF
!-- Couple CWM, F_ice, & F_rain arrays
            CWM(I,J,K)=QCW+QRain+QCI+QSnow+QGraup
            F_ICE(I,J,K)=0.
            F_RAIN(I,J,K)=0.
            IF(CWM(I,J,K)>EPSQ) F_ICE(I,J,K)=(QCI+QSnow+QGraup)/CWM(I,J,K)
            IF(QRain>EPSQ) F_RAIN(I,J,K)=QRain/(QCW+QRain)
          ENDIF  sas_test
!
        ENDDO

!
!***  UPDATE PRECIPITATION
!
        PCPCOL=RAINCV(I,J)*1.E-3*NCNVC
        PREC(I,J)=PREC(I,J)+PCPCOL
        ACPREC(I,J)=ACPREC(I,J)+PCPCOL
        CUPREC(I,J)=CUPREC(I,J)+PCPCOL
        CUPPT(I,J)=CUPPT(I,J)+PCPCOL
        CPRATE(I,J)=PCPCOL
!
!***  SAVE CLOUD TOP AND BOTTOM FOR RADIATION (HTOP/HBOT) AND
!***  FOR OUTPUT (CNVTOP/CNVBOT, HTOPS/HBOTS, HTOPD/HBOTD) ARRAYS.
!***  MUST BE TREATED SEPARATELY FROM EACH OTHER.
!
        NCUTOP=NINT(CUTOP(I,J))
        NCUBOT=NINT(CUBOT(I,J))
!
        IF(NCUTOP>1.AND.NCUTOP<LM+1)THEN
          HTOP(I,J)=MAX(CUTOP(I,J),HTOP(I,J))
          CNVTOP(I,J)=MAX(CUTOP(I,J),CNVTOP(I,J))
          IF(PCPCOL>0.)THEN
            HTOPD(I,J)=MAX(CUTOP(I,J),HTOPD(I,J))
          ELSE
            HTOPS(I,J)=MAX(CUTOP(I,J),HTOPS(I,J))
          ENDIF
        ENDIF
        IF(NCUBOT>0.AND.NCUBOT<LM+1)THEN
          HBOT(I,J)=MIN(CUBOT(I,J),HBOT(I,J))
          CNVBOT(I,J)=MIN(CUBOT(I,J),CNVBOT(I,J))
          IF(PCPCOL>0.)THEN
            HBOTD(I,J)=MIN(CUBOT(I,J),HBOTD(I,J))
          ELSE
            HBOTS(I,J)=MIN(CUBOT(I,J),HBOTS(I,J))
          ENDIF
        ENDIF
!
      ENDDO
      ENDDO
!.......................................................................
!zj$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CUCNVC
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_CONVECTION
!
!-----------------------------------------------------------------------
