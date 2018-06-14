!-----------------------------------------------------------------------
!
      MODULE MODULE_CU_SAS
!
!     12-10-2010  Created by Weiguo Wang
!-----------------------------------------------------------------------
!
!***  THE CONVECTION DRIVERS AND PACKAGES
!
!-----------------------------------------------------------------------
!
      USE MODULE_KINDS
!
      USE MODULE_CONSTANTS,ONLY : g99 => g, CP, ELWV,EPSQ
      use machine , only : kind_phys
      use funcphys , only : fpvs, gpvs

!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
      REAL,    PARAMETER ::    XLV=ELWV
!
      PUBLIC :: SASDRV 
      PUBLIC :: SAS_INIT
!
!-----------------------------------------------------------------------
       CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE SASDRV( &
                        IMS,IME,JMS,JME,LM &
                       ,DT,NTSD,NCNVC &
                       ,TH,T,SICE,SHEAT,LHEAT,PBLH,U,V &
                       ,Q,QC,QR,QI,QS,QG &
                       ,F_QC,F_QR,F_QI,F_QS,F_QG &
                       ,PHINT,PHMID,exner,RR,DZ &
                       ,XLAND,CU_ACT_FLAG &
                       ,VVL &
                       ,RAINCV,CUTOP,CUBOT &   !! out below
                       ,DUDT,DVDT &
                      ! optional
                       ,RTHCUTEN,RQCUTEN &
                       ,RQCCUTEN,RQRCUTEN &
                       ,RQICUTEN,RQSCUTEN &
                       ,RQGCUTEN &
                       )
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER(kind=kint),INTENT(IN):: &
       IMS,IME,JMS,JME,LM
!
      INTEGER(kind=kint),INTENT(IN) :: ntsd,NCNVC
      REAL(kind=kfpt),   INTENT(IN) :: DT
!
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(IN):: &
       XLAND,SICE,PBLH,SHEAT,LHEAT
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:lm),INTENT(IN):: &
       dz,exner,phmid,rr,t,th,U,V 
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:lm),INTENT(IN):: VVL
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:lm+1),INTENT(IN):: &
       phint
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:lm),INTENT(IN):: Q
      REAL(kind=kfpt),DIMENSION(:,:,:),INTENT(IN):: QC,QR,QI,QS,QG
      LOGICAL,INTENT(IN) :: F_QC,F_QR,F_QI,F_QS,F_QG
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:lm),optional,intent(inout):: &
       RQCUTEN,RTHCUTEN &
      ,RQCCUTEN,RQRCUTEN &
      ,RQSCUTEN,RQICUTEN &
      ,RQGCUTEN 
! 
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT):: &
       RAINCV
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(OUT):: &
       CUBOT,CUTOP
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:lm),INTENT(OUT):: &
       DUDT,DVDT 
!
      LOGICAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT):: &
       CU_ACT_FLAG
!
      LOGICAL DEEP, SHALLOW
!
!-----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
!-----------------------------------------------------------------------
!
      INTEGER :: I,J,K,ICLDCK,KFLIP,idbg,jdbg

! For SAS
      INTEGER :: KM,NUM_ICE,NSHAL,NDEEP
      INTEGER, PARAMETER :: IX=1, IM=1, ncloud=1
      INTEGER :: jcap, kcnv(IX), KBOT(IX), KTOP(IX), islimsk(IX)
      REAL(kind=kind_phys), DIMENSION(IX,lm) :: delp, prsl,phil,q1,t1,u1,v1,VVEL,     &
                                ud_mf,dd_mf,dt_mf, q0,t0,u0,v0,cnvc,cnvw
      REAL(kind=kind_phys), DIMENSION(IX) :: psp,cldwrk,rn,hpbl,hflx,evap
      REAL(kind=kind_phys), DIMENSION(IX,lm,2) :: CLW, CLW0  !! 1-ice  2-liquid 
      REAL(kind=kind_phys) :: TMP, DELT, RDELT, landmask
      REAL(kind=kind_phys), PARAMETER :: H1=1., H0=0.,    &
                         mommix=1.0    !HWRF uses this to adjust/tune moment mixing
      REAL, DIMENSION(lm+1)    :: ZF
      LOGICAL, PARAMETER :: LPR=.FALSE.  !- Set to .TRUE. for debugging
      LOGICAL :: MULTI_ICE
!
!------------------------------------------------------------------------
!
      KM = LM
!
      jcap = 126
!      print *,'in module_CU_SAS, sasdriver,jcap=',jcap
!
      NUM_ICE=0
      IF(F_QI) NUM_ICE=1
      IF(F_QS) NUM_ICE=NUM_ICE+1
      IF(F_QG) NUM_ICE=NUM_ICE+1
      IF(NUM_ICE>1) THEN
        MULTI_ICE=.TRUE.
      ELSE
        MULTI_ICE=.FALSE.
      ENDIF
!.......................................................................
!$omp parallel do                &
!$omp     private(k,j,i)
!.......................................................................
      DO K=1,lm
        DO J=JMS,JME
          DO I=IMS,IME
            DUDT(I,J,K) = 0.0
            DVDT(I,J,K) = 0.0
          ENDDO
        ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do              
!.......................................................................
!.......................................................................
!$omp parallel do                &
!$omp     private(k,j,i)
!.......................................................................
      DO K=1,lm
        DO J=JMS,JME
          DO I=IMS,IME
            RTHCUTEN(I,J,K) = 0.0
            RQCUTEN(I,J,K) = 0.0
            RQCCUTEN(I,J,K) = 0.0
            IF(F_QR) RQRCUTEN(I,J,K) = 0.0
            IF(F_QI) RQICUTEN(I,J,K) = 0.0
            IF(F_QS) RQSCUTEN(I,J,K) = 0.0
            IF(F_QG) RQGCUTEN(I,J,K) = 0.0
          ENDDO
        ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do                
!.......................................................................
!
      DELT=DT*NCNVC
      RDELT=1./DELT
dbg1: IF(LPR) THEN
        write(0,*)'delt,rdelt=',delt,rdelt
        idbg=(ims+ime+1)/2   !- or set to fixed "I"
        jdbg=(jms+jme+1)/2   !- or set to fixed "J"
        NSHAL=0
        NDEEP=0
      ENDIF dbg1
!
!-----------------------------------------------------------------------
!
!***  PREPARE TO CALL SAS CONVECTION SCHEME
!
!-----------------------------------------------------------------------
!
!***  CHECK TO SEE IF THIS IS A CONVECTION TIMESTEP
!                                                                        
      ICLDCK=MOD(ntsd,NCNVC)                                              
      IF(ICLDCK/=0) RETURN
!
!-----------------------------------------------------------------------
!                                                                      
!***  COMPUTE CONVECTION EVERY NCNVC*DT/60.0 MINUTES
!                                                                     
      DO J=JMS,JME
        DO I=IMS,IME
          CU_ACT_FLAG(I,J)=.TRUE.
        ENDDO
      ENDDO
!
!
!.......................................................................
!$omp parallel do                &
!$omp     private(j,i,k,landmask,islimsk,zf,kflip,psp,prsl,delp,phil,u1,    &
!$omp             v1,t1,q1,clw,ud_mf,dd_mf,dt_mf,cldwrk,vvel,hflx,evap,hpbl,&
!$omp             kcnv,kbot,ktop,u0,v0,t0,q0,clw0,cnvc,cnvw,tmp,rn,jcap)
!.......................................................................
      DO J=JMS,JME  
        DO I=IMS,IME
          RAINCV(I,J)=0.
!
!***  CONVERT TO BMJ LAND MASK (1.0 FOR SEA; 0.0 FOR LAND)
!
          LANDMASK=XLAND(I,J)-1.
          ISLIMSK(1) = 1. - LANDMASK
          IF(SICE(I,J) > 0.5) ISLIMSK(1) = 2     !! 0-sea; 1-land; 2-ice 
!
!***  FILL 1-D VERTICAL ARRAYS 
!
          ZF(1) = 0.0
          DO K=2,LM+1 
            KFLIP = LM + 1 + 1 -K
            ZF(K) = ZF(K-1) + DZ(I,J,KFLIP)
          ENDDO
          PSP(1) = PHINT(I,J,lm+1)        ! Surface pressure, Pa
vloop1:   DO K=1,lm
            kflip = LM + 1 -K
            prsl(1,K) = phmid(I,J,KFLIP)
            delp(1,K) = RR(I,J,KFLIP)*g99*DZ(I,J,KFLIP) 
            phil(1,K) = 0.5*(ZF(K) + ZF(K+1) )*g99              
            u1(1,K) = U(I,J,KFLIP)
            v1(1,K) = V(I,J,KFLIP)
            t1(1,K) = T(I,J,KFLIP)
            q1(1,K) = MAX(EPSQ,Q(I,J,KFLIP)) 
            clw(1,K,1) = 0.0
            if (f_qc) clw(1,K,1) = QC(I,J,KFLIP)
            if (f_qr) clw(1,K,1) = clw(1,K,1) + QR(I,J,KFLIP)
            clw(1,K,2) = 0.0
            if (f_qi) clw(1,K,2) = QI(I,J,KFLIP)
            if (f_qs) clw(1,K,2) = clw(1,K,2) + QS(I,J,KFLIP)
            if (f_qg) clw(1,K,2) = clw(1,K,2) + QG(I,J,KFLIP)
            ud_mf(1,K) = 0.0
            dd_mf(1,K) = 0.0
            dt_mf(1,K) = 0.0 
            cnvc(1,K) = 0.  !-- convective cloud cover (new, not yet used)
            cnvw(1,K) = 0.  !-- convective cloud water (new, not yet used)
            cldwrk(1) = 0.0
            VVEL(1,K) = 0.
            if(kflip-1 <= lm-1 .and. kflip-1 >= 1 )     &
              VVEL(1,K)=VVL(I,J,KFLIP-1)
          ENDDO  vloop1
          hflx(1) = SHEAT(I,J)/RR(I,J,LM)/CP            ! W/m2 to K m/s
          evap(1) = LHEAT(I,J)/RR(I,J,LM)/XLV
          hpbl(1) = PBLH(I,J)
          KCNV(1) = 0     
          KBOT(1) = KM 
          KTOP(1) = 1       
          u0 = u1
          v0 = v1
          t0 = t1
          q0 = q1
          clw0 = clw   
!
!---  CALL CONVECTION
!
!          print *,'in module_CU_SAS, sasdriver, call sascnvn,jcap=',jcap
          CALL sascnvn(im,ix,km,jcap,delt,delp,prsl,psp,phil,clw,       &
               q1,t1,u1,v1,cldwrk,rn,kbot,ktop,kcnv,islimsk,            &
               VVEL,ncloud,ud_mf,dd_mf,dt_mf,cnvc,cnvw)
          IF(KCNV(1)>0) THEN
            DEEP=.TRUE.
            SHALLOW=.FALSE.
          ELSE
            DEEP=.FALSE.
            SHALLOW=.TRUE.
          ENDIF
!
          IF(SHALLOW) THEN
!            print *,'call shalcnv,jcap=',jcap,'delt=',delt
            CALL shalcnv(im,ix,km,jcap,delt,delp,prsl,psp,phil,clw,     &
                q1,t1,u1,v1,rn,kbot,ktop,kcnv,islimsk,                  &
                VVEL,ncloud,hpbl,hflx,evap,ud_mf,dt_mf,cnvc,cnvw)
            IF(KTOP(1)<1) SHALLOW=.FALSE.
          ENDIF
!
          CUTOP(I,J) = REAL(KTOP(1))
          CUBOT(I,J) = REAL(KBOT(1))
          RAINCV(I,J) = RN(1)*1.E3/NCNVC
!
!-- Consistency checks
!
          IF(DEEP .OR. SHALLOW) THEN
            IF(KTOP(1)<1) write(0,*)'WARNING: KTOP,DEEP,SHALLOW=',      &
              KTOP(1),DEEP,SHALLOW
            IF(KBOT(1)>LM) write(0,*)'WARNING: KBOT,DEEP,SHALLOW=',     &
              KBOT(1),DEEP,SHALLOW
          ENDIF
          IF(.NOT.DEEP .AND. .NOT.SHALLOW) THEN
            IF(RN(1)>EPSQ) write(0,*)'WARNING: RAIN,DEEP,SHALLOW=',     &
              RN(1),DEEP,SHALLOW
          ENDIF
!
!*** COMPUTE HEATING, MOISTENING, AND MOMENTUM TENDENCIES
!
convect:  IF (SHALLOW .OR. DEEP) THEN
vloop2:     DO K=1,LM
              KFLIP = LM+1-K
              RTHCUTEN(I,J,KFLIP)=(t1(1,K)-t0(1,K))*RDELT/exner(I,J,KFLIP)
              RQCUTEN(I,J,KFLIP)=(q1(1,K)-q0(1,K))*RDELT
              DUDT(I,J,KFLIP) = mommix*(u1(1,K)-u0(1,K))*RDELT
              DVDT(I,J,KFLIP) = mommix*(v1(1,K)-v0(1,K))*RDELT
!
              TMP = (CLW(1,K,1)-CLW0(1,K,1))*RDELT   !- DETRAINED LIQUID WATER
              RQCCUTEN(I,J,KFLIP) = TMP             
              IF(CLW0(1,K,1)>EPSQ) THEN
                RQCCUTEN(I,J,KFLIP)=TMP*QC(I,J,KFLIP)/CLW0(1,K,1)
                IF(F_QR) RQRCUTEN(I,J,KFLIP)=TMP*QR(I,J,KFLIP)/CLW0(1,K,1)
              ENDIF
dbg2:         if(abs(TMP)>0.1) then
                write(0,*)'WARNING: DETRAINED LIQUID IS TOO LARGE. TMP=',TMP
                write(0,*)'i,j,k,kflip,exner=',i,j,k,kflip,exner(I,J,kflip)
                write(0,*)'t1,t0,rthcuten=',t1(1,k),t0(1,k),rthcuten(i,j,kflip)
                write(0,*)'q1,q0,rqcuten=',q1(1,k),q0(1,k),rqcuten(i,j,kflip)
                write(0,*)'clw,clw0=',clw(1,k,1),clw0(1,k,1)
                write(0,*)'qc,rqccuten=',qc(i,j,kflip),rqccuten(i,j,kflip)
                IF(F_QR) write(0,*)'qr,rqrcuten=',qr(i,j,kflip),rqrcuten(i,j,kflip)
              endif dbg2
!
              TMP = (CLW(1,K,2)-CLW0(1,K,2))/DELT    !- DETRAINED ICE
              IF(F_QI) THEN
                RQICUTEN(I,J,KFLIP) = TMP
              ELSE IF(F_QS) THEN
                RQSCUTEN(I,J,KFLIP) = TMP
              ENDIF
              IF(MULTI_ICE .AND. CLW0(1,K,2)>EPSQ) THEN
                IF(F_QI) RQICUTEN(I,J,KFLIP)=TMP*QI(I,J,KFLIP)/CLW0(1,K,2)
                IF(F_QS) RQSCUTEN(I,J,KFLIP)=TMP*QS(I,J,KFLIP)/CLW0(1,K,2)
                IF(F_QG) RQGCUTEN(I,J,KFLIP)=TMP*QG(I,J,KFLIP)/CLW0(1,K,2)
              ENDIF
dbg3:         if(abs(TMP)>0.1) then
                write(0,*)'WARNING: DETRAINED ICE IS TOO LARGE. TMP=',TMP
                write(0,*)'i,j,k,kflip,exner=',i,j,k,kflip,exner(I,J,kflip)
                write(0,*)'t1,t0,rthcuten=',t1(1,k),t0(1,k),rthcuten(i,j,kflip)
                write(0,*)'q1,q0,rqcuten=',q1(1,k),q0(1,k),rqcuten(i,j,kflip)
                write(0,*)'clw,clw0,delt=',clw(1,k,1),clw0(1,k,1),delt
                IF(F_QI) write(0,*)'qi,rqicuten=',qi(i,j,kflip),rqicuten(i,j,kflip)
                IF(F_QS) write(0,*)'qs,rqscuten=',qs(i,j,kflip),rqscuten(i,j,kflip)
                IF(F_QG) write(0,*)'qg,rqgcuten=',qg(i,j,kflip),rqgcuten(i,j,kflip)
              endif  dbg3
            ENDDO  vloop2
!
!-----------------------------------------------------------------------
!
dbg4:       IF(LPR) THEN
              IF(DEEP) NDEEP=NDEEP+1
              IF(SHALLOW) NSHAL=NSHAL+1
              if(i==idbg .and. j==jdbg) then
                write(0,*)'Conv rain=',rn(1)*1.E3/NCNVC 
                write(0,*)'kbot,ktop,hpbl=,',kbot(1),ktop(1),hpbl(1)
                write(0,*)'shallow,deep=,',shallow,deep
                write(0,*)'islimsk,psp=,',islimsk(1),psp
                write(0,*)'prsl=,',prsl
                write(0,*)'delp=,',delp
                write(0,*)'phil=,',phil
                write(0,*)'exner=',exner(i,j,:)
                write(0,*)'w=,',vvel
                write(0,*)'u0=,',u0
                write(0,*)'u1=,',u1
                write(0,*)'dudt=',dudt(i,j,:)
                write(0,*)'v0=,',v0
                write(0,*)'v1=,',v1
                write(0,*)'dvdt=',dvdt(i,j,:)
                write(0,*)'t0=,',t0
                write(0,*)'t1=,',t1
                write(0,*)'dthdt=',rthcuten(i,j,:)
                write(0,*)'q0=,',q0
                write(0,*)'q1=,',q1
                write(0,*)'dqdt=',rqcuten(i,j,:)
                if(F_QC) THEN
                  write(0,*)'qc=',qc(i,j,:)
                  write(0,*)'dqcdt=',rqccuten(i,j,:)
                endif
                if(F_QR) THEN
                  write(0,*)'qr=',qr(i,j,:)
                  write(0,*)'dqrdt=',rqrcuten(i,j,:)
                endif
                if(F_QI) THEN
                  write(0,*)'qi=',qi(i,j,:)
                  write(0,*)'dqidt=',rqicuten(i,j,:)
                endif
                if(F_QS) THEN
                  write(0,*)'qs=',qs(i,j,:)
                  write(0,*)'dqsdt=',rqscuten(i,j,:)
                endif
                if(F_QG) THEN
                  write(0,*)'qg=',qg(i,j,:)
                  write(0,*)'dqgdt=',rqgcuten(i,j,:)
                endif
                write(0,*)'clw0(1)=',clw0(1,:,1)
                write(0,*)'clw0(2)=',clw0(1,:,2)
              endif           
            ENDIF  dbg4
!
          ENDIF  convect
!
        ENDDO  !- i loop
      ENDDO    !- j loop
!
!.......................................................................
!$omp end parallel do              
!.......................................................................
!
dbg5: IF(LPR) THEN
        write(0,*)'NSHAL,NDEEP=',NSHAL,NDEEP
        write(0,*)'max sfc_rain,location=', maxval(raincv),maxloc(raincv)
        write(0,*)'max dudt,location=', maxval(abs(dudt)),maxloc(abs(dudt))
        write(0,*)'max dvdt,location=', maxval(abs(dvdt)),maxloc(abs(dvdt))
        write(0,*)'max dthdt,location=', maxval(abs(rthcuten)),maxloc(abs(rthcuten))
        write(0,*)'max dqdt,location=', maxval(abs(rqcuten)),maxloc(abs(rqcuten))
        if(F_QC) &
          write(0,*)'max dqcdt,location=', maxval(abs(rqccuten)),maxloc(abs(rqccuten))
        if(F_QR) &
          write(0,*)'max dqrdt,location=', maxval(abs(rqrcuten)),maxloc(abs(rqrcuten))
        if(F_QI) &
          write(0,*)'max dqidt,location=', maxval(abs(rqicuten)),maxloc(abs(rqicuten))
        if(F_QS) &
           write(0,*)'max dqsdt,location=', maxval(abs(rqscuten)),maxloc(abs(rqscuten))
        if(F_QG) &
           write(0,*)'max dqgdt,location=', maxval(abs(rqgcuten)),maxloc(abs(rqgcuten))
      ENDIF  dbg5
!
!
      END SUBROUTINE SASDRV
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        SUBROUTINE SAS_INIT
          CALL GPVS
        END SUBROUTINE SAS_INIT

      END MODULE MODULE_CU_SAS
!
!-----------------------------------------------------------------------
