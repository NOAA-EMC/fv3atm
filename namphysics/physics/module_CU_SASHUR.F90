!-----------------------------------------------------------------------
!
      MODULE MODULE_CU_SASHUR
!
!     HISTORY LOG
!     2014-06-18  Created by Weiguo Wang, move CU_SAS in HWRF to NMMB,
!     TUNED FOR HURRICANE APPLICATIONS
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
      PUBLIC :: SASDRV_HUR
      PUBLIC :: SASHUR_INIT
!
!-----------------------------------------------------------------------
       CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE SASDRV_HUR( &
                        IMS,IME,JMS,JME,LM &
                       ,DT,NTSD,NCNVC &
                       ,TH,T,SICE,VVL,SHEAT,LHEAT,PBLH,U,V &
                       ,Q,QC,QR,QI,QS,QG &
                       ,F_QC,F_QR,F_QI,F_QS,F_QG &
                       ,PHINT,PHMID,exner,RR,DZ &
                       ,XLAND,CU_ACT_FLAG &
                       ,MOMMIX,PGCON,SAS_MASS_FLUX   &   ! hwrf in
                       ,SHALCONV,SHAL_PGCON          &   ! hwrf in
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
      integer(kind=kint),INTENT(IN):: &
       IMS,IME,JMS,JME,LM
!
      integer(kind=kint),INTENT(IN) :: ntsd,NCNVC
      real(kind=kfpt),   INTENT(IN) :: DT
!
!
      real(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(IN):: &
       XLAND,SICE,PBLH,SHEAT,LHEAT
!
      real(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:lm),INTENT(IN):: &
       dz,exner,phmid,rr,t,th,U,V
      real(kind=kfpt),DIMENSION(IMS:IME,1:lm),INTENT(IN):: &
       VVL
!
      real(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:lm+1),INTENT(IN):: &
       phint
!
      real(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:lm),INTENT(IN):: Q,QC,QR,QI,QS,QG
      LOGICAL,INTENT(IN) :: F_QC,F_QR,F_QI,F_QS,F_QG
!
      real(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:lm),optional,intent(inout):: &
       RQCUTEN,RTHCUTEN &
      ,RQCCUTEN,RQRCUTEN &
      ,RQSCUTEN,RQICUTEN &
      ,RQGCUTEN 
! 
      real(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT):: &
       RAINCV
      !real(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: PRATEC
      real(kind=kfpt),DIMENSION(IMS:IME,JMS:JME) :: PRATEC
!
      real(kind=kfpt),DIMENSION(IMS:IME,JMS:JME),INTENT(OUT):: &
       CUBOT,CUTOP
!
      real(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:lm),INTENT(OUT):: &
       DUDT,DVDT 
!
      LOGICAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT):: &
       CU_ACT_FLAG

       real(kind=kfpt), OPTIONAL, INTENT(IN) ::    PGCON,sas_mass_flux  &
                                       ,shal_pgcon,shalconv
     !  integer(kind=kint), OPTIONAL, INTENT(IN) :: shalconv
       real(kind=kfpt), OPTIONAL,   INTENT(IN) ::    MOMMIX
!
      LOGICAL DEEP, SHALLOW
!
!-----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
!-----------------------------------------------------------------------
!      INTEGER :: LBOT,LPBL,LTOP
! 
!      REAL,DIMENSION(1:lm) :: DPCOL,DQDT,DTDT,PCOL,QCOL,TCOL
!
      INTEGER :: I,J,K,ICLDCK,KFLIP

! For SAS
      INTEGER :: KM
      INTEGER, PARAMETER :: IX=1, IM=1, ncloud=1
      INTEGER :: jcap, kcnv(IX), KBOT(IX), KTOP(IX)
      REAL(kind=kind_phys), DIMENSION(IX,lm) :: delp, prsl,phil,q1,t1,  &
                                                u1,v1,VVEL,             &
                                ud_mf,dd_mf,dt_mf, q0,t0,u0,v0
      REAL(kind=kind_phys), DIMENSION(IX) :: psp,cldwrk,rn,slimsk,hpbl, &
                                             hflx,evap,rcs, ps_kpa
      REAL(kind=kind_phys), DIMENSION(IX,lm,2) :: CLW, CLW0  !! 1-ice  2-liquid 
      REAL(kind=kind_phys) :: triggerpert(im)
      REAL(kind=kind_phys) :: fract, tmp, delt, landmask, DTCNVC
      REAL, DIMENSION(lm+1)    :: ZF

      REAL(kind=kind_phys)       :: PGCON_USE,SHAL_PGCON_USE,massf

      LOGICAL :: lpr
      LOGICAL :: MESOSAS
       lpr=.true.
       lpr=.false.
       mesosas=.false.

       !if (.not. lpr) then
       if (lpr) then
         write(0,*)'namelist options used'
         write(0,*)'mommix=',mommix
         write(0,*)'pgcon=',pgcon
         write(0,*)'sas_mass_flux=',sas_mass_flux
         write(0,*)'SHALCONV=',SHALCONV
         write(0,*)'SHAL_PGCON=',SHAL_PGCON
       endif


      DEEP = .TRUE.
      SHALLOW= .FALSE.
      if ( present(shalconv) ) then
      SHALLOW = .TRUE.
      endif
      KM = lm
!  input from arguments
!       mommix = 1.0    !!! HWRF uses this to adjust/tune moment mixing
    
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
            RQRCUTEN(I,J,K) = 0.0
            RQICUTEN(I,J,K) = 0.0
            RQSCUTEN(I,J,K) = 0.0
            RQGCUTEN(I,J,K) = 0.0
         ENDDO
        ENDDO
       ENDDO
!.......................................................................
!$omp end parallel do                
!.......................................................................
      IF ( (.NOT. DEEP) .AND. (.NOT. SHALLOW) ) RETURN

!-----------------------------------------------------------------------
!
!***  PREPARE TO CALL SAS CONVECTION SCHEME
!
  
      if(present(pgcon)) then
         pgcon_use  = pgcon
      else
!        pgcon_use  = 0.7     ! Gregory et al. (1997, QJRMS)
         pgcon_use  = 0.55    ! Zhang & Wu (2003,JAS), used in GFS (25km res spectral)
!        pgcon_use  = 0.2     ! HWRF, for model tuning purposes
!        pgcon_use  = 0.3     ! GFDL, or so I am told

         ! For those attempting to tune pgcon:

         ! The value of 0.55 comes from an observational study of
         ! synoptic-scale deep convection and 0.7 came from an
         ! incorrect fit to the same data.  That value is likely
         ! correct for deep convection at gridscales near that of GFS,
         ! but is questionable in shallow convection, or for scales
         ! much finer than synoptic scales.

         ! Then again, the assumptions of SAS break down when the
         ! gridscale is near the convection scale anyway.  In a large
         ! storm such as a hurricane, there is often no environment to
         ! detrain into since adjancent gridsquares are also undergoing
         ! active convection.  Each gridsquare will no longer have many
         ! updrafts and downdrafts.  At sub-convective timescales, you
         ! will find unstable columns for many (say, 5 second length)
         ! timesteps in a real atmosphere during a convection cell's
         ! lifetime, so forcing it to be neutrally stable is unphysical.

         ! Hence, in scales near the convection scale (cells have
         ! ~0.5-4km diameter in hurricanes), this parameter is more of a
         ! tuning parameter to get a scheme that is inappropriate for
         ! that resolution to do a reasonable job.

         ! Your mileage might vary.

         ! - Sam Trahan
      endif

      if(present(sas_mass_flux)) then
         massf=sas_mass_flux
         ! Use this to reduce the fluxes added by SAS to prevent
         ! computational instability as a result of large fluxes.
      else
         massf=9e9 ! large number to disable check
      endif

      if(present(shal_pgcon)) then
         if(shal_pgcon>=0) then
            shal_pgcon_use  = shal_pgcon
         else
            ! shal_pgcon<0 means use deep pgcon
            shal_pgcon_use  = pgcon_use
         endif
      else
         ! Default: Same as deep convection pgcon
         shal_pgcon_use  = pgcon_use
         ! Read the warning above though.  It may be advisable for
         ! these to be different.  
      endif

!
!-----------------------------------------------------------------------
!
!***  CHECK TO SEE IF THIS IS A CONVECTION TIMESTEP
!                                                                        
      ICLDCK=MOD(ntsd,NCNVC)                                              
!-----------------------------------------------------------------------
!                                                                      
!***  COMPUTE CONVECTION EVERY NCNVC*DT/60.0 MINUTES
!                                                                     

      IF(ICLDCK==0.OR.ntsd==0)THEN                       !!! call convection
!
        DO J=JMS,JME
        DO I=IMS,IME
          CU_ACT_FLAG(I,J)=.TRUE.
        ENDDO
        ENDDO
!
        DTCNVC=DT*NCNVC
!
!.......................................................................
!$omp parallel do                &
!$omp     private(j,i,k,landmask,slimsk,zf,kflip,delt,psp,prsl,delp,phil,u1,&
!$omp             v1,t1,q1,clw,ud_mf,dd_mf,dt_mf,cldwrk,vvel,hflx,evap,hpbl,&
!$omp             kcnv,kbot,ktop,u0,v0,t0,q0,clw0,tmp,fract,rn,jcap)
!.......................................................................
        DO J=JMS,JME  
        DO I=IMS,IME
          triggerpert(1) = 0.0
          RAINCV(I,J)=0.
          PRATEC(I,J)=0.0
!
!***  CONVERT TO BMJ LAND MASK (1.0 FOR SEA; 0.0 FOR LAND)
!
          LANDMASK=XLAND(I,J)-1.
          SLIMSK(1) = 1. - LANDMASK
          IF(SICE(I,J) > 0.5) SLIMSK(1) = 2     !! 0-sea; 1-land; 2-ice 
          RCS(1) = 1.0
!
!***  FILL 1-D VERTICAL ARRAYS 
!
          ZF(1) = 0.0
          DO K=2,LM+1 
           KFLIP = LM + 1 + 1 -K
           ZF(K) = ZF(K-1) + DZ(I,J,KFLIP)
          ENDDO
           delt = 2.0 * DTCNVC
           PSP(1) = PHINT(I,J,lm+1)        ! Surface pressure, Pa
           PS_kpa(1) = 0.001*PSP(1)
          DO K=1,lm
           kflip = LM + 1 -K
           prsl(1,K)  = phmid(I,J,KFLIP)
           delp(1,K)  = RR(I,J,KFLIP)*g99*DZ(I,J,KFLIP) 
           phil(1,K)  = 0.5*(ZF(K) + ZF(K+1) )*g99              
!           u1(1,K)    = (U(I,J  ,KFLIP)+U(I-1,J  ,KFLIP)                       & 
!                        +U(I,J-1,KFLIP)+U(I-1,J-1,KFLIP))*0.25
!           v1(1,K)    = (V(I,J  ,KFLIP)+V(I-1,J  ,KFLIP)                       &
!                        +V(I,J-1,KFLIP)+V(I-1,J-1,KFLIP))*0.25
            u1(1,K) = U(I,J  ,KFLIP)       ! now, input is already at phy point.
            v1(1,K) = V(I,J,KFLIP)

           t1(1,K)    = T(I,J,KFLIP)
           q1(1,K)    = MAX(EPSQ,Q(I,J,KFLIP)) 
           clw(1,K,1) = 0.0
        !   clw(1,K,1) = QC(I,J,KFLIP)+QR(I,J,KFLIP)                 ! Liquid
           if (f_qc) clw(1,K,1) = clw(1,K,1) + QC(I,J,KFLIP)
           if (f_qr) clw(1,K,1) = clw(1,K,1) + QR(I,J,KFLIP)
        !   clw(1,K,2) = QI(I,J,KFLIP)+QS(I,J,KFLIP)+QG(I,J,KFLIP)   ! ICE
           clw(1,K,2) = 0.0
           if (f_qi) clw(1,K,2) = clw(1,K,2) + QI(I,J,KFLIP)
           if (f_qs) clw(1,K,2) = clw(1,K,2) + QS(I,J,KFLIP)
           if (f_qg) clw(1,K,2) = clw(1,K,2) + QG(I,J,KFLIP)
           ud_mf(1,K) = 0.0
           dd_mf(1,K) = 0.0
           dt_mf(1,K) = 0.0 
           cldwrk(1) = 0.0
           VVEL(1,K)    = VVL(i,k)*0.001       !! kpa/s  

          ENDDO
            hflx(1) = SHEAT(I,J)/RR(I,J,LM)/CP            ! W/m2 to K m/s
            evap(1) = LHEAT(I,J)/RR(I,J,LM)/XLV
            hpbl(1) = PBLH(I,J)
             if ( lpr .and. i == 10 .and. j == 10 ) then
              write(0,*)'SHEAT,hflx=',SHEAT(I,J),hflx(1)
              write(0,*)'LHEAT,evap=',LHEAT(I,J),evap(1)
              write(0,*)'hpbl(1)=',hpbl(1)
             endif

           KCNV(1)  = 0     
           KBOT(1)  = KM 
           KTOP(1)  = 1       
           u0 = u1
           v0 = v1
           t0 = t1
           q0 = q1
           clw0 = clw   


!
!-----------------------------------------------------------------------
!***
!***  CALL CONVECTION
!***
      IF(DEEP) THEN                      !! DEEP

     !  CALL sascnvn(im,ix,km,jcap,delt,delp,prsl,psp,phil,clw,           &
     !     q1,t1,u1,v1,cldwrk,rn,kbot,ktop,kcnv,slimsk,                     &
     !     VVEL,ncloud,ud_mf,dd_mf,dt_mf,triggerpert)
!hwrf

! use kpa for delp, prsl, ps
      IF(MESOSAS) THEN
      CALL  sascnvn_meso(im,ix,km,jcap,delt,delp*0.001,prsl*0.001,psp*0.001,phil,clw,   & 
     &     q1,t1,u1,v1,rcs,cldwrk,rn,kbot,ktop,kcnv,slimsk,        &
     &     vvel,ncloud,pgcon_use,massf)    
      ELSE
      CALL  sascnvn_hur(im,ix,km,jcap,delt,delp*0.001,prsl*0.001,psp*0.001,phil,clw,   & 
     &     q1,t1,u1,v1,rcs,cldwrk,rn,kbot,ktop,kcnv,slimsk,        &
     &     vvel,ncloud,pgcon_use,massf)    
      ENDIF                     
!hwrf

!***  CONVECTIVE CLOUD TOP AND BOTTOM FROM THIS CALL
!
       !   CUTOP(I,J) = REAL( lm+1-KTOP(1) )   !BMJ
       !   CUBOT(I,J) = REAL( lm+1-KBOT(1) )   !BMJ
          CUTOP(I,J) = KTOP(1)
          CUBOT(I,J) = KBOT(1)

!***  ALL UNITS IN BMJ SCHEME ARE MKS, THUS CONVERT PRECIP FROM METERS
!***  TO MILLIMETERS PER STEP FOR OUTPUT.
!
          if(lpr .and. i == 20 .and. j == 10)write(0,*)'deep rain=',0.5*rn(1)*1e3/ncnvc 
          
          RAINCV(I,J)=RAINCV(I,J) + 0.5 * rn(1)*1.E3/NCNVC    !! Rain from Deep conv 
          PRATEC(I,J)=PRATEC(I,J) + 0.5 * rn(1)*1.E3/NCNVC/DT
!
      ENDIF                             !! DEEP

      IF (SHALLOW) THEN                 !! Shallow
!      CALL shalcnv(im,ix,km,jcap,delt,delp,prsl,psp,phil,clw,          &
!                   q1,t1,u1,v1,rn,kbot,ktop,kcnv,slimsk,                          &
!                   VVEL,ncloud,hpbl,hflx,evap,ud_mf,dt_mf)

! use kpa for delp, prsl, ps
       CALL shalcnv_hur(im,ix,km,jcap,delt,delp*0.001,prsl*0.001,psp*0.001,phil,clw,   &
     &     q1,t1,u1,v1,rcs,rn,kbot,ktop,kcnv,slimsk,               &
     &     VVEL,ncloud,hpbl,hflx,evap,shal_pgcon_use)


          if(lpr .and.i == 20 .and. j == 10)write(0,*)'shallow rain=',0.5*rn(1)*1.E3/NCNVC 

          RAINCV(I,J)=RAINCV(I,J) + 0.5 * rn(1)*1.E3/NCNVC   !! Rain from shallow conv
          PRATEC(I,J)=PRATEC(I,J) + 0.5 * rn(1)*1.E3/NCNVC/DT

      ENDIF                             !! Shallow

!   compute tendency , either shallow or deep happens. only one of them happens
!***  COMPUTE HEATING AND MOISTENING TENDENCIES
!
              DO K=1,LM
                KFLIP = LM+1-K
         !       DUDT(I,J,KFLIP) = mommix*(u1(1,K)-u0(1,K))/delt
         !       DVDT(I,J,KFLIP) = mommix*(v1(1,K)-v0(1,K))/delt
         ! in HWRF 3.5,  mommix is not actually used.
                DUDT(I,J,KFLIP) = (u1(1,K)-u0(1,K))/delt
                DVDT(I,J,KFLIP) = (v1(1,K)-v0(1,K))/delt
              ENDDO

            IF(PRESENT(RTHCUTEN).AND.PRESENT(RQCUTEN))THEN
              DO K=1,lm
                KFLIP = LM+1-K
                RTHCUTEN(I,J,KFLIP)=(t1(1,K)-t0(1,K))/delt/exner(I,J,KFLIP)
                RQCUTEN(I,J,KFLIP)=(q1(1,K)-q0(1,K))/DELT
              ENDDO
            ENDIF
            IF(    PRESENT(RQCCUTEN).OR.PRESENT(RQRCUTEN)     &
               .OR.PRESENT(RQICUTEN).OR.PRESENT(RQSCUTEN)     &
               .OR.PRESENT(RQGCUTEN))THEN
                 DO K=1,LM                           !! K
                   KFLIP=LM+1-K
                   tmp   = (CLW(1,K,1)-CLW0(1,K,1))/DELT
              ! IF liquid water=0 at t0, then change is assigned to QC tendency
                   RQCCUTEN(I,J,KFLIP) = tmp             
                    IF(CLW0(1,K,1) .GT. EPSQ ) THEN
                       fract = QC(I,J,KFLIP)/CLW0(1,K,1)
                       RQCCUTEN(I,J,KFLIP) = tmp*fract
                       RQRCUTEN(I,J,KFLIP) = tmp*(1.0-fract)
                           
                          if(abs(rqccuten(i,j,kflip)) .gt. 0.1) then
                            write(0,*)'i=,j=',i,j,kflip
                            write(0,*)'qc=',qc(i,j,kflip)
                            write(0,*)'qr=',qr(i,j,kflip)
                            write(0,*)'clw,clw0=',clw(1,k,1),clw0(1,k,1)
                            write(0,*)'rqccuten=',rqccuten(i,j,kflip)
                            write(0,*)'delt=',delt
                            write(0,*)'q1,q0=',q1(1,k),q0(1,k)
                            write(0,*)'t1,t0=',t1(1,k),t0(1,k)
                            stop
                          endif
                    ENDIF

                   tmp   = (CLW(1,K,2)-CLW0(1,K,2))/DELT 
                   RQICUTEN(I,J,KFLIP) = tmp             
                    IF(CLW0(1,K,2) .GT. EPSQ ) THEN

                       RQICUTEN(I,J,KFLIP) = 0.0
                       IF (F_QI) RQICUTEN(I,J,KFLIP) = tmp*QI(I,J,KFLIP)/CLW0(1,K,2)

                       RQSCUTEN(I,J,KFLIP) = 0.0
                       IF (F_QS) RQSCUTEN(I,J,KFLIP) = tmp*QS(I,J,KFLIP)/CLW0(1,K,2)

                       RQGCUTEN(I,J,KFLIP) = 0.0
                       IF (F_QG) RQGCUTEN(I,J,KFLIP) = tmp*QG(I,J,KFLIP)/CLW0(1,K,2)

                    ENDIF
                   
                 ENDDO                              !! K 
            ENDIF 


!
!-----------------------------------------------------------------------
!
        IF(LPR) THEN
          if(i == 20 .and. j == 10) then
           write(0,*)'u1=,',u1
           write(0,*)'v1=,',v1
           write(0,*)'q1=,',q1
           write(0,*)'W=,',vvel
           write(0,*)'psp=,',psp
           write(0,*)'prsl=,',prsl
           write(0,*)'delp=,',delp
           write(0,*)'phil=,',phil
           write(0,*)'dudt=',dudt(i,j,:)
           write(0,*)'dvdt=',dvdt(i,j,:)
           write(0,*)'dthdt=',rthcuten(i,j,:)
           write(0,*)'dqdt=',rqcuten(i,j,:)
           write(0,*)'dqcdt=',rqccuten(i,j,:)
           write(0,*)'dqidt=',rqicuten(i,j,:)
           write(0,*)'dqsdt=',rqscuten(i,j,:)
           write(0,*)'dqgdt=',rqgcuten(i,j,:)
           write(0,*)'max dqdt,location=', maxval(abs(rqcuten)),maxloc(abs(rqcuten))
           write(0,*)'max dqcdt,location=', maxval(abs(rqccuten)),maxloc(abs(rqccuten))
           write(0,*)'max dqrdt,location=', maxval(abs(rqrcuten)),maxloc(abs(rqrcuten))
           write(0,*)'max dqidt,location=', maxval(abs(rqicuten)),maxloc(abs(rqicuten))
           write(0,*)'max dqsdt,location=', maxval(abs(rqscuten)),maxloc(abs(rqscuten))
           write(0,*)'max dqgdt,location=', maxval(abs(rqgcuten)),maxloc(abs(rqgcuten))
           write(0,*)'clw0(1,2)=',clw0(1,10,1),clw0(1,10,2)
           write(0,*)'water(p_qc,p_qr)=',qc(i,j,lm+1-10),qr(i,j,lm+1-10)
           write(0,*)'water(p_qi,p_qs,p_qg)=',qi(i,j,lm+1-10),qs(i,j,lm+1-10),qg(i,j,lm+1-10)
           write(0,*)'EXNER=',exner(i,j,:)
           write(0,*)'hPBL=',hpbl(1)
           write(0,*)'SICE=',sice(i,j)
           write(0,*)'RAIN=,',raincv(I,J)
           write(0,*)'kbot=,',kbot
           write(0,*)'ktop=,',ktop
          endif           

        ENDIF
        ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do              
!.......................................................................
!
      ENDIF                                               !! end of convection
!
      END SUBROUTINE SASDRV_HUR
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        SUBROUTINE SASHUR_INIT
          CALL GPVS
        END SUBROUTINE SASHUR_INIT

! ------------------------------------------------------------------------


! ------------------------------------------------------------------------


!-----------------------------------------------------------------------
      SUBROUTINE TRIDI2T3(L,N,CL,CM,CU,R1,R2,AU,A1,A2)
!yt      INCLUDE DBTRIDI2;
!!
      USE MACHINE , ONLY : kind_phys
      implicit none
      integer             k,n,l,i
      real(kind=kind_phys) fk
!!
      real(kind=kind_phys)                                              &
     &          CL(L,2:N),CM(L,N),CU(L,N-1),R1(L,N),R2(L,N),            &
     &          AU(L,N-1),A1(L,N),A2(L,N)
!-----------------------------------------------------------------------
      DO I=1,L
        FK=1./CM(I,1)
        AU(I,1)=FK*CU(I,1)
        A1(I,1)=FK*R1(I,1)
        A2(I,1)=FK*R2(I,1)
      ENDDO
      DO K=2,N-1
        DO I=1,L
          FK=1./(CM(I,K)-CL(I,K)*AU(I,K-1))
          AU(I,K)=FK*CU(I,K)
          A1(I,K)=FK*(R1(I,K)-CL(I,K)*A1(I,K-1))
          A2(I,K)=FK*(R2(I,K)-CL(I,K)*A2(I,K-1))
        ENDDO
      ENDDO
      DO I=1,L
        FK=1./(CM(I,N)-CL(I,N)*AU(I,N-1))
        A1(I,N)=FK*(R1(I,N)-CL(I,N)*A1(I,N-1))
        A2(I,N)=FK*(R2(I,N)-CL(I,N)*A2(I,N-1))
      ENDDO
      DO K=N-1,1,-1
        DO I=1,L
          A1(I,K)=A1(I,K)-AU(I,K)*A1(I,K+1)
          A2(I,K)=A2(I,K)-AU(I,K)*A2(I,K+1)
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE TRIDI2T3
!-----------------------------------------------------------------------

      SUBROUTINE MSTADBT3(IM,KM,K1,K2,PRSL,PRSLK,TENV,QENV,             &
     &                  KLCL,KBOT,KTOP,TCLD,QCLD)
!yt      INCLUDE DBMSTADB;
!!
      USE MACHINE, ONLY : kind_phys
      USE FUNCPHYS, ONLY : FTDP, FTHE, FTLCL, STMA
      USE PHYSCONS, EPS => con_eps, EPSM1 => con_epsm1, FV => con_FVirt

      implicit none
!!
!     include 'constant.h'
!!
      integer              k,k1,k2,km,i,im
      real(kind=kind_phys) pv,qma,slklcl,tdpd,thelcl,tlcl
      real(kind=kind_phys) tma,tvcld,tvenv
!!
      real(kind=kind_phys) PRSL(IM,KM), PRSLK(IM,KM), TENV(IM,KM),      &
     &                     QENV(IM,KM), TCLD(IM,KM),  QCLD(IM,KM)
      INTEGER              KLCL(IM),    KBOT(IM),      KTOP(IM)
!  LOCAL ARRAYS
      real(kind=kind_phys) SLKMA(IM), THEMA(IM)
!-----------------------------------------------------------------------
!  DETERMINE WARMEST POTENTIAL WET-BULB TEMPERATURE BETWEEN K1 AND K2.
!  COMPUTE ITS LIFTING CONDENSATION LEVEL.
!
      DO I=1,IM
        SLKMA(I) = 0.
        THEMA(I) = 0.
      ENDDO
      DO K=K1,K2
        DO I=1,IM
          PV   = 1000.0 * PRSL(I,K)*QENV(I,K)/(EPS-EPSM1*QENV(I,K))
          TDPD = TENV(I,K)-FTDP(PV)
          IF(TDPD.GT.0.) THEN
            TLCL   = FTLCL(TENV(I,K),TDPD)
            SLKLCL = PRSLK(I,K)*TLCL/TENV(I,K)
          ELSE
            TLCL   = TENV(I,K)
            SLKLCL = PRSLK(I,K)
          ENDIF
          THELCL=FTHE(TLCL,SLKLCL)
          IF(THELCL.GT.THEMA(I)) THEN
            SLKMA(I) = SLKLCL
            THEMA(I) = THELCL
          ENDIF
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!  SET CLOUD TEMPERATURES AND HUMIDITIES WHEREVER THE PARCEL LIFTED UP
!  THE MOIST ADIABAT IS BUOYANT WITH RESPECT TO THE ENVIRONMENT.
      DO I=1,IM
        KLCL(I)=KM+1
        KBOT(I)=KM+1
        KTOP(I)=0
      ENDDO
      DO K=1,KM
        DO I=1,IM
          TCLD(I,K)=0.
          QCLD(I,K)=0.
        ENDDO
      ENDDO
      DO K=K1,KM
        DO I=1,IM
          IF(PRSLK(I,K).LE.SLKMA(I)) THEN
            KLCL(I)=MIN(KLCL(I),K)
            CALL STMA(THEMA(I),PRSLK(I,K),TMA,QMA)
!           TMA=FTMA(THEMA(I),PRSLK(I,K),QMA)
            TVCLD=TMA*(1.+FV*QMA)
            TVENV=TENV(I,K)*(1.+FV*QENV(I,K))
            IF(TVCLD.GT.TVENV) THEN
              KBOT(I)=MIN(KBOT(I),K)
              KTOP(I)=MAX(KTOP(I),K)
              TCLD(I,K)=TMA-TENV(I,K)
              QCLD(I,K)=QMA-QENV(I,K)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE MSTADBT3

      subroutine sascnvn_hur(im,ix,km,jcap,delt,del,prsl,ps,phil,ql,   & 
     &     q1,t1,u1,v1,rcs,cldwrk,rn,kbot,ktop,kcnv,slimsk,        &
     &     dot,ncloud,pgcon,sas_mass_flux)                         
!     &     dot,ncloud,ud_mf,dd_mf,dt_mf)                         
!    &     dot,ncloud,ud_mf,dd_mf,dt_mf,me)
!
!      use machine , only : kind_phys
!      use funcphys , only : fpvs
!      use physcons, grav => con_g, cp => con_cp, hvap => con_hvap  &
      USE MACHINE, ONLY : kind_phys
      USE FUNCPHYS, ONLY : fpvs
      USE PHYSCONS, grav => con_g, cp => con_cp         &
     &,             hvap => con_hvap                               &
     &,             rv => con_rv, fv => con_fvirt, t0c => con_t0c  &
     &,             cvap => con_cvap, cliq => con_cliq             &
     &,             eps => con_eps, epsm1 => con_epsm1
      implicit none
!
      integer            im, ix,  km, jcap, ncloud,                &
     &                   kbot(im), ktop(im), kcnv(im) 
!    &,                  me
      real(kind=kind_phys) delt,sas_mass_flux
      real(kind=kind_phys) ps(im),     del(ix,km),  prsl(ix,km),   &
     &                     ql(ix,km,2),q1(ix,km),   t1(ix,km),     &
     &                     u1(ix,km),  v1(ix,km),   rcs(im),       &
     &                     cldwrk(im), rn(im),      slimsk(im),    &
     &                     dot(ix,km), phil(ix,km)
! hchuang code change mass flux output
!     &,                    ud_mf(im,km),dd_mf(im,km),dt_mf(im,km)
!
      integer              i, j, indx, jmn, k, kk, latd, lond, km1
!
      real(kind=kind_phys) clam, cxlamu, xlamde, xlamdd
! 
      real(kind=kind_phys) adw,     aup,     aafac,                &
     &                     beta,    betal,   betas,                &
     &                     c0,      cpoel,   dellat,  delta,       &
     &                     desdt,   deta,    detad,   dg,          &
     &                     dh,      dhh,     dlnsig,  dp,          &
     &                     dq,      dqsdp,   dqsdt,   dt,          &
     &                     dt2,     dtmax,   dtmin,   dv1h,        &
     &                     dv1q,    dv2h,    dv2q,    dv1u,        &
     &                     dv1v,    dv2u,    dv2v,    dv3q,        &
     &                     dv3h,    dv3u,    dv3v,                 &
     &                     dz,      dz1,     e1,      edtmax,      &
     &                     edtmaxl, edtmaxs, el2orc,  elocp,       &
     &                     es,      etah,    cthk,    dthk,        &
     &                     evef,    evfact,  evfactl, fact1,       &
     &                     fact2,   factor,  fjcap,   fkm,         &
     &                     g,       gamma,   pprime,               &
     &                     qlk,     qrch,    qs,      c1,          &
     &                     rain,    rfact,   shear,   tem1,        &
     &                     tem2,    terr,    val,     val1,        &
     &                     val2,    w1,      w1l,     w1s,         &
     &                     w2,      w2l,     w2s,     w3,          &
     &                     w3l,     w3s,     w4,      w4l,         &
     &                     w4s,     xdby,    xpw,     xpwd,        &
     &                     xqrch,   mbdt,    tem,                  &
     &                     ptem,    ptem1
!
      real(kind=kind_phys), intent(in) :: pgcon

      integer              kb(im), kbcon(im), kbcon1(im),          &
     &                     ktcon(im), ktcon1(im),                  &
     &                     jmin(im), lmin(im), kbmax(im),          &
     &                     kbm(im), kmax(im)
!
      real(kind=kind_phys) aa1(im),     acrt(im),   acrtfct(im),   &
     &                     delhbar(im), delq(im),   delq2(im),     &
     &                     delqbar(im), delqev(im), deltbar(im),   &
     &                     deltv(im),   dtconv(im), edt(im),       &
     &                     edto(im),    edtx(im),   fld(im),       &
     &                     hcdo(im,km), hmax(im),   hmin(im),      &
     &                     ucdo(im,km), vcdo(im,km),aa2(im),       &
     &                     pbcdif(im),  pdot(im),   po(im,km),     &
     &                     pwavo(im),   pwevo(im),  xlamud(im),    &
     &                     qcdo(im,km), qcond(im),  qevap(im),     &
     &                     rntot(im),   vshear(im), xaa0(im),      &
     &                     xk(im),      xlamd(im),                 &
     &                     xmb(im),     xmbmax(im), xpwav(im),     &
     &                     xpwev(im),   delubar(im),delvbar(im)
!cj
      real(kind=kind_phys) cincr, cincrmax, cincrmin
      real(kind=kind_phys) xmbmx1
!cj
!c  physical parameters
      parameter(g=grav)
      parameter(cpoel=cp/hvap,elocp=hvap/cp,                       &
     &          el2orc=hvap*hvap/(rv*cp))
      parameter(terr=0.,c0=.002,c1=.002,delta=fv)
      parameter(fact1=(cvap-cliq)/rv,fact2=hvap/rv-fact1*t0c)
      parameter(cthk=150.,cincrmax=180.,cincrmin=120.,dthk=25.)
!c  local variables and arrays
      real(kind=kind_phys) pfld(im,km),to(im,km), qo(im,km),       &
     &                     uo(im,km),  vo(im,km), qeso(im,km)
!c  cloud water
      real(kind=kind_phys)qlko_ktcon(im),dellal(im,km),tvo(im,km), &
     &                dbyo(im,km), zo(im,km),    xlamue(im,km),    &
     &                fent1(im,km),fent2(im,km), frh(im,km),       &
     &                heo(im,km),  heso(im,km),                    &
     &                qrcd(im,km), dellah(im,km), dellaq(im,km),   &
     &                dellau(im,km),dellav(im,km), hcko(im,km),    &
     &                ucko(im,km), vcko(im,km),   qcko(im,km),     &
     &                eta(im,km),  etad(im,km),   zi(im,km),       &
     &                qrcdo(im,km),pwo(im,km),    pwdo(im,km),     &
     &                tx1(im),     sumx(im)
!    &,               rhbar(im)
!
      logical totflg, cnvflg(im), flg(im)
!
      real(kind=kind_phys) pcrit(15), acritt(15), acrit(15)
!     save pcrit, acritt
      data pcrit/850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,&
     &           350.,300.,250.,200.,150./
      data acritt/.0633,.0445,.0553,.0664,.075,.1082,.1521,.2216,  &
     &           .3151,.3677,.41,.5255,.7663,1.1686,1.6851/
!c  gdas derived acrit
!c     data acritt/.203,.515,.521,.566,.625,.665,.659,.688,
!c    &            .743,.813,.886,.947,1.138,1.377,1.896/
      real(kind=kind_phys) tf, tcr, tcrf
      parameter (tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf))
!
!c-----------------------------------------------------------------------
!

      km1 = km - 1
!c
!c  initialize arrays
!c
      do i=1,im
        kcnv(i)=0
        cnvflg(i) = .true.
        rn(i)=0.
        kbot(i)=km+1
        ktop(i)=0
        kbcon(i)=km
        ktcon(i)=1
        dtconv(i) = 3600.
        cldwrk(i) = 0.
        pdot(i) = 0.
        pbcdif(i)= 0.
        lmin(i) = 1
        jmin(i) = 1
        qlko_ktcon(i) = 0.
        edt(i)  = 0.
        edto(i) = 0.
        edtx(i) = 0.
        acrt(i) = 0.
        acrtfct(i) = 1.
        aa1(i)  = 0.
        aa2(i)  = 0.
        xaa0(i) = 0.
        pwavo(i)= 0.
        pwevo(i)= 0.
        xpwav(i)= 0.
        xpwev(i)= 0.
        vshear(i) = 0.
      enddo
! hchuang code change
!      do k = 1, km
!        do i = 1, im
!          ud_mf(i,k) = 0.
!          dd_mf(i,k) = 0.
!          dt_mf(i,k) = 0.
!        enddo
!      enddo
!c
      do k = 1, 15
        acrit(k) = acritt(k) * (975. - pcrit(k))
      enddo
      dt2 = delt
      val   =         1200.
      dtmin = max(dt2, val )
      val   =         3600.
      dtmax = max(dt2, val )
!c  model tunable parameters are all here
      mbdt    = 10.
      edtmaxl = .3
      edtmaxs = .3
      clam    = .1
      aafac   = .1
!     betal   = .15
!     betas   = .15
      betal   = .05
      betas   = .05
!c     evef    = 0.07
      evfact  = 0.3
      evfactl = 0.3
#if ( EM_CORE == 1 )
!  HAWAII TEST - ZCX
      BETAl   = .05
      betas   = .05
      evfact  = 0.5
      evfactl = 0.5
#endif
!
      cxlamu  = 1.0e-4
      xlamde  = 1.0e-4
      xlamdd  = 1.0e-4
!
      fjcap   = (float(jcap) / 126.) ** 2
      val     =           1.
      fjcap   = max(fjcap,val)
      fkm     = (float(km) / 28.) ** 2
      fkm     = max(fkm,val)
      w1l     = -8.e-3 
      w2l     = -4.e-2
      w3l     = -5.e-3 
      w4l     = -5.e-4
      w1s     = -2.e-4
      w2s     = -2.e-3
      w3s     = -1.e-3
      w4s     = -2.e-5
!c
!c  define top layer for search of the downdraft originating layer
!c  and the maximum thetae for updraft
!c
      do i=1,im
        kbmax(i) = km
        kbm(i)   = km
        kmax(i)  = km
        tx1(i)   = 1.0 / ps(i)
      enddo
!     
      do k = 1, km
        do i=1,im
          IF (prSL(I,K)*tx1(I) .GT. 0.04) KMAX(I)  = MIN(KM,K + 1)
!2011bugfix          if (prsl(i,k)*tx1(i) .gt. 0.04) kmax(i)  = k + 1
          if (prsl(i,k)*tx1(i) .gt. 0.45) kbmax(i) = k + 1
          if (prsl(i,k)*tx1(i) .gt. 0.70) kbm(i)   = k + 1
        enddo
      enddo
      do i=1,im
        kbmax(i) = min(kbmax(i),kmax(i))
        kbm(i)   = min(kbm(i),kmax(i))
      enddo
!c
!c  hydrostatic height assume zero terr and initially assume
!c    updraft entrainment rate as an inverse function of height 
!c
      do k = 1, km
        do i=1,im
          zo(i,k) = phil(i,k) / g
        enddo
      enddo
      do k = 1, km1
        do i=1,im
          zi(i,k) = 0.5*(zo(i,k)+zo(i,k+1))
          xlamue(i,k) = clam / zi(i,k)
        enddo
      enddo
!c
!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c   convert surface pressure to mb from cb
!c
      do k = 1, km
        do i = 1, im
          if (k .le. kmax(i)) then
            pfld(i,k) = prsl(i,k) * 10.0
            eta(i,k)  = 1.
            fent1(i,k)= 1.
            fent2(i,k)= 1.
            frh(i,k)  = 0.
            hcko(i,k) = 0.
            qcko(i,k) = 0.
            ucko(i,k) = 0.
            vcko(i,k) = 0.
            etad(i,k) = 1.
            hcdo(i,k) = 0.
            qcdo(i,k) = 0.
            ucdo(i,k) = 0.
            vcdo(i,k) = 0.
            qrcd(i,k) = 0.
            qrcdo(i,k)= 0.
            dbyo(i,k) = 0.
            pwo(i,k)  = 0.
            pwdo(i,k) = 0.
            dellal(i,k) = 0.
            to(i,k)   = t1(i,k)
            qo(i,k)   = q1(i,k)
            uo(i,k)   = u1(i,k) * rcs(i)
            vo(i,k)   = v1(i,k) * rcs(i)
          endif
        enddo
      enddo
!c
!c  column variables
!c  p is pressure of the layer (mb)
!c  t is temperature at t-dt (k)..tn
!c  q is mixing ratio at t-dt (kg/kg)..qn
!c  to is temperature at t+dt (k)... this is after advection and turbulan
!c  qo is mixing ratio at t+dt (kg/kg)..q1
!c
      do k = 1, km
        do i=1,im
          if (k .le. kmax(i)) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
!           tvo(i,k)  = to(i,k) + delta * to(i,k) * qo(i,k)
          endif
        enddo
      enddo
!c
!c  compute moist static energy
!c
      do k = 1, km
        do i=1,im
          if (k .le. kmax(i)) then
!           tem       = g * zo(i,k) + cp * to(i,k)
            tem       = phil(i,k) + cp * to(i,k)
            heo(i,k)  = tem  + hvap * qo(i,k)
            heso(i,k) = tem  + hvap * qeso(i,k)
!c           heo(i,k)  = min(heo(i,k),heso(i,k))
          endif
        enddo
      enddo
!c
!c  determine level with largest moist static energy
!c  this is the level where updraft starts
!c
      do i=1,im
        hmax(i) = heo(i,1)
        kb(i)   = 1
      enddo
      do k = 2, km
        do i=1,im
          if (k .le. kbm(i)) then
            if(heo(i,k).gt.hmax(i)) then
              kb(i)   = k
              hmax(i) = heo(i,k)
            endif
          endif
        enddo
      enddo
!c
      do k = 1, km1
        do i=1,im
          if (k .le. kmax(i)-1) then
            dz      = .5 * (zo(i,k+1) - zo(i,k))
            dp      = .5 * (pfld(i,k+1) - pfld(i,k))
            es      = 0.01 * fpvs(to(i,k+1))      ! fpvs is in pa
            pprime  = pfld(i,k+1) + epsm1 * es
            qs      = eps * es / pprime
            dqsdp   = - qs / pprime
            desdt   = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
            dqsdt   = qs * pfld(i,k+1) * desdt / (es * pprime)
            gamma   = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
            dt      = (g * dz + hvap * dqsdp * dp) / (cp * (1. + gamma))
            dq      = dqsdt * dt + dqsdp * dp
            to(i,k) = to(i,k+1) + dt
            qo(i,k) = qo(i,k+1) + dq
            po(i,k) = .5 * (pfld(i,k) + pfld(i,k+1))
          endif
        enddo
      enddo
!
      do k = 1, km1
        do i=1,im
          if (k .le. kmax(i)-1) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (po(i,k) + epsm1*qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
            val1      = 1.0
            frh(i,k)  = 1. - min(qo(i,k)/qeso(i,k), val1)
            heo(i,k)  = .5 * g * (zo(i,k) + zo(i,k+1)) +      &
     &                  cp * to(i,k) + hvap * qo(i,k)
            heso(i,k) = .5 * g * (zo(i,k) + zo(i,k+1)) +      &
     &                  cp * to(i,k) + hvap * qeso(i,k)
            uo(i,k)   = .5 * (uo(i,k) + uo(i,k+1))
            vo(i,k)   = .5 * (vo(i,k) + vo(i,k+1))
          endif
        enddo
      enddo
!c
!c  look for the level of free convection as cloud base
!c
      do i=1,im
        flg(i)   = .true.
        kbcon(i) = kmax(i)
      enddo
      do k = 1, km1
        do i=1,im
          if (flg(i).and.k.le.kbmax(i)) then
            if(k.gt.kb(i).and.heo(i,kb(i)).gt.heso(i,k)) then
              kbcon(i) = k
              flg(i)   = .false.
            endif
          endif
        enddo
      enddo
!c
      do i=1,im
        if(kbcon(i).eq.kmax(i)) cnvflg(i) = .false.
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  determine critical convective inhibition
!c  as a function of vertical velocity at cloud base.
!c
      do i=1,im
        if(cnvflg(i)) then
          pdot(i)  = 10.* dot(i,kbcon(i))
        endif
      enddo
      do i=1,im
        if(cnvflg(i)) then
          if(slimsk(i).eq.1.) then
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          endif
          if(pdot(i).le.w4) then
            tem = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i).ge.-w4) then
            tem = - (pdot(i) + w4) / (w4 - w3)
          else
            tem = 0.
          endif
          val1    =             -1.
          tem = max(tem,val1)
          val2    =             1.
          tem = min(tem,val2)
          tem = 1. - tem
          tem1= .5*(cincrmax-cincrmin)
          cincr = cincrmax - tem * tem1
          pbcdif(i) = pfld(i,kb(i)) - pfld(i,kbcon(i))
          if(pbcdif(i).gt.cincr) then
             cnvflg(i) = .false.
          endif
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  assume that updraft entrainment rate above cloud base is
!c    same as that at cloud base
!c
      do k = 2, km1
        do i=1,im
          if(cnvflg(i).and.                            &
     &      (k.gt.kbcon(i).and.k.lt.kmax(i))) then
              xlamue(i,k) = xlamue(i,kbcon(i))
          endif
        enddo
      enddo
!c
!c  assume the detrainment rate for the updrafts to be same as
!c  the entrainment rate at cloud base
!c
      do i = 1, im
        if(cnvflg(i)) then
          xlamud(i) = xlamue(i,kbcon(i))
        endif
      enddo
!c
!c  functions rapidly decreasing with height, mimicking a cloud ensemble
!c    (Bechtold et al., 2008)
!c
      do k = 2, km1
        do i=1,im
          if(cnvflg(i).and.                          &
     &      (k.gt.kbcon(i).and.k.lt.kmax(i))) then
              tem = qeso(i,k)/qeso(i,kbcon(i))
              fent1(i,k) = tem**2
              fent2(i,k) = tem**3
          endif
        enddo
      enddo
!c
!c  final entrainment rate as the sum of turbulent part and organized entrainment
!c    depending on the environmental relative humidity
!c    (Bechtold et al., 2008)
!c
      do k = 2, km1
        do i=1,im
          if(cnvflg(i).and.                         &
     &      (k.ge.kbcon(i).and.k.lt.kmax(i))) then
              tem = cxlamu * frh(i,k) * fent2(i,k)
              xlamue(i,k) = xlamue(i,k)*fent1(i,k) + tem
          endif
        enddo
      enddo
!c
!c  determine updraft mass flux for the subcloud layers
!c
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.lt.kbcon(i).and.k.ge.kb(i)) then
              dz       = zi(i,k+1) - zi(i,k)
              ptem     = 0.5*(xlamue(i,k)+xlamue(i,k+1))-xlamud(i)
              eta(i,k) = eta(i,k+1) / (1. + ptem * dz)
            endif
          endif
        enddo
      enddo
!c
!c  compute mass flux above cloud base
!c
      do k = 2, km1
        do i = 1, im
         if(cnvflg(i))then
           if(k.gt.kbcon(i).and.k.lt.kmax(i)) then
              dz       = zi(i,k) - zi(i,k-1)
              ptem     = 0.5*(xlamue(i,k)+xlamue(i,k-1))-xlamud(i)
              eta(i,k) = eta(i,k-1) * (1 + ptem * dz)
           endif
         endif
        enddo
      enddo
!c
!c  compute updraft cloud properties
!c
      do i = 1, im
        if(cnvflg(i)) then
          indx         = kb(i)
          hcko(i,indx) = heo(i,indx)
          ucko(i,indx) = uo(i,indx)
          vcko(i,indx) = vo(i,indx)
          pwavo(i)     = 0.
        endif
      enddo
!c
!c  cloud property is modified by the entrainment process
!c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.kmax(i)) then
              dz   = zi(i,k) - zi(i,k-1)
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              ptem = 0.5 * tem + pgcon
              ptem1= 0.5 * tem - pgcon
              hcko(i,k) = ((1.-tem1)*hcko(i,k-1)+tem*0.5*     &
     &                     (heo(i,k)+heo(i,k-1)))/factor
              ucko(i,k) = ((1.-tem1)*ucko(i,k-1)+ptem*uo(i,k) &
     &                     +ptem1*uo(i,k-1))/factor
              vcko(i,k) = ((1.-tem1)*vcko(i,k-1)+ptem*vo(i,k) &
     &                     +ptem1*vo(i,k-1))/factor
              dbyo(i,k) = hcko(i,k) - heso(i,k)
            endif
          endif
        enddo
      enddo
!c
!c   taking account into convection inhibition due to existence of
!c    dry layers below cloud base
!c
      do i=1,im
        flg(i) = cnvflg(i)
        kbcon1(i) = kmax(i)
      enddo
      do k = 2, km1
      do i=1,im
        if (flg(i).and.k.lt.kmax(i)) then
          if(k.ge.kbcon(i).and.dbyo(i,k).gt.0.) then
            kbcon1(i) = k
            flg(i)    = .false.
          endif
        endif
      enddo
      enddo
      do i=1,im
        if(cnvflg(i)) then
          if(kbcon1(i).eq.kmax(i)) cnvflg(i) = .false.
        endif
      enddo
      do i=1,im
        if(cnvflg(i)) then
          tem = pfld(i,kbcon(i)) - pfld(i,kbcon1(i))
          if(tem.gt.dthk) then
             cnvflg(i) = .false.
          endif
        endif
      enddo
!!
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  determine first guess cloud top as the level of zero buoyancy
!c
      do i = 1, im
        flg(i) = cnvflg(i)
        ktcon(i) = 1
      enddo
      do k = 2, km1
      do i = 1, im
        if (flg(i).and.k .lt. kmax(i)) then
          if(k.gt.kbcon1(i).and.dbyo(i,k).lt.0.) then
             ktcon(i) = k
             flg(i)   = .false.
          endif
        endif
      enddo
      enddo
!c
      do i = 1, im
        if(cnvflg(i)) then
          tem = pfld(i,kbcon(i))-pfld(i,ktcon(i))
          if(tem.lt.cthk) cnvflg(i) = .false.
        endif
      enddo
!!
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  search for downdraft originating level above theta-e minimum
!c
      do i = 1, im
        if(cnvflg(i)) then
           hmin(i) = heo(i,kbcon1(i))
           lmin(i) = kbmax(i)
           jmin(i) = kbmax(i)
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i) .and. k .le. kbmax(i)) then
            if(k.gt.kbcon1(i).and.heo(i,k).lt.hmin(i)) then
               lmin(i) = k + 1
               hmin(i) = heo(i,k)
            endif
          endif
        enddo
      enddo
!c
!c  make sure that jmin(i) is within the cloud
!c
      do i = 1, im
        if(cnvflg(i)) then
          jmin(i) = min(lmin(i),ktcon(i)-1)
          jmin(i) = max(jmin(i),kbcon1(i)+1)
          if(jmin(i).ge.ktcon(i)) cnvflg(i) = .false.
        endif
      enddo
!c
!c  specify upper limit of mass flux at cloud base
!c
      do i = 1, im
        if(cnvflg(i)) then
!         xmbmax(i) = .1
!
          k = kbcon(i)
          dp = 1000. * del(i,k)
          xmbmax(i) = dp / (g * dt2)
          xmbmax(i) = min(sas_mass_flux,xmbmax(i))
!
!         tem = dp / (g * dt2)
!         xmbmax(i) = min(tem, xmbmax(i))
        endif
      enddo
!c
!c  compute cloud moisture property and precipitation
!c
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = 0.
          qcko(i,kb(i)) = qo(i,kb(i))
!         rhbar(i) = 0.
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k)                             &
     &             + gamma * dbyo(i,k) / (hvap * (1. + gamma))
!cj
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5*  &
     &                     (qo(i,k)+qo(i,k-1)))/factor
!cj
              dq = eta(i,k) * (qcko(i,k) - qrch)
!c
!             rhbar(i) = rhbar(i) + qo(i,k) / qeso(i,k)
!c
!c  check if there is excess moisture to release latent heat
!c
              if(k.ge.kbcon(i).and.dq.gt.0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                if(ncloud.gt.0..and.k.gt.jmin(i)) then
                  dp = 1000. * del(i,k)
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                endif
                aa1(i) = aa1(i) - dz * g * qlk
                qcko(i,k) = qlk + qrch
                pwo(i,k) = etah * c0 * dz * qlk
                pwavo(i) = pwavo(i) + pwo(i,k)
              endif
            endif
          endif
        enddo
      enddo
!c
!     do i = 1, im
!       if(cnvflg(i)) then
!         indx = ktcon(i) - kb(i) - 1
!         rhbar(i) = rhbar(i) / float(indx)
!       endif
!     enddo
!c
!c  calculate cloud work function
!c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.ge.kbcon(i).and.k.lt.ktcon(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma            &
     &                 * to(i,k) / hvap
              aa1(i) = aa1(i) +                           &
     &                 dz1 * (g / (cp * to(i,k)))         &
     &                 * dbyo(i,k) / (1. + gamma)         &
     &                 * rfact
              val = 0.
              aa1(i)=aa1(i)+                              &
     &                 dz1 * g * delta *                  &
     &                 max(val,(qeso(i,k) - qo(i,k)))
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i).and.aa1(i).le.0.) cnvflg(i) = .false.
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  estimate the onvective overshooting as the level 
!c    where the [aafac * cloud work function] becomes zero,
!c    which is the final cloud top
!c
      do i = 1, im
        if (cnvflg(i)) then
          aa2(i) = aafac * aa1(i)
        endif
      enddo
!c
      do i = 1, im
        flg(i) = cnvflg(i)
        ktcon1(i) = kmax(i) - 1
      enddo
      do k = 2, km1
        do i = 1, im
          if (flg(i)) then
            if(k.ge.ktcon(i).and.k.lt.kmax(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma          &
     &                 * to(i,k) / hvap
              aa2(i) = aa2(i) +                         &
     &                 dz1 * (g / (cp * to(i,k)))       &
     &                 * dbyo(i,k) / (1. + gamma)       &
     &                 * rfact
              if(aa2(i).lt.0.) then
                ktcon1(i) = k
                flg(i) = .false.
              endif
            endif
          endif
        enddo
      enddo
!c
!c  compute cloud moisture property, detraining cloud water 
!c    and precipitation in overshooting layers 
!c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.ge.ktcon(i).and.k.lt.ktcon1(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k)                              &
     &             + gamma * dbyo(i,k) / (hvap * (1. + gamma))
!cj
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5*   &
     &                     (qo(i,k)+qo(i,k-1)))/factor
!cj
              dq = eta(i,k) * (qcko(i,k) - qrch)
!c
!c  check if there is excess moisture to release latent heat
!c
              if(dq.gt.0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                if(ncloud.gt.0.) then
                  dp = 1000. * del(i,k)
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                endif
                qcko(i,k) = qlk + qrch
                pwo(i,k) = etah * c0 * dz * qlk
                pwavo(i) = pwavo(i) + pwo(i,k)
              endif
            endif
          endif
        enddo
      enddo
!c
!c exchange ktcon with ktcon1
!c
      do i = 1, im
        if(cnvflg(i)) then
          kk = ktcon(i)
          ktcon(i) = ktcon1(i)
          ktcon1(i) = kk
        endif
      enddo
!c
!c  this section is ready for cloud water
!c
      if(ncloud.gt.0) then
!c
!c  compute liquid and vapor separation at cloud top
!c
      do i = 1, im
        if(cnvflg(i)) then
          k = ktcon(i) - 1
          gamma = el2orc * qeso(i,k) / (to(i,k)**2)
          qrch = qeso(i,k)                              &
     &         + gamma * dbyo(i,k) / (hvap * (1. + gamma))
          dq = qcko(i,k) - qrch
!c
!c  check if there is excess moisture to release latent heat
!c
          if(dq.gt.0.) then
            qlko_ktcon(i) = dq
            qcko(i,k) = qrch
          endif
        endif
      enddo
      endif
!c
!ccccc if(lat.eq.latd.and.lon.eq.lond.and.cnvflg(i)) then
!ccccc   print *, ' aa1(i) before dwndrft =', aa1(i)
!ccccc endif
!c
!c------- downdraft calculations
!c
!c--- compute precipitation efficiency in terms of windshear
!c
      do i = 1, im
        if(cnvflg(i)) then
          vshear(i) = 0.
        endif
      enddo
      do k = 2, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              shear= sqrt((uo(i,k)-uo(i,k-1)) ** 2      &
     &                  + (vo(i,k)-vo(i,k-1)) ** 2)
              vshear(i) = vshear(i) + shear
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          vshear(i) = 1.e3 * vshear(i) / (zi(i,ktcon(i))-zi(i,kb(i)))
          e1=1.591-.639*vshear(i)                       &
     &       +.0953*(vshear(i)**2)-.00496*(vshear(i)**3)
          edt(i)=1.-e1
          val =         .9
          edt(i) = min(edt(i),val)
          val =         .0
          edt(i) = max(edt(i),val)
          edto(i)=edt(i)
          edtx(i)=edt(i)
        endif
      enddo
!c
!c  determine detrainment rate between 1 and kbcon
!c
      do i = 1, im
        if(cnvflg(i)) then
          sumx(i) = 0.
        endif
      enddo
      do k = 1, km1
      do i = 1, im
        if(cnvflg(i).and.k.ge.1.and.k.lt.kbcon(i)) then
          dz = zi(i,k+1) - zi(i,k)
          sumx(i) = sumx(i) + dz
        endif
      enddo
      enddo
      do i = 1, im
        beta = betas
        if(slimsk(i).eq.1.) beta = betal
        if(cnvflg(i)) then
          dz  = (sumx(i)+zi(i,1))/float(kbcon(i))
          tem = 1./float(kbcon(i))
          xlamd(i) = (1.-beta**tem)/dz
        endif
      enddo
!c
!c  determine downdraft mass flux
!c
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)-1) then
           if(k.lt.jmin(i).and.k.ge.kbcon(i)) then
              dz        = zi(i,k+1) - zi(i,k)
              ptem      = xlamdd - xlamde
              etad(i,k) = etad(i,k+1) * (1. - ptem * dz)
           else if(k.lt.kbcon(i)) then
              dz        = zi(i,k+1) - zi(i,k)
              ptem      = xlamd(i) + xlamdd - xlamde
              etad(i,k) = etad(i,k+1) * (1. - ptem * dz)
           endif
          endif
        enddo
      enddo
!c
!c--- downdraft moisture properties
!c
      do i = 1, im
        if(cnvflg(i)) then
          jmn = jmin(i)
          hcdo(i,jmn) = heo(i,jmn)
          qcdo(i,jmn) = qo(i,jmn)
          qrcdo(i,jmn)= qeso(i,jmn)
          ucdo(i,jmn) = uo(i,jmn)
          vcdo(i,jmn) = vo(i,jmn)
          pwevo(i) = 0.
        endif
      enddo
!cj
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i) .and. k.lt.jmin(i)) then
              dz = zi(i,k+1) - zi(i,k)
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = 0.5 * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = 0.5 * (xlamd(i)+xlamdd) * dz
              endif
              factor = 1. + tem - tem1
              ptem = 0.5 * tem - pgcon
              ptem1= 0.5 * tem + pgcon
              hcdo(i,k) = ((1.-tem1)*hcdo(i,k+1)+tem*0.5*       &
     &                     (heo(i,k)+heo(i,k+1)))/factor
              ucdo(i,k) = ((1.-tem1)*ucdo(i,k+1)+ptem*uo(i,k+1) &
     &                     +ptem1*uo(i,k))/factor
              vcdo(i,k) = ((1.-tem1)*vcdo(i,k+1)+ptem*vo(i,k+1) &
     &                     +ptem1*vo(i,k))/factor
              dbyo(i,k) = hcdo(i,k) - heso(i,k)
          endif
        enddo
      enddo
!c
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i).and.k.lt.jmin(i)) then
              gamma      = el2orc * qeso(i,k) / (to(i,k)**2)
              qrcdo(i,k) = qeso(i,k)+                          &
     &                (1./hvap)*(gamma/(1.+gamma))*dbyo(i,k)
!             detad      = etad(i,k+1) - etad(i,k)
!cj
              dz = zi(i,k+1) - zi(i,k)
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = 0.5 * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = 0.5 * (xlamd(i)+xlamdd) * dz
              endif
              factor = 1. + tem - tem1
              qcdo(i,k) = ((1.-tem1)*qcdo(i,k+1)+tem*0.5*     &
     &                     (qo(i,k)+qo(i,k+1)))/factor
!cj
!             pwdo(i,k)  = etad(i,k+1) * qcdo(i,k+1) -
!    &                     etad(i,k) * qrcdo(i,k)
!             pwdo(i,k)  = pwdo(i,k) - detad *
!    &                    .5 * (qrcdo(i,k) + qrcdo(i,k+1))
!cj
              pwdo(i,k)  = etad(i,k+1) * (qcdo(i,k) - qrcdo(i,k))
              qcdo(i,k)  = qrcdo(i,k)
              pwevo(i)   = pwevo(i) + pwdo(i,k)
          endif
        enddo
      enddo
!c
!c--- final downdraft strength dependent on precip
!c--- efficiency (edt), normalized condensate (pwav), and
!c--- evaporate (pwev)
!c
      do i = 1, im
        edtmax = edtmaxl
        if(slimsk(i).eq.0.) edtmax = edtmaxs
        if(cnvflg(i)) then
          if(pwevo(i).lt.0.) then
            edto(i) = -edto(i) * pwavo(i) / pwevo(i)
            edto(i) = min(edto(i),edtmax)
          else
            edto(i) = 0.
          endif
        endif
      enddo
!c
!c--- downdraft cloudwork functions
!c
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i) .and. k .lt. jmin(i)) then
              gamma = el2orc * qeso(i,k) / to(i,k)**2
              dhh=hcdo(i,k)
              dt=to(i,k)
              dg=gamma
              dh=heso(i,k)
              dz=-1.*(zo(i,k+1)-zo(i,k))
              aa1(i)=aa1(i)+edto(i)*dz*(g/(cp*dt))*((dhh-dh)/(1.+dg)) &
     &               *(1.+delta*cp*dg*dt/hvap)
              val=0.
              aa1(i)=aa1(i)+edto(i)*                    &
     &        dz*g*delta*max(val,(qeso(i,k)-qo(i,k)))
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i).and.aa1(i).le.0.) then
           cnvflg(i) = .false.
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c--- what would the change be, that a cloud with unit mass
!c--- will do to the environment?
!c
      do k = 1, km
        do i = 1, im
          if(cnvflg(i) .and. k .le. kmax(i)) then
            dellah(i,k) = 0.
            dellaq(i,k) = 0.
            dellau(i,k) = 0.
            dellav(i,k) = 0.
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          dp = 1000. * del(i,1)
          dellah(i,1) = edto(i) * etad(i,1) * (hcdo(i,1)     &
     &                   - heo(i,1)) * g / dp
          dellaq(i,1) = edto(i) * etad(i,1) * (qcdo(i,1)     &
     &                   - qo(i,1)) * g / dp
          dellau(i,1) = edto(i) * etad(i,1) * (ucdo(i,1)     &
     &                   - uo(i,1)) * g / dp
          dellav(i,1) = edto(i) * etad(i,1) * (vcdo(i,1)     &
     &                   - vo(i,1)) * g / dp
        endif
      enddo
!c
!c--- changed due to subsidence and entrainment
!c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i).and.k.lt.ktcon(i)) then
              aup = 1.
              if(k.le.kb(i)) aup = 0.
              adw = 1.
              if(k.gt.jmin(i)) adw = 0.
              dp = 1000. * del(i,k)
              dz = zi(i,k) - zi(i,k-1)
!c
              dv1h = heo(i,k)
              dv2h = .5 * (heo(i,k) + heo(i,k-1))
              dv3h = heo(i,k-1)
              dv1q = qo(i,k)
              dv2q = .5 * (qo(i,k) + qo(i,k-1))
              dv3q = qo(i,k-1)
              dv1u = uo(i,k)
              dv2u = .5 * (uo(i,k) + uo(i,k-1))
              dv3u = uo(i,k-1)
              dv1v = vo(i,k)
              dv2v = .5 * (vo(i,k) + vo(i,k-1))
              dv3v = vo(i,k-1)
!c
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1))
              tem1 = xlamud(i)
!c
              if(k.le.kbcon(i)) then
                ptem  = xlamde
                ptem1 = xlamd(i)+xlamdd
              else
                ptem  = xlamde
                ptem1 = xlamdd
              endif
!cj
              dellah(i,k) = dellah(i,k) +                           &
     &     ((aup*eta(i,k)-adw*edto(i)*etad(i,k))*dv1h               &
     &    - (aup*eta(i,k-1)-adw*edto(i)*etad(i,k-1))*dv3h           &
     &    - (aup*tem*eta(i,k-1)+adw*edto(i)*ptem*etad(i,k))*dv2h*dz &
     &    +  aup*tem1*eta(i,k-1)*.5*(hcko(i,k)+hcko(i,k-1))*dz      &
     &    +  adw*edto(i)*ptem1*etad(i,k)*.5*(hcdo(i,k)+hcdo(i,k-1)) &
     &         *dz) *g/dp
!cj
              dellaq(i,k) = dellaq(i,k) +                             &
     &     ((aup*eta(i,k)-adw*edto(i)*etad(i,k))*dv1q                 &
     &    - (aup*eta(i,k-1)-adw*edto(i)*etad(i,k-1))*dv3q             &
     &    - (aup*tem*eta(i,k-1)+adw*edto(i)*ptem*etad(i,k))*dv2q*dz   &
     &    +  aup*tem1*eta(i,k-1)*.5*(qcko(i,k)+qcko(i,k-1))*dz        &
     &    +  adw*edto(i)*ptem1*etad(i,k)*.5*(qrcdo(i,k)+qrcdo(i,k-1)) &
     &         *dz) *g/dp
!23456789012345678901234567890123456789012345678901234567890123456789012
!cj
              dellau(i,k) = dellau(i,k) +                             &
     &     ((aup*eta(i,k)-adw*edto(i)*etad(i,k))*dv1u                 &
     &    - (aup*eta(i,k-1)-adw*edto(i)*etad(i,k-1))*dv3u             &
     &    - (aup*tem*eta(i,k-1)+adw*edto(i)*ptem*etad(i,k))*dv2u*dz   &
     &    +  aup*tem1*eta(i,k-1)*.5*(ucko(i,k)+ucko(i,k-1))*dz        &
     &    + adw*edto(i)*ptem1*etad(i,k)*.5*(ucdo(i,k)+ucdo(i,k-1))*dz &
     &    -  pgcon*(aup*eta(i,k-1)-adw*edto(i)*etad(i,k))*(dv1u-dv3u) &
     &         ) *g/dp
!cj
              dellav(i,k) = dellav(i,k) +                             &
     &     ((aup*eta(i,k)-adw*edto(i)*etad(i,k))*dv1v                 &
     &    - (aup*eta(i,k-1)-adw*edto(i)*etad(i,k-1))*dv3v             &
     &    - (aup*tem*eta(i,k-1)+adw*edto(i)*ptem*etad(i,k))*dv2v*dz   &
     &    +  aup*tem1*eta(i,k-1)*.5*(vcko(i,k)+vcko(i,k-1))*dz        &
     &    + adw*edto(i)*ptem1*etad(i,k)*.5*(vcdo(i,k)+vcdo(i,k-1))*dz &
     &    -  pgcon*(aup*eta(i,k-1)-adw*edto(i)*etad(i,k))*(dv1v-dv3v) &
     &         ) *g/dp
!cj
          endif
        enddo
      enddo
!c
!c------- cloud top
!c
      do i = 1, im
        if(cnvflg(i)) then
          indx = ktcon(i)
          dp = 1000. * del(i,indx)
          dv1h = heo(i,indx-1)
          dellah(i,indx) = eta(i,indx-1) *                    &
     &                     (hcko(i,indx-1) - dv1h) * g / dp
          dv1q = qo(i,indx-1)
          dellaq(i,indx) = eta(i,indx-1) *                    &
     &                     (qcko(i,indx-1) - dv1q) * g / dp
          dv1u = uo(i,indx-1)
          dellau(i,indx) = eta(i,indx-1) *                    &
     &                     (ucko(i,indx-1) - dv1u) * g / dp
          dv1v = vo(i,indx-1)
          dellav(i,indx) = eta(i,indx-1) *                    &
     &                     (vcko(i,indx-1) - dv1v) * g / dp
!c
!c  cloud water
!c
          dellal(i,indx) = eta(i,indx-1) *                    &
     &                     qlko_ktcon(i) * g / dp
        endif
      enddo
!c
!c------- final changed variable per unit mass flux
!c
      do k = 1, km
        do i = 1, im
          if (cnvflg(i).and.k .le. kmax(i)) then
            if(k.gt.ktcon(i)) then
              qo(i,k) = q1(i,k)
              to(i,k) = t1(i,k)
            endif
            if(k.le.ktcon(i)) then
              qo(i,k) = dellaq(i,k) * mbdt + q1(i,k)
              dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cp
              to(i,k) = dellat * mbdt + t1(i,k)
              val   =           1.e-10
              qo(i,k) = max(qo(i,k), val  )
            endif
          endif
        enddo
      enddo
!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c
!c--- the above changed environment is now used to calulate the
!c--- effect the arbitrary cloud (with unit mass flux)
!c--- would have on the stability,
!c--- which then is used to calculate the real mass flux,
!c--- necessary to keep this change in balance with the large-scale
!c--- destabilization.
!c
!c--- environmental conditions again, first heights
!c
      do k = 1, km
        do i = 1, im
          if(cnvflg(i) .and. k .le. kmax(i)) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k)+epsm1*qeso(i,k))
            val       =             1.e-8
            qeso(i,k) = max(qeso(i,k), val )
!           tvo(i,k)  = to(i,k) + delta * to(i,k) * qo(i,k)
          endif
        enddo
      enddo
!c
!c--- moist static energy
!c
      do k = 1, km1
        do i = 1, im
          if(cnvflg(i) .and. k .le. kmax(i)-1) then
            dz = .5 * (zo(i,k+1) - zo(i,k))
            dp = .5 * (pfld(i,k+1) - pfld(i,k))
            es = 0.01 * fpvs(to(i,k+1))      ! fpvs is in pa
            pprime = pfld(i,k+1) + epsm1 * es
            qs = eps * es / pprime
            dqsdp = - qs / pprime
            desdt = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
            dqsdt = qs * pfld(i,k+1) * desdt / (es * pprime)
            gamma = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
            dt = (g * dz + hvap * dqsdp * dp) / (cp * (1. + gamma))
            dq = dqsdt * dt + dqsdp * dp
            to(i,k) = to(i,k+1) + dt
            qo(i,k) = qo(i,k+1) + dq
            po(i,k) = .5 * (pfld(i,k) + pfld(i,k+1))
          endif
        enddo
      enddo
      do k = 1, km1
        do i = 1, im
          if(cnvflg(i) .and. k .le. kmax(i)-1) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (po(i,k) + epsm1 * qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
            heo(i,k)   = .5 * g * (zo(i,k) + zo(i,k+1)) +     &
     &                    cp * to(i,k) + hvap * qo(i,k)
            heso(i,k) = .5 * g * (zo(i,k) + zo(i,k+1)) +      &
     &                  cp * to(i,k) + hvap * qeso(i,k)
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          k = kmax(i)
          heo(i,k) = g * zo(i,k) + cp * to(i,k) + hvap * qo(i,k)
          heso(i,k) = g * zo(i,k) + cp * to(i,k) + hvap * qeso(i,k)
!c         heo(i,k) = min(heo(i,k),heso(i,k))
        endif
      enddo
!c
!c**************************** static control
!c
!c------- moisture and cloud work functions
!c
      do i = 1, im
        if(cnvflg(i)) then
          xaa0(i) = 0.
          xpwav(i) = 0.
        endif
      enddo
!c
      do i = 1, im
        if(cnvflg(i)) then
          indx = kb(i)
          hcko(i,indx) = heo(i,indx)
          qcko(i,indx) = qo(i,indx)
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              dz = zi(i,k) - zi(i,k-1)
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              hcko(i,k) = ((1.-tem1)*hcko(i,k-1)+tem*0.5*    &
     &                     (heo(i,k)+heo(i,k-1)))/factor
            endif
          endif
        enddo
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              dz = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              xdby = hcko(i,k) - heso(i,k)
              xqrch = qeso(i,k)                             &
     &              + gamma * xdby / (hvap * (1. + gamma))
!cj
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5*   &
     &                     (qo(i,k)+qo(i,k-1)))/factor
!cj
              dq = eta(i,k) * (qcko(i,k) - xqrch)
!c
              if(k.ge.kbcon(i).and.dq.gt.0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                if(ncloud.gt.0..and.k.gt.jmin(i)) then
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                endif
                if(k.lt.ktcon1(i)) then
                  xaa0(i) = xaa0(i) - dz * g * qlk
                endif
                qcko(i,k) = qlk + xqrch
                xpw = etah * c0 * dz * qlk
                xpwav(i) = xpwav(i) + xpw
              endif
            endif
            if(k.ge.kbcon(i).and.k.lt.ktcon1(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma          &
     &                 * to(i,k) / hvap
              xaa0(i) = xaa0(i)                         &
     &                + dz1 * (g / (cp * to(i,k)))      &
     &                * xdby / (1. + gamma)             &
     &                * rfact
              val=0.
              xaa0(i)=xaa0(i)+                          &
     &                 dz1 * g * delta *                &
     &                 max(val,(qeso(i,k) - qo(i,k)))
            endif
          endif
        enddo
      enddo
!c
!c------- downdraft calculations
!c
!c--- downdraft moisture properties
!c
      do i = 1, im
        if(cnvflg(i)) then
          jmn = jmin(i)
          hcdo(i,jmn) = heo(i,jmn)
          qcdo(i,jmn) = qo(i,jmn)
          qrcd(i,jmn) = qeso(i,jmn)
          xpwev(i) = 0.
        endif
      enddo
!cj
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i) .and. k.lt.jmin(i)) then
              dz = zi(i,k+1) - zi(i,k)
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = 0.5 * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = 0.5 * (xlamd(i)+xlamdd) * dz
              endif
              factor = 1. + tem - tem1
              hcdo(i,k) = ((1.-tem1)*hcdo(i,k+1)+tem*0.5*      &
     &                     (heo(i,k)+heo(i,k+1)))/factor
          endif
        enddo
      enddo
!cj
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i) .and. k .lt. jmin(i)) then
              dq = qeso(i,k)
              dt = to(i,k)
              gamma    = el2orc * dq / dt**2
              dh       = hcdo(i,k) - heso(i,k)
              qrcd(i,k)=dq+(1./hvap)*(gamma/(1.+gamma))*dh
!             detad    = etad(i,k+1) - etad(i,k)
!cj
              dz = zi(i,k+1) - zi(i,k)
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = 0.5 * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = 0.5 * (xlamd(i)+xlamdd) * dz
              endif
              factor = 1. + tem - tem1
              qcdo(i,k) = ((1.-tem1)*qcdo(i,k+1)+tem*0.5*     &
     &                     (qo(i,k)+qo(i,k+1)))/factor
!cj
!             xpwd     = etad(i,k+1) * qcdo(i,k+1) -
!    &                   etad(i,k) * qrcd(i,k)
!             xpwd     = xpwd - detad *
!    &                 .5 * (qrcd(i,k) + qrcd(i,k+1))
!cj
              xpwd     = etad(i,k+1) * (qcdo(i,k) - qrcd(i,k))
              qcdo(i,k)= qrcd(i,k)
              xpwev(i) = xpwev(i) + xpwd
          endif
        enddo
      enddo
!c
      do i = 1, im
        edtmax = edtmaxl
        if(slimsk(i).eq.0.) edtmax = edtmaxs
        if(cnvflg(i)) then
          if(xpwev(i).ge.0.) then
            edtx(i) = 0.
          else
            edtx(i) = -edtx(i) * xpwav(i) / xpwev(i)
            edtx(i) = min(edtx(i),edtmax)
          endif
        endif
      enddo
!c
!c
!c--- downdraft cloudwork functions
!c
!c
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i) .and. k.lt.jmin(i)) then
              gamma = el2orc * qeso(i,k) / to(i,k)**2
              dhh=hcdo(i,k)
              dt= to(i,k)
              dg= gamma
              dh= heso(i,k)
              dz=-1.*(zo(i,k+1)-zo(i,k))
              xaa0(i)=xaa0(i)+edtx(i)*dz*(g/(cp*dt))*((dhh-dh)/(1.+dg)) &
     &                *(1.+delta*cp*dg*dt/hvap)
              val=0.
              xaa0(i)=xaa0(i)+edtx(i)*                         &
     &        dz*g*delta*max(val,(qeso(i,k)-qo(i,k)))
          endif
        enddo
      enddo
!c
!c  calculate critical cloud work function
!c
      do i = 1, im
        if(cnvflg(i)) then
          if(pfld(i,ktcon(i)).lt.pcrit(15))then
            acrt(i)=acrit(15)*(975.-pfld(i,ktcon(i)))          &
     &              /(975.-pcrit(15))
          else if(pfld(i,ktcon(i)).gt.pcrit(1))then
            acrt(i)=acrit(1)
          else
            k =  int((850. - pfld(i,ktcon(i)))/50.) + 2
            k = min(k,15)
            k = max(k,2)
            acrt(i)=acrit(k)+(acrit(k-1)-acrit(k))*            &
     &           (pfld(i,ktcon(i))-pcrit(k))/(pcrit(k-1)-pcrit(k))
          endif
        endif
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          if(slimsk(i).eq.1.) then
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          endif
!c
!c  modify critical cloud workfunction by cloud base vertical velocity
!c
          if(pdot(i).le.w4) then
            acrtfct(i) = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i).ge.-w4) then
            acrtfct(i) = - (pdot(i) + w4) / (w4 - w3)
          else
            acrtfct(i) = 0.
          endif
          val1    =             -1.
          acrtfct(i) = max(acrtfct(i),val1)
          val2    =             1.
          acrtfct(i) = min(acrtfct(i),val2)
          acrtfct(i) = 1. - acrtfct(i)
!c
!c  modify acrtfct(i) by colume mean rh if rhbar(i) is greater than 80 percent
!c
!c         if(rhbar(i).ge..8) then
!c           acrtfct(i) = acrtfct(i) * (.9 - min(rhbar(i),.9)) * 10.
!c         endif
!c
!c  modify adjustment time scale by cloud base vertical velocity
!c
          val1=0.
          dtconv(i) = dt2 + max((1800. - dt2),val1) *         &
     &                (pdot(i) - w2) / (w1 - w2)
!c         dtconv(i) = max(dtconv(i), dt2)
!c         dtconv(i) = 1800. * (pdot(i) - w2) / (w1 - w2)
          dtconv(i) = max(dtconv(i),dtmin)
          dtconv(i) = min(dtconv(i),dtmax)
!c
        endif
      enddo
!c
!c--- large scale forcing
!c
      xmbmx1=-1.e20
      do i= 1, im
        if(cnvflg(i)) then
          fld(i)=(aa1(i)-acrt(i)* acrtfct(i))/dtconv(i)
          if(fld(i).le.0.) cnvflg(i) = .false.
        endif
        if(cnvflg(i)) then
!c         xaa0(i) = max(xaa0(i),0.)
          xk(i) = (xaa0(i) - aa1(i)) / mbdt
          if(xk(i).ge.0.) cnvflg(i) = .false.
        endif
!c
!c--- kernel, cloud base mass flux
!c
        if(cnvflg(i)) then
          xmb(i) = -fld(i) / xk(i)
          xmb(i) = min(xmb(i),xmbmax(i))
          xmbmx1=max(xmbmx1,xmb(i))
        endif
      enddo
!      if(xmbmx1.gt.0.4)print*,'qingfu test xmbmx1=',xmbmx1
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  restore to,qo,uo,vo to t1,q1,u1,v1 in case convection stops
!c
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)) then
            to(i,k) = t1(i,k)
            qo(i,k) = q1(i,k)
            uo(i,k) = u1(i,k)
            vo(i,k) = v1(i,k)
            qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val     =             1.e-8
            qeso(i,k) = max(qeso(i,k), val )
          endif
        enddo
      enddo
!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c
!c--- feedback: simply the changes from the cloud with unit mass flux
!c---           multiplied by  the mass flux necessary to keep the
!c---           equilibrium with the larger-scale.
!c
      do i = 1, im
        delhbar(i) = 0.
        delqbar(i) = 0.
        deltbar(i) = 0.
        delubar(i) = 0.
        delvbar(i) = 0.
        qcond(i) = 0.
      enddo
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)) then
            if(k.le.ktcon(i)) then
              dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cp
              t1(i,k) = t1(i,k) + dellat * xmb(i) * dt2
              q1(i,k) = q1(i,k) + dellaq(i,k) * xmb(i) * dt2
              tem = 1./rcs(i)
              u1(i,k) = u1(i,k) + dellau(i,k) * xmb(i) * dt2 * tem
              v1(i,k) = v1(i,k) + dellav(i,k) * xmb(i) * dt2 * tem
              dp = 1000. * del(i,k)
              delhbar(i) = delhbar(i) + dellah(i,k)*xmb(i)*dp/g
              delqbar(i) = delqbar(i) + dellaq(i,k)*xmb(i)*dp/g
              deltbar(i) = deltbar(i) + dellat*xmb(i)*dp/g
              delubar(i) = delubar(i) + dellau(i,k)*xmb(i)*dp/g
              delvbar(i) = delvbar(i) + dellav(i,k)*xmb(i)*dp/g
            endif
          endif
        enddo
      enddo
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)) then
            if(k.le.ktcon(i)) then
              qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
              qeso(i,k) = eps * qeso(i,k)/(pfld(i,k) + epsm1*qeso(i,k))
              val     =             1.e-8
              qeso(i,k) = max(qeso(i,k), val )
            endif
          endif
        enddo
      enddo
!c
      do i = 1, im
        rntot(i) = 0.
        delqev(i) = 0.
        delq2(i) = 0.
        flg(i) = cnvflg(i)
      enddo
      do k = km, 1, -1
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)) then
            if(k.lt.ktcon(i)) then
              aup = 1.
              if(k.le.kb(i)) aup = 0.
              adw = 1.
              if(k.ge.jmin(i)) adw = 0.
              rain =  aup * pwo(i,k) + adw * edto(i) * pwdo(i,k)
              rntot(i) = rntot(i) + rain * xmb(i) * .001 * dt2
            endif
          endif
        enddo
      enddo
      do k = km, 1, -1
        do i = 1, im
          if (k .le. kmax(i)) then
            deltv(i) = 0.
            delq(i) = 0.
            qevap(i) = 0.
            if(cnvflg(i).and.k.lt.ktcon(i)) then
              aup = 1.
              if(k.le.kb(i)) aup = 0.
              adw = 1.
              if(k.ge.jmin(i)) adw = 0.
              rain =  aup * pwo(i,k) + adw * edto(i) * pwdo(i,k)
              rn(i) = rn(i) + rain * xmb(i) * .001 * dt2
            endif
            if(flg(i).and.k.lt.ktcon(i)) then
              evef = edt(i) * evfact
              if(slimsk(i).eq.1.) evef=edt(i) * evfactl
!             if(slimsk(i).eq.1.) evef=.07
!c             if(slimsk(i).ne.1.) evef = 0.
              qcond(i) = evef * (q1(i,k) - qeso(i,k))     &
     &                 / (1. + el2orc * qeso(i,k) / t1(i,k)**2)
              dp = 1000. * del(i,k)
              if(rn(i).gt.0..and.qcond(i).lt.0.) then
                qevap(i) = -qcond(i) * (1.-exp(-.32*sqrt(dt2*rn(i))))
                qevap(i) = min(qevap(i), rn(i)*1000.*g/dp)
                delq2(i) = delqev(i) + .001 * qevap(i) * dp / g
              endif
              if(rn(i).gt.0..and.qcond(i).lt.0..and.      &
     &           delq2(i).gt.rntot(i)) then
                qevap(i) = 1000.* g * (rntot(i) - delqev(i)) / dp
                flg(i) = .false.
              endif
              if(rn(i).gt.0..and.qevap(i).gt.0.) then
                q1(i,k) = q1(i,k) + qevap(i)
                t1(i,k) = t1(i,k) - elocp * qevap(i)
                rn(i) = rn(i) - .001 * qevap(i) * dp / g
                deltv(i) = - elocp*qevap(i)/dt2
                delq(i) =  + qevap(i)/dt2
                delqev(i) = delqev(i) + .001*dp*qevap(i)/g
              endif
              dellaq(i,k) = dellaq(i,k) + delq(i) / xmb(i)
              delqbar(i) = delqbar(i) + delq(i)*dp/g
              deltbar(i) = deltbar(i) + deltv(i)*dp/g
            endif
          endif
        enddo
      enddo
!cj
!     do i = 1, im
!     if(me.eq.31.and.cnvflg(i)) then
!     if(cnvflg(i)) then
!       print *, ' deep delhbar, delqbar, deltbar = ',
!    &             delhbar(i),hvap*delqbar(i),cp*deltbar(i)
!       print *, ' deep delubar, delvbar = ',delubar(i),delvbar(i)
!       print *, ' precip =', hvap*rn(i)*1000./dt2
!       print*,'pdif= ',pfld(i,kbcon(i))-pfld(i,ktcon(i))
!     endif
!     enddo
!c
!c  precipitation rate converted to actual precip
!c  in unit of m instead of kg
!c
      do i = 1, im
        if(cnvflg(i)) then
!c
!c  in the event of upper level rain evaporation and lower level downdraft
!c    moistening, rn can become negative, in this case, we back out of the
!c    heating and the moistening
!c
          if(rn(i).lt.0..and..not.flg(i)) rn(i) = 0.
          if(rn(i).le.0.) then
            rn(i) = 0.
          else
            ktop(i) = ktcon(i)
            kbot(i) = kbcon(i)
            kcnv(i) = 1
            cldwrk(i) = aa1(i)
          endif
        endif
      enddo
!c
!c  cloud water
!c
      if (ncloud.gt.0) then
!
      val1=1.0
      val2=0.0
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. rn(i).gt.0.) then
            if (k.gt.kb(i).and.k.le.ktcon(i)) then
              tem  = dellal(i,k) * xmb(i) * dt2
              tem1 = max(val2, min(val1, (tcr-t1(i,k))*tcrf))
              if (ql(i,k,2) .gt. -999.0) then
                ql(i,k,1) = ql(i,k,1) + tem * tem1            ! ice
                ql(i,k,2) = ql(i,k,2) + tem *(1.0-tem1)       ! water
              else
                ql(i,k,1) = ql(i,k,1) + tem
              endif
            endif
          endif
        enddo
      enddo
!
      endif
!c
      do k = 1, km
        do i = 1, im
          if(cnvflg(i).and.rn(i).le.0.) then
            if (k .le. kmax(i)) then
              t1(i,k) = to(i,k)
              q1(i,k) = qo(i,k)
              u1(i,k) = uo(i,k)
              v1(i,k) = vo(i,k)
            endif
          endif
        enddo
      enddo
!
! hchuang code change
!
!      do k = 1, km
!        do i = 1, im
!          if(cnvflg(i).and.rn(i).gt.0.) then
!            if(k.ge.kb(i) .and. k.lt.ktop(i)) then
!              ud_mf(i,k) = eta(i,k) * xmb(i) * dt2
!            endif
!          endif
!        enddo
!      enddo
!      do i = 1, im
!        if(cnvflg(i).and.rn(i).gt.0.) then
!           k = ktop(i)-1
!           dt_mf(i,k) = ud_mf(i,k)
!        endif
!      enddo
!      do k = 1, km
!        do i = 1, im
!          if(cnvflg(i).and.rn(i).gt.0.) then
!            if(k.ge.1 .and. k.le.jmin(i)) then
!              dd_mf(i,k) = edto(i) * etad(i,k) * xmb(i) * dt2
!            endif
!          endif
!        enddo
!      enddo
!!
      return
      end subroutine sascnvn_hur      

      subroutine shalcnv_hur(im,ix,km,jcap,delt,del,prsl,ps,phil,ql,   &
     &     q1,t1,u1,v1,rcs,rn,kbot,ktop,kcnv,slimsk,               &
     &     dot,ncloud,hpbl,heat,evap,pgcon)
!
      use machine , only : kind_phys
      use funcphys , only : fpvs
      use physcons, grav => con_g, cp => con_cp, hvap => con_hvap         &
!      use MODULE_GFS_machine , only : kind_phys
!      use MODULE_GFS_funcphys , only : fpvs
!      use MODULE_GFS_physcons, grav => con_g, cp => con_cp, hvap => con_hvap         &
     &,             rv => con_rv, fv => con_fvirt, t0c => con_t0c         &
     &,             rd => con_rd, cvap => con_cvap, cliq => con_cliq      &
     &,             eps => con_eps, epsm1 => con_epsm1
      implicit none
!
      integer            im, ix,  km, jcap, ncloud,                       &
     &                   kbot(im), ktop(im), kcnv(im)                   
      real(kind=kind_phys) delt
      real(kind=kind_phys) ps(im),     del(ix,km),  prsl(ix,km),          &
     &                     ql(ix,km,2),q1(ix,km),   t1(ix,km),            &
     &                     u1(ix,km),  v1(ix,km),   rcs(im),              &
     &                     rn(im),     slimsk(im),                        &
     &                     dot(ix,km), phil(ix,km), hpbl(im),             &
     &                     heat(im),   evap(im)                           
!     &,                    ud_mf(im,km),dt_mf(im,km)

      real  ud_mf(im,km),dt_mf(im,km)
!
      integer              i,j,indx, jmn, k, kk, latd, lond, km1
      integer              kpbl(im)
!
      real(kind=kind_phys) c0,      cpoel,   dellat,  delta,        &
     &                     desdt,   deta,    detad,   dg,           &
     &                     dh,      dhh,     dlnsig,  dp,           &
     &                     dq,      dqsdp,   dqsdt,   dt,           &
     &                     dt2,     dtmax,   dtmin,   dv1h,         &
     &                     dv1q,    dv2h,    dv2q,    dv1u,         &
     &                     dv1v,    dv2u,    dv2v,    dv3q,         &
     &                     dv3h,    dv3u,    dv3v,    clam,         &
     &                     dz,      dz1,     e1,                    &
     &                     el2orc,  elocp,   aafac,   cthk,         &
     &                     es,      etah,    h1,      dthk,         &
     &                     evef,    evfact,  evfactl, fact1,        &
     &                     fact2,   factor,  fjcap,                 &
     &                     g,       gamma,   pprime,  betaw,        &
     &                     qlk,     qrch,    qs,      c1,           &
     &                     rain,    rfact,   shear,   tem1,         &
     &                     tem2,    terr,    val,     val1,         &
     &                     val2,    w1,      w1l,     w1s,          &
     &                     w2,      w2l,     w2s,     w3,           &
     &                     w3l,     w3s,     w4,      w4l,          &
     &                     w4s,     tem,     ptem,    ptem1,        &
     &                     pgcon
!
      integer              kb(im), kbcon(im), kbcon1(im),           &
     &                     ktcon(im), ktcon1(im),                   &
     &                     kbm(im), kmax(im)
!
      real(kind=kind_phys) aa1(im),                                 &
     &                     delhbar(im), delq(im),   delq2(im),      &
     &                     delqbar(im), delqev(im), deltbar(im),    &
     &                     deltv(im),   edt(im),                    &
     &                     wstar(im),   sflx(im),                   &
     &                     pdot(im),    po(im,km),                  &
     &                     qcond(im),   qevap(im),  hmax(im),       &
     &                     rntot(im),   vshear(im),                 &
     &                     xlamud(im),  xmb(im),    xmbmax(im),     &
     &                     delubar(im), delvbar(im)
!c
      real(kind=kind_phys) cincr, cincrmax, cincrmin
!cc
!c  physical parameters
      parameter(g=grav)
      parameter(cpoel=cp/hvap,elocp=hvap/cp,                            &
     &          el2orc=hvap*hvap/(rv*cp))
      parameter(terr=0.,c0=.002,c1=5.e-4,delta=fv)
      parameter(fact1=(cvap-cliq)/rv,fact2=hvap/rv-fact1*t0c)
      parameter(cthk=50.,cincrmax=180.,cincrmin=120.,dthk=25.)
      parameter(h1=0.33333333)
!c  local variables and arrays
      real(kind=kind_phys) pfld(im,km),    to(im,km),     qo(im,km),    &
     &                     uo(im,km),      vo(im,km),     qeso(im,km)
!c  cloud water
      real(kind=kind_phys) qlko_ktcon(im), dellal(im,km),                   &
     &                     dbyo(im,km),    zo(im,km),     xlamue(im,km),    &
     &                     heo(im,km),     heso(im,km),                     &
     &                     dellah(im,km),  dellaq(im,km),                   &
     &                     dellau(im,km),  dellav(im,km), hcko(im,km),      &
     &                     ucko(im,km),    vcko(im,km),   qcko(im,km),      &
     &                     eta(im,km),     zi(im,km),     pwo(im,km),       &
     &                     tx1(im)
!
      logical totflg, cnvflg(im), flg(im)
!
      real(kind=kind_phys) tf, tcr, tcrf
      parameter (tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf))
!
!c-----------------------------------------------------------------------
!
      km1 = km - 1
!c
!c  compute surface buoyancy flux
!c
      do i=1,im
        sflx(i) = heat(i)+fv*t1(i,1)*evap(i)
      enddo
!c
!c  initialize arrays
!c
      do i=1,im
        cnvflg(i) = .true.
        if(kcnv(i).eq.1) cnvflg(i) = .false.
        if(sflx(i).le.0.) cnvflg(i) = .false.
        if(cnvflg(i)) then
          kbot(i)=km+1
          ktop(i)=0
        endif
        rn(i)=0.
        kbcon(i)=km
        ktcon(i)=1
        kb(i)=km
        pdot(i) = 0.
        qlko_ktcon(i) = 0.
        edt(i)  = 0.
        aa1(i)  = 0.
        vshear(i) = 0.
      enddo
! hchuang code change
      do k = 1, km
        do i = 1, im
          ud_mf(i,k) = 0.
          dt_mf(i,k) = 0.
        enddo
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
      dt2   = delt
      val   =         1200.
      dtmin = max(dt2, val )
      val   =         3600.
      dtmax = max(dt2, val )
!c  model tunable parameters are all here
      clam    = .3
      aafac   = .1
      betaw   = .03
!c     evef    = 0.07
      evfact  = 0.3
      evfactl = 0.3
!
      fjcap   = (float(jcap) / 126.) ** 2
      val     =           1.
      fjcap   = max(fjcap,val)
      w1l     = -8.e-3 
      w2l     = -4.e-2
      w3l     = -5.e-3 
      w4l     = -5.e-4
      w1s     = -2.e-4
      w2s     = -2.e-3
      w3s     = -1.e-3
      w4s     = -2.e-5
!c
!c  define top layer for search of the downdraft originating layer
!c  and the maximum thetae for updraft
!c
      do i=1,im
        kbm(i)   = km
        kmax(i)  = km
        tx1(i)   = 1.0 / ps(i)
      enddo
!     
      do k = 1, km
        do i=1,im
          if (prsl(i,k)*tx1(i) .gt. 0.70) kbm(i)   = k + 1
          if (prsl(i,k)*tx1(i) .gt. 0.60) kmax(i)  = k + 1
        enddo
      enddo
      do i=1,im
        kbm(i)   = min(kbm(i),kmax(i))
      enddo
!c
!!c  hydrostatic height assume zero terr and compute
!c  updraft entrainment rate as an inverse function of height
!c
      do k = 1, km
        do i=1,im
          zo(i,k) = phil(i,k) / g
        enddo
      enddo
      do k = 1, km1
        do i=1,im
          zi(i,k) = 0.5*(zo(i,k)+zo(i,k+1))
          xlamue(i,k) = clam / zi(i,k)
        enddo
      enddo
      do i=1,im
        xlamue(i,km) = xlamue(i,km1)
      enddo
!c
!c  pbl height
!c
      do i=1,im
        flg(i) = cnvflg(i)
        kpbl(i)= 1
      enddo
      do k = 2, km1
        do i=1,im
          if (flg(i).and.zo(i,k).le.hpbl(i)) then
            kpbl(i) = k
          else
            flg(i) = .false.
          endif
        enddo
      enddo
      do i=1,im
        kpbl(i)= min(kpbl(i),kbm(i))
      enddo
!c
!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c   convert surface pressure to mb from cb
!c
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)) then
            pfld(i,k) = prsl(i,k) * 10.0
            eta(i,k)  = 1.
            hcko(i,k) = 0.
            qcko(i,k) = 0.
            ucko(i,k) = 0.
            vcko(i,k) = 0.
            dbyo(i,k) = 0.
            pwo(i,k)  = 0.
            dellal(i,k) = 0.
            to(i,k)   = t1(i,k)
            qo(i,k)   = q1(i,k)
            uo(i,k)   = u1(i,k) * rcs(i)
            vo(i,k)   = v1(i,k) * rcs(i)
          endif
        enddo
      enddo
!c
!c  column variables
!c  p is pressure of the layer (mb)
!c  t is temperature at t-dt (k)..tn
!c  q is mixing ratio at t-dt (kg/kg)..qn
!c  to is temperature at t+dt (k)... this is after advection and turbulan
!c  qo is mixing ratio at t+dt (kg/kg)..q1
!c
      do k = 1, km
        do i=1,im
          if (cnvflg(i) .and. k .le. kmax(i)) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
!           tvo(i,k)  = to(i,k) + delta * to(i,k) * qo(i,k)
          endif
        enddo
      enddo
!c
!c  compute moist static energy
!c
      do k = 1, km
        do i=1,im
          if (cnvflg(i) .and. k .le. kmax(i)) then
!           tem       = g * zo(i,k) + cp * to(i,k)
            tem       = phil(i,k) + cp * to(i,k)
            heo(i,k)  = tem  + hvap * qo(i,k)
            heso(i,k) = tem  + hvap * qeso(i,k)
!c           heo(i,k)  = min(heo(i,k),heso(i,k))
          endif
        enddo
      enddo
!c
!c  determine level with largest moist static energy within pbl
!c  this is the level where updraft starts
!c
      do i=1,im
         if (cnvflg(i)) then
            hmax(i) = heo(i,1)
            kb(i) = 1
         endif
      enddo
      do k = 2, km
        do i=1,im
          if (cnvflg(i).and.k.le.kpbl(i)) then
            if(heo(i,k).gt.hmax(i)) then
              kb(i)   = k
              hmax(i) = heo(i,k)
            endif
          endif
        enddo
      enddo
!c
      do k = 1, km1
        do i=1,im
          if (cnvflg(i) .and. k .le. kmax(i)-1) then
            dz      = .5 * (zo(i,k+1) - zo(i,k))
            dp      = .5 * (pfld(i,k+1) - pfld(i,k))
            es      = 0.01 * fpvs(to(i,k+1))      ! fpvs is in pa
            pprime  = pfld(i,k+1) + epsm1 * es
            qs      = eps * es / pprime
            dqsdp   = - qs / pprime
            desdt   = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
            dqsdt   = qs * pfld(i,k+1) * desdt / (es * pprime)
            gamma   = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
            dt      = (g * dz + hvap * dqsdp * dp) / (cp * (1. + gamma))
            dq      = dqsdt * dt + dqsdp * dp
            to(i,k) = to(i,k+1) + dt
            qo(i,k) = qo(i,k+1) + dq
            po(i,k) = .5 * (pfld(i,k) + pfld(i,k+1))
          endif
        enddo
      enddo
!
      do k = 1, km1
        do i=1,im
          if (cnvflg(i) .and. k .le. kmax(i)-1) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (po(i,k) + epsm1*qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
            heo(i,k)  = .5 * g * (zo(i,k) + zo(i,k+1)) +                  &
     &                  cp * to(i,k) + hvap * qo(i,k)
            heso(i,k) = .5 * g * (zo(i,k) + zo(i,k+1)) +                  &
     &                  cp * to(i,k) + hvap * qeso(i,k)
            uo(i,k)   = .5 * (uo(i,k) + uo(i,k+1))
            vo(i,k)   = .5 * (vo(i,k) + vo(i,k+1))
          endif
        enddo
      enddo
!c
!c  look for the level of free convection as cloud base
!c
      do i=1,im
        flg(i)   = cnvflg(i)
        if(flg(i)) kbcon(i) = kmax(i)
      enddo
      do k = 2, km1
        do i=1,im
          if (flg(i).and.k.lt.kbm(i)) then
            if(k.gt.kb(i).and.heo(i,kb(i)).gt.heso(i,k)) then
              kbcon(i) = k
              flg(i)   = .false.
            endif
          endif
        enddo
      enddo
!c
      do i=1,im
        if(cnvflg(i)) then
          if(kbcon(i).eq.kmax(i)) cnvflg(i) = .false.
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  determine critical convective inhibition
!c  as a function of vertical velocity at cloud base.
!c
      do i=1,im
        if(cnvflg(i)) then
          pdot(i)  = 10.* dot(i,kbcon(i))
        endif
      enddo
      do i=1,im
        if(cnvflg(i)) then
          if(slimsk(i).eq.1.) then
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          endif
          if(pdot(i).le.w4) then
            ptem = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i).ge.-w4) then
            ptem = - (pdot(i) + w4) / (w4 - w3)
          else
            ptem = 0.
          endif
          val1    =             -1.
          ptem = max(ptem,val1)
          val2    =             1.
          ptem = min(ptem,val2)
          ptem = 1. - ptem
          ptem1= .5*(cincrmax-cincrmin)
          cincr = cincrmax - ptem * ptem1
          tem1 = pfld(i,kb(i)) - pfld(i,kbcon(i))
          if(tem1.gt.cincr) then
             cnvflg(i) = .false.
          endif
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  assume the detrainment rate for the updrafts to be same as 
!c  the entrainment rate at cloud base
!c
      do i = 1, im
        if(cnvflg(i)) then
          xlamud(i) = xlamue(i,kbcon(i))
        endif
      enddo
!c
!c  determine updraft mass flux for the subcloud layers
!c
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.lt.kbcon(i).and.k.ge.kb(i)) then
              dz       = zi(i,k+1) - zi(i,k)
              ptem     = 0.5*(xlamue(i,k)+xlamue(i,k+1))-xlamud(i)
              eta(i,k) = eta(i,k+1) / (1. + ptem * dz)
            endif
          endif
        enddo
      enddo
!c
!c  compute mass flux above cloud base
!c
      do k = 2, km1
        do i = 1, im
         if(cnvflg(i))then
           if(k.gt.kbcon(i).and.k.lt.kmax(i)) then
              dz       = zi(i,k) - zi(i,k-1)
              ptem     = 0.5*(xlamue(i,k)+xlamue(i,k-1))-xlamud(i)
              eta(i,k) = eta(i,k-1) * (1 + ptem * dz)
           endif
         endif
        enddo
      enddo
!c
!c  compute updraft cloud property
!c
      do i = 1, im
        if(cnvflg(i)) then
          indx         = kb(i)
          hcko(i,indx) = heo(i,indx)
          ucko(i,indx) = uo(i,indx)
          vcko(i,indx) = vo(i,indx)
        endif
      enddo
!c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.kmax(i)) then
              dz   = zi(i,k) - zi(i,k-1)
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              ptem = 0.5 * tem + pgcon
              ptem1= 0.5 * tem - pgcon
              hcko(i,k) = ((1.-tem1)*hcko(i,k-1)+tem*0.5*                        &
     &                     (heo(i,k)+heo(i,k-1)))/factor
              ucko(i,k) = ((1.-tem1)*ucko(i,k-1)+ptem*uo(i,k)                    &
     &                     +ptem1*uo(i,k-1))/factor
              vcko(i,k) = ((1.-tem1)*vcko(i,k-1)+ptem*vo(i,k)                    &
     &                     +ptem1*vo(i,k-1))/factor
              dbyo(i,k) = hcko(i,k) - heso(i,k)
            endif
          endif
        enddo
      enddo
!c
!c   taking account into convection inhibition due to existence of
!c    dry layers below cloud base
!c
      do i=1,im
        flg(i) = cnvflg(i)
        kbcon1(i) = kmax(i)
      enddo
      do k = 2, km1
      do i=1,im
        if (flg(i).and.k.lt.kbm(i)) then
          if(k.ge.kbcon(i).and.dbyo(i,k).gt.0.) then
            kbcon1(i) = k
            flg(i)    = .false.
          endif
        endif
      enddo
      enddo
      do i=1,im
        if(cnvflg(i)) then
          if(kbcon1(i).eq.kmax(i)) cnvflg(i) = .false.
        endif
      enddo
      do i=1,im
        if(cnvflg(i)) then
          tem = pfld(i,kbcon(i)) - pfld(i,kbcon1(i))
          if(tem.gt.dthk) then
             cnvflg(i) = .false.
          endif
        endif
      enddo
!!
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  determine first guess cloud top as the level of zero buoyancy
!c    limited to the level of sigma=0.7
!c
      do i = 1, im
        flg(i) = cnvflg(i)
        if(flg(i)) ktcon(i) = kbm(i)
      enddo
      do k = 2, km1
      do i=1,im
        if (flg(i).and.k .lt. kbm(i)) then
          if(k.gt.kbcon1(i).and.dbyo(i,k).lt.0.) then
             ktcon(i) = k
             flg(i)   = .false.
          endif
        endif
      enddo
      enddo
!c
!c  turn off shallow convection if cloud top is less than pbl top
!c
     do i=1,im
       if(cnvflg(i)) then
         kk = kpbl(i)+1
         if(ktcon(i).le.kk) cnvflg(i) = .false.
       endif
     enddo
! c
! c  turn off shallow convection if cloud depth is less than
! c    a threshold value (cthk)
! c
       do i = 1, im
         if(cnvflg(i)) then
           tem = pfld(i,kbcon(i))-pfld(i,ktcon(i))
           if(tem.lt.cthk) cnvflg(i) = .false.
         endif
       enddo
!!
     totflg = .true.
     do i = 1, im
       totflg = totflg .and. (.not. cnvflg(i))
     enddo
     if(totflg) return
!!
!c
!c  specify upper limit of mass flux at cloud base
!c
      do i = 1, im
        if(cnvflg(i)) then
!         xmbmax(i) = .1
!
          k = kbcon(i)
          dp = 1000. * del(i,k)
          xmbmax(i) = dp / (g * dt2)
!
!         tem = dp / (g * dt2)
!         xmbmax(i) = min(tem, xmbmax(i))
        endif
      enddo
!c
!c  compute cloud moisture property and precipitation
!c
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = 0.
          qcko(i,kb(i)) = qo(i,kb(i))
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k)                                      &
     &             + gamma * dbyo(i,k) / (hvap * (1. + gamma))
!cj
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5*           &
     &                     (qo(i,k)+qo(i,k-1)))/factor
!cj
              dq = eta(i,k) * (qcko(i,k) - qrch)
!c
!             rhbar(i) = rhbar(i) + qo(i,k) / qeso(i,k)
!c
!c  below lfc check if there is excess moisture to release latent heat
!c
              if(k.ge.kbcon(i).and.dq.gt.0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                if(ncloud.gt.0.) then
                  dp = 1000. * del(i,k)
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                endif
                aa1(i) = aa1(i) - dz * g * qlk
                qcko(i,k)= qlk + qrch
                pwo(i,k) = etah * c0 * dz * qlk
              endif
            endif
          endif
        enddo
      enddo
!c
!c  calculate cloud work function
!c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.ge.kbcon(i).and.k.lt.ktcon(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma                 &
     &                 * to(i,k) / hvap
              aa1(i) = aa1(i) +                                &
     &                 dz1 * (g / (cp * to(i,k)))              &
     &                 * dbyo(i,k) / (1. + gamma)              &
     &                 * rfact
              val = 0.
              aa1(i)=aa1(i)+                                   &
     &                 dz1 * g * delta *                       &
     &                 max(val,(qeso(i,k) - qo(i,k)))
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i).and.aa1(i).le.0.) cnvflg(i) = .false.
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  estimate the onvective overshooting as the level
!c    where the [aafac * cloud work function] becomes zero,
!c    which is the final cloud top
!c    limited to the level of sigma=0.7
!c
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = aafac * aa1(i)
        endif
      enddo
!c
      do i = 1, im
        flg(i) = cnvflg(i)
        ktcon1(i) = kbm(i)
      enddo
      do k = 2, km1
        do i = 1, im
          if (flg(i)) then
            if(k.ge.ktcon(i).and.k.lt.kbm(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma                            &
     &                 * to(i,k) / hvap
              aa1(i) = aa1(i) +                                           &
     &                 dz1 * (g / (cp * to(i,k)))                         &
     &                 * dbyo(i,k) / (1. + gamma)                         &
     &                 * rfact
              if(aa1(i).lt.0.) then
                ktcon1(i) = k
                flg(i) = .false.
              endif
            endif
          endif
        enddo
      enddo
!c
!c  compute cloud moisture property, detraining cloud water
!c    and precipitation in overshooting layers
!c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.ge.ktcon(i).and.k.lt.ktcon1(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k)                                            &
     &             + gamma * dbyo(i,k) / (hvap * (1. + gamma))
!cj
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5*                  &
     &                     (qo(i,k)+qo(i,k-1)))/factor
!cj
              dq = eta(i,k) * (qcko(i,k) - qrch)
!c
!c  check if there is excess moisture to release latent heat
!c
              if(dq.gt.0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                if(ncloud.gt.0.) then
                  dp = 1000. * del(i,k)
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                endif
                qcko(i,k) = qlk + qrch
                pwo(i,k) = etah * c0 * dz * qlk
              endif
            endif
          endif
        enddo
      enddo
!c
!c exchange ktcon with ktcon1
!c
      do i = 1, im
        if(cnvflg(i)) then
          kk = ktcon(i)
          ktcon(i) = ktcon1(i)
          ktcon1(i) = kk
        endif
      enddo
!c
!c  this section is ready for cloud water
!c
      if(ncloud.gt.0) then
!c
!c  compute liquid and vapor separation at cloud top
!c
      do i = 1, im
        if(cnvflg(i)) then
          k = ktcon(i) - 1
          gamma = el2orc * qeso(i,k) / (to(i,k)**2)
          qrch = qeso(i,k)                                             &
     &         + gamma * dbyo(i,k) / (hvap * (1. + gamma))
          dq = qcko(i,k) - qrch
!c
!c  check if there is excess moisture to release latent heat
!c
          if(dq.gt.0.) then
            qlko_ktcon(i) = dq
            qcko(i,k) = qrch
          endif
        endif
      enddo
      endif
!!c
!c--- compute precipitation efficiency in terms of windshear
!c
      do i = 1, im
        if(cnvflg(i)) then
          vshear(i) = 0.
        endif
      enddo
      do k = 2, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              shear= sqrt((uo(i,k)-uo(i,k-1)) ** 2                       &
     &                  + (vo(i,k)-vo(i,k-1)) ** 2)
              vshear(i) = vshear(i) + shear
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          vshear(i) = 1.e3 * vshear(i) / (zi(i,ktcon(i))-zi(i,kb(i)))
          e1=1.591-.639*vshear(i)                                               &
     &       +.0953*(vshear(i)**2)-.00496*(vshear(i)**3)
          edt(i)=1.-e1
          val =         .9
          edt(i) = min(edt(i),val)
          val =         .0
          edt(i) = max(edt(i),val)
        endif
      enddo
!c
!c--- what would the change be, that a cloud with unit mass
!c--- will do to the environment?
!c
      do k = 1, km
        do i = 1, im
          if(cnvflg(i) .and. k .le. kmax(i)) then
            dellah(i,k) = 0.
            dellaq(i,k) = 0.
            dellau(i,k) = 0.
            dellav(i,k) = 0.
          endif
        enddo
      enddo
!c
!c--- changed due to subsidence and entrainment
!c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              dp = 1000. * del(i,k)
              dz = zi(i,k) - zi(i,k-1)
!c
              dv1h = heo(i,k)
              dv2h = .5 * (heo(i,k) + heo(i,k-1))
              dv3h = heo(i,k-1)
              dv1q = qo(i,k)
              dv2q = .5 * (qo(i,k) + qo(i,k-1))
              dv3q = qo(i,k-1)
              dv1u = uo(i,k)
              dv2u = .5 * (uo(i,k) + uo(i,k-1))
              dv3u = uo(i,k-1)
              dv1v = vo(i,k)
              dv2v = .5 * (vo(i,k) + vo(i,k-1))
              dv3v = vo(i,k-1)
!c
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1))
              tem1 = xlamud(i)
!cj
              dellah(i,k) = dellah(i,k) +                        &
     &     ( eta(i,k)*dv1h - eta(i,k-1)*dv3h                     &
     &    -  tem*eta(i,k-1)*dv2h*dz                              &
     &    +  tem1*eta(i,k-1)*.5*(hcko(i,k)+hcko(i,k-1))*dz       &
     &         ) *g/dp
!cj
              dellaq(i,k) = dellaq(i,k) +                        &
     &     ( eta(i,k)*dv1q - eta(i,k-1)*dv3q                     &
     &    -  tem*eta(i,k-1)*dv2q*dz                              &
     &    +  tem1*eta(i,k-1)*.5*(qcko(i,k)+qcko(i,k-1))*dz       &
     &         ) *g/dp
!cj
              dellau(i,k) = dellau(i,k) +                        &
     &     ( eta(i,k)*dv1u - eta(i,k-1)*dv3u                     &
     &    -  tem*eta(i,k-1)*dv2u*dz                              &
     &    +  tem1*eta(i,k-1)*.5*(ucko(i,k)+ucko(i,k-1))*dz       &
     &    -  pgcon*eta(i,k-1)*(dv1u-dv3u)                        &
     &         ) *g/dp
!cj
              dellav(i,k) = dellav(i,k) +                        &
     &     ( eta(i,k)*dv1v - eta(i,k-1)*dv3v                     &
     &    -  tem*eta(i,k-1)*dv2v*dz                              &
     &    +  tem1*eta(i,k-1)*.5*(vcko(i,k)+vcko(i,k-1))*dz       &
     &    -  pgcon*eta(i,k-1)*(dv1v-dv3v)                        &
     &         ) *g/dp
!cj
            endif
          endif
        enddo
      enddo
!c
!c------- cloud top
!c
      do i = 1, im
        if(cnvflg(i)) then
          indx = ktcon(i)
          dp = 1000. * del(i,indx)
          dv1h = heo(i,indx-1)
          dellah(i,indx) = eta(i,indx-1) *                      &
     &                     (hcko(i,indx-1) - dv1h) * g / dp
          dv1q = qo(i,indx-1)
          dellaq(i,indx) = eta(i,indx-1) *                      &
     &                     (qcko(i,indx-1) - dv1q) * g / dp
          dv1u = uo(i,indx-1)
          dellau(i,indx) = eta(i,indx-1) *                      &
     &                     (ucko(i,indx-1) - dv1u) * g / dp
          dv1v = vo(i,indx-1)
          dellav(i,indx) = eta(i,indx-1) *                      &
     &                     (vcko(i,indx-1) - dv1v) * g / dp
!c
!c  cloud water
!c
          dellal(i,indx) = eta(i,indx-1) *                      &
     &                     qlko_ktcon(i) * g / dp
        endif
      enddo
!c
!c  mass flux at cloud base for shallow convection
!c  (Grant, 2001)
!c
      do i= 1, im
        if(cnvflg(i)) then
          k = kbcon(i)
!         ptem = g*sflx(i)*zi(i,k)/t1(i,1)
          ptem = g*sflx(i)*hpbl(i)/t1(i,1)
          wstar(i) = ptem**h1
          tem = po(i,k)*100. / (rd*t1(i,k))
          xmb(i) = betaw*tem*wstar(i)
          xmb(i) = min(xmb(i),xmbmax(i))
        endif
      enddo
!c
!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)) then
            qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val     =             1.e-8
            qeso(i,k) = max(qeso(i,k), val )
          endif
        enddo
      enddo
!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c
      do i = 1, im
        delhbar(i) = 0.
        delqbar(i) = 0.
        deltbar(i) = 0.
        delubar(i) = 0.
        delvbar(i) = 0.
        qcond(i) = 0.
      enddo
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cp
              t1(i,k) = t1(i,k) + dellat * xmb(i) * dt2
              q1(i,k) = q1(i,k) + dellaq(i,k) * xmb(i) * dt2
              tem = 1./rcs(i)
              u1(i,k) = u1(i,k) + dellau(i,k) * xmb(i) * dt2 * tem
              v1(i,k) = v1(i,k) + dellav(i,k) * xmb(i) * dt2 * tem
              dp = 1000. * del(i,k)
              delhbar(i) = delhbar(i) + dellah(i,k)*xmb(i)*dp/g
              delqbar(i) = delqbar(i) + dellaq(i,k)*xmb(i)*dp/g
              deltbar(i) = deltbar(i) + dellat*xmb(i)*dp/g
              delubar(i) = delubar(i) + dellau(i,k)*xmb(i)*dp/g
              delvbar(i) = delvbar(i) + dellav(i,k)*xmb(i)*dp/g
            endif
          endif
        enddo
      enddo
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
              qeso(i,k) = eps * qeso(i,k)/(pfld(i,k) + epsm1*qeso(i,k))
              val     =             1.e-8
              qeso(i,k) = max(qeso(i,k), val )
            endif
          endif
        enddo
      enddo
!c
      do i = 1, im
        rntot(i) = 0.
        delqev(i) = 0.
        delq2(i) = 0.
        flg(i) = cnvflg(i)
      enddo
      do k = km, 1, -1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.lt.ktcon(i).and.k.gt.kb(i)) then
              rntot(i) = rntot(i) + pwo(i,k) * xmb(i) * .001 * dt2
            endif
          endif
        enddo
      enddo
!c
!c evaporating rain
!c
      do k = km, 1, -1
        do i = 1, im
          if (k .le. kmax(i)) then
            deltv(i) = 0.
            delq(i) = 0.
            qevap(i) = 0.
            if(cnvflg(i)) then
              if(k.lt.ktcon(i).and.k.gt.kb(i)) then
                rn(i) = rn(i) + pwo(i,k) * xmb(i) * .001 * dt2
              endif
            endif
            if(flg(i).and.k.lt.ktcon(i)) then
              evef = edt(i) * evfact
              if(slimsk(i).eq.1.) evef=edt(i) * evfactl
!             if(slimsk(i).eq.1.) evef=.07
!c             if(slimsk(i).ne.1.) evef = 0.
              qcond(i) = evef * (q1(i,k) - qeso(i,k))                            &
     &                 / (1. + el2orc * qeso(i,k) / t1(i,k)**2)
              dp = 1000. * del(i,k)
              if(rn(i).gt.0..and.qcond(i).lt.0.) then
                qevap(i) = -qcond(i) * (1.-exp(-.32*sqrt(dt2*rn(i))))
                qevap(i) = min(qevap(i), rn(i)*1000.*g/dp)
                delq2(i) = delqev(i) + .001 * qevap(i) * dp / g
              endif
              if(rn(i).gt.0..and.qcond(i).lt.0..and.                            &
     &           delq2(i).gt.rntot(i)) then
                qevap(i) = 1000.* g * (rntot(i) - delqev(i)) / dp
                flg(i) = .false.
              endif
              if(rn(i).gt.0..and.qevap(i).gt.0.) then
                tem  = .001 * dp / g
                tem1 = qevap(i) * tem
                if(tem1.gt.rn(i)) then
                  qevap(i) = rn(i) / tem
                  rn(i) = 0.
                else
                  rn(i) = rn(i) - tem1
                endif
                q1(i,k) = q1(i,k) + qevap(i)
                t1(i,k) = t1(i,k) - elocp * qevap(i)
                deltv(i) = - elocp*qevap(i)/dt2
                delq(i) =  + qevap(i)/dt2
                delqev(i) = delqev(i) + .001*dp*qevap(i)/g
              endif
              dellaq(i,k) = dellaq(i,k) + delq(i) / xmb(i)
              delqbar(i) = delqbar(i) + delq(i)*dp/g
              deltbar(i) = deltbar(i) + deltv(i)*dp/g
            endif
          endif
        enddo
      enddo
!cj
!     do i = 1, im
!     if(me.eq.31.and.cnvflg(i)) then
!     if(cnvflg(i)) then
!       print *, ' shallow delhbar, delqbar, deltbar = ',
!    &             delhbar(i),hvap*delqbar(i),cp*deltbar(i)
!       print *, ' shallow delubar, delvbar = ',delubar(i),delvbar(i)
!       print *, ' precip =', hvap*rn(i)*1000./dt2
!       print*,'pdif= ',pfld(i,kbcon(i))-pfld(i,ktcon(i))
!     endif
!     enddo
!cj
      do i = 1, im
        if(cnvflg(i)) then
          if(rn(i).lt.0..or..not.flg(i)) rn(i) = 0.
          ktop(i) = ktcon(i)
          kbot(i) = kbcon(i)
          kcnv(i) = 0
        endif
      enddo
!c
!c  cloud water
!c
      if (ncloud.gt.0) then
!
      val1 = 1.0
      val2 = 0.
      do k = 1, km1
        do i = 1, im
          if (cnvflg(i)) then
            if (k.gt.kb(i).and.k.le.ktcon(i)) then
              tem  = dellal(i,k) * xmb(i) * dt2
!             tem1 = max(0.0,  min(1.0,  (tcr-t1(i,k))*tcrf))
              tem1 = max(val2, min(val1, (tcr-t1(i,k))*tcrf))
              if (ql(i,k,2) .gt. -999.0) then
                ql(i,k,1) = ql(i,k,1) + tem * tem1            ! ice
                ql(i,k,2) = ql(i,k,2) + tem *(1.0-tem1)       ! water
              else
                ql(i,k,1) = ql(i,k,1) + tem
              endif
            endif
          endif
        enddo
      enddo
!
      endif
!
! hchuang code change
!
      do k = 1, km
        do i = 1, im
          if(cnvflg(i)) then
            if(k.ge.kb(i) .and. k.lt.ktop(i)) then
              ud_mf(i,k) = eta(i,k) * xmb(i) * dt2
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
           k = ktop(i)-1
           dt_mf(i,k) = ud_mf(i,k)
        endif
      enddo
!!
      return
    end subroutine shalcnv_hur


!!!!MESO SAS
      subroutine sascnvn_meso(im,ix,km,jcap,delt,del,prsl,ps,phil,ql, &
     &     q1,t1,u1,v1,rcs,cldwrk,rn,kbot,ktop,kcnv,slimsk,        &
     &     dot,ncloud,pgcon,sas_mass_flux)
!     &     dot,ncloud,pgcon,sas_mass_flux,sigma,jqfliu)
!     &     dot,ncloud,sigma,pgcon,sas_mass_flux)
!    &     dot,ncloud,ud_mf,dd_mf,dt_mf,me)
!
! Version 20120809
!  Modified on 20120803 to add dbyod, include definition of heotd to jmin level, and fix bug
!   on the calculation of qotd
!  Modified on 20120807 to fix bug in the dhdt calculation
!
!  Adding in consistency with the pwo, pwdo so rain is consistent with heating and drying
!  after the tilda terms are computed.
!
!  20120822
!   Turns off SAS when sigma is greater than .9
!   Correct cloud top cloud water detrainment
!
!      use machine , only : kind_phys
!      use funcphys , only : fpvs
!      use physcons, grav => con_g, cp => con_cp, hvap => con_hvap  &
      USE MACHINE, ONLY : kind_phys
      USE FUNCPHYS, ONLY : fpvs
      USE PHYSCONS, grav => con_g, cp => con_cp         &
     &,             hvap => con_hvap                               &
     &,             rv => con_rv, fv => con_fvirt, t0c => con_t0c  &
     &,             cvap => con_cvap, cliq => con_cliq             &
     &,             eps => con_eps, epsm1 => con_epsm1             &
     &,             rd => con_rd
      implicit none
!
      integer            im, ix,  km, jcap, ncloud,                &
     &                   kbot(im), ktop(im), kcnv(im),jqfliu
!    &,                  me
      real(kind=kind_phys) delt, sas_mass_flux
      real(kind=kind_phys) ps(im),     del(ix,km),  prsl(ix,km),   &
     &                     ql(ix,km,2),q1(ix,km),   t1(ix,km),     &
     &                     u1(ix,km),  v1(ix,km),   rcs(im),       &
     &                     cldwrk(im), rn(im),      slimsk(im),    &
     &                     dot(ix,km), phil(ix,km)                 &
! hchuang code change mass flux output
     &,                    ud_mf(im,km),dd_mf(im,km),dt_mf(im,km)
!
      integer              i, j, indx, jmn, k, kk, latd, lond, km1
!
      real(kind=kind_phys) clam, cxlamu, xlamde, xlamdd
!
      real(kind=kind_phys) adw,     aup,     aafac,                &
     &                     beta,    betal,   betas,                &
     &                     c0,      cpoel,   dellat,  delta,       &
     &                     desdt,   deta,    detad,   dg,          &
     &                     dh,      dhh,     dlnsig,  dp,          &
     &                     dq,      dqsdp,   dqsdt,   dt,          &
     &                     dt2,     dtmax,   dtmin,   dv1h,        &
     &                     dv1q,    dv2h,    dv2q,    dv1u,        &
     &                     dv1v,    dv2u,    dv2v,    dv3q,        &
     &                     dv3h,    dv3u,    dv3v,                 &
     &                     dv1hd,   dv1qd,   dv2hd,   dv2qd,       &
     &                     dv1ud,   dv1vd,   dv2ud,   dv2vd,       &
     &                     dv3hd,   dv3qd,   dv3ud,   dv3vd,       &
     &                     dz,      dz1,     e1,      edtmax,      &
     &                     edtmaxl, edtmaxs, el2orc,  elocp,       &
     &                     es,      etah,    cthk,    dthk,        &
     &                     evef,    evfact,  evfactl, fact1,       &
     &                     fact2,   factor,  fjcap,   fkm,         &
     &                     g,       gamma,   pprime,               &
     &                     qlk,     qrch,    qs,      c1,          &
     &                     rain,    rfact,   shear,   tem1,        &
     &                     tem2,    terr,    val,     val1,        &
     &                     val2,    w1,      w1l,     w1s,         &
     &                     w2,      w2l,     w2s,     w3,          &
     &                     w3l,     w3s,     w4,      w4l,         &
     &                     w4s,     xdby,    xpw,     xpwd,        &
     &                     xqrch,   armb,    ardt,    mbdt,        &
     &                     delhz,   delqz,   deluz,   delvz,       &
     &                     tem,     ptem,    ptem1
!
      real(kind=kind_phys), intent(in) :: pgcon
!
      integer              kb(im), kbcon(im), kbcon1(im),          &
     &                     ktcon(im), ktcon1(im),                  &
     &                     jmin(im), lmin(im), kbmax(im),          &
     &                     kbm(im), kmax(im)
!
      real(kind=kind_phys) aa1(im),     acrt(im),   acrtfct(im),   &
     &                     delhbar(im), delq(im),   delq2(im),     &
     &                     delqbar(im), delqev(im), deltbar(im),   &
     &                     deltv(im),   dtconv(im), edt(im),       &
     &                     edto(im),    edtx(im),   fld(im),       &
     &                     hcdo(im,km), hmax(im),   hmin(im),      &
     &                     ucdo(im,km), vcdo(im,km),aa2(im),       &
     &                     pbcdif(im),  pdot(im),   po(im,km),     &
     &                     pwavo(im),   pwevo(im),  xlamud(im),    &
     &                     qcdo(im,km), qcond(im),  qevap(im),     &
     &                     rntot(im),   vshear(im), xaa0(im),      &
     &                     xk(im),      xlamd(im),                 &
     &                     xmb(im),     xmbmax(im), xpwav(im),     &
     &                     xpwev(im),   delubar(im),delvbar(im)
!cj
      real(kind=kind_phys) cincr, cincrmax, cincrmin
!cj
!c  physical parameters
      parameter(g=grav)
      parameter(cpoel=cp/hvap,elocp=hvap/cp,                       &
     &          el2orc=hvap*hvap/(rv*cp))
      parameter(terr=0.,c0=.002,c1=.002,delta=fv)
      parameter(fact1=(cvap-cliq)/rv,fact2=hvap/rv-fact1*t0c)
      parameter(cthk=150.,cincrmax=180.,cincrmin=120.,dthk=25.)
! Qingfu modified
!      parameter(cthk=150.,cincrmax=160.,cincrmin=100.,dthk=25.)
!c  local variables and arrays
      real(kind=kind_phys) pfld(im,km),to(im,km), qo(im,km),       &
     &                     uo(im,km),  vo(im,km), qeso(im,km)
!c  cloud water
      real(kind=kind_phys)qlko_ktcon(im),dellal(im,km),tvo(im,km), &
     &                dbyo(im,km),  zo(im,km),     xlamue(im,km),  &
     &                fent1(im,km), fent2(im,km),  frh(im,km),     &
     &                heo(im,km),   heso(im,km),   doto(im,km-1),  &
!c  heotu and qeotu are the environmental mean h and q for updraft
!c  heotd and qeotd are the environmental mean h and q for downdraft
     &                heotu(im,km), qotu(im,km),   uotu(im,km),    &
     &                votu(im,km),  heotd(im,km),  qotd(im,km),    &
     &                uotd(im,km),  votd(im,km),                   &
     &                delhx(im,km), delqx(im,km),                  &
     &                delux(im,km), delvx(im,km),                  &
     &                qrcd(im,km),  dellah(im,km), dellaq(im,km),  &
     &                dellau(im,km),dellav(im,km), hcko(im,km),    &
     &                ucko(im,km),  vcko(im,km),   qcko(im,km),    &
     &                eta(im,km),   etad(im,km),   zi(im,km),      &
     &                qrcdo(im,km), pwo(im,km),    pwdo(im,km),    &
     &                wc(im),       wbar(im),                      &
     &                sigma(im),    sigi1(im),     sigi2(im),      &
     &                tx1(im),      sumx(im),      dbyod(im,km)
!    &,               rhbar(im)
!
      logical totflg, cnvflg(im), cnvdflg(im), flg(im)
!
      real(kind=kind_phys) pcrit(15), acritt(15), acrit(15)
!     save pcrit, acritt
      data pcrit/850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,  &
     &           350.,300.,250.,200.,150./
      data acritt/.0633,.0445,.0553,.0664,.075,.1082,.1521,.2216,    &
     &           .3151,.3677,.41,.5255,.7663,1.1686,1.6851/
!c  gdas derived acrit
!c     data acritt/.203,.515,.521,.566,.625,.665,.659,.688,
!c    &            .743,.813,.886,.947,1.138,1.377,1.896/
      real(kind=kind_phys) tf, tcr, tcrf
      parameter (tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf))

      real    sigma_sum
!
!c-----------------------------------------------------------------------
!
      km1 = km - 1
!c
!c  initialize arrays
!c
      do i=1,im
        kcnv(i)=0
        cnvflg(i) = .true.
        rn(i)=0.
        kbot(i)=km+1
        ktop(i)=0
        kbcon(i)=km
        ktcon(i)=1
        dtconv(i) = 3600.
        cldwrk(i) = 0.
        pdot(i) = 0.
        pbcdif(i)= 0.
        lmin(i) = 1
        jmin(i) = 1
        qlko_ktcon(i) = 0.
        edt(i)  = 0.
        edto(i) = 0.
        edtx(i) = 0.
        acrt(i) = 0.
        acrtfct(i) = 1.
        aa1(i)  = 0.
        aa2(i)  = 0.
        xaa0(i) = 0.
        pwavo(i)= 0.
        pwevo(i)= 0.
        xpwav(i)= 0.
        xpwev(i)= 0.
        vshear(i) = 0.
        wc(i) = 0.
        wbar(i) = 0.
        xmb(i) = 0.
      enddo
! hchuang code change
      do k = 1, km
        do i = 1, im
          ud_mf(i,k) = 0.
          dd_mf(i,k) = 0.
          dt_mf(i,k) = 0.
        enddo
      enddo
!c
      do k = 1, 15
        acrit(k) = acritt(k) * (975. - pcrit(k))
      enddo
      dt2 = delt
      val   =         1200.
      dtmin = max(dt2, val )
      val   =         3600.
      dtmax = max(dt2, val )
!c  model tunable parameters are all here
!      mbdt    = 10.
      armb    = 1.             ! arbitrary cloud base mass flux
      ardt    = 10.            ! arbitrary time step
      mbdt    = armb * ardt
      mbdt    = min(mbdt, dt2)
      edtmaxl = .3
      edtmaxs = .3
      clam    = .1
      aafac   = .1
!     betal   = .15
!     betas   = .15
      betal   = .05
      betas   = .05
!c     evef    = 0.07
      evfact  = 0.3
      evfactl = 0.3
!
      cxlamu  = 1.0e-4
      xlamde  = 1.0e-4
      xlamdd  = 1.0e-4
!
      fjcap   = (float(jcap) / 126.) ** 2
      val     =           1.
      fjcap   = max(fjcap,val)
      fkm     = (float(km) / 28.) ** 2
      fkm     = max(fkm,val)
      w1l     = -8.e-3
      w2l     = -4.e-2
      w3l     = -5.e-3
      w4l     = -5.e-4
      w1s     = -2.e-4
      w2s     = -2.e-3
      w3s     = -1.e-3
      w4s     = -2.e-5
!c
!c  define top layer for search of the downdraft originating layer
!c  and the maximum thetae for updraft
!c
      do i=1,im
        kbmax(i) = km
        kbm(i)   = km
        kmax(i)  = km
        tx1(i)   = 1.0 / ps(i)
      enddo
!
      do k = 1, km1
        do i=1,im
          if (prsl(i,k)*tx1(i) .gt. 0.04) kmax(i)  = k + 1
          if (prsl(i,k)*tx1(i) .gt. 0.45) kbmax(i) = k + 1
          if (prsl(i,k)*tx1(i) .gt. 0.70) kbm(i)   = k + 1
        enddo
      enddo
      do i=1,im
        kbmax(i) = min(kbmax(i),kmax(i))
        kbm(i)   = min(kbm(i),kmax(i))
        kmax(i) = min(km,kmax(i))
      enddo
!c
!c  hydrostatic height assume zero terr and initially assume
!c    updraft entrainment rate as an inverse function of height
!c
      do k = 1, km
        do i=1,im
          zo(i,k) = phil(i,k) / g
        enddo
      enddo
      do k = 1, km1
        do i=1,im
          zi(i,k) = 0.5*(zo(i,k)+zo(i,k+1))
          xlamue(i,k) = clam / zi(i,k)
        enddo
      enddo
!c
!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c   convert surface pressure to mb from cb
!c
      do k = 1, km
        do i = 1, im
          if (k .le. kmax(i)) then
            pfld(i,k) = prsl(i,k) * 10.0
            eta(i,k)  = 1.
            fent1(i,k)= 1.
            fent2(i,k)= 1.
            frh(i,k)  = 0.
            hcko(i,k) = 0.
            qcko(i,k) = 0.
            ucko(i,k) = 0.
            vcko(i,k) = 0.
            etad(i,k) = 1.
            hcdo(i,k) = 0.
            qcdo(i,k) = 0.
            ucdo(i,k) = 0.
            vcdo(i,k) = 0.
            qrcd(i,k) = 0.
            qrcdo(i,k)= 0.
            dbyo(i,k) = 0.
            dbyod(i,k) = 0.
            pwo(i,k)  = 0.
            pwdo(i,k) = 0.
            dellal(i,k) = 0.
            to(i,k)   = t1(i,k)
            qo(i,k)   = q1(i,k)
            uo(i,k)   = u1(i,k) * rcs(i)
            vo(i,k)   = v1(i,k) * rcs(i)
            delhx(i,k) = 0.
            delqx(i,k) = 0.
            delux(i,k) = 0.
            delvx(i,k) = 0.
          endif
        enddo
      enddo
!c
!c  column variables
!c  p is pressure of the layer (mb)
!c  t is temperature at t-dt (k)..tn
!c  q is mixing ratio at t-dt (kg/kg)..qn
!c  to is temperature at t+dt (k)... this is after advection and turbulan
!c  qo is mixing ratio at t+dt (kg/kg)..q1
!c
      do k = 1, km
        do i=1,im
          if (k .le. kmax(i)) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
!           tvo(i,k)  = to(i,k) + delta * to(i,k) * qo(i,k)
          endif
        enddo
      enddo
!c
!c  compute moist static energy
!c
      do k = 1, km
        do i=1,im
          if (k .le. kmax(i)) then
!           tem       = g * zo(i,k) + cp * to(i,k)
            tem       = phil(i,k) + cp * to(i,k)
            heo(i,k)  = tem  + hvap * qo(i,k)
            heso(i,k) = tem  + hvap * qeso(i,k)
!c           heo(i,k)  = min(heo(i,k),heso(i,k))
          endif
        enddo
      enddo
!c
!c  determine level with largest moist static energy
!c  this is the level where updraft starts
!c
      do i=1,im
        hmax(i) = heo(i,1)
        kb(i)   = 1
      enddo
      do k = 2, km
        do i=1,im
          if (k .le. kbm(i)) then
            if(heo(i,k).gt.hmax(i)) then
              kb(i)   = k
              hmax(i) = heo(i,k)
            endif
          endif
        enddo
      enddo
!c
      do k = 1, km1
        do i=1,im
          if (k .le. kmax(i)-1) then
            dz      = .5 * (zo(i,k+1) - zo(i,k))
            dp      = .5 * (pfld(i,k+1) - pfld(i,k))
            es      = 0.01 * fpvs(to(i,k+1))      ! fpvs is in pa
            pprime  = pfld(i,k+1) + epsm1 * es
            qs      = eps * es / pprime
            dqsdp   = - qs / pprime
            desdt   = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
            dqsdt   = qs * pfld(i,k+1) * desdt / (es * pprime)
            gamma   = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
            dt      = (g * dz + hvap * dqsdp * dp) / (cp * (1. + gamma))
            dq      = dqsdt * dt + dqsdp * dp
            to(i,k) = to(i,k+1) + dt
            qo(i,k) = qo(i,k+1) + dq
            po(i,k) = .5 * (pfld(i,k) + pfld(i,k+1))
          endif
        enddo
      enddo
!
      do k = 1, km1
        do i=1,im
          if (k .le. kmax(i)-1) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (po(i,k) + epsm1*qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
            val1      = 1.
            frh(i,k)  = 1. - min(qo(i,k)/qeso(i,k), val1)
            heo(i,k)  = .5 * g * (zo(i,k) + zo(i,k+1)) +            &
     &                  cp * to(i,k) + hvap * qo(i,k)
            heso(i,k) = .5 * g * (zo(i,k) + zo(i,k+1)) +            &
     &                  cp * to(i,k) + hvap * qeso(i,k)
            uo(i,k)   = .5 * (uo(i,k) + uo(i,k+1))
            vo(i,k)   = .5 * (vo(i,k) + vo(i,k+1))
            doto(i,k) = 1000. * dot(i,k+1)          ! pa/s
          endif
        enddo
      enddo
!c
!c  initialize environmental property as grid mean value
!c
      do k = 1, km
        do i=1,im
          if (k .le. kmax(i)) then
            heotu(i,k) = heo(i,k)
            qotu(i,k) = qo(i,k)
            uotu(i,k) = uo(i,k)
            votu(i,k) = vo(i,k)
            heotd(i,k) = heo(i,k)
            qotd(i,k) = qo(i,k)
            uotd(i,k) = uo(i,k)
            votd(i,k) = vo(i,k)
          endif
        enddo
      enddo
!c
!c  look for the level of free convection as cloud base
!c
      do i=1,im
        flg(i)   = .true.
        kbcon(i) = kmax(i)
      enddo
      do k = 1, km1
        do i=1,im
          if (flg(i).and.k.le.kbmax(i)) then
            if(k.gt.kb(i).and.heo(i,kb(i)).gt.heso(i,k)) then
              kbcon(i) = k
              flg(i)   = .false.
            endif
          endif
        enddo
      enddo
!c
      do i=1,im
        if(kbcon(i).eq.kmax(i)) cnvflg(i) = .false.
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  determine critical convective inhibition
!c  as a function of vertical velocity at cloud base.
!c
      do i=1,im
        if(cnvflg(i)) then
          pdot(i)  = 10.* dot(i,kbcon(i))
        endif
      enddo
      do i=1,im
        if(cnvflg(i)) then
          if(slimsk(i).eq.1.) then
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          endif
          if(pdot(i).le.w4) then
            tem = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i).ge.-w4) then
            tem = - (pdot(i) + w4) / (w4 - w3)
          else
            tem = 0.
          endif
          val1    =             -1.
          tem = max(tem,val1)
          val2    =             1.
          tem = min(tem,val2)
          tem = 1. - tem
          tem1= .5*(cincrmax-cincrmin)
          cincr = cincrmax - tem * tem1
          pbcdif(i) = pfld(i,kb(i)) - pfld(i,kbcon(i))
          if(pbcdif(i).gt.cincr) then
             cnvflg(i) = .false.
          endif
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  assume that updraft entrainment rate above cloud base is
!c    same as that at cloud base
!c
      do k = 2, km1
        do i=1,im
          if(cnvflg(i).and.                                 &
     &      (k.gt.kbcon(i).and.k.lt.kmax(i))) then
              xlamue(i,k) = xlamue(i,kbcon(i))
          endif
        enddo
      enddo
!c
!c  assume the detrainment rate for the updrafts to be same as
!c  the entrainment rate at cloud base
!c
      do i = 1, im
        if(cnvflg(i)) then
          xlamud(i) = xlamue(i,kbcon(i))
        endif
      enddo
!c
!c  functions rapidly decreasing with height, mimicking a cloud ensemble
!c    (Bechtold et al., 2008)
!c
      val1=1.0
      do k = 2, km1
        do i=1,im
          if(cnvflg(i).and.                                &
     &      (k.gt.kbcon(i).and.k.lt.kmax(i))) then
!              tem = qeso(i,k)/qeso(i,kbcon(i))
              tem = min(val1,qeso(i,k)/qeso(i,kbcon(i)))
              fent1(i,k) = tem**2
              fent2(i,k) = tem**3
          endif
        enddo
      enddo
!c
!c  final entrainment rate as the sum of turbulent part and organized entrainment
!c    depending on the environmental relative humidity
!c    (Bechtold et al., 2008)
!c
      do k = 2, km1
        do i=1,im
          if(cnvflg(i).and.                                &
     &      (k.ge.kbcon(i).and.k.lt.kmax(i))) then
              tem = cxlamu * frh(i,k) * fent2(i,k)
              xlamue(i,k) = xlamue(i,k)*fent1(i,k) + tem
          endif
        enddo
      enddo
!c
!c  determine updraft mass flux for the subcloud layers
!c
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.lt.kbcon(i).and.k.ge.kb(i)) then
              dz       = zi(i,k+1) - zi(i,k)
              ptem     = 0.5*(xlamue(i,k)+xlamue(i,k+1))-xlamud(i)
              eta(i,k) = eta(i,k+1) / (1. + ptem * dz)
            endif
          endif
        enddo
      enddo
!c
!c  compute mass flux above cloud base
!c
      do k = 2, km1
        do i = 1, im
         if(cnvflg(i))then
           if(k.gt.kbcon(i).and.k.lt.kmax(i)) then
              dz       = zi(i,k) - zi(i,k-1)
              ptem     = 0.5*(xlamue(i,k)+xlamue(i,k-1))-xlamud(i)
              eta(i,k) = eta(i,k-1) * (1 + ptem * dz)
           endif
         endif
        enddo
      enddo
!c
!c  compute updraft cloud properties
!c
      do i = 1, im
        if(cnvflg(i)) then
          indx         = kb(i)
          hcko(i,indx) = heo(i,indx)
          ucko(i,indx) = uo(i,indx)
          vcko(i,indx) = vo(i,indx)
          pwavo(i)     = 0.
        endif
      enddo
!c
!c  cloud property is modified by the entrainment process
!c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.kmax(i)) then
              dz   = zi(i,k) - zi(i,k-1)
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              ptem = 0.5 * tem + pgcon
              ptem1= 0.5 * tem - pgcon
              hcko(i,k) = ((1.-tem1)*hcko(i,k-1)+tem*0.5*       &
     &                     (heo(i,k)+heo(i,k-1)))/factor
              ucko(i,k) = ((1.-tem1)*ucko(i,k-1)+ptem*uo(i,k)   &
     &                     +ptem1*uo(i,k-1))/factor
              vcko(i,k) = ((1.-tem1)*vcko(i,k-1)+ptem*vo(i,k)   &
     &                     +ptem1*vo(i,k-1))/factor
              dbyo(i,k) = hcko(i,k) - heso(i,k)
            endif
          endif
        enddo
      enddo
!c
!c   taking account into convection inhibition due to existence of
!c    dry layers below cloud base
!c
      do i=1,im
        flg(i) = cnvflg(i)
        kbcon1(i) = kmax(i)
      enddo
      do k = 2, km1
      do i=1,im
        if (flg(i).and.k.lt.kmax(i)) then
          if(k.ge.kbcon(i).and.dbyo(i,k).gt.0.) then
            kbcon1(i) = k
            flg(i)    = .false.
          endif
        endif
      enddo
      enddo
      do i=1,im
        if(cnvflg(i)) then
          if(kbcon1(i).eq.kmax(i)) cnvflg(i) = .false.
        endif
      enddo
      do i=1,im
        if(cnvflg(i)) then
          tem = pfld(i,kbcon(i)) - pfld(i,kbcon1(i))
          if(tem.gt.dthk) then
             cnvflg(i) = .false.
          endif
        endif
      enddo
!!
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  determine first guess cloud top as the level of zero buoyancy
!c
      do i = 1, im
        flg(i) = cnvflg(i)
        ktcon(i) = 1
      enddo
      do k = 2, km1
      do i = 1, im
        if (flg(i).and.k .lt. kmax(i)) then
          if(k.gt.kbcon1(i).and.dbyo(i,k).lt.0.) then
             ktcon(i) = k
             flg(i)   = .false.
          endif
        endif
      enddo
      enddo
!c
      do i = 1, im
        if(cnvflg(i)) then
          tem = pfld(i,kbcon(i))-pfld(i,ktcon(i))
          if(tem.lt.cthk) cnvflg(i) = .false.
        endif
      enddo
!!
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  search for downdraft originating level above theta-e minimum
!c
      do i = 1, im
        if(cnvflg(i)) then
           hmin(i) = heo(i,kbcon1(i))
           lmin(i) = kbmax(i)
           jmin(i) = kbmax(i)
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i) .and. k .le. kbmax(i)) then
            if(k.gt.kbcon1(i).and.heo(i,k).lt.hmin(i)) then
               lmin(i) = k + 1
               hmin(i) = heo(i,k)
            endif
          endif
        enddo
      enddo
!c
!c  make sure that jmin(i) is within the cloud
!c
      do i = 1, im
        if(cnvflg(i)) then
          jmin(i) = min(lmin(i),ktcon(i)-1)
          jmin(i) = max(jmin(i),kbcon1(i)+1)
          if(jmin(i).ge.ktcon(i)) cnvflg(i) = .false.
        endif
      enddo
!c
!c  specify upper limit of mass flux at cloud base
!c
      do i = 1, im
        if(cnvflg(i)) then
!         xmbmax(i) = .1
!
          k = kbcon(i)
          dp = 1000. * del(i,k)
          xmbmax(i) = dp / (g * dt2)
          xmbmax(i) = min(sas_mass_flux,xmbmax(i))
!
!         tem = dp / (g * dt2)
!         xmbmax(i) = min(tem, xmbmax(i))
        endif
      enddo
!c
!c  compute cloud moisture property and precipitation
!c
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = 0.
          indx = kb(i)
          qcko(i,indx) = qo(i,indx)
          qcko(i,1) = qcko(i,indx)
!         rhbar(i) = 0.
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k)                                  &
     &             + gamma * dbyo(i,k) / (hvap * (1. + gamma))
!cj
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5*       &
     &                     (qo(i,k)+qo(i,k-1)))/factor
!cj
              dq = eta(i,k) * (qcko(i,k) - qrch)
!c
!             rhbar(i) = rhbar(i) + qo(i,k) / qeso(i,k)
!c
!c  check if there is excess moisture to release latent heat
!c
              if(k.ge.kbcon(i).and.dq.gt.0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                if(ncloud.gt.0.and.k.gt.jmin(i)) then
                  dp = 1000. * del(i,k)
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                endif
                aa1(i) = aa1(i) - dz * g * qlk
                qcko(i,k) = qlk + qrch
                pwo(i,k) = etah * c0 * dz * qlk
                pwavo(i) = pwavo(i) + pwo(i,k)
              endif
            endif
          endif
        enddo
      enddo
!c
!     do i = 1, im
!       if(cnvflg(i)) then
!         indx = ktcon(i) - kb(i) - 1
!         rhbar(i) = rhbar(i) / float(indx)
!       endif
!     enddo
!c
!c  calculate cloud work function
!c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.ge.kbcon(i).and.k.lt.ktcon(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma                &
     &                 * to(i,k) / hvap
              aa1(i) = aa1(i) +                               &
     &                 dz1 * (g / (cp * to(i,k)))             &
     &                 * dbyo(i,k) / (1. + gamma)             &
     &                 * rfact
              val = 0.
              aa1(i)=aa1(i)+                                  &
     &                 dz1 * g * delta *                      &
     &                 max(val,(qeso(i,k) - qo(i,k)))
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i).and.aa1(i).le.0.) cnvflg(i) = .false.
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  calculate the updraft area sigma as a function of the updraft speed wc=sqrt(2*aa1 + wbar**2 (at cloud base))
!c  and the area mean vertical wind speed wbar = -omega / (rho * g)
!c
!c  po is in the unit of mb and dot in the unit of cb/sec
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.ge.kbcon(i).and.k.lt.ktcon(i)) then
              tem = po(i,k) / (rd * to(i,k) * (1. + delta * qo(i,k)))
              tem1 = - 10. * dot(i,k) / (tem * g)
              wbar(i) = max(wbar(i),tem1)
            endif
          endif
        enddo
      enddo
!
!
!c   cloud base updraft speed is added here. For the time being, we use the same wbar as above. This guarantee that
!c   the calculated sigma never exceeds one.
!
      do i = 1, im
        if(cnvflg(i)) then
!          k = kbcon(i)
!          tem = po(i,k) / (rd * to(i,k) * (1. + delta * qo(i,k)))
!          tem1 = - 10. * dot(i,k) / (tem * g)
!          wc(i) = sqrt(tem1*tem1+2.*aa1(i))
          wc(i) = sqrt(wbar(i) * wbar(i) + 2. * aa1(i))
        endif
      enddo
      sigma_sum=0.
      val1=0.09
!      val1=0.0
      do i = 1, im
        if(cnvflg(i).and.wc(i).gt.0.) then
!
!  Scale sigma assuming magnitude of w_tilda to be .1 w_c
!
          sigma(i) = .91 * wbar(i) / (wc(i) + 1.E-20) + .09
!          sigma(i) = wbar(i) / (wc(i) + 1.E-20)
          sigma(i) = max(sigma(i),val1)
          if(sigma(i).gt.0.5.and.wbar(i).lt.10.)sigma(i)=0.5
          if(sigma(i).gt.0.9) then
            sigma(i)=0.9
            cnvflg(i)=.false.
          end if
        endif
        if(sigma_sum.lt.sigma(i))sigma_sum=sigma(i)
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  turn off downdraft if sigma is larger than 0.5
!c
      do i = 1, im
        cnvdflg(i) = cnvflg(i)
        if(cnvflg(i).and.sigma(i).gt.0.5) then
          cnvdflg(i) = .false.
        endif
      enddo
!c
!c  estimate the onvective overshooting as the level
!c    where the [aafac * cloud work function] becomes zero,
!c    which is the final cloud top
!c
      do i = 1, im
        if (cnvflg(i)) then
          aa2(i) = aafac * aa1(i)
        endif
      enddo
!c
      do i = 1, im
        flg(i) = cnvflg(i)
        ktcon1(i) = kmax(i) - 1
      enddo
      do k = 2, km1
        do i = 1, im
          if (flg(i)) then
            if(k.ge.ktcon(i).and.k.lt.kmax(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma           &
     &                 * to(i,k) / hvap
              aa2(i) = aa2(i) +                          &
     &                 dz1 * (g / (cp * to(i,k)))        &
     &                 * dbyo(i,k) / (1. + gamma)        &
     &                 * rfact
              if(aa2(i).lt.0.) then
                ktcon1(i) = k
                flg(i) = .false.
              endif
            endif
          endif
        enddo
      enddo
!c
!c  compute cloud moisture property, detraining cloud water
!c    and precipitation in overshooting layers
!c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.ge.ktcon(i).and.k.lt.ktcon1(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k)                               &
     &             + gamma * dbyo(i,k) / (hvap * (1. + gamma))
!cj
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5*    &
     &                     (qo(i,k)+qo(i,k-1)))/factor
!cj
              dq = eta(i,k) * (qcko(i,k) - qrch)
!c
!c  check if there is excess moisture to release latent heat
!c
              if(dq.gt.0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                if(ncloud.gt.0) then
                  dp = 1000. * del(i,k)
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                endif
                qcko(i,k) = qlk + qrch
                pwo(i,k) = etah * c0 * dz * qlk
                pwavo(i) = pwavo(i) + pwo(i,k)
              endif
            endif
          endif
        enddo
      enddo
!c
!c exchange ktcon with ktcon1
!c
      do i = 1, im
        if(cnvflg(i)) then
          kk = ktcon(i)
          ktcon(i) = ktcon1(i)
          ktcon1(i) = kk
        endif
      enddo
!c
!c  this section is ready for cloud water
!c
      if(ncloud.gt.0) then
!c
!c  compute liquid and vapor separation at cloud top
!c
      do i = 1, im
        if(cnvflg(i)) then
          k = ktcon(i) - 1
          gamma = el2orc * qeso(i,k) / (to(i,k)**2)
          qrch = qeso(i,k)                               &
     &         + gamma * dbyo(i,k) / (hvap * (1. + gamma))
          dq = qcko(i,k) - qrch
!c
!c  check if there is excess moisture to release latent heat
!c
          if(dq.gt.0.) then
            qlko_ktcon(i) = dq
            qcko(i,k) = qrch
          endif
        endif
      enddo
      endif
!c
!ccccc if(lat.eq.latd.and.lon.eq.lond.and.cnvflg(i)) then
!ccccc   print *, ' aa1(i) before dwndrft =', aa1(i)
!ccccc endif
!c
!c------- downdraft calculations
!c
!c--- compute precipitation efficiency in terms of windshear
!c
      do i = 1, im
        if(cnvdflg(i)) then
          vshear(i) = 0.
        endif
      enddo
      do k = 2, km
        do i = 1, im
          if (cnvdflg(i)) then
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              shear= sqrt((uo(i,k)-uo(i,k-1)) ** 2      &
     &                  + (vo(i,k)-vo(i,k-1)) ** 2)
              vshear(i) = vshear(i) + shear
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvdflg(i)) then
          vshear(i) = 1.e3 * vshear(i) / (zi(i,ktcon(i))-zi(i,kb(i)))
          e1=1.591-.639*vshear(i)                       &
     &       +.0953*(vshear(i)**2)-.00496*(vshear(i)**3)
          edt(i)=1.-e1
          val =         .9
          edt(i) = min(edt(i),val)
          val =         .0
          edt(i) = max(edt(i),val)
          edto(i)=edt(i)
          edtx(i)=edt(i)
        endif
      enddo
!c
!c  determine detrainment rate between 1 and kbcon
!c
      do i = 1, im
        if(cnvdflg(i)) then
          sumx(i) = 0.
        endif
      enddo
      do k = 1, km1
      do i = 1, im
        if(cnvdflg(i).and.k.ge.1.and.k.lt.kbcon(i)) then
          dz = zi(i,k+1) - zi(i,k)
          sumx(i) = sumx(i) + dz
        endif
      enddo
      enddo
      do i = 1, im
        beta = betas
        if(slimsk(i).eq.1.) beta = betal
        if(cnvdflg(i)) then
          dz  = (sumx(i)+zi(i,1))/float(kbcon(i))
          tem = 1./float(kbcon(i))
          xlamd(i) = (1.-beta**tem)/dz
        endif
      enddo
!c
!c  determine downdraft mass flux
!c
      do k = km1, 1, -1
        do i = 1, im
          if (cnvdflg(i) .and. k .le. kmax(i)-1) then
           if(k.lt.jmin(i).and.k.ge.kbcon(i)) then
              dz        = zi(i,k+1) - zi(i,k)
              ptem      = xlamdd - xlamde
              etad(i,k) = etad(i,k+1) * (1. - ptem * dz)
           else if(k.lt.kbcon(i)) then
              dz        = zi(i,k+1) - zi(i,k)
              ptem      = xlamd(i) + xlamdd - xlamde
              etad(i,k) = etad(i,k+1) * (1. - ptem * dz)
           endif
          endif
        enddo
      enddo
!c
!c--- downdraft moisture properties
!c
      do i = 1, im
        if(cnvdflg(i)) then
          jmn = jmin(i)
          hcdo(i,jmn) = heo(i,jmn)
          qcdo(i,jmn) = qo(i,jmn)
          qrcdo(i,jmn)= qeso(i,jmn)
          ucdo(i,jmn) = uo(i,jmn)
          vcdo(i,jmn) = vo(i,jmn)
          pwevo(i) = 0.
        endif
      enddo
!cj
      do k = km1, 1, -1
        do i = 1, im
          if (cnvdflg(i) .and. k.lt.jmin(i)) then
              dz = zi(i,k+1) - zi(i,k)
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = 0.5 * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = 0.5 * (xlamd(i)+xlamdd) * dz
              endif
              factor = 1. + tem - tem1
              ptem = 0.5 * tem - pgcon
              ptem1= 0.5 * tem + pgcon
              hcdo(i,k) = ((1.-tem1)*hcdo(i,k+1)+tem*0.5*        &
     &                     (heo(i,k)+heo(i,k+1)))/factor
              ucdo(i,k) = ((1.-tem1)*ucdo(i,k+1)+ptem*uo(i,k+1)  &
     &                     +ptem1*uo(i,k))/factor
              vcdo(i,k) = ((1.-tem1)*vcdo(i,k+1)+ptem*vo(i,k+1)  &
     &                     +ptem1*vo(i,k))/factor
              dbyod(i,k) = hcdo(i,k) - heso(i,k)
          endif
        enddo
      enddo
!c
      do k = km1, 1, -1
        do i = 1, im
          if (cnvdflg(i).and.k.lt.jmin(i)) then
              gamma      = el2orc * qeso(i,k) / (to(i,k)**2)
              qrcdo(i,k) = qeso(i,k)+                            &
     &                (1./hvap)*(gamma/(1.+gamma))*dbyod(i,k)
!             detad      = etad(i,k+1) - etad(i,k)
!cj
              dz = zi(i,k+1) - zi(i,k)
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = 0.5 * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = 0.5 * (xlamd(i)+xlamdd) * dz
              endif
              factor = 1. + tem - tem1
              qcdo(i,k) = ((1.-tem1)*qcdo(i,k+1)+tem*0.5*        &
     &                     (qo(i,k)+qo(i,k+1)))/factor
!cj
!             pwdo(i,k)  = etad(i,k+1) * qcdo(i,k+1) -
!    &                     etad(i,k) * qrcdo(i,k)
!             pwdo(i,k)  = pwdo(i,k) - detad *
!    &                    .5 * (qrcdo(i,k) + qrcdo(i,k+1))
!cj
              pwdo(i,k)  = etad(i,k+1) * (qcdo(i,k) - qrcdo(i,k))
              qcdo(i,k)  = qrcdo(i,k)
              pwevo(i)   = pwevo(i) + pwdo(i,k)
          endif
        enddo
      enddo
!c
!c--- final downdraft strength dependent on precip
!c--- efficiency (edt), normalized condensate (pwav), and
!c--- evaporate (pwev)
!c
      do i = 1, im
        edtmax = edtmaxl
        if(slimsk(i).eq.0.) edtmax = edtmaxs
        if(cnvdflg(i)) then
          if(pwevo(i).lt.0.) then
            edto(i) = -edto(i) * pwavo(i) / pwevo(i)
            edto(i) = min(edto(i),edtmax)
          else
            edto(i) = 0.
          endif
        endif
      enddo
!c
!c--- downdraft cloudwork functions
!c
      do k = km1, 1, -1
        do i = 1, im
          if (cnvdflg(i) .and. k .lt. jmin(i)) then
              gamma = el2orc * qeso(i,k) / to(i,k)**2
              dhh=hcdo(i,k)
              dt=to(i,k)
              dg=gamma
              dh=heso(i,k)
              dz=-1.*(zo(i,k+1)-zo(i,k))
              aa1(i)=aa1(i)+edto(i)*dz*(g/(cp*dt))*((dhh-dh)/(1.+dg))  &
     &               *(1.+delta*cp*dg*dt/hvap)
              val=0.
              aa1(i)=aa1(i)+edto(i)*                                   &
     &        dz*g*delta*max(val,(qeso(i,k)-qo(i,k)))
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvdflg(i).and.aa1(i).le.0.) then
           cnvdflg(i) = .false.
           cnvflg(i) = .false.
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  calculate environmental values of heo, qo, uo, and vo
!c
!c   updraft
      sigma_sum=0.
      do i = 1, im
        if(cnvflg(i)) then
           tem = 1. - sigma(i)
           sigi1(i) = 1. / tem
           sigi2(i) = sigma(i) / tem
           if(sigma_sum.lt.sigi1(i))sigma_sum=sigi1(i)
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              heotu(i,k)=sigi1(i)*heo(i,k)-sigi2(i)*hcko(i,k)
              qotu(i,k) =sigi1(i)*qo(i,k) -sigi2(i)*qcko(i,k)
              uotu(i,k) =sigi1(i)*uo(i,k) -sigi2(i)*ucko(i,k)
              votu(i,k) =sigi1(i)*vo(i,k) -sigi2(i)*vcko(i,k)
            endif
          endif
        enddo
      enddo
!c   downdraft
      do i = 1, im
        if(cnvdflg(i)) then
           tem = 1. - edto(i)*sigma(i)
           sigi1(i) = 1. / tem
           sigi2(i) = edto(i)*sigma(i) / tem
        endif
      enddo
      do k = 1, km1
        do i = 1, im
          if (cnvdflg(i)) then
!
!   we need to define heotd at jmin level
!
            if(k.le.jmin(i)) then
!            if(k.lt.jmin(i)) then
              heotd(i,k)=sigi1(i)*heo(i,k)-sigi2(i)*hcdo(i,k)
              qotd(i,k) =sigi1(i)*qo(i,k) -sigi2(i)*qcdo(i,k)
              uotd(i,k) =sigi1(i)*uo(i,k) -sigi2(i)*ucdo(i,k)
              votd(i,k) =sigi1(i)*vo(i,k) -sigi2(i)*vcdo(i,k)
            endif
          endif
        enddo
      enddo

!      GO TO 659
! 
! Do iteration to the cloud property using the environmental properties now
!
!c
!c  compute updraft cloud properties
!c
      do i = 1, im
        if(cnvflg(i)) then
          indx         = kb(i)
          hcko(i,indx) = heo(i,indx)
          ucko(i,indx) = uo(i,indx)
          vcko(i,indx) = vo(i,indx)
          pwavo(i)     = 0.
        endif
      enddo
!c
!c  cloud property is modified by the entrainment process
!c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.kmax(i)) then
              dz   = zi(i,k) - zi(i,k-1)
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              ptem = 0.5 * tem + pgcon
              ptem1= 0.5 * tem - pgcon
              hcko(i,k) = ((1.-tem1)*hcko(i,k-1)+tem*0.5*       &
     &                     (heotu(i,k)+heotu(i,k-1)))/factor
              ucko(i,k) = ((1.-tem1)*ucko(i,k-1)+ptem*uotu(i,k)   &
     &                     +ptem1*uotu(i,k-1))/factor
              vcko(i,k) = ((1.-tem1)*vcko(i,k-1)+ptem*votu(i,k)   &
     &                     +ptem1*votu(i,k-1))/factor
            endif
          endif
        enddo
      enddo
!c
!c  compute cloud moisture property and precipitation
!c
      do i = 1, im
        if (cnvflg(i)) then
          indx = kb(i)
          qcko(i,indx) = qo(i,indx)
          qcko(i,1) = qcko(i,indx)
          pwavo(i) = 0.
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k)                                  &
     &             + gamma * dbyo(i,k) / (hvap * (1. + gamma))
!cj
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5*       &
     &                     (qotu(i,k)+qotu(i,k-1)))/factor
!cj
              dq = eta(i,k) * (qcko(i,k) - qrch)
!c
!c
!c  check if there is excess moisture to release latent heat
!c
              if(k.ge.kbcon(i).and.dq.gt.0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                if(ncloud.gt.0.and.k.gt.jmin(i)) then
                  dp = 1000. * del(i,k)
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                endif
                qcko(i,k) = qlk + qrch
                pwo(i,k) = etah * c0 * dz * qlk
                pwavo(i) = pwavo(i) + pwo(i,k)
              endif
            endif
          endif
        enddo
      enddo
!c
!c  this section is ready for cloud water
!c
      if(ncloud.gt.0) then
!c
!c  compute liquid and vapor separation at cloud top
!c
      do i = 1, im
        if(cnvflg(i)) then
          k = ktcon(i) - 1
          gamma = el2orc * qeso(i,k) / (to(i,k)**2)
          qrch = qeso(i,k)                               &
     &         + gamma * dbyo(i,k) / (hvap * (1. + gamma))
          dq = qcko(i,k) - qrch
!c
!c  check if there is excess moisture to release latent heat
!c
          if(dq.gt.0.) then
            qlko_ktcon(i) = dq
            qcko(i,k) = qrch
          endif
        endif
      enddo
      endif
!c
!c--- downdraft moisture properties
!c
      do i = 1, im
        if(cnvdflg(i)) then
          jmn = jmin(i)
          hcdo(i,jmn) = heo(i,jmn)
          qcdo(i,jmn) = qo(i,jmn)
          qrcdo(i,jmn)= qeso(i,jmn)
          ucdo(i,jmn) = uo(i,jmn)
          vcdo(i,jmn) = vo(i,jmn)
          pwevo(i) = 0.
        endif
      enddo
!cj
      do k = km1, 1, -1
        do i = 1, im
          if (cnvdflg(i) .and. k.lt.jmin(i)) then
              dz = zi(i,k+1) - zi(i,k)
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = 0.5 * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = 0.5 * (xlamd(i)+xlamdd) * dz
              endif
              factor = 1. + tem - tem1
              ptem = 0.5 * tem - pgcon
              ptem1= 0.5 * tem + pgcon
              hcdo(i,k) = ((1.-tem1)*hcdo(i,k+1)+tem*0.5*        &
     &                     (heotd(i,k)+heotd(i,k+1)))/factor
              ucdo(i,k) = ((1.-tem1)*ucdo(i,k+1)+ptem*uotd(i,k+1)  &
     &                     +ptem1*uotd(i,k))/factor
              vcdo(i,k) = ((1.-tem1)*vcdo(i,k+1)+ptem*votd(i,k+1)  &
     &                     +ptem1*votd(i,k))/factor
          endif
        enddo
      enddo
!c
      do k = km1, 1, -1
        do i = 1, im
          if (cnvdflg(i).and.k.lt.jmin(i)) then
              gamma      = el2orc * qeso(i,k) / (to(i,k)**2)
              qrcdo(i,k) = qeso(i,k)+                            &
     &                (1./hvap)*(gamma/(1.+gamma))*dbyod(i,k)
!             detad      = etad(i,k+1) - etad(i,k)
!cj
              dz = zi(i,k+1) - zi(i,k)
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = 0.5 * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = 0.5 * (xlamd(i)+xlamdd) * dz
              endif
              factor = 1. + tem - tem1
              qcdo(i,k) = ((1.-tem1)*qcdo(i,k+1)+tem*0.5*        &
     &                     (qotd(i,k)+qotd(i,k+1)))/factor
              pwdo(i,k)  = etad(i,k+1) * (qcdo(i,k) - qrcdo(i,k))
              qcdo(i,k)  = qrcdo(i,k)
              pwevo(i)   = pwevo(i) + pwdo(i,k)
          endif
        enddo
      enddo

 659  continue

!c
!c--- what would the change be, that a cloud with unit mass
!c--- will do to the environment?
!c
      do k = 1, km
        do i = 1, im
          if(cnvflg(i) .and. k .le. kmax(i)) then
            dellah(i,k) = 0.
            dellaq(i,k) = 0.
            dellau(i,k) = 0.
            dellav(i,k) = 0.
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvdflg(i)) then
          dp = 1000. * del(i,1)
          dellah(i,1) = edto(i) * etad(i,1) * (hcdo(i,1)          &
     &                   - heotd(i,1)) * g / dp
!
          dellaq(i,1) = edto(i) * etad(i,1) * (qcdo(i,1)          &
     &                   - qotd(i,1)) * g / dp
!
          dellau(i,1) = edto(i) * etad(i,1) * (ucdo(i,1)          &
     &                   - uotd(i,1)) * g / dp
!
          dellav(i,1) = edto(i) * etad(i,1) * (vcdo(i,1)          &
     &                   - votd(i,1)) * g / dp
        endif
      enddo
!c
!c--- changed due to subsidence and entrainment
!c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i).and.k.lt.ktcon(i)) then
              aup = 1.
              if(k.le.kb(i)) aup = 0.
              adw = 1.
              if(k.gt.jmin(i).or..not.cnvdflg(i)) adw = 0.
              dp = 1000. * del(i,k)
              dz = zi(i,k) - zi(i,k-1)
!c
              dv1h = heotu(i,k)
              dv2h = .5 * (heotu(i,k) + heotu(i,k-1))
              dv3h = heotu(i,k-1)
              dv1q = qotu(i,k)
              dv2q = .5 * (qotu(i,k) + qotu(i,k-1))
              dv3q = qotu(i,k-1)
              dv1u = uotu(i,k)
              dv2u = .5 * (uotu(i,k) + uotu(i,k-1))
              dv3u = uotu(i,k-1)
              dv1v = votu(i,k)
              dv2v = .5 * (votu(i,k) + votu(i,k-1))
              dv3v = votu(i,k-1)
!c
              dv1hd = heotd(i,k)
              dv2hd = .5 * (heotd(i,k) + heotd(i,k-1))
              dv3hd = heotd(i,k-1)
              dv1qd = qotd(i,k)
              dv2qd = .5 * (qotd(i,k) + qotd(i,k-1))
              dv3qd = qotd(i,k-1)
              dv1ud = uotd(i,k)
              dv2ud = .5 * (uotd(i,k) + uotd(i,k-1))
              dv3ud = uotd(i,k-1)
              dv1vd = votd(i,k)
              dv2vd = .5 * (votd(i,k) + votd(i,k-1))
              dv3vd = votd(i,k-1)
!c
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1))
              tem1 = xlamud(i)
!c
              if(k.le.kbcon(i)) then
                ptem  = xlamde
                ptem1 = xlamd(i)+xlamdd
              else
                ptem  = xlamde
                ptem1 = xlamdd
              endif
!cj
              dellah(i,k) = dellah(i,k) +                               &
     &     (aup*(eta(i,k)*heotu(i,k)-eta(i,k-1)*heotu(i,k-1))           &
     &    - adw*edto(i)*(etad(i,k)*heotd(i,k)-etad(i,k-1)*heotd(i,k-1)) &
!     &    - (aup*tem*eta(i,k-1)+adw*edto(i)*ptem*etad(i,k))*dv2h*dz     &
     &    - (aup*tem*eta(i,k-1))*dv2h*dz-(adw*edto(i)*ptem*etad(i,k))*dv2hd*dz     &
     &    +  aup*tem1*eta(i,k-1)*.5*(hcko(i,k)+hcko(i,k-1))*dz          &
     &    +  adw*edto(i)*ptem1*etad(i,k)*.5*(hcdo(i,k)+hcdo(i,k-1))*dz  &
     &         ) *g/dp

!
!c delhx = -(g/dp) * (rho * wbar) * del(htilda-hbar) 
!c rho * g * wbar is replaced by -omega_bar which is doto in Pa/sec
!c dp is in Pa
!
              delhx(i,k) = ( aup*(doto(i,k)*(heotu(i,k)-heo(i,k))           &
     &                          - doto(i,k-1)*(heotu(i,k-1)-heo(i,k-1)))     &
     &                     ) / dp
!cj
              dellaq(i,k) = dellaq(i,k) +                               &
     &     (aup*(eta(i,k)*qotu(i,k)-eta(i,k-1)*qotu(i,k-1))             &
     &    - adw*edto(i)*(etad(i,k)*qotd(i,k)-etad(i,k-1)*qotd(i,k-1))   &
     &    - (aup*tem*eta(i,k-1))*dv2q*dz-(adw*edto(i)*ptem*etad(i,k))*dv2qd*dz     &
     &    +  aup*tem1*eta(i,k-1)*.5*(qcko(i,k)+qcko(i,k-1))*dz          &
     &    +  adw*edto(i)*ptem1*etad(i,k)*.5*(qrcdo(i,k)+qrcdo(i,k-1))*dz  &
     &         ) *g/dp
!

              delqx(i,k) = ( aup*(doto(i,k)*(qotu(i,k)-qo(i,k))            &
     &                          - doto(i,k-1)*(qotu(i,k-1)-qo(i,k-1)))       &
     &                     ) / dp
!cj
              dellau(i,k) = dellau(i,k) +                               &
     &     (aup*(eta(i,k)*uotu(i,k)-eta(i,k-1)*uotu(i,k-1))             &
     &    - adw*edto(i)*(etad(i,k)*uotd(i,k)-etad(i,k-1)*uotd(i,k-1))   &
     &    - (aup*tem*eta(i,k-1))*dv2u*dz-(adw*edto(i)*ptem*etad(i,k))*dv2ud*dz     &
     &    +  aup*tem1*eta(i,k-1)*.5*(ucko(i,k)+ucko(i,k-1))*dz          &
     &    +  adw*edto(i)*ptem1*etad(i,k)*.5*(ucdo(i,k)+ucdo(i,k-1))*dz  &
     &    -  pgcon*(aup*eta(i,k-1)*(dv1u-dv3u)-adw*edto(i)*etad(i,k)*(dv1ud-dv3ud))   &
     &         ) *g/dp
!
              delux(i,k) = ( aup*(doto(i,k)*(uotu(i,k)-uo(i,k))            &
     &                          - doto(i,k-1)*(uotu(i,k-1)-uo(i,k-1)))       &
     &                     ) / dp
!cj
              dellav(i,k) = dellav(i,k) +                               &
     &     (aup*(eta(i,k)*votu(i,k)-eta(i,k-1)*votu(i,k-1))             &
     &    - adw*edto(i)*(etad(i,k)*votd(i,k)-etad(i,k-1)*votd(i,k-1))   &
     &    - (aup*tem*eta(i,k-1))*dv2v*dz-(adw*edto(i)*ptem*etad(i,k))*dv2vd*dz     &
     &    +  aup*tem1*eta(i,k-1)*.5*(vcko(i,k)+vcko(i,k-1))*dz          &
     &    +  adw*edto(i)*ptem1*etad(i,k)*.5*(vcdo(i,k)+vcdo(i,k-1))*dz  &
     &    -  pgcon*(aup*eta(i,k-1)*(dv1v-dv3v)-adw*edto(i)*etad(i,k)*(dv1vd-dv3vd))   &
     &         ) *g/dp
!

              delvx(i,k) = ( aup*(doto(i,k)*(votu(i,k)-vo(i,k))            &
     &                          - doto(i,k-1)*(votu(i,k-1)-vo(i,k-1)))       &
     &                     ) / dp
!          if(abs(delvx(i,k)).gt.1.0)print*,'qingfu test999=',        &
!              i,k,delvx(i,k),aup,adw,doto(i,k),(votu(i,k)-vo(i,k))      &
!              ,(votu(i,k-1)-vo(i,k-1)),dp

!cj
          endif
        enddo
      enddo
!c
!c------- cloud top
!c
      do i = 1, im
        if(cnvflg(i)) then
          indx = ktcon(i)
          dp = 1000. * del(i,indx)
          dv1h = heo(i,indx-1)
          dellah(i,indx) = eta(i,indx-1) *                              &
     &              (hcko(i,indx-1) - heotu(i,indx-1)) * g / dp
!          delhx(i,indx) = doto(i,indx-1)*dv1h / dp
          delhx(i,indx) = doto(i,indx-1)*(dv1h-heotu(i,indx-1)) / dp
!
          dv1q = qo(i,indx-1)
          dellaq(i,indx) = eta(i,indx-1) *                              &
     &              (qcko(i,indx-1) - qotu(i,indx-1)) * g / dp
!          delqx(i,indx) = doto(i,indx-1)*dv1q / dp
          delqx(i,indx) = doto(i,indx-1)*(dv1q-qotu(i,indx-1)) / dp
!
          dv1u = uo(i,indx-1)
          dellau(i,indx) = eta(i,indx-1) *                              &
     &              (ucko(i,indx-1) - uotu(i,indx-1)) * g / dp
!          delux(i,indx) = doto(i,indx-1)*dv1u / dp
          delux(i,indx) = doto(i,indx-1)*(dv1u-uotu(i,indx-1)) / dp
!
          dv1v = vo(i,indx-1)
          dellav(i,indx) = eta(i,indx-1) *                              &
     &              (vcko(i,indx-1) - votu(i,indx-1)) * g / dp
!          delvx(i,indx) = doto(i,indx-1)*dv1v / dp
          delvx(i,indx) = doto(i,indx-1)*(dv1v-votu(i,indx-1)) / dp
!          if(abs(delvx(i,indx)).gt.5.0)print*,'qingfu test888=',      &
!              i,indx,delvx(i,indx),doto(i,indx-1),dv1v,votu(i,indx-1),dp
!c
!c  cloud water
!c
          dellal(i,indx) = eta(i,indx-1) *                              &
     &                     qlko_ktcon(i) * g / dp
        endif
      enddo
!c
!c------- final changed variable per unit mass flux
!c
      do k = 1, km
        do i = 1, im
          if (cnvflg(i).and.k .le. kmax(i)) then
            if(k.gt.ktcon(i)) then
              qo(i,k) = q1(i,k)
              to(i,k) = t1(i,k)
            endif
            if(k.le.ktcon(i)) then
!
!c   We need to scale the w-bar contribution (delhx and delqx) by rho * wc
!c   po is in mb but rho (tem) is now in standard unit, wc is in m/sec
!c   tem1 is wbar in m/sec, doto is in pa/sec
!
              tem = po(i,k) / (rd * to(i,k) * (1. + delta * qo(i,k))) * 100.
              tem1 = -doto(i,k) / (tem * g)
              tem1 = tem1 / (wc(i) + 1.E-20)
              tem1 = max(tem1,real(0.,kind=kind_phys))
              tem1 = min(tem1,real(1.,kind=kind_phys))
!              delqz = dellaq(i,k)* (1.-tem1) *mbdt / (1. - sigma)
!     &              + delqx(i,k)*ardt / ((1. - sigma(i)) * tem * wc(i))
              delqz = dellaq(i,k) * mbdt
              qo(i,k) = q1(i,k) + delqz
!              delhz = dellah(i,k)* (1.-tem1) *mbdt / (1. - sigma)
!     &              + delhx(i,k)*ardt / ((1. - sigma(i)) * tem * wc(i))
              delhz = dellah(i,k) * mbdt
              dellat = (delhz - hvap * delqz) / cp
              to(i,k) = t1(i,k) + dellat
              val   =           1.e-10
              qo(i,k) = max(qo(i,k), val  )
            endif
          endif
        enddo
      enddo
!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c
!c--- the above changed environment is now used to calulate the
!c--- effect the arbitrary cloud (with unit mass flux)
!c--- would have on the stability,
!c--- which then is used to calculate the real mass flux,
!c--- necessary to keep this change in balance with the large-scale
!c--- destabilization.
!c
!c--- environmental conditions again, first heights
!c
      do k = 1, km
        do i = 1, im
          if(cnvflg(i) .and. k .le. kmax(i)) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k)+epsm1*qeso(i,k))
            val       =             1.e-8
            qeso(i,k) = max(qeso(i,k), val )
!           tvo(i,k)  = to(i,k) + delta * to(i,k) * qo(i,k)
          endif
        enddo
      enddo
!c
!c--- moist static energy
!c
      do k = 1, km1
        do i = 1, im
          if(cnvflg(i) .and. k .le. kmax(i)-1) then
            dz = .5 * (zo(i,k+1) - zo(i,k))
            dp = .5 * (pfld(i,k+1) - pfld(i,k))
            es = 0.01 * fpvs(to(i,k+1))      ! fpvs is in pa
            pprime = pfld(i,k+1) + epsm1 * es
            qs = eps * es / pprime
            dqsdp = - qs / pprime
            desdt = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
            dqsdt = qs * pfld(i,k+1) * desdt / (es * pprime)
            gamma = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
            dt = (g * dz + hvap * dqsdp * dp) / (cp * (1. + gamma))
            dq = dqsdt * dt + dqsdp * dp
            to(i,k) = to(i,k+1) + dt
            qo(i,k) = qo(i,k+1) + dq
            po(i,k) = .5 * (pfld(i,k) + pfld(i,k+1))
          endif
        enddo
      enddo
      do k = 1, km1
        do i = 1, im
          if(cnvflg(i) .and. k .le. kmax(i)-1) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (po(i,k) + epsm1 * qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
            heo(i,k)   = .5 * g * (zo(i,k) + zo(i,k+1)) +           &
     &                    cp * to(i,k) + hvap * qo(i,k)
            heso(i,k) = .5 * g * (zo(i,k) + zo(i,k+1)) +            &
     &                  cp * to(i,k) + hvap * qeso(i,k)
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          k = kmax(i)
          heo(i,k) = g * zo(i,k) + cp * to(i,k) + hvap * qo(i,k)
          heso(i,k) = g * zo(i,k) + cp * to(i,k) + hvap * qeso(i,k)
!c         heo(i,k) = min(heo(i,k),heso(i,k))
        endif
      enddo
!c
!c**************************** static control
!c
!c------- moisture and cloud work functions
!c
      do i = 1, im
        if(cnvflg(i)) then
          xaa0(i) = 0.
          xpwav(i) = 0.
        endif
      enddo
!c
      do i = 1, im
        if(cnvflg(i)) then
          indx = kb(i)
          hcko(i,indx) = heo(i,indx)
          qcko(i,indx) = qo(i,indx)
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              dz = zi(i,k) - zi(i,k-1)
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              hcko(i,k) = ((1.-tem1)*hcko(i,k-1)+tem*0.5*          &
     &                     (heo(i,k)+heo(i,k-1)))/factor
            endif
          endif
        enddo
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              dz = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              xdby = hcko(i,k) - heso(i,k)
              xqrch = qeso(i,k)                                    &
     &              + gamma * xdby / (hvap * (1. + gamma))
!cj
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5*          &
     &                     (qo(i,k)+qo(i,k-1)))/factor
!cj
              dq = eta(i,k) * (qcko(i,k) - xqrch)
!c
              if(k.ge.kbcon(i).and.dq.gt.0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                if(ncloud.gt.0.and.k.gt.jmin(i)) then
                  qlk = dq / (eta(i,k) + etah * (c0 + c1) * dz)
                else
                  qlk = dq / (eta(i,k) + etah * c0 * dz)
                endif
                if(k.lt.ktcon1(i)) then
                  xaa0(i) = xaa0(i) - dz * g * qlk
                endif
                qcko(i,k) = qlk + xqrch
                xpw = etah * c0 * dz * qlk
                xpwav(i) = xpwav(i) + xpw
              endif
            endif
            if(k.ge.kbcon(i).and.k.lt.ktcon1(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma                    &
     &                 * to(i,k) / hvap
              xaa0(i) = xaa0(i)                                   &
     &                + dz1 * (g / (cp * to(i,k)))                &
     &                * xdby / (1. + gamma)                       &
     &                * rfact
              val=0.
              xaa0(i)=xaa0(i)+                                    &
     &                 dz1 * g * delta *                          &
     &                 max(val,(qeso(i,k) - qo(i,k)))
            endif
          endif
        enddo
      enddo
!c
!c------- downdraft calculations
!c
!c--- downdraft moisture properties
!c
      do i = 1, im
        if(cnvdflg(i)) then
          jmn = jmin(i)
          hcdo(i,jmn) = heo(i,jmn)
          qcdo(i,jmn) = qo(i,jmn)
          qrcd(i,jmn) = qeso(i,jmn)
          xpwev(i) = 0.
        endif
      enddo
!cj
      do k = km1, 1, -1
        do i = 1, im
          if (cnvdflg(i) .and. k.lt.jmin(i)) then
              dz = zi(i,k+1) - zi(i,k)
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = 0.5 * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = 0.5 * (xlamd(i)+xlamdd) * dz
              endif
              factor = 1. + tem - tem1
              hcdo(i,k) = ((1.-tem1)*hcdo(i,k+1)+tem*0.5*        &
     &                     (heo(i,k)+heo(i,k+1)))/factor
          endif
        enddo
      enddo
!cj
      do k = km1, 1, -1
        do i = 1, im
          if (cnvdflg(i) .and. k .lt. jmin(i)) then
              dq = qeso(i,k)
              dt = to(i,k)
              gamma    = el2orc * dq / dt**2
              dh       = hcdo(i,k) - heso(i,k)
              qrcd(i,k)=dq+(1./hvap)*(gamma/(1.+gamma))*dh
!             detad    = etad(i,k+1) - etad(i,k)
!cj
              dz = zi(i,k+1) - zi(i,k)
              if(k.ge.kbcon(i)) then
                 tem  = xlamde * dz
                 tem1 = 0.5 * xlamdd * dz
              else
                 tem  = xlamde * dz
                 tem1 = 0.5 * (xlamd(i)+xlamdd) * dz
              endif
              factor = 1. + tem - tem1
              qcdo(i,k) = ((1.-tem1)*qcdo(i,k+1)+tem*0.5*        &
     &                     (qo(i,k)+qo(i,k+1)))/factor
!cj
!             xpwd     = etad(i,k+1) * qcdo(i,k+1) -
!    &                   etad(i,k) * qrcd(i,k)
!             xpwd     = xpwd - detad *
!    &                 .5 * (qrcd(i,k) + qrcd(i,k+1))
!cj
              xpwd     = etad(i,k+1) * (qcdo(i,k) - qrcd(i,k))
              qcdo(i,k)= qrcd(i,k)
              xpwev(i) = xpwev(i) + xpwd
          endif
        enddo
      enddo
!c
      do i = 1, im
        edtmax = edtmaxl
        if(slimsk(i).eq.0.) edtmax = edtmaxs
        if(cnvdflg(i)) then
          if(xpwev(i).ge.0.) then
            edtx(i) = 0.
          else
            edtx(i) = -edtx(i) * xpwav(i) / xpwev(i)
            edtx(i) = min(edtx(i),edtmax)
          endif
        endif
      enddo
!c
!c
!c--- downdraft cloudwork functions
!c
!c
      do k = km1, 1, -1
        do i = 1, im
          if (cnvdflg(i) .and. k.lt.jmin(i)) then
              gamma = el2orc * qeso(i,k) / to(i,k)**2
              dhh=hcdo(i,k)
              dt= to(i,k)
              dg= gamma
              dh= heso(i,k)
              dz=-1.*(zo(i,k+1)-zo(i,k))
              xaa0(i)=xaa0(i)+edtx(i)*dz*(g/(cp*dt))*((dhh-dh)/(1.+dg))  &
     &                *(1.+delta*cp*dg*dt/hvap)
              val=0.
              xaa0(i)=xaa0(i)+edtx(i)*                             &
     &        dz*g*delta*max(val,(qeso(i,k)-qo(i,k)))
          endif
        enddo
      enddo
!c
!c  calculate critical cloud work function
!c
      do i = 1, im
        if(cnvflg(i)) then
          if(pfld(i,ktcon(i)).lt.pcrit(15))then
            acrt(i)=acrit(15)*(975.-pfld(i,ktcon(i)))              &
     &              /(975.-pcrit(15))
          else if(pfld(i,ktcon(i)).gt.pcrit(1))then
            acrt(i)=acrit(1)
          else
            k =  int((850. - pfld(i,ktcon(i)))/50.) + 2
            k = min(k,15)
            k = max(k,2)
            acrt(i)=acrit(k)+(acrit(k-1)-acrit(k))*                &
     &           (pfld(i,ktcon(i))-pcrit(k))/(pcrit(k-1)-pcrit(k))
          endif
        endif
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          if(slimsk(i).eq.1.) then
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          endif
!c
!c  modify critical cloud workfunction by cloud base vertical velocity
!c
          if(pdot(i).le.w4) then
            acrtfct(i) = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i).ge.-w4) then
            acrtfct(i) = - (pdot(i) + w4) / (w4 - w3)
          else
            acrtfct(i) = 0.
          endif
          val1    =             -1.
          acrtfct(i) = max(acrtfct(i),val1)
          val2    =             1.
          acrtfct(i) = min(acrtfct(i),val2)
          acrtfct(i) = 1. - acrtfct(i)
!c
!c  modify acrtfct(i) by colume mean rh if rhbar(i) is greater than 80 percent
!c
!c         if(rhbar(i).ge..8) then
!c           acrtfct(i) = acrtfct(i) * (.9 - min(rhbar(i),.9)) * 10.
!c         endif
!c
!c  modify adjustment time scale by cloud base vertical velocity
!c
          val1=0.0
          dtconv(i) = dt2 + max((1800. - dt2),val1) *             &
     &                (pdot(i) - w2) / (w1 - w2)
!c         dtconv(i) = max(dtconv(i), dt2)
!c         dtconv(i) = 1800. * (pdot(i) - w2) / (w1 - w2)
          dtconv(i) = max(dtconv(i),dtmin)
          dtconv(i) = min(dtconv(i),dtmax)
!c
        endif
      enddo
!c
!c--- large scale forcing
!c
      do i= 1, im
        if(cnvflg(i)) then
          fld(i)=(aa1(i)-acrt(i)* acrtfct(i))/dtconv(i)
          if(fld(i).le.0.) cnvflg(i) = .false.
        endif
        if(cnvflg(i)) then
!c         xaa0(i) = max(xaa0(i),0.)
          xk(i) = (xaa0(i) - aa1(i)) / mbdt
          if(xk(i).ge.0.) cnvflg(i) = .false.
        endif
!c
!c--- kernel, cloud base mass flux
!c
        if(cnvflg(i)) then
          xmb(i) = -fld(i) / xk(i)
          xmb(i) = min(xmb(i),xmbmax(i))
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!c
!c  restore to,qo,uo,vo to t1,q1,u1,v1 in case convection stops
!c
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)) then
            to(i,k) = t1(i,k)
            qo(i,k) = q1(i,k)
            uo(i,k) = u1(i,k)
            vo(i,k) = v1(i,k)
            qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val     =             1.e-8
            qeso(i,k) = max(qeso(i,k), val )
          endif
        enddo
      enddo
!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!c
!c--- feedback: simply the changes from the cloud with unit mass flux
!c---           multiplied by  the mass flux necessary to keep the
!c---           equilibrium with the larger-scale.
!c
      do i = 1, im
        delhbar(i) = 0.
        delqbar(i) = 0.
        deltbar(i) = 0.
        delubar(i) = 0.
        delvbar(i) = 0.
        qcond(i) = 0.
      enddo
!c
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)) then
            if(k.le.ktcon(i)) then
              delhz = dellah(i,k)*xmb(i) + delhx(i,k)
              delqz = dellaq(i,k)*xmb(i) + delqx(i,k)
              deluz = dellau(i,k)*xmb(i) + delux(i,k)
              delvz = dellav(i,k)*xmb(i) + delvx(i,k)
              dellat = (delhz - hvap * delqz) / cp
              t1(i,k) = t1(i,k) + dellat * dt2
              q1(i,k) = q1(i,k) + delqz * dt2
              tem = 1./rcs(i)
              u1(i,k) = u1(i,k) + deluz * dt2 * tem
              v1(i,k) = v1(i,k) + delvz * dt2 * tem
              dp = 1000. * del(i,k)
              delhbar(i) = delhbar(i) + delhz * dp / g
              delqbar(i) = delqbar(i) + delqz * dp / g
              deltbar(i) = deltbar(i) + dellat * dp / g
              delubar(i) = delubar(i) + deluz * dp / g
              delvbar(i) = delvbar(i) + delvz * dp / g
            endif
          endif
        enddo
      enddo
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)) then
            if(k.le.ktcon(i)) then
              qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
              qeso(i,k) = eps * qeso(i,k)/(pfld(i,k) + epsm1*qeso(i,k))
              val     =             1.e-8
              qeso(i,k) = max(qeso(i,k), val )
            endif
          endif
        enddo
      enddo
!c
      do i = 1, im
        rntot(i) = 0.
        delqev(i) = 0.
        delq2(i) = 0.
        flg(i) = cnvflg(i)
      enddo
      do k = km, 1, -1
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)) then
            if(k.lt.ktcon(i)) then
              aup = 1.
              if(k.le.kb(i)) aup = 0.
              adw = 1.
              if(k.ge.jmin(i).or..not.cnvdflg(i)) adw = 0.
              rain =  aup * pwo(i,k) + adw * edto(i) * pwdo(i,k)
              rntot(i) = rntot(i) + rain * xmb(i) * .001 * dt2
            endif
          endif
        enddo
      enddo
      do k = km, 1, -1
        do i = 1, im
          if (k .le. kmax(i)) then
            deltv(i) = 0.
            delq(i) = 0.
            qevap(i) = 0.
            if(cnvflg(i).and.k.lt.ktcon(i)) then
              aup = 1.
              if(k.le.kb(i)) aup = 0.
              adw = 1.
              if(k.ge.jmin(i).or..not.cnvdflg(i)) adw = 0.
              rain =  aup * pwo(i,k) + adw * edto(i) * pwdo(i,k)
              rn(i) = rn(i) + rain * xmb(i) * .001 * dt2
            endif
            if(flg(i).and.k.lt.ktcon(i)) then
              evef = edt(i) * evfact
              if(slimsk(i).eq.1.) evef=edt(i) * evfactl
!             if(slimsk(i).eq.1.) evef=.07
!c             if(slimsk(i).ne.1.) evef = 0.
              qcond(i) = evef * (q1(i,k) - qeso(i,k))                &
     &                 / (1. + el2orc * qeso(i,k) / t1(i,k)**2)
              dp = 1000. * del(i,k)
              if(rn(i).gt.0..and.qcond(i).lt.0.) then
                qevap(i) = -qcond(i) * (1.-exp(-.32*sqrt(dt2*rn(i))))
                qevap(i) = min(qevap(i), rn(i)*1000.*g/dp)
                delq2(i) = delqev(i) + .001 * qevap(i) * dp / g
              endif
              if(rn(i).gt.0..and.qcond(i).lt.0..and.                 &
     &           delq2(i).gt.rntot(i)) then
                qevap(i) = 1000.* g * (rntot(i) - delqev(i)) / dp
                flg(i) = .false.
              endif
              if(rn(i).gt.0..and.qevap(i).gt.0.) then
                q1(i,k) = q1(i,k) + qevap(i)
                t1(i,k) = t1(i,k) - elocp * qevap(i)
                rn(i) = rn(i) - .001 * qevap(i) * dp / g
                deltv(i) = - elocp*qevap(i)/dt2
                delq(i) =  + qevap(i)/dt2
                delqev(i) = delqev(i) + .001*dp*qevap(i)/g
              endif
              dellaq(i,k) = dellaq(i,k) + delq(i) / xmb(i)
              delqbar(i) = delqbar(i) + delq(i)*dp/g
              deltbar(i) = deltbar(i) + deltv(i)*dp/g
            endif
          endif
        enddo
      enddo
!cj
!     do i = 1, im
!     if(me.eq.31.and.cnvflg(i)) then
!     if(cnvflg(i)) then
!       print *, ' deep delhbar, delqbar, deltbar = ',
!    &             delhbar(i),hvap*delqbar(i),cp*deltbar(i)
!       print *, ' deep delubar, delvbar = ',delubar(i),delvbar(i)
!       print *, ' precip =', hvap*rn(i)*1000./dt2
!       print*,'pdif= ',pfld(i,kbcon(i))-pfld(i,ktcon(i))
!     endif
!     enddo
!c
!c  precipitation rate converted to actual precip
!c  in unit of m instead of kg
!c
      do i = 1, im
        if(cnvflg(i)) then
!c
!c  in the event of upper level rain evaporation and lower level downdraft
!c    moistening, rn can become negative, in this case, we back out of the
!c    heating and the moistening
!c
          if(rn(i).lt.0..and..not.flg(i)) rn(i) = 0.
          if(rn(i).le.0.) then
            rn(i) = 0.
          else
            ktop(i) = ktcon(i)
            kbot(i) = kbcon(i)
            kcnv(i) = 1
            cldwrk(i) = aa1(i)
          endif
        endif
      enddo
!c
!c  cloud water
!c
      if (ncloud.gt.0) then
!
      val1=0.0
      val2=1.0
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. rn(i).gt.0.) then
            if (k.gt.kb(i).and.k.le.ktcon(i)) then
              tem  = dellal(i,k) * xmb(i) * dt2
              tem1 = max(val1, min(val2, (tcr-t1(i,k))*tcrf))
              if (ql(i,k,2) .gt. -999.0) then
                ql(i,k,1) = ql(i,k,1) + tem * tem1            ! ice
                ql(i,k,2) = ql(i,k,2) + tem *(1.0-tem1)       ! water
              else
                ql(i,k,1) = ql(i,k,1) + tem
              endif
            endif
          endif
        enddo
      enddo
!
      endif
!c
      do k = 1, km
        do i = 1, im
          if(cnvflg(i).and.rn(i).le.0.) then
            if (k .le. kmax(i)) then
              t1(i,k) = to(i,k)
              q1(i,k) = qo(i,k)
              u1(i,k) = uo(i,k)
              v1(i,k) = vo(i,k)
            endif
          endif
        enddo
      enddo
!
! hchuang code change
!
      do k = 1, km
        do i = 1, im
          if(cnvflg(i).and.rn(i).gt.0.) then
            if(k.ge.kb(i) .and. k.lt.ktop(i)) then
              ud_mf(i,k) = eta(i,k) * xmb(i) * dt2
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i).and.rn(i).gt.0.) then
           k = ktop(i)-1
           dt_mf(i,k) = ud_mf(i,k)
        endif
      enddo
      do k = 1, km
        do i = 1, im
          if(cnvflg(i).and.rn(i).gt.0.) then
            if(k.ge.1 .and. k.le.jmin(i)) then
              dd_mf(i,k) = edto(i) * etad(i,k) * xmb(i) * dt2
            endif
          endif
        enddo
      enddo
!!

      sigma_sum=0.
      do I=1,im
         sigma_sum=sigma_sum+abs(sigma(I))
      end do
!      if(sigma_sum.gt.0.1)then
!        print*,'qliu test sigma_c='
!        write(*,333)sigma
!      end if
!333   format(1x,'inside sascnvn_h sigma_c=',9F10.3)

      return
      end subroutine sascnvn_meso

      END MODULE MODULE_CU_SASHUR
!
!-----------------------------------------------------------------------
