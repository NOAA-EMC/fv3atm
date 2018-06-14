!----------------------------------------------------------------------
!
      MODULE MODULE_BL_GFSPBL
!
!-----------------------------------------------------------------------
!
!***  THE GFS PBL SCHEME
!
!     2010-0910    Created by Weiguo Wang to use GFS PBL 
!     2014-0703    Weiguo Wang, use F_QI, in some MP, QI is not defined/used.
!     2014-11-21   Brad Ferrier, removed F_QR, F_QS, F_QG, etc. dependencies
!-----------------------------------------------------------------------
!
      USE MODULE_KINDS
!
       USE MODULE_CONSTANTS,ONLY : cp99 => CP,ELWV,                &
     &                            g99 =>G, rd99 => r_d,            &
     &                            ep199 => ep_1, CAPPA, PQ0,ELIV,  &
                                  AA2 =>A2, AA3=>A3, AA4=>A4, P608
       use machine     , only : kind_phys

!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: GFSPBL
!
!-----------------------------------------------------------------------
!
      REAL,PARAMETER :: VKARMAN=0.4
      REAL,PARAMETER :: XLV=ELWV, RLIVWV=ELIV/XLV
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
! REFERENCES:  moved from GFS 
!
! ABSTRACT:
!
!-----------------------------------------------------------------------
        SUBROUTINE GFSPBL(DT,NPHS,DP,AIRDEN                              &
     &                    ,RIB                                            &
     &                    ,PHMID,PHINT,T,ZINT                             &
     &                    ,Q,QC,QI                                        &
     &                    ,F_QC,F_QI                                      &
     &                    ,U,V                                            &
     &                    ,USTAR                                          &
     &                    ,SHEAT, LHEAT                                   &
                          ,XLAND                                          &
     &                    ,AKHS,AKMS                                      &
                          ,THZ0,QZ0                                       &
                          ,QSFC                                           &
                          ,TSK,SNOW,SICE,CHKLOWQ                          &
                          ,factrs,rswtt,rlwtt                             &
     &                    ,PBLH,PBLK                                      &     !! out below
     &                    ,MIXHT                                          &
     &                    ,RUBLTEN                                        &
     &                    ,RVBLTEN                                        &
     &                    ,RTHBLTEN                                       &
     &                    ,RQBLTEN                                        &
     &                    ,RQCBLTEN                                       &
     &                    ,RQIBLTEN                                       &
     &                   ,IMS,IME,JMS,JME,LM)

!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
       REAL,PARAMETER :: A2S=17.2693882,A3S=273.16,A4S=35.86
       REAL,PARAMETER :: SEAFC=0.98,PQ0SEA=PQ0*SEAFC

      INTEGER,INTENT(IN) :: IMS,IME,JMS,JME,LM
!
      INTEGER,INTENT(IN) :: NPHS
!
!
      REAL,INTENT(IN) :: DT
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: SHEAT, LHEAT,RIB,USTAR ,XLAND
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: AKHS,AKMS,factrs
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: SICE, SNOW, TSK ,CHKLOWQ
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: QZ0, THZ0
!
      REAL,DIMENSION(IMS:IME,JMS:JME,LM),INTENT(IN)::  RSWTT, RLWTT

      REAL,DIMENSION(IMS:IME,JMS:JME,LM),INTENT(IN) :: DP,PHMID,AIRDEN
      REAL,DIMENSION(IMS:IME,JMS:JME,LM),INTENT(IN) :: PHINT,ZINT
      REAL,DIMENSION(IMS:IME,JMS:JME,LM),INTENT(IN) :: U,V,T

      REAL,DIMENSION(IMS:IME,JMS:JME,LM),INTENT(IN):: Q,QC,QI
      LOGICAL,INTENT(IN) :: F_QC,F_QI
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: MIXHT,PBLH,QSFC
      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: PBLK
!
!
      REAL,DIMENSION(IMS:IME,JMS:JME,LM)                               &
     &    ,INTENT(OUT) ::                                              &
     &                                         RQCBLTEN,RQIBLTEN       &
     &                                        ,RUBLTEN,RVBLTEN         &
     &                                        ,RTHBLTEN,RQBLTEN

!
!
!
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: I,J,K,KFLIP,K1,nvdiff,ntcw, levs , ipr
      INTEGER, DIMENSION(1) :: kpbl, kinver
      REAL(kind=kind_phys) :: dtp, surface,xkzm_m,xkzm_h,xkzm_s
      REAL(kind=kind_phys), DIMENSION(1,LM-1) :: dkt
      REAL(kind=kind_phys), DIMENSION(1,LM) :: dvdt,dudt,dtdt,ugrs,vgrs,tgrs,swh,  &
                                   hlw,del,prsl,prslk,phil
      REAL(kind=kind_phys), DIMENSION(1,LM+1) :: prsi,phii,prsik
      REAL(kind=kind_phys), DIMENSION(1)       :: xmu, psk,rb,ffmm,ffhh,tsea,qss,hflx,& 
                                   evap,stress,wind, &
                                   dusfc1,dvsfc1,dtsfc1,dqsfc1,hpbl,gamt,gamq,pii
      REAL(kind=kind_phys), DIMENSION(1,LM,3) :: dqdt ,qgrs
       real a96, a97, temp1, plow, tz0,seamask,qz0ss
      REAL :: QKLOW, CWMKLOW,RHOKLOW,QFC1,EXNSFC, PSFC, THSK, zmid1
      LOGICAL :: lpr, lprnt, dspheat
!----------------------------------------------------------------------
    !  lpr=.true.  !.false.
      lpr=.false.  !.false.
      ipr=0
      lprnt=.false.
      dspheat=.false.   !-- Dissipative heating (supported in newer version of GFS PBL)
    !  dtp=DT*float(NPHS)
!
!-- Time step is reduced by 0.5 because it is doubled (2x) in subroutine moninq (BSF, 4 Dec 2014)
!   
      dtp=0.5*DT*float(NPHS)

      levs = LM
      nvdiff = 3 
      ntcw  = 2         !-- Mix only cloud water; mixing of other species may not be robust
   !
      kinver(1) = levs          !! temp
      xkzm_m = 3.0
     !! xkzm_h = 1.0
      xkzm_h = 0.05  ! 0.1  !0.0  !1.0 !0.1  !0.2 !#0.5
      xkzm_s = 0.2              !! background diffusivity, see compns_physics.f in gfs/phys
!--------------------------------------------------------------
!.......................................................................
!$omp parallel do                &
!$omp     private(k,j,i)
!.......................................................................
      DO K=1,LM
      DO J=JMS,JME
      DO I=IMS,IME
      RQBLTEN(I,J,K) = 0.0
      RQCBLTEN(I,J,K) = 0.0
      RQIBLTEN(I,J,K) = 0.0
      RTHBLTEN(I,J,K) = 0.0
      RUBLTEN(I,J,K) = 0.0
      RVBLTEN(I,J,K) = 0.0
      ENDDO
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................

!.......................................................................
!$omp parallel do                &
!$omp     private(j,i)
!.......................................................................
      DO J=JMS,JME
      DO I=IMS,IME
       PBLH(I,J) = 0.9
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! THIS PART FOLLOWS MYJPBL TO UPDATE QSFC AND QZ0
!!!  NOTE: THIS PART IS SUPPOSED TO BE DONE CORRECTLY IN JSFC.F90, BUT IT IS NOT. THIS IS 
!!!  WHY WE UPDATE HERE. BE CAREFUL ABOUT DIFFERENCES BETWEEN SPECIFIC HUMIDITY AND 
!!!  MIXING RATIO.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
!--------------------------------------------------------------
!.......................................................................
!$omp parallel do                &
!$omp     private(j,i,k,qklow,cwmklow,rhoklow,thsk,seamask,qfc1,psfc,exnsfc)
!.......................................................................
      DO J=JMS,JME
      DO I=IMS,IME
           K=LM
           QKLOW=Q(I,J,K)
           CWMKLOW=QC(I,J,K)
           IF(F_QI)CWMKLOW=QC(I,J,K)+QI(I,J,K)
           RHOKLOW=PHMID(I,J,K)/(RD99*T(I,J,K)*(1.+P608*QKLOW-CWMKLOW))
           THSK=TSK(I,J)*(1.E5/PHINT(I,J,LM+1))**CAPPA

  !
  !***  COUNTING DOWNWARD FROM THE TOP, THE EXCHANGE COEFFICIENTS AKH
  !***  ARE DEFINED ON THE BOTTOMS OF THE LAYERS KTS TO KTE-1.  THESE COEFFICIENTS
  !***  ARE ALSO MULTIPLIED BY THE DENSITY AT THE BOTTOM INTERFACE LEVEL.
  !
  !
  !
                SEAMASK=XLAND(I,J)-1.
                THZ0(I,J)=(1.-SEAMASK)*THSK+SEAMASK*THZ0(I,J)
  !
                IF(SEAMASK<0.5)THEN
                  QFC1=XLV*CHKLOWQ(I,J)*AKHS(I,J)*RHOKLOW
  !
                  IF(SNOW(I,J)>0..OR.SICE(I,J)>0.5)THEN
                    QFC1=QFC1*RLIVWV
                  ENDIF
  !
                  IF(QFC1>0.)THEN
                    QSFC(I,J)=QKLOW+LHEAT(I,J)/QFC1
                  ENDIF
  !
                ELSE
                  PSFC=PHINT(I,J,LM+1)
                  EXNSFC=(1.E5/PSFC)**CAPPA
  
                 QSFC(I,J)=PQ0SEA/PSFC                                      &
          &         *EXP(AA2*(THSK-AA3*EXNSFC)/(THSK-AA4*EXNSFC))
               ENDIF
 !
               QZ0 (I,J)=(1.-SEAMASK)*QSFC(I,J)+SEAMASK*QZ0 (I,J)
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
      if(lpr) write(0,*)'new qsfc,qz0,thz0',qsfc(35,17),qz0(35,17),thz0(35,17)         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF UPDATING QSFC, QZ0 !!!!!!!!!!!!!!!!!!!!!!!!!

!--------------------------------------------------------------
!.......................................................................
!$omp parallel do                                                                     &
!$omp     private(j,i,k,dvdt,dudt,dtdt,k1,dqdt,kflip,prsi,prsik,phii,ugrs,vgrs,       &
!$omp            tgrs,qgrs,del,prsl,prslk,swh,hlw,zmid1,phil,xmu,surface,seamask,plow,&
!$omp             tz0,hflx,evap,rb,ffmm,ffhh,tsea,qss,wind,stress,kpbl,dusfc1,dvsfc1, &
!$omp             dqsfc1,hpbl,gamt,gamq,dkt,pii)
!.......................................................................
      DO J=JMS,JME
      DO I=IMS,IME

           DO K=1,LM  
            dvdt(1,K) = 0.0
            dudt(1,K) = 0.0
            dtdt(1,K) = 0.0
             DO k1=1,nvdiff 
              dqdt(1,K,K1) = 0.0
             ENDDO
           ENDDO
         
           DO K=1,LM+1
              KFLIP = LM+1+1-K
              prsi(1,K)  = PHINT(I,J,KFLIP)   !! pa
              prsik(1,K) = (prsi(1,K)*1.e-5)**CAPPA
              phii(1,K)  = ZINT(I,J,KFLIP)*G99
           ENDDO
          DO K=1,LM
             KFLIP = LM+1-K
              ugrs(1,K) = U(I,J,KFLIP)
              vgrs(1,K) = V(I,J,KFLIP)
              tgrs(1,K) = T(I,J,KFLIP)
              qgrs(1,K,1) = Q(I,J,KFLIP)
              qgrs(1,K,2) = QC(I,J,KFLIP)
              if(F_QI) qgrs(1,K,3) = QI(I,J,KFLIP)

              del(1,K)  = DP(I,J,KFLIP)     !! pa
              prsl(1,K) = PHMID(I,J,KFLIP)  !! pa
              prslk(1,K)= (prsl(1,K)*1.0e-5)**CAPPA
             !! phil(1,K) = G99*0.5*(ZINT(I,J,KFLIP)+ZINT(I,J,KFLIP+1))
             !! phil(1,K) = 0.5*(phii(1,K)+phii(1,K+1)) 
              swh(1,K) = RSWTT(I,J,KFLIP)    !!0.0  
              hlw(1,K) = RLWTT(I,J,KFLIP)    !!0.0 
                 zmid1=zint(i,j,kflip+1)+phmid(i,j,kflip)/airden(i,j,kflip)/g99 &
                          *alog(phint(i,j,kflip+1)/phmid(i,j,kflip))
              !!write(0,*)'K=',phil(1,k)/g99
                 phil(1,K)=zmid1*G99 
              !!write(0,*)'K=,new',phil(1,k)/g99
           ENDDO

              xmu(1) = factrs(I,J) 

              surface=phii(1,1)
              phil(1,:)=phil(1,:)-surface  !!phii(1,1)  !! surface=0
              phii(1,:)=phii(1,:)-surface  !!phii(1,1)

           seamask = xland(i,j) - 1.0
           plow    = phint(i,j,LM+1)
           tz0     = thz0(i,j)*(plow*1.0e-05)**CAPPA
           hflx(1) = SHEAT(I,J)/AIRDEN(I,J,LM)/CP99            ! W/m2 to K m/s
           evap(1) = LHEAT(I,J)/AIRDEN(I,J,LM)/XLV
           

           rb(1) = max(RIB(I,J),-5000.0)
           ffmm(1) =USTAR(I,J)*VKARMAN/AKMS(I,J)
           ffhh(1) =USTAR(I,J)*VKARMAN/AKHS(I,J)
           !  ffmm(1) = alog(phil(1,1)/9.8/0.05)
           !  ffhh(1) = ffmm(1)
           tsea(1) = 0.0    ! not in use
           qss(1)  = 0.0    ! not in use
           WIND(1)   = SQRT(ugrs(1,1)*ugrs(1,1)+vgrs(1,1)*vgrs(1,1))
           WIND(1)   = max(WIND(1),1.0d0)
           stress(1) = USTAR(I,J)*USTAR(I,J)
           KPBL(1)  = 1.0
           dusfc1(1) = 0.0
           dvsfc1(1) = 0.0
           dtsfc1(1) = 0.0
           dqsfc1(1) = 0.0
           hpbl(1)   = phil(1,1)/9.8 !! 0.0
           gamt(1)   = 0.0
           gamq(1)   = 0.0
           dkt(1,:)    = 0.0 
           pii(1)    = 1.0 
          call moninq(1,1,levs,nvdiff,ntcw,dvdt,dudt,dtdt,dqdt,         &
     &     ugrs,vgrs,tgrs,qgrs,swh,hlw,xmu,                             &
  !!   &     prsik(1,1),rb,ffmm,ffhh,tsea,qss,hflx,evap,stress,wind,kpbl, &
     &     pii,rb,ffmm,ffhh,tsea,qss,hflx,evap,stress,wind,kpbl, &
     &     prsi,del,prsl,prslk,phii,phil,dtp,dspheat,                   &
     &     dusfc1,dvsfc1,dtsfc1,dqsfc1,hpbl,gamt,gamq,dkt,              &
     &     kinver,xkzm_m,xkzm_h,xkzm_s                                  &    
     &    ,lprnt,ipr)
!! AFTER CALLING, flip Z , then back to NEMS variables                 
          
           DO K=1,LM
             KFLIP = LM+1-K
             RUBLTEN(I,J,K)  = dudt(1,KFLIP) 
             RVBLTEN(I,J,K)  = dvdt(1,KFLIP)                      
             RTHBLTEN(I,J,K) = dtdt(1,KFLIP)/prslk(1,KFLIP)  !! /EXNER(I,J,K)
        !!     RTHBLTEN(I,J,K) = dtdt(1,KFLIP)*prsik(1,1)/prslk(1,KFLIP)  !! /EXNER(I,J,K)
             RQBLTEN(I,J,K) = dqdt(1,KFLIP,1)
             RQCBLTEN(I,J,K) = dqdt(1,KFLIP,2)
             if(F_QI) RQIBLTEN(I,J,K) = dqdt(1,KFLIP,3)
          ENDDO

             PBLH(I,J)  = hpbl(1)
             PBLK(I,J)  = LM+1-kpbl(1)
             MIXHT(I,J)  = hpbl(1)
 
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................

!----------------------------------------------------------------------
      END SUBROUTINE GFSPBL
!----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END MODULE MODULE_BL_GFSPBL
!
!-----------------------------------------------------------------------
