!-----------------------------------------------------------------------------
     MODULE MODULE_MP_GFS
!-----------------------------------------------------------------------------
!
!    12-10-2010   Created by Weiguo Wang
!-- 7 May 2013 changes:
!   1) rh00=0.95 rather than 0.85 (critical threshold for the onset of condensation)
!   2) Update SR (snow ratio) array from the GFS microphysics and pass to the rest of the model
!
       USE MACHINE , ONLY : kind_phys
       USE FUNCPHYS , ONLY : gpvs, fpvs
       IMPLICIT NONE
       REAL, Private, PARAMETER ::                   & 
                                g99=9.80665,         &
                                t0c=273.15,          &
                                t_ice=-40.0+t0c

      PUBLIC :: GFSMP, GFSMP_INIT
      CONTAINS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE GFSMP    (DT,                                         &
                           dz8w,rho_phy,p_phy,pi_phy,th_phy,           &
                           SR,QT,F_ICE_phy,                            &
                           RAINNC,RAINNCV,                             &
                           Q,QC,QI,F_QC,F_QI,                          &
                           TP1,QP1,PSP1,                               &
                           ims,ime, jms,jme, lm)
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------

      INTEGER,INTENT(IN) :: IMS,IME,JMS,JME,LM

      REAL, INTENT(IN)   :: DT
      REAL, INTENT(IN),     DIMENSION(ims:ime, jms:jme, lm)::      &
                           dz8w,p_phy,pi_phy,rho_phy
      REAL, INTENT(INOUT),  DIMENSION(ims:ime, jms:jme, lm)::      &
                           th_phy,F_ICE_phy, QT
      REAL, INTENT(INOUT),  DIMENSION(ims:ime,jms:jme)           ::     &
                                                         RAINNC,RAINNCV
      REAL, INTENT(OUT),    DIMENSION(ims:ime,jms:jme):: SR
      REAL, INTENT(INOUT),  DIMENSION(IMS:IME,JMS:JME,LM):: Q,QC,QI
      LOGICAL,INTENT(IN) :: F_QC,F_QI
      REAL, INTENT(INOUT),  DIMENSION(ims:ime,jms:jme)           :: PSP1
      REAL, INTENT(INOUT),  DIMENSION(ims:ime,jms:jme,LM)     :: TP1,QP1

!-----------------------------------------------------------------------
!     LOCAL VARS
!-----------------------------------------------------------------------

       REAL(kind=kind_phys), PARAMETER    :: rh00 = 0.950  ! was 0.850
! Zhao scheme default opr value

       REAL(kind=kind_phys) psautco, prautco, evpco, wminco(2)
       data psautco /4.0E-4/, prautco /1.0E-4/
       data evpco   /2.0E-5/, wminco  /1.0E-5, 1.0E-5/

!      REAL(kind=kind_phys), PARAMETER    :: psautco = 4.0E-4    &  ! Zhao scheme default opr value
!                                           ,prautco = 1.0E-4    &  ! Zhao scheme default opr value
!                                           ,evpco  = 2.0E-5

!     TLATGS_PHY,TRAIN_PHY,APREC,PREC,ACPREC,SR are not directly related 
!     the microphysics scheme. Instead, they will be used by Eta precip 
!     assimilation.

      REAL,  DIMENSION( ims:ime, jms:jme,lm ) ::                  &
            TLATGS_PHY,TRAIN_PHY
      REAL,  DIMENSION(ims:ime,jms:jme):: APREC,PREC,ACPREC
      REAL,  DIMENSION( ims:ime, jms:jme,lm ) :: t_phy
      INTEGER :: I,J,K, KFLIP, L, KM
!!!LOCAL VARS FOR GFS MICROPHY
       Integer, parameter :: IX=1, IM=1, ipr=1
       REAL(kind=kind_phys), DIMENSION(IX,LM) :: PRSL, Q_COL, CWM_COL,T_COL, RHC, DELP &
                                 ,TP1_COL,QP1_COL,TP2_COL,QP2_COL
       REAL(kind=kind_phys), DIMENSION(IX) :: PS ,rain1, psp1_1,psp2_1, sratio
       REAL(kind=kind_phys) :: dtp,frain,fice
       INTEGER, dimension(ix,lm) :: IW            !! ice flag
       logical lprnt
       REAL(kind=kind_phys), DIMENSION(IX,LM) :: RAINP            ! not in use
       logical diag
        lprnt = .false.
        diag = .false.  !.true.
!------------------------------------------------------------------------
!**********************************************************************
        KM = LM
!
!-- Because the NMMB does not use leapfrog time differencing, only arrays from
!   the previous time step are needed (TP1, QP1, PSP1) and the time step (DTP)
!   used for physics rates does *NOT* need to be doubled (BSF, 03-30-2011).
!
        DTP = DT   !- was 2.0*DT
        frain = DT/DTP
               DO k = 1,lm
                rhc(ix,k)=rh00
               enddo
!.......................................................................
!$omp parallel do                &
!$omp     private(k,j,i)
!.......................................................................

               DO k = 1,lm
               DO j = jms,jme
               DO i = ims,ime
                t_phy(i,j,k) = th_phy(i,j,k)*pi_phy(i,j,k)
                TLATGS_PHY (i,j,k)=0.
                TRAIN_PHY  (i,j,k)=0.
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
      DO j = jms,jme
      DO i = ims,ime
        ACPREC(i,j)=0.
        APREC (i,j)=0.
        PREC  (i,j)=0.
        SR    (i,j)=0.
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!-- Start of original driver for EGCP01COLUMN
!-----------------------------------------------------------------------
!.......................................................................
!$omp parallel do                &
!$omp     private(j,i,k,kflip,delp,prsl,q_col,cwm_col,t_col,tp1_col, &
!$omp             qp1_col,tp2_col,qp2_col,psp1_1,psp2_1,iw,ps,rain1, &
!$omp             fice,sratio)
!.......................................................................
       DO J=JMS,JME    
        DO I=IMS,IME  
          DO K=1,LM
            KFLIP = LM + 1 -K

            DELP(IX,KFLIP)=RHO_PHY(I,J,K)*g99*dz8w(I,J,K)
            PRSL(IX,KFLIP)=P_phy(I,J,K)
            Q_COL(IX,KFLIP) = Q(I,J,K)
            CWM_COL(IX,KFLIP)=QC(I,J,K)+QI(I,J,K) 
            T_COL(IX,KFLIP) = t_phy(i,j,k) 
!
!-- The original GFS uses a leapfrog time-differencing scheme that goes back
!   2 time steps, represented by arrays with the names TP1, QP1, & PSP1.
!   The arrays with the names TP2, QP2, & PSP2 are associated with the previous
!   time step.  But because the NMMB does not use leapfrog time differencing,
!   the TP1, QP1, PSP1 set of arrays now represent the previous time step,
!   and the TP2, QP2, PSP2 set of arrays are treated as dummy input values to
!   subroutine gscond, and they will also store the T, Q, & P values for the
!   previous time step (BSF, 03-30-2011).
!
!-- The situation is a little different after leaving subroutine gscond, so 
!   please read the next set of comments after the "300 continue" line below.
!   
            TP1_COL(IX,K) = tp1(i,j,k)
            QP1_COL(IX,K) = qp1(i,j,k)
            TP2_COL(IX,K) = tp1(i,j,k)   !-- was tp2(i,j,k)
            QP2_COL(IX,K) = qp1(i,j,k)   !-- was qp2(i,j,k)
            psp1_1(IX)   = psp1(i,j)
            psp2_1(IX)   = psp1(i,j)     !-- was psp2(i,j,k)
            iw(ix,KFLIP) = 0 !F_ICE_phy(I,J,K) 
          ENDDO

            PS(IX) = DELP(IX,1)*0.5+PRSL(IX,1)
            rain1(ix) = 0.0
            sratio(ix) = 0.0

               IF(DIAG .and. i == 10 .and. j == 10 ) THEN
                  write(0,*)'before calling MICRO'
                  write(0,*)'tp1=',tp1(i,j,:)
                  write(0,*)'qp1=',qp1(i,j,:)
                  write(0,*)'psp1=',psp1(i,j)
                  write(0,*)'qt=',qt(i,j,:)
                  write(0,*)'cwm=',cwm_col(1,:)
                  write(0,*)'qc=',qc(i,j,:)
                  write(0,*)'qi=',qi(i,j,:)
               ENDIF    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  CALL MICROPHY                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           call gscond(im, ix, KM, dtp, dtp, prsl, ps,                &
                       Q_COL, CWM_COL, T_COL,                           &
                       TP1_COL, QP1_COL,  PSP1_1,                     &
                       TP2_COL, QP2_COL,  PSP2_1,                     &
                       rhc,lprnt, ipr)
           call precpd_nmmb(im, ix, KM, dtp, delp, prsl, ps,           &
                       Q_COL, CWM_COL, T_COL, rain1, sratio,           &
                       rainp, rhc, psautco, prautco, evpco, wminco,    &
                       lprnt, ipr)
 300      continue 

!
!-- Before exiting gscond, the arrays TP1_COL, QP1_COL, PSP1_1 are set
!   to the previous time step, while TP2_COL, QP2_COL, PSP2_1 are set
!   to values at the current time step.  The arrays below will be set
!   to values associated with TP2_COL, QP2_COL, PSP2_1, while the
!   TP1_COL, QP1_COL, PSP1_1 arrays will not be used (BSF, 03-30-2011).
!   
            DO K=1,lm
             tp1(i,j,k) = tp2_COL(ix,k)   !- was tp1_COL(ix,k)
             qp1(i,j,k) = qp2_COL(ix,k)   !- was qp1_COL(ix,k)
             psp1(i,j)  = psp2_1(ix)      !- was psp1_1(ix)
            ENDDO
!#######################################################################
!
!--- Update storage arrays
!
          DO L=1,KM
            KFLIP = KM + 1 - L
            TRAIN_phy(I,J,L)= (T_col(1,KFLIP)-T_phy(I,J,L))/DTp
            TLATGS_phy(I,J,L)=(T_col(1,KFLIP)-T_phy(I,J,L))*frain
          ENDDO
          DO K=1,KM
            KFLIP=KM + 1 - K 
            T_phy(I,J,K)=T_col(1,KFLIP)
            Q(I,J,K)= Q_col(1,KFLIP)
              fice=1.0
              IF(T_COL(1,KFLIP) .GT. t_ice .and.                   &
                 T_COL(1,KFLIP) .LE. t0c ) THEN
                 fice = 1.0 - (T_COL(1,KFLIP)-t_ice)/(t0c-t_ice)
              ENDIF
              IF(T_COL(1,KFLIP) .GT. t0c ) fice=0.0
      !      fice = float( IW(1,KFLIP) )
            QC(I,J,K) = CWM_COL(1,KFLIP)*(1.0-fice)
            QI(I,J,K) = CWM_COL(1,KFLIP)*fice
            QT(I,J,K) = CWM_COL(1,KFLIP)
            F_ICE_phy(I,J,K) = fice
          ENDDO
!
!--- Update accumulated precipitation statistics
!
!--- Surface precipitation statistics; SR is fraction of surface 
!    precipitation (if >0) associated with snow
!
        APREC(I,J)=rain1(1)*frain       ! Accumulated surface precip (depth in m)  !<--- Ying
        PREC(I,J)=PREC(I,J)+APREC(I,J)
        ACPREC(I,J)=ACPREC(I,J)+APREC(I,J)
        SR(I,J)=sratio(ix)
!
       IF(DIAG .and. i == 10 .and. j == 10 ) THEN
          write(0,*)'RAIN=',APREC(I,J)
          write(0,*)'DELP=',DELP  
          write(0,*)'PS=',PS
          write(0,*)'tp1=',tp1(i,j,:)
          write(0,*)'qp1=',qp1(i,j,:)
          write(0,*)'psp1=',psp1(i,j)
          write(0,*)'p,cwm,T,Fice,Q'
           do k=1,km
            write(0,100)prsl(1,k),cwm_col(1,k),t_col(1,k),f_ice_phy(i,j,km+1-K),q_col(1,k),fpvs(t_col(1,k))
           enddo
          write(0,*)'Max,min T',maxval(tp1),minval(tp1)
       ENDIF
100    format(F10.1,E12.3,F6.1,F4.1,2E12.3)
    enddo                          ! End "I" loop
    enddo                          ! End "J" loop
!.......................................................................
!$omp end parallel do
!.......................................................................
!.......................................................................
!$omp parallel do                &
!$omp     private(k,j,i)
!.......................................................................
     DO k = 1,lm
       DO j = jms,jme
       DO i = ims,ime
         th_phy(i,j,k) = t_phy(i,j,k)/pi_phy(i,j,k)
       ENDDO   !- i
       ENDDO   !- j
     ENDDO   !- k
!.......................................................................
!$omp end parallel do
!.......................................................................
! 
!- Update rain (convert from m to kg/m**2, which is also equivalent to mm depth)
! 
       DO j=jms,jme
       DO i=ims,ime
          RAINNC(i,j)=APREC(i,j)*1000.+RAINNC(i,j)
          RAINNCV(i,j)=APREC(i,j)*1000.
       ENDDO
       ENDDO
!
!-----------------------------------------------------------------------
!
  END SUBROUTINE GFSMP
!
!-----------------------------------------------------------------------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE GFSMP_INIT
          CALL GPVS
        END SUBROUTINE GFSMP_INIT    

!
      END MODULE module_mp_gfs
