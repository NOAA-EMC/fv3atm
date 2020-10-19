 !----------------------------------------------------------------------
!
      MODULE MODULE_DIAGNOSE
!
!----------------------------------------------------------------------
!
      USE MODULE_KINDS
      USE MODULE_CONSTANTS, ONLY : R_D,R_V,CPV,CP,G,CLIQ,PSAT,P608      &
     & ,XLV,TIW,EPSQ,DBZmin
      IMPLICIT NONE
!
      REAL, PRIVATE, PARAMETER :: Cice=1.634e13    &  !-- For dry ice (T<0C)
     , Cwet=1./.189      &   !-- Wet ice spheres at >=0C (Smith, JCAM, 1984, p. 1259, eq. 10)
     , Cboth=Cice*Cwet   &   !-- Rain + wet ice at >0C
     , CU_A=300, CU_B=1.4    &   !-- For convective precipitation reflectivity
     , TFRZ=TIW, TTP=TIW+0.01, Zmin=0.01                                &
     , EPSILON=R_D/R_V, ONE_MINUS_EPSILON=1.-EPSILON                    &
     , R_FACTOR=1./EPSILON-1., CP_FACTOR=CPV/CP-1., RCP=R_D/CP          &
     , P00_INV=1.E-5, XA=(CLIQ-CPV)/R_V, XB=XA+XLV/(R_V*TTP)
!
!
!----------------------------------------------------------------------
!
      CONTAINS
!
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
!
      SUBROUTINE CALMICT(P1D,T1D,Q1D,C1D,FI1D,FR1D,FS1D,CUREFL, &
                       DBZ1,I,J, Ilook,Jlook)

      USE MODULE_MP_ETANEW, ONLY : FERRIER_INIT, FPVS,RQR_DRmin,  &
                                   RQR_DRmax,MASSI,CN0R0,         &
                                   CN0r_DMRmin,CN0r_DMRmax

      IMPLICIT NONE

      
      REAL, PARAMETER :: DMImin=.05e-3, DMImax=1.e-3, &
                         XMImin=1.e6*DMImin, XMImax=1.e6*DMImax
      INTEGER, PARAMETER :: MDImin=XMImin, MDImax=XMImax
!
!-----------------------------------------------------------------------
!
!--- Mean rain drop diameters vary from 50 microns to 1000 microns (1 mm)
!
      REAL, PARAMETER :: DMRmin=.05E-3, DMRmax=1.E-3, DelDMR=1.E-6,        &
         XMRmin=1.E6*DMRmin, XMRmax=1.E6*DMRmax
      INTEGER, PARAMETER :: MDRmin=XMRmin, MDRmax=XMRmax

      INTEGER INDEXS, INDEXR, I, J, Ilook, Jlook

      REAL ::  NLICE, N0r,Ztot,Zrain,Zice,Zconv
      REAL ::  P1D,T1D,Q1D,C1D,                                            &
               FI1D,FR1D,FS1D,CUREFL,                                      &
               QW1,QI1,QR1,QS1,                                            &
               DBZ1

      REAL :: FLARGE, FSMALL, WV, ESAT, TC, WC, RHO,                       &
              RRHO, RQR, Fice, Frain, Rimef , XLI, QICE, DRmm,             &
              DLI, XSIMASS, XLIMASS, DUM, WVQW, QSIGRD, QLICE, FLIMASS

      REAL, PARAMETER ::                                                   &
     &  RHgrd=1.                                                           &
     & ,T_ICE=-40.                                                         &
     & ,NLImax=5.E3                                                        &
     & ,NLImin=1.E3                                                        &
     & ,N0r0=8.E6                                                          &
     & ,N0rmin=1.E4

! ---------
      DBZ1=DBZmin
      IF (C1D<=EPSQ) THEN
!
!--- Skip rest of calculatiions if no condensate is present
!
        RETURN
      ELSE
        WC=C1D
      ENDIF
!
!--- Code below is from GSMDRIVE for determining:
!    QI1 - total ice (cloud ice & snow) mixing ratio
!    QW1 - cloud water mixing ratio
!    QR1 - rain mixing ratio
!
      Zrain=0.            !--- Radar reflectivity from rain
      Zice=0.             !--- Radar reflectivity from ice
      Zconv=CUREFL   !--- Radar reflectivity from convection
      QW1=0.
      QI1=0.
      QR1=0.
      QS1=0.
      TC=T1D-TFRZ
      Fice=FI1D
      Frain=FR1D
      IF (TC.LE.T_ICE .OR. Fice.GE.1.) THEN
        QI1=WC
      ELSE IF (Fice .LE. 0.) THEN
        QW1=WC
      ELSE
        QI1=Fice*WC
        QW1=WC-QI1
      ENDIF
!
      IF (QW1>0. .AND. Frain>0.) THEN
        IF (Frain>=1.) THEN
          QR1=QW1
          QW1=0.
        ELSE
          QR1=Frain*QW1
          QW1=QW1-QR1
        ENDIF
      ENDIF
      WV=Q1D/(1.-Q1D)
      RHO=P1D/(R_D*T1D*(1.+P608*Q1D))
      RRHO=1./RHO
  !
  !--- Based on code from GSMCOLUMN in model to determine reflectivity from rain
  !
rain_dbz: IF (QR1>EPSQ) THEN
        RQR=RHO*QR1
        IF (RQR<=RQR_DRmin) THEN
          N0r=MAX(N0rmin, CN0r_DMRmin*RQR)
          INDEXR=MDRmin
        ELSE IF (RQR>=RQR_DRmax) THEN
          N0r=CN0r_DMRmax*RQR
          INDEXR=MDRmax
        ELSE
          N0r=N0r0
          INDEXR=MAX( XMRmin, MIN(CN0r0*RQR**.25, XMRmax) )
        ENDIF
  !
  !--- INDEXR is the mean drop size in microns; convert to mm
  !
        DRmm=1.e-3*REAL(INDEXR)
        Zrain=0.72*N0r*DRmm*DRmm*DRmm*DRmm*DRmm*DRmm*DRmm
      ENDIF   rain_dbz     !--- End IF (QR1 .GT. EPSQ)
!
!--- Based on code from GSMCOLUMN in model to determine partition of
!    total ice into cloud ice & snow (precipitation ice)
!
ice_dbz: IF (QI1 .GT. EPSQ) THEN
        QICE=QI1
        RHO=P1D/(R_D*T1D*(1.+P608*Q1D))
        RRHO=1./RHO
!- FPVS - saturation vapor pressure w/r/t water ( >=0C ) or ice ( <0C ) in kPa
        ESAT=1000.*FPVS(T1D)      !-- saturation w/r/t ice at <0C in Pa
        QSIgrd=RHgrd*EPSILON*ESAT/(P1D-ESAT)
        WVQW=WV+QW1
!
! * FLARGE  - ratio of number of large ice to total (large & small) ice
! * FSMALL  - ratio of number of small ice crystals to large ice particles
!  ->  Small ice particles are assumed to have a mean diameter of 50 microns.
!  * XSIMASS - used for calculating small ice mixing ratio
!  * XLIMASS - used for calculating large ice mixing ratio
!  * INDEXS  - mean size of snow to the nearest micron (units of microns)
!  * RimeF   - Rime Factor, which is the mass ratio of total (unrimed &
!              rimed) ice mass to the unrimed ice mass (>=1)
!  * FLIMASS - mass fraction of large ice
!  * QTICE   - time-averaged mixing ratio of total ice
!  * QLICE   - time-averaged mixing ratio of large ice
!  * NLICE   - time-averaged number concentration of large ice
!
        IF (TC>=0. .OR. WVQW<QSIgrd) THEN
          FLARGE=1.
        ELSE
          FLARGE=.2
          IF (TC>=-8. .AND. TC<=-3.) FLARGE=.5*FLARGE
        ENDIF
        FSMALL=(1.-FLARGE)/FLARGE
        XSIMASS=RRHO*MASSI(MDImin)*FSMALL
!
        DUM=XMImax*EXP(.0536*TC)
        INDEXS=MIN(MDImax, MAX(MDImin, INT(DUM) ) )
        RimeF=AMAX1(1., FS1D )
        XLIMASS=RRHO*RimeF*MASSI(INDEXS)
        FLIMASS=XLIMASS/(XLIMASS+XSIMASS)
        QLICE=FLIMASS*QICE
        NLICE=QLICE/XLIMASS
new_nlice: IF (NLICE<NLImin .OR. NLICE>NLImax) THEN
!
!--- Force NLICE to be between NLImin and NLImax
!
          DUM=MAX(NLImin, MIN(NLImax, NLICE) )
          XLI=RHO*(QICE/DUM-XSIMASS)/RimeF
          IF (XLI<=MASSI(MDImin) ) THEN
            INDEXS=MDImin
          ELSE IF (XLI<=MASSI(450) ) THEN
            DLI=9.5885E5*XLI**.42066         ! DLI in microns
            INDEXS=MIN(MDImax, MAX(MDImin, INT(DLI) ) )
          ELSE IF (XLI<=MASSI(MDImax) ) THEN
            DLI=3.9751E6*XLI**.49870         ! DLI in microns
            INDEXS=MIN(MDImax, MAX(MDImin, INT(DLI) ) )
          ELSE
            INDEXS=MDImax
!
!--- 8/22/01: Increase density of large ice if maximum limits
!    are reached for number concentration (NLImax) and mean size
!    (MDImax).  Done to increase fall out of ice.
!
            IF (DUM>=NLImax) THEN
              RimeF=RHO*(QICE/NLImax-XSIMASS)/MASSI(INDEXS)
            ENDIF
          ENDIF             ! End IF (XLI .LE. MASSI(MDImin) )
          XLIMASS=RRHO*RimeF*MASSI(INDEXS)
          FLIMASS=XLIMASS/(XLIMASS+XSIMASS)
          QLICE=FLIMASS*QICE
          NLICE=QLICE/XLIMASS
        ENDIF  new_nlice
        QS1=AMIN1(QI1, QLICE)
        QI1=AMAX1(0., QI1-QS1)
   !
   !--- Equation (C.8) in Ferrier (1994, JAS, p. 272), which when
   !    converted from cgs units to mks units results in the same
   !    value for Cice, which is equal to the {} term below:
   !
   !    Zi={.224*720*(10**18)/[(PI*RHOL)**2]}*(RHO*QLICE)**2/NLICE,
   !    where RHOL=1000 kg/m**3 is the density of liquid water
   !
   !--- Valid only for exponential ice distributions
   !
         IF (NLICE>0. .AND. QLICE>0.) THEN
           Zice=Cice*RHO*RHO*QLICE*QLICE/NLICE
         ENDIF
      ENDIF  ice_dbz         ! End IF (QI1 .GT. EPSQ)
!
!---  Calculate total (convective + grid-scale) radar reflectivity
!
      Ztot=Zrain+Zice+Zconv
      IF (Ztot>Zmin) DBZ1= 10.*ALOG10(Ztot)
      RETURN
      END SUBROUTINE CALMICT
!
!-----------------------------------------------------------------------
!
      SUBROUTINE MAX_FIELDS_driver(T,Q,U              &
                           ,V,CW                      &
                           ,F_RAIN,F_ICE              &
                           ,F_RIMEF                   &
                           ,Z_DAM,W,REFL_10CM         &
                           ,QR,QS,QG                  &
                           ,PINT,PREC                 &
                           ,CPRATE,HTOP               &
                           ,T2,U10,V10                &
                           ,PSHLTR,TSHLTR,QSHLTR      &
                           ,PMID                      &
                           ,REFDMAX,PRATEMAX          &
                           ,FPRATEMAX,SR              &
                           ,UPVVELMAX,DNVVELMAX       &
                           ,TLMAX,TLMIN               &
                           ,T02MAX,T02MIN             &
                           ,RH02MAX,RH02MIN           &
                           ,U10MAX,V10MAX,TH10,T10    &
                           ,SPD10MAX,T10AVG,PSFCAVG   &
                           ,AKHS,AKMS                 &
                           ,AKHSAVG,AKMSAVG           &
                           ,SNO,SNOAVG                &
                           ,UPHLMAX                   &
                           ,DT,NPHS,NTSD              &
                           ,NSTEPS_PER_RESET,FIS      &
                           ,IMS,IME                   &
                           ,LM,NCOUNT,MICROPHYSICS)

      USE MODULE_CONSTANTS ,ONLY : DBZmin

      IMPLICIT NONE

      INTEGER,PARAMETER :: JMS=1,JME=1

      CHARACTER *25,INTENT(IN) :: MICROPHYSICS

      INTEGER,INTENT(IN) :: IMS,IME,LM,NTSD,NSTEPS_PER_RESET
      INTEGER,INTENT(IN) :: NPHS

      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: NCOUNT

      REAL, DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: T, Q, U, V, CW   &
                                              ,F_RAIN,F_ICE,F_RIMEF        &
                                              ,QR,QS,QG,W,Z_DAM,REFL_10CM

      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM+1),INTENT(IN) :: PINT

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: PREC,CPRATE,HTOP       &
                                                    ,T2,U10,V10            &
                                                    ,PSHLTR,TSHLTR,QSHLTR  &
                                                    ,TH10,AKHS             &
                                                    ,AKMS,SNO,FIS,SR

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT):: REFDMAX,PRATEMAX      &
                                                   ,FPRATEMAX               &
                                                   ,UPVVELMAX,DNVVELMAX     &
                                                   ,TLMAX,TLMIN             &
                                                   ,T02MAX,T02MIN           &
                                                   ,RH02MAX,RH02MIN         &
                                                   ,U10MAX,V10MAX           &
                                                   ,SPD10MAX,T10AVG,PSFCAVG &
                                                   ,UPHLMAX,T10,AKHSAVG     &
                                                   ,AKMSAVG,SNOAVG

      REAL, DIMENSION(IMS:IME,LM),INTENT(IN) :: PMID

      REAL, INTENT(IN) :: DT

      INTEGER:: I,J
!
!-----------------------------------------------------------------------
!***  Start here
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!***  At the appropriate times, reset the various min/max/average
!***  diagnostic fields to begin accumulating for the next period
!-----------------------------------------------------------------------
!
      IF(NTSD == 0 .or. MOD(NTSD,NSTEPS_PER_RESET)==0) THEN
        DO J=JMS,JME
        DO I=IMS,IME
          TLMAX(I,J)     = -999.
          TLMIN(I,J)     =  999.
          T02MAX(I,J)    = -999.
          T02MIN(I,J)    =  999.
          RH02MAX(I,J)   = -999.
          RH02MIN(I,J)   =  999.
          SPD10MAX(I,J)  = -999.
          U10MAX(I,J)    = -999.
          V10MAX(I,J)    = -999.
          UPVVELMAX(I,J) = -999.
          DNVVELMAX(I,J) =  999.
          t10avg(I,J)    =    0.
          T10(I,J)       =    0.
          psfcavg(I,J)   =    0.
          akhsavg(I,J)   =    0.
          akmsavg(I,J)   =    0.
          snoavg(I,J)    =    0.
          REFDMAX(I,J)   = DBZmin
          PRATEMAX(I,J)  =    0
          FPRATEMAX(I,J) =    0
          UPHLMAX(I,J)   = -999.
          NCOUNT(I,J)    =    0
        ENDDO
        ENDDO
      ENDIF
!
      IF (TRIM(MICROPHYSICS) == 'fer') THEN
        CALL MAX_FIELDS(T,Q,U                     &
                       ,V,CW                      &
                       ,F_RAIN,F_ICE              &
                       ,F_RIMEF                   &
                       ,Z_DAM,W,PINT,PREC         &
                       ,CPRATE,HTOP               &
                       ,T2,U10,V10                &
                       ,PSHLTR,TSHLTR,QSHLTR      &
                       ,PMID                      &
                       ,REFDMAX,PRATEMAX          &
                       ,FPRATEMAX,SR              &
                       ,UPVVELMAX,DNVVELMAX       &
                       ,TLMAX,TLMIN               &
                       ,T02MAX,T02MIN             &
                       ,RH02MAX,RH02MIN           &
                       ,U10MAX,V10MAX,TH10,T10    &
                       ,SPD10MAX,T10AVG,PSFCAVG   &
                       ,AKHS,AKMS                 &
                       ,AKHSAVG,AKMSAVG           &
                       ,SNO,SNOAVG                &
                       ,UPHLMAX                   &
                       ,DT,NPHS,NTSD              &
                       ,FIS                       &
                       ,IMS,IME,JMS,JME           &
                       ,LM,NCOUNT)
      ELSEIF (TRIM(MICROPHYSICS) == 'fer_hires') THEN
        CALL MAX_FIELDS_HR(T,Q,U                  &
                       ,V,CW                      &
                       ,F_RAIN,F_ICE              &
                       ,F_RIMEF                   &
                       ,Z_DAM,W,REFL_10CM         &
                       ,PINT,PREC                 &
                       ,CPRATE,HTOP               &
                       ,T2,U10,V10                &
                       ,PSHLTR,TSHLTR,QSHLTR      &
                       ,PMID                      &
                       ,REFDMAX,PRATEMAX          &
                       ,FPRATEMAX,SR              &
                       ,UPVVELMAX,DNVVELMAX       &
                       ,TLMAX,TLMIN               &
                       ,T02MAX,T02MIN             &
                       ,RH02MAX,RH02MIN           &
                       ,U10MAX,V10MAX,TH10,T10    &
                       ,SPD10MAX,T10AVG,PSFCAVG   &
                       ,AKHS,AKMS                 &
                       ,AKHSAVG,AKMSAVG           &
                       ,SNO,SNOAVG                &
                       ,UPHLMAX                   &
                       ,DT,NPHS,NTSD              &
                       ,FIS                       &
                       ,IMS,IME,JMS,JME           &
                       ,LM,NCOUNT)
      ELSEIF (TRIM(MICROPHYSICS) == 'wsm6') THEN
        CALL MAX_FIELDS_w6(T,Q,U,V,Z_DAM,W        &
                       ,QR,QS,QG,PINT,PREC        &
                       ,CPRATE,HTOP               &
                       ,T2,U10,V10                &
                       ,PSHLTR,TSHLTR,QSHLTR      &
                       ,PMID                      &
                       ,REFDMAX,PRATEMAX          &
                       ,FPRATEMAX,SR              &
                       ,UPVVELMAX,DNVVELMAX       &
                       ,TLMAX,TLMIN               &
                       ,T02MAX,T02MIN             &
                       ,RH02MAX,RH02MIN           &
                       ,U10MAX,V10MAX,TH10,T10    &
                       ,SPD10MAX,T10AVG,PSFCAVG   &
                       ,AKHS,AKMS                 &
                       ,AKHSAVG,AKMSAVG           &
                       ,SNO,SNOAVG                &
                       ,UPHLMAX                   &
                       ,DT,NPHS,NTSD              &
                       ,FIS                       &
                       ,IMS,IME,JMS,JME           &
                       ,LM,NCOUNT)
      ELSEIF (TRIM(MICROPHYSICS) == 'thompson') THEN
        CALL MAX_FIELDS_THO(T,Q,U,V,Z_DAM,W       &
                       ,REFL_10CM,PINT,PREC       &
                       ,CPRATE,HTOP               &
                       ,T2,U10,V10                &
                       ,PSHLTR,TSHLTR,QSHLTR      &
                       ,PMID                      &
                       ,REFDMAX,PRATEMAX          &
                       ,FPRATEMAX,SR              &
                       ,UPVVELMAX,DNVVELMAX       &
                       ,TLMAX,TLMIN               &
                       ,T02MAX,T02MIN             &
                       ,RH02MAX,RH02MIN           &
                       ,U10MAX,V10MAX,TH10,T10    &
                       ,SPD10MAX,T10AVG,PSFCAVG   &
                       ,AKHS,AKMS                 &
                       ,AKHSAVG,AKMSAVG           &
                       ,SNO,SNOAVG                &
                       ,UPHLMAX                   &
                       ,DT,NPHS,NTSD              &
                       ,FIS                       &
                       ,IMS,IME,JMS,JME           &
                       ,LM,NCOUNT)
      ENDIF
 
      END SUBROUTINE MAX_FIELDS_driver
!
!----------------------------------------------------------------------
!
      SUBROUTINE MAX_FIELDS(T,Q,U                     &
                           ,V,CW                      &
                           ,F_RAIN,F_ICE              &
                           ,F_RIMEF                   &
                           ,Z_DAM,W,PINT,PREC         &
                           ,CPRATE,HTOP               &
                           ,T2,U10,V10                &
                           ,PSHLTR,TSHLTR,QSHLTR      &
                           ,PMID                      &
                           ,REFDMAX,PRATEMAX          &
                           ,FPRATEMAX,SR              &
                           ,UPVVELMAX,DNVVELMAX       &
                           ,TLMAX,TLMIN               &
                           ,T02MAX,T02MIN             &
                           ,RH02MAX,RH02MIN           &
                           ,U10MAX,V10MAX,TH10,T10    &
                           ,SPD10MAX,T10AVG,PSFCAVG   &
                           ,AKHS,AKMS                 &
                           ,AKHSAVG,AKMSAVG           &
                           ,SNO,SNOAVG                &
                           ,UPHLMAX                   &
                           ,DT,NPHS,NTSD              &
                           ,FIS                       &
                           ,IMS,IME,JMS,JME           &
                           ,LM,NCOUNT)

      USE MODULE_MP_ETANEW, ONLY : FERRIER_INIT, FPVS0

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: IMS,IME,JMS,JME,LM,NTSD
      INTEGER,INTENT(IN) :: NPHS

      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: NCOUNT

      REAL, DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: T, Q, U, V, CW   &
                                              ,F_RAIN,F_ICE,F_RIMEF        &
                                              ,W,Z_DAM

      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM+1),INTENT(IN) :: PINT

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: PREC,CPRATE,HTOP       &
                                                    ,T2,U10,V10            &
                                                    ,PSHLTR,TSHLTR,QSHLTR  &
                                                    ,TH10,AKHS             &
                                                    ,AKMS,SNO,FIS,SR

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT):: REFDMAX,PRATEMAX      &
                                                   ,FPRATEMAX               &
                                                   ,UPVVELMAX,DNVVELMAX     &
                                                   ,TLMAX,TLMIN             &
                                                   ,T02MAX,T02MIN           &
                                                   ,RH02MAX,RH02MIN         &
                                                   ,U10MAX,V10MAX           &
                                                   ,SPD10MAX,T10AVG,PSFCAVG &
                                                   ,UPHLMAX,T10,AKHSAVG     &
                                                   ,AKMSAVG,SNOAVG

      REAL, DIMENSION(IMS:IME,LM),INTENT(IN) :: PMID

      REAL, INTENT(IN) :: DT

      INTEGER :: UPINDX(IMS:IME,JMS:JME)
      REAL, DIMENSION(IMS:IME,JMS:JME) :: P10

      REAL, DIMENSION(IMS:IME,JMS:JME) :: ZINTSFC
      REAL, DIMENSION(IMS:IME,JMS:JME,LM) :: Z

      REAL :: PLOW, PUP,WGTa,WGTb,ZMIDloc,ZMIDP1
      REAL :: P1Da,P1Db,P1D(2)
      REAL :: T1Da,T1Db,T1D(2),fact
      REAL :: Q1Da,Q1Db,Q1D(2)
      REAL :: C1Da,C1Db,C1D(2)
      REAL :: FR1Da,FR1Db,FR1D(2)
      REAL :: FI1Da,FI1Db,FI1D(2)
      REAL :: FS1Da,FS1Db,FS1D(2),DBZ1(2)

      REAL :: CUPRATE,CUREFL,CUREFL_I,ZFRZ,DBZ1avg,FCTR,DELZ,Z1KM,ZCTOP
      REAL :: T02, RH02, TERM
      REAL :: CAPPA_MOIST, VAPOR_PRESS, SAT_VAPOR_PRESS
      REAL, SAVE:: DTPHS, RDTPHS
      REAL :: MAGW2

      INTEGER :: LCTOP
      INTEGER :: I,J,L,LL, RC, Ilook,Jlook


!***  COMPUTE AND SAVE THE FACTORS IN R AND CP TO ACCOUNT FOR
!***  WATER VAPOR IN THE AIR.
!***
!***  RECALL: R  = Rd * (1. + Q * (1./EPSILON - 1.))
!***          CP = CPd * (1. + Q * (CPv/CPd - 1.))

      Ilook=99
      Jlook=275

      DTPHS=DT*NPHS
      RDTPHS=3.6e6/DTPHS
!
      DO J=JMS,JME
       DO I=IMS,IME
         ZINTSFC(I,J)=FIS(I,J)/g
         DO L=1,LM
           Z(I,J,L)=Z_DAM(I,J,L)/g !rv FV3 is in gpdam, not gpm
         ENDDO
       ENDDO
      ENDDO
!
!     WON'T BOTHER TO REBUILD HEIGHTS AS IS DONE IN POST.
!     THE NONHYDROSTATIC MID-LAYER Z VALUES MATCH CLOSELY ENOUGH
!     AT 1000 m AGL
!
      DO J=JMS,JME
       DO I=IMS,IME
 L_LOOP: DO L=1,LM-1
          PLOW= PMID(I,L+1)
          PUP=  PMID(I,L)
          IF (PLOW .ge. 40000. .and. PUP .le. 40000.) THEN
            UPINDX(I,J)=L
            exit L_LOOP
          ENDIF
         ENDDO L_LOOP
       ENDDO
      ENDDO
!
      DO J=JMS,JME
       DO I=IMS,IME
  vloop: DO L=8,LM-1
          IF ( (Z(I,J,L+1)-ZINTSFC(I,J)) .LE. 1000.                &
          .AND.(Z(I,J,L)-ZINTSFC(I,J))   .GE. 1000.)  THEN
            ZMIDP1=Z(I,J,L)
            ZMIDloc=Z(I,J,L+1)
            P1D(1)=PMID(I,L)
            P1D(2)=PMID(I,L+1)
            T1D(1)=T(I,J,L)
            T1D(2)=T(I,J,L+1)
            Q1D(1)=Q(I,J,L)
            Q1D(2)=Q(I,J,L+1)
            C1D(1)=CW(I,J,L)
            C1D(2)=CW(I,J,L+1)
            FR1D(1)=F_RAIN(I,J,L)
            FR1D(2)=F_RAIN(I,J,L+1)
            FI1D(1)=F_ICE(I,J,L)
            FI1D(2)=F_ICE(I,J,L+1)
            FS1D(1)=F_RIMEF(I,J,L)
            FS1D(2)=F_RIMEF(I,J,L+1)
            EXIT vloop
          ENDIF
        ENDDO vloop
!
!!! INITIAL CUREFL VALUE WITHOUT REDUCTION ABOVE FREEZING LEVEL
!
        CUPRATE=RDTPHS*CPRATE(I,J)
        CUREFL=0.
        IF (CUPRATE>0.) CUREFL=CU_A*CUPRATE**CU_B
        ZFRZ=Z(I,J,LM)
!
 culoop: IF (CUREFL>0. .AND. NINT(HTOP(I,J)) > 0) THEN
 vloop2:  DO L=1,LM
            IF (T(I,J,L) >= TFRZ) THEN
              ZFRZ=Z(I,J,L)
              EXIT vloop2
            ENDIF
          ENDDO vloop2
!
          LCTOP=NINT(HTOP(I,J))
          ZCTOP=Z(I,J,LCTOP)
          Z1KM=ZINTSFC(I,J)+1000.
          FCTR=0.
vloop3:   IF (ZCTOP >= Z1KM) THEN
            DELZ=Z1KM-ZFRZ
            IF (DELZ <= 0.) THEN
              FCTR=1.        !-- Below the highest freezing level
            ELSE
!
!--- Reduce convective radar reflectivity above freezing level
!
              CUREFL_I=-2./MAX(1000.,ZCTOP-ZFRZ)
              FCTR=10.**(CUREFL_I*DELZ)
            ENDIF
          ENDIF  vloop3
          CUREFL=FCTR*CUREFL
        ENDIF culoop
!
         DO LL=1,2
          IF (C1D(LL) .GE. 1.e-12 .OR. CUREFL .GT. 0.) then
           CALL CALMICT(P1D(LL),T1D(LL),Q1D(LL),C1D(LL), &
                        FI1D(LL),FR1D(LL),FS1D(LL),CUREFL, &
                        DBZ1(LL), I, J, Ilook, Jlook)
          ELSE
           DBZ1(LL)=-20.
          ENDIF
         ENDDO
         FACT=(1000.+ZINTSFC(I,J)-ZMIDloc)/(ZMIDloc-ZMIDP1)
         DBZ1avg=DBZ1(2)+(DBZ1(2)-DBZ1(1))*FACT
         REFDMAX(I,J)=max(REFDMAX(I,J),DBZ1avg)
       ENDDO
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO L=1,LM
       DO J=JMS,JME
        DO I=IMS,IME
         IF (L >= UPINDX(I,J)) THEN
           UPVVELMAX(I,J)=max(UPVVELMAX(I,J),W(I,J,L))
           DNVVELMAX(I,J)=min(DNVVELMAX(I,J),W(I,J,L))
         ENDIF
        ENDDO
       ENDDO
      ENDDO
!
      DO J=JMS,JME
      DO I=IMS,IME
        TLMAX(I,J)=MAX(TLMAX(I,J),T(I,J,LM))  !<--- Hourly max lowest layer T
        TLMIN(I,J)=MIN(TLMIN(I,J),T(I,J,LM))  !<--- Hourly min lowest layer T
        IF (NTSD > 0) THEN
          CAPPA_MOIST=RCP*(1.+QSHLTR(I,J)*R_FACTOR)/(1.+QSHLTR(I,J)*CP_FACTOR)
          T02=TSHLTR(I,J)*(P00_INV*PSHLTR(I,J))**CAPPA_MOIST
          T02MAX(I,J)=MAX(T02MAX(I,J),T02)  !<--- Hourly max 2m T
          T02MIN(I,J)=MIN(T02MIN(I,J),T02)  !<--- Hourly min 2m T
!
          VAPOR_PRESS=PSHLTR(I,J)*QSHLTR(I,J)/                          &
                     (EPSILON+QSHLTR(I,J)*ONE_MINUS_EPSILON)
!- FPVS0 - saturation w/r/t liquid water at all temperatures for RH w/r/t water
          SAT_VAPOR_PRESS=1.E3*FPVS0(T02)
          RH02=MIN(VAPOR_PRESS/SAT_VAPOR_PRESS,0.99)
!
          RH02MAX(I,J)=MAX(RH02MAX(I,J),RH02)     !<--- Hourly max shelter RH
          RH02MIN(I,J)=MIN(RH02MIN(I,J),RH02)     !<--- Hourly min shelter RH
!
          MAGW2=(U10(I,J)**2.+V10(I,J)**2.)
          IF (MAGW2 .gt. SPD10MAX(I,J)) THEN
            U10MAX(I,J)=U10(I,J)                 !<--- U assoc with Hrly max 10m wind speed
            V10MAX(I,J)=V10(I,J)                 !<--- V assoc with Hrly max 10m wind speed
            SPD10MAX(I,J)=MAGW2
          ENDIF
        ENDIF
      ENDDO
      ENDDO

      DO J=JMS,JME
       DO I=IMS,IME
        NCOUNT(I,J)=NCOUNT(I,J)+1
        TERM=-0.273133/T2(I,J)
        P10(I,J)=PSHLTR(I,J)*exp(TERM)
        T10(I,J)=TH10(I,J)*(P10(I,J)/1.e5)**RCP
        T10AVG(I,J)=T10AVG(I,J)*(NCOUNT(I,J)-1)+T10(I,J)
        T10AVG(I,J)=T10AVG(I,J)/NCOUNT(I,J)
        PSFCAVG(I,J)=PSFCAVG(I,J)*(NCOUNT(I,J)-1)+PINT(I,J,LM+1)
        PSFCAVG(I,J)=PSFCAVG(I,J)/NCOUNT(I,J)
        AKHSAVG(I,J)=AKHSAVG(I,J)*(NCOUNT(I,J)-1)+AKHS(I,J)
        AKHSAVG(I,J)=AKHSAVG(I,J)/NCOUNT(I,J)
        AKMSAVG(I,J)=AKMSAVG(I,J)*(NCOUNT(I,J)-1)+AKMS(I,J)
        AKMSAVG(I,J)=AKMSAVG(I,J)/NCOUNT(I,J)
        IF (SNO(I,J) > 0.) THEN
         SNOAVG(I,J)=SNOAVG(I,J)*(NCOUNT(I,J)-1)+1
        ELSE
         SNOAVG(I,J)=SNOAVG(I,J)*(NCOUNT(I,J)-1)
        ENDIF
        SNOAVG(I,J)=SNOAVG(I,J)/NCOUNT(I,J)
       ENDDO
      ENDDO

!-- Maximum precipitation rate (total, frozen)

      DO J=JMS,JME
       DO I=IMS,IME
        PRATEMAX(I,J)=MAX(PRATEMAX(I,J),RDTPHS*PREC(I,J) )
        FPRATEMAX(I,J)=MAX(FPRATEMAX(I,J),RDTPHS*SR(I,J)*PREC(I,J) )
       ENDDO
      ENDDO

      END SUBROUTINE MAX_FIELDS
!
!----------------------------------------------------------------------
!
      SUBROUTINE MAX_FIELDS_HR(T,Q,U                  &
                           ,V,CW                      &
                           ,F_RAIN,F_ICE              &
                           ,F_RIMEF                   &
                           ,Z_DAM,W,REFL_10CM         &
                           ,PINT,PREC                 &
                           ,CPRATE,HTOP               &
                           ,T2,U10,V10                &
                           ,PSHLTR,TSHLTR,QSHLTR      &
                           ,PMID                      &
                           ,REFDMAX,PRATEMAX          &
                           ,FPRATEMAX,SR              &
                           ,UPVVELMAX,DNVVELMAX       &
                           ,TLMAX,TLMIN               &
                           ,T02MAX,T02MIN             &
                           ,RH02MAX,RH02MIN           &
                           ,U10MAX,V10MAX,TH10,T10    &
                           ,SPD10MAX,T10AVG,PSFCAVG   &
                           ,AKHS,AKMS                 &
                           ,AKHSAVG,AKMSAVG           &
                           ,SNO,SNOAVG                &
                           ,UPHLMAX                   &
                           ,DT,NPHS,NTSD              &
                           ,FIS                       &
                           ,IMS,IME,JMS,JME           &
                           ,LM,NCOUNT)

      USE MODULE_MP_FER_HIRES, ONLY : FERRIER_INIT_HR, FPVS0

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: IMS,IME,JMS,JME,LM,NTSD
      INTEGER,INTENT(IN) :: NPHS

      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: NCOUNT

      REAL, DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: T, Q, U, V, CW   &
                                              ,F_RAIN,F_ICE,F_RIMEF        &
                                              ,W,Z_DAM,REFL_10CM

      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM+1),INTENT(IN) :: PINT

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: PREC,CPRATE,HTOP       &
                                                    ,T2,U10,V10            &
                                                    ,PSHLTR,TSHLTR,QSHLTR  &
                                                    ,TH10,AKHS             &
                                                    ,AKMS,SNO,FIS,SR

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT):: REFDMAX,PRATEMAX      &
                                                   ,FPRATEMAX               &
                                                   ,UPVVELMAX,DNVVELMAX     &
                                                   ,TLMAX,TLMIN             &
                                                   ,T02MAX,T02MIN           &
                                                   ,RH02MAX,RH02MIN         &
                                                   ,U10MAX,V10MAX           &
                                                   ,SPD10MAX,T10AVG,PSFCAVG &
                                                   ,UPHLMAX,T10,AKHSAVG     &
                                                   ,AKMSAVG,SNOAVG

      REAL, DIMENSION(IMS:IME,LM),INTENT(IN) :: PMID

      REAL, INTENT(IN) :: DT

      INTEGER :: UPINDX(IMS:IME,JMS:JME)
      REAL, DIMENSION(IMS:IME,JMS:JME) :: P10

      REAL, DIMENSION(IMS:IME,JMS:JME) :: ZINTSFC
      REAL, DIMENSION(IMS:IME,JMS:JME,LM) :: Z

      REAL :: PLOW, PUP,WGTa,WGTb,ZMIDloc,ZMIDP1
      REAL :: P1Da,P1Db,P1D(2)
      REAL :: T1Da,T1Db,T1D(2),fact
      REAL :: Q1Da,Q1Db,Q1D(2)
      REAL :: C1Da,C1Db,C1D(2)
      REAL :: FR1Da,FR1Db,FR1D(2)
      REAL :: FI1Da,FI1Db,FI1D(2)
      REAL :: FS1Da,FS1Db,FS1D(2),DBZ1(2),REFL

      REAL :: CUPRATE,CUREFL,CUREFL_I,ZFRZ,DBZ1avg,FCTR,DELZ,Z1KM,ZCTOP
      REAL :: T02, RH02, TERM
      REAL,SAVE :: CAPPA_MOIST, VAPOR_PRESS, SAT_VAPOR_PRESS
      REAL, SAVE:: DTPHS, RDTPHS
      REAL :: MAGW2

      INTEGER :: LCTOP
      INTEGER :: I,J,L,LL, RC, Ilook,Jlook


!***  COMPUTE AND SAVE THE FACTORS IN R AND CP TO ACCOUNT FOR
!***  WATER VAPOR IN THE AIR.
!***
!***  RECALL: R  = Rd * (1. + Q * (1./EPSILON - 1.))
!***          CP = CPd * (1. + Q * (CPv/CPd - 1.))

      Ilook=99
      Jlook=275

      DTPHS=DT*NPHS
      RDTPHS=3.6e6/DTPHS
!
      DO J=JMS,JME
       DO I=IMS,IME
         ZINTSFC(I,J)=FIS(I,J)/g
         DO L=1,LM
           Z(I,J,L)=Z_DAM(I,J,L)/g !rv FV3 is in gpdam, not gpm
         ENDDO
       ENDDO
      ENDDO
!
!     WON'T BOTHER TO REBUILD HEIGHTS AS IS DONE IN POST.
!     THE NONHYDROSTATIC MID-LAYER Z VALUES MATCH CLOSELY ENOUGH
!     AT 1000 m AGL
!
      DO J=JMS,JME
       DO I=IMS,IME
 L_LOOP: DO L=1,LM-1
          PLOW= PMID(I,L+1)
          PUP=  PMID(I,L)
          IF (PLOW .ge. 40000. .and. PUP .le. 40000.) THEN
            UPINDX(I,J)=L
            exit L_LOOP
          ENDIF
         ENDDO L_LOOP
       ENDDO
      ENDDO
!
      DO J=JMS,JME
       DO I=IMS,IME
  vloop: DO L=8,LM-1
          IF ( (Z(I,J,L+1)-ZINTSFC(I,J)) .LE. 1000.                &
          .AND.(Z(I,J,L)-ZINTSFC(I,J))   .GE. 1000.)  THEN
            ZMIDP1=Z(I,J,L)
            ZMIDloc=Z(I,J,L+1)
            P1D(1)=PMID(I,L)
            P1D(2)=PMID(I,L+1)
            T1D(1)=T(I,J,L)
            T1D(2)=T(I,J,L+1)
            Q1D(1)=Q(I,J,L)
            Q1D(2)=Q(I,J,L+1)
            DBZ1(1)=REFL_10CM(I,J,L)   !- dBZ (not Z) values
            DBZ1(2)=REFL_10CM(I,J,L+1) !- dBZ values
            EXIT vloop
          ENDIF
        ENDDO vloop
!
!!! INITIAL CUREFL VALUE WITHOUT REDUCTION ABOVE FREEZING LEVEL
!
         CUREFL=0.
         IF (CPRATE(I,J)>0.) THEN
           CUPRATE=RDTPHS*CPRATE(I,J)
           CUREFL=CU_A*CUPRATE**CU_B
         ENDIF
!
!-- Ignore convective vertical profile effects when the freezing
!   level is below 1000 m AGL, approximate using the surface value
!
         DO LL=1,2
           REFL=0.
           IF (DBZ1(LL)>DBZmin) REFL=10.**(0.1*DBZ1(LL))
           DBZ1(LL)=CUREFL+REFL    !- in Z units
         ENDDO
!-- Vertical interpolation of Z (units of mm**6/m**3)
         FACT=(1000.+ZINTSFC(I,J)-ZMIDloc)/(ZMIDloc-ZMIDP1)
         DBZ1avg=DBZ1(2)+(DBZ1(2)-DBZ1(1))*FACT
!-- Convert to dBZ (10*logZ) as the last step
         IF (DBZ1avg>ZMIN) THEN
           DBZ1avg=10.*ALOG10(DBZ1avg)
         ELSE
           DBZ1avg=DBZmin
         ENDIF
         REFDMAX(I,J)=max(REFDMAX(I,J),DBZ1avg)
       ENDDO
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO L=1,LM
       DO J=JMS,JME
        DO I=IMS,IME
         IF (L >= UPINDX(I,J)) THEN
           UPVVELMAX(I,J)=max(UPVVELMAX(I,J),W(I,J,L))
           DNVVELMAX(I,J)=min(DNVVELMAX(I,J),W(I,J,L))
         ENDIF
        ENDDO
       ENDDO
      ENDDO
!
      DO J=JMS,JME
      DO I=IMS,IME
        TLMAX(I,J)=MAX(TLMAX(I,J),T(I,J,LM))  !<--- Hourly max lowest layer T
        TLMIN(I,J)=MIN(TLMIN(I,J),T(I,J,LM))  !<--- Hourly min lowest layer T
        IF (NTSD > 0) THEN
          CAPPA_MOIST=RCP*(1.+QSHLTR(I,J)*R_FACTOR)/(1.+QSHLTR(I,J)*CP_FACTOR)
          T02=TSHLTR(I,J)*(P00_INV*PSHLTR(I,J))**CAPPA_MOIST
          T02MAX(I,J)=MAX(T02MAX(I,J),T02)  !<--- Hourly max 2m T
          T02MIN(I,J)=MIN(T02MIN(I,J),T02)  !<--- Hourly min 2m T
!
          VAPOR_PRESS=PSHLTR(I,J)*QSHLTR(I,J)/                          &
                     (EPSILON+QSHLTR(I,J)*ONE_MINUS_EPSILON)
!- FPVS0 - saturation w/r/t liquid water at all temperatures
          SAT_VAPOR_PRESS=1.E3*FPVS0(T02)
          RH02=MIN(VAPOR_PRESS/SAT_VAPOR_PRESS,0.99)
!
          RH02MAX(I,J)=MAX(RH02MAX(I,J),RH02)     !<--- Hourly max shelter RH
          RH02MIN(I,J)=MIN(RH02MIN(I,J),RH02)     !<--- Hourly min shelter RH
!
          MAGW2=(U10(I,J)**2.+V10(I,J)**2.)
          IF (MAGW2 .gt. SPD10MAX(I,J)) THEN
            U10MAX(I,J)=U10(I,J)                 !<--- U assoc with Hrly max 10m wind speed
            V10MAX(I,J)=V10(I,J)                 !<--- V assoc with Hrly max 10m wind speed
            SPD10MAX(I,J)=MAGW2
          ENDIF
        ENDIF
      ENDDO
      ENDDO

      DO J=JMS,JME
       DO I=IMS,IME
        NCOUNT(I,J)=NCOUNT(I,J)+1
        TERM=-0.273133/T2(I,J)
        P10(I,J)=PSHLTR(I,J)*exp(TERM)
        T10(I,J)=TH10(I,J)*(P10(I,J)/1.e5)**RCP
        T10AVG(I,J)=T10AVG(I,J)*(NCOUNT(I,J)-1)+T10(I,J)
        T10AVG(I,J)=T10AVG(I,J)/NCOUNT(I,J)
        PSFCAVG(I,J)=PSFCAVG(I,J)*(NCOUNT(I,J)-1)+PINT(I,J,LM+1)
        PSFCAVG(I,J)=PSFCAVG(I,J)/NCOUNT(I,J)
        AKHSAVG(I,J)=AKHSAVG(I,J)*(NCOUNT(I,J)-1)+AKHS(I,J)
        AKHSAVG(I,J)=AKHSAVG(I,J)/NCOUNT(I,J)
        AKMSAVG(I,J)=AKMSAVG(I,J)*(NCOUNT(I,J)-1)+AKMS(I,J)
        AKMSAVG(I,J)=AKMSAVG(I,J)/NCOUNT(I,J)
        IF (SNO(I,J) > 0.) THEN
         SNOAVG(I,J)=SNOAVG(I,J)*(NCOUNT(I,J)-1)+1
        ELSE
         SNOAVG(I,J)=SNOAVG(I,J)*(NCOUNT(I,J)-1)
        ENDIF
        SNOAVG(I,J)=SNOAVG(I,J)/NCOUNT(I,J)
       ENDDO
      ENDDO

!-- Maximum precipitation rate (total, frozen)

      DO J=JMS,JME
       DO I=IMS,IME
        PRATEMAX(I,J)=MAX(PRATEMAX(I,J),RDTPHS*PREC(I,J) )
        FPRATEMAX(I,J)=MAX(FPRATEMAX(I,J),RDTPHS*SR(I,J)*PREC(I,J) )
       ENDDO
      ENDDO

      END SUBROUTINE MAX_FIELDS_HR
!
!----------------------------------------------------------------------
!

      SUBROUTINE MAX_FIELDS_w6(T,Q,U,V,Z_DAM,W        &
                           ,QR,QS,QG,PINT,PREC        &
                           ,CPRATE,HTOP               &
                           ,T2,U10,V10                &
                           ,PSHLTR,TSHLTR,QSHLTR      &
                           ,PMID                      &
                           ,REFDMAX,PRATEMAX          &
                           ,FPRATEMAX,SR              &
                           ,UPVVELMAX,DNVVELMAX       &
                           ,TLMAX,TLMIN               &
                           ,T02MAX,T02MIN             &
                           ,RH02MAX,RH02MIN           &
                           ,U10MAX,V10MAX,TH10,T10    &
                           ,SPD10MAX,T10AVG,PSFCAVG   &
                           ,AKHS,AKMS                 &
                           ,AKHSAVG,AKMSAVG           &
                           ,SNO,SNOAVG                &
                           ,UPHLMAX                   &
                           ,DT,NPHS,NTSD              &
                           ,FIS                       &
                           ,IMS,IME,JMS,JME           &
                           ,LM,NCOUNT)

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: IMS,IME,JMS,JME,LM,NTSD
      INTEGER,INTENT(IN) :: NPHS

      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: NCOUNT

      REAL, DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: T,Q,U,V,Z_DAM,W,QR,QS,QG

      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM+1),INTENT(IN) :: PINT

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: PREC,CPRATE,HTOP       &
                                                    ,T2,U10,V10            &
                                                    ,PSHLTR,TSHLTR,QSHLTR  &
                                                    ,TH10,AKHS             &
                                                    ,AKMS,SNO,FIS,SR

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT):: REFDMAX,PRATEMAX      &
                                                   ,FPRATEMAX               &
                                                   ,UPVVELMAX,DNVVELMAX     &
                                                   ,TLMAX,TLMIN             &
                                                   ,T02MAX,T02MIN           &
                                                   ,RH02MAX,RH02MIN         &
                                                   ,U10MAX,V10MAX           &
                                                   ,SPD10MAX,T10AVG,PSFCAVG &
                                                   ,UPHLMAX,T10,AKHSAVG     &
                                                   ,AKMSAVG,SNOAVG

      REAL, DIMENSION(IMS:IME,LM),INTENT(IN) :: PMID

      REAL, INTENT(IN) :: DT

      INTEGER :: UPINDX(IMS:IME,JMS:JME)
      REAL, DIMENSION(IMS:IME,JMS:JME) :: P10

      REAL, DIMENSION(IMS:IME,JMS:JME) :: ZINTSFC
      REAL, DIMENSION(IMS:IME,JMS:JME,LM) :: Z

      REAL :: PLOW, PUP,WGTa,WGTb,ZMIDloc,ZMIDP1
      REAL :: P1Da,P1Db,P1D(2)
      REAL :: T1Da,T1Db,T1D(2),fact
      REAL :: Q1Da,Q1Db,Q1D(2)
      REAL :: QQR(2),QQS(2),QQG(2),QPCP,DENS,N0S
      REAL :: DBZR,DBZS,DBZG,DBZ1(2)
      REAL, PARAMETER :: N0S0=2.E6,N0Smax=1.E11,ALPHA=0.12   &
             ,N0G=4.E6,RHOS=100.,RHOG=500.,ZRADR=3.631E9     &
             ,DBZmin=-20.

      REAL :: CUPRATE, CUREFL, CUREFL_I, ZFRZ, DBZ1avg, FCTR, DELZ
      REAL :: T02, RH02, TERM, TREF
      REAL,SAVE :: CAPPA_MOIST, VAPOR_PRESS, SAT_VAPOR_PRESS
      REAL, SAVE:: DTPHS, RDTPHS, ZRADS,ZRADG,ZMIN
      REAL :: MAGW2

      INTEGER :: LCTOP
      INTEGER :: I,J,L,LL, RC, Ilook,Jlook


!***  COMPUTE AND SAVE THE FACTORS IN R AND CP TO ACCOUNT FOR
!***  WATER VAPOR IN THE AIR.
!***
!***  RECALL: R  = Rd * (1. + Q * (1./EPSILON - 1.))
!***          CP = CPd * (1. + Q * (CPv/CPd - 1.))

      Ilook=99
      Jlook=275

      DTPHS=DT*NPHS
      RDTPHS=3.6e6/DTPHS
!-- For calculating radar reflectivity
      ZRADS=2.17555E13*RHOS**0.25
      ZRADG=2.17555E13*RHOG**0.25/N0G**0.75
      ZMIN=10.**(0.1*DBZmin)
!
      DO J=JMS,JME
       DO I=IMS,IME
         ZINTSFC(I,J)=FIS(I,J)/g
         DO L=1,LM
           Z(I,J,L)=Z_DAM(I,J,L)/g !rv FV3 is in gpdam, not gpm
         ENDDO
       ENDDO
      ENDDO
!
!     WON'T BOTHER TO REBUILD HEIGHTS AS IS DONE IN POST.
!     THE NONHYDROSTATIC MID-LAYER Z VALUES MATCH CLOSELY ENOUGH
!     AT 1000 m AGL
!
      DO J=JMS,JME
       DO I=IMS,IME
 L_LOOP: DO L=1,LM-1
          PLOW= PMID(I,L+1)
          PUP=  PMID(I,L)
          IF (PLOW .ge. 40000. .and. PUP .le. 40000.) THEN
            UPINDX(I,J)=L
            exit L_LOOP
          ENDIF
         ENDDO L_LOOP
       ENDDO
      ENDDO
!
      DO J=JMS,JME
       DO I=IMS,IME
  vloop: DO L=8,LM-1
          IF ( (Z(I,J,L+1)-ZINTSFC(I,J)) .LE. 1000.                &
          .AND.(Z(I,J,L)-ZINTSFC(I,J))   .GE. 1000.)  THEN
            ZMIDP1=Z(I,J,L)
            ZMIDloc=Z(I,J,L+1)
            P1D(1)=PMID(I,L)
            P1D(2)=PMID(I,L+1)
            T1D(1)=T(I,J,L)
            T1D(2)=T(I,J,L+1)
            Q1D(1)=Q(I,J,L)
            Q1D(2)=Q(I,J,L+1)
            QQR(1)=QR(I,J,L)
            QQR(2)=QR(I,J,L+1)
            QQS(1)=QS(I,J,L)
            QQS(2)=QS(I,J,L+1)
            QQG(1)=QG(I,J,L)
            QQG(2)=QG(I,J,L+1)
            EXIT vloop
          ENDIF
        ENDDO vloop
!
!!! INITIAL CUREFL VALUE WITHOUT REDUCTION ABOVE FREEZING LEVEL
!
        CUPRATE=RDTPHS*CPRATE(I,J)
        CUREFL=0.
        IF (CUPRATE>0.) CUREFL=CU_A*CUPRATE**CU_B
!
!-- Ignore convective vertical profile effects when the freezing
!   level is below 1000 m AGL, approximate using the surface value
!
         DO LL=1,2
           DBZ1(LL)=CUREFL
           QPCP=QQR(LL)+QQS(LL)+QQG(LL)
!-- A higher threshold can be used for calculating radar reflectivities
!   above DBZmin=-20 dBZ; note the DBZ arrays below are actually in
!   Z units of mm**6/m**3
           IF (QPCP>1.E-8) THEN
              DBZR=0.
              DBZS=0.
              DBZG=0.
              DENS=P1D(LL)/(R_D*T1D(LL)*(Q1D(LL)*P608+1.0))
              IF(QQR(LL)>1.E-8) DBZR=ZRADR*((QQR(LL)*DENS)**1.75)
              IF(QQS(LL)>1.E-8) THEN
                 N0S=N0S0*MAX(1., EXP(ALPHA*(TIW-T1D(LL) ) ) )
                 N0S=MIN(N0S, N0Smax)
                 DBZS=ZRADS*((QQS(LL)*DENS)**1.75)/N0S**0.75
              ENDIF
              IF(QQG(LL)>1.E-8) DBZG=ZRADG*((QQG(LL)*DENS)**1.75)
              DBZ1(LL)=DBZ1(LL)+DBZR+DBZS+DBZG
           ENDIF
         ENDDO
!-- Vertical interpolation of Z (units of mm**6/m**3)
         FACT=(1000.+ZINTSFC(I,J)-ZMIDloc)/(ZMIDloc-ZMIDP1)
         DBZ1avg=DBZ1(2)+(DBZ1(2)-DBZ1(1))*FACT
!-- Convert to dBZ (10*logZ) as the last step
         IF (DBZ1avg>ZMIN) THEN
            DBZ1avg=10.*ALOG10(DBZ1avg)
         ELSE
            DBZ1avg=DBZmin
         ENDIF
         REFDMAX(I,J)=max(REFDMAX(I,J),DBZ1avg)
       ENDDO
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO L=1,LM
       DO J=JMS,JME
        DO I=IMS,IME
         IF (L >= UPINDX(I,J)) THEN
           UPVVELMAX(I,J)=max(UPVVELMAX(I,J),W(I,J,L))
           DNVVELMAX(I,J)=min(DNVVELMAX(I,J),W(I,J,L))
         ENDIF
        ENDDO
       ENDDO
      ENDDO
!
      DO J=JMS,JME
      DO I=IMS,IME
        TLMAX(I,J)=MAX(TLMAX(I,J),T(I,J,LM))  !<--- Hourly max lowest layer T
        TLMIN(I,J)=MIN(TLMIN(I,J),T(I,J,LM))  !<--- Hourly min lowest layer T
        IF (NTSD > 0) THEN
          CAPPA_MOIST=RCP*(1.+QSHLTR(I,J)*R_FACTOR)/(1.+QSHLTR(I,J)*CP_FACTOR)
          T02=TSHLTR(I,J)*(P00_INV*PSHLTR(I,J))**CAPPA_MOIST
          T02MAX(I,J)=MAX(T02MAX(I,J),T02)  !<--- Hourly max 2m T
          T02MIN(I,J)=MIN(T02MIN(I,J),T02)  !<--- Hourly min 2m T
!
          VAPOR_PRESS=PSHLTR(I,J)*QSHLTR(I,J)/                          &
                     (EPSILON+QSHLTR(I,J)*ONE_MINUS_EPSILON)

!-- Adapted from WSM6 code:
          TREF=TTP/T02
          SAT_VAPOR_PRESS=PSAT*EXP(LOG(TREF)*(XA))*EXP(XB*(1.-TREF))

          RH02=MIN(VAPOR_PRESS/SAT_VAPOR_PRESS,0.99)
!
          RH02MAX(I,J)=MAX(RH02MAX(I,J),RH02)     !<--- Hourly max shelter RH
          RH02MIN(I,J)=MIN(RH02MIN(I,J),RH02)     !<--- Hourly min shelter RH
!
          MAGW2=(U10(I,J)**2.+V10(I,J)**2.)
          IF (MAGW2 .gt. SPD10MAX(I,J)) THEN
            U10MAX(I,J)=U10(I,J)                 !<--- U assoc with Hrly max 10m wind speed
            V10MAX(I,J)=V10(I,J)                 !<--- V assoc with Hrly max 10m wind speed
            SPD10MAX(I,J)=MAGW2
          ENDIF
        ENDIF
      ENDDO
      ENDDO

      DO J=JMS,JME
       DO I=IMS,IME
        NCOUNT(I,J)=NCOUNT(I,J)+1
        TERM=-0.273133/T2(I,J)
        P10(I,J)=PSHLTR(I,J)*exp(TERM)
        T10(I,J)=TH10(I,J)*(P10(I,J)/1.e5)**RCP
        T10AVG(I,J)=T10AVG(I,J)*(NCOUNT(I,J)-1)+T10(I,J)
        T10AVG(I,J)=T10AVG(I,J)/NCOUNT(I,J)
        PSFCAVG(I,J)=PSFCAVG(I,J)*(NCOUNT(I,J)-1)+PINT(I,J,LM+1)
        PSFCAVG(I,J)=PSFCAVG(I,J)/NCOUNT(I,J)
        AKHSAVG(I,J)=AKHSAVG(I,J)*(NCOUNT(I,J)-1)+AKHS(I,J)
        AKHSAVG(I,J)=AKHSAVG(I,J)/NCOUNT(I,J)
        AKMSAVG(I,J)=AKMSAVG(I,J)*(NCOUNT(I,J)-1)+AKMS(I,J)
        AKMSAVG(I,J)=AKMSAVG(I,J)/NCOUNT(I,J)
        IF (SNO(I,J) > 0.) THEN
         SNOAVG(I,J)=SNOAVG(I,J)*(NCOUNT(I,J)-1)+1
        ELSE
         SNOAVG(I,J)=SNOAVG(I,J)*(NCOUNT(I,J)-1)
        ENDIF
        SNOAVG(I,J)=SNOAVG(I,J)/NCOUNT(I,J)
       ENDDO
      ENDDO

!-- Maximum precipitation rate (total, frozen)

      DO J=JMS,JME
       DO I=IMS,IME
        PRATEMAX(I,J)=MAX(PRATEMAX(I,J),RDTPHS*PREC(I,J) )
        FPRATEMAX(I,J)=MAX(FPRATEMAX(I,J),RDTPHS*SR(I,J)*PREC(I,J) )
       ENDDO
      ENDDO

      END SUBROUTINE MAX_FIELDS_W6
!
!----------------------------------------------------------------------
!
      SUBROUTINE MAX_FIELDS_THO(T,Q,U,V,Z_DAM,W       &
                           ,REFL_10CM,PINT,PREC       &
                           ,CPRATE,HTOP               &
                           ,T2,U10,V10                &
                           ,PSHLTR,TSHLTR,QSHLTR      &
                           ,PMID                      &
                           ,REFDMAX,PRATEMAX          &
                           ,FPRATEMAX,SR              &
                           ,UPVVELMAX,DNVVELMAX       &
                           ,TLMAX,TLMIN               &
                           ,T02MAX,T02MIN             &
                           ,RH02MAX,RH02MIN           &
                           ,U10MAX,V10MAX,TH10,T10    &
                           ,SPD10MAX,T10AVG,PSFCAVG   &
                           ,AKHS,AKMS                 &
                           ,AKHSAVG,AKMSAVG           &
                           ,SNO,SNOAVG                &
                           ,UPHLMAX                   &
                           ,DT,NPHS,NTSD              &
                           ,FIS                       &
                           ,IMS,IME,JMS,JME           &
                           ,LM,NCOUNT)

      USE MODULE_MP_THOMPSON, ONLY : RSLF

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: IMS,IME,JMS,JME,LM,NTSD
      INTEGER,INTENT(IN) :: NPHS

      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: NCOUNT

      REAL, DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: T,Q,U,V,Z_DAM,W

      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: REFL_10CM

      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM+1),INTENT(IN) :: PINT

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: PREC,CPRATE,HTOP       &
                                                    ,T2,U10,V10            &
                                                    ,PSHLTR,TSHLTR,QSHLTR  &
                                                    ,TH10,AKHS             &
                                                    ,AKMS,SNO,FIS,SR

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT):: REFDMAX,PRATEMAX      &
                                                   ,FPRATEMAX               &
                                                   ,UPVVELMAX,DNVVELMAX     &
                                                   ,TLMAX,TLMIN             &
                                                   ,T02MAX,T02MIN           &
                                                   ,RH02MAX,RH02MIN         &
                                                   ,U10MAX,V10MAX           &
                                                   ,SPD10MAX,T10AVG,PSFCAVG &
                                                   ,UPHLMAX,T10,AKHSAVG     &
                                                   ,AKMSAVG,SNOAVG

      REAL, DIMENSION(IMS:IME,LM),INTENT(IN) :: PMID

      REAL, INTENT(IN) :: DT

      INTEGER :: UPINDX(IMS:IME,JMS:JME)
      REAL, DIMENSION(IMS:IME,JMS:JME) :: P10

      REAL, DIMENSION(IMS:IME,JMS:JME) :: ZINTSFC
      REAL, DIMENSION(IMS:IME,JMS:JME,LM) :: Z

      REAL :: PLOW, PUP,WGTa,WGTb,ZMIDloc,ZMIDP1
      REAL :: P1Da,P1Db,P1D(2)
      REAL :: T1Da,T1Db,T1D(2),fact
      REAL :: Q1Da,Q1Db,Q1D(2)
      REAL :: QQR(2),QQS(2),QQG(2),QPCP,DENS,N0S
      REAL :: REFL,DBZ1(2)
      REAL, PARAMETER :: N0S0=2.E6,N0Smax=1.E11,ALPHA=0.12   &
             ,N0G=4.E6,RHOS=100.,RHOG=500.,ZRADR=3.631E9     &
             ,DBZmin=-20.

      REAL :: CUPRATE, CUREFL, CUREFL_I, ZFRZ, DBZ1avg, FCTR, DELZ
      REAL :: T02, RH02, TERM, TREF
      REAL,SAVE :: CAPPA_MOIST, QVSHLTR, QVSAT
      REAL, SAVE:: DTPHS, RDTPHS, ZRADS,ZRADG,ZMIN
      REAL :: MAGW2

      INTEGER :: LCTOP
      INTEGER :: I,J,L,LL, RC, Ilook,Jlook


!***  COMPUTE AND SAVE THE FACTORS IN R AND CP TO ACCOUNT FOR
!***  WATER VAPOR IN THE AIR.
!***
!***  RECALL: R  = Rd * (1. + Q * (1./EPSILON - 1.))
!***          CP = CPd * (1. + Q * (CPv/CPd - 1.))

      Ilook=99
      Jlook=275

      DTPHS=DT*NPHS
      RDTPHS=3.6e6/DTPHS
!-- For calculating radar reflectivity
      ZMIN=10.**(0.1*DBZmin)
!
      DO J=JMS,JME
       DO I=IMS,IME
         ZINTSFC(I,J)=FIS(I,J)/g
         DO L=1,LM
           Z(I,J,L)=Z_DAM(I,J,L)/g !rv FV3 is in gpdam, not gpm
         ENDDO
       ENDDO
      ENDDO
!
!     WON'T BOTHER TO REBUILD HEIGHTS AS IS DONE IN POST.
!     THE NONHYDROSTATIC MID-LAYER Z VALUES MATCH CLOSELY ENOUGH
!     AT 1000 m AGL
!
      DO J=JMS,JME
       DO I=IMS,IME
 L_LOOP: DO L=1,LM-1
          PLOW= PMID(I,L+1)
          PUP=  PMID(I,L)
          IF (PLOW .ge. 40000. .and. PUP .le. 40000.) THEN
            UPINDX(I,J)=L
            exit L_LOOP
          ENDIF
         ENDDO L_LOOP
       ENDDO
      ENDDO
!
      DO J=JMS,JME
       DO I=IMS,IME
  vloop: DO L=8,LM-1
          IF ( (Z(I,J,L+1)-ZINTSFC(I,J)) .LE. 1000.                &
          .AND.(Z(I,J,L)-ZINTSFC(I,J))   .GE. 1000.)  THEN
            ZMIDP1=Z(I,J,L)
            ZMIDloc=Z(I,J,L+1)
            P1D(1)=PMID(I,L)
            P1D(2)=PMID(I,L+1)
            T1D(1)=T(I,J,L)
            T1D(2)=T(I,J,L+1)
            Q1D(1)=Q(I,J,L)
            Q1D(2)=Q(I,J,L+1)
            DBZ1(1)=REFL_10CM(I,J,L)   !- dBZ (not Z) values
            DBZ1(2)=REFL_10CM(I,J,L+1) !- dBZ values
            EXIT vloop
          ENDIF
         ENDDO vloop
!
!!! INITIAL CUREFL VALUE WITHOUT REDUCTION ABOVE FREEZING LEVEL
!
         CUREFL=0.
         IF (CPRATE(I,J)>0.) THEN
           CUPRATE=RDTPHS*CPRATE(I,J)
           CUREFL=CU_A*CUPRATE**CU_B
         ENDIF
!
!-- Ignore convective vertical profile effects when the freezing
!   level is below 1000 m AGL, approximate using the surface value
!
         DO LL=1,2
           REFL=0.
           IF (DBZ1(LL)>DBZmin) REFL=10.**(0.1*DBZ1(LL))
           DBZ1(LL)=CUREFL+REFL    !- in Z units
         ENDDO
!-- Vertical interpolation of Z (units of mm**6/m**3)
         FACT=(1000.+ZINTSFC(I,J)-ZMIDloc)/(ZMIDloc-ZMIDP1)
         DBZ1avg=DBZ1(2)+(DBZ1(2)-DBZ1(1))*FACT
!-- Convert to dBZ (10*logZ) as the last step
         IF (DBZ1avg>ZMIN) THEN
           DBZ1avg=10.*ALOG10(DBZ1avg)
         ELSE
           DBZ1avg=DBZmin
         ENDIF
         REFDMAX(I,J)=max(REFDMAX(I,J),DBZ1avg)
       ENDDO
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO L=1,LM
       DO J=JMS,JME
        DO I=IMS,IME
         IF (L >= UPINDX(I,J)) THEN
           UPVVELMAX(I,J)=max(UPVVELMAX(I,J),W(I,J,L))
           DNVVELMAX(I,J)=min(DNVVELMAX(I,J),W(I,J,L))
         ENDIF
        ENDDO
       ENDDO
      ENDDO
!
      DO J=JMS,JME
      DO I=IMS,IME
        TLMAX(I,J)=MAX(TLMAX(I,J),T(I,J,LM))  !<--- Hourly max lowest layer T
        TLMIN(I,J)=MIN(TLMIN(I,J),T(I,J,LM))  !<--- Hourly min lowest layer T
        IF (NTSD > 0) THEN
          CAPPA_MOIST=RCP*(1.+QSHLTR(I,J)*R_FACTOR)/(1.+QSHLTR(I,J)*CP_FACTOR)
          T02=TSHLTR(I,J)*(P00_INV*PSHLTR(I,J))**CAPPA_MOIST
          T02MAX(I,J)=MAX(T02MAX(I,J),T02)  !<--- Hourly max 2m T
          T02MIN(I,J)=MIN(T02MIN(I,J),T02)  !<--- Hourly min 2m T
!
          QVSHLTR=QSHLTR(I,J)/(1.-QSHLTR(I,J))  !<-- 2-m water vapor mixing ratio
!
!-- Adapted from Thompson code:
!
          TREF=T02-273.15
          QVSAT=RSLF(PSHLTR(I,J),TREF)
!
          RH02=MIN(QVSHLTR/QVSAT,0.99)
!
          RH02MAX(I,J)=MAX(RH02MAX(I,J),RH02)     !<--- Hourly max shelter RH
          RH02MIN(I,J)=MIN(RH02MIN(I,J),RH02)     !<--- Hourly min shelter RH
!
          MAGW2=(U10(I,J)**2.+V10(I,J)**2.)
          IF (MAGW2 .gt. SPD10MAX(I,J)) THEN
            U10MAX(I,J)=U10(I,J)                 !<--- U assoc with Hrly max 10m wind speed
            V10MAX(I,J)=V10(I,J)                 !<--- V assoc with Hrly max 10m wind speed
            SPD10MAX(I,J)=MAGW2
          ENDIF
        ENDIF
      ENDDO
      ENDDO

      DO J=JMS,JME
       DO I=IMS,IME
        NCOUNT(I,J)=NCOUNT(I,J)+1
        TERM=-0.273133/T2(I,J)
        P10(I,J)=PSHLTR(I,J)*exp(TERM)
        T10(I,J)=TH10(I,J)*(P10(I,J)/1.e5)**RCP
        T10AVG(I,J)=T10AVG(I,J)*(NCOUNT(I,J)-1)+T10(I,J)
        T10AVG(I,J)=T10AVG(I,J)/NCOUNT(I,J)
        PSFCAVG(I,J)=PSFCAVG(I,J)*(NCOUNT(I,J)-1)+PINT(I,J,LM+1)
        PSFCAVG(I,J)=PSFCAVG(I,J)/NCOUNT(I,J)
        AKHSAVG(I,J)=AKHSAVG(I,J)*(NCOUNT(I,J)-1)+AKHS(I,J)
        AKHSAVG(I,J)=AKHSAVG(I,J)/NCOUNT(I,J)
        AKMSAVG(I,J)=AKMSAVG(I,J)*(NCOUNT(I,J)-1)+AKMS(I,J)
        AKMSAVG(I,J)=AKMSAVG(I,J)/NCOUNT(I,J)
        IF (SNO(I,J) > 0.) THEN
         SNOAVG(I,J)=SNOAVG(I,J)*(NCOUNT(I,J)-1)+1
        ELSE
         SNOAVG(I,J)=SNOAVG(I,J)*(NCOUNT(I,J)-1)
        ENDIF
        SNOAVG(I,J)=SNOAVG(I,J)/NCOUNT(I,J)
       ENDDO
      ENDDO

!-- Maximum precipitation rate (total, frozen)

      DO J=JMS,JME
       DO I=IMS,IME
        PRATEMAX(I,J)=MAX(PRATEMAX(I,J),RDTPHS*PREC(I,J) )
        FPRATEMAX(I,J)=MAX(FPRATEMAX(I,J),RDTPHS*SR(I,J)*PREC(I,J) )
       ENDDO
      ENDDO

      END SUBROUTINE MAX_FIELDS_THO
!
!----------------------------------------------------------------------
!
      END MODULE MODULE_DIAGNOSE
      SUBROUTINE LA2GA
! stub
      END SUBROUTINE LA2GA
!
!----------------------------------------------------------------------
