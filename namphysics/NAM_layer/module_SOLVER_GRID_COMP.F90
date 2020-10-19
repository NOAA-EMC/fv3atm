!-----------------------------------------------------------------------
!
      MODULE module_SOLVER_GRID_COMP
!
!-----------------------------------------------------------------------
!
      USE MODULE_KINDS
      USE MODULE_DIAGNOSE         ,ONLY : MAX_FIELDS_driver
      USE MODULE_RADIATION        ,ONLY : RADIATION
      USE MODULE_RA_GFDL          ,ONLY : RDTEMP,TIME_MEASURE
      USE MODULE_TURBULENCE
      USE MODULE_CONVECTION
      USE MODULE_MICROPHYSICS_NMM ,ONLY : GSMDRIVE,UPDATE_WATER,TQADJUST

      use GFS_typedefs,      only: GFS_statein_type, GFS_stateout_type, &
                                   GFS_sfcprop_type, GFS_coupling_type, &
                                   GFS_control_type, GFS_grid_type, &
                                   GFS_tbd_type,     GFS_cldprop_type, &
                                   GFS_radtend_type, GFS_diag_type
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: SOLVER_RUN
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE SOLVER_RUN &
         (Model, Statein, Stateout, Sfcprop, Coupling,  &
          Grid, Tbd, Cldprop, Radtend, Diag)
!
      type(GFS_control_type),         intent(in)    :: Model
      type(GFS_statein_type),         intent(inout) :: Statein
      type(GFS_stateout_type),        intent(inout) :: Stateout
      type(GFS_sfcprop_type),         intent(inout) :: Sfcprop
      type(GFS_coupling_type),        intent(inout) :: Coupling
      type(GFS_grid_type),            intent(in)    :: Grid
      type(GFS_tbd_type),             intent(inout) :: Tbd
      type(GFS_cldprop_type),         intent(inout) :: Cldprop
      type(GFS_radtend_type),         intent(inout) :: Radtend
      type(GFS_diag_type),            intent(inout) :: Diag
!
      INTEGER(kind=KINT),DIMENSION(size(Grid%xlon,1))            :: ISLTYP,IVGTYP
      REAL(kind=KFPT)   ,DIMENSION(size(Grid%xlon,1),Model%levs) :: DUDT,DVDT,DTDT
!
      LOGICAL(kind=KLOG) :: LISS_RESTART
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: IMS=1,IME
      INTEGER(kind=KINT) :: I,K,L,NTIMESTEP,NTIMESTEP_RAD,MYPE,LM
!
      INTEGER(kind=KINT) :: JULDAY,JULYR
!
      REAL(kind=KFPT)    :: JULIAN,XTIME,PI
!
      LOGICAL(kind=KLOG) :: CALL_LONGWAVE                               &
                           ,CALL_SHORTWAVE                              &
                           ,CALL_TURBULENCE                             &
                           ,CALL_PRECIP
!
!-----------------------------------------------------------------------
!***  Start here
!-----------------------------------------------------------------------
!
      IME=size(Grid%xlon,1)
      LM=Model%levs
      NTIMESTEP=Model%kdt-1
      MYPE=Model%me
      PI=2.*asin(1.)
!
!-----------------------------------------------------------------------
!***  Call  radiation so that updated fields are written to the
!***  history files after 0 hours.
!-----------------------------------------------------------------------
!
      IF(NTIMESTEP==0)THEN
        NTIMESTEP_RAD=NTIMESTEP
      ELSE
        NTIMESTEP_RAD=NTIMESTEP+1
      ENDIF
!
      CALL TIME_MEASURE(Model%idat(1),Model%idat(2),Model%idat(3)       &
                       ,Model%idat(5),Model%idat(6),Model%idat(7)       &
                       ,NTIMESTEP_rad,Model%dtp                         &
                       ,JULDAY,JULYR,JULIAN,XTIME)
!
!-----------------------------------------------------------------------
!*** Set tendencies to zero
!*** Initialize next step's values
!-----------------------------------------------------------------------
!
        DTDT(:,:)=0.
        DUDT(:,:)=0.
        DVDT(:,:)=0.
        Stateout%gu0 = Statein%ugrs
        Stateout%gv0 = Statein%vgrs
        Stateout%gq0 = Statein%qgrs
!
!-----------------------------------------------------------------------
!***  This part should be in INIT
!-----------------------------------------------------------------------
!
      IF(NTIMESTEP==0) THEN
        do i=ims,ime
          Sfcprop%sst(i)=Sfcprop%tsfc(i)
          Sfcprop%THS(i)=Sfcprop%tsfc(i)
          Radtend%sfalb(i)=(Sfcprop%alvwf(i)+Sfcprop%alnwf(i))*0.5
          Sfcprop%ALBASE(i)=Radtend%sfalb(i)
          Sfcprop%sncovr(I)=0.0; if(Sfcprop%sncovr(I) > 0.0) Sfcprop%sncovr(I)=0.98
          Sfcprop%canopy(i)=Sfcprop%canopy(i)*0.001
          Sfcprop%oro(i)=Sfcprop%oro(i)*9.81
          if(Grid%xlon(i)>PI) Grid%xlon(i)=Grid%xlon(i)-2.*PI
          Sfcprop%sice(i)=int(Sfcprop%slmsk(i)*0.5)
          Sfcprop%sm(i)=1.; if(Sfcprop%slmsk(i) > 0.5 ) Sfcprop%sm(i)=0.
          Sfcprop%zorl(i)=Sfcprop%zorl(i)*0.01
        enddo
      ENDIF
!
      IVGTYP(:)=Sfcprop%vtype(:)
      ISLTYP(:)=Sfcprop%stype(:)
!
!-----------------------------------------------------------------------
!
        IF (mod(NTIMESTEP,Model%NSTEPS_PER_CHECK) == 0 ) THEN
!
          CALL MAX_FIELDS_driver(Statein%tgrs(:,:)                      &
                         ,Stateout%gq0(:,:,1),Statein%ugrs(:,:)         &
                         ,Statein%vgrs(:,:),Stateout%gq0(:,:,2)         & !rv CW vs. QC
                         ,Statein%f_rain(:,:),Statein%f_ice(:,:)        &
                         ,Statein%f_rimef(:,:),Statein%phil(:,:)        &
                         ,Statein%vvl(:,:),Statein%refl_10cm(:,:)       &
                         ,Stateout%gq0(:,:,Model%ntrw)                  &
                         ,Stateout%gq0(:,:,Model%ntsw)                  &
                         ,Stateout%gq0(:,:,Model%ntgl)                  &
                         ,Statein%prsi(:,:),Sfcprop%tprcp(:)            &
                         ,Diag%rainc(:),Diag%htop(:)                    &
                         ,Sfcprop%t2m(:),Diag%u10m(:),Diag%v10m(:)      &
                         ,Diag%pshltr(:),Diag%tshltr(:)                 &
                         ,Diag%qshltr(:),Statein%prsl(:,:)              &
                         ,Diag%REFDMAX(:),Diag%PRATEMAX(:)              &
                         ,Diag%FPRATEMAX(:),Diag%sr(:)                  &
                         ,Diag%UPVVELMAX(:),Diag%DNVVELMAX(:)           &
                         ,Diag%TLMAX(:),Diag%TLMIN(:)                   &
                         ,Diag%T02MAX(:),Diag%T02MIN(:)                 &
                         ,Diag%RH02MAX(:),Diag%RH02MIN(:)               &
                         ,Diag%U10MAX(:),Diag%V10MAX(:)                 &
                         ,Diag%TH10(:),Diag%T10(:)                      &
                         ,Diag%SPD10MAX(:),Diag%t10avg(:)               &
                         ,Diag%psfcavg(:),Diag%chh(:),Diag%cmm(:)       &
                         ,Diag%akhsavg(:),Diag%akmsavg(:)               &
                         ,Sfcprop%weasd(:),Diag%snoavg(:)               &
                         ,Diag%UPHLMAX(:)                               &
                         ,Model%dtp,Model%NPHS,NTIMESTEP                &
                         ,Model%NSTEPS_PER_RESET,Sfcprop%oro(:)         &
                         ,ims,ime                                       &
                         ,LM,Diag%NCOUNT(:),Model%MICROPHYSICS)
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Set logical switches for calling each of the Physics schemes.
!-----------------------------------------------------------------------
!
        CALL_SHORTWAVE  = MOD(NTIMESTEP_RAD,Model%NRADS)==0
        CALL_LONGWAVE   = MOD(NTIMESTEP_RAD,Model%NRADL)==0
        CALL_TURBULENCE = MOD(NTIMESTEP,Model%NPHS)==0
        CALL_PRECIP     = MOD(NTIMESTEP,Model%NPRECIP)==0
!
!-----------------------------------------------------------------------
!***  Update WATER array from CWM, F_ICE, F_RAIN for Ferrier
!***  microphysics but only if any of the Physics subroutines
!***  are called
!***  Expanded to also update CWM, F_ICE, F_RAIN, F_RIMEF for non-Ferrier
!***  microphysics.
!-----------------------------------------------------------------------
!
        IF((Model%MICROPHYSICS=='fer'                                   &
                       .OR.                                             &
            Model%MICROPHYSICS=='fer_hires'                             &
                       .OR.                                             &
            Model%MICROPHYSICS=='gfs'                                   &
                       .OR.                                             &
            Model%MICROPHYSICS=='wsm6'                                  &
                       .OR.                                             &
            Model%MICROPHYSICS=='thompson')                             &
                       .AND.                                            &
           (CALL_SHORTWAVE .OR. CALL_LONGWAVE .OR.                      &
            CALL_TURBULENCE .OR. CALL_PRECIP) ) THEN
!
          CALL UPDATE_WATER(Stateout%gq0(:,:,2),Statein%f_ice(:,:)      &
                           ,Statein%f_rain(:,:),Statein%f_rimef(:,:)    &
                           ,Statein%tgrs(:,:)                           &
                           ,Stateout%gq0(:,:,Model%ntcw)                &
                           ,Stateout%gq0(:,:,Model%ntrw)                &
                           ,Stateout%gq0(:,:,Model%ntsw)                &
                           ,Stateout%gq0(:,:,Model%ntiw)                &
                           ,Stateout%gq0(:,:,Model%ntgl)                &
                           ,Model%MICROPHYSICS,Model%SPEC_ADV           &
                           ,NTIMESTEP,LM,ims,ime)
!
        ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Radiation
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
        IF(CALL_SHORTWAVE.OR.CALL_LONGWAVE)THEN
          DO I=IMS,IME
            CALL RADIATION(NTIMESTEP_RAD                                &
                          ,Model%dtp,JULDAY,JULYR,XTIME,JULIAN          &
                          ,Model%idat(5),Model%NPHS                     &
                          ,Grid%xlat(I:I),Grid%xlon(I:I)                &
                          ,Model%NRADS,Model%NRADL                      &
                          ,Statein%prsi(I:I,:),Statein%prsl(I:I,:)      &
                          ,Statein%tgrs(I:I,:),Stateout%gq0(I:I,:,1)    &
                          ,Sfcprop%THS(I:I),Radtend%sfalb(I:I)          &
                          ,Stateout%gq0(i:i,:,Model%ntcw)               &
                          ,Stateout%gq0(i:i,:,Model%ntrw)               &
                          ,Stateout%gq0(i:i,:,Model%ntiw)               &
                          ,Stateout%gq0(i:i,:,Model%ntsw)               &
                          ,Stateout%gq0(i:i,:,Model%ntgl)               &
                          ,Stateout%gq0(i:i,:,Model%ntinc)              &
                          ,Model%F_QC,Model%F_QR,Model%F_QI             &
                          ,Model%F_QS,Model%F_QG,Model%F_NI             &
                          ,Sfcprop%sm(I:I),Diag%cldfra(I:I,:)           &
                          ,Radtend%lwhc(I:I,:),Radtend%swhc(I:I,:)      &
                          ,Diag%RLWIN(I:I),Diag%DSWSFCI(I:I)            &
                          ,Diag%fluxr(i:i,32),Diag%fluxr(i:i,30)        &
                          ,Diag%fluxr(i:i,31),Diag%fluxr(i:i,33)        &
                          ,Diag%fluxr(i:i,28),Diag%fluxr(i:i,29)        &
                          ,Diag%USWSFCI(I:I)                            &
                          ,Diag%fluxr(i:i,1),Diag%fluxr(i:i,2)          &
                          ,Radtend%coszen(I:I),Diag%DLWSFCI(I:I)        &
                          ,Diag%CFRACL(I:I),Diag%CFRACM(I:I)            &
                          ,Diag%CFRACH(I:I)                             &
                          ,Diag%ACFRST(I:I),Diag%ACFRCV(I:I)            &
                          ,Diag%cuppt(I:I),Sfcprop%weasd(I:I)           &
                          ,Diag%htop(I:I),Diag%hbot(I:I)                &
                          ,Model%SHORTWAVE,Model%LONGWAVE               &
                          ,Model%CLDFRACTION,Grid%dx(i:i)               &
!---- RRTM part ---------------------------------------------------------
                          ,Model%jdat                                   &
                          ,Stateout%gq0(I:I,:,2),Stateout%gq0(I:I,:,3)  &
                          ,Statein%f_ice(I:I,:),Statein%f_rain(I:I,:)   &
                          ,Statein%f_rimef(I:I,:)                       &
                          ,Sfcprop%snowd(I:I),Sfcprop%tsfc(I:I)         &
                          ,Sfcprop%zorl(I:I),Sfcprop%sice(I:I)          &
                          ,Sfcprop%snoalb(I:I)                          &
                          ,Sfcprop%hprim(I:I),Statein%vvl(I:I,:)        &
                          ,Sfcprop%alvsf(I:I),Sfcprop%alnsf(I:I)        &  ! vis+uv & near IR beam albedos
                          ,Sfcprop%alvwf(I:I),Sfcprop%alnwf(I:I)        &  ! vis+uv & near IR diffuse albedos
                          ,Sfcprop%sncovr(I:I)                          &
!------------------------------------------------------------------------
                          ,LM,i,i,mype)
          ENDDO
        ENDIF
!
!-----------------------------------------------------------------------
!***  Empty the ACFRST and ACFRCV accumulation arrays if it is time
!***  to do so prior to their being updated by the radiation.
!-----------------------------------------------------------------------
!
        IF(MOD(NTIMESTEP,Model%NCLOD)==0)THEN
          Diag%ACFRST(:)=0.
          Diag%ACFRCV(:)=0.
        ENDIF
!
!-----------------------------------------------------------------------
!***  Update the temperature with the radiative tendency.
!-----------------------------------------------------------------------
!
        CALL RDTEMP(NTIMESTEP,Model%dtp,JULDAY,JULYR,Model%idat(5)      &
                   ,Grid%xlat(:),Grid%xlon(:)                           &
                   ,Radtend%CZEN(:),Radtend%coszen(:)                   &
                   ,Statein%tgrs(:,:)                                   &
                   ,Radtend%swhc(:,:),Radtend%lwhc(:,:)                 &
                   ,ims,ime,LM)
!
!-----------------------------------------------------------------------
!*** Initialize next step's temperature (new temperature after radiation)
!-----------------------------------------------------------------------
!
        Stateout%gt0 = Statein%tgrs
!
!-----------------------------------------------------------------------
!***  Empty the accumulators of sfc energy flux and sfc hydrology if
!***  it is time to do so prior to their being updated by turbulence.
!-----------------------------------------------------------------------
!
        IF(MOD(NTIMESTEP,Model%NRDLW)==0)THEN
          Diag%DLWSFC(:) =0.
          Diag%ULWSFC(:)=0.
          Diag%ARDLW =0.
        ENDIF
!
        IF(MOD(NTIMESTEP,Model%NRDSW)==0)THEN
          Diag%fluxr(:,4)=0.
          Diag%fluxr(:,3)=0.
          Diag%fluxr(:,23)=0.
          Diag%ARDSW =0.
        ENDIF
!
        IF(MOD(NTIMESTEP,Model%NSRFC)==0)THEN
          Diag%dtsfc(:)=0.
          Diag%dqsfc(:)=0.
          Diag%gflux(:)=0.
          Sfcprop%snopcx(:)=0.
          Diag%POTFLX(:)=0.
          Diag%ASRFC =0.
        ENDIF
!
        IF(MOD(NTIMESTEP,Model%NPREC)==0)THEN
          Sfcprop%acsnow(:)=0.
          Sfcprop%acsnom(:)=0.
          Diag%srunoff(:)=0.
          Diag%runoff(:)=0.
          Diag%POTEVP(:)=0.
        ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Turbulence, Sfc Layer, and Land Surface
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
        IF(CALL_TURBULENCE)THEN
!
          CALL TURBL(NTIMESTEP,Model%dtp,Model%NPHS                     &
                    ,Statein%ugrs(:,:),Statein%vgrs(:,:)                &
                    ,Statein%prsi(:,:),Statein%prsl(:,:)                &
                    ,Sfcprop%sm(:),Radtend%CZEN(:),Radtend%coszen(:)    &
                    ,Diag%DLWSFCI(:),Diag%RLWIN(:),Diag%DSWSFCI(:)      &
                    ,Diag%ULWSFCI(:) ,Radtend%lwhc(:,:)                 &
                    ,Radtend%swhc(:,:),Statein%tgrs(:,:)                &
                    ,Stateout%gq0(:,:,1),Stateout%gq0(:,:,2)            &
                    ,Statein%f_ice(:,:)                                 &
                    ,Statein%f_rain(:,:),Statein%f_rimef(:,:)           &
                    ,Diag%sr(:),Diag%Q2(:,:),DTDT(:,:),DUDT(:,:)        &
                    ,DVDT(:,:),Sfcprop%THS(:),Sfcprop%tsfc(:)           &
                    ,Sfcprop%SST(:),Sfcprop%tprcp(:),Sfcprop%weasd(:)   &
                    ,Sfcprop%sncovr(:)                                  &
                    ,Stateout%gq0(:,:,Model%ntcw)                       &
                    ,Stateout%gq0(:,:,Model%ntrw)                       &
                    ,Stateout%gq0(:,:,Model%ntiw)                       &
                    ,Stateout%gq0(:,:,Model%ntsw)                       &
                    ,Stateout%gq0(:,:,Model%ntgl)                       &
                    ,Model%F_QC,Model%F_QR,Model%F_QI                   &
                    ,Model%F_QS,Model%F_QG                              &
                    ,Sfcprop%oro(:),Sfcprop%zorl(:),Sfcprop%Z0BASE(:)   &
                    ,Sfcprop%uustar(:),Diag%hpbl(:),Diag%XLEN_MIX(:,:)  &
                    ,Diag%RMOL(:),Diag%chh(:),Diag%cmm(:)               &
                    ,Diag%AKHS_OUT(:),Diag%AKMS_OUT(:)                  &
                    ,Sfcprop%THZ0(:),Sfcprop%QZ0(:),Sfcprop%UZ0(:)      &
                    ,Sfcprop%VZ0(:),Sfcprop%QSH(:),Sfcprop%stc(:,:)     &
                    ,Sfcprop%smc(:,:),Sfcprop%canopy(:),Diag%SMSTAV(:)  &
                    ,Diag%soilm(:),Diag%srunoff(:),Diag%runoff(:)       &
                    ,IVGTYP(:),ISLTYP(:),Sfcprop%vfrac(:),Diag%gfluxi(:)&
                    ,Diag%SFCEXC(:),Sfcprop%acsnow(:)                   &
                    ,Sfcprop%acsnom(:),Sfcprop%snopcx(:)                &
                    ,Sfcprop%sice(:),Sfcprop%tg3(:),Diag%SOILTB(:)      &
                    ,Sfcprop%ALBASE(:),Sfcprop%snoalb(:)                &
                    ,Radtend%sfalb(:),Sfcprop%slc(:,:),Sfcprop%snowd(:) &
                    ,Diag%EPSR(:),Diag%u10m(:),Diag%v10m(:)             &
                    ,Diag%TH10(:),Diag%Q10(:),Diag%tshltr(:)            &
                    ,Diag%qshltr(:),Diag%pshltr(:),Diag%psurf(:)        &
                    ,Sfcprop%t2m(:),Diag%dtsfci(:),Diag%dqsfci(:)       &
                    ,Diag%dtsfc(:),Diag%dqsfc(:),Diag%EP(:),Diag%EPI(:) &
                    ,Diag%POTEVP(:),Diag%POTFLX(:),Diag%gflux(:)        &
                    ,Diag%evbsa(:),Diag%evcwa(:),Diag%APHTIM            &
                    ,Diag%ARDSW,Diag%ARDLW,Diag%ASRFC,Sfcprop%CROT(:)   &
                    ,Sfcprop%SROT(:),Diag%MIXHT(:)                      &
                    ,Sfcprop%hprime(:, 1),Sfcprop%hprime(:, 2)          &
                    ,Sfcprop%hprime(:, 3),Sfcprop%hprime(:, 4)          &
                    ,Sfcprop%hprime(:, 5),Sfcprop%hprime(:, 6)          &
                    ,Sfcprop%hprime(:, 7),Sfcprop%hprime(:, 8)          &
                    ,Sfcprop%hprime(:, 9),Sfcprop%hprime(:,10)          &
                    ,Sfcprop%hprime(:,11),Sfcprop%hprime(:,12)          &
                    ,Sfcprop%hprime(:,13),Sfcprop%hprime(:,14)          &
                    ,Model%CDMB,Model%CLEFF,Model%SIGFAC                &
                    ,Model%FACTOP,Model%RLOLEV,Model%DPMIN              &
                    ,Diag%USWSFCI(:),Diag%fluxr(:,2),Diag%fluxr(:,1)    &
                    ,Diag%fluxr(:,4),Diag%fluxr(:,3),Diag%fluxr(:,23)   &
                    ,Diag%DLWSFC(:),Diag%ULWSFC(:)                      &
                    ,Model%GWDFLG,.false.,Diag%DDATA(:),Model%UCMCALL   &
                    ,Model%IVEGSRC,Model%TURBULENCE,Model%SFC_LAYER     &
                    ,Model%LAND_SURFACE,Model%MICROPHYSICS,LISS_RESTART &
                    ,Model%VAR_RIC,Model%COEF_RIC_L,Model%COEF_RIC_S    &
                    ,Model%DISHEAT,Model%ALPHA,Model%SFENTH             &
                    ,ims,ime,LM)
!
!-----------------------------------------------------------------------
!***  Update temperature & wind with the turbulence tendencies
!-----------------------------------------------------------------------
          Stateout%gt0(:,:) = Stateout%gt0(:,:) + DTDT(:,:) * Model%NPHS * Model%dtp
          Stateout%gu0(:,:) = Stateout%gu0(:,:) + DUDT(:,:) * Model%NPHS * Model%dtp
          Stateout%gv0(:,:) = Stateout%gv0(:,:) + DVDT(:,:) * Model%NPHS * Model%dtp
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Empty the accumulators of precipitation and latent heating if is
!***  is time prior to their being updated by convection/microphysics.
!-----------------------------------------------------------------------
!
        IF(MOD(NTIMESTEP,Model%NPREC)==0)THEN
          Diag%totprcp(:)=0.
          Diag%cnvprcp(:)=0.
        ENDIF
!
        IF(MOD(NTIMESTEP,Model%NHEAT)==0)THEN
          Diag%AVCNVC=0.
          Diag%AVRAIN=0.
!
          DO L=1,LM
            Diag%TRAIN(:,L)=0.
            Diag%TCUCN(:,L)=0.
          ENDDO
        ENDIF
!
!-----------------------------------------------------------------------
!***  1 of 3 calls to CLTEND, save Told array before convection & microphysics
!-----------------------------------------------------------------------
!
        cld_tend1: IF(CALL_PRECIP .AND. Model%NPRECIP>1) THEN
          DO K=1,LM
            Statein%Told(:,K)=Stateout%gt0(:,K)
          ENDDO
        ENDIF cld_tend1
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Convection
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
        IF(CALL_PRECIP.AND.Model%CONVECTION/='none')THEN
!
          CALL CUCNVC(NTIMESTEP,Model%dtp,Model%NPRECIP,Model%NRADS     &
                     ,Model%NRADL,Model%MINUTES_HISTORY,Model%ENTRAIN   &
                     ,Model%NEWALL,Model%NEWSWAP,Model%NEWUPUP          &
                     ,Model%NODEEP,Model%FRES,Model%FR,Model%FSL        &
                     ,Model%FSS,Diag%CLDEFI(:)                          &
                     ,Stateout%gu0(:,:),Stateout%gv0(:,:)               &
                     ,Statein%f_ice(:,:),Statein%f_rain(:,:)            &
                     ,Stateout%gq0(:,:,Model%ntcw)                      &
                     ,Stateout%gq0(:,:,Model%ntrw)                      &
                     ,Stateout%gq0(:,:,Model%ntiw)                      &
                     ,Stateout%gq0(:,:,Model%ntsw)                      &
                     ,Stateout%gq0(:,:,Model%ntgl)                      &
                     ,Model%F_QC,Model%F_QR,Model%F_QI,Model%F_QS       &
                     ,Model%F_QG,Statein%prsi(:,:),Statein%prsl(:,:)    &
                     ,Grid%area(:),Stateout%gt0(:,:)                    &
                     ,Stateout%gq0(:,:,1)                               &
                     ,Stateout%gq0(:,:,2),Diag%TCUCN(:,:)               &
                     ,Statein%vvl(:,:),Sfcprop%oro(:)                   &
                     ,Sfcprop%tprcp(:),Diag%totprcp(:),Diag%cnvprcp(:)  &
                     ,Diag%cuppt(:),Diag%rainc(:),Diag%CNVBOT(:)        &
                     ,Diag%CNVTOP(:),Sfcprop%sm(:)                      &
                     ,Diag%htop(:),Diag%htopd(:),Diag%htops(:)          &
                     ,Diag%hbot(:),Diag%hbotd(:),Diag%hbots(:)          &
                     ,Diag%AVCNVC,Diag%ACUTIM                           &
                     ,Diag%DSWSFCI(:),Diag%USWSFCI(:)                   &
                     ,Model%CONVECTION,Model%MICROPHYSICS               &
                     ,Sfcprop%sice(:),Diag%dqsfci(:),Diag%dtsfci(:)     &
                     ,Diag%hpbl(:),DTDT(:,:),DUDT(:,:),DVDT(:,:)        &
                     ,Model%SAS_MOMMIX,Model%SAS_PGCON                  &
                     ,Model%SAS_MASS_FLUX                               &
                     ,Model%SAS_SHALCONV,Model%SAS_SHAL_PGCON           &
                     ,ims,ime,LM)
!
!-----------------------------------------------------------------------
!***  Update temperature & wind with the convection tendencies
!-----------------------------------------------------------------------
          Stateout%gt0(:,:) = Stateout%gt0(:,:) + DTDT(:,:) * Model%NPRECIP * Model%dtp
          Stateout%gu0(:,:) = Stateout%gu0(:,:) + DUDT(:,:) * Model%NPRECIP * Model%dtp
          Stateout%gv0(:,:) = Stateout%gv0(:,:) + DVDT(:,:) * Model%NPRECIP * Model%dtp
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Microphysics
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
        IF(CALL_PRECIP)THEN
!
          CALL GSMDRIVE(NTIMESTEP,Model%dtp,Model%NPRECIP               &
                       ,Sfcprop%sm(:),Sfcprop%oro(:)                    &
                       ,Statein%prsi(:,:),Statein%prsl(:,:)             &
                       ,Stateout%gt0(:,:)                               &
                       ,Stateout%gq0(:,:,1),Stateout%gq0(:,:,2)         & !rv CW vs. QC
                       ,Diag%TRAIN(:,:),Diag%sr(:)                      &
                       ,Statein%f_ice(:,:),Statein%f_rain(:,:)          &
                       ,Statein%f_rimef(:,:)                            &
                       ,Stateout%gq0(:,:,Model%ntcw)                    &
                       ,Stateout%gq0(:,:,Model%ntrw)                    &
                       ,Stateout%gq0(:,:,Model%ntiw)                    &
                       ,Stateout%gq0(:,:,Model%ntsw)                    &
                       ,Stateout%gq0(:,:,Model%ntgl)                    &
                       ,Stateout%gq0(:,:,Model%ntinc)                   &
                       ,Stateout%gq0(:,:,Model%ntrnc)                   &
                       ,Model%F_QC,Model%F_QR,Model%F_QI,Model%F_QS     &
                       ,Model%F_QG,Model%F_NI,Model%F_NR                &
                       ,Model%has_reqc,Model%has_reqi,Model%has_reqs    &
                       ,Sfcprop%tprcp(:),Diag%totprcp(:)                &
                       ,Diag%AVRAIN,Statein%refl_10cm(:,:)              &
                       ,Statein%re_cloud(:,:)                           &
                       ,Statein%re_ice(:,:),Statein%re_snow(:,:)        &
                       ,Model%MICROPHYSICS,Model%RHGRD,Diag%TP1(:,:)    &
                       ,Diag%QP1(:,:),Diag%PSP1(:)                      &
                       ,ims,ime,LM)
!
!-----------------------------------------------------------------------
!***  2 of 3 calls to CLTEND, calculate Tadj and replace T with Told
!-----------------------------------------------------------------------
!
        cld_tend2: IF(Model%NPRECIP>1) THEN
          DO K=1,LM
            Statein%Tadj(:,K)=(Stateout%gt0(:,K)-Statein%Told(:,K))/REAL(Model%NPRECIP)
            Stateout%gt0(:,K)=Statein%Told(:,K)
          ENDDO
        ENDIF  cld_tend2
!
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  3 of 3 calls to CLTEND, incremental updates of T using Told & Tadj
!-----------------------------------------------------------------------
!
        cld_tend3: IF(Model%NPRECIP>1) THEN
          DO K=1,LM
            Stateout%gt0(:,K)=Stateout%gt0(:,K)+Statein%Tadj(:,K)
          ENDDO
        ENDIF  cld_tend3
!
!-----------------------------------------------------------------------
!***  Prevent supersaturation w/r/t water and smooth temperature profiles
!     if lapse rates are steeper than dry adiabatic above lowest levels.
!-----------------------------------------------------------------------
!
        CALL TQADJUST(Stateout%gt0(:,:),Stateout%gq0(:,:,1)             &
                     ,Stateout%gq0(:,:,Model%ntcw)                      &
                     ,Stateout%gq0(:,:,2),Statein%f_ice(:,:)            &
                     ,Statein%f_rain(:,:)                               &
                     ,Statein%prsi(:,:),Statein%prsl(:,:)               &
                     ,Model%SPEC_ADV,Model%RHGRD                        &
                     ,LM,ims,ime)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SOLVER_RUN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_SOLVER_GRID_COMP
!
!-----------------------------------------------------------------------
