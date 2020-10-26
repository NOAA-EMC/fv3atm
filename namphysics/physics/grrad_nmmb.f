!!!!!  ==========================================================  !!!!!
!!!!!             'module_radiation_driver' descriptions           !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!   this is the radiation driver module.  it prepares atmospheric      !
!   profiles and invokes main radiation calculations.                  !
!                                                                      !
!   in module 'module_radiation_driver' there are twe externally       !
!   callable subroutine:                                               !
!                                                                      !
!      'radinit'    -- initialization routine                          !
!         input:                                                       !
!           ( si, NLAY, iflip, idate, jdate, ICTM, ISOL, ICO2,         !
!             IAER, IALB, IEMS, ICWP, NP3D, isubcsw, isubclw,          !
!             iovrsw, iovrlw, me )                                     !
!         output:                                                      !
!           ( none )                                                   !
!                                                                      !
!      'grrad'      -- setup and invoke main radiation calls           !
!         input:                                                       !
!          ( prsi,prsl,prslk,tgrs,qgrs,oz,vvl,slmsk,                   !
!            xlon,xlat,tsfc,snowd,sncovr,snoalb,zorl,hprim,            !
!            alvsf,alnsf,alvwf,alnwf,facsf,facwf,fice,tisfc,           !
!            solcon,coszen,coszdg,k1oz,k2oz,facoz,                     !
!            cv,cvt,cvb,iovrsw,iovrlw,fcice,frain,rrime,flgmin,        !
!            icsdsw,icsdlw, np3d,ntoz, NTRAC,NFXR,                     !
!            dtlw,dtsw, lsswr,lslwr,lssav,sashal,norad_precip,         !
!            crick_proof, ccnorm,                                      !
!            IX, IM, LM, iflip, me, lprnt, ipt, kdt,                   !
!         output:                                                      !
!            htrsw,topfsw,sfcfsw,sfalb,                                !
!            htrlw,topflw,sfcflw,tsflw,semis,cldcov,cldsa              !
!         input/output:                                                !
!            fluxr                                                     !
!         optional output:                                             !
!            HTRSWB,HTRLWB)                                            !
!                                                                      !
!                                                                      !
!   external modules referenced:                                       !
!       'module machine'                    in 'machine.f'             !
!       'module funcphys'                   in 'funcphys.f'            !
!       'module physcons'                   in 'physcons.f             !
!                                                                      !
!       'module module_radiation_gases'     in 'radiation_gases.f'     !
!       'module module_radiation_aerosols'  in 'radiation_aerosols.f'  !
!       'module module_radiation_surface'   in 'radiation_surface.f'   !
!       'module module_radiation_clouds'    in 'radiation_clouds.f'    !
!                                                                      !
!       'module module_radsw_cntr_para'     in 'radsw_xxxx_param.f'    !
!       'module module_radsw_parameters'    in 'radsw_xxxx_param.f'    !
!       'module module_radsw_main'          in 'radsw_xxxx_main.f'     !
!                                                                      !
!       'module module_radlw_cntr_para'     in 'radlw_xxxx_param.f'    !
!       'module module_radlw_parameters'    in 'radlw_xxxx_param.f'    !
!       'module module_radlw_main'          in 'radlw_xxxx_main.f'     !
!                                                                      !
!    where xxxx may vary according to different scheme selection       !
!                                                                      !
!                                                                      !
!   program history log:                                               !
!     mm-dd-yy    ncep         - created program grrad                 !
!     08-12-03    yu-tai hou   - re-written for modulized radiations   !
!     11-06-03    yu-tai hou   - modified                              !
!     01-18-05    s. moorthi   - NOAH/ICE model changes added          !
!     05-10-05    yu-tai hou   - modified module structure             !
!     12-xx-05    s. moorthi   - sfc lw flux adj by mean temperature   !
!     02-20-06    yu-tai hou   - add time variation for co2 data, and  !
!                                solar const. add sfc emiss change     !
!     03-21-06    s. Moorthi   - added surface temp over ice           !
!     07-28-06    yu-tai hou   - add stratospheric vocanic aerosols    !
!     03-14-07    yu-tai hou   - add generalized spectral band interp  !
!                                for aerosol optical prop. (sw and lw) !
!     04-10-07    yu-tai hou   - spectral band sw/lw heating rates     !
!     05-04-07    yu-tai hou   - make options for clim based and modis !
!                                based (h. wei and c. marshall) albedo !
!     09-05-08    yu-tai hou   - add the initial date and time 'idate' !
!                    and control param 'ICTM' to the passing param list!
!                    to handel different time/date requirements for    !
!                    external data (co2, aeros, solcon, ...)           !
!     10-10-08    yu-tai hou   - add the ICTM=-2 option for combining  !
!                    initial condition data with seasonal cycle from   !
!                    climatology.                                      !
!     03-12-09    yu-tai hou   - use two time stamps to keep tracking  !
!                    dates for init cond and fcst time. remove volcanic!
!                    aerosols data in climate hindcast (ICTM=-2).      !
!     03-16-09    yu-tai hou   - included sub-column clouds approx.    !
!                    control flags isubcsw/isubclw in initialization   !
!                    subroutine. passed auxiliary cloud control arrays !
!                    icsdsw/icsdlw (if isubcsw/isubclw =2, it will be  !
!                    the user provided permutation seeds) to the sw/lw !
!                    radiation calculation programs. also moved cloud  !
!                    overlapping control flags iovrsw/iovrlw from main !
!                    radiation routines to the initialization routines.!
!     04-02-09    yu-tai hou   - modified surface control flag iems to !
!                    have additional function of if the surface-air    !
!                    interface have the same or different temperature  !
!                    for radiation calculations.                       !
!     04-03-09    yu-tai hou   - modified to add lw surface emissivity !
!                    as output variable. changed the sign of sfcnsw to !
!                    be positive value denote solar flux goes into the !
!                    ground (this is needed to reduce sign confusion   !
!                    in other part of model)                           !
!     04-20-09    carlos perez - prepare driver for nmmb.  added option!
!                    of run the gfs's radiation on nmmb                !
!     09-09-09    fanglin yang (thru s.moorthi) added QME5 QME6 to E-20!
!     01-09-10    sarah lu     - added gocart option, revised grrad for!
!                    gocart coupling. calling argument modifed: ldiag3 !
!                    removed; cldcov/fluxr sequence changed; cldcov is !
!                    changed from accumulative to instant field and    !
!                    from input/output to output field                 !
!     01-24-10    sarah lu     - added aod to fluxr, added prslk and   !
!                    oz to setaer input argument (for gocart coupling),!
!                    added tau_gocart to setaer output argument (for,  !
!                    aerosol diag)                                     !
!     07-08-10    s.moorthi - updated the NEMS version for new physics !
!     07-28-10    yu-tai hou   - changed grrad interface to allow all  !
!                    components of sw/lw toa/sfc instantaneous values  !
!                    being passed to the calling program. moved the    !
!                    computaion of sfc net sw flux (sfcnsw) to the     !
!                    calling program. merged carlos' nmmb modification.!
!     07-30-10    s. moorthi - corrected some errors associated with   !
!                    unit changes                                      !
!     12-02-10    s. moorthi/y. hou - removed the use of aerosol flags !
!                    'iaersw' 'iaerlw' from radiations and replaced    !
!                    them by using the runtime variable iaerflg and    !
!                    laswflg defined in module radiation_aerosols.     !
!                    also replaced param nspc in grrad with the use of !
!                    max_num_gridcomp in module radiation_aerosols.    !
!     01-03-11    y. hou     - added sea/land madk 'slmsk' to the      !
!                    argument list of subrotine setaer call for the    !
!                    newly modified horizontal bi-linear interpolation !
!                    in climatological aerosols schem.                 !
!     09-30-11    H.M. Lin   - rename "grrad" to "grrad_nmmb" for using!
!                    in regional to avoid the confliction with global  !
!                    when compile the nems. Also change names of the   !
!                    related subroutines to "astronomy_nmmb" &         !
!                    "solinit_nmmb"                                    !
!     03-15-12    H.M. Lin   - changed grrad interface to allow 'cldsa'!
!                    (fraction of clouds for low, mid, hi, tot, bl)    !
!                     passed to the calling program.                   !
!     09-25-12    y. hou     - added optional extra top layer 'LTP' on !
!                    top of low ceiling models                         !
!                    (** to cover the gases atop, expecially O3)       !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!



!=============================================!
      module module_radiation_driver_nmmb     !
!.............................................!
!
      use physparam
      use physcons,                 only : con_eps, con_epsm1,          &
     &                                     con_fvirt, PI=>con_pi
      use funcphys,                 only : fpkapx

      use MODULE_CONSTANTS,         only : EPSQ

      use module_radiation_astronomy_nmmb,only : sol_init_nmmb,         &
     &                                     sol_update_nmmb, coszmn_nmmb
      use module_radiation_gases_nmmb,only : NF_VGAS, getgases, getozn, &
     &                                     gas_init, gas_update
      use module_radiation_aerosols_nmmb,only : NF_AESW, NF_AELW,setaer,&
     &                                     aer_init, aer_update
!    &,                                    NSPC1                        ! optn for aod output
      use module_radiation_surface_nmmb, only : NF_ALBD,sfc_init,setalb,&
     &                                     setemis
      use module_radiation_clouds_nmmb,  only : NF_CLDS, cld_init

      use module_radsw_parameters,  only : topfsw_type, sfcfsw_type,    &
     &                                     profsw_type,cmpfsw_type,NBDSW
      use module_radsw_main_nmmb,   only : rswinit,  swrad

      use module_radlw_parameters,  only : topflw_type, sfcflw_type,    &
     &                                     proflw_type, NBDLW
      use module_radlw_main_nmmb,   only : rlwinit,  lwrad

!
      implicit   none
!
      private

!  ---  version tag and last revision date
      character(40), parameter ::                                       &
     &   VTAGRAD='NCEP-Radiation_driver    v5.2  Jan 2013 '
!    &   VTAGRAD='NCEP-Radiation_driver    v5.1  Nov 2012 '
!    &   VTAGRAD='NCEP-Radiation_driver    v5.0  Aug 2012 '

!  ---  constant values
      real (kind=kind_phys) :: QMIN, QME5, QME6
!     parameter (QMIN=1.0e-10, QME5=1.0e-5,  QME6=1.0e-6,  EPSQ=1.0e-12)
      parameter (QMIN=1.0e-10, QME5=1.0e-7,  QME6=1.0e-7)
!     parameter (QMIN=1.0e-10, QME5=1.0e-20, QME6=1.0e-20, EPSQ=1.0e-12)
      real, parameter :: prsmin = 1.0e-6 ! toa pressure minimum value in mb (hpa)

!  ---  control flags set in subr radinit:
      integer :: itsfc  =0             ! flag for lw sfc air/ground interface temp setting

!  ---  data input control variables set in subr radupdate:
      integer :: month0=0,   iyear0=0,   monthd=0
      logical :: loz1st =.true.       ! first-time clim ozone data read flag

!  ---  optional extra top layer on top of low ceiling models (Y. Hou, 2012-09-25)
!     integer, parameter :: LTP = 0   ! do no add an extra top layer
      integer, parameter :: LTP = 1   ! add an extra top layer
      logical :: lextop = (LTP > 0)
!
!  ---  publicly accessible module programs:

      public radinit_nmmb, radupdate_nmmb, grrad_nmmb


! =================
      contains
! =================


!-----------------------------------
      subroutine radinit_nmmb                                           &
!...................................

!  ---  inputs:
     &     ( si, NLAY, me )
!  ---  outputs:
!          ( none )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:   radinit     initialization of radiation calculations    !
!                                                                       !
!                                                                       !
! program history log:                                                  !
!   08-14-2003   yu-tai hou   created                                   !
!                                                                       !
! usage:        call radinit                                            !
!                                                                       !
! attributes:                                                           !
!   language:  fortran 90                                               !
!   machine:   ibm sp                                                   !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input parameters:                                                     !
!   si               : model vertical sigma interface                   !
!   NLAY             : number of model vertical layers                  !
!   iflip            : control flag for direction of vertical index     !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   idate(8)         : ncep absolute date and time of initial condition !
!                      (yr, mon, day, t-zone, hr, min, sec, mil-sec)    !
!   jdate(8)         : ncep absolute date and time at fcst time         !
!                      (yr, mon, day, t-zone, hr, min, sec, mil-sec)    !
!   ICTM             :=yyyy#, external data time/date control flag      !
!                     =   -2: same as 0, but superimpose seasonal cycle !
!                             from climatology data set.                !
!                     =   -1: use user provided external data for the   !
!                             forecast time, no extrapolation.          !
!                     =    0: use data at initial cond time, if not     !
!                             available, use latest, no extrapolation.  !
!                     =    1: use data at the forecast time, if not     !
!                             available, use latest and extrapolation.  !
!                     =yyyy0: use yyyy data for the forecast time,      !
!                             no further data extrapolation.            !
!                     =yyyy1: use yyyy data for the fcst. if needed, do !
!                             extrapolation to match the fcst time.     !
!   ISOL             :=0: use a fixed solar constant value              !
!                     =1: use 11-year cycle solar constant table        !
!   ICO2             :=0: use prescribed global mean co2 (old  oper)    !
!                     =1: use observed co2 annual mean value only       !
!                     =2: use obs co2 monthly data with 2-d variation   !
!   IAER             : 3-digit aerosol flag (for volc, lw, sw)          !
!                     =  0: turn all aeros effects off (sw,lw,volc)     !
!                     =  1: use clim tropspheric aerosol for sw only    !
!                     = 10: use clim tropspheric aerosol for lw only    !
!                     = 11: use clim tropspheric aerosol for both sw/lw !
!                     =100: volc aerosol only for both sw and lw        !
!                     =101: volc and clim trops aerosol for sw only     !
!                     =110: volc and clim trops aerosol for lw only     !
!                     =111: volc and clim trops aerosol for both sw/lw  !
!                     = 2: gocart prognostic, without volc forcing      !
!                     =12: gocart prognostic, with volcanic forcing     !
!   IALB             : control flag for surface albedo schemes          !
!                     =0: climatology, based on surface veg types       !
!                     =1: modis retrieval based surface albedo scheme   !
!                     =2: use externally provided albedoes directly.    !
!   IEMS             : ab 2-digit control flag                          !
!                      a =0 set sfc air/ground t same for lw radiation  !
!                        =1 set sfc air/ground t diff for lw radiation  !
!                      b =0 use fixed sfc emissivity=1.0 (black-body)   !
!                        =1 use varying climtology sfc emiss (veg based)!
!                        =2 future development (not yet)                !
!   ICWP             : control flag for cloud generation schemes        !
!                     =-1: use diagnostic cloud scheme (GFDL type       !
!                     =0 : use diagnostic cloud scheme                  !
!                     =1 : use prognostic cloud scheme (default)        !
!   NP3D             :=3: ferrier's microphysics cloud scheme           !
!                     =4: zhao/carr/sundqvist microphysics cloud        !
!                     =5: nmmb ferrier+bmj microphysics scheme          !
!                     =8: Thompson microphysics scheme                  !
!   isubcsw/isubclw  : sub-column cloud approx control flag (sw/lw rad) !
!                     =0: with out sub-column cloud approximation       !
!                     =1: mcica sub-col approx. prescribed random seed  !
!                     =2: mcica sub-col approx. provided random seed    !
!   iovrsw/iovrlw    : control flag for cloud overlap (sw/lw rad)       !
!                     =0: random overlapping clouds                     !
!                     =1: max/ran overlapping clouds                    !
!   me               : print control flag                               !
!                                                                       !
!  outputs: (none)                                                      !
!                                                                       !
!  usage:       call radinit                                            !
!                                                                       !
!  subroutines called:    cldinit, aerinit, rlwinit, rswinit, gasinit   !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NLAY, me

      real (kind=kind_phys), intent(in) :: si(:)

!  ---  outputs: (none, to module variables)

!  ---  locals:
      integer :: RICWP, icldflg_org

!
!===> ...  begin here
!
!  --- ...  set up control variables
      itsfc  = iemsflg / 10             ! sfc air/ground temp control
      loz1st = (ioznflg == 0)           ! first-time clim ozone data read flag

      if (me == 0) then
!       print *,' NEW RADIATION PROGRAM STRUCTURES -- SEP 01 2004'
        print *,' NEW RADIATION PROGRAM STRUCTURES BECAME OPER. ',      &
     &          '  May 01 2007'
        print *, VTAGRAD                !print out version tag
        print *,' - Selected Control Flag settings: ICTMflg=',ictmflg,  &
     &    ' ISOLar =',isolar, ' ICO2flg=',ico2flg,' IAERflg=',iaerflg,  &
     &    ' IALBflg=',ialbflg,' IEMSflg=',iemsflg,' ICLDflg=',icldflg,  &
     &    ' ICMPHYS=',icmphys,' IOZNflg=',ioznflg
        print *,' IVFLIP=',ivflip,' IOVRSW=',iovrsw,' IOVRLW=',iovrlw,  &
     &    ' ISUBCSW=',isubcsw,' ISUBCLW=',isubclw
        print *,' LSASHAL=',lsashal,' LCRICK=',lcrick,' LCNORM=',lcnorm,&
     &    ' LNOPREC=',lnoprec
        print *,' LTP =',LTP,', add extra top layer =',lextop

        if ( ictmflg==0 .or. ictmflg==-2 ) then
          print *,'   Data usage is limited by initial condition!'
          print *,'   No volcanic aerosols'
        endif

        if ( isubclw == 0 ) then
          print *,' - ISUBCLW=',isubclw,' No McICA, use grid ',         &
     &            'averaged cloud in LW radiation'
        elseif ( isubclw == 1 ) then
          print *,' - ISUBCLW=',isubclw,' Use McICA with fixed',        &
     &            'permutation seeds for LW random number generator'
        elseif ( isubclw == 2 ) then
          print *,' - ISUBCLW=',isubclw,' Use McICA with random ',      &
     &            'permutation seeds for LW random number generator'
        else
          print *,' - ERROR!!! ISUBCLW=',isubclw,' is not a ',          &
     &            'valid option '
          stop
        endif

        if ( isubcsw == 0 ) then
          print *,' - ISUBCSW=',isubcsw,' No McICA, use grid ',         &
     &            'averaged cloud in SW radiation'
        elseif ( isubcsw == 1 ) then
          print *,' - ISUBCSW=',isubcsw,' Use McICA with fixed',        &
     &            'permutation seeds for SW random number generator'
        elseif ( isubcsw == 2 ) then
          print *,' - ISUBCSW=',isubcsw,' Use McICA with random ',      &
     &            'permutation seeds for SW random number generator'
        else
          print *,' - ERROR!!! ISUBCSW=',isubcsw,' is not a ',          &
     &            'valid option '
          stop
        endif

        if ( isubcsw /= isubclw ) then
          print *,' - *** Notice *** ISUBCSW /= ISUBCLW !!!',           &
     &            isubcsw, isubclw
        endif
      endif

!  --- ...  call astronomy initialization routine

      call sol_init_nmmb (  me )

!  --- ...  call aerosols initialization routine

      call aer_init ( NLAY, me )

!  --- ...  call co2 and other gases initialization routine

      call gas_init ( me )

!  --- ...  call surface initialization routine

      call sfc_init ( me )

!  --- ...  call cloud initialization routine

      call cld_init ( si, NLAY, me)


!==========================================================
! The following use for NP3D=5 GFDL cloud in radiation
!==========================================================

      if ( icmphys==5 ) then
         ricwp = 0
         icldflg_org = icldflg
         icldflg = ricwp
      endif

!  --- ...  call lw radiation initialization routine

      call rlwinit ( me )

!  --- ...  call sw radiation initialization routine

      call rswinit ( me )
!

!==========================================================
! The following use for NP3D=5 GFDL cloud in radiation
! return icldflg to initial setting
!==========================================================

      if ( icmphys==5 ) then
         icldflg = icldflg_org
      endif

      
      return
!...................................
      end subroutine radinit_nmmb
!-----------------------------------


!-----------------------------------
      subroutine radupdate_nmmb                                              &
!...................................

!  ---  inputs:
     &     ( idate, jdate, deltsw, deltim, lsswr, me,                   &
!  ---  outputs:
     &       slag,sdec,cdec,solcon                                      &
     &     )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:   radupdate   calls many update subroutines to check and  !
!   update radiation required but time varying data sets and module     !
!   variables.                                                          !
!                                                                       !
! usage:        call radupdate                                          !
!                                                                       !
! attributes:                                                           !
!   language:  fortran 90                                               !
!   machine:   ibm sp                                                   !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input parameters:                                                     !
!   idate(8)       : ncep absolute date and time of initial condition   !
!                    (yr, mon, day, t-zone, hr, min, sec, mil-sec)      !
!   jdate(8)       : ncep absolute date and time at fcst time           !
!                    (yr, mon, day, t-zone, hr, min, sec, mil-sec)      !
!   deltsw         : sw radiation calling frequency in seconds          !
!   deltim         : model timestep in seconds                          !
!   lsswr          : logical flags for sw radiation calculations        !
!   me             : print control flag                                 !
!                                                                       !
!  outputs:                                                             !
!   slag           : equation of time in radians                        !
!   sdec, cdec     : sin and cos of the solar declination angle         !
!   solcon         : sun-earth distance adjusted solar constant (w/m2)  !
!                                                                       !
!  module variables:                                                    !
!   isolar   : solar constant cntrl  (see ISOL description)             !
!   ictmfl   : flag for initial condition gh-gas data source            !
!   loz1st   : first-time clim ozone data read flag                     !
!                                                                       !
!  subroutines called: sol_update, aer_update, gas_update               !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: idate(:), jdate(:), me
      logical, intent(in) :: lsswr

      real (kind=kind_phys), intent(in) :: deltsw, deltim

!  ---  outputs:
      real (kind=kind_phys), intent(out) :: slag, sdec, cdec, solcon

!  ---  locals:
      integer :: iyear, imon, iday, ihour
      integer :: kyear, kmon, kday, khour

      logical :: lmon_chg       ! month change flag
      logical :: lco2_chg       ! cntrl flag for updating co2 data
      logical :: lsol_chg       ! cntrl flag for updating solar constant
!
!===> ...  begin here
!

!  --- ...  time stamp at fcst time

      iyear = jdate(1)
      imon  = jdate(2)
      iday  = jdate(3)
      ihour = jdate(5)

!  --- ...  set up time stamp used for green house gases (** currently co2 only)

      if ( ictmflg==0 .or. ictmflg==-2 ) then  ! get external data at initial condition time
        kyear = idate(1)
        kmon  = idate(2)
        kday  = idate(3)
        khour = idate(5)
      else                           ! get external data at fcst or specified time
        kyear = iyear
        kmon  = imon
        kday  = iday
        khour = ihour
      endif   ! end if_ictmflg_block

      if ( month0 /= imon ) then
        lmon_chg = .true.
        month0 = imon
      else
        lmon_chg = .false.
      endif

!  --- ...  call astronomy update routine, yearly update, no time interpolation

      if (lsswr) then

        if  ( isolar == 0 .or. isolar == 10 ) then
          lsol_chg = .false.
        elseif ( iyear0/=iyear ) then
          lsol_chg = .true.
        else
          lsol_chg = ( isolar==4 .and. lmon_chg )
        endif
        iyear0 = iyear

        call sol_update_nmmb                                            &
!  ---  input:
     &     ( jdate,kyear,deltsw,deltim,lsol_chg, me,                    &
!  ---  outputs:
     &       slag,sdec,cdec,solcon                                      &
     &     )

      endif  ! end_if_lsswr_block

!  --- ...  call aerosols update routine, monthly update, no time interpolation

      if ( lmon_chg ) then
        call aer_update ( iyear, imon, me )
      endif

!  --- ...  call co2 and other gases update routine

      if ( monthd /= kmon ) then
        monthd = kmon
        lco2_chg = .true.
      else
        lco2_chg = .false.
      endif

      call gas_update ( kyear,kmon,kday,khour,loz1st,lco2_chg, me )

      if ( loz1st ) loz1st = .false.

!  --- ...  call surface update routine (currently not needed)
!     call sfc_update ( iyear, imon, me )

!  --- ...  call clouds update routine (currently not needed)
!     call cld_update ( iyear, imon, me )
!
      return
!...................................
      end subroutine radupdate_nmmb
!-----------------------------------



!-----------------------------------
      subroutine grrad_nmmb                                             &
!...................................

!  ---  inputs:
     &     ( prsi,plyr1,plvl1,prslk1,tgrs,qgrs,tracer1,slmsk,rhly1,     &
     &       qlyr1,tvly1,                                               &
     &       xlon,xlat,tsfc,snowd,sncovr,snoalb,zorl,hprim,             &
     &       alvsf,alnsf,alvwf,alnwf,facsf,facwf,                       &  ! processed inside
     &       SALBEDO,SM,                                                &  ! input for albedo cal.
     &       fice,tisfc,                                                &
     &       sinlat,coslat,solhr,jdate,solcon,                          &
     &       dtswav,nrads,                                              &  ! extra input
     &       icsdsw,icsdlw, ntoz, NTRAC,NFXR,                           &
     &       CPATHFAC4LW,                                               &  ! enhance factor of cloud depth for LW
     &       dtlw,dtsw, lsswr,lslwr,lssav,                              &
     &       IX, IM, LM, me, lprnt, ipt, kdt,                           &
     &       clouds1,cldsa,mtopa,mbota,                                 &  !! clouds properties for radiation
!  ---  outputs:
     &       htrsw,topfsw,sfcfsw,sfalb,coszen,coszdg,                   &
     &       htrlw,topflw,sfcflw,tsflw,semis,                           &
!  ---  input/output:
     &       fluxr                                                      &
!! ---  optional outputs:
     &,      HTRSWB,HTRLWB                                              &
     &     )

! =================   subprogram documentation block   ================ !
!                                                                       !
!    this program is the driver of radiation calculation subroutines. * !
!    It sets up profile variables for radiation input, including      * !
!    clouds, surface albedos, atmospheric aerosols, ozone, etc.       * !
!                                                                     * !
!                                                                     * !
!    usage:        call grrad                                         * !
!                                                                     * !
!    subprograms called:                                              * !
!                  setalb, setemis, setaer, getozn, getgases,         * !
!                  swrad, lwrad, fpvs                                 * !
!                                                                     * !
!    attributes:                                                      * !
!      language:   fortran 90                                         * !
!      machine:    ibm-sp, sgi                                        * !
!                                                                     * !
!                                                                     * !
!  ====================  defination of variables  ====================  !
!                                                                       !
!    input variables:                                                   !
!      prsi  (IX,LM+1) : model level pressure in Pa                     !
!      prsl  (IX,LM)   : model layer mean pressure in Pa                !
!      prslk (IX,LM)   : Exner function                                 !
!      tgrs  (IX,LM)   : model layer mean temperature in k              !
!      qgrs  (IX,LM)   : layer specific humidity in gm/gm               !
!      oz  (IX,LM,NTRAC):layer ozone mass mixing ratio                  !
!      vvl   (IX,LM)   : layer mean vertical velocity in Pa/sec         !
!      slmsk (IM)      : sea/land mask array (sea:0,land:1,sea-ice:2)   !
!      xlon,xlat (IM)  : grid longitude/latitude in radians             !
!      tsfc  (IM)      : surface temperature in k                       !
!      snowd (IM)      : snow depth water equivalent in mm              !
!      sncovr(IM)      : snow cover in fraction                         !
!      snoalb(IM)      : maximum snow albedo in fraction                !
!      zorl  (IM)      : surface roughness in cm                        !
!      hprim (IM)      : topographic standard deviation in m            !
!      alvsf (IM)      : mean vis albedo with strong cosz dependency    !
!      alnsf (IM)      : mean nir albedo with strong cosz dependency    !
!      alvwf (IM)      : mean vis albedo with weak cosz dependency      !
!      alnwf (IM)      : mean nir albedo with weak cosz dependency      !
!      facsf (IM)      : fractional coverage with strong cosz dependen  !
!      facwf (IM)      : fractional coverage with weak cosz dependency  !
!      fice  (IM)      : ice fraction over open water grid              !
!      tisfc (IM)      : surface temperature over ice fraction          !
!      solcon          : solar constant (sun-earth distant adjusted)    !
!      coszen(IM)      : mean cos of zenith angle over rad call period  !
!      coszdg(IM)      : daytime mean cosz over rad call period         !
!      k1oz,k2oz,facoz : parameters for climatological ozone            !
!      cv    (IM)      : fraction of convective cloud                   !
!      cvt, cvb (IM)   : convective cloud top/bottom pressure in cb     !
!      iovrsw/iovrlw   : control flag for cloud overlap (sw/lw rad)     !
!                        =0 random overlapping clouds                   !
!                        =1 max/ran overlapping clouds                  !
!      fcice           : fraction of cloud ice  (in ferrier scheme)     !
!      frain           : fraction of rain water (in ferrier scheme)     !
!      rrime           : mass ratio of total to unrimed ice ( >= 1 )    !
!      flgmin          : minimim large ice fraction                     !
!      icsdsw/icsdlw   : auxiliary cloud control arrays passed to main  !
!           (IM)         radiations. if isubcsw/isubclw (input to init) !
!                        are set to 2, the arrays contains provided     !
!                        random seeds for sub-column clouds generators  !
!      np3d            : =3 brad ferrier microphysics scheme            !
!                        =4 zhao/carr/sundqvist microphysics scheme     !
!                        =5 external microphysics scheme provided bulk/ !
!                           grey quantities of cloud fields. (nmmb vars)!
!                        =8 Thompson microphysics scheme                !
!      ntoz            : =0 climatological ozone profile                !
!                        >0 interactive ozone profile                   !
!      NTRAC           : dimension veriable for array oz                !
!      NFXR            : second dimension of input/output array fluxr   !
!      dtlw, dtsw      : time duration for lw/sw radiation call in sec  !
!      lsswr, lslwr    : logical flags for sw/lw radiation calls        !
!      lssav           : logical flag for store 3-d cloud field         !
!      sashal          : logical flag for Jongil's shallow convection   !
!      norad_precip    : logical flag for not using precip in radiation !
!      crick_proof     : logical flag for eliminating CRICK             !
!      ccnorm          : logical flag for incloud condensate mixing ratio!
!      IX,IM           : horizontal dimention and num of used points    !
!      LM              : vertical layer dimension                       !
!      iflip           : control flag for in/out vertical indexing      !
!                        =0 index from toa to surface                   !
!                        =1 index from surface to toa                   !
!      me              : control flag for parallel process              !
!      lprnt           : control flag for diagnostic print out          !
!      ipt             : index for diagnostic printout point            !
!      kdt             : time-step number                               !
!                                                                       !
!    output variables:                                                  !
!      htrsw (IX,LM)   : total sky sw heating rate in k/sec             !
!      topfsw(IM)      : sw radiation fluxes at toa, components:        !
!                      (check module_radsw_parameters for definition)   !
!       %upfxc           - total sky upward sw flux at toa (w/m**2)     !
!       %dnflx           - total sky downward sw flux at toa (w/m**2)   !
!       %upfx0           - clear sky upward sw flux at toa (w/m**2)     !
!      sfcfsw(IM)      : sw radiation fluxes at sfc, components:        !
!                      (check module_radsw_parameters for definition)   !
!       %upfxc           - total sky upward sw flux at sfc (w/m**2)     !
!       %dnfxc           - total sky downward sw flux at sfc (w/m**2)   !
!       %upfx0           - clear sky upward sw flux at sfc (w/m**2)     !
!       %dnfx0           - clear sky downward sw flux at sfc (w/m**2)   !
!      sfalb (IM)      : mean surface diffused sw albedo                !
!      cldcov(IX,LM)   : 3-d cloud fraction                             !
!      cldsa(IX,5)     : fraction of clouds for low, mid, hi, tot, bl   !
!      htrlw (IX,LM)   : total sky lw heating rate in k/sec             !
!      topflw(IM)      : lw radiation fluxes at top, component:         !
!                        (check module_radlw_paramters for definition)  !
!       %upfxc           - total sky upward lw flux at toa (w/m**2)     !
!       %upfx0           - clear sky upward lw flux at toa (w/m**2)     !
!      sfcflw(IM)      : lw radiation fluxes at sfc, component:         !
!                        (check module_radlw_paramters for definition)  !
!       %upfxc           - total sky upward lw flux at sfc (w/m**2)     !
!       %upfx0           - clear sky upward lw flux at sfc (w/m**2)     !
!       %dnfxc           - total sky downward lw flux at sfc (w/m**2)   !
!       %dnfx0           - clear sky downward lw flux at sfc (w/m**2)   !
!      semis (IM)      : surface lw emissivity in fraction              !
!      tsflw (IM)      : surface air temp during lw calculation in k    !
!                                                                       !
!    input and output variables:                                        !
!      fluxr (IX,NFXR) : to save time accumulated 2-d fields defined as:!
!                 1      - toa total sky upwd lw radiation flux         !
!                 2      - toa total sky upwd sw radiation flux         !
!                 3      - sfc total sky upwd sw radiation flux         !
!                 4      - sfc total sky dnwd sw radiation flux         !
!                 5      - high domain cloud fraction                   !
!                 6      - mid  domain cloud fraction                   !
!                 7      - low  domain cloud fraction                   !
!                 8      - high domain mean cloud top pressure          !
!                 9      - mid  domain mean cloud top pressure          !
!                10      - low  domain mean cloud top pressure          !
!                11      - high domain mean cloud base pressure         !
!                12      - mid  domain mean cloud base pressure         !
!                13      - low  domain mean cloud base pressure         !
!                14      - high domain mean cloud top temperature       !
!                15      - mid  domain mean cloud top temperature       !
!                16      - low  domain mean cloud top temperature       !
!                17      - total cloud fraction                         !
!                18      - boundary layer domain cloud fraction         !
!                19      - sfc total sky dnwd lw radiation flux         !
!                20      - sfc total sky upwd lw radiation flux         !
!                21      - sfc total sky dnwd sw uv-b radiation flux    !
!                22      - sfc clear sky dnwd sw uv-b radiation flux    !
!                23      - toa incoming solar radiation flux            !
!                24      - sfc vis beam dnwd sw radiation flux          !
!                25      - sfc vis diff dnwd sw radiation flux          !
!                26      - sfc nir beam dnwd sw radiation flux          !
!                27      - sfc nir diff dnwd sw radiation flux          !
!                28      - toa clear sky upwd lw radiation flux         !
!                29      - toa clear sky upwd sw radiation flux         !
!                30      - sfc clear sky dnwd lw radiation flux         !
!                31      - sfc clear sky upwd sw radiation flux         !
!                32      - sfc clear sky dnwd sw radiation flux         !
!                33      - sfc clear sky upwd lw radiation flux         !
!optional        34      - aeros opt depth at 550nm (all components)    !
!               ....     - optional for test and future use             !
!                                                                       !
!    optional output variables:                                         !
!      htrswb(IX,LM,NBDSW) : spectral band total sky sw heating rate    !
!      htrlwb(IX,LM,NBDLW) : spectral band total sky lw heating rate    !
!                                                                       !
!                                                                       !
!    definitions of internal variable arrays:                           !
!                                                                       !
!     1. fixed gases:         (defined in 'module_radiation_gases')     !
!          gasvmr(:,:,1)  -  co2 volume mixing ratio                    !
!          gasvmr(:,:,2)  -  n2o volume mixing ratio                    !
!          gasvmr(:,:,3)  -  ch4 volume mixing ratio                    !
!          gasvmr(:,:,4)  -  o2  volume mixing ratio                    !
!          gasvmr(:,:,5)  -  co  volume mixing ratio                    !
!          gasvmr(:,:,6)  -  cf11 volume mixing ratio                   !
!          gasvmr(:,:,7)  -  cf12 volume mixing ratio                   !
!          gasvmr(:,:,8)  -  cf22 volume mixing ratio                   !
!          gasvmr(:,:,9)  -  ccl4 volume mixing ratio                   !
!                                                                       !
!     2. cloud profiles:      (defined in 'module_radiation_clouds')    !
!                ---  for  prognostic cloud  ---                        !
!          clouds(:,:,1)  -  layer total cloud fraction                 !
!          clouds(:,:,2)  -  layer cloud liq water path                 !
!          clouds(:,:,3)  -  mean effective radius for liquid cloud     !
!          clouds(:,:,4)  -  layer cloud ice water path                 !
!          clouds(:,:,5)  -  mean effective radius for ice cloud        !
!          clouds(:,:,6)  -  layer rain drop water path                 !
!          clouds(:,:,7)  -  mean effective radius for rain drop        !
!          clouds(:,:,8)  -  layer snow flake water path                !
!          clouds(:,:,9)  -  mean effective radius for snow flake       !
!                ---  for  diagnostic cloud  ---                        !
!          clouds(:,:,1)  -  layer total cloud fraction                 !
!          clouds(:,:,2)  -  layer cloud optical depth                  !
!          clouds(:,:,3)  -  layer cloud single scattering albedo       !
!          clouds(:,:,4)  -  layer cloud asymmetry factor               !
!                                                                       !
!     3. surface albedo:      (defined in 'module_radiation_surface')   !
!          sfcalb( :,1 )  -  near ir direct beam albedo                 !
!          sfcalb( :,2 )  -  near ir diffused albedo                    !
!          sfcalb( :,3 )  -  uv+vis direct beam albedo                  !
!          sfcalb( :,4 )  -  uv+vis diffused albedo                     !
!                                                                       !
!     4. sw aerosol profiles: (defined in 'module_radiation_aerosols')  !
!          faersw(:,:,:,1)-  sw aerosols optical depth                  !
!          faersw(:,:,:,2)-  sw aerosols single scattering albedo       !
!          faersw(:,:,:,3)-  sw aerosols asymmetry parameter            !
!                                                                       !
!     5. lw aerosol profiles: (defined in 'module_radiation_aerosols')  !
!          faerlw(:,:,:,1)-  lw aerosols optical depth                  !
!          faerlw(:,:,:,2)-  lw aerosols single scattering albedo       !
!          faerlw(:,:,:,3)-  lw aerosols asymmetry parameter            !
!                                                                       !
!     6. sw fluxes at toa:    (defined in 'module_radsw_main')          !
!        (topfsw_type -- derived data type for toa rad fluxes)          !
!          topfsw(:)%upfxc  -  total sky upward flux at toa             !
!          topfsw(:)%dnfxc  -  total sky downward flux at toa           !
!          topfsw(:)%upfx0  -  clear sky upward flux at toa             !
!                                                                       !
!     7. lw fluxes at toa:    (defined in 'module_radlw_main')          !
!        (topflw_type -- derived data type for toa rad fluxes)          !
!          topflw(:)%upfxc  -  total sky upward flux at toa             !
!          topflw(:)%upfx0  -  clear sky upward flux at toa             !
!                                                                       !
!     8. sw fluxes at sfc:    (defined in 'module_radsw_main')          !
!        (sfcfsw_type -- derived data type for sfc rad fluxes)          !
!          sfcfsw(:)%upfxc  -  total sky upward flux at sfc             !
!          sfcfsw(:)%dnfxc  -  total sky downward flux at sfc           !
!          sfcfsw(:)%upfx0  -  clear sky upward flux at sfc             !
!          sfcfsw(:)%dnfx0  -  clear sky downward flux at sfc           !
!                                                                       !
!     9. lw fluxes at sfc:    (defined in 'module_radlw_main')          !
!        (sfcflw_type -- derived data type for sfc rad fluxes)          !
!          sfcflw(:)%upfxc  -  total sky upward flux at sfc             !
!          sfcflw(:)%dnfxc  -  total sky downward flux at sfc           !
!          sfcflw(:)%dnfx0  -  clear sky downward flux at sfc           !
!                                                                       !
!! optional radiation outputs:                                          !
!!   10. sw flux profiles:    (defined in 'module_radsw_main')          !
!!       (profsw_type -- derived data type for rad vertical profiles)   !
!!         fswprf(:,:)%upfxc - total sky upward flux                    !
!!         fswprf(:,:)%dnfxc - total sky downward flux                  !
!!         fswprf(:,:)%upfx0 - clear sky upward flux                    !
!!         fswprf(:,:)%dnfx0 - clear sky downward flux                  !
!!                                                                      !
!!   11. lw flux profiles:    (defined in 'module_radlw_main')          !
!!       (proflw_type -- derived data type for rad vertical profiles)   !
!!         flwprf(:,:)%upfxc - total sky upward flux                    !
!!         flwprf(:,:)%dnfxc - total sky downward flux                  !
!!         flwprf(:,:)%upfx0 - clear sky upward flux                    !
!!         flwprf(:,:)%dnfx0 - clear sky downward flux                  !
!!                                                                      !
!!   12. sw sfc components:   (defined in 'module_radsw_main')          !
!!       (cmpfsw_type -- derived data type for component sfc fluxes)    !
!!         scmpsw(:)%uvbfc  -  total sky downward uv-b flux at sfc      !
!!         scmpsw(:)%uvbf0  -  clear sky downward uv-b flux at sfc      !
!!         scmpsw(:)%nirbm  -  total sky sfc downward nir direct flux   !
!!         scmpsw(:)%nirdf  -  total sky sfc downward nir diffused flux !
!!         scmpsw(:)%visbm  -  total sky sfc downward uv+vis direct flx !
!!         scmpsw(:)%visdf  -  total sky sfc downward uv+vis diff flux  !
!                                                                       !
!  ======================  end of definations  =======================  !
!
      implicit none
 
!  ---  constant parameter
!       The following need to be commented out after the newer version 
!       (Hsin-mu Lin, 2011-04-25)

 !     integer, parameter :: NSPC = 6


!  ---  inputs: (for rank>1 arrays, horizontal dimensioned by IX)
      integer,  intent(in) :: IX,IM, LM, NTRAC, NFXR, me,               &
     &                        ntoz, ipt, kdt,                           &
     &                        nrads
      logical              :: isday(IX)  ! true = day, false = night
      integer,  intent(in) :: icsdsw(IX), icsdlw(IX), jdate(8)
      integer,  intent(in) :: mbota(IX,3), mtopa(IX,3)

      real (kind=kind_phys), intent(in) :: clouds1(IX,LM,NF_CLDS)

      logical,  intent(in) :: lsswr, lslwr, lssav, lprnt

      real (kind=kind_phys), dimension(IX,LM+1), intent(in) :: prsi,    &
     &       plvl1

      real (kind=kind_phys), dimension(IX,LM), intent(in) :: plyr1,     &
     &       prslk1, tgrs, qgrs, rhly1, qlyr1, tvly1

      real (kind=kind_phys), dimension(IX),      intent(in) ::  slmsk,  &
     &       xlon, xlat, tsfc, snowd, zorl, hprim,                      &
     &       fice, tisfc,                                               &
     &       sncovr, snoalb, sinlat, coslat

      real (kind=kind_phys), dimension(IX),    intent(inout) ::         &
     &       alvsf, alnsf, alvwf, alnwf, facsf, facwf

      real (kind=kind_phys), intent(in) :: solcon, dtlw, dtsw, solhr,   &
     &       tracer1(IX,LM,NTRAC)                                       &
     &     , dtswav, SALBEDO(IX), SM(IX) ! extra input

! --- outputs: (horizontal dimensioned by IX)
      real (kind=kind_phys), dimension(IX,LM),intent(out):: htrsw, htrlw

      real (kind=kind_phys), dimension(IX,5),intent(in):: cldsa    

      real (kind=kind_phys), dimension(IX),   intent(out):: tsflw,      &
     &       sfalb, semis, coszen, coszdg

!rv   real, intent(in) :: CPATHFAC4LW
      real(kind=kind_phys), intent(in) :: CPATHFAC4LW

      type (topfsw_type), dimension(IX), intent(out) :: topfsw
      type (sfcfsw_type), dimension(IX), intent(out) :: sfcfsw
      type (topflw_type), dimension(IX), intent(out) :: topflw
      type (sfcflw_type), dimension(IX), intent(out) :: sfcflw

!  ---  variables are for both input and output:
      real (kind=kind_phys), intent(inout) :: fluxr(IX,NFXR)

!! ---  optional outputs:
      real (kind=kind_phys), dimension(IX,LM,NBDSW), optional, &
     &                       intent(out) :: htrswb
      real (kind=kind_phys), dimension(IX,LM,NBDLW), optional, &
     &                       intent(out) :: htrlwb

      real (kind=kind_phys), dimension(im,LM+1+LTP)::plvl, tlvl

      real (kind=kind_phys), dimension(im,LM+LTP)  ::                   &
     &       plyr, tlyr, qlyr,                                          &
     &       olyr, rhly, vvel, prslk, tem2da, tem2db, tvly

      real (kind=kind_phys), dimension(im,LM+LTP)  ::                   &
     &       qc2d, qi2d, qs2d, ni2d

      real (kind=kind_phys), dimension(im) :: tsfa, sfcemis,    &
     &       tsfg, tskn, tem1d

      real (kind=kind_phys), dimension(im) :: CZMEAN

      real (kind=kind_phys), dimension(im,LM+LTP,NF_CLDS) :: clouds
      real (kind=kind_phys), dimension(im,LM+LTP,NF_CLDS) :: dummys
      real (kind=kind_phys), dimension(im,LM+LTP,NF_VGAS) :: gasvmr
      real (kind=kind_phys), dimension(im,       NF_ALBD) :: sfcalb
!     real (kind=kind_phys), dimension(im,       NSPC1)   :: aerodp      ! optn for aod output
      real (kind=kind_phys), dimension(im,LM+LTP,NTRAC)   :: tracer

      real (kind=kind_phys), dimension(im,LM+LTP,NBDSW,NF_AESW):: faersw
      real (kind=kind_phys), dimension(im,LM+LTP,NBDLW,NF_AELW):: faerlw

 !     real (kind=kind_phys), dimension(im,LM+LTP,NSPC-1) :: tau_gocart

      real (kind=kind_phys), dimension(im,LM+LTP) :: htswc
      real (kind=kind_phys), dimension(im,LM+LTP) :: htlwc

      real (kind=kind_phys), dimension(im,LM+LTP) :: gcice, grain, grime

!! ---  may be used for optional sw/lw outputs:
!!      take out "!!" as needed
      type (cmpfsw_type),    dimension(im)          :: scmpsw
      real (kind=kind_phys), dimension(im,LM+LTP,NBDSW) :: htswb

      real (kind=kind_phys), dimension(im,LM+LTP,NBDLW) :: htlwb

      real (kind=kind_phys)::raddt,qs, delt, tem0d                       &
     & ,SFCALBEDO(im), SMX(im)

      integer j2
      integer :: i, j, k, k1, lv, icec, itop, ibtc, nday, idxday(im),    &
     &      LP1, nb, LMK, LMP, kd, lla, llb, lya, lyb, kt, kb, np3d

!  ---  for debug test use
!     real (kind=kind_phys) :: temlon, temlat, alon, alat
!     integer :: ipt
!     logical :: lprnt1

!!==========================================================================
!  --- Albedo (ETA era) calculation from 2012 RRTM  (Lin, 20130225)
!
      INTEGER :: IQ,JX

      INTEGER,EXTERNAL :: omp_get_thread_num

 !     REAL (kind=kind_phys) :: ALBD0, ALVD1, ALND1
      REAL (kind=kind_phys) :: ZEN, DZEN, ALB1, ALB2

      REAL (kind=kind_phys), PARAMETER :: TWENTY=20.0,        &
                                          HP537=0.537,        &
                                          ONE=1.,             &
                                          DEGRAD1=180.0/PI,   &
                                          H74E1=74.0,         &
                                          HAF=0.5,            &
                                          HNINETY=90.,        &
                                          FIFTY=50.,          &
                                          QUARTR=0.25,        &
                                          HNINE=9.0,          &
                                          HP1=0.1,            &
                                          H15E1=15.0

      REAL (kind=kind_phys), PARAMETER :: coszenmin =1.0E-4

      REAL (kind=kind_phys), DIMENSION(IM) :: ALVB,ALNB,ALVD,ALND

      REAL (kind=kind_phys), DIMENSION(20) :: ZA
      REAL (kind=kind_phys), DIMENSION(19) :: DZA
      REAL (kind=kind_phys), DIMENSION(21) :: TRN
      REAL (kind=kind_phys), DIMENSION(21,20) :: ALBD

      DATA TRN/.00,.05,.10,.15,.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,  &
               .70,.75,.80,.85,.90,.95,1.00/

      DATA  ALBD/.061,.062,.072,.087,.115,.163,.235,.318,.395,.472,.542, &
       .604,.655,.693,.719,.732,.730,.681,.581,.453,.425,.061,.062,.070, &
       .083,.108,.145,.198,.263,.336,.415,.487,.547,.595,.631,.656,.670, &
       .652,.602,.494,.398,.370,.061,.061,.068,.079,.098,.130,.174,.228, &
       .290,.357,.424,.498,.556,.588,.603,.592,.556,.488,.393,.342,.325, &
       .061,.061,.065,.073,.086,.110,.150,.192,.248,.306,.360,.407,.444, &
       .469,.480,.474,.444,.386,.333,.301,.290,.061,.061,.065,.070,.082, &
       .101,.131,.168,.208,.252,.295,.331,.358,.375,.385,.377,.356,.320, &
       .288,.266,.255,.061,.061,.063,.068,.077,.092,.114,.143,.176,.210, &
       .242,.272,.288,.296,.300,.291,.273,.252,.237,.266,.220,.061,.061, &
       .062,.066,.072,.084,.103,.127,.151,.176,.198,.219,.236,.245,.250, &
       .246,.235,.222,.211,.205,.200,                                    &
                 .061,.061,.061,.065,.071,.079,.094,.113,.134,.154,.173, &
       .185,.190,.193,.193,.190,.188,.185,.182,.180,.178,.061,.061,.061, &
       .064,.067,.072,.083,.099,.117,.135,.150,.160,.164,.165,.164,.162, &
       .160,.159,.158,.157,.157,.061,.061,.061,.062,.065,.068,.074,.084, &
       .097,.111,.121,.127,.130,.131,.131,.130,.129,.127,.126,.125,.122, &
       .061,.061,.061,.061,.062,.064,.070,.076,.085,.094,.101,.105,.107, &
       .106,.103,.100,.097,.096,.095,.095,.095,.061,.061,.061,.060,.061, &
       .062,.065,.070,.075,.081,.086,.089,.090,.088,.084,.080,.077,.075, &
       .074,.074,.074,.061,.061,.060,.060,.060,.061,.063,.065,.068,.072, &
       .076,.077,.076,.074,.071,.067,.064,.062,.061,.061,.061,.061,.061, &
       .060,.060,.060,.060,.061,.062,.065,.068,.069,.069,.068,.065,.061, &
       .058,.055,.054,.053,.052,.052,                                    &
                 .061,.061,.060,.060,.060,.060,.060,.060,.062,.065,.065, &
       .063,.060,.057,.054,.050,.047,.046,.045,.044,.044,.061,.061,.060, &
       .060,.060,.059,.059,.059,.059,.059,.058,.055,.051,.047,.043,.039, &
       .035,.033,.032,.031,.031,.061,.061,.060,.060,.060,.059,.059,.058, &
       .057,.056,.054,.051,.047,.043,.039,.036,.033,.030,.028,.027,.026, &
       .061,.061,.060,.060,.060,.059,.059,.058,.057,.055,.052,.049,.045, &
       .040,.036,.032,.029,.027,.026,.025,.025,.061,.061,.060,.060,.060, &
       .059,.059,.058,.056,.053,.050,.046,.042,.038,.034,.031,.028,.026, &
       .025,.025,.025,.061,.061,.060,.060,.059,.058,.058,.057,.055,.053, &
       .050,.046,.042,.038,.034,.030,.028,.029,.025,.025,.025/
!
      DATA ZA/90.,88.,86.,84.,82.,80.,78.,76.,74.,70.,66.,62.,58.,54., &
              50.,40.,30.,20.,10.,0.0/
!
      DATA DZA/8*2.0,6*4.0,5*10.0/
!
!  --- end of Albedo (ETA era) calculation from 2012 RRTM
!!==========================================================================

!
!===> ...  begin here
!
      np3d = icmphys

      LP1 = LM + 1               ! num of in/out levels

!===========================================================================
!  --- ...  set local /level/layer indexes corresponding to in/out variables
!           ( for optional extra top layer)

      LMK = LM + LTP             ! num of local layer
      LMP = LMK + 1              ! num of local level

      if ( lextop ) then
        if ( ivflip == 1 ) then   ! vertical from sfc upward
          kd = 0                   ! index diff between in/out and local
          kt = 1                   ! index diff between lyr and upper bound
          kb = 0                   ! index diff between lyr and lower bound
          lla = LMK                ! local index at the 2nd level from top
          llb = LMP                ! local index at toa level
          lya = LM                 ! local index for the 2nd layer from top
          lyb = LP1                ! local index for the top layer
        else                     ! vertical from toa downward
          kd = 1                   ! index diff between in/out and local
          kt = 0                   ! index diff between lyr and upper bound
          kb = 1                   ! index diff between lyr and lower bound
          lla = 2                  ! local index at the 2nd level from top
          llb = 1                  ! local index at toa level
          lya = 2                  ! local index for the 2nd layer from top
          lyb = 1                  ! local index for the top layer
        endif                    ! end if_ivflip_block
      else
        kd = 0
        if ( ivflip == 1 ) then  ! vertical from sfc upward
          kt = 1                   ! index diff between lyr and upper bound
          kb = 0                   ! index diff between lyr and lower bound
        else                     ! vertical from toa downward
          kt = 0                   ! index diff between lyr and upper bound
          kb = 1                   ! index diff between lyr and lower bound
        endif                    ! end if_ivflip_block
      endif   ! end if_lextop_block
!
!  --- End of index setting for optional extra top layer
!===========================================================================

      raddt = min(dtsw, dtlw)

!  --- ...  for debug test
!     alon = 120.0
!     alat = 29.5
!     ipt = 0
!     do i = 1, IM
!       temlon = xlon(i) * 57.29578
!       if (temlon < 0.0) temlon = temlon + 360.0
!       temlat = xlat(i) * 57.29578
!       lprnt1 = abs(temlon-alon) < 1.1 .and. abs(temlat-alat) < 1.1
!       if ( lprnt1 ) then
!         ipt = i
!         exit
!       endif
!     enddo

!     print *,' in grrad : raddt=',raddt

!  --- ...  setup surface ground temp and ground/air skin temp if required

      if ( itsfc == 0 ) then            ! use same sfc skin-air/ground temp
        do i = 1, IM
          tskn(i) = tsfc(i)
          tsfg(i) = tsfc(i)
        enddo
      else                              ! use diff sfc skin-air/ground temp
        do i = 1, IM
!!        tskn(i) = ta  (i)               ! not yet
!!        tsfg(i) = tg  (i)               ! not yet
          tskn(i) = tsfc(i)
          tsfg(i) = tsfc(i)
        enddo
      endif

!=============================================================

      do k = 1, LM
        k1 = k + kd      ! kd: index diff between in/out and local

!DEC$ SIMD
        do i = 1, IM
          plvl(i,k1) = plvl1(i,k)
          plyr(i,k1) = plyr1(i,k)
          tlyr(i,k1) = tgrs(i,k)
          prslk(i,k1) = prslk1(i,k)
        ENDDO
!DEC$ SIMD

        DO i = 1, IM
          rhly(i,k1) = rhly1(i,k)
        enddo

        do j = 1, NF_CLDS
          do i = 1, IM
             clouds(i,k1,j) = clouds1(i,k,j)
          enddo
        enddo

        do j = 1, NTRAC
          do i = 1, IM
             tracer(i,k1,j) = tracer1(i,k,j)
          enddo
        enddo
      enddo

      do i = 1, IM
        plvl(i,LP1+kd) = plvl1(i,LP1)
      enddo

!======================================================================
!  --- ...  values for extra top layer

      if ( lextop ) then                 ! values for extra top layer
!DEC$ SIMD
        do i = 1, IM
          plvl(i,llb) = prsmin
          if ( plvl(i,lla)<=prsmin ) plvl(i,lla) = 2.0*prsmin
          plyr(i,lyb) = 0.5 * plvl(i,lla)
          tlyr(i,lyb) = tlyr(i,lya)
          rhly(i,lyb) = rhly(i,lya)
        enddo

! loop won't vectorize because of call to fpkapx
        do i = 1, IM
          prslk(i,lyb) = fpkapx(plyr(i,lyb)*100.0) ! fpkapx in pa
        enddo

!  ---  note: may need to take care the top layer mount

         do i = 1, IM
           tracer(i,lyb,1) = tracer(i,lya,1)    ! 1st tracer element is ozone
         enddo

         do j = 2, NTRAC     ! This loop variable should now be correct
           do i = 1, IM
             tracer(i,lyb,j) = 0.       !  enforce no cloud above domain top
           enddo
         enddo

         if (np3d == 5) then         ! nmmb ferrier+bmj  (GFDL type diagnostic)
           do i = 1, IM
             clouds(i,lyb,1) = clouds(i,lya,1)
             clouds(i,lyb,2) = clouds(i,lya,2)
             clouds(i,lyb,3) = clouds(i,lya,3)
             clouds(i,lyb,4) = clouds(i,lya,4)
           enddo
         else
           do j = 1, NF_CLDS
             do i = 1, IM
               clouds(i,lyb,j) = 0.
             enddo
           enddo
         endif

      endif
!
!======================================================================

!  --- ...  get layer ozone mass mixing ratio

      if (ntoz > 0) then            ! interactive ozone generation

        do k = 1, LMK
          do i = 1, IM
            olyr(i,k) = max( QMIN, tracer(i,k,ntoz) )
          enddo
        enddo

      else                          ! climatological ozone

!     print *,' in grrad : calling getozn'

        call getozn                                                     &
!  ---  inputs:
     &     ( prslk,xlat,                                                &
     &       IM, LMK,                                                   &
!  ---  outputs:
     &       olyr                                                       &
     &     )

      endif                            ! end_if_ntoz

!  --- ...  compute cosin of zenith angle

      call coszmn_nmmb                                                  &
!  ---  inputs:
     &     ( xlon,sinlat,coslat,solhr,IM, me,                           &
     &       dtswav,nrads ,                                             &  ! Extra input
!  ---  outputs:
     &       coszen, coszdg                                             &
     &      )

!  --- ...  set up non-prognostic gas volume mixing ratioes

      call getgases                                                     &
!  ---  inputs:
     &    ( plvl, xlon, xlat,                                           &
     &      IM, LMK,                                                    &
!  ---  outputs:
     &      gasvmr                                                      &
     &     )

!  --- ...  get temperature at layer interface, and layer moisture

      do k = 2, LMK
        do i = 1, IM
          tem2da(i,k) = log( plyr(i,k) )
          tem2db(i,k) = log( plvl(i,k) )
        enddo
      enddo

      if (ivflip == 0) then               ! input data from toa to sfc

        do i = 1, IM
          tem1d (i)   = QME6
          tem2da(i,1) = log( plyr(i,1) )

        !---------------------------------------------
        ! take care for the model top (20130410, Lin)
        !---------------------------------------------

          if ( plvl(i,1) == 0.0 ) then
             tem2db(i,1) = -14.0
          else
             tem2db(i,1) = log( plvl(i,1) )
          endif

         ! tem2db(i,1) = 1.0

          tsfa  (i)   = tlyr(i,LMK)                  ! sfc layer air temp
          tlvl(i,1)   = tlyr(i,1)
          tlvl(i,LMP) = tskn(i)
        enddo

        do k = 1, LM
          k1 = k + kd      ! kd: index diff between in/out and local

          do i = 1, IM
            qlyr(i,k1) = max( tem1d(i), qgrs(i,k) )
            tem1d(i)   = min( QME5, qlyr(i,k1) )
            tvly(i,k1) = tgrs(i,k) * (1.0 + con_fvirt*qlyr(i,k1))          ! virtual temp in K

           ! qlyr(i,k1) = qlyr1(i,k)
           ! tvly(i,k1) = tvly1(i,k)          ! virtual temp in K
          enddo
        enddo

        !============================================
        !  --- ...  values for extra top layer

        if ( lextop ) then
          do i = 1, IM
            qlyr(i,lyb) = qlyr(i,lya)
            tvly(i,lyb) = tvly(i,lya)
          enddo
        endif

        !============================================

        do k = 2, LMK
          do i = 1, IM
            tlvl(i,k) = tlyr(i,k) + (tlyr(i,k-1) - tlyr(i,k))           &
     &                * (tem2db(i,k)   - tem2da(i,k))                   &
     &                / (tem2da(i,k-1) - tem2da(i,k))
          enddo
        enddo

      else                               ! input data from sfc to toa

        do i = 1, IM
          tem1d (i)   = QME6
          tem2da(i,1) = log( plyr(i,1) )
          tem2db(i,1) = log( plvl(i,1) )
          tsfa  (i)   = tlyr(i,1)                    ! sfc layer air temp
          tlvl(i,1)   = tskn(i)
          tlvl(i,LMP) = tlyr(i,LMK)
        enddo

        do k = LM, 1, -1
          do i = 1, IM
            qlyr(i,k) = max( tem1d(i), qgrs(i,k) )
            tem1d(i)  = min( QME5, qlyr(i,k) )
            tvly(i,k) = tgrs(i,k) * (1.0 + con_fvirt*qlyr(i,k))          ! virtual temp in K

           ! qlyr(i,k) = qlyr1(i,k)
           ! tvly(i,k) = tvly1(i,k)        ! virtual temp in K
          enddo
        enddo

        !============================================
        !  --- ...  values for extra top layer

        if ( lextop ) then
          do i = 1, IM
            qlyr(i,lyb) = qlyr(i,lya)
            tvly(i,lyb) = tvly(i,lya)
          enddo
        endif

        !============================================

        do k = 1, LMK-1
          do i = 1, IM
            tlvl(i,k+1) = tlyr(i,k) + (tlyr(i,k+1) - tlyr(i,k))         &
     &                  * (tem2db(i,k+1) - tem2da(i,k))                 &
     &                  / (tem2da(i,k+1) - tem2da(i,k))
          enddo
        enddo

      endif                              ! end_if_ivflip


!  --- ...  setup aerosols property profile for radiation

!check  print *,' in grrad : calling setaer '

      call setaer                                                     &
!  ---  inputs:
     &   ( plvl,plyr,prslk,tvly,rhly,slmsk,tracer,xlon,xlat,          &
     &     IM,LMK,LMP, lsswr,lslwr,                                   &
!  ---   extra 2 vars for old code that has been replaced by 'tvly'
!        can be eliminated after aggrements being reached
     &     tlyr,qlyr,                                                 &
!  ---  outputs:
     &     faersw,faerlw                                              &
 !   &,    tau_gocart                                                 &
     &   )


!==========================================================================
!  Albedo (ETA era) calculation from 2012 RRTM  (Lin, 20130225)
!==========================================================================
      IF (ialbflg == 2) THEN

         do i = 1, IM

            SMX(i) = SM(i)
            SFCALBEDO(i) = SALBEDO(i)

!..... THE FOLLOWING CODE GETS ALBEDO FROM PAYNE,1972 TABLES IF
!         1) OPEN SEA POINT (SLMSK=1);  2) KALB=0

            IQ=INT(TWENTY*HP537+ONE)
            CZMEAN(i)=coszen(i)

            IF(CZMEAN(i).GT.0.0 .AND. SMX(i).GT.0.5) THEN
               ZEN=DEGRAD1*ACOS(MAX(CZMEAN(i),0.0))

               IF(ZEN.GE.H74E1) JX=INT(HAF*(HNINETY-ZEN)+ONE)
               IF(ZEN.LT.H74E1 .AND. ZEN.GE.FIFTY) &
                 JX=INT(QUARTR*(H74E1-ZEN)+HNINE)
               IF(ZEN .LT. FIFTY) JX=INT(HP1*(FIFTY-ZEN)+H15E1)

               DZEN=-(ZEN-ZA(JX))/DZA(JX)

               ALB1=ALBD(IQ,JX)+DZEN*(ALBD(IQ,JX+1)-ALBD(IQ,JX))
               ALB2=ALBD(IQ+1,JX)+DZEN*(ALBD(IQ+1,JX+1)-ALBD(IQ+1,JX))

               SFCALBEDO(I)=ALB1+TWENTY*(ALB2-ALB1)*(HP537-TRN(IQ)) ! BSF
            ENDIF

!.....     VISIBLE AND NEAR IR DIFFUSE ALBEDO
            ALVD(I) = SFCALBEDO(I) ! BSF
            ALND(I) = SFCALBEDO(I) ! BSF

!.....     VISIBLE AND NEAR IR DIRECT BEAM ALBEDO
            ALVB(I) = SFCALBEDO(I) ! BSF
            ALNB(I) = SFCALBEDO(I) ! BSF
!
!--- Remove diurnal variation of land surface albedos (Ferrier, 6/28/05)
!--- Turn back on to mimic NAM 8/17/05
!
!.....     VISIBLE AND NEAR IR DIRECT BEAM ALBEDO, IF NOT OCEAN NOR SNOW
!        ..FUNCTION OF COSINE SOLAR ZENITH ANGLE..
!
!=== The following line commented out by the "fixed" are used instead of the
!    original version (09/2012)
!
!fixed            IF (SMX .LT. 0.5) THEN
!fixed             IF (SFCALBEDO .LE. 0.5) THEN
!fixed             ALBD0=-18.0 * (0.5 - ACOS(CZMEAN(I))/PI)
!fixed             ALBD0=EXP (ALBD0)
!fixed             ALVD1=(ALVD(I) - 0.054313) / 0.945687
!fixed             ALND1=(ALND(I) - 0.054313) / 0.945687
!fixed             ALVB(I)=ALVD1 + (1.0 - ALVD1) * ALBD0
!fixed             ALNB(I)=ALND1 + (1.0 - ALND1) * ALBD0
!fixed! !-- Put in an upper limit on beam albedos
!fixed             ALVB(I)=MIN(0.5,ALVB(I))
!fixed             ALNB(I)=MIN(0.5,ALNB(I))
!fixed            ENDIF
!fixed           ENDIF

!!! WE INTRODUCE HERE DIRECT AND DIFFUSE ALBEDO... FOR THIS OPTION, THERE IS A CHANGE IN GRRAD.f

            alvsf(i) = ALVB(I) !For this option ALVSF1 is direct visible albedo
            alnsf(i) = ALNB(I) !For this option ALNSF1 is direct nir albedo
            alvwf(i) = ALVD(I) !For this option ALVWF1 is diffuse visible albedo
            alnwf(i) = ALND(I) !For this option ALNWF1 is diffuse nir albedo
            facsf(i) = 0.      !not used with this option
            facwf(i) = 0.      !not used for this option

         ENDDO
      ENDIF    !---- end of IALB=2  GFDL TYPE RADIATION

!
! --- End of Albedo (ETA era) calculation from 2012 RRTM ------------------
!==========================================================================



!  --- ...  start radiation calculations 
!           remember to set heating rate unit to k/sec!

      if (lsswr) then

!  ---  setup surface albedo for sw radiation, incl xw (nov04) sea-ice

        call setalb                                                     &
!  ---  inputs:
     &     ( slmsk,snowd,sncovr,snoalb,zorl,coszen,tsfg,tsfa,hprim,     &
     &       alvsf,alnsf,alvwf,alnwf,facsf,facwf,fice,tisfc,            &
     &       IM,                                                        &
!  ---  outputs:
     &       sfcalb                                                     &
     &     )

!  --- lu [+4L]: derive SFALB from vis- and nir- diffuse surface albedo

        do i = 1, IM
          sfalb(i) = max(0.01, 0.5 * (sfcalb(i,2) + sfcalb(i,4)))
        enddo

!  ---  check for daytime points

        nday = 0

        do i = 1, IM
          isday(i) = (coszen(i) >= coszenmin)    ! day
          if ( isday(i) ) then
            nday = nday + 1
            idxday(nday) = i
          endif
        enddo

! ===+===================================================================
!
        if ( nday > 0 ) then

!----------------------------------------------------------------

!     print *,' in grrad : calling swrad'

          if ( present(htrswb) ) then

            call swrad                                                  &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,icsdsw,faersw,sfcalb,                               &
     &       coszen,solcon, nday,idxday,                                &
     &       IM, LMK, LMP, lprnt,                                       &
!  ---  outputs:
     &       htswc,topfsw,sfcfsw                                        &
!! ---  optional:
!!   &,      HSW0=htsw0,FLXPRF=fswprf                                   &
     &,      HSWB=htswb,FDNCMP=scmpsw                                   &
     &     )

            do k = 1, LM
              k1 = k + kd      ! kd: index diff between in/out and local

              do j = 1, NBDSW
                do i = 1, IM
                  if ( isday(i) ) then
                    htrswb(i,k,j) = htswb(i,k1,j)
                  else
                    htrswb(i,k,j) = 0.0     ! night time ==> 0
                  endif
                enddo
              enddo
            enddo

          else

            call swrad                                                  &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,icsdsw,faersw,sfcalb,                               &
     &       coszen,solcon, nday,idxday,                                &
     &       IM, LMK, LMP, lprnt,                                       &
!  ---  outputs:
     &       htswc,topfsw,sfcfsw                                        &
!! ---  optional:
!!   &,      HSW0=htsw0,FLXPRF=fswprf,HSWB=htswb                        &
     &,      FDNCMP=scmpsw                                              &
     &     )

          endif

          do k = 1, LM
            k1 = k + kd      ! kd: index diff between in/out and local

            do i = 1, IM
              if ( isday(i) ) then
                htrsw(i,k) = htswc(i,k1)
              else
                htrsw(i,k) = 0.0      ! night time ==> 0
              endif
            enddo
          enddo

          !---  night time ==> 0

          do i = 1, IM
            if ( .not.  isday(i) ) then
               sfcfsw(i) = sfcfsw_type( 0.0, 0.0, 0.0, 0.0 )
               topfsw(i) = topfsw_type( 0.0, 0.0, 0.0 )
               scmpsw(i) = cmpfsw_type( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )
            endif
          enddo

        else                   ! if_nday_block, for all night section

          do k = 1, LM
            do i = 1, IM
              htrsw(i,k) = 0.0
            enddo
          enddo

          sfcfsw = sfcfsw_type( 0.0, 0.0, 0.0, 0.0 )
          topfsw = topfsw_type( 0.0, 0.0, 0.0 )
          scmpsw = cmpfsw_type( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )

!! ---  optional:
!!        fswprf= profsw_type( 0.0, 0.0, 0.0, 0.0 )

          if ( present(htrswb) ) then
            do j = 1, NBDSW
              do k = 1, LM
                do i = 1, IM
                  htrswb(i,k,j) = 0.0
                enddo
              enddo
            enddo
          endif

        endif                  ! end_if_nday

      endif                                ! end_if_lsswr

! ======= end of SW process ==================================


! ======= process LW =========================================

      if (lslwr) then

!  ---  setup surface emissivity for lw radiation

        call setemis                                                    &
!  ---  inputs:
     &     ( xlon,xlat,slmsk,snowd,sncovr,zorl,tsfg,tsfa,hprim,         &
     &       IM,                                                        &
!  ---  outputs:
     &       sfcemis                                                    &
     &     )

!     print *,' in grrad : calling lwrad'

!==========================================================================
!  Special for lwrad to enhence the emissivity
!  it is similar to *CPATHFAC4LW to odcld in radlw (Hsin-Mu Lin, 20140520)
!==========================================================================

        if (np3d == 3) then          ! ferrier's microphysics
          do k = 1, LMK
            do i = 1, IM
              dummys(i,k,2) = clouds(i,k,2)
              dummys(i,k,4) = clouds(i,k,4)
              clouds(i,k,2) = dummys(i,k,2)*CPATHFAC4LW  ! cloud liquid path 
              clouds(i,k,4) = dummys(i,k,4)*CPATHFAC4LW  ! cloud ice path
            enddo
          enddo
        endif

!==========================================================================

        if ( present(htrlwb) ) then

          call lwrad                                                    &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,icsdlw,faerlw,sfcemis,tsfg,                         &
     &       IM, LMK, LMP, lprnt,                                       &
!  ---  outputs:
     &       htlwc,topflw,sfcflw                                        &
!! ---  optional:
!!   &,      HLW0=htlw0,FLXPRF=flwprf                                   &
     &,      HLWB=htlwb                                                 &
     &     )

          do k = 1, LM
            k1 = k + kd      ! kd: index diff between in/out and local

            do j = 1, NBDLW
              do i = 1, IM
                htrlwb(i,k,j) = htlwb(i,k1,j)
              enddo
            enddo
          enddo

        else

          call lwrad                                                    &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,icsdlw,faerlw,sfcemis,tsfg,                         &
     &       IM, LMK, LMP, lprnt,                                       &
!  ---  outputs 
     &       htlwc,topflw,sfcflw                                        &
!! ---  optional:
!!   &,      HLW0=htlw0,FLXPRF=flwprf,HLWB=htlwb                        &
     &     )

        endif

        do i = 1, IM
          semis (i) = sfcemis(i)
!  ---  save surface air temp for diurnal adjustment at model t-steps
          tsflw (i) = tsfa(i)
        enddo

        do k = 1, LM
          k1 = k + kd      ! kd: index diff between in/out and local

          do i = 1, IM
            htrlw(i,k) = htlwc(i,k1)
          enddo
        enddo

!==========================================================================
!   return *1.5 clouds to original value  (Hsin-Mu Lin, 20140520)
!==========================================================================

        if (np3d == 3) then          ! ferrier's microphysics
          do k = 1, LMK
            do i = 1, IM
              clouds(i,k,2) = dummys(i,k,2)  ! cloud liquid path
              clouds(i,k,4) = dummys(i,k,4)  ! cloud ice path
            enddo
          enddo
        endif

! ======= end of LW process ==================================

      endif                                ! end_if_lslwr


!==========================================================================

!  --- ...  collect the fluxr data for wrtsfc

      if (lssav) then

!  ---  save lw toa and sfc fluxes

        if (lslwr) then
          do i = 1, IM
!  ---  lw total-sky fluxes
            fluxr(i,1 ) = fluxr(i,1 ) + dtlw * topflw(i)%upfxc   ! total sky top lw up
            fluxr(i,19) = fluxr(i,19) + dtlw * sfcflw(i)%dnfxc   ! total sky sfc lw dn
            fluxr(i,20) = fluxr(i,20) + dtlw * sfcflw(i)%upfxc   ! total sky sfc lw up
!  ---  lw clear-sky fluxes
            fluxr(i,28) = fluxr(i,28) + dtlw * topflw(i)%upfx0   ! clear sky top lw up
            fluxr(i,30) = fluxr(i,30) + dtlw * sfcflw(i)%dnfx0   ! clear sky sfc lw dn
            fluxr(i,33) = fluxr(i,33) + dtlw * sfcflw(i)%upfx0   ! clear sky sfc lw up
          enddo
        endif

!  ---  save sw toa and sfc fluxes with proper diurnal sw wgt. coszen=mean cosz over daylight
!       part of sw calling interval, while coszdg= mean cosz over entire interval

        if (lsswr) then
          do i = 1, IM
            if (coszen(i) > 0.) then
!  ---  sw total-sky fluxes
              tem0d = dtsw * coszdg(i) / coszen(i)
              fluxr(i,2 ) = fluxr(i,2)  + topfsw(i)%upfxc * tem0d  ! total sky top sw up
              fluxr(i,3 ) = fluxr(i,3)  + sfcfsw(i)%upfxc * tem0d  ! total sky sfc sw up
              fluxr(i,4 ) = fluxr(i,4)  + sfcfsw(i)%dnfxc * tem0d  ! total sky sfc sw dn
!  ---  sw uv-b fluxes
              fluxr(i,21) = fluxr(i,21) + scmpsw(i)%uvbfc * tem0d  ! total sky uv-b sw dn
              fluxr(i,22) = fluxr(i,22) + scmpsw(i)%uvbf0 * tem0d  ! clear sky uv-b sw dn
!  ---  sw toa incoming fluxes
              fluxr(i,23) = fluxr(i,23) + topfsw(i)%dnfxc * tem0d  ! top sw dn
!  ---  sw sfc flux components
              fluxr(i,24) = fluxr(i,24) + scmpsw(i)%visbm * tem0d  ! uv/vis beam sw dn
              fluxr(i,25) = fluxr(i,25) + scmpsw(i)%visdf * tem0d  ! uv/vis diff sw dn
              fluxr(i,26) = fluxr(i,26) + scmpsw(i)%nirbm * tem0d  ! nir beam sw dn        !?
              fluxr(i,27) = fluxr(i,27) + scmpsw(i)%nirdf * tem0d  ! nir diff sw dn        !?
!  ---  sw clear-sky fluxes
              fluxr(i,29) = fluxr(i,29) + topfsw(i)%upfx0 * tem0d  ! clear sky top sw up
              fluxr(i,31) = fluxr(i,31) + sfcfsw(i)%upfx0 * tem0d  ! clear sky sfc sw up
              fluxr(i,32) = fluxr(i,32) + sfcfsw(i)%dnfx0 * tem0d  ! clear sky sfc sw dn
            endif
          enddo
        endif

!  ---  save total cloud and bl cloud

        if (lsswr .or. lslwr) then
          do i = 1, IM
            fluxr(i,17) = fluxr(i,17) + raddt * cldsa(i,4)             !?
            fluxr(i,18) = fluxr(i,18) + raddt * cldsa(i,5)             !?
          enddo

!  ---  save cld frac,toplyr,botlyr and top temp, note that the order
!       of h,m,l cloud is reversed for the fluxr output.
!  ---  save interface pressure (Pa) of top/bot

          do j = 1, 3
            do i = 1, IM
              tem0d = raddt * cldsa(i,j)
              itop  = mtopa(i,j) - kd
              ibtc  = mbota(i,j) - kd

            !==========================================
            !=== Zavisa mentioned on 01/15/2013
            !    for our current run, there are no problem for np3d=3
            !    without this modification, while np3d=5 need this

              itop=min(max(itop,1),lm) !zj quick fix
              ibtc=min(max(ibtc,1),lm) !zj quick fix
 
            !==========================================

              fluxr(i, 8-j) = fluxr(i, 8-j) + tem0d
              fluxr(i,11-j) = fluxr(i,11-j) + prsi(i,itop+kt) * tem0d
              fluxr(i,14-j) = fluxr(i,14-j) + prsi(i,ibtc+kb) * tem0d
              fluxr(i,17-j) = fluxr(i,17-j) + tgrs(i,itop) * tem0d
            enddo
          enddo
        endif

!  ---  save optional vertically integrated aerosol optical depth at
!       wavelenth of 550nm aerodp(:,1), and other optional aod for
!       individual species aerodp(:,2:NSPC1)

!       if ( laswflg ) then
!         if ( NFXR > 33 ) then
!           do i = 1, IM
!             fluxr(i,34) = fluxr(i,34) + dtsw*aerodp(i,1)  ! total aod at 550nm (all species)
!           enddo

!           if ( lspcodp ) then
!             do j = 2, NSPC1
!               k = 33 + j

!               do i = 1, IM
!                 fluxr(i,k) = fluxr(i,k) + dtsw*aerodp(i,j) ! aod at 550nm for indiv species
!               enddo
!             enddo
!           endif     ! end_if_lspcodp
!         else
!           print *,'  !Error! Need to increase array fluxr size NFXR ',&
!    &              ' to be able to output aerosol optical depth'
!           stop
!         endif     ! end_if_nfxr
!       endif       ! end_if_laswflg

      endif                                ! end_if_lssav
!
      return
!...................................
      end subroutine grrad_nmmb
!-----------------------------------


!
!.............................................!
      end module module_radiation_driver_nmmb !
!=============================================!
