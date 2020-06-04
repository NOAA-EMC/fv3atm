!-----------------------------------
      subroutine rad_initialize                                         &
!...................................
!  ---  inputs:
     &     ( si,levr,ictm,isol,ico2,iaer,ialb,iems,ntcw, num_p2d,       &
     &       num_p3d,npdf3d,ntoz,iovr_sw,iovr_lw,isubc_sw,isubc_lw,     &
     &       icliq_sw,crick_proof,ccnorm,                               &                  
     &       imp_physics,norad_precip,idate,iflip,me )
!  ---  outputs: ( none )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:   rad_initialize - a subprogram to initialize radiation   !
!                                                                       !
! usage:        call rad_initialize                                     !
!                                                                       !
! attributes:                                                           !
!   language:  fortran 90                                               !
!                                                                       !
! program history:                                                      !
!   mar 2012  - yu-tai hou   create the program to initialize fixed     !
!                 control variables for radiaion processes.  this       !
!                 subroutine is called at the start of model run.       !
!   nov 2012  - yu-tai hou   modified control parameter through         !
!                 module 'physparam'.                                   !
!   mar 2014  - sarah lu  iaermdl is determined from iaer               !        
!   jul 2014  - s moorthi add npdf3d for pdf clouds                     !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input parameters:                                                     !
!   si               : model vertical sigma interface or equivalence    !
!   levr             : number of model vertical layers                  !
!   ictm             :=yyyy#, external data time/date control flag      !
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
!   isol             := 0: use the old fixed solar constant in "physcon"!
!                     =10: use the new fixed solar constant in "physcon"!
!                     = 1: use noaa ann-mean tsi tbl abs-scale data tabl!
!                     = 2: use noaa ann-mean tsi tbl tim-scale data tabl!
!                     = 3: use cmip5 ann-mean tsi tbl tim-scale data tbl!
!                     = 4: use cmip5 mon-mean tsi tbl tim-scale data tbl!
!   ico2             :=0: use prescribed global mean co2 (old  oper)    !
!                     =1: use observed co2 annual mean value only       !
!                     =2: use obs co2 monthly data with 2-d variation   !
!   iaer             : 4-digit aerosol flag (dabc for aermdl,volc,lw,sw)!
!                     d: =0 or none, opac-climatology aerosol scheme    !                
!                        =1 use gocart climatology aerosol scheme       !  
!                        =2 use gocart progostic aerosol scheme         !  
!                     a: =0 use background stratospheric aerosol        !
!                        =1 incl stratospheric vocanic aeros            !
!                     b: =0 no topospheric aerosol in lw radiation      !
!                        =1 include tropspheric aerosols for lw         !
!                     c: =0 no topospheric aerosol in sw radiation      !
!                        =1 include tropspheric aerosols for sw         !
!   ialb             : control flag for surface albedo schemes          !
!                     =0: climatology, based on surface veg types       !
!                     =1: modis retrieval based surface albedo scheme   !
!   iems             : ab 2-digit control flag                          !
!                     a: =0 set sfc air/ground t same for lw radiation  !
!                        =1 set sfc air/ground t diff for lw radiation  !
!                     b: =0 use fixed sfc emissivity=1.0 (black-body)   !
!                        =1 use varying climtology sfc emiss (veg based)!
!                        =2 future development (not yet)                !
!   ntcw             :=0 no cloud condensate calculated                 !
!                     >0 array index location for cloud condensate      !
!   num_p3d          :=3: ferrier's microphysics cloud scheme           !
!                     =4: zhao/carr/sundqvist microphysics cloud        !
!   npdf3d            =0 no pdf clouds                                  !
!                     =3 (when num_p3d=4) pdf clouds with zhao/carr/    !
!                        sundqvist scheme                               !
!   ntoz             : ozone data control flag                          !
!                     =0: use climatological ozone profile              !
!                     >0: use interactive ozone profile                 !
!   icliq_sw         : sw optical property for liquid clouds            !
!                     =0:input cld opt depth, ignoring iswcice setting  !
!                     =1:cloud optical property scheme based on Hu and  !
!                        Stamnes(1993) \cite hu_and_stamnes_1993 method !
!                     =2:cloud optical property scheme based on Hu and  !
!                        Stamnes(1993) -updated                         !
!   iovr_sw/iovr_lw  : control flag for cloud overlap (sw/lw rad)       !
!                     =0: random overlapping clouds                     !
!                     =1: max/ran overlapping clouds                    !
!                     =2: maximum overlap clouds       (mcica only)     !
!                     =3: decorrelation-length overlap (mcica only)     !
!   isubc_sw/isubc_lw: sub-column cloud approx control flag (sw/lw rad) !
!                     =0: with out sub-column cloud approximation       !
!                     =1: mcica sub-col approx. prescribed random seed  !
!                     =2: mcica sub-col approx. provided random seed    !
!   crick_proof      : control flag for eliminating CRICK               !
!   ccnorm           : control flag for in-cloud condensate mixing ratio!
!   norad_precip     : control flag for not using precip in radiation   !
!   idate(4)         : ncep absolute date and time of initial condition !
!                      (hour, month, day, year)                         !
!   iflip            : control flag for direction of vertical index     !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   me               : print control flag                               !
!                                                                       !
!  subroutines called: radinit                                          !
!                                                                       !
!  ===================================================================  !
!
      use physparam, only : isolar , ictmflg, ico2flg, ioznflg, iaerflg,&
!    &             iaermdl, laswflg, lalwflg, lavoflg, icldflg,         &
     &             iaermdl,                            icldflg,         &
     &             iovrsw , iovrlw , lcrick , lcnorm , lnoprec,         &
     &             ialbflg, iemsflg, isubcsw, isubclw, ivflip , ipsd0,  &
     &             iswcliq,                                             &
     &             kind_phys

      use module_radiation_driver, only : radinit
!
      implicit   none

!  ---  input:
      integer,  intent(in) :: levr, ictm, isol, ico2, iaer, num_p2d,    &
     &       ntcw, ialb, iems, num_p3d, npdf3d, ntoz, iovr_sw, iovr_lw, &
     &       isubc_sw, isubc_lw, icliq_sw, iflip, me, idate(4)

      real (kind=kind_phys), intent(in) :: si(levr+1)
      integer, intent(in) :: imp_physics

      logical, intent(in) :: crick_proof, ccnorm, norad_precip

!  ---  output: ( none )

!  ---  local:
!
!===> ...  start here
!
!  ---  set up parameters for radiation initialization

      isolar = isol                     ! solar constant control flag

      ictmflg= ictm                     ! data ic time/date control flag
      ico2flg= ico2                     ! co2 data source control flag
      ioznflg= ntoz                     ! ozone data source control flag

      if ( ictm==0 .or. ictm==-2 ) then
        iaerflg = mod(iaer, 100)        ! no volcanic aerosols for clim hindcast
      else
        iaerflg = mod(iaer, 1000)   
      endif
      iaermdl = iaer/1000               ! control flag for aerosol scheme selection                              
      if ( iaermdl < 0 .or.  (iaermdl>2 .and. iaermdl/=5) ) then
         print *, ' Error -- IAER flag is incorrect, Abort'
         stop 7777
      endif

!     if ( ntcw > 0 ) then
        icldflg = 1                     ! prognostic cloud optical prop scheme
!     else
!       icldflg = 0                     ! no support for diag cloud opt prop scheme
!     endif

      iswcliq = icliq_sw                ! optical property for liquid clouds for sw

      iovrsw = iovr_sw                  ! cloud overlapping control flag for sw
      iovrlw = iovr_lw                  ! cloud overlapping control flag for lw

      lcrick  = crick_proof             ! control flag for eliminating CRICK 
      lcnorm  = ccnorm                  ! control flag for in-cld condensate 
      lnoprec = norad_precip            ! precip effect on radiation flag (ferrier microphysics)
      isubcsw = isubc_sw                ! sub-column cloud approx flag in sw radiation
      isubclw = isubc_lw                ! sub-column cloud approx flag in lw radiation

      ialbflg= ialb                     ! surface albedo control flag
      iemsflg= iems                     ! surface emissivity control flag

      ivflip = iflip                    ! vertical index direction control flag

!  ---  assign initial permutation seed for mcica cloud-radiation
      if ( isubc_sw>0 .or. isubc_lw>0 ) then
!       ipsd0 = 17*idate(1)+43*idate(2)+37*idate(3)+23*idate(4) + ipsd0
        ipsd0 = 17*idate(1)+43*idate(2)+37*idate(3)+23*idate(4)
      endif

      if ( me == 0 ) then
        print *,'  In rad_initialize, before calling radinit'
        print *,' si =',si
        print *,' levr=',levr,' ictm=',ictm,' isol=',isol,' ico2=',ico2,&
     &          ' iaer=',iaer,' ialb=',ialb,' iems=',iems,' ntcw=',ntcw
        print *,' np3d=',num_p3d,' ntoz=',ntoz,' iovr_sw=',iovr_sw,     &
     &          ' iovr_lw=',iovr_lw,' isubc_sw=',isubc_sw,              &
     &          ' isubc_lw=',isubc_lw,' icliq_sw=',icliq_sw,            &
     &          ' iflip=',iflip,'  me=',me
        print *,' crick_proof=',crick_proof,                            &
     &          ' ccnorm=',ccnorm,' norad_precip=',norad_precip
      endif

      call radinit                                                      &
!  ---  inputs:
     &     ( si, levr, imp_physics,  me )
!  ---  outputs:
!          ( none )

      if ( me == 0 ) then
        print *,'  Radiation sub-cloud initial seed =',ipsd0,           &
     &          ' IC-idate =',idate
        print *,' return from rad_initialize - after calling radinit'
      endif
!
      return
!...................................
      end subroutine rad_initialize
!-----------------------------------
