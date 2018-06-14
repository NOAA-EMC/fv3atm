!-----------------------------------
      subroutine rad_initialize_nmmb                                    &
!...................................
!  ---  inputs:
     &     ( si,levr,ictm,isol,ico2,iaer,iaer_mdl,ialb,iems,ntcw,       &
     &       num_p3d,ntoz,iovr_sw,iovr_lw,isubc_sw,isubc_lw,            &
     &       icliq_sw, icice_sw, icliq_lw, icice_lw,                    &
     &       sashal,crick_proof,ccnorm,norad_precip, iflip, me )
!  ---  outputs:
!          ( none )

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
!   iaer             : 3-digit aerosol flag (abc for volc, lw, sw)      !
!                     a: =0 use background stratospheric aerosol        !
!                        =1 incl stratospheric vocanic aeros            !
!                     b: =0 no topospheric aerosol in lw radiation      !
!                        =1 comp tropspheric aero in 1 broad band for lw!
!                        =2 comp tropspheric aero in multi bands for lw !
!                     c: =0 no topospheric aerosol in sw radiation      !
!                        =1 include tropspheric aerosols for sw         !
!   iaer_mdl         : control flag for topospheric aerosol models      !
!                     =0: opac-climatology aerosol model                !
!                     =1: gocart-climatology aerosol model              !
!                     =2: gocart-prognostic aerosol model               !
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
!                     =8: Thompson microphysics cloud                   !
!   ntoz             : ozone data control flag                          !
!                     =0: use climatological ozone profile              !
!                     >0: use interactive ozone profile                 !
!   iovr_sw/iovr_lw  : control flag for cloud overlap (sw/lw rad)       !
!                     =0: random overlapping clouds                     !
!                     =1: max/ran overlapping clouds                    !
!   isubc_sw/isubc_lw: sub-column cloud approx control flag (sw/lw rad) !
!                     =0: with out sub-column cloud approximation       !
!                     =1: mcica sub-col approx. prescribed random seed  !
!                     =2: mcica sub-col approx. provided random seed    !
!   sashal           : control flag for shallow convection              !
!   crick_proof      : control flag for eliminating CRICK               !
!   ccnorm           : control flag for in-cloud condensate mixing ratio!
!   norad_precip     : control flag for not using precip in radiation   !
!   iflip            : control flag for direction of vertical index     !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   me               : print control flag                               !
!                                                                       !
!  outputs: (none)                                                      !
!                                                                       !
!  subroutines called: radinit                                          !
!                                                                       !
!  ===================================================================  !
!
      use physparam,only : isolar , ictmflg, ico2flg, ioznflg, iaerflg, &
     &            iaermdl, laswflg, lalwflg, lavoflg, icldflg, icmphys, &
     &            iovrsw , iovrlw , lsashal, lcrick , lcnorm , lnoprec, &
     &            ialbflg, iemsflg, isubcsw, isubclw, ivflip , ipsd0,   &
     &            ilwcice, ilwcliq, iswcice, iswcliq, kind_phys

      use module_radiation_driver_nmmb, only : radinit_nmmb
!
      implicit   none

!  ---  input:

      integer,  intent(in) :: levr, ictm, isol, ico2, iaer, iaer_mdl,   &
     &       ntcw, ialb, iems, num_p3d, ntoz, iovr_sw, iovr_lw,         &
     &       icliq_sw, icice_sw, icliq_lw, icice_lw,                    &
     &       isubc_sw, isubc_lw, iflip, me

      real (kind=kind_phys), intent(in) :: si(levr+1)
      ! real , intent(in) :: si(levr+1)

      logical, intent(in) :: sashal, crick_proof, ccnorm, norad_precip

!  ---  output: ( none )

!  ---  local:
      integer :: icld
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
        iaerflg = iaer
      endif
      laswflg= (mod(iaerflg,10) > 0)    ! control flag for sw tropospheric aerosol
      lalwflg= (mod(iaerflg/10,10) > 0) ! control flag for lw tropospheric aerosol
      lavoflg= (iaerflg >= 100)         ! control flag for stratospheric volcanic aeros
      iaermdl= iaer_mdl                 ! control flag for aerosol scheme selections

      if ( ntcw > 0 ) then
        icldflg = 1         ! prognostic cloud scheme
      else
        icldflg = 0         ! diagnostic cloud scheme
      endif

  !    if ( num_p3d == 4 ) then
  !      icmphys = 1                     ! zhao/moorthi's prognostic cloud scheme
  !    elseif ( num_p3d == 3 ) then
  !      icmphys = 2                     ! ferrier's microphysics
  !    endif

      icmphys = num_p3d

      iovrsw = iovr_sw                  ! cloud overlapping control flag for sw
      iovrlw = iovr_lw                  ! cloud overlapping control flag for lw

      lsashal = sashal                  ! shallow convection scheme flag
      lcrick  = crick_proof             ! control flag for eliminating CRICK
      lcnorm  = ccnorm                  ! control flag for in-cld condensate
      lnoprec = norad_precip            ! precip effect on radiation flag (ferrier microphysics)
      isubcsw = isubc_sw                ! sub-column cloud approx flag in sw radiation
      isubclw = isubc_lw                ! sub-column cloud approx flag in lw radiation

      iswcliq = icliq_sw                ! sw optical property for liquid clouds
      iswcice = icice_sw                ! sw optical property for ice clouds (only iswcliq>0)
      ilwcliq = icliq_lw                ! lw optical property for liquid clouds
      ilwcice = icice_lw                ! lw optical property for ice clouds (only ilwcliq>0)

      ialbflg = ialb                    ! surface albedo control flag
      iemsflg = iems                    ! surface emissivity control flag

      ivflip = iflip                    ! vertical index direction control flag

! ! ---  assign initial permutation seed for mcica cloud-radiation
  !    if ( isubc_sw>0 .or. isubc_lw>0 ) then
! !      ipsd0 = 17*idate(1)+43*idate(2)+37*idate(3)+23*idate(4) + ipsd0
  !      ipsd0 = 17*idate(1)+43*idate(2)+37*idate(3)+23*idate(4)
  !    endif

      if ( me == 0 ) then
        print *,'  In rad_initialize, before calling radinit'
        print *,' si =',si
        print *,' levr=',levr,' ictm=',ictm,' isol=',isol,' ico2=',ico2,&
     &          ' iaer=',iaer,' ialb=',ialb,' iems=',iems,' ntcw=',ntcw
        print *,' np3d=',num_p3d,' ntoz=',ntoz,' iovrsw=',iovr_sw,      &
     &          ' iovrlw=',iovr_lw,' isubcsw=',isubc_sw,                &
     &          ' isubclw=',isubc_lw,' iflip=',iflip,'  me=',me
        print *,' sashal=',sashal,' crick=',crick_proof,                &
     &          ' ccnorm=',ccnorm,' norad=',norad_precip

        print *,'=== lw optical property for clouds ==='
        print *,' ilwcice=',ilwcice, ' ilwcliq=',ilwcliq

        print *,'=== sw optical property for clouds ==='
        print *,' iswcice=',iswcice, ' iswcliq=',iswcliq

      endif

      call radinit_nmmb                                                 &
!  ---  inputs:
     &     ( si, levr, me )
!  ---  outputs:
!          ( none )

  !    if ( me == 0 ) then
  !      print *,'  Radiation sub-cloud initial seed =',ipsd0,           &
  !   &          ' IC-idate =',idate
  !      print *,' return from rad_initialize - after calling radinit'
  !    endif
!
      return
!...................................
      end subroutine rad_initialize_nmmb
!-----------------------------------
