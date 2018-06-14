!!!!!  ==========================================================  !!!!!
!!!!!            'module_radiation_surface' description            !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!    this module sets up surface albedo for sw radiation and surface   !
!    emissivity for lw radiation.                                      !
!                                                                      !
!                                                                      !
!    in the module, the externally callabe subroutines are :           !
!                                                                      !
!      'sfc_init'   -- initialization radiation surface data           !
!         inputs:                                                      !
!           ( IALB, IEMS, me )                                         !
!         outputs:                                                     !
!           (none)                                                     !
!                                                                      !
!      'setalb'     -- set up four-component surface albedoes          !
!         inputs:                                                      !
!           (slmsk,snowf,sncovr,snoalb,zorlf,coszf,tsknf,tairf,hprif,  !
!            alvsf,alnsf,alvwf,alnwf,facsf,facwf,fice,tisfc            !
!            IMAX)                                                     !
!         outputs:                                                     !
!           (sfcalb)                                                   !
!                                                                      !
!      'setemis'    -- set up surface emissivity for lw radiation      !
!         inputs:                                                      !
!           (xlon,xlat,slmsk,snowf,sncovr,zorlf,tsknf,tairf,hprif,     !
!            IMAX)                                                     !
!         outputs:                                                     !
!           (sfcemis)                                                  !
!                                                                      !
!    external modules referenced:                                      !
!                                                                      !
!       'module machine'             in 'machine.f'                    !
!       'module physcons'            in 'physcons.f'                   !
!       'module module_iounitdef'    in 'iounitdef.f'                  !
!                                                                      !
!                                                                      !
!    program history log:                                              !
!           1995   y.t. hou     - created albaer.f (include albedo     !
!                                 and aerosols calculations)           !
!      nov  1997   y.t. hou     - modified snow albedo                 !
!      jan  1998   y.t. hou     - included grumbine's sea-ice scheme   !
!      feb  1998   h.l. pan     - seasonal interpolation in cycle      !
!      mar  2000   y.t. hou     - modified to use opac aerosol data    !
!      apr  2003   y.t. hou     - seperate albedo and aerosols into    !
!                    two subroutines, rewritten in f90 modulized form  !
!      jan  2005   s. moorthi   - xingren's sea-ice fraction added     !
!      apr  2005   y.t. hou     - revised module structure             !
!      feb  2006   y.t. hou     - add varying surface emissivity,      !
!                    modified sfc albedo structure for modis shceme    !
!      Mar  2006   s. moorthi   - added surface temp over ice fraction !
!      mar  2007   c. marshall & h. wei                                !
!                               - added modis based sfc albedo scheme  !
!      may  2007   y. hou & s. moorthi                                 !
!                               - fix bug in modis albedo over ocean   !
!      aug  2007   h. wei & s. moorthi                                 !
!                               - fix bug in modis albedo over sea-ice !
!      aug  2007   y. hou       - fix bug in emissivity over ocean in  !
!                                 the modis scheme option              !
!           2009   f. yang      - replaced orig zenith angle treatment !
!                    for the direct beam part of surface albedo with a !
!                    better fit to obs formula over land. yang et al.  !
!                    (2008, jamc).                                     !
!      jul  2010   y. hou       - added option (ialb=2) for sw surface !
!                    albedo that uses input values directly.           !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!


!=============================================!
      module module_radiation_surface_nmmb    !
!.....................................=====...!
!
      USE ESMF

      use physparam,         only : ialbflg, iemsflg, semis_file,       &
     &                              kind_phys

      use physcons,          only : con_t0c, con_ttp, con_pi, con_tice
      use module_iounitdef,  only : NIRADSF
!
      implicit   none
!
      private

!  ---  version tag and last revision date
      character(40), parameter ::                                       &
     &   VTAGSFC='NCEP-Radiation_surface   v5.1  Nov 2012 '
!    &   VTAGSFC='NCEP-Radiation_surface   v5.0  Aug 2012 '

!  ---  constant parameters

      integer, parameter, public :: NF_ALBD = 4  ! num of sfc albedo components

      integer, parameter, public :: IMXEMS = 360 ! num of lon-pts in glb emis-type map
      integer, parameter, public :: JMXEMS = 180 ! num of lat-pts in glb emis-type map

      real (kind=kind_phys), parameter :: f_zero = 0.0
      real (kind=kind_phys), parameter :: f_one  = 1.0

      real (kind=kind_phys), parameter :: rad2dg= 180.0 / con_pi

!  ---  global surface emissivity index array set up in 'sfc_init'

      integer, allocatable  ::  idxems(:,:)
      integer :: iemslw = 0
!
      public  sfcinit, sfc_init, setalb, setemis

! =================
      contains
! =================


!-----------------------------------
      subroutine sfc_init                                               &
!...................................

!  ---  inputs:
     &     ( me )
!  ---  outputs: ( none )

!  ===================================================================  !
!                                                                       !
!  this program is the initialization program for surface radiation     !
!  related quantities (albedo, emissivity, etc.)                        !
!                                                                       !
! program history log:                                                  !
!      feb.  2006   yu-tai hou     - created                            !
!                                                                       !
! usage:         call sfc_init                                          !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                              !
!     IALB          - control flag for surface albedo schemes           !
!                     =0: climatology, based on surface veg types       !
!                     =1: modis retrieval based scheme from boston univ.!
!                     =2: use externally provided albedoes directly.    !
!     IEMS          - control flag for surface emissivity schemes       !
!                     =0: fixed value of 1.0                            !
!                     =1: varying value based on surface veg types      !
!      me           - print control flag                                !
!                                                                       !
!  outputs: (none) to module variables only                             !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: me

!  ---  outputs: ( none )

!  ---  locals:
      integer    :: i, k, ia, ja
      logical    :: file_exist
      character  :: cline*80
!
!===> ...  begin here
!
      if ( me == 0 ) print *, VTAGSFC   ! print out version tag

!  --- ...  initialization of surface albedo section

      if ( ialbflg == 0 ) then

        if ( me == 0 ) then
          print *,' - Using climatology surface albedo scheme for sw'
        endif

      else if ( ialbflg == 1 ) then

        if ( me == 0 ) then
          print *,' - Using MODIS based land surface albedo for sw'
        endif

      else if ( ialbflg == 2 ) then

        if ( me == 0 ) then
          print *,' - Using externally provided surface albedo for sw'
        endif

      else
        print *,' !! ERROR in Albedo Scheme Setting, IALB=',ialbflg
        ! stop
        CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif    ! end if_ialbflg_block

!  --- ...  initialization of surface emissivity section

      iemslw = mod(iemsflg, 10)          ! emissivity control

      if ( iemslw == 0 ) then           ! fixed sfc emis at 1.0

        if ( me == 0 ) then
          print *,' - Using Fixed Surface Emissivity = 1.0 for lw'
        endif

      elseif ( iemslw == 1 ) then       ! input sfc emiss type map

!  ---  allocate data space
        if ( .not. allocated(idxems) ) then
          allocate ( idxems(IMXEMS,JMXEMS)    )
        endif

!  ---  check to see if requested emissivity data file existed

        inquire (file=semis_file, exist=file_exist)

        if ( .not. file_exist ) then
          if ( me == 0 ) then
            print *,' - Using Varying Surface Emissivity for lw'
            print *,'   Requested data file "',semis_file,'" not found!'
            print *,'   Change to fixed surface emissivity = 1.0 !'
          endif

          iemslw = 0
        else
          close(NIRADSF)
          open (NIRADSF,file=semis_file,form='formatted',status='old')
          rewind NIRADSF

          read (NIRADSF,12) cline
  12      format(a80)

          read (NIRADSF,14) idxems
  14      format(80i1)

          if ( me == 0 ) then
            print *,' - Using Varying Surface Emissivity for lw'
            print *,'   Opened data file: ', semis_file
            print *, cline
!check      print *,' CHECK: Sample emissivity index data'
!           ia = IMXEMS / 5
!           ja = JMXEMS / 5
!           print *, idxems(1:IMXEMS:ia,1:JMXEMS:ja)
          endif

          close(NIRADSF)
        endif    ! end if_file_exist_block

      else
        print *,' !! ERROR in Emissivity Scheme Setting, IEMS=',iemsflg
        ! stop
        CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif   ! end if_iemslw_block

!
      return
!...................................
      end subroutine sfc_init
!-----------------------------------


!-----------------------------------
      subroutine sfcinit                                                &
!...................................

!  ---  inputs:
     &     ( LMAX, iflip, IALB, IEMS, me )
!  ---  outputs: ( none )

!  ===================================================================  !
!                                                                       !
!  this program is the initialization program for surface radiation     !
!  related quantities (albedo, emissivity, etc.)                        !
!                                                                       !
! program history log:                                                  !
!      feb.  2006   yu-tai hou     - created                            !
!                                                                       !
! usage:         call sfcinit                                           !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                              !
!     LMAX          - model vertical layer dimension                    !
!     iflip         - control flag for direction of vertical index      !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!     IALB          - control flag for surface albedo schemes           !
!                     =0: climatology, based on surface veg types       !
!                     =1: modis retrieval based scheme from boston univ.!
!                     =2: use externally provided albedoes directly.    !
!     IEMS          - control flag for surface emissivity schemes       !
!                     =0: fixed value of 1.0                            !
!                     =1: varying value based on surface veg types      !
!      me           - print control flag                                !
!                                                                       !
!  outputs: (none) to module variables only                             !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: LMAX, iflip, IALB, IEMS, me

!  ---  outputs: ( none )

!  ---  locals:
      integer    :: i, k, ia, ja
      logical    :: file_exist
      character  :: cline*80, cfems*24

      data  cfems / 'sfc_emissivity_idx.txt  ' /

!
!===> ...  begin here
!

!  --- ...  initialization of surface albedo section

      ialbflg = IALB

      if ( IALB == 0 ) then

        if ( me == 0 ) then
          print *,' - Using climatology surface albedo scheme for sw'
        endif

      else if ( IALB == 1 ) then

        if ( me == 0 ) then
          print *,' - Using MODIS based land surface albedo for sw'
        endif

      else if ( IALB == 2 ) then

        if ( me == 0 ) then
          print *,' - Using externally provided surface albedo for sw'
        endif

      else
        print *,' !! ERROR in Albedo Scheme Setting, IALB=',IALB
        ! stop
        CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif    ! end if_IALB_block

!  --- ...  initialization of surface emissivity section

      iemsflg = IEMS

      if ( IEMS == 0 ) then           ! fixed sfc emis at 1.0

        if ( me == 0 ) then
          print *,' - Using Fixed Surface Emissivity = 1.0 for lw'
        endif

      elseif ( IEMS == 1 ) then       ! input sfc emiss type map

!  ---  allocate data space
        if ( .not. allocated(idxems) ) then
          allocate ( idxems(IMXEMS,JMXEMS)    )
        endif

!  ---  check to see if requested emissivity data file existed

        inquire (file=cfems, exist=file_exist)

        if ( .not. file_exist ) then
          if ( me == 0 ) then
            print *,' - Using Varying Surface Emissivity for lw'
            print *,'   Requested data file "',cfems,'" not found!'
            print *,'   Change to fixed surface emissivity = 1.0 !'
          endif

          iemsflg = 0
        else
          open (NIRADSF,file=cfems,form='formatted',status='old')
          rewind NIRADSF

          read (NIRADSF,12) cline
  12      format(a80)

          read (NIRADSF,14) idxems
  14      format(80i1)

          if ( me == 0 ) then
            print *,' - Using Varying Surface Emissivity for lw'
            print *,'   Opened data file: ',cfems
            print *, cline
!check      print *,' CHECK: Sample emissivity index data'
!           ia = IMXEMS / 5
!           ja = JMXEMS / 5
!           print *, idxems(1:IMXEMS:ia,1:JMXEMS:ja)
          endif
        endif    ! end if_file_exist_block

      else
        print *,' !! ERROR in Emissivity Scheme Setting, IEMS=',IEMS
        ! stop
        CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif   ! end if_IEMS_block

!
      return
!...................................
      end subroutine sfcinit
!-----------------------------------



!-----------------------------------
      subroutine setalb                                                 &
!...................................

!  ---  inputs:
     &     ( slmsk,snowf,sncovr,snoalb,zorlf,coszf,tsknf,tairf,hprif,   &
     &       alvsf,alnsf,alvwf,alnwf,facsf,facwf,fice,tisfc,            &
     &       IMAX,                                                      &
!  ---  outputs:
     &       sfcalb                                                     &
     &     )

!  ===================================================================  !
!                                                                       !
!  this program computes four components of surface albedos (i.e.       !
!  vis-nir, direct-diffused) according to controflag ialbflg.           !
!   1) climatological surface albedo scheme (briegleb 1992)             !
!   2) modis retrieval based scheme from boston univ.                   !
!                                                                       !
!                                                                       !
! usage:         call setalb                                            !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                              !
!     slmsk (IMAX)  - sea(0),land(1),ice(2) mask on fcst model grid     !
!     snowf (IMAX)  - snow depth water equivalent in mm                 !
!     sncovr(IMAX)  - ialgflg=0: not used                               !
!                     ialgflg=1: snow cover over land in fraction       !
!     snoalb(IMAX)  - ialbflg=0: not used                               !
!                     ialgflg=1: max snow albedo over land in fraction  !
!     zorlf (IMAX)  - surface roughness in cm                           !
!     coszf (IMAX)  - cosin of solar zenith angle                       !
!     tsknf (IMAX)  - ground surface temperature in k                   !
!     tairf (IMAX)  - lowest model layer air temperature in k           !
!     hprif (IMAX)  - topographic sdv in m                              !
!           ---  for ialbflg=0 climtological albedo scheme  ---         !
!     alvsf (IMAX)  - 60 degree vis albedo with strong cosz dependency  !
!     alnsf (IMAX)  - 60 degree nir albedo with strong cosz dependency  !
!     alvwf (IMAX)  - 60 degree vis albedo with weak cosz dependency    !
!     alnwf (IMAX)  - 60 degree nir albedo with weak cosz dependency    !
!           ---  for ialbflg=1 modis based land albedo scheme ---       !
!     alvsf (IMAX)  - visible black sky albedo at zenith 60 degree      !
!     alnsf (IMAX)  - near-ir black sky albedo at zenith 60 degree      !
!     alvwf (IMAX)  - visible white sky albedo                          !
!     alnwf (IMAX)  - near-ir white sky albedo                          !
!           ---  for ialbflg=2 externally provided surface albedo  ---  !
!     alvsf (IMAX)  - visible direct beam albedo                        !
!     alnsf (IMAX)  - near-ir direct beam albedo                        !
!     alvwf (IMAX)  - visible diffused albedo                           !
!     alnwf (IMAX)  - near-ir diffused albedo                           !
!                                                                       !
!     facsf (IMAX)  - fractional coverage with strong cosz dependency   !
!     facwf (IMAX)  - fractional coverage with weak cosz dependency     !
!     fice  (IMAX)  - sea-ice fraction                                  !
!     tisfc (IMAX)  - sea-ice surface temperature                       !
!     IMAX          - array horizontal dimension                        !
!                                                                       !
!  outputs:                                                             !
!     sfcalb(IMAX,NF_ALBD)                                              !
!           ( :, 1) -     near ir direct beam albedo                    !
!           ( :, 2) -     near ir diffused albedo                       !
!           ( :, 3) -     uv+vis direct beam albedo                     !
!           ( :, 4) -     uv+vis diffused albedo                        !
!                                                                       !
!  module internal control variables:                                   !
!     ialbflg       - =0 use the default climatology surface albedo     !
!                     =1 use modis retrieved albedo and input snow cover!
!                        for land areas                                 !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer, intent(in) :: IMAX

      real (kind=kind_phys), dimension(:), intent(in) ::                &
     &       slmsk, snowf, zorlf, coszf, tsknf, tairf, hprif,           &
     &       alvsf, alnsf, alvwf, alnwf, facsf, facwf, fice, tisfc,     &
     &       sncovr, snoalb

!  ---  outputs
      real (kind=kind_phys), dimension(:,:), intent(out) :: sfcalb

!  ---  locals:
      real (kind=kind_phys) :: asnvb, asnnb, asnvd, asnnd, asevb        &
     &,     asenb, asevd, asend, fsno,  fsea,  rfcs,  rfcw,  flnd       &
     &,     asnow, argh,  hrgh,  fsno0, fsno1, flnd0, fsea0, csnow      &
     &,     a1, a2, b1, b2, b3

      real (kind=kind_phys) ffw, dtgd

      integer :: i, k

!
!===> ...  begin here
!

      if ( ialbflg == 0 ) then   ! use climatological albedo scheme

        do i = 1, IMAX

! --- modified snow albedo scheme - units convert to m
!     (originally snowf in mm; zorlf in cm)

         asnow = 0.02*snowf(i)
         argh  = min(0.50, max(.025, 0.01*zorlf(i)))
         hrgh  = min(1.0, max(0.20, 1.0577-1.1538e-3*hprif(i) ) )
         fsno0 = asnow / (argh + asnow) * hrgh
         if (nint(slmsk(i)) == 0 .and. tsknf(i) > con_tice) fsno0 = 0.0
         fsno1 = 1.0 - fsno0
         flnd0 = min(1.0, facsf(i)+facwf(i))
         fsea0 = max(0.0, 1.0 - flnd0)
         fsno  = fsno0
         fsea  = fsea0 * fsno1
         flnd  = flnd0 * fsno1

! --- diffused sea surface albedo

         if (tsknf(i) >= 271.5) then
            asevd = 0.06
            asend = 0.06
         elseif (tsknf(i) < 271.1) then
            asevd = 0.70
            asend = 0.65
         else
            a1 = (tsknf(i) - 271.1)**2
            asevd = 0.7 - 4.0*a1
            asend = 0.65 - 3.6875*a1
         endif

! --- diffused snow albedo

         if (nint(slmsk(i)) == 2) then
            ffw   = 1.0 - fice(i)
            if (ffw < 1.0) then
               dtgd = max(0.0, min(5.0, (con_ttp-tisfc(i)) ))
               b1   = 0.03 * dtgd
            else
               b1 = 0.0
            endif

            b3   = 0.06 * ffw
            asnvd = (0.70 + b1) * fice(i) + b3
            asnnd = (0.60 + b1) * fice(i) + b3
            asevd = 0.70        * fice(i) + b3
            asend = 0.60        * fice(i) + b3
         else
            asnvd = 0.90
            asnnd = 0.75
         endif

! --- direct snow albedo

         if (coszf(i) < 0.5) then
            csnow = 0.5 * (3.0 / (1.0+4.0*coszf(i)) - 1.0)
            asnvb = min( 0.98, asnvd+(1.0-asnvd)*csnow )
            asnnb = min( 0.98, asnnd+(1.0-asnnd)*csnow )
         else
            asnvb = asnvd
            asnnb = asnnd
         endif

! --- direct sea surface albedo

         if (coszf(i) > 0.0001) then
!           rfcs = 1.4 / (1.0 + 0.8*coszf(i))     ! briegleb's zenith angle treatment
!           rfcw = 1.3 / (1.0 + 0.6*coszf(i))
            rfcs = 2.14 / (1.0 + 1.48*coszf(i))   ! f. yang's zenith angle treatment
            rfcw = rfcs

            if (tsknf(i) >= con_t0c) then
              asevb = max(asevd, 0.026/(coszf(i)**1.7+0.065)            &
     &              + 0.15 * (coszf(i)-0.1) * (coszf(i)-0.5)            &
     &              * (coszf(i)-1.0))
              asenb = asevb
            else
              asevb = asevd
              asenb = asend
            endif
         else
            rfcs  = 1.0
            rfcw  = 1.0
            asevb = asevd
            asenb = asend
         endif

         a1   = alvsf(i) * facsf(i)
         b1   = alvwf(i) * facwf(i)
         a2   = alnsf(i) * facsf(i)
         b2   = alnwf(i) * facwf(i)
         sfcalb(i,1) = (a2*rfcs+b2*rfcw)*flnd + asenb*fsea + asnnb*fsno
         sfcalb(i,2) = (a2 + b2) * 0.96 *flnd + asend*fsea + asnnd*fsno
         sfcalb(i,3) = (a1*rfcs+b1*rfcw)*flnd + asevb*fsea + asnvb*fsno
         sfcalb(i,4) = (a1 + b1) * 0.96 *flnd + asevd*fsea + asnvd*fsno

        enddo    ! end_do_i_loop

      elseif ( ialbflg == 1 ) then ! use modis based albedo for land area

        do i = 1, IMAX

! --- snow cover input directly from land model, no conversion needed

         fsno0 = sncovr(i)

         if (nint(slmsk(i)) == 0 .and. tsknf(i) > con_tice) fsno0 = 0.0

         if (nint(slmsk(i)) == 2) then
           asnow = 0.02*snowf(i)
           argh  = min(0.50, max(.025, 0.01*zorlf(i)))
           hrgh  = min(1.0, max(0.20, 1.0577-1.1538e-3*hprif(i) ) )
           fsno0 = asnow / (argh + asnow) * hrgh
         endif

         fsno1 = 1.0 - fsno0
         flnd0 = min(1.0, facsf(i)+facwf(i))
         fsea0 = max(0.0, 1.0 - flnd0)
         fsno  = fsno0
         fsea  = fsea0 * fsno1
         flnd  = flnd0 * fsno1

! --- diffused sea surface albedo

         if (tsknf(i) >= 271.5) then
            asevd = 0.06
            asend = 0.06
         elseif (tsknf(i) < 271.1) then
            asevd = 0.70
            asend = 0.65
         else
            a1 = (tsknf(i) - 271.1)**2
            asevd = 0.7 - 4.0*a1
            asend = 0.65 - 3.6875*a1
         endif

! --- diffused snow albedo, land area use input max snow albedo

         if (nint(slmsk(i)) == 2) then
            ffw   = 1.0 - fice(i)
            if (ffw < 1.0) then
               dtgd = max(0.0, min(5.0, (con_ttp-tisfc(i)) ))
               b1   = 0.03 * dtgd
            else
               b1 = 0.0
            endif

            b3   = 0.06 * ffw
            asnvd = (0.70 + b1) * fice(i) + b3
            asnnd = (0.60 + b1) * fice(i) + b3
            asevd = 0.70        * fice(i) + b3
            asend = 0.60        * fice(i) + b3
         else
            asnvd = snoalb(i)
            asnnd = snoalb(i)
         endif

! --- direct snow albedo

         if (nint(slmsk(i)) == 2) then
           if (coszf(i) < 0.5) then
              csnow = 0.5 * (3.0 / (1.0+4.0*coszf(i)) - 1.0)
              asnvb = min( 0.98, asnvd+(1.0-asnvd)*csnow )
              asnnb = min( 0.98, asnnd+(1.0-asnnd)*csnow )
           else
             asnvb = asnvd
             asnnb = asnnd
           endif
         else
           asnvb = snoalb(i)
           asnnb = snoalb(i)
         endif

! --- direct sea surface albedo

         if (coszf(i) > 0.0001) then

!           rfcs = 1.89 - 3.34*coszf(i) + 4.13*coszf(i)*coszf(i)        &
!    &           - 2.02*coszf(i)*coszf(i)*coszf(i)
            rfcs = 1.775/(1.0+1.55*coszf(i))           ! f. yang's zenith angle treatment

            if (tsknf(i) >= con_t0c) then
              asevb = max(asevd, 0.026/(coszf(i)**1.7+0.065)            &
     &              + 0.15 * (coszf(i)-0.1) * (coszf(i)-0.5)            &
     &              * (coszf(i)-1.0))
              asenb = asevb
            else
              asevb = asevd
              asenb = asend
            endif
         else
            rfcs  = 1.0
            asevb = asevd
            asenb = asend
         endif

         sfcalb(i,1) = alnsf(i)*rfcs*flnd + asenb*fsea + asnnb*fsno
         sfcalb(i,2) = alnwf(i)     *flnd + asend*fsea + asnnd*fsno
         sfcalb(i,3) = alvsf(i)*rfcs*flnd + asevb*fsea + asnvb*fsno
         sfcalb(i,4) = alvwf(i)     *flnd + asevd*fsea + asnvd*fsno

        enddo    ! end_do_i_loop

      else                         ! use externally provided albedoes

        do i = 1, IMAX
          sfcalb(i,1) = alnsf(i)
          sfcalb(i,2) = alnwf(i)
          sfcalb(i,3) = alvsf(i)
          sfcalb(i,4) = alvwf(i)
        enddo

      endif   ! end if_ialbflg
!
      return
!...................................
      end subroutine setalb
!-----------------------------------


!-----------------------------------
      subroutine setemis                                                &
!...................................

!  ---  inputs:
     &     ( xlon,xlat,slmsk,snowf,sncovr,zorlf,tsknf,tairf,hprif,      &
     &       IMAX,                                                      &
!  ---  outputs:
     &       sfcemis                                                    &
     &     )

!  ===================================================================  !
!                                                                       !
!  this program computes surface emissivity for lw radiation.           !
!                                                                       !
!  usage:         call setemis                                          !
!                                                                       !
!  subprograms called:  none                                            !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                              !
!     xlon  (IMAX)  - longitude in radiance                             !
!     xlat  (IMAX)  - latitude  in radiance                             !
!     slmsk (IMAX)  - sea(0),land(1),ice(2) mask on fcst model grid     !
!     snowf (IMAX)  - snow depth water equivalent in mm                 !
!     sncovr(IMAX)  - ialbflg=1: snow cover over land in fraction       !
!     zorlf (IMAX)  - surface roughness in cm                           !
!     tsknf (IMAX)  - ground surface temperature in k                   !
!     tairf (IMAX)  - lowest model layer air temperature in k           !
!     hprif (IMAX)  - topographic sdv in m                              !
!     IMAX          - array horizontal dimension                        !
!                                                                       !
!  outputs:                                                             !
!     sfcemis(IMAX) - surface emissivity                                !
!                                                                       !
!  -------------------------------------------------------------------  !
!                                                                       !
!  surface type definations:                                            !
!     1. open water                   2. grass/wood/shrub land          !
!     3. tundra/bare soil             4. sandy desert                   !
!     5. rocky desert                 6. forest                         !
!     7. ice                          8. snow                           !
!                                                                       !
!  input index data lon from 0 towards east, lat from n to s            !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer, intent(in) :: IMAX

      real (kind=kind_phys), dimension(:), intent(in) ::                &
     &       xlon,xlat, slmsk, snowf,sncovr, zorlf, tsknf, tairf, hprif

!  ---  outputs
      real (kind=kind_phys), dimension(:), intent(out) :: sfcemis

!  ---  locals:
      integer :: i, i1, i2, j1, j2, idx

      real (kind=kind_phys) :: dltg, hdlt, tmp1, tmp2,                  &
     &      asnow, argh, hrgh, fsno, fsno0, fsno1

!  ---  reference emiss value for diff surface emiss index
!       1-open water, 2-grass/shrub land, 3-bare soil, tundra,
!       4-sandy desert, 5-rocky desert, 6-forest, 7-ice, 8-snow

      real (kind=kind_phys) ::  emsref(8)
      data  emsref / 0.97, 0.95, 0.94, 0.90, 0.93, 0.96, 0.96, 0.99 /

!
!===> ...  begin here
!

      if ( iemslw == 0 ) then        ! sfc emiss default to 1.0

        sfcemis(:) = 1.0
        return

      else                           ! emiss set by sfc type and condition

        dltg = 360.0 / float(IMXEMS)
        hdlt = 0.5 * dltg

!  --- ...  mapping input data onto model grid
!           note: this is a simple mapping method, an upgrade is needed if
!           the model grid is much corcer than the 1-deg data resolution

        lab_do_IMAX : do i = 1, IMAX

          if ( nint(slmsk(i)) == 0 ) then          ! sea point

            sfcemis(i) = emsref(1)

          else if ( nint(slmsk(i)) == 2 ) then     ! sea-ice

            sfcemis(i) = emsref(7)

          else                                     ! land

!  ---  map grid in longitude direction
            i2 = 1
            j2 = 1

            tmp1 = xlon(i) * rad2dg
            if (tmp1 < f_zero) tmp1 = tmp1 + 360.0

            lab_do_IMXEMS : do i1 = 1, IMXEMS
              tmp2 = dltg * (i1 - 1) + hdlt

              if (abs(tmp1-tmp2) <= hdlt) then
               i2 = i1
                exit lab_do_IMXEMS
              endif
            enddo  lab_do_IMXEMS

!  ---  map grid in latitude direction
            tmp1 = xlat(i) * rad2dg

            lab_do_JMXEMS : do j1 = 1, JMXEMS
              tmp2 = 90.0 - dltg * (j1 - 1)

              if (abs(tmp1-tmp2) <= hdlt) then
                j2 = j1
                exit lab_do_JMXEMS
              endif
            enddo  lab_do_JMXEMS


            idx = max( 2, idxems(i2,j2) )
            if ( idx >= 7 ) idx = 2
            sfcemis(i) = emsref(idx)

          endif   ! end if_slmsk_block

!  ---  check for snow covered area

          if ( ialbflg==1 .and. nint(slmsk(i))==1 ) then ! input land area snow cover

            fsno0 = sncovr(i)
            fsno1 = f_one - fsno0
            sfcemis(i) = sfcemis(i)*fsno1 + emsref(8)*fsno0

          else                                           ! compute snow cover from snow depth
            if ( snowf(i) > f_zero ) then
              asnow = 0.02*snowf(i)
              argh  = min(0.50, max(.025, 0.01*zorlf(i)))
              hrgh  = min(1.0, max(0.20, 1.0577-1.1538e-3*hprif(i) ) )
              fsno0 = asnow / (argh + asnow) * hrgh
              if (nint(slmsk(i)) == 0 .and. tsknf(i) > 271.2)           &
     &                               fsno0=f_zero
              fsno1 = f_one - fsno0
              sfcemis(i) = sfcemis(i)*fsno1 + emsref(8)*fsno0
            endif

          endif                                          ! end if_ialbflg

        enddo  lab_do_IMAX

      endif   ! end if_iemslw_block

!chk  print *,' In setemis, iemsflg, sfcemis =',iemsflg,sfcemis

!
      return
!...................................
      end subroutine setemis
!-----------------------------------

!
!..............................................!
      end module module_radiation_surface_nmmb !
!==============================================!
