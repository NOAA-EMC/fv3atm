!-- Modifications:
!   1) recwat_def=5 microns
!   2) Assumed TNW=200 cm^-3 is fixed
!   3) Corrected comments concerning cloud fractions
!   4) Effectively treat all ice as cloud ice, eliminate snow in rsipath_tmp
!   5) Parameterize effective radius of (cloud) ice following eq. (5) from
!      McFarquhar & Heymsfield (1996), discussion concerning area cross
!      section of ice particles (Ac) on p. 2437, and eqs. (3.10) and (3.12)
!      in Fu (1996).
!   
!!!!!              module_radiation_clouds description             !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!    the 'radiation_clouds.f' contains:                                !
!                                                                      !
!       'module_radiation_clouds' ---  compute cloud related quantities!
!                for radiation computations                            !
!                                                                      !
!    the following are the externally accessable subroutines:          !
!                                                                      !
!       'cld_init'            --- initialization routine               !
!          inputs:                                                     !
!           (si, NLAY, iflip, np3d, icld, me)                          !
!            sashal,crick_proof,ccnorm,norad_precip, me )              !
!          outputs:                                                    !
!           ( none )                                                   !
!                                                                      !
!       'progcld1'           --- zhao/moorthi prognostic cloud scheme  !
!          inputs:                                                     !
!           (plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                   !
!            xlat,xlon,slmsk,                                          !
!            IX, NLAY, NLP1,                                           !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot)                                    !
!                                                                      !
!       'progcld2'           --- ferrier prognostic cloud microphysics !
!          inputs:                                                     !
!           (plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                   !
!            xlat,xlon,slmsk, f_ice,f_rain,r_rime,flgmin,              !
!            IX, NLAY, NLP1,                                           !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot)                                    !
!                                                                      !
!       'progcld8'           --- Thompson prognostic cloud microphysics!
!          inputs:                                                     !
!           (plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,                       !
!            clw, qi, qs, ni, ! cloud water, cloud ice, snow, ice numb !
!            IX, NLAY, NLP1,                                           !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot)                                    !
!                                                                      !
!       'diagcld1'           --- diagnostic cloud calc routine         !
!          inputs:                                                     !
!           (plyr,plvl,tlyr,rhly,vvel,cv,cvt,cvb,                      !
!            xlat,xlon,slmsk,                                          !
!            IX, NLAY, NLP1,                                           !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot)                                    !
!                                                                      !
!    internal accessable only subroutines:                             !
!       'gethml'             --- get diagnostic hi, mid, low clouds    !
!                                                                      !
!       'rhtable'            --- rh lookup table for diag cloud scheme !
!                                                                      !
!                                                                      !
!    cloud array description:                                          !
!                ---  for prognostic cloud: icld=1  ---                !
!          clouds(:,:,1)  -  layer total cloud fraction                !
!          clouds(:,:,2)  -  layer cloud liq water path                !
!          clouds(:,:,3)  -  mean effective radius for liquid cloud    !
!          clouds(:,:,4)  -  layer cloud ice water path                !
!          clouds(:,:,5)  -  mean effective radius for ice cloud       !
!          clouds(:,:,6)  -  layer rain drop water path                !
!          clouds(:,:,7)  -  mean effective radius for rain drop       !
!   **     clouds(:,:,8)  -  layer snow flake water path               !
!          clouds(:,:,9)  -  mean effective radius for snow flake      !
!   ** fu's scheme need to be normalized by snow density (g/m**3/1.0e6)!
!                ---  for diagnostic cloud: icld=0  ---                !
!          clouds(:,:,1)  -  layer total cloud fraction                !
!          clouds(:,:,2)  -  layer cloud optical depth                 !
!          clouds(:,:,3)  -  layer cloud single scattering albedo      !
!          clouds(:,:,4)  -  layer cloud asymmetry factor              !
!                                                                      !
!    external modules referenced:                                      !
!                                                                      !
!       'module machine'             in 'machine.f'                    !
!       'module physcons'            in 'physcons.f'                   !
!       'module module_microphysics' in 'module_bfmicrophysics.f'      !
!                                                                      !
!                                                                      !
!    modification history log:                                         !
!                                                                      !
!       apr 2003,  yu-tai hou                                          !
!                  created 'module_rad_clouds' from combining the      !
!                  original subroutine 'cldjms', 'cldprp', and 'gcljms'!
!       may 2004,  yu-tai hou                                          !
!                  incorporate ferrier's cloud microphysics scheme.    !
!       apr 2005,  yu-tai hou                                          !
!                  modified cloud array and module structures.         !
!       dec 2008,  yu-tai hou                                          !
!                  changed low-cld calc, now cantains clds from sfc    !
!                  layer and upward to the low/mid boundary (include   !
!                  bl-cld). h,m,l clds domain boundaries are adjusted  !
!                  for better agreement with observations.             !
!       jan 2011,  yu-tai hou                                          !
!                  changed virtual temperature as input variable       !
!                  instead of originally computed inside the two       !
!                  prognostic cld schemes 'progcld1' and 'progcld2'.   !
!       apr 2012, yu-tai hou                                           !
!                  modified subroutine cld_init to pass all fixed      !
!                  control variables at the start. and set their       !
!                  correponding module variables for accessing by other!
!                  module subroutines.                                 !
!       May 2012,  brad ferrier, hsin-mu lin                           !
!                  combine the radiation impact of snow and ice to be  !
!                  only ice for progcld2.                              !
!                  The layer cloud ice path is summation of ice & snow.!
!                  The mean effective radius for cloud ice is weighted !
!                  from ice & snow                                     !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!


!=============================================!
      module module_radiation_clouds_nmmb     !
!.............................................!
!
      USE ESMF

      use physparam,           only : icldflg, icmphys, iovrsw, iovrlw, &
     &                                lcrick, lcnorm, lnoprec, lsashal, &
     &                                ivflip, kind_phys, kind_io4
      use physcons,            only : con_fvirt, con_ttp, con_rocp,     &
     &                                con_t0c, con_pi, con_g, con_rd

      use module_mp_thompson,  only : Nt_c, calc_effectRad

!BSF 20120319      use module_microphysics,    only : rsipath2
      use module_iounitdef,    only : NICLTUN
!
      implicit   none
!
      private

!  ---  version tag and last revision date
      character(40), parameter ::                                       &
     &   VTAGCLD='NCEP-Radiation_clouds    v5.1  Nov 2012 '
!    &   VTAGCLD='NCEP-Radiation_clouds    v5.0  Aug 2012 '


!  ---  set constant parameters

      real (kind=kind_phys), parameter :: gfac=1.0e5/con_g              &
     &,                                   gord=con_g/con_rd
      integer, parameter, public :: NF_CLDS = 9   ! number of fields in cloud array
      integer, parameter, public :: NK_CLDS = 3   ! number of cloud vertical domains

!  ---  pressure limits of cloud domain interfaces (low,mid,high) in mb (0.1kPa)
      real (kind=kind_phys), save :: ptopc(NK_CLDS+1,2)

!org  data ptopc / 1050., 642., 350., 0.0,  1050., 750., 500., 0.0 /
      data ptopc / 1050., 650., 400., 0.0,  1050., 750., 500., 0.0 /

!     real (kind=kind_phys), parameter :: climit = 0.01
      real (kind=kind_phys), parameter :: climit = 0.001, climit2=0.05
      real (kind=kind_phys), parameter :: ovcst  = 1.0 - 1.0e-8

!  ---  set default quantities as parameters (for prognostic cloud)

      real (kind=kind_phys), parameter :: reliq_def = 10.0    ! default liq radius to 10 micron
      real (kind=kind_phys), parameter :: reice_def = 50.0    ! default ice radius to 50 micron
      real (kind=kind_phys), parameter :: rrain_def = 1000.0  ! default rain radius to 1000 micron
      real (kind=kind_phys), parameter :: rsnow_def = 250.0   ! default snow radius to 250 micron

!--- set default quantities for new version of progcld2

      real (kind=kind_phys), parameter :: cclimit = 0.001, cclimit2=0.05
      real (kind=kind_phys), parameter :: recwat_def = 5.0    ! default liq radius to 5 microns at <0C
      real (kind=kind_phys), parameter :: recice_def = 10.0   ! default ice radius to 10 microns
      real (kind=kind_phys), parameter :: rerain_def = 100.0  ! default rain radius to 100 microns
      real (kind=kind_phys), parameter :: resnow_def = 50.0   ! default snow radius to 50 microns

!  ---  set look-up table dimensions and other parameters (for diagnostic cloud)

      integer, parameter :: NBIN=100     ! rh in one percent interval
      integer, parameter :: NLON=2       ! =1,2 for eastern and western hemispheres
      integer, parameter :: NLAT=4       ! =1,4 for 60n-30n,30n-equ,equ-30s,30s-60s
      integer, parameter :: MCLD=4       ! =1,4 for bl,low,mid,hi cld type
      integer, parameter :: NSEAL=2      ! =1,2 for land,sea

      real (kind=kind_phys), parameter :: cldssa_def = 0.99   ! default cld single scat albedo
      real (kind=kind_phys), parameter :: cldasy_def = 0.84   ! default cld asymmetry factor

!  ---  xlabdy: lat bndry between tuning regions, +/- xlim for transition
!       xlobdy: lon bndry between tuning regions
      real (kind=kind_phys), parameter :: xlim=5.0
      real (kind=kind_phys)            :: xlabdy(3), xlobdy(3)

      data xlabdy / 30.0,  0.0, -30.0 /,  xlobdy / 0.0, 180., 360. /

!  ---  low cloud vertical velocity adjustment boundaries in mb/sec
      real (kind=kind_phys), parameter :: vvcld1= 0.0003e0
      real (kind=kind_phys), parameter :: vvcld2=-0.0005e0

!  ---  those data will be set up by "cld_init"
!       rhcl : tuned rh relation table for diagnostic cloud scheme

      real (kind=kind_phys) :: rhcl(NBIN,NLON,NLAT,MCLD,NSEAL)

      integer  :: llyr   = 2           ! upper limit of boundary layer clouds
      integer  :: iovr   = 1           ! cloud over lapping method for diagnostic 3-domain
                                       ! output calc (see iovrsw/iovrlw description)

      public progcld1, progcld2, progcld8, diagcld1, cld_init


! =================
      contains
! =================


!-----------------------------------
      subroutine cld_init                                                &
!...................................

!  ---  inputs:
     &     ( si, NLAY, me )
!  ---  outputs:
!          ( none )

!  ===================================================================  !
!                                                                       !
! abstract: cld_init is an initialization program for cloud-radiation   !
!   calculations. it sets up boundary layer cloud top.                  !
!                                                                       !
!                                                                       !
! inputs:                                                               !
!   si (L+1)        : model vertical sigma layer interface              !
!   NLAY            : vertical layer number                             !
!   icld            : cloud computation method flag                     !
!                     =0: model use diagnostic cloud method             !
!                     =1: model use prognostic cloud method             !
!   icmphys         : =3: ferrier microphysics cloud scheme             !
!                     =4: zhao/carr/sundqvist microphysics cloud        !
!                     =8: Thompson microphysics cloud                   !
!   iovrsw/iovrlw   : sw/lw control flag for cloud overlap              !
!                     =0 random overlapping clouds                      !
!                     =1 max/ran overlapping clouds                     !
!   iflip           : control flag for direction of vertical index      !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   sashal          : control flag for shallow convection               !
!   crick_proof     : control flag for eliminating CRICK                !
!   ccnorm          : control flag for in-cloud condensate mixing ratio !
!   norad_precip    : control flag for not using precip in radiation    !
!   me              : print control flag                                !
!                                                                       !
!  outputs: (none)                                                      !
!           to module variables                                         !
!                                                                       !
!  usage:       call cld_init                                            !
!                                                                       !
!  subroutines called:    rhtable                                       !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NLAY, me

      real (kind=kind_phys), intent(in) :: si(:)

!  ---  outputs: (none)

!  ---  locals:
      integer :: k, kl, ier

!
!===> ...  begin here
!
!  ---  set up module variables

      iovr    = max( iovrsw, iovrlw )    !cld ovlp used for diag HML cld output

      if (me == 0) print *, VTAGCLD      !print out version tag

      if ( icldflg == 0 ) then
        if (me == 0) print *,' - Using Diagnostic Cloud Method'

!  ---  set up tuned rh table

        call rhtable( me, ier )

        if (ier < 0) then
          write(6,99) ier
  99      format(3x,' *** Error in finding tuned RH table ***'          &
     &,         /3x,'     STOP at calling subroutine RHTABLE !!'/)
          stop 99
        endif
      else
        if (me == 0) then
          if (icmphys == 3) then
             print *,' - Using Prognostic Cloud Method'
             print *, '   --- Ferrier cloud microphysics'
          elseif (icmphys == 4) then
             print *,' - Using Prognostic Cloud Method'
             print *, '   --- Zhao/Carr/Sundqvist microphysics'
          elseif (icmphys == 5) then
             print *,' - Using Diagnostic Cloud Method'
             print *, '   --- Old GFDL cloud'
          elseif (icmphys == 8) then
             print *,' - Using Prognostic Cloud Method'
             print *, '   --- Thompson microphysics'
          else
            print *,'  !!! ERROR in cloud microphysc specification!!!', &
     &              '  icmphys (NP3D) =',icmphys
            ! stop
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          endif
        endif
      endif

!  ---  compute llyr - the top of bl cld and is topmost non cld(low) layer
!       for stratiform (at or above lowest 0.1 of the atmosphere)

      if ( ivflip == 0 ) then      ! data from toa to sfc
        lab_do_k0 : do k = NLAY, 2, -1
          kl = k
          if (si(k) < 0.9e0) exit lab_do_k0
        enddo  lab_do_k0

        llyr = kl 

      else                      ! data from sfc to top
        lab_do_k1 : do k = 2, NLAY
          kl = k
          if (si(k) < 0.9e0) exit lab_do_k1
        enddo  lab_do_k1

        llyr = kl - 1

      endif                     ! end_if_ivflip

!
      return
!...................................
      end subroutine cld_init
!-----------------------------------


!-----------------------------------
      subroutine progcld1                                               &
!...................................

!  ---  inputs:
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                    &
     &       xlat,xlon,slmsk,                                           &
     &       IX, NLAY, NLP1,                                            &
!  ---  outputs:
     &       clouds,clds,mtop,mbot                                      &
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld1    computes cloud related quantities using    !
!   zhao/moorthi's prognostic cloud microphysics scheme.                !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                       !
!                                                                       !
! program history log:                                                  !
!      11-xx-1992   y.h., k.a.c, a.k. - cloud parameterization          !
!         'cldjms' patterned after slingo and slingo's work (jgr,       !
!         1992), stratiform clouds are allowed in any layer except      !
!         the surface and upper stratosphere. the relative humidity     !
!         criterion may cery in different model layers.                 !
!      10-25-1995   kenneth campana   - tuned cloud rh curves           !
!         rh-cld relation from tables created using mitchell-hahn       !
!         tuning technique on airforce rtneph observations.             !
!      11-02-1995   kenneth campana   - the bl relationships used       !
!         below llyr, except in marine stratus regions.                 !
!      04-11-1996   kenneth campana   - save bl cld amt in cld(,5)      !
!      12-29-1998   s. moorthi        - prognostic cloud method         !
!      04-15-2003   yu-tai hou        - rewritten in frotran 90         !
!         modulized form, seperate prognostic and diagnostic methods    !
!         into two packages.                                            !
!                                                                       !
! usage:         call progcld1                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   tvly  (IX,NLAY) : model layer virtual temperature in k              !
!   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
!   qstl  (IX,NLAY) : layer saturate humidity in gm/gm                  !
!   rhly  (IX,NLAY) : layer relative humidity (=qlyr/qstl)              !
!   clw   (IX,NLAY) : layer cloud condensate amount                     !
!   xlat  (IX)      : grid latitude in radians                          !
!   xlon  (IX)      : grid longitude in radians  (not used)             !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!   iflip           : control flag for in/out vertical indexing         !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   iovr            : control flag for cloud overlap                    !
!                     =0 random overlapping clouds                      !
!                     =1 max/ran overlapping clouds                     !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lsashal         : control flag for shallow convection               !
!   lcrick          : control flag for eliminating CRICK                !
!   lcnorm          : control flag for in-cld condensate                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, tvly, qlyr, qstl, rhly, clw

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldtot, cldcnv,      &
     &       cwp, cip, crp, csp, rew, rei, res, rer, delp, tem2d, clwf

      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1)

      real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,       &
     &       tem1, tem2, tem3


      integer :: i, k, id, nf


!
!===> ... begin here
!
      do nf=1,nf_clds
        do k=1,nlay
          do i=1,ix
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo
!     clouds(:,:,:) = 0.0

      do k = 1, NLAY
        do i = 1, IX
          cldtot(i,k) = 0.0
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          rew   (i,k) = reliq_def            ! default liq radius
          rei   (i,k) = reice_def            ! default ice radius
          rer   (i,k) = rrain_def            ! default rain radius
          res   (i,k) = rsnow_def            ! default snow radius
          tem2d (i,k) = min( 1.0, max( 0.0, (con_ttp-tlyr(i,k))*0.05 ) )
          clwf(i,k)   = 0.0
        enddo
      enddo
!
      if ( lcrick ) then
        do i = 1, IX
          clwf(i,1)    = 0.75*clw(i,1)    + 0.25*clw(i,2)
          clwf(i,nlay) = 0.75*clw(i,nlay) + 0.25*clw(i,nlay-1)
        enddo
        do k = 2, NLAY-1
          do i = 1, IX
            clwf(i,K) = 0.25*clw(i,k-1) + 0.5*clw(i,k) + 0.25*clw(i,k+1)
          enddo
        enddo
      else
        do k = 1, NLAY
          do i = 1, IX
            clwf(i,k) = clw(i,k)
          enddo
        enddo
      endif

!  ---  find top pressure for each cloud domain for given latitude
!       ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          tem2 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
!         tem2 = 0.5 - xlat(i)/con_pi    ! if xlat in 0 -> pi range

          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*abs(tem2)-1.0 )
        enddo
      enddo

!  ---  compute liquid/ice condensate path in g/m**2

      if ( ivflip == 0 ) then          ! input data from toa to sfc
        do k = 1, NLAY
          do i = 1, IX
            delp(i,k) = plvl(i,k+1) - plvl(i,k)
            clwt     = max(0.0, clwf(i,k)) * gfac * delp(i,k)
            cip(i,k) = clwt * tem2d(i,k)
            cwp(i,k) = clwt - cip(i,k)
          enddo
        enddo
      else                             ! input data from sfc to toa
        do k = 1, NLAY
          do i = 1, IX
            delp(i,k) = plvl(i,k) - plvl(i,k+1)
            clwt     = max(0.0, clwf(i,k)) * gfac * delp(i,k)
            cip(i,k) = clwt * tem2d(i,k)
            cwp(i,k) = clwt - cip(i,k)
          enddo
        enddo
      endif                            ! end_if_ivflip

!  ---  effective liquid cloud droplet radius over land

      do i = 1, IX
        if (nint(slmsk(i)) == 1) then
          do k = 1, NLAY
            rew(i,k) = 5.0 + 5.0 * tem2d(i,k)
          enddo
        endif
      enddo

!  ---  layer cloud fraction

      if ( ivflip == 0 ) then              ! input data from toa to sfc

        clwmin = 0.0
        if (.not. lsashal) then
        do k = NLAY, 1, -1
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)
!           clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then

              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

              tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
              tem1  = 2000.0 / tem1
!             tem1  = 1000.0 / tem1

              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
        enddo
        else
        do k = NLAY, 1, -1
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)
!           clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then

              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

!             tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
!             tem1  = 2000.0 / tem1

              tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)  !jhan
              tem1  = 100.0 / tem1
!
!             tem1  = 2000.0 / tem1
!             tem1  = 1000.0 / tem1
!

              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )
              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
        enddo
        endif

      else                                 ! input data from sfc to toa

        clwmin = 0.0
        if (.not. lsashal) then
        do k = 1, NLAY
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)
!           clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then

              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

              tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
              tem1  = 2000.0 / tem1

!             tem1  = 1000.0 / tem1

              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
        enddo
        else
                do k = 1, NLAY
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)
!           clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then

              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

!             tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
!             tem1  = 2000.0 / tem1

              tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)  !jhan
              tem1  = 100.0 / tem1
!
!             tem1  = 2000.0 / tem1
!             tem1  = 1000.0 / tem1
!

              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
        enddo
        endif

      endif                                ! end_if_flip

      do k = 1, NLAY
        do i = 1, IX
          if (cldtot(i,k) < climit) then
            cldtot(i,k) = 0.0
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo
!     where (cldtot < climit)
!       cldtot = 0.0
!       cwp    = 0.0
!       cip    = 0.0
!       crp    = 0.0
!       csp    = 0.0
!     endwhere
!
      if ( lcnorm ) then
        do k = 1, NLAY
          do i = 1, IX
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!  ---  effective ice cloud droplet radius

      do k = 1, NLAY
        do i = 1, IX
          tem2 = tlyr(i,k) - con_ttp

          if (cip(i,k) > 0.0) then
            tem3 = gord * cip(i,k) * (plyr(i,k)/delp(i,k)) / tvly(i,k)

            if (tem2 < -50.0) then
              rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
            elseif (tem2 < -40.0) then
              rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
            elseif (tem2 < -30.0) then
              rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
            else
              rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
            endif
            rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(10.0, min(rei(i,k), 100.0))
          endif
        enddo
      enddo

!
      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
!         clouds(i,k,6) = 0.0
          clouds(i,k,7) = rer(i,k)
!         clouds(i,k,8) = 0.0
          clouds(i,k,9) = rei(i,k)
        enddo
      enddo


!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which is
!  ---  also used by the lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       IX,NLAY,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld1
!-----------------------------------


!-----------------------------------
      subroutine progcld2                                               &
!...................................

!  ---  inputs:
     &     ( plyr, xlat, IX, NLAY, cldtot,                              &
!  ---  outputs:
     &       clds,mtop,mbot                                             &
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld2    computes cloud related quantities using    !
!   ferrier's prognostic cloud microphysics scheme.                     !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                       !
!                                                                       !
! program history log:                                                  !
!        -  -       brad ferrier      - original development            !
!        -  -2003   s. moorthi        - adapted to ncep gfs model       !
!      05-05-2004   yu-tai hou        - rewritten as a separated        !
!                   program in the cloud module.                        !
!                                                                       !
! usage:         call progcld2                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   tvly  (IX,NLAY) : model layer virtual temperature in k              !
!   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
!   qstl  (IX,NLAY) : layer saturate humidity in gm/gm                  !
!   rhly  (IX,NLAY) : layer relative humidity (=qlyr/qstl)              !
!   clw   (IX,NLAY) : layer cloud condensate amount                     !
!   f_ice (IX,NLAY) : fraction of layer cloud ice  (ferrier micro-phys) !
!   f_rain(IX,NLAY) : fraction of layer rain water (ferrier micro-phys) !
!   r_rime(IX,NLAY) : mass ratio of total ice to unrimed ice (>=1)      !
!   xlat  (IX)      : grid latitude in radians                          !
!   xlon  (IX)      : grid longitude in radians  (not used)             !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!   iflip           : control flag for in/out vertical indexing         !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   iovr            : control flag for cloud overlap                    !
!                     =0 random overlapping clouds                      !
!                     =1 max/ran overlapping clouds                     !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         (g/m**2)      !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        (g/m**2)      !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lsashal         : control flag for shallow convection               !
!   lcrick          : control flag for eliminating CRICK                !
!   lcnorm          : control flag for in-cld condensate                !
!   lnoprec         : control flag for no precip in rad                 !
!                                                                       !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs

      real (kind=kind_phys), dimension(:,:), intent(in) :: plyr, cldtot
      real (kind=kind_phys), dimension(:),   intent(in) :: xlat

      integer,  intent(in) :: IX, NLAY

!  ---  outputs

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds
      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:

      real (kind=kind_phys) :: cldcnv(IX,NLAY), ptop1(IX,NK_CLDS+1)
      real (kind=kind_phys) :: tem1, tem2, tem3

      integer :: i, k, id

!
!===> ... begin here
!
!
      do k = 1, NLAY
        do i = 1, IX
          cldcnv(i,k) = 0.0
        enddo
      enddo


!  ---  find top pressure for each cloud domain for given latitude
!       ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          tem3 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
          tem2 = max( 0.0, 4.0*abs(tem3)-1.0 )
          ptop1(i,id) = ptopc(id,1) + tem1*tem2
        enddo
      enddo

!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       IX,NLAY,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld2
!-----------------------------------


!+---+-----------------------------------------------------------------+

      subroutine progcld8 (plyr, plvl, tlyr, qlyr, qc, qi, qs, ni,      &
     &                     xlat,IX, NLAY, NLP1,                         &
     &                     cldfk1, cld_fraction,                        &
     &                     clouds,clds,mtop,mbot)


! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld8    computes cloud related quantities using    !
!   Thompson's prognostic cloud microphysics scheme.                    !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice/snow radiative effective radius. !
!                                                                       !
! program history log:                                                  !
!      22Jan2014    G. Thompson                                         !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
!   qc    (IX,NLAY) : layer cloud water mixing ratio (kg/kg)            !
!   qi    (IX,NLAY) : layer cloud ice mixing ratio (kg/kg)              !
!   qs    (IX,NLAY) : layer snow mixing ratio (kg/kg)                   !
!   ni    (IX,NLAY) : layer cloud ice number concentration (#/kg)       !
!   xlat  (IX)      : grid latitude in radians                          !
!   IX              : horizontal dimension                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!   iflip           : control flag for in/out vertical indexing         !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   iovr            : control flag for cloud overlap                    !
!                     =0 random overlapping clouds                      !
!                     =1 max/ran overlapping clouds                     !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(I,K,1) - layer total cloud fraction                       !
!      clouds(I,K,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(I,K,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(I,K,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(I,K,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(I,K,6) - layer rain drop water path         (g/m**2)      !
!      clouds(I,K,7) - mean eff radius for rain drop      (micron)      !
!      clouds(I,K,8) - layer snow flake water path        (g/m**2)      !
!      clouds(I,K,9) - mean eff radius for snow flake     (micron)      !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lsashal         : control flag for shallow convection               !
!   lcrick          : control flag for eliminating CRICK                !
!   lcnorm          : control flag for in-cld condensate                !
!   lnoprec         : control flag for no precip in rad                 !
!                                                                       !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  constants
      real (kind=kind_phys), parameter :: EPSQ = 1.0e-12

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1, cld_fraction

      real (kind=kind_phys), dimension(:,:), intent(in) :: plyr, plvl,  &
     &       tlyr, qlyr, qc, qi, qs, ni, cldfk1
      real (kind=kind_phys), dimension(:),   intent(in) :: xlat

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds
      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds
      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      integer:: i, k, nf, id
      real (kind=kind_phys), dimension(IX,NLAY) :: cldtot, cldcnv
      double precision, dimension(NLAY):: re_qc1d, re_qi1d, re_qs1d
      double precision, dimension(NLAY):: t1d, p1d, qv1d, qc1d,         &
     &                        nc1d, qi1d, ni1d, qs1d
      real (kind=kind_phys):: tem1, tem2
      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1)
      real (kind=kind_phys) :: rho_k, cpath, snow_mass_factor
      logical:: is_any_clouds

!+---+

      do nf=1,NF_CLDS
        do k = 1, NLAY
          do i = 1, IX
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo
      do k = 1, 5
        do i = 1, IX
          clds(i,k) = 0.0
        enddo
      enddo
      do k = 1, NK_CLDS
        do i = 1, IX
          mbot(i,k) = 0.0
          mtop(i,k) = 0.0
        enddo
      enddo

! --- choice of cloud fraction

      if (cld_fraction==0) then   ! non Thompson cloud fraction (default)
        do k = 1, NLAY
          do i = 1, IX
            cldtot(i,k) = 0.0          ! original
          enddo
        enddo
      else                        ! Thompson cloud fraction
        do k = 1, NLAY
          do i = 1, IX
            cldtot(i,k) = cldfk1(i,k)   ! for Thompson cloud fraction
          enddo
        enddo
      endif

      do k = 1, NLAY
        do i = 1, IX
          cldcnv(i,k) = 0.0
        enddo
      enddo

      do i = 1, IX
        is_any_clouds = .false.
        do k = 1, NLAY
          t1d(k) = tlyr(i,k)
          p1d(k) = plyr(i,k)*100.
          qv1d(k) = MAX(1.E-8, qlyr(i,k)/(1.-qlyr(i,k)))
          qc1d(k) = qc(i,k)
          nc1d(k) = Nt_c
          qi1d(k) = qi(i,k)
          ni1d(k) = ni(i,k)
          qs1d(k) = qs(i,k)
          re_qc1d(k) = 2.501D-6
          re_qi1d(k) = 10.01D-6
          re_qs1d(k) = 25.01D-6
          if ((qc1d(k)+qi1d(k)+qs1d(k)).gt.1.E-9) then
            is_any_clouds = .true.
            cldtot(i,k) = 1.0
          endif
        enddo

        if (is_any_clouds) then
          call calc_effectRad (t1d, p1d, qv1d, qc1d, nc1d, qi1d,        &
     &                ni1d, qs1d, re_qc1d, re_qi1d, re_qs1d, 1, NLAY)

!     do k = 1, NLAY
!       if(cldtot(i,k).gt.0.) then
!        write(*,*) ' DEBUG_GT, K, qc, re_qc, qi, re_qi, qs, re_qs ',   &
!    &   k,qc1d(k)*1000.,re_qc1d(k)*1.E6,qi1d(k)*1000.,re_qi1d(k)*1.E6,qs1d(k)*1000.,re_qs1d(k)*1.E6
!       endif
!     enddo

          do k = 1, NLAY
            clouds(i,k,1) = cldtot(i,k)
            rho_k = 0.622*p1d(k)/(287.05*t1d(k)*(qv1d(k)+0.622))
            cpath = ABS(plvl(i,k+1)-plvl(i,k))*(100000.0/(rho_k*9.8))
            clouds(i,k,2) = max(0.0, qc(i,k)) * cpath
            clouds(i,k,3) = max(2.501, re_qc1d(k)*1.E6)
            clouds(i,k,4) = max(0.0, qi(i,k)) * cpath
            clouds(i,k,5) = max(10.01, re_qi1d(k)*1.E6)
            clouds(i,k,8) = max(0.0, qs(i,k)) * cpath
            clouds(i,k,9) = max(25.01, re_qs1d(k)*1.E6)
            if (clouds(i,k,8).gt.0.0 .AND. clouds(i,k,9).gt.130.) then
              snow_mass_factor = (130.0E-6/re_qs1d(k))*(130.0E-6/re_qs1d(k))
              clouds(i,k,8) = snow_mass_factor*qs(i,k) * cpath
              clouds(i,k,9) = 130.0
            endif
          enddo
        endif
      enddo


!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

!  ---  find top pressure for each cloud domain for given latitude
!       ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)
!       data ptopc / 1050., 650., 400., 0.0,  1050., 750., 500., 0.0 /

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          tem2 = max( 0.0, 4.0*abs(xlat(i))/con_pi-1.0 )
          ptop1(i,id) = ptopc(id,1) + tem1*tem2
        enddo
      enddo

      call gethml (plyr, ptop1, cldtot, cldcnv, IX,NLAY,                &
     &       clds, mtop, mbot)


      return

      end subroutine progcld8

!+---+-----------------------------------------------------------------+


!-----------------------------------
      subroutine diagcld1                                               &
!...................................

!  ---  inputs:
     &     ( plyr,plvl,tlyr,rhly,vvel,cv,cvt,cvb,                       &
     &       xlat,xlon,slmsk,                                           &
     &       IX, NLAY, NLP1,                                            &
!  ---  outputs:
     &       clouds,clds,mtop,mbot                                      &
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    diagcld1    computes cloud fractions for radiation     !
!   calculations.                                                       !
!                                                                       !
! abstract:  clouds are diagnosed from layer relative humidity, and     !
!   estimate cloud optical depth from temperature and layer thickness.  !
!   then computes the low, mid, high, total and boundary layer cloud    !
!   fractions and the vertical indices of low, mid, and high cloud top  !
!   and base.  the three vertical cloud domains are set up in the       !
!   initial subroutine "cld_init".                                       !
!                                                                       !
! program history log:                                                  !
!      11-xx-1992   y.h., k.a.c, a.k. - cloud parameterization          !
!         'cldjms' patterned after slingo and slingo's work (jgr,       !
!         1992), stratiform clouds are allowed in any layer except      !
!         the surface and upper stratosphere. the relative humidity     !
!         criterion may cery in different model layers.                 !
!      10-25-1995   kenneth campana   - tuned cloud rh curves           !
!         rh-cld relation from tables created using mitchell-hahn       !
!         tuning technique on airforce rtneph observations.             !
!      11-02-1995   kenneth campana   - the bl relationships used       !
!         below llyr, except in marine stratus regions.                 !
!      04-11-1996   kenneth campana   - save bl cld amt in cld(,5)      !
!      12-29-1998   s. moorthi        - prognostic cloud method         !
!      04-15-2003   yu-tai hou        - rewritten in frotran 90         !
!         modulized form, seperate prognostic and diagnostic methods    !
!         into two packages.                                            !
!                                                                       !
! usage:         call diagcld1                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   rhly  (IX,NLAY) : layer relative humidity                           !
!   vvel  (IX,NLAY) : layer mean vertical velocity in mb/sec            !
!   clw   (IX,NLAY) : layer cloud condensate amount         (not used)  !
!   xlat  (IX)      : grid latitude in radians                          !
!   xlon  (IX)      : grid longitude in radians                         !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   cv    (IX)      : fraction of convective cloud                      !
!   cvt, cvb (IX)   : conv cloud top/bottom pressure in mb              !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!   iflip           : control flag for in/out vertical indexing         !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   iovr            : control flag for cloud overlap                    !
!                     =0 random overlapping clouds                      !
!                     =1 max/ran overlapping clouds                     !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud optical depth                        !
!      clouds(:,:,3) - layer cloud single scattering albedo             !
!      clouds(:,:,4) - layer cloud asymmetry factor                     !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1   ! , iflip, iovr

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, rhly, vvel

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk, cv, cvt, cvb

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldtot, cldcnv,      &
     &       cldtau, taufac, dthdp, tem2d

      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1)

      real (kind=kind_phys) :: cr1, cr2, crk, pval, cval, omeg, value,  &
     &       tem1, tem2

      integer, dimension(IX):: idom, kcut

!  ---  for rh-cl calculation
      real (kind=kind_phys) :: xlatdg, xlondg, xlnn, xlss, xrgt, xlft,  &
     &       rhcla(NBIN,NLON,MCLD,NSEAL),  rhcld(IX,NBIN,MCLD)

      integer :: ireg, ib, ic, id, id1, il, is, nhalf

      integer :: i, j, k, klowt, klowb

      logical :: notstop

!
!===> ... begin here
!
      clouds(:,:,:) = 0.0

      tem1 = 180.0 / con_pi

      lab_do_i_IX : do i = 1, IX

        xlatdg = xlat(i) * tem1                    ! if xlat in pi/2 -> -pi/2 range
!       xlatdg = 90.0 - xlat(i)*tem1               ! if xlat in 0 -> pi range

        xlondg = xlon(i) * tem1
        if (xlondg < 0.0) xlondg = xlondg + 360.0  ! if in -180->180, chg to 0->360 range

        ireg = 4

!  ---  get rh-cld relation for this lat

        lab_do_j : do j = 1, 3
          if (xlatdg > xlabdy(j)) then
            ireg = j
            exit lab_do_j
          endif
        enddo  lab_do_j

        do is = 1, NSEAL
          do ic = 1, MCLD
            do il = 1, NLON
              do ib = 1, NBIN
                rhcla(ib,il,ic,is) = rhcl(ib,il,ireg,ic,is)
              enddo
            enddo
          enddo
        enddo

!  ---  linear transition between latitudinal regions...
        do j = 1, 3
          xlnn = xlabdy(j) + xlim
          xlss = xlabdy(j) - xlim

          if (xlatdg < xlnn .and. xlatdg > xlss) then
            do is = 1, NSEAL
              do ic = 1, MCLD
                do il = 1, NLON
                  do ib = 1, NBIN
                    rhcla(ib,il,ic,is) = rhcl(ib,il,j+1,ic,is)          &
     &                + (rhcl(ib,il,j,ic,is)-rhcl(ib,il,j+1,ic,is))     &
     &                * (xlatdg-xlss) / (xlnn-xlss)
                  enddo
                enddo
              enddo
            enddo
          endif

        enddo        ! end_j_loop

!  ---  get rh-cld relationship for each grid point, interpolating
!       longitudinally between regions if necessary..

        if (slmsk(i) < 1.0) then
          is = 2
        else
          is = 1
        endif

!  ---  which hemisphere (e,w)

        if (xlondg > 180.e0) then
          il = 2
        else
          il = 1
        endif

        do ic = 1, MCLD
          do ib = 1, NBIN
            rhcld(i,ib,ic) = rhcla(ib,il,ic,is)
          enddo

          lab_do_k : do k = 1, 3
            tem2 = abs(xlondg - xlobdy(k))

            if (tem2 < xlim) then
              id = il
              id1= id + 1
              if (id1 > NLON) id1 = 1

              xlft = xlobdy(k) - xlim
              xrgt = xlobdy(k) + xlim

              do ib = 1, NBIN
                rhcld(i,ib,ic) = rhcla(ib,id1,ic,is)                    &
     &            + (rhcla(ib,id,ic,is) - rhcla(ib,id1,ic,is))          &
     &            * (xlondg-xrgt)/(xlft-xrgt)
              enddo
              exit lab_do_k
            endif

          enddo  lab_do_k

        enddo   ! end_do_ic_loop
      enddo  lab_do_i_IX

!  ---  find top pressure for each cloud domain

      do j = 1, 4
        tem1 = ptopc(j,2) - ptopc(j,1)

        do i = 1, IX
          tem2 = max( 0.0, 4.0*abs(xlat(i))/con_pi-1.0 )
          ptop1(i,j) = ptopc(j,1) + tem1*tem2
        enddo
      enddo

!  ---  stratiform cloud optical depth

      do k = 1, NLAY
        do i = 1, IX
          tem1 = tlyr(i,k) - con_ttp
          if (tem1 <= -10.0) then
            cldtau(i,k) = max( 0.1e-3, 2.0e-6*(tem1+82.5)**2 )
          else
            cldtau(i,k) = min( 0.08, 6.949e-3*tem1+0.08 )
          endif
        enddo
      enddo

!  ---  potential temperature and its lapse rate

      do k = 1, NLAY
        do i = 1, IX
          cldtot(i,k) = 0.0
          cldcnv(i,k) = 0.0
          tem1        = (plyr(i,k)*0.001) ** (-con_rocp)
          tem2d(i,k)  = tem1 * tlyr(i,k)
        enddo
      enddo

      do k = 1, NLAY-1
        do i = 1, IX
          dthdp(i,k) = (tem2d(i,k+1)-tem2d(i,k))/(plyr(i,k+1)-plyr(i,k))
        enddo
      enddo
!
!===> ... diagnostic method to find cloud amount cldtot, cldcnv
!

      if ( ivflip == 0 ) then                 ! input data from toa to sfc

!  ---  find the lowest low cloud top sigma level, computed for each lat cause
!       domain definition changes with latitude...

!       klowb = 1
        klowt = 1
        do k = 1, NLAY
          do i = 1, IX
!           if (plvl(i,k) < ptop1(i,2))  klowb = k
            if (plvl(i,k) < ptop1(i,2))  klowt = max(klowt,k)
            taufac(i,k) = plvl(i,k+1) - plvl(i,k)
          enddo
        enddo

        do i = 1, IX

!  ---  find the stratosphere cut off layer for high cloud (about 250mb).
!       it is assumed to be above the layer with dthdp less than -0.25 in
!       the high cloud domain

          kcut(i) = 2
          lab_do_kcut0 : do k = klowt-1, 2, -1
            if (plyr(i,k) <= ptop1(i,3) .and.                           &
     &          dthdp(i,k) < -0.25e0) then
              kcut(i) = k
              exit lab_do_kcut0
            endif
          enddo  lab_do_kcut0

!  ---  put convective cloud into 'cldcnv', no merge at this point..

          if (cv(i) >= climit .and. cvt(i) < cvb(i)) then
            id  = NLAY
            id1 = NLAY

            lab_do_k_cvt0 : do k = 2, NLAY
              if (cvt(i) <= plyr(i,k)) then
                id = k - 1
                exit lab_do_k_cvt0
              endif
            enddo  lab_do_k_cvt0

            lab_do_k_cvb0 : do k = NLAY-1, 1, -1
              if (cvb(i) >= plyr(i,k)) then
                id1 = k + 1
                exit lab_do_k_cvb0
              endif
            enddo  lab_do_k_cvb0

            tem1 = plyr(i,id1) - plyr(i,id)
            do k = id, id1
              cldcnv(i,k) = cv(i)
              taufac(i,k) = taufac(i,k) * max( 0.25, 1.0-0.125*tem1 )
              cldtau(i,k) = 0.06
            enddo
          endif

        enddo                ! end_do_i_loop

!  ---  calculate stratiform cloud and put into array 'cldtot' using
!       the cloud-rel.humidity relationship from table look-up..where
!       tables obtained using k.mitchell frequency distribution tuning
!bl       (observations are daily means from us af rtneph).....k.a.c.
!bl       tables created without lowest 10 percent of atmos.....k.a.c.
!      (observations are synoptic using -6,+3 window from rtneph)
!       tables are created with lowest 10-percent-of-atmos, and are
!  ---  now used..  25 october 1995 ... kac.

        do k = NLAY-1, 2, -1

          if (k < llyr) then
            do i = 1, IX
              idom(i) = 0
            enddo

            do i = 1, IX
              lab_do_ic0 : do ic = 2, 4
                if(plyr(i,k) >= ptop1(i,ic)) then
                  idom(i) = ic
                  exit lab_do_ic0
                endif
              enddo  lab_do_ic0
            enddo
          else
            do i = 1, IX
              idom(i) = 1
            enddo
          endif

          do i = 1, IX
            id = idom(i)
            nhalf = (NBIN + 1) / 2

            if (id <= 0 .or. k < kcut(i)) then
              cldtot(i,k) = 0.0
            elseif (rhly(i,k) <= rhcld(i,1,id)) then
              cldtot(i,k) = 0.0
            elseif (rhly(i,k) >= rhcld(i,NBIN,id)) then
              cldtot(i,k) = 1.0
            else
              ib = nhalf
              crk = rhly(i,k)

              notstop = .true.
              do while ( notstop )
                nhalf = (nhalf + 1) / 2
                cr1 = rhcld(i,ib,  id)
                cr2 = rhcld(i,ib+1,id)

                if (crk <= cr1) then
                  ib = max( ib-nhalf, 1 )
                elseif (crk > cr2) then
                  ib = min( ib+nhalf, NBIN-1 )
                else
                  cldtot(i,k) = 0.01 * (ib + (crk - cr1)/(cr2 - cr1))
                  notstop = .false.
                endif
              enddo      ! end_do_while
            endif
          enddo          ! end_do_i_loop

        enddo            ! end_do_k_loop

! --- vertical velocity adjustment on low clouds

        value = vvcld1 - vvcld2
        do k = klowt, llyr+1
          do i = 1, IX

            omeg = vvel(i,k)
            cval = cldtot(i,k)
            pval = plyr(i,k)

! --- vertical velocity adjustment on low clouds

            if (cval >= climit .and. pval >= ptop1(i,2)) then
              if (omeg >= vvcld1) then
                cldtot(i,k) = 0.0
              elseif (omeg > vvcld2) then
                tem1 = (vvcld1 - omeg) / value
                cldtot(i,k) = cldtot(i,k) * sqrt(tem1)
              endif
            endif

          enddo     ! end_do_i_loop
        enddo       ! end_do_k_loop

      else                                    ! input data from sfc to toa

!  ---  find the lowest low cloud top sigma level, computed for each lat cause
!       domain definition changes with latitude...

!       klowb = NLAY
        klowt = NLAY
        do k = NLAY, 1, -1
          do i = 1, IX
!           if (plvl(i,k) < ptop1(i,2))  klowb = k
            if (plvl(i,k) < ptop1(i,2))  klowt = min(klowt,k)
            taufac(i,k) = plvl(i,k) - plvl(i,k+1)       ! dp for later cal cldtau use
          enddo
        enddo

        do i = 1, IX

!  ---  find the stratosphere cut off layer for high cloud (about 250mb).
!       it is assumed to be above the layer with dthdp less than -0.25 in
!       the high cloud domain

          kcut(i) = NLAY - 1
          lab_do_kcut1 : do k = klowt+1, NLAY-1
            if (plyr(i,k) <= ptop1(i,3) .and.                           &
     &          dthdp(i,k) < -0.25e0) then
              kcut(i) = k
              exit lab_do_kcut1
            endif
          enddo  lab_do_kcut1

!  ---  put convective cloud into 'cldcnv', no merge at this point..

          if (cv(i) >= climit .and. cvt(i) < cvb(i)) then
            id  = 1
            id1 = 1

            lab_do_k_cvt : do k = NLAY-1, 1, -1
              if (cvt(i) <= plyr(i,k)) then
                id = k + 1
                exit lab_do_k_cvt
              endif
            enddo  lab_do_k_cvt

            lab_do_k_cvb : do k = 2, NLAY
              if (cvb(i) >= plyr(i,k)) then
                id1 = k - 1
                exit lab_do_k_cvb
              endif
            enddo  lab_do_k_cvb

            tem1 = plyr(i,id1) - plyr(i,id)
            do k = id1, id
              cldcnv(i,k) = cv(i)
              taufac(i,k) = taufac(i,k) * max( 0.25, 1.0-0.125*tem1 )
              cldtau(i,k) = 0.06
            enddo
          endif

        enddo     ! end_do_i_loop

!  ---  calculate stratiform cloud and put into array 'cldtot' using
!       the cloud-rel.humidity relationship from table look-up..where
!       tables obtained using k.mitchell frequency distribution tuning
!bl       (observations are daily means from us af rtneph).....k.a.c.
!bl       tables created without lowest 10 percent of atmos.....k.a.c.
!      (observations are synoptic using -6,+3 window from rtneph)
!       tables are created with lowest 10-percent-of-atmos, and are
!  ---  now used..  25 october 1995 ... kac.

        do k = 2, NLAY-1

          if (k > llyr) then
            do i = 1, IX
              idom(i) = 0
            enddo

            do i = 1, IX
              lab_do_ic1 : do ic = 2, 4
                if(plyr(i,k) >= ptop1(i,ic)) then
                  idom(i) = ic
                  exit lab_do_ic1
                endif
              enddo  lab_do_ic1
            enddo
          else
            do i = 1, IX
              idom(i) = 1
            enddo
          endif

          do i = 1, IX
            id = idom(i)
            nhalf = (NBIN + 1) / 2

            if (id <= 0 .or. k > kcut(i)) then
              cldtot(i,k) = 0.0
            elseif (rhly(i,k) <= rhcld(i,1,id)) then
              cldtot(i,k) = 0.0
            elseif (rhly(i,k) >= rhcld(i,NBIN,id)) then
              cldtot(i,k) = 1.0
            else
              ib = nhalf
              crk = rhly(i,k)

              notstop = .true.
              do while ( notstop )
                nhalf = (nhalf + 1) / 2
                cr1 = rhcld(i,ib,  id)
                cr2 = rhcld(i,ib+1,id)

                if (crk <= cr1) then
                  ib = max( ib-nhalf, 1 )
                elseif (crk > cr2) then
                  ib = min( ib+nhalf, NBIN-1 )
                else
                  cldtot(i,k) = 0.01 * (ib + (crk - cr1)/(cr2 - cr1))
                  notstop = .false.
                endif
              enddo      ! end_do_while
            endif
          enddo          ! end_do_i_loop

        enddo            ! end_do_k_loop

! --- vertical velocity adjustment on low clouds

        value = vvcld1 - vvcld2
        do k = llyr-1, klowt
          do i = 1, IX

            omeg = vvel(i,k)
            cval = cldtot(i,k)
            pval = plyr(i,k)

! --- vertical velocity adjustment on low clouds

            if (cval >= climit .and. pval >= ptop1(i,2)) then
              if (omeg >= vvcld1) then
                cldtot(i,k) = 0.0
              elseif (omeg > vvcld2) then
                tem1 = (vvcld1 - omeg) / value
                cldtot(i,k) = cldtot(i,k) * sqrt(tem1)
              endif
            endif

          enddo     ! end_do_i_loop
        enddo       ! end_do_k_loop

      endif                                   ! end_if_ivflip

!  ---  diagnostic cloud optical depth
!     cldtau = cldtau * taufac

      where (cldtot < climit)
        cldtot = 0.0
      endwhere
      where (cldcnv < climit)
        cldcnv = 0.0
      endwhere

      where (cldtot < climit .and. cldcnv < climit)
        cldtau = 0.0
      endwhere

      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = max(cldtot(i,k), cldcnv(i,k))
          clouds(i,k,2) = cldtau(i,k) * taufac(i,k)
          clouds(i,k,3) = cldssa_def
          clouds(i,k,4) = cldasy_def
        enddo
      enddo

!
!===> ... compute low, mid, high, total, and boundary layer cloud fractions
!         and clouds top/bottom layer indices for low, mid, and high clouds.
!         the three cloud domain boundaries are defined by ptopc.  the cloud
!         overlapping method is defined by control flag 'iovr', which is
!         also used by the lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       IX, NLAY,                                                  & 
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )

!
      return
!...................................
      end subroutine diagcld1
!-----------------------------------


!-----------------------------------                                    !
      subroutine gethml                                                 &
!...................................                                    !

!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       IX, NLAY,                                                  &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )

!  ===================================================================  !
!                                                                       !
! abstract: compute high, mid, low, total, and boundary cloud fractions !
!   and cloud top/bottom layer indices for model diagnostic output.     !
!   the three cloud domain boundaries are defined by ptopc.  the cloud  !
!   overlapping method is defined by control flag 'iovr', which is also !
!   used by lw and sw radiation programs.                               !
!                                                                       !
! program history log:                                                  !
!      04-29-2004   yu-tai hou        - separated to become individule  !
!         subprogram to calculate averaged h,m,l,bl cloud amounts.      !
!                                                                       !
! usage:         call gethml                                            !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   ptop1 (IX,4)    : pressure limits of cloud domain interfaces        !
!                     (sfc,low,mid,high) in mb (100Pa)                  !
!   cldtot(IX,NLAY) : total or straiform cloud profile in fraction      !
!   cldcnv(IX,NLAY) : convective cloud (for diagnostic scheme only)     !
!   IX              : horizontal dimention                              !
!   NLAY            : vertical layer dimensions                         !
!   iflip           : control flag for in/out vertical indexing         !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   iovr            : control flag for cloud overlap                    !
!                     =0 random overlapping clouds                      !
!                     =1 max/ran overlapping clouds                     !
!                                                                       !
! output variables:                                                     !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   iovr            : control flag for cloud overlap                    !
!                     =0 random overlapping clouds                      !
!                     =1 max/ran overlapping clouds                     !
!                                                                       !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none!

!  ---  inputs:
      integer, intent(in) :: IX, NLAY

      real (kind=kind_phys), dimension(:,:), intent(in) :: plyr, ptop1, &
     &       cldtot, cldcnv

!  ---  outputs
      real (kind=kind_phys), dimension(:,:), intent(out) :: clds

      integer,               dimension(:,:), intent(out) :: mtop, mbot

!  ---  local variables:
      real (kind=kind_phys) :: cl1(IX), cl2(IX)

      real (kind=kind_phys) :: pcur, pnxt, ccur, cnxt

      integer, dimension(IX):: idom, kbt1, kth1, kbt2, kth2

      integer :: i, k, id, id1, kstr, kend, kinc, n_clds

!
!===> ... begin here
!
      do i = 1, IX
         do n_clds = 1, 5
            clds(i, n_clds) = 0.0
         enddo
      enddo

      ! clds(:,:) = 0.0

      do i = 1, IX
        cl1 (i) = 1.0
        cl2 (i) = 1.0
      enddo

!  ---  total and bl clouds, where cl1, cl2 are fractions of clear-sky view
!       layer processed from surface and up

      if ( ivflip == 0 ) then                   ! input data from toa to sfc
        kstr = NLAY
        kend = 1
        kinc = -1
      else                                      ! input data from sfc to toa
        kstr = 1
        kend = NLAY
        kinc = 1
      endif                                     ! end_if_ivflip

      if ( iovr == 0 ) then                     ! random overlap

        do k = kstr, kend, kinc
          do i = 1, IX
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
            if (ccur >= climit) cl1(i) = cl1(i) * (1.0 - ccur)
          enddo

          if (k == llyr) then
            do i = 1, IX
              clds(i,5) = 1.0 - cl1(i)          ! save bl cloud
            enddo
          endif
        enddo

        do i = 1, IX
          clds(i,4) = 1.0 - cl1(i)              ! save total cloud
        enddo

      else                                      ! max/ran overlap

        do k = kstr, kend, kinc
          do i = 1, IX
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
            if (ccur >= climit) then             ! cloudy layer
              cl2(i) = min( cl2(i), (1.0 - ccur) )
            else                                ! clear layer
              cl1(i) = cl1(i) * cl2(i)
              cl2(i) = 1.0
            endif
          enddo

          if (k == llyr) then
            do i = 1, IX
              clds(i,5) = 1.0 - cl1(i) * cl2(i) ! save bl cloud
            enddo
          endif
        enddo

        do i = 1, IX
          clds(i,4) = 1.0 - cl1(i) * cl2(i)     ! save total cloud
        enddo

      endif                                     ! end_if_iovr

!  ---  high, mid, low clouds, where cl1, cl2 are cloud fractions
!       layer processed from one layer below llyr and up
!  ---  change! layer processed from surface to top, so low clouds will
!       contains both bl and low clouds.

      if ( ivflip == 0 ) then                   ! input data from toa to sfc

        do i = 1, IX
          cl1 (i) = 0.0
          cl2 (i) = 0.0
          kbt1(i) = NLAY
          kbt2(i) = NLAY
          kth1(i) = 0
          kth2(i) = 0
          idom(i) = 1
          mbot(i,1) = NLAY
          mtop(i,1) = NLAY
        enddo

!org    do k = llyr-1, 1, -1
        do k = NLAY, 1, -1
          do i = 1, IX
            id = idom(i)
            id1= id + 1

            pcur = plyr(i,k)
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))

            if (k > 1) then
              pnxt = plyr(i,k-1)
              cnxt = min( ovcst, max( cldtot(i,k-1), cldcnv(i,k-1) ))
            else
              pnxt = -1.0
              cnxt = 0.0
            endif

            if (pcur < ptop1(i,id1)) then
              id = id + 1
              id1= id1 + 1
              idom(i) = id
            endif

            if (ccur >= climit) then
              if (kth2(i) == 0) kbt2(i) = k
              kth2(i) = kth2(i) + 1

              if ( iovr == 0 ) then
                cl2(i) = cl2(i) + ccur - cl2(i)*ccur
              else
                cl2(i) = max( cl2(i), ccur )
              endif

              if (cnxt < climit .or. pnxt < ptop1(i,id1)) then
                kbt1(i) = nint( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i) )      &
     &                  / (cl1(i) + cl2(i)) )
                kth1(i) = nint( (cl1(i)*kth1(i) + cl2(i)*kth2(i) )      &
     &                  / (cl1(i) + cl2(i)) )
                cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)

                kbt2(i) = k - 1
                kth2(i) = 0
                cl2 (i) = 0.0
              endif   ! end_if_cnxt_or_pnxt
            endif     ! end_if_ccur

            if (pnxt < ptop1(i,id1)) then
              clds(i,id) = cl1(i)
              mtop(i,id) = min( kbt1(i), kbt1(i)-kth1(i)+1 )
              mbot(i,id) = kbt1(i)

              cl1 (i) = 0.0
              kbt1(i) = k - 1
              kth1(i) = 0

              if (id1 <= NK_CLDS) then
                mbot(i,id1) = kbt1(i)
                mtop(i,id1) = kbt1(i)
              endif
            endif     ! end_if_pnxt

          enddo       ! end_do_i_loop
        enddo         ! end_do_k_loop

      else                                      ! input data from sfc to toa

        do i = 1, IX
          cl1 (i) = 0.0
          cl2 (i) = 0.0
          kbt1(i) = 1
          kbt2(i) = 1
          kth1(i) = 0
          kth2(i) = 0
          idom(i) = 1
          mbot(i,1) = 1
          mtop(i,1) = 1
        enddo

!org    do k = llyr+1, NLAY
        do k = 1, NLAY
          do i = 1, IX
            id = idom(i)
            id1= id + 1

            pcur = plyr(i,k)
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))

            if (k < NLAY) then
              pnxt = plyr(i,k+1)
              cnxt = min( ovcst, max( cldtot(i,k+1), cldcnv(i,k+1) ))
            else
              pnxt = -1.0
              cnxt = 0.0
            endif

            if (pcur < ptop1(i,id1)) then
              id = id + 1
              id1= id1 + 1
              idom(i) = id
            endif

            if (ccur >= climit) then
              if (kth2(i) == 0) kbt2(i) = k
              kth2(i) = kth2(i) + 1

              if ( iovr == 0 ) then
                cl2(i) = cl2(i) + ccur - cl2(i)*ccur
              else
                cl2(i) = max( cl2(i), ccur )
              endif

              if (cnxt < climit .or. pnxt < ptop1(i,id1)) then
                kbt1(i) = nint( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i))       &
     &                  / (cl1(i) + cl2(i)) )
                kth1(i) = nint( (cl1(i)*kth1(i) + cl2(i)*kth2(i))       &
     &                  / (cl1(i) + cl2(i)) )
                cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)

                kbt2(i) = k + 1
                kth2(i) = 0
                cl2 (i) = 0.0
              endif     ! end_if_cnxt_or_pnxt
            endif       ! end_if_ccur

            if (pnxt < ptop1(i,id1)) then
              clds(i,id) = cl1(i)
              mtop(i,id) = max( kbt1(i), kbt1(i)+kth1(i)-1 )
              mbot(i,id) = kbt1(i)

              cl1 (i) = 0.0
              kbt1(i) = k + 1
              kth1(i) = 0

              if (id1 <= NK_CLDS) then
                 mbot(i,id1) = kbt1(i)
                 mtop(i,id1) = kbt1(i)
              endif
            endif     ! end_if_pnxt

          enddo       ! end_do_i_loop
        enddo         ! end_do_k_loop

      endif                                     ! end_if_ivflip

!
      return
!...................................
      end subroutine gethml
!-----------------------------------


!-----------------------------------                                    !
      subroutine rhtable                                                &
!...................................                                    !

!  ---  inputs:
     &     ( me                                                         &
!  ---  outputs:
     &,      ier )

!  ===================================================================  !
!                                                                       !
! abstract: cld-rh relations obtained from mitchell-hahn procedure,     !
!   here read cld/rh tuning tables for day 0,1,...,5 and merge into 1   !
!   file.                                                               !
!                                                                       !
! program history log:                                                  !
!   03-xx-1993     kenneth campana     - created original crhtab        !
!   02-xx-1994     kenneth campana     - use only one table for all     !
!                                        forecast hours                 !
!   08-xx-1997     kenneth campana     - smooth out last bunch of       !
!                                        bins of the tables             !
!   04-21-2003     yu-tai hou          - seperate prognostic and        !
!                         diagnostic cloud schemes, re-write into f90   !
!                         modulized form.                               !
!                                                                       !
!                                                                       !
! inputs:                                                               !
!   me              : check print control flag                          !
!                                                                       !
! outputs:                                                              !
!   ier             : error flag                                        !
!                                                                       !
!  ===================================================================  !
!
      implicit none!

!  ---  inputs:
      integer, intent(in) :: me

!  ---  output:
      integer, intent(out) :: ier

!  ---  locals:
      real (kind=kind_phys), dimension(NBIN,NLON,NLAT,MCLD,NSEAL) ::    &
     &      rhfd, rtnfd, rhcf, rtncf, rhcla

      real (kind=kind_io4), dimension(NBIN,NLON,NLAT,MCLD,NSEAL) ::     &
     &      rhfd4, rtnfd4

      real(kind=kind_io4)  :: fhour

      real(kind=kind_phys) :: binscl, cfrac, clsat, rhsat, cstem

      integer, dimension(NLON,NLAT,MCLD,NSEAL) :: kpts, kkpts

      integer :: icdays(15), idate(4), nbdayi, isat

      integer :: i, i1, j, k, l, m, id, im, iy

!
!===> ...  begin here
!

      ier = 1

      rewind NICLTUN

      binscl = 1.0 / NBIN

!  ---  array initializations

      do m=1,NSEAL
       do l=1,MCLD
        do k=1,NLAT
         do j=1,NLON
          do i=1,NBIN
            rhcf (i,j,k,l,m) = 0.0
            rtncf(i,j,k,l,m) = 0.0
            rhcla(i,j,k,l,m) = -0.1
          enddo
         enddo
        enddo
       enddo
      enddo

      kkpts = 0

!  ---  read the data off the rotating file

      read (NICLTUN,ERR=998,END=999) nbdayi, icdays

      if (me == 0) print 11, nbdayi
  11  format('   from rhtable DAYS ON FILE =',i5)

      do i = 1, nbdayi
       id = icdays(i) / 10000
       im = (icdays(i)-id*10000) / 100
       iy = icdays(i)-id*10000-im*100
       if (me == 0) print 51, id,im,iy
  51   format('   from rhtable ARCHV DATA FROM DA,MO,YR=',3i4)
      enddo

      read (NICLTUN,ERR=998,END=999) fhour,idate

      do i1 = 1, nbdayi
        read (NICLTUN) rhfd4
        rhfd = rhfd4

        read (NICLTUN) rtnfd4
        rtnfd = rtnfd4

        read (NICLTUN) kpts

        do m=1,NSEAL
         do l=1,MCLD
          do k=1,NLAT
           do j=1,NLON
            do i=1,NBIN
              rhcf (i,j,k,l,m) = rhcf (i,j,k,l,m) + rhfd (i,j,k,l,m)
              rtncf(i,j,k,l,m) = rtncf(i,j,k,l,m) + rtnfd(i,j,k,l,m)
            enddo
           enddo
          enddo
         enddo
        enddo

        kkpts = kkpts + kpts

      enddo     ! end_do_i1_loop

      do m = 1, NSEAL
       do l = 1, MCLD
        do k = 1, NLAT
         do j = 1, NLON

!  ---  compute the cumulative frequency distribution

           do i = 2, NBIN
             rhcf (i,j,k,l,m) = rhcf (i-1,j,k,l,m) + rhcf (i,j,k,l,m)
             rtncf(i,j,k,l,m) = rtncf(i-1,j,k,l,m) + rtncf(i,j,k,l,m)
           enddo   ! end_do_i_loop

           if (kkpts(j,k,l,m) > 0) then
             do i = 1, NBIN
               rhcf (i,j,k,l,m)= rhcf (i,j,k,l,m)/kkpts(j,k,l,m)
               rtncf(i,j,k,l,m)=min(1., rtncf(i,j,k,l,m)/kkpts(j,k,l,m))
             enddo

!  ---  cause we mix calculations of rh retune with cray and ibm words
!       the last value of rhcf is close to but ne 1.0,
!  ---  so we reset it in order that the 360 loop gives complete tabl

             rhcf(NBIN,j,k,l,m) = 1.0

             do i = 1, NBIN
               lab_do_i1 : do i1 = 1, NBIN
                 if (rhcf(i1,j,k,l,m) >= rtncf(i,j,k,l,m)) then
                   rhcla(i,j,k,l,m) = i1 * binscl
                   exit  lab_do_i1
                 endif
               enddo  lab_do_i1
             enddo

           else                   ! if_kkpts
!  ---  no critical rh

             do i = 1, NBIN
               rhcf (i,j,k,l,m) = -0.1
               rtncf(i,j,k,l,m) = -0.1
             enddo

             if (me == 0) then
               print 210, k,j,m
 210           format('  NO CRIT RH FOR LAT=',I3,' AND LON BAND=',I3,   &
     &                ' LAND(=1) SEA=',I3//'  MODEL RH ',' OBS RTCLD')
               do i = 1, NBIN
                 print 203, rhcf(i,j,k,l,m), rtncf(i,j,k,l,m)
 203             format(2f10.2)
               enddo
             endif

           endif               ! if_kkpts

         enddo    ! end_do_j_loop
        enddo     ! end_do_k_loop
       enddo      ! end_do_l_loop
      enddo       ! end_do_m_loop

      do m = 1, NSEAL
       do l = 1, MCLD
        do k = 1, NLAT
         do j = 1, NLON

           isat = 0
           do i = 1, NBIN-1
             cfrac = binscl * (i - 1)

             if (rhcla(i,j,k,l,m) < 0.0) then
               print 1941, i,m,l,k,j
 1941          format('  NEG RHCLA FOR IT,NSL,NC,LAT,LON=',5I4          &
     &,               '...STOPPP..')
               ! stop
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
             endif

             if (rtncf(i,j,k,l,m) >= 1.0) then
               if (isat <= 0) then
                 isat  = i
                 rhsat = rhcla(i,j,k,l,m)
                 clsat = cfrac
               endif

               rhcla(i,j,k,l,m) = rhsat + (1.0 - rhsat)                 &
     &                         * (cfrac - clsat) / (1.0 - clsat)
             endif
           enddo

           rhcla(NBIN,j,k,l,m) = 1.0

         enddo    ! end_do_j_loop
        enddo     ! end_do_k_loop
       enddo      ! end_do_l_loop
      enddo       ! end_do_m_loop

!  ---  smooth out the table as it reaches rh=1.0, via linear interpolation
!       between location of rh ge .98 and the NBIN bin (where rh=1.0)
!       previously rh=1.0 occurred for many of the latter bins in the
!  ---  table, thereby giving a cloud value of less then 1.0 for rh=1.0

      rhcl = rhcla

      do m = 1, NSEAL
       do l = 1, MCLD
        do k = 1, NLAT
         do j = 1, NLON

           lab_do_i : do i = 1, NBIN - 2
             cfrac = binscl * i

             if (rhcla(i,j,k,l,m) >= 0.98) then
               do i1 = i, NBIN
                 cstem = binscl * i1

                 rhcl(i1,j,k,l,m) = rhcla(i,j,k,l,m)                    &
     &                    + (rhcla(NBIN,j,k,l,m) - rhcla(i,j,k,l,m))    &
     &                    * (cstem - cfrac) / (1.0 - cfrac)
               enddo
               exit  lab_do_i
             endif
           enddo  lab_do_i

         enddo    ! end_do_j_loop
        enddo     ! end_do_k_loop
       enddo      ! end_do_l_loop
      enddo       ! end_do_m_loop

      if (me == 0) then
        print *,'completed rhtable for cloud tuninig tables'
      endif
      return

 998  print 988
 988  format(' from rhtable ERROR READING TABLES')
      ier = -1
      return

 999  print 989
 989  format(' from rhtable E.O.F READING TABLES')
      ier = -1
      return

!...................................
      end subroutine rhtable
!-----------------------------------

!
!.............................................!
      end module module_radiation_clouds_nmmb !
!=============================================!

