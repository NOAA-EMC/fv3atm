      subroutine getcp_idea(im,ix,levs,ntrac,adr,xcp,                   
     &                      thermodyn_id,gen_coord_hybrid)
!
      use tracer_const              ! NENS/gsmphys: see tracer_const_h.f:      module tracer_const
!     USE multi_gases_mod, only: rilist=>ri, cpilist=>cpi
!
      implicit none
      integer, intent(in) :: im     ! number of data points in adr (first dim)
      integer, intent(in) :: ix     ! max data points in adr (first dim)
      integer, intent(in) :: levs   ! number of pressure levels
      integer, intent(in) :: ntrac  ! number of tracer
!
      real, intent(in)    :: adr(ix,levs,ntrac)    ! tracer kg/kg
      real, intent(out)   :: xcp(ix,levs)          !CP (J/kg/k)
      integer thermodyn_id
      logical gen_coord_hybrid
!
! local
!
      real sumq(ix,levs),work1
      integer i,j,k
      integer  :: ntb ! no any other configurations for WAM with NTRAC=5
      sumq = 0.0
      xcp  = 0.0
      if( gen_coord_hybrid .and. thermodyn_id == 3 ) then
        ntb = 1
      elseif (ntrac >= 4) then
        ntb = 4
      else
        return
      endif
!
!VAY-2016
!                             ! reflecting computation for 3 major tracers O-O2-N2 (cpi[0] for N2
      do i=ntb,ntrac
        if( cpi(i) /= 0.0 ) then
          do k=1,levs
            do j=1,im
              work1     = adr(j,k,i)
              sumq(j,k) = sumq(j,k) + work1
              xcp(j,k)  = xcp(j,k)  + work1*cpi(i)
            enddo
           enddo
        endif
      enddo
      do k=1,levs
        do j=1,im
          xcp(j,k) = xcp(j,k) + (1.-sumq(j,k))*cpi(0)
        enddo
      enddo
      return
      end subroutine getcp_idea
!
!
      subroutine rad_merge(im,ix,levs,hlw,swh,prsi,prsl,wtot,           
     &                 xmu,dtrad, dtco2c,dtco2h,dth2oh,dth2oc,dto3)
!
! VAY-201702
!
      use idea_solar,only :  xb, xt, rdx, xlogps
      implicit none
      integer, intent(in) :: im               ! number of data points in hlw,dt..(first dim)
      integer, intent(in) :: ix               ! max data points in hlw,... (first dim)
      integer, intent(in) :: levs             ! number of pressure levels
!      real, parameter     :: xb=7.5, xt=8.5   ! for Hp = 7 km:  52.5 km < Z_logp < 59.5 km 
!      real, parameter     :: rdx=1./(xt-xb)
!      real, parameter     :: xlogps = 11.5129 ! alog(1.e5=Ps_in_Pa)
      real, intent(in)    :: hlw(ix,levs)     ! GFS lw rad (K/s)
      real, intent(in)    :: swh(ix,levs)     ! GFS sw rad (K/s)
      real, intent(in)    :: prsi(ix,levs+1)  ! pressure
      real, intent(in)    :: prsl(ix,levs)    ! pressure
!
      real, intent(in)    :: xmu(im)          ! im-1D array GSM-rules  normalized(??) cos zenith angle
!
      real, intent(in)    :: dtrad (ix,levs)  ! idea EUV-UV(SRB-SRC-Lya) MLT (K/s)
      real, intent(in)    :: dtco2c(ix,levs)  ! idea co2 cooling(K/s)
      real, intent(in)    :: dtco2h(ix,levs)  ! idea co2 heating(K/s)
      real, intent(in)    :: dth2oc(ix,levs)  ! idea h2o cooling(K/s)
      real, intent(in)    :: dth2oh(ix,levs)  ! idea h2o heating(K/s)
      real, intent(in)    :: dto3(ix,levs)    ! idea o3 heating(K/s)
      real, intent(out)   :: wtot(ix,levs)    ! GFS idea combined  rad
     
!     local
      real xk,wl,wh
      integer i,k,j
!
      do k=1,levs
        do i=1,im
!
!vay-2016          xk = log(prsi(i,1)/prsl(i,k))
!   fixing incorrect IDEA-2013/14 dep-ce on the surface pressure in the "Isobaric" Stratosphere
!
          xk  = xlogps - alog(prsl(i,k))
          wh = dtco2c(i,k)+dth2oc(i,k)+dtco2h(i,k)+dth2oh(i,k)+dto3(i,k)+dtrad(i,k)
          wl = hlw(i,k)+swh(i,k)*xmu(i)    ! see analog in gbphys.f: dt3dt(i,k) = dt3dt(i,k) + swh(i,k)*dtf*xmu(i)
          if(xk < xb) then
             wtot(i,k) = wl
          elseif(xk >= xb .and. xk <= xt) then
             wtot(i,k) = (wl*(xt-xk) + wh*(xk-xb))*rdx     ! 52.5 km <= Z_logp <= 59.5 km
          else
             wtot(i,k) = wh
          endif
        enddo
      enddo
      return
      end subroutine rad_merge
!
!======================= Hydrostatic equation: Zgeo and Gravity(Zgeo)
!
      subroutine phi2z(im,ix,levs,phi,soro,z,grav)

! Subroutine to calculate geometric height and gravity from geopotential
! in a hydrostatic atmosphere, assuming a spherically symmetric planet
! and Newton's gravity.

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! File history

! Feb 26, 2010: Rashid Akmaev
! Loosely based on Hojun Wang's phi2hgt but generalized to rid of
! recursive calculations, include surface orography, and calculate 
! gravity.

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Define constants
! - Earth radius (m) and 
! - gravity at sea level (m/s**2) 

! If used with GFS/WAM codes "use" this module
      use physcons, only: re => con_rerth, g0 => con_g

      implicit none

! If the module is not available, comment out the "use" line above and
! uncomment this line
!      real, parameter:: re = 6.3712e+6, g0 = 9.80665e+0

      real, parameter:: g0re = g0*re, g0re2 = g0*re*re

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Subroutine parameters
! INPUTS
! - array dimensions (following GFS conventios): first actual, first 
! maximum, number of levels

      integer, intent(in):: im,ix,levs

! - geopotential (m**2/s**2)
! - surface orography (m)

      real, intent(in):: phi(ix,levs)
      real, intent(in):: soro(im)

! OUTPUTS
! - height (m)
      real, intent(out):: z(ix,levs)
! - gravity (m/s**2)
      real, intent(out):: grav(ix,levs)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Local variables

      integer:: i,l
      real:: phis(im)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Calculate surface geopotential

      do i = 1,im
        phis(i) = g0re*soro(i)/(re+soro(i))
      enddo

! Calculate height

      do l = 1,levs
        do i = 1,im
          z(i,l) = re*(phis(i)+phi(i,l))/(g0re-(phis(i)+phi(i,l)))
        enddo
      enddo

! calculate gravity

         do l = 1,levs
           do i = 1,im
             grav(i,l) = g0re2/((re+z(i,l))*(re+z(i,l)))
           enddo
         enddo

      end subroutine phi2z
!----------------------------------------------------------------------------
      subroutine gravco2(levs,phi,soro,gg)
!----------------------------------------------------------------------------
! VAY OCT-2016:   gloopb.f: call gravco2(levs,philco2,sfc_fld%oro(1,lanlat1),gg1
! TODO this should be revized because of "the first gaussian point" =/= global 
!                better to put soro =0.
!
! Subroutine is "lously" modified from phi2z above to compute gravity for co2cin,
! the first gaussian point is chosen to represent the whole data domain

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! File history

! Dec 26, 2012: Jun Wang        modified from phi2z from Rashid

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Define constants
! - Earth radius (m) and
! - gravity at sea level (m/s**2)

! If used with GFS/WAM codes "use" this module
      use physcons, only: re => con_rerth, g0 => con_g

      implicit none

! If the module is not available, comment out the "use" line above and
      real, parameter:: g0re = g0*re, g0re2 = g0*re**2

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Subroutine parameters
! INPUTS

      integer, intent(in):: levs

! - geopotential (m**2/s**2)
! - surface orography (m)

      real, intent(in) :: phi(levs)
      real, intent(in) :: soro

! OUTPUTS   gravity (m/s**2)

      real, intent(out):: gg(levs)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Local variables

      integer:: i,l
      real:: phis
      real:: z(levs)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Calculate surface geopotential

      phis = g0re*soro/(re+soro)
!      print *,'in grevco2 phis=',phis,'phi=',phi(1:100:10),'soro=',soro,
!     &   're=',re

! Calculate height

      do l = 1,levs
        z(l) = re*(phis+phi(l))/(g0re-(phis+phi(l)))
      enddo

! calculate gravity

      do l = 1,levs
         gg(l) = g0re2/((re+z(l))*(re+z(l)))
      enddo
!      print *,'in grevco2 gg=',gg(1:100:10)
!
      end subroutine gravco2
!----------------------------------------------------------------------------
      subroutine getphilvl(levs,ntrac,ps,t,q,dp,gen_coord_hybrid,
     &  thermodyn_id,phil,prsi)
!
!..same thing as with gravco2: it copmputes PHI in gloopb.f for init-n of WAM 
!
!gloopb.f:              call getphilvl
!
!              call getphilvl(levs,ntrac, grid_fld%ps(1,lanlat1),
!     &                   grid_fld%t(1,lanlat1,1:levs),qtrac,
!     &                   grid_fld%dp(1,lanlat1,1:levs),gen_coord_hybrid,
!     &                   thermodyn_id,philco2,prsilvl1)
!                                       phil,    prsi)
! change prsi from cb to pascal
!              prsilvl1 = prsilvl1*1000. in [Pa]
!              call gravco2(levs,philco2,sfc_fld%oro(1,lanlat1),gg1)
!              call ideaca_init(prsilvl,levs+1)
! vay 05/2015  it copmputes PHI in gloopb.f 
!                                 prsi in [Pa]  not in cpbs as in Tprsi 
! Subroutine computes phi on a single point on model levels from p,tmp,
!  and tracers for general
! 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! File history

! Dec 26, 2010: Jun Wang

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      use tracer_const, only : ri
!     USE multi_gases_mod, only: rilist=>ri
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Subroutine parameters
! INPUTS

      integer, intent(in):: levs,ntrac
      logical, intent(in):: gen_coord_hybrid
      integer, intent(in):: thermodyn_id
!
! Local variables
      real,parameter :: pa2cb=0.001, zero=0.0
! - sfc pressure  (pascal)
! - pressure  thickness (pascal)
! - tmp (k)
! - tracers

      real, intent(in):: ps,t(levs),dp(levs),q(levs,ntrac)

! OUTPUTS
! 

      real, intent(out):: phil(levs),prsi(levs+1)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Local variables
      real:: tem,dphi,phii,sumq(levs),xr(levs)
      integer :: k,n

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! init



! Calculate enthalpy
!
!       print *,'in getphilvl,thermodyn_id=',thermodyn_id,
!     &    thermodyn_id.eq.3

      if( gen_coord_hybrid ) then
        if( thermodyn_id == 3 ) then           ! Enthalpy case
!get r
          sumq = zero
          xr   = zero
          do n=1,ntrac
            if( ri(n) > 0.0 ) then
              do k=1,levs
                xr(k)   = xr(k)   + q(k,n) * ri(n)
                sumq(k) = sumq(k) + q(k,n)
              enddo
            endif
          enddo
!
! add N2 with ri(0)
!
          do k=1,levs
            xr(k)    = (1.-sumq(k))*ri(0)  + xr(k)
          enddo
!
! pressure interfaces, options get "dP*0.001"
!
          prsi(levs+1) = 0.
          do k=levs,1,-1
            prsi(k) = prsi(k+1) + dp(k)*pa2cb
          enddo
!          print *,'in getphilvl,prsi=',prsi(1:100:10)
!
          phii = zero               ! should be surface Phis assumes sea level
                                    ! no consistency with gravco2
          do k = 1,levs
            tem   = xr(k) * T(k)
            dphi  = tem*(prsi(k) - prsi(k+1))/(prsi(k) + prsi(k+1)) ! RT*dln(Pi)/Pmid
            phil(k)   = phii + dphi       ! dphi =1/2*DPHI for mid-points
!
!    OK                         but dphi = 1/2 of DPHI between interfaces              
            phii      = phil(k) + dphi    ! Phii(k-1)+ 2*dphi
          enddo
!
        else
!          print *,'ERROR: No phil is compute, this routine is ',
!     &          'for gen-hybrid  with enthalpy'
!
        print *, '  thermodyn_id  in getphilvl  ', thermodyn_id 
        endif
      endif
!
      end subroutine getphilvl
!------------------------------------------------------------------
      subroutine getmax(ain,n1,n,m,rmin,j1,rmax,j2)
      real ain(n1,m)
      rmin =  1.e36
      rmax = -1.e36
      i1 = 500
      j1 = 500
      i2 = 500
      j2 = 500
      do j=1,m
        do i=1,n
          if(rmin > ain(i,j)) then
            rmin = ain(i,j)
            i1 = i
            j1 = j
          endif
          if(rmax < ain(i,j)) then
            rmax = ain(i,j)
            i2 = i
            j2 = j
          endif
        enddo
      enddo
      return
      end subroutine getmax
!
      subroutine getmax2(ain,ain1,n1,n,m,rmax,j2)
      real ain(n1,m),ain1(n1,m)
      rmax = -1.e36
      i1   = 500
      j1   = 500
      i2   = 500
      j2   = 500
      do j=1,m
        do i=1,n
          sq = sqrt(ain(i,j)*ain(i,j) + ain1(i,j)*ain1(i,j))
          if(rmax < sq) then
            rmax = sq
            i2   = i
            j2   = j
          endif
        enddo
      enddo
      RETURN
      END subroutine getmax2
!
!vay Oct/2015
!      CALL get_exner_phi_zgrav(ix, im, levs, ntrac,
!     & adr, adt,
!     & prsl,  prsi, rdelp, prsik, prslk,phii,phil, 
!     & thermodyn_id, sfcpress_id, gen_coord_hybrid,                                    
!     & oro,zg, grav, exner, delpa) 
!
       SUBROUTINE  get_exner_phi_zgrav(ix, im, levs, ntrac,
     & adr, adt, 
     & prsl, prsi, rdelp, prsik, prslk, phii,phil, 
     & thermodyn_id, sfcpress_id, gen_coord_hybrid,                                    
     & oro,zg, grav, exner, delPa)
!
      IMPLICIT NONE
      integer, intent(IN) :: ix, im, levs, ntrac    ! 
      real, intent(in)    :: prsl(ix, levs)         ! Pa
      real, intent(in)    :: prsi(ix, levs+1)       ! Pa
      real, intent(in)    :: adT(ix, levs)
      real, intent(in)    :: adR(ix, levs, ntrac)

      real, intent(in)    :: oro(iM)
      integer,intent(in)  :: thermodyn_id, sfcpress_id
      logical,intent(in)  :: gen_coord_hybrid
! output
!
      real, intent(out)    :: rdelp(ix, levs), delpa(ix, levs)
      real, intent(out)    :: phil(ix, levs)
      real, intent(out)    :: prslk(ix, levs)
      real, intent(out)    :: zg(ix, levs)
      real, intent(out)    :: grav(ix, levs)
      real, intent(out)    :: exner(ix, levs)
!

      real, intent(out)    :: prsiK(ix, levs+1)
      real, intent(out)    :: phii(ix, levs+1)
!
!locals
      real, parameter     ::  p00i=1.0e-5
!
      real :: xr(ix, levs), xcp(ix,levs)               ! values before tracer updates 
      real :: kappa(ix, levs), tem
      integer :: i,k 
! 
!  COMPUTE phil, phii, zg, grav, prsik,prslk, exner
!
      call GET_CPR(im,ix,levs,ntrac,adr,xcp,xr)        ! see get_prs.f
!
      call GET_PHI_WAM_PinPa(im,ix,levs, adt,                                                              
     &             prsi, prsl, xr, phii,phil)
!
! get height in meters "zg"
!
      call phi2z(im,ix,levs,phil,oro,zg,grav)
!
! get exner prsik / prslk
!
            do k=1,levs
              do i=1,im         
              DelPa(i,k) = prsi(i,k) -prsi(i,k+1) 
              rdelp(i,k) =1./DelPa(i,k)
                kappa(i,k) = xr(i,k)/xcp(i,k)
                prslk(i,k)  = (prsl(i,k)*p00i) ** kappa(i,k)  ! T = prkl*PT
                exner(i,k) =1./prslk(i,k)                     ! PT= exner*T
              enddo
            enddo 
            do k=2,levs
              do i=1,im
               tem = 0.5 * (kappa(i,k) + kappa(i,k-1))
               prsik(i,k-1) = (prsi(i,k)*p00i) ** tem 
              enddo
            enddo
            do i=1,im
              prsik(i,1) = (prsi(i,1)*p00i) ** kappa(i,1)
              prsik(i,levs+1) = prsik(i, levs)                 ! make non-zero value at the top for PT<=> T transforms
            enddo
            k = levs + 1

            if (prsi(1,k) .gt. 0.0) then
              do i=1,im
                prsik(i,k) = (prsi(i,k)*p00i) ** kappa(i,levs) ! make "realistic" value at the top for PT<=> T transforms
              enddo
            endif
! 
      END SUBROUTINE  get_exner_phi_zgrav
!
      subroutine GET_PHI_WAM_PinPa(im,ix,levs,t,
     &                   prsi, prsl, XR, phii,phil)
!======================================================
!   purpose: get phi using preesure in Pascals
!======================================================
      USE MACHINE ,              ONLY : kind_phys
!
      implicit none
!
      integer, intent(in) ::  im, ix, levs

!
      real(kind=kind_phys), intent(in) :: prsi(ix,levs+1)     ! pressure in SI-units Pascal !!!!!
      real(kind=kind_phys), intent(in) :: prsl(ix,levs)
      real(kind=kind_phys), intent(in) :: T(ix,levs)          ! temp-re in Kelvins
      real(kind=kind_phys), intent(in) :: XR(ix,levs)    
!
      real(kind=kind_phys),intent(out) :: phii(ix,levs+1)
      real(kind=kind_phys),intent(out) :: phil(ix,levs)        ! g*Z = m2/s2
    
      real(kind=kind_phys) tem, dphi
      real (kind=kind_phys), parameter :: zero=0.0
      integer i, k, n
!
      do i=1,im
        phii(i,1)   = zero                         ! Ignoring topography height here
      enddo
! 
          DO k=1,levs
            do i=1,im
              TEM         = xr(i,k) * T(i,k)   !
              DPHI        = (PRSI(i,k) - PRSI(i,k+1)) * TEM
     &                     /(PRSI(i,k) + PRSI(i,k+1))
              phil(i,k)   = phii(i,k) + DPHI
              phii(i,k+1) = phil(i,k) + DPHI
!
            ENDDO
          ENDDO
!

      return
      end subroutine GET_PHI_WAM_PinPa
