! Apr 06 2012    Henry Juang, initial implement for nems
! Oct 20 2015    Weiyu Yang,  add f10.7 inputted data.
! Oct 15 2016    VAY, put f10.7 as a parameter
! Nov    2016    Correction of JO2 and oxygen chemistry
! September 2017 Weiyu Yang, add IPE back coipling WAM code.
! October 2017   Rashid Akmaev, eddy mixing, corrections and clean-up
!-----------------------------------------------------------------------

      module idea_tracer_mod
! hold jprofile-old; oh-ho2 global profiles

      implicit none
      real, allocatable::  jj(:)
      real, allocatable::  oh(:), ho2(:)
      real, allocatable:: vmr_glob(:,:)     !(levs, nvmr)
      end module idea_tracer_mod
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine idea_tracer_init(levs)

      use idea_tracer_mod,    only :   jj, oh, ho2, vmr_glob
      use idea_tracers_input, only :
     &    jprofile, hprofile, init_tracer_constants,  WAM_GLOBAL_TRACERS
      use IDEA_MPI_def,       only : mpi_WAM_quit

      implicit none
      real, parameter :: f107=100.  ! any value of f107 to init 1D-JJs
                                    ! now updated with time...JO2-3D
      integer, parameter :: nvmr=15 !  15-global vertical arrays returned by WAM_GLOBAL_TRACERS
      integer, intent(in):: levs    !number of pres levels

      allocate (jj(levs))
      allocate (oh(levs), ho2(levs))
      allocate (vmr_glob(levs, nvmr))

      call jprofile(levs,f107,jj)    ! old 1D-Jo2 profile
      call hprofile(levs,oh,ho2)     ! 1D [oh-ho2]  profiles 
      call init_tracer_constants     

      call WAM_GLOBAL_TRACERS(levs, nvmr, vmr_glob)
      print *, ' VAY WAM_GLOBAL_TRACERS INIT'
      return
      end subroutine idea_tracer_init
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine idea_tracer(im,ix,levs,ntrac,ntrac_i,grav,prsi,prsl,   
     &  adt,q,dtp,n1,n2,ozn, n3,n,rho,am, am29, 
     &  cospass,dayno,zg,f107,f107d,me, go2dr, plow, phigh, xpk_low,
     &  xpk_high)

      use physcons, only          : avgd => con_avgd             
      use idea_composition, only  : bz,amo,amn2, amo2, amo3, amh2o
      use idea_composition, only  : rbz, rmo, rmo2, rmn2, rmh2o,rmo3
      use idea_tracer_mod,   only : jj

      use module_IPE_to_WAM, only: lowst_ipe_level, ipe_to_wam_coupling

      implicit none
! Argument
      integer, intent(in) :: me              ! current PE
      integer, intent(in) :: im              ! number of long data points
      integer, intent(in) :: ix              ! max data points in fields
      integer, intent(in) :: levs            ! number of pressure levels
      integer, intent(in) :: ntrac           ! number of tracer (total)
      integer, intent(in) :: ntrac_i         ! number of tracer add by IDEA
      integer, intent(in) :: dayno           ! calendar day
      real, intent(in)    :: cospass(im)     ! cos zenith angle
      real, intent(in)    :: zg(ix,levs)     !layer height (m)
      real, intent(in)    :: f107, f107d     ! variable F107 to recomput JJ
      real, intent(in)    :: prsi(ix,levs+1) ! interface pressure in Pa
      real, intent(in)    :: prsl(ix,levs)   ! layer pressure in Pa
      real, intent(in)    :: grav(ix,levs)   ! (m/s2)
      real, intent(in)    :: adt(ix,levs)    ! input  temp 
      real, intent(in)    :: dtp             ! time step in second
      real, intent(in), dimension(ix)    :: xpk_low, xpk_high
      integer,intent(in), dimension(ix) :: plow, phigh
      real, intent(in)    :: go2dr(ix, lowst_ipe_level:levs)
      real, intent(inout) :: q(ix,levs,ntrac)   ! input output tracer
      real, intent(out)   :: n1(ix,levs)     ! number density of o (/m3)
      real, intent(out)   :: n2(ix,levs)     ! number density of o2 (/m3)
      real, intent(out)   :: ozn(ix,levs)    ! number density of o2 (/m3)
      real, intent(out)   :: n3(ix,levs)     ! number density of n2 (/m3)
      real, intent(out)   :: n(ix,levs)      ! total number density (/m3)
      real, intent(out)   :: rho(ix,levs)    ! density of  (kg/m3)
      real, intent(out)   :: am(ix,levs)     ! avg mass of mix  (kg)
      real, intent(out)   :: am29(ix,levs)   ! avg mass of air 28.84 => 16.
! local argument
      real ::  dq1(ix,levs,ntrac_i),dq2(ix,levs,ntrac_i),
     &         dq3(ix,levs,ntrac_i)
      real ::  mh2o,mo3, mo,mo2,mn2        
      real ::  qin(ix,levs,ntrac_i)
      real ::  qn2
      integer i,k,in
       real :: Jrates_O2(ix, levs)            ! O2 dissociation rate

!     Two tracers added by IDEA, O and O2
      do in=1,ntrac_i
        do i=1,im
          do k=1,levs
            qin(i,k,in)=q(i,k,ntrac-ntrac_i+in)
          enddo
        enddo
      enddo

! mean mass, mass and number densities
!     here n,n1,n2 in /m3 , rho in kg/m3
      do i=1,im
        do k=1,levs
           qn2=1.-q(i,k,1)-q(i,k,2)-q(i,k,4)-q(i,k,5)

! mean molecular mass of gaseous tracers
           am(i,k)=1./(qin(i,k,1)*rmo+qin(i,k,2)*rmo2+q(i,k,1)*rmh2o+       
     &          q(i,k,2)*rmo3+qn2*rmn2)

! total number density and mass density
          n(i,k)=rbz*prsl(i,k)/adt(i,k)
          rho(i,k)=am(i,k)*n(i,k)

! partial number densities for radiation and chemistry,
!     make sure non-negative
          n1(i,k)=max(qin(i,k,1)*rho(i,k)*rmo,0.)
          n2(i,k)=max(qin(i,k,2)*rho(i,k)*rmo2,0.)
          n3(i,k) =max(qn2*rho(i,k)*rmn2,0.)
          ozn(i,k)= max(q(i,k,2)*rho(i,k)*rmo3,0.)
        enddo
      enddo

! Eddy mixing
      call idea_tracer_eddy(im,ix,levs,ntrac_i,grav,prsi,prsl,rho,dtp,
     &     qin,dayno,dq3)

! Mutual molecular diffusion of major thermospheric species O, O2, N2
      call idea_tracer_m(im,ix,levs,ntrac_i,grav,prsi,prsl,adt,dtp,     
     &     qin,am,dq1)

! O2 dissociation rate
      call idea_dissociation_jo2(im,ix,levs,adt,cospass,n1,n2, ozn, n3,
     &          dayno,zg,grav, f107, f107d, Jrates_O2)

! Merge the  IPE back coupling WAM variable arrays into WAM.
      IF (ipe_to_wam_coupling) THEN
        call idea_merge_ipe_to_wam(GO2DR, jrates_O2,
     &    im, ix, levs, lowst_ipe_level, prsl, plow, phigh,
     &    xpk_low, xpk_high)
      END IF

! Oxygen chemistry + ad hoc HOx sinks of O
      call idea_tracer_c(im,ix,levs,ntrac_i,adt,dtp,Jrates_O2,
     &     n1,n2,n,rho, qin,dq2)

! Update mmr of O and O2 due to mixing, diffusion, and chemistry
      do in=1,ntrac_i
        do i=1,im
          do k=1,levs
            q(i,k,in+ntrac-ntrac_i)=q(i,k,in+ntrac-ntrac_i)+            
     &        dq1(i,k,in)  +dq2(i,k,in) + dq3(i,k,in)
          enddo
        enddo
      enddo

! Update number densities (mass density is conserved) and mean mass
      do i=1,im
         do k=1,levs
            qn2=1.-q(i,k,1)-q(i,k,2)-q(i,k,4)-q(i,k,5)
            am(i,k)=1./(qin(i,k,1)*rmo+qin(i,k,2)*rmo2+q(i,k,1)*rmh2o+       
     &           q(i,k,2)*rmo3+qn2*rmn2)
            n1(i,k)=max(qin(i,k,1)*rho(i,k)*rmo,0.)
            n2(i,k)=max(qin(i,k,2)*rho(i,k)*rmo2,0.)
            n3(i,k) =max(qn2*rho(i,k)*rmn2,0.)
            ozn(i,k)= max(q(i,k,2)*rho(i,k)*rmo3,0.)
         enddo
      enddo
         am29 = am * 1.e3 * avgd

      return
      end

      subroutine idea_tracer_m(im,ix,levs,ntrac_i,grav,prsi,prsl,adt,   
     &dtp,qin,am,dq)

! Calculate tracer changes by mutual molecular diffusion of
!     major thermospheric species O, O2, N2
! 2006 Rashid Akmaev and Fei Wu
!-----------------------------------------------------------------------
      use physcons,  only :rgas=>con_rgas,            
     &               avgd => con_avgd
      use machine, only : kind_phys
      use idea_composition, only:  amo, amo2, amn2, bz, rbz
      implicit none
! Argument
      integer, intent(in) :: im    ! number of long data points in fields
      integer, intent(in) :: ix    ! max data points in fields
      integer, intent(in) :: levs  ! number of pressure levels
      integer, intent(in) :: ntrac_i ! number of tracer add by IDEA
      real,    intent(in) :: dtp   ! time step in second
      real, intent(in)    :: prsi(ix,levs+1) ! interface pressure in KPa
      real, intent(in)    :: prsl(ix,levs)   ! layer pressure in KPa
      real, intent(in)    :: grav(ix,levs)   ! (m/s2)
      real, intent(in) :: adt(ix,levs)   ! input  temp at dt=0
      real, intent(in) :: qin(ix,levs,ntrac_i)   ! input tracer
      real, intent(in)   :: am(ix,levs)   ! avg mass of mix  (kg)
      real, intent(out):: dq(ix,levs,ntrac_i) ! output tracer changes
!local  variables
      real n1_i(levs+1),n2_i(levs+1),n3_i(levs+1),n_i(levs+1)
      real t_i(levs+1),am_i(levs+1),qout(ix,levs,ntrac_i)
      real beta(2,2,levs+1),a(2,2,levs),b(2,2,levs),c(2,2,levs)
      real ggg(2,2),ee(2,2,levs+1),f(2,levs+1),                         
     &     d12,d13,d23,a12,a13,a23,s12,s13,s23,mo,mo2,mn2,              
     &     dp1(levs),dp1_i(levs+1)
      real partb_i(levs+1),parta(levs),hold1,dtp1,hold2
      integer k,i,kk,kk1,in
! change unit from g/mol to kg
      mo=amo*1.e-3/avgd
      mo2=amo2*1.e-3/avgd
      mn2=amn2*1.e-3/avgd
! some constants
      a12=9.69e18
      a13=9.69e18
      a23=8.3e18

      s12=0.774
      s13=0.774
      s23=0.724
! set boundary
      beta(1:2,1:2,1)=0.
      beta(1:2,1:2,levs+1)=0.
      a(1:2,1:2,1)=0.
      c(1:2,1:2,levs)=0.
      ee(1:2,1:2,levs+1)=0.
      f(1:2,levs+1)=0.
!
      dtp1=1./dtp
      t_i=0.
      am_i=0.
      n_i=0.
      n1_i=0.
      n2_i=0.
      n3_i=0.
!
! for each longitude
!
      do i=1,im
! calculate temp in interface pressure levels
! get compositions at interface pressure levels
        do k=2,levs
          t_i(k)=(adt(i,k-1)+adt(i,k))*.5
          am_i(k)=.5*(am(i,k-1)+am(i,k))
          n_i(k)=prsi(i,k)/bz/t_i(k) 
          n1_i(k)=max(0.,
     &         .5*(qin(i,k,1)+qin(i,k-1,1))*am_i(k)*n_i(k)/mo)
          n2_i(k)=max(0.,
     &         .5*(qin(i,k,2)+qin(i,k-1,2))*am_i(k)*n_i(k)/mo2)
          n3_i(k)=max(0.,n_i(k)-n1_i(k)-n2_i(k))
        enddo

! calculate beta at interface pressure
        do k=2,levs
          d12=a12*t_i(k)**(s12)
          d13=a13*t_i(k)**(s13)
          d23=a23*t_i(k)**(s23)
          hold1=1./(n1_i(k)*d23+n2_i(k)*d13+n3_i(k)*d12)
          beta(1,1,k)=hold1*d13*mo*(n1_i(k)*mn2*d23+                    
     &            (n2_i(k)*mo2+n3_i(k)*mn2)*d12)
          beta(2,2,k)=hold1*d23*mo2*(n2_i(k)*mn2*d13+                   
     &            (n1_i(k)*mo+n3_i(k)*mn2)*d12)
          beta(1,2,k)=hold1*d23*mo*n1_i(k)*(mn2*d13-mo2*d12)
          beta(2,1,k)=hold1*d13*mo2*n2_i(k)*(mn2*d23-mo*d12)
        enddo

! solve tridiagonal problem
        do k=1,levs
          dp1(k)=1./(prsi(i,k)-prsi(i,k+1))
          parta(k)=dtp*grav(i,k)*dp1(k)/bz
        enddo
        do k=2,levs
          dp1_i(k)=1./(prsl(i,k-1)-prsl(i,k))
          partb_i(k)=.5*(grav(i,k)+grav(i,k-1))/t_i(k)
        enddo
        do k=2,levs
          hold1=parta(k)*partb_i(k)
          hold2=am(i,k-1)*prsl(i,k-1)*dp1_i(k)
          a(1,1,k)=hold1*beta(1,1,k)*(hold2/mo-.5)
          a(1,2,k)=hold1*beta(1,2,k)*(hold2/mo2-.5)
          a(2,1,k)=hold1*beta(2,1,k)*(hold2/mo-.5)
          a(2,2,k)=hold1*beta(2,2,k)*(hold2/mo2-.5)
         enddo
        do k=1,levs-1
          hold1=parta(k)*partb_i(k+1)
          hold2=am(i,k+1)*prsl(i,k+1)*dp1_i(k+1)
          c(1,1,k)=hold1*beta(1,1,k+1)*(hold2/mo+.5)
          c(1,2,k)=hold1*beta(1,2,k+1)*(hold2/mo2+.5)
          c(2,1,k)=hold1*beta(2,1,k+1)*(hold2/mo+.5)
          c(2,2,k)=hold1*beta(2,2,k+1)*(hold2/mo2+.5)
         enddo
        do k=2,levs-1
          hold1=am(i,k)*prsl(i,k)*dp1_i(k+1)
          hold2=am(i,k)*prsl(i,k)*dp1_i(k)
      b(1,1,k)=1.+parta(k)*(partb_i(k+1)*beta(1,1,k+1)*(hold1/mo-.5)    
     &                    +partb_i(k)*beta(1,1,k)*(hold2/mo+.5))
      b(2,2,k)=1.+parta(k)*(partb_i(k+1)*beta(2,2,k+1)*(hold1/mo2-.5)   
     &                    +partb_i(k)*beta(2,2,k)*(hold2/mo2+.5))
      b(1,2,k)=parta(k)*(partb_i(k+1)*beta(1,2,k+1)*(hold1/mo2-.5)      
     &                    +partb_i(k)*beta(1,2,k)*(hold2/mo2+.5))
      b(2,1,k)=parta(k)*(partb_i(k+1)*beta(2,1,k+1)*(hold1/mo-.5)       
     &                    +partb_i(k)*beta(2,1,k)*(hold2/mo+.5))
        enddo
          hold1=am(i,1)*prsl(i,1)*dp1_i(2)
      b(1,1,1)=1.+parta(1)*partb_i(2)*beta(1,1,2)*(hold1/mo-.5)
      b(2,2,1)=1.+parta(1)*partb_i(2)*beta(2,2,2)*(hold1/mo2-.5)
      b(1,2,1)=parta(1)*partb_i(2)*beta(1,2,2)*(hold1/mo2-.5)
      b(2,1,1)=parta(1)*partb_i(2)*beta(2,1,2)*(hold1/mo-.5)
          hold2=am(i,levs)*prsl(i,levs)*dp1_i(levs)
      b(1,1,levs)=1.+parta(levs)*partb_i(levs)*beta(1,1,levs)*          
     &(hold2/mo+.5)
      b(2,2,levs)=1.+parta(levs)*partb_i(levs)*beta(2,2,levs)*          
     &(hold2/mo2+.5)
      b(1,2,levs)=parta(levs)*partb_i(levs)*beta(1,2,levs)*             
     &(hold2/mo2+.5)
      b(2,1,levs)=parta(levs)*partb_i(levs)*beta(2,1,levs)*             
     &(hold2/mo+.5)
       do k=levs,1,-1
         ggg(1,1)=b(2,2,k)-c(2,1,k)*ee(1,2,k+1)-c(2,2,k)*ee(2,2,k+1)
         ggg(2,2)=b(1,1,k)-c(1,1,k)*ee(1,1,k+1)-c(1,2,k)*ee(2,1,k+1)
         ggg(1,2)=-1.*b(1,2,k)+c(1,1,k)*ee(1,2,k+1)+c(1,2,k)*ee(2,2,k+1)
         ggg(2,1)=-1.*b(2,1,k)+c(2,1,k)*ee(1,1,k+1)+c(2,2,k)*ee(2,1,k+1)
         hold1=1./(ggg(1,1)*ggg(2,2)-ggg(1,2)*ggg(2,1))
         ggg=ggg*hold1
         ee(1,1,k)=ggg(1,1)*a(1,1,k)+ggg(1,2)*a(2,1,k)       
         ee(1,2,k)=ggg(1,1)*a(1,2,k)+ggg(1,2)*a(2,2,k)       
         ee(2,1,k)=ggg(2,1)*a(1,1,k)+ggg(2,2)*a(2,1,k)       
         ee(2,2,k)=ggg(2,1)*a(1,2,k)+ggg(2,2)*a(2,2,k)       
      f(1,k)=ggg(1,1)*(qin(i,k,1)+c(1,1,k)*f(1,k+1)                      
     &+c(1,2,k)*f(2,k+1))+ggg(1,2)*(qin(i,k,2)+c(2,1,k)*f(1,k+1)         
     &+c(2,2,k)*f(2,k+1))
      f(2,k)=ggg(2,1)*(qin(i,k,1)+c(1,1,k)*f(1,k+1)                      
     &+c(1,2,k)*f(2,k+1))+ggg(2,2)*(qin(i,k,2)+c(2,1,k)*f(1,k+1)         
     &+c(2,2,k)*f(2,k+1))
        enddo
        do in=1,ntrac_i
          qout(i,1,in)=f(in,1)
          dq(i,1,in)=qout(i,1,in)-qin(i,1,in)
        enddo
        do k=2,levs
          qout(i,k,1)=ee(1,1,k)*qout(i,k-1,1)+ee(1,2,k)*qout(i,k-1,2)+  
     &              f(1,k)
          qout(i,k,2)=ee(2,1,k)*qout(i,k-1,1)+ee(2,2,k)*qout(i,k-1,2)+  
     &              f(2,k)
          do in=1,ntrac_i
            dq(i,k,in)=qout(i,k,in)-qin(i,k,in)
          enddo
        enddo
      enddo !i
      return
      end subroutine idea_tracer_m

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!   
      subroutine idea_tracer_c(im,ix,levs,ntrac_i,adt,dtp,jo2,n1,n2,     
     &        n,rho,qin,dq)
!-----------------------------------------------------------------------
! calculate tracer changes caused by chemistry reaction
! 2006 Rashid Akmaev and Fei Wu
!      Tim Fuller-Rowell HOx oxygen sinks
!-----------------------------------------------------------------------
      use physcons, only : rgas=>con_rgas
      use physcons, only : avgd => con_avgd
      use machine,  only : kind_phys
      use idea_composition, only : amo, amn2, amo2
      use idea_tracer_mod,  only : oh, ho2, vmr_glob

      implicit none

! Argument
      integer, intent(in) :: im          ! number of long data points
      integer, intent(in) :: ix          ! max data points in fields
      integer, intent(in) :: levs        ! number of pressure levels
      integer, intent(in) :: ntrac_i     ! number of tracer add by IDEA
      real,    intent(in) :: dtp         ! time step in second
      real, intent(in) :: adt(ix,levs)   ! input  temp at dt=0

      real, intent(in) :: qin(ix,levs,ntrac_i)   ! input tracer

      real, intent(in) :: jo2(ix, levs)   ! input photo diss rate
      real, intent(in) :: n1(ix,levs)     ! number density of o
      real, intent(in) :: n2(ix,levs)     ! number density of o2
      real, intent(in) :: n(ix,levs)      ! number density of mixture
      real, intent(in) :: rho(ix,levs)    ! density of mixture
      real, intent(out):: dq(ix,levs,ntrac_i) ! output
! Local variables
      real, dimension(levs) :: noh, nh, nho2, natom, ndens
      real k1,k2,p1,p2,L1,L2
      real :: k3, k4, k5
      real :: mo,mo2,mn2
      integer k,i
      real :: q1new, q2new

      mo=amo*1.e-3/avgd
      mo2=amo2*1.e-3/avgd
      mn2=amn2*1.e-3/avgd

! O2 +hv =JJ => O + O
! O + O+ M  =k1 => O2+M 
! O + O2    =k2=>  O3
! O +OH  =k3=> O2 +H 
! O + O     =k4=>  O2 ?
! O +HO2  =k5=> O2 +OH 
! coefficents
        k3=4.2e-17
        k4=1.0e-26
        k5=3.5e-17

      do k=1,levs
         do i=1,im
            k1=4.7e-45*(300./adt(i,k))**2
            k2=6.e-46*(300./adt(i,k))**(2.4)

            p1=2.*Jo2(i,k)*n2(i,k) * mo/rho(i,k) !O-production mmr
            L1=2.*k1*n1(i,k)*n(i,k)+k2*n2(i,k)*n(i,k) ! O-loss
     &           +k3*oh(k)+k5*ho2(k)+2.*k4*n1(i,k)
            q1new=(qin(i,k,1)+p1*dtp)/(1.+L1*dtp)

            dq(i,k,1)=q1new-qin(i,k,1)
! conserve oxygen
            dq(i,k,2)= -dq(i,k,1)

        enddo
      enddo

      return
      end subroutine idea_tracer_c

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!   
      subroutine idea_tracer_eddy(im,ix,levs,ntrac_i,grav,prsi,prsl,
     &     rho,dtp,qin,dayno,dq)

! Calculate major species changes by eddy mixing O, O2, and
!     (indirectly) N2
! October 2017 Rashid Akmaev

      implicit none
! Arguments
      integer, intent(in) :: im    ! number of long data points in fields
      integer, intent(in) :: ix    ! max data points in fields
      integer, intent(in) :: levs  ! number of pressure levels
      integer, intent(in) :: ntrac_i ! number of addtl tracers in IDEA
      integer, intent(in) :: dayno ! for semiannual variation
      real, intent(in) :: dtp      ! time step in second
      real, intent(in) :: prsi(ix,levs+1) ! interface pressure in Pa
      real, intent(in) :: prsl(ix,levs)   ! layer pressure in Pa
      real, intent(in) :: grav(ix,levs)   ! (m/s**2)
      real, intent(in) :: rho(ix,levs)   ! mass density (kg/m**3)
      real, intent(in) :: qin(ix,levs,ntrac_i)   ! input tracers
      real, intent(out):: dq(ix,levs,ntrac_i) ! output tracer changes
! Locals
      real alpha(levs+1),beta(levs),qout(ntrac_i)
      real a(levs),b(levs),c(levs)
      real e(levs),d(levs,ntrac_i)
      integer i,k

!-----------------------------------------------------------------------
! Calculate Keddy (m**2/s) (move this to init subroutine/module later)
! Keddy parameters: mean, width in scale heights, height of max

      real, parameter:: pi = 3.141592653
!     real, parameter:: kmax = 120.
      real, parameter:: kmax = 140.
! semiannual amp
      real, parameter:: kampsa = 60.
!      real, parameter:: dkeddy = 2.
      real, parameter:: dkeddy = 0.
      real, parameter:: xmax = 15.
      real keddy(levs+1),x


      if(dkeddy <= 1e-10) then
!         keddy(:) = kmax
! Add semiannual variation
         keddy(:) = kmax + kampsa*(cos(4.*pi*(dayno+9.)/365.))
      else
         do k=1,levs+1
! height in scale heights
            x = alog(1e5/prsi(1,k))
            keddy(k)= kmax*exp(-((x-xmax)/dkeddy)**2)
         enddo
      endif
!      print *, 'kampsa=',kampsa,'keddy=', keddy(135)
!-----------------------------------------------------------------------
! Boundary conditions
      a(1) = 0.
      c(levs) = 0.

      do i = 1,im
! Auxiliary arrays
! at interfaces
         do k = 2,levs
            alpha(k) = keddy(k)*(.5*(rho(i,k-1)+rho(i,k)))**2*
     &           (.5*(grav(i,k-1)+grav(i,k)))/(prsl(i,k-1) -
     &           prsl(i,k))
         enddo

! in layers
         do k = 1,levs
            beta(k) = dtp*grav(i,k)/(prsi(i,k) - prsi(i,k+1))
         enddo

! Coefficients a(k), c(k) and b(k)
         do k = 2,levs
            a(k) = beta(k)*alpha(k)
         enddo
         do k = 1,levs-1
            c(k) = beta(k)*alpha(k+1)
         enddo
         do k = 1,levs
            b(k) = 1. + a(k) + c(k)
         enddo

! Solve tridiagonal problem for each tracer
! boundary conditions
         e(levs) = a(levs)/b(levs)
         d(levs,:) = qin(i,levs,:)/b(levs)

! go down, find e(k) and d(k)
         do k = levs-1,1,-1
            e(k) = a(k)/(b(k) - c(k)*e(k+1))
            d(k,:) = (c(k)*d(k+1,:) + qin(i,k,:))/(b(k) - c(k)*e(k+1))
         enddo
         
! go up, find solution
         qout(:) = d(1,:)
         dq(i,1,:) = qout(:) - qin(i,1,:)
         do k = 2,levs
            qout(:) = e(k)*qout(:) + d(k,:)
            dq(i,k,:) = qout(:) - qin(i,k,:)
         enddo

      enddo
      end subroutine idea_tracer_eddy

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!   
