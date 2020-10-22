      subroutine idea_getcoef(levs,ntrac,q,plyr,mur,lam,d12)

! calculate global avg viscosity, conductivity, and diffusion coeffs
! Apr 06 2012    Henry Juang, initial implement
!

      use physcons, only:  amo2 => con_amo2,
     &                     avgd => con_avgd, p0 => con_p0,
     &                     amh2o =>con_amw, amo3 =>con_amo3
      use tracer_const, only: ri, cpi
      use idea_composition, only: am=>amgm
      implicit none

! subroutine params

      integer, intent(in) :: levs   ! number of pressure model layers
      integer, intent(in) :: ntrac   ! number of tracers
      real,    intent(in) :: q(levs,0:ntrac) ! temp*cp,h2o,o3,cld,o,o2 tracer
      real,    intent(in) :: plyr(levs) ! global mean presure at model layer
      real,    intent(out):: mur(levs) ! mu/rho (m2/s)
      real,    intent(out):: lam(levs) ! lambda/rho/cp (m2/s)
      real,    intent(out):: d12(levs) ! d12/n

! local params

      real amo,amn2,muo,muo2,mun2,lao,lao2,lan2
      parameter (amo=15.9994, amn2=28.013) !g/mol
      parameter (muo=3.9e-7, muo2=4.03e-7, mun2=3.43e-7) !kg/m/s
      parameter (lao=75.9e-5, lao2=56.e-5, lan2=56.e-5) !kg/m/s         
      real, parameter:: bz=1.3806505e-23 ! Boltzmann constant 
      real, parameter:: s12=0.774,a12=9.69e18 ! O-O2 diffusion params

! local variables

      real hold1,la,rho,o_n,o2_n,cp,t69,mu,n,n2_n,q_n2,tem,rn
      integer k
      logical, save :: skprnt = .false.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     print *,'in idea_getcoef,cp=',cpi(1:5)
      do k=1,levs

! get global mean pressure

!hmhj    plyr=(ak5(k)+ak5(k+1)+p0*1.e-3*(bk5(k)+bk5(k+1)))*500.

! get n2 kg/kg

       if( ntrac.eq.5 ) then
         q_n2=1.-q(k,1)-q(k,2)-q(k,4)-q(k,5)
       else
         q_n2=1.
         do n=1,ntrac
           if(cpi(n).ne.0.0) q_n2=q_n2-q(k,n)
         enddo
       endif

! get cp, recover temp from enthalpy

       if( ntrac.eq.5 ) then
         cp=cpi(1)*q(k,1)+cpi(2)*q(k,2)+cpi(4)*q(k,4)+cpi(5)*q(k,5)+
     &        cpi(0)*q_n2
       else
         cp=cpi(0)*q_n2
         do n=1,ntrac
           if(cpi(n).ne.0.0) cp=cp+cpi(n)*q(k,n)
         enddo
       endif
       if(abs(cp)<1.0e-10) then
        print *,'in idea_getcoef,k=',k,'q(k1)=',q(k,1),q(k,2),
     &   q(k,3),q(k,4),q(k,5),'q_n2=',q_n2
       endif
         tem=q(k,0)/cp    ! K

! get molecular mass, number densities

       if( ntrac.eq.5 ) then
         am(k)=1./(q(k,4)/amo+q(k,5)/amo2+q_n2/amn2+q(k,1)/amh2o+
     &        q(k,2)/amo3)                      ! g/mol
         n=plyr(k)/tem                          ! 1/m3/bz
         rn=avgd*bz
         rho=1e-3*am(k)*n/rn                     ! kg/m3
         o_n=q(k,4)*am(k)/amo            
         o2_n=q(k,5)*am(k)/amo2         
         n2_n=q_n2*am(k)/amn2          
         hold1=tem/plyr(k)
!! get d12

         d12(k)=(a12*bz)*tem**(s12)*hold1      !d12 

! calculate mass nu=mu/rho,lamda/(rho*cp)

         mu=o_n*muo+o2_n*muo2+n2_n*mun2
         la=o_n*lao+o2_n*lao2+n2_n*lan2

! now use tem**0.69

         t69=tem**(0.69)

         mur(k)=mu*t69/rho
         lam(k)=la*t69/rho/cp

!      print*,'ntrac=5,k=',k,plyr(k),mur(k),lam(k),d12(k),
!     &  'q(k,0)=',q(k,0),'q(k,1)=',q(k,1),'q(k,2)=',q(k,2),
!     &  'q(k,3)=',q(k,3),'q(k,4)=',q(k,4),'q(k,5)=',q(k,5),
!     &  'cp=',cp
!
!------------------------------------------------------
       else  ! assume ntrac=3
         am(k)=1./(q_n2/amn2+q(k,1)/amh2o+q(k,2)/amo3)
         o_n=0.0
         o2_n=0.0
!jw
         n=plyr(k)/tem                          ! 1/m3/bz
         rn=bz*avgd                          ! 1/m3
         rho=1e-3*am(k)*n/rn                   ! kg/m3
         n2_n=q_n2*am(k)/amn2              !
         hold1=tem/plyr(k)

!       print *,'ntrac=3,k am rho n2_n ',k,am(k),rho,n2_n,n,rn,tem,      &
!   &  'q=',q_n2,q(k,1:ntrac),cp,'cpi=',cpi(0:ntrac),q(k,0)

! get d12

         d12(k)=(a12*bz)*tem**(s12)*hold1   

! calculate mass nu=mu/rho,lamda/(rho*cp)

         mu=n2_n*mun2
         la=n2_n*lan2

! now use tem**0.69

         t69=tem**(0.69)
         
         mur(k)=mu*t69/rho
         lam(k)=la*t69/rho/cp

       endif

      enddo
!SK2020Oct5
      if (skprnt) print *,' skprnt:in idea_getcoef,am=',am(1:levs)
      skprnt = .false.
!ske
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      end subroutine idea_getcoef
