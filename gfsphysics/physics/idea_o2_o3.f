
      subroutine idea_o2_o3(im,ix,levs,cosz,adt,o2_n,o3_n,rho,cp,       
     &zg,grav,dth)
!
! Apr 06 2012  Henry Juang, initial implement for nems
! Jan 02 2013  Jun Wang,    move o3ini out of column physics
! Nov 29 2016  VAY review/suggestions and SRB/SRC-out
      use physcons, only:    avgd=>con_avgd  !, amo3=> con_amo3 
      use idea_composition
      use idea_solar, only : eff_hugg, eff_chap, eff_herz, eff_hart
!
      implicit none
! Argument
      integer, intent(in) :: im               ! number of data points in adt (first dim)
      integer, intent(in) :: ix               ! max data points in adt (first dim)
      integer, intent(in) :: levs             ! number of pressure levels
      real, intent(in)    :: cosz(im)         ! cos zenith angle
      real, intent(in)    :: adt(ix,levs)     !temp(k) 
      real, intent(in)    :: o2_n(ix,levs)    ! /m3
      real, intent(in)    :: o3_n(ix,levs)    ! /m3
      real, intent(in)    :: rho(ix,levs)     ! kg/m3
      real, intent(in)    :: cp(ix,levs)      ! J/kg/k
      real, intent(in)    :: zg(ix,levs)      ! height (m)
      real, intent(in)    :: grav(ix,levs)    ! (m/s2)
      real, intent(out)   :: dth(ix,levs)     ! heating rate k/s
!
      real hc,fc,dc,hha,fha,dha,hhu,i1,i2,m,dhu,lams,laml               
     &,hhz,fhz,dhzo2,dhzo3,hsrb,fsrb,dsrb,ysrb,h1,rodfac
      real clmo2(levs),clmo3(levs) 
      real  :: rdzg
      real  :: tg_vay
      real  :: fo2_vay, fo3_vay
      real  :: fc_dc, fha_dha
      real  :: hu_exp1, hu_exp2, rm, i2_i1
      integer i,k
!
! very ....."dangerous" games with constants !!!!
!
      fc=370.     !J/m2/s
      dc=2.85E-25 !m2
      fc_dc = fc*dc
      fha=5.13  !J/m2/s
      dha=8.7E-22 !m2
      fha_dha = fha*dha
      i1=0.07   !J/m2/s/A
      i2=0.05
      m=0.01273   !/A
      rm = 1./m

      lams=2805.
      laml=3015.
      dhu=1.15e-6    !m2
      fhz=1.5        !J/m2/s
      dhzo2=6.e-28  !m2
      dhzo3=4.e-22  !m2
      fsrb=0.0128    !J/m2/s
      dsrb=2.07e-24  !m2
      ysrb=0.0152

       hu_exp1 = exp(-m*laml)*dhu
       hu_exp2 = exp(-m*lams)*dhu
       i2_i1 = -0.02
!
! rewrite and compute constants in the init_solar (?)
!   exp(-1.*m*laml) in loops
!   mmr*dp => <n>dZ
      dth(:,:)=0.
      fo3_vay = 1000.*avgd/amo3*bz
      fo2_vay = 1000.*avgd/amo2*bz
      do i=1,im
        if(cosz(i).ge.0.) then
          rodfac=35./sqrt(1224.*cosz(i)**2+1.)
          tg_vay = adt(i,levs)/grav(i,levs)
!         clmo2(levs)=1.e3*o2_n(i,levs)*bz*adt(i,levs)*avgd/(grav(i,levs)*amo2)
          clmo2(levs)=o2_n(i,levs)*tg_vay*fo2_vay
          clmo3(levs)=o3_n(i,levs)*tg_vay*fo3_vay
          do k=levs-1,1,-1
           rdzg =  .5*(zg(i,k+1)-zg(i,k))
          clmo2(k)=clmo2(k+1)+(o2_n(i,k+1)+o2_n(i,k))*rdzg                                     
          clmo3(k)=clmo3(k+1)+(o3_n(i,k+1)+o3_n(i,k))*rdzg              
!
          enddo
          clmo2=clmo2*rodfac   !rad path
          clmo3=clmo3*rodfac
!
!very ....."dangerous" games with constants !!!!exp(-1.*m*laml)  i2/i1 etc....  
!            can be restricted to MLT
!
          do k=1,levs
!
!Chappius acc-cy 2% according Strobel 1978
!
            hc=fc_dc*exp(-dc*clmo3(k))
! Hartley
            hha=fha_dha*exp(-dha*clmo3(k))
! Huggins
            hhu=rm*(i1+i2_i1*exp(-clmo3(k)*hu_exp1)       
     &      -i2*exp(-clmo3(k)*hu_exp2)) /clmo3(k)
! Herzberg
            hhz=fhz*(dhzo2*o2_n(i,k)+dhzo3*o3_n(i,k))*         
     &        exp(-dhzo2*clmo2(k)-dhzo3*clmo3(k))
!
! VAY-2016: hsrb 2 times above P > 0.02 orchetrating with SOLAR in idea_solar_heating
!           moved to idea_solar_heating.f
! Comments for Zhu's updates
!vay            h1=sqrt(1.+4.*dsrb*clmo2(k)/(pi*ysrb))
!vay            hsrb=fsrb*dsrb*o2_n(i,k)*exp(-.5*pi*ysrb*(h1-1.))/h1
!           dth(i,k)=((hc+hha+hhu)*o3_n(i,k)+hhz+hsrb)/                 
!    &             (cp(i,k)*rho(i,k))
        dth(i,k)=((hc+hha*eff_hart(k)+hhu)*o3_n(i,k)+hhz)/           
     &             (cp(i,k)*rho(i,k))
          enddo
        else
          dth(i,1:levs)=0.
        endif
      enddo
      return
      end
      subroutine o3pro(im,ix,levs,ntrac,adr,am,n,o3_n)
!
!      use physcons, amo3=> con_amo3, avgd=> con_avgd 
      use idea_composition
!
      implicit none
! IN/OUT
      integer, intent(in) :: im                    ! number of data points in PE
      integer, intent(in) :: ix                    ! max data points reserved & fixed ngptc
      integer, intent(in) :: levs                  ! number of pressure levels
      integer, intent(in) :: ntrac                 ! number of tracer
      real, intent(in)    :: adr(ix,levs,ntrac)    ! tracers in "mmr"
      real, intent(in)    :: am(ix,levs)           ! mixture mol weight kg
      real, intent(in)    :: n(ix,levs)            ! number density  /m3
      real, intent(out)   :: o3_n(ix,levs)         ! /m3
! locals
      real rate
      integer i,k
!
!      mo3=amo3*1.e-3/avgd
!      rmo3 = 1000.*avgd/amo3
      do i=1,im
        do k=1,k71-1
          o3_n(i,k)=adr(i,k,2)*am(i,k)*n(i,k)*rmo3
        enddo
          rate=adr(i,k71,2)/o3ra(k71)
        do k=k71,levs
          o3_n(i,k)=o3ra(k)*rate*am(i,k)*n(i,k)*rmo3
        enddo
      enddo
      return
      end
      subroutine o3ini(levs)
!
      use idea_composition
!
      implicit none
      integer,intent(in) :: levs  ! number of pressure levels
      integer i
      real c0(2),c1(2),c2(2),c3(2),logp,x
!=====================================
!MS-93
!     data c0/0.66965,0.92621/
!     data c1/-0.009682,0.13396/
!     data c2/0.033093,-0.076863/
!     data c3/0.017938,0.006897/
!?? Why small corrections
      data c0/0.66965,0.932363/
      data c1/-0.009682,0.139425/
      data c2/0.033093,-0.076863/
      data c3/0.017938,0.005075/
!
      allocate (ef(levs))
!
      do i=1,levs
        logp=log10(pr_idea(i))
          if(logp.ge.0.) then
            ef(i)=c0(2)+c1(2)+c2(2)+c3(2)      ! vay should be ~1 or 0.99
          elseif(logp.ge.-2) then
            x=1.+logp
            ef(i)=c0(2)+c1(2)*x+c2(2)*x**2+c3(2)*x**3  
          elseif(logp.ge.-4) then
            x=3.+logp
            ef(i)=c0(1)+c1(1)*x+c2(1)*x**2+c3(1)*x**3  
          else
            ef(i)=c0(1)-c1(1)+c2(1)-c3(1)  
          endif
      enddo
      return
      end
!
!vay-2016
!
      subroutine  heat_uveff(levs, eff_hart, eff_src)
!
      use idea_composition, only : pr_idea
!
      implicit none
      integer,intent(in) :: levs             ! number of pressure levels
      real ,  intent(out) :: eff_hart(levs)  ! heat-eff by ML-93 for Hartley
      real ,  intent(out) :: eff_src(levs)  !

      integer i
      real c0(3),c1(3),c2(3),c3(3),lp,x
      real EFF_SRC17(17), EFF_EUV17(17)
!=====================================
!MS-93
      data c0/0.66965,0.92621,    0.75349/
      data c1/-0.009682,0.13396,  0.0036/
      data c2/0.033093,-0.076863, 0.0595/
      data c3/0.017938,0.006897, -0.0228/
!
! pr_idea in mb
!
!     DATA EFF_SRC17/5*.28,.29,.32,.38,.4,.4,.4,.39,.34,.26,.19,.17,.16/         
      do i=1,levs
        lp=log10(pr_idea(i))       
 
         if(lp.ge.0.) then
            eff_hart(i)=c0(2)+c1(2)+c2(2)+c3(2)      ! vay should be ~1 or 0.99
          elseif(lp.ge.-2) then
            x=1.+lp
            eff_hart(i)=c0(2)+c1(2)*x+c2(2)*x*x+c3(2)*x*x*x  
          elseif(lp.ge.-4) then
            x=3.+lp
            eff_hart(i)=c0(1)+c1(1)*x+c2(1)*x*x+c3(1)*x*x*x  
            eff_src(i)=c0(3)+c1(3)*x+c2(3)*x*x+c3(3)*x*x*x  
          else
            eff_hart(i)=c0(1)-c1(1)+c2(1)-c3(1) 
            eff_src(i)=c0(3)-c1(3)+c2(3)-c3(3)  
          endif
       enddo
      return
      end  subroutine heat_uveff
