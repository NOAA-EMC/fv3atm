
                        module module_control
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
use module_kinds
use module_constants
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---look-up tables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint),parameter :: &
 itb=201 &                   ! convection tables, dimension 1
,jtb=601 &                   ! convection tables, dimension 2
,kexm=10001 &                ! size of exponentiation table
,kztm=10001                  ! size of surface layer stability function table

real(kind=kfpt),parameter :: &
 ph=105000. &                ! upper bound of pressure range
,thh=350. &                  ! upper bound of potential temperature range
,thl=200.                    ! upper bound of potential temperature range

integer(kind=kint):: &
 kexm2 &                     ! internal points of exponentiation table
,kztm2                       ! internal pts. of the stability function table

real(kind=kfpt) :: &
 dex &                       ! exponentiation table step
,dzeta1 &                    ! sea table z/L step
,dzeta2 &                    ! land table z/L step
,fh01 &                      ! prandtl number for sea stability functions
,fh02 &                      ! prandtl number for land stability functions
,pl &                        ! lower bound of pressure range
,rdp &                       ! scaling factor for pressure
,rdq &                       ! scaling factor for humidity
,rdth &                      ! scaling factor for potential temperature
,rdthe &                     ! scaling factor for equivalent pot. temperature
,rdex &                      ! exponentiation table scaling factor
,xmax &                      ! upper bound for exponent in the table
,xmin &                      ! lower bound for exponent in the table
,ztmax1 &                    ! upper bound for z/L for sea stab. functions
,ztmin1 &                    ! lower bound for z/L for sea stab. functions
,ztmax2 &                    ! upper bound for z/L for land stab. functions
,ztmin2                      ! lower bound for z/L for land stab. functions

real(kind=kfpt),dimension(1:itb):: &
 sthe &                      ! range for equivalent potential temperature
,the0                        ! base for equivalent potential temperature

real(kind=kfpt),dimension(1:jtb):: &
 qs0 &                       ! base for saturation specific humidity
,sqs                         ! range for saturation specific humidity

real(kind=kfpt),dimension(1:kexm):: &
 expf                        ! exponentiation table

real(kind=kfpt),dimension(1:kztm):: &
 psih1 &                     ! sea heat stability function
,psim1 &                     ! sea momentum stability function
,psih2 &                     ! land heat stability function
,psim2                       ! land momentum stability function

real(kind=kfpt),dimension(1:itb,1:jtb):: &
 ptbl                        ! saturation pressure table

real(kind=kfpt),dimension(1:jtb,1:itb):: &
 ttbl                        ! temperature table
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---miscellaneous control parameters------------------------------------
!-----------------------------------------------------------------------
character(64):: &
 infile

logical(kind=klog):: &
 readbc &                    ! read regional boundary conditions
,runbc                       ! boundary data ready, start run

integer(kind=kint):: &
 ihr &                       ! current forecast hour
,ihrbc &                     ! boundary condition hour
,ihrstbc &                   ! boundary conditions starting time
,nbc                         ! boundary data logical unit

integer(kind=kint),dimension(1:3):: &
 idatbc(3)                   ! date of boundary data, day, month, year

integer(kind=kint):: &
 lnsbc                       ! # of boundary lines with enhanced diffusion
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
       contains
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine consts (pt)
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(kind=kfpt),intent(in) :: pt  ! Pressure at top of domain (Pa)
!
!-----------------------------------------------------------------------
      call exptbl
      pl=max(pt,1000.)
      call tablep
      call tablet
      call psitbl
!-----------------------------------------------------------------------
!
                        end subroutine consts
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
                        subroutine exptbl
!     ******************************************************************
!     *                                                                *
!     *               exponential function table                       *
!     *               responsible person: z.janjic                     *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 k                           ! index

real(kind=kfpt):: &
 x &                         ! argument
,xrng                        ! argument range
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      xmax= 30.
      xmin=-30.
!
      kexm2=kexm-2
      xrng=xmax-xmin
!
      dex=xrng/(kexm-1)
      rdex=1./dex
!--------------function definition loop---------------------------------
      x=xmin-dex
      do k=1,kexm
        x=x+dex
        expf(k)=exp(x)
      enddo
!-----------------------------------------------------------------------
!
                        end subroutine exptbl
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        function zjexp(x)
!     ******************************************************************
!     *                                                                *
!     *               exponential function table                       *
!     *               responsible person: z.janjic                     *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
real(kind=kfpt):: &
 zjexp

!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 k                           ! index

real(kind=kfpt):: &
 ak &                        ! position in table
,x                           ! argument
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      ak=(x-xmin)*rdex
      k=int(ak)
      k=max(k,0)
      k=min(k,kexm2)
!
      zjexp=(expf(k+2)-expf(k+1))*(ak-real(k))+expf(k+1)
!-----------------------------------------------------------------------
!
                        end function zjexp
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        subroutine tablep
!     ******************************************************************
!     *                                                                *
!     *    generates the table for finding pressure from               *
!     *    saturation specific humidity and potential temperature      *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 eps=1.e-10
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 kth &                       ! index
,kp                          ! index

real(kind=kfpt):: &
 ape &                       ! exner function
,dth &                       ! potential temperature step
,dp &                        ! pressure step
,dqs &                       ! saturation specific humidity step
,p &                         ! pressure
,qs0k &                      ! base value for saturation humidity
,sqsk &                      ! saturation spec. humidity range
,th                          ! potential temperature

real(kind=kfpt),dimension(1:itb):: &
 app &                       ! temporary
,aqp &                       ! temporary
,pnew &                      ! new pressures
,pold &                      ! old pressure
,qsnew &                     ! new saturation spec. humidity
,qsold &                     ! old saturation spec. humidity
,y2p                         ! temporary
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      dth=(thh-thl)/real(jtb-1)
      dp=(ph-pl)/real(itb-1)
      rdth=1./dth
!-----------------------------------------------------------------------
      th=thl-dth
      do kth=1,jtb
        th=th+dth
        p=pl-dp
        do kp=1,itb
          p=p+dp
          ape=(100000./p)**cappa
          qsold(kp)=pq0/p*exp(a2*(th-a3*ape)/(th-a4*ape))
          pold(kp)=p
        enddo
!
        qs0k=qsold(1)
        sqsk=qsold(itb)-qsold(1)
        qsold(1)=0.
        qsold(itb)=1.
!
        do kp=2,itb-1
          qsold(kp)=(qsold(kp)-qs0k)/sqsk
!wwwwwwwwwwwwww fix due to 32 bit precision limitation wwwwwwwwwwwwwwwww
          if((qsold(kp)-qsold(kp-1)).lt.eps) then
            qsold(kp)=qsold(kp-1)+eps
          endif
!mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
        enddo
!
        qs0(kth)=qs0k
        sqs(kth)=sqsk
!
        qsnew(1)=0.
        qsnew(itb)=1.
        dqs=1./real(itb-1)
        rdq=1./dqs
!
        do kp=2,itb-1
          qsnew(kp)=qsnew(kp-1)+dqs
        enddo
!
        y2p(1)=0.
        y2p(itb)=0.
!
        call spline(jtb,itb,qsold,pold,y2p,itb,qsnew,pnew,app,aqp)
!
        do kp=1,itb
          ptbl(kp,kth)=pnew(kp)
        enddo
!-----------------------------------------------------------------------
      enddo
!-----------------------------------------------------------------------
!
                        end subroutine tablep
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
                        subroutine tablet
!     ******************************************************************
!     *                                                                *
!     *    generates the table for finding temperature from            *
!     *    pressure and equivalent potential temperature               *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 eps=1.e-10
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 kth &                       ! index
,kp                          ! index

real(kind=kfpt):: &
 ape &                       ! exner function
,dth &                       ! potential temperature step
,dp &                        ! pressure step
,dthe &                      ! equivalent pot. temperature step
,p &                         ! pressure
,qs &                        ! saturation specific humidity
,the0k &                     ! base value for equivalent pot. temperature
,sthek &                     ! equivalent pot. temperature range
,th                          ! potential temperature

real(kind=kfpt),dimension(1:jtb):: &
 apt &                       ! temporary
,aqt &                       ! temporary
,tnew &                      ! new temperature
,told &                      ! old temperature
,thenew &                    ! new equivalent potential temperature
,theold &                    ! old equivalent potential temperature
,y2t                         ! temporary
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      dth=(thh-thl)/real(jtb-1)
      dp=(ph-pl)/real(itb-1)
      rdp=1./dp
!-----------------------------------------------------------------------
      p=pl-dp
      do kp=1,itb
        p=p+dp
        th=thl-dth
        do kth=1,jtb
          th=th+dth
          ape=(100000./p)**cappa
          qs=pq0/p*exp(a2*(th-a3*ape)/(th-a4*ape))
          told(kth)=th/ape
          theold(kth)=th*exp(eliwv*qs/(cp*told(kth)))
        enddo
!
        the0k=theold(1)
        sthek=theold(jtb)-theold(1)
        theold(1)=0.
        theold(jtb)=1.
!
        do kth=2,jtb-1
          theold(kth)=(theold(kth)-the0k)/sthek
!wwwwwwwwwwwwww fix due to 32 bit precision limitation wwwwwwwwwwwwwwwww
          if((theold(kth)-theold(kth-1)).lt.eps) then
            theold(kth)=theold(kth-1)+eps
          endif
!mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
        enddo
!
        the0(kp)=the0k
        sthe(kp)=sthek
!
        thenew(1)=0.
        thenew(jtb)=1.
        dthe=1./real(jtb-1)
        rdthe=1./dthe
!
        do kth=2,jtb-1
          thenew(kth)=thenew(kth-1)+dthe
        enddo
!
        y2t(1)=0.
        y2t(jtb)=0.
!
        call spline(jtb,jtb,theold,told,y2t,jtb,thenew,tnew,apt,aqt)
!
        do kth=1,jtb
          ttbl(kth,kp)=tnew(kth)
        enddo
!-----------------------------------------------------------------------
      enddo
!-----------------------------------------------------------------------
!
                        end subroutine tablet
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      subroutine spline(jtbx,nold,xold,yold,y2,nnew,xnew,ynew,p,q)
!   ********************************************************************
!   *                                                                  *
!   *  this is a one-dimensional cubic spline fitting routine          *
!   *  programed for a small scalar machine.                           *
!   *                                                                  *
!   *  programer z. janjic                                             *
!   *                                                                  *
!   *  nold - number of given values of the function.  must be ge 3.   *
!   *  xold - locations of the points at which the values of the       *
!   *         function are given.  must be in ascending order.         *
!   *  yold - the given values of the function at the points xold.     *
!   *  y2   - the second derivatives at the points xold.  if natural   *
!   *         spline is fitted y2(1)=0. and y2(nold)=0. must be        *
!   *         specified.                                               *
!   *  nnew - number of values of the function to be calculated.       *
!   *  xnew - locations of the points at which the values of the       *
!   *         function are calculated.  xnew(k) must be ge xold(1)     *
!   *         and le xold(nold).                                       *
!   *  ynew - the values of the function to be calculated.             *
!   *  p, q - auxiliary vectors of the length nold-2.                  *
!   *                                                                  *
!   ********************************************************************
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer(kind=kint),intent(in):: &
       jtbx,nnew,nold

      real(kind=kfpt),dimension(jtbx),intent(in):: &
       xnew,xold,yold

      real(kind=kfpt),dimension(jtbx),intent(inout):: &
       p,q,y2

      real(kind=kfpt),dimension(jtbx),intent(out):: &
       ynew
!
      integer(kind=kint):: &
       k,k1,k2,kold,noldm1

      real(kind=kfpt):: &
       ak,bk,ck,den,dx,dxc,dxl,dxr,dydxl,dydxr &
      ,rdx,rtdxc,x,xk,xsq,y2k,y2kp1
!-----------------------------------------------------------------------
      noldm1=nold-1
!
      dxl=xold(2)-xold(1)
      dxr=xold(3)-xold(2)
      dydxl=(yold(2)-yold(1))/dxl
      dydxr=(yold(3)-yold(2))/dxr
      rtdxc=0.5/(dxl+dxr)
!
      p(1)= rtdxc*(6.*(dydxr-dydxl)-dxl*y2(1))
      q(1)=-rtdxc*dxr
!
      if(nold==3)go to 150
!-----------------------------------------------------------------------
      k=3
!
  100 dxl=dxr
      dydxl=dydxr
      dxr=xold(k+1)-xold(k)
      dydxr=(yold(k+1)-yold(k))/dxr
      dxc=dxl+dxr
      den=1./(dxl*q(k-2)+dxc+dxc)
!
      p(k-1)= den*(6.*(dydxr-dydxl)-dxl*p(k-2))
      q(k-1)=-den*dxr
!
      k=k+1
      if(k<nold)go to 100
!-----------------------------------------------------------------------
  150 k=noldm1
!
  200 y2(k)=p(k-1)+q(k-1)*y2(k+1)
!
      k=k-1
      if(k>1)go to 200
!-----------------------------------------------------------------------
      k1=1
!
  300 xk=xnew(k1)
!
      do 400 k2=2,nold
!
      if(xold(k2)>xk)then
        kold=k2-1
        go to 450
      endif
!
  400 continue
!
      ynew(k1)=yold(nold)
      go to 600
!
  450 if(k1==1)go to 500
      if(k==kold)go to 550
!
  500 k=kold
!
      y2k=y2(k)
      y2kp1=y2(k+1)
      dx=xold(k+1)-xold(k)
      rdx=1./dx
!
      ak=.1666667*rdx*(y2kp1-y2k)
      bk=0.5*y2k
      ck=rdx*(yold(k+1)-yold(k))-.1666667*dx*(y2kp1+y2k+y2k)
!
  550 x=xk-xold(k)
      xsq=x*x
!
      ynew(k1)=ak*xsq*x+bk*xsq+ck*x+yold(k)
!
  600 k1=k1+1
      if(k1<=nnew)go to 300
!-----------------------------------------------------------------------
      end subroutine spline
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine psitbl
!     ******************************************************************
!     *                                                                *
!     *               surface layer integral functions                 *
!     *               responsible person: z.janjic                     *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 eps=0.000001
!--local variables------------------------------------------------------
integer(kind=kint):: &
 k                           ! index

real(kind=kfpt):: &
 x &                         ! temporary
,zeta1 &                     ! z/L, sea
,zeta2 &                     ! z/L, land
,zrng1 &                     ! z/L range, sea
,zrng2                       ! z/L range, land
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      kztm2=kztm-2
!
      fh01=1.
      fh02=1.
!
      ztmin1=-5.0
      ztmax1= 1.0
!
      ztmin2=-5.0
      ztmax2= 1.0
!
      zrng1=ztmax1-ztmin1
      zrng2=ztmax2-ztmin2
!
      dzeta1=zrng1/(kztm-1)
      dzeta2=zrng2/(kztm-1)
!--------------function definition loop---------------------------------
      zeta1=ztmin1
      zeta2=ztmin2
      do k=1,kztm
!--------------unstable range-------------------------------------------
        if(zeta1.lt.0.)then
!--------------paulson 1970 functions-----------------------------------
          x=sqrt(sqrt(1.-16.*zeta1))
          psim1(k)=-2.*log((x+1.)/2.)-log((x*x+1.)/2.)+2.*atan(x)-pihf
          psih1(k)=-2.*log((x*x+1.)/2.)
!--------------stable range---------------------------------------------
        else
!--------------holtslag and de bruin 1988-------------------------------
          psim1(k)=0.7*zeta1+0.75*zeta1*(6.-0.35*zeta1)*exp(-0.35*zeta1)
          psih1(k)=0.7*zeta1+0.75*zeta1*(6.-0.35*zeta1)*exp(-0.35*zeta1)
!-----------------------------------------------------------------------
        endif
!--------------unstable range-------------------------------------------
        if(zeta2.lt.0.)then
!--------------paulson 1970 functions-----------------------------------
          x=sqrt(sqrt(1.-16.*zeta2))
          psim2(k)=-2.*log((x+1.)/2.)-log((x*x+1.)/2.)+2.*atan(x)-pihf
          psih2(k)=-2.*log((x*x+1.)/2.)
!--------------stable range---------------------------------------------
        else
!--------------holtslag and de bruin 1988-------------------------------
          psim2(k)=0.7*zeta2+0.75*zeta2*(6.-0.35*zeta2)*exp(-0.35*zeta2)
          psih2(k)=0.7*zeta2+0.75*zeta2*(6.-0.35*zeta2)*exp(-0.35*zeta2)
!-----------------------------------------------------------------------
        endif
!-----------------------------------------------------------------------
        if(k.eq.kztm)then
          ztmax1=zeta1
          ztmax2=zeta2
        endif
!
        zeta1=zeta1+dzeta1
        zeta2=zeta2+dzeta2
!-----------------------------------------------------------------------
      enddo
!-----------------------------------------------------------------------
      ztmax1=ztmax1-eps
      ztmax2=ztmax2-eps
!-----------------------------------------------------------------------
!
                        endsubroutine psitbl
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
                       end module module_control
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
