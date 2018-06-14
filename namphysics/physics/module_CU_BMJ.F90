!-----------------------------------------------------------------------
!
      module module_cu_bmj
!
!-----------------------------------------------------------------------
!
!***  the convection drivers and packages
!
!-----------------------------------------------------------------------
!
      use module_kinds
      use module_constants ,only : &
      a2,a3,a4,cappa,cp,eliv,elwv,epsq,g,p608,pq0,r_d,tiw
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      private
!
      public :: bmj_init,bmjdrv
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  for bmj convection
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      real(kind=kfpt),parameter:: &
       dttop=0.,efifc=5.0,efimn=0.20 &
      ,efmntl=0.70,efmnts=0.70 &
      ,eliwv=2.683e6,enplo=98500.,enpup=95000. &
      ,epsdn=1.05,epsdt=0. &
      ,epsntp=.0001,epsntt=.0001,epspr=1.e-7 &
      ,epsup=1.00 &
      ,fup=1./200000. &
      ,pbm=13000.,pfrz=15000.,pno=1000. &
      ,pone=2500.,pqm=20000. &
      ,psh=20000.,pshu=45000. &
      ,rendp=1./(enplo-enpup) &
      ,rhlsc=0.00,rhhsc=1.00 &
      ,row=1.e3 &
      ,stabdf=0.90,stabds=0.90 &
      ,stabs=1.0,stresh=1.10 &
      ,dtshal=-1.0,trel=2400.
!
      real(kind=kfpt),parameter:: &
       dttrigr=-0.0 &
      ,dtptrigr=dttrigr*pone  !<-- average parcel virtual temperature deficit over depth pone
                              !<-- note: capetrigr is scaled by the cloud base temperature 
                              !    (see below)
!
      real(kind=kfpt),parameter:: &
       dspbfl_base=-3875. &
      ,dsp0fl_base=-5875. &
      ,dsptfl_base=-1875. &
      ,dspbfs_base=-3875. &
      ,dsp0fs_base=-5875. &
      ,dsptfs_base=-1875.
!
      real(kind=kfpt),parameter:: &
       pl=2500.,plq=70000.,ph=105000. &
      ,thl=210.,thh=365.,thhq=325.
!
      integer(kind=kint),parameter:: &
       itb=76,jtb=134,itbq=152,jtbq=440
!
      integer(kind=kint),parameter:: &
       itrefi_max=3
!
!***  arrays for lookup tables
!
      real(kind=kfpt),dimension(itb),private,save:: &
       sthe,the0

      real(kind=kfpt),dimension(jtb),private,save:: &
       qs0,sqs

      real(kind=kfpt),dimension(itbq),private,save:: &
       stheq,the0q

      real(kind=kfpt),dimension(itb,jtb),private,save:: &
       ptbl

      real(kind=kfpt),dimension(jtb,itb),private,save:: &
       ttbl

      real(kind=kfpt),dimension(jtbq,itbq),private,save:: &
       ttblq
                         
!***  share copies for module_bl_myjpbl
!
      real(kind=kfpt),dimension(jtb):: &
       qs0_exp,sqs_exp

      real(kind=kfpt),dimension(itb,jtb):: &
       ptbl_exp
!
      real(kind=kfpt),parameter:: &
       rdp=(itb-1.)/(ph-pl),rdpq=(itbq-1.)/(ph-plq)   &
      ,rdq=itb-1,rdth=(jtb-1.)/(thh-thl) &
      ,rdthe=jtb-1.,rdtheq=jtbq-1. &
      ,rsfcp=1./101300.
!
      real(kind=kfpt),parameter:: &
       avgefi=(efimn+1.)*0.5 &
      ,stefi=1.
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      subroutine bmjdrv( &
                        ims,ime,jms,jme,lm &
                       ,entrain,newall,newswap,newupup,nodeep &
                       ,fres,fr,fsl,fss &
                       ,dt,ntsd,ncnvc &
                       ,raincv,cutop,cubot &
                       ,th,t,q,u_phy,v_phy,dudt_phy,dvdt_phy &
                       ,phint,phmid,exner &
                       ,cldefi,xland,cu_act_flag &
                     ! optional
                       ,rthcuten,rqcuten &
                       )
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
      logical(kind=klog),intent(in):: &
       entrain,newall,newswap,newupup,nodeep
!
      logical(kind=klog),dimension(ims:ime,jms:jme),intent(inout):: &
       cu_act_flag
!
      integer(kind=kint),intent(in):: &
       ims,ime,jms,jme,lm &
      ,ntsd,ncnvc
!
      real(kind=kfpt),intent(in):: &
       dt,fres,fr,fsl,fss
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
       xland
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
       q &
      ,exner,phmid,t,th,u_phy,v_phy
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm+1),intent(in):: &
       phint
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm) &
                     ,optional &
                     ,intent(inout):: &
       rqcuten,rthcuten 
! 
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout):: &
       cldefi,raincv
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(out):: &
       cubot,cutop

      real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out):: &
       dudt_phy,dvdt_phy
!
!-----------------------------------------------------------------------
!***
!***  local variables
!***
!-----------------------------------------------------------------------
      integer(kind=kint):: &
       i,icldck,ierr,j,k,lbot,lmh,ltop
!
      real(kind=kfpt):: &
       dspbfl,dsp0fl,dsptfl,dspbfs,dsp0fs,dsptfs &
      ,dspbsl,dsp0sl,dsptsl,dspbss,dsp0ss,dsptss &
      ,slopbl,slop0l,sloptl,slopbs,slop0s,slopts &
      ,dtcnvc,seamask,pcpcol,psfc,ptop,qnew,qold
! 
      real(kind=kfpt),dimension(1:lm):: &
       dpcol,dqdt,dtdt,dudt,dvdt,pcol,qcol,tcol,exnercol,ucol,vcol
!
!***  begin debugging convection
      real(kind=kfpt) :: delq,delt,plyr
      integer(kind=kint) :: imd,jmd
      logical(kind=klog) :: print_diag
!***  end debugging convection
!
!-----------------------------------------------------------------------
!*********************************************************************** 
!-----------------------------------------------------------------------
!
!***  prepare to call bmj convection scheme
!
!-----------------------------------------------------------------------
!
!***  check to see if this is a convection timestep
!                                                                        
      icldck=mod(ntsd,ncnvc)                                              
!-----------------------------------------------------------------------
!                                                                      
!***  compute convection every ncnvc*dt/60.0 minutes
!                                                                     
!***  begin debugging convection
      imd=(ims+ime)/2
      jmd=(jms+jme)/2
      print_diag=.false.
!***  end debugging convection

      if(icldck==0.or.ntsd==0)then
!
        do j=jms,jme
        do i=ims,ime
          cu_act_flag(i,j)=.true.
        enddo
        enddo
!
        dtcnvc=dt*ncnvc
!.......................................................................
!zj$omp parallel do &
!zj$omp private (j,i,dqdt,dtdt,dudt,dvdt &
!zj$omp         ,dxh,pcpcol,psfc,ptop,seamask,k &
!zj$omp         ,tcol,pcol,dpcol,qcol,ucol,vcol &
!zj$omp         ,lmh,delt,delq,plyr,lbot,ltop)
!.......................................................................
        do j=jms,jme
          do i=ims,ime
!--set up deficit saturation pressures depending on resolution----------
            dspbfl=dspbfl_base*fres*fr
            dsp0fl=dsp0fl_base*fres*fr
            dsptfl=dsptfl_base*fres*fr
            dspbfs=dspbfs_base*fres
            dsp0fs=dsp0fs_base*fres
            dsptfs=dsptfs_base*fres
!
            dspbsl=dspbfl*fsl
            dsp0sl=dsp0fl*fsl
            dsptsl=dsptfl*fsl
            dspbss=dspbfs*fss
            dsp0ss=dsp0fs*fss
            dsptss=dsptfs*fss
!
            slopbl=(dspbfl-dspbsl)/(1.-efimn)
            slop0l=(dsp0fl-dsp0sl)/(1.-efimn)
            sloptl=(dsptfl-dsptsl)/(1.-efimn)
            slopbs=(dspbfs-dspbss)/(1.-efimn)
            slop0s=(dsp0fs-dsp0ss)/(1.-efimn)
            slopts=(dsptfs-dsptss)/(1.-efimn)
!------------------------------------------------------------------------
          do k=1,lm
            dqdt(k)=0.
            dtdt(k)=0.
            dudt(k)=0.
            dvdt(k)=0.
          enddo
!
          raincv(i,j)=0.
          pcpcol=0.
          psfc=phint(i,j,lm+1)
          ptop=phint(i,j,1)
!
!***  convert to bmj land mask (1.0 for sea; 0.0 for land)
!
          seamask=xland(i,j)-1.
!
!***  fill 1-d vertical arrays 
!
          do k=1,lm
!
            ucol    (k)=u_phy(i,j,k)
            vcol    (k)=v_phy(i,j,k)
            qcol    (k)=max(epsq,q(i,j,k))
            tcol    (k)=t(i,j,k)
            pcol    (k)=phmid(i,j,k)
            exnercol(k)=exner(i,j,k)
            dpcol   (k)=phint(i,j,k+1)-phint(i,j,k)
          enddo
!
!***  lowest layer above ground must also be flipped
!
          lmh=lm
!-----------------------------------------------------------------------
!***
!***  call convection
!***
!-----------------------------------------------------------------------
!***  begin debugging convection
!         print_diag=.false.
!         if(i==imd.and.j==jmd)print_diag=.true.
!***  end debugging convection
!-----------------------------------------------------------------------
          call bmj( &
                   entrain,newall,newswap,newupup,nodeep,print_diag &
                  ,i,j,lm,lmh,ntsd &
                  ,lbot,ltop &
                  ,cldefi(i,j),dtcnvc &
                  ,psfc,ptop,seamask &
                  ,dspbfl,dsp0fl,dsptfl,dspbfs,dsp0fs,dsptfs &
                  ,dspbsl,dsp0sl,dsptsl,dspbss,dsp0ss,dsptss &
                  ,slopbl,slop0l,sloptl,slopbs,slop0s,slopts &
                  ,dpcol,pcol,qcol,tcol,exnercol,ucol,vcol &
                  ,dqdt,dtdt,dudt,dvdt,pcpcol &
                  )
!-----------------------------------------------------------------------
!***  compute momentum tendencies
!
          do k=1,lm
            dudt_phy(i,j,k)=dudt(k)
            dvdt_phy(i,j,k)=dvdt(k)
          enddo
!
!***  compute heating and moistening tendencies
!
          if(present(rthcuten).and.present(rqcuten))then
            do k=1,lm
              rthcuten(i,j,k)=dtdt(k)/exner(i,j,k)
              rqcuten(i,j,k)=dqdt(k)
            enddo
          endif
!
!***  all units in bmj scheme are mks, thus convert precip from meters
!***  to millimeters per step for output.
!
          raincv(i,j)=pcpcol*1.e3/ncnvc
!
!***  convective cloud top and bottom from this call
!
          cutop(i,j)=float(lm+1-ltop)
          cubot(i,j)=float(lm+1-lbot)
!
!-----------------------------------------------------------------------
!***  begin debugging convection
          if(print_diag)then
            delt=0.
            delq=0.
            plyr=0.
            if(lbot>0.and.ltop<lbot)then
              do k=ltop,lbot
                plyr=plyr+dpcol(k)
                delq=delq+dpcol(k)*dtcnvc*abs(dqdt(k))
                delt=delt+dpcol(k)*dtcnvc*abs(dtdt(k))
              enddo
              delq=delq/plyr
              delt=delt/plyr
            endif
!
            write(6,"(2a,2i4,3e12.4,f7.2,4i3)") &
                 '{cu3 i,j,pcpcol,dtavg,dqavg,plyr,'   &
                 ,'lbot,ltop,cubot,cutop = '   &
                 ,i,j, pcpcol,delt,1000.*delq,.01*plyr   &
                 ,lbot,ltop,nint(cubot(i,j)),nint(cutop(i,j))
!
            if(plyr> 0.)then
              do k=lbot,ltop,-1
                write(6,"(a,i3,2e12.4,f7.2)") '{cu3a k,dt,dq,dp = ' &
                     ,k,dtcnvc*dtdt(k),1000.*dtcnvc*dqdt(k),.01*dpcol(k)
              enddo
            endif
          endif
!***  end debugging convection
!-----------------------------------------------------------------------
!
        enddo
        enddo
!.......................................................................
!zj$omp end parallel do
!.......................................................................
!
      endif
!
      end subroutine bmjdrv
!-----------------------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!-----------------------------------------------------------------------
                          subroutine bmj                                 &
!-----------------------------------------------------------------------
        ( &
         entrain,newall,newswap,newupup,nodeep,print_diag &
        ,i,j,lm,lmh,ntsd &
        ,lbot,ltop &
        ,cldefi,dtcnvc &
        ,psfc,pt,sm &
        ,dspbfl,dsp0fl,dsptfl,dspbfs,dsp0fs,dsptfs &
        ,dspbsl,dsp0sl,dsptsl,dspbss,dsp0ss,dsptss &
        ,slopbl,slop0l,sloptl,slopbs,slop0s,slopts &
        ,dprs,prsmid,q,t,exner,u,v &
        ,dqdt,dtdt,dudt,dvdt,pcpcol &
        )
!-----------------------------------------------------------------------
!zj  new shallow cloud added in june 2008 to address swap point
!zj  convection and shallow convection transporting both heat and moisture
!zj  up.  'soft' version of the cloud is implemented.
!zj  reintroduced entrainment on 11/02/2008
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
      logical(kind=klog),intent(in):: &
       entrain,newall,newswap,newupup,nodeep
!
      integer(kind=kint),intent(in):: &
       i,j,lm,lmh,ntsd

      integer(kind=kint),intent(out):: &
       lbot,ltop
!
      real(kind=kfpt),intent(in):: &
       dtcnvc &
      ,psfc,pt,sm &
      ,dspbfl,dsp0fl,dsptfl,dspbfs,dsp0fs,dsptfs &
      ,dspbsl,dsp0sl,dsptsl,dspbss,dsp0ss,dsptss &
      ,slopbl,slop0l,sloptl,slopbs,slop0s,slopts
!
      real(kind=kfpt),dimension(1:lm),intent(in):: &
       dprs,prsmid,q,t,exner,u,v
!
      real(kind=kfpt),intent(inout):: &
       cldefi,pcpcol
!
      real(kind=kfpt),dimension(1:lm),intent(inout):: &
       dqdt,dtdt,dudt,dvdt
!
!-----------------------------------------------------------------------
!***  define local variables
!-----------------------------------------------------------------------
      logical(kind=klog):: &
       deep,mmntdeep,mmntshal1,mmntshal2,plume,shallow,overshoot   !jun04
!
      real(kind=kfpt):: &
       dum,dvm,facuv,uvscald,uvscals1,uvscals2,ubar,vbar

      real(kind=kfpt),dimension(1:lm):: &
       rxnerk,rxnersk,el,fpk &
      ,pk,psk,qbte,qbtk,qk,qrefk,qsatk &
      ,therk,thesp,thevrf,thsk &
      ,thvmod,thvref,tk,trefk,urefk,vrefk,wcld
!
      real(kind=kfpt),dimension(1:lm):: &
       rxner,difq,dift,thee,thes,tref
!
      real(kind=kfpt),dimension(1:lm):: &
       cpe,cpecnv,dtv,dtvcnv,thescnv    !<-- cpe for shallow convection buoyancy check (24 aug 2006)
!
      real(kind=kfpt),dimension(1:lm):: &
       rhk,thmak,thvmk
!
!***  begin debugging convection
      logical(kind=klog) :: print_diag
!***  end debugging convection
!
!-----------------------------------------------------------------------
!***
!***  local scalars
!***
!-----------------------------------------------------------------------
      integer(kind=kint):: &
       iq,iqtb,iswap,it,iter,itrefi,ittb,ittbk &
      ,kb,knumh,knuml &
      ,l,l0,l0m1,lb,lbm1,lcor,lpt1 &
      ,lqm,lshu,ltp1,ltp2,ltsh, lbotcnv,ltopcnv,lmid
!
      real(kind=kfpt):: &
       avrgt,avrgtl,bq,bqk,bqs00k,bqs10k &
      ,cup,den,dentpy,depmin,depth &
      ,depwl,dhdt,difql,diftl,dp,dpkl,dplo,dpmix,dpot &
      ,dpup,dqref,drhdp,drheat,dsp &
      ,dsp0,dsp0k,dspb,dspbk,dspt,dsptk &
      ,dsq,dst,dstq,dthem,dtdp,efi &
      ,fefi,ffup,fprs,fptk,fup,hcorr &
      ,otsum,p,p00k,p01k,p10k,p11k &
      ,part1,part2,part3,pbot,pbotfc,pbtk &
      ,pk0,pkb,pkl,pkt,pklo,pkhi &
      ,plmh,pelevfc,pbtmx,plo,potsum,pp1,ppk,preck &
      ,presk,psp,psum,pthrs,ptop,ptpk,pup &
      ,qbt,qkl,qnew,qotsum,qq1,qqk,qrfkl &
      ,qrftp,qsp,qsum,qup,rdp0t &
      ,rdpsum,rdtcnvc,rhh,rhl,rhmax,rotsum,rtbar,rhavg &
      ,rxnerp,rxnerlo,rxnerhi,rxners,rxnersts &
      ,sm1,smix,sq,sqk,sqs00k,sqs10k,stabdl,sumde,sumdp &
      ,sumdt,tauk,tauksc,tcorr,thbt,therkx,therky &
      ,thskl,thtpk,thvmkl,tkl,tlo,tnew &
      ,tq,tqk,treflo,trfkl,trmlo,trmup,tskl,tsp,tth &
      ,tthk,tup &
      ,capecnv,pspcnv,thbtcnv,capetrigr,cape &
      ,tlev2,qsat1,qsat2,rhshmax,rh_deep      !mar11
!-----------------------------------------------------------------------
!
      real(kind=kfpt),parameter:: &
       elevfc=0.6
!
      real(kind=kfpt),parameter:: &
       slopst=(stabdf-stabds)/(1.-efimn) &
      ,slopel=(1.-efmntl)/(1.-efimn) &
      ,slopes=(1.-efmnts)/(1.-efimn)
!
      real(kind=kfpt),parameter:: &
!       wdry=0.50 &
       wdry=0.75 &
!       wdry=0.90 &
      ,deftop=.95 !soft2
!zj      ,deftop=1.0 ! hard
!
      real(kind=kfpt):: &
       a23m4l,cprlg,elocp,rcp,qwat &
      ,a11,a12,a21,a22,ama,aqs,arh,au,av,avm &
      ,b1qsat,b1rh,b1thma,b1thvm,b1u,b1v &
      ,b2qsat,b2rh,b2thma,b2thvm,b2u,b2v &
      ,bma,bqs,brh,bu,bvm,bv &
      ,qcorr,rden,rhmean,rhref,sumdq,sumdu,sumdv,sumrh &
      ,ucorr,vcorr &
      ,adef,fk
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      cprlg=cp/(row*g*elwv)
      elocp=eliwv/cp
      rcp=1./cp
      a23m4l=a2*(a3-a4)*elwv
!
      uvscald=0.01
      uvscals1=0.01
      uvscals2=0.01
!
      rdtcnvc=1./dtcnvc
      depmin=psh*psfc*rsfcp
!-----------------------------------------------------------------------
      deep=.false.
      plume=.false.
      shallow=.false.
!
!mar11 changes start
      if(nodeep) then
        rh_deep=0.90
      else
        rh_deep=0.75
      endif
!mar11 changes end
!
      mmntdeep=.false. !.true. !.false. !.true.
      mmntshal1=.true. !.false.
      mmntshal2=.true. !.false.
!-----------------------------------------------------------------------
      tauk  =dtcnvc/(trel*01.0)
      tauksc=dtcnvc/(trel*01.0)
!-----------------------------------------------------------------------
!-----------------------------preparations------------------------------
!-----------------------------------------------------------------------
      iswap=0
!
      cup=0.
      dsp0=0.
      dspb=0.
      dspt=0.
      pcpcol=0.
      tref(1)=t(1)
!
      do l=1,lmh
        tk(l)=t(l)
        qk(l)=q(l)
        rxner(l)=1./exner(l)
        cpecnv(l)=0.
        dtvcnv(l)=0.
        thes(l)=0.
        thesp(l)=0.
        thescnv(l)=0.
        dudt(l)=0.
        dvdt(l)=0.
        qsatk(l)=pq0/prsmid(l)*exp(a2*(tk(l)-a3)/(tk(l)-a4))
        rhk(l)=qk(l)/qsatk(l)
      enddo
!-----------------------------------------------------------------------
!----------------search for maximum buoyancy level----------------------
!-----------------------------------------------------------------------
      pelevfc=prsmid(lmh)*elevfc
      pbtmx  =prsmid(lmh)-pone
      capecnv=0.
      pspcnv =0.
      thbtcnv=0.
      lbot=lmh
      ltop=lmh
      lbotcnv=lbot
      ltopcnv=lbot
!-----------------------------------------------------------------------
!----------------trial maximum buoyancy level variables-----------------
!-----------------------------------------------------------------------
      prep_loop: do kb=lmh,1,-1
!-----------------------------------------------------------------------
        if(prsmid(kb).lt.pelevfc.and..not.entrain) exit
!---preparation for search for max cape---------------------------------
        qbt=q(kb)
        thbt=t(kb)*rxner(kb)
        tth=(thbt-thl)*rdth
        qq1=tth-aint(tth)
        ittb=int(tth)+1
!---keeping indices within the table------------------------------------
        if(ittb.lt.1)then
          ittb=1
          qq1=0.
        else if(ittb.ge.jtb)then
          ittb=jtb-1
          qq1=0.
        endif
!---base and scaling factor for spec. humidity--------------------------
        ittbk=ittb
        bqs00k=qs0(ittbk)
        sqs00k=sqs(ittbk)
        bqs10k=qs0(ittbk+1)
        sqs10k=sqs(ittbk+1)
!--------------scaling spec. humidity & table index---------------------
        bq=(bqs10k-bqs00k)*qq1+bqs00k
        sq=(sqs10k-sqs00k)*qq1+sqs00k
        tq=(qbt-bq)/sq*rdq
        pp1=tq-aint(tq)
        iqtb=int(tq)+1
!----------------keeping indices within the table-----------------------
        if(iqtb.lt.1)then
          iqtb=1
          pp1=0.
        else if(iqtb.ge.itb)then
          iqtb=itb-1
          pp1=0.
        endif
!--------------saturation pressure at four surrounding table pts.-------
        iq=iqtb
        it=ittb
        p00k=ptbl(iq  ,it  )
        p10k=ptbl(iq+1,it  )
        p01k=ptbl(iq  ,it+1)
        p11k=ptbl(iq+1,it+1)
!--------------saturation point variables at the bottom-----------------
        psp=p00k+(p10k-p00k)*pp1+(p01k-p00k)*qq1 &
      &     +(p00k-p10k-p01k+p11k)*pp1*qq1
        rxners=(1.e5/psp)**cappa
        thesp  (kb)=thbt*exp(elocp*qbt*rxners/thbt)
        psk    (kb)=psp
        rxnersk(kb)=rxners
!-----------------------------------------------------------------------
      enddo prep_loop

!write(0,*)'thesp',thesp
!write(0,*)'psk',psk
!write(0,*)'rxnersk',rxnersk

!-----------------------------------------------------------------------
      max_buoy_loop: do kb=lmh,1,-1
!---choose cloud base as model level just below psp---------------------
!-----------------------------------------------------------------------
        if(prsmid(kb).lt.pelevfc) exit
!---search over a scaled depth to find the parcel with the max cape-----
        qbt=q(kb)
        thbt=t(kb)*rxner(kb)
        psp=psk(kb)
        rxners=rxnersk(kb)
!
        do l=1,lmh-1
          p=prsmid(l)
          if(p.lt.psp.and.p.ge.pqm) lbot=l+1
        enddo
!***
!*** warning: lbot must not be .gt. lm-1 in shallow convection
!*** make sure cloud base is at least pone above the surface
!***
        pbot=prsmid(lbot)
        if(pbot.ge.pbtmx.or.lbot.ge.lmh)then
          do l=1,lmh-1
            p=prsmid(l)
            if(p.lt.pbtmx)lbot=l
          enddo
          pbot=prsmid(lbot)
        endif
!---cloud top computation-----------------------------------------------
        ltop=lbot
        ptop=pbot
!---entrainment during parcel ascent------------------------------------
        if(entrain) then
!-----------------------------------------------------------------------
          do l=lmh,kb,-1
            thes(l)=thesp(kb)
            qbtk(l)=qk   (kb)
          enddo
          do l=kb,1,-1
            thes(l)=thesp(l)
            qbtk(l)=qk   (l)
          enddo
!
          do l=1,lmh
            thee(l)=thes(l)
            qbte(l)=qbtk(l)
          enddo
!
          if(sm.gt.0.5) then
            fefi=(cldefi-efimn)*slopes+efmnts
          else
            fefi=(cldefi-efimn)*slopel+efmntl
          endif
!
          ffup=fup/(fefi*fefi)
!
          if(pbot.gt.enplo)then
            fprs=1.
          elseif(pbot.gt.enpup)then
            fprs=(pbot-enpup)*rendp
          else
            fprs=0.
          endif
!
          ffup=ffup*fprs*fprs*0.5
          dpup=dprs(kb)
!
          do l=kb-1,1,-1
            dplo=dpup
            dpup=dprs(l)
!
            thes(l)=((-ffup*dplo+1.)*thes(l+1) &
                    +(thee(l)*dpup+thee(l+1)*dplo)*ffup) &
                   /(ffup*dpup+1.)
            qbtk(l)=((-ffup*dplo+1.)*qbtk(l+1) &
                    +(qbte(l)*dpup+qbte(l+1)*dplo)*ffup) &
                   /(ffup*dpup+1.)
          enddo
!---no entrainment------------------------------------------------------
        else
!-----------------------------------------------------------------------
          do l=lmh,1,-1
            thes(l)=thesp(kb)
            qbtk(l)=qk   (kb)
          enddo
!-----------------------------------------------------------------------
        endif ! end of entrainment/no entrainment
!------------------first entropy check----------------------------------
        do l=1,lmh
          cpe(l)=0.
          dtv(l)=0.
        enddo
!-----------------------------------------------------------------------
!-- begin: buoyancy check including deep convection (24 aug 2006)
!-----------------------------------------------------------------------
        if(kb.lt.lbot) go to 170
!-----------------------------------------------------------------------
        l=kb
        plo=prsmid(l)
        tlo=thbt*exner(l)
        trmlo=((qbt*p608+1.)*tlo-(q(l)*p608+1.)*t(l))*r_d/plo
        capetrigr=dtptrigr/t(lbot)
!
!--- below cloud base
!
        if(kb-lbot.ge.1) then
          do l=kb-1,lbot,-1
            pup=prsmid(l)
            tup=thbt*exner(l)
            dp=plo-pup
            trmup=((qbt*p608+1.)*tup-(q(l)*p608+1.)*t(l))*r_d/pup
!
            dtv(l)=(trmlo+trmup)*dp*0.5
            cpe(l)=dtv(l)+cpe(l+1)
!
            if(cpe(l).lt.capetrigr) go to 170
!
            plo=pup
            trmlo=trmup
          enddo
        endif
!
!--- cloud base layer
!
        l=lbot-1
        pup=psp
        tup=thbt/rxners
        tsp=((t(l+1)-t(l))/(plo-prsmid(l)))*(pup-pbot)+t(l)
        qsp=((q(l+1)-q(l))/(plo-prsmid(l)))*(pup-pbot)+q(l)
        dp=plo-pup
        trmup=((qbt*p608+1.)*tup-(qsp*p608+1.)*tsp)*r_d/pup
!
        dtv(l)=(trmlo+trmup)*dp*0.5
!
        plo=pup
        trmlo=trmup
!
!--- above saturation pressure updraft t and q along moist adiabat
!
        pup=prsmid(l)
!
!--- calculate updraft temperature along moist adiabat (tup)
!
        if(pup<plq)then
          call ttblex(itb,jtb,pl,pup,rdp,rdthe &
                     ,sthe,the0,thes(l),ttbl,tup)
        else
          call ttblex(itbq,jtbq,plq,pup,rdpq,rdtheq &
                     ,stheq,the0q,thes(l),ttblq,tup)
        endif
!
        qup=pq0/pup*exp(a2*(tup-a3)/(tup-a4))
        qwat=qbt-qup  !-- water loading effects, reversible adiabat
        dp=plo-pup
        trmup=((qup*p608+1.-qwat)*tup-(q(l)*p608+1.)*t(l))*r_d/pup
!
        dtv(l)=(trmlo+trmup)*dp*0.5+dtv(l)
        cpe(l)=dtv(l)+cpe(l+1)
!
        if(cpe(l).lt.capetrigr) go to 170
!
        plo=pup
        trmlo=trmup
!
!-----------------------------------------------------------------------
!--- in cloud above cloud base
!-----------------------------------------------------------------------
!
        do l=lbot-2,1,-1
          pup=prsmid(l)
!
!--- calculate updraft temperature along moist adiabat (tup)
!
          if(pup<plq)then
            call ttblex(itb,jtb,pl,pup,rdp,rdthe &
                       ,sthe,the0,thes(l),ttbl,tup)
          else
            call ttblex(itbq,jtbq,plq,pup,rdpq,rdtheq &
                       ,stheq,the0q,thes(l),ttblq,tup)
          endif
!
          qup=pq0/pup*exp(a2*(tup-a3)/(tup-a4))
          qwat=qbt-qup  !-- water loading effects, reversible adiabat
          dp=plo-pup
          trmup=((qup*p608+1.-qwat)*tup-(q(l)*p608+1.)*t(l))*r_d/pup
!
          dtv(l)=(trmlo+trmup)*dp*0.5
          cpe(l)=dtv(l)+cpe(l+1)
!
          if(cpe(l).lt.capetrigr) go to 170
!
          plo=pup
          trmlo=trmup
        enddo
!
!-----------------------------------------------------------------------
!
170     ltp1=kb
        cape=0.
!
!-----------------------------------------------------------------------
!--- cloud top level (ltop) is located where cape is a maximum
!--- exit cloud-top search if cape < capetrigr
!-   Also exit cloud-top search if newswap=T and RH<rh_deep  !jun04
!-----------------------------------------------------------------------
!
          do l=kb,1,-1
            if (cpe(l) < capetrigr) then
              exit
            else if (cpe(l) > cape) then
              ltp1=l
              cape=cpe(l)
            endif
            if(newswap .and. rhk(l)<rh_deep) exit   !jun04
          enddo      !-- end do l=kb,1,-1
!
          ltop=max(min(ltp1,lbot),1)
!
!-----------------------------------------------------------------------
!--------------- check for maximum instability  ------------------------
!-----------------------------------------------------------------------
          if(cape > capecnv) then
            capecnv=cape
            pspcnv=psp
            thbtcnv=thbt
            lbotcnv=lbot
            ltopcnv=ltop
            do l=lmh,1,-1
              cpecnv(l)=cpe(l)
              dtvcnv(l)=dtv(l)
              thescnv(l)=thes(l)
            enddo
          endif    ! end if(cape > capecnv) then
!
!-----------------------------------------------------------------------
!
      enddo max_buoy_loop
!
!-----------------------------------------------------------------------
!------------------------  maximum instability  ------------------------
!-----------------------------------------------------------------------
!
      if(capecnv > 0.) then
        psp=pspcnv
        thbt=thbtcnv
        lbot=lbotcnv
        ltop=ltopcnv
        pbot=prsmid(lbot)
        ptop=prsmid(ltop)
!
        do l=lmh,1,-1
          cpe(l)=cpecnv(l)
          dtv(l)=dtvcnv(l)
          thes(l)=thescnv(l)
        enddo
!
      endif
!
!-----------------------------------------------------------------------
!-----  quick exit if cloud is too thin or no cape is present  ---------
!-----------------------------------------------------------------------
!
      if(ptop>pbot-pno.or.ltop>lbot-2.or.capecnv<=0.)then
        lbot=0
        ltop=lm
        pbot=prsmid(lmh)
        ptop=pbot
        cldefi=avgefi*sm+stefi*(1.-sm)
        return
      endif
!
!***  depth of cloud required to make the point a deep convection point
!***  is a scaled value of psfc.
!
      depth=pbot-ptop
!
!zj      if(depth.ge.depmin*0.50) then
!zj      if(depth.ge.depmin*0.25) then
!zj      if(depth.ge.depmin*0.625) then
      if(depth.ge.depmin*0.75) then
!zj      if(depth.ge.depmin*0.85) then
        plume=.true.
      endif
!
      if(depth>=depmin) then
        deep=.true.
      else
        shallow=.true.
        go to 600
      endif
!
!-----------------------------------------------------------------------
!dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
!dcdcdcdcdcdcdcdcdcdcdc    deep convection   dcdcdcdcdcdcdcdcdcdcdcdcdcd
!dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
!-----------------------------------------------------------------------
!
  300 continue
!
      lb =lbot
      efi=cldefi
!-----------------------------------------------------------------------
!--------------initialize variables in the convective column------------
!-----------------------------------------------------------------------
!***
!***  one should note that the values assigned to the array trefk
!***  in the following loop are really only relevant in anchoring the
!***  reference temperature profile at level lb.  when building the
!***  reference profile from cloud base, then assigning the
!***  ambient temperature to trefk is acceptable.  however, when
!***  building the reference profile from some other level (such as
!***  one level above the ground), then trefk should be filled with
!***  the temperatures in tref(l) which are the temperatures of
!***  the moist adiabat through cloud base.  by the time the line
!***  numbered 450 has been reached, trefk actually does hold the
!***  reference temperature profile.
!***
      do l=1,lmh
        dift(l)=0.
        difq(l)=0.
        tkl=t(l)
        tk(l)=tkl
        trefk(l)=tkl
        qkl=q(l)
        qk(l)=qkl
        qrefk(l)=qkl
        pkl=prsmid(l)
        pk(l)=pkl
        psk(l)=pkl
        rxnerp=rxner(l)
        rxnerk(l)=rxnerp
!
!--- calculate temperature along moist adiabat (tref)
!
        if(pkl<plq)then
          call ttblex(itb,jtb,pl,pkl,rdp,rdthe                           &
      &               ,sthe,the0,thes(l),ttbl,tref(l))
        else
          call ttblex(itbq,jtbq,plq,pkl,rdpq,rdtheq                      &
      &               ,stheq,the0q,thes(l),ttblq,tref(l))
        endif
        therk (l)=tref(l)*rxnerp
      enddo
!
!------------deep convection reference temperature profile------------
      ltp1=ltop+1
      lbm1=lb-1
      pkb=pk(lb)
      pkt=pk(ltop)
      stabdl=(efi-efimn)*slopst+stabds
!------------temperature reference profile below freezing level-------
      el(lb) = elwv
      l0=lb
      pk0=pk(lb)
      treflo=trefk(lb)
      therkx=therk(lb)
      rxnerlo=rxnerk(lb)
      therky=therk(lbm1)
      rxnerhi=rxnerk(lbm1)
!
      do l=lbm1,ltop,-1
        if(t(l+1)<tiw)go to 430
        treflo=((therky-therkx)*stabdl &
               +treflo*rxnerlo)/rxnerhi
        trefk(l)=treflo
        el(l)=elwv
        rxnerlo=rxnerhi
        therkx=therky
        rxnerhi=rxnerk(l-1)
        therky=therk(l-1)
        l0=l
        pk0=pk(l0)
      enddo
!--------------freezing level at or above the cloud top-----------------
      go to 450
!--------------temperature reference profile above freezing level-------
  430 l0m1=l0-1
      rdp0t=1./(pk0-pkt)
      dthem=therk(l0)-trefk(l0)*rxnerk(l0)
!
      do l=ltop,l0m1
        trefk(l)=(therk(l)-(pk(l)-pkt)*dthem*rdp0t)/rxnerk(l)
!        el(l)=elwv
        el(l)=eliv
      enddo
!
!-----------------------------------------------------------------------
!--------------deep convection reference humidity profile---------------
!-----------------------------------------------------------------------
!
!***  depwl is the pressure difference between cloud base and
!***  the freezing level
!
  450 continue
      depwl=pkb-pk0
      depth=pfrz*psfc*rsfcp
      sm1=1.-sm
      pbotfc=1.
!
!-------------first adjustment of temperature profile-------------------
!!
!      sumdt=0.
!      sumdp=0.
!!
!      do l=ltop,lb
!        sumdt=(tk(l)-trefk(l))*dprs(l)+sumdt
!        sumdp=sumdp+dprs(l)
!      enddo
!!
!      tcorr=sumdt/sumdp
!!
!      do l=ltop,lb
!        trefk(l)=trefk(l)+tcorr
!      enddo
!!
!-----------------------------------------------------------------------
!--------------- iteration loop for cloud efficiency -------------------
!-----------------------------------------------------------------------
!
      cloud_efficiency : do itrefi=1,itrefi_max
!
!-----------------------------------------------------------------------
        dspbk=((efi-efimn)*slopbs+dspbss*pbotfc)*sm &
             +((efi-efimn)*slopbl+dspbsl*pbotfc)*sm1
        dsp0k=((efi-efimn)*slop0s+dsp0ss*pbotfc)*sm &
             +((efi-efimn)*slop0l+dsp0sl*pbotfc)*sm1
        dsptk=((efi-efimn)*slopts+dsptss*pbotfc)*sm &
             +((efi-efimn)*sloptl+dsptsl*pbotfc)*sm1
!
!-----------------------------------------------------------------------
!
        do l=ltop,lb
!
!***
!***  saturation pressure difference
!***
          if(depwl>=depth)then
            if(l<l0)then
              dsp=((pk0-pk(l))*dsptk+(pk(l)-pkt)*dsp0k)/(pk0-pkt)
            else
              dsp=((pkb-pk(l))*dsp0k+(pk(l)-pk0)*dspbk)/(pkb-pk0)
            endif
          else
            dsp=dsp0k
            if(l<l0)then
              dsp=((pk0-pk(l))*dsptk+(pk(l)-pkt)*dsp0k)/(pk0-pkt)
            endif
          endif
!***
!***  humidity profile
!***
          if(pk(l)>pqm)then
            psk(l)=pk(l)+dsp
            rxnersk(l)=(1.e5/psk(l))**cappa
            thsk(l)=trefk(l)*rxnerk(l)
            qrefk(l)=pq0/psk(l)*exp(a2*(thsk(l)-a3*rxnersk(l)) &
                                      /(thsk(l)-a4*rxnersk(l)))
          else
            qrefk(l)=qk(l)
          endif
!
        enddo
!-----------------------------------------------------------------------
!***
!***  enthalpy conservation integral
!***
!-----------------------------------------------------------------------
        enthalpy_conservation : do iter=1,2
!
          sumde=0.
          sumdp=0.
!
          do l=ltop,lb
            sumde=((tk(l)-trefk(l))*cp+(qk(l)-qrefk(l))*el(l))*dprs(l)   &
      &            +sumde
            sumdp=sumdp+dprs(l)
          enddo
!
          hcorr=sumde/(sumdp-dprs(ltop))
          lcor=ltop+1
!***
!***  find lqm
!***
          lqm=1
          do l=1,lb
            if(pk(l)<=pqm)lqm=l
          enddo
!***
!***  above lqm correct temperature only
!***
          if(lcor<=lqm)then
            do l=lcor,lqm
              trefk(l)=trefk(l)+hcorr*rcp
            enddo
            lcor=lqm+1
          endif
!***
!***  below lqm correct both temperature and moisture
!***
          do l=lcor,lb
            tskl=trefk(l)*rxnerk(l)/rxnersk(l)
            dhdt=qrefk(l)*a23m4l/(tskl-a4)**2+cp
            trefk(l)=hcorr/dhdt+trefk(l)
            thskl=trefk(l)*rxnerk(l)
            qrefk(l)=pq0/psk(l)*exp(a2*(thskl-a3*rxnersk(l)) &
                                      /(thskl-a4*rxnersk(l)))
          enddo
!
        enddo  enthalpy_conservation
!-----------------------------------------------------------------------
!
!***  heating, moistening, precipitation
!
!-----------------------------------------------------------------------
        avrgt=0.
        preck=0.
        dsq=0.
        dst=0.
!
        do l=ltop,lb
          tkl=tk(l)
          diftl=(trefk(l)-tkl  )*tauk
          difql=(qrefk(l)-qk(l))*tauk
          avrgtl=(tkl+tkl+diftl)
          dpot=dprs(l)/avrgtl
          dst=diftl*dpot+dst
          dsq=difql*el(l)*dpot+dsq
          avrgt=avrgtl*dprs(l)+avrgt
          preck=diftl*dprs(l)+preck
          dift(l)=diftl
          difq(l)=difql
        enddo
!
        dst=(dst+dst)*cp
        dsq=dsq+dsq
        dentpy=dst+dsq
        avrgt=avrgt/(sumdp+sumdp)
!
!        drheat=preck*cp/avrgt
        drheat=(preck*sm+max(1.e-7,preck)*(1.-sm))*cp/avrgt
        drheat=max(drheat,1.e-20)
        efi=efifc*dentpy/drheat
!-----------------------------------------------------------------------
        efi=min(efi,1.)
        efi=max(efi,efimn)
!-----------------------------------------------------------------------
!
      enddo  cloud_efficiency
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!---------------------- deep convection --------------------------------
!-----------------------------------------------------------------------
      if(dentpy>=epsntp.and.preck>epspr.and..not.nodeep) then
!
        iswap=0 ! deep convection, no swap
        cldefi=efi
!
        if(sm.gt.0.5) then
          fefi=(cldefi-efimn)*slopes+efmnts
        else
          fefi=(cldefi-efimn)*slopel+efmntl
        endif
!
        fefi=(dentpy-epsntp)*fefi/dentpy
        preck=preck*fefi
!
!***  update precipitation and tendencies of temperature and moisture
!
        cup=preck*cprlg
        pcpcol=cup
!
        do l=ltop,lb
          dtdt(l)=dift(l)*fefi*rdtcnvc
          dqdt(l)=difq(l)*fefi*rdtcnvc
        enddo
!-----------------------------------------------------------------------
        if(mmntdeep) then
          facuv=fefi*rdtcnvc*uvscald
          if(l0.gt.ltop.and.l0.lt.lb) then
            ubar=0.
            vbar=0.
            sumdp=0.
!
            do l=l0,lb
              ubar=u(l)*dprs(l)+ubar
              vbar=v(l)*dprs(l)+vbar
              sumdp=dprs(l)+sumdp
            enddo
!
            rdpsum=1./sumdp
            ubar=ubar*rdpsum
            vbar=vbar*rdpsum
!
            do l=l0,lb
              dudt(l)=(ubar-u(l))*facuv
              dvdt(l)=(vbar-v(l))*facuv
            enddo
!
            dum=ubar-u(l0)
            dvm=vbar-v(l0)
!
            do l=ltop,l0-1
              dudt(l)=(pk(l)-pkt)*dum*rdp0t*facuv
              dvdt(l)=(pk(l)-pkt)*dvm*rdp0t*facuv
            enddo
          else
            ubar=0.
            vbar=0.
            sumdp=0.
!
            do l=ltop,lb
              ubar=u(l)*dprs(l)+ubar
              vbar=v(l)*dprs(l)+vbar
              sumdp=dprs(l)+sumdp
            enddo
!
            rdpsum=1./sumdp
            ubar=ubar*rdpsum
            vbar=vbar*rdpsum
!
            do l=ltop,lb
              dudt(l)=(ubar-u(l))*facuv
              dvdt(l)=(vbar-v(l))*facuv
            enddo
          endif
        endif
!-----------------------------------------------------------------------
      else
!-----------------------------------------------------------------------
!***  reduce the cloud top
!-----------------------------------------------------------------------
!
!        ltop=ltop+3           !iterate cloud top
!        ptop=prsmid(ltop)     !iterate cloud top
!        depmin=psh*psfc*rsfcp !iterate cloud top
!        depth=pbot-ptop       !iterate cloud top
!***
!***  iterate deep convection procedure if needed
!***
!        if(depth>=depmin)then !iterate cloud top
!          go to 300           !iterate cloud top
!        endif                 !iterate cloud top
!
!         cldefi=avgefi
         cldefi=efimn*sm+stefi*(1.-sm)
!***
!***  search for shallow cloud top
!***
!        ltsh=lbot
!        lbm1=lbot-1
!        pbtk=pk(lbot)
!        depmin=psh*psfc*rsfcp
!        ptpk=pbtk-depmin
        ptpk=max(pshu, pk(lbot)-depmin)
!***
!***  cloud top is the level just below pbtk-psh or just below pshu
!***
      if(.not.newswap) then   !jun04
        do l=1,lmh
          if(pk(l)<=ptpk)ltop=l+1
        enddo
      endif    !jun04
!
!        ptpk=pk(ltop)
!!***
!!***  highest level allowed is level just below pshu
!!***
!        if(ptpk<=pshu)then
!!
!          do l=1,lmh
!            if(pk(l)<=pshu)lshu=l+1
!          enddo
!!
!          ltop=lshu
!          ptpk=pk(ltop)
!        endif
!
!        if(ltop>=lbot)then
!!!!!!     lbot=0
!          ltop=lmh
!!!!!!     pbot=pk(lbot)
!          ptop=pk(ltop)
!          pbot=ptop
!          go to 600
!        endif
!
!        ltp1=ltop+1
!        ltp2=ltop+2
!!
!        do l=ltop,lbot
!          qsatk(l)=pq0/pk(l)*exp(a2*(tk(l)-a3)/(tk(l)-a4))
!        enddo
!!
!        rhh=qk(ltop)/qsatk(ltop)
!        rhmax=0.
!        ltsh=ltop
!!
!        do l=ltp1,lbm1
!          rhl=qk(l)/qsatk(l)
!          drhdp=(rhh-rhl)/(pk(l-1)-pk(l))
!!
!          if(drhdp>rhmax)then
!            ltsh=l-1
!            rhmax=drhdp
!          endif
!!
!          rhh=rhl
!        enddo
!
!-----------------------------------------------------------------------
!-- make shallow cloud top a function of virtual temperature excess (dtv)
!-----------------------------------------------------------------------
!
      if(.not.newswap) then   !jun04
        ltp1=lbot
        do l=lbot-1,ltop,-1
          if (dtv(l) > 0.) then
            ltp1=l
          else
            exit
          endif
        enddo
        ltop=min(ltp1,lbot)
      endif   !jun04
!***
!***  cloud must be at least two layers thick
!***
!        if(lbot-ltop<2)ltop=lbot-2  (eliminate this criterion)
!
!-- end: buoyancy check (24 aug 2006)
!
        iswap=1 ! failed deep convection, shallow swap point
        ptop=pk(ltop)
        shallow=.true.
        deep=.false.
!
      endif
!dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
!dcdcdcdcdcdcdc          end of deep convection            dcdcdcdcdcdcd
!dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  600 continue
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!----------------gather shallow convection points-----------------------
!
!      if(ptop<=pbot-pno.and.ltop<=lbot-2)then
!         depmin=psh*psfc*rsfcp
!!
!!        if(lpbl<lbot)lbot=lpbl
!!        if(lbot>lmh-1)lbot=lmh-1
!!        pbot=prsmid(lbot)
!!
!         if(ptop+1.>=pbot-depmin)shallow=.true.
!      else
!         lbot=0
!         ltop=lm
!      endif
!
!***********************************************************************
!-----------------------------------------------------------------------
!***  begin debugging convection
      if(print_diag)then
        write(6,"(a,2i3,l2,3e12.4)")   &
             '{cu2a lbot,ltop,shallow,pbot,ptop,depmin = ' &
             ,lbot,ltop,shallow,pbot,ptop,depmin
      endif
!***  end debugging convection
!-----------------------------------------------------------------------
!
      if(.not.shallow)return
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!scscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
!scscscscscscsc         shallow convection          cscscscscscscscscscs
!scscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
!-----------------------------------------------------------------------
      do l=1,lmh
        pk(l)=prsmid(l)
        tk(l)=t(l)
        qk(l)=q(l)
        trefk(l)=t(l)
        qrefk(l)=q(l)
        qsatk(l)=pq0/pk(l)*exp(a2*(tk(l)-a3)/(tk(l)-a4))
        rxnerk(l)=rxner(l)
        thvref(l)=tk(l)*rxnerk(l)*(qk(l)*p608+1.)
!
        if(tk(l)>=tiw)then
          el(l)=elwv
        else
          el(l)=eliv
        endif
      enddo
!
!-----------------------------------------------------------------------
!-- begin: raise cloud top if avg rh>rhshmax and cape>0 (dtv>0 if newswap=T)   !jun04
!   rhshmax=rh at cloud base associated with a dsp of pone
!-----------------------------------------------------------------------
!
      tlev2=t(lbot)*((pk(lbot)-pone)/pk(lbot))**cappa
      qsat1=pq0/pk(lbot)*exp(a2*(t(lbot)-a3)/(tk(lbot)-a4))
      qsat2=pq0/(pk(lbot)-pone)*exp(a2*(tlev2-a3)/(tlev2-a4))
      rhshmax=qsat2/qsat1
      sumdp=0.
      rhavg=0.
!
      do l=lbot,ltop,-1
        rhavg=rhavg+dprs(l)*qk(l)/qsatk(l)
        sumdp=sumdp+dprs(l)
      enddo
!
      if (rhavg/sumdp > rhshmax) then
        ltsh=ltop
!-- overshoot: allow cloud top to be 1 level above neutral/positive buoyancy   !apr21
        if(newswap) overshoot=.true.   !jun04
        do l=ltop-1,1,-1
          rhavg=rhavg+dprs(l)*qk(l)/qsatk(l)
          sumdp=sumdp+dprs(l)
          if(.not.newswap) then   !jun04
            if (cpe(l) > 0.) then
              ltsh=l
            else
              exit
            endif
          else    !jun04 begin
            if (dtv(l) <= 0.) then   !-- More strict positive buoyancy criterion
              if (cpe(l) <= 0.) exit
              if (overshoot) then
                overshoot=.false.
              else
                exit
              endif
            endif
            ltsh=l   !jun04 end
          endif
          if (rhavg/sumdp <= rhshmax) exit
          if (pk(l) <= pshu) exit
        enddo
        ltop=ltsh
!swapnomoist        iswap=0 ! old cloud for moist clouds
      endif
!
!-- end: raise cloud top if avg rh>rhshmax and cape>0 (dtv>0 if newswap=T)   !jun04
!
!---------------------------shallow cloud top---------------------------
      lbm1=lbot-1
      ptpk=ptop
      ltp1=ltop-1
      depth=pbot-ptop
!-----------------------------------------------------------------------
!zj      if(depth.ge.depmin*0.50) then
!zj      if(depth.ge.depmin*0.25) then
!zj      if(depth.ge.depmin*0.625) then
!      if(depth.ge.depmin*0.75) then
!zj      if(depth.ge.depmin*0.85) then
!        plume=.true.
!      endif
!-----------------------------------------------------------------------
!***  begin debugging convection
      if(print_diag)then
        write(6,"(a,4e12.4)") '{cu2b pbot,ptop,depth,depmin= ' &
             ,pbot,ptop,depth,depmin
      endif
!***  end debugging convection
!-----------------------------------------------------------------------
!
!bsf      if(depth<depmin)then
!bsf        return
!bsf      endif
!-----------------------------------------------------------------------
      if(ptop>pbot-pno.or.ltop>lbot-2)then
        lbot=0
        ltop=lm
        ptop=pbot
        return
      endif
!-----------------------------------------------------------------------
!***  new cloud at all shallow points
!-----------------------------------------------------------------------
!zj      if(newall) go to 810 ! new cloud at all shallow points

!zj      if(newall.and.sm.lt.0.5) go to 810 ! new cloud at land points
      if(newall.and.plume) go to 810 ! new cloud at plume points
!zj      if(newall.and.plume.and.sm.lt.0.5) go to 810 ! new cloud at plume land points
!-----------------------------------------------------------------------
!***  new cloud at swap shallow points
!-----------------------------------------------------------------------
!zj      if(newswap.and.iswap.gt.0) go to 810 ! new cloud only at swap pts.

!zj      if(newswap.and.iswap.gt.0.and.sm.lt.0.5) go to 810 ! new cloud only at swap pts.
!zj      if(newswap.and.iswap.gt.0.and.plume) go to 810 ! new cloud if plume at swap pts.
!zj      if(newswap.and.iswap.gt.0.and.plume.and.sm.lt.0.5) go to 810 ! new cloud only at swap pts.
!-----------------------------------------------------------------------
!
!--------------scaling potential temperature & table index at top-------
!
      thtpk=t(ltp1)*rxner(ltp1)
!
      tthk=(thtpk-thl)*rdth
      qqk =tthk-aint(tthk)
      it  =int(tthk)+1
!
      if(it<1)then
        it=1
        qqk=0.
      endif
!
      if(it>=jtb)then
        it=jtb-1
        qqk=0.
      endif
!
!--------------base and scaling factor for spec. humidity at top--------
!
      bqs00k=qs0(it)
      sqs00k=sqs(it)
      bqs10k=qs0(it+1)
      sqs10k=sqs(it+1)
!
!--------------scaling spec. humidity & table index at top--------------
!
      bqk=(bqs10k-bqs00k)*qqk+bqs00k
      sqk=(sqs10k-sqs00k)*qqk+sqs00k
!
!     tqk=(q(ltop)-bqk)/sqk*rdq
      tqk=(q(ltp1)-bqk)/sqk*rdq
!
      ppk=tqk-aint(tqk)
      iq =int(tqk)+1
!
      if(iq<1)then
        iq=1
        ppk=0.
      endif
!
      if(iq>=itb)then
        iq=itb-1
        ppk=0.
      endif
!
!----------------cloud top saturation point pressure--------------------
      part1=(ptbl(iq+1,it)-ptbl(iq,it))*ppk
      part2=(ptbl(iq,it+1)-ptbl(iq,it))*qqk
      part3=(ptbl(iq  ,it  )-ptbl(iq+1,it  ) &
            -ptbl(iq  ,it+1)+ptbl(iq+1,it+1))*ppk*qqk
      ptpk=ptbl(iq,it)+part1+part2+part3
!-----------------------------------------------------------------------
      dpmix=ptpk-psp
      if(abs(dpmix).lt.3000.)dpmix=-3000.
!----------------temperature profile slope------------------------------
      smix=(thtpk-thbt)/dpmix*stabs
!
      treflo=trefk(lbot+1)
      pklo=pk(lbot+1)
      pkhi=pk(lbot)
      rxnerlo=rxnerk(lbot+1)
      rxnerhi=rxnerk(lbot)
!
      lmid=.5*(lbot+ltop)
!
      do l=lbot,ltop,-1
        treflo=((pkhi-pklo)*smix+treflo*rxnerlo)/rxnerhi
        trefk(l)=treflo
        if(l<=lmid) trefk(l)=max(trefk(l),tk(l)+dtshal)
        rxnerlo=rxnerhi
        pklo=pkhi
        rxnerhi=rxnerk(l-1)
        pkhi=pk(l-1)
      enddo
!----------------temperature reference profile correction---------------
      sumdt=0.
      sumdp=0.
!
      do l=ltop,lbot
        sumdt=(tk(l)-trefk(l))*dprs(l)+sumdt
        sumdp=sumdp+dprs(l)
      enddo
!
      rdpsum=1./sumdp
      fpk(lbot)=trefk(lbot)
!
      tcorr=sumdt*rdpsum
!
      do l=ltop,lbot
        trfkl   =trefk(l)+tcorr
        trefk(l)=trfkl
        fpk  (l)=trfkl
      enddo
!----------------humidity profile equations-----------------------------
      psum  =0.
      qsum  =0.
      potsum=0.
      qotsum=0.
      otsum =0.
      dst   =0.
      fptk  =fpk(ltop)
!
      do l=ltop,lbot
        dpkl  =fpk(l)-fptk
        psum  =dpkl *dprs(l)+psum
        qsum  =qk(l)*dprs(l)+qsum
        rtbar =2./(trefk(l)+tk(l))
        otsum =dprs(l)*rtbar+otsum
        potsum=dpkl   *rtbar*dprs(l)+potsum
        qotsum=qk(l)  *rtbar*dprs(l)+qotsum
        dst   =(trefk(l)-tk(l))*rtbar*dprs(l)/el(l)+dst
      enddo
!
      psum  =psum*rdpsum
      qsum  =qsum*rdpsum
      rotsum=1./otsum
      potsum=potsum*rotsum
      qotsum=qotsum*rotsum
      dst   =dst*rotsum*cp
!
!-----------------------------------------------------------------------
!***  begin debugging convection
      if(print_diag)then
        write(6,"(a,5e12.4)") '{cu2c dst,psum,qsum,potsum,qotsum = ' &
             ,dst,psum,qsum,potsum,qotsum
      endif
!***  end debugging convection
!-----------------------------------------------------------------------
!***  if upward transport of temperature go to new cloud
!-----------------------------------------------------------------------
!zj      if(newupup.and.dst.gt.0.) go to 810 ! new shallow cloud for both heat and moisture up

!zj      if(newupup.and.dst.gt.0..and.sm.lt.0.5) go to 810 ! new shallow cloud for both heat and moisture up
       if(newupup.and.dst.gt.0..and.plume) go to 810 ! new shallow cloud if plume for both heat and moisture up
!zj      if(newupup.and.dst.gt.0..and.plume.and.sm.lt.0.5) go to 810 ! new shallow cloud for both heat and moisture up
!-----------------------------------------------------------------------
!*** otherwise old cloud
!-----------------------------------------------------------------------
      if(dst.gt.0.) then 
        if (newswap) go to 810 ! new shallow cloud for both heat and moisture up   !jun04
        lbot=0          
        ltop=lm     
        ptop=pbot   
        return 
      endif
!-----------------------------------------------------------------------
!***  otherwise continue with old cloud
!----------------ensure positive entropy change-------------------------
      dstq=dst*epsdn
!----------------check for isothermal atmosphere------------------------
      den=potsum-psum
!
      if(-den/psum<5.e-5)then
        if (newswap) go to 810 ! new shallow cloud for both heat and moisture up   !jun04
        lbot=0
        ltop=lm
        ptop=pbot
        return
!----------------slope of the reference humidity profile----------------
!
      else
        dqref=(qotsum-dstq-qsum)/den
      endif
!
!-------------- humidity does not increase with height------------------
!
      if(dqref<0.)then
        if (newswap) go to 810 ! new shallow cloud for both heat and moisture up   !jun04
        lbot=0
        ltop=lm
        ptop=pbot
        return
      endif
!
!----------------humidity at the cloud top------------------------------
!
      qrftp=qsum-dqref*psum
!
!----------------humidity profile---------------------------------------
!
      do l=ltop,lbot
        qrfkl=(fpk(l)-fptk)*dqref+qrftp
!
!***  too dry clouds not allowed
!
        tnew=(trefk(l)-tk(l))*tauksc+tk(l)
        qsatk(l)=pq0/pk(l)*exp(a2*(tnew-a3)/(tnew-a4))
        qnew=(qrfkl-qk(l))*tauksc+qk(l)
!
        if(qnew<qsatk(l)*rhlsc)then
          if (newswap) go to 810 ! new shallow cloud for both heat and moisture up   !jun04
          lbot=0
          ltop=lm
          ptop=pbot
          return
        endif
!
!-------------too moist clouds not allowed------------------------------
!
        if(qnew>qsatk(l)*rhhsc)then
          if (newswap) go to 810 ! new shallow cloud for both heat and moisture up   !jun04
          lbot=0
          ltop=lm
          ptop=pbot
          return
        endif

!
        thvref(l)=trefk(l)*rxnerk(l)*(qrfkl*p608+1.)
        qrefk(l)=qrfkl
      enddo
!
!------------------ eliminate clouds with bottoms too dry --------------
!!
!      qnew=(qrefk(lbot)-qk(lbot))*tauksc+qk(lbot)
!!
!      if(qnew<qk(lbot+1)*stresh)then  !!?? stresh too large!!
!        lbot=0
!        ltop=lm
!        ptop=pbot
!        return
!      endif
!!
!-------------- eliminate impossible slopes (betts,dtheta/dq)------------
!
      do l=ltop,lbot
        dtdp=(thvref(l-1)-thvref(l))/(prsmid(l)-prsmid(l-1))
!
        if(dtdp<epsdt)then
          if (newswap) go to 810 ! new shallow cloud for both heat and moisture up   !jun04
          lbot=0
          ltop=lm
          ptop=pbot
          return
        endif
!
      enddo
!-----------------------------------------------------------------------
!***  relaxation to reference profiles
!-----------------------------------------------------------------------
      if(mmntshal1) then
        facuv=tauksc*rdtcnvc*uvscals1
!
        ubar=0.
        vbar=0.
        sumdp=0.
!
        do l=ltop,lbot
          ubar=u(l)*dprs(l)+ubar
          vbar=v(l)*dprs(l)+vbar
          sumdp=dprs(l)+sumdp
        enddo
!
        rdpsum=1./sumdp
        ubar=ubar*rdpsum
        vbar=vbar*rdpsum
!
        do l=ltop,lbot
          dudt(l)=(ubar-u(l))*facuv
          dvdt(l)=(vbar-v(l))*facuv
        enddo
      endif
!-----------------------------------------------------------------------
              go to 820 ! relaxation
!-----------------------------------------------------------------------
!***  new cloud starts here
!-----------------------------------------------------------------------
 810  do l=1,lmh
        wcld(l)=0.
        rhk(l)=qk(l)/qsatk(l)
        thvmk(l)=tk(l)*rxnerk(l) !zj *(qk(l)*0p608+1.)
!----calculate updraft temperature along moist adiabat tref(l)----------
        if(prsmid(l).lt.plq) then
          call ttblex(itb,jtb,pl,pk(l),rdp,rdthe &
                     ,sthe,the0,thescnv(l),ttbl,tref(l))
        else
          call ttblex(itbq,jtbq,plq,pk(l),rdpq,rdtheq &
                     ,stheq,the0q,thescnv(l),ttblq,tref(l))
        endif
!-----------------------------------------------------------------------
        thmak(l)=tref(l)*rxnerk(l)
      enddo
!-------------mean rh and slopes within cloud---------------------------
      sumdp=0.
      sumrh=0.
      a11=0.
      a12=0.
      b1qsat=0.
      b2qsat=0.
      b1thvm=0.
      b1thma=0.
      b2thvm=0.
      b2thma=0.
      b1rh  =0.
      b2rh  =0.
!
      do l=ltop,lbot
        sumdp=dprs(l)+sumdp
        sumrh=rhk(l)*dprs(l)+sumrh
        a11=prsmid(l)**2*dprs(l)+a11
        a12=prsmid(l)*dprs(l)+a12
        b1qsat=qsatk(l)*prsmid(l)*dprs(l)+b1qsat
        b1thvm=thvmk(l)*prsmid(l)*dprs(l)+b1thvm
        b1thma=thmak(l)*prsmid(l)*dprs(l)+b1thma
        b1rh  =rhk  (l)*prsmid(l)*dprs(l)+b1rh
        b2qsat=qsatk(l)*dprs(l)+b2qsat
        b2thvm=thvmk(l)*dprs(l)+b2thvm
        b2thma=thmak(l)*dprs(l)+b2thma
        b2rh  =rhk  (l)*dprs(l)+b2rh
      enddo
!
      rhmean=sumrh/sumdp
!-------------no shallow convection if the cloud is saturated-----------
      if(rhmean.gt.0.95) then
        lbot=0
        ltop=lm
        ptop=pbot
        return
      endif
!-----------------------------------------------------------------------
      a21=a12
      a22=sumdp
!
      rden=1./(a11*a22-a12*a21)
!
      aqs=(b1qsat*a22-a12*b2qsat)*rden
      avm=(b1thvm*a22-a12*b2thvm)*rden
      ama=(b1thma*a22-a12*b2thma)*rden
      arh=(b1rh  *a22-a12*b2rh  )*rden
!
      bqs=(a11*b2qsat-b1qsat*a21)*rden
      bvm=(a11*b2thvm-b1thvm*a21)*rden
      bma=(a11*b2thma-b1thma*a21)*rden
      brh=(a11*b2rh  -b1rh  *a21)*rden
!-------------no shallow convection if the cloud moister on top---------
      if(arh.lt.0.) then !soft2
        lbot=0           !soft2
        ltop=lm          !soft2
        ptop=pbot        !soft2
        return        !soft2
      endif              !soft2
!-------------first guess t & q profiles--------------------------------
      adef=(1.-deftop)*2./(pk(lbot)-pk(ltop)) !soft2
!
      do l=ltop,lbot
        fk=(pk(l)-pk(ltop))*adef+deftop !soft2
        rhref=rhmean*fk                 !soft2
!
        wcld(l)=(1.-wdry)*rhref/(1.-wdry*rhref)
        trefk(l)=((1.-wcld(l))*(avm*pk(l)+bvm) &
                 +    wcld(l) *(ama*pk(l)+bma))/rxnerk(l)
        if (nodeep .and. abs(trefk(l)-tk(l))>5.) then
          lbot=0
          ltop=lm
          ptop=pbot
          return
        endif
        qrefk(l)=rhref*(aqs*prsmid(l)+bqs)
      enddo
!-------------enthalpy conservation-------------------------------------
      sumdp=0.
      sumdt=0.
      sumdq=0.
!
      do l=ltop,lbot
        sumdp=dprs(l)+sumdp
        sumdt=(tk(l)-trefk(l))*dprs(l)+sumdt
        sumdq=(qk(l)-qrefk(l))*dprs(l)+sumdq
      enddo
!
      rdpsum=1./sumdp
!
      tcorr=sumdt*rdpsum
      qcorr=sumdq*rdpsum
!
      do l=ltop,lbot
        trefk(l)=trefk(l)+tcorr
        qrefk(l)=qrefk(l)+qcorr
      enddo
!-----------------------------------------------------------------------
      dsq=0.
      dst=0.
!
      do l=ltop,lbot
        tkl=tk(l)
        diftl=(trefk(l)-tkl  )*tauksc
        difql=(qrefk(l)-qk(l))*tauksc
        dpot=dprs(l)/(tkl+tkl+diftl)
        dst=diftl      *dpot+dst
        dsq=difql*el(l)*dpot+dsq
      enddo
!
      dst=(dst+dst)*cp
      dsq=dsq+dsq
      dentpy=dst+dsq
!
      if(dentpy.lt.0.) then
        lbot=0
        ltop=lm
        ptop=pbot
        return
      endif
!-----------------------------------------------------------------------
      if(mmntshal2) then
!-------------mean momentum and profile slopes--------------------------
go to 8888
        b1u=0.
        b1v=0.
        b2u=0.
        b2v=0.
!
        do l=ltop,lbot
          b1u=u(l)*prsmid(l)*dprs(l)+b1u
          b1v=v(l)*prsmid(l)*dprs(l)+b1v
          b2u=u(l)*dprs(l)+b2u
          b2v=v(l)*dprs(l)+b2v
        enddo
!-----------------------------------------------------------------------
        au=(b1u*a22-a12*b2u)*rden
        av=(b1v*a22-a12*b2v)*rden
        bu=(a11*b2u-b1u*a21)*rden
        bv=(a11*b2v-b1v*a21)*rden
!-------------first guess u & v profiles--------------------------------
        do l=ltop,lbot
          urefk(l)=(1.-wcld(l))*u(l)+wcld(l)*(au*pk(l)+bu)
          vrefk(l)=(1.-wcld(l))*v(l)+wcld(l)*(av*pk(l)+bv)
        enddo
!-------------momentum conservation-------------------------------------
        sumdu=0.
        sumdv=0.
!
        do l=ltop,lbot
          sumdu=(u(l)-urefk(l))*dprs(l)+sumdu
          sumdv=(v(l)-vrefk(l))*dprs(l)+sumdv
        enddo
!
        ucorr=sumdu*rdpsum
        vcorr=sumdv*rdpsum
!
        do l=ltop,lbot
          urefk(l)=urefk(l)+ucorr
          vrefk(l)=vrefk(l)+vcorr
        enddo
!-----------------------------------------------------------------------
        facuv=tauksc*rdtcnvc*uvscals2

        do l=ltop,lbot
          dudt(l)=(urefk(l)-u(l))*facuv
          dvdt(l)=(vrefk(l)-v(l))*facuv
        enddo
8888 continue

!go to 7777
        ubar=0.
        vbar=0.
        sumdp=0.
!
        do l=ltop,lbot
          ubar=u(l)*dprs(l)+ubar
          vbar=v(l)*dprs(l)+vbar
          sumdp=dprs(l)+sumdp
        enddo
!
        rdpsum=1./sumdp
        ubar=ubar*rdpsum
        vbar=vbar*rdpsum
!
        facuv=tauksc*rdtcnvc*uvscals2
!
        do l=ltop,lbot
          dudt(l)=(ubar-u(l))*facuv
          dvdt(l)=(vbar-v(l))*facuv
        enddo
7777 continue
      endif
!--------------relaxation towards reference profiles--------------------
 820  do l=ltop,lbot
        dtdt(l)=(trefk(l)-tk(l))*tauksc*rdtcnvc
        dqdt(l)=(qrefk(l)-qk(l))*tauksc*rdtcnvc
      enddo
!-----------------------------------------------------------------------
!***  begin debugging convection
      if(print_diag)then
        do l=lbot,ltop,-1
          write(6,"(a,i3,4e12.4)") '{cu2 kflip,dt,dtdt,dq,dqdt = ' &
               ,lm+1-l,trefk(l)-tk(l),dtdt(l),qrefk(l)-qk(l),dqdt(l)
        enddo
      endif
!***  end debugging convection
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!scscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
!scscscscscscsc         end of shallow convection        scscscscscscscs
!scscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
!-----------------------------------------------------------------------
      end subroutine bmj
!-----------------------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!-----------------------------------------------------------------------
                           subroutine ttblex &
      (itbx,jtbx,plx,prsmid,rdpx,rdthex,sthe &
      ,the0,thesp,ttbl,tref)
!-----------------------------------------------------------------------
!     ******************************************************************
!     *                                                                *
!     *           extract temperature of the moist adiabat from        *
!     *                      the appropriate ttbl                      *
!     *                                                                *
!     ******************************************************************
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer(kind=kint),intent(in):: &
       itbx,jtbx
!
      real(kind=kfpt),intent(in):: &
       plx,prsmid,rdpx,rdthex,thesp
!
      real(kind=kfpt),dimension(itbx),intent(in):: &
       sthe,the0
!
      real(kind=kfpt),dimension(jtbx,itbx),intent(in):: &
       ttbl
!
      real(kind=kfpt),intent(out):: &
       tref
!-----------------------------------------------------------------------
      integer(kind=kint):: &
       iptb,ithtb
!
      real(kind=kfpt):: &
       bthe00k,bthe10k,bthk,pk,pp,qq,sthe00k,sthe10k,sthk &
      ,t00k,t01k,t10k,t11k,tpk,tthk
!-----------------------------------------------------------------------
!----------------scaling pressure & tt table index----------------------
!-----------------------------------------------------------------------
      pk=prsmid
      tpk=(pk-plx)*rdpx
      qq=tpk-aint(tpk)
      iptb=int(tpk)+1
!----------------keeping indices within the table-----------------------
      if(iptb<1)then
        iptb=1
        qq=0.
      endif
!
      if(iptb>=itbx)then
        iptb=itbx-1
        qq=0.
      endif
!----------------base and scaling factor for thetae---------------------
      bthe00k=the0(iptb)
      sthe00k=sthe(iptb)
      bthe10k=the0(iptb+1)
      sthe10k=sthe(iptb+1)
!----------------scaling the & tt table index---------------------------
      bthk=(bthe10k-bthe00k)*qq+bthe00k
      sthk=(sthe10k-sthe00k)*qq+sthe00k
      tthk=(thesp-bthk)/sthk*rdthex
      pp=tthk-aint(tthk)
      ithtb=int(tthk)+1
!----------------keeping indices within the table-----------------------
      if(ithtb<1)then
        ithtb=1
        pp=0.
      endif
!
      if(ithtb>=jtbx)then
        ithtb=jtbx-1
        pp=0.
      endif
!----------------temperature at four surrounding tt table pts.----------
      t00k=ttbl(ithtb,iptb)
      t10k=ttbl(ithtb+1,iptb)
      t01k=ttbl(ithtb,iptb+1)
      t11k=ttbl(ithtb+1,iptb+1)
!-----------------------------------------------------------------------
!----------------parcel temperature-------------------------------------
!-----------------------------------------------------------------------
      tref=(t00k+(t10k-t00k)*pp+(t01k-t00k)*qq                           &
      &    +(t00k-t10k-t01k+t11k)*pp*qq)
!-----------------------------------------------------------------------
      end subroutine ttblex
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      subroutine bmj_init
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      real(kind=kfpt),parameter:: &
       eliwv=2.683e6,eps=1.e-9
!
      integer(kind=kint):: &
       kth,kthm,kthm1,kp,kpm,kpm1
!
      real(kind=kfpt):: &
       rxner,dp,dqs,dth,dthe,p,qs,qs0k,sqsk,sthek &
      ,th,the0k,denom
!
      real(kind=kfpt), dimension(jtb):: &
       app,apt,aqp,aqt,pnew,pold,qsnew,qsold &
      ,thenew,theold,tnew,told,y2p,y2t
!
      real(kind=kfpt),dimension(jtbq):: &
       aptq,aqtq,thenewq,theoldq &
      ,tnewq,toldq,y2tq
!
!-----------------------------------------------------------------------
!----------------coarse look-up table for saturation point--------------
!-----------------------------------------------------------------------
!
      kthm=jtb
      kpm=itb
      kthm1=kthm-1
      kpm1=kpm-1
!
      dth=(thh-thl)/float(kthm-1)
      dp =(ph -pl )/float(kpm -1)
!
      th=thl-dth
!-----------------------------------------------------------------------
!
      do 100 kth=1,kthm
!
      th=th+dth
      p=pl-dp
!
      do kp=1,kpm
        p=p+dp
        rxner=(100000./p)**cappa
        denom=th-a4*rxner
        if (denom>eps) then
           qsold(kp)=pq0/p*exp(a2*(th-a3*rxner)/denom)
        else
           qsold(kp)=0.
        endif
        pold(kp)=p
      enddo
!
      qs0k=qsold(1)
      sqsk=qsold(kpm)-qsold(1)
      qsold(1  )=0.
      qsold(kpm)=1.
!
      do kp=2,kpm1
        qsold(kp)=(qsold(kp)-qs0k)/sqsk
        if((qsold(kp)-qsold(kp-1))<eps)qsold(kp)=qsold(kp-1)+eps
      enddo
!
      qs0(kth)=qs0k
      qs0_exp(kth)=qs0k
      sqs(kth)=sqsk
      sqs_exp(kth)=sqsk
!-----------------------------------------------------------------------
      qsnew(1  )=0.
      qsnew(kpm)=1.
      dqs=1./float(kpm-1)
!
      do kp=2,kpm1
        qsnew(kp)=qsnew(kp-1)+dqs
      enddo
!
      y2p(1   )=0.
      y2p(kpm )=0.
!
      call spline(jtb,kpm,qsold,pold,y2p,kpm,qsnew,pnew,app,aqp)
!
      do kp=1,kpm
        ptbl(kp,kth)=pnew(kp)
        ptbl_exp(kp,kth)=pnew(kp)
      enddo
!-----------------------------------------------------------------------
  100 continue
!-----------------------------------------------------------------------
!------------coarse look-up table for t(p) from constant the------------
!-----------------------------------------------------------------------
      p=pl-dp
!
      do 200 kp=1,kpm
!
      p=p+dp
      th=thl-dth
!
      do kth=1,kthm
        th=th+dth
        rxner=(1.e5/p)**cappa
        denom=th-a4*rxner
        if (denom>eps) then
           qs=pq0/p*exp(a2*(th-a3*rxner)/denom)
        else
           qs=0.
        endif
!        qs=pq0/p*exp(a2*(th-a3*rxner)/(th-a4*rxner))
        told(kth)=th/rxner
        theold(kth)=th*exp(eliwv*qs/(cp*told(kth)))
      enddo
!
      the0k=theold(1)
      sthek=theold(kthm)-theold(1)
      theold(1   )=0.
      theold(kthm)=1.
!
      do kth=2,kthm1
        theold(kth)=(theold(kth)-the0k)/sthek
        if((theold(kth)-theold(kth-1)).lt.eps) &
      &      theold(kth)=theold(kth-1)  +  eps
      enddo
!
      the0(kp)=the0k
      sthe(kp)=sthek
!-----------------------------------------------------------------------
      thenew(1  )=0.
      thenew(kthm)=1.
      dthe=1./float(kthm-1)
!
      do kth=2,kthm1
        thenew(kth)=thenew(kth-1)+dthe
      enddo
!
      y2t(1   )=0.
      y2t(kthm)=0.
!
      call spline(jtb,kthm,theold,told,y2t,kthm,thenew,tnew,apt,aqt)
!
      do kth=1,kthm
        ttbl(kth,kp)=tnew(kth)
      enddo
!-----------------------------------------------------------------------
  200 continue
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!---------------fine look-up table for saturation point-----------------
!-----------------------------------------------------------------------
      kthm=jtbq
      kpm=itbq
      kthm1=kthm-1
      kpm1=kpm-1
!
      dth=(thhq-thl)/float(kthm-1)
      dp=(ph-plq)/float(kpm-1)
!
      th=thl-dth
      p=plq-dp
!-----------------------------------------------------------------------
!---------------fine look-up table for t(p) from constant the-----------
!-----------------------------------------------------------------------
      do 300 kp=1,kpm
!
      p=p+dp
      th=thl-dth
!
      do kth=1,kthm
        th=th+dth
        rxner=(1.e5/p)**cappa
        denom=th-a4*rxner
        if (denom>eps) then
           qs=pq0/p*exp(a2*(th-a3*rxner)/denom)
        else
           qs=0.
        endif
!        qs=pq0/p*exp(a2*(th-a3*rxner)/(th-a4*rxner))
        toldq(kth)=th/rxner
        theoldq(kth)=th*exp(eliwv*qs/(cp*toldq(kth)))
      enddo
!
      the0k=theoldq(1)
      sthek=theoldq(kthm)-theoldq(1)
      theoldq(1   )=0.
      theoldq(kthm)=1.
!
      do kth=2,kthm1
        theoldq(kth)=(theoldq(kth)-the0k)/sthek
        if((theoldq(kth)-theoldq(kth-1))<eps) &
             theoldq(kth)=theoldq(kth-1)+eps
      enddo
!
      the0q(kp)=the0k
      stheq(kp)=sthek
!-----------------------------------------------------------------------
      thenewq(1  )=0.
      thenewq(kthm)=1.
      dthe=1./float(kthm-1)
!
      do kth=2,kthm1
        thenewq(kth)=thenewq(kth-1)+dthe
      enddo
!
      y2tq(1   )=0.
      y2tq(kthm)=0.
!
      call spline(jtbq,kthm,theoldq,toldq,y2tq,kthm &
                  ,thenewq,tnewq,aptq,aqtq)
!
      do kth=1,kthm
        ttblq(kth,kp)=tnewq(kth)
      enddo
!-----------------------------------------------------------------------
  300 continue
!-----------------------------------------------------------------------
      end subroutine bmj_init
!-----------------------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
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
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      end module module_cu_bmj
!
!-----------------------------------------------------------------------
