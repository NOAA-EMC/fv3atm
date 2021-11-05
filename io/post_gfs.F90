!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
module post_gfs

  use module_fv3_io_def,    only : wrttasks_per_group,filename_base
  use write_internal_state, only : wrt_internal_state

  implicit none

  include 'mpif.h'

  integer mype, nbdl
  logical setvar_atmfile, setvar_sfcfile, read_postcntrl
  public  post_run_gfs, post_getattr_gfs

  contains

  subroutine post_run_gfs(wrt_int_state,mypei,mpicomp,lead_write,      &
             mynfhr,mynfmin,mynfsec)
!
!  revision history:
!     Jul 2019    J. Wang      create interface to run inline post for FV3
!     Apr 2021    R. Sun       Added variables for Thomspon MP
!
!-----------------------------------------------------------------------
!*** run post on write grid comp
!-----------------------------------------------------------------------
!
      use ctlblk_mod, only : komax,ifhr,ifmin,modelname,datapd,fld_info, &
                             npset,grib,gocart_on,icount_calmict, jsta,  &
                             jend,im, nsoil, filenameflat
      use gridspec_mod, only : maptype, gridtype
      use grib2_module, only : gribit2,num_pset,nrecout,first_grbtbl
      use xml_perl_data,only : paramset
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(wrt_internal_state),intent(in)       :: wrt_int_state
      integer,intent(in)                        :: mypei
      integer,intent(in)                        :: mpicomp
      integer,intent(in)                        :: lead_write
      integer,intent(in)                        :: mynfhr
      integer,intent(in)                        :: mynfmin
      integer,intent(in)                        :: mynfsec
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      integer n,nwtpg,ieof,lcntrl,ierr,i,j,k,jts,jte,mynsoil
      integer,allocatable  :: jstagrp(:),jendgrp(:)
      integer,save         :: kpo,kth,kpv
      logical,save         :: log_postalct=.false.
      real,dimension(komax),save :: po, th, pv
      logical        :: Log_runpost
      character(255) :: post_fname*255

      integer,save :: iostatusD3D=-1
!
      real(kind=8)   :: btim0, btim1, btim2, btim3,btim4,btim5,btim6,btim7
!
!      print *,'in post_run start'
!-----------------------------------------------------------------------
!*** set up dimensions
!-----------------------------------------------------------------------
!
      btim0 = MPI_Wtime()

      modelname = "GFS"
      grib      = "grib2"
      gridtype  = "A"
      nsoil     = 4
      mype      = mypei
      nwtpg     = wrt_int_state%petcount
      jts       = wrt_int_state%lat_start              !<-- Starting J of this write task's subsection
      jte       = wrt_int_state%lat_end                !<-- Ending J of this write task's subsection
      maptype   = wrt_int_state%post_maptype
      nbdl      = wrt_int_state%FBCount

      if(mype==0) print *,'in post_run,jts=',jts,'jte=',jte,'nwtpg=',nwtpg,'nwtpg=',nwtpg, &
        'jts=',jts,'jte=',jte,'maptype=',maptype,'nbdl=',nbdl,'log_postalct=',log_postalct

!
!-----------------------------------------------------------------------
!*** set up fields to run post
!-----------------------------------------------------------------------
!
      if (.not.log_postalct) then
!
        allocate(jstagrp(nwtpg),jendgrp(nwtpg))
!
        do n=0,nwtpg-1
          jstagrp(n+1) = wrt_int_state%lat_start_wrtgrp(n+1)
          jendgrp(n+1) = wrt_int_state%lat_end_wrtgrp  (n+1)
        enddo
        if(mype==0) print *,'in post_run,jstagrp=',jstagrp,'jendgrp=',jendgrp

!-----------------------------------------------------------------------
!*** read namelist for pv,th,po
!-----------------------------------------------------------------------
!
        call read_postnmlt(kpo,kth,kpv,po,th,pv,wrt_int_state%post_nlunit, &
                           wrt_int_state%post_namelist)
!
!-----------------------------------------------------------------------
!*** allocate post variables
!-----------------------------------------------------------------------
!
!     if(mype==0) print *,'in post_run,be post_alctvars, dim=',wrt_int_state%im, &
!       wrt_int_state%jm, wrt_int_state%lm,'mype=',mype,'wrttasks_per_group=', &
!       wrttasks_per_group,'lead_write=',lead_write,'jts=',jts,'jte=',jte,   &
!       'jstagrp=',jstagrp,'jendgrp=',jendgrp
        call post_alctvars(wrt_int_state%im,wrt_int_state%jm,        &
          wrt_int_state%lm,mype,wrttasks_per_group,lead_write,    &
          mpicomp,jts,jte,jstagrp,jendgrp)
!
!-----------------------------------------------------------------------
!*** read namelist for pv,th,po
!-----------------------------------------------------------------------
!
        log_postalct = .true.
        first_grbtbl = .true.
        read_postcntrl = .true.
!
      ENDIF
!
!-----------------------------------------------------------------------
!*** fill post variables with values from forecast results
!-----------------------------------------------------------------------
!
      ifhr  = mynfhr
      ifmin = mynfmin
      if (ifhr == 0 ) ifmin = 0
      if(mype==0) print *,'bf set_postvars,ifmin=',ifmin,'ifhr=',ifhr
      setvar_atmfile=.false.
      setvar_sfcfile=.false.
      call set_postvars_gfs(wrt_int_state,mpicomp,setvar_atmfile,   &
           setvar_sfcfile)

!       print *,'af set_postvars,setvar_atmfile=',setvar_atmfile,  &
!        'setvar_sfcfile=',setvar_sfcfile
!
      if (setvar_atmfile.and.setvar_sfcfile) then
! 20190807 no need to call microinit for GFDLMP
!        call MICROINIT
!
        if(grib=="grib2" .and. read_postcntrl) then
          if (ifhr == 0) then
            filenameflat = 'postxconfig-NT_FH00.txt'
            call read_xml()
            if(mype==0) print *,'af read_xml at fh00,name=',trim(filenameflat)
          else if(ifhr > 0) then
            filenameflat = 'postxconfig-NT.txt'
            if(associated(paramset)) then
              if( size(paramset)>0) then
                do i=1,size(paramset)
                  if (associated(paramset(i)%param)) then
                    if (size(paramset(i)%param)>0) then
                      deallocate(paramset(i)%param)
                      nullify(paramset(i)%param)
                    endif
                  endif
                enddo
                deallocate(paramset)
                nullify(paramset)
              endif
            endif
            num_pset = 0
            call read_xml()
            if(mype==0) print *,'af read_xml,name=',trim(filenameflat),'ifhr=',ifhr
            read_postcntrl = .false.
          endif
        endif
!
        IEOF  = 0
        npset = 0
        icount_calmict = 0
        do while( IEOF == 0)
!
          if(grib == "grib2") then
            npset = npset + 1
            call set_outflds(kth,th,kpv,pv)
            if(allocated(datapd))deallocate(datapd)
            allocate(datapd(wrt_int_state%im,jte-jts+1,nrecout+100))
!$omp parallel do default(none),private(i,j,k),shared(nrecout,jend,jsta,im,datapd)
            do k=1,nrecout+100
              do j=1,jend+1-jsta
                do i=1,im
                  datapd(i,j,k) = 0.
                enddo
              enddo
            enddo
            call get_postfilename(post_fname)
            if (mype==0) write(0,*)'post_fname=',trim(post_fname)
!
            if ( ieof == 0) call process(kth,kpv,th(1:kth),pv(1:kpv),iostatusD3D)
!
            call mpi_barrier(mpicomp,ierr)
            call gribit2(post_fname)
            if(allocated(datapd))deallocate(datapd)
            if(allocated(fld_info))deallocate(fld_info)
            if(npset >= num_pset) exit

          endif
!
        enddo
!
      endif

    end subroutine post_run_gfs
!
!-----------------------------------------------------------------------
!
    subroutine post_getattr_gfs(wrt_int_state)
!
      use esmf
      use ctlblk_mod,           only: im, jm, mpi_comm_comp
      use masks,                only: gdlat, gdlon, dx, dy
      use gridspec_mod,         only: latstart, latlast, lonstart,    &
                                      lonlast, cenlon, cenlat
!
      implicit none
!
      type(wrt_internal_state),intent(inout)    :: wrt_int_state
!
! local variable
      integer i,j,k,n,kz, attcount, nfb
      integer ni,naryi,nr4,nr8,rc
      integer aklen,varival
      real(4) varr4val
      real(8) varr8val
      character(80) attName, hydrostatics, fldname
      type(ESMF_TypeKind_Flag)           :: typekind
      real(4), dimension(:), allocatable :: ak4,bk4
      real(8), dimension(:), allocatable :: ak8,bk8
      type(ESMF_FieldBundle)             :: fldbundle
!
! field bundle
     do nfb=1, wrt_int_state%FBcount
       fldbundle = wrt_int_state%wrtFB(nfb) 

! look at the field bundle attributes
      call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
        attnestflag=ESMF_ATTNEST_OFF, Count=attcount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,file=__FILE__))return  ! bail out
!
      aklen=0.
      do i=1, attCount

        call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
          attnestflag=ESMF_ATTNEST_OFF, attributeIndex=i, name=attName, &
          typekind=typekind, itemCount=n,  rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,file=__FILE__))return  ! bail out

        if (typekind==ESMF_TYPEKIND_I4 ) then
          if(n==1) then
            call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
              name=trim(attName), value=varival, rc=rc)
            if (trim(attName) == 'ncnsto') wrt_int_state%ntrac=varival
            if (trim(attName) == 'ncld')   wrt_int_state%ncld=varival
            if (trim(attName) == 'nsoil')  wrt_int_state%nsoil=varival
            if (trim(attName) == 'fhzero') wrt_int_state%fhzero=varival
            if (trim(attName) == 'imp_physics') wrt_int_state%imp_physics=varival
          endif
        else if (typekind==ESMF_TYPEKIND_R4) then
          if(n==1) then
            call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
              name=trim(attName), value=varr4val, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=__FILE__)) return  ! bail out
            if (trim(attName) == 'dtp')   then
               wrt_int_state%dtp=varr4val
            endif
          else if(n>1) then
            if(trim(attName) =="ak") then
              if(allocated(wrt_int_state%ak)) deallocate(wrt_int_state%ak)
              allocate(wrt_int_state%ak(n))
              call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                name=trim(attName), valueList=wrt_int_state%ak, rc=rc)
              wrt_int_state%lm = n-1
            else if(trim(attName) =="bk") then
              if(allocated(wrt_int_state%bk)) deallocate(wrt_int_state%bk)
              allocate(wrt_int_state%bk(n))
              call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                name=trim(attName), valueList=wrt_int_state%bk, rc=rc)
            endif
          endif
        else if (typekind==ESMF_TYPEKIND_R8) then
          if(n==1) then
            call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
              name=trim(attName), value=varr8val, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=__FILE__)) return  ! bail out
            if (trim(attName) == 'dtp')   then
               wrt_int_state%dtp=varr8val
            endif
          else if(n>1) then
            if(trim(attName) =="ak") then
              if(allocated(wrt_int_state%ak)) deallocate(wrt_int_state%ak)
              allocate(wrt_int_state%ak(n))
              call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
              name=trim(attName), valueList=wrt_int_state%ak, rc=rc)
              wrt_int_state%lm = n-1
            else if(trim(attName) =="bk") then
              if(allocated(wrt_int_state%bk)) deallocate(wrt_int_state%bk)
              allocate(wrt_int_state%bk(n))
              call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
              name=trim(attName), valueList=wrt_int_state%bk, rc=rc)
            endif
            wrt_int_state%lm = size(wrt_int_state%ak) - 1
          endif
        endif
!
      enddo
!
      enddo !end nfb
!      print *,'in post_getattr, dtp=',wrt_int_state%dtp
!
    end subroutine post_getattr_gfs
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
    subroutine set_postvars_gfs(wrt_int_state,mpicomp,setvar_atmfile,   &
                                setvar_sfcfile)
!
!  revision history:
!     Jul 2019    J. Wang      Initial code
!
!-----------------------------------------------------------------------
!*** set up post fields from nmint_state
!-----------------------------------------------------------------------
!
      use esmf
      use vrbls3d,     only: t, q, uh, vh, wh, alpint, dpres, zint, zmid, o3,  &
                             qqr, qqs, cwm, qqi, qqw, qqg, omga, cfr, pmid,    &
                             q2, rlwtt, rswtt, tcucn, tcucns, train, el_pbl,   &
                             pint, exch_h, ref_10cm, qqni,qqnr,qqnwfa,qqnifa
      use vrbls2d,     only: f, pd, sigt4, fis, pblh, ustar, z0, ths, qs, twbs,&
                             qwbs, avgcprate, cprate, avgprec, prec, lspa, sno,&
                             cldefi, th10, q10, tshltr, pshltr, tshltr, albase,&
                             avgalbedo, avgtcdc, czen, czmean, mxsnal, radot,  &
                             cfrach, cfracl, cfracm, avgcfrach, qshltr,        &
                             avgcfracl, avgcfracm, cnvcfr, islope, cmc, grnflx,&
                             vegfrc, acfrcv, ncfrcv, acfrst, ncfrst, ssroff,   &
                             bgroff, rlwin,      &
                             rlwtoa, cldwork, alwin, alwout, alwtoa, rswin,    &
                             rswinc, rswout, aswin, auvbin, auvbinc, aswout,   &
                             aswtoa, sfcshx, sfclhx, subshx, snopcx, sfcux,    &
                             sfcvx, sfcuvx, gtaux, gtauy, potevp, u10, v10,    &
                             smstav, smstot, ivgtyp, isltyp, sfcevp, sfcexc,   &
                             acsnow, acsnom, sst, thz0, qz0, uz0, vz0, ptop,   &
                             htop, pbot, hbot, ptopl, pbotl, ttopl, ptopm,     &
                             pbotm, ttopm, ptoph, pboth, pblcfr, ttoph, runoff,&
                             maxtshltr, mintshltr, maxrhshltr, minrhshltr,     &
                             dzice, smcwlt, suntime, fieldcapa, htopd, hbotd,  &
                             htops, hbots, aswintoa, maxqshltr, minqshltr,     &
                             acond, sr, u10h, v10h, avgedir, avgecan,          &
                             avgetrans, avgesnow, avgprec_cont, avgcprate_cont,&
                             avisbeamswin, avisdiffswin, airbeamswin, airdiffswin, &
                             alwoutc, alwtoac, aswoutc, aswtoac, alwinc, aswinc,& 
                             avgpotevp, snoavg, ti, si, cuppt
      use soil,        only: sldpth, sh2o, smc, stc
      use masks,       only: lmv, lmh, htm, vtm, gdlat, gdlon, dx, dy, hbm2, sm, sice
      use ctlblk_mod,  only: im, jm, lm, lp1, jsta, jend, jsta_2l, jend_2u, jsta_m,jend_m, &
                             lsm, pt, imp_physics, spval, mpi_comm_comp, gdsdegr,  &
                             tprec, tclod, trdlw, trdsw, tsrfc, tmaxmin, theat, &
                             ardlw, ardsw, asrfc, avrain, avcnvc, iSF_SURFACE_PHYSICS,&
                             td3d, idat, sdat, ifhr, ifmin, dt, nphs, dtq2, pt_tbl, &
                             alsl, spl, ihrst 
      use params_mod,  only: erad, dtr, capa, p1000
      use gridspec_mod,only: latstart, latlast, lonstart, lonlast, cenlon, cenlat
      use lookup_mod,  only: thl, plq, ptbl, ttbl, rdq, rdth, rdp, rdthe, pl,   &
                             qs0, sqs, sthe, ttblq, rdpq, rdtheq, stheq, the0q, the0
      use physcons,    only: grav => con_g, fv => con_fvirt, rgas => con_rd,    &
                             eps => con_eps, epsm1 => con_epsm1
      use rqstfld_mod
!
!      use write_internal_state, only: wrt_internal_state
!
!-----------------------------------------------------------------------
!
      implicit none
!
      include 'mpif.h'
!
!-----------------------------------------------------------------------
!
      type(wrt_internal_state),intent(in) :: wrt_int_state
      integer,intent(in)                  :: mpicomp
      logical,intent(inout)               :: setvar_atmfile,setvar_sfcfile
!
!-----------------------------------------------------------------------
!
      integer i, ip1, j, l, k, n, iret, ibdl, rc, kstart, kend
      integer ista,iend,fieldDimCount,gridDimCount,ncount_field
      integer jdate(8)
      logical foundland, foundice, found
      real(4) rinc(5)
      real    tlmh,RADI,TMP,ES,TV,RHOAIR,tem,tstart,dtp
      real, dimension(:),allocatable    :: ak5, bk5
      real(4),dimension(:,:),pointer    :: arrayr42d
      real(8),dimension(:,:),pointer    :: arrayr82d
      real(4),dimension(:,:,:),pointer  :: arrayr43d
      real(8),dimension(:,:,:),pointer  :: arrayr83d
      real,dimension(:),    allocatable :: slat,qstl
      real,external::FPVSNEW
      real,dimension(:,:),allocatable :: dummy, p2d, t2d, q2d,  qs2d,  &
                             cw2d, cfr2d
      character(len=80)              :: fieldname, wrtFBName
      type(ESMF_Grid)                :: wrtGrid
      type(ESMF_Field)               :: theField
      type(ESMF_Field), allocatable  :: fcstField(:)
      type(ESMF_TypeKind_Flag)       :: typekind
!
!-----------------------------------------------------------------------
!***  INTEGER SCALAR/1D HISTORY VARIABLES
!-----------------------------------------------------------------------
!
      imp_physics = wrt_int_state%imp_physics       !set GFS mp physics to 99 for Zhao scheme
      dtp         = wrt_int_state%dtp
      iSF_SURFACE_PHYSICS = 2
      spval = 9.99e20

!
! nems gfs has zhour defined
      tprec   = float(wrt_int_state%fhzero)
      tclod   = tprec
      trdlw   = tprec
      trdsw   = tprec
      tsrfc   = tprec
      tmaxmin = tprec
      td3d    = tprec
      if(mype==0)print*,'MP_PHYSICS= ',imp_physics,'nbdl=',nbdl, 'tprec=',tprec,'tclod=',tclod, &
       'dtp=',dtp,'tmaxmin=',tmaxmin

!      write(6,*) 'maptype and gridtype is ', maptype,gridtype
!
!$omp parallel do default(shared),private(i,j)
      do j=jsta,jend
        do  i=1,im
          gdlat(i,j) = wrt_int_state%latPtr(i,j)
          gdlon(i,j) = wrt_int_state%lonPtr(i,j)
        enddo
      enddo
!
      lonstart = nint(wrt_int_state%lonstart*gdsdegr)
      lonlast  = nint(wrt_int_state%lonlast*gdsdegr)
      latstart = nint(wrt_int_state%latstart*gdsdegr)
      latlast  = nint(wrt_int_state%latlast*gdsdegr)
!      print*,'latstart,latlast B bcast= ',latstart,latlast
!      print*,'lonstart,lonlast B bcast= ',lonstart,lonlast

!$omp parallel do default(none),private(i,j,ip1), &
!$omp&  shared(jsta,jend_m,im,dx,gdlat,gdlon,dy)
      do j = jsta, jend_m
        do i = 1, im
          ip1 = i + 1
          if (ip1 > im) ip1 = ip1 - im
          dx(i,j) = erad*cos(gdlat(i,j)*dtr)*(gdlon(ip1,j)-gdlon(i,j))*dtr
          dy(i,j)  = erad*(gdlat(i,j)-gdlat(i,j+1))*dtr  ! like A*DPH
        end do
      end do
!
      if(.not. allocated(ak5)) allocate(ak5(lm+1),bk5(lm+1))
      do i=1,lm+1
        ak5(i) = wrt_int_state%ak(i)
        bk5(i) = wrt_int_state%bk(i)
      enddo

!$omp parallel do default(none) private(i,j) shared(jsta,jend,im,f,gdlat)
      do j=jsta,jend
        do i=1,im
          f(I,J) = 1.454441e-4*sin(gdlat(i,j)*dtr)   ! 2*omeg*sin(phi)
        end do
      end do
!
! GFS does not output PD
      pt    = ak5(1)

! GFS may not have model derived radar ref.
!                        TKE
!                        cloud amount
!$omp parallel do default(none),private(i,j,l), &
!$omp& shared(lm,jsta,jend,im,spval,ref_10cm,q2,cfr)
      do l=1,lm
        do j=jsta,jend
          do i=1,im
            ref_10cm(i,j,l) = SPVAL
            q2(i,j,l) = SPVAL
            cfr(i,j,l) = SPVAL
          enddo
        enddo
      enddo

! GFS does not have surface specific humidity
!                   inst sensible heat flux
!                   inst latent heat flux
!$omp parallel do default(none),private(i,j),shared(jsta,jend,im,spval,qs,twbs,qwbs,ths)
      do j=jsta,jend
        do i=1,im
          qs(i,j) = SPVAL
          twbs(i,j) = SPVAL
          qwbs(i,j) = SPVAL
          ths(i,j) = SPVAL
        enddo
      enddo

! GFS set up DT to compute accumulated fields, set it to one
      dtq2 = wrt_int_state%dtp
      nphs = 2.
      dt   = dtq2/nphs
!
! GFS does not have convective cloud efficiency
!                   similated precip
!                   10 m theta
!                   10 m humidity
!                   snow free albedo
!$omp parallel do default(none), private(i,j), shared(jsta,jend,im,spval), &
!$omp& shared(cldefi,lspa,th10,q10,albase)
      do j=jsta,jend
        do i=1,im
          cldefi(i,j) = SPVAL
          lspa(i,j) = SPVAL
          th10(i,j) = SPVAL
          q10(i,j) = SPVAL
          albase(i,j) = SPVAL
        enddo
      enddo

! GFS does not have convective precip
!$omp parallel do default(none) private(i,j) shared(jsta,jend,im,cprate)
      do j=jsta,jend
        do i=1,im
          cprate(i,j) = 0.
        enddo
      enddo

! GFS probably does not use zenith angle, czen, czmean
!                       inst surface outgoing longwave, radot
!                       inst cloud fraction for high, middle, and low cloud,
!                            cfrach
!                       inst ground heat flux, grnflx
!$omp parallel do default(none) private(i,j) shared(jsta,jend,im,spval), &
!$omp& shared(czen,czmean,radot,cfrach,cfracl,cfracm,grnflx)
      do j=jsta,jend
        do i=1,im
          czen(i,j)   = SPVAL
          czmean(i,j) = SPVAL
          radot(i,j)  = SPVAL
          cfrach(i,j) = SPVAL
          cfracl(i,j) = SPVAL
          cfracm(i,j) = SPVAL
          grnflx(i,j) = SPVAL
        enddo
      enddo
!
! GFS doesn not yet output soil layer thickness, assign SLDPTH to be the same as nam
      sldpth(1) = 0.10
      sldpth(2) = 0.3
      sldpth(3) = 0.6
      sldpth(4) = 1.0

! GFS does not output time averaged convective and strat cloud fraction, set acfrcv to spval, n
! cfrcv to 1
!                     time averaged cloud fraction, set acfrst to spval, ncfrst to 1
!                     UNDERGROUND RUNOFF, bgroff
!                     inst incoming sfc longwave, rlwin
!                     inst model top outgoing longwave,rlwtoa
!                     inst incoming sfc shortwave, rswin
!                     inst incoming clear sky sfc shortwave, rswinc
!                     inst outgoing sfc shortwave, rswout
!                     snow phase change heat flux, snopcx
! GFS does not use total momentum flux,sfcuvx
!$omp parallel do default(none),private(i,j),shared(jsta,jend,im,spval), &
!$omp& shared(acfrcv,ncfrcv,acfrst,ncfrst,bgroff,rlwin,rlwtoa,rswin,rswinc,rswout,snopcx,sfcuvx)
      do j=jsta,jend
        do i=1,im
          acfrcv(i,j) = spval
          ncfrcv(i,j) = 1.0
          acfrst(i,j) = spval
          ncfrst(i,j) = 1.0
          bgroff(i,j) = spval
          rlwin(i,j)  = spval
          rlwtoa(i,j) = spval
          rswin(i,j)  = spval
          rswinc(i,j) = spval
          rswout(i,j) = spval
          snopcx(i,j) = spval
          sfcuvx(i,j) = spval
        enddo
      enddo

! GFS incoming sfc longwave has been averaged over 6 hr bucket, set ARDLW to 1
      ardlw  = 1.0
! GFS incoming sfc longwave has been averaged, set ARDLW to 1
      ardsw = 1.0
! GFS surface flux has been averaged, set  ASRFC to 1
      asrfc = 1.0

! GFS does not have temperature tendency due to long wave radiation
!                   temperature tendency due to short wave radiation
!                   temperature tendency due to latent heating from convection
!                   temperature tendency due to latent heating from grid scale
      do l=1,lm
!$omp parallel do default(none),private(i,j),shared(jsta_2l,jend_2u,im,spval,l), &
!$omp& shared(rlwtt,rswtt,tcucn,tcucns,train)
        do j=jsta_2l,jend_2u
          do i=1,im
            rlwtt(i,j,l) = spval
            rswtt(i,j,l)  = spval
            tcucn(i,j,l)  = spval
            tcucns(i,j,l) = spval
            train(i,j,l)  = spval
          enddo
        enddo
      enddo

! set avrain to 1
      avrain = 1.0
      avcnvc = 1.0
      theat  = 6.0 ! just in case GFS decides to output T tendency

! GFS does not have temperature tendency due to latent heating from grid scale
      train  = spval

! GFS does not have soil moisture availability, smstav
!                   accumulated surface evaporatio, sfcevp
!                   averaged accumulated snow, acsnow
!                   snow melt,acsnom
!                   humidity at roughness length, qz0
!                   u at roughness length, uz0
!                   v at roughness length, vz0
!                   shelter rh max, maxrhshltr
!                   shelter rh min, minrhshltr
!$omp parallel do default(none),private(i,j),shared(jsta_2l,jend_2u,im,spval), &
!$omp& shared(smstav,sfcevp,acsnow,acsnom,qz0,uz0,vz0,maxrhshltr,minrhshltr)
      do j=jsta_2l,jend_2u
        do i=1,im
          smstav(i,j) = spval
          sfcevp(i,j) = spval
          acsnow(i,j) = spval
          acsnom(i,j) = spval
          qz0(i,j)    = spval
          uz0(i,j)    = spval
          vz0(i,j)    = spval
          maxrhshltr(i,j) = SPVAL
          minrhshltr(i,j) = SPVAL
        enddo
      enddo

! GFS does not have mixing length,el_pbl
!                   exchange coefficient, exch_h
      do l=1,lm
!$omp parallel do default(none),private(i,j),shared(jsta_2l,jend_2u,im,l,spval,el_pbl,exch_h)
        do j=jsta_2l,jend_2u
          do i=1,im
            el_pbl(i,j,l) = spval
            exch_h(i,j,l) = spval
          enddo
        enddo
      enddo

! GFS does not have deep convective cloud top and bottom fields
!$omp parallel do default(none),private(i,j),shared(jsta_2l,jend_2u,im,spval), &
!$omp& shared(htopd,hbotd,htops,hbots,cuppt)
      do j=jsta_2l,jend_2u
        do i=1,im
          htopd(i,j) = SPVAL
          hbotd(i,j) = SPVAL
          htops(i,j) = SPVAL
          hbots(i,j) = SPVAL
          cuppt(i,j) = SPVAL
        enddo
      enddo
!
! get inital date
      sdat(1)  = wrt_int_state%idate(2)   !month
      sdat(2)  = wrt_int_state%idate(3)   !day
      sdat(3)  = wrt_int_state%idate(1)   !year
      ihrst    = wrt_int_state%idate(4)   !hour

      idat(1)  = wrt_int_state%fdate(2)
      idat(2)  = wrt_int_state%fdate(3)
      idat(3)  = wrt_int_state%fdate(1)
      idat(4)  = wrt_int_state%fdate(4)
      idat(5)  = wrt_int_state%fdate(5)
!
      if(mype==0) print *,'idat=',idat,'sdat=',sdat,'ihrst=',ihrst
!      CALL W3DIFDAT(JDATE,IDATE,0,RINC)
!
!      if(mype==0)print *,' rinc=',rinc
!      ifhr = nint(rinc(2)+rinc(1)*24.)
!      if(mype==0)print *,' ifhr=',ifhr
!      ifmin = nint(rinc(3))
!      if(ifhr /= nint(fhour))print*,'find wrong Grib file';stop
!      if(mype==0)print*,' in INITPOST ifhr ifmin =',ifhr,ifmin
!
      tstart = 0.
!
!** initialize cloud water and ice mixing ratio
!$omp parallel do default(none),private(i,j,l),shared(lm,jsta,jend,im), &
!$omp& shared(qqw,qqr,qqs,qqi)
      do l = 1,lm
        do j = jsta, jend
          do i = 1,im
            qqw(i,j,l) = 0.
            qqr(i,j,l) = 0.
            qqs(i,j,l) = 0.
            qqi(i,j,l) = 0.
          enddo
        enddo
      enddo
!
!-----------------------------------------------------------------------------
! get post fields
!-----------------------------------------------------------------------------
!
     foundland = .false.
     foundice = .false.
     get_lsmsk: do ibdl=1, wrt_int_state%FBCount

! find lans sea mask
        found = .false.
        call ESMF_FieldBundleGet(wrt_int_state%wrtFB(ibdl),fieldName='land',isPresent=found, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
!        if(mype==0) print *,'ibdl=',ibdl,'land, found=',found
        if (found) then
          call ESMF_FieldBundleGet(wrt_int_state%wrtFB(ibdl),'land',field=theField, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return  ! bail out
          call ESMF_FieldGet(theField, localDe=0, farrayPtr=arrayr42d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return  ! bail out
          ista = lbound(arrayr42d,1)
          iend = ubound(arrayr42d,1)
          !$omp parallel do default(none),private(i,j),shared(jsta,jend,ista,iend,spval,arrayr42d,sm)
          do j=jsta, jend
            do i=ista, iend
              if (arrayr42d(i,j) /= spval) sm(i,j) = 1.- arrayr42d(i,j)
            enddo
          enddo
          foundland = .true.
        endif

! find ice fraction
        found = .false.
        call ESMF_FieldBundleGet(wrt_int_state%wrtFB(ibdl),'icec',isPresent=found, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
!        if(mype==0) print *,'ibdl=',ibdl,'ice, found=',found
        if (found) then
          call ESMF_FieldBundleGet(wrt_int_state%wrtFB(ibdl),'icec',field=theField, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return  ! bail out
          call ESMF_FieldGet(theField, localDe=0, farrayPtr=arrayr42d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return  ! bail out
          ista = lbound(arrayr42d,1)
          iend = ubound(arrayr42d,1)
          !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sice,arrayr42d,sm)
          do j=jsta, jend
            do i=ista, iend
              sice(i,j) = arrayr42d(i,j)
              if (sm(i,j) /= spval .and. sm(i,j) == 0.0) sice(i,j) = 0.0
            enddo
          enddo
          foundice = .true.
        endif

     enddo get_lsmsk
     if (.not.foundland .or. .not.foundice) then
       rc=999
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return  ! bail out
     endif
     if(mype==0) print *,'after find sm and sice,imp_physics=',imp_physics,'nbdl=',wrt_int_state%FBCount
!
     file_loop_all: do ibdl=1, wrt_int_state%FBCount
!
! get grid dimension count
!       if(mype==0) print *,'in setvar, read field, ibdl=',ibdl,'idim=',   &
!         ista,iend,'jdim=',jsta,jend
       call ESMF_FieldBundleGet(wrt_int_state%wrtFB(ibdl), grid=wrtGrid,  &
         fieldCount=ncount_field, name=wrtFBName,rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
       call ESMF_GridGet(wrtgrid, dimCount=gridDimCount, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out

!       if(mype==0) print *,'in setvar, allocate fcstField,ibdl=',ibdl,'count=',ncount_field,'wrtFBname=',trim(wrtFBName)
       allocate(fcstField(ncount_field))
       call ESMF_FieldBundleGet(wrt_int_state%wrtFB(ibdl),           &
         fieldList=fcstField, itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out

!       if(mype==0) print *,'in setvar, read field, ibdl=',ibdl, 'nfield=',ncount_field
       do n=1, ncount_field
!
          call ESMF_FieldGet(fcstField(n),typekind=typekind, name=fieldname, &
            dimCount=fieldDimCount,rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return  ! bail out
          if (index(trim(fieldname),"vector") >0) cycle
!
!** for 2D fields
          if (fieldDimCount == 2) then

            if (typekind == ESMF_TYPEKIND_R4) then
              call ESMF_FieldGet(fcstField(n), localDe=0, farrayPtr=arrayr42d, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out
            else if (typekind == ESMF_TYPEKIND_R8) then
              call ESMF_FieldGet(fcstField(n), localDe=0, farrayPtr=arrayr82d, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out
              allocate( arrayr42d(ista:iend,jsta:jend))
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,arrayr42d,arrayr82d)
              do j=jsta, jend
                do i=ista, iend
                  arrayr42d(i,j) = arrayr82d(i,j)
                enddo
              enddo
            endif

            ! Terrain height (*G later)
            if(trim(fieldname)=='hgtsfc') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,fis,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  fis(i,j)=arrayr42d(i,j)
                enddo
              enddo
            endif

            ! Surface pressure
!            if(trim(fieldname)=='pressfc') then
!              !$omp parallel do private(i,j)
!              do j=jsta,jend
!                do i=ista, iend
!                  pint(i,j,lp1)=arrayr42d(i,j)
!                enddo
!              enddo
!            endif

            ! PBL height using nemsio
            if(trim(fieldname)=='hpbl') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,pblh,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  pblh(i,j)=arrayr42d(i,j)
                enddo
              enddo
            endif

            ! frictional velocity
            if(trim(fieldname)=='fricv') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,ustar,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  ustar(i,j)=arrayr42d(i,j)
                enddo
              enddo
            endif

            ! roughness length
            if(trim(fieldname)=='sfcr') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,z0,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  z0(i,j)=arrayr42d(i,j)
                enddo
              enddo
            endif

            ! sfc exchange coeff
            if(trim(fieldname)=='sfexc') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,sfcexc,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  sfcexc(i,j)=arrayr42d(i,j)
                enddo
              enddo
            endif

            ! aerodynamic conductance
            if(trim(fieldname)=='acond') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,acond,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  acond(i,j)=arrayr42d(i,j)
                enddo
              enddo
            endif

            ! surface potential T
            if(trim(fieldname)=='tmpsfc') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,arrayr42d,ths)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) /= spval) then
                    ths(i,j) = arrayr42d(i,j)
                  endif
                enddo
              enddo
            endif

            ! convective precip in m per physics time step
            if(trim(fieldname)=='cpratb_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,dtq2,arrayr42d,avgcprate)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) /= spval)   &
                    avgcprate(i,j) = arrayr42d(i,j) * (dtq2*0.001)
                enddo
              enddo
            endif

            ! continuous bucket convective precip in m per physics time step
            if(trim(fieldname)=='cprat_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,dtq2,arrayr42d,avgcprate_cont)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) /= spval) then
                    avgcprate_cont(i,j) = arrayr42d(i,j) * (dtq2*0.001)
                  endif
                enddo
              enddo
            endif

            ! time averaged bucketed precip rate
            if(trim(fieldname)=='prateb_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,dtq2,arrayr42d,avgprec)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) /= spval) then
                    avgprec(i,j) = arrayr42d(i,j) * (dtq2*0.001)
                  endif
                enddo
              enddo
            endif

            ! time averaged continuous precip rate in m per physics time step
            if(trim(fieldname)=='prate_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,dtq2,arrayr42d,avgprec_cont)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) /= spval) then
                    avgprec_cont(i,j) = arrayr42d(i,j) * (dtq2*0.001)
                  endif
                enddo
              enddo
            endif

            ! precip rate in m per physics time step
            if(trim(fieldname)=='tprcp') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,dtq2,dtp,arrayr42d,prec)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) /= spval) then
                    prec(i,j) = arrayr42d(i,j) * (dtq2*0.001) * 1000./dtp
                  endif
                enddo
              enddo
            endif

            ! convective precip rate in m per physics time step
            if(trim(fieldname)=='cnvprcp') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,dtq2,dtp,arrayr42d,cprate)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) /= spval) then
                    cprate(i,j) = max(0.,arrayr42d(i,j)) * (dtq2*0.001) * 1000./dtp
                  endif
                enddo
              enddo
            endif

            ! inst snow water eqivalent
            if(trim(fieldname)=='weasd') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sno,arrayr42d,sm,sice)
              do j=jsta,jend
                do i=ista, iend
                  sno(i,j) = arrayr42d(i,j)
                  if (sm(i,j) == 1.0 .and. sice(i,j)==0.)sno(i,j) = spval
                enddo
              enddo
            endif

            ! ave snow cover
            if(trim(fieldname)=='snowc_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,snoavg,arrayr42d,sm,sice)
              do j=jsta,jend
                do i=ista, iend
                  snoavg(i,j) = arrayr42d(i,j)
                  if (sm(i,j)==1.0 .and. sice(i,j)==0.) snoavg(i,j) = spval
                  if (snoavg(i,j) /= spval) snoavg(i,j) = snoavg(i,j)/100.
                enddo
              enddo
            endif

            ! snow depth in mm
            if(trim(fieldname)=='snod') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,si,arrayr42d,sm,sice)
              do j=jsta,jend
                do i=ista, iend
                  si(i,j) = arrayr42d(i,j)
                  if (sm(i,j)==1.0 .and. sice(i,j)==0.) si(i,j)=spval
                  if (si(i,j) /= spval) si(i,j) = si(i,j) * 1000.0
                enddo
              enddo
            endif

            ! 2m potential T (computed later)
            if(trim(fieldname)=='tmp2m') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,tshltr,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  tshltr(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! surface potential T
            if(trim(fieldname)=='spfh2m') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,qshltr,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  qshltr(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! mid day avg albedo in fraction
            if(trim(fieldname)=='albdo_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgalbedo,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  avgalbedo(i,j) = arrayr42d(i,j)
                  if (arrayr42d(i,j) /= spval) then
                    avgalbedo(i,j) = avgalbedo(i,j) * 0.01
                  endif
                enddo
              enddo
            endif

            ! time averaged column cloud fraction
            if(trim(fieldname)=='tcdc_aveclm') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgtcdc,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  avgtcdc(i,j) = arrayr42d(i,j)
                  if (arrayr42d(i,j) /= spval) then
                    avgtcdc(i,j) = avgtcdc(i,j) * 0.01
                  endif
                enddo
              enddo
            endif

            ! maximum snow albedo in fraction
            if(trim(fieldname)=='snoalb') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,mxsnal,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  mxsnal(i,j) = arrayr42d(i,j)
                  if (arrayr42d(i,j) /= spval) then
                    mxsnal(i,j) = mxsnal(i,j) * 0.01
                  endif
                enddo
              enddo
            endif

            ! ave high cloud fraction
            if(trim(fieldname)=='tcdc_avehcl') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgcfrach,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  avgcfrach(i,j) = arrayr42d(i,j)
                  if (arrayr42d(i,j) /= spval) then
                    avgcfrach(i,j) = avgcfrach(i,j) * 0.01
                  endif
                enddo
              enddo
            endif

            ! ave low cloud fraction
            if(trim(fieldname)=='tcdc_avelcl') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgcfracl,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  avgcfracl(i,j) = arrayr42d(i,j)
                  if (arrayr42d(i,j) /= spval) then
                    avgcfracl(i,j) = avgcfracl(i,j) * 0.01
                  endif
                enddo
              enddo
            endif

            ! ave middle cloud fraction
            if(trim(fieldname)=='tcdc_avemcl') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgcfracm,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  avgcfracm(i,j) = arrayr42d(i,j)
                  if (arrayr42d(i,j) /= spval) then
                    avgcfracm(i,j) = avgcfracm(i,j) * 0.01
                  endif
                enddo
              enddo
            endif

            ! inst convective cloud fraction
            if(trim(fieldname)=='tcdccnvcl') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,cnvcfr,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  cnvcfr(i,j) = arrayr42d(i,j)
                  if (arrayr42d(i,j) /= spval) then
                    cnvcfr(i,j) = cnvcfr(i,j) * 0.01
                  endif
                enddo
              enddo
            endif

            ! slope type
            if(trim(fieldname)=='sltyp') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,arrayr42d,islope)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) < spval) then
                    islope(i,j) = nint(arrayr42d(i,j))
                  else
                    islope(i,j) = 0
                  endif
                enddo
              enddo
            endif

            ! time averaged column cloud fraction
            if(trim(fieldname)=='cnwat') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,cmc,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  cmc(i,j) = arrayr42d(i,j)
                  if (arrayr42d(i,j) /= spval) cmc(i,j) = cmc(i,j) * 0.001
                  if (sm(i,j) /= 0.0) cmc(i,j) = spval
                enddo
              enddo
            endif

            ! frozen precip fraction
            if(trim(fieldname)=='cpofp') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,arrayr42d,sr)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) /= spval) then
                  !set range within (0,1)
                    sr(i,j) = min(1.,max(0.,arrayr42d(i,j)))
                  else
                    sr(i,j) = spval
                  endif
                enddo
              enddo
            endif

            ! sea ice skin temperature
            if(trim(fieldname)=='tisfc') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sice,arrayr42d,ti)
              do j=jsta,jend
                do i=ista,iend
                  if (arrayr42d(i,j) /= spval) then
                    ti(i,j) = arrayr42d(i,j)
                    if (sice(i,j) == spval .or. sice(i,j) == 0.) ti(i,j)=spval
                  else
                    ti(i,j) = spval
                  endif
                enddo
              enddo
            endif

            ! vegetation fraction
            if(trim(fieldname)=='veg') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,vegfrc,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  vegfrc(i,j) = arrayr42d(i,j)
                  if (arrayr42d(i,j) /= spval) then
                    vegfrc(i,j) = vegfrc(i,j) * 0.01
                  else
                    vegfrc(i,j) = 0.0
                  endif
                  if (sm(i,j) /= 0.0) vegfrc(i,j) = spval
                enddo
              enddo
            endif

            ! liquid volumetric soil mpisture in fraction
            if(trim(fieldname)=='soill1') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sh2o,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  sh2o(i,j,1) = arrayr42d(i,j)
                  if (sm(i,j) /= 0.0) sh2o(i,j,1) = spval
                enddo
              enddo
            endif

            ! liquid volumetric soil mpisture in fraction
            if(trim(fieldname)=='soill2') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sh2o,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  sh2o(i,j,2) = arrayr42d(i,j)
                  if (sm(i,j) /= 0.0) sh2o(i,j,2) = spval
                enddo
              enddo
            endif

            ! liquid volumetric soil mpisture in fraction
            if(trim(fieldname)=='soill3') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sh2o,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  sh2o(i,j,3) = arrayr42d(i,j)
                  if (sm(i,j) /= 0.0) sh2o(i,j,3) = spval
                enddo
              enddo
            endif

            ! liquid volumetric soil mpisture in fraction
            if(trim(fieldname)=='soill4') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sh2o,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  sh2o(i,j,4) = arrayr42d(i,j)
                  if (sm(i,j) /= 0.0) sh2o(i,j,4) = spval
                enddo
              enddo
            endif

            ! volumetric soil moisture
            if(trim(fieldname)=='soilw1') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,smc,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  smc(i,j,1) = arrayr42d(i,j)
                  if (sm(i,j) /= 0.0) smc(i,j,1) = spval
                enddo
              enddo
            endif

            ! volumetric soil moisture
            if(trim(fieldname)=='soilw2') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,smc,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  smc(i,j,2) = arrayr42d(i,j)
                  if (sm(i,j) /= 0.0) smc(i,j,2) = spval
                enddo
              enddo
            endif

            ! volumetric soil moisture
            if(trim(fieldname)=='soilw3') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,smc,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  smc(i,j,3) = arrayr42d(i,j)
                  if (sm(i,j) /= 0.0) smc(i,j,3) = spval
                enddo
              enddo
            endif

            ! volumetric soil moisture
            if(trim(fieldname)=='soilw4') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,smc,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  smc(i,j,4) = arrayr42d(i,j)
                  if (sm(i,j) /= 0.0) smc(i,j,4) = spval
                enddo
              enddo
            endif

            ! soil temperature
            if(trim(fieldname)=='soilt1') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,stc,arrayr42d,sm,sice)
              do j=jsta,jend
                do i=ista, iend
                  stc(i,j,1) = arrayr42d(i,j)
                  !mask open water areas, combine with sea ice tmp
                  if (sm(i,j) /= 0.0 .and. sice(i,j) ==0.) stc(i,j,1) = spval
                enddo
              enddo
            endif

            ! soil temperature
            if(trim(fieldname)=='soilt2') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,stc,arrayr42d,sm,sice)
              do j=jsta,jend
                do i=ista, iend
                  stc(i,j,2) = arrayr42d(i,j)
                  !mask open water areas, combine with sea ice tmp
                  if (sm(i,j) /= 0.0 .and. sice(i,j) ==0.) stc(i,j,2) = spval
                enddo
              enddo
            endif

            ! soil temperature
            if(trim(fieldname)=='soilt3') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,stc,arrayr42d,sm,sice)
              do j=jsta,jend
                do i=ista, iend
                  stc(i,j,3) = arrayr42d(i,j)
                  !mask open water areas, combine with sea ice tmp
                  if (sm(i,j) /= 0.0 .and. sice(i,j) ==0.) stc(i,j,3) = spval
                enddo
              enddo
            endif

            ! soil temperature
            if(trim(fieldname)=='soilt4') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,stc,arrayr42d,sm,sice)
              do j=jsta,jend
                do i=ista, iend
                  stc(i,j,4) = arrayr42d(i,j)
                  !mask open water areas, combine with sea ice tmp
                  if (sm(i,j) /= 0.0 .and. sice(i,j) ==0.) stc(i,j,4) = spval
                enddo
              enddo
            endif

            ! time averaged incoming sfc longwave
            if(trim(fieldname)=='dlwrf_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,alwin,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  alwin(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! inst incoming sfc longwave
            if(trim(fieldname)=='dlwrf') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,rlwin,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  rlwin(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged outgoing sfc longwave, CLDRAD puts a minus sign
            if(trim(fieldname)=='ulwrf_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,alwout,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  alwout(i,j) = arrayr42d(i,j)
                  if (alwout(i,j) /= spval) alwout(i,j) = -alwout(i,j)
                enddo
              enddo
            endif

            ! inst outgoing sfc longwave
            if(trim(fieldname)=='ulwrf') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,radot,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  radot(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged outgoing model top longwave
            if(trim(fieldname)=='ulwrf_avetoa') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,alwtoa,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  alwtoa(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged incoming sfc shortwave
            if(trim(fieldname)=='dswrf_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,aswin,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  aswin(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! inst incoming sfc shortwave
            if(trim(fieldname)=='dswrf') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,rswin,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  rswin(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged incoming sfc uv-b
            if(trim(fieldname)=='duvb_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,auvbin,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  auvbin(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged incoming sfc clear sky uv-b
            if(trim(fieldname)=='cduvb_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,auvbinc,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  auvbinc(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged outgoing sfc shortwave,CLDRAD puts a minus sign
            if(trim(fieldname)=='uswrf_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,aswout,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  aswout(i,j) = arrayr42d(i,j)
                  if (aswout(i,j) /= spval) aswout(i,j) = -aswout(i,j)
                enddo
              enddo
            endif

            ! inst outgoing sfc shortwave
            if(trim(fieldname)=='uswrf') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,rswout,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  rswout(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged model top incoming shortwave
            if(trim(fieldname)=='dswrf_avetoa') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,aswintoa,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  aswintoa(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! ime averaged model top outgoing shortwave
            if(trim(fieldname)=='uswrf_avetoa') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,aswtoa,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  aswtoa(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged surface sensible heat flux, multiplied by -1 because
            ! wrf model fluxhas reversed sign convention using gfsio
            if(trim(fieldname)=='shtfl_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sfcshx,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  sfcshx(i,j) = arrayr42d(i,j)
                  if (sfcshx(i,j) /= spval) sfcshx(i,j) = -sfcshx(i,j)
                enddo
              enddo
            endif

            ! inst surface sensible heat flux
            if(trim(fieldname)=='shtfl') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,twbs,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  twbs(i,j) = arrayr42d(i,j)
                  if (twbs(i,j) /= spval) twbs(i,j) = -twbs(i,j)
                enddo
              enddo
            endif

            ! time averaged surface latent heat flux, multiplied by -1 because
            ! wrf model flux has reversed sign vonvention using gfsio
            if(trim(fieldname)=='lhtfl_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sfclhx,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  sfclhx(i,j) = arrayr42d(i,j)
                  if (sfclhx(i,j) /= spval) sfclhx(i,j) = -sfclhx(i,j)
                enddo
              enddo
            endif

            ! inst surface latent heat flux
            if(trim(fieldname)=='lhtfl') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,qwbs,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  qwbs(i,j) = arrayr42d(i,j)
                  if (qwbs(i,j) /= spval) qwbs(i,j) = -qwbs(i,j)
                enddo
              enddo
            endif

            ! time averaged ground heat flux
            if(trim(fieldname)=='gflux_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,subshx,arrayr42d,sm,sice)
              do j=jsta,jend
                do i=ista, iend
                  subshx(i,j) = arrayr42d(i,j)
                  if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) subshx(i,j) = spval
                enddo
              enddo
            endif

            ! inst ground heat flux
            if(trim(fieldname)=='gflux') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,grnflx,arrayr42d,sm,sice)
              do j=jsta,jend
                do i=ista, iend
                  grnflx(i,j) = arrayr42d(i,j)
                  if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) grnflx(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged zonal momentum flux
            if(trim(fieldname)=='uflx_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,sfcux,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  sfcux(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged meridional momentum flux
            if(trim(fieldname)=='vflx_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,sfcvx,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  sfcvx(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged zonal gravity wave stress
            if(trim(fieldname)=='u-gwd_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,gtaux,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  gtaux(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged meridional gravity wave stress
            if(trim(fieldname)=='v-gwd_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,gtauy,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  gtauy(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged accumulated potential evaporation
            if(trim(fieldname)=='pevpr_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgpotevp,arrayr42d,sm,sice)
              do j=jsta,jend
                do i=ista, iend
                  avgpotevp(i,j) = arrayr42d(i,j)
                  if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) avgpotevp(i,j) = spval
                enddo
              enddo
            endif

            ! inst potential evaporation
            if(trim(fieldname)=='pevpr') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,potevp,arrayr42d,sm,sice)
              do j=jsta,jend
                do i=ista, iend
                  potevp(i,j) = arrayr42d(i,j)
                  if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) potevp(i,j) = spval
                enddo
              enddo
            endif

            ! 10 m u
            if(trim(fieldname)=='ugrd10m') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,u10,arrayr42d,u10h)
              do j=jsta,jend
                do i=ista, iend
                  u10(i,j) = arrayr42d(i,j)
                  u10h(i,j) = u10(i,j)
                enddo
              enddo
            endif

            ! 10 m v
            if(trim(fieldname)=='vgrd10m') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,v10,arrayr42d,v10h)
              do j=jsta,jend
                do i=ista, iend
                  v10(i,j) = arrayr42d(i,j)
                  v10h(i,j) = v10(i,j)
                enddo
              enddo
            endif

            ! vegetation type
            if(trim(fieldname)=='vtype') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,arrayr42d,ivgtyp)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) < spval) then
                    ivgtyp(i,j) = nint(arrayr42d(i,j))
                  else
                    ivgtyp(i,j) = 0
                  endif
                enddo
              enddo
            endif

            ! soil type
            if(trim(fieldname)=='sotyp') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,arrayr42d,isltyp)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) < spval) then
                    isltyp(i,j) = nint(arrayr42d(i,j))
                  else
                    isltyp(i,j) = 0
                  endif
                enddo
              enddo
            endif

            ! inst  cloud top pressure
            if(trim(fieldname)=='prescnvclt') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,ptop,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  ptop(i,j) = arrayr42d(i,j)
                  if(ptop(i,j) <= 0.0) ptop(i,j) = spval
                enddo
              enddo
            endif

            ! inst cloud bottom pressure
            if(trim(fieldname)=='prescnvclb') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,pbot,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  pbot(i,j) = arrayr42d(i,j)
                  if(pbot(i,j) <= 0.0) pbot(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged low cloud top pressure
            if(trim(fieldname)=='pres_avelct') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,ptopl,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  ptopl(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged low cloud bottom pressure
            if(trim(fieldname)=='pres_avelcb') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,pbotl,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  pbotl(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged low cloud top temperature
            if(trim(fieldname)=='tmp_avelct') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,ttopl,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  ttopl(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged middle cloud top pressure
            if(trim(fieldname)=='pres_avemct') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,ptopm,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  ptopm(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged middle cloud bottom pressure
            if(trim(fieldname)=='pres_avemcb') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,pbotm,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  pbotm(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged middle cloud top temperature
            if(trim(fieldname)=='tmp_avemct') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,ttopm,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  ttopm(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged high cloud top pressure
            if(trim(fieldname)=='pres_avehct') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,ptoph,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  ptoph(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged high cloud bottom pressure
            if(trim(fieldname)=='pres_avehcb') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,pboth,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  pboth(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged high cloud top temperature
            if(trim(fieldname)=='tmp_avehct') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,ttoph,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  ttoph(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged boundary layer cloud cover
            if(trim(fieldname)=='tcdc_avebndcl') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,pblcfr,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  pblcfr(i,j) = arrayr42d(i,j)
                  if (pblcfr(i,j) < spval) pblcfr(i,j) = pblcfr(i,j) * 0.01
                enddo
              enddo
            endif

            ! cloud work function
            if(trim(fieldname)=='cwork_aveclm') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,cldwork,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  cldwork(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! water runoff
            if(trim(fieldname)=='watr_acc') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,runoff,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  runoff(i,j) = arrayr42d(i,j)
                  if (sm(i,j) /= 0.0) runoff(i,j) = spval
                enddo
              enddo
            endif

            ! shelter max temperature
            if(trim(fieldname)=='tmax_max2m') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,maxtshltr,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  maxtshltr(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! shelter min temperature
            if(trim(fieldname)=='tmin_min2m') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,mintshltr,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  mintshltr(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! ice thickness
            if(trim(fieldname)=='icetk') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,dzice,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  dzice(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! wilting point
            if(trim(fieldname)=='wilt') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend, spval,smcwlt,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  smcwlt(i,j) = arrayr42d(i,j)
                  if (sm(i,j) /= 0.0) smcwlt(i,j) = spval
                enddo
              enddo
            endif

            ! sunshine duration
            if(trim(fieldname)=='sunsd_acc') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,suntime,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  suntime(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! field capacity
            if(trim(fieldname)=='fldcp') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,fieldcapa,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  fieldcapa(i,j) = arrayr42d(i,j)
                  if (sm(i,j) /= 0.0) fieldcapa(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged surface visible beam downward solar flux
            if(trim(fieldname)=='vbdsf_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,avisbeamswin,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  avisbeamswin(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged surface visible diffuse downward solar flux
            if(trim(fieldname)=='vddsf_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,avisdiffswin,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  avisdiffswin(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged surface near IR beam downward solar flux
            if(trim(fieldname)=='nbdsf_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,airbeamswin,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  airbeamswin(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged surface near IR diffuse downward solar flux
            if(trim(fieldname)=='nddsf_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,airdiffswin,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  airdiffswin(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged surface clear sky outgoing LW
            if(trim(fieldname)=='csulf') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,alwoutc,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  alwoutc(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged TOA clear sky outgoing LW
            if(trim(fieldname)=='csulftoa') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,alwtoac,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  alwtoac(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged surface clear sky outgoing SW
            if(trim(fieldname)=='csusf') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,aswoutc,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  aswoutc(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged TOA clear sky outgoing SW
            if(trim(fieldname)=='csusftoa') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,aswtoac,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  aswtoac(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged surface clear sky incoming LW
            if(trim(fieldname)=='csdlf') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,alwinc,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  alwinc(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! time averaged surface clear sky incoming SW
            if(trim(fieldname)=='csdsf') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,aswinc,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  aswinc(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! shelter max specific humidity
            if(trim(fieldname)=='spfhmax_max2m') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,maxqshltr,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  maxqshltr(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! shelter min temperature
            if(trim(fieldname)=='spfhmin_min2m') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,minqshltr,arrayr42d)
              do j=jsta,jend
                do i=ista, iend
                  minqshltr(i,j) = arrayr42d(i,j)
                enddo
              enddo
            endif

            ! storm runoffs
            if(trim(fieldname)=='ssrun_acc') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,ssroff,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  ssroff(i,j) = arrayr42d(i,j)
                  if (sm(i,j) /= 0.0) ssroff(i,j) = spval
                enddo
              enddo
            endif

            !  direct soil evaporation
            if(trim(fieldname)=='evbs_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgedir,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  avgedir(i,j) = arrayr42d(i,j)
                  if (sm(i,j) /= 0.0) avgedir(i,j) = spval
                enddo
              enddo
            endif

            ! canopy water evap
            if(trim(fieldname)=='evcw_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgecan,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  avgecan(i,j) = arrayr42d(i,j)
                  if (sm(i,j) /= 0.0) avgecan(i,j) = spval
                enddo
              enddo
            endif

            ! plant transpiration
            if(trim(fieldname)=='trans_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgetrans,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  avgetrans(i,j) = arrayr42d(i,j)
                  if (sm(i,j) /= 0.0) avgetrans(i,j) = spval
                enddo
              enddo
            endif

            ! snow sublimation
            if(trim(fieldname)=='sbsno_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgesnow,arrayr42d,sm,sice)
              do j=jsta,jend
                do i=ista, iend
                  avgesnow(i,j) = arrayr42d(i,j)
                  if (sm(i,j)==1.0 .and. sice(i,j)==0.) avgesnow(i,j)=spval
                enddo
              enddo
            endif

            ! total soil moisture
            if(trim(fieldname)=='soilm') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,smstot,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  smstot(i,j) = arrayr42d(i,j)
                  if (sm(i,j) /= 0.0) smstot(i,j) = spval
                enddo
              enddo
            endif

            ! snow phase change heat flux
            if(trim(fieldname)=='snohf') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,snopcx,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  snopcx(i,j) = arrayr42d(i,j)
                  if (sm(i,j) /= 0.0) snopcx(i,j) = spval
                enddo
              enddo
            endif

!          else if (fieldDimCount > gridDimCount) then
          else if (fieldDimCount ==3) then
            if (typekind == ESMF_TYPEKIND_R4) then
              call ESMF_FieldGet(fcstField(n), localDe=0, farrayPtr=arrayr43d, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out
            else if (typekind == ESMF_TYPEKIND_R8) then
              call ESMF_FieldGet(fcstField(n), localDe=0, farrayPtr=arrayr83d, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out
              allocate(arrayr43d(ista:iend,jsta:jend,kstart:kend))
              arrayr43d = 0.
              do k=kstart,kend
             !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,k,arrayr43d,arrayr83d)
                do j=jsta,jend
                  do i=ista,iend
                    arrayr43d(i,j,k) = arrayr83d(i,j,k)
                  enddo
                enddo
              enddo
            endif

            ! model level T
            if(trim(fieldname)=='tmp') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,t,arrayr43d)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    t(i,j,l)=arrayr43d(i,j,l)
                  enddo
                enddo
              enddo

              !! sig4
              !$omp parallel do default(none) private(i,j,tlmh) shared(lm,jsta,jend,ista,iend,t,sigt4)
              do j=jsta,jend
                do i=ista, iend
                  tlmh = t(i,j,lm) * t(i,j,lm)
                  sigt4(i,j) = 5.67E-8 * tlmh * tlmh
                enddo
              enddo
            endif

            ! model level spfh
            if(trim(fieldname)=='spfh') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,q,arrayr43d)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    q(i,j,l)=arrayr43d(i,j,l)
                  enddo
                enddo
              enddo
            endif

            ! model level u wind
            if(trim(fieldname)=='ugrd') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,uh,arrayr43d)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    uh(i,j,l)=arrayr43d(i,j,l)
                  enddo
                enddo
              enddo
            endif

            ! model level v wind
            if(trim(fieldname)=='vgrd') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,vh,arrayr43d)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    vh(i,j,l)=arrayr43d(i,j,l)
                  enddo
                enddo
              enddo
            endif

            ! model level pressure thinkness
            if(trim(fieldname)=='dpres') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,dpres,arrayr43d)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    dpres(i,j,l)=arrayr43d(i,j,l)
                  enddo
                enddo
              enddo
            endif

            ! model level gh thinkness, model output negative delz
            if(trim(fieldname)=='delz') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,zint,arrayr43d)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    zint(i,j,l)=-1.*arrayr43d(i,j,l)
                  enddo
                enddo
              enddo
            endif

            ! model level w
            if(trim(fieldname)=='dzdt') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,wh,arrayr43d)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    wh(i,j,l)=arrayr43d(i,j,l)
                  enddo
                enddo
              enddo
            endif

            ! model level ozone mixing ratio
#ifdef MULTI_GASES
            if(trim(fieldname)=='spo3') then
#else
            if(trim(fieldname)=='o3mr') then
#endif
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,o3,arrayr43d)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    o3(i,j,l)=arrayr43d(i,j,l)
                  enddo
                enddo
              enddo
            endif

! for GFDL MP or Thompson MP 
            if (imp_physics == 11 .or. imp_physics == 8) then
              ! model level cloud water mixing ratio
              if(trim(fieldname)=='clwmr') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqw,arrayr43d)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqw(i,j,l)=arrayr43d(i,j,l)
                    enddo
                  enddo
                enddo
              endif

              ! model level ice mixing ratio
              if(trim(fieldname)=='icmr') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqi,arrayr43d)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqi(i,j,l)=arrayr43d(i,j,l)
                    enddo
                  enddo
                enddo
              endif

              ! model level rain water mixing ratio
              if(trim(fieldname)=='rwmr') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqr,arrayr43d)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqr(i,j,l)=arrayr43d(i,j,l)
                    enddo
                  enddo
                enddo
              endif

              ! model level snow mixing ratio
              if(trim(fieldname)=='snmr') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqs,arrayr43d)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqs(i,j,l)=arrayr43d(i,j,l)
                    enddo
                  enddo
                enddo
              endif

              ! model level rain water mixing ratio
              if(trim(fieldname)=='grle') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqg,arrayr43d)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqg(i,j,l)=arrayr43d(i,j,l)
                    enddo
                  enddo
                enddo
              endif

              if(imp_physics == 8) then
              ! model level rain number
              if(trim(fieldname)=='ncrain') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqnr,arrayr43d)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqnr(i,j,l)=arrayr43d(i,j,l)
                    enddo
                  enddo
                enddo
              endif

              ! model level rain number
              if(trim(fieldname)=='ncice') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqni,arrayr43d)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqni(i,j,l)=arrayr43d(i,j,l)
                    enddo
                  enddo
                enddo
              endif

              ! model level rain number
              if(trim(fieldname)=='nwfa') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqnwfa,arrayr43d)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqnwfa(i,j,l)=arrayr43d(i,j,l)
                    enddo
                  enddo
                enddo
              endif

              ! model level rain number
              if(trim(fieldname)=='nifa') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqnifa,arrayr43d)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqnifa(i,j,l)=arrayr43d(i,j,l)
                    enddo
                  enddo
                enddo
              endif
              endif !if(imp_physics == 8) then
!gfdlmp
            endif

            ! model level cloud amount
            if(trim(fieldname)=='cld_amt') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,cfr,arrayr43d)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    cfr(i,j,l)=arrayr43d(i,j,l)
                  enddo
                enddo
              enddo
            endif

            ! model level ref3d
            if(trim(fieldname)=='ref3D') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,ref_10cm,arrayr43d)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    ref_10cm(i,j,l)=arrayr43d(i,j,l)
                  enddo
                enddo
              enddo
!              print *,'in gfs_post, get ref_10cm=',maxval(ref_10cm), minval(ref_10cm)
            endif

            ! model level ref3d
            if(trim(fieldname)=='tke') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,q2,arrayr43d)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    q2(i,j,l)=arrayr43d(i,j,l)
                  enddo
                enddo
              enddo
            endif
!3d fields
          endif

! end loop ncount_field
        enddo

        if ( index(trim(wrt_int_state%wrtFB_names(ibdl)),trim(filename_base(1))) > 0)   &
          setvar_atmfile = .true.
        if ( index(trim(wrt_int_state%wrtFB_names(ibdl)),trim(filename_base(2))) > 0)   &
          setvar_sfcfile = .true.
        deallocate(fcstField)

! end file_loop_all
      enddo file_loop_all

! recompute full layer of zint
!$omp parallel do default(none) private(i,j) shared(jsta,jend,im,lp1,spval,zint,fis)
      do j=jsta,jend
        do i=1,im
          if (fis(i,j) /= spval) then
            zint(i,j,lp1) = fis(i,j)
            fis(i,j)      = fis(i,j) * grav
          endif
        enddo
      enddo

      do l=lm,1,-1
!$omp parallel do default(none) private(i,j) shared(l,jsta,jend,im,omga,wh,dpres,zint)
        do j=jsta,jend
          do i=1,im
            omga(i,j,l) = (-1.) * wh(i,j,l) * dpres(i,j,l)/zint(i,j,l)
            zint(i,j,l) = zint(i,j,l) + zint(i,j,l+1)
          enddo
        enddo
      enddo

! compute pint from top down
!$omp parallel do default(none) private(i,j) shared(jsta,jend,im,ak5,pint)
      do j=jsta,jend
        do i=1,im
          pint(i,j,1) = ak5(1)
        end do
      end do

      do l=2,lp1
!$omp parallel do default(none) private(i,j) shared(l,jsta,jend,im,pint,dpres)
        do j=jsta,jend
          do i=1,im
            pint(i,j,l) = pint(i,j,l-1) + dpres(i,j,l-1)
          enddo
        enddo
      end do

!compute pmid from averaged two layer pint
      do l=lm,1,-1
!$omp parallel do default(none) private(i,j) shared(l,jsta,jend,im,pmid,pint)
        do j=jsta,jend
          do i=1,im
            pmid(i,j,l) = 0.5*(pint(i,j,l)+pint(i,j,l+1))
          enddo
        enddo
      enddo

!$omp parallel do default(none) private(i,j) shared(jsta,jend,im,spval,pt,pd,pint)
      do j=jsta,jend
        do i=1,im
          pd(i,j)     = spval
          pint(i,j,1) = pt
        end do
      end do
!      print *,'in setvar, pt=',pt,'ak5(lp1)=', ak5(lp1),'ak5(1)=',ak5(1)

! compute alpint
      do l=lp1,1,-1
!$omp parallel do default(none) private(i,j) shared(l,jsta,jend,im,alpint,pint)
        do j=jsta,jend
          do i=1,im
            alpint(i,j,l)=log(pint(i,j,l))
          end do
        end do
      end do

! compute zmid  
      do l=lm,1,-1
!$omp parallel do default(none) private(i,j) shared(l,jsta,jend,im,zmid,zint,pmid,alpint)
        do j=jsta,jend
          do i=1,im
            zmid(i,j,l)=zint(i,j,l+1)+(zint(i,j,l)-zint(i,j,l+1))* &
                    (log(pmid(i,j,l))-alpint(i,j,l+1))/ &
                    (alpint(i,j,l)-alpint(i,j,l+1))
          end do
        end do
      end do
!        print *,'in post_gfs,zmid=',maxval(zmid(1:im,jsta:jend,1)), &
!          minval(zmid(1:im,jsta:jend,1)),maxloc(zmid(1:im,jsta:jend,1)), &
!          'zint=',maxval(zint(1:im,jsta:jend,2)),minval(zint(1:im,jsta:jend,1)),  &
!          'pmid=',maxval(pmid(1:im,jsta:jend,1)),minval(pmid(1:im,jsta:jend,1)),  &
!          'alpint=',maxval(alpint(1:im,jsta:jend,2)),minval(alpint(1:im,jsta:jend,2))
!        print *,'in post_gfs,alpint=',maxval(alpint(1:im,jsta:jend,1)), &
!          minval(alpint(1:im,jsta:jend,1))

! surface potential T, and potential T at roughness length
!$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,lp1,sm,ths,sst,thz0,pint)
      do j=jsta,jend
        do i=ista, iend
          !assign sst
          if (sm(i,j) /= 0.0 .and. ths(i,j) /= spval) then
             sst(i,j) = ths(i,j)
          else
             sst(i,j) = spval
          endif
          if (ths(i,j) /= spval) then
            ths(i,j)  = ths(i,j)* (p1000/pint(i,j,lp1))**capa
            thz0(i,j) = ths(i,j)
          endif
        enddo
      enddo

! compute cwm for gfdlmp or Thompson 
      if(  imp_physics == 11 .or. imp_physics == 8) then
        do l=1,lm
!$omp parallel do default(none) private(i,j) shared(l,jsta,jend,ista,iend,cwm,qqg,qqs,qqr,qqi,qqw)
          do j=jsta,jend
            do i=ista,iend
              cwm(i,j,l)=qqg(i,j,l)+qqs(i,j,l)+qqr(i,j,l)+qqi(i,j,l)+qqw(i,j,l)
            enddo
          enddo
        enddo
      endif

! estimate 2m pres and convert t2m to theta
!$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,lm,pshltr,pint,tshltr)
      do j=jsta,jend
        do i=ista, iend
          pshltr(I,J)=pint(i,j,lm+1)*EXP(-0.068283/tshltr(i,j))
          tshltr(i,j)= tshltr(i,j)*(p1000/pshltr(I,J))**CAPA
        enddo
      enddo

!htop
      do j=jsta,jend
        do i=1,im
          htop(i,j) = spval
          if(ptop(i,j) < spval)then
            do l=1,lm
              if(ptop(i,j) <= pmid(i,j,l))then
                htop(i,j)=l
                exit
              end if
            end do
          end if
        end do
      end do

! hbot
      do j=jsta,jend
        do i=1,im
          hbot(i,j) = spval
          if(pbot(i,j) < spval)then
            do l=lm,1,-1
              if(pbot(i,j) >= pmid(i,j,l)) then
                hbot(i,j) = l
                exit
              end if
            end do
          end if
        end do
      end do

! generate look up table for lifted parcel calculations
      thl    = 210.
      plq    = 70000.
      pt_tbl = 10000.          ! this is for 100 hPa added by Moorthi

      call table(ptbl,ttbl,pt_tbl,                                     &
                 rdq,rdth,rdp,rdthe,pl,thl,qs0,sqs,sthe,the0)

      call tableq(ttblq,rdpq,rdtheq,plq,thl,stheq,the0q)

      if(mype == 0)then
        write(6,*)'  SPL (POSTED PRESSURE LEVELS) BELOW: '
        write(6,51) (SPL(L),L=1,LSM)
   50   format(14(F4.1,1X))
   51   format(8(F8.1,1X))
      endif
!
!$omp parallel do default(none) private(l) shared(lsm,alsl,spl)
      do l = 1,lsm
         alsl(l) = log(spl(l))
      end do
!
!      print *,'in gfs_post, end ref_10cm=',maxval(ref_10cm), minval(ref_10cm)
!!! above is fv3 change
!
!more fields need to be computed
!
    end subroutine set_postvars_gfs


    end module post_gfs
