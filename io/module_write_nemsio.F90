module module_write_nemsio

  use esmf
  use nemsio_module
 
  implicit none
 
  include 'mpif.h'

  private
  logical :: first_nemsio_call
  integer :: im,jm,lm, idate(7),nmeta, nsoil,ncld, idrt, ntrac
  integer :: mp_physi, CU_PHYSICS
  integer :: mype, ntasks, mpi_comm, nbdl
  logical :: hydrostatic
  integer,dimension(200,100)      :: nfldlev
  character(16),dimension(3000,5) :: recname,reclevtyp 
  integer,dimension(3000,5)       :: reclev 

  integer,dimension(:), allocatable    :: nrec
  integer,dimension(:), allocatable    :: fieldcount
  integer, dimension(:), allocatable :: idisp, irecv
!
  integer,dimension(:), allocatable :: nmetavari,nmetavarc, nmetavarr4,nmetavarr8
  integer,dimension(:), allocatable :: nmetaaryi,nmetaaryc, nmetaaryr4,nmetaaryr8
  character(16),dimension(:,:),allocatable :: variname, varcname, varr4name, varr8name
  integer, dimension(:,:), allocatable   :: varival
  real(4), dimension(:,:), allocatable   :: varr4val
  real(8), dimension(:,:), allocatable   :: varr8val
  character(16), dimension(:,:), allocatable   :: varcval
  logical, dimension(:), allocatable :: extrameta

  public nemsio_first_call, write_nemsio
  
  contains

  subroutine nemsio_first_call(fieldbundle, imo, jmo, &
             wrt_mype, wrt_ntasks, wrt_mpi_comm, wrt_nbdl, mybdl, inidate, rc)
    type(ESMF_FieldBundle), intent(in)     :: fieldbundle
    integer, intent(in)                    :: imo, jmo
    integer, intent(in)                    :: wrt_mype, wrt_ntasks, wrt_mpi_comm
    integer, intent(in)                    :: wrt_nbdl, mybdl
    integer, intent(in)                    :: inidate(7)
    integer, optional,intent(out)          :: rc

!** local vars
    integer i,j, nfld
    integer fieldDimCount,gridDimCount
    character(100) :: fieldname
    type(ESMF_GRID)            :: wrtgrid
    type(ESMF_TypeKind_Flag)   :: typekind

    integer, dimension(:), allocatable :: ungriddedLBound, ungriddedUBound
    type(ESMF_Field), allocatable  :: fcstField(:)
    
!-------------------------------------------------------------------
!
    im    = imo
    jm    = jmo
    nmeta = 5
    idrt  = 4
    nsoil = 4
    ntrac = 3
    ncld  = 1
    idate(1:7) = inidate(1:7)
    mype  = wrt_mype
    ntasks= wrt_ntasks
    nbdl  = wrt_nbdl
    mpi_comm = wrt_mpi_comm
    if(.not.allocated(idisp)) allocate(idisp(ntasks),irecv(ntasks))
    if(.not.allocated(fieldCount)) allocate(fieldCount(nbdl))
    if(.not.allocated(nrec)) allocate(nrec(nbdl))
!
!** get attibute info from fieldbundle
    call get_global_attr(fieldbundle, mybdl, rc=rc)

!** get meta info from fieldbundle
    call ESMF_FieldBundleGet(fieldbundle, fieldCount=fieldCount(mybdl), &
         grid=wrtGrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return  ! bail out
!
    call ESMF_GridGet(wrtgrid, dimCount=gridDimCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!
    allocate(fcstField(fieldCount(mybdl)))
    call ESMF_FieldBundleGet(fieldbundle, fieldList=fcstField,     &
         itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    nrec(mybdl)=0
    lm=1
    do i=1,fieldcount(mybdl)

       call ESMF_FieldGet(fcstField(i), typekind=typekind, &
            dimCount=fieldDimCount, grid=wrtGrid, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
       if (fieldDimCount > gridDimCount) then
         allocate(ungriddedLBound(fieldDimCount-gridDimCount))
         allocate(ungriddedUBound(fieldDimCount-gridDimCount))
         call ESMF_FieldGet(fcstField(i), ungriddedLBound=ungriddedLBound, &
                    ungriddedUBound=ungriddedUBound, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, &
           file=__FILE__)) &
           return  ! bail out
         nfldlev(i,mybdl) = ungriddedUBound(fieldDimCount-gridDimCount) - &
                        ungriddedLBound(fieldDimCount-gridDimCount) + 1
         nrec(mybdl) = nrec(mybdl) + nfldlev(i,mybdl)
         lm   = nfldlev(i,mybdl)
         deallocate(ungriddedLBound)
         deallocate(ungriddedUBound)
       else if(fieldDimCount == 2) then
         nfldlev(i,mybdl) = 1
         nrec(mybdl) = nrec(mybdl) + 1
       endif
       
    enddo
!
    nfld = 1
    do i=1,fieldcount(mybdl)
      call ESMF_FieldGet(fcstField(i),name=fieldName,rc=rc)
      if( nfldlev(i,mybdl) == 1) then
        recname(nfld,mybdl)   = trim(fieldName)
        reclevtyp(nfld,mybdl) = "sfc"
        reclev(nfld,mybdl)    = 1
        nfld = nfld + 1
      else
        do j = 1,nfldlev(i,mybdl)
          recname(nfld,mybdl)   = trim(fieldName)
          reclevtyp(nfld,mybdl) = "isobaric_sfc"
          reclev(nfld,mybdl)    = j
          nfld = nfld + 1
        enddo
      endif
    enddo
!
      
  end subroutine nemsio_first_call

!----------------------------------------------------------------------------------------
  subroutine write_nemsio(fieldbundle, filename, nf_hours, &
             nf_minutes, nf_seconds, mybdl, rc)
!
    type(ESMF_FieldBundle), intent(in)     :: fieldbundle
    character(*), intent(in)               :: filename
    integer, intent(in)                    :: nf_hours, nf_minutes, nf_seconds
    integer, intent(in)                    :: mybdl
    integer, optional,intent(out)          :: rc
!
!** local vars
    integer i,j,m,n,k, jrec
    integer istart, iend, jstart, jend, kstart, kend, nlen
    real(4),dimension(:),allocatable    :: tmp
    real(4),dimension(:,:),allocatable  :: arrayr4
    real(4),dimension(:,:),pointer      :: arrayr42d
    real(8),dimension(:,:),pointer      :: arrayr82d
    real(4),dimension(:,:,:),pointer    :: arrayr43d
    real(8),dimension(:,:,:),pointer    :: arrayr83d
    type(ESMF_Field), allocatable  :: fcstField(:)
    type(ESMF_TypeKind_Flag)       :: typekind
    type(nemsio_gfile) ::  nemsiofile
!

!** init nemsio
    call nemsio_init(iret=rc)
!
!**  OPEN NEMSIO FILE
!
!    print *,'in write_nemsio,bf nemsio_open, filename=',trim(filename), &
!      'idate=',idate,'nfour=',NF_HOURS,NF_MINUTES,NF_SECONDS, 'mybdl=',mybdl,&
!      'dim=',im,jm,lm,'nmeta=',nmeta,'idrt=',idrt,'nsoil=',nsoil, &
!      'ntrac=',ntrac,'nrec=',nrec(mybdl),'extrameta=',extrameta(mybdl), &
!      'nmetavari=',nmetavari(mybdl),'nmetavarc=',nmetavarc(mybdl)
!    if(nmetavari(mybdl)>0) print *,'in write_nemsio,bf nemsio_open,nmetavari=', &
!      nmetavari(mybdl),'varival=',trim(variname(1,mybdl)),varival(1,mybdl)
!    if(nmetavarc(mybdl)>0) print *,'in write_nemsio,bf nemsio_open,nmetavarc=', &
!      nmetavarc(mybdl),'varcval=',trim(varcname(1,mybdl)),varcval(1,mybdl)
 
    if(mype==0) then
      call nemsio_open(nemsiofile,trim(FILENAME),'write',rc,    &
        modelname="FV3", gdatatype="bin4",                      &
        idate=idate,nfhour=NF_HOURS, nfminute=NF_MINUTES,       &
        nfsecondn=NF_SECONDS*100, nfsecondd=100,          &
        dimx=im,dimy=jm,dimz=lm, nmeta=nmeta,idrt=idrt,         &
        nsoil=nsoil,ntrac=ntrac,nrec=nrec(mybdl), ncldt=ncld,          &
        extrameta=extrameta(mybdl),recname=RECNAME(1:nrec(mybdl),mybdl), &
        reclevtyp=RECLEVTYP(1:nrec(mybdl),mybdl),  &
        reclev=RECLEV(1:nrec(mybdl),mybdl),        &
        nmetavari=nmetavari(mybdl), nmetavarc=nmetavarc(mybdl),  &
        variname=variname(1:nmetavari(mybdl),mybdl), &
        varival=varival(1:nmetavari(mybdl),mybdl),   &
        varcname=varcname(1:nmetavarc(mybdl),mybdl), &
        varcval=varcval(1:nmetavarc(mybdl),mybdl) )
       if(rc/=0) print *,'nemsio_open, file=',trim(filename),' iret=',rc
    endif
!
!** collect data to first pe and write out data
!
    allocate(arrayr4(im,jm),tmp(im*jm))
    allocate(fcstField(fieldCount(mybdl)))
    call ESMF_FieldBundleGet(fieldbundle, fieldList=fcstField,     &
         itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) & return  ! bail out

    jrec = 1
    do i=1, fieldcount(mybdl)
!
       call ESMF_FieldGet(fcstField(i),typekind=typekind, rc=rc)

       if( nfldlev(i,mybdl) == 1) then
         if( typekind == ESMF_TYPEKIND_R4) then
           call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=arrayr42d, rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) return  ! bail out
           istart = lbound(arrayr42d,1)
           iend   = ubound(arrayr42d,1)
           jstart = lbound(arrayr42d,2)
           jend   = ubound(arrayr42d,2)
           nlen   = (iend-istart+1) * (jend-jstart+1)
         elseif( typekind == ESMF_TYPEKIND_R8) then
           call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=arrayr82d, rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) return  ! bail out
           istart = lbound(arrayr82d,1)
           iend   = ubound(arrayr82d,1)
           jstart = lbound(arrayr82d,2)
           jend   = ubound(arrayr82d,2)
           nlen   = (iend-istart+1) * (jend-jstart+1)
           allocate( arrayr42d(istart:iend,jstart:jend))
           do n=jstart,jend
             do m=istart,iend
               arrayr42d(m,n) = arrayr82d(m,n)
             enddo
           enddo
         endif
! send data to task 0
         call mpi_gather(nlen, 1, MPI_INTEGER, irecv(:), 1, MPI_INTEGER, 0, mpi_comm, rc)
         if(mype == 0) then
           idisp(1) = 0
           do n=1,ntasks-1
             idisp(n+1) = idisp(n) + irecv(n)
           enddo
!           if(mype==0) print *,' collect data, idisp=',idisp(:)
!           if(mype==0) print *,' collect data, irecv=',irecv(:)
         endif
!         if( trim(recname(jrec,mybdl))=="HGTsfc" .and. trim(recname(jrec,mybdl))=="sfc") then
!           print *,'in write nemsio,fb=',i,' write jrec=',jrec,' val=',maxval(arrayr42d(istart:iend,jstart:jend)),  &
!             minval(arrayr42d(istart:iend,jstart:jend)),maxloc(arrayr42d(istart:iend,jstart:jend)), &
!             minloc(arrayr42d(istart:iend,jstart:jend))
!         endif
         call mpi_gatherv(arrayr42d,nlen,MPI_REAL, arrayr4,irecv,idisp(:), MPI_REAL, &
              0, mpi_comm, rc)
         if(mype==0) then
!           print *,'in write nemsio, value=',maxval(arrayr4(1:im,1:jm)), &
!              minval(arrayr4(1:im,1:jm)),maxloc(arrayr4(1:im,1:jm)),minloc(arrayr4(1:im,1:jm))
           tmp = reshape(arrayr4, (/im*jm/))
           call nemsio_writerec(nemsiofile, jrec, tmp, iret=rc)
!           print *,'in write nemsio,fb=',i,' write jrec=',jrec,'fld is',  &
!             trim(recname(jrec,mybdl)), 'rc=', &
!             rc, 'value=',maxval(tmp),minval(tmp),maxloc(tmp),minloc(tmp)
         endif
         jrec = jrec + 1

      elseif (nfldlev(i,mybdl) > 1) then

         if( typekind == ESMF_TYPEKIND_R4) then
           call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=arrayr43d, rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) return  ! bail out

           istart = lbound(arrayr43d,1)
           iend   = ubound(arrayr43d,1)
           jstart = lbound(arrayr43d,2)
           jend   = ubound(arrayr43d,2)
           kstart = lbound(arrayr43d,3)
           kend   = ubound(arrayr43d,3)
           nlen   = (iend-istart+1) * (jend-jstart+1)
           lm     = kend - kstart + 1
         elseif( typekind == ESMF_TYPEKIND_R8) then
           call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=arrayr83d, rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) return  ! bail out
           istart = lbound(arrayr83d,1)
           iend   = ubound(arrayr83d,1)
           jstart = lbound(arrayr83d,2)
           jend   = ubound(arrayr83d,2)
           kstart = lbound(arrayr43d,3)
           kend   = ubound(arrayr43d,3)
           nlen   = (iend-istart+1) * (jend-jstart+1)
           lm     = kend - kstart + 1
         endif

         ! send data to task 0
         call mpi_gather(nlen, 1, MPI_INTEGER, irecv, 1, MPI_INTEGER, 0, mpi_comm, rc)
         if(mype == 0) then
           idisp(1) = 0
           do n=1,ntasks-1
             idisp(n+1) = idisp(n) + irecv(n)
           enddo
         endif
! write out all levels
         allocate(arrayr42d(istart:iend,jstart:jend))
         do k=kstart,kend
           if (typekind == ESMF_TYPEKIND_R4) then
             do n=jstart,jend
               do m=istart,iend
                  arrayr42d(m,n)=arrayr43d(m,n,k)    
               enddo
             enddo
           elseif (typekind == ESMF_TYPEKIND_R8) then
             do n=jstart,jend
               do m=istart,iend
                  arrayr42d(m,n)=arrayr83d(m,n,k)    
               enddo
             enddo
           endif
!
           call mpi_gatherv(arrayr42d, nlen, MPI_REAL, arrayr4, irecv,idisp, MPI_REAL, &
              0, mpi_comm, rc)
           if(mype==0) then
             tmp = reshape(arrayr4, (/im*jm/))
             call nemsio_writerec(nemsiofile, jrec, tmp, iret=rc)
             jrec = jrec + 1
           endif
         enddo
         deallocate(arrayr42d)
!
      endif
    enddo
!
    deallocate(tmp)
    deallocate(arrayr4)
    deallocate(fcstField)
!
!** close nemsio file
    call nemsio_close(nemsiofile, iret=rc)
!
    call nemsio_finalize()

  end subroutine write_nemsio

!----------------------------------------------------------------------------------------

  subroutine write_nemaio_final()

!**
    deallocate(irecv)
    deallocate(idisp)
    deallocate(fieldcount)

  end subroutine write_nemaio_final
!
!----------------------------------------------------------------------------------------

  subroutine get_global_attr(fldbundle, mybdl, rc)
    type(ESMF_FieldBundle), intent(in) :: fldbundle
    integer, intent(in)                :: mybdl
    integer, intent(out)               :: rc

! local variable
   integer i,j, k,n, attcount
   integer ni,nr4,nr8, nc
   character(80) attName, hydrostatics, fldname
   type(ESMF_TypeKind_Flag)                :: typekind
!
! look at the field bundle attributes 
    call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
      attnestflag=ESMF_ATTNEST_OFF, Count=attcount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

! first loop over all the attributes to find the count for integer attr, real
! attr, etc
    j=1
    k=1
    if (.not. allocated(nmetavari))  then
      allocate(nmetavari(nbdl),nmetavarr4(nbdl),nmetavarr8(nbdl),nmetavarc(nbdl))
      allocate(nmetaaryi(nbdl),nmetaaryr4(nbdl),nmetaaryr8(nbdl),nmetaaryc(nbdl))
      allocate(extrameta(nbdl))
    endif
    nmetavari(mybdl)=0; nmetavarr4(mybdl)=0; nmetavarr8(mybdl)=0; nmetavarc(mybdl)=0
    nmetaaryi(mybdl)=0; nmetaaryr4(mybdl)=0; nmetaaryr8(mybdl)=0; nmetaaryc(mybdl)=0
    do i=1, attCount

      call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
        attnestflag=ESMF_ATTNEST_OFF, attributeIndex=i,name=attName, typekind=typekind, &
        itemCount=n, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

! add this attribute to the list of transfers
      if (typekind==ESMF_TYPEKIND_CHARACTER) then
        if( n == 1) then
          nmetavarc(mybdl) = nmetavarc(mybdl) + 1
        else if (n > 1) then
          nmetaaryc(mybdl) = nmetaaryc(mybdl) + 1
        endif
      else if (typekind==ESMF_TYPEKIND_I4) then
        if( n == 1) then
          nmetavari(mybdl) = nmetavari(mybdl) + 1
        else if (n > 1) then
          nmetaaryi(mybdl) = nmetaaryi(mybdl) + 1
        endif
      else if (typekind==ESMF_TYPEKIND_R4) then
        if( n == 1) then
          nmetavarr4(mybdl) = nmetavarr4(mybdl) + 1
        else if (n > 1) then
          nmetaaryr4(mybdl) = nmetaaryr4(mybdl) + 1
        endif
      else if (typekind==ESMF_TYPEKIND_R8) then
        if( n == 1) then
          nmetavarr8(mybdl) = nmetavarr8(mybdl) + 1
        else if (n > 1) then
          nmetaaryr8(mybdl) = nmetaaryr8(mybdl) + 1
        endif
      endif
    enddo
!    print *,'in get _global_attr, nmetavarc=',nmetavarc(mybdl),'nmetaaryc=',nmetaaryc(mybdl), &
!      'nmetavari=',nmetavari(mybdl),'nmetaaryi=',nmetaaryi(mybdl),'nmetavarr4=',nmetavarr4(mybdl), &
!      'nmetavarr8=',nmetavarr8(mybdl) 
!
! get value:
    if (nmetavari(mybdl) > 0) then 
      if(.not.allocated(variname))   allocate(variname(100,nbdl),varival(100,nbdl))
    endif
    if (nmetavarr4(mybdl) > 0) then
      if(.not.allocated(varr4name))  allocate(varr4name(100,nbdl),varr4val(100,nbdl))
    endif
    if (nmetavarr8(mybdl) > 0) then
      if(.not.allocated(varr8name))  allocate(varr8name(100,nbdl),varr8val(100,nbdl))
    endif
    if (nmetavarc(mybdl) > 0)  then
      if(.not.allocated(varcname))  allocate(varcname(100,nbdl),varcval(100,nbdl))
    endif
!
    ni=0; nr4=0; nr8=0; nc=0
    do i=1, attCount

      call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
        attnestflag=ESMF_ATTNEST_OFF, attributeIndex=i, name=attName, &
        typekind=typekind, rc=rc)
      
      if (typekind==ESMF_TYPEKIND_I4 ) then
         ni = ni + 1
         variname(ni,mybdl) = trim(attName)
         call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
           name=trim(variname(ni,mybdl)), value=varival(ni,mybdl), rc=rc)
         if (trim(variname(ni,mybdl)) == 'ncnsto') ntrac=varival(ni,mybdl)
      else if (typekind==ESMF_TYPEKIND_R4) then
         nr4 = nr4 + 1
         varr4name(nr4,mybdl) = trim(attName)
         call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
           name=trim(variname(nr4,mybdl)), value=varival(nr4,mybdl), rc=rc)
      else if (typekind==ESMF_TYPEKIND_R8) then
         nr8 = nr8 + 1
         varr8name(nr8,mybdl) = trim(attName)
         call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
           name=trim(variname(nr8,mybdl)), value=varival(nr8,mybdl), rc=rc)
      else if (typekind==ESMF_TYPEKIND_CHARACTER) then
         nc = nc + 1
         varcname(nc,mybdl) = trim(attName)
         call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
           name=trim(varcname(nc,mybdl)), value=varcval(nc,mybdl), rc=rc)
      endif

!      if(nmetavari(mybdl)>0) print *,'variname=',variname(1,mybdl),'varival=',varival(1,mybdl)
!      if(nmetavarc(mybdl)>0) print *,'varcname=',varcname(1,mybdl),'varcval=',varcval(1,mybdl)
!
      if( nmetavari(mybdl)>0 .or. nmetavarc(mybdl)>0 .or. nmetavarr4(mybdl) >0 .or. nmetavarr8(mybdl)>0) then
        extrameta(mybdl) = .true.
      else
        extrameta(mybdl) = .false.
      endif

    enddo

  end subroutine get_global_attr
!
!----------------------------------------------------------------------------------------

end module module_write_nemsio
