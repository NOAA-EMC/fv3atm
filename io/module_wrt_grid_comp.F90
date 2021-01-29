!-----------------------------------------------------------------------
!
    module module_wrt_grid_comp
!
!-----------------------------------------------------------------------
!***  This module includes the funcationailty of write gridded component.
!-----------------------------------------------------------------------
!***  At intialization step, write grid is defined. The firecast field
!***  bundle is mirrored and output field information inside the field 
!***  bundle is used to create ESMF field on the write grid and added in
!***  the mirrror field bundle on write grid component. Also the IO_BaseTime
!***  is set to the intial clock time.
!***  At the run step, output time is set from the write grid comp clock 
!***  the ESMF field bundles that contains the data on write grid are 
!***  writen out through ESMF field bundle write to netcdf files.
!***  The ESMF field bundle write uses parallel write, so if output grid
!***  is cubed sphere grid, the  six tiles file will be written out at
!***  same time.
!-----------------------------------------------------------------------
!***
!***  Revision history
!***
!     Jul 2017:  J. Wang/G. Theurich  - initial code for fv3 write grid component
!     Aug 2017:  J. Wang              - add nemsio binary output for Gaussian grid
!     Mar 2018:  S  Moorthi           - changing cfhour to accomodate up to 99999 hours
!     Aug 2019:  J. Wang              - add inline post
!
!---------------------------------------------------------------------------------
!
      use fms_io_mod,         only: field_exist, read_data

      use esmf
      use write_internal_state
      use module_fv3_io_def,   only : num_pes_fcst,lead_wrttask, last_wrttask,  &
                                      n_group, num_files, app_domain,           &
                                      filename_base, output_grid, output_file,  &
                                      imo,jmo,ichunk2d,jchunk2d,write_nemsioflip,&
                                      ichunk3d,jchunk3d,kchunk3d,nbits,         &
                                      nsout => nsout_io,                        &
                                      cen_lon, cen_lat,                         &
                                      lon1, lat1, lon2, lat2, dlon, dlat,       &
                                      stdlat1, stdlat2, dx, dy, iau_offset
      use module_write_nemsio, only : nemsio_first_call, write_nemsio
      use module_write_netcdf, only : write_netcdf
      use physcons,            only : pi => con_pi
      use inline_post,         only : inline_post_run, inline_post_getattr
      use module_write_netcdf_parallel, only : write_netcdf_parallel
!
!-----------------------------------------------------------------------
!
      implicit none
!
      include 'mpif.h'
!
!-----------------------------------------------------------------------

      private
!
!-----------------------------------------------------------------------
!

      real, parameter   :: rdgas=287.04, grav=9.80
      real, parameter   :: stndrd_atmos_ps = 101325.
      real, parameter   :: stndrd_atmos_lapse = 0.0065
!
      integer,save      :: lead_write_task                                !<-- Rank of the first write task in the write group
      integer,save      :: last_write_task                                !<-- Rank of the last write task in the write group
      integer,save      :: ntasks                                         !<-- # of write tasks in the current group 

      integer,save      :: mytile                                         !<-- the tile number in write task
      integer,save      :: wrt_mpi_comm                                   !<-- the mpi communicator in the write comp
      integer,save      :: idate(7)
      logical,save      :: first_init=.false.
      logical,save      :: first_run=.false.
      logical,save      :: first_getlatlon=.true.
      logical,save      :: first_getmaskwrt=.true.                        !<-- for mask the output grid of the write comp
      logical,save      :: change_wrtidate=.false.
!
!-----------------------------------------------------------------------
!
      type(wrt_internal_state),pointer :: wrt_int_state                 ! The internal state pointer.
      type(ESMF_FieldBundle)           :: gridFB
      integer                          :: FBcount
      character(len=esmf_maxstr),allocatable    :: fcstItemNameList(:)
      real(ESMF_KIND_R4), dimension(:,:), allocatable  :: maskwrt
!
!-----------------------------------------------------------------------
      REAL(KIND=8)             :: btim,btim0
      REAL(KIND=8),PUBLIC,SAVE :: write_init_tim, write_run_tim  
      REAL(KIND=8), parameter  :: radi=180.0d0/pi
!-----------------------------------------------------------------------
!
      public SetServices
!
      interface splat
        module procedure splat4
        module procedure splat8
      end interface splat
!
      type optimizeT
        type(ESMF_State)                   :: state
        type(ESMF_GridComp), allocatable   :: comps(:)
      end type

      contains
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine SetServices(wrt_comp, rc)
        type(ESMF_GridComp)  :: wrt_comp
        integer, intent(out) :: rc

        rc = ESMF_SUCCESS

        call ESMF_GridCompSetEntryPoint(wrt_comp, ESMF_METHOD_INITIALIZE, &
                                        userRoutine=wrt_initialize, rc=rc)
        if(rc/=0) write(*,*)'Error: write grid comp, initial'
!
        call ESMF_GridCompSetEntryPoint(wrt_comp, ESMF_METHOD_RUN, &
                                        userRoutine=wrt_run, rc=rc)
        if(rc/=0) write(*,*)'Error: write grid comp, run'
!
        call ESMF_GridCompSetEntryPoint(wrt_comp, ESMF_METHOD_FINALIZE, &
                                        userRoutine=wrt_finalize, rc=rc)
        if(rc/=0) write(*,*)'Error: write grid comp, run'

      end subroutine SetServices
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine wrt_initialize(wrt_comp, imp_state_write, exp_state_write, clock, rc)
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE WRITE GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
      type(esmf_GridComp)               :: wrt_comp
      type(ESMF_State)                  :: imp_state_write, exp_state_write
      type(esmf_Clock)                  :: clock
      integer,intent(out)               :: rc
!
!***  LOCAL VARIABLES
!
      TYPE(ESMF_VM)                           :: VM
      type(write_wrap)                        :: WRAP
      type(wrt_internal_state),pointer        :: wrt_int_state

      integer                                 :: ISTAT, tl, i, j, n, k
      integer,dimension(2,6)                  :: decomptile
      integer,dimension(2)                    :: regDecomp !define delayout for the nest grid
      integer                                 :: fieldCount
      integer                                 :: vm_mpi_comm
      character(40)                           :: fieldName, axesname,longname
      type(ESMF_Config)                       :: cf
      type(ESMF_DELayout)                     :: delayout
      type(ESMF_Grid)                         :: wrtGrid, fcstGrid
      type(ESMF_Array)                        :: array_work, array
      type(ESMF_FieldBundle)                  :: fieldbdl_work
      type(ESMF_Field)                        :: field_work, field
      type(ESMF_Decomp_Flag)                  :: decompflagPTile(2,6)

      character(len=80)                       :: attrValueSList(2)
      type(ESMF_StateItem_Flag), allocatable  :: fcstItemTypeList(:)
      type(ESMF_FieldBundle)                  :: fcstFB, wrtFB, fieldbundle
      type(ESMF_Field),          allocatable  :: fcstField(:)
      type(ESMF_TypeKind_Flag)                :: typekind
      character(len=80),         allocatable  :: fieldnamelist(:)
      integer                                 :: fieldDimCount, gridDimCount
      integer,                   allocatable  :: petMap(:)
      integer,                   allocatable  :: gridToFieldMap(:)
      integer,                   allocatable  :: ungriddedLBound(:)
      integer,                   allocatable  :: ungriddedUBound(:)
      character(len=80)                       :: attName
      character(len=80),         allocatable  :: attNameList(:),attNameList2(:)
      type(ESMF_TypeKind_Flag),  allocatable  :: typekindList(:)
      character(len=80)                       :: valueS
      integer                                 :: valueI4
      real(ESMF_KIND_R4)                      :: valueR4
      real(ESMF_KIND_R8)                      :: valueR8

      integer :: attCount, axeslen, jidx, idx, noutfile
      character(19)  :: newdate
      character(128) :: FBlist_outfilename(100), outfile_name
      character(128),dimension(:,:), allocatable    :: outfilename
      real(8), dimension(:),         allocatable    :: slat
      real(8), dimension(:),         allocatable    :: lat, lon
      real(ESMF_KIND_R8), dimension(:,:), pointer   :: lonPtr, latPtr
      real(ESMF_KIND_R8)                            :: rot_lon, rot_lat
      real(ESMF_KIND_R8)                            :: geo_lon, geo_lat
      real(ESMF_KIND_R8)                            :: lon1_r8, lat1_r8
      real(ESMF_KIND_R8)                            :: x1, y1, x, y, delat
      type(ESMF_TimeInterval)                       :: IAU_offsetTI
      type(ESMF_DataCopy_Flag) :: copyflag=ESMF_DATACOPY_REFERENCE
!     real(8),parameter :: PI=3.14159265358979d0

      character(256)                          :: gridfile

!
      logical,save                            :: first=.true.
      logical                                 :: lprnt
!test
      integer myattCount
      real(ESMF_KIND_R8),dimension(:,:), pointer :: glatPtr, glonPtr
!
!----------------------------------------------------------------------- 
!*********************************************************************** 
!----------------------------------------------------------------------- 
!
      rc = ESMF_SUCCESS
!
!----------------------------------------------------------------------- 
!***  initialize the write component timers.
!----------------------------------------------------------------------- 
!
      write_init_tim = 0.
      write_run_tim  = 0.
      btim0 = MPI_Wtime()
!
!----------------------------------------------------------------------- 
!***  set the write component's internal state.
!----------------------------------------------------------------------- 
!
      allocate(wrt_int_state,stat=RC)
      wrap%write_int_state => wrt_int_state
      call ESMF_GridCompSetInternalState(wrt_comp, wrap, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      allocate(wrt_int_state%wrtFB(num_files))
!
      call ESMF_VMGetCurrent(vm=VM,rc=RC)        
      call ESMF_VMGet(vm=VM, localPet=wrt_int_state%mype,               &
                      petCount=wrt_int_state%petcount,mpiCommunicator=vm_mpi_comm,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call mpi_comm_dup(vm_mpi_comm, wrt_mpi_comm, rc)

      ntasks = wrt_int_state%petcount
      jidx   = wrt_int_state%petcount/6
      lead_write_task = 0
      last_write_task = ntasks -1
      lprnt = lead_write_task == wrt_int_state%mype

!      print *,'in wrt, lead_write_task=', &
!         lead_write_task,'last_write_task=',last_write_task, &
!         'mype=',wrt_int_state%mype,'jidx=',jidx,' comm=',wrt_mpi_comm
!
!-----------------------------------------------------------------------
!*** get configuration variables
!-----------------------------------------------------------------------
!
      call esmf_GridCompGet(gridcomp=wrt_comp,config=CF,rc=RC)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out
! variables for post
      call ESMF_ConfigGetAttribute(config=CF,value=wrt_int_state%output_history,default=.true., &
                                   label='output_history:',rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_ConfigGetAttribute(config=CF,value=wrt_int_state%write_dopost,default=.false., &
                                   label='write_dopost:',rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      if( wrt_int_state%write_dopost ) then
#ifdef NO_INLINE_POST
        rc = ESMF_RC_NOT_IMPL
        print *,'inline post not available on this machine'
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
#endif
        call esmf_configgetattribute(cf,wrt_int_state%post_nlunit,default=777,label='nlunit:',rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ConfigGetAttribute(config=CF,value=wrt_int_state%post_namelist,default='itag', &
                                     label ='post_namelist:',rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      endif
!
!-----------------------------------------------------------------------
!*** Create the cubed sphere grid with field on PETs
!*** first try: Create cubed sphere grid from file
!-----------------------------------------------------------------------
!
      if ( trim(output_grid) == 'cubed_sphere_grid' ) then

        mytile = mod(wrt_int_state%mype,ntasks)+1
        if ( trim(app_domain) == 'global' ) then
          do tl=1,6
            decomptile(1,tl) = 1
            decomptile(2,tl) = jidx
            decompflagPTile(:,tl) = (/ESMF_DECOMP_SYMMEDGEMAX,ESMF_DECOMP_SYMMEDGEMAX/)
          enddo
          call ESMF_AttributeGet(imp_state_write, convention="NetCDF", purpose="FV3", &
                                 name="gridfile", value=gridfile, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          CALL ESMF_LogWrite("wrtComp: gridfile:"//trim(gridfile),ESMF_LOGMSG_INFO,rc=rc)
          wrtgrid = ESMF_GridCreateMosaic(filename="INPUT/"//trim(gridfile),                                 &
                                          regDecompPTile=decomptile,tileFilePath="INPUT/",                   &
                                          decompflagPTile=decompflagPTile,                                   &
                                          staggerlocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
                                          name='wrt_grid', rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        else
          if(trim(app_domain) == 'nested') then
            gridfile='grid.nest02.tile7.nc'
          else if(trim(app_domain) == 'regional') then
            gridfile='grid.tile7.halo0.nc'
          endif
          regDecomp(1) = 1
          regDecomp(2) = ntasks
          allocate(petMap(ntasks))
          do i=1, ntasks
            petMap(i) = i-1
          enddo
          delayout = ESMF_DELayoutCreate(petMap=petMap, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          ! create the nest Grid by reading it from file but use DELayout
          wrtGrid = ESMF_GridCreate(filename="INPUT/"//trim(gridfile),                                       &
                                    fileformat=ESMF_FILEFORMAT_GRIDSPEC, regDecomp=regDecomp,                &
                                    decompflag=(/ESMF_DECOMP_SYMMEDGEMAX,ESMF_DECOMP_SYMMEDGEMAX/),          &
                                    delayout=delayout, isSphere=.false., indexflag=ESMF_INDEX_DELOCAL,       &
                                    rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          print *,'in nested/regional cubed_sphere grid, regDecomp=',regDecomp,' PetMap=',petMap(1),petMap(ntasks), &
            'gridfile=',trim(gridfile)
          deallocate(petMap)
        endif
      else if ( trim(output_grid) == 'gaussian_grid') then
        
        wrtgrid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/),                             &
                                          maxIndex=(/imo,jmo/), regDecomp=(/1,ntasks/), &
                                          indexflag=ESMF_INDEX_GLOBAL,                  &
                                          name='wrt_grid',rc=rc)
!                                         indexflag=ESMF_INDEX_GLOBAL, coordSys=ESMF_COORDSYS_SPH_DEG

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_GridAddCoord(wrtgrid, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_GridGetCoord(wrtgrid, coordDim=1, farrayPtr=lonPtr, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_GridGetCoord(wrtgrid, coordDim=2, farrayPtr=latPtr, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
        allocate(slat(jmo), lat(jmo), lon(imo))
        call splat(4, jmo, slat)
        if(write_nemsioflip) then
          do j=1,jmo
            lat(j) = asin(slat(j)) * radi
          enddo
        else
          do j=1,jmo
            lat(jmo-j+1) = asin(slat(j)) * radi
          enddo
        endif
        wrt_int_state%latstart = lat(1)
        wrt_int_state%latlast  = lat(jmo)
        do j=1,imo
          lon(j) = 360.d0/real(imo,8) *real(j-1,8)
        enddo
        wrt_int_state%lonstart = lon(1)
        wrt_int_state%lonlast  = lon(imo)
        do j=lbound(latPtr,2),ubound(latPtr,2)
          do i=lbound(lonPtr,1),ubound(lonPtr,1)
            lonPtr(i,j) = 360.d0/real(imo,8) * real(i-1,8)
            latPtr(i,j) = lat(j)
          enddo
        enddo 
!        print *,'aft wrtgrd, Gaussian, dimi,i=',lbound(lonPtr,1),ubound(lonPtr,1), &
!         ' j=',lbound(lonPtr,2),ubound(lonPtr,2),'imo=',imo,'jmo=',jmo
!       if(wrt_int_state%mype==0) print *,'aft wrtgrd, lon=',lonPtr(1:5,1), &
!        'lat=',latPtr(1,1:5),'imo,jmo=',imo,jmo
!        lonPtr(lbound(lonPtr,1),ubound(lonPtr,2)),'lat=',latPtr(lbound(lonPtr,1),lbound(lonPtr,2)), &
!        latPtr(lbound(lonPtr,1),ubound(lonPtr,2))
        wrt_int_state%lat_start = lbound(latPtr,2)
        wrt_int_state%lat_end   = ubound(latPtr,2)
        wrt_int_state%lon_start = lbound(lonPtr,1)
        wrt_int_state%lon_end   = ubound(lonPtr,1)
        allocate( wrt_int_state%lat_start_wrtgrp(wrt_int_state%petcount))
        allocate( wrt_int_state%lat_end_wrtgrp  (wrt_int_state%petcount))
        call mpi_allgather(wrt_int_state%lat_start,1,MPI_INTEGER,    &
                           wrt_int_state%lat_start_wrtgrp, 1, MPI_INTEGER, wrt_mpi_comm, rc)
        call mpi_allgather(wrt_int_state%lat_end,  1,MPI_INTEGER,    &
                           wrt_int_state%lat_end_wrtgrp,   1, MPI_INTEGER, wrt_mpi_comm, rc)
        if( lprnt ) print *,'aft wrtgrd, Gaussian, dimj_start=',wrt_int_state%lat_start_wrtgrp, &
          'dimj_end=',wrt_int_state%lat_end_wrtgrp, 'wrt_group=',n_group
        allocate( wrt_int_state%latPtr(wrt_int_state%lon_start:wrt_int_state%lon_end, &
                  wrt_int_state%lat_start:wrt_int_state%lat_end))
        allocate( wrt_int_state%lonPtr(wrt_int_state%lon_start:wrt_int_state%lon_end, &
                  wrt_int_state%lat_start:wrt_int_state%lat_end))
        do j=wrt_int_state%lat_start,wrt_int_state%lat_end
          do i=wrt_int_state%lon_start,wrt_int_state%lon_end
            wrt_int_state%latPtr(i,j) = latPtr(i,j)
            wrt_int_state%lonPtr(i,j) = lonPtr(i,j)
          enddo
        enddo
        wrt_int_state%im = imo
        wrt_int_state%jm = jmo
        wrt_int_state%post_maptype = 4

        deallocate(slat)
      else if ( trim(output_grid) == 'global_latlon') then
        wrtgrid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/),                             &
                                          maxIndex=(/imo,jmo/), regDecomp=(/1,ntasks/), &
                                          indexflag=ESMF_INDEX_GLOBAL, name='wrt_grid',rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_GridAddCoord(wrtgrid, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_GridGetCoord(wrtgrid, coordDim=1, farrayPtr=lonPtr, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_GridGetCoord(wrtgrid, coordDim=2, farrayPtr=latPtr, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
        allocate(lat(jmo), lon(imo))
        if (mod(jmo,2) == 0) then
          ! if jmo even, lats do not include poles and equator
          delat = 180.d0/real(jmo,8)
          if(write_nemsioflip) then
            do j=1,jmo
              lat(j) = 90.d0 - 0.5*delat - real(j-1,8)*delat
            enddo
          else
            do j=1,jmo
              lat(j) = -90.d0 + 0.5*delat + real(j-1,8)*delat
            enddo
          endif
        else
          ! if jmo odd, lats include poles and equator
          delat = 180.d0/real(jmo-1,8)
          if(write_nemsioflip) then
            do j=1,jmo
              lat(j) = 90.d0 - real(j-1,8)*delat
            enddo
          else
            do j=1,jmo
              lat(j) = -90.d0 + real(j-1,8)*delat
            enddo
          endif
        endif
        wrt_int_state%latstart = lat(1)
        wrt_int_state%latlast  = lat(jmo)
        do j=1,imo
          lon(j) = 360.d0/real(imo,8) *real(j-1,8)
        enddo
        wrt_int_state%lonstart = lon(1)
        wrt_int_state%lonlast  = lon(imo)
        do j=lbound(latPtr,2),ubound(latPtr,2)
          do i=lbound(lonPtr,1),ubound(lonPtr,1)
            lonPtr(i,j) = 360.d0/real(imo,8) * real(i-1,8)
            latPtr(i,j) = lat(j)
          enddo
        enddo
        wrt_int_state%im = imo
        wrt_int_state%jm = jmo
        wrt_int_state%post_maptype = 0

      else if ( trim(output_grid) == 'regional_latlon' .or. &
                trim(output_grid) == 'rotated_latlon'  .or. &
                trim(output_grid) == 'lambert_conformal' ) then

        wrtgrid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/),                             &
                                          maxIndex=(/imo,jmo/), regDecomp=(/1,ntasks/), &
                                          indexflag=ESMF_INDEX_GLOBAL,                  &
                                          name='wrt_grid',rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_GridAddCoord(wrtgrid, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_GridGetCoord(wrtgrid, coordDim=1, farrayPtr=lonPtr, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_GridGetCoord(wrtgrid, coordDim=2, farrayPtr=latPtr, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        wrt_int_state%im = imo
        wrt_int_state%jm = jmo
        if ( trim(output_grid) == 'regional_latlon' ) then
            do j=lbound(lonPtr,2),ubound(lonPtr,2)
            do i=lbound(lonPtr,1),ubound(lonPtr,1)
              lonPtr(i,j) = lon1 + (lon2-lon1)/(imo-1) * (i-1)
              latPtr(i,j) = lat1 + (lat2-lat1)/(jmo-1) * (j-1)
            enddo
            enddo
            wrt_int_state%post_maptype = 0
        else if ( trim(output_grid) == 'rotated_latlon' ) then
            do j=lbound(lonPtr,2),ubound(lonPtr,2)
            do i=lbound(lonPtr,1),ubound(lonPtr,1)
              rot_lon = lon1 + (lon2-lon1)/(imo-1) * (i-1)
              rot_lat = lat1 + (lat2-lat1)/(jmo-1) * (j-1)
              call rtll(rot_lon, rot_lat, geo_lon, geo_lat, dble(cen_lon), dble(cen_lat))
              if (geo_lon < 0.0) geo_lon = geo_lon + 360.0
              lonPtr(i,j) = geo_lon
              latPtr(i,j) = geo_lat
            enddo
            enddo
            wrt_int_state%post_maptype = 207
        else if ( trim(output_grid) == 'lambert_conformal' ) then
            lon1_r8 = dble(lon1)
            lat1_r8 = dble(lat1)
            call lambert(dble(stdlat1),dble(stdlat2),dble(cen_lat),dble(cen_lon), &
                         lon1_r8,lat1_r8,x1,y1, 1)
            do j=lbound(lonPtr,2),ubound(lonPtr,2)
            do i=lbound(lonPtr,1),ubound(lonPtr,1)
              x = x1 + dx * (i-1)
              y = y1 + dy * (j-1)
              call lambert(dble(stdlat1),dble(stdlat2),dble(cen_lat),dble(cen_lon), &
                           geo_lon,geo_lat,x,y,-1)
              if (geo_lon <0.0) geo_lon = geo_lon + 360.0
              lonPtr(i,j) = geo_lon
              latPtr(i,j) = geo_lat
            enddo
            enddo
            wrt_int_state%post_maptype = 1
        endif

        wrt_int_state%lat_start = lbound(latPtr,2)
        wrt_int_state%lat_end   = ubound(latPtr,2)
        wrt_int_state%lon_start = lbound(lonPtr,1)
        wrt_int_state%lon_end   = ubound(lonPtr,1)
        allocate( wrt_int_state%lat_start_wrtgrp(wrt_int_state%petcount))
        allocate( wrt_int_state%lat_end_wrtgrp  (wrt_int_state%petcount))
        call mpi_allgather(wrt_int_state%lat_start,1,MPI_INTEGER,    &
                       wrt_int_state%lat_start_wrtgrp, 1, MPI_INTEGER, wrt_mpi_comm, rc)
        call mpi_allgather(wrt_int_state%lat_end,  1,MPI_INTEGER,    &
                           wrt_int_state%lat_end_wrtgrp,   1, MPI_INTEGER, wrt_mpi_comm, rc)
        allocate( wrt_int_state%latPtr(wrt_int_state%lon_start:wrt_int_state%lon_end, &
                  wrt_int_state%lat_start:wrt_int_state%lat_end))
        allocate( wrt_int_state%lonPtr(wrt_int_state%lon_start:wrt_int_state%lon_end, &
                  wrt_int_state%lat_start:wrt_int_state%lat_end))
        do j=wrt_int_state%lat_start,wrt_int_state%lat_end
        do i=wrt_int_state%lon_start,wrt_int_state%lon_end
          wrt_int_state%latPtr(i,j) = latPtr(i,j)
          wrt_int_state%lonPtr(i,j) = lonPtr(i,j)
        enddo
        enddo

      else

        write(0,*)"wrt_initialize: Unknown output_grid ", trim(output_grid)
        call ESMF_LogWrite("wrt_initialize: Unknown output_grid "//trim(output_grid),ESMF_LOGMSG_ERROR,rc=RC)
        call ESMF_Finalize(endflag=ESMF_END_ABORT)

      endif
!
!-----------------------------------------------------------------------
!***  get write grid component initial time from clock
!-----------------------------------------------------------------------
!
      call ESMF_ClockGet(clock    =CLOCK                                &  !<-- The ESMF Clock
                        ,startTime=wrt_int_state%IO_BASETIME            &  !<-- The Clock's starting time
                        ,rc       =RC)

      call ESMF_TimeGet(time=wrt_int_state%IO_BASETIME,yy=idate(1),mm=idate(2),dd=idate(3),h=idate(4), &
                        m=idate(5),s=idate(6),rc=rc)
!     if (lprnt) write(0,*) 'in wrt initial, io_baseline time=',idate,'rc=',rc
      idate(7) = 1
      wrt_int_state%idate = idate
      wrt_int_state%fdate = idate
! update IO-BASETIME and idate on write grid comp when IAU is enabled
      if(iau_offset > 0 ) then
        call ESMF_TimeIntervalSet(IAU_offsetTI, h=iau_offset, rc=rc)
        wrt_int_state%IO_BASETIME = wrt_int_state%IO_BASETIME + IAU_offsetTI
        call ESMF_TimeGet(time=wrt_int_state%IO_BASETIME,yy=idate(1),mm=idate(2),dd=idate(3),h=idate(4), &
                          m=idate(5),s=idate(6),rc=rc)
        wrt_int_state%idate = idate
        change_wrtidate = .true.
       if (lprnt) print *,'in wrt initial, with iau, io_baseline time=',idate,'rc=',rc
      endif
!
! Create field bundle
!-------------------------------------------------------------------
!
!---  check grid dim count first
      call ESMF_GridGet(wrtgrid, dimCount=gridDimCount, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
!--- Look at the incoming FieldBundles in the imp_state_write, and mirror them
!
      call ESMF_StateGet(imp_state_write, itemCount=FBCount, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      wrt_int_state%FBCount = FBCount

!      if (lprnt) write(0,*) 'in wrt,fcst FBCount=',FBCount
!      if (lprnt) write(0,*) 'in wrt,fcst wrt_int_state%FBCount=',wrt_int_state%FBCount

      allocate(fcstItemNameList(FBCount), fcstItemTypeList(FBCount))
      allocate(wrt_int_state%wrtFB_names(FBCount))
      allocate(wrt_int_state%ncount_fields(FBCount),wrt_int_state%ncount_attribs(FBCount))
      allocate(wrt_int_state%field_names(2000,FBCount))
      allocate(outfilename(2000,FBcount))
      outfilename = ''

      call ESMF_StateGet(imp_state_write, itemNameList=fcstItemNameList, &
                         itemTypeList=fcstItemTypeList, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!loop over all items in the imp_state_write and collect all FieldBundles
      do i=1, FBcount

        if (fcstItemTypeList(i) == ESMF_STATEITEM_FIELDBUNDLE) then

          call ESMF_StateGet(imp_state_write, itemName=fcstItemNameList(i), &
                             fieldbundle=fcstFB, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! create a mirror FieldBundle and add it to importState
          fieldbundle = ESMF_FieldBundleCreate(name="mirror_"//trim(fcstItemNameList(i)), rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_StateAdd(imp_state_write, (/fieldbundle/), rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! copy the fcstFB Attributes to the mirror FieldBundle
          call ESMF_AttributeCopy(fcstFB, fieldbundle, attcopy=ESMF_ATTCOPY_REFERENCE, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
          call ESMF_FieldBundleGet(fcstFB, fieldCount=fieldCount, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if (fieldCount > 0) then

            wrt_int_state%ncount_fields(i) = fieldCount

            allocate(fcstField(fieldCount))
            call ESMF_FieldBundleGet(fcstFB, fieldList=fcstField,     &
                                     itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            do j=1, fieldCount

              call ESMF_FieldGet(fcstField(j), typekind=typekind, &
                                 dimCount=fieldDimCount, name=fieldName, grid=fcstGrid, rc=rc)

              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

              allocate(gridToFieldMap(gridDimCount))
              allocate(ungriddedLBound(fieldDimCount-gridDimCount))
              allocate(ungriddedUBound(fieldDimCount-gridDimCount))

              call ESMF_FieldGet(fcstField(j), gridToFieldMap=gridToFieldMap,                      &
                                 ungriddedLBound=ungriddedLBound, ungriddedUBound=ungriddedUBound, &
                                 rc=rc)
              CALL ESMF_LogWrite("after field create on wrt comp",ESMF_LOGMSG_INFO,rc=RC)

!             if (lprnt) print *,'in wrt,fcstfld,fieldname=',                                         &
!                        trim(fieldname),'fieldDimCount=',fieldDimCount,'gridDimCount=',gridDimCount, &
!                        'gridToFieldMap=',gridToFieldMap,'ungriddedLBound=',ungriddedLBound,         &
!                        'ungriddedUBound=',ungriddedUBound,'rc=',rc

              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! create the mirror field

              CALL ESMF_LogWrite("call field create on wrt comp",ESMF_LOGMSG_INFO,rc=RC)
              field_work = ESMF_FieldCreate(wrtGrid, typekind, name=fieldName, &
                                            gridToFieldMap=gridToFieldMap,     &
                                            ungriddedLBound=ungriddedLBound,   &
                                            ungriddedUBound=ungriddedUBound, rc=rc)
              CALL ESMF_LogWrite("aft call field create on wrt comp",ESMF_LOGMSG_INFO,rc=RC)

              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

              wrt_int_state%field_names(j,i) = trim(fieldName)

              call ESMF_AttributeCopy(fcstField(j), field_work, attcopy=ESMF_ATTCOPY_REFERENCE, rc=rc)

              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
! get output file name
              call ESMF_AttributeGet(fcstField(j), convention="NetCDF", purpose="FV3", &
                                     name="output_file", value=outfile_name, rc=rc)

              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              CALL ESMF_LogWrite("bf fcstfield, get output_file "//trim(outfile_name)//" "//trim(fieldName),ESMF_LOGMSG_INFO,rc=RC)
              if (trim(outfile_name) /= '') then
                outfilename(j,i) = trim(outfile_name)
              endif
              CALL ESMF_LogWrite("af fcstfield, get output_file",ESMF_LOGMSG_INFO,rc=RC)

!             if (lprnt) print *,' i=',i,' j=',j,' outfilename=',trim(outfilename(j,i))

! add the mirror field to the mirror FieldBundle
              call ESMF_FieldBundleAdd(fieldbundle, (/field_work/), rc=rc)

              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! local garbage collection
              deallocate(gridToFieldMap, ungriddedLBound, ungriddedUBound)
            enddo
!
            call ESMF_AttributeCopy(fcstGrid, wrtGrid, &
                                    attcopy=ESMF_ATTCOPY_REFERENCE, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            deallocate(fcstField)

          endif !if (fieldCount > 0) then

        else  ! anything but a FieldBundle in the state is unexpected here
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
                                msg="Only FieldBundles supported in fcstState.", line=__LINE__, file=__FILE__)
          return
        endif


!end FBCount
      enddo
!
!loop over all items in the imp_state_write and count output FieldBundles

      call get_outfile(FBcount, outfilename,FBlist_outfilename,noutfile)
      wrt_int_state%FBCount = noutfile

!
!create output field bundles
      do i=1, wrt_int_state%FBcount

        wrt_int_state%wrtFB_names(i) = trim(FBlist_outfilename(i))
        wrt_int_state%wrtFB(i) = ESMF_FieldBundleCreate(name=trim(wrt_int_state%wrtFB_names(i)), rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        do n=1, FBcount

          call ESMF_StateGet(imp_state_write, itemName="mirror_"//trim(fcstItemNameList(n)), &
                             fieldbundle=fcstFB, rc=rc)

          if( index(trim(fcstItemNameList(n)),trim(FBlist_outfilename(i))) > 0 ) then
!
! copy the mirror fcstfield bundle Attributes to the output field bundle
            call ESMF_AttributeCopy(fcstFB,  wrt_int_state%wrtFB(i), &
                                    attcopy=ESMF_ATTCOPY_REFERENCE, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
            call ESMF_FieldBundleGet(fcstFB, fieldCount=fieldCount, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            allocate(fcstField(fieldCount),fieldnamelist(fieldCount))
            call ESMF_FieldBundleGet(fcstFB, fieldList=fcstField, fieldNameList=fieldnamelist,   &
                                     itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            do j=1, fieldCount

              call ESMF_AttributeGet(fcstField(j),convention="NetCDF", purpose="FV3", &
                                     name='output_file',value=outfile_name, rc=rc)

              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!             if (lprnt) print *,'in wrt,add field,i=',i,'n=',n,' j=',j, &
!                        'fieldname=',trim(fieldnamelist(j)), ' outfile_name=',trim(outfile_name), &
!                       ' field bundle name, FBlist_outfilename(i)=',trim(FBlist_outfilename(i))

              if( trim(outfile_name) == trim(FBlist_outfilename(i))) then
                call ESMF_FieldBundleAdd(wrt_int_state%wrtFB(i), (/fcstField(j)/), rc=rc)

                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              endif

            enddo
            deallocate(fcstField, fieldnamelist)

          endif

! add output grid related attributes

            call ESMF_AttributeAdd(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                   attrList=(/"source","grid  "/), rc=rc)
            call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                   name="source", value="FV3GFS", rc=rc)

            if (trim(output_grid) == 'cubed_sphere_grid') then

              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="grid", value="cubed_sphere", rc=rc)

            else if (trim(output_grid) == 'gaussian_grid') then

              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="grid", value="gaussian", rc=rc)
              call ESMF_AttributeAdd(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     attrList=(/"im","jm"/), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="im", value=imo, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="jm", value=jmo, rc=rc)

            else if (trim(output_grid) == 'global_latlon') then

              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="grid", value="latlon", rc=rc)
              call ESMF_AttributeAdd(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     attrList=(/"lonstart","latstart","lonlast ","latlast "/), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="lonstart", value=wrt_int_state%lonstart, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="latstart", value=wrt_int_state%latstart, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="lonlast", value=wrt_int_state%lonlast, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="latlast", value=wrt_int_state%latlast, rc=rc)

            else if (trim(output_grid) == 'regional_latlon') then

              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="grid", value="latlon", rc=rc)
              call ESMF_AttributeAdd(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     attrList=(/"lon1","lat1","lon2","lat2","dlon","dlat"/), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="lon1", value=lon1, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="lat1", value=lat1, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="lon2", value=lon2, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="lat2", value=lat2, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="dlon", value=dlon, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="dlat", value=dlat, rc=rc)

            else if (trim(output_grid) == 'rotated_latlon') then

              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="grid", value="rotated_latlon", rc=rc)
              call ESMF_AttributeAdd(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     attrList=(/"cen_lon",&
                                                "cen_lat",&
                                                "lon1   ",&
                                                "lat1   ",&
                                                "lon2   ",&
                                                "lat2   ",&
                                                "dlon   ",&
                                                "dlat   "/), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="cen_lon", value=cen_lon, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="cen_lat", value=cen_lat, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="lon1", value=lon1, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="lat1", value=lat1, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="lon2", value=lon2, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="lat2", value=lat2, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="dlon", value=dlon, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="dlat", value=dlat, rc=rc)

            else if (trim(output_grid) == 'lambert_conformal') then

              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="grid", value="lambert_conformal", rc=rc)
              call ESMF_AttributeAdd(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     attrList=(/"cen_lon",&
                                                "cen_lat",&
                                                "stdlat1",&
                                                "stdlat2",&
                                                "nx     ",&
                                                "ny     ",&
                                                "lon1   ",&
                                                "lat1   ",&
                                                "dx     ",&
                                                "dy     "/), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="cen_lon", value=cen_lon, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="cen_lat", value=cen_lat, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="stdlat1", value=stdlat1, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="stdlat2", value=stdlat2, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="nx", value=imo, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="ny", value=jmo, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="lat1", value=lat1, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="lon1", value=lon1, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="dx", value=dx, rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="dy", value=dy, rc=rc)

            end if

        enddo ! end FBcount
      enddo ! end wrt_int_state%FBcount
!
! add time Attribute
! look at the importState attributes and copy those starting with "time"

      call ESMF_AttributeGet(imp_state_write, convention="NetCDF", purpose="FV3", &
                             attnestflag=ESMF_ATTNEST_OFF, count=attCount, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_LogWrite("Write component AttributeGet, attCount ", ESMF_LOGMSG_INFO, rc=rc)

! prepare the lists needed to transfer attributes

      allocate(attNameList(attCount), attNameList2(attCount))
      allocate(typekindList(attCount))

! loop over all the attributes on importState within AttPack
      j = 1
      k = 1
      do i=1, attCount
        call ESMF_AttributeGet(imp_state_write, convention="NetCDF", purpose="FV3",          &
                               attnestflag=ESMF_ATTNEST_OFF, attributeIndex=i, name=attName, &
                               typekind=typekind, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! test for name starting with "time"
        if (index(trim(attName), "time") == 1) then

! add this attribute to the list of transfers
          attNameList(j)  = attName
          typekindList(j) = typekind
          j = j + 1
          if (index(trim(attName), "time:") == 1) then
          ! store names of attributes starting with "time:" for later use
            attNameList2(k) = attName
            k = k+1
          endif
        endif
      enddo

! add the transfer attributes from importState to grid
    call ESMF_AttributeAdd(wrtgrid, convention="NetCDF", purpose="FV3", &
                           attrList=attNameList(1:j-1), rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! loop over the added attributes, access the value (only scalar allowed),
! and set them on the grid
    do i=1, j-1
      if (typekindList(i) == ESMF_TYPEKIND_CHARACTER) then
        call ESMF_AttributeGet(imp_state_write,                    &
                               convention="NetCDF", purpose="FV3", &
                               name=trim(attNameList(i)), value=valueS, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! update the time:units when idate on write grid component is changed 
        if ( index(trim(attNameList(i)),'time:units')>0) then
          if ( change_wrtidate ) then
            idx = index(trim(valueS),' since ')
            if(lprnt) print *,'in write grid comp, time:unit=',trim(valueS)
            write(newdate,'(I4.4,a,I2.2,a,I2.2,a,I2.2,a,I2.2,a,I2.2)') idate(1),'-',   &
              idate(2),'-',idate(3),' ',idate(4),':',idate(5),':',idate(6)
            valueS = valueS(1:idx+6)//newdate
            if(lprnt) print *,'in write grid comp, new time:unit=',trim(valueS)
          endif
        endif
        call ESMF_AttributeSet(wrtgrid, convention="NetCDF", purpose="FV3", &
                               name=trim(attNameList(i)), value=valueS, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      else if (typekindList(i) == ESMF_TYPEKIND_I4) then
        call ESMF_AttributeGet(imp_state_write,                    &
                               convention="NetCDF", purpose="FV3", &
                               name=trim(attNameList(i)), value=valueI4, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_AttributeSet(wrtgrid, convention="NetCDF", purpose="FV3", &
                               name=trim(attNameList(i)), value=valueI4, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      else if (typekindList(i) == ESMF_TYPEKIND_R4) then
        call ESMF_AttributeGet(imp_state_write,                    &
                               convention="NetCDF", purpose="FV3", &
                               name=trim(attNameList(i)), value=valueR4, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_AttributeSet(wrtgrid, convention="NetCDF", purpose="FV3", &
                               name=trim(attNameList(i)), value=valueR4, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      else if (typekindList(i) == ESMF_TYPEKIND_R8) then
        call ESMF_AttributeGet(imp_state_write,                    &
                               convention="NetCDF", purpose="FV3", &
                               name=trim(attNameList(i)), value=valueR8, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_AttributeSet(wrtgrid, convention="NetCDF", purpose="FV3", &
                               name=trim(attNameList(i)), value=valueR8, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      endif
    enddo

! Add special attribute that holds names of "time" related attributes
! for faster access during Run().
    call ESMF_AttributeAdd(wrtgrid, convention="NetCDF", purpose="FV3", &
                           attrList=(/"TimeAttributes"/), rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeSet(wrtgrid, convention="NetCDF", purpose="FV3", &
                           name="TimeAttributes", valueList=attNameList2(1:k-1), rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    deallocate(attNameList, attNameList2, typekindList)

!
!*** create temporary field bundle for  axes information
! write the Grid coordinate arrays into the output files via temporary FB

      gridFB = ESMF_FieldBundleCreate(rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeGet(wrtGrid, convention="NetCDF", purpose="FV3", &
                             name="ESMF:gridded_dim_labels", valueList=attrValueSList, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_GridGetCoord(wrtGrid, coordDim=1, &
                             staggerloc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!      write(0,*) 'create gridFB,fieldname=',trim(attrValueSList(1)),trim(attrValueSList(2)), &
!         'lon value=',array(1:5)

      field = ESMF_FieldCreate(wrtGrid, array, name=trim(attrValueSList(1)), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!add attribute info
! long name
      call ESMF_AttributeAdd(field,convention="NetCDF",purpose="FV3",  &
                             attrList=(/'long_name'/), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(field,convention="NetCDF",purpose="FV3",name='long_name', &
                             value="T-cell longitude", rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
! units
      call ESMF_AttributeAdd(field,convention="NetCDF",purpose="FV3",  &
                             attrList=(/'units'/), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(field,convention="NetCDF",purpose="FV3",name='units', &
                             value="degrees_E", rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
! cartesian_axis
      call ESMF_AttributeAdd(field,convention="NetCDF",purpose="FV3",  &
                            attrList=(/'cartesian_axis'/), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(field,convention="NetCDF",purpose="FV3",name='cartesian_axis', &
                             value="X", rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! add field to bundle
      call ESMF_FieldBundleAdd(gridFB, (/field/), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
! get 2nd dimension
      call ESMF_GridGetCoord(wrtGrid, coordDim=2, &
                             staggerloc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!      write(0,*) 'create gridFB,fieldname=',trim(attrValueSList(1)),trim(attrValueSList(2)), &
!         'lat value=',array(1:5,1),array(1,1:5)

      field = ESMF_FieldCreate(wrtGrid, array, name=trim(attrValueSList(2)), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!add attribute info
! long name
      call ESMF_AttributeAdd(field,convention="NetCDF",purpose="FV3",  &
                             attrList=(/'long_name'/), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(field,convention="NetCDF",purpose="FV3",name='long_name', &
                              value="T-cell latitude", rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
! units
      call ESMF_AttributeAdd(field,convention="NetCDF",purpose="FV3",  &
                             attrList=(/'units'/), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(field,convention="NetCDF",purpose="FV3",name='units', &
                             value="degrees_N", rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! cartesian_axis
      call ESMF_AttributeAdd(field,convention="NetCDF",purpose="FV3",  &
                             attrList=(/'cartesian_axis'/), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(field,convention="NetCDF",purpose="FV3",name='cartesian_axis', &
                             value="Y", rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_FieldBundleAdd(gridFB, (/field/), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
!-----------------------------------------------------------------------
!***  SET THE FIRST HISTORY FILE'S TIME INDEX.
!-----------------------------------------------------------------------
!
      wrt_int_state%NFHOUR = 0
!
!-----------------------------------------------------------------------
!***  Initialize for POST
!-----------------------------------------------------------------------
!
      call ESMF_LogWrite("before initialize for POST", ESMF_LOGMSG_INFO, rc=rc)
      if (lprnt) print *,'in wrt grid comp, dopost=',wrt_int_state%write_dopost
      if( wrt_int_state%write_dopost ) then
        call inline_post_getattr(wrt_int_state)
      endif
!
!-----------------------------------------------------------------------
!***  Initialize for nemsio file
!-----------------------------------------------------------------------
!
      call ESMF_LogWrite("before initialize for nemsio file", ESMF_LOGMSG_INFO, rc=rc)
      do i= 1, wrt_int_state%FBcount
         if (trim(output_grid) == 'gaussian_grid' .and. trim(output_file(i)) == 'nemsio') then
!           if (lprnt) write(0,*) 'in wrt initial, befnemsio_first_call wrt_int_state%FBcount=',wrt_int_state%FBcount
             call nemsio_first_call(wrt_int_state%wrtFB(i), imo, jmo,         &
                                   wrt_int_state%mype, ntasks, wrt_mpi_comm, &
                                   wrt_int_state%FBcount, i, idate, lat, lon, rc) 
         endif
      enddo
      call ESMF_LogWrite("after initialize for nemsio file", ESMF_LOGMSG_INFO, rc=rc)
!
!-----------------------------------------------------------------------
!
      IF(RC /= ESMF_SUCCESS) THEN
        WRITE(0,*)"FAIL: Write_Initialize."
!      ELSE
!        WRITE(0,*)"PASS: Write_Initialize."
      ENDIF
!
!      write_init_tim = MPI_Wtime() - btim0
!
!----------------------------------------------------------------------- 
!
      end subroutine wrt_initialize
!
!----------------------------------------------------------------------- 
!####################################################################### 
!----------------------------------------------------------------------- 
!
      subroutine wrt_run(wrt_comp, imp_state_write, exp_state_write,clock,rc)
!
!----------------------------------------------------------------------- 
!***  the run step for the write gridded component.  
!----------------------------------------------------------------------- 
!
!
      type(ESMF_GridComp)            :: wrt_comp
      type(ESMF_State)               :: imp_state_write, exp_state_write
      type(ESMF_Clock)               :: clock
      integer,intent(out)            :: rc 
!
!-----------------------------------------------------------------------
!***  local variables
!
      TYPE(ESMF_VM)                         :: VM
      type(ESMF_FieldBundle)                :: file_bundle 
      type(ESMF_Time)                       :: currtime
      type(ESMF_TypeKind_Flag)              :: datatype
      type(ESMF_Field)                      :: field_work
      type(ESMF_Grid)                       :: grid_work, fbgrid, wrtgrid
      type(ESMF_Array)                      :: array_work
      type(ESMF_State),save                 :: stateGridFB
      type(optimizeT), save                 :: optimize(4)
      type(ESMF_GridComp), save, allocatable   :: compsGridFB(:)
!
      type(write_wrap)                      :: wrap
      type(wrt_internal_state),pointer      :: wrt_int_state
!
      integer,dimension(:),allocatable,save :: ih_int, ih_real
!
      INTEGER,SAVE                          :: NPOSN_1,NPOSN_2
!
      integer                               :: i,j,n,mype,nolog
!
      integer                               :: nf_hours,nf_seconds, nf_minutes,     &
                                               nseconds,nseconds_num,nseconds_den
!
      integer                               :: id
      integer                               :: nbdl, idx, date(6), ndig
      integer                               :: step=1
!
      REAL                                  :: DEGRAD
!
      logical                               :: opened
      logical                               :: lmask_fields
      logical,save                          :: first=.true.
      logical,save                          :: file_first=.true.
!
      character(esmf_maxstr)            :: filename,compname,bundle_name
      character(40)                         :: cfhour, cform
      character(10)                         :: stepString
      character(80)                         :: attrValueS
      integer                               :: attrValueI
      real                                  :: attrValueR
      real(ESMF_KIND_R8)                    :: time

!
!-----------------------------------------------------------------------
!
      real(kind=8)  :: wait_time, MPI_Wtime
      real(kind=8)  :: times,times2,etim
      character(10) :: timeb
      real(kind=8)  :: tbeg,tend
      real(kind=8)  :: wbeg,wend

      integer fieldcount, dimCount
      real(kind=ESMF_KIND_R8), dimension(:,:,:), pointer   :: datar8
      real(kind=ESMF_KIND_R8), dimension(:,:),   pointer   :: datar82d
!
      integer myattCount
      logical lprnt
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      tbeg = MPI_Wtime()
      rc   = esmf_success
!
!-----------------------------------------------------------------------
!***  get the current write grid comp name, id, and internal state
!
      call ESMF_GridCompGet(wrt_comp, name=compname, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!      print *,'in wrt run. compname=',trim(compname),' rc=',rc

! instance id from name
      read(compname(10:11),"(I2)") id

! Provide log message indicating which wrtComp is active
      call ESMF_LogWrite("Write component activated: "//trim(compname), &
                          ESMF_LOGMSG_INFO, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! access the internal state
      call ESMF_GridCompGetInternalState(wrt_Comp, wrap, rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      wrt_int_state => wrap%write_int_state

      call ESMF_VMGetCurrent(VM,rc=RC)

      mype = wrt_int_state%mype
      lprnt = mype == lead_write_task
!    print *,'in wrt run, mype=',mype,'lead_write_task=',lead_write_task
!
!-----------------------------------------------------------------------
!*** get current time and elapsed forecast time

      call ESMF_ClockGet(clock=CLOCK, currTime=CURRTIME, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_TimeGet(time=currTime,yy=date(1),mm=date(2),dd=date(3),h=date(4), &
                        m=date(5),s=date(6),rc=rc)
      wrt_int_state%fdate(7) = 1
      wrt_int_state%fdate(1:6) = date(1:6)

       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!    if(mype == lead_write_task) print *,'in wrt run, curr time=',date
!
      call ESMF_TimeGet(time=wrt_int_state%IO_BASETIME,yy=date(1),mm=date(2),dd=date(3),h=date(4), &
                        m=date(5),s=date(6),rc=rc)

       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!    print *,'in wrt run, io_baseline time=',date
!
      wrt_int_state%IO_CURRTIMEDIFF = CURRTIME-wrt_int_state%IO_BASETIME
!
      call ESMF_TimeIntervalGet(timeinterval=wrt_int_state%IO_CURRTIMEDIFF &
                                   ,h           =nf_hours               &  !<-- Hours of elapsed time
                                   ,m           =nf_minutes             &  !<-- Minutes of elapsed time
                                   ,s           =nseconds               &  !<-- Seconds of elapsed time
                                   ,sN          =nseconds_num           &  !<-- Numerator of fractional elapsed seconds
                                   ,sD          =nseconds_den           &  !<-- denominator of fractional elapsed seconds
                                   ,rc          =RC)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!       if (lprnt) print *,'in wrt run, nf_hours=',nf_hours,nf_minutes,nseconds, &
!         'nseconds_num=',nseconds_num,nseconds_den,'mype=',mype
!
      nf_seconds = nf_hours*3600+nf_minuteS*60+nseconds+real(nseconds_num)/real(nseconds_den)
      wrt_int_state%nfhour = nf_seconds/3600.
      nf_hours   = int(nf_seconds/3600.)
      if(mype == lead_write_task) print *,'in write grid comp, nf_hours=',nf_hours
      ! if iau_offset > nf_hours, don't write out anything
      if (nf_hours < 0) return

      nf_minutes = int((nf_seconds-nf_hours*3600.)/60.)
      nseconds   = int(nf_seconds-nf_hours*3600.-nf_minutes*60.)
!      if (nf_seconds-nf_hours*3600 > 0 .and. nsout > 0) then
      if (nsout > 0) then
        ndig = max(log10(nf_hours+0.5)+1., 3.)
        write(cform, '("(I",I1,".",I1,",A1,I2.2,A1,I2.2)")') ndig, ndig
        write(cfhour, cform) nf_hours,':',nf_minutes,':',nseconds
      else
        ndig = max(log10(nf_hours+0.5)+1., 3.)
        write(cform, '("(I",I1,".",I1,")")') ndig, ndig
        write(cfhour, cform) nf_hours
      endif
!
       if(lprnt) print *,'in wrt run, nf_hours=',nf_hours,nf_minutes,nseconds, &
                'nseconds_num=',nseconds_num,nseconds_den,' FBCount=',FBCount,' cfhour=',trim(cfhour)

!    if(lprnt) print *,'in wrt run, cfhour=',cfhour, &
!     print *,'in wrt run, cfhour=',cfhour, &
!        ' nf_seconds=',nf_seconds,wrt_int_state%nfhour

! access the time Attribute which is updated by the driver each time
      call ESMF_LogWrite("before Write component get time", ESMF_LOGMSG_INFO, rc=rc)
      call ESMF_AttributeGet(imp_state_write, convention="NetCDF", purpose="FV3", &
                              name="time", value=time, rc=rc)

       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_LogWrite("before Write component af get time", ESMF_LOGMSG_INFO, rc=rc)
!
!-----------------------------------------------------------------------
!*** loop on the files that need to write out
!-----------------------------------------------------------------------

      do i=1, FBCount
        call ESMF_LogWrite("before Write component get mirror file bundle", ESMF_LOGMSG_INFO, rc=rc)
        call ESMF_StateGet(imp_state_write, itemName="mirror_"//trim(fcstItemNameList(i)), &
                           fieldbundle=file_bundle, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        call ESMF_LogWrite("before Write component af get mirror file bundle", ESMF_LOGMSG_INFO, rc=rc)
!recover fields from cartesian vector and sfc pressure
        call recover_fields(file_bundle,rc)
      enddo
!
!-----------------------------------------------------------------------
!*** do post
!-----------------------------------------------------------------------
      lmask_fields = .false.
      if( wrt_int_state%write_dopost ) then
!
        wbeg = MPI_Wtime()
        if (trim(output_grid) == 'regional_latlon' .or. &
            trim(output_grid) == 'rotated_latlon'  .or. &
            trim(output_grid) == 'lambert_conformal') then

            !mask fields according to sfc pressure
            do nbdl=1, wrt_int_state%FBCount
              call ESMF_LogWrite("before mask_fields for wrt field bundle", ESMF_LOGMSG_INFO, rc=rc)
              call mask_fields(wrt_int_state%wrtFB(nbdl),rc)
              call ESMF_LogWrite("after mask_fields for wrt field bundle", ESMF_LOGMSG_INFO, rc=rc)
            enddo
            lmask_fields = .true.
        endif

        call inline_post_run(wrt_int_state, mype, wrt_mpi_comm, lead_write_task, &
                          nf_hours, nf_minutes,nseconds)
        wend = MPI_Wtime()
        if (lprnt) then
          write(*,'(A,F10.5,A,I4.2,A,I2.2)')' actual    inline post Time is ',wend-wbeg &
                     ,' at Fcst ',nf_hours,':',nf_minutes
            endif

      endif
!
!-----------------------------------------------------------------------
! ** now loop through output field bundle
!-----------------------------------------------------------------------

      if ( wrt_int_state%output_history ) then

        file_loop_all: do nbdl=1, wrt_int_state%FBCount
!
          if(step == 1) then
            file_bundle = wrt_int_state%wrtFB(nbdl)
          endif

          ! set default chunksizes for netcdf output
          ! (use MPI decomposition size).
          ! if chunksize parameter set to negative value,
          ! netcdf library default is used.
          if (output_file(nbdl)(1:6) == 'netcdf') then 
             if (ichunk2d == 0) then
                if( wrt_int_state%mype == 0 ) &
                  ichunk2d = wrt_int_state%lon_end-wrt_int_state%lon_start+1
                call mpi_bcast(ichunk2d,1,mpi_integer,0,wrt_mpi_comm,rc)
             endif
             if (jchunk2d == 0) then
                if( wrt_int_state%mype == 0 ) &
                  jchunk2d = wrt_int_state%lat_end-wrt_int_state%lat_start+1
                call mpi_bcast(jchunk2d,1,mpi_integer,0,wrt_mpi_comm,rc)
             endif
             if (ichunk3d == 0) then
                if( wrt_int_state%mype == 0 ) &
                  ichunk3d = wrt_int_state%lon_end-wrt_int_state%lon_start+1
                call mpi_bcast(ichunk3d,1,mpi_integer,0,wrt_mpi_comm,rc)
             endif
             if (jchunk3d == 0) then
                if( wrt_int_state%mype == 0 ) &
                  jchunk3d = wrt_int_state%lat_end-wrt_int_state%lat_start+1
                call mpi_bcast(jchunk3d,1,mpi_integer,0,wrt_mpi_comm,rc)
             endif
             if (kchunk3d == 0 .and. nbdl == 1) then
                if( wrt_int_state%mype == 0 )  then
                  call ESMF_FieldBundleGet(wrt_int_state%wrtFB(nbdl), grid=wrtgrid)
                  call ESMF_AttributeGet(wrtgrid, convention="NetCDF", purpose="FV3", &
                          attnestflag=ESMF_ATTNEST_OFF, name='pfull', &
                          itemCount=kchunk3d, rc=rc)
                endif
                call mpi_bcast(kchunk3d,1,mpi_integer,0,wrt_mpi_comm,rc)
             endif
             if (wrt_int_state%mype == 0) then
                print *,'ichunk2d,jchunk2d',ichunk2d,jchunk2d
                print *,'ichunk3d,jchunk3d,kchunk3d',ichunk3d,jchunk3d,kchunk3d
             endif
          endif

          if ( trim(output_file(nbdl)) == 'nemsio' ) then
             filename = trim(wrt_int_state%wrtFB_names(nbdl))//'f'//trim(cfhour)//'.nemsio'
          else
             filename = trim(wrt_int_state%wrtFB_names(nbdl))//'f'//trim(cfhour)//'.nc'
          endif
!           if(mype == lead_write_task) print *,'in wrt run,filename=',trim(filename)

!
! set the time Attribute on the grid to carry it into the lower levels
          call ESMF_FieldBundleGet(file_bundle, grid=fbgrid, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_AttributeSet(fbgrid, convention="NetCDF", purpose="FV3", &
                               name="time", value=real(wrt_int_state%nfhour,ESMF_KIND_R8), rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!*** write out grid bundle:
! Provide log message indicating which wrtComp is active
          call ESMF_LogWrite("before Write component before gridFB ", ESMF_LOGMSG_INFO, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if (trim(output_grid) == 'cubed_sphere_grid') then

            wbeg = MPI_Wtime()
            call ESMFproto_FieldBundleWrite(gridFB, filename=trim(filename),      &
                                            convention="NetCDF", purpose="FV3",   &
                                            status=ESMF_FILESTATUS_REPLACE,       &
                                            state=stateGridFB, comps=compsGridFB,rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            call ESMFproto_FieldBundleWrite(wrt_int_state%wrtFB(nbdl),                    &
                                           filename=trim(filename), convention="NetCDF",  &
                                           purpose="FV3", status=ESMF_FILESTATUS_OLD,     &
                                           timeslice=step, state=optimize(nbdl)%state,    &
                                           comps=optimize(nbdl)%comps, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            wend = MPI_Wtime()
            if (lprnt) then
              write(*,'(A,F10.5,A,I4.2,A,I2.2)')' actual    netcdf Write Time is ',wend-wbeg &
                     ,' at Fcst ',NF_HOURS,':',NF_MINUTES
            endif

          else if (trim(output_grid) == 'gaussian_grid') then

            if (trim(output_file(nbdl)) == 'nemsio') then

              wbeg = MPI_Wtime()
              call write_nemsio(file_bundle,trim(filename),nf_hours, nf_minutes, &
                                nseconds, nseconds_num, nseconds_den,nbdl, rc)
              wend = MPI_Wtime()
              if (lprnt) then
                write(*,'(A,F10.5,A,I4.2,A,I2.2)')' nemsio      Write Time is ',wend-wbeg  &
                      ,' at Fcst ',NF_HOURS,':',NF_MINUTES
              endif

            else if (trim(output_file(nbdl)) == 'netcdf') then

              wbeg = MPI_Wtime()
              call write_netcdf(file_bundle,wrt_int_state%wrtFB(nbdl),trim(filename), &
                               wrt_mpi_comm,wrt_int_state%mype,imo,jmo,&
                               ichunk2d,jchunk2d,ichunk3d,jchunk3d,kchunk3d,rc)
              wend = MPI_Wtime()
              if (lprnt) then
                write(*,'(A,F10.5,A,I4.2,A,I2.2)')' netcdf      Write Time is ',wend-wbeg  &
                        ,' at Fcst ',NF_HOURS,':',NF_MINUTES
              endif

            else if (trim(output_file(nbdl)) == 'netcdf_parallel') then

#ifdef NO_PARALLEL_NETCDF
              rc = ESMF_RC_NOT_IMPL
              print *,'netcdf_parallel not available on this machine'
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=__FILE__)) return
#endif
              wbeg = MPI_Wtime()
              call write_netcdf_parallel(file_bundle,wrt_int_state%wrtFB(nbdl),   &
                                         trim(filename), wrt_mpi_comm,wrt_int_state%mype,imo,jmo,&
                                         ichunk2d,jchunk2d,ichunk3d,jchunk3d,kchunk3d,rc)
              wend = MPI_Wtime()
              if (lprnt) then
                write(*,'(A,F10.5,A,I4.2,A,I2.2)')' parallel netcdf      Write Time is ',wend-wbeg  &
                        ,' at Fcst ',NF_HOURS,':',NF_MINUTES
              endif

            else if (trim(output_file(nbdl)) == 'netcdf_esmf') then

              wbeg = MPI_Wtime()
              call ESMFproto_FieldBundleWrite(gridFB, filename=trim(filename),    &
                                              convention="NetCDF", purpose="FV3", &
                                              status=ESMF_FILESTATUS_REPLACE, state=stateGridFB, comps=compsGridFB,rc=rc)

              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

              call ESMFproto_FieldBundleWrite(wrt_int_state%wrtFB(nbdl),                   &
                                             filename=trim(filename), convention="NetCDF", &
                                             purpose="FV3", status=ESMF_FILESTATUS_OLD,    &
                                             timeslice=step, state=optimize(nbdl)%state,   &
                                             comps=optimize(nbdl)%comps, rc=rc)

              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

              wend = MPI_Wtime()
              if (lprnt) then
                write(*,'(A,F10.5,A,I4.2,A,I2.2)')' netcdf_esmf Write Time is ',wend-wbeg &
                        ,' at Fcst ',NF_HOURS,':',NF_MINUTES
              endif
            endif

          else if (trim(output_grid) == 'global_latlon') then

            if (trim(output_file(nbdl)) == 'netcdf') then

              wbeg = MPI_Wtime()
              call write_netcdf(file_bundle,wrt_int_state%wrtFB(nbdl),trim(filename), &
                               wrt_mpi_comm,wrt_int_state%mype,imo,jmo,&
                               ichunk2d,jchunk2d,ichunk3d,jchunk3d,kchunk3d,rc)
              wend = MPI_Wtime()
              if (lprnt) then
                write(*,'(A,F10.5,A,I4.2,A,I2.2)')' netcdf      Write Time is ',wend-wbeg  &
                        ,' at Fcst ',NF_HOURS,':',NF_MINUTES
              endif

            else if (trim(output_file(nbdl)) == 'netcdf_parallel') then

#ifdef NO_PARALLEL_NETCDF
              rc = ESMF_RC_NOT_IMPL
              print *,'netcdf_parallel not available on this machine'
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=__FILE__)) return
#endif
              wbeg = MPI_Wtime()
              call write_netcdf_parallel(file_bundle,wrt_int_state%wrtFB(nbdl), &
                                         trim(filename), wrt_mpi_comm,wrt_int_state%mype,imo,jmo,&
                                         ichunk2d,jchunk2d,ichunk3d,jchunk3d,kchunk3d,rc)
              wend = MPI_Wtime()
              if (lprnt) then
                write(*,'(A,F10.5,A,I4.2,A,I2.2)')' parallel netcdf      Write Time is ',wend-wbeg  &
                        ,' at Fcst ',NF_HOURS,':',NF_MINUTES
              endif

            else ! unknown output_file

              call ESMF_LogWrite("wrt_run: Unknown output_file",ESMF_LOGMSG_ERROR,rc=RC)
              call ESMF_Finalize(endflag=ESMF_END_ABORT)

            endif

          else if (trim(output_grid) == 'regional_latlon' .or. &
                   trim(output_grid) == 'rotated_latlon'  .or. &
                 trim(output_grid) == 'lambert_conformal') then

            !mask fields according to sfc pressure
            !if (mype == lead_write_task) print *,'before mask_fields'
            if( .not. lmask_fields ) then
              wbeg = MPI_Wtime()
              call ESMF_LogWrite("before mask_fields for wrt field bundle", ESMF_LOGMSG_INFO, rc=rc)
              !call mask_fields(wrt_int_state%wrtFB(nbdl),rc)
              call mask_fields(file_bundle,rc)
              !if (mype == lead_write_task) print *,'after mask_fields'
              call ESMF_LogWrite("after mask_fields for wrt field bundle", ESMF_LOGMSG_INFO, rc=rc)
              wend = MPI_Wtime()
              if (mype == lead_write_task) then
                write(*,'(A,F10.5,A,I4.2,A,I2.2)')' mask_fields time is ',wend-wbeg
              endif
            endif

            if (trim(output_file(nbdl)) == 'netcdf' .and. nbits==0) then

              wbeg = MPI_Wtime()
              call write_netcdf(file_bundle,wrt_int_state%wrtFB(nbdl),trim(filename), &
                                wrt_mpi_comm,wrt_int_state%mype,imo,jmo,&
                                ichunk2d,jchunk2d,ichunk3d,jchunk3d,kchunk3d,rc)
              wend = MPI_Wtime()
              if (mype == lead_write_task) then
                write(*,'(A,F10.5,A,I4.2,A,I2.2)')' netcdf      Write Time is ',wend-wbeg  &
                        ,' at Fcst ',NF_HOURS,':',NF_MINUTES
              endif

            else if (trim(output_file(nbdl)) == 'netcdf_parallel' .and. nbits==0) then

#ifdef NO_PARALLEL_NETCDF
              rc = ESMF_RC_NOT_IMPL
              print *,'netcdf_parallel not available on this machine'
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=__FILE__)) return
#endif
              wbeg = MPI_Wtime()
              call write_netcdf_parallel(file_bundle,wrt_int_state%wrtFB(nbdl), &
                                         trim(filename), wrt_mpi_comm,wrt_int_state%mype,imo,jmo,&
                                         ichunk2d,jchunk2d,ichunk3d,jchunk3d,kchunk3d,rc)
              wend = MPI_Wtime()
              if (lprnt) then
                write(*,'(A,F10.5,A,I4.2,A,I2.2)')' parallel netcdf      Write Time is ',wend-wbeg  &
                        ,' at Fcst ',NF_HOURS,':',NF_MINUTES
              endif
            else ! unknown output_file

              if( nbits /= 0) then
                call ESMF_LogWrite("wrt_run: lossy compression is not supported for regional grids",ESMF_LOGMSG_ERROR,rc=RC)
                call ESMF_Finalize(endflag=ESMF_END_ABORT)
              else 
                call ESMF_LogWrite("wrt_run: Unknown output_file",ESMF_LOGMSG_ERROR,rc=RC)
                call ESMF_Finalize(endflag=ESMF_END_ABORT)
              endif

            endif

          else ! unknown output_grid

            call ESMF_LogWrite("wrt_run: Unknown output_grid",ESMF_LOGMSG_ERROR,rc=RC)
            call ESMF_Finalize(endflag=ESMF_END_ABORT)

        endif

      enddo file_loop_all

! end output history
    endif
!
!** write out log file
!
    if(mype == lead_write_task) then
      do n=701,900
        inquire(n,opened=OPENED)
        if(.not.opened)then
          nolog = n
          exit
        endif
      enddo
!
      open(nolog,file='logf'//trim(cfhour),form='FORMATTED')
        write(nolog,100)wrt_int_state%nfhour,idate(1:6)
100     format(' completed fv3gfs fhour=',f10.3,2x,6(i4,2x))
      close(nolog)
    endif
!
!
!-----------------------------------------------------------------------
!
      call ESMF_VMBarrier(VM, rc=rc)
!
      write_run_tim = MPI_Wtime() - tbeg
!
      IF (lprnt) THEN
        WRITE(*,'(A,F10.5,A,I4.2,A,I2.2)')' total            Write Time is ',write_run_tim  &
                 ,' at Fcst ',NF_HOURS,':',NF_MINUTES
      ENDIF
!
      IF(RC /= ESMF_SUCCESS) THEN
        WRITE(0,*)"FAIL: WRITE_RUN"
!     ELSE
!       WRITE(0,*)"PASS: WRITE_RUN"
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE wrt_run
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      subroutine wrt_finalize(wrt_comp, imp_state_write, exp_state_write, clock, rc)
!
!-----------------------------------------------------------------------
!***  finalize the write gridded component.
!-----------------------------------------------------------------------
!
!***  HISTORY
!       Feb 2017:  J. Wang  - deallocate for fv3
!
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp)            :: wrt_comp
      type(ESMF_State)               :: imp_state_write, exp_state_write
      type(ESMF_Clock)               :: clock
      integer,intent(out)            :: rc
!
!***  local variables
!
      integer :: stat
!
      type(write_wrap)               :: wrap
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  retrieve the write component's esmf internal state(used later for 
!***  post finalization)
!-----------------------------------------------------------------------
!
      call ESMF_GridCompGetInternalState(wrt_comp, wrap, rc)  
      deallocate(wrap%write_int_state,stat=stat)  
!
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
          msg="Deallocation of internal state memory failed.", &
          line=__LINE__, file=__FILE__)) return
!
!-----------------------------------------------------------------------
!
      IF(RC /= ESMF_SUCCESS)THEN
        WRITE(0,*)'FAIL: Write_Finalize.'
!      ELSE
!        WRITE(0,*)'PASS: Write_Finalize.'
      ENDIF
!
!-----------------------------------------------------------------------
!
    end subroutine wrt_finalize
!
!-----------------------------------------------------------------------
!
   subroutine recover_fields(file_bundle,rc)

     type(ESMF_FieldBundle), intent(in)              :: file_bundle
     integer,                intent(out),   optional :: rc
!
     integer i,j,k,ifld,fieldCount,nstt,nend,fieldDimCount,gridDimCount
     integer istart,iend,jstart,jend,kstart,kend,km
     logical uPresent, vPresent
     type(ESMF_Grid)  fieldGrid
     type(ESMF_Field)  ufield, vfield
     type(ESMF_TypeKind_Flag) typekind
     character(100) fieldName,uwindname,vwindname
     type(ESMF_Field),   allocatable  :: fcstField(:)
     real(ESMF_KIND_R8), dimension(:,:),     pointer  :: lon, lat
     real(ESMF_KIND_R8), dimension(:,:),     pointer  :: lonloc, latloc
     real(ESMF_KIND_R4), dimension(:,:),     pointer  :: pressfc
     real(ESMF_KIND_R4), dimension(:,:),     pointer  :: uwind2dr4,vwind2dr4
     real(ESMF_KIND_R4), dimension(:,:,:),   pointer  :: uwind3dr4,vwind3dr4
     real(ESMF_KIND_R4), dimension(:,:,:),   pointer  :: cart3dPtr2dr4
     real(ESMF_KIND_R4), dimension(:,:,:,:), pointer  :: cart3dPtr3dr4
     real(ESMF_KIND_R8), dimension(:,:,:,:), pointer  :: cart3dPtr3dr8
     save lonloc, latloc
     real(ESMF_KIND_R8) :: coslon, sinlon, sinlat
!
! get filed count
     call ESMF_FieldBundleGet(file_bundle, fieldCount=fieldCount, &
                              grid=fieldGrid, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
     CALL ESMF_LogWrite("call recover field on wrt comp",ESMF_LOGMSG_INFO,rc=RC)
     call ESMF_GridGet(fieldgrid, dimCount=gridDimCount, rc=rc)

     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

     if( first_getlatlon ) then
       CALL ESMF_LogWrite("call recover field get coord 1",ESMF_LOGMSG_INFO,rc=RC)

       call ESMF_GridGetCoord(fieldgrid, coordDim=1, farrayPtr=lon, rc=rc)

       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

       allocate(lonloc(lbound(lon,1):ubound(lon,1),lbound(lon,2):ubound(lon,2)))
       istart = lbound(lon,1)
       iend   = ubound(lon,1)
       jstart = lbound(lon,2)
       jend   = ubound(lon,2)
!$omp parallel do default(none) shared(lon,lonloc,jstart,jend,istart,iend) &
!$omp             private(i,j)
       do j=jstart,jend
        do i=istart,iend
          lonloc(i,j) = lon(i,j) * pi/180.
        enddo
       enddo

       CALL ESMF_LogWrite("call recover field get coord 2",ESMF_LOGMSG_INFO,rc=RC)

       call ESMF_GridGetCoord(fieldgrid, coordDim=2, farrayPtr=lat, rc=rc)

       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

       allocate(latloc(lbound(lat,1):ubound(lat,1),lbound(lat,2):ubound(lat,2)))
       istart = lbound(lat,1)
       iend   = ubound(lat,1)
       jstart = lbound(lat,2)
       jend   = ubound(lat,2)
!$omp parallel do default(none) shared(lat,latloc,jstart,jend,istart,iend) &
!$omp             private(i,j)
       do j=jstart,jend
        do i=istart,iend
          latloc(i,j) = lat(i,j) * pi/180.d0
        enddo
       enddo
       first_getlatlon = .false.
     endif
!
     allocate(fcstField(fieldCount))
     CALL ESMF_LogWrite("call recover field get fcstField",ESMF_LOGMSG_INFO,rc=RC)
     call ESMF_FieldBundleGet(file_bundle, fieldList=fcstField, itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
!
     do ifld=1,fieldCount

       CALL ESMF_LogWrite("call recover field get fieldname, type dimcount",ESMF_LOGMSG_INFO,rc=RC)
       call ESMF_FieldGet(fcstField(ifld),name=fieldName,typekind=typekind,dimCount=fieldDimCount, rc=rc)

! convert back wind
       if(index(trim(fieldName),"vector")>0) then
         nstt = index(trim(fieldName),"wind")+4
         nend = index(trim(fieldName),"vector")-1
         if( nend>nstt ) then
           uwindname = 'u'//fieldName(nstt+1:nend)
           vwindname = 'v'//fieldName(nstt+1:nend)
         else
           uwindname = 'ugrd'
           vwindname = 'vgrd'
         endif
!         print *,'in get 3D vector wind, uwindname=',trim(uwindname),' v=', trim(vwindname),' fieldname=',trim(fieldname)
! get u , v wind
         CALL ESMF_LogWrite("call recover field get u, v field",ESMF_LOGMSG_INFO,rc=RC)
         call ESMF_FieldBundleGet(file_bundle,trim(uwindname),field=ufield,isPresent=uPresent,rc=rc)
         call ESMF_FieldBundleGet(file_bundle,trim(vwindname),field=vfield,isPresent=vPresent,rc=rc)
         if(.not. uPresent .or. .not.vPresent) then
           rc=990
           print *,' ERROR ,the local wind is not present! rc=', rc
           exit
         endif

! get field data
         if ( typekind == ESMF_TYPEKIND_R4 ) then
           if( fieldDimCount > gridDimCount+1 ) then
             CALL ESMF_LogWrite("call recover field get 3d card wind farray",ESMF_LOGMSG_INFO,rc=RC)
             call ESMF_FieldGet(fcstField(ifld), localDe=0, farrayPtr=cart3dPtr3dr4, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
             if( ubound(cart3dPtr3dr4,1)-lbound(cart3dPtr3dr4,1)+1/=3) then
               rc=991
               print *,'ERROR, 3D the vector dimension /= 3, rc=',rc
               exit
             endif
             iend   = ubound(cart3dPtr3dr4,2)
             istart = lbound(cart3dPtr3dr4,2)
             iend   = ubound(cart3dPtr3dr4,2)
             jstart = lbound(cart3dPtr3dr4,3)
             jend   = ubound(cart3dPtr3dr4,3)
             kstart = lbound(cart3dPtr3dr4,4)
             kend   = ubound(cart3dPtr3dr4,4)
             call ESMF_FieldGet(ufield, localDe=0, farrayPtr=uwind3dr4,rc=rc)
             call ESMF_FieldGet(vfield, localDe=0, farrayPtr=vwind3dr4,rc=rc)
! update u , v wind
!$omp parallel do default(shared) private(i,j,k,coslon,sinlon,sinlat)
             do k=kstart,kend
!$omp parallel do default(none) shared(uwind3dr4,vwind3dr4,lonloc,latloc,cart3dPtr3dr4,jstart,jend,istart,iend,k) &
!$omp             private(i,j,coslon,sinlon,sinlat)
               do j=jstart, jend
                 do i=istart, iend
                  coslon = cos(lonloc(i,j))
                  sinlon = sin(lonloc(i,j))
                  sinlat = sin(latloc(i,j))
                  uwind3dr4(i,j,k) = cart3dPtr3dr4(1,i,j,k) * coslon           &
                                   + cart3dPtr3dr4(2,i,j,k) * sinlon
                  vwind3dr4(i,j,k) =-cart3dPtr3dr4(1,i,j,k) * sinlat*sinlon    &
                                   + cart3dPtr3dr4(2,i,j,k) * sinlat*coslon    &
                                   + cart3dPtr3dr4(3,i,j,k) * cos(latloc(i,j))
                 enddo
               enddo
             enddo
           else
             call ESMF_FieldGet(fcstField(ifld), localDe=0, farrayPtr=cart3dPtr2dr4, rc=rc)
               if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
             if( ubound(cart3dPtr2dr4,1)-lbound(cart3dPtr2dr4,1)+1 /= 3) then
               rc=991
               print *,'ERROR, 2D the vector dimension /= 3, rc=',rc
               exit
             endif
             istart = lbound(cart3dPtr2dr4,2)
             iend   = ubound(cart3dPtr2dr4,2)
             jstart = lbound(cart3dPtr2dr4,3)
             jend   = ubound(cart3dPtr2dr4,3)

             call ESMF_FieldGet(ufield, localDe=0, farrayPtr=uwind2dr4,rc=rc)
             call ESMF_FieldGet(vfield, localDe=0, farrayPtr=vwind2dr4,rc=rc)
              ! update u , v wind
!$omp parallel do default(none) shared(uwind2dr4,vwind2dr4,lonloc,latloc,cart3dPtr2dr4,jstart,jend,istart,iend) &
!$omp             private(i,j,k,coslon,sinlon,sinlat)
             do j=jstart, jend
               do i=istart, iend
                  coslon = cos(lonloc(i,j))
                  sinlon = sin(lonloc(i,j))
                  sinlat = sin(latloc(i,j))
                  uwind2dr4(i,j) = cart3dPtr2dr4(1,i,j) * coslon         &
                                 + cart3dPtr2dr4(2,i,j) * sinlon
                  vwind2dr4(i,j) =-cart3dPtr2dr4(1,i,j) * sinlat*sinlon  &
                                 + cart3dPtr2dr4(2,i,j) * sinlat*coslon  &
                                 + cart3dPtr2dr4(3,i,j) * cos(latloc(i,j))
               enddo
             enddo
           endif
         endif
!     print *,'in 3DCartesian2wind, uwindname=', trim(uwindname),'uPresent =',uPresent, 'vPresent=',vPresent,'fieldDimCount=', &
!       fieldDimCount,'gridDimCount=',gridDimCount

! convert back surface pressure
       else if(index(trim(fieldName),"pressfc")>0) then
         call ESMF_FieldGet(fcstField(ifld),localDe=0, farrayPtr=pressfc, rc=rc)
         istart = lbound(pressfc,1)
         iend   = ubound(pressfc,1)
         jstart = lbound(pressfc,2)
         jend   = ubound(pressfc,2)
!$omp parallel do default(none) shared(pressfc,jstart,jend,istart,iend) private(i,j)
         do j=jstart, jend
           do i=istart, iend
             pressfc(i,j) = pressfc(i,j)**(grav/(rdgas*stndrd_atmos_lapse))*stndrd_atmos_ps
           enddo
         enddo
       endif
     enddo
!
     deallocate(fcstField)
     rc = 0

   end subroutine recover_fields
!
!-----------------------------------------------------------------------
!
   subroutine mask_fields(file_bundle,rc)

     type(ESMF_FieldBundle), intent(in)              :: file_bundle
     integer,                intent(out),   optional :: rc
!
     integer i,j,k,ifld,fieldCount,nstt,nend,fieldDimCount,gridDimCount
     integer istart,iend,jstart,jend,kstart,kend,km
     type(ESMF_Grid)  fieldGrid
     type(ESMF_TypeKind_Flag) typekind
     type(ESMF_TypeKind_Flag) attTypeKind
     character(len=ESMF_MAXSTR) fieldName
     type(ESMF_Field),   allocatable                  :: fcstField(:)
     real(ESMF_KIND_R4), dimension(:,:),     pointer  :: var2dPtr2dr4
     real(ESMF_KIND_R4), dimension(:,:,:),   pointer  :: var3dPtr3dr4
     real(ESMF_KIND_R4), dimension(:,:,:),   pointer  :: vect3dPtr2dr4
     real(ESMF_KIND_R4), dimension(:,:,:,:), pointer  :: vect4dPtr3dr4
     real(ESMF_KIND_R4), dimension(:,:), allocatable  :: maskwrt

     logical :: mvispresent=.false.
     real(ESMF_KIND_R4) :: missing_value_r4=-1.e+10
     real(ESMF_KIND_R8) :: missing_value_r8=9.99e20
     character(len=ESMF_MAXSTR) :: msg

     save maskwrt

     call ESMF_LogWrite("call mask field on wrt comp",ESMF_LOGMSG_INFO,rc=RC)

! get fieldCount
     call ESMF_FieldBundleGet(file_bundle, fieldCount=fieldCount, &
         grid=fieldGrid, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out
! get gridDimCount
     call ESMF_GridGet(fieldgrid, dimCount=gridDimCount, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out

     allocate(fcstField(fieldCount))
     call ESMF_LogWrite("call mask field get fcstField",ESMF_LOGMSG_INFO,rc=RC)
     call ESMF_FieldBundleGet(file_bundle, fieldList=fcstField, itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)

! generate the maskwrt according to surface pressure
     if( first_getmaskwrt ) then

     do ifld=1,fieldCount
       !call ESMF_LogWrite("call mask field get fieldname, type dimcount",ESMF_LOGMSG_INFO,rc=RC)
       call ESMF_FieldGet(fcstField(ifld),name=fieldName,typekind=typekind,dimCount=fieldDimCount, rc=rc)
       !write(msg,*) 'fieldName,typekind,fieldDimCount=',trim(fieldName),typekind,fieldDimCount
       !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
       if (.not. allocated(maskwrt)) then
         if ( typekind == ESMF_TYPEKIND_R4 .and. fieldDimCount == gridDimCount) then
           call ESMF_FieldGet(fcstField(ifld),localDe=0, farrayPtr=var2dPtr2dr4, rc=rc)
           istart = lbound(var2dPtr2dr4,1)
           iend   = ubound(var2dPtr2dr4,1)
           jstart = lbound(var2dPtr2dr4,2)
           jend   = ubound(var2dPtr2dr4,2)
           allocate(maskwrt(istart:iend,jstart:jend))
           maskwrt(istart:iend,jstart:jend)=1.0
         endif
       endif
       if(index(trim(fieldName),"pressfc")>0) then
         call ESMF_FieldGet(fcstField(ifld),localDe=0, farrayPtr=var2dPtr2dr4, rc=rc)
         istart = lbound(var2dPtr2dr4,1)
         iend   = ubound(var2dPtr2dr4,1)
         jstart = lbound(var2dPtr2dr4,2)
         jend   = ubound(var2dPtr2dr4,2)
         if (.not. allocated(maskwrt)) then
           allocate(maskwrt(istart:iend,jstart:jend))
           maskwrt(istart:iend,jstart:jend)=1.0
         endif
!$omp parallel do default(shared) private(i,j)
         do j=jstart, jend
           do i=istart, iend
             if(abs(var2dPtr2dr4(i,j)-0.) < 1.0e-6) maskwrt(i,j)=0.
           enddo
         enddo
         call ESMF_LogWrite("call mask field pressfc found, maskwrt generated",ESMF_LOGMSG_INFO,rc=RC)
         exit
       endif
     enddo
     first_getmaskwrt = .false.

     endif !first_getmaskwrt

! loop to mask all fields according to maskwrt
     do ifld=1,fieldCount
       !call ESMF_LogWrite("call mask field get fieldname, type dimcount",ESMF_LOGMSG_INFO,rc=RC)
       call ESMF_FieldGet(fcstField(ifld),name=fieldName,typekind=typekind,dimCount=fieldDimCount, rc=rc)
       !write(msg,*) 'fieldName,typekind,fieldDimCount=',trim(fieldName),typekind,fieldDimCount
       !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
       ! For vector fields
       if(index(trim(fieldName),"vector")>0) then
         ! Only work on ESMF_TYPEKIND_R4 fields for now
         if ( typekind == ESMF_TYPEKIND_R4 ) then
           ! 3-d vector fields with 4-d arrays
           if( fieldDimCount > gridDimCount+1 ) then
             !call ESMF_LogWrite("call mask field get vector 3d farray",ESMF_LOGMSG_INFO,rc=RC)
             call ESMF_FieldGet(fcstField(ifld), localDe=0, farrayPtr=vect4dPtr3dr4, rc=rc)
               if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                 line=__LINE__, file=__FILE__)) return  ! bail out
             if( ubound(vect4dPtr3dr4,1)-lbound(vect4dPtr3dr4,1)+1/=3 ) then
               rc=991
               print *,'ERROR, 3D the vector dimension /= 3, rc=',rc
               exit
             endif
             ! Get the _FillValue from the field attribute if exists
             call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                 name="_FillValue", typekind=attTypeKind, isPresent=mvispresent, rc=rc)
             !write(msg,*) 'fieldName,attTypeKind,isPresent=',trim(fieldName),attTypeKind,mvispresent
             !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
             if ( mvispresent ) then
               if (attTypeKind==ESMF_TYPEKIND_R4) then
                 call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                        name="_FillValue", value=missing_value_r4, isPresent=mvispresent, rc=rc)
                 !write(msg,*) 'fieldName,_FillValue,isPresent=',trim(fieldName),missing_value_r4,mvispresent
                 !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
               else if (attTypeKind==ESMF_TYPEKIND_R8) then
                 call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                        name="_FillValue", value=missing_value_r8, isPresent=mvispresent, rc=rc)
                 !write(msg,*) 'fieldName,_FillValue,isPresent=',trim(fieldName),missing_value_r8,mvispresent
                 !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
               endif
               istart = lbound(vect4dPtr3dr4,2)
               iend   = ubound(vect4dPtr3dr4,2)
               jstart = lbound(vect4dPtr3dr4,3)
               jend   = ubound(vect4dPtr3dr4,3)
               kstart = lbound(vect4dPtr3dr4,4)
               kend   = ubound(vect4dPtr3dr4,4)
!$omp parallel do default(shared) private(i,j,k)
               do k=kstart,kend
                 do j=jstart, jend
                   do i=istart, iend
                     if (maskwrt(i,j)<1.0 .and. attTypeKind==ESMF_TYPEKIND_R4) vect4dPtr3dr4(:,i,j,k)=missing_value_r4
                     if (maskwrt(i,j)<1.0 .and. attTypeKind==ESMF_TYPEKIND_R8) vect4dPtr3dr4(:,i,j,k)=missing_value_r8
                   enddo
                 enddo
               enddo
             endif !mvispresent
           ! 2-d vector fields with 3-d arrays
           else
             call ESMF_FieldGet(fcstField(ifld), localDe=0, farrayPtr=vect3dPtr2dr4, rc=rc)
               if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                 line=__LINE__, file=__FILE__)) return  ! bail out
             if( ubound(vect3dPtr2dr4,1)-lbound(vect3dPtr2dr4,1)+1 /= 3 ) then
               rc=991
               print *,'ERROR, 2D the vector dimension /= 3, rc=',rc
               exit
             endif
             ! Get the _FillValue from the field attribute if exists
             call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                 name="_FillValue", typekind=attTypeKind, isPresent=mvispresent, rc=rc)
             !write(msg,*) 'fieldName,attTypeKind,isPresent=',trim(fieldName),attTypeKind,mvispresent
             !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
             if ( mvispresent ) then
               if (attTypeKind==ESMF_TYPEKIND_R4) then
                 call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                        name="_FillValue", value=missing_value_r4, isPresent=mvispresent, rc=rc)
                 !write(msg,*) 'fieldName,_FillValue,isPresent=',trim(fieldName),missing_value_r4,mvispresent
                 !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
               else if (attTypeKind==ESMF_TYPEKIND_R8) then
                 call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                        name="_FillValue", value=missing_value_r8, isPresent=mvispresent, rc=rc)
                 !write(msg,*) 'fieldName,_FillValue,isPresent=',trim(fieldName),missing_value_r8,mvispresent
                 !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
               endif
               istart = lbound(vect3dPtr2dr4,2)
               iend   = ubound(vect3dPtr2dr4,2)
               jstart = lbound(vect3dPtr2dr4,3)
               jend   = ubound(vect3dPtr2dr4,3)
!$omp parallel do default(shared) private(i,j)
               do j=jstart, jend
                 do i=istart, iend
                   if (maskwrt(i,j)<1.0 .and. attTypeKind==ESMF_TYPEKIND_R4) vect3dPtr2dr4(:,i,j)=missing_value_r4
                   if (maskwrt(i,j)<1.0 .and. attTypeKind==ESMF_TYPEKIND_R8) vect3dPtr2dr4(:,i,j)=missing_value_r8
                 enddo
               enddo
             endif !mvispresent
           endif
         endif
! For non-vector fields
       else
         ! Only work on ESMF_TYPEKIND_R4 fields for now
         if ( typekind == ESMF_TYPEKIND_R4 ) then
           ! 2-d fields
           if(fieldDimCount == gridDimCount) then
             call ESMF_FieldGet(fcstField(ifld),localDe=0, farrayPtr=var2dPtr2dr4, rc=rc)
             ! Get the _FillValue from the field attribute if exists
             call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                 name="_FillValue", typekind=attTypeKind, isPresent=mvispresent, rc=rc)
             !write(msg,*) 'fieldName,attTypeKind,isPresent=',trim(fieldName),attTypeKind,mvispresent
             !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
             if ( mvispresent ) then
               if (attTypeKind==ESMF_TYPEKIND_R4) then
                 call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                        name="_FillValue", value=missing_value_r4, isPresent=mvispresent, rc=rc)
                 !write(msg,*) 'fieldName,_FillValue,isPresent=',trim(fieldName),missing_value_r4,mvispresent
                 !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
               else if (attTypeKind==ESMF_TYPEKIND_R8) then
                 call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                        name="_FillValue", value=missing_value_r8, isPresent=mvispresent, rc=rc)
                 !write(msg,*) 'fieldName,_FillValue,isPresent=',trim(fieldName),missing_value_r8,mvispresent
                 !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
               endif
               istart = lbound(var2dPtr2dr4,1)
               iend   = ubound(var2dPtr2dr4,1)
               jstart = lbound(var2dPtr2dr4,2)
               jend   = ubound(var2dPtr2dr4,2)
!$omp parallel do default(shared) private(i,j)
               do j=jstart, jend
                 do i=istart, iend
                   if (maskwrt(i,j)<1.0 .and. attTypeKind==ESMF_TYPEKIND_R4) var2dPtr2dr4(i,j)=missing_value_r4
                   if (maskwrt(i,j)<1.0 .and. attTypeKind==ESMF_TYPEKIND_R8) var2dPtr2dr4(i,j)=missing_value_r8
                 enddo
               enddo
             endif !mvispresent
           ! 3-d fields
           else if(fieldDimCount == gridDimCount+1) then
             call ESMF_FieldGet(fcstField(ifld),localDe=0, farrayPtr=var3dPtr3dr4, rc=rc)
             ! Get the _FillValue from the field attribute if exists
             call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                 name="_FillValue", typekind=attTypeKind, isPresent=mvispresent, rc=rc)
             !write(msg,*) 'fieldName,attTypeKind,isPresent=',trim(fieldName),attTypeKind,mvispresent
             !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
             if ( mvispresent ) then
               if (attTypeKind==ESMF_TYPEKIND_R4) then
                 call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                        name="_FillValue", value=missing_value_r4, isPresent=mvispresent, rc=rc)
                 !write(msg,*) 'fieldName,_FillValue,isPresent=',trim(fieldName),missing_value_r4,mvispresent
                 !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
               else if (attTypeKind==ESMF_TYPEKIND_R8) then
                 call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                        name="_FillValue", value=missing_value_r8, isPresent=mvispresent, rc=rc)
                 !write(msg,*) 'fieldName,_FillValue,isPresent=',trim(fieldName),missing_value_r8,mvispresent
                 !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
               endif
               istart = lbound(var3dPtr3dr4,1)
               iend   = ubound(var3dPtr3dr4,1)
               jstart = lbound(var3dPtr3dr4,2)
               jend   = ubound(var3dPtr3dr4,2)
               kstart = lbound(var3dPtr3dr4,3)
               kend   = ubound(var3dPtr3dr4,3)
!$omp parallel do default(shared) private(i,j,k)
               do k=kstart,kend
                 do j=jstart, jend
                   do i=istart, iend
                     if (maskwrt(i,j)<1.0 .and. attTypeKind==ESMF_TYPEKIND_R4) var3dPtr3dr4(i,j,k)=missing_value_r4
                     if (maskwrt(i,j)<1.0 .and. attTypeKind==ESMF_TYPEKIND_R8) var3dPtr3dr4(i,j,k)=missing_value_r8
                   enddo
                 enddo
               enddo
             endif !mvispresent
           endif
         endif
       endif
     enddo
!
     deallocate(fcstField)
     rc = 0

   end subroutine mask_fields
!
!-----------------------------------------------------------------------
!

  subroutine ESMFproto_FieldBundleWrite(fieldbundle, fileName, &
    convention, purpose, status, timeslice, state, comps, rc)
    type(ESMF_FieldBundle),     intent(in)              :: fieldbundle
    character(*),               intent(in)              :: fileName
    character(*),               intent(in),    optional :: convention
    character(*),               intent(in),    optional :: purpose
    type(ESMF_FileStatus_Flag), intent(in),    optional :: status
    integer,                    intent(in),    optional :: timeslice
    type(ESMF_State),           intent(inout), optional :: state
    type(ESMF_GridComp), allocatable, intent(inout), optional :: comps(:)
    integer,                    intent(out),   optional :: rc
    
    ! Prototype multi-tile implementation for FieldBundleWrite().
    ! Produces as many output files as there are tiles. The naming of the 
    ! output files is such that the string in the fileName argument is used
    ! as the basis. If fileName ends with ".nc", then this suffix is replaced
    ! by ".tileN.nc", where "N" is the tile number. If fileName does not 
    ! end in ".nc", then ".tileN.nc" will simply be appended.
    !
    ! Restrictions:
    !   - All Fields in the FieldBundle must have the same tileCount
    
    integer                             :: i, j, ind
    integer                             :: fieldCount, tileCount, itemCount
    type(ESMF_Field), allocatable       :: fieldList(:), tileFieldList(:)
    type(ESMF_Grid)                     :: grid
    type(ESMF_Array)                    :: array
    type(ESMF_DistGrid)                 :: distgrid
    type(ESMF_DELayout)                 :: delayout
    type(ESMF_FieldBundle)              :: wrtTileFB
    type(ESMF_FieldBundle), allocatable :: wrtTileFBList(:)
    character(len=80), allocatable      :: itemNameList(:)
    logical                             :: stateIsEmpty
    character(len=80)                   :: tileFileName
    character(len=80)                   :: statusStr
    integer, allocatable                :: petList(:)
    character(1024)                     :: msgString
    type(ESMF_State), allocatable       :: ioState(:)
    integer                             :: timesliceOpt
    integer                             :: urc

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_VMLogMemInfo("Entering ESMFproto_FieldBundleWrite",rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! query number of fields in fieldbundle
    call ESMF_FieldBundleGet(fieldbundle, fieldCount=fieldCount, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! early successful exit if there are no fields present
    if (fieldCount==0) return
    ! obtain list of fields in the fieldbundle

    allocate(fieldList(fieldCount), tileFieldList(fieldCount))
    call ESMF_FieldBundleGet(fieldbundle, fieldList=fieldList, &
                             itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    ! determine tileCount by looking at first field

    call ESMF_FieldGet(fieldList(1), array=array, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_ArrayGet(array, tileCount=tileCount, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      
    ! deal with optional state argument
    stateIsEmpty = .true.
    if (present(state)) then
      if (.not.ESMF_StateIsCreated(state, rc=rc)) then
        state = ESMF_StateCreate(rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      endif
      call ESMF_StateGet(state, itemCount=itemCount, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if (itemCount /= 0) then
        stateIsEmpty = .false.
        if (itemCount /= tileCount) then
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
                                msg="Number of items in state must match number of tiles.", &
                                line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
        endif
        allocate(itemNameList(itemCount),wrtTileFBList(itemCount))
        call ESMF_StateGet(state, itemNameList=itemNameList, &
                           itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        do i=1, itemCount
          call ESMF_StateGet(state, itemName=itemNameList(i), &
                             fieldbundle=wrtTileFBList(i), rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        enddo
      endif
    endif
    
    ! loop over all the tiles and construct a tile specific fieldbundle
    call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() before tileCount-loop",&
      ESMF_LOGMSG_INFO, rc=rc)
      
    !TODO: remove this once comps is hidden within state
    if (present(comps)) then
      allocate(ioState(tileCount))
      if (.not.allocated(comps)) then
        ! first-write
        allocate(comps(tileCount))
        do i=1, tileCount
          call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() before "// &
                             "ESMFproto_FieldMakeSingleTile() w/ petList", &
                              ESMF_LOGMSG_INFO, rc=rc)
          do j=1, fieldCount
            ! access only tile specific part of field
            call ESMFproto_FieldMakeSingleTile(fieldList(j), tile=i, &
                                               tileField=tileFieldList(j), petList=petList, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          enddo
!          write(msgString, *) petList
!          call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() after "// &
!            "ESMFproto_FieldMakeSingleTile(), petList:"//trim(msgString), &
!            ESMF_LOGMSG_INFO, rc=rc)
          ! create component to handle this tile I/O
          call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() before "// &
                             "tile-component creation", ESMF_LOGMSG_INFO, rc=rc)
          comps(i) = ESMF_GridCompCreate(petList=petList, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          ! convention
          call ESMF_AttributeSet(comps(i), name="convention", &
                                 value=convention, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          ! purpose
          call ESMF_AttributeSet(comps(i), name="purpose", &
                                 value=purpose, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          ! timeslice
          timesliceOpt = -1 ! init
          if (present(timeslice)) timesliceOpt = timeslice
          call ESMF_AttributeSet(comps(i), name="timeslice", &
                                 value=timesliceOpt, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          ! status
          statusStr = "ESMF_FILESTATUS_UNKNOWN" ! default
          if (present(status)) then
            if (status==ESMF_FILESTATUS_UNKNOWN) then
              statusStr="ESMF_FILESTATUS_UNKNOWN" ! default
            else if (status==ESMF_FILESTATUS_NEW) then
              statusStr="ESMF_FILESTATUS_NEW" ! default
            else if (status==ESMF_FILESTATUS_OLD) then
              statusStr="ESMF_FILESTATUS_OLD" ! default
            else if (status==ESMF_FILESTATUS_REPLACE) then
              statusStr="ESMF_FILESTATUS_REPLACE" ! default
            endif
          endif
          call ESMF_AttributeSet(comps(i), name="status", &
                                 value=statusStr, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_GridCompSetServices(comps(i), ioCompSS, userRc=urc, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() after "// &
                             "tile-component creation", ESMF_LOGMSG_INFO, rc=rc)
        enddo
      endif
    endif
      
    do i=1, tileCount
      if (stateIsEmpty) then
        ! loop over all the fields and add tile specific part to fieldbundle
        call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() before "// &
                           "ESMFproto_FieldMakeSingleTile()", ESMF_LOGMSG_INFO, rc=rc)
        do j=1, fieldCount
          ! access only tile specific part of field
          call ESMFproto_FieldMakeSingleTile(fieldList(j), tile=i, &
                                             tileField=tileFieldList(j), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        enddo
        call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() after "// &
                           "ESMFproto_FieldMakeSingleTile()", ESMF_LOGMSG_INFO, rc=rc)
        ! create tile specific fieldbundle
        wrtTileFB = ESMF_FieldBundleCreate(fieldList=tileFieldList, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        ! ensure global attributes on the fieldbundle are passed on by reference
        call ESMF_AttributeCopy(fieldbundle, wrtTileFB, &
                                attcopy=ESMF_ATTCOPY_REFERENCE, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        ! store this fieldbundle in state if present
        if (present(state)) then
          call ESMF_StateAdd(state, fieldbundleList=(/wrtTileFB/), rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        endif
      else
        ! state brought in existing fieldbundles
        wrtTileFB = wrtTileFBList(i)
      endif
      ! write out the tile specific fieldbundle
      if(tileCount>1) then
        ind=min(index(trim(fileName),".nc",.true.)-1,len_trim(fileName))
        write(tileFileName, "(A,A,I1,A)") fileName(1:ind), ".tile", i, ".nc"
      else
        tileFileName=trim(fileName)
      endif
      if (present(comps)) then
        ioState(i) = ESMF_StateCreate(rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_StateAdd(ioState(i), (/wrtTileFB/), rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_AttributeSet(comps(i), name="tileFileName", value=tileFileName, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_GridCompRun(comps(i), importState=ioState(i), userRc=urc, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_StateDestroy(ioState(i), rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__))  return
      else
        call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() before "//      &
                           "ESMF_FieldBundleWrite(): "//trim(tileFileName), &
                            ESMF_LOGMSG_INFO, rc=rc)
        call ESMF_FieldBundleWrite(fieldbundle=wrtTileFB, fileName=tileFileName,          &
                                   convention=convention, purpose=purpose, status=status, &
                                   timeslice=timeslice, overwrite=.true., rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() after "// &
                           "ESMF_FieldBundleWrite()", ESMF_LOGMSG_INFO, rc=rc)
      endif
      if (.not.present(state)) then
        ! local garbage collection of fields
        do j=1, fieldCount
          call ESMF_FieldGet(tileFieldList(j), array=array, grid=grid, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_FieldDestroy(tileFieldList(j), noGarbage=.true., rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_GridDestroy(grid, noGarbage=.true., rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_ArrayGet(array, distgrid=distgrid, delayout=delayout, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_ArrayDestroy(array, noGarbage=.true., rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_DistGridDestroy(distgrid, noGarbage=.true., rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_DELayoutDestroy(delayout, noGarbage=.true., rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        enddo
        ! destroy tile specific fieldbundle
        call ESMF_FieldBundleDestroy(wrtTileFB, noGarbage=.true., rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      endif
    enddo
    if (present(comps)) then
      deallocate(ioState)
    endif
    call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() after tileCount-loop",&
                        ESMF_LOGMSG_INFO, rc=rc)
    
    ! deallocate temporary lists
    deallocate(fieldList, tileFieldList)

    call ESMF_VMLogMemInfo("Exiting ESMFproto_FieldBundleWrite",rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine ESMFproto_FieldBundleWrite

  !-----------------------------------------------------------------------------

  subroutine ioCompSS(comp, rc)
    type(ESMF_GridComp)   :: comp
    integer, intent(out)  :: rc

    rc = ESMF_SUCCESS
  
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
                                    userRoutine=ioCompRun, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine ioCompRun(comp, importState, exportState, clock, rc)
    use netcdf
    type(ESMF_GridComp)   :: comp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc

    type(ESMF_FieldBundle)           :: wrtTileFB
    character(len=80)                :: tileFileName
    character(len=80)                :: convention
    character(len=80)                :: purpose
    integer                          :: timeslice
    character(len=80)                :: statusStr
    type(ESMF_FileStatus_Flag)       :: status
    character(len=80)                :: itemNameList(1)

    integer                          :: localPet, i, j, k, ind
    type(ESMF_Grid)                  :: grid
    real(ESMF_KIND_R4), allocatable  :: valueListr4(:)
    real(ESMF_KIND_R8), allocatable  :: valueListr8(:)
    integer                          :: valueCount, fieldCount, udimCount
    character(80),      allocatable  :: udimList(:)
    integer                          :: ncerr, ncid, dimid, varid
    type(ESMF_Field),   allocatable  :: fieldList(:)
    type(ESMF_Field)                 :: field
    logical                          :: isPresent
    real(ESMF_KIND_R8)               :: time
    integer                          :: itemCount, attCount
    character(len=80),  allocatable  :: attNameList(:)
    character(len=80)                :: attName
    type(ESMF_TypeKind_Flag)         :: typekind
    character(len=80)                :: valueS
    integer                          :: valueI4
    real(ESMF_KIND_R4)               :: valueR4
    real(ESMF_KIND_R8)               :: valueR8
    logical                          :: thereAreVerticals

    rc = ESMF_SUCCESS

    ! Access the FieldBundle
    call ESMF_StateGet(importState, itemNameList=itemNameList, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_StateGet(importState, itemName=itemNameList(1), &
                       fieldBundle=wrtTileFB, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Access attributes on the component and use as parameters for Write()
    call ESMF_AttributeGet(comp, name="tileFileName", value=tileFileName, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeGet(comp, name="convention", value=convention, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeGet(comp, name="purpose", value=purpose, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeGet(comp, name="timeslice", value=timeslice, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeGet(comp, name="status", value=statusStr, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (trim(statusStr) == "ESMF_FILESTATUS_UNKNOWN") then
      status = ESMF_FILESTATUS_UNKNOWN
    else if (trim(statusStr) == "ESMF_FILESTATUS_NEW") then
      status = ESMF_FILESTATUS_NEW
    else if (trim(statusStr) == "ESMF_FILESTATUS_OLD") then
      status=ESMF_FILESTATUS_OLD
    else if (trim(statusStr) == "ESMF_FILESTATUS_REPLACE") then
      status = ESMF_FILESTATUS_REPLACE
    endif

    call ESMF_LogWrite("In ioCompRun() before writing to: "// &
                       trim(tileFileName), ESMF_LOGMSG_INFO, rc=rc)

    if (status == ESMF_FILESTATUS_OLD) then
      ! This writes the vectical coordinates and the time dimension into the
      ! file. Doing this before the large data sets are written, assuming that
      ! the first time coming into ioCompRun() with this tileFileName, only
      ! the grid info is written. Second time in, with ESMF_FILESTATUS_OLD,
      ! the large data sets are written. That is when vertical and time info
      ! is also written.
      ! Hoping for better performance because file exists, but is still small.

      call ESMF_GridCompGet(comp, localPet=localPet, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      if (localPet==0) then
        ! do this work only on the root pet
        call ESMF_FieldBundleGet(wrtTileFB, grid=grid, fieldCount=fieldCount, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        allocate(fieldList(fieldCount))
        call ESMF_FieldBundleGet(wrtTileFB, fieldList=fieldList, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        ! open this tile's NetCDF file
        ncerr = nf90_open(tileFileName, NF90_WRITE, ncid=ncid)

        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
        ! loop over all the fields in the bundle and handle their vectical dims

        thereAreVerticals = .false.
        do i=1, fieldCount
          field = fieldList(i)
          call ESMF_AttributeGetAttPack(field, convention="NetCDF", purpose="FV3", &
                                        isPresent=isPresent, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if (.not.isPresent) cycle ! field does not have the AttPack
          call ESMF_AttributeGet(field, convention="NetCDF", purpose="FV3", &
                                 name="ESMF:ungridded_dim_labels", isPresent=isPresent, &
            itemCount=udimCount, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if (udimCount==0 .or. .not.isPresent) cycle ! nothing there to do

          thereAreVerticals = .true.
          allocate(udimList(udimCount))
          call ESMF_AttributeGet(field, convention="NetCDF", purpose="FV3", &
                                 name="ESMF:ungridded_dim_labels", valueList=udimList, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          ! loop over all ungridded dimension labels
          do k=1, udimCount
            call write_out_ungridded_dim_atts(dimLabel=trim(udimList(k)), rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          enddo
          deallocate(udimList)
        enddo ! fieldCount
        deallocate(fieldList)
        if (thereAreVerticals) then
          ! see if the vertical_dim_labels attribute exists on the grid, and
          ! if so access it and write out vecticals accordingly
          call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                 name="vertical_dim_labels", isPresent=isPresent, &
                                 itemCount=udimCount, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if (isPresent .and. (udimCount>0) ) then
            allocate(udimList(udimCount))
            call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                   name="vertical_dim_labels", valueList=udimList, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ! loop over all ungridded dimension labels
            do k=1, udimCount
              call write_out_ungridded_dim_atts(dimLabel=trim(udimList(k)), rc=rc)

              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            enddo
            deallocate(udimList)
          endif
        endif
        ! inquire if NetCDF file already contains the "time" variable
        ncerr = nf90_inq_varid(ncid, "time", varid=varid)
        if (ncerr /= NF90_NOERR) then
          ! the variable does not exist in the NetCDF file yet -> add it
          ! access the "time" attribute on the grid
          call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                 name="time", value=time, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          ncerr = nf90_redef(ncid=ncid)
          if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          ncerr = nf90_inq_dimid(ncid, "time", dimid=dimid)
          if (ncerr /= NF90_NOERR) then
            ! "time" dimension does not yet exist, define as unlimited dim
            ncerr = nf90_def_dim(ncid, "time", NF90_UNLIMITED, dimid=dimid)

            if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
          endif
          ncerr = nf90_def_var(ncid, "time", NF90_DOUBLE, &
                               dimids=(/dimid/), varid=varid)
          if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          ncerr = nf90_enddef(ncid=ncid)

          if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          ncerr = nf90_put_var(ncid, varid, values=time)

          if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          ! loop over all the grid attributes that start with "time:", and
          ! put them on the "time" variable in the NetCDF file

          call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                 name="TimeAttributes", itemCount=itemCount, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if (itemCount > 0) then
            ncerr = nf90_redef(ncid=ncid)

            if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

            allocate(attNameList(itemCount))
            call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                   name="TimeAttributes", valueList=attNameList, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            do i=1, itemCount
              attName = attNameList(i)
              call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                     name=trim(attNameList(i)), typekind=typekind, rc=rc)
!                print *,'in esmf call, att name=',trim(attNameList(i))

              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

              if (typekind==ESMF_TYPEKIND_CHARACTER) then
                call ESMF_AttributeGet(grid,                               &
                                       convention="NetCDF", purpose="FV3", &
                                       name=trim(attNameList(i)), value=valueS, rc=rc)
!                print *,'in esmf call, att string value=',trim(valueS)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

                ncerr = nf90_put_att(ncid, varid, &
                                     trim(attName(6:len(attName))), values=valueS)

                if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

              else if (typekind==ESMF_TYPEKIND_I4) then

                call ESMF_AttributeGet(grid,                               &
                                       convention="NetCDF", purpose="FV3", &
                                       name=trim(attNameList(i)), value=valueI4, rc=rc)
!                print *,'in esmf call, att I4 value=',valueR8
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
                ncerr = nf90_put_att(ncid, varid, &
                                     trim(attName(6:len(attName))), values=valueI4)

                if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

              else if (typekind==ESMF_TYPEKIND_R4) then
                call ESMF_AttributeGet(grid,                               &
                                       convention="NetCDF", purpose="FV3", &
                                       name=trim(attNameList(i)), value=valueR4, rc=rc)
!                print *,'in esmf call, att r4 value=',valueR8

                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

                ncerr = nf90_put_att(ncid, varid, &
                                     trim(attName(6:len(attName))), values=valueR4)

                if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

              else if (typekind==ESMF_TYPEKIND_R8) then
                call ESMF_AttributeGet(grid,                               &
                                       convention="NetCDF", purpose="FV3", &
                                       name=trim(attNameList(i)), value=valueR8, rc=rc)
!                print *,'in esmf call, att r8 value=',valueR8

                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
                ncerr = nf90_put_att(ncid, varid, &
                                     trim(attName(6:len(attName))), values=valueR8)

                if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
              endif
            enddo
            deallocate(attNameList)
            ncerr = nf90_enddef(ncid=ncid)

            if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
          endif
        endif
        ! close the NetCDF file
        ncerr = nf90_close(ncid=ncid)

        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      endif
      call ESMF_LogWrite("In ioCompRun() after "// &
                         "writing vectical and time dimensions.", ESMF_LOGMSG_INFO, rc=rc)
    endif

    !TODO: remove this block once the ESMF_FieldBundleWrite() below allows to
    !TODO: specify the NetCDF 64-bit-offset file format.
    if (status==ESMF_FILESTATUS_REPLACE) then
      ! First time in with this filename (therefore 'replace'), create the
      ! file with 64bit-offset format in order to accommodate larger data
      ! volume.
      call ESMF_GridCompGet(comp, localPet=localPet, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      if (localPet==0) then
        ! only single PET to deal with NetCDF
        ncerr = nf90_create(tileFileName, &
                            cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncid)

        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        ncerr = nf90_close(ncid=ncid)

        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      endif
      status = ESMF_FILESTATUS_OLD  ! switch status to 'OLD' to not overwrite
      call ESMF_LogWrite("In ioCompRun() after creating the NetCDF file", &
                         ESMF_LOGMSG_INFO, rc=rc)
    endif

    if (timeslice==-1) then
      call ESMF_FieldBundleWrite(fieldbundle=wrtTileFB, fileName=tileFileName,          &
                                 convention=convention, purpose=purpose, status=status, &
                                 overwrite=.true., rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    else
      call ESMF_FieldBundleWrite(fieldbundle=wrtTileFB, fileName=tileFileName, &
                                 convention=convention, purpose=purpose, status=status, &
                                 timeslice=timeslice, overwrite=.true., rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif
    call ESMF_LogWrite("In ioCompRun() after "// &
                       "ESMF_FieldBundleWrite()", ESMF_LOGMSG_INFO, rc=rc)

  contains

    subroutine write_out_ungridded_dim_atts(dimLabel, rc)
      character(len=*)      :: dimLabel
      integer, intent(out)  :: rc

      ! inquire if NetCDF file already contains this ungridded dimension
      ncerr = nf90_inq_varid(ncid, trim(dimLabel), varid=varid)
      if (ncerr == NF90_NOERR) return
      ! the variable does not exist in the NetCDF file yet -> add it
      ! access the undistributed dimension attribute on the grid
      call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                             name=trim(dimLabel), itemCount=valueCount, typekind=typekind, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if( typekind == ESMF_TYPEKIND_R4 ) then
        allocate(valueListr4(valueCount))
        call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                               name=trim(dimLabel), valueList=valueListr4, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      else if ( typekind == ESMF_TYPEKIND_R8) then
        allocate(valueListr8(valueCount))
        call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                              name=trim(dimLabel), valueList=valueListr8, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      endif
      ! now add it to the NetCDF file
      ncerr = nf90_redef(ncid=ncid)
      if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      ncerr = nf90_inq_dimid(ncid, trim(dimLabel), dimid=dimid)
      if (ncerr /= NF90_NOERR) then
        ! dimension does not yet exist, and must be defined
        ncerr = nf90_def_dim(ncid, trim(dimLabel), valueCount, dimid=dimid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      endif
      if( typekind == ESMF_TYPEKIND_R4 ) then
        ncerr = nf90_def_var(ncid, trim(dimLabel), NF90_FLOAT, &
                             dimids=(/dimid/), varid=varid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        ncerr = nf90_enddef(ncid=ncid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        ncerr = nf90_put_var(ncid, varid, values=valueListr4)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        deallocate(valueListr4)
      else if(typekind == ESMF_TYPEKIND_R8) then
        ncerr = nf90_def_var(ncid, trim(dimLabel), NF90_DOUBLE, &
                             dimids=(/dimid/), varid=varid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        ncerr = nf90_enddef(ncid=ncid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        ncerr = nf90_put_var(ncid, varid, values=valueListr8)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
        deallocate(valueListr8)
      endif
      ! add attributes to this vertical variable
      call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                             attnestflag=ESMF_ATTNEST_OFF, count=attCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if (attCount>0) then
        ncerr = nf90_redef(ncid=ncid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      endif
      ! loop over all the attributes
      do j=1, attCount
        call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3",       &
                               attnestflag=ESMF_ATTNEST_OFF, attributeIndex=j, &
                               name=attName, typekind=typekind, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        ! test for name starting with trim(dimLabel)":"
        if (index(trim(attName), trim(dimLabel)//":") == 1) then
          ind = len(trim(dimLabel)//":")
          ! found a matching attributes
          if (typekind == ESMF_TYPEKIND_CHARACTER) then
            call ESMF_AttributeGet(grid, &
                                   convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=valueS, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ncerr = nf90_put_att(ncid, varid, &
                                 trim(attName(ind+1:len(attName))), values=valueS)

            if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          else if (typekind == ESMF_TYPEKIND_I4) then
            call ESMF_AttributeGet(grid, &
                                   convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=valueI4, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ncerr = nf90_put_att(ncid, varid, &
                                 trim(attName(ind+1:len(attName))), values=valueI4)

            if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          else if (typekind == ESMF_TYPEKIND_R4) then
            call ESMF_AttributeGet(grid, &
                                   convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=valueR4, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ncerr = nf90_put_att(ncid, varid, &
                                 trim(attName(ind+1:len(attName))), values=valueR4)

            if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          else if (typekind == ESMF_TYPEKIND_R8) then
            call ESMF_AttributeGet(grid, &
                                   convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=valueR8, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ncerr = nf90_put_att(ncid, varid, &
                                 trim(attName(ind+1:len(attName))), values=valueR8)

            if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
          endif
        endif
      enddo
      if (attCount>0) then
        ncerr = nf90_enddef(ncid=ncid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      endif
    end subroutine

  end subroutine ioCompRun

  !-----------------------------------------------------------------------------

  subroutine ESMFproto_FieldMakeSingleTile(field, tile, tileField, petList, rc)
    type(ESMF_Field),     intent(in)              :: field
    integer,              intent(in)              :: tile
    type(ESMF_Field),     intent(out)             :: tileField
    integer, allocatable, intent(inout), optional :: petList(:)
    integer,              intent(out),   optional :: rc
    
    ! Take in a field on a multi-tile grid and return a field that only 
    ! references a single tile.
    
    ! This routine only works with references, no data copies are being 
    ! made. The single tile field that is returned points to the original
    ! field allocation. 
    
    ! The original field passed in remains valid.

    type(ESMF_TypeKind_Flag)                :: typekind
    type(ESMF_Index_Flag)                   :: indexflag
    type(ESMF_Grid)                         :: grid, tileGrid
    type(ESMF_Array)                        :: array
    type(ESMF_DistGrid)                     :: distgrid
    type(ESMF_DELayout)                     :: delayout
    character(40)                           :: fieldName
    integer                                 :: fieldDimCount, gridDimCount
    integer                                 :: undistDims
    integer                                 :: localDeCount, deCount, tileCount
    integer,                   allocatable  :: gridToFieldMap(:)
    integer,                   allocatable  :: ungriddedLBound(:)
    integer,                   allocatable  :: ungriddedUBound(:)
    integer,                   allocatable  :: localDeToDeMap(:), deToTileMap(:)
    type(ESMF_LocalArray),     allocatable  :: lArrayList(:)
    integer                                 :: i
    integer                                 :: tileDeCount, tileLocalDeCount
    integer,                   allocatable  :: distgridToArrayMap(:)
    integer,                   allocatable  :: undistLBound(:), undistUBound(:)
    integer,                   allocatable  :: petMap(:), tilePetMap(:)
    integer,                   allocatable  :: minIndexPDe(:,:)
    integer,                   allocatable  :: maxIndexPDe(:,:)
    integer,                   allocatable  :: minIndexPTile(:,:)
    integer,                   allocatable  :: maxIndexPTile(:,:)
    integer,                   allocatable  :: deBlockList(:,:,:)

    character(800)                          :: msg

    if (present(rc)) rc = ESMF_SUCCESS

    ! access information from the incoming field
    call ESMF_FieldGet(field, array=array, typekind=typekind, &
                       dimCount=fieldDimCount, name=fieldName, grid=grid, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! access information from the associated grid
    call ESMF_GridGet(grid, dimCount=gridDimCount, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
#if 0
    write(msg,*) "fieldDimCount=",fieldDimCount,"gridDimCount=",gridDimCount
    call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)
#endif
    ! access list type information from the incoming field
    allocate(gridToFieldMap(gridDimCount))
    undistDims = fieldDimCount-gridDimCount
    if (undistDims < 0) undistDims = 0  ! this supports replicated dimensions
    allocate(ungriddedLBound(undistDims))
    allocate(ungriddedUBound(undistDims))
    call ESMF_FieldGet(field, gridToFieldMap=gridToFieldMap, &
                       ungriddedLBound=ungriddedLBound, ungriddedUBound=ungriddedUBound, &
                       rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! access information from associated array
    call ESMF_ArrayGet(array, distgrid=distgrid, delayout=delayout, &
                       indexflag=indexflag, localDeCount=localDeCount, deCount=deCount, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! access list type information from associated array
    allocate(localDeToDeMap(localDeCount), deToTileMap(deCount))
    allocate(distgridToArrayMap(gridDimCount))
    allocate(undistLBound(undistDims))
    allocate(undistUBound(undistDims))
    call ESMF_ArrayGet(array, tileCount=tileCount,                             &
                       localDeToDeMap=localDeToDeMap, deToTileMap=deToTileMap, &
                       distgridToArrayMap=distgridToArrayMap,                  &
                       undistLBound=undistLBound, undistUBound=undistUBound, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! access list type information from associated distgrid
    allocate(minIndexPDe(gridDimCount,deCount))
    allocate(maxIndexPDe(gridDimCount,deCount))
    allocate(minIndexPTile(gridDimCount,tileCount))
    allocate(maxIndexPTile(gridDimCount,tileCount))
    call ESMF_DistGridGet(distgrid,                                                 &
                          minIndexPDe=minIndexPDe, maxIndexPDe=maxIndexPDe,         &
                          minIndexPTile=minIndexPTile, maxIndexPTile=maxIndexPTile, &
                          rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! access list type information from associated delayout
    allocate(petMap(deCount))
    call ESMF_DELayoutGet(delayout, petMap=petMap, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! construct data structures selecting specific tile
    allocate(lArrayList(localDeCount))
    tileLocalDeCount = 0
    do i=1, localDeCount
      if (deToTileMap(localDeToDeMap(i)+1) == tile) then
        ! localDe is on tile
        tileLocalDeCount = tileLocalDeCount + 1
        call ESMF_ArrayGet(array, localDe=i-1, &
                           localarray=lArrayList(tileLocalDeCount), rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      endif
    enddo
    allocate(tilePetMap(deCount))
    allocate(deBlockList(gridDimCount,2,deCount))
    tileDeCount=0
    do i=1, deCount
      if (deToTileMap(i) == tile) then
        ! DE is on tile
        tileDeCount = tileDeCount + 1
        tilePetMap(tileDeCount) = petMap(i)
        deBlockList(:,1,tileDeCount) = minIndexPDe(:,i)
        deBlockList(:,2,tileDeCount) = maxIndexPDe(:,i)
      endif
    enddo
    if (present(petList)) then
      if (.not.allocated(petList)) then
        allocate(petList(tileDeCount))
      else if (size(petList)/=tileDeCount) then
        deallocate(petList)
        allocate(petList(tileDeCount))
      endif
      petList(:) = tilePetMap(1:tileDeCount)
    endif
    ! create DELayout and DistGrid that only contain the single tile
    delayout = ESMF_DELayoutCreate(tilePetMap(1:tileDeCount), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    distgrid = ESMF_DistGridCreate(minIndex=minIndexPTile(:,tile),                    &
                                   maxIndex=maxIndexPTile(:,tile), delayout=delayout, &
                                   deBlockList=deBlockList(:,:,1:tileDeCount), rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! create an Array that only holds tile specific allocations
    if (tileLocalDeCount>0) then
      array = ESMF_ArrayCreate(distgrid, lArrayList(1:tileLocalDeCount),                         &
                               indexflag=indexflag,                                              &
                               distgridToArrayMap=distgridToArrayMap, undistLBound=undistLBound, &
                               undistUBound=undistUBound, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    else
      array = ESMF_ArrayCreate(distgrid, typekind, indexflag=indexflag,                          &
                               distgridToArrayMap=distgridToArrayMap, undistLBound=undistLBound, &
                               undistUBound=undistUBound, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif

    ! create a grid on the new distgrid
    tileGrid = ESMF_GridCreate(distgrid, indexflag=indexflag, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! alias the Attributes on grid level
    call ESMF_AttributeCopy(grid, tileGrid, attcopy=ESMF_ATTCOPY_REFERENCE, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! create the tile specific field from the array
    tileField = ESMF_FieldCreate(tileGrid, array=array, name=fieldName, &
                                 gridToFieldMap=gridToFieldMap, ungriddedLBound=ungriddedLBound, &
                                 ungriddedUBound=ungriddedUBound, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! alias the Attributes on field level
    call ESMF_AttributeCopy(field, tileField, attcopy=ESMF_ATTCOPY_REFERENCE, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! local garbage collection
    deallocate(localDeToDeMap, deToTileMap)
    deallocate(petMap,         tilePetMap)
    deallocate(minIndexPDe,    maxIndexPDe)
    deallocate(minIndexPTile,  maxIndexPTile)
    deallocate(gridToFieldMap, ungriddedLBound, ungriddedUBound)
    deallocate(deBlockList)

  end subroutine ESMFproto_FieldMakeSingleTile

!
!-----------------------------------------------------------------------
  subroutine splat4(idrt,jmax,aslat)

      implicit none
      integer,intent(in)  :: idrt,jmax
      real(4),intent(out) :: aslat(jmax)
!
      integer,parameter   :: KD=SELECTED_REAL_KIND(15,45)
      real(kind=KD)       :: pk(jmax/2),pkm1(jmax/2),pkm2(jmax/2)
      real(kind=KD)       :: aslatd(jmax/2),sp,spmax,eps=10.d0*epsilon(sp)
      integer,PARAMETER   :: JZ=50
      real(8) bz(jz)
      data bz        / 2.4048255577d0,  5.5200781103d0, &
       8.6537279129d0, 11.7915344391d0, 14.9309177086d0, 18.0710639679d0, &
      21.2116366299d0, 24.3524715308d0, 27.4934791320d0, 30.6346064684d0, &
      33.7758202136d0, 36.9170983537d0, 40.0584257646d0, 43.1997917132d0, &
      46.3411883717d0, 49.4826098974d0, 52.6240518411d0, 55.7655107550d0, &
      58.9069839261d0, 62.0484691902d0, 65.1899648002d0, 68.3314693299d0, &
      71.4729816036d0, 74.6145006437d0, 77.7560256304d0, 80.8975558711d0, &
      84.0390907769d0, 87.1806298436d0, 90.3221726372d0, 93.4637187819d0, &
      96.6052679510d0, 99.7468198587d0, 102.888374254d0, 106.029930916d0, &
      109.171489649d0, 112.313050280d0, 115.454612653d0, 118.596176630d0, &
      121.737742088d0, 124.879308913d0, 128.020877005d0, 131.162446275d0, &
      134.304016638d0, 137.445588020d0, 140.587160352d0, 143.728733573d0, &
      146.870307625d0, 150.011882457d0, 153.153458019d0, 156.295034268d0 /
      real(8)           :: dlt,d1=1.d0
      integer           :: jhe,jho,j0=0
!     real(8),parameter :: PI=3.14159265358979d0,C=(1.d0-(2.d0/PI)**2)*0.25d0
      real(8),parameter ::                       C=(1.d0-(2.d0/PI)**2)*0.25d0
      real(8) r
      integer jh,js,n,j
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  GAUSSIAN LATITUDES
      IF(IDRT.EQ.4) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        R=1.d0/SQRT((JMAX+0.5d0)**2+C)
        DO J=1,MIN(JH,JZ)
          ASLATD(J)=COS(BZ(J)*R)
        ENDDO
        DO J=JZ+1,JH
          ASLATD(J)=COS((BZ(JZ)+(J-JZ)*PI)*R)
        ENDDO
        SPMAX=1.d0
        DO WHILE(SPMAX.GT.EPS)
          SPMAX=0.d0
          DO J=1,JH
            PKM1(J)=1.d0
            PK(J)=ASLATD(J)
          ENDDO
          DO N=2,JMAX
            DO J=1,JH
              PKM2(J)=PKM1(J)
              PKM1(J)=PK(J)
              PK(J)=((2*N-1)*ASLATD(J)*PKM1(J)-(N-1)*PKM2(J))/N
            ENDDO
          ENDDO
          DO J=1,JH
            SP=PK(J)*(1.d0-ASLATD(J)**2)/(JMAX*(PKM1(J)-ASLATD(J)*PK(J)))
            ASLATD(J)=ASLATD(J)-SP
            SPMAX=MAX(SPMAX,ABS(SP))
          ENDDO
        ENDDO
!
        DO J=1,JH
          ASLAT(J)=ASLATD(J)
          ASLAT(JMAX+1-J)=-ASLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          ASLAT(JHE)=0.d0
        ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  EQUALLY-SPACED LATITUDES INCLUDING POLES
      ELSEIF(IDRT.EQ.0) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        JHO=JHE-1
        DLT=PI/(JMAX-1)
        ASLAT(1)=1.d0
        DO J=2,JH
          ASLAT(J)=COS((J-1)*DLT)
        ENDDO
!
        DO J=1,JH
          ASLAT(JMAX+1-J)=-ASLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          ASLAT(JHE)=0.d0
        ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  EQUALLY-SPACED LATITUDES EXCLUDING POLES
      ELSEIF(IDRT.EQ.256) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        JHO=JHE
        DLT=PI/JMAX
        ASLAT(1)=1.d0
        DO J=1,JH
          ASLAT(J)=COS((J-0.5)*DLT)
        ENDDO

        DO J=1,JH
          ASLAT(JMAX+1-J)=-ASLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          ASLAT(JHE)=0.d0
        ENDIF
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     end subroutine splat4
!----------------------------------------------------------------------
     subroutine splat8(idrt,jmax,aslat)
!$$$
      implicit none
      integer,intent(in)  :: idrt,jmax
      real(8),intent(out) :: aslat(jmax)
!
      integer,parameter   :: KD=SELECTED_REAL_KIND(15,45)
      real(kind=KD)       :: pk(jmax/2),pkm1(jmax/2),pkm2(jmax/2)
      real(kind=KD)       :: aslatd(jmax/2),sp,spmax,eps=10.d0*epsilon(sp)
      integer,parameter   :: jz=50
      real(8) bz(jz)
      data bz        / 2.4048255577d0,  5.5200781103d0, &
       8.6537279129d0, 11.7915344391d0, 14.9309177086d0, 18.0710639679d0, &
      21.2116366299d0, 24.3524715308d0, 27.4934791320d0, 30.6346064684d0, &
      33.7758202136d0, 36.9170983537d0, 40.0584257646d0, 43.1997917132d0, &
      46.3411883717d0, 49.4826098974d0, 52.6240518411d0, 55.7655107550d0, &
      58.9069839261d0, 62.0484691902d0, 65.1899648002d0, 68.3314693299d0, &
      71.4729816036d0, 74.6145006437d0, 77.7560256304d0, 80.8975558711d0, &
      84.0390907769d0, 87.1806298436d0, 90.3221726372d0, 93.4637187819d0, &
      96.6052679510d0, 99.7468198587d0, 102.888374254d0, 106.029930916d0, &
      109.171489649d0, 112.313050280d0, 115.454612653d0, 118.596176630d0, &
      121.737742088d0, 124.879308913d0, 128.020877005d0, 131.162446275d0, &
      134.304016638d0, 137.445588020d0, 140.587160352d0, 143.728733573d0, &
      146.870307625d0, 150.011882457d0, 153.153458019d0, 156.295034268d0 /
      real(8)           :: dlt,d1=1.d0
      integer(4)        :: jhe,jho,j0=0
!     real(8),parameter :: PI=3.14159265358979d0,C=(1.d0-(2.d0/PI)**2)*0.25d0
      real(8),parameter ::                       C=(1.d0-(2.d0/PI)**2)*0.25d0
      real(8) r
      integer jh,js,n,j
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  GAUSSIAN LATITUDES
      IF(IDRT.EQ.4) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        R=1.d0/SQRT((JMAX+0.5d0)**2+C)
        DO J=1,MIN(JH,JZ)
          ASLATD(J)=COS(BZ(J)*R)
        ENDDO
        DO J=JZ+1,JH
          ASLATD(J)=COS((BZ(JZ)+(J-JZ)*PI)*R)
        ENDDO
        SPMAX=1.d0
        DO WHILE(SPMAX.GT.EPS)
          SPMAX=0.d0
          DO J=1,JH
            PKM1(J)=1.d0
            PK(J)=ASLATD(J)
          ENDDO
          DO N=2,JMAX
            DO J=1,JH
              PKM2(J)=PKM1(J)
              PKM1(J)=PK(J)
              PK(J)=((2*N-1)*ASLATD(J)*PKM1(J)-(N-1)*PKM2(J))/N
            ENDDO
          ENDDO
          DO J=1,JH
            SP=PK(J)*(1.d0-ASLATD(J)**2)/(JMAX*(PKM1(J)-ASLATD(J)*PK(J)))
            ASLATD(J)=ASLATD(J)-SP
            SPMAX=MAX(SPMAX,ABS(SP))
          ENDDO
        ENDDO
!
        DO J=1,JH
          ASLAT(J)=ASLATD(J)
          ASLAT(JMAX+1-J)=-ASLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          ASLAT(JHE)=0.d0
        ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  EQUALLY-SPACED LATITUDES INCLUDING POLES
      ELSEIF(IDRT.EQ.0) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        JHO=JHE-1
        DLT=PI/(JMAX-1)
        ASLAT(1)=1.d0
        DO J=2,JH
          ASLAT(J)=COS((J-1)*DLT)
        ENDDO
        DO J=1,JH
          ASLAT(JMAX+1-J)=-ASLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          ASLAT(JHE)=0.d0
        ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  EQUALLY-SPACED LATITUDES EXCLUDING POLES
      ELSEIF(IDRT.EQ.256) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        JHO=JHE
        DLT=PI/JMAX
        ASLAT(1)=1.d0
        DO J=1,JH
          ASLAT(J)=COS((J-0.5d0)*DLT)
        ENDDO
!
        DO J=1,JH
          ASLAT(JMAX+1-J)=-ASLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          ASLAT(JHE)=0.d0
        ENDIF
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     end subroutine splat8
!
!
   subroutine rtll(tlmd,tphd,almd,aphd,tlm0d,tph0d)
!-------------------------------------------------------------------------------
      real(ESMF_KIND_R8), intent(in) :: tlmd, tphd
      real(ESMF_KIND_R8), intent(out) :: almd, aphd
      real(ESMF_KIND_R8), intent(in) :: tph0d, tlm0d
!-------------------------------------------------------------------------------
!     real(ESMF_KIND_R8), parameter :: pi=3.14159265358979323846
      real(ESMF_KIND_R8), parameter :: dtr=pi/180.0
!
      real(ESMF_KIND_R8) :: tph0, ctph0, stph0, tlm, tph, stph, ctph, ctlm, stlm, aph, cph
      real(ESMF_KIND_R8) :: xx, yy
!-------------------------------------------------------------------------------
!
      tph0=tph0d*dtr
      ctph0=cos(tph0)
      stph0=sin(tph0)
!
      tlm=tlmd*dtr
      tph=tphd*dtr
      stph=sin(tph)
      ctph=cos(tph)
      ctlm=cos(tlm)
      stlm=sin(tlm)
!
      xx=stph0*ctph*ctlm+ctph0*stph
      xx=max(xx,-1.0)
      xx=min(xx, 1.0)
      aph=asin(xx)
      cph=cos(aph)
!
      xx=(ctph0*ctph*ctlm-stph0*stph)/cph
      xx=max(xx,-1.0)
      xx=min(xx, 1.0)
      xx=acos(xx)/dtr
      yy=ctph*stlm/cph
      xx=sign(xx,yy)
      almd=tlm0d+xx

      aphd=aph/dtr
!
      if (almd > 180.0) then
         almd=almd-360.0
      end if
      if (almd < -180.0)  then
         almd=almd+360.0
      end if
!
      return
!
     end subroutine rtll
!
!-----------------------------------------------------------------------
!
     subroutine lambert(stlat1,stlat2,c_lat,c_lon,glon,glat,x,y,inv)

!-------------------------------------------------------------------------------
      real(ESMF_KIND_R8),      intent(in)  :: stlat1,stlat2,c_lat,c_lon
      real(ESMF_KIND_R8),      intent(inout)  :: glon, glat
      real(ESMF_KIND_R8),      intent(inout) :: x, y
      integer,                 intent(in)  :: inv
!-------------------------------------------------------------------------------
!     real(ESMF_KIND_R8), parameter :: pi=3.14159265358979323846
      real(ESMF_KIND_R8), parameter :: dtor=pi/180.0
      real(ESMF_KIND_R8), parameter :: rtod=180.0/pi
      real(ESMF_KIND_R8), parameter :: a = 6371200.0
!-------------------------------------------------------------------------------
! inv == 1     (glon,glat) ---> (x,y)    lat/lon to grid
! inv == -1    (x,y) ---> (glon,glat)    grid to lat/lon

      real(ESMF_KIND_R8) :: en,f,rho,rho0, dlon, theta, xp, yp

      IF (stlat1 == stlat2) THEN
         en=sin(stlat1*dtor)
      ELSE
         en=log(cos(stlat1*dtor)/cos(stlat2*dtor))/ &
            log(tan((45+0.5*stlat2)*dtor)/tan((45+0.5*stlat1)*dtor))
      ENDIF

      f=(cos(stlat1*dtor)*tan((45+0.5*stlat1)*dtor)**en)/en
      rho0=a*f/(tan((45+0.5*c_lat)*dtor)**en)

      if (inv == 1) then          ! FORWARD TRANSFORMATION
            rho=a*f/(tan((45+0.5*glat)*dtor)**en)
            dlon=modulo(glon-c_lon+180+3600,360.)-180.D0
            theta=en*dlon*dtor
            x=rho*sin(theta)
            y=rho0-rho*cos(theta)
      else if (inv == -1) then    ! INVERSE TRANSFORMATION
            y=rho0-y
            rho = sqrt(x*x+y*y)
            theta=atan2(x,y)
            glon=c_lon+(theta/en)*rtod
            glon=modulo(glon+180+3600,360.)-180.D0
!            glat=(2.0*atan((a*f/rho)**(1.0/en))-0.5*pi)*rtod
            glat=(0.5*pi-2.0*atan((rho/(a*f))**(1.0/en)))*rtod
      else
        write (unit=*,fmt=*) " lambert: unknown inv argument"
        return
      end if

      return
     end subroutine lambert
!
!-----------------------------------------------------------------------
!
     subroutine get_outfile(nfl, filename, outfile_name,noutfile)
       integer, intent(in)          :: nfl
       character(*), intent(in)     :: filename(:,:)
       character(*), intent(inout)  :: outfile_name(:)
       integer, intent(inout)       :: noutfile

       integer        :: i,j,n,idx
       logical        :: found
!
       noutfile = 0
       do i=1,nfl

         loopj: do j=1, 2000

           if( trim(filename(j,i)) == '') exit loopj
           if( trim(filename(j,i)) == 'none') cycle

           found = .false.
           loopn: do n=1, noutfile
             if(trim(filename(j,i)) == trim(outfile_name(n))) then
               found = .true.
               exit loopn
             endif
           enddo loopn

           if (.not.found) then
             noutfile = noutfile + 1
             outfile_name(noutfile) = trim(filename(j,i))
!             print *,'in get outfile,noutfile=', noutfile,' outfile=',trim(filename(j,i))
           endif
         enddo loopj
!
       enddo

     end subroutine get_outfile
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
    end module  module_wrt_grid_comp
!
!-----------------------------------------------------------------------
