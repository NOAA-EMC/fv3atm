#define ESMF_ERR_RETURN(rc) \
    if (ESMF_LogFoundError(rc, msg="Breaking out of subroutine", line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

#define NC_ERR_STOP(status) \
    if (status /= nf90_noerr) write(0,*) "file: ", __FILE__, " line: ", __LINE__, trim(nf90_strerror(status)); \
    if (status /= nf90_noerr) call ESMF_Finalize(endflag=ESMF_END_ABORT)

module module_write_netcdf

  use esmf
  use netcdf
  use module_fv3_io_def,only : ideflate, nbits, &
                               ichunk2d,jchunk2d,ichunk3d,jchunk3d,kchunk3d, &
                               output_grid,dx,dy,lon1,lat1,lon2,lat2
  use mpi

  implicit none
  private
  public write_netcdf

  logical :: par

  interface quantize_array
     module procedure quantize_array_3d
     module procedure quantize_array_4d
  end interface

  contains

!----------------------------------------------------------------------------------------
  subroutine write_netcdf(wrtfb, filename, &
                          use_parallel_netcdf, mpi_comm, mype, &
                          grid_id, rc)
!
    type(ESMF_FieldBundle), intent(in) :: wrtfb
    character(*), intent(in)           :: filename
    logical, intent(in)                :: use_parallel_netcdf
    integer, intent(in)                :: mpi_comm
    integer, intent(in)                :: mype
    integer, intent(in)                :: grid_id
    integer, optional,intent(out)      :: rc
!
!** local vars
    integer :: i,j,t, istart,iend,jstart,jend
    integer :: im, jm, lm

    integer, dimension(:), allocatable              :: fldlev

    real(ESMF_KIND_R4), dimension(:,:), pointer     :: array_r4
    real(ESMF_KIND_R4), dimension(:,:,:), pointer   :: array_r4_cube
    real(ESMF_KIND_R4), dimension(:,:,:), pointer   :: array_r4_3d
    real(ESMF_KIND_R4), dimension(:,:,:,:), pointer :: array_r4_3d_cube

    real(ESMF_KIND_R8), dimension(:,:), pointer     :: array_r8
    real(ESMF_KIND_R8), dimension(:,:,:), pointer   :: array_r8_cube
    real(ESMF_KIND_R8), dimension(:,:,:), pointer   :: array_r8_3d
    real(ESMF_KIND_R8), dimension(:,:,:,:), pointer :: array_r8_3d_cube

    real(8), dimension(:), allocatable :: x,y
    integer :: fieldCount, fieldDimCount, gridDimCount
    integer, dimension(:), allocatable   :: ungriddedLBound, ungriddedUBound
    integer, dimension(:), allocatable   :: start_idx

    type(ESMF_Field), allocatable        :: fcstField(:)
    type(ESMF_TypeKind_Flag)             :: typekind
    type(ESMF_TypeKind_Flag)             :: attTypeKind
    type(ESMF_Grid)                      :: wrtgrid
    type(ESMF_Array)                     :: array
    type(ESMF_DistGrid)                  :: distgrid

    integer :: attCount
    character(len=ESMF_MAXSTR) :: attName, fldName

    integer :: varival
    real(4) :: varr4val, dataMin, dataMax
    real(4), allocatable, dimension(:) :: compress_err
    real(8) :: varr8val
    character(len=ESMF_MAXSTR) :: varcval

    integer :: ncerr,ierr
    integer :: ncid
    integer :: oldMode
    integer :: im_dimid, jm_dimid, tile_dimid, pfull_dimid, phalf_dimid, time_dimid
    integer :: im_varid, jm_varid, tile_varid, lon_varid, lat_varid
    integer, dimension(:), allocatable :: dimids_2d, dimids_3d
    integer, dimension(:), allocatable :: varids
    logical shuffle

    logical :: is_cubed_sphere
    integer :: rank, deCount, localDeCount, dimCount, tileCount
    integer :: my_tile, start_i, start_j
    integer, dimension(:,:), allocatable :: minIndexPDe, maxIndexPDe
    integer, dimension(:,:), allocatable :: minIndexPTile, maxIndexPTile
    integer, dimension(:), allocatable :: deToTileMap, localDeToDeMap
    logical :: do_io
    integer :: par_access
!
    is_cubed_sphere = .false.
    tileCount = 0
    my_tile = 0
    start_i = -10000000
    start_j = -10000000

    par = use_parallel_netcdf
    do_io = par .or. (mype==0)

    call ESMF_FieldBundleGet(wrtfb, fieldCount=fieldCount, rc=rc); ESMF_ERR_RETURN(rc)

    allocate(compress_err(fieldCount)); compress_err=-999.
    allocate(fldlev(fieldCount)) ; fldlev = 0
    allocate(fcstField(fieldCount))
    allocate(varids(fieldCount))

    call ESMF_FieldBundleGet(wrtfb, fieldList=fcstField, grid=wrtGrid, &
!                             itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                             rc=rc); ESMF_ERR_RETURN(rc)

    call ESMF_GridGet(wrtgrid, dimCount=gridDimCount, rc=rc); ESMF_ERR_RETURN(rc)

    do i=1,fieldCount
       call ESMF_FieldGet(fcstField(i), dimCount=fieldDimCount, array=array, rc=rc); ESMF_ERR_RETURN(rc)

       if (fieldDimCount > 3) then
          write(0,*)"write_netcdf: Only 2D and 3D fields are supported!"
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if

       ! use first field to determine tile number, grid size, start index etc.
       if (i == 1) then
          call ESMF_ArrayGet(array, &
                             distgrid=distgrid, &
                             dimCount=dimCount, &
                             deCount=deCount, &
                             localDeCount=localDeCount, &
                             tileCount=tileCount, &
                             rc=rc); ESMF_ERR_RETURN(rc)

          allocate(minIndexPDe(dimCount,deCount))
          allocate(maxIndexPDe(dimCount,deCount))
          allocate(minIndexPTile(dimCount, tileCount))
          allocate(maxIndexPTile(dimCount, tileCount))
          call ESMF_DistGridGet(distgrid, &
                                minIndexPDe=minIndexPDe, maxIndexPDe=maxIndexPDe, &
                                minIndexPTile=minIndexPTile, maxIndexPTile=maxIndexPTile, &
                                rc=rc); ESMF_ERR_RETURN(rc)

          allocate(deToTileMap(deCount))
          allocate(localDeToDeMap(localDeCount))
          call ESMF_ArrayGet(array, &
                             deToTileMap=deToTileMap, &
                             localDeToDeMap=localDeToDeMap, &
                             rc=rc); ESMF_ERR_RETURN(rc)

          is_cubed_sphere = (tileCount == 6)
          my_tile = deToTileMap(localDeToDeMap(1)+1)
          im = maxIndexPTile(1,1)
          jm = maxIndexPTile(2,1)
          start_i = minIndexPDe(1,localDeToDeMap(1)+1)
          start_j = minIndexPDe(2,localDeToDeMap(1)+1)
          if (.not. par) then
             start_i = 1
             start_j = 1
          end if
       end if

       if (fieldDimCount > gridDimCount) then
         allocate(ungriddedLBound(fieldDimCount-gridDimCount))
         allocate(ungriddedUBound(fieldDimCount-gridDimCount))
         call ESMF_FieldGet(fcstField(i), &
                            ungriddedLBound=ungriddedLBound, &
                            ungriddedUBound=ungriddedUBound, rc=rc); ESMF_ERR_RETURN(rc)
         fldlev(i) = ungriddedUBound(fieldDimCount-gridDimCount) - &
                     ungriddedLBound(fieldDimCount-gridDimCount) + 1
         deallocate(ungriddedLBound)
         deallocate(ungriddedUBound)
       else if (fieldDimCount == 2) then
         fldlev(i) = 1
       end if
    end do

    lm = maxval(fldlev(:))

    ! for serial output allocate 'global' arrays
    if (.not. par) then
       allocate(array_r4(im,jm))
       allocate(array_r8(im,jm))
       allocate(array_r4_3d(im,jm,lm))
       allocate(array_r8_3d(im,jm,lm))
       if (is_cubed_sphere) then
          allocate(array_r4_cube(im,jm,tileCount))
          allocate(array_r8_cube(im,jm,tileCount))
          allocate(array_r4_3d_cube(im,jm,lm,tileCount))
          allocate(array_r8_3d_cube(im,jm,lm,tileCount))
       end if
    end if

    ! create netcdf file and enter define mode
    if (do_io) then

       if (par) then
          ncerr = nf90_create(trim(filename),&
                  cmode=IOR(IOR(NF90_CLOBBER,NF90_NETCDF4),NF90_CLASSIC_MODEL),&
                  comm=mpi_comm, info = MPI_INFO_NULL, ncid=ncid); NC_ERR_STOP(ncerr)
       else
          ncerr = nf90_create(trim(filename),&
                  cmode=IOR(IOR(NF90_CLOBBER,NF90_NETCDF4),NF90_CLASSIC_MODEL),&
                  ncid=ncid); NC_ERR_STOP(ncerr)
       end if

       ! disable auto filling.
       ncerr = nf90_set_fill(ncid, NF90_NOFILL, oldMode); NC_ERR_STOP(ncerr)

       ! define dimensions [grid_xt, grid_yta ,(pfull/phalf), (tile), time]
       ncerr = nf90_def_dim(ncid, "grid_xt", im, im_dimid); NC_ERR_STOP(ncerr)
       ncerr = nf90_def_dim(ncid, "grid_yt", jm, jm_dimid); NC_ERR_STOP(ncerr)
       if (lm > 1) then
         call add_dim(ncid, "pfull", pfull_dimid, wrtgrid, rc)
         call add_dim(ncid, "phalf", phalf_dimid, wrtgrid, rc)
       end if
       if (is_cubed_sphere) then
          ncerr = nf90_def_dim(ncid, "tile", tileCount, tile_dimid); NC_ERR_STOP(ncerr)
       end if
       call add_dim(ncid, "time", time_dimid, wrtgrid, rc)

       ! define coordinate variables
       ncerr = nf90_def_var(ncid, "grid_xt", NF90_DOUBLE, im_dimid, im_varid); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, im_varid, "cartesian_axis", "X"); NC_ERR_STOP(ncerr)
       ncerr = nf90_def_var(ncid, "grid_yt", NF90_DOUBLE, jm_dimid, jm_varid); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, jm_varid, "cartesian_axis", "Y"); NC_ERR_STOP(ncerr)
       if (is_cubed_sphere) then
          ncerr = nf90_def_var(ncid, "tile", NF90_INT, tile_dimid, tile_varid); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, tile_varid, "long_name", "cubed-spehere face"); NC_ERR_STOP(ncerr)
       end if

       ! coordinate variable attributes based on output_grid type
       if (trim(output_grid(grid_id)) == 'gaussian_grid' .or. &
           trim(output_grid(grid_id)) == 'global_latlon' .or. &
           trim(output_grid(grid_id)) == 'regional_latlon') then
          ncerr = nf90_put_att(ncid, im_varid, "long_name", "T-cell longitude"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, im_varid, "units", "degrees_E"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, jm_varid, "long_name", "T-cell latiitude"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, jm_varid, "units", "degrees_N"); NC_ERR_STOP(ncerr)
       else if (trim(output_grid(grid_id)) == 'rotated_latlon') then
          ncerr = nf90_put_att(ncid, im_varid, "long_name", "rotated T-cell longiitude"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, im_varid, "units", "degrees"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, jm_varid, "long_name", "rotated T-cell latiitude"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, jm_varid, "units", "degrees"); NC_ERR_STOP(ncerr)
       else if (trim(output_grid(grid_id)) == 'lambert_conformal') then
          ncerr = nf90_put_att(ncid, im_varid, "long_name", "x-coordinate of projection"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, im_varid, "units", "meters"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, jm_varid, "long_name", "y-coordinate of projection"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, jm_varid, "units", "meters"); NC_ERR_STOP(ncerr)
       end if

       ! define longitude variable
       if (is_cubed_sphere) then
          ncerr = nf90_def_var(ncid, "lon", NF90_DOUBLE, [im_dimid,jm_dimid,tile_dimid], lon_varid); NC_ERR_STOP(ncerr)
       else
          ncerr = nf90_def_var(ncid, "lon", NF90_DOUBLE, [im_dimid,jm_dimid           ], lon_varid); NC_ERR_STOP(ncerr)
       end if
       ncerr = nf90_put_att(ncid, lon_varid, "long_name", "T-cell longitude"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, lon_varid, "units", "degrees_E"); NC_ERR_STOP(ncerr)

       ! define latitude variable
       if (is_cubed_sphere) then
          ncerr = nf90_def_var(ncid, "lat", NF90_DOUBLE, [im_dimid,jm_dimid,tile_dimid], lat_varid); NC_ERR_STOP(ncerr)
       else
          ncerr = nf90_def_var(ncid, "lat", NF90_DOUBLE, [im_dimid,jm_dimid           ], lat_varid); NC_ERR_STOP(ncerr)
       end if
       ncerr = nf90_put_att(ncid, lat_varid, "long_name", "T-cell latitude"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, lat_varid, "units", "degrees_N"); NC_ERR_STOP(ncerr)

       if (par) then
          ncerr = nf90_var_par_access(ncid, im_varid, NF90_INDEPENDENT); NC_ERR_STOP(ncerr)
          ncerr = nf90_var_par_access(ncid, lon_varid, NF90_INDEPENDENT); NC_ERR_STOP(ncerr)
          ncerr = nf90_var_par_access(ncid, jm_varid, NF90_INDEPENDENT); NC_ERR_STOP(ncerr)
          ncerr = nf90_var_par_access(ncid, lat_varid, NF90_INDEPENDENT); NC_ERR_STOP(ncerr)
          if (is_cubed_sphere) then
             ncerr = nf90_var_par_access(ncid, tile_varid, NF90_INDEPENDENT); NC_ERR_STOP(ncerr)
          end if
       end if


       call get_global_attr(wrtfb, ncid, rc)


       ! define variables (fields)
       if (is_cubed_sphere) then
          allocate(dimids_2d(4))
          allocate(dimids_3d(5))
          dimids_2d =             [im_dimid,jm_dimid,            tile_dimid,time_dimid]
          if (lm > 1) dimids_3d = [im_dimid,jm_dimid,pfull_dimid,tile_dimid,time_dimid]
       else
          allocate(dimids_2d(3))
          allocate(dimids_3d(4))
          dimids_2d =             [im_dimid,jm_dimid,                       time_dimid]
          if (lm > 1) dimids_3d = [im_dimid,jm_dimid,pfull_dimid,           time_dimid]
       end if

       do i=1, fieldCount
         call ESMF_FieldGet(fcstField(i), name=fldName, rank=rank, typekind=typekind, rc=rc); ESMF_ERR_RETURN(rc)

         par_access = NF90_INDEPENDENT
         ! define variables
         if (rank == 2) then
           if (typekind == ESMF_TYPEKIND_R4) then
             if (ideflate(grid_id) > 0) then
               if (ichunk2d(grid_id) < 0 .or. jchunk2d(grid_id) < 0) then
                  ! let netcdf lib choose chunksize
                  ! shuffle filter on for 2d fields (lossless compression)
                  ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
                          dimids_2d, varids(i), &
                          shuffle=.true.,deflate_level=ideflate(grid_id)); NC_ERR_STOP(ncerr)
               else
                  if (is_cubed_sphere) then
                  ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
                                       dimids_2d, varids(i), &
                                       shuffle=.true.,deflate_level=ideflate(grid_id),&
                                       chunksizes=[ichunk2d(grid_id),jchunk2d(grid_id),tileCount,1]); NC_ERR_STOP(ncerr)
                  else
                  ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
                                       dimids_2d, varids(i), &
                                       shuffle=.true.,deflate_level=ideflate(grid_id),&
                                       chunksizes=[ichunk2d(grid_id),jchunk2d(grid_id),          1]); NC_ERR_STOP(ncerr)
                  end if
               end if
               ! compression filters require collective access.
               par_access = NF90_COLLECTIVE
             else
               ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
                                    dimids_2d, varids(i)); NC_ERR_STOP(ncerr)
             end if
           else if (typekind == ESMF_TYPEKIND_R8) then
             ncerr = nf90_def_var(ncid, trim(fldName), NF90_DOUBLE, &
                                  dimids_2d, varids(i)); NC_ERR_STOP(ncerr)
           else
             write(0,*)'Unsupported typekind ', typekind
             call ESMF_Finalize(endflag=ESMF_END_ABORT)
           end if
         else if (rank == 3) then
           if (typekind == ESMF_TYPEKIND_R4) then
             if (ideflate(grid_id) > 0) then
               ! shuffle filter off for 3d fields using lossy compression
               if (nbits(grid_id) > 0) then
                   shuffle=.false.
               else
                   shuffle=.true.
               end if
               if (ichunk3d(grid_id) < 0 .or. jchunk3d(grid_id) < 0 .or. kchunk3d(grid_id) < 0) then
                  ! let netcdf lib choose chunksize
                  ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
                                       dimids_3d, varids(i), &
                                       shuffle=shuffle,deflate_level=ideflate(grid_id)); NC_ERR_STOP(ncerr)
               else
                  if (is_cubed_sphere) then
                  ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
                                       dimids_3d, varids(i), &
                                       shuffle=shuffle,deflate_level=ideflate(grid_id),&
                                       chunksizes=[ichunk3d(grid_id),jchunk3d(grid_id),kchunk3d(grid_id),tileCount,1]); NC_ERR_STOP(ncerr)
                  else
                  ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
                                       dimids_3d, varids(i), &
                                       shuffle=shuffle,deflate_level=ideflate(grid_id),&
                                       chunksizes=[ichunk3d(grid_id),jchunk3d(grid_id),kchunk3d(grid_id),          1]); NC_ERR_STOP(ncerr)
                  end if
               end if
               ! compression filters require collective access.
               par_access = NF90_COLLECTIVE
             else
               ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
                       dimids_3d, varids(i)); NC_ERR_STOP(ncerr)
             end if
           else if (typekind == ESMF_TYPEKIND_R8) then
             ncerr = nf90_def_var(ncid, trim(fldName), NF90_DOUBLE, &
                                  dimids_3d, varids(i)); NC_ERR_STOP(ncerr)
           else
             write(0,*)'Unsupported typekind ', typekind
             call ESMF_Finalize(endflag=ESMF_END_ABORT)
           end if
         else
           write(0,*)'Unsupported rank ', rank
           call ESMF_Finalize(endflag=ESMF_END_ABORT)
         end if
         if (par) then
             ncerr = nf90_var_par_access(ncid, varids(i), par_access); NC_ERR_STOP(ncerr)
         end if

         ! define variable attributes
         call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                                attnestflag=ESMF_ATTNEST_OFF, count=attCount, &
                                rc=rc); ESMF_ERR_RETURN(rc)

         do j=1,attCount
           call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                                  attnestflag=ESMF_ATTNEST_OFF, attributeIndex=j, &
                                  name=attName, typekind=attTypeKind, &
                                  rc=rc); ESMF_ERR_RETURN(rc)

           if (index(trim(attName),"ESMF") /= 0) then
              cycle
           end if

           if (attTypeKind==ESMF_TYPEKIND_I4) then
              call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                                     name=trim(attName), value=varival, &
                                     rc=rc); ESMF_ERR_RETURN(rc)
              ncerr = nf90_put_att(ncid, varids(i), trim(attName), varival); NC_ERR_STOP(ncerr)

           else if (attTypeKind==ESMF_TYPEKIND_R4) then
              call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                                     name=trim(attName), value=varr4val, &
                                     rc=rc); ESMF_ERR_RETURN(rc)
              ncerr = nf90_put_att(ncid, varids(i), trim(attName), varr4val); NC_ERR_STOP(ncerr)

           else if (attTypeKind==ESMF_TYPEKIND_R8) then
              call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                                     name=trim(attName), value=varr8val, &
                                     rc=rc); ESMF_ERR_RETURN(rc)
              if (trim(attName) /= '_FillValue') then
                 ! FIXME:  _FillValue must be cast to var type when using NF90_NETCDF4
                 ncerr = nf90_put_att(ncid, varids(i), trim(attName), varr8val); NC_ERR_STOP(ncerr)
              end if

           else if (attTypeKind==ESMF_TYPEKIND_CHARACTER) then
              call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                                     name=trim(attName), value=varcval, &
                                     rc=rc); ESMF_ERR_RETURN(rc)
              ncerr = nf90_put_att(ncid, varids(i), trim(attName), trim(varcval)); NC_ERR_STOP(ncerr)

           end if

         end do ! j=1,attCount

         if (is_cubed_sphere) then
            ncerr = nf90_put_att(ncid, varids(i), 'coordinates', 'lon lat'); NC_ERR_STOP(ncerr)
            ncerr = nf90_put_att(ncid, varids(i), 'grid_mapping', 'cubed_sphere'); NC_ERR_STOP(ncerr)
         end if

       end do   ! i=1,fieldCount

       ncerr = nf90_enddef(ncid); NC_ERR_STOP(ncerr)
    end if
    ! end of define mode

    !
    ! write dimension variables and lon,lat variables
    !
    if (allocated(start_idx)) deallocate(start_idx)
    if (is_cubed_sphere) then
       allocate(start_idx(3))
       start_idx = [start_i, start_j, my_tile]
    else
       allocate(start_idx(2))
       start_idx = [start_i, start_j]
    end if

    ! write lon (lon_varid)
    if (par) then
       call ESMF_GridGetCoord(wrtGrid, coordDim=1, farrayPtr=array_r8, rc=rc); ESMF_ERR_RETURN(rc)
       ncerr = nf90_put_var(ncid, lon_varid, values=array_r8, start=start_idx); NC_ERR_STOP(ncerr)
    else
       call ESMF_GridGetCoord(wrtGrid, coordDim=1, array=array, rc=rc); ESMF_ERR_RETURN(rc)
       if (is_cubed_sphere) then
          do t=1,tileCount
             call ESMF_ArrayGather(array, array_r8_cube(:,:,t), rootPet=0, tile=t, rc=rc); ESMF_ERR_RETURN(rc)
          end do
          if (do_io) then
             ncerr = nf90_put_var(ncid, lon_varid, values=array_r8_cube, start=start_idx); NC_ERR_STOP(ncerr)
          end if
       else
          call ESMF_ArrayGather(array, array_r8, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
          if (do_io) then
             ncerr = nf90_put_var(ncid, lon_varid, values=array_r8, start=start_idx); NC_ERR_STOP(ncerr)
          end if
       endif
    end if

    istart = lbound(array_r8,1); iend   = ubound(array_r8,1)
    jstart = lbound(array_r8,2); jend   = ubound(array_r8,2)

    ! write grid_xt (im_varid)
    if (do_io) then
       allocate (x(im))
       if (trim(output_grid(grid_id)) == 'gaussian_grid' .or. &
           trim(output_grid(grid_id)) == 'global_latlon' .or. &
           trim(output_grid(grid_id)) == 'regional_latlon') then
          ncerr = nf90_put_var(ncid, im_varid, values=array_r8(:,jstart), start=[istart], count=[iend-istart+1]); NC_ERR_STOP(ncerr)
       else if (trim(output_grid(grid_id)) == 'rotated_latlon') then
          do i=1,im
             x(i) = lon1(grid_id) + (lon2(grid_id)-lon1(grid_id))/(im-1) * (i-1)
          end do
          ncerr = nf90_put_var(ncid, im_varid, values=x); NC_ERR_STOP(ncerr)
       else if (trim(output_grid(grid_id)) == 'lambert_conformal') then
          do i=1,im
             x(i) = dx(grid_id) * (i-1)
          end do
          ncerr = nf90_put_var(ncid, im_varid, values=x); NC_ERR_STOP(ncerr)
       else if (trim(output_grid(grid_id)) == 'cubed_sphere_grid') then
          do i=1,im
             x(i) = i
          end do
          ncerr = nf90_put_var(ncid, im_varid, values=x); NC_ERR_STOP(ncerr)
       else
          write(0,*)'unknown output_grid ', trim(output_grid(grid_id))
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
    end if

    ! write lat (lat_varid)
    if (par) then
       call ESMF_GridGetCoord(wrtGrid, coordDim=2, farrayPtr=array_r8, rc=rc); ESMF_ERR_RETURN(rc)
       ncerr = nf90_put_var(ncid, lat_varid, values=array_r8, start=start_idx); NC_ERR_STOP(ncerr)
    else
       call ESMF_GridGetCoord(wrtGrid, coordDim=2, array=array, rc=rc); ESMF_ERR_RETURN(rc)
       if (is_cubed_sphere) then
          do t=1,tileCount
             call ESMF_ArrayGather(array, array_r8_cube(:,:,t), rootPet=0, tile=t, rc=rc); ESMF_ERR_RETURN(rc)
          end do
          if (do_io) then
             ncerr = nf90_put_var(ncid, lat_varid, values=array_r8_cube, start=start_idx); NC_ERR_STOP(ncerr)
          end if
       else
          call ESMF_ArrayGather(array, array_r8, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
          if (do_io) then
             ncerr = nf90_put_var(ncid, lat_varid, values=array_r8, start=start_idx); NC_ERR_STOP(ncerr)
          end if
       endif
    end if

    ! write grid_yt (jm_varid)
    if (do_io) then
       allocate (y(jm))
       if (trim(output_grid(grid_id)) == 'gaussian_grid' .or. &
           trim(output_grid(grid_id)) == 'global_latlon' .or. &
           trim(output_grid(grid_id)) == 'regional_latlon') then
          ncerr = nf90_put_var(ncid, jm_varid, values=array_r8(istart,:), start=[jstart], count=[jend-jstart+1]); NC_ERR_STOP(ncerr)
       else if (trim(output_grid(grid_id)) == 'rotated_latlon') then
          do j=1,jm
             y(j) = lat1(grid_id) + (lat2(grid_id)-lat1(grid_id))/(jm-1) * (j-1)
          end do
          ncerr = nf90_put_var(ncid, jm_varid, values=y); NC_ERR_STOP(ncerr)
       else if (trim(output_grid(grid_id)) == 'lambert_conformal') then
          do j=1,jm
             y(j) = dy(grid_id) * (j-1)
          end do
          ncerr = nf90_put_var(ncid, jm_varid, values=y); NC_ERR_STOP(ncerr)
       else if (trim(output_grid(grid_id)) == 'cubed_sphere_grid') then
          do j=1,jm
             y(j) = j
          end do
          ncerr = nf90_put_var(ncid, jm_varid, values=y); NC_ERR_STOP(ncerr)
       else
          write(0,*)'unknown output_grid ', trim(output_grid(grid_id))
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
    end if

    ! write tile (tile_varid)
    if (do_io .and. is_cubed_sphere) then
       ncerr = nf90_put_var(ncid, tile_varid, values=[1,2,3,4,5,6]); NC_ERR_STOP(ncerr)
    end if

    ! write variables (fields)
    do i=1, fieldCount

       call ESMF_FieldGet(fcstField(i),name=fldName,rank=rank,typekind=typekind, rc=rc); ESMF_ERR_RETURN(rc)

       if (rank == 2) then

         if (allocated(start_idx)) deallocate(start_idx)
         if (is_cubed_sphere) then
            allocate(start_idx(4))
            start_idx = [start_i,start_j,my_tile,1]
         else
            allocate(start_idx(3))
            start_idx = [start_i,start_j,        1]
         end if

         if (typekind == ESMF_TYPEKIND_R4) then
            if (par) then
               call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=array_r4, rc=rc); ESMF_ERR_RETURN(rc)
               ncerr = nf90_put_var(ncid, varids(i), values=array_r4, start=start_idx); NC_ERR_STOP(ncerr)
            else
               if (is_cubed_sphere) then
                  call ESMF_FieldGet(fcstField(i), array=array, rc=rc); ESMF_ERR_RETURN(rc)
                  do t=1,tileCount
                     call ESMF_ArrayGather(array, array_r4_cube(:,:,t), rootPet=0, tile=t, rc=rc); ESMF_ERR_RETURN(rc)
                  end do
                  if (do_io) then
                     ncerr = nf90_put_var(ncid, varids(i), values=array_r4_cube, start=start_idx); NC_ERR_STOP(ncerr)
                  end if
               else
                  call ESMF_FieldGather(fcstField(i), array_r4, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
                  if (do_io) then
                     ncerr = nf90_put_var(ncid, varids(i), values=array_r4, start=start_idx); NC_ERR_STOP(ncerr)
                  end if
               end if
            end if
         else if (typekind == ESMF_TYPEKIND_R8) then
            if (par) then
               call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=array_r8, rc=rc); ESMF_ERR_RETURN(rc)
               ncerr = nf90_put_var(ncid, varids(i), values=array_r8, start=start_idx); NC_ERR_STOP(ncerr)
            else
               if (is_cubed_sphere) then
                  call ESMF_FieldGet(fcstField(i), array=array, rc=rc); ESMF_ERR_RETURN(rc)
                  do t=1,tileCount
                     call ESMF_ArrayGather(array, array_r8_cube(:,:,t), rootPet=0, tile=t, rc=rc); ESMF_ERR_RETURN(rc)
                  end do
                  if (do_io) then
                     ncerr = nf90_put_var(ncid, varids(i), values=array_r8_cube, start=start_idx); NC_ERR_STOP(ncerr)
                  end if
               else
                  call ESMF_FieldGather(fcstField(i), array_r8, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
                  if (do_io) then
                     ncerr = nf90_put_var(ncid, varids(i), values=array_r8, start=start_idx); NC_ERR_STOP(ncerr)
                  end if
               end if
            end if
         end if

      else if (rank == 3) then

         if (allocated(start_idx)) deallocate(start_idx)
         if (is_cubed_sphere) then
            allocate(start_idx(5))
            start_idx = [start_i,start_j,1,my_tile,1]
         else
            allocate(start_idx(4))
            start_idx = [start_i,start_j,1,        1]
         end if

         if (typekind == ESMF_TYPEKIND_R4) then
            if (par) then
               call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=array_r4_3d, rc=rc); ESMF_ERR_RETURN(rc)
               if (ideflate(grid_id) > 0 .and. nbits(grid_id) > 0) then
                  dataMax = maxval(array_r4_3d)
                  dataMin = minval(array_r4_3d)
                  call mpi_allreduce(mpi_in_place,dataMax,1,mpi_real4,mpi_max,mpi_comm,ierr)
                  call mpi_allreduce(mpi_in_place,dataMin,1,mpi_real4,mpi_min,mpi_comm,ierr)
                  call quantize_array(array_r4_3d, dataMin, dataMax, nbits(grid_id), compress_err(i))
                  call mpi_allreduce(mpi_in_place,compress_err(i),1,mpi_real4,mpi_max,mpi_comm,ierr)
               end if
               ncerr = nf90_put_var(ncid, varids(i), values=array_r4_3d, start=start_idx); NC_ERR_STOP(ncerr)
            else
               if (is_cubed_sphere) then
                  call ESMF_FieldGet(fcstField(i), array=array, rc=rc); ESMF_ERR_RETURN(rc)
                  do t=1,tileCount
                     call ESMF_ArrayGather(array, array_r4_3d_cube(:,:,:,t), rootPet=0, tile=t, rc=rc); ESMF_ERR_RETURN(rc)
                  end do
                  if (mype==0) then
                     if (ideflate(grid_id) > 0 .and. nbits(grid_id) > 0) then
                        call quantize_array(array_r4_3d_cube, minval(array_r4_3d_cube), maxval(array_r4_3d_cube), nbits(grid_id), compress_err(i))
                     end if
                     ncerr = nf90_put_var(ncid, varids(i), values=array_r4_3d_cube, start=start_idx); NC_ERR_STOP(ncerr)
                  end if
               else
                  call ESMF_FieldGather(fcstField(i), array_r4_3d, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
                  if (mype==0) then
                     if (ideflate(grid_id) > 0 .and. nbits(grid_id) > 0) then
                        call quantize_array(array_r4_3d, minval(array_r4_3d), maxval(array_r4_3d), nbits(grid_id), compress_err(i))
                     end if
                     ncerr = nf90_put_var(ncid, varids(i), values=array_r4_3d, start=start_idx); NC_ERR_STOP(ncerr)
                  end if
               end if
            end if
         else if (typekind == ESMF_TYPEKIND_R8) then
            if (par) then
               call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=array_r8_3d, rc=rc); ESMF_ERR_RETURN(rc)
               ncerr = nf90_put_var(ncid, varids(i), values=array_r8_3d, start=start_idx); NC_ERR_STOP(ncerr)
            else
               if (is_cubed_sphere) then
                  call ESMF_FieldGet(fcstField(i), array=array, rc=rc); ESMF_ERR_RETURN(rc)
                  do t=1,tileCount
                     call ESMF_ArrayGather(array, array_r8_3d_cube(:,:,:,t), rootPet=0, tile=t, rc=rc); ESMF_ERR_RETURN(rc)
                  end do
                  if (mype==0) then
                     ncerr = nf90_put_var(ncid, varids(i), values=array_r8_3d_cube, start=start_idx); NC_ERR_STOP(ncerr)
                  end if
               else
                  call ESMF_FieldGather(fcstField(i), array_r8_3d, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
                  if (mype==0) then
                     ncerr = nf90_put_var(ncid, varids(i), values=array_r8_3d, start=start_idx); NC_ERR_STOP(ncerr)
                  end if
               end if
            end if
         end if ! end typekind

      else

         write(0,*)'Unsupported rank ', rank
         call ESMF_Finalize(endflag=ESMF_END_ABORT)

      end if ! end rank

    end do ! end fieldCount

    if (ideflate(grid_id) > 0 .and. nbits(grid_id) > 0 .and. do_io) then
       ncerr = nf90_redef(ncid=ncid); NC_ERR_STOP(ncerr)
       do i=1, fieldCount
          if (compress_err(i) > 0) then
             ncerr = nf90_put_att(ncid, varids(i), 'max_abs_compression_error', compress_err(i)); NC_ERR_STOP(ncerr)
             ncerr = nf90_put_att(ncid, varids(i), 'nbits', nbits(grid_id)); NC_ERR_STOP(ncerr)
          end if
       end do
       ncerr = nf90_enddef(ncid=ncid); NC_ERR_STOP(ncerr)
    end if

    if (.not. par) then
       deallocate(array_r4)
       deallocate(array_r8)
       deallocate(array_r4_3d)
       deallocate(array_r8_3d)
       if (is_cubed_sphere) then
          deallocate(array_r4_cube)
          deallocate(array_r8_cube)
          deallocate(array_r4_3d_cube)
          deallocate(array_r8_3d_cube)
       end if
    end if

    if (do_io) then
       deallocate(dimids_2d)
       deallocate(dimids_3d)
    end if

    deallocate(fcstField)
    deallocate(varids)
    deallocate(compress_err)

    if (do_io) then
       ncerr = nf90_close(ncid=ncid); NC_ERR_STOP(ncerr)
    end if

  end subroutine write_netcdf

!----------------------------------------------------------------------------------------
  subroutine get_global_attr(fldbundle, ncid, rc)
    type(ESMF_FieldBundle), intent(in) :: fldbundle
    integer, intent(in)                :: ncid
    integer, intent(out)               :: rc

! local variable
    integer :: i, attCount
    integer :: ncerr
    character(len=ESMF_MAXSTR) :: attName
    type(ESMF_TypeKind_Flag)   :: typekind

    integer :: varival
    real(ESMF_KIND_R4), dimension(:), allocatable :: varr4list
    real(ESMF_KIND_R8), dimension(:), allocatable :: varr8list
    integer :: itemCount
    character(len=ESMF_MAXSTR) :: varcval
!
    call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                           attnestflag=ESMF_ATTNEST_OFF, count=attCount, &
                           rc=rc); ESMF_ERR_RETURN(rc)

    do i=1,attCount

      call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                             attnestflag=ESMF_ATTNEST_OFF, attributeIndex=i, name=attName, &
                             typekind=typekind, itemCount=itemCount, rc=rc); ESMF_ERR_RETURN(rc)

      if (typekind==ESMF_TYPEKIND_I4) then
         call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                                name=trim(attName), value=varival, rc=rc); ESMF_ERR_RETURN(rc)
         ncerr = nf90_put_att(ncid, NF90_GLOBAL, trim(attName), varival); NC_ERR_STOP(ncerr)

      else if (typekind==ESMF_TYPEKIND_R4) then
         allocate (varr4list(itemCount))
         call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                                name=trim(attName), valueList=varr4list, rc=rc); ESMF_ERR_RETURN(rc)
         ncerr = nf90_put_att(ncid, NF90_GLOBAL, trim(attName), varr4list); NC_ERR_STOP(ncerr)
         deallocate(varr4list)

      else if (typekind==ESMF_TYPEKIND_R8) then
         allocate (varr8list(itemCount))
         call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                                name=trim(attName), valueList=varr8list, rc=rc); ESMF_ERR_RETURN(rc)
         ncerr = nf90_put_att(ncid, NF90_GLOBAL, trim(attName), varr8list); NC_ERR_STOP(ncerr)
         deallocate(varr8list)

      else if (typekind==ESMF_TYPEKIND_CHARACTER) then
         call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                                name=trim(attName), value=varcval, rc=rc); ESMF_ERR_RETURN(rc)
         ncerr = nf90_put_att(ncid, NF90_GLOBAL, trim(attName), trim(varcval)); NC_ERR_STOP(ncerr)

      end if

    end do

  end subroutine get_global_attr

!----------------------------------------------------------------------------------------
  subroutine get_grid_attr(grid, prefix, ncid, varid, rc)
    type(ESMF_Grid), intent(in)  :: grid
    character(len=*), intent(in) :: prefix
    integer, intent(in)          :: ncid
    integer, intent(in)          :: varid
    integer, intent(out)         :: rc

! local variable
    integer :: i, attCount, n, ind
    integer :: ncerr
    character(len=ESMF_MAXSTR) :: attName
    type(ESMF_TypeKind_Flag)   :: typekind

    integer :: varival
    real(ESMF_KIND_R4) :: varr4val
    real(ESMF_KIND_R8) :: varr8val
    character(len=ESMF_MAXSTR) :: varcval
!
    call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                           attnestflag=ESMF_ATTNEST_OFF, count=attCount, &
                           rc=rc); ESMF_ERR_RETURN(rc)

    do i=1,attCount

      call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                             attnestflag=ESMF_ATTNEST_OFF, attributeIndex=i, name=attName, &
                             typekind=typekind, itemCount=n, rc=rc); ESMF_ERR_RETURN(rc)

      if (index(trim(attName), trim(prefix)//":")==1) then
         ind = len(trim(prefix)//":")

         if (typekind==ESMF_TYPEKIND_I4) then
            call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=varival, rc=rc); ESMF_ERR_RETURN(rc)
            ncerr = nf90_put_att(ncid, varid, trim(attName(ind+1:len(attName))), varival); NC_ERR_STOP(ncerr)

         else if (typekind==ESMF_TYPEKIND_R4) then
            call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=varr4val, rc=rc); ESMF_ERR_RETURN(rc)
            ncerr = nf90_put_att(ncid, varid, trim(attName(ind+1:len(attName))), varr4val); NC_ERR_STOP(ncerr)

         else if (typekind==ESMF_TYPEKIND_R8) then
            call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=varr8val, rc=rc); ESMF_ERR_RETURN(rc)
            if (trim(attName) /= '_FillValue') then
              ! FIXME:  _FillValue must be cast to var type when using
              ! NF90_NETCDF4. Until this is fixed, using netCDF default _FillValue.
              ncerr = nf90_put_att(ncid, varid, trim(attName(ind+1:len(attName))), varr8val); NC_ERR_STOP(ncerr)
            end if

         else if (typekind==ESMF_TYPEKIND_CHARACTER) then
            call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=varcval, rc=rc); ESMF_ERR_RETURN(rc)
            ncerr = nf90_put_att(ncid, varid, trim(attName(ind+1:len(attName))), trim(varcval)); NC_ERR_STOP(ncerr)

         end if

      end if

    end do

  end subroutine get_grid_attr

!----------------------------------------------------------------------------------------
  subroutine add_dim(ncid, dim_name, dimid, grid, rc)
    integer, intent(in)             :: ncid
    character(len=*), intent(in)    :: dim_name
    integer, intent(inout) :: dimid
    type(ESMF_Grid), intent(in)     :: grid
    integer, intent(out)            :: rc

! local variable
    integer :: n, dim_varid
    integer :: ncerr
    type(ESMF_TypeKind_Flag)   :: typekind

    real(ESMF_KIND_R4), allocatable  :: valueListR4(:)
    real(ESMF_KIND_R8), allocatable  :: valueListR8(:)
!
    call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                           attnestflag=ESMF_ATTNEST_OFF, name=dim_name, &
                           typekind=typekind, itemCount=n, rc=rc); ESMF_ERR_RETURN(rc)

    if (trim(dim_name) == "time") then
    ! using an unlimited dim requires collective mode (NF90_COLLECTIVE)
    ! for parallel writes, which seems to slow things down on hera.
    !ncerr = nf90_def_dim(ncid, trim(dim_name), NF90_UNLIMITED, dimid); NC_ERR_STOP(ncerr)
    ncerr = nf90_def_dim(ncid, trim(dim_name), 1, dimid); NC_ERR_STOP(ncerr)
    else
    ncerr = nf90_def_dim(ncid, trim(dim_name), n, dimid); NC_ERR_STOP(ncerr)
    end if

    if (typekind==ESMF_TYPEKIND_R8) then
       ncerr = nf90_def_var(ncid, dim_name, NF90_REAL8, dimids=[dimid], varid=dim_varid); NC_ERR_STOP(ncerr)
       allocate(valueListR8(n))
       call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                              name=trim(dim_name), valueList=valueListR8, rc=rc); ESMF_ERR_RETURN(rc)
       ncerr = nf90_enddef(ncid=ncid); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_var(ncid, dim_varid, values=valueListR8); NC_ERR_STOP(ncerr)
       ncerr = nf90_redef(ncid=ncid); NC_ERR_STOP(ncerr)
       deallocate(valueListR8)
     else if (typekind==ESMF_TYPEKIND_R4) then
       ncerr = nf90_def_var(ncid, dim_name, NF90_REAL4, dimids=[dimid], varid=dim_varid); NC_ERR_STOP(ncerr)
       allocate(valueListR4(n))
       call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                              name=trim(dim_name), valueList=valueListR4, rc=rc); ESMF_ERR_RETURN(rc)
       ncerr = nf90_enddef(ncid=ncid); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_var(ncid, dim_varid, values=valueListR4); NC_ERR_STOP(ncerr)
       ncerr = nf90_redef(ncid=ncid); NC_ERR_STOP(ncerr)
       deallocate(valueListR4)
     else
        write(0,*)'Error in module_write_netcdf.F90(add_dim) unknown typekind for ',trim(dim_name)
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    if (par) then
       ncerr = nf90_var_par_access(ncid, dim_varid, NF90_INDEPENDENT); NC_ERR_STOP(ncerr)
    end if

    call get_grid_attr(grid, dim_name, ncid, dim_varid, rc)

  end subroutine add_dim

!----------------------------------------------------------------------------------------
  subroutine quantize_array_3d(array, dataMin, dataMax, nbits, compress_err)

    real(4), dimension(:,:,:), intent(inout)   :: array
    real(4), intent(in)                        :: dataMin, dataMax
    integer, intent(in)                        :: nbits
    real(4), intent(out)                       :: compress_err

    real(4) :: scale_fact, offset
    real(4), dimension(:,:,:), allocatable     :: array_save
    ! Lossy compression if nbits>0.
    ! The floating point data is quantized to improve compression
    ! See doi:10.5194/gmd-10-413-2017.  The method employed
    ! here is identical to the 'scaled linear packing' method in
    ! that paper, except that the data are scaling into an arbitrary
    ! range (2**nbits-1 not just 2**16-1) and are stored as
    ! re-scaled floats instead of short integers.
    ! The zlib algorithm does almost as
    ! well packing the re-scaled floats as it does the scaled
    ! integers, and this avoids the need for the client to apply the
    ! rescaling (plus it allows the ability to adjust the packing
    ! range).
    scale_fact = (dataMax - dataMin) / (2**nbits-1)
    offset = dataMin
    if (scale_fact > 0.) then
       allocate(array_save, source=array)
       array = scale_fact*(nint((array_save - offset) / scale_fact)) + offset
       ! compute max abs compression error
       compress_err = maxval(abs(array_save-array))
       deallocate(array_save)
    else
       ! field is constant
       compress_err = 0.
    end if
  end subroutine quantize_array_3d

  subroutine quantize_array_4d(array, dataMin, dataMax, nbits, compress_err)

    real(4), dimension(:,:,:,:), intent(inout) :: array
    real(4), intent(in)                        :: dataMin, dataMax
    integer, intent(in)                        :: nbits
    real(4), intent(out)                       :: compress_err

    real(4) :: scale_fact, offset
    real(4), dimension(:,:,:,:), allocatable   :: array_save

    ! Lossy compression if nbits>0.
    ! The floating point data is quantized to improve compression
    ! See doi:10.5194/gmd-10-413-2017.  The method employed
    ! here is identical to the 'scaled linear packing' method in
    ! that paper, except that the data are scaling into an arbitrary
    ! range (2**nbits-1 not just 2**16-1) and are stored as
    ! re-scaled floats instead of short integers.
    ! The zlib algorithm does almost as
    ! well packing the re-scaled floats as it does the scaled
    ! integers, and this avoids the need for the client to apply the
    ! rescaling (plus it allows the ability to adjust the packing
    ! range).
    scale_fact = (dataMax - dataMin) / (2**nbits-1)
    offset = dataMin
    if (scale_fact > 0.) then
       allocate(array_save, source=array)
       array = scale_fact*(nint((array_save - offset) / scale_fact)) + offset
       ! compute max abs compression error
       compress_err = maxval(abs(array_save-array))
       deallocate(array_save)
    else
       ! field is constant
       compress_err = 0.
    end if
  end subroutine quantize_array_4d

!----------------------------------------------------------------------------------------
end module module_write_netcdf
