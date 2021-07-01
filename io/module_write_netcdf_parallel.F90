#define ESMF_ERR_RETURN(rc) if (ESMF_LogFoundError(rc, msg="Breaking out of subroutine", line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

#define NC_ERR_STOP(status) \
    if (status /= nf90_noerr) write(0,*) "line ", __LINE__, trim(nf90_strerror(status)); \
    if (status /= nf90_noerr) call ESMF_Finalize(endflag=ESMF_END_ABORT)

module module_write_netcdf_parallel

  use esmf
  use netcdf
  use module_fv3_io_def,only : ideflate, nbits, &
                               output_grid,dx,dy,lon1,lat1,lon2,lat2
  use mpi

  implicit none
  private
  public write_netcdf_parallel

  contains

#ifdef NO_PARALLEL_NETCDF
!----------------------------------------------------------------------------------------
  subroutine write_netcdf_parallel(fieldbundle, wrtfb, filename, mpi_comm, mype, im, jm, ichunk2d, jchunk2d, ichunk3d, jchunk3d, kchunk3d, rc)
    type(ESMF_FieldBundle), intent(in) :: fieldbundle
    type(ESMF_FieldBundle), intent(in) :: wrtfb
    character(*), intent(in)           :: filename
    integer, intent(in)                :: mpi_comm
    integer, intent(in)                :: mype
    integer, intent(in)                :: im, jm, ichunk2d, jchunk2d, &
                                          ichunk3d, jchunk3d, kchunk3d
    integer, optional,intent(out)      :: rc
    print *,'in stub write_netcdf_parallel - model not built with parallel netcdf support, return'
  end subroutine write_netcdf_parallel
#else
!----------------------------------------------------------------------------------------
  subroutine write_netcdf_parallel(fieldbundle, wrtfb, filename, mpi_comm, mype, im, jm, ichunk2d, jchunk2d, ichunk3d, jchunk3d, kchunk3d, rc)
!
    type(ESMF_FieldBundle), intent(in) :: fieldbundle
    type(ESMF_FieldBundle), intent(in) :: wrtfb
    character(*), intent(in)           :: filename
    integer, intent(in)                :: mpi_comm
    integer, intent(in)                :: mype
    integer, intent(in)                :: im, jm, ichunk2d, jchunk2d, &
                                          ichunk3d, jchunk3d, kchunk3d
    integer, optional,intent(out)      :: rc
!
!** local vars
    integer :: i,j,m,n,k,istart,iend,jstart,jend,i1,i2,j1,j2,k1,k2
    integer :: lm

    integer, dimension(:), allocatable     :: fldlev
    real(ESMF_KIND_R4), dimension(:,:), pointer   :: arrayr4
    real(ESMF_KIND_R8), dimension(:,:), pointer   :: arrayr8
    real(ESMF_KIND_R4), dimension(:,:,:), pointer :: arrayr4_3d,arrayr4_3d_save
    real(ESMF_KIND_R8), dimension(:,:,:), pointer :: arrayr8_3d

    real(8) x(im),y(jm)
    integer :: fieldCount, fieldDimCount, gridDimCount
    integer, dimension(:), allocatable   :: ungriddedLBound, ungriddedUBound

    type(ESMF_Field), allocatable        :: fcstField(:)
    type(ESMF_TypeKind_Flag)             :: typekind
    type(ESMF_TypeKind_Flag)             :: attTypeKind
    type(ESMF_Grid)                      :: wrtgrid
    type(ESMF_Array)                     :: array

    integer :: attcount
    character(len=ESMF_MAXSTR) :: attName, fldName
    integer :: totalLBound2d(2),totalUBound2d(2),totalLBound3d(3),totalUBound3d(3)

    integer :: varival
    real(4) :: varr4val, scale_fact, offset, dataMin, dataMax
    real(4), allocatable, dimension(:) :: compress_err
    real(8) :: varr8val
    character(len=ESMF_MAXSTR) :: varcval

    character(128) :: time_units

    integer :: ncerr,ierr
    integer :: ncid
    integer :: oldMode
    integer :: im_dimid, jm_dimid, pfull_dimid, phalf_dimid, time_dimid
    integer :: im_varid, jm_varid, lm_varid, time_varid, lon_varid, lat_varid
    integer, dimension(:), allocatable :: varids
    logical shuffle
!
    call ESMF_FieldBundleGet(fieldbundle, fieldCount=fieldCount, rc=rc); ESMF_ERR_RETURN(rc)

    allocate(compress_err(fieldCount)); compress_err=-999.
    allocate(fldlev(fieldCount)) ; fldlev = 0
    allocate(fcstField(fieldCount))
    allocate(varids(fieldCount))

    call ESMF_FieldBundleGet(fieldbundle, fieldList=fcstField, grid=wrtGrid, &
!                             itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                             rc=rc); ESMF_ERR_RETURN(rc)

    call ESMF_GridGet(wrtgrid, dimCount=gridDimCount, rc=rc); ESMF_ERR_RETURN(rc)

    do i=1,fieldCount
       call ESMF_FieldGet(fcstField(i), dimCount=fieldDimCount, rc=rc); ESMF_ERR_RETURN(rc)
       if (fieldDimCount > 3) then
          write(0,*)"write_netcdf: Only 2D and 3D fields are supported!"
          stop
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

! create netcdf file for parallel access

    ncerr = nf90_create(trim(filename),&
            cmode=IOR(IOR(NF90_CLOBBER,NF90_NETCDF4),NF90_CLASSIC_MODEL),&
            comm=mpi_comm, info = MPI_INFO_NULL, ncid=ncid); NC_ERR_STOP(ncerr)
! disable auto filling.
    ncerr = nf90_set_fill(ncid, NF90_NOFILL, oldMode); NC_ERR_STOP(ncerr)

    ! define dimensions
    ncerr = nf90_def_dim(ncid, "grid_xt", im, im_dimid); NC_ERR_STOP(ncerr)
    ncerr = nf90_def_dim(ncid, "grid_yt", jm, jm_dimid); NC_ERR_STOP(ncerr)
    ! define coordinate variables
    ncerr = nf90_def_var(ncid, "grid_xt", NF90_DOUBLE, im_dimid, im_varid); NC_ERR_STOP(ncerr)
    ncerr = nf90_var_par_access(ncid, im_varid, NF90_INDEPENDENT)
    ncerr = nf90_def_var(ncid, "lon", NF90_DOUBLE, (/im_dimid,jm_dimid/), lon_varid); NC_ERR_STOP(ncerr)
    !ncerr = nf90_var_par_access(ncid, lon_varid, NF90_INDEPENDENT)
    ncerr = nf90_put_att(ncid, lon_varid, "long_name", "T-cell longitude"); NC_ERR_STOP(ncerr)
    ncerr = nf90_put_att(ncid, lon_varid, "units", "degrees_E"); NC_ERR_STOP(ncerr)
    ncerr = nf90_put_att(ncid, im_varid, "cartesian_axis", "X"); NC_ERR_STOP(ncerr)
    ncerr = nf90_def_var(ncid, "grid_yt", NF90_DOUBLE, jm_dimid, jm_varid); NC_ERR_STOP(ncerr)
    ncerr = nf90_var_par_access(ncid, jm_varid, NF90_INDEPENDENT)
    ncerr = nf90_def_var(ncid, "lat", NF90_DOUBLE, (/im_dimid,jm_dimid/), lat_varid); NC_ERR_STOP(ncerr)
    ncerr = nf90_var_par_access(ncid, lat_varid, NF90_INDEPENDENT)
    ncerr = nf90_put_att(ncid, lat_varid, "long_name", "T-cell latitude"); NC_ERR_STOP(ncerr)
    ncerr = nf90_put_att(ncid, lat_varid, "units", "degrees_N"); NC_ERR_STOP(ncerr)
    ncerr = nf90_put_att(ncid, jm_varid, "cartesian_axis", "Y"); NC_ERR_STOP(ncerr)

    if (lm > 1) then
      call add_dim(ncid, "pfull", pfull_dimid, wrtgrid, rc)
      call add_dim(ncid, "phalf", phalf_dimid, wrtgrid, rc)
    end if

    call add_dim(ncid, "time", time_dimid, wrtgrid, rc)

    call get_global_attr(wrtfb, ncid, rc)

    do i=1, fieldCount
      call ESMF_FieldGet(fcstField(i), name=fldName, typekind=typekind, rc=rc); ESMF_ERR_RETURN(rc)

      ! define variables
      if (fldlev(i) == 1) then
        if (typekind == ESMF_TYPEKIND_R4) then
          if (ideflate > 0) then
            if (ichunk2d < 0 .or. jchunk2d < 0) then
               ! let netcdf lib choose chunksize
               ! shuffle filter on for 2d fields (lossless compression)
               ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
                       (/im_dimid,jm_dimid,time_dimid/), varids(i), &
                       shuffle=.true.,deflate_level=ideflate); NC_ERR_STOP(ncerr)
            else
               ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
                       (/im_dimid,jm_dimid,time_dimid/), varids(i), &
                       shuffle=.true.,deflate_level=ideflate,&
                       chunksizes=(/ichunk2d,jchunk2d,1/)); NC_ERR_STOP(ncerr)
            endif
            ! compression filters require collective access.
            ncerr = nf90_var_par_access(ncid, varids(i), NF90_COLLECTIVE) 
          else
            ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
            (/im_dimid,jm_dimid,time_dimid/), varids(i)); NC_ERR_STOP(ncerr)
            ncerr = nf90_var_par_access(ncid, varids(i), NF90_INDEPENDENT) 
          endif
        else if (typekind == ESMF_TYPEKIND_R8) then
          ncerr = nf90_def_var(ncid, trim(fldName), NF90_DOUBLE, &
                               (/im_dimid,jm_dimid,time_dimid/), varids(i)); NC_ERR_STOP(ncerr)
          ncerr = nf90_var_par_access(ncid, varids(i), NF90_INDEPENDENT) 
        else
          write(0,*)'Unsupported typekind ', typekind
          stop
        end if
      else if (fldlev(i) > 1) then
        if (typekind == ESMF_TYPEKIND_R4) then
          if (ideflate > 0) then
            ! shuffle filter off for 3d fields using lossy compression
            if (nbits > 0) then
                shuffle=.false.
            else
                shuffle=.true.
            endif
            if (ichunk3d < 0 .or. jchunk3d < 0 .or. kchunk3d < 0) then
               ! let netcdf lib choose chunksize
               ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
                       (/im_dimid,jm_dimid,pfull_dimid,time_dimid/), varids(i), &
                       shuffle=shuffle,deflate_level=ideflate); NC_ERR_STOP(ncerr)
            else
               ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
                       (/im_dimid,jm_dimid,pfull_dimid,time_dimid/), varids(i), &
                       shuffle=shuffle,deflate_level=ideflate,&
                       chunksizes=(/ichunk3d,jchunk3d,kchunk3d,1/)); NC_ERR_STOP(ncerr)
            endif
            ! compression filters require collective access.
            ncerr = nf90_var_par_access(ncid, varids(i), NF90_COLLECTIVE) 
          else
            ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, &
            (/im_dimid,jm_dimid,pfull_dimid,time_dimid/), varids(i)); NC_ERR_STOP(ncerr)
            ncerr = nf90_var_par_access(ncid, varids(i), NF90_INDEPENDENT) 
          endif
        else if (typekind == ESMF_TYPEKIND_R8) then
          ncerr = nf90_def_var(ncid, trim(fldName), NF90_DOUBLE, &
                                (/im_dimid,jm_dimid,pfull_dimid,time_dimid/), varids(i)); NC_ERR_STOP(ncerr)
          ncerr = nf90_var_par_access(ncid, varids(i), NF90_INDEPENDENT) 
        else
          write(0,*)'Unsupported typekind ', typekind
          stop
        end if
      end if

      ! define variable attributes
      call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                             attnestflag=ESMF_ATTNEST_OFF, Count=attcount, &
                             rc=rc); ESMF_ERR_RETURN(rc)

      do j=1,attCount
        call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                               attnestflag=ESMF_ATTNEST_OFF, attributeIndex=j, &
                               name=attName, typekind=attTypeKind, itemCount=n, &
                               rc=rc); ESMF_ERR_RETURN(rc)

        if ( index(trim(attName),"ESMF") /= 0 ) then
           cycle
        endif

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
           endif

        else if (attTypeKind==ESMF_TYPEKIND_CHARACTER) then
           call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                                  name=trim(attName), value=varcval, &
                                  rc=rc); ESMF_ERR_RETURN(rc)
           ncerr = nf90_put_att(ncid, varids(i), trim(attName), trim(varcval)); NC_ERR_STOP(ncerr)

        end if

      end do ! j=1,attCount

    end do   ! i=1,fieldCount

    ! write grid_xt, grid_yt attributes
    if (trim(output_grid) == 'gaussian_grid' .or. &
        trim(output_grid) == 'global_latlon' .or. &
        trim(output_grid) == 'regional_latlon') then
       ncerr = nf90_put_att(ncid, im_varid, "long_name", "T-cell longitude"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, im_varid, "units", "degrees_E"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, jm_varid, "long_name", "T-cell latiitude"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, jm_varid, "units", "degrees_N"); NC_ERR_STOP(ncerr)
    else if (trim(output_grid) == 'rotated_latlon') then
       ncerr = nf90_put_att(ncid, im_varid, "long_name", "rotated T-cell longiitude"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, im_varid, "units", "degrees"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, jm_varid, "long_name", "rotated T-cell latiitude"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, jm_varid, "units", "degrees"); NC_ERR_STOP(ncerr)
    else if (trim(output_grid) == 'lambert_conformal') then
       ncerr = nf90_put_att(ncid, im_varid, "long_name", "x-coordinate of projection"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, im_varid, "units", "meters"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, jm_varid, "long_name", "y-coordinate of projection"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, jm_varid, "units", "meters"); NC_ERR_STOP(ncerr)
    endif

    ncerr = nf90_enddef(ncid); NC_ERR_STOP(ncerr)

! end of define mode

    ! write grid_xt, grid_yt values
    call ESMF_GridGetCoord(wrtGrid, coordDim=1, farrayPtr=arrayr8, rc=rc); ESMF_ERR_RETURN(rc)
    istart = lbound(arrayr8,1); iend   = ubound(arrayr8,1)
    jstart = lbound(arrayr8,2); jend   = ubound(arrayr8,2)
    !print *,'in write netcdf mpi dim 1',istart,iend,jstart,jend,shape(arrayr8),minval(arrayr8(:,jstart)),maxval(arrayr8(:,jstart))

    if (trim(output_grid) == 'gaussian_grid' .or. &
        trim(output_grid) == 'global_latlon' .or. &
        trim(output_grid) == 'regional_latlon') then
      ncerr = nf90_put_var(ncid, im_varid, values=arrayr8(:,jstart),start=(/istart/), count=(/iend-istart+1/)); NC_ERR_STOP(ncerr)
    else if (trim(output_grid) == 'rotated_latlon') then
      do i=1,im
         x(i) = lon1 + (lon2-lon1)/(im-1) * (i-1)
      enddo
      ncerr = nf90_put_var(ncid, im_varid, values=x  ); NC_ERR_STOP(ncerr)
    else if (trim(output_grid) == 'lambert_conformal') then
      do i=1,im
         x(i) = dx * (i-1)
      enddo
      ncerr = nf90_put_var(ncid, im_varid, values=x  ); NC_ERR_STOP(ncerr)
    endif
    ncerr = nf90_put_var(ncid, lon_varid, values=arrayr8, start=(/istart,jstart/)); NC_ERR_STOP(ncerr)

    call ESMF_GridGetCoord(wrtGrid, coordDim=2, farrayPtr=arrayr8, rc=rc); ESMF_ERR_RETURN(rc)
    !print *,'in write netcdf mpi dim 2',istart,iend,jstart,jend,shape(arrayr8),minval(arrayr8(istart,:)),maxval(arrayr8(istart,:))
    if (trim(output_grid) == 'gaussian_grid' .or. &
        trim(output_grid) == 'global_latlon' .or. &
        trim(output_grid) == 'regional_latlon') then
          ncerr = nf90_put_var(ncid, jm_varid, values=arrayr8(istart,:),start=(/jstart/),count=(/jend-jstart+1/)); NC_ERR_STOP(ncerr)
    else if (trim(output_grid) == 'rotated_latlon') then
          do j=1,jm
             y(j) = lat1 + (lat2-lat1)/(jm-1) * (j-1)
          enddo
          ncerr = nf90_put_var(ncid, jm_varid, values=y  ); NC_ERR_STOP(ncerr)
    else if (trim(output_grid) == 'lambert_conformal') then
          do j=1,jm
             y(j) = dy * (j-1)
          enddo
          ncerr = nf90_put_var(ncid, jm_varid, values=y  ); NC_ERR_STOP(ncerr)
    endif
    ncerr = nf90_put_var(ncid, lat_varid, values=arrayr8, start=(/istart,jstart/)); NC_ERR_STOP(ncerr)

    do i=1, fieldCount

       call ESMF_FieldGet(fcstField(i),name=fldName,typekind=typekind, rc=rc); ESMF_ERR_RETURN(rc)

       if (fldlev(i) == 1) then
         if (typekind == ESMF_TYPEKIND_R4) then
           call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=arrayr4, totalLBound=totalLBound2d, totalUBound=totalUBound2d,rc=rc); ESMF_ERR_RETURN(rc)
           !print *,'field name=',trim(fldName),'bound=',totalLBound2d,'ubound=',totalUBound2d
           ncerr = nf90_put_var(ncid, varids(i), values=arrayr4, start=(/totalLBound2d(1),totalLBound2d(2),1/)); NC_ERR_STOP(ncerr)
         else if (typekind == ESMF_TYPEKIND_R8) then
           call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=arrayr8, totalLBound=totalLBound2d, totalUBound=totalUBound2d,rc=rc); ESMF_ERR_RETURN(rc)
           ncerr = nf90_put_var(ncid, varids(i), values=arrayr8, start=(/totalLBound2d(1),totalLBound2d(2),1/)); NC_ERR_STOP(ncerr)
         end if
      else if (fldlev(i) > 1) then
         if (typekind == ESMF_TYPEKIND_R4) then
           call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=arrayr4_3d, totalLBound=totalLBound3d, totalUBound=totalUBound3d,rc=rc); ESMF_ERR_RETURN(rc)
           if (ideflate > 0 .and. nbits > 0) then
              i1=totalLBound3d(1);i2=totalUBound3d(1)
              j1=totalLBound3d(2);j2=totalUBound3d(2)
              k1=totalLBound3d(3);k2=totalUBound3d(3)
              dataMax = maxval(arrayr4_3d(i1:i2,j1:j2,k1:k2))
              dataMin = minval(arrayr4_3d(i1:i2,j1:j2,k1:k2))
              call mpi_allreduce(mpi_in_place,dataMax,1,mpi_real4,mpi_max,mpi_comm,ierr)
              call mpi_allreduce(mpi_in_place,dataMin,1,mpi_real4,mpi_min,mpi_comm,ierr)
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
              ! range)
              scale_fact = (dataMax - dataMin) / (2**nbits-1); offset = dataMin
              if (scale_fact > 0.) then
                allocate(arrayr4_3d_save(i1:i2,j1:j2,k1:k2))
                arrayr4_3d_save(i1:i2,j1:j2,k1:k2)=arrayr4_3d(i1:i2,j1:j2,k1:k2)
                arrayr4_3d = scale_fact*(nint((arrayr4_3d_save - offset) / scale_fact)) + offset
                ! compute max abs compression error.
                compress_err(i) = &
                maxval(abs(arrayr4_3d_save(i1:i2,j1:j2,k1:k2)-arrayr4_3d(i1:i2,j1:j2,k1:k2)))
                deallocate(arrayr4_3d_save)
                call mpi_allreduce(mpi_in_place,compress_err(i),1,mpi_real4,mpi_max,mpi_comm,ierr)
                !print *,'field name=',trim(fldName),dataMin,dataMax,compress_err(i)
              else
                ! field is constant
                compress_err(i) = 0.
              endif
           endif
           ncerr = nf90_put_var(ncid, varids(i), values=arrayr4_3d, start=(/totalLBound3d(1),totalLBound3d(2),totalLBound3d(3),1/)); NC_ERR_STOP(ncerr)
         else if (typekind == ESMF_TYPEKIND_R8) then
           call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=arrayr8_3d, totalLBound=totalLBound3d, totalUBound=totalUBound3d,rc=rc); ESMF_ERR_RETURN(rc)
           !print *,'field name=',trim(fldName),'bound=',totalLBound3d,'ubound=',totalUBound3d
           ncerr = nf90_put_var(ncid, varids(i), values=arrayr8_3d, start=(/totalLBound3d(1),totalLBound3d(2),totalLBound3d(3),1/)); NC_ERR_STOP(ncerr)
         end if

      end if  !end fldlev(i)

    end do  ! end fieldCount

    if (ideflate > 0 .and. nbits > 0) then
       ncerr = nf90_redef(ncid=ncid); NC_ERR_STOP(ncerr)
       do i=1, fieldCount
          if (compress_err(i) > 0) then
             ncerr = nf90_put_att(ncid, varids(i), 'max_abs_compression_error', compress_err(i)); NC_ERR_STOP(ncerr)
             ncerr = nf90_put_att(ncid, varids(i), 'nbits', nbits); NC_ERR_STOP(ncerr)
          endif
       enddo
       ncerr = nf90_enddef(ncid=ncid); NC_ERR_STOP(ncerr)
    endif

    deallocate(fcstField)
    deallocate(varids)
    deallocate(compress_err)

    ncerr = nf90_close(ncid=ncid); NC_ERR_STOP(ncerr)
    !call mpi_barrier(mpi_comm,ierr)
    !print *,'netcdf parallel close, finished write_netcdf_parallel'

  end subroutine write_netcdf_parallel
#endif

!----------------------------------------------------------------------------------------
  subroutine get_global_attr(fldbundle, ncid, rc)
    type(ESMF_FieldBundle), intent(in) :: fldbundle
    integer, intent(in)                :: ncid
    integer, intent(out)               :: rc

! local variable
    integer :: i, attcount
    integer :: ncerr
    character(len=ESMF_MAXSTR) :: attName
    type(ESMF_TypeKind_Flag)   :: typekind

    integer :: varival
    real(ESMF_KIND_R4) :: varr4val
    real(ESMF_KIND_R4), dimension(:), allocatable :: varr4list
    real(ESMF_KIND_R8) :: varr8val
    real(ESMF_KIND_R8), dimension(:), allocatable :: varr8list
    integer :: itemCount
    character(len=ESMF_MAXSTR) :: varcval
!
    call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                           attnestflag=ESMF_ATTNEST_OFF, Count=attcount, &
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
!
!----------------------------------------------------------------------------------------
  subroutine get_grid_attr(grid, prefix, ncid, varid, rc)
    type(ESMF_Grid), intent(in)  :: grid
    character(len=*), intent(in) :: prefix
    integer, intent(in)          :: ncid
    integer, intent(in)          :: varid
    integer, intent(out)         :: rc

! local variable
    integer :: i, attcount, n, ind
    integer :: ncerr
    character(len=ESMF_MAXSTR) :: attName
    type(ESMF_TypeKind_Flag)   :: typekind

    integer :: varival
    real(ESMF_KIND_R4) :: varr4val
    real(ESMF_KIND_R8) :: varr8val
    character(len=ESMF_MAXSTR) :: varcval
!
    call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                           attnestflag=ESMF_ATTNEST_OFF, Count=attcount, &
                           rc=rc); ESMF_ERR_RETURN(rc)

    !write(0,*)'grid attcount = ', attcount
    do i=1,attCount

      call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                             attnestflag=ESMF_ATTNEST_OFF, attributeIndex=i, name=attName, &
                             typekind=typekind, itemCount=n, rc=rc); ESMF_ERR_RETURN(rc)
      !write(0,*)'grid att = ',i,trim(attName), ' itemCount = ' , n

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
            endif

         else if (typekind==ESMF_TYPEKIND_CHARACTER) then
            call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=varcval, rc=rc); ESMF_ERR_RETURN(rc)
            ncerr = nf90_put_att(ncid, varid, trim(attName(ind+1:len(attName))), trim(varcval)); NC_ERR_STOP(ncerr)

         end if

      end if

    end do

  end subroutine get_grid_attr

  subroutine add_dim(ncid, dim_name, dimid, grid, rc)
    integer, intent(in)             :: ncid
    character(len=*), intent(in)    :: dim_name
    integer, intent(inout) :: dimid
    type(ESMF_Grid), intent(in)     :: grid
    integer, intent(out)            :: rc

! local variable
    integer :: i, attcount, n, dim_varid
    integer :: ncerr
    character(len=ESMF_MAXSTR) :: attName
    type(ESMF_TypeKind_Flag)   :: typekind

    integer, allocatable  :: valueListI(:)
    real(ESMF_KIND_R4), allocatable  :: valueListR4(:)
    real(ESMF_KIND_R8), allocatable  :: valueListR8(:)
    character(len=ESMF_MAXSTR), allocatable  :: valueListC(:)
!
    call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                           attnestflag=ESMF_ATTNEST_OFF, name=dim_name, &
                           typekind=typekind, itemCount=n, rc=rc); ESMF_ERR_RETURN(rc)

    if ( trim(dim_name) == "time" ) then
    ! using an unlimited dim requires collective mode (NF90_COLLECTIVE)
    ! for parallel writes, which seems to slow things down on hera.
    !ncerr = nf90_def_dim(ncid, trim(dim_name), NF90_UNLIMITED, dimid); NC_ERR_STOP(ncerr)
    ncerr = nf90_def_dim(ncid, trim(dim_name), 1, dimid); NC_ERR_STOP(ncerr)
    else
    ncerr = nf90_def_dim(ncid, trim(dim_name), n, dimid); NC_ERR_STOP(ncerr)
    end if

    if (typekind==ESMF_TYPEKIND_R8) then
       ncerr = nf90_def_var(ncid, dim_name, NF90_REAL8, dimids=(/dimid/), varid=dim_varid); NC_ERR_STOP(ncerr)
       ncerr = nf90_var_par_access(ncid, dim_varid, NF90_INDEPENDENT)
       allocate(valueListR8(n))
       call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                              name=trim(dim_name), valueList=valueListR8, rc=rc); ESMF_ERR_RETURN(rc)
       ncerr = nf90_enddef(ncid=ncid); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_var(ncid, dim_varid, values=valueListR8 ); NC_ERR_STOP(ncerr)
       ncerr = nf90_redef(ncid=ncid); NC_ERR_STOP(ncerr)
       deallocate(valueListR8)
     else if (typekind==ESMF_TYPEKIND_R4) then
       ncerr = nf90_def_var(ncid, dim_name, NF90_REAL4, dimids=(/dimid/), varid=dim_varid); NC_ERR_STOP(ncerr)
       ncerr = nf90_var_par_access(ncid, dim_varid, NF90_INDEPENDENT)
       allocate(valueListR4(n))
       call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                              name=trim(dim_name), valueList=valueListR4, rc=rc); ESMF_ERR_RETURN(rc)
       ncerr = nf90_enddef(ncid=ncid); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_var(ncid, dim_varid, values=valueListR4 ); NC_ERR_STOP(ncerr)
       ncerr = nf90_redef(ncid=ncid); NC_ERR_STOP(ncerr)
       deallocate(valueListR4)
     else
        write(0,*)'Error in module_write_netcdf.F90(add_dim) unknown typekind for ',trim(dim_name)
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    call get_grid_attr(grid, dim_name, ncid, dim_varid, rc)

  end subroutine add_dim
!
!----------------------------------------------------------------------------------------
  subroutine nccheck(status)
    use netcdf
    implicit none
    integer, intent (in) :: status

    if (status /= nf90_noerr) then
      write(0,*) status, trim(nf90_strerror(status))
      stop "stopped"
    end if
  end subroutine nccheck
 
end module module_write_netcdf_parallel
