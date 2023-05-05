module FV3GFS_common_io
  use GFS_typedefs, only: kind_phys

  implicit none
  private

  public :: copy_from_GFS_Data, copy_to_GFS_Data
  public :: copy_from_GFS_Data_2d_phys2phys, copy_from_GFS_Data_3d_phys2phys, &
       copy_from_GFS_Data_2d_int2phys, copy_from_GFS_Data_3d_int2phys, &
       copy_from_GFS_Data_2d_stack_phys2phys, copy_to_GFS_Data_3d_slice_phys2phys, &
       copy_to_GFS_Data_2d_phys2phys, copy_to_GFS_Data_3d_phys2phys, &
       copy_to_GFS_Data_2d_int2phys, copy_to_GFS_Data_3d_int2phys

  public :: create_2d_field_and_add_to_bundle
  public :: create_3d_field_and_add_to_bundle
  public :: add_zaxis_to_field

  interface copy_from_GFS_Data
    module procedure copy_from_GFS_Data_2d_phys2phys, &
         copy_from_GFS_Data_3d_phys2phys, &
         copy_from_GFS_Data_2d_int2phys, &
         copy_from_GFS_Data_3d_int2phys, &
         copy_from_GFS_Data_2d_stack_phys2phys
  end interface

  interface copy_to_GFS_Data
    module procedure copy_to_GFS_Data_2d_phys2phys, &
         copy_to_GFS_Data_3d_phys2phys, &
         copy_to_GFS_Data_2d_int2phys, &
         copy_to_GFS_Data_3d_int2phys, &
         copy_to_GFS_Data_3d_slice_phys2phys
  end interface copy_to_GFS_Data

contains

   pure subroutine copy_from_GFS_Data_2d_phys2phys(ii1,jj1,isc,jsc,nt,var2d,var_block)
     implicit none
     integer, intent(in) :: ii1(:), jj1(:), isc, jsc
     integer, intent(inout) :: nt
     real(kind=kind_phys), intent(in) :: var_block(:)
     real(kind=kind_phys), intent(out) :: var2d(:,:,:)
     integer ix

     nt=nt+1
     do ix=1,size(var_block)
       var2d(ii1(ix),jj1(ix),nt) = var_block(ix)
     enddo
   end subroutine copy_from_GFS_Data_2d_phys2phys

   pure subroutine copy_from_GFS_Data_3d_phys2phys(ii1,jj1,isc,jsc,nt,var3d,var_block)
     implicit none
     integer, intent(in) :: ii1(:), jj1(:), isc, jsc
     integer, intent(inout) :: nt
     real(kind=kind_phys), intent(in) :: var_block(:,:)
     real(kind=kind_phys), intent(out) :: var3d(:,:,:,:)
     integer ix, k

     nt=nt+1
     do k=lbound(var_block,2),ubound(var_block,2)
       do ix=1,size(var_block,1)
         var3d(ii1(ix),jj1(ix),k,nt) = var_block(ix,k)
       enddo
     enddo
   end subroutine copy_from_GFS_Data_3d_phys2phys

   pure subroutine copy_from_GFS_Data_2d_int2phys(ii1,jj1,isc,jsc,nt,var2d,var_block)
     implicit none
     integer, intent(in) :: ii1(:), jj1(:), isc, jsc, var_block(:)
     integer, intent(inout) :: nt
     real(kind=kind_phys), intent(out) :: var2d(:,:,:)
     integer ix

     nt=nt+1
     do ix=1,size(var_block)
       var2d(ii1(ix),jj1(ix),nt) = var_block(ix)
     enddo
   end subroutine copy_from_GFS_Data_2d_int2phys

   pure subroutine copy_from_GFS_Data_2d_stack_phys2phys(ii1,jj1,isc,jsc,nt,var3d,var_block)
     ! For copying phy_f2d and phy_fctd
     implicit none
     integer, intent(in) :: ii1(:), jj1(:), isc, jsc
     integer, intent(inout) :: nt
     real(kind=kind_phys), intent(in) :: var_block(:,:)
     real(kind=kind_phys), intent(out) :: var3d(:,:,:)
     integer ix, k

     nt=nt+1
     do k=lbound(var_block,2),ubound(var_block,2)
       do ix=1,size(var_block,1)
         var3d(ii1(ix),jj1(ix),nt) = var_block(ix,k)
       enddo
     enddo
   end subroutine copy_from_GFS_Data_2d_stack_phys2phys

   pure subroutine copy_from_GFS_Data_3d_int2phys(ii1,jj1,isc,jsc,nt,var3d,var_block)
     implicit none
     integer, intent(in) :: ii1(:), jj1(:), var_block(:,:), isc, jsc
     integer, intent(inout) :: nt
     real(kind=kind_phys), intent(out) :: var3d(:,:,:,:)
     integer ix, k

     nt=nt+1
     do k=lbound(var_block,2),ubound(var_block,2)
       do ix=1,size(var_block,1)
         var3d(ii1(ix),jj1(ix),k,nt) = real(var_block(ix,k),kind_phys)
       enddo
     enddo
   end subroutine copy_from_GFS_Data_3d_int2phys

   pure subroutine copy_to_GFS_Data_2d_phys2phys(ii1,jj1,isc,jsc,nt,var2d,var_block)
     implicit none
     integer, intent(in) :: ii1(:), jj1(:), isc, jsc
     integer, intent(inout) :: nt
     real(kind=kind_phys), intent(out) :: var_block(:)
     real(kind=kind_phys), intent(in) :: var2d(:,:,:)
     integer ix

     nt=nt+1
     do ix=1,size(var_block)
       var_block(ix) = var2d(ii1(ix),jj1(ix),nt)
     enddo
   end subroutine copy_to_GFS_Data_2d_phys2phys

   pure subroutine copy_to_GFS_Data_3d_phys2phys(ii1,jj1,isc,jsc,nt,var3d,var_block)
     implicit none
     integer, intent(in) :: ii1(:), jj1(:), isc, jsc
     integer, intent(inout) :: nt
     real(kind=kind_phys), intent(out) :: var_block(:,:)
     real(kind=kind_phys), intent(in) :: var3d(:,:,:,:)
     integer ix, k

     nt=nt+1
     do k=lbound(var_block,2),ubound(var_block,2)
       do ix=1,size(var_block,1)
         var_block(ix,k) = var3d(ii1(ix),jj1(ix),k,nt)
       enddo
     enddo
   end subroutine copy_to_GFS_Data_3d_phys2phys

   pure subroutine copy_to_GFS_Data_3d_slice_phys2phys(ii1,jj1,isc,jsc,nt,k1,k2,var3d,var_block)
     implicit none
     integer, intent(in) :: ii1(:), jj1(:), isc, jsc, k1, k2
     integer, intent(inout) :: nt
     real(kind=kind_phys), intent(out) :: var_block(:,:)
     real(kind=kind_phys), intent(in) :: var3d(:,:,:,:)
     integer ix, k

     nt=nt+1
     do k=k1,k2
       do ix=1,size(var_block,1)
         var_block(ix,k) = var3d(ii1(ix),jj1(ix),k,nt)
       enddo
     enddo
   end subroutine copy_to_GFS_Data_3d_slice_phys2phys

   pure subroutine copy_to_GFS_Data_2d_int2phys(ii1,jj1,isc,jsc,nt,var2d,var_block)
     implicit none
     integer, intent(in) :: ii1(:), jj1(:), isc, jsc
     integer, intent(inout) :: nt
     integer, intent(out) :: var_block(:)
     real(kind=kind_phys), intent(in) :: var2d(:,:,:)
     integer ix

     nt=nt+1
     do ix=1,size(var_block)
       var_block(ix) = int(var2d(ii1(ix),jj1(ix),nt))
     enddo
   end subroutine copy_to_GFS_Data_2d_int2phys

   pure subroutine copy_to_GFS_Data_3d_int2phys(ii1,jj1,isc,jsc,nt,var3d,var_block)
     implicit none
     integer, intent(in) :: ii1(:), jj1(:), isc, jsc
     integer, intent(inout) :: nt
     integer, intent(out) :: var_block(:,:)
     real(kind=kind_phys), intent(in) :: var3d(:,:,:,:)
     integer ix

     nt=nt+1
     do ix=1,size(var_block,1)
       var_block(ix,:) = int(var3d(ii1(ix),jj1(ix),:,nt))
     enddo
   end subroutine copy_to_GFS_Data_3d_int2phys

 subroutine create_2d_field_and_add_to_bundle(temp_r2d, field_name, outputfile, grid, bundle)

   use esmf

   implicit none

   real(kind_phys), dimension(:,:),   pointer, intent(in)    :: temp_r2d
   character(len=*),                           intent(in)    :: field_name
   character(len=*),                           intent(in)    :: outputfile
   type(ESMF_Grid),                            intent(in)    :: grid
   type(ESMF_FieldBundle),                     intent(inout) :: bundle

   type(ESMF_Field) :: field

   integer :: rc, i

   field = ESMF_FieldCreate(grid, temp_r2d, datacopyflag=ESMF_DATACOPY_REFERENCE, &
                            name=trim(field_name), indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, file=__FILE__)) &
   call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", name='output_file', value=trim(outputfile), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_FieldBundleAdd(bundle, (/field/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  end subroutine create_2d_field_and_add_to_bundle

 subroutine create_3d_field_and_add_to_bundle(temp_r3d, field_name, axis_name, num_levels, outputfile, grid, bundle)

   use esmf

   implicit none

   real(kind_phys), dimension(:,:,:), pointer, intent(in)    :: temp_r3d
   character(len=*),                           intent(in)    :: field_name
   character(len=*),                           intent(in)    :: axis_name
   integer,                                    intent(in)    :: num_levels
   character(len=*),                           intent(in)    :: outputfile
   type(ESMF_Grid),                            intent(in)    :: grid
   type(ESMF_FieldBundle),                     intent(inout) :: bundle

   type(ESMF_Field) :: field

   integer :: rc, i

   field = ESMF_FieldCreate(grid, temp_r3d, datacopyflag=ESMF_DATACOPY_REFERENCE, &
                            name=trim(field_name), indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, file=__FILE__)) &
   call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", name='output_file', value=trim(outputfile), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call add_zaxis_to_field(field, axis_name, num_levels)

   call ESMF_FieldBundleAdd(bundle, (/field/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  end subroutine create_3d_field_and_add_to_bundle

 subroutine add_zaxis_to_field(field, axis_name, num_levels)

   use esmf

   implicit none

   type(ESMF_Field), intent(inout) :: field
   character(len=*), intent(in)    :: axis_name
   integer,          intent(in)    :: num_levels

   real(kind_phys), allocatable, dimension(:) :: buffer
   integer :: rc, i

   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
                          name="ESMF:ungridded_dim_labels", valueList=(/trim(axis_name)/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

   allocate( buffer(num_levels) )
   do i=1, num_levels
      buffer(i)=i
   end do
   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3-dim", &
                          name=trim(axis_name), valueList=buffer, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
   deallocate(buffer)

   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3-dim", &
                          name=trim(axis_name)//"cartesian_axis", value="Z", rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

 end subroutine add_zaxis_to_field

end module FV3GFS_common_io
