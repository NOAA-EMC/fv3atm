module FV3GFS_restart_io_mod

  use esmf
  use block_control_mod,  only: block_control_type
  use GFS_typedefs,       only: GFS_sfcprop_type, GFS_control_type, kind_phys
  use GFS_restart,        only: GFS_restart_type
  use FV3GFS_sfc_io
  use FV3GFS_common_io,   only: copy_from_GFS_Data

  implicit none
  private

  real(kind_phys), parameter:: zero = 0.0, one = 1.0
  integer, parameter :: r8 = kind_phys

  integer :: nvar2d, nvar3d, npz
  real(kind=kind_phys), allocatable, target, dimension(:,:,:)   :: phy_var2
  real(kind=kind_phys), allocatable, target, dimension(:,:,:,:) :: phy_var3
  character(len=32),dimension(:),allocatable :: phy_var2_names, phy_var3_names

  type(Sfc_io_data_type) :: sfc
  type(ESMF_FieldBundle) :: phy_bundle, sfc_bundle

  public FV3GFS_restart_register

  public fv_phy_restart_output
  public fv_phy_restart_bundle_setup

  public fv_sfc_restart_output
  public fv_sfc_restart_bundle_setup

 contains

 subroutine FV3GFS_restart_register (Sfcprop, GFS_restart, Atm_block, Model)

   ! this subroutine must allocate all data buffers and set the variable names
   ! for both 'phy' and 'sfc' restart bundles

   implicit none

   type(GFS_sfcprop_type),      intent(in) :: Sfcprop(:)
   type(GFS_restart_type),      intent(in) :: GFS_Restart
   type(block_control_type),    intent(in) :: Atm_block
   type(GFS_control_type),      intent(in) :: Model

   logical was_changed
   integer :: isc, iec, jsc, jec, nx, ny
   integer :: num

   isc = Atm_block%isc
   iec = Atm_block%iec
   jsc = Atm_block%jsc
   jec = Atm_block%jec
   npz = Atm_block%npz
   nx  = (iec - isc + 1)
   ny  = (jec - jsc + 1)

   !--------------- phy
   nvar2d = GFS_Restart%num2d
   nvar3d = GFS_Restart%num3d

   allocate (phy_var2(nx,ny,nvar2d), phy_var2_names(nvar2d))
   allocate (phy_var3(nx,ny,npz,nvar3d), phy_var3_names(nvar3d))
   phy_var2 = zero
   phy_var3 = zero
   do num = 1,nvar2d
     phy_var2_names(num) = trim(GFS_Restart%name2d(num))
   enddo
   do num = 1,nvar3d
     phy_var3_names(num) = trim(GFS_Restart%name3d(num))
   enddo

   !--------------- sfc
   was_changed = sfc%allocate_arrays(Model, Atm_block, .false., .true.)
   call sfc%fill_2d_names(Model, .true.)
   call sfc%fill_3d_names(Model, .true.)

 end subroutine FV3GFS_restart_register

 subroutine fv_phy_restart_output(GFS_Restart, Atm_block)

   implicit none

   type(GFS_restart_type),      intent(in) :: GFS_Restart
   type(block_control_type),    intent(in) :: Atm_block

!*** local variables
   integer :: i, j, k, n
   integer :: nb, ix, num
   integer :: isc, iec, jsc, jec, npz, nx, ny
   integer(8) :: rchk

   isc = Atm_block%isc
   iec = Atm_block%iec
   jsc = Atm_block%jsc
   jec = Atm_block%jec
   npz = Atm_block%npz
   nx  = (iec - isc + 1)
   ny  = (jec - jsc + 1)

   !--- register the restart fields
   if (.not. allocated(phy_var2)) then
      write(0,*)'phy_var2 must be allocated'
   endif
   if (.not. allocated(phy_var3)) then
      write(0,*)'phy_var3 must be allocated'
   endif

   !--- 2D variables
    do num = 1,nvar2d
      do nb = 1,Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          i = Atm_block%index(nb)%ii(ix) - isc + 1
          j = Atm_block%index(nb)%jj(ix) - jsc + 1
          phy_var2(i,j,num) = GFS_Restart%data(nb,num)%var2p(ix)
        enddo
      enddo
    enddo

    !--- 3D variables
    do num = 1,nvar3d
      do nb = 1,Atm_block%nblks
        do k=1,npz
          do ix = 1, Atm_block%blksz(nb)
            i = Atm_block%index(nb)%ii(ix) - isc + 1
            j = Atm_block%index(nb)%jj(ix) - jsc + 1
            phy_var3(i,j,k,num) = GFS_Restart%data(nb,num)%var3p(ix,k)
          enddo
        enddo
      enddo
    enddo

 end subroutine fv_phy_restart_output

 subroutine fv_sfc_restart_output(Sfcprop, Atm_block, Model)
    !--- interface variable definitions
   implicit none

    type(GFS_sfcprop_type),      intent(in) :: Sfcprop(:)
    type(block_control_type),    intent(in) :: Atm_block
    type(GFS_control_type),      intent(in) :: Model

    call sfc%copy_from_grid(Model, Atm_block, Sfcprop)

 end subroutine fv_sfc_restart_output

 subroutine fv_phy_restart_bundle_setup(bundle, grid, rc)
!
!-------------------------------------------------------------
!*** set esmf bundle for phys restart fields
!------------------------------------------------------------
!
   use esmf

   implicit none

   type(ESMF_FieldBundle),intent(inout)        :: bundle
   type(ESMF_Grid),intent(inout)               :: grid
   integer,intent(out)                         :: rc

!*** local variables
   integer i, j, k, n
   character(128)    :: bdl_name
   type(ESMF_Field)  :: field
   character(128)    :: outputfile
   real(kind_phys),dimension(:,:),pointer   :: temp_r2d
   real(kind_phys),dimension(:,:,:),pointer   :: temp_r3d
   integer :: num

   if (.not. allocated(phy_var2)) then
     write(0,*)'ERROR phy_var2, NOT allocated'
   endif
   if (.not. allocated(phy_var3)) then
     write(0,*)'ERROR phy_var3 NOT allocated'
   endif

   phy_bundle = bundle

   call ESMF_FieldBundleGet(bundle, name=bdl_name,rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

   outputfile = trim(bdl_name)

!*** add esmf fields

   do num = 1,nvar2d
       temp_r2d => phy_var2(:,:,num)
       call create_2d_field_and_add_to_bundle(temp_r2d, trim(phy_var2_names(num)), trim(outputfile), grid, bundle)
   enddo

   do num = 1,nvar3d
       temp_r3d => phy_var3(:,:,:,num)
       call create_3d_field_and_add_to_bundle(temp_r3d, trim(phy_var3_names(num)), "zaxis_1", npz, trim(outputfile), grid, bundle)
   enddo

 end subroutine fv_phy_restart_bundle_setup

 subroutine fv_sfc_restart_bundle_setup(bundle, grid, Model, rc)
!
!-------------------------------------------------------------
!*** set esmf bundle for sfc restart fields
!------------------------------------------------------------
!
   use esmf

   implicit none

   type(ESMF_FieldBundle),intent(inout)        :: bundle
   type(ESMF_Grid),intent(inout)               :: grid
   type(GFS_control_type),          intent(in) :: Model
   integer,intent(out)                         :: rc

!*** local variables
   integer i, j, k, n
   character(128)    :: sfcbdl_name
   type(ESMF_Field)  :: field
   character(128)    :: outputfile
   real(kind_phys),dimension(:,:),pointer   :: temp_r2d
   real(kind_phys),dimension(:,:,:),pointer :: temp_r3d

   integer :: num

   if (.not. associated(sfc%var2)) then
     write(0,*)'ERROR sfc%var2, NOT associated'
   endif
   if (.not. associated(sfc%var3)) then
     write(0,*)'ERROR sfc%var3 NOT associated'
   endif

   sfc_bundle = bundle

   call ESMF_FieldBundleGet(bundle, name=sfcbdl_name,rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

   outputfile = trim(sfcbdl_name)

!*** add esmf fields

   do num = 1,sfc%nvar2m
      temp_r2d => sfc%var2(:,:,num)
      call create_2d_field_and_add_to_bundle(temp_r2d, trim(sfc%name2(num)), outputfile, grid, bundle)
   enddo

   if (Model%nstf_name(1) > 0) then
      do num = sfc%nvar2m+1,sfc%nvar2m+sfc%nvar2o
         temp_r2d => sfc%var2(:,:,num)
         call create_2d_field_and_add_to_bundle(temp_r2d, trim(sfc%name2(num)), outputfile, grid, bundle)
      enddo
   endif

   if (Model%lsm == Model%lsm_ruc) then ! sfc%nvar2mp =0
      do num = sfc%nvar2m+sfc%nvar2o+1, sfc%nvar2m+sfc%nvar2o+sfc%nvar2r
         temp_r2d => sfc%var2(:,:,num)
         call create_2d_field_and_add_to_bundle(temp_r2d, trim(sfc%name2(num)), outputfile, grid, bundle)
      enddo
   else if (Model%lsm == Model%lsm_noahmp) then ! sfc%nvar2r =0
      do num = sfc%nvar2m+sfc%nvar2o+1,sfc%nvar2m+sfc%nvar2o+sfc%nvar2mp
         temp_r2d => sfc%var2(:,:,num)
         call create_2d_field_and_add_to_bundle(temp_r2d, trim(sfc%name2(num)), outputfile, grid, bundle)
      enddo
   endif

   temp_r3d => sfc%var3ice(:,:,:)
   call create_3d_field_and_add_to_bundle(temp_r3d, trim(sfc%name3(0)), "zaxis_1", Model%kice, trim(outputfile), grid, bundle)

   if(Model%lsm == Model%lsm_ruc) then
      do num = 1,sfc%nvar3
         temp_r3d => sfc%var3(:,:,:,num)
         call create_3d_field_and_add_to_bundle(temp_r3d, trim(sfc%name3(num)), "zaxis_1", Model%kice, trim(outputfile), grid, bundle)
      enddo
   else
      do num = 1,sfc%nvar3
         temp_r3d => sfc%var3(:,:,:,num)
         call create_3d_field_and_add_to_bundle(temp_r3d, trim(sfc%name3(num)), "zaxis_2", Model%lsoil, trim(outputfile), grid, bundle)
      enddo
   endif

   if (Model%lsm == Model%lsm_noahmp) then
      do num = sfc%nvar3+1,sfc%nvar3+3
         temp_r3d => sfc%var3sn(:,:,:,num)
         call create_3d_field_and_add_to_bundle(temp_r3d, trim(sfc%name3(num)), "zaxis_3", 3, trim(outputfile), grid, bundle)
      enddo

      temp_r3d => sfc%var3eq(:,:,:,7)
      call create_3d_field_and_add_to_bundle(temp_r3d, trim(sfc%name3(7)), "zaxis_2", Model%lsoil, trim(outputfile), grid, bundle)

      temp_r3d => sfc%var3zn(:,:,:,8)
      call create_3d_field_and_add_to_bundle(temp_r3d, trim(sfc%name3(8)), "zaxis_4", 7, trim(outputfile), grid, bundle)
   endif ! lsm = lsm_noahmp

 end subroutine fv_sfc_restart_bundle_setup

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

end module FV3GFS_restart_io_mod
