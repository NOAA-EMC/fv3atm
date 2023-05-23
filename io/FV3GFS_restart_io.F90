module FV3GFS_restart_io_mod

  use esmf
  use block_control_mod,  only: block_control_type
  use GFS_typedefs,       only: GFS_sfcprop_type, GFS_control_type, kind_phys
  use GFS_restart,        only: GFS_restart_type
  use FV3GFS_sfc_io
  use FV3GFS_common_io,   only: create_2d_field_and_add_to_bundle, &
       create_3d_field_and_add_to_bundle, add_zaxis_to_field
  use FV3GFS_rrfs_sd_io
  use FV3GFS_clm_lake_io

  implicit none
  private

  real(kind_phys), parameter:: zero = 0.0, one = 1.0
  integer, parameter :: r8 = kind_phys

  integer :: nvar2d, nvar3d, npz
  real(kind=kind_phys), allocatable, target, dimension(:,:,:)   :: phy_var2
  real(kind=kind_phys), allocatable, target, dimension(:,:,:,:) :: phy_var3
  character(len=32),dimension(:),allocatable :: phy_var2_names, phy_var3_names

  type(rrfs_sd_state_type) :: rrfs_sd
  type(clm_lake_data_type) :: clm_lake
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

   if(Model%iopt_lake == 2 .and. Model%lkm > 0) then
     call clm_lake%allocate_data(Model)
   endif

   if(Model%rrfs_sd) then
     call rrfs_sd%allocate_data(Model)
   endif

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
    if(Model%iopt_lake == 2 .and. Model%lkm > 0) then
      call clm_lake%copy_from_grid(Model, Atm_block, Sfcprop)
    endif
    if(Model%rrfs_sd) then
      call rrfs_sd%copy_from_grid(Model, Atm_block, Sfcprop)
    endif

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

   sfc_bundle = bundle

   call ESMF_FieldBundleGet(bundle, name=sfcbdl_name,rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

   outputfile = trim(sfcbdl_name)

!*** add esmf fields

   call sfc%bundle_2d_fields(bundle, grid, Model, outputfile)
   call sfc%bundle_3d_fields(bundle, grid, Model, outputfile)

   if(Model%iopt_lake == 2 .and. Model%lkm > 0) then
     call clm_lake%bundle_fields(bundle, grid, Model, outputfile)
   endif
   if(Model%rrfs_sd) then
     call rrfs_sd%bundle_fields(bundle, grid, Model, outputfile)
   endif

 end subroutine fv_sfc_restart_bundle_setup

end module FV3GFS_restart_io_mod
