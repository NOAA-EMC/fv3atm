module FV3GFS_restart_io_mod

  use esmf
  use block_control_mod,  only: block_control_type
  use GFS_typedefs,       only: GFS_sfcprop_type, GFS_control_type, kind_phys
  use GFS_restart,        only: GFS_restart_type

  implicit none
  private

  real(kind_phys), parameter:: zero = 0.0, one = 1.0
  integer, parameter :: r8 = kind_phys

  integer :: nvar2d, nvar3d, npz
  real(kind=kind_phys), allocatable, target, dimension(:,:,:)   :: phy_var2
  real(kind=kind_phys), allocatable, target, dimension(:,:,:,:) :: phy_var3
  character(len=32),dimension(:),allocatable :: phy_var2_names, phy_var3_names

  integer :: nvar2m, nvar2o, nvar3, nvar2r, nvar2mp, nvar3mp
  real(kind=kind_phys), allocatable, target, dimension(:,:,:)   :: sfc_var2
  real(kind=kind_phys), allocatable, target, dimension(:,:,:)   :: sfc_var3ice
  real(kind=kind_phys), allocatable, target, dimension(:,:,:,:) :: sfc_var3, sfc_var3sn,sfc_var3eq,sfc_var3zn
  character(len=32),allocatable,dimension(:) :: sfc_name2, sfc_name3

  type(ESMF_FieldBundle) :: phy_bundle, sfc_bundle

  public FV3GFS_restart_register

  public fv_phy_restart_output
  public fv_phy_restart_bundle_setup

  public fv_sfc_restart_output
  public fv_sfc_restart_bundle_setup

  interface copy_from_GFS_Data
    module procedure copy_from_GFS_Data_2d_phys2phys, &
         copy_from_GFS_Data_3d_phys2phys, &
         copy_from_GFS_Data_2d_int2phys, &
         copy_from_GFS_Data_3d_int2phys, &
         copy_from_GFS_Data_2d_stack_phys2phys
  end interface

 contains

 subroutine FV3GFS_restart_register (Sfcprop, GFS_restart, Atm_block, Model)

   ! this subroutine must allocate all data buffers and set the variable names
   ! for both 'phy' and 'sfc' restart bundles

   implicit none

   type(GFS_sfcprop_type),      intent(in) :: Sfcprop(:)
   type(GFS_restart_type),      intent(in) :: GFS_Restart
   type(block_control_type),    intent(in) :: Atm_block
   type(GFS_control_type),      intent(in) :: Model

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
   nvar2m = 48
   if (Model%use_cice_alb .or. Model%lsm == Model%lsm_ruc) then
     nvar2m = nvar2m + 4
!    nvar2m = nvar2m + 5
   endif
   if (Model%cplwav) nvar2m = nvar2m + 1
   nvar2o = 18
   if (Model%lsm == Model%lsm_ruc) then
     if (Model%rdlai) then
       nvar2r = 13
     else
       nvar2r = 12
     endif
     nvar3  = 5
    else
     nvar2r = 0
     nvar3  = 3
   endif
   nvar2mp = 0
   nvar3mp = 0
   if (Model%lsm == Model%lsm_noahmp) then
     nvar2mp = 29
     nvar3mp = 5
   endif

   !--- allocate the various containers needed for restarts
   allocate(sfc_name2(nvar2m+nvar2o+nvar2mp+nvar2r))
   allocate(sfc_var2(nx,ny,nvar2m+nvar2o+nvar2mp+nvar2r))
   allocate(sfc_name3(0:nvar3+nvar3mp))
   allocate(sfc_var3ice(nx,ny,Model%kice))

   if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp) then
     allocate(sfc_var3(nx,ny,Model%lsoil,nvar3))
   elseif (Model%lsm == Model%lsm_ruc) then
     allocate(sfc_var3(nx,ny,Model%lsoil_lsm,nvar3))
   endif

   sfc_var2   = -9999.0_r8
   sfc_var3   = -9999.0_r8
   sfc_var3ice= -9999.0_r8

   if (Model%lsm == Model%lsm_noahmp) then
     allocate(sfc_var3sn(nx,ny,-2:0,4:6))
     allocate(sfc_var3eq(nx,ny,1:4,7:7))
     allocate(sfc_var3zn(nx,ny,-2:4,8:8))

     sfc_var3sn = -9999.0_r8
     sfc_var3eq = -9999.0_r8
     sfc_var3zn = -9999.0_r8
   endif

   call fill_Sfcprop_names(Model,sfc_name2,sfc_name3,nvar2m,.true.)

   sfc_name3(0) = 'tiice'

   if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp) then
      sfc_name3(1) = 'stc'
      sfc_name3(2) = 'smc'
      sfc_name3(3) = 'slc'
      if (Model%lsm == Model%lsm_noahmp) then
         sfc_name3(4) = 'snicexy'
         sfc_name3(5) = 'snliqxy'
         sfc_name3(6) = 'tsnoxy'
         sfc_name3(7) = 'smoiseq'
         sfc_name3(8) = 'zsnsoxy'
      endif
   else if (Model%lsm == Model%lsm_ruc) then
      sfc_name3(1) = 'tslb'
      sfc_name3(2) = 'smois'
      sfc_name3(3) = 'sh2o'
      sfc_name3(4) = 'smfr'
      sfc_name3(5) = 'flfr'
   end if

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

    integer :: i, j, k, nb, ix, lsoil, num, nt
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer, allocatable :: ii1(:), jj1(:)
    real(kind_phys) :: ice
    integer :: is, ie

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)

!$omp parallel do default(shared) private(i, j, nb, ix, nt, ii1, jj1, lsoil, k, ice)
    block_loop: do nb = 1, Atm_block%nblks
       allocate(ii1(Atm_block%blksz(nb)))
       allocate(jj1(Atm_block%blksz(nb)))
       ii1=Atm_block%index(nb)%ii - isc + 1
       jj1=Atm_block%index(nb)%jj - jsc + 1

       nt=0

       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%slmsk) !--- slmsk
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tsfco) !--- tsfc (tsea in sfc file)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%weasd) !--- weasd (sheleg in sfc file)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tg3)   !--- tg3
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%zorl)  !--- zorl
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%alvsf) !--- alvsf
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%alvwf) !--- alvwf
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%alnsf) !--- alnsf
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%alnwf) !--- alnwf
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%facsf) !--- facsf
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%facwf) !--- facwf
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%vfrac) !--- vfrac
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%canopy)!--- canopy
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%f10m)  !--- f10m
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%t2m)   !--- t2m
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%q2m)   !--- q2m

       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%vtype) !--- vtype
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%stype) !--- stype
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%uustar)!--- uustar
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%ffmm)  !--- ffmm
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%ffhh)  !--- ffhh
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%hice)  !--- hice
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%fice)  !--- fice
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tisfc) !--- tisfc
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tprcp) !--- tprcp
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%srflag)!--- srflag
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%snowd) !--- snowd (snwdph in the file)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%shdmin)!--- shdmin
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%shdmax)!--- shdmax
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%slope) !--- slope
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%snoalb)!--- snoalb
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%sncovr) !--- sncovr
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%snodl)  !--- snodl (snowd on land)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%weasdl) !--- weasdl (weasd on land)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tsfc)   !--- tsfc composite
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tsfcl)  !--- tsfcl (temp on land)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%zorlw)  !--- zorl (zorl on water)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%zorll)  !--- zorll (zorl on land)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%zorli)  !--- zorli (zorl on ice)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%albdirvis_lnd)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%albdirnir_lnd)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%albdifvis_lnd)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%albdifnir_lnd)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%emis_lnd)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%emis_ice)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%sncovr_ice)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%snodi)  !--- snodi (snowd on ice)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%weasdi) !--- weasdi (weasd on ice)
       if (Model%use_cice_alb .or. Model%lsm == Model%lsm_ruc) then
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%albdirvis_ice)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%albdifvis_ice)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%albdirnir_ice)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%albdifnir_ice)
!        sfc_var2(i,j,53) = Sfcprop(nb)%sfalb_ice(ix)
       endif
       if (Model%cplwav) then
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%zorlwav) !--- zorlwav (zorl from wav)
       endif
       !--- NSSTM variables
       if (Model%nstf_name(1) > 0) then
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tref)   !--- nsstm tref
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%z_c)    !--- nsstm z_c
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%c_0)    !--- nsstm c_0
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%c_d)    !--- nsstm c_d
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%w_0)    !--- nsstm w_0
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%w_d)    !--- nsstm w_d
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xt)     !--- nsstm xt
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xs)     !--- nsstm xs
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xu)     !--- nsstm xu
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xv)     !--- nsstm xv
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xz)     !--- nsstm xz
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%zm)     !--- nsstm zm
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xtts)   !--- nsstm xtts
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xzts)   !--- nsstm xzts
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%d_conv) !--- nsstm d_conv
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%ifd)    !--- nsstm ifd
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%dt_cool)!--- nsstm dt_cool
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%qrain)  !--- nsstm qrain

         ! FIXME convert negative zero (-0.0) to zero (0.0)
         do j=1,ny
         do i=1,nx
            if(sfc_var2(i,j,nt) == 0.0) sfc_var2(i,j,nt) = 0.0
         end do
         end do
       endif

       if (Model%lsm == Model%lsm_ruc) then
         !--- Extra RUC variables
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%wetness)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%clw_surf_land)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%clw_surf_ice)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%qwv_surf_land)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%qwv_surf_ice)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tsnow_land)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tsnow_ice)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%snowfallac_land)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%snowfallac_ice)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%sfalb_lnd)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%sfalb_lnd_bck)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%sfalb_ice)
         if (Model%rdlai) then
           call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xlaixy)
         endif
       else if (Model%lsm == Model%lsm_noahmp) then
         !--- Extra Noah MP variables
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%snowxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tvxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tgxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%canicexy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%canliqxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%eahxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tahxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%cmxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%chxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%fwetxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%sneqvoxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%alboldxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%qsnowxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%wslakexy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%zwtxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%waxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%wtxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%lfmassxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%rtmassxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%stmassxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%woodxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%stblcpxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%fastcpxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xsaixy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xlaixy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%taussxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%smcwtdxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%deeprechxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%rechxy)
       endif

       do k = 1,Model%kice
         do ix = 1, Atm_block%blksz(nb)
           ice=Sfcprop(nb)%tiice(ix,k)
           if(ice<one) then
             sfc_var3ice(ii1(ix),jj1(ix),k) = zero
           else
             sfc_var3ice(ii1(ix),jj1(ix),k) = ice
           endif
         enddo
       enddo

       if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp) then
          !--- 3D variables
          nt=0
          call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var3,Sfcprop(nb)%stc)
          call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var3,Sfcprop(nb)%smc)
          call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var3,Sfcprop(nb)%slc)
! 5 Noah MP 3D
          if (Model%lsm == Model%lsm_noahmp) then

            ! These arrays use bizarre indexing, which does not pass
            ! through function calls in Fortran, so we use loops here:

            nt=nt+1
            do lsoil = -2,0
              do ix = 1, Atm_block%blksz(nb)
                sfc_var3sn(ii1(ix),jj1(ix),lsoil,nt) = Sfcprop(nb)%snicexy(ix,lsoil)
              enddo
            enddo

            nt=nt+1
            do lsoil = -2,0
              do ix = 1, Atm_block%blksz(nb)
                sfc_var3sn(ii1(ix),jj1(ix),lsoil,nt) = Sfcprop(nb)%snliqxy(ix,lsoil)
              enddo
            enddo

            nt=nt+1
            do lsoil = -2,0
              do ix = 1, Atm_block%blksz(nb)
                sfc_var3sn(ii1(ix),jj1(ix),lsoil,nt) = Sfcprop(nb)%tsnoxy(ix,lsoil)
              enddo
            enddo

            nt=nt+1
            do lsoil = 1,Model%lsoil
              do ix = 1, Atm_block%blksz(nb)
                sfc_var3eq(ii1(ix),jj1(ix),lsoil,nt)  = Sfcprop(nb)%smoiseq(ix,lsoil)
              enddo
            enddo

            nt=nt+1
            do lsoil = -2,4
              do ix = 1, Atm_block%blksz(nb)
                sfc_var3zn(ii1(ix),jj1(ix),lsoil,nt)  = Sfcprop(nb)%zsnsoxy(ix,lsoil)
              enddo
            enddo

          endif  ! Noah MP
        else if (Model%lsm == Model%lsm_ruc) then
          !--- 3D variables
          nt=0
          call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var3,Sfcprop(nb)%tslb)
          call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var3,Sfcprop(nb)%smois)
          call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var3,Sfcprop(nb)%sh2o)
          call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var3,Sfcprop(nb)%keepsmfr)
          call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var3,Sfcprop(nb)%flag_frsoil)
        end if

       deallocate(ii1,jj1)
    enddo block_loop

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

   if (.not. allocated(sfc_var2)) then
     write(0,*)'ERROR sfc_var2, NOT allocated'
   endif
   if (.not. allocated(sfc_var3)) then
     write(0,*)'ERROR sfc_var3 NOT allocated'
   endif

   sfc_bundle = bundle

   call ESMF_FieldBundleGet(bundle, name=sfcbdl_name,rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

   outputfile = trim(sfcbdl_name)

!*** add esmf fields

   do num = 1,nvar2m
      temp_r2d => sfc_var2(:,:,num)
      call create_2d_field_and_add_to_bundle(temp_r2d, trim(sfc_name2(num)), outputfile, grid, bundle)
   enddo

   if (Model%nstf_name(1) > 0) then
      do num = nvar2m+1,nvar2m+nvar2o
         temp_r2d => sfc_var2(:,:,num)
         call create_2d_field_and_add_to_bundle(temp_r2d, trim(sfc_name2(num)), outputfile, grid, bundle)
      enddo
   endif

   if (Model%lsm == Model%lsm_ruc) then ! nvar2mp =0
      do num = nvar2m+nvar2o+1, nvar2m+nvar2o+nvar2r
         temp_r2d => sfc_var2(:,:,num)
         call create_2d_field_and_add_to_bundle(temp_r2d, trim(sfc_name2(num)), outputfile, grid, bundle)
      enddo
   else if (Model%lsm == Model%lsm_noahmp) then ! nvar2r =0
      do num = nvar2m+nvar2o+1,nvar2m+nvar2o+nvar2mp
         temp_r2d => sfc_var2(:,:,num)
         call create_2d_field_and_add_to_bundle(temp_r2d, trim(sfc_name2(num)), outputfile, grid, bundle)
      enddo
   endif

   temp_r3d => sfc_var3ice(:,:,:)
   call create_3d_field_and_add_to_bundle(temp_r3d, trim(sfc_name3(0)), "zaxis_1", Model%kice, trim(outputfile), grid, bundle)

   if(Model%lsm == Model%lsm_ruc) then
      do num = 1,nvar3
         temp_r3d => sfc_var3(:,:,:,num)
         call create_3d_field_and_add_to_bundle(temp_r3d, trim(sfc_name3(num)), "zaxis_1", Model%kice, trim(outputfile), grid, bundle)
      enddo
   else
      do num = 1,nvar3
         temp_r3d => sfc_var3(:,:,:,num)
         call create_3d_field_and_add_to_bundle(temp_r3d, trim(sfc_name3(num)), "zaxis_2", Model%lsoil, trim(outputfile), grid, bundle)
      enddo
   endif

   if (Model%lsm == Model%lsm_noahmp) then
      do num = nvar3+1,nvar3+3
         temp_r3d => sfc_var3sn(:,:,:,num)
         call create_3d_field_and_add_to_bundle(temp_r3d, trim(sfc_name3(num)), "zaxis_3", 3, trim(outputfile), grid, bundle)
      enddo

      temp_r3d => sfc_var3eq(:,:,:,7)
      call create_3d_field_and_add_to_bundle(temp_r3d, trim(sfc_name3(7)), "zaxis_2", Model%lsoil, trim(outputfile), grid, bundle)

      temp_r3d => sfc_var3zn(:,:,:,8)
      call create_3d_field_and_add_to_bundle(temp_r3d, trim(sfc_name3(8)), "zaxis_4", 7, trim(outputfile), grid, bundle)
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

   pure subroutine fill_Sfcprop_names(Model,sfc_name2,sfc_name3,nvar_s2m,warm_start)
     implicit none
     type(GFS_control_type),    intent(in) :: Model
     integer, intent(in) :: nvar_s2m
     character(len=32),intent(out) :: sfc_name2(:), sfc_name3(:)
     logical, intent(in) :: warm_start
     integer :: nt

      !--- names of the 2D variables to save
      nt=0
      nt=nt+1 ; sfc_name2(nt) = 'slmsk'
      nt=nt+1 ; sfc_name2(nt) = 'tsea'    !tsfc
      nt=nt+1 ; sfc_name2(nt) = 'sheleg'  !weasd
      nt=nt+1 ; sfc_name2(nt) = 'tg3'
      nt=nt+1 ; sfc_name2(nt) = 'zorl'
      nt=nt+1 ; sfc_name2(nt) = 'alvsf'
      nt=nt+1 ; sfc_name2(nt) = 'alvwf'
      nt=nt+1 ; sfc_name2(nt) = 'alnsf'
      nt=nt+1 ; sfc_name2(nt) = 'alnwf'
      nt=nt+1 ; sfc_name2(nt) = 'facsf'
      nt=nt+1 ; sfc_name2(nt) = 'facwf'
      nt=nt+1 ; sfc_name2(nt) = 'vfrac'
      nt=nt+1 ; sfc_name2(nt) = 'canopy'
      nt=nt+1 ; sfc_name2(nt) = 'f10m'
      nt=nt+1 ; sfc_name2(nt) = 't2m'
      nt=nt+1 ; sfc_name2(nt) = 'q2m'
      nt=nt+1 ; sfc_name2(nt) = 'vtype'
      nt=nt+1 ; sfc_name2(nt) = 'stype'
      nt=nt+1 ; sfc_name2(nt) = 'uustar'
      nt=nt+1 ; sfc_name2(nt) = 'ffmm'
      nt=nt+1 ; sfc_name2(nt) = 'ffhh'
      nt=nt+1 ; sfc_name2(nt) = 'hice'
      nt=nt+1 ; sfc_name2(nt) = 'fice'
      nt=nt+1 ; sfc_name2(nt) = 'tisfc'
      nt=nt+1 ; sfc_name2(nt) = 'tprcp'
      nt=nt+1 ; sfc_name2(nt) = 'srflag'
      nt=nt+1 ; sfc_name2(nt) = 'snwdph'  !snowd
      nt=nt+1 ; sfc_name2(nt) = 'shdmin'
      nt=nt+1 ; sfc_name2(nt) = 'shdmax'
      nt=nt+1 ; sfc_name2(nt) = 'slope'
      nt=nt+1 ; sfc_name2(nt) = 'snoalb'
      !--- variables below here are optional
      nt=nt+1 ; sfc_name2(nt) = 'sncovr'
      nt=nt+1 ; sfc_name2(nt) = 'snodl' !snowd on land portion of a cell
      nt=nt+1 ; sfc_name2(nt) = 'weasdl'!weasd on land portion of a cell
      nt=nt+1 ; sfc_name2(nt) = 'tsfc'  !tsfc composite
      nt=nt+1 ; sfc_name2(nt) = 'tsfcl' !temp on land portion of a cell
      nt=nt+1 ; sfc_name2(nt) = 'zorlw' !zorl on water portion of a cell
      nt=nt+1 ; sfc_name2(nt) = 'zorll' !zorl on land portion of a cell
      nt=nt+1 ; sfc_name2(nt) = 'zorli' !zorl on ice portion of a cell
      nt=nt+1 ; sfc_name2(nt) = 'albdirvis_lnd'
      nt=nt+1 ; sfc_name2(nt) = 'albdirnir_lnd'
      nt=nt+1 ; sfc_name2(nt) = 'albdifvis_lnd'
      nt=nt+1 ; sfc_name2(nt) = 'albdifnir_lnd'
      nt=nt+1 ; sfc_name2(nt) = 'emis_lnd'
      nt=nt+1 ; sfc_name2(nt) = 'emis_ice'
      nt=nt+1 ; sfc_name2(nt) = 'sncovr_ice'
      nt=nt+1 ; sfc_name2(nt) = 'snodi' ! snowd on ice portion of a cell
      nt=nt+1 ; sfc_name2(nt) = 'weasdi'! weasd on ice portion of a cell

      if (Model%use_cice_alb .or. Model%lsm == Model%lsm_ruc) then
        nt=nt+1 ; sfc_name2(nt) = 'albdirvis_ice'
        nt=nt+1 ; sfc_name2(nt) = 'albdifvis_ice'
        nt=nt+1 ; sfc_name2(nt) = 'albdirnir_ice'
        nt=nt+1 ; sfc_name2(nt) = 'albdifnir_ice'
!        nt=nt+1 ; sfc_name2(nt) = 'sfalb_ice'
      endif

      if(Model%cplwav) then
        sfc_name2(nvar_s2m) = 'zorlwav' !zorl from wave component
      endif

      nt = nvar_s2m ! next variable will be at nvar_s2m

      !--- NSSTM inputs only needed when (nstf_name(1) > 0) .and. (nstf_name(2)) == 0)
      nt=nt+1 ; sfc_name2(nt) = 'tref'
      nt=nt+1 ; sfc_name2(nt) = 'z_c'
      nt=nt+1 ; sfc_name2(nt) = 'c_0'
      nt=nt+1 ; sfc_name2(nt) = 'c_d'
      nt=nt+1 ; sfc_name2(nt) = 'w_0'
      nt=nt+1 ; sfc_name2(nt) = 'w_d'
      nt=nt+1 ; sfc_name2(nt) = 'xt'
      nt=nt+1 ; sfc_name2(nt) = 'xs'
      nt=nt+1 ; sfc_name2(nt) = 'xu'
      nt=nt+1 ; sfc_name2(nt) = 'xv'
      nt=nt+1 ; sfc_name2(nt) = 'xz'
      nt=nt+1 ; sfc_name2(nt) = 'zm'
      nt=nt+1 ; sfc_name2(nt) = 'xtts'
      nt=nt+1 ; sfc_name2(nt) = 'xzts'
      nt=nt+1 ; sfc_name2(nt) = 'd_conv'
      nt=nt+1 ; sfc_name2(nt) = 'ifd'
      nt=nt+1 ; sfc_name2(nt) = 'dt_cool'
      nt=nt+1 ; sfc_name2(nt) = 'qrain'
!
! Only needed when Noah MP LSM is used - 29 2D
!
      if (Model%lsm == Model%lsm_noahmp) then
        nt=nt+1 ; sfc_name2(nt) = 'snowxy'
        nt=nt+1 ; sfc_name2(nt) = 'tvxy'
        nt=nt+1 ; sfc_name2(nt) = 'tgxy'
        nt=nt+1 ; sfc_name2(nt) = 'canicexy'
        nt=nt+1 ; sfc_name2(nt) = 'canliqxy'
        nt=nt+1 ; sfc_name2(nt) = 'eahxy'
        nt=nt+1 ; sfc_name2(nt) = 'tahxy'
        nt=nt+1 ; sfc_name2(nt) = 'cmxy'
        nt=nt+1 ; sfc_name2(nt) = 'chxy'
        nt=nt+1 ; sfc_name2(nt) = 'fwetxy'
        nt=nt+1 ; sfc_name2(nt) = 'sneqvoxy'
        nt=nt+1 ; sfc_name2(nt) = 'alboldxy'
        nt=nt+1 ; sfc_name2(nt) = 'qsnowxy'
        nt=nt+1 ; sfc_name2(nt) = 'wslakexy'
        nt=nt+1 ; sfc_name2(nt) = 'zwtxy'
        nt=nt+1 ; sfc_name2(nt) = 'waxy'
        nt=nt+1 ; sfc_name2(nt) = 'wtxy'
        nt=nt+1 ; sfc_name2(nt) = 'lfmassxy'
        nt=nt+1 ; sfc_name2(nt) = 'rtmassxy'
        nt=nt+1 ; sfc_name2(nt) = 'stmassxy'
        nt=nt+1 ; sfc_name2(nt) = 'woodxy'
        nt=nt+1 ; sfc_name2(nt) = 'stblcpxy'
        nt=nt+1 ; sfc_name2(nt) = 'fastcpxy'
        nt=nt+1 ; sfc_name2(nt) = 'xsaixy'
        nt=nt+1 ; sfc_name2(nt) = 'xlaixy'
        nt=nt+1 ; sfc_name2(nt) = 'taussxy'
        nt=nt+1 ; sfc_name2(nt) = 'smcwtdxy'
        nt=nt+1 ; sfc_name2(nt) = 'deeprechxy'
        nt=nt+1 ; sfc_name2(nt) = 'rechxy'
      else if (Model%lsm == Model%lsm_ruc .and. warm_start) then
        nt=nt+1 ; sfc_name2(nt) = 'wetness'
        nt=nt+1 ; sfc_name2(nt) = 'clw_surf_land'
        nt=nt+1 ; sfc_name2(nt) = 'clw_surf_ice'
        nt=nt+1 ; sfc_name2(nt) = 'qwv_surf_land'
        nt=nt+1 ; sfc_name2(nt) = 'qwv_surf_ice'
        nt=nt+1 ; sfc_name2(nt) = 'tsnow_land'
        nt=nt+1 ; sfc_name2(nt) = 'tsnow_ice'
        nt=nt+1 ; sfc_name2(nt) = 'snowfall_acc_land'
        nt=nt+1 ; sfc_name2(nt) = 'snowfall_acc_ice'
        nt=nt+1 ; sfc_name2(nt) = 'sfalb_lnd'
        nt=nt+1 ; sfc_name2(nt) = 'sfalb_lnd_bck'
        nt=nt+1 ; sfc_name2(nt) = 'sfalb_ice'
        if (Model%rdlai) then
          nt=nt+1 ; sfc_name2(nt) = 'lai'
        endif
      else if (Model%lsm == Model%lsm_ruc .and. Model%rdlai) then
        nt=nt+1 ; sfc_name2(nt) = 'lai'
      endif
   end subroutine fill_sfcprop_names

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

end module FV3GFS_restart_io_mod
