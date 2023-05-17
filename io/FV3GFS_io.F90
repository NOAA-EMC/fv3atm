module FV3GFS_io_mod

!-----------------------------------------------------------------------
!    gfs_physics_driver_mod defines the GFS physics routines used by
!    the GFDL FMS system to obtain tendencies and boundary fluxes due
!    to the physical parameterizations and processes that drive
!    atmospheric time tendencies for use by other components, namely
!    the atmospheric dynamical core.
!
!    NOTE: This module currently supports only the operational GFS
!          parameterizations as of September 2015.  Further development
!          is needed to support the full suite of physical
!          parameterizations present in the GFS physics package.
!-----------------------------------------------------------------------
!
!--- FMS/GFDL modules
  use block_control_mod,  only: block_control_type
  use mpp_mod,            only: mpp_error,  mpp_pe, mpp_root_pe, &
                                mpp_chksum, NOTE,   FATAL
  use fms_mod,            only: stdout
  use fms2_io_mod,        only: FmsNetcdfDomainFile_t, unlimited,      &
                                open_file, close_file,                 &
                                register_axis, register_restart_field, &
                                register_variable_attribute, register_field, &
                                read_restart, write_restart, write_data,     &
                                get_global_io_domain_indices, variable_exists
  use mpp_domains_mod,    only: domain1d, domain2d, domainUG
  use time_manager_mod,   only: time_type
  use diag_manager_mod,   only: register_diag_field, send_data
  use diag_axis_mod,      only: get_axis_global_length, get_diag_axis, &
                                get_diag_axis_name
  use diag_data_mod,      only: output_fields, max_output_fields
  use diag_util_mod,      only: find_input_field
  use constants_mod,      only: grav, rdgas
  use physcons,           only: con_tice          !saltwater freezing temp (K)
  use clm_lake_io,        only: clm_lake_data_type
!
!--- GFS_typedefs
  use GFS_typedefs,       only: GFS_sfcprop_type, GFS_control_type, &
                                GFS_data_type, kind_phys
  use GFS_restart,        only: GFS_restart_type
  use GFS_diagnostics,    only: GFS_externaldiag_type

  use FV3GFS_common_io,   only: copy_from_GFS_Data, copy_to_GFS_Data
  use FV3GFS_sfc_io
  use FV3GFS_oro_io

!
!-----------------------------------------------------------------------
  implicit none
  private

  !--- public interfaces ---
  public  FV3GFS_restart_read, FV3GFS_restart_write
  public  FV3GFS_GFS_checksum
  public  fv3gfs_diag_register, fv3gfs_diag_output
#ifdef use_WRTCOMP
  public  fv_phys_bundle_setup
#endif

  !--- GFDL filenames
  character(len=32)  :: fn_oro    = 'oro_data.nc'
  character(len=32)  :: fn_oro_ls = 'oro_data_ls.nc'
  character(len=32)  :: fn_oro_ss = 'oro_data_ss.nc'
  character(len=32)  :: fn_srf    = 'sfc_data.nc'
  character(len=32)  :: fn_phy    = 'phy_data.nc'
  character(len=32)  :: fn_dust12m= 'dust12m_data.nc'
  character(len=32)  :: fn_emi    = 'emi_data.nc'
  character(len=32)  :: fn_rrfssd = 'SMOKE_RRFS_data.nc'

  !--- GFDL FMS netcdf restart data types defined in fms2_io
  type(FmsNetcdfDomainFile_t) :: Oro_restart, Sfc_restart, Phy_restart, dust12m_restart, emi_restart, rrfssd_restart
  type(FmsNetcdfDomainFile_t) :: Oro_ls_restart, Oro_ss_restart

  !--- GFDL FMS restart containers
  real(kind=kind_phys), allocatable, target, dimension(:,:,:,:) :: phy_var3
  character(len=32),    allocatable,         dimension(:)       :: oro_ls_ss_name
  real(kind=kind_phys), allocatable, target, dimension(:,:,:)   :: oro_ls_var, oro_ss_var
  real(kind=kind_phys), allocatable, target, dimension(:,:,:)   :: phy_var2

  !--- Noah MP restart containers
  real(kind=kind_phys) :: zhour
!
  integer, parameter :: r8 = kind_phys
  integer :: tot_diag_idx = 0
  integer :: total_outputlevel = 0
  integer :: isco,ieco,jsco,jeco,levo,num_axes_phys
  integer :: fhzero, ncld, nsoil, imp_physics, landsfcmdl
  real(4) :: dtp
  logical :: lprecip_accu
  character(len=64)  :: Sprecip_accu
  integer,dimension(:),        allocatable         :: nstt, nstt_vctbl, all_axes
  character(20),dimension(:),  allocatable         :: axis_name, axis_name_vert
  real(4), dimension(:,:,:),   allocatable, target :: buffer_phys_bl, buffer_phys_nb
  real(4), dimension(:,:,:,:), allocatable, target :: buffer_phys_windvect
  real(kind=kind_phys),dimension(:,:),allocatable  :: lon, lat, uwork
  real(kind=kind_phys),dimension(:,:,:),allocatable:: uwork3d
  logical                    :: uwork_set = .false.
  character(128)             :: uwindname
  integer, parameter, public :: DIAG_SIZE = 800
  real, parameter :: missing_value = 9.99e20_r8
  real, parameter:: stndrd_atmos_ps = 101325.0_r8
  real, parameter:: stndrd_atmos_lapse = 0.0065_r8
  real, parameter:: drythresh = 1.e-4_r8, zero = 0.0_r8, one = 1.0_r8

  type(Sfc_io_data_type) :: sfc
  type(Oro_io_data_type) :: oro

!--- miscellaneous other variables
  logical :: use_wrtgridcomp_output = .FALSE.
  logical :: module_is_initialized  = .FALSE.

  CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!--------------------
! FV3GFS_restart_read
!--------------------
  subroutine FV3GFS_restart_read (GFS_Data, GFS_Restart, Atm_block, Model, fv_domain, warm_start, ignore_rst_cksum)
    type(GFS_data_type),      intent(inout) :: GFS_Data(:)
    type(GFS_restart_type),   intent(inout) :: GFS_Restart
    type(block_control_type), intent(in)    :: Atm_block
    type(GFS_control_type),   intent(inout) :: Model
    type(domain2d),           intent(in)    :: fv_domain
    logical,                  intent(in)    :: warm_start
    logical,                  intent(in)    :: ignore_rst_cksum

    !--- read in surface data from chgres
    call sfc_prop_restart_read (GFS_Data%Sfcprop, Atm_block, Model, fv_domain, warm_start, ignore_rst_cksum)

    !--- read in physics restart data
    call phys_restart_read (GFS_Restart, Atm_block, Model, fv_domain, ignore_rst_cksum)

  end subroutine FV3GFS_restart_read

!---------------------
! FV3GFS_restart_write
!---------------------
  subroutine FV3GFS_restart_write (GFS_Data, GFS_Restart, Atm_block, Model, fv_domain, timestamp)
    type(GFS_data_type),         intent(inout) :: GFS_Data(:)
    type(GFS_restart_type),      intent(inout) :: GFS_Restart
    type(block_control_type),    intent(in)    :: Atm_block
    type(GFS_control_type),      intent(in)    :: Model
    type(domain2d),              intent(in)    :: fv_domain
    character(len=32), optional, intent(in)    :: timestamp

    !--- write surface data from chgres
    call sfc_prop_restart_write (GFS_Data%Sfcprop, Atm_block, Model, fv_domain, timestamp)

    !--- write physics restart data
    call phys_restart_write (GFS_Restart, Atm_block, Model, fv_domain, timestamp)

  end subroutine FV3GFS_restart_write


!--------------------
! FV3GFS_GFS_checksum
!--------------------
 subroutine FV3GFS_GFS_checksum (Model, GFS_Data, Atm_block)
   !--- interface variables
   type(GFS_control_type),    intent(in) :: Model
   type(GFS_data_type),       intent(in) :: GFS_Data(:)
   type (block_control_type), intent(in) :: Atm_block
   !--- local variables
   integer :: outunit, j, i, ix, nb, isc, iec, jsc, jec, lev, ct, l, ntr, k
   integer :: nsfcprop2d, idx_opt, nt
   real(kind=kind_phys), allocatable :: temp2d(:,:,:)
   real(kind=kind_phys), allocatable :: temp3d(:,:,:,:)
   real(kind=kind_phys), allocatable :: temp3dlevsp1(:,:,:,:)
   integer, allocatable :: ii1(:), jj1(:)
   character(len=32) :: name

   isc = Model%isc
   iec = Model%isc+Model%nx-1
   jsc = Model%jsc
   jec = Model%jsc+Model%ny-1
   lev = Model%levs

   ntr = size(GFS_Data(1)%Statein%qgrs,3)

     nsfcprop2d = 93
   if (Model%lsm == Model%lsm_noahmp) then
     nsfcprop2d = nsfcprop2d + 49
     if (Model%use_cice_alb) then
       nsfcprop2d = nsfcprop2d + 4
     endif
   elseif (Model%lsm == Model%lsm_ruc) then
     nsfcprop2d = nsfcprop2d + 4 + 12
     if (Model%rdlai) then
       nsfcprop2d = nsfcprop2d + 1
     endif
   else
     if (Model%use_cice_alb) then
       nsfcprop2d = nsfcprop2d + 4
     endif
   endif

   if (Model%nstf_name(1) > 0) then
     nsfcprop2d = nsfcprop2d + 16
   endif

   if(Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_flake) then
     nsfcprop2d = nsfcprop2d + 10
   endif

   allocate (temp2d(isc:iec,jsc:jec,nsfcprop2d+Model%ntot2d+Model%nctp))
   allocate (temp3d(isc:iec,jsc:jec,1:lev,14+Model%ntot3d+2*ntr))
   allocate (temp3dlevsp1(isc:iec,jsc:jec,1:lev+1,3))

   temp2d = zero
   temp3d = zero
   temp3dlevsp1 = zero

!$omp parallel do default(shared) private(i, j, k, nb, ix, nt, ii1, jj1)
    block_loop: do nb = 1, Atm_block%nblks
       allocate(ii1(Atm_block%blksz(nb)))
       allocate(jj1(Atm_block%blksz(nb)))
       ii1=Atm_block%index(nb)%ii - isc + 1
       jj1=Atm_block%index(nb)%jj - jsc + 1

       ! Copy into temp2d
       nt=0

       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Statein%pgr)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%slmsk)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tsfc)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tisfc)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%snowd)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%zorl)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%fice)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%hprime(:,1))
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%sncovr)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%snoalb)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%alvsf)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%alnsf)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%alvwf)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%alnwf)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%facsf)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%facwf)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%slope)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%shdmin)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%shdmax)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tg3)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%vfrac)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%vtype)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%stype)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%uustar)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%oro)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%oro_uf)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%hice)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%weasd)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%canopy)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%ffmm)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%ffhh)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%f10m)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tprcp)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%srflag)
       lsm_choice: if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp) then
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%slc)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%smc)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%stc)
       elseif (Model%lsm == Model%lsm_ruc) then
         do k=1,3
           call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%sh2o(:,k))
         enddo
         ! Combine levels 4 to lsoil_lsm (9 for RUC) into one
         nt=nt+1
         do ix=1,Atm_block%blksz(nb)
           temp2d(ii1(ix),jj1(ix),nt) = sum(GFS_Data(nb)%Sfcprop%sh2o(ix,4:Model%lsoil_lsm))
         enddo
         do k=1,3
           call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%smois(:,k))
         enddo
         ! Combine levels 4 to lsoil_lsm (9 for RUC) into one
         nt=nt+1
         do ix=1,Atm_block%blksz(nb)
           temp2d(ii1(ix),jj1(ix),nt) = sum(GFS_Data(nb)%Sfcprop%smois(ix,4:Model%lsoil_lsm))
         enddo
         do k=1,3
           call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tslb(:,k))
         enddo
         ! Combine levels 4 to lsoil_lsm (9 for RUC) into one
         nt=nt+1
         do ix=1,Atm_block%blksz(nb)
           temp2d(ii1(ix),jj1(ix),nt) = sum(GFS_Data(nb)%Sfcprop%tslb(ix,4:Model%lsoil_lsm))
         enddo
       endif lsm_choice

       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%t2m)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%q2m)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%nirbmdi)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%nirdfdi)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%visbmdi)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%visdfdi)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%nirbmui)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%nirdfui)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%visbmui)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%visdfui)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%sfcdsw)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%sfcnsw)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%sfcdlw)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Grid%xlon)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Grid%xlat)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Grid%xlat_d)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Grid%sinlat)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Grid%coslat)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Grid%area)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Grid%dx)
       if (Model%ntoz > 0) then
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Grid%ddy_o3)
       endif
       if (Model%h2o_phys) then
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Grid%ddy_h)
       endif
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Cldprop%cv)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Cldprop%cvt)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Cldprop%cvb)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Radtend%sfalb)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Radtend%coszen)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Radtend%tsflw)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Radtend%semis)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Radtend%coszdg)

       ! Radtend%sfcfsw is an array of derived type, so we copy all
       ! eight elements of the type in one loop
       do ix=1,Atm_block%blksz(nb)
         temp2d(ii1(ix),jj1(ix),nt+1) = GFS_Data(nb)%Radtend%sfcfsw(ix)%upfxc
         temp2d(ii1(ix),jj1(ix),nt+2) = GFS_Data(nb)%Radtend%sfcfsw(ix)%upfx0
         temp2d(ii1(ix),jj1(ix),nt+3) = GFS_Data(nb)%Radtend%sfcfsw(ix)%dnfxc
         temp2d(ii1(ix),jj1(ix),nt+4) = GFS_Data(nb)%Radtend%sfcfsw(ix)%dnfx0
         temp2d(ii1(ix),jj1(ix),nt+5) = GFS_Data(nb)%Radtend%sfcflw(ix)%upfxc
         temp2d(ii1(ix),jj1(ix),nt+6) = GFS_Data(nb)%Radtend%sfcflw(ix)%upfx0
         temp2d(ii1(ix),jj1(ix),nt+7) = GFS_Data(nb)%Radtend%sfcflw(ix)%dnfxc
         temp2d(ii1(ix),jj1(ix),nt+8) = GFS_Data(nb)%Radtend%sfcflw(ix)%dnfx0
       enddo
       nt = nt + 8

       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tiice(:,1))
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tiice(:,2))
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%albdirvis_lnd)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%albdirnir_lnd)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%albdifvis_lnd)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%albdifnir_lnd)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%emis_lnd)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%emis_ice)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%sncovr_ice)

       if (Model%use_cice_alb .or. Model%lsm == Model%lsm_ruc) then
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%albdirvis_ice)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%albdirnir_ice)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%albdifvis_ice)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%albdifnir_ice)
       endif

       lsm_choice_2: if (Model%lsm == Model%lsm_noahmp) then
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%snowxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tvxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tgxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%canicexy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%canliqxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%eahxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tahxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%cmxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%chxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%fwetxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%sneqvoxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%alboldxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%qsnowxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%wslakexy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%zwtxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%waxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%wtxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%lfmassxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%rtmassxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%stmassxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%woodxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%stblcpxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%fastcpxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%xsaixy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%xlaixy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%taussxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%smcwtdxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%deeprechxy)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%rechxy)

         ! These five arrays use bizarre indexing, so we use loops:
         do k=-2,0
           call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%snicexy(:,k))
         enddo

         do k=-2,0
           call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%snliqxy(:,k))
         enddo

         do k=-2,0
           call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tsnoxy(:,k))
         enddo

         do k=1,4
           call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%smoiseq(:,k))
         enddo

         do k=-2,4
           call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%zsnsoxy(:,k))
         enddo
       elseif (Model%lsm == Model%lsm_ruc) then
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%wetness)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%clw_surf_land)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%clw_surf_ice)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%qwv_surf_land)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%qwv_surf_ice)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tsnow_land)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tsnow_ice)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%snowfallac_land)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%snowfallac_ice)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%sfalb_lnd)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%sfalb_lnd_bck)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%sfalb_ice)
         if (Model%rdlai) then
           call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%xlaixy)
         endif
       endif lsm_choice_2

       nstf_name_choice: if (Model%nstf_name(1) > 0) then
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tref)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%z_c)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%c_0)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%c_d)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%w_0)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%w_d)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%xt)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%xs)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%xu)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%xz)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%zm)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%xtts)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%xzts)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%ifd)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%dt_cool)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%qrain)
       endif nstf_name_choice

! Flake
       if (Model%lkm > 0 .and. Model%iopt_lake==Model%iopt_lake_flake) then
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%T_snow)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%T_ice)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%h_ML)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%t_ML)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%t_mnw)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%h_talb)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%t_talb)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%t_bot1)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%t_bot2)
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%c_t)
       endif
       
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Tbd%phy_f2d)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Tbd%phy_fctd)

       ! Copy to temp3dlevsp1
       nt=0

       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3dlevsp1, GFS_Data(nb)%Statein%phii)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3dlevsp1, GFS_Data(nb)%Statein%prsi)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3dlevsp1, GFS_Data(nb)%Statein%prsik)

       ! Copy to temp3d
       nt=0

       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Statein%phil)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Statein%prsl)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Statein%prslk)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Statein%ugrs)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Statein%vgrs)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Statein%vvl)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Statein%tgrs)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Stateout%gu0)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Stateout%gv0)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Stateout%gt0)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Radtend%htrsw)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Radtend%htrlw)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Radtend%swhc)
       call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Radtend%lwhc)
       do l = 1,Model%ntot3d
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Tbd%phy_f3d(:,:,l))
       enddo
       do l = 1,ntr
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Statein%qgrs(:,:,l))
         call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Stateout%gq0(:,:,l))
       enddo
   enddo block_loop


   outunit = stdout()
   do i = 1,nsfcprop2d+Model%ntot2d+Model%nctp
     write (name, '(i3.3,3x,4a)') i, ' 2d '
     write(outunit,100) name, mpp_chksum(temp2d(:,:,i:i))
   enddo
   do i = 1,3
     write (name, '(i2.2,3x,4a)') i, ' 3d levsp1'
     write(outunit,100) name, mpp_chksum(temp3dlevsp1(:,:,:,i:i))
   enddo
   do i = 1,14+Model%ntot3d+2*ntr
     write (name, '(i2.2,3x,4a)') i, ' 3d levs'
     write(outunit,100) name, mpp_chksum(temp3d(:,:,:,i:i))
   enddo
100 format("CHECKSUM::",A32," = ",Z20)

   deallocate(temp2d)
   deallocate(temp3d)
   deallocate(temp3dlevsp1)
   end subroutine FV3GFS_GFS_checksum

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!----------------------------------------------------------------------
! sfc_prop_restart_read
!----------------------------------------------------------------------
!    creates and populates a data type which is then used to "register"
!    restart variables with the GFDL FMS restart subsystem.
!    calls a GFDL FMS routine to restore the data from a restart file.
!    calculates sncovr if it is not present in the restart file.
!
!    calls:  register_restart_field, restart_state, free_restart
!
!    opens:  oro_data.tile?.nc, sfc_data.tile?.nc
!
!----------------------------------------------------------------------
  subroutine sfc_prop_restart_read (Sfcprop, Atm_block, Model, fv_domain, warm_start, ignore_rst_cksum)
    use rrfs_sd_io
    implicit none
    !--- interface variable definitions
    type(GFS_sfcprop_type),    intent(inout) :: Sfcprop(:)
    type (block_control_type), intent(in)    :: Atm_block
    type(GFS_control_type),    intent(inout) :: Model
    type (domain2d),           intent(in)    :: fv_domain
    logical,                   intent(in)    :: warm_start
    logical,                   intent(in)    :: ignore_rst_cksum
    !--- local variables
    integer :: i, j, k, ix, lsoil, num, nb, i_start, j_start, i_end, j_end, nt, n
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: id_restart
    integer :: nvar_o2
    integer :: nvar_oro_ls_ss
    integer :: nvar_vegfr, nvar_soilfr
    integer :: isnow
    integer, allocatable :: ii1(:), jj1(:)
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p  => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p  => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p1 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p2 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p3 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_fr => NULL()
    !--- local variables for sncovr calculation
    integer :: vegtyp
    logical :: mand
    real(kind=kind_phys) :: rsnow, tem, tem1
    !--- directory of the input files
    character(5)  :: indir='INPUT'
    character(37) :: infile
    !--- fms2_io file open logic
    logical :: amiopen
    logical :: override_frac_grid

    type(clm_lake_data_type) :: clm_lake
    type(rrfs_sd_state_type) :: rrfs_sd_state
    type(rrfs_sd_emissions_type) :: rrfs_sd_emis

    nvar_oro_ls_ss = 10

    nvar_vegfr  = Model%nvegcat
    nvar_soilfr = Model%nsoilcat

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)

    !--- OROGRAPHY FILE

    !--- open file
    infile=trim(indir)//'/'//trim(fn_oro)
    amiopen=open_file(Oro_restart, trim(infile), 'read', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
    if (.not.amiopen) call mpp_error( FATAL, 'Error with opening file '//trim(infile) )

    call oro%register(Model,Oro_restart,Atm_block)

   !--- read the orography restart/data
   call mpp_error(NOTE,'reading topographic/orographic information from INPUT/oro_data.tile*.nc')
   call read_restart(Oro_restart, ignore_checksum=ignore_rst_cksum)
   call close_file(Oro_restart)

   !--- copy data into GFS containers
   call oro%copy(Model, Sfcprop, Atm_block)

    if_smoke: if(Model%rrfs_sd) then  ! for RRFS-SD

      !--- Dust input FILE
      !--- open file
      infile=trim(indir)//'/'//trim(fn_dust12m)
      amiopen=open_file(dust12m_restart, trim(infile), 'read', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
      if (.not.amiopen) call mpp_error( FATAL, 'Error with opening file'//trim(infile) )
      
      !--- Register axes and variables, allocate memory:
      call rrfs_sd_emis%register_dust12m(dust12m_restart, Atm_block)
      
      !--- read new GSL created dust12m restart/data
      call mpp_error(NOTE,'reading dust12m information from INPUT/dust12m_data.tile*.nc')
      call read_restart(dust12m_restart)
      call close_file(dust12m_restart)
      
      !--- Copy to Sfcprop and free temporary arrays:
      call rrfs_sd_emis%copy_dust12m(Sfcprop, Atm_block)

      !----------------------------------------------

      !--- open anthropogenic emission file
      infile=trim(indir)//'/'//trim(fn_emi)
      amiopen=open_file(emi_restart, trim(infile), 'read', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
      if (.not.amiopen) call mpp_error( FATAL, 'Error with opening file'//trim(infile) )

      ! Register axes and variables, allocate memory
      call rrfs_sd_emis%register_emi(emi_restart, Atm_block)
      
      !--- read anthropogenic emi restart/data
      call mpp_error(NOTE,'reading emi information from INPUT/emi_data.tile*.nc')
      call read_restart(emi_restart)
      call close_file(emi_restart)
      
      !--- Copy to Sfcprop and free temporary arrays:
      call rrfs_sd_emis%copy_emi(Sfcprop, Atm_block)

      !----------------------------------------------
    
      !--- Dust input FILE
      !--- open file
      infile=trim(indir)//'/'//trim(fn_rrfssd)
      amiopen=open_file(rrfssd_restart, trim(infile), 'read', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
      if (.not.amiopen) call mpp_error( FATAL, 'Error with opening file'//trim(infile) )

      ! Register axes and variables, allocate memory
      call rrfs_sd_emis%register_fire(rrfssd_restart, Atm_block)

      !--- read new GSL created rrfssd restart/data
      call mpp_error(NOTE,'reading rrfssd information from INPUT/SMOKE_RRFS_data.nc')
      call read_restart(rrfssd_restart)
      call close_file(rrfssd_restart)

      !--- Copy to Sfcprop and free temporary arrays:
      call rrfs_sd_emis%copy_fire(Sfcprop, Atm_block)

    endif if_smoke  ! RRFS_SD

    !--- Modify/read-in additional orographic static fields for GSL drag suite
    if (Model%gwd_opt==3 .or. Model%gwd_opt==33 .or. &
        Model%gwd_opt==2 .or. Model%gwd_opt==22 ) then

      !--- open restart file
      infile=trim(indir)//'/'//trim(fn_oro_ls)
      amiopen=open_file(Oro_ls_restart, trim(infile), 'read', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
      if( .not.amiopen ) call mpp_error( FATAL, 'Error with opening file '//trim(infile) )

      !--- open restart file
      infile=trim(indir)//'/'//trim(fn_oro_ss)
      amiopen=open_file(Oro_ss_restart, trim(infile), 'read', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
      if( .not.amiopen ) call mpp_error( FATAL, 'Error with opening file '//trim(infile) )

      if (.not. allocated(oro_ls_ss_name)) then
      !--- allocate the various containers needed for orography data
        allocate(oro_ls_ss_name(nvar_oro_ls_ss))
        allocate(oro_ls_var(nx,ny,nvar_oro_ls_ss))
        allocate(oro_ss_var(nx,ny,nvar_oro_ls_ss))

        oro_ls_ss_name(1)  = 'stddev'
        oro_ls_ss_name(2)  = 'convexity'
        oro_ls_ss_name(3)  = 'oa1'
        oro_ls_ss_name(4)  = 'oa2'
        oro_ls_ss_name(5)  = 'oa3'
        oro_ls_ss_name(6)  = 'oa4'
        oro_ls_ss_name(7)  = 'ol1'
        oro_ls_ss_name(8)  = 'ol2'
        oro_ls_ss_name(9)  = 'ol3'
        oro_ls_ss_name(10) = 'ol4'

        call register_axis(Oro_ls_restart, "lon", 'X')
        call register_axis(Oro_ls_restart, "lat", 'Y')
        call register_axis(Oro_ss_restart, "lon", 'X')
        call register_axis(Oro_ss_restart, "lat", 'Y')

        do num = 1,nvar_oro_ls_ss
          var2_p => oro_ls_var(:,:,num)
          call register_restart_field(Oro_ls_restart, oro_ls_ss_name(num), var2_p, dimensions=(/'lon','lat'/))
        enddo
        nullify(var2_p)
        do num = 1,nvar_oro_ls_ss
          var2_p => oro_ss_var(:,:,num)
          call register_restart_field(Oro_ss_restart, oro_ls_ss_name(num), var2_p, dimensions=(/'lon','lat'/))
        enddo
        nullify(var2_p)
      end if

      !--- read new GSL created orography restart/data
      call mpp_error(NOTE,'reading topographic/orographic information from &
           &INPUT/oro_data_ls.tile*.nc')
      call read_restart(Oro_ls_restart, ignore_checksum=ignore_rst_cksum)
      call close_file(Oro_ls_restart)
      call mpp_error(NOTE,'reading topographic/orographic information from &
           &INPUT/oro_data_ss.tile*.nc')
      call read_restart(Oro_ss_restart, ignore_checksum=ignore_rst_cksum)
      call close_file(Oro_ss_restart)


      do nb = 1, Atm_block%nblks
        !--- 2D variables
        do ix = 1, Atm_block%blksz(nb)
          i = Atm_block%index(nb)%ii(ix) - isc + 1
          j = Atm_block%index(nb)%jj(ix) - jsc + 1
          ! Replace hprime(1:10) with GSL oro stat data only when using GSL
          ! drag suite with large scale GWD and blocking as part of unified drag
          ! suite. Otherwise, original oro stat data is used.
          if ( (Model%gwd_opt==3 .or. Model%gwd_opt==33) .or.    &
               ( (Model%gwd_opt==2 .or. Model%gwd_opt==22) .and. &
                  Model%do_gsl_drag_ls_bl ) ) then
            !--- assign hprime(1:10) and hprime(15:24) with new oro stat data
            Sfcprop(nb)%hprime(ix,1)  = oro_ls_var(i,j,1)
            Sfcprop(nb)%hprime(ix,2)  = oro_ls_var(i,j,2)
            Sfcprop(nb)%hprime(ix,3)  = oro_ls_var(i,j,3)
            Sfcprop(nb)%hprime(ix,4)  = oro_ls_var(i,j,4)
            Sfcprop(nb)%hprime(ix,5)  = oro_ls_var(i,j,5)
            Sfcprop(nb)%hprime(ix,6)  = oro_ls_var(i,j,6)
            Sfcprop(nb)%hprime(ix,7)  = oro_ls_var(i,j,7)
            Sfcprop(nb)%hprime(ix,8)  = oro_ls_var(i,j,8)
            Sfcprop(nb)%hprime(ix,9)  = oro_ls_var(i,j,9)
            Sfcprop(nb)%hprime(ix,10)  = oro_ls_var(i,j,10)
          endif
          Sfcprop(nb)%hprime(ix,15)  = oro_ss_var(i,j,1)
          Sfcprop(nb)%hprime(ix,16)  = oro_ss_var(i,j,2)
          Sfcprop(nb)%hprime(ix,17)  = oro_ss_var(i,j,3)
          Sfcprop(nb)%hprime(ix,18)  = oro_ss_var(i,j,4)
          Sfcprop(nb)%hprime(ix,19)  = oro_ss_var(i,j,5)
          Sfcprop(nb)%hprime(ix,20)  = oro_ss_var(i,j,6)
          Sfcprop(nb)%hprime(ix,21)  = oro_ss_var(i,j,7)
          Sfcprop(nb)%hprime(ix,22)  = oro_ss_var(i,j,8)
          Sfcprop(nb)%hprime(ix,23)  = oro_ss_var(i,j,9)
          Sfcprop(nb)%hprime(ix,24)  = oro_ss_var(i,j,10)
        enddo
      enddo

   end if

   !--- SURFACE FILE

   !--- open file
   infile=trim(indir)//'/'//trim(fn_srf)
   amiopen=open_file(Sfc_restart, trim(infile), "read", domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
   if( .not.amiopen ) call mpp_error(FATAL, 'Error opening file'//trim(infile))

   if(sfc%allocate_arrays(Model, Atm_block, .true., warm_start)) then
      call sfc%fill_2d_names(Model, warm_start)
      call sfc%register_axes(Model, Sfc_restart, .true., warm_start)

      ! Tell CLM Lake to allocate data, and register its axes and fields
      if(Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_clm) then
        call clm_lake%allocate_data(Model)
        call clm_lake%copy_to_temporaries(Model,Sfcprop,Atm_block)
        call clm_lake%register_axes(Model, Sfc_restart)
        call clm_lake%register_fields(Sfc_restart)
      endif

      if(Model%rrfs_sd) then
        call rrfs_sd_state%allocate_data(Model)
        call rrfs_sd_state%fill_data(Model, Sfcprop, Atm_block)
        call rrfs_sd_state%register_axis(Model, Sfc_restart)
        call rrfs_sd_state%register_fields(Sfc_restart)
      endif

      call sfc%register_2d_fields(Model,Sfc_restart,.true.,warm_start)
   endif  ! if not allocated

   call sfc%fill_3d_names(Model,warm_start)
   call sfc%register_3d_fields(Model,Sfc_restart,.true.,warm_start)
   call sfc%init_fields(Model)

    !--- read the surface restart/data
    call mpp_error(NOTE,'reading surface properties data from INPUT/sfc_data.tile*.nc')
    call read_restart(Sfc_restart, ignore_checksum=ignore_rst_cksum)
    call close_file(Sfc_restart)

    ! Tell clm_lake to copy data to temporary arrays
    if(Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_clm) then
      call clm_lake%copy_from_temporaries(Model,Sfcprop,Atm_block)
    endif

    if(Model%rrfs_sd) then
      call rrfs_sd_state%copy_from_temporaries(Model,Sfcprop,Atm_block)
    end if

!   write(0,*)' stype read in min,max=',minval(sfc%var2(:,:,35)),maxval(sfc%var2(:,:,35)),' sfc%name2=',sfc%name2(35)
!   write(0,*)' stype read in min,max=',minval(sfc%var2(:,:,18)),maxval(sfc%var2(:,:,18))
!   write(0,*)' sfc%var2=',sfc%var2(:,:,12)

    !--- place the data into the block GFS containers
    override_frac_grid=Model%frac_grid
    call sfc%copy_to_grid(Model, Atm_block, Sfcprop, warm_start, override_frac_grid)
    Model%frac_grid=override_frac_grid

    call mpp_error(NOTE, 'gfs_driver:: - after put to container ')

    call sfc%apply_safeguards(Model, Atm_block, Sfcprop)

    ! A standard-compliant Fortran 2003 compiler will call clm_lake_final and rrfs_sd_final here.

  end subroutine sfc_prop_restart_read


!----------------------------------------------------------------------
! sfc_prop_restart_write
!----------------------------------------------------------------------
!    routine to write out GFS surface restarts via the GFDL FMS restart
!    subsystem.
!    takes an optional argument to append timestamps for intermediate
!    restarts.
!
!    calls:  register_restart_field, save_restart
!----------------------------------------------------------------------
  subroutine sfc_prop_restart_write (Sfcprop, Atm_block, Model, fv_domain, timestamp)
    use rrfs_sd_io
    implicit none
    !--- interface variable definitions
    type(GFS_sfcprop_type),      intent(in) :: Sfcprop(:)
    type(block_control_type),    intent(in) :: Atm_block
    type(GFS_control_type),      intent(in) :: Model
    type(domain2d),              intent(in) :: fv_domain
    character(len=32), optional, intent(in) :: timestamp
    !--- local variables
    integer :: i, j, k, nb, ix, lsoil, num, nt
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: id_restart
    logical :: mand
    integer, allocatable :: ii1(:), jj1(:)
    character(len=32) :: fn_srf = 'sfc_data.nc'
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p  => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p  => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p1 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p2 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p3 => NULL()
    real(kind_phys) :: ice
    !--- directory of the input files
    character(7)  :: indir='RESTART'
    character(72) :: infile
    !--- fms2_io file open logic
    logical :: amiopen
    !--- variables used for fms2_io register axis
    integer :: is, ie
    integer, allocatable, dimension(:) :: buffer
    type(clm_lake_data_type), target :: clm_lake
    !--- temporary variables for storing rrfs_sd fields
    type(rrfs_sd_state_type) :: rrfs_sd_state

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)

    !--- set filename
    infile=trim(indir)//'/'//trim(fn_srf)
    if( present(timestamp) ) infile=trim(indir)//'/'//trim(timestamp)//'.'//trim(fn_srf)

    !--- register axis
    amiopen=open_file(Sfc_restart, trim(infile), 'overwrite', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
    if_amiopen: if( amiopen ) then
      call sfc%register_axes(Model, Sfc_restart, .false., .true.)
      call sfc%write_axes(Model, Sfc_restart)
    else
      call mpp_error(FATAL, 'Error in opening file'//trim(infile) )
    end if if_amiopen
    
    ! Tell clm_lake to allocate data, register its axes, and call write_data for each axis's variable
    if(Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_clm) then
      call clm_lake%allocate_data(Model)
      call clm_lake%register_axes(Model, Sfc_restart)
      call clm_lake%write_axes(Model, Sfc_restart)
    endif

    if(Model%rrfs_sd) then
      call rrfs_sd_state%allocate_data(Model)
      call rrfs_sd_state%register_axis(Model,Sfc_restart)
      call rrfs_sd_state%write_axis(Model,Sfc_restart)
    end if

    if (sfc%allocate_arrays(Model, Atm_block, .false., .true.)) then
      call sfc%fill_2d_names(Model,.true.)
    end if

   if(Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_clm) then
       ! Tell clm_lake to register all of its fields
       call clm_lake%register_fields(Sfc_restart)
   endif

   if(Model%rrfs_sd) then
     call rrfs_sd_state%register_fields(Sfc_restart)
   endif

   ! Register 2D surface property fields (except lake, smoke, and dust)
   call sfc%register_2d_fields(Model, Sfc_restart, .false., .true.)

   ! Determine list of 3D surface property fields names:
   call sfc%fill_3d_names(Model, .true.)

   ! Register 3D surface property fields (except lake, smoke, and dust)
   call sfc%register_3d_fields(Model, Sfc_restart, .false., .true.)

    ! Tell clm_lake to copy Sfcprop data to its internal temporary arrays.
    if(Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_clm) then
      call clm_lake%copy_to_temporaries(Model,Sfcprop,Atm_block)
    endif

    if(Model%rrfs_sd) then
     call rrfs_sd_state%copy_to_temporaries(Model,Sfcprop,Atm_block)
    endif

    call sfc%copy_from_grid(Model, Atm_block, Sfcprop)

    call write_restart(Sfc_restart)
    call close_file(Sfc_restart)

    ! A standard-compliant Fortran 2003 compiler will call rrfs_sd_final and clm_lake_final here

  end subroutine sfc_prop_restart_write

!----------------------------------------------------------------------
! phys_restart_read
!----------------------------------------------------------------------
!    creates and populates a data type which is then used to "register"
!    restart variables with the GFDL FMS restart subsystem.
!    calls a GFDL FMS routine to restore the data from a restart file.
!    calculates sncovr if it is not present in the restart file.
!
!    calls:  register_restart_field, restart_state, free_restart
!
!    opens:  phys_data.tile?.nc
!
!----------------------------------------------------------------------
  subroutine phys_restart_read (GFS_Restart, Atm_block, Model, fv_domain, ignore_rst_cksum)
    !--- interface variable definitions
    type(GFS_restart_type),      intent(in) :: GFS_Restart
    type(block_control_type),    intent(in) :: Atm_block
    type(GFS_control_type),      intent(in) :: Model
    type(domain2d),              intent(in) :: fv_domain
    logical,                     intent(in) :: ignore_rst_cksum
    !--- local variables
    integer :: i, j, k, nb, ix, num
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: id_restart
    integer :: nvar2d, nvar3d, fdiag, ldiag
    character(len=64) :: fname
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p => NULL()
    !--- directory of the input files
    character(5)  :: indir='INPUT'
    logical :: amiopen

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)

    nvar2d = GFS_Restart%num2d
    nvar3d = GFS_Restart%num3d
    fdiag  = GFS_Restart%fdiag
    ldiag  = GFS_Restart%ldiag

    !--- open restart file and register axes
    fname = trim(indir)//'/'//trim(fn_phy)
    amiopen=open_file(Phy_restart, trim(fname), 'read', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
    if( amiopen ) then
       call register_axis(Phy_restart, 'xaxis_1', 'X')
       call register_axis(Phy_restart, 'yaxis_1', 'Y')
       call register_axis(Phy_restart, 'zaxis_1', npz)
       call register_axis(Phy_restart, 'Time', unlimited)
    else
       call mpp_error(NOTE,'No physics restarts - cold starting physical parameterizations')
       return
    endif

    !--- register the restart fields
    if (.not. allocated(phy_var2)) then
      allocate (phy_var2(nx,ny,nvar2d))
      allocate (phy_var3(nx,ny,npz,nvar3d))
      phy_var2 = zero
      phy_var3 = zero

      do num = 1,nvar2d
        var2_p => phy_var2(:,:,num)
        call register_restart_field(Phy_restart, trim(GFS_Restart%name2d(num)), var2_p, dimensions=(/'xaxis_1','yaxis_1','Time   '/),&
                                   &is_optional=.true.)
      enddo
      do num = 1,nvar3d
        var3_p => phy_var3(:,:,:,num)
        call register_restart_field(Phy_restart, trim(GFS_restart%name3d(num)), var3_p, dimensions=(/'xaxis_1','yaxis_1','zaxis_1','Time   '/), is_optional=.true.)
      enddo
      nullify(var2_p)
      nullify(var3_p)
    endif

    !--- read the surface restart/data
    call mpp_error(NOTE,'reading physics restart data from INPUT/phy_data.tile*.nc')
    call read_restart(Phy_restart, ignore_checksum=ignore_rst_cksum)
    call close_file(Phy_restart)

    !--- place the data into the block GFS containers
    !--- phy_var* variables
!$omp parallel do default(shared) private(i, j, nb, ix)
    do num = 1,nvar2d
      do nb = 1,Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          i = Atm_block%index(nb)%ii(ix) - isc + 1
          j = Atm_block%index(nb)%jj(ix) - jsc + 1
          GFS_Restart%data(nb,num)%var2p(ix) = phy_var2(i,j,num)
        enddo
      enddo
    enddo
    !-- if restart from init time, reset accumulated diag fields
    if( Model%phour < 1.e-7) then
      do num = fdiag,ldiag
!$omp parallel do default(shared) private(i, j, nb, ix)
        do nb = 1,Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            i = Atm_block%index(nb)%ii(ix) - isc + 1
            j = Atm_block%index(nb)%jj(ix) - jsc + 1
            GFS_Restart%data(nb,num)%var2p(ix) = zero
          enddo
        enddo
      enddo
    endif
    do num = 1,nvar3d
!$omp parallel do default(shared) private(i, j, k, nb, ix)
      do nb = 1,Atm_block%nblks
        do k=1,npz
          do ix = 1, Atm_block%blksz(nb)
            i = Atm_block%index(nb)%ii(ix) - isc + 1
            j = Atm_block%index(nb)%jj(ix) - jsc + 1
            GFS_Restart%data(nb,num)%var3p(ix,k) = phy_var3(i,j,k,num)
          enddo
        enddo
      enddo
    enddo

  end subroutine phys_restart_read


!----------------------------------------------------------------------
! phys_restart_write
!----------------------------------------------------------------------
!    routine to write out GFS surface restarts via the GFDL FMS restart
!    subsystem.
!    takes an optional argument to append timestamps for intermediate
!    restarts.
!
!    calls:  register_restart_field, save_restart
!----------------------------------------------------------------------
  subroutine phys_restart_write (GFS_Restart, Atm_block, Model, fv_domain, timestamp)
    !--- interface variable definitions
    type(GFS_restart_type),      intent(in) :: GFS_Restart
    type(block_control_type),    intent(in) :: Atm_block
    type(GFS_control_type),      intent(in) :: Model
    type(domain2d),              intent(in) :: fv_domain
    character(len=32), optional, intent(in) :: timestamp
    !--- local variables
    integer :: i, j, k, nb, ix, num
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: id_restart
    integer :: nvar2d, nvar3d
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p => NULL()
    !--- used for axis data for fms2_io
    integer :: is, ie
    integer, allocatable, dimension(:) :: buffer
    character(7) :: indir='RESTART'
    character(72) :: infile
    logical :: amiopen

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)
    nvar2d = GFS_Restart%num2d
    nvar3d = GFS_Restart%num3d

    !--- set file name
    infile=trim(indir)//'/'//trim(fn_phy)
    if( present(timestamp) ) infile=trim(indir)//'/'//trim(timestamp)//'.'//trim(fn_phy)
    !--- register axis
    amiopen=open_file(Phy_restart, trim(infile), 'overwrite', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
    if( amiopen ) then
      call register_axis(Phy_restart, 'xaxis_1', 'X')
      call register_field(Phy_restart, 'xaxis_1', 'double', (/'xaxis_1'/))
      call register_variable_attribute(Phy_restart, 'xaxis_1', 'cartesian_axis', 'X', str_len=1)
      call get_global_io_domain_indices(Phy_restart, 'xaxis_1', is, ie, indices=buffer)
      call write_data(Phy_restart, "xaxis_1", buffer)
      deallocate(buffer)

      call register_axis(Phy_restart, 'yaxis_1', 'Y')
      call register_field(Phy_restart, 'yaxis_1', 'double', (/'yaxis_1'/))
      call register_variable_attribute(Phy_restart, 'yaxis_1', 'cartesian_axis', 'Y', str_len=1)
      call get_global_io_domain_indices(Phy_restart, 'yaxis_1', is, ie, indices=buffer)
      call write_data(Phy_restart, "yaxis_1", buffer)
      deallocate(buffer)

      call register_axis(Phy_restart, 'zaxis_1', npz)
      call register_field(Phy_restart, 'zaxis_1', 'double', (/'zaxis_1'/))
      call register_variable_attribute(Phy_restart, 'zaxis_1', 'cartesian_axis', 'Z', str_len=1)
      allocate( buffer(npz) )
      do i=1, npz
         buffer(i)=i
      end do
      call write_data(Phy_restart, "zaxis_1", buffer)
      deallocate(buffer)

      call register_axis(Phy_restart, 'Time', unlimited)
      call register_field(Phy_restart, 'Time', 'double', (/'Time'/))
      call register_variable_attribute(Phy_restart, 'Time', 'cartesian_axis', 'T', str_len=1)
      call write_data(Phy_restart, "Time", 1)
    else
      call mpp_error(FATAL, 'Error opening file '//trim(infile))
    end if

    !--- register the restart fields
    if (.not. allocated(phy_var2)) then
      allocate (phy_var2(nx,ny,nvar2d))
      allocate (phy_var3(nx,ny,npz,nvar3d))
      phy_var2 = zero
      phy_var3 = zero
    endif

    do num = 1,nvar2d
       var2_p => phy_var2(:,:,num)
       call register_restart_field(Phy_restart, trim(GFS_Restart%name2d(num)), var2_p, dimensions=(/'xaxis_1','yaxis_1','Time   '/),&
                                  &is_optional=.true.)
    enddo
    do num = 1,nvar3d
       var3_p => phy_var3(:,:,:,num)
       call register_restart_field(Phy_restart, trim(GFS_Restart%name3d(num)), var3_p, dimensions=(/'xaxis_1','yaxis_1','zaxis_1','Time   '/),&
                                  &is_optional=.true.)
    enddo
    nullify(var2_p)
    nullify(var3_p)

    !--- 2D variables
!$omp parallel do default(shared) private(i, j, num, nb, ix)
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
!$omp parallel do default(shared) private(i, j, k, num, nb, ix)
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

    call write_restart(Phy_restart)
    call close_file(Phy_restart)

  end subroutine phys_restart_write

!-------------------------------------------------------------------------
!--- gfdl_diag_register ---
!-------------------------------------------------------------------------
!    creates and populates a data type which is then used to "register"
!    GFS physics diagnostic variables with the GFDL FMS diagnostic manager.
!    includes short & long names, units, conversion factors, etc.
!    there is no copying of data, but instead a clever use of pointers.
!    calls a GFDL FMS routine to register diagnositcs and compare against
!    the diag_table to determine what variables are to be output.
!
!    calls:  register_diag_field
!-------------------------------------------------------------------------
  subroutine fv3gfs_diag_register(Diag, Time, Atm_block, Model, xlon, xlat, axes)
    use physcons,  only: con_g
!--- subroutine interface variable definitions
    type(GFS_externaldiag_type),       intent(inout) :: Diag(:)
    type(time_type),           intent(in)    :: Time
    type (block_control_type), intent(in)    :: Atm_block
    type(GFS_control_type),    intent(in)    :: Model
    real(kind=kind_phys),      intent(in)    :: xlon(:,:)
    real(kind=kind_phys),      intent(in)    :: xlat(:,:)
    integer, dimension(4),     intent(in)    :: axes
!--- local variables
    integer :: idx, nrgst_bl, nrgst_nb, nrgst_vctbl

    isco   = Atm_block%isc
    ieco   = Atm_block%iec
    jsco   = Atm_block%jsc
    jeco   = Atm_block%jec
    levo   = model%levs
    fhzero = nint(Model%fhzero)
!   ncld   = Model%ncld
    ncld   = Model%imp_physics
    nsoil  = Model%lsoil
    dtp    = Model%dtp
    imp_physics  = Model%imp_physics
    landsfcmdl  = Model%lsm
!    print *,'in fv3gfs_diag_register,ncld=',Model%ncld,Model%lsoil,Model%imp_physics, &
!      ' dtp=',dtp,' landsfcmdl=',Model%lsm
!
!save lon/lat for vector interpolation
    allocate(lon(isco:ieco,jsco:jeco))
    allocate(lat(isco:ieco,jsco:jeco))
    lon = xlon
    lat = xlat

    do idx = 1,DIAG_SIZE
      if (trim(Diag(idx)%name) == '') exit
      tot_diag_idx = idx
    enddo

    if (tot_diag_idx == DIAG_SIZE) then
      call mpp_error(fatal, 'FV3GFS_io::fv3gfs_diag_register - need to increase parameter DIAG_SIZE')
    endif

    allocate(nstt(tot_diag_idx), nstt_vctbl(tot_diag_idx))
    nstt          = 0
    nstt_vctbl    = 0
    nrgst_bl      = 0
    nrgst_nb      = 0
    nrgst_vctbl   = 0
    num_axes_phys = 2
    do idx = 1,tot_diag_idx
      if (diag(idx)%axes == -99) then
        call mpp_error(fatal, 'gfs_driver::gfs_diag_register - attempt to register an undefined variable')
      endif
      Diag(idx)%id = register_diag_field (trim(Diag(idx)%mod_name), trim(Diag(idx)%name),  &
                                          axes(1:Diag(idx)%axes), Time, trim(Diag(idx)%desc), &
                                          trim(Diag(idx)%unit), missing_value=real(missing_value))
      if(Diag(idx)%id > 0) then
        if (Diag(idx)%axes == 2) then
           if( index(trim(Diag(idx)%intpl_method),'bilinear') > 0 ) then
             nrgst_bl = nrgst_bl + 1
             nstt(idx) = nrgst_bl
           else if (trim(Diag(idx)%intpl_method) == 'nearest_stod' ) then
             nrgst_nb = nrgst_nb + 1
             nstt(idx) = nrgst_nb
           endif
           if(trim(Diag(idx)%intpl_method) == 'vector_bilinear') then
             if(Diag(idx)%name(1:1) == 'v' .or. Diag(idx)%name(1:1) == 'V') then
               nrgst_vctbl = nrgst_vctbl + 1
               nstt_vctbl(idx) = nrgst_vctbl
!             print *,'in phy_setup, vector_bilinear, name=', trim(Diag(idx)%name),' nstt_vctbl=', nstt_vctbl(idx), 'idx=',idx
             endif
           endif
        else if (diag(idx)%axes == 3) then
           if( index(trim(diag(idx)%intpl_method),'bilinear') > 0 ) then
             nstt(idx) = nrgst_bl + 1
             nrgst_bl  = nrgst_bl + levo
           else if (trim(diag(idx)%intpl_method) == 'nearest_stod' ) then
             nstt(idx) = nrgst_nb + 1
             nrgst_nb  = nrgst_nb + levo
           endif
           if(trim(diag(idx)%intpl_method) == 'vector_bilinear') then
             if(diag(idx)%name(1:1) == 'v' .or. diag(idx)%name(1:1) == 'V') then
               nstt_vctbl(idx) = nrgst_vctbl + 1
               nrgst_vctbl = nrgst_vctbl + levo
!             print *,'in phy_setup, vector_bilinear, name=', trim(diag(idx)%name),' nstt_vctbl=', nstt_vctbl(idx), 'idx=',idx
             endif
           endif
           num_axes_phys = 3
        endif
      endif

    enddo

    total_outputlevel = nrgst_bl + nrgst_nb
    allocate(buffer_phys_bl(isco:ieco,jsco:jeco,nrgst_bl))
    allocate(buffer_phys_nb(isco:ieco,jsco:jeco,nrgst_nb))
    allocate(buffer_phys_windvect(3,isco:ieco,jsco:jeco,nrgst_vctbl))
    buffer_phys_bl = zero
    buffer_phys_nb = zero
    buffer_phys_windvect = zero
    if(mpp_pe() == mpp_root_pe()) print *,'in fv3gfs_diag_register, nrgst_bl=',nrgst_bl,' nrgst_nb=',nrgst_nb, &
       ' nrgst_vctbl=',nrgst_vctbl, 'isco=',isco,ieco,'jsco=',jsco,jeco,' num_axes_phys=', num_axes_phys

  end subroutine fv3gfs_diag_register
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!--- gfs_diag_output ---
!-------------------------------------------------------------------------
!    routine to transfer the diagnostic data to the gfdl fms diagnostic
!    manager for eventual output to the history files.
!
!    calls:  send_data
!-------------------------------------------------------------------------
  subroutine fv3gfs_diag_output(time, diag, atm_block, nx, ny, levs, ntcw, ntoz, &
                                dt, time_int, time_intfull, time_radsw, time_radlw)
!--- subroutine interface variable definitions
    type(time_type),           intent(in) :: time
    type(GFS_externaldiag_type),       intent(in) :: diag(:)
    type (block_control_type), intent(in) :: atm_block
    integer,                   intent(in) :: nx, ny, levs, ntcw, ntoz
    real(kind=kind_phys),      intent(in) :: dt
    real(kind=kind_phys),      intent(in) :: time_int
    real(kind=kind_phys),      intent(in) :: time_intfull
    real(kind=kind_phys),      intent(in) :: time_radsw
    real(kind=kind_phys),      intent(in) :: time_radlw
!--- local variables
    integer :: i, j, k, idx, nblks, nb, ix, ii, jj
    integer :: is_in, js_in, isc, jsc
    character(len=2) :: xtra
#ifdef CCPP_32BIT
    real, dimension(nx*ny)      :: var2p
    real, dimension(nx*ny,levs) :: var3p
    real, dimension(nx,ny)      :: var2
    real, dimension(nx,ny,levs) :: var3
#else
    real(kind=kind_phys), dimension(nx*ny)      :: var2p
    real(kind=kind_phys), dimension(nx*ny,levs) :: var3p
    real(kind=kind_phys), dimension(nx,ny)      :: var2
    real(kind=kind_phys), dimension(nx,ny,levs) :: var3
#endif
    real(kind=kind_phys) :: rdt, rtime_int, rtime_intfull, lcnvfac
    real(kind=kind_phys) :: rtime_radsw, rtime_radlw
    logical :: used

     nblks         = atm_block%nblks
     rdt           = one/dt
     rtime_int     = one/time_int
     rtime_intfull = one/time_intfull
     rtime_radsw   = one/time_radsw
     rtime_radlw   = one/time_radlw

     isc   = atm_block%isc
     jsc   = atm_block%jsc
     is_in = atm_block%isc
     js_in = atm_block%jsc

!     if(mpp_pe()==mpp_root_pe())print *,'in,fv3gfs_io. time avg, time_int=',time_int
     do idx = 1,tot_diag_idx
       if (diag(idx)%id > 0) then
         lcnvfac = diag(idx)%cnvfac
         if (diag(idx)%time_avg) then
           if ( trim(diag(idx)%time_avg_kind) == 'full' ) then
             lcnvfac = lcnvfac*rtime_intfull
!             if(mpp_pe()==mpp_root_pe())print *,'in,fv3gfs_io. full time avg, field=',trim(Diag(idx)%name),' time=',time_intfull
           else if ( trim(diag(idx)%time_avg_kind) == 'rad_lw' ) then
             lcnvfac = lcnvfac*min(rtime_radlw,rtime_int)
!             if(mpp_pe()==mpp_root_pe())print *,'in,fv3gfs_io. rad longwave avg, field=',trim(Diag(idx)%name),' time=',time_radlw
           else if ( trim(diag(idx)%time_avg_kind) == 'rad_sw' ) then
             lcnvfac = lcnvfac*min(rtime_radsw,rtime_int)
!             if(mpp_pe()==mpp_root_pe())print *,'in,fv3gfs_io. rad shortwave avg, field=',trim(Diag(idx)%name),' time=',time_radsw
           else if ( trim(diag(idx)%time_avg_kind) == 'rad_swlw_min' ) then
             lcnvfac = lcnvfac*min(max(rtime_radsw,rtime_radlw),rtime_int)
!             if(mpp_pe()==mpp_root_pe())print *,'in,fv3gfs_io. rad swlw min avg, field=',trim(Diag(idx)%name),' time=',time_radlw,time_radsw,time_int
           else
             lcnvfac = lcnvfac*rtime_int
           endif
         endif
         if (diag(idx)%axes == 2) then
           ! Integer data
           int_or_real: if (associated(Diag(idx)%data(1)%int2)) then
             if (trim(Diag(idx)%intpl_method) == 'nearest_stod') then
               var2(1:nx,1:ny) = 0._kind_phys
               do j = 1, ny
                 jj = j + jsc -1
                 do i = 1, nx
                   ii = i + isc -1
                   nb = Atm_block%blkno(ii,jj)
                   ix = Atm_block%ixp(ii,jj)
                   var2(i,j) = real(Diag(idx)%data(nb)%int2(ix), kind=kind_phys)
                 enddo
               enddo
               call store_data(Diag(idx)%id, var2, Time, idx, Diag(idx)%intpl_method, Diag(idx)%name)
             else
               call mpp_error(FATAL, 'Interpolation method ' // trim(Diag(idx)%intpl_method) // ' for integer variable ' &
                                    // trim(Diag(idx)%name) // ' not supported.')
             endif
           ! Real data
           else ! int_or_real
             if (trim(diag(idx)%mask) == 'positive_flux') then
               !--- albedos are actually a ratio of two radiation surface properties
               var2(1:nx,1:ny) = 0._kind_phys
               do j = 1, ny
                 jj = j + jsc -1
                 do i = 1, nx
                   ii = i + isc -1
                   nb = Atm_block%blkno(ii,jj)
                   ix = Atm_block%ixp(ii,jj)
                   if (Diag(idx)%data(nb)%var21(ix) > 0._kind_phys) &
                     var2(i,j) = max(0._kind_phys,min(1._kind_phys,Diag(idx)%data(nb)%var2(ix)/Diag(idx)%data(nb)%var21(ix)))*lcnvfac
                 enddo
               enddo
             elseif (trim(Diag(idx)%mask) == 'land_ice_only') then
               !--- need to "mask" gflux to output valid data over land/ice only
               var2(1:nx,1:ny) = missing_value
               do j = 1, ny
                 jj = j + jsc -1
                 do i = 1, nx
                   ii = i + isc -1
                   nb = Atm_block%blkno(ii,jj)
                   ix = Atm_block%ixp(ii,jj)
                    if (Diag(idx)%data(nb)%var21(ix) /= 0) var2(i,j) = Diag(idx)%data(nb)%var2(ix)*lcnvfac
                 enddo
               enddo
             elseif (trim(Diag(idx)%mask) == 'land_only') then
               !--- need to "mask" soilm to have value only over land
               var2(1:nx,1:ny) = missing_value
               do j = 1, ny
                 jj = j + jsc -1
                 do i = 1, nx
                   ii = i + isc -1
                   nb = Atm_block%blkno(ii,jj)
                   ix = Atm_block%ixp(ii,jj)
                   if (Diag(idx)%data(nb)%var21(ix) == 1) var2(i,j) = Diag(idx)%data(nb)%var2(ix)*lcnvfac
                 enddo
               enddo
             elseif (trim(Diag(idx)%mask) == 'cldmask') then
               !--- need to "mask" soilm to have value only over land
               var2(1:nx,1:ny) = missing_value
               do j = 1, ny
                 jj = j + jsc -1
                 do i = 1, nx
                   ii = i + isc -1
                   nb = Atm_block%blkno(ii,jj)
                   ix = Atm_block%ixp(ii,jj)
                   if (Diag(idx)%data(nb)%var21(ix)*100. > 0.5) var2(i,j) = Diag(idx)%data(nb)%var2(ix)*lcnvfac
                 enddo
               enddo
             elseif (trim(Diag(idx)%mask) == 'cldmask_ratio') then
               !--- need to "mask" soilm to have value only over land
               var2(1:nx,1:ny) = missing_value
               do j = 1, ny
                 jj = j + jsc -1
                 do i = 1, nx
                   ii = i + isc -1
                   nb = Atm_block%blkno(ii,jj)
                   ix = Atm_block%ixp(ii,jj)
                   if (Diag(idx)%data(nb)%var21(ix)*100.*lcnvfac > 0.5) var2(i,j) = Diag(idx)%data(nb)%var2(ix)/ &
                       Diag(idx)%data(nb)%var21(ix)
                 enddo
               enddo
             elseif (trim(Diag(idx)%mask) == 'pseudo_ps') then
               if ( use_wrtgridcomp_output ) then
                 do j = 1, ny
                   jj = j + jsc -1
                   do i = 1, nx
                     ii = i + isc -1
                     nb = Atm_block%blkno(ii,jj)
                     ix = Atm_block%ixp(ii,jj)
                     var2(i,j) = (Diag(idx)%data(nb)%var2(ix)/stndrd_atmos_ps)**(rdgas/grav*stndrd_atmos_lapse)
                   enddo
                 enddo
               else
                 do j = 1, ny
                   jj = j + jsc -1
                   do i = 1, nx
                     ii = i + isc -1
                     nb = Atm_block%blkno(ii,jj)
                     ix = Atm_block%ixp(ii,jj)
                     var2(i,j) = Diag(idx)%data(nb)%var2(ix)
                   enddo
                 enddo
               endif
             elseif (trim(Diag(idx)%mask) == '') then
               do j = 1, ny
                 jj = j + jsc -1
                 do i = 1, nx
                   ii = i + isc -1
                   nb = Atm_block%blkno(ii,jj)
                   ix = Atm_block%ixp(ii,jj)
                   var2(i,j) = Diag(idx)%data(nb)%var2(ix)*lcnvfac
                 enddo
               enddo
             endif
           endif int_or_real
!           used=send_data(Diag(idx)%id, var2, Time)
!           print *,'in phys, after store_data, idx=',idx,' var=', trim(Diag(idx)%name)
           call store_data(Diag(idx)%id, var2, Time, idx, Diag(idx)%intpl_method, Diag(idx)%name)
!           if(trim(Diag(idx)%name) == 'totprcp_ave' ) print *,'in gfs_io, totprcp=',Diag(idx)%data(1)%var2(1:3), &
!             ' lcnvfac=', lcnvfac
         elseif (Diag(idx)%axes == 3) then
         !---
         !--- skipping other 3D variables with the following else statement
         !---
!         if(mpp_pe()==mpp_root_pe())print *,'in,fv3gfs_io. 3D fields, idx=',idx,'varname=',trim(diag(idx)%name), &
!             'lcnvfac=',lcnvfac, 'levo=',levo,'nx=',nx,'ny=',ny
           do k=1, levo
             do j = 1, ny
               jj = j + jsc -1
               do i = 1, nx
                 ii = i + isc -1
                 nb = Atm_block%blkno(ii,jj)
                 ix = Atm_block%ixp(ii,jj)
!         if(mpp_pe()==mpp_root_pe())print *,'in,fv3gfs_io,sze(Diag(idx)%data(nb)%var3)=',  &
!             size(Diag(idx)%data(nb)%var3,1),size(Diag(idx)%data(nb)%var3,2)
                 var3(i,j,k) = Diag(idx)%data(nb)%var3(ix,levo-k+1)*lcnvfac
               enddo
             enddo
           enddo
           call store_data3D(Diag(idx)%id, var3, Time, idx, Diag(idx)%intpl_method, Diag(idx)%name)
#ifdef JUNK
         else
           !--- temperature tendency
           if (trim(Diag(idx)%name) == 'dtemp_dt') then
             var3(1:nx,1:ny,1:levs) =  RESHAPE(Statein%tgrs(1:ngptc,levs:1:-1), (/nx,ny,levs/))
             var3(1:nx,1:ny,1:levs) = (RESHAPE(Stateout%gt0(1:ngptc,levs:1:-1), (/nx,ny,levs/))  &
                                        - var3(1:nx,1:ny,1:levs))*rdt
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1)
           endif
           !--- horizontal wind component tendency
           if (trim(Diag(idx)%name) == 'du_dt') then
             var3(1:nx,1:ny,1:levs) =  RESHAPE(Statein%ugrs(1:ngptc,levs:1:-1), (/nx,ny,levs/))
             var3(1:nx,1:ny,1:levs) = (RESHAPE(Stateout%gu0(1:ngptc,levs:1:-1), (/nx,ny,levs/))  &
                                        - var3(1:nx,1:ny,1:levs))*rdt
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1)
           endif
           !--- meridional wind component tendency
           if (trim(Diag(idx)%name) == 'dv_dt') then
             var3(1:nx,1:ny,1:levs) =  RESHAPE(Statein%vgrs(1:ngptc,levs:1:-1), (/nx,ny,levs/))
             var3(1:nx,1:ny,1:levs) = (RESHAPE(Stateout%gv0(1:ngptc,levs:1:-1), (/nx,ny,levs/))  &
                                        - var3(1:nx,1:ny,1:levs))*rdt
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1)
           endif
           !--- specific humidity tendency
           if (trim(Diag(idx)%name) == 'dsphum_dt') then
             var3(1:nx,1:ny,1:levs) =  RESHAPE(Statein%qgrs(1:ngptc,levs:1:-1,1:1), (/nx,ny,levs/))
             var3(1:nx,1:ny,1:levs) = (RESHAPE(Stateout%gq0(1:ngptc,levs:1:-1,1:1), (/nx,ny,levs/))  &
                                        - var3(1:nx,1:ny,1:levs))*rdt
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1)
           endif
           !--- cloud water mixing ration tendency
           if (trim(Diag(idx)%name) == 'dclwmr_dt') then
             var3(1:nx,1:ny,1:levs) =  RESHAPE(Statein%qgrs(1:ngptc,levs:1:-1,ntcw:ntcw), (/nx,ny,levs/))
             var3(1:nx,1:ny,1:levs) = (RESHAPE(Stateout%gq0(1:ngptc,levs:1:-1,ntcw:ntcw), (/nx,ny,levs/))  &
                                        - var3(1:nx,1:ny,1:levs))*rdt
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1)
           endif
           !--- ozone mixing ration tendency
#ifdef MULTI_GASES
           if (trim(Diag(idx)%name) == 'dspo3_dt') then
#else
           if (trim(Diag(idx)%name) == 'do3mr_dt') then
#endif
             var3(1:nx,1:ny,1:levs) =  RESHAPE(Statein%qgrs(1:ngptc,levs:1:-1,ntoz:ntoz), (/nx,ny,levs/))
             var3(1:nx,1:ny,1:levs) = (RESHAPE(Stateout%gq0(1:ngptc,levs:1:-1,ntoz:ntoz), (/nx,ny,levs/))  &
                                        - var3(1:nx,1:ny,1:levs))*rdt
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1)
           endif
#endif
         endif
       endif
     enddo


  end subroutine fv3gfs_diag_output
!
!-------------------------------------------------------------------------
  subroutine store_data(id, work, Time, idx, intpl_method, fldname)
    integer, intent(in)                 :: id
    integer, intent(in)                 :: idx
#ifdef CCPP_32BIT
    real, intent(in)                    :: work(:,:)
#else
    real(kind=kind_phys), intent(in)    :: work(ieco-isco+1,jeco-jsco+1)
#endif
    type(time_type), intent(in)         :: Time
    character(*), intent(in)            :: intpl_method
    character(*), intent(in)            :: fldname
!
    real(kind=kind_phys)                :: sinlat, sinlon, coslon
    integer k,j,i,kb,nv,i1,j1
    logical used
!
    if( id > 0 ) then
      if( use_wrtgridcomp_output ) then
        if( trim(intpl_method) == 'bilinear') then
!$omp parallel do default(shared) private(i,j)
          do j= jsco,jeco
            do i= isco,ieco
              buffer_phys_bl(i,j,nstt(idx)) = work(i-isco+1,j-jsco+1)
            enddo
          enddo
        else if(trim(intpl_method) == 'nearest_stod') then
!$omp parallel do default(shared) private(i,j)
          do j= jsco,jeco
            do i= isco,ieco
              buffer_phys_nb(i,j,nstt(idx)) = work(i-isco+1,j-jsco+1)
            enddo
          enddo
        else if(trim(intpl_method) == 'vector_bilinear') then
!first save the data
!$omp parallel do default(shared) private(i,j)
          do j= jsco,jeco
            do i= isco,ieco
              buffer_phys_bl(i,j,nstt(idx)) = work(i-isco+1,j-jsco+1)
            enddo
          enddo
          if( fldname(1:1) == 'u' .or. fldname(1:1) == 'U') then
            if(.not.allocated(uwork)) allocate(uwork(isco:ieco,jsco:jeco))
!$omp parallel do default(shared) private(i,j)
            do j= jsco,jeco
              do i= isco,ieco
                uwork(i,j) = work(i-isco+1,j-jsco+1)
              enddo
            enddo
            uwindname = fldname
            uwork_set = .true.
          endif
          if( fldname(1:1) == 'v' .or. fldname(1:1) == 'V') then
!set up wind vector
            if( uwork_set .and. trim(uwindname(2:)) == trim(fldname(2:))) then
              nv = nstt_vctbl(idx)
!$omp parallel do default(shared) private(i,j,i1,j1,sinlat,sinlon,coslon)
              do j= jsco,jeco
                j1 = j-jsco+1
                do i= isco,ieco
                  i1 = i-isco+1
                  sinlat = sin(lat(i,j))
                  sinlon = sin(lon(i,j))
                  coslon = cos(lon(i,j))
                  buffer_phys_windvect(1,i,j,nv) = uwork(i,j)*coslon - work(i1,j1)*sinlat*sinlon
                  buffer_phys_windvect(2,i,j,nv) = uwork(i,j)*sinlon + work(i1,j1)*sinlat*coslon
                  buffer_phys_windvect(3,i,j,nv) =                     work(i1,j1)*cos(lat(i,j))
                enddo
              enddo
            endif
            uwork     = zero
            uwindname = ''
            uwork_set = .false.
          endif

        endif
      else
        used = send_data(id, work, Time)
      endif
    endif
!
 end subroutine store_data
!
!-------------------------------------------------------------------------
!
  subroutine store_data3D(id, work, Time, idx, intpl_method, fldname)
    integer, intent(in)                 :: id
    integer, intent(in)                 :: idx
#ifdef CCPP_32BIT
    real, intent(in)                    :: work(:,:,:)
#else
    real(kind=kind_phys), intent(in)    :: work(ieco-isco+1,jeco-jsco+1,levo)
#endif
    type(time_type), intent(in)         :: Time
    character(*), intent(in)            :: intpl_method
    character(*), intent(in)            :: fldname
!
    real(kind=kind_phys), allocatable, dimension(:,:) :: sinlon, coslon, sinlat, coslat
    integer k,j,i,kb,nv,i1,j1
    logical used
!
    if( id > 0 ) then
      if( use_wrtgridcomp_output ) then
        if( trim(intpl_method) == 'bilinear') then
!$omp parallel do default(shared) private(i,j,k)
          do k= 1,levo
            do j= jsco,jeco
              do i= isco,ieco
                buffer_phys_bl(i,j,nstt(idx)+k-1) = work(i-isco+1,j-jsco+1,k)
              enddo
            enddo
          enddo
        else if(trim(intpl_method) == 'nearest_stod') then
!$omp parallel do default(shared) private(i,j,k)
          do k= 1,levo
            do j= jsco,jeco
              do i= isco,ieco
                buffer_phys_nb(i,j,nstt(idx)+k-1) = work(i-isco+1,j-jsco+1,k)
              enddo
            enddo
          enddo
        else if(trim(intpl_method) == 'vector_bilinear') then
!first save the data
!$omp parallel do default(shared) private(i,j,k)
          do k= 1,levo
            do j= jsco,jeco
              do i= isco,ieco
                buffer_phys_bl(i,j,nstt(idx)+k-1) = work(i-isco+1,j-jsco+1,k)
              enddo
            enddo
          enddo
          if( fldname(1:1) == 'u' .or. fldname(1:1) == 'U') then
            if(.not.allocated(uwork3d)) allocate(uwork3d(isco:ieco,jsco:jeco,levo))
!$omp parallel do default(shared) private(i,j,k)
            do k= 1, levo
              do j= jsco,jeco
                do i= isco,ieco
                  uwork3d(i,j,k) = work(i-isco+1,j-jsco+1,k)
                enddo
              enddo
            enddo
            uwindname = fldname
            uwork_set = .true.
          endif
          if( fldname(1:1) == 'v' .or. fldname(1:1) == 'V') then
!set up wind vector
            if( uwork_set .and. trim(uwindname(2:)) == trim(fldname(2:))) then
              allocate (sinlon(isco:ieco,jsco:jeco), coslon(isco:ieco,jsco:jeco), &
                        sinlat(isco:ieco,jsco:jeco), coslat(isco:ieco,jsco:jeco))
!$omp parallel do default(shared) private(i,j)
              do j= jsco,jeco
                do i= isco,ieco
                  sinlon(i,j) = sin(lon(i,j))
                  coslon(i,j) = cos(lon(i,j))
                  sinlat(i,j) = sin(lat(i,j))
                  coslat(i,j) = cos(lat(i,j))
                enddo
              enddo
!$omp parallel do default(shared) private(i,j,k,nv,i1,j1)
              do k= 1, levo
                nv = nstt_vctbl(idx)+k-1
                do j= jsco,jeco
                  j1 = j-jsco+1
                  do i= isco,ieco
                    i1 = i-isco+1
                    buffer_phys_windvect(1,i,j,nv) = uwork3d(i,j,k)*coslon(i,j) &
                                                   - work(i1,j1,k)*sinlat(i,j)*sinlon(i,j)
                    buffer_phys_windvect(2,i,j,nv) = uwork3d(i,j,k)*sinlon(i,j) &
                                                   + work(i1,j1,k)*sinlat(i,j)*coslon(i,j)
                    buffer_phys_windvect(3,i,j,nv) = work(i1,j1,k)*coslat(i,j)
                  enddo
                enddo
              enddo
              deallocate (sinlon, coslon, sinlat, coslat)
            endif
            uwork3d   = zero
            uwindname = ''
            uwork_set = .false.
          endif

        endif
      else
        used = send_data(id, work, Time)
      endif
    endif
!
 end subroutine store_data3D
!
!-------------------------------------------------------------------------
!
#ifdef use_WRTCOMP

 subroutine fv_phys_bundle_setup(Diag, axes, phys_bundle, fcst_grid, quilting, nbdlphys, rc)
!
!-------------------------------------------------------------
!*** set esmf bundle for phys output fields
!------------------------------------------------------------
!
   use esmf
   use diag_data_mod, ONLY:  diag_atttype
!
   implicit none
!
   type(GFS_externaldiag_type),intent(in)      :: Diag(:)
   integer, intent(in)                         :: axes(:)
   type(ESMF_FieldBundle),intent(inout)        :: phys_bundle(:)
   type(ESMF_Grid),intent(inout)               :: fcst_grid
   logical,intent(in)                          :: quilting
   integer, intent(in)                         :: nbdlphys
   integer,intent(out)                         :: rc

!
!*** local variables
   integer i, j, k, n, idx, ibdl, nbdl
   integer id, axis_length, direction, edges, axis_typ
   integer num_attributes, num_field_dyn
   integer currdate(6)
   character(2) axis_id
   character(255)    :: units, long_name, cart_name, axis_direct, edgesS
   character(128)    :: output_name, physbdl_name, outputfile1
   logical           :: lput2physbdl, loutputfile, l2dvector
   type(domain1d)    :: Domain
   type(domainUG)    :: DomainU
   type(ESMF_Field)  :: field
   real,dimension(:),allocatable               :: axis_data
   character(128),dimension(:), allocatable    :: bdl_intplmethod, outputfile
   type(diag_atttype),dimension(:),allocatable :: attributes
   real(4),dimension(:,:),pointer              :: dataPtr2d
!
   logical isPresent
   integer udimCount
   character(80),dimension(:),allocatable :: udimList
!
!------------------------------------------------------------
!--- use wrte grid component for output
   use_wrtgridcomp_output = quilting
!   if(mpp_pe()==mpp_root_pe())print *,'in fv_phys bundle,use_wrtgridcomp_output=',use_wrtgridcomp_output, &
!   print *,'in fv_phys bundle,use_wrtgridcomp_output=',use_wrtgridcomp_output, &
!       'isco=',isco,ieco,'jsco=',jsco,jeco,'tot_diag_idx=',tot_diag_idx
!
!------------------------------------------------------------
!*** add attributes to the bundle such as subdomain limtis,
!*** axes, output time, etc
!------------------------------------------------------------
!
   allocate(bdl_intplmethod(nbdlphys), outputfile(nbdlphys))
   if(mpp_pe()==mpp_root_pe()) print *,'in fv_phys bundle,nbdl=',nbdlphys
   do ibdl = 1, nbdlphys
     loutputfile = .false.
     call ESMF_FieldBundleGet(phys_bundle(ibdl), name=physbdl_name,rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
     idx = index(physbdl_name,'_bilinear')
     if(idx > 0) then
       outputfile(ibdl)      = physbdl_name(1:idx-1)
       bdl_intplmethod(ibdl) = 'bilinear'
       loutputfile           = .true.
     endif
     idx = index(physbdl_name,'_nearest_stod')
     if(idx > 0) then
       outputfile(ibdl)      = physbdl_name(1:idx-1)
       bdl_intplmethod(ibdl) = 'nearest_stod'
       loutputfile           = .true.
     endif
     if( .not. loutputfile) then
       outputfile(ibdl)      = 'phy'
       bdl_intplmethod(ibdl) = 'nearest_stod'
     endif
!    print *,'in fv_phys bundle,i=',ibdl,'outputfile=',trim(outputfile(ibdl)), &
!      'bdl_intplmethod=',trim(bdl_intplmethod(ibdl))

     call ESMF_AttributeAdd(phys_bundle(ibdl), convention="NetCDF", purpose="FV3", &
                            attrList=(/"fhzero     ", "ncld       ", "nsoil      ",&
                                       "imp_physics", "dtp        ", "landsfcmdl "/), rc=rc)

     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

     call ESMF_AttributeSet(phys_bundle(ibdl), convention="NetCDF", purpose="FV3", &
                            name="fhzero", value=fhzero, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

     call ESMF_AttributeSet(phys_bundle(ibdl), convention="NetCDF", purpose="FV3", &
                            name="ncld", value=ncld, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

     call ESMF_AttributeSet(phys_bundle(ibdl), convention="NetCDF", purpose="FV3", &
                            name="nsoil", value=nsoil, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

     call ESMF_AttributeSet(phys_bundle(ibdl), convention="NetCDF", purpose="FV3", &
                            name="imp_physics", value=imp_physics, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

     call ESMF_AttributeSet(phys_bundle(ibdl), convention="NetCDF", purpose="FV3", &
                            name="dtp", value=dtp, rc=rc)
!     print *,'in fcst gfdl diag, dtp=',dtp,' ibdl=',ibdl
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

     call ESMF_AttributeSet(phys_bundle(ibdl), convention="NetCDF", purpose="FV3", &
                            name="landsfcmdl", value=landsfcmdl, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!end ibdl
   enddo
!
!*** get axis names
   allocate(axis_name(num_axes_phys))
   do id = 1,num_axes_phys
     call get_diag_axis_name( axes(id), axis_name(id))
   enddo
   isPresent = .false.
   if( num_axes_phys>2 ) then
     allocate(axis_name_vert(num_axes_phys-2))
     do id=3,num_axes_phys
       axis_name_vert(id-2) = axis_name(id)
     enddo
!
     call ESMF_AttributeGet(fcst_grid, convention="NetCDF", purpose="FV3", &
                            name="vertical_dim_labels", isPresent=isPresent, &
                            itemCount=udimCount, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

     if (isPresent .and. (udimCount>num_axes_phys-2) ) then
       allocate(udimList(udimCount))
       call ESMF_AttributeGet(fcst_grid, convention="NetCDF", purpose="FV3", &
                              name="vertical_dim_labels", valueList=udimList, rc=rc)
!       if(mpp_pe()==mpp_root_pe()) print *,'in fv3gfsio, vertical
!       list=',udimList(1:udimCount),'rc=',rc

       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

     else

       if(mpp_pe()==mpp_root_pe()) print *,'in fv_dyn bundle,axis_name_vert=',axis_name_vert
       call ESMF_AttributeAdd(fcst_grid, convention="NetCDF", purpose="FV3",  &
                              attrList=(/"vertical_dim_labels"/), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
                              name="vertical_dim_labels", valueList=axis_name_vert, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
     endif
     deallocate(axis_name_vert)
   endif

!*** add attributes
   if(allocated(all_axes)) deallocate(all_axes)
   allocate(all_axes(num_axes_phys))
   all_axes(1:num_axes_phys) = axes(1:num_axes_phys)
   if (.not. isPresent .or. (udimCount<num_axes_phys-2) ) then
     do id = 1,num_axes_phys
       axis_length =  get_axis_global_length(axes(id))
       allocate(axis_data(axis_length))
       call get_diag_axis( axes(id), axis_name(id), units, long_name, cart_name, &
                         direction, edges, Domain, DomainU, axis_data,           &
                         num_attributes=num_attributes, attributes=attributes)
!
       edgesS = ''
       do i = 1,num_axes_phys
         if(axes(i) == edges) edgesS=axis_name(i)
       enddo
! Add vertical dimension Attributes to Grid
       if( id>2 ) then
!      if(mpp_pe()==mpp_root_pe()) print *,' in dyn add grid, axis_name=',     &
!         trim(axis_name(id)),'axis_data=',axis_data
         if(trim(edgesS)/='') then
           call ESMF_AttributeAdd(fcst_grid, convention="NetCDF", purpose="FV3",  &
             attrList=(/trim(axis_name(id)),trim(axis_name(id))//":long_name",    &
                    trim(axis_name(id))//":units", trim(axis_name(id))//":cartesian_axis", &
                    trim(axis_name(id))//":positive", trim(axis_name(id))//":edges"/), rc=rc)
         else
           call ESMF_AttributeAdd(fcst_grid, convention="NetCDF", purpose="FV3",  &
             attrList=(/trim(axis_name(id)),trim(axis_name(id))//":long_name",    &
                    trim(axis_name(id))//":units", trim(axis_name(id))//":cartesian_axis", &
                    trim(axis_name(id))//":positive"/), rc=rc)
         endif
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

         call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
                                name=trim(axis_name(id)), valueList=axis_data, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

         call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
                                name=trim(axis_name(id))//":long_name", value=trim(long_name), rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

         call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
                                name=trim(axis_name(id))//":units", value=trim(units), rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

         call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
                                name=trim(axis_name(id))//":cartesian_axis", value=trim(cart_name), rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

         if(direction > 0) then
           axis_direct = "up"
         else
           axis_direct = "down"
         endif
         call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
                                name=trim(axis_name(id))//":positive", value=trim(axis_direct), rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

         if(trim(edgesS)/='') then
           call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
                                  name=trim(axis_name(id))//":edges", value=trim(edgesS), rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
         endif

        endif
!
        deallocate(axis_data)
      enddo
    endif
!   print *,'in setup fieldbundle_phys, num_axes_phys=',num_axes_phys,'tot_diag_idx=',tot_diag_idx, &
!       'nbdlphys=',nbdlphys
!
!-----------------------------------------------------------------------------------------
!*** add esmf fields
!
   do idx= 1,tot_diag_idx

     lput2physbdl = .false.
     do ibdl = 1, nbdlphys

       if( index(trim(Diag(idx)%intpl_method),trim(bdl_intplmethod(ibdl))) > 0) then
         lput2physbdl = .true.
         if( Diag(idx)%id > 0 ) then
           call find_output_name(trim(Diag(idx)%mod_name),trim(Diag(idx)%name),output_name)

!add origin field
           call add_field_to_phybundle(trim(output_name),trim(Diag(idx)%desc),trim(Diag(idx)%unit), "time: point",         &
                                       axes(1:Diag(idx)%axes), fcst_grid, nstt(idx), phys_bundle(ibdl), outputfile(ibdl),  &
                                       bdl_intplmethod(ibdl), rcd=rc)
!           if( mpp_pe() == mpp_root_pe()) print *,'phys, add field,',trim(Diag(idx)%name),'idx=',idx,'ibdl=',ibdl
!
           if( index(trim(Diag(idx)%intpl_method), "vector") > 0) then
             l2dvector = .true.
             if (nstt_vctbl(idx) > 0) then
               output_name = 'wind'//trim(output_name)//'vector'
               outputfile1 = 'none'
               call add_field_to_phybundle(trim(output_name),trim(Diag(idx)%desc),trim(Diag(idx)%unit), "time: point",       &
                                          axes(1:Diag(idx)%axes), fcst_grid, nstt_vctbl(idx),phys_bundle(ibdl), outputfile1, &
                                          bdl_intplmethod(ibdl),l2dvector=l2dvector,  rcd=rc)
!               if( mpp_pe() == mpp_root_pe()) print *,'in phys, add vector field,',trim(Diag(idx)%name),' idx=',idx,' ibdl=',ibdl
             endif
           endif

         endif
       endif
     enddo
     if( .not. lput2physbdl ) then
         if( mpp_pe() == mpp_root_pe()) print *,'WARNING: not matching interpolation method, field ',trim(Diag(idx)%name), &
           ' is not added to phys bundle '
     endif

   enddo
   deallocate(axis_name)
   deallocate(all_axes)

 end subroutine fv_phys_bundle_setup
!
!-----------------------------------------------------------------------------------------
 subroutine add_field_to_phybundle(var_name,long_name,units,cell_methods, axes,phys_grid, &
                                   kstt,phys_bundle,output_file,intpl_method,range,l2dvector,rcd)
!
   use esmf
!
   implicit none

   character(*), intent(in)             :: var_name, long_name, units, cell_methods
   character(*), intent(in)             :: output_file, intpl_method
   integer, intent(in)                  :: axes(:)
   type(esmf_grid), intent(in)          :: phys_grid
   integer, intent(in)                  :: kstt
   type(esmf_fieldbundle),intent(inout) :: phys_bundle
   real, intent(in), optional           :: range(2)
   logical, intent(in), optional        :: l2dvector
   integer, intent(out), optional       :: rcd
!
!*** local variable
   type(ESMF_Field)         :: field
   type(ESMF_DataCopy_Flag) :: copyflag=ESMF_DATACOPY_REFERENCE
   integer rc, i, j, idx
   real(4),dimension(:,:),pointer   :: temp_r2d
   real(4),dimension(:,:,:),pointer :: temp_r3d
   logical :: l2dvector_local
!
   ! fix for non-standard compilers (e.g. PGI)
   l2dvector_local = .false.
   if (present(l2dvector)) then
     if (l2dvector) then
         l2dvector_local = .true.
     end if
   end if
!
!*** create esmf field
   if (l2dvector_local .and. size(axes)==2) then
     temp_r3d => buffer_phys_windvect(1:3,isco:ieco,jsco:jeco,kstt)
!     if( mpp_root_pe() == 0) print *,'phys, create wind vector esmf field'
     call ESMF_LogWrite('bf create winde vector esmf field '//trim(var_name), ESMF_LOGMSG_INFO, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!datacopyflag=ESMF_DATACOPY_VALUE, &
     field = ESMF_FieldCreate(phys_grid, temp_r3d, datacopyflag=ESMF_DATACOPY_REFERENCE,          &
                            gridToFieldMap=(/2,3/), ungriddedLBound=(/1/), ungriddedUBound=(/3/), &
                            name=var_name, indexFlag=ESMF_INDEX_DELOCAL, rc=rc)

     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
     call ESMF_LogWrite('af winde vector esmf field create '//trim(var_name), ESMF_LOGMSG_INFO, rc=rc)

     call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
                            attrList=(/"output_file"/), rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
       line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

     call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
                            name='output_file',value=trim(output_file),rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
       line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

     call ESMF_LogWrite('before winde vector esmf field add output_file', ESMF_LOGMSG_INFO, rc=rc)

!     if( mpp_root_pe() == 0)print *,'phys, aftercreate wind vector esmf field'
     call ESMF_FieldBundleAdd(phys_bundle,(/field/), rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

     if( present(rcd)) rcd=rc
     call ESMF_LogWrite('aft winde vector esmf field add to fieldbundle'//trim(var_name), ESMF_LOGMSG_INFO, rc=rc)
     return
   else if( trim(intpl_method) == 'nearest_stod' ) then
     if(size(axes) == 2) then
       temp_r2d => buffer_phys_nb(isco:ieco,jsco:jeco,kstt)
       field = ESMF_FieldCreate(phys_grid, temp_r2d, datacopyflag=copyflag, &
                              name=var_name, indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,       &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

     else if(size(axes) == 3) then
       temp_r3d => buffer_phys_nb(isco:ieco,jsco:jeco,kstt:kstt+levo-1)
       field = ESMF_FieldCreate(phys_grid, temp_r3d, datacopyflag=copyflag, &
                              name=var_name, indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,       &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

     if( mpp_pe() == mpp_root_pe()) print *,'add 3D field to after nearest_stod, fld=', trim(var_name)
     endif
   else if( trim(intpl_method) == 'bilinear' ) then
     if(size(axes) == 2) then
       temp_r2d => buffer_phys_bl(isco:ieco,jsco:jeco,kstt)
       field = ESMF_FieldCreate(phys_grid, temp_r2d, datacopyflag=copyflag, &
                            name=var_name, indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,       &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
     else if(size(axes) == 3) then
       temp_r3d => buffer_phys_bl(isco:ieco,jsco:jeco,kstt:kstt+levo-1)
       field = ESMF_FieldCreate(phys_grid, temp_r3d, datacopyflag=copyflag, &
                            name=var_name, indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,       &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
       if( mpp_pe() == mpp_root_pe()) print *,'add field to after bilinear, fld=', trim(var_name)
     endif
   endif
!
!*** add field attributes
   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
                          attrList=(/"long_name"/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
                          name='long_name',value=trim(long_name),rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
                          attrList=(/"units"/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
                          name='units',value=trim(units),rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
                          attrList=(/"missing_value"/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
                          name='missing_value',value=real(missing_value,kind=4),rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
                          attrList=(/"_FillValue"/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
                          name='_FillValue',value=real(missing_value,kind=4),rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
                          attrList=(/"cell_methods"/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
                          name='cell_methods',value=trim(cell_methods),rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
                          attrList=(/"output_file"/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     line=__LINE__, file=__FILE__))  call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
                          name='output_file',value=trim(output_file),rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!
!*** add vertical coord attribute:
   if( size(axes) > 2) then
     do i=3,size(axes)
       idx=0
       do j=1,size(all_axes)
         if (axes(i)==all_axes(j)) then
           idx=j
           exit
         endif
       enddo
       if (idx>0) then
         call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
                                attrList=(/"ESMF:ungridded_dim_labels"/), rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
         call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
                                name="ESMF:ungridded_dim_labels", valueList=(/trim(axis_name(idx))/), rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       endif
     enddo
   endif

!*** add field into bundle
   call ESMF_FieldBundleAdd(phys_bundle,(/field/), rc=rc)
   if( present(rcd)) rcd=rc
!
   call ESMF_LogWrite('phys field add to fieldbundle '//trim(var_name), ESMF_LOGMSG_INFO, rc=rc)

 end subroutine add_field_to_phybundle
!
!
 subroutine find_output_name(module_name,field_name,output_name)
   character(*), intent(in)     :: module_name
   character(*), intent(in)     :: field_name
   character(*), intent(out)    :: output_name
!
   integer i,in_num, out_num
   integer tile_count
!
   tile_count = 1
   in_num = find_input_field(module_name, field_name, tile_count)
!
   output_name = ''
   do i=1, max_output_fields
     if(output_fields(i)%input_field == in_num) then
       output_name = output_fields(i)%output_name
       exit
     endif
   enddo
   if(output_name == '') then
     print *,'Error, cant find out put name'
   endif

 end subroutine find_output_name
#endif
!-------------------------------------------------------------------------

end module FV3GFS_io_mod
