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

!
!--- GFS_typedefs
  use GFS_typedefs,       only: GFS_sfcprop_type, GFS_control_type, &
                                GFS_data_type, kind_phys
  use GFS_restart,        only: GFS_restart_type
  use GFS_diagnostics,    only: GFS_externaldiag_type

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
  character(len=32),    allocatable,         dimension(:)       :: oro_name2, sfc_name2, sfc_name3
  real(kind=kind_phys), allocatable, target, dimension(:,:,:)   :: oro_var2, sfc_var2, phy_var2, sfc_var3ice
  character(len=32),    allocatable,         dimension(:)       :: oro_ls_ss_name
  real(kind=kind_phys), allocatable, target, dimension(:,:,:)   :: oro_ls_var, oro_ss_var, oro_var3v, oro_var3s
  real(kind=kind_phys), allocatable, target, dimension(:,:,:,:) :: sfc_var3, phy_var3
  character(len=32),    allocatable,         dimension(:)       :: dust12m_name, emi_name, rrfssd_name
  real(kind=kind_phys), allocatable, target, dimension(:,:,:,:) :: rrfssd_var
  real(kind=kind_phys), allocatable, target, dimension(:,:,:,:) :: dust12m_var
  real(kind=kind_phys), allocatable, target, dimension(:,:,:,:) :: emi_var
  !--- Noah MP restart containers
  real(kind=kind_phys), allocatable, target, dimension(:,:,:,:) :: sfc_var3sn,sfc_var3eq,sfc_var3zn

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
  integer, parameter, public :: DIAG_SIZE = 500
  real, parameter :: missing_value = 9.99e20_r8
  real, parameter:: stndrd_atmos_ps = 101325.0_r8
  real, parameter:: stndrd_atmos_lapse = 0.0065_r8
  real, parameter:: drythresh = 1.e-4_r8, zero = 0.0_r8, one = 1.0_r8
  real, parameter:: min_lake_orog = 200.0_r8
  real(kind=kind_phys), parameter :: timin = 173.0_r8  ! minimum temperature allowed for snow/ice

!--- miscellaneous other variables
  logical :: use_wrtgridcomp_output = .FALSE.
  logical :: module_is_initialized  = .FALSE.

  type rrfs_sd_data_type
    ! The smoke_data_type stores temporary arrays used to read or
    ! write RRFS-SD restart and axis variables.

    real(kind_phys), pointer, private, dimension(:,:) :: & ! i,j variables
         emdust=>null(), emseas=>null(), emanoc=>null(), fhist=>null(), coef_bb_dc=>null()

    real(kind_phys), pointer, private, dimension(:,:,:) :: &
         fire_in=>null() ! i, j, fire_aux_data_levels

  contains
    procedure, public :: register_axis => rrfs_sd_register_axis ! register fire_aux_data_levels axis
    procedure, public :: write_axis => rrfs_sd_write_axis ! write fire_aux_data_levels variable
    procedure, public :: allocate_data => rrfs_sd_allocate_data ! allocate all pointers
    procedure, public :: fill_data => rrfs_sd_fill_data ! fill data with default values
    procedure, public :: register_fields => rrfs_sd_register_fields ! register rrfs_sd fields
    procedure, public :: deallocate_data => rrfs_sd_deallocate_data ! deallocate pointers
    procedure, public :: copy_to_temporaries => rrfs_sd_copy_to_temporaries ! Copy Sfcprop to arrays
    procedure, public :: copy_from_temporaries => rrfs_sd_copy_from_temporaries ! Copy arrays to Sfcprop
    final :: rrfs_sd_final ! Destructor; calls deallocate_data
  end type rrfs_sd_data_type

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
      endif

      if(Model%cplwav) then
        nt=nt+1 ; sfc_name2(nvar_s2m) = 'zorlwav' !zorl from wave component
      endif

      if (Model%nstf_name(1) > 0) then
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
      endif
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
    integer :: nvar_o2, nvar_s2m, nvar_s2o, nvar_s3
    integer :: nvar_oro_ls_ss
    integer :: nvar_vegfr, nvar_soilfr
    integer :: nvar_s2r, nvar_s2mp, nvar_s3mp, isnow
    integer :: nvar_emi, nvar_dust12m, nvar_rrfssd
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
    logical :: is_lsoil

    type(rrfs_sd_data_type) :: rrfs_sd_data

    nvar_o2  = 19
    nvar_oro_ls_ss = 10

    nvar_vegfr  = Model%nvegcat
    nvar_soilfr = Model%nsoilcat

    if (Model%nstf_name(1) > 0) then
      nvar_s2o = 18
    else
      nvar_s2o = 0
    endif
    if(Model%rrfs_sd) then
      nvar_dust12m = 5
      nvar_rrfssd  = 3
      nvar_emi     = 1
    else
      nvar_dust12m = 0
      nvar_rrfssd  = 0
      nvar_emi     = 0
    endif

    if (Model%lsm == Model%lsm_ruc .and. warm_start) then
      if(Model%rdlai) then
        nvar_s2r = 13
      else
        nvar_s2r = 12
      end if
      nvar_s3  = 5
    else
      if(Model%rdlai) then
       nvar_s2r = 1
      else
       nvar_s2r = 0
      endif
      nvar_s3  = 3
    endif

    if (Model%lsm == Model%lsm_noahmp) then
      nvar_s2mp = 29       !mp 2D
      nvar_s3mp = 5        !mp 3D
    else
      nvar_s2mp = 0        !mp 2D
      nvar_s3mp = 0        !mp 3D
    endif

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

    if (.not. allocated(oro_name2)) then
    !--- allocate the various containers needed for orography data
      allocate(oro_name2(nvar_o2))
      allocate(oro_var2(nx,ny,nvar_o2))

      allocate(oro_var3v(nx,ny,nvar_vegfr))
      allocate(oro_var3s(nx,ny,nvar_soilfr))

      oro_var2 = -9999._kind_phys

      num = 1       ; oro_name2(num)  = 'stddev'     ! hprime(ix,1)
      num = num + 1 ; oro_name2(num)  = 'convexity'  ! hprime(ix,2)
      num = num + 1 ; oro_name2(num)  = 'oa1'        ! hprime(ix,3)
      num = num + 1 ; oro_name2(num)  = 'oa2'        ! hprime(ix,4)
      num = num + 1 ; oro_name2(num)  = 'oa3'        ! hprime(ix,5)
      num = num + 1 ; oro_name2(num)  = 'oa4'        ! hprime(ix,6)
      num = num + 1 ; oro_name2(num)  = 'ol1'        ! hprime(ix,7)
      num = num + 1 ; oro_name2(num)  = 'ol2'        ! hprime(ix,8)
      num = num + 1 ; oro_name2(num)  = 'ol3'        ! hprime(ix,9)
      num = num + 1 ; oro_name2(num) = 'ol4'        ! hprime(ix,10)
      num = num + 1 ; oro_name2(num) = 'theta'      ! hprime(ix,11)
      num = num + 1 ; oro_name2(num) = 'gamma'      ! hprime(ix,12)
      num = num + 1 ; oro_name2(num) = 'sigma'      ! hprime(ix,13)
      num = num + 1 ; oro_name2(num) = 'elvmax'     ! hprime(ix,14)
      num = num + 1 ; oro_name2(num) = 'orog_filt'  ! oro
      num = num + 1 ; oro_name2(num) = 'orog_raw'   ! oro_uf
      num = num + 1 ; oro_name2(num) = 'land_frac'  ! land fraction [0:1]
      !--- variables below here are optional
      num = num + 1 ; oro_name2(num) = 'lake_frac'  ! lake fraction [0:1]
      num = num + 1 ; oro_name2(num) = 'lake_depth' ! lake depth(m)

      !--- register axis
      call register_axis( Oro_restart, "lon", 'X' )
      call register_axis( Oro_restart, "lat", 'Y' )
      !--- register the 2D fields
      do n = 1,num
         var2_p => oro_var2(:,:,n)
         if (trim(oro_name2(n)) == 'lake_frac' .or. trim(oro_name2(n)) == 'lake_depth' ) then
            call register_restart_field(Oro_restart, oro_name2(n), var2_p, dimensions=(/'lat','lon'/), is_optional=.true.)
         else
            call register_restart_field(Oro_restart, oro_name2(n), var2_p, dimensions=(/'lat','lon'/))
         endif
      enddo
      nullify(var2_p)

     !--- register 3D vegetation and soil fractions
      var3_fr => oro_var3v(:,:,:)
      call register_restart_field(Oro_restart, 'vegetation_type_pct', var3_fr, dimensions=(/'num_veg_cat','lat        ','lon        '/) , is_optional=.true.)
      var3_fr => oro_var3s(:,:,:)
      call register_restart_field(Oro_restart, 'soil_type_pct', var3_fr, dimensions=(/'num_soil_cat','lat         ','lon         '/) , is_optional=.true.)
      nullify(var3_fr)

   endif

   !--- read the orography restart/data
   call mpp_error(NOTE,'reading topographic/orographic information from INPUT/oro_data.tile*.nc')
   call read_restart(Oro_restart, ignore_checksum=ignore_rst_cksum)
   call close_file(Oro_restart)


   !--- copy data into GFS containers

!$omp parallel do default(shared) private(i, j, nb, ix, num)
    do nb = 1, Atm_block%nblks
      !--- 2D variables
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1
        !--- stddev
!       Sfcprop(nb)%hprim(ix)     = oro_var2(i,j,1)
        !--- hprime(1:14)
        num = 1       ; Sfcprop(nb)%hprime(ix,num)  = oro_var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num)  = oro_var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num)  = oro_var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num)  = oro_var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num)  = oro_var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num)  = oro_var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num)  = oro_var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num)  = oro_var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num)  = oro_var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num) = oro_var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num) = oro_var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num) = oro_var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num) = oro_var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num) = oro_var2(i,j,num)
        !--- oro
        num = num + 1 ; Sfcprop(nb)%oro(ix)       = oro_var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%oro_uf(ix)    = oro_var2(i,j,num)

        Sfcprop(nb)%landfrac(ix)  = -9999.0
        Sfcprop(nb)%lakefrac(ix)  = -9999.0

        num = num + 1 ; Sfcprop(nb)%landfrac(ix)  = oro_var2(i,j,num) !land frac [0:1]
        num = num + 1 ; Sfcprop(nb)%lakefrac(ix)  = oro_var2(i,j,num) !lake frac [0:1]
        num = num + 1 ; Sfcprop(nb)%lakedepth(ix) = oro_var2(i,j,num) !lake depth [m]    !YWu

        Sfcprop(nb)%vegtype_frac(ix,:)  =  -9999.0
        Sfcprop(nb)%soiltype_frac(ix,:) =  -9999.0

        Sfcprop(nb)%vegtype_frac(ix,:)  = oro_var3v(i,j,:) ! vegetation type fractions, [0:1]
        Sfcprop(nb)%soiltype_frac(ix,:) = oro_var3s(i,j,:) ! soil type fractions, [0:1]

        !do n=1,nvar_vegfr
        !  if (Sfcprop(nb)%vegtype_frac(ix,n) > 0.) print *,'Sfcprop(nb)%vegtype_frac(ix,n)',Sfcprop(nb)%vegtype_frac(ix,n),n
        !enddo
        !do n=1,nvar_soilfr
        !  if (Sfcprop(nb)%soiltype_frac(ix,n) > 0.) print *,'Sfcprop(nb)%soiltype_frac(ix,n)',Sfcprop(nb)%soiltype_frac(ix,n),n
        !enddo

      enddo
    enddo

    nvar_s2m = 48
    if (Model%use_cice_alb .or. Model%lsm == Model%lsm_ruc) then
      nvar_s2m = nvar_s2m + 4
!     nvar_s2m = nvar_s2m + 5
    endif
    if (Model%cplwav) then
      nvar_s2m = nvar_s2m + 1
    endif

    !--- deallocate containers and free restart container
    deallocate(oro_name2, oro_var2)
    deallocate(oro_var3v)
    deallocate(oro_var3s)

    if_smoke: if(Model%rrfs_sd) then  ! for RRFS-SD

    !--- Dust input FILE
    !--- open file
    infile=trim(indir)//'/'//trim(fn_dust12m)
    amiopen=open_file(dust12m_restart, trim(infile), 'read', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
    if (.not.amiopen) call mpp_error( FATAL, 'Error with opening file'//trim(infile) )

    if (.not. allocated(dust12m_name)) then
    !--- allocate the various containers needed for fengsha dust12m data
      allocate(dust12m_name(nvar_dust12m))
      allocate(dust12m_var(nx,ny,12,nvar_dust12m))

      dust12m_name(1)  = 'clay'
      dust12m_name(2)  = 'rdrag'
      dust12m_name(3)  = 'sand'
      dust12m_name(4)  = 'ssm'
      dust12m_name(5)  = 'uthr'

      !--- register axis
      call register_axis(dust12m_restart, 'lon', 'X')
      call register_axis(dust12m_restart, 'lat', 'Y')
      call register_axis(dust12m_restart, 'time', 12)
      !--- register the 3D fields
      do num = 1,nvar_dust12m
        var3_p2 => dust12m_var(:,:,:,num)
        call register_restart_field(dust12m_restart, dust12m_name(num), var3_p2, dimensions=(/'time', 'lat ', 'lon '/),&
                                  &is_optional=.not.mand)
      enddo
      nullify(var3_p2)
    endif

    !--- read new GSL created dust12m restart/data
    call mpp_error(NOTE,'reading dust12m information from INPUT/dust12m_data.tile*.nc')
    call read_restart(dust12m_restart)
    call close_file(dust12m_restart)

    do nb = 1, Atm_block%nblks
      !--- 3D variables
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1
        do k = 1, 12
          Sfcprop(nb)%dust12m_in(ix,k,1)  = dust12m_var(i,j,k,1)
          Sfcprop(nb)%dust12m_in(ix,k,2)  = dust12m_var(i,j,k,2)
          Sfcprop(nb)%dust12m_in(ix,k,3)  = dust12m_var(i,j,k,3)
          Sfcprop(nb)%dust12m_in(ix,k,4)  = dust12m_var(i,j,k,4)
          Sfcprop(nb)%dust12m_in(ix,k,5)  = dust12m_var(i,j,k,5)
        enddo
      enddo
    enddo

    deallocate(dust12m_name,dust12m_var)

    read_emi: if(nvar_emi>0) then
    !--- open anthropogenic emission file
    infile=trim(indir)//'/'//trim(fn_emi)
    amiopen=open_file(emi_restart, trim(infile), 'read', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
    if (.not.amiopen) call mpp_error( FATAL, 'Error with opening file'//trim(infile) )

    !if (.not. allocated(emi_name)) then
     !--- allocate the various containers needed for anthropogenic emission data
      if(allocated(emi_name)) deallocate(emi_name)
      if(allocated(emi_var)) deallocate(emi_var)
      allocate(emi_name(nvar_emi))
      allocate(emi_var(nx,ny,1,nvar_emi))

      emi_name(1)  = 'e_oc'
      !--- register axis
      call register_axis( emi_restart, 'time', 1) ! only read first time level, even if multiple are present
      call register_axis( emi_restart, "grid_xt", 'X' )
      call register_axis( emi_restart, "grid_yt", 'Y' )
      !--- register the 2D fields
      do num = 1,nvar_emi
        var3_p2 => emi_var(:,:,:,num)
        call register_restart_field(emi_restart, emi_name(num), var3_p2, dimensions=(/'time   ','grid_yt','grid_xt'/))
      enddo
      nullify(var3_p2)
    !endif

    !--- read anthropogenic emi restart/data
    call mpp_error(NOTE,'reading emi information from INPUT/emi_data.tile*.nc')
    call read_restart(emi_restart)
    call close_file(emi_restart)

    do num=1,nvar_emi
    do nb = 1, Atm_block%nblks
      !--- 2D variables
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1
        Sfcprop(nb)%emi_in(ix,num)  = emi_var(i,j,1,num)
      enddo
    enddo
    enddo

    !--- deallocate containers and free restart container
    deallocate(emi_name, emi_var)
    endif read_emi

    !--- Dust input FILE
    !--- open file
    infile=trim(indir)//'/'//trim(fn_rrfssd)
    amiopen=open_file(rrfssd_restart, trim(infile), 'read', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
    if (.not.amiopen) call mpp_error( FATAL, 'Error with opening file'//trim(infile) )

    if (.not. allocated(rrfssd_name)) then
      !--- allocate the various containers needed for rrfssd fire data
      allocate(rrfssd_name(nvar_rrfssd))
      allocate(rrfssd_var(nx,ny,24,nvar_rrfssd))

      rrfssd_name(1)  = 'ebb_smoke_hr'
      rrfssd_name(2)  = 'frp_avg_hr'
      rrfssd_name(3)  = 'frp_std_hr'

      !--- register axis
      call register_axis(rrfssd_restart, 'lon', 'X')
      call register_axis(rrfssd_restart, 'lat', 'Y')
      call register_axis(rrfssd_restart, 't', 24)
      !--- register the 3D fields
      mand = .false.
      do num = 1,nvar_rrfssd
       var3_p2 => rrfssd_var(:,:,:,num)
       call register_restart_field(rrfssd_restart, rrfssd_name(num), var3_p2, dimensions=(/'t  ', 'lat', 'lon'/),&
                                  &is_optional=.not.mand)
      enddo
      nullify(var3_p2)
    endif

    !--- read new GSL created rrfssd restart/data
    call mpp_error(NOTE,'reading rrfssd information from INPUT/SMOKE_RRFS_data.nc')
    call read_restart(rrfssd_restart)
    call close_file(rrfssd_restart)

    do nb = 1, Atm_block%nblks
      !--- 3D variables
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1
        !--- assign hprime(1:10) and hprime(15:24) with new oro stat data
        do k = 1, 24
          Sfcprop(nb)%smoke_RRFS(ix,k,1)  = rrfssd_var(i,j,k,1)
          Sfcprop(nb)%smoke_RRFS(ix,k,2)  = rrfssd_var(i,j,k,2)
          Sfcprop(nb)%smoke_RRFS(ix,k,3)  = rrfssd_var(i,j,k,3)
        enddo
      enddo
    enddo

    deallocate(rrfssd_name, rrfssd_var)
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

    if (.not. allocated(sfc_name2)) then
      !--- allocate the various containers needed for restarts
      allocate(sfc_name2(nvar_s2m+nvar_s2o+nvar_s2mp+nvar_s2r))
      allocate(sfc_name3(0:nvar_s3+nvar_s3mp))
      allocate(sfc_var2(nx,ny,nvar_s2m+nvar_s2o+nvar_s2mp+nvar_s2r))
      ! Note that this may cause problems with RUC LSM for coldstart runs from GFS data
      ! if the initial conditions do contain this variable, because Model%kice is 9 for
      ! RUC LSM, but tiice in the initial conditions will only have two vertical layers
      allocate(sfc_var3ice(nx,ny,Model%kice))

      if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp .or. (.not.warm_start)) then
        allocate(sfc_var3(nx,ny,Model%lsoil,nvar_s3))
      else if (Model%lsm == Model%lsm_ruc) then
        allocate(sfc_var3(nx,ny,Model%lsoil_lsm,nvar_s3))
      end if

      sfc_var2   = -9999.0_r8
      sfc_var3   = -9999.0_r8
      sfc_var3ice= -9999.0_r8
!
      if (Model%lsm == Model%lsm_noahmp) then
        allocate(sfc_var3sn(nx,ny,-2:0,4:6))
        allocate(sfc_var3eq(nx,ny,1:4,7:7))
        allocate(sfc_var3zn(nx,ny,-2:4,8:8))
        sfc_var3sn = -9999.0_r8
        sfc_var3eq = -9999.0_r8
        sfc_var3zn = -9999.0_r8
      end if

      call fill_Sfcprop_names(Model,sfc_name2,sfc_name3,nvar_s2m,warm_start)

      is_lsoil=.false.
      if ( .not. warm_start ) then
        if( variable_exists(Sfc_restart,"lsoil") ) then
          is_lsoil=.true.
          call register_axis(Sfc_restart, 'lon', 'X')
          call register_axis(Sfc_restart, 'lat', 'Y')
          call register_axis(Sfc_restart, 'lsoil', dimension_length=Model%lsoil)
       else
          call register_axis(Sfc_restart, 'xaxis_1', 'X')
          call register_axis(Sfc_restart, 'yaxis_1', 'Y')
          call register_axis(Sfc_restart, 'zaxis_1', dimension_length=4)
          call register_axis(Sfc_restart, 'Time', 1)
        end if
      else
        call register_axis(Sfc_restart, 'xaxis_1', 'X')
        call register_axis(Sfc_restart, 'yaxis_1', 'Y')
        call register_axis(Sfc_restart, 'zaxis_1', dimension_length=Model%kice)

        if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp) then
          call register_axis(Sfc_restart, 'zaxis_2', dimension_length=Model%lsoil)
        else if(Model%lsm == Model%lsm_ruc) then
          call register_axis(Sfc_restart, 'zaxis_2', dimension_length=Model%lsoil_lsm)
        end if
        if(Model%lsm == Model%lsm_noahmp) then
          call register_axis(Sfc_restart, 'zaxis_3', dimension_length=3)
          call register_axis(Sfc_restart, 'zaxis_4', dimension_length=7)
        end if
        call register_axis(Sfc_restart, 'Time', unlimited)
      end if

      if(Model%rrfs_sd) then
        call rrfs_sd_data%allocate_data(Model)
        call rrfs_sd_data%fill_data(Model, Sfcprop, Atm_block)
        call rrfs_sd_data%register_axis(Model)
        call rrfs_sd_data%register_fields
      endif

      !--- register the 2D fields
      do num = 1,nvar_s2m
        var2_p => sfc_var2(:,:,num)
        if (trim(sfc_name2(num)) == 'sncovr'.or. trim(sfc_name2(num)) == 'tsfcl' .or. trim(sfc_name2(num)) == 'zorll'   &
                                            .or. trim(sfc_name2(num)) == 'zorli' .or. trim(sfc_name2(num)) == 'zorlwav' &
                                            .or. trim(sfc_name2(num)) == 'snodl' .or. trim(sfc_name2(num)) == 'weasdl'  &
                                            .or. trim(sfc_name2(num)) == 'snodi' .or. trim(sfc_name2(num)) == 'weasdi'  &
                                            .or. trim(sfc_name2(num)) == 'tsfc'  .or. trim(sfc_name2(num)) ==  'zorlw'  &
                                            .or. trim(sfc_name2(num)) == 'albdirvis_lnd' .or. trim(sfc_name2(num)) == 'albdirnir_lnd' &
                                            .or. trim(sfc_name2(num)) == 'albdifvis_lnd' .or. trim(sfc_name2(num)) == 'albdifnir_lnd' &
                                            .or. trim(sfc_name2(num)) == 'albdirvis_ice' .or. trim(sfc_name2(num)) == 'albdirnir_ice' &
                                            .or. trim(sfc_name2(num)) == 'albdifvis_ice' .or. trim(sfc_name2(num)) == 'albdifnir_ice' &
                                            .or. trim(sfc_name2(num)) == 'emis_lnd'      .or. trim(sfc_name2(num)) == 'emis_ice'      &
                                            .or. trim(sfc_name2(num)) == 'sncovr_ice') then
           if(is_lsoil) then
              call register_restart_field(Sfc_restart, sfc_name2(num), var2_p, dimensions=(/'lat','lon'/), is_optional=.true.)
           else
              call register_restart_field(Sfc_restart, sfc_name2(num), var2_p, dimensions=(/'Time   ','yaxis_1','xaxis_1'/),&
                                         &is_optional=.true.)
           end if
        else
           if(is_lsoil) then
              call register_restart_field(Sfc_restart,sfc_name2(num),var2_p, dimensions=(/'lat','lon'/))
           else
              call register_restart_field(Sfc_restart,sfc_name2(num),var2_p, dimensions=(/'Time   ','yaxis_1','xaxis_1'/))
           end if
        endif
     enddo

      if (Model%nstf_name(1) > 0) then
         mand = .false.
         if (Model%nstf_name(2) == 0) mand = .true.
         do num = nvar_s2m+1,nvar_s2m+nvar_s2o
            var2_p => sfc_var2(:,:,num)
            if(is_lsoil) then
               call register_restart_field(Sfc_restart, sfc_name2(num), var2_p, dimensions=(/'lat','lon'/), is_optional=.not.mand)
            else
               call register_restart_field(Sfc_restart, sfc_name2(num), var2_p, dimensions=(/'Time   ','yaxis_1','xaxis_1'/), &
                                          &is_optional=.not.mand)
            endif
         enddo
      endif

      if (Model%lsm == Model%lsm_ruc) then ! nvar_s2mp = 0
         do num = nvar_s2m+nvar_s2o+1, nvar_s2m+nvar_s2o+nvar_s2r
            var2_p => sfc_var2(:,:,num)
            if(is_lsoil) then
               call register_restart_field(Sfc_restart, sfc_name2(num), var2_p, dimensions=(/'lat','lon'/) )
            else
               call register_restart_field(Sfc_restart, sfc_name2(num), var2_p, dimensions=(/'Time   ','yaxis_1','xaxis_1'/) )
            end if
         enddo
      endif ! mp/ruc


! Noah MP register only necessary only lsm = 2, not necessary has values
      if (nvar_s2mp > 0) then
         mand = .false.
         do num = nvar_s2m+nvar_s2o+1,nvar_s2m+nvar_s2o+nvar_s2mp
            var2_p => sfc_var2(:,:,num)
            if(is_lsoil) then
               call register_restart_field(Sfc_restart, sfc_name2(num), var2_p, dimensions=(/'lat','lon'/), is_optional=.not.mand)
            else
               call register_restart_field(Sfc_restart, sfc_name2(num), var2_p, dimensions=(/'Time   ','yaxis_1','xaxis_1'/), &
                                          &is_optional=.not.mand)
            end if
         enddo
      endif ! noahmp
      nullify(var2_p)
   endif  ! if not allocated


    if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp .or. (.not.warm_start)) then
      !--- names of the 3D variables to save
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
      !--- names of the 2D variables to save
      sfc_name3(1) = 'tslb'
      sfc_name3(2) = 'smois'
      sfc_name3(3) = 'sh2o'
      sfc_name3(4) = 'smfr'
      sfc_name3(5) = 'flfr'
    endif

      !--- register the 3D fields
    sfc_name3(0) = 'tiice'
    var3_p => sfc_var3ice(:,:,:)
    call register_restart_field(Sfc_restart, sfc_name3(0), var3_p, dimensions=(/'xaxis_1', 'yaxis_1', 'zaxis_1', 'Time   '/),&
                              &is_optional=.true.)

    do num = 1,nvar_s3
       var3_p => sfc_var3(:,:,:,num)
       if ( warm_start ) then
          call register_restart_field(Sfc_restart, sfc_name3(num), var3_p, dimensions=(/'xaxis_1', 'yaxis_1', 'lsoil  ', 'Time   '/),&
                                     &is_optional=.true.)
       else
          if(is_lsoil) then
             call register_restart_field(Sfc_restart, sfc_name3(num), var3_p, dimensions=(/'lat  ', 'lon  ', 'lsoil'/), is_optional=.true.)
          else
             call register_restart_field(Sfc_restart, sfc_name3(num), var3_p, dimensions=(/'xaxis_1','yaxis_1','zaxis_1','Time   '/),&
                                        &is_optional=.true.)
          end if
       end if
    enddo

    if (Model%lsm == Model%lsm_noahmp) then
       mand = .false.
       do num = nvar_s3+1,nvar_s3+3
          var3_p1 => sfc_var3sn(:,:,:,num)
          call register_restart_field(Sfc_restart, sfc_name3(num), var3_p1, dimensions=(/'xaxis_1', 'yaxis_1','zaxis_2', 'Time   '/),&
                                     &is_optional=.not.mand)
       enddo

       var3_p2 => sfc_var3eq(:,:,:,7)
       call register_restart_field(Sfc_restart, sfc_name3(7), var3_p2, dimensions=(/'xaxis_1', 'yaxis_1', 'zaxis_3', 'Time   '/),&
                                  &is_optional=.not.mand)

       var3_p3 => sfc_var3zn(:,:,:,8)
       call register_restart_field(Sfc_restart, sfc_name3(8), var3_p3, dimensions=(/'xaxis_1', 'yaxis_1', 'zaxis_4', 'Time   '/),&
                                  &is_optional=.not.mand)

       nullify(var3_p1)
       nullify(var3_p2)
       nullify(var3_p3)
    endif   !mp

    nullify(var3_p)

!--- Noah MP define arbitrary value (number layers of snow) to indicate
!coldstart(sfcfile doesn't include noah mp fields) or not

    if (Model%lsm == Model%lsm_noahmp) then
      sfc_var2(1,1,nvar_s2m+19) = -66666.0_r8
    endif

    !--- read the surface restart/data
    call mpp_error(NOTE,'reading surface properties data from INPUT/sfc_data.tile*.nc')
    call read_restart(Sfc_restart, ignore_checksum=ignore_rst_cksum)
    call close_file(Sfc_restart)

    if(Model%rrfs_sd) then
      call rrfs_sd_data%copy_from_temporaries(Model,Sfcprop,Atm_block)
    end if

!   write(0,*)' stype read in min,max=',minval(sfc_var2(:,:,35)),maxval(sfc_var2(:,:,35)),' sfc_name2=',sfc_name2(35)
!   write(0,*)' stype read in min,max=',minval(sfc_var2(:,:,18)),maxval(sfc_var2(:,:,18))
!   write(0,*)' sfc_var2=',sfc_var2(:,:,12)

    !--- place the data into the block GFS containers

!$omp parallel do default(shared) private(i, j, nb, ix, nt, ii1, jj1, lsoil)
    block_loop: do nb = 1, Atm_block%nblks
       allocate(ii1(Atm_block%blksz(nb)))
       allocate(jj1(Atm_block%blksz(nb)))
       ii1=Atm_block%index(nb)%ii - isc + 1
       jj1=Atm_block%index(nb)%jj - jsc + 1

       nt=0

!--- 2D variables
!    ------------
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%slmsk)   !--- slmsk
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tsfco)   !--- tsfc (tsea in sfc file)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%weasd)   !--- weasd (sheleg in sfc file)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tg3)     !--- tg3
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%zorl)    !--- zorl composite
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%alvsf)   !--- alvsf
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%alvwf)   !--- alvwf
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%alnsf)   !--- alnsf
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%alnwf)   !--- alnwf
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%facsf)   !--- facsf
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%facwf)   !--- facwf
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%vfrac)   !--- vfrac
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%canopy)  !--- canopy
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%f10m)    !--- f10m
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%t2m)     !--- t2m
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%q2m)     !--- q2m
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%vtype)   !--- vtype
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%stype)   !--- stype
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%uustar)  !--- uustar
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%ffmm)    !--- ffmm
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%ffhh)    !--- ffhh
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%hice)    !--- hice
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%fice)    !--- fice
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tisfc)   !--- tisfc
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tprcp)   !--- tprcp
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%srflag)  !--- srflag
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%snowd)   !--- snowd (snwdph in the file)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%shdmin)  !--- shdmin
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%shdmax)  !--- shdmax
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%slope)   !--- slope
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%snoalb)  !--- snoalb
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%sncovr)  !--- sncovr
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%snodl)   !--- snodl (snowd on land  portion of a cell)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%weasdl)  !--- weasdl (weasd on land  portion of a cell)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tsfc)    !--- tsfc composite
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tsfcl)   !--- tsfcl  (temp on land portion of a cell)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%zorlw)   !--- zorlw (zorl on water portion of a cell)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%zorll)   !--- zorll (zorl on land portion of a cell)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%zorli)   !--- zorli (zorl on ice  portion of a cell)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%albdirvis_lnd)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%albdirnir_lnd)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%albdifvis_lnd)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%albdifnir_lnd)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%emis_lnd)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%emis_ice)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%sncovr_ice)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%snodi)   !--- snodi (snowd on ice  portion of a cell)
        call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%weasdi)  !--- weasdi (weasd on ice  portion of a cell)
        if (Model%use_cice_alb .or. Model%lsm == Model%lsm_ruc) then
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%albdirvis_ice)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%albdifvis_ice)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%albdirnir_ice)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%albdifnir_ice)
!         call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%sfalb_ice)
        endif
        if(Model%cplwav) then
          !tgs - the following line is a bug. It should be nt = nt
          !nt = nvar_s2m-1 ! Next item will be at nvar_s2m
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%zorlwav) !--- (zorl from wave model)
        else
          Sfcprop(nb)%zorlwav  = Sfcprop(nb)%zorlw
        endif

        do_lsi_fractions: do ix = 1, Atm_block%blksz(nb)
         if (Sfcprop(nb)%stype(ix) == 14 .or. Sfcprop(nb)%stype(ix) <= 0) then
          Sfcprop(nb)%landfrac(ix) = zero
          Sfcprop(nb)%stype(ix) = 0
          if (Sfcprop(nb)%lakefrac(ix) > zero) then
            Sfcprop(nb)%lakefrac(ix) = one
          endif
         endif

         if_frac_grid: if (Model%frac_grid) then
          if (Sfcprop(nb)%landfrac(ix) > -999.0_r8) then
            Sfcprop(nb)%slmsk(ix) = ceiling(Sfcprop(nb)%landfrac(ix)-1.0e-6)
            if (Sfcprop(nb)%slmsk(ix) == 1 .and. Sfcprop(nb)%stype(ix) == 14) &
              Sfcprop(nb)%slmsk(ix) = 0
            if (Sfcprop(nb)%lakefrac(ix) > zero) then
              Sfcprop(nb)%oceanfrac(ix) = zero ! lake & ocean don't coexist in a cell
              if (nint(Sfcprop(nb)%slmsk(ix)) /= 1) then
                if(Sfcprop(nb)%fice(ix) >= Model%min_lakeice) then
                  Sfcprop(nb)%slmsk(ix) = 2
                else
                  Sfcprop(nb)%slmsk(ix) = 0
                endif
              endif
            else
              Sfcprop(nb)%lakefrac(ix)  = zero
              Sfcprop(nb)%oceanfrac(ix) = one - Sfcprop(nb)%landfrac(ix)
              if (nint(Sfcprop(nb)%slmsk(ix)) /= 1) then
                if (Sfcprop(nb)%fice(ix) >= Model%min_seaice) then
                  Sfcprop(nb)%slmsk(ix) = 2
                else
                  Sfcprop(nb)%slmsk(ix) = 0
                endif
              endif
            endif
          else
            Model%frac_grid = .false.
            if (nint(Sfcprop(nb)%slmsk(ix)) == 1) then
              Sfcprop(nb)%landfrac(ix)  = one
              Sfcprop(nb)%lakefrac(ix)  = zero
              Sfcprop(nb)%oceanfrac(ix) = zero
            else
              if (Sfcprop(nb)%slmsk(ix) < 0.1_r8 .or. Sfcprop(nb)%slmsk(ix) > 1.9_r8) then
                Sfcprop(nb)%landfrac(ix) = zero
                if (Sfcprop(nb)%oro_uf(ix) > min_lake_orog) then   ! lakes
                  Sfcprop(nb)%lakefrac(ix)  = one
                  Sfcprop(nb)%oceanfrac(ix) = zero
                else                                               ! ocean
                  Sfcprop(nb)%lakefrac(ix)  = zero
                  Sfcprop(nb)%oceanfrac(ix) = one
                endif
              endif
            endif
          endif
         else                                             ! not a fractional grid
          if (Sfcprop(nb)%landfrac(ix) > -999.0_r8) then
            if (Sfcprop(nb)%lakefrac(ix) > zero) then
              Sfcprop(nb)%oceanfrac(ix) = zero
              Sfcprop(nb)%landfrac(ix)  = zero
              Sfcprop(nb)%lakefrac(ix)  = one
              Sfcprop(nb)%slmsk(ix)     = zero
              if (Sfcprop(nb)%fice(ix) >= Model%min_lakeice) Sfcprop(nb)%slmsk(ix) = 2.0
            else
              Sfcprop(nb)%slmsk(ix) = nint(Sfcprop(nb)%landfrac(ix))
              if (Sfcprop(nb)%stype(ix) <= 0 .or. Sfcprop(nb)%stype(ix) == 14) &
                Sfcprop(nb)%slmsk(ix) = zero
              if (nint(Sfcprop(nb)%slmsk(ix)) == 0) then
                Sfcprop(nb)%oceanfrac(ix) = one
                Sfcprop(nb)%landfrac(ix)  = zero
                Sfcprop(nb)%lakefrac(ix)  = zero
                if (Sfcprop(nb)%fice(ix) >= Model%min_seaice) Sfcprop(nb)%slmsk(ix) = 2.0
              else
                Sfcprop(nb)%landfrac(ix)  = one
                Sfcprop(nb)%lakefrac(ix)  = zero
                Sfcprop(nb)%oceanfrac(ix) = zero
              endif
            endif
          else
            if (nint(Sfcprop(nb)%slmsk(ix)) == 1 .and. Sfcprop(nb)%stype(ix) > 0      &
                                                 .and. Sfcprop(nb)%stype(ix) /= 14) then
              Sfcprop(nb)%landfrac(ix)  = one
              Sfcprop(nb)%lakefrac(ix)  = zero
              Sfcprop(nb)%oceanfrac(ix) = zero
            else
              Sfcprop(nb)%slmsk(ix)    = zero
              Sfcprop(nb)%landfrac(ix) = zero
              if (Sfcprop(nb)%oro_uf(ix) > min_lake_orog) then   ! lakes
                Sfcprop(nb)%lakefrac(ix) = one
                Sfcprop(nb)%oceanfrac(ix) = zero
                if (Sfcprop(nb)%fice(ix) > Model%min_lakeice) Sfcprop(nb)%slmsk(ix) = 2.0
              else                                       ! ocean
                Sfcprop(nb)%lakefrac(ix)  = zero
                Sfcprop(nb)%oceanfrac(ix) = one
                if (Sfcprop(nb)%fice(ix) > Model%min_seaice) Sfcprop(nb)%slmsk(ix) = 2.0
              endif
            endif
          endif
         endif if_frac_grid
        enddo do_lsi_fractions

        if (warm_start .and. Model%kdt > 1) then
          do ix = 1, Atm_block%blksz(nb)
            Sfcprop(nb)%slmsk(ix)  = sfc_var2(ii1(ix),jj1(ix),1)    !--- slmsk
          enddo
        endif

        !
        !--- NSSTM variables
        !tgs - the following line is a bug that will show if(Model%cplwav) = true
        !nt = nvar_s2m 
        if (Model%nstf_name(1) > 0) then
          if (Model%nstf_name(2) == 1) then             ! nsst spinup
          !--- nsstm tref
            nt = nt + 18
            Sfcprop(nb)%tref    = Sfcprop(nb)%tsfco
            Sfcprop(nb)%z_c     = zero
            Sfcprop(nb)%c_0     = zero
            Sfcprop(nb)%c_d     = zero
            Sfcprop(nb)%w_0     = zero
            Sfcprop(nb)%w_d     = zero
            Sfcprop(nb)%xt      = zero
            Sfcprop(nb)%xs      = zero
            Sfcprop(nb)%xu      = zero
            Sfcprop(nb)%xv      = zero
            Sfcprop(nb)%xz      = 20.0_r8
            Sfcprop(nb)%zm      = zero
            Sfcprop(nb)%xtts    = zero
            Sfcprop(nb)%xzts    = zero
            Sfcprop(nb)%d_conv  = zero
            Sfcprop(nb)%ifd     = zero
            Sfcprop(nb)%dt_cool = zero
            Sfcprop(nb)%qrain   = zero
          elseif (Model%nstf_name(2) == 0) then         ! nsst restart
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tref)  !--- nsstm tref
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%z_c)  !--- nsstm z_c
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%c_0)  !--- nsstm c_0
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%c_d)  !--- nsstm c_d
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%w_0)  !--- nsstm w_0
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%w_d)  !--- nsstm w_d
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xt)  !--- nsstm xt
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xs)  !--- nsstm xs
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xu)  !--- nsstm xu
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xv) !--- nsstm xv
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xz) !--- nsstm xz
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%zm) !--- nsstm zm
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xtts) !--- nsstm xtts
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xzts) !--- nsstm xzts
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%d_conv) !--- nsstm d_conv
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%ifd) !--- nsstm ifd
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%dt_cool) !--- nsstm dt_cool
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%qrain) !--- nsstm qrain
          endif
        endif

        if (Model%lsm == Model%lsm_ruc .and. warm_start) then
          !--- Extra RUC variables
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%wetness)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%clw_surf_land)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%clw_surf_ice)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%qwv_surf_land)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%qwv_surf_ice)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tsnow_land)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tsnow_ice)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%snowfallac_land)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%snowfallac_ice)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%sfalb_lnd)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%sfalb_lnd_bck)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%sfalb_ice)
          if (Model%rdlai) then
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xlaixy)
          endif
        else if (Model%lsm == Model%lsm_ruc) then
          ! Initialize RUC snow cover on ice from snow cover
          Sfcprop(nb)%sncovr_ice = Sfcprop(nb)%sncovr
          if (Model%rdlai) then
            call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xlaixy)
          end if
        elseif (Model%lsm == Model%lsm_noahmp) then
          !--- Extra Noah MP variables
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%snowxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tvxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tgxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%canicexy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%canliqxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%eahxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%tahxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%cmxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%chxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%fwetxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%sneqvoxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%alboldxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%qsnowxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%wslakexy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%zwtxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%waxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%wtxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%lfmassxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%rtmassxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%stmassxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%woodxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%stblcpxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%fastcpxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xsaixy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%xlaixy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%taussxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%smcwtdxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%deeprechxy)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,sfc_var2,Sfcprop(nb)%rechxy)
        endif

        if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp .or. (.not.warm_start)) then
          !--- 3D variables
          nt=0
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,1,Model%lsoil,sfc_var3,Sfcprop(nb)%stc)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,1,Model%lsoil,sfc_var3,Sfcprop(nb)%smc)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,1,Model%lsoil,sfc_var3,Sfcprop(nb)%slc)

          if (Model%lsm == Model%lsm_noahmp) then
            ! These use weird indexing which is lost during a Fortran subroutine call, so we use loops instead:
            nt=nt+1
            do lsoil = -2, 0
              do ix = 1, Atm_block%blksz(nb)
                Sfcprop(nb)%snicexy(ix,lsoil) = sfc_var3sn(ii1(ix),jj1(ix),lsoil,nt)
              enddo
            enddo
            
            nt=nt+1
            do lsoil = -2, 0
              do ix = 1, Atm_block%blksz(nb)
                Sfcprop(nb)%snliqxy(ix,lsoil) = sfc_var3sn(ii1(ix),jj1(ix),lsoil,nt)
              enddo
            enddo

            nt=nt+1
            do lsoil = -2, 0
              do ix = 1, Atm_block%blksz(nb)
                Sfcprop(nb)%tsnoxy(ix,lsoil)  = sfc_var3sn(ii1(ix),jj1(ix),lsoil,nt)
              enddo
            enddo

            nt=nt+1
            do lsoil = 1, 4
              do ix = 1, Atm_block%blksz(nb)
                Sfcprop(nb)%smoiseq(ix,lsoil)  = sfc_var3eq(ii1(ix),jj1(ix),lsoil,nt)
              enddo
            enddo

            nt=nt+1
            do lsoil = -2, 4
              do ix = 1, Atm_block%blksz(nb)
                Sfcprop(nb)%zsnsoxy(ix,lsoil)  = sfc_var3zn(ii1(ix),jj1(ix),lsoil,nt)
              enddo
            enddo
          endif

        else if (Model%lsm == Model%lsm_ruc) then
          !--- 3D variables
          nt=0
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,1,Model%lsoil_lsm,sfc_var3,Sfcprop(nb)%tslb)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,1,Model%lsoil_lsm,sfc_var3,Sfcprop(nb)%smois)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,1,Model%lsoil_lsm,sfc_var3,Sfcprop(nb)%sh2o)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,1,Model%lsoil_lsm,sfc_var3,Sfcprop(nb)%keepsmfr)
          call copy_to_GFS_Data(ii1,jj1,isc,jsc,nt,1,Model%lsoil_lsm,sfc_var3,Sfcprop(nb)%flag_frsoil)
        endif

        do k = 1,Model%kice
         do ix = 1, Atm_block%blksz(nb)
          Sfcprop(nb)%tiice(ix,k) = sfc_var3ice(ii1(ix),jj1(ix),k)   !--- internal ice temp
         enddo
        enddo

        deallocate(ii1,jj1)

      end do block_loop
    call mpp_error(NOTE, 'gfs_driver:: - after put to container ')

! so far: At cold start everything is 9999.0, warm start snowxy has values
!         but the 3D of snow fields are not available because not allocated yet.
!         ix,nb loops may be consolidate with the Noah MP isnowxy init
!         restore traditional vars first,we need some of them to init snow fields
!         snow depth to actual snow layers; so we can allocate and register
!         note zsnsoxy is from -2:4 - isnowxy is from 0:-2, but we need
!         exact snow layers to pass 3D fields correctly, snow layers are
!         different fro grid to grid, we have to init point by point/grid.
!         It has to be done after the weasd is available
!         sfc_var2(1,1,32) is the first; we need this to allocate snow related fields

    i = Atm_block%index(1)%ii(1) - isc + 1
    j = Atm_block%index(1)%jj(1) - jsc + 1

    if (sfc_var2(i,j,33) < -9990.0_r8) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing snodl')
!$omp parallel do default(shared) private(nb, ix, tem)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          if (Sfcprop(nb)%landfrac(ix) > zero) then
            tem = one / (Sfcprop(nb)%fice(ix)*(one-Sfcprop(nb)%landfrac(ix))+Sfcprop(nb)%landfrac(ix))
            Sfcprop(nb)%snodl(ix)  = Sfcprop(nb)%snowd(ix) * tem
          else
            Sfcprop(nb)%snodl(ix)  = zero
          endif
        enddo
      enddo
    endif

    if (sfc_var2(i,j,34) < -9990.0_r8) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing weasdl')
!$omp parallel do default(shared) private(nb, ix, tem)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          if (Sfcprop(nb)%landfrac(ix) > zero) then
            tem = one / (Sfcprop(nb)%fice(ix)*(one-Sfcprop(nb)%landfrac(ix))+Sfcprop(nb)%landfrac(ix))
            Sfcprop(nb)%weasdl(ix) = Sfcprop(nb)%weasd(ix) * tem
          else
            Sfcprop(nb)%weasdl(ix) = zero
          endif
        enddo
      enddo
    endif

    if (sfc_var2(i,j,36) < -9990.0_r8) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing tsfcl')
!$omp parallel do default(shared) private(nb, ix)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          Sfcprop(nb)%tsfcl(ix) = Sfcprop(nb)%tsfco(ix) !--- compute tsfcl from existing variables
        enddo
      enddo
    endif

    if (sfc_var2(i,j,37) < -9990.0_r8) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing zorlw')
!$omp parallel do default(shared) private(nb, ix)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          if (Sfcprop(nb)%landfrac(ix) < one .and. Sfcprop(nb)%fice(ix) < one) then
          Sfcprop(nb)%zorlw(ix) = min(Sfcprop(nb)%zorl(ix), 0.317)
          endif
        enddo
      enddo
    endif

    if (sfc_var2(i,j,38) < -9990.0_r8) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing zorll')
!$omp parallel do default(shared) private(nb, ix)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          Sfcprop(nb)%zorll(ix) = Sfcprop(nb)%zorl(ix) !--- compute zorll from existing variables
        enddo
      enddo
    endif

    if (sfc_var2(i,j,39) < -9990.0_r8) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing zorli')
!$omp parallel do default(shared) private(nb, ix)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          if (Sfcprop(nb)%fice(ix)*(one-Sfcprop(nb)%landfrac(ix)) > zero) then
            Sfcprop(nb)%zorli(ix) = one
          endif
        enddo
      enddo
    endif

    if (sfc_var2(i,j,45) < -9990.0_r8) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing emis_ice')
!$omp parallel do default(shared) private(nb, ix)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          Sfcprop(nb)%emis_ice(ix) = 0.96
        enddo
      enddo
    endif

    if (sfc_var2(i,j,46) < -9990.0_r8 .and. Model%lsm /= Model%lsm_ruc) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing sncovr_ice')
!$omp parallel do default(shared) private(nb, ix)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
!         Sfcprop(nb)%sncovr_ice(ix) = Sfcprop(nb)%sncovr(ix)
          Sfcprop(nb)%sncovr_ice(ix) = zero
        enddo
      enddo
    endif

    if (sfc_var2(i,j,47) < -9990.0_r8) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing snodi')
!$omp parallel do default(shared) private(nb, ix, tem)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          if (Sfcprop(nb)%fice(ix) > zero) then
            tem = one / (Sfcprop(nb)%fice(ix)*(one-Sfcprop(nb)%landfrac(ix))+Sfcprop(nb)%landfrac(ix))
            Sfcprop(nb)%snodi(ix)  = min(Sfcprop(nb)%snowd(ix) * tem, 3.0)
          else
            Sfcprop(nb)%snodi(ix)  = zero
          endif
        enddo
      enddo
    endif

    if (sfc_var2(i,j,48) < -9990.0_r8) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing weasdi')
!$omp parallel do default(shared) private(nb, ix, tem)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          if (Sfcprop(nb)%fice(ix) > zero) then
            tem = one / (Sfcprop(nb)%fice(ix)*(one-Sfcprop(nb)%landfrac(ix))+Sfcprop(nb)%landfrac(ix))
            Sfcprop(nb)%weasdi(ix)  = Sfcprop(nb)%weasd(ix)*tem
          else
            Sfcprop(nb)%weasdi(ix)  = zero
          endif
        enddo
      enddo
    endif

    if (Model%use_cice_alb) then
      if (sfc_var2(i,j,49) < -9990.0_r8) then
!$omp parallel do default(shared) private(nb, ix)
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            if (Sfcprop(nb)%oceanfrac(ix) > zero .and. &
                Sfcprop(nb)%fice(ix) >= Model%min_seaice) then
              Sfcprop(nb)%albdirvis_ice(ix) = 0.6_kind_phys
              Sfcprop(nb)%albdifvis_ice(ix) = 0.6_kind_phys
              Sfcprop(nb)%albdirnir_ice(ix) = 0.6_kind_phys
              Sfcprop(nb)%albdifnir_ice(ix) = 0.6_kind_phys
            endif
          enddo
        enddo
      endif

    endif

      ! Fill in composite tsfc for coldstart runs - must happen after tsfcl is computed
    compute_tsfc_for_colstart: if (sfc_var2(i,j,35) < -9990.0_r8) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing composite tsfc')
      if(Model%frac_grid) then ! 3-way composite
!$omp parallel do default(shared) private(nb, ix, tem, tem1)
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            Sfcprop(nb)%tsfco(ix) = max(con_tice, Sfcprop(nb)%tsfco(ix)) ! this may break restart reproducibility
            tem1 = one - Sfcprop(nb)%landfrac(ix)
            tem  = tem1 * Sfcprop(nb)%fice(ix) ! tem = ice fraction wrt whole cell
            Sfcprop(nb)%tsfc(ix) = Sfcprop(nb)%tsfcl(ix) * Sfcprop(nb)%landfrac(ix) &
                                 + Sfcprop(nb)%tisfc(ix) * tem                      &
                                 + Sfcprop(nb)%tsfco(ix) * (tem1-tem)
          enddo
        enddo
      else
!$omp parallel do default(shared) private(nb, ix, tem)
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            if (Sfcprop(nb)%slmsk(ix) == 1) then
              Sfcprop(nb)%tsfc(ix) = Sfcprop(nb)%tsfcl(ix)
            else
              tem = one - Sfcprop(nb)%fice(ix)
              Sfcprop(nb)%tsfc(ix) = Sfcprop(nb)%tisfc(ix) * Sfcprop(nb)%fice(ix) &
                                   + Sfcprop(nb)%tsfco(ix) * tem
            endif
          enddo
        enddo
      endif
    endif compute_tsfc_for_colstart

    if (sfc_var2(i,j,nvar_s2m) < -9990.0_r8) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing zorlwav')
!$omp parallel do default(shared) private(nb, ix)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          Sfcprop(nb)%zorlwav(ix) = Sfcprop(nb)%zorl(ix) !--- compute zorlwav from existing variables
        enddo
      enddo
    endif

    if (nint(sfc_var3ice(1,1,1)) == -9999) then    !--- initialize internal ice temp from layer 1 and 2 soil temp
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing tiice')
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          Sfcprop(nb)%tiice(ix,1) = max(timin, min(con_tice, Sfcprop(nb)%stc(ix,1)))
          Sfcprop(nb)%tiice(ix,2) = max(timin, min(con_tice, Sfcprop(nb)%stc(ix,2)))
        enddo
      enddo
    endif

    ! A standard-compliant Fortran 2003 compiler will call rrfs_sd_final here

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
    integer :: nvar2m, nvar2o, nvar3
    integer :: nvar2r, nvar2mp, nvar3mp
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
    !--- temporary variables for storing rrfs_sd fields
    type(rrfs_sd_data_type) :: rrfs_sd_data

    nvar2m = 48
    if (Model%use_cice_alb .or. Model%lsm == Model%lsm_ruc) then
      nvar2m = nvar2m + 4
!     nvar2m = nvar2m + 5
    endif
    if (Model%cplwav) nvar2m = nvar2m + 1
    if (Model%nstf_name(1) > 0) then
      nvar2o = 18
    else
      nvar2o = 0
    endif
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

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)

    if (Model%lsm == Model%lsm_ruc) then
      if (allocated(sfc_name2)) then
        ! Re-allocate if one or more of the dimensions don't match
        if (size(sfc_name2).ne.nvar2m+nvar2o+nvar2mp+nvar2r .or. &
            size(sfc_name3).ne.nvar3+nvar3mp .or.                &
            size(sfc_var3,dim=3).ne.Model%lsoil_lsm) then
          !--- deallocate containers and free restart container
          deallocate(sfc_name2)
          deallocate(sfc_name3)
          deallocate(sfc_var2)
          deallocate(sfc_var3)
       end if
      end if
    end if

    !--- set filename
    infile=trim(indir)//'/'//trim(fn_srf)
    if( present(timestamp) ) infile=trim(indir)//'/'//trim(timestamp)//'.'//trim(fn_srf)

    !--- register axis
    amiopen=open_file(Sfc_restart, trim(infile), 'overwrite', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
    if_amiopen: if( amiopen ) then
      call register_axis(Sfc_restart, 'xaxis_1', 'X')
      call register_field(Sfc_restart, 'xaxis_1', 'double', (/'xaxis_1'/))
      call register_variable_attribute(Sfc_restart, 'xaxis_1', 'cartesian_axis', 'X', str_len=1)
      call get_global_io_domain_indices(Sfc_restart, 'xaxis_1', is, ie, indices=buffer)
      call write_data(Sfc_restart, "xaxis_1", buffer)
      deallocate(buffer)

      call register_axis(Sfc_restart, 'yaxis_1', 'Y')
      call register_field(Sfc_restart, 'yaxis_1', 'double', (/'yaxis_1'/))
      call register_variable_attribute(Sfc_restart, 'yaxis_1', 'cartesian_axis', 'Y', str_len=1)
      call get_global_io_domain_indices(Sfc_restart, 'yaxis_1', is, ie, indices=buffer)
      call write_data(Sfc_restart, "yaxis_1", buffer)
      deallocate(buffer)

      call register_axis(Sfc_restart, 'zaxis_1', dimension_length=Model%kice)
      call register_field(Sfc_restart, 'zaxis_1', 'double', (/'zaxis_1'/))
      call register_variable_attribute(Sfc_restart, 'zaxis_1', 'cartesian_axis', 'Z', str_len=1)
      allocate( buffer(Model%kice) )
      do i=1, Model%kice
         buffer(i) = i
      end do
      call write_data(Sfc_restart, 'zaxis_1', buffer)
      deallocate(buffer)

      if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp) then
        call register_axis(Sfc_restart, 'zaxis_2', dimension_length=Model%lsoil)
        call register_field(Sfc_restart, 'zaxis_2', 'double', (/'zaxis_2'/))
        call register_variable_attribute(Sfc_restart, 'zaxis_2', 'cartesian_axis', 'Z', str_len=1)
        allocate( buffer(Model%lsoil) )
        do i=1, Model%lsoil
          buffer(i)=i
        end do
        call write_data(Sfc_restart, 'zaxis_2', buffer)
        deallocate(buffer)
      endif

      if(Model%lsm == Model%lsm_noahmp) then
        call register_axis(Sfc_restart, 'zaxis_3', dimension_length=3)
        call register_field(Sfc_restart, 'zaxis_3', 'double', (/'zaxis_3'/))
        call register_variable_attribute(Sfc_restart, 'zaxis_3', 'cartesian_axis', 'Z', str_len=1)
        allocate(buffer(3))
        do i=1, 3
           buffer(i) = i
        end do
        call write_data(Sfc_restart, 'zaxis_3', buffer)
        deallocate(buffer)

        call register_axis(Sfc_restart, 'zaxis_4', dimension_length=7)
        call register_field(Sfc_restart, 'zaxis_4', 'double', (/'zaxis_4'/))
        call register_variable_attribute(Sfc_restart, 'zaxis_4', 'cartesian_axis' ,'Z', str_len=1)
        allocate(buffer(7))
        do i=1, 7
           buffer(i)=i
        end do
        call write_data(Sfc_restart, 'zaxis_4', buffer)
        deallocate(buffer)
      end if
      call register_axis(Sfc_restart, 'Time', unlimited)
      call register_field(Sfc_restart, 'Time', 'double', (/'Time'/))
      call register_variable_attribute(Sfc_restart, 'Time', 'cartesian_axis', 'T', str_len=1)
      call write_data( Sfc_restart, 'Time', 1)
    else
      call mpp_error(FATAL, 'Error in opening file'//trim(infile) )
    end if if_amiopen

    if(Model%rrfs_sd) then
      call rrfs_sd_data%allocate_data(Model)
      call rrfs_sd_data%register_axis(Model)
      call rrfs_sd_data%write_axis(Model)
    end if

    if (.not. allocated(sfc_name2)) then
      !--- allocate the various containers needed for restarts
      allocate(sfc_name2(nvar2m+nvar2o+nvar2mp+nvar2r))
      allocate(sfc_name3(0:nvar3+nvar3mp))
      allocate(sfc_var2(nx,ny,nvar2m+nvar2o+nvar2mp+nvar2r))
      if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp) then
        allocate(sfc_var3(nx,ny,Model%lsoil,nvar3))
      elseif (Model%lsm == Model%lsm_ruc) then
        allocate(sfc_var3(nx,ny,Model%lsoil_lsm,nvar3))
      endif
      sfc_var2   = -9999.0_r8
      sfc_var3   = -9999.0_r8
      if (Model%lsm == Model%lsm_noahmp) then
        allocate(sfc_var3sn(nx,ny,-2:0,4:6))
        allocate(sfc_var3eq(nx,ny,1:4,7:7))
        allocate(sfc_var3zn(nx,ny,-2:4,8:8))

        sfc_var3sn = -9999.0_r8
        sfc_var3eq = -9999.0_r8
        sfc_var3zn = -9999.0_r8
      endif
      call fill_Sfcprop_names(Model,sfc_name2,sfc_name3,nvar2m,.true.)
   end if

   if(Model%rrfs_sd) then
     call rrfs_sd_data%register_fields
   endif

   !--- register the 2D fields
   do num = 1,nvar2m
      var2_p => sfc_var2(:,:,num)
      if (trim(sfc_name2(num)) == 'sncovr' .or. trim(sfc_name2(num)) == 'tsfcl' .or.trim(sfc_name2(num))  == 'zorll'   &
           .or. trim(sfc_name2(num)) == 'zorli' .or.trim(sfc_name2(num))  == 'zorlwav' &
           .or. trim(sfc_name2(num)) == 'snodl' .or. trim(sfc_name2(num)) == 'weasdl'  &
           .or. trim(sfc_name2(num)) == 'snodi' .or. trim(sfc_name2(num)) == 'weasdi'  &
           .or. trim(sfc_name2(num)) == 'tsfc'  .or. trim(sfc_name2(num)) ==  'zorlw'  &
           .or. trim(sfc_name2(num)) == 'albdirvis_lnd' .or. trim(sfc_name2(num)) == 'albdirnir_lnd' &
           .or. trim(sfc_name2(num)) == 'albdifvis_lnd' .or. trim(sfc_name2(num)) == 'albdifnir_lnd' &
           .or. trim(sfc_name2(num)) == 'albdirvis_ice' .or. trim(sfc_name2(num)) == 'albdirnir_ice' &
           .or. trim(sfc_name2(num)) == 'albdifvis_ice' .or. trim(sfc_name2(num)) == 'albdifnir_ice' &
           .or. trim(sfc_name2(num)) == 'emis_lnd'      .or. trim(sfc_name2(num)) == 'emis_ice'      &
           .or. trim(sfc_name2(num)) == 'sncovr_ice' ) then
         call register_restart_field(Sfc_restart, sfc_name2(num), var2_p, dimensions=(/'xaxis_1','yaxis_1','Time   '/), is_optional=.true.)
      else
         call register_restart_field(Sfc_restart, sfc_name2(num), var2_p, dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/) )
      endif
   enddo
   if (Model%nstf_name(1) > 0) then
      mand = .false.
      if (Model%nstf_name(2) ==0) mand = .true.
      do num = nvar2m+1,nvar2m+nvar2o
         var2_p => sfc_var2(:,:,num)
         call register_restart_field(Sfc_restart, sfc_name2(num), var2_p, dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/),&
                                    &is_optional=.not.mand)
      enddo
   endif

   if (Model%lsm == Model%lsm_ruc) then ! nvar2mp =0
      do num = nvar2m+nvar2o+1, nvar2m+nvar2o+nvar2r
         var2_p => sfc_var2(:,:,num)
         call register_restart_field(Sfc_restart, sfc_name2(num), var2_p, dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/))
      enddo
   else if (Model%lsm == Model%lsm_noahmp) then ! nvar2r =0
      mand = .true.                  ! actually should be true since it is after cold start
      do num = nvar2m+nvar2o+1,nvar2m+nvar2o+nvar2mp
         var2_p => sfc_var2(:,:,num)
         call register_restart_field(Sfc_restart, sfc_name2(num), var2_p, dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/),&
                                    &is_optional=.not.mand)
      enddo
   endif
   nullify(var2_p)

   if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp) then
      !--- names of the 3D variables to save
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
      !--- names of the 3D variables to save
      sfc_name3(1) = 'tslb'
      sfc_name3(2) = 'smois'
      sfc_name3(3) = 'sh2o'
      sfc_name3(4) = 'smfr'
      sfc_name3(5) = 'flfr'
   end if

   !--- register the 3D fields
   !     if (Model%frac_grid) then
   sfc_name3(0) = 'tiice'
   var3_p => sfc_var3ice(:,:,:)
   call register_restart_field(Sfc_restart, sfc_name3(0), var3_p, dimensions=(/'xaxis_1', 'yaxis_1', 'zaxis_1', 'Time   '/))
   !     endif

   if(Model%lsm == Model%lsm_ruc) then
      do num = 1,nvar3
         var3_p => sfc_var3(:,:,:,num)
         call register_restart_field(Sfc_restart, sfc_name3(num), var3_p, dimensions=(/'xaxis_1', 'yaxis_1', 'zaxis_1', 'Time   '/))
      enddo
      nullify(var3_p)
   else
      do num = 1,nvar3
         var3_p => sfc_var3(:,:,:,num)
         call register_restart_field(Sfc_restart, sfc_name3(num), var3_p, dimensions=(/'xaxis_1', 'yaxis_1', 'zaxis_2', 'Time   '/))
      enddo
      nullify(var3_p)
   endif

   if (Model%lsm == Model%lsm_noahmp) then
      mand = .true.
      do num = nvar3+1,nvar3+3
         var3_p1 => sfc_var3sn(:,:,:,num)
         call register_restart_field(Sfc_restart, sfc_name3(num), var3_p1, dimensions=(/'xaxis_1', 'yaxis_1', 'zaxis_3', 'Time   '/),&
                                    &is_optional=.not.mand)
      enddo

      var3_p2 => sfc_var3eq(:,:,:,7)
      call register_restart_field(Sfc_restart, sfc_name3(7), var3_p2, dimensions=(/'xaxis_1', 'yaxis_1', 'zaxis_2', 'Time   '/),&
                                    &is_optional=.not.mand)

      var3_p3 => sfc_var3zn(:,:,:,8)
      call register_restart_field(Sfc_restart, sfc_name3(8), var3_p3, dimensions=(/'xaxis_1', 'yaxis_1', 'zaxis_4', 'Time   '/),&
                                 &is_optional=.not.mand)

      nullify(var3_p1)
      nullify(var3_p2)
      nullify(var3_p3)
   endif ! lsm = lsm_noahmp

   if(Model%rrfs_sd) then
     call rrfs_sd_data%copy_to_temporaries(Model,Sfcprop,Atm_block)
    endif

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

    call write_restart(Sfc_restart)
    call close_file(Sfc_restart)

    ! A standard-compliant Fortran 2003 compiler will call rrfs_sd_final here

  end subroutine sfc_prop_restart_write

  subroutine rrfs_sd_register_axis(data,Model)
    implicit none
    class(rrfs_sd_data_type) :: data
    type(GFS_control_type),      intent(in) :: Model
    call register_axis(Sfc_restart, 'fire_aux_data_levels', &
         dimension_length=Model%fire_aux_data_levels)
  end subroutine rrfs_sd_register_axis

  subroutine rrfs_sd_write_axis(data,Model)
    implicit none
    class(rrfs_sd_data_type) :: data
    type(GFS_control_type),      intent(in) :: Model
    real(kind_phys) :: fire_aux_data_levels(Model%fire_aux_data_levels)
    integer :: i

    call register_field(Sfc_restart, 'fire_aux_data_levels', 'double', (/'fire_aux_data_levels'/))
    call register_variable_attribute(Sfc_restart, 'fire_aux_data_levels', 'cartesian_axis' ,'Z', str_len=1)

    do i=1,Model%fire_aux_data_levels
      fire_aux_data_levels(i) = i
    enddo

    call write_data(Sfc_restart, 'fire_aux_data_levels', fire_aux_data_levels)
  end subroutine rrfs_sd_write_axis

  subroutine rrfs_sd_allocate_data(data,Model)
    implicit none
    class(rrfs_sd_data_type) :: data
    type(GFS_control_type),   intent(in) :: Model
    integer :: nx, ny

    call data%deallocate_data

    nx=Model%nx
    ny=Model%ny

    allocate(data%emdust(nx,ny))
    allocate(data%emseas(nx,ny))
    allocate(data%emanoc(nx,ny))
    allocate(data%fhist(nx,ny))
    allocate(data%coef_bb_dc(nx,ny))

    allocate(data%fire_in(nx,ny,Model%fire_aux_data_levels))

  end subroutine rrfs_sd_allocate_data

  subroutine rrfs_sd_fill_data(data, Model, Sfcprop, Atm_block)
    ! Fills all temporary variables with default values.
    ! Terrible things will happen if you don't call data%allocate_data first.
    ! IMPORTANT: This must match the corresponding code in sfcprop_create in
    ! GFS_typedefs.F90
    implicit none
    class(rrfs_sd_data_type) :: data
    type(GFS_sfcprop_type),   intent(in) :: Sfcprop(:)
    type(GFS_control_type),   intent(in) :: Model
    type(block_control_type), intent(in) :: Atm_block

    integer :: nb, ix, isc, jsc, i, j

    isc = Model%isc
    jsc = Model%jsc

!$omp parallel do default(shared) private(i, j, nb, ix)
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1

        data%emdust(i,j) = 0
        data%emseas(i,j) = 0
        data%emanoc(i,j) = 0
        data%fhist(i,j) = 1.
        data%coef_bb_dc(i,j) = 0

        data%fire_in(i,j,:) = 0
      end do
    end do
  end subroutine rrfs_sd_fill_data

  subroutine rrfs_sd_register_fields(data)
    ! Registers all restart fields needed by the RRFS-SD
    ! Terrible things will happen if you don't call data%allocate_data
    ! and data%register_axes first.
    implicit none
    class(rrfs_sd_data_type) :: data

    ! Register 2D fields
    call register_restart_field(Sfc_restart, 'emdust', data%emdust, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)
    call register_restart_field(Sfc_restart, 'emseas', data%emseas, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)
    call register_restart_field(Sfc_restart, 'emanoc', data%emanoc, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)
    call register_restart_field(Sfc_restart, 'fhist', data%fhist, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)
    call register_restart_field(Sfc_restart, 'coef_bb_dc', data%coef_bb_dc, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), is_optional=.true.)

    ! Register 3D field
    call register_restart_field(Sfc_restart, 'fire_in', data%fire_in, &
         dimensions=(/'xaxis_1             ', 'yaxis_1             ', &
                      'fire_aux_data_levels', 'Time                '/), &
         is_optional=.true.)
  end subroutine rrfs_sd_register_fields

  subroutine rrfs_sd_final(data)
    ! Final routine for rrfs_sd_data_type, called automatically when
    ! an object of that type goes out of scope.  This is a wrapper
    ! around data%deallocate_data() with necessary syntactic
    ! differences.
    implicit none
    type(rrfs_sd_data_type) :: data
    call rrfs_sd_deallocate_data(data)
  end subroutine rrfs_sd_final

  subroutine rrfs_sd_deallocate_data(data)
    ! Deallocates all data used, and nullifies the pointers. The data
    ! object can safely be used again after this call. This is also
    ! the implementation of the rrfs_sd_deallocate_data final routine.
    implicit none
    class(rrfs_sd_data_type) :: data

    ! This #define reduces code length by a lot
#define IF_ASSOC_DEALLOC_NULL(var) \
    if(associated(data%var)) then ; \
      deallocate(data%var) ; \
      nullify(data%var) ; \
    endif

    IF_ASSOC_DEALLOC_NULL(emdust)
    IF_ASSOC_DEALLOC_NULL(emseas)
    IF_ASSOC_DEALLOC_NULL(emanoc)
    IF_ASSOC_DEALLOC_NULL(fhist)
    IF_ASSOC_DEALLOC_NULL(coef_bb_dc)

    IF_ASSOC_DEALLOC_NULL(fire_in)

    ! Undefine this to avoid cluttering the cpp scope:
#undef IF_ASSOC_DEALLOC_NULL
  end subroutine rrfs_sd_deallocate_data

  subroutine rrfs_sd_copy_from_temporaries(data, Model, Sfcprop, Atm_block)
    implicit none
    class(rrfs_sd_data_type) :: data
    type(GFS_sfcprop_type),   intent(in) :: Sfcprop(:)
    type(GFS_control_type),   intent(in) :: Model
    type(block_control_type), intent(in) :: Atm_block

    integer :: nb, ix, isc, jsc, i, j

    isc = Model%isc
    jsc = Model%jsc

!$omp parallel do default(shared) private(i, j, nb, ix)
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1

        Sfcprop(nb)%emdust(ix) = data%emdust(i,j)
        Sfcprop(nb)%emseas(ix) = data%emseas(i,j)
        Sfcprop(nb)%emanoc(ix) = data%emanoc(i,j)
        Sfcprop(nb)%fhist(ix) = data%fhist(i,j)
        Sfcprop(nb)%coef_bb_dc(ix) = data%coef_bb_dc(i,j)

        Sfcprop(nb)%fire_in(ix,:) = data%fire_in(i,j,:)
      enddo
    enddo
  end subroutine rrfs_sd_copy_from_temporaries

  subroutine rrfs_sd_copy_to_temporaries(data, Model, Sfcprop, Atm_block)
    implicit none
    class(rrfs_sd_data_type) :: data
    type(GFS_sfcprop_type),   intent(in) :: Sfcprop(:)
    type(GFS_control_type),   intent(in) :: Model
    type(block_control_type), intent(in) :: Atm_block

    integer :: nb, ix, isc, jsc, i, j

    isc = Model%isc
    jsc = Model%jsc

!$omp parallel do default(shared) private(i, j, nb, ix)
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1

        data%emdust(i,j) = Sfcprop(nb)%emdust(ix)
        data%emseas(i,j) = Sfcprop(nb)%emseas(ix)
        data%emanoc(i,j) = Sfcprop(nb)%emanoc(ix)
        data%fhist(i,j) = Sfcprop(nb)%fhist(ix)
        data%coef_bb_dc(i,j) = Sfcprop(nb)%coef_bb_dc(ix)

        data%fire_in(i,j,:) = Sfcprop(nb)%fire_in(ix,:)
      enddo
    enddo
  end subroutine rrfs_sd_copy_to_temporaries

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
