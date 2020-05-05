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
  use fms_mod,            only: file_exist, stdout
  use fms_io_mod,         only: restart_file_type, free_restart_type, &
                                register_restart_field,               &
                                restore_state, save_restart
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
!--- GFS physics modules
!#ifndef CCPP
!--- variables needed for calculating 'sncovr'
  use namelist_soilveg,   only: salp_data, snupx
!#endif

!
! --- variables needed for Noah MP init
!
  use noahmp_tables,      only: laim_table,saim_table,sla_table,      &
                                bexp_table,smcmax_table,smcwlt_table, &
                                dwsat_table,dksat_table,psisat_table, &
                                isurban_table,isbarren_table,         &
                                isice_table,iswater_table

!
!--- GFS_typedefs
!rab  use GFS_typedefs,       only: GFS_sfcprop_type, GFS_diag_type, &
!rab                                GFS_cldprop_type, GFS_grid_type
  use GFS_typedefs,       only: GFS_sfcprop_type
!
!--- IPD typdefs
  use IPD_typedefs,       only: IPD_control_type, IPD_data_type, &
                                IPD_restart_type, IPD_diag_type, &
                                kind_phys => IPD_kind_phys
!
!-----------------------------------------------------------------------
  implicit none
  private
 
  !--- public interfaces ---
  public  FV3GFS_restart_read, FV3GFS_restart_write
  public  FV3GFS_IPD_checksum
  public  fv3gfs_diag_register, fv3gfs_diag_output
#ifdef use_WRTCOMP
  public  fv_phys_bundle_setup
#endif

  !--- GFDL filenames
  character(len=32)  :: fn_oro = 'oro_data.nc'
  character(len=32)  :: fn_srf = 'sfc_data.nc'
  character(len=32)  :: fn_phy = 'phy_data.nc'

  !--- GFDL FMS netcdf restart data types
  type(restart_file_type) :: Oro_restart, Sfc_restart, Phy_restart
 
  !--- GFDL FMS restart containers
  character(len=32),    allocatable,         dimension(:)       :: oro_name2, sfc_name2, sfc_name3
  real(kind=kind_phys), allocatable, target, dimension(:,:,:)   :: oro_var2, sfc_var2, phy_var2
  real(kind=kind_phys), allocatable, target, dimension(:,:,:,:) :: sfc_var3, phy_var3
  !--- Noah MP restart containers
  real(kind=kind_phys), allocatable, target, dimension(:,:,:,:) :: sfc_var3sn,sfc_var3eq,sfc_var3zn

  real(kind=kind_phys) :: zhour
!
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
  real, parameter :: missing_value = 9.99e20
  real, parameter:: stndrd_atmos_ps = 101325.
  real, parameter:: stndrd_atmos_lapse = 0.0065
  real, parameter:: drythresh = 1.e-4
 
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
#ifdef CCPP
  subroutine FV3GFS_restart_read (IPD_Data, IPD_Restart, Atm_block, Model, fv_domain, warm_start)
#else
  subroutine FV3GFS_restart_read (IPD_Data, IPD_Restart, Atm_block, Model, fv_domain)
#endif
    type(IPD_data_type),      intent(inout) :: IPD_Data(:)
    type(IPD_restart_type),   intent(inout) :: IPD_Restart
    type(block_control_type), intent(in)    :: Atm_block
    type(IPD_control_type),   intent(inout) :: Model
    type(domain2d),           intent(in)    :: fv_domain
#ifdef CCPP
    logical,                  intent(in)    :: warm_start
#endif
 
    !--- read in surface data from chgres 
#ifdef CCPP
    call sfc_prop_restart_read (IPD_Data%Sfcprop, Atm_block, Model, fv_domain, warm_start)
#else
    call sfc_prop_restart_read (IPD_Data%Sfcprop, Atm_block, Model, fv_domain)
#endif

    !--- read in physics restart data
    call phys_restart_read (IPD_Restart, Atm_block, Model, fv_domain)

  end subroutine FV3GFS_restart_read

!---------------------
! FV3GFS_restart_write
!---------------------
  subroutine FV3GFS_restart_write (IPD_Data, IPD_Restart, Atm_block, Model, fv_domain, timestamp)
    type(IPD_data_type),         intent(inout) :: IPD_Data(:)
    type(IPD_restart_type),      intent(inout) :: IPD_Restart
    type(block_control_type),    intent(in)    :: Atm_block
    type(IPD_control_type),      intent(in)    :: Model
    type(domain2d),              intent(in)    :: fv_domain
    character(len=32), optional, intent(in)    :: timestamp
 
    !--- write surface data from chgres 
    call sfc_prop_restart_write (IPD_Data%Sfcprop, Atm_block, Model, fv_domain, timestamp)
 
    !--- write physics restart data
    call phys_restart_write (IPD_Restart, Atm_block, Model, fv_domain, timestamp)

  end subroutine FV3GFS_restart_write


!--------------------
! FV3GFS_IPD_checksum
!--------------------
 subroutine FV3GFS_IPD_checksum (Model, IPD_Data, Atm_block)
   !--- interface variables
   type(IPD_control_type),    intent(in) :: Model
   type(IPD_data_type),       intent(in) :: IPD_Data(:)
   type (block_control_type), intent(in) :: Atm_block
   !--- local variables
   integer :: outunit, j, i, ix, nb, isc, iec, jsc, jec, lev, ct, l, ntr
   integer :: nsfcprop2d, idx_opt
   real(kind=kind_phys), allocatable :: temp2d(:,:,:)
   real(kind=kind_phys), allocatable :: temp3d(:,:,:,:)
   real(kind=kind_phys), allocatable :: temp3dlevsp1(:,:,:,:)
   character(len=32) :: name

   isc = Model%isc
   iec = Model%isc+Model%nx-1
   jsc = Model%jsc
   jec = Model%jsc+Model%ny-1
   lev = Model%levs

   ntr = size(IPD_Data(1)%Statein%qgrs,3)

   if(Model%lsm == Model%lsm_noahmp) then
     nsfcprop2d = 149  
   else
     nsfcprop2d = 100
   endif

   allocate (temp2d(isc:iec,jsc:jec,nsfcprop2d+Model%ntot3d+Model%nctp))
   allocate (temp3d(isc:iec,jsc:jec,1:lev,14+Model%ntot3d+2*ntr))
   allocate (temp3dlevsp1(isc:iec,jsc:jec,1:lev+1,3))

   temp2d = 0.
   temp3d = 0.
   temp3dlevsp1 = 0.

   do j=jsc,jec
     do i=isc,iec
       nb = Atm_block%blkno(i,j) 
       ix = Atm_block%ixp(i,j) 
       !--- statein pressure
       temp2d(i,j, 1) = IPD_Data(nb)%Statein%pgr(ix)
       temp2d(i,j, 2) = IPD_Data(nb)%Sfcprop%slmsk(ix)
       temp2d(i,j, 3) = IPD_Data(nb)%Sfcprop%tsfc(ix)
       temp2d(i,j, 4) = IPD_Data(nb)%Sfcprop%tisfc(ix)
       temp2d(i,j, 5) = IPD_Data(nb)%Sfcprop%snowd(ix)
       temp2d(i,j, 6) = IPD_Data(nb)%Sfcprop%zorl(ix)
       temp2d(i,j, 7) = IPD_Data(nb)%Sfcprop%fice(ix)
       temp2d(i,j, 8) = IPD_Data(nb)%Sfcprop%hprime(ix,1)
       temp2d(i,j, 9) = IPD_Data(nb)%Sfcprop%sncovr(ix)
       temp2d(i,j,10) = IPD_Data(nb)%Sfcprop%snoalb(ix)
       temp2d(i,j,11) = IPD_Data(nb)%Sfcprop%alvsf(ix)
       temp2d(i,j,12) = IPD_Data(nb)%Sfcprop%alnsf(ix)
       temp2d(i,j,13) = IPD_Data(nb)%Sfcprop%alvwf(ix)
       temp2d(i,j,14) = IPD_Data(nb)%Sfcprop%alnwf(ix)
       temp2d(i,j,15) = IPD_Data(nb)%Sfcprop%facsf(ix)
       temp2d(i,j,16) = IPD_Data(nb)%Sfcprop%facwf(ix)
       temp2d(i,j,17) = IPD_Data(nb)%Sfcprop%slope(ix)
       temp2d(i,j,18) = IPD_Data(nb)%Sfcprop%shdmin(ix)
       temp2d(i,j,19) = IPD_Data(nb)%Sfcprop%shdmax(ix)
       temp2d(i,j,20) = IPD_Data(nb)%Sfcprop%tg3(ix)
       temp2d(i,j,21) = IPD_Data(nb)%Sfcprop%vfrac(ix)
       temp2d(i,j,22) = IPD_Data(nb)%Sfcprop%vtype(ix)
       temp2d(i,j,23) = IPD_Data(nb)%Sfcprop%stype(ix)
       temp2d(i,j,24) = IPD_Data(nb)%Sfcprop%uustar(ix)
       temp2d(i,j,25) = IPD_Data(nb)%Sfcprop%oro(ix)
       temp2d(i,j,26) = IPD_Data(nb)%Sfcprop%oro_uf(ix)
       temp2d(i,j,27) = IPD_Data(nb)%Sfcprop%hice(ix)
       temp2d(i,j,28) = IPD_Data(nb)%Sfcprop%weasd(ix)
       temp2d(i,j,29) = IPD_Data(nb)%Sfcprop%canopy(ix)
       temp2d(i,j,30) = IPD_Data(nb)%Sfcprop%ffmm(ix)
       temp2d(i,j,31) = IPD_Data(nb)%Sfcprop%ffhh(ix)
       temp2d(i,j,32) = IPD_Data(nb)%Sfcprop%f10m(ix)
       temp2d(i,j,33) = IPD_Data(nb)%Sfcprop%tprcp(ix)
       temp2d(i,j,34) = IPD_Data(nb)%Sfcprop%srflag(ix)
#ifdef CCPP
     if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp) then
#endif
       temp2d(i,j,35) = IPD_Data(nb)%Sfcprop%slc(ix,1)
       temp2d(i,j,36) = IPD_Data(nb)%Sfcprop%slc(ix,2)
       temp2d(i,j,37) = IPD_Data(nb)%Sfcprop%slc(ix,3)
       temp2d(i,j,38) = IPD_Data(nb)%Sfcprop%slc(ix,4)
       temp2d(i,j,39) = IPD_Data(nb)%Sfcprop%smc(ix,1)
       temp2d(i,j,40) = IPD_Data(nb)%Sfcprop%smc(ix,2)
       temp2d(i,j,41) = IPD_Data(nb)%Sfcprop%smc(ix,3)
       temp2d(i,j,42) = IPD_Data(nb)%Sfcprop%smc(ix,4)
       temp2d(i,j,43) = IPD_Data(nb)%Sfcprop%stc(ix,1)
       temp2d(i,j,44) = IPD_Data(nb)%Sfcprop%stc(ix,2)
       temp2d(i,j,45) = IPD_Data(nb)%Sfcprop%stc(ix,3)
       temp2d(i,j,46) = IPD_Data(nb)%Sfcprop%stc(ix,4)
#ifdef CCPP
     elseif (Model%lsm == Model%lsm_ruc) then
       temp2d(i,j,35) = IPD_Data(nb)%Sfcprop%sh2o(ix,1)
       temp2d(i,j,36) = IPD_Data(nb)%Sfcprop%sh2o(ix,2)
       temp2d(i,j,37) = IPD_Data(nb)%Sfcprop%sh2o(ix,3)
       ! Combine levels 4 to lsoil_lsm (9 for RUC) into one
       temp2d(i,j,38) = sum(IPD_Data(nb)%Sfcprop%sh2o(ix,4:Model%lsoil_lsm))
       temp2d(i,j,39) = IPD_Data(nb)%Sfcprop%smois(ix,1)
       temp2d(i,j,40) = IPD_Data(nb)%Sfcprop%smois(ix,2)
       temp2d(i,j,41) = IPD_Data(nb)%Sfcprop%smois(ix,3)
       ! Combine levels 4 to lsoil_lsm (9 for RUC) into one
       temp2d(i,j,42) = sum(IPD_Data(nb)%Sfcprop%smois(ix,4:Model%lsoil_lsm))
       temp2d(i,j,43) = IPD_Data(nb)%Sfcprop%tslb(ix,1)
       temp2d(i,j,44) = IPD_Data(nb)%Sfcprop%tslb(ix,2)
       temp2d(i,j,45) = IPD_Data(nb)%Sfcprop%tslb(ix,3)
       ! Combine levels 4 to lsoil_lsm (9 for RUC) into one
       temp2d(i,j,46) = sum(IPD_Data(nb)%Sfcprop%tslb(ix,4:Model%lsoil_lsm))
     endif ! LSM choice
#endif
       temp2d(i,j,47) = IPD_Data(nb)%Sfcprop%t2m(ix)
       temp2d(i,j,48) = IPD_Data(nb)%Sfcprop%q2m(ix)
       temp2d(i,j,49) = IPD_Data(nb)%Coupling%nirbmdi(ix)
       temp2d(i,j,50) = IPD_Data(nb)%Coupling%nirdfdi(ix)
       temp2d(i,j,51) = IPD_Data(nb)%Coupling%visbmdi(ix)
       temp2d(i,j,52) = IPD_Data(nb)%Coupling%visdfdi(ix)
       temp2d(i,j,53) = IPD_Data(nb)%Coupling%nirbmui(ix)
       temp2d(i,j,54) = IPD_Data(nb)%Coupling%nirdfui(ix)
       temp2d(i,j,55) = IPD_Data(nb)%Coupling%visbmui(ix)
       temp2d(i,j,56) = IPD_Data(nb)%Coupling%visdfui(ix)
       temp2d(i,j,57) = IPD_Data(nb)%Coupling%sfcdsw(ix)
       temp2d(i,j,58) = IPD_Data(nb)%Coupling%sfcnsw(ix)
       temp2d(i,j,59) = IPD_Data(nb)%Coupling%sfcdlw(ix)
       temp2d(i,j,60) = IPD_Data(nb)%Grid%xlon(ix)
       temp2d(i,j,61) = IPD_Data(nb)%Grid%xlat(ix)
       temp2d(i,j,62) = IPD_Data(nb)%Grid%xlat_d(ix)
       temp2d(i,j,63) = IPD_Data(nb)%Grid%sinlat(ix)
       temp2d(i,j,64) = IPD_Data(nb)%Grid%coslat(ix)
       temp2d(i,j,65) = IPD_Data(nb)%Grid%area(ix)
       temp2d(i,j,66) = IPD_Data(nb)%Grid%dx(ix)
       if (Model%ntoz > 0) then
         temp2d(i,j,67) = IPD_Data(nb)%Grid%ddy_o3(ix)
       endif
       if (Model%h2o_phys) then
         temp2d(i,j,68) = IPD_Data(nb)%Grid%ddy_h(ix)
       endif
       temp2d(i,j,69) = IPD_Data(nb)%Cldprop%cv(ix)
       temp2d(i,j,70) = IPD_Data(nb)%Cldprop%cvt(ix)
       temp2d(i,j,71) = IPD_Data(nb)%Cldprop%cvb(ix)
       temp2d(i,j,72) = IPD_Data(nb)%Radtend%sfalb(ix)
       temp2d(i,j,73) = IPD_Data(nb)%Radtend%coszen(ix)
       temp2d(i,j,74) = IPD_Data(nb)%Radtend%tsflw(ix)
       temp2d(i,j,75) = IPD_Data(nb)%Radtend%semis(ix)
       temp2d(i,j,76) = IPD_Data(nb)%Radtend%coszdg(ix)
       temp2d(i,j,77) = IPD_Data(nb)%Radtend%sfcfsw(ix)%upfxc
       temp2d(i,j,78) = IPD_Data(nb)%Radtend%sfcfsw(ix)%upfx0
       temp2d(i,j,79) = IPD_Data(nb)%Radtend%sfcfsw(ix)%dnfxc
       temp2d(i,j,80) = IPD_Data(nb)%Radtend%sfcfsw(ix)%dnfx0
       temp2d(i,j,81) = IPD_Data(nb)%Radtend%sfcflw(ix)%upfxc
       temp2d(i,j,82) = IPD_Data(nb)%Radtend%sfcflw(ix)%upfx0
       temp2d(i,j,83) = IPD_Data(nb)%Radtend%sfcflw(ix)%dnfxc
       temp2d(i,j,84) = IPD_Data(nb)%Radtend%sfcflw(ix)%dnfx0

        idx_opt = 85 
       if (Model%lsm == Model%lsm_noahmp) then
        temp2d(i,j,idx_opt) = IPD_Data(nb)%Sfcprop%snowxy(ix)
        temp2d(i,j,idx_opt+1) = IPD_Data(nb)%Sfcprop%tvxy(ix)
        temp2d(i,j,idx_opt+2) = IPD_Data(nb)%Sfcprop%tgxy(ix)
        temp2d(i,j,idx_opt+3) = IPD_Data(nb)%Sfcprop%canicexy(ix)
        temp2d(i,j,idx_opt+4) = IPD_Data(nb)%Sfcprop%canliqxy(ix)
        temp2d(i,j,idx_opt+5) = IPD_Data(nb)%Sfcprop%eahxy(ix)
        temp2d(i,j,idx_opt+6) = IPD_Data(nb)%Sfcprop%tahxy(ix)
        temp2d(i,j,idx_opt+7) = IPD_Data(nb)%Sfcprop%cmxy(ix)
        temp2d(i,j,idx_opt+8) = IPD_Data(nb)%Sfcprop%chxy(ix)
        temp2d(i,j,idx_opt+9) = IPD_Data(nb)%Sfcprop%fwetxy(ix)
        temp2d(i,j,idx_opt+10) = IPD_Data(nb)%Sfcprop%sneqvoxy(ix)
        temp2d(i,j,idx_opt+11) = IPD_Data(nb)%Sfcprop%alboldxy(ix)
        temp2d(i,j,idx_opt+12) = IPD_Data(nb)%Sfcprop%qsnowxy(ix)
        temp2d(i,j,idx_opt+13) = IPD_Data(nb)%Sfcprop%wslakexy(ix)
        temp2d(i,j,idx_opt+14) = IPD_Data(nb)%Sfcprop%zwtxy(ix)
        temp2d(i,j,idx_opt+15) = IPD_Data(nb)%Sfcprop%waxy(ix)
        temp2d(i,j,idx_opt+16) = IPD_Data(nb)%Sfcprop%wtxy(ix)
        temp2d(i,j,idx_opt+17) = IPD_Data(nb)%Sfcprop%lfmassxy(ix)
        temp2d(i,j,idx_opt+18) = IPD_Data(nb)%Sfcprop%rtmassxy(ix)
        temp2d(i,j,idx_opt+19) = IPD_Data(nb)%Sfcprop%stmassxy(ix)
        temp2d(i,j,idx_opt+20) = IPD_Data(nb)%Sfcprop%woodxy(ix)
        temp2d(i,j,idx_opt+21) = IPD_Data(nb)%Sfcprop%stblcpxy(ix)
        temp2d(i,j,idx_opt+22) = IPD_Data(nb)%Sfcprop%fastcpxy(ix)
        temp2d(i,j,idx_opt+23) = IPD_Data(nb)%Sfcprop%xsaixy(ix)
        temp2d(i,j,idx_opt+24) = IPD_Data(nb)%Sfcprop%xlaixy(ix)
        temp2d(i,j,idx_opt+25) = IPD_Data(nb)%Sfcprop%taussxy(ix)
        temp2d(i,j,idx_opt+26) = IPD_Data(nb)%Sfcprop%smcwtdxy(ix)
        temp2d(i,j,idx_opt+27) = IPD_Data(nb)%Sfcprop%deeprechxy(ix)
        temp2d(i,j,idx_opt+28) = IPD_Data(nb)%Sfcprop%rechxy(ix)

        temp2d(i,j,idx_opt+29) = IPD_Data(nb)%Sfcprop%snicexy(ix,-2)
        temp2d(i,j,idx_opt+30) = IPD_Data(nb)%Sfcprop%snicexy(ix,-1)
        temp2d(i,j,idx_opt+31) = IPD_Data(nb)%Sfcprop%snicexy(ix,0)
        temp2d(i,j,idx_opt+32) = IPD_Data(nb)%Sfcprop%snliqxy(ix,-2)
        temp2d(i,j,idx_opt+33) = IPD_Data(nb)%Sfcprop%snliqxy(ix,-1)
        temp2d(i,j,idx_opt+34) = IPD_Data(nb)%Sfcprop%snliqxy(ix,0)
        temp2d(i,j,idx_opt+35) = IPD_Data(nb)%Sfcprop%tsnoxy(ix,-2)
        temp2d(i,j,idx_opt+36) = IPD_Data(nb)%Sfcprop%tsnoxy(ix,-1)
        temp2d(i,j,idx_opt+37) = IPD_Data(nb)%Sfcprop%tsnoxy(ix,0)
        temp2d(i,j,idx_opt+38) = IPD_Data(nb)%Sfcprop%smoiseq(ix,1)
        temp2d(i,j,idx_opt+39) = IPD_Data(nb)%Sfcprop%smoiseq(ix,2)
        temp2d(i,j,idx_opt+40) = IPD_Data(nb)%Sfcprop%smoiseq(ix,3)
        temp2d(i,j,idx_opt+41) = IPD_Data(nb)%Sfcprop%smoiseq(ix,4)
        temp2d(i,j,idx_opt+42) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,-2)
        temp2d(i,j,idx_opt+43) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,-1)
        temp2d(i,j,idx_opt+44) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,0)
        temp2d(i,j,idx_opt+45) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,1)
        temp2d(i,j,idx_opt+46) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,2)
        temp2d(i,j,idx_opt+47) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,3)
        temp2d(i,j,idx_opt+48) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,4)
        idx_opt = 134
       endif

       if (Model%nstf_name(1) > 0) then
         temp2d(i,j,idx_opt) = IPD_Data(nb)%Sfcprop%tref(ix)
         temp2d(i,j,idx_opt+1) = IPD_Data(nb)%Sfcprop%z_c(ix)
         temp2d(i,j,idx_opt+2) = IPD_Data(nb)%Sfcprop%c_0(ix)
         temp2d(i,j,idx_opt+3) = IPD_Data(nb)%Sfcprop%c_d(ix)
         temp2d(i,j,idx_opt+4) = IPD_Data(nb)%Sfcprop%w_0(ix)
         temp2d(i,j,idx_opt+5) = IPD_Data(nb)%Sfcprop%w_d(ix)
         temp2d(i,j,idx_opt+6) = IPD_Data(nb)%Sfcprop%xt(ix)
         temp2d(i,j,idx_opt+7) = IPD_Data(nb)%Sfcprop%xs(ix)
         temp2d(i,j,idx_opt+8) = IPD_Data(nb)%Sfcprop%xu(ix)
         temp2d(i,j,idx_opt+9) = IPD_Data(nb)%Sfcprop%xz(ix)
         temp2d(i,j,idx_opt+10) = IPD_Data(nb)%Sfcprop%zm(ix)
         temp2d(i,j,idx_opt+11) = IPD_Data(nb)%Sfcprop%xtts(ix)
         temp2d(i,j,idx_opt+12) = IPD_Data(nb)%Sfcprop%xzts(ix)
         temp2d(i,j,idx_opt+13) = IPD_Data(nb)%Sfcprop%ifd(ix)
         temp2d(i,j,idx_opt+14) = IPD_Data(nb)%Sfcprop%dt_cool(ix)
         temp2d(i,j,idx_opt+15) = IPD_Data(nb)%Sfcprop%qrain(ix)
       endif

       do l = 1,Model%ntot2d
         temp2d(i,j,nsfcprop2d+l) = IPD_Data(nb)%Tbd%phy_f2d(ix,l)
       enddo

       do l = 1,Model%nctp
         temp2d(i,j,nsfcprop2d+Model%ntot2d+l) = IPD_Data(nb)%Tbd%phy_fctd(ix,l)
       enddo

       temp3dlevsp1(i,j,:, 1) = IPD_Data(nb)%Statein%phii(ix,:)
       temp3dlevsp1(i,j,:, 2) = IPD_Data(nb)%Statein%prsi(ix,:)
       temp3dlevsp1(i,j,:, 3) = IPD_Data(nb)%Statein%prsik(ix,:)

       temp3d(i,j,:, 1) = IPD_Data(nb)%Statein%phil(ix,:)
       temp3d(i,j,:, 2) = IPD_Data(nb)%Statein%prsl(ix,:)
       temp3d(i,j,:, 3) = IPD_Data(nb)%Statein%prslk(ix,:)
       temp3d(i,j,:, 4) = IPD_Data(nb)%Statein%ugrs(ix,:)
       temp3d(i,j,:, 5) = IPD_Data(nb)%Statein%vgrs(ix,:)
       temp3d(i,j,:, 6) = IPD_Data(nb)%Statein%vvl(ix,:)
       temp3d(i,j,:, 7) = IPD_Data(nb)%Statein%tgrs(ix,:)
       temp3d(i,j,:, 8) = IPD_Data(nb)%Stateout%gu0(ix,:)
       temp3d(i,j,:, 9) = IPD_Data(nb)%Stateout%gv0(ix,:)
       temp3d(i,j,:,10) = IPD_Data(nb)%Stateout%gt0(ix,:)
       temp3d(i,j,:,11) = IPD_Data(nb)%Radtend%htrsw(ix,:)
       temp3d(i,j,:,12) = IPD_Data(nb)%Radtend%htrlw(ix,:)
       temp3d(i,j,:,13) = IPD_Data(nb)%Radtend%swhc(ix,:)
       temp3d(i,j,:,14) = IPD_Data(nb)%Radtend%lwhc(ix,:)
       do l = 1,Model%ntot3d
         temp3d(i,j,:,14+l) = IPD_Data(nb)%Tbd%phy_f3d(ix,:,l)
       enddo
       do l = 1,ntr
         temp3d(i,j,:,14+Model%ntot3d+l)     = IPD_Data(nb)%Statein%qgrs(ix,:,l)
         temp3d(i,j,:,14+Model%ntot3d+ntr+l) = IPD_Data(nb)%Stateout%gq0(ix,:,l)
       enddo
     enddo
   enddo

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
   end subroutine FV3GFS_IPD_checksum

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
#ifdef CCPP
  subroutine sfc_prop_restart_read (Sfcprop, Atm_block, Model, fv_domain, warm_start)
#else
  subroutine sfc_prop_restart_read (Sfcprop, Atm_block, Model, fv_domain)
#endif
    !--- interface variable definitions
    type(GFS_sfcprop_type),    intent(inout) :: Sfcprop(:)
    type (block_control_type), intent(in)    :: Atm_block
    type(IPD_control_type),    intent(inout) :: Model
    type (domain2d),           intent(in)    :: fv_domain
#ifdef CCPP
    logical,                   intent(in)    :: warm_start
#endif
    !--- local variables
    integer :: i, j, k, ix, lsoil, num, nb
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: id_restart
    integer :: nvar_o2, nvar_s2m, nvar_s2o, nvar_s3
    integer :: nvar_s2mp, nvar_s3mp,isnow
#ifdef CCPP
    integer :: nvar_s2r
#endif
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p  => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p  => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p1 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p2 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p3 => NULL()
    !--- local variables for sncovr calculation
    integer :: vegtyp
    logical :: mand
    real(kind=kind_phys) :: rsnow, tem
    !--- Noah MP
    integer              :: soiltyp,ns,imon,iter,imn
    real(kind=kind_phys) :: masslai, masssai,snd
    real(kind=kind_phys) :: ddz,expon,aa,bb,smc,func,dfunc,dx
    real(kind=kind_phys) :: bexp, smcmax, smcwlt,dwsat,dksat,psisat

    real(kind=kind_phys), dimension(-2:0) :: dzsno
    real(kind=kind_phys), dimension(-2:4) :: dzsnso

    real(kind=kind_phys), dimension(4), save :: zsoil,dzs
    data dzs   /0.1,0.3,0.6,1.0/
    data zsoil /-0.1,-0.4,-1.0,-2.0/

    
    if (Model%cplflx) then ! needs more variables
      nvar_s2m = 34
    else
      nvar_s2m = 32
    endif
    nvar_o2  = 19
    nvar_s2o = 18
#ifdef CCPP
    if (Model%lsm == Model%lsm_ruc .and. warm_start) then
      if(Model%rdlai) then
        nvar_s2r = 7
      else
        nvar_s2r = 6
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
#else
    nvar_s3  = 3
#endif

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
    if (.not. allocated(oro_name2)) then
    !--- allocate the various containers needed for orography data
      allocate(oro_name2(nvar_o2))
      allocate(oro_var2(nx,ny,nvar_o2))
      oro_var2 = -9999._kind_phys

      oro_name2(1)  = 'stddev'     ! hprime(ix,1)
      oro_name2(2)  = 'convexity'  ! hprime(ix,2)
      oro_name2(3)  = 'oa1'        ! hprime(ix,3)
      oro_name2(4)  = 'oa2'        ! hprime(ix,4)
      oro_name2(5)  = 'oa3'        ! hprime(ix,5)
      oro_name2(6)  = 'oa4'        ! hprime(ix,6)
      oro_name2(7)  = 'ol1'        ! hprime(ix,7)
      oro_name2(8)  = 'ol2'        ! hprime(ix,8)
      oro_name2(9)  = 'ol3'        ! hprime(ix,9)
      oro_name2(10) = 'ol4'        ! hprime(ix,10)
      oro_name2(11) = 'theta'      ! hprime(ix,11)
      oro_name2(12) = 'gamma'      ! hprime(ix,12)
      oro_name2(13) = 'sigma'      ! hprime(ix,13)
      oro_name2(14) = 'elvmax'     ! hprime(ix,14)
      oro_name2(15) = 'orog_filt'  ! oro
      oro_name2(16) = 'orog_raw'   ! oro_uf
      oro_name2(17) = 'land_frac'  ! land fraction [0:1]
      !--- variables below here are optional
      oro_name2(18) = 'lake_frac'  ! lake fraction [0:1]
      oro_name2(19) = 'lake_depth' ! lake depth(m)
      !--- register the 2D fields
      do num = 1,nvar_o2
        var2_p => oro_var2(:,:,num)
        if (trim(oro_name2(num)) == 'lake_frac' .or. trim(oro_name2(num)) == 'lake_depth') then
          id_restart = register_restart_field(Oro_restart, fn_oro, oro_name2(num), var2_p, domain=fv_domain, mandatory=.false.)
        else
          id_restart = register_restart_field(Oro_restart, fn_oro, oro_name2(num), var2_p, domain=fv_domain)
        endif
      enddo
      nullify(var2_p)
    endif

    !--- read the orography restart/data
    call mpp_error(NOTE,'reading topographic/orographic information from INPUT/oro_data.tile*.nc')
    call restore_state(Oro_restart)

    Model%frac_grid = .false.
    !--- copy data into GFS containers
    do nb = 1, Atm_block%nblks
      !--- 2D variables
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1
        !--- stddev
!       Sfcprop(nb)%hprim(ix)     = oro_var2(i,j,1)
        !--- hprime(1:14)
        Sfcprop(nb)%hprime(ix,1)  = oro_var2(i,j,1)
        Sfcprop(nb)%hprime(ix,2)  = oro_var2(i,j,2)
        Sfcprop(nb)%hprime(ix,3)  = oro_var2(i,j,3)
        Sfcprop(nb)%hprime(ix,4)  = oro_var2(i,j,4)
        Sfcprop(nb)%hprime(ix,5)  = oro_var2(i,j,5)
        Sfcprop(nb)%hprime(ix,6)  = oro_var2(i,j,6)
        Sfcprop(nb)%hprime(ix,7)  = oro_var2(i,j,7)
        Sfcprop(nb)%hprime(ix,8)  = oro_var2(i,j,8)
        Sfcprop(nb)%hprime(ix,9)  = oro_var2(i,j,9)
        Sfcprop(nb)%hprime(ix,10) = oro_var2(i,j,10)
        Sfcprop(nb)%hprime(ix,11) = oro_var2(i,j,11)
        Sfcprop(nb)%hprime(ix,12) = oro_var2(i,j,12)
        Sfcprop(nb)%hprime(ix,13) = oro_var2(i,j,13)
        Sfcprop(nb)%hprime(ix,14) = oro_var2(i,j,14)
        !--- oro
        Sfcprop(nb)%oro(ix)       = oro_var2(i,j,15)
        !--- oro_uf
        Sfcprop(nb)%oro_uf(ix)    = oro_var2(i,j,16)
        Sfcprop(nb)%landfrac(ix)  = oro_var2(i,j,17) !land frac [0:1]
        Sfcprop(nb)%lakefrac(ix)  = oro_var2(i,j,18) !lake frac [0:1]
      enddo
    enddo
 
    if (nint(oro_var2(1,1,18)) == -9999._kind_phys) then ! lakefrac doesn't exist in the restart, need to create it
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - will computing lakefrac') 
      Model%frac_grid = .false.
    elseif (Model%frac_grid_off) then
      Model%frac_grid = .false.
    else
      Model%frac_grid = .true.
    endif

    if (Model%me == Model%master ) write(0,*)' resetting Model%frac_grid=',Model%frac_grid

    !--- deallocate containers and free restart container
    deallocate(oro_name2, oro_var2)
    call free_restart_type(Oro_restart)
 
    !--- SURFACE FILE
    if (.not. allocated(sfc_name2)) then
      !--- allocate the various containers needed for restarts
#ifdef CCPP
      allocate(sfc_name2(nvar_s2m+nvar_s2o+nvar_s2mp+nvar_s2r))
      allocate(sfc_name3(nvar_s3+nvar_s3mp))

      allocate(sfc_var2(nx,ny,nvar_s2m+nvar_s2o+nvar_s2mp+nvar_s2r))
      if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp .or. (.not.warm_start)) then
        allocate(sfc_var3(nx,ny,Model%lsoil,nvar_s3))
      else if (Model%lsm == Model%lsm_ruc) then
        allocate(sfc_var3(nx,ny,Model%lsoil_lsm,nvar_s3))
      end if
#else
      allocate(sfc_name2(nvar_s2m+nvar_s2o+nvar_s2mp))
      allocate(sfc_name3(nvar_s3+nvar_s3mp))

      allocate(sfc_var2(nx,ny,nvar_s2m+nvar_s2o+nvar_s2mp))
      allocate(sfc_var3(nx,ny,Model%lsoil,nvar_s3))
#endif
      sfc_var2   = -9999._kind_phys
      sfc_var3   = -9999._kind_phys
!
      if (Model%lsm == Model%lsm_noahmp) then
        allocate(sfc_var3sn(nx,ny,-2:0,4:6))
        allocate(sfc_var3eq(nx,ny,1:4,7:7))
        allocate(sfc_var3zn(nx,ny,-2:4,8:8))
        sfc_var3sn = -9999._kind_phys
        sfc_var3eq = -9999._kind_phys
        sfc_var3zn = -9999._kind_phys
      end if

      !--- names of the 2D variables to save
      sfc_name2(1)  = 'slmsk'
      sfc_name2(2)  = 'tsea'    !tsfc
      sfc_name2(3)  = 'sheleg'  !weasd
      sfc_name2(4)  = 'tg3'
      sfc_name2(5)  = 'zorl'
      sfc_name2(6)  = 'alvsf'
      sfc_name2(7)  = 'alvwf'
      sfc_name2(8)  = 'alnsf'
      sfc_name2(9)  = 'alnwf'
      sfc_name2(10) = 'facsf'
      sfc_name2(11) = 'facwf'
      sfc_name2(12) = 'vfrac'
      sfc_name2(13) = 'canopy'
      sfc_name2(14) = 'f10m'
      sfc_name2(15) = 't2m'
      sfc_name2(16) = 'q2m'
      sfc_name2(17) = 'vtype'
      sfc_name2(18) = 'stype'
      sfc_name2(19) = 'uustar'
      sfc_name2(20) = 'ffmm'
      sfc_name2(21) = 'ffhh'
      sfc_name2(22) = 'hice'
      sfc_name2(23) = 'fice'
      sfc_name2(24) = 'tisfc'
      sfc_name2(25) = 'tprcp'
      sfc_name2(26) = 'srflag'
      sfc_name2(27) = 'snwdph'  !snowd
      sfc_name2(28) = 'shdmin'
      sfc_name2(29) = 'shdmax'
      sfc_name2(30) = 'slope'
      sfc_name2(31) = 'snoalb'
      !--- variables below here are optional
      sfc_name2(32) = 'sncovr'
      if(Model%cplflx) then
        sfc_name2(33) = 'tsfcl'   !temp on land portion of a cell
        sfc_name2(34) = 'zorll'   !zorl on land portion of a cell
      end if

      !--- NSSTM inputs only needed when (nstf_name(1) > 0) .and. (nstf_name(2)) == 0) 
      sfc_name2(nvar_s2m+1)  = 'tref'
      sfc_name2(nvar_s2m+2)  = 'z_c'
      sfc_name2(nvar_s2m+3)  = 'c_0'
      sfc_name2(nvar_s2m+4)  = 'c_d'
      sfc_name2(nvar_s2m+5)  = 'w_0'
      sfc_name2(nvar_s2m+6)  = 'w_d'
      sfc_name2(nvar_s2m+7)  = 'xt'
      sfc_name2(nvar_s2m+8)  = 'xs'
      sfc_name2(nvar_s2m+9)  = 'xu'
      sfc_name2(nvar_s2m+10) = 'xv'
      sfc_name2(nvar_s2m+11) = 'xz'
      sfc_name2(nvar_s2m+12) = 'zm'
      sfc_name2(nvar_s2m+13) = 'xtts'
      sfc_name2(nvar_s2m+14) = 'xzts'
      sfc_name2(nvar_s2m+15) = 'd_conv'
      sfc_name2(nvar_s2m+16) = 'ifd'
      sfc_name2(nvar_s2m+17) = 'dt_cool'
      sfc_name2(nvar_s2m+18) = 'qrain'
!
! Only needed when Noah MP LSM is used - 29 2D
!
      if (Model%lsm == Model%lsm_noahmp) then
        sfc_name2(nvar_s2m+19) = 'snowxy'
        sfc_name2(nvar_s2m+20) = 'tvxy'
        sfc_name2(nvar_s2m+21) = 'tgxy'
        sfc_name2(nvar_s2m+22) = 'canicexy'
        sfc_name2(nvar_s2m+23) = 'canliqxy'
        sfc_name2(nvar_s2m+24) = 'eahxy'
        sfc_name2(nvar_s2m+25) = 'tahxy'
        sfc_name2(nvar_s2m+26) = 'cmxy'
        sfc_name2(nvar_s2m+27) = 'chxy'
        sfc_name2(nvar_s2m+28) = 'fwetxy'
        sfc_name2(nvar_s2m+29) = 'sneqvoxy'
        sfc_name2(nvar_s2m+30) = 'alboldxy'
        sfc_name2(nvar_s2m+31) = 'qsnowxy'
        sfc_name2(nvar_s2m+32) = 'wslakexy'
        sfc_name2(nvar_s2m+33) = 'zwtxy'
        sfc_name2(nvar_s2m+34) = 'waxy'
        sfc_name2(nvar_s2m+35) = 'wtxy'
        sfc_name2(nvar_s2m+36) = 'lfmassxy'
        sfc_name2(nvar_s2m+37) = 'rtmassxy'
        sfc_name2(nvar_s2m+38) = 'stmassxy'
        sfc_name2(nvar_s2m+39) = 'woodxy'
        sfc_name2(nvar_s2m+40) = 'stblcpxy'
        sfc_name2(nvar_s2m+41) = 'fastcpxy'
        sfc_name2(nvar_s2m+42) = 'xsaixy'
        sfc_name2(nvar_s2m+43) = 'xlaixy'
        sfc_name2(nvar_s2m+44) = 'taussxy'
        sfc_name2(nvar_s2m+45) = 'smcwtdxy'
        sfc_name2(nvar_s2m+46) = 'deeprechxy'
        sfc_name2(nvar_s2m+47) = 'rechxy'
#ifdef CCPP
      else if (Model%lsm == Model%lsm_ruc .and. warm_start) then
        sfc_name2(nvar_s2m+19) = 'wetness'
        sfc_name2(nvar_s2m+20) = 'clw_surf'
        sfc_name2(nvar_s2m+21) = 'qwv_surf'
        sfc_name2(nvar_s2m+22) = 'tsnow'
        sfc_name2(nvar_s2m+23) = 'snowfall_acc'
        sfc_name2(nvar_s2m+24) = 'swe_snowfall_acc'
        if (Model%rdlai) then
          sfc_name2(nvar_s2m+25) = 'lai'
        endif
      else if (Model%lsm == Model%lsm_ruc .and. Model%rdlai) then
        sfc_name2(nvar_s2m+19) = 'lai'
#endif
      endif

      !--- register the 2D fields
      do num = 1,nvar_s2m
        var2_p => sfc_var2(:,:,num)
        if (trim(sfc_name2(num)) == 'sncovr'.or.trim(sfc_name2(num)) == 'tsfcl'.or.trim(sfc_name2(num)) == 'zorll') then
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=fv_domain, mandatory=.false.)
        else
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=fv_domain)
        endif
      enddo


      if (Model%nstf_name(1) > 0) then
        mand = .false.
        if (Model%nstf_name(2) == 0) mand = .true.
        do num = nvar_s2m+1,nvar_s2m+nvar_s2o
          var2_p => sfc_var2(:,:,num)
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=fv_domain, mandatory=mand)
        enddo
      endif
#ifdef CCPP
      if (Model%lsm == Model%lsm_ruc) then ! nvar_s2mp = 0
        do num = nvar_s2m+nvar_s2o+1, nvar_s2m+nvar_s2o+nvar_s2r
          var2_p => sfc_var2(:,:,num)
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=fv_domain)
        enddo
      endif ! mp/ruc
#endif
! Noah MP register only necessary only lsm = 2, not necessary has values
      if (nvar_s2mp > 0) then
        mand = .false.
        do num = nvar_s2m+nvar_s2o+1,nvar_s2m+nvar_s2o+nvar_s2mp
          var2_p => sfc_var2(:,:,num)
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=fv_domain, mandatory=mand)
        enddo
      endif ! noahmp

      nullify(var2_p)
    endif  ! if not allocated

 
#ifdef CCPP
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
#else
      !--- names of the 3D variables to save
    sfc_name3(1) = 'stc'
    sfc_name3(2) = 'smc'
    sfc_name3(3) = 'slc'
      !--- Noah MP
    if (Model%lsm == Model%lsm_noahmp) then
      sfc_name3(4) = 'snicexy'
      sfc_name3(5) = 'snliqxy'
      sfc_name3(6) = 'tsnoxy'
      sfc_name3(7) = 'smoiseq'
      sfc_name3(8) = 'zsnsoxy'
    endif
#endif
      !--- register the 3D fields
    do num = 1,nvar_s3
      var3_p => sfc_var3(:,:,:,num)
      id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name3(num), var3_p, domain=fv_domain)
    enddo
    if (Model%lsm == Model%lsm_noahmp) then
      mand = .false.
      do num = nvar_s3+1,nvar_s3+3
        var3_p1 => sfc_var3sn(:,:,:,num)
        id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name3(num), var3_p1, domain=fv_domain,mandatory=mand)
      enddo

      var3_p2 => sfc_var3eq(:,:,:,7)
      id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name3(7), var3_p2, domain=fv_domain,mandatory=mand)

      var3_p3 => sfc_var3zn(:,:,:,8)
      id_restart = register_restart_fIeld(Sfc_restart, fn_srf, sfc_name3(8), var3_p3, domain=fv_domain,mandatory=mand)

      nullify(var3_p1)
      nullify(var3_p2)
      nullify(var3_p3)

    endif   !mp

    nullify(var3_p)

!--- Noah MP define arbitrary value (number layers of snow) to indicate
!coldstart(sfcfile doesn't include noah mp fields) or not

    if (Model%lsm == Model%lsm_noahmp) then
      sfc_var2(1,1,nvar_s2m+19) = -66666.
    endif

    !--- read the surface restart/data
    call mpp_error(NOTE,'reading surface properties data from INPUT/sfc_data.tile*.nc')
    call restore_state(Sfc_restart)

!   write(0,*)' stype read in min,max=',minval(sfc_var2(:,:,18)),maxval(sfc_var2(:,:,18))
!   write(0,*)' sfc_var2=',sfc_var2(:,:,12)

    !--- place the data into the block GFS containers
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1

!--- 2D variables
!    ------------
        Sfcprop(nb)%slmsk(ix)  = sfc_var2(i,j,1)    !--- slmsk
        Sfcprop(nb)%tsfco(ix)  = sfc_var2(i,j,2)    !--- tsfc (tsea in sfc file)
        Sfcprop(nb)%weasd(ix)  = sfc_var2(i,j,3)    !--- weasd (sheleg in sfc file)
        Sfcprop(nb)%tg3(ix)    = sfc_var2(i,j,4)    !--- tg3
        Sfcprop(nb)%zorlo(ix)  = sfc_var2(i,j,5)    !--- zorl on ocean
        Sfcprop(nb)%alvsf(ix)  = sfc_var2(i,j,6)    !--- alvsf
        Sfcprop(nb)%alvwf(ix)  = sfc_var2(i,j,7)    !--- alvwf
        Sfcprop(nb)%alnsf(ix)  = sfc_var2(i,j,8)    !--- alnsf
        Sfcprop(nb)%alnwf(ix)  = sfc_var2(i,j,9)    !--- alnwf
        Sfcprop(nb)%facsf(ix)  = sfc_var2(i,j,10)   !--- facsf
        Sfcprop(nb)%facwf(ix)  = sfc_var2(i,j,11)   !--- facwf
        Sfcprop(nb)%vfrac(ix)  = sfc_var2(i,j,12)   !--- vfrac
        Sfcprop(nb)%canopy(ix) = sfc_var2(i,j,13)   !--- canopy
        Sfcprop(nb)%f10m(ix)   = sfc_var2(i,j,14)   !--- f10m
        Sfcprop(nb)%t2m(ix)    = sfc_var2(i,j,15)   !--- t2m
        Sfcprop(nb)%q2m(ix)    = sfc_var2(i,j,16)   !--- q2m
        Sfcprop(nb)%vtype(ix)  = sfc_var2(i,j,17)   !--- vtype
        Sfcprop(nb)%stype(ix)  = sfc_var2(i,j,18)   !--- stype
        Sfcprop(nb)%uustar(ix) = sfc_var2(i,j,19)   !--- uustar
        Sfcprop(nb)%ffmm(ix)   = sfc_var2(i,j,20)   !--- ffmm
        Sfcprop(nb)%ffhh(ix)   = sfc_var2(i,j,21)   !--- ffhh
        Sfcprop(nb)%hice(ix)   = sfc_var2(i,j,22)   !--- hice
        Sfcprop(nb)%fice(ix)   = sfc_var2(i,j,23)   !--- fice
        Sfcprop(nb)%tisfc(ix)  = sfc_var2(i,j,24)   !--- tisfc
        Sfcprop(nb)%tprcp(ix)  = sfc_var2(i,j,25)   !--- tprcp
        Sfcprop(nb)%srflag(ix) = sfc_var2(i,j,26)   !--- srflag
        Sfcprop(nb)%snowd(ix)  = sfc_var2(i,j,27)   !--- snowd (snwdph in the file)
        Sfcprop(nb)%shdmin(ix) = sfc_var2(i,j,28)   !--- shdmin
        Sfcprop(nb)%shdmax(ix) = sfc_var2(i,j,29)   !--- shdmax
        Sfcprop(nb)%slope(ix)  = sfc_var2(i,j,30)   !--- slope
        Sfcprop(nb)%snoalb(ix) = sfc_var2(i,j,31)   !--- snoalb
        Sfcprop(nb)%sncovr(ix) = sfc_var2(i,j,32)   !--- sncovr
        if(Model%cplflx) then
          Sfcprop(nb)%tsfcl(ix)  = sfc_var2(i,j,33) !--- sfcl  (temp on land portion of a cell)
          Sfcprop(nb)%zorll(ix)  = sfc_var2(i,j,34) !--- zorll (zorl on land portion of a cell)
        end if

        if(Model%frac_grid) then ! obtain slmsk from landfrac
!! next 5 lines are temporary till lake model is available
          if (Sfcprop(nb)%lakefrac(ix) > 0.0) then
            Sfcprop(nb)%lakefrac(ix) = nint(Sfcprop(nb)%lakefrac(ix))
            Sfcprop(nb)%landfrac(ix) = 1.-Sfcprop(nb)%lakefrac(ix)
            if (Sfcprop(nb)%lakefrac(ix) == 0) Sfcprop(nb)%fice(ix)=0.
          end if 
          Sfcprop(nb)%slmsk(ix) = ceiling(Sfcprop(nb)%landfrac(ix))
          if (Sfcprop(nb)%fice(ix) > 0. .and. Sfcprop(nb)%landfrac(ix)==0.) Sfcprop(nb)%slmsk(ix) = 2 ! land dominates ice if co-exist
        else ! obtain landfrac from slmsk
          if (Sfcprop(nb)%slmsk(ix) > 1.9) then
            Sfcprop(nb)%landfrac(ix) = 0.0
          else
            Sfcprop(nb)%landfrac(ix) = Sfcprop(nb)%slmsk(ix)
          endif
        end if

        if (Sfcprop(nb)%lakefrac(ix) > 0.0) then
          Sfcprop(nb)%oceanfrac(ix) = 0.0 ! lake & ocean don't coexist in a cell
          if (Sfcprop(nb)%fice(ix) < Model%min_lakeice) Sfcprop(nb)%fice(ix) = 0.
        else
          Sfcprop(nb)%oceanfrac(ix) = 1.0 - Sfcprop(nb)%landfrac(ix)
          if (Sfcprop(nb)%fice(ix) < Model%min_seaice) Sfcprop(nb)%fice(ix) = 0.
        endif
        !
        !--- NSSTM variables
        if ((Model%nstf_name(1) > 0) .and. (Model%nstf_name(2) == 1)) then
          !--- nsstm tref
          Sfcprop(nb)%tref(ix)    = Sfcprop(nb)%tsfco(ix)
          Sfcprop(nb)%xz(ix)      = 30.0d0
        endif
        if ((Model%nstf_name(1) > 0) .and. (Model%nstf_name(2) == 0)) then
          Sfcprop(nb)%tref(ix)    = sfc_var2(i,j,nvar_s2m+1)  !--- nsstm tref
          Sfcprop(nb)%z_c(ix)     = sfc_var2(i,j,nvar_s2m+2)  !--- nsstm z_c
          Sfcprop(nb)%c_0(ix)     = sfc_var2(i,j,nvar_s2m+3)  !--- nsstm c_0
          Sfcprop(nb)%c_d(ix)     = sfc_var2(i,j,nvar_s2m+4)  !--- nsstm c_d
          Sfcprop(nb)%w_0(ix)     = sfc_var2(i,j,nvar_s2m+5)  !--- nsstm w_0
          Sfcprop(nb)%w_d(ix)     = sfc_var2(i,j,nvar_s2m+6)  !--- nsstm w_d
          Sfcprop(nb)%xt(ix)      = sfc_var2(i,j,nvar_s2m+7)  !--- nsstm xt
          Sfcprop(nb)%xs(ix)      = sfc_var2(i,j,nvar_s2m+8)  !--- nsstm xs
          Sfcprop(nb)%xu(ix)      = sfc_var2(i,j,nvar_s2m+9)  !--- nsstm xu
          Sfcprop(nb)%xv(ix)      = sfc_var2(i,j,nvar_s2m+10) !--- nsstm xv
          Sfcprop(nb)%xz(ix)      = sfc_var2(i,j,nvar_s2m+11) !--- nsstm xz
          Sfcprop(nb)%zm(ix)      = sfc_var2(i,j,nvar_s2m+12) !--- nsstm zm
          Sfcprop(nb)%xtts(ix)    = sfc_var2(i,j,nvar_s2m+13) !--- nsstm xtts
          Sfcprop(nb)%xzts(ix)    = sfc_var2(i,j,nvar_s2m+14) !--- nsstm xzts
          Sfcprop(nb)%d_conv(ix)  = sfc_var2(i,j,nvar_s2m+15) !--- nsstm d_conv
          Sfcprop(nb)%ifd(ix)     = sfc_var2(i,j,nvar_s2m+16) !--- nsstm ifd
          Sfcprop(nb)%dt_cool(ix) = sfc_var2(i,j,nvar_s2m+17) !--- nsstm dt_cool
          Sfcprop(nb)%qrain(ix)   = sfc_var2(i,j,nvar_s2m+18) !--- nsstm qrain
        endif
#ifdef CCPP
        if (Model%lsm == Model%lsm_ruc .and. warm_start) then
          !--- Extra RUC variables
          Sfcprop(nb)%wetness(ix)    = sfc_var2(i,j,nvar_s2m+19)
          Sfcprop(nb)%clw_surf(ix)   = sfc_var2(i,j,nvar_s2m+20)
          Sfcprop(nb)%qwv_surf(ix)   = sfc_var2(i,j,nvar_s2m+21)
          Sfcprop(nb)%tsnow(ix)      = sfc_var2(i,j,nvar_s2m+22)
          Sfcprop(nb)%snowfallac(ix) = sfc_var2(i,j,nvar_s2m+23)
          Sfcprop(nb)%acsnow(ix)     = sfc_var2(i,j,nvar_s2m+24)
          if (Model%rdlai) then
            Sfcprop(nb)%xlaixy(ix)   = sfc_var2(i,j,nvar_s2m+25)
          endif
        else if (Model%lsm == Model%lsm_ruc .and. Model%rdlai) then
          Sfcprop(nb)%xlaixy(ix)     = sfc_var2(i,j,nvar_s2m+19)
        elseif (Model%lsm == Model%lsm_noahmp) then
          !--- Extra Noah MP variables
#else
! Noah MP
! -------
        if (Model%lsm == Model%lsm_noahmp) then
#endif
          Sfcprop(nb)%snowxy(ix)     = sfc_var2(i,j,nvar_s2m+19)
          Sfcprop(nb)%tvxy(ix)       = sfc_var2(i,j,nvar_s2m+20)
          Sfcprop(nb)%tgxy(ix)       = sfc_var2(i,j,nvar_s2m+21)
          Sfcprop(nb)%canicexy(ix)   = sfc_var2(i,j,nvar_s2m+22)
          Sfcprop(nb)%canliqxy(ix)   = sfc_var2(i,j,nvar_s2m+23)
          Sfcprop(nb)%eahxy(ix)      = sfc_var2(i,j,nvar_s2m+24)
          Sfcprop(nb)%tahxy(ix)      = sfc_var2(i,j,nvar_s2m+25)
          Sfcprop(nb)%cmxy(ix)       = sfc_var2(i,j,nvar_s2m+26)
          Sfcprop(nb)%chxy(ix)       = sfc_var2(i,j,nvar_s2m+27)
          Sfcprop(nb)%fwetxy(ix)     = sfc_var2(i,j,nvar_s2m+28)
          Sfcprop(nb)%sneqvoxy(ix)   = sfc_var2(i,j,nvar_s2m+29)
          Sfcprop(nb)%alboldxy(ix)   = sfc_var2(i,j,nvar_s2m+30)
          Sfcprop(nb)%qsnowxy(ix)    = sfc_var2(i,j,nvar_s2m+31)
          Sfcprop(nb)%wslakexy(ix)   = sfc_var2(i,j,nvar_s2m+32)
          Sfcprop(nb)%zwtxy(ix)      = sfc_var2(i,j,nvar_s2m+33)
          Sfcprop(nb)%waxy(ix)       = sfc_var2(i,j,nvar_s2m+34)
          Sfcprop(nb)%wtxy(ix)       = sfc_var2(i,j,nvar_s2m+35)
          Sfcprop(nb)%lfmassxy(ix)   = sfc_var2(i,j,nvar_s2m+36)
          Sfcprop(nb)%rtmassxy(ix)   = sfc_var2(i,j,nvar_s2m+37)
          Sfcprop(nb)%stmassxy(ix)   = sfc_var2(i,j,nvar_s2m+38)
          Sfcprop(nb)%woodxy(ix)     = sfc_var2(i,j,nvar_s2m+39)
          Sfcprop(nb)%stblcpxy(ix)   = sfc_var2(i,j,nvar_s2m+40)
          Sfcprop(nb)%fastcpxy(ix)   = sfc_var2(i,j,nvar_s2m+41)
          Sfcprop(nb)%xsaixy(ix)     = sfc_var2(i,j,nvar_s2m+42)
          Sfcprop(nb)%xlaixy(ix)     = sfc_var2(i,j,nvar_s2m+43)
          Sfcprop(nb)%taussxy(ix)    = sfc_var2(i,j,nvar_s2m+44)
          Sfcprop(nb)%smcwtdxy(ix)   = sfc_var2(i,j,nvar_s2m+45)
          Sfcprop(nb)%deeprechxy(ix) = sfc_var2(i,j,nvar_s2m+46)
          Sfcprop(nb)%rechxy(ix)     = sfc_var2(i,j,nvar_s2m+47)
        endif

#ifdef CCPP
        if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp .or. (.not.warm_start)) then
          !--- 3D variables
          do lsoil = 1,Model%lsoil
            Sfcprop(nb)%stc(ix,lsoil) = sfc_var3(i,j,lsoil,1)   !--- stc
            Sfcprop(nb)%smc(ix,lsoil) = sfc_var3(i,j,lsoil,2)   !--- smc
            Sfcprop(nb)%slc(ix,lsoil) = sfc_var3(i,j,lsoil,3)   !--- slc
          enddo

          if (Model%lsm == Model%lsm_noahmp) then
            do lsoil = -2, 0
              Sfcprop(nb)%snicexy(ix,lsoil) = sfc_var3sn(i,j,lsoil,4)
              Sfcprop(nb)%snliqxy(ix,lsoil) = sfc_var3sn(i,j,lsoil,5)
              Sfcprop(nb)%tsnoxy(ix,lsoil)  = sfc_var3sn(i,j,lsoil,6)
            enddo 

            do lsoil = 1, 4
              Sfcprop(nb)%smoiseq(ix,lsoil)  = sfc_var3eq(i,j,lsoil,7)
            enddo 

            do lsoil = -2, 4
              Sfcprop(nb)%zsnsoxy(ix,lsoil)  = sfc_var3zn(i,j,lsoil,8)
            enddo 
          endif

        else if (Model%lsm == Model%lsm_ruc) then
          !--- 3D variables
          do lsoil = 1,Model%lsoil_lsm
            Sfcprop(nb)%tslb(ix,lsoil)        = sfc_var3(i,j,lsoil,1) !--- tslb
            Sfcprop(nb)%smois(ix,lsoil)       = sfc_var3(i,j,lsoil,2) !--- smois
            Sfcprop(nb)%sh2o(ix,lsoil)        = sfc_var3(i,j,lsoil,3) !--- sh2o
            Sfcprop(nb)%keepsmfr(ix,lsoil)    = sfc_var3(i,j,lsoil,4) !--- keepsmfr
            Sfcprop(nb)%flag_frsoil(ix,lsoil) = sfc_var3(i,j,lsoil,5) !--- flag_frsoil
          enddo
        end if
#else
        !--- 3D variables
        do lsoil = 1,Model%lsoil
          Sfcprop(nb)%stc(ix,lsoil) = sfc_var3(i,j,lsoil,1)   !--- stc
          Sfcprop(nb)%smc(ix,lsoil) = sfc_var3(i,j,lsoil,2)   !--- smc
          Sfcprop(nb)%slc(ix,lsoil) = sfc_var3(i,j,lsoil,3)   !--- slc
        enddo

        if (Model%lsm == Model%lsm_noahmp) then
          do lsoil = -2, 0
            Sfcprop(nb)%snicexy(ix,lsoil) = sfc_var3sn(i,j,lsoil,4)
            Sfcprop(nb)%snliqxy(ix,lsoil) = sfc_var3sn(i,j,lsoil,5)
            Sfcprop(nb)%tsnoxy(ix,lsoil)  = sfc_var3sn(i,j,lsoil,6)
          enddo 

          do lsoil = 1, 4
            Sfcprop(nb)%smoiseq(ix,lsoil)  = sfc_var3eq(i,j,lsoil,7)
          enddo 

          do lsoil = -2, 4
            Sfcprop(nb)%zsnsoxy(ix,lsoil)  = sfc_var3zn(i,j,lsoil,8)
          enddo 
        endif
#endif

      enddo   !ix
    enddo    !nb
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

#ifdef CCPP
    ! Calculating sncovr does NOT belong into an I/O routine!
    ! TODO: move to physics and stop building namelist_soilveg/set_soilveg
    ! in the FV3/non-CCPP physics when the CCPP-enabled executable is built.
#endif
!#ifndef CCPP
    !--- if sncovr does not exist in the restart, need to create it
    if (nint(sfc_var2(1,1,32)) == -9999) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing sncovr') 
      !--- compute sncovr from existing variables
      !--- code taken directly from read_fix.f
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          Sfcprop(nb)%sncovr(ix) = 0.0
          if (Sfcprop(nb)%landfrac(ix) >= drythresh .or. Sfcprop(nb)%fice(ix) >= Model%min_seaice) then
            vegtyp = Sfcprop(nb)%vtype(ix)
            if (vegtyp == 0) vegtyp = 7
            rsnow  = 0.001*Sfcprop(nb)%weasd(ix)/snupx(vegtyp)
            if (0.001*Sfcprop(nb)%weasd(ix) < snupx(vegtyp)) then
              Sfcprop(nb)%sncovr(ix) = 1.0 - (exp(-salp_data*rsnow) - rsnow*exp(-salp_data))
            else
              Sfcprop(nb)%sncovr(ix) = 1.0
            endif
          endif
        enddo
      enddo
    endif

    if(Model%cplflx .or. Model%frac_grid) then
      if (nint(sfc_var2(1,1,33)) == -9999) then
        if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing tsfcl')
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            Sfcprop(nb)%tsfcl(ix) = Sfcprop(nb)%tsfco(ix) !--- compute tsfcl from existing variables
          enddo
        enddo
      endif

      if (nint(sfc_var2(1,1,34)) == -9999) then
        if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing zorll')
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            Sfcprop(nb)%zorll(ix) = Sfcprop(nb)%zorlo(ix) !--- compute zorll from existing variables
          enddo
        enddo
      endif
    endif

!#endif

    if(Model%frac_grid) then ! 3-way composite
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          Sfcprop(nb)%tsfco(ix) = max(con_tice, Sfcprop(nb)%tsfco(ix))
          tem = (1.-Sfcprop(nb)%landfrac(ix)) * Sfcprop(nb)%fice(ix) ! tem = ice fraction wrt whole cell
          Sfcprop(nb)%zorl(ix) = Sfcprop(nb)%zorll(ix) * Sfcprop(nb)%landfrac(ix) &
                               + Sfcprop(nb)%zorll(ix) * tem &     !zorl ice = zorl land
                               + Sfcprop(nb)%zorlo(ix) * (1.-Sfcprop(nb)%landfrac(ix)-tem)
          Sfcprop(nb)%tsfc(ix) = Sfcprop(nb)%tsfcl(ix) * Sfcprop(nb)%landfrac(ix) &
                               + Sfcprop(nb)%tisfc(ix) * tem &
                               + Sfcprop(nb)%tsfco(ix) * (1.-Sfcprop(nb)%landfrac(ix)-tem)
        enddo
      enddo
    else
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
      !--- specify tsfcl/zorll from existing variable tsfco/zorlo
          Sfcprop(nb)%tsfcl(ix) = Sfcprop(nb)%tsfco(ix)
          Sfcprop(nb)%zorll(ix) = Sfcprop(nb)%zorlo(ix)
          Sfcprop(nb)%zorl(ix)  = Sfcprop(nb)%zorlo(ix)
          Sfcprop(nb)%tsfc(ix)  = Sfcprop(nb)%tsfco(ix)
        enddo
      enddo
    endif ! if (Model%frac_grid)

    if (Model%lsm == Model%lsm_noahmp) then 
      if (nint(sfc_var2(1,1,nvar_s2m+19)) == -66666) then
        if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver:: - Cold start Noah MP ')

        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)

            Sfcprop(nb)%tvxy(ix)     = missing_value
            Sfcprop(nb)%tgxy(ix)     = missing_value
            Sfcprop(nb)%tahxy(ix)    = missing_value
            Sfcprop(nb)%canicexy(ix) = missing_value
            Sfcprop(nb)%canliqxy(ix) = missing_value
            Sfcprop(nb)%eahxy(ix)    = missing_value
            Sfcprop(nb)%cmxy(ix)     = missing_value
            Sfcprop(nb)%chxy(ix)     = missing_value
            Sfcprop(nb)%fwetxy(ix)   = missing_value
            Sfcprop(nb)%sneqvoxy(ix) = missing_value
            Sfcprop(nb)%alboldxy(ix) = missing_value
            Sfcprop(nb)%qsnowxy(ix)  = missing_value
            Sfcprop(nb)%wslakexy     = missing_value
            Sfcprop(nb)%taussxy      = missing_value
            Sfcprop(nb)%waxy(ix)     = missing_value
            Sfcprop(nb)%wtxy(ix)     = missing_value
            Sfcprop(nb)%zwtxy(ix)    = missing_value
            Sfcprop(nb)%xlaixy(ix)   = missing_value
            Sfcprop(nb)%xsaixy(ix)   = missing_value

            Sfcprop(nb)%lfmassxy(ix) = missing_value
            Sfcprop(nb)%stmassxy(ix) = missing_value
            Sfcprop(nb)%rtmassxy(ix) = missing_value
            Sfcprop(nb)%woodxy(ix)   = missing_value
            Sfcprop(nb)%stblcpxy(ix) = missing_value
            Sfcprop(nb)%fastcpxy(ix) = missing_value
            Sfcprop(nb)%smcwtdxy(ix) = missing_value
            Sfcprop(nb)%deeprechxy(ix) = missing_value
            Sfcprop(nb)%rechxy(ix)     = missing_value

            Sfcprop(nb)%snowxy (ix)   = missing_value
            Sfcprop(nb)%snicexy(ix, -2:0) = missing_value
            Sfcprop(nb)%snliqxy(ix, -2:0) = missing_value
            Sfcprop(nb)%tsnoxy (ix, -2:0) = missing_value
            Sfcprop(nb)%smoiseq(ix,  1:4) = missing_value
            Sfcprop(nb)%zsnsoxy(ix, -2:4) = missing_value

            if (Sfcprop(nb)%landfrac(ix) >= drythresh) then

              Sfcprop(nb)%tvxy(ix)     = Sfcprop(nb)%tsfcl(ix)
              Sfcprop(nb)%tgxy(ix)     = Sfcprop(nb)%tsfcl(ix)
              Sfcprop(nb)%tahxy(ix)    = Sfcprop(nb)%tsfcl(ix)

              if (Sfcprop(nb)%snowd(ix) > 0.01 .and. Sfcprop(nb)%tsfcl(ix) > 273.15 ) Sfcprop(nb)%tvxy  = 273.15
              if (Sfcprop(nb)%snowd(ix) > 0.01 .and. Sfcprop(nb)%tsfcl(ix) > 273.15 ) Sfcprop(nb)%tgxy  = 273.15
              if (Sfcprop(nb)%snowd(ix) > 0.01 .and. Sfcprop(nb)%tsfcl(ix) > 273.15 ) Sfcprop(nb)%tahxy = 273.15
  
              Sfcprop(nb)%canicexy(ix) = 0.0
              Sfcprop(nb)%canliqxy(ix) = Sfcprop(nb)%canopy(ix)

              Sfcprop(nb)%eahxy(ix)    = 2000.0

!      eahxy = psfc*qv/(0.622+qv); qv is mixing ratio, converted from sepcific
!      humidity specific humidity /(1.0 - specific humidity)

              Sfcprop(nb)%cmxy(ix)     = 0.0
              Sfcprop(nb)%chxy(ix)     = 0.0
              Sfcprop(nb)%fwetxy(ix)   = 0.0
              Sfcprop(nb)%sneqvoxy(ix) = Sfcprop(nb)%weasd(ix)     ! mm
              Sfcprop(nb)%alboldxy(ix) = 0.65
              Sfcprop(nb)%qsnowxy(ix)  = 0.0

!           if (Sfcprop(nb)%srflag(ix) > 0.001) Sfcprop(nb)%qsnowxy(ix) = Sfcprop(nb)%tprcp(ix)/Model%dtp
! already set to 0.0
              Sfcprop(nb)%wslakexy     = 0.0
              Sfcprop(nb)%taussxy      = 0.0


              Sfcprop(nb)%waxy(ix)     = 4900.0
              Sfcprop(nb)%wtxy(ix)     = Sfcprop(nb)%waxy(ix)
              Sfcprop(nb)%zwtxy(ix)    = (25.0 + 2.0) - Sfcprop(nb)%waxy(ix) / 1000.0 /0.2
!
              vegtyp                   = Sfcprop(nb)%vtype(ix)
              if (vegtyp == 0) vegtyp = 7
              imn                      = Model%idate(2)

              if ((vegtyp == isbarren_table) .or. (vegtyp == isice_table) .or.  (vegtyp == isurban_table) .or. (vegtyp == iswater_table)) then

                Sfcprop(nb)%xlaixy(ix)   = 0.0
                Sfcprop(nb)%xsaixy(ix)   = 0.0

                Sfcprop(nb)%lfmassxy(ix) = 0.0
                Sfcprop(nb)%stmassxy(ix) = 0.0
                Sfcprop(nb)%rtmassxy(ix) = 0.0

                Sfcprop(nb)%woodxy   (ix) = 0.0       
                Sfcprop(nb)%stblcpxy (ix) = 0.0      
                Sfcprop(nb)%fastcpxy (ix) = 0.0     

              else

!             print *, 'vegtyp', vegtyp
!             print *, 'imn', imn
!             print *, 'xlaixy', Sfcprop(nb)%xlaixy(ix) 

                Sfcprop(nb)%xlaixy(ix)   = max(laim_table(vegtyp, imn),0.05)
!             Sfcprop(nb)%xsaixy(ix)   = max(saim_table(vegtyp, imn),0.05)
                Sfcprop(nb)%xsaixy(ix)   = max(Sfcprop(nb)%xlaixy(ix)*0.1,0.05)

                masslai                  = 1000.0 / max(sla_table(vegtyp),1.0)
                Sfcprop(nb)%lfmassxy(ix) = Sfcprop(nb)%xlaixy(ix)*masslai
                masssai                  = 1000.0 / 3.0
                Sfcprop(nb)%stmassxy(ix) = Sfcprop(nb)%xsaixy(ix)* masssai

                Sfcprop(nb)%rtmassxy(ix) = 500.0      

                Sfcprop(nb)%woodxy  (ix) = 500.0       
                Sfcprop(nb)%stblcpxy(ix) = 1000.0      
                Sfcprop(nb)%fastcpxy(ix) = 1000.0     

              endif  ! non urban ...

              if ( vegtyp == isice_table )  then
                do lsoil = 1,Model%lsoil
                  Sfcprop(nb)%stc(ix,lsoil) = min(Sfcprop(nb)%stc(ix,lsoil),min(Sfcprop(nb)%tg3(ix),263.15))
                  Sfcprop(nb)%smc(ix,lsoil) = 1
                  Sfcprop(nb)%slc(ix,lsoil) = 0
                enddo
              endif

              snd   = Sfcprop(nb)%snowd(ix)/1000.0  ! go to m from snwdph

              if (Sfcprop(nb)%weasd(ix) /= 0.0 .and. snd == 0.0 ) then
                snd = Sfcprop(nb)%weasd(ix)/1000.0
              endif

              if (vegtyp == 15) then                      ! land ice in MODIS/IGBP
                if ( Sfcprop(nb)%weasd(ix) < 0.1) then
                  Sfcprop(nb)%weasd(ix) = 0.1
                  snd                   = 0.01
                endif
              endif

            if (snd < 0.025 ) then
              Sfcprop(nb)%snowxy(ix)   = 0.0
              dzsno(-2:0)              = 0.0
            elseif (snd >= 0.025 .and. snd <= 0.05 ) then
              Sfcprop(nb)%snowxy(ix)   = -1.0
              dzsno(0)                 = snd
            elseif (snd > 0.05 .and. snd <= 0.10 ) then
              Sfcprop(nb)%snowxy(ix)   = -2.0
              dzsno(-1)                = 0.5*snd
              dzsno(0)                 = 0.5*snd
            elseif (snd > 0.10 .and. snd <= 0.25 ) then
              Sfcprop(nb)%snowxy(ix)   = -2.0
              dzsno(-1)                = 0.05
              dzsno(0)                 = snd - 0.05
            elseif (snd > 0.25 .and. snd <= 0.45 ) then
              Sfcprop(nb)%snowxy(ix)   = -3.0
              dzsno(-2)                = 0.05
              dzsno(-1)                = 0.5*(snd-0.05)
              dzsno(0)                 = 0.5*(snd-0.05)
            elseif (snd > 0.45) then 
              Sfcprop(nb)%snowxy(ix)   = -3.0
              dzsno(-2)                = 0.05
              dzsno(-1)                = 0.20
              dzsno(0)                 = snd - 0.05 - 0.20
            else
              call mpp_error(FATAL, 'problem with the logic assigning snow layers.') 
            endif

! Now we have the snowxy field
! snice + snliq + tsno allocation and compute them from what we have
             
!
              Sfcprop(nb)%tsnoxy(ix,-2:0)  = 0.0
              Sfcprop(nb)%snicexy(ix,-2:0) = 0.0
              Sfcprop(nb)%snliqxy(ix,-2:0) = 0.0
              Sfcprop(nb)%zsnsoxy(ix,-2:4) = 0.0

              isnow = nint(Sfcprop(nb)%snowxy(ix))+1    ! snowxy <=0.0, dzsno >= 0.0

              do ns = isnow , 0
                Sfcprop(nb)%tsnoxy(ix,ns)  = Sfcprop(nb)%tgxy(ix)
                Sfcprop(nb)%snliqxy(ix,ns) = 0.0
                Sfcprop(nb)%snicexy(ix,ns) = 1.00 * dzsno(ns) * Sfcprop(nb)%weasd(ix)/snd
              enddo
!
!zsnsoxy, all negative ?
!
              do ns = isnow, 0
                dzsnso(ns) = -dzsno(ns)
              enddo

              do ns = 1 , 4
                dzsnso(ns) = -dzs(ns)
              enddo
!
! Assign to zsnsoxy
!
              Sfcprop(nb)%zsnsoxy(ix,isnow) = dzsnso(isnow)
              do ns = isnow+1,4
                Sfcprop(nb)%zsnsoxy(ix,ns) = Sfcprop(nb)%zsnsoxy(ix,ns-1) + dzsnso(ns)
              enddo
 
!
! smoiseq
! Init water table related quantities here
!
              soiltyp  = Sfcprop(nb)%stype(ix)

              if (soiltyp /= 0) then
                bexp   = bexp_table(soiltyp)
                smcmax = smcmax_table(soiltyp)
                smcwlt = smcwlt_table(soiltyp)
                dwsat  = dwsat_table(soiltyp)
                dksat  = dksat_table(soiltyp)
                psisat = -psisat_table(soiltyp)
              endif

              if (vegtyp == isurban_table) then
                smcmax = 0.45
                smcwlt = 0.40
              endif

              if ((bexp > 0.0) .and. (smcmax > 0.0) .and. (-psisat > 0.0 )) then
                do ns = 1, Model%lsoil          
                  if ( ns == 1 )then
                    ddz = -zsoil(ns+1) * 0.5
                  elseif ( ns < Model%lsoil ) then
                    ddz = ( zsoil(ns-1) - zsoil(ns+1) ) * 0.5
                  else
                    ddz = zsoil(ns-1) - zsoil(ns)
                  endif
!
! Use newton-raphson method to find eq soil moisture
!
                expon = bexp + 1.
                aa    = dwsat / ddz
                bb    = dksat / smcmax ** expon

                smc = 0.5 * smcmax

                do iter = 1, 100
                  func  = (smc - smcmax) * aa +  bb * smc ** expon
                  dfunc = aa + bb * expon * smc ** bexp
                  dx    = func / dfunc
                  smc   = smc - dx
                  if ( abs (dx) < 1.e-6) exit
                enddo                               ! iteration
                Sfcprop(nb)%smoiseq(ix,ns) = min(max(smc,1.e-4),smcmax*0.99)
              enddo                                 ! ddz soil layer
            else                                    ! bexp <= 0.0 
              Sfcprop(nb)%smoiseq(ix,1:4) = smcmax
              endif                                   ! end the bexp condition
!
              Sfcprop(nb)%smcwtdxy(ix)   = smcmax
              Sfcprop(nb)%deeprechxy(ix) = 0.0
              Sfcprop(nb)%rechxy(ix)     = 0.0
 
            endif !end if slmsk>0.01 (land only)

          enddo ! ix
        enddo  ! nb
      endif
    endif   !if Noah MP cold start ends

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
    type(IPD_control_type),      intent(in) :: Model
    type(domain2d),              intent(in) :: fv_domain
    character(len=32), optional, intent(in) :: timestamp
    !--- local variables
    integer :: i, j, k, nb, ix, lsoil, num
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: id_restart
    integer :: nvar2m, nvar2o, nvar3
    integer :: nvar2mp, nvar3mp
#ifdef CCPP
    integer :: nvar2r
#endif
    logical :: mand
    character(len=32) :: fn_srf = 'sfc_data.nc'
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p  => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p  => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p1 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p2 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p3 => NULL()

    if (Model%cplflx) then ! needs more variables
      nvar2m = 34
    else
      nvar2m = 32
    endif
    nvar2o = 18
#ifdef CCPP
    if (Model%lsm == Model%lsm_ruc) then
      if (Model%rdlai) then
        nvar2r = 7
      else
        nvar2r = 6
      endif
      nvar3  = 5
    else
      nvar2r = 0
      nvar3  = 3
    endif
#else
    nvar3  = 3
#endif
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

#ifdef CCPP
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
          call free_restart_type(Sfc_restart)
        end if
      end if
    end if
#endif

    if (.not. allocated(sfc_name2)) then
      !--- allocate the various containers needed for restarts
#ifdef CCPP
      allocate(sfc_name2(nvar2m+nvar2o+nvar2mp+nvar2r))
      allocate(sfc_name3(nvar3+nvar3mp))
      allocate(sfc_var2(nx,ny,nvar2m+nvar2o+nvar2mp+nvar2r))
      if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp) then
        allocate(sfc_var3(nx,ny,Model%lsoil,nvar3))
      elseif (Model%lsm == Model%lsm_ruc) then
        allocate(sfc_var3(nx,ny,Model%lsoil_lsm,nvar3))
      endif
#else
      allocate(sfc_name2(nvar2m+nvar2o+nvar2mp))
      allocate(sfc_name3(nvar3+nvar3mp))
      allocate(sfc_var2(nx,ny,nvar2m+nvar2o+nvar2mp))
      allocate(sfc_var3(nx,ny,Model%lsoil,nvar3))
#endif
      sfc_var2   = -9999._kind_phys
      sfc_var3   = -9999._kind_phys
      if (Model%lsm == Model%lsm_noahmp) then
        allocate(sfc_var3sn(nx,ny,-2:0,4:6))
        allocate(sfc_var3eq(nx,ny,1:4,7:7))
        allocate(sfc_var3zn(nx,ny,-2:4,8:8))

        sfc_var3sn = -9999._kind_phys
        sfc_var3eq = -9999._kind_phys
        sfc_var3zn = -9999._kind_phys
      endif


    !--- names of the 2D variables to save
      sfc_name2(1)  = 'slmsk'
      sfc_name2(2)  = 'tsea'    !tsfc
      sfc_name2(3)  = 'sheleg'  !weasd
      sfc_name2(4)  = 'tg3'
      sfc_name2(5)  = 'zorl'
      sfc_name2(6)  = 'alvsf'
      sfc_name2(7)  = 'alvwf'
      sfc_name2(8)  = 'alnsf'
      sfc_name2(9)  = 'alnwf'
      sfc_name2(10) = 'facsf'
      sfc_name2(11) = 'facwf'
      sfc_name2(12) = 'vfrac'
      sfc_name2(13) = 'canopy'
      sfc_name2(14) = 'f10m'
      sfc_name2(15) = 't2m'
      sfc_name2(16) = 'q2m'
      sfc_name2(17) = 'vtype'
      sfc_name2(18) = 'stype'
      sfc_name2(19) = 'uustar'
      sfc_name2(20) = 'ffmm'
      sfc_name2(21) = 'ffhh'
      sfc_name2(22) = 'hice'
      sfc_name2(23) = 'fice'
      sfc_name2(24) = 'tisfc'
      sfc_name2(25) = 'tprcp'
      sfc_name2(26) = 'srflag'
      sfc_name2(27) = 'snwdph'  !snowd
      sfc_name2(28) = 'shdmin'
      sfc_name2(29) = 'shdmax'
      sfc_name2(30) = 'slope'
      sfc_name2(31) = 'snoalb'
    !--- variables below here are optional
      sfc_name2(32) = 'sncovr'
      if (Model%cplflx) then
        sfc_name2(33) = 'tsfcl'   !temp on land portion of a cell
        sfc_name2(34) = 'zorll'   !zorl on land portion of a cell
      end if
    !--- NSSTM inputs only needed when (nstf_name(1) > 0) .and. (nstf_name(2)) == 0)
      sfc_name2(nvar2m+1)  = 'tref'
      sfc_name2(nvar2m+2)  = 'z_c'
      sfc_name2(nvar2m+3)  = 'c_0'
      sfc_name2(nvar2m+4)  = 'c_d'
      sfc_name2(nvar2m+5)  = 'w_0'
      sfc_name2(nvar2m+6)  = 'w_d'
      sfc_name2(nvar2m+7)  = 'xt'
      sfc_name2(nvar2m+8)  = 'xs'
      sfc_name2(nvar2m+9)  = 'xu'
      sfc_name2(nvar2m+10) = 'xv'
      sfc_name2(nvar2m+11) = 'xz'
      sfc_name2(nvar2m+12) = 'zm'
      sfc_name2(nvar2m+13) = 'xtts'
      sfc_name2(nvar2m+14) = 'xzts'
      sfc_name2(nvar2m+15) = 'd_conv'
      sfc_name2(nvar2m+16) = 'ifd'
      sfc_name2(nvar2m+17) = 'dt_cool'
      sfc_name2(nvar2m+18) = 'qrain'
#ifdef CCPP
      if (Model%lsm == Model%lsm_ruc) then
        sfc_name2(nvar2m+19) = 'wetness'
        sfc_name2(nvar2m+20) = 'clw_surf'
        sfc_name2(nvar2m+21) = 'qwv_surf'
        sfc_name2(nvar2m+22) = 'tsnow'
        sfc_name2(nvar2m+23) = 'snowfall_acc'
        sfc_name2(nvar2m+24) = 'swe_snowfall_acc'
        if (Model%rdlai) then
          sfc_name2(nvar2m+25) = 'lai'
        endif
      else if(Model%lsm == Model%lsm_noahmp) then
#else
! Only needed when Noah MP LSM is used - 29 2D
      if(Model%lsm == Model%lsm_noahmp) then
#endif
        sfc_name2(nvar2m+19) = 'snowxy'
        sfc_name2(nvar2m+20) = 'tvxy'
        sfc_name2(nvar2m+21) = 'tgxy'
        sfc_name2(nvar2m+22) = 'canicexy'
        sfc_name2(nvar2m+23) = 'canliqxy'
        sfc_name2(nvar2m+24) = 'eahxy'
        sfc_name2(nvar2m+25) = 'tahxy'
        sfc_name2(nvar2m+26) = 'cmxy'
        sfc_name2(nvar2m+27) = 'chxy'
        sfc_name2(nvar2m+28) = 'fwetxy'
        sfc_name2(nvar2m+29) = 'sneqvoxy'
        sfc_name2(nvar2m+30) = 'alboldxy'
        sfc_name2(nvar2m+31) = 'qsnowxy'
        sfc_name2(nvar2m+32) = 'wslakexy'
        sfc_name2(nvar2m+33) = 'zwtxy'
        sfc_name2(nvar2m+34) = 'waxy'
        sfc_name2(nvar2m+35) = 'wtxy'
        sfc_name2(nvar2m+36) = 'lfmassxy'
        sfc_name2(nvar2m+37) = 'rtmassxy'
        sfc_name2(nvar2m+38) = 'stmassxy'
        sfc_name2(nvar2m+39) = 'woodxy'
        sfc_name2(nvar2m+40) = 'stblcpxy'
        sfc_name2(nvar2m+41) = 'fastcpxy'
        sfc_name2(nvar2m+42) = 'xsaixy'
        sfc_name2(nvar2m+43) = 'xlaixy'
        sfc_name2(nvar2m+44) = 'taussxy'
        sfc_name2(nvar2m+45) = 'smcwtdxy'
        sfc_name2(nvar2m+46) = 'deeprechxy'
        sfc_name2(nvar2m+47) = 'rechxy'
      endif
 
    !--- register the 2D fields
      do num = 1,nvar2m
        var2_p => sfc_var2(:,:,num)
        if (trim(sfc_name2(num)) == 'sncovr'.or.trim(sfc_name2(num)) == 'tsfcl'.or.trim(sfc_name2(num)) == 'zorll') then
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=fv_domain, mandatory=.false.)
        else
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=fv_domain)
        endif
      enddo
      if (Model%nstf_name(1) > 0) then
        mand = .false.
        if (Model%nstf_name(2) ==0) mand = .true.
        do num = nvar2m+1,nvar2m+nvar2o
          var2_p => sfc_var2(:,:,num)
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=fv_domain, mandatory=mand)
        enddo
      endif
#ifdef CCPP
      if (Model%lsm == Model%lsm_ruc) then ! nvar2mp =0 
        do num = nvar2m+nvar2o+1, nvar2m+nvar2o+nvar2r
          var2_p => sfc_var2(:,:,num)
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=fv_domain)
        enddo
      else if (Model%lsm == Model%lsm_noahmp) then ! nvar2r =0 
#else
      if (Model%lsm == Model%lsm_noahmp) then
#endif
        mand = .true.                  ! actually should be true since it is after cold start
        do num = nvar2m+nvar2o+1,nvar2m+nvar2o+nvar2mp
          var2_p => sfc_var2(:,:,num)
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=fv_domain, mandatory=mand)
        enddo
      endif
      nullify(var2_p)

#ifdef CCPP
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
#else
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
#endif

      !--- register the 3D fields
      do num = 1,nvar3
        var3_p => sfc_var3(:,:,:,num)
        id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name3(num), var3_p, domain=fv_domain)
      enddo
      nullify(var3_p)

      if (Model%lsm == Model%lsm_noahmp) then
        mand = .true.
        do num = nvar3+1,nvar3+3
          var3_p1 => sfc_var3sn(:,:,:,num)
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name3(num), var3_p1, domain=fv_domain,mandatory=mand)
        enddo

        var3_p2 => sfc_var3eq(:,:,:,7)
        id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name3(7), var3_p2, domain=fv_domain,mandatory=mand)

        var3_p3 => sfc_var3zn(:,:,:,8)
        id_restart = register_restart_fIeld(Sfc_restart, fn_srf, sfc_name3(8), var3_p3, domain=fv_domain,mandatory=mand)

        nullify(var3_p1)
        nullify(var3_p2)
        nullify(var3_p3)
      endif ! lsm = lsm_noahmp
    endif

   
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        !--- 2D variables
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1
        sfc_var2(i,j,1)  = Sfcprop(nb)%slmsk(ix) !--- slmsk
        sfc_var2(i,j,2)  = Sfcprop(nb)%tsfc(ix)  !--- tsfc (tsea in sfc file)
        sfc_var2(i,j,3)  = Sfcprop(nb)%weasd(ix) !--- weasd (sheleg in sfc file)
        sfc_var2(i,j,4)  = Sfcprop(nb)%tg3(ix)   !--- tg3
        sfc_var2(i,j,5)  = Sfcprop(nb)%zorl(ix)  !--- zorl
        sfc_var2(i,j,6)  = Sfcprop(nb)%alvsf(ix) !--- alvsf
        sfc_var2(i,j,7)  = Sfcprop(nb)%alvwf(ix) !--- alvwf
        sfc_var2(i,j,8)  = Sfcprop(nb)%alnsf(ix) !--- alnsf
        sfc_var2(i,j,9)  = Sfcprop(nb)%alnwf(ix) !--- alnwf
        sfc_var2(i,j,10) = Sfcprop(nb)%facsf(ix) !--- facsf
        sfc_var2(i,j,11) = Sfcprop(nb)%facwf(ix) !--- facwf
        sfc_var2(i,j,12) = Sfcprop(nb)%vfrac(ix) !--- vfrac
        sfc_var2(i,j,13) = Sfcprop(nb)%canopy(ix)!--- canopy
        sfc_var2(i,j,14) = Sfcprop(nb)%f10m(ix)  !--- f10m
        sfc_var2(i,j,15) = Sfcprop(nb)%t2m(ix)   !--- t2m
        sfc_var2(i,j,16) = Sfcprop(nb)%q2m(ix)   !--- q2m
        sfc_var2(i,j,17) = Sfcprop(nb)%vtype(ix) !--- vtype
        sfc_var2(i,j,18) = Sfcprop(nb)%stype(ix) !--- stype
        sfc_var2(i,j,19) = Sfcprop(nb)%uustar(ix)!--- uustar
        sfc_var2(i,j,20) = Sfcprop(nb)%ffmm(ix)  !--- ffmm
        sfc_var2(i,j,21) = Sfcprop(nb)%ffhh(ix)  !--- ffhh
        sfc_var2(i,j,22) = Sfcprop(nb)%hice(ix)  !--- hice
        sfc_var2(i,j,23) = Sfcprop(nb)%fice(ix)  !--- fice
        sfc_var2(i,j,24) = Sfcprop(nb)%tisfc(ix) !--- tisfc
        sfc_var2(i,j,25) = Sfcprop(nb)%tprcp(ix) !--- tprcp
        sfc_var2(i,j,26) = Sfcprop(nb)%srflag(ix)!--- srflag
        sfc_var2(i,j,27) = Sfcprop(nb)%snowd(ix) !--- snowd (snwdph in the file)
        sfc_var2(i,j,28) = Sfcprop(nb)%shdmin(ix)!--- shdmin
        sfc_var2(i,j,29) = Sfcprop(nb)%shdmax(ix)!--- shdmax
        sfc_var2(i,j,30) = Sfcprop(nb)%slope(ix) !--- slope
        sfc_var2(i,j,31) = Sfcprop(nb)%snoalb(ix)!--- snoalb
        sfc_var2(i,j,32) = Sfcprop(nb)%sncovr(ix)!--- sncovr
        if (Model%cplflx) then
          sfc_var2(i,j,33) = Sfcprop(nb)%tsfcl(ix) !--- tsfcl (temp on land)
          sfc_var2(i,j,34) = Sfcprop(nb)%zorll(ix) !--- zorll (zorl on land)
        end if
        !--- NSSTM variables
        if (Model%nstf_name(1) > 0) then
          sfc_var2(i,j,nvar2m+1) = Sfcprop(nb)%tref(ix)    !--- nsstm tref
          sfc_var2(i,j,nvar2m+2) = Sfcprop(nb)%z_c(ix)     !--- nsstm z_c
          sfc_var2(i,j,nvar2m+3) = Sfcprop(nb)%c_0(ix)     !--- nsstm c_0
          sfc_var2(i,j,nvar2m+4) = Sfcprop(nb)%c_d(ix)     !--- nsstm c_d
          sfc_var2(i,j,nvar2m+5) = Sfcprop(nb)%w_0(ix)     !--- nsstm w_0
          sfc_var2(i,j,nvar2m+6) = Sfcprop(nb)%w_d(ix)     !--- nsstm w_d
          sfc_var2(i,j,nvar2m+7) = Sfcprop(nb)%xt(ix)      !--- nsstm xt
          sfc_var2(i,j,nvar2m+8) = Sfcprop(nb)%xs(ix)      !--- nsstm xs
          sfc_var2(i,j,nvar2m+9) = Sfcprop(nb)%xu(ix)      !--- nsstm xu
          sfc_var2(i,j,nvar2m+10) = Sfcprop(nb)%xv(ix)     !--- nsstm xv
          sfc_var2(i,j,nvar2m+11) = Sfcprop(nb)%xz(ix)     !--- nsstm xz
          sfc_var2(i,j,nvar2m+12) = Sfcprop(nb)%zm(ix)     !--- nsstm zm
          sfc_var2(i,j,nvar2m+13) = Sfcprop(nb)%xtts(ix)   !--- nsstm xtts
          sfc_var2(i,j,nvar2m+14) = Sfcprop(nb)%xzts(ix)   !--- nsstm xzts
          sfc_var2(i,j,nvar2m+15) = Sfcprop(nb)%d_conv(ix) !--- nsstm d_conv
          sfc_var2(i,j,nvar2m+16) = Sfcprop(nb)%ifd(ix)    !--- nsstm ifd
          sfc_var2(i,j,nvar2m+17) = Sfcprop(nb)%dt_cool(ix)!--- nsstm dt_cool
          sfc_var2(i,j,nvar2m+18) = Sfcprop(nb)%qrain(ix)  !--- nsstm qrain
        endif
#ifdef CCPP
        if (Model%lsm == Model%lsm_ruc) then
          !--- Extra RUC variables
          sfc_var2(i,j,nvar2m+19) = Sfcprop(nb)%wetness(ix)
          sfc_var2(i,j,nvar2m+20) = Sfcprop(nb)%clw_surf(ix)
          sfc_var2(i,j,nvar2m+21) = Sfcprop(nb)%qwv_surf(ix)
          sfc_var2(i,j,nvar2m+22) = Sfcprop(nb)%tsnow(ix)
          sfc_var2(i,j,nvar2m+23) = Sfcprop(nb)%snowfallac(ix)
          sfc_var2(i,j,nvar2m+24) = Sfcprop(nb)%acsnow(ix)
          if (Model%rdlai) then
            sfc_var2(i,j,nvar2m+25) = Sfcprop(nb)%xlaixy(ix)
          endif
        else if (Model%lsm == Model%lsm_noahmp) then

#else
! Noah MP
        if (Model%lsm == Model%lsm_noahmp) then
#endif

          sfc_var2(i,j,nvar2m+19) = Sfcprop(nb)%snowxy(ix)
          sfc_var2(i,j,nvar2m+20) = Sfcprop(nb)%tvxy(ix)
          sfc_var2(i,j,nvar2m+21) = Sfcprop(nb)%tgxy(ix)
          sfc_var2(i,j,nvar2m+22) = Sfcprop(nb)%canicexy(ix)
          sfc_var2(i,j,nvar2m+23) = Sfcprop(nb)%canliqxy(ix)
          sfc_var2(i,j,nvar2m+24) = Sfcprop(nb)%eahxy(ix)
          sfc_var2(i,j,nvar2m+25) = Sfcprop(nb)%tahxy(ix)
          sfc_var2(i,j,nvar2m+26) = Sfcprop(nb)%cmxy(ix)
          sfc_var2(i,j,nvar2m+27) = Sfcprop(nb)%chxy(ix)
          sfc_var2(i,j,nvar2m+28) = Sfcprop(nb)%fwetxy(ix)
          sfc_var2(i,j,nvar2m+29) = Sfcprop(nb)%sneqvoxy(ix)
          sfc_var2(i,j,nvar2m+30) = Sfcprop(nb)%alboldxy(ix)
          sfc_var2(i,j,nvar2m+31) = Sfcprop(nb)%qsnowxy(ix)
          sfc_var2(i,j,nvar2m+32) = Sfcprop(nb)%wslakexy(ix)
          sfc_var2(i,j,nvar2m+33) = Sfcprop(nb)%zwtxy(ix)
          sfc_var2(i,j,nvar2m+34) = Sfcprop(nb)%waxy(ix)
          sfc_var2(i,j,nvar2m+35) = Sfcprop(nb)%wtxy(ix)
          sfc_var2(i,j,nvar2m+36) = Sfcprop(nb)%lfmassxy(ix)
          sfc_var2(i,j,nvar2m+37) = Sfcprop(nb)%rtmassxy(ix)
          sfc_var2(i,j,nvar2m+38) = Sfcprop(nb)%stmassxy(ix)
          sfc_var2(i,j,nvar2m+39) = Sfcprop(nb)%woodxy(ix)
          sfc_var2(i,j,nvar2m+40) = Sfcprop(nb)%stblcpxy(ix)
          sfc_var2(i,j,nvar2m+41) = Sfcprop(nb)%fastcpxy(ix)
          sfc_var2(i,j,nvar2m+42) = Sfcprop(nb)%xsaixy(ix)
          sfc_var2(i,j,nvar2m+43) = Sfcprop(nb)%xlaixy(ix)
          sfc_var2(i,j,nvar2m+44) = Sfcprop(nb)%taussxy(ix)
          sfc_var2(i,j,nvar2m+45) = Sfcprop(nb)%smcwtdxy(ix)
          sfc_var2(i,j,nvar2m+46) = Sfcprop(nb)%deeprechxy(ix)
          sfc_var2(i,j,nvar2m+47) = Sfcprop(nb)%rechxy(ix)
        endif

#ifdef CCPP
        if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp) then
          !--- 3D variables
          do lsoil = 1,Model%lsoil
            sfc_var3(i,j,lsoil,1) = Sfcprop(nb)%stc(ix,lsoil) !--- stc
            sfc_var3(i,j,lsoil,2) = Sfcprop(nb)%smc(ix,lsoil) !--- smc
            sfc_var3(i,j,lsoil,3) = Sfcprop(nb)%slc(ix,lsoil) !--- slc
          enddo
! 5 Noah MP 3D
          if (Model%lsm == Model%lsm_noahmp) then

             do lsoil = -2,0
              sfc_var3sn(i,j,lsoil,4) = Sfcprop(nb)%snicexy(ix,lsoil)
              sfc_var3sn(i,j,lsoil,5) = Sfcprop(nb)%snliqxy(ix,lsoil)
              sfc_var3sn(i,j,lsoil,6) = Sfcprop(nb)%tsnoxy(ix,lsoil)
            enddo

            do lsoil = 1,Model%lsoil
              sfc_var3eq(i,j,lsoil,7)  = Sfcprop(nb)%smoiseq(ix,lsoil)
            enddo

            do lsoil = -2,4
              sfc_var3zn(i,j,lsoil,8)  = Sfcprop(nb)%zsnsoxy(ix,lsoil)
            enddo

          endif  ! Noah MP
        else if (Model%lsm == Model%lsm_ruc) then
          !--- 3D variables
          do lsoil = 1,Model%lsoil_lsm
            sfc_var3(i,j,lsoil,1) = Sfcprop(nb)%tslb(ix,lsoil)         !--- tslb  = stc
            sfc_var3(i,j,lsoil,2) = Sfcprop(nb)%smois(ix,lsoil)        !--- smois = smc
            sfc_var3(i,j,lsoil,3) = Sfcprop(nb)%sh2o(ix,lsoil)         !--- sh2o  = slc
            sfc_var3(i,j,lsoil,4) = Sfcprop(nb)%keepsmfr(ix,lsoil)     !--- keepsmfr
            sfc_var3(i,j,lsoil,5) = Sfcprop(nb)%flag_frsoil(ix,lsoil)  !--- flag_frsoil
          enddo
        end if
#else
        !--- 3D variables
        do lsoil = 1,Model%lsoil
          sfc_var3(i,j,lsoil,1) = Sfcprop(nb)%stc(ix,lsoil) !--- stc
          sfc_var3(i,j,lsoil,2) = Sfcprop(nb)%smc(ix,lsoil) !--- smc
          sfc_var3(i,j,lsoil,3) = Sfcprop(nb)%slc(ix,lsoil) !--- slc
        enddo
! 5 Noah MP 3D
        if (Model%lsm == Model%lsm_noahmp) then

          do lsoil = -2,0
            sfc_var3sn(i,j,lsoil,4) = Sfcprop(nb)%snicexy(ix,lsoil)
            sfc_var3sn(i,j,lsoil,5) = Sfcprop(nb)%snliqxy(ix,lsoil)
            sfc_var3sn(i,j,lsoil,6) = Sfcprop(nb)%tsnoxy(ix,lsoil)
          enddo

          do lsoil = 1,Model%lsoil
            sfc_var3eq(i,j,lsoil,7)  = Sfcprop(nb)%smoiseq(ix,lsoil)
          enddo

          do lsoil = -2,4
            sfc_var3zn(i,j,lsoil,8)  = Sfcprop(nb)%zsnsoxy(ix,lsoil)
          enddo

        endif  ! Noah MP
#endif
      enddo
    enddo

    call save_restart(Sfc_restart, timestamp)

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
  subroutine phys_restart_read (IPD_Restart, Atm_block, Model, fv_domain)
    !--- interface variable definitions
    type(IPD_restart_type),      intent(in) :: IPD_Restart
    type(block_control_type),    intent(in) :: Atm_block
    type(IPD_control_type),      intent(in) :: Model
    type(domain2d),              intent(in) :: fv_domain
    !--- local variables
    integer :: i, j, k, nb, ix, num
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: id_restart
    integer :: nvar2d, nvar3d, fdiag, ldiag
    character(len=64) :: fname
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p => NULL()


    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)
    nvar2d = IPD_Restart%num2d
    nvar3d = IPD_Restart%num3d
    fdiag  = IPD_Restart%fdiag
    ldiag  = IPD_Restart%ldiag
 
    !--- register the restart fields
    if (.not. allocated(phy_var2)) then
      allocate (phy_var2(nx,ny,nvar2d))
      allocate (phy_var3(nx,ny,npz,nvar3d))
      phy_var2 = 0.0_kind_phys
      phy_var3 = 0.0_kind_phys
      
      do num = 1,nvar2d
        var2_p => phy_var2(:,:,num)
        id_restart = register_restart_field (Phy_restart, fn_phy, trim(IPD_Restart%name2d(num)), &
                                             var2_p, domain=fv_domain, mandatory=.false.)
      enddo
      do num = 1,nvar3d
        var3_p => phy_var3(:,:,:,num)
        id_restart = register_restart_field (Phy_restart, fn_phy, trim(IPD_restart%name3d(num)), &
                                             var3_p, domain=fv_domain, mandatory=.false.)
      enddo
      nullify(var2_p)
      nullify(var3_p)
    endif

    fname = 'INPUT/'//trim(fn_phy)
    if (file_exist(fname)) then
      !--- read the surface restart/data
      call mpp_error(NOTE,'reading physics restart data from INPUT/phy_data.tile*.nc')
      call restore_state(Phy_restart)
    else
      call mpp_error(NOTE,'No physics restarts - cold starting physical parameterizations')
      return
    endif
 
    !--- place the data into the block GFS containers
    !--- phy_var* variables
    do num = 1,nvar2d
      do nb = 1,Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)            
          i = Atm_block%index(nb)%ii(ix) - isc + 1
          j = Atm_block%index(nb)%jj(ix) - jsc + 1
          IPD_Restart%data(nb,num)%var2p(ix) = phy_var2(i,j,num)
        enddo
      enddo
    enddo
    !-- if restart from init time, reset accumulated diag fields
    if( Model%phour < 1.e-7) then
      do num = fdiag,ldiag
        do nb = 1,Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            i = Atm_block%index(nb)%ii(ix) - isc + 1
            j = Atm_block%index(nb)%jj(ix) - jsc + 1
            IPD_Restart%data(nb,num)%var2p(ix) = 0.
          enddo
        enddo 
      enddo
    endif
    do num = 1,nvar3d
      do nb = 1,Atm_block%nblks
        do k=1,npz
          do ix = 1, Atm_block%blksz(nb)            
            i = Atm_block%index(nb)%ii(ix) - isc + 1
            j = Atm_block%index(nb)%jj(ix) - jsc + 1
            IPD_Restart%data(nb,num)%var3p(ix,k) = phy_var3(i,j,k,num)
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
  subroutine phys_restart_write (IPD_Restart, Atm_block, Model, fv_domain, timestamp)
    !--- interface variable definitions
    type(IPD_restart_type),      intent(in) :: IPD_Restart
    type(block_control_type),    intent(in) :: Atm_block
    type(IPD_control_type),      intent(in) :: Model
    type(domain2d),              intent(in) :: fv_domain
    character(len=32), optional, intent(in) :: timestamp
    !--- local variables
    integer :: i, j, k, nb, ix, num
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: id_restart
    integer :: nvar2d, nvar3d
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p => NULL()


    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)
    nvar2d = IPD_Restart%num2d
    nvar3d = IPD_Restart%num3d

    !--- register the restart fields 
    if (.not. allocated(phy_var2)) then
      allocate (phy_var2(nx,ny,nvar2d))
      allocate (phy_var3(nx,ny,npz,nvar3d))
      phy_var2 = 0.0_kind_phys
      phy_var3 = 0.0_kind_phys
      
      do num = 1,nvar2d
        var2_p => phy_var2(:,:,num)
        id_restart = register_restart_field (Phy_restart, fn_phy, trim(IPD_Restart%name2d(num)), &
                                             var2_p, domain=fv_domain, mandatory=.false.)
      enddo
      do num = 1,nvar3d
        var3_p => phy_var3(:,:,:,num)
        id_restart = register_restart_field (Phy_restart, fn_phy, trim(IPD_restart%name3d(num)), &
                                             var3_p, domain=fv_domain, mandatory=.false.)
      enddo
      nullify(var2_p)
      nullify(var3_p)
    endif

    !--- 2D variables
    do num = 1,nvar2d
      do nb = 1,Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)            
          i = Atm_block%index(nb)%ii(ix) - isc + 1
          j = Atm_block%index(nb)%jj(ix) - jsc + 1
          phy_var2(i,j,num) = IPD_Restart%data(nb,num)%var2p(ix)
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
            phy_var3(i,j,k,num) = IPD_Restart%data(nb,num)%var3p(ix,k)
          enddo
        enddo
      enddo
    enddo

    call save_restart(Phy_restart, timestamp)

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
    type(IPD_diag_type),       intent(inout) :: Diag(:)
    type(time_type),           intent(in)    :: Time
    type (block_control_type), intent(in)    :: Atm_block
    type(IPD_control_type),    intent(in)    :: Model
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
    ncld   = Model%ncld
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
    buffer_phys_bl = 0.
    buffer_phys_nb = 0.
    buffer_phys_windvect = 0.
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
    type(IPD_diag_type),       intent(in) :: diag(:)
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
    real(kind=kind_phys), dimension(nx*ny)      :: var2p
    real(kind=kind_phys), dimension(nx*ny,levs) :: var3p
    real(kind=kind_phys), dimension(nx,ny)      :: var2
    real(kind=kind_phys), dimension(nx,ny,levs) :: var3
    real(kind=kind_phys) :: rdt, rtime_int, rtime_intfull, lcnvfac
    real(kind=kind_phys) :: rtime_radsw, rtime_radlw
    logical :: used

     nblks         = atm_block%nblks
     rdt           = 1.0d0/dt
     rtime_int     = 1.0d0/time_int
     rtime_intfull = 1.0d0/time_intfull
     rtime_radsw   = 1.0d0/time_radsw
     rtime_radlw   = 1.0d0/time_radlw

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
!           used=send_data(Diag(idx)%id, var2, Time)
!           print *,'in phys, after store_data, idx=',idx,' var=', trim(Diag(idx)%name)
           call store_data(Diag(idx)%id, var2, Time, idx, Diag(idx)%intpl_method, Diag(idx)%name)
!           if(trim(Diag(idx)%name) == 'totprcp_ave' ) print *,'in gfs_io, totprcp=',Diag(idx)%data(1)%var2(1:3), &
!             ' lcnvfac=', lcnvfac
         elseif (Diag(idx)%axes == 3) then
         !---
         !--- skipping other 3D variables with the following else statement
         !---
         if(mpp_pe()==mpp_root_pe())print *,'in,fv3gfs_io. 3D fields, idx=',idx,'varname=',trim(diag(idx)%name), &
             'lcnvfac=',lcnvfac, 'levo=',levo,'nx=',nx,'ny=',ny
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
           !--- dt3dt variables
           do num = 1,6
             write(xtra,'(i1)') num
             if (trim(Diag(idx)%name) == 'dt3dt_'//trim(xtra)) then
               var3(1:nx,1:ny,1:levs) = RESHAPE(Gfs_diag%dt3dt(1:ngptc,levs:1:-1,num:num), (/nx,ny,levs/))
               used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1) 
             endif
           enddo
           !--- dq3dt variables
           do num = 1,5+Mdl_parms%pl_coeff
             write(xtra,'(i1)') num
             if (trim(Diag(idx)%name) == 'dq3dt_'//trim(xtra)) then
               var3(1:nx,1:ny,1:levs) = RESHAPE(Gfs_diag%dq3dt(1:ngptc,levs:1-1,num:num), (/nx,ny,levs/))
               used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1) 
             endif
           enddo
           !--- du3dt and dv3dt variables
           do num = 1,4
             write(xtra,'(i1)') num
             if (trim(Diag(idx)%name) == 'du3dt_'//trim(xtra)) then
               var3(1:nx,1:ny,1:levs) = RESHAPE(Gfs_diag%du3dt(1:ngptc,levs:1:-1,num:num), (/nx,ny,levs/))
               used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1) 
             endif
             if (trim(Diag(idx)%name) == 'dv3dt_'//trim(xtra)) then
               var3(1:nx,1:ny,1:levs) = RESHAPE(Gfs_diag%dv3dt(1:ngptc,levs:1:-1,num:num), (/nx,ny,levs/))
               used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1) 
             endif
           enddo
           if (trim(Diag(idx)%name) == 'dqdt_v') then
             var3(1:nx,1:ny,1:levs) = RESHAPE(Gfs_diag%dqdt_v(1:ngptc,levs:1:-1), (/nx,ny,levs/))
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1) 
           endif
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
           if (trim(Diag(idx)%name) == 'do3mr_dt') then
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
    real(kind=kind_phys), intent(in)    :: work(ieco-isco+1,jeco-jsco+1)
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
            uwork     = 0.0
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
    real(kind=kind_phys), intent(in)    :: work(ieco-isco+1,jeco-jsco+1,levo)
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
            uwork3d   = 0.
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

 subroutine fv_phys_bundle_setup(Diag, axes, phys_bundle, fcst_grid, quilting, nbdlphys)
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
   type(IPD_diag_type),intent(in)              :: Diag(:)
   integer, intent(in)                         :: axes(:)
   type(ESMF_FieldBundle),intent(inout)        :: phys_bundle(:)
   type(ESMF_Grid),intent(inout)               :: fcst_grid
   logical,intent(in)                          :: quilting
   integer, intent(in)                         :: nbdlphys
!
!*** local variables
   integer i, j, k, n, rc, idx, ibdl, nbdl
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
   if(mpp_pe()==mpp_root_pe())print *,'in fv_phys bundle,nbdl=',nbdlphys
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
!       if(mpp_pe()==mpp_root_pe())print *,'in fv3gfsio, vertical
!       list=',udimList(1:udimCount),'rc=',rc

       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

     else

       if(mpp_pe()==mpp_root_pe())print *,'in fv_dyn bundle,axis_name_vert=',axis_name_vert
       call ESMF_AttributeAdd(fcst_grid, convention="NetCDF", purpose="FV3",  &
                              attrList=(/"vertical_dim_labels"/), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
                              name="vertical_dim_labels", valueList=axis_name_vert, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
     endif
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
       edgesS=''
       do i = 1,num_axes_phys
         if(axes(i) == edges) edgesS=axis_name(i)
       enddo
! Add vertical dimension Attributes to Grid
       if( id>2 ) then
!      if(mpp_pe()==mpp_root_pe())print *,' in dyn add grid, axis_name=',     &
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
        name='missing_value',value=missing_value,rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
        attrList=(/"_FillValue"/), rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
        name='_FillValue',value=missing_value,rc=rc)
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
