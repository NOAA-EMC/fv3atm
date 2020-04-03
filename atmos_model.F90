!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of fvGFS.                                       *
!*                                                                     *
!* fvGFS is free software; you can redistribute it and/or modify it    *
!* and are expected to follow the terms of the GNU General Public      *
!* License as published by the Free Software Foundation; either        *
!* version 2 of the License, or (at your option) any later version.    *
!*                                                                     *
!* fvGFS is distributed in the hope that it will be useful, but        *
!* WITHOUT ANY WARRANTY; without even the implied warranty of          *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   *
!* General Public License for more details.                            *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************
module atmos_model_mod
!-----------------------------------------------------------------------
!<OVERVIEW>
!  Driver for the atmospheric model, contains routines to advance the
!  atmospheric model state by one time step.
!</OVERVIEW>

!<DESCRIPTION>
!     This version of atmos_model_mod has been designed around the implicit
!     version diffusion scheme of the GCM. It requires two routines to advance
!     the atmospheric model one time step into the future. These two routines
!     correspond to the down and up sweeps of the standard tridiagonal solver.
!     Most atmospheric processes (dynamics,radiation,etc.) are performed
!     in the down routine. The up routine finishes the vertical diffusion
!     and computes moisture related terms (convection,large-scale condensation,
!     and precipitation).

!     The boundary variables needed by other component models for coupling
!     are contained in a derived data type. A variable of this derived type
!     is returned when initializing the atmospheric model. It is used by other
!     routines in this module and by coupling routines. The contents of
!     this derived type should only be modified by the atmospheric model.

!</DESCRIPTION>

use mpp_mod,            only: mpp_pe, mpp_root_pe, mpp_clock_id, mpp_clock_begin
use mpp_mod,            only: mpp_clock_end, CLOCK_COMPONENT, MPP_CLOCK_SYNC
use mpp_mod,            only: FATAL, mpp_min, mpp_max, mpp_error, mpp_chksum
use mpp_domains_mod,    only: domain2d
use mpp_mod,            only: mpp_get_current_pelist_name
#ifdef INTERNAL_FILE_NML
use mpp_mod,            only: input_nml_file
#else
use fms_mod,            only: open_namelist_file
#endif
use fms_mod,            only: file_exist, error_mesg
use fms_mod,            only: close_file, write_version_number, stdlog, stdout
use fms_mod,            only: clock_flag_default
use fms_mod,            only: check_nml_error
use diag_manager_mod,   only: diag_send_complete_instant
use time_manager_mod,   only: time_type, get_time, get_date, &
                              operator(+), operator(-),real_to_time_type
use field_manager_mod,  only: MODEL_ATMOS
use tracer_manager_mod, only: get_number_tracers, get_tracer_names, &
                              get_tracer_index, NO_TRACER
use xgrid_mod,          only: grid_box_type
use atmosphere_mod,     only: atmosphere_init
use atmosphere_mod,     only: atmosphere_restart
use atmosphere_mod,     only: atmosphere_end
use atmosphere_mod,     only: atmosphere_state_update
use atmosphere_mod,     only: atmos_phys_driver_statein
use atmosphere_mod,     only: atmosphere_control_data
use atmosphere_mod,     only: atmosphere_resolution, atmosphere_domain
use atmosphere_mod,     only: atmosphere_grid_bdry, atmosphere_grid_ctr
use atmosphere_mod,     only: atmosphere_dynamics, atmosphere_diag_axes
use atmosphere_mod,     only: atmosphere_etalvls, atmosphere_hgt
!rab use atmosphere_mod,     only: atmosphere_tracer_postinit
use atmosphere_mod,     only: atmosphere_diss_est, atmosphere_nggps_diag
use atmosphere_mod,     only: atmosphere_scalar_field_halo
use atmosphere_mod,     only: atmosphere_get_bottom_layer
use atmosphere_mod,     only: set_atmosphere_pelist
use atmosphere_mod,     only: Atm, mytile
use block_control_mod,  only: block_control_type, define_blocks_packed
use DYCORE_typedefs,    only: DYCORE_data_type, DYCORE_diag_type
#ifdef CCPP
use IPD_typedefs,       only: IPD_init_type, IPD_diag_type,    &
                              IPD_restart_type, IPD_kind_phys, &
                              IPD_func0d_proc, IPD_func1d_proc
#else
use IPD_typedefs,       only: IPD_init_type, IPD_control_type, &
                              IPD_data_type, IPD_diag_type,    &
                              IPD_restart_type, IPD_kind_phys, &
                              IPD_func0d_proc, IPD_func1d_proc
#endif

#ifdef CCPP
use CCPP_data,          only: ccpp_suite,                      &
                              IPD_control => GFS_control,      &
                              IPD_data => GFS_data,            &
                              IPD_interstitial => GFS_interstitial
use IPD_driver,         only: IPD_initialize, IPD_initialize_rst
use CCPP_driver,        only: CCPP_step, non_uniform_blocks
#else
use IPD_driver,         only: IPD_initialize, IPD_initialize_rst, IPD_step
use physics_abstraction_layer, only: time_vary_step, radiation_step1, physics_step1, physics_step2
#endif

use stochastic_physics, only: init_stochastic_physics,         &
                              run_stochastic_physics
use stochastic_physics_sfc, only: run_stochastic_physics_sfc

use FV3GFS_io_mod,      only: FV3GFS_restart_read, FV3GFS_restart_write, &
                              FV3GFS_IPD_checksum,                       &
                              FV3GFS_diag_register, FV3GFS_diag_output,  &
                              DIAG_SIZE
use fv_iau_mod,         only: iau_external_data_type,getiauforcing,iau_initialize
use module_fv3_config,  only: output_1st_tstep_rst, first_kdt, nsout,    &
                              frestart, restart_endfcst

!-----------------------------------------------------------------------

implicit none
private

public update_atmos_radiation_physics
public update_atmos_model_state
public update_atmos_model_dynamics
public atmos_model_init, atmos_model_end, atmos_data_type
public atmos_model_exchange_phase_1, atmos_model_exchange_phase_2
public atmos_model_restart
public get_atmos_model_ungridded_dim
public addLsmask2grid
!-----------------------------------------------------------------------

!<PUBLICTYPE >
 type atmos_data_type
     integer                       :: axes(4)            ! axis indices (returned by diag_manager) for the atmospheric grid 
                                                         ! (they correspond to the x, y, pfull, phalf axes)
     integer, pointer              :: pelist(:) =>null() ! pelist where atmosphere is running.
     integer                       :: layout(2)          ! computer task laytout
     logical                       :: regional           ! true if domain is regional
     logical                       :: nested             ! true if there is a nest
     integer                       :: mlon, mlat
     integer                       :: iau_offset         ! iau running window length
     logical                       :: pe                 ! current pe.
     real(kind=8),             pointer, dimension(:)     :: ak, bk
     real,                     pointer, dimension(:,:)   :: lon_bnd  => null() ! local longitude axis grid box corners in radians.
     real,                     pointer, dimension(:,:)   :: lat_bnd  => null() ! local latitude axis grid box corners in radians.
     real(kind=IPD_kind_phys), pointer, dimension(:,:)   :: lon      => null() ! local longitude axis grid box centers in radians.
     real(kind=IPD_kind_phys), pointer, dimension(:,:)   :: lat      => null() ! local latitude axis grid box centers in radians.
     real(kind=IPD_kind_phys), pointer, dimension(:,:)   :: dx, dy
     real(kind=8),             pointer, dimension(:,:)   :: area
     real(kind=8),             pointer, dimension(:,:,:) :: layer_hgt, level_hgt
     type(domain2d)                :: domain             ! domain decomposition
     type(time_type)               :: Time               ! current time
     type(time_type)               :: Time_step          ! atmospheric time step.
     type(time_type)               :: Time_init          ! reference time.
     type(grid_box_type)           :: grid               ! hold grid information needed for 2nd order conservative flux exchange 
     type(IPD_diag_type), pointer, dimension(:) :: Diag
 end type atmos_data_type
                                                         ! to calculate gradient on cubic sphere grid.
!</PUBLICTYPE >

integer :: fv3Clock, getClock, updClock, setupClock, radClock, physClock

!-----------------------------------------------------------------------
integer :: blocksize    = 1
logical :: chksum_debug = .false.
logical :: dycore_only  = .false.
logical :: debug        = .false.
!logical :: debug        = .true.
logical :: sync         = .false.
integer, parameter     :: maxhr = 4096
real, dimension(maxhr) :: fdiag = 0.
real                   :: fhmax=384.0, fhmaxhf=120.0, fhout=3.0, fhouthf=1.0,avg_max_length=3600.
#ifdef CCPP
namelist /atmos_model_nml/ blocksize, chksum_debug, dycore_only, debug, sync, fdiag, fhmax, fhmaxhf, fhout, fhouthf, ccpp_suite, avg_max_length
#else
namelist /atmos_model_nml/ blocksize, chksum_debug, dycore_only, debug, sync, fdiag, fhmax, fhmaxhf, fhout, fhouthf, avg_max_length
#endif

type (time_type) :: diag_time, diag_time_fhzero

!--- concurrent and decoupled radiation and physics variables
!-------------------
!  DYCORE containers
!-------------------
type(DYCORE_data_type),    allocatable :: DYCORE_Data(:)  ! number of blocks
type(DYCORE_diag_type)                 :: DYCORE_Diag(25)

!----------------
!  IPD containers
!----------------
#ifndef CCPP
type(IPD_control_type)              :: IPD_Control
type(IPD_data_type),    allocatable :: IPD_Data(:)  ! number of blocks
type(IPD_diag_type),    target      :: IPD_Diag(DIAG_SIZE)
type(IPD_restart_type)              :: IPD_Restart
#else
! IPD_Control and IPD_Data are coming from CCPP_data
type(IPD_diag_type),    target      :: IPD_Diag(DIAG_SIZE)
type(IPD_restart_type)              :: IPD_Restart
#endif

!--------------
! IAU container
!--------------
type(iau_external_data_type)        :: IAU_Data ! number of blocks

!-----------------
!  Block container
!-----------------
type (block_control_type), target   :: Atm_block

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

#ifdef NAM_phys
  logical,parameter :: flip_vc = .false.
#else
  logical,parameter :: flip_vc = .true.
#endif

  real(kind=IPD_kind_phys), parameter :: zero = 0.0_IPD_kind_phys, &
                                         one  = 1.0_IPD_kind_phys

contains

!#######################################################################
! <SUBROUTINE NAME="update_radiation_physics">
!
!<DESCRIPTION>
!   Called every time step as the atmospheric driver to compute the
!   atmospheric tendencies for dynamics, radiation, vertical diffusion of
!   momentum, tracers, and heat/moisture.  For heat/moisture only the
!   downward sweep of the tridiagonal elimination is performed, hence
!   the name "_down". 
!</DESCRIPTION>

!   <TEMPLATE>
!     call  update_atmos_radiation_physics (Atmos)
!   </TEMPLATE>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </INOUT>

subroutine update_atmos_radiation_physics (Atmos)
#ifdef OPENMP
    use omp_lib
#endif
!-----------------------------------------------------------------------
  type (atmos_data_type), intent(in) :: Atmos
!--- local variables---
    integer :: nb, jdat(8), rc
    procedure(IPD_func0d_proc), pointer :: Func0d => NULL()
    procedure(IPD_func1d_proc), pointer :: Func1d => NULL()
    integer :: nthrds
#ifdef CCPP
    integer :: ierr
#endif

#ifdef OPENMP
    nthrds = omp_get_max_threads()
#else
    nthrds = 1
#endif

    if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "statein driver"
!--- get atmospheric state from the dynamic core
    call set_atmosphere_pelist()
    call mpp_clock_begin(getClock)
    if (IPD_control%do_skeb) call atmosphere_diss_est (IPD_control%skeb_npass) !  do smoothing for SKEB
    call atmos_phys_driver_statein (IPD_data, Atm_block, flip_vc)
    call mpp_clock_end(getClock)

!--- if dycore only run, set up the dummy physics output state as the input state
    if (dycore_only) then
      do nb = 1,Atm_block%nblks
        IPD_Data(nb)%Stateout%gu0 = IPD_Data(nb)%Statein%ugrs
        IPD_Data(nb)%Stateout%gv0 = IPD_Data(nb)%Statein%vgrs
        IPD_Data(nb)%Stateout%gt0 = IPD_Data(nb)%Statein%tgrs
        IPD_Data(nb)%Stateout%gq0 = IPD_Data(nb)%Statein%qgrs
      enddo
    else
      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "setup step"

!--- update IPD_Control%jdat(8)
      jdat(:) = 0
      call get_date (Atmos%Time, jdat(1), jdat(2), jdat(3),  &
                                 jdat(5), jdat(6), jdat(7))
      IPD_Control%jdat(:) = jdat(:)

!--- execute the IPD atmospheric setup step
      call mpp_clock_begin(setupClock)
#ifdef CCPP
      call CCPP_step (step="time_vary", nblks=Atm_block%nblks, ierr=ierr)
      if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP time_vary step failed')
#else
      Func1d => time_vary_step
      call IPD_step (IPD_Control, IPD_Data(:), IPD_Diag, IPD_Restart, IPD_func1d=Func1d)
#endif

!--- call stochastic physics pattern generation / cellular automata
    if (IPD_Control%do_sppt .OR. IPD_Control%do_shum .OR. IPD_Control%do_skeb .OR. IPD_Control%do_sfcperts) then
       call run_stochastic_physics(IPD_Control, IPD_Data(:)%Grid, IPD_Data(:)%Coupling, nthrds)
    end if

    if(IPD_Control%do_ca)then
       ! DH* The current implementation of cellular_automata assumes that all blocksizes are the
       ! same, this is tested in the initialization call to cellular_automata, no need to redo *DH
       call cellular_automata(IPD_Control%kdt, IPD_Data(:)%Statein, IPD_Data(:)%Coupling, IPD_Data(:)%Intdiag, &
                              Atm_block%nblks, IPD_Control%levs, IPD_Control%nca, IPD_Control%ncells,          &
                              IPD_Control%nlives, IPD_Control%nfracseed, IPD_Control%nseed,                    &
                              IPD_Control%nthresh, IPD_Control%ca_global, IPD_Control%ca_sgs,                  &
                              IPD_Control%iseed_ca, IPD_Control%ca_smooth, IPD_Control%nspinup,                &
                              Atm_block%blksz(1))
    endif

!--- if coupled, assign coupled fields
      if( IPD_Control%cplflx .or. IPD_Control%cplwav ) then
!        print *,'in atmos_model,nblks=',Atm_block%nblks
!        print *,'in atmos_model,IPD_Data size=',size(IPD_Data)
!        print *,'in atmos_model,tsfc(1)=',IPD_Data(1)%sfcprop%tsfc(1)
!        print *,'in atmos_model, tsfc size=',size(IPD_Data(1)%sfcprop%tsfc)
        call assign_importdata(rc)
!        print *,'in atmos_model, after assign_importdata, rc=',rc
      endif

      call mpp_clock_end(setupClock)

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "radiation driver"

!--- execute the IPD atmospheric radiation subcomponent (RRTM)

      call mpp_clock_begin(radClock)
#ifdef CCPP
      ! Performance improvement. Only enter if it is time to call the radiation physics.
      if (IPD_Control%lsswr .or. IPD_Control%lslwr) then
        call CCPP_step (step="radiation", nblks=Atm_block%nblks, ierr=ierr)
        if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP radiation step failed')
      endif
#else
      Func0d => radiation_step1
!$OMP parallel do default (none)       &
!$OMP            schedule (dynamic,1), &
!$OMP            shared   (Atm_block, IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, Func0d) &
!$OMP            private  (nb)
      do nb = 1,Atm_block%nblks
        call IPD_step (IPD_Control, IPD_Data(nb:nb), IPD_Diag, IPD_Restart, IPD_func0d=Func0d)
      enddo
#endif
      call mpp_clock_end(radClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'RADIATION STEP  ', IPD_Control%kdt, IPD_Control%fhour
        call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
      endif

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "physics driver"

!--- execute the IPD atmospheric physics step1 subcomponent (main physics driver)

      call mpp_clock_begin(physClock)
#ifdef CCPP
      call CCPP_step (step="physics", nblks=Atm_block%nblks, ierr=ierr)
      if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP physics step failed')
#else
      Func0d => physics_step1
!$OMP parallel do default (none) &
!$OMP            schedule (dynamic,1), &
!$OMP            shared   (Atm_block, IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, Func0d) &
!$OMP            private  (nb)
      do nb = 1,Atm_block%nblks
        call IPD_step (IPD_Control, IPD_Data(nb:nb), IPD_Diag, IPD_Restart, IPD_func0d=Func0d)
      enddo
#endif
      call mpp_clock_end(physClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'PHYSICS STEP1   ', IPD_Control%kdt, IPD_Control%fhour
        call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
      endif

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "stochastic physics driver"

!--- execute the IPD atmospheric physics step2 subcomponent (stochastic physics driver)

      call mpp_clock_begin(physClock)
#ifdef CCPP
      call CCPP_step (step="stochastics", nblks=Atm_block%nblks, ierr=ierr)
      if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP stochastics step failed')
#else
      Func0d => physics_step2
!$OMP parallel do default (none) &
!$OMP            schedule (dynamic,1), &
!$OMP            shared   (Atm_block, IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, Func0d) &
!$OMP            private  (nb)
      do nb = 1,Atm_block%nblks
        call IPD_step (IPD_Control, IPD_Data(nb:nb), IPD_Diag, IPD_Restart, IPD_func0d=Func0d)
      enddo
#endif
      call mpp_clock_end(physClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'PHYSICS STEP2   ', IPD_Control%kdt, IPD_Control%fhour
        call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
      endif
      call getiauforcing(IPD_Control,IAU_data)
      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "end of radiation and physics step"
    endif

#ifdef CCPP
    ! Update flag for first time step of time integration
    IPD_Control%first_time_step = .false.
#endif
!-----------------------------------------------------------------------
 end subroutine update_atmos_radiation_physics
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="atmos_model_init">
!
! <OVERVIEW>
! Routine to initialize the atmospheric model
! </OVERVIEW>

subroutine atmos_model_init (Atmos, Time_init, Time, Time_step)

#ifdef OPENMP
  use omp_lib
#endif
#ifdef CCPP
  use fv_mp_mod, only: commglobal
#endif
  use mpp_mod, only: mpp_npes

  type (atmos_data_type), intent(inout) :: Atmos
  type (time_type), intent(in) :: Time_init, Time, Time_step
!--- local variables ---
  integer :: unit, ntdiag, ntfamily, i, j, k
  integer :: mlon, mlat, nlon, nlat, nlev, sec, dt
  integer :: ierr, io, logunit
  integer :: idx, tile_num
  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: blk, ibs, ibe, jbs, jbe
  real(kind=IPD_kind_phys) :: dt_phys
  real, allocatable    :: q(:,:,:,:), p_half(:,:,:)
  character(len=80)    :: control
  character(len=64)    :: filename, filename2, pelist_name
  character(len=132)   :: text
  logical              :: p_hydro, hydro, fexist
  logical, save        :: block_message = .true.
  type(IPD_init_type)  :: Init_parm
  integer              :: bdat(8), cdat(8)
  integer              :: ntracers, maxhf, maxh
  character(len=32), allocatable, target :: tracer_names(:)
  integer :: nthrds

!-----------------------------------------------------------------------

!---- set the atmospheric model time ------

   Atmos % Time_init = Time_init
   Atmos % Time      = Time
   Atmos % Time_step = Time_step
   call get_time (Atmos % Time_step, sec)
   dt_phys = real(sec)      ! integer seconds

   logunit = stdlog()

!-----------------------------------------------------------------------
! initialize atmospheric model -----

#ifndef CCPP
!---------- initialize atmospheric dynamics -------
   call atmosphere_init (Atmos%Time_init, Atmos%Time, Atmos%Time_step,&
                         Atmos%grid, Atmos%area)
#endif

   IF ( file_exist('input.nml')) THEN
#ifdef INTERNAL_FILE_NML
      read(input_nml_file, nml=atmos_model_nml, iostat=io)
      ierr = check_nml_error(io, 'atmos_model_nml')
#else
      unit = open_namelist_file ( )
      ierr=1
      do while (ierr /= 0)
         read  (unit, nml=atmos_model_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'atmos_model_nml')
      enddo
 10     call close_file (unit)
#endif
   endif

#ifdef CCPP
!---------- initialize atmospheric dynamics after reading the namelist -------
!---------- (need name of CCPP suite definition file from input.nml) ---------
   call atmosphere_init (Atmos%Time_init, Atmos%Time, Atmos%Time_step,&
                         Atmos%grid, Atmos%area)
#endif

!-----------------------------------------------------------------------
   call atmosphere_resolution (nlon, nlat, global=.false.)
   call atmosphere_resolution (mlon, mlat, global=.true.)
   call alloc_atmos_data_type (nlon, nlat, Atmos)
   call atmosphere_domain (Atmos%domain, Atmos%layout, Atmos%regional, Atmos%nested, Atmos%pelist)
   call atmosphere_diag_axes (Atmos%axes)
   call atmosphere_etalvls (Atmos%ak, Atmos%bk, flip=flip_vc)
   call atmosphere_grid_bdry (Atmos%lon_bnd, Atmos%lat_bnd, global=.false.)
   call atmosphere_grid_ctr (Atmos%lon, Atmos%lat)
   call atmosphere_hgt (Atmos%layer_hgt, 'layer', relative=.false., flip=flip_vc)
   call atmosphere_hgt (Atmos%level_hgt, 'level', relative=.false., flip=flip_vc)

   Atmos%mlon = mlon
   Atmos%mlat = mlat
!-----------------------------------------------------------------------
!--- before going any further check definitions for 'blocks'
!-----------------------------------------------------------------------
   call atmosphere_control_data (isc, iec, jsc, jec, nlev, p_hydro, hydro, tile_num)
   call define_blocks_packed ('atmos_model', Atm_block, isc, iec, jsc, jec, nlev, &
                              blocksize, block_message)
   
   allocate(DYCORE_Data(Atm_block%nblks))
   allocate(IPD_Data(Atm_block%nblks))

#ifdef OPENMP
   nthrds = omp_get_max_threads()
#else
   nthrds = 1
#endif

#ifdef CCPP
   ! This logic deals with non-uniform block sizes for CCPP.
   ! When non-uniform block sizes are used, it is required
   ! that only the last block has a different (smaller)
   ! size than all other blocks. This is the standard in
   ! FV3. If this is the case, set non_uniform_blocks (a
   ! variable imported from CCPP_driver) to .true. and
   ! allocate nthreads+1 elements of the interstitial array.
   ! The extra element will be used by the thread that
   ! runs over the last, smaller block.
   if (minval(Atm_block%blksz)==maxval(Atm_block%blksz)) then
      non_uniform_blocks = .false.
      allocate(IPD_Interstitial(nthrds))
   else if (all(minloc(Atm_block%blksz)==(/size(Atm_block%blksz)/))) then
      non_uniform_blocks = .true.
      allocate(IPD_Interstitial(nthrds+1))
   else
      call mpp_error(FATAL, 'For non-uniform blocksizes, only the last element ' // &
                            'in Atm_block%blksz can be different from the others')
   end if

#endif

!--- update IPD_Control%jdat(8)
   bdat(:) = 0
   call get_date (Time_init, bdat(1), bdat(2), bdat(3),  &
                             bdat(5), bdat(6), bdat(7))
   cdat(:) = 0
   call get_date (Time,      cdat(1), cdat(2), cdat(3),  &
                             cdat(5), cdat(6), cdat(7))
   call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
   allocate (tracer_names(ntracers))
   do i = 1, ntracers
     call get_tracer_names(MODEL_ATMOS, i, tracer_names(i))
   enddo
!--- setup IPD Init_parm
   Init_parm%me              =  mpp_pe()
   Init_parm%master          =  mpp_root_pe()
   Init_parm%tile_num        =  tile_num
   Init_parm%isc             =  isc
   Init_parm%jsc             =  jsc
   Init_parm%nx              =  nlon
   Init_parm%ny              =  nlat
   Init_parm%levs            =  nlev
   Init_parm%cnx             =  mlon
   Init_parm%cny             =  mlat
   Init_parm%gnx             =  Init_parm%cnx*4
   Init_parm%gny             =  Init_parm%cny*2
   Init_parm%nlunit          =  9999
   Init_parm%logunit         =  logunit
   Init_parm%bdat(:)         =  bdat(:)
   Init_parm%cdat(:)         =  cdat(:)
   Init_parm%dt_dycore       =  dt_phys
   Init_parm%dt_phys         =  dt_phys
   Init_parm%iau_offset      =  Atmos%iau_offset
   Init_parm%blksz           => Atm_block%blksz
   Init_parm%ak              => Atmos%ak
   Init_parm%bk              => Atmos%bk
   Init_parm%xlon            => Atmos%lon
   Init_parm%xlat            => Atmos%lat
   Init_parm%area            => Atmos%area
   Init_parm%tracer_names    => tracer_names
#ifdef CCPP
   Init_parm%restart         = Atm(mytile)%flagstruct%warm_start
   Init_parm%hydrostatic     = Atm(mytile)%flagstruct%hydrostatic
#endif

#ifdef INTERNAL_FILE_NML
   Init_parm%input_nml_file  => input_nml_file
   Init_parm%fn_nml='using internal file'
#else
   pelist_name=mpp_get_current_pelist_name()
   Init_parm%fn_nml='input_'//trim(pelist_name)//'.nml'
   inquire(FILE=Init_parm%fn_nml, EXIST=fexist)
   if (.not. fexist ) then
      Init_parm%fn_nml='input.nml'
   endif
#endif

#ifdef CCPP
   call IPD_initialize (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, &
                        IPD_Interstitial, commglobal, mpp_npes(), Init_parm)
#else
   call IPD_initialize (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, Init_parm)
#endif

   if (IPD_Control%do_sppt .OR. IPD_Control%do_shum .OR. IPD_Control%do_skeb .OR. IPD_Control%do_sfcperts) then
      ! Initialize stochastic physics
      call init_stochastic_physics(IPD_Control, Init_parm, mpp_npes(), nthrds)
      if(IPD_Control%me == IPD_Control%master) print *,'do_skeb=',IPD_Control%do_skeb
   end if

   Atmos%Diag => IPD_Diag

   if (IPD_Control%do_sfcperts) then
      ! Get land surface perturbations here (move to GFS_time_vary
      ! step if wanting to update each time-step)
      call run_stochastic_physics_sfc(IPD_Control, IPD_Data(:)%Grid, IPD_Data(:)%Coupling)
   end if

   ! Initialize cellular automata
   if(IPD_Control%do_ca)then
      ! DH* The current implementation of cellular_automata assumes that all blocksizes are the
      ! same - abort if this is not the case, otherwise proceed with Atm_block%blksz(1) below
      if (.not. minval(Atm_block%blksz)==maxval(Atm_block%blksz)) then
         call mpp_error(FATAL, 'Logic errror: cellular_automata not compatible with non-uniform blocksizes')
      end if
      ! *DH
      call cellular_automata(IPD_Control%kdt, IPD_Data(:)%Statein, IPD_Data(:)%Coupling, IPD_Data(:)%Intdiag, &
                             Atm_block%nblks, IPD_Control%levs, IPD_Control%nca, IPD_Control%ncells,          &
                             IPD_Control%nlives, IPD_Control%nfracseed, IPD_Control%nseed,                    &
                             IPD_Control%nthresh, IPD_Control%ca_global, IPD_Control%ca_sgs,                  &
                             IPD_Control%iseed_ca, IPD_Control%ca_smooth, IPD_Control%nspinup,                &
                             Atm_block%blksz(1))
   endif

   Atm(mytile)%flagstruct%do_skeb = IPD_Control%do_skeb

!  initialize the IAU module
   call iau_initialize (IPD_Control,IAU_data,Init_parm)

   Init_parm%blksz           => null()
   Init_parm%ak              => null()
   Init_parm%bk              => null()
   Init_parm%xlon            => null()
   Init_parm%xlat            => null()
   Init_parm%area            => null()
   Init_parm%tracer_names    => null()
   deallocate (tracer_names)

   !--- update tracers in FV3 with any initialized during the physics/radiation init phase
!rab   call atmosphere_tracer_postinit (IPD_Data, Atm_block)

   call atmosphere_nggps_diag (Time, init=.true.)
   call FV3GFS_diag_register (IPD_Diag, Time, Atm_block, IPD_Control, Atmos%lon, Atmos%lat, Atmos%axes)
   call IPD_initialize_rst (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, Init_parm)
#ifdef CCPP
   call FV3GFS_restart_read (IPD_Data, IPD_Restart, Atm_block, IPD_Control, Atmos%domain, Atm(mytile)%flagstruct%warm_start)
#else
   call FV3GFS_restart_read (IPD_Data, IPD_Restart, Atm_block, IPD_Control, Atmos%domain)
#endif

#ifdef CCPP
   ! Populate the IPD_Data%Statein container with the prognostic state
   ! in Atm_block, which contains the initial conditions/restart data.
   call atmos_phys_driver_statein (IPD_data, Atm_block, flip_vc)
   ! Initialize the CCPP framework
   call CCPP_step (step="init", nblks=Atm_block%nblks, ierr=ierr)
   if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP init step failed')
   ! Initialize the CCPP physics
   call CCPP_step (step="physics_init", nblks=Atm_block%nblks, ierr=ierr)
   if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP physics_init step failed')
#endif

   !--- set the initial diagnostic timestamp
   diag_time = Time 
   if (output_1st_tstep_rst) then
     diag_time = Time - real_to_time_type(mod(int((first_kdt - 1)*dt_phys/3600.),6)*3600.0)
   endif
   if (Atmos%iau_offset > zero) then
     diag_time = Atmos%Time_init
     diag_time_fhzero = Atmos%Time
   endif

   !---- print version number to logfile ----

   call write_version_number ( version, tagname )
   !--- write the namelist to a log file
   if (mpp_pe() == mpp_root_pe()) then
      unit = stdlog( )
      write (unit, nml=atmos_model_nml)
      call close_file (unit)
   endif

   !--- get fdiag
#ifdef GFS_PHYS
!--- check fdiag to see if it is an interval or a list
   if (nint(fdiag(2)) == 0) then
     if (fhmaxhf > 0) then
       maxhf = fhmaxhf / fhouthf
       maxh  = maxhf + (fhmax-fhmaxhf) / fhout
       fdiag(1) = fhouthf
       do i=2,maxhf
        fdiag(i) = fdiag(i-1) + fhouthf
       enddo
       do i=maxhf+1,maxh
         fdiag(i) = fdiag(i-1) + fhout
       enddo
     else
       maxh  = fhmax / fhout
       do i = 2, maxh
         fdiag(i) = fdiag(i-1) + fhout
       enddo
     endif
   endif
   if (mpp_pe() == mpp_root_pe()) write(6,*) "---fdiag",fdiag(1:40)
#endif

   setupClock = mpp_clock_id( 'GFS Step Setup        ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   radClock   = mpp_clock_id( 'GFS Radiation         ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   physClock  = mpp_clock_id( 'GFS Physics           ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   getClock   = mpp_clock_id( 'Dynamics get state    ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   updClock   = mpp_clock_id( 'Dynamics update state ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   if (sync) then
     fv3Clock = mpp_clock_id( 'FV3 Dycore            ', flags=clock_flag_default+MPP_CLOCK_SYNC, grain=CLOCK_COMPONENT )
   else
     fv3Clock = mpp_clock_id( 'FV3 Dycore            ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   endif

#ifdef CCPP
   ! Set flag for first time step of time integration
   IPD_Control%first_time_step = .true.
#endif
!-----------------------------------------------------------------------
end subroutine atmos_model_init
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_dynamics"
!
! <OVERVIEW>
subroutine update_atmos_model_dynamics (Atmos)
! run the atmospheric dynamics to advect the properties
  type (atmos_data_type), intent(in) :: Atmos

    call set_atmosphere_pelist()
    call mpp_clock_begin(fv3Clock)
    call atmosphere_dynamics (Atmos%Time)
    call mpp_clock_end(fv3Clock)

end subroutine update_atmos_model_dynamics
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="atmos_model_exchange_phase_1"
!
! <OVERVIEW>
!   Perform data exchange with coupled components in run phase 1
! </OVERVIEW>
!
! <DESCRIPTION>
!  This subroutine currently exports atmospheric fields and tracers
!  to the chemistry component during the model's run phase 1, i.e.
!  before chemistry is run.
! </DESCRIPTION>

subroutine atmos_model_exchange_phase_1 (Atmos, rc)

  use ESMF

  type (atmos_data_type), intent(inout) :: Atmos
  integer, optional,      intent(out)   :: rc
!--- local variables
  integer :: localrc

    !--- begin
    if (present(rc)) rc = ESMF_SUCCESS

    !--- if coupled, exchange coupled fields
    if( IPD_Control%cplchm ) then
      ! -- export fields to chemistry
      call update_atmos_chemistry('export', rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
    endif

 end subroutine atmos_model_exchange_phase_1
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="atmos_model_exchange_phase_2"
!
! <OVERVIEW>
!   Perform data exchange with coupled components in run phase 2
! </OVERVIEW>
!
! <DESCRIPTION>
!  This subroutine currently imports fields updated by the coupled
!  chemistry component back into the atmospheric model during run
!  phase 2.
! </DESCRIPTION>

subroutine atmos_model_exchange_phase_2 (Atmos, rc)

  use ESMF

  type (atmos_data_type), intent(inout) :: Atmos
  integer, optional,      intent(out)   :: rc
!--- local variables
  integer :: localrc

    !--- begin
    if (present(rc)) rc = ESMF_SUCCESS

    !--- if coupled, exchange coupled fields
    if( IPD_Control%cplchm ) then
      ! -- import fields from chemistry
      call update_atmos_chemistry('import', rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
    endif

 end subroutine atmos_model_exchange_phase_2
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_state"
!
! <OVERVIEW>
subroutine update_atmos_model_state (Atmos)
! to update the model state after all concurrency is completed
  type (atmos_data_type), intent(inout) :: Atmos
!--- local variables
  integer :: isec, seconds, isec_fhzero
  integer :: rc
  real(kind=IPD_kind_phys) :: time_int, time_intfull
!
    call set_atmosphere_pelist()
    call mpp_clock_begin(fv3Clock)
    call mpp_clock_begin(updClock)
    call atmosphere_state_update (Atmos%Time, IPD_Data, IAU_Data, Atm_block, flip_vc)
    call mpp_clock_end(updClock)
    call mpp_clock_end(fv3Clock)

    if (chksum_debug) then
      if (mpp_pe() == mpp_root_pe()) print *,'UPDATE STATE    ', IPD_Control%kdt, IPD_Control%fhour
      if (mpp_pe() == mpp_root_pe()) print *,'in UPDATE STATE    ', size(IPD_Data(1)%SfcProp%tsfc),'nblks=',Atm_block%nblks
      call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
    endif

    !--- advance time ---
    Atmos % Time = Atmos % Time + Atmos % Time_step

    call get_time (Atmos%Time - diag_time, isec)
    call get_time (Atmos%Time - Atmos%Time_init, seconds)
    call atmosphere_nggps_diag(Atmos%Time,ltavg=.true.,avg_max_length=avg_max_length)
    if (ANY(nint(fdiag(:)*3600.0) == seconds) .or. (IPD_Control%kdt == first_kdt) .or. nsout > 0) then
      if (mpp_pe() == mpp_root_pe()) write(6,*) "---isec,seconds",isec,seconds
      time_int = real(isec)
      if(Atmos%iau_offset > zero) then
        if( time_int - Atmos%iau_offset*3600. > zero ) then
          time_int = time_int - Atmos%iau_offset*3600.
        else if(seconds == Atmos%iau_offset*3600) then
          call get_time (Atmos%Time - diag_time_fhzero, isec_fhzero)
          time_int = real(isec_fhzero)
          if (mpp_pe() == mpp_root_pe()) write(6,*) "---iseczero",isec_fhzero
        endif
      endif
      time_intfull = real(seconds)
      if(Atmos%iau_offset > zero) then
        if( time_intfull - Atmos%iau_offset*3600. > zero) then
          time_intfull = time_intfull - Atmos%iau_offset*3600.
        endif
      endif
      if (mpp_pe() == mpp_root_pe()) write(6,*) ' gfs diags time since last bucket empty: ',time_int/3600.,'hrs'
      call atmosphere_nggps_diag(Atmos%Time)
      call FV3GFS_diag_output(Atmos%Time, IPD_DIag, Atm_block, IPD_Control%nx, IPD_Control%ny, &
                            IPD_Control%levs, 1, 1, 1.d0, time_int, time_intfull,              &
                            IPD_Control%fhswr, IPD_Control%fhlwr)
      if (nint(IPD_Control%fhzero) > 0) then 
        if (mod(isec,3600*nint(IPD_Control%fhzero)) == 0) diag_time = Atmos%Time
      else
        if (mod(isec,nint(3600*IPD_Control%fhzero)) == 0) diag_time = Atmos%Time
      endif
      call diag_send_complete_instant (Atmos%Time)
    endif

    !--- this may not be necessary once write_component is fully implemented
    !!!call diag_send_complete_extra (Atmos%Time)

    !--- get bottom layer data from dynamical core for coupling
    call atmosphere_get_bottom_layer (Atm_block, DYCORE_Data) 

    !if in coupled mode, set up coupled fields
    if (IPD_Control%cplflx .or. IPD_Control%cplwav) then
!     if (mpp_pe() == mpp_root_pe()) print *,'COUPLING: IPD layer'
!jw       call setup_exportdata(IPD_Control, IPD_Data, Atm_block)
      call setup_exportdata(rc)
    endif

 end subroutine update_atmos_model_state
! </SUBROUTINE>



!#######################################################################
! <SUBROUTINE NAME="atmos_model_end">
!
! <OVERVIEW>
!  termination routine for atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!  Call once to terminate this module and any other modules used.
!  This routine writes a restart file and deallocates storage
!  used by the derived-type variable atmos_boundary_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_model_end (Atmos)
! </TEMPLATE>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </INOUT>

subroutine atmos_model_end (Atmos)
  type (atmos_data_type), intent(inout) :: Atmos
!---local variables
  integer :: idx, seconds
#ifdef CCPP
  integer :: ierr
#endif

!-----------------------------------------------------------------------
!---- termination routine for atmospheric model ----
                                              
    call atmosphere_end (Atmos % Time, Atmos%grid, restart_endfcst)
    if(restart_endfcst) then
      call FV3GFS_restart_write (IPD_Data, IPD_Restart, Atm_block, &
                                 IPD_Control, Atmos%domain)
    endif

#ifdef CCPP
!   Fast physics (from dynamics) are finalized in atmosphere_end above;
!   standard/slow physics (from IPD) are finalized in CCPP_step 'finalize'.
!   The CCPP framework for all cdata structures is finalized in CCPP_step 'finalize'.
    call CCPP_step (step="finalize", nblks=Atm_block%nblks, ierr=ierr)
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP finalize step failed')
#endif

end subroutine atmos_model_end

! </SUBROUTINE>
!#######################################################################
! <SUBROUTINE NAME="atmos_model_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine atmos_model_restart(Atmos, timestamp)
  type (atmos_data_type),   intent(inout) :: Atmos
  character(len=*),  intent(in)           :: timestamp

    call atmosphere_restart(timestamp)
    call FV3GFS_restart_write (IPD_Data, IPD_Restart, Atm_block, &
                               IPD_Control, Atmos%domain, timestamp)

end subroutine atmos_model_restart
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="get_atmos_model_ungridded_dim">
!
! <DESCRIPTION>
!  Retrieve ungridded dimensions of atmospheric model arrays
! </DESCRIPTION>

subroutine get_atmos_model_ungridded_dim(nlev, nsoillev, ntracers,     &
  num_diag_sfc_emis_flux, num_diag_down_flux, num_diag_type_down_flux, &
  num_diag_burn_emis_flux, num_diag_cmass)

  integer, optional, intent(out) :: nlev, nsoillev, ntracers,            &
    num_diag_sfc_emis_flux, num_diag_down_flux, num_diag_type_down_flux, &
    num_diag_burn_emis_flux, num_diag_cmass

  !--- number of atmospheric vertical levels
  if (present(nlev)) nlev = Atm_block%npz

  !--- number of soil levels
  if (present(nsoillev)) then
    nsoillev = 0
    if (allocated(IPD_Data)) then
      if (associated(IPD_Data(1)%Sfcprop%slc)) &
        nsoillev = size(IPD_Data(1)%Sfcprop%slc, dim=2)
    end if
  end if

  !--- total number of atmospheric tracers
  if (present(ntracers)) call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)

  !--- number of tracers used in chemistry diagnostic output
  if (present(num_diag_down_flux)) then
    num_diag_down_flux = 0
    if (associated(IPD_Data(1)%IntDiag%sedim)) &
      num_diag_down_flux = size(IPD_Data(1)%IntDiag%sedim, dim=2)
    if (present(num_diag_type_down_flux)) then
      num_diag_type_down_flux = 0
      if (associated(IPD_Data(1)%IntDiag%sedim))  &
        num_diag_type_down_flux = num_diag_type_down_flux + 1
      if (associated(IPD_Data(1)%IntDiag%drydep)) &
        num_diag_type_down_flux = num_diag_type_down_flux + 1
      if (associated(IPD_Data(1)%IntDiag%wetdpl)) &
        num_diag_type_down_flux = num_diag_type_down_flux + 1
      if (associated(IPD_Data(1)%IntDiag%wetdpc)) &
        num_diag_type_down_flux = num_diag_type_down_flux + 1
    end if
  end if

  !--- number of bins for chemistry diagnostic output
  if (present(num_diag_sfc_emis_flux)) then
    num_diag_sfc_emis_flux = 0
    if (associated(IPD_Data(1)%IntDiag%duem)) &
      num_diag_sfc_emis_flux = size(IPD_Data(1)%IntDiag%duem, dim=2)
    if (associated(IPD_Data(1)%IntDiag%ssem)) &
      num_diag_sfc_emis_flux = &
        num_diag_sfc_emis_flux + size(IPD_Data(1)%IntDiag%ssem, dim=2)
  end if

  !--- number of tracers used in emission diagnostic output
  if (present(num_diag_burn_emis_flux)) then
    num_diag_burn_emis_flux = 0
    if (associated(IPD_Data(1)%IntDiag%abem)) &
      num_diag_burn_emis_flux = size(IPD_Data(1)%IntDiag%abem, dim=2)
  end if

  !--- number of tracers used in column mass density diagnostics
  if (present(num_diag_cmass)) then
    num_diag_cmass = 0
    if (associated(IPD_Data(1)%IntDiag%aecm)) &
      num_diag_cmass = size(IPD_Data(1)%IntDiag%aecm, dim=2)
  end if

end subroutine get_atmos_model_ungridded_dim
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="update_atmos_chemistry">
! <DESCRIPTION>
!  Populate exported chemistry fields with current atmospheric state
!  data (state='export'). Update tracer concentrations for atmospheric
!  chemistry with values from chemistry component (state='import').
!  Fields should be exported/imported from/to the atmospheric state
!  after physics calculations.
!
!  NOTE: It is assumed that all the chemical tracers follow the standard
!  atmospheric tracers, which end with ozone. The order of the chemical
!  tracers must match their order in the chemistry component.
!
!  Requires:
!         IPD_Data
!         Atm_block
! </DESCRIPTION>
subroutine update_atmos_chemistry(state, rc)

  use ESMF
  use module_cplfields,   only: cplFieldGet

  character(len=*),  intent(in)  :: state
  integer, optional, intent(out) :: rc

  !--- local variables
  integer :: localrc
  integer :: ni, nj, nk, nt, ntb, nte
  integer :: nb, ix, i, j, k, it
  integer :: ib, jb

  real(ESMF_KIND_R8), dimension(:,:,:),   pointer :: prsl, phil,  &
                                                     prsi, phii,  &
                                                     temp, dqdt,  &
                                                     ua, va, vvl, &
                                                     dkt, slc,    &
                                                     qb, qm, qu
  real(ESMF_KIND_R8), dimension(:,:,:,:), pointer :: qd, q

  real(ESMF_KIND_R8), dimension(:,:), pointer :: hpbl, area, stype, rainc, &
    uustar, rain, sfcdsw, slmsk, tsfc, shfsfc, snowd, vtype, vfrac, zorl

! logical, parameter :: diag = .true.

  ! -- begin
  if (present(rc)) rc = ESMF_SUCCESS

  ni  = Atm_block%iec - Atm_block%isc + 1
  nj  = Atm_block%jec - Atm_block%jsc + 1
  nk  = Atm_block%npz

  !--- get total number of tracers
  call get_number_tracers(MODEL_ATMOS, num_tracers=nt)

  select case (trim(state))
    case ('import')
      !--- retrieve references to allocated memory for each field
      call cplFieldGet(state,'inst_tracer_mass_frac', farrayPtr4d=q, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      call cplFieldGet(state,'inst_tracer_up_surface_flx', &
        farrayPtr3d=qu, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      call cplFieldGet(state,'inst_tracer_down_surface_flx', &
        farrayPtr4d=qd, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      call cplFieldGet(state,'inst_tracer_clmn_mass_dens', &
        farrayPtr3d=qm, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      call cplFieldGet(state,'inst_tracer_anth_biom_flx', &
        farrayPtr3d=qb, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      !--- do not import tracer concentrations by default
      ntb = nt + 1
      nte = nt

      !--- if chemical tracers are present, set bounds appropriately
      if (IPD_Control%ntchm > 0) then
        if (IPD_Control%ntchs /= NO_TRACER) then
          ntb = IPD_Control%ntchs
          nte = IPD_Control%ntchm + ntb - 1
        end if
      end if

      !--- tracer concentrations
      do it = ntb, nte
!$OMP parallel do default (none) &
!$OMP             shared  (it, nk, nj, ni, Atm_block, IPD_Data, q)  &
!$OMP             private (k, j, jb, i, ib, nb, ix)
        do k = 1, nk
          do j = 1, nj
            jb = j + Atm_block%jsc - 1
            do i = 1, ni
              ib = i + Atm_block%isc - 1
              nb = Atm_block%blkno(ib,jb)
              ix = Atm_block%ixp(ib,jb)
              IPD_Data(nb)%Stateout%gq0(ix,k,it) = q(i,j,k,it)
            enddo
          enddo
        enddo
      enddo

      !--- tracer diagnostics
      !--- (a) column mass densities
      do it = 1, size(qm, dim=3)
!$OMP parallel do default (none) &
!$OMP             shared  (it, nj, ni, Atm_block, IPD_Data, qm)  &
!$OMP             private (j, jb, i, ib, nb, ix)
        do j = 1, nj
          jb = j + Atm_block%jsc - 1
          do i = 1, ni
            ib = i + Atm_block%isc - 1
            nb = Atm_block%blkno(ib,jb)
            ix = Atm_block%ixp(ib,jb)
            IPD_Data(nb)%IntDiag%aecm(ix,it) = qm(i,j,it)
          enddo
        enddo
      enddo

      !--- (b) dust and sea salt emissions
      ntb = size(IPD_Data(1)%IntDiag%duem, dim=2)
      nte = size(qu, dim=3)
      do it = 1, min(ntb, nte)
        do j = 1, nj
          jb = j + Atm_block%jsc - 1
          do i = 1, ni
            ib = i + Atm_block%isc - 1
            nb = Atm_block%blkno(ib,jb)
            ix = Atm_block%ixp(ib,jb)
            IPD_Data(nb)%IntDiag%duem(ix,it) = qu(i,j,it)
          enddo
        enddo
      enddo

      nte = nte - ntb
      do it = 1, min(size(IPD_Data(1)%IntDiag%ssem, dim=2), nte)
        do j = 1, nj
          jb = j + Atm_block%jsc - 1
          do i = 1, ni
            ib = i + Atm_block%isc - 1
            nb = Atm_block%blkno(ib,jb)
            ix = Atm_block%ixp(ib,jb)
            IPD_Data(nb)%IntDiag%ssem(ix,it) = qu(i,j,it+ntb)
          enddo
        enddo
      enddo

      !--- (c) sedimentation and dry/wet deposition
      do it = 1, size(qd, dim=3)
!$OMP parallel do default (none) &
!$OMP             shared  (it, nj, ni, Atm_block, IPD_Data, qd)  &
!$OMP             private (j, jb, i, ib, nb, ix)
        do j = 1, nj
          jb = j + Atm_block%jsc - 1
          do i = 1, ni
            ib = i + Atm_block%isc - 1
            nb = Atm_block%blkno(ib,jb)
            ix = Atm_block%ixp(ib,jb)
            IPD_Data(nb)%IntDiag%sedim (ix,it) = qd(i,j,it,1)
            IPD_Data(nb)%IntDiag%drydep(ix,it) = qd(i,j,it,2)
            IPD_Data(nb)%IntDiag%wetdpl(ix,it) = qd(i,j,it,3)
            IPD_Data(nb)%IntDiag%wetdpc(ix,it) = qd(i,j,it,4)
          enddo
        enddo
      enddo

      !--- (d) anthropogenic and biomass burning emissions
      do it = 1, size(qb, dim=3)
!$OMP parallel do default (none) &
!$OMP             shared  (it, nj, ni, Atm_block, IPD_Data, qb)  &
!$OMP             private (j, jb, i, ib, nb, ix)
        do j = 1, nj
          jb = j + Atm_block%jsc - 1
          do i = 1, ni
            ib = i + Atm_block%isc - 1
            nb = Atm_block%blkno(ib,jb)
            ix = Atm_block%ixp(ib,jb)
            IPD_Data(nb)%IntDiag%abem(ix,it) = qb(i,j,it)
          enddo
        enddo
      enddo

      if (IPD_Control%debug) then
        write(6,'("update_atmos: ",a,": qgrs - min/max/avg",3g16.6)') &
          trim(state), minval(q), maxval(q), sum(q)/size(q)
        write(6,'("update_atmos: ",a,": qup  - min/max/avg",3g16.6)') &
          trim(state), minval(qu), maxval(qu), sum(qu)/size(qu)
        write(6,'("update_atmos: ",a,": qdwn - min/max/avg",3g16.6)') &
          trim(state), minval(qd), maxval(qd), sum(qd)/size(qd)
        write(6,'("update_atmos: ",a,": qcmd - min/max/avg",3g16.6)') &
          trim(state), minval(qm), maxval(qm), sum(qm)/size(qm)
        write(6,'("update_atmos: ",a,": qabb - min/max/avg",3g16.6)') &
          trim(state), minval(qb), maxval(qb), sum(qb)/size(qb)
      end if

    case ('export')
      !--- retrieve references to allocated memory for each field
      call cplFieldGet(state,'inst_pres_interface', farrayPtr3d=prsi, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_pres_levels', farrayPtr3d=prsl, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_geop_interface', farrayPtr3d=phii, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_geop_levels', farrayPtr3d=phil, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_temp_levels', farrayPtr3d=temp, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_zonal_wind_levels', farrayPtr3d=ua, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_merid_wind_levels', farrayPtr3d=va, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_omega_levels', farrayPtr3d=vvl, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_spec_humid_conv_tendency_levels', &
                       farrayPtr3d=dqdt, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_tracer_mass_frac', farrayPtr4d=q, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_soil_moisture_content', &
                       farrayPtr3d=slc, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'soil_type', farrayPtr2d=stype, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_pbl_height', farrayPtr2d=hpbl, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'surface_cell_area', farrayPtr2d=area, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_convective_rainfall_amount', &
                       farrayPtr2d=rainc, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_exchange_coefficient_heat_levels', &
                       farrayPtr3d=dkt, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_friction_velocity', farrayPtr2d=uustar, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_rainfall_amount', farrayPtr2d=rain, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_down_sw_flx', farrayPtr2d=sfcdsw, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_land_sea_mask', farrayPtr2d=slmsk, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_temp_height_surface', farrayPtr2d=tsfc, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_up_sensi_heat_flx', farrayPtr2d=shfsfc, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_lwe_snow_thickness', farrayPtr2d=snowd, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'vegetation_type', farrayPtr2d=vtype, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_vegetation_area_frac', farrayPtr2d=vfrac, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_surface_roughness', farrayPtr2d=zorl, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      !--- handle all three-dimensional variables
!$OMP parallel do default (none) &
!$OMP             shared  (nk, nj, ni, Atm_block, IPD_Data, prsi, phii, prsl, phil, temp, ua, va, vvl, dkt, dqdt)  &
!$OMP             private (k, j, jb, i, ib, nb, ix)
      do k = 1, nk
        do j = 1, nj
          jb = j + Atm_block%jsc - 1
          do i = 1, ni
            ib = i + Atm_block%isc - 1
            nb = Atm_block%blkno(ib,jb)
            ix = Atm_block%ixp(ib,jb)
            !--- interface values
            prsi(i,j,k) = IPD_Data(nb)%Statein%prsi(ix,k)
            phii(i,j,k) = IPD_Data(nb)%Statein%phii(ix,k)
            !--- layer values
            prsl(i,j,k) = IPD_Data(nb)%Statein%prsl(ix,k)
            phil(i,j,k) = IPD_Data(nb)%Statein%phil(ix,k)
            temp(i,j,k) = IPD_Data(nb)%Stateout%gt0(ix,k)
            ua  (i,j,k) = IPD_Data(nb)%Stateout%gu0(ix,k)
            va  (i,j,k) = IPD_Data(nb)%Stateout%gv0(ix,k)
            vvl (i,j,k) = IPD_Data(nb)%Statein%vvl (ix,k)
            dkt (i,j,k) = IPD_Data(nb)%Coupling%dkt(ix,k)
            dqdt(i,j,k) = IPD_Data(nb)%Coupling%dqdti(ix,k)
          enddo
        enddo
      enddo

      !--- top interface values
      k = nk+1
!$omp parallel do default(shared) private(i,j,jb,ib,nb,ix)
      do j = 1, nj
        jb = j + Atm_block%jsc - 1
        do i = 1, ni
          ib = i + Atm_block%isc - 1
          nb = Atm_block%blkno(ib,jb)
          ix = Atm_block%ixp(ib,jb)
          prsi(i,j,k) = IPD_Data(nb)%Statein%prsi(ix,k)
          phii(i,j,k) = IPD_Data(nb)%Statein%phii(ix,k)
        enddo
      enddo

      !--- tracers quantities
      do it = 1, nt
!$OMP parallel do default (none) &
!$OMP             shared  (it, nk, nj, ni, Atm_block, IPD_Data, q)  &
!$OMP             private (k, j, jb, i, ib, nb, ix)
        do k = 1, nk
          do j = 1, nj
            jb = j + Atm_block%jsc - 1
            do i = 1, ni
              ib = i + Atm_block%isc - 1
              nb = Atm_block%blkno(ib,jb)
              ix = Atm_block%ixp(ib,jb)
              q(i,j,k,it) = IPD_Data(nb)%Stateout%gq0(ix,k,it)
            enddo
          enddo
        enddo
      enddo

!$OMP parallel do default (none) &
!$OMP             shared  (nj, ni, Atm_block, IPD_Data, &
!$OMP                      hpbl, area, stype, rainc, rain, uustar, sfcdsw, &
!$OMP                      slmsk, snowd, tsfc, shfsfc, vtype, vfrac, zorl, slc) &
!$OMP             private (j, jb, i, ib, nb, ix)
      do j = 1, nj
        jb = j + Atm_block%jsc - 1
        do i = 1, ni
          ib = i + Atm_block%isc - 1
          nb = Atm_block%blkno(ib,jb)
          ix = Atm_block%ixp(ib,jb)
          hpbl(i,j)   = IPD_Data(nb)%Tbd%hpbl(ix)
          area(i,j)   = IPD_Data(nb)%Grid%area(ix)
          stype(i,j)  = IPD_Data(nb)%Sfcprop%stype(ix)
          rainc(i,j)  = IPD_Data(nb)%Coupling%rainc_cpl(ix)
          rain(i,j)   = IPD_Data(nb)%Coupling%rain_cpl(ix)  &
                      + IPD_Data(nb)%Coupling%snow_cpl(ix)
          uustar(i,j) = IPD_Data(nb)%Sfcprop%uustar(ix)
          sfcdsw(i,j) = IPD_Data(nb)%Coupling%sfcdsw(ix)
          slmsk(i,j)  = IPD_Data(nb)%Sfcprop%slmsk(ix)
          snowd(i,j)  = IPD_Data(nb)%Sfcprop%snowd(ix)
          tsfc(i,j)   = IPD_Data(nb)%Sfcprop%tsfc(ix)
          shfsfc(i,j) = IPD_Data(nb)%Coupling%ushfsfci(ix)
          vtype(i,j)  = IPD_Data(nb)%Sfcprop%vtype(ix)
          vfrac(i,j)  = IPD_Data(nb)%Sfcprop%vfrac(ix)
          zorl(i,j)   = IPD_Data(nb)%Sfcprop%zorl(ix)
          slc(i,j,:)  = IPD_Data(nb)%Sfcprop%slc(ix,:)
        enddo
      enddo

      ! -- zero out accumulated fields
!$OMP parallel do default (none) &
!$OMP             shared  (nj, ni, Atm_block, IPD_Control, IPD_Data) &
!$OMP             private (j, jb, i, ib, nb, ix)
      do j = 1, nj
        jb = j + Atm_block%jsc - 1
        do i = 1, ni
          ib = i + Atm_block%isc - 1
          nb = Atm_block%blkno(ib,jb)
          ix = Atm_block%ixp(ib,jb)
          IPD_Data(nb)%coupling%rainc_cpl(ix)  = zero
          if (.not.IPD_Control%cplflx) then
            IPD_Data(nb)%coupling%rain_cpl(ix) = zero
            IPD_Data(nb)%coupling%snow_cpl(ix) = zero
          end if
        enddo
      enddo

      if (IPD_Control%debug) then
        ! -- diagnostics
        write(6,'("update_atmos: prsi   - min/max/avg",3g16.6)') minval(prsi),   maxval(prsi),   sum(prsi)/size(prsi)
        write(6,'("update_atmos: phii   - min/max/avg",3g16.6)') minval(phii),   maxval(phii),   sum(phii)/size(phii)
        write(6,'("update_atmos: prsl   - min/max/avg",3g16.6)') minval(prsl),   maxval(prsl),   sum(prsl)/size(prsl)
        write(6,'("update_atmos: phil   - min/max/avg",3g16.6)') minval(phil),   maxval(phil),   sum(phil)/size(phil)
        write(6,'("update_atmos: tgrs   - min/max/avg",3g16.6)') minval(temp),   maxval(temp),   sum(temp)/size(temp)
        write(6,'("update_atmos: ugrs   - min/max/avg",3g16.6)') minval(ua),     maxval(ua),     sum(ua)/size(ua)
        write(6,'("update_atmos: vgrs   - min/max/avg",3g16.6)') minval(va),     maxval(va),     sum(va)/size(va)
        write(6,'("update_atmos: vvl    - min/max/avg",3g16.6)') minval(vvl),    maxval(vvl),    sum(vvl)/size(vvl)
        write(6,'("update_atmos: dqdt   - min/max/avg",3g16.6)') minval(dqdt),   maxval(dqdt),   sum(dqdt)/size(dqdt)
        write(6,'("update_atmos: qgrs   - min/max/avg",3g16.6)') minval(q),      maxval(q),      sum(q)/size(q)

        write(6,'("update_atmos: hpbl   - min/max/avg",3g16.6)') minval(hpbl),   maxval(hpbl),   sum(hpbl)/size(hpbl)
        write(6,'("update_atmos: rainc  - min/max/avg",3g16.6)') minval(rainc),  maxval(rainc),  sum(rainc)/size(rainc)
        write(6,'("update_atmos: rain   - min/max/avg",3g16.6)') minval(rain),   maxval(rain),   sum(rain)/size(rain)
        write(6,'("update_atmos: shfsfc - min/max/avg",3g16.6)') minval(shfsfc), maxval(shfsfc), sum(shfsfc)/size(shfsfc)
        write(6,'("update_atmos: sfcdsw - min/max/avg",3g16.6)') minval(sfcdsw), maxval(sfcdsw), sum(sfcdsw)/size(sfcdsw)
        write(6,'("update_atmos: slmsk  - min/max/avg",3g16.6)') minval(slmsk),  maxval(slmsk),  sum(slmsk)/size(slmsk)
        write(6,'("update_atmos: snowd  - min/max/avg",3g16.6)') minval(snowd),  maxval(snowd),  sum(snowd)/size(snowd)
        write(6,'("update_atmos: tsfc   - min/max/avg",3g16.6)') minval(tsfc),   maxval(tsfc),   sum(tsfc)/size(tsfc)
        write(6,'("update_atmos: vtype  - min/max/avg",3g16.6)') minval(vtype),  maxval(vtype),  sum(vtype)/size(vtype)
        write(6,'("update_atmos: vfrac  - min/max/avg",3g16.6)') minval(vfrac),  maxval(vfrac),  sum(vfrac)/size(vfrac)
        write(6,'("update_atmos: area   - min/max/avg",3g16.6)') minval(area),   maxval(area),   sum(area)/size(area)
        write(6,'("update_atmos: stype  - min/max/avg",3g16.6)') minval(stype),  maxval(stype),  sum(stype)/size(stype)
        write(6,'("update_atmos: zorl   - min/max/avg",3g16.6)') minval(zorl),   maxval(zorl),   sum(zorl)/size(zorl)
        write(6,'("update_atmos: slc    - min/max/avg",3g16.6)') minval(slc),    maxval(slc),    sum(slc)/size(slc)
      end if

    case default
      ! -- do nothing
  end select

end subroutine update_atmos_chemistry
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="atmos_data_type_chksum">
!
! <OVERVIEW>
!  Print checksums of the various fields in the atmos_data_type.
! </OVERVIEW>

! <DESCRIPTION>
!  Routine to print checksums of the various fields in the atmos_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_data_type_chksum(id, timestep, atm)
! </TEMPLATE>

! <IN NAME="Atm" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields in the atmos_data_type.
! </INOUT>
!
! <IN NAME="id" TYPE="character">
!   Label to differentiate where this routine in being called from.
! </IN>
!
! <IN NAME="timestep" TYPE="integer">
!   An integer to indicate which timestep this routine is being called for.
! </IN>
!
subroutine atmos_data_type_chksum(id, timestep, atm)
type(atmos_data_type), intent(in) :: atm 
    character(len=*),  intent(in) :: id
    integer         ,  intent(in) :: timestep
    integer :: n, outunit

100 format("CHECKSUM::",A32," = ",Z20)
101 format("CHECKSUM::",A16,a,'%',a," = ",Z20)

  outunit = stdout()
  write(outunit,*) 'BEGIN CHECKSUM(Atmos_data_type):: ', id, timestep
  write(outunit,100) ' atm%lon_bnd                ', mpp_chksum(atm%lon_bnd)
  write(outunit,100) ' atm%lat_bnd                ', mpp_chksum(atm%lat_bnd)
  write(outunit,100) ' atm%lon                    ', mpp_chksum(atm%lon)
  write(outunit,100) ' atm%lat                    ', mpp_chksum(atm%lat)

end subroutine atmos_data_type_chksum

! </SUBROUTINE>

  subroutine alloc_atmos_data_type (nlon, nlat, Atmos)
   integer, intent(in) :: nlon, nlat
   type(atmos_data_type), intent(inout) :: Atmos
    allocate ( Atmos % lon_bnd  (nlon+1,nlat+1), &
               Atmos % lat_bnd  (nlon+1,nlat+1), &
               Atmos % lon      (nlon,nlat),     &
               Atmos % lat      (nlon,nlat)      )

  end subroutine alloc_atmos_data_type

  subroutine dealloc_atmos_data_type (Atmos)
   type(atmos_data_type), intent(inout) :: Atmos
    deallocate (Atmos%lon_bnd, &
                Atmos%lat_bnd, &
                Atmos%lon,     &
                Atmos%lat      )
  end subroutine dealloc_atmos_data_type

  subroutine assign_importdata(rc)

    use module_cplfields,  only: importFields, nImportFields, QueryFieldList, &
                                 ImportFieldsList, importFieldsValid
    use ESMF
!
    implicit none
    integer, intent(out) :: rc

    !--- local variables
    integer :: n, j, i, ix, nb, isc, iec, jsc, jec, dimCount, findex
    character(len=128) :: impfield_name, fldname
    type(ESMF_TypeKind_Flag)                           :: datatype
    real(kind=ESMF_KIND_R4),  dimension(:,:), pointer  :: datar42d
    real(kind=ESMF_KIND_R8),  dimension(:,:), pointer  :: datar82d
    real(kind=IPD_kind_phys), dimension(:,:), pointer  :: datar8
    logical found, isFieldCreated, lcpl_fice
!
!------------------------------------------------------------------------------
!
! set up local dimension
    rc  = -999
    isc = IPD_control%isc
    iec = IPD_control%isc+IPD_control%nx-1
    jsc = IPD_control%jsc
    jec = IPD_control%jsc+IPD_control%ny-1
    lcpl_fice = .false.

    allocate(datar8(isc:iec,jsc:jec))

!   if (mpp_pe() == mpp_root_pe() .and. debug) print *,'in cplImp,dim=',isc,iec,jsc,jec
!   if (mpp_pe() == mpp_root_pe() .and. debug) print *,'in cplImp,IPD_Data, size', size(IPD_Data)
!   if (mpp_pe() == mpp_root_pe() .and. debug) print *,'in cplImp,tsfc, size', size(IPD_Data(1)%sfcprop%tsfc)
!   if (mpp_pe() == mpp_root_pe() .and. debug) print *,'in cplImp,tsfc, min_seaice', IPD_Control%min_seaice

    do n=1,nImportFields ! Each import field is only available if it was connected in the import state.

      found = .false.

      isFieldCreated = ESMF_FieldIsCreated(importFields(n), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if (isFieldCreated) then ! put the data from local cubed sphere grid to column grid for phys

        datar8 = -99999.0
        call ESMF_FieldGet(importFields(n), dimCount=dimCount ,typekind=datatype, &
                           name=impfield_name, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if ( dimCount == 2) then
          if ( datatype == ESMF_TYPEKIND_R8) then
            call ESMF_FieldGet(importFields(n),farrayPtr=datar82d,localDE=0, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            datar8 = datar82d
            if (mpp_pe() == mpp_root_pe() .and. debug) print *,'in cplIMP,atmos gets ',trim(impfield_name),' datar8=', &
                                                               datar8(isc,jsc), maxval(datar8), minval(datar8)
            found = .true.
! gfs physics runs with r8
!          else
!            call ESMF_FieldGet(importFields(n),farrayPtr=datar42d,localDE=0, rc=rc)
!            datar8 = datar42d
          endif
        endif
!
        if (found .and. datar8(isc,jsc) > -99998.0) then
!
        ! get sea land mask: in order to update the coupling fields over the ocean/ice
!        fldname = 'land_mask'
!        if (trim(impfield_name) == trim(fldname)) then
!          findex = QueryFieldList(ImportFieldsList,fldname)
!          if (importFieldsValid(findex)) then
!!$omp parallel do default(shared) private(i,j,nb,ix)
!            do j=jsc,jec
!              do i=isc,iec
!                nb = Atm_block%blkno(i,j)
!                ix = Atm_block%ixp(i,j)
!                IPD_Data(nb)%Coupling%slimskin_cpl(ix) = datar8(i,j)
!              enddo
!            enddo
!            if( mpp_pe()==mpp_root_pe()) print *,'get land mask from mediator'
!          endif
!        endif

! get sea ice surface temperature
!--------------------------------
          fldname = 'sea_ice_surface_temperature'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = QueryFieldList(ImportFieldsList,fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  IPD_Data(nb)%Coupling%tisfcin_cpl(ix) = datar8(i,j)
                enddo
              enddo
            endif
          endif

! get sst:  sst needs to be adjusted by land sea mask before passing to fv3
!--------------------------------------------------------------------------
          fldname = 'sea_surface_temperature'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = QueryFieldList(ImportFieldsList,fldname)
!       if (mpp_pe() == mpp_root_pe() .and. debug)  print *,' for sst', &
!    ' fldname=',fldname,' findex=',findex,' importFieldsValid=',importFieldsValid(findex)

            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (IPD_Data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    IPD_Data(nb)%Coupling%tseain_cpl(ix) = datar8(i,j)
                    IPD_Data(nb)%Sfcprop%tsfco(ix)       = datar8(i,j)
!                   IPD_Data(nb)%Sfcprop%tsfc(ix)        = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'get sst from mediator'
            endif
          endif

! get sea ice fraction:  fice or sea ice concentration from the mediator
!-----------------------------------------------------------------------
          fldname = 'ice_fraction'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = QueryFieldList(ImportFieldsList,fldname)
            if (importFieldsValid(findex)) then
            lcpl_fice = .true.
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  IPD_Data(nb)%Coupling%ficein_cpl(ix)   = zero
                  if (IPD_Data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    if (datar8(i,j) >= IPD_control%min_seaice*IPD_Data(nb)%Sfcprop%oceanfrac(ix)) then
                      IPD_Data(nb)%Coupling%ficein_cpl(ix) = max(zero, min(datar8(i,j),one))
!                     if (IPD_Data(nb)%Sfcprop%oceanfrac(ix) == one) IPD_Data(nb)%Sfcprop%slmsk(ix) = 2. !slmsk=2 crashes in gcycle on partial land points
                      IPD_Data(nb)%Sfcprop%slmsk(ix)         = 2.                                        !slmsk=2 crashes in gcycle on partial land points
                      IPD_Data(nb)%Coupling%slimskin_cpl(ix) = 4.
                    else
                      if (IPD_Data(nb)%Sfcprop%oceanfrac(ix) == one) IPD_Data(nb)%Sfcprop%slmsk(ix) = zero
                      IPD_Data(nb)%Coupling%slimskin_cpl(ix) = zero
                    endif
                  else
                    IPD_Data(nb)%Sfcprop%slmsk(ix)         = one
                    IPD_Data(nb)%Coupling%slimskin_cpl(ix) = one
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get fice from mediator'
            endif
          endif

! get upward LW flux:  for sea ice covered area
!----------------------------------------------
          fldname = 'mean_up_lw_flx'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = QueryFieldList(ImportFieldsList,fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
!               do i=isc,iec
!                 nb = Atm_block%blkno(i,j)
!                 ix = Atm_block%ixp(i,j)
!                if (IPD_Data(nb)%Sfcprop%slmsk(ix) < 0.1 .or. IPD_Data(nb)%Sfcprop%slmsk(ix) > 1.9) then
!                   IPD_Data(nb)%Coupling%ulwsfcin_cpl(ix) = -datar8(i,j)
!                 endif
!               enddo
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (IPD_Data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    IPD_Data(nb)%Coupling%ulwsfcin_cpl(ix) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get lwflx from mediator'
            endif
          endif

! get latent heat flux:  for sea ice covered area
!------------------------------------------------
          fldname = 'mean_laten_heat_flx'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = QueryFieldList(ImportFieldsList,fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (IPD_Data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    IPD_Data(nb)%Coupling%dqsfcin_cpl(ix) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get laten_heat from mediator'
            endif
          endif

! get sensible heat flux:  for sea ice covered area
!--------------------------------------------------
          fldname = 'mean_sensi_heat_flx'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = QueryFieldList(ImportFieldsList,fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (IPD_Data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    IPD_Data(nb)%Coupling%dtsfcin_cpl(ix) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get sensi_heat from mediator'
            endif
          endif

! get zonal compt of momentum flux:  for sea ice covered area
!------------------------------------------------------------
          fldname = 'mean_zonal_moment_flx'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = QueryFieldList(ImportFieldsList,fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (IPD_Data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    IPD_Data(nb)%Coupling%dusfcin_cpl(ix) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get zonal_moment_flx from mediator'
            endif
          endif

! get meridional compt of momentum flux:  for sea ice covered area
!-----------------------------------------------------------------
          fldname = 'mean_merid_moment_flx'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = QueryFieldList(ImportFieldsList,fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (IPD_Data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    IPD_Data(nb)%Coupling%dvsfcin_cpl(ix) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get merid_moment_flx from mediator'
            endif
          endif

! get sea ice volume:  for sea ice covered area
!----------------------------------------------
          fldname = 'mean_ice_volume'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = QueryFieldList(ImportFieldsList,fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (IPD_Data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    IPD_Data(nb)%Coupling%hicein_cpl(ix) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug) print *,'fv3 assign_import: get ice_volume from mediator'
            endif
          endif

! get snow volume:  for sea ice covered area
!-------------------------------------------
          fldname = 'mean_snow_volume'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = QueryFieldList(ImportFieldsList,fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (IPD_Data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    IPD_Data(nb)%Coupling%hsnoin_cpl(ix) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get snow_volume from mediator'
            endif
          endif

        endif ! if (datar8(isc,jsc) > -99999.0) then
      endif   ! if (isFieldCreated) then
    enddo
!
    deallocate(datar8)

! update sea ice related fields:
    if( lcpl_fice ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          if (IPD_Data(nb)%Sfcprop%oceanfrac(ix) > zero) then
!if it is ocean or ice get surface temperature from mediator
            if(IPD_Data(nb)%Coupling%ficein_cpl(ix) >= IPD_control%min_seaice*IPD_Data(nb)%Sfcprop%oceanfrac(ix)) then
              IPD_Data(nb)%Sfcprop%tisfc(ix) = IPD_Data(nb)%Coupling%tisfcin_cpl(ix)
              IPD_Data(nb)%Sfcprop%fice(ix)  = IPD_Data(nb)%Coupling%ficein_cpl(ix)
              IPD_Data(nb)%Sfcprop%hice(ix)  = IPD_Data(nb)%Coupling%hicein_cpl(ix)
              IPD_Data(nb)%Sfcprop%snowd(ix) = IPD_Data(nb)%Coupling%hsnoin_cpl(ix)
            else 
              IPD_Data(nb)%Sfcprop%fice(ix)  = zero
              IPD_Data(nb)%Sfcprop%hice(ix)  = zero
              IPD_Data(nb)%Sfcprop%snowd(ix) = zero
!
              IPD_Data(nb)%Coupling%dtsfcin_cpl(ix)  = -99999.0 ! over open water - should not be used in ATM
              IPD_Data(nb)%Coupling%dqsfcin_cpl(ix)  = -99999.0 !                 ,,
              IPD_Data(nb)%Coupling%dusfcin_cpl(ix)  = -99999.0 !                 ,,
              IPD_Data(nb)%Coupling%dvsfcin_cpl(ix)  = -99999.0 !                 ,,
              IPD_Data(nb)%Coupling%dtsfcin_cpl(ix)  = -99999.0 !                 ,,
              IPD_Data(nb)%Coupling%ulwsfcin_cpl(ix) = -99999.0 !                 ,,
              if (IPD_Data(nb)%Sfcprop%oceanfrac(ix) == one) IPD_Data(nb)%Sfcprop%slmsk(ix) = zero ! 100% open water
            endif
          endif
        enddo
      enddo
    endif

    rc=0
!
    if (mpp_pe() == mpp_root_pe()) print *,'end of assign_importdata'
  end subroutine assign_importdata

!
  subroutine setup_exportdata (rc)

    use module_cplfields,  only: exportData, nExportFields, exportFieldsList, &
                                 queryFieldList, fillExportFields

    implicit none

!------------------------------------------------------------------------------

    !--- interface variables
    integer, intent(out) :: rc

    !--- local variables
    integer                :: j, i, ix, nb, isc, iec, jsc, jec, idx
    real(IPD_kind_phys)    :: rtime, rtimek
!
!   if (mpp_pe() == mpp_root_pe()) print *,'enter setup_exportdata'

    isc = IPD_control%isc
    iec = IPD_control%isc+IPD_control%nx-1
    jsc = IPD_control%jsc
    jec = IPD_control%jsc+IPD_control%ny-1

    rtime  = one / IPD_control%dtp
    rtimek = IPD_control%rho_h2o * rtime
!    print *,'in cplExp,dim=',isc,iec,jsc,jec,'nExportFields=',nExportFields
!    print *,'in cplExp,IPD_Data, size', size(IPD_Data)
!    print *,'in cplExp,u10micpl, size', size(IPD_Data(1)%coupling%u10mi_cpl)

    if(.not.allocated(exportData)) then
      allocate(exportData(isc:iec,jsc:jec,nExportFields))
    endif

    ! set cpl fields to export Data

    if (IPD_Control%cplflx .or. IPD_Control%cplwav) then 
    ! Instantaneous u wind (m/s) 10 m above ground
    idx = queryfieldlist(exportFieldsList,'inst_zonal_wind_height10m')
    if (idx > 0 ) then
      if (mpp_pe() == mpp_root_pe() .and. debug) print *,'cpl, in get u10mi_cpl'
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%u10mi_cpl(ix)
        enddo
      enddo
    endif

    ! Instantaneous v wind (m/s) 10 m above ground
    idx = queryfieldlist(exportFieldsList,'inst_merid_wind_height10m')
    if (idx > 0 ) then
      if (mpp_pe() == mpp_root_pe() .and. debug) print *,'cpl, in get v10mi_cpl'
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%v10mi_cpl(ix)
        enddo
      enddo
      if (mpp_pe() == mpp_root_pe() .and. debug) print *,'cpl, get v10mi_cpl, exportData=',exportData(isc,jsc,idx),'idx=',idx
    endif

    endif !if cplflx or cplwav 

    if (IPD_Control%cplflx) then
    ! MEAN Zonal compt of momentum flux (N/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_zonal_moment_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dusfc_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN Merid compt of momentum flux (N/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_merid_moment_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dvsfc_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN Sensible heat flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_sensi_heat_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dtsfc_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN Latent heat flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_laten_heat_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dqsfc_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN Downward LW heat flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_down_lw_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dlwsfc_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN Downward SW heat flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_down_sw_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dswsfc_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN precipitation rate (kg/m2/s)
    idx = queryfieldlist(exportFieldsList,'mean_prec_rate')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%rain_cpl(ix) * rtimek
        enddo
      enddo
    endif

    ! Instataneous Zonal compt of momentum flux (N/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_zonal_moment_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dusfci_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Merid compt of momentum flux (N/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_merid_moment_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dvsfci_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Sensible heat flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_sensi_heat_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dtsfci_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Latent heat flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_laten_heat_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dqsfci_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Downward long wave radiation flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_down_lw_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dlwsfci_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Downward solar radiation flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_down_sw_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dswsfci_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Temperature (K) 2 m above ground
    idx = queryfieldlist(exportFieldsList,'inst_temp_height2m')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%t2mi_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Specific humidity (kg/kg) 2 m above ground
    idx = queryfieldlist(exportFieldsList,'inst_spec_humid_height2m')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%q2mi_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Temperature (K) at surface
    idx = queryfieldlist(exportFieldsList,'inst_temp_height_surface')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%tsfci_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Pressure (Pa) land and sea surface
    idx = queryfieldlist(exportFieldsList,'inst_pres_height_surface')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%psurfi_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous Surface height (m)
    idx = queryfieldlist(exportFieldsList,'inst_surface_height')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%oro_cpl(ix)
        enddo
      enddo
    endif

    ! MEAN NET long wave radiation flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_net_lw_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nlwsfc_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN NET solar radiation flux over the ocean (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_net_sw_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nswsfc_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! Instataneous NET long wave radiation flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_net_lw_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nlwsfci_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous NET solar radiation flux over the ocean (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_net_sw_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nswsfci_cpl(ix)
        enddo
      enddo
    endif

    ! MEAN sfc downward nir direct flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_down_sw_ir_dir_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dnirbm_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN sfc downward nir diffused flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_down_sw_ir_dif_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dnirdf_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN sfc downward uv+vis direct flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_down_sw_vis_dir_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dvisbm_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN sfc downward uv+vis diffused flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_down_sw_vis_dif_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dvisdf_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! Instataneous sfc downward nir direct flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_down_sw_ir_dir_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dnirbmi_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous sfc downward nir diffused flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_down_sw_ir_dif_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dnirdfi_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous sfc downward uv+vis direct flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_down_sw_vis_dir_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dvisbmi_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous sfc downward uv+vis diffused flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_down_sw_vis_dif_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%dvisdfi_cpl(ix)
        enddo
      enddo
    endif

    ! MEAN NET sfc nir direct flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_net_sw_ir_dir_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nnirbm_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN NET sfc nir diffused flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_net_sw_ir_dif_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nnirdf_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN NET sfc uv+vis direct flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_net_sw_vis_dir_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nvisbm_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! MEAN NET sfc uv+vis diffused flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'mean_net_sw_vis_dif_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nvisdf_cpl(ix) * rtime
        enddo
      enddo
    endif

    ! Instataneous net sfc nir direct flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_net_sw_ir_dir_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nnirbmi_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous net sfc nir diffused flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_net_sw_ir_dif_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nnirdfi_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous net sfc uv+vis direct flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_net_sw_vis_dir_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nvisbmi_cpl(ix)
        enddo
      enddo
    endif

    ! Instataneous net sfc uv+vis diffused flux (W/m**2)
    idx = queryfieldlist(exportFieldsList,'inst_net_sw_vis_dif_flx')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%nvisdfi_cpl(ix)
        enddo
      enddo
    endif

    ! Land/Sea mask (sea:0,land:1)
    idx = queryfieldlist(exportFieldsList,'inst_land_sea_mask')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%slmsk_cpl(ix)
        enddo
      enddo
    endif

! Data from DYCORE:

    ! bottom layer temperature (t)
    idx = queryfieldlist(exportFieldsList,'inst_temp_height_lowest')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          if (associated(DYCORE_Data(nb)%coupling%t_bot)) then 
            exportData(i,j,idx) = DYCORE_Data(nb)%coupling%t_bot(ix)
          else 
            exportData(i,j,idx) = zero
          endif 
        enddo
      enddo
    endif

    ! bottom layer specific humidity (q)
    !!! CHECK if tracer 1 is for specific humidity !!!
    idx = queryfieldlist(exportFieldsList,'inst_spec_humid_height_lowest')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          if (associated(DYCORE_Data(nb)%coupling%tr_bot)) then
            exportData(i,j,idx) = DYCORE_Data(nb)%coupling%tr_bot(ix,1)
          else 
            exportData(i,j,idx) = zero
          endif 
        enddo
      enddo
    endif

    ! bottom layer zonal wind (u)
    idx = queryfieldlist(exportFieldsList,'inst_zonal_wind_height_lowest')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          if (associated(DYCORE_Data(nb)%coupling%u_bot)) then
            exportData(i,j,idx) = DYCORE_Data(nb)%coupling%u_bot(ix)
          else
            exportData(i,j,idx) = zero
          endif 
        enddo
      enddo
    endif

    ! bottom layer meridionalw wind (v)
    idx = queryfieldlist(exportFieldsList,'inst_merid_wind_height_lowest')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          if (associated(DYCORE_Data(nb)%coupling%v_bot)) then
            exportData(i,j,idx) = DYCORE_Data(nb)%coupling%v_bot(ix)
          else 
            exportData(i,j,idx) = zero 
          endif 
        enddo
      enddo
    endif

    ! bottom layer pressure (p)
    idx = queryfieldlist(exportFieldsList,'inst_pres_height_lowest')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          if (associated(DYCORE_Data(nb)%coupling%p_bot)) then
            exportData(i,j,idx) = DYCORE_Data(nb)%coupling%p_bot(ix)
          else 
            exportData(i,j,idx) = zero
          endif 
        enddo
      enddo
    endif

    ! bottom layer height (z)
    idx = queryfieldlist(exportFieldsList,'inst_height_lowest')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          if (associated(DYCORE_Data(nb)%coupling%z_bot)) then
            exportData(i,j,idx) = DYCORE_Data(nb)%coupling%z_bot(ix)
          else 
            exportData(i,j,idx) = zero 
          endif 
        enddo
      enddo
    endif

! END Data from DYCORE.

    ! MEAN snow precipitation rate (kg/m2/s)
    idx = queryfieldlist(exportFieldsList,'mean_fprec_rate')
    if (idx > 0 ) then
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          exportData(i,j,idx) = IPD_Data(nb)%coupling%snow_cpl(ix) * rtimek
        enddo
      enddo
    endif
    endif !cplflx 

!---
! Fill the export Fields for ESMF/NUOPC style coupling
    call fillExportFields(exportData)

!---
    if (IPD_Control%cplflx) then 
! zero out accumulated fields
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          IPD_Data(nb)%coupling%dusfc_cpl(ix)  = zero
          IPD_Data(nb)%coupling%dvsfc_cpl(ix)  = zero
          IPD_Data(nb)%coupling%dtsfc_cpl(ix)  = zero
          IPD_Data(nb)%coupling%dqsfc_cpl(ix)  = zero
          IPD_Data(nb)%coupling%dlwsfc_cpl(ix) = zero
          IPD_Data(nb)%coupling%dswsfc_cpl(ix) = zero
          IPD_Data(nb)%coupling%rain_cpl(ix)   = zero
          IPD_Data(nb)%coupling%nlwsfc_cpl(ix) = zero
          IPD_Data(nb)%coupling%nswsfc_cpl(ix) = zero
          IPD_Data(nb)%coupling%dnirbm_cpl(ix) = zero
          IPD_Data(nb)%coupling%dnirdf_cpl(ix) = zero
          IPD_Data(nb)%coupling%dvisbm_cpl(ix) = zero
          IPD_Data(nb)%coupling%dvisdf_cpl(ix) = zero
          IPD_Data(nb)%coupling%nnirbm_cpl(ix) = zero
          IPD_Data(nb)%coupling%nnirdf_cpl(ix) = zero
          IPD_Data(nb)%coupling%nvisbm_cpl(ix) = zero
          IPD_Data(nb)%coupling%nvisdf_cpl(ix) = zero
          IPD_Data(nb)%coupling%snow_cpl(ix)   = zero
        enddo
      enddo
      if (mpp_pe() == mpp_root_pe()) print *,'zeroing coupling fields at kdt= ',IPD_Control%kdt
    endif !cplflx
!   if (mpp_pe() == mpp_root_pe()) print *,'end of setup_exportdata'

  end subroutine setup_exportdata

  subroutine addLsmask2grid(fcstgrid, rc)

    use ESMF
!
    implicit none
    type(ESMF_Grid)      :: fcstgrid
    integer, optional, intent(out) :: rc
!
!  local vars
    integer isc, iec, jsc, jec
    integer i, j, nb, ix
!    integer CLbnd(2), CUbnd(2), CCount(2), TLbnd(2), TUbnd(2), TCount(2)
    type(ESMF_StaggerLoc) :: staggerloc
    integer, allocatable :: lsmask(:,:)
    integer(kind=ESMF_KIND_I4), pointer  :: maskPtr(:,:)
!
    isc = IPD_control%isc
    iec = IPD_control%isc+IPD_control%nx-1
    jsc = IPD_control%jsc
    jec = IPD_control%jsc+IPD_control%ny-1
    allocate(lsmask(isc:iec,jsc:jec))
!
!$omp parallel do default(shared) private(i,j,nb,ix)
    do j=jsc,jec
      do i=isc,iec
        nb = Atm_block%blkno(i,j)
        ix = Atm_block%ixp(i,j)
! use land sea mask: land:1, ocean:0
        lsmask(i,j) = floor(IPD_Data(nb)%SfcProp%landfrac(ix))
      enddo
    enddo
!
! Get mask
    call ESMF_GridAddItem(fcstgrid, itemflag=ESMF_GRIDITEM_MASK,   &
                          staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!    call ESMF_GridGetItemBounds(fcstgrid, itemflag=ESMF_GRIDITEM_MASK,   &
!         staggerloc=ESMF_STAGGERLOC_CENTER, computationalLBound=ClBnd,  &
!         computationalUBound=CUbnd, computationalCount=Ccount,  &
!         totalLBound=TLbnd, totalUBound=TUbnd, totalCount=Tcount, rc=rc)
!    print *,'in set up grid, aft add esmfgridadd item mask, rc=',rc, &
!     'ClBnd=',ClBnd,'CUbnd=',CUbnd,'Ccount=',Ccount, &
!     'TlBnd=',TlBnd,'TUbnd=',TUbnd,'Tcount=',Tcount
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridGetItem(fcstgrid, itemflag=ESMF_GRIDITEM_MASK,   &
                          staggerloc=ESMF_STAGGERLOC_CENTER,farrayPtr=maskPtr, rc=rc)
!    print *,'in set up grid, aft get maskptr, rc=',rc, 'size=',size(maskPtr,1),size(maskPtr,2), &
!      'bound(maskPtr)=', LBOUND(maskPtr,1),LBOUND(maskPtr,2),UBOUND(maskPtr,1),UBOUND(maskPtr,2)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!    
!$omp parallel do default(shared) private(i,j)
    do j=jsc,jec
      do i=isc,iec
        maskPtr(i-isc+1,j-jsc+1) = lsmask(i,j)
      enddo
    enddo
!      print *,'in set set lsmask, maskPtr=', maxval(maskPtr), minval(maskPtr)
!
    deallocate(lsmask)  

  end subroutine addLsmask2grid
!------------------------------------------------------------------------------

end module atmos_model_mod
