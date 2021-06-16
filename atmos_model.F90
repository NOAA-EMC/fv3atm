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
use atmosphere_mod,     only: Atm, mygrid
use block_control_mod,  only: block_control_type, define_blocks_packed
use DYCORE_typedefs,    only: DYCORE_data_type, DYCORE_diag_type

use GFS_typedefs,       only: GFS_init_type, GFS_kind_phys => kind_phys
use GFS_restart,        only: GFS_restart_type, GFS_restart_populate
use GFS_diagnostics,    only: GFS_externaldiag_type, &
                              GFS_externaldiag_populate
use CCPP_data,          only: ccpp_suite, GFS_control, &
                              GFS_data, GFS_interstitial
use GFS_init,           only: GFS_initialize
use CCPP_driver,        only: CCPP_step, non_uniform_blocks

use stochastic_physics_wrapper_mod, only: stochastic_physics_wrapper,stochastic_physics_wrapper_end

use FV3GFS_io_mod,      only: FV3GFS_restart_read, FV3GFS_restart_write, &
                              FV3GFS_GFS_checksum,                       &
                              FV3GFS_diag_register, FV3GFS_diag_output,  &
                              DIAG_SIZE
use fv_iau_mod,         only: iau_external_data_type,getiauforcing,iau_initialize
use module_fv3_config,  only: output_1st_tstep_rst, first_kdt, nsout,    &
                              restart_endfcst
use module_block_data

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
public setup_exportdata
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
     real(kind=GFS_kind_phys), pointer, dimension(:,:)   :: lon      => null() ! local longitude axis grid box centers in radians.
     real(kind=GFS_kind_phys), pointer, dimension(:,:)   :: lat      => null() ! local latitude axis grid box centers in radians.
     real(kind=GFS_kind_phys), pointer, dimension(:,:)   :: dx, dy
     real(kind=8),             pointer, dimension(:,:)   :: area
     real(kind=8),             pointer, dimension(:,:,:) :: layer_hgt, level_hgt
     type(domain2d)                :: domain             ! domain decomposition
     type(time_type)               :: Time               ! current time
     type(time_type)               :: Time_step          ! atmospheric time step.
     type(time_type)               :: Time_init          ! reference time.
     type(grid_box_type)           :: grid               ! hold grid information needed for 2nd order conservative flux exchange
     type(GFS_externaldiag_type), pointer, dimension(:) :: Diag
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
namelist /atmos_model_nml/ blocksize, chksum_debug, dycore_only, debug, sync, fdiag, fhmax, fhmaxhf, fhout, fhouthf, ccpp_suite, avg_max_length

type (time_type) :: diag_time, diag_time_fhzero

!--- concurrent and decoupled radiation and physics variables
!-------------------
!  DYCORE containers
!-------------------
type(DYCORE_data_type),    allocatable :: DYCORE_Data(:)  ! number of blocks
type(DYCORE_diag_type)                 :: DYCORE_Diag(25)

!----------------
!  GFS containers
!----------------
type(GFS_externaldiag_type), target :: GFS_Diag(DIAG_SIZE)
type(GFS_restart_type)              :: GFS_restart_var

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

  real(kind=GFS_kind_phys), parameter :: zero    = 0.0_GFS_kind_phys,     &
                                         one     = 1.0_GFS_kind_phys,     &
                                         epsln   = 1.0e-10_GFS_kind_phys, &
                                         zorlmin = 1.0e-7_GFS_kind_phys

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
!-----------------------------------------------------------------------
  type (atmos_data_type), intent(in) :: Atmos
!--- local variables---
    integer :: nb, jdat(8), rc, ierr

    if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "statein driver"
!--- get atmospheric state from the dynamic core
    call set_atmosphere_pelist()
    call mpp_clock_begin(getClock)
    if (GFS_control%do_skeb) call atmosphere_diss_est (GFS_control%skeb_npass) !  do smoothing for SKEB
    call atmos_phys_driver_statein (GFS_data, Atm_block, flip_vc)
    call mpp_clock_end(getClock)

!--- if dycore only run, set up the dummy physics output state as the input state
    if (dycore_only) then
      do nb = 1,Atm_block%nblks
        GFS_data(nb)%Stateout%gu0 = GFS_data(nb)%Statein%ugrs
        GFS_data(nb)%Stateout%gv0 = GFS_data(nb)%Statein%vgrs
        GFS_data(nb)%Stateout%gt0 = GFS_data(nb)%Statein%tgrs
        GFS_data(nb)%Stateout%gq0 = GFS_data(nb)%Statein%qgrs
      enddo
    else
      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "setup step"

!--- update GFS_control%jdat(8)
      jdat(:) = 0
      call get_date (Atmos%Time, jdat(1), jdat(2), jdat(3),  &
                                 jdat(5), jdat(6), jdat(7))
      GFS_control%jdat(:) = jdat(:)

!--- execute the atmospheric setup step
      call mpp_clock_begin(setupClock)
      call CCPP_step (step="timestep_init", nblks=Atm_block%nblks, ierr=ierr)
      if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP timestep_init step failed')

!--- call stochastic physics pattern generation / cellular automata
      call stochastic_physics_wrapper(GFS_control, GFS_data, Atm_block, ierr)
      if (ierr/=0)  call mpp_error(FATAL, 'Call to stochastic_physics_wrapper failed')

!--- if coupled, assign coupled fields

      if (.not. GFS_control%cplchm) then
        call assign_importdata(rc)
      endif

      ! Calculate total non-physics tendencies by substracting old GFS Stateout
      ! variables from new/updated GFS Statein variables (gives the tendencies
      ! due to anything else than physics)
      if (GFS_control%ldiag3d) then
        do nb = 1,Atm_block%nblks
          GFS_data(nb)%Intdiag%du3dt(:,:,8)  = GFS_data(nb)%Intdiag%du3dt(:,:,8)  &
                                              + (GFS_data(nb)%Statein%ugrs - GFS_data(nb)%Stateout%gu0)
          GFS_data(nb)%Intdiag%dv3dt(:,:,8)  = GFS_data(nb)%Intdiag%dv3dt(:,:,8)  &
                                              + (GFS_data(nb)%Statein%vgrs - GFS_data(nb)%Stateout%gv0)
          GFS_data(nb)%Intdiag%dt3dt(:,:,11) = GFS_data(nb)%Intdiag%dt3dt(:,:,11) &
                                              + (GFS_data(nb)%Statein%tgrs - GFS_data(nb)%Stateout%gt0)
        enddo
        if (GFS_control%qdiag3d) then
          do nb = 1,Atm_block%nblks
            GFS_data(nb)%Intdiag%dq3dt(:,:,12) = GFS_data(nb)%Intdiag%dq3dt(:,:,12) &
                  + (GFS_data(nb)%Statein%qgrs(:,:,GFS_control%ntqv) - GFS_data(nb)%Stateout%gq0(:,:,GFS_control%ntqv))
            GFS_data(nb)%Intdiag%dq3dt(:,:,13) = GFS_data(nb)%Intdiag%dq3dt(:,:,13) &
                  + (GFS_data(nb)%Statein%qgrs(:,:,GFS_control%ntoz) - GFS_data(nb)%Stateout%gq0(:,:,GFS_control%ntoz))
          enddo
        endif
      endif

      call mpp_clock_end(setupClock)

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "radiation driver"

!--- execute the atmospheric radiation subcomponent (RRTM)

      call mpp_clock_begin(radClock)
      ! Performance improvement. Only enter if it is time to call the radiation physics.
      if (GFS_control%lsswr .or. GFS_control%lslwr) then
        call CCPP_step (step="radiation", nblks=Atm_block%nblks, ierr=ierr)
        if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP radiation step failed')
      endif
      call mpp_clock_end(radClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'RADIATION STEP  ', GFS_control%kdt, GFS_control%fhour
        call FV3GFS_GFS_checksum(GFS_control, GFS_data, Atm_block)
      endif

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "physics driver"

!--- execute the atmospheric physics step1 subcomponent (main physics driver)

      call mpp_clock_begin(physClock)
      call CCPP_step (step="physics", nblks=Atm_block%nblks, ierr=ierr)
      if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP physics step failed')
      call mpp_clock_end(physClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'PHYSICS STEP1   ', GFS_control%kdt, GFS_control%fhour
        call FV3GFS_GFS_checksum(GFS_control, GFS_data, Atm_block)
      endif

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "stochastic physics driver"

!--- execute the atmospheric physics step2 subcomponent (stochastic physics driver)

      call mpp_clock_begin(physClock)
      call CCPP_step (step="stochastics", nblks=Atm_block%nblks, ierr=ierr)
      if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP stochastics step failed')
      call mpp_clock_end(physClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'PHYSICS STEP2   ', GFS_control%kdt, GFS_control%fhour
        call FV3GFS_GFS_checksum(GFS_control, GFS_data, Atm_block)
      endif
      call getiauforcing(GFS_control,IAU_data)
      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "end of radiation and physics step"

!--- execute the atmospheric timestep finalize step
      call mpp_clock_begin(setupClock)
      call CCPP_step (step="timestep_finalize", nblks=Atm_block%nblks, ierr=ierr)
      if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP timestep_finalize step failed')
      call mpp_clock_end(setupClock)

    endif

    ! Update flag for first time step of time integration
    GFS_control%first_time_step = .false.

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

#ifdef _OPENMP
  use omp_lib
#endif
  use fv_mp_mod, only: commglobal
  use mpp_mod, only: mpp_npes
  use update_ca, only: read_ca_restart

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
  real(kind=GFS_kind_phys) :: dt_phys
  real, allocatable    :: q(:,:,:,:), p_half(:,:,:)
  character(len=80)    :: control
  character(len=64)    :: filename, filename2, pelist_name
  character(len=132)   :: text
  logical              :: p_hydro, hydro, fexist
  logical, save        :: block_message = .true.
  type(GFS_init_type)  :: Init_parm
  integer              :: bdat(8), cdat(8)
  integer              :: ntracers, maxhf, maxh
  character(len=32), allocatable, target :: tracer_names(:)
  integer,           allocatable, target :: tracer_types(:)
  integer :: nthrds, nb

!-----------------------------------------------------------------------

!---- set the atmospheric model time ------

   Atmos % Time_init = Time_init
   Atmos % Time      = Time
   Atmos % Time_step = Time_step
   call get_time (Atmos % Time_step, sec)
   dt_phys = real(sec)      ! integer seconds

   logunit = stdlog()

!---------- initialize atmospheric dynamics after reading the namelist -------
!---------- (need name of CCPP suite definition file from input.nml) ---------
   call atmosphere_init (Atmos%Time_init, Atmos%Time, Atmos%Time_step,&
                         Atmos%grid, Atmos%area)

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

!----------------------------------------------------------------------------------------------
! initialize atmospheric model - must happen AFTER atmosphere_init so that nests work correctly

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

!-----------------------------------------------------------------------
!--- before going any further check definitions for 'blocks'
!-----------------------------------------------------------------------
   call atmosphere_control_data (isc, iec, jsc, jec, nlev, p_hydro, hydro, tile_num)
   call define_blocks_packed ('atmos_model', Atm_block, isc, iec, jsc, jec, nlev, &
                              blocksize, block_message)

   allocate(DYCORE_Data(Atm_block%nblks))
   allocate(GFS_data(Atm_block%nblks))

#ifdef _OPENMP
   nthrds = omp_get_max_threads()
#else
   nthrds = 1
#endif

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
      allocate(GFS_interstitial(nthrds))
   else if (all(minloc(Atm_block%blksz)==(/size(Atm_block%blksz)/))) then
      non_uniform_blocks = .true.
      allocate(GFS_interstitial(nthrds+1))
   else
      call mpp_error(FATAL, 'For non-uniform blocksizes, only the last element ' // &
                            'in Atm_block%blksz can be different from the others')
   end if

!--- update GFS_control%jdat(8)
   bdat(:) = 0
   call get_date (Time_init, bdat(1), bdat(2), bdat(3),  &
                             bdat(5), bdat(6), bdat(7))
   cdat(:) = 0
   call get_date (Time,      cdat(1), cdat(2), cdat(3),  &
                             cdat(5), cdat(6), cdat(7))
   call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
   allocate (tracer_names(ntracers), tracer_types(ntracers))
   do i = 1, ntracers
     call get_tracer_names(MODEL_ATMOS, i, tracer_names(i))
   enddo
   call get_atmos_tracer_types(tracer_types)
!--- setup Init_parm
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
   Init_parm%nwat            = Atm(mygrid)%flagstruct%nwat
   Init_parm%tracer_names    => tracer_names
   Init_parm%tracer_types    => tracer_types
   Init_parm%restart         = Atm(mygrid)%flagstruct%warm_start
   Init_parm%hydrostatic     = Atm(mygrid)%flagstruct%hydrostatic

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

   call GFS_initialize (GFS_control, GFS_data%Statein, GFS_data%Stateout, GFS_data%Sfcprop,     &
                        GFS_data%Coupling, GFS_data%Grid, GFS_data%Tbd, GFS_data%Cldprop, GFS_data%Radtend, &
                        GFS_data%Intdiag, GFS_interstitial, commglobal, mpp_npes(), Init_parm)

   !--- populate/associate the Diag container elements
   call GFS_externaldiag_populate (GFS_Diag, GFS_Control, GFS_Data%Statein, GFS_Data%Stateout,   &
                                             GFS_Data%Sfcprop, GFS_Data%Coupling, GFS_Data%Grid, &
                                             GFS_Data%Tbd, GFS_Data%Cldprop, GFS_Data%Radtend,   &
                                             GFS_Data%Intdiag, Init_parm)

   Atmos%Diag => GFS_Diag

   Atm(mygrid)%flagstruct%do_skeb = GFS_control%do_skeb

!  initialize the IAU module
   call iau_initialize (GFS_control,IAU_data,Init_parm)

   Init_parm%blksz           => null()
   Init_parm%ak              => null()
   Init_parm%bk              => null()
   Init_parm%xlon            => null()
   Init_parm%xlat            => null()
   Init_parm%area            => null()
   Init_parm%tracer_names    => null()
   deallocate (tracer_names)
   deallocate (tracer_types)

   !--- update tracers in FV3 with any initialized during the physics/radiation init phase
!rab   call atmosphere_tracer_postinit (GFS_data, Atm_block)

   call atmosphere_nggps_diag (Time, init=.true.)
   call FV3GFS_diag_register (GFS_Diag, Time, Atm_block, GFS_control, Atmos%lon, Atmos%lat, Atmos%axes)
   call GFS_restart_populate (GFS_restart_var, GFS_control, GFS_data%Statein, GFS_data%Stateout, GFS_data%Sfcprop, &
                              GFS_data%Coupling, GFS_data%Grid, GFS_data%Tbd, GFS_data%Cldprop, GFS_data%Radtend, &
                              GFS_data%IntDiag, Init_parm, GFS_Diag)
   call FV3GFS_restart_read (GFS_data, GFS_restart_var, Atm_block, GFS_control, Atmos%domain, Atm(mygrid)%flagstruct%warm_start)
   if(GFS_control%ca_sgs)then
      call read_ca_restart (Atmos%domain,GFS_control%scells)
   endif
   ! Populate the GFS_data%Statein container with the prognostic state
   ! in Atm_block, which contains the initial conditions/restart data.
   call atmos_phys_driver_statein (GFS_data, Atm_block, flip_vc)

   ! When asked to calculate 3-dim. tendencies, set Stateout variables to
   ! Statein variables here in order to capture the first call to dycore
    if (GFS_control%ldiag3d) then
      do nb = 1,Atm_block%nblks
        GFS_data(nb)%Stateout%gu0 = GFS_data(nb)%Statein%ugrs
        GFS_data(nb)%Stateout%gv0 = GFS_data(nb)%Statein%vgrs
        GFS_data(nb)%Stateout%gt0 = GFS_data(nb)%Statein%tgrs
        GFS_data(nb)%Stateout%gq0 = GFS_data(nb)%Statein%qgrs
      enddo
    endif

   ! Initialize the CCPP framework
   call CCPP_step (step="init", nblks=Atm_block%nblks, ierr=ierr)
   if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP init step failed')
   ! Initialize the CCPP physics
   call CCPP_step (step="physics_init", nblks=Atm_block%nblks, ierr=ierr)
   if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP physics_init step failed')

!--- Initialize stochastic physics pattern generation / cellular automata for first time step
   call stochastic_physics_wrapper(GFS_control, GFS_data, Atm_block, ierr)
   if (ierr/=0)  call mpp_error(FATAL, 'Call to stochastic_physics_wrapper failed')

   !--- set the initial diagnostic timestamp
   diag_time = Time
   if (output_1st_tstep_rst) then
     diag_time = Time - real_to_time_type(mod(int((first_kdt - 1)*dt_phys/3600.),6)*3600.0)
   endif
   if (Atmos%iau_offset > zero) then
     call get_time (Atmos%Time - Atmos%Time_init, sec)
     if (sec < Atmos%iau_offset*3600) then
       diag_time = Atmos%Time_init
       diag_time_fhzero = Atmos%Time
     endif
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

!--- get bottom layer data from dynamical core for coupling
   call atmosphere_get_bottom_layer (Atm_block, DYCORE_Data)

   ! Set flag for first time step of time integration
   GFS_control%first_time_step = .true.

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
    if( GFS_control%cplchm ) then
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
    if( GFS_control%cplchm ) then
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
subroutine update_atmos_model_state (Atmos, rc)
! to update the model state after all concurrency is completed
  use ESMF
  type (atmos_data_type), intent(inout) :: Atmos
  integer, optional,      intent(out)   :: rc
!--- local variables
  integer :: localrc
  integer :: isec, seconds, isec_fhzero
  real(kind=GFS_kind_phys) :: time_int, time_intfull
!
    if (present(rc)) rc = ESMF_SUCCESS

    call set_atmosphere_pelist()
    call mpp_clock_begin(fv3Clock)
    call mpp_clock_begin(updClock)
    call atmosphere_state_update (Atmos%Time, GFS_data, IAU_Data, Atm_block, flip_vc)
    call mpp_clock_end(updClock)
    call mpp_clock_end(fv3Clock)

    if (chksum_debug) then
      if (mpp_pe() == mpp_root_pe()) print *,'UPDATE STATE    ', GFS_control%kdt, GFS_control%fhour
      if (mpp_pe() == mpp_root_pe()) print *,'in UPDATE STATE    ', size(GFS_data(1)%SfcProp%tsfc),'nblks=',Atm_block%nblks
      call FV3GFS_GFS_checksum(GFS_control, GFS_data, Atm_block)
    endif

    !--- advance time ---
    Atmos % Time = Atmos % Time + Atmos % Time_step

    call get_time (Atmos%Time - diag_time, isec)
    call get_time (Atmos%Time - Atmos%Time_init, seconds)
    call atmosphere_nggps_diag(Atmos%Time,ltavg=.true.,avg_max_length=avg_max_length)
    if (ANY(nint(fdiag(:)*3600.0) == seconds) .or. (GFS_control%kdt == first_kdt) .or. nsout > 0) then
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
      call FV3GFS_diag_output(Atmos%Time, GFS_Diag, Atm_block, GFS_control%nx, GFS_control%ny, &
                            GFS_control%levs, 1, 1, 1.0_GFS_kind_phys, time_int, time_intfull,              &
                            GFS_control%fhswr, GFS_control%fhlwr)
      if (nint(GFS_control%fhzero) > 0) then
        if (mod(isec,3600*nint(GFS_control%fhzero)) == 0) diag_time = Atmos%Time
      else
        if (mod(isec,nint(3600*GFS_control%fhzero)) == 0) diag_time = Atmos%Time
      endif
      call diag_send_complete_instant (Atmos%Time)
    endif

    !--- this may not be necessary once write_component is fully implemented
    !!!call diag_send_complete_extra (Atmos%Time)

    !--- get bottom layer data from dynamical core for coupling
    call atmosphere_get_bottom_layer (Atm_block, DYCORE_Data)

    !--- if in coupled mode, set up coupled fields
    call setup_exportdata(rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return

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
  use get_stochy_pattern_mod, only: write_stoch_restart_atm
  use update_ca, only: write_ca_restart
  type (atmos_data_type), intent(inout) :: Atmos
!---local variables
  integer :: idx, seconds, ierr

!-----------------------------------------------------------------------
!---- termination routine for atmospheric model ----

    call atmosphere_end (Atmos % Time, Atmos%grid, restart_endfcst)

    if(restart_endfcst) then
      call FV3GFS_restart_write (GFS_data, GFS_restart_var, Atm_block, &
                                 GFS_control, Atmos%domain)
      call write_stoch_restart_atm('RESTART/atm_stoch.res.nc')
      if(GFS_control%ca_sgs)then
         call write_ca_restart(Atmos%domain,GFS_control%scells)
      endif
    endif
    call stochastic_physics_wrapper_end(GFS_control)

!   Fast physics (from dynamics) are finalized in atmosphere_end above;
!   standard/slow physics (from CCPP) are finalized in CCPP_step 'finalize'.
!   The CCPP framework for all cdata structures is finalized in CCPP_step 'finalize'.
    call CCPP_step (step="finalize", nblks=Atm_block%nblks, ierr=ierr)
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP finalize step failed')

end subroutine atmos_model_end

! </SUBROUTINE>
!#######################################################################
! <SUBROUTINE NAME="atmos_model_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine atmos_model_restart(Atmos, timestamp)
  use update_ca, only: write_ca_restart
  type (atmos_data_type),   intent(inout) :: Atmos
  character(len=*),  intent(in)           :: timestamp

    call atmosphere_restart(timestamp)
    call FV3GFS_restart_write (GFS_data, GFS_restart_var, Atm_block, &
                               GFS_control, Atmos%domain, timestamp)
    if(GFS_control%ca_sgs)then
       call write_ca_restart(Atmos%domain,GFS_control%scells,timestamp)
    endif
end subroutine atmos_model_restart
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="get_atmos_model_ungridded_dim">
!
! <DESCRIPTION>
!  Retrieve ungridded dimensions of atmospheric model arrays
! </DESCRIPTION>

subroutine get_atmos_model_ungridded_dim(nlev, nsoillev, ntracers)

  integer, optional, intent(out) :: nlev, nsoillev, ntracers

  !--- number of atmospheric vertical levels
  if (present(nlev)) nlev = Atm_block%npz

  !--- number of soil levels
  if (present(nsoillev)) then
    nsoillev = 0
    if (allocated(GFS_data)) then
      if (associated(GFS_data(1)%Sfcprop%slc)) &
        nsoillev = size(GFS_data(1)%Sfcprop%slc, dim=2)
    end if
  end if

  !--- total number of atmospheric tracers
  if (present(ntracers)) call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)

end subroutine get_atmos_model_ungridded_dim
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="get_atmos_tracer_types">
! <DESCRIPTION>
!  Identify and return usage and type id of atmospheric tracers.
!  Ids are defined as:
!    0 = generic tracer
!    1 = chemistry - prognostic
!    2 = chemistry - diagnostic
!
!  Tracers are identified via the additional 'tracer_usage' keyword and
!  their optional 'type' qualifier. A tracer is assumed prognostic if
!  'type' is not provided. See examples from the field_table file below:
!
!  Prognostic tracer:
!  ------------------
!  "TRACER", "atmos_mod",    "so2"
!            "longname",     "so2 mixing ratio"
!            "units",        "ppm"
!            "tracer_usage", "chemistry"
!            "profile_type", "fixed", "surface_value=5.e-6" /
!
!  Diagnostic tracer:
!  ------------------
!  "TRACER", "atmos_mod",    "pm25"
!            "longname",     "PM2.5"
!            "units",        "ug/m3"
!            "tracer_usage", "chemistry", "type=diagnostic"
!            "profile_type", "fixed", "surface_value=5.e-6" /
!
!  For atmospheric chemistry, the order of both prognostic and diagnostic
!  tracers is validated against the model's internal assumptions.
!
! </DESCRIPTION>
subroutine get_atmos_tracer_types(tracer_types)

  use field_manager_mod,  only: parse
  use tracer_manager_mod, only: query_method

  integer, intent(out) :: tracer_types(:)

  !--- local variables
  logical :: found
  integer :: n, num_tracers, num_types
  integer :: id_max, id_min, id_num, ip_max, ip_min, ip_num
  character(len=32)  :: tracer_usage
  character(len=128) :: control, tracer_type

  !--- begin

  !--- validate array size
  call get_number_tracers(MODEL_ATMOS, num_tracers=num_tracers)

  if (size(tracer_types) < num_tracers) &
    call mpp_error(FATAL, 'insufficient size of tracer type array')

  !--- initialize tracer indices
  id_min = num_tracers + 1
  id_max = -id_min
  ip_min = id_min
  ip_max = id_max
  id_num = 0
  ip_num = 0

  do n = 1, num_tracers
    tracer_types(n) = 0
    found = query_method('tracer_usage',MODEL_ATMOS,n,tracer_usage,control)
    if (found) then
      if (trim(tracer_usage) == 'chemistry') then
        !--- set default to prognostic
        tracer_type = 'prognostic'
        num_types = parse(control, 'type', tracer_type)
        select case (trim(tracer_type))
          case ('diagnostic')
            tracer_types(n) = 2
            id_num = id_num + 1
            id_max = n
            if (id_num == 1) id_min = n
          case ('prognostic')
            tracer_types(n) = 1
            ip_num = ip_num + 1
            ip_max = n
            if (ip_num == 1) ip_min = n
        end select
      end if
    end if
  end do

  if (ip_num > 0) then
    !--- check if prognostic tracers are contiguous
    if (ip_num > ip_max - ip_min + 1) &
      call mpp_error(FATAL, 'prognostic chemistry tracers must be contiguous')
  end if

  if (id_num > 0) then
    !--- check if diagnostic tracers are contiguous
    if (id_num > id_max - id_min + 1) &
      call mpp_error(FATAL, 'diagnostic chemistry tracers must be contiguous')
  end if

  !--- prognostic tracers must precede diagnostic ones
  if (ip_max > id_min) &
    call mpp_error(FATAL, 'diagnostic chemistry tracers must follow prognostic ones')

end subroutine get_atmos_tracer_types
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
!         GFS_data
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
  integer :: nb, ix, i, j, k, k1, it
  integer :: ib, jb

  real(ESMF_KIND_R8), dimension(:,:,:),   pointer :: prsl, phil,   &
                                                     prsi, phii,   &
                                                     temp, cldfra, &
                                                     pflls, pfils, &
                                                     ua, va, slc
  real(ESMF_KIND_R8), dimension(:,:,:,:), pointer :: q

  real(ESMF_KIND_R8), dimension(:,:), pointer :: hpbl, area, rainc, &
    uustar, rain, slmsk, tsfc, shfsfc, zorl, focn, flake, fice,     &
    fsnow, u10m, v10m, swet

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

      !--- do not import tracer concentrations by default
      ntb = nt + 1
      nte = nt

      !--- if chemical tracers are present, set bounds appropriately
      if (GFS_control%ntchm > 0) then
        ntb = GFS_control%ntchs
        nte = GFS_control%ntche
      end if

      !--- prognostic tracer concentrations
      do it = ntb, nte
!$OMP parallel do default (none) &
!$OMP             shared  (it, nk, nj, ni, Atm_block, GFS_data, q)  &
!$OMP             private (k, j, jb, i, ib, nb, ix)
        do k = 1, nk
          do j = 1, nj
            jb = j + Atm_block%jsc - 1
            do i = 1, ni
              ib = i + Atm_block%isc - 1
              nb = Atm_block%blkno(ib,jb)
              ix = Atm_block%ixp(ib,jb)
              GFS_data(nb)%Stateout%gq0(ix,k,it) = q(i,j,k,it)
            enddo
          enddo
        enddo
      enddo

      !--- diagnostic tracers
      !--- set tracer concentrations in the atmospheric state directly
      !--- since the atmosphere's driver cannot perform this step while
      !--- updating the state
      if (GFS_control%ndchm > 0) then
        ntb = GFS_control%ndchs
        nte = GFS_control%ndche
!$OMP parallel do default (none) &
!$OMP             shared  (mygrid, nk, ntb, nte, Atm, Atm_block, q) &
!$OMP             private (i, ib, ix, j, jb, k, k1, nb)
        do nb = 1, Atm_block%nblks
          do k = 1, nk
            if(flip_vc) then
              k1 = nk+1-k !reverse the k direction
            else
              k1 = k
            endif
            do ix = 1, Atm_block%blksz(nb)
              ib = Atm_block%index(nb)%ii(ix)
              jb = Atm_block%index(nb)%jj(ix)
              i = ib - Atm_block%isc + 1
              j = jb - Atm_block%jsc + 1
              Atm(mygrid)%q(ib,jb,k1,ntb:nte) = q(i,j,k,ntb:nte)
            enddo
          end do
        end do
      end if

      if (GFS_control%debug) then
        write(6,'("update_atmos: ",a,": qgrs - min/max/avg",3g16.6)') &
          trim(state), minval(q), maxval(q), sum(q)/size(q)
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

      call cplFieldGet(state,'inst_tracer_mass_frac', farrayPtr4d=q, rc=localrc)
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

      call cplFieldGet(state,'inst_friction_velocity', farrayPtr2d=uustar, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_rainfall_amount', farrayPtr2d=rain, rc=localrc)
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

      call cplFieldGet(state,'inst_surface_roughness', farrayPtr2d=zorl, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_soil_moisture_content', farrayPtr3d=slc, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_liq_nonconv_tendency_levels', &
                       farrayPtr3d=pflls, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_ice_nonconv_tendency_levels', &
                       farrayPtr3d=pfils, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_cloud_frac_levels', farrayPtr3d=cldfra, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_zonal_wind_height10m', farrayPtr2d=u10m, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_merid_wind_height10m', farrayPtr2d=v10m, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_surface_soil_wetness', farrayPtr2d=swet, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'ice_fraction_in_atm', farrayPtr2d=fice, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'lake_fraction', farrayPtr2d=flake, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'ocean_fraction', farrayPtr2d=focn, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'surface_snow_area_fraction', farrayPtr2d=fsnow, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      !--- handle all three-dimensional variables
!$OMP parallel do default (none) &
!$OMP             shared  (nk, nj, ni, Atm_block, GFS_Data, GFS_Control, &
!$OMP                      cldfra, pfils, pflls, prsi, phii, prsl, phil, &
!$OMP                      temp, ua, va)  &
!$OMP             private (k, j, jb, i, ib, nb, ix)
      do k = 1, nk
        do j = 1, nj
          jb = j + Atm_block%jsc - 1
          do i = 1, ni
            ib = i + Atm_block%isc - 1
            nb = Atm_block%blkno(ib,jb)
            ix = Atm_block%ixp(ib,jb)
            !--- interface values
            prsi(i,j,k) = GFS_data(nb)%Statein%prsi(ix,k)
            phii(i,j,k) = GFS_data(nb)%Statein%phii(ix,k)
            !--- layer values
            prsl(i,j,k) = GFS_Data(nb)%Statein%prsl(ix,k)
            phil(i,j,k) = GFS_Data(nb)%Statein%phil(ix,k)
            temp(i,j,k) = GFS_Data(nb)%Stateout%gt0(ix,k)
            ua  (i,j,k) = GFS_Data(nb)%Stateout%gu0(ix,k)
            va  (i,j,k) = GFS_Data(nb)%Stateout%gv0(ix,k)
            cldfra(i,j,k) = GFS_Data(nb)%IntDiag%cldfra(ix,k)
            pfils (i,j,k) = GFS_Data(nb)%Coupling%pfi_lsan(ix,k)
            pflls (i,j,k) = GFS_Data(nb)%Coupling%pfl_lsan(ix,k)
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
          prsi(i,j,k) = GFS_data(nb)%Statein%prsi(ix,k)
          phii(i,j,k) = GFS_data(nb)%Statein%phii(ix,k)
        enddo
      enddo

      !--- tracers quantities
      do it = 1, nt
!$OMP parallel do default (none) &
!$OMP             shared  (it, nk, nj, ni, Atm_block, GFS_data, q)  &
!$OMP             private (k, j, jb, i, ib, nb, ix)
        do k = 1, nk
          do j = 1, nj
            jb = j + Atm_block%jsc - 1
            do i = 1, ni
              ib = i + Atm_block%isc - 1
              nb = Atm_block%blkno(ib,jb)
              ix = Atm_block%ixp(ib,jb)
              q(i,j,k,it) = GFS_data(nb)%Stateout%gq0(ix,k,it)
            enddo
          enddo
        enddo
      enddo

!$OMP parallel do default (none) &
!$OMP             shared  (nj, ni, Atm_block, GFS_data, GFS_Control, &
!$OMP                      hpbl, area, rainc, rain, uustar,      &
!$OMP                      fice, flake, focn, fsnow, u10m, v10m, &
!$OMP                      slmsk, tsfc, shfsfc, zorl, slc, swet) &
!$OMP             private (j, jb, i, ib, nb, ix)
      do j = 1, nj
        jb = j + Atm_block%jsc - 1
        do i = 1, ni
          ib = i + Atm_block%isc - 1
          nb = Atm_block%blkno(ib,jb)
          ix = Atm_block%ixp(ib,jb)
          hpbl(i,j)   = GFS_Data(nb)%Tbd%hpbl(ix)
          area(i,j)   = GFS_Data(nb)%Grid%area(ix)
          rainc(i,j)  = GFS_Data(nb)%Coupling%rainc_cpl(ix)
          rain(i,j)   = GFS_Data(nb)%Coupling%rain_cpl(ix)  &
                      + GFS_Data(nb)%Coupling%snow_cpl(ix)
          uustar(i,j) = GFS_Data(nb)%Sfcprop%uustar(ix)
          slmsk(i,j)  = GFS_Data(nb)%Sfcprop%slmsk(ix)
          shfsfc(i,j) = GFS_Data(nb)%Coupling%ushfsfci(ix)
          tsfc(i,j)   = GFS_Data(nb)%Coupling%tsfci_cpl(ix)
          zorl(i,j)   = GFS_Data(nb)%Sfcprop%zorl(ix)
          slc(i,j,:)  = GFS_Data(nb)%Sfcprop%slc(ix,:)
          u10m(i,j)   = GFS_Data(nb)%Coupling%u10mi_cpl(ix)
          v10m(i,j)   = GFS_Data(nb)%Coupling%v10mi_cpl(ix)
          focn(i,j)   = GFS_Data(nb)%Sfcprop%oceanfrac(ix)
          flake(i,j)  = max(zero, GFS_Data(nb)%Sfcprop%lakefrac(ix))
          fice(i,j)   = GFS_Data(nb)%Sfcprop%fice(ix)
          fsnow(i,j)  = GFS_Data(nb)%Sfcprop%sncovr(ix)
          if (GFS_Control%lsm == GFS_Control%lsm_ruc) then
            swet(i,j) = GFS_Data(nb)%Sfcprop%wetness(ix)
          else
            swet(i,j) = GFS_Data(nb)%IntDiag%wet1(ix)
          end if
        enddo
      enddo

      ! -- zero out accumulated fields
!$OMP parallel do default (none) &
!$OMP             shared  (nj, ni, Atm_block, GFS_control, GFS_data) &
!$OMP             private (j, jb, i, ib, nb, ix)
      do j = 1, nj
        jb = j + Atm_block%jsc - 1
        do i = 1, ni
          ib = i + Atm_block%isc - 1
          nb = Atm_block%blkno(ib,jb)
          ix = Atm_block%ixp(ib,jb)
          GFS_data(nb)%coupling%rainc_cpl(ix)  = zero
          if (.not.GFS_control%cplflx) then
            GFS_data(nb)%coupling%rain_cpl(ix) = zero
            GFS_data(nb)%coupling%snow_cpl(ix) = zero
          end if
        enddo
      enddo

      if (GFS_control%debug) then
        ! -- diagnostics
        write(6,'("update_atmos: prsi   - min/max/avg",3g16.6)') minval(prsi),   maxval(prsi),   sum(prsi)/size(prsi)
        write(6,'("update_atmos: phii   - min/max/avg",3g16.6)') minval(phii),   maxval(phii),   sum(phii)/size(phii)
        write(6,'("update_atmos: prsl   - min/max/avg",3g16.6)') minval(prsl),   maxval(prsl),   sum(prsl)/size(prsl)
        write(6,'("update_atmos: phil   - min/max/avg",3g16.6)') minval(phil),   maxval(phil),   sum(phil)/size(phil)
        write(6,'("update_atmos: tgrs   - min/max/avg",3g16.6)') minval(temp),   maxval(temp),   sum(temp)/size(temp)
        write(6,'("update_atmos: ugrs   - min/max/avg",3g16.6)') minval(ua),     maxval(ua),     sum(ua)/size(ua)
        write(6,'("update_atmos: vgrs   - min/max/avg",3g16.6)') minval(va),     maxval(va),     sum(va)/size(va)
        write(6,'("update_atmos: qgrs   - min/max/avg",3g16.6)') minval(q),      maxval(q),      sum(q)/size(q)

        write(6,'("update_atmos: hpbl   - min/max/avg",3g16.6)') minval(hpbl),   maxval(hpbl),   sum(hpbl)/size(hpbl)
        write(6,'("update_atmos: rainc  - min/max/avg",3g16.6)') minval(rainc),  maxval(rainc),  sum(rainc)/size(rainc)
        write(6,'("update_atmos: rain   - min/max/avg",3g16.6)') minval(rain),   maxval(rain),   sum(rain)/size(rain)
        write(6,'("update_atmos: shfsfc - min/max/avg",3g16.6)') minval(shfsfc), maxval(shfsfc), sum(shfsfc)/size(shfsfc)
        write(6,'("update_atmos: slmsk  - min/max/avg",3g16.6)') minval(slmsk),  maxval(slmsk),  sum(slmsk)/size(slmsk)
        write(6,'("update_atmos: tsfc   - min/max/avg",3g16.6)') minval(tsfc),   maxval(tsfc),   sum(tsfc)/size(tsfc)
        write(6,'("update_atmos: area   - min/max/avg",3g16.6)') minval(area),   maxval(area),   sum(area)/size(area)
        write(6,'("update_atmos: zorl   - min/max/avg",3g16.6)') minval(zorl),   maxval(zorl),   sum(zorl)/size(zorl)
        write(6,'("update_atmos: slc    - min/max/avg",3g16.6)') minval(slc),    maxval(slc),    sum(slc)/size(slc)
        write(6,'("update_atmos: cldfra - min/max/avg",3g16.6)') minval(cldfra), maxval(cldfra), sum(cldfra)/size(cldfra)
        write(6,'("update_atmos: fice   - min/max/avg",3g16.6)') minval(fice),   maxval(fice),   sum(fice)/size(fice)
        write(6,'("update_atmos: flake  - min/max/avg",3g16.6)') minval(flake),  maxval(flake),  sum(flake)/size(flake)
        write(6,'("update_atmos: focn   - min/max/avg",3g16.6)') minval(focn),   maxval(focn),   sum(focn)/size(focn)
        write(6,'("update_atmos: pfils  - min/max/avg",3g16.6)') minval(pfils),  maxval(pfils),  sum(pfils)/size(pfils)
        write(6,'("update_atmos: pflls  - min/max/avg",3g16.6)') minval(pflls),  maxval(pflls),  sum(pflls)/size(pflls)
        write(6,'("update_atmos: swet   - min/max/avg",3g16.6)') minval(swet),   maxval(swet),   sum(swet)/size(swet)
        write(6,'("update_atmos: u10m   - min/max/avg",3g16.6)') minval(u10m),   maxval(u10m),   sum(u10m)/size(u10m)
        write(6,'("update_atmos: v10m   - min/max/avg",3g16.6)') minval(v10m),   maxval(v10m),   sum(v10m)/size(v10m)
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

    use module_cplfields,  only: importFields, nImportFields, queryImportFields, &
                                 importFieldsValid
    use ESMF
!
    implicit none
    integer, intent(out) :: rc

    !--- local variables
    integer :: n, j, i, k, ix, nb, isc, iec, jsc, jec, nk, dimCount, findex
    integer :: sphum, liq_wat, ice_wat, o3mr
    character(len=128) :: impfield_name, fldname
    type(ESMF_TypeKind_Flag)                           :: datatype
    real(kind=ESMF_KIND_R4),  dimension(:,:), pointer  :: datar42d
    real(kind=ESMF_KIND_R8),  dimension(:,:), pointer  :: datar82d
    real(kind=ESMF_KIND_R8),  dimension(:,:,:), pointer:: datar83d
    real(kind=GFS_kind_phys), dimension(:,:), pointer  :: datar8
    real(kind=GFS_kind_phys)                           :: tem, ofrac
    logical found, isFieldCreated, lcpl_fice
    real (kind=GFS_kind_phys), parameter :: z0ice=1.1    !  (in cm)
!
!------------------------------------------------------------------------------
!
    rc  = -999

! set up local dimension
    isc = GFS_control%isc
    iec = GFS_control%isc+GFS_control%nx-1
    jsc = GFS_control%jsc
    jec = GFS_control%jsc+GFS_control%ny-1
    nk  = Atm_block%npz
    lcpl_fice = .false.

    allocate(datar8(isc:iec,jsc:jec))

!   if (mpp_pe() == mpp_root_pe() .and. debug) print *,'in cplImp,dim=',isc,iec,jsc,jec
!   if (mpp_pe() == mpp_root_pe() .and. debug) print *,'in cplImp,GFS_data, size', size(GFS_data)
!   if (mpp_pe() == mpp_root_pe() .and. debug) print *,'in cplImp,tsfc, size', size(GFS_data(1)%sfcprop%tsfc)
!   if (mpp_pe() == mpp_root_pe() .and. debug) print *,'in cplImp,tsfc, min_seaice', GFS_control%min_seaice

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

        else if( dimCount == 3) then
          if ( datatype == ESMF_TYPEKIND_R8) then
            call ESMF_FieldGet(importFields(n),farrayPtr=datar83d,localDE=0, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            found = .true.
          endif
        endif
!
        if (found) then
         if (datar8(isc,jsc) > -99998.0) then
!
        ! get sea land mask: in order to update the coupling fields over the ocean/ice
!        fldname = 'land_mask'
!        if (trim(impfield_name) == trim(fldname)) then
!          findex = queryImportFields(fldname)
!          if (importFieldsValid(findex)) then
!!$omp parallel do default(shared) private(i,j,nb,ix)
!            do j=jsc,jec
!              do i=isc,iec
!                nb = Atm_block%blkno(i,j)
!                ix = Atm_block%ixp(i,j)
!                GFS_data(nb)%Coupling%slimskin_cpl(ix) = datar8(i,j)
!              enddo
!            enddo
!            if( mpp_pe()==mpp_root_pe()) print *,'get land mask from mediator'
!          endif
!        endif


! get sea-state dependent surface roughness (if cplwav2atm=true)
!----------------------------
          fldname = 'wave_z0_roughness_length'
          if (trim(impfield_name) == trim(fldname)) then
            findex = queryImportFields(fldname)
            if (importFieldsValid(findex) .and. GFS_control%cplwav2atm) then
!$omp parallel do default(shared) private(i,j,nb,ix,tem)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (GFS_data(nb)%Sfcprop%oceanfrac(ix) > zero .and.  datar8(i,j) > zorlmin) then
                    tem = 100.0_GFS_kind_phys * min(0.1_GFS_kind_phys, datar8(i,j))
!                   GFS_data(nb)%Coupling%zorlwav_cpl(ix) = tem
                    GFS_data(nb)%Sfcprop%zorlw(ix)        = tem
                    GFS_data(nb)%Sfcprop%zorlwav(ix)      = tem
                  else
                    GFS_data(nb)%Sfcprop%zorlwav(ix) = -999.0_GFS_kind_phys

                  endif
                enddo
              enddo
            endif
          endif

! get sea ice surface temperature
!--------------------------------
          fldname = 'sea_ice_surface_temperature'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (GFS_data(nb)%Sfcprop%oceanfrac(ix) > zero .and.  datar8(i,j) > 150.0) then
!                   GFS_data(nb)%Coupling%tisfcin_cpl(ix) = datar8(i,j)
                    GFS_data(nb)%Sfcprop%tisfc(ix)       = datar8(i,j)
                  endif
                enddo
              enddo
            endif
          endif

! get sst:  sst needs to be adjusted by land sea mask before passing to fv3
!--------------------------------------------------------------------------
          fldname = 'sea_surface_temperature'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (GFS_data(nb)%Sfcprop%oceanfrac(ix) > zero .and. datar8(i,j) > 150.0) then
!                   GFS_data(nb)%Coupling%tseain_cpl(ix) = datar8(i,j)
                    GFS_data(nb)%Sfcprop%tsfco(ix)       = datar8(i,j)
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
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
              lcpl_fice = .true.
!$omp parallel do default(shared) private(i,j,nb,ix,ofrac)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)

                  GFS_data(nb)%Coupling%slimskin_cpl(ix) = GFS_data(nb)%Sfcprop%slmsk(ix)
                  ofrac = GFS_data(nb)%Sfcprop%oceanfrac(ix)
                  if (ofrac > zero) then
                    GFS_data(nb)%Sfcprop%fice(ix) = max(zero, min(one, datar8(i,j)/ofrac)) !LHS: ice frac wrt water area
                    if (GFS_data(nb)%Sfcprop%fice(ix) >= GFS_control%min_seaice) then
                      if (GFS_data(nb)%Sfcprop%fice(ix) > one-epsln) GFS_data(nb)%Sfcprop%fice(ix) = one
                      if (abs(one-ofrac) < epsln) GFS_data(nb)%Sfcprop%slmsk(ix) = 2.0_GFS_kind_phys !slmsk=2 crashes in gcycle on partial land points
!                     GFS_data(nb)%Sfcprop%slmsk(ix)         = 2.0_GFS_kind_phys
                      GFS_data(nb)%Coupling%slimskin_cpl(ix) = 4.0_GFS_kind_phys
                    else
                      GFS_data(nb)%Sfcprop%fice(ix) = zero
                      if (abs(one-ofrac) < epsln) then
                        GFS_data(nb)%Sfcprop%slmsk(ix)         = zero
                        GFS_data(nb)%Coupling%slimskin_cpl(ix) = zero
                      endif
                    endif
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get fice from mediator'
            endif
          endif

! get upward LW flux:  for sea ice covered area
!----------------------------------------------
          fldname = 'mean_up_lw_flx_ice'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
!               do i=isc,iec
!                 nb = Atm_block%blkno(i,j)
!                 ix = Atm_block%ixp(i,j)
!                if (GFS_data(nb)%Sfcprop%slmsk(ix) < 0.1 .or. GFS_data(nb)%Sfcprop%slmsk(ix) > 1.9) then
!                   GFS_data(nb)%Coupling%ulwsfcin_cpl(ix) = -datar8(i,j)
!                 endif
!               enddo
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (GFS_data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    GFS_data(nb)%Coupling%ulwsfcin_cpl(ix) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get lwflx from mediator'
            endif
          endif

! get latent heat flux:  for sea ice covered area
!------------------------------------------------
          fldname = 'mean_laten_heat_flx_atm_into_ice'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (GFS_data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    GFS_data(nb)%Coupling%dqsfcin_cpl(ix) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get laten_heat from mediator'
            endif
          endif

! get sensible heat flux:  for sea ice covered area
!--------------------------------------------------
          fldname = 'mean_sensi_heat_flx_atm_into_ice'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (GFS_data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    GFS_data(nb)%Coupling%dtsfcin_cpl(ix) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get sensi_heat from mediator'
            endif
          endif

! get zonal compt of momentum flux:  for sea ice covered area
!------------------------------------------------------------
          fldname = 'stress_on_air_ice_zonal'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (GFS_data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    GFS_data(nb)%Coupling%dusfcin_cpl(ix) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get zonal_moment_flx from mediator'
            endif
          endif

! get meridional compt of momentum flux:  for sea ice covered area
!-----------------------------------------------------------------
          fldname = 'stress_on_air_ice_merid'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (GFS_data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    GFS_data(nb)%Coupling%dvsfcin_cpl(ix) = -datar8(i,j)
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
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (GFS_data(nb)%Sfcprop%oceanfrac(ix) > zero) then
!                   GFS_data(nb)%Coupling%hicein_cpl(ix) = datar8(i,j)
                    GFS_data(nb)%Sfcprop%hice(ix)        = datar8(i,j)
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
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (GFS_data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    GFS_data(nb)%Coupling%hsnoin_cpl(ix) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get snow_volume from mediator'
            endif
          endif
!
! get instantaneous near IR albedo for diffuse radiation: for sea ice covered area
!---------------------------------------------------------------------------------
          fldname = 'inst_ice_ir_dif_albedo'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (GFS_data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    GFS_data(nb)%Coupling%sfc_alb_nir_dif_cpl(ix) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get sfc_alb_nir_dif_cpl from mediator'
            endif
          endif
!
! get instantaneous near IR albedo for direct radiation: for sea ice covered area
!---------------------------------------------------------------------------------
          fldname = 'inst_ice_ir_dir_albedo'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (GFS_data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    GFS_data(nb)%Coupling%sfc_alb_nir_dir_cpl(ix) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get sfc_alb_nir_dir_cpl from mediator'
            endif
          endif
!
! get instantaneous visible albedo for diffuse radiation: for sea ice covered area
!---------------------------------------------------------------------------------
          fldname = 'inst_ice_vis_dif_albedo'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (GFS_data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    GFS_data(nb)%Coupling%sfc_alb_vis_dif_cpl(ix) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get sfc_alb_vis_dif_cpl from mediator'
            endif
          endif

!
! get instantaneous visible IR albedo for direct radiation: for sea ice covered area
!---------------------------------------------------------------------------------
          fldname = 'inst_ice_vis_dir_albedo'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  if (GFS_data(nb)%Sfcprop%oceanfrac(ix) > zero) then
                    GFS_data(nb)%Coupling%sfc_alb_vis_dir_cpl(ix) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get inst_ice_vis_dir_albedo from mediator'
            endif
          endif


        endif ! if (datar8(isc,jsc) > -99999.0) then

!-------------------------------------------------------

       ! For JEDI

        sphum   = get_tracer_index(MODEL_ATMOS, 'sphum')
        liq_wat = get_tracer_index(MODEL_ATMOS, 'liq_wat')
        ice_wat = get_tracer_index(MODEL_ATMOS, 'ice_wat')
        o3mr    = get_tracer_index(MODEL_ATMOS, 'o3mr')

        fldname = 'u'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%u(i,j,k) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 'v'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%v(i,j,k) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 'ua'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%ua(i,j,k) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 'va'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%va(i,j,k) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 't'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%pt(i,j,k) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 'delp'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%delp(i,j,k) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 'sphum'
        if (trim(impfield_name) == trim(fldname) .and. sphum > 0) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%q(i,j,k,sphum) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 'ice_wat'
        if (trim(impfield_name) == trim(fldname) .and. ice_wat > 0) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%q(i,j,k,ice_wat) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 'liq_wat'
        if (trim(impfield_name) == trim(fldname) .and. liq_wat > 0) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%q(i,j,k,sphum) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 'o3mr'
        if (trim(impfield_name) == trim(fldname) .and. o3mr > 0) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%q(i,j,k,o3mr) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 'phis'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j)
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%phis(i,j) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'u_srf'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j)
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%u_srf(i,j) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'v_srf'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j)
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%v_srf(i,j) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        ! physics
        fldname = 'slmsk'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                GFS_data(nb)%Sfcprop%slmsk(ix) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'weasd'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                GFS_data(nb)%Sfcprop%weasd(ix) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'tsea'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                GFS_data(nb)%Sfcprop%tsfco(ix) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'vtype'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                GFS_data(nb)%Sfcprop%vtype(ix) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'stype'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                GFS_data(nb)%Sfcprop%stype(ix) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'vfrac'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                GFS_data(nb)%Sfcprop%vfrac(ix) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'stc'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                GFS_data(nb)%Sfcprop%stc(ix,:) = datar83d(i-isc+1,j-jsc+1,:)
              enddo
            enddo
          endif
        endif

        fldname = 'smc'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                GFS_data(nb)%Sfcprop%smc(ix,:) = datar83d(i-isc+1,j-jsc+1,:)
              enddo
            enddo
          endif
        endif

        fldname = 'snwdph'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                GFS_data(nb)%Sfcprop%snowd(ix) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'f10m'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                GFS_data(nb)%Sfcprop%f10m(ix) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'zorl'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                GFS_data(nb)%Sfcprop%zorl(ix) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 't2m'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                GFS_data(nb)%Sfcprop%t2m(ix) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        endif ! if (found) then
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
          if (GFS_data(nb)%Sfcprop%oceanfrac(ix) > zero) then
!if it is ocean or ice get surface temperature from mediator
            if (GFS_data(nb)%Sfcprop%fice(ix) >= GFS_control%min_seaice) then

!           if(GFS_data(nb)%Coupling%ficein_cpl(ix) >= GFS_control%min_seaice) then
!             GFS_data(nb)%Sfcprop%tisfc(ix)       = GFS_data(nb)%Coupling%tisfcin_cpl(ix)
!             GFS_data(nb)%Sfcprop%fice(ix)        = GFS_data(nb)%Coupling%ficein_cpl(ix)
!             GFS_data(nb)%Sfcprop%hice(ix)        = GFS_data(nb)%Coupling%hicein_cpl(ix)
!             GFS_data(nb)%Sfcprop%snowd(ix)       = GFS_data(nb)%Coupling%hsnoin_cpl(ix)

              GFS_data(nb)%Coupling%hsnoin_cpl(ix) = GFS_data(nb)%Coupling%hsnoin_cpl(ix) &
                                                   / max(0.01_GFS_kind_phys, GFS_data(nb)%Sfcprop%fice(ix))
!                                                  / max(0.01_GFS_kind_phys, GFS_data(nb)%Coupling%ficein_cpl(ix))
              GFS_data(nb)%Sfcprop%zorli(ix)       = z0ice
            else
!             GFS_data(nb)%Sfcprop%tisfc(ix)       = GFS_data(nb)%Coupling%tseain_cpl(ix)
              GFS_data(nb)%Sfcprop%tisfc(ix)       = GFS_data(nb)%Sfcprop%tsfco(ix)
              GFS_data(nb)%Sfcprop%fice(ix)        = zero
              GFS_data(nb)%Sfcprop%hice(ix)        = zero
!             GFS_data(nb)%Sfcprop%snowd(ix)       = zero
              GFS_data(nb)%Coupling%hsnoin_cpl(ix) = zero
!
              GFS_data(nb)%Coupling%dtsfcin_cpl(ix)  = -99999.0 ! over open water - should not be used in ATM
              GFS_data(nb)%Coupling%dqsfcin_cpl(ix)  = -99999.0 !                 ,,
              GFS_data(nb)%Coupling%dusfcin_cpl(ix)  = -99999.0 !                 ,,
              GFS_data(nb)%Coupling%dvsfcin_cpl(ix)  = -99999.0 !                 ,,
              GFS_data(nb)%Coupling%dtsfcin_cpl(ix)  = -99999.0 !                 ,,
              GFS_data(nb)%Coupling%ulwsfcin_cpl(ix) = -99999.0 !                 ,,
              if (abs(one-GFS_data(nb)%Sfcprop%oceanfrac(ix)) < epsln) then !  100% open water
                GFS_data(nb)%Coupling%slimskin_cpl(ix) = zero
                GFS_data(nb)%Sfcprop%slmsk(ix)         = zero
              endif
            endif
          endif
        enddo
      enddo
    endif
!
!-------------------------------------------------------------------------------
!   do j=jsc,jec
!     do i=isc,iec
!       nb = Atm_block%blkno(i,j)
!       ix = Atm_block%ixp(i,j)
!       if (abs(GFS_data(nb)%Grid%xlon_d(ix)-2.89) < 0.1 .and. &
!           abs(GFS_data(nb)%Grid%xlat_d(ix)+58.99) < 0.1) then
!         write(0,*)' in assign tisfc=',GFS_data(nb)%Sfcprop%tisfc(ix),     &
!          ' oceanfrac=',GFS_data(nb)%Sfcprop%oceanfrac(ix),' i=',i,' j=',j,&
!!         ' tisfcin=',GFS_data(nb)%Coupling%tisfcin_cpl(ix),               &
!          ' tisfcin=',GFS_data(nb)%Sfcprop%tisfc(ix),                      &
!          ' fice=',GFS_data(nb)%Sfcprop%fice(ix)
!       endif
!     enddo
!   enddo
!-------------------------------------------------------------------------------
!

    rc=0
!
    if (mpp_pe() == mpp_root_pe()) print *,'end of assign_importdata'
  end subroutine assign_importdata

!
  subroutine setup_exportdata(rc)

    use ESMF

    use module_cplfields, only: exportFields

    !--- arguments
    integer, optional, intent(out) :: rc

    !--- local variables
    integer                :: i, j, k, idx, ix
    integer                :: isc, iec, jsc, jec
    integer                :: ib, jb, nb, nsb, nk
    integer                :: sphum, liq_wat, ice_wat, o3mr
    real(GFS_kind_phys)    :: rtime, rtimek

    integer                                     :: localrc
    integer                                     :: n,rank
    logical                                     :: isFound
    type(ESMF_TypeKind_Flag)                    :: datatype
    character(len=ESMF_MAXSTR)                  :: fieldName
    real(kind=ESMF_KIND_R4), dimension(:,:), pointer   :: datar42d
    real(kind=ESMF_KIND_R8), dimension(:,:), pointer   :: datar82d
    real(kind=ESMF_KIND_R8), dimension(:,:,:), pointer :: datar83d

    !--- local parameters
    real(kind=ESMF_KIND_R8), parameter :: zeror8 = 0._ESMF_KIND_R8

    !--- begin
    if (present(rc)) rc = ESMF_SUCCESS

    !--- disable if coupling with chemistry
    if (GFS_control%cplchm) return

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    nk  = Atm_block%npz
    nsb = Atm_block%blkno(isc,jsc)

    rtime  = one / GFS_control%dtp
    rtimek = GFS_control%rho_h2o * rtime

    do n=1, size(exportFields)

      datar42d => null()
      datar82d => null()
      datar83d => null()

      isFound = ESMF_FieldIsCreated(exportFields(n), rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      if (isFound) then
        call ESMF_FieldGet(exportFields(n), name=fieldname, rank=rank, typekind=datatype, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
        if (datatype == ESMF_TYPEKIND_R8) then
           select case (rank)
             case (2)
               call ESMF_FieldGet(exportFields(n),farrayPtr=datar82d,localDE=0, rc=localrc)
               if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
             case (3)
               call ESMF_FieldGet(exportFields(n),farrayPtr=datar83d,localDE=0, rc=localrc)
               if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
             case default
               !--- skip field
               isFound = .false.
           end select
        else if (datatype == ESMF_TYPEKIND_R4) then
           select case (rank)
             case (2)
               call ESMF_FieldGet(exportFields(n),farrayPtr=datar42d,localDE=0, rc=localrc)
               if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
             case default
               !--- skip field
               isFound = .false.
           end select
        else
          !--- skip field
          isFound = .false.
        end if
      end if

      if (isFound) then
!$omp parallel do default(shared) private(nb) reduction(max:localrc)
        do nb = 1, Atm_block%nblks
          select case (trim(fieldname))
            !--- Instantaneous quantities
            ! Instantaneous u wind (m/s) 10 m above ground
            case ('inst_zonal_wind_height10m')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%u10mi_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous v wind (m/s) 10 m above ground
            case ('inst_merid_wind_height10m')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%v10mi_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous Zonal compt of momentum flux (N/m**2)
            case ('inst_zonal_moment_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dusfci_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous Merid compt of momentum flux (N/m**2)
            case ('inst_merid_moment_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dvsfci_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous Sensible heat flux (W/m**2)
            case ('inst_sensi_heat_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dtsfci_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous Latent heat flux (W/m**2)
            case ('inst_laten_heat_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dqsfci_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous Downward long wave radiation flux (W/m**2)
            case ('inst_down_lw_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dlwsfci_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous Downward solar radiation flux (W/m**2)
            case ('inst_down_sw_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dswsfci_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous Temperature (K) 2 m above ground
            case ('inst_temp_height2m')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%t2mi_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous Specific humidity (kg/kg) 2 m above ground
            case ('inst_spec_humid_height2m')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%q2mi_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous Temperature (K) at surface
            case ('inst_temp_height_surface')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%tsfci_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous Pressure (Pa) land and sea surface
            case ('inst_pres_height_surface')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%psurfi_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous Surface height (m)
            case ('inst_surface_height')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%oro_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous NET long wave radiation flux (W/m**2)
            case ('inst_net_lw_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%nlwsfci_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous NET solar radiation flux over the ocean (W/m**2)
            case ('inst_net_sw_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%nswsfci_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous sfc downward nir direct flux (W/m**2)
            case ('inst_down_sw_ir_dir_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dnirbmi_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous sfc downward nir diffused flux (W/m**2)
            case ('inst_down_sw_ir_dif_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dnirdfi_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous sfc downward uv+vis direct flux (W/m**2)
            case ('inst_down_sw_vis_dir_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dvisbmi_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous sfc downward uv+vis diffused flux (W/m**2)
            case ('inst_down_sw_vis_dif_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dvisdfi_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous net sfc nir direct flux (W/m**2)
            case ('inst_net_sw_ir_dir_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%nnirbmi_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous net sfc nir diffused flux (W/m**2)
            case ('inst_net_sw_ir_dif_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%nnirdfi_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous net sfc uv+vis direct flux (W/m**2)
            case ('inst_net_sw_vis_dir_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%nvisbmi_cpl, Atm_block, nb, rc=localrc)
            ! Instantaneous net sfc uv+vis diffused flux (W/m**2)
            case ('inst_net_sw_vis_dif_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%nvisdfi_cpl, Atm_block, nb, rc=localrc)
            ! Land/Sea mask (sea:0,land:1)
            case ('inst_land_sea_mask', 'slmsk')
              call block_data_copy(datar82d, GFS_data(nb)%sfcprop%slmsk, Atm_block, nb, rc=localrc)
            !--- Mean quantities
            ! MEAN Zonal compt of momentum flux (N/m**2)
            case ('mean_zonal_moment_flx_atm')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dusfc_cpl, Atm_block, nb, scale_factor=rtime, rc=localrc)
            ! MEAN Merid compt of momentum flux (N/m**2)
            case ('mean_merid_moment_flx_atm')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dvsfc_cpl, Atm_block, nb, scale_factor=rtime, rc=localrc)
            ! MEAN Sensible heat flux (W/m**2)
            case ('mean_sensi_heat_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dtsfc_cpl, Atm_block, nb, scale_factor=rtime, rc=localrc)
            ! MEAN Latent heat flux (W/m**2)
            case ('mean_laten_heat_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dqsfc_cpl, Atm_block, nb, scale_factor=rtime, rc=localrc)
            ! MEAN Downward LW heat flux (W/m**2)
            case ('mean_down_lw_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dlwsfc_cpl, Atm_block, nb, scale_factor=rtime, rc=localrc)
            ! MEAN Downward SW heat flux (W/m**2)
            case ('mean_down_sw_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dswsfc_cpl, Atm_block, nb, scale_factor=rtime, rc=localrc)
            ! MEAN NET long wave radiation flux (W/m**2)
            case ('mean_net_lw_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%nlwsfc_cpl, Atm_block, nb, scale_factor=rtime, rc=localrc)
            ! MEAN NET solar radiation flux over the ocean (W/m**2)
            case ('mean_net_sw_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%nswsfc_cpl, Atm_block, nb, scale_factor=rtime, rc=localrc)
            ! MEAN sfc downward nir direct flux (W/m**2)
            case ('mean_down_sw_ir_dir_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dnirbm_cpl, Atm_block, nb, scale_factor=rtime, rc=localrc)
            ! MEAN sfc downward nir diffused flux (W/m**2)
            case ('mean_down_sw_ir_dif_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dnirdf_cpl, Atm_block, nb, scale_factor=rtime, rc=localrc)
            ! MEAN sfc downward uv+vis direct flux (W/m**2)
            case ('mean_down_sw_vis_dir_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dvisbm_cpl, Atm_block, nb, scale_factor=rtime, rc=localrc)
            ! MEAN sfc downward uv+vis diffused flux (W/m**2)
            case ('mean_down_sw_vis_dif_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%dvisdf_cpl, Atm_block, nb, scale_factor=rtime, rc=localrc)
            ! MEAN NET sfc nir direct flux (W/m**2)
            case ('mean_net_sw_ir_dir_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%nnirbm_cpl, Atm_block, nb, scale_factor=rtime, rc=localrc)
            ! MEAN NET sfc nir diffused flux (W/m**2)
            case ('mean_net_sw_ir_dif_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%nnirdf_cpl, Atm_block, nb, scale_factor=rtime, rc=localrc)
            ! MEAN NET sfc uv+vis direct flux (W/m**2)
            case ('mean_net_sw_vis_dir_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%nvisbm_cpl, Atm_block, nb, scale_factor=rtime, rc=localrc)
            ! MEAN NET sfc uv+vis diffused flux (W/m**2)
            case ('mean_net_sw_vis_dif_flx')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%nvisdf_cpl, Atm_block, nb, scale_factor=rtime, rc=localrc)
            ! MEAN precipitation rate (kg/m2/s)
            case ('mean_prec_rate')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%rain_cpl, Atm_block, nb, scale_factor=rtimek, rc=localrc)
            ! MEAN snow precipitation rate (kg/m2/s)
            case ('mean_fprec_rate')
              call block_data_copy(datar82d, GFS_data(nb)%coupling%snow_cpl, Atm_block, nb, scale_factor=rtimek, rc=localrc)
            ! oceanfrac used by atm to calculate fluxes
            case ('openwater_frac_in_atm')
              call block_data_combine_fractions(datar82d, GFS_data(nb)%sfcprop%oceanfrac, GFS_Data(nb)%sfcprop%fice, Atm_block, nb, rc=localrc)
            !--- Dycore quantities
            ! bottom layer temperature (t)
            case('inst_temp_height_lowest')
              call block_data_copy_or_fill(datar82d, DYCORE_data(nb)%coupling%t_bot, zeror8, Atm_block, nb, rc=localrc)
            ! bottom layer specific humidity (q)
            !    !    ! CHECK if tracer 1 is for specific humidity     !    !    !
            case('inst_spec_humid_height_lowest')
              call block_data_copy_or_fill(datar82d, DYCORE_data(nb)%coupling%tr_bot, 1, zeror8, Atm_block, nb, rc=localrc)
            ! bottom layer zonal wind (u)
            case('inst_zonal_wind_height_lowest')
              call block_data_copy_or_fill(datar82d, DYCORE_data(nb)%coupling%u_bot, zeror8, Atm_block, nb, rc=localrc)
            ! bottom layer meridionalw wind (v)
            case('inst_merid_wind_height_lowest')
              call block_data_copy_or_fill(datar82d, DYCORE_data(nb)%coupling%v_bot, zeror8, Atm_block, nb, rc=localrc)
            ! bottom layer pressure (p)
            case('inst_pres_height_lowest')
              call block_data_copy_or_fill(datar82d, DYCORE_data(nb)%coupling%p_bot, zeror8, Atm_block, nb, rc=localrc)
            ! bottom layer height (z)
            case('inst_height_lowest')
              call block_data_copy_or_fill(datar82d, DYCORE_data(nb)%coupling%z_bot, zeror8, Atm_block, nb, rc=localrc)
            !--- JEDI fields
            case ('u')
              call block_atmos_copy(datar83d, Atm(mygrid)%u, Atm_block, nb, rc=localrc)
            case ('v')
              call block_atmos_copy(datar83d, Atm(mygrid)%v, Atm_block, nb, rc=localrc)
            case ('ua')
              call block_atmos_copy(datar83d, Atm(mygrid)%ua, Atm_block, nb, rc=localrc)
            case ('va')
              call block_atmos_copy(datar83d, Atm(mygrid)%va, Atm_block, nb, rc=localrc)
            case ('t')
              call block_atmos_copy(datar83d, Atm(mygrid)%pt, Atm_block, nb, rc=localrc)
            case ('delp')
              call block_atmos_copy(datar83d, Atm(mygrid)%delp, Atm_block, nb, rc=localrc)
            case ('sphum')
              sphum = get_tracer_index(MODEL_ATMOS, 'sphum')
              call block_atmos_copy(datar83d, Atm(mygrid)%q, sphum, Atm_block, nb, rc=localrc)
            case ('ice_wat')
              ice_wat = get_tracer_index(MODEL_ATMOS, 'ice_wat')
              call block_atmos_copy(datar83d, Atm(mygrid)%q, ice_wat, Atm_block, nb, rc=localrc)
            case ('liq_wat')
              liq_wat = get_tracer_index(MODEL_ATMOS, 'liq_wat')
              call block_atmos_copy(datar83d, Atm(mygrid)%q, liq_wat, Atm_block, nb, rc=localrc)
            case ('o3mr')
              o3mr = get_tracer_index(MODEL_ATMOS, 'o3mr')
              call block_atmos_copy(datar83d, Atm(mygrid)%q, o3mr, Atm_block, nb, rc=localrc)
            case ('phis')
              call block_atmos_copy(datar82d, Atm(mygrid)%phis, Atm_block, nb, rc=localrc)
            case ('u_srf')
              call block_atmos_copy(datar82d, Atm(mygrid)%u_srf, Atm_block, nb, rc=localrc)
            case ('v_srf')
              call block_atmos_copy(datar82d, Atm(mygrid)%v_srf, Atm_block, nb, rc=localrc)
            case ('weasd')
              call block_data_copy(datar82d, GFS_data(nb)%sfcprop%weasd, Atm_block, nb, rc=localrc)
            case ('tsea')
              call block_data_copy(datar82d, GFS_data(nb)%sfcprop%tsfco, Atm_block, nb, rc=localrc)
            case ('vtype')
              call block_data_copy(datar82d, GFS_data(nb)%sfcprop%vtype, Atm_block, nb, rc=localrc)
            case ('stype')
              call block_data_copy(datar82d, GFS_data(nb)%sfcprop%stype, Atm_block, nb, rc=localrc)
            case ('vfrac')
              call block_data_copy(datar82d, GFS_data(nb)%sfcprop%vfrac, Atm_block, nb, rc=localrc)
            case ('stc')
              call block_data_copy(datar83d, GFS_data(nb)%sfcprop%stc, Atm_block, nb, rc=localrc)
            case ('smc')
              call block_data_copy(datar83d, GFS_data(nb)%sfcprop%smc, Atm_block, nb, rc=localrc)
            case ('snwdph')
              call block_data_copy(datar82d, GFS_data(nb)%sfcprop%snowd, Atm_block, nb, rc=localrc)
            case ('f10m')
              call block_data_copy(datar82d, GFS_data(nb)%sfcprop%f10m, Atm_block, nb, rc=localrc)
            case ('zorl')
              call block_data_copy(datar82d, GFS_data(nb)%sfcprop%zorl, Atm_block, nb, rc=localrc)
            case ('t2m')
              call block_data_copy(datar82d, GFS_data(nb)%sfcprop%t2m, Atm_block, nb, rc=localrc)
            case default
              localrc = ESMF_RC_NOT_FOUND
          end select
        enddo
        if (ESMF_LogFoundError(rcToCheck=localrc, msg="Failure to populate exported field: "//trim(fieldname), &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      endif
    enddo ! exportFields

!---
    if (GFS_control%cplflx) then
! zero out accumulated fields
!$omp parallel do default(shared) private(i,j,nb,ix)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          GFS_data(nb)%coupling%dusfc_cpl(ix)  = zero
          GFS_data(nb)%coupling%dvsfc_cpl(ix)  = zero
          GFS_data(nb)%coupling%dtsfc_cpl(ix)  = zero
          GFS_data(nb)%coupling%dqsfc_cpl(ix)  = zero
          GFS_data(nb)%coupling%dlwsfc_cpl(ix) = zero
          GFS_data(nb)%coupling%dswsfc_cpl(ix) = zero
          GFS_data(nb)%coupling%rain_cpl(ix)   = zero
          GFS_data(nb)%coupling%nlwsfc_cpl(ix) = zero
          GFS_data(nb)%coupling%nswsfc_cpl(ix) = zero
          GFS_data(nb)%coupling%dnirbm_cpl(ix) = zero
          GFS_data(nb)%coupling%dnirdf_cpl(ix) = zero
          GFS_data(nb)%coupling%dvisbm_cpl(ix) = zero
          GFS_data(nb)%coupling%dvisdf_cpl(ix) = zero
          GFS_data(nb)%coupling%nnirbm_cpl(ix) = zero
          GFS_data(nb)%coupling%nnirdf_cpl(ix) = zero
          GFS_data(nb)%coupling%nvisbm_cpl(ix) = zero
          GFS_data(nb)%coupling%nvisdf_cpl(ix) = zero
          GFS_data(nb)%coupling%snow_cpl(ix)   = zero
        enddo
      enddo
      if (mpp_pe() == mpp_root_pe()) print *,'zeroing coupling accumulated fields at kdt= ',GFS_control%kdt
    endif !cplflx

  end subroutine setup_exportdata

  subroutine addLsmask2grid(fcstGrid, rc)

    use ESMF
!
    implicit none
    type(ESMF_Grid)      :: fcstGrid
    integer, optional, intent(out) :: rc
!
!  local vars
    integer isc, iec, jsc, jec
    integer i, j, nb, ix
!    integer CLbnd(2), CUbnd(2), CCount(2), TLbnd(2), TUbnd(2), TCount(2)
    type(ESMF_StaggerLoc) :: staggerloc
    integer, allocatable  :: lsmask(:,:)
    integer(kind=ESMF_KIND_I4), pointer  :: maskPtr(:,:)
!
    isc = GFS_control%isc
    iec = GFS_control%isc+GFS_control%nx-1
    jsc = GFS_control%jsc
    jec = GFS_control%jsc+GFS_control%ny-1
    allocate(lsmask(isc:iec,jsc:jec))
!
!$omp parallel do default(shared) private(i,j,nb,ix)
    do j=jsc,jec
      do i=isc,iec
        nb = Atm_block%blkno(i,j)
        ix = Atm_block%ixp(i,j)
! use land sea mask: land:1, ocean:0
        lsmask(i,j) = floor(one + epsln - GFS_data(nb)%SfcProp%oceanfrac(ix))
      enddo
    enddo
!
! Get mask
    call ESMF_GridAddItem(fcstGrid, itemflag=ESMF_GRIDITEM_MASK,   &
                          staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!    call ESMF_GridGetItemBounds(fcstGrid, itemflag=ESMF_GRIDITEM_MASK,   &
!         staggerloc=ESMF_STAGGERLOC_CENTER, computationalLBound=ClBnd,  &
!         computationalUBound=CUbnd, computationalCount=Ccount,  &
!         totalLBound=TLbnd, totalUBound=TUbnd, totalCount=Tcount, rc=rc)
!    print *,'in set up grid, aft add esmfgridadd item mask, rc=',rc, &
!     'ClBnd=',ClBnd,'CUbnd=',CUbnd,'Ccount=',Ccount, &
!     'TlBnd=',TlBnd,'TUbnd=',TUbnd,'Tcount=',Tcount
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridGetItem(fcstGrid, itemflag=ESMF_GRIDITEM_MASK,   &
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
