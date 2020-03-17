#undef MULTI_GASES

module GFS_typedefs

       use machine,                  only: kind_phys
#ifdef CCPP
       use physcons,                 only: con_cp, con_fvirt, con_g,                       &
                                           con_hvap, con_hfus, con_pi, con_rd, con_rv,     &
                                           con_t0c, con_cvap, con_cliq, con_eps, con_epsq, &
                                           con_epsm1, con_ttp, rlapse, con_jcal, con_rhw0, &
                                           con_sbc, con_tice, cimin, con_p0, rhowater
       use module_radsw_parameters,  only: topfsw_type, sfcfsw_type, cmpfsw_type, NBDSW
       use module_radlw_parameters,  only: topflw_type, sfcflw_type, NBDLW
#else
       use physcons,                 only: rhowater
       use module_radsw_parameters,  only: topfsw_type, sfcfsw_type
       use module_radlw_parameters,  only: topflw_type, sfcflw_type
       use ozne_def,                 only: levozp,      oz_coeff
       use h2o_def,                  only: levh2o,      h2o_coeff
       use aerclm_def,               only: ntrcaer,     ntrcaerm
#endif

       implicit none

#ifdef CCPP
      ! To ensure that these values match what's in the physics,
      ! array sizes are compared during model init in GFS_rrtmg_setup_init()
      private :: NF_AESW, NF_AELW, NSPC, NSPC1, NF_CLDS, NF_VGAS, NF_ALBD, ntrcaerm
      ! from module_radiation_aerosols
      integer, parameter :: NF_AESW = 3
      integer, parameter :: NF_AELW = 3
      integer, parameter :: NSPC    = 5
      integer, parameter :: NSPC1   = NSPC + 1
      ! from module_radiation_clouds
      integer, parameter :: NF_CLDS = 9
      ! from module_radiation_gases
      integer, parameter :: NF_VGAS = 10
      ! from module_radiation_surface
      integer, parameter :: NF_ALBD = 4
      ! from aerclm_def
      integer, parameter :: ntrcaerm = 15

      ! These will be set later in GFS_Control%initialize,
      ! since they depend on the runtime config (e.g. Model%ntoz, Model%h2o_phys, Model%aero_in)
      private :: levozp, oz_coeff, levh2o, h2o_coeff, ntrcaer
      integer :: levozp, oz_coeff, levh2o, h2o_coeff, ntrcaer
#endif

!> \section arg_table_GFS_typedefs
!! \htmlinclude GFS_typedefs.html
!!

       !--- version of physics
       character(len=64) :: phys_version = 'v2018 FV3GFS BETA VERSION PHYSICS'

       !--- parameter constants used for default initializations
       real(kind=kind_phys), parameter :: zero      = 0.0_kind_phys
       real(kind=kind_phys), parameter :: huge      = 9.9692099683868690E36 ! NetCDF float FillValue
       real(kind=kind_phys), parameter :: clear_val = zero
      !real(kind=kind_phys), parameter :: clear_val = -9.9999e80
       real(kind=kind_phys), parameter :: rann_init = 0.6_kind_phys
       real(kind=kind_phys), parameter :: cn_one    = 1._kind_phys
       real(kind=kind_phys), parameter :: cn_100    = 100._kind_phys
       real(kind=kind_phys), parameter :: cn_th     = 1000._kind_phys
       real(kind=kind_phys), parameter :: cn_hr     = 3600._kind_phys
#ifdef CCPP
       ! optional extra top layer on top of low ceiling models
       ! this parameter was originally defined in the radiation driver
       ! (and is still for standard non-CCPP builds), but is required
       ! here for CCPP to allocate arrays used for the interstitial
       ! calculations previously in GFS_{physics,radiation}_driver.F90
       ! LTP=0: no extra top layer
       integer, parameter :: LTP = 0   ! no extra top layer
       !integer, parameter :: LTP = 1   ! add an extra top layer
#endif

!----------------
! Data Containers
!----------------
!    !--- GFS external initialization type
!    GFS_init_type
!    !--- GFS Derived Data Types (DDTs)
!    GFS_statein_type        !< prognostic state data in from dycore
!    GFS_stateout_type       !< prognostic state or tendencies return to dycore
!    GFS_sfcprop_type        !< surface fields
!    GFS_coupling_type       !< fields to/from coupling with other components (e.g. land/ice/ocean/etc.)
!    !---GFS specific containers
!    GFS_control_type        !< model control parameters 
!    GFS_grid_type           !< grid and interpolation related data
!    GFS_tbd_type            !< to be determined data that doesn't fit in any one container
!    GFS_cldprop_type        !< cloud fields needed by radiation from physics
!    GFS_radtend_type        !< radiation tendencies needed in physics
!    GFS_diag_type           !< fields targetted for diagnostic output
#ifdef CCPP
!    GFS_interstitial_type   !< fields required to replace interstitial code in GFS_{physics,radiation}_driver.F90 in CCPP
!    GFS_data_type           !< combined type of all of the above except GFS_control_type and GFS_interstitial_type
#endif

!--------------------------------------------------------------------------------
! GFS_init_type
!--------------------------------------------------------------------------------
!   This container is the minimum set of data required from the dycore/atmosphere
!   component to allow proper initialization of the GFS physics
!--------------------------------------------------------------------------------
!! \section arg_table_GFS_init_type
!! \htmlinclude GFS_init_type.html
!!
  type GFS_init_type
    integer :: me                                !< my MPI-rank
    integer :: master                            !< master MPI-rank
    integer :: tile_num                          !< tile number for this MPI rank
    integer :: isc                               !< starting i-index for this MPI-domain
    integer :: jsc                               !< starting j-index for this MPI-domain
    integer :: nx                                !< number of points in i-dir for this MPI rank
    integer :: ny                                !< number of points in j-dir for this MPI rank
    integer :: levs                              !< number of vertical levels
    integer :: cnx                               !< number of points in i-dir for this cubed-sphere face
                                                 !< equal to gnx for lat-lon grids
    integer :: cny                               !< number of points in j-dir for this cubed-sphere face
                                                 !< equal to gny for lat-lon grids
    integer :: gnx                               !< number of global points in x-dir (i) along the equator
    integer :: gny                               !< number of global points in y-dir (j) along any meridian
    integer :: nlunit                            !< fortran unit number for file opens
    integer :: logunit                           !< fortran unit number for writing logfile
    integer :: bdat(8)                           !< model begin date in GFS format   (same as idat)
    integer :: cdat(8)                           !< model current date in GFS format (same as jdat)
    integer :: iau_offset                        !< iau running window length
    real(kind=kind_phys) :: dt_dycore            !< dynamics time step in seconds
    real(kind=kind_phys) :: dt_phys              !< physics  time step in seconds
#ifdef CCPP
!--- restart information
    logical :: restart                           !< flag whether this is a coldstart (.false.) or a warmstart/restart (.true.)
!--- hydrostatic/non-hydrostatic flag
    logical :: hydrostatic                       !< flag whether this is a hydrostatic or non-hydrostatic run
#endif
!--- blocking data
    integer, pointer :: blksz(:)                 !< for explicit data blocking
                                                 !< default blksz(1)=[nx*ny]
!--- ak/bk for pressure level calculations
    real(kind=kind_phys), pointer :: ak(:)       !< from surface (k=1) to TOA (k=levs)
    real(kind=kind_phys), pointer :: bk(:)       !< from surface (k=1) to TOA (k=levs)
!--- grid metrics
    real(kind=kind_phys), pointer :: xlon(:,:)   !< column longitude for MPI rank
    real(kind=kind_phys), pointer :: xlat(:,:)   !< column latitude  for MPI rank
    real(kind=kind_phys), pointer :: area(:,:)   !< column area for length scale calculations

    character(len=32), pointer :: tracer_names(:) !< tracers names to dereference tracer id
                                                  !< based on name location in array
    character(len=65) :: fn_nml                   !< namelist filename
    character(len=256), pointer :: input_nml_file(:) !< character string containing full namelist
                                                     !< for use with internal file reads
  end type GFS_init_type


!----------------------------------------------------------------
! GFS_statein_type
!   prognostic state variables with layer and level specific data
!----------------------------------------------------------------
!! \section arg_table_GFS_statein_type
!! \htmlinclude GFS_statein_type.html
!!
  type GFS_statein_type

!--- level geopotential and pressures 
    real (kind=kind_phys), pointer :: phii  (:,:) => null()   !< interface geopotential height
    real (kind=kind_phys), pointer :: prsi  (:,:) => null()   !< model level pressure in Pa
    real (kind=kind_phys), pointer :: prsik (:,:) => null()   !< Exner function at interface

!--- layer geopotential and pressures
    real (kind=kind_phys), pointer :: phil  (:,:) => null()   !< layer geopotential height
    real (kind=kind_phys), pointer :: prsl  (:,:) => null()   !< model layer mean pressure Pa
    real (kind=kind_phys), pointer :: prslk (:,:) => null()   !< exner function = (p/p0)**rocp

!--- prognostic variables
    real (kind=kind_phys), pointer :: pgr  (:)     => null()  !< surface pressure (Pa) real
    real (kind=kind_phys), pointer :: ugrs (:,:)   => null()  !< u component of layer wind
    real (kind=kind_phys), pointer :: vgrs (:,:)   => null()  !< v component of layer wind
    real (kind=kind_phys), pointer :: vvl  (:,:)   => null()  !< layer mean vertical velocity in pa/sec
    real (kind=kind_phys), pointer :: tgrs (:,:)   => null()  !< model layer mean temperature in k
    real (kind=kind_phys), pointer :: qgrs (:,:,:) => null()  !< layer mean tracer concentration
! dissipation estimate
    real (kind=kind_phys), pointer :: diss_est(:,:)   => null()  !< model layer mean temperature in k
    ! soil state variables - for soil SPPT - sfc-perts, mgehne
    real (kind=kind_phys), pointer :: smc (:,:)   => null()  !< soil moisture content
    real (kind=kind_phys), pointer :: stc (:,:)   => null()  !< soil temperature content
    real (kind=kind_phys), pointer :: slc (:,:)   => null()  !< soil liquid water content
 
    contains
      procedure :: create  => statein_create  !<   allocate array data
  end type GFS_statein_type


!------------------------------------------------------------------
! GFS_stateout_type
!   prognostic state or tendencies after physical parameterizations
!------------------------------------------------------------------
!! \section arg_table_GFS_stateout_type
!! \htmlinclude GFS_stateout_type.html
!!
  type GFS_stateout_type

    !-- Out (physics only)
    real (kind=kind_phys), pointer :: gu0 (:,:)   => null()  !< updated zonal wind
    real (kind=kind_phys), pointer :: gv0 (:,:)   => null()  !< updated meridional wind
    real (kind=kind_phys), pointer :: gt0 (:,:)   => null()  !< updated temperature
    real (kind=kind_phys), pointer :: gq0 (:,:,:) => null()  !< updated tracers

    contains
      procedure :: create  => stateout_create  !<   allocate array data
  end type GFS_stateout_type


!---------------------------------------------------------------------------------------
! GFS_sfcprop_type
!   surface properties that may be read in and/or updated by climatology or observations
!---------------------------------------------------------------------------------------
!! \section arg_table_GFS_sfcprop_type
!! \htmlinclude GFS_sfcprop_type.html
!!
  type GFS_sfcprop_type

!--- In (radiation and physics)
    real (kind=kind_phys), pointer :: slmsk  (:)   => null()  !< sea/land mask array (sea:0,land:1,sea-ice:2)
    real (kind=kind_phys), pointer :: oceanfrac(:) => null()  !< ocean fraction [0:1]
    real (kind=kind_phys), pointer :: landfrac(:)  => null()  !< land  fraction [0:1]
    real (kind=kind_phys), pointer :: lakefrac(:)  => null()  !< lake  fraction [0:1]
    real (kind=kind_phys), pointer :: tsfc   (:)   => null()  !< surface air temperature in K
                                                              !< [tsea in gbphys.f]
    real (kind=kind_phys), pointer :: tsfco  (:)   => null()  !< sst in K
    real (kind=kind_phys), pointer :: tsfcl  (:)   => null()  !< surface land temperature in K
    real (kind=kind_phys), pointer :: tisfc  (:)   => null()  !< surface temperature over ice fraction 
    real (kind=kind_phys), pointer :: snowd  (:)   => null()  !< snow depth water equivalent in mm ; same as snwdph
    real (kind=kind_phys), pointer :: zorl   (:)   => null()  !< composite surface roughness in cm 
    real (kind=kind_phys), pointer :: zorlo  (:)   => null()  !< ocean surface roughness in cm 
    real (kind=kind_phys), pointer :: zorll  (:)   => null()  !< land surface roughness in cm 
    real (kind=kind_phys), pointer :: fice   (:)   => null()  !< ice fraction over open water grid 
!   real (kind=kind_phys), pointer :: hprim  (:)   => null()  !< topographic standard deviation in m
    real (kind=kind_phys), pointer :: hprime (:,:) => null()  !< orographic metrics

!--- In (radiation only)
    real (kind=kind_phys), pointer :: sncovr (:)   => null()  !< snow cover in fraction
    real (kind=kind_phys), pointer :: snoalb (:)   => null()  !< maximum snow albedo in fraction
    real (kind=kind_phys), pointer :: alvsf  (:)   => null()  !< mean vis albedo with strong cosz dependency
    real (kind=kind_phys), pointer :: alnsf  (:)   => null()  !< mean nir albedo with strong cosz dependency
    real (kind=kind_phys), pointer :: alvwf  (:)   => null()  !< mean vis albedo with weak cosz dependency
    real (kind=kind_phys), pointer :: alnwf  (:)   => null()  !< mean nir albedo with weak cosz dependency
    real (kind=kind_phys), pointer :: facsf  (:)   => null()  !< fractional coverage with strong cosz dependency
    real (kind=kind_phys), pointer :: facwf  (:)   => null()  !< fractional coverage with   weak cosz dependency

!--- In (physics only)
    real (kind=kind_phys), pointer :: slope  (:)   => null()  !< sfc slope type for lsm
    real (kind=kind_phys), pointer :: shdmin (:)   => null()  !< min fractional coverage of green veg
    real (kind=kind_phys), pointer :: shdmax (:)   => null()  !< max fractnl cover of green veg (not used)
    real (kind=kind_phys), pointer :: tg3    (:)   => null()  !< deep soil temperature
    real (kind=kind_phys), pointer :: vfrac  (:)   => null()  !< vegetation fraction
    real (kind=kind_phys), pointer :: vtype  (:)   => null()  !< vegetation type
    real (kind=kind_phys), pointer :: stype  (:)   => null()  !< soil type
    real (kind=kind_phys), pointer :: uustar (:)   => null()  !< boundary layer parameter
    real (kind=kind_phys), pointer :: oro    (:)   => null()  !< orography 
    real (kind=kind_phys), pointer :: oro_uf (:)   => null()  !< unfiltered orography

    real (kind=kind_phys), pointer :: qss_ice(:)       => null()  !<
    real (kind=kind_phys), pointer :: qss_land(:)      => null()  !<
    real (kind=kind_phys), pointer :: qss_ocean(:)     => null()  !<
    real (kind=kind_phys), pointer :: snowd_ice(:)     => null()  !<
    real (kind=kind_phys), pointer :: snowd_land(:)    => null()  !<
    real (kind=kind_phys), pointer :: snowd_ocean(:)    => null() !<
    real (kind=kind_phys), pointer :: tsfc_ice(:)      => null()  !<
    real (kind=kind_phys), pointer :: tsfc_land(:)     => null()  !<
    real (kind=kind_phys), pointer :: tsfc_ocean(:)    => null()  !<
    real (kind=kind_phys), pointer :: tsurf_ice(:)     => null()  !<
    real (kind=kind_phys), pointer :: tsurf_land(:)    => null()  !<
    real (kind=kind_phys), pointer :: tsurf_ocean(:)   => null()  !<
    real (kind=kind_phys), pointer :: uustar_ice(:)    => null()  !<
    real (kind=kind_phys), pointer :: uustar_land(:)   => null()  !<
    real (kind=kind_phys), pointer :: uustar_ocean(:)  => null()  !<
    real (kind=kind_phys), pointer :: zorl_ice(:)      => null()  !<
    real (kind=kind_phys), pointer :: zorl_land(:)     => null()  !<
    real (kind=kind_phys), pointer :: zorl_ocean(:)    => null()  !<
    real (kind=kind_phys), pointer :: evap(:)          => null()  !<
    real (kind=kind_phys), pointer :: hflx(:)          => null()  !<

!-- In/Out
#ifdef CCPP
    real (kind=kind_phys), pointer :: conv_act(:)  => null()  !< convective activity counter hli 09/2017
#endif
    real (kind=kind_phys), pointer :: hice   (:)   => null()  !< sea ice thickness
    real (kind=kind_phys), pointer :: weasd  (:)   => null()  !< water equiv of accumulated snow depth (kg/m**2)
                                                              !< over land and sea ice
    real (kind=kind_phys), pointer :: canopy (:)   => null()  !< canopy water
    real (kind=kind_phys), pointer :: ffmm   (:)   => null()  !< fm parameter from PBL scheme
    real (kind=kind_phys), pointer :: ffhh   (:)   => null()  !< fh parameter from PBL scheme
    real (kind=kind_phys), pointer :: f10m   (:)   => null()  !< fm at 10m - Ratio of sigma level 1 wind and 10m wind
    real (kind=kind_phys), pointer :: tprcp  (:)   => null()  !< sfc_fld%tprcp - total precipitation
    real (kind=kind_phys), pointer :: srflag (:)   => null()  !< sfc_fld%srflag - snow/rain flag for precipitation
    real (kind=kind_phys), pointer :: slc    (:,:) => null()  !< liquid soil moisture
    real (kind=kind_phys), pointer :: smc    (:,:) => null()  !< total soil moisture
    real (kind=kind_phys), pointer :: stc    (:,:) => null()  !< soil temperature

!--- Out
    real (kind=kind_phys), pointer :: t2m    (:)   => null()  !< 2 meter temperature
#ifdef CCPP
    real (kind=kind_phys), pointer :: th2m   (:)   => null()  !< 2 meter potential temperature
#endif
    real (kind=kind_phys), pointer :: q2m    (:)   => null()  !< 2 meter humidity

! -- In/Out for Noah MP 
    real (kind=kind_phys), pointer  :: snowxy (:)  => null() !
    real (kind=kind_phys), pointer :: tvxy    (:)  => null()  !< veg temp
    real (kind=kind_phys), pointer :: tgxy    (:)  => null()  !< ground temp
    real (kind=kind_phys), pointer :: canicexy(:)  => null()  !<
    real (kind=kind_phys), pointer :: canliqxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: eahxy   (:)  => null()  !<
    real (kind=kind_phys), pointer :: tahxy   (:)  => null()  !<
    real (kind=kind_phys), pointer :: cmxy    (:)  => null()  !<
    real (kind=kind_phys), pointer :: chxy    (:)  => null()  !<
    real (kind=kind_phys), pointer :: fwetxy  (:)  => null()  !<
    real (kind=kind_phys), pointer :: sneqvoxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: alboldxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: qsnowxy (:)  => null()  !<
    real (kind=kind_phys), pointer :: wslakexy(:)  => null()  !<
    real (kind=kind_phys), pointer :: zwtxy   (:)  => null()  !<
    real (kind=kind_phys), pointer :: waxy    (:)  => null()  !<
    real (kind=kind_phys), pointer :: wtxy    (:)  => null()  !<
    real (kind=kind_phys), pointer :: lfmassxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: rtmassxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: stmassxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: woodxy  (:)  => null()  !<
    real (kind=kind_phys), pointer :: stblcpxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: fastcpxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: xsaixy  (:)  => null()  !<
    real (kind=kind_phys), pointer :: xlaixy  (:)  => null()  !<
    real (kind=kind_phys), pointer :: taussxy (:)  => null()  !<
    real (kind=kind_phys), pointer :: smcwtdxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: deeprechxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: rechxy  (:)  => null()  !<

    real (kind=kind_phys), pointer :: snicexy   (:,:) => null()  !<
    real (kind=kind_phys), pointer :: snliqxy   (:,:) => null()  !<
    real (kind=kind_phys), pointer :: tsnoxy    (:,:) => null()  !<
    real (kind=kind_phys), pointer :: smoiseq   (:,:) => null()  !<
    real (kind=kind_phys), pointer :: zsnsoxy   (:,:) => null()  !<



!--- NSSTM variables  (only allocated when [Model%nstf_name(1) > 0])
    real (kind=kind_phys), pointer :: tref   (:)   => null()  !< nst_fld%Tref - Reference Temperature
    real (kind=kind_phys), pointer :: z_c    (:)   => null()  !< nst_fld%z_c - Sub layer cooling thickness
    real (kind=kind_phys), pointer :: c_0    (:)   => null()  !< nst_fld%c_0 - coefficient1 to calculate d(Tz)/d(Ts)
    real (kind=kind_phys), pointer :: c_d    (:)   => null()  !< nst_fld%c_d - coefficient2 to calculate d(Tz)/d(Ts)
    real (kind=kind_phys), pointer :: w_0    (:)   => null()  !< nst_fld%w_0 - coefficient3 to calculate d(Tz)/d(Ts)
    real (kind=kind_phys), pointer :: w_d    (:)   => null()  !< nst_fld%w_d - coefficient4 to calculate d(Tz)/d(Ts)
    real (kind=kind_phys), pointer :: xt     (:)   => null()  !< nst_fld%xt      heat content in DTL
    real (kind=kind_phys), pointer :: xs     (:)   => null()  !< nst_fld%xs      salinity  content in DTL
    real (kind=kind_phys), pointer :: xu     (:)   => null()  !< nst_fld%xu      u current content in DTL
    real (kind=kind_phys), pointer :: xv     (:)   => null()  !< nst_fld%xv      v current content in DTL
    real (kind=kind_phys), pointer :: xz     (:)   => null()  !< nst_fld%xz      DTL thickness
    real (kind=kind_phys), pointer :: zm     (:)   => null()  !< nst_fld%zm      MXL thickness
    real (kind=kind_phys), pointer :: xtts   (:)   => null()  !< nst_fld%xtts    d(xt)/d(ts)
    real (kind=kind_phys), pointer :: xzts   (:)   => null()  !< nst_fld%xzts    d(xz)/d(ts)
    real (kind=kind_phys), pointer :: d_conv (:)   => null()  !< nst_fld%d_conv  thickness of Free Convection Layer (FCL)
    real (kind=kind_phys), pointer :: ifd    (:)   => null()  !< nst_fld%ifd     index to start DTM run or not
    real (kind=kind_phys), pointer :: dt_cool(:)   => null()  !< nst_fld%dt_cool Sub layer cooling amount
    real (kind=kind_phys), pointer :: qrain  (:)   => null()  !< nst_fld%qrain   sensible heat flux due to rainfall (watts)

#ifdef CCPP
    ! Soil properties for RUC LSM (number of levels different from NOAH 4-layer model)
    real (kind=kind_phys), pointer :: wetness(:)       => null()  !< normalized soil wetness for lsm
    real (kind=kind_phys), pointer :: sh2o(:,:)        => null()  !< volume fraction of unfrozen soil moisture for lsm
    real (kind=kind_phys), pointer :: keepsmfr(:,:)    => null()  !< RUC LSM: frozen moisture in soil
    real (kind=kind_phys), pointer :: smois(:,:)       => null()  !< volumetric fraction of soil moisture for lsm
    real (kind=kind_phys), pointer :: tslb(:,:)        => null()  !< soil temperature for land surface model
    real (kind=kind_phys), pointer :: flag_frsoil(:,:) => null()  !< RUC LSM: flag for frozen soil physics
    !
    real (kind=kind_phys), pointer :: zs(:)            => null()  !< depth of soil levels for land surface model
    real (kind=kind_phys), pointer :: clw_surf(:)      => null()  !< RUC LSM: moist cloud water mixing ratio at surface
    real (kind=kind_phys), pointer :: qwv_surf(:)      => null()  !< RUC LSM: water vapor mixing ratio at surface
    real (kind=kind_phys), pointer :: cndm_surf(:)     => null()  !< RUC LSM: surface condensation mass
    real (kind=kind_phys), pointer :: rhofr(:)         => null()  !< RUC LSM: density of frozen precipitation
    real (kind=kind_phys), pointer :: tsnow(:)         => null()  !< RUC LSM: snow temperature at the bottom of the first soil layer
    real (kind=kind_phys), pointer :: snowfallac(:)    => null()  !< ruc lsm diagnostics
    real (kind=kind_phys), pointer :: acsnow(:)        => null()  !< ruc lsm diagnostics

    !  MYNN surface layer
    real (kind=kind_phys), pointer :: ustm (:)         => null()  !u* including drag
    real (kind=kind_phys), pointer :: zol(:)           => null()  !surface stability parameter
    real (kind=kind_phys), pointer :: mol(:)           => null()  !theta star
    real (kind=kind_phys), pointer :: rmol(:)          => null()  !reciprocal of obukhov length
    real (kind=kind_phys), pointer :: flhc(:)          => null()  !drag coeff for heat
    real (kind=kind_phys), pointer :: flqc(:)          => null()  !drag coeff for moisture
    real (kind=kind_phys), pointer :: chs2(:)          => null()  !exch coeff for heat at 2m
    real (kind=kind_phys), pointer :: cqs2(:)          => null()  !exch coeff for moisture at 2m
    real (kind=kind_phys), pointer :: lh(:)            => null()  !latent heating at the surface
#endif

    !---- precipitation amounts from previous time step for RUC LSM/NoahMP LSM
    real (kind=kind_phys), pointer :: raincprv  (:)    => null()  !< explicit rainfall from previous timestep
    real (kind=kind_phys), pointer :: rainncprv (:)    => null()  !< convective_precipitation_amount from previous timestep
    real (kind=kind_phys), pointer :: iceprv    (:)    => null()  !< ice amount from previous timestep
    real (kind=kind_phys), pointer :: snowprv   (:)    => null()  !< snow amount from previous timestep
    real (kind=kind_phys), pointer :: graupelprv(:)    => null()  !< graupel amount from previous timestep

    !---- precipitation rates from previous time step for NoahMP LSM
    real (kind=kind_phys), pointer :: draincprv  (:)    => null()  !< convective precipitation rate from previous timestep
    real (kind=kind_phys), pointer :: drainncprv (:)    => null()  !< explicit rainfall rate from previous timestep
    real (kind=kind_phys), pointer :: diceprv    (:)    => null()  !< ice precipitation rate from previous timestep
    real (kind=kind_phys), pointer :: dsnowprv   (:)    => null()  !< snow precipitation rate from previous timestep
    real (kind=kind_phys), pointer :: dgraupelprv(:)    => null()  !< graupel precipitation rate from previous timestep

    contains
      procedure :: create  => sfcprop_create  !<   allocate array data
  end type GFS_sfcprop_type


!---------------------------------------------------------------------
! GFS_coupling_type
!   fields to/from other coupled components (e.g. land/ice/ocean/etc.)
!---------------------------------------------------------------------
!! \section arg_table_GFS_coupling_type
!! \htmlinclude GFS_coupling_type.html
!!
  type GFS_coupling_type

!--- Out (radiation only)
    real (kind=kind_phys), pointer :: nirbmdi(:)     => null()   !< sfc nir beam sw downward flux (w/m2) 
    real (kind=kind_phys), pointer :: nirdfdi(:)     => null()   !< sfc nir diff sw downward flux (w/m2)
    real (kind=kind_phys), pointer :: visbmdi(:)     => null()   !< sfc uv+vis beam sw downward flux (w/m2)
    real (kind=kind_phys), pointer :: visdfdi(:)     => null()   !< sfc uv+vis diff sw downward flux (w/m2)
    real (kind=kind_phys), pointer :: nirbmui(:)     => null()   !< sfc nir beam sw upward flux (w/m2)
    real (kind=kind_phys), pointer :: nirdfui(:)     => null()   !< sfc nir diff sw upward flux (w/m2)
    real (kind=kind_phys), pointer :: visbmui(:)     => null()   !< sfc uv+vis beam sw upward flux (w/m2)
    real (kind=kind_phys), pointer :: visdfui(:)     => null()   !< sfc uv+vis diff sw upward flux (w/m2)

    !--- In (physics only)
    real (kind=kind_phys), pointer :: sfcdsw(:)      => null()   !< total sky sfc downward sw flux ( w/m**2 )
                                                                 !< GFS_radtend_type%sfcfsw%dnfxc
    real (kind=kind_phys), pointer :: sfcnsw(:)      => null()   !< total sky sfc netsw flx into ground(w/m**2)
                                                                 !< difference of dnfxc & upfxc from GFS_radtend_type%sfcfsw
    real (kind=kind_phys), pointer :: sfcdlw(:)      => null()   !< total sky sfc downward lw flux ( w/m**2 )
                                                                 !< GFS_radtend_type%sfclsw%dnfxc

!--- incoming quantities
    real (kind=kind_phys), pointer :: dusfcin_cpl(:) => null()   !< aoi_fld%dusfcin(item,lan)
    real (kind=kind_phys), pointer :: dvsfcin_cpl(:) => null()   !< aoi_fld%dvsfcin(item,lan)
    real (kind=kind_phys), pointer :: dtsfcin_cpl(:) => null()   !< aoi_fld%dtsfcin(item,lan)
    real (kind=kind_phys), pointer :: dqsfcin_cpl(:) => null()   !< aoi_fld%dqsfcin(item,lan)
    real (kind=kind_phys), pointer :: ulwsfcin_cpl(:)=> null()   !< aoi_fld%ulwsfcin(item,lan)
    real (kind=kind_phys), pointer :: tseain_cpl(:)  => null()   !< aoi_fld%tseain(item,lan)
    real (kind=kind_phys), pointer :: tisfcin_cpl(:) => null()   !< aoi_fld%tisfcin(item,lan)
    real (kind=kind_phys), pointer :: ficein_cpl(:)  => null()   !< aoi_fld%ficein(item,lan)
    real (kind=kind_phys), pointer :: hicein_cpl(:)  => null()   !< aoi_fld%hicein(item,lan)
    real (kind=kind_phys), pointer :: hsnoin_cpl(:)  => null()   !< aoi_fld%hsnoin(item,lan)
    !--- only variable needed for cplwav=.TRUE.
    !--- also needed for ice/ocn coupling - Xingren
    real (kind=kind_phys), pointer :: slimskin_cpl(:)=> null()   !< aoi_fld%slimskin(item,lan)

!--- outgoing accumulated quantities
    real (kind=kind_phys), pointer :: rain_cpl  (:)  => null()   !< total rain precipitation
    real (kind=kind_phys), pointer :: rainc_cpl (:)  => null()   !< convective rain precipitation
    real (kind=kind_phys), pointer :: snow_cpl  (:)  => null()   !< total snow precipitation  
    real (kind=kind_phys), pointer :: dusfc_cpl (:)  => null()   !< sfc u momentum flux
    real (kind=kind_phys), pointer :: dvsfc_cpl (:)  => null()   !< sfc v momentum flux
    real (kind=kind_phys), pointer :: dtsfc_cpl (:)  => null()   !< sfc sensible heat flux
    real (kind=kind_phys), pointer :: dqsfc_cpl (:)  => null()   !< sfc   latent heat flux
    real (kind=kind_phys), pointer :: dlwsfc_cpl(:)  => null()   !< sfc downward lw flux (w/m**2)
    real (kind=kind_phys), pointer :: dswsfc_cpl(:)  => null()   !< sfc downward sw flux (w/m**2)
    real (kind=kind_phys), pointer :: dnirbm_cpl(:)  => null()   !< sfc nir beam downward sw flux (w/m**2)
    real (kind=kind_phys), pointer :: dnirdf_cpl(:)  => null()   !< sfc nir diff downward sw flux (w/m**2)
    real (kind=kind_phys), pointer :: dvisbm_cpl(:)  => null()   !< sfc uv+vis beam dnwd sw flux (w/m**2)
    real (kind=kind_phys), pointer :: dvisdf_cpl(:)  => null()   !< sfc uv+vis diff dnwd sw flux (w/m**2)
    real (kind=kind_phys), pointer :: nlwsfc_cpl(:)  => null()   !< net downward lw flux (w/m**2)
    real (kind=kind_phys), pointer :: nswsfc_cpl(:)  => null()   !< net downward sw flux (w/m**2)
    real (kind=kind_phys), pointer :: nnirbm_cpl(:)  => null()   !< net nir beam downward sw flux (w/m**2)
    real (kind=kind_phys), pointer :: nnirdf_cpl(:)  => null()   !< net nir diff downward sw flux (w/m**2)
    real (kind=kind_phys), pointer :: nvisbm_cpl(:)  => null()   !< net uv+vis beam downward sw rad flux (w/m**2)
    real (kind=kind_phys), pointer :: nvisdf_cpl(:)  => null()   !< net uv+vis diff downward sw rad flux (w/m**2)

!--- outgoing instantaneous quantities
    real (kind=kind_phys), pointer :: dusfci_cpl (:) => null()   !< instantaneous sfc u momentum flux
    real (kind=kind_phys), pointer :: dvsfci_cpl (:) => null()   !< instantaneous sfc v momentum flux
    real (kind=kind_phys), pointer :: dtsfci_cpl (:) => null()   !< instantaneous sfc sensible heat flux
    real (kind=kind_phys), pointer :: dqsfci_cpl (:) => null()   !< instantaneous sfc   latent heat flux
    real (kind=kind_phys), pointer :: dlwsfci_cpl(:) => null()   !< instantaneous sfc downward lw flux
    real (kind=kind_phys), pointer :: dswsfci_cpl(:) => null()   !< instantaneous sfc downward sw flux
    real (kind=kind_phys), pointer :: dnirbmi_cpl(:) => null()   !< instantaneous sfc nir beam downward sw flux
    real (kind=kind_phys), pointer :: dnirdfi_cpl(:) => null()   !< instantaneous sfc nir diff downward sw flux
    real (kind=kind_phys), pointer :: dvisbmi_cpl(:) => null()   !< instantaneous sfc uv+vis beam downward sw flux
    real (kind=kind_phys), pointer :: dvisdfi_cpl(:) => null()   !< instantaneous sfc uv+vis diff downward sw flux
    real (kind=kind_phys), pointer :: nlwsfci_cpl(:) => null()   !< instantaneous net sfc downward lw flux
    real (kind=kind_phys), pointer :: nswsfci_cpl(:) => null()   !< instantaneous net sfc downward sw flux
    real (kind=kind_phys), pointer :: nnirbmi_cpl(:) => null()   !< instantaneous net nir beam sfc downward sw flux
    real (kind=kind_phys), pointer :: nnirdfi_cpl(:) => null()   !< instantaneous net nir diff sfc downward sw flux
    real (kind=kind_phys), pointer :: nvisbmi_cpl(:) => null()   !< instantaneous net uv+vis beam downward sw flux
    real (kind=kind_phys), pointer :: nvisdfi_cpl(:) => null()   !< instantaneous net uv+vis diff downward sw flux
    real (kind=kind_phys), pointer :: t2mi_cpl   (:) => null()   !< instantaneous T2m
    real (kind=kind_phys), pointer :: q2mi_cpl   (:) => null()   !< instantaneous Q2m
    real (kind=kind_phys), pointer :: u10mi_cpl  (:) => null()   !< instantaneous U10m
    real (kind=kind_phys), pointer :: v10mi_cpl  (:) => null()   !< instantaneous V10m
    real (kind=kind_phys), pointer :: tsfci_cpl  (:) => null()   !< instantaneous sfc temperature
    real (kind=kind_phys), pointer :: psurfi_cpl (:) => null()   !< instantaneous sfc pressure

    !--- topography-based information for the coupling system
    real (kind=kind_phys), pointer :: oro_cpl    (:) => null()   !< orography          (  oro from GFS_sfcprop_type)
    real (kind=kind_phys), pointer :: slmsk_cpl  (:) => null()   !< Land/Sea/Ice mask  (slmsk from GFS_sfcprop_type)

    !--- cellular automata
    real (kind=kind_phys), pointer :: tconvtend(:,:) => null()
    real (kind=kind_phys), pointer :: qconvtend(:,:) => null()
    real (kind=kind_phys), pointer :: uconvtend(:,:) => null()
    real (kind=kind_phys), pointer :: vconvtend(:,:) => null()
    real (kind=kind_phys), pointer :: ca_out   (:)   => null() !
    real (kind=kind_phys), pointer :: ca_deep  (:)   => null() !
    real (kind=kind_phys), pointer :: ca_turb  (:)   => null() !
    real (kind=kind_phys), pointer :: ca_shal  (:)   => null() !
    real (kind=kind_phys), pointer :: ca_rad   (:)   => null() !
    real (kind=kind_phys), pointer :: ca_micro (:)   => null() !
    real (kind=kind_phys), pointer :: cape     (:)   => null() !
    
    !--- stochastic physics
    real (kind=kind_phys), pointer :: shum_wts  (:,:) => null()  !
    real (kind=kind_phys), pointer :: sppt_wts  (:,:) => null()  !
    real (kind=kind_phys), pointer :: skebu_wts (:,:) => null()  !
    real (kind=kind_phys), pointer :: skebv_wts (:,:) => null()  !
    real (kind=kind_phys), pointer :: sfc_wts   (:,:) => null()  ! mg, sfc-perts
    integer              :: nsfcpert=6                             !< number of sfc perturbations

    !--- aerosol surface emissions for Thompson microphysics
    real (kind=kind_phys), pointer :: nwfa2d  (:)     => null()  !< instantaneous water-friendly sfc aerosol source
    real (kind=kind_phys), pointer :: nifa2d  (:)     => null()  !< instantaneous ice-friendly sfc aerosol source

    !--- instantaneous quantities for GSDCHEM coupling
    real (kind=kind_phys), pointer :: dqdti   (:,:)   => null()  !< instantaneous total moisture tendency (kg/kg/s)
    real (kind=kind_phys), pointer :: ushfsfci(:)     => null()  !< instantaneous upward sensible heat flux (w/m**2)
    real (kind=kind_phys), pointer :: dkt     (:,:)   => null()  !< instantaneous dkt diffusion coefficient for temperature (m**2/s)
    real (kind=kind_phys), pointer :: qci_conv(:,:)   => null()  !< convective cloud condesate after rainout


    contains
      procedure :: create  => coupling_create  !<   allocate array data
  end type GFS_coupling_type


!----------------------------------------------------------------------------------
! GFS_control_type
!   model control parameters input from a namelist and/or derived from others
!   list of those that can be modified during the run are at the bottom of the list 
!----------------------------------------------------------------------------------
!! \section arg_table_GFS_control_type
!! \htmlinclude GFS_control_type.html
!!
  type GFS_control_type

    integer              :: me              !< MPI rank designator
    integer              :: master          !< MPI rank of master atmosphere processor
#ifdef CCPP
    integer              :: communicator    !< MPI communicator
    integer              :: ntasks          !< MPI size in communicator
    integer              :: nthreads        !< OpenMP threads available for physics
#endif
    integer              :: nlunit          !< unit for namelist
    character(len=64)    :: fn_nml          !< namelist filename for surface data cycling
    character(len=256), pointer :: input_nml_file(:) !< character string containing full namelist
                                                   !< for use with internal file reads
#ifdef CCPP
    integer              :: input_nml_file_length
    integer              :: logunit
#endif
    real(kind=kind_phys) :: fhzero          !< hours between clearing of diagnostic buckets
    logical              :: ldiag3d         !< flag for 3d diagnostic fields
    logical              :: lssav           !< logical flag for storing diagnostics
    real(kind=kind_phys) :: fhcyc           !< frequency for surface data cycling (hours)
    integer              :: thermodyn_id    !< valid for GFS only for get_prs/phi
    integer              :: sfcpress_id     !< valid for GFS only for get_prs/phi
    logical              :: gen_coord_hybrid!< for Henry's gen coord

!--- set some grid extent parameters
    integer              :: isc             !< starting i-index for this MPI-domain
    integer              :: jsc             !< starting j-index for this MPI-domain
    integer              :: nx              !< number of points in the i-dir for this MPI-domain
    integer              :: ny              !< number of points in the j-dir for this MPI-domain
    integer              :: levs            !< number of vertical levels
#ifdef CCPP
    integer              :: levsp1          !< number of vertical levels plus one
    integer              :: levsm1          !< number of vertical levels minus one
    !--- ak/bk for pressure level calculations
    real(kind=kind_phys), pointer :: ak(:)  !< from surface (k=1) to TOA (k=levs)
    real(kind=kind_phys), pointer :: bk(:)  !< from surface (k=1) to TOA (k=levs)
#endif
    integer              :: cnx             !< number of points in the i-dir for this cubed-sphere face
    integer              :: cny             !< number of points in the j-dir for this cubed-sphere face
    integer              :: lonr            !< number of global points in x-dir (i) along the equator
    integer              :: latr            !< number of global points in y-dir (j) along any meridian
    integer              :: tile_num
    integer              :: nblks           !< for explicit data blocking: number of blocks
    integer,     pointer :: blksz(:)        !< for explicit data blocking: block sizes of all blocks
#ifdef CCPP
    integer,     pointer :: blksz2(:)       !< for explicit data blocking: block sizes of all blocks (duplicate)
#endif

!--- coupling parameters
    logical              :: cplflx          !< default no cplflx collection
    logical              :: cplwav          !< default no cplwav collection
    logical              :: cplchm          !< default no cplchm collection

!--- integrated dynamics through earth's atmosphere
    logical              :: lsidea         

!vay 2018  GW physics switches

    logical              :: ldiag_ugwp
    logical              :: do_ugwp         ! do mesoscale UGWP + TOFD + RF
    logical              :: do_tofd         ! tofd flag in gwdps.f
    logical              :: do_gwd          ! logical for gravity wave drag (gwd)
    logical              :: do_cnvgwd       ! logical for convective gwd

!--- calendars and time parameters and activation triggers
    real(kind=kind_phys) :: dtp             !< physics timestep in seconds
    real(kind=kind_phys) :: dtf             !< dynamics timestep in seconds
    integer              :: nscyc           !< trigger for surface data cycling
    integer              :: nszero          !< trigger for zeroing diagnostic buckets
    integer              :: idat(1:8)       !< initialization date and time
                                            !< (yr, mon, day, t-zone, hr, min, sec, mil-sec)
    integer              :: idate(4)        !< initial date with different size and ordering
                                            !< (hr, mon, day, yr)
!--- radiation control parameters
    real(kind=kind_phys) :: fhswr           !< frequency for shortwave radiation (secs)
    real(kind=kind_phys) :: fhlwr           !< frequency for longwave radiation (secs)
    integer              :: nsswr           !< integer trigger for shortwave radiation
    integer              :: nslwr           !< integer trigger for longwave  radiation
#ifdef CCPP
    integer              :: nhfrad          !< number of timesteps for which to call radiation on physics timestep (coldstarts)
#endif
    integer              :: levr            !< number of vertical levels for radiation calculations
#ifdef CCPP
    integer              :: levrp1          !< number of vertical levels for radiation calculations plus one
#endif
    integer              :: nfxr            !< second dimension for fluxr diagnostic variable (radiation)
    logical              :: aero_in         !< flag for initializing aerosol data
#ifdef CCPP
    integer              :: ntrcaer         !< number of aerosol tracers for Morrison-Gettelman microphysics
#endif
    logical              :: lmfshal         !< parameter for radiation
    logical              :: lmfdeep2        !< parameter for radiation
    integer              :: nrcm            !< second dimension of random number stream for RAS
    integer              :: iflip           !< iflip - is not the same as flipv
    integer              :: isol            !< use prescribed solar constant
    integer              :: ico2            !< prescribed global mean value (old opernl)
    integer              :: ialb            !< use climatology alb, based on sfc type
                                            !< 1 => use modis based alb
    integer              :: iems            !< use fixed value of 1.0
    integer              :: iaer            !< default aerosol effect in sw only
    integer              :: icliq_sw        !< sw optical property for liquid clouds
    integer              :: iovr_sw         !< sw: max-random overlap clouds
    integer              :: iovr_lw         !< lw: max-random overlap clouds
    integer              :: ictm            !< ictm=0 => use data at initial cond time, if not
                                            !<           available; use latest; no extrapolation.
                                            !< ictm=1 => use data at the forecast time, if not
                                            !<           available; use latest; do extrapolation.
                                            !< ictm=yyyy0 => use yyyy data for the forecast time;
                                            !<           no extrapolation.
                                            !< ictm=yyyy1 = > use yyyy data for the fcst. If needed, 
                                            !<           do extrapolation to match the fcst time.
                                            !< ictm=-1 => use user provided external data for
                                            !<           the fcst time; no extrapolation.
                                            !< ictm=-2 => same as ictm=0, but add seasonal cycle
                                            !<           from climatology; no extrapolation.
    integer              :: isubc_sw        !< sw clouds without sub-grid approximation
    integer              :: isubc_lw        !< lw clouds without sub-grid approximation
                                            !< =1 => sub-grid cloud with prescribed seeds
                                            !< =2 => sub-grid cloud with randomly generated
                                            !< seeds
    logical              :: crick_proof     !< CRICK-Proof cloud water
    logical              :: ccnorm          !< Cloud condensate normalized by cloud cover 
    logical              :: norad_precip    !< radiation precip flag for Ferrier/Moorthi
    logical              :: lwhtr           !< flag to output lw heating rate (Radtend%lwhc)
    logical              :: swhtr           !< flag to output sw heating rate (Radtend%swhc)

!--- microphysical switch
    integer              :: ncld            !< choice of cloud scheme
    !--- new microphysical switch
    integer              :: imp_physics                    !< choice of microphysics scheme
    integer              :: imp_physics_gfdl = 11          !< choice of GFDL     microphysics scheme
    integer              :: imp_physics_thompson = 8       !< choice of Thompson microphysics scheme
    integer              :: imp_physics_wsm6 = 6           !< choice of WSMG     microphysics scheme
    integer              :: imp_physics_zhao_carr = 99     !< choice of Zhao-Carr microphysics scheme
    integer              :: imp_physics_zhao_carr_pdf = 98 !< choice of Zhao-Carr microphysics scheme with PDF clouds
    integer              :: imp_physics_mg = 10            !< choice of Morrison-Gettelman microphysics scheme
    integer              :: imp_physics_fer_hires = 15     !< choice of Ferrier-Aligo microphysics scheme
    !--- Z-C microphysical parameters
    real(kind=kind_phys) :: psautco(2)         !< [in] auto conversion coeff from ice to snow
    real(kind=kind_phys) :: prautco(2)         !< [in] auto conversion coeff from cloud to rain
    real(kind=kind_phys) :: evpco              !< [in] coeff for evaporation of largescale rain
    real(kind=kind_phys) :: wminco(2)          !< [in] water and ice minimum threshold for Zhao
    real(kind=kind_phys) :: avg_max_length     !< reset time in seconds for max hourly fields
    !--- M-G microphysical parameters
    integer              :: fprcp              !< no prognostic rain and snow (MG)
    integer              :: pdfflag            !< pdf flag for MG macrophysics
    real(kind=kind_phys) :: mg_dcs             !< Morrison-Gettelman microphysics parameters
    real(kind=kind_phys) :: mg_qcvar      
    real(kind=kind_phys) :: mg_ts_auto_ice(2)  !< ice auto conversion time scale
#ifdef CCPP
    real(kind=kind_phys) :: mg_rhmini          !< relative humidity threshold parameter for nucleating ice
#endif

    real(kind=kind_phys) :: mg_ncnst           !< constant droplet num concentration (m-3)
    real(kind=kind_phys) :: mg_ninst           !< constant ice num concentration (m-3)
    real(kind=kind_phys) :: mg_ngnst           !< constant graupel/hail num concentration (m-3)
    real(kind=kind_phys) :: mg_berg_eff_factor !< berg efficiency factor
    real(kind=kind_phys) :: mg_alf             !< tuning factor for alphs in MG macrophysics
    real(kind=kind_phys) :: mg_qcmin(2)        !< min liquid and ice mixing ratio in Mg macro clouds
    character(len=16)    :: mg_precip_frac_method ! type of precipitation fraction method
#ifdef CCPP
    real(kind=kind_phys) :: tf
    real(kind=kind_phys) :: tcr
    real(kind=kind_phys) :: tcrf
#endif
!
    logical              :: effr_in            !< eg to turn on ffective radii for MG
    logical              :: microp_uniform
    logical              :: do_cldliq
    logical              :: do_cldice
    logical              :: hetfrz_classnuc

    logical              :: mg_nccons
    logical              :: mg_nicons
    logical              :: mg_ngcons
    logical              :: sed_supersat
    logical              :: do_sb_physics
    logical              :: mg_do_graupel
    logical              :: mg_do_hail
    logical              :: mg_do_ice_gmao
    logical              :: mg_do_liq_liu

    real(kind=kind_phys) :: shoc_parm(5)    !< critical pressure in Pa for tke dissipation in shoc
    integer              :: ncnd            !< number of cloud condensate types

    !--- Thompson's microphysical parameters
    logical              :: ltaerosol       !< flag for aerosol version
    logical              :: lradar          !< flag for radar reflectivity 
    real(kind=kind_phys) :: ttendlim        !< temperature tendency limiter per time step in K/s

    !--- GFDL microphysical paramters
    logical              :: lgfdlmprad      !< flag for GFDL mp scheme and radiation consistency 

    !--- Thompson,GFDL mp parameter
    logical              :: lrefres          !< flag for radar reflectivity in restart file

    !--- land/surface model parameters
    integer              :: lsm             !< flag for land surface model lsm=1 for noah lsm
    integer              :: lsm_noah=1      !< flag for NOAH land surface model
    integer              :: lsm_noahmp=2    !< flag for NOAH land surface model
    integer              :: lsm_ruc=3       !< flag for RUC land surface model
    integer              :: lsoil           !< number of soil layers
#ifdef CCPP
    integer              :: lsoil_lsm       !< number of soil layers internal to land surface model
    integer              :: lsnow_lsm       !< maximum number of snow layers internal to land surface model
    integer              :: lsnow_lsm_lbound!< lower bound for snow arrays, depending on lsnow_lsm
    logical              :: rdlai
#endif
    integer              :: ivegsrc         !< ivegsrc = 0   => USGS, 
                                            !< ivegsrc = 1   => IGBP (20 category)
                                            !< ivegsrc = 2   => UMD  (13 category)
    integer              :: isot            !< isot = 0   => Zobler soil type  ( 9 category)
                                            !< isot = 1   => STATSGO soil type (19 category)
    ! -- the Noah MP options

    integer              :: iopt_dveg ! 1-> off table lai 2-> on 3-> off;4->off;5 -> on
    integer              :: iopt_crs  !canopy stomatal resistance (1-> ball-berry; 2->jarvis)
    integer              :: iopt_btr  !soil moisture factor for stomatal resistance (1-> noah; 2-> clm; 3-> ssib)
    integer              :: iopt_run  !runoff and groundwater (1->simgm; 2->simtop; 3->schaake96; 4->bats)
    integer              :: iopt_sfc  !surface layer drag coeff (ch & cm) (1->m-o; 2->chen97)
    integer              :: iopt_frz  !supercooled liquid water (1-> ny06; 2->koren99)
    integer              :: iopt_inf  !frozen soil permeability (1-> ny06; 2->koren99)
    integer              :: iopt_rad  !radiation transfer (1->gap=f(3d,cosz); 2->gap=0; 3->gap=1-fveg)
    integer              :: iopt_alb  !snow surface albedo (1->bats; 2->class)
    integer              :: iopt_snf  !rainfall & snowfall (1-jordan91; 2->bats; 3->noah)
    integer              :: iopt_tbot !lower boundary of soil temperature (1->zero-flux; 2->noah)
    integer              :: iopt_stc  !snow/soil temperature time scheme (only layer 1)

    logical              :: use_ufo         !< flag for gcycle surface option

!--- tuning parameters for physical parameterizations
    logical              :: ras             !< flag for ras convection scheme
    logical              :: flipv           !< flag for vertical direction flip (ras)
                                            !< .true. implies surface at k=1
    logical              :: trans_trac      !< flag for convective transport of tracers (RAS, CS, or SAMF)
    logical              :: old_monin       !< flag for diff monin schemes
    logical              :: cnvgwd          !< flag for conv gravity wave drag
#ifdef CCPP
    integer              :: gwd_opt         !< gwd_opt = 1  => original GFS gwd (gwdps.f)
                                            !< gwd_opt = 2  => unified GWD (placeholder)
                                            !< gwd_opt = 3  => GSD drag suite
                                            !< gwd_opt = 33 => GSD drag suite with extra output
#endif
    logical              :: mstrat          !< flag for moorthi approach for stratus
    logical              :: moist_adj       !< flag for moist convective adjustment
    logical              :: cscnv           !< flag for Chikira-Sugiyama convection
    logical              :: cal_pre         !< flag controls precip type algorithm
#ifdef CCPP
    real(kind=kind_phys) :: rhgrd           !< fer_hires microphysics only
    logical              :: spec_adv        !< flag for individual cloud species advected
#endif
    logical              :: do_aw           !< AW scale-aware option in cs convection
    logical              :: do_awdd         !< AW scale-aware option in cs convection
    logical              :: flx_form        !< AW scale-aware option in cs convection
    logical              :: do_shoc         !< flag for SHOC
    logical              :: shocaftcnv      !< flag for SHOC
    logical              :: shoc_cld        !< flag for clouds
    logical              :: uni_cld         !< flag for clouds in grrad
#ifdef CCPP
    logical              :: oz_phys         !< flag for old (2006) ozone physics
    logical              :: oz_phys_2015    !< flag for new (2015) ozone physics
#endif
    logical              :: h2o_phys        !< flag for stratosphere h2o
    logical              :: pdfcld          !< flag for pdfcld
    logical              :: shcnvcw         !< flag for shallow convective cloud
    logical              :: redrag          !< flag for reduced drag coeff. over sea
    logical              :: hybedmf         !< flag for hybrid edmf pbl scheme
    logical              :: satmedmf        !< flag for scale-aware TKE-based moist edmf
                                            !< vertical turbulent mixing scheme
    logical              :: shinhong        !< flag for scale-aware Shinhong vertical turbulent mixing scheme
    logical              :: do_ysu          !< flag for YSU turbulent mixing scheme
    logical              :: dspheat         !< flag for tke dissipative heating
    logical              :: lheatstrg       !< flag for canopy heat storage parameterization
    logical              :: cnvcld        
    logical              :: random_clds     !< flag controls whether clouds are random
    logical              :: shal_cnv        !< flag for calling shallow convection
    logical              :: do_deep         !< whether to do deep convection
    integer              :: imfshalcnv      !< flag for mass-flux shallow convection scheme
                                            !<     1: July 2010 version of mass-flux shallow conv scheme
                                            !<         current operational version as of 2016
                                            !<     2: scale- & aerosol-aware mass-flux shallow conv scheme (2017)
                                            !<     3: scale- & aerosol-aware Grell-Freitas scheme (GSD)
                                            !<     4: New Tiedtke scheme (CAPS)
                                            !<     0: modified Tiedtke's eddy-diffusion shallow conv scheme
                                            !<    -1: no shallow convection used
#ifdef CCPP
    integer              :: imfshalcnv_sas      = 1 !< flag for SAS mass-flux shallow convection scheme
    integer              :: imfshalcnv_samf     = 2 !< flag for SAMF scale- & aerosol-aware mass-flux shallow convection scheme
    integer              :: imfshalcnv_gf       = 3 !< flag for scale- & aerosol-aware Grell-Freitas scheme (GSD)
    integer              :: imfshalcnv_ntiedtke = 4 !< flag for new Tiedtke scheme (CAPS)
#endif
    integer              :: imfdeepcnv      !< flag for mass-flux deep convection scheme
                                            !<     1: July 2010 version of SAS conv scheme
                                            !<           current operational version as of 2016
                                            !<     2: scale- & aerosol-aware mass-flux deep conv scheme (2017)
                                            !<     3: scale- & aerosol-aware Grell-Freitas scheme (GSD)
                                            !<     4: New Tiedtke scheme (CAPS)
                                            !<     0: old SAS Convection scheme before July 2010
#ifdef CCPP
    integer              :: imfdeepcnv_sas      = 1 !< flag for SAS mass-flux deep convection scheme
    integer              :: imfdeepcnv_samf     = 2 !< flag for SAMF scale- & aerosol-aware mass-flux deep convection scheme
    integer              :: imfdeepcnv_gf       = 3 !< flag for scale- & aerosol-aware Grell-Freitas scheme (GSD)
    integer              :: imfdeepcnv_ntiedtke = 4 !< flag for new Tiedtke scheme (CAPS)
#endif
    integer              :: isatmedmf       !< flag for scale-aware TKE-based moist edmf scheme
                                            !<     0: initial version of satmedmf (Nov. 2018)
                                            !<     1: updated version of satmedmf (as of May 2019)
#ifdef CCPP
    integer              :: isatmedmf_vdif  = 0 !< flag for initial version of satmedmf (Nov. 2018)
    integer              :: isatmedmf_vdifq = 1 !< flag for updated version of satmedmf (as of May 2019)
#endif
    integer              :: nmtvr           !< number of topographic variables such as variance etc
                                            !< used in the GWD parameterization
    integer              :: jcap            !< number of spectral wave trancation used only by sascnv shalcnv
    real(kind=kind_phys) :: cs_parm(10)     !< tunable parameters for Chikira-Sugiyama convection
    real(kind=kind_phys) :: flgmin(2)       !< [in] ice fraction bounds
    real(kind=kind_phys) :: cgwf(2)         !< multiplication factor for convective GWD
    real(kind=kind_phys) :: ccwf(2)         !< multiplication factor for critical cloud
                                            !< workfunction for RAS
    real(kind=kind_phys) :: cdmbgwd(4)      !< multiplication factors for cdmb, gwd and NS gwd, tke based enhancement
    real(kind=kind_phys) :: sup             !< supersaturation in pdf cloud when t is very low
    real(kind=kind_phys) :: ctei_rm(2)      !< critical cloud top entrainment instability criteria 
                                            !< (used if mstrat=.true.)
    real(kind=kind_phys) :: crtrh(3)        !< critical relative humidity at the surface
                                            !< PBL top and at the top of the atmosphere
    real(kind=kind_phys) :: dlqf(2)         !< factor for cloud condensate detrainment 
                                            !< from cloud edges for RAS
    real(kind=kind_phys) :: psauras(2)      !< [in] auto conversion coeff from ice to snow in ras
    real(kind=kind_phys) :: prauras(2)      !< [in] auto conversion coeff from cloud to rain in ras
    real(kind=kind_phys) :: wminras(2)      !< [in] water and ice minimum threshold for ras

    integer              :: seed0           !< random seed for radiation

    real(kind=kind_phys) :: rbcr            !< Critical Richardson Number in the PBL scheme
#ifdef CCPP
    !--- MYNN parameters/switches
    logical              :: do_mynnedmf
    logical              :: do_mynnsfclay
    ! DH* TODO - move this to MYNN namelist section
    integer              :: grav_settling      !< flag for initalizing fist time step
    integer              :: bl_mynn_tkebudget  !< flag for activating TKE budget
    logical              :: bl_mynn_tkeadvect  !< activate computation of TKE advection (not yet in use for FV3)
    integer              :: bl_mynn_cloudpdf   !< flag to determine which cloud PDF to use
    integer              :: bl_mynn_mixlength  !< flag for different version of mixing length formulation
    integer              :: bl_mynn_edmf       !< flag to activate the mass-flux scheme
    integer              :: bl_mynn_edmf_mom   !< flag to activate the transport of momentum
    integer              :: bl_mynn_edmf_tke   !< flag to activate the transport of TKE
    integer              :: bl_mynn_edmf_part  !< flag to partitioning og the MF and ED areas
    integer              :: bl_mynn_cloudmix   !< flag to activate mixing of cloud species
    integer              :: bl_mynn_mixqt      !< flag to mix total water or individual species
    integer              :: icloud_bl          !< flag for coupling sgs clouds to radiation
    ! *DH
    ! MYJ switches
    logical              :: do_myjsfc          !< flag for MYJ surface layer scheme
    logical              :: do_myjpbl          !< flag for MYJ PBL scheme
#endif

!--- Rayleigh friction
    real(kind=kind_phys) :: prslrd0         !< pressure level from which Rayleigh Damping is applied
    real(kind=kind_phys) :: ral_ts          !< time scale for Rayleigh damping in days

!--- mass flux deep convection
    real(kind=kind_phys) :: clam_deep       !< c_e for deep convection (Han and Pan, 2011, eq(6))
    real(kind=kind_phys) :: c0s_deep        !< convective rain conversion parameter
    real(kind=kind_phys) :: c1_deep         !< conversion parameter of detrainment from liquid water into grid-scale cloud water
    real(kind=kind_phys) :: betal_deep      !< fraction factor of downdraft air mass reaching ground surface over land
    real(kind=kind_phys) :: betas_deep      !< fraction factor of downdraft air mass reaching ground surface over sea
    real(kind=kind_phys) :: evfact_deep     !< evaporation factor from convective rain
    real(kind=kind_phys) :: evfactl_deep    !< evaporation factor from convective rain over land
    real(kind=kind_phys) :: pgcon_deep      !< reduction factor in momentum transport due to convection induced pressure gradient force
                                            !< 0.7 : Gregory et al. (1997, QJRMS)
                                            !< 0.55: Zhang & Wu (2003, JAS)
    real(kind=kind_phys) :: asolfac_deep    !< aerosol-aware parameter based on Lim (2011)
                                            !< asolfac= cx / c0s(=.002)
                                            !< cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s)
                                            !< Nccn: CCN number concentration in cm^(-3)
                                            !< Until a realistic Nccn is provided, Nccns are assumed
                                            !< as Nccn=100 for sea and Nccn=1000 for land

!--- mass flux shallow convection
    real(kind=kind_phys) :: clam_shal       !< c_e for shallow convection (Han and Pan, 2011, eq(6))
    real(kind=kind_phys) :: c0s_shal        !< convective rain conversion parameter
    real(kind=kind_phys) :: c1_shal         !< conversion parameter of detrainment from liquid water into grid-scale cloud water
    real(kind=kind_phys) :: pgcon_shal      !< reduction factor in momentum transport due to convection induced pressure gradient force
                                            !< 0.7 : Gregory et al. (1997, QJRMS)
                                            !< 0.55: Zhang & Wu (2003, JAS)
    real(kind=kind_phys) :: asolfac_shal    !< aerosol-aware parameter based on Lim (2011)
                                            !< asolfac= cx / c0s(=.002)
                                            !< cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s)
                                            !< Nccn: CCN number concentration in cm^(-3)
                                            !< Until a realistic Nccn is provided, Nccns are assumed
                                            !< as Nccn=100 for sea and Nccn=1000 for land 

!--- near surface temperature model
    logical              :: nst_anl         !< flag for NSSTM analysis in gcycle/sfcsub
    integer              :: lsea
    integer              :: nstf_name(5)    !< flag 0 for no nst  1 for uncoupled nst  and 2 for coupled NST
                                            !< nstf_name contains the NSST related parameters
                                            !< nstf_name(1) : 0 = NSSTM off, 1 = NSSTM on but uncoupled
                                            !<                2 = NSSTM on and coupled
                                            !< nstf_name(2) : 1 = NSSTM spin up on, 0 = NSSTM spin up off
                                            !< nstf_name(3) : 1 = NSST analysis on, 0 = NSSTM analysis off
                                            !< nstf_name(4) : zsea1 in mm
                                            !< nstf_name(5) : zsea2 in mm
!--- fractional grid
    logical              :: frac_grid       !< flag for fractional grid
    real(kind=kind_phys) :: min_lakeice     !< minimum lake ice value
    real(kind=kind_phys) :: min_seaice      !< minimum sea  ice value
    real(kind=kind_phys) :: rho_h2o         !< density of fresh water

!--- surface layer z0 scheme
    integer              :: sfc_z0_type     !< surface roughness options over ocean: 
                                            !< 0=no change
                                            !< 6=areodynamical roughness over water with input 10-m wind
                                            !< 7=slightly decrease Cd for higher wind speed compare to 6

!--- background vertical diffusion
    real(kind=kind_phys) :: xkzm_m          !< [in] bkgd_vdif_m  background vertical diffusion for momentum
    real(kind=kind_phys) :: xkzm_h          !< [in] bkgd_vdif_h  background vertical diffusion for heat q
    real(kind=kind_phys) :: xkzm_s          !< [in] bkgd_vdif_s  sigma threshold for background mom. diffusion
    real(kind=kind_phys) :: xkzminv         !< diffusivity in inversion layers
    real(kind=kind_phys) :: moninq_fac      !< turbulence diffusion coefficient factor
    real(kind=kind_phys) :: dspfac          !< tke dissipative heating factor
    real(kind=kind_phys) :: bl_upfr         !< updraft fraction in boundary layer mass flux scheme
    real(kind=kind_phys) :: bl_dnfr         !< downdraft fraction in boundary layer mass flux scheme

 !---cellular automata control parameters
    integer              :: nca             !< number of independent cellular automata
    integer              :: nlives          !< cellular automata lifetime
    integer              :: ncells          !< cellular automata finer grid
    real(kind=kind_phys) :: nfracseed       !< cellular automata seed probability
    integer              :: nseed           !< cellular automata seed frequency
    logical              :: do_ca           !< cellular automata main switch
    logical              :: ca_sgs          !< switch for sgs ca
    logical              :: ca_global       !< switch for global ca
    logical              :: ca_smooth       !< switch for gaussian spatial filter
    logical              :: isppt_deep      !< switch for combination with isppt_deep. OBS! Switches off SPPT on other tendencies!
    integer              :: iseed_ca        !< seed for random number generation in ca scheme
    integer              :: nspinup         !< number of iterations to spin up the ca
    real(kind=kind_phys) :: nthresh         !< threshold used for perturbed vertical velocity

!--- stochastic physics control parameters
    logical              :: do_sppt
    logical              :: use_zmtnblck
    logical              :: do_shum
    logical              :: do_skeb
    integer              :: skeb_npass
    logical              :: do_sfcperts
    integer              :: nsfcpert=6
    real(kind=kind_phys) :: pertz0(5)          ! mg, sfc-perts
    real(kind=kind_phys) :: pertzt(5)          ! mg, sfc-perts
    real(kind=kind_phys) :: pertshc(5)         ! mg, sfc-perts
    real(kind=kind_phys) :: pertlai(5)         ! mg, sfc-perts
    real(kind=kind_phys) :: pertalb(5)         ! mg, sfc-perts
    real(kind=kind_phys) :: pertvegf(5)        ! mg, sfc-perts
!--- tracer handling
    character(len=32), pointer :: tracer_names(:) !< array of initialized tracers from dynamic core
    integer              :: ntrac           !< number of tracers
#ifdef CCPP
    integer              :: ntracp1         !< number of tracers plus one
    integer              :: ntqv            !< tracer index for water vapor (specific humidity)
    integer              :: nqrimef         !< tracer index for mass weighted rime factor
#endif
    integer              :: ntoz            !< tracer index for ozone mixing ratio
    integer              :: ntcw            !< tracer index for cloud condensate (or liquid water)
    integer              :: ntiw            !< tracer index for ice water
    integer              :: ntrw            !< tracer index for rain water
    integer              :: ntsw            !< tracer index for snow water
    integer              :: ntgl            !< tracer index for graupel
    integer              :: ntclamt         !< tracer index for cloud amount
    integer              :: ntlnc           !< tracer index for liquid number concentration
    integer              :: ntinc           !< tracer index for ice    number concentration
    integer              :: ntrnc           !< tracer index for rain   number concentration
    integer              :: ntsnc           !< tracer index for snow   number concentration
    integer              :: ntgnc           !< tracer index for graupel number concentration
    integer              :: ntke            !< tracer index for kinetic energy
    integer              :: nto             !< tracer index for oxygen ion
    integer              :: nto2            !< tracer index for oxygen
    integer              :: ntwa            !< tracer index for water friendly aerosol 
    integer              :: ntia            !< tracer index for ice friendly aerosol
    integer              :: ntchm           !< number of chemical tracers
    integer              :: ntchs           !< tracer index for first chemical tracer
    logical, pointer     :: ntdiag(:) => null() !< array to control diagnostics for chemical tracers
    real(kind=kind_phys), pointer :: fscav(:)  => null() !< array of aerosol scavenging coefficients

    !--- derived totals for phy_f*d
    integer              :: ntot2d          !< total number of variables for phyf2d
    integer              :: ntot3d          !< total number of variables for phyf3d
    integer              :: indcld          !< location of cloud fraction in phyf3d (used only for SHOC or MG)
    integer              :: num_p2d         !< number of 2D arrays needed for microphysics
    integer              :: num_p3d         !< number of 3D arrays needed for microphysics
    integer              :: nshoc_2d        !< number of 2d fields for SHOC
    integer              :: nshoc_3d        !< number of 3d fields for SHOC
    integer              :: ncnvcld3d       !< number of convective 3d clouds fields
    integer              :: npdf3d          !< number of 3d arrays associated with pdf based clouds/microphysics
    integer              :: nctp            !< number of cloud types in Chikira-Sugiyama scheme
    integer              :: ncnvw           !< the index of cnvw in phy_f3d
    integer              :: ncnvc           !< the index of cnvc in phy_f3d
    integer              :: nleffr          !< the index of cloud liquid water effective radius in phy_f3d
    integer              :: nieffr          !< the index of ice effective radius in phy_f3d
    integer              :: nreffr          !< the index of rain effective radius in phy_f3d
    integer              :: nseffr          !< the index of snow effective radius in phy_f3d
    integer              :: ngeffr          !< the index of graupel effective radius in phy_f3d
#ifdef CCPP
    integer              :: nkbfshoc        !< the index of upward kinematic buoyancy flux from SHOC in phy_f3d
    integer              :: nahdshoc        !< the index of diffusivity for heat from from SHOC in phy_f3d
    integer              :: nscfshoc        !< the index of subgrid-scale cloud fraction from from SHOC in phy_f3d
#endif

!--- debug flag
    logical              :: debug         
    logical              :: pre_rad         !< flag for testing purpose

!--- variables modified at each time step
    integer              :: ipt             !< index for diagnostic printout point
    logical              :: lprnt           !< control flag for diagnostic print out
    logical              :: lsswr           !< logical flags for sw radiation calls
    logical              :: lslwr           !< logical flags for lw radiation calls
    real(kind=kind_phys) :: solhr           !< hour time after 00z at the t-step
    real(kind=kind_phys) :: solcon          !< solar constant (sun-earth distant adjusted)  [set via radupdate]
    real(kind=kind_phys) :: slag            !< equation of time ( radian )                  [set via radupdate]
    real(kind=kind_phys) :: sdec            !< sin of the solar declination angle           [set via radupdate]
    real(kind=kind_phys) :: cdec            !< cos of the solar declination angle           [set via radupdate]
    real(kind=kind_phys) :: clstp           !< index used by cnvc90 (for convective clouds) 
                                            !< legacy stuff - does not affect forecast
    real(kind=kind_phys) :: phour           !< previous forecast hour
    real(kind=kind_phys) :: fhour           !< current forecast hour
    real(kind=kind_phys) :: zhour           !< previous hour diagnostic buckets emptied
    integer              :: kdt             !< current forecast iteration
#ifdef CCPP
    logical              :: first_time_step !< flag signaling first time step for time integration routine
    logical              :: restart         !< flag whether this is a coldstart (.false.) or a warmstart/restart (.true.)
    logical              :: hydrostatic     !< flag whether this is a hydrostatic or non-hydrostatic run
#endif
    integer              :: jdat(1:8)       !< current forecast date and time
                                            !< (yr, mon, day, t-zone, hr, min, sec, mil-sec)
    integer              :: imn             !< current forecast month
    real(kind=kind_phys) :: julian          !< current forecast julian date
    integer              :: yearlen         !< current length of the year
!
    logical              :: iccn            !< using IN CCN forcing for MG2/3
#ifdef CCPP
    real(kind=kind_phys)          :: sec    !< seconds since model initialization
    real(kind=kind_phys), pointer :: si(:)  !< vertical sigma coordinate for model initialization
#endif

!--- IAU
    integer              :: iau_offset
    real(kind=kind_phys) :: iau_delthrs     ! iau time interval (to scale increments) in hours
    character(len=240)   :: iau_inc_files(7)! list of increment files
    real(kind=kind_phys) :: iaufhrs(7)      ! forecast hours associated with increment files
    logical :: iau_filter_increments

#ifdef CCPP
    ! From physcons.F90, updated/set in control_initialize
    real(kind=kind_phys) :: dxinv           ! inverse scaling factor for critical relative humidity, replaces dxinv in physcons.F90
    real(kind=kind_phys) :: dxmax           ! maximum scaling factor for critical relative humidity, replaces dxmax in physcons.F90
    real(kind=kind_phys) :: dxmin           ! minimum scaling factor for critical relative humidity, replaces dxmin in physcons.F90
    real(kind=kind_phys) :: rhcmax          ! maximum critical relative humidity, replaces rhc_max in physcons.F90
#endif

    contains
      procedure :: init  => control_initialize
      procedure :: print => control_print
  end type GFS_control_type


!--------------------------------------------------------------------
! GFS_grid_type
!   grid data needed for interpolations and length-scale calculations
!--------------------------------------------------------------------
!! \section arg_table_GFS_grid_type
!! \htmlinclude GFS_grid_type.html
!!
  type GFS_grid_type
 
    real (kind=kind_phys), pointer :: xlon   (:)    => null()   !< grid longitude in radians, ok for both 0->2pi
                                                                !! or -pi -> +pi ranges
    real (kind=kind_phys), pointer :: xlat   (:)    => null()   !< grid latitude in radians, default to pi/2 ->
                                                                !! -pi/2 range, otherwise adj in subr called 
    real (kind=kind_phys), pointer :: xlat_d (:)    => null()   !< grid latitude in degrees, default to 90 ->
                                                                !! -90 range, otherwise adj in subr called
    real (kind=kind_phys), pointer :: xlon_d (:)    => null()   !< grid longitude in degrees, default to 0 ->
                                                                !! 360 range, otherwise adj in subr called
    real (kind=kind_phys), pointer :: sinlat (:)    => null()   !< sine of the grids corresponding latitudes
    real (kind=kind_phys), pointer :: coslat (:)    => null()   !< cosine of the grids corresponding latitudes
    real (kind=kind_phys), pointer :: area   (:)    => null()   !< area of the grid cell
    real (kind=kind_phys), pointer :: dx     (:)    => null()   !< relative dx for the grid cell

!--- grid-related interpolation data for prognostic ozone
    real (kind=kind_phys), pointer :: ddy_o3    (:) => null()   !< interpolation     weight for ozone
    integer,               pointer :: jindx1_o3 (:) => null()   !< interpolation  low index for ozone
    integer,               pointer :: jindx2_o3 (:) => null()   !< interpolation high index for ozone

!--- grid-related interpolation data for stratosphere water
    real (kind=kind_phys), pointer :: ddy_h     (:) => null()   !< interpolation     weight for h2o
    integer,               pointer :: jindx1_h  (:) => null()   !< interpolation  low index for h2o
    integer,               pointer :: jindx2_h  (:) => null()   !< interpolation high index for h2o

!--- grid-related interpolation data for prognostic iccn
    real (kind=kind_phys), pointer :: ddy_ci    (:) => null()   !< interpolation     weight for iccn
    integer,               pointer :: jindx1_ci (:) => null()   !< interpolation  low index for iccn
    integer,               pointer :: jindx2_ci (:) => null()   !< interpolation high index for iccn
    real (kind=kind_phys), pointer :: ddx_ci    (:) => null()   !< interpolation     weight for iccn
    integer,               pointer :: iindx1_ci (:) => null()   !< interpolation  low index for iccn
    integer,               pointer :: iindx2_ci (:) => null()   !< interpolation high index for iccn

!--- grid-related interpolation data for prescribed aerosols
    real (kind=kind_phys), pointer :: ddy_aer    (:) => null()   !< interpolation     weight for iaerclm
    integer,               pointer :: jindx1_aer (:) => null()   !< interpolation  low index for iaerclm
    integer,               pointer :: jindx2_aer (:) => null()   !< interpolation high index for iaerclm
    real (kind=kind_phys), pointer :: ddx_aer    (:) => null()   !< interpolation     weight for iaerclm
    integer,               pointer :: iindx1_aer (:) => null()   !< interpolation  low index for iaerclm
    integer,               pointer :: iindx2_aer (:) => null()   !< interpolation high index for iaerclm
    contains
      procedure :: create   => grid_create   !<   allocate array data
  end type GFS_grid_type


!-----------------------------------------------
! GFS_tbd_type
!   data not yet assigned to a defined container
!-----------------------------------------------
!! \section arg_table_GFS_tbd_type
!! \htmlinclude GFS_tbd_type.html
!!
  type GFS_tbd_type

!--- radiation random seeds
    integer,               pointer :: icsdsw   (:)     => null()  !< (rad. only) auxiliary cloud control arrays passed to main
    integer,               pointer :: icsdlw   (:)     => null()  !< (rad. only) radiations. if isubcsw/isubclw (input to init)
                                                                  !< (rad. only) are set to 2, the arrays contains provided
                                                                  !< (rad. only) random seeds for sub-column clouds generators

!--- In
    real (kind=kind_phys), pointer :: ozpl     (:,:,:) => null()  !< ozone forcing data
    real (kind=kind_phys), pointer :: h2opl    (:,:,:) => null()  !< water forcing data
    real (kind=kind_phys), pointer :: in_nm    (:,:)   => null()  !< IN number concentration
    real (kind=kind_phys), pointer :: ccn_nm   (:,:)   => null()  !< CCN number concentration
    real (kind=kind_phys), pointer :: aer_nm   (:,:,:) => null()  !< GOCART aerosol climo

    !--- active when ((.not. newsas .or. cal_pre) .and. random_clds)
#ifdef CCPP
    integer,               pointer :: imap     (:)     => null()  !< map of local index ix to global index i for this block
    integer,               pointer :: jmap     (:)     => null()  !< map of local index ix to global index j for this block
#endif
    real (kind=kind_phys), pointer :: rann     (:,:)   => null()  !< random number array (0-1)

!--- In/Out
    real (kind=kind_phys), pointer :: acv      (:)     => null()  !< array containing accumulated convective clouds
    real (kind=kind_phys), pointer :: acvb     (:)     => null()  !< arrays used by cnvc90 bottom
    real (kind=kind_phys), pointer :: acvt     (:)     => null()  !< arrays used by cnvc90 top (cnvc90.f)

!--- Stochastic physics properties calculated in physics_driver
    real (kind=kind_phys), pointer :: dtdtr     (:,:)  => null()  !< temperature change due to radiative heating per time step (K)
    real (kind=kind_phys), pointer :: dtotprcp  (:)    => null()  !< change in totprcp  (diag_type)
    real (kind=kind_phys), pointer :: dcnvprcp  (:)    => null()  !< change in cnvprcp  (diag_type)
    real (kind=kind_phys), pointer :: drain_cpl (:)    => null()  !< change in rain_cpl (coupling_type)
    real (kind=kind_phys), pointer :: dsnow_cpl (:)    => null()  !< change in show_cpl (coupling_type)

!--- phy_f*d variables needed for seamless restarts and moving data between grrad and gbphys
    real (kind=kind_phys), pointer :: phy_fctd (:,:)   => null()  !< cloud base mass flux for CS convection
    real (kind=kind_phys), pointer :: phy_f2d  (:,:)   => null()  !< 2d arrays saved for restart
    real (kind=kind_phys), pointer :: phy_f3d  (:,:,:) => null()  !< 3d arrays saved for restart

!--- Diagnostic that needs to be carried over to the next time step (removed from diag_type)
    real (kind=kind_phys), pointer :: hpbl     (:)     => null()  !< Planetary boundary layer height

#ifndef CCPP
!--- for explicit data blocking
    integer                        :: blkno                       !< block number of this block
#endif

#ifdef CCPP
    !--- radiation variables that need to be carried over from radiation to physics
    real (kind=kind_phys), pointer :: htlwc(:,:)       => null()  !<
    real (kind=kind_phys), pointer :: htlw0(:,:)       => null()  !<
    real (kind=kind_phys), pointer :: htswc(:,:)       => null()  !<
    real (kind=kind_phys), pointer :: htsw0(:,:)       => null()  !<

    !--- dynamical forcing variables for Grell-Freitas convection
    real (kind=kind_phys), pointer :: forcet (:,:)     => null()  !<
    real (kind=kind_phys), pointer :: forceq (:,:)     => null()  !<
    real (kind=kind_phys), pointer :: prevst (:,:)     => null()  !<
    real (kind=kind_phys), pointer :: prevsq (:,:)     => null()  !<
    integer,               pointer :: cactiv   (:)     => null()  !< convective activity memory contour

    !--- MYNN prognostic variables that can't be in the Intdiag or Interstitial DDTs
    real (kind=kind_phys), pointer :: CLDFRA_BL  (:,:)   => null()  !
    real (kind=kind_phys), pointer :: QC_BL      (:,:)   => null()  !
    real (kind=kind_phys), pointer :: el_pbl     (:,:)   => null()  !
    real (kind=kind_phys), pointer :: Sh3D       (:,:)   => null()  !
    real (kind=kind_phys), pointer :: qke        (:,:)   => null()  !
    real (kind=kind_phys), pointer :: tsq        (:,:)   => null()  !
    real (kind=kind_phys), pointer :: qsq        (:,:)   => null()  !
    real (kind=kind_phys), pointer :: cov        (:,:)   => null()  !

    !--- MYJ schemes saved variables (from previous time step)
    real (kind=kind_phys), pointer :: phy_myj_qsfc(:)    => null()  ! 
    real (kind=kind_phys), pointer :: phy_myj_thz0(:)    => null()  ! 
    real (kind=kind_phys), pointer :: phy_myj_qz0(:)     => null()  ! 
    real (kind=kind_phys), pointer :: phy_myj_uz0(:)     => null()  ! 
    real (kind=kind_phys), pointer :: phy_myj_vz0(:)     => null()  ! 
    real (kind=kind_phys), pointer :: phy_myj_z0base(:)  => null()  ! 
    real (kind=kind_phys), pointer :: phy_myj_akhs(:)    => null()  ! 
    real (kind=kind_phys), pointer :: phy_myj_akms(:)    => null()  ! 
    real (kind=kind_phys), pointer :: phy_myj_chkqlm(:)  => null()  ! 
    real (kind=kind_phys), pointer :: phy_myj_elflx(:)   => null()  ! 
    real (kind=kind_phys), pointer :: phy_myj_a1u(:)     => null()  ! 
    real (kind=kind_phys), pointer :: phy_myj_a1t(:)     => null()  ! 
    real (kind=kind_phys), pointer :: phy_myj_a1q(:)     => null()  ! 
#endif

    contains
      procedure :: create  => tbd_create  !<   allocate array data
  end type GFS_tbd_type


!------------------------------------------------------------------
! GFS_cldprop_type
!  cloud properties and tendencies needed by radiation from physics 
!------------------------------------------------------------------
!! \section arg_table_GFS_cldprop_type
!! \htmlinclude GFS_cldprop_type.html
!!
  type GFS_cldprop_type

!--- In     (radiation)
!--- In/Out (physics)
    real (kind=kind_phys), pointer :: cv  (:)     => null()  !< fraction of convective cloud ; phys
    real (kind=kind_phys), pointer :: cvt (:)     => null()  !< convective cloud top pressure in pa ; phys
    real (kind=kind_phys), pointer :: cvb (:)     => null()  !< convective cloud bottom pressure in pa ; phys, cnvc90

    contains
      procedure :: create  => cldprop_create  !<   allocate array data
  end type GFS_cldprop_type


!-----------------------------------------
! GFS_radtend_type
!   radiation tendencies needed by physics
!-----------------------------------------
!! \section arg_table_GFS_radtend_type
!! \htmlinclude GFS_radtend_type.html
!!
  type GFS_radtend_type

    type (sfcfsw_type),    pointer :: sfcfsw(:)   => null()   !< sw radiation fluxes at sfc
                                                              !< [dim(im): created in grrad.f], components:
                                                              !!     (check module_radsw_parameters for definition)
                                                              !!\n   %upfxc - total sky upward sw flux at sfc (w/m**2)
                                                              !!\n   %upfx0 - clear sky upward sw flux at sfc (w/m**2)
                                                              !!\n   %dnfxc - total sky downward sw flux at sfc (w/m**2)
                                                              !!\n   %dnfx0 - clear sky downward sw flux at sfc (w/m**2)

    type (sfcflw_type),    pointer :: sfcflw(:)    => null()  !< lw radiation fluxes at sfc
                                                              !< [dim(im): created in grrad.f], components:
                                                              !!     (check module_radlw_paramters for definition)
                                                              !!\n   %upfxc - total sky upward lw flux at sfc (w/m**2)
                                                              !!\n   %upfx0 - clear sky upward lw flux at sfc (w/m**2)
                                                              !!\n   %dnfxc - total sky downward lw flux at sfc (w/m**2)
                                                              !!\n   %dnfx0 - clear sky downward lw flux at sfc (w/m**2)

!--- Out (radiation only)
    real (kind=kind_phys), pointer :: htrsw (:,:)  => null()  !< swh  total sky sw heating rate in k/sec
    real (kind=kind_phys), pointer :: htrlw (:,:)  => null()  !< hlw  total sky lw heating rate in k/sec
    real (kind=kind_phys), pointer :: sfalb (:)    => null()  !< mean surface diffused sw albedo 

    real (kind=kind_phys), pointer :: coszen(:)    => null()  !< mean cos of zenith angle over rad call period
    real (kind=kind_phys), pointer :: tsflw (:)    => null()  !< surface air temp during lw calculation in k 
    real (kind=kind_phys), pointer :: semis (:)    => null()  !< surface lw emissivity in fraction

!--- In/Out (???) (radiaition only)
    real (kind=kind_phys), pointer :: coszdg(:)    => null()  !< daytime mean cosz over rad call period

!--- In/Out (???) (physics only)
    real (kind=kind_phys), pointer :: swhc (:,:)   => null()  !< clear sky sw heating rates ( k/s )
    real (kind=kind_phys), pointer :: lwhc (:,:)   => null()  !< clear sky lw heating rates ( k/s )
    real (kind=kind_phys), pointer :: lwhd (:,:,:) => null()  !< idea sky lw heating rates ( k/s )

    contains
      procedure :: create  => radtend_create   !<   allocate array data
  end type GFS_radtend_type

!----------------------------------------------------------------
! GFS_diag_type
!  internal diagnostic type used as arguments to gbphys and grrad 
!----------------------------------------------------------------
!! \section arg_table_GFS_diag_type
!! \htmlinclude GFS_diag_type.html
!!
  type GFS_diag_type

!! Input/Output only in radiation
    real (kind=kind_phys), pointer :: fluxr(:,:)     => null()   !< to save time accumulated 2-d fields defined as:!
                                                                 !< hardcoded field indices, opt. includes aerosols!
    type (topfsw_type),    pointer :: topfsw(:)      => null()   !< sw radiation fluxes at toa, components:
                                               !       %upfxc    - total sky upward sw flux at toa (w/m**2)
                                               !       %dnfxc    - total sky downward sw flux at toa (w/m**2)
                                               !       %upfx0    - clear sky upward sw flux at toa (w/m**2)
    type (topflw_type),    pointer :: topflw(:)      => null()   !< lw radiation fluxes at top, component:
                                               !        %upfxc    - total sky upward lw flux at toa (w/m**2)
                                               !        %upfx0    - clear sky upward lw flux at toa (w/m**2)

! Input/output - used by physics
    real (kind=kind_phys), pointer :: srunoff(:)     => null()   !< surface water runoff (from lsm)
    real (kind=kind_phys), pointer :: evbsa  (:)     => null()   !< noah lsm diagnostics
    real (kind=kind_phys), pointer :: evcwa  (:)     => null()   !< noah lsm diagnostics
    real (kind=kind_phys), pointer :: snohfa (:)     => null()   !< noah lsm diagnostics
    real (kind=kind_phys), pointer :: transa (:)     => null()   !< noah lsm diagnostics
    real (kind=kind_phys), pointer :: sbsnoa (:)     => null()   !< noah lsm diagnostics
    real (kind=kind_phys), pointer :: snowca (:)     => null()   !< noah lsm diagnostics
    real (kind=kind_phys), pointer :: soilm  (:)     => null()   !< soil moisture
    real (kind=kind_phys), pointer :: tmpmin (:)     => null()   !< min temperature at 2m height (k)
    real (kind=kind_phys), pointer :: tmpmax (:)     => null()   !< max temperature at 2m height (k)
    real (kind=kind_phys), pointer :: dusfc  (:)     => null()   !< u component of surface stress
    real (kind=kind_phys), pointer :: dvsfc  (:)     => null()   !< v component of surface stress
    real (kind=kind_phys), pointer :: dtsfc  (:)     => null()   !< sensible heat flux (w/m2)
    real (kind=kind_phys), pointer :: dqsfc  (:)     => null()   !< latent heat flux (w/m2)
    real (kind=kind_phys), pointer :: totprcp(:)     => null()   !< accumulated total precipitation (kg/m2)
    real (kind=kind_phys), pointer :: totprcpb(:)    => null()   !< accumulated total precipitation in bucket(kg/m2)
    real (kind=kind_phys), pointer :: gflux  (:)     => null()   !< groud conductive heat flux
    real (kind=kind_phys), pointer :: dlwsfc (:)     => null()   !< time accumulated sfc dn lw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: ulwsfc (:)     => null()   !< time accumulated sfc up lw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: suntim (:)     => null()   !< sunshine duration time (s)
    real (kind=kind_phys), pointer :: runoff (:)     => null()   !< total water runoff
    real (kind=kind_phys), pointer :: ep     (:)     => null()   !< potential evaporation
    real (kind=kind_phys), pointer :: cldwrk (:)     => null()   !< cloud workfunction (valid only with sas)
    real (kind=kind_phys), pointer :: dugwd  (:)     => null()   !< vertically integrated u change by OGWD
    real (kind=kind_phys), pointer :: dvgwd  (:)     => null()   !< vertically integrated v change by OGWD
    real (kind=kind_phys), pointer :: psmean (:)     => null()   !< surface pressure (kPa)
    real (kind=kind_phys), pointer :: cnvprcp(:)     => null()   !< accumulated convective precipitation (kg/m2)
    real (kind=kind_phys), pointer :: cnvprcpb(:)    => null()   !< accumulated convective precipitation in bucket (kg/m2)
    real (kind=kind_phys), pointer :: spfhmin(:)     => null()   !< minimum specific humidity
    real (kind=kind_phys), pointer :: spfhmax(:)     => null()   !< maximum specific humidity
    real (kind=kind_phys), pointer :: u10mmax(:)     => null()   !< maximum u-wind
    real (kind=kind_phys), pointer :: v10mmax(:)     => null()   !< maximum v-wind
    real (kind=kind_phys), pointer :: wind10mmax(:)  => null()   !< maximum wind speed
    real (kind=kind_phys), pointer :: u10max(:)      => null()   !< maximum u-wind used with avg_max_length
    real (kind=kind_phys), pointer :: v10max(:)      => null()   !< maximum v-wind used with avg_max_length
    real (kind=kind_phys), pointer :: spd10max(:)    => null()   !< maximum wind speed used with avg_max_length
    real (kind=kind_phys), pointer :: rain   (:)     => null()   !< total rain at this time step
    real (kind=kind_phys), pointer :: rainc  (:)     => null()   !< convective rain at this time step
    real (kind=kind_phys), pointer :: ice    (:)     => null()   !< ice fall at this time step
    real (kind=kind_phys), pointer :: snow   (:)     => null()   !< snow fall at this time step
    real (kind=kind_phys), pointer :: graupel(:)     => null()   !< graupel fall at this time step
    real (kind=kind_phys), pointer :: totice (:)     => null()   !< accumulated ice precipitation (kg/m2)
    real (kind=kind_phys), pointer :: totsnw (:)     => null()   !< accumulated snow precipitation (kg/m2)
    real (kind=kind_phys), pointer :: totgrp (:)     => null()   !< accumulated graupel precipitation (kg/m2)
    real (kind=kind_phys), pointer :: toticeb(:)     => null()   !< accumulated ice precipitation in bucket (kg/m2)
    real (kind=kind_phys), pointer :: totsnwb(:)     => null()   !< accumulated snow precipitation in bucket (kg/m2)
    real (kind=kind_phys), pointer :: totgrpb(:)     => null()   !< accumulated graupel precipitation in bucket (kg/m2)

#ifdef CCPP
    !--- MYNN variables                                              
    real (kind=kind_phys), pointer :: edmf_a     (:,:)   => null()  !
    real (kind=kind_phys), pointer :: edmf_w     (:,:)   => null()  !
    real (kind=kind_phys), pointer :: edmf_qt    (:,:)   => null()  !
    real (kind=kind_phys), pointer :: edmf_thl   (:,:)   => null()  !
    real (kind=kind_phys), pointer :: edmf_ent   (:,:)   => null()  !
    real (kind=kind_phys), pointer :: edmf_qc    (:,:)   => null()  !
    real (kind=kind_phys), pointer :: maxMF       (:)    => null()  !
    integer, pointer               :: nupdraft    (:)    => null()  !
    integer, pointer               :: ktop_shallow (:)   => null()  !
    real (kind=kind_phys), pointer :: exch_h     (:,:)   => null()  !
    real (kind=kind_phys), pointer :: exch_m     (:,:)   => null()  !

    !--- Drag Suite variables
    real (kind=kind_phys), pointer :: dtaux2d_ls  (:,:)   => null()  !
    real (kind=kind_phys), pointer :: dtauy2d_ls  (:,:)   => null()  !
    real (kind=kind_phys), pointer :: dtaux2d_bl  (:,:)   => null()  !
    real (kind=kind_phys), pointer :: dtauy2d_bl  (:,:)   => null()  !
    real (kind=kind_phys), pointer :: dtaux2d_ss  (:,:)   => null()  !
    real (kind=kind_phys), pointer :: dtauy2d_ss  (:,:)   => null()  !
    real (kind=kind_phys), pointer :: dtaux2d_fd  (:,:)   => null()  !
    real (kind=kind_phys), pointer :: dtauy2d_fd  (:,:)   => null()  !
    real (kind=kind_phys), pointer :: dusfc_ls    (:)     => null()  !
    real (kind=kind_phys), pointer :: dvsfc_ls    (:)     => null()  !
    real (kind=kind_phys), pointer :: dusfc_bl    (:)     => null()  !
    real (kind=kind_phys), pointer :: dvsfc_bl    (:)     => null()  !
    real (kind=kind_phys), pointer :: dusfc_ss    (:)     => null()  !
    real (kind=kind_phys), pointer :: dvsfc_ss    (:)     => null()  !
    real (kind=kind_phys), pointer :: dusfc_fd    (:)     => null()  !
    real (kind=kind_phys), pointer :: dvsfc_fd    (:)     => null()  !
#endif

! Output - only in physics
    real (kind=kind_phys), pointer :: u10m   (:)     => null()   !< 10 meter u/v wind speed
    real (kind=kind_phys), pointer :: v10m   (:)     => null()   !< 10 meter u/v wind speed
    real (kind=kind_phys), pointer :: dpt2m  (:)     => null()   !< 2 meter dew point temperature
    real (kind=kind_phys), pointer :: zlvl   (:)     => null()   !< layer 1 height (m)
    real (kind=kind_phys), pointer :: psurf  (:)     => null()   !< surface pressure (Pa)
    real (kind=kind_phys), pointer :: pwat   (:)     => null()   !< precipitable water
    real (kind=kind_phys), pointer :: t1     (:)     => null()   !< layer 1 temperature (K)
    real (kind=kind_phys), pointer :: q1     (:)     => null()   !< layer 1 specific humidity (kg/kg)
    real (kind=kind_phys), pointer :: u1     (:)     => null()   !< layer 1 zonal wind (m/s)
    real (kind=kind_phys), pointer :: v1     (:)     => null()   !< layer 1 merdional wind (m/s)
    real (kind=kind_phys), pointer :: chh    (:)     => null()   !< thermal exchange coefficient
    real (kind=kind_phys), pointer :: cmm    (:)     => null()   !< momentum exchange coefficient
    real (kind=kind_phys), pointer :: dlwsfci(:)     => null()   !< instantaneous sfc dnwd lw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: ulwsfci(:)     => null()   !< instantaneous sfc upwd lw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: dswsfci(:)     => null()   !< instantaneous sfc dnwd sw flux ( w/m**2 )
#ifdef CCPP
    real (kind=kind_phys), pointer :: nswsfci(:)     => null()   !< instantaneous sfc net dnwd sw flux ( w/m**2 )
#endif
    real (kind=kind_phys), pointer :: uswsfci(:)     => null()   !< instantaneous sfc upwd sw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: dusfci (:)     => null()   !< instantaneous u component of surface stress
    real (kind=kind_phys), pointer :: dvsfci (:)     => null()   !< instantaneous v component of surface stress
    real (kind=kind_phys), pointer :: dtsfci (:)     => null()   !< instantaneous sfc sensible heat flux
    real (kind=kind_phys), pointer :: dqsfci (:)     => null()   !< instantaneous sfc latent heat flux
    real (kind=kind_phys), pointer :: gfluxi (:)     => null()   !< instantaneous sfc ground heat flux
    real (kind=kind_phys), pointer :: epi    (:)     => null()   !< instantaneous sfc potential evaporation
    real (kind=kind_phys), pointer :: smcwlt2(:)     => null()   !< wilting point (volumetric)
    real (kind=kind_phys), pointer :: smcref2(:)     => null()   !< soil moisture threshold (volumetric)
    real (kind=kind_phys), pointer :: wet1   (:)     => null()   !< normalized soil wetness
    real (kind=kind_phys), pointer :: sr     (:)     => null()   !< snow ratio : ratio of snow to total precipitation
    real (kind=kind_phys), pointer :: tdomr  (:)     => null()   !< dominant accumulated rain type
    real (kind=kind_phys), pointer :: tdomzr (:)     => null()   !< dominant accumulated freezing rain type
    real (kind=kind_phys), pointer :: tdomip (:)     => null()   !< dominant accumulated sleet type
    real (kind=kind_phys), pointer :: tdoms  (:)     => null()   !< dominant accumulated snow type

    real (kind=kind_phys), pointer :: ca_out  (:)    => null()   !< cellular automata fraction
    real (kind=kind_phys), pointer :: ca_deep  (:)   => null()   !< cellular automata fraction
    real (kind=kind_phys), pointer :: ca_turb  (:)   => null()   !< cellular automata fraction
    real (kind=kind_phys), pointer :: ca_shal  (:)   => null()   !< cellular automata fraction
    real (kind=kind_phys), pointer :: ca_rad   (:)   => null()   !< cellular automata fraction
    real (kind=kind_phys), pointer :: ca_micro (:)   => null()   !< cellular automata fraction

    real (kind=kind_phys), pointer :: skebu_wts(:,:) => null()   !< 10 meter u wind speed
    real (kind=kind_phys), pointer :: skebv_wts(:,:) => null()   !< 10 meter v wind speed
    real (kind=kind_phys), pointer :: sppt_wts(:,:)  => null()   !<
    real (kind=kind_phys), pointer :: shum_wts(:,:)  => null()   !<
    real (kind=kind_phys), pointer :: zmtnblck(:)    => null()   !<mountain blocking evel
    real (kind=kind_phys), pointer :: du3dt (:,:,:)  => null()   !< u momentum change due to physics
    real (kind=kind_phys), pointer :: dv3dt (:,:,:)  => null()   !< v momentum change due to physics
    real (kind=kind_phys), pointer :: dt3dt (:,:,:)  => null()   !< temperature change due to physics
    real (kind=kind_phys), pointer :: dq3dt (:,:,:)  => null()   !< moisture change due to physics
    real (kind=kind_phys), pointer :: refdmax (:)    => null()   !< max hourly 1-km agl reflectivity
    real (kind=kind_phys), pointer :: refdmax263k(:) => null()   !< max hourly -10C reflectivity
    real (kind=kind_phys), pointer :: t02max  (:)    => null()   !< max hourly 2m T
    real (kind=kind_phys), pointer :: t02min  (:)    => null()   !< min hourly 2m T
    real (kind=kind_phys), pointer :: rh02max (:)    => null()   !< max hourly 2m RH
    real (kind=kind_phys), pointer :: rh02min (:)    => null()   !< min hourly 2m RH
!--- accumulated quantities for 3D diagnostics
    real (kind=kind_phys), pointer :: upd_mf (:,:)   => null()  !< instantaneous convective updraft mass flux
    real (kind=kind_phys), pointer :: dwn_mf (:,:)   => null()  !< instantaneous convective downdraft mass flux
    real (kind=kind_phys), pointer :: det_mf (:,:)   => null()  !< instantaneous convective detrainment mass flux
    real (kind=kind_phys), pointer :: cldcov (:,:)   => null()  !< instantaneous 3D cloud fraction
!--- F-A MP scheme
#ifdef CCPP
    real (kind=kind_phys), pointer :: TRAIN  (:,:)   => null()  !< accumulated stratiform T tendency (K s-1)
#endif

    !--- MP quantities for 3D diagnositics 
    real (kind=kind_phys), pointer :: refl_10cm(:,:) => null()  !< instantaneous refl_10cm 
!
!---vay-2018 UGWP-diagnostics daily mean
!
    real (kind=kind_phys), pointer :: dudt_tot (:,:) => null()  !< daily aver GFS_phys tend for WE-U
    real (kind=kind_phys), pointer :: dvdt_tot (:,:) => null()  !< daily aver GFS_phys tend for SN-V
    real (kind=kind_phys), pointer :: dtdt_tot (:,:) => null()  !< daily aver GFS_phys tend for Temp-re
!
    real (kind=kind_phys), pointer :: du3dt_pbl(:,:) => null()  !< daily aver GFS_phys tend for WE-U pbl
    real (kind=kind_phys), pointer :: dv3dt_pbl(:,:) => null()  !< daily aver GFS_phys tend for SN-V pbl
    real (kind=kind_phys), pointer :: dt3dt_pbl(:,:) => null()  !< daily aver GFS_phys tend for Temp pbl
!
    real (kind=kind_phys), pointer :: du3dt_ogw(:,:) => null()  !< daily aver GFS_phys tend for WE-U OGW
    real (kind=kind_phys), pointer :: dv3dt_ogw(:,:) => null()  !< daily aver GFS_phys tend for SN-V OGW
    real (kind=kind_phys), pointer :: dt3dt_ogw(:,:) => null()  !< daily aver GFS_phys tend for Temp OGW
!
    real (kind=kind_phys), pointer :: du3dt_mtb(:,:) => null()  !< daily aver GFS_phys tend for WE-U MTB
    real (kind=kind_phys), pointer :: dv3dt_mtb(:,:) => null()  !< daily aver GFS_phys tend for SN-V MTB
    real (kind=kind_phys), pointer :: dt3dt_mtb(:,:) => null()  !< daily aver GFS_phys tend for Temp MTB
!
    real (kind=kind_phys), pointer :: du3dt_tms(:,:) => null()  !< daily aver GFS_phys tend for WE-U TMS
    real (kind=kind_phys), pointer :: dv3dt_tms(:,:) => null()  !< daily aver GFS_phys tend for SN-V TMS
    real (kind=kind_phys), pointer :: dt3dt_tms(:,:) => null()  !< daily aver GFS_phys tend for Temp TMS
!
    real (kind=kind_phys), pointer :: du3dt_ngw(:,:) => null()  !< daily aver GFS_phys tend for WE-U NGW
    real (kind=kind_phys), pointer :: dv3dt_ngw(:,:) => null()  !< daily aver GFS_phys tend for SN-V NGW
    real (kind=kind_phys), pointer :: dt3dt_ngw(:,:) => null()  !< daily aver GFS_phys tend for Temp NGW
!
    real (kind=kind_phys), pointer :: du3dt_cgw(:,:) => null()  !< daily aver GFS_phys tend for WE-U NGW
    real (kind=kind_phys), pointer :: dv3dt_cgw(:,:) => null()  !< daily aver GFS_phys tend for SN-V NGW
    real (kind=kind_phys), pointer :: dt3dt_cgw(:,:) => null()  !< daily aver GFS_phys tend for Temp NGW
!
    real (kind=kind_phys), pointer :: du3dt_moist(:,:) => null()  !< daily aver GFS_phys tend for WE-U MOIST
    real (kind=kind_phys), pointer :: dv3dt_moist(:,:) => null()  !< daily aver GFS_phys tend for SN-V MOIST
    real (kind=kind_phys), pointer :: dt3dt_moist(:,:) => null()  !< daily aver GFS_phys tend for Temp MOIST
!
!--- Instantaneous UGWP-diagnostics  16-variables
!       Diag%gwp_ax, Diag%gwp_axo, Diag%gwp_axc, Diag%gwp_axf,       &
!       Diag%gwp_ay, Diag%gwp_ayo, Diag%gwp_ayc, Diag%gwp_ayf,       &
!       Diag%gwp_dtdt,   Diag%gwp_kdis, Diag%gwp_okw, Diag%gwp_fgf,  &
!       Diag%gwp_dcheat, Diag%gwp_precip, Diag%gwp_klevs,            &
!       Diag%gwp_scheat

    real (kind=kind_phys), pointer :: gwp_scheat(:,:) => null()  ! instant shal-conv heat tendency
    real (kind=kind_phys), pointer :: gwp_dcheat(:,:) => null()  ! instant deep-conv heat tendency
    real (kind=kind_phys), pointer :: gwp_precip(:) => null()    ! total precip rates
    integer , pointer              :: gwp_klevs(:,:)=> null()    ! instant levels for GW-launches
    real (kind=kind_phys), pointer :: gwp_fgf(:)    => null()    ! fgf triggers
    real (kind=kind_phys), pointer :: gwp_okw(:)    => null()    ! okw triggers

    real (kind=kind_phys), pointer :: gwp_ax(:,:)   => null()    ! instant total UGWP tend m/s/s EW
    real (kind=kind_phys), pointer :: gwp_ay(:,:)   => null()    ! instant total UGWP tend m/s/s NS
    real (kind=kind_phys), pointer :: gwp_dtdt(:,:) => null()    ! instant total heat tend   K/s
    real (kind=kind_phys), pointer :: gwp_kdis(:,:) => null()    ! instant total eddy mixing m2/s
    real (kind=kind_phys), pointer :: gwp_axc(:,:)   => null()   ! instant con-UGWP tend m/s/s EW
    real (kind=kind_phys), pointer :: gwp_ayc(:,:)   => null()   ! instant con-UGWP tend m/s/s NS
    real (kind=kind_phys), pointer :: gwp_axo(:,:)   => null()   ! instant oro-UGWP tend m/s/s EW
    real (kind=kind_phys), pointer :: gwp_ayo(:,:)   => null()   ! instant oro-UGWP tend m/s/s NS
    real (kind=kind_phys), pointer :: gwp_axf(:,:)   => null()   ! instant jet-UGWP tend m/s/s EW
    real (kind=kind_phys), pointer :: gwp_ayf(:,:)   => null()   ! instant jet-UGWP tend m/s/s NS

    real (kind=kind_phys), pointer :: uav_ugwp(:,:)   => null()   ! aver  wind UAV from physics
    real (kind=kind_phys), pointer :: tav_ugwp(:,:)   => null()   ! aver  temp UAV from physics
    real (kind=kind_phys), pointer :: du3dt_dyn(:,:)  => null()   ! U Tend-dynamics "In"-"PhysOut"

!--- COODRE ORO diagnostics
    real (kind=kind_phys), pointer :: zmtb(:)         => null()   !
    real (kind=kind_phys), pointer :: zogw(:)         => null()   !
    real (kind=kind_phys), pointer :: zlwb(:)         => null()   !!
    real (kind=kind_phys), pointer :: tau_ogw(:)      => null()   !!
    real (kind=kind_phys), pointer :: tau_ngw(:)      => null()   !!
    real (kind=kind_phys), pointer :: tau_mtb(:)      => null()   !
    real (kind=kind_phys), pointer :: tau_tofd(:)     => null()   !
!---vay-2018 UGWP-diagnostics

    !--- Output diagnostics for coupled chemistry
#ifdef CCPP
    integer                        :: ndust                    !< number of dust bins for diagnostics
    integer                        :: nseasalt                 !< number of seasalt bins for diagnostics
    integer                        :: ntchmdiag                !< number of chemical tracers for diagnostics
#endif
    real (kind=kind_phys), pointer :: duem  (:,:) => null()    !< instantaneous dust emission flux                             ( kg/m**2/s )
    real (kind=kind_phys), pointer :: ssem  (:,:) => null()    !< instantaneous sea salt emission flux                         ( kg/m**2/s )
    real (kind=kind_phys), pointer :: sedim (:,:) => null()    !< instantaneous sedimentation                                  ( kg/m**2/s )
    real (kind=kind_phys), pointer :: drydep(:,:) => null()    !< instantaneous dry deposition                                 ( kg/m**2/s )
    real (kind=kind_phys), pointer :: wetdpl(:,:) => null()    !< instantaneous large-scale wet deposition                     ( kg/m**2/s )
    real (kind=kind_phys), pointer :: wetdpc(:,:) => null()    !< instantaneous convective-scale wet deposition                ( kg/m**2/s )
    real (kind=kind_phys), pointer :: abem  (:,:) => null()    !< instantaneous anthopogenic and biomass burning emissions
                                                               !< for black carbon, organic carbon, and sulfur dioxide         ( ug/m**2/s )
    real (kind=kind_phys), pointer :: aecm  (:,:) => null()    !< instantaneous aerosol column mass densities for
                                                               !< pm2.5, black carbon, organic carbon, sulfate, dust, sea salt ( g/m**2 )
    contains
      procedure :: create    => diag_create
      procedure :: rad_zero  => diag_rad_zero
      procedure :: phys_zero => diag_phys_zero
      procedure :: chem_init => diag_chem_init
  end type GFS_diag_type

#ifdef CCPP
!---------------------------------------------------------------------
! GFS_interstitial_type
!   fields required for interstitial code in CCPP schemes, previously
!   in GFS_{physics,radiation}_driver.F90
!---------------------------------------------------------------------
!! \section arg_table_GFS_interstitial_type
!! \htmlinclude GFS_interstitial_type.html
!!
  type GFS_interstitial_type

    real (kind=kind_phys), pointer      :: adjsfculw_land(:)  => null()  !<
    real (kind=kind_phys), pointer      :: adjsfculw_ice(:)   => null()  !<
    real (kind=kind_phys), pointer      :: adjsfculw_ocean(:) => null()  !<
    real (kind=kind_phys), pointer      :: adjnirbmd(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjnirbmu(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjnirdfd(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjnirdfu(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjvisbmd(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjvisbmu(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjvisdfu(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjvisdfd(:)       => null()  !<
    real (kind=kind_phys), pointer      :: aerodp(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: alb1d(:)           => null()  !<
    real (kind=kind_phys), pointer      :: bexp1d(:)          => null()  !<
    real (kind=kind_phys), pointer      :: cd(:)              => null()  !<
    real (kind=kind_phys), pointer      :: cd_ice(:)          => null()  !<
    real (kind=kind_phys), pointer      :: cd_land(:)         => null()  !<
    real (kind=kind_phys), pointer      :: cd_ocean(:)        => null()  !<
    real (kind=kind_phys), pointer      :: cdq(:)             => null()  !<
    real (kind=kind_phys), pointer      :: cdq_ice(:)         => null()  !<
    real (kind=kind_phys), pointer      :: cdq_land(:)        => null()  !<
    real (kind=kind_phys), pointer      :: cdq_ocean(:)       => null()  !<
    real (kind=kind_phys), pointer      :: cf_upi(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: chh_ice(:)         => null()  !<
    real (kind=kind_phys), pointer      :: chh_land(:)        => null()  !<
    real (kind=kind_phys), pointer      :: chh_ocean(:)       => null()  !<
    real (kind=kind_phys), pointer      :: clcn(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: cldf(:)            => null()  !<
    real (kind=kind_phys), pointer      :: cldsa(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: cldtaulw(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: cldtausw(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: cld1d(:)           => null()  !<
    real (kind=kind_phys), pointer      :: clouds(:,:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: clw(:,:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: clw_surf(:)        => null()  !<
    real (kind=kind_phys), pointer      :: clx(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: cmm_ice(:)         => null()  !<
    real (kind=kind_phys), pointer      :: cmm_land(:)        => null()  !<
    real (kind=kind_phys), pointer      :: cmm_ocean(:)       => null()  !<
    real (kind=kind_phys), pointer      :: cndm_surf(:)       => null()  !<
    real (kind=kind_phys), pointer      :: cnv_dqldt(:,:)     => null()  !<
    real (kind=kind_phys), pointer      :: cnv_fice(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: cnv_mfd(:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: cnv_ndrop(:,:)     => null()  !<
    real (kind=kind_phys), pointer      :: cnv_nice(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: cnvc(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: cnvw(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: ctei_r(:)          => null()  !<
    real (kind=kind_phys), pointer      :: ctei_rml(:)        => null()  !<
    real (kind=kind_phys), pointer      :: cumabs(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dd_mf(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: de_lgth(:)         => null()  !<
    real (kind=kind_phys), pointer      :: del(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: del_gz(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: delr(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: dkt(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: dlength(:)         => null()  !<
    real (kind=kind_phys), pointer      :: dqdt(:,:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: dqsfc1(:)          => null()  !<
    real (kind=kind_phys), pointer      :: drain(:)           => null()  !<
    real (kind=kind_phys), pointer      :: dtdt(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: dtdtc(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: dtsfc1(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dtzm(:)            => null()  !<
    real (kind=kind_phys), pointer      :: dt_mf(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: dudt(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: dusfcg(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dusfc1(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dvdftra(:,:,:)     => null()  !<
    real (kind=kind_phys), pointer      :: dvdt(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: dvsfcg(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dvsfc1(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dzlyr(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: elvmax(:)          => null()  !<
    real (kind=kind_phys), pointer      :: ep1d(:)            => null()  !<
    real (kind=kind_phys), pointer      :: ep1d_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: ep1d_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: ep1d_ocean(:)      => null()  !<
    real (kind=kind_phys), pointer      :: evap_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: evap_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: evap_ocean(:)      => null()  !<
    real (kind=kind_phys), pointer      :: evbs(:)            => null()  !<
    real (kind=kind_phys), pointer      :: evcw(:)            => null()  !<
    real (kind=kind_phys), pointer      :: faerlw(:,:,:,:)    => null()  !<
    real (kind=kind_phys), pointer      :: faersw(:,:,:,:)    => null()  !<
    real (kind=kind_phys), pointer      :: ffhh_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: ffhh_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: ffhh_ocean(:)      => null()  !<
    real (kind=kind_phys), pointer      :: fh2(:)             => null()  !<
    real (kind=kind_phys), pointer      :: fh2_ice(:)         => null()  !<
    real (kind=kind_phys), pointer      :: fh2_land(:)        => null()  !<
    real (kind=kind_phys), pointer      :: fh2_ocean(:)       => null()  !<
    logical,               pointer      :: flag_cice(:)       => null()  !<
    logical,               pointer      :: flag_guess(:)      => null()  !<
    logical,               pointer      :: flag_iter(:)       => null()  !<
    real (kind=kind_phys), pointer      :: ffmm_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: ffmm_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: ffmm_ocean(:)      => null()  !<
    real (kind=kind_phys), pointer      :: fm10(:)            => null()  !<
    real (kind=kind_phys), pointer      :: fm10_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: fm10_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: fm10_ocean(:)      => null()  !<
    real (kind=kind_phys)               :: frain                         !<
    real (kind=kind_phys), pointer      :: frland(:)          => null()  !<
    real (kind=kind_phys), pointer      :: fscav(:)           => null()  !<
    real (kind=kind_phys), pointer      :: fswtr(:)           => null()  !<
    real (kind=kind_phys), pointer      :: gabsbdlw(:)        => null()  !<
    real (kind=kind_phys), pointer      :: gabsbdlw_ice(:)    => null()  !<
    real (kind=kind_phys), pointer      :: gabsbdlw_land(:)   => null()  !<
    real (kind=kind_phys), pointer      :: gabsbdlw_ocean(:)  => null()  !<
    real (kind=kind_phys), pointer      :: gamma(:)           => null()  !<
    real (kind=kind_phys), pointer      :: gamq(:)            => null()  !<
    real (kind=kind_phys), pointer      :: gamt(:)            => null()  !<
    real (kind=kind_phys), pointer      :: gasvmr(:,:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: gflx(:)            => null()  !<
    real (kind=kind_phys), pointer      :: gflx_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: gflx_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: gflx_ocean(:)      => null()  !<
    real (kind=kind_phys), pointer      :: graupelmp(:)       => null()  !<
    real (kind=kind_phys), pointer      :: gwdcu(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: gwdcv(:,:)         => null()  !<
    integer                             :: h2o_coeff                     !<
    real (kind=kind_phys), pointer      :: h2o_pres(:)        => null()  !<
    real (kind=kind_phys), pointer      :: hflx_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: hflx_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: hflx_ocean(:)      => null()  !<
    real (kind=kind_phys), pointer      :: icemp(:)           => null()  !<
    logical,               pointer      :: dry(:)             => null()  !<
    integer,               pointer      :: idxday(:)          => null()  !<
    logical,               pointer      :: icy(:)             => null()  !<
    logical,               pointer      :: lake(:)            => null()  !<
    logical,               pointer      :: ocean(:)           => null()  !<
    integer                             :: ipr                           !<
    integer,               pointer      :: islmsk(:)          => null()  !<
    integer,               pointer      :: islmsk_cice(:)     => null()  !<
    integer                             :: itc                           !<
    logical,               pointer      :: wet(:)             => null()  !<
    integer                             :: kb                            !<
    integer,               pointer      :: kbot(:)            => null()  !<
    integer,               pointer      :: kcnv(:)            => null()  !<
    integer                             :: kd                            !<
    integer,               pointer      :: kinver(:)          => null()  !<
    integer,               pointer      :: kpbl(:)            => null()  !<
    integer                             :: kt                            !<
    integer,               pointer      :: ktop(:)            => null()  !<
    integer                             :: latidxprnt                    !<
    integer                             :: levi                          !<
    integer                             :: levh2o                        !<
    integer                             :: levozp                        !<
    integer                             :: lmk                           !<
    integer                             :: lmp                           !<
    integer,               pointer      :: mbota(:,:)         => null()  !<
    logical                             :: mg3_as_mg2                    !<
    integer,               pointer      :: mtopa(:,:)         => null()  !<
    integer                             :: nbdlw                         !<
    integer                             :: nbdsw                         !<
    real (kind=kind_phys), pointer      :: ncgl(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: ncpi(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: ncpl(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: ncpr(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: ncps(:,:)          => null()  !<
    integer                             :: ncstrac                       !<
    integer                             :: nday                          !<
    integer                             :: nf_aelw                       !<
    integer                             :: nf_aesw                       !<
    integer                             :: nn                            !<
    integer                             :: nncl                          !<
    integer                             :: nsamftrac                     !<
    integer                             :: nscav                         !<
    integer                             :: nspc1                         !<
    integer                             :: ntiwx                         !<
    integer                             :: ntk                           !<
    integer                             :: ntkev                         !<
    integer                             :: nvdiff                        !<
    real (kind=kind_phys), pointer      :: oa4(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: oc(:)              => null()  !<
    real (kind=kind_phys), pointer      :: olyr(:,:)          => null()  !<
    logical              , pointer      :: otspt(:,:)         => null()  !<
    integer                             :: oz_coeff                      !<
    integer                             :: oz_coeffp5                    !<
    real (kind=kind_phys), pointer      :: oz_pres(:)         => null()  !<
    logical                             :: phys_hydrostatic              !<
    real (kind=kind_phys), pointer      :: plvl(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: plyr(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: prcpmp(:)          => null()  !<
    real (kind=kind_phys), pointer      :: prnum(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: q2mp(:)            => null()  !<
    real (kind=kind_phys), pointer      :: qgl(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: qicn(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: qlcn(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: qlyr(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: qrn(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: qsnw(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: qss(:)             => null()  !<
    real (kind=kind_phys)               :: raddt                         !<
    real (kind=kind_phys), pointer      :: rainmp(:)          => null()  !<
    real (kind=kind_phys), pointer      :: raincd(:)          => null()  !<
    real (kind=kind_phys), pointer      :: raincs(:)          => null()  !<
    real (kind=kind_phys), pointer      :: rainmcadj(:)       => null()  !<
    real (kind=kind_phys), pointer      :: rainp(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: rb(:)              => null()  !<
    real (kind=kind_phys), pointer      :: rb_ice(:)          => null()  !<
    real (kind=kind_phys), pointer      :: rb_land(:)         => null()  !<
    real (kind=kind_phys), pointer      :: rb_ocean(:)        => null()  !<
    logical                             :: reset                         !<
    real (kind=kind_phys), pointer      :: rhc(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: runoff(:)          => null()  !<
    real (kind=kind_phys), pointer      :: save_q(:,:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: save_t(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: save_tcp(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: save_u(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: save_v(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: sbsno(:)           => null()  !<
    type (cmpfsw_type),    pointer      :: scmpsw(:)          => null()  !<
    real (kind=kind_phys), pointer      :: semis_ice(:)       => null()  !<
    real (kind=kind_phys), pointer      :: semis_land(:)      => null()  !<
    real (kind=kind_phys), pointer      :: semis_ocean(:)     => null()  !<
    real (kind=kind_phys), pointer      :: sfcalb(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: sigma(:)           => null()  !<
    real (kind=kind_phys), pointer      :: sigmaf(:)          => null()  !<
    real (kind=kind_phys), pointer      :: sigmafrac(:,:)     => null()  !<
    real (kind=kind_phys), pointer      :: sigmatot(:,:)      => null()  !<
    logical                             :: skip_macro                    !<
    integer, pointer                    :: slopetype(:)       => null()  !<
    real (kind=kind_phys), pointer      :: snowc(:)           => null()  !<
    real (kind=kind_phys), pointer      :: snohf(:)           => null()  !<
    real (kind=kind_phys), pointer      :: snowmp(:)          => null()  !<
    real (kind=kind_phys), pointer      :: snowmt(:)          => null()  !<
    integer, pointer                    :: soiltype(:)        => null()  !<
    real (kind=kind_phys), pointer      :: stress(:)          => null()  !<
    real (kind=kind_phys), pointer      :: stress_ice(:)      => null()  !<
    real (kind=kind_phys), pointer      :: stress_land(:)     => null()  !<
    real (kind=kind_phys), pointer      :: stress_ocean(:)    => null()  !<
    real (kind=kind_phys), pointer      :: t2mmp(:)           => null()  !<
    real (kind=kind_phys), pointer      :: theta(:)           => null()  !<
    real (kind=kind_phys), pointer      :: tice(:)            => null()  !<
    real (kind=kind_phys), pointer      :: tlvl(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: tlyr(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: tprcp_ice(:)       => null()  !<
    real (kind=kind_phys), pointer      :: tprcp_land(:)      => null()  !<
    real (kind=kind_phys), pointer      :: tprcp_ocean(:)     => null()  !<
    integer                             :: tracers_start_index           !<
    integer                             :: tracers_total                 !<
    integer                             :: tracers_water                 !<
    logical                             :: trans_aero                    !<
    real (kind=kind_phys), pointer      :: trans(:)           => null()  !<
    real (kind=kind_phys), pointer      :: tseal(:)           => null()  !<
    real (kind=kind_phys), pointer      :: tsfa(:)            => null()  !<
    real (kind=kind_phys), pointer      :: tsfg(:)            => null()  !<
    real (kind=kind_phys), pointer      :: tsnow(:)           => null()  !<
    real (kind=kind_phys), pointer      :: tsurf(:)           => null()  !<
    real (kind=kind_phys), pointer      :: ud_mf(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: ulwsfc_cice(:)     => null()  !<
    real (kind=kind_phys), pointer      :: dusfc_cice(:)      => null()  !<
    real (kind=kind_phys), pointer      :: dvsfc_cice(:)      => null()  !<
    real (kind=kind_phys), pointer      :: dqsfc_cice(:)      => null()  !<
    real (kind=kind_phys), pointer      :: dtsfc_cice(:)      => null()  !<
    real (kind=kind_phys), pointer      :: vdftra(:,:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: vegf1d(:)          => null()  !<
    integer, pointer                    :: vegtype(:)         => null()  !<
    real (kind=kind_phys), pointer      :: w_upi(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: wcbmax(:)          => null()  !<
    real (kind=kind_phys), pointer      :: weasd_ocean(:)     => null()  !<
    real (kind=kind_phys), pointer      :: weasd_land(:)      => null()  !<
    real (kind=kind_phys), pointer      :: weasd_ice(:)       => null()  !<
    real (kind=kind_phys), pointer      :: wind(:)            => null()  !<
    real (kind=kind_phys), pointer      :: work1(:)           => null()  !<
    real (kind=kind_phys), pointer      :: work2(:)           => null()  !<
    real (kind=kind_phys), pointer      :: work3(:)           => null()  !<
    real (kind=kind_phys), pointer      :: xcosz(:)           => null()  !<
    real (kind=kind_phys), pointer      :: xlai1d(:)          => null()  !<
    real (kind=kind_phys), pointer      :: xmu(:)             => null()  !<
    real (kind=kind_phys), pointer      :: z01d(:)            => null()  !<
    real (kind=kind_phys), pointer      :: zt1d(:)            => null()  !<
    real (kind=kind_phys), pointer      :: gw_dudt(:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: gw_dvdt(:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: gw_dtdt(:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: gw_kdis(:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: tau_tofd(:)        => null()  !< instantaneous momentum flux due to TOFD
    real (kind=kind_phys), pointer      :: tau_mtb(:)         => null()  !< instantaneous momentum flux due to mountain blocking drag
    real (kind=kind_phys), pointer      :: tau_ogw(:)         => null()  !< instantaneous momentum flux due to orographic gravity wave drag
    real (kind=kind_phys), pointer      :: tau_ngw(:)         => null()  !< instantaneous momentum flux due to nonstationary gravity waves
    real (kind=kind_phys), pointer      :: zmtb(:)            => null()  !< mountain blocking height
    real (kind=kind_phys), pointer      :: zlwb(:)            => null()  !< low level wave breaking height
    real (kind=kind_phys), pointer      :: zogw(:)            => null()  !< height of drag due to orographic gravity wave
    real (kind=kind_phys), pointer      :: dudt_mtb(:,:)      => null()  !< daily aver u-wind tend due to mountain blocking drag
    real (kind=kind_phys), pointer      :: dudt_ogw(:,:)      => null()  !< daily aver u-wind tend due to orographic gravity wave drag
    real (kind=kind_phys), pointer      :: dudt_tms(:,:)      => null()  !< daily aver u-wind tend due to TMS

    !-- HWRF physics: dry mixing ratios
    real (kind=kind_phys), pointer :: qv_r(:,:)               => null()  !<
    real (kind=kind_phys), pointer :: qc_r(:,:)               => null()  !<
    real (kind=kind_phys), pointer :: qi_r(:,:)               => null()  !<
    real (kind=kind_phys), pointer :: qr_r(:,:)               => null()  !<
    real (kind=kind_phys), pointer :: qs_r(:,:)               => null()  !<
    real (kind=kind_phys), pointer :: qg_r(:,:)               => null()  !<


    !-- Ferrier-Aligo MP scheme
    real (kind=kind_phys), pointer :: f_rain     (:,:)   => null()  !<
    real (kind=kind_phys), pointer :: f_ice      (:,:)   => null()  !<
    real (kind=kind_phys), pointer :: f_rimef    (:,:)   => null()  !<
    real (kind=kind_phys), pointer :: cwm        (:,:)   => null()  !<


    contains
      procedure :: create      => interstitial_create     !<   allocate array data
      procedure :: rad_reset   => interstitial_rad_reset  !<   reset array data for radiation
      procedure :: phys_reset  => interstitial_phys_reset !<   reset array data for physics
      procedure :: mprint      => interstitial_print      !<   print array data

  end type GFS_interstitial_type
#endif

!-------------------------
! GFS sub-containers
!-------------------------
#ifdef CCPP
!------------------------------------------------------------------------------------
! combined type of all of the above except GFS_control_type and GFS_interstitial_type
!------------------------------------------------------------------------------------
!! \section arg_table_GFS_data_type
!! \htmlinclude GFS_data_type.html
!!
  type GFS_data_type
     type(GFS_statein_type)  :: Statein
     type(GFS_stateout_type) :: Stateout
     type(GFS_sfcprop_type)  :: Sfcprop
     type(GFS_coupling_type) :: Coupling
     type(GFS_grid_type)     :: Grid
     type(GFS_tbd_type)      :: Tbd
     type(GFS_cldprop_type)  :: Cldprop
     type(GFS_radtend_type)  :: Radtend
     type(GFS_diag_type)     :: Intdiag
  end type GFS_data_type
#endif

!----------------
! PUBLIC ENTITIES
!----------------
  public GFS_init_type
  public GFS_statein_type,  GFS_stateout_type, GFS_sfcprop_type, &
         GFS_coupling_type
  public GFS_control_type,  GFS_grid_type,     GFS_tbd_type, &
         GFS_cldprop_type,  GFS_radtend_type,  GFS_diag_type
#ifdef CCPP
  public GFS_interstitial_type
#endif

!*******************************************************************************************
  CONTAINS

!------------------------
! GFS_statein_type%create
!------------------------
  subroutine statein_create (Statein, IM, Model) 
    implicit none

    class(GFS_statein_type)             :: Statein
    integer,                 intent(in) :: IM
    type(GFS_control_type),  intent(in) :: Model

    !--- level geopotential and pressures
    allocate (Statein%phii  (IM,Model%levs+1))
    allocate (Statein%prsi  (IM,Model%levs+1))
    allocate (Statein%prsik (IM,Model%levs+1))

    Statein%phii  = clear_val
    Statein%prsi  = clear_val
    Statein%prsik = clear_val

    !--- layer geopotential and pressures
    allocate (Statein%phil  (IM,Model%levs))
    allocate (Statein%prsl  (IM,Model%levs))
    allocate (Statein%prslk (IM,Model%levs))

    Statein%phil  = clear_val
    Statein%prsl  = clear_val
    Statein%prslk = clear_val

    !--- shared radiation and physics variables
    allocate (Statein%vvl  (IM,Model%levs))
    allocate (Statein%tgrs (IM,Model%levs))

    Statein%vvl  = clear_val
    Statein%tgrs = clear_val
! stochastic physics SKEB variable
    allocate (Statein%diss_est(IM,Model%levs))
    Statein%diss_est= clear_val
    !--- physics only variables
    allocate (Statein%pgr    (IM))
    allocate (Statein%ugrs   (IM,Model%levs))
    allocate (Statein%vgrs   (IM,Model%levs))
    allocate (Statein%qgrs   (IM,Model%levs,Model%ntrac))

    Statein%qgrs   = clear_val
    Statein%pgr    = clear_val
    Statein%ugrs   = clear_val
    Statein%vgrs   = clear_val

    !--- soil state variables - for soil SPPT - sfc-perts, mgehne
    allocate (Statein%smc  (IM,Model%lsoil))
    allocate (Statein%stc  (IM,Model%lsoil))
    allocate (Statein%slc  (IM,Model%lsoil))

    Statein%smc   = clear_val
    Statein%stc   = clear_val
    Statein%slc   = clear_val

  end subroutine statein_create


!-------------------------
! GFS_stateout_type%create
!-------------------------
  subroutine stateout_create (Stateout, IM, Model)

    implicit none

    class(GFS_stateout_type)           :: Stateout
    integer,                intent(in) :: IM
    type(GFS_control_type), intent(in) :: Model

    allocate (Stateout%gu0 (IM,Model%levs))
    allocate (Stateout%gv0 (IM,Model%levs))
    allocate (Stateout%gt0 (IM,Model%levs))
    allocate (Stateout%gq0 (IM,Model%levs,Model%ntrac))

    Stateout%gu0 = clear_val
    Stateout%gv0 = clear_val
    Stateout%gt0 = clear_val
    Stateout%gq0 = clear_val

 end subroutine stateout_create


!------------------------
! GFS_sfcprop_type%create
!------------------------
  subroutine sfcprop_create (Sfcprop, IM, Model)

    implicit none

    class(GFS_sfcprop_type)            :: Sfcprop
    integer,                intent(in) :: IM
    type(GFS_control_type), intent(in) :: Model

    !--- physics and radiation
    allocate (Sfcprop%slmsk    (IM))
    allocate (Sfcprop%oceanfrac(IM))
    allocate (Sfcprop%landfrac (IM))
    allocate (Sfcprop%lakefrac (IM))
    allocate (Sfcprop%tsfc     (IM))
    allocate (Sfcprop%tsfco    (IM))
    allocate (Sfcprop%tsfcl    (IM))
    allocate (Sfcprop%tisfc    (IM))
    allocate (Sfcprop%snowd    (IM))
    allocate (Sfcprop%zorl     (IM))
    allocate (Sfcprop%zorlo    (IM))
    allocate (Sfcprop%zorll    (IM))
    allocate (Sfcprop%fice     (IM))
!   allocate (Sfcprop%hprim    (IM))
    allocate (Sfcprop%hprime   (IM,Model%nmtvr))

    Sfcprop%slmsk     = clear_val
    Sfcprop%oceanfrac = clear_val
    Sfcprop%landfrac  = clear_val
    Sfcprop%lakefrac  = clear_val
    Sfcprop%tsfc      = clear_val
    Sfcprop%tsfco     = clear_val
    Sfcprop%tsfcl     = clear_val
    Sfcprop%tisfc     = clear_val
    Sfcprop%snowd     = clear_val
    Sfcprop%zorl      = clear_val
    Sfcprop%zorlo     = clear_val
    Sfcprop%zorll     = clear_val
    Sfcprop%fice      = clear_val
!   Sfcprop%hprim     = clear_val
    Sfcprop%hprime    = clear_val

!--- In (radiation only)
    allocate (Sfcprop%sncovr (IM))
    allocate (Sfcprop%snoalb (IM))
    allocate (Sfcprop%alvsf  (IM))
    allocate (Sfcprop%alnsf  (IM))
    allocate (Sfcprop%alvwf  (IM))
    allocate (Sfcprop%alnwf  (IM))
    allocate (Sfcprop%facsf  (IM))
    allocate (Sfcprop%facwf  (IM))

    Sfcprop%sncovr = clear_val
    Sfcprop%snoalb = clear_val
    Sfcprop%alvsf  = clear_val
    Sfcprop%alnsf  = clear_val
    Sfcprop%alvwf  = clear_val
    Sfcprop%alnwf  = clear_val
    Sfcprop%facsf  = clear_val
    Sfcprop%facwf  = clear_val

!--- physics surface props
!--- In
    allocate (Sfcprop%slope   (IM))
    allocate (Sfcprop%shdmin  (IM))
    allocate (Sfcprop%shdmax  (IM))
    allocate (Sfcprop%snoalb  (IM))
    allocate (Sfcprop%tg3     (IM))
    allocate (Sfcprop%vfrac   (IM))
    allocate (Sfcprop%vtype   (IM))
    allocate (Sfcprop%stype   (IM))
    allocate (Sfcprop%uustar  (IM))
    allocate (Sfcprop%oro     (IM))
    allocate (Sfcprop%oro_uf  (IM))

    Sfcprop%slope   = clear_val
    Sfcprop%shdmin  = clear_val
    Sfcprop%shdmax  = clear_val
    Sfcprop%snoalb  = clear_val
    Sfcprop%tg3     = clear_val
    Sfcprop%vfrac   = clear_val
    Sfcprop%vtype   = clear_val
    Sfcprop%stype   = clear_val
    Sfcprop%uustar  = clear_val
    Sfcprop%oro     = clear_val
    Sfcprop%oro_uf  = clear_val

    allocate (Sfcprop%qss_ice     (IM))
    allocate (Sfcprop%qss_land    (IM))
    allocate (Sfcprop%qss_ocean   (IM))
    allocate (Sfcprop%snowd_ice   (IM))
    allocate (Sfcprop%snowd_land  (IM))
    allocate (Sfcprop%snowd_ocean (IM))
    allocate (Sfcprop%tsfc_ice    (IM))
    allocate (Sfcprop%tsfc_land   (IM))
    allocate (Sfcprop%tsfc_ocean  (IM))
    allocate (Sfcprop%tsurf_ice   (IM))
    allocate (Sfcprop%tsurf_land  (IM))
    allocate (Sfcprop%tsurf_ocean (IM))
    allocate (Sfcprop%uustar_ice  (IM))
    allocate (Sfcprop%uustar_land (IM))
    allocate (Sfcprop%uustar_ocean(IM))
    allocate (Sfcprop%zorl_ice    (IM))
    allocate (Sfcprop%zorl_land   (IM))
    allocate (Sfcprop%zorl_ocean  (IM))
    allocate (Sfcprop%evap        (IM))
    allocate (Sfcprop%hflx        (IM))

    Sfcprop%qss_ice       = huge
    Sfcprop%qss_land      = huge
    Sfcprop%qss_ocean     = huge
    Sfcprop%snowd_ice     = huge
    Sfcprop%snowd_land    = huge
    Sfcprop%snowd_ocean   = huge
    Sfcprop%tsfc_ice      = huge
    Sfcprop%tsfc_land     = huge
    Sfcprop%tsfc_ocean    = huge
    Sfcprop%tsurf_ice     = huge
    Sfcprop%tsurf_land    = huge
    Sfcprop%tsurf_ocean   = huge
    Sfcprop%uustar_ice    = huge
    Sfcprop%uustar_land   = huge
    Sfcprop%uustar_ocean  = huge
    Sfcprop%zorl_ice      = huge
    Sfcprop%zorl_land     = huge
    Sfcprop%zorl_ocean    = huge
    Sfcprop%evap          = clear_val
    Sfcprop%hflx          = clear_val

!--- In/Out
    allocate (Sfcprop%hice   (IM))
    allocate (Sfcprop%weasd  (IM))
    allocate (Sfcprop%sncovr (IM))
    allocate (Sfcprop%canopy (IM))
    allocate (Sfcprop%ffmm   (IM))
    allocate (Sfcprop%ffhh   (IM))
    allocate (Sfcprop%f10m   (IM))
    allocate (Sfcprop%tprcp  (IM))
    allocate (Sfcprop%srflag (IM))
    allocate (Sfcprop%slc    (IM,Model%lsoil))
    allocate (Sfcprop%smc    (IM,Model%lsoil))
    allocate (Sfcprop%stc    (IM,Model%lsoil))

    Sfcprop%hice   = clear_val
    Sfcprop%weasd  = clear_val
    Sfcprop%sncovr = clear_val
    Sfcprop%canopy = clear_val
    Sfcprop%ffmm   = clear_val
    Sfcprop%ffhh   = clear_val
    Sfcprop%f10m   = clear_val
    Sfcprop%tprcp  = clear_val
    Sfcprop%srflag = clear_val
    Sfcprop%slc    = clear_val
    Sfcprop%smc    = clear_val
    Sfcprop%stc    = clear_val

!--- Out
    allocate (Sfcprop%t2m (IM))
#ifdef CCPP
    allocate (Sfcprop%th2m(IM))
#endif
    allocate (Sfcprop%q2m (IM))

    Sfcprop%t2m = clear_val
#ifdef CCPP
    Sfcprop%th2m = clear_val
#endif
    Sfcprop%q2m = clear_val

    if (Model%nstf_name(1) > 0) then
      allocate (Sfcprop%tref   (IM))
      allocate (Sfcprop%z_c    (IM))
      allocate (Sfcprop%c_0    (IM))
      allocate (Sfcprop%c_d    (IM))
      allocate (Sfcprop%w_0    (IM))
      allocate (Sfcprop%w_d    (IM))
      allocate (Sfcprop%xt     (IM))
      allocate (Sfcprop%xs     (IM))
      allocate (Sfcprop%xu     (IM))
      allocate (Sfcprop%xv     (IM))
      allocate (Sfcprop%xz     (IM))
      allocate (Sfcprop%zm     (IM))
      allocate (Sfcprop%xtts   (IM))
      allocate (Sfcprop%xzts   (IM))
      allocate (Sfcprop%d_conv (IM))
      allocate (Sfcprop%ifd    (IM))
      allocate (Sfcprop%dt_cool(IM))
      allocate (Sfcprop%qrain  (IM))

      Sfcprop%tref    = zero
      Sfcprop%z_c     = zero
      Sfcprop%c_0     = zero
      Sfcprop%c_d     = zero
      Sfcprop%w_0     = zero
      Sfcprop%w_d     = zero
      Sfcprop%xt      = zero
      Sfcprop%xs      = zero
      Sfcprop%xu      = zero
      Sfcprop%xv      = zero
      Sfcprop%xz      = zero
      Sfcprop%zm      = zero
      Sfcprop%xtts    = zero
      Sfcprop%xzts    = zero
      Sfcprop%d_conv  = zero
      Sfcprop%ifd     = zero
      Sfcprop%dt_cool = zero
      Sfcprop%qrain   = zero
    endif
    if (Model%lsm == Model%lsm_ruc .or. Model%lsm == Model%lsm_noahmp) then
      allocate(Sfcprop%raincprv  (IM))
      allocate(Sfcprop%rainncprv (IM))
      allocate(Sfcprop%iceprv    (IM))
      allocate(Sfcprop%snowprv   (IM))
      allocate(Sfcprop%graupelprv(IM))
      Sfcprop%raincprv   = clear_val
      Sfcprop%rainncprv  = clear_val
      Sfcprop%iceprv     = clear_val
      Sfcprop%snowprv    = clear_val
      Sfcprop%graupelprv = clear_val
    end if
! Noah MP allocate and init when used
!
    if (Model%lsm == Model%lsm_noahmp ) then

    allocate (Sfcprop%snowxy   (IM))
    allocate (Sfcprop%tvxy     (IM))
    allocate (Sfcprop%tgxy     (IM))
    allocate (Sfcprop%canicexy (IM))
    allocate (Sfcprop%canliqxy (IM))
    allocate (Sfcprop%eahxy    (IM))
    allocate (Sfcprop%tahxy    (IM))
    allocate (Sfcprop%cmxy     (IM))
    allocate (Sfcprop%chxy     (IM))
    allocate (Sfcprop%fwetxy   (IM))
    allocate (Sfcprop%sneqvoxy (IM))
    allocate (Sfcprop%alboldxy (IM))
    allocate (Sfcprop%qsnowxy  (IM))
    allocate (Sfcprop%wslakexy (IM))
    allocate (Sfcprop%zwtxy    (IM))
    allocate (Sfcprop%waxy     (IM))
    allocate (Sfcprop%wtxy     (IM))
    allocate (Sfcprop%lfmassxy (IM))
    allocate (Sfcprop%rtmassxy (IM))
    allocate (Sfcprop%stmassxy (IM))
    allocate (Sfcprop%woodxy   (IM))
    allocate (Sfcprop%stblcpxy (IM))
    allocate (Sfcprop%fastcpxy (IM))
    allocate (Sfcprop%xsaixy   (IM))
    allocate (Sfcprop%xlaixy   (IM))
    allocate (Sfcprop%taussxy  (IM))
    allocate (Sfcprop%smcwtdxy (IM))
    allocate (Sfcprop%deeprechxy (IM))
    allocate (Sfcprop%rechxy    (IM))
#ifdef CCPP
    allocate (Sfcprop%snicexy    (IM, Model%lsnow_lsm_lbound:0))
    allocate (Sfcprop%snliqxy    (IM, Model%lsnow_lsm_lbound:0))
    allocate (Sfcprop%tsnoxy     (IM, Model%lsnow_lsm_lbound:0))
    allocate (Sfcprop%smoiseq    (IM, Model%lsoil_lsm))
    allocate (Sfcprop%zsnsoxy    (IM, Model%lsnow_lsm_lbound:Model%lsoil_lsm))
#else
    allocate (Sfcprop%snicexy    (IM,-2:0))
    allocate (Sfcprop%snliqxy    (IM,-2:0))
    allocate (Sfcprop%tsnoxy     (IM,-2:0))
    allocate (Sfcprop%smoiseq    (IM, 1:4))
    allocate (Sfcprop%zsnsoxy    (IM,-2:4))
#endif

    Sfcprop%snowxy     = clear_val
    Sfcprop%tvxy       = clear_val
    Sfcprop%tgxy       = clear_val
    Sfcprop%canicexy   = clear_val
    Sfcprop%canliqxy   = clear_val
    Sfcprop%eahxy      = clear_val
    Sfcprop%tahxy      = clear_val
    Sfcprop%cmxy       = clear_val
    Sfcprop%chxy       = clear_val
    Sfcprop%fwetxy     = clear_val
    Sfcprop%sneqvoxy   = clear_val
    Sfcprop%alboldxy   = clear_val
    Sfcprop%qsnowxy    = clear_val
    Sfcprop%wslakexy   = clear_val
    Sfcprop%zwtxy      = clear_val
    Sfcprop%waxy       = clear_val
    Sfcprop%wtxy       = clear_val
    Sfcprop%lfmassxy   = clear_val
    Sfcprop%rtmassxy   = clear_val
    Sfcprop%stmassxy   = clear_val
    Sfcprop%woodxy     = clear_val
    Sfcprop%stblcpxy   = clear_val
    Sfcprop%fastcpxy   = clear_val
    Sfcprop%xsaixy     = clear_val
    Sfcprop%xlaixy     = clear_val
    Sfcprop%taussxy    = clear_val
    Sfcprop%smcwtdxy   = clear_val
    Sfcprop%deeprechxy = clear_val
    Sfcprop%rechxy     = clear_val

    Sfcprop%snicexy    = clear_val
    Sfcprop%snliqxy    = clear_val
    Sfcprop%tsnoxy     = clear_val
    Sfcprop%smoiseq    = clear_val
    Sfcprop%zsnsoxy    = clear_val
    
    allocate(Sfcprop%draincprv  (IM))
    allocate(Sfcprop%drainncprv (IM))
    allocate(Sfcprop%diceprv    (IM))
    allocate(Sfcprop%dsnowprv   (IM))
    allocate(Sfcprop%dgraupelprv(IM))

    Sfcprop%draincprv   = clear_val
    Sfcprop%drainncprv  = clear_val
    Sfcprop%diceprv     = clear_val
    Sfcprop%dsnowprv    = clear_val
    Sfcprop%dgraupelprv = clear_val
    
   endif

#ifdef CCPP
    if (Model%lsm == Model%lsm_ruc) then
       ! For land surface models with different numbers of levels than the four NOAH levels
       allocate (Sfcprop%wetness     (IM))
       allocate (Sfcprop%sh2o        (IM,Model%lsoil_lsm))
       allocate (Sfcprop%keepsmfr    (IM,Model%lsoil_lsm))
       allocate (Sfcprop%smois       (IM,Model%lsoil_lsm))
       allocate (Sfcprop%tslb        (IM,Model%lsoil_lsm))
       allocate (Sfcprop%flag_frsoil (IM,Model%lsoil_lsm))
       allocate (Sfcprop%zs          (Model%lsoil_lsm))
       allocate (Sfcprop%clw_surf    (IM))
       allocate (Sfcprop%qwv_surf    (IM))
       allocate (Sfcprop%cndm_surf   (IM))
       allocate (Sfcprop%rhofr       (IM))
       allocate (Sfcprop%tsnow       (IM))
       allocate (Sfcprop%snowfallac  (IM))
       allocate (Sfcprop%acsnow      (IM))
       !
       Sfcprop%wetness     = clear_val
       Sfcprop%sh2o        = clear_val
       Sfcprop%keepsmfr    = clear_val
       Sfcprop%smois       = clear_val
       Sfcprop%tslb        = clear_val
       Sfcprop%zs          = clear_val
       Sfcprop%clw_surf    = clear_val
       Sfcprop%qwv_surf    = clear_val
       Sfcprop%cndm_surf   = clear_val
       Sfcprop%flag_frsoil = clear_val
       Sfcprop%rhofr       = clear_val
       Sfcprop%tsnow       = clear_val
       Sfcprop%snowfallac  = clear_val
       Sfcprop%acsnow      = clear_val
       !
       if (Model%rdlai) then
          allocate (Sfcprop%xlaixy (IM))
          Sfcprop%xlaixy = clear_val
       end if
              
    end if
    if (Model%do_mynnsfclay) then
    ! For MYNN surface layer scheme
       !print*,"Allocating all MYNN-sfclay variables"
       allocate (Sfcprop%ustm   (IM ))
       allocate (Sfcprop%zol    (IM ))
       allocate (Sfcprop%mol    (IM ))
       allocate (Sfcprop%rmol   (IM ))
       allocate (Sfcprop%flhc   (IM ))
       allocate (Sfcprop%flqc   (IM ))
       allocate (Sfcprop%chs2   (IM ))
       allocate (Sfcprop%cqs2   (IM ))
       allocate (Sfcprop%lh     (IM ))
       !
       !print*,"Initializing all MYNN-SfcLay variables with ",clear_val
       Sfcprop%ustm        = clear_val
       Sfcprop%zol         = clear_val
       Sfcprop%mol         = clear_val
       Sfcprop%rmol        = clear_val
       Sfcprop%flhc        = clear_val
       Sfcprop%flqc        = clear_val
       Sfcprop%chs2        = clear_val
       Sfcprop%cqs2        = clear_val
       Sfcprop%lh          = clear_val
    end if
    if (Model%imfdeepcnv == Model%imfdeepcnv_gf) then
        allocate (Sfcprop%conv_act(IM))
        Sfcprop%conv_act = zero
    end if
    
#endif

  end subroutine sfcprop_create


!-------------------------
! GFS_coupling_type%create
!-------------------------
  subroutine coupling_create (Coupling, IM, Model)

    implicit none

    class(GFS_coupling_type)           :: Coupling
    integer,                intent(in) :: IM
    type(GFS_control_type), intent(in) :: Model

    !--- radiation out
    !--- physics in
    allocate (Coupling%nirbmdi (IM))
    allocate (Coupling%nirdfdi (IM))
    allocate (Coupling%visbmdi (IM))   
    allocate (Coupling%visdfdi (IM))   
    allocate (Coupling%nirbmui (IM))   
    allocate (Coupling%nirdfui (IM))   
    allocate (Coupling%visbmui (IM))   
    allocate (Coupling%visdfui (IM))   

    Coupling%nirbmdi = clear_val
    Coupling%nirdfdi = clear_val
    Coupling%visbmdi = clear_val
    Coupling%visdfdi = clear_val
    Coupling%nirbmui = clear_val
    Coupling%nirdfui = clear_val
    Coupling%visbmui = clear_val
    Coupling%visdfui = clear_val

    allocate (Coupling%sfcdsw (IM))
    allocate (Coupling%sfcnsw (IM))
    allocate (Coupling%sfcdlw (IM))

    Coupling%sfcdsw = clear_val
    Coupling%sfcnsw = clear_val
    Coupling%sfcdlw = clear_val

    if (Model%cplflx .or. Model%do_sppt .or. Model%cplchm) then
      allocate (Coupling%rain_cpl (IM))
      allocate (Coupling%snow_cpl (IM))
      Coupling%rain_cpl = clear_val
      Coupling%snow_cpl = clear_val
    endif

    if (Model%cplflx .or. Model%cplwav) then
      !--- instantaneous quantities 
      allocate (Coupling%u10mi_cpl (IM))
      allocate (Coupling%v10mi_cpl (IM))

      Coupling%u10mi_cpl = clear_val
      Coupling%v10mi_cpl = clear_val
    endif 

    if (Model%cplflx) then
      !--- incoming quantities
      allocate (Coupling%slimskin_cpl (IM))
      allocate (Coupling%dusfcin_cpl  (IM))
      allocate (Coupling%dvsfcin_cpl  (IM))
      allocate (Coupling%dtsfcin_cpl  (IM))
      allocate (Coupling%dqsfcin_cpl  (IM))
      allocate (Coupling%ulwsfcin_cpl (IM))
      allocate (Coupling%tseain_cpl   (IM))
      allocate (Coupling%tisfcin_cpl  (IM))
      allocate (Coupling%ficein_cpl   (IM))
      allocate (Coupling%hicein_cpl   (IM))
      allocate (Coupling%hsnoin_cpl   (IM))

      Coupling%slimskin_cpl = clear_val
      Coupling%dusfcin_cpl  = clear_val
      Coupling%dvsfcin_cpl  = clear_val
      Coupling%dtsfcin_cpl  = clear_val
      Coupling%dqsfcin_cpl  = clear_val
      Coupling%ulwsfcin_cpl = clear_val
      Coupling%tseain_cpl   = clear_val
      Coupling%tisfcin_cpl  = clear_val
      Coupling%ficein_cpl   = clear_val
      Coupling%hicein_cpl   = clear_val
      Coupling%hsnoin_cpl   = clear_val

      !--- accumulated quantities
      allocate (Coupling%dusfc_cpl  (IM))
      allocate (Coupling%dvsfc_cpl  (IM))
      allocate (Coupling%dtsfc_cpl  (IM))
      allocate (Coupling%dqsfc_cpl  (IM))
      allocate (Coupling%dlwsfc_cpl (IM))
      allocate (Coupling%dswsfc_cpl (IM))
      allocate (Coupling%dnirbm_cpl (IM))
      allocate (Coupling%dnirdf_cpl (IM))
      allocate (Coupling%dvisbm_cpl (IM))
      allocate (Coupling%dvisdf_cpl (IM))
      allocate (Coupling%nlwsfc_cpl (IM))
      allocate (Coupling%nswsfc_cpl (IM))
      allocate (Coupling%nnirbm_cpl (IM))
      allocate (Coupling%nnirdf_cpl (IM))
      allocate (Coupling%nvisbm_cpl (IM))
      allocate (Coupling%nvisdf_cpl (IM))

      Coupling%dusfc_cpl  = clear_val
      Coupling%dvsfc_cpl  = clear_val
      Coupling%dtsfc_cpl  = clear_val
      Coupling%dqsfc_cpl  = clear_val
      Coupling%dlwsfc_cpl = clear_val
      Coupling%dswsfc_cpl = clear_val
      Coupling%dnirbm_cpl = clear_val
      Coupling%dnirdf_cpl = clear_val
      Coupling%dvisbm_cpl = clear_val
      Coupling%dvisdf_cpl = clear_val
      Coupling%nlwsfc_cpl = clear_val
      Coupling%nswsfc_cpl = clear_val
      Coupling%nnirbm_cpl = clear_val
      Coupling%nnirdf_cpl = clear_val
      Coupling%nvisbm_cpl = clear_val
      Coupling%nvisdf_cpl = clear_val

      !--- instantaneous quantities
      allocate (Coupling%dusfci_cpl  (IM))
      allocate (Coupling%dvsfci_cpl  (IM))
      allocate (Coupling%dtsfci_cpl  (IM))
      allocate (Coupling%dqsfci_cpl  (IM))
      allocate (Coupling%dlwsfci_cpl (IM))
      allocate (Coupling%dswsfci_cpl (IM))
      allocate (Coupling%dnirbmi_cpl (IM))
      allocate (Coupling%dnirdfi_cpl (IM))
      allocate (Coupling%dvisbmi_cpl (IM))
      allocate (Coupling%dvisdfi_cpl (IM))
      allocate (Coupling%nlwsfci_cpl (IM))
      allocate (Coupling%nswsfci_cpl (IM))
      allocate (Coupling%nnirbmi_cpl (IM))
      allocate (Coupling%nnirdfi_cpl (IM))
      allocate (Coupling%nvisbmi_cpl (IM))
      allocate (Coupling%nvisdfi_cpl (IM))
      allocate (Coupling%t2mi_cpl    (IM))
      allocate (Coupling%q2mi_cpl    (IM))
      allocate (Coupling%tsfci_cpl   (IM))
      allocate (Coupling%psurfi_cpl  (IM))
      allocate (Coupling%oro_cpl     (IM))
      allocate (Coupling%slmsk_cpl   (IM))

      Coupling%dusfci_cpl  = clear_val
      Coupling%dvsfci_cpl  = clear_val
      Coupling%dtsfci_cpl  = clear_val
      Coupling%dqsfci_cpl  = clear_val
      Coupling%dlwsfci_cpl = clear_val
      Coupling%dswsfci_cpl = clear_val
      Coupling%dnirbmi_cpl = clear_val
      Coupling%dnirdfi_cpl = clear_val
      Coupling%dvisbmi_cpl = clear_val
      Coupling%dvisdfi_cpl = clear_val
      Coupling%nlwsfci_cpl = clear_val
      Coupling%nswsfci_cpl = clear_val
      Coupling%nnirbmi_cpl = clear_val
      Coupling%nnirdfi_cpl = clear_val
      Coupling%nvisbmi_cpl = clear_val
      Coupling%nvisdfi_cpl = clear_val
      Coupling%t2mi_cpl    = clear_val
      Coupling%q2mi_cpl    = clear_val
      Coupling%tsfci_cpl   = clear_val
      Coupling%psurfi_cpl  = clear_val
      Coupling%oro_cpl     = clear_val  !< pointer to sfcprop%oro
      Coupling%slmsk_cpl   = clear_val  !< pointer to sfcprop%slmsk
    endif

   !-- cellular automata
    if (Model%do_ca) then
      allocate (Coupling%tconvtend (IM,Model%levs))
      allocate (Coupling%qconvtend (IM,Model%levs))
      allocate (Coupling%uconvtend (IM,Model%levs))
      allocate (Coupling%vconvtend (IM,Model%levs))
      allocate (Coupling%cape     (IM))
      allocate (Coupling%ca_out   (IM))
      allocate (Coupling%ca_deep  (IM))
      allocate (Coupling%ca_turb  (IM))
      allocate (Coupling%ca_shal  (IM))
      allocate (Coupling%ca_rad   (IM))
      allocate (Coupling%ca_micro (IM))
      Coupling%ca_out    = clear_val
      Coupling%ca_deep   = clear_val
      Coupling%ca_turb   = clear_val
      Coupling%ca_shal   = clear_val
      Coupling%ca_rad    = clear_val
      Coupling%ca_micro  = clear_val   
      Coupling%cape      = clear_val
      Coupling%tconvtend = clear_val
      Coupling%qconvtend = clear_val
      Coupling%uconvtend = clear_val
      Coupling%vconvtend = clear_val
    endif

    ! -- GSDCHEM coupling options
    if (Model%cplchm) then
      !--- outgoing instantaneous quantities
      allocate (Coupling%ushfsfci  (IM))
      allocate (Coupling%dkt       (IM,Model%levs))
      allocate (Coupling%dqdti     (IM,Model%levs))
      !--- accumulated convective rainfall
      allocate (Coupling%rainc_cpl (IM))

      Coupling%rainc_cpl = clear_val
      Coupling%ushfsfci  = clear_val
      Coupling%dkt       = clear_val
      Coupling%dqdti     = clear_val
    endif

    !--- stochastic physics option
    if (Model%do_sppt) then
      allocate (Coupling%sppt_wts  (IM,Model%levs))
      Coupling%sppt_wts = clear_val
    endif

    !--- stochastic shum option
    if (Model%do_shum) then
      allocate (Coupling%shum_wts  (IM,Model%levs))
      Coupling%shum_wts = clear_val
    endif

    !--- stochastic skeb option
    if (Model%do_skeb) then
      allocate (Coupling%skebu_wts (IM,Model%levs))
      allocate (Coupling%skebv_wts (IM,Model%levs))

      Coupling%skebu_wts = clear_val
      Coupling%skebv_wts = clear_val
    endif

    !--- stochastic physics option
    if (Model%do_sfcperts) then
      allocate (Coupling%sfc_wts  (IM,Model%nsfcpert))
      Coupling%sfc_wts = clear_val
    endif

    !--- needed for Thompson's aerosol option
    if(Model%imp_physics == Model%imp_physics_thompson .and. Model%ltaerosol) then 
      allocate (Coupling%nwfa2d (IM))
      allocate (Coupling%nifa2d (IM))
      Coupling%nwfa2d   = clear_val
      Coupling%nifa2d   = clear_val
    endif

#ifdef CCPP
    if (Model%imfdeepcnv == Model%imfdeepcnv_gf) then
      allocate (Coupling%qci_conv (IM,Model%levs))
      Coupling%qci_conv   = clear_val
    endif
#endif

  end subroutine coupling_create


!----------------------
! GFS_control_type%init
!----------------------
  subroutine control_initialize (Model, nlunit, fn_nml, me, master, &
                                 logunit, isc, jsc, nx, ny, levs,   &
                                 cnx, cny, gnx, gny, dt_dycore,     &
                                 dt_phys, iau_offset, idat, jdat,   &
                                 tracer_names,                      &
                                 input_nml_file, tile_num, blksz    &
#ifdef CCPP
                                ,ak, bk, restart, hydrostatic,      &
                                 communicator, ntasks, nthreads     &
#endif
                                 )

!--- modules
#ifdef CCPP
    use physcons,         only: con_rerth, con_pi
#else
    use physcons,         only: dxmax, dxmin, dxinv, con_rerth, con_pi, rhc_max
#endif
    use mersenne_twister, only: random_setseed, random_number
#ifndef CCPP
    use module_ras,       only: nrcmax
#endif
    use parse_tracers,    only: get_tracer_index
#ifndef CCPP
    use wam_f107_kp_mod,  only: f107_kp_size, f107_kp_interval,     &
                                f107_kp_skip_size, f107_kp_data_size
#endif
    implicit none

!--- interface variables
    class(GFS_control_type)            :: Model
    integer,                intent(in) :: nlunit
    character(len=64),      intent(in) :: fn_nml
    integer,                intent(in) :: me
    integer,                intent(in) :: master
    integer,                intent(in) :: logunit
    integer,                intent(in) :: tile_num
    integer,                intent(in) :: isc
    integer,                intent(in) :: jsc
    integer,                intent(in) :: nx
    integer,                intent(in) :: ny
    integer,                intent(in) :: levs
    integer,                intent(in) :: cnx
    integer,                intent(in) :: cny
    integer,                intent(in) :: gnx
    integer,                intent(in) :: gny
    real(kind=kind_phys),   intent(in) :: dt_dycore
    real(kind=kind_phys),   intent(in) :: dt_phys
    integer,                intent(in) :: iau_offset
    integer,                intent(in) :: idat(8)
    integer,                intent(in) :: jdat(8)
    character(len=32),      intent(in) :: tracer_names(:)
    character(len=256),     intent(in), pointer :: input_nml_file(:)
    integer,                intent(in) :: blksz(:)
#ifdef CCPP
    real(kind=kind_phys), dimension(:), intent(in) :: ak
    real(kind=kind_phys), dimension(:), intent(in) :: bk
    logical,                intent(in) :: restart
    logical,                intent(in) :: hydrostatic
    integer,                intent(in) :: communicator
    integer,                intent(in) :: ntasks
    integer,                intent(in) :: nthreads
#endif
    !--- local variables
    integer :: i, j, n
    integer :: ios
    integer :: seed0
    logical :: exists
    real(kind=kind_phys) :: tem
    real(kind=kind_phys) :: rinc(5)
    real(kind=kind_phys) :: wrk(1)
    real(kind=kind_phys), parameter :: con_hr = 3600.

!--- BEGIN NAMELIST VARIABLES
    real(kind=kind_phys) :: fhzero         = 0.0             !< hours between clearing of diagnostic buckets
    logical              :: ldiag3d        = .false.         !< flag for 3d diagnostic fields
    logical              :: lssav          = .false.         !< logical flag for storing diagnostics

    real(kind=kind_phys) :: fhcyc          = 0.              !< frequency for surface data cycling (hours)
    integer              :: thermodyn_id   =  1              !< valid for GFS only for get_prs/phi
    integer              :: sfcpress_id    =  1              !< valid for GFS only for get_prs/phi

    !--- coupling parameters
    logical              :: cplflx         = .false.         !< default no cplflx collection
    logical              :: cplwav         = .false.         !< default no cplwav collection
    logical              :: cplchm         = .false.         !< default no cplchm collection

!--- integrated dynamics through earth's atmosphere
    logical              :: lsidea         = .false.

!--- radiation parameters
    real(kind=kind_phys) :: fhswr          = 3600.           !< frequency for shortwave radiation (secs)
    real(kind=kind_phys) :: fhlwr          = 3600.           !< frequency for longwave radiation (secs)
#ifdef CCPP
    integer              :: nhfrad         = 0               !< number of timesteps for which to call radiation on physics timestep (coldstarts)
#endif
    integer              :: levr           = -99             !< number of vertical levels for radiation calculations
    integer              :: nfxr           = 39+6            !< second dimension of input/output array fluxr   
    logical              :: aero_in        = .false.         !< flag for initializing aero data 
    logical              :: iccn           = .false.         !< logical to use IN CCN forcing for MG2/3
    integer              :: iflip          =  1              !< iflip - is not the same as flipv
    integer              :: isol           =  0              !< use prescribed solar constant
    integer              :: ico2           =  0              !< prescribed global mean value (old opernl)
    integer              :: ialb           =  0              !< use climatology alb, based on sfc type
                                                             !< 1 => use modis based alb
    integer              :: iems           =  0              !< use fixed value of 1.0
    integer              :: iaer           =  1              !< default aerosol effect in sw only
    integer              :: icliq_sw       =  1              !< sw optical property for liquid clouds
    integer              :: iovr_sw        =  1              !< sw: max-random overlap clouds
    integer              :: iovr_lw        =  1              !< lw: max-random overlap clouds
    integer              :: ictm           =  1              !< ictm=0 => use data at initial cond time, if not
                                                             !<           available; use latest; no extrapolation.
                                                             !< ictm=1 => use data at the forecast time, if not
                                                             !<           available; use latest; do extrapolation.
                                                             !< ictm=yyyy0 => use yyyy data for the forecast time;
                                                             !<           no extrapolation.
                                                             !< ictm=yyyy1 = > use yyyy data for the fcst. If needed, 
                                                             !<           do extrapolation to match the fcst time.
                                                             !< ictm=-1 => use user provided external data for
                                                             !<           the fcst time; no extrapolation.
                                                             !< ictm=-2 => same as ictm=0, but add seasonal cycle
                                                             !<           from climatology; no extrapolation.
    integer              :: isubc_sw          =  0           !< sw clouds without sub-grid approximation
    integer              :: isubc_lw          =  0           !< lw clouds without sub-grid approximation
                                                             !< =1 => sub-grid cloud with prescribed seeds
                                                             !< =2 => sub-grid cloud with randomly generated
                                                             !< seeds
    logical              :: crick_proof       = .false.      !< CRICK-Proof cloud water
    logical              :: ccnorm            = .false.      !< Cloud condensate normalized by cloud cover 
    logical              :: norad_precip      = .false.      !< radiation precip flag for Ferrier/Moorthi
    logical              :: lwhtr             = .true.       !< flag to output lw heating rate (Radtend%lwhc)
    logical              :: swhtr             = .true.       !< flag to output sw heating rate (Radtend%swhc)

!--- Z-C microphysical parameters
    integer              :: ncld              =  1                 !< choice of cloud scheme
    integer              :: imp_physics       =  99                !< choice of cloud scheme
    real(kind=kind_phys) :: psautco(2)        = (/6.0d-4,3.0d-4/)  !< [in] auto conversion coeff from ice to snow
    real(kind=kind_phys) :: prautco(2)        = (/1.0d-4,1.0d-4/)  !< [in] auto conversion coeff from cloud to rain
    real(kind=kind_phys) :: evpco             = 2.0d-5             !< [in] coeff for evaporation of largescale rain
    real(kind=kind_phys) :: wminco(2)         = (/1.0d-5,1.0d-5/)  !< [in] water and ice minimum threshold for Zhao
!---Max hourly
    real(kind=kind_phys) :: avg_max_length = 3600.              !< reset value in seconds for max hourly.
!--- Ferrier-Aligo microphysical parameters
#ifdef CCPP
    real(kind=kind_phys) :: rhgrd             = 0.98               !< fer_hires microphysics only     
    logical              :: spec_adv          = .true.            !< Individual cloud species advected
#endif
!--- M-G microphysical parameters
    integer              :: fprcp             =  0                 !< no prognostic rain and snow (MG)
    integer              :: pdfflag           =  4                 !< pdf flag for MG macro physics
    real(kind=kind_phys) :: mg_dcs            = 200.0              !< Morrison-Gettelman microphysics parameters
    real(kind=kind_phys) :: mg_qcvar          = 1.0
    real(kind=kind_phys) :: mg_ts_auto_ice(2) = (/180.0,180.0/)    !< ice auto conversion time scale
#ifdef CCPP
    real(kind=kind_phys) :: mg_rhmini         = 1.01               !< relative humidity threshold parameter for nucleating ice
#endif
    real(kind=kind_phys) :: mg_ncnst          = 100.e6             !< constant droplet num concentration (m-3)
    real(kind=kind_phys) :: mg_ninst          = 0.15e6             !< constant ice num concentration (m-3)
    real(kind=kind_phys) :: mg_ngnst          = 0.10e6             !< constant graupel/hail num concentration (m-3) = 0.1e6_r8
    real(kind=kind_phys) :: mg_alf            = 1.0                !< tuning factor for alphs in MG macrophysics
    real(kind=kind_phys) :: mg_qcmin(2)       = (/1.0d-9,1.0d-9/)  !< min liquid and ice mixing ratio in Mg macro clouds
    real(kind=kind_phys) :: mg_berg_eff_factor = 2.0               !< berg efficiency factor
    character(len=16)    :: mg_precip_frac_method = 'max_overlap'  !< type of precipitation fraction method
#ifdef CCPP
    real(kind=kind_phys) :: tf              = 258.16d0
    real(kind=kind_phys) :: tcr             = 273.16d0
#endif
!
    logical              :: effr_in         = .false.              !< flag to use effective radii of cloud species in radiation
    logical              :: microp_uniform  = .true.
    logical              :: do_cldliq       = .true.
    logical              :: do_cldice       = .true.
    logical              :: hetfrz_classnuc = .false.
    logical              :: mg_nccons       = .false.           !< set .true. to specify constant cloud droplet number
    logical              :: mg_nicons       = .false.           !< set .true. to specify constant cloud ice number
    logical              :: mg_ngcons       = .false.           !< set .true. to specify constant graupel/hail number
    logical              :: sed_supersat    = .true.
    logical              :: do_sb_physics   = .true.
    logical              :: mg_do_graupel   = .true.            !< set .true. to turn on prognostic grapuel (with fprcp=2)
    logical              :: mg_do_hail      = .false.           !< set .true. to turn on prognostic hail (with fprcp=2)
    logical              :: mg_do_ice_gmao  = .false.           !< set .true. to turn on gmao ice formulation
    logical              :: mg_do_liq_liu   = .true.            !< set .true. to turn on liu liquid treatment


    !--- Thompson microphysical parameters
    logical              :: ltaerosol      = .false.            !< flag for aerosol version
    logical              :: lradar         = .false.            !< flag for radar reflectivity 
    real(kind=kind_phys) :: ttendlim       = -999.0             !< temperature tendency limiter, set to <0 to deactivate

    !--- GFDL microphysical parameters
    logical              :: lgfdlmprad     = .false.            !< flag for GFDLMP radiation interaction 

    !--- Thompson,GFDL microphysical parameter
    logical              :: lrefres        = .false.            !< flag for radar reflectivity in restart file

    !--- land/surface model parameters
    integer              :: lsm            =  1              !< flag for land surface model to use =0  for osu lsm; =1  for noah lsm; =2  for noah mp lsm; =3  for RUC lsm
    integer              :: lsoil          =  4              !< number of soil layers
#ifdef CCPP
    integer              :: lsoil_lsm      =  -1             !< number of soil layers internal to land surface model; -1 use lsoil
    integer              :: lsnow_lsm      =  3              !< maximum number of snow layers internal to land surface model
    logical              :: rdlai          = .false.
#endif
    integer              :: ivegsrc        =  2              !< ivegsrc = 0   => USGS,
                                                             !< ivegsrc = 1   => IGBP (20 category)
                                                             !< ivegsrc = 2   => UMD  (13 category)
    integer              :: isot           =  0              !< isot = 0   => Zobler soil type  ( 9 category)
                                                             !< isot = 1   => STATSGO soil type (19 category)
    ! -- to use Noah MP, lsm needs to be set to 2 and both ivegsrc and isot are set
    ! to 1 - MODIS IGBP and STATSGO - the defaults are the same as in the
    ! scripts;change from namelist

    integer              :: iopt_dveg      =  4  ! 4 -> off (use table lai; use maximum vegetation fraction)
    integer              :: iopt_crs       =  1  !canopy stomatal resistance (1-> ball-berry; 2->jarvis)
    integer              :: iopt_btr       =  1  !soil moisture factor for stomatal resistance (1-> noah; 2-> clm; 3-> ssib)
    integer              :: iopt_run       =  3  !runoff and groundwater (1->simgm; 2->simtop; 3->schaake96; 4->bats)
    integer              :: iopt_sfc       =  1  !surface layer drag coeff (ch & cm) (1->m-o; 2->chen97)
    integer              :: iopt_frz       =  1  !supercooled liquid water (1-> ny06; 2->koren99)
    integer              :: iopt_inf       =  1  !frozen soil permeability (1-> ny06; 2->koren99)
    integer              :: iopt_rad       =  3  !radiation transfer (1->gap=f(3d,cosz); 2->gap=0; 3->gap=1-fveg)
    integer              :: iopt_alb       =  2  !snow surface albedo (1->bats; 2->class)
    integer              :: iopt_snf       =  1  !rainfall & snowfall (1-jordan91; 2->bats; 3->noah)
    integer              :: iopt_tbot      =  2  !lower boundary of soil temperature (1->zero-flux; 2->noah)
    integer              :: iopt_stc       =  1  !snow/soil temperature time scheme (only layer 1)

    logical              :: use_ufo        = .false.         !< flag for gcycle surface option

!--- tuning parameters for physical parameterizations
    logical              :: ras            = .false.                  !< flag for ras convection scheme
    logical              :: flipv          = .true.                   !< flag for vertical direction flip (ras)
                                                                      !< .true. implies surface at k=1
    logical              :: trans_trac     = .false.                  !< flag for convective transport of tracers (RAS, CS, or SAMF)
    logical              :: old_monin      = .false.                  !< flag for diff monin schemes
    logical              :: cnvgwd         = .false.                  !< flag for conv gravity wave drag
    integer              :: gwd_opt        =  1                       !< flag for configuring gwd scheme
                                                                      !< gwd_opt = 3 : GSDdrag suite
                                                                      !< gwd_opt = 33: GSDdrag suite with extra output
!--- vay-2018
    logical              :: ldiag_ugwp     = .false.                  !< flag for UGWP diag fields
    logical              :: do_ugwp        = .false.                  !< flag do UGWP+RF
    logical              :: do_tofd        = .false.                  !< flag do Turb oro Form Drag

    logical              :: do_gwd         = .false.                  !< flag for running gravity wave drag
    logical              :: do_cnvgwd      = .false.                  !< flag for running conv gravity wave drag

    logical              :: mstrat         = .false.                  !< flag for moorthi approach for stratus
    logical              :: moist_adj      = .false.                  !< flag for moist convective adjustment
    logical              :: cscnv          = .false.                  !< flag for Chikira-Sugiyama convection
    logical              :: cal_pre        = .false.                  !< flag controls precip type algorithm
    logical              :: do_aw          = .false.                  !< AW scale-aware option in cs convection
    logical              :: do_awdd        = .false.                  !< AW scale-aware option in cs convection
    logical              :: flx_form       = .false.                  !< AW scale-aware option in cs convection
    logical              :: do_shoc        = .false.                  !< flag for SHOC
    logical              :: shocaftcnv     = .false.                  !< flag for SHOC
    logical              :: shoc_cld       = .false.                  !< flag for SHOC in grrad
#ifdef CCPP
    logical              :: oz_phys        = .true.                   !< flag for old (2006) ozone physics
    logical              :: oz_phys_2015   = .false.                  !< flag for new (2015) ozone physics
#endif
    logical              :: h2o_phys       = .false.                  !< flag for stratosphere h2o
    logical              :: pdfcld         = .false.                  !< flag for pdfcld
    logical              :: shcnvcw        = .false.                  !< flag for shallow convective cloud
    logical              :: redrag         = .false.                  !< flag for reduced drag coeff. over sea
    logical              :: hybedmf        = .false.                  !< flag for hybrid edmf pbl scheme
    logical              :: satmedmf       = .false.                  !< flag for scale-aware TKE-based moist edmf
                                                                      !< vertical turbulent mixing scheme
    logical              :: shinhong       = .false.                  !< flag for scale-aware Shinhong vertical turbulent mixing scheme
    logical              :: do_ysu         = .false.                  !< flag for YSU vertical turbulent mixing scheme
    logical              :: dspheat        = .false.                  !< flag for tke dissipative heating
    logical              :: lheatstrg      = .false.                  !< flag for canopy heat storage parameterization
    logical              :: cnvcld         = .false.
    logical              :: random_clds    = .false.                  !< flag controls whether clouds are random
    logical              :: shal_cnv       = .false.                  !< flag for calling shallow convection
    integer              :: imfshalcnv     =  1                       !< flag for mass-flux shallow convection scheme
                                                                      !<     1: July 2010 version of mass-flux shallow conv scheme
                                                                      !<         current operational version as of 2016
                                                                      !<     2: scale- & aerosol-aware mass-flux shallow conv scheme (2017)
                                                                      !<     3: scale- & aerosol-aware Grell-Freitas scheme (GSD)
                                                                      !<     4: New Tiedtke scheme (CAPS)
                                                                      !<     0: modified Tiedtke's eddy-diffusion shallow conv scheme
                                                                      !<    -1: no shallow convection used
    integer              :: imfdeepcnv     =  1                       !< flag for mass-flux deep convection scheme
                                                                      !<     1: July 2010 version of SAS conv scheme
                                                                      !<           current operational version as of 2016
                                                                      !<     2: scale- & aerosol-aware mass-flux deep conv scheme (2017)
                                                                      !<     3: scale- & aerosol-aware Grell-Freitas scheme (GSD)
                                                                      !<     4: New Tiedtke scheme (CAPS)
    integer              :: isatmedmf      =  0                       !< flag for scale-aware TKE-based moist edmf scheme
                                                                      !<     0: initial version of satmedmf (Nov. 2018)
                                                                      !<     1: updated version of satmedmf (as of May 2019)
    logical              :: do_deep        = .true.                   !< whether to do deep convection
#ifdef CCPP
    logical              :: do_mynnedmf       = .false.               !< flag for MYNN-EDMF
    logical              :: do_mynnsfclay     = .false.               !< flag for MYNN Surface Layer Scheme
    ! DH* TODO - move to MYNN namelist section
    integer              :: grav_settling     = 0
    integer              :: bl_mynn_tkebudget = 0
    logical              :: bl_mynn_tkeadvect = .false.
    integer              :: bl_mynn_cloudpdf  = 2
    integer              :: bl_mynn_mixlength = 2
    integer              :: bl_mynn_edmf      = 0
    integer              :: bl_mynn_edmf_mom  = 1
    integer              :: bl_mynn_edmf_tke  = 0
    integer              :: bl_mynn_edmf_part = 0
    integer              :: bl_mynn_cloudmix  = 1
    integer              :: bl_mynn_mixqt     = 0
    integer              :: icloud_bl         = 1
    ! *DH
    logical              :: do_myjsfc         = .false.               !< flag for MYJ surface layer scheme
    logical              :: do_myjpbl         = .false.               !< flag for MYJ PBL scheme
#endif
    integer              :: nmtvr          = 14                       !< number of topographic variables such as variance etc
                                                                      !< used in the GWD parameterization
    integer              :: jcap           =  1              !< number of spectral wave trancation used only by sascnv shalcnv
!   real(kind=kind_phys) :: cs_parm(10) = (/5.0,2.5,1.0e3,3.0e3,20.0,-999.,-999.,0.,0.,0./)
    real(kind=kind_phys) :: cs_parm(10) = (/8.0,4.0,1.0e3,3.5e3,20.0,1.0,-999.,1.,0.6,0./)
    real(kind=kind_phys) :: flgmin(2)      = (/0.180,0.220/)          !< [in] ice fraction bounds
    real(kind=kind_phys) :: cgwf(2)        = (/0.5d0,0.05d0/)         !< multiplication factor for convective GWD
    real(kind=kind_phys) :: ccwf(2)        = (/1.0d0,1.0d0/)          !< multiplication factor for critical cloud
                                                                      !< workfunction for RAS
    real(kind=kind_phys) :: cdmbgwd(4)     = (/2.0d0,0.25d0,1.0d0,1.0d0/)   !< multiplication factors for cdmb, gwd, and NS gwd, tke based enhancement
    real(kind=kind_phys) :: sup            = 1.0                      !< supersaturation in pdf cloud (IMP_physics=98) when t is very low
                                                                      !< or ice super saturation in SHOC (when do_shoc=.true.)
    real(kind=kind_phys) :: ctei_rm(2)     = (/10.0d0,10.0d0/)        !< critical cloud top entrainment instability criteria 
                                                                      !< (used if mstrat=.true.)
    real(kind=kind_phys) :: crtrh(3)       = (/0.90d0,0.90d0,0.90d0/) !< critical relative humidity at the surface
                                                                      !< PBL top and at the top of the atmosphere
    real(kind=kind_phys) :: dlqf(2)        = (/0.0d0,0.0d0/)          !< factor for cloud condensate detrainment 
                                                                      !< from cloud edges for RAS
    real(kind=kind_phys) :: psauras(2)     = (/1.0d-3,1.0d-3/)        !< [in] auto conversion coeff from ice to snow in ras
    real(kind=kind_phys) :: prauras(2)     = (/2.0d-3,2.0d-3/)        !< [in] auto conversion coeff from cloud to rain in ras
    real(kind=kind_phys) :: wminras(2)     = (/1.0d-5,1.0d-5/)        !< [in] water and ice minimum threshold for ras

    real(kind=kind_phys) :: rbcr           = 0.25                     !< Critical Richardson Number in PBL scheme
    real(kind=kind_phys) :: shoc_parm(5)   = (/7000.0,1.0,4.2857143,0.7,-999.0/)  !< some tunable parameters for shoc

!--- Rayleigh friction
    real(kind=kind_phys) :: prslrd0        = 0.0d0           !< pressure level from which Rayleigh Damping is applied
    real(kind=kind_phys) :: ral_ts         = 0.0d0           !< time scale for Rayleigh damping in days

!--- mass flux deep convection
    real(kind=kind_phys) :: clam_deep      = 0.1             !< c_e for deep convection (Han and Pan, 2011, eq(6))
    real(kind=kind_phys) :: c0s_deep       = 0.002           !< convective rain conversion parameter
    real(kind=kind_phys) :: c1_deep        = 0.002           !< conversion parameter of detrainment from liquid water into grid-scale cloud water
    real(kind=kind_phys) :: betal_deep     = 0.05            !< fraction factor of downdraft air mass reaching ground surface over land
    real(kind=kind_phys) :: betas_deep     = 0.05            !< fraction factor of downdraft air mass reaching ground surface over sea
    real(kind=kind_phys) :: evfact_deep    = 0.3             !< evaporation factor from convective rain
    real(kind=kind_phys) :: evfactl_deep   = 0.3             !< evaporation factor from convective rain over land
    real(kind=kind_phys) :: pgcon_deep     = 0.55            !< reduction factor in momentum transport due to convection induced pressure gradient force
                                                             !< 0.7 : Gregory et al. (1997, QJRMS)
                                                             !< 0.55: Zhang & Wu (2003, JAS)
    real(kind=kind_phys) :: asolfac_deep   = 0.958           !< aerosol-aware parameter based on Lim (2011)
                                                             !< asolfac= cx / c0s(=.002)
                                                             !< cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s)
                                                             !< Nccn: CCN number concentration in cm^(-3)
                                                             !< Until a realistic Nccn is provided, Nccns are assumed
                                                             !< as Nccn=100 for sea and Nccn=1000 for land 

!--- mass flux shallow convection
    real(kind=kind_phys) :: clam_shal      = 0.3             !< c_e for shallow convection (Han and Pan, 2011, eq(6))
    real(kind=kind_phys) :: c0s_shal       = 0.002           !< conversion parameter of detrainment from liquid water into convetive precipitaiton
    real(kind=kind_phys) :: c1_shal        = 5.e-4           !< conversion parameter of detrainment from liquid water into grid-scale cloud water
    real(kind=kind_phys) :: pgcon_shal     = 0.55            !< reduction factor in momentum transport due to convection induced pressure gradient force
                                                             !< 0.7 : Gregory et al. (1997, QJRMS)
                                                             !< 0.55: Zhang & Wu (2003, JAS)
    real(kind=kind_phys) :: asolfac_shal   = 0.958           !< aerosol-aware parameter based on Lim (2011)
                                                             !< asolfac= cx / c0s(=.002)
                                                             !< cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s)
                                                             !< Nccn: CCN number concentration in cm^(-3)
                                                             !< Until a realistic Nccn is provided, Nccns are assumed
                                                             !< as Nccn=100 for sea and Nccn=1000 for land 

!--- near surface sea temperature model
    logical              :: nst_anl        = .false.         !< flag for NSSTM analysis in gcycle/sfcsub
    integer              :: lsea           = 0
    integer              :: nstf_name(5)   = (/0,0,1,0,5/)   !< flag 0 for no nst  1 for uncoupled nst  and 2 for coupled NST
                                                             !< nstf_name(1) : 0 = NSSTM off, 1 = NSSTM on but uncoupled
                                                             !<                2 = NSSTM on and coupled
                                                             !< nstf_name(2) : 1 = NSSTM spin up on, 0 = NSSTM spin up off
                                                             !< nstf_name(3) : 1 = NSSTM analysis on, 0 = NSSTM analysis off
                                                             !< nstf_name(4) : zsea1 in mm
                                                             !< nstf_name(5) : zsea2 in mm
!--- fractional grid
    logical              :: frac_grid      = .false.         !< flag for fractional grid
    real(kind=kind_phys) :: min_lakeice    = 0.15d0          !< minimum lake ice value
    real(kind=kind_phys) :: min_seaice     = 1.0d-6          !< minimum sea  ice value
    real(kind=kind_phys) :: rho_h2o        = rhowater        !< fresh water density

!--- surface layer z0 scheme
    integer              :: sfc_z0_type    = 0               !< surface roughness options over ocean
                                                             !< 0=no change
                                                             !< 6=areodynamical roughness over water with input 10-m wind
                                                             !< 7=slightly decrease Cd for higher wind speed compare to 6

!--- background vertical diffusion
    real(kind=kind_phys) :: xkzm_m         = 1.0d0           !< [in] bkgd_vdif_m  background vertical diffusion for momentum  
    real(kind=kind_phys) :: xkzm_h         = 1.0d0           !< [in] bkgd_vdif_h  background vertical diffusion for heat q  
    real(kind=kind_phys) :: xkzm_s         = 1.0d0           !< [in] bkgd_vdif_s  sigma threshold for background mom. diffusion  
    real(kind=kind_phys) :: xkzminv        = 0.3             !< diffusivity in inversion layers
    real(kind=kind_phys) :: moninq_fac     = 1.0             !< turbulence diffusion coefficient factor
    real(kind=kind_phys) :: dspfac         = 1.0             !< tke dissipative heating factor
    real(kind=kind_phys) :: bl_upfr        = 0.13            !< updraft fraction in boundary layer mass flux scheme
    real(kind=kind_phys) :: bl_dnfr        = 0.1             !< downdraft fraction in boundary layer mass flux scheme

 
!---Cellular automaton options
    integer              :: nca            = 1
    integer              :: ncells         = 5
    integer              :: nlives         = 10
    real(kind=kind_phys) :: nfracseed      = 0.5
    integer              :: nseed          = 100000
    integer              :: iseed_ca       = 0
    integer              :: nspinup        = 1
    logical              :: do_ca          = .false.
    logical              :: ca_sgs         = .false. 
    logical              :: ca_global      = .false.
    logical              :: ca_smooth      = .false.
    logical              :: isppt_deep     = .false.
    real(kind=kind_phys) :: nthresh        = 0.0
  

!--- IAU options
    real(kind=kind_phys)  :: iau_delthrs      = 0           !< iau time interval (to scale increments)
    character(len=240)    :: iau_inc_files(7) = ''          !< list of increment files
    real(kind=kind_phys)  :: iaufhrs(7)       = -1          !< forecast hours associated with increment files
    logical  :: iau_filter_increments         = .false.     !< filter IAU increments

!--- debug flag
    logical              :: debug          = .false.
    logical              :: pre_rad        = .false.         !< flag for testing purpose

!  max and min lon and lat for critical relative humidity
    integer :: max_lon=5000, max_lat=2000, min_lon=192, min_lat=94
    real(kind=kind_phys) :: rhcmax = 0.9999999               !< max critical rel. hum.

!--- stochastic physics control parameters
    logical :: do_sppt      = .false.
    logical :: use_zmtnblck = .false.
    logical :: do_shum      = .false.
    logical :: do_skeb      = .false.
    integer :: skeb_npass = 11
    logical :: do_sfcperts = .false.   ! mg, sfc-perts
    integer :: nsfcpert    =  6        ! mg, sfc-perts
    real(kind=kind_phys) :: pertz0   = -999.
    real(kind=kind_phys) :: pertzt   = -999.
    real(kind=kind_phys) :: pertshc  = -999.
    real(kind=kind_phys) :: pertlai  = -999.
    real(kind=kind_phys) :: pertalb  = -999.
    real(kind=kind_phys) :: pertvegf = -999.

!--- aerosol scavenging factors
    character(len=20) :: fscav_aero(20) = 'default'

!--- END NAMELIST VARIABLES

    NAMELIST /gfs_physics_nml/                                                              &
                          !--- general parameters
                               fhzero, ldiag3d, lssav, fhcyc,                               &
                               thermodyn_id, sfcpress_id,                                   &
                          !--- coupling parameters
                               cplflx, cplwav, cplchm, lsidea,                              &
                          !--- radiation parameters
                               fhswr, fhlwr, levr, nfxr, aero_in, iflip, isol, ico2, ialb,  &
                               isot, iems, iaer, icliq_sw, iovr_sw, iovr_lw, ictm, isubc_sw,&
                               isubc_lw, crick_proof, ccnorm, lwhtr, swhtr,                 &
#ifdef CCPP
                               nhfrad,                                                      &
#endif
                          ! IN CCN forcing
                               iccn,                                                        &
                          !--- microphysical parameterizations
                               ncld, imp_physics, psautco, prautco, evpco, wminco,          &
#ifdef CCPP
                               fprcp, pdfflag, mg_dcs, mg_qcvar, mg_ts_auto_ice, mg_rhmini, &
                               effr_in, tf, tcr,                                            &
#else
                               fprcp, pdfflag, mg_dcs, mg_qcvar, mg_ts_auto_ice, effr_in,   &
#endif
                               microp_uniform, do_cldice, hetfrz_classnuc,                  &
                               mg_do_graupel, mg_do_hail, mg_nccons, mg_nicons, mg_ngcons,  &
                               mg_ncnst, mg_ninst, mg_ngnst, sed_supersat, do_sb_physics,   &
                               mg_alf,   mg_qcmin, mg_do_ice_gmao, mg_do_liq_liu,           &
                               ltaerosol, lradar, lrefres, ttendlim, lgfdlmprad,            &
                          !--- max hourly
                               avg_max_length,                                              &
                          !--- land/surface model control
#ifdef CCPP
                               lsm, lsoil, lsoil_lsm, lsnow_lsm, rdlai,                     &
                               nmtvr, ivegsrc, use_ufo,                                     &
#else
                               lsm, lsoil, nmtvr, ivegsrc, use_ufo,                         &
#endif
                          !    Noah MP options
                               iopt_dveg,iopt_crs,iopt_btr,iopt_run,iopt_sfc, iopt_frz,     &
                               iopt_inf, iopt_rad,iopt_alb,iopt_snf,iopt_tbot,iopt_stc,     &
                          !--- physical parameterizations
                               ras, trans_trac, old_monin, cnvgwd, mstrat, moist_adj,       &
                               cscnv, cal_pre, do_aw, do_shoc, shocaftcnv, shoc_cld,        &
#ifdef CCPP
                               oz_phys, oz_phys_2015,                                       &
                               do_mynnedmf, do_mynnsfclay,                                  &
                               ! DH* TODO - move to MYNN namelist section
                               bl_mynn_cloudpdf, bl_mynn_edmf, bl_mynn_edmf_mom,            &
                               bl_mynn_edmf_tke, bl_mynn_edmf_part, bl_mynn_cloudmix,       &
                               bl_mynn_mixqt, icloud_bl, bl_mynn_tkeadvect, gwd_opt,        &
                               ! *DH
                               do_myjsfc, do_myjpbl,                                        &
#endif
                               h2o_phys, pdfcld, shcnvcw, redrag, hybedmf, satmedmf,        &
                               shinhong, do_ysu, dspheat, lheatstrg, cnvcld,                &
                               random_clds, shal_cnv, imfshalcnv, imfdeepcnv, isatmedmf,    &
                               do_deep, jcap,                                               &
                               cs_parm, flgmin, cgwf, ccwf, cdmbgwd, sup, ctei_rm, crtrh,   &
                               dlqf, rbcr, shoc_parm, psauras, prauras, wminras,            &
                               do_sppt, do_shum, do_skeb, do_sfcperts,                      &
                          !--- Rayleigh friction
                               prslrd0, ral_ts,  ldiag_ugwp, do_ugwp, do_tofd,              &
                          ! --- Ferrier-Aligo
#ifdef CCPP
                               spec_adv, rhgrd,                                             &
#endif
                          !--- mass flux deep convection
                               clam_deep, c0s_deep, c1_deep, betal_deep,                    &
                               betas_deep, evfact_deep, evfactl_deep, pgcon_deep,           &
                               asolfac_deep,                                                &
                          !--- mass flux shallow convection
                               clam_shal, c0s_shal, c1_shal, pgcon_shal, asolfac_shal,      &
                          !--- near surface sea temperature model
                               nst_anl, lsea, nstf_name,                                    &
                               frac_grid, min_lakeice, min_seaice,                          &
                               frac_grid,                                                   &
                          !--- surface layer
                               sfc_z0_type,                                                 &
                          !    background vertical diffusion
                               xkzm_m, xkzm_h, xkzm_s, xkzminv, moninq_fac, dspfac,         &
                               bl_upfr, bl_dnfr,                                            &
                          !--- cellular automata
                               nca, ncells, nlives, nfracseed,nseed, nthresh, do_ca,        &
                               ca_sgs, ca_global,iseed_ca,ca_smooth,isppt_deep,nspinup,     &
                          !--- IAU
                               iau_delthrs,iaufhrs,iau_inc_files,iau_filter_increments,     &
                          !--- debug options
                               debug, pre_rad,                                              &
                          !--- parameter range for critical relative humidity
                               max_lon, max_lat, min_lon, min_lat, rhcmax,                  &
                               phys_version,                                                &
                          !--- aerosol scavenging factors ('name:value' string array)
                               fscav_aero

!--- other parameters 
    integer :: nctp    =  0                !< number of cloud types in CS scheme
    logical :: gen_coord_hybrid = .false.  !< for Henry's gen coord

!--- SHOC parameters
    integer :: nshoc_2d  = 0  !< number of 2d fields for SHOC
    integer :: nshoc_3d  = 0  !< number of 3d fields for SHOC

!--- convective clouds
    integer :: ncnvcld3d = 0       !< number of convective 3d clouds fields


!--- read in the namelist
#ifdef INTERNAL_FILE_NML
    Model%input_nml_file => input_nml_file
    read(Model%input_nml_file, nml=gfs_physics_nml)
#ifdef CCPP
    ! Set length (number of lines) in namelist for internal reads
    Model%input_nml_file_length = size(Model%input_nml_file)
#endif
#else
    inquire (file=trim(fn_nml), exist=exists)
    if (.not. exists) then
      write(6,*) 'GFS_namelist_read:: namelist file: ',trim(fn_nml),' does not exist'
      stop
    else
      open (unit=nlunit, file=fn_nml, action='READ', status='OLD', iostat=ios)
    endif
    rewind(nlunit)
    read (nlunit, nml=gfs_physics_nml)
    close (nlunit)
#ifdef CCPP
    ! Set length (number of lines) in namelist for internal reads
    Model%input_nml_file_length = 0
#endif
#endif
!--- write version number and namelist to log file ---
    if (me == master) then
      write(logunit, '(a80)') '================================================================================'
      write(logunit, '(a64)') phys_version
      write(logunit, nml=gfs_physics_nml)
    endif

!--- MPI parameters
    Model%me               = me
    Model%master           = master
#ifdef CCPP
    Model%communicator     = communicator
    Model%ntasks           = ntasks
    Model%nthreads         = nthreads
#endif
    Model%nlunit           = nlunit
    Model%fn_nml           = fn_nml
#ifdef CCPP
    Model%logunit          = logunit
#endif
    Model%fhzero           = fhzero
    Model%ldiag3d          = ldiag3d
!
!VAY-ugwp  --- set some GW-related switches
!
    Model%ldiag_ugwp       = ldiag_ugwp
    Model%do_ugwp          = do_ugwp
    Model%do_tofd          = do_tofd

    Model%lssav            = lssav
    Model%fhcyc            = fhcyc
    Model%thermodyn_id     = thermodyn_id
    Model%sfcpress_id      = sfcpress_id
    Model%gen_coord_hybrid = gen_coord_hybrid

    !--- set some grid extent parameters
    Model%tile_num         = tile_num
    Model%isc              = isc
    Model%jsc              = jsc
    Model%nx               = nx
    Model%ny               = ny
    Model%levs             = levs
#ifdef CCPP
    Model%levsp1           = Model%levs + 1
    Model%levsm1           = Model%levs - 1
    allocate(Model%ak(1:size(ak)))
    allocate(Model%bk(1:size(bk)))
    Model%ak               = ak
    Model%bk               = bk
#endif
    Model%cnx              = cnx
    Model%cny              = cny
    Model%lonr             = gnx         ! number longitudinal points
    Model%latr             = gny         ! number of latitudinal points from pole to pole
    Model%nblks            = size(blksz)
    allocate(Model%blksz(1:Model%nblks))
    Model%blksz            = blksz
#ifdef CCPP
    allocate(Model%blksz2(1:Model%nblks))
    Model%blksz2           = blksz
#endif

!--- coupling parameters
    Model%cplflx           = cplflx
    Model%cplwav           = cplwav
    Model%cplchm           = cplchm

!--- integrated dynamics through earth's atmosphere
    Model%lsidea           = lsidea

!--- calendars and time parameters and activation triggers
    Model%dtp              = dt_phys
    Model%dtf              = dt_dycore
    Model%nscyc            = nint(Model%fhcyc*con_hr/Model%dtp)
    Model%nszero           = nint(Model%fhzero*con_hr/Model%dtp)
    Model%idat(1:8)        = idat(1:8)
    Model%idate            = 0
    Model%idate(1)         = Model%idat(5)
    Model%idate(2)         = Model%idat(2)
    Model%idate(3)         = Model%idat(3)
    Model%idate(4)         = Model%idat(1)
    Model%iau_offset       = iau_offset

!--- radiation control parameters
    Model%fhswr            = fhswr
    Model%fhlwr            = fhlwr
    Model%nsswr            = nint(fhswr/Model%dtp)
    Model%nslwr            = nint(fhlwr/Model%dtp)
#ifdef CCPP
    if (restart) then
      Model%nhfrad         = 0
      if (Model%me == Model%master .and. nhfrad>0) &
        write(*,'(a)') 'Disable high-frequency radiation calls for restart run'
    else
      Model%nhfrad         = nhfrad
      if (Model%me == Model%master .and. nhfrad>0) &
        write(*,'(a,i0)') 'Number of high-frequency radiation calls for coldstart run: ', nhfrad
    endif
#endif
    if (levr < 0) then
      Model%levr           = levs
    else
      Model%levr           = levr
    endif
#ifdef CCPP
    Model%levrp1           = Model%levr + 1
#endif
    Model%nfxr             = nfxr
    Model%aero_in          = aero_in
    if (Model%aero_in) then
      ntrcaer = ntrcaerm
    else
      ntrcaer = 1
    endif
#ifdef CCPP
    Model%ntrcaer          = ntrcaer
#endif
    Model%iccn             = iccn
    if (Model%aero_in) Model%iccn = .false.
    ! further down: set Model%iccn to .false.
    ! for all microphysics schemes except
    ! MG2/3 (these are the only ones using ICCN)
    Model%iflip            = iflip
    Model%isol             = isol
    Model%ico2             = ico2
    Model%ialb             = ialb
    Model%iems             = iems
    Model%iaer             = iaer
    Model%icliq_sw         = icliq_sw
    Model%iovr_sw          = iovr_sw
    Model%iovr_lw          = iovr_lw
    Model%ictm             = ictm
    Model%isubc_sw         = isubc_sw
    Model%isubc_lw         = isubc_lw
    Model%crick_proof      = crick_proof
    Model%ccnorm           = ccnorm
    Model%lwhtr            = lwhtr
    Model%swhtr            = swhtr
#ifdef CCPP
    ! The CCPP versions of the RRTMG lw/sw schemes are configured
    ! such that lw and sw heating rate are output, i.e. they rely
    ! on the corresponding arrays to be allocated.
    if (.not.lwhtr .or. .not.swhtr) then
      write(0,*) "Logic error, the CCPP version of RRTMG lwrad/swrad require the output" // &
             " of the lw/sw heating rates to be turned on (namelist options lwhtr and swhtr)"
      stop
    end if
#endif

!--- microphysical switch
    Model%ncld             = ncld
    Model%imp_physics      = imp_physics
    ! turn off ICCN interpolation when MG2/3 are not used
    if (.not. Model%imp_physics==Model%imp_physics_mg) Model%iccn = .false.
!--- Zhao-Carr MP parameters
    Model%psautco          = psautco
    Model%prautco          = prautco
    Model%evpco            = evpco
    Model%wminco           = wminco
!--- Max hourly
    Model%avg_max_length   = avg_max_length
!--- Morrison-Gettelman MP parameters
    Model%fprcp            = fprcp
    Model%pdfflag          = pdfflag
    Model%mg_dcs           = mg_dcs
    Model%mg_qcvar         = mg_qcvar
    Model%mg_ts_auto_ice   = mg_ts_auto_ice
#ifdef CCPP
    Model%mg_rhmini        = mg_rhmini
#endif
    Model%mg_alf           = mg_alf
    Model%mg_qcmin         = mg_qcmin
    Model%effr_in          = effr_in
    Model%microp_uniform   = microp_uniform
    Model%do_cldice        = do_cldice
    Model%hetfrz_classnuc  = hetfrz_classnuc
    Model%mg_do_graupel    = mg_do_graupel
    Model%mg_do_hail       = mg_do_hail
    Model%mg_do_ice_gmao   = mg_do_ice_gmao
    Model%mg_do_liq_liu    = mg_do_liq_liu
    Model%mg_nccons        = mg_nccons
    Model%mg_nicons        = mg_nicons
    Model%mg_ngcons        = mg_ngcons
    Model%mg_ncnst         = mg_ncnst
    Model%mg_ninst         = mg_ninst
    Model%mg_ngnst         = mg_ngnst
    Model%sed_supersat     = sed_supersat
    Model%do_sb_physics    = do_sb_physics
    Model%mg_precip_frac_method  = mg_precip_frac_method
    Model%mg_berg_eff_factor     = mg_berg_eff_factor
#ifdef CCPP
    Model%tf               = tf
    Model%tcr              = tcr
    Model%tcrf             = 1.0/(tcr-tf)
#endif

!--- Thompson MP parameters
    Model%ltaerosol        = ltaerosol
    Model%lradar           = lradar
    Model%ttendlim         = ttendlim
!--- F-A MP parameters
#ifdef CCPP
    Model%rhgrd            = rhgrd
    Model%spec_adv         = spec_adv
#endif 

!--- gfdl  MP parameters
    Model%lgfdlmprad       = lgfdlmprad
!--- Thompson,GFDL MP parameter
    Model%lrefres          = lrefres

!--- land/surface model parameters
    Model%lsm              = lsm
    Model%lsoil            = lsoil
#ifdef CCPP
    ! Consistency check for RUC LSM
    if (Model%lsm == Model%lsm_ruc .and. Model%nscyc>0) then
      write(0,*) 'Logic error: RUC LSM cannot be used with surface data cycling at this point (fhcyc>0)'
      stop
    end if
    ! Flag to read leaf area index from input files (initial conditions)
    Model%rdlai = rdlai
    if (Model%rdlai .and. .not. Model%lsm == Model%lsm_ruc) then
      write(0,*) 'Logic error: rdlai = .true. only works with RUC LSM'
      stop
    end if
    ! Set surface layers for CCPP physics
    if (lsoil_lsm==-1) then
      Model%lsoil_lsm      = lsoil
    else
      Model%lsoil_lsm      = lsoil_lsm
    end if
    if (lsnow_lsm /= 3) then
      write(0,*) 'Logic error: NoahMP expects the maximum number of snow layers to be exactly 3 (see sfc_noahmp_drv.f)'
      stop
    else
      Model%lsnow_lsm        = lsnow_lsm
      ! Set lower bound for LSM model, runs from negative (above surface) to surface (zero)
      Model%lsnow_lsm_lbound = -Model%lsnow_lsm+1
    end if
#endif
    Model%ivegsrc          = ivegsrc
    Model%isot             = isot
    Model%use_ufo          = use_ufo

! Noah MP options from namelist
!
    Model%iopt_dveg        = iopt_dveg
    Model%iopt_crs         = iopt_crs
    Model%iopt_btr         = iopt_btr
    Model%iopt_run         = iopt_run
    Model%iopt_sfc         = iopt_sfc
    Model%iopt_frz         = iopt_frz
    Model%iopt_inf         = iopt_inf
    Model%iopt_rad         = iopt_rad
    Model%iopt_alb         = iopt_alb
    Model%iopt_snf         = iopt_snf
    Model%iopt_tbot        = iopt_tbot
    Model%iopt_stc         = iopt_stc

!--- tuning parameters for physical parameterizations
    Model%ras              = ras
    Model%flipv            = flipv
    Model%trans_trac       = trans_trac
    Model%old_monin        = old_monin
    Model%cnvgwd           = cnvgwd
    Model%mstrat           = mstrat
    Model%moist_adj        = moist_adj
    Model%cscnv            = cscnv
    Model%cal_pre          = cal_pre
    Model%do_aw            = do_aw
    Model%cs_parm          = cs_parm
    Model%do_shoc          = do_shoc
#ifdef CCPP
    if (Model%do_shoc) then
      print *, "Error, update of SHOC from May 22 2019 not yet in CCPP"
      stop
    end if
#endif
    Model%shoc_parm        = shoc_parm
    Model%shocaftcnv       = shocaftcnv
    Model%shoc_cld         = shoc_cld
#ifdef CCPP
    if (oz_phys .and. oz_phys_2015) then
       write(*,*) 'Logic error: can only use one ozone physics option (oz_phys or oz_phys_2015), not both. Exiting.'
       stop
    end if
    Model%oz_phys          = oz_phys
    Model%oz_phys_2015     = oz_phys_2015
#endif
    Model%h2o_phys         = h2o_phys
#ifdef CCPP
    ! To ensure that these values match what's in the physics,
    ! array sizes are compared during model init in GFS_phys_time_vary_init()
    !
    ! from module h2ointerp
    if (h2o_phys) then
       levh2o    = 72
       h2o_coeff = 3
    else
       levh2o    = 1
       h2o_coeff = 1
    end if
#endif
    Model%pdfcld            = pdfcld
    Model%shcnvcw           = shcnvcw
    Model%redrag            = redrag
    Model%hybedmf           = hybedmf
    Model%satmedmf          = satmedmf
    Model%shinhong          = shinhong
    Model%do_ysu            = do_ysu
    Model%dspheat           = dspheat
    Model%lheatstrg         = lheatstrg
    Model%cnvcld            = cnvcld
    Model%random_clds       = random_clds
    Model%shal_cnv          = shal_cnv
    Model%imfshalcnv        = imfshalcnv
    Model%imfdeepcnv        = imfdeepcnv
    Model%isatmedmf         = isatmedmf
    Model%do_deep           = do_deep
    Model%nmtvr             = nmtvr
    Model%jcap              = jcap
    Model%flgmin            = flgmin
    Model%cgwf              = cgwf
    Model%ccwf              = ccwf
    Model%cdmbgwd           = cdmbgwd
    Model%sup               = sup
    Model%ctei_rm           = ctei_rm
    Model%crtrh             = crtrh
    Model%dlqf              = dlqf
    Model%psauras           = psauras
    Model%prauras           = prauras
    Model%wminras           = wminras
    Model%rbcr              = rbcr
    Model%do_gwd            = maxval(Model%cdmbgwd) > 0.0
    Model%do_cnvgwd         = Model%cnvgwd .and. maxval(Model%cdmbgwd(3:4)) == 0.0
#ifdef CCPP
    Model%do_mynnedmf       = do_mynnedmf
    Model%do_mynnsfclay     = do_mynnsfclay
    ! DH* TODO - move to MYNN namelist section
    Model%bl_mynn_cloudpdf  = bl_mynn_cloudpdf
    Model%bl_mynn_mixlength = bl_mynn_mixlength
    Model%bl_mynn_edmf      = bl_mynn_edmf
    Model%bl_mynn_edmf_mom  = bl_mynn_edmf_mom
    Model%bl_mynn_edmf_tke  = bl_mynn_edmf_tke
    Model%bl_mynn_cloudmix  = bl_mynn_cloudmix
    Model%bl_mynn_mixqt     = bl_mynn_mixqt
    Model%bl_mynn_edmf_part = bl_mynn_edmf_part
    Model%bl_mynn_tkeadvect = bl_mynn_tkeadvect
    Model%grav_settling     = grav_settling
    Model%icloud_bl         = icloud_bl
    ! *DH
    Model%gwd_opt           = gwd_opt
    Model%do_myjsfc         = do_myjsfc
    Model%do_myjpbl         = do_myjpbl
#endif

!--- Rayleigh friction
    Model%prslrd0          = prslrd0
    Model%ral_ts           = ral_ts

!--- mass flux deep convection
    Model%clam_deep        = clam_deep
    Model%c0s_deep         = c0s_deep
    Model%c1_deep          = c1_deep
    Model%betal_deep       = betal_deep
    Model%betas_deep       = betas_deep
    Model%evfact_deep      = evfact_deep
    Model%evfactl_deep     = evfactl_deep
    Model%pgcon_deep       = pgcon_deep
    Model%asolfac_deep     = asolfac_deep

!--- mass flux shallow convection
    Model%clam_shal        = clam_shal
    Model%c0s_shal         = c0s_shal
    Model%c1_shal          = c1_shal
    Model%pgcon_shal       = pgcon_shal
    Model%asolfac_shal     = asolfac_shal

!--- near surface sea temperature model
    Model%nst_anl          = nst_anl
    Model%lsea             = lsea
    Model%nstf_name        = nstf_name

!--- fractional grid
    Model%frac_grid        = frac_grid
#ifdef CCPP
    if (Model%frac_grid) then
      write(0,*) "ERROR: CCPP has not been tested with fractional landmask turned on"
      stop
    end if
#endif
    Model%min_lakeice      = min_lakeice
    Model%min_seaice       = min_seaice
    Model%rho_h2o          = rho_h2o

!--- surface layer
    Model%sfc_z0_type      = sfc_z0_type

!--- backgroud vertical diffusion
    Model%xkzm_m           = xkzm_m
    Model%xkzm_h           = xkzm_h
    Model%xkzm_s           = xkzm_s
    Model%xkzminv          = xkzminv
    Model%moninq_fac       = moninq_fac
    Model%dspfac           = dspfac
    Model%bl_upfr          = bl_upfr
    Model%bl_dnfr          = bl_dnfr

!--- stochastic physics options
    ! do_sppt, do_shum, do_skeb and do_sfcperts are namelist variables in group
    ! physics that are parsed here and then compared in init_stochastic_physics
    ! to the stochastic physics namelist parametersto ensure consistency.
    Model%do_sppt          = do_sppt
    Model%use_zmtnblck     = use_zmtnblck
    Model%do_shum          = do_shum
    Model%do_skeb          = do_skeb
    Model%do_sfcperts      = do_sfcperts ! mg, sfc-perts
    Model%nsfcpert         = nsfcpert    ! mg, sfc-perts
    Model%pertz0           = pertz0
    Model%pertzt           = pertzt
    Model%pertshc          = pertshc
    Model%pertlai          = pertlai
    Model%pertalb          = pertalb
    Model%pertvegf         = pertvegf

    !--- cellular automata options
    Model%nca              = nca
    Model%ncells           = ncells
    Model%nlives           = nlives
    Model%nfracseed        = nfracseed
    Model%nseed            = nseed
    Model%ca_global        = ca_global
    Model%do_ca            = do_ca
    Model%ca_sgs           = ca_sgs
    Model%iseed_ca         = iseed_ca
    Model%ca_smooth        = ca_smooth
    Model%isppt_deep       = isppt_deep
    Model%nspinup          = nspinup  
    Model%nthresh          = nthresh 

    ! IAU flags
    !--- iau parameters
    Model%iaufhrs         = iaufhrs
    Model%iau_inc_files   = iau_inc_files
    Model%iau_delthrs     = iau_delthrs
    Model%iau_filter_increments = iau_filter_increments
    if(Model%me==0) print *,' model init,iaufhrs=',Model%iaufhrs

!--- tracer handling
    Model%ntrac            = size(tracer_names)
#ifdef CCPP
    Model%ntracp1          = Model%ntrac + 1
#endif
    allocate (Model%tracer_names(Model%ntrac))
    Model%tracer_names(:)  = tracer_names(:)
#ifdef CCPP
    Model%ntqv             = 1
#endif
#ifdef MULTI_GASES
    Model%nto              = get_tracer_index(Model%tracer_names, 'spfo',        Model%me, Model%master, Model%debug)
    Model%nto2             = get_tracer_index(Model%tracer_names, 'spfo2',       Model%me, Model%master, Model%debug)
    Model%ntoz             = get_tracer_index(Model%tracer_names, 'spfo3',       Model%me, Model%master, Model%debug)
#else
    Model%ntoz             = get_tracer_index(Model%tracer_names, 'o3mr',       Model%me, Model%master, Model%debug)
#endif
    Model%ntcw             = get_tracer_index(Model%tracer_names, 'liq_wat',    Model%me, Model%master, Model%debug)
    Model%ntiw             = get_tracer_index(Model%tracer_names, 'ice_wat',    Model%me, Model%master, Model%debug)
    Model%ntrw             = get_tracer_index(Model%tracer_names, 'rainwat',    Model%me, Model%master, Model%debug)
    Model%ntsw             = get_tracer_index(Model%tracer_names, 'snowwat',    Model%me, Model%master, Model%debug)
    Model%ntgl             = get_tracer_index(Model%tracer_names, 'graupel',    Model%me, Model%master, Model%debug)
    Model%ntclamt          = get_tracer_index(Model%tracer_names, 'cld_amt',    Model%me, Model%master, Model%debug)
    Model%ntlnc            = get_tracer_index(Model%tracer_names, 'water_nc',   Model%me, Model%master, Model%debug)
    Model%ntinc            = get_tracer_index(Model%tracer_names, 'ice_nc',     Model%me, Model%master, Model%debug)
    Model%ntrnc            = get_tracer_index(Model%tracer_names, 'rain_nc',    Model%me, Model%master, Model%debug)
    Model%ntsnc            = get_tracer_index(Model%tracer_names, 'snow_nc',    Model%me, Model%master, Model%debug)
    Model%ntgnc            = get_tracer_index(Model%tracer_names, 'graupel_nc', Model%me, Model%master, Model%debug)
    Model%ntke             = get_tracer_index(Model%tracer_names, 'sgs_tke',    Model%me, Model%master, Model%debug)
#ifdef CCPP
    Model%nqrimef          = get_tracer_index(Model%tracer_names, 'q_rimef',    Model%me, Model%master, Model%debug)
#endif
    Model%ntwa             = get_tracer_index(Model%tracer_names, 'liq_aero',   Model%me, Model%master, Model%debug)
    Model%ntia             = get_tracer_index(Model%tracer_names, 'ice_aero',   Model%me, Model%master, Model%debug)
    Model%ntchm            = 0
    Model%ntchs            = get_tracer_index(Model%tracer_names, 'so2',        Model%me, Model%master, Model%debug)
    if (Model%ntchs > 0) then
      Model%ntchm          = get_tracer_index(Model%tracer_names, 'pp10',       Model%me, Model%master, Model%debug)
      if (Model%ntchm > 0) then
        Model%ntchm = Model%ntchm - Model%ntchs + 1
        allocate(Model%ntdiag(Model%ntchm))
        ! -- turn on all tracer diagnostics to .true. by default, except for so2
        Model%ntdiag(1)  = .false.
        Model%ntdiag(2:) = .true.
        ! -- turn off diagnostics for DMS
        n = get_tracer_index(Model%tracer_names, 'DMS', Model%me, Model%master, Model%debug) - Model%ntchs + 1
        if (n > 0) Model%ntdiag(n) = .false.
        ! -- turn off diagnostics for msa
        n = get_tracer_index(Model%tracer_names, 'msa', Model%me, Model%master, Model%debug) - Model%ntchs + 1
        if (n > 0) Model%ntdiag(n) = .false.
      endif
    endif

    ! -- setup aerosol scavenging factors
    allocate(Model%fscav(Model%ntchm))
    if (Model%ntchm > 0) then
      ! -- initialize to default
      Model%fscav = 0.6_kind_phys
      n = get_tracer_index(Model%tracer_names, 'seas1', Model%me, Model%master, Model%debug) - Model%ntchs + 1
      if (n > 0) Model%fscav(n) = 1.0_kind_phys
      n = get_tracer_index(Model%tracer_names, 'seas2', Model%me, Model%master, Model%debug) - Model%ntchs + 1
      if (n > 0) Model%fscav(n) = 1.0_kind_phys
      n = get_tracer_index(Model%tracer_names, 'seas3', Model%me, Model%master, Model%debug) - Model%ntchs + 1
      if (n > 0) Model%fscav(n) = 1.0_kind_phys
      n = get_tracer_index(Model%tracer_names, 'seas4', Model%me, Model%master, Model%debug) - Model%ntchs + 1
      if (n > 0) Model%fscav(n) = 1.0_kind_phys
      n = get_tracer_index(Model%tracer_names, 'seas5', Model%me, Model%master, Model%debug) - Model%ntchs + 1
      if (n > 0) Model%fscav(n) = 1.0_kind_phys
      ! -- read factors from namelist
      do i = 1, size(fscav_aero)
        j = index(fscav_aero(i),":")
        if (j > 1) then
          read(fscav_aero(i)(j+1:), *, iostat=ios) tem
          if (ios /= 0) cycle
          if (adjustl(fscav_aero(i)(:j-1)) == "*") then
            Model%fscav = tem
            exit
          else
            n = get_tracer_index(Model%tracer_names, adjustl(fscav_aero(i)(:j-1)), Model%me, Model%master, Model%debug) &
                - Model%ntchs + 1
            if (n > 0) Model%fscav(n) = tem
          endif
        endif
      enddo
    endif

#ifdef CCPP
    ! To ensure that these values match what's in the physics,
    ! array sizes are compared during model init in GFS_phys_time_vary_init()
    !
    ! from module ozinterp
    if (Model%ntoz>0) then
       if (Model%oz_phys) then
          levozp   = 80
          oz_coeff = 4
       else if (Model%oz_phys_2015) then
          levozp   = 53
          oz_coeff = 6
       else
          write(*,*) 'Logic error, ntoz>0 but no ozone physics selected'
          stop
       end if
    else
       if (Model%oz_phys .or. Model%oz_phys_2015) then
          write(*,*) 'Logic error, ozone physics are selected, but ntoz<=0'
          stop
       else
          levozp   = 1
          oz_coeff = 0
       end if
    end if
#endif

!--- quantities to be used to derive phy_f*d totals
    Model%nshoc_2d         = nshoc_2d
    Model%nshoc_3d         = nshoc_3d
    Model%ncnvcld3d        = ncnvcld3d
    Model%nctp             = nctp

!--- debug flag
    Model%debug            = debug
    Model%pre_rad          = pre_rad

!--- set initial values for time varying properties
    Model%ipt              = 1
    Model%lprnt            = .false.
    Model%lsswr            = .false.
    Model%lslwr            = .false.
    Model%solhr            = -9999.
    Model%solcon           = -9999.
    Model%slag             = -9999.
    Model%sdec             = -9999.
    Model%cdec             = -9999.
    Model%clstp            = -9999
    rinc(1:5)              = 0 
    call w3difdat(jdat,idat,4,rinc)
    Model%phour            = rinc(4)/con_hr
    Model%fhour            = (rinc(4) + Model%dtp)/con_hr
    Model%zhour            = mod(Model%phour,Model%fhzero)
    Model%kdt              = 0
#ifdef CCPP
    Model%first_time_step  = .true.
    Model%restart          = restart
    Model%hydrostatic      = hydrostatic
#endif
    Model%jdat(1:8)        = jdat(1:8)
#ifdef CCPP
    Model%sec              = 0
    if (Model%lsm == Model%lsm_noahmp) then
      Model%yearlen          = 365
      Model%julian           = -9999.
    endif
    ! DH* what happens if LTP>0? Does this have to change? 
    ! A conversation with Yu-Tai suggests that we can probably
    ! eliminate LTP altogether *DH
    allocate(Model%si(Model%levr+1))
    !--- Define sigma level for radiation initialization
    !--- The formula converting hybrid sigma pressure coefficients to sigma coefficients follows Eckermann (2009, MWR)
    !--- ps is replaced with p0. The value of p0 uses that in http://www.emc.ncep.noaa.gov/officenotes/newernotes/on461.pdf
    !--- ak/bk have been flipped from their original FV3 orientation and are defined sfc -> toa
    Model%si = (ak + bk * con_p0 - ak(Model%levr+1)) / (con_p0 - ak(Model%levr+1))
#endif

#ifndef CCPP
    ! Beware! The values set here reside in wam_f107_kp_mod and determine sizes of arrays
    ! inside that module. These arrays get used later in modules idea_tracer.f, idea_ion.f,
    ! idea_solar_heating.f, efield.f, and idea_composition.f.
    ! Since in wam_f107_kp_mod no default values are assigned to the four integers below, not
    ! setting them here can lead to memory corruption that is hard to detect.
!--- stored in wam_f107_kp module
    f107_kp_size      = 56
    f107_kp_skip_size = 0
    f107_kp_data_size = 56
    f107_kp_interval  = 10800
#endif

!--- BEGIN CODE FROM GFS_PHYSICS_INITIALIZE
!--- define physcons module variables
    tem     = con_rerth*con_rerth*(con_pi+con_pi)*con_pi
#ifdef CCPP
    Model%dxmax  = log(tem/(max_lon*max_lat))
    Model%dxmin  = log(tem/(min_lon*min_lat))
    Model%dxinv  = 1.0d0 / (Model%dxmax-Model%dxmin)
    Model%rhcmax = rhcmax
    if (Model%me == Model%master) write(*,*)' dxmax=',Model%dxmax,' dxmin=',Model%dxmin,' dxinv=',Model%dxinv, &
       'max_lon=',max_lon,' max_lat=',max_lat,' min_lon=',min_lon,' min_lat=',min_lat,       &
       ' rhc_max=',Model%rhcmax
#else
    dxmax   = log(tem/(max_lon*max_lat))
    dxmin   = log(tem/(min_lon*min_lat))
    dxinv   = 1.0d0 / (dxmax-dxmin)
    rhc_max = rhcmax
    if (Model%me == Model%master) write(*,*)' dxmax=',dxmax,' dxmin=',dxmin,' dxinv=',dxinv, &
       'max_lon=',max_lon,' max_lat=',max_lat,' min_lon=',min_lon,' min_lat=',min_lat,       &
       ' rhc_max=',rhc_max
#endif

!--- set nrcm 

#ifndef CCPP
    if (Model%ras) then
      Model%nrcm = min(nrcmax, Model%levs-1) * (Model%dtp/1200.d0) + 0.10001d0
    else
      Model%nrcm = 2
    endif
#else
    Model%nrcm = 2
#endif

!--- cal_pre
    if (Model%cal_pre) then
      Model%random_clds = .true.
    endif
!--- END CODE FROM GFS_PHYSICS_INITIALIZE


!--- BEGIN CODE FROM COMPNS_PHYSICS
!--- shoc scheme
    if (do_shoc) then
      Model%nshoc_3d   = 3
      Model%nshoc_2d   = 0
      Model%shal_cnv   = .false.
      Model%imfshalcnv = -1
      Model%hybedmf    = .false.
      Model%satmedmf   = .false.
      if (Model%me == Model%master) print *,' Simplified Higher Order Closure Model used for', &
                                            ' Boundary layer and Shallow Convection',          &
                                            ' nshoc_3d=',Model%nshoc_3d,                       &
                                            ' nshoc_2d=',Model%nshoc_2d,                       &
                                            ' ntke=',Model%ntke,' shoc_parm=',shoc_parm
    endif

#ifdef CCPP
    !--- mynn-edmf scheme
    if (Model%do_mynnedmf) then
      if (Model%do_shoc .or. Model%hybedmf .or. Model%satmedmf) then
          print *,' Logic error: MYNN EDMF cannot be run with SHOC, HEDMF or SATMEDMF'
          stop
      end if
!      Model%shal_cnv   = .false.
!      Model%imfshalcnv = -1
      ! DH* substitute for MYNN namelist section
      Model%icloud_bl         = 1
      !Model%bl_mynn_tkeadvect = .true.
      Model%bl_mynn_edmf      = 1
      !Model%bl_mynn_edmf_mom  = 1
      ! *DH
      if (Model%me == Model%master) print *,' MYNN-EDMF scheme is used for both',                &
                                            ' boundary layer turbulence and shallow convection', &
                                            ' bl_mynn_cloudpdf=',Model%bl_mynn_cloudpdf,         &
                                            ' bl_mynn_mixlength=',Model%bl_mynn_mixlength,       &
                                            ' bl_mynn_edmf=',Model%bl_mynn_edmf
    endif
#endif

!--- set number of cloud types
    if (Model%cscnv) then
      Model%nctp = nint(Model%cs_parm(5))
      Model%nctp = max(Model%nctp,10)
      if (Model%cs_parm(7) < 0.0) Model%cs_parm(7) = Model%dtp
      Model%do_awdd  = Model%do_aw .and. Model%cs_parm(6) > 0.0
!     Model%flx_form = Model%do_aw .and. Model%cs_parm(8) > 0.0
      Model%flx_form = Model%cs_parm(8) > 0.0
    endif
!   Model%nctp = max(Model%nctp,1)

!--- output information about the run
    if (Model%me == Model%master) then
      if (Model%lsm == 1) then
        print *,' NOAH Land Surface Model used'
      elseif (Model%lsm == 0) then
        print *,' OSU no longer supported - job aborted'
        stop
      elseif (Model%lsm == Model%lsm_noahmp) then
        if (Model%ivegsrc /= 1) then
          print *,'Vegetation type must be IGBP if Noah MP is used'
          stop
        elseif (Model%isot /= 1) then
          print *,'Soil type must be STATSGO if Noah MP is used'
          stop
        endif
        print *, 'New Noah MP Land Surface Model will be used'
        print *, 'The Physics options are'

        print *,'iopt_dveg  =  ', Model%iopt_dveg
        print *,'iopt_crs   =  ', Model%iopt_crs
        print *,'iopt_btr   =  ', Model%iopt_btr
        print *,'iopt_run   =  ', Model%iopt_run
        print *,'iopt_sfc   =  ', Model%iopt_sfc
        print *,'iopt_frz   =  ', Model%iopt_frz
        print *,'iopt_inf   =  ', Model%iopt_inf
        print *,'iopt_rad   =  ', Model%iopt_rad
        print *,'iopt_alb   =  ', Model%iopt_alb
        print *,'iopt_snf   =  ', Model%iopt_snf
        print *,'iopt_tbot   =  ',Model%iopt_tbot
        print *,'iopt_stc   =  ', Model%iopt_stc
#ifdef CCPP
      elseif (Model%lsm == Model%lsm_ruc) then
        print *,' RUC Land Surface Model used'
#else
      elseif (Model%lsm == Model%lsm_ruc) then
        print *,' RUC Land Surface Model only available through CCPP - job aborted'
        stop
#endif
      else
        print *,' Unsupported LSM type - job aborted - lsm=',Model%lsm
        stop
      endif

      if (Model%lsm == Model%lsm_noahmp .and. Model%iopt_snf == 4) then
        if (Model%imp_physics /= Model%imp_physics_gfdl) stop 'iopt_snf == 4 must use GFDL MP'
      endif

      print *,' nst_anl=',Model%nst_anl,' use_ufo=',Model%use_ufo,' frac_grid=',Model%frac_grid
      print *,' min_lakeice=',Model%min_lakeice,' min_seaice=',Model%min_seaice
      if (Model%nstf_name(1) > 0 ) then
        print *,' NSSTM is active '
        print *,' nstf_name(1)=',Model%nstf_name(1)
        print *,' nstf_name(2)=',Model%nstf_name(2)
        print *,' nstf_name(3)=',Model%nstf_name(3)
        print *,' nstf_name(4)=',Model%nstf_name(4)
        print *,' nstf_name(5)=',Model%nstf_name(5)
      endif
      if (Model%do_deep) then
#ifdef CCPP
        ! Consistency check for NTDK convection: deep and shallow convection are bundled
        ! and cannot be combined with any other deep or shallow convection scheme
        if ( (Model%imfdeepcnv == Model%imfdeepcnv_ntiedtke .or. Model%imfshalcnv == Model%imfshalcnv_ntiedtke) .and. &
            .not. (Model%imfdeepcnv == Model%imfdeepcnv_ntiedtke .and. Model%imfshalcnv == Model%imfshalcnv_ntiedtke) ) then
            write(0,*) "Logic error: if NTDK deep convection is used, must also use NTDK shallow convection (and vice versa)"
            stop
        end if
#else
        if (Model%imfdeepcnv == 3 .or. Model%imfshalcnv == 3) then
            write(0,*) "Error, GF convection scheme only available through CCPP"
            stop
        else if (Model%imfdeepcnv == 4 .or. Model%imfshalcnv == 4) then
            write(0,*) "Error, NTDK convection scheme only available through CCPP"
            stop
        end if
#endif
        if (.not. Model%cscnv) then
          if (Model%ras) then
            print *,' RAS Convection scheme used with ccwf=',Model%ccwf
            Model%imfdeepcnv = -1
          else
            if (Model%imfdeepcnv == 0) then
               print *,' old SAS Convection scheme before July 2010 used'
#ifdef CCPP
            elseif(Model%imfdeepcnv == Model%imfdeepcnv_sas) then
               print *,' July 2010 version of SAS conv scheme used'
            elseif(Model%imfdeepcnv == Model%imfdeepcnv_samf) then
               print *,' scale & aerosol-aware mass-flux deep conv scheme'
            elseif(Model%imfdeepcnv == Model%imfdeepcnv_gf) then
               print *,' Grell-Freitas scale & aerosol-aware mass-flux deep conv scheme'
            elseif(Model%imfdeepcnv == Model%imfdeepcnv_ntiedtke) then
               print *,' New Tiedtke cumulus scheme'
#else
            elseif(Model%imfdeepcnv == 1) then
               print *,' July 2010 version of SAS conv scheme used'
            elseif(Model%imfdeepcnv == 2) then
               print *,' scale & aerosol-aware mass-flux deep conv scheme'
#endif
            endif
          endif
        else
          if (Model%do_aw) then
            print *,'Chikira-Sugiyama convection scheme with Arakawa-Wu'&
     &,                ' unified parameterization used'
          else
              print *,'Chikira-Sugiyama convection scheme used'
          endif
          print *,' cs_parm=',Model%cs_parm,' nctp=',Model%nctp
        endif
      else
        print*, ' Deep convection scheme disabled'
      endif
      if (Model%satmedmf) then
#ifdef CCPP
        if (Model%isatmedmf == Model%isatmedmf_vdif) then
          print *,' initial version (Nov 2018) of sale-aware TKE-based moist EDMF scheme used'
        elseif(Model%isatmedmf == Model%isatmedmf_vdifq) then
          print *,' update version (May 2019) of sale-aware TKE-based moist EDMF scheme used'
        endif
#else
        if (Model%isatmedmf == 0) then
          print *,' initial version (Nov 2018) of sale-aware TKE-based moist EDMF scheme used'
        elseif(Model%isatmedmf == 1) then
          print *,' update version (May 2019) of sale-aware TKE-based moist EDMF scheme used'
        endif
#endif
      elseif (Model%hybedmf) then
        print *,' scale-aware hybrid edmf PBL scheme used'
      elseif (Model%old_monin) then
        print *,' old (old_monin) PBL scheme used'
#ifdef CCPP
      elseif (Model%do_mynnedmf) then
        print *,' MYNN PBL scheme used'
      elseif (Model%do_myjpbl)then
        print *,' MYJ PBL scheme used'
#endif
      endif
      if (.not. Model%shal_cnv) then
        Model%imfshalcnv = -1
        print *,' No shallow convection used'
      else
        if (Model%imfshalcnv == 0) then
          print *,' modified Tiedtke eddy-diffusion shallow conv scheme used'
#ifdef CCPP
        elseif (Model%imfshalcnv == Model%imfshalcnv_sas) then
          print *,' July 2010 version of mass-flux shallow conv scheme used'
        elseif (Model%imfshalcnv == Model%imfshalcnv_samf) then
          print *,' scale- & aerosol-aware mass-flux shallow conv scheme (2017)'
        elseif (Model%imfshalcnv == Model%imfshalcnv_gf) then
          print *,' Grell-Freitas scale- & aerosol-aware mass-flux shallow conv scheme (2013)'
        elseif (Model%imfshalcnv == Model%imfshalcnv_ntiedtke) then
          print *,' New Tiedtke cumulus scheme'
#else
        elseif (Model%imfshalcnv == 1) then
          print *,' July 2010 version of mass-flux shallow conv scheme used'
        elseif (Model%imfshalcnv == 2) then
          print *,' scale- & aerosol-aware mass-flux shallow conv scheme (2017)'
#endif
        else
          print *,' unknown mass-flux scheme in use - defaulting to no shallow convection'
          Model%imfshalcnv = -1
        endif
      endif
      if (Model%do_gwd) then
        if (Model%do_ugwp) then
          print *,' Unified gravity wave drag parameterization used'
        else
          print *,' Original mountain blocking and oragraphic  gravity wave drag parameterization used'
          if (cdmbgwd(3) > 0.0) print *,' non-statioary gravity wave drag parameterization used'
        endif
          print *,' do_gwd=',Model%do_gwd
      endif
      if (Model%do_cnvgwd) then
        print *,' Convective GWD parameterization used, do_cnvgwd=',do_cnvgwd
      endif
      if (Model%crick_proof) print *,' CRICK-Proof cloud water used in radiation '
      if (Model%ccnorm)      print *,' Cloud condensate normalized by cloud cover for radiation'

      print *,' Radiative heating calculated at',Model%levr, ' layers'
      if (Model%iovr_sw == 0) then
        print *,' random cloud overlap for Shortwave IOVR_SW=',Model%iovr_sw
      else
        print *,' max-random cloud overlap for Shortwave IOVR_SW=',Model%iovr_sw
      endif
      if (Model%iovr_lw == 0) then
        print *,' random cloud overlap for Longwave IOVR_LW=',Model%iovr_lw
      else
        print *,' max-random cloud overlap for Longwave IOVR_LW=',Model%iovr_lw
      endif
      if (Model%isubc_sw == 0) then
        print *,' no sub-grid cloud for Shortwave ISUBC_SW=',Model%isubc_sw
      else
        print *,' sub-grid cloud for Shortwave ISUBC_SW=',Model%isubc_sw
      endif
      if (Model%isubc_lw == 0) then
        print *,' no sub-grid cloud for Longwave ISUBC_LW=',Model%isubc_lw
      else
        print *,' sub-grid cloud for Longwave ISUBC_LW=',Model%isubc_lw
      endif
    endif

!--- set up cloud schemes and tracer elements
    Model%nleffr = -999
    Model%nieffr = -999
    Model%nreffr = -999
    Model%nseffr = -999
    Model%ngeffr = -999
    if (Model%imp_physics == Model%imp_physics_zhao_carr) then
      Model%npdf3d  = 0
      Model%num_p3d = 4
      Model%num_p2d = 3
      Model%shcnvcw = .false.
      Model%ncnd    = 1                   ! ncnd is the number of cloud condensate types
      if (Model%me == Model%master) print *,' Using Zhao/Carr/Sundqvist Microphysics'

    elseif (Model%imp_physics == Model%imp_physics_zhao_carr_pdf) then !Zhao Microphysics with PDF cloud
      Model%npdf3d  = 3
      Model%num_p3d = 4
      Model%num_p2d = 3
      Model%ncnd    = 1
      if (Model%me == Model%master) print *,'Using Zhao/Carr/Sundqvist Microphysics with PDF Cloud'

    !else if (Model%imp_physics == 5) then        ! F-A goes here
    !  print *,' Ferrier Microphysics scheme has been deprecated - job aborted'
    !  stop
    else if (Model%imp_physics == Model%imp_physics_fer_hires) then     ! Ferrier-Aligo scheme
      Model%npdf3d  = 0
      Model%num_p3d = 3
      Model%num_p2d = 1
      Model%pdfcld  = .false.
      Model%shcnvcw = .false.
      Model%ncnd    = 5
      Model%nleffr = 1
      Model%nieffr = 2
      Model%nseffr = 3
      if (Model%me == Model%master) print *,' Using Ferrier-Aligo MP scheme', &
                                          ' microphysics', &
                                          ' lradar =',Model%lradar


    elseif (Model%imp_physics == Model%imp_physics_wsm6) then !WSM6 microphysics
      Model%npdf3d  = 0
      Model%num_p3d = 3
      Model%num_p2d = 1
      Model%pdfcld  = .false.
      Model%shcnvcw = .false.
      Model%ncnd    = 5
      Model%nleffr  = 1
      Model%nieffr  = 2
      Model%nseffr  = 3
      if (Model%me == Model%master) print *,' Using wsm6 microphysics'

    elseif (Model%imp_physics == Model%imp_physics_thompson) then !Thompson microphysics
      Model%npdf3d  = 0
      Model%num_p3d = 3
      Model%num_p2d = 1
      Model%pdfcld  = .false.
      Model%shcnvcw = .false.
      Model%ncnd    = 5
      Model%nleffr = 1
      Model%nieffr = 2
      Model%nseffr = 3
      if (Model%me == Model%master) print *,' Using Thompson double moment', &
                                          ' microphysics',' ltaerosol = ',Model%ltaerosol, &
                                          ' ttendlim =',Model%ttendlim, &
                                          ' lradar =',Model%lradar,Model%num_p3d,Model%num_p2d

    else if (Model%imp_physics == Model%imp_physics_mg) then        ! Morrison-Gettelman Microphysics
      Model%npdf3d  = 0
      Model%num_p3d = 5
      Model%num_p2d = 1
      Model%pdfcld  = .false.
      Model%shcnvcw = .false.
      Model%ncnd    = 2
      Model%nleffr  = 2
      Model%nieffr  = 3
      Model%nreffr  = 4
      Model%nseffr  = 5
      if (abs(Model%fprcp) == 1) then
        Model%ncnd  = 4
      elseif (Model%fprcp >= 2) then
        Model%ncnd  = 4
        if (Model%mg_do_graupel .or. Model%mg_do_hail) then
          Model%ncnd = 5
        endif
        Model%num_p3d = 6
        Model%ngeffr = 6
      endif
      if (Model%me == Model%master)                                                                 &
         print *,' Using Morrison-Gettelman double moment microphysics',                            &
                 ' aero_in=',         Model%aero_in,         ' iccn=',          Model%iccn,         &
                 ' mg_dcs=',          Model%mg_dcs,          ' mg_qcvar=',      Model%mg_qcvar,     &
                 ' mg_ts_auto_ice=',  Model%mg_ts_auto_ice,  ' pdfflag=',       Model%pdfflag,      &
                 ' mg_do_graupel=',   Model%mg_do_graupel,   ' mg_do_hail=',    Model%mg_do_hail,   &
                 ' mg_nccons=',       Model%mg_nccons,       ' mg_nicon=',      Model%mg_nicons,    &
                 ' mg_ngcons=',       Model%mg_ngcons ,      ' mg_ncnst=',      Model%mg_ncnst,     &
                 ' mg_ninst=',        Model%mg_ninst ,       ' mg_ngnst=',      Model%mg_ngnst,     &
                 ' sed_supersat=',    Model%sed_supersat ,   ' do_sb_physics=', Model%do_sb_physics,&
                 ' microp_uniform=',  Model%microp_uniform,  ' do_cldice=',     Model%do_cldice,    &
                 ' hetfrz_classnuc=', Model%hetfrz_classnuc, ' ncnd=',          Model%ncnd,         &
                 ' mg_alf=',          Model%mg_alf,          ' mg_qcmin=',      Model%mg_qcmin,     &
                 ' mg_do_ice_gmao=',  Model%mg_do_ice_gmao,  ' mg_do_liq_liu=', Model%mg_do_liq_liu

    elseif (Model%imp_physics == Model%imp_physics_gfdl) then !GFDL microphysics
      Model%npdf3d  = 0
      if(Model%effr_in) then
        Model%num_p3d = 5
        Model%nleffr = 1
        Model%nieffr = 2
        Model%nreffr = 3
        Model%nseffr = 4
        Model%ngeffr = 5
      else
        Model%num_p3d = 1
        ! Effective radii not used, point to valid index in dummy phy_f3d array
        Model%nleffr = 1
        Model%nieffr = 1
        Model%nreffr = 1
        Model%nseffr = 1
        Model%ngeffr = 1
      end if
      Model%num_p2d = 1
      Model%pdfcld  = .false.
      Model%shcnvcw = .false.
      Model%ncnd    = 5
      if (Model%me == Model%master) print *,' avg_max_length=',Model%avg_max_length
      if (Model%me == Model%master) print *,' Using GFDL Cloud Microphysics'
    else
      if (Model%me == Model%master) print *,'Wrong imp_physics value. Job abort.'
      stop
    endif

    if(Model%ras     .or. Model%cscnv)  Model%cnvcld = .false.
#ifdef CCPP
    if(Model%do_shoc .or. Model%pdfcld .or. Model%do_mynnedmf) Model%cnvcld = .false.
#else
    if(Model%do_shoc .or. Model%pdfcld) Model%cnvcld = .false.
#endif
    if(Model%cnvcld) Model%ncnvcld3d = 1

!--- get cnvw and cnvc indices in phy_f3d
    Model%ncnvw = -999
    Model%ncnvc = -999
    if ((Model%npdf3d == 3) .and. (Model%num_p3d == 4)) then
      Model%ncnvw = Model%num_p3d + 2
      Model%ncnvc = Model%ncnvw + 1
    elseif ((Model%npdf3d == 0) .and. (Model%ncnvcld3d == 1)) then
      Model%ncnvw = Model%num_p3d + 1
    endif
 
!--- derived totals for phy_f*d
    Model%ntot2d = Model%num_p2d + Model%nshoc_2d
    Model%ntot3d = Model%num_p3d + Model%nshoc_3d + Model%npdf3d + Model%ncnvcld3d
!
!   Unified cloud for SHOC and/or MG3
    Model%uni_cld = .false.
    Model%indcld  = -1
!   if (Model%shoc_cld .or. Model%ncld == 2 .or. Model%ntclamt > 0) then
    if (Model%imp_physics == Model%imp_physics_mg) then
      Model%uni_cld = .true.
      Model%indcld  = 1
    elseif (Model%shoc_cld) then
      Model%uni_cld = .true.
      Model%indcld  = Model%ntot3d - 2
    endif

#ifdef CCPP
    if (Model%do_shoc) then
      Model%nkbfshoc = Model%ntot3d   !< the index of upward kinematic buoyancy flux from SHOC in phy_f3d
      Model%nahdshoc = Model%ntot3d-1 !< the index of diffusivity for heat from from SHOC in phy_f3d
      Model%nscfshoc = Model%ntot3d-2 !< the index of subgrid-scale cloud fraction from from SHOC in phy_f3d
    else
      Model%nkbfshoc = -999
      Model%nahdshoc = -999
      Model%nscfshoc = -999
    endif
#endif

    if (me == Model%master)                                                     &
      write(0,*) ' num_p3d=',   Model%num_p3d,   ' num_p2d=',  Model%num_p2d,   &
                 ' crtrh=',     Model%crtrh,     ' npdf3d=',   Model%npdf3d,    &
                 ' pdfcld=',    Model%pdfcld,    ' shcnvcw=',  Model%shcnvcw,   &
                 ' cnvcld=',    Model%cnvcld,    ' ncnvcld3d=',Model%ncnvcld3d, &
                 ' do_shoc=',   Model%do_shoc,   ' nshoc3d=',  Model%nshoc_3d,  &
                 ' nshoc_2d=',  Model%nshoc_2d,  ' shoc_cld=', Model%shoc_cld,  &
#ifdef CCPP
                 ' nkbfshoc=',  Model%nkbfshoc,  ' nahdshoc=', Model%nahdshoc,  &
                 ' nscfshoc=',  Model%nscfshoc,                                 &
#endif
                 ' uni_cld=',   Model%uni_cld,                                  &
                 ' ntot3d=',    Model%ntot3d,    ' ntot2d=',   Model%ntot2d,    &
                 ' shocaftcnv=',Model%shocaftcnv,' indcld=',   Model%indcld,    &
                 ' shoc_parm=', Model%shoc_parm,                                &
                 ' ncnvw=',     Model%ncnvw,     ' ncnvc=',     Model%ncnvc

!--- END CODE FROM COMPNS_PHYSICS


!--- BEGIN CODE FROM GLOOPR
!--- set up parameters for Xu & Randell's cloudiness computation (Radiation)

    Model%lmfshal  = (Model%shal_cnv .and. Model%imfshalcnv > 0)
#ifdef CCPP
    Model%lmfdeep2 = (Model%imfdeepcnv == Model%imfdeepcnv_samf         &
                      .or. Model%imfdeepcnv == Model%imfdeepcnv_gf      &
                      .or. Model%imfdeepcnv == Model%imfdeepcnv_ntiedtke)
#else
    Model%lmfdeep2 = (Model%imfdeepcnv == 2)
#endif
!--- END CODE FROM GLOOPR

!--- BEGIN CODE FROM GLOOPB
!--- set up random number seed needed for RAS and old SAS and when cal_pre=.true.
!    Model%imfdeepcnv < 0 when Model%ras = .true.

    if (Model%imfdeepcnv <= 0 .or. Model%cal_pre ) then
      if (Model%random_clds) then
        seed0 = Model%idate(1) + Model%idate(2) + Model%idate(3) + Model%idate(4)
        call random_setseed(seed0)
        call random_number(wrk)
        Model%seed0 = seed0 + nint(wrk(1)*1000.0d0)
      endif
    endif
!--- END CODE FROM GLOOPB

    call Model%print ()

  end subroutine control_initialize


!------------------
! GFS_control%print
!------------------
  subroutine control_print(Model)

    implicit none

!--- interface variables
    class(GFS_control_type) :: Model
 
    if (Model%me == Model%master) then
      print *, ' '
      print *, 'basic control parameters'
      print *, ' me                : ', Model%me
      print *, ' master            : ', Model%master
#ifdef CCPP
      print *, ' communicator      : ', Model%communicator
#endif
      print *, ' nlunit            : ', Model%nlunit
      print *, ' fn_nml            : ', trim(Model%fn_nml)
      print *, ' fhzero            : ', Model%fhzero
      print *, ' ldiag3d           : ', Model%ldiag3d
      print *, ' lssav             : ', Model%lssav
      print *, ' fhcyc             : ', Model%fhcyc
      print *, ' thermodyn_id      : ', Model%thermodyn_id
      print *, ' sfcpress_id       : ', Model%sfcpress_id
      print *, ' gen_coord_hybrid  : ', Model%gen_coord_hybrid
      print *, ' '
      print *, 'grid extent parameters'
      print *, ' isc               : ', Model%isc
      print *, ' jsc               : ', Model%jsc
      print *, ' nx                : ', Model%nx
      print *, ' ny                : ', Model%ny
      print *, ' levs              : ', Model%levs
      print *, ' cnx               : ', Model%cnx
      print *, ' cny               : ', Model%cny
      print *, ' lonr              : ', Model%lonr
      print *, ' latr              : ', Model%latr
      print *, ' blksz(1)          : ', Model%blksz(1)
      print *, ' blksz(nblks)      : ', Model%blksz(Model%nblks)
      print *, ' '
      print *, 'coupling parameters'
      print *, ' cplflx            : ', Model%cplflx
      print *, ' cplwav            : ', Model%cplwav
      print *, ' cplchm            : ', Model%cplchm
      print *, ' '
      print *, 'integrated dynamics through earth atmosphere'
      print *, ' lsidea            : ', Model%lsidea
      print *, ' '
      print *, 'calendars and time parameters and activation triggers'
      print *, ' dtp               : ', Model%dtp
      print *, ' dtf               : ', Model%dtf
      print *, ' nscyc             : ', Model%nscyc
      print *, ' nszero            : ', Model%nszero
      print *, ' idat              : ', Model%idat
      print *, ' idate             : ', Model%idate
      print *, ' '
      print *, 'radiation control parameters'
      print *, ' fhswr             : ', Model%fhswr
      print *, ' fhlwr             : ', Model%fhlwr
      print *, ' nsswr             : ', Model%nsswr
      print *, ' nslwr             : ', Model%nslwr
#ifdef CCPP
      print *, ' nhfrad            : ', Model%nhfrad
#endif
      print *, ' levr              : ', Model%levr
      print *, ' nfxr              : ', Model%nfxr
      print *, ' aero_in           : ', Model%aero_in
#ifdef CCPP
      print *, ' ntrcaer           : ', Model%ntrcaer
#endif
      print *, ' lmfshal           : ', Model%lmfshal
      print *, ' lmfdeep2          : ', Model%lmfdeep2
      print *, ' nrcm              : ', Model%nrcm
      print *, ' iflip             : ', Model%iflip
      print *, ' isol              : ', Model%isol
      print *, ' ico2              : ', Model%ico2
      print *, ' ialb              : ', Model%ialb
      print *, ' iems              : ', Model%iems
      print *, ' iaer              : ', Model%iaer
      print *, ' icliq_sw          : ', Model%icliq_sw
      print *, ' iovr_sw           : ', Model%iovr_sw
      print *, ' iovr_lw           : ', Model%iovr_lw
      print *, ' ictm              : ', Model%ictm
      print *, ' isubc_sw          : ', Model%isubc_sw
      print *, ' isubc_lw          : ', Model%isubc_lw
      print *, ' crick_proof       : ', Model%crick_proof
      print *, ' ccnorm            : ', Model%ccnorm
      print *, ' norad_precip      : ', Model%norad_precip
      print *, ' lwhtr             : ', Model%lwhtr
      print *, ' swhtr             : ', Model%swhtr
      print *, ' '
      print *, 'microphysical switch'
      print *, ' ncld              : ', Model%ncld
      print *, ' imp_physics       : ', Model%imp_physics
      print *, ' '

      if (Model%imp_physics == Model%imp_physics_zhao_carr .or. Model%imp_physics == Model%imp_physics_zhao_carr_pdf) then
        print *, ' Z-C microphysical parameters'
        print *, ' psautco           : ', Model%psautco
        print *, ' prautco           : ', Model%prautco
        print *, ' evpco             : ', Model%evpco
        print *, ' wminco            : ', Model%wminco
        print *, ' '
      endif
      if (Model%imp_physics == Model%imp_physics_wsm6 .or. Model%imp_physics == Model%imp_physics_thompson) then
        print *, ' Thompson microphysical parameters'
        print *, ' ltaerosol         : ', Model%ltaerosol
        print *, ' lradar            : ', Model%lradar
        print *, ' lrefres           : ', Model%lrefres
        print *, ' ttendlim          : ', Model%ttendlim
        print *, ' '
      endif
      if (Model%imp_physics == Model%imp_physics_mg) then
        print *, ' M-G microphysical parameters'
        print *, ' fprcp             : ', Model%fprcp
        print *, ' mg_dcs            : ', Model%mg_dcs
        print *, ' mg_qcvar          : ', Model%mg_qcvar
        print *, ' mg_ts_auto_ice    : ', Model%mg_ts_auto_ice
        print *, ' mg_alf            : ', Model%mg_alf
        print *, ' mg_qcmin          : ', Model%mg_qcmin
        print *, ' pdfflag           : ', Model%pdfflag
        print *, ' '
      endif
      if (Model%imp_physics == Model%imp_physics_gfdl) then
        print *, ' GFDL microphysical parameters'
        print *, ' GFDL MP radiation inter: ', Model%lgfdlmprad
        print *, ' lrefres                : ', Model%lrefres
        print *, ' '
      endif
#ifdef CCPP
      if (Model%imp_physics == Model%imp_physics_fer_hires) then
        print *, ' Ferrier-Aligo microphysical parameters'
        print *, ' spec_adv          : ', Model%spec_adv
        print *, ' rhgrd             : ', Model%rhgrd
        print *, ' '
      endif
#endif
      print *, 'land/surface model parameters'
      print *, ' lsm               : ', Model%lsm
      print *, ' lsoil             : ', Model%lsoil
#ifdef CCPP
      print *, ' rdlai             : ', Model%rdlai
      print *, ' lsoil_lsm         : ', Model%lsoil_lsm
      print *, ' lsnow_lsm         : ', Model%lsnow_lsm
#endif
      print *, ' ivegsrc           : ', Model%ivegsrc
      print *, ' isot              : ', Model%isot

      if (Model%lsm == Model%lsm_noahmp) then
      print *, ' Noah MP LSM is used, the options are'
      print *, ' iopt_dveg         : ', Model%iopt_dveg
      print *, ' iopt_crs          : ', Model%iopt_crs
      print *, ' iopt_btr          : ', Model%iopt_btr
      print *, ' iopt_run          : ', Model%iopt_run
      print *, ' iopt_sfc          : ', Model%iopt_sfc
      print *, ' iopt_frz          : ', Model%iopt_frz
      print *, ' iopt_inf          : ', Model%iopt_inf
      print *, ' iopt_rad          : ', Model%iopt_rad
      print *, ' iopt_alb          : ', Model%iopt_alb
      print *, ' iopt_snf          : ', Model%iopt_snf
      print *, ' iopt_tbot         : ', Model%iopt_tbot
      print *, ' iopt_stc          : ', Model%iopt_stc

     endif

      print *, ' use_ufo           : ', Model%use_ufo
      print *, ' '
      print *, 'tuning parameters for physical parameterizations'
      print *, ' ras               : ', Model%ras
      if (Model%ras) then
        print *, ' psauras           : ', Model%psauras
        print *, ' prauras           : ', Model%prauras
        print *, ' wminras           : ', Model%wminras
      endif
      print *, ' flipv             : ', Model%flipv
      print *, ' trans_trac        : ', Model%trans_trac
      print *, ' old_monin         : ', Model%old_monin
      print *, ' do_gwd            : ', Model%do_gwd
      print *, ' cnvgwd            : ', Model%cnvgwd
      print *, ' do_cnvgwd         : ', Model%do_cnvgwd
      print *, ' mstrat            : ', Model%mstrat
      print *, ' moist_adj         : ', Model%moist_adj
      print *, ' cscnv             : ', Model%cscnv
      print *, ' cal_pre           : ', Model%cal_pre
      print *, ' do_aw             : ', Model%do_aw
      print *, ' flx_form          : ', Model%flx_form
      print *, ' do_shoc           : ', Model%do_shoc
      print *, ' shoc_parm         : ', Model%shoc_parm
      print *, ' shocaftcnv        : ', Model%shocaftcnv
      print *, ' shoc_cld          : ', Model%shoc_cld
      print *, ' uni_cld           : ', Model%uni_cld
      print *, ' h2o_phys          : ', Model%h2o_phys
      print *, ' pdfcld            : ', Model%pdfcld
      print *, ' shcnvcw           : ', Model%shcnvcw
      print *, ' redrag            : ', Model%redrag
      print *, ' hybedmf           : ', Model%hybedmf
      print *, ' satmedmf          : ', Model%satmedmf
      print *, ' isatmedmf         : ', Model%isatmedmf
      print *, ' shinhong          : ', Model%shinhong
      print *, ' do_ysu            : ', Model%do_ysu
      print *, ' dspheat           : ', Model%dspheat
      print *, ' lheatstrg         : ', Model%lheatstrg
      print *, ' cnvcld            : ', Model%cnvcld
      print *, ' random_clds       : ', Model%random_clds
      print *, ' shal_cnv          : ', Model%shal_cnv
      print *, ' imfshalcnv        : ', Model%imfshalcnv
      print *, ' imfdeepcnv        : ', Model%imfdeepcnv
      print *, ' do_deep           : ', Model%do_deep
      print *, ' nmtvr             : ', Model%nmtvr
      print *, ' jcap              : ', Model%jcap
      print *, ' cs_parm           : ', Model%cs_parm
      print *, ' flgmin            : ', Model%flgmin
      print *, ' cgwf              : ', Model%cgwf
      print *, ' ccwf              : ', Model%ccwf
      print *, ' cdmbgwd           : ', Model%cdmbgwd
      print *, ' sup               : ', Model%sup
      print *, ' ctei_rm           : ', Model%ctei_rm
      print *, ' crtrh             : ', Model%crtrh
      print *, ' dlqf              : ', Model%dlqf
      print *, ' seed0             : ', Model%seed0
      print *, ' rbcr              : ', Model%rbcr
#ifdef CCPP
      print *, ' do_mynnedmf       : ', Model%do_mynnedmf
      print *, ' do_mynnsfclay     : ', Model%do_mynnsfclay
      print *, ' do_myjsfc         : ', Model%do_myjsfc
      print *, ' do_myjpbl         : ', Model%do_myjpbl
      print *, ' gwd_opt           : ', Model%gwd_opt
#endif
      print *, ' '
      print *, 'Rayleigh friction'
      print *, ' prslrd0           : ', Model%prslrd0
      print *, ' ral_ts            : ', Model%ral_ts
      print *, ' '
      if (Model%imfdeepcnv >= 0) then
        print *, 'mass flux deep convection'
        print *, ' clam_deep         : ', Model%clam_deep
        print *, ' c0s_deep          : ', Model%c0s_deep
        print *, ' c1_deep           : ', Model%c1_deep
        print *, ' betal_deep        : ', Model%betal_deep
        print *, ' betas_deep        : ', Model%betas_deep
        print *, ' evfact_deep       : ', Model%evfact_deep
        print *, ' evfactl_deep      : ', Model%evfactl_deep
        print *, ' pgcon_deep        : ', Model%pgcon_deep
        print *, ' asolfac_deep      : ', Model%asolfac_deep
        print *, ' '
      endif
      if (Model%imfshalcnv >= 0) then
        print *, 'mass flux shallow convection'
        print *, ' clam_shal         : ', Model%clam_shal
        print *, ' c0s_shal          : ', Model%c0s_shal
        print *, ' c1_shal           : ', Model%c1_shal
        print *, ' pgcon_shal        : ', Model%pgcon_shal
        print *, ' asolfac_shal      : ', Model%asolfac_shal
      endif
      print *, ' '
      print *, 'near surface sea temperature model'
      print *, ' nst_anl           : ', Model%nst_anl
      print *, ' nstf_name         : ', Model%nstf_name
      print *, ' lsea              : ', Model%lsea
      print *, ' '
      print *, 'surface layer options'
      print *, ' sfc_z0_type       : ', Model%sfc_z0_type
      print *, ' '
      print *, 'background vertical diffusion coefficients'
      print *, ' xkzm_m            : ', Model%xkzm_m
      print *, ' xkzm_h            : ', Model%xkzm_h
      print *, ' xkzm_s            : ', Model%xkzm_s
      print *, ' xkzminv           : ', Model%xkzminv
      print *, ' moninq_fac        : ', Model%moninq_fac
      print *, ' dspfac            : ', Model%dspfac
      print *, ' bl_upfr           : ', Model%bl_upfr
      print *, ' bl_dnfr           : ', Model%bl_dnfr
      print *, ' '
      print *, 'stochastic physics'
      print *, ' do_sppt           : ', Model%do_sppt
      print *, ' do_shum           : ', Model%do_shum
      print *, ' do_skeb           : ', Model%do_skeb
      print *, ' do_sfcperts       : ', Model%do_sfcperts
      print *, ' '
      print *, 'cellular automata'
      print *, ' nca               : ', Model%ncells
      print *, ' ncells            : ', Model%ncells
      print *, ' nlives            : ', Model%nlives
      print *, ' nfracseed         : ', Model%nfracseed
      print *, ' nseed             : ', Model%nseed
      print *, ' ca_global         : ', Model%ca_global
      print *, ' ca_sgs            : ', Model%ca_sgs
      print *, ' do_ca             : ', Model%do_ca
      print *, ' iseed_ca          : ', Model%iseed_ca
      print *, ' ca_smooth         : ', Model%ca_smooth
      print *, ' isppt_deep        : ', Model%isppt_deep
      print *, ' nspinup           : ', Model%nspinup
      print *, ' nthresh           : ', Model%nthresh
      print *, ' '
      print *, 'tracers'
      print *, ' tracer_names      : ', Model%tracer_names
      print *, ' ntrac             : ', Model%ntrac
#ifdef CCPP
      print *, ' ntqv              : ', Model%ntqv
      print *, ' nqrimef           : ', Model%nqrimef
#endif
      print *, ' ntoz              : ', Model%ntoz
      print *, ' ntcw              : ', Model%ntcw
      print *, ' ntiw              : ', Model%ntiw
      print *, ' ntrw              : ', Model%ntrw
      print *, ' ntsw              : ', Model%ntsw
      print *, ' ntgl              : ', Model%ntgl
      print *, ' ntclamt           : ', Model%ntclamt
      print *, ' ntlnc             : ', Model%ntlnc
      print *, ' ntinc             : ', Model%ntinc
      print *, ' ntrnc             : ', Model%ntrnc
      print *, ' ntsnc             : ', Model%ntsnc
      print *, ' ntgnc             : ', Model%ntgnc
      print *, ' ntke              : ', Model%ntke
      print *, ' nto               : ', Model%nto
      print *, ' nto2              : ', Model%nto2
      print *, ' ntwa              : ', Model%ntwa
      print *, ' ntia              : ', Model%ntia
      print *, ' ntchm             : ', Model%ntchm
      print *, ' ntchs             : ', Model%ntchs
      print *, ' fscav             : ', Model%fscav
      print *, ' '
      print *, 'derived totals for phy_f*d'
      print *, ' ntot2d            : ', Model%ntot2d
      print *, ' ntot3d            : ', Model%ntot3d
      print *, ' num_p2d           : ', Model%num_p2d
      print *, ' num_p3d           : ', Model%num_p3d
      print *, ' nshoc_2d          : ', Model%nshoc_2d
      print *, ' nshoc_3d          : ', Model%nshoc_3d
      print *, ' ncnvcld3d         : ', Model%ncnvcld3d
      print *, ' npdf3d            : ', Model%npdf3d
      print *, ' nctp              : ', Model%nctp
#ifdef CCPP
      print *, ' nkbfshoc          : ', Model%nkbfshoc
      print *, ' nahdshoc          : ', Model%nahdshoc
      print *, ' nscfshoc          : ', Model%nscfshoc
#endif
      print *, ' '
      print *, 'debug flags'
      print *, ' debug             : ', Model%debug 
      print *, ' pre_rad           : ', Model%pre_rad
      print *, ' '
      print *, 'variables modified at each time step'
      print *, ' ipt               : ', Model%ipt
      print *, ' lprnt             : ', Model%lprnt
      print *, ' lsswr             : ', Model%lsswr
      print *, ' lslwr             : ', Model%lslwr
      print *, ' solhr             : ', Model%solhr
      print *, ' solcon            : ', Model%solcon
      print *, ' slag              : ', Model%slag
      print *, ' sdec              : ', Model%sdec
      print *, ' cdec              : ', Model%cdec
      print *, ' clstp             : ', Model%clstp
      print *, ' phour             : ', Model%phour
      print *, ' fhour             : ', Model%fhour
      print *, ' zhour             : ', Model%zhour
      print *, ' kdt               : ', Model%kdt
      print *, ' jdat              : ', Model%jdat
#ifdef CCPP
      print *, ' sec               : ', Model%sec
      print *, ' si                : ', Model%si
      print *, ' first_time_step   : ', Model%first_time_step
      print *, ' restart           : ', Model%restart
      print *, ' hydrostatic       : ', Model%hydrostatic
#endif
    endif

  end subroutine control_print


!----------------
! GFS_grid%create
!----------------
  subroutine grid_create (Grid, IM, Model)

    implicit none

    class(GFS_grid_type)              :: Grid
    integer,                intent(in) :: IM
    type(GFS_control_type), intent(in) :: Model

    allocate (Grid%xlon   (IM))
    allocate (Grid%xlat   (IM))
    allocate (Grid%xlat_d (IM))
    allocate (Grid%xlon_d (IM))
    allocate (Grid%sinlat (IM))
    allocate (Grid%coslat (IM))
    allocate (Grid%area   (IM))
    allocate (Grid%dx     (IM))

    Grid%xlon   = clear_val
    Grid%xlat   = clear_val
    Grid%xlat_d = clear_val
    Grid%xlon_d = clear_val
    Grid%sinlat = clear_val
    Grid%coslat = clear_val
    Grid%area   = clear_val
    Grid%dx     = clear_val

!--- ozone active
    if ( Model%ntoz > 0 ) then
      allocate (Grid%ddy_o3    (IM))
      allocate (Grid%jindx1_o3 (IM))
      allocate (Grid%jindx2_o3 (IM))
    endif

!--- stratosphere h2o active
    if ( Model%h2o_phys ) then
      allocate (Grid%ddy_h    (IM))
      allocate (Grid%jindx1_h (IM))
      allocate (Grid%jindx2_h (IM))
    endif

!--- iccn active
    if ( Model%iccn ) then
      allocate (Grid%ddy_ci    (IM))
      allocate (Grid%jindx1_ci (IM))
      allocate (Grid%jindx2_ci (IM))
      allocate (Grid%ddx_ci    (IM))
      allocate (Grid%iindx1_ci (IM))
      allocate (Grid%iindx2_ci (IM))
    endif

!--- iaerclm active
    if ( Model%aero_in ) then
      allocate (Grid%ddy_aer   (IM))
      allocate (Grid%jindx1_aer(IM))
      allocate (Grid%jindx2_aer(IM))
      allocate (Grid%ddx_aer   (IM))
      allocate (Grid%iindx1_aer(IM))
      allocate (Grid%iindx2_aer(IM))
    endif
 end subroutine grid_create


!--------------------
! GFS_tbd_type%create
!--------------------
#ifndef CCPP
  subroutine tbd_create (Tbd, IM, BLKNO, Model)
#else
  subroutine tbd_create (Tbd, IM, Model)
#endif

    implicit none

    class(GFS_tbd_type)                :: Tbd
    integer,                intent(in) :: IM
#ifndef CCPP
    integer,                intent(in) :: BLKNO
#endif
    type(GFS_control_type), intent(in) :: Model

!--- In
!--- sub-grid cloud radiation
    if ( Model%isubc_lw == 2 .or. Model%isubc_sw == 2 ) then
      allocate (Tbd%icsdsw (IM))
      allocate (Tbd%icsdlw (IM))
    endif

!--- ozone and stratosphere h2o needs
    ! DH* oz_coeff is set to zero if both ozphys options are false,
    ! better to use conditional allocations here for ozpl (and h2opl)? *DH
    allocate (Tbd%ozpl  (IM,levozp,oz_coeff))
    allocate (Tbd%h2opl (IM,levh2o,h2o_coeff))
    Tbd%ozpl  = clear_val
    Tbd%h2opl = clear_val

!--- ccn and in needs
    ! DH* allocate only for MG? *DH
    allocate (Tbd%in_nm  (IM,Model%levs))
    allocate (Tbd%ccn_nm (IM,Model%levs))
    Tbd%in_nm  = clear_val
    Tbd%ccn_nm = clear_val

!--- aerosol fields
    ! DH* allocate only for MG? *DH
    allocate (Tbd%aer_nm  (IM,Model%levs,ntrcaer))
    Tbd%aer_nm = clear_val

#ifdef CCPP
! DH* TODO - MOVE THIS TO a block-vector dependent structure in GFS_control?
! e.g. GFS_Control%imap(blk), GFS_Control%jmap(blk), or ii instead if imap etc? *DH
!--- maps of local index ix to global indices i and j for this block
    allocate (Tbd%imap (IM))
    allocate (Tbd%jmap (IM))
    Tbd%imap = 0
    Tbd%jmap = 0
#endif

    allocate (Tbd%rann (IM,Model%nrcm))
    Tbd%rann = rann_init

!--- In/Out
    allocate (Tbd%acv  (IM))
    allocate (Tbd%acvb (IM))
    allocate (Tbd%acvt (IM))

    Tbd%acv  = clear_val
    Tbd%acvb = clear_val
    Tbd%acvt = clear_val

    if (Model%do_sppt) then
      allocate (Tbd%dtdtr     (IM,Model%levs))
      allocate (Tbd%dtotprcp  (IM))
      allocate (Tbd%dcnvprcp  (IM))
      allocate (Tbd%drain_cpl (IM))
      allocate (Tbd%dsnow_cpl (IM))

      Tbd%dtdtr     = clear_val
      Tbd%dtotprcp  = clear_val
      Tbd%dcnvprcp  = clear_val
      Tbd%drain_cpl = clear_val
      Tbd%dsnow_cpl = clear_val
    endif

    allocate (Tbd%phy_f2d  (IM,Model%ntot2d))
    allocate (Tbd%phy_f3d  (IM,Model%levs,Model%ntot3d))
    Tbd%phy_f2d  = clear_val
    Tbd%phy_f3d  = clear_val

    if (Model%nctp > 0 .and. Model%cscnv) then
      allocate (Tbd%phy_fctd (IM,Model%nctp))
      Tbd%phy_fctd = clear_val
    endif

!   if (Model%do_shoc) Tbd%phy_f3d(:,1,Model%ntot3d-1) = 3.0
!   if (Model%do_shoc) Tbd%phy_f3d(:,:,Model%ntot3d-1) = 1.0

    allocate (Tbd%hpbl (IM))
    Tbd%hpbl     = clear_val

#ifndef CCPP
    Tbd%blkno = BLKNO
#endif

#ifdef CCPP
    allocate (Tbd%htlwc (IM,Model%levr+LTP))
    allocate (Tbd%htlw0 (IM,Model%levr+LTP))
    allocate (Tbd%htswc (IM,Model%levr+LTP))
    allocate (Tbd%htsw0 (IM,Model%levr+LTP))

    Tbd%htlwc = clear_val
    Tbd%htlw0 = clear_val
    Tbd%htswc = clear_val
    Tbd%htsw0 = clear_val

    if (Model%imfdeepcnv == Model%imfdeepcnv_gf .or. Model%imfdeepcnv == Model%imfdeepcnv_ntiedtke) then
       allocate(Tbd%forcet(IM, Model%levs))
       allocate(Tbd%forceq(IM, Model%levs))
       allocate(Tbd%prevst(IM, Model%levs))
       allocate(Tbd%prevsq(IM, Model%levs))
       Tbd%forcet = clear_val
       Tbd%forceq = clear_val
       Tbd%prevst = clear_val
       Tbd%prevsq = clear_val
   end if

   if (Model%imfdeepcnv == Model%imfdeepcnv_gf) then
      allocate(Tbd%cactiv(IM))
      Tbd%cactiv = zero
   end if

    !--- MYNN variables:
    if (Model%do_mynnedmf) then
       !print*,"Allocating all MYNN-EDMF variables:"
       allocate (Tbd%cldfra_bl (IM,Model%levs))
       allocate (Tbd%qc_bl     (IM,Model%levs))
       allocate (Tbd%el_pbl    (IM,Model%levs))
       allocate (Tbd%sh3d      (IM,Model%levs))
       allocate (Tbd%qke       (IM,Model%levs))
       allocate (Tbd%tsq       (IM,Model%levs))
       allocate (Tbd%qsq       (IM,Model%levs))
       allocate (Tbd%cov       (IM,Model%levs))
       !print*,"Allocating all MYNN-EDMF variables:"
       Tbd%cldfra_bl     = clear_val
       Tbd%qc_bl         = clear_val
       Tbd%el_pbl        = clear_val
       Tbd%sh3d          = clear_val
       Tbd%qke           = zero
       Tbd%tsq           = clear_val
       Tbd%qsq           = clear_val
       Tbd%cov           = clear_val
    end if

    ! MYJ variables
    if (Model%do_myjsfc.or.Model%do_myjpbl) then
       !print*,"Allocating all MYJ surface variables:"
       allocate (Tbd%phy_myj_qsfc   (IM))
       allocate (Tbd%phy_myj_thz0   (IM)) 
       allocate (Tbd%phy_myj_qz0    (IM)) 
       allocate (Tbd%phy_myj_uz0    (IM)) 
       allocate (Tbd%phy_myj_vz0    (IM)) 
       allocate (Tbd%phy_myj_z0base (IM)) 
       allocate (Tbd%phy_myj_akhs   (IM)) 
       allocate (Tbd%phy_myj_akms   (IM)) 
       allocate (Tbd%phy_myj_chkqlm (IM)) 
       allocate (Tbd%phy_myj_elflx  (IM)) 
       allocate (Tbd%phy_myj_a1u    (IM)) 
       allocate (Tbd%phy_myj_a1t    (IM)) 
       allocate (Tbd%phy_myj_a1q    (IM))
       !print*,"Allocating all MYJ schemes variables:"
       Tbd%phy_myj_qsfc   = clear_val 
       Tbd%phy_myj_thz0   = clear_val 
       Tbd%phy_myj_qz0    = clear_val 
       Tbd%phy_myj_uz0    = clear_val 
       Tbd%phy_myj_vz0    = clear_val 
       Tbd%phy_myj_z0base = clear_val 
       Tbd%phy_myj_akhs   = clear_val 
       Tbd%phy_myj_akms   = clear_val 
       Tbd%phy_myj_chkqlm = clear_val 
       Tbd%phy_myj_elflx  = clear_val 
       Tbd%phy_myj_a1u    = clear_val 
       Tbd%phy_myj_a1t    = clear_val 
       Tbd%phy_myj_a1q    = clear_val 
    end if
#endif

  end subroutine tbd_create


!------------------------
! GFS_cldprop_type%create
!------------------------
  subroutine cldprop_create (Cldprop, IM, Model)

    implicit none

    class(GFS_cldprop_type)            :: Cldprop
    integer,                intent(in) :: IM
    type(GFS_control_type), intent(in) :: Model

    allocate (Cldprop%cv  (IM))
    allocate (Cldprop%cvt (IM)) 
    allocate (Cldprop%cvb (IM))
    
    Cldprop%cv  = clear_val
    Cldprop%cvt = clear_val
    Cldprop%cvb = clear_val
  
  end subroutine cldprop_create


!******************************************
! GFS_radtend_type%create
!******************************************
  subroutine radtend_create (Radtend, IM, Model)
                               
    implicit none
       
    class(GFS_radtend_type)            :: Radtend
    integer,                intent(in) :: IM
    type(GFS_control_type), intent(in) :: Model

    !--- Out (radiation only) 
    allocate (Radtend%sfcfsw (IM))
    allocate (Radtend%sfcflw (IM))

    Radtend%sfcfsw%upfxc = clear_val
    Radtend%sfcfsw%upfx0 = clear_val
    Radtend%sfcfsw%dnfxc = clear_val
    Radtend%sfcfsw%dnfx0 = clear_val
    Radtend%sfcflw%upfxc = clear_val
    Radtend%sfcflw%upfx0 = clear_val
    Radtend%sfcflw%dnfxc = clear_val
    Radtend%sfcflw%dnfx0 = clear_val
         
    allocate (Radtend%htrsw  (IM,Model%levs))
    allocate (Radtend%htrlw  (IM,Model%levs))
    allocate (Radtend%sfalb  (IM))
    allocate (Radtend%coszen (IM))
    allocate (Radtend%tsflw  (IM))
    allocate (Radtend%semis  (IM))

    Radtend%htrsw  = clear_val
    Radtend%htrlw  = clear_val
    Radtend%sfalb  = clear_val
    Radtend%coszen = clear_val
    Radtend%tsflw  = clear_val
    Radtend%semis  = clear_val
             
!--- In/Out (???) (radiation only)
    allocate (Radtend%coszdg (IM))

    Radtend%coszdg = clear_val
             
!--- In/Out (???) (physics only)
    allocate (Radtend%swhc  (IM,Model%levs))
    allocate (Radtend%lwhc  (IM,Model%levs))
    allocate (Radtend%lwhd  (IM,Model%levs,6))

    Radtend%lwhd  = clear_val
    Radtend%lwhc  = clear_val
    Radtend%swhc  = clear_val

  end subroutine radtend_create


!----------------
! GFS_diag%create
!----------------
  subroutine diag_create (Diag, IM, Model)
    class(GFS_diag_type)               :: Diag
    integer,                intent(in) :: IM
    type(GFS_control_type), intent(in) :: Model

!
    logical, save :: linit

    !--- Radiation
    allocate (Diag%fluxr   (IM,Model%nfxr))
    allocate (Diag%topfsw  (IM))
    allocate (Diag%topflw  (IM))
!--- Physics
!--- In/Out
    allocate (Diag%srunoff (IM))
    allocate (Diag%evbsa   (IM))
    allocate (Diag%evcwa   (IM))
    allocate (Diag%snohfa  (IM))
    allocate (Diag%transa  (IM))
    allocate (Diag%sbsnoa  (IM))
    allocate (Diag%snowca  (IM))
    allocate (Diag%soilm   (IM))
    allocate (Diag%tmpmin  (IM))
    allocate (Diag%tmpmax  (IM))
    allocate (Diag%dusfc   (IM))
    allocate (Diag%dvsfc   (IM))
    allocate (Diag%dtsfc   (IM))
    allocate (Diag%dqsfc   (IM))
    allocate (Diag%totprcp (IM))
    allocate (Diag%totprcpb(IM))
    allocate (Diag%gflux   (IM))
    allocate (Diag%dlwsfc  (IM))
    allocate (Diag%ulwsfc  (IM))
    allocate (Diag%suntim  (IM))
    allocate (Diag%runoff  (IM))
    allocate (Diag%ep      (IM))
    allocate (Diag%cldwrk  (IM))
    allocate (Diag%dugwd   (IM))
    allocate (Diag%dvgwd   (IM))
    allocate (Diag%psmean  (IM))
    allocate (Diag%cnvprcp (IM))
    allocate (Diag%cnvprcpb(IM))
    allocate (Diag%spfhmin (IM))
    allocate (Diag%spfhmax (IM))
    allocate (Diag%u10mmax (IM))
    allocate (Diag%v10mmax (IM))
    allocate (Diag%wind10mmax (IM))
    allocate (Diag%u10max (IM))
    allocate (Diag%v10max (IM))
    allocate (Diag%spd10max (IM))
    allocate (Diag%rain    (IM))
    allocate (Diag%rainc   (IM))
    allocate (Diag%ice     (IM))
    allocate (Diag%snow    (IM))
    allocate (Diag%graupel (IM))
    allocate (Diag%totice  (IM))
    allocate (Diag%totsnw  (IM))
    allocate (Diag%totgrp  (IM))
    allocate (Diag%toticeb (IM))
    allocate (Diag%totsnwb (IM))
    allocate (Diag%totgrpb (IM))
    allocate (Diag%u10m    (IM))
    allocate (Diag%v10m    (IM))
    allocate (Diag%dpt2m   (IM))
    allocate (Diag%zlvl    (IM))
    allocate (Diag%psurf   (IM))
    allocate (Diag%pwat    (IM))
    allocate (Diag%t1      (IM))
    allocate (Diag%q1      (IM))
    allocate (Diag%u1      (IM))
    allocate (Diag%v1      (IM))
    allocate (Diag%chh     (IM))
    allocate (Diag%cmm     (IM))
    allocate (Diag%dlwsfci (IM))
    allocate (Diag%ulwsfci (IM))
    allocate (Diag%dswsfci (IM))
#ifdef CCPP
    allocate (Diag%nswsfci (IM))
#endif
    allocate (Diag%uswsfci (IM))
    allocate (Diag%dusfci  (IM))
    allocate (Diag%dvsfci  (IM))
    allocate (Diag%dtsfci  (IM))
    allocate (Diag%dqsfci  (IM))
    allocate (Diag%gfluxi  (IM))
    allocate (Diag%epi     (IM))
    allocate (Diag%smcwlt2 (IM))
    allocate (Diag%smcref2 (IM))
    if (.not. Model%lsm == Model%lsm_ruc) then
      allocate (Diag%wet1    (IM))
    end if
    allocate (Diag%sr      (IM))
    allocate (Diag%tdomr   (IM))
    allocate (Diag%tdomzr  (IM))
    allocate (Diag%tdomip  (IM))
    allocate (Diag%tdoms   (IM))
    allocate (Diag%skebu_wts(IM,Model%levs))
    allocate (Diag%skebv_wts(IM,Model%levs))
    allocate (Diag%sppt_wts(IM,Model%levs))
    allocate (Diag%shum_wts(IM,Model%levs))
    allocate (Diag%zmtnblck(IM))    

    ! F-A MP scheme
#ifdef CCPP
    if (Model%imp_physics == Model%imp_physics_fer_hires) then
     allocate (Diag%TRAIN     (IM,Model%levs))
    end if
#endif

    allocate (Diag%ca_out  (IM))
    allocate (Diag%ca_deep  (IM))
    allocate (Diag%ca_turb  (IM))
    allocate (Diag%ca_shal  (IM))
    allocate (Diag%ca_rad (IM))
    allocate (Diag%ca_micro  (IM))
    
    !--- 3D diagnostics
    if (Model%ldiag3d) then
      allocate (Diag%du3dt  (IM,Model%levs,4))
      allocate (Diag%dv3dt  (IM,Model%levs,4))
      allocate (Diag%dt3dt  (IM,Model%levs,7))
      allocate (Diag%dq3dt  (IM,Model%levs,9))
!      allocate (Diag%dq3dt  (IM,Model%levs,oz_coeff+5))
!--- needed to allocate GoCart coupling fields
!      allocate (Diag%upd_mf (IM,Model%levs))
!      allocate (Diag%dwn_mf (IM,Model%levs))
!      allocate (Diag%det_mf (IM,Model%levs))
!      allocate (Diag%cldcov (IM,Model%levs))
    endif

!vay-2018
    if (Model%ldiag_ugwp) then
      allocate (Diag%du3dt_dyn  (IM,Model%levs) )

      allocate (Diag%du3dt_pbl  (IM,Model%levs) )
      allocate (Diag%dv3dt_pbl  (IM,Model%levs) )
      allocate (Diag%dt3dt_pbl  (IM,Model%levs) )

      allocate (Diag%du3dt_ogw  (IM,Model%levs) )
      allocate (Diag%dv3dt_ogw  (IM,Model%levs) )
      allocate (Diag%dt3dt_ogw  (IM,Model%levs) )

      allocate (Diag%du3dt_mtb  (IM,Model%levs) )
      allocate (Diag%dv3dt_mtb  (IM,Model%levs) )
      allocate (Diag%dt3dt_mtb  (IM,Model%levs) )

      allocate (Diag%du3dt_tms  (IM,Model%levs) )
      allocate (Diag%dv3dt_tms  (IM,Model%levs) )
      allocate (Diag%dt3dt_tms  (IM,Model%levs) )

      allocate (Diag%du3dt_ngw  (IM,Model%levs) )
      allocate (Diag%dv3dt_ngw  (IM,Model%levs) )
      allocate (Diag%dt3dt_ngw  (IM,Model%levs) )

      allocate (Diag%du3dt_cgw  (IM,Model%levs) )
      allocate (Diag%dv3dt_cgw  (IM,Model%levs) )
      allocate (Diag%dt3dt_moist  (IM,Model%levs) )

      allocate (Diag%dudt_tot  (IM,Model%levs) )
      allocate (Diag%dvdt_tot  (IM,Model%levs) )
      allocate (Diag%dtdt_tot  (IM,Model%levs) )

       allocate (Diag%uav_ugwp  (IM,Model%levs) )
       allocate (Diag%tav_ugwp  (IM,Model%levs) )
    endif

       allocate (Diag%zmtb      (IM) )
       allocate (Diag%zogw      (IM) )
       allocate (Diag%zlwb      (IM) )
       allocate (Diag%tau_ogw   (IM) )
       allocate (Diag%tau_ngw   (IM) )
       allocate (Diag%tau_mtb   (IM) )
       allocate (Diag%tau_tofd  (IM) )
!   endif

!
!ugwp - instant
!
    if (Model%do_ugwp) then
      allocate (Diag%gwp_ax  (IM,Model%levs) )
      allocate (Diag%gwp_ay  (IM,Model%levs) )
      allocate (Diag%gwp_dtdt(IM,Model%levs) )
      allocate (Diag%gwp_kdis(IM,Model%levs) )

      allocate (Diag%gwp_axo  (IM,Model%levs) )
      allocate (Diag%gwp_ayo  (IM,Model%levs) )
      allocate (Diag%gwp_axc  (IM,Model%levs) )
      allocate (Diag%gwp_ayc  (IM,Model%levs) )
      allocate (Diag%gwp_axf  (IM,Model%levs) )
      allocate (Diag%gwp_ayf  (IM,Model%levs) )
!GW-sources
      allocate (Diag%gwp_dcheat(IM,Model%levs) )
      allocate (Diag%gwp_scheat(IM,Model%levs) )
      allocate (Diag%gwp_fgf  (IM            ) )
      allocate (Diag%gwp_okw  (IM            ) )

      allocate (Diag%gwp_precip(IM) )
      allocate (Diag%gwp_klevs (IM, 3) )

    endif

    !--- 3D diagnostics for Thompson MP / GFDL MP
    allocate (Diag%refl_10cm(IM,Model%levs))

    !--  New max hourly diag.
    allocate (Diag%refdmax(IM))
    allocate (Diag%refdmax263k(IM))
    allocate (Diag%t02max(IM))
    allocate (Diag%t02min(IM))
    allocate (Diag%rh02max(IM))
    allocate (Diag%rh02min(IM))

#ifdef CCPP
    !--- MYNN variables:
    if (Model%do_mynnedmf) then
      !print*,"Allocating all MYNN-EDMF variables:"
      allocate (Diag%edmf_a    (IM,Model%levs))
      allocate (Diag%edmf_w    (IM,Model%levs))
      allocate (Diag%edmf_qt   (IM,Model%levs))
      allocate (Diag%edmf_thl  (IM,Model%levs))
      allocate (Diag%edmf_ent  (IM,Model%levs))
      allocate (Diag%edmf_qc   (IM,Model%levs))
      allocate (Diag%nupdraft  (IM))
      allocate (Diag%maxmf     (IM))
      allocate (Diag%ktop_shallow(IM))
      allocate (Diag%exch_h    (IM,Model%levs))
      allocate (Diag%exch_m    (IM,Model%levs))
      !print*,"Initializing all MYNN-EDMF variables with ",clear_val
      Diag%edmf_a        = clear_val
      Diag%edmf_w        = clear_val
      Diag%edmf_qt       = clear_val
      Diag%edmf_thl      = clear_val
      Diag%edmf_ent      = clear_val
      Diag%edmf_qc       = clear_val
      Diag%nupdraft      = 0
      Diag%maxmf         = clear_val
      Diag%ktop_shallow  = 0
      Diag%exch_h        = clear_val
      Diag%exch_m        = clear_val
    endif

    !--- Drag Suite variables:
    if (Model%gwd_opt == 33) then
      !print*,"Allocating all Drag Suite variables:"
      allocate (Diag%dtaux2d_ls  (IM,Model%levs))
      allocate (Diag%dtauy2d_ls  (IM,Model%levs))
      allocate (Diag%dtaux2d_bl  (IM,Model%levs))
      allocate (Diag%dtauy2d_bl  (IM,Model%levs))
      allocate (Diag%dtaux2d_ss  (IM,Model%levs))
      allocate (Diag%dtauy2d_ss  (IM,Model%levs))
      allocate (Diag%dtaux2d_fd  (IM,Model%levs))
      allocate (Diag%dtauy2d_fd  (IM,Model%levs))
      Diag%dtaux2d_ls    = clear_val
      Diag%dtauy2d_ls    = clear_val
      Diag%dtaux2d_bl    = clear_val
      Diag%dtauy2d_bl    = clear_val
      Diag%dtaux2d_ss    = clear_val
      Diag%dtauy2d_ss    = clear_val
      Diag%dtaux2d_fd    = clear_val
      Diag%dtauy2d_fd    = clear_val
      allocate (Diag%dusfc_ls    (IM))
      allocate (Diag%dvsfc_ls    (IM))
      allocate (Diag%dusfc_bl    (IM))
      allocate (Diag%dvsfc_bl    (IM))
      allocate (Diag%dusfc_ss    (IM))
      allocate (Diag%dvsfc_ss    (IM))
      allocate (Diag%dusfc_fd    (IM))
      allocate (Diag%dvsfc_fd    (IM))
      Diag%dusfc_ls      = 0
      Diag%dvsfc_ls      = 0
      Diag%dusfc_bl      = 0
      Diag%dvsfc_bl      = 0
      Diag%dusfc_ss      = 0
      Diag%dvsfc_ss      = 0
      Diag%dusfc_fd      = 0
      Diag%dvsfc_fd      = 0
    endif
#endif

    !--- diagnostics for coupled chemistry
    if (Model%cplchm) call Diag%chem_init(IM,Model)

    call Diag%rad_zero  (Model)
!    if(Model%me==0) print *,'in diag_create, call rad_zero'
    linit = .true.
    call Diag%phys_zero (Model, linit=linit)
!    if(Model%me==0) print *,'in diag_create, call phys_zero'
    linit = .false.

  end subroutine diag_create

!-----------------------
! GFS_diag%rad_zero
!-----------------------
  subroutine diag_rad_zero(Diag, Model)
    class(GFS_diag_type)               :: Diag
    type(GFS_control_type), intent(in) :: Model

    Diag%fluxr        = zero
    Diag%topfsw%upfxc = zero
    Diag%topfsw%dnfxc = zero
    Diag%topfsw%upfx0 = zero
    Diag%topflw%upfxc = zero
    Diag%topflw%upfx0 = zero
    if (Model%ldiag3d) then
      Diag%cldcov     = zero
    endif

  end subroutine diag_rad_zero

!------------------------
! GFS_diag%phys_zero
!------------------------
  subroutine diag_phys_zero (Diag, Model, linit, iauwindow_center)
    class(GFS_diag_type)               :: Diag
    type(GFS_control_type), intent(in) :: Model
    logical,optional, intent(in)       :: linit, iauwindow_center

    logical set_totprcp

    !--- In/Out
    Diag%srunoff    = zero
    Diag%evbsa      = zero
    Diag%evcwa      = zero
    Diag%snohfa     = zero
    Diag%transa     = zero
    Diag%sbsnoa     = zero
    Diag%snowca     = zero
    Diag%soilm      = zero
    Diag%tmpmin     = huge
    Diag%tmpmax     = zero
    Diag%dusfc      = zero
    Diag%dvsfc      = zero
    Diag%dtsfc      = zero
    Diag%dqsfc      = zero
    Diag%gflux      = zero
    Diag%dlwsfc     = zero
    Diag%ulwsfc     = zero
    Diag%suntim     = zero
    Diag%runoff     = zero
    Diag%ep         = zero
    Diag%cldwrk     = zero
    Diag%dugwd      = zero
    Diag%dvgwd      = zero
    Diag%psmean     = zero
    Diag%spfhmin    = huge
    Diag%spfhmax    = zero
    Diag%u10mmax    = zero
    Diag%v10mmax    = zero
    Diag%wind10mmax = zero
    Diag%u10max     = zero
    Diag%v10max     = zero
    Diag%spd10max   = zero
    Diag%rain       = zero
    Diag%rainc      = zero
    Diag%ice        = zero
    Diag%snow       = zero
    Diag%graupel    = zero

    !--- Out
    Diag%u10m       = zero
    Diag%v10m       = zero
    Diag%dpt2m      = zero
    Diag%zlvl       = zero
    Diag%psurf      = zero
    Diag%pwat       = zero
    Diag%t1         = zero
    Diag%q1         = zero
    Diag%u1         = zero
    Diag%v1         = zero
    Diag%chh        = zero
    Diag%cmm        = zero
    Diag%dlwsfci    = zero
    Diag%ulwsfci    = zero
    Diag%dswsfci    = zero
#ifdef CCPP
    Diag%nswsfci    = zero
#endif
    Diag%uswsfci    = zero
    Diag%dusfci     = zero
    Diag%dvsfci     = zero
    Diag%dtsfci     = zero
    Diag%dqsfci     = zero
    Diag%gfluxi     = zero
    Diag%epi        = zero
    Diag%smcwlt2    = zero
    Diag%smcref2    = zero
    if (.not. Model%lsm == Model%lsm_ruc) then
      Diag%wet1       = zero
    end if
    Diag%sr         = zero
    Diag%tdomr      = zero
    Diag%tdomzr     = zero
    Diag%tdomip     = zero
    Diag%tdoms      = zero
    Diag%skebu_wts  = zero
    Diag%skebv_wts  = zero
    Diag%sppt_wts   = zero
    Diag%shum_wts   = zero
    Diag%zmtnblck   = zero

#ifdef CCPP
    if (Model%imp_physics == Model%imp_physics_fer_hires) then
       Diag%TRAIN      = zero
    end if
#endif
    Diag%totprcpb   = zero
    Diag%cnvprcpb   = zero
    Diag%toticeb    = zero
    Diag%totsnwb    = zero
    Diag%totgrpb    = zero
!
    if (Model%do_ca) then
      Diag%ca_out   = zero
      Diag%ca_deep  = zero
      Diag%ca_turb  = zero
      Diag%ca_shal  = zero
      Diag%ca_rad   = zero
      Diag%ca_micro = zero
    endif
!    if(Model%me == Model%master) print *,'in diag_phys_zero, totprcpb set to 0,kdt=',Model%kdt

    if (Model%ldiag3d) then
      Diag%du3dt    = zero
      Diag%dv3dt    = zero
      Diag%dt3dt    = zero
!     Diag%dq3dt    = zero
!     Diag%upd_mf   = zero
!     Diag%dwn_mf   = zero
!     Diag%det_mf   = zero
    endif

!
!-----------------------------
    if (Model%ldiag_ugwp)   then
      if(Model%me == Model%master) print *,'VAY in diag_phys_zero at kdt=',Model%kdt, Model%ldiag_ugwp
      Diag%du3dt_pbl   = zero
      Diag%dv3dt_pbl   = zero
      Diag%dt3dt_pbl   = zero
!
      Diag%du3dt_ogw   = zero
      Diag%dv3dt_ogw   = zero
      Diag%dt3dt_ogw   = zero

      Diag%du3dt_mtb   = zero
      Diag%dv3dt_mtb   = zero
      Diag%dt3dt_mtb   = zero

      Diag%du3dt_tms   = zero
      Diag%dv3dt_tms   = zero
      Diag%dt3dt_tms   = zero

      Diag%du3dt_ngw   = zero
      Diag%dv3dt_ngw   = zero
      Diag%dt3dt_ngw   = zero

      Diag%du3dt_moist = zero
      Diag%dv3dt_moist = zero
      Diag%dt3dt_moist = zero

      Diag%dudt_tot    = zero
      Diag%dvdt_tot    = zero
      Diag%dtdt_tot    = zero

      Diag%uav_ugwp    = zero
      Diag%tav_ugwp    = zero
!COORDE
      Diag%du3dt_dyn   = zero
      Diag%zmtb        = zero
      Diag%zogw        = zero
      Diag%zlwb        = zero

      Diag%tau_mtb     = zero
      Diag%tau_ogw     = zero
      Diag%tau_ngw     = zero
      Diag%tau_tofd    = zero
    endif
!
    if (Model%do_ugwp)   then
      Diag%gwp_ax     = zero
      Diag%gwp_ay     = zero
      Diag%gwp_dtdt   = zero
      Diag%gwp_kdis   = zero
      Diag%gwp_axo    = zero
      Diag%gwp_ayo    = zero
      Diag%gwp_axc    = zero
      Diag%gwp_ayc    = zero
      Diag%gwp_axf    = zero
      Diag%gwp_ayf    = zero
      Diag%gwp_dcheat = zero
      Diag%gwp_scheat = zero
      Diag%gwp_precip = zero
      Diag%gwp_klevs  = -99
      Diag%gwp_fgf    = zero
      Diag%gwp_okw    = zero
    endif
!-----------------------------

! max hourly diagnostics
    Diag%refl_10cm   = zero
    Diag%refdmax     = -35.
    Diag%refdmax263k = -35.
    Diag%t02max      = -999.
    Diag%t02min      = 999.
    Diag%rh02max     = -999.
    Diag%rh02min     = 999.
    set_totprcp      = .false.
    if (present(linit) ) set_totprcp = linit
    if (present(iauwindow_center) ) set_totprcp = iauwindow_center
    if (set_totprcp) then
      Diag%totprcp = zero
      Diag%cnvprcp = zero
      Diag%totice  = zero
      Diag%totsnw  = zero
      Diag%totgrp  = zero
    endif
  end subroutine diag_phys_zero

!-----------------------
! GFS_diag%chem_init
!-----------------------
  subroutine diag_chem_init(Diag, IM, Model)

    use parse_tracers,    only: get_tracer_index, NO_TRACER

    class(GFS_diag_type)               :: Diag
    integer,                intent(in) :: IM
    type(GFS_control_type), intent(in) :: Model

    ! -- local variables
    integer :: n

    ! -- initialize diagnostic variables depending on
    ! -- specific chemical tracers
    if (Model%ntchm > 0) then
      ! -- retrieve number of dust bins
      n = get_number_bins('dust')
#ifdef CCPP
      Diag%ndust = n
#endif
      if (n > 0) then
        allocate (Diag%duem(IM,n))
        Diag%duem = zero
      end if

      ! -- retrieve number of sea salt bins
      n = get_number_bins('seas')
#ifdef CCPP
      Diag%nseasalt = n
#endif
      if (n > 0) then
        allocate (Diag%ssem(IM,n))
        Diag%ssem = zero
      end if
    end if

    ! -- sedimentation and dry/wet deposition diagnostics
    if (associated(Model%ntdiag)) then
      ! -- get number of tracers with enabled diagnostics
      n = count(Model%ntdiag)
#ifdef CCPP
      Diag%ntchmdiag = n
#endif
      ! -- initialize sedimentation
      allocate (Diag%sedim(IM,n))
      Diag%sedim = zero

      ! -- initialize dry deposition
      allocate (Diag%drydep(IM,n))
      Diag%drydep = zero

      ! -- initialize large-scale wet deposition
      allocate (Diag%wetdpl(IM,n))
      Diag%wetdpl = zero

      ! -- initialize convective-scale wet deposition
      allocate (Diag%wetdpc(IM,n))
      Diag%wetdpc = zero
    end if

    ! -- initialize anthropogenic and biomass
    ! -- burning emission diagnostics for
    ! -- (in order): black carbon,
    ! -- organic carbon, and sulfur dioxide
    allocate (Diag%abem(IM,6))
    Diag%abem = zero

    ! -- initialize column burden diagnostics
    ! -- for aerosol species (in order): pm2.5
    ! -- black carbon, organic carbon, sulfate,
    ! -- dust, sea salt
    allocate (Diag%aecm(IM,6))
    Diag%aecm = zero

  contains

    integer function get_number_bins(tracer_type)
      character(len=*), intent(in) :: tracer_type

      logical :: next
      integer :: n
      character(len=5) :: name

      get_number_bins = 0

      n = 0
      next = .true.
      do while (next)
        n = n + 1
        write(name,'(a,i1)') tracer_type, n + 1
        next = get_tracer_index(Model%tracer_names, name, &
          Model%me, Model%master, Model%debug) /= NO_TRACER
      end do

      get_number_bins = n

    end function get_number_bins

  end subroutine diag_chem_init

#ifdef CCPP
  !-------------------------
  ! GFS_interstitial_type%create
  !-------------------------
  subroutine interstitial_create (Interstitial, IM, Model)
    !
    implicit none
    !
    class(GFS_interstitial_type)       :: Interstitial
    integer,                intent(in) :: IM
    type(GFS_control_type), intent(in) :: Model
    !
    allocate (Interstitial%otspt      (Model%ntracp1,2))
    ! Set up numbers of tracers for PBL, convection, etc: sets
    ! Interstitial%{nncl,nvdiff,mg3_as_mg2,nn,tracers_total,ntiwx,ntk,ntkev,otspt,nsamftrac,ncstrac,nscav}
    call interstitial_setup_tracers(Interstitial, Model)
    ! Allocate arrays
    allocate (Interstitial%adjsfculw_land  (IM))
    allocate (Interstitial%adjsfculw_ice   (IM))
    allocate (Interstitial%adjsfculw_ocean (IM))
    allocate (Interstitial%adjnirbmd       (IM))
    allocate (Interstitial%adjnirbmu       (IM))
    allocate (Interstitial%adjnirdfd       (IM))
    allocate (Interstitial%adjnirdfu       (IM))
    allocate (Interstitial%adjvisbmd       (IM))
    allocate (Interstitial%adjvisbmu       (IM))
    allocate (Interstitial%adjvisdfu       (IM))
    allocate (Interstitial%adjvisdfd       (IM))
    allocate (Interstitial%aerodp          (IM,NSPC1))
    allocate (Interstitial%alb1d           (IM))
    allocate (Interstitial%bexp1d          (IM))
    allocate (Interstitial%cd              (IM))
    allocate (Interstitial%cd_ice          (IM))
    allocate (Interstitial%cd_land         (IM))
    allocate (Interstitial%cd_ocean        (IM))
    allocate (Interstitial%cdq             (IM))
    allocate (Interstitial%cdq_ice         (IM))
    allocate (Interstitial%cdq_land        (IM))
    allocate (Interstitial%cdq_ocean       (IM))
    allocate (Interstitial%chh_ice         (IM))
    allocate (Interstitial%chh_land        (IM))
    allocate (Interstitial%chh_ocean       (IM))
    allocate (Interstitial%cldf            (IM))
    allocate (Interstitial%cldsa           (IM,5))
    allocate (Interstitial%cldtaulw        (IM,Model%levr+LTP))
    allocate (Interstitial%cldtausw        (IM,Model%levr+LTP))
    allocate (Interstitial%cld1d           (IM))
    allocate (Interstitial%clouds          (IM,Model%levr+LTP,NF_CLDS))
    allocate (Interstitial%clw             (IM,Model%levs,Interstitial%nn))
    allocate (Interstitial%clx             (IM,4))
    allocate (Interstitial%cmm_ice         (IM))
    allocate (Interstitial%cmm_land        (IM))
    allocate (Interstitial%cmm_ocean       (IM))
    allocate (Interstitial%cnvc            (IM,Model%levs))
    allocate (Interstitial%cnvw            (IM,Model%levs))
    allocate (Interstitial%ctei_r          (IM))
    allocate (Interstitial%ctei_rml        (IM))
    allocate (Interstitial%cumabs          (IM))
    allocate (Interstitial%dd_mf           (IM,Model%levs))
    allocate (Interstitial%de_lgth         (IM))
    allocate (Interstitial%del             (IM,Model%levs))
    allocate (Interstitial%del_gz          (IM,Model%levs+1))
    allocate (Interstitial%delr            (IM,Model%levr+LTP))
    allocate (Interstitial%dkt             (IM,Model%levs-1))
    allocate (Interstitial%dlength         (IM))
    allocate (Interstitial%dqdt            (IM,Model%levs,Model%ntrac))
    allocate (Interstitial%dqsfc1          (IM))
    allocate (Interstitial%drain           (IM))
    allocate (Interstitial%dtdt            (IM,Model%levs))
    allocate (Interstitial%dtdtc           (IM,Model%levs))
    allocate (Interstitial%dtsfc1          (IM))
    allocate (Interstitial%dt_mf           (IM,Model%levs))
    allocate (Interstitial%dtzm            (IM))
    allocate (Interstitial%dudt            (IM,Model%levs))
    allocate (Interstitial%dusfcg          (IM))
    allocate (Interstitial%dusfc1          (IM))
    allocate (Interstitial%dvdt            (IM,Model%levs))
    allocate (Interstitial%dvsfcg          (IM))
    allocate (Interstitial%dvsfc1          (IM))
    allocate (Interstitial%dvdftra         (IM,Model%levs,Interstitial%nvdiff))
    allocate (Interstitial%dzlyr           (IM,Model%levr+LTP))
    allocate (Interstitial%elvmax          (IM))
    allocate (Interstitial%ep1d            (IM))
    allocate (Interstitial%ep1d_ice        (IM))
    allocate (Interstitial%ep1d_land       (IM))
    allocate (Interstitial%ep1d_ocean      (IM))
    allocate (Interstitial%evap_ice        (IM))
    allocate (Interstitial%evap_land       (IM))
    allocate (Interstitial%evap_ocean      (IM))
    allocate (Interstitial%evbs            (IM))
    allocate (Interstitial%evcw            (IM))
    allocate (Interstitial%faerlw          (IM,Model%levr+LTP,NBDLW,NF_AELW))
    allocate (Interstitial%faersw          (IM,Model%levr+LTP,NBDSW,NF_AESW))
    allocate (Interstitial%ffhh_ice        (IM))
    allocate (Interstitial%ffhh_land       (IM))
    allocate (Interstitial%ffhh_ocean      (IM))
    allocate (Interstitial%fh2             (IM))
    allocate (Interstitial%fh2_ice         (IM))
    allocate (Interstitial%fh2_land        (IM))
    allocate (Interstitial%fh2_ocean       (IM))
    allocate (Interstitial%flag_cice       (IM))
    allocate (Interstitial%flag_guess      (IM))
    allocate (Interstitial%flag_iter       (IM))
    allocate (Interstitial%ffmm_ice        (IM))
    allocate (Interstitial%ffmm_land       (IM))
    allocate (Interstitial%ffmm_ocean      (IM))
    allocate (Interstitial%fm10            (IM))
    allocate (Interstitial%fm10_ice        (IM))
    allocate (Interstitial%fm10_land       (IM))
    allocate (Interstitial%fm10_ocean      (IM))
    allocate (Interstitial%frland          (IM))
    allocate (Interstitial%fscav           (Interstitial%nscav))
    allocate (Interstitial%fswtr           (Interstitial%nscav))
    allocate (Interstitial%gabsbdlw        (IM))
    allocate (Interstitial%gabsbdlw_ice    (IM))
    allocate (Interstitial%gabsbdlw_land   (IM))
    allocate (Interstitial%gabsbdlw_ocean  (IM))
    allocate (Interstitial%gamma           (IM))
    allocate (Interstitial%gamq            (IM))
    allocate (Interstitial%gamt            (IM))
    allocate (Interstitial%gasvmr          (IM,Model%levr+LTP,NF_VGAS))
    allocate (Interstitial%gflx            (IM))
    allocate (Interstitial%gflx_ice        (IM))
    allocate (Interstitial%gflx_land       (IM))
    allocate (Interstitial%gflx_ocean      (IM))
    allocate (Interstitial%gwdcu           (IM,Model%levs))
    allocate (Interstitial%gwdcv           (IM,Model%levs))
    allocate (Interstitial%h2o_pres        (levh2o))
    allocate (Interstitial%hflx_ice        (IM))
    allocate (Interstitial%hflx_land       (IM))
    allocate (Interstitial%hflx_ocean      (IM))
    allocate (Interstitial%dry             (IM))
    allocate (Interstitial%idxday          (IM))
    allocate (Interstitial%icy             (IM))
    allocate (Interstitial%lake            (IM))
    allocate (Interstitial%ocean           (IM))
    allocate (Interstitial%islmsk          (IM))
    allocate (Interstitial%islmsk_cice      (IM))
    allocate (Interstitial%wet             (IM))
    allocate (Interstitial%kbot            (IM))
    allocate (Interstitial%kcnv            (IM))
    allocate (Interstitial%kinver          (IM))
    allocate (Interstitial%kpbl            (IM))
    allocate (Interstitial%ktop            (IM))
    allocate (Interstitial%mbota           (IM,3))
    allocate (Interstitial%mtopa           (IM,3))
    allocate (Interstitial%oa4             (IM,4))
    allocate (Interstitial%oc              (IM))
    allocate (Interstitial%olyr            (IM,Model%levr+LTP))
    allocate (Interstitial%oz_pres         (levozp))
    allocate (Interstitial%plvl            (IM,Model%levr+1+LTP))
    allocate (Interstitial%plyr            (IM,Model%levr+LTP))
    allocate (Interstitial%prnum           (IM,Model%levs))
    allocate (Interstitial%qlyr            (IM,Model%levr+LTP))
    allocate (Interstitial%prcpmp          (IM))
    allocate (Interstitial%qss             (IM))
    allocate (Interstitial%raincd          (IM))
    allocate (Interstitial%raincs          (IM))
    allocate (Interstitial%rainmcadj       (IM))
    allocate (Interstitial%rainp           (IM,Model%levs))
    allocate (Interstitial%rb              (IM))
    allocate (Interstitial%rb_ice          (IM))
    allocate (Interstitial%rb_land         (IM))
    allocate (Interstitial%rb_ocean        (IM))
    allocate (Interstitial%rhc             (IM,Model%levs))
    allocate (Interstitial%runoff          (IM))
    allocate (Interstitial%save_q          (IM,Model%levs,Model%ntrac))
    allocate (Interstitial%save_t          (IM,Model%levs))
    allocate (Interstitial%save_tcp        (IM,Model%levs))
    allocate (Interstitial%save_u          (IM,Model%levs))
    allocate (Interstitial%save_v          (IM,Model%levs))
    allocate (Interstitial%sbsno           (IM))
    allocate (Interstitial%scmpsw          (IM))
    allocate (Interstitial%semis_ice       (IM))
    allocate (Interstitial%semis_land      (IM))
    allocate (Interstitial%semis_ocean     (IM))
    allocate (Interstitial%sfcalb          (IM,NF_ALBD))
    allocate (Interstitial%sigma           (IM))
    allocate (Interstitial%sigmaf          (IM))
    allocate (Interstitial%sigmafrac       (IM,Model%levs))
    allocate (Interstitial%sigmatot        (IM,Model%levs))
    allocate (Interstitial%slopetype       (IM))
    allocate (Interstitial%snowc           (IM))
    allocate (Interstitial%snohf           (IM))
    allocate (Interstitial%snowmt          (IM))
    allocate (Interstitial%soiltype        (IM))
    allocate (Interstitial%stress          (IM))
    allocate (Interstitial%stress_ice      (IM))
    allocate (Interstitial%stress_land     (IM))
    allocate (Interstitial%stress_ocean    (IM))
    allocate (Interstitial%theta           (IM))
    allocate (Interstitial%tice            (IM))
    allocate (Interstitial%tlvl            (IM,Model%levr+1+LTP))
    allocate (Interstitial%tlyr            (IM,Model%levr+LTP))
    allocate (Interstitial%tprcp_ice       (IM))
    allocate (Interstitial%tprcp_land      (IM))
    allocate (Interstitial%tprcp_ocean     (IM))
    allocate (Interstitial%trans           (IM))
    allocate (Interstitial%tseal           (IM))
    allocate (Interstitial%tsfa            (IM))
    allocate (Interstitial%tsfg            (IM))
    allocate (Interstitial%tsurf           (IM))
    allocate (Interstitial%ud_mf           (IM,Model%levs))
    allocate (Interstitial%ulwsfc_cice     (IM))
    allocate (Interstitial%dusfc_cice      (IM))
    allocate (Interstitial%dvsfc_cice      (IM))
    allocate (Interstitial%dtsfc_cice      (IM))
    allocate (Interstitial%dqsfc_cice      (IM))
    allocate (Interstitial%vdftra          (IM,Model%levs,Interstitial%nvdiff))  !GJF first dimension was set as 'IX' in GFS_physics_driver
    allocate (Interstitial%vegf1d          (IM))
    allocate (Interstitial%vegtype         (IM))
    allocate (Interstitial%wcbmax          (IM))
    allocate (Interstitial%weasd_ice       (IM))
    allocate (Interstitial%weasd_land      (IM))
    allocate (Interstitial%weasd_ocean     (IM))
    allocate (Interstitial%wind            (IM))
    allocate (Interstitial%work1           (IM))
    allocate (Interstitial%work2           (IM))
    allocate (Interstitial%work3           (IM))
    allocate (Interstitial%xcosz           (IM))
    allocate (Interstitial%xlai1d          (IM))
    allocate (Interstitial%xmu             (IM))
    allocate (Interstitial%z01d            (IM))
    allocate (Interstitial%zt1d            (IM))
! CIRES UGWP v0
    allocate (Interstitial%gw_dudt         (IM,Model%levs))
    allocate (Interstitial%gw_dvdt         (IM,Model%levs))
    allocate (Interstitial%gw_dtdt         (IM,Model%levs))
    allocate (Interstitial%gw_kdis         (IM,Model%levs))
    allocate (Interstitial%tau_mtb         (IM))
    allocate (Interstitial%tau_ogw         (IM))
    allocate (Interstitial%tau_tofd        (IM))
    allocate (Interstitial%tau_ngw         (IM))
    allocate (Interstitial%zmtb            (IM))
    allocate (Interstitial%zlwb            (IM))
    allocate (Interstitial%zogw            (IM))
    allocate (Interstitial%dudt_mtb        (IM,Model%levs))
    allocate (Interstitial%dudt_ogw        (IM,Model%levs))
    allocate (Interstitial%dudt_tms        (IM,Model%levs))
!
    ! Allocate arrays that are conditional on physics choices
    if (Model%imp_physics == Model%imp_physics_gfdl .or. Model%imp_physics == Model%imp_physics_thompson) then
       allocate (Interstitial%graupelmp  (IM))
       allocate (Interstitial%icemp      (IM))
       allocate (Interstitial%rainmp     (IM))
       allocate (Interstitial%snowmp     (IM))
    else if (Model%imp_physics == Model%imp_physics_mg) then
       allocate (Interstitial%ncgl       (IM,Model%levs))
       allocate (Interstitial%ncpr       (IM,Model%levs))
       allocate (Interstitial%ncps       (IM,Model%levs))
       allocate (Interstitial%qgl        (IM,Model%levs))
       allocate (Interstitial%qrn        (IM,Model%levs))
       allocate (Interstitial%qsnw       (IM,Model%levs))
       allocate (Interstitial%qlcn       (IM,Model%levs))
       allocate (Interstitial%qicn       (IM,Model%levs))
       allocate (Interstitial%w_upi      (IM,Model%levs))
       allocate (Interstitial%cf_upi     (IM,Model%levs))
       allocate (Interstitial%cnv_mfd    (IM,Model%levs))
       allocate (Interstitial%cnv_dqldt  (IM,Model%levs))
       allocate (Interstitial%clcn       (IM,Model%levs))
       allocate (Interstitial%cnv_fice   (IM,Model%levs))
       allocate (Interstitial%cnv_ndrop  (IM,Model%levs))
       allocate (Interstitial%cnv_nice   (IM,Model%levs))
    end if
    if (Model%imp_physics == Model%imp_physics_fer_hires) then
    !--- if HWRF physics?
       allocate (Interstitial%qv_r        (IM,Model%levs))
       allocate (Interstitial%qc_r        (IM,Model%levs))
       allocate (Interstitial%qi_r        (IM,Model%levs))
       allocate (Interstitial%qr_r        (IM,Model%levs))
       allocate (Interstitial%qs_r        (IM,Model%levs))
       allocate (Interstitial%qg_r        (IM,Model%levs))

    !--- Ferrier-Aligo MP scheme
       allocate (Interstitial%f_ice       (IM,Model%levs))
       allocate (Interstitial%f_rain      (IM,Model%levs))
       allocate (Interstitial%f_rimef     (IM,Model%levs))
       allocate (Interstitial%cwm         (IM,Model%levs))
    end if
    if (Model%do_shoc) then
       if (.not. associated(Interstitial%qrn))  allocate (Interstitial%qrn  (IM,Model%levs))
       if (.not. associated(Interstitial%qsnw)) allocate (Interstitial%qsnw (IM,Model%levs))
       ! DH* updated version of shoc from May 22 2019 (not yet in CCPP) doesn't use qgl? remove?
       if (.not. associated(Interstitial%qgl))  allocate (Interstitial%qgl  (IM,Model%levs))
       ! *DH
       allocate (Interstitial%ncpi (IM,Model%levs))
       allocate (Interstitial%ncpl (IM,Model%levs))
    end if
    if (Model%lsm == Model%lsm_noahmp) then
       allocate (Interstitial%t2mmp (IM))
       allocate (Interstitial%q2mp  (IM))
    end if
    !
    ! Set components that do not change
    Interstitial%frain            = Model%dtf/Model%dtp
    Interstitial%ipr              = min(IM,10)
    Interstitial%latidxprnt       = 1
    Interstitial%levi             = Model%levs+1
    Interstitial%levh2o           = levh2o
    Interstitial%levozp           = levozp
    Interstitial%lmk              = Model%levr+LTP
    Interstitial%lmp              = Model%levr+1+LTP
    Interstitial%h2o_coeff        = h2o_coeff
    Interstitial%nbdlw            = NBDLW
    Interstitial%nbdsw            = NBDSW
    Interstitial%nf_aelw          = NF_AELW
    Interstitial%nf_aesw          = NF_AESW
    Interstitial%nspc1            = NSPC1
    Interstitial%oz_coeff         = oz_coeff
    Interstitial%oz_coeffp5       = oz_coeff+5
    ! h2o_pres and oz_pres do not change during the run, but
    ! need to be set later in GFS_phys_time_vary_init (after
    ! h2o_pres/oz_pres are read in read_h2odata/read_o3data)
    Interstitial%h2o_pres         = clear_val
    Interstitial%oz_pres          = clear_val
    !
    Interstitial%skip_macro       = .false.
    ! The value phys_hydrostatic from dynamics does not match the
    ! hardcoded value for calling GFDL MP in GFS_physics_driver.F90,
    ! which is set to .true.
    Interstitial%phys_hydrostatic = .true.
    !
    ! Reset all other variables
    call Interstitial%rad_reset (Model)
    call Interstitial%phys_reset (Model)
    !
  end subroutine interstitial_create

  subroutine interstitial_setup_tracers(Interstitial, Model)
    !
    implicit none
    !
    class(GFS_interstitial_type)       :: Interstitial
    type(GFS_control_type), intent(in) :: Model
    integer :: n, tracers

    !first, initialize the values (in case the values don't get initialized within if statements below)
    Interstitial%nncl             = Model%ncld
    Interstitial%nvdiff           = Model%ntrac
    Interstitial%mg3_as_mg2       = .false.
    Interstitial%nn               = Model%ntrac + 1
    Interstitial%itc              = 0
    Interstitial%ntk              = 0
    Interstitial%ntkev            = 0
    Interstitial%tracers_total    = 0
    Interstitial%otspt(:,:)       = .true.
    Interstitial%nsamftrac        = 0
    Interstitial%ncstrac          = 0
    Interstitial%nscav            = Model%ntrac-Model%ncld+2

    ! perform aerosol convective transport and PBL diffusion
    Interstitial%trans_aero = Model%cplchm .and. Model%trans_trac

    if (Model%imp_physics == Model%imp_physics_thompson) then
      if (Model%ltaerosol) then
        Interstitial%nvdiff = 12
      else
        Interstitial%nvdiff = 9
      endif
      if (Model%satmedmf) Interstitial%nvdiff = Interstitial%nvdiff + 1
      Interstitial%nncl = 5
    elseif (Model%imp_physics == Model%imp_physics_wsm6) then
      Interstitial%nvdiff = Model%ntrac -3
      if (Model%satmedmf) Interstitial%nvdiff = Interstitial%nvdiff + 1
      Interstitial%nncl = 5
    elseif (Model%ntclamt > 0) then             ! for GFDL MP don't diffuse cloud amount
      Interstitial%nvdiff = Model%ntrac - 1
    endif

    if (Model%imp_physics == Model%imp_physics_gfdl) then
      Interstitial%nncl = 5
    endif

    if (Model%imp_physics == Model%imp_physics_mg) then
      if (abs(Model%fprcp) == 1) then
        Interstitial%nncl = 4                          ! MG2 with rain and snow
        Interstitial%mg3_as_mg2 = .false.
      elseif (Model%fprcp >= 2) then
        if(Model%ntgl > 0 .and. (Model%mg_do_graupel .or. Model%mg_do_hail)) then
          Interstitial%nncl = 5                        ! MG3 with rain and snow and grapuel/hail
          Interstitial%mg3_as_mg2 = .false.
        else                              ! MG3 code run without graupel/hail i.e. as MG2
          Interstitial%nncl = 4
          Interstitial%mg3_as_mg2 = .true.
        endif
      endif
    endif

    ! DH* STILL VALID GIVEN THE CHANGES BELOW FOR CPLCHM?
    if (Interstitial%nvdiff == Model%ntrac) then
      Interstitial%ntiwx = Model%ntiw
    else
      if (Model%imp_physics == Model%imp_physics_wsm6) then
        Interstitial%ntiwx = 3
      elseif (Model%imp_physics == Model%imp_physics_thompson) then
        if(Model%ltaerosol) then
          Interstitial%ntiwx = 3
        else
          Interstitial%ntiwx = 3
        endif
      elseif (Model%imp_physics == Model%imp_physics_gfdl) then
        Interstitial%ntiwx = 3
      ! F-A MP scheme
      elseif (Model%imp_physics == Model%imp_physics_fer_hires) then
        Interstitial%ntiwx = 3 ! total ice or total condensate
      elseif (Model%imp_physics == Model%imp_physics_mg) then
        Interstitial%ntiwx = 3
      else
        Interstitial%ntiwx = 0
      endif
    endif
    ! *DH

    if (Model%cplchm) then
      ! Only Zhao/Carr/Sundqvist and GFDL microphysics schemes are supported
      ! when coupling with chemistry. PBL diffusion of aerosols is only supported
      ! for GFDL microphysics and MG microphysics.
      if (Model%imp_physics == Model%imp_physics_zhao_carr) then
        Interstitial%nvdiff = 3
      elseif (Model%imp_physics == Model%imp_physics_mg) then
        if (Model%ntgl > 0) then
          Interstitial%nvdiff = 12
        else
          Interstitial%nvdiff = 10
        endif
      elseif (Model%imp_physics == Model%imp_physics_gfdl) then
        Interstitial%nvdiff = 7
      else
        write(0,*) "Only Zhao/Carr/Sundqvist and GFDL microphysics schemes are supported when coupling with chemistry"
        stop
      endif
      if (Interstitial%trans_aero) Interstitial%nvdiff = Interstitial%nvdiff + Model%ntchm
      if (Model%ntke > 0) Interstitial%nvdiff = Interstitial%nvdiff + 1    ! adding tke to the list
    endif

    Interstitial%ntkev = Interstitial%nvdiff

    if (Model%ntiw > 0) then
      if (Model%ntclamt > 0) then
        Interstitial%nn = Model%ntrac - 2
      else
        Interstitial%nn = Model%ntrac - 1
      endif
    elseif (Model%ntcw > 0) then
      Interstitial%nn = Model%ntrac
    else
      Interstitial%nn = Model%ntrac + 1
    endif

    if (Model%cscnv .or. Model%satmedmf .or. Model%trans_trac ) then
      Interstitial%otspt(:,:)   = .true.     ! otspt is used only for cscnv
      Interstitial%otspt(1:3,:) = .false.    ! this is for sp.hum, ice and liquid water
      tracers = 2
      do n=2,Model%ntrac
        if ( n /= Model%ntcw  .and. n /= Model%ntiw  .and. n /= Model%ntclamt .and. &
             n /= Model%ntrw  .and. n /= Model%ntsw  .and. n /= Model%ntrnc   .and. &
             n /= Model%ntsnc .and. n /= Model%ntgl  .and. n /= Model%ntgnc) then
          tracers = tracers + 1
          if (Model%ntke  == n ) then
            Interstitial%otspt(tracers+1,1) = .false.
            Interstitial%ntk = tracers
          endif
          if (Model%ntlnc == n .or. Model%ntinc == n .or. Model%ntrnc == n .or. Model%ntsnc == n .or. Model%ntgnc == n)    &
!           if (ntlnc == n .or. ntinc == n .or. ntrnc == n .or. ntsnc == n .or.&
!               ntrw  == n .or. ntsw  == n .or. ntgl  == n)                    &
                  Interstitial%otspt(tracers+1,1) = .false.
          if (Interstitial%trans_aero .and. Model%ntchs == n) Interstitial%itc = tracers
        endif
      enddo
      Interstitial%tracers_total = tracers - 2
    endif   ! end if_ras or cfscnv or samf
    if(.not. Model%satmedmf .and. .not. Model%trans_trac) then
       Interstitial%nsamftrac = 0
    else
       Interstitial%nsamftrac = Interstitial%tracers_total
    endif
    Interstitial%ncstrac = Interstitial%tracers_total + 3

  end subroutine interstitial_setup_tracers

  subroutine interstitial_rad_reset (Interstitial, Model)
    !
    implicit none
    !
    class(GFS_interstitial_type) :: Interstitial
    type(GFS_control_type), intent(in) :: Model
    !
    Interstitial%aerodp       = clear_val
    Interstitial%alb1d        = clear_val
    Interstitial%cldsa        = clear_val
    Interstitial%cldtaulw     = clear_val
    Interstitial%cldtausw     = clear_val
    Interstitial%clouds       = clear_val
    Interstitial%de_lgth      = clear_val
    Interstitial%delr         = clear_val
    Interstitial%dzlyr        = clear_val
    Interstitial%faerlw       = clear_val
    Interstitial%faersw       = clear_val
    Interstitial%gasvmr       = clear_val
    Interstitial%idxday       = 0
    Interstitial%kb           = 0
    Interstitial%kd           = 0
    Interstitial%kt           = 0
    Interstitial%mbota        = 0
    Interstitial%mtopa        = 0
    Interstitial%nday         = 0
    Interstitial%olyr         = clear_val
    Interstitial%plvl         = clear_val
    Interstitial%plyr         = clear_val
    Interstitial%qlyr         = clear_val
    Interstitial%raddt        = clear_val
    Interstitial%scmpsw%uvbfc = clear_val
    Interstitial%scmpsw%uvbf0 = clear_val
    Interstitial%scmpsw%nirbm = clear_val
    Interstitial%scmpsw%nirdf = clear_val
    Interstitial%scmpsw%visbm = clear_val
    Interstitial%scmpsw%visdf = clear_val
    Interstitial%sfcalb       = clear_val
    Interstitial%tlvl         = clear_val
    Interstitial%tlyr         = clear_val
    Interstitial%tsfa         = clear_val
    Interstitial%tsfg         = clear_val

! F-A scheme
    if (Model%imp_physics == Model%imp_physics_fer_hires) then
         Interstitial%qv_r       = clear_val
         Interstitial%qc_r       = clear_val
         Interstitial%qi_r       = clear_val
         Interstitial%qr_r       = clear_val
         Interstitial%qs_r       = clear_val
         Interstitial%qg_r       = clear_val
       if(Model%spec_adv) then
         Interstitial%f_ice     = clear_val
         Interstitial%f_rain    = clear_val
         Interstitial%f_rimef   = clear_val
         Interstitial%cwm       = clear_val
       end if
    end if

    !
  end subroutine interstitial_rad_reset

  subroutine interstitial_phys_reset (Interstitial, Model)
    !
    implicit none
    !
    class(GFS_interstitial_type) :: Interstitial
    type(GFS_control_type), intent(in) :: Model
    !
    Interstitial%adjsfculw_land  = clear_val
    Interstitial%adjsfculw_ice   = clear_val
    Interstitial%adjsfculw_ocean = clear_val
    Interstitial%adjnirbmd       = clear_val
    Interstitial%adjnirbmu       = clear_val
    Interstitial%adjnirdfd       = clear_val
    Interstitial%adjnirdfu       = clear_val
    Interstitial%adjvisbmd       = clear_val
    Interstitial%adjvisbmu       = clear_val
    Interstitial%adjvisdfu       = clear_val
    Interstitial%adjvisdfd       = clear_val
    Interstitial%bexp1d          = clear_val
    Interstitial%cd              = clear_val
    Interstitial%cd_ice          = huge
    Interstitial%cd_land         = huge
    Interstitial%cd_ocean        = huge
    Interstitial%cdq             = clear_val
    Interstitial%cdq_ice         = huge
    Interstitial%cdq_land        = huge
    Interstitial%cdq_ocean       = huge
    Interstitial%chh_ice         = huge
    Interstitial%chh_land        = huge
    Interstitial%chh_ocean       = huge
    Interstitial%cld1d           = clear_val
    Interstitial%cldf            = clear_val
    Interstitial%clw             = clear_val
    Interstitial%clw(:,:,2)      = -999.9
    Interstitial%clx             = clear_val
    Interstitial%cmm_ice         = huge
    Interstitial%cmm_land        = huge
    Interstitial%cmm_ocean       = huge
    Interstitial%cnvc            = clear_val
    Interstitial%cnvw            = clear_val
    Interstitial%ctei_r          = clear_val
    Interstitial%ctei_rml        = clear_val
    Interstitial%cumabs          = clear_val
    Interstitial%dd_mf           = clear_val
    Interstitial%del             = clear_val
    Interstitial%del_gz          = clear_val
    Interstitial%dkt             = clear_val
    Interstitial%dlength         = clear_val
    Interstitial%dqdt            = clear_val
    Interstitial%dqsfc1          = clear_val
    Interstitial%drain           = clear_val
    Interstitial%dt_mf           = clear_val
    Interstitial%dtdt            = clear_val
    Interstitial%dtdtc           = clear_val
    Interstitial%dtsfc1          = clear_val
    Interstitial%dtzm            = clear_val
    Interstitial%dudt            = clear_val
    Interstitial%dusfcg          = clear_val
    Interstitial%dusfc1          = clear_val
    Interstitial%dvdftra         = clear_val
    Interstitial%dvdt            = clear_val
    Interstitial%dvsfcg          = clear_val
    Interstitial%dvsfc1          = clear_val
    Interstitial%elvmax          = clear_val
    Interstitial%ep1d            = clear_val
    Interstitial%ep1d_ice        = huge
    Interstitial%ep1d_land       = huge
    Interstitial%ep1d_ocean      = huge
    Interstitial%evap_ice        = huge
    Interstitial%evap_land       = huge
    Interstitial%evap_ocean      = huge
    Interstitial%evbs            = clear_val
    Interstitial%evcw            = clear_val
    Interstitial%ffhh_ice        = huge
    Interstitial%ffhh_land       = huge
    Interstitial%ffhh_ocean      = huge
    Interstitial%fh2             = clear_val
    Interstitial%fh2_ice         = huge
    Interstitial%fh2_land        = huge
    Interstitial%fh2_ocean       = huge
    Interstitial%flag_cice       = .false.
    Interstitial%flag_guess      = .false.
    Interstitial%flag_iter       = .true.
    Interstitial%ffmm_ice        = huge
    Interstitial%ffmm_land       = huge
    Interstitial%ffmm_ocean      = huge
    Interstitial%fm10            = clear_val
    Interstitial%fm10_ice        = huge
    Interstitial%fm10_land       = huge
    Interstitial%fm10_ocean      = huge
    Interstitial%frland          = clear_val
    Interstitial%fscav           = clear_val
    Interstitial%fswtr           = clear_val
    Interstitial%gabsbdlw        = clear_val
    Interstitial%gabsbdlw_ice    = clear_val
    Interstitial%gabsbdlw_land   = clear_val
    Interstitial%gabsbdlw_ocean  = clear_val
    Interstitial%gamma           = clear_val
    Interstitial%gamq            = clear_val
    Interstitial%gamt            = clear_val
    Interstitial%gflx            = clear_val
    Interstitial%gflx_ice        = zero
    Interstitial%gflx_land       = zero
    Interstitial%gflx_ocean      = zero
    Interstitial%gwdcu           = clear_val
    Interstitial%gwdcv           = clear_val
    Interstitial%hflx_ice        = huge
    Interstitial%hflx_land       = huge
    Interstitial%hflx_ocean      = huge
    Interstitial%dry             = .false.
    Interstitial%icy             = .false.
    Interstitial%lake            = .false.
    Interstitial%ocean           = .false.
    Interstitial%islmsk          = 0
    Interstitial%islmsk_cice     = 0
    Interstitial%wet             = .false.
    Interstitial%kbot            = Model%levs
    Interstitial%kcnv            = 0
    Interstitial%kinver          = Model%levs
    Interstitial%kpbl            = 0
    Interstitial%ktop            = 1
    Interstitial%oa4             = clear_val
    Interstitial%oc              = clear_val
    Interstitial%prcpmp          = clear_val
    Interstitial%prnum           = clear_val
    Interstitial%qss             = clear_val
    Interstitial%raincd          = clear_val
    Interstitial%raincs          = clear_val
    Interstitial%rainmcadj       = clear_val
    Interstitial%rainp           = clear_val
    Interstitial%rb              = clear_val
    Interstitial%rb_ice          = huge
    Interstitial%rb_land         = huge
    Interstitial%rb_ocean        = huge
    Interstitial%rhc             = clear_val
    Interstitial%runoff          = clear_val
    Interstitial%save_q          = clear_val
    Interstitial%save_t          = clear_val
    Interstitial%save_tcp        = clear_val
    Interstitial%save_u          = clear_val
    Interstitial%save_v          = clear_val
    Interstitial%sbsno           = clear_val
    Interstitial%semis_ice       = clear_val
    Interstitial%semis_land      = clear_val
    Interstitial%semis_ocean     = clear_val
    Interstitial%sigma           = clear_val
    Interstitial%sigmaf          = clear_val
    Interstitial%sigmafrac       = clear_val
    Interstitial%sigmatot        = clear_val
    Interstitial%slopetype       = 0
    Interstitial%snowc           = clear_val
    Interstitial%snohf           = clear_val
    Interstitial%snowmt          = clear_val
    Interstitial%soiltype        = 0
    Interstitial%stress          = clear_val
    Interstitial%stress_ice      = huge
    Interstitial%stress_land     = huge
    Interstitial%stress_ocean    = huge
    Interstitial%theta           = clear_val
    Interstitial%tice            = clear_val
    Interstitial%tprcp_ice       = huge
    Interstitial%tprcp_land      = huge
    Interstitial%tprcp_ocean     = huge
    Interstitial%trans           = clear_val
    Interstitial%tseal           = clear_val
    Interstitial%tsurf           = clear_val
    Interstitial%ud_mf           = clear_val
    Interstitial%ulwsfc_cice     = clear_val
    Interstitial%dusfc_cice      = clear_val
    Interstitial%dvsfc_cice      = clear_val
    Interstitial%dtsfc_cice      = clear_val
    Interstitial%dqsfc_cice      = clear_val
    Interstitial%vdftra          = clear_val
    Interstitial%vegf1d          = clear_val
    Interstitial%vegtype         = 0
    Interstitial%wcbmax          = clear_val
    Interstitial%weasd_ice       = huge
    Interstitial%weasd_land      = huge
    Interstitial%weasd_ocean     = huge
    Interstitial%wind            = huge
    Interstitial%work1           = clear_val
    Interstitial%work2           = clear_val
    Interstitial%work3           = clear_val
    Interstitial%xcosz           = clear_val
    Interstitial%xlai1d          = clear_val
    Interstitial%xmu             = clear_val
    Interstitial%z01d            = clear_val
    Interstitial%zt1d            = clear_val
! CIRES UGWP v0
    Interstitial%gw_dudt         = clear_val
    Interstitial%gw_dvdt         = clear_val
    Interstitial%gw_dtdt         = clear_val
    Interstitial%gw_kdis         = clear_val
    Interstitial%tau_mtb         = clear_val
    Interstitial%tau_ogw         = clear_val
    Interstitial%tau_tofd        = clear_val
    Interstitial%tau_ngw         = clear_val
    Interstitial%zmtb            = clear_val
    Interstitial%zlwb            = clear_val
    Interstitial%zogw            = clear_val
    Interstitial%dudt_mtb        = clear_val
    Interstitial%dudt_ogw        = clear_val
    Interstitial%dudt_tms        = clear_val
!
    ! Reset fields that are conditional on physics choices
    if (Model%imp_physics == Model%imp_physics_gfdl .or. Model%imp_physics == Model%imp_physics_thompson) then
       Interstitial%graupelmp = clear_val
       Interstitial%icemp     = clear_val
       Interstitial%rainmp    = clear_val
       Interstitial%snowmp    = clear_val
    else if (Model%imp_physics == Model%imp_physics_mg) then
       Interstitial%ncgl      = clear_val
       Interstitial%ncpr      = clear_val
       Interstitial%ncps      = clear_val
       Interstitial%qgl       = clear_val
       Interstitial%qrn       = clear_val
       Interstitial%qsnw      = clear_val
       Interstitial%qlcn      = clear_val
       Interstitial%qicn      = clear_val
       Interstitial%w_upi     = clear_val
       Interstitial%cf_upi    = clear_val
       Interstitial%cnv_mfd   = clear_val
       Interstitial%cnv_dqldt = clear_val
       Interstitial%clcn      = clear_val
       Interstitial%cnv_fice  = clear_val
       Interstitial%cnv_ndrop = clear_val
       Interstitial%cnv_nice  = clear_val
    end if
    if (Model%imp_physics == Model%imp_physics_fer_hires .and. Model%spec_adv) then
       Interstitial%f_ice     = clear_val
       Interstitial%f_rain    = clear_val
       Interstitial%f_rimef   = clear_val
       Interstitial%cwm       = clear_val
    end if
    if (Model%do_shoc) then
       Interstitial%qrn       = clear_val
       Interstitial%qsnw      = clear_val
       ! DH* updated version of shoc from May 22 2019 doesn't use qgl? remove?
       Interstitial%qgl       = clear_val
       ! *DH
       Interstitial%ncpi      = clear_val
       Interstitial%ncpl      = clear_val
    end if
    if (Model%lsm == Model%lsm_noahmp) then
       Interstitial%t2mmp     = clear_val
       Interstitial%q2mp      = clear_val
    end if
    !
    ! Set flag for resetting maximum hourly output fields
    Interstitial%reset = mod(Model%kdt-1, nint(Model%avg_max_length/Model%dtp)) == 0
    !
  end subroutine interstitial_phys_reset

  subroutine interstitial_print(Interstitial, Model, mpirank, omprank, blkno)
    !
    implicit none
    !
    class(GFS_interstitial_type) :: Interstitial
    type(GFS_control_type), intent(in) :: Model
    integer, intent(in) :: mpirank, omprank, blkno
    !
    ! Print static variables
    write (0,'(a,3i6)') 'Interstitial_print for mpirank, omprank, blkno: ', mpirank, omprank, blkno
    write (0,*) 'Interstitial_print: values that do not change'
    write (0,*) 'Interstitial%h2o_coeff         = ', Interstitial%h2o_coeff
    write (0,*) 'sum(Interstitial%h2o_pres)     = ', sum(Interstitial%h2o_pres)
    write (0,*) 'Interstitial%ipr               = ', Interstitial%ipr
    write (0,*) 'Interstitial%itc               = ', Interstitial%itc
    write (0,*) 'Interstitial%latidxprnt        = ', Interstitial%latidxprnt
    write (0,*) 'Interstitial%levi              = ', Interstitial%levi
    write (0,*) 'Interstitial%levh2o            = ', Interstitial%levh2o
    write (0,*) 'Interstitial%levozp            = ', Interstitial%levozp
    write (0,*) 'Interstitial%lmk               = ', Interstitial%lmk
    write (0,*) 'Interstitial%lmp               = ', Interstitial%lmp
    write (0,*) 'Interstitial%nbdlw             = ', Interstitial%nbdlw
    write (0,*) 'Interstitial%nbdsw             = ', Interstitial%nbdsw
    write (0,*) 'Interstitial%nf_aelw           = ', Interstitial%nf_aelw
    write (0,*) 'Interstitial%nf_aesw           = ', Interstitial%nf_aesw
    write (0,*) 'Interstitial%nsamftrac         = ', Interstitial%nsamftrac
    write (0,*) 'Interstitial%nscav             = ', Interstitial%nscav
    write (0,*) 'Interstitial%nspc1             = ', Interstitial%nspc1
    write (0,*) 'Interstitial%ntiwx             = ', Interstitial%ntiwx
    write (0,*) 'Interstitial%nvdiff            = ', Interstitial%nvdiff
    write (0,*) 'Interstitial%oz_coeff          = ', Interstitial%oz_coeff
    write (0,*) 'sum(Interstitial%oz_pres)      = ', sum(Interstitial%oz_pres)
    write (0,*) 'Interstitial%phys_hydrostatic  = ', Interstitial%phys_hydrostatic
    write (0,*) 'Interstitial%skip_macro        = ', Interstitial%skip_macro
    write (0,*) 'Interstitial%trans_aero        = ', Interstitial%trans_aero
    ! Print all other variables
    write (0,*) 'Interstitial_print: values that change'
    write (0,*) 'sum(Interstitial%adjsfculw_land  ) = ', sum(Interstitial%adjsfculw_land  )
    write (0,*) 'sum(Interstitial%adjsfculw_ice   ) = ', sum(Interstitial%adjsfculw_ice   )
    write (0,*) 'sum(Interstitial%adjsfculw_ocean ) = ', sum(Interstitial%adjsfculw_ocean )
    write (0,*) 'sum(Interstitial%adjnirbmd       ) = ', sum(Interstitial%adjnirbmd       )
    write (0,*) 'sum(Interstitial%adjnirbmu       ) = ', sum(Interstitial%adjnirbmu       )
    write (0,*) 'sum(Interstitial%adjnirdfd       ) = ', sum(Interstitial%adjnirdfd       )
    write (0,*) 'sum(Interstitial%adjnirdfu       ) = ', sum(Interstitial%adjnirdfu       )
    write (0,*) 'sum(Interstitial%adjvisbmd       ) = ', sum(Interstitial%adjvisbmd       )
    write (0,*) 'sum(Interstitial%adjvisbmu       ) = ', sum(Interstitial%adjvisbmu       )
    write (0,*) 'sum(Interstitial%adjvisdfu       ) = ', sum(Interstitial%adjvisdfu       )
    write (0,*) 'sum(Interstitial%adjvisdfd       ) = ', sum(Interstitial%adjvisdfd       )
    write (0,*) 'sum(Interstitial%aerodp          ) = ', sum(Interstitial%aerodp          )
    write (0,*) 'sum(Interstitial%alb1d           ) = ', sum(Interstitial%alb1d           )
    write (0,*) 'sum(Interstitial%bexp1d          ) = ', sum(Interstitial%bexp1d          )
    write (0,*) 'sum(Interstitial%cd              ) = ', sum(Interstitial%cd              )
    write (0,*) 'sum(Interstitial%cd_ice          ) = ', sum(Interstitial%cd_ice          )
    write (0,*) 'sum(Interstitial%cd_land         ) = ', sum(Interstitial%cd_land         )
    write (0,*) 'sum(Interstitial%cd_ocean        ) = ', sum(Interstitial%cd_ocean        )
    write (0,*) 'sum(Interstitial%cdq             ) = ', sum(Interstitial%cdq             )
    write (0,*) 'sum(Interstitial%cdq_ice         ) = ', sum(Interstitial%cdq_ice         )
    write (0,*) 'sum(Interstitial%cdq_land        ) = ', sum(Interstitial%cdq_land        )
    write (0,*) 'sum(Interstitial%cdq_ocean       ) = ', sum(Interstitial%cdq_ocean       )
    write (0,*) 'sum(Interstitial%chh_ice         ) = ', sum(Interstitial%chh_ice         )
    write (0,*) 'sum(Interstitial%chh_land        ) = ', sum(Interstitial%chh_land        )
    write (0,*) 'sum(Interstitial%chh_ocean       ) = ', sum(Interstitial%chh_ocean       )
    write (0,*) 'sum(Interstitial%cldf            ) = ', sum(Interstitial%cldf            )
    write (0,*) 'sum(Interstitial%cldsa           ) = ', sum(Interstitial%cldsa           )
    write (0,*) 'sum(Interstitial%cldtaulw        ) = ', sum(Interstitial%cldtaulw        )
    write (0,*) 'sum(Interstitial%cldtausw        ) = ', sum(Interstitial%cldtausw        )
    write (0,*) 'sum(Interstitial%cld1d           ) = ', sum(Interstitial%cld1d           )
    write (0,*) 'sum(Interstitial%clw             ) = ', sum(Interstitial%clw             )
    write (0,*) 'sum(Interstitial%clx             ) = ', sum(Interstitial%clx             )
    write (0,*) 'sum(Interstitial%clouds          ) = ', sum(Interstitial%clouds          )
    write (0,*) 'sum(Interstitial%cmm_ice         ) = ', sum(Interstitial%cmm_ice         )
    write (0,*) 'sum(Interstitial%cmm_land        ) = ', sum(Interstitial%cmm_land        )
    write (0,*) 'sum(Interstitial%cmm_ocean       ) = ', sum(Interstitial%cmm_ocean       )
    write (0,*) 'sum(Interstitial%cnvc            ) = ', sum(Interstitial%cnvc            )
    write (0,*) 'sum(Interstitial%cnvw            ) = ', sum(Interstitial%cnvw            )
    write (0,*) 'sum(Interstitial%ctei_r          ) = ', sum(Interstitial%ctei_r          )
    write (0,*) 'sum(Interstitial%ctei_rml        ) = ', sum(Interstitial%ctei_rml        )
    write (0,*) 'sum(Interstitial%cumabs          ) = ', sum(Interstitial%cumabs          )
    write (0,*) 'sum(Interstitial%dd_mf           ) = ', sum(Interstitial%dd_mf           )
    write (0,*) 'sum(Interstitial%de_lgth         ) = ', sum(Interstitial%de_lgth         )
    write (0,*) 'sum(Interstitial%del             ) = ', sum(Interstitial%del             )
    write (0,*) 'sum(Interstitial%del_gz          ) = ', sum(Interstitial%del_gz          )
    write (0,*) 'sum(Interstitial%delr            ) = ', sum(Interstitial%delr            )
    write (0,*) 'sum(Interstitial%dkt             ) = ', sum(Interstitial%dkt             )
    write (0,*) 'sum(Interstitial%dlength         ) = ', sum(Interstitial%dlength         )
    write (0,*) 'sum(Interstitial%dqdt            ) = ', sum(Interstitial%dqdt            )
    write (0,*) 'sum(Interstitial%dqsfc1          ) = ', sum(Interstitial%dqsfc1          )
    write (0,*) 'sum(Interstitial%drain           ) = ', sum(Interstitial%drain           )
    write (0,*) 'sum(Interstitial%dtdt            ) = ', sum(Interstitial%dtdt            )
    write (0,*) 'sum(Interstitial%dtdtc           ) = ', sum(Interstitial%dtdtc           )
    write (0,*) 'sum(Interstitial%dtsfc1          ) = ', sum(Interstitial%dtsfc1          )
    write (0,*) 'sum(Interstitial%dtzm            ) = ', sum(Interstitial%dtzm            )
    write (0,*) 'sum(Interstitial%dt_mf           ) = ', sum(Interstitial%dt_mf           )
    write (0,*) 'sum(Interstitial%dudt            ) = ', sum(Interstitial%dudt            )
    write (0,*) 'sum(Interstitial%dusfcg          ) = ', sum(Interstitial%dusfcg          )
    write (0,*) 'sum(Interstitial%dusfc1          ) = ', sum(Interstitial%dusfc1          )
    write (0,*) 'sum(Interstitial%dvdftra         ) = ', sum(Interstitial%dvdftra         )
    write (0,*) 'sum(Interstitial%dvdt            ) = ', sum(Interstitial%dvdt            )
    write (0,*) 'sum(Interstitial%dvsfcg          ) = ', sum(Interstitial%dvsfcg          )
    write (0,*) 'sum(Interstitial%dvsfc1          ) = ', sum(Interstitial%dvsfc1          )
    write (0,*) 'sum(Interstitial%dzlyr           ) = ', sum(Interstitial%dzlyr           )
    write (0,*) 'sum(Interstitial%elvmax          ) = ', sum(Interstitial%elvmax          )
    write (0,*) 'sum(Interstitial%ep1d            ) = ', sum(Interstitial%ep1d            )
    write (0,*) 'sum(Interstitial%ep1d_ice        ) = ', sum(Interstitial%ep1d_ice        )
    write (0,*) 'sum(Interstitial%ep1d_land       ) = ', sum(Interstitial%ep1d_land       )
    write (0,*) 'sum(Interstitial%ep1d_ocean      ) = ', sum(Interstitial%ep1d_ocean      )
    write (0,*) 'sum(Interstitial%evap_ice        ) = ', sum(Interstitial%evap_ice        )
    write (0,*) 'sum(Interstitial%evap_land       ) = ', sum(Interstitial%evap_land       )
    write (0,*) 'sum(Interstitial%evap_ocean      ) = ', sum(Interstitial%evap_ocean      )
    write (0,*) 'sum(Interstitial%evbs            ) = ', sum(Interstitial%evbs            )
    write (0,*) 'sum(Interstitial%evcw            ) = ', sum(Interstitial%evcw            )
    write (0,*) 'sum(Interstitial%faerlw          ) = ', sum(Interstitial%faerlw          )
    write (0,*) 'sum(Interstitial%faersw          ) = ', sum(Interstitial%faersw          )
    write (0,*) 'sum(Interstitial%ffhh_ice        ) = ', sum(Interstitial%ffhh_ice        )
    write (0,*) 'sum(Interstitial%ffhh_land       ) = ', sum(Interstitial%ffhh_land       )
    write (0,*) 'sum(Interstitial%ffhh_ocean      ) = ', sum(Interstitial%ffhh_ocean      )
    write (0,*) 'sum(Interstitial%fh2             ) = ', sum(Interstitial%fh2             )
    write (0,*) 'sum(Interstitial%fh2_ice         ) = ', sum(Interstitial%fh2_ice         )
    write (0,*) 'sum(Interstitial%fh2_land        ) = ', sum(Interstitial%fh2_land        )
    write (0,*) 'sum(Interstitial%fh2_ocean       ) = ', sum(Interstitial%fh2_ocean       )
    write (0,*) 'Interstitial%flag_cice(1)          = ', Interstitial%flag_cice(1)
    write (0,*) 'Interstitial%flag_guess(1)         = ', Interstitial%flag_guess(1)
    write (0,*) 'Interstitial%flag_iter(1)          = ', Interstitial%flag_iter(1)
    write (0,*) 'sum(Interstitial%ffmm_ice        ) = ', sum(Interstitial%ffmm_ice        )
    write (0,*) 'sum(Interstitial%ffmm_land       ) = ', sum(Interstitial%ffmm_land       )
    write (0,*) 'sum(Interstitial%ffmm_ocean      ) = ', sum(Interstitial%ffmm_ocean      )
    write (0,*) 'sum(Interstitial%fm10            ) = ', sum(Interstitial%fm10            )
    write (0,*) 'sum(Interstitial%fm10_ice        ) = ', sum(Interstitial%fm10_ice        )
    write (0,*) 'sum(Interstitial%fm10_land       ) = ', sum(Interstitial%fm10_land       )
    write (0,*) 'sum(Interstitial%fm10_ocean      ) = ', sum(Interstitial%fm10_ocean      )
    write (0,*) 'Interstitial%frain                 = ', Interstitial%frain
    write (0,*) 'sum(Interstitial%frland          ) = ', sum(Interstitial%frland          )
    write (0,*) 'sum(Interstitial%fscav           ) = ', sum(Interstitial%fscav           )
    write (0,*) 'sum(Interstitial%fswtr           ) = ', sum(Interstitial%fswtr           )
    write (0,*) 'sum(Interstitial%gabsbdlw        ) = ', sum(Interstitial%gabsbdlw        )
    write (0,*) 'sum(Interstitial%gabsbdlw_ice    ) = ', sum(Interstitial%gabsbdlw_ice    )
    write (0,*) 'sum(Interstitial%gabsbdlw_land   ) = ', sum(Interstitial%gabsbdlw_land   )
    write (0,*) 'sum(Interstitial%gabsbdlw_ocean  ) = ', sum(Interstitial%gabsbdlw_ocean  )
    write (0,*) 'sum(Interstitial%gamma           ) = ', sum(Interstitial%gamma           )
    write (0,*) 'sum(Interstitial%gamq            ) = ', sum(Interstitial%gamq            )
    write (0,*) 'sum(Interstitial%gamt            ) = ', sum(Interstitial%gamt            )
    write (0,*) 'sum(Interstitial%gasvmr          ) = ', sum(Interstitial%gasvmr          )
    write (0,*) 'sum(Interstitial%gflx            ) = ', sum(Interstitial%gflx            )
    write (0,*) 'sum(Interstitial%gflx_ice        ) = ', sum(Interstitial%gflx_ice        )
    write (0,*) 'sum(Interstitial%gflx_land       ) = ', sum(Interstitial%gflx_land       )
    write (0,*) 'sum(Interstitial%gflx_ocean      ) = ', sum(Interstitial%gflx_ocean      )
    write (0,*) 'sum(Interstitial%gwdcu           ) = ', sum(Interstitial%gwdcu           )
    write (0,*) 'sum(Interstitial%gwdcv           ) = ', sum(Interstitial%gwdcv           )
    write (0,*) 'sum(Interstitial%hflx_ice        ) = ', sum(Interstitial%hflx_ice        )
    write (0,*) 'sum(Interstitial%hflx_land       ) = ', sum(Interstitial%hflx_land       )
    write (0,*) 'sum(Interstitial%hflx_ocean      ) = ', sum(Interstitial%hflx_ocean      )
    write (0,*) 'Interstitial%dry(:)==.true.        = ', count(Interstitial%dry(:)        )
    write (0,*) 'sum(Interstitial%idxday          ) = ', sum(Interstitial%idxday          )
    write (0,*) 'Interstitial%icy(:)==.true.        = ', count(Interstitial%icy(:)        )
    write (0,*) 'Interstitial%lake(:)==.true.       = ', count(Interstitial%lake(:)       )
    write (0,*) 'Interstitial%ocean(:)==.true.      = ', count(Interstitial%ocean(:)      )
    write (0,*) 'sum(Interstitial%islmsk          ) = ', sum(Interstitial%islmsk          )
    write (0,*) 'sum(Interstitial%islmsk_cice     ) = ', sum(Interstitial%islmsk_cice     )
    write (0,*) 'Interstitial%wet(:)==.true.        = ', count(Interstitial%wet(:)        )
    write (0,*) 'Interstitial%kb                    = ', Interstitial%kb
    write (0,*) 'sum(Interstitial%kbot            ) = ', sum(Interstitial%kbot            )
    write (0,*) 'sum(Interstitial%kcnv            ) = ', sum(Interstitial%kcnv            )
    write (0,*) 'Interstitial%kd                    = ', Interstitial%kd
    write (0,*) 'sum(Interstitial%kinver          ) = ', sum(Interstitial%kinver          )
    write (0,*) 'sum(Interstitial%kpbl            ) = ', sum(Interstitial%kpbl            )
    write (0,*) 'Interstitial%kt                    = ', Interstitial%kt
    write (0,*) 'sum(Interstitial%ktop            ) = ', sum(Interstitial%ktop            )
    write (0,*) 'sum(Interstitial%mbota           ) = ', sum(Interstitial%mbota           )
    write (0,*) 'sum(Interstitial%mtopa           ) = ', sum(Interstitial%mtopa           )
    write (0,*) 'Interstitial%nday                  = ', Interstitial%nday
    write (0,*) 'sum(Interstitial%oa4             ) = ', sum(Interstitial%oa4             )
    write (0,*) 'sum(Interstitial%oc              ) = ', sum(Interstitial%oc              )
    write (0,*) 'sum(Interstitial%olyr            ) = ', sum(Interstitial%olyr            )
    write (0,*) 'sum(Interstitial%plvl            ) = ', sum(Interstitial%plvl            )
    write (0,*) 'sum(Interstitial%plyr            ) = ', sum(Interstitial%plyr            )
    write (0,*) 'sum(Interstitial%prcpmp          ) = ', sum(Interstitial%prcpmp          )
    write (0,*) 'sum(Interstitial%prnum           ) = ', sum(Interstitial%prnum           )
    write (0,*) 'sum(Interstitial%qlyr            ) = ', sum(Interstitial%qlyr            )
    write (0,*) 'sum(Interstitial%qss             ) = ', sum(Interstitial%qss             )
    write (0,*) 'Interstitial%raddt                 = ', Interstitial%raddt
    write (0,*) 'sum(Interstitial%raincd          ) = ', sum(Interstitial%raincd          )
    write (0,*) 'sum(Interstitial%raincs          ) = ', sum(Interstitial%raincs          )
    write (0,*) 'sum(Interstitial%rainmcadj       ) = ', sum(Interstitial%rainmcadj       )
    write (0,*) 'sum(Interstitial%rainp           ) = ', sum(Interstitial%rainp           )
    write (0,*) 'sum(Interstitial%rb              ) = ', sum(Interstitial%rb              )
    write (0,*) 'sum(Interstitial%rb_ice          ) = ', sum(Interstitial%rb_ice          )
    write (0,*) 'sum(Interstitial%rb_land         ) = ', sum(Interstitial%rb_land         )
    write (0,*) 'sum(Interstitial%rb_ocean        ) = ', sum(Interstitial%rb_ocean        )
    write (0,*) 'Interstitial%reset                 = ', Interstitial%reset
    write (0,*) 'sum(Interstitial%rhc             ) = ', sum(Interstitial%rhc             )
    write (0,*) 'sum(Interstitial%runoff          ) = ', sum(Interstitial%runoff          )
    write (0,*) 'sum(Interstitial%save_q          ) = ', sum(Interstitial%save_q          )
    write (0,*) 'sum(Interstitial%save_t          ) = ', sum(Interstitial%save_t          )
    write (0,*) 'sum(Interstitial%save_tcp        ) = ', sum(Interstitial%save_tcp        )
    write (0,*) 'sum(Interstitial%save_u          ) = ', sum(Interstitial%save_u          )
    write (0,*) 'sum(Interstitial%save_v          ) = ', sum(Interstitial%save_v          )
    write (0,*) 'sum(Interstitial%sbsno           ) = ', sum(Interstitial%sbsno           )
    write (0,*) 'sum(Interstitial%scmpsw%uvbfc    ) = ', sum(Interstitial%scmpsw%uvbfc    )
    write (0,*) 'sum(Interstitial%scmpsw%uvbf0    ) = ', sum(Interstitial%scmpsw%uvbf0    )
    write (0,*) 'sum(Interstitial%scmpsw%nirbm    ) = ', sum(Interstitial%scmpsw%nirbm    )
    write (0,*) 'sum(Interstitial%scmpsw%nirdf    ) = ', sum(Interstitial%scmpsw%nirdf    )
    write (0,*) 'sum(Interstitial%scmpsw%visbm    ) = ', sum(Interstitial%scmpsw%visbm    )
    write (0,*) 'sum(Interstitial%scmpsw%visdf    ) = ', sum(Interstitial%scmpsw%visdf    )
    write (0,*) 'sum(Interstitial%semis_ice       ) = ', sum(Interstitial%semis_ice       )
    write (0,*) 'sum(Interstitial%semis_land      ) = ', sum(Interstitial%semis_land      )
    write (0,*) 'sum(Interstitial%semis_ocean     ) = ', sum(Interstitial%semis_ocean     )
    write (0,*) 'sum(Interstitial%sfcalb          ) = ', sum(Interstitial%sfcalb          )
    write (0,*) 'sum(Interstitial%sigma           ) = ', sum(Interstitial%sigma           )
    write (0,*) 'sum(Interstitial%sigmaf          ) = ', sum(Interstitial%sigmaf          )
    write (0,*) 'sum(Interstitial%sigmafrac       ) = ', sum(Interstitial%sigmafrac       )
    write (0,*) 'sum(Interstitial%sigmatot        ) = ', sum(Interstitial%sigmatot        )
    write (0,*) 'sum(Interstitial%slopetype       ) = ', sum(Interstitial%slopetype       )
    write (0,*) 'sum(Interstitial%snowc           ) = ', sum(Interstitial%snowc           )
    write (0,*) 'sum(Interstitial%snohf           ) = ', sum(Interstitial%snohf           )
    write (0,*) 'sum(Interstitial%snowmt          ) = ', sum(Interstitial%snowmt          )
    write (0,*) 'sum(Interstitial%soiltype        ) = ', sum(Interstitial%soiltype        )
    write (0,*) 'sum(Interstitial%stress          ) = ', sum(Interstitial%stress          )
    write (0,*) 'sum(Interstitial%stress_ice      ) = ', sum(Interstitial%stress_ice      )
    write (0,*) 'sum(Interstitial%stress_land     ) = ', sum(Interstitial%stress_land     )
    write (0,*) 'sum(Interstitial%stress_ocean    ) = ', sum(Interstitial%stress_ocean    )
    write (0,*) 'sum(Interstitial%theta           ) = ', sum(Interstitial%theta           )
    write (0,*) 'sum(Interstitial%tice            ) = ', sum(Interstitial%tice            )
    write (0,*) 'sum(Interstitial%tlvl            ) = ', sum(Interstitial%tlvl            )
    write (0,*) 'sum(Interstitial%tlyr            ) = ', sum(Interstitial%tlyr            )
    write (0,*) 'sum(Interstitial%tprcp_ice       ) = ', sum(Interstitial%tprcp_ice       )
    write (0,*) 'sum(Interstitial%tprcp_land      ) = ', sum(Interstitial%tprcp_land      )
    write (0,*) 'sum(Interstitial%tprcp_ocean     ) = ', sum(Interstitial%tprcp_ocean     )
    write (0,*) 'sum(Interstitial%trans           ) = ', sum(Interstitial%trans           )
    write (0,*) 'sum(Interstitial%tseal           ) = ', sum(Interstitial%tseal           )
    write (0,*) 'sum(Interstitial%tsfa            ) = ', sum(Interstitial%tsfa            )
    write (0,*) 'sum(Interstitial%tsfg            ) = ', sum(Interstitial%tsfg            )
    write (0,*) 'sum(Interstitial%tsurf           ) = ', sum(Interstitial%tsurf           )
    write (0,*) 'sum(Interstitial%ud_mf           ) = ', sum(Interstitial%ud_mf           )
    write (0,*) 'sum(Interstitial%ulwsfc_cice     ) = ', sum(Interstitial%ulwsfc_cice     )
    write (0,*) 'sum(Interstitial%dusfc_cice      ) = ', sum(Interstitial%dusfc_cice      )
    write (0,*) 'sum(Interstitial%dvsfc_cice      ) = ', sum(Interstitial%dvsfc_cice      )
    write (0,*) 'sum(Interstitial%dtsfc_cice      ) = ', sum(Interstitial%dtsfc_cice      )
    write (0,*) 'sum(Interstitial%dqsfc_cice      ) = ', sum(Interstitial%dqsfc_cice      )
    write (0,*) 'sum(Interstitial%vdftra          ) = ', sum(Interstitial%vdftra          )
    write (0,*) 'sum(Interstitial%vegf1d          ) = ', sum(Interstitial%vegf1d          )
    write (0,*) 'sum(Interstitial%vegtype         ) = ', sum(Interstitial%vegtype         )
    write (0,*) 'sum(Interstitial%wcbmax          ) = ', sum(Interstitial%wcbmax          )
    write (0,*) 'sum(Interstitial%weasd_ice       ) = ', sum(Interstitial%weasd_ice       )
    write (0,*) 'sum(Interstitial%weasd_land      ) = ', sum(Interstitial%weasd_land      )
    write (0,*) 'sum(Interstitial%weasd_ocean     ) = ', sum(Interstitial%weasd_ocean     )
    write (0,*) 'sum(Interstitial%wind            ) = ', sum(Interstitial%wind            )
    write (0,*) 'sum(Interstitial%work1           ) = ', sum(Interstitial%work1           )
    write (0,*) 'sum(Interstitial%work2           ) = ', sum(Interstitial%work2           )
    write (0,*) 'sum(Interstitial%work3           ) = ', sum(Interstitial%work3           )
    write (0,*) 'sum(Interstitial%xcosz           ) = ', sum(Interstitial%xcosz           )
    write (0,*) 'sum(Interstitial%xlai1d          ) = ', sum(Interstitial%xlai1d          )
    write (0,*) 'sum(Interstitial%xmu             ) = ', sum(Interstitial%xmu             )
    write (0,*) 'sum(Interstitial%z01d            ) = ', sum(Interstitial%z01d            )
    write (0,*) 'sum(Interstitial%zt1d            ) = ', sum(Interstitial%zt1d            )
! CIRES UGWP v0
    write (0,*) 'sum(Interstitial%gw_dudt         ) = ', sum(Interstitial%gw_dudt         )
    write (0,*) 'sum(Interstitial%gw_dvdt         ) = ', sum(Interstitial%gw_dvdt         )
    write (0,*) 'sum(Interstitial%gw_dtdt         ) = ', sum(Interstitial%gw_dtdt         )
    write (0,*) 'sum(Interstitial%gw_kdis         ) = ', sum(Interstitial%gw_kdis         )
    write (0,*) 'sum(Interstitial%tau_mtb         ) = ', sum(Interstitial%tau_mtb         )
    write (0,*) 'sum(Interstitial%tau_ogw         ) = ', sum(Interstitial%tau_ogw         )
    write (0,*) 'sum(Interstitial%tau_tofd        ) = ', sum(Interstitial%tau_tofd        )
    write (0,*) 'sum(Interstitial%tau_ngw         ) = ', sum(Interstitial%tau_ngw         )
    write (0,*) 'sum(Interstitial%zmtb            ) = ', sum(Interstitial%zmtb            )
    write (0,*) 'sum(Interstitial%zlwb            ) = ', sum(Interstitial%zlwb            )
    write (0,*) 'sum(Interstitial%zogw            ) = ', sum(Interstitial%zogw            )
    write (0,*) 'sum(Interstitial%dudt_mtb        ) = ', sum(Interstitial%dudt_mtb        )
    write (0,*) 'sum(Interstitial%dudt_ogw        ) = ', sum(Interstitial%dudt_ogw        )
    write (0,*) 'sum(Interstitial%dudt_tms        ) = ', sum(Interstitial%dudt_tms        )
!
    ! Print arrays that are conditional on physics choices
    if (Model%imp_physics == Model%imp_physics_gfdl .or. Model%imp_physics == Model%imp_physics_thompson) then
       write (0,*) 'Interstitial_print: values specific to GFDL/Thompson microphysics'
       write (0,*) 'sum(Interstitial%graupelmp    ) = ', sum(Interstitial%graupelmp       )
       write (0,*) 'sum(Interstitial%icemp        ) = ', sum(Interstitial%icemp           )
       write (0,*) 'sum(Interstitial%rainmp       ) = ', sum(Interstitial%rainmp          )
       write (0,*) 'sum(Interstitial%snowmp       ) = ', sum(Interstitial%snowmp          )
    !F-A scheme
    else if (Model%imp_physics == Model%imp_physics_fer_hires) then
       write (0,*) 'Interstitial_print: values specific to F-A microphysics'
       write (0,*) 'sum(Interstitial%f_ice        ) = ', sum(Interstitial%f_ice           )
       write (0,*) 'sum(Interstitial%f_rain       ) = ', sum(Interstitial%f_rain          )
       write (0,*) 'sum(Interstitial%f_rimef      ) = ', sum(Interstitial%f_rimef         )
       write (0,*) 'sum(Interstitial%cwm          ) = ', sum(Interstitial%cwm             )
    else if (Model%imp_physics == Model%imp_physics_mg) then
       write (0,*) 'Interstitial_print: values specific to MG microphysics'
       write (0,*) 'sum(Interstitial%ncgl         ) = ', sum(Interstitial%ncgl            )
       write (0,*) 'sum(Interstitial%ncpr         ) = ', sum(Interstitial%ncpr            )
       write (0,*) 'sum(Interstitial%ncps         ) = ', sum(Interstitial%ncps            )
       write (0,*) 'sum(Interstitial%qgl          ) = ', sum(Interstitial%qgl             )
       write (0,*) 'sum(Interstitial%qrn          ) = ', sum(Interstitial%qrn             )
       write (0,*) 'sum(Interstitial%qsnw         ) = ', sum(Interstitial%qsnw            )
       write (0,*) 'sum(Interstitial%qlcn         ) = ', sum(Interstitial%qlcn            )
       write (0,*) 'sum(Interstitial%qicn         ) = ', sum(Interstitial%qicn            )
       write (0,*) 'sum(Interstitial%w_upi        ) = ', sum(Interstitial%w_upi           )
       write (0,*) 'sum(Interstitial%cf_upi       ) = ', sum(Interstitial%cf_upi          )
       write (0,*) 'sum(Interstitial%cnv_mfd      ) = ', sum(Interstitial%cnv_mfd         )
       write (0,*) 'sum(Interstitial%cnv_dqldt    ) = ', sum(Interstitial%cnv_dqldt       )
       write (0,*) 'sum(Interstitial%clcn         ) = ', sum(Interstitial%clcn            )
       write (0,*) 'sum(Interstitial%cnv_fice     ) = ', sum(Interstitial%cnv_fice        )
       write (0,*) 'sum(Interstitial%cnv_ndrop    ) = ', sum(Interstitial%cnv_ndrop       )
       write (0,*) 'sum(Interstitial%cnv_nice     ) = ', sum(Interstitial%cnv_nice        )
    end if
    if (Model%do_shoc) then
       write (0,*) 'Interstitial_print: values specific to SHOC'
       write (0,*) 'sum(Interstitial%ncgl         ) = ', sum(Interstitial%ncgl            )
       write (0,*) 'sum(Interstitial%qrn          ) = ', sum(Interstitial%qrn             )
       write (0,*) 'sum(Interstitial%qsnw         ) = ', sum(Interstitial%qsnw            )
       write (0,*) 'sum(Interstitial%qgl          ) = ', sum(Interstitial%qgl             )
       write (0,*) 'sum(Interstitial%ncpi         ) = ', sum(Interstitial%ncpi            )
       write (0,*) 'sum(Interstitial%ncpl         ) = ', sum(Interstitial%ncpl            )
    end if
    if (Model%lsm == Model%lsm_noahmp) then
       write (0,*) 'sum(Interstitial%t2mmp        ) = ', sum(Interstitial%t2mmp           )
       write (0,*) 'sum(Interstitial%q2mp         ) = ', sum(Interstitial%q2mp            )
    end if
    write (0,*) 'Interstitial_print: end'
    !
  end subroutine interstitial_print
#endif

end module GFS_typedefs
