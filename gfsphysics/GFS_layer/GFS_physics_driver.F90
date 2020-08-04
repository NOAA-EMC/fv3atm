module module_physics_driver

  use machine,               only: kind_phys
  use physcons,              only: con_cp, con_fvirt, con_g, con_rd,    &
                                   con_rv, con_hvap, con_hfus,          &
                                   con_rerth, con_pi, rhc_max, dxmin,   &
                                   dxinv, pa2mb, rlapse, con_eps,       &
                                   con_epsm1, PQ0, A2A, A3, A4, RHmin,  &
                                   tgice => con_tice

  use cs_conv,               only: cs_convr
  use ozne_def,              only: levozp,  oz_coeff, oz_pres
  use h2o_def,               only: levh2o, h2o_coeff, h2o_pres
  use gfs_fv3_needs,         only: get_prs_fv3, get_phi_fv3
  use module_nst_water_prop, only: get_dtzm_2d
  use GFS_typedefs,          only: GFS_statein_type, GFS_stateout_type, &
                                   GFS_sfcprop_type, GFS_coupling_type, &
                                   GFS_control_type, GFS_grid_type,     &
                                   GFS_tbd_type,     GFS_cldprop_type,  &
!                                  GFS_radtend_type, GFS_diag_type
                                   GFS_radtend_type, GFS_diag_type, huge
  use gfdl_cloud_microphys_mod, only: gfdl_cloud_microphys_driver,      &
                                      cloud_diagnosis
  use module_mp_thompson,    only: mp_gt_driver
  use module_mp_wsm6,        only: wsm6
  use funcphys,              only: ftdp
  use surface_perturbation,  only: cdfnor

  use module_sfc_diff,       only: sfc_diff
  use module_sfc_ocean,      only: sfc_ocean
  use module_sfc_drv,        only: sfc_drv
  use module_sfc_sice,       only: sfc_sice
  use module_sfc_cice,       only: sfc_cice
  use module_sfc_nst,        only: sfc_nst
  use module_sfc_diag,       only: sfc_diag
!
!vay-2018
!
  use cires_ugwp_module,     only:  cires_ugwp_driver, knob_ugwp_version
!

  implicit none


  !--- CONSTANT PARAMETERS
  real(kind=kind_phys), parameter :: hocp    = con_hvap/con_cp
  real(kind=kind_phys), parameter :: epsln   = 1.0e-10_kind_phys
  real(kind=kind_phys), parameter :: qmin    = 1.0e-8_kind_phys
  real(kind=kind_phys), parameter :: qsmall  = 1.0e-20_kind_phys
  real(kind=kind_phys), parameter :: rainmin = 1.0e-13_kind_phys
  real(kind=kind_phys), parameter :: p850    = 85000.0_kind_phys
  real(kind=kind_phys), parameter :: epsq    = 1.0e-20_kind_phys
  real(kind=kind_phys), parameter :: hsub    = con_hvap+con_hfus
  real(kind=kind_phys), parameter :: czmin   = 0.0001_kind_phys      ! cos(89.994)
  real(kind=kind_phys), parameter :: zero    = 0.0_kind_phys, one   = 1.0_kind_phys,  &
                                     half    = 0.5_kind_phys, onebg = one/con_g
  real(kind=kind_phys), parameter :: albdf   = 0.06_kind_phys
  real(kind=kind_phys), parameter :: tf=258.16_kind_phys, tcr=273.16_kind_phys, tcrf=one/(tcr-tf)
  real(kind=kind_phys), parameter :: con_p001= 0.001_kind_phys
  real(kind=kind_phys), parameter :: con_day = 86400.0_kind_phys
  real(kind=kind_phys), parameter :: rad2dg  = 180.0_kind_phys/con_pi
  real(kind=kind_phys), parameter :: omz1    = 2.0_kind_phys
! real(kind=kind_phys), parameter :: huge    = 0.0_kind_phys

!> GFS Physics Implementation Layer
!> @brief Layer that invokes individual GFS physics routines
!> @{
!at tune step===========================================================!
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call GFS_physics_driver                                            !
!                                                                       !
!  ---  interface variables                                             !
!     type(GFS_control_type),         intent(in)    :: Model            !
!     type(GFS_statein_type),         intent(inout) :: Statein          !
!     type(GFS_stateout_type),        intent(inout) :: Stateout         !
!     type(GFS_sfcprop_type),         intent(inout) :: Sfcprop          !
!     type(GFS_coupling_type),        intent(inout) :: Coupling         !
!     type(GFS_grid_type),            intent(in)    :: Grid             !
!     type(GFS_tbd_type),             intent(inout  :: Tbd              !
!     type(GFS_cldprop_type),         intent(inout) :: Cldprop          !
!     type(GFS_radtend_type),         intent(inout) :: Radtend          !
!     type(GFS_diag_type),            intent(inout) :: Diag             !
!                                                                       !
!  subprograms called:                                                  !
!                                                                       !
!     get_prs,  dcyc2t2_pre_rad (testing),    dcyc2t3,  sfc_diff,       !
!     sfc_ocean,sfc_drv,  sfc_sice, sfc_cice, sfc_diag, moninp1,        !
!     moninp,   moninq1,  moninq,   satmedmfvdif,                       !
!     gwdps,    ozphys,   get_phi,                                      !
!     sascnv,   sascnvn,  samfdeepcnv, rascnv,   cs_convr, gwdc,        !
!     shalcvt3, shalcv,   samfshalcnv,                                  !
!     shalcnv,  cnvc90,   lrgscl,   gsmdrive, gscond,   precpd,         !
!     progt2.                                                           !
!                                                                       !
!                                                                       !
!  program history log:                                                 !
!           19xx  - ncep mrf/gfs                                        !
!           2002  - s. moorthi  modify and restructure and add Ferrier  !
!                               microphysics as an option               !
!           200x  - h. juang    modify (what?)                          !
!      nov  2004  - x. wu       modify sea-ice model                    !
!      may  2005  - s. moorthi  modify and restructure                  !
!           2005  - s. lu       modify to include noah lsm              !
!      oct  2006  - h. wei      modify lsm options to include both      !
!                               noah and osu lsms.                      !
!           2006  - s. moorthi  added a. johansson's convective gravity !
!                               wave parameterization code              !
!           2007  - s. moorthi  added j. han's modified pbl/sas options !
!      dec  2007  - xu li       modified the operational version for    !
!                               nst model                               !
!           2008  - s. moorthi  applied xu li's nst model to new gfs    !
!      mar  2008  - y.-t. hou   added sunshine duration var (suntim) as !
!                     an input/output argument.                         !
!           2008  - jun wang    added spfhmax/spfhmin as input/output.  !
!      apr  2008  - y.-t. hou   added lw sfc emissivity var (sfcemis),  !
!                     define the lw sfc dn/up fluxes in two forms: atmos!
!                     and ground. also changed sw sfc net flux direction!
!                     (positive) from ground -> atmos to the direction  !
!                     of atmos -> ground. recode the program and add    !
!                     program documentation block.
!           2008/ - s. moorthi and y.t. hou upgraded the code to more   !
!           2009    modern form and changed all the inputs to MKS units.!
!      feb  2009  - s. moorthi  upgraded to add Hochun's gocart changes !
!      jul  2009  - s. moorthi  added rqtk for sela's semi-lagrangian   !
!      aug  2009  - s. moorthi  added j. han and h. pan updated shallow !
!                               convection package                      !
!      sep  2009  - s. moorthi  updated for the mcica (rrtm3) radiation !
!      feb  2011  - sarah lu    add the option to update surface diag   !
!                               fields (t2m,q2m,u10m,v10m) at the end   !
!      Jun  2011  - s. moorthi and Xu Li - updated the nst model        !
!                               !
!      sep  2011  - sarah lu    correct dqdt_v calculations             !
!      apr  2012  - henry juang add idea                                !
!      sep  2012  - s. moorthi  merge with operational version          !
!      Mar  2013  - Jun Wang    set idea heating rate to tmp tendency   !
!      May  2013  - Jun Wang    tmp updated after idea phys             !
!      Jun  2013  - s. moorthi  corrected a bug in 3d diagnostics for T !
!      Aug  2013  - s. moorthi updating J. Whitekar's changes related   !
!                              to stochastic physics perturnbation      !
!      Oct  2013  - Xingren Wu  add dusfci/dvsfci                       !
!      Mar  2014  - Xingren Wu  add "_cpl" for coupling                 !
!      Mar  2014  - Xingren Wu  add "nir/vis beam and nir/vis diff"     !
!      Apr  2014  - Xingren Wu  add "NET LW/SW including nir/vis"       !
!      Jan  2014  - Jun Wang    merge Moorthi's gwdc change and H.Juang !
!                               and F. Yang's energy conversion from GWD!
!      jan  2014  - y-t hou     revised sw sfc spectral component fluxes!
!                     for coupled mdl, added estimation of ocean albedo !
!                     without ice contamination.                        !
!      Jun  2014  - Xingren Wu  update net SW fluxes over the ocean     !
!                               (no ice contamination)                  !
!      Jul  2014  - Xingren Wu  add Sea/Land/Ice Mask - slmsk_cpl       !
!      Jul  2014  - s. moorthi  merge with standalone gfs and cleanup   !
!      Aug  2014  - s. moorthi  add tracer fixer                        !
!      Sep  2014  - Sarah Lu    disable the option to compute tracer    !
!                               scavenging in GFS phys (set fscav=0.)   !
!      Dec  2014  - Jun Wang    add cnvqc_v for gocart                  !

!  ====================  defination of variables  ====================  !
!      ---  2014  - D. Dazlich  Added Chikira-Sugiyama (CS) convection  !
!                               as an option in opr GFS.                !
!      Apr  2015    S. Moorthi  Added CS scheme to NEMS/GSM             !
!      Jun  2015    S. Moorthi  Added SHOC  to NEMS/GSM                 !
!      Aug  2015  - Xu  Li      change nst_fcst to be nstf_name         !
!                               and introduce depth mean SST            !
!      Sep  2015  - Xingren Wu  remove oro_cpl & slmsk_cpl              !
!      Sep  2015  - Xingren Wu  add sfc_cice                            !
!      Sep  2015  - Xingren Wu  connect CICE output to sfc_cice         !
!      Jan  2016  - P. Tripp    NUOPC/GSM merge                         !
!      Mar  2016  - J. Han  -   add ncnvcld3d integer                   !
!                               for convective cloudiness enhancement   !
!      Mar  2016  - J. Han  -   change newsas & sashal to imfdeepcnv    !
!                            &  imfshalcnv, respectively                !
!      Mar  2016    F. Yang     add pgr to rayleigh damping call        !
!      Mar  2016    S. Moorthi  add ral_ts                              !
!      Mar  2016    Ann Cheng   add morrison 2m microphysics (gsfc)     !
!      May  2016    S. Moorthi  cleanup 2m microphysics implementation  !
!      Jun  2016    X. Li       change all nst_fld as inout             !
!      jul  2016    S. Moorthi  fix some bugs in shoc/2m microphysics   !
!      au-nv2016a   S. Moorthi  CS with AW and connect with shoc/2m     !
!      Dec  2016    Anning C.   Add prognostic rain and snow with 2M    !
!      Oct  2017    S. Moorthi  fix tracers to account for ice, snow etc!
!                               with this RAS and CSAW advect condensates!
!      Mar  2017    Ruiyu S.    Add Thompson's 2M aerosol MP            !
!      May  2017    Ruiyu S.    Add WSM6 MP                             !
!      Dec  2017    S. Moorthi  Merge/update Ruiyu's update on vertical !
!                               diffusion of tracers for all monins     !
!      Jan 04 2018  S. Moorthi  fix a bug in rhc for use in MG          !
!                               macrophysics and replace ntrac by nvdiff!
!                               in call to moninshoc                    !
!      Jun  2018    J. Han      Add scal-aware TKE-based moist EDMF     !
!                               vertical turbulent mixng scheme         !
!      Nov  2018    J. Han      Add canopy heat storage parameterization!
!      Feb  2019    Ruiyu S.    Add an alternate method to use          ! 
!				hydrometeors from GFDL MP in radiation  !
!      Mar  2019    Rongqian &Helin    Add Noah MP LSM                  ! 
!      Mar  2019    S. Moorthi  update slflag for MG3 and update        !
!                               rain/snow over sea-ice.  Update sfc_sice!
!                               sfc_cice calls                          !
!
!      Apr 22 2019  S. Moorthi  Porting Unified Gravitiy Wave drag      !
!                               parameterrizaion package from V. Yudin, !
!                               J. Alpert, T. Fuller-Rowll and R. Akmaev! 
!      May  2019    J. Han      Add updated scal-aware TKE-based moist  !
!                               EDMF vertical turbulent mixng scheme    !
!      july 2019    S. Moorthi  Move original GWD to inside of UGW such !
!                               that it can be called along with non-   !
!                               stationary GWD and make this part a     !
!                               function of precip or TKE.              !
!      Jul  2019    Weiguo Wang Update PBL scheme for HAFS              !
!
!  ====================    end of description    =====================
!  ====================  definition of variables  ====================  !

!> @details This subroutine is the suite driver for the GFS atmospheric physics and surface.
!! It is responsible for calculating and applying tendencies of the atmospheric state
!! variables due to the atmospheric physics and due to the surface layer scheme. In addition,
!! this routine applies radiative heating rates that were calculated during the
!! antecedent call to the radiation scheme. Code within this subroutine is executed on the
!! physics sub-timestep. The sub-timestep loop is executed in the subroutine gloopb.
!!
!!  \section general General Algorithm
!!  -# Prepare input variables for calling individual parameterizations.
!!  -# Using a two-iteration loop, calculate the state variable tendencies for the surface layer.
!!  -# Calculate the state variable tendencies due to the PBL (vertical diffusion) scheme.
!!  -# Calculate the state variable tendencies due to orographic gravity wave drag and Rayleigh damping.
!!  -# Apply tendencies to the state variables calculated so far:
!!   - for temperature: radiation, surface, PBL, oro. GWD, Rayleigh damping
!!   - for momentum: surface, PBL, oro. GWD, Rayleigh damping
!!   - for water vapor: surface, PBL
!!  -# Calculate and apply the tendency of ozone.
!!  -# Prepare input variables for physics routines that update the state variables within their subroutines.
!!  -# If SHOC is active and is supposed to be called before convection, call it and update the state variables within.
!!  -# Calculate and apply the state variable tendencies (within the subroutine) due to deep convection.
!!  -# Calculate the state variable tendencies due to convective gravity wave drag and apply them afterwards.
!!  -# Calculate and apply the state variable tendencies (within the subroutine) due to shallow convection.
!!  -# If SHOC is active and is supposed to be called after convection, call it and update the state variables within.
!!  -# Prepare for microphysics call by calculating preliminary variables.
!!  -# If necessary, call the moist convective adjustment subroutine and update the state temperature and moisture variable within.
!!  -# Calculate and apply the state variable tendencies (within the subroutine) due to microphysics.
!!  -# Determine the precipitation type and update land surface properties if necessary.
!!  -# Fill the output variables from the local variables as necessary and deallocate allocatable arrays.
!!  \section detailed Detailed Algorithm
!!  ## Prepare input variables for calling individual parameterizations.
!!  Before calling any parameterizations, there is a section at the beginning of the subroutine for
!!  preparing input arguments to the various schemes based on general input to the driver and initializing
!!  variables used throughout the driver.
!!  - General initialization:
!!   - set a flag for running in debug mode and the horizontal index of the column to print
!!   - calculate the pressure at layer centers, the exner function at layer centers and interfaces,
!!     geopotential at layer centers and interfaces, and the layer-centered pressure difference
!!   - calculate the ratio of dynamics time step to physics time step for applying tendencies
!!   - initialize local tendency arrays to zero
!!  - Radiation:
!!   - adjust radiative fluxes and heating rates to the shorter physics time step (from the longer radiation time step),
!!    unless idealized physics is true (lsidea) where radiative heating rates are set to 0
!!   - compute diagnostics from the radiation scheme needed for other schemes (e.g., downward longwave flux absorbed by the surface)
!!   - accumulate the upward and downward longwave fluxes at the surface
!!  - Surface:
!!   - set NOAH and OSU scheme variables from gbphys input arguments and initialize local soil moisture variables
!!   - set local sea ice variables from gbphys arguments
!!   - set up A/O/I coupling variables from gbphys arguments
!!  - PBL:
!!   - set the number of tracers that are diffused vertically
!!  - SHOC:
!!   - determine the index of TKE (ntk) in the convectively transported tracer array (clw)
!!   - allocate precipitation mixing ratio cloud droplet number concentration arrays
!!  - Deep Convection:
!!   - determine which tracers in the tracer input array undergo convective transport (valid for the RAS and Chikira-Sugiyama, and SAMF schemes) and allocate a local convective transported tracer array (clw)
!!   - apply an adjustment to the tracers from the dynamics
!!   - calculate horizontal grid-related parameters needed for some parameterizations
!!   - calculate the maxiumum cloud base updraft speed for the Chikira-Sugiyama scheme
!!   - allocate array for cloud water and cloud cover (for non-RAS and non-Chikira-Sugiyama deep convective schemes)
!!  - Shallow Convection:
!!   - when using the Tiedtke shallow convection scheme with the stratus modifications, find the lowest
!!     model level where a temperature inversion exists in the absence of CTEI
!!  - Microphysics:
!!   - for the Morrison (MGB) scheme, calculate 'FRLAND' if the grid point is over land
!!   - allocate arrays associated with the Morrison scheme
!!   - assign the local critical relative humidity variables from the gbphys arguments
!!  - Gravity Wave Drag:
!!   - calculate the deep convective cloud fraction at cloud top for the convective GWD scheme
!!  .
!!  ## Using a two-iteration loop, calculate the state variable tendencies for the surface layer.
!!   - Each iteration of the loop calls the following:
!!    - 'sfc_diff' to calculate surface exchange coefficients and near-surface wind
!!    - surface energy balances routines are called regardless of surface type; the surface type is checked within each to determine whether the routine is "active"
!!    - for the surface energy balance over the ocean, call 'sfc_nst' if NSST is on, otherwise, call 'sfc_ocean'
!!    - for the surface energy balance over the land, call 'sfc_drv' for the NOAH model and 'sfc_land' for the OSU model
!!    - for the surface energy balance over sea ice, call sfc_sice; if A/O/I coupling, call sfc_cice
!!   - The initial iteration has flag_guess = F unless wind < 2 m/s; flag_iter = T
!!   - After the initial iteration, flag_guess = F and flag_iter = F (unless wind < 2 m/s and over a land surface or an ocean surface with NSST on)
!!   - The following actions are performed after the iteration to calculate surface energy balance:
!!    - set surface output variables from their local values
!!    - call 'sfc_diag' to calculate state variable values at 2 and 10 m as appropriate from near-surface model levels and the surface exchange coefficients
!!    - if A/O/I coupling, set coupling variables from local variables and calculate the open water albedo
!!    - finally, accumulate surface-related diagnostics and calculate the max/min values of T and q at 2 m height.
!!  .
!!  ## Calculate the state variable tendencies due to the PBL (vertical diffusion) scheme.
!!   - Call the vertical diffusion scheme (PBL) based on the following logical flags: do_shoc, hybedmf, satmedmf, old_monin, mstrat
!!    - the PBL scheme is expected to return tendencies of the state variables
!!   - If A/O/I coupling and the surface is sea ice, overwrite some surface-related variables to their states before PBL was called
!!   - For diagnostics, do the following:
!!    - accumulate surface state variable tendencies and set the instantaneous values for output
!!    - accumulate the temperature tendency due to the PBL scheme in dt3dt(:,:,3), subtracting out the radiative heating rate if necessary
!!    - accumulate the u, v tendencies due to the PBL in du3dt(:,:,1:2) and dv3dt(:,:,1:2)
!!    - accumulate the water vapor tendency due to the PBL in dq3dt(:,:,1)
!!    - accumulate the ozone tendency in dq3dt(:,:,5)
!!  .
!!  ## Calculate the state variable tendencies due to orographic gravity wave drag and Rayleigh damping.
!!   - Based on the variable nmtvr, unpack orographic gravity wave varibles from the hprime array
!!   - Call 'gwdps' to calculate tendencies of u, v, T, and surface stress
!!   - Accumulate gravity wave drag surface stresses.
!!   - Accumulate change in u, v, and T due to oro. gravity wave drag in du3dt(:,:,2), dv3dt(:,:,2), and dt3dt(:,:,2)
!!   - Call 'rayleigh_damp' to calculate tendencies to u, v, and T due to Rayleigh friction
!!  .
!!  ## Apply tendencies to the state variables calculated so far.
!!  ## Calculate and apply the tendency of ozone.
!!   - Call the convective adjustment scheme for IDEA
!!   - Call 'ozphys_2015' or 'ozphys' depending on the value of pl_coeff, updating the ozone tracer within and outputing the tendency of ozone in dq3dt(:,:,6)
!!   - Call 'h2ophys' if necessary ("adaptation of NRL H2O phys for stratosphere and mesophere")
!!  .
!!  ## Prepare input variables for physics routines that update the state variables within their subroutines.
!!  - If diagnostics is active, save the updated values of the state variables in 'dudt', 'dvdt', 'dTdt', and 'dqdt(:,:,1)'
!!  - Call 'get_phi' to calculate geopotential from p, q, T
!!  - Initialize the cloud water and ice portions of the convectively transported tracer array (clw) and (if the deep convective scheme is not RAS or Chikira-Sugiyama) the convective cloud water and cloud cover.
!!  - If the dep convective scheme is RAS or Chikira-Sugiyama, fill the 'clw' array with tracers to be transported by convection
!!  - Initialize 'ktop' and 'kbot' (to be modified by all convective schemes except Chikira-Sugiyama)
!!  - Prepare for microphysics call (if cloud condensate is in the input tracer array):
!!   - all schemes: calculate critical relative humidity
!!   - Morrison et al. scheme (occasionally denoted MGB) (when ncld==2): set clw(:,:,1) to cloud ice and clw(:,:,2) to cloud liquid water
!!   - Ferrier scheme (num_p3d==3): set the cloud water variable and separate hydrometeors into cloud ice, cloud water, and rain; set clw(:,:,1) to cloud ice and clw(:,:,2) to cloud liquid water
!!   - Zhao-Carr scheme (num_p3d==4): calculate autoconversion coefficients from input constants and grid info; set set clw(:,:,1) to cloud liquid water
!!   - otherwise: set autoconversion parameters like in Zhao-Carr and set critical relative humidity to 1
!!  .
!!  ##  If SHOC is active and is supposed to be called before convection, call it and update the state variables within.
!!   - Prior to calling SHOC, prepare some microphysics variables:
!!    - if Morrison et al. scheme: set 'skip_macro', fill clw(:,:,1,2) with cloud ice, liquid from the tracer array, and fill cloud droplet number concentration arrays from the input tracer array
!!    - if Zhao-Carr scheme: initialize precip. mixing ratios to 0, fill clw(:,:,1,2) with cloud ice, liquid from the tracer array (as a function of temperature)
!!   - Call 'shoc' (modifies state variables within the subroutine)
!!   - Afterward, set updated cloud number concentrations in the tracer array from the updated 'ncpl' and 'ncpi'
!!  .
!!  ## Calculate and apply the state variable tendencies (within the subroutine) due to deep convection.
!!   - Call deep convective scheme according to the parameter 'imfdeepcnv', 'ras', and 'cscnv'.
!!    - if imfdeepcnv == 0, 1, or 2, no special processing is needed
!!    - if the Chikira-Sugiyama scheme (cscnv), convert rain rate to accumulated rain (rain1)
!!    - if RAS, initialize 'ccwfac', 'dlqfac', and revap before the call to 'rascnv'
!!   - Zero out 'cld1d' (cloud work function calculated in non-RAS, non-Chikira-Sugiyama schemes)
!!   - Update tracers in the tracer array (gq0) due to convective transport (RAS, CS only) from the 'clw' array
!!   - Calculate accumulated surface convective precip. for this physics time step (rainc)
!!   - If necessary, accumulate cloud work function, convective precipitation, and convective mass fluxes; accumulate dt3dt(:,:,4), dq3dt(:,:,2), du3dt(:,:,3), dv3dt(:,:,3) as change in state variables due to deep convection
!!   - If PDF-based clouds are active and Zhao-Carr microphysics, save convective cloud cover and water in 'phy_f3d' array
!!    - otherwise, if non-PDF-based clouds and the "convective cloudiness enhancement" is active, save convective cloud water in 'phy_f3d' array
!!  .
!!  ## Calculate the state variable tendencies due to convective gravity wave drag and apply them afterwards.
!!   - Calculate the average deep convective heating rate in the column to pass into 'gwdc'
!!   - Call 'gwdc' to calculate tendencies of u, v due to convective GWD
!!   - For diagnostics, accumulate the vertically-integrated change in u, v due to conv. GWD; accumulate change in u, v, due to conv. GWD in du3dt(:,:,4) and dv3dt(:,:,4)
!!   - Calculate updated values of u, v, T using conv. GWD tendencies
!!  .
!!  ## Calculate and apply the state variable tendencies (within the subroutine) due to shallow convection.
!!   - If diagnostics are active, set 'dtdt' and 'dqdt' to updated values of T and q before shallow convection
!!   - If SHOC is not active, do the following:
!!    - for the mass-flux shallow convection scheme (imfshalcnv == 1), call 'shalcnv'
!!    - for the scale- and aerosol-aware scheme (imfshalcnv == 2), call 'samfshalcnv'
!!    - for either of the first two schemes, perform the following after the call:
!!     - if Zhao-Carr microphysics with PDF-based clouds, save convective cloud water an cover in 'phy_f3d'
!!     - if non-PDF-based clouds and convective cloudiness enhancement is active, save convective cloud water in 'phy_f3d'
!!     - calculate shallow convective precip. and add it to convective precip; accumulate convective precip.
!!    - for the Tiedtke scheme (imfshalcnv == 0), find the top level where shallow convection must stratosphere
!!     - if using Moorthi's approach to stratus, call 'shalcv'
!!     - otherwise, call 'shalcvt3'
!!     - save the change in T and q due to shallow convection in dt3dt(:,:,5) and dq3dt(:,:,3); reset dtdt and dqdt to the updated values of T, q after shallow Convection
!!     - if 'clw' is not partitioned into ice/water, set 'clw(ice)' to zero
!!   - If SHOC is active (and shocaftcnv)
!!    - if Morrison et al. scheme: set 'skip_macro' and fill cloud droplet number concentration arrays from the input tracer array
!!    - initialize precip. mixing ratios to 0
!!    - call 'shoc' (modifies state variables within the subroutine)
!!    - afterward, set updated cloud number concentrations in the tracer array from the updated 'ncpl' and 'ncpi'
!!  .
!!  ## Prepare for microphysics call by calculating preliminary variables.
!!   - For Morrison et al. microphysics, set cloud water and ice arrays to the convecitvely transported values
!!   - For Ferrier microphysics, combine convectively transported cloud water and ice with column rain and save in cloud water array
!!    - calculate and save ice fraction and rain fraction in phy_f3d(1),(2)
!!   - For Zhao-Carr, combine convectively transported cloud water and ice into the cloud water array
!!   - Otherwise, combine convectively transported cloud water and ice into the convectively transported cloud water
!!   - Call 'cnvc90'; a "legacy routine which determines convective clouds"; outputs 'acv','acvb','acvt','cv','cvb','cvt'
!!  .
!!  ## If necessary, call the moist convective adjustment subroutine and update the state temperature and moisture variable within.
!!   - Updates T, q, 'rain1', cloud water array
!!   - Accumulate convective precip
!!   - For diagnostics, accumulate the change in T, q due to moist convective adjustment; reset 'dtdt' and 'dqdt' to updated T, q before call to microphysics
!!  .
!!  ## Calculate and apply the state variable tendencies (within the subroutine) due to microphysics.
!!   - If no cloud microphysics (ncld == 0), call 'lrgscl' to update T, q and output large scale precipitation and cloud water
!!   - Otherwise, a more advanced microphysics scheme is called (which scheme depends on values of 'num_p3d','do_shoc',and 'ncld')
!!   - Ferrier scheme (num_p3d == 3):
!!    - calculate droplet number concentration and minimum large ice fraction
!!    - call 'gsmdrive' (modifies T, q, cloud water, 'f_ice', 'f_rain', 'f_rimef', 'rain1')
!!   - Zhao-Carr-Sundqvist scheme (num_p3d == 4):
!!    - if non-PDF-based clouds:
!!     - if 'do_shoc', call 'precpd_shoc' (precpd modified for SHOC)
!!     - else, call 'gscond' (grid-scale condensation/evaporation); updates water vapor, cloud water, temperature
!!      - call 'precpd'; updates water vapor, cloud water, temperature and outputs precip., snow ratio, and rain water path
!!    - for PDF-based clouds:
!!     - call 'gscondp' followed by 'precpdp' (similar arguments to gscond, precpd above)
!!   - Morrison et al. scheme (ncld = 2):
!!    - if 'do_shoc', set clw(1),(2) from updated values; set phy_f3d(:,:,1) from phy_f3d(:,:,ntot3d-2)
!!    - else, set clw(1),(2) from updated values; set phy_f3d(:,:,1) to cloud cover from previous time step + convective cloud water from convective scheme
!!    - call 'm_micro_driver'; updates water vapor, temperature, droplet number concentrations, cloud cover
!!   - Combine large scale and convective precip.
!!   - For diagnostics, accumulate total surface precipitation and accumulate change in T and q due to microphysics in dt3dt(:,:,6) and dq3dt(:,:,4)
!!  .
!!  ## Determine the precipitation type and update land surface properties if necessary.
!!   - If 'cal_pre', diagnose the surface precipitation type
!!    - call 'calpreciptype'; set snow flag to 1 if snow or sleet, 0 otherwise
!!   - For rain/snow decision, calculate temperature at 850 hPa (\f$T_{850}\f$)
!!    - If not 'cal_pre', set snow flag to 1 if \f$T_{850}\f$ is below freezing
!!   - For coupling, accumulate rain if \f$T_{850}\f$ is above freezing, otherwise accumulate snow
!!   - If using the OSU land model, accumulate surface snow depth if \f$T_{850}\f$ is below freezing and not over an ocean surface
!!    - call 'progt2' (canopy and soil moisture?) and set the soil liquid water equal to soil total water
!!   - If necessary (lssav), update the 2m max/min values of T and q
!!   - If necessary (lssav), accumulate total runoff and surface runoff.
!!  .
!!  ## Fill the output variables from the local variables as necessary and deallocate allocatable arrays.
!!   - Set global sea ice thickness and concentration as well as the temperature of the sea ice
!!   - Set global soil moisture variables
!!   - Calculate precipitable water and water vapor mass change due to all physics for the column
!!   - Deallocate arrays for SHOC scheme, deep convective scheme, and Morrison et al. microphysics


  public GFS_physics_driver

  CONTAINS
!*******************************************************************************************

    subroutine GFS_physics_driver                                 &
                   (Model, Statein, Stateout, Sfcprop, Coupling,  &
                   Grid, Tbd, Cldprop, Radtend, Diag)

      implicit none
!
!  ---  interface variables
      type(GFS_control_type),         intent(inout) :: Model
      type(GFS_statein_type),         intent(inout) :: Statein
      type(GFS_stateout_type),        intent(inout) :: Stateout
      type(GFS_sfcprop_type),         intent(inout) :: Sfcprop
      type(GFS_coupling_type),        intent(inout) :: Coupling
      type(GFS_grid_type),            intent(inout) :: Grid
      type(GFS_tbd_type),             intent(inout) :: Tbd
      type(GFS_cldprop_type),         intent(inout) :: Cldprop
      type(GFS_radtend_type),         intent(inout) :: Radtend
      type(GFS_diag_type),            intent(inout) :: Diag
!
!## CCPP ## Note: Variables defined locally in this file for temporary calculations
! or transfer of data between schemes are defined in gfsphysics/GFS_layer/GFS_typedefs.F90
! in the GFS_interstitial_type datatype. Type-bound procedures create, rad_reset, 
! phys_reset, and mprint exist to allocate memory, to reset variables used in GFS_radiation_driver.F90,
! to reset variables used in GFS_physics_driver.F90, and to print the contents of the
! data type to the console

!  ---  local variables

!--- INTEGER VARIABLES
      integer :: me,   ipr,  ix,   im, levs, ntrac, nvdiff, kdt,        &
                 ntoz, ntcw, ntiw, ncld,ntke,ntkev, ntlnc, ntinc, lsoil,&
                 ntrw, ntsw, ntrnc, ntsnc, ntot3d, ntgl, ntgnc, ntclamt,&
                 ims, ime, kms, kme, its, ite, kts, kte, imp_physics,   &
                 ntwa, ntia, nmtvr

      integer :: i, kk, ic, itc, k, n, k1, iter, levshcm, tracers,      &
                 tottracer, nsamftrac, num2, num3, nshocm, nshoc, ntk,  &
                 nn, nncl, ntiwx, seconds

      integer, dimension(size(Grid%xlon,1)) ::                          &
           kbot, ktop, kcnv, soiltyp, vegtype, kpbl, slopetyp, kinver,  &
           levshc, islmsk,                                              &
!--- coupling inputs for physics
           islmsk_cice

!--- LOGICAL VARIABLES
      logical :: lprnt, revap, mg3_as_mg2, skip_macro, trans_aero

      logical, dimension(size(Grid%xlon,1)) ::                          &
           flag_iter, flag_guess, invrsn,                               &
!--- coupling inputs for physics
           flag_cice

      logical, dimension(Model%ntrac+1,2) :: otspt

      real(kind=kind_phys), dimension(Model%ntrac+2) :: trcmin

!--- REAL VARIABLES
      real(kind=kind_phys) ::                                           &
           dtf, dtp,  frain, tem,   tem1, tem2,                         &
           xcosz_loc, zsea1, zsea2, eng0, eng1, dpshc,                  &
           txl, txi, txo, dt_warm,                                      &
!--- experimental for shoc sub-stepping
           dtshoc,                                                      &
!--- GFDL Cloud microphysics
           crain, csnow, total_precip

      real(kind=kind_phys) :: rho


      real(kind=kind_phys), dimension(Model%ntrac-Model%ncld+2) ::      &
           fscav, fswtr

      real(kind=kind_phys), dimension(size(Grid%xlon,1))  ::            &
           ccwfac, garea, dlength, cumabs, fice, zice, tice,       gflx,&
           rain1,         snowmt, cd, cdq, qss, dusfcg, dvsfcg, dusfc1, &
           dvsfc1,  dtsfc1, dqsfc1, rb, drain,  cld1d, evap, hflx,      &
           stress, t850, ep1d, gamt, gamq, sigmaf,                      &
           wind, work1, work2, work3, work4, runof, xmu, fm10, fh2,     &
                   tx1, tx2, tx3, tx4, ctei_r, evbs, evcw, trans, sbsno,&
           snowc, frland, adjsfcdsw, adjsfcnsw, adjsfcdlw, adjsfculw,   &
           adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu, adjnirbmd,       &
           adjnirdfd, adjvisbmd, adjvisdfd,           xcosz, tseal,     &
!          adjnirdfd, adjvisbmd, adjvisdfd, gabsbdlw, xcosz, tseal,     &
           snohf, dlqfac, ctei_rml, cldf, domr, domzr, domip,           &
           doms, psautco_l, prautco_l, ocalnirbm_cpl, ocalnirdf_cpl,    &
           ocalvisbm_cpl, ocalvisdf_cpl, dtzm, temrain1, t2mmp, q2mp,   &
           psaur_l, praur_l,                                            &
!--- for CS-convection
           wcbmax

!  1 - land, 2 - ice, 3 - ocean
      real(kind=kind_phys), dimension(size(Grid%xlon,1),3)  ::           &
             zorl3, cd3, cdq3, rb3, stress3, ffmm3, ffhh3, uustar3,      &
             fm103, fh23, qss3, cmm3, chh3, gflx3, evap3, hflx3, ep1d3,  &
             weasd3, snowd3, tprcp3, tsfc3, tsurf3, adjsfculw3, semis3,  &
             gabsbdlw3

      logical, dimension(size(Grid%xlon,1))                ::           &
           wet, dry,              icy

      real(kind=kind_phys), dimension(size(Grid%xlon,1),1) ::           &
          area, land, rain0, snow0, ice0, graupel0

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%lsoil) :: &
          smsoil, stsoil, slsoil

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs) ::  &
          del, rhc, dtdt, dudt, dvdt, dtdtc,                            &
          ud_mf, dd_mf, dt_mf, prnum, dkt
!         ud_mf, dd_mf, dt_mf, prnum, dkt, sigmatot, sigmafrac, txa
      real(kind=kind_phys), allocatable, dimension(:,:) :: sigmatot,    &
          gwdcu, gwdcv, rainp, sigmafrac, tke


!--- GFDL modification for FV3 

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs+1) ::&
           del_gz
      real(kind=kind_phys), allocatable, dimension(:,:,:) ::            &
           delp, dz, uin, vin, pt, qv1, ql1, qr1, qg1, qa1, qn1, qi1,   &
           qs1, pt_dt, qa_dt, udt, vdt, w, qv_dt, ql_dt, qr_dt, qi_dt,  &
           qs_dt, qg_dt, p123, refl
!
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs,Model%ntrac) :: &
           dqdt

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs,oz_coeff+5) ::  &
           dq3dt_loc

!  mg, sfc perts
      real (kind=kind_phys), dimension(size(Grid%xlon,1)) :: &
         z01d, zt1d, bexp1d, xlai1d, alb1d, vegf1d
      real(kind=kind_phys) :: cdfz
!--- ALLOCATABLE ELEMENTS
      !--- in clw, the first two varaibles are cloud water and ice.
      !--- from third to ntrac are convective transportable tracers,
      !--- third being the ozone, when ntrac=3 (valid with ras, csaw, or samf)
      !--- Anning Cheng 9/21/2016 leave a hook here for diagnosed snow,
      !--- rain, and their numbers
      real(kind=kind_phys), allocatable ::                              &
           clw(:,:,:), qrn(:,:),  qsnw(:,:), ncpl(:,:), ncpi(:,:),      &
           ncpr(:,:),  ncps(:,:), cnvc(:,:), cnvw(:,:),                 &
           qgl(:,:),   ncgl(:,:)
!--- for 2 M microphysics
!     real(kind=kind_phys), allocatable, dimension(:) ::                &
!            cn_prc, cn_snr
      real(kind=kind_phys), allocatable, dimension(:,:) ::              &
!            qlcn, qicn, w_upi, cf_upi, CNV_MFD, CNV_PRC3, CNV_DQLDT,   &
             qlcn, qicn, w_upi, cf_upi, CNV_MFD,           CNV_DQLDT,   &
             CLCN, CNV_FICE, CNV_NDROP, CNV_NICE
!      real(kind=kind_phys),parameter :: slope_mg = 0.02, slope_upmg = 0.02,  &
!      real(kind=kind_phys),parameter :: slope_mg = 0.02, slope_upmg = 0.04,  &
!                         turnrhcrit = 0.900, turnrhcrit_upper = 0.150
! in the following inverse of slope_mg and slope_upmg are specified
       real(kind=kind_phys),parameter :: slope_mg   = 50.0_kind_phys,   &
                                         slope_upmg = 25.0_kind_phys
!
      !--- for 2 M Thompson MP
      real(kind=kind_phys), allocatable, dimension(:,:,:) ::            &
            vdftra, dvdftra
      real(kind=kind_phys), allocatable, dimension(:,:)   ::            &
            ice00, liq0
!     real(kind=kind_phys), allocatable, dimension(:) ::  nwfa2d
      real(kind=kind_phys), parameter :: liqm = 4./3.*con_pi*1.e-12,    &
                              icem = 4./3.*con_pi*3.2768*1.e-14*890.
!===============================================================================
!
! vay ---  local variables Local PdXdt after each Physics chain
!                TdXdt total Tendency for X due to ALL GFS_physics except
!                radiance
! vay-2018 PROCESS-oriented diagnostics for 3D-fields in UGWP for COORDE
!
!          New 2D-process oriented arrays for Daily mean (6-hr aver) diagnostics
!          Diag%dXdT_pbl  Diag%dXdT_ogw  Diag%dXdT_congw Diag%dXdT_moist
!          Diag%dXdT_total
!          Additional 2D/3D diagnostic containers and arrays
!
     logical   :: ldiag_ugwp

!    real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs) ::  &
     real(kind=kind_phys)                                              &
                             Pdtdt, Pdudt, Pdvdt
!    real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs) ::  &
!                            Tdtdt, Tdudt, Tdvdt
!-----------------------------------------
! ugwp: oro-stationary + non-stationary
!-----------------------------------------
      real(kind=kind_phys), dimension(size(Grid%xlon,1))  :: hprime,   &
                                        sigma, elvmax, oc, theta, gamma
      real(kind=kind_phys), dimension(size(Grid%xlon,1),4) :: oa4, clx
      real(kind=kind_phys), dimension(size(Grid%xlon,1))   :: sgh30      !proxy for small-scale turb oro
!
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs) :: &
                               gw_dudt,  gw_dvdt, gw_dtdt, gw_kdis
!
      real(kind=kind_phys)  :: ftausec, fdaily, fwindow
      integer               :: master

! COODRE-averaged diagnostics
!
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs) :: ax_mtb, &
                            ax_ogw, ax_tms, ax_ngw
      real(kind=kind_phys), dimension(size(Grid%xlon,1))  :: &
                            tau_tms, tau_mtb, tau_ogw, tau_ngw
      real(kind=kind_phys), dimension(size(Grid%xlon,1))  :: &
                            zm_mtb, zm_ogw, zm_ngw, zm_lwb
!------------------------------------------------------
!  parameters for canopy heat storage parametrization
!------------------------------------------------------
      real(kind=kind_phys), dimension(size(Grid%xlon,1))  :: &
                                 hflxq, evapq, hffac, hefac
      real (kind=kind_phys), parameter :: z0min=0.2, z0max=1.0
      real (kind=kind_phys), parameter :: u10min=2.5, u10max=7.5
      real (kind=kind_phys), parameter :: z0ice=1.1
!
!===============================================================================

      real, allocatable, dimension(:) :: refd, REFD263K
      integer :: kdtminus1
      logical :: reset
! For computing saturation vapor pressure and rh at 2m
      real    :: pshltr,QCQ,rh02
      real(kind=kind_phys), allocatable, dimension(:,:) :: den

      !! Initialize local variables (for debugging purposes only,
      !! because the corresponding variables Interstitial(nt)%...
      !! are reset to zero every time).
      !snowmt = 0.
      !gamq   = 0.
      !gamt   = 0.
      !gflx   = 0.
      !hflx   = 0.

      !! Strictly speaking, this is not required. But when
      !! hunting for bit-for-bit differences, doing the same as
      !! in GFS_suite_stateout_reset makes life a lot easier.
      !Stateout%gt0(:,:)   = Statein%tgrs(:,:)
      !Stateout%gu0(:,:)   = Statein%ugrs(:,:)
      !Stateout%gv0(:,:)   = Statein%vgrs(:,:)
      !Stateout%gq0(:,:,:) = Statein%qgrs(:,:,:)

!## CCPP ## Note: Setting local variables from the Model DDT (without additional 
! logic attached) is not necessary with the CCPP interstitial schemes with exceptions
! noted below.

!===> ...  begin here
      ldiag_ugwp = Model%ldiag_ugwp
!
!===>
      master = Model%master

      me      = Model%me
      ix     = size(Grid%xlon,1) !## CCPP ## set in GFS_typedefs.F90/interstitial_create
      im     = size(Grid%xlon,1) !## CCPP ## set in GFS_typedefs.F90/interstitial_create
      ipr    = min(im,10)        !## CCPP ## set in GFS_typedefs.F90/interstitial_create
      levs   = Model%levs
      lsoil  = Model%lsoil
      ntrac  = Model%ntrac
      dtf    = Model%dtf
      dtp    = Model%dtp

!## CCPP ##* this block not yet in CCPP
!-------
! For COORDE-2019 averaging with fwindow, it was done before
! 3Diag fixes and averaging ingested using "fdaily"-factor
!
      ftausec = 86400.0
      fdaily  = dtp / ftausec
      if (Model%fhzero /= 0) then
        ftausec = Model%fhzero*3600
        fwindow = dtp/ftausec
        fdaily  = fwindow
      else
        print *, 'VAY Model%fhzero = 0., Bad Averaged-diagnostics '
      endif
!-------
!*## CCPP ##

      kdt     = Model%kdt
      lprnt   = Model%lprnt
!## CCPP ## see GFS_typedefs.F90/interstitial_setup_tracers for logic for setting nvdiff 
      nvdiff  = ntrac           ! vertical diffusion of all tracers!
      ntcw    = Model%ntcw
      ntoz    = Model%ntoz
      ntiw    = Model%ntiw
      ncld    = Model%ncld
      ntke    = Model%ntke
!
      ntlnc   = Model%ntlnc
      ntinc   = Model%ntinc
      ntrw    = Model%ntrw
      ntsw    = Model%ntsw
      ntrnc   = Model%ntrnc
      ntsnc   = Model%ntsnc
      ntgl    = Model%ntgl
      ntgnc   = Model%ntgnc
      ntclamt = Model%ntclamt
      ntot3d  = Model%ntot3d
      ntwa    = Model%ntwa
      ntia    = Model%ntia
      nmtvr   = Model%nmtvr

      imp_physics = Model%imp_physics

!## CCPP ##* GFS_typedefs.F90/interstitial_setup_tracers
      nncl = ncld

      ! perform aerosol convective transport and PBL diffusion
      trans_aero = Model%cplchm .and. Model%trans_trac

      if (imp_physics == Model%imp_physics_thompson) then
        if (Model%ltaerosol) then
          nvdiff = 8
        else
          nvdiff = 5
        endif
        if (Model%satmedmf) nvdiff = nvdiff + 1
        nncl = 5
      elseif (imp_physics == Model%imp_physics_wsm6) then
        nvdiff = ntrac -3
        if (Model%satmedmf) nvdiff = nvdiff + 1
        nncl = 5
      elseif (ntclamt > 0) then             ! for GFDL MP don't diffuse cloud amount
        nvdiff = ntrac - 1
      endif

      if (imp_physics == Model%imp_physics_gfdl) then
        nncl = 5
      endif

      if (imp_physics == Model%imp_physics_mg) then
        if (abs(Model%fprcp) == 1) then
          nncl = 4                          ! MG2 with rain and snow
          mg3_as_mg2 = .false.
        elseif (Model%fprcp >= 2) then
          if (ntgl > 0 .and. (Model%mg_do_graupel .or. Model%mg_do_hail)) then
            nncl = 5                        ! MG3 with rain and snow and grapuel/hail
            mg3_as_mg2 = .false.
          else                              ! MG3 code run without graupel/hail i.e. as MG2
            nncl = 4
            mg3_as_mg2 = .true.
          endif
        endif
      endif
!
      if (Model%cplchm) then
        ! Only Zhao/Carr/Sundqvist and GFDL microphysics schemes are supported
        ! when coupling with chemistry. PBL diffusion of aerosols is only  supported
        ! Adding MG microphysics - Moorthi
        if (imp_physics == Model%imp_physics_zhao_carr) then
          nvdiff = 3
        elseif (imp_physics == Model%imp_physics_mg) then
          if (ntgl > 0) then
            nvdiff = 12
          else
            nvdiff = 10
          endif
        elseif (imp_physics == Model%imp_physics_gfdl) then
          nvdiff = 7
        endif
        if (trans_aero) nvdiff = nvdiff + Model%ntchm
        if (ntke > 0)   nvdiff = nvdiff + 1    ! adding tke to the list
      endif
!*## CCPP ##
!
!## CCPP ##* GFS_typedefs.F90/interstitial_phys_reset
      kdtminus1 = kdt - 1
      reset     = mod(kdtminus1, nint(Model%avg_max_length/dtp)) == 0
!*## CCPP ##

!
!-------------------------------------------------------------------------------------------
!     lprnt   = .false.

!     do i=1,im
!       lprnt = Model%me == 23 .and. i == 25
!       lprnt = Model%me == 127 .and. i == 11
!       lprnt = kdt >=  20 .and. abs(grid%xlon(i)*rad2dg-102.65) < 0.101  &
!                          .and. abs(grid%xlat(i)*rad2dg-0.12) < 0.201
!       lprnt = kdt >=  20 .and. abs(grid%xlon(i)*rad2dg-184.00) < 0.301  &
!                          .and. abs(grid%xlat(i)*rad2dg-83.23) < 0.301
!       lprnt = kdt >=   7 .and. abs(grid%xlon(i)*rad2dg-216.20) < 0.101  &
!                          .and. abs(grid%xlat(i)*rad2dg-81.23) < 0.101
!       lprnt = kdt >=   7 .and. abs(grid%xlon(i)*rad2dg-28.800) < 0.101  &
!                          .and. abs(grid%xlat(i)*rad2dg+2.45) < 0.101
!       lprnt = kdt >=   1 .and. abs(grid%xlon(i)*rad2dg-293.91) < 0.101  &
!                          .and. abs(grid%xlat(i)*rad2dg+72.02) < 0.101
!       lprnt = kdt >=   1 .and. abs(grid%xlon(i)*rad2dg-113.48) < 0.101  &
!                          .and. abs(grid%xlat(i)*rad2dg-21.07) < 0.101
!       lprnt = kdt >=   1 .and. abs(grid%xlon(i)*rad2dg-169.453) < 0.501  &
!                          .and. abs(grid%xlat(i)*rad2dg-72.96) < 0.501
!       if (kdt == 1) &
!         write(2000+me,*)' i=',i,' xlon=',grid%xlon(i)*rad2dg,          &
!                       ' xlat=',grid%xlat(i)*rad2dg,' me=',me
!       if (lprnt) then
!         ipr = i
!         write(0,*)' ipr=',ipr,'xlon=',grid%xlon(i)*rad2dg,' xlat=',grid%xlat(i)*rad2dg,' me=',me
!         exit
!       endif
!     enddo
!     if (lprnt) then
!       if (Model%cplflx) then
!         write(0,*)' sfcprop%tisfc=',Sfcprop%tisfc(ipr),' kdt=',kdt,   &
!    ' fice=',Sfcprop%fice(ipr),' ulw=',Coupling%ulwsfcin_cpl(ipr),     &
!    ' tsfc=',Sfcprop%tsfc(ipr)
!        else
!          write(0,*)' sfcprop%tisfc=',Sfcprop%tisfc(ipr),' kdt=',kdt,  &
!    ' fice=',Sfcprop%fice(ipr), ' tsfc=',Sfcprop%tsfc(ipr),            &
!    'tsfcl=',Sfcprop%tsfcl(ipr),' tsfco=',Sfcprop%tsfco(ipr)
!        endif
!       if (Model%nstf_name(1) > 0) then
!          write(0,*)' begin sfcprop%tref=',Sfcprop%tref(ipr),' kdt=',kdt, &
!                    ' landfrac=',Sfcprop%landfrac(ipr)
!       endif
!      endif
!-------------------------------------------------------------------------------------------
!
!     if (lprnt) then
!       write(0,*)' in phydrv tgrs=',Statein%tgrs(ipr,:)
!       write(0,*)' in phydrv ugrs=',Statein%ugrs(ipr,:)
!       write(0,*)' in phydrv vgrs=',Statein%vgrs(ipr,:)
!       write(0,*)' in phydrv qgrs=',Statein%qgrs(ipr,:,1)*1000.0
!       write(0,*)' in phydrv tke=',Statein%qgrs(ipr,:,ntke)
!       write(0,*)' in phydrv phii=',Statein%phii(ipr,:)
!     endif
!
!  --- ...  frain=factor for centered difference scheme correction of rain amount.

      frain = dtf / dtp

!## CCPP ##* GFS_typedefs.F90/interstitial_create
      skip_macro = .false.
!*## CCPP ##
!## CCPP ##* GFS_typedefs.F90/interstitial_setup_tracers
      if (ntiw > 0) then
        if (ntclamt > 0) then
          nn = ntrac - 2
        else
          nn = ntrac - 1
        endif
      elseif (ntcw > 0) then
        nn = ntrac
      else
        nn = ntrac + 1
      endif
!*## CCPP ##
!## CCPP ##* GFS_typedefs.F90/interstitial_create
      allocate (clw(ix,levs,nn))
!*## CCPP ##
!## CCPP ##* GFS_typedefs.F90/interstitial_create Note: cnvc and cnvw are always allocated and initialized regardless of test condition
      if (Model%imfdeepcnv >= 0 .or.  Model%imfshalcnv > 0  .or. &
         (Model%npdf3d == 3     .and. Model%num_p3d   == 4) .or. &
         (Model%npdf3d == 0     .and. Model%ncnvcld3d == 1) ) then
        allocate (cnvc(ix,levs), cnvw(ix,levs))
        do k=1,levs
          do i=1,im
            cnvc(i,k) = zero
            cnvw(i,k) = zero
          enddo
        enddo
!*## CCPP ##
!## CCPP ##* GFS_typedefs.F90/control_initialize Note: these are calculated regardless of test condition
        if (Model%npdf3d == 3 .and. Model%num_p3d == 4) then
          num2 = Model%num_p3d + 2
          num3 = num2 + 1
        elseif (Model%npdf3d == 0 .and. Model%ncnvcld3d == 1) then
          num2 = Model%num_p3d + 1
        endif
        !CCPP: num2 = Model%ncnvw
        !CCPP: num3 = Model%ncnvc
      endif
!*## CCPP ##

!## CCPP ##* GFS_surface_generic.F90/GFS_surface_generic_pre_run
!  ---  set initial quantities for stochastic physics deltas
      if (Model%do_sppt .or. Model%ca_global)then
        Tbd%dtdtr = zero
      endif

! mg, sfc-perts
!  ---  scale random patterns for surface perturbations with perturbation size
!  ---  turn vegetation fraction pattern into percentile pattern
!## CCPP ##* Note: initialzations to zero are not needed in GFS_surface_generic.F90/GFS_surface_generic_pre_run
! since this function occurs in GFS_typedefs.F90/interstitial_phys_reset
      do i=1,im
         z01d(i)   = zero
         zt1d(i)   = zero
         bexp1d(i) = zero
         xlai1d(i) = zero
!        alb1d(i)  = zero
         vegf1d(i) = zero
      enddo
      if (Model%do_sfcperts) then
        if (Model%pertz0(1) > zero) then
          z01d(:) = Model%pertz0(1) * Coupling%sfc_wts(:,1)
!          if (me == 0) print*,'Coupling%sfc_wts(:,1) min and max',minval(Coupling%sfc_wts(:,1)),maxval(Coupling%sfc_wts(:,1))
!          if (me == 0) print*,'z01d min and max ',minval(z01d),maxval(z01d)
        endif
        if (Model%pertzt(1) > zero) then
          zt1d(:) = Model%pertzt(1) * Coupling%sfc_wts(:,2)
        endif
        if (Model%pertshc(1) > zero) then
          bexp1d(:) = Model%pertshc(1) * Coupling%sfc_wts(:,3)
        endif
        if (Model%pertlai(1) > zero) then
          xlai1d(:) = Model%pertlai(1) * Coupling%sfc_wts(:,4)
        endif
! --- do the albedo percentile calculation in GFS_radiation_driver instead --- !
!        if (Model%pertalb(1) > zero) then
!          do i=1,im
!            call cdfnor(Coupling%sfc_wts(i,5),cdfz)
!            alb1d(i) = cdfz
!          enddo
!        endif
        if (Model%pertvegf(1) > zero) then
          do i=1,im
            call cdfnor(Coupling%sfc_wts(i,6),cdfz)
            vegf1d(i) = cdfz
          enddo
        endif
      endif
!*## CCPP ##
!
!## CCPP ##* GFS_typedefs.F90/interstitial_create
      if (Model%do_shoc) then
        allocate (qrn(im,levs),   qsnw(im,levs), &
                  ncpl(im,levs),  ncpi(im,levs))
        do k=1,levs
          do i=1,im
            ncpl(i,k) = zero
            ncpi(i,k) = zero
            qrn(i,k)  = zero
            qsnw(i,k) = zero
          enddo
        enddo
      endif
!## CCPP ##* GFS_typedefs.F90/coupling_create ##
      if (imp_physics == Model%imp_physics_thompson) then
        if(Model%ltaerosol) then
          allocate(ice00(im,levs))
          allocate(liq0(im,levs))
!         allocate(nwfa2d(im))
        else
          allocate(ice00(im,levs))
        endif
      endif
!*## CCPP ##
!## CCPP ##* allocated in GFS_typedefs.F90/interstitial_create; initialized in GFS_typedefs.F90/interstitial_phys_reset
      if (imp_physics == Model%imp_physics_mg) then         ! For MGB double moment microphysics
        allocate (qlcn(im,levs),      qicn(im,levs),    w_upi(im,levs),     &
                  cf_upi(im,levs),    CNV_MFD(im,levs),                     &
!                 cf_upi(im,levs),    CNV_MFD(im,levs), CNV_PRC3(im,levs),  &
                  CNV_DQLDT(im,levs), clcn(im,levs),    cnv_fice(im,levs),  &
                  cnv_ndrop(im,levs), cnv_nice(im,levs))
!       allocate (cn_prc(im),    cn_snr(im))
        allocate (ncpr(im,levs), ncps(im,levs), ncgl(im,levs))
        if (.not. allocated(qrn))  allocate (qrn(im,levs))
        if (.not. allocated(qsnw)) allocate (qsnw(im,levs))
        if (.not. allocated(qgl))  allocate (qgl(im,levs))
        do k=1,levs
          do i=1,im
            qrn(i,k)  = zero
            qsnw(i,k) = zero
            qgl(i,k)  = zero
            ncpr(i,k) = zero
            ncps(i,k) = zero
            ncgl(i,k) = zero
          enddo
        enddo
!*## CCPP ##
!## CCPP ##* These variables are currently being allocated fully (im,levs) in GFS_typedefs.F90/interstitial_create
      else
        allocate (qlcn(1,1),    qicn(1,1),     w_upi(1,1),    cf_upi(1,1),  &
                  CNV_MFD(1,1),                CNV_DQLDT(1,1),              &
!                 CNV_MFD(1,1), CNV_PRC3(1,1), CNV_DQLDT(1,1),              &
                  clcn(1,1),    cnv_fice(1,1), cnv_ndrop(1,1), cnv_nice(1,1))
!## CCPP ##* The following variables are local to gfdl_cloud_microphys.F90/gfdl_cloud_microphys_run
        if (imp_physics == Model%imp_physics_gfdl) then       ! GFDL MP
          allocate (delp(im,1,levs),  dz(im,1,levs),    uin(im,1,levs),                    &
                    vin(im,1,levs),   pt(im,1,levs),    qv1(im,1,levs),   ql1(im,1,levs),  &
                    qr1(im,1,levs),   qg1(im,1,levs),   qa1(im,1,levs),   qn1(im,1,levs),  &
                    qi1(im,1,levs),   qs1(im,1,levs),   pt_dt(im,1,levs), qa_dt(im,1,levs),&
                    udt(im,1,levs),   vdt(im,1,levs),   w(im,1,levs),     qv_dt(im,1,levs),&
                    ql_dt(im,1,levs), qr_dt(im,1,levs), qi_dt(im,1,levs), qs_dt(im,1,levs),&
                    qg_dt(im,1,levs), p123(im,1,levs),  refl(im,1,levs),  den(im,levs))
        endif
      endif
!*## CCPP ##

!## CCPP ## Only get_prs_fv3.F90/get_prs_fv3_run is a scheme (GFS_HYDRO is assumed to be undefined)
#ifdef GFS_HYDRO
      call get_prs(im, ix, levs, ntrac, Statein%tgrs, Statein%qgrs,     &
                   Model%thermodyn_id, Model%sfcpress_id,               &
                   Model%gen_coord_hybrid, Statein%prsi, Statein%prsik, &
                   Statein%prsl, Statein%prslk, Statein%phii, Statein%phil, del)
#else
!GFDL   Adjust the geopotential height hydrostatically in a way consistent with FV3 discretization
!  if (lprnt) write(0,*)'bef get_prs_fv3 phii=',Statein%phii(ipr,:),' kdt=',kdt

      call get_prs_fv3 (ix, levs, ntrac, Statein%phii, Statein%prsi,    &
                        Statein%tgrs, Statein%qgrs, del, del_gz)
#endif
!  if (lprnt) write(0,*)'aft get_prs_fv3 phii=',Statein%phii(ipr,:)
!  if (lprnt) write(0,*)'aft get_prs_fv3 del_gz=',del_gz(ipr,:)
!*## CCPP ##

!## CCPP ##* GFS_surface_generic.F90/GFS_surface_generic_pre_run
      do i = 1, IM
        sigmaf(i) = max( Sfcprop%vfrac(i),0.01_kind_phys )
        islmsk(i) = nint(Sfcprop%slmsk(i))
        islmsk_cice(i) = islmsk(i)

        if (islmsk(i) == 2) then
          if (Model%isot == 1) then
            soiltyp(i)  = 16
          else
            soiltyp(i)  = 9
          endif
          if (Model%ivegsrc == 1) then
            vegtype(i)  = 15
          elseif(Model%ivegsrc == 2) then
            vegtype(i)  = 13
          endif
          slopetyp(i) = 9
        else
          soiltyp(i)  = int( Sfcprop%stype(i)+half )
          vegtype(i)  = int( Sfcprop%vtype(i)+half )
          slopetyp(i) = int( Sfcprop%slope(i)+half )    !! clu: slope -> slopetyp
          if (soiltyp(i)  < 1) soiltyp(i)  = 14
          if (vegtype(i)  < 1) vegtype(i)  = 17
          if (slopetyp(i) < 1) slopetyp(i) = 1
        endif
!*## CCPP ##
!  --- ...  xw: transfer ice thickness & concentration from global to local variables
!## CCPP ## global to local variable transfer not necessary for these two
        zice(i)    = Sfcprop%hice(i)
        fice(i)    = Sfcprop%fice(i)
!*## CCPP ##* 
!## CCPP ##* GFS_surface_composites.F90/GFS_surface_composites_pre_run
        tice(i)    = Sfcprop%tisfc(i)
!
!GFDL   work1(i)   = (log(coslat(i) / (nlons(i)*latr)) - dxmin) * dxinv
!GFS         Moorthi thinks this should be area and not dx
!       work1(i)   = (log(Grid%dx(i)) - dxmin) * dxinv
        work1(i)   = (log(Grid%area(i)) - dxmin) * dxinv
        work1(i)   = max(zero, min(one, work1(i)))
        work2(i)   = one - work1(i)
        Diag%psurf(i) = Statein%pgr(i)
!*## CCPP ##
!## CCPP ##* GFS_surface_generic.F90/GFS_surface_generic_pre_run
        work3(i)   = Statein%prsik(i,1) / Statein%prslk(i,1)
!*## CCPP ##
!GFDL   tem1       = con_rerth * (con_pi+con_pi)*coslat(i)/nlons(i)
!GFDL   tem2       = con_rerth * con_pi / latr
!GFDL   garea(i)   = tem1 * tem2
!## CCPP ## global to local variable transfer not necessary for these variables
        tem1       = Grid%dx(i)
        tem2       = Grid%dx(i)
        garea(i)   = Grid%area(i)
!*## CCPP ##
!## CCPP ##* gwdc.f/gwdc_pre_run
        dlength(i) = sqrt( tem1*tem1+tem2*tem2 )
        cldf(i)    = Model%cgwf(1)    * work1(i) + Model%cgwf(2)    * work2(i)
!*## CCPP ##
!## CCPP ##* cs_conv.F90/cs_conv_pre_run
        wcbmax(i)  = Model%cs_parm(1) * work1(i) + Model%cs_parm(2) * work2(i)
!*## CCPP ##
!## CCPP ##* GFS_typedefs.F90/interstitial_phys_reset
        dry(i)     = .false.
        icy(i)     = .false.
        wet(i)     = .false.
        flag_cice(i) = .false.
!*## CCPP ##
      enddo
!
!## CCPP ##* GFS_surface_generic.F90/GFS_surface_generic_pre_run
      if (Model%cplflx) then
        do i=1,im
          islmsk_cice(i) = nint(Coupling%slimskin_cpl(i))
          flag_cice(i)   = (islmsk_cice(i) == 4)
        enddo
      endif
!*## CCPP ##

!## CCPP ##* GFS_surface_composites.F90/GFS_surface_composites_pre
      if (Model%frac_grid) then
        do i = 1, IM
          frland(i) = Sfcprop%landfrac(i)
          if (frland(i) > zero) dry(i) = .true.
          if (frland(i) < one) then
            if (flag_cice(i)) then
              if (fice(i) >= Model%min_seaice) then
                icy(i)  = .true.
                if (fice(i) < one) wet(i) = .true. ! some open ocean/lake water exists
              else
                fice(i)        = zero
                flag_cice(i)   = .false.
                islmsk_cice(i) = 0
!               islmsk(i)      = 0
                wet(i) = .true. ! some open ocean/lake water exists
              endif
            else
              if (fice(i) >= Model%min_lakeice) then
                icy(i) = .true.
                if (fice(i) < one) wet(i) = .true. ! some open ocean/lake water exists
                islmsk(i) = 2
              else
                fice(i)   = zero
!               islmsk(i) = 0
                wet(i) = .true. ! some open ocean/lake water exists
              endif
            endif
            if (wet(i) .and. .not. Model%cplflx) then
              if (Sfcprop%oceanfrac(i) > zero) then
                Sfcprop%tsfco(i) = max(Sfcprop%tsfco(i), Sfcprop%tisfc(i), tgice)
              elseif (icy(i)) then
                 Sfcprop%tsfco(i) = max(Sfcprop%tisfc(i), tgice)
              endif
            endif
          else
            fice(i) = zero
          endif
        enddo
      else
        do i = 1, IM
          if (islmsk(i) == 1) then
!           Sfcprop%tsfcl(i) = Sfcprop%tsfc(i)
            dry(i)    = .true.
            frland(i) = one
            fice(i)   = zero
          else
            frland(i) = zero
            if (flag_cice(i)) then
              if (fice(i) > Model%min_seaice) then
                icy(i) = .true.
              else
                fice(i)        = zero
                flag_cice(i)   = .false.
                islmsk_cice(i) = 0
                islmsk(i)      = 0
              endif
            else
              if (fice(i) > Model%min_lakeice) then
                icy(i) = .true.
              else
                fice(i)   = zero
                islmsk(i) = 0
              endif
            endif
            if (fice(i) < one) then
              wet(i)=.true. ! some open ocean/lake water exists
              if (.not. Model%cplflx .and. icy(i))                                   &
                 Sfcprop%tsfco(i) = max(Sfcprop%tisfc(i), tgice)
            endif
          endif
        enddo
      endif
!
      do k=1,3
        do i=1,im
             cd3(i,k) = huge
            cdq3(i,k) = huge
             rb3(i,k) = huge
         stress3(i,k) = huge
           ffmm3(i,k) = huge
           ffhh3(i,k) = huge
           fm103(i,k) = huge
            fh23(i,k) = huge
            qss3(i,k) = huge
            cmm3(i,k) = huge
            chh3(i,k) = huge
           gflx3(i,k) = zero
!          gflx3(i,k) = huge
           evap3(i,k) = huge
           hflx3(i,k) = huge
           ep1d3(i,k) = huge
         uustar3(i,k) = huge
          weasd3(i,k) = huge
          snowd3(i,k) = huge
          tprcp3(i,k) = Sfcprop%tprcp(i)
           tsfc3(i,k) = huge
          tsurf3(i,k) = huge
           zorl3(i,k) = huge
!          oro3(i,k)  = Sfcprop%oro(i)
!       oro_uf3(i,k)  = Sfcprop%oro_uf(i)
     adjsfculw3(i,k)  = zero
      gabsbdlw3(i,k)  = zero
        enddo
      enddo
      zorl3(:,2) = z0ice

!     if (.not. Model%cplflx .or. .not. Model%frac_grid) then
!       if (Model%cplwav2atm) then
!         do i=1,im
!           Sfcprop%zorll(i) = Sfcprop%zorl(i)
!         enddo
!       else
!         do i=1,im
!           Sfcprop%zorll(i) = Sfcprop%zorl(i)
!           Sfcprop%zorlo(i) = Sfcprop%zorl(i)
!         enddo
!       endif
!     endif
!     if (lprnt) write(0,*)' dry=',dry(ipr),' wet=',wet(ipr),' icy=',icy(ipr) ,&
!       ' tsfco=',Sfcprop%tsfco(ipr)
      do i=1,im
        if(wet(i)) then                    ! Water
            zorl3(i,3) = Sfcprop%zorlo(i)
            tsfc3(i,3) = Sfcprop%tsfco(i)
           tsurf3(i,3) = Sfcprop%tsfco(i)
!          weasd3(i,3) = Sfcprop%weasd(i)
!          snowd3(i,3) = Sfcprop%snowd(i)
           snowd3(i,3) = zero
           weasd3(i,3) = zero
           semis3(i,3) = 0.984_kind_phys
        endif
!
        if (dry(i)) then                   ! Land
          uustar3(i,1) = Sfcprop%uustar(i)
           weasd3(i,1) = Sfcprop%weasd(i)
            zorl3(i,1) = Sfcprop%zorll(i)
            tsfc3(i,1) = Sfcprop%tsfcl(i)
           tsurf3(i,1) = Sfcprop%tsfcl(i)
           snowd3(i,1) = Sfcprop%snowd(i)
           semis3(i,1) = Radtend%semis(i)
        endif
!
        if (icy(i)) then                   ! Ice
          uustar3(i,2) = Sfcprop%uustar(i)
           weasd3(i,2) = Sfcprop%weasd(i)
            zorl3(i,2) = Sfcprop%zorli(i)
            tsfc3(i,2) = Sfcprop%tisfc(i)
           tsurf3(i,2) = Sfcprop%tisfc(i)
           snowd3(i,2) = Sfcprop%snowd(i)
            ep1d3(i,2) = zero
            gflx3(i,2) = zero
           semis3(i,2) = 0.95_kind_phys
        endif
      enddo
!*## CCPP ##

!## CCPP ## global to local variable transfer not necessary for these variables
!  --- ...  transfer soil moisture and temperature from global to local variables
      do k=1,lsoil
        do i=1,im
          smsoil(i,k) = Sfcprop%smc(i,k)
          stsoil(i,k) = Sfcprop%stc(i,k)
          slsoil(i,k) = Sfcprop%slc(i,k)          !! clu: slc -> slsoil
        enddo
      enddo
!*## CCPP ##

      do k=1,levs
        do i=1,im
          dudt(i,k)  = zero
          dvdt(i,k)  = zero
          dtdt(i,k)  = zero
          dtdtc(i,k) = zero

!## CCPP ##* GFS_typedefs.F90/interstitial_phys_reset
!vay-2018
! Pure tendency arrays w/o accumulation of Phys-tendencies from each
!      chain of GFS-physics (later add container for species)
!
!         Pdudt(i,k)  = zero
!         Pdvdt(i,k)  = zero
!         Pdtdt(i,k)  = zero

!
!ugwp-marked can be later accumulated as Pdudt Pdvdt Pdtdt
!
          gw_dudt(i,k) = zero
          gw_dvdt(i,k) = zero
          gw_dtdt(i,k) = zero
          gw_kdis(i,k) = zero
!*## CCPP ##
        enddo
      enddo
!## CCPP ##* GFS_suite_interstitial.F90/GFS_suite_interstitial_1_run
      do n=1,ntrac
        do k=1,levs
          do i=1,im
            dqdt(i,k,n) = zero
          enddo
        enddo
      enddo
!*## CCPP ##

!## CCPP ##* This block is not yet in CCPP.
!-----------------------------------------------
!vay-2018-19 ORO/UGWP process-oriented diagnostics
!
      if (ldiag_ugwp) then
        do i=1,im
          tau_tms(i) = zero  ; tau_mtb(i) = zero
          tau_ogw(i) = zero  ; tau_ngw(i) = zero
          zm_mtb(i)  = zero  ;  zm_lwb(i) = zero
          zm_ogw(i)  = zero  ;  zm_ngw(i) = zero
        enddo
        do k=1,levs
          do i=1,im
            ax_mtb(i,k)  = zero ; ax_ogw(i,k)  = zero
            ax_tms(i,k)  = zero ; ax_ngw(i,k)  = zero
          enddo
        enddo
      endif

      if (mod((kdt-1)*dtp, ftausec) == zero) then
        do i=1,im
          Diag%tau_tofd(i) = zero
          Diag%tau_mtb(i)  = zero
          Diag%tau_ogw(i)  = zero
          Diag%tau_ngw(i)  = zero
          Diag%zmtb(i)     = zero
          Diag%zlwb(i)     = zero
          Diag%zogw(i)     = zero
!         Diag%dugwd(i)    = zero
!         Diag%dvgwd(i)    = zero
        enddo
      endif
!===========================
! can be taken out by "call Diag%zero"  => call Diag(nb)%phys_zero (Model)
!       in GFS_driver.F90
!  It can be also done by hands w/o
!  relying on  FV3GFS_io_mod
!=================================
      if (ldiag_ugwp) then
!       do k=1,levs
!         do i=1,im
!           Diag%du3dt_pbl(i,k)   = zero
!           Diag%dv3dt_pbl(i,k)   = zero
!           Diag%dt3dt_pbl(i,k)   = zero
!
!           Diag%du3dt_ogw(i,k)   = zero
!           Diag%dv3dt_ogw(i,k)   = zero
!           Diag%dt3dt_ogw(i,k)   = zero

!           Diag%du3dt_mtb(i,k)   = zero
!           Diag%dv3dt_mtb(i,k)   = zero
!           Diag%dt3dt_mtb(i,k)   = zero

!           Diag%du3dt_tms(i,k)   = zero
!           Diag%dv3dt_tms(i,k)   = zero
!           Diag%dt3dt_tms(i,k)   = zero

!           Diag%du3dt_ngw(i,k)   = zero
!           Diag%dv3dt_ngw(i,k)   = zero
!           Diag%dt3dt_ngw(i,k)   = zero
!
! employed for "storage" of State%out to compute DyCore_Tendencies
!!          Diag%du3dt_cgw(i,k)   = zero
!!          Diag%dv3dt_cgw(i,k)   = zero
!!          Diag%dt3dt_cgw(i,k)   = zero

!           Diag%du3dt_moist(i,k) = zero
!           Diag%dv3dt_moist(i,k) = zero
!           Diag%dt3dt_moist(i,k) = zero

!           Diag%dudt_tot(i,k)    = zero
!           Diag%dvdt_tot(i,k)    = zero
!           Diag%dtdt_tot(i,k)    = zero

!           Diag%uav_ugwp(i,k)    = zero
!           Diag%tav_ugwp(i,k)    = zero

!
!           Tdudt(i,k)  = zero
!           Tdvdt(i,k)  = zero
!           Tdtdt(i,k)  = zero
!         enddo
!       enddo
!
        if (kdt > 1) then
          do k=1,levs
            do i=1,im
!
!---- dycore_tend =  Statein - Stateout , assuming that Statein-after Dycore and out-after Physics
!     Statein%ugrs-- "Stateout%gu0 = Diag%du3dt_cgw"
!
              Diag%dudt_tot(i,k) = (Statein%ugrs(i,k) - Diag%du3dt_cgw(i,k))*fdaily &
                                 + Diag%dudt_tot(i,k) !
              Diag%dtdt_tot(i,k) = (Statein%tgrs(i,k) - Diag%dt3dt_cgw(i,k))*fdaily &
                                 + Diag%dtdt_tot(i,k)
            enddo
          enddo
          if (kdt == -2) then
             print *, maxval(Statein%ugrs), maxval(Diag%du3dt_cgw), ' max Uin-out'
             print *, minval(Statein%ugrs), minval(Diag%du3dt_cgw), ' min Uin-out'
             print *, maxval(Statein%tgrs), maxval(Diag%dt3dt_cgw), ' max Tin-out'
             print *, minval(Statein%tgrs), minval(Diag%dt3dt_cgw), ' min Tin-out'
          endif
        endif
      endif
!===========================Above Phys-tend Diag for COORDE ======================
!*## CCPP ##

!  --- ...  initialize dtdt with heating rate from dcyc2

!  --- ...  adjust mean radiation fluxes and heating rates to fit for
!           faster model time steps.
!      sw:  using cos of zenith angle as scaling factor
!      lw:  using surface air skin temperature as scaling factor
!## CCPP ##* This is not in the CCPP yet.
      if (Model%pre_rad) then
        call dcyc2t3_pre_rad                                                &
!  ---  inputs:
           ( Model%solhr, Model%slag, Model%sdec, Model%cdec, Grid%sinlat,  &
             Grid%coslat, Grid%xlon, Radtend%coszen, Sfcprop%tsfc,          &
             Statein%tgrs(1,1), Statein%tgrs(1,1), Coupling%sfcdsw,         &
             Coupling%sfcnsw, Coupling%sfcdlw, Radtend%htrsw, Radtend%htrlw,&
             Coupling%nirbmui, Coupling%nirdfui, Coupling%visbmui,          &
             Coupling%visdfui, Coupling%nirbmdi, Coupling%nirdfdi,          &
             Coupling%visbmdi, Coupling%visdfdi, ix, im, levs,              &
!  ---  input/output:
             dtdt,                                                          &
!  ---  outputs:
             adjsfcdsw, adjsfcnsw, adjsfcdlw, adjsfculw, xmu, xcosz,        &
             adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,                    &
             adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd                     &
           )

      else
!*## CCPP ##
!** CCPP ## dcyc2.f/dcyc2t3_run Note: Check for Model%pre_rad was omitted, so this option is broken in CCPP
        call dcyc2t3                                                        &
!  ---  inputs:
           ( Model%solhr, Model%slag, Model%sdec, Model%cdec, Grid%sinlat,  &
             Grid%coslat, Grid%xlon, Radtend%coszen, tsfc3,                 &
!            Statein%tgrs(1,1), Radtend%tsflw,  Radtend%semis,              &
             Statein%tgrs(1,1), Radtend%tsflw,  semis3,                     &
             Coupling%sfcdsw,  Coupling%sfcnsw, Coupling%sfcdlw,            &
             Radtend%htrsw,    Radtend%swhc,    Radtend%htrlw, Radtend%lwhc,&
             Coupling%nirbmui, Coupling%nirdfui, Coupling%visbmui,          &
             Coupling%visdfui, Coupling%nirbmdi, Coupling%nirdfdi,          &
             Coupling%visbmdi, Coupling%visdfdi, ix, im, levs, dtf,         &
             Model%fhswr, dry, icy, wet,                                    &
!            lprnt, ipr,                                                    &
!  ---  input/output:
             dtdt, dtdtc,                                                   &
!  ---  outputs:
             adjsfcdsw, adjsfcnsw, adjsfcdlw, adjsfculw3, xmu, xcosz,       &
             adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,                    &
             adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd                     &
           )
!*## CCPP ##
!
! save temp change due to radiation - need for sttp stochastic physics
!---------------------------------------------------------------------
      endif
!
!## CCPP ##* This is not in the CCPP yet.
      if (Model%lsidea) then                       !idea jw
        dtdt(:,:) = zero
      endif
!*## CCPP ##

!  ---  convert lw fluxes for land/ocean/sea-ice models
!  note: for sw: adjsfcdsw and adjsfcnsw are zenith angle adjusted downward/net fluxes.
!        for lw: adjsfcdlw is (sfc temp adjusted) downward fluxe with no emiss effect.
!                adjsfculw is (sfc temp adjusted) upward fluxe including emiss effect.
!        one needs to be aware that that the absorbed downward lw flux (used by land/ocean
!        models as downward flux) is not the same as adjsfcdlw but a value reduced by
!        the factor of emissivity.  however, the net effects are the same when seeing
!        it either above the surface interface or below.
!
!   - flux above the interface used by atmosphere model:
!        down: adjsfcdlw;    up: adjsfculw = sfcemis*sigma*T**4 + (1-sfcemis)*adjsfcdlw
!        net = up - down = sfcemis * (sigma*T**4 - adjsfcdlw)
!   - flux below the interface used by lnd/oc/ice models:
!        down: sfcemis*adjsfcdlw;  up: sfcemis*sigma*T**4
!        net = up - down = sfcemis * (sigma*T**4 - adjsfcdlw)

!  --- ...  define the downward lw flux absorbed by ground

!## CCPP ##* GFS_surface_composites.F90/GFS_surface_composites_pre_run
      do i=1,im
        if (dry(i)) gabsbdlw3(i,1) = semis3(i,1) * adjsfcdlw(i)
        if (icy(i)) gabsbdlw3(i,2) = semis3(i,2) * adjsfcdlw(i)
        if (wet(i)) gabsbdlw3(i,3) = semis3(i,3) * adjsfcdlw(i)
      enddo
!*## CCPP ##

!## CCPP ##* GFS_suite_interstitial.F90/GFS_suite_interstitial_2_run
      if (Model%lssav) then      !  --- ...  accumulate/save output variables

!  --- ...  sunshine duration time is defined as the length of time (in mdl output
!           interval) that solar radiation falling on a plane perpendicular to the
!           direction of the sun >= 120 w/m2

        do i=1,im
          if ( xcosz(i) >= czmin ) then   ! zenth angle > 89.994 deg
            tem1 = adjsfcdsw(i) / xcosz(i)
            if ( tem1 >= 120.0_kind_phys) then
              Diag%suntim(i) = Diag%suntim(i) + dtf
            endif
          endif
        enddo

!  --- ...  sfc lw fluxes used by atmospheric model are saved for output

        if (Model%frac_grid) then
          do i=1,im
            tem = (one - frland(i)) * fice(i) ! tem = ice fraction wrt whole cell
            if (flag_cice(i)) then
              adjsfculw(i) = adjsfculw3(i,1) * frland(i)                 &
                           + Coupling%ulwsfcin_cpl(i) * tem              &
                           + adjsfculw3(i,3) * (one - frland(i) - tem)
            else
              adjsfculw(i) = adjsfculw3(i,1) * frland(i)                 &
                           + adjsfculw3(i,2) * tem                       &
                           + adjsfculw3(i,3) * (one - frland(i) - tem)
            endif
          enddo
        else
          do i=1,im
            if (dry(i)) then                     ! all land
              adjsfculw(i) = adjsfculw3(i,1)
            elseif (icy(i)) then                 ! ice (and water)
              tem = one - fice(i)
              if (flag_cice(i)) then
                if (wet(i) .and. adjsfculw3(i,3) /= huge) then
                  adjsfculw(i) = Coupling%ulwsfcin_cpl(i)*fice(i) + adjsfculw3(i,3)*tem
                else
                  adjsfculw(i) = Coupling%ulwsfcin_cpl(i)
                endif
              else
                if (wet(i) .and. adjsfculw3(i,3) /= huge) then
                  adjsfculw(i) = adjsfculw3(i,2)*fice(i) + adjsfculw3(i,3)*tem
                else
                  adjsfculw(i) = adjsfculw3(i,2)
                endif
              endif
            else                                 ! all water
              adjsfculw(i) = adjsfculw3(i,3)
            endif
          enddo
        endif
!     if (lprnt) write(0,*)' kdt=',kdt,' tsfc=',Sfcprop%tsfc(ipr),' adjsfculw=',adjsfculw(ipr),&
!       ' adjsfculw3=',adjsfculw3(ipr,:),' icefr=',fice(ipr),' tsfc3=',tsfc3(ipr,:)
!
        do i=1,im
          Diag%dlwsfc(i) = Diag%dlwsfc(i) +   adjsfcdlw(i)*dtf
          Diag%ulwsfc(i) = Diag%ulwsfc(i) +   adjsfculw(i)*dtf
          Diag%psmean(i) = Diag%psmean(i) + Statein%pgr(i)*dtf        ! mean surface pressure
        enddo

        if (Model%ldiag3d) then
          if (Model%lsidea) then
            do k=1,levs
              do i=1,im
                Diag%dt3dt(i,k,1) = Diag%dt3dt(i,k,1) + Radtend%lwhd(i,k,1)*dtf
                Diag%dt3dt(i,k,2) = Diag%dt3dt(i,k,2) + Radtend%lwhd(i,k,2)*dtf
                Diag%dt3dt(i,k,3) = Diag%dt3dt(i,k,3) + Radtend%lwhd(i,k,3)*dtf
                Diag%dt3dt(i,k,4) = Diag%dt3dt(i,k,4) + Radtend%lwhd(i,k,4)*dtf
                Diag%dt3dt(i,k,5) = Diag%dt3dt(i,k,5) + Radtend%lwhd(i,k,5)*dtf
                Diag%dt3dt(i,k,6) = Diag%dt3dt(i,k,6) + Radtend%lwhd(i,k,6)*dtf
              enddo
            enddo
          else
            do k=1,levs
              do i=1,im
                Diag%dt3dt(i,k,1) = Diag%dt3dt(i,k,1) + Radtend%htrlw(i,k)*dtf
                Diag%dt3dt(i,k,2) = Diag%dt3dt(i,k,2) + Radtend%htrsw(i,k)*dtf*xmu(i)
              enddo
            enddo
          endif
        endif
      endif    ! end if_lssav_block

      do i=1,im
        kcnv(i)   = 0     !## CCPP ## GFS_typedefs.F90/interstitial_phys_reset
        kinver(i) = levs  !## CCPP ## GFS_typedefs.F90/interstitial_phys_reset
        invrsn(i) = .false.
        tx1(i)    = zero
        tx2(i)    = 10.0_kind_phys
        ctei_r(i) = 10.0_kind_phys
      enddo

!    Only used for old shallow convection with mstrat=.true.

      if ((((Model%imfshalcnv == 0 .and. Model%shal_cnv) .or. Model%old_monin)        &
                                   .and. Model%mstrat)   .or. Model%do_shoc) then
        ctei_rml(:) = Model%ctei_rm(1)*work1(:) + Model%ctei_rm(2)*work2(:)
        do k=1,levs/2
          do i=1,im
            if (Statein%prsi(i,1)-Statein%prsi(i,k+1) < 0.35_kind_phys*Statein%prsi(i,1)       &
                .and. (.not. invrsn(i))) then
              tem = (Statein%tgrs(i,k+1) - Statein%tgrs(i,k))  &
                  / (Statein%prsl(i,k)   - Statein%prsl(i,k+1))

              if (((tem > 0.00010_kind_phys) .and. (tx1(i) < zero)) .or.  &
                  ((tem-abs(tx1(i)) > zero) .and. (tx2(i) < zero))) then
                invrsn(i) = .true.

                if (Statein%qgrs(i,k,1) > Statein%qgrs(i,k+1,1)) then
                  tem1 = Statein%tgrs(i,k+1) + hocp*max(Statein%qgrs(i,k+1,1),qmin)
                  tem2 = Statein%tgrs(i,k)   + hocp*max(Statein%qgrs(i,k,1),qmin)

                  tem1 = tem1 / Statein%prslk(i,k+1) - tem2 / Statein%prslk(i,k)

!  --- ...  (cp/l)(deltathetae)/(deltatwater) > ctei_rm -> conditon for CTEI
                  ctei_r(i) = (one/hocp)*tem1/(Statein%qgrs(i,k+1,1)-Statein%qgrs(i,k,1)  &
                            + Statein%qgrs(i,k+1,ntcw)-Statein%qgrs(i,k,ntcw))
                else
                  ctei_r(i) = 10.0_kind_phys
                endif

                if ( ctei_rml(i) > ctei_r(i) ) then
                  kinver(i) = k
                else
                  kinver(i) = levs
                endif
              endif

              tx2(i) = tx1(i)
              tx1(i) = tem
            endif
          enddo
        enddo
      endif
!*## CCPP ##

!  --- ...  lu: initialize flag_guess, flag_iter, tsurf

!## CCPP ##* These initializations are done in GFS_typedefs.F90/interstitial_phys_reset except for as noted below
      do i=1,im
!       tsurf(i)        = Sfcprop%tsfc(i)
        flag_guess(i)   = .false.
        flag_iter(i)    = .true.
        drain(i)        = zero
        ep1d(i)         = zero
        gflx(i)         = zero
        runof(i)        = zero
        hflx(i)         = zero
        evap(i)         = zero
        evbs(i)         = zero
        evcw(i)         = zero
        trans(i)        = zero
        sbsno(i)        = zero
        snowc(i)        = zero
        snohf(i)        = zero
        !## CCPP ##* GFS_surface_generic.F90/GFS_surface_generic_pre_run
        Diag%zlvl(i)    = Statein%phil(i,1) * onebg
        Diag%smcwlt2(i) = zero
        Diag%smcref2(i) = zero
        wind(i)         = max(sqrt(Statein%ugrs(i,1)*Statein%ugrs(i,1) + &
                                   Statein%vgrs(i,1)*Statein%vgrs(i,1))  &
                        + max(zero, min(Tbd%phy_f2d(i,Model%num_p2d), 30.0_kind_phys)), one)
        !*## CCPP ##
      enddo
!*## CCPP ##

!  --- ...  lu: iter-loop over (sfc_diff,sfc_drv,sfc_ocean,sfc_sice)
!## CCPP ##* This loop is implemented using the subcycle/iteration capability in the CCPP SDF
      do iter=1,2
!*## CCPP ##

!  --- ...  surface exchange coefficients
!
!     if (lprnt) write(0,*)' tsfc=',Sfcprop%tsfc(ipr),'iter=', &
!           iter ,'wet=',wet(ipr),'dry=',dry(ipr),' icy=',icy(ipr),&
!           ' zorl3=',zorl3(ipr,:),' uustar3=',uustar3(ipr,:)

!## CCPP ##* sfc_diff.f/sfc_diff_run
        call sfc_diff                                                   &
!  ---  inputs:
          (im, Statein%pgr,                                             &
           Statein%tgrs(:,1), Statein%qgrs(:,1,1), Diag%zlvl, wind,     &
           Statein%prsl(:,1), work3,                                    &
           sigmaf, vegtype,   Sfcprop%shdmax, Model%ivegsrc,            &
           z01d, zt1d,                                                  & ! mg, sfc-perts
           flag_iter, Model%redrag,                                     &
           Diag%u10m,    Diag%v10m,  Model%sfc_z0_type,                 &
           wet, dry, icy, tsfc3, tsurf3, snowd3,                        &
!  ---  input/output:
           zorl3, Sfcprop%zorlw, uustar3,                               &
!  ---  outputs:
           cd3, cdq3, rb3, stress3, ffmm3, ffhh3, fm103, fh23)
!          cd3, cdq3, rb3, stress3, ffmm3, ffhh3, fm103, fh23, wind, lprnt, ipr)
!
!     if (lprnt) write(0,*)' aft sfc_diff cd3=',cd3(ipr,:),' cdq3=',cdq3(ipr,:),'iter=', iter, &
!           ' zorl3=',zorl3(ipr,:),' uustar3=',uustar3(ipr,:)
!  --- ...  lu: update flag_guess
!*## CCPP ##
!## CCPP ##* GFS_surface_loop_control/GFS_surface_loop_control_part1_run
        do i=1,im
          if (iter == 1 .and. wind(i) < 2.0_kind_phys) then
            flag_guess(i) = .true.
          endif
        enddo
!*## CCPP ##
!## CCPP ##* sfc_nst.f/sfc_nst_pre_run Note: the conditional is not included in the CCPP scheme, so calling 
! this code is controlled by its presence in the active CCPP SDF
        if (Model%nstf_name(1) > 0) then
          do i=1,im
            if (wet(i)) then
!             tem         = (Sfcprop%oro(i)-Sfcprop%oro_uf(i)) * rlapse
              tem         = zero
              tseal(i)    = tsfc3(i,3)  + tem
              tsurf3(i,3) = tsurf3(i,3) + tem
            endif
          enddo
          if (Model%cplflx) then       ! apply only at ocean points
            call get_dtzm_2d (Sfcprop%xt,  Sfcprop%xz, Sfcprop%dt_cool,  &
                              Sfcprop%z_c, wet, zero, omz1, im, 1, dtzm)
            do i=1,im
              if (wet(i) .and. Sfcprop%oceanfrac(i) > zero) then
                Sfcprop%tref(i) = Sfcprop%tsfco(i) - dtzm(i)    ! update Tf with T1 and NSST T-Profile
                if (abs(Sfcprop%xz(i)) > zero) then
                  tem2 = one / Sfcprop%xz(i)
                else
                  tem2 = zero
                endif
                tseal(i)    = Sfcprop%tref(i) + (Sfcprop%xt(i)+Sfcprop%xt(i)) * tem2 &
                                              - Sfcprop%dt_cool(i)
                tsurf3(i,3) = tseal(i)
              endif
            enddo
          endif

!     if (lprnt) write(0,*)' bef nst tseal=',tseal(ipr) &
!     ,' tsfc3=',tsfc3(ipr,3),' tsurf3=',tsurf3(ipr,3), &
!           iter ,'wet=',wet(ipr),'dry=',dry(ipr),' icy=',icy(ipr),&
!      ' tref=',Sfcprop%tref(ipr),' tgrs=',Statein%tgrs(ipr,1),' qgrs=',Statein%qgrs(ipr,1,1), &
!      ' prsl=',Statein%prsl(ipr,1),' cd3=',cd3(ipr,3),' cdq3=',cdq3(ipr,3),' work3=',    &
!      work3(ipr),' semis3=',semis3(ipr,3),' gabsbdlw3=',gabsbdlw3(ipr,3),' adjsfcnsw=',    &
!      adjsfcnsw(ipr),' wind=',wind(ipr),' tseal=',tseal(ipr),' xcosz=',xcosz(ipr)
!*## CCPP ##
!## CCPP ##* sfc_nst.f/sfc_nst_run
          call sfc_nst                                                  &
!  ---  inputs:
            (im, Statein%pgr, Statein%ugrs(:,1), Statein%vgrs(:,1),     &
             Statein%tgrs(:,1), Statein%qgrs(:,1,1),                    &
             Sfcprop%tref, cd3(:,3), cdq3(:,3), Statein%prsl(:,1),      &
             work3, wet, Grid%xlon, Grid%sinlat, stress3(:,3),          &
             semis3(:,3), gabsbdlw3(:,3), adjsfcnsw, tprcp3(:,3),       &
             dtf, kdt, Model%solhr, xcosz,                              &
             wind, flag_iter,                                           &
             flag_guess, Model%nstf_name, lprnt, ipr,                   &
!  ---  input/output
             tseal, tsurf3(:,3), Sfcprop%xt, Sfcprop%xs,                &
             Sfcprop%xu,  Sfcprop%xv,   Sfcprop%xz, Sfcprop%zm,         &
             Sfcprop%xtts,Sfcprop%xzts, Sfcprop%dt_cool,                &
             Sfcprop%z_c, Sfcprop%c_0,  Sfcprop%c_d,                    &
             Sfcprop%w_0, Sfcprop%w_d,  Sfcprop%d_conv,                 &
             Sfcprop%ifd, Sfcprop%qrain,                                &
!  ---  outputs:
             qss3(:,3),  gflx3(:,3), cmm3(:,3), chh3(:,3), evap3(:,3),  &
             hflx3(:,3), ep1d3(:,3))
!*## CCPP ##
!## CCPP ##* sfc_nst.f/sfc_nst_post_run

!         do i=1,im
!!          if (wet(i) .and. .not.icy(i)) then
!!          if (wet(i) .and. (Model%frac_grid .or. .not. icy(i))) then
!           if (wet(i)) then
!             tsurf3(i,3) = tsurf3(i,3)                                    &
!                         - (Sfcprop%oro(i)-Sfcprop%oro_uf(i)) * rlapse
!           endif
!         enddo

!  --- ...  run nsst model  ... ---

          if (Model%nstf_name(1) > 1) then
            zsea1 = 0.001_kind_phys*real(Model%nstf_name(4))
            zsea2 = 0.001_kind_phys*real(Model%nstf_name(5))
            call get_dtzm_2d (Sfcprop%xt,  Sfcprop%xz, Sfcprop%dt_cool, &
                              Sfcprop%z_c, wet, zsea1, zsea2, im, 1, dtzm)
            do i=1,im
!             if (wet(i) .and. .not.icy(i)) then
!             if (wet(i) .and. (Model%frac_grid .or. .not. icy(i))) then
              if (wet(i)) then
                tsfc3(i,3) = max(tgice,Sfcprop%tref(i) + dtzm(i))
!               tsfc3(i,3) = max(271.2,Sfcprop%tref(i) + dtzm(i)) -  &
!                               (Sfcprop%oro(i)-Sfcprop%oro_uf(i))*rlapse
              endif
            enddo
          endif

!     if (lprnt) write(0,*)' aft nst tref=',Sfcprop%tref(ipr) &
!     ,' tsfc3=',tsfc3(ipr,3),' dtzm=',dtzm(ipr),' hflx33=',hflx3(ipr,3)
!*## CCPP ##
!         if (lprnt) print *,' tseaz2=',Sfcprop%tsfc(ipr),' tref=',tref(ipr),   &
!    &    ' dt_cool=',dt_cool(ipr),' dt_warm=',dt_warm(ipr),' kdt=',kdt

!## CCPP ## Note: This conditional is replaced by whether the sfc_ocean scheme is in the CCPP SDF
        else

!  --- ...  surface energy balance over ocean
!## CCPP ##* sfc_ocean.F/sfc_ocean_run
          call sfc_ocean                                                &
!  ---  inputs:
           (im, Statein%pgr,                                            &
            Statein%tgrs(:,1), Statein%qgrs(:,1,1), tsfc3(:,3),         &
            cd3(:,3), cdq3(:,3), Statein%prsl(:,1), work3, wet,         &
            wind, flag_iter,                                            &
!  ---  outputs:
            qss3(:,3), cmm3(:,3), chh3(:,3), gflx3(:,3), evap3(:,3),    &
            hflx3(:,3), ep1d3(:,3))
!*## CCPP ##

        endif       ! if nstf_name(1) > 0

!       if (lprnt) write(0,*)' sfalb=',Radtend%sfalb(ipr),' ipr=',ipr   &
!     ,   ' weasd=',Sfcprop%weasd(ipr)                                  &
!     ,   ' tprcp=',Sfcprop%tprcp(ipr),' kdt=',kdt,' iter=',iter        &
!     ,' tseabefland=',Sfcprop%tsfc(ipr)

!  --- ...  surface energy balance over land
!
!## CCPP ##* Note: the conditional is not included in the CCPP, so calling 
! the LSM scheme is controlled by its presence in the active CCPP SDF
        if (Model%lsm == Model%lsm_noah) then                          ! noah lsm call
!*## CCPP ##

!     if (lprnt) write(0,*)' tseal=',tseal(ipr),' tsurf=',tsurf(ipr),iter &
!     ,' stsoil0=',stsoil(ipr,:)
!    &,' pgr=',pgr(ipr),' sfcemis=',sfcemis(ipr)

!## CCPP ##* sfc_drv.f/lsm_noah_run
          call sfc_drv                                                   &
!  ---  inputs:
           (im, lsoil, Statein%pgr,                                      &
            Statein%tgrs(:,1), Statein%qgrs(:,1,1), soiltyp, vegtype,    &
            sigmaf, semis3(:,1), gabsbdlw3(:,1), adjsfcdsw, adjsfcnsw, dtf,&
!           sigmaf, Radtend%semis, gabsbdlw, adjsfcdsw, adjsfcnsw, dtf,  &
            Sfcprop%tg3, cd3(:,1), cdq3(:,1), Statein%prsl(:,1), work3,  &
            Diag%zlvl, dry, wind, slopetyp,                              &
            Sfcprop%shdmin, Sfcprop%shdmax, Sfcprop%snoalb,              &
            Radtend%sfalb, flag_iter, flag_guess, Model%lheatstrg,       &
            Model%isot, Model%ivegsrc,                                   &
            bexp1d, xlai1d, vegf1d, Model%pertvegf,                      &
!  ---  input/output:
            weasd3(:,1), snowd3(:,1), tsfc3(:,1), tprcp3(:,1),           &
            Sfcprop%srflag, smsoil, stsoil, slsoil, Sfcprop%canopy,      &
            trans, tsurf3(:,1), zorl3(:,1),                              &
!  ---  outputs:
            Sfcprop%sncovr, qss3(:,1), gflx3(:,1), drain, evap3(:,1),    &
            hflx3(:,1), ep1d3(:,1), runof,                               &
            cmm3(:,1),  chh3(:,1), evbs, evcw, sbsno, snowc, Diag%soilm, &
            snohf, Diag%smcwlt2, Diag%smcref2, Diag%wet1)
!*## CCPP ##

!     if (lprnt) write(0,*)' tseae=',tseal(ipr),' tsurf=',tsurf(ipr),iter
!                         ,' phy_f2d=',phy_f2d(ipr,num_p2d)

!     if (lprnt) write(0,*)' hflx3=',hflx3(ipr,:),' evap3=',evap3(ipr,:)

!## CCPP ##* sfc_noahmp_drv.f/noahmpdrv_run
! Noah MP call
!
       elseif (Model%lsm == Model%lsm_noahmp) then
          call noahmpdrv                                               &
!  ---  inputs:
           (im, lsoil,kdt, Statein%pgr,  Statein%ugrs, Statein%vgrs,   &
            Statein%tgrs,  Statein%qgrs, soiltyp, vegtype, sigmaf,     &
            semis3(:,1),   gabsbdlw3(:,1), adjsfcdsw, adjsfcnsw, dtf,  &
!           Radtend%semis, gabsbdlw,     adjsfcdsw,  adjsfcnsw, dtf,   &
            Sfcprop%tg3, cd3(:,1), cdq3(:,1), Statein%prsl(:,1), work3,&
            Diag%zlvl, dry,   wind, slopetyp,                          &
            Sfcprop%shdmin,   Sfcprop%shdmax,  Sfcprop%snoalb,         &
            Radtend%sfalb,    flag_iter,       flag_guess,             &
            Model%iopt_dveg,  Model%iopt_crs,  Model%iopt_btr,         &
            Model%iopt_run,   Model%iopt_sfc,  Model%iopt_frz,         &
            Model%iopt_inf,   Model%iopt_rad,  Model%iopt_alb,         &
            Model%iopt_snf,   Model%iopt_tbot, Model%iopt_stc,         &
            grid%xlat, xcosz, Model%yearlen,   Model%julian, Model%imn,&
            Sfcprop%drainncprv, Sfcprop%draincprv, Sfcprop%dsnowprv,   &
            Sfcprop%dgraupelprv, Sfcprop%diceprv,                      &
!  ---  in/outs:
            weasd3(:,1), snowd3(:,1), tsfc3(:,1), tprcp3(:,1),         &
            Sfcprop%srflag, smsoil, stsoil, slsoil, Sfcprop%canopy,    &
            trans, tsurf3(:,1), zorl3(:,1),                            &
!
            Sfcprop%snowxy,   Sfcprop%tvxy,    Sfcprop%tgxy,  Sfcprop%canicexy, &
            Sfcprop%canliqxy, Sfcprop%eahxy,   Sfcprop%tahxy, Sfcprop%cmxy,     &
            Sfcprop%chxy,     Sfcprop%fwetxy,  Sfcprop%sneqvoxy,                &
            Sfcprop%alboldxy, Sfcprop%qsnowxy, Sfcprop%wslakexy,                &
            Sfcprop%zwtxy,    Sfcprop%waxy,    Sfcprop%wtxy, Sfcprop%tsnoxy,    &
            Sfcprop%zsnsoxy,  Sfcprop%snicexy, Sfcprop%snliqxy,                 &
            Sfcprop%lfmassxy, Sfcprop%rtmassxy,                                 &
            Sfcprop%stmassxy, Sfcprop%woodxy,  Sfcprop%stblcpxy,                &
            Sfcprop%fastcpxy, Sfcprop%xlaixy,  Sfcprop%xsaixy,                  &
            Sfcprop%taussxy,  Sfcprop%smoiseq, Sfcprop%smcwtdxy,                &
            Sfcprop%deeprechxy, Sfcprop%rechxy,                                 &
!  ---  outputs:
            Sfcprop%sncovr, qss3(:,1), gflx3(:,1), drain, evap3(:,1),           &
            hflx3(:,1), ep1d3(:,1), runof,                                      &
            cmm3(:,1), chh3(:,1), evbs, evcw, sbsno, snowc, Diag%soilm,         &
            snohf, Diag%smcwlt2, Diag%smcref2, Diag%wet1, t2mmp, q2mp)

!     if (lprnt) write(0,*)' tseae=',tsea(ipr),' tsurf=',tsurf(ipr),iter &
!    &,' phy_f2d=',phy_f2d(ipr,num_p2d)
!*## CCPP ##

        elseif (Model%lsm == Model%lsm_ruc) then
           write (0,*) 'RUC LSM is available only in CCPP'
           stop

        endif !lsm

        !! Strictly speaking, this is not required. But when
        !! hunting for bit-for-bit differences, updating the
        !! subsurface variables in the Sfcprop DDT makes
        !! life a lot easier
        !if (Model%frac_grid) then
        !  do k=1,lsoil
        !    do i=1,im
        !      if (dry(i)) then
        !        Sfcprop%smc(i,k) = smsoil(i,k)
        !        Sfcprop%stc(i,k) = stsoil(i,k)
        !        Sfcprop%slc(i,k) = slsoil(i,k)
        !      endif
        !    enddo
        !  enddo
        !else
        !  do k=1,lsoil
        !    do i=1,im
        !      Sfcprop%smc(i,k) = smsoil(i,k)
        !      Sfcprop%stc(i,k) = stsoil(i,k)
        !      Sfcprop%slc(i,k) = slsoil(i,k)
        !    enddo
        !  enddo
        !endif

!       if (lprnt) write(0,*)' tseabeficemodel =',Sfcprop%tsfc(ipr),' me=',me   &
!    &,   ' kdt=',kdt,' tsfc32=',tsfc3(ipr,2),' fice=',fice(ipr)                &
!    &,' stsoil=',stsoil(ipr,:),' tsfc33=',tsfc3(ipr,3),' islmsk=',islmsk(ipr)

!  --- ...  surface energy balance over seaice
!## CCPP ##* sfc_sice.f/sfc_sice_run (local adjustment to avoid resetting islmsk after call to sfc_sice_run)
        if (Model%cplflx) then
          do i=1,im
            if (flag_cice(i)) then
               islmsk(i) = islmsk_cice(i)
            endif
          enddo
!*## CCPP ##

!## CCPP ##* sfc_cice.f/sfc_cice_run
! call sfc_cice for sea ice points in the coupled model (i.e. islmsk=4)
!
          call sfc_cice                                                  &
!  ---  inputs:
           (im, Statein%tgrs(:,1),                                       &
            Statein%qgrs(:,1,1),  cd3(:,2), cdq3(:,2),                   &
            Statein%prsl(:,1),    wind,                                  &
            flag_cice, flag_iter,                                        &
            Coupling%dqsfcin_cpl, Coupling%dtsfcin_cpl,                  &
            Coupling%dusfcin_cpl, Coupling%dvsfcin_cpl,                  &
            Coupling%hsnoin_cpl,                                         &
!  ---  outputs:
            qss3(:,2), cmm3(:,2), chh3(:,2), evap3(:,2), hflx3(:,2),     &
            stress3(:,2), weasd3(:,2), snowd3(:,2), ep1d3(:,2))
        endif
!*## CCPP ##

!
! call sfc_sice for lake ice and for the uncoupled case, sea ice (i.e. islmsk=2)
!
        if (Model%frac_grid) then
          do i=1,im
            if (icy(i) .and. islmsk(i) < 2) then
              if (Sfcprop%oceanfrac(i) > zero) then
                tem = Model%min_seaice
              else
                tem = Model%min_lakeice
              endif
              if (fice(i) > tem) then
                islmsk(i) = 2
                tsfc3(i,2) = Sfcprop%tisfc(i)
              endif
            endif
          enddo
        endif
!## CCPP ##* sfc_sice.f/sfc_sice_run
        call sfc_sice                                                            &
!  ---  inputs:
           (im, lsoil, Statein%pgr,                                              &
            Statein%tgrs(:,1), Statein%qgrs(:,1,1), dtf, semis3(:,2),            &
!           Statein%tgrs(:,1), Statein%qgrs(:,1,1), dtf, Radtend%semis,          &
            gabsbdlw3(:,2), adjsfcnsw, adjsfcdsw, Sfcprop%srflag,                &
            cd3(:,2), cdq3(:,2),                                                 &
            Statein%prsl(:,1), work3, islmsk, wind,                              &
            flag_iter, lprnt, ipr, Model%min_lakeice,                            &
!  ---  input/output:
            zice, fice, tice, weasd3(:,2), tsfc3(:,2), tprcp3(:,2),              &
            stsoil, ep1d3(:,2),                                                  &
!  ---  outputs:
            snowd3(:,2), qss3(:,2), snowmt, gflx3(:,2), cmm3(:,2), chh3(:,2),    &
            evap3(:,2),  hflx3(:,2))
!*## CCPP ##
!## CCPP ##* This section is not needed for CCPP.
        if (Model%frac_grid) then
          do i = 1, im
            if (islmsk(i) == 2 .and. fice(i) < one) then
              wet(i) = .true.
              tsfc3(i,3) = max(Sfcprop%tisfc(i), tgice)
            endif
          enddo
        endif
        if (Model%cplflx) then
          do i = 1, im
            if (flag_cice(i)) then
               islmsk(i) = nint(Sfcprop%slmsk(i))
            endif
          enddo
        endif
!*## CCPP ##

!       if (lprnt) write(0,*)' tseaafticemodel =',tsfc3(ipr,:),' me=',me &
!    &,   ' kdt=',kdt,' iter=',iter,' fice=',fice(ipr),' wet=',wet(ipr),' icy=',icy(ipr)&
!    &,' dry=',dry(ipr)

!  --- ...  lu: update flag_iter and flag_guess
!## CCPP ##* GFS_surface_loop_control.F90/GFS_surface_loop_control_part_2
        do i=1,im
          flag_iter(i)  = .false.
          flag_guess(i) = .false.

          if (iter == 1 .and. wind(i) < 2.0_kind_phys) then
!           if (dry(i) .or. (wet(i) .and. .not.icy(i)                   &
            if (dry(i) .or. (wet(i) .and. Model%nstf_name(1) > 0)) then
              flag_iter(i) = .true.
            endif
          endif

        enddo
!*## CCPP ##

      enddo   ! end iter_loop


! --- generate ocean/land/ice composites

!## CCPP ##* GFS_surface_compoistes.F90/GFS_surface_composites_post_run
      if (Model%frac_grid) then
        do i=1, im
!
! Three-way composites (fields from sfc_diff)
          txl = frland(i)
          txi = fice(i)*(one - frland(i)) ! txi = ice fraction wrt whole cell
          txo = max(zero, one - txl - txi)

!      if (i == ipr .and. lprnt) write(0,*)' txl=',txl,' fice=',fice(i),' txi=',txi,&
!         ' txo=',txo,' dry=',dry(i),' wet=',wet(i),' icy=',icy(i),' oceanfrac=',&
!           Sfcprop%oceanfrac(i),' frland=',frland(i)

          Sfcprop%zorl(i)   = txl*zorl3(i,1)   + txi*zorl3(i,2)   + txo*zorl3(i,3)
          cd(i)             = txl*cd3(i,1)     + txi*cd3(i,2)     + txo*cd3(i,3)
          cdq(i)            = txl*cdq3(i,1)    + txi*cdq3(i,2)    + txo*cdq3(i,3)
          rb(i)             = txl*rb3(i,1)     + txi*rb3(i,2)     + txo*rb3(i,3)
          stress(i)         = txl*stress3(i,1) + txi*stress3(i,2) + txo*stress3(i,3)
          Sfcprop%ffmm(i)   = txl*ffmm3(i,1)   + txi*ffmm3(i,2)   + txo*ffmm3(i,3)
          Sfcprop%ffhh(i)   = txl*ffhh3(i,1)   + txi*ffhh3(i,2)   + txo*ffhh3(i,3)
          Sfcprop%uustar(i) = txl*uustar3(i,1) + txi*uustar3(i,2) + txo*uustar3(i,3)
          fm10(i)           = txl*fm103(i,1)   + txi*fm103(i,2)   + txo*fm103(i,3)
          fh2(i)            = txl*fh23(i,1)    + txi*fh23(i,2)    + txo*fh23(i,3)
!         tsurf(i)          = txl*tsurf3(i,1)  + txi*tice(i)      + txo*tsurf3(i,3)
!         tsurf(i)          = txl*tsurf3(i,1)  + txi*tsurf3(i,2)  + txo*tsurf3(i,3)  ! not used again! Moorthi
          Diag%cmm(i)       = txl*cmm3(i,1)    + txi*cmm3(i,2)    + txo*cmm3(i,3)
          Diag%chh(i)       = txl*chh3(i,1)    + txi*chh3(i,2)    + txo*chh3(i,3)
!         gflx(i)           = txl*gflx3(i,1)   + txi*gflx3(i,2)   + txo*gflx3(i,3)
          ep1d(i)           = txl*ep1d3(i,1)   + txi*ep1d3(i,2)   + txo*ep1d3(i,3)
!         Sfcprop%weasd(i)  = txl*weasd3(i,1)  + txi*weasd3(i,2)  + txo*weasd3(i,3)
!         Sfcprop%snowd(i)  = txl*snowd3(i,1)  + txi*snowd3(i,2)  + txo*snowd3(i,3)
          Sfcprop%weasd(i)  = txl*weasd3(i,1)  + txi*weasd3(i,2)
          Sfcprop%snowd(i)  = txl*snowd3(i,1)  + txi*snowd3(i,2)
!         Sfcprop%tprcp(i)  = txl*tprcp3(i,1)  + txi*tprcp3(i,2)  + txo*tprcp3(i,3)

          if (.not. flag_cice(i) .and. islmsk(i) == 2) then
            tem             = one - txl
            evap(i)         = txl*evap3(i,1)   + tem*evap3(i,2)
            hflx(i)         = txl*hflx3(i,1)   + tem*hflx3(i,2)
            qss(i)          = txl*qss3(i,1)    + tem*qss3(i,2)
            gflx(i)         = txl*gflx3(i,1)   + tem*gflx3(i,2)
          else
            evap(i)         = txl*evap3(i,1)   + txi*evap3(i,2)   + txo*evap3(i,3)
            hflx(i)         = txl*hflx3(i,1)   + txi*hflx3(i,2)   + txo*hflx3(i,3)
            qss(i)          = txl*qss3(i,1)    + txi*qss3(i,2)    + txo*qss3(i,3)
            gflx(i)         = txl*gflx3(i,1)   + txi*gflx3(i,2)   + txo*gflx3(i,3)
          endif
          Sfcprop%tsfc(i)   = txl*tsfc3(i,1)   + txi*tice(i)      + txo*tsfc3(i,3)
!         Sfcprop%tsfc(i)   = txl*tsfc3(i,1)   + txi*tsfc3(i,2)   + txo*tsfc3(i,3)

!     if (i == ipr .and. lprnt) then
!       write(0,*)' tsfc=',Sfcprop%tsfc(i),' txl=',txl,' txi=',txi,' txo=',txo, &
!        ' tsfc3=',tsfc3(i,:),' evap3=',evap3(i,:),' evap=',evap(i),' tice=',tice(i),&
!        'Sfcprop%zorl=',Sfcprop%zorl(ipr)
!     endif

!         Diag%cmm(i)       = txl*cmm3(i,1)    + txi*cmm3(i,2)    + txo*cmm3(i,3)
!         Diag%chh(i)       = txl*chh3(i,1)    + txi*chh3(i,2)    + txo*chh3(i,3)

          Sfcprop%zorll(i) = zorl3(i,1)
          Sfcprop%zorli(i) = zorl3(i,2)
          Sfcprop%zorlo(i) = zorl3(i,3)

          if (dry(i)) then
            Sfcprop%tsfcl(i) = tsfc3(i,1)      ! over land
          elseif (wet(i)) then
            Sfcprop%tsfcl(i) = tsfc3(i,3)      ! over land
          else
            Sfcprop%tsfcl(i) = tice(i)         ! over land
          endif
          if (wet(i)) then
            Sfcprop%tsfco(i) = tsfc3(i,3)      ! over lake or ocean when uncoupled
          elseif (icy(i)) then
            Sfcprop%tsfco(i) = tice(i)         ! over lake or ocean when uncoupled
          else
            Sfcprop%tsfco(i) = tsfc3(i,1)      ! over lake or ocean when uncoupled
          endif
          if (icy(i)) then
            Sfcprop%tisfc(i) = tice(i)         ! over lake or ocean when uncoupled
!           if (Sfcprop%zorll(i) > 1000.0) Sfcprop%zorll(i) = zorl3(i,2)
          elseif (wet(i)) then
            Sfcprop%tisfc(i) = tsfc3(i,3)       ! over lake or ocean when uncoupled
          else
            Sfcprop%tisfc(i) = tsfc3(i,1)       ! over lake or ocean when uncoupled
          endif
                                                         ! for coupled model ocean will replace this
!         if (icy(i)) Sfcprop%tisfc(i) = tsfc3(i,2)      ! over ice when uncoupled
!         if (icy(i)) Sfcprop%tisfc(i) = tice(i)         ! over ice when uncoupled

!         if (wet(i) .and. .not. Model%cplflx) then
!           Sfcprop%tsfco(i) = tsfc3(i,3)                ! over lake or ocean when uncoupled
!           Sfcprop%tisfc(i) = tsfc3(i,2)                ! over ice when uncoupled
!         endif

          if (.not. flag_cice(i)) then
!           if (islmsk(i) == 2) then                     ! return updated lake ice thickness & concentration to global array
            if (icy(i)) then                             ! return updated lake ice thickness & concentration to global array
              Sfcprop%hice(i)  = zice(i)
              Sfcprop%fice(i)  = fice(i) 
              Sfcprop%tisfc(i) = tice(i)
            else                                         ! this would be over open ocean or land (no ice fraction)
              Sfcprop%hice(i)  = zero
              Sfcprop%fice(i)  = zero
              Sfcprop%tisfc(i) = Sfcprop%tsfc(i)
            endif
          endif
        enddo
      else
        do i=1,im
          if (flag_cice(i) .and. wet(i) .and. fice(i) < Model%min_seaice) then 
              islmsk(i) = 0
              fice(i) = zero
           endif
          if (islmsk(i) == 1) then
            k = 1
            Sfcprop%tsfcl(i) = tsfc3(i,1)  ! over land
            stress(i)        = stress3(i,1)
!           Sfcprop%tprcp(i) = tprcp3(i,1)
            Sfcprop%tsfco(i) = tsfc3(i,1)
            Sfcprop%tisfc(i) = tsfc3(i,1)
          elseif (islmsk(i) == 0) then
            k = 3
            Sfcprop%tsfco(i) = tsfc3(i,3)  ! over lake (and ocean when uncoupled)
            stress(i)        = stress3(i,3)
!           Sfcprop%tprcp(i) = tprcp3(i,3)
            Sfcprop%tisfc(i) = tsfc3(i,3)
            Sfcprop%tsfcl(i) = tsfc3(i,3)
          else
            k = 2
            stress(i)        = stress3(i,2)
!           Sfcprop%tprcp(i) = fice(i)*tprcp3(i,2)  + (one-fice(i))*tprcp3(i,3)
          endif
          Sfcprop%zorl(i)   = zorl3(i,k)
          cd(i)             = cd3(i,k)
          cdq(i)            = cdq3(i,k)
          rb(i)             = rb3(i,k)
          Sfcprop%ffmm(i)   = ffmm3(i,k)
          Sfcprop%ffhh(i)   = ffhh3(i,k)
          Sfcprop%uustar(i) = uustar3(i,k)
          fm10(i)           = fm103(i,k)
          fh2(i)            = fh23(i,k)
!         tsurf(i)          = tsurf3(i,k)
          Diag%cmm(i)       = cmm3(i,k)
          Diag%chh(i)       = chh3(i,k)
          gflx(i)           = gflx3(i,k)
          ep1d(i)           = ep1d3(i,k)
          Sfcprop%weasd(i)  = weasd3(i,k)
          Sfcprop%snowd(i)  = snowd3(i,k)
          evap(i)           = evap3(i,k)
          hflx(i)           = hflx3(i,k)
          qss(i)            = qss3(i,k)
          Sfcprop%tsfc(i)   = tsfc3(i,k)

          Sfcprop%zorll(i)  = zorl3(i,1)
          Sfcprop%zorli(i)  = zorl3(i,2)
          Sfcprop%zorlo(i)  = zorl3(i,3)

          if (flag_cice(i)) then
            if (wet(i) .and. fice(i) > Model%min_seaice) then  ! this was already done for lake ice in sfc_sice
              txi = fice(i)
              txo = one - txi
              evap(i)         = txi * evap3(i,2)   + txo * evap3(i,3)
              hflx(i)         = txi * hflx3(i,2)   + txo * hflx3(i,3)
              Sfcprop%tsfc(i) = txi * tsfc3(i,2)   + txo * tsfc3(i,3)
              stress(i)       = txi  *stress3(i,2) + txo * stress3(i,3)
              qss(i)          = txi * qss3(i,2)    + txo * qss3(i,3)
              ep1d(i)         = txi * ep1d3(i,2)   + txo * ep1d3(i,3)
              Sfcprop%zorl(i) = txi*zorl3(i,2)     + txo*zorl3(i,3)
            endif
          elseif (islmsk(i) == 2) then  ! return updated lake ice thickness & concentration to global array
            Sfcprop%tisfc(i) = tice(i)  ! over lake ice (and sea ice when uncoupled)
            Sfcprop%hice(i)  = zice(i)
            Sfcprop%fice(i)  = fice(i)  ! fice is fraction of lake area that is frozen
            Sfcprop%zorl(i)  = fice(i)*zorl3(i,2) + (one-fice(i))*zorl3(i,3)
          else                          ! this would be over open ocean or land (no ice fraction)
            Sfcprop%hice(i)  = zero
            Sfcprop%fice(i)  = zero
            Sfcprop%tisfc(i) = Sfcprop%tsfc(i)
            icy(i)           = .false.
          endif
          Sfcprop%tsfcl(i) = Sfcprop%tsfc(i)
          if (wet(i)) then
            Sfcprop%tsfco(i) = tsfc3(i,3)
          else
            Sfcprop%tsfco(i) =Sfcprop%tsfc(i)
          endif
          do k=1,Model%kice ! store tiice in stc to reduce output in the nonfrac grid case
            Sfcprop%stc(i,k) = Sfcprop%tiice(i,k)
          enddo
        enddo
      endif       ! if (Model%frac_grid)
!*## CCPP ##

! --- compositing done

!     if (lprnt) write(0,*) 'tisfc=',Sfcprop%tisfc(ipr),'tice=',tice(ipr),' kdt=',kdt

      do i=1,im
!## CCPP ##* GFS_surface_generic.F90/GFS_surface_generic_post_run
        Diag%epi(i)     = ep1d(i)
!*## CCPP ##
        Diag%dlwsfci(i) = adjsfcdlw(i)
        Diag%ulwsfci(i) = adjsfculw(i)
!## CCPP ##* GFS_surface_composites.F90/GFS_surface_composites_inter_run
        Diag%uswsfci(i) = adjsfcdsw(i) - adjsfcnsw(i)
!*## CCPP ##
        Diag%dswsfci(i) = adjsfcdsw(i)
!## CCPP ##* GFS_surface_generic.F90/GFS_surface_generic_post_run
        Diag%gfluxi(i)  = gflx(i)
        Diag%t1(i)      = Statein%tgrs(i,1)
        Diag%q1(i)      = Statein%qgrs(i,1,1)
        Diag%u1(i)      = Statein%ugrs(i,1)
        Diag%v1(i)      = Statein%vgrs(i,1)
!*## CCPP ##
      enddo

!  --- ...  update near surface fields

!## CCPP ##* sfc_diag.f/sfc_diag_run
      call sfc_diag (im, Statein%pgr, Statein%ugrs(:,1), Statein%vgrs(:,1),       &
                     Statein%tgrs(:,1), Statein%qgrs(:,1,1), work3, evap,         &
                     Sfcprop%ffmm, Sfcprop%ffhh, fm10, fh2, Sfcprop%tsfc, qss,    &
                     Sfcprop%f10m, Diag%u10m, Diag%v10m, Sfcprop%t2m, Sfcprop%q2m)
!*## CCPP ##

!## CCPP ##* This block is not in the CCPP
      Tbd%phy_f2d(:,Model%num_p2d) = zero
!*## CCPP ##

      if (Model%lsm == Model%lsm_noahmp) then
        do i=1,im
         if (dry(i)) then
          Sfcprop%t2m(i) = t2mmp(i)
          Sfcprop%q2m(i) = q2mp(i)
         endif
        enddo
      endif ! if (Model%lsm == Model%lsm_noahmp)
! *DH

      if (Model%cplflx .or. Model%cplwav) then
        do i=1,im
          Coupling%u10mi_cpl   (i) = Diag%u10m(i)
          Coupling%v10mi_cpl   (i) = Diag%v10m(i)
        enddo
      endif

      if (Model%cplflx) then
        do i=1,im
          Coupling%dlwsfci_cpl (i) = adjsfcdlw(i)
          Coupling%dswsfci_cpl (i) = adjsfcdsw(i)
          Coupling%dlwsfc_cpl  (i) = Coupling%dlwsfc_cpl(i) + adjsfcdlw(i)*dtf
          Coupling%dswsfc_cpl  (i) = Coupling%dswsfc_cpl(i) + adjsfcdsw(i)*dtf
          Coupling%dnirbmi_cpl (i) = adjnirbmd(i)
          Coupling%dnirdfi_cpl (i) = adjnirdfd(i)
          Coupling%dvisbmi_cpl (i) = adjvisbmd(i)
          Coupling%dvisdfi_cpl (i) = adjvisdfd(i)
          Coupling%dnirbm_cpl  (i) = Coupling%dnirbm_cpl(i) + adjnirbmd(i)*dtf
          Coupling%dnirdf_cpl  (i) = Coupling%dnirdf_cpl(i) + adjnirdfd(i)*dtf
          Coupling%dvisbm_cpl  (i) = Coupling%dvisbm_cpl(i) + adjvisbmd(i)*dtf
          Coupling%dvisdf_cpl  (i) = Coupling%dvisdf_cpl(i) + adjvisdfd(i)*dtf
          Coupling%nlwsfci_cpl (i) = adjsfcdlw(i)           - adjsfculw(i)
          if (wet(i)) then
            Coupling%nlwsfci_cpl(i) = adjsfcdlw(i)          - adjsfculw3(i,3)
          endif
          Coupling%nlwsfc_cpl  (i) = Coupling%nlwsfc_cpl(i) + Coupling%nlwsfci_cpl(i)*dtf
          Coupling%t2mi_cpl    (i) = Sfcprop%t2m(i)
          Coupling%q2mi_cpl    (i) = Sfcprop%q2m(i)
          Coupling%tsfci_cpl   (i) = Sfcprop%tsfc(i)
!         Coupling%tsfci_cpl   (i) = tsfc3(i,3)
          Coupling%psurfi_cpl  (i) = Statein%pgr(i)
        enddo

!  ---  estimate mean albedo for ocean point without ice cover and apply
!       them to net SW heat fluxes

        do i=1,im
!          if (Sfcprop%landfrac(i) < one) then ! Not 100% land
           if (wet(i)) then                    ! some open water
!  ---  compute open water albedo
            xcosz_loc = max( zero, min( one, xcosz(i) ))
            ocalnirdf_cpl(i) = 0.06_kind_phys
            ocalnirbm_cpl(i) = max(albdf, 0.026_kind_phys/(xcosz_loc**1.7_kind_phys+0.065_kind_phys)  &
     &                       + 0.15_kind_phys * (xcosz_loc-0.1_kind_phys) * (xcosz_loc-0.5_kind_phys) &
     &                       * (xcosz_loc-one))
            ocalvisdf_cpl(i) = 0.06
            ocalvisbm_cpl(i) = ocalnirbm_cpl(i)

            Coupling%nnirbmi_cpl(i) = adjnirbmd(i) * (one-ocalnirbm_cpl(i))
            Coupling%nnirdfi_cpl(i) = adjnirdfd(i) * (one-ocalnirdf_cpl(i))
            Coupling%nvisbmi_cpl(i) = adjvisbmd(i) * (one-ocalvisbm_cpl(i))
            Coupling%nvisdfi_cpl(i) = adjvisdfd(i) * (one-ocalvisdf_cpl(i))
          else
            Coupling%nnirbmi_cpl(i) = adjnirbmd(i) - adjnirbmu(i)
            Coupling%nnirdfi_cpl(i) = adjnirdfd(i) - adjnirdfu(i)
            Coupling%nvisbmi_cpl(i) = adjvisbmd(i) - adjvisbmu(i)
            Coupling%nvisdfi_cpl(i) = adjvisdfd(i) - adjvisdfu(i)
          endif
          Coupling%nswsfci_cpl(i) = Coupling%nnirbmi_cpl(i) + Coupling%nnirdfi_cpl(i)   &
                                  + Coupling%nvisbmi_cpl(i) + Coupling%nvisdfi_cpl(i)
          Coupling%nswsfc_cpl(i)  = Coupling%nswsfc_cpl(i)  + Coupling%nswsfci_cpl(i)*dtf
          Coupling%nnirbm_cpl(i)  = Coupling%nnirbm_cpl(i)  + Coupling%nnirbmi_cpl(i)*dtf
          Coupling%nnirdf_cpl(i)  = Coupling%nnirdf_cpl(i)  + Coupling%nnirdfi_cpl(i)*dtf
          Coupling%nvisbm_cpl(i)  = Coupling%nvisbm_cpl(i)  + Coupling%nvisbmi_cpl(i)*dtf
          Coupling%nvisdf_cpl(i)  = Coupling%nvisdf_cpl(i)  + Coupling%nvisdfi_cpl(i)*dtf
        enddo
      endif
      if (Model%lssav) then
        do i=1,im
          Diag%gflux(i)   = Diag%gflux(i)  + gflx(i)  * dtf
          Diag%evbsa(i)   = Diag%evbsa(i)  + evbs(i)  * dtf
          Diag%evcwa(i)   = Diag%evcwa(i)  + evcw(i)  * dtf
          Diag%transa(i)  = Diag%transa(i) + trans(i) * dtf
          Diag%sbsnoa(i)  = Diag%sbsnoa(i) + sbsno(i) * dtf
          Diag%snowca(i)  = Diag%snowca(i) + snowc(i) * dtf
          Diag%snohfa(i)  = Diag%snohfa(i) + snohf(i) * dtf
          Diag%ep(i)      = Diag%ep(i)     + ep1d(i)  * dtf
!*## CCPP ##
!## CCPP ##* sfc_diag_post.F90/sfc_diag_post_run
          Diag%tmpmax(i)  = max(Diag%tmpmax(i), Sfcprop%t2m(i))
          Diag%tmpmin(i)  = min(Diag%tmpmin(i), Sfcprop%t2m(i))

          Diag%spfhmax(i) = max(Diag%spfhmax(i), Sfcprop%q2m(i))
          Diag%spfhmin(i) = min(Diag%spfhmin(i), Sfcprop%q2m(i))
        enddo

        do i=1, im
! find max wind speed then decompose
           tem = sqrt(Diag%u10m(i)*Diag%u10m(i) + Diag%v10m(i)*Diag%v10m(i))
           if (tem > Diag%wind10mmax(i)) then
              Diag%wind10mmax(i) = tem
              Diag%u10mmax(i)    = Diag%u10m(i)
              Diag%v10mmax(i)    = Diag%v10m(i)
           endif

! Compute dew point, first using vapor pressure
           tem = max(Statein%pgr(i) * Sfcprop%q2m(i) / ( con_eps - con_epsm1 * Sfcprop%q2m(i)), qmin)
           Diag%dpt2m(i) = 243.5 / ( ( 17.67 / log(tem/611.2) ) - one) + 273.14
        enddo

      endif
!*## CCPP ##

!!!!!!!!!!!!!!!!!Commented by Moorthi on July 18, 2012 !!!!!!!!!!!!!!!!!!!
!     do i=1,im
!  --- ...  compute coefficient of evaporation in evapc
!
!       if (evapc(i) > one) evapc(i) = one
!  --- ...  over snow cover or ice or sea, coef of evap =one
!       if (weasd(i) > zero .or. slmsk(i) /= one) evapc(i) = one
!     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!## CCPP ##* GFS_surface_generic.F90/GFS_surface_generic_post_run
!  --- ...  Boundary Layer and Free atmospheic turbulence parameterization
!
!  in order to achieve heat storage within canopy layer, in the canopy heat
!    storage parameterization the kinematic sensible and latent heat fluxes
!    (hflx & evap) as surface boundary forcings to the pbl scheme are
!    reduced as a function of surface roughness
!
      do i=1,im
        hflxq(i) = hflx(i)
        evapq(i) = evap(i)
        hffac(i) = one
        hefac(i) = one
      enddo
      if (Model%lheatstrg) then
        do i=1,im
          tem = 0.01_kind_phys * Sfcprop%zorl(i)     ! change unit from cm to m
          tem1 = (tem - z0min) / (z0max - z0min)
          hffac(i) = Model%z0fac * min(max(tem1, zero), one)
          tem = sqrt(Diag%u10m(i)*Diag%u10m(i)+Diag%v10m(i)*Diag%v10m(i))
          tem1 = (tem - u10min) / (u10max - u10min)
          tem2 = one - min(max(tem1, zero), one)
          hffac(i) = tem2 * hffac(i)
          hefac(i) = one + Model%e0fac * hffac(i)
          hffac(i) = one + hffac(i)
          hflxq(i) = hflx(i) / hffac(i)
          evapq(i) = evap(i) / hefac(i)
        enddo
      endif
!*## CCPP ##
!
!     if (lprnt) write(0,*)' tsea3=',Sfcprop%tsfc(ipr),' slmsk=',Sfcprop%slmsk(ipr)     &
!    &, ' kdt=',kdt,' evap=',evapq(ipr)
!     if (lprnt)  write(0,*)' dtdtb=',(dtdt(ipr,k),k=1,15)

!     do i=1,im
!       if (islmsk(i) == 0) then
!         oro_land(i) = zero
!       else
!         oro_land(i) = oro(i)
!       endif
!     enddo

!     write(0,*)' before monin clstp=',clstp,' kdt=',kdt,' lat=',lat
!  if (lprnt) write(0,*)'befmonshoc phii=',Statein%phii(ipr,:)
!  if (lprnt) write(0,*)'befmonshoc=',Statein%tgrs(ipr,:)
!  if (lprnt) write(0,*)'befmonshocdtdt=',dtdt(ipr,1:10)
!  if (lprnt) write(0,*)'befmonshoctkh=',Tbd%phy_f3d(ipr,1:10,ntot3d-1)
!  if (lprnt) write(0,*)'befmonshochflx=',hflxq(ipr),' tsea=',Sfcprop%tsfc(ipr),&
!      ' evap=',evapq(ipr)
!  if (lprnt) write(0,*)'befmonshocq=',Statein%qgrs(ipr,:,1)
!  if (lprnt) write(0,*)'befmonice=',Statein%qgrs(ipr,:,ntiw)
!  if (lprnt) write(0,*)'befmonwat=',Statein%qgrs(ipr,:,ntcw)
!  if (lprnt) write(0,*)'befmonshoctke=',Statein%qgrs(ipr,:,ntke)

!     write(0,*)' before monsho hflx=',hflxq,' me=',me
!     write(0,*)' before monsho evap=',evapq,' me=',me

!## CCPP ##* Note: In the CCPP, the vdftra array is prepared in GFS_PBL_generic.F90/GFS_PBL_generic_pre_run
! regardless of the following conditions. Therefore, this block is redundant in the CCPP and is not included.

      if (nvdiff == ntrac .or. Model%do_ysu .or. Model%shinhong) then
!
        ntiwx = 0

        if (Model%do_shoc) then
          call moninshoc(ix, im, levs, nvdiff, ntcw, nncl, dvdt, dudt, dtdt, dqdt, &
                         Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs,   &
                         Tbd%phy_f3d(1,1,ntot3d-1), prnum, ntke,                   &
                         Statein%prsik(1,1), rb, Sfcprop%zorl, Diag%u10m,          &
                         Diag%v10m, Sfcprop%ffmm, Sfcprop%ffhh, Sfcprop%tsfc,hflxq,&
                         evapq,stress, wind, kpbl, Statein%prsi, del, Statein%prsl,&
                         Statein%prslk, Statein%phii, Statein%phil, dtp, dusfc1,   &
                         dvsfc1, dtsfc1, dqsfc1, dkt, Tbd%hpbl, kinver,           &
                         Model%xkzm_m, Model%xkzm_h, Model%xkzm_s, Model%xkzminv,  &
                         lprnt, ipr, me)
!         if (lprnt) then
!           write(0,*)' aftpbl phii=',Statein%phii(ipr,:)
!           write(0,*)' aftpbl dtdt=',dtdt(ipr,:)
!          write(0,*)' aftpbl dqdtv=',dqdt(ipr,:,1)
!           write(0,*)'aftmonshoc=',Statein%tgrs(ipr,:)
!           write(0,*)'aftmonshocq=',Statein%qgrs(ipr,:,1)
!           write(0,*)'aftmonshoctke=',Statein%qgrs(ipr,:,ntke)
!           write(0,*)'aftmonice=',Statein%qgrs(ipr,:,ntiw)
!           write(0,*)'aftmonwat=',Statein%qgrs(ipr,:,ntcw)
!           write(0,*)'aftmonshocdtdt=',dtdt(ipr,1:10)
!         endif

        else
          if (Model%satmedmf) then
             if (Model%isatmedmf == 0) then   ! initial version of satmedmfvdif (Nov 2018)
!## CCPP ##* satmedmfvdif.F/satmedmfvdif_run Note: The conditional above is checked in satmedmfvdif_init
                call satmedmfvdif(ix, im, levs, nvdiff, ntcw, ntiw, ntke,           &
                       dvdt, dudt, dtdt, dqdt,                                      &
                       Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs,      &
                       Radtend%htrsw, Radtend%htrlw, xmu, garea,                    &
                       Statein%prsik(1,1), rb, Sfcprop%zorl, Diag%u10m, Diag%v10m,  &
                       Sfcprop%ffmm, Sfcprop%ffhh, Sfcprop%tsfc, hflxq, evapq,      &
                       stress, wind, kpbl, Statein%prsi, del, Statein%prsl,         &
                       Statein%prslk, Statein%phii, Statein%phil, dtp,              &
                       Model%dspheat, dusfc1, dvsfc1, dtsfc1, dqsfc1, Tbd%hpbl,    &
                       kinver, Model%xkzm_m, Model%xkzm_h, Model%xkzm_s)
!*## CCPP ##
             elseif (Model%isatmedmf == 1) then   ! updated version of satmedmfvdif (May 2019)
!## CCPP ##* satmedmfvdifq.F/satmedmfvdifq_run Note: The conditional above is checked in satmedmfvdifq_init
                call satmedmfvdifq(ix, im, levs, nvdiff, ntcw, ntiw, ntke,          &
                       dvdt, dudt, dtdt, dqdt,                                      &
                       Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs,      &
                       Radtend%htrsw, Radtend%htrlw, xmu, garea, islmsk, snowd3,    &
                       Statein%prsik(1,1), rb, Sfcprop%zorl, Diag%u10m, Diag%v10m,  &
                       Sfcprop%ffmm, Sfcprop%ffhh, Sfcprop%tsfc, hflxq, evapq,      &
                       stress, wind, kpbl, Statein%prsi, del, Statein%prsl,         &
                       Statein%prslk, Statein%phii, Statein%phil, dtp,              &
                       Model%dspheat, dusfc1, dvsfc1, dtsfc1, dqsfc1, Tbd%hpbl,    &
                       kinver, Model%xkzm_m, Model%xkzm_h, Model%xkzm_s,            &
                       Model%dspfac, Model%bl_upfr, Model%bl_dnfr)
!*## CCPP ##
             endif
          elseif (Model%hybedmf) then
            if (Model%moninq_fac > 0) then
              call moninedmf(ix, im, levs, nvdiff, ntcw, dvdt, dudt, dtdt, dqdt,    &
                           Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs,  &
                           Radtend%htrsw, Radtend%htrlw, xmu, Statein%prsik(1,1),   &
                           rb, Sfcprop%zorl, Diag%u10m, Diag%v10m, Sfcprop%ffmm,    &
                           Sfcprop%ffhh, Sfcprop%tsfc, qss, hflxq, evapq, stress,   &
                           wind, kpbl, Statein%prsi, del, Statein%prsl,             &
                           Statein%prslk, Statein%phii, Statein%phil, dtp,          &
                           Model%dspheat, dusfc1, dvsfc1, dtsfc1, dqsfc1, Tbd%hpbl,&
                           gamt, gamq, dkt, kinver, Model%xkzm_m, Model%xkzm_h,     &
                           Model%xkzm_s, lprnt, ipr,                                &
                           Model%xkzminv, Model%moninq_fac)
            else
              call moninedmf_hafs(ix, im, levs, nvdiff, ntcw, dvdt, dudt, dtdt, dqdt,&
                           Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs,  &
                           Radtend%htrsw, Radtend%htrlw, xmu, Statein%prsik(1,1),   &
                           rb, Sfcprop%zorl, Diag%u10m, Diag%v10m, Sfcprop%ffmm,    &
                           Sfcprop%ffhh, Sfcprop%tsfc, qss, hflxq, evapq, stress,   &
                           wind, kpbl, Statein%prsi, del, Statein%prsl,             &
                           Statein%prslk, Statein%phii, Statein%phil, dtp,          &
                           Model%dspheat, dusfc1, dvsfc1, dtsfc1, dqsfc1, Tbd%hpbl,&
                           gamt, gamq, dkt, kinver, Model%xkzm_m, Model%xkzm_h,     &
                           Model%xkzm_s, lprnt, ipr,                                &
                           Model%xkzminv, Model%moninq_fac,islmsk)
            endif
!     if (lprnt)  write(0,*)' dtdtm=',(dtdt(ipr,k),k=1,15)
!     if (lprnt)  write(0,*)' dqdtm=',(dqdt(ipr,k,1),k=1,15)
          !elseif (Model%do_ysu) then
          !  if (Model%me==0) then
          !      write(0,*) 'Error, ysuvdif only available through CCPP'
          !      stop
          !  end if
          !elseif (Model%shinhong) then
          !  if (Model%me==0) then
          !      write(0,*) 'Error, shinhongvdif only available through CCPP'
          !      stop
          !  end if
          elseif (.not. Model%old_monin) then
            call moninq(ix, im, levs, nvdiff, ntcw, dvdt, dudt, dtdt, dqdt,         &
                        Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs,     &
                        Radtend%htrsw, Radtend%htrlw, xmu, Statein%prsik(1,1), rb,  &
                        Sfcprop%ffmm, Sfcprop%ffhh, Sfcprop%tsfc, qss, hflxq, evapq,&
                        stress, wind, kpbl, Statein%prsi, del, Statein%prsl,        &
                        Statein%prslk, Statein%phii, Statein%phil, dtp,             &
                        Model%dspheat, dusfc1, dvsfc1, dtsfc1, dqsfc1, Tbd%hpbl,   &
                        gamt, gamq, dkt, kinver, Model%xkzm_m, Model%xkzm_h,        &
                        Model%xkzm_s, lprnt, ipr,                                   &
                        Model%xkzminv, Model%moninq_fac, Model%rbcr)
          else
            if (Model%mstrat) then
              call moninp1(ix, im, levs, nvdiff, dvdt, dudt, dtdt, dqdt,            &
                           Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs,  &
                           Statein%prsik(1,1), rb, Sfcprop%ffmm, Sfcprop%ffhh,      &
                           Sfcprop%tsfc, qss, hflxq, evapq, stress, wind, kpbl,     &
                           Statein%prsi, del, Statein%prsl, Statein%prslk,          &
                           Statein%phii, Statein%phil, dtp, dusfc1, dvsfc1,         &
                           dtsfc1, dqsfc1, Tbd%hpbl, gamt, gamq, dkt, kinver,      &
                           Model%xkzm_m, Model%xkzm_h)
            else
              call moninp(ix, im, levs, nvdiff, dvdt, dudt, dtdt, dqdt,             &
                           Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs,  &
                           Statein%prsik(1,1), rb, Sfcprop%ffmm, Sfcprop%ffhh,      &
                           Sfcprop%tsfc, qss, hflxq, evapq, stress, wind, kpbl,     &
                           Statein%prsi, del, Statein%prsl, Statein%phii,           &
                           Statein%phil, dtp, dusfc1, dvsfc1, dtsfc1, dqsfc1,       &
                           Tbd%hpbl, gamt, gamq, dkt, Model%xkzm_m, Model%xkzm_h)
            endif

          endif   ! end if_hybedmf
        endif     ! end if_do_shoc
      else
!*## CCPP ##
!## CCPP ## These variables are allocated in GFS_typedefs.F90/interstitial_create and 
! initialized in GFS_typedefs.F90/interstitial_phys_reset; ntiwx is set in 
! GFS_typedef.F90/interstitial_setup_tracers
        allocate(vdftra(ix,levs,nvdiff), dvdftra(im,levs,nvdiff))
        dvdftra(:,:,:) = zero
        ntiwx = 0
!
!## CCPP ##* GFS_PBL_generic.F90/GFS_PBL_generic_pre_run (ntiwx is set in GFS_typedef.F90/interstitial_setup_tracers)
        if (imp_physics == Model%imp_physics_wsm6) then
! WSM6
          do k=1,levs
            do i=1,im
              vdftra(i,k,1) = Statein%qgrs(i,k,1)
              vdftra(i,k,2) = Statein%qgrs(i,k,ntcw)
              vdftra(i,k,3) = Statein%qgrs(i,k,ntiw)
              vdftra(i,k,4) = Statein%qgrs(i,k,ntoz)
            enddo
          enddo
          kk    = 4
          ntiwx = 3
        elseif (imp_physics == Model%imp_physics_thompson) then
! Thompson
          if(Model%ltaerosol) then
            do k=1,levs
              do i=1,im
                vdftra(i,k,1) = Statein%qgrs(i,k,1)
                vdftra(i,k,2) = Statein%qgrs(i,k,ntcw)
                vdftra(i,k,3) = Statein%qgrs(i,k,ntiw)
                vdftra(i,k,4) = Statein%qgrs(i,k,ntlnc)
                vdftra(i,k,5) = Statein%qgrs(i,k,ntinc)
                vdftra(i,k,6) = Statein%qgrs(i,k,ntoz)
                vdftra(i,k,7) = Statein%qgrs(i,k,ntwa)
                vdftra(i,k,8) = Statein%qgrs(i,k,ntia)
              enddo
            enddo
            kk    = 8
            ntiwx = 3
          else
            do k=1,levs
              do i=1,im
                vdftra(i,k,1) = Statein%qgrs(i,k,1)
                vdftra(i,k,2) = Statein%qgrs(i,k,ntcw)
                vdftra(i,k,3) = Statein%qgrs(i,k,ntiw)
                vdftra(i,k,4) = Statein%qgrs(i,k,ntinc)
                vdftra(i,k,5) = Statein%qgrs(i,k,ntoz)
              enddo
            enddo
            kk    = 5
            ntiwx = 3
          endif
        elseif (imp_physics == Model%imp_physics_mg) then  ! MG3/2
          if (ntgl > 0) then                               ! MG3
            do k=1,levs
              do i=1,im
                vdftra(i,k,1)  = Statein%qgrs(i,k,1)
                vdftra(i,k,2)  = Statein%qgrs(i,k,ntcw)
                vdftra(i,k,3)  = Statein%qgrs(i,k,ntiw)
                vdftra(i,k,4)  = Statein%qgrs(i,k,ntrw)
                vdftra(i,k,5)  = Statein%qgrs(i,k,ntsw)
                vdftra(i,k,6)  = Statein%qgrs(i,k,ntgl)
                vdftra(i,k,7)  = Statein%qgrs(i,k,ntlnc)
                vdftra(i,k,8)  = Statein%qgrs(i,k,ntinc)
                vdftra(i,k,9)  = Statein%qgrs(i,k,ntrnc)
                vdftra(i,k,10) = Statein%qgrs(i,k,ntsnc)
                vdftra(i,k,11) = Statein%qgrs(i,k,ntgnc)
                vdftra(i,k,12) = Statein%qgrs(i,k,ntoz)
              enddo
            enddo
            kk = 12
          else                                             ! MG2
            do k=1,levs
              do i=1,im
                vdftra(i,k,1)  = Statein%qgrs(i,k,1)
                vdftra(i,k,2)  = Statein%qgrs(i,k,ntcw)
                vdftra(i,k,3)  = Statein%qgrs(i,k,ntiw)
                vdftra(i,k,4)  = Statein%qgrs(i,k,ntrw)
                vdftra(i,k,5)  = Statein%qgrs(i,k,ntsw)
                vdftra(i,k,6)  = Statein%qgrs(i,k,ntlnc)
                vdftra(i,k,7)  = Statein%qgrs(i,k,ntinc)
                vdftra(i,k,8)  = Statein%qgrs(i,k,ntrnc)
                vdftra(i,k,9)  = Statein%qgrs(i,k,ntsnc)
                vdftra(i,k,10) = Statein%qgrs(i,k,ntoz)
              enddo
            enddo
            kk = 10
          endif
          ntiwx = 3
!
        elseif (imp_physics == Model%imp_physics_gfdl) then! GFDL MP
          do k=1,levs
            do i=1,im
              vdftra(i,k,1) = Statein%qgrs(i,k,1)
              vdftra(i,k,2) = Statein%qgrs(i,k,ntcw)
              vdftra(i,k,3) = Statein%qgrs(i,k,ntiw)
              vdftra(i,k,4) = Statein%qgrs(i,k,ntrw)
              vdftra(i,k,5) = Statein%qgrs(i,k,ntsw)
              vdftra(i,k,6) = Statein%qgrs(i,k,ntgl)
              vdftra(i,k,7) = Statein%qgrs(i,k,ntoz)
            enddo
          enddo
          kk    = 7
          ntiwx = 3
        elseif (imp_physics == Model%imp_physics_zhao_carr) then  ! Zhao/Carr/Sundqvist
          do k=1,levs
            do i=1,im
              vdftra(i,k,1) = Statein%qgrs(i,k,1)
              vdftra(i,k,2) = Statein%qgrs(i,k,ntcw)
              vdftra(i,k,3) = Statein%qgrs(i,k,ntoz)
            enddo
          enddo
          kk = 3
        endif
!
        if (trans_aero) then
          k1 = kk
          do n=Model%ntchs,Model%ntchm+Model%ntchs-1
            k1 = k1 + 1
            do k=1,levs
              do i=1,im
                vdftra(i,k,k1) = Statein%qgrs(i,k,n)
              enddo
            enddo
          enddo
        endif
!
        if (ntke > 0) then                                 ! prognostic TKE
          ntkev = nvdiff
          do k=1,levs
            do i=1,im
              vdftra(i,k,ntkev) = Statein%qgrs(i,k,ntke)
            enddo
          enddo
        endif
!*## CCPP ##
!       for SHOC nvdiff=ntrac, so the following is not needed unless cplchm is true
!       -----------------------------------------------------
        if (Model%do_shoc) then
!## CCPP ##* moninshoc.f/moninshoc_run Note: The conditional above is not checked in the CCPP scheme;
! therefore the use of this scheme is controlled via the CCPP SDF
            call moninshoc(ix, im, levs, nvdiff, ntcw, nncl, dvdt, dudt, dtdt, dvdftra, &
                           Statein%ugrs, Statein%vgrs, Statein%tgrs, vdftra,            &
                           Tbd%phy_f3d(1,1,ntot3d-1), prnum, ntkev,                     &
                           Statein%prsik(1,1), rb, Sfcprop%zorl, Diag%u10m,             &
                           Diag%v10m, Sfcprop%ffmm, Sfcprop%ffhh, Sfcprop%tsfc, hflxq,  &
                           evapq, stress, wind, kpbl, Statein%prsi, del, Statein%prsl,  &
                           Statein%prslk, Statein%phii, Statein%phil, dtp, dusfc1,      &
                           dvsfc1, dtsfc1, dqsfc1, dkt, Tbd%hpbl, kinver,              &
                           Model%xkzm_m, Model%xkzm_h, Model%xkzm_s, Model%xkzminv,     &
                           lprnt, ipr, me)
!*## CCPP ##
        else
          if (Model%satmedmf) then
             if (Model%isatmedmf == 0) then   ! initial version of satmedmfvdif (Nov 2018)
!## CCPP ##* satmedmfvdif.F/satmedmfvdif_run Note: The conditional above is checked in satmedmfvdif_init
                call satmedmfvdif(ix, im, levs, nvdiff, ntcw, ntiwx, ntkev,           &
                         dvdt, dudt, dtdt, dvdftra,                                   &
                         Statein%ugrs, Statein%vgrs, Statein%tgrs, vdftra,            &
                         Radtend%htrsw, Radtend%htrlw, xmu, garea,                    &
                         Statein%prsik(1,1), rb, Sfcprop%zorl, Diag%u10m, Diag%v10m,  &
                         Sfcprop%ffmm, Sfcprop%ffhh, Sfcprop%tsfc, hflxq, evapq,      &
                         stress, wind, kpbl, Statein%prsi, del, Statein%prsl,         &
                         Statein%prslk, Statein%phii, Statein%phil, dtp,              &
                         Model%dspheat, dusfc1, dvsfc1, dtsfc1, dqsfc1, Tbd%hpbl,    &
                         kinver, Model%xkzm_m, Model%xkzm_h, Model%xkzm_s)
!*## CCPP ##
             elseif (Model%isatmedmf == 1) then   ! updated version of satmedmfvdif (May 2019)
!## CCPP ##* satmedmfvdifq.F/satmedmfvdifq_run Note: The conditional above is checked in satmedmfvdifq_init
                call satmedmfvdifq(ix, im, levs, nvdiff, ntcw, ntiwx, ntkev,          &
                         dvdt, dudt, dtdt, dvdftra,                                   &
                         Statein%ugrs, Statein%vgrs, Statein%tgrs, vdftra,            &
                         Radtend%htrsw, Radtend%htrlw, xmu, garea, islmsk, snowd3,    &
                         Statein%prsik(1,1), rb, Sfcprop%zorl, Diag%u10m, Diag%v10m,  &
                         Sfcprop%ffmm, Sfcprop%ffhh, Sfcprop%tsfc, hflxq, evapq,      &
                         stress, wind, kpbl, Statein%prsi, del, Statein%prsl,         &
                         Statein%prslk, Statein%phii, Statein%phil, dtp,              &
                         Model%dspheat, dusfc1, dvsfc1, dtsfc1, dqsfc1, Tbd%hpbl,    &
                         kinver, Model%xkzm_m, Model%xkzm_h, Model%xkzm_s,            &
                         Model%dspfac, Model%bl_upfr, Model%bl_dnfr)
!*## CCPP ##
             endif
          elseif (Model%hybedmf) then
!## CCPP ## moninedmf.f/hedmf_run Note: The conditional above is not checked in the CCPP scheme;
! therefore the use of this scheme is controlled via the CCPP SDF
           if ( Model%moninq_fac > 0 ) then 
            call moninedmf(ix, im, levs, nvdiff, ntcw, dvdt, dudt, dtdt, dvdftra,       &
                           Statein%ugrs, Statein%vgrs, Statein%tgrs, vdftra,            &
                           Radtend%htrsw, Radtend%htrlw, xmu, Statein%prsik(1,1),       &
                           rb, Sfcprop%zorl, Diag%u10m, Diag%v10m, Sfcprop%ffmm,        &
                           Sfcprop%ffhh, Sfcprop%tsfc, qss, hflxq, evapq, stress,       &
                           wind, kpbl, Statein%prsi, del, Statein%prsl,                 &
                           Statein%prslk, Statein%phii, Statein%phil, dtp,              &
                           Model%dspheat, dusfc1, dvsfc1, dtsfc1, dqsfc1, Tbd%hpbl,    &
                           gamt, gamq, dkt, kinver, Model%xkzm_m, Model%xkzm_h,         &
                           Model%xkzm_s, lprnt, ipr,                                    &
                           Model%xkzminv, Model%moninq_fac)
!*## CCPP ##                                                                                                                                                                                                                               
!## CCPP ##* The following schemes are not in the CCPP yet.
           else
            call moninedmf_hafs(ix, im, levs, nvdiff, ntcw, dvdt, dudt, dtdt, dvdftra,  &
                           Statein%ugrs, Statein%vgrs, Statein%tgrs, vdftra,            &
                           Radtend%htrsw, Radtend%htrlw, xmu, Statein%prsik(1,1),       &
                           rb, Sfcprop%zorl, Diag%u10m, Diag%v10m, Sfcprop%ffmm,        &
                           Sfcprop%ffhh, Sfcprop%tsfc, qss, hflxq, evapq, stress,       &
                           wind, kpbl, Statein%prsi, del, Statein%prsl,                 &
                           Statein%prslk, Statein%phii, Statein%phil, dtp,              &
                           Model%dspheat, dusfc1, dvsfc1, dtsfc1, dqsfc1, Tbd%hpbl,    &
                           gamt, gamq, dkt, kinver, Model%xkzm_m, Model%xkzm_h,         &
                           Model%xkzm_s, lprnt, ipr,                                    &
                           Model%xkzminv, Model%moninq_fac,islmsk)
           endif
          elseif (.not. Model%old_monin) then
            call moninq(ix, im, levs, nvdiff, ntcw, dvdt, dudt, dtdt, dvdftra,          &
                        Statein%ugrs, Statein%vgrs, Statein%tgrs, vdftra,               &
                        Radtend%htrsw, Radtend%htrlw, xmu, Statein%prsik(1,1), rb,      &
                        Sfcprop%ffmm, Sfcprop%ffhh, Sfcprop%tsfc, qss, hflxq, evapq,    &
                        stress, wind, kpbl, Statein%prsi, del, Statein%prsl,            &
                        Statein%prslk, Statein%phii, Statein%phil, dtp,                 &
                        Model%dspheat, dusfc1, dvsfc1, dtsfc1, dqsfc1, Tbd%hpbl,       &
                        gamt, gamq, dkt, kinver, Model%xkzm_m, Model%xkzm_h,            &
                        Model%xkzm_s, lprnt, ipr,                                       &
                        Model%xkzminv, Model%moninq_fac, Model%rbcr)
          else
            if (Model%mstrat) then
              call moninp1(ix, im, levs, nvdiff, dvdt, dudt, dtdt, dvdftra,             &
                           Statein%ugrs, Statein%vgrs, Statein%tgrs, vdftra,            &
                           Statein%prsik(1,1), rb, Sfcprop%ffmm, Sfcprop%ffhh,          &
                           Sfcprop%tsfc, qss, hflxq, evapq, stress, wind, kpbl,         &
                           Statein%prsi, del, Statein%prsl, Statein%prslk,              &
                           Statein%phii, Statein%phil, dtp, dusfc1, dvsfc1,             &
                           dtsfc1, dqsfc1, Tbd%hpbl, gamt, gamq, dkt, kinver,          &
                           Model%xkzm_m, Model%xkzm_h)
            else
              call moninp(ix, im, levs, nvdiff, dvdt, dudt, dtdt, dvdftra,              &
                          Statein%ugrs, Statein%vgrs, Statein%tgrs, vdftra,             &
                          Statein%prsik(1,1), rb, Sfcprop%ffmm, Sfcprop%ffhh,           &
                          Sfcprop%tsfc, qss, hflxq, evapq, stress, wind, kpbl,          &
                          Statein%prsi, del, Statein%prsl, Statein%phii,                &
                          Statein%phil, dtp, dusfc1, dvsfc1, dtsfc1, dqsfc1,            &
                          Tbd%hpbl, gamt, gamq, dkt, Model%xkzm_m, Model%xkzm_h)
            endif

          endif   ! end if_satmedmf
        endif     ! end if_do_shoc
!*## CCPP ##
!## CCPP ## GFS_PBL_generic.F90/GFS_PBL_generic_post_run
        if (ntke > 0) then
          do k=1,levs
            do i=1,im
              dqdt(i,k,ntke)  = dvdftra(i,k,ntkev)
            enddo
          enddo
        endif
        if (trans_aero) then
          k1 = kk
          do n=Model%ntchs,Model%ntchm+Model%ntchs-1
            k1 = k1 + 1
            do k=1,levs
              do i=1,im
                dqdt(i,k,n) = dvdftra(i,k,k1)
              enddo
            enddo
          enddo
        endif
        if (imp_physics == Model%imp_physics_wsm6) then         ! WSM6
          do k=1,levs
            do i=1,im
              dqdt(i,k,1)     = dvdftra(i,k,1)
              dqdt(i,k,ntcw)  = dvdftra(i,k,2)
              dqdt(i,k,ntiw)  = dvdftra(i,k,3)
              dqdt(i,k,ntoz)  = dvdftra(i,k,4)
            enddo
          enddo
        elseif (imp_physics == Model%imp_physics_thompson) then ! Thompson
          if(Model%ltaerosol) then
            do k=1,levs
              do i=1,im
                dqdt(i,k,1)     = dvdftra(i,k,1)
                dqdt(i,k,ntcw)  = dvdftra(i,k,2)
                dqdt(i,k,ntiw)  = dvdftra(i,k,3)
                dqdt(i,k,ntlnc) = dvdftra(i,k,4)
                dqdt(i,k,ntinc) = dvdftra(i,k,5)
                dqdt(i,k,ntoz)  = dvdftra(i,k,6)
                dqdt(i,k,ntwa)  = dvdftra(i,k,7)
                dqdt(i,k,ntia)  = dvdftra(i,k,8)
              enddo
            enddo
          else
            do k=1,levs
              do i=1,im
                dqdt(i,k,1)     = dvdftra(i,k,1)
                dqdt(i,k,ntcw)  = dvdftra(i,k,2)
                dqdt(i,k,ntiw)  = dvdftra(i,k,3)
                dqdt(i,k,ntinc) = dvdftra(i,k,4)
                dqdt(i,k,ntoz)  = dvdftra(i,k,5)
              enddo
            enddo
          endif
        elseif (imp_physics == Model%imp_physics_mg) then    ! MG3/2
          if (ntgl > 0) then                                 ! MG
            do k=1,levs
              do i=1,im
                dqdt(i,k,1)     = dvdftra(i,k,1)
                dqdt(i,k,ntcw)  = dvdftra(i,k,2)
                dqdt(i,k,ntiw)  = dvdftra(i,k,3)
                dqdt(i,k,ntrw)  = dvdftra(i,k,4)
                dqdt(i,k,ntsw)  = dvdftra(i,k,5)
                dqdt(i,k,ntgl)  = dvdftra(i,k,6)
                dqdt(i,k,ntlnc) = dvdftra(i,k,7)
                dqdt(i,k,ntinc) = dvdftra(i,k,8)
                dqdt(i,k,ntrnc) = dvdftra(i,k,9)
                dqdt(i,k,ntsnc) = dvdftra(i,k,10)
                dqdt(i,k,ntgnc) = dvdftra(i,k,11)
                dqdt(i,k,ntoz)  = dvdftra(i,k,12)
              enddo
            enddo
          else                                               ! MG2
            do k=1,levs
              do i=1,im
                dqdt(i,k,1)     = dvdftra(i,k,1)
                dqdt(i,k,ntcw)  = dvdftra(i,k,2)
                dqdt(i,k,ntiw)  = dvdftra(i,k,3)
                dqdt(i,k,ntrw)  = dvdftra(i,k,4)
                dqdt(i,k,ntsw)  = dvdftra(i,k,5)
                dqdt(i,k,ntlnc) = dvdftra(i,k,6)
                dqdt(i,k,ntinc) = dvdftra(i,k,7)
                dqdt(i,k,ntrnc) = dvdftra(i,k,8)
                dqdt(i,k,ntsnc) = dvdftra(i,k,9)
                dqdt(i,k,ntoz)  = dvdftra(i,k,10)
              enddo
            enddo
          endif
!
        elseif (imp_physics == Model%imp_physics_gfdl) then  ! GFDL MP
          do k=1,levs
            do i=1,im
              dqdt(i,k,1)    = dvdftra(i,k,1)
              dqdt(i,k,ntcw) = dvdftra(i,k,2)
              dqdt(i,k,ntiw) = dvdftra(i,k,3)
              dqdt(i,k,ntrw) = dvdftra(i,k,4)
              dqdt(i,k,ntsw) = dvdftra(i,k,5)
              dqdt(i,k,ntgl) = dvdftra(i,k,6)
              dqdt(i,k,ntoz) = dvdftra(i,k,7)
            enddo
          enddo

        elseif (imp_physics == Model%imp_physics_zhao_carr) then   !  Zhao/Carr/Sundqvist
          do k=1,levs
            do i=1,im
              dqdt(i,k,1)    = dvdftra(i,k,1)
              dqdt(i,k,ntcw) = dvdftra(i,k,2)
              dqdt(i,k,ntoz) = dvdftra(i,k,3)
            enddo
          enddo
        endif
!
        deallocate(vdftra, dvdftra)

      endif
!*## CCPP ##

!## CCPP ##* GFS_PBL_generic.F90/GFS_PBL_generic_post_run
      if (Model%cplchm) then
        do i = 1, im
          tem1 = max(Diag%q1(i), qmin)
          tem  = Statein%prsl(i,1) / (con_rd*Diag%t1(i)*(one+con_fvirt*tem1))
          Coupling%ushfsfci(i) = -con_cp * tem * hflx(i) ! upward sensible heat flux
        enddo
        Coupling%dkt     (:,:) = dkt (:,:)
      endif

!     if (lprnt) then
!       write(0,*) ' dusfc1=',dusfc1(ipr),' kdt=',kdt
!       write(0,*) ' dvsfc1=',dvsfc1(ipr),' kdt=',kdt
!       write(0,*)' dtsfc1=',dtsfc1(ipr)*hffac(ipr)
!       write(0,*)' dqsfc1=',dqsfc1(ipr)*hefac(ipr)
!       write(0,*)' dtdtc=',(dtdt(ipr,k),k=1,15)
!       write(0,*)' dqdtc=',(dqdt(ipr,k,1),k=1,15)
!       print *,' dudtm=',dudt(ipr,:)
!     endif

!  --- ...  coupling insertion

      if (Model%cplflx) then
        do i=1,im
          if (Sfcprop%oceanfrac(i) > zero) then               ! Ocean only, NO LAKES
            if (Sfcprop%fice(i) > one - epsln) then ! no open water, thus use results from CICE
              Coupling%dusfci_cpl(i) = Coupling%dusfcin_cpl(i)
              Coupling%dvsfci_cpl(i) = Coupling%dvsfcin_cpl(i)
              Coupling%dtsfci_cpl(i) = Coupling%dtsfcin_cpl(i)
              Coupling%dqsfci_cpl(i) = Coupling%dqsfcin_cpl(i)
            elseif (icy(i) .or. dry(i)) then ! use stress_ocean from sfc_diff for opw component at mixed point
              tem1 = max(Diag%q1(i), qmin)
              rho = Statein%prsl(i,1) / (con_rd*Diag%t1(i)*(one+con_fvirt*tem1))
              if (wind(i) > zero) then
                tem = - rho * stress3(i,3) / wind(i)
                Coupling%dusfci_cpl(i) = tem * Statein%ugrs(i,1)   ! U-momentum flux
                Coupling%dvsfci_cpl(i) = tem * Statein%vgrs(i,1)   ! V-momentum flux
              else
                Coupling%dusfci_cpl(i) = zero
                Coupling%dvsfci_cpl(i) = zero
              endif
              Coupling%dtsfci_cpl(i) = con_cp   * rho * hflx3(i,3) ! sensible heat flux over open ocean
              Coupling%dqsfci_cpl(i) = con_hvap * rho * evap3(i,3) ! latent heat flux over open ocean
            else                                                   ! use results from PBL scheme for 100% open ocean
              Coupling%dusfci_cpl(i) = dusfc1(i)
              Coupling%dvsfci_cpl(i) = dvsfc1(i)
              Coupling%dtsfci_cpl(i) = dtsfc1(i)*hffac(i)
              Coupling%dqsfci_cpl(i) = dqsfc1(i)*hefac(i)
            endif

            Coupling%dusfc_cpl (i) = Coupling%dusfc_cpl(i) + Coupling%dusfci_cpl(i) * dtf
            Coupling%dvsfc_cpl (i) = Coupling%dvsfc_cpl(i) + Coupling%dvsfci_cpl(i) * dtf
            Coupling%dtsfc_cpl (i) = Coupling%dtsfc_cpl(i) + Coupling%dtsfci_cpl(i) * dtf
            Coupling%dqsfc_cpl (i) = Coupling%dqsfc_cpl(i) + Coupling%dqsfci_cpl(i) * dtf
!
          else
            Coupling%dusfc_cpl(i) = huge
            Coupling%dvsfc_cpl(i) = huge
            Coupling%dtsfc_cpl(i) = huge
            Coupling%dqsfc_cpl(i) = huge
          endif ! Ocean only, NO LAKES
        enddo
      endif
!*## CCPP ##
!-------------------------------------------------------lssav if loop ----------
!## CCPP ## GFS_PBL_generic.F90/GFS_PBL_generic_post_run
      if (Model%lssav) then
        do i=1,im
          Diag%dusfc (i) = Diag%dusfc(i) + dusfc1(i)*dtf
          Diag%dvsfc (i) = Diag%dvsfc(i) + dvsfc1(i)*dtf
          Diag%dtsfc (i) = Diag%dtsfc(i) + dtsfc1(i)*hffac(i)*dtf
          Diag%dqsfc (i) = Diag%dqsfc(i) + dqsfc1(i)*hefac(i)*dtf
          Diag%dusfci(i) = dusfc1(i)
          Diag%dvsfci(i) = dvsfc1(i)
          Diag%dtsfci(i) = dtsfc1(i)*hffac(i)
          Diag%dqsfci(i) = dqsfc1(i)*hefac(i)
        enddo
!       if (lprnt) then
!         write(0,*)' dusfc=',dusfc(ipr),' dusfc1=',dusfc1(ipr),' dtf=',
!    &     dtf,' kdt=',kdt,' lat=',lat
!       endif

        if (Model%ldiag3d) then
          if (Model%lsidea) then
            Diag%dt3dt(1:im,:,3) = Diag%dt3dt(1:im,:,3) + dtdt(1:im,:)*dtf
          else
            do k=1,levs
              do i=1,im
                tem  = dtdt(i,k) - (Radtend%htrlw(i,k)+Radtend%htrsw(i,k)*xmu(i))
                Diag%dt3dt(i,k,3) = Diag%dt3dt(i,k,3) + tem*dtf
              enddo
            enddo
          endif
          do k=1,levs
            do i=1,im
              Diag%du3dt(i,k,1) = Diag%du3dt(i,k,1) + dudt(i,k) * dtf
              Diag%du3dt(i,k,2) = Diag%du3dt(i,k,2) - dudt(i,k) * dtf
              Diag%dv3dt(i,k,1) = Diag%dv3dt(i,k,1) + dvdt(i,k) * dtf
              Diag%dv3dt(i,k,2) = Diag%dv3dt(i,k,2) - dvdt(i,k) * dtf
            enddo
          enddo
        endif

      endif   ! end if_lssav
!*## CCPP ##

!## CCPP ##* This block not yet in CCPP.
!
      if (ldiag_ugwp) then
!
! here for COORDE-2018 clean way to store averaged du3dt_pbl
!
        do k=1,levs
          do i=1,im
            Diag%du3dt_pbl(i,k) = Diag%du3dt_pbl(i,k) + dUdt(i,k) * fdaily
            Diag%dv3dt_pbl(i,k) = Diag%dv3dt_pbl(i,k) + dVdt(i,k) * fdaily
            Diag%dt3dt_pbl(i,k) = Diag%dt3dt_pbl(i,k) + dTdt(i,k) * fdaily
!           Tdudt(i,k) = Tdudt(i,k) + dUdt(i,k) * fdaily
!           Tdvdt(i,k) = Tdvdt(i,k) + dVdt(i,k) * fdaily
!           Tdtdt(i,k) = Tdtdt(i,k) + dTdt(i,k) * fdaily
          enddo
        enddo
      endif

      if (Model%lssav) then
        if (Model%ldiag3d) then
          do k=1,levs
            do i=1,im
              Diag%dt3dt(i,k,7) = Diag%dt3dt(i,k,7) - dtdt(i,k)*dtf
            enddo
         enddo
        endif
      endif
!*## CCPP ##

!=============================================================  GW-physics start
!
!            Orographic gravity wave drag parameterization
!            ---------------------------------------------

!## CCPP ##* GFS_GWD_generic.F90/GFS_GWD_generic_pre_run
      if (nmtvr == 14) then         ! current operational - as of 2014
        do i=1,im
! vay-2018
! copy to the separate container to avoid "use" of Sfcprop as "static" field
! sgh30 for TOFD
!
          sgh30(i)  = abs(Sfcprop%oro(i) - Sfcprop%oro_uf(i))

          oc(i)     = Sfcprop%hprime(i,2)
          oa4(i,1)  = Sfcprop%hprime(i,3)
          oa4(i,2)  = Sfcprop%hprime(i,4)
          oa4(i,3)  = Sfcprop%hprime(i,5)
          oa4(i,4)  = Sfcprop%hprime(i,6)
          clx(i,1)  = Sfcprop%hprime(i,7)
          clx(i,2)  = Sfcprop%hprime(i,8)
          clx(i,3)  = Sfcprop%hprime(i,9)
          clx(i,4)  = Sfcprop%hprime(i,10)
          theta(i)  = Sfcprop%hprime(i,11)
          gamma(i)  = Sfcprop%hprime(i,12)
          sigma(i)  = Sfcprop%hprime(i,13)
          elvmax(i) = Sfcprop%hprime(i,14)
        enddo
      elseif (nmtvr == 10) then
        do i=1,im
          oc(i)     = Sfcprop%hprime(i,2)
          oa4(i,1)  = Sfcprop%hprime(i,3)
          oa4(i,2)  = Sfcprop%hprime(i,4)
          oa4(i,3)  = Sfcprop%hprime(i,5)
          oa4(i,4)  = Sfcprop%hprime(i,6)
          clx(i,1)  = Sfcprop%hprime(i,7)
          clx(i,2)  = Sfcprop%hprime(i,8)
          clx(i,3)  = Sfcprop%hprime(i,9)
          clx(i,4)  = Sfcprop%hprime(i,10)
        enddo
      elseif (nmtvr == 6) then
        do i=1,im
          oc(i)     = Sfcprop%hprime(i,2)
          oa4(i,1)  = Sfcprop%hprime(i,3)
          oa4(i,2)  = Sfcprop%hprime(i,4)
          oa4(i,3)  = Sfcprop%hprime(i,5)
          oa4(i,4)  = Sfcprop%hprime(i,6)
          clx(i,1)  = zero
          clx(i,2)  = zero
          clx(i,3)  = zero
          clx(i,4)  = zero
        enddo
      else
!
!  no-oro effects
!
        sgh30(:) = zero
        oc    = zero ; oa4   = zero ; clx    = zero ; theta = zero
        gamma = zero ; sigma = zero ; elvmax = zero

      endif   ! end if_nmtvr
!*## CCPP ##

!## CCPP ##* cires_ugwp.F90/cires_ugwp_run - only V0 is implemented
!
!===== UGWP-start: two versions V0 (knob_ugwp_version=0) and V1(knob_ugwp_version=1)
!
!
      if (Model%do_gwd) then
        if (knob_ugwp_version == 1 ) then
          if (kdt < 2  .and.  me == master) then
            print *, ' VAY-attention UGWP-V1 cires_ugwp_driver '
            print *, ' Only Test-mode by developers '
            stop ' cires_ugwp_driver Test-mode Jan 2019 '
          endif

          call cires_ugwp_driver                                         &
              (im, levs, dtp, kdt, me, lprnt,  Model%lonr,               &
               Model%prslrd0, Model%ral_ts,  Model%cdmbgwd,              &
               Grid%xlat, Grid%xlat_d, Grid%sinlat,  Grid%coslat,        &
               Statein%ugrs, Statein%vgrs, Statein%tgrs,                 &
               Statein%qgrs(1:im,1:levs,1), Statein%prsi, Statein%prsl,  &
               Statein%prslk, Statein%phii, Statein%phil,                &
               del, Sfcprop%hprime, kpbl,                                &
               dusfcg, dvsfcg, gw_dudt,  gw_dvdt, gw_dtdt, gw_kdis,      &
!diagnostics
               Diag%gwp_ax, Diag%gwp_axo, Diag%gwp_axc, Diag%gwp_axf,    &
               Diag%gwp_ay, Diag%gwp_ayo, Diag%gwp_ayc, Diag%gwp_ayf,    &
               Diag%gwp_dtdt, Diag%gwp_kdis, Diag%gwp_okw, Diag%gwp_fgf, &
               Diag%gwp_dcheat, Diag%gwp_precip, Diag%gwp_klevs,         &
               Diag%zmtb,   Diag%gwp_scheat, dlength, cldf,              &
!COORDE-2019 diagnostics    without 3d-fluxes:  tauz_ogw, tauz_ngw ....
               Diag%tau_tofd, Diag%tau_mtb, Diag%tau_ogw, Diag%tau_ngw,  &
               Diag%zmtb, Diag%zlwb, Diag%zogw, Diag%du3dt_mtb,          &
               Diag%du3dt_ogw, Diag%du3dt_tms )

!           do k=1,levs
!             do i=1,im
!               Pdtdt(i,k) = gw_dtdt(i,k)
!               Pdudt(i,k) = gw_dudt(i,k)
!               Pdvdt(i,k) = gw_dvdt(i,k)
!             enddo
!           enddo

        else
!
!knob_ugwp_version == o
!
          if (kdt < 2 .and.  me == master) then
             print *, ' VAY-attention UGWP-V0, Jan 2019 '
          endif
!
          allocate (tke(im,levs))
          if (ntke > 0) then
            tke(1:im,:) = Statein%qgrs(1:im,:,ntke) + dqdt(1:im,:,ntke) * dtp
          else
            tke(:,:) = -9999.0_kind_phys
          endif
!
! tendency without PBL-accumulations
!
          call cires_ugwp_driver_v0                                           &
              (me, master, im, levs, nmtvr, dtp, kdt, Model%lonr,             &
               Model%do_ugwp, Model%do_tofd, Model%cdmbgwd,                   &
               Grid%xlat,    Grid%xlat_d,                                     &
               Grid%sinlat,  Grid%coslat,  Grid%area,                         &
               Statein%ugrs, Statein%vgrs, Statein%tgrs,  Statein%qgrs(1,1,1),&
               Statein%prsi, Statein%prsl, Statein%prslk, Statein%phii,       &
               Statein%phil, del, Sfcprop%hprime(:,1), oc, oa4, clx, theta,   &
               gamma,   sigma,   elvmax, sgh30, kpbl,                         &
               dusfcg,  dvsfcg,  gw_dudt,  gw_dvdt, gw_dtdt, gw_kdis,         &
               tau_tms, tau_mtb, tau_ogw,  tau_ngw,                           &
               zm_mtb,  zm_lwb,  zm_ogw,   ax_mtb, ax_ogw, ax_tms,            &
               Diag%zmtnblck, Diag%rain, ntke, tke, lprnt, ipr)


!Diag for COORDE-2019....... for cires_ugwp_driver_v0

          if (ldiag_ugwp) then
            do i=1,im
              Diag%zmtb(i)     = Diag%zmtb(i)     + fdaily * zm_mtb(i)
              Diag%zlwb(i)     = Diag%zlwb(i)     + fdaily * zm_lwb(i)
              Diag%zogw(i)     = Diag%zogw(i)     + fdaily * zm_ogw(i)
              Diag%tau_tofd(i) = Diag%tau_tofd(i) + fdaily * tau_tms(i)
              Diag%tau_mtb(i)  = Diag%tau_mtb(i)  + fdaily * tau_mtb(i)
              Diag%tau_ogw(i)  = Diag%tau_ogw(i)  + fdaily * tau_ogw(i)
              Diag%tau_ngw(i)  = Diag%tau_ngw(i)  + fdaily * tau_ngw(i)
            enddo
            do k=1,levs
              do i=1,im
                Diag%du3dt_mtb(i,k) = Diag%du3dt_mtb(i,k) + fdaily * ax_mtb(i,k)
                Diag%du3dt_tms(i,k) = Diag%du3dt_tms(i,k) + fdaily * ax_tms(i,k)
                Diag%du3dt_ogw(i,k) = Diag%du3dt_ogw(i,k) + fdaily * ax_ogw(i,k)
                Diag%du3dt_ngw(i,k) = Diag%du3dt_ngw(i,k) + fdaily * gw_dudt(i,k)
                Diag%dv3dt_ngw(i,k) = Diag%dv3dt_ngw(i,k) + fdaily * gw_dvdt(i,k)

!               Tdudt(i,k) = Tdudt(i,k) + gw_dudt(i,k)* fdaily
!               Tdvdt(i,k) = Tdvdt(i,k) + gw_dvdt(i,k)* fdaily
!               Tdtdt(i,k) = Tdtdt(i,k) + gw_dvdt(i,k)* fdaily
              enddo
            enddo
          endif
!
          do k=1,levs
            do i=1,im
              dtdt(i,k) = dtdt(i,k) + gw_dtdt(i,k)
              dudt(i,k) = dudt(i,k) + gw_dudt(i,k)
              dvdt(i,k) = dvdt(i,k) + gw_dvdt(i,k)
            enddo
          enddo

        endif ! if (knob_ugwp_version == 1 ) then
      endif    ! if (do_gwd) then

! *DH UGWD not yet in CCPP
!
!===== UGWP-end ===== ===== =====
!

      if (Model%lssav) then
        do i=1,im
          Diag%dugwd(i) = Diag%dugwd(i) + dusfcg(i)*dtf
          Diag%dvgwd(i) = Diag%dvgwd(i) + dvsfcg(i)*dtf
        enddo
        if (Model%ldiag3d) then
          do k=1,levs
            do i=1,im
              Diag%du3dt(i,k,2) = Diag%du3dt(i,k,2) + dudt(i,k) * dtf
              Diag%dv3dt(i,k,2) = Diag%dv3dt(i,k,2) + dvdt(i,k) * dtf
              Diag%dt3dt(i,k,7) = Diag%dt3dt(i,k,7) + dtdt(i,k) * dtf
            enddo
          enddo
        endif
      endif

!
!===============================================
!
!!    if (ldiag_ugwp) then
!!      do k=1,levs
!!        do i=1,im
!!          Tdudt(i,k) =      Tdudt(i,k) + PdUdt(i,k) * fdaily
!!          Tdvdt(i,k) =      Tdvdt(i,k) + PdVdt(i,k) * fdaily
!!          Tdtdt(i,k) =      Tdtdt(i,k) + PdTdt(i,k) * fdaily
!
!!        enddo
!!      enddo
!!    endif

!## CCPP ##* rayleigh_damp.f/rayleigh_damp_run Note: Conditional IS checked
! within the scheme (returns from scheme if condition is not met)
!    Rayleigh damping  near the model top
      if( .not. Model%lsidea .and. Model%ral_ts > zero) then
        call rayleigh_damp(im, ix, im, levs, dvdt, dudt, dtdt,      &
                           Statein%ugrs, Statein%vgrs, dtp, con_cp, &
                           Model%levr, Statein%pgr, Statein%prsl,   &
                           Model%prslrd0, Model%ral_ts)
      endif

!     if (lprnt) then
!       write(0,*)' tgrs1=',(Statein%tgrs(ipr,k),k=1,10)
!       write(0,*)' dtdt=',(dtdt(ipr,k),k=1,10)
!     endif

! Standard accum-Update before "moist physics" by "PBL + GWP + RF" as in GFS/GSM
!

!## CCPP ##* GFS_suite_interstitial.F90/GFS_suite_stateout_update Note: Terms 
! containing gw_* are related to the CIRES UGWP code and are not currently in
! this scheme.
      do k=1,levs
        do i=1,im
          Stateout%gt0(i,k)  = Statein%tgrs(i,k) + dtdt(i,k) * dtp
          Stateout%gu0(i,k)  = Statein%ugrs(i,k) + dudt(i,k) * dtp
          Stateout%gv0(i,k)  = Statein%vgrs(i,k) + dvdt(i,k) * dtp
        enddo
      enddo
      Stateout%gq0(1:im,:,:) = Statein%qgrs(1:im,:,:) + dqdt(1:im,:,:) * dtp
!*## CCPP ##

!## CCPP ##* This is not in the CCPP yet.
!================================================================================
!     above: updates of the state by UGWP oro-GWS and RF-damp
!  Diag%tav_ugwp & Diag%uav_ugwp(i,k)-Updated U-T state before moist/micro !  physics
!================================================================================

      if (ldiag_ugwp) then
        do k=1,levs
          do i=1,im
            Diag%tav_ugwp(i,k) = Diag%tav_ugwp(i,k) + Stateout%gt0(i,k) * fdaily
            Diag%uav_ugwp(i,k) = Diag%uav_ugwp(i,k) + Stateout%gu0(i,k) * fdaily
!           Diag%vav_ogw(i,k)  = Diag%vav_ogw(i,k)  + Stateout%gv0(i,k) * fdaily
          enddo
        enddo
      endif
!*## CCPP ##

!================================================================================
! It is not clear Do we need it, "ideaca_up", having stability check inside UGWP-module
!## CCPP ##* This is not in the CCPP yet.
      if (Model%lsidea) then            ! idea convective adjustment
        call ideaca_up(Statein%prsi,Stateout%gt0,ix,im,levs+1)
      endif
!*## CCPP ##

!  --- ...  ozone physics

      if (ntoz > 0 .and. ntrac >= ntoz) then
        if (oz_coeff > 4) then
!## CCPP ##* ozphys_2015.f/ozphys_2015_run Note: The conditionals above are not
! checked in the scheme. The scheme's use is controlled by its presense in the
! CCPP SDF
          call ozphys_2015 (ix, im, levs, levozp, dtp,               &
                            Stateout%gq0(1,1,ntoz),                  &
                            Stateout%gq0(1,1,ntoz),                  &
                            Stateout%gt0, oz_pres, Statein%prsl,     &
                            Tbd%ozpl, oz_coeff, del, Model%ldiag3d,  &
                            dq3dt_loc(1,1,6), me)
!*## CCPP ##
!          if (Model%ldiag3d) then
!            do k=1,levs
!              do i=1,im
!                Diag%dq3dt(i,k,6) = dq3dt_loc(i,k,6)
!                Diag%dq3dt(i,k,7) = dq3dt_loc(i,k,7)
!                Diag%dq3dt(i,k,8) = dq3dt_loc(i,k,8)
!                Diag%dq3dt(i,k,9) = dq3dt_loc(i,k,9)
!              enddo
!            enddo
!          endif
        else
!## CCPP ##* ozphys.f/ozphys_run
          call ozphys (ix, im, levs, levozp, dtp,                 &
                       Stateout%gq0(1,1,ntoz),                    &
                       Stateout%gq0(1,1,ntoz),                    &
                       Stateout%gt0, oz_pres, Statein%prsl,       &
                       Tbd%ozpl, oz_coeff, del, Model%ldiag3d,    &
                       dq3dt_loc(1,1,6), me)
!*## CCPP ##
!          if (Model%ldiag3d) then
!            do k=1,levs
!              do i=1,im
!                Diag%dq3dt(i,k,6) = dq3dt_loc(i,k,6)
!                Diag%dq3dt(i,k,7) = dq3dt_loc(i,k,7)
!                Diag%dq3dt(i,k,8) = dq3dt_loc(i,k,8)
!                Diag%dq3dt(i,k,9) = dq3dt_loc(i,k,9)
!              enddo
!            enddo
!          endif
        endif
      endif

      if (Model%h2o_phys) then
!## CCPP ## h2ophys.f/h2ophys_run Note: The conditional is not checked within 
! the scheme. The scheme's use is controlled via the CCPP SDF.
        call h2ophys (ix, im, levs, levh2o, dtp, Stateout%gq0(1,1,1),  &
                      Stateout%gq0(1,1,1), h2o_pres, Statein%prsl,     &
                      Tbd%h2opl, h2o_coeff, Model%ldiag3d,             &
                      dq3dt_loc(1,1,1), me)
!*## CCPP ##
      endif

!  --- ...  to side-step the ozone physics

!      if (ntrac >= 2) then
!        do k = 1, levs
!          gq0(k,ntoz) = qgrs(k,ntoz)
!        enddo
!      endif

!     if (lprnt) then
!       write(0,*) ' levs=',levs,' jcap=',jcap,' dtp',dtp               &
!    &,  ' slmsk=',slmsk(ilon,ilat),' kdt=',kdt
!       print *,' rann=',rann,' ncld=',ncld,' iq=',iq,' lat=',lat
!       print *,' pgr=',pgr
!       print *,' del=',del(ipr,:)
!       print *,' prsl=',prsl(ipr,:)
!       print *,' prslk=',prslk(ipr,:)
!       print *,' rann=',rann(ipr,1)
!       write(0,*)' gt0=',Stateout%gt0(ipr,:)                    &
!    &,         ' kdt=',kdt,' xlon=',grid%xlon(ipr),' xlat=',grid%xlat(ipr)
!       print *,' dtdt=',dtdt(ipr,:)
!       print *,' gu0=',gu0(ipr,:)
!       print *,' gv0=',gv0(ipr,:)
!       write(0,*) ' gt0=',(Stateout%gt0(ipr,k),k=1,levs),' kdt=',kdt
!       write(0,*)' gq0=',(Stateout%gq0(ipr,k,1),k=1,levs)
!       write(0,*)' gq0i2=',(Stateout%gq0(ipr,k,ntiw),k=1,levs)
!       write(0,*)' gq1=',(Stateout%gq0(ipr,k,ntcw),k=1,levs)
!       print *,' vvel=',vvel
!     endif
!     if (lprnt) write(7000,*)' bef convection gu0=',gu0(ipr,:)
!    &,' lat=',lat,' kdt=',kdt,' me=',me
!     if (lprnt) write(7000,*)' bef convection gv0=',gv0(ipr,:)

!## CCPP ## GFS_DCNV_generic.F90/GFS_DCNV_generic_pre_run
      if (Model%ldiag3d) then
        do k=1,levs
          do i=1,im
            dtdt(i,k) = Stateout%gt0(i,k)
            dudt(i,k) = Stateout%gu0(i,k)
            dvdt(i,k) = Stateout%gv0(i,k)
          enddo
        enddo
      elseif (Model%do_cnvgwd) then
        dtdt(1:im,:) = Stateout%gt0(1:im,:)
      endif   ! end if_ldiag3d/cnvgwd

      if (Model%ldiag3d .or. Model%cplchm) then
        dqdt(1:im,:,1) = Stateout%gq0(1:im,:,1)
      endif   ! end if_ldiag3d/cplchm

      if (Model%cplchm) then
        Coupling%dqdti(1:im,:) = zero
      endif   ! end if_cplchm
!*## CCPP ##

!## CCPP ## Only get_prs_fv3.F90/get_phi_fv3_run is a scheme (GFS_HYDRO is assumed to be undefined)
#ifdef GFS_HYDRO
      call get_phi(im, ix, levs, ntrac, Stateout%gt0, Stateout%gq0,    &
                   Model%thermodyn_id, Model%sfcpress_id,              &
                   Model%gen_coord_hybrid Statein%prsi, Statein%prsik, &
                   Statein%prsl, Statein%prslk, Statein%phii, Statein%phil)
#else
!  if (lprnt) write(0,*)'bef get_phi_fv3 gt0=',Stateout%gt0(ipr,:),' kdt=',kdt
!  if (lprnt) write(0,*)'bef get_phi_fv3 gq0=',Stateout%gq0(ipr,:,1),' kdt=',kdt
!  if (lprnt) write(0,*)'bef get_phi_fv3 phii=',Statein%phii(ipr,:),' kdt=',kdt

!GFDL   Adjust the height hydrostatically in a way consistent with FV3 discretization
      call get_phi_fv3 (ix, levs, ntrac, Stateout%gt0, Stateout%gq0, &
                        del_gz, Statein%phii, Statein%phil)

!  if (lprnt) write(0,*)'aft get_phi_fv3 phii=',Statein%phii(ipr,:),' kdt=',kdt
#endif
!*## CCPP ## 

!## CCPP ## These variables are initialized every physics time step through 
! GFS_typedefs.F90/interstitial_phys_reset
      do k=1,levs
        do i=1,im
          clw(i,k,1) = zero
          clw(i,k,2) = -999.9_kind_phys
        enddo
      enddo

      if(imp_physics == Model%imp_physics_thompson) then
        if(Model%ltaerosol) then
          ice00 (:,:) = zero
          liq0  (:,:) = zero
        else
          ice00 (:,:) = zero
        endif
      endif
!*## CCPP ##

!  --- ...  for convective tracer transport (while using ras, csaw, or samf)
!           (the code here implicitly assumes that ntiw=ntcw+1)

!## CCPP ## Most of this code block is in GFS_typedefs.F90/interstitial_setup_tracers except
! for code that needs to be executed every time step (noted below). For those lines,
! they are in GFS_suite_interstitial.F90/GFS_suite_interstitial_3_run.
      ntk       = 0
      tottracer = 0
      if (Model%cscnv .or. Model%satmedmf .or. Model%trans_trac ) then
        otspt(:,:)   = .true.     ! otspt is used only for cscnv
        otspt(1:3,:) = .false.    ! this is for sp.hum, ice and liquid water
!## CCPP ##* GFS_suite_interstitial.F90/GFS_suite_interstitial_3_run
        tracers = 2
        do n=2,ntrac
          if ( n /= ntcw  .and. n /= ntiw  .and. n /= ntclamt .and. &
               n /= ntrw  .and. n /= ntsw  .and. n /= ntrnc   .and. &
!              n /= ntlnc .and. n /= ntinc .and.                    &
               n /= ntsnc .and. n /= ntgl  .and. n /= ntgnc) then
            tracers = tracers + 1
            do k=1,levs
              do i=1,im
                clw(i,k,tracers) = Stateout%gq0(i,k,n)
              enddo
            enddo
!*## CCPP ##
            if (ntke  == n ) then
              otspt(tracers+1,1) = .false.
              ntk = tracers
            endif
            if (ntlnc == n .or. ntinc == n .or. ntrnc == n .or. ntsnc == n .or. ntgnc == n)    &
!           if (ntlnc == n .or. ntinc == n .or. ntrnc == n .or. ntsnc == n .or.&
!               ntrw  == n .or. ntsw  == n .or. ntgl  == n)                    &
                    otspt(tracers+1,1) = .false.
            if (trans_aero .and. Model%ntchs == n) itc = tracers
          endif
        enddo
        tottracer = tracers - 2
      endif   ! end if_ras or cfscnv or samf
!*## CCPP ##

!    if (kdt == 1 .and. me == 0)                                       &
!        write(0,*)' trans_trac=',Model%trans_trac,' tottracer=',      &
!    &               tottracer,' kdt=',kdt,' ntk=',ntk
!## CCPP ##* These variables are initialized in GFS_typedefs.F90/interstitial_phys_reset
      do i=1,im
        ktop(i) = 1
        kbot(i) = levs
      enddo
!*## CCPP ##

!  --- ...  calling condensation/precipitation processes
!           --------------------------------------------
!## CCPP ## GFS_suite_interstitial.F90/GFS_suite_interstitial_3_run
      if (ntcw > 0) then
!       if (imp_physics == Model%imp_physics_mg .and. .not. Model%do_shoc) then ! compute rhc for GMAO macro physics cloud pdf
        if (imp_physics == Model%imp_physics_mg .and. Model%crtrh(2) < half) then ! compute rhc for GMAO macro physics cloud pdf
          do i=1,im
            tx1(i) = one / Statein%prsi(i,1)
            tx2(i) = one - rhc_max*work1(i) - Model%crtrh(1)*work2(i)
            kk     = min(kinver(i), max(2,kpbl(i)))
            tx3(i) = Statein%prsi(i,kk)*tx1(i)
            tx4(i) = Model%crtrh(2) - Model%crtrh(3)*abs(cos(Grid%xlat(i)))
          enddo
          do k = 1, levs
            do i = 1, im
              tem  = Statein%prsl(i,k) * tx1(i)
              tem1 = min(max((tem-tx3(i))*slope_mg, -20.0_kind_phys), 20.0_kind_phys)
!
!     Using crtrh(2) and crtrh(3) from the namelist instead of 0.3 and 0.2
!     and crtrh(1) represents pbl top critical relative humidity
              tem2 = min(max((tx4(i)-tem)*slope_upmg, -20.0_kind_phys), 20.0_kind_phys)

              if (islmsk(i) > 0) then
                tem1 = one / (one+exp(tem1+tem1))
              else
                tem1 = 2.0_kind_phys / (one+exp(tem1+tem1))
              endif
              tem2 = one / (one+exp(tem2))

              rhc(i,k) = min(rhc_max, max(0.7_kind_phys, one-tx2(i)*tem1*tem2))
            enddo
          enddo
        else
          do k=1,levs
            do i=1,im
              kk = max(10,kpbl(i))
              if (k < kk) then
                tem    = Model%crtrh(1) - (Model%crtrh(1)-Model%crtrh(2))     &
                                        * (one-Statein%prslk(i,k)) / (one-Statein%prslk(i,kk))
              else
                tem    = Model%crtrh(2) - (Model%crtrh(2)-Model%crtrh(3))     &
                                        * (Statein%prslk(i,kk)-Statein%prslk(i,k)) / Statein%prslk(i,kk)
              endif
              if (rhc_max > tem) tem = rhc_max * work1(i) + tem * work2(i)
              rhc(i,k) = max(zero, min(one, tem))
            enddo
          enddo
        endif
      endif      ! ntcw > 0
!*## CCPP ##
!
      if (imp_physics == Model%imp_physics_zhao_carr .or. &
          imp_physics == Model%imp_physics_zhao_carr_pdf) then   ! zhao-carr microphysics
!## CCPP ##* precpd.f/zhaocarr_precpd_run
        do i=1,im
          psautco_l(i) = Model%psautco(1)*work1(i) + Model%psautco(2)*work2(i)
          prautco_l(i) = Model%prautco(1)*work1(i) + Model%prautco(2)*work2(i)
        enddo
!*## CCPP ##
!## CCPP ##* GFS_suite_interstitial.F90/GFS_suite_interstitial_3_run
        do k=1,levs
          do i=1,im
            clw(i,k,1) = Stateout%gq0(i,k,ntcw)
          enddo
        enddo
      elseif (imp_physics == Model%imp_physics_gfdl) then
        clw(1:im,:,1) = Stateout%gq0(1:im,:,ntcw)
      elseif (imp_physics == Model%imp_physics_thompson) then
        do k=1,levs
          do i=1,im
            clw(i,k,1) = Stateout%gq0(i,k,ntiw)                    ! ice
            clw(i,k,2) = Stateout%gq0(i,k,ntcw)                    ! water
          enddo
        enddo
        if(Model%ltaerosol) then
          ice00(:,:) = clw(:,:,1)
          liq0(:,:)  = clw(:,:,2)
        else
          ice00(:,:) = clw(:,:,1)
        endif
      elseif (imp_physics == Model%imp_physics_wsm6 .or. &
              imp_physics == Model%imp_physics_mg) then
        do k=1,levs
          do i=1,im
            clw(i,k,1) = Stateout%gq0(i,k,ntiw)                    ! ice
            clw(i,k,2) = Stateout%gq0(i,k,ntcw)                    ! water
          enddo
        enddo
!*## CCPP ##
!## CCPP ## These lines are not in the CCPP since it appeared that they were 
! not needed. These variables are only ever used if (imp_physics == 99 .or. imp_physics == 98)
! which is handled by the first if statement.
      else
        do i=1,im
          psautco_l(i) = Model%psautco(1)*work1(i) + Model%psautco(2)*work2(i)
          prautco_l(i) = Model%prautco(1)*work1(i) + Model%prautco(2)*work2(i)
        enddo
!*## CCPP ##
!## CCPP ##* GFS_suite_interstitial.F90/GFS_suite_interstitial_3_run
        rhc(:,:) = one
!*## CCPP ##
      endif

!     if (lprnt) write(0,*)' clwice=',clw(ipr,:,1)
!     if (lprnt) write(0,*)' clwwat=',clw(ipr,:,2)
!     if (lprnt) write(0,*)' rhc=',rhc(ipr,:)

!
!        Call SHOC if do_shoc is true and shocaftcnv is false
!
!## CCPP ##* gcm_shoc.F90/shoc_run Note: do_shoc is not checked in the scheme, so 
! using this scheme is controlled via the CCPP SDF.
      if (Model%do_shoc .and. .not. Model%shocaftcnv) then
        if (imp_physics == Model%imp_physics_mg) then
          do k=1,levs
            do i=1,im
!## CCPP ##* These lines are commented out in gcm_shoc.F90/shoc_run since they are 
! previously executed in GFS_suite_interstitial.F90/GFS_suite_interstitial_3_run
!             clw(i,k,1) = Stateout%gq0(i,k,ntiw)                    ! ice
!             clw(i,k,2) = Stateout%gq0(i,k,ntcw)                    ! water
!*## CCPP ##
!## CCPP ##* These lines are commented out in gcm_shoc.F90/shoc_run since it is 
! not necessary to copy global variables to local variables
              ncpl(i,k)  = Stateout%gq0(i,k,ntlnc)
              ncpi(i,k)  = Stateout%gq0(i,k,ntinc)
!*## CCPP ##
            enddo
          enddo
          if (abs(Model%fprcp) == 1 .or. mg3_as_mg2) then
            do k=1,levs
              do i=1,im
                qrn(i,k)  = Stateout%gq0(i,k,ntrw)
                qsnw(i,k) = Stateout%gq0(i,k,ntsw)
              enddo
            enddo
          elseif (Model%fprcp > 1) then
            do k=1,levs
              do i=1,im
                qrn(i,k)  = Stateout%gq0(i,k,ntrw)
                qsnw(i,k) = Stateout%gq0(i,k,ntsw) + Stateout%gq0(i,k,ntgl)
!               clw(i,k,1) = clw(i,k,1) + Stateout%gq0(i,k,ntgl)
              enddo
            enddo
          endif
        elseif (imp_physics == Model%imp_physics_gfdl) then  ! GFDL MP - needs modify for condensation
          do k=1,levs
            do i=1,im
              clw(i,k,1) = Stateout%gq0(i,k,ntiw)                    ! ice
              clw(i,k,2) = Stateout%gq0(i,k,ntcw)                    ! water
              qrn(i,k)   = Stateout%gq0(i,k,ntrw)
              qsnw(i,k)  = Stateout%gq0(i,k,ntsw)
            enddo
          enddo
        elseif (imp_physics == Model%imp_physics_zhao_carr .or. &
                imp_physics == Model%imp_physics_zhao_carr_pdf) then
          do k=1,levs
            do i=1,im
              if (abs(Stateout%gq0(i,k,ntcw)) < epsq) then
                Stateout%gq0(i,k,ntcw) = zero
              endif
              tem = Stateout%gq0(i,k,ntcw)              &
     &            * max(zero, MIN(one, (TCR-Stateout%gt0(i,k))*TCRF))
              clw(i,k,1) = tem                              ! ice
              clw(i,k,2) = Stateout%gq0(i,k,ntcw) - tem              ! water
            enddo
          enddo
        endif

!  if (lprnt) write(0,*)'gt01=',Stateout%gt0(ipr,:)
!  if (lprnt) write(0,*)'gq01=',Stateout%gq0(ipr,:,1)
!  if (lprnt) write(0,*)'phii=',Statein%phii(ipr,:),' kdt=',kdt
!  if (lprnt) write(0,*)'clwi=',clw(ipr,:,1)
!  if (lprnt) write(0,*)'clwl=',clw(ipr,:,2)
!  if (lprnt) write(0,*)'befncpi=',ncpi(ipr,:)
!  if (lprnt) write(0,*)'tkh=',Tbd%phy_f3d(ipr,:,ntot3d-1)
!  if (lprnt) write(0,*) ' befshoc hflx=',hflxq(ipr),' evap=',evapq(ipr),&
!    ' stress=',stress(ipr)
!       dtshoc = 60.0
!       dtshoc = 120.0
!       dtshoc = dtp
!       dtshoc = min(dtp, 300.0)
!       nshocm = max(1, nint(dtp/dtshoc))
!       dtshoc = dtp / nshocm
!       do nshoc=1,nshocm
!      if (lprnt) write(0,*)' before shoc tke=',clw(ipr,1:45,ntk), &
!    &' kdt=',kdt,'xlon=',grid%xlon(ipr),' xlat=',grid%xlat(ipr)

!     phy_f3d(1,1,ntot3d-2) - shoc determined sgs clouds
!     phy_f3d(1,1,ntot3d-1) - shoc determined diffusion coefficients
!     phy_f3d(1,1,ntot3d  ) - shoc determined  w'theta'
!
!     dqdt(1:im,:,1) = Stateout%gq0(1:im,:,1)
!     dqdt(1:im,:,2) = Stateout%gq0(1:im,:,ntiw)
!     dqdt(1:im,:,3) = Stateout%gq0(1:im,:,ntcw)
!GFDL lat has no meaning inside of shoc - changed to "1"
!GFDL     call shoc(ix, im, 1, levs, levs+1, dtp, me, lat,
!         call shoc (ix, im, 1, levs, levs+1, dtp, me, 1, Statein%prsl(1,1),  &
!         call shoc (ix, im, 1, levs, levs+1, dtshoc, me, 1, Statein%prsl(1,1),  &
!         call shoc (ix, im, 1, levs, levs+1, dtp, me, 1, Staotein%prsl(1,1),  &
!     write(0,*)' before shoc hflx=',hflxq, ' me=',me
!     write(0,*)' before shoc evap=',evapq,' me=',me
          call shoc (ix, im, levs, levs+1, dtp, me, 1, Statein%prsl(1,1), del,&
                     Statein%phii(1,1), Statein%phil(1,1), Stateout%gu0(1,1), &
                     Stateout%gv0(1,1), Statein%vvl(1,1), Stateout%gt0(1,1),  &
                     Stateout%gq0(1,1,1), clw(1,1,1), clw(1,1,2), qsnw, qrn,  &
                     rhc, Model%sup, Model%shoc_parm(1), Model%shoc_parm(2),  &
                     Model%shoc_parm(3), Model%shoc_parm(4),                  &
                     Model%shoc_parm(5), Tbd%phy_f3d(1,1,ntot3d-2),           &
                     clw(1,1,ntk), hflxq, evapq, prnum,                       &
                     Tbd%phy_f3d(1,1,ntot3d-1), Tbd%phy_f3d(1,1,ntot3d),      &
                     lprnt, ipr, imp_physics, ncpl, ncpi)


!  if (lprnt) write(0,*)'aftncpi=',ncpi(ipr,:)
!       enddo
!         if (imp_physics == Model%imp_physics_mg .and. Model%fprcp > 1) then
!           do k=1,levs
!             do i=1,im
!               clw(i,k,1) = clw(i,k,1) - Stateout%gq0(i,k,ntgl)
!             enddo
!           enddo
!         endif

!     if (lprnt) write(0,*)'aftshocgt0=',Stateout%gt0(ipr,:)
!     if (lprnt) write(0,*)'aftshocgq0=',Stateout%gq0(ipr,:,1)
!     if (lprnt) write(0,*)' aft shoc tke=',clw(ipr,1:25,ntk), &
!    &' kdt=',kdt,'xlon=',grid%xlon(ipr),' xlat=',grid%xlat(ipr)
!     if (lprnt) write(0,*)' aftshoccld=',Tbd%phy_f3d(ipr,:,ntot3d-2)*100
!     if (lprnt) write(0,*)' aftshocice=',clw(ipr,:,1)
!     if (lprnt) write(0,*)' aftshocwat=',clw(ipr,:,2)
!     write(1000+me,*)' at latitude = ',lat
!     rain1 = zero
!     call moist_bud(im,im,ix,levs,me,kdt,con_g,dtp,del,rain1
!    &,              dqdt(1,1,1), dqdt(1,1,2), dqdt(1,1,3)
!    &,              gq0(1,1,1),clw(1,1,2),clw(1,1,1),'shoc      ')
!     tem = 1000.0
!     call moist_bud(im,im,ix,levs,me,kdt,con_g,tem,del,rain1               &
!    &,              dqdt(1,1,1), dqdt(1,1,2), dqdt(1,1,3)                  &
!    &,              Stateout%gq0(1:ix,1:levs,1),clw(1,1,2),clw(1,1,1)      &
!    &,              '   shoc   ', grid%xlon(1:im), grid%xlat(1:im))

!## CCPP ## this is in CCPP's gcm_shoc (but commented out because not needed)
          if (imp_physics == Model%imp_physics_mg) then
            do k=1,levs
              do i=1,im
                Stateout%gq0(i,k,ntlnc) = ncpl(i,k)
                Stateout%gq0(i,k,ntinc) = ncpi(i,k)
              enddo
            enddo
          endif
!*## CCPP ##
!       do k=1,levs
!         do i=1,im
!           sgs_cld(i,k) = sgs_cld(i,k) + shoc_cld(i,k)
!         enddo
!       enddo
!     if (lprnt) write(0,*)' gt03=',gt0(ipr,1:10)
!     if (lprnt) write(0,*)' tke=',clw(ipr,1:10,ntk)

!      if (lprnt) write(1000+me,*)' after shoc tke=',clw(1,:,ntk),
!    &' kdt=',kdt
!       enddo
!
!      do k=1,levs
!      write(1000+me,*)' maxcld=',maxval(sgs_cld(1:im,k)),
!      write(1000+me,*)' maxtkh=',maxval(phy_f3d(1:im,k,ntot3d-1)),
!    &' k=',k,' kdt=',kdt,' lat=',lat
!      enddo

!     write(0,*)' aft shoc gt0=',gt0(1,:),' lat=',lat
!     write(0,*)' aft shoc gq0=',gq0(1,:,1),' lat=',lat
!     write(0,*)' aft shoc gu0=',gu0(1,:),' lat=',lat
!
!*## CCPP ##
      endif   ! if(do_shoc)

!
!  --- ...  calling convective parameterization
!           -----------------------------------
      if (Model%do_deep) then
 
!## CCPP ## GFS_DCNV_generic.F90/GFS_DCNV_generic_pre_run Note: The conditional
! above is not checked within the scheme, so the execution of the code below 
! is controlled via its presence in the CCPP SDF.

!*## CCPP ##
        if (.not. Model%ras .and. .not. Model%cscnv) then

          if (Model%imfdeepcnv == 1) then             ! no random cloud top
!## CCPP ##* sascnvn.F/sascnvn_run
            call sascnvn (im, ix, levs, Model%jcap, dtp, del,                    &
                          Statein%prsl, Statein%pgr, Statein%phil, clw(:,:,1:2), &
                          Stateout%gq0(:,:,1), Stateout%gt0, Stateout%gu0,       &
                          Stateout%gv0, cld1d, rain1, kbot, ktop, kcnv,          &
                          islmsk, Statein%vvl, ncld, ud_mf, dd_mf,               &
                          dt_mf, cnvw, cnvc,                                     &
                          QLCN, QICN, w_upi,cf_upi, CNV_MFD,                     &
!                         QLCN, QICN, w_upi,cf_upi, CNV_MFD, CNV_PRC3,           &
                          CNV_DQLDT,CLCN,CNV_FICE,CNV_NDROP,CNV_NICE,imp_physics,&
                          Model%clam_deep,   Model%c0s_deep,                     &
                          Model%c1_deep,     Model%betal_deep, Model%betas_deep, &
                          Model%evfact_deep, Model%evfactl_deep,                 &
                          Model%pgcon_deep)
!*## CCPP ## 
          elseif (Model%imfdeepcnv == 2) then
!## CCPP ##* GFS_typedefs.f90/interstitial_setup_tracers
            if(.not. Model%satmedmf .and. .not. Model%trans_trac) then
               nsamftrac = 0
            else
               nsamftrac = tottracer
            endif
!*## CCPP ##
!## CCPP ##* samfdeepcnv.f/samfdeepcnv_run
            call samfdeepcnv(im, ix, levs, dtp, itc, Model%ntchm, ntk, nsamftrac,  &
                             del, Statein%prsl, Statein%pgr, Statein%phil, clw,    &
                             Stateout%gq0(:,:,1), Stateout%gt0,                    &
                             Stateout%gu0, Stateout%gv0, Model%fscav,              &
                             cld1d, rain1, kbot, ktop, kcnv,                       &
                             islmsk, garea,                                        &
                             Statein%vvl, ncld, ud_mf, dd_mf, dt_mf, cnvw, cnvc,   &
                             QLCN, QICN, w_upi,cf_upi, CNV_MFD,                    &
!                            QLCN, QICN, w_upi,cf_upi, CNV_MFD, CNV_PRC3,          &
                             CNV_DQLDT,CLCN,CNV_FICE,CNV_NDROP,CNV_NICE,           &
                             imp_physics,                                          &
                             Model%clam_deep,   Model%c0s_deep,                    &
                             Model%c1_deep,  Model%betal_deep, Model%betas_deep,   &
                             Model%evfact_deep, Model%evfactl_deep,                &
                             Model%pgcon_deep,  Model%asolfac_deep,                &
                             Model%do_ca, Model%ca_closure, Model%ca_entr,         & 
                             Model%ca_trigger, Model%nthresh, Coupling%ca_deep,    &
                             Coupling%condition)
!*## CCPP ##
!           if (lprnt) print *,' rain1=',rain1(ipr)
          !elseif (Model%imfdeepcnv == 3) then
          !  if (Model%me==0) then
          !      write(0,*) 'Error, GF convection scheme only available through CCPP'
          !      stop
          !  end if
          !elseif (Model%imfdeepcnv == 4) then
          !  if (Model%me==0) then
          !      write(0,*) 'Error, New Tiedtke convection scheme only available through CCPP'
          !      stop
          !  end if
          elseif (Model%imfdeepcnv == 0) then         ! random cloud top
!## CCPP ##* This is not in the CCPP yet.
            call sascnv (im, ix, levs, Model%jcap, dtp, del,                     &
                         Statein%prsl, Statein%pgr, Statein%phil, clw(:,:,1:2),  &
                         Stateout%gq0(:,:,1), Stateout%gt0, Stateout%gu0,        &
                         Stateout%gv0, cld1d, rain1, kbot, ktop, kcnv,           &
                         islmsk, Statein%vvl, Tbd%rann, ncld,                    &
                         ud_mf, dd_mf, dt_mf, cnvw, cnvc,                        &
                         QLCN, QICN, w_upi,cf_upi, CNV_MFD,                      &
!                        QLCN, QICN, w_upi,cf_upi, CNV_MFD, CNV_PRC3,            &
                         CNV_DQLDT,CLCN,CNV_FICE,CNV_NDROP,CNV_NICE,imp_physics )
!           if (lprnt) print *,' rain1=',rain1(ipr),' rann=',rann(ipr,1)
!*## CCPP ##
          endif

!## CCPP ##* GFS_DCNV_generic.F90/GFS_DCNV_generic_post_run
          if (Model%npdf3d == 3 .and. Model%num_p3d == 4) then
            do k=1,levs
              do i=1,im
                Tbd%phy_f3d(i,k,num2) = cnvw(i,k)
                Tbd%phy_f3d(i,k,num3) = cnvc(i,k)
                cnvw(i,k)             = zero
                cnvc(i,k)             = zero
              enddo
            enddo
          elseif (Model%npdf3d == 0 .and. Model%ncnvcld3d == 1) then
            do k=1,levs
              do i=1,im
                Tbd%phy_f3d(i,k,num2) = cnvw(i,k)
                cnvw(i,k)             = zero
              enddo
            enddo
          endif
!*## CCPP ##
!
        else        ! ras or cscnv
!## CCPP ## cs_conv.F90/cs_conv_pre_run
          fscav(:) = zero
          if (Model%cscnv) then    ! Chikira-Sugiyama  convection scheme (via CSU)

           fswtr(:) = zero
!     write(0,*)' bef cs_cconv phii=',phii(ipr,:)
!    &,' sizefsc=',size(fscav)
!     write(0,*)' bef cs_cconv otspt=',otspt,' kdt=',kdt,' me=',me
!           do k=1,levs
!             do i=1,im
!               dqdt(i,k,1) = Stateout%gq0(i,k,1)
!               dqdt(i,k,2) = clw(i,k,2)
!               dqdt(i,k,3) = clw(i,k,1)
!             enddo
!           enddo
!*## CCPP ##

!
! JLS NOTE:  The convective mass fluxes (dt_mf, dd_mf and ud_mf) passed in and out of cs_conv have not been multiplied by
!            the timestep (i.e, the are in kg/m2/sec) as they are in all other convective schemes.  EMC is aware of this problem, 
!            and in the future will be fixing this discrepancy.  In the meantime, CCPP will use the same mass flux standard_name
!            and long_name as the other convective schemes, where the units are in kg/m2. (Aug 2018)
!
!           if (lprnt) write(0,*)'befcsgt0=',Stateout%gt0(ipr,:)
!           if (lprnt) write(0,*)'befcstke=',clw(ipr,1:25,ntk)

            allocate (sigmatot(im,levs), sigmafrac(im,levs))

! JLS NOTE:  The variable rain1 output from cs_convr (called prec inside the subroutine) is a precipitation flux (kg/m2/sec),
!            not meters LWE like the other schemes.  It is converted to m after the call to cs_convr.

!## CCPP ## cs_conv.F90/cs_conv_run
            call cs_convr (ix, im, levs, ntrac+1, nn, tottracer+3,          &
                           Model%nctp, otspt(1:ntrac+1,1:2), 1,             &
                           kdt, Stateout%gt0, Stateout%gq0(:,:,1:1), rain1, &
                           clw, Statein%phil, Statein%phii, Statein%prsl,   &
                           Statein%prsi, dtp, dtf, ud_mf, dd_mf, dt_mf,     &
                           Stateout%gu0, Stateout%gv0, fscav, fswtr,        &
                           Tbd%phy_fctd, me, wcbmax, Model%cs_parm(3),      &
                           Model%cs_parm(4), Model%cs_parm(9), sigmatot,    &
!                          Model%cs_parm(4), sigmai, sigmatot, vverti,      &
                           Model%do_aw, Model%do_awdd, Model%flx_form,      &
                           lprnt, ipr, kcnv, QLCN, QICN,                    &
                           w_upi, cf_upi, CNV_MFD,           CNV_DQLDT,     &
!                          w_upi, cf_upi, CNV_MFD, CNV_PRC3, CNV_DQLDT,     &
                           CLCN, CNV_FICE, CNV_NDROP, CNV_NICE, imp_physics)
!*## CCPP ##
!           if (lprnt) write(0,*)'aftcsgt0=',Stateout%gt0(ipr,:)
!           if (lprnt) write(0,*)'aftcstke=',clw(ipr,1:25,ntk)

!     write(1000+me,*)' at latitude = ',lat
!     call moist_bud(im,im,ix,levs,me,kdt,con_g,dtp,del,rain1
!    &,                    dqdt(1,1,1), dqdt(1,1,2), dqdt(1,1,3)
!    &,                    gq0(1,1,1),clw(1,1,2),clw(1,1,1),' cs_conv')
!     tem = 1000.0
!     call moist_bud(im,im,ix,levs,me,kdt,con_g,tem,del,rain1               &
!    &,              dqdt(1,1,1), dqdt(1,1,2), dqdt(1,1,3)                  &
!    &,              Stateout%gq0(1:ix,1:levs,1),clw(1,1,2),clw(1,1,1)      &
!    &,              '   cs_conv', grid%xlon(1:im), grid%xlat(1:im))

!## CCPP ##* Not in the CCPP. TODO: Does this need to be in cs_conv_post_run?
            rain1(:) = rain1(:) * (dtp*con_p001)
!## CCPP ##* cs_conv.F90/cs_conv_post_run
            if (Model%do_aw) then
              do k=1,levs
                kk = min(k+1,levs)  ! assuming no cloud top reaches the model top
                do i=1,im                                               !DD
                  sigmafrac(i,k) = half * (sigmatot(i,k)+sigmatot(i,kk))
                enddo
              enddo
            endif
!*## CCPP ##
!           if (lprnt) then
!             write(0,*)' gt01=',stateout%gt0(ipr,:),' kdt=',kdt
!             write(0,*)' gq01=',stateout%gq0(ipr,:,1),' kdt=',kdt
!             write(0,*)' clw1=',clw(ipr,:,1),' kdt=',kdt
!             write(0,*)' clw2=',clw(ipr,:,1),' kdt=',kdt
!             write(0,*)' aft cs rain1=',rain1(ipr)*86400
!             write(0,*)' aft cs rain1=',rain1(ipr)
!           endif

          else      ! ras version 2

!## CCPP ##* This code not yet in CCPP Note: Likely belongs in rascnv_pre.
            if (Model%ccwf(1) >= zero .or. Model%ccwf(2) >= 0) then
              do i=1,im
                ccwfac(i)  = Model%ccwf(1)*work1(i)    + Model%ccwf(2)*work2(i)
                dlqfac(i)  = Model%dlqf(1)*work1(i)    + Model%dlqf(2)*work2(i)
                psaur_l(i) = Model%psauras(1)*work1(i) + Model%psauras(2)*work2(i)
                praur_l(i) = Model%prauras(1)*work1(i) + Model%prauras(2)*work2(i)
              enddo
            else
              do i=1,im
                ccwfac(i)  = -999.0_kind_phys
                dlqfac(i)  = zero
                psaur_l(i) = Model%psauras(1)*work1(i) + Model%psauras(2)*work2(i)
                praur_l(i) = Model%prauras(1)*work1(i) + Model%prauras(2)*work2(i)
              enddo
            endif

!           if  (lprnt) write(0,*) ' calling ras for kdt=',kdt,' me=',me    &
!    &,                            ' lprnt=',lprnt,' ccwfac=',ccwfac(ipr)

!           do k=1,levs
!             do i=1,im
!               dqdt(i,k,1) = Stateout%gq0(i,k,1)
!               dqdt(i,k,2) = clw(i,k,2)
!               dqdt(i,k,3) = clw(i,k,1)
!             enddo
!           enddo

            revap = .true.
!           if (ncld ==2) revap = .false.
            trcmin(:)     = -999999.0_kind_phys
            if (ntk-2 > 0) trcmin(ntk-2) = 1.0e-4_kind_phys
!*## CCPP ##
!           if (lprnt) write(0,*)' gt04bras=',Stateout%gt0(ipr,:)
!           if (lprnt) write(0,*)' gq04bras=',Stateout%gq0(ipr,:,1)
!           if (lprnt) write(0,*)'befrasclw1=',clw(ipr,:,1)
!           if (lprnt) write(0,*)'befrasclw2=',clw(ipr,:,2)
!           if (lprnt) write(0,*)'befrastke=',clw(ipr,:,ntk)
!           if (lprnt) write(0,*)'trcmin=',trcmin(ntk-2),' ntk=',ntk
!## CCPP ## Not in CCPP yet.
            call rascnv (im, ix, levs, dtp, dtf, Tbd%rann, Stateout%gt0,     &
                         Stateout%gq0, Stateout%gu0, Stateout%gv0, clw,      &
                         tottracer, fscav, Statein%prsi, Statein%prsl,       &
                         Statein%prsik, Statein%prslk, Statein%phil,         &
                         Statein%phii, kpbl, cd, rain1, kbot, ktop, kcnv,    &
                         Tbd%phy_f2d(1,Model%num_p2d), Model%flipv, pa2mb,   &
                         me, garea, ccwfac, Model%nrcm, rhc, ud_mf,          &
                         dd_mf, dt_mf, praur_l, Model%wminras(1),            &
                         psaur_l, Model%wminras(2), dlqfac,                  &
                         lprnt, ipr, kdt, revap, QLCN,                       &
                         QICN, w_upi, cf_upi, CNV_MFD,           CNV_DQLDT,  &
!                        QICN, w_upi, cf_upi, CNV_MFD, CNV_PRC3, CNV_DQLDT,  &
                         CLCN, CNV_FICE, CNV_NDROP, CNV_NICE, imp_physics,   &
!                        trcmin)
                         trcmin, ntk)
!*## CCPP ##

!          if (lprnt) write(0,*)' gt04=',Stateout%gt0(ipr,:)
!          if (lprnt) write(0,*)' gq04=',Stateout%gq0(ipr,:,1)
!          if (lprnt) write(0,*)'aftrasclw1=',clw(ipr,:,1)
!          if (lprnt) write(0,*)'aftrasclw2=',clw(ipr,:,2)
!          if (lprnt) write(0,*)'aftrastke=',clw(ipr,:,ntk)

          endif

!     write(1000+me,*)' at latitude = ',lat
!     tem = 1000.0
!     call moist_bud(im,im,ix,levs,me,kdt,con_g,tem,del,rain1               &
!    &,              dqdt(1,1,1), dqdt(1,1,2), dqdt(1,1,3)                  &
!    &,              Stateout%gq0(1:ix,1:levs,1),clw(1,1,2),clw(1,1,1)      &
!    &,              '  ras_conv', grid%xlon(1:im), grid%xlat(1:im))
!     if(lprnt) write(0,*)' after ras rain1=',rain1(ipr),' me=',me,' kdt=',kdt
!    &,' cnv_prc3sum=',sum(cnv_prc3(ipr,1:levs))
!     if (lprnt) write(0,*)' gt04=',gt0(ipr,1:10)

!## CCPP ##* Not in CCPP yet.
          cld1d = 0
!*## CCPP ##

        endif   ! end if_not_ras

!## CCPP ##* GFS_DCNV_generic.F90/GFS_DCNV_generic_post

!*## CCPP ##
      else      ! no parameterized deep convection
!## CCPP ##* GFS_typedefs.F90/interstitial_phys_reset Note: These are only zeroed out 
! initially, prior to calling physics.
        cld1d = zero
        rain1 = zero
        ud_mf = zero
        dd_mf = zero
        dt_mf = zero
!*## CCPP ##
      endif

!     if (lprnt) then
!       write(0,*)' aftcnvgt0=',stateout%gt0(ipr,:),' kdt=',kdt
!       write(0,*)' aftcnvgq0=',(stateout%gq0(ipr,k,1),k=1,levs)
!       write(0,*)' gq0i2=',(stateout%gq0(ipr,k,ntiw),k=1,levs)
!       write(0,*)' aftcnvgq1=',(stateout%gq0(ipr,k,ntcw),k=1,levs)
!     endif
!
!## CCPP ## GFS_DCNV_generic.F90/GFS_DCNV_generic_post_run
      do i=1,im
        Diag%rainc(i) = frain * rain1(i)
      enddo
!
      if (Model%lssav) then
        do i=1,im
          Diag%cldwrk (i)  = Diag%cldwrk (i)  + cld1d(i) * dtf
        enddo

        if (Model%ldiag3d) then
          do k=1,levs
            do i=1,im
              Diag%dt3dt(i,k,4) = Diag%dt3dt(i,k,4) + (Stateout%gt0(i,k)-dtdt(i,k)) * frain
!             Diag%dq3dt(i,k,2) = Diag%dq3dt(i,k,2) + (Stateout%gq0(i,k,1)-dqdt(i,k,1)) * frain
              Diag%du3dt(i,k,3) = Diag%du3dt(i,k,3) + (Stateout%gu0(i,k)-dudt(i,k)) * frain
              Diag%dv3dt(i,k,3) = Diag%dv3dt(i,k,3) + (Stateout%gv0(i,k)-dvdt(i,k)) * frain

!             Diag%upd_mf(i,k)  = Diag%upd_mf(i,k)  + ud_mf(i,k) * (con_g*frain)
!             Diag%dwn_mf(i,k)  = Diag%dwn_mf(i,k)  + dd_mf(i,k) * (con_g*frain)
!             Diag%det_mf(i,k)  = Diag%det_mf(i,k)  + dt_mf(i,k) * (con_g*frain)
            enddo
          enddo
        endif ! if (ldiag3d)

      endif   ! end if_lssav
!*## CCPP ##
!
!## CCPP ##* This block not yet in CCPP.
      if (ldiag_ugwp) then
        tem = frain/dtp
        do k=1,levs
          do i=1,im
!
! frain = dtf / dtp = 1
!
            PdUdt = (Stateout%gu0(i,k)-dudt(i,k)) * tem
            PdVdt = (Stateout%gv0(i,k)-dVdt(i,k)) * tem
            PdTdt = (Stateout%gt0(i,k)-dTdt(i,k)) * tem

            Diag%du3dt_moist(i,k) = Diag%du3dt_moist(i,k) + PdUdt
            Diag%dv3dt_moist(i,k) = Diag%dv3dt_moist(i,k) + PdVdt
            Diag%dt3dt_moist(i,k) = Diag%dt3dt_moist(i,k) + PdTdt
!
! Attention : frain and increments
!
!           Tdudt(i,k) = Tdudt(i,k) + PdUdt * fdaily
!           Tdvdt(i,k) = Tdvdt(i,k) + PdVdt * fdaily
!           Tdtdt(i,k) = Tdtdt(i,k) + PdTdt * fdaily
          enddo
        enddo
      endif
!     if (Model%do_ugwp) then
!
! Put in the instantaneous "Diag%-arrays" to drive UGWP-convective triggers
!     from previous time step we need:  LH-release + cld_top/bot + precip
!
!     endif
!*## CCPP ##

!     if (lprnt) write(7000,*)' bef cnvgwd gu0=',gu0(ipr,:)
!    &,' lat=',lat,' kdt=',kdt,' me=',me
!     if (lprnt) write(7000,*)' bef cnvgwd gv0=',gv0(ipr,:)
!
!----------------Convective gravity wave drag parameterization starting --------

!     if (lprnt) then
!       write(0,*) ' befgwgt0=',Stateout%gt0(ipr,:)
!       write(0,*) ' befgwgq0=',Stateout%gq0(ipr,:,1)
!       write(0,*) ' do_cnvgwd=',Model%do_cnvgwd
!     endif

! DH* this block is in gwdc_pre
!## CCPP ##* gwdc.f/gwdc_pre Note: The conditional above is not in the scheme, so 
! the execution of the code below is controlled by its presence in the CCPP SDF
!  --- ...  calculate maximum convective heating rate 
      if (Model%do_cnvgwd) then         !        call convective gravity wave drag

        allocate(gwdcu(im,levs), gwdcv(im,levs))
!  --- ...  calculate maximum convective heating rate
!           cuhr = temperature change due to deep convection

        do i=1,im
          cumabs(i) = zero
          work4 (i) = zero
        enddo
        do k=1,levs
          do i=1,im
            if (k >= kbot(i) .and. k <= ktop(i)) then
              cumabs(i) = cumabs(i) + (Stateout%gt0(i,k)-dtdt(i,k)) * del(i,k)
              work4(i)  = work4(i)  + del(i,k)
            endif
          enddo
        enddo
        do i=1,im
          if (work4(i) > zero) cumabs(i) = cumabs(i) / (dtp*work4(i))
        enddo
!*## CCPP ##

! DH* 20180817 - note: the above non-CCPP code modifies work3, which until then was defined
! as the ratio of the exner function between midlayer and interface at lowest model layer:
!    work3(i) = Statein%prsik(i,1) / Statein%prslk(i,1)
! This does not happen for the CCPP code, because gwdc_pre uses an internal array
! work3 (maybe not a good name, given that we have work1/2/3 in GFS_physics_driver and
! in the GFS_Interstitial DDT). Therefore, work3 is different from here on until the end
! of GFS_physics_driver. This is ok as long as Model%lgocart is set to .false. - if
! Model%lgocart is set to .true., sfc_diag is called again, which uses work3 as input.
! This work3 used in sfc_diag should be the ratio of the exner function, not the modified
! value derived in the non-CCPP code above. If we get different results for the surface
! diagnstics with Model%lgocart=.true., then the CCPP code is correct! *DH 20180817

!       do i = 1, im
!         do k = kbot(i), ktop(i)
!           do k1 = kbot(i), k
!             cumchr(i,k) = cuhr(i,k1) + cumchr(i,k)
!           enddo
!           cumchr(i,k) = cumchr(i,k) / cumabs(i)
!         enddo
!       enddo

!  --- ...  begin check print ******************************************

!       if (lprnt) then
!         if (kbot(ipr) <= ktop(ipr)) then
!           write(*,*) 'kbot <= ktop     for (lat,lon) = ',             &
!    &            xlon(ipr)*rad2dg,xlat(ipr)*rad2dg
!           write(*,*) 'kcnv kbot ktop dlength  ',kcnv(ipr),            &
!    &            kbot(ipr),ktop(ipr),dlength(ipr)
!           write(*,9000) kdt
!9000       format(/,3x,'k',5x,'cuhr(k)',4x,'cumchr(k)',5x,             &
!    &            'at kdt = ',i4,/)

!           do k = ktop(ipr), kbot(ipr),-1
!             write(*,9010) k,(86400.*cuhr(ipr,k)),(100.*cumchr(ipr,k))
!9010         format(2x,i2,2x,f8.2,5x,f6.0)
!           enddo
!         endif

!         if (fhour >= fhourpr) then
!           print *,' before gwdc in gbphys start print'
!           write(*,*) 'fhour ix im levs = ',fhour,ix,im,levs
!           print *,'dtp  dtf  = ',dtp,dtf

!           write(*,9100)
!9100       format(//,14x,'pressure levels',//                          &
!    &             ' ilev',7x,'prsi',8x,'prsl',8x,'delp',/)

!           k = levs + 1
!           write(*,9110) k,(10.*prsi(ipr,k))
!9110       format(i4,2x,f10.3)

!           do k = levs, 1, -1
!             write(*,9120) k,(10.*prsl(ipr,k)),(10.*del(ipr,k))
!             write(*,9110) k,(10.*prsi(ipr,k))
!           enddo
!9120       format(i4,12x,2(2x,f10.3))

!           write(*,9130)
!9130       format(//,10x,'before gwdc in gbphys',//,' ilev',6x,        &
!    &             'ugrs',9x,'gu0',8x,'vgrs',9x,'gv0',8x,               &
!    &             'tgrs',9x,'gt0',8x,'gt0b',8x,'dudt',8x,'dvdt',/)

!           do k = levs, 1, -1
!             write(*,9140) k,ugrs(ipr,k),gu0(ipr,k),                   &
!    &                        vgrs(ipr,k),gv0(ipr,k),                   &
!    &                        tgrs(ipr,k),gt0(ipr,k),dtdt(ipr,k),       &
!    &                        dudt(ipr,k),dvdt(ipr,k)
!           enddo
!9140       format(i4,9(2x,f10.3))

!           print *,' before gwdc in gbphys end print'
!         endif
!       endif   ! end if_lprnt

!  --- ...  end check print ********************************************

!## CCPP ##* gwdc.f/gwdc_run Note: The conditional above is not in the scheme, so 
! the execution of the code below is controlled by its presence in the CCPP SDF

!GFDL replacing lat with "1"
!       call gwdc(im, ix, im, levs, lat, gu0, gv0, gt0, gq0, dtp,        &
        call gwdc (im, ix, im, levs, 1, Statein%ugrs, Statein%vgrs,      &
                   Statein%tgrs, Statein%qgrs(1,1,1), dtp, Statein%prsl, &
                   Statein%prsi, del, cumabs, ktop, kbot, kcnv, cldf,    &
                   con_g, con_cp, con_rd, con_fvirt, con_pi, dlength,    &
                   lprnt, ipr, Model%fhour, gwdcu, gwdcv, dusfcg, dvsfcg)

!       if (lprnt) then
!         if (fhour >= fhourpr) then
!           print *,' after gwdc in gbphys start print'

!           write(*,9131)
!9131       format(//,10x,'after gwdc in gbphys',//,' ilev',6x,         &
!    &             'ugrs',9x,'gu0',8x,'vgrs',9x,'gv0',8x,               &
!    &             'tgrs',9x,'gt0',8x,'gt0b',7x,'gwdcu',7x,'gwdcv',/)

!           do k = levs, 1, -1
!             write(*,9141) k,ugrs(ipr,k),gu0(ipr,k),                   &
!    &                        vgrs(ipr,k),gv0(ipr,k),                   &
!    &                        tgrs(ipr,k),gt0(ipr,k),dtdt(ipr,k),       &
!    &                        gwdcu(ipr,k),gwdcv(ipr,k)
!           enddo
!9141       format(i4,9(2x,f10.3))

!           print *,' after gwdc in gbphys end print'
!         endif
!       endif

!  --- ...  write out cloud top stress and wind tendencies
!## CCPP ## gwdc.f/gwdc_post_run
        if (Model%lssav) then
          do i=1,im
            Diag%dugwd(i) = Diag%dugwd(i) + dusfcg(i)*dtf
            Diag%dvgwd(i) = Diag%dvgwd(i) + dvsfcg(i)*dtf
          enddo

          if (Model%ldiag3d) then
            do k=1,levs
              do i=1,im
                Diag%du3dt(i,k,4) = Diag%du3dt(i,k,4) + gwdcu(i,k)  * dtf
                Diag%dv3dt(i,k,4) = Diag%dv3dt(i,k,4) + gwdcv(i,k)  * dtf
              enddo
            enddo
          endif
        endif   ! end if_lssav

!  --- ...  update the wind components with  gwdc tendencies

        do k=1,levs
          do i=1,im
            eng0               = half*(Stateout%gu0(i,k)*Stateout%gu0(i,k)+Stateout%gv0(i,k)*Stateout%gv0(i,k))
            Stateout%gu0(i,k)  = Stateout%gu0(i,k) + gwdcu(i,k) * dtp
            Stateout%gv0(i,k)  = Stateout%gv0(i,k) + gwdcv(i,k) * dtp
            eng1               = half*(Stateout%gu0(i,k)*Stateout%gu0(i,k)+Stateout%gv0(i,k)*Stateout%gv0(i,k))
            Stateout%gt0(i,k)  = Stateout%gt0(i,k) + (eng0-eng1)/(dtp*con_cp)
          enddo
!         if (lprnt) write(7000,*)' gu0=',gu0(ipr,k),' gwdcu=',
!    &gwdcu(ipr,k), ' gv0=', gv0(ipr,k),' gwdcv=',gwdcv(ipr,k)
!    &,' k=',k
        enddo
!*## CCPP ##

!       if (lprnt) then
!         if (fhour >= fhourpr) then
!           print *,' after tendency gwdc in gbphys start print'

!           write(*,9132)
!9132       format(//,10x,'after tendency gwdc in gbphys',//,' ilev',6x,&
!    &             'ugrs',9x,'gu0',8x,'vgrs',9x,'gv0',8x,               &
!    &             'tgrs',9x,'gt0',8x,'gt0b',7x,'gwdcu',7x,'gwdcv',/)

!           do k = levs, 1, -1
!             write(*,9142) k,ugrs(ipr,k),gu0(ipr,k),vgrs(ipr,k),       &
!    &              gv0(ipr,k),tgrs(ipr,k),gt0(ipr,k),dtdt(ipr,k),      &
!    &              gwdcu(ipr,k),gwdcv(ipr,k)
!           enddo
!9142       format(i4,9(2x,f10.3))

!           print *,' after tendency gwdc in gbphys end print'
!         endif
!       endif

        deallocate(gwdcu, gwdcv)
      endif   ! end if_cnvgwd (convective gravity wave drag)

!     if (lprnt) then
!       write(0,*) ' befgwegt0=',Stateout%gt0(ipr,:)
!       write(0,*) ' befgwegq0=',Stateout%gq0(ipr,:,1)
!     endif

!     if (lprnt) write(7000,*)' aft cnvgwd gu0=',gu0(ipr,:)
!     if (lprnt) write(7000,*)' aft cnvgwd gv0=',gv0(ipr,:)
!    &,' lat=',lat,' kdt=',kdt,' me=',me
!----------------Convective gravity wave drag parameterization over --------

!## CCPP ## GFS_SCNV_generic.F90/GFS_SCNV_generic_pre_run
      if (Model%ldiag3d) then
        do k=1,levs
          do i=1,im
            dtdt(i,k)   = Stateout%gt0(i,k)
          enddo
        enddo
      endif
!*## CCPP ##

      if (.not. Model%do_shoc) then

        if (Model%shal_cnv) then               ! Shallow convection parameterizations
!                                               --------------------------------------
          if (Model%imfshalcnv == 1) then      ! opr option now at 2014
                                               !-----------------------
!## CCPP ##* shalcnv.F/shalcnv_run
            call shalcnv (im, ix, levs, Model%jcap, dtp, del, Statein%prsl, &
                          Statein%pgr, Statein%phil, clw, Stateout%gq0,     &
                          Stateout%gt0, Stateout%gu0, Stateout%gv0, rain1,  &
                          kbot, ktop, kcnv, islmsk, Statein%vvl, ncld,      &
                          Tbd%hpbl, hflxq, evapq, ud_mf, dt_mf, cnvw, cnvc, &
                          Model%clam_shal, Model%c0s_shal, Model%c1_shal,   &
                          Model%pgcon_shal)
!*## CCPP ##

!## CCPP ##* GFS_SCNV_generic.F90/GFS_SCNV_generic_post_run
            do i=1,im
              Diag%rainc(i) = Diag%rainc(i) + frain * rain1(i)
            enddo
! in shalcnv,  'cnvw' and 'cnvc' are not set to zero
            if (Model%shcnvcw .and. Model%num_p3d == 4 .and. Model%npdf3d == 3) then
              do k=1,levs
                do i=1,im
                  Tbd%phy_f3d(i,k,num2) = Tbd%phy_f3d(i,k,num2) + cnvw(i,k)
                  Tbd%phy_f3d(i,k,num3) = Tbd%phy_f3d(i,k,num3) + cnvc(i,k)
                enddo
              enddo
            elseif (Model%npdf3d == 0 .and. Model%ncnvcld3d == 1) then
              do k=1,levs
                do i=1,im
                  Tbd%phy_f3d(i,k,num2) = Tbd%phy_f3d(i,k,num2) + cnvw(i,k)
                enddo
              enddo
            endif
!*## CCPP ##

          elseif (Model%imfshalcnv == 2) then
!## CCPP ##* GFS_typedef.F90/interstitial_setup_tracers
            if(.not. Model%satmedmf .and. .not. Model%trans_trac) then
               nsamftrac = 0
            else
               nsamftrac = tottracer
            endif
!       if (lprnt) then
!         write(0,*) ' befshgt0=',Stateout%gt0(ipr,:)
!         write(0,*) ' befshgq0=',Stateout%gq0(ipr,:,1)
!       endif
!*## CCPP ##
!## CCPP ##* samfshalcnv.f/samfshalcnv_run
            call samfshalcnv (im, ix, levs, dtp, itc, Model%ntchm, ntk, nsamftrac, &
                              del, Statein%prsl, Statein%pgr, Statein%phil, clw,   &
                              Stateout%gq0(:,:,1), Stateout%gt0,                   &
                              Stateout%gu0, Stateout%gv0, Model%fscav,             &
                              rain1, kbot, ktop, kcnv, islmsk, garea,              &
                              Statein%vvl, ncld, Tbd%hpbl, ud_mf,                 &
                              dt_mf, cnvw, cnvc,                                   &
                              Model%clam_shal,  Model%c0s_shal, Model%c1_shal,     &
                              Model%pgcon_shal, Model%asolfac_shal)
!*## CCPP ##
!## CCPP ##* GFS_SCNV_generic.F90/GFS_SCNV_generic_post_run
            do i=1,im
              Diag%rainc(i) = Diag%rainc(i) + frain * rain1(i)
            enddo

! in  mfshalcnv,  'cnvw' and 'cnvc' are set to zero before computation starts:
            if (Model%shcnvcw .and. Model%num_p3d == 4 .and. Model%npdf3d == 3) then
              do k=1,levs
                do i=1,im
                  Tbd%phy_f3d(i,k,num2) = Tbd%phy_f3d(i,k,num2) + cnvw(i,k)
                  Tbd%phy_f3d(i,k,num3) = Tbd%phy_f3d(i,k,num3) + cnvc(i,k)
                enddo
              enddo
            elseif (Model%npdf3d == 0 .and. Model%ncnvcld3d == 1) then
              do k=1,levs
                do i=1,im
                  Tbd%phy_f3d(i,k,num2) = Tbd%phy_f3d(i,k,num2) +  cnvw(i,k)
                enddo
              enddo
            endif
!*## CCPP ##

          !elseif (Model%imfshalcnv == 3) then
          !if (Model%me==0) write(0,*) "CCPP DEBUG: shallow convection of GF is called in gf_driver"

          !elseif (Model%imfshalcnv == 4) then
          !if (Model%me==0) write(0,*) "CCPP DEBUG: shallow convection of New Tiedtke is called in cu_tiedtke"

          elseif (Model%imfshalcnv == 0) then    ! modified Tiedtke Shallow convecton
                                                 !-----------------------------------
!## CCPP ## This block is not in the CCPP yet.
            levshc(:) = 0
            do k=2,levs
              do i=1,im
                dpshc = 0.3_kind_phys * Statein%prsi(i,1)
                if (Statein%prsi(i,1)-Statein%prsi(i,k) <= dpshc) levshc(i) = k
              enddo
            enddo
            levshcm = 1
            do i=1,im
              levshcm = max(levshcm, levshc(i))
            enddo

!           if (lprnt) print *,' levshcm=',levshcm,' gt0sh=',gt0(ipr,:)
!    &,    ' lat=',lat

            if (Model%mstrat) then             !  As in CFSv2
              call shalcv (im, ix, levshcm, dtp, del, Statein%prsi,        &
                           Statein%prsl, Statein%prslk,kcnv, Stateout%gq0, &
                           Stateout%gt0, levshc, Statein%phil, kinver,     &
                           ctei_r, ctei_rml, lprnt, ipr)
            else
              call shalcvt3 (im, ix, levshcm, dtp, del, Statein%prsi, &
                             Statein%prsl, Statein%prslk, kcnv,       &
                             Stateout%gq0, Stateout%gt0)
            endif
!           if (lprnt) print *,' levshcm=',levshcm,' gt0sha=',gt0(ipr,:)

          endif   ! end if_imfshalcnv

!*## CCPP ##
        endif     ! end if_shal_cnv

!## CCPP ## GFS_SCNV_generic.F90/GFS_SCNV_generic_post_run
        if (Model%lssav) then
          if (Model%ldiag3d) then
            do k=1,levs
              do i=1,im
                Diag%dt3dt(i,k,5) = Diag%dt3dt(i,k,5) + (Stateout%gt0(i,k)  -dtdt(i,k))   * frain
!               Diag%dq3dt(i,k,3) = Diag%dq3dt(i,k,3) + (Stateout%gq0(i,k,1)-dqdt(i,k,1)) * frain
              enddo
            enddo
          endif
        endif   ! end if_lssav

        if (Model%cplchm) then
          do k=1,levs
            do i=1,im
              tem  = (Stateout%gq0(i,k,1)-dqdt(i,k,1)) * frain
              Coupling%dqdti(i,k) = Coupling%dqdti(i,k)  + tem
            enddo
          enddo
        endif
!
        do k=1,levs
          do i=1,im
            if (clw(i,k,2) <= -999.0_kind_phys) clw(i,k,2) = zero
          enddo
        enddo
!*## CCPP ##
!       if (lprnt) then
!         write(0,*)' prsl=',prsl(ipr,:)
!         write(0,*) ' del=',del(ipr,:)
!         write(0,*) ' befshgt0=',gt0(ipr,:)
!         write(0,*) ' befshgq0=',gq0(ipr,:,1)
!       endif

!## CCPP ##* gcm_shoc.F90/shoc_run
      elseif (Model%shocaftcnv) then ! if do_shoc is true and shocaftcnv is true call shoc
        if (imp_physics == Model%imp_physics_mg) then
          do k=1,levs
            do i=1,im
              ncpl(i,k)  = Stateout%gq0(i,k,ntlnc)
              ncpi(i,k)  = Stateout%gq0(i,k,ntinc)
            enddo
          enddo

!       else
!         if (clw(1,1,2) < -999.0) then ! if clw is not partitioned to ice and water
!           do k=1,levs
!             do i=1,im
!               tem = gq0(i,k,ntcw)                                     &
!    &              * max(0.0, MIN(one, (TCR-gt0(i,k))*TCRF))
!               clw(i,k,1) = tem                              ! ice
!               clw(i,k,2) = gq0(i,k,ntcw) - tem              ! water
!             enddo
!           enddo
!         endif     ! Anning ncld ==2
          if (abs(Model%fprcp) == 1 .or. mg3_as_mg2) then
            do k=1,levs
              do i=1,im
                qrn(i,k)  = Stateout%gq0(i,k,ntrw)
                qsnw(i,k) = Stateout%gq0(i,k,ntsw)
              enddo
            enddo
          elseif (Model%fprcp > 1) then
            do k=1,levs
              do i=1,im
                qrn(i,k)  = Stateout%gq0(i,k,ntrw)
                qsnw(i,k) = Stateout%gq0(i,k,ntsw) + Stateout%gq0(i,k,ntgl)
              enddo
            enddo
          endif
        endif

!       dtshoc = 60.0
!       dtshoc = min(dtp, 300.0)
!       nshocm = max(1, nint(dtp/dtshoc))
!       dtshoc = dtp / nshocm
!       do nshoc=1,nshocm
!       call shoc(im, 1, levs, levs+1, dtp, me, lat,        &
!!       call shoc(im, 1, levs, levs+1, dtshoc, me, lat, &
!    &                       prsl(1:im,:), phii (1:im,:),  phil(1:im,:),&
!    &          gu0(1:im,:),gv0(1:im,:), vvl(1:im,:), gt0(1:im,:),     &
!    &                                                   gq0(1:im,:,1), &
!    &          clw(1:im,:,1), clw(1:im,:,2), qsnw, qrn, sgs_cld(1:im,:)&
!    &,         gq0(1:im,:,ntke),                                       &
!    &          phy_f3d(1:im,:,ntot3d-1), phy_f3d(1:im,:,ntot3d),       &
!    &          lprnt, ipr,                                             &
!    &          con_cp, con_g, con_hvap, con_hfus, con_hvap+con_hfus,   &
!    &          con_rv, con_rd, con_pi, con_fvirt)

!GFDL  replace lat with "1:
!       call shoc(ix, im, 1, levs, levs+1, dtshoc, me, lat,             &
!       call shoc (ix, im, 1, levs, levs+1, dtp, me, 1, Statein%prsl(1,1),    &
        call shoc (ix, im, levs, levs+1, dtp, me, 1, Statein%prsl(1,1), del,  &
                   Statein%phii(1,1), Statein%phil(1,1), Stateout%gu0(1,1),   &
                   Stateout%gv0(1,1), Statein%vvl(1,1), Stateout%gt0(1,1),    &
                   Stateout%gq0(1,1,1), clw(1,1,1), clw(1,1,2), qsnw, qrn,    &
                   rhc, Model%sup, Model%shoc_parm(1), Model%shoc_parm(2),    &
                   Model%shoc_parm(3), Model%shoc_parm(4),                    &
                   Model%shoc_parm(5), Tbd%phy_f3d(1,1,ntot3d-2),             &
                   clw(1,1,ntk), hflxq, evapq, prnum,                         &
                   Tbd%phy_f3d(1,1,ntot3d-1), Tbd%phy_f3d(1,1,ntot3d),        &
                   lprnt, ipr, imp_physics, ncpl, ncpi)
!       enddo

        if (imp_physics == Model%imp_physics_mg) then
          do k=1,levs
            do i=1,im
              Stateout%gq0(i,k,ntlnc) = ncpl(i,k)
              Stateout%gq0(i,k,ntinc) = ncpi(i,k)
            enddo
          enddo
        endif

!
!      do k=1,levs
!      write(1000+me,*)' maxtkh=',maxval(phy_f3d(1:im,k,ntot3d-1)), &
!     ' k=',k,' kdt=',kdt,' lat=',lat
!      enddo

!     write(0,*)' aft shoc gt0=',gt0(1,:),' lat=',lat
!     write(0,*)' aft shoc gq0=',gq0(1,:,1),' lat=',lat
!     write(0,*)' aft shoc gu0=',gu0(1,:),' lat=',lat
!
!*## CCPP ##
      endif   ! if( .not. do_shoc)
!
!       if (lprnt) then
!         write(0,*)' prsl=',prsl(ipr,:)
!         write(0,*) ' del=',del(ipr,:)
!         write(0,*) ' aftshgt0=',gt0(ipr,:)
!         write(0,*) ' aftshgq0=',gq0(ipr,:,1)
!       endif
!
!## CCPP ##* GFS_suite_interstitial.F90/GFS_suite_interstitial_4_run
!------------------------------------------------------------------------------
!  --- update the tracers due to deep & shallow cumulus convective transport
!           (except for suspended water and ice)
!
      if (tottracer > 0) then
        tracers = 2
        do n=2,ntrac
!         if ( n /= ntcw .and. n /= ntiw .and. n /= ntclamt) then
          if ( n /= ntcw  .and. n /= ntiw  .and. n /= ntclamt .and. &
               n /= ntrw  .and. n /= ntsw  .and. n /= ntrnc   .and. &
!              n /= ntlnc .and. n /= ntinc .and.                    &
               n /= ntsnc .and. n /= ntgl  .and. n /= ntgnc ) then
              tracers = tracers + 1
            do k=1,levs
              do i=1,im
                Stateout%gq0(i,k,n) = clw(i,k,tracers)
              enddo
            enddo
          endif
        enddo
      endif
!-------------------------------------------------------------------------------
!
      if (ntcw > 0) then

!  for microphysics

        if (imp_physics == Model%imp_physics_zhao_carr     .or. &
            imp_physics == Model%imp_physics_zhao_carr_pdf .or. &
            imp_physics == Model%imp_physics_gfdl) then
          Stateout%gq0(1:im,:,ntcw) = clw(1:im,:,1) + clw(1:im,:,2)
        elseif (ntiw > 0) then
          do k=1,levs
            do i=1,im
              Stateout%gq0(i,k,ntiw) = clw(i,k,1)                     ! ice
              Stateout%gq0(i,k,ntcw) = clw(i,k,2)                     ! water
            enddo
          enddo
          if (imp_physics == Model%imp_physics_thompson) then
            if (Model%ltaerosol) then
              do k=1,levs
                do i=1,im
                  Stateout%gq0(i,k,ntlnc) = Stateout%gq0(i,k,ntlnc)  &
                           +  max(zero, (clw(i,k,2)-liq0(i,k))) / liqm
                  Stateout%gq0(i,k,ntinc) = Stateout%gq0(i,k,ntinc)  &
                           +  max(zero, (clw(i,k,1)-ice00(i,k))) / icem
                enddo
              enddo
            else
              do k=1,levs
                do i=1,im
                  Stateout%gq0(i,k,ntinc) = Stateout%gq0(i,k,ntinc)  &
                           +  max(zero, (clw(i,k,1)-ice00(i,k))) / icem
                enddo
              enddo
            endif
          endif
        else
          do k=1,levs
            do i=1,im
              Stateout%gq0(i,k,ntcw) = clw(i,k,1) + clw(i,k,2)
            enddo
          enddo
        endif   ! end if_ntiw
      else
        do k=1,levs
          do i=1,im
            clw(i,k,1) = clw(i,k,1) + clw(i,k,2)
          enddo
        enddo
      endif   ! end if_ntcw
!*## CCPP ##

!       if (lprnt) then
!         write(0,*)' aft shallow physics kdt=',kdt
!         write(0,*)'qt0s=',Stateout%gt0(ipr,:)
!         write(0,*)'qq0s=',Stateout%gq0(ipr,:,1)
!         write(0,*)'qq0ws=',Stateout%gq0(ipr,:,ntcw)
!         write(0,*)'qq0is=',Stateout%gq0(ipr,:,ntiw)
!         write(0,*)'qq0ntic=',Stateout%gq0(ipr,:,ntinc)
!         write(0,*)'qq0os=',Stateout%gq0(ipr,:,ntoz)
!       endif

!  Legacy routine which determines convectve clouds - should be removed at some point
!## CCPP ## cnvc90.f/cnvc90_run
      call cnvc90 (Model%clstp, im, ix, Diag%rainc, kbot, ktop, levs, Statein%prsi,  &
                   Tbd%acv, Tbd%acvb, Tbd%acvt, Cldprop%cv, Cldprop%cvb, Cldprop%cvt)
!*## CCPP ##

!## CCPP ##* This is not in the CCPP yet.
      if (Model%moist_adj) then       ! moist convective adjustment
!                                     ---------------------------
!
!       To call moist convective adjustment
!
!       if (lprnt) then
!         print *,' prsl=',prsl(ipr,:)
!         print *,' del=',del(ipr,:)
!         print *,' gt0b=',gt0(ipr,:)
!         print *,' gq0b=',gq0(ipr,:,1)
!       endif

        call mstcnv (im, ix, levs, dtp, Stateout%gt0, Stateout%gq0, &
                     Statein%prsl,del, Statein%prslk, rain1,        &
                     Stateout%gq0(1,1,ntcw), rhc, lprnt, ipr)

!       if (lprnt) then
!         print *,' rain1=',rain1(ipr),' rainc=',rainc(ipr)
!         print *,' gt0a=',gt0(ipr,:)
!         print *,' gq0a=',gq0(ipr,:,1)
!       endif
        do i=1,im
          Diag%rainc(i) = Diag%rainc(i) + frain * rain1(i)
        enddo

!       if(Model%lssav) then
! update dqdt_v to include moisture tendency due to surface processes
! dqdt_v : instaneous moisture tendency (kg/kg/sec)
!          if (Model%ldiag3d) then
!            do k=1,levs
!              do i=1,im
!                Diag%dt3dt(i,k,8) = Diag%dt3dt(i,k,8) + (Stateout%gt0(i,k)  -dtdt(i,k)  ) * frain
!!                Diag%dq3dt(i,k,2) = Diag%dq3dt(i,k,2) + (Stateout%gq0(i,k,1)-dqdt(i,k,1)) * frain
!              enddo
!            enddo
!          endif
!       endif
      endif               !       moist convective adjustment over
!*## CCPP ##
!## CCPP ## GFS_MP_generic.F90/GFS_MP_generic_pre_run
      if (Model%ldiag3d .or. Model%do_aw) then
        do k=1,levs
          do i=1,im
            dtdt(i,k)   = Stateout%gt0(i,k)
            dqdt(i,k,1) = Stateout%gq0(i,k,1)
          enddo
        enddo
        do n=ntcw,ntcw+nncl-1
          dqdt(1:im,:,n) = Stateout%gq0(1:im,:,n)
        enddo
      endif
!*## CCPP ##
! dqdt_v : instaneous moisture tendency (kg/kg/sec)
!## CCPP ##* GFS_suite_interstitial.F90/GFS_suite_interstitial_4_run 
! Note: ( these lines are relevant for shallow and deep convection)
      if (Model%cplchm) then
        do k=1,levs
          do i=1,im
            Coupling%dqdti(i,k) = Coupling%dqdti(i,k) * (one / dtf)
          enddo
        enddo
      endif
!*## CCPP ##
!
!     grid-scale condensation/precipitations and microphysics parameterization
!     ------------------------------------------------------------------------
!## CCPP ##* This is not in the CCPP yet.

      if (ncld == 0) then                   ! no cloud microphysics

        call lrgscl (ix, im, levs, dtp, Stateout%gt0, Stateout%gq0, &
                     Statein%prsl, del, Statein%prslk, rain1, clw)
!*## CCPP ##
      else                                  ! all microphysics

        if (imp_physics == Model%imp_physics_zhao_carr) then  ! call zhao/carr/sundqvist microphysics
                                                              ! ------------

          allocate(rainp(im,levs))
!         if (lprnt) then
!           write(0,*)' prsl=',prsl(ipr,:)
!           write(0,*) ' del=',del(ipr,:)
!           write(0,*) ' beflsgt0=',gt0(ipr,:),' kdt=',kdt
!           write(0,*) ' beflsgq0=',gq0(ipr,:,1),' kdt=',kdt
!           write(0,*) ' beflsgw0=',gq0(ipr,:,3),' kdt=',kdt
!         endif
                                              ! ------------------
!## CCPP ##* This is not in the CCPP yet.
          if (Model%do_shoc) then
            call precpd_shoc (im, ix, levs, dtp, del, Statein%prsl,            &
                              Stateout%gq0(1,1,1), Stateout%gq0(1,1,ntcw),     &
                              Stateout%gt0, rain1, Diag%sr, rainp, rhc,        &
                              psautco_l, prautco_l, Model%evpco, Model%wminco, &
                              Tbd%phy_f3d(1,1,ntot3d-2), lprnt, ipr)
!*## CCPP ##
          else
!## CCPP ##* gscond.f/zhaocarr_gscond_run
            call gscond (im, ix, levs, dtp, dtf, Statein%prsl, Statein%pgr,    &
                         Stateout%gq0(1,1,1), Stateout%gq0(1,1,ntcw),          &
                         Stateout%gt0, Tbd%phy_f3d(1,1,1), Tbd%phy_f3d(1,1,2), &
                         Tbd%phy_f2d(1,1), Tbd%phy_f3d(1,1,3),                 &
                         Tbd%phy_f3d(1,1,4), Tbd%phy_f2d(1,2), rhc,lprnt, ipr)
!*## CCPP ##
!## CCPP ##* precpd.f/zhaocarr_precpd_run
            call precpd (im, ix, levs, dtp, del, Statein%prsl,                 &
                        Stateout%gq0(1,1,1), Stateout%gq0(1,1,ntcw),           &
                        Stateout%gt0, rain1, Diag%sr, rainp, rhc, psautco_l,   &
                        prautco_l, Model%evpco, Model%wminco, lprnt, ipr)
!*## CCPP ##
          endif
!         if (lprnt) then
!           write(0,*)' prsl=',prsl(ipr,:)
!           write(0,*) ' del=',del(ipr,:)
!           write(0,*) ' aftlsgt0=',gt0(ipr,:),' kdt=',kdt
!           write(0,*) ' aftlsgq0=',gq0(ipr,:,1),' kdt=',kdt
!           write(0,*) ' aftlsgw0=',gq0(ipr,:,3),' kdt=',kdt
!           write(0,*)' aft precpd rain1=',rain1(1:3),' lat=',lat
!           endif

          deallocate(rainp)
        elseif (imp_physics == Model%imp_physics_zhao_carr_pdf) then ! with pdf clouds
!## CCPP ##* These schemes are not in the CCPP yet.
          allocate(rainp(im,levs))
          call gscondp (im, ix, levs, dtp, dtf, Statein%prsl,        &
                        Statein%pgr, Stateout%gq0(1,1,1),            &
                        Stateout%gq0(1,1,ntcw), Stateout%gt0,        &
                        Tbd%phy_f3d(1,1,1), Tbd%phy_f3d(1,1,2),      &
                        Tbd%phy_f2d(1,1),   Tbd%phy_f3d(1,1,3),      &
                        Tbd%phy_f3d(1,1,4), Tbd%phy_f2d(1,2), rhc,   &
                        Tbd%phy_f3d(1,1,Model%num_p3d+1), Model%sup, &
                        lprnt, ipr, kdt)

          call precpdp (im, ix, levs,  dtp, del, Statein%prsl,       &
                        Statein%pgr, Stateout%gq0(1,1,1),            &
                        Stateout%gq0(1,1,ntcw), Stateout%gt0,        &
                        rain1, Diag%sr, rainp, rhc,                  &
                        Tbd%phy_f3d(1,1,Model%num_p3d+1), psautco_l, &
                        prautco_l, Model%evpco, Model%wminco, lprnt, ipr)
          deallocate(rainp)
!*## CCPP ##

!     if (lprnt) write(0,*) ' rain1=',rain1(ipr),' rainc=',rainc(ipr),' lat=',lat

        elseif (imp_physics == Model%imp_physics_thompson) then      !  Thompson MP
                                            ! ------------
!## CCPP ##* mp_thompson.F90/mp_thompson_run
          ims = 1 ; ime = ix ; kms = 1 ; kme = levs ; its = 1 ; ite = ix ; kts = 1 ; kte = levs

          if (Model%ltaerosol) then
            print*,'aerosol version of the Thompson scheme is not included'

!           call mp_gt_driver(ims,ime,kms,kme,its,ite,kts,kte,                             &
!              Stateout%gq0(1:im,1:levs,1),                                                &
!              Stateout%gq0(1:im,1:levs,Model%ntcw), Stateout%gq0(1:im,1:levs,Model%ntrw), &
!              Stateout%gq0(1:im,1:levs,Model%ntiw), Stateout%gq0(1:im,1:levs,Model%ntsw), &
!              Stateout%gq0(1:im,1:levs,Model%ntgl), Stateout%gq0(1:im,1:levs,Model%ntinc),&
!              Stateout%gq0(1:im,1:im,Model%ntrnc),                                        &
!              Stateout%gt0, Statein%prsl, Statein%vvl, del, dtp, kdt,                     &
!              rain1,                                                                      &
!              Diag%sr,                                                                    &
!!             Diag%refl_10cm, Model%lradar,                                               &
!!             Tbd%phy_f3d(:,:,1),Tbd%phy_f3d(:,:,2),Tbd%phy_f3d(:,:,3),                   & !has_reqc, has_reqi, has_reqs,
!!             ims,ime,kms,kme,its,ite,kts,kte)
!              Tbd%phy_f3d(:,:,1),Tbd%phy_f3d(:,:,2),Tbd%phy_f3d(:,:,3),me,                &
!              nc=Stateout%gq0(1:im,1:levs,Model%ntlnc),                                   &
!              nwfa=Stateout%gq0(1:im,1:levs,Model%ntwa),                                  &
!              nifa=Stateout%gq0(1:im,1:levs,Model%ntia),                                  &
!!             nwfa2d=Sfcprop%nwfa2d(1:im))
!              nwfa2d=Coupling%nwfa2d(1:im))
          else
            call mp_gt_driver(ims,ime,kms,kme,its,ite,kts,kte,                             &
               Stateout%gq0(1:im,1:levs,1),                                                &
               Stateout%gq0(1:im,1:levs,Model%ntcw), Stateout%gq0(1:im,1:levs,Model%ntrw), &
               Stateout%gq0(1:im,1:levs,Model%ntiw), Stateout%gq0(1:im,1:levs,Model%ntsw), &
               Stateout%gq0(1:im,1:levs,Model%ntgl), Stateout%gq0(1:im,1:levs,Model%ntinc),&
               Stateout%gq0(1:im,1:levs,Model%ntrnc),                                      &
!2014v         Stateout%gt0, Statein%prsl, Statein%vvl, del, dtp, kdt,                     &
               Stateout%gt0, Statein%prsl, del, dtp, kdt,                                  &
               rain1,                                                                      &
               Diag%sr,                                                                    &
               islmsk,                                                                     &
               Diag%refl_10cm, Model%lradar,                                               &
               Tbd%phy_f3d(:,:,1),Tbd%phy_f3d(:,:,2),Tbd%phy_f3d(:,:,3),me,Statein%phii)
          endif 
!*## CCPP ##
        elseif (imp_physics == Model%imp_physics_wsm6) then      ! WSM6
!## CCPP ##* This is not in the CCPP yet.                                            ! -----
          ims = 1 ; ime = ix ; kms = 1 ; kme = levs ; its = 1 ; ite = ix ; kts = 1 ; kte = levs

           call wsm6(Stateout%gt0, Statein%phii(1:im,1:levs+1),                                 &
                                Stateout%gq0(1:im,1:levs,1),                                    &
                                Stateout%gq0(1:im,1:levs,Model%ntcw),                           &
                                Stateout%gq0(1:im,1:levs,Model%ntrw),                           &
                                Stateout%gq0(1:im,1:levs,Model%ntiw),                           &
                                Stateout%gq0(1:im,1:levs,Model%ntsw),                           &
                                Stateout%gq0(1:im,1:levs,Model%ntgl),                           &
                                Statein%prsl, del, dtp, rain1,                                  &
                                Diag%sr,                                                        &
                                islmsk,                                                         & 
                                Tbd%phy_f3d(:,:,1),Tbd%phy_f3d(:,:,2),Tbd%phy_f3d(:,:,3),       &
                                ims,ime, kms,kme,                                               &
                                its,ite, kts,kte)
!
!*## CCPP ##
        elseif (imp_physics == Model%imp_physics_mg) then       ! MGB double-moment microphysics
                                                                ! ------------------------------
!## CCPP ##* GFS_typedefs.F90/control_initialize
          kk = 5
          if (Model%fprcp >= 2) kk = 6
!*## CCPP ##
!       Acheng used clw here for other code to run smoothly and minimum change
!       to make the code work. However, the nc and clw should be treated
!       in other procceses too.  August 28/2015; Hope that can be done next
!       year. I believe this will make the physical interaction more reasonable
!       Anning 12/5/2015 changed ntcw hold liquid only
!## CCPP ##* m_micro_insterstitial.F90/m_micro_pre_run
          if (Model%do_shoc) then
            skip_macro = Model%do_shoc
            if (Model%fprcp == 0) then
              do k=1,levs
                do i=1,im
                  clw(i,k,1) = Stateout%gq0(i,k,ntiw)             ! ice
                  clw(i,k,2) = Stateout%gq0(i,k,ntcw)             ! water
                  Tbd%phy_f3d(i,k,1) = Tbd%phy_f3d(i,k,ntot3d-2) ! clouds from shoc
                enddo
              enddo
            elseif (abs(Model%fprcp) == 1 .or. mg3_as_mg2) then
              do k=1,levs
                do i=1,im
                  clw(i,k,1) = Stateout%gq0(i,k,ntiw)             ! ice
                  clw(i,k,2) = Stateout%gq0(i,k,ntcw)             ! water
                  qrn(i,k)   = Stateout%gq0(i,k,ntrw)
                  qsnw(i,k)  = Stateout%gq0(i,k,ntsw)
                  ncpr(i,k)  = Stateout%gq0(i,k,ntrnc)
                  ncps(i,k)  = Stateout%gq0(i,k,ntsnc)
                  Tbd%phy_f3d(i,k,1) = Tbd%phy_f3d(i,k,ntot3d-2) ! clouds from shoc
                enddo
              enddo
            else
              do k=1,levs
                do i=1,im
                  clw(i,k,1) = Stateout%gq0(i,k,ntiw)             ! ice
                  clw(i,k,2) = Stateout%gq0(i,k,ntcw)             ! water
                  qrn(i,k)   = Stateout%gq0(i,k,ntrw)
                  qsnw(i,k)  = Stateout%gq0(i,k,ntsw)
                  qgl(i,k)   = Stateout%gq0(i,k,ntgl)
                  ncpr(i,k)  = Stateout%gq0(i,k,ntrnc)
                  ncps(i,k)  = Stateout%gq0(i,k,ntsnc)
                  ncgl(i,k)  = Stateout%gq0(i,k,ntgnc)
                  Tbd%phy_f3d(i,k,1) = Tbd%phy_f3d(i,k,ntot3d-2) ! clouds from shoc
                enddo
              enddo

            endif

          else
                                                     ! clouds from t-dt and cnvc
            if (Model%fprcp == 0 ) then
              do k=1,levs
                do i=1,im
                  clw(i,k,1) = Stateout%gq0(i,k,ntiw)             ! ice
                  clw(i,k,2) = Stateout%gq0(i,k,ntcw)             ! water
                enddo
              enddo
            elseif (abs(Model%fprcp) == 1 .or. mg3_as_mg2) then
              do k=1,levs
                do i=1,im
                  clw(i,k,1) = Stateout%gq0(i,k,ntiw)             ! ice
                  clw(i,k,2) = Stateout%gq0(i,k,ntcw)             ! water
                  qrn(i,k)   = Stateout%gq0(i,k,ntrw)
                  qsnw(i,k)  = Stateout%gq0(i,k,ntsw)
                  ncpr(i,k)  = Stateout%gq0(i,k,ntrnc)
                  ncps(i,k)  = Stateout%gq0(i,k,ntsnc)
                enddo
              enddo
            else
              do k=1,levs
                do i=1,im
                  clw(i,k,1) = Stateout%gq0(i,k,ntiw)             ! ice
                  clw(i,k,2) = Stateout%gq0(i,k,ntcw)             ! water
                  qrn(i,k)   = Stateout%gq0(i,k,ntrw)
                  qsnw(i,k)  = Stateout%gq0(i,k,ntsw)
                  qgl(i,k)   = Stateout%gq0(i,k,ntgl)
                  ncpr(i,k)  = Stateout%gq0(i,k,ntrnc)
                  ncps(i,k)  = Stateout%gq0(i,k,ntsnc)
                  ncgl(i,k)  = Stateout%gq0(i,k,ntgnc)
                enddo
              enddo
            endif
          endif
! add convective cloud fraction
          do k = 1,levs
            do i = 1,im
              Tbd%phy_f3d(i,k,1) = min(one, Tbd%phy_f3d(i,k,1) + clcn(i,k))
            enddo
          enddo
!*## CCPP ##

!       notice clw ix instead of im
!       call m_micro_driver(im,ix,levs,flipv,del,dtp,prsl,prsi,
!    &    prslk,prsik,pgr,vvl,clw(1,1,2), QLCN, clw(1,1,1),QICN,
!       if (lprnt) write(0,*)' cnv_mfdbef=',cnv_mfd(ipr,:),' flipv=',flipv
!       if(lprnt) write(0,*) ' befgt0=',Stateout%gt0(ipr,:),' kdt=',kdt
!       if(lprnt) write(0,*) ' befgq0=',Stateout%gq0(ipr,:,1),' kdt=',kdt
!       if(lprnt) write(0,*) ' befntlnc=',Stateout%gq0(ipr,:,ntlnc),' kdt=',kdt
!       if(lprnt) write(0,*) ' befntinc=',Stateout%gq0(ipr,:,ntinc),' kdt=',kdt
!       if (lprnt) write(0,*)' clw1bef=',clw(ipr,:,1),' kdt=',kdt
!       if (lprnt) write(0,*)' clw2bef=',clw(ipr,:,2),' kdt=',kdt
!       if (lprnt) write(0,*)' qrnb=',qrn(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qsnwb=',qsnw(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qglb=',qgl(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' rhc=',rhc(ipr,:),' kdt=',kdt,' kk=',kk
!       if (lprnt) write(0,*)' cloudsb=',Tbd%phy_f3d(ipr,:,1)*100,' kdt=',kdt
!       if (lprnt) write(0,*)' cloudsb=',Tbd%phy_f3d(ipr,:,1)*100,' kdt=',kdt
!       if (lprnt) write(0,*)' clcn=',clcn(ipr,:)*100,' kdt=',kdt
!       txa(:,:) = Stateout%gq0(:,:,1)
!       do k=1,levs
!       write(1000+me,*)' maxwatncb=',maxval(Stateout%gq0(1:im,k,ntlnc)),' k=',k,' kdt',kdt
!       enddo

!## CCPP ##* m_micro.F90/m_micro_run
        call m_micro_driver (im, ix, levs, Model%flipv, dtp,  Statein%prsl,      &
                             Statein%prsi, Statein%phil, Statein%phii,           &
                             Statein%vvl, clw(1,1,2), QLCN, clw(1,1,1), QICN,    &
                             Radtend%htrlw, Radtend%htrsw, w_upi, cf_upi,        &
                             FRLAND, Tbd%Hpbl, CNV_MFD,           CNV_DQLDT,    &
!                            FRLAND, Tbd%Hpbl, CNV_MFD, CNV_PRC3, CNV_DQLDT,    &
                             CLCN, Stateout%gu0, Stateout%gv0, Diag%dusfc,       &
                             Diag%dvsfc, dusfc1, dvsfc1, dusfc1, dvsfc1,         &
                             CNV_FICE, CNV_NDROP, CNV_NICE, Stateout%gq0(1,1,1), &
                             Stateout%gq0(1,1,ntcw),                             &
                             Stateout%gq0(1,1,ntiw), Stateout%gt0, rain1,        &
                             Diag%sr, Stateout%gq0(1,1,ntlnc),                   &
                             Stateout%gq0(1,1,ntinc), Model%fprcp, qrn,          &
                             qsnw, qgl, ncpr, ncps, ncgl,                        &
                             Tbd%phy_f3d(1,1,1),  kbot,                          &
                             Tbd%phy_f3d(1,1,2),  Tbd%phy_f3d(1,1,3),            &
                             Tbd%phy_f3d(1,1,4),  Tbd%phy_f3d(1,1,5),            &
                             Tbd%phy_f3d(1,1,kk), Tbd%aer_nm,                    &
                             Tbd%in_nm, Tbd%ccn_nm, Model%iccn,                  &
                             skip_macro,                 lprnt,                  &
!                            skip_macro, cn_prc, cn_snr, lprnt,                  &
!                            ipr, kdt, Grid%xlat, Grid%xlon)
                             Model%mg_alf, Model%mg_qcmin, Model%pdfflag,        &
                             ipr, kdt, Grid%xlat, Grid%xlon, rhc)
!*## CCPP ##
!     do k=1,levs
!     write(1000+me,*)' maxwatnca=',maxval(Stateout%gq0(1:im,k,ntlnc)),' k=',k,' kdt=',kdt
!     enddo
!     write(1000+me,*)' at kdt = ',kdt
!     tem = 1000.0

!     call moist_bud2(im,ix,ix,levs,me,kdt,con_g,tem,del,rain1 &
!    &,               txa, clw(1,1,2), clw(1,1,1)         &
!    &,           Stateout%gq0(1:ix,1:levs,ntrw),Stateout%gq0(1:ix,1:levs,ntsw)&
!    &,           Stateout%gq0(1:ix,1:levs,ntgl)                       &
!    &,           Stateout%gq0(1:ix,1:levs,1),Stateout%gq0(1:ix,1:levs,ntcw)   &
!    &,           Stateout%gq0(1:ix,1:levs,ntiw)                       &
!    &,           qrn, qsnw, qgl, ' m_micro  ', grid%xlon(1:im), grid%xlat(1:im))

!       if (lprnt) write(0,*) ' rain1=',rain1(ipr)*86400.0, &
!    &' rainc=',diag%rainc(ipr)*86400.0
!    &,' cn_prc=',cn_prc(ipr),' cn_snr=',cn_snr(ipr),' kdt=',kdt
!       if(lprnt) write(0,*) ' aftgt0=',Stateout%gt0(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*) ' aftlsgq0=',stateout%gq0(ipr,:,1),' kdt=',kdt
!       if (lprnt) write(0,*)' cli1aft=',stateout%gq0(ipr,:,ntiw),' kdt=',kdt
!       if (ntgl > 0 .and. lprnt)  &
!                  write(0,*)' cgw1aft=',stateout%gq0(ipr,:,ntgl),' kdt=',kdt
!       if (lprnt) write(0,*)' cloudsm=',Tbd%phy_f3d(ipr,:,1)*100,' kdt=',kdt
!       if (lprnt) write(0,*)' clw2aft=',stateout%gq0(ipr,:,ntcw),' kdt=',kdt
!       if (lprnt) write(0,*)' qrna=',qrn(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qsnwa=',qsnw(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qglba',qgl(ipr,:),' kdt=',kdt

!## CCPP ##* m_micro_interstitial.F90/m_micro_post_run

          tem = dtp * con_p001 / con_day
          if (abs(Model%fprcp) == 1 .or. mg3_as_mg2) then
            do k=1,levs
              do i=1,im
                if (abs(qrn(i,k))  < qsmall) qrn(i,k)  = zero
                if (abs(qsnw(i,k)) < qsmall) qsnw(i,k) = zero
                Stateout%gq0(i,k,ntrw)  = qrn(i,k)
                Stateout%gq0(i,k,ntsw)  = qsnw(i,k)
                Stateout%gq0(i,k,ntrnc) = ncpr(i,k)
                Stateout%gq0(i,k,ntsnc) = ncps(i,k)
              enddo
            enddo
            do i=1,im
              Diag%ice(i)  = tem * Stateout%gq0(i,1,ntiw) 
              Diag%snow(i) = tem * qsnw(i,1) 
            enddo
          elseif (Model%fprcp > 1) then
            do k=1,levs
              do i=1,im
                if (abs(qrn(i,k))  < qsmall) qrn(i,k)  = zero
                if (abs(qsnw(i,k)) < qsmall) qsnw(i,k) = zero
                if (abs(qgl(i,k))  < qsmall) qgl(i,k)  = zero
                Stateout%gq0(i,k,ntrw)  = qrn(i,k)
                Stateout%gq0(i,k,ntsw)  = qsnw(i,k)
                Stateout%gq0(i,k,ntgl)  = qgl(i,k)
                Stateout%gq0(i,k,ntrnc) = ncpr(i,k)
                Stateout%gq0(i,k,ntsnc) = ncps(i,k)
                Stateout%gq0(i,k,ntgnc) = ncgl(i,k)
              enddo
            enddo
            do i=1,im
              Diag%ice(i)     = tem * Stateout%gq0(i,1,ntiw)
              Diag%snow(i)    = tem * qsnw(i,1)
              Diag%graupel(i) = tem * qgl(i,1)
            enddo

          endif
!*## CCPP ##

!       if (lprnt) write(0,*)' cloudsm=',Tbd%phy_f3d(ipr,:,1)*100,' kdt=',kdt
!       if (lprnt) write(0,*)' clw2aft=',stateout%gq0(ipr,:,ntcw),' kdt=',kdt
!       if (lprnt) write(0,*)' qrna=',qrn(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qsnwa=',qsnw(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qglba',qgl(ipr,:),' kdt=',kdt
!

        elseif (imp_physics == Model%imp_physics_gfdl) then     ! GFDL MP
                                                                ! -------
          do i = 1, im
!## CCPP ##* Not necessary in the CCPP.
            land     (i,1)   = frland(i)
            area     (i,1)   = Grid%area(i)
!*## CCPP ##
!## CCPP ##* gfdl_cloud_microphys.F90/gfdl_cloud_microphys_run
            rain0    (i,1)   = zero
            snow0    (i,1)   = zero
            ice0     (i,1)   = zero
            graupel0 (i,1)   = zero
          enddo

          do k = 1, levs
            kk = levs-k+1
            do i = 1, im
              qn1  (i,1,k) = zero
              qv_dt(i,1,k) = zero
              ql_dt(i,1,k) = zero
              qr_dt(i,1,k) = zero
              qi_dt(i,1,k) = zero
              qs_dt(i,1,k) = zero
              qg_dt(i,1,k) = zero
              qa_dt(i,1,k) = zero
              pt_dt(i,1,k) = zero
              udt  (i,1,k) = zero
              vdt  (i,1,k) = zero
!
              qv1  (i,1,k) =  Stateout%gq0(i,kk,1)
              ql1  (i,1,k) =  Stateout%gq0(i,kk,ntcw)
              qr1  (i,1,k) =  Stateout%gq0(i,kk,ntrw)
              qi1  (i,1,k) =  Stateout%gq0(i,kk,ntiw)
              qs1  (i,1,k) =  Stateout%gq0(i,kk,ntsw)
              qg1  (i,1,k) =  Stateout%gq0(i,kk,ntgl)
              qa1  (i,1,k) =  Stateout%gq0(i,kk,ntclamt)
              pt   (i,1,k) =  Stateout%gt0(i,kk)
              w    (i,1,k) = -Statein%vvl(i,kk)*(one+con_fvirt*qv1(i,1,k))   &
                           *  Stateout%gt0(i,kk) / Statein%prsl(i,kk) * (con_rd*onebg)
              uin  (i,1,k) =  Stateout%gu0(i,kk)
              vin  (i,1,k) =  Stateout%gv0(i,kk)
              delp (i,1,k) =  del(i,kk)
              dz   (i,1,k) = (Statein%phii(i,kk)-Statein%phii(i,kk+1)) * onebg
              p123 (i,1,k) = Statein%prsl(i,kk)
              refl (i,1,k) = Diag%refl_10cm(i,kk)
            enddo
          enddo


          call gfdl_cloud_microphys_driver(qv1, ql1, qr1, qi1, qs1, qg1, qa1, &
                                           qn1, qv_dt, ql_dt, qr_dt, qi_dt,   &
                                           qs_dt, qg_dt, qa_dt, pt_dt, pt, w, &
                                           uin, vin, udt, vdt, dz, delp,      &
                                           area, dtp, land, rain0, snow0,     &
                                           ice0, graupel0, .false., .true.,   &
                                           1, im, 1, 1, 1, levs, 1, levs,     &
                                           seconds,p123,Model%lradar,refl,    &
                                           reset)
          tem = dtp * con_p001 / con_day
          do i = 1, im
!           rain0(i,1)     = max(zero, rain0(i,1))
!           snow0(i,1)     = max(zero, snow0(i,1))
!           ice0(i,1)      = max(zero, ice0(i,1))
!           graupel0(i,1)  = max(zero, graupel0(i,1))
            if (rain0(i,1)*tem < rainmin) then
               rain0(i,1) = zero
            endif
            if (ice0(i,1)*tem < rainmin) then
              ice0(i,1) = zero
            endif
            if (snow0(i,1)*tem < rainmin) then
              snow0(i,1) = zero
            endif
            if (graupel0(i,1)*tem < rainmin) then
              graupel0(i,1) = zero
            endif

            rain1(i)        = (rain0(i,1)+snow0(i,1)+ice0(i,1)+graupel0(i,1)) * tem
            Diag%ice(i)     = ice0    (i,1) * tem
            Diag%snow(i)    = snow0   (i,1) * tem
            Diag%graupel(i) = graupel0(i,1) * tem
            if ( rain1(i) > rainmin ) then
              Diag%sr(i) = (snow0(i,1) + ice0(i,1)  + graupel0(i,1)) &
                         / (rain0(i,1) + snow0(i,1) + ice0(i,1) + graupel0(i,1))
            else
              Diag%sr(i)  = zero
            endif
          enddo
#ifdef REPRO
          ! Convert rain0, ice0, graupel0 and snow0 from mm/day to m/physics-timestep
          ! for later use (approx. lines 7970, calculation of srflag)
          rain0 = tem*rain0
          ice0  = tem*ice0
          snow0 = tem*snow0
          graupel0 = tem*graupel0
#endif
          do k = 1, levs
            kk = levs-k+1
            do i=1,im
              Stateout%gq0(i,k,1   )    = qv1(i,1,kk)       + qv_dt(i,1,kk) * dtp
              Stateout%gq0(i,k,ntcw)    = ql1(i,1,kk)       + ql_dt(i,1,kk) * dtp
              Stateout%gq0(i,k,ntrw)    = qr1(i,1,kk)       + qr_dt(i,1,kk) * dtp
              Stateout%gq0(i,k,ntiw)    = qi1(i,1,kk)       + qi_dt(i,1,kk) * dtp
              Stateout%gq0(i,k,ntsw)    = qs1(i,1,kk)       + qs_dt(i,1,kk) * dtp
              Stateout%gq0(i,k,ntgl)    = qg1(i,1,kk)       + qg_dt(i,1,kk) * dtp
              Stateout%gq0(i,k,ntclamt) = qa1(i,1,kk)       + qa_dt(i,1,kk) * dtp
              Stateout%gt0(i,k)         = Stateout%gt0(i,k) + pt_dt(i,1,kk) * dtp
              Stateout%gu0(i,k)         = Stateout%gu0(i,k) + udt  (i,1,kk) * dtp
              Stateout%gv0(i,k)         = Stateout%gv0(i,k) + vdt  (i,1,kk) * dtp
              Diag%refl_10cm(i,k)       = refl(i,1,kk)
            enddo


            if (Model%effr_in) then
              do i =1, im
                den(i,k) = 0.622_kind_phys*Statein%prsl(i,k) / &
                          (con_rd*Stateout%gt0(i,k)*(Stateout%gq0(i,k,1)+0.622_kind_phys))
              enddo
            endif
          enddo
!*## CCPP ##
!## CCPP ##* maximum_hourly_diagnostics.F90/maximum_hourly_diagnsostics_run
!Calculate hourly max 1-km agl and -10C reflectivity
          if (Model%lradar .and.                              &
             (imp_physics == Model%imp_physics_gfdl .or.      &
              imp_physics == Model%imp_physics_thompson)) then
            allocate(refd(im))
            allocate(refd263k(im))
            call max_fields(Statein%phil,Diag%refl_10cm,con_g,im,levs,refd,Stateout%gt0,refd263k)
            if (reset) then
              do i=1,im
                Diag%refdmax(I)     = -35.0_kind_phys
                Diag%refdmax263k(I) = -35.0_kind_phys
              enddo
            endif
            do i=1,im
              Diag%refdmax(i)     = max(Diag%refdmax(i),refd(i))
              Diag%refdmax263k(i) = max(Diag%refdmax263k(i),refd263k(i))
            enddo
            deallocate (refd) 
            deallocate (refd263k)
          endif
!*## CCPP ##
!## CCPP ##* gfdl_cloud_microphys.F90/gfdl_cloud_microphys_run
          if(Model%effr_in) then 
            call cloud_diagnosis (1, im, 1, levs, den(1:im,1:levs),            &
                                  del(1:im,1:levs), islmsk(1:im),              &
               Stateout%gq0(1:im,1:levs,ntcw), Stateout%gq0(1:im,1:levs,ntiw), &
               Stateout%gq0(1:im,1:levs,ntrw),                                 &
               Stateout%gq0(1:im,1:levs,ntsw)+Stateout%gq0(1:im,1:levs,ntgl),  &
               Stateout%gq0(1:im,1:levs,ntgl)*0.0, Stateout%gt0(1:im,1:levs),  &
               Tbd%phy_f3d(1:im,1:levs,1),     Tbd%phy_f3d(1:im,1:levs,2),     &
               Tbd%phy_f3d(1:im,1:levs,3),     Tbd%phy_f3d(1:im,1:levs,4),     &
               Tbd%phy_f3d(1:im,1:levs,5))

!*## CCPP ##
!           do k = 1, levs
!             do i=1,im
!               if(Model%me==0) then
!                 if(Tbd%phy_f3d(i,k,1) > 5.) then
!                   write(6,*) 'phy driver:cloud radii:',Model%kdt, i,k,        &
!                              Tbd%phy_f3d(i,k,1)
!                 endif
!                if(Tbd%phy_f3d(i,k,3)> zero) then
!                   write(6,*) 'phy driver:rain radii:',Model%kdt, i,k,         &
!                              Tbd%phy_f3d(i,k,3)
!                 endif
!
!               endif
!             enddo
!           enddo

          endif

        endif  ! end of if(Model%imp_physics)
      endif    ! end if_ncld

!     if (lprnt) write(0,*)' rain1 after ls=',rain1(ipr)
!
      if (Model%cscnv .and. Model%do_aw) then
!  Arakawa-Wu adjustment of large-scale microphysics tendencies:
!  reduce by factor of (1-sigma)
!  these are microphysics increments. We want to keep (1-sigma) of the increment,
!  we will remove sigma*increment from final values
!         fsigma = zero  ! don't apply any AW correction, in addition comment next line
!         fsigma = sigmafrac

!  adjust sfc rainrate for conservation
!  vertically integrate reduction of water increments, reduce precip by that amount
!## CCPP ##* cs_conv_aw_adj.F90/cs_conv_aw_adj_run Note: The conditional above 
! is not checked in the scheme, so the control of the code below is through its 
! inclusion in a CCPP SDF

        temrain1(:) = zero
        do k = 1,levs
          do i = 1,im
            tem1               = sigmafrac(i,k)
            Stateout%gt0(i,k)  = Stateout%gt0(i,k) - tem1 * (Stateout%gt0(i,k)-dtdt(i,k))
            tem2                = tem1 * (Stateout%gq0(i,k,1)-dqdt(i,k,1))
            Stateout%gq0(i,k,1) = Stateout%gq0(i,k,1) - tem2
            temrain1(i) = temrain1(i) - (Statein%prsi(i,k)-Statein%prsi(i,k+1)) &
                                      * tem2 * onebg
          enddo
        enddo
! add convective clouds if shoc is true and not MG microphysics
        if (Model%do_shoc .and. imp_physics /= Model%imp_physics_mg) then
          do k = 1,levs
            do i = 1,im
              Tbd%phy_f3d(i,k,ntot3d-2) = min(one, Tbd%phy_f3d(i,k,ntot3d-2)    &
     &                                             + sigmafrac(i,k))
            enddo
          enddo
        endif

!     if (lprnt) write(0,*)' gt0aftpraw=',Stateout%gt0(ipr,:),' kdt=',kdt,'me=',me
        do n=ntcw,ntcw+nncl-1
          do k = 1,levs
            do i = 1,im
              tem1                = sigmafrac(i,k) * (Stateout%gq0(i,k,n)-dqdt(i,k,n))
              Stateout%gq0(i,k,n) = Stateout%gq0(i,k,n) - tem1
              temrain1(i)         = temrain1(i) - (Statein%prsi(i,k)-Statein%prsi(i,k+1)) &
                                        * tem1 * onebg
            enddo
          enddo
        enddo
!     write(1000+me,*)' rain1=',rain1(4),' temrain1=',temrain1(i)*0.001
        do i = 1,im
          rain1(i) = max(rain1(i) - temrain1(i)*con_p001, zero)
        enddo
      endif

!*## CCPP ##
!## CCPP ##* GFS_MP_generic.F90/GFS_MP_generic_post_run
      Diag%rain(:) = Diag%rainc(:) + frain * rain1(:)  ! total rain per timestep

!  ---  get the amount of different precip type for Noah MP
!  ---  convert from m/dtp to mm/s
      if (Model%lsm==Model%lsm_noahmp) then
        if (Model%imp_physics == Model%imp_physics_mg .or. &
            Model%imp_physics == Model%imp_physics_gfdl) then
          !GJF: Should all precipitation rates have the same denominator below? 
          ! It appears that Diag%rain and Diag%rainc are on the dynamics time step,
          ! but Diag%snow,graupel,ice are on the physics time step? This doesn't
          ! matter as long as dtp=dtf (frain=1).
          tem = one / (dtp*con_p001)
          Sfcprop%draincprv(:)   = tem * Diag%rainc(:)
          Sfcprop%drainncprv(:)  = tem * (frain * rain1(:))
          Sfcprop%dsnowprv(:)    = tem * Diag%snow(:)
          Sfcprop%dgraupelprv(:) = tem * Diag%graupel(:)
          Sfcprop%diceprv(:)     = tem * Diag%ice(:)
        else
          Sfcprop%draincprv(:)   = zero
          Sfcprop%drainncprv(:)  = zero
          Sfcprop%dsnowprv(:)    = zero
          Sfcprop%dgraupelprv(:) = zero
          Sfcprop%diceprv(:)     = zero
        endif
      end if !  if (Model%lsm == Model%lsm_noahmp)
      
      if (Model%cal_pre) then       ! hchuang: add dominant precipitation type algorithm
!
        call calpreciptype (kdt, Model%nrcm, im, ix, levs, levs+1,           &
                            Tbd%rann, Grid%xlat, Grid%xlon, Stateout%gt0,    &
                            Stateout%gq0, Statein%prsl, Statein%prsi,        &
                            Diag%rain, Statein%phii, Sfcprop%tsfc,           &  !input
                            domr, domzr, domip, doms)                           ! output
!
!        if (lprnt) print*,'debug calpreciptype: DOMR,DOMZR,DOMIP,DOMS '
!     &,DOMR(ipr),DOMZR(ipr),DOMIP(ipr),DOMS(ipr)
!        do i=1,im
!         if (abs(xlon(i)*rad2dg-114.0) .lt. 0.2  .and.
!     &       abs(xlat(i)*rad2dg- 40.0) .lt. 0.2)
!     &    print*,'debug calpreciptype: DOMR,DOMZR,DOMIP,DOMS ',
!     &    DOMR(i),DOMZR(i),DOMIP(i),DOMS(i)
!       enddo
!       HCHUANG: use new precipitation type to decide snow flag for LSM snow accumulation

        if (Model%imp_physics /= Model%imp_physics_gfdl) then
          do i=1,im
            Sfcprop%tprcp(i)  = max(zero, Diag%rain(i) )
            if(doms(i) > zero .or. domip(i) > zero) then
              Sfcprop%srflag(i) = one
            else
              Sfcprop%srflag(i) = zero
            endif
          enddo
        endif
        if (Model%lssav) then
          do i=1,im
            Diag%tdomr(i)  = Diag%tdomr(i)  + domr(i)  * dtf
            Diag%tdomzr(i) = Diag%tdomzr(i) + domzr(i) * dtf
            Diag%tdomip(i) = Diag%tdomip(i) + domip(i) * dtf
            Diag%tdoms(i)  = Diag%tdoms(i)  + doms(i)  * dtf
          enddo
        endif

      endif

!*## CCPP ##
!## CCPP ##* this block not yet in CCPP
!--------------------------------
! vay-2018 for Dycore-Tendencies save Stateout%X => Diag%dX3dt_cgw
!
      if (ldiag_ugwp) then
        Diag%dt3dt_cgw = Stateout%gt0
        Diag%dv3dt_cgw = Stateout%gv0
        Diag%du3dt_cgw = Stateout%gu0
      endif
!--------------------------------
!*## CCPP ##
!## CCPP ##* GFS_MP_generic.F90/GFS_MP_generic_post_run

!  --- ...  estimate t850 for rain-snow decision

      t850(1:im) = Stateout%gt0(1:im,1)

      do k = 1, levs-1
        do i = 1, im
          if (Statein%prsl(i,k) > p850 .and. Statein%prsl(i,k+1) <= p850) then
            t850(i) = Stateout%gt0(i,k) - (Statein%prsl(i,k)-p850) / &
                     (Statein%prsl(i,k)-Statein%prsl(i,k+1)) *       &
                     (Stateout%gt0(i,k)-Stateout%gt0(i,k+1))
          endif
        enddo
      enddo

      if (Model%imp_physics == Model%imp_physics_gfdl) then         ! GFDL microphysics
                                                                    ! -----------------
! determine convective rain/snow by surface temperature
! determine large-scale rain/snow by rain/snow coming out directly from MP
        tem = dtp * con_p001 / con_day
        do i = 1, im
          Sfcprop%tprcp(i)  = max(zero, Diag%rain(i) )! clu: rain -> tprcp
          Sfcprop%srflag(i) = zero                    ! clu: default srflag as 'rain' (i.e. 0)
          if (Sfcprop%tsfc(i) >= 273.15_kind_phys) then
            crain = Diag%rainc(i)
            csnow = zero
          else
            crain = zero
            csnow = Diag%rainc(i)
          endif
!         if (snow0(i,1)+ice0(i,1)+graupel0(i,1)+csnow > rain0(i,1)+crain) then
!          if (snow0(i,1)+ice0(i,1)+graupel0(i,1)+csnow > zero) then
!            Sfcprop%srflag(i) = one                   ! clu: set srflag to 'snow' (i.e. 1)
!          endif
! compute fractional srflag
#ifdef REPRO
          ! For bit-for-bit identical results with CCPP code, snow0/ice0/graupel0/rain0
          ! were converted from mm per day to m per physics timestep previously in the code
          total_precip = snow0(i,1)+ice0(i,1)+graupel0(i,1)+rain0(i,1)+Diag%rainc(i)
          if (total_precip > rainmin) then
            Sfcprop%srflag(i) = (snow0(i,1)+ice0(i,1)+graupel0(i,1)+csnow)/total_precip
          endif
#else
          tem1 = snow0(i,1)+ice0(i,1)+graupel0(i,1)
          total_precip = (tem1+rain0(i,1)) * tem + Diag%rainc(i)
          if (total_precip > rainmin) then
            Sfcprop%srflag(i) = (tem1*tem+csnow) / total_precip
          endif
#endif
        enddo
      elseif( .not. Model%cal_pre) then
        if (Model%imp_physics == Model%imp_physics_mg) then         ! MG microphysics
                                                                    ! ---------------
          do i=1,im
            if (Diag%rain(i) > rainmin) then
              tem1 = max(zero, (Diag%rain(i)-Diag%rainc(i))) * Diag%sr(i)
              tem2 = one / Diag%rain(i)
              if (t850(i) > 273.16_kind_phys) then
                Sfcprop%srflag(i) = max(zero, min(one, tem1*tem2))
              else
                Sfcprop%srflag(i) = max(zero, min(one, (tem1+Diag%rainc(i))*tem2))
              endif
            else
              Sfcprop%srflag(i) = zero
              Diag%rain(i)      = zero
              Diag%rainc(i)     = zero
            endif
            Sfcprop%tprcp(i)  = max(zero, Diag%rain(i))
          enddo
        else                                                        ! not GFDL or MG microphysics
                                                                    ! ---------------------------
          do i = 1, im
            Sfcprop%tprcp(i)  = max(zero, Diag%rain(i))
            Sfcprop%srflag(i) = Diag%sr(i)
          enddo
        endif
      endif

      if (Model%lssav) then
!        if (Model%me == 0) print *,'in phys drive, kdt=',Model%kdt, &
!          'totprcpb=', Diag%totprcpb(1),'totprcp=',Diag%totprcp(1), &
!          'rain=',Diag%rain(1)
        do i=1,im
          Diag%cnvprcp(i)  = Diag%cnvprcp(i)  + Diag%rainc(i)
          Diag%totprcp (i) = Diag%totprcp (i) + Diag%rain(i)
          Diag%totice  (i) = Diag%totice  (i) + Diag%ice(i)
          Diag%totsnw  (i) = Diag%totsnw  (i) + Diag%snow(i)
          Diag%totgrp  (i) = Diag%totgrp  (i) + Diag%graupel(i)
!
          Diag%cnvprcpb(i) = Diag%cnvprcpb(i) + Diag%rainc(i)
          Diag%totprcpb(i) = Diag%totprcpb(i) + Diag%rain(i)
          Diag%toticeb (i) = Diag%toticeb (i) + Diag%ice(i)
          Diag%totsnwb (i) = Diag%totsnwb (i) + Diag%snow(i)
          Diag%totgrpb (i) = Diag%totgrpb (i) + Diag%graupel(i)
        enddo

        if (Model%ldiag3d) then
          do k=1,levs
            do i=1,im
              Diag%dt3dt(i,k,6) = Diag%dt3dt(i,k,6) + (Stateout%gt0(i,k)-dtdt(i,k)) * frain
!             Diag%dq3dt(i,k,4) = Diag%dq3dt(i,k,4) + (Stateout%gq0(i,k,1)-dqdt(i,k,1)) * frain
            enddo
          enddo
        endif
      endif

!  --- ...  coupling insertion

      if (Model%cplflx .or. Model%cplchm) then
        do i = 1, im
          Tbd%dsnow_cpl(i)= max(zero, Diag%rain(i) * Sfcprop%srflag(i))
          Tbd%drain_cpl(i)= max(zero, Diag%rain(i) - Tbd%dsnow_cpl(i))
          Coupling%rain_cpl(i) = Coupling%rain_cpl(i) + Tbd%drain_cpl(i)
          Coupling%snow_cpl(i) = Coupling%snow_cpl(i) + Tbd%dsnow_cpl(i)
        enddo
      endif

      if (Model%cplchm) then
        do i = 1, im
          Coupling%rainc_cpl(i) = Coupling%rainc_cpl(i) + Diag%rainc(i)
        enddo
      endif
!*## CCPP ##
!  --- ...  end coupling insertion

!## CCPP ##* GFS_surface_generic.F90/GFS_surface_generic_post_run
!  --- ...  total runoff is composed of drainage into water table and
!           runoff at the surface and is accumulated in unit of meters
      if (Model%lssav) then
        do i=1,im
          Diag%runoff(i)  = Diag%runoff(i)  + (drain(i)+runof(i)) * dtf
          Diag%srunoff(i) = Diag%srunoff(i) + runof(i) * dtf
        enddo
      endif
!*## CCPP ##
!## CCPP ##* This block is not in the CCPP because data transfer between global 
! and local variables is not necessary in the CCPP.
!  --- ...  return updated smsoil and stsoil to global arrays
      if (Model%frac_grid) then
        do k=1,lsoil
          do i=1,im
            if (dry(i)) then
              Sfcprop%smc(i,k) = smsoil(i,k)
              Sfcprop%stc(i,k) = stsoil(i,k)
              Sfcprop%slc(i,k) = slsoil(i,k)
            endif
          enddo
        enddo
      else
        do k=1,lsoil
          do i=1,im
            Sfcprop%smc(i,k) = smsoil(i,k)
            Sfcprop%stc(i,k) = stsoil(i,k)
            Sfcprop%slc(i,k) = slsoil(i,k)
          enddo
        enddo
      endif
!*## CCPP ##

!  --- ...  calculate column precipitable water "pwat"
      Diag%pwat(:) = zero
      do k = 1, levs
        do i=1,im
          work1(i) = zero
        enddo
        if (ncld > 0) then
          do ic = ntcw, ntcw+nncl-1
            do i=1,im
              work1(i) = work1(i) + Stateout%gq0(i,k,ic)
            enddo
          enddo
        endif
        do i=1,im
          Diag%pwat(i) = Diag%pwat(i) + del(i,k)*(Stateout%gq0(i,k,1)+work1(i))
        enddo
!     if (lprnt .and. i == ipr) write(0,*)' gq0=',
!    &gq0(i,k,1),' qgrs=',qgrs(i,k,1),' work2=',work2(i),' k=',k
      enddo
      do i=1,im
        Diag%pwat(i) = Diag%pwat(i) * onebg
      enddo

!     tem = dtf * 0.03456 / 86400.0
!       write(1000+me,*)' pwat=',pwat(i),'i=',i,',
!    &' rain=',rain(i)*1000.0,' dqsfc1=',dqsfc1(i)*hefac(i)*tem,' kdt=',kdt
!    &,' e-p=',dqsfc1(i)*hefac(i)*tem-rain(i)*1000.0
!     if (lprnt) write(0,*)' pwat=',pwat(ipr),',
!    &' rain=',rain(ipr)*1000.0,' dqsfc1=',dqsfc1(ipr)*hefac(ipr)*tem,' kdt=',kdt
!    &,' e-p=',dqsfc1(ipr)*hefac(ipr)*tem-rain(ipr)*1000.0

!
!     if (lprnt .and. rain(ipr) > 5) call mpi_quit(5678)
!     if (lat == 45) write(1000+me,*)' pwat=',pwat(1),' kdt=',kdt
!       if (lprnt) then
!         write(7000,*) ' endgu0=',gu0(ipr,:),' kdt=',kdt
!         write(7000,*) ' endgv0=',gv0(ipr,:),' kdt=',kdt,' nnp=',nnp
!         write(0,*) ' endgt0=',Stateout%gt0(ipr,:),' kdt=',kdt
!         write(0,*) ' endgq0=',Stateout%gq0(ipr,:,1),' kdt=',kdt
!         write(0,*) ' endgw0=',gq0(ipr,:,3),' kdt=',kdt,' lat=',lat
!         write(0,*) ' endzorl=',Sfcprop%zorl(ipr),' kdt=',kdt
!       endif

      if (Model%do_sppt .or. Model%ca_global)then
!--- radiation heating rate
        Tbd%dtdtr(1:im,:) = Tbd%dtdtr(1:im,:) + dtdtc(1:im,:)*dtf
      endif
!*## CCPP ##
!## CCPP ##* This block is not in the CCPP since it is not needed in the CCPP.
                            deallocate (clw)
      if (allocated(cnvc))  deallocate(cnvc)
      if (allocated(cnvw))  deallocate(cnvw)
      if (allocated(qrn))   deallocate(qrn)
      if (allocated(qsnw))  deallocate(qsnw)
      if (allocated(qgl))   deallocate(qgl)
      if (allocated(ncpl))  deallocate(ncpl)
      if (allocated(ncpi))  deallocate(ncpi)
      if (allocated(ncpr))  deallocate(ncpr)
      if (allocated(ncps))  deallocate(ncps)
      if (allocated(ncgl))  deallocate(ncgl)

      if (allocated(liq0))  deallocate(liq0)
      if (allocated(ice00)) deallocate(ice00)


!     deallocate (fscav, fswtr)
!
!     if (lprnt) write(0,*)' end of gbphys maxu=',
!    &maxval(gu0(1:im,1:levs)),' minu=',minval(gu0(1:im,1:levs))
!    &,' maxv=',maxval(gv0(1:im,1:levs)),' minv=',
!    & minval(gv0(1:im,1:levs)),' kdt=',kdt,' lat=',lat,' nnp=',nnp
!     if (lprnt) write(0,*)' end of gbphys gv0=',gv0(:,120:128)
!     if (lprnt) write(0,*)' end of gbphys at kdt=',kdt,&
!    &' rain=',rain(ipr),' rainc=',rainc(ipr)
!     if (lprnt) call mpi_quit(7)
!     if (kdt > 2 ) call mpi_quit(70)
!    if (lprnt) then
!      write(0,*)' at the end of physics kdt=',kdt
!      write(0,*)' end rain=',diag%rain(ipr),' rainc=',diag%rainc(ipr)
!      write(0,*)'qt0out=',Stateout%gt0(ipr,:)
!      write(0,*)'qq0outv=',Stateout%gq0(ipr,:,1)
!      write(0,*)'qq0outw=',Stateout%gq0(ipr,:,ntcw)
!      write(0,*)'qq0outi=',Stateout%gq0(ipr,:,ntiw)
!      write(0,*)'qq0outo=',Stateout%gq0(ipr,:,ntoz)
!    endif
!    if (lprnt) write(0,*)'gq0outtke=',Stateout%gq0(ipr,1:25,ntke)    &
!      ,'xlon=',grid%xlon(ipr)*rad2dg,' xlat=',grid%xlat(ipr)*rad2dg
!     if (lprnt) write(0,*)' clouddriverend=',Tbd%phy_f3d(ipr,:,1)*100,' kdt=',kdt

!     deallocate (qlcn, qicn, w_upi, cf_upi, CNV_MFD, CNV_PRC3, &
      deallocate (qlcn, qicn, w_upi, cf_upi, CNV_MFD,           &
                  CNV_DQLDT, clcn, cnv_fice, cnv_ndrop, cnv_nice)
      if (imp_physics == Model%imp_physics_gfdl) then
        deallocate (delp,  dz,    uin,   vin,   pt,    qv1,   ql1, qr1,        &
                    qg1,   qa1,   qn1,   qi1,   qs1,   pt_dt, qa_dt, udt, vdt, &
                    w,     qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt,p123,refl)
        deallocate (den)
      endif
!*## CCPP ##
!## CCPP ##* maximum_hourly_diagnostics.F90/maximum_hourly_diagnostics_run
      if (allocated(tke)) deallocate (tke)
      if (Model%cscnv) then
        deallocate (sigmatot, sigmafrac)
      endif
!
      if (reset) then
        do i=1, im
! find max hourly wind speed then decompose
          Diag%spd10max(i) = -999.0_kind_phys
          Diag%u10max(i)   = -999.0_kind_phys
          Diag%v10max(i)   = -999.0_kind_phys
          Diag%t02max(i)   = -999.0_kind_phys
          Diag%t02min(i)   =  999.0_kind_phys
          Diag%rh02max(i)  = -999.0_kind_phys
          Diag%rh02min(i)  =  999.0_kind_phys
        enddo
      endif
      do i=1, im
! find max hourly wind speed then decompose
        tem = sqrt(Diag%u10m(i)*Diag%u10m(i) + Diag%v10m(i)*Diag%v10m(i))
        if (tem > Diag%spd10max(i)) then
          Diag%spd10max(i) = tem
          Diag%u10max(i)   = Diag%u10m(i)
          Diag%v10max(i)   = Diag%v10m(i)
        endif
        pshltr = Statein%pgr(i)*exp(-0.068283_kind_phys/Stateout%gt0(i,1))
        QCQ    = PQ0/pshltr*EXP(A2A*(Sfcprop%t2m(i)-A3)/(Sfcprop%t2m(i)-A4))
        rh02   = Sfcprop%q2m(i) / QCQ
        IF (rh02 > one) THEN
          rh02 = one
        ENDIF
        IF (rh02 < RHmin) THEN  !use smaller RH limit for stratosphere
          rh02 = RHmin
        ENDIF
        Diag%rh02max(i) = max(Diag%rh02max(i), rh02)
        Diag%rh02min(i) = min(Diag%rh02min(i), rh02)
        Diag%T02MAX(I)  = MAX(Diag%T02MAX(I), Sfcprop%t2m(i))  !<--- Hourly max 2m T
        Diag%T02MIN(I)  = MIN(Diag%T02MIN(I), Sfcprop%t2m(i))  !<--- Hourly min 2m T
      enddo
!*## CCPP ##
!     if (kdt > 2 ) stop

!       if (Model%nstf_name(1) > 0) then
!          if (lprnt) write(0,*)' end driver sfcprop%tref=',Sfcprop%tref(ipr),' kdt=',kdt
!       endif
!       if (Model%frac_grid) then
!          if (lprnt) write(0,*)' end driver sfcprop%tsfcl=',Sfcprop%tsfcl(ipr),' kdt=',kdt
!          if (lprnt) write(0,*)' end driver sfcprop%tsfco=',Sfcprop%tsfco(ipr),' kdt=',kdt
!          if (lprnt) write(0,*)' end driver sfcprop%tisfc=',Sfcprop%tisfc(ipr),' kdt=',kdt
!       endif

      return
!...................................
      end subroutine GFS_physics_driver
!-----------------------------------


      subroutine max_fields(phil,ref3D,grav,im,levs,refd,tk,refd263k)
      use machine, only : kind_phys
      integer, intent(in)               :: im,levs
      real (kind=kind_phys), intent(in) :: grav
      real (kind=kind_phys), intent(in),dimension(im,levs)  :: phil,ref3D,tk
      integer                  :: i,k,ll,ipt,kpt
      real                     :: dbz1avg,zmidp1,zmidloc,refl,fact
      real, dimension(im,levs) :: z
      real, dimension(im)      :: zintsfc
      real, dimension(im), intent(inout) :: refd,refd263k
      REAL :: dbz1(2),dbzk,dbzk1
      logical counter
      do i=1,im
         do k=1,levs
            z(i,k) = phil(i,k)/grav
         enddo
      enddo
      do i=1,im
         refd(I) = -35.
  vloop:  do k=1,levs-1
            if ( z(i,k+1) >= 1000. .and. z(i,k) <= 1000.)  then
               zmidp1  = z(i,k+1)
               zmidLOC = z(i,k)
               dbz1(1) = ref3d(i,k+1)   !- dBZ (not Z) values
               dbz1(2) = ref3d(i,k)     !- dBZ values
               exit vloop
            endif
         enddo vloop

!!! Initial curefl value without reduction above freezing level
!
!         curefl=0.
!         if (cprate(i,j)>0.) then
!           cuprate=rdtphs*cprate(i,j)
!           curefl=cu_a*cuprate**cu_b
!         endif
         do ll=1,2
           refl=0.
           if (dbz1(ll)>-35.) refl=10.**(0.1*dbz1(ll))
!           dbz1(l)=curefl+refl    !- in Z units
             dbz1(ll)=refl
         enddo
!-- Vertical interpolation of Z (units of mm**6/m**3)
         fact=(1000.-zmidloc)/(zmidloc-zmidp1)
         dbz1avg=dbz1(2)+(dbz1(2)-dbz1(1))*fact
!-- Convert to dBZ (10*logZ) as the last step
         if (dbz1avg>0.01) then
           dbz1avg=10.*alog10(dbz1avg)
         else
           dbz1avg=-35.
         endif
         refd(I)=max(refd(I),dbz1avg)
      enddo

!-- refl at -10C
      do i=1,im
         dbz1(1) = -35.
         dbz1(2) = -35.
  vloopm10:  do k=1,levs-1
            if (tk(i,k+1) .le. 263.15 .and. tk(i,k) .ge. 263.15)  then     
               dbz1(1)=ref3d(i,k+1)   !- dBZ (not Z) values
               dbz1(2)=ref3d(i,k) !- dBZ values
               exit vloopm10
            endif
         enddo vloopm10
         
         do ll=1,2
           refl=0.
           if (dbz1(ll)>-35.) refl=10.**(0.1*dbz1(ll))
!           dbz1(l)=curefl+refl    !- in Z units
             dbz1(ll)=refl
         enddo
!-- Take max of bounding reflectivity values 
         dbz1avg=maxval(dbz1)
!-- Convert to dBZ (10*logZ) as the last step
         if (dbz1avg>0.01) then
           dbz1avg=10.*alog10(dbz1avg)
         else
           dbz1avg=-35.
         endif
         refd263K(I)=dbz1avg
      enddo
      end subroutine max_fields

 subroutine moist_bud(im,ix,ix2,levs,me,kdt,grav,dtp,delp,rain, &
                           qv0,ql0,qi0,qv1,ql1,qi1,comp, xlon, xlat)
!  nov 2016 - S. Moorthi - routine to compute local moisture budget
      use machine, only : kind_phys
      implicit none
      character*10          :: comp
      integer               :: im,ix,ix2,levs,me,kdt
      real (kind=kind_phys) :: grav, rain(im), dtp, xlon(im), xlat(im)
      real (kind=kind_phys), dimension(ix,levs)  :: qv0,ql0,qi0,delp
      real (kind=kind_phys), dimension(ix2,levs) :: qv1,ql1,qi1
      REAL (kind=kind_phys), dimension(im) :: sumq, sumqv, sumql, sumqi
      integer               :: i, k
!
      do i=1,im
        sumqv(i) = 0.0_kind_phys
        sumql(i) = 0.0_kind_phys
        sumqi(i) = 0.0_kind_phys
        sumq (i) = 0.0_kind_phys
      enddo
      do k=1,levs
        do i=1,im
          sumqv(i) = sumqv(i) + (qv1(i,k) - qv0(i,k)) * delp(i,k)
          sumql(i) = sumql(i) + (ql1(i,k) - ql0(i,k)) * delp(i,k)
          sumqi(i) = sumqi(i) + (qi1(i,k) - qi0(i,k)) * delp(i,k)
        enddo
      enddo
      do i=1,im
        sumqv(i) = - sumqv(i) * (1.0_kind_phys/grav)
        sumql(i) = - sumql(i) * (1.0_kind_phys/grav)
        sumqi(i) = - sumqi(i) * (1.0_kind_phys/grav)
        sumq (i) =  sumqv(i) + sumql(i) + sumqi(i)
      enddo
      do i=1,im
        write(2000+me,*)' in moist_bud:',' i=',i,' sumq=',sumq(i), &
       ' sumqv=',sumqv(i),' sumql=',sumql(i),' sumqi=',sumqi(i),   &
       ' rain=',rain(i)*dtp,' kdt=',kdt,' component=',trim(comp),  &
       ' qv:=',qv1(i,1),qv0(i,1),' ql=',ql1(i,1),ql0(i,1),         &
       ' qi=',qi1(i,1), qi0(i,1),' xlon=',xlon(i),' xlat=',xlat(i)
      enddo
      return

      end subroutine moist_bud


      subroutine moist_bud2(im,ix,ix2,levs,me,kdt,grav,dtp,delp,rain, &
                            qv0,ql0,qi0,qr0,qs0,qg0,                  &
                            qv1,ql1,qi1,qr1,qs1,qg1,comp,xlon,xlat)
!  aug 2018 - S. Moorthi - routine to compute local moisture budget
      use machine, only : kind_phys
      implicit none
      character*10          :: comp
      integer               :: im,ix,ix2,levs,me,kdt
      real (kind=kind_phys) :: grav, rain(im), dtp, oneog, xlon(im), xlat(im)
      real (kind=kind_phys), dimension(ix,levs)  :: qv0,ql0,qi0,delp, &
                                                    qr0,qs0,qg0
      real (kind=kind_phys), dimension(ix2,levs) :: qv1,ql1,qi1,      &
                                                    qr1,qs1,qg1
      REAL (kind=kind_phys), dimension(im) :: sumq, sumqv, sumql, sumqi, &
                                              sumqr, sumqs, sumqg
      integer               :: i, k
!
      do i=1,im
        sumqv(i) = 0.0_kind_phys
        sumql(i) = 0.0_kind_phys
        sumqi(i) = 0.0_kind_phys
        sumqr(i) = 0.0_kind_phys
        sumqs(i) = 0.0_kind_phys
        sumqg(i) = 0.0_kind_phys
        sumq (i) = 0.0_kind_phys
      enddo
      do k=1,levs
        do i=1,im
          sumqv(i) = sumqv(i) + (qv1(i,k) - qv0(i,k)) * delp(i,k)
          sumql(i) = sumql(i) + (ql1(i,k) - ql0(i,k)) * delp(i,k)
          sumqi(i) = sumqi(i) + (qi1(i,k) - qi0(i,k)) * delp(i,k)
          sumqr(i) = sumqr(i) + (qr1(i,k) - qr0(i,k)) * delp(i,k)
          sumqs(i) = sumqs(i) + (qs1(i,k) - qs0(i,k)) * delp(i,k)
          sumqg(i) = sumqg(i) + (qg1(i,k) - qg0(i,k)) * delp(i,k)
        enddo
      enddo
      oneog = 1.0_kind_phys / grav
      do i=1,im
        sumqv(i) = - sumqv(i) * oneog
        sumql(i) = - sumql(i) * oneog
        sumqi(i) = - sumqi(i) * oneog
        sumqr(i) = - sumqr(i) * oneog
        sumqs(i) = - sumqs(i) * oneog
        sumqg(i) = - sumqg(i) * oneog
        sumq (i) =  sumqv(i) + sumql(i) + sumqi(i) + sumqr(i) &
                 +  sumqs(i) + sumqg(i)
      enddo
      do i=1,im
        write(1000+me,*)' in moist_bud:',' i=',i,' sumq=',sumq(i), &
       ' sumqv=',sumqv(i),' sumql=',sumql(i),' sumqi=',sumqi(i),   &
       ' sumqr=',sumqr(i),' sumqs=',sumqs(i),' sumqg=',sumqg(i),   &
       ' rain=',rain(i)*dtp,' kdt=',kdt,' component=',trim(comp),  &
       ' qv:=',qv1(i,1),qv0(i,1),' ql=',ql1(i,1),ql0(i,1),         &
       ' qi=',qi1(i,1), qi0(i,1),' qr=',qr1(i,1),qr0(i,1),         &
       ' qs=',qs1(i,1), qs0(i,1),' qg=',qg1(i,1),qg0(i,1),         &
       ' xlon=',xlon(i),' xlat=',xlat(i)
      enddo
      return

      end subroutine moist_bud2


!> @}

end module module_physics_driver

