      module idea_solar
!---------------------------------------------------------------------------
! hold effuv,effeuv ro (density (kg/m3)),nps (start pressure levels index)
! Apr 06 2012   Henry Juang, initial implement for nems
! Oct 20 2015   Weiyu Yang - add the f107 and kp inputted data.
!---------------------------------------------------------------------------
      use machine,          only : kind_phys
      use idea_composition, only : pr=> pr_idea
      use idea_composition, only : amo, amo2, amn2, amno, pi

      use physcons, only         : rgas => con_rgas
      use physcons, only         : avgd => con_avgd
!      use wam_f107_kp_mod, only: f107, kp, kdt_3h
      implicit none
!
      real, parameter   ::  PCC=1.985E-25    !  wavelength/energy conversion SI E=hc/lamda
      real, parameter   ::  eccentric=1.    !  eccentricity of Earth's orbit
! 
!     real   :: f107 = 100., f107a = f107
!
      integer,  parameter :: NWAVES = 37
      integer,  parameter :: NWAVES_EUV = 22
      integer,  parameter :: NWAVES_SRC =NWAVES-NWAVES
      integer,  parameter :: lyman_a_num  = NWAVES+1-12   
      integer,  parameter :: nsp_euv = NWAVES
      real                :: euv37(nsp_euv)
!
!TIEGCM
!
      real, dimension(nwaves) :: sfmin, afac,rlmeuv
      real, dimension(nwaves) :: csao, csao2, csan2, csao3
      real, dimension(nwaves) :: csio, csio2, csin2
      real, dimension(nwaves) :: csdo2, csdeo2
      real, dimension(nwaves) :: rwpcc
!
! Efficiencies
!
      real (kind=kind_phys), allocatable ::
     &       SRBEFF(:), o2_scale_factor(:) ! Jo2-TIEGCM
!
! middle & upper atm-re heating efficiencies, effuv, effeuv(SRC)
!
      real (kind=kind_phys), allocatable :: effuv(:), effeuv(:)
      real (kind=kind_phys), allocatable :: eff_hart(:), eff_hugg(:)
      real (kind=kind_phys), allocatable :: eff_chap(:), eff_herz(:)
      real (kind=kind_phys), allocatable :: eff_srb(:), eff_src(:)
      real (kind=kind_phys), allocatable :: eff_lya(:)
!
! merging scheme for heating rates
!
      real, parameter     :: xb=7.5, xt=8.5            ! for Hp = 7 km:  52.5 km < Z_logp < 59.5 km 
      real, parameter     :: xbl= 0.99*xb              ! xbl < xb
      real, parameter     :: rdx=1./(xt-xb)
      real, parameter     :: xlogps = 11.5129          ! alog(1.e5=Ps_in_Pa)
      real, parameter     :: prdot02 = 1.e3*exp(-xbl)  ! mbars because pr = pr_idea in (mb)
      integer             :: nps                       ! layer where Pressure < 0.02   (2Pa)
                                                       ! nps-pressure index to start WAM-solar/photo <= 52.5 km 
      end module idea_solar
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine idea_solar_init(mpi_id, master, levs)
!SK   subroutine idea_solar_init(levs)
!----------------------------------------------------------------------------
!
! define: effuv,effeuv, solar fluxes (37), wlengths, cross-sections,
!     other eff_hart, eff_hugg, eff_chap, eff_herz, eff_srb, eff_src, eff_lya
!     init solar_calendar and F107/Kp datasets
!
!----------------------------------------------------------------------------    
      use machine,            only : kind_phys
      use idea_wam_control,   only : SPW_DRIVERS
      use IDEA_IO_UNITS    ,  only : nml_solar, nlun_solar, ch100
!
      use idea_solarno_input, only : solar_readno_snoewx

      use idea_solar_input,  only : idea_solar_fix, itheia
      use idea_solar_input,  only : wf107_s,   wkp_s
      use idea_solar_input,  only : wf107a_s,  wap_s
      use idea_solar_input,  only : weuv_s,  nwafix

      use idea_solar_input,  only : solar_read_namelist
      use idea_solar_input,  only : solar_read_myy1947_2016, 
     &                              solar_read_wam_init
      use idea_solar_input,  only : solar_waccmx_advance
      use idea_solar_input,  only : solar_wamstep_advance

!      use idea_solar_input,  only : f107 => wf107_s,   kp => wkp_s
!      use idea_solar_input,  only : f107a => wf107a_s, Ap => wap_s
!      use idea_solar_input,  only : EUV37 => weuv_s,  nsp_euv => nwafix
!
!SK   use idea_mpi_def,       only : mpi_id, mpi_err, MPI_COMM_ALL,info
!
      use idea_solar,         only : nps,  effuv, effeuv, pr, prdot02
      use idea_solar,         only : o2_scale_factor,  SRBEFF
!
      use idea_solar,      only : eff_hart, eff_hugg, eff_chap, eff_herz
      use idea_solar,      only : eff_srb, eff_src, eff_lya
!
      use idea_solar,       only  : rwpcc, pcc
      use idea_solar,       only  : csao, csao2, csao3, csan2
      use idea_solar,       only  : csio, csio2, csin2
      use idea_solar,       only  : csdo2, csdeo2
      use idea_solar,       only  : eccentric
      use idea_solar,       only  : nwaves,NWAVES_EUV,lyman_a_num,
     &       nwaves_src
      use idea_solar,       only  : sfmin, afac,rlmeuv
! 
      use wam_date_calendar, only : COPY_IDAT_NEMS_2_WAM, idat_wam,
     &            irhour_wam
      use wam_date_calendar, only : CURRENT_NCEP_JDAT, ndwam  
      use wam_date_calendar, only : curday_wam, curmonth_wam, curddd_wam 
      use wam_date_calendar, only : curyear_wam, curutsec_wam 
!
      use wam_f107_kp_mod,    ONLY: read_wam_f107_kp_txt, 
     &                              f107_wy, kp_wy, f107_kp_size,
     &            fix_spweather_data, kpa_wy, nhp_wy, nhpi_wy, f107d_wy,
     &                              shp_wy, shpi_wy,
     &                    swbt_wy, swvel_wy, swang_wy, swbz_wy, swden_wy
!

!
      implicit none
!
!SK   include 'mpif.h'
! define some constants
      integer, intent(in)  :: mpi_id         !=me  <-- my MPI-rank
      integer, intent(in)  :: master         !<-- master MPI-rank
      integer, intent(in)  :: levs           !number of pressure level
! Local variables  !SK
      logical :: same_MPI_rank
!
      character(len=256)   :: noeof_file
!
      integer, parameter        :: nz_euv = 17
      real, dimension(nz_euv)   :: effeuv17,effuv17, p17, z17
!
      real                :: dz, dz1 
      real                :: z(levs)
      integer             :: kup  ,kdw  
      integer             :: i, k,  j , kref, jinv
      real                :: unit_conv = 1.e-18        ! Xsections into cm^2 

! O2 scale factor for 15 pressure level 
      integer, parameter  ::  nz_Jo2sf = 15
      real                ::  JJ_scale_factor(nz_Jo2sf)
      real                ::  z15(nz_Jo2sf)
! TIEGCM SRB-efficiency factors on 63-pressures
      integer, parameter  ::  nz_63 = 63
      real                ::  pres63(nz_63), SRBEFF63(nz_63)
      real                ::  z63(nz_63)
!
! solar calendar and data files
      character(ch100)    :: nc_file, ncfile_fpath
      integer             :: Mjdat(ndwam)         ! IDAT_WAM +IrHOUR_WAM
      real                :: Hcur                 ! IDAT_WAM = Mjdatfor init-stage
!
!  wavelength (cm)
!
!     TIEGCM data
!
      real                    :: dsfmin(nwaves), dafac(nwaves)
      real                    :: drlmeuv(nwaves) ! wavelengths (cm)
!                                           ! absoption cross sections (x1e18cm^2)
      real, dimension(nwaves) :: sigeuv_o,sigeuv_o2,sigeuv_N2, sigeuv_o3
!                                           !ionization branching ratios (off absorption)
      real, dimension(nwaves) ::  BphotonI_O, BphotonI_O2,BphotonI_N2
!                                           ! O2 dsociation branching ratios
      real, dimension(nwaves) ::  bro2DPh(nwaves)
!                                           !O2 Photo-e dsociation branching ratios
      real, dimension(nwaves) ::  bro2Del(nwaves)

      DATA  drlmeuv/1.725e-05, 1.675e-05, 1.625e-05, 1.575e-05,
     &             1.525e-05, 1.475e-05, 1.425e-05, 1.375e-05,
     &             1.325e-05, 1.275e-05, 1.225e-05, 1.216e-05,
     &             1.175e-05, 1.125e-05, 1.075e-05, 1.038e-05,
     &             1.007e-05, 9.810e-06, 9.440e-06, 9.440e-06,
     &             9.440e-06, 8.555e-06, 8.555e-06, 8.555e-06,
     &             7.240e-06, 7.240e-06, 5.950e-06, 4.300e-06,
     &             3.050e-06, 2.570e-06, 1.895e-06, 1.125e-06,
     &             5.100e-07, 2.500e-07, 1.300e-07, 6.000e-08,
     &             2.250e-08/

! O2 absorption coefficient: (1,2,3,4)= (O2,O, N2,O3) 
! Cross Section cm^2 after scaling by 10-18 (unit_conv)
      DATA  sigeuv_O2 /                                                 
     &           5.00e-01, 1.50e+00, 3.40e+00, 6.00e+00, 1.00e+01,
     &           1.30e+01, 1.50e+01, 1.20e+01, 2.20e+00, 4.00e-01,
     &           1.30e+01, 1.00e-02, 1.40e+00, 4.00e-01, 1.00e+00,
     &           1.15e+00, 1.63e+00, 1.87e+01, 3.25e+01, 1.44e+01,
     &           1.34e+01, 1.33e+01, 1.09e+01, 1.05e+01, 2.49e+01,
     &           2.36e+01, 2.70e+01, 2.03e+01, 1.68e+01, 1.32e+01,
     &           7.63e+00, 2.63e+00, 6.46e-01, 2.10e-01, 2.25e-01,
     &           3.40e-02, 4.54e-03/

! O absorption coefficient:
      DATA  sigeuv_O /                                                  
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 3.79e+00, 4.10e+00, 3.00e+00, 4.79e+00,
     &           8.52e+00, 1.31e+01, 1.07e+01, 7.72e+00, 6.02e+00,
     &           3.78e+00, 1.32e+00, 3.25e-01, 1.05e-01, 1.13e-01,
     &           1.70e-02, 2.27e-03/

! N2 absorption coefficient:
      DATA  sigeuv_N2 /                                                 
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 2.55e+00, 1.15e+02, 1.44e+01,
     &           2.18e+00, 7.17e+01, 1.31e+01, 2.14e+00, 5.45e+01,
     &           2.30e+01, 2.31e+01, 1.97e+01, 1.17e+01, 9.94e+00,
     &           5.09e+00, 1.53e+00, 3.46e-01, 1.14e+00, 1.41e-01,
     &           2.01e-02, 2.53e-03/
!
! O3 absorption cross section
!
      sigeuv_o3 = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 1.25e+01, 9.20e+00,
     |          9.20e+00, 9.20e+00, 9.16e+00, 9.50e+00, 9.50e+00,
     |          9.50e+00, 1.47e+01, 1.47e+01, 1.47e+01, 2.74e+01,
     |          2.02e+01, 3.33e+01, 7.74e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00/)

 
! The three major species' ionization branching ratio (off absorption):
! O2
      DATA  BPhotonI_O2 /                                               
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 6.13e-01, 8.30e-01, 6.20e-01, 7.86e-01,
     &           7.56e-01, 5.34e-01, 5.74e-01, 5.49e-01, 4.76e-01,
     &           6.73e-01, 9.83e-01, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 1.00e+00/
! O       use idea_solar_input,  only : itheia

      DATA  BPhotonI_O /                                                
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 1.00e+00/
! N2
      DATA  BPhotonI_N2 /                                               
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 4.29e-01,
     &           6.80e-01, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 1.00e+00/
! O2 photon dsociation branching ratio
      DATA  bro2DPh  /                                                  
     &           1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 3.87e-01, 1.70e-01, 3.80e-01, 2.14e-01,
     &           2.44e-01, 4.66e-01, 4.26e-01, 4.51e-01, 5.24e-01,
     &           3.27e-01, 1.74e-02, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00/
! O2 photoelectron dsociation branching ratio
      DATA  bro2DEl  /                                                  
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 1.10e-02, 6.53e-01, 7.62e-01, 9.96e-01,
     &           1.27e+00, 2.04e+00, 4.11e+00, 5.70e+01, 1.78e+01,
     &           2.03e+01, 8.79e+01/
! Solar spectrum based on EUVAC and glow for wave length less than 1050 A
! and Woods for wavelength greater than 1050 A
! solar minimum flux (when P_index=80, unit:photon cm^-2 S^-1)
      DATA  dsfmin /3.397e+11, 1.998e+11, 1.055e+11, 7.260e+10,
     &            5.080e+10, 2.802e+10, 1.824e+10, 1.387e+10,
     &            2.659e+10, 7.790e+09, 1.509e+10, 3.940e+11,
     &            8.399e+09, 3.200e+09, 3.298e+09, 4.235e+09,
     &            4.419e+09, 4.482e+09, 7.156e+08, 1.028e+09,
     &            3.818e+08, 8.448e+08, 3.655e+09, 2.364e+09,
     &            1.142e+09, 1.459e+09, 4.830e+09, 2.861e+09,
     &            8.380e+09, 4.342e+09, 5.612e+09, 1.270e+09,
     &            5.326e+08, 2.850e+07, 2.000e+06, 1.000e+04,
     &            5.010e+01/
! scaling factor A as defined in EUVAC model
      DATA  dafac /5.937e-04, 6.089e-04, 1.043e-03, 1.125e-03,
     &           1.531e-03, 1.202e-03, 1.873e-03, 2.632e-03,
     &           2.877e-03, 2.610e-03, 3.739e-03, 4.230e-03,
     &           2.541e-03, 2.099e-03, 3.007e-03, 4.825e-03,
     &           5.021e-03, 3.950e-03, 4.422e-03, 4.955e-03,
     &           4.915e-03, 5.437e-03, 5.261e-03, 5.310e-03,
     &           3.680e-03, 5.719e-03, 5.857e-03, 1.458e-02,
     &           7.059e-03, 2.575e-02, 1.433e-02, 9.182e-03,
     &           1.343e-02, 6.247e-02, 2.000e-01, 3.710e-01,
     &           6.240e-01/
!============================================================
!
       DATA JJ_scale_factor/ 1.0, 4.465536, 4.365480, 3.904985,         
     &                    3.367959, 3.202786, 2.378429, 1.636311,       
     &                    1.423021, 1.452178, 1.588099, 1.714328,       
     &                    1.811639, 1.907779, 1.987971/

! Log(pres) of efficiencies grid
       DATA pres63/6.90775528,  6.57442194,                             
     &   6.24108859,  5.90775525,  5.57442191,                          
     &   5.24108856,  4.90775522,  4.57442188,                          
     &   4.24108853,  3.90775519,  3.57442185,                          
     &   3.2410885,   2.90775516,  2.57442182,                          
     &   2.24108847,  1.90775513,  1.57442179,                          
     &   1.24108844,  0.9077551,   0.574421757,                         
     &   0.241088414, -0.0922449296,-0.425578273,                       
     &   -0.758911616,-1.09224496,  -1.4255783,                         
     &   -1.75891165, -2.09224499,  -2.42557833,
     &   -2.75891168, -3.09224502,  -3.42557836,
     &   -3.75891171, -4.09224505,  -4.42557839,
     &   -4.75891174, -5.09224508,  -5.42557842,
     &   -5.75891177, -6.09224511,  -6.42557845,
     &   -6.75891179, -7.09224514,  -7.42557848,
     &   -7.75891182, -8.09224517,  -8.42557851,
     &   -8.75891185, -9.0922452,   -9.42557854,
     &   -9.75891188, -10.0922452,  -10.4255786,
     &   -10.7589119, -11.0922453,  -11.4255786,
     &   -11.7589119, -12.0922453,  -12.4255786,
     &   -12.758912,  -13.0922453,  -13.4255787,
     &   -13.758912/

! SRB heating by O2 (don't need efficiencies if self constently
! calculating chemical heating.
! Mylynczak and Solomon bulk O2 heating efficienies
!
        DATA SRBEFF63/1.000,1.000,1.000,1.000,1.000,                    
     &   1.000,1.000,1.000,1.000,1.000,1.000,                           
     &   1.000,1.000,1.000,.980,.983,.982,                              
     &   .970,.945,.913,.880,.852,.832,.820,                            
     &   .817,.816,.810,.795,.777,.765,.764,
     &   .759,.730,.664,.579,.500,.446,.416,
     &   .400,.393,.390,.390,.391,.391,.390,
     &   .388,.384,.380,.375,.366,.350,.324,
     &   .291,.260,.234,.214,.200,.190,.184,
     &   .180,.176,.173,.170/    
!
!  ** EUV and UV heating efficiency on 17 pressure levels
!
      DATA EFFEUV17/8*1.0,.75,.6,.62,.54,.49,.41,.33,.30,.30/
!
! OLD     DATA EFFUV17/5*.28,.29,.32,.38,.4,.4,.4,.39,.34,.26,.19,.17,.16/
!
! Tab-EFF       P(pa)        EUV          UV           JO2
!       1      5.22850      1.00000     0.280000      1.00000
!       2      1.92346      1.00000     0.280000      1.00000
!       3     0.707601      1.00000     0.280000      1.00000
!       4     0.260312      1.00000     0.280000      4.42724
!       5    0.0957633      1.00000     0.280000      4.18921
!       6    0.0352294      1.00000     0.290000      3.69942
!       7    0.0129602      1.00000     0.320000      3.30473
!       8   0.00476777      1.00000     0.380000      2.88723
!       9   0.00175397     0.750000     0.400000      2.09436
!      10  0.000645248     0.600000     0.400000      1.55467
!      11  0.000237374     0.620000     0.400000      1.43418
!      12  8.73248e-05     0.540000     0.390000      1.50421
!      13  3.21250e-05     0.490000     0.340000      1.63642
!      14  1.18181e-05     0.410000     0.260000      1.75158
!      15  4.34765e-06     0.330000     0.190000      1.84844
!      16  1.59941e-06     0.300000     0.170000      1.92998
!      17  5.88390e-07     0.300000     0.160000      1.98797
!
! NEW  8-points
!
       DATA EFFUV17/0.59, 0.59, 0.58, 0.57, 0.56, 0.52, 0.48, 0.43,
     &.4,.4,.4,.39,.34,.26,.19,.17,.16/
 
        same_MPI_rank = mpi_id == master   !SK2020
        afac  = dafac
        sfmin = dsfmin
        rlmeuv =drlmeuv
!======================================================================
!
! VAY 12/2026: start with selection of SPW_DRIVERS and read of 
!       "GLOBAL" year-to-year & day-to-day variable solar-geo inputs
!
! Below we read data on each PE separately, namelts/global F107/Kp/NO 
!
! should be done on the designated IO-PE and MPI_BCASTED => to all PEs
!
! th was done on  Zeus, but initial issues with mpif.h on Theia
! enforced to read on all PEs.....
!
!======================================================================

!SK   call MPI_COMM_RANK (MPI_COMM_WORLD, mpi_id, mpi_err)
      CALL solar_read_namelist(nml_solar, nlun_solar, ncfile_fpath, 
     &     noeof_file,  mpi_id)
!
!SNOE-NO
!
      noeof_file='snoe_eof.nc'
!
      call solar_readno_snoewx(noeof_file, mpi_id)
!
! F107/Kp
!
! TODO Reserve 3-paths Swpc_fst, Sair_2012, WX_YYYYs, 
!                      now only WACCM-X and SWPC-3day forecast files
!                      many thanks to Dan Marsh, ACOM/NCAR
!  '
!
       CALL  COPY_IDAT_NEMS_2_WAM      !------------------------------------------------------------------        
                                       ! copy NEMS-cal[Y4-D3-M2-H1] to idat_wam (Y1-M2-D3-H4), irhour_wam
                                       ! idate ...... in sighdr & sfchdr
                                       ! -----------------------------------------------------------------
!
       CALL  CURRENT_NCEP_JDAT(idat_wam, irhour_wam, Mjdat, Hcur)
  
       if (mpi_id ==0 ) print *, 'idea_solar_init, idat_wam ',  idat_wam
!      print *, 'idea_solar_init, idat_wam ', idat_wam, ' me=', mpi_id
!      if (same_MPI_rank) print *,
!    &      'idea_solar_init, idat_wam ',  idat_wam
!
!     data_swpc make a final decion keep it in "gloopb.f" or in "solar_init.f"
!



        if (trim(SPW_DRIVERS)=='swpc_fst') then 

          IF(.NOT. ALLOCATED(f107_wy))  ALLOCATE(f107_wy (f107_kp_size))
          IF(.NOT. ALLOCATED(kp_wy))    ALLOCATE(kp_wy   (f107_kp_size))
          IF(.NOT. ALLOCATED(f107d_wy)) ALLOCATE(f107d_wy(f107_kp_size))
          IF(.NOT. ALLOCATED(kpa_wy))   ALLOCATE(kpa_wy  (f107_kp_size))
          IF(.NOT. ALLOCATED(nhp_wy))   ALLOCATE(nhp_wy  (f107_kp_size))
          IF(.NOT. ALLOCATED(nhpi_wy))  ALLOCATE(nhpi_wy (f107_kp_size))
          IF(.NOT. ALLOCATED(shp_wy))   ALLOCATE(shp_wy  (f107_kp_size))
          IF(.NOT. ALLOCATED(shpi_wy))  ALLOCATE(shpi_wy (f107_kp_size))
          IF(.NOT. ALLOCATED(swbt_wy))  ALLOCATE(swbt_wy (f107_kp_size))
          IF(.NOT. ALLOCATED(swang_wy)) ALLOCATE(swang_wy(f107_kp_size))
          IF(.NOT. ALLOCATED(swvel_wy)) ALLOCATE(swvel_wy(f107_kp_size))
          IF(.NOT. ALLOCATED(swbz_wy))  ALLOCATE(swbz_wy (f107_kp_size))
          IF(.NOT. ALLOCATED(swden_wy)) ALLOCATE(swden_wy(f107_kp_size))
          call read_wam_f107_kp_txt
!         if (same_MPI_rank) write(6,*) 
          if (mpi_id ==0 ) write(6,*) 
     & ' SPW_DRIVERS => swpc_fst, 3-day forecasts:', trim(SPW_DRIVERS)
        endif
!
!  data only for 2012
!
      if (trim(SPW_DRIVERS)=='sair_wam'.and.idea_solar_fix.le.1 ) then
           CALL solar_read_wam_init(ncfile_fpath, mpi_id) 
           CALL solar_wamstep_advance(mpi_id, Mjdat, Hcur)
       endif
!
! data_wx
!
       if (trim(SPW_DRIVERS)=='cires_wam'.and.idea_solar_fix.le.1) then
           CALL solar_read_myy1947_2016(ncfile_fpath, mpi_id )
           CALL solar_waccmx_advance(mpi_id, Mjdat, Hcur, 1)
       endif
!
!climate_fix TODO  from solar_in
!
       if (trim(SPW_DRIVERS)=='wam_climate') then
          wap_s = 3.0
          wkp_s = 1.0
          wf107_s =100.
          wf107a_s =100.
       endif

!     if (same_MPI_rank) then
      if (mpi_id == 0) then
        write(6,*) ' VAY-end of idea_solar_init in idea_solar_heating.f'
      endif
!
      allocate (effeuv(levs), effuv(levs))
! efficiencies suggested by Mlynchack and Solomon for UV-heating rates 40-110 km
!
      allocate (eff_hart(levs), eff_hugg(levs), eff_chap(levs), 
     &  eff_herz(levs))
      allocate (eff_srb(levs), eff_src(levs), eff_lya(levs))

      eff_srb(1:levs) =1.00
      eff_lya(1:levs) =0.95
      eff_src(1:levs) =0.85   ! should be modified to introduce vertical profile with min at ~90km
      eff_hugg(1:levs)=1.00
      eff_chap(1:levs)=1.00
      eff_herz(1:levs)=1.00
      eff_hart(1:levs)=1.00   ! should be modified to suppress ozone mesospheric heating
!
      allocate (o2_scale_factor(levs),  SRBEFF(levs) )  

!
!SK 2020Jun30
      if (same_MPI_rank) print *, 'Returned from idea_solar_init'
      RETURN
      end subroutine idea_solar_init

      subroutine idea_solar_pre(mpi_id, master, levs)
!----------------------------------------------------------------------------
! define: effuv,effeuv, solar fluxes (37), wlengths, cross-sections,
!     other eff_hart, eff_hugg, eff_chap, eff_herz, eff_srb, eff_src, eff_lya
!     init solar_calendar and F107/Kp datasets
!
!----------------------------------------------------------------------------    
      use machine,            only : kind_phys
      use idea_wam_control,   only : SPW_DRIVERS
      use IDEA_IO_UNITS    ,  only : nml_solar, nlun_solar, ch100
!
      use idea_solarno_input, only : solar_readno_snoewx

      use idea_solar_input,  only : idea_solar_fix, itheia
      use idea_solar_input,  only : wf107_s,   wkp_s
      use idea_solar_input,  only : wf107a_s,  wap_s
      use idea_solar_input,  only : weuv_s,  nwafix

      use idea_solar_input,  only : solar_read_namelist
      use idea_solar_input,  only : solar_read_myy1947_2016, 
     &                              solar_read_wam_init
      use idea_solar_input,  only : solar_waccmx_advance
      use idea_solar_input,  only : solar_wamstep_advance

!      use idea_solar_input,  only : f107 => wf107_s,   kp => wkp_s
!      use idea_solar_input,  only : f107a => wf107a_s, Ap => wap_s
!      use idea_solar_input,  only : EUV37 => weuv_s,  nsp_euv => nwafix
!
!SK   use idea_mpi_def,       only : mpi_id, mpi_err, MPI_COMM_ALL,info
!
      use idea_solar,         only : nps,  effuv, effeuv, pr, prdot02
      use idea_solar,         only : o2_scale_factor,  SRBEFF
!
      use idea_solar,      only : eff_hart, eff_hugg, eff_chap, eff_herz
      use idea_solar,      only : eff_srb, eff_src, eff_lya
!
     
       use idea_solar,       only  : rwpcc, pcc
       use idea_solar,       only  : csao, csao2, csao3, csan2
       use idea_solar,       only  : csio, csio2, csin2
       use idea_solar,       only  : csdo2, csdeo2
       use idea_solar,       only  : eccentric
       use idea_solar,       only  : nwaves,NWAVES_EUV,lyman_a_num,
     &       nwaves_src
       use idea_solar,       only  : sfmin, afac,rlmeuv
! 
      use wam_date_calendar, only : COPY_IDAT_NEMS_2_WAM, idat_wam,
     &            irhour_wam
      use wam_date_calendar, only : CURRENT_NCEP_JDAT, ndwam  
      use wam_date_calendar, only : curday_wam, curmonth_wam, curddd_wam 
      use wam_date_calendar, only : curyear_wam, curutsec_wam 
!
      use wam_f107_kp_mod,    ONLY: read_wam_f107_kp_txt, 
     &                              f107_wy, kp_wy, f107_kp_size,
     &            fix_spweather_data, kpa_wy, nhp_wy, nhpi_wy, f107d_wy,
     &                              shp_wy, shpi_wy,
     &                    swbt_wy, swvel_wy, swang_wy, swbz_wy, swden_wy
!

!
      implicit none
!
!SK   include 'mpif.h'
! define some constants
      integer, intent(in)  :: mpi_id         !=me  <-- my MPI-rank
      integer, intent(in)  :: master         !<-- master MPI-rank
      integer, intent(in)  :: levs           !number of pressure level
! Local variables  !SK
      logical :: same_MPI_rank
!
      character(len=256)   :: noeof_file
!
      integer, parameter        :: nz_euv = 17
      real, dimension(nz_euv)   :: effeuv17,effuv17, p17, z17
!
      real                :: dz, dz1 
      real                :: z(levs)
      integer             :: kup  ,kdw  
      integer             :: i, k,  j , kref, jinv
      real                :: unit_conv = 1.e-18        ! Xsections into cm^2 

! O2 scale factor for 15 pressure level 
      integer, parameter  ::  nz_Jo2sf = 15
      real                ::  JJ_scale_factor(nz_Jo2sf)
      real                ::  z15(nz_Jo2sf)
! TIEGCM SRB-efficiency factors on 63-pressures
      integer, parameter  ::  nz_63 = 63
      real                ::  pres63(nz_63), SRBEFF63(nz_63)
      real                ::  z63(nz_63)
!
! solar calendar and data files
      character(ch100)    :: nc_file, ncfile_fpath
      integer             :: Mjdat(ndwam)         ! IDAT_WAM +IrHOUR_WAM
      real                :: Hcur                 ! IDAT_WAM = Mjdatfor init-stage
!
!  wavelength (cm)
!
!     TIEGCM data
!
      real                    :: dsfmin(nwaves), dafac(nwaves)
      real                    :: drlmeuv(nwaves) ! wavelengths (cm)
!                                           ! absoption cross sections (x1e18cm^2)
      real, dimension(nwaves) :: sigeuv_o,sigeuv_o2,sigeuv_N2, sigeuv_o3
!                                           !ionization branching ratios (off absorption)
      real, dimension(nwaves) ::  BphotonI_O, BphotonI_O2,BphotonI_N2
!                                           ! O2 dsociation branching ratios
      real, dimension(nwaves) ::  bro2DPh(nwaves)
!                                           !O2 Photo-e dsociation branching ratios
      real, dimension(nwaves) ::  bro2Del(nwaves)

      DATA  drlmeuv/1.725e-05, 1.675e-05, 1.625e-05, 1.575e-05,
     &             1.525e-05, 1.475e-05, 1.425e-05, 1.375e-05,
     &             1.325e-05, 1.275e-05, 1.225e-05, 1.216e-05,
     &             1.175e-05, 1.125e-05, 1.075e-05, 1.038e-05,
     &             1.007e-05, 9.810e-06, 9.440e-06, 9.440e-06,
     &             9.440e-06, 8.555e-06, 8.555e-06, 8.555e-06,
     &             7.240e-06, 7.240e-06, 5.950e-06, 4.300e-06,
     &             3.050e-06, 2.570e-06, 1.895e-06, 1.125e-06,
     &             5.100e-07, 2.500e-07, 1.300e-07, 6.000e-08,
     &             2.250e-08/

! O2 absorption coefficient: (1,2,3,4)= (O2,O, N2,O3) 
! Cross Section cm^2 after scaling by 10-18 (unit_conv)
      DATA  sigeuv_O2 /                                                 
     &           5.00e-01, 1.50e+00, 3.40e+00, 6.00e+00, 1.00e+01,
     &           1.30e+01, 1.50e+01, 1.20e+01, 2.20e+00, 4.00e-01,
     &           1.30e+01, 1.00e-02, 1.40e+00, 4.00e-01, 1.00e+00,
     &           1.15e+00, 1.63e+00, 1.87e+01, 3.25e+01, 1.44e+01,
     &           1.34e+01, 1.33e+01, 1.09e+01, 1.05e+01, 2.49e+01,
     &           2.36e+01, 2.70e+01, 2.03e+01, 1.68e+01, 1.32e+01,
     &           7.63e+00, 2.63e+00, 6.46e-01, 2.10e-01, 2.25e-01,
     &           3.40e-02, 4.54e-03/

! O absorption coefficient:
      DATA  sigeuv_O /                                                  
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 3.79e+00, 4.10e+00, 3.00e+00, 4.79e+00,
     &           8.52e+00, 1.31e+01, 1.07e+01, 7.72e+00, 6.02e+00,
     &           3.78e+00, 1.32e+00, 3.25e-01, 1.05e-01, 1.13e-01,
     &           1.70e-02, 2.27e-03/

! N2 absorption coefficient:
      DATA  sigeuv_N2 /                                                 
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 2.55e+00, 1.15e+02, 1.44e+01,
     &           2.18e+00, 7.17e+01, 1.31e+01, 2.14e+00, 5.45e+01,
     &           2.30e+01, 2.31e+01, 1.97e+01, 1.17e+01, 9.94e+00,
     &           5.09e+00, 1.53e+00, 3.46e-01, 1.14e+00, 1.41e-01,
     &           2.01e-02, 2.53e-03/
!
! O3 absorption cross section
!
      sigeuv_o3 = (/
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 1.25e+01, 9.20e+00,
     |          9.20e+00, 9.20e+00, 9.16e+00, 9.50e+00, 9.50e+00,
     |          9.50e+00, 1.47e+01, 1.47e+01, 1.47e+01, 2.74e+01,
     |          2.02e+01, 3.33e+01, 7.74e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     |          0.00e+00, 0.00e+00/)

 
! The three major species' ionization branching ratio (off absorption):
! O2
      DATA  BPhotonI_O2 /                                               
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 6.13e-01, 8.30e-01, 6.20e-01, 7.86e-01,
     &           7.56e-01, 5.34e-01, 5.74e-01, 5.49e-01, 4.76e-01,
     &           6.73e-01, 9.83e-01, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 1.00e+00/
! O       use idea_solar_input,  only : itheia

      DATA  BPhotonI_O /                                                
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 1.00e+00/
! N2
      DATA  BPhotonI_N2 /                                               
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 4.29e-01,
     &           6.80e-01, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 1.00e+00/
! O2 photon dsociation branching ratio
      DATA  bro2DPh  /                                                  
     &           1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00,
     &           1.00e+00, 3.87e-01, 1.70e-01, 3.80e-01, 2.14e-01,
     &           2.44e-01, 4.66e-01, 4.26e-01, 4.51e-01, 5.24e-01,
     &           3.27e-01, 1.74e-02, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00/
! O2 photoelectron dsociation branching ratio
      DATA  bro2DEl  /                                                  
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     &           0.00e+00, 1.10e-02, 6.53e-01, 7.62e-01, 9.96e-01,
     &           1.27e+00, 2.04e+00, 4.11e+00, 5.70e+01, 1.78e+01,
     &           2.03e+01, 8.79e+01/
! Solar spectrum based on EUVAC and glow for wave length less than 1050 A
! and Woods for wavelength greater than 1050 A
! solar minimum flux (when P_index=80, unit:photon cm^-2 S^-1)
      DATA  dsfmin /3.397e+11, 1.998e+11, 1.055e+11, 7.260e+10,
     &            5.080e+10, 2.802e+10, 1.824e+10, 1.387e+10,
     &            2.659e+10, 7.790e+09, 1.509e+10, 3.940e+11,
     &            8.399e+09, 3.200e+09, 3.298e+09, 4.235e+09,
     &            4.419e+09, 4.482e+09, 7.156e+08, 1.028e+09,
     &            3.818e+08, 8.448e+08, 3.655e+09, 2.364e+09,
     &            1.142e+09, 1.459e+09, 4.830e+09, 2.861e+09,
     &            8.380e+09, 4.342e+09, 5.612e+09, 1.270e+09,
     &            5.326e+08, 2.850e+07, 2.000e+06, 1.000e+04,
     &            5.010e+01/
! scaling factor A as defined in EUVAC model
      DATA  dafac /5.937e-04, 6.089e-04, 1.043e-03, 1.125e-03,
     &           1.531e-03, 1.202e-03, 1.873e-03, 2.632e-03,
     &           2.877e-03, 2.610e-03, 3.739e-03, 4.230e-03,
     &           2.541e-03, 2.099e-03, 3.007e-03, 4.825e-03,
     &           5.021e-03, 3.950e-03, 4.422e-03, 4.955e-03,
     &           4.915e-03, 5.437e-03, 5.261e-03, 5.310e-03,
     &           3.680e-03, 5.719e-03, 5.857e-03, 1.458e-02,
     &           7.059e-03, 2.575e-02, 1.433e-02, 9.182e-03,
     &           1.343e-02, 6.247e-02, 2.000e-01, 3.710e-01,
     &           6.240e-01/
!============================================================
!
       DATA JJ_scale_factor/ 1.0, 4.465536, 4.365480, 3.904985,         
     &                    3.367959, 3.202786, 2.378429, 1.636311,       
     &                    1.423021, 1.452178, 1.588099, 1.714328,       
     &                    1.811639, 1.907779, 1.987971/

! Log(pres) of efficiencies grid
       DATA pres63/6.90775528,  6.57442194,                             
     &   6.24108859,  5.90775525,  5.57442191,                          
     &   5.24108856,  4.90775522,  4.57442188,                          
     &   4.24108853,  3.90775519,  3.57442185,                          
     &   3.2410885,   2.90775516,  2.57442182,                          
     &   2.24108847,  1.90775513,  1.57442179,                          
     &   1.24108844,  0.9077551,   0.574421757,                         
     &   0.241088414, -0.0922449296,-0.425578273,                       
     &   -0.758911616,-1.09224496,  -1.4255783,                         
     &   -1.75891165, -2.09224499,  -2.42557833,
     &   -2.75891168, -3.09224502,  -3.42557836,
     &   -3.75891171, -4.09224505,  -4.42557839,
     &   -4.75891174, -5.09224508,  -5.42557842,
     &   -5.75891177, -6.09224511,  -6.42557845,
     &   -6.75891179, -7.09224514,  -7.42557848,
     &   -7.75891182, -8.09224517,  -8.42557851,
     &   -8.75891185, -9.0922452,   -9.42557854,
     &   -9.75891188, -10.0922452,  -10.4255786,
     &   -10.7589119, -11.0922453,  -11.4255786,
     &   -11.7589119, -12.0922453,  -12.4255786,
     &   -12.758912,  -13.0922453,  -13.4255787,
     &   -13.758912/

! SRB heating by O2 (don't need efficiencies if self constently
! calculating chemical heating.
! Mylynczak and Solomon bulk O2 heating efficienies
!
        DATA SRBEFF63/1.000,1.000,1.000,1.000,1.000,                    
     &   1.000,1.000,1.000,1.000,1.000,1.000,                           
     &   1.000,1.000,1.000,.980,.983,.982,                              
     &   .970,.945,.913,.880,.852,.832,.820,                            
     &   .817,.816,.810,.795,.777,.765,.764,
     &   .759,.730,.664,.579,.500,.446,.416,
     &   .400,.393,.390,.390,.391,.391,.390,
     &   .388,.384,.380,.375,.366,.350,.324,
     &   .291,.260,.234,.214,.200,.190,.184,
     &   .180,.176,.173,.170/    
!
!  ** EUV and UV heating efficiency on 17 pressure levels
!
      DATA EFFEUV17/8*1.0,.75,.6,.62,.54,.49,.41,.33,.30,.30/
!
! OLD     DATA EFFUV17/5*.28,.29,.32,.38,.4,.4,.4,.39,.34,.26,.19,.17,.16/
!
! Tab-EFF       P(pa)        EUV          UV           JO2
!       1      5.22850      1.00000     0.280000      1.00000
!       2      1.92346      1.00000     0.280000      1.00000
!       3     0.707601      1.00000     0.280000      1.00000
!       4     0.260312      1.00000     0.280000      4.42724
!       5    0.0957633      1.00000     0.280000      4.18921
!       6    0.0352294      1.00000     0.290000      3.69942
!       7    0.0129602      1.00000     0.320000      3.30473
!       8   0.00476777      1.00000     0.380000      2.88723
!       9   0.00175397     0.750000     0.400000      2.09436
!      10  0.000645248     0.600000     0.400000      1.55467
!      11  0.000237374     0.620000     0.400000      1.43418
!      12  8.73248e-05     0.540000     0.390000      1.50421
!      13  3.21250e-05     0.490000     0.340000      1.63642
!      14  1.18181e-05     0.410000     0.260000      1.75158
!      15  4.34765e-06     0.330000     0.190000      1.84844
!      16  1.59941e-06     0.300000     0.170000      1.92998
!      17  5.88390e-07     0.300000     0.160000      1.98797
!
! NEW  8-points
!
       DATA EFFUV17/0.59, 0.59, 0.58, 0.57, 0.56, 0.52, 0.48, 0.43,
     &.4,.4,.4,.39,.34,.26,.19,.17,.16/
 
!=====================================================================================
!
!
! Below initilization of IDEA_SOLAR_HEAT, cross sections/efficiencies/NO-snoe
!
!     also a good place to "init constants" for Ozone/O2  MA heating rates
!
!=====================================================================================
      same_MPI_rank = mpi_id == master   !SK2020
!     find nps
      do k=1,levs
        if(pr(k).le.prdot02) then
          nps=k
          exit
        endif
      enddo
      if (same_MPI_rank) print *, 
     &   'idea_solar_init, VAY solar-init nps', nps
!SK    print *, 'idea_solar_init, VAY solar-init nps', nps
!
! get effuv,effeuv from interplating effuv17, effeuv17 to 150 levs 
! get no from interplating no17 to 150 levs  
      do k=1, nz_euv
        p17(k)=5.2285*exp(1.-k)
        z17(k)=-log(p17(k))
      enddo
      do k=1, nz_Jo2sf
        z15(k)=-log(1.0376*exp(1.-k))
      enddo
      do k=1, nz_63
        z63(k) = -pres63(k)
      enddo
      do k=1,levs
        z(k)= -log(pr(k)*100.)
      enddo
!
!SK-> allocate (effeuv(levs), effuv(levs))
! efficiencies suggested by Mlynchack and Solomon for UV-heating rates 40-110 km
!
!     allocate (eff_hart(levs), eff_hugg(levs), eff_chap(levs), 
!    &  eff_herz(levs))
!     allocate (eff_srb(levs), eff_src(levs), eff_lya(levs))

!     eff_srb(1:levs) =1.00
!     eff_lya(1:levs) =0.95
!     eff_src(1:levs) =0.85   ! should be modified to introduce vertical profile with min at ~90km
!     eff_hugg(1:levs)=1.00
!     eff_chap(1:levs)=1.00
!     eff_herz(1:levs)=1.00
!     eff_hart(1:levs)=1.00   ! should be modified to suppress ozone mesospheric heating
!<-SK
      call heat_uveff(levs, eff_hart, eff_src)
!========================================================================
!
!SK   allocate (o2_scale_factor(levs),  SRBEFF(levs) )  
!
! (17) => (levs)
!
!SK2020Sep8 the 8th-argument kdw is bogus, hence removed.
!SK   CALL interpol_wamz(nz_euv, z17, effuv17,  levs, Z, effuv, kup,kdw)
!SK   CALL interpol_wamz(nz_euv, z17, effeuv17, levs, Z, effeuv,kup,kdw)
      CALL interpol_wamz(nz_euv, z17, effuv17,  levs, Z, effuv, kup)
      CALL interpol_wamz(nz_euv, z17, effeuv17, levs, Z, effeuv,kup)
      CALL interpol_wamz_down(nz_63, z63, SRBEFF63, levs, Z,
     & SRBEFF, 1.0) 
      CALL  interpol_wamz_down(nz_Jo2sf, z15,JJ_scale_factor , levs, Z,
     & o2_scale_factor , 1.0) 
!===============================================================VAY-2016
!
!      before taking out J-eff & srb-eff put them as 1.0
!
!       o2_scale_factor(1:levs) =1.0
!       SRBEFF(1:levs) =1.0
!
!===============================================================VAY-2016

!
!
!  faulty/dirty vert. interpolations are replaced by "interpol_wamz"-
!      should be applied for all WAM-physics with external inputs
!      call z15toz(levs,JJ_scale_factor, Z, o2_scale_factor,0.)
!      call z63toz(levs,     SRBEFF63,   Z,   SRBEFF,       0.)
!      do k=1, levs
!       print *, 'VAY-EF',Z(k)*7.-35., SRBEFF(k), o2_scale_factor(k)
!      enddo
!
!      print *, 'VAY-SRBEF', maxval(SRBEFF63), minval(SRBEFF63)
!      print *, 'VAY-SRBEF', maxval(SRBEFF), minval(SRBEFF)
! compute:  RWPCC
! VAY
        do J = 1, NWAVES  
        jinv = nwaves+1-J
!         WAVELS(J) = rlmeuv(jinv)*1.e-2              ! (convert to m)
!        if(wavels(J) == 1.216e-07) lyman_a_num  = J ! danger for real
         CSAO(J)   = sigeuv_O(jinv)*unit_conv
         CSAO2(J)  = sigeuv_O2(jinv)*unit_conv
         CSAO3(J)  = sigeuv_O3(jinv)*unit_conv
         CSAN2(J)  = sigeuv_N2(jinv)*unit_conv
         CSIO(J)   = CSAO(J)*BphotonI_O(jinv)
         CSIO2(J)  = CSAO2(J)*BphotonI_O2(jinv)
         CSIN2(J)  = CSAN2(J)*BphotonI_N2(jinv)
         CSDO2(J)  = CSAO2(J)*bro2Dph(jinv)
         CSDeO2(J) = CSIO2(J)*bro2Del(jinv)
         RWPCC(J)  = PCC/(rlmeuv(jinv)*1.e-2)      
         enddo
      RETURN
      end subroutine idea_solar_pre
