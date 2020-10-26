      module namelist_wamphysics_def
!
!VAY-20161101  WAM-control Namelist for Solar-Geo drivers
!spw_drivers  = 'cires_wam' ! 'swpc_fst', 'cires_wam', 'sair_wam', 'wam_climate'
      implicit none
!====================================================
! defines only wam-idea control variables
! that allow to switch WAM-space weather applications 
! 1) Space weather drivers
! 2) GW-physics GSM vs WAM
! 3) Solar_in
! 4) Ion_in
! 5) Das_in
! 6) Smet_in
! 7) Nc_inout
! 8) Tides_inout
! 9) Pws_inout
!======================================================
! Three Main types of WAM-predictions
! (1) Climate (2) SWPC-FKP-forecats (3) CIRES real data
!======================================================
      logical wam_climate
      logical wam_swpc_3day
      logical wam_cires_rdata
      logical wam_sair2012
!  in case wam_swpc_3day='true', and nowcast
      logical wam_swin
! wam-climate/perturb cases
      logical wam_smin
      logical wam_smax
      logical wam_saver
      logical wam_geostorm
!
! wam-physics
!
      logical wam_gwphys
      logical wam_solar_in
      logical wam_ion_in

      real JH0, JH_tanh, JH_semiann, JH_ann
      real skeddy0, skeddy_semiann, skeddy_ann
      real tkeddy0, tkeddy_semiann, tkeddy_ann
!
! DAS & Nudging on the fly
!
      logical wam_das_in
      logical wam_smet_in
!io-wam
      logical wam_netcdf_inout
!diagnostics
      logical wam_tides_diag
      logical wam_pws_diag
      logical wam_gws_diag
      public wam_control_default, select_case4_wamclimate
      contains
!
      subroutine wam_control_default
      character(len=100) :: spw_drivers
! 3 types of WAM predictions
      wam_climate     =.false.   ! simulations with arbitrary solar/wea forcing
      wam_cires_rdata =.true.    ! simulations with realistic daily F107/Kp YMD-WX-file
      wam_swpc_3day   =.false.   ! simulations with SWPC-predicted F107/Kp
      wam_sair2012    =.false.   ! simulations with realistic daily F107/Kp/EUV 2012-file
      spw_drivers  = 'cires_wam' ! 'swpc_fst', 'cires_wam', 'sair_wam', 'wam_climate'
      
      wam_swin   =.false.   ! simulations with kp f107, not solar wind
!
! wam_climate state:  Solar Min, Solar Max Solar Aver and perturbed case GEOSTORM
!
      wam_saver =.true.          ! defines  specific sets of F107 and Kp
      wam_smin  =.false.         !
      wam_smax  =.false.         !
      wam_geostorm  =.false.     !
!
! wam-physics
!
      wam_gwphys=.false.         ! GW-physics of GSM w/o vertical eddy dif/cond/visc
                                 ! .true.  all NGWs are in idea_phys with eddy effects
      wam_solar_in=.false.       ! no specific data files for EUV/F107/KP
                                 !  .true.  run like wam_cires_rdata with YMD-files
      wam_ion_in=.false.         ! no specific data files for aurora and ion-physics
                                 ! .true.   run with data files (imf, aurora etc...)
      JH0  = 1.75
      JH_tanh = 0.5
      JH_semiann = 0.5
      JH_ann  = 0.
      skeddy0 = 140.
      skeddy_semiann = 60.
      skeddy_ann     = 0.
      tkeddy0 = 280.
      tkeddy_semiann = 0.
      tkeddy_ann     = 0.

! WAM with "DAS & Nudging" on the fly the LA-drivers w/o GDAS/NOAA
!
      wam_das_in=.false.         ! UFO for SABER & MLS with EKF corrections
      wam_smet_in=.false.        ! Specified Meteorology Option like WX-GEOS5
!
! Additional WAM-netcdf output
      wam_netcdf_inout=.true.    ! NC-outs along with sigio/nemsio from DYN-core
! Extra diagnostics on the fly   ! in the NC-formats
      wam_tides_diag =.true.     ! Acos
      wam_pws_diag   =.true.
      wam_gws_diag   =.true.                  
!
      end subroutine wam_control_default
!
      subroutine select_case4_wamclimate
       write(6,*) ' select_case4_wamclimate'
!
! define F107/F107d/Kp indices for max/min/aver
!
      end subroutine select_case4_wamclimate
!
      end module namelist_wamphysics_def

!==================================================
      module idea_wam_control

      character(len=80) :: SPW_DRIVERS     ! four options a) 'swpc_fst' & b) 'cires_wam'
      character(len=80) :: swin_drivers    !  c) 'sair_wam', d) 'wam_climate;
      integer, parameter :: nlun_con = 133 ! unit for nam_wam_control
      character(len=80)  :: nml_control='wam_control_in'    ! in $RUNDIR
      end module idea_wam_control

      subroutine idea_wamcontrol_init(mpi_id)
!
!SK   use module_physics_driver, only: is_master
!     use idea_wam_control, only : SPW_DRIVERS, nlun_con, nml_control
      use idea_wam_control, only : SPW_DRIVERS, swin_drivers, nlun_con,
     & nml_control
      use namelist_wamphysics_def                            ! ./gsmphys/namelist_wamphysics_def.f
!
!SK   use idea_mpi_def,    only :   mpi_WAM_quit             !(iret, message)
      implicit NONE
      integer, intent(in) :: mpi_id
      integer :: ierr
      namelist /nam_wam_control/ 
     & wam_climate, wam_swpc_3day, wam_cires_rdata,wam_sair2012,
     & wam_swin, 
     & wam_smin, wam_smax,
     & wam_saver, wam_geostorm,
     & wam_gwphys, wam_solar_in, wam_ion_in, wam_das_in, wam_smet_in,
     & wam_netcdf_inout, wam_tides_diag, wam_pws_diag, wam_gws_diag,
     & JH0, JH_tanh, JH_semiann, JH_ann,
     & skeddy0, skeddy_semiann, skeddy_ann,
     & tkeddy0, tkeddy_semiann, tkeddy_ann
     

      open(nlun_con, file=trim(nml_control), status='old' )
      read(nlun_con, nam_wam_control, iostat=ierr)   
      close(nlun_con)
      if (ierr .ne. 0) then 
        print *, ' error in nam_wam_control '
        print *, ' file of nam_wam_control ', trim(nml_control)
        write(6, nam_wam_control)
!SK   call mpi_WAM_quit(23999, 'wam_control namelist_wamphysics_def.f')
      endif
       if (wam_climate)   SPW_DRIVERS = 'climate_wam'
       if (wam_swpc_3day) SPW_DRIVERS = 'swpc_fst'
       if (wam_cires_rdata) SPW_DRIVERS = 'cires_wam'
       if (wam_sair2012)    SPW_DRIVERS = 'sair_wam'
   
       if(wam_swin)   swin_drivers = 'swin_wam'

!SK    if (is_master) then
       if (mpi_id == 0) then
       print *, ' VAY idea_wamcontrol_init '
       print *, ' VAY SPW_DRIVERS ', SPW_DRIVERS
       print *, ' WAM_SWIN ', wam_swin
       print *,' nluncon=',nlun_con,'nam_wam_control =',
     &                 trim(nml_control)
       write(6, nam_wam_control)    
       endif
       end subroutine idea_wamcontrol_init

