! VAY 05/09  2015 First cut for WAM solar/geo drivers
! VAY Oct 18 2016 Updated to the NEMS-Legacy/WAM
!
!                solar input read ....F107 KP AP EUV(37)
!
! in "solar_in" - WAM namelist file under ../wam_trunk/job/solar_in
!===============================================================================
! idea_solar_fix = 2     True fixed values of EUV, KP, F107
! idea_solar_fix = 1     True fixed values of EUV and Time variable KP, F107
! idea_solar_fix = 0     Time variable all : EUV, KP, F107
!===============================================================================
       MODULE IDEA_SOLAR_INPUT
!
       use IDEA_IO_UNITS, only : iulog
       implicit none
       public :: solar_wam_get_feuv
       public :: solar_read_wam_init
       public :: solar_wamstep_advance
       public :: solar_waccmx_advance
       public :: solar_read_myy1947_2016
       public :: solar_readno_snoewx
       public :: solar_read_namelist
       public :: dealloc_solar
!
       save 
!
! no-snoe input Marsh et al. (2004)
!
       integer, parameter  ::  no_ny33=33
       integer, parameter  ::  no_nz16=16
       integer, parameter  ::  no_neofs=7       ! total NO-snoe model has 7 modes
       real, allocatable   ::  no_eof(:, : ,:)
       real , allocatable  ::  no_m(:,:)
       real, allocatable   ::  no_zkm(:), no_mlat(:)
!------------
       integer, parameter  :: itheia = 1     !vay-aug THEIA no-mpi broadcast 
!
       integer :: idea_solar_fix    ! flag for solar-variable (0) or fixed (1) solar conditions
                                    ! read_in from namelist "solar_in"
!
! Current time step geo-solar variables in WAM
! 
       integer, parameter :: nwafix= 37         ! dimension of EUV-bands
       real :: wf107_s, wf107a_s, wap_s, wkp_s, weuv_s(nwafix)
!
! Geo-solar variables from WAM-solar namelist
!
       real :: f107_fix, kp_fix, ap_fix, f107a_fix
       real :: Euv_fix(nwafix)
!
! solar calendar in the format  idat(yyyy, mm, dd) + real(hour)
!
       integer, parameter      ::  ndi = 4 ! idat and jdat -array dimensions
       integer, dimension(ndi) :: idatc    ! initial idat[ymdh]

       integer, dimension(ndi) :: jdat1    ! current closest data ymdh_d < ymdh
       integer, dimension(ndi) :: jdatc    ! current model ymdh
       integer, dimension(ndi) :: jdat2    ! current closest data ymdh_d > ymdh

       real :: hr1, hrc, hr2               ! corresponding real hours of jdat1-jdatc-jdat2

      
      integer ::  ndays
      integer ::  nwaves
      integer ::  ntimes
      integer ::  ntimes_wx                ! dimension (in days) of WACCM-X file
 
      integer,  allocatable :: dates(:)    ! for Kp and F107, days
      integer,  allocatable :: dfhours(:)  ! hours for Kp-sampling

 
      real, allocatable :: times(:)         ! data-based netcdf-file time-grid for data 20090101.frac_of_day
                                            ! it is used by WAM-SAIR for inday-cadence 1-hr, 3-hr etc...
                                            ! daily files like WACCMX has "0" UT, 20090101.00
                                            ! two readers of data
      real, allocatable :: Af107(:)
      real, allocatable :: Af107a(:)
      real, allocatable :: Akp(:)
      real, allocatable :: Aap(:)
!
      real, allocatable :: AEUV(:,:)


!
! Values for the current model step of WAM-physics      times(tim_ndx1) <= time_wam <= times(tim_ndx2)
!
       integer               :: tim_ndx1
       integer               :: tim_ndx2
       real                  :: w_ndx1 
       real                  :: w_ndx2 

       CONTAINS
!
      subroutine solar_readno_snoewx(file, mpi_id)
       use netcdf      
       use module_physics_driver, only : is_master    !SK 2020
!SK    use idea_mpi_def, ONLY:  info, mpi_comm_all
       implicit none
!SK    include 'mpif.h'
       character(len=*), intent(in) :: file
       integer, intent(in) :: mpi_id   
!locals
       integer :: nz, ny, neofs
       integer :: istat, ierr, astat
       integer :: dim_id, var_id 
       integer :: ncid, vid, iernc
       character(len=256) :: locfn
       integer, dimension(nf90_max_var_dims) :: dimidT  
!----------------------------------------------------------------------
!	... open the netcdf file
!----------------------------------------------------------------------
!SK    if(mpi_id.eq.0) then
        if (is_master) then
          write(iulog,*)file        
          write(iulog,*) 'solar_readno_snoewx: opening file for readno',
     &                    trim(file)
        endif
       ierNC=NF90_OPEN(trim(File), nf90_nowrite, ncid)   
       if (iernc /=0) write(iulog,*) ncid, 'ncid ', iernc, ' iernc ' 
!----------------------------------------------------------------------
!	... read the snoe dimensions
!----------------------------------------------------------------------
      iernc=nf90_inq_varid( ncid, 'EOF', vid )
         ierNC=nf90_inquire_variable(ncid, vid, dimids=dimidT) 
         iernc = nf90_inquire_dimension(ncid, dimidT(3), len=neofs)
         iernc = nf90_inquire_dimension(ncid, dimidT(2), len=nz)
         iernc = nf90_inquire_dimension(ncid, dimidT(1), len=ny)
!SK      if(mpi_id.eq.0) then
         if (is_master) then
          write(iulog,*) neofs, nz, ny, ' ne-nz-ny of NO-EOFs VAY'
!SK       write(iulog,*) no_neofs, no_nz16, no_ny33,' ne-nz-ny of NO-EOFs VAY'
!SK       if (nz .ne. no_nz16 .or. ny.ne.no_ny33 .or. neofs.ne.no_neofs) then
!SK       write(iulog,*)'snoe_rdeof: failed to read expected neofs=nz=ny'
          write(iulog,*) no_neofs, no_nz16, no_ny33,
     &                   ' ne-nz-ny of NO-EOFs VAY'
      if (nz .ne. no_nz16 .or. ny.ne.no_ny33 .or.neofs.ne.no_neofs) then
         write(iulog,*)'snoe_rdeof: failed to read expected neofs=nz=ny'
!SK      call mpi_quit(23901)
      endif
         endif
         
!----------------------------------------------------------------------
!	... allocate snoe variables
!----------------------------------------------------------------------
 
       allocate( no_mlat(ny),  no_zkm(nz),stat=astat )  
       allocate( no_m(ny, nz),stat=astat )
       allocate( no_eof(ny, nz, neofs),stat=astat )
       if( astat /= 0 ) then
       write(iulog,*) ' alloc_err in read_no_snoe no_eof ' 
!SK    write(iulog,*) 'snoe_rdeof: failed to allocate eofs; error = ',astat
       write(iulog,*) 'snoe_rdeof: failed to allocate eofs; error = ',
     &                astat
       end if 

!----------------------------------------------------------------------
!	... read the snoe variables
!----------------------------------------------------------------------
        iernc=nf90_inq_varid( ncid, 'lat', vid )
        iernc= nf90_get_var( ncid, vid, no_mlat)
        iernc=nf90_inq_varid( ncid, 'z', vid )
        iernc= nf90_get_var( ncid, vid, no_zkm)

      iernc = nf90_inq_varid( ncid, 'NO', var_id )
      ierr  = nf90_get_var( ncid, var_id, no_m )

      iernc = nf90_inq_varid( ncid, 'EOF', var_id )
      iernc = nf90_get_var( ncid, var_id, no_eof )  !(/1,1,1/), (/ny, nz, neofs/), no_eof )

!----------------------------------------------------------------------
!	... close the netcdf file
!----------------------------------------------------------------------
        iernc=nf90_close(ncid)     
!SK     if(mpi_id.eq.0) then
        if (is_master) then
         write(iulog,*) ' VAYsnoe ZKM:', no_zkm(1), ': ', no_zkm(nz)
         write(iulog,*) ' VAYsnoe MLT:', no_mlat(1), ': ', no_mlat(ny)
         write(iulog,*) ' VAYsnoe NO:', maxval(no_m), minval(no_m)
        endif
      end subroutine solar_readno_snoewx
!
       subroutine dealloc_solar(mpi_id)
       use module_physics_driver, only : is_master    !SK 2020
       integer :: mpi_id

       deallocate(af107, af107a, akp, aap)
       if( allocated(dates))          deallocate( dates)
       if( allocated(times))          deallocate( times)
       if( allocated(aeuv))           deallocate( aeuv)
       if( allocated(dfhours))        deallocate(dfhours)

!SK    if (mpi_id.eq.0) then
       if (is_master) then
       write(iulog, *)  'subroutine dealloc_solar: free memory '
       endif

       end subroutine dealloc_solar
!
       SUBROUTINE Found_jdates(mpi_id, ind_sol)
!
! find   [idat1, hr1, tim_ndx1], [idat2, hr2, tim_ndx2] for the start day
! update [idat1, hr1, tim_ndx1], [idat2, hr2, tim_ndx2] for 'the next call'
!
! call weights_time_interp(ndi, idat1, hr1, idat2, hr2, idatc, hrc, w_ndx1, w_ndx2)
!
        use wam_date_calendar, only : wam_split_ymd
        use module_physics_driver, only : is_master     !SK 2020
        implicit none
        integer, intent(in)    :: mpi_id
        integer, intent(inout) :: ind_sol
        integer :: wrk_date  
        real    :: wrk_time
        integer :: n, nk
        integer :: ymd1, ymd2, hh1, hh2, hhc
        hhc = Jdatc(4)
        wrk_date = 10000*Jdatc(1) + 100*Jdatc(2) + Jdatc(3) 
        wrk_time = float(wrk_date) + hrc/24.                    !wrk_time = flt_date( wrk_date, 0 )
!SK    if (mpi_id.eq.0) then        
       if (is_master) then        
        if( wrk_time < times(1) .or. wrk_time > times(ntimes) ) then
         write(iulog,*) 'solar_files: model time is out of-range F107'
!SK       write(iulog,*)   times(1) ,   times(ntimes) ,' times(start -/- end '
          write(iulog,*)   times(1) ,   times(ntimes) ,
     &                     ' times(start -/- end '
          write(iulog,*) wrk_time, wrk_date, ' wrk_time, wrk_date '
!         call endrun
        end if
       ENDIF

! time is growing
        nk =ind_sol
        do n = ind_sol,ntimes
           if( wrk_time <= times(n) ) then
              nk = n
           exit
          end if
       end do
!tim_ndx1
       ind_sol  = nk-1
       tim_ndx1 = nk-1
       tim_ndx2 = nk

       ymd1 = dates(tim_ndx1)
       ymd2 = dates(tim_ndx2)
       hh1 = 0               ! WACCMX data for UT=0
       hh2 = 0    
       hr1 = float(hh1)
       hr2 = float(hh2)
       call wam_split_ymd(ymd1, hh1, jdat1, ndi)
       call wam_split_ymd(ymd2, hh2, jdat2, ndi)
!
      end subroutine Found_Jdates
!
      subroutine solar_waccmx_advance(mpi_id, Mjdat_cur, hour_cur, kdt)

!      update calendar find "tim_ndx" and compute "weights-interp"
!      compute  "f107_s, f107a_s, ap_s, kp_s, euv_s(37)"
!      values for 'the current step of WAM"
!      times(1:NTIMES) - fixed "data-calendar" like "20130101.75"
!
       use wam_date_calendar, only : weights_time_interp
        use module_physics_driver, only : is_master     !SK 2020
       integer, intent(in) ::   kdt   
       integer, intent(in) ::   mpi_id          ! current mpi-PE
       integer, intent(in) ::   Mjdat_cur(ndi)  ! Mjdat for current step
       real ,   intent(in) ::   Hour_cur        ! Real fractional HOUR of timestep
       integer, save :: ind_dsol                ! the last index from which we search "SOLAR" file
!
! jdatc ( yr, mon, day, day_fraction) 

! define     [Jdat1, hr1, Jdat2, hr2, Jdatc, hrc]
!             next 3-lines

        if (kdt == 1)  ind_dsol =1     ! or postion for year of 20090101
        Jdatc = Mjdat_cur
        hrc   = hour_cur
        CALL Found_jdates(mpi_id, ind_dsol)         ! "Found_jdates" to "bound" Jdatc by Jdat1/Jdat2
!                                               
!       
!      "weights_time_interp" works with Julian Days 
!
!SK     CALL weights_time_interp(ndi, Jdat1, hr1, Jdat2, hr2, Jdatc, hrc, w_ndx1, w_ndx2)
        CALL weights_time_interp(ndi, Jdat1, hr1, Jdat2, hr2, Jdatc, 
     &                           hrc, w_ndx1, w_ndx2)
!
!       
        if (idea_solar_fix.eq.1) CALL solar_wam_get_f107kp    !( f107_s, f107a_s, ap_s, kp_s)
!
!SK    if (mpi_id.eq.0) then
       if (is_master) then
!SK    write(iulog,*) ' sol_adv indices/Weights: ', tim_ndx1, tim_ndx2, w_ndx1, w_ndx2
!SK    write(iulog,*) ' VAY-sol_adv dates:Jdat1', Jdat1(1), Jdat1(2), Jdat1(3)
!SK    write(iulog,*) ' VAY-sol_adv dates:Jdat1', Jdat2(1), Jdat2(2), Jdat2(3)  
       write(iulog,*) ' sol_adv indices/Weights: ', tim_ndx1, tim_ndx2,
     &                w_ndx1, w_ndx2
       write(iulog,*) ' VAY-sol_adv dates:Jdat1', Jdat1(1), Jdat1(2),
     &                Jdat1(3)
       write(iulog,*) ' VAY-sol_adv dates:Jdat1', Jdat2(1), Jdat2(2),
     &                Jdat2(3)  
       write(iulog,*) 'solar_waccmx_adv f107_s,  kp_s',  wf107_s, wkp_s
       endif
       end subroutine solar_waccmx_advance
!
       subroutine solar_wamstep_advance(mpi_id, Mjdat_cur, hour_cur)

!      update calendar find "tim_ndx" and compute "weights-interp"
!      compute  "f107_s, f107a_s, ap_s, kp_s, euv_s(37)"
!      values for 'the current step of WAM"
!      times(1:NTIMES) - fixed "data-calendar" like "20130101.75"
!
       use wam_date_calendar, only : weights_time_interp
       use module_physics_driver, only : is_master    !SK 2020
       integer ::  mpi_id             ! current mpi-PE
       integer ::  Mjdat_cur(ndi)     ! Mjdat for current step
       real    ::  Hour_cur           ! Real fractional HOUR of timestep
! jdatc ( yr, mon, day, day_fraction) 

! define     [Jdat1, hr1, Jdat2, hr2, Jdatc, hrc]
!             next 3-lines

        Jdatc = Mjdat_cur
        hrc   = hour_cur
        CALL start_jdates(mpi_id)      ! more adequate name "update_jdates" to "bound" Jdatc by Jdat1/Jdat2
!                                      ! for same YEAR YYYY you can use simple LINTERP in time
                                       ! for NY and OY you can use Julian Day Transform   
!       after update_jdates compute "w_ndx1, w_ndx2" and call  "solar_wam_get_feuv" or "solar_wam_get_f107kp"
!
!        "weights_time_interp" works with Julian Days 
!
!SK     CALL weights_time_interp(ndi, Jdat1, hr1, Jdat2, hr2, Jdatc, hrc, w_ndx1, w_ndx2)
        CALL weights_time_interp(ndi, Jdat1, hr1, Jdat2, hr2, Jdatc, 
     &                           hrc, w_ndx1, w_ndx2)
!
!vay 08/2015  idea_solar_fix -FLAG for time interpolation
!       
       if (idea_solar_fix.eq.0) CALL solar_wam_get_feuv      !( f107_s, f107a_s, ap_s, kp_s, euv_s)
       if (idea_solar_fix.eq.1) CALL solar_wam_get_f107kp    !( f107_s, f107a_s, ap_s, kp_s)
!
!SK    if (mpi_id.eq.0) then
       if (is_master) then
!SK    write(iulog,*) ' sol_adv indices/Weights: ', tim_ndx1, tim_ndx2, w_ndx1, w_ndx2
       write(iulog,*) ' sol_adv indices/Weights: ', tim_ndx1, tim_ndx2,
     &                w_ndx1, w_ndx2
       write(iulog,*) ' f107_s,  kp_s', wf107_s, wkp_s
       endif
       end subroutine solar_wamstep_advance
!
       SUBROUTINE start_jdates(mpi_id)
!
! find   [idat1, hr1, tim_ndx1], [idat2, hr2, tim_ndx2] for the start day
! update [idat1, hr1, tim_ndx1], [idat2, hr2, tim_ndx2] for 'the next call'
!       call weights_time_interp(ndi, idat1, hr1, idat2, hr2, idatc, hrc, w_ndx1, w_ndx2)
!
!        integer :: idatc(4)
        use wam_date_calendar, only : wam_split_ymd
        use module_physics_driver, only : is_master    !SK 2020
        implicit none
        integer :: mpi_id
        integer :: wrk_date  
        real    :: wrk_time
        integer :: n, nk
        integer :: ymd1, ymd2, hh1, hh2, hhc
        hhc = Jdatc(4)
        wrk_date = 10000*Jdatc(1) + 100*Jdatc(2) + Jdatc(3) 
        wrk_time = float(wrk_date) + hrc/24.                    !wrk_time = flt_date( wrk_date, 0 )
!SK    if (mpi_id.eq.0) then        
       if (is_master) then     !SK        
        if( wrk_time < times(1) .or. wrk_time > times(ntimes) ) then
         write(iulog,*) 'solar_files: model time is out of-range F107'
!SK       write(iulog,*)   times(1) ,   times(ntimes) ,' times(start -/- end '
          write(iulog,*)   times(1) ,   times(ntimes) ,
     &                    ' times(start -/- end '
          write(iulog,*) wrk_time, wrk_date, ' wrk_time, wrk_date '
!         call endrun
        end if
       ENDIF

! time is growing
        nk =2
        do n = 2,ntimes
           if( wrk_time <= times(n) ) then
              nk = n
           exit
          end if
       end do
!tim_ndx1
       tim_ndx1 = nk-1
       tim_ndx2 = nk
!SK    w_ndx2=(wrk_time-times(tim_ndx1))/(times(tim_ndx2)-times(tim_ndx1))
       w_ndx2=(wrk_time-times(tim_ndx1))/
     &        (times(tim_ndx2)-times(tim_ndx1))
       w_ndx1 = 1. -w_ndx2 
!
! we need updates of hr1 & hr2
!       hr1 = (times(tim_ndx1)-float(yddds(tim_ndx1)))*24.
!       hr2 = (times(tim_ndx2)-float(yddds(tim_ndx2)))*24.

       ymd1 = dates(tim_ndx1)
       ymd2 = dates(tim_ndx2)
       if (hr2.lt.hrc.and.(ymd1 == ymd2)) hr2 = hrc
       hh1 = nint(hr1)
       hh2 = nint(hr2)    

       call wam_split_ymd(ymd1, hh1, jdat1, ndi)
       call wam_split_ymd(ymd2, hh2, jdat2, ndi)
!
       end subroutine start_Jdates
!

       subroutine solar_read_namelist(nml_solar, nlun_solar, 
     &                                ncfile_fpath, file_no, mpi_id)
!SK    subroutine solar_read_namelist(nml_solar, nlun_solar, ncfile_fpath, file_no, mpi_id)
!
! read name-list
!
      use module_physics_driver, only : is_master
      integer :: mpi_id
      character(len=*),intent(in)   :: NmL_solar
      integer, intent(in)           :: nlun_solar
      character(len=*), intent(out) :: ncfile_fpath, file_no

      integer, parameter :: ch100 = 256
      integer, parameter :: me =0
      integer, parameter :: masterproc =0

      INTEGER :: k, i,j
      integer :: unitn, ierr  
      integer          :: isolar_file   
      character(ch100) :: Dir_swpc       !='/scratch1/portfolios/NCEPDEV/swpc/save/'
      character(ch100) :: Dir_uid        !='Valery.Yudin/NEMS/wam_april_2015/data_euv/'
      character(ch100) :: Dir_solar      !=Dir_swpc//Dir_uid
      character(ch100) :: Solar_file     ! name of file w/o the full-path
      character(ch100) :: noeof_file     ! name of file w/o the full-path
      character(ch100) :: wxdan_file     !
      character(ch100) :: wam2012_file   !
!
      namelist /solar_parms_nl/ idea_solar_fix, solar_file, Dir_swpc,
     &          Dir_uid, noeof_file, wxdan_file, wam2012_file,
     &          isolar_file,f107_fix, f107a_fix, kp_fix, ap_fix, Euv_fix
!SK   namelist /solar_parms_nl/ idea_solar_fix, solar_file, Dir_swpc, Dir_uid, 
!SK  &        noeof_file, wxdan_file, wam2012_file, isolar_file,
!SK  &        f107_fix, f107a_fix, kp_fix, ap_fix, Euv_fix
!=========================================================================
! solar_in should be copied to $RUNDIR along with other namelists
!
! additional 
! ideaphys_in will be next to "make" controls for WAM-physics 
!=========================================================================

      isolar_file = 2016       ! all years
      wxdan_file='wasolar_dan_20161019.nc'
      noeof_file='snoe_eof.nc'
! put defaults above
!
      open(nlun_solar, file=trim(nml_solar), status='old' )
      read(nlun_solar, solar_parms_nl, iostat=ierr)   
      close(nlun_solar)    

      Dir_solar=trim(Dir_swpc)//trim(Dir_uid)
      ncfile_fpath= trim(Dir_solar)//trim(solar_file)
!
! select b-n   wxdan_file & wam2012_file
!
      ncfile_fpath= trim(wxdan_file)
      if (isolar_file == 2012) ncfile_fpath= trim(wam2012_file)
!SK   if (mpi_id.eq.0) then     
      if (is_master) then     
      write(iulog,*) idea_solar_fix, 'idea_solar_fix - flag (1-fixed) '
      write(iulog,*)  f107_fix, f107a_fix, ' F107 '
      write(iulog,*)  kp_fix, ap_fix,  ' Kp-Ap '
      write(iulog,*) Euv_fix(1), Euv_fix(37), ' EUV (1:37) '

      write(iulog,*)    '+++++++++++ ncfile_fpath '
      write(iulog,133) ncfile_fpath
      endif
133   format(A122)
!                                    time-invariant solar and geo- parameters
      if (idea_solar_fix == 2) then
       wf107_s  = f107_fix
       wf107a_s = f107a_fix
       wkp_s    = kp_fix
       wap_s    = ap_fix
       weuv_s   = Euv_fix 
!SK    if (mpi_id.eq.0) then    
       if (is_master) then    
!SK        write(iulog,*) idea_solar_fix, ' Time-invariant  Solar-Geo Inputs'
           write(iulog,*) idea_solar_fix, 
     &         ' Time-invariant  Solar-Geo Inputs'
           write(iulog,*) wf107_s, ' F107-fix '
           write(iulog,*) wKp_s, ' Kp-fix '
           write(iulog,*) wAp_s, ' Ap-fix '
       endif
      endif
      if (idea_solar_fix == 1) then
       weuv_s   = Euv_fix 
       if (is_master) then    
!SK    if (mpi_id.eq.0) then    
!SK        write(iulog,*) idea_solar_fix,' Time-variable F107/KP; fixed EUV-inputs'
           write(iulog,*) idea_solar_fix,
     &    ' Time-variable F107/KP; fixed EUV-inputs'
!           write(iulog,*) wf107_s, ' F107-fix '
!           write(iulog,*) wKp_s, ' Kp-fix '
!           write(iulog,*) wAp_s, ' Ap-fix '
       endif
      endif
      if (idea_solar_fix == 0) then
       weuv_s   = Euv_fix 
       if (is_master) then
!SK    if (mpi_id.eq.0) then    
!SK        write(iulog,*) idea_solar_fix,' Time-variable EUV-inputs with F107/KP'
           write(iulog,*) idea_solar_fix,
     &   ' Time-variable EUV-inputs with F107/KP'
!
       endif
      endif
!
! full path to NO-eof file
!
      file_NO = trim(noeof_file)
!
      end subroutine solar_read_namelist
!
!
!       
       subroutine solar_wam_get_feuv   !( f107_s, f107a_s, ap_s, kp_s, euv_s )
       implicit none
!local
        integer  :: k

        wf107_s  =  af107(tim_ndx1)*w_ndx1 + af107(tim_ndx2)*w_ndx2
!
         wf107a_s  = af107a(tim_ndx1)*w_ndx1 + af107a(tim_ndx2)*w_ndx2
 
!  
         wkp_s  =  akp(tim_ndx1)*w_ndx1 + akp(tim_ndx2)*w_ndx2 
! 
         wap_s  =  aap(tim_ndx1)*w_ndx1 + aap(tim_ndx2)*w_ndx2  
!   
        do k=1, nwafix
!SK        weuv_s(k)  =  Aeuv(k, tim_ndx1)*w_ndx1 + Aeuv(k, tim_ndx2)*w_ndx2
           weuv_s(k)  =  Aeuv(k, tim_ndx1)*w_ndx1 + Aeuv(k, tim_ndx2)
     &                  *w_ndx2
        enddo
!
       end subroutine solar_wam_get_feuv
!
       subroutine solar_wam_get_f107kp   !( f107_s, f107a_s, ap_s, kp_s)
       implicit none
!local
        integer  :: k

        wf107_s  =  af107(tim_ndx1)*w_ndx1 + af107(tim_ndx2)*w_ndx2
!
         wf107a_s  = af107a(tim_ndx1)*w_ndx1 + af107a(tim_ndx2)*w_ndx2
 
!  
         wkp_s  =  akp(tim_ndx1)*w_ndx1 + akp(tim_ndx2)*w_ndx2 
! 
         wap_s  =  aap(tim_ndx1)*w_ndx1 + aap(tim_ndx2)*w_ndx2  

!
       end subroutine solar_wam_get_f107kp
!
!
!
       subroutine solar_read_wam_init ( file, mpi_id )
!
! called from idea_solar_heating.f:CALL solar_read_wam_init(ncfile_fpath, mpi_id)   
!
       use netcdf      
       use module_physics_driver, only : is_master
!SK    use idea_mpi_def, ONLY:  info, mpi_comm_all
     
       implicit none
       include 'mpif.h'
!input
        integer :: mpi_id
        character(len=*) :: file
!
!locals
        integer ::  ierr
        integer ::  ncid, vid, ierNC
        integer  :: astat        
             
        real     :: wrk_time   
        integer  :: wrk_date
        integer  :: yr, mon, day, day_fraction
        integer  :: dimid    
        integer, dimension(nf90_max_var_dims) :: dimidT         
        integer :: n
        integer :: masterproc

!SK     if(mpi_id.eq.0) then
        if(is_master) then
           write(iulog,*)file        
           write(iulog,*) 'SOLAR_PARMS: opening file ', trim(file) 
        endif
!
       ierNC=NF90_OPEN(trim(File), nf90_nowrite, ncid)   
!SK    if (iernc /=0) write(iulog,*) ncid, 'ncid ', iernc, ' iernc '
       if (is_master.and.iernc /=0) 
     &     write(iulog,*) ncid, 'ncid ', iernc, ' iernc '

!       call cam_pio_openfile ( ncid, locfn, PIO_NOWRITE)    
!       ierr = pio_inq_dimid( ncid, 'time', dimid )
!       ierr = pio_inq_dimlen( ncid, dimid, ntimes )
!        ierr = pio_inq_varid( ncid, 'date', varid )
!        ierr = pio_get_var( ncid, varid, dates )
!

         iernc=nf90_inq_varid( ncid, 'EUV', vid )
!SK      if (iernc /=0) write(iulog,*) ncid, 'ncid ', iernc, ' iernc '
         if (is_master.and.iernc /=0) write(iulog,*) ncid, 'ncid ',
     &     iernc, ' iernc '
         ierNC=nf90_inquire_variable(ncid, vid, dimids=dimidT) 
         iernc = nf90_inquire_dimension(ncid, dimidT(2), len=ntimes)
         iernc = nf90_inquire_dimension(ncid, dimidT(1), len=nwaves)

!SK      if(mpi_id.eq.0) then
         if(is_master) then
              write(iulog,*) ntimes, nwaves, ' nt-nw  idea_solar_input'
         endif
         
 !   
       allocate( dates(ntimes),  times(ntimes),stat=astat )  
       allocate( dfhours(ntimes),stat=astat )  
       if( is_master.and.astat /= 0 ) then
!SK    if( astat /= 0 ) then
!SK    write(iulog,*) ' alloc_err in read_waccm_solar for dates,times', ntimes 
       write(iulog,*) ' alloc_err in read_waccm_solar for dates,times',
     &   ntimes 
       end if    


!         iernc = nf90_inquire_dimension(ncid, dimidT(1), len=nlev)
!         iernc = nf90_inquire_dimension(ncid, dimidT(2), len=npix)  

        iernc=nf90_inq_varid( ncid, 'date', vid )
        iernc= nf90_get_var( ncid, vid, dates)

          if(is_master) then
!SK       if(mpi_id.eq.0) then
            write(iulog,*) dates, ' dates ' 
          endif

        do n = 1,ntimes
           dfhours(n) = 0          ! integer .......current for daily  12UT
           times(n) = float(dates(n))  + dfhours(n)/24.    
        end do
!
! init hr1 & hr2
          hr1 = dfhours(1)
          hr2 = dfhours(2)
    !---------------------------------------------------------------
    !	... allocate and read solar parms ..... ALL-time series
    !   call dealloc_solar(mpi_id) in the End of WAM-RUN
    !   we do not put these data-sets in the restart files
    !---------------------------------------------------------------
       allocate(  Af107(ntimes), Af107a(ntimes),stat=astat )
       allocate(  Akp(ntimes),  Aap(ntimes),    stat=astat )
       allocate(  AEUV(nwaves,ntimes), stat=astat )
       if( is_master.and.astat /= 0 ) then
!SK    if( astat /= 0 ) then
         write(iulog,*) ' alloc_err( astat, f107 ... ap ', ntimes 
       end if

        iernc=nf90_inq_varid( ncid, 'f107', vid )
        iernc= nf90_get_var( ncid, vid, Af107)

        iernc=nf90_inq_varid( ncid, 'f107a', vid )
        iernc= nf90_get_var( ncid, vid, Af107a)

        iernc=nf90_inq_varid( ncid, 'EUV', vid )
        iernc= nf90_get_var( ncid, vid, AEUV)


        iernc=nf90_inq_varid( ncid, 'kp', vid )
        iernc= nf90_get_var( ncid, vid, Akp)

        iernc=nf90_inq_varid( ncid, 'ap', vid )
        iernc= nf90_get_var( ncid, vid, Aap)
        iernc=nf90_close(ncid)     
!
!
          if(is_master) then
!SK       if(mpi_id.eq.0) then
          write(iulog,*) '  read_wam_solar: ntimes  ', ntimes   
          write(iulog,*)     maxval(af107),   minval(af107), ' F107 '
          write(iulog,*)     maxval(af107a),   minval(af107), ' F107a ' 
          write(iulog,*) maxval(AEUV), minval(AEUV), ' EUV ', nwaves

          write(iulog,*)     maxval(aKp),   minval(aKp), ' Kp-daily ' 
          write(iulog,*)     maxval(aAp),   minval(aAp), ' Aap '      
          write(iulog,*)            ' mpi_bcast in solar_read_wam_init'
          write(iulog,*)  ' VAY completed solar_read_wam_init'
         endif

        RETURN    ! Here RETURN is a temporary FIX of mpif.h MPI_REAL8/mpi_integer for THEIA
                  ! ALL PEs read nc-file
!
! result of solar_read_wam_init(ncfile_fpath): ALLOCATE, READ and BROADCAST
! call MPI_BCAST (data_to_be_sent, send_count, send_type, broadcasting_process_ID, comm, ierr)

!SK    call mpi_bcast(ntimes,1,mpi_integer,0,mpi_comm_all,info)
!      call mpi_bcast(nwaves,1,mpi_integer,0,mpi_comm_all,info)

!      call mpi_bcast(dates,ntimes,mpi_integer,0,mpi_comm_all,info)
!      call mpi_bcast(dfhours,ntimes, mpi_integer,0,MPI_COMM_ALL,info)

!      call mpi_bcast(times  ,ntimes,MPI_REAL8,0,MPI_COMM_ALL,info)

!      call mpi_bcast(af107  ,ntimes,MPI_REAL8,0,MPI_COMM_ALL,info)
!      call mpi_bcast(af107a  ,ntimes,MPI_REAL8,0,MPI_COMM_ALL,info)

!      call mpi_bcast(akp  ,ntimes,MPI_REAL8,0,MPI_COMM_ALL,info)
!      call mpi_bcast(aap  ,ntimes,MPI_REAL8,0,MPI_COMM_ALL,info)
!      call mpi_bcast(aEUV,nwaves*ntimes,MPI_REAL8,0,MPI_COMM_ALL,info)

!       call mpi_barrier(mpi_comm_all,info)         
!      if(mpi_id.eq.0) then
!         write(iulog,*)  ' VAY completed solar_read_wam_init'
!      endif
       end  subroutine solar_read_wam_init
!
!
       SUBROUTINE solar_read_myy1947_2016( file, mpi_id )
!========================================================
! Oct 20 2016:  VAY include WACCM-solar data file
! data span from 19470410 to 20160723
! /scratch3/NCEPDEV/swpc/save/Valery.Yudin/BASE_SVN/BASE_WAM_DATA/WACCM_DAILY_2016
! file: wasolar_dan_20161019.nc created by Daniel Marsh, Oct 19/2016 
!	time = UNLIMITED ; // (24671 currently)
!variables:
!	float f107(time) ;
!		f107:long_name = "10.7 cm solar radio flux (F10.7)" ;
!		f107:units = "10^-22 W m^-2 Hz^-1" ;
!	float f107a(time) ;
!		f107a:long_name = "81-day centered mean of 10.7 cm solar radio flux (F10.7)" ;
!	float kp(time) ;
!		kp:long_name = "Daily planetary K index" ;
!	short ap(time) ;
!		ap:long_name = "Daily planetary a index" ;
!		ap:units = "nanoTeslas" ;
!	short isn(time) ;
!		isn:long_name = "International Sunspot Number" ;
!	int date(time) ;
!		date:long_name = "current date (YYYYMMDD)" ;
!=======================================================
       use netcdf      
!SK    use idea_mpi_def, ONLY:  info, mpi_comm_all     
       use module_physics_driver, only : is_master
       implicit none
!SK    include 'mpif.h'
!input
        integer :: mpi_id
        character(len=*) :: file
!
!locals
!
        integer ::  ierr
        integer ::  ncid, vid, ierNC
        integer  :: astat      
!        integer  :: ntimes_wx                ! data-line of FILE, # of days  as a pert of module
             
        real     :: wrk_time   
        real, parameter     :: r24 = 1./24.
        integer  :: wrk_date
        integer  :: yr, mon, day, day_fraction
        integer  :: dimid    
        integer, dimension(nf90_max_var_dims) :: dimidT         
        integer  :: n
        integer  :: masterproc
        integer(2), allocatable  :: Apshort(:)

        if(is_master) then
!SK     if(mpi_id.eq.0) then
           write(iulog,*)file        
           write(iulog,*) 'VAY SOLAR_MULTI-YEARS: opening file ',
     &      trim(file) 
        endif
!
       ierNC=NF90_OPEN(trim(File), nf90_nowrite, ncid)   
!SK    if (iernc /= 0) write(iulog,*) ncid, 'ncid ', iernc, ' iernc '
       if (is_master.and.iernc /= 0) 
     &    write(iulog,*) ncid, 'ncid ', iernc, ' iernc '

!       call cam_pio_openfile ( ncid, locfn, PIO_NOWRITE)    
!       ierr = pio_inq_dimid( ncid, 'time', dimid )
!       ierr = pio_inq_dimlen( ncid, dimid, ntimes )
!        ierr = pio_inq_varid( ncid, 'date', varid )
!        ierr = pio_get_var( ncid, varid, dates )
!

         iernc=nf90_inq_varid( ncid, 'f107', vid )
!SK      if (iernc /=0) write(iulog,*) 'err ind_varid ', iernc, ' f107 '
         if (is_master.and.iernc /=0) 
     &      write(iulog,*) 'err ind_varid ', iernc, ' f107 '

         ierNC=nf90_inquire_variable(ncid, vid, dimids=dimidT) 
         iernc = nf90_inquire_dimension(ncid, dimidT(1), len=ntimes_wx)
!         iernc = nf90_inquire_dimension(ncid, dimidT(1), len=nwaves)
         ntimes = ntimes_wx
         if(is_master) then
!SK      if(mpi_id.eq.0) then
!SK           write(iulog,*) ntimes_wx, ' ntimes_wx  idea_solar_input_myy'
            write(iulog,*) ntimes_wx, ' ntimes_wx  idea_solar_input_myy'
         endif
         
 !   
       allocate( dates(ntimes_wx),  times(ntimes_wx),stat=astat )  
       allocate( dfhours(ntimes_wx),stat=astat )  
       if( is_master.and.astat /= 0 ) then
!SK    if( astat /= 0 ) then
!SK    write(iulog,*) ' alloc_err in read_waccm_solar for dates/times', ntimes_wx 
       write(iulog,*) ' alloc_err in read_waccm_solar for dates/times',
     &     ntimes_wx 
       end if    

 

        iernc=nf90_inq_varid( ncid, 'date', vid )
        iernc= nf90_get_var( ncid, vid, dates)

          if(is_master) then
!SK       if(mpi_id.eq.0) then
            write(iulog,*) ' dates-last 365 days ' 
            write(iulog,*) dates(ntimes-365:ntimes)
            write(iulog,*) ' dates-last 365 days ' 
          endif

        do n = 1,ntimes_wx
           dfhours(n) = 0                              ! integer .......current for daily  12UT
           times(n) = float(dates(n))  + dfhours(n)*r24   
        end do
!
! init hr1 & hr2
          hr1 = dfhours(1)
          hr2 = dfhours(2)
    !---------------------------------------------------------------
    !	... allocate and read solar parms ..... ALL-time series
    !   call dealloc_solar(mpi_id) in the End of WAM-RUN
    !   we do not put these data-sets in the restart files
    !---------------------------------------------------------------
       allocate(  Af107(ntimes_wx), Af107a(ntimes_wx),stat=astat )
       allocate(  Akp(ntimes_wx),  Aap(ntimes_wx),    stat=astat )
       allocate(  Apshort(ntimes_wx)) 

       if( is_master.and.astat /= 0 ) then
!SK    if( astat /= 0 ) then
         write(iulog,*) ' alloc_err( astat, f107 ... ap ', ntimes_wx 
       end if

        iernc=nf90_inq_varid( ncid, 'f107', vid )
        iernc= nf90_get_var( ncid, vid, Af107)

        iernc=nf90_inq_varid( ncid, 'f107a', vid )
        iernc= nf90_get_var( ncid, vid, Af107a)

        iernc=nf90_inq_varid( ncid, 'kp', vid )
        iernc= nf90_get_var( ncid, vid, Akp)

        iernc=nf90_inq_varid( ncid, 'ap', vid )
        iernc= nf90_get_var( ncid, vid, Apshort)
        iernc=nf90_close(ncid)     
!
!
          Aap = float(Apshort)
          if(is_master) then
!SK       if(mpi_id.eq.0) then
          write(iulog,*) '  solar_read_myyread, ntimes_wx:  ', ntimes_wx
          write(iulog,*)     maxval(af107), minval(af107), ' F107 '
          write(iulog,*)     maxval(af107a),minval(af107a), ' F107a ' 
          write(iulog,*)     maxval(aKp),   minval(aKp), ' Kp-daily ' 
          write(iulog,*)     maxval(aAp),   minval(aAp), ' Aap-daily ' 
          write(iulog,*)            ' mpi_bcast in solar_read_wam_init'
          write(iulog,*)  ' VAY completed solar_read_WACCMX_init'
         endif
        
       if(is_master) then
!SK    if(mpi_id.eq.0) then
!SK       write(iulog,*)  ' VAY completed solar_read_myy1947_2016 ntimes', ntimes
          write(iulog,*)' VAY completed solar_read_myy1947_2016 ntimes',
     &                   ntimes
       endif
        RETURN    ! Here RETURN is a temporary FIX of mpif.h MPI_REAL8/mpi_integer for THEIA
                  ! ALL PEs read nc-file
!
! result of solar_read_wam_init(ncfile_fpath): ALLOCATE, READ and BROADCAST
! call MPI_BCAST (data_to_be_sent, send_count, send_type, broadcasting_process_ID, comm, ierr)
      
!SK    call mpi_bcast(ntimes,1,mpi_integer,0,mpi_comm_all,info)
!      call mpi_bcast(nwaves,1,mpi_integer,0,mpi_comm_all,info)

!      call mpi_bcast(dates,ntimes,mpi_integer,0,mpi_comm_all,info)
!      call mpi_bcast(dfhours,ntimes, mpi_integer,0,MPI_COMM_ALL,info)

!      call mpi_bcast(times  ,ntimes,MPI_REAL8,0,MPI_COMM_ALL,info)

!      call mpi_bcast(af107  ,ntimes,MPI_REAL8,0,MPI_COMM_ALL,info)
!      call mpi_bcast(af107a  ,ntimes,MPI_REAL8,0,MPI_COMM_ALL,info)

!      call mpi_bcast(akp  ,ntimes,MPI_REAL8,0,MPI_COMM_ALL,info)
!      call mpi_bcast(aap  ,ntimes,MPI_REAL8,0,MPI_COMM_ALL,info)
 
!       call mpi_barrier(mpi_comm_all,info)         
!      if(mpi_id.eq.0) then
!         write(iulog,*)  ' VAY completed solar_read_myy1947_2016 '
!SK    endif
!
!
       END SUBROUTINE solar_read_myy1947_2016
!
       END MODULE IDEA_SOLAR_INPUT      
