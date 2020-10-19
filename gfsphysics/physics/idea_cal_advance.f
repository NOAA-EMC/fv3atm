      MODULE wam_date_calendar
!================================================================== 
! 
! vay-2015/16 WAM YYYYMMDD_HH calendar for SOLAR-GEO inputs
!
!==================================================================   
! Collections of WAM date-calendar subroutines
!   to work properly with the time-dependent 
!   SOLAR and METEO input/output fiels
! 
! Accepted DATE-TIME format IDAT(4) + FHOUR ( from IC-start
!                           JDAT(4) + CHOUR ( CHOUR -digital from current day)
!
!      subroutine test_wam_calendar
!
!      SUBROUTINE CURRENT_NCEP_JDAT(idat,hr_after_idat,jdat, hr_jdat)
!
!      subroutine wam_split_ymd(ymd, hh, jdat, ndi)
!
!      subroutine wam_calendar(ndi, jdat1,hr1, jdat2, hr2, jdatc, hrc)
!
!      subroutine weights_time_interp(ndi, Jdat1, hr1, Jdat2, hr2, Jdatc, hrc, w1, w2)
!
!      subroutine wam_julday_doy(idat, jdoy,jday)
!
!      SUBROUTINE W3FS26_VY(JLDAYN,IYEAR,MONTH,IDAY,IDAYWK,IDAYYR)
!      
!      SUBROUTINE JULDAY_2_YMD(JLDAYN,IYEAR,MONTH,IDAY)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Still in write_phys_nems idate[hr, month, day, year]
!
!      iyr     = idate(4)
!      imo     = idate(2)
!      ida     = idate(3)
!      ihr     = idate(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! idea_solar_heating.f:      use date_def
! in idea_solar_heating.f pre....idat(8) for w3lib
!                                idate(4) from date_def
!      idat=0
!      idat(1)=idate(4) year
!      idat(2)=idate(2) mont
!      idat(3)=idate(3) day
!      idat(5)=idate(1) ihr
!      utsec=solhr*3600   ...... "solhr" input for presolar
!====================================================================
!    
! idea_ion.f:!     use date_def
! idea_ion.f:      use date_def
!      subroutine idea_geteb(im,ix,dayno,utsec,f107,kp,maglat,maglon,    &
!     &essa,ee1,ee2)
!
!        iday = dayno                   ! day of year .... variable
!        imo   = idate(2) ..............! initial day of integration
!        iday_m= idate(3)               ! initial month of integration
!        iyear = 1995                   !
!        f107d=f107
!
! WAM/NEMS employ   "use date_def" to get a time-date for solar_subs
!====================================================================
     
      use idea_io_units, only : iulog
!==================================================================     
       implicit none

      integer, parameter :: ndwam = 4
      integer            :: curdate_wam(ndwam)  ! idat-format
      real               :: curhour_wam         ! real hour 1.35 = 1hour+ 35/100
      integer            :: curyear_wam         ! current model year
      integer            :: curmonth_wam        ! current model month
      integer            :: curddd_wam          ! current day of year (1:366)
      integer            :: curday_wam          ! current day oo month (1:31)
      integer            :: curjulday_wam       ! current Julian day      ! 
      integer            :: curutsec_wam        ! current # of seconds   nint(curhour_wam)*3600 ! 
      integer            :: curnsteps_wam       ! current number of timesteps
      real               :: curdtime_wam        ! current model time step in sec
      integer            :: idat_wam(4)         ! idat-format (yyyy-month-dd-hh)
      real               :: irhour_wam          ! initial real hour 1.35 = 1hour+ 35/100
      CONTAINS
      SUBROUTINE COPY_IDAT_NEMS_2_WAM
!
! call from "sub-ne idea_solar_init" of "idea_solar_heating.f" in gloopb.f
!     idat_nems = Idate_header(4-year;3-day;2-month;1-ihour)
!     idar_wam  =             [1-year/2-month/3-day/4-hour]
!
      use date_def, only : idat_nems => IDATE   ! nems_format (hr-month-dd-YYYY) according to HEADER-file
      use date_def, only : fhour                ! nems digital hour
                                                ! check it out
      idat_wam(1)=idat_nems(4)                  ! year
      idat_wam(2)=idat_nems(2)                  ! month
      idat_wam(3)=idat_nems(3)                  ! day of month
      idat_wam(4)=idat_nems(1)                  ! hour of initial day
      irhour_wam =   fhour     

!       print *, 'vay-cal year',   idat_wam(1), idat_nems(4)  
!       print *, 'vay-cal mth',    idat_wam(2), idat_nems(2) 
!       print *, 'vay-cal day',    idat_wam(3), idat_nems(3)  
         
      END SUBROUTINE COPY_IDAT_NEMS_2_WAM
!
      SUBROUTINE CURRENT_NCEP_JDAT(idat,hr_after_idat,jdat, hr_jdat)
!
!     refresh all current calendar-related "cur..._wam scalars"
! 
!    [Fhour, IDAT(4)] =>  [hr_jdat, JDAT(4)] current YMD + Real(Hour)
!
!     transform (idat + hr_after_idat) => jdat (current day-format)
!
      real :: hr_after_idat      ! how many hours from IC or IDAT
      integer :: idat(4)
      integer :: jdat(4)         ! current YY MM DD HH
      real :: hr_jdat            ! real HH witg time-step precision 
!
      integer :: hrin, hri_dres, iddd
      integer ::jday, jday0101, jday_ycur0101
      hrin = floor(hr_after_idat)
      hri_dres = modulo(hrin,24)
      hr_jdat = hri_dres+ (hr_after_idat - hrin)

      iddd = (hrin - hri_dres)/24
      JDAY0101 = iw3jdn_VY(idat(1),1,1) 
      JDAY = iw3jdn_VY(idat(1),idat(2),idat(3)) + iddd 
      jdat(4) = hri_dres 
!
! Back from Julian to Normal
!   
      CALL JULDAY_2_YMD(JDAY,JDAT(1),JDAT(2),JDAT(3)) 
!
! Populate the PUBLIC CURRENT_WAM INTERFACE  
!      
      curdate_wam = jdat
      curhour_wam = hr_jdat 
      curyear_wam = jdat(1)
      curmonth_wam = jdat(2)
      curday_wam =  jdat(3)
      JDAY_ycur0101=iw3jdn_VY(Jdat(1), 1, 1) 

      curddd_wam =  JDAY - JDAY_ycur0101 +1
      curjulday_wam  = JDAY
      curutsec_wam = nint(hr_jdat*3600.) 
!
!     curnsteps_wam 
!     curdtime_wam 
!
      END SUBROUTINE CURRENT_NCEP_JDAT

      SUBROUTINE wam_split_ymd(ymd, hh, jdat, ndi)
!
! purpose: transform (ymd, hh) => jdat[4]: yy, mm, dd, hh
!
      integer :: ymd
      integer :: hh
      integer :: ndi
      integer :: mmdd
      integer :: jdat(ndi)
      jdat(4) = hh
      jdat(1) = (ymd - mod(ymd, 10000))/10000
      mmdd = ymd -Jdat(1)*10000
      jdat(2) = (mmdd - mod(mmdd, 100))/100
      jdat(3) =  mmdd -jdat(2)*100
      END subroutine wam_split_ymd
!
      SUBROUTINE wam_calendar(ndi, jdat1,hr1, jdat2, hr2, jdatc, hrc)  
!

      implicit NONE
      integer :: ndi
      integer, dimension(ndi) :: jdat1, jdat2, jdatc
      real :: hr1, hr2, hrc
      real :: w1, w2
      call weights_time_interp(ndi,jdat1,hr1,jdat2,hr2,jdatc,hrc,w1,w2)

      write(iulog, *)  jdat1(1),jdat1(2),jdat1(3), ' ymd_1 '
      write(iulog, *)  jdatc(1),jdatc(2),jdatc(3), ' ymd_c '
      write(iulog, *)  jdat2(1),jdat2(2),jdat2(3), ' ymd_2 '   
       write(iulog, *)  w1, ' weight w1 for ymd_1 '
       write(iulog, *)  w2, ' weight w2 for ymd_2 '

      end subroutine wam_calendar
!
      SUBROUTINE weights_time_interp(ndi, Jdat1, hr1, Jdat2, hr2, Jdatc,
     &           hrc, w1, w2)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Current & Input data "DATE" representations
!   idate_ymd(3): (yyyy, mm, dd) + fraction of day = ncsec/86400.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      use idea_io_units, only : iulog
!SK   use idea_mpi_def, only  : mpi_WAM_quit
      implicit NONE
      integer :: ndi
      integer, dimension(ndi) :: Jdat1, Jdat2     ! two consecutive data records in YMDH-format
      integer, dimension(ndi) :: Jdatc            ! current YMDH
      real                    :: hr1, hr2         ! real hours of data with digits sy 2.55 hours
      real                    ::  hrc             ! real hours of model for given day
      real, intent(out) :: w1, w2                 ! weights for interpolation in time for (jdatc, hrc)
      real :: ydh1
      real :: ydh2
      real :: ydhc
      integer :: jd1,jd2, jdc
      integer, parameter :: iret = 23901 
!
      real, parameter :: f24 =1./24.
!
      JD1 = iw3jdn_vy(Jdat1(1),Jdat1(2),Jdat1(3))    ! compute Julian day from YYYYMMDD-format
      JD2 = iw3jdn_vy(Jdat2(1),Jdat2(2),Jdat2(3))
      JDc = iw3jdn_vy(Jdatc(1),Jdatc(2),Jdatc(3))
!
! add fractional day Hr/24. to the Julian Day
!
      ydh1 = float(jd1)+ hr1*f24
      ydh2 = float(jd2)+ hr2*f24
      ydhc = float(jdc)+ hrc*f24
         
!
!  Year does not matter for JULIAN date-line format
!
      if (ydhc.gt.ydh2.or.ydhc.lt.ydh1) then
      write(iulog, *) ' check yd1-ydc-yd2 in idea_cal_advance.f '
      write(iulog, *)  Jdat1(1),Jdat1(2),Jdat1(3), ' ymd_1 '
      write(iulog, *)  Jdatc(1),Jdatc(2),Jdatc(3), ' ymd_c '
      write(iulog, *)  Jdat2(1),Jdat2(2),Jdat2(3), ' ymd_2 '
      write(iulog, *) hr1, hrc, hr2,  ' hr1 < hrc < hr2  '
      write(iulog, *) jd1, jdc, jd2,  ' Julians-VAY 1-C-2 ' 
!SK   CALL mpi_WAM_quit(iret,'weights_time_interp in <idea_cal_adv.f>' )
      endif
      if(ydh2.ne.ydh1) w1 = (ydh2-ydhc)/(ydh2-ydh1)
      if(ydh2.eq.ydh1) w1 =0.5
      w2 =1.-w1 
!
      end  subroutine weights_time_interp
!
      subroutine wam_julday_doy(idat, jdoy,jday)
!
! it RETURNS two INTEGERS , THE DAY OF YEAR, AND JULIAN DAY
!                            from IDAT(4)
! Here we use only YMDH-format 
!         adding also the real HOUR with "digits" after "." like 0.55 or 23.55
!
!      IDAT(4) rather than  INTEGER (8) of the NCEP ABSOLUTE DATE AND TIME
!      (YEAR, MONTH, DAY, "TIME ZONE", HOUR,  MINUTE, SECOND, MILLISECOND)
!
        integer, intent(in)  :: idat(4) !(YEAR, MONTH, DAY, HOUR)
        integer, intent(out) :: jday    !      JDAY       INTEGER JULIAN DAY (DAY NUMBER FROM JAN. 1,4713 B.C.)
        integer, intent(out) :: jdoy    !      JDOY       INTEGER DAY OF YEAR (1-366, WHERE 1 IS JANUARY 1)
!local
        integer :: jdow  !      JDOW       INTEGER DAY OF WEEK (1-7, WHERE 1 IS SUNDAY)
!JDAY
      JDAY = iw3jdn_VY(idat(1),idat(2),idat(3))

!      JDOY       INTEGER DAY OF YEAR (1-366, WHERE 1 IS JANUARY 1)
      CALL w3fs26_VY(jday, idat(1),idat(2),idat(3),jdow,jdoy)

      END subroutine wam_julday_doy
!
!
      FUNCTION IW3JDN_VY(IYEAR,MONTH,IDAY)
!
!     compute Julian day from YYYYMMDD-format
!
!     IYEAR  ARG LIST  INTEGER   YEAR           ( 4 DIGITS)
!     MONTH  ARG LIST  INTEGER   MONTH OF YEAR   (1 - 12)
!     IDAY   ARG LIST  INTEGER   DAY OF MONTH    (1 - 31)
  
      INTEGER :: IW3JDN_VY 
!      INTEGER :: YYYYMMDD
!locals
      INTEGER :: IYEAR,MONTH,IDAY
       IW3JDN_VY  =    IDAY - 32075
     &            + 1461 * (IYEAR + 4800 + (MONTH - 14) / 12) / 4
     &            + 367 * (MONTH - 2 - (MONTH -14) / 12 * 12) / 12
     &            - 3 * ((IYEAR + 4900 + (MONTH - 14) / 12) / 100) / 4
       RETURN
       END FUNCTION IW3JDN_VY
!
       SUBROUTINE W3FS26_VY(JLDAYN,IYEAR,MONTH,IDAY,IDAYWK,IDAYYR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! compute from Julian day: (YYYY, Month, Day_of_week(?) and Day of year
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       integer :: JLDAYN
       integer :: IYEAR
       integer :: MONTH
       integer :: IDAY
       integer :: IDAYWK
       integer :: IDAYYR
!$$$   SUBPROGRAM  DOCUMENTATION  BLOCK
!
! SUBPROGRAM: W3FS26         YEAR, MONTH, DAY FROM JULIAN DAY NUMBER
!     IYEAR  ARG LIST  INTEGER   YEAR  (4 DIGITS)
!     MONTH  ARG LIST  INTEGER   MONTH
!     IDAY   ARG LIST  INTEGER   DAY
!     IDAYWK ARG LIST  INTEGER   DAY OF WEEK (1 IS SUNDAY, 7 IS SAT)
!     IDAYYR ARG LIST  INTEGER   DAY OF YEAR (1 TO 366)
       integer :: L, N, I, J

       L      = JLDAYN + 68569
       N      = 4 * L / 146097
       L      = L - (146097 * N + 3) / 4
       I      = 4000 * (L + 1) / 1461001
       L      = L - 1461 * I / 4 + 31
       J      = 80 * L / 2447
       IDAY   = L - 2447 * J / 80
       L      = J / 11
       MONTH  = J + 2 - 12 * L
       IYEAR  = 100 * (N - 49) + I + L
       IDAYWK = MOD((JLDAYN + 1),7) + 1
!
       IDAYYR = JLDAYN -
     &  (-31739 +1461 * (IYEAR+4799) / 4 - 3 * ((IYEAR+4899)/100)/4)
!
       RETURN
       END SUBROUTINE W3FS26_VY
!

!
      SUBROUTINE JULDAY_2_YMD(JLDAYN,IYEAR,MONTH,IDAY)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! compute from Julian day: (YYYY, Month, Day_of_week(?) and Day of year
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       integer :: JLDAYN
       integer :: IYEAR
       integer :: MONTH
       integer :: IDAY
!
! SUBPROGRAM: JULDAY_2_YMD FROM JULIAN DAY NUMBER =>YEAR, MONTH, DAY
!     IYEAR  ARG LIST  INTEGER   YEAR  (4 DIGITS)
!     MONTH  ARG LIST  INTEGER   MONTH
!     IDAY   ARG LIST  INTEGER   DAY
!     IDAYWK ARG LIST  INTEGER   DAY OF WEEK (1 IS SUNDAY, 7 IS SAT)
!     IDAYYR ARG LIST  INTEGER   DAY OF YEAR (1 TO 366)
       integer :: L, N, I, J

       L      = JLDAYN + 68569
       N      = 4 * L / 146097
       L      = L - (146097 * N + 3) / 4
       I      = 4000 * (L + 1) / 1461001
       L      = L - 1461 * I / 4 + 31
       J      = 80 * L / 2447
       IDAY   = L - 2447 * J / 80
       L      = J / 11
       MONTH  = J + 2 - 12 * L
       IYEAR  = 100 * (N - 49) + I + L
!
!       IDAYYR = JLDAYN -
!     &  (-31739 +1461 * (IYEAR+4799) / 4 - 3 * ((IYEAR+4899)/100)/4)
!
       RETURN
       END SUBROUTINE JULDAY_2_YMD
!
! below test-package
!
      subroutine test_wam_calendar
      integer, parameter :: ndi=4
      integer, dimension(ndi) :: jdat, idat
      character(len=4), dimension(ndi) :: Cjdat
      integer :: ymd, hh
      integer :: k
      real :: hr_a_idat, hr_jdat
      Cjdat(1) = 'year'
      Cjdat(2) = 'mont'
      Cjdat(3) = ' day'
      Cjdat(4) = 'hour'
      hh  = 0
      ymd = 20131231
      call wam_split_ymd(ymd, hh, Idat, ndi)
      do k=1, ndi
      write(iulog,*)  idat(k), ' ', Cjdat(k) 
      enddo
      hr_a_idat = 366.*24.+ 0.5
      CALL CURRENT_NCEP_JDAT(idat,hr_a_idat,jdat, hr_jdat)
      do k=1, ndi
       write(iulog,*) jdat(k), ' ', Cjdat(k) 
      enddo 

           write(iulog,*)  hr_jdat, ' hr_jdat  0.5 '

      end subroutine test_wam_calendar
!
!========================================================
! Brief Summary on the NCEP/GFS-NEMS calendar
!                      w3lib - std-library
! NCEP-calendar INTEGER (8)
!
!      IDAT   INTEGER (8) NCEP ABSOLUTE DATE AND TIME
!             (YEAR, MONTH, DAY, TIME ZONE,
!              HOUR, MINUTE, SECOND, MILLISECOND)
! DO simple wam_calendar subs yyddd.fraction
!                             yymmdd (20010101)
!     handling also leap-year calendars
!     returning seconds in the curr-day
! grep -i idate ../gsm/phys/gfs_physics_start_time_get_mod.f
!      use date_def,     only: idate,idate7
!
!        idate  = head%idate  (reverse order from file-header)
!        yy     = idate(4)
!        mm     = idate(2)
!        dd     = idate(3)
!        hh     = idate(1)
!
!        compute weihgts for the interpolation in time
!
! gbphys.f:       real(kind=kind_phys), intent(in) ::  dtp,  dtf, fhour, solhr
! date_def.f:     real(kind=kind_evod),target      ::  fhour
!=======================================================================
!   USE of calendar in the IDEA-physics
! based on "date_def.f" :      integer,target :: idate(4),idate7(7)
!           idea_ion.f:      use date_def
! idea_solar_heating.f:      use date_def
!=======================================================================
!idea_ion.f:        imo=idate(2)
!idea_ion.f:        iday_m=idate(3)
!idea_solar_heating.f:      INTEGER   i,idat(8),jdat(8),jdow,jday
!idea_solar_heating.f:      idat=0
!idea_solar_heating.f:      idat(1)=idate(4) YY
!idea_solar_heating.f:      idat(2)=idate(2) MM
!idea_solar_heating.f:      idat(3)=idate(3) DD
!idea_solar_heating.f:      idat(5)=idate(1) HH
!idea_solar_heating.f:      call w3movdat(rinc,idat,jdat)
      END MODULE wam_date_calendar     

