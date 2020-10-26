      MODULE IDEA_IMF_INPUT
!===============================================================================
! VAY Oct 2015  Implementation with T.Matsuo's wam_nems_imf_2012.nc
!
! file_imf =/scratch3/NCEPDEV/swpc/sair/data/input/20120101/wam_nems_imf_2012.nc
!
!	double ydh(time) ;   
!       2012001.00347222, 2012001.00694444......2012001.99652778
!       FORMAT:    YYYY_DDD . FRACT_of_day
!================================================================================
       use IDEA_IO_UNITS, only : iulog
       implicit none
       public :: imf_wam_get_4fields
       public :: imf_read_wam_init
       public :: imf_wamstep_advance
       public :: dealloc_imf
!
       save 
! 
       integer, parameter  :: itheia = 1     !vay-aug THEIA no-mpi broadcast 
!
       integer :: idea_imf_fix      ! flag for imf-variable input (0) or fixed (1) in nml idea_ion

       integer ::  ntimes_imf
       real, allocatable :: times_imf(:)

      real, allocatable :: ydh_imf(:)
      real, allocatable :: BZ_imf(:)
      real, allocatable :: BY_imf(:)
      real, allocatable :: Swden(:)
      real, allocatable :: Swvel(:)
!
       integer, parameter      ::  ndi = 4 ! idat and jdat -array dimensions
       integer, dimension(ndi) :: idatc    ! initial idat[ymdh]

       integer, dimension(ndi) :: jdat1    ! current closest data ymdh_d < ymdh
       integer, dimension(ndi) :: jdatc    ! current model ymdh
       integer, dimension(ndi) :: jdat2    ! current closest data ymdh_d > ymdh

       real :: hr1, hr2, hrc               ! corresponding real hours
!
! Values for the current model step of WAM-physics
!
       integer               :: tim_ndx1_imf
       integer               :: tim_ndx2_imf
       real :: w_ndx1_imf 
       real :: w_ndx2_imf
!
       real :: Bz_s, By_s, Swden_s, Swvel_s           ! input => "idea_ion"
       real :: Bz_fix, By_fix, Swden_fix, Swvel_fix   ! fixed values => "idea_ion"
!
       CONTAINS
!
       subroutine dealloc_imf(mpi_id)
       integer :: mpi_id

       deallocate(ydh_imf, times_imf, Bz_imf, By_imf, Swden, Swvel)
 
       if (mpi_id.eq.0) then
       write(iulog, *)  'subroutine dealloc_imf: free memory '
       endif

       end subroutine dealloc_imf
!       
       subroutine imf_wam_get_4fields  !( f107_s, f107a_s, ap_s, kp_s)
       implicit none
!local

        Bz_s  =  Bz_imf(tim_ndx1_imf)*w_ndx1_imf
     &        + Bz_imf(tim_ndx2_imf)*w_ndx2_imf
        By_s  =  By_imf(tim_ndx1_imf)*w_ndx1_imf
     &        + By_imf(tim_ndx2_imf)*w_ndx2_imf
!
        Swden_s = Swden(tim_ndx1_imf)*w_ndx1_imf
     &          + Swden(tim_ndx2_imf)*w_ndx2_imf
        Swvel_s = Swvel(tim_ndx1_imf)*w_ndx1_imf
     &          + Swvel(tim_ndx2_imf)*w_ndx2_imf
!
       end subroutine imf_wam_get_4fields  !
!
       subroutine imf_wamstep_advance(mpi_id, Mjdat_cur, hour_cur)

!      update calendar find "tim_ndx" and compute "weights-interp"
!      compute  "f107_s, f107a_s, ap_s, kp_s, euv_s(37)"
!      values for 'the current step of WAM"
!      times_imf(1:NTIMES_imf) - fixed "data-calendar" like "20130101.75"
!
       use wam_date_calendar, only : weights_time_interp
       integer ::  mpi_id
       integer ::  Mjdat_cur(ndi)
       real    ::  Hour_cur      
!
! jdatc ( yr, mon, day, day_fraction) 
!      define
!     [Jdat1, hr1, Jdat2, hr2, Jdatc, hrc]
!

       Jdatc = Mjdat_cur
       hrc  = hour_cur
       call IMF_update_jdates(mpi_id)     ! for given year "IMF_update_jdates" computes w_ndx1_imf, w_ndx2_imf
      

!     after update_jdates anf hours => compute weights
!     General for Y-2-Y transition ..... [Jdat1, hr1, Jdat2, hr2, Jdatc, hrc]
!
!     call imf_weights_time_interp(ndi, Jdat1, hr1, Jdat2, hr2, Jdatc, hrc, w_ndx1_imf, w_ndx2_imf)
!
!       
       if (idea_imf_fix.eq.0) CALL imf_wam_get_4fields ! Bz_s, By_s, Swden_s, SwVel_s
 
!
!       if (mpi_id.eq.0) then
!       write(iulog,*) ' imf_adv IND + Weights: ', tim_ndx1_imf, tim_ndx2_imf, w_ndx1_imf, w_ndx2_imf
!       write(iulog,*) ' Bz_s,  By_s', Bz_s, By_s
!       endif
!
       end subroutine imf_wamstep_advance
!
       SUBROUTINE IMF_update_jdates(mpi_id)
!
! find   [idat1, hr1, tim_ndx1], [idat2, hr2, tim_ndx2] for the start day
! update [idat1, hr1, tim_ndx1], [idat2, hr2, tim_ndx2] for
!       call weights_time_interp(ndi, idat1, hr1, idat2, hr2, idatc, hrc, w_ndx1, w_ndx2)
!
!        integer :: idatc(4)
        use wam_date_calendar, only : wam_split_ymd, curddd_wam
        use module_physics_driver, only: is_master
        implicit none
        integer :: mpi_id
        integer :: wrk_date_ymd,  wrk_date_yddd 
        real    :: wrk_time
        integer :: n, nk
        integer :: iymd1, iymd2, iyd1, iyd2, im1, im2, id1, id2
        integer :: iy1, iy2
        integer :: ihh1, ihh2, ihhc
        ihhc = Jdatc(4)
        wrk_date_ymd = 10000*Jdatc(1) + 100*Jdatc(2) + Jdatc(3) 
!
!       CALL wam_julday_doy(Jdatc, YDDD, JUL_day)
!
        wrk_date_yddd = 1000*Jdatc(1)+curddd_wam
        wrk_time = float(wrk_date_yddd) + hrc/24.                !wrk_time = flt_date( wrk_date, 0 )

       if (is_master) then     !SK 2020      
!SK    if (mpi_id.eq.0) then      
        write(iulog,*) 
     & ' vay wrk_time:',wrk_time, hrc, ' hour ', wrk_date_yddd, ' yddd '

        if( wrk_time < times_imf(1) .or. 
     &      wrk_time > times_imf(ntimes_imf) ) then

         write(iulog,*) 'imf_files: model time is out of-range Bz-By'
         write(iulog,*)   times_imf(1) ,   times_imf(ntimes_imf) ,
     &    ' times_imf(start -/- end '
         write(iulog,*) wrk_time, wrk_date_yddd, 
     &   ' wrk_time, wrk_ddd, wrk_ymd ', wrk_date_ymd
! 
        end if
       ENDIF

! time is growing  FIND INDEX
        nk = 2
        do n = 2,ntimes_imf  
           if( wrk_time <= times_imf(n) ) then
           nk = n
           exit
          end if
       end do
!
       tim_ndx1_imf = nk-1
       tim_ndx2_imf = nk
!
       w_ndx2_imf=(wrk_time-times_imf(nk-1))/
     &            (times_imf(nk)-times_imf(nk-1))
       w_ndx1_imf = 1. -w_ndx2_imf 

       RETURN          ! the rest is "under development"
!
!  the rest is "under development" for handling Year1-to-Year2 transition
!
       iyd1 = int( times_imf(nk-1) )         !dates(tim_ndx1)
       iyd2 = int( times_imf(nk-1) )         !dates(tim_ndx2)

       hr1 = (times_imf(nk-1)-float(iyd1))*24.
       hr2 = (times_imf(nk)  -float(iyd2))*24.
       ihh1 = nint(hr1)
       ihh2 = nint(hr2)  
!
!       call yddd_ymd(iyd1, iy1, im1, id1, iymd1)
!       call yddd_ymd(iyd2, iy2, im2, id2, iymd2)
!
       iy1 = int( times_imf(nk-1) )/1000          !(iyd1 - mod(iyd1, 1000))/1000 
       iy2 = int( times_imf(nk  ) )/1000  
       im1 = 1
       im2 = 1
       id1 = 1
       id2 = 1
       iymd1 = iy1*10000+im1*100+id1                       !dates(tim_ndx1)
       iymd2 = iy2*10000+im2*100+id2                       !dates(tim_ndx2)

       if (hr2.lt.hrc.and.(iymd1 == iymd2)) hr2 = hrc
!
! lines below should be activated when we have TRANSITION from Year1 => Year2
!       and  Julian days Jdat1(4) & Jdat2(4) can be used
!       for  now we will use "growing" times-array .....
!    
! fill in Jdat(1:4) Y-M-D-H
!
       call wam_split_ymd(iymd1, ihh1, jdat1, ndi)
       call wam_split_ymd(iymd2, ihh2, jdat2, ndi)
!
       end subroutine IMF_update_Jdates
!
       subroutine IMF_read_wam_init(file , mpi_id)   ! take it out from  "idea_solar_input"
!
! called from idea_solar_heating.f:CALL imf_read_wam_init(ncIMF_fpath, mpi_id)   
!
       use netcdf      
       use module_physics_driver, only: is_master
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

        if(is_master) then      !SK 2020
!SK     if(mpi_id.eq.0) then
           write(iulog,*)file        
           write(iulog,*) 'IMF_PARMS: opening file ', trim(file) 
        endif
!--------------------------------------------------------------------------------------
!dimensions:	time = 105410 ;
!	double by(time) ;
!		by:units = "nT" ;
!		by:long_name = "1AU IP By, GSM" ;
!	double bz(time) ;
!		bz:units = "nT" ;
!		bz:long_name = "1AU IP Bz, GSM" ;
!	double swden(time) ;
!		swden:units = "per cc" ;
!		swden:long_name = "1AU IP N (ion)" ;
!	double swvel(time) ;
!		swvel:units = "Km/s" ;
!		swvel:long_name = "1AU IP Plasma Speed" ;
!	double ydh(time) ;
!		ydh:long_name = "year-day plus fractional day: yyyyddd.frac" ;
!---------------------------------------------------------------------------------------
       ierNC=NF90_OPEN(trim(File), nf90_nowrite, ncid)   
       if (iernc /=0) write(iulog,*) ncid, 'ncid ', iernc, ' iernc '
!
         iernc=nf90_inq_varid( ncid, 'ydh', vid )
         if (iernc /=0) write(iulog,*) ncid, 'ncid ', iernc, ' iernc '
         ierNC=nf90_inquire_variable(ncid, vid, dimids=dimidT) 
         iernc = nf90_inquire_dimension(ncid, dimidT(1), len=ntimes_imf)
         

         if(is_master) then      !SK 2020
!SK      if(mpi_id.eq.0) then
              write(iulog,*) ntimes_imf, ' nt-nw  idea_imf_input'
         endif
         
 !   
       allocate( times_imf(ntimes_imf),stat=astat )  
 ! 
       if( astat /= 0 ) then
       write(iulog,*) ' alloc_err in read_imf for ntimes', ntimes_imf 
       end if     

        iernc=nf90_inq_varid( ncid, 'ydh', vid )
        iernc= nf90_get_var( ncid, vid, times_imf)

          if(is_master) then      !SK 2020
!SK       if(mpi_id.eq.0) then
            write(iulog,*) times_imf(1), times_imf(ntimes_imf), 
     &                    ' ydh-imf ' 
          endif

!        do n = 1,ntimes
!           dfhours_imf(n) = 0.0          ! current for daily  12UT
!           times_imf(n) = float(dates(n))  + dfhours(n)/24.    
!        end do
!
! init hr1 & hr2
!          hr1 = dfhours(1)
!          hr2 = dfhours(2)
    !---------------------------------------------------------------
    !	... allocate and read solar parms ..... ALL-time series
    !   call dealloc_solar(mpi_id) in the End of WAM-RUN
    !   we do not put these data-sets in the restart files
    !---------------------------------------------------------------
       allocate(  BY_imf(ntimes_imf), BZ_imf(ntimes_imf),stat=astat )
       allocate(  SWDEN(ntimes_imf),  SWVEL(ntimes_imf), stat=astat )
      
       if( astat /= 0 ) then
         write(iulog,*) ' alloc_err( astat, BY ...SWVEL ', ntimes_imf
       end if

        iernc=nf90_inq_varid( ncid, 'by', vid )
        iernc= nf90_get_var( ncid, vid, by_imf)

        iernc=nf90_inq_varid( ncid, 'bz', vid )
        iernc= nf90_get_var( ncid, vid, bz_imf)

        iernc=nf90_inq_varid( ncid, 'swden', vid )
        iernc= nf90_get_var( ncid, vid, SWDEN)


        iernc=nf90_inq_varid( ncid, 'swvel', vid )
        iernc= nf90_get_var( ncid, vid, SWVEL)

        iernc=nf90_close(ncid)     
!
     
!
          if(is_master) then      !SK 2020
!SK       if(mpi_id.eq.0) then
          write(iulog,*) '  read_wam_IMF: ntimes  ', ntimes_imf   
          write(iulog,*) maxval(bz_imf),   minval(bz_imf), ' BZ_imf '
          write(iulog,*) maxval(bz_imf),   minval(bz_imf), ' BY_IMF ' 
        

          write(iulog,*) maxval(swden),   minval(swden), ' SW-density ' 
          write(iulog,*) maxval(swvel),   minval(swvel), ' SW-velocity'      
          write(iulog,*)            ' mpi_bcast in IMF_read_wam_init'
          write(iulog,*)  ' completed IMF_read_wam_init'
         endif

        RETURN    ! Here RETURN is a temporary FIX of mpif.h MPI_REAL8/mpi_integer for THEIA
                  ! ALL PEs read nc-file
!
! result of solar_read_wam_init(ncfile_fpath): ALLOCATE, READ and BROADCAST
! call MPI_BCAST (data_to_be_sent, send_count, send_type, broadcasting_process_ID, comm, ierr)

!SK    call mpi_bcast(ntimes_imf,1,mpi_integer,0,mpi_comm_all,info)

!       call mpi_bcast(dates,ntimes,mpi_integer,0,mpi_comm_all,info)
!       call mpi_bcast(dfhours,ntimes,MPI_REAL8,0,MPI_COMM_ALL,info)
!SK    call mpi_bcast(times_imf,ntimes_imf,
!SK  &                MPI_REAL8,0,MPI_COMM_ALL,info)

!SK    call mpi_bcast(bz_imf  ,ntimes_imf,MPI_REAL8,0,MPI_COMM_ALL,info)
!SK    call mpi_bcast(by_imf  ,ntimes_imf,MPI_REAL8,0,MPI_COMM_ALL,info)

!SK    call mpi_bcast(SWDEN  ,ntimes_imf,MPI_REAL8,0,MPI_COMM_ALL,info)
!SK    call mpi_bcast(SWVEL  ,ntimes_imf,MPI_REAL8,0,MPI_COMM_ALL,info)


!       call mpi_barrier(mpi_comm_all,info)         
!SK    if(mpi_id.eq.0) then
!SK       write(iulog,*)  ' completed IMF_read_wam_init'
!SK    endif
!
!
       end  subroutine IMF_read_wam_init
!
!
!
       END MODULE IDEA_IMF_INPUT
