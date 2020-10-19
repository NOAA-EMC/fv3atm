!=======================================================================
! vay Oct 31 2015 first implementation for
!                 initial fixed wam-ion input fields
!                 transfer input files & data in the 
!                 netcdf files to avoid "danger" of program. bugs 
!                 by 1000-es lines of "data"-operators
! /scratch3/NCEPDEV/swpc/save/Tzu-Wei.Fang/Tiros:  ionprof  tiros_spectra
! nc-files see namelist :  ion_in
!========================================================================
      Module IDEA_ION_INPUT
!=====================================================================
!idat1 real cormag(20,91),btot(20,91),dipang(20,91),glat(91),glon(20)
!
      use IDEA_IO_UNITS, only : iulog
!
      integer, parameter :: tiros_switch =1  ! 0 -ncfile; 1-ascii
      integer, parameter :: NYMAG = 91       ! dimension of glat
      integer, parameter :: NXMAG = 20       ! dimension of glon
      real               :: glat(NYMAG)
      real               :: glon(NXMAG)
      real, dimension(NXMAG, NYMAG) ::    cormag, btot, dipang
!=====================================================================
      real  :: emaps1(21,20,7),cmaps1(21,20,7), djspectra1(15,21)
!
      integer, parameter :: NT_21 = 21   ! TIROS-21  1st dimension
      integer, parameter :: NT_20 = 20   ! TIROS-20  2nd dimension
      integer, parameter :: NT_7  =  7   ! TIROS-7   3rd dimension
!
      real, dimension(NT_21, NT_20, NT_7) :: EMAPS, CMAPS     ! TIROS

      integer, parameter :: N_FLX    = 15               ! TIROS
      integer, parameter :: N_BND    = 21               ! TIROS
      real               :: djspectra(N_FLX, N_BND)     ! TIROS
!=====================================================================
      integer, parameter :: jmaxwell = 6
      real, parameter    :: width_maxwell = 0.050
      real               :: en_maxwell(jmaxwell)
!
      real    :: GW_fixnam                   ! the number of gigawatts of auroral power
      integer :: tiros_activity_fixnam       ! tiros_activity_level = 7 in idea_ion.f     
!
!  
!     real, parameter    ::      pi = 3.145926   
      real, parameter    ::      bz = 1.38e-23
      real, parameter    ::      gscon = 8.314e3
      real, parameter    ::      mo = 16.
      real, parameter    ::      mo2 = 32.
      real, parameter    ::      mn2 = 28.
      real, parameter    ::      mh = 1.
      real, parameter    ::      mhe = 4.
      real, parameter    ::      amu = 1.661e-27
      real, parameter    ::      E0=0.035
!23456
      real :: width(N_FLX), en(N_FLX)            ! defined by data-stat
      real :: TE11(N_BND),TE15(N_BND)            !precomp_iondata_fixed
!
       real ::  RATIO(21)                        !precomp_iondata_fixed
       real ::  ionchr(21), RLAM(21)  
       real ::  ion_recomb(8),lognpres(8)        ! defined by data-stat
!
      data en/.37,.6,.92,1.37,2.01,2.91,4.19,6.,8.56,12.18,
     &       17.3,24.49,36.66,54.77,81.82/
      data width/.158,.315,.315,.63,.631,1.261,1.26,2.522,
     &       2.522,5.043,5.043,10.,14.81,22.13,33.06/



      data RLAM/1.49,1.52,1.51,1.48,1.43,1.37,1.30,1.22,
     &      1.12,1.01,0.895,0.785,0.650,0.540,0.415,0.320,0.225,
     &      0.14,0.08,0.04,0.0/
!
      DATA ionchr/.378 , .458 , .616 , .773 , .913 , 1.088 , 1.403 ,
     &     1.718 , 2.033 , 2.349 , 2.979 , 3.610 , 4.250 , 4.780 ,
     &     6.130 , 7.392 , 8.653 , 9.914 , 12.436 , 14.957 , 17.479/
!
       
      data lognpres/-3.425,-4.835,-5.918,-7.066,-7.784,-8.366,-9.314,
     &     -10.507/

      data ion_recomb/3.20e-13,3.20e-13,2.75e-13,1.45e-13,1.13e-13,
     &      8.30e-14,3.70e-14,2.00e-14/


!      CHARACTER (LEN=*), parameter :: filename3=
!     & '/scratch3/NCEPDEV/swpc/save/Tzu-Wei.Fang/Tiros/ionprof'
!       CHARACTER (LEN=*), parameter :: filename7=
!     & '/scratch3/NCEPDEV/swpc/save/Tzu-Wei.Fang/Tiros/tiros_spectra'
!
      CONTAINS
!
      subroutine precomp_iondata_fixed
!
!
! precompute 4-fixed: ratio(m=1:21
!                    te15(m=1:21) te11(m=1:21) en_maxwell(1:6)
!  
      integer :: m, iband, j
      do  m=1,21
       ratio(m) = (m-1)*0.05
      enddo
!
      do iband=1,21
      te15(iband)=0.0
      te11(iband)=0.0
! the ionization rates will need to be normalized to TE11, which is
! the energy flux between 300eV and 20keV, which is provided by the
! TIROS energy influx maps emaps, rather than the energy from
! 300eV to 100keV, which is what the spectra were normalized to
!
! check the energy influx is normalized to 1 erg/cm2/s
      do m=1,15
       TE15(iband)=TE15(iband)+djspectra(m,iband)*en(m)*width(m)*1.6E-06
      enddo                !normalize with the energy influx 300eV to 20keV
      do  m=1,11
      TE11(iband)=TE11(iband)+djspectra(m,iband)*en(m)*width(m)*1.6E-06
      enddo 
                           !print *, iband, TE11(iband), TE15(iband)
      enddo ! iband=1,21
!   
      do j = 1,jmaxwell
         en_maxwell(j) = j*0.05 - 0.025
      enddo    
!
      end subroutine precomp_iondata_fixed
!       
      subroutine ion_read_namelist(nml_ion, nlun_ion, 
     & ncfile_fpath, nctiros_fpath, ncimf_fpath, mpi_id)
!
! read name-list
!
!      USE IDEA_IMF_INPUT, only : idea_imf_fix    ! data-based Bz-By-Swden Swvel
!
      integer                       :: idea_imf_fix=1
!
      integer ,intent(in)           :: mpi_id
      character(len=*),intent(in)   :: NmL_ion
      integer, intent(in)           :: nlun_ion
      character(len=*), intent(out) :: ncfile_fpath
      character(len=*), intent(out) :: nctiros_fpath
      character(len=*), intent(out) :: ncimf_fpath

      integer, parameter :: ch100 = 256
      integer, parameter :: me =0
      integer, parameter :: masterproc =0


      INTEGER :: k, i,j
      integer :: unitn, ierr     
      character(ch100) :: Dir_swpc       !='/scratch1/portfolios/NCEPDEV/swpc/save/'
      character(ch100) :: Dir_uid        !='Valery.Yudin/NEMS/wam_april_2015/data_euv/'
      character(ch100) :: Dir_ion        !=Dir_swpc//Dir_uid
      character(ch100) :: Dir_imf        !=Dir_swpc//Dir_uid
      character(ch100) :: dation_file    ! name of file w/o the full-path
      character(ch100) :: tiros_file     ! name of file w/o the full-path
      character(ch100) :: imf_file       ! name of file w/o the full-path for IMF
!
!
      real    :: GW_fix                  ! the number of gigawatts of auroral power
      integer ::tiros_activity_fix       ! tiros_activity_level = 7 in idea_ion.f
!
      namelist /ion_parms_nl/ GW_fix, tiros_activity_fix, idea_imf_fix,
     &  dation_file, tiros_file, imf_file, Dir_swpc, Dir_uid, Dir_imf
!=========================================================================
!
! "ion_in"-namleist file
!      should be copied to $RUNDIR along with other namelists
!
! module wam_phys_control
!-----------------------------------------------------------------------
! will
! Provide a control interface to WAM physics packages
!
! ....VDIFF(mol+eddy)-GWs-SOLAR-RAD/NLTE-CHEM/NEUT-TRACERS
! 
! read "nml_ion" by defauly "ion_in"
!=========================================================================
      open(nlun_ion, file=trim(nml_ion), status='old' )
      read(nlun_ion, ion_parms_nl, iostat=ierr)   
      close(nlun_ion)    
!
!
!
       Dir_ion=trim(Dir_swpc)//trim(Dir_uid)
      
      ncfile_fpath= trim(Dir_ion)//trim(dation_file) 
      nctiros_fpath= trim(Dir_ion)//trim(tiros_file) 
      ncimf_fpath= trim(Dir_imf)//trim(imf_file) 

      if (mpi_id.eq.0) then     
      write(iulog,*) GW_fix, 'idea_ion_input GW_fix '
      write(iulog,*) tiros_activity_fix, ' tiros_activity_level   '


      write(iulog,*)    '+++++++++++ ncfile_fpath '
      write(iulog,133) ncfile_fpath
      write(iulog,*)    '+++++++++++  '
!                                                       write(iulog,*) idea_IMF_fix, 'idea_IMF_fix - flag (1-fixed) '
      write(iulog,*)    '+++++++++++ nctiros_fpath '
      write(iulog,133)  nctiros_fpath
      write(iulog,*)    '+++++++++++  '

      write(iulog,*)    '+++++++++++ ncimf_fpath '
      write(iulog,133)   ncimf_fpath
      write(iulog,*)    '+++++++++++ '
      write(iulog,*) idea_imf_fix, 'idea_ion_input ... idea_imf_fix ' 
      endif
133   format(A122)
!                                  
      GW_fixnam = GW_fix
      tiros_activity_fixnam =tiros_activity_fix

      end subroutine ion_read_namelist
!
!
!
      subroutine tiros_read_wam_init(file , mpi_id)
!
! called from idea_ion.f:CALL tiros_read_wam_init(nctiros_fpath, mpi_id)   
!
       use netcdf      
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
             
!        integer  :: dimid    
!        integer, dimension(nf90_max_var_dims) :: dimidT         
 

        if(mpi_id.eq.0) then
           write(iulog,*)file        
           write(iulog,*) 'TIROS_READ...: opening file ', trim(file) 
        endif
!--------------------------------------------------------------------------------------
!         double emaps(nt_21, nt_20, nt_7) ;
!         double cmaps(nt_21, nt_20, nt_7) ;
!         double djspectra(n_flx, n_bnd) 
!---------------------------------------------------------------------------------------
       ierNC=NF90_OPEN(trim(File), nf90_nowrite, ncid)   
       if (iernc /=0) write(iulog,*) ncid, 'ncid ', iernc, ' iernc '
!
        iernc=nf90_inq_varid( ncid, 'emaps', vid )
        iernc= nf90_get_var( ncid, vid, emaps)
        iernc=nf90_inq_varid( ncid, 'cmaps', vid )
        iernc= nf90_get_var( ncid, vid, cmaps)
        iernc=nf90_inq_varid( ncid, 'djspectra', vid )
        iernc= nf90_get_var( ncid, vid, djspectra)
        iernc=nf90_close(ncid)    

!        iernc=nf90_inq_varid( ncid, 'ydh', vid )!
!         if (iernc /=0) write(iulog,*) ncid, 'ncid ', iernc, ' iernc '
!         ierNC=nf90_inquire_variable(ncid, vid, dimids=dimidT) 
!         iernc = nf90_inquire_dimension(ncid, dimidT(1), len=ntimes_imf)
!
 !      allocate( times_imf(ntimes),stat=astat )  
 !      allocate( dfhours_imf(ntimes),stat=astat )  
 !      if( astat /= 0 ) then
 !      write(iulog,*) ' alloc_err in read_waccm_solar for dates,times', ntimes 
 !      end if     

      if(mpi_id.eq.0) then
      write(iulog,*) '  tiros_read_wam_init '   
      write(iulog,*) maxval(emaps),   minval(emaps), ' EMAPS-tiros '
      write(iulog,*) maxval(cmaps),   minval(cmaps), ' CMAPS-tiros ' 
      write(iulog,*)  maxval(djspectra),   minval(djspectra), 
     &                    ' DJSPECT-tiros '      
          write(iulog,*)  ' completed tiros_read_wam_init'
      endif

      RETURN    ! Here RETURN is a temporary FIX of mpif.h MPI_REAL8/mpi_integer for THEIA
                  ! ALL PEs read nc-file on THEIA 
!
! result of solar_read_wam_init(ncfile_fpath): ALLOCATE, READ and BROADCAST
! call MPI_BCAST (data_to_be_sent, send_count, send_type, broadcasting_process_ID, comm, ierr)

!       call mpi_bcast(ntimes_imf,1,mpi_integer,0,mpi_comm_all,info)
!       call mpi_bcast(times_imf  ,ntimes_imf,MPI_REAL8,0,MPI_COMM_ALL,info)

!       call mpi_bcast(bz_imf  ,ntimes_imf,MPI_REAL8,0,MPI_COMM_ALL,info)
!       call mpi_bcast(by_imf  ,ntimes_imf,MPI_REAL8,0,MPI_COMM_ALL,info)

!       call mpi_bcast(SWDEN  ,ntimes_imf,MPI_REAL8,0,MPI_COMM_ALL,info)
!       call mpi_bcast(SWVEL  ,ntimes_imf,MPI_REAL8,0,MPI_COMM_ALL,info)
!       call mpi_barrier(mpi_comm_all,info)         
!      if(mpi_id.eq.0) then
!         write(iulog,*)  ' completed IMF_read_wam_init'
!      endif
!
!
       end  subroutine tiros_read_wam_init
!
!
!
       subroutine ion_read_wam_init(file , mpi_id)
!
! called from idea_ion.f:CALL tiros_read_wam_init(nctiros_fpath, mpi_id)   
!
       use netcdf      
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
             
!        integer  :: dimid    
!        integer, dimension(nf90_max_var_dims) :: dimidT         
 

        if(mpi_id.eq.0) then
           write(iulog,*)  file        
           write(iulog,*) 'TIROS_READ...: opening file ', trim(file) 
        endif
!--------------------------------------------------------------------------------------
!        double cormag(nxmag, nymag) ;
!        double btot(nxmag, nymag) ;
!        double dipang(nxmag, nymag) ;
!        double glon(nxmag) ;
!        double glat(nymag)
!---------------------------------------------------------------------------------------
       ierNC=NF90_OPEN(trim(File), nf90_nowrite, ncid)   
       if (iernc /=0) write(iulog,*) ncid, 'ncid ', iernc, ' iernc '
!
        iernc=nf90_inq_varid( ncid, 'glon', vid )
        iernc= nf90_get_var( ncid, vid, glon)
        iernc=nf90_inq_varid( ncid, 'glat', vid )
        iernc= nf90_get_var( ncid, vid, glat)
        iernc=nf90_inq_varid( ncid, 'btot', vid )
        iernc= nf90_get_var( ncid, vid, btot)
        iernc=nf90_inq_varid( ncid, 'cormag', vid )
        iernc= nf90_get_var( ncid, vid, cormag)
        iernc=nf90_inq_varid( ncid, 'dipang', vid )
        iernc= nf90_get_var( ncid, vid, dipang)
        iernc=nf90_close(ncid)    

      if(mpi_id.eq.0) then
       write(iulog,*) '  ion_read_wam_init '   
       write(iulog,*) maxval(glon),   minval(glon), ' GLON-ion '
       write(iulog,*) maxval(glat),   minval(glat), ' GLAT-ion ' 
        

       write(iulog,*)  maxval(btot),   minval(btot), ' BTOT-ion '    
       write(iulog,*)  maxval(cormag),minval(cormag),' CORMAG-ion '   
       write(iulog,*)  maxval(dipang),minval(dipang),' DIPANG-ion '  
       write(iulog,*)  ' completed ION_read_wam_init'
      endif

      RETURN    ! Here RETURN is a temporary FIX of mpif.h MPI_REAL8/mpi_integer for THEIA
                  ! ALL PEs read nc-file on THEIA
      end  subroutine ion_read_wam_init
!23456
       END MODULE IDEA_ION_INPUT
!
!Ascii reader of TIROS data
!
      SUBROUTINE tiros_init(emaps,cmaps,djspectra)
!
! Data location: /scratch3/NCEPDEV/swpc/save/Valery.Yudin/BASE_SVN/BASE_WAM_DATA
!   
!      CHARACTER (LEN=*), parameter :: filename3=
!     & '/scratch3/NCEPDEV/swpc/save/Valery.Yudin/BASE_SVN/BASE_WAM_DATA/Tiros/ionprof'
!       CHARACTER (LEN=*), parameter :: filename7=
!     & '/scratch3/NCEPDEV/swpc/save/Valery.Yudin/BASE_SVN/BASE_WAM_DATA/Tiros/tiros_spectra'

      CHARACTER (LEN=*), parameter :: filename3='ionprof'
      CHARACTER (LEN=*), parameter :: filename7='tiros_spectra'
!
      CHARACTER (LEN=100) :: string_dum
      INTEGER :: istat, i
      INTEGER, parameter :: UNIT3=103
      INTEGER, parameter :: UNIT7=107
      real  :: emaps(21,20,7),cmaps(21,20,7),djspectra(15,21)

      open(UNIT=UNIT3,file=trim(filename3),status='old',
     &form='formatted', iostat=istat)
         READ (UNIT3,99001) emaps
         READ (UNIT3,99001) cmaps
       CLOSE(UNIT=UNIT3)
99001 FORMAT (1x,6E13.6)
      
      open(UNIT=UNIT7,file=trim(filename7),status='old',
     &form='formatted', iostat=istat)

         read(UNIT=UNIT7,fmt=*) string_dum
         read(UNIT=UNIT7,fmt=*) string_dum
         read(UNIT7,fmt=*)
         do iband=1,21
         read(UNIT7,fmt=*) string_dum
         read(UNIT7,fmt=*)
         read(UNIT=UNIT7,fmt='(1X,5e10.4)')(djspectra(iflux,iband),
     &   iflux=1,15)
         read(UNIT7,fmt=*)
         enddo
       CLOSE(UNIT=UNIT7)
         return
         end subroutine tiros_init
