!vay-2015
!               iounitdef.f
!
      MODULE IDEA_IO_UNITS
!
! list of XXX-units in open(unit=XXX, file=IDEA_XXXX, ......
!  
! names of namelist files . in WAM-physics
!                              NmL_solar='solar_in'
!                              NmL_ion='ion_in'
!                              NmL_gwp='gwp_in'
!
! suggestion IDEA_UNITS from [111 to 141]
! don't use with two-digital UNITS of GFS/NEMS 
!
      integer, parameter :: ch100 = 256
      
      character(100), parameter :: NmL_solar='solar_in'
      character(100), parameter :: NmL_ion  ='ion_in'
      character(100), parameter :: NmL_gwp  ='gwp_in'
!
!    efield.f:      open(unit=unit,file=trim(locfn)   3-times
! read-units for WAM-physics  nml + files
!
      integer, parameter :: iulog = 6
      integer, parameter :: nlun_solar = 111
      integer, parameter :: nlun_ion = 112   
      integer, parameter :: nlun_gwp = 114
! efield.f ....... open(unit=unit, where unit =10, 11, 600)  
      integer, parameter :: lun_ef       = 10
      integer, parameter :: lun_efld11   = 11
      integer, parameter :: lun_efld600  =600

      integer, parameter :: lun1_co2 = 30
      integer, parameter :: lun2_co2 = 31
!
! efield.f ....... open(unit=unit, where unit =11)
!
!     efield_lflux_file='global_idea_coeff_lflux.dat',
!     efield_hflux_file='global_idea_coeff_hflux.dat',
!     efield_wei96_file='global_idea_wei96.cofcnts'

!co2hc.f:        open(30,file='global_idea_coeff_lte.150',status = 'OLD')
!co2hc.f:        open(31,file='global_idea_coeff_lte.360',status = 'OLD')
!co2hc.f:        open(32,file='global_idea_coeff_lte.540',status = 'OLD')
!co2hc.f:        open(33,file='global_idea_coeff_lte.720',status = 'OLD')
! h2oc.f         OPEN(11,FILE='global_idea_ggww_in4.par',STATUS='OLD')
!                OPEN(71,FILE='global_idea_h2ort_kg7t.par',STATUS='OLD')
!                OPEN(11,FILE='global_idea_ggww_in1.par',STATUS='OLD')
!                OPEN(71,FILE='global_idea_h2ovb_kg7t.par',STATUS='OLD')
      END MODULE IDEA_IO_UNITS

!  --- ...  input units
!      integer, parameter :: NIMICPH = 1
!      integer, parameter :: NISIGI  = 11
!      integer, parameter :: NISIGI2 = 12
!      integer, parameter :: NISFCI  = 14
!      integer, parameter :: NICO2TR = 15
!      integer, parameter :: NIOFRAD = 16
!      integer, parameter :: NIMTNVR = 24
!      integer, parameter :: NIDTBTH = 27
!      integer, parameter :: NIO3PRD = 28
!      integer, parameter :: NIO3LOS = 29
!      integer, parameter :: NINAMSF = 35
!      integer, parameter :: NICLTUN = 43
!      integer, parameter :: NIO3CLM = 48
!      integer, parameter :: NOSIGR1 = 51
!      integer, parameter :: NOSIGR2 = 52
!      integer, parameter :: NOSFCR  = 53
!      integer, parameter :: NOSIGF  = 61
!      integer, parameter :: NOSFCF  = 62
!      integer, parameter :: NOFLXF  = 63
!      integer, parameter :: NOD3DF  = 64
!      integer, parameter :: NOAERF  = 65    ! for g2d_fld
!      integer, parameter :: NOG3DF  = 69
!      integer, parameter :: NOLOGF  = 99
!      integer, parameter :: NISFCYC = 101
!      integer, parameter :: NIAERCM = 102
!      integer, parameter :: NIRADSF = 102!
!      integer, parameter :: NICO2CN = 102
!  --- ... output units
