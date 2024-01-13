! This program provides a unit test to verify that the variables in
! module_fv3_io_def.F90 are accessible.

program test_fv3_io_def

  use module_fv3_io_def

  implicit none

  num_pes_fcst = 2
  if (num_pes_fcst .ne. 2) stop 1
  wrttasks_per_group = 4
  if (wrttasks_per_group .ne. 4) stop 2
  write_groups = 8
  if (write_groups .ne. 8) stop 3
  n_group = 16
  if (n_group .ne. 16) stop 4
  num_files = 32
  if (num_files .ne. 32) stop 5
  nbdlphys = 64
  if (nbdlphys .ne. 64) stop 6
  iau_offset = 128
  if (iau_offset .ne. 128) stop 7

  lflname_fulltime = .true.
  if (.not. lflname_fulltime) stop 8
  time_unlimited = .true.
  if (.not. time_unlimited) stop 9

  allocate(filename_base(1))
  filename_base(1)="abc"
  if (trim(filename_base(1)).ne."abc") stop 10
  allocate(output_file(1))
  output_file(1)="def"
  if (trim(output_file(1)).ne."def") stop 11

  allocate(lead_wrttask(10))
  lead_wrttask = 42
  if (.not.all(lead_wrttask(1:10).eq.42)) stop 12
  allocate(last_wrttask(12))
  last_wrttask = 43
  if (.not.all(last_wrttask(1:12).eq.43)) stop 13

  allocate(output_grid(2))
  output_grid(1) = "ABC"
  output_grid(2) = "DEF"
  if (trim(output_grid(1)).ne."ABC") stop 14
  if (trim(output_grid(2)).ne."DEF") stop 15

  allocate(imo(6))
  imo(1:6) = 17
  if (.not.all(imo.eq.17)) stop 16
  allocate(jmo(6))
  jmo(1:6) = 18
  if (.not.all(jmo.eq.18)) stop 17

  allocate(cen_lon(8))
  cen_lon(1:8) = 2.1
  if(.not.all(cen_lon.eq.2.1)) stop 18
  allocate(cen_lat(8))
  cen_lat(1:8) = 3.2
  if(.not.all(cen_lat.eq.3.2)) stop 19

  allocate(lon1(10))
  lon1(1:10) = 4.3
  if(.not.all(lon1.eq.4.3)) stop 20
  allocate(lat1(10))
  lat1(1:10) = 5.4
  if(.not.all(lat1.eq.5.4)) stop 21

  allocate(lon2(12))
  lon2(1:12) = 6.5
  if(.not.all(lon2.eq.6.5)) stop 22
  allocate(lat2(12))
  lat2(1:12) = 7.6
  if(.not.all(lat2.eq.7.6)) stop 23

  allocate(dlon(15))
  dlon(1:15) = 8.7
  if(.not.all(dlon.eq.8.7)) stop 24
  allocate(dlat(15))
  dlat(1:15) = 9.8
  if(.not.all(dlat.eq.9.8)) stop 25

  allocate(stdlat1(18))
  stdlat1(1:18) = 10.9
  if(.not.all(stdlat1.eq.10.9)) stop 26
  allocate(stdlat2(18))
  stdlat2(1:18) = 11.0
  if(.not.all(stdlat2.eq.11.0)) stop 27

  allocate(dx(20))
  dx(1:20) = 12.1
  if(.not.all(dx.eq.12.1)) stop 28
  allocate(dy(20))
  dy(1:20) = 13.2
  if(.not.all(dy.eq.13.2)) stop 29

  allocate(ideflate(25))
  ideflate(1:25) = 123
  if(.not.all(ideflate.eq.123)) stop 30
  allocate(quantize_nsd(25))
  quantize_nsd(1:25) = 123
  if(.not.all(quantize_nsd.eq.123)) stop 31
  allocate(zstandard_level(25))
  zstandard_level(1:25) = 123
  if(.not.all(zstandard_level.eq.123)) stop 32

  allocate(quantize_mode(5))
  quantize_mode(1:5) = "xyz"
  if(.not.trim(quantize_mode(5)).eq."xyz") stop 33

  allocate(ichunk2d(30))
  ichunk2d(1:30) = 50
  if(.not.all(ichunk2d.eq.50)) stop 34
  allocate(jchunk2d(30))
  jchunk2d(1:30) = 51
  if(.not.all(jchunk2d.eq.51)) stop 35
  allocate(ichunk3d(30))
  ichunk3d(1:30) = 52
  if(.not.all(ichunk3d.eq.52)) stop 36
  allocate(jchunk3d(30))
  jchunk3d(1:30) = 53
  if(.not.all(jchunk3d.eq.53)) stop 37
  allocate(kchunk3d(30))
  kchunk3d(1:30) = 54
  if(.not.all(kchunk3d.eq.54)) stop 38

end program test_fv3_io_def
