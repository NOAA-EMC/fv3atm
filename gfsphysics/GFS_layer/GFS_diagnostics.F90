module GFS_diagnostics

!-----------------------------------------------------------------------
!    GFS_diagnostics_mod defines a data type and contains the routine 
!    to populate said type with diagnostics from the GFS physics for
!    use by the modeling system for output
!-----------------------------------------------------------------------
 
  use machine,            only: kind_phys

  !--- GFS_typedefs ---
  use GFS_typedefs,       only: GFS_control_type,  GFS_statein_type,  &
                                GFS_stateout_type, GFS_sfcprop_type,  &
                                GFS_coupling_type, GFS_grid_type,     &
                                GFS_tbd_type,      GFS_cldprop_type,  &
                                GFS_radtend_type,  GFS_diag_type,  &
                                GFS_init_type
  implicit none
  private

  !--- private data type definition ---
  type data_subtype
    real(kind=kind_phys), dimension(:),   pointer :: var2  => NULL()
    real(kind=kind_phys), dimension(:),   pointer :: var21 => NULL()
    real(kind=kind_phys), dimension(:,:), pointer :: var3  => NULL()
  end type data_subtype

  !--- data type definition for use with GFDL FMS diagnostic manager until write component is working
  type GFS_externaldiag_type
    integer :: id
    integer :: axes
    logical :: time_avg
    character(len=64)    :: time_avg_kind
    character(len=64)    :: mod_name
    character(len=64)    :: name
    character(len=128)   :: desc
    character(len=64)    :: unit
    character(len=64)    :: mask
    character(len=64)    :: intpl_method
    real(kind=kind_phys) :: cnvfac
    type(data_subtype), dimension(:), allocatable :: data
   end type GFS_externaldiag_type

  !--- public data type ---
  public  GFS_externaldiag_type

  !--- public interfaces ---
  public  GFS_externaldiag_populate
 
  CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-------------------------------------------------------------------------      
!--- GFS_externaldiag_populate ---
!-------------------------------------------------------------------------      
!    creates and populates a data type with GFS physics diagnostic 
!    variables which is then handed off to the IPD for use by the model
!    infrastructure layer to output as needed.  The data type includes 
!    names, units, conversion factors, etc.  There is no copying of data, 
!    but instead pointers are associated to the internal representation 
!    of each individual physics diagnostic.
!-------------------------------------------------------------------------      
  subroutine GFS_externaldiag_populate (ExtDiag, Model, Statein, Stateout, Sfcprop, Coupling,  &
                                        Grid, Tbd, Cldprop, Radtend, IntDiag, Init_parm)
!---------------------------------------------------------------------------------------------!
!   DIAGNOSTIC_METADATA                                                                       !
!     ExtDiag%id                   [integer ]   switch to turn on/off variable output         !
!     ExtDiag%axes                 [integer ]   dimensionality of variable (2 or 3)           !
!     ExtDiag%time_avg             [logical ]   bucketed accumulation time average            !
!     ExtDiag%time_avg_kind        [char*64 ]   time average period                           !
!     ExtDiag%mod_name             [char*64 ]   classification of the variable                !
!     ExtDiag%name                 [char*64 ]   output name for variable                      !
!     ExtDiag%desc                 [char*128]   long description of field                     !
!     ExtDiag%unit                 [char*64 ]   units associated with field                   !
!     ExtDiag%mask                 [char*64 ]   description of mask-type                      !
!     ExtDiag%intpl_method         [char*64 ]   method to use for interpolation               !
!     ExtDiag%cnvfac               [real*8  ]   conversion factor to output specified units   !
!     ExtDiag%data(nb)%var2(:)     [real*8  ]   pointer to 2D data [=> null() for a 3D field] !
!     ExtDiag%data(nb)%var21(:)    [real*8  ]   pointer to 2D data for ratios                 !
!     ExtDiag%data(nb)%var3(:,:)   [real*8  ]   pointer to 3D data [=> null() for a 2D field] !
!---------------------------------------------------------------------------------------------!

      implicit none
!
!  ---  interface variables
    type(GFS_externaldiag_type),  intent(inout) :: ExtDiag(:)
    type(GFS_control_type),       intent(in)    :: Model
    type(GFS_statein_type),       intent(in)    :: Statein(:)
    type(GFS_stateout_type),      intent(in)    :: Stateout(:)
    type(GFS_sfcprop_type),       intent(in)    :: Sfcprop(:)
    type(GFS_coupling_type),      intent(in)    :: Coupling(:)
    type(GFS_grid_type),          intent(in)    :: Grid(:)
    type(GFS_tbd_type),           intent(in)    :: Tbd(:)
    type(GFS_cldprop_type),       intent(in)    :: Cldprop(:)
    type(GFS_radtend_type),       intent(in)    :: Radtend(:)
    type(GFS_diag_type),          intent(in)    :: IntDiag(:)
    type(GFS_init_type),          intent(in)    :: Init_parm

!--- local variables
    integer :: idt, idx, num, nb, nblks, NFXR
    character(len=2) :: xtra
    real(kind=kind_phys), parameter :: cn_one = 1._kind_phys
    real(kind=kind_phys), parameter :: cn_100 = 100._kind_phys
    real(kind=kind_phys), parameter :: cn_th  = 1000._kind_phys
    real(kind=kind_phys), parameter :: cn_hr  = 3600._kind_phys

    NFXR = Model%NFXR
    nblks = size(Statein)

    ExtDiag(:)%id = -99
    ExtDiag(:)%axes = -99
    ExtDiag(:)%cnvfac = cn_one
    ExtDiag(:)%time_avg = .FALSE.
    ExtDiag(:)%time_avg_kind = ''
    ExtDiag(:)%mask = ''
    ExtDiag(:)%name = ''
    ExtDiag(:)%intpl_method = 'nearest_stod'

    idx = 0 

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ALBDO_ave'
    ExtDiag(idx)%desc = 'surface albedo'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_100
    ExtDiag(idx)%mask = 'positive_flux'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2  => IntDiag(nb)%fluxr(:,3)
      ExtDiag(idx)%data(nb)%var21 => IntDiag(nb)%fluxr(:,4)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'DLWRF'
    ExtDiag(idx)%desc = 'surface downward longwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%dlwsfc(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'DLWRFI'
    ExtDiag(idx)%desc = 'instantaneous surface downward longwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%dlwsfci(:)
    enddo


    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ULWRF'
    ExtDiag(idx)%desc = 'surface upward longwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%ulwsfc(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ULWRFI'
    ExtDiag(idx)%desc = 'instantaneous surface upward longwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%ulwsfci(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'DSWRF'
    ExtDiag(idx)%desc = 'averaged surface downward shortwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,4)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'DSWRFI'
    ExtDiag(idx)%desc = 'instantaneous surface downward shortwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%dswsfci(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'USWRF'
    ExtDiag(idx)%desc = 'averaged surface upward shortwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,3)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'USWRFI'
    ExtDiag(idx)%desc = 'instantaneous surface upward shortwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%uswsfci(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'duvb_ave'
    ExtDiag(idx)%desc = 'UV-B Downward Solar Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,21)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'cduvb_ave'
    ExtDiag(idx)%desc = 'Clear sky UV-B Downward Solar Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,22)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'vbdsf_ave'
    ExtDiag(idx)%desc = 'Visible Beam Downward Solar Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,24)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'vddsf_ave'
    ExtDiag(idx)%desc = 'Visible Diffuse Downward Solar Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,25)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'nbdsf_ave'
    ExtDiag(idx)%desc = 'Near IR Beam Downward Solar Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,26)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'nddsf_ave'
    ExtDiag(idx)%desc = 'Near IR Diffuse Downward Solar Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,27)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'csulf_avetoa'
    ExtDiag(idx)%desc = 'Clear Sky Upward Long Wave Flux at toa'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_lw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,28)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'csusf_avetoa'
    ExtDiag(idx)%desc = 'Clear Sky Upward Short Wave Flux at toa'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,29)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'csdlf_ave'
    ExtDiag(idx)%desc = 'Clear Sky Downward Long Wave Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_lw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,30)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'csusf_ave'
    ExtDiag(idx)%desc = 'Clear Sky Upward Short Wave Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,31)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'csdsf_ave'
    ExtDiag(idx)%desc = 'Clear Sky Downward Short Wave Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,32)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'csulf_ave'
    ExtDiag(idx)%desc = 'Clear Sky Upward Long Wave Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_lw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,33)
    enddo


    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'DSWRFtoa'
    ExtDiag(idx)%desc = 'top of atmos downward shortwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,23)
    enddo


    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'USWRFtoa'
    ExtDiag(idx)%desc = 'top of atmos upward shortwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,2)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ULWRFtoa'
    ExtDiag(idx)%desc = 'top of atmos upward longwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_lw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,1)
    enddo
!    if(mpp_pe()==mpp_root_pe())print *,'in gfdl_diag_register,bf ULWRFtoa,idx=',idx

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'TCDC_aveclm'
    ExtDiag(idx)%desc = 'atmos column total cloud cover'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_100
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,17)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'TCDC_avebndcl'
    ExtDiag(idx)%desc = 'boundary layer cloud layer total cloud cover'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_100
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,18)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'TCDCcnvcl'
    ExtDiag(idx)%desc = 'convective cloud layer total cloud cover'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_100
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Cldprop(nb)%cv(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'PREScnvclt'
    ExtDiag(idx)%desc = 'pressure at convective cloud top level'
    ExtDiag(idx)%unit = 'pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%mask = 'cldmask'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Cldprop(nb)%cvt(:)
      ExtDiag(idx)%data(nb)%var21 => Cldprop(nb)%cv(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'PREScnvclb'
    ExtDiag(idx)%desc = 'pressure at convective cloud bottom level'
    ExtDiag(idx)%unit = 'pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%mask = 'cldmask'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Cldprop(nb)%cvb(:)
      ExtDiag(idx)%data(nb)%var21 => Cldprop(nb)%cv(:)
    enddo
!    if(mpp_pe()==mpp_root_pe())print *,'in gfdl_diag_register,af PREScnvclb,idx=',idx


    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'TCDC_avehcl'
    ExtDiag(idx)%desc = 'high cloud level total cloud cover'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_100
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,5)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'PRES_avehct'
    ExtDiag(idx)%desc = 'pressure high cloud top level'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    ExtDiag(idx)%mask = "cldmask_ratio"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,8)
      ExtDiag(idx)%data(nb)%var21 => IntDiag(nb)%fluxr(:,5)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'PRES_avehcb'
    ExtDiag(idx)%desc = 'pressure high cloud bottom level'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    ExtDiag(idx)%mask = "cldmask_ratio"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,11)
      ExtDiag(idx)%data(nb)%var21 => IntDiag(nb)%fluxr(:,5)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'TEMP_avehct'
    ExtDiag(idx)%desc = 'temperature high cloud top level'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    ExtDiag(idx)%mask = "cldmask_ratio"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,14)
      ExtDiag(idx)%data(nb)%var21 => IntDiag(nb)%fluxr(:,5)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'TCDC_avemcl'
    ExtDiag(idx)%desc = 'mid cloud level total cloud cover'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_100
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,6)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'PRES_avemct'
    ExtDiag(idx)%desc = 'pressure middle cloud top level'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    ExtDiag(idx)%mask = "cldmask_ratio"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,9)
      ExtDiag(idx)%data(nb)%var21 => IntDiag(nb)%fluxr(:,6)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'PRES_avemcb'
    ExtDiag(idx)%desc = 'pressure middle cloud bottom level'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    ExtDiag(idx)%mask = "cldmask_ratio"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,12)
      ExtDiag(idx)%data(nb)%var21 => IntDiag(nb)%fluxr(:,6)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'TEMP_avemct'
    ExtDiag(idx)%desc = 'temperature middle cloud top level'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    ExtDiag(idx)%mask = "cldmask_ratio"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,15)
      ExtDiag(idx)%data(nb)%var21 => IntDiag(nb)%fluxr(:,6)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'TCDC_avelcl'
    ExtDiag(idx)%desc = 'low cloud level total cloud cover'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_100
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,7)
    enddo

   idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'PRES_avelct'
    ExtDiag(idx)%desc = 'pressure low cloud top level'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    ExtDiag(idx)%mask = "cldmask_ratio"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,10)
      ExtDiag(idx)%data(nb)%var21 => IntDiag(nb)%fluxr(:,7)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'PRES_avelcb'
    ExtDiag(idx)%desc = 'pressure low cloud bottom level'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    ExtDiag(idx)%mask = "cldmask_ratio"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,13)
      ExtDiag(idx)%data(nb)%var21 => IntDiag(nb)%fluxr(:,7)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'TEMP_avelct'
    ExtDiag(idx)%desc = 'temperature low cloud top level'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    ExtDiag(idx)%mask = "cldmask_ratio"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,16)
      ExtDiag(idx)%data(nb)%var21 => IntDiag(nb)%fluxr(:,7)
    enddo
!    if(mpp_pe()==mpp_root_pe())print *,'in gfdl_diag_register,af TEMP_avelct,idx=',idx

!
!--- accumulated diagnostics ---
    do num = 1,NFXR
      write (xtra,'(I2.2)') num 
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'fluxr_'//trim(xtra)
      ExtDiag(idx)%desc = 'fluxr diagnostic '//trim(xtra)//' - GFS radiation'
      ExtDiag(idx)%unit = 'XXX'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%fluxr(:,num)
      enddo
    enddo

!--- the next two appear to be appear to be coupling fields in gloopr
!--- each has four elements
!rab    do num = 1,4
!rab      write (xtra,'(I1)') num 
!rab      idx = idx + 1
!rab      ExtDiag(idx)%axes = 2
!rab      ExtDiag(idx)%name = 'dswcmp_'//trim(xtra)
!rab      ExtDiag(idx)%desc = 'dswcmp dagnostic '//trim(xtra)//' - GFS radiation'
!rab      ExtDiag(idx)%unit = 'XXX'
!rab      ExtDiag(idx)%mod_name = 'gfs_phys'
!rab      do nb = 1,nblks
!rab        ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%dswcmp(:,num)
!rab      enddo
!rab    enddo
!rab
!rab    do num = 1,4
!rab      write (xtra,'(I1)') num 
!rab      idx = idx + 1
!rab      ExtDiag(idx)%axes = 2
!rab      ExtDiag(idx)%name = 'uswcmp_'//trim(xtra)
!rab      ExtDiag(idx)%desc = 'uswcmp dagnostic '//trim(xtra)//' - GFS radiation'
!rab      ExtDiag(idx)%unit = 'XXX'
!rab      ExtDiag(idx)%mod_name = 'gfs_phys'
!rab      do nb = 1,nblks
!rab        ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%uswcmp(:,num)
!rab      enddo
!rab    enddo

! DH gfortran cannot point to members of arrays of derived types such
! as IntDiag(nb)%topfsw(:)%upfxc (the compilation succeeds, but the
! pointers do not reference the correct data and the output either
! contains garbage (Inf, NaN), or the netCDF I/O layer crashes.
#ifndef __GFORTRAN__
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'sw_upfxc'
    ExtDiag(idx)%desc = 'total sky upward sw flux at toa - GFS radiation'
    ExtDiag(idx)%unit = 'w/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%topfsw(:)%upfxc
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'sw_dnfxc'
    ExtDiag(idx)%desc = 'total sky downward sw flux at toa - GFS radiation'
    ExtDiag(idx)%unit = 'w/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%topfsw(:)%dnfxc
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'sw_upfx0'
    ExtDiag(idx)%desc = 'clear sky upward sw flux at toa - GFS radiation'
    ExtDiag(idx)%unit = 'w/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%topfsw(:)%upfx0
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'lw_upfxc'
    ExtDiag(idx)%desc = 'total sky upward lw flux at toa - GFS radiation'
    ExtDiag(idx)%unit = 'w/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%topflw(:)%upfxc
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'lw_upfx0'
    ExtDiag(idx)%desc = 'clear sky upward lw flux at toa - GFS radiation'
    ExtDiag(idx)%unit = 'w/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%topflw(:)%upfx0
    enddo
#endif

!--- physics accumulated diagnostics ---
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ssrun_acc'
    ExtDiag(idx)%desc = 'surface storm water runoff - GFS lsm'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%srunoff(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'evbs_ave'
    ExtDiag(idx)%desc = 'Direct Evaporation from Bare Soil - GFS lsm'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%evbsa(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'evcw_ave'
    ExtDiag(idx)%desc = 'Canopy water evaporation - GFS lsm'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%evcwa(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'snohf'
    ExtDiag(idx)%desc = 'Snow Phase Change Heat Flux - GFS lsm'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%snohfa(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'trans_ave'
    ExtDiag(idx)%desc = 'transpiration - GFS lsm'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%transa(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'sbsno_ave'
    ExtDiag(idx)%desc = 'Sublimation (evaporation from snow) - GFS lsm'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%sbsnoa(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'snowc_ave'
    ExtDiag(idx)%desc = 'snow cover - GFS lsm'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%cnvfac = cn_100
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%snowca(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'soilm'
    ExtDiag(idx)%desc = 'total column soil moisture content'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%mask = "land_only"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2  => IntDiag(nb)%soilm(:)
      ExtDiag(idx)%data(nb)%var21 => Sfcprop(nb)%slmsk(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tmpmin2m'
    ExtDiag(idx)%desc = 'min temperature at 2m height'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%tmpmin(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tmpmax2m'
    ExtDiag(idx)%desc = 'max temperature at 2m height'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%tmpmax(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dusfc'
    ExtDiag(idx)%desc = 'surface zonal momentum flux'
    ExtDiag(idx)%unit = 'N/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%dusfc(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dvsfc'
    ExtDiag(idx)%desc = 'surface meridional momentum flux'
    ExtDiag(idx)%unit = 'N/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%dvsfc(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'shtfl_ave'
    ExtDiag(idx)%desc = 'surface sensible heat flux'
    ExtDiag(idx)%unit = 'w/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%dtsfc(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'lhtfl_ave'
    ExtDiag(idx)%desc = 'surface latent heat flux'
    ExtDiag(idx)%unit = 'w/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%dqsfc(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'totprcp_ave'
    ExtDiag(idx)%desc = 'surface precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'full'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%totprcp(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'totprcpb_ave'
    ExtDiag(idx)%desc = 'bucket surface precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%totprcpb(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'gflux_ave'
    ExtDiag(idx)%desc = 'surface ground heat flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
!    ExtDiag(idx)%mask = "land_ice_only"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2  => IntDiag(nb)%gflux(:)
!      ExtDiag(idx)%data(nb)%var21 => Sfcprop(nb)%slmsk(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dlwsfc'
    ExtDiag(idx)%desc = 'time accumulated downward lw flux at surface- GFS physics'
    ExtDiag(idx)%unit = 'w/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%dlwsfc(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ulwsfc'
    ExtDiag(idx)%desc = 'time accumulated upward lw flux at surface- GFS physics'
    ExtDiag(idx)%unit = 'w/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%ulwsfc(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'sunsd_acc'
    ExtDiag(idx)%desc = 'Sunshine Duration'
    ExtDiag(idx)%unit = 's'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%suntim(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'watr_acc'
    ExtDiag(idx)%desc = 'total water runoff'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%runoff(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'pevpr_ave'
    ExtDiag(idx)%desc = 'averaged potential evaporation rate'
    ExtDiag(idx)%unit = 'W/M**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%ep(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'cwork_ave'
    ExtDiag(idx)%desc = 'cloud work function (valid only with sas)'
    ExtDiag(idx)%unit = 'J/kg'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%cldwrk(:)
    enddo


    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'u-gwd_ave'
    ExtDiag(idx)%desc = 'surface zonal gravity wave stress'
    ExtDiag(idx)%unit = 'N/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%dugwd(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'v-gwd_ave'
    ExtDiag(idx)%desc = 'surface meridional gravity wave stress'
    ExtDiag(idx)%unit = 'N/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%dvgwd(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'psmean'
    ExtDiag(idx)%desc = 'surface pressure'
    ExtDiag(idx)%unit = 'kPa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%psmean(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'cnvprcp_ave'
    ExtDiag(idx)%desc = 'averaged surface convective precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'full'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%cnvprcp(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'cnvprcpb_ave'
    ExtDiag(idx)%desc = 'averaged bucket surface convective precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%cnvprcpb(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'cnvprcp'
    ExtDiag(idx)%desc = 'surface convective precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%cnvprcp(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'spfhmin2m'
    ExtDiag(idx)%desc = 'minimum specific humidity'
    ExtDiag(idx)%unit = 'kg/kg'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%spfhmin(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'spfhmax2m'
    ExtDiag(idx)%desc = 'maximum specific humidity'
    ExtDiag(idx)%unit = 'kg/kg'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%spfhmax(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'u10mmax'
    ExtDiag(idx)%desc = 'maximum (magnitude) u-wind'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'vector_bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%u10mmax(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'v10mmax'
    ExtDiag(idx)%desc = 'maximum (magnitude) v-wind'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'vector_bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%v10mmax(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'wind10mmax'
    ExtDiag(idx)%desc = 'maximum wind speed'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%wind10mmax(:)
    enddo
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'u10max'
    ExtDiag(idx)%desc = 'hourly maximum (magnitude) u-wind'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'vector_bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%u10max(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'v10max'
    ExtDiag(idx)%desc = 'hourly maximum (magnitude) v-wind'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'vector_bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%v10max(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'spd10max'
    ExtDiag(idx)%desc = 'hourly maximum wind speed'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%spd10max(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 't02max'
    ExtDiag(idx)%desc = 'max hourly 2m Temperature'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%t02max(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 't02min'
    ExtDiag(idx)%desc = 'min hourly 2m Temperature'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%t02min(:)
    enddo
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'rh02max'
    ExtDiag(idx)%desc = 'max hourly 2m RH'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%rh02max(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'rh02min'
    ExtDiag(idx)%desc = 'min hourly 2m RH'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%rh02min(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'rain'
    ExtDiag(idx)%desc = 'total rain at this time step'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%rain(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'rainc'
    ExtDiag(idx)%desc = 'convective rain at this time step'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%rainc(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ice'
    ExtDiag(idx)%desc = 'ice fall at this time step'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%ice(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'snow'
    ExtDiag(idx)%desc = 'snow fall at this time step'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%snow(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'graupel'
    ExtDiag(idx)%desc = 'graupel fall at this time step'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%graupel(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'totice_ave'
    ExtDiag(idx)%desc = 'surface ice precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'full'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%totice(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'toticeb_ave'
    ExtDiag(idx)%desc = 'bucket surface ice precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%toticeb(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'totsnw_ave'
    ExtDiag(idx)%desc = 'surface snow precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'full'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%totsnw(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'totsnwb_ave'
    ExtDiag(idx)%desc = 'bucket surface snow precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%totsnwb(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'totgrp_ave'
    ExtDiag(idx)%desc = 'surface graupel precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'full'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%totgrp(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'totgrpb_ave'
    ExtDiag(idx)%desc = 'bucket surface graupel precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%totgrpb(:)
    enddo

!    if(mpp_pe()==mpp_root_pe())print *,'in gfdl_diag_register,af totgrp,idx=',idx

!--- physics instantaneous diagnostics ---
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'u10m'
    ExtDiag(idx)%desc = '10 meter u wind'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'vector_bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%u10m(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'v10m'
    ExtDiag(idx)%desc = '10 meter v wind'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'vector_bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%v10m(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dpt2m'
    ExtDiag(idx)%desc = '2 meter dew point temperature'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%dpt2m(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'hgt_hyblev1'
    ExtDiag(idx)%desc = 'layer 1 height'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%zlvl(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'psurf'
    ExtDiag(idx)%desc = 'surface pressure'
    ExtDiag(idx)%unit = 'Pa'
    ExtDiag(idx)%mask = 'pseudo_ps'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%psurf(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'hpbl'
    ExtDiag(idx)%desc = 'surface planetary boundary layer height'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Tbd(nb)%hpbl(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'pwat'
    ExtDiag(idx)%desc = 'atmos column precipitable water'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%pwat(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tmp_hyblev1'
    ExtDiag(idx)%desc = 'layer 1 temperature'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%t1(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'spfh_hyblev1'
    ExtDiag(idx)%desc = 'layer 1 specific humidity'
    ExtDiag(idx)%unit = 'kg/kg'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%q1(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ugrd_hyblev1'
    ExtDiag(idx)%desc = 'layer 1 zonal wind'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'vector_bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%u1(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'vgrd_hyblev1'
    ExtDiag(idx)%desc = 'layer 1 meridional wind'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'vector_bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%v1(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'sfexc'
    ExtDiag(idx)%desc = 'Exchange Coefficient'
    ExtDiag(idx)%unit = 'kg/m2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%chh(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'acond'
    ExtDiag(idx)%desc = 'Aerodynamic conductance'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%cmm(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dlwsfci'
    ExtDiag(idx)%desc = 'instantaneous sfc downward lw flux'
    ExtDiag(idx)%unit = 'w/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%dlwsfci(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ulwsfci'
    ExtDiag(idx)%desc = 'instantaneous sfc upward lw flux'
    ExtDiag(idx)%unit = 'w/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%ulwsfci(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dswsfci'
    ExtDiag(idx)%desc = 'instantaneous sfc downward sw flux'
    ExtDiag(idx)%unit = 'w/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%dswsfci(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'uswsfci'
    ExtDiag(idx)%desc = 'instantaneous sfc upward sw flux'
    ExtDiag(idx)%unit = 'w/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%uswsfci(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dusfci'
    ExtDiag(idx)%desc = 'instantaneous u component of surface stress'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%dusfci(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dvsfci'
    ExtDiag(idx)%desc = 'instantaneous v component of surface stress'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%dvsfci(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'shtfl'
    ExtDiag(idx)%desc = 'instantaneous surface sensible heat net flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%dtsfci(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'lhtfl'
    ExtDiag(idx)%desc = 'instantaneous surface latent heat net flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%dqsfci(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'gfluxi'
    ExtDiag(idx)%desc = 'instantaneous surface ground heat flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%gfluxi(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'pevpr'
    ExtDiag(idx)%desc = 'instantaneous surface potential evaporation'
    ExtDiag(idx)%unit = 'W/M**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%epi(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'wilt'
    ExtDiag(idx)%desc = 'wiltimg point (volumetric)'
    ExtDiag(idx)%unit = 'Proportion'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%smcwlt2(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'fldcp'
    ExtDiag(idx)%desc = 'Field Capacity (volumetric)'
    ExtDiag(idx)%unit = 'fraction'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%smcref2(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'wet1'
    ExtDiag(idx)%desc = 'normalized soil wetness'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
#ifdef CCPP
    if (Model%lsm==Model%lsm_ruc) then
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%wetness(:)
      enddo
    else
#endif
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%wet1(:)
      enddo
#ifdef CCPP
    endif
#endif

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'cpofp'
    ExtDiag(idx)%desc = 'Percent frozen precipitation'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%sr(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'crain_ave'
    ExtDiag(idx)%desc = 'averaged categorical rain'
    ExtDiag(idx)%unit = 'number'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%tdomr(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'csnow_ave'
    ExtDiag(idx)%desc = 'averaged categorical snow'
    ExtDiag(idx)%unit = 'number'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%tdoms(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'cfrzr_ave'
    ExtDiag(idx)%desc = 'averaged categorical freezing rain'
    ExtDiag(idx)%unit = 'number'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%tdomzr(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'cicep_ave'
    ExtDiag(idx)%desc = 'averaged categorical sleet'
    ExtDiag(idx)%unit = 'number'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%tdomip(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'refl_10cm'
    ExtDiag(idx)%desc = 'Radar reflectivity'
    ExtDiag(idx)%unit = 'dBz'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%refl_10cm(:,:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'cnvw'
    ExtDiag(idx)%desc = 'subgrid scale convective cloud water'
    ExtDiag(idx)%unit = 'kg/kg'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    if( Model%ncnvw > 0 ) then
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Tbd(nb)%phy_f3d(:,:,Model%ncnvw)
      enddo
    endif

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'skebu_wts'
    ExtDiag(idx)%desc = 'perturbation velocity'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%skebu_wts(:,:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'skebv_wts'
    ExtDiag(idx)%desc = 'perturbation velocity'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%skebv_wts(:,:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'zmtnblck'
    ExtDiag(idx)%desc = 'level of dividing streamline'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%zmtnblck(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'refdmax'
    ExtDiag(idx)%desc = 'max hourly 1-km agl reflectivity'
    ExtDiag(idx)%unit = 'dBZ'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%refdmax(:)
    enddo
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'refdmax263k'
    ExtDiag(idx)%desc = 'max hourly -10C reflectivity'
    ExtDiag(idx)%unit = 'dBZ'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%refdmax263k(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'sppt_wts'
    ExtDiag(idx)%desc = 'perturbation velocity'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%sppt_wts(:,:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'shum_wts'
    ExtDiag(idx)%desc = 'perturbation velocity'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%shum_wts(:,:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ca_out'
    ExtDiag(idx)%desc = 'Cellular Automata'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%ca_out(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ca_deep'
    ExtDiag(idx)%desc = 'CA deep conv'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%ca_deep(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ca_turb'
    ExtDiag(idx)%desc = 'CA turbulence'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%ca_turb(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ca_shal'
    ExtDiag(idx)%desc = 'CA shallow conv'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%ca_shal(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ca_rad'
    ExtDiag(idx)%desc = 'CA radiation'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%ca_rad(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ca_micro'
    ExtDiag(idx)%desc = 'CA microphys'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%ca_micro(:)
    enddo

  if (Model%ldiag_ugwp) THEN
!
! VAY-2018: Momentum and Temp-re tendencies
! du3dt_pbl dv3dt_pbl  dT3dt_pbl
!
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'zmtb'
    ExtDiag(idx)%desc = 'height of dividing streamline'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%zmtb(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'zogw'
    ExtDiag(idx)%desc = 'height of OGW-launch'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%zogw(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'zlwb'
    ExtDiag(idx)%desc = 'height of LWB-level'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%zlwb(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tau_ogw'
    ExtDiag(idx)%desc = ' OGW vertical MF at launch level'
    ExtDiag(idx)%unit = 'N/m2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%tau_ogw(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tau_mtb'
    ExtDiag(idx)%desc = ' ORO-MTB integrated flux from surface'
    ExtDiag(idx)%unit = 'N/m2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%tau_mtb(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tau_tofd'
    ExtDiag(idx)%desc = ' ORO-TOFD integrated flux from surface'
    ExtDiag(idx)%unit = 'N/m2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%tau_tofd(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tau_ngw'
    ExtDiag(idx)%desc = ' NGW momentum flux at launch level '
    ExtDiag(idx)%unit = 'N/m2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%tau_ngw(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'du3dt_pbl_ugwp'
    ExtDiag(idx)%desc = 'U-tendency due to PBL physics'
    ExtDiag(idx)%unit = 'm/s/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%du3dt_pbl(:,:)
    enddo
!    
! dv3dt_pbl     
    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dv3dt_pbl_ugwp'
    ExtDiag(idx)%desc = 'V-tendency due to PBL physics'
    ExtDiag(idx)%unit = 'm/s/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%dv3dt_pbl(:,:)
    enddo
!    
! dt3dt_pbl     
    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dt3dt_pbl_ugwp'
    ExtDiag(idx)%desc = 'T-tendency due to PBL physics'
    ExtDiag(idx)%unit = 'K/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%dt3dt_pbl(:,:)
    enddo
!
! uav_ugwp 
! 
    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'uav_ugwp'
    ExtDiag(idx)%desc = 'U-daily mean for UGWP'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%uav_ugwp(:,:)
    enddo
!
! tav_ugwp 
! 
    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'tav_ugwp'
    ExtDiag(idx)%desc = 'T-daily mean for UGWP'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%tav_ugwp(:,:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'du3dt_ogw'
    ExtDiag(idx)%desc = 'axz_oro averaged E-W OROGW-tendency'
    ExtDiag(idx)%unit = 'm/s/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%du3dt_ogw(:,:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'du3dt_ngw'
    ExtDiag(idx)%desc = 'axz_oro averaged E-W GWALL-tendency'
    ExtDiag(idx)%unit = 'm/s/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%du3dt_ngw(:,:)
    enddo
!
!    
    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'du3dt_mtb'
    ExtDiag(idx)%desc = 'axz_oro averaged E-W MTB-tendency'
    ExtDiag(idx)%unit = 'm/s/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%du3dt_mtb(:,:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'du3dt_tms'
    ExtDiag(idx)%desc = 'axz_oro averaged E-W TMS-tendency'
    ExtDiag(idx)%unit = 'm/s/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%du3dt_tms(:,:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dudt_tot'
    ExtDiag(idx)%desc = ' dudt_tot averaged E-W dycore-tendency'
    ExtDiag(idx)%unit = 'm/s/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%dudt_tot(:,:)
    enddo

     idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dtdt_tot'
    ExtDiag(idx)%desc = ' dtdt_tot averaged Temp dycore-tendency'
    ExtDiag(idx)%unit = 'Ks'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%dtdt_tot(:,:)
    enddo
   ENDIF



!    if(mpp_pe()==mpp_root_pe())print *,'in gfdl_diag_register,af shum_wts,idx=',idx

!--- three-dimensional variables that need to be handled special when writing 
    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dt3dt_lw'
    ExtDiag(idx)%desc = 'temperature change due to long wave radiation'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%dt3dt(:,:,1)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dt3dt_sw'
    ExtDiag(idx)%desc = 'temperature change due to short wave radiation'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%dt3dt(:,:,2)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dt3dt_pbl'
    ExtDiag(idx)%desc = 'temperature change due to pbl'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%dt3dt(:,:,3)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dt3dt_deepcnv'
    ExtDiag(idx)%desc = 'temperature change due to deep convection'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%dt3dt(:,:,4)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dt3dt_shlwcnv'
    ExtDiag(idx)%desc = 'temperature change due to shallow convection'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%dt3dt(:,:,5)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dt3dt_microphy'
    ExtDiag(idx)%desc = 'temperature change due to micro-physics'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%dt3dt(:,:,6)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dt3dt_orogwd'
    ExtDiag(idx)%desc = 'temperature change due to orographic gravity wave drag'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%dt3dt(:,:,7)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'du3dt_pbl'
    ExtDiag(idx)%desc = 'u momentum change due to PBL'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%du3dt(:,:,1)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dv3dt_pbl'
    ExtDiag(idx)%desc = 'v momentum change due to PBL'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%dv3dt(:,:,1)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'du3dt_orogwd'
    ExtDiag(idx)%desc = 'u momentum change due to orographic gravity wave drag'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%du3dt(:,:,2)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dv3dt_orogwd'
    ExtDiag(idx)%desc = 'v momentum change due to orographic gravity wave drag'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%dv3dt(:,:,2)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'du3dt_deepcnv'
    ExtDiag(idx)%desc = 'u momentum change due to deep convection'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%du3dt(:,:,3)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dv3dt_deepcnv'
    ExtDiag(idx)%desc = 'v momentum change due to deep convection'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%dv3dt(:,:,3)
    enddo


    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'du3dt_cnvgwd'
    ExtDiag(idx)%desc = 'u momentum change due to convective gravity wave drag'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%du3dt(:,:,4)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dv3dt_cnvgwd'
    ExtDiag(idx)%desc = 'v momentum change due to convective gravity wave drag'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag(nb)%dv3dt(:,:,4)
    enddo


!rab
!rab    do num = 1,5+Mdl_parms%pl_coeff
!rab      write (xtra,'(I1)') num 
!rab      idx = idx + 1
!rab      ExtDiag(idx)%axes = 3
!rab      ExtDiag(idx)%name = 'dq3dt_'//trim(xtra)
!rab      ExtDiag(idx)%desc = 'moisture change due to physics '//trim(xtra)//''
!rab      ExtDiag(idx)%unit = 'XXX'
!rab      ExtDiag(idx)%mod_name = 'gfs_phys'
!rab    enddo
!rab
!rab    idx = idx + 1
!rab    ExtDiag(idx)%axes = 3
!rab    !Requires lgocart = .T.
!rab    ExtDiag(idx)%name = 'dqdt_v'
!rab    ExtDiag(idx)%desc = 'instantaneous total moisture tendency'
!rab    ExtDiag(idx)%unit = 'XXX'
!rab    ExtDiag(idx)%mod_name = 'gfs_phys'

!--- Surface diagnostics in gfs_sfc
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'alnsf'
    ExtDiag(idx)%desc = 'mean nir albedo with strong cosz dependency'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%alnsf(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'alnwf'
    ExtDiag(idx)%desc = 'mean nir albedo with weak cosz dependency'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%alnwf(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'alvsf'
    ExtDiag(idx)%desc = 'mean vis albedo with strong cosz dependency'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%alvsf(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'alvwf'
    ExtDiag(idx)%desc = 'mean vis albedo with weak cosz dependency'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%alvwf(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'canopy'
    ExtDiag(idx)%desc = 'canopy water (cnwat in gfs data)'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%canopy(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'f10m'
    ExtDiag(idx)%desc = '10-meter wind speed divided by lowest model wind speed'
    ExtDiag(idx)%unit = 'N/A'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%f10m(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'facsf'
    ExtDiag(idx)%desc = 'fractional coverage with strong cosz dependency'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%facsf(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'facwf'
    ExtDiag(idx)%desc = 'fractional coverage with weak cosz dependency'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%facwf(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ffhh'
    ExtDiag(idx)%desc = 'fh parameter from PBL scheme'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%ffhh(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ffmm'
    ExtDiag(idx)%desc = 'fm parameter from PBL scheme'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%ffmm(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'uustar'
    ExtDiag(idx)%desc = 'uustar surface frictional wind'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%uustar(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'slope'
    ExtDiag(idx)%desc = 'surface slope type'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%slope(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'fice'
    ExtDiag(idx)%desc = 'surface ice concentration (ice=1; no ice=0)'
    ExtDiag(idx)%unit = 'fraction'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%fice(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'hice'
    ExtDiag(idx)%desc = 'sea ice thickness (icetk in gfs_data)'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%hice(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'snoalb'
    ExtDiag(idx)%desc = 'maximum snow albedo in fraction'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%snoalb(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'shdmax'
    ExtDiag(idx)%desc = 'maximum fractional coverage of green vegetation'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%shdmax(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'shdmin'
    ExtDiag(idx)%desc = 'minimum fractional coverage of green vegetation'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%shdmin(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'snowd'
    ExtDiag(idx)%desc = 'surface snow depth'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    ExtDiag(idx)%cnvfac = cn_one/cn_th
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%snowd(:)
    enddo

#ifdef CCPP
    if (Model%lsm == Model%lsm_ruc) then
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'snowfall_acc'
      ExtDiag(idx)%desc = 'total accumulated frozen precipitation'
      ExtDiag(idx)%unit = 'kg m-2'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%snowfallac(:)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'swe_snowfall_acc'
      ExtDiag(idx)%desc = 'accumulated water equivalent of frozen precipitation'
      ExtDiag(idx)%unit = 'kg m-2'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%acsnow(:)
      enddo
    endif
#endif

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'crain'
    ExtDiag(idx)%desc = 'instantaneous categorical rain'
    ExtDiag(idx)%unit = 'number'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    ExtDiag(idx)%cnvfac = cn_one
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%srflag(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'stype'
    ExtDiag(idx)%desc = 'soil type in integer 1-9'
    ExtDiag(idx)%unit = 'number'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%stype(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'q2m'
    ExtDiag(idx)%desc = '2m specific humidity'
    ExtDiag(idx)%unit = 'kg/kg'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%q2m(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 't2m'
    ExtDiag(idx)%desc = '2m temperature'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%t2m(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tsfc'
    ExtDiag(idx)%desc = 'surface temperature'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%tsfc(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tg3'
    ExtDiag(idx)%desc = 'deep soil temperature'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%tg3(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tisfc'
    ExtDiag(idx)%desc = 'surface temperature over ice fraction'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%tisfc(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tprcp'
    ExtDiag(idx)%desc = 'total precipitation'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%tprcp(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'vtype'
    ExtDiag(idx)%desc = 'vegetation type in integer 1-13'
    ExtDiag(idx)%unit = 'number'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%vtype(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'weasd'
    ExtDiag(idx)%desc = 'surface snow water equivalent'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%weasd(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'hgtsfc'
    ExtDiag(idx)%desc = 'surface geopotential height'
    ExtDiag(idx)%unit = 'gpm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    ExtDiag(idx)%cnvfac = cn_one
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%oro(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'slmsksfc'
    ExtDiag(idx)%desc = 'sea-land-ice mask (0-sea, 1-land, 2-ice)'
    ExtDiag(idx)%unit = 'numerical'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%slmsk(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'zorlsfc'
    ExtDiag(idx)%desc = 'surface roughness'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    ExtDiag(idx)%cnvfac = cn_one/cn_100
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%zorl(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'vfracsfc'
    ExtDiag(idx)%desc = 'vegetation fraction'
    ExtDiag(idx)%unit = 'fraction'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    ExtDiag(idx)%cnvfac = cn_100
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%vfrac(:)
    enddo

#ifdef CCPP
  if (Model%rdlai) then
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'xlaixy'
    ExtDiag(idx)%desc = 'leaf area index'
    ExtDiag(idx)%unit = 'number'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%xlaixy(:)
    enddo
  endif

  if (Model%lsm == Model%lsm_ruc) then
    do num = 1,Model%lsoil_lsm
      write (xtra,'(i1)') num
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'slc_'//trim(xtra)
      ExtDiag(idx)%desc = 'liquid soil moisture ' // trim(soil_layer_depth(Model%lsm, Model%lsm_ruc, Model%lsm_noah, num))
      ExtDiag(idx)%unit = 'm**3/m**3'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%sh2o(:,num)
      enddo
    enddo
  else
    do num = 1,Model%lsoil_lsm
      write (xtra,'(i1)') num
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'slc_'//trim(xtra)
! DH* Can't use correct unit/description because of the way
! bit for bit tests are conducted (using cmp -> test fails)
#if 0
      ExtDiag(idx)%desc = 'liquid soil moisture ' // trim(soil_layer_depth(Model%lsm, Model%lsm_ruc, Model%lsm_noah, num))
      ExtDiag(idx)%unit = 'm**3/m**3'
#else
      ExtDiag(idx)%desc = 'liquid soil mositure at layer-'//trim(xtra)
      ExtDiag(idx)%unit = 'xxx'
#endif
! *DH
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%slc(:,num)
      enddo
    enddo
  endif
#else
    do num = 1,4
      write (xtra,'(i1)') num
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'slc_'//trim(xtra)
      ExtDiag(idx)%desc = 'liquid soil mositure at layer-'//trim(xtra)
      ExtDiag(idx)%unit = 'xxx'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%slc(:,num)
      enddo
    enddo
#endif

#ifdef CCPP
  if (Model%lsm == Model%lsm_ruc) then
    do num = 1,Model%lsoil_lsm
      write (xtra,'(i1)') num
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'soilw'//trim(xtra)
      ExtDiag(idx)%desc = 'volumetric soil moisture ' // trim(soil_layer_depth(Model%lsm, Model%lsm_ruc, Model%lsm_noah, num))
      ExtDiag(idx)%unit = 'fraction'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%smois(:,num)
      enddo
    enddo
  else
    do num = 1,Model%lsoil_lsm
      write (xtra,'(i1)') num
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'soilw'//trim(xtra)
      ExtDiag(idx)%desc = 'volumetric soil moisture ' // trim(soil_layer_depth(Model%lsm, Model%lsm_ruc, Model%lsm_noah, num))
      ExtDiag(idx)%unit = 'fraction'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%smc(:,num)
      enddo
    enddo
  endif
#else
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'soilw1'
    ExtDiag(idx)%desc = 'volumetric soil moisture 0-10cm'
    ExtDiag(idx)%unit = 'fraction'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%smc(:,1)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'soilw2'
    ExtDiag(idx)%desc = 'volumetric soil moisture 10-40cm'
    ExtDiag(idx)%unit = 'fraction'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%smc(:,2)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'soilw3'
    ExtDiag(idx)%desc = 'volumetric soil moisture 40-100cm'
    ExtDiag(idx)%unit = 'fraction'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%smc(:,3)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'soilw4'
    ExtDiag(idx)%desc = 'volumetric soil moisture 100-200cm'
    ExtDiag(idx)%unit = 'fraction'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%smc(:,4)
    enddo
#endif

#ifdef CCPP
    if (Model%lsm == Model%lsm_ruc) then
      do num = 1,Model%lsoil_lsm
        write (xtra,'(i1)') num
        idx = idx + 1
        ExtDiag(idx)%axes = 2
        ExtDiag(idx)%name = 'soilt'//trim(xtra)
        ExtDiag(idx)%desc = 'soil temperature ' // trim(soil_layer_depth(Model%lsm, Model%lsm_ruc, Model%lsm_noah, num))
        ExtDiag(idx)%unit = 'K'
        ExtDiag(idx)%mod_name = 'gfs_sfc'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%tslb(:,num)
        enddo
      enddo
    else
      do num = 1,Model%lsoil_lsm
        write (xtra,'(i1)') num
        idx = idx + 1
        ExtDiag(idx)%axes = 2
        ExtDiag(idx)%name = 'soilt'//trim(xtra)
        ExtDiag(idx)%desc = 'soil temperature ' // trim(soil_layer_depth(Model%lsm, Model%lsm_ruc, Model%lsm_noah, num))
        ExtDiag(idx)%unit = 'K'
        ExtDiag(idx)%mod_name = 'gfs_sfc'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%stc(:,num)
        enddo
      enddo
    endif
#else
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'soilt1'
    ExtDiag(idx)%desc = 'soil temperature 0-10cm'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%stc(:,1)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'soilt2'
    ExtDiag(idx)%desc = 'soil temperature 10-40cm'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%stc(:,2)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'soilt3'
    ExtDiag(idx)%desc = 'soil temperature 40-100cm'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%stc(:,3)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'soilt4'
    ExtDiag(idx)%desc = 'soil temperature 100-200cm'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%stc(:,4)
    enddo
#endif

!--------------------------nsst variables
  if (model%nstf_name(1) > 0) then
!--------------------------nsst variables

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tref'
    ExtDiag(idx)%desc = 'nsst reference or foundation temperature'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%tref(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'z_c'
    ExtDiag(idx)%desc = 'nsst sub-layer cooling thickness'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%z_c(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'c_0'
    ExtDiag(idx)%desc = 'nsst coefficient1 to calculate d(tz)/d(ts)'
    ExtDiag(idx)%unit = 'numerical'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%c_0(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'c_d'
    ExtDiag(idx)%desc = 'nsst coefficient2 to calculate d(tz)/d(ts)'
    ExtDiag(idx)%unit = 'n/a'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%c_d(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'w_0'
    ExtDiag(idx)%desc = 'nsst coefficient3 to calculate d(tz)/d(ts)'
    ExtDiag(idx)%unit = 'n/a'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%w_0(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'w_d'
    ExtDiag(idx)%desc = 'nsst coefficient4 to calculate d(tz)/d(ts)'
    ExtDiag(idx)%unit = 'n/a'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%w_d(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'xt'
    ExtDiag(idx)%desc = 'nsst heat content in diurnal thermocline layer'
    ExtDiag(idx)%unit = 'k*m'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%xt(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'xs'
    ExtDiag(idx)%desc = 'nsst salinity content in diurnal thermocline layer'
    ExtDiag(idx)%unit = 'n/a'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%xs(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'xu'
    ExtDiag(idx)%desc = 'nsst u-current content in diurnal thermocline layer'
    ExtDiag(idx)%unit = 'm2/s'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%xu(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'xv'
    ExtDiag(idx)%desc = 'nsst v-current content in diurnal thermocline layer'
    ExtDiag(idx)%unit = 'm2/s'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%xv(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'xz'
    ExtDiag(idx)%desc = 'nsst diurnal thermocline layer thickness'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%xz(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'zm'
    ExtDiag(idx)%desc = 'nsst mixed layer thickness'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%zm(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'xtts'
    ExtDiag(idx)%desc = 'nsst d(xt)/d(ts)'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%xtts(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'xzts'
    ExtDiag(idx)%desc = 'nsst d(xt)/d(ts)'
    ExtDiag(idx)%unit = 'm/k'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%xzts(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'd_conv'
    ExtDiag(idx)%desc = 'nsst thickness of free convection layer'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%d_conv(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ifd'
    ExtDiag(idx)%desc = 'nsst index to start dtlm run or not'
    ExtDiag(idx)%unit = 'n/a'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%ifd(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dt_cool'
    ExtDiag(idx)%desc = 'nsst sub-layer cooling amount'
    ExtDiag(idx)%unit = 'k'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%dt_cool(:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'qrain'
    ExtDiag(idx)%desc = 'nsst sensible heat flux due to rainfall'
    ExtDiag(idx)%unit = 'w/m2'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop(nb)%qrain(:)
    enddo
!--------------------------nsst variables
  endif

!--------------------------aerosols
#ifdef CCPP
    if (Model%ntwa>0) then
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'nwfa'
      ExtDiag(idx)%desc = 'number concentration of water-friendly aerosols'
      ExtDiag(idx)%unit = 'kg-1'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Statein(nb)%qgrs(:,:,Model%ntwa)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'nwfa2d'
      ExtDiag(idx)%desc = 'water-friendly surface aerosol source'
      ExtDiag(idx)%unit = 'kg-1 s-1'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling(nb)%nwfa2d
      enddo
    endif

    if (Model%ntia>0) then
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'nifa'
      ExtDiag(idx)%desc = 'number concentration of ice-friendly aerosols'
      ExtDiag(idx)%unit = 'kg-1'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Statein(nb)%qgrs(:,:,Model%ntia)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'nifa2d'
      ExtDiag(idx)%desc = 'ice-friendly surface aerosol source'
      ExtDiag(idx)%unit = 'kg-1 s-1'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling(nb)%nifa2d
      enddo
    endif

    !! Cloud effective radii from Microphysics
    !if (Model%imp_physics == Model%imp_physics_thompson .or. Model%imp_physics == Model%imp_physics_wsm6) then
    !  idx = idx + 1
    !  ExtDiag(idx)%axes = 3
    !  ExtDiag(idx)%name = 'cleffr'
    !  ExtDiag(idx)%desc = 'effective radius of cloud liquid water particle'
    !  ExtDiag(idx)%unit = 'um'
    !  ExtDiag(idx)%mod_name = 'gfs_phys'
    !  allocate (ExtDiag(idx)%data(nblks))
    !  do nb = 1,nblks
    !    ExtDiag(idx)%data(nb)%var3 => Tbd(nb)%phy_f3d(:,:,Model%nleffr)
    !  enddo
    !  idx = idx + 1
    !  ExtDiag(idx)%axes = 3
    !  ExtDiag(idx)%name = 'cieffr'
    !  ExtDiag(idx)%desc = 'effective radius of stratiform cloud ice particle in um'
    !  ExtDiag(idx)%unit = 'um'
    !  ExtDiag(idx)%mod_name = 'gfs_phys'
    !  allocate (ExtDiag(idx)%data(nblks))
    !  do nb = 1,nblks
    !    ExtDiag(idx)%data(nb)%var3 => Tbd(nb)%phy_f3d(:,:,Model%nieffr)
    !  enddo
    !  idx = idx + 1
    !  ExtDiag(idx)%axes = 3
    !  ExtDiag(idx)%name = 'cseffr'
    !  ExtDiag(idx)%desc = 'effective radius of stratiform cloud snow particle in um'
    !  ExtDiag(idx)%unit = 'um'
    !  ExtDiag(idx)%mod_name = 'gfs_phys'
    !  allocate (ExtDiag(idx)%data(nblks))
    !  do nb = 1,nblks
    !    ExtDiag(idx)%data(nb)%var3 => Tbd(nb)%phy_f3d(:,:,Model%nseffr)
    !  enddo
    !endif

    !MYNN
    if (Model%do_mynnedmf) then

      !idx = idx + 1
      !ExtDiag(idx)%axes = 2
      !ExtDiag(idx)%name = 'ktop_shallow'
      !ExtDiag(idx)%desc = 'k-level of plume top'
      !ExtDiag(idx)%unit = 'n/a'
      !ExtDiag(idx)%mod_name = 'gfs_sfc'
      !allocate (ExtDiag(idx)%data(nblks))
      !do nb = 1,nblks
      !  ExtDiag(idx)%data(nb)%var2 => real(IntDiag(nb)%ktop_shallow(:),kind=kind_phys)
      !enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'maxmf'
      ExtDiag(idx)%desc = 'maximum mass-flux in column'
      ExtDiag(idx)%unit = 'm s-1'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%maxmf(:)
      enddo

      !idx = idx + 1
      !ExtDiag(idx)%axes = 2
      !ExtDiag(idx)%name = 'nupdraft'
      !ExtDiag(idx)%desc = 'number of plumes in grid column'
      !ExtDiag(idx)%unit = 'n/a'
      !ExtDiag(idx)%mod_name = 'gfs_sfc'
      !allocate (ExtDiag(idx)%data(nblks))
      !do nb = 1,nblks
      !  ExtDiag(idx)%data(nb)%var2 => real(IntDiag(nb)%nupdraft(:),kind=kind_phys)
      !enddo
    endif

    if (Model%do_mynnsfclay) then
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'zol'
      ExtDiag(idx)%desc = 'monin obukhov surface stability parameter'
      ExtDiag(idx)%unit = 'n/a'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%zol(:)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'flhc'
      ExtDiag(idx)%desc = 'surface exchange coefficient for heat'
      ExtDiag(idx)%unit = 'W m-2 K-1'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%flhc(:)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'flqc'
      ExtDiag(idx)%desc = 'surface exchange coefficient for moisture'
      ExtDiag(idx)%unit = 'kg m-2 s-1'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop(nb)%flqc(:)
      enddo
    endif

    if (Model%do_mynnedmf) then
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'CLDFRA_BL'
      ExtDiag(idx)%desc = 'subgrid cloud fraction'
      ExtDiag(idx)%unit = 'frac'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Tbd(nb)%CLDFRA_BL(:,:)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'QC_BL'
      ExtDiag(idx)%desc = 'subgrid cloud mixing ratio'
      ExtDiag(idx)%unit = 'frac'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Tbd(nb)%QC_BL(:,:)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'EL_PBL'
      ExtDiag(idx)%desc = 'turbulent mixing length'
      ExtDiag(idx)%unit = 'm'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Tbd(nb)%el_pbl(:,:)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'QKE'
      ExtDiag(idx)%desc = '2 X TKE (from mynn)'
      ExtDiag(idx)%unit = 'm2 s-2'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Tbd(nb)%QKE(:,:)
      enddo
    endif
#endif

!  print *,'in gfdl_diag_register,af all extdiag, idx=',idx

! -- chemistry diagnostic variables
  if (Model%cplchm) then

    if (Model%ntchm > 0) then

      if (associated(IntDiag(1)%duem)) then
        do num = 1, size(IntDiag(1)%duem, dim=2)
          idx = idx + 1
          ExtDiag(idx)%axes = 2
          write(ExtDiag(idx)%name,'("duem",i3.3)') num
          write(ExtDiag(idx)%desc,'("Dust Emission Bin ",i0)') num
          ExtDiag(idx)%unit = 'kg/m2/s'
          ExtDiag(idx)%mod_name = 'gfs_phys'
          allocate (ExtDiag(idx)%data(nblks))
          do nb = 1,nblks
            ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%duem(:,num)
          enddo
        enddo
      endif

      if (associated(IntDiag(1)%ssem)) then
        do num = 1, size(IntDiag(1)%ssem, dim=2)
          idx = idx + 1
          ExtDiag(idx)%axes = 2
          write(ExtDiag(idx)%name,'("ssem",i3.3)') num
          write(ExtDiag(idx)%desc,'("Seasalt Emission Bin ",i0)') num
          ExtDiag(idx)%unit = 'kg/m2/s'
          ExtDiag(idx)%mod_name = 'gfs_phys'
          allocate (ExtDiag(idx)%data(nblks))
          do nb = 1,nblks
            ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%ssem(:,num)
          enddo
        enddo
      endif

      if (associated(Model%ntdiag)) then
        idt = 0
        do num = Model%ntchs, Model%ntchm + Model%ntchs - 1
          if (Model%ntdiag(num-Model%ntchs+1)) then
            idt = idt + 1
            idx = idx + 1
            ExtDiag(idx)%axes = 2
            ExtDiag(idx)%name = trim(Model%tracer_names(num)) // 'sd'
            ExtDiag(idx)%desc = trim(Model%tracer_names(num)) // ' Sedimentation'
            ExtDiag(idx)%unit = 'kg/m2/s'
            ExtDiag(idx)%mod_name = 'gfs_phys'
            allocate (ExtDiag(idx)%data(nblks))
            do nb = 1,nblks
              ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%sedim(:,idt)
            enddo

            idx = idx + 1
            ExtDiag(idx)%axes = 2
            ExtDiag(idx)%name = trim(Model%tracer_names(num)) // 'dp'
            ExtDiag(idx)%desc = trim(Model%tracer_names(num)) // ' Dry Deposition'
            ExtDiag(idx)%unit = 'kg/m2/s'
            ExtDiag(idx)%mod_name = 'gfs_phys'
            allocate (ExtDiag(idx)%data(nblks))
            do nb = 1,nblks
              ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%drydep(:,idt)
            enddo

            idx = idx + 1
            ExtDiag(idx)%axes = 2
            ExtDiag(idx)%name = trim(Model%tracer_names(num)) // 'wtl'
            ExtDiag(idx)%desc = trim(Model%tracer_names(num)) // ' Large-Scale Wet Deposition'
            ExtDiag(idx)%unit = 'kg/m2/s'
            ExtDiag(idx)%mod_name = 'gfs_phys'
            allocate (ExtDiag(idx)%data(nblks))
            do nb = 1,nblks
              ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%wetdpl(:,idt)
            enddo

            idx = idx + 1
            ExtDiag(idx)%axes = 2
            ExtDiag(idx)%name = trim(Model%tracer_names(num)) // 'wtc'
            ExtDiag(idx)%desc = trim(Model%tracer_names(num)) // ' Convective-Scale Wet Deposition'
            ExtDiag(idx)%unit = 'kg/m2/s'
            ExtDiag(idx)%mod_name = 'gfs_phys'
            allocate (ExtDiag(idx)%data(nblks))
            do nb = 1,nblks
              ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%wetdpc(:,idt)
            enddo
          endif
        enddo
      endif

    endif

    num = size(IntDiag(1)%abem, dim=2)
    do num = 1, size(IntDiag(1)%abem, dim=2)
      idx = idx + 1
      select case (mod(num,3))
        case (0)
          ExtDiag(idx)%name = 'bcem'
          ExtDiag(idx)%desc = 'Black Carbon'
        case (1)
          ExtDiag(idx)%name = 'ocem'
          ExtDiag(idx)%desc = 'Organic Carbon'
        case (2)
          ExtDiag(idx)%name = 'so2em'
          ExtDiag(idx)%desc = 'SO2'
      end select

      if (num > 3) then
        ExtDiag(idx)%name = trim(ExtDiag(idx)%name) // 'bb'
        ExtDiag(idx)%desc = trim(ExtDiag(idx)%desc) // ' Biomass Burning Emissions'
      else
        ExtDiag(idx)%name = trim(ExtDiag(idx)%name) // 'an'
        ExtDiag(idx)%desc = trim(ExtDiag(idx)%desc) // ' Anthropogenic Emissions'
      end if

      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%unit = 'ug/m2/s'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%abem(:,num)
      enddo
    end do

    do num = 1, size(IntDiag(1)%aecm, dim=2)
      idx = idx + 1
      select case (num)
        case(1)
          ExtDiag(idx)%name = 'aecmass'
          ExtDiag(idx)%desc = 'Aerosol Column Mass Density (PM2.5)'
        case(2)
          ExtDiag(idx)%name = 'bccmass'
          ExtDiag(idx)%desc = 'Black Carbon Column Mass Density'
        case(3)
          ExtDiag(idx)%name = 'occmass'
          ExtDiag(idx)%desc = 'Organic Carbon Column Mass Density'
        case(4)
          ExtDiag(idx)%name = 'sucmass'
          ExtDiag(idx)%desc = 'Sulfate Column Mass Density'
        case(5)
          ExtDiag(idx)%name = 'ducmass'
          ExtDiag(idx)%desc = 'Dust Column Mass Density'
        case(6)
          ExtDiag(idx)%name = 'sscmass'
          ExtDiag(idx)%desc = 'Seasalt Column Mass Density'
      end select

      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%unit = 'g/m2'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => IntDiag(nb)%aecm(:,num)
      enddo
    end do

  endif

!--- prognostic variable tendencies (t, u, v, sph, clwmr, o3)
!rab    idx = idx + 1
!rab    ExtDiag(idx)%axes = 3
!rab    ExtDiag(idx)%name = 'dtemp_dt'
!rab    ExtDiag(idx)%desc = 'gfs radiation/physics temperature tendency'
!rab    ExtDiag(idx)%unit = 'k/s'
!rab    ExtDiag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    ExtDiag(idx)%axes = 3
!rab    ExtDiag(idx)%name = 'du_dt'
!rab    ExtDiag(idx)%desc = 'gfs radiation/physics horizontal wind component tendency'
!rab    ExtDiag(idx)%unit = 'm/s/s'
!rab    ExtDiag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    ExtDiag(idx)%axes = 3
!rab    ExtDiag(idx)%name = 'dv_dt'
!rab    ExtDiag(idx)%desc = 'gfs radiation/physics meridional wind component tendency'
!rab    ExtDiag(idx)%unit = 'm/s/s'
!rab    ExtDiag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    ExtDiag(idx)%axes = 3
!rab    ExtDiag(idx)%name = 'dsphum_dt'
!rab    ExtDiag(idx)%desc = 'gfs radiation/physics specific humidity tendency'
!rab    ExtDiag(idx)%unit = 'kg/kg/s'
!rab    ExtDiag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    ExtDiag(idx)%axes = 3
!rab    ExtDiag(idx)%name = 'dclwmr_dt'
!rab    ExtDiag(idx)%desc = 'gfs radiation/radiation cloud water mixing ratio tendency'
!rab    ExtDiag(idx)%unit = 'kg/kg/s'
!rab    ExtDiag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    ExtDiag(idx)%axes = 3
!rab    ExtDiag(idx)%name = 'do3mr_dt'
!rab    ExtDiag(idx)%desc = 'gfs radiation/radiation ozone mixing ratio tendency'
!rab    ExtDiag(idx)%unit = 'kg/kg/s'
!rab    ExtDiag(idx)%mod_name = 'gfs_phys'

  end subroutine GFS_externaldiag_populate

#ifdef CCPP
  function soil_layer_depth(lsm, lsm_ruc, lsm_noah, layer) result(layer_depth)
     character(len=30)   :: layer_depth
     integer, intent(in) :: lsm, lsm_ruc, lsm_noah, layer
     !
     continue
     !
     if (lsm==lsm_ruc) then
        select case (layer)
           case (1)
              layer_depth = 'at 0 cm depth'
           case (2)
              layer_depth = 'at 5 cm depth'
           case (3)
              layer_depth = 'at 20 cm depth'
           case (4)
              layer_depth = 'at 40 cm depth'
           case (5)
              layer_depth = 'at 60 cm depth'
           case (6)
              layer_depth = 'at 100 cm depth'
           case (7)
              layer_depth = 'at 160 cm depth'
           case (8)
              layer_depth = 'at 220 cm depth'
           case (9)
              layer_depth = 'at 300 cm depth'
           case default
              write (layer_depth,'(a,i0)') 'invalid layer ', layer
        end select
     else if (lsm==lsm_noah) then
        select case (layer)
           case (1)
              layer_depth = '0-10cm'
           case (2)
              layer_depth = '10-40cm'
           case (3)
              layer_depth = '40-100cm'
           case (4)
              layer_depth = '100-200cm'
           case default
              write (layer_depth,'(a,i0)') 'invalid layer ', layer
        end select
     else
        write (layer_depth,'(a,i0)') 'unknown layer ', layer
     end if
     !
     return
     !
  end function soil_layer_depth
#endif

!-------------------------------------------------------------------------      

end module GFS_diagnostics
