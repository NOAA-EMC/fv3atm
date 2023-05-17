module FV3GFS_oro_io

  use block_control_mod,  only: block_control_type
  use fms2_io_mod,        only: FmsNetcdfDomainFile_t, unlimited, &
       register_axis, register_restart_field,       &
       register_variable_attribute, register_field, &
       read_restart, write_restart, write_data,     &
       get_global_io_domain_indices, variable_exists
  use FV3GFS_common_io,   only: get_nx_ny_from_atm
  use GFS_typedefs,       only: GFS_sfcprop_type, GFS_control_type, kind_phys
  use GFS_restart,        only: GFS_restart_type
  use mpp_mod,            only: mpp_error,  mpp_pe, mpp_root_pe, &
       mpp_chksum, NOTE,   FATAL
  use physcons,           only: con_tice          !saltwater freezing temp (K)

  implicit none
  private

  public :: Oro_io_data_type, Oro_io_register, Oro_io_copy, Oro_io_final

  type Oro_io_data_type
    character(len=32),    pointer, private, dimension(:)       :: name2 => null()
    real(kind=kind_phys), pointer, private, dimension(:,:,:)   :: var2  => null()
    real(kind=kind_phys), pointer, private, dimension(:,:,:)   :: var3v => null()
    real(kind=kind_phys), pointer, private, dimension(:,:,:)   :: var3s => null()
  contains
    procedure, public :: register => Oro_io_register
    procedure, public :: copy => Oro_io_copy
    final :: Oro_io_final
  end type Oro_io_data_type

contains

  subroutine Oro_io_register(oro, Model, Oro_restart, Atm_block)
    implicit none
    class(Oro_io_data_type) :: oro
    type(GFS_control_type),      intent(in) :: Model
    type(FmsNetcdfDomainFile_t) :: Oro_restart
    type(block_control_type), intent(in) :: Atm_block

    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p  => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_fr => NULL()
    integer :: nx, ny

    integer :: nvar_o2, nvar_vegfr, nvar_soilfr, n, num

    ! This #define reduces code length by a lot
#define WARN_DISASSOCIATE(name) \
    if(associated(name)) then ; \
      write(0,*) 'Internal error. Called oro%copy twice. Will try to keep going anyway.' ; \
      deallocate(name); \
      nullify(name) ; \
    endif

    call get_nx_ny_from_atm(Atm_block, nx, ny)

    WARN_DISASSOCIATE(oro%name2)
    WARN_DISASSOCIATE(oro%var2)
    WARN_DISASSOCIATE(oro%var3v)
    WARN_DISASSOCIATE(oro%var3s)
#undef WARN_DISASSOCIATE

    nvar_o2  = 19
    nvar_vegfr  = Model%nvegcat
    nvar_soilfr = Model%nsoilcat

    allocate(oro%name2(nvar_o2))
    allocate(oro%var2(nx,ny,nvar_o2))

    allocate(oro%var3v(nx,ny,nvar_vegfr))
    allocate(oro%var3s(nx,ny,nvar_soilfr))

    oro%var2 = -9999._kind_phys

    num = 1       ; oro%name2(num)  = 'stddev'     ! hprime(ix,1)
    num = num + 1 ; oro%name2(num)  = 'convexity'  ! hprime(ix,2)
    num = num + 1 ; oro%name2(num)  = 'oa1'        ! hprime(ix,3)
    num = num + 1 ; oro%name2(num)  = 'oa2'        ! hprime(ix,4)
    num = num + 1 ; oro%name2(num)  = 'oa3'        ! hprime(ix,5)
    num = num + 1 ; oro%name2(num)  = 'oa4'        ! hprime(ix,6)
    num = num + 1 ; oro%name2(num)  = 'ol1'        ! hprime(ix,7)
    num = num + 1 ; oro%name2(num)  = 'ol2'        ! hprime(ix,8)
    num = num + 1 ; oro%name2(num)  = 'ol3'        ! hprime(ix,9)
    num = num + 1 ; oro%name2(num) = 'ol4'        ! hprime(ix,10)
    num = num + 1 ; oro%name2(num) = 'theta'      ! hprime(ix,11)
    num = num + 1 ; oro%name2(num) = 'gamma'      ! hprime(ix,12)
    num = num + 1 ; oro%name2(num) = 'sigma'      ! hprime(ix,13)
    num = num + 1 ; oro%name2(num) = 'elvmax'     ! hprime(ix,14)
    num = num + 1 ; oro%name2(num) = 'orog_filt'  ! oro
    num = num + 1 ; oro%name2(num) = 'orog_raw'   ! oro_uf
    num = num + 1 ; oro%name2(num) = 'land_frac'  ! land fraction [0:1]
    !--- variables below here are optional
    num = num + 1 ; oro%name2(num) = 'lake_frac'  ! lake fraction [0:1]
    num = num + 1 ; oro%name2(num) = 'lake_depth' ! lake depth(m)

    !--- register axis
    call register_axis( Oro_restart, "lon", 'X' )
    call register_axis( Oro_restart, "lat", 'Y' )
    !--- register the 2D fields
    do n = 1,num
      var2_p => oro%var2(:,:,n)
      if (trim(oro%name2(n)) == 'lake_frac' .or. trim(oro%name2(n)) == 'lake_depth' ) then
        call register_restart_field(Oro_restart, oro%name2(n), var2_p, dimensions=(/'lat','lon'/), is_optional=.true.)
      else
        call register_restart_field(Oro_restart, oro%name2(n), var2_p, dimensions=(/'lat','lon'/))
      endif
    enddo

    !--- register 3D vegetation and soil fractions
    var3_fr => oro%var3v(:,:,:)
    call register_restart_field(Oro_restart, 'vegetation_type_pct', var3_fr, dimensions=(/'num_veg_cat','lat        ','lon        '/) , is_optional=.true.)
    var3_fr => oro%var3s(:,:,:)
    call register_restart_field(Oro_restart, 'soil_type_pct', var3_fr, dimensions=(/'num_soil_cat','lat         ','lon         '/) , is_optional=.true.)

  end subroutine Oro_io_register

  subroutine Oro_io_copy(oro, Model, Sfcprop, Atm_block)
    implicit none
    class(Oro_io_data_type) :: oro
    type(GFS_control_type),      intent(in) :: Model
    type(GFS_sfcprop_type)                  :: Sfcprop(:)
    type(FmsNetcdfDomainFile_t) :: Oro_restart
    type(block_control_type), intent(in) :: Atm_block

    integer :: i,j,nb,ix,num

    !$omp parallel do default(shared) private(i, j, nb, ix, num)
    do nb = 1, Atm_block%nblks
      !--- 2D variables
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - Atm_block%isc + 1
        j = Atm_block%index(nb)%jj(ix) - Atm_block%jsc + 1
        !--- stddev
        !       Sfcprop(nb)%hprim(ix)     = oro%var2(i,j,1)
        !--- hprime(1:14)
        num = 1       ; Sfcprop(nb)%hprime(ix,num)  = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num)  = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num)  = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num)  = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num)  = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num)  = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num)  = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num)  = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num)  = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num) = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num) = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num) = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num) = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%hprime(ix,num) = oro%var2(i,j,num)
        !--- oro
        num = num + 1 ; Sfcprop(nb)%oro(ix)       = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop(nb)%oro_uf(ix)    = oro%var2(i,j,num)

        Sfcprop(nb)%landfrac(ix)  = -9999.0
        Sfcprop(nb)%lakefrac(ix)  = -9999.0

        num = num + 1 ; Sfcprop(nb)%landfrac(ix)  = oro%var2(i,j,num) !land frac [0:1]
        if (Model%lkm > 0  ) then
          if(oro%var2(i,j,num+1)>Model%lakefrac_threshold .and. &
               oro%var2(i,j,num+2)>Model%lakedepth_threshold) then
            Sfcprop(nb)%lakefrac(ix)  = oro%var2(i,j,num+1) !lake frac [0:1]
            Sfcprop(nb)%lakedepth(ix) = oro%var2(i,j,num+2) !lake depth [m]    !YWu
          else
            Sfcprop(nb)%lakefrac(ix)  = 0
            Sfcprop(nb)%lakedepth(ix) = -9999
          endif
        else
          Sfcprop(nb)%lakefrac(ix)  = oro%var2(i,j,num+1) !lake frac [0:1]
          Sfcprop(nb)%lakedepth(ix) = oro%var2(i,j,num+2) !lake depth [m]    !YWu
        endif
        num = num + 2 ! To account for lakefrac and lakedepth

        Sfcprop(nb)%vegtype_frac(ix,:)  =  -9999.0
        Sfcprop(nb)%soiltype_frac(ix,:) =  -9999.0

        Sfcprop(nb)%vegtype_frac(ix,:)  = oro%var3v(i,j,:) ! vegetation type fractions, [0:1]
        Sfcprop(nb)%soiltype_frac(ix,:) = oro%var3s(i,j,:) ! soil type fractions, [0:1]

      enddo
    enddo

    !--- deallocate containers and free restart container
    deallocate(oro%name2)
    deallocate(oro%var2)
    deallocate(oro%var3v)
    deallocate(oro%var3s)

    nullify(oro%name2)
    nullify(oro%var2)
    nullify(oro%var3v)
    nullify(oro%var3s)

  end subroutine Oro_io_copy

  subroutine Oro_io_final(oro)
    implicit none
    type(Oro_io_data_type) :: oro
    
    ! This #define reduces code length by a lot
#define IF_ASSOC_DEALLOC_NULL(var) \
    if(associated(oro%var)) then ; \
      deallocate(oro%var) ; \
      nullify(oro%var) ; \
    endif

    IF_ASSOC_DEALLOC_NULL(name2)
    IF_ASSOC_DEALLOC_NULL(var2)
    IF_ASSOC_DEALLOC_NULL(var3s)
    IF_ASSOC_DEALLOC_NULL(var3v)

#undef IF_ASSOC_DEALLOC_NULL
  end subroutine Oro_io_final

end module FV3GFS_oro_io
