module CCPP_typedefs

!> \section arg_table_CCPP_typedefs Argument Table
!! \htmlinclude CCPP_typedefs.html
!!

    use machine, only: kind_grid, kind_dyn

    implicit none

    private

    public CCPP_interstitial_type

!! \section arg_table_CCPP_interstitial_type Argument Table
!! \htmlinclude CCPP_interstitial_type.html
!!
  type CCPP_interstitial_type

     real(kind_dyn)                      :: akap
     real(kind_dyn)                      :: bdt
     real(kind_dyn), pointer             :: cappa(:,:,:)
     logical                             :: do_qa
     real(kind_dyn), pointer             :: dtdt(:,:,:)
     logical                             :: fast_mp_consv
     integer                             :: kmp
     logical                             :: last_step
     real(kind_dyn)                      :: mdt
     logical                             :: out_dt
     real(kind_dyn), pointer             :: pfull(:)
     real(kind_dyn), pointer             :: te0_2d(:,:) ! called te_2d in fv_dynamics, te0_2d in Lagrangian_to_Eulerian, te0_2d in fv_sat_adj
     real(kind_dyn), pointer             :: te0(:,:,:)  ! called dp1 in fv_dynamics, te in Lagrangian_to_Eulerian, te0 in fv_sat_adj
     real(kind_dyn)                      :: zvir
     logical                             :: do_sat_adj
     integer                             :: is
     integer                             :: ie
     integer                             :: isd
     integer                             :: ied
     integer                             :: js
     integer                             :: je
     integer                             :: jsd
     integer                             :: jed
     integer                             :: ng
     integer                             :: npz
     integer                             :: npzp1
     integer                             :: npzdelz
     integer                             :: npzq_con
     integer                             :: npzcappa
     real(kind_dyn),  pointer            :: delp(:,:,:)
     real(kind_dyn),  pointer            :: delz(:,:,:)
     real(kind_grid), pointer            :: area(:,:)
     real(kind_dyn),  pointer            :: peln(:,:,:)
     real(kind_dyn),  pointer            :: phis(:,:)
     real(kind_dyn),  pointer            :: pkz(:,:,:)
     real(kind_dyn),  pointer            :: pt(:,:,:)
     real(kind_dyn),  pointer            :: qvi(:,:,:,:)
     real(kind_dyn),  pointer            :: qv(:,:,:)
     real(kind_dyn),  pointer            :: ql(:,:,:)
     real(kind_dyn),  pointer            :: qi(:,:,:)
     real(kind_dyn),  pointer            :: qr(:,:,:)
     real(kind_dyn),  pointer            :: qs(:,:,:)
     real(kind_dyn),  pointer            :: qg(:,:,:)
     real(kind_dyn),  pointer            :: qc(:,:,:)
     real(kind_dyn),  pointer            :: q_con(:,:,:)
     integer                             :: nthreads
     logical                             :: hydrostatic
     integer                             :: nwat
     integer                             :: ngas
     real(kind_dyn),  pointer            :: rilist(:)
     real(kind_dyn),  pointer            :: cpilist(:)
     integer                             :: mpirank
     integer                             :: mpiroot

  contains

    procedure :: create  => interstitial_create     !<   allocate array data
    procedure :: reset   => interstitial_reset      !<   reset array data
    procedure :: mprint  => interstitial_print      !<   print array data

  end type CCPP_interstitial_type

contains

!-----------------------------
! CCPP_interstitial_type
!-----------------------------
  subroutine interstitial_create (Interstitial, is, ie, isd, ied, js, je, jsd, jed, npz, ng, &
                                  dt_atmos, p_split, k_split, zvir, p_ref, ak, bk,           &
                                  do_ql, do_qi, do_qr, do_qs, do_qg, do_qa,                  &
                                  kappa, hydrostatic, do_sat_adj,                            &
                                  delp, delz, area, peln, phis, pkz, pt,                     &
                                  qvi, qv, ql, qi, qr, qs, qg, qc, q_con,                    &
                                  nthreads, nwat, ngas, rilist, cpilist, mpirank, mpiroot)
    !
    implicit none
    !
    class(CCPP_interstitial_type) :: Interstitial
    integer, intent(in) :: is
    integer, intent(in) :: ie
    integer, intent(in) :: isd
    integer, intent(in) :: ied
    integer, intent(in) :: js
    integer, intent(in) :: je
    integer, intent(in) :: jsd
    integer, intent(in) :: jed
    integer, intent(in) :: npz
    integer, intent(in) :: ng
    real(kind_dyn),    intent(in) :: dt_atmos
    integer, intent(in) :: p_split
    integer, intent(in) :: k_split
    real(kind_dyn),    intent(in) :: zvir
    real(kind_dyn),    intent(in) :: p_ref
    real(kind_dyn),    intent(in) :: ak(:)
    real(kind_dyn),    intent(in) :: bk(:)
    logical, intent(in) :: do_ql
    logical, intent(in) :: do_qi
    logical, intent(in) :: do_qr
    logical, intent(in) :: do_qs
    logical, intent(in) :: do_qg
    logical, intent(in) :: do_qa
    real(kind_dyn),    intent(in) :: kappa
    logical, intent(in) :: hydrostatic
    logical, intent(in) :: do_sat_adj
    real(kind_dyn),  target, intent(in) :: delp(:,:,:)
    real(kind_dyn),  target, intent(in) :: delz(:,:,:)
    real(kind_grid), target, intent(in) :: area(:,:)
    real(kind_dyn),  target, intent(in) :: peln(:,:,:)
    real(kind_dyn),  target, intent(in) :: phis(:,:)
    real(kind_dyn),  target, intent(in) :: pkz(:,:,:)
    real(kind_dyn),  target, intent(in) :: pt(:,:,:)
    real(kind_dyn),  target, intent(in) :: qvi(:,:,:,:)
    real(kind_dyn),  target, intent(in) :: qv(:,:,:)
    real(kind_dyn),  target, intent(in) :: ql(:,:,:)
    real(kind_dyn),  target, intent(in) :: qi(:,:,:)
    real(kind_dyn),  target, intent(in) :: qr(:,:,:)
    real(kind_dyn),  target, intent(in) :: qs(:,:,:)
    real(kind_dyn),  target, intent(in) :: qg(:,:,:)
    real(kind_dyn),  target, intent(in) :: qc(:,:,:)
    real(kind_dyn),  target, intent(in) :: q_con(:,:,:)
    integer, intent(in) :: nthreads
    ! For multi-gases physics
    integer,        intent(in)           :: nwat
    integer,        intent(in), optional :: ngas
    real(kind_dyn), intent(in), optional :: rilist(:)
    real(kind_dyn), intent(in), optional :: cpilist(:)
    integer,        intent(in)           :: mpirank
    integer,        intent(in)           :: mpiroot
    !
#ifdef MOIST_CAPPA
    Interstitial%npzcappa = npz
    allocate (Interstitial%cappa  (isd:ied, jsd:jed, 1:npz) )
#else
    Interstitial%npzcappa = 1
    allocate (Interstitial%cappa  (isd:isd, jsd:jsd, 1)     )
#endif
    allocate (Interstitial%dtdt   (is:ie, js:je, 1:npz)     )
    allocate (Interstitial%pfull  (1:npz)                   )
    allocate (Interstitial%te0_2d (is:ie, js:je)            )
    allocate (Interstitial%te0    (isd:ied, jsd:jed, 1:npz) )
    !
    ! Initialize variables to default values
#ifdef SW_DYNAMICS
    Interstitial%akap = 1.
#else
    Interstitial%akap = kappa
#endif
    Interstitial%bdt       = dt_atmos/real(abs(p_split))
    Interstitial%do_qa     = do_qa
    Interstitial%mdt       = Interstitial%bdt/real(k_split)
    !
    ! Flag for hydrostatic/non-hydrostatic physics
    Interstitial%hydrostatic = hydrostatic
    if (hydrostatic) then
       Interstitial%npzdelz = 1
    else
       Interstitial%npzdelz = npz
    end if
    !
    Interstitial%zvir      = zvir
    !
    Interstitial%do_sat_adj =  do_sat_adj
    Interstitial%is         =  is
    Interstitial%ie         =  ie
    Interstitial%isd        =  isd
    Interstitial%ied        =  ied
    Interstitial%js         =  js
    Interstitial%je         =  je
    Interstitial%jsd        =  jsd
    Interstitial%jed        =  jed
    Interstitial%ng         =  ng
    Interstitial%npz        =  npz
    Interstitial%npzp1      =  npz+1
    ! Set up links from CCPP_interstitial DDT to ATM DDT
    Interstitial%delp       => delp
    Interstitial%delz       => delz
    Interstitial%area       => area
    Interstitial%peln       => peln
    Interstitial%phis       => phis
    Interstitial%pkz        => pkz
    Interstitial%pt         => pt
    Interstitial%qvi        => qvi
    Interstitial%qv         => qv
    if (do_ql) Interstitial%ql => ql
    if (do_qi) Interstitial%qi => qi
    if (do_qr) Interstitial%qr => qr
    if (do_qs) Interstitial%qs => qs
    if (do_qg) Interstitial%qg => qg
    if (do_qa) Interstitial%qc => qc
#ifdef USE_COND
    Interstitial%npzq_con = npz
#else
    Interstitial%npzq_con = 1
#endif
    Interstitial%q_con      => q_con
    !
    ! Number of OpenMP threads available for schemes
    Interstitial%nthreads   = nthreads
    !
    ! For multi-gases physics
    Interstitial%nwat  = nwat
    ! If ngas, rilist and cpilist are present, then
    ! multi-gases physics are used. If not, set ngas=0
    ! (safe value), allocate rilist/cpilist and set to zero
    if(present(ngas)) then
      Interstitial%ngas  = ngas
    else
      Interstitial%ngas  = 0
    end if
    allocate(Interstitial%rilist(0:Interstitial%ngas))
    allocate(Interstitial%cpilist(0:Interstitial%ngas))
    if (present(rilist)) then
      Interstitial%rilist  = rilist(0:Interstitial%ngas)
      Interstitial%cpilist = cpilist(0:Interstitial%ngas)
    else
      Interstitial%rilist  = 0.0
      Interstitial%cpilist = 0.0
    end if
    !
    Interstitial%mpirank = mpirank
    Interstitial%mpiroot = mpiroot
    !
    ! Calculate vertical pressure levels
    call interstitital_calculate_pressure_levels(Interstitial, npz, p_ref, ak, bk)
    !
    ! Reset all other variables
    call Interstitial%reset()
    !
  end subroutine interstitial_create

  subroutine interstitital_calculate_pressure_levels(Interstitial, npz, p_ref, ak, bk)

     implicit none

     class(CCPP_interstitial_type)       :: Interstitial
     integer, intent(in) :: npz
     real(kind_dyn), intent(in) :: p_ref
     real(kind_dyn), intent(in) :: ak(:)
     real(kind_dyn), intent(in) :: bk(:)

     real(kind_dyn) :: ph1
     real(kind_dyn) :: ph2
     integer :: k

#ifdef SW_DYNAMICS
      Interstitial%pfull(1) = 0.5*p_ref
#else
      do k=1,npz
         ph1 = ak(k  ) + bk(k  )*p_ref
         ph2 = ak(k+1) + bk(k+1)*p_ref
         Interstitial%pfull(k) = (ph2 - ph1) / log(ph2/ph1)
      enddo
#endif
      ! DH* This is copied from fv_mapz.F90, does it work with SW_DYNAMICS?
      do k=1,npz
         Interstitial%kmp = k
         if ( Interstitial%pfull(k) > 10.E2 ) exit
      enddo
  end subroutine interstitital_calculate_pressure_levels

  subroutine interstitial_reset (Interstitial)
    !
    implicit none
    !
    class(CCPP_interstitial_type) :: Interstitial
    !
    Interstitial%cappa         = 0.0
    Interstitial%dtdt          = 0.0
    Interstitial%fast_mp_consv = .false.
    Interstitial%last_step     = .false.
    Interstitial%out_dt        = .false.
    Interstitial%te0_2d        = 0.0
    Interstitial%te0           = 0.0
    !
  end subroutine interstitial_reset

  subroutine interstitial_print(Interstitial)
    !
    implicit none
    !
    class(CCPP_interstitial_type) :: Interstitial
    !
    ! Print static variables
    write (0,'(a)') 'Interstitial_print'
    write (0,*) 'Interstitial_print: values that do not change'
    write (0,*) 'Interstitial%akap              = ', Interstitial%akap
    write (0,*) 'Interstitial%bdt               = ', Interstitial%bdt
    write (0,*) 'Interstitial%kmp               = ', Interstitial%kmp
    write (0,*) 'Interstitial%mdt               = ', Interstitial%mdt
    write (0,*) 'sum(Interstitial%pfull)        = ', sum(Interstitial%pfull)
    write (0,*) 'Interstitial%zvir              = ', Interstitial%zvir
    write (0,*) 'Interstitial%do_qa             = ', Interstitial%do_qa
    write (0,*) 'Interstitial%do_sat_adj        = ', Interstitial%do_sat_adj
    write (0,*) 'Interstitial%is                = ', Interstitial%is
    write (0,*) 'Interstitial%ie                = ', Interstitial%ie
    write (0,*) 'Interstitial%isd               = ', Interstitial%isd
    write (0,*) 'Interstitial%ied               = ', Interstitial%ied
    write (0,*) 'Interstitial%js                = ', Interstitial%js
    write (0,*) 'Interstitial%je                = ', Interstitial%je
    write (0,*) 'Interstitial%jsd               = ', Interstitial%jsd
    write (0,*) 'Interstitial%jed               = ', Interstitial%jed
    write (0,*) 'sum(Interstitial%area)         = ', sum(Interstitial%area)
    write (0,*) 'Interstitial%ng                = ', Interstitial%ng
    write (0,*) 'Interstitial%npz               = ', Interstitial%npz
    write (0,*) 'Interstitial%npzp1             = ', Interstitial%npzp1
    write (0,*) 'Interstitial%npzdelz           = ', Interstitial%npzdelz
    write (0,*) 'Interstitial%npzq_con          = ', Interstitial%npzq_con
    write (0,*) 'Interstitial%npzcappa          = ', Interstitial%npzcappa
    ! Print all other variables
    write (0,*) 'Interstitial_print: values that change'
    write (0,*) 'sum(Interstitial%cappa)        = ', sum(Interstitial%cappa)
    write (0,*) 'sum(Interstitial%dtdt)         = ', sum(Interstitial%dtdt)
    write (0,*) 'Interstitial%fast_mp_consv     = ', Interstitial%fast_mp_consv
    write (0,*) 'Interstitial%last_step         = ', Interstitial%last_step
    write (0,*) 'Interstitial%out_dt            = ', Interstitial%out_dt
    write (0,*) 'sum(Interstitial%te0_2d)       = ', sum(Interstitial%te0_2d)
    write (0,*) 'sum(Interstitial%te0)          = ', sum(Interstitial%te0)
    write (0,*) 'sum(Interstitial%delp)         = ', Interstitial%delp
    write (0,*) 'sum(Interstitial%delz)         = ', Interstitial%delz
    write (0,*) 'sum(Interstitial%peln)         = ', Interstitial%peln
    write (0,*) 'sum(Interstitial%phis)         = ', Interstitial%phis
    write (0,*) 'sum(Interstitial%pkz)          = ', Interstitial%pkz
    write (0,*) 'sum(Interstitial%pt)           = ', Interstitial%pt
    write (0,*) 'sum(Interstitial%qvi)          = ', Interstitial%qvi
    write (0,*) 'sum(Interstitial%qv)           = ', Interstitial%qv
    if (associated(Interstitial%ql)) write (0,*) 'sum(Interstitial%ql)           = ', Interstitial%ql
    if (associated(Interstitial%qi)) write (0,*) 'sum(Interstitial%qi)           = ', Interstitial%qi
    if (associated(Interstitial%qr)) write (0,*) 'sum(Interstitial%qr)           = ', Interstitial%qr
    if (associated(Interstitial%qs)) write (0,*) 'sum(Interstitial%qs)           = ', Interstitial%qs
    if (associated(Interstitial%qg)) write (0,*) 'sum(Interstitial%qg)           = ', Interstitial%qg
    if (associated(Interstitial%qc)) write (0,*) 'sum(Interstitial%qc)           = ', Interstitial%qc
    write (0,*) 'sum(Interstitial%q_con)        = ', Interstitial%q_con
    write (0,*) 'Interstitial%hydrostatic       = ', Interstitial%hydrostatic
    write (0,*) 'Interstitial%nwat              = ', Interstitial%nwat
    write (0,*) 'Interstitial%ngas              = ', Interstitial%ngas
    write (0,*) 'Interstitial%rilist            = ', Interstitial%rilist
    write (0,*) 'Interstitial%cpilist           = ', Interstitial%cpilist
    write (0,*) 'Interstitial%mpirank           = ', Interstitial%mpirank
    write (0,*) 'Interstitial%mpiroot           = ', Interstitial%mpiroot
    write (0,*) 'Interstitial%nthreads          = ', Interstitial%nthreads
    write (0,*) 'Interstitial_print: end'
    !
  end subroutine interstitial_print

end module CCPP_typedefs

