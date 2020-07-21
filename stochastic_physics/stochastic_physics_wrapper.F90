module stochastic_physics_wrapper_mod

  use machine, only: kind_phys

  implicit none

  ! For stochastic physics pattern generation
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: xlat
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: xlon
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: sppt_wts
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: shum_wts
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: skebu_wts
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: skebv_wts
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: sfc_wts

  ! For cellular automata
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: ugrs
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: qgrs
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: pgr
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: vvl
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: prsl
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: condition
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: ca_deep_cpl, ca_turb_cpl, ca_shal_cpl
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: ca_deep_diag,ca_turb_diag,ca_shal_diag
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: ca1_cpl, ca2_cpl, ca3_cpl
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: ca1_diag,ca2_diag,ca3_diag


!----------------
! Public Entities
!----------------
! functions
  public stochastic_physics_wrapper

  contains

!-------------------------------
!  CCPP step
!-------------------------------
  subroutine stochastic_physics_wrapper (GFS_Control, GFS_Data, Atm_block)

#ifdef OPENMP
    use omp_lib
#endif

    use GFS_typedefs,       only: GFS_control_type, GFS_data_type
    use mpp_mod,            only: FATAL, mpp_error
    use block_control_mod,  only: block_control_type
    use atmosphere_mod,     only: Atm, mygrid

    use stochastic_physics,           only: init_stochastic_physics, run_stochastic_physics
    use stochastic_physics_sfc,       only: run_stochastic_physics_sfc
    use cellular_automata_global_mod, only: cellular_automata_global
    use cellular_automata_sgs_mod,    only: cellular_automata_sgs

    implicit none

    type(GFS_control_type),   intent(inout) :: GFS_Control
    type(GFS_data_type),      intent(inout) :: GFS_Data(:)
    type(block_control_type), intent(inout) :: Atm_block

    integer :: nthreads, nb

#ifdef OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

    ! Initialize
    initalize_stochastic_physics: if (GFS_Control%kdt==0) then

      if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb .OR. GFS_Control%do_sfcperts) then
        ! Initialize stochastic physics
        call init_stochastic_physics(GFS_Control%levs, GFS_Control%blksz, GFS_Control%dtp,                                               &
            GFS_Control%input_nml_file, GFS_Control%fn_nml, GFS_Control%nlunit, GFS_Control%do_sppt, GFS_Control%do_shum,                &
            GFS_Control%do_skeb, GFS_Control%do_sfcperts, GFS_Control%use_zmtnblck, GFS_Control%skeb_npass, GFS_Control%nsfcpert,        &
            GFS_Control%pertz0, GFS_Control%pertzt, GFS_Control%pertshc, GFS_Control%pertlai, GFS_Control%pertalb, GFS_Control%pertvegf, &
            GFS_Control%ak, GFS_Control%bk, nthreads, GFS_Control%master, GFS_Control%communicator)
      end if

      ! Get land surface perturbations here (move to "else" block below if wanting to update each time-step)
      if (GFS_Control%do_sfcperts) then
         ! Copy blocked data into contiguous arrays; no need to copy sfc_wts in (intent out)
         allocate(xlat(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
         allocate(xlon(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
         allocate(sfc_wts(1:Atm_block%nblks,maxval(GFS_Control%blksz),1:GFS_Control%levs))
         do nb=1,Atm_block%nblks
            xlat(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Grid%xlat(:)
            xlon(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Grid%xlon(:)
         end do
         call run_stochastic_physics_sfc(GFS_Control%blksz, xlat=xlat, xlon=xlon, sfc_wts=sfc_wts)
         ! Copy contiguous data back; no need to copy xlat/xlon, these are intent(in) - just deallocate
         do nb=1,Atm_block%nblks
            GFS_Data(nb)%Coupling%sfc_wts(:,:) = sfc_wts(nb,1:GFS_Control%blksz(nb),:)
         end do
         deallocate(xlat)
         deallocate(xlon)
         deallocate(sfc_wts)
      end if

      ! Consistency check for cellular automata
      if(GFS_Control%do_ca)then
        ! DH* The current implementation of cellular_automata assumes that all blocksizes are the
        ! same - abort if this is not the case, otherwise proceed with Atm_block%blksz(1) below
        if (.not. minval(Atm_block%blksz)==maxval(Atm_block%blksz)) then
           call mpp_error(FATAL, 'Logic errror: cellular_automata not compatible with non-uniform blocksizes')
        end if
        ! *DH
      endif

    else initalize_stochastic_physics

      if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb) then
         ! Copy blocked data into contiguous arrays; no need to copy weights in (intent(out))
         allocate(xlat(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
         allocate(xlon(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
         do nb=1,Atm_block%nblks
            xlat(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Grid%xlat(:)
            xlon(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Grid%xlon(:)
         end do
         if (GFS_Control%do_sppt) then
            allocate(sppt_wts(1:Atm_block%nblks,maxval(GFS_Control%blksz),1:GFS_Control%levs))
         end if
         if (GFS_Control%do_shum) then
            allocate(shum_wts(1:Atm_block%nblks,maxval(GFS_Control%blksz),1:GFS_Control%levs))
         end if
         if (GFS_Control%do_skeb) then
            allocate(skebu_wts(1:Atm_block%nblks,maxval(GFS_Control%blksz),1:GFS_Control%levs))
            allocate(skebv_wts(1:Atm_block%nblks,maxval(GFS_Control%blksz),1:GFS_Control%levs))
         end if
         call run_stochastic_physics(GFS_Control%levs, GFS_Control%kdt, GFS_Control%phour, GFS_Control%blksz, xlat=xlat, xlon=xlon, &
                                     sppt_wts=sppt_wts, shum_wts=shum_wts, skebu_wts=skebu_wts, skebv_wts=skebv_wts, nthreads=nthreads)
         ! Copy contiguous data back; no need to copy xlat/xlon, these are intent(in) - just deallocate
         deallocate(xlat)
         deallocate(xlon)
         if (GFS_Control%do_sppt) then
            do nb=1,Atm_block%nblks
                GFS_Data(nb)%Coupling%sppt_wts(:,:) = sppt_wts(nb,1:GFS_Control%blksz(nb),:)
            end do
            deallocate(sppt_wts)
         end if
         if (GFS_Control%do_shum) then
            do nb=1,Atm_block%nblks
                GFS_Data(nb)%Coupling%shum_wts(:,:) = shum_wts(nb,1:GFS_Control%blksz(nb),:)
            end do
            deallocate(shum_wts)
         end if
         if (GFS_Control%do_skeb) then
            do nb=1,Atm_block%nblks
                GFS_Data(nb)%Coupling%skebu_wts(:,:) = skebu_wts(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Coupling%skebv_wts(:,:) = skebv_wts(nb,1:GFS_Control%blksz(nb),:)
            end do
            deallocate(skebu_wts)
            deallocate(skebv_wts)
         end if
      end if

    endif initalize_stochastic_physics

    ! Cellular automata code is identical for initialization (kstep=0) and time integration (kstep>0)
    if(GFS_Control%do_ca)then
       if(GFS_Control%ca_sgs)then
         ! Allocate contiguous arrays; copy in as needed
         allocate(ugrs        (1:Atm_block%nblks,maxval(GFS_Control%blksz),1:GFS_Control%levs))
         allocate(qgrs        (1:Atm_block%nblks,maxval(GFS_Control%blksz),1:GFS_Control%levs))
         allocate(pgr         (1:Atm_block%nblks,maxval(GFS_Control%blksz)                   ))
         allocate(vvl         (1:Atm_block%nblks,maxval(GFS_Control%blksz),1:GFS_Control%levs))
         allocate(prsl        (1:Atm_block%nblks,maxval(GFS_Control%blksz),1:GFS_Control%levs))
         allocate(ca_deep_diag(1:Atm_block%nblks,maxval(GFS_Control%blksz)                   ))
         allocate(ca_turb_diag(1:Atm_block%nblks,maxval(GFS_Control%blksz)                   ))
         allocate(ca_shal_diag(1:Atm_block%nblks,maxval(GFS_Control%blksz)                   ))
         allocate(condition   (1:Atm_block%nblks,maxval(GFS_Control%blksz)                   ))
         allocate(ca_deep_cpl (1:Atm_block%nblks,maxval(GFS_Control%blksz)                   ))
         allocate(ca_turb_cpl (1:Atm_block%nblks,maxval(GFS_Control%blksz)                   ))
         allocate(ca_shal_cpl (1:Atm_block%nblks,maxval(GFS_Control%blksz)                   ))
         do nb=1,Atm_block%nblks
             ugrs       (nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Statein%ugrs(:,:)
             qgrs       (nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Statein%qgrs(:,:,1)
             pgr        (nb,1:GFS_Control%blksz(nb))   = GFS_Data(nb)%Statein%pgr(:)
             vvl        (nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Statein%vvl(:,:)
             prsl       (nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Statein%prsl(:,:)
             condition  (nb,1:GFS_Control%blksz(nb))   = GFS_Data(nb)%Coupling%condition(:)
             ca_deep_cpl(nb,1:GFS_Control%blksz(nb))   = GFS_Data(nb)%Coupling%ca_deep(:)
             ca_turb_cpl(nb,1:GFS_Control%blksz(nb))   = GFS_Data(nb)%Coupling%ca_turb(:)
             ca_shal_cpl(nb,1:GFS_Control%blksz(nb))   = GFS_Data(nb)%Coupling%ca_shal(:)
         enddo
         call cellular_automata_sgs(GFS_Control%kdt,ugrs,qgrs,pgr,vvl,prsl,condition,ca_deep_cpl,ca_turb_cpl,ca_shal_cpl, &
            ca_deep_diag,ca_turb_diag,ca_shal_diag,Atm(mygrid)%domain_for_coupler,Atm_block%nblks,                        &
            Atm_block%isc,Atm_block%iec,Atm_block%jsc,Atm_block%jec,Atm(mygrid)%npx,Atm(mygrid)%npy, GFS_Control%levs,    &
            GFS_Control%nca,GFS_Control%ncells,GFS_Control%nlives,GFS_Control%nfracseed,                                  &
            GFS_Control%nseed,GFS_Control%nthresh,GFS_Control%ca_global,GFS_Control%ca_sgs,GFS_Control%iseed_ca,          &
            GFS_Control%ca_smooth,GFS_Control%nspinup,Atm_block%blksz(1),GFS_Control%master,GFS_Control%communicator)
         ! Copy contiguous data back as needed
         do nb=1,Atm_block%nblks
             GFS_Data(nb)%Intdiag%ca_deep(:)  = ca_deep_diag(nb,1:GFS_Control%blksz(nb))
             GFS_Data(nb)%Intdiag%ca_turb(:)  = ca_turb_diag(nb,1:GFS_Control%blksz(nb))
             GFS_Data(nb)%Intdiag%ca_shal(:)  = ca_shal_diag(nb,1:GFS_Control%blksz(nb))
             GFS_Data(nb)%Coupling%ca_deep(:) = ca_deep_cpl (nb,1:GFS_Control%blksz(nb))
             GFS_Data(nb)%Coupling%ca_turb(:) = ca_turb_cpl (nb,1:GFS_Control%blksz(nb))
             GFS_Data(nb)%Coupling%ca_shal(:) = ca_shal_cpl (nb,1:GFS_Control%blksz(nb))
         enddo
         deallocate(ugrs        )
         deallocate(qgrs        )
         deallocate(pgr         )
         deallocate(vvl         )
         deallocate(prsl        )
         deallocate(condition   )
         deallocate(ca_deep_cpl )
         deallocate(ca_turb_cpl )
         deallocate(ca_shal_cpl )
         deallocate(ca_deep_diag)
         deallocate(ca_turb_diag)
         deallocate(ca_shal_diag)
       endif
       if(GFS_Control%ca_global)then
          ! Allocate contiguous arrays; no need to copy in (intent out)
          allocate(ca1_cpl (1:Atm_block%nblks,maxval(GFS_Control%blksz)))
          allocate(ca2_cpl (1:Atm_block%nblks,maxval(GFS_Control%blksz)))
          allocate(ca3_cpl (1:Atm_block%nblks,maxval(GFS_Control%blksz)))
          allocate(ca1_diag(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
          allocate(ca2_diag(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
          allocate(ca3_diag(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
          call cellular_automata_global(GFS_Control%kdt,ca1_cpl,ca2_cpl,ca3_cpl,ca1_diag,ca2_diag,ca3_diag,Atm(mygrid)%domain_for_coupler, &
            Atm_block%nblks,Atm_block%isc,Atm_block%iec,Atm_block%jsc,Atm_block%jec,Atm(mygrid)%npx,Atm(mygrid)%npy,GFS_Control%levs,      &
            GFS_Control%nca_g,GFS_Control%ncells_g,GFS_Control%nlives_g,GFS_Control%nfracseed,GFS_Control%nseed_g,GFS_Control%nthresh,     &
            GFS_Control%ca_global,GFS_Control%ca_sgs,GFS_Control%iseed_ca,GFS_Control%ca_smooth,GFS_Control%nspinup,Atm_block%blksz(1),    &
            GFS_Control%nsmooth,GFS_Control%ca_amplitude,GFS_Control%master,GFS_Control%communicator)
          ! Copy contiguous data back
          do nb=1,Atm_block%nblks
             GFS_Data(nb)%Coupling%ca1(:) = ca1_cpl (nb,1:GFS_Control%blksz(nb))
             GFS_Data(nb)%Coupling%ca2(:) = ca2_cpl (nb,1:GFS_Control%blksz(nb))
             GFS_Data(nb)%Coupling%ca3(:) = ca3_cpl (nb,1:GFS_Control%blksz(nb))
             GFS_Data(nb)%Intdiag%ca1(:)  = ca1_diag(nb,1:GFS_Control%blksz(nb))
             GFS_Data(nb)%Intdiag%ca2(:)  = ca2_diag(nb,1:GFS_Control%blksz(nb))
             GFS_Data(nb)%Intdiag%ca3(:)  = ca3_diag(nb,1:GFS_Control%blksz(nb))
          enddo
          deallocate(ca1_cpl )
          deallocate(ca2_cpl )
          deallocate(ca3_cpl )
          deallocate(ca1_diag)
          deallocate(ca2_diag)
          deallocate(ca3_diag)
       endif
    endif

  end subroutine stochastic_physics_wrapper

end module stochastic_physics_wrapper_mod
