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

  integer, save :: lsoil = -999
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: smc
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: stc
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: slc
  !
  real(kind=kind_phys), dimension(:,:), allocatable, save :: vfrac
  !albedo
  real(kind=kind_phys), dimension(:,:), allocatable, save :: snoalb
  real(kind=kind_phys), dimension(:,:), allocatable, save :: alvsf
  real(kind=kind_phys), dimension(:,:), allocatable, save :: alnsf
  real(kind=kind_phys), dimension(:,:), allocatable, save :: alvwf
  real(kind=kind_phys), dimension(:,:), allocatable, save :: alnwf
  real(kind=kind_phys), dimension(:,:), allocatable, save :: facsf
  real(kind=kind_phys), dimension(:,:), allocatable, save :: facwf
  !emissivity
  real(kind=kind_phys), dimension(:,:), allocatable, save :: semis
  !roughness length for land
  real(kind=kind_phys), dimension(:,:), allocatable, save :: zorll

  real(kind=kind_phys), dimension(:,:), allocatable, save :: stype

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
  subroutine stochastic_physics_wrapper (GFS_Control, GFS_Data, Atm_block, ierr)

#ifdef _OPENMP
    use omp_lib
#endif

    use GFS_typedefs,       only: GFS_control_type, GFS_data_type
    use mpp_mod,            only: FATAL, mpp_error
    use block_control_mod,  only: block_control_type
    use atmosphere_mod,     only: Atm, mygrid

    use stochastic_physics,           only: init_stochastic_physics, run_stochastic_physics
    use cellular_automata_global_mod, only: cellular_automata_global
    use cellular_automata_sgs_mod,    only: cellular_automata_sgs
    use lndp_apply_perts_mod, only: lndp_apply_perts

    implicit none

    type(GFS_control_type),   intent(inout) :: GFS_Control
    type(GFS_data_type),      intent(inout) :: GFS_Data(:)
    type(block_control_type), intent(inout) :: Atm_block
    integer,                  intent(out)   :: ierr

    integer :: nthreads, nb
    logical :: param_update_flag

#ifdef _OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif
    ierr = 0

    ! Initialize
    initalize_stochastic_physics: if (GFS_Control%kdt==0) then

      if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb .OR. (GFS_Control%lndp_type .GT. 0) ) then
         allocate(xlat(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
         allocate(xlon(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
         do nb=1,Atm_block%nblks
            xlat(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Grid%xlat(:)
            xlon(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Grid%xlon(:)
         end do
        ! Initialize stochastic physics
        call init_stochastic_physics(GFS_Control%levs, GFS_Control%blksz, GFS_Control%dtp, GFS_Control%sppt_amp,                         &
            GFS_Control%input_nml_file, GFS_Control%fn_nml, GFS_Control%nlunit, xlon, xlat, GFS_Control%do_sppt, GFS_Control%do_shum,                &
            GFS_Control%do_skeb, GFS_Control%lndp_type, GFS_Control%n_var_lndp, GFS_Control%use_zmtnblck, GFS_Control%skeb_npass, &
            GFS_Control%lndp_var_list, GFS_Control%lndp_prt_list,    &
            GFS_Control%ak, GFS_Control%bk, nthreads, GFS_Control%master, GFS_Control%communicator, ierr)
            if (ierr/=0)  then
                    write(6,*) 'call to init_stochastic_physics failed'
                    return
            endif
      end if
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
      if ( GFS_Control%lndp_type .EQ. 2 ) then ! this scheme updates through forecast
         allocate(sfc_wts(1:Atm_block%nblks,maxval(GFS_Control%blksz),1:GFS_Control%n_var_lndp))
      end if
      if (GFS_Control%lndp_type .EQ. 2) then ! save wts, and apply lndp scheme
          if (GFS_Control%lsm == GFS_Control%lsm_noah) then
            lsoil = GFS_Control%lsoil
          elseif (GFS_Control%lsm == GFS_Control%lsm_ruc) then
            lsoil = GFS_Control%lsoil_lsm
          endif
          allocate(smc(1:Atm_block%nblks,maxval(GFS_Control%blksz),lsoil))
          allocate(slc(1:Atm_block%nblks,maxval(GFS_Control%blksz),lsoil))
          allocate(stc(1:Atm_block%nblks,maxval(GFS_Control%blksz),lsoil))
          allocate(stype(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
          allocate(vfrac(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
          allocate(snoalb(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
          allocate(alvsf(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
          allocate(alnsf(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
          allocate(alvwf(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
          allocate(alnwf(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
          allocate(facsf(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
          allocate(facwf(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
          allocate(semis(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
          allocate(zorll(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
      endif


      if ( GFS_Control%lndp_type .EQ. 1 ) then ! this scheme sets perts once
         allocate(sfc_wts(1:Atm_block%nblks,maxval(GFS_Control%blksz),GFS_Control%n_var_lndp))
         call run_stochastic_physics(GFS_Control%levs, GFS_Control%kdt, GFS_Control%fhour, GFS_Control%blksz,                     &
                                 sppt_wts=sppt_wts, shum_wts=shum_wts, skebu_wts=skebu_wts, skebv_wts=skebv_wts, sfc_wts=sfc_wts, &
                                 nthreads=nthreads)
         ! Copy contiguous data back
         do nb=1,Atm_block%nblks
            GFS_Data(nb)%Coupling%sfc_wts(:,:) = sfc_wts(nb,1:GFS_Control%blksz(nb),:)
         end do
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
      if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb .OR. (GFS_Control%lndp_type .EQ. 2) ) then
         call run_stochastic_physics(GFS_Control%levs, GFS_Control%kdt, GFS_Control%fhour, GFS_Control%blksz, &
                                 sppt_wts=sppt_wts, shum_wts=shum_wts, skebu_wts=skebu_wts, skebv_wts=skebv_wts, sfc_wts=sfc_wts, &
                                 nthreads=nthreads)
         ! Copy contiguous data back
         if (GFS_Control%do_sppt) then
            do nb=1,Atm_block%nblks
                GFS_Data(nb)%Coupling%sppt_wts(:,:) = sppt_wts(nb,1:GFS_Control%blksz(nb),:)
            end do
         end if
         if (GFS_Control%do_shum) then
            do nb=1,Atm_block%nblks
                GFS_Data(nb)%Coupling%shum_wts(:,:) = shum_wts(nb,1:GFS_Control%blksz(nb),:)
            end do
         end if
         if (GFS_Control%do_skeb) then
            do nb=1,Atm_block%nblks
                GFS_Data(nb)%Coupling%skebu_wts(:,:) = skebu_wts(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Coupling%skebv_wts(:,:) = skebv_wts(nb,1:GFS_Control%blksz(nb),:)
            end do
         end if
         if (GFS_Control%lndp_type .EQ. 2) then ! save wts, and apply lndp scheme
             do nb=1,Atm_block%nblks
                GFS_Data(nb)%Coupling%sfc_wts(:,:) = sfc_wts(nb,1:GFS_Control%blksz(nb),:)
             end do

             do nb=1,Atm_block%nblks
                stype(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%stype(:)
                vfrac(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%vfrac(:)
                snoalb(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Sfcprop%snoalb(:)
                alvsf(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%alvsf(:)
                alnsf(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%alnsf(:)
                alvwf(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%alvwf(:)
                alnwf(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%alnwf(:)
                facsf(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%facsf(:)
                facwf(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%facwf(:)
                semis(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Radtend%semis(:)
                zorll(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%zorll(:)
             end do

             if (GFS_Control%lsm == GFS_Control%lsm_noah) then
               do nb=1,Atm_block%nblks
                 smc(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Sfcprop%smc(:,:)
                 slc(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Sfcprop%slc(:,:)
                 stc(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Sfcprop%stc(:,:)
               end do
             elseif (GFS_Control%lsm == GFS_Control%lsm_ruc) then
               do nb=1,Atm_block%nblks
                 smc(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Sfcprop%smois(:,:)
                 slc(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Sfcprop%sh2o(:,:)
                 stc(nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Sfcprop%tslb(:,:)
               end do
             endif

             ! determine whether land paramaters have been over-written to
             ! trigger applying perturbations (logic copied from GFS_driver),
             ! or if perturbations should be applied at every time step
             if (mod(GFS_Control%kdt,GFS_Control%nscyc) == 1 ) then
               param_update_flag = .true.
             else
               param_update_flag = .false.
             endif

             call lndp_apply_perts(GFS_Control%blksz, GFS_Control%lsm, GFS_Control%lsm_noah, GFS_Control%lsm_ruc, lsoil,      &
                               GFS_Control%dtf, GFS_Control%kdt, GFS_Control%lndp_each_step,                                  &
                               GFS_Control%n_var_lndp, GFS_Control%lndp_var_list, GFS_Control%lndp_prt_list,                  &
                               sfc_wts, xlon, xlat, stype, GFS_Control%pores, GFS_Control%resid,param_update_flag,            &
                               smc, slc, stc, vfrac, alvsf, alnsf, alvwf, alnwf, facsf, facwf, snoalb, semis, zorll, ierr)

             if (ierr/=0)  then
                    write(6,*) 'call to GFS_apply_lndp failed'
                    return
             endif

             do nb=1,Atm_block%nblks
               GFS_Data(nb)%Sfcprop%vfrac(:)  = vfrac(nb,1:GFS_Control%blksz(nb))
               GFS_Data(nb)%Sfcprop%snoalb(:) = snoalb(nb,1:GFS_Control%blksz(nb))
               GFS_Data(nb)%Sfcprop%alvsf(:)  = alvsf(nb,1:GFS_Control%blksz(nb))
               GFS_Data(nb)%Sfcprop%alnsf(:)  = alnsf(nb,1:GFS_Control%blksz(nb))
               GFS_Data(nb)%Sfcprop%alvwf(:)  = alvwf(nb,1:GFS_Control%blksz(nb))
               GFS_Data(nb)%Sfcprop%alnwf(:)  = alnwf(nb,1:GFS_Control%blksz(nb))
               GFS_Data(nb)%Sfcprop%facsf(:)  = facsf(nb,1:GFS_Control%blksz(nb))
               GFS_Data(nb)%Sfcprop%facwf(:)  = facwf(nb,1:GFS_Control%blksz(nb))
               GFS_Data(nb)%Radtend%semis(:)  = semis(nb,1:GFS_Control%blksz(nb))
               GFS_Data(nb)%Sfcprop%zorll(:)  = zorll(nb,1:GFS_Control%blksz(nb))
             enddo

             if (GFS_Control%lsm == GFS_Control%lsm_noah) then
               do nb=1,Atm_block%nblks
                   GFS_Data(nb)%Sfcprop%smc(:,:) = smc(nb,1:GFS_Control%blksz(nb),:)
                   GFS_Data(nb)%Sfcprop%slc(:,:) = slc(nb,1:GFS_Control%blksz(nb),:)
                   GFS_Data(nb)%Sfcprop%stc(:,:) = stc(nb,1:GFS_Control%blksz(nb),:)
               enddo
             elseif (GFS_Control%lsm == GFS_Control%lsm_ruc) then
               do nb=1,Atm_block%nblks
                   GFS_Data(nb)%Sfcprop%smois(:,:) = smc(nb,1:GFS_Control%blksz(nb),:)
                   GFS_Data(nb)%Sfcprop%sh2o(:,:)  = slc(nb,1:GFS_Control%blksz(nb),:)
                   GFS_Data(nb)%Sfcprop%tslb(:,:)  = stc(nb,1:GFS_Control%blksz(nb),:)
               enddo
             endif

         endif ! lndp block
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


  subroutine stochastic_physics_wrapper_end (GFS_Control)

  use GFS_typedefs,       only: GFS_control_type, GFS_data_type
  use stochastic_physics, only: finalize_stochastic_physics

  implicit none

  type(GFS_control_type),   intent(inout) :: GFS_Control

  if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb .OR. (GFS_Control%lndp_type .GT. 0) ) then
      if (allocated(xlat)) deallocate(xlat)
      if (allocated(xlon)) deallocate(xlon)
      if (GFS_Control%do_sppt) then
         if (allocated(sppt_wts)) deallocate(sppt_wts)
      end if
      if (GFS_Control%do_shum) then
         if (allocated(shum_wts)) deallocate(shum_wts)
      end if
      if (GFS_Control%do_skeb) then
         if (allocated(skebu_wts)) deallocate(skebu_wts)
         if (allocated(skebv_wts)) deallocate(skebv_wts)
      end if
      if ( GFS_Control%lndp_type .EQ. 2 ) then ! this scheme updates through forecast
         lsoil = -999
         if (allocated(sfc_wts)) deallocate(sfc_wts)
      end if
      if (GFS_Control%lndp_type .EQ. 2) then ! save wts, and apply lndp scheme
          if (allocated(smc)) deallocate(smc)
          if (allocated(slc)) deallocate(slc)
          if (allocated(stc)) deallocate(stc)
          if (allocated(stype)) deallocate(stype)
          if (allocated(vfrac)) deallocate(vfrac)
          if (allocated(snoalb)) deallocate(snoalb)
          if (allocated(alvsf)) deallocate(alvsf)
          if (allocated(alnsf)) deallocate(alnsf)
          if (allocated(alvwf)) deallocate(alvwf)
          if (allocated(alnwf)) deallocate(alnwf)
          if (allocated(facsf)) deallocate(facsf)
          if (allocated(facwf)) deallocate(facwf)
          if (allocated(semis)) deallocate(semis)
          if (allocated(zorll)) deallocate(zorll)
      endif
      call finalize_stochastic_physics()
   endif
  end subroutine stochastic_physics_wrapper_end

end module stochastic_physics_wrapper_mod
