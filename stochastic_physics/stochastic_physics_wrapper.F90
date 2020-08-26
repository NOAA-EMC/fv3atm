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

  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: smc 
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: stc 
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: slc 
  real(kind=kind_phys), dimension(:,:), allocatable, save :: vfrac
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

#ifdef OPENMP
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
    use namelist_soilveg, only: maxsmc

    implicit none

    type(GFS_control_type),   intent(inout) :: GFS_Control
    type(GFS_data_type),      intent(inout) :: GFS_Data(:)
    type(block_control_type), intent(inout) :: Atm_block
    integer,                  intent(out)   :: ierr

    integer :: nthreads, nb
    logical :: param_update_flag

#ifdef OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif
    ierr = 0 

    ! Initialize
    initalize_stochastic_physics: if (GFS_Control%kdt==0) then

      if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb .OR. (GFS_Control%lndp_type .GT. 0) ) then
        ! Initialize stochastic physics
        call init_stochastic_physics(GFS_Control%levs, GFS_Control%blksz, GFS_Control%dtp,                                               &
            GFS_Control%input_nml_file, GFS_Control%fn_nml, GFS_Control%nlunit, GFS_Control%do_sppt, GFS_Control%do_shum,                &
            GFS_Control%do_skeb, GFS_Control%lndp_type, GFS_Control%n_var_lndp, GFS_Control%use_zmtnblck, GFS_Control%skeb_npass, & 
            GFS_Control%lndp_var_list, GFS_Control%lndp_prt_list,    &
            GFS_Control%ak, GFS_Control%bk, nthreads, GFS_Control%master, GFS_Control%communicator, ierr)
            if (ierr/=0)  then 
                    write(6,*) 'call to init_stochastic_physics failed'
                    return
            endif
      end if

      if ( GFS_Control%lndp_type .EQ. 1 ) then ! this scheme sets perts once
         ! Copy blocked data into contiguous arrays; no need to copy sfc_wts in (intent out)
         allocate(xlat(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
         allocate(xlon(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
         allocate(sfc_wts(1:Atm_block%nblks,maxval(GFS_Control%blksz),GFS_Control%n_var_lndp))
         do nb=1,Atm_block%nblks
            xlat(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Grid%xlat(:)
            xlon(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Grid%xlon(:)
         end do
         call run_stochastic_physics(GFS_Control%levs, GFS_Control%kdt, GFS_Control%phour, GFS_Control%blksz, xlat=xlat, xlon=xlon, &
                                 sppt_wts=sppt_wts, shum_wts=shum_wts, skebu_wts=skebu_wts, skebv_wts=skebv_wts, sfc_wts=sfc_wts, &
                                 nthreads=nthreads)
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

      if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb .OR. (GFS_Control%lndp_type .EQ. 2) ) then
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
         if ( GFS_Control%lndp_type .EQ. 2 ) then ! this scheme updates through forecast
            allocate(sfc_wts(1:Atm_block%nblks,maxval(GFS_Control%blksz),1:GFS_Control%n_var_lndp))
         end if 

         call run_stochastic_physics(GFS_Control%levs, GFS_Control%kdt, GFS_Control%phour, GFS_Control%blksz, xlat=xlat, xlon=xlon, &
                                 sppt_wts=sppt_wts, shum_wts=shum_wts, skebu_wts=skebu_wts, skebv_wts=skebv_wts, sfc_wts=sfc_wts, &
                                 nthreads=nthreads)
         ! Copy contiguous data back; no need to copy xlat/xlon, these are intent(in) - just deallocate
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
         if (GFS_Control%lndp_type .EQ. 2) then ! save wts, and apply lndp scheme 
             do nb=1,Atm_block%nblks
                GFS_Data(nb)%Coupling%sfc_wts(:,:) = sfc_wts(nb,1:GFS_Control%blksz(nb),:)
             end do

             allocate(smc(1:Atm_block%nblks,maxval(GFS_Control%blksz),GFS_Control%lsoil))
             allocate(slc(1:Atm_block%nblks,maxval(GFS_Control%blksz),GFS_Control%lsoil))
             allocate(stc(1:Atm_block%nblks,maxval(GFS_Control%blksz),GFS_Control%lsoil))
             allocate(stype(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
             allocate(vfrac(1:Atm_block%nblks,maxval(GFS_Control%blksz)))
             do nb=1,Atm_block%nblks
                stype(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%stype(:)
                smc(nb,1:GFS_Control%blksz(nb),:)  = GFS_Data(nb)%Sfcprop%smc(:,:) 
                slc(nb,1:GFS_Control%blksz(nb),:)  = GFS_Data(nb)%Sfcprop%slc(:,:) 
                stc(nb,1:GFS_Control%blksz(nb),:)  = GFS_Data(nb)%Sfcprop%stc(:,:) 
                vfrac(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%vfrac(:) 
             end do

             ! determine whether land paramaters have been over-written 
             if (mod(GFS_Control%kdt,GFS_Control%nscyc) == 1)  then ! logic copied from GFS_driver
                    param_update_flag = .true. 
             else
                    param_update_flag = .false.
             endif 
             call lndp_apply_perts( GFS_Control%blksz, GFS_Control%lsm,  GFS_Control%lsoil, GFS_Control%dtf, & 
                             GFS_Control%n_var_lndp, GFS_Control%lndp_var_list, GFS_Control%lndp_prt_list, & 
                             sfc_wts, xlon, xlat, stype, maxsmc,param_update_flag, smc, slc,stc, vfrac, ierr) 
             if (ierr/=0)  then 
                    write(6,*) 'call to GFS_apply_lndp failed'
                    return
             endif
             deallocate(stype) 
             deallocate(sfc_wts) 
             do nb=1,Atm_block%nblks
                 GFS_Data(nb)%Sfcprop%smc(:,:) =  smc(nb,1:GFS_Control%blksz(nb),:)
                 GFS_Data(nb)%Sfcprop%slc(:,:) =  slc(nb,1:GFS_Control%blksz(nb),:)
                 GFS_Data(nb)%Sfcprop%stc(:,:) =  stc(nb,1:GFS_Control%blksz(nb),:)
                 GFS_Data(nb)%Sfcprop%vfrac(:) =  vfrac(nb,1:GFS_Control%blksz(nb))
             enddo
             deallocate(smc) 
             deallocate(slc) 
             deallocate(stc) 
             deallocate(vfrac) 
         endif ! lndp block
         deallocate(xlat)
         deallocate(xlon)
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
