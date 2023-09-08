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
  real(kind=kind_phys), dimension(:,:,:,:), allocatable, save :: spp_wts

  logical, save :: is_initialized = .false.
  integer, save :: lsoil = -999
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: smc
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: stc
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: slc
  !
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: vfrac
  !albedo
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: snoalb
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: alnsf
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: alnwf
  !emissivity
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: semis
  !roughness length for land
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: zorll

  !real(kind=kind_phys), dimension(:,:),   allocatable, save :: stype
  integer, dimension(:,:),   allocatable, save :: stype

  ! For cellular automata
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: sst
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: lmsk
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: lake
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: uwind
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: vwind
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: height
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: dx
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: condition
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: ca_deep_cpl, ca_turb_cpl, ca_shal_cpl
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: ca1_cpl, ca2_cpl, ca3_cpl


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
    use lndp_apply_perts_mod,         only: lndp_apply_perts

    implicit none

    type(GFS_control_type),   intent(inout) :: GFS_Control
    type(GFS_data_type),      intent(inout) :: GFS_Data(:)
    type(block_control_type), intent(inout) :: Atm_block
    integer,                  intent(out)   :: ierr

    integer :: nthreads, nb, levs, maxblk, nblks, n, v
    logical :: param_update_flag

#ifdef _OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif
    ierr = 0

    levs   = GFS_Control%levs
    maxblk = maxval(GFS_Control%blksz)
    nblks  = Atm_block%nblks

    ! Initialize

    initalize_stochastic_physics: if (.not. is_initialized) then

      if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb .OR. (GFS_Control%lndp_type > 0) .OR. GFS_Control%do_spp) then
         allocate(xlat(1:nblks,maxblk))
         allocate(xlon(1:nblks,maxblk))
         do nb=1,nblks
            xlat(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Grid%xlat(:)
            xlon(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Grid%xlon(:)
         end do
        ! Initialize stochastic physics
        call init_stochastic_physics(levs, GFS_Control%blksz, GFS_Control%dtp, GFS_Control%sppt_amp,                                  &
            GFS_Control%input_nml_file, GFS_Control%fn_nml, GFS_Control%nlunit, xlon, xlat, GFS_Control%do_sppt, GFS_Control%do_shum, &
            GFS_Control%do_skeb, GFS_Control%lndp_type, GFS_Control%n_var_lndp, GFS_Control%use_zmtnblck, GFS_Control%skeb_npass,     &
            GFS_Control%lndp_var_list, GFS_Control%lndp_prt_list,    &
            GFS_Control%n_var_spp, GFS_Control%spp_var_list, GFS_Control%spp_prt_list, GFS_Control%spp_stddev_cutoff, GFS_Control%do_spp,                            &
            GFS_Control%ak, GFS_Control%bk, nthreads, GFS_Control%master, GFS_Control%communicator, ierr)
            if (ierr/=0)  then
                    write(6,*) 'call to init_stochastic_physics failed'
                    return
            endif
      end if
      if (GFS_Control%do_sppt) then
         allocate(sppt_wts(1:nblks,maxblk,1:levs))
      end if
      if (GFS_Control%do_shum) then
         allocate(shum_wts(1:nblks,maxblk,1:levs))
      end if
      if (GFS_Control%do_skeb) then
         allocate(skebu_wts(1:nblks,maxblk,1:levs))
         allocate(skebv_wts(1:nblks,maxblk,1:levs))
      end if
      if ( GFS_Control%do_spp ) then
         allocate(spp_wts(1:nblks,maxblk,1:levs,1:GFS_Control%n_var_spp))
         do n=1,GFS_Control%n_var_spp
           select case (trim(GFS_Control%spp_var_list(n)))
           case('pbl')
             GFS_Control%spp_pbl = 1
           case('sfc')
             GFS_Control%spp_sfc = 1
           case('mp')
             GFS_Control%spp_mp = 7
           case('rad')
             GFS_Control%spp_rad = 1
           case('gwd')
             GFS_Control%spp_gwd = 1
           end select
         end do
      end if
      if ( GFS_Control%lndp_type == 2 ) then
          allocate(sfc_wts(1:nblks,maxblk,1:GFS_Control%n_var_lndp))
          if ( (GFS_Control%lsm == GFS_Control%lsm_noah) .or. (GFS_Control%lsm == GFS_Control%lsm_noahmp)) then
            lsoil = GFS_Control%lsoil
          elseif (GFS_Control%lsm == GFS_Control%lsm_ruc) then
            lsoil = GFS_Control%lsoil_lsm
          endif
          allocate(smc   (1:nblks, maxblk, lsoil))
          do v = 1,GFS_Control%n_var_lndp
            select case (trim(GFS_Control%lndp_var_list(v)))
            case('smc')
              allocate(slc   (1:nblks, maxblk, lsoil))
              allocate(stype (1:nblks, maxblk))
            case('stc')
              allocate(stc   (1:nblks, maxblk, lsoil))
            case('vgf') 
              allocate(vfrac (1:nblks, maxblk))
            case('alb') 
              allocate(alnsf (1:nblks, maxblk))
              allocate(alnwf (1:nblks, maxblk))
            case('sal') 
              allocate(snoalb(1:nblks, maxblk))
            case('emi') 
              allocate(semis (1:nblks, maxblk))
            case('zol') 
              allocate(zorll (1:nblks, maxblk))
            endselect 
          enddo 
      endif


      if ( GFS_Control%lndp_type == 1 ) then ! this scheme sets perts once
         allocate(sfc_wts(1:nblks, maxblk, GFS_Control%n_var_lndp))
         call run_stochastic_physics(levs, GFS_Control%kdt, GFS_Control%fhour, GFS_Control%blksz,       &
                                     sppt_wts=sppt_wts, shum_wts=shum_wts, skebu_wts=skebu_wts,         &
                                     skebv_wts=skebv_wts, sfc_wts=sfc_wts,                              &
                                     spp_wts=spp_wts, nthreads=nthreads)
         ! Copy contiguous data back
         do nb=1,nblks
            GFS_Data(nb)%Coupling%sfc_wts(:,:) = sfc_wts(nb,1:GFS_Control%blksz(nb),:)
         end do
         deallocate(sfc_wts)
      end if
      ! Consistency check for cellular automata
      if(GFS_Control%do_ca)then
        if(GFS_Control%ca_sgs)then
           allocate(sst         (1:nblks, maxblk))
           allocate(lmsk        (1:nblks, maxblk))
           allocate(lake        (1:nblks, maxblk))
           allocate(uwind       (1:nblks, maxblk, 1:levs))
           allocate(vwind       (1:nblks, maxblk, 1:levs))
           allocate(height      (1:nblks, maxblk, 1:levs))
           allocate(condition   (1:nblks, maxblk))
           allocate(dx          (1:nblks, maxblk))
           allocate(ca_deep_cpl (1:nblks, maxblk))
           allocate(ca_turb_cpl (1:nblks, maxblk))
           allocate(ca_shal_cpl (1:nblks, maxblk))
        endif
        if(GFS_Control%ca_global)then
          ! Allocate contiguous arrays; no need to copy in (intent out)
          allocate(ca1_cpl (1:nblks, maxblk))
          allocate(ca2_cpl (1:nblks, maxblk))
          allocate(ca3_cpl (1:nblks, maxblk))
        endif
      endif

      is_initialized = .true.

    else initalize_stochastic_physics
      if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb .OR. (GFS_Control%lndp_type == 2) .OR. GFS_Control%do_spp) then
         call run_stochastic_physics(levs, GFS_Control%kdt, GFS_Control%fhour, GFS_Control%blksz, &
                                 sppt_wts=sppt_wts, shum_wts=shum_wts, skebu_wts=skebu_wts, skebv_wts=skebv_wts, sfc_wts=sfc_wts, &
                                 spp_wts=spp_wts, nthreads=nthreads)
         ! Copy contiguous data back
         if (GFS_Control%do_sppt) then
            do nb=1,nblks
                GFS_Data(nb)%Coupling%sppt_wts(:,:) = sppt_wts(nb,1:GFS_Control%blksz(nb),:)
            end do
         end if
         if (GFS_Control%do_shum) then
            do nb=1,nblks
                GFS_Data(nb)%Coupling%shum_wts(:,:) = shum_wts(nb,1:GFS_Control%blksz(nb),:)
            end do
         end if
         if (GFS_Control%do_skeb) then
            do nb=1,nblks
                GFS_Data(nb)%Coupling%skebu_wts(:,:) = skebu_wts(nb,1:GFS_Control%blksz(nb),:)
                GFS_Data(nb)%Coupling%skebv_wts(:,:) = skebv_wts(nb,1:GFS_Control%blksz(nb),:)
            end do
         end if
         if (GFS_Control%do_spp) then
            do n=1,GFS_Control%n_var_spp
               select case (trim(GFS_Control%spp_var_list(n)))
               case('pbl')
                 do nb=1,Atm_block%nblks
                     GFS_Data(nb)%Coupling%spp_wts_pbl(:,:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               case('sfc')
                 do nb=1,Atm_block%nblks
                     GFS_Data(nb)%Coupling%spp_wts_sfc(:,:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               case('mp')
                 do nb=1,Atm_block%nblks
                     GFS_Data(nb)%Coupling%spp_wts_mp(:,:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               case('gwd')
                 do nb=1,Atm_block%nblks
                     GFS_Data(nb)%Coupling%spp_wts_gwd(:,:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               case('rad')
                 do nb=1,Atm_block%nblks
                     GFS_Data(nb)%Coupling%spp_wts_rad(:,:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               end select
            end do
         end if

         if (GFS_Control%lndp_type == 2) then ! save wts, and apply lndp scheme
             do nb=1,nblks
                GFS_Data(nb)%Coupling%sfc_wts(:,:) = sfc_wts(nb,1:GFS_Control%blksz(nb),:)
             end do
 
             do nb=1,nblks
                do v = 1,GFS_Control%n_var_lndp
                  ! used to identify locations with land model (=soil) 
                  if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                     smc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Data(nb)%Sfcprop%smois(1:GFS_Control%blksz(nb),1:lsoil)
                  else  ! noah or noah-MP
                     smc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Data(nb)%Sfcprop%smc(1:GFS_Control%blksz(nb),1:lsoil)
                  endif

                  select case (trim(GFS_Control%lndp_var_list(v)))
                  case('smc')
                      ! stype used to fetch soil params
                      stype(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%stype(1:GFS_Control%blksz(nb))
                      if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                         slc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Data(nb)%Sfcprop%sh2o(1:GFS_Control%blksz(nb),1:lsoil)
                      else  ! noah or noah-MP
                         slc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Data(nb)%Sfcprop%slc(1:GFS_Control%blksz(nb),1:lsoil)
                      endif
                  case('stc')
                      if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                         stc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Data(nb)%Sfcprop%tslb(1:GFS_Control%blksz(nb),1:lsoil)
                      else ! noah or noah-MP 
                         stc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Data(nb)%Sfcprop%stc(1:GFS_Control%blksz(nb),1:lsoil)
                      endif 
                  case('vgf')
                      if ( (GFS_Control%lsm == GFS_Control%lsm_noahmp) ) then 
                         ! assumes iopt_dveg = 4 (will be checked later)
                         vfrac(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%shdmax(1:GFS_Control%blksz(nb))
                      else ! ruc or noah-MP
                         vfrac(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%vfrac(1:GFS_Control%blksz(nb))
                      endif
                  case('alb')
                      alnsf(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%alnsf(1:GFS_Control%blksz(nb))
                      alnwf(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%alnwf(1:GFS_Control%blksz(nb))
                  case('sal')
                      snoalb(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Sfcprop%snoalb(1:GFS_Control%blksz(nb))
                  case('emi')
                      semis(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Radtend%semis(1:GFS_Control%blksz(nb))
                  case('zol')
                      zorll(nb,1:GFS_Control%blksz(nb))  = GFS_Data(nb)%Sfcprop%zorll(1:GFS_Control%blksz(nb))
                  endselect
              enddo
             enddo


             param_update_flag = .false.
             ! noah and noah-MP treated differently, as global cycle doesn't overwrite shdmax for Noah-MP 
             ! determine whether land paramaters have been over-written to
             ! trigger applying perturbations (logic copied from GFS_driver),
             if ( (GFS_Control%lsm == GFS_Control%lsm_noah) .and. GFS_Control%nscyc >  0) then
                 if (mod(GFS_Control%kdt,GFS_Control%nscyc) == 1 ) then
                   param_update_flag = .true.
                 endif
             endif
             if ( ( GFS_Control%nscyc ==  0 .or. GFS_Control%lsm == GFS_Control%lsm_noahmp) .and. GFS_Control%first_time_step ) then
             ! call once at start of the forecast.
                    param_update_flag = .true.
             endif
              
             call lndp_apply_perts(GFS_Control%blksz, GFS_Control%lsm, GFS_Control%lsm_noah, GFS_Control%lsm_ruc,             &
                               GFS_Control%lsm_noahmp, GFS_Control%iopt_dveg, lsoil, GFS_Control%dtp, GFS_Control%kdt,        &
                               GFS_Control%n_var_lndp, GFS_Control%lndp_var_list, GFS_Control%lndp_prt_list,                  &
                               sfc_wts, xlon, xlat, stype, GFS_Control%pores, GFS_Control%resid,param_update_flag,            &
                               smc, slc, stc, vfrac, alnsf, alnwf, snoalb, semis, zorll, ierr)

             if (ierr/=0)  then
                    write(6,*) 'call to GFS_apply_lndp failed'
                    return
             endif

             do nb=1,nblks
                do v = 1,GFS_Control%n_var_lndp

                  select case (trim(GFS_Control%lndp_var_list(v)))
                  case('smc')
                      if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                           GFS_Data(nb)%Sfcprop%smois(1:GFS_Control%blksz(nb),1:lsoil) = smc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                           GFS_Data(nb)%Sfcprop%sh2o(1:GFS_Control%blksz(nb),1:lsoil)  = slc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                      else  ! noah or noah-MP
                           GFS_Data(nb)%Sfcprop%smc(1:GFS_Control%blksz(nb),1:lsoil) = smc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                           GFS_Data(nb)%Sfcprop%slc(1:GFS_Control%blksz(nb),1:lsoil) = slc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                      endif
                  case('stc')
                      if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                           GFS_Data(nb)%Sfcprop%tslb(1:GFS_Control%blksz(nb),1:lsoil)  = stc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                      else ! noah or noah-MP 
                           GFS_Data(nb)%Sfcprop%stc(1:GFS_Control%blksz(nb),1:lsoil) = stc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                      endif 
                  case('vgf')
                      if ( (GFS_Control%lsm == GFS_Control%lsm_noahmp) ) then 
                        GFS_Data(nb)%Sfcprop%shdmax(1:GFS_Control%blksz(nb))  = vfrac(nb,1:GFS_Control%blksz(nb))
                      else 
                        GFS_Data(nb)%Sfcprop%vfrac(1:GFS_Control%blksz(nb))  = vfrac(nb,1:GFS_Control%blksz(nb))
                      endif
                  case('alb')
                       GFS_Data(nb)%Sfcprop%alnsf(1:GFS_Control%blksz(nb))  = alnsf(nb,1:GFS_Control%blksz(nb))
                       GFS_Data(nb)%Sfcprop%alnwf(1:GFS_Control%blksz(nb))  = alnwf(nb,1:GFS_Control%blksz(nb))
                  case('sal')
                        GFS_Data(nb)%Sfcprop%snoalb(1:GFS_Control%blksz(nb)) = snoalb(nb,1:GFS_Control%blksz(nb))
                  case('emi')
                        GFS_Data(nb)%Radtend%semis(1:GFS_Control%blksz(nb))  = semis(nb,1:GFS_Control%blksz(nb))
                  case('zol')
                        GFS_Data(nb)%Sfcprop%zorll(1:GFS_Control%blksz(nb))  = zorll(nb,1:GFS_Control%blksz(nb))
                  end select   
                enddo 
            enddo
         endif ! lndp block
      endif ! if do* block

      if (GFS_Control%do_ca) then

       if(GFS_Control%ca_sgs)then
         do nb=1,nblks
             sst        (nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Sfcprop%tsfco(:)
             lmsk       (nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Sfcprop%slmsk(:)
             lake       (nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Sfcprop%lakefrac(:)
             uwind      (nb,1:GFS_Control%blksz(nb),:) = GFS_Data(nb)%Statein%ugrs(:,:)
             vwind      (nb,1:GFS_Control%blksz(nb),:) =  GFS_Data(nb)%Statein%vgrs(:,:)
             height     (nb,1:GFS_Control%blksz(nb),:) =  GFS_Data(nb)%Statein%phil(:,:)
             dx         (nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Grid%dx(:)
             condition  (nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Coupling%condition(:)
             ca_deep_cpl(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Coupling%ca_deep(:)
             ca_turb_cpl(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Coupling%ca_turb(:)
             ca_shal_cpl(nb,1:GFS_Control%blksz(nb)) = GFS_Data(nb)%Coupling%ca_shal(:)
         enddo
         call cellular_automata_sgs(GFS_Control%kdt,GFS_control%dtp,GFS_control%restart,GFS_Control%first_time_step,              &
            sst,lmsk,lake,uwind,vwind,height,dx,condition,ca_deep_cpl,ca_turb_cpl,ca_shal_cpl, Atm(mygrid)%domain_for_coupler,nblks,      &
            Atm_block%isc,Atm_block%iec,Atm_block%jsc,Atm_block%jec,Atm(mygrid)%npx,Atm(mygrid)%npy, levs,                        &
            GFS_Control%nthresh,GFS_Control%tile_num,GFS_Control%nca,GFS_Control%ncells,GFS_Control%nlives,                       &
            GFS_Control%nfracseed, GFS_Control%nseed,GFS_Control%iseed_ca,GFS_Control%ca_advect,                                  &
            GFS_Control%nspinup,GFS_Control%ca_trigger,Atm_block%blksz(1),GFS_Control%master,GFS_Control%communicator)
         ! Copy contiguous data back as needed
         do nb=1,nblks
             GFS_Data(nb)%Coupling%ca_deep(:) = ca_deep_cpl (nb,1:GFS_Control%blksz(nb))
             GFS_Data(nb)%Coupling%ca_turb(:) = ca_turb_cpl (nb,1:GFS_Control%blksz(nb))
             GFS_Data(nb)%Coupling%ca_shal(:) = ca_shal_cpl (nb,1:GFS_Control%blksz(nb))
         enddo
       endif
       if(GFS_Control%ca_global)then
          call cellular_automata_global(GFS_Control%kdt,GFS_control%restart,GFS_Control%first_time_step,ca1_cpl,ca2_cpl,ca3_cpl,                &
            Atm(mygrid)%domain_for_coupler, nblks,Atm_block%isc,Atm_block%iec,Atm_block%jsc,Atm_block%jec,Atm(mygrid)%npx,Atm(mygrid)%npy,levs, &
            GFS_Control%nca_g,GFS_Control%ncells_g,GFS_Control%nlives_g,GFS_Control%nfracseed,GFS_Control%nseed_g,                              &
            GFS_Control%iseed_ca,GFS_control%tile_num,GFS_Control%ca_smooth,GFS_Control%nspinup,Atm_block%blksz(1),                             &
            GFS_Control%nsmooth,GFS_Control%ca_amplitude,GFS_Control%master,GFS_Control%communicator)
          ! Copy contiguous data back
          do nb=1,nblks
             GFS_Data(nb)%Coupling%ca1(:) = ca1_cpl(nb,1:GFS_Control%blksz(nb))
             GFS_Data(nb)%Coupling%ca2(:) = ca2_cpl(nb,1:GFS_Control%blksz(nb))
             GFS_Data(nb)%Coupling%ca3(:) = ca3_cpl(nb,1:GFS_Control%blksz(nb))
          enddo
       endif

      endif !do_ca

    endif initalize_stochastic_physics

  end subroutine stochastic_physics_wrapper


  subroutine stochastic_physics_wrapper_end (GFS_Control)

  use GFS_typedefs,       only: GFS_control_type, GFS_data_type
  use stochastic_physics, only: finalize_stochastic_physics

  implicit none

  type(GFS_control_type),   intent(inout) :: GFS_Control

  if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb .OR. (GFS_Control%lndp_type > 0) .OR. GFS_Control%do_spp) then
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
      if (GFS_Control%do_spp) then
         if (allocated(spp_wts)) deallocate(spp_wts)
      end if
      if ( GFS_Control%lndp_type == 2 ) then
         lsoil = -999
         if (allocated(sfc_wts)) deallocate(sfc_wts)
      end if
      if (GFS_Control%lndp_type == 2) then
          if (allocated(smc))    deallocate(smc)
          if (allocated(slc))    deallocate(slc)
          if (allocated(stc))    deallocate(stc)
          if (allocated(stype))  deallocate(stype)
          if (allocated(vfrac))  deallocate(vfrac)
          if (allocated(snoalb)) deallocate(snoalb)
          if (allocated(alnsf))  deallocate(alnsf)
          if (allocated(alnwf))  deallocate(alnwf)
          if (allocated(semis))  deallocate(semis)
          if (allocated(zorll))  deallocate(zorll)
      endif
      call finalize_stochastic_physics()
   endif
   if(GFS_Control%do_ca)then
        if(GFS_Control%ca_sgs)then
           deallocate(sst         )
           deallocate(lmsk        )
           deallocate(lake        )
           deallocate(uwind       )
           deallocate(vwind       )
           deallocate(height      )
           deallocate(dx          )
           deallocate(condition   )
           deallocate(ca_deep_cpl )
           deallocate(ca_turb_cpl )
           deallocate(ca_shal_cpl )
        endif
        if(GFS_Control%ca_global)then
            deallocate(ca1_cpl )
            deallocate(ca2_cpl )
            deallocate(ca3_cpl )
        endif
   endif
  end subroutine stochastic_physics_wrapper_end

end module stochastic_physics_wrapper_mod
