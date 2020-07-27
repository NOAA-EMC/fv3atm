module GFS_land_perts

    use machine,                  only: kind_phys

    implicit none

    private

    public :: GFS_apply_lndp

    contains

!====================================================================
! GFS_apply_lndp
!====================================================================
! Driver for applying perturbations to sprecified land states or parameters
! Draper, July 2020. 
! Note on location: requires access to namelist_soilveg

    subroutine GFS_apply_lndp(Model, Coupling, Grid,param_update_flag, Sfcprop, ierr) 

        use GFS_typedefs,       only: GFS_control_type, GFS_grid_type, & 
                                      GFS_coupling_type, GFS_sfcprop_type

        use namelist_soilveg ! needed for maxsmc

        implicit none

        ! intent(in) 
        type(GFS_control_type),         intent(in) :: Model
        type(GFS_coupling_type),        intent(in) :: Coupling(:)
        type(GFS_grid_type),            intent(in) :: Grid(:)
        logical,                        intent(in) ::  param_update_flag    
                                        ! true =  parameters have been updated, apply perts
        ! intent(inout) 
        type(GFS_sfcprop_type),         intent(inout) :: Sfcprop(:)

        ! intent(out) 
        integer,                        intent(out) :: ierr

        ! local
        integer         :: nblks, print_i, print_nb, i, nb
        integer         ::  this_im, v, soiltyp, k
        logical         :: print_flag 

        real(kind=kind_phys) :: p, min_bound, max_bound, tmp_sic,  pert

        ! decrease in applied pert with depth
        real(kind=kind_phys), dimension(4), parameter  :: smc_vertscale = (/1.0,0.5,0.25,0.125/)
        real(kind=kind_phys), dimension(4), parameter  :: stc_vertscale = (/1.0,0.5,0.25,0.125/)

        ! model-dependent values, hard-wired in noah code.
        real(kind=kind_phys), dimension(4), parameter  :: zs_noah = (/0.1, 0.3, 0.6, 1.0/)
        real(kind=kind_phys), parameter                :: minsmc = 0.02

        ierr = 0 

        nblks = size(Model%blksz)

        call  set_printing_nb_i(nblks,Grid,print_i,print_nb)

        do nb =1,nblks
           this_im = size(Sfcprop(nb)%smc(:,1))
           do i = 1, this_im

             if ( nint(Sfcprop(nb)%slmsk(i)) .NE. 1) cycle ! skip if not land
             soiltyp  = int( Sfcprop(nb)%stype(i)+0.5 )  ! also need for maxsmc

             if ( ((Model%isot == 1) .and. (soiltyp == 16)) &
               .or.( (Model%isot == 0) .and. (soiltyp  == 9 )) ) cycle ! skip if land-ice

             ! set printing
             if ( (i==print_i)  .and. (nb==print_nb) ) then 
                print_flag = .true.
             else 
                print_flag=.false. 
             endif

             do v = 1,Model%n_var_lndp 
                select case (trim(Model%lndp_var_list(v)))
                !=================================================================
                ! State updates - performed every cycle
                !=================================================================
                case('smc') 
                    p=5. 
                    min_bound = minsmc
                    max_bound = maxsmc(soiltyp)

                    do k=1,Model%lsoil
                         !store frozen soil moisture
                         tmp_sic= Sfcprop(nb)%smc(i,k)  - Sfcprop(nb)%slc(i,k)

                         ! perturb total soil moisture 
                         ! factor of sldepth*1000 converts from mm to m3/m3
                         pert = Coupling(nb)%sfc_wts(i,v)*smc_vertscale(k)*Model%lndp_prt_list(v)/(zs_noah(k)*1000.)                    
                         pert = pert*Model%dtf/3600. ! lndp_prt_list input is per hour, convert to per timestep 
                                                     ! (necessary for state vars only)
                         call apply_pert('smc',pert,print_flag, Sfcprop(nb)%smc(i,k),ierr,p,min_bound, max_bound)

                         ! assign all of applied pert to the liquid soil moisture 
                         Sfcprop(nb)%slc(i,k)  =  Sfcprop(nb)%smc(i,k) -  tmp_sic
                    enddo

                case('stc') 

                    do k=1,Model%lsoil
                         pert = Coupling(nb)%sfc_wts(i,v)*stc_vertscale(k)*Model%lndp_prt_list(v)
                         pert = pert*Model%dtf/3600. ! lndp_prt_list input is per hour, convert to per timestep
                                                     ! (necessary for state vars only)
                         call apply_pert('stc',pert,print_flag, Sfcprop(nb)%stc(i,k),ierr)
                    enddo
                !=================================================================
                ! Parameter updates - only if param_update_flag = TRUE
                !=================================================================
                case('vgf')  ! vegetation fraction
                     if (param_update_flag) then 
                         p =5. 
                         min_bound=0.
                         max_bound=1.

                         pert = Coupling(nb)%sfc_wts(i,v)*Model%lndp_prt_list(v)
                         call apply_pert ('vfrac',pert,print_flag, Sfcprop(nb)%vfrac(i), ierr,p,min_bound, max_bound)
                     endif
                case default 
                    print*, &
                     'unrecognised lndp_prt_list option in GFS_apply_lndp, exiting', trim(Model%lndp_var_list(v)) 
                    ierr = 10 
                    return 
                end select 
                !if (param_update_flag) then 
                !endif
             enddo 
           enddo 
        enddo
    end subroutine  GFS_apply_lndp

!====================================================================
! apply_perts
!====================================================================
! Apply perturbations to selected land states or parameters

  subroutine apply_pert(vname,pert,print_flag, state,ierr,p,vmin, vmax)

   ! intent in
    logical, intent(in)                 :: print_flag
    real(kind=kind_phys), intent(in)    :: pert
    character(len=*), intent(in)        :: vname ! name of variable being perturbed

    real(kind=kind_phys), optional, intent(in)    :: p ! flat-top paramater, 0 = no flat-top
                                                       ! flat-top function is used for bounded variables 
                                                       ! to reduce the magnitude of perturbations  near boundaries.
    real(kind=kind_phys), optional, intent(in) :: vmin, vmax ! min,max bounds of variable being perturbed

    ! intent (inout)
    real(kind=kind_phys), intent(inout) :: state
    
    ! intent out 
    integer                             :: ierr

    !local
    real(kind=kind_phys) :: z

       if ( print_flag ) then
              write(*,*) 'LNDP - applying lndp to ',vname, ', initial value', state
       endif

       ! apply perturbation
       if (present(p) ) then 
           if ( .not. (present(vmin) .and. present(vmax) )) then
              ierr=20 
              print*, 'error, flat-top function requires min & max to be specified'
           endif

           z = -1. + 2*(state  - vmin)/(vmax - vmin) ! flat-top function
           state =  state  + pert*(1-abs(z**p))
       else
          state =  state  + pert
       endif

       if (present(vmax)) state =  min( state , vmax )
       if (present(vmin)) state =  max( state , vmin )
       !state = max( min( state , vmax ), vmin )

       if ( print_flag ) then
              write(*,*) 'LNDP - applying lndp to ',vname, ', final value', state
       endif

  end subroutine apply_pert
       

!====================================================================
! set_printing_nb_i 
!====================================================================
! routine to turn on print flag for selected location
! 
    subroutine set_printing_nb_i(nblks,Grid,print_i,print_nb)

        use GFS_typedefs,       only: GFS_grid_type

        implicit none 

        ! intent (in)
        integer, intent(in) :: nblks
        type(GFS_grid_type),      intent(in) :: Grid(nblks) 

        ! intent (out)
        integer, intent(out) :: print_i, print_nb

        ! local
        integer :: nb,i,this_im
        real, parameter :: plon_trunc =  114.9
        real, parameter :: plat_trunc =  -26.6
        real, parameter  :: delta  = 1.

        print_i = -9
        print_nb = -9
        do nb = 1,nblks
         this_im = size(Grid(nb)%xlon(:)) 
         do i = 1,this_im
        if ( (Grid(nb)%xlon(i)*57.29578 > plon_trunc) .and.  (Grid(nb)%xlon(i)*57.29578 < plon_trunc+delta ) .and. &
           (Grid(nb)%xlat(i)*57.29578 >  plat_trunc ) .and.  (Grid(nb)%xlat(i)*57.29578 < plat_trunc+delta ) ) then
                      print_i=i
                      print_nb=nb
                      write(*,*) 'LNDP -print flag is on', Grid(nb)%xlon(i)*57.29578, Grid(nb)%xlat(i)*57.29578, nb, i
                      return  
         endif
         enddo
        enddo

    end subroutine set_printing_nb_i

end module GFS_land_perts


