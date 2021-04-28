module GFS_restart

  use machine,          only: kind_phys
  use GFS_typedefs,     only: GFS_control_type,  GFS_statein_type,  &
                              GFS_stateout_type, GFS_sfcprop_type,  &
                              GFS_coupling_type, GFS_grid_type,     &
                              GFS_tbd_type,      GFS_cldprop_type,  &
                              GFS_radtend_type,  GFS_diag_type,     &
                              GFS_init_type
  use GFS_diagnostics,  only: GFS_externaldiag_type

  type var_subtype
    real(kind=kind_phys), pointer :: var2p(:)   => null()  !< 2D data saved in packed format [dim(ix)]
    real(kind=kind_phys), pointer :: var3p(:,:) => null()  !< 3D data saved in packed format [dim(ix,levs)]
  end type var_subtype

  type GFS_restart_type
    integer           :: num2d                    !< current number of registered 2D restart variables
    integer           :: num3d                    !< current number of registered 3D restart variables
    integer           :: fdiag                    !< index of first diagnostic field in restart file
    integer           :: ldiag                    !< index of last diagnostic field in restart file

    character(len=32), allocatable :: name2d(:)   !< variable name as it will appear in the restart file
    character(len=32), allocatable :: name3d(:)   !< variable name as it will appear in the restart file
    type(var_subtype), allocatable :: data(:,:)   !< holds pointers to data in packed format (allocated to (nblks,max(2d/3dfields))
  end type GFS_restart_type

  public GFS_restart_type, GFS_restart_populate

  CONTAINS
!*******************************************************************************************

!---------------------
! GFS_restart_populate
!---------------------
  subroutine GFS_restart_populate (Restart, Model, Statein, Stateout, Sfcprop, &
                                   Coupling, Grid, Tbd, Cldprop, Radtend, IntDiag, Init_parm, ExtDiag)
!----------------------------------------------------------------------------------------!
!   RESTART_METADATA                                                                         !
!     Restart%num2d          [int*4  ]  number of 2D variables to output             !
!     Restart%num3d          [int*4  ]  number of 3D variables to output             !
!     Restart%name2d         [char=32]  variable name in restart file                !
!     Restart%name3d         [char=32]  variable name in restart file                !
!     Restart%fld2d(:,:,:)   [real*8 ]  pointer to 2D data (im,nblks,MAX_RSTRT)      !
!     Restart%fld3d(:,:,:,:) [real*8 ]  pointer to 3D data (im,levs,nblks,MAX_RSTRT) !
!----------------------------------------------------------------------------------------!
    type(GFS_restart_type),     intent(inout) :: Restart
    type(GFS_control_type),     intent(in)    :: Model
    type(GFS_statein_type),     intent(in)    :: Statein(:)
    type(GFS_stateout_type),    intent(in)    :: Stateout(:)
    type(GFS_sfcprop_type),     intent(in)    :: Sfcprop(:)
    type(GFS_coupling_type),    intent(in)    :: Coupling(:)
    type(GFS_grid_type),        intent(in)    :: Grid(:)
    type(GFS_tbd_type),         intent(in)    :: Tbd(:)
    type(GFS_cldprop_type),     intent(in)    :: Cldprop(:)
    type(GFS_radtend_type),     intent(in)    :: Radtend(:)
    type(GFS_diag_type),        intent(in)    :: IntDiag(:)
    type(GFS_init_type),        intent(in)    :: Init_parm
    type(GFS_externaldiag_type),intent(in)    :: ExtDiag(:)

    !--- local variables
    integer :: idx, ndiag_rst
    integer :: ndiag_idx(20)
    integer :: nblks, num, nb, max_rstrt, offset 
    character(len=2) :: c2 = ''
    
    nblks = size(Init_parm%blksz)
    max_rstrt = size(Restart%name2d)

    !--- check if continuous accumulated total precip and total cnvc precip are
    !requested in output 
    ndiag_rst = 0
    ndiag_idx(1:20) = 0
    do idx=1, size(ExtDiag)
      if( ExtDiag(idx)%id > 0) then
        if( trim(ExtDiag(idx)%name) == 'totprcp_ave') then
          ndiag_rst = ndiag_rst +1
          ndiag_idx(ndiag_rst) = idx
        else if( trim(ExtDiag(idx)%name) == 'cnvprcp_ave') then
          ndiag_rst = ndiag_rst +1
          ndiag_idx(ndiag_rst) = idx
        else if( trim(ExtDiag(idx)%name) == 'totice_ave') then
          ndiag_rst = ndiag_rst +1
          ndiag_idx(ndiag_rst) = idx
        else if( trim(ExtDiag(idx)%name) == 'totsnw_ave') then
          ndiag_rst = ndiag_rst +1
          ndiag_idx(ndiag_rst) = idx
        else if( trim(ExtDiag(idx)%name) == 'totgrp_ave') then
          ndiag_rst = ndiag_rst +1
          ndiag_idx(ndiag_rst) = idx
        endif
      endif
    enddo

    ! Store first and last index of diagnostic fields:
    Restart%fdiag = 3 + Model%ntot2d + Model%nctp + 1
    Restart%ldiag = 3 + Model%ntot2d + Model%nctp + ndiag_rst
    Restart%num2d = 3 + Model%ntot2d + Model%nctp + ndiag_rst

    ! GF
    if (Model%imfdeepcnv == 3) then
      Restart%num2d = Restart%num2d + 1
    endif
    ! NoahMP
    if (Model%lsm == Model%lsm_noahmp) then
      Restart%num2d = Restart%num2d + 10
    endif
    ! RUC 
    if (Model%lsm == Model%lsm_ruc) then
      Restart%num2d = Restart%num2d + 5
    endif
    ! MYNN SFC
    if (Model%do_mynnsfclay) then
      Restart%num2d = Restart%num2d + 1
    endif
    ! Thompson aerosol-aware
    if (Model%imp_physics == Model%imp_physics_thompson .and. Model%ltaerosol) then
      Restart%num2d = Restart%num2d + 2
    endif

    Restart%num3d = Model%ntot3d
    if(Model%lrefres) then
       Restart%num3d = Model%ntot3d+1
    endif
    ! GF
    if (Model%imfdeepcnv == 3) then
      Restart%num3d = Restart%num3d + 3
    endif
    ! MYNN PBL 
    if (Model%do_mynnedmf) then
      Restart%num3d = Restart%num3d + 9
    endif

    allocate (Restart%name2d(Restart%num2d))
    allocate (Restart%name3d(Restart%num3d))
    allocate (Restart%data(nblks,max(Restart%num2d,Restart%num3d)))

    Restart%name2d(:) = ' '
    Restart%name3d(:) = ' '

    !--- Cldprop variables
    Restart%name2d(1) = 'cv'
    Restart%name2d(2) = 'cvt'
    Restart%name2d(3) = 'cvb'
    do nb = 1,nblks
      Restart%data(nb,1)%var2p => Cldprop(nb)%cv(:)
      Restart%data(nb,2)%var2p => Cldprop(nb)%cvt(:)
      Restart%data(nb,3)%var2p => Cldprop(nb)%cvb(:)
    enddo

    !--- phy_f2d variables
    offset = 3
    do num = 1,Model%ntot2d
       !--- set the variable name
      write(c2,'(i2.2)') num
      Restart%name2d(num+offset) = 'phy_f2d_'//c2
      do nb = 1,nblks
        Restart%data(nb,num+offset)%var2p => Tbd(nb)%phy_f2d(:,num)
      enddo
    enddo
    offset = offset + Model%ntot2d

    !--- phy_fctd variables
    if (Model%nctp > 0) then
      do num = 1, Model%nctp
       !--- set the variable name
        write(c2,'(i2.2)') num
        Restart%name2d(num+offset) = 'phy_fctd_'//c2
        do nb = 1,nblks
          Restart%data(nb,num+offset)%var2p => Tbd(nb)%phy_fctd(:,num)
        enddo
      enddo
      offset = offset + Model%nctp
    endif

    !--- Diagnostic variables
    do idx = 1,ndiag_rst
      if( ndiag_idx(idx) > 0 ) then
        Restart%name2d(offset+idx) = trim(ExtDiag(ndiag_idx(idx))%name)
        do nb = 1,nblks
          Restart%data(nb,offset+idx)%var2p => ExtDiag(ndiag_idx(idx))%data(nb)%var2
        enddo
      endif
!      print *,'in restart 2d field, Restart%name2d(',offset+idx,')=',trim(Restart%name2d(offset+idx))
    enddo

    !--- RAP/HRRR-specific variables, 2D
    num = offset + ndiag_rst
    ! GF
    if (Model%imfdeepcnv == Model%imfdeepcnv_gf) then
      num = num + 1
      Restart%name2d(num) = 'gf_2d_conv_act'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop(nb)%conv_act(:)
      enddo
    endif
    ! NoahMP
    if (Model%lsm == Model%lsm_noahmp) then
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_raincprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop(nb)%raincprv(:)
      enddo
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_rainncprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop(nb)%rainncprv(:)
      enddo
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_iceprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop(nb)%iceprv(:)
      enddo
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_snowprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop(nb)%snowprv(:)
      enddo
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_graupelprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop(nb)%graupelprv(:)
      enddo
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_draincprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop(nb)%draincprv(:)
      enddo
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_drainncprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop(nb)%drainncprv(:)
      enddo
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_diceprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop(nb)%diceprv(:)
      enddo
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_dsnowprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop(nb)%dsnowprv(:)
      enddo
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_dgraupelprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop(nb)%dgraupelprv(:)
      enddo
    endif
    ! RUC 
    if (Model%lsm == Model%lsm_ruc) then
      num = num + 1
      Restart%name2d(num) = 'ruc_2d_raincprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop(nb)%raincprv(:)
      enddo
      num = num + 1
      Restart%name2d(num) = 'ruc_2d_rainncprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop(nb)%rainncprv(:)
      enddo
      num = num + 1
      Restart%name2d(num) = 'ruc_2d_iceprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop(nb)%iceprv(:)
      enddo
      num = num + 1
      Restart%name2d(num) = 'ruc_2d_snowprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop(nb)%snowprv(:)
      enddo
      num = num + 1
      Restart%name2d(num) = 'ruc_2d_graupelprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop(nb)%graupelprv(:)
      enddo
    endif
    ! MYNN SFC
    if (Model%do_mynnsfclay) then
        num = num + 1
        Restart%name2d(num) = 'mynn_2d_uustar'
        do nb = 1,nblks
          Restart%data(nb,num)%var2p => Sfcprop(nb)%uustar(:)
        enddo
    endif
    ! Thompson aerosol-aware
    if (Model%imp_physics == Model%imp_physics_thompson .and. Model%ltaerosol) then
      num = num + 1
      Restart%name2d(num) = 'thompson_2d_nwfa2d'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Coupling(nb)%nwfa2d(:)
      enddo
      num = num + 1
      Restart%name2d(num) = 'thompson_2d_nifa2d'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Coupling(nb)%nifa2d(:)
      enddo
    endif

    !--- phy_f3d variables
    do num = 1,Model%ntot3d
       !--- set the variable name
      write(c2,'(i2.2)') num
      Restart%name3d(num) = 'phy_f3d_'//c2
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd(nb)%phy_f3d(:,:,num)
      enddo
    enddo
    if (Model%lrefres) then
      num = Model%ntot3d+1
      restart%name3d(num) = 'ref_f3d'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => IntDiag(nb)%refl_10cm(:,:)
      enddo
    endif

    if (Model%lrefres) then
       num = Model%ntot3d+1
    else
       num = Model%ntot3d
    endif
    !--- RAP/HRRR-specific variables, 3D
    ! GF
    if (Model%imfdeepcnv == Model%imfdeepcnv_gf) then
      num = num + 1
      Restart%name3d(num) = 'gf_3d_prevst'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd(nb)%prevst(:,:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'gf_3d_prevsq'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd(nb)%prevsq(:,:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'gf_3d_qci_conv'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Coupling(nb)%qci_conv(:,:)
      enddo
    endif
    ! MYNN PBL
    if (Model%do_mynnedmf) then
      num = num + 1
      Restart%name3d(num) = 'mynn_3d_cldfra_bl'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd(nb)%cldfra_bl(:,:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'mynn_3d_qc_bl'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd(nb)%qc_bl(:,:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'mynn_3d_qi_bl'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd(nb)%qi_bl(:,:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'mynn_3d_el_pbl'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd(nb)%el_pbl(:,:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'mynn_3d_sh3d'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd(nb)%sh3d(:,:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'mynn_3d_qke'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd(nb)%qke(:,:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'mynn_3d_tsq'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd(nb)%tsq(:,:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'mynn_3d_qsq'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd(nb)%qsq(:,:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'mynn_3d_cov'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd(nb)%cov(:,:)
      enddo
    endif

  end subroutine GFS_restart_populate

end module GFS_restart
