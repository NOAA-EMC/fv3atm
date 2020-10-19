      MODULE IDEA_solarno_input


       use IDEA_IO_UNITS, only : iulog
       implicit none
!       public :: solar_wam_get_feuv
!       public :: solar_read_wam_init
!       public :: solar_wamstep_advance
!       public :: solar_waccmx_advance
!       public :: solar_read_myy1947_2016
       public :: solar_readno_snoewx
!       public :: solar_read_namelist
!       public :: dealloc_solar
!
       save 
!
! no-snoe input Marsh et al. (2004)
!
       integer, parameter  ::  no_ny33=33
       integer, parameter  ::  no_nz16=16
       integer, parameter  ::  no_neofs=7       ! total NO-snoe model has 7 modes
       real, allocatable   ::  no_eof(:, : ,:)
       real , allocatable  ::  no_m(:,:)
       real, allocatable   ::  no_zkm(:), no_mlat(:)
      contains
!
      subroutine solar_readno_snoewx(file, mpi_id)
       use netcdf      
       use module_physics_driver, only : is_master
!SK    use idea_mpi_def, ONLY:  info, mpi_comm_all
       implicit none
!SK    include 'mpif.h'
       character(len=*), intent(in) :: file
       integer, intent(in) :: mpi_id   
!locals
       integer :: nz, ny, neofs
       integer :: istat, ierr, astat
       integer :: dim_id, var_id 
       integer :: ncid, vid, iernc
       character(len=256) :: locfn
       integer, dimension(nf90_max_var_dims) :: dimidT  
!----------------------------------------------------------------------
!	... open the netcdf file
!----------------------------------------------------------------------
       if(is_master) then
!SK    if(mpi_id.eq.0) then
         write(iulog,*)file        
         write(iulog,*) ' solar_readno_snoewx: opening file for readno',
     &  trim(file)
        endif
       ierNC=NF90_OPEN(trim(File), nf90_nowrite, ncid)   
       if (iernc /=0) write(iulog,*) ncid, 'ncid ', iernc, ' iernc '
!SK    if (is_master.and.iernc /=0) 
!SK  &    write(iulog,*) ncid, 'ncid ', iernc, ' iernc '
!----------------------------------------------------------------------
!	... read the snoe dimensions
!----------------------------------------------------------------------
      iernc=nf90_inq_varid( ncid, 'EOF', vid )
         ierNC=nf90_inquire_variable(ncid, vid, dimids=dimidT) 
         iernc = nf90_inquire_dimension(ncid, dimidT(3), len=neofs)
         iernc = nf90_inquire_dimension(ncid, dimidT(2), len=nz)
         iernc = nf90_inquire_dimension(ncid, dimidT(1), len=ny)
         if(is_master) then
!SK      if(mpi_id.eq.0) then
          write(iulog,*) neofs, nz, ny, ' ne-nz-ny of NO-EOFs VAY'
          write(iulog,*) no_neofs, no_nz16, no_ny33,
     &                 ' ne-nz-ny of NO-EOFs VAY'
          if (nz.ne.no_nz16.or.ny.ne.no_ny33.or.neofs.ne.no_neofs) then
         write(iulog,*)'snoe_rdeof: failed to read expected neofs=nz=ny'
!SK        call mpi_quit(23901)
          endif
         endif
         
!----------------------------------------------------------------------
!	... allocate snoe variables
!----------------------------------------------------------------------
 
       allocate( no_mlat(ny),  no_zkm(nz),stat=astat )  
       allocate( no_m(ny, nz),stat=astat )
       allocate( no_eof(ny, nz, neofs),stat=astat )
!SK    if( is_master .and. astat /= 0 ) then
       if( astat /= 0 ) then
       write(iulog,*) ' alloc_err in read_no_snoe no_eof ' 
       write(iulog,*) 'snoe_rdeof: failed to allocate eofs; error = ',
     &                astat
       end if 

!----------------------------------------------------------------------
!	... read the snoe variables
!----------------------------------------------------------------------
        iernc=nf90_inq_varid( ncid, 'lat', vid )
        iernc= nf90_get_var( ncid, vid, no_mlat)
        iernc=nf90_inq_varid( ncid, 'z', vid )
        iernc= nf90_get_var( ncid, vid, no_zkm)

      iernc = nf90_inq_varid( ncid, 'NO', var_id )
      ierr  = nf90_get_var( ncid, var_id, no_m )

      iernc = nf90_inq_varid( ncid, 'EOF', var_id )
      iernc = nf90_get_var( ncid, var_id, no_eof )  !(/1,1,1/), (/ny, nz, neofs/), no_eof )

!----------------------------------------------------------------------
!	... close the netcdf file
!----------------------------------------------------------------------
        iernc=nf90_close(ncid)     
        if(is_master) then
!SK     if(mpi_id.eq.0) then
         write(iulog,*) ' VAYsnoe ZKM:', no_zkm(1), ': ', no_zkm(nz)
         write(iulog,*) ' VAYsnoe MLT:', no_mlat(1), ': ', no_mlat(ny)
         write(iulog,*) ' VAYsnoe NO:', maxval(no_m), minval(no_m)
        endif
      end subroutine solar_readno_snoewx
!
      END MODULE IDEA_solarno_input
!=====================================================
      SUBROUTINE WAM_COOLNO1(np,nps,T,O,NO,QNO)   
!
!vay Oct 1  2015 Clean-up
!    Oct 28 2016 for new trunk
!      use   physcons, only : bz =>  con_boltz             
!-------------------------------------------------------------------------
! calculate NO cooling
!-------------------------------------------------------------------------
!  **
!  input:
!  T temperature profile K
!  O atomic oxygen number density profile m-3
!  NO nitric oxide number density profile m-3
!  output:
!  QNO: NO cooling rate J/m-3
!  **
      implicit none
      integer, intent(in):: np           ! numer of pressure levels
      integer, intent(in):: nps          ! pressure index to start
      real, intent(in)   :: O(np),NO(np) ! number density/m3        
      real, intent(in)   :: T(np)        ! temp (K)   
!out  
      real, intent(out)  :: QNO(np)    
!locals     
      real :: K10,HV,A10, G
      real :: A1,A2,A3,OM1,OM     
      real :: A10gHV, HVBZ, A23   
      real :: BZ
      integer i
!
!vay-2015/16
!                           ! HV=phot_e  = 3.726e-13_r8   at 5.3 mum (erg)
!
       A10=13.3             !  trans_prob = 13.3_r8       
       BZ=1.38E-23          ! boltzman
       K10=2.7E-17          ! k10*[O3P] vs 2.7 in WACCM  K10=3.6E-17 
                            ! decreased because cm3/s = 10(-11)
                            !                       m3/s    10(-17)                                       
       HV=3.726E-20         ! in Joules   1 erg is equal to 1.0E-7 joule.
       HVBZ = HV/BZ                           
       G=1.0                                 
!                                                     
       A2=5.4E-6*(1./(EXP(HVBZ/5800.)-1.))    
       A3=0.5*EXP(-HVBZ/247.5)
       A23 = A2+A3
       A10=13.3     
       A10gHV = A10*G*HV
!-------------------------------------------------
! K10/A10, HBVZ, A10gHV, A23 => idear_solar_init
!=================================================
       QNO(1:np) =0.
      do i=nps,np                                                   
        OM1=K10*O(i)          
        OM=OM1/(OM1+A10)
        A1=EXP(-HVBZ/T(i)) 
        QNO(i)=A10gHV*NO(i)*OM*(A1-A23)      ! should be in "J/kg/s" ~"J/s*[1/m3-NO]"
      enddo
!====================================
! deactivation term trans_prob = 13.3_r8                      HVBZ=2700
!
!   nocool(i,k) = -1.e-4_r8 * phot_e * trans_prob* exp(-2700._r8/t(i,k)) &                     
!                  * no_conc * (no_deact / (no_deact + trans_prob)) 
!
! [NO] 1/m3
! [OM] dimensionless
! [A10] -dimensionless
!      A1*A10*HV*[NO]/CP = [K/sec]
!      A1 -dimens
! units:     Hv*[NO]/Cp
!====================================
      return                                                            
      end  SUBROUTINE WAM_COOLNO1
!
      subroutine getno1d(levs,f107in,kpain,mlatrad,doy,alt,pr,n,am,no)
!
! Oct 2016 VAY: new interface  getno1d ... mozaic of bugs in 2012-14 versions
!               (a) rad/deg (b) pr(cb) but expected in Pa (c) no-NO for lat > 80deg
!               (d) dangerous interp-n and DATA for eofs etc...
!
! Mar 2018 Zhuxiao Li and Tzu-Wei Fang,read in the 24hr ave kp (kpa)
! from driving parameters file instead of reading kp.
!
!       use idea_solar_input, only : F107 => wf107_s
!       use idea_solar_input, only : kp => wkp_s
!       use idea_solar, only       : amno
       use idea_solarno_input, only : eof => no_eof, nom => no_m, z16 =>
     &                                no_zkm , lat33 => no_mlat
       use idea_solarno_input, only : no_ny33, no_nz16
       use idea_composition, only : r2d => R_2_D, d2r => dtr, twopi =>
     &                              pi2    !180/!pi & !pi/180.
       use idea_composition, only : con_nzero, amno
!
      implicit none
! input
      integer, intent(in)  :: levs         !number of pressure level
      integer, intent(in)  :: doy          ! day of year from 1 to 365   

      real,    intent(in)  :: f107in       !F10.7 index
      real,    intent(in)  :: kpain         ! 24hr average kp index
      real,    intent(in)  :: mlatrad      ! magnetic latitude in radians

      real,    intent(in)  :: alt(levs)    ! in km WAM
      real,    intent(in)  :: pr(levs)     ! in Pa
      real,    intent(in)  :: am(levs)     ! avg mass g/mol
      real,    intent(in)  :: n(levs)      !/m3 number density
! out
      real,    intent(out) :: no(levs)     ! number density of NO (/m3) 
!
! locals
!  we have         33, 16,7 now...
!      real :: eof(33,16,3),nom(33,16),z16(16), lat33(33)
      real :: mlat          ! degrees like in snoe data
      real ::  zm(no_nz16)
      real ::  dz(levs)
      real ::  f107, kpa
      real ::  dx, dl,m1,m2,m3,theta0,dec
!
      integer :: iref,kref(levs)
      integer :: i,k,il,k1,k2
      integer :: kup
!
       kpa =kpain
       f107 = f107in

       mlat = mlatrad *r2d

!       print *, kpa, f107, ' getno1d-VAY-mlat-deg ', mlat
       if(kpa .lt. 0.7)   kpa=0.7
       if(f107.lt.70.0) f107 = 70.0
!
! find interp latitude lat33 in degrees...
!     
! VAY-2016: mlat-bug
        do i=1,no_ny33-1
          if(mlat.gt.lat33(i).and.mlat.lt.lat33(i+1)) then
            iref=i
            dl=(mlat-lat33(i))/(lat33(i+1)-lat33(i))
            exit
          endif
        enddo
!vay oct-2016 add abs(mlat) > 80.
         if (mlat.le.lat33(1)) then 
            dl = 0.
            iref =1
         endif
         if (mlat.ge.lat33(no_ny33)) then 
            dl = 1.
            iref =no_ny33-1
         endif
!
!  snoe NO interpolated to model grid (molecules/cm^3)
!... eof1 - kpa
!     m1 =  kpa * 0.689254 - 1.53366
      m1 =  kpa * 0.785760 - 1.94262           !waccmx-2015
!... eof2 - declination
      theta0 = twopi/365.*float(doy - 1)
      dec = 0.006918                                                    
     &    - 0.399912 * cos(theta0)   + 0.070257 * sin(theta0)           
     &    - 0.006758 * cos(2.*theta0) + 0.000907 * sin(2.*theta0)       
     &    - 0.002697 * cos(3.*theta0) + 0.001480 * sin(3.*theta0)
      dec = dec * r2d   !180./3.1415927
      m2 = -.319782 + dec*(.0973109 + dec*(.00048981 - dec*.000103608))
!      m2 = -0.31978                                                     
!     &   + dec    * 0.097309                                            
!     &   + dec**2 * 0.00048979                                          
!     &   - dec**3 * 0.00010360
!... eof3 - f107
!      m3 =  alog10(f107) * 6.35777 - 13.8163
       m3 =  alog10(f107) * 6.44069 - 13.9832                 !waccmx-2015
!
!... zonal mean distrib. is sum of mean and eofs
! (dl)*no(i+1) + (1-dl)*no(i) interpolation to WAM "mlat"
!  WX-code     do k = 1,nlev
!  zm(:,k) = no_mean(:,k) - m1 * eofs(:,k,1) + m2 * eofs(:,k,2) - m3 * eofs(:,k,3) 
!              end do
    
        do k=1,no_nz16
          zm(k) =dl*(nom(iref+1,k)                              
     &        - m1 * eof(iref+1,k,1)                                
     &        + m2 * eof(iref+1,k,2)                                
     &        - m3 * eof(iref+1,k,3))+                              
     &         (1.-dl)*(nom(iref,k)                             
     &        - m1 * eof(iref,k,1)                                  
     &        + m2 * eof(iref,k,2)                                  
     &        - m3 * eof(iref,k,3))
          if (zm(k).le. 0.0) zm(k)=con_nzero
        enddo
        zm = zm*1.e6 ! zm transform fom cm-3 to m-3  OK... due to data UNITS
!
! vertical interp, from k1 to k2-1, extend k2 to levs, keep 
! cons 1 to k1-1
!

!interpolate
 
!        print *, ' getno-ZmNO ', maxval(zM), minval(zM)
 
        no(1:levs)=con_nzero
        CALL  interpol_wamz( no_nz16, z16, zm, levs, alt, NO, kup ) 

!        print *, ' getno-ZmNOi ', maxval(NO), minval(NO)
       
!
!extrapolate 
!
        do k=kup+1,levs
           dx =log(pr(k-1))-log(pr(k))
           no(k)=no(k-1)*n(k)/n(k-1)*                        
     &    exp(dx*(1.-.5*amno*(1./am(k-1)+1./am(k))))
        enddo
!
! extra-check for positive no
!
        do k=1,levs
          no(k)=max(no(k), con_nzero)
        enddo
!
! incease NO by a factor when kpa gt 5.0, same with CTIPe, based on CHAMP
! neutral density obervations 
!  
        do k=1,levs
          if (kpa.gt.5.0.and.kpa.le.6.0) then
            no(k) = no(k)*1.5
          elseif (kpa.gt.6.0.and.kpa.le.7.0) then
            no(k) = no(k)*2.5
          elseif (kpa.gt.7.0.and.kpa.le.8.0) then
            no(k) = no(k)*3.5
          elseif (kpa.gt.8.0.and.kpa.le.9.0) then
            no(k) = no(k)*4.5
          endif
        enddo
!
!
! check no
!       print *, ' getno-ZmNOi2 ', maxval(NO), minval(NO)
!        print *, 'Kup getno, iref', kup, iref, con_nzero
      return
      end subroutine getno1d
!
      subroutine interpol_wamz( nin, xin, yin, nout, xout, yout, kup )
!-----------------------------------------------------------------------
!	... linear interpolation in vertical
!           does not extrapolate, but repeats edge values
!-----------------------------------------------------------------------

      implicit none
!-----------------------------------------------------------------------
      integer,  intent(in) :: nin, nout
      real, intent(in)    :: xin(nin)   ! reverse 150 => 100
      real, intent(in)    :: yin(nin)
      real, intent(in)    :: xout(nout)
      real, intent(out)   :: yout(nout)
      integer,intent(out) ::  kup      ! zout(kout) > zin(1) = 150 km. zin(nin)=100.
      integer             :: dxin      ! top => bot dxin < 0; bot => top dxin > 0
!-----------------------------------------------------------------------
!	... local variables 
!-----------------------------------------------------------------------
      integer :: i, j
      kup = nout                       ! for extrap-n up

      dxin  = xin(2)-xin(1)
      IF (dxin < 0) THEN               ! snoe-no
      do j = 1,nout
       if( xout(j) .ge. xin(1) ) then
          yout(j) = yin(1)
          kup = min(j, kup)
       else   if( xout(j) .le. xin(nin) ) then 
               yout(j) = yin(nin)
       else
         do i = 1, nin-1
         if ((xout(j) >= xin(i+1)) .and. (xout(j) < xin(i)) )
     &   yout(j) =  yin(i) + (yin(i+1) - yin(i)) * (xout(j) - xin(i)) /
     &   (xin(i+1) - xin(i))
         end do
       end if
      end do
! dxin < 0
      ELSE
! dxin > 0
     
       do j = 1,nout
        if( xout(j) .le. xin(1) )   then
           yout(j) = yin(1)
        else  if( xout(j) .ge. xin(nin) ) then
            yout(j) = yin(nin)
            kup = min(j, kup)
        else
         do i = 1, nin-1
         if ((xout(j) .gt. xin(i)) .and. (xout(j) .lt. xin(i+1) ))
     & yout(j) = yin(i) + (yin(i+1) - yin(i)) * (xout(j) - xin(i)) / 
     & (xin(i+1) - xin(i))
       
         end do
        end if
      end do
!
      ENDIF
      end subroutine interpol_wamz
      subroutine interpol_wamz_down(nin,xin,yin,nout, xout, yout, down )
!-----------------------------------------------------------------------
!	... linear interpolation in vertical
!           does not extrapolate, but repeats edge values
!-----------------------------------------------------------------------

      implicit none
!-----------------------------------------------------------------------
      integer,  intent(in) :: nin, nout
      real, intent(in)    :: xin(nin)   ! reverse 150 => 100
      real, intent(in)    :: yin(nin)
      real, intent(in)    :: xout(nout)
      real, intent(in)    :: down      ! zeroes below
      real, intent(out)   :: yout(nout)
      integer             ::  kup      ! zout(kout) > zin(1) = 150 km. zin(nin)=100.
      integer             :: dxin      ! top => bot dxin < 0; bot => top dxin > 0
!-----------------------------------------------------------------------
!	... local variables 
!-----------------------------------------------------------------------
      integer :: i, j
      kup = nout                       ! for extrap-n up

      dxin  = xin(2)-xin(1)
      IF (dxin < 0) THEN               ! snoe-no
      do j = 1,nout
       if( xout(j) .gt. xin(1) ) then
          yout(j) = yin(1)
          kup = min(j, kup)
       else   if( xout(j) .lt. xin(nin) ) then 
               yout(j) = down
       else
         do i = 1, nin-1
         if ((xout(j) >= xin(i+1)) .and. (xout(j) <= xin(i)) )
     &   yout(j) =  yin(i) + (yin(i+1) - yin(i)) * (xout(j) - xin(i)) / 
     &   (xin(i+1) - xin(i))
         end do
       end if
      end do
! dxin < 0
      ELSE
! dxin > 0
     
       do j = 1,nout
        if( xout(j) .lt. xin(1) )   then
           yout(j) = down
        else  if( xout(j) .gt. xin(nin) ) then
            yout(j) = yin(nin)
            kup = min(j, kup)
        else
         do i = 1, nin-1
         if ((xout(j) .ge. xin(i)) .and. (xout(j) .le. xin(i+1) ))
     & yout(j) = yin(i) + (yin(i+1) - yin(i)) * (xout(j) - xin(i)) /
     & (xin(i+1) - xin(i))
       
         end do
        end if
      end do
!
      ENDIF
      end subroutine interpol_wamz_down
!
      subroutine z15toz(ain,levs, Z, aout,down)
! interpolate 15 pressure levels (from Tim's grid) to
! idea pressure grid pr(levs)
      use idea_composition, only : pr=> pr_idea
      implicit none
      integer, parameter :: np=15      !number of pressure levels of input
      integer, intent(in) :: levs      !number of pressure levels of output 
      real,    intent(in) :: Z(levs)   !model -log(pressure) grid
      real,    intent(in) :: ain(np)   !input field in 15 pressure grid
      real,    intent(in) :: down      !field value below 1.0376Pa => surface
      real,    intent(out):: aout(levs)!output in levs pressure grid
!local variable
      real p15(np),z15(np),dz
      integer kref,k,i
!
      do k=1,np
        p15(k)=1.0376*exp(1.-k)
        z15(k)=-log(p15(k))
      enddo
!
! grids
! interpolation
!
      do k=1,levs
        do i=1,np-1
          if(z(k).ge.z15(i).and.z(k).le.z15(i+1)) then
          dz=(z(k)-z15(i))/(z15(i+1)-z15(i))*(ain(i+1)-ain(i))
          aout(k)=ain(i) +dz
          endif
        enddo   
        if(z(k).lt.z15(1))  aout(k)=down    ! zero
        if(z(k).gt.z15(15)) aout(k)=ain(15)
    
      enddo
      return
      end
      subroutine z63toz(ain, levs, Z, aout,down)
! interpolate 63 pressure levels (from Tim's grid) to
! idea pressure grid pr(levs)
      use idea_composition, only : pr=> pr_idea
      implicit none
      integer, parameter  :: np=63     !number of pressure levels of input
      integer, intent(in) :: levs      !number of pressure levels of output 
      real,    intent(in) :: Z(levs)   ! model grid
      real,    intent(in) :: ain(np)   !input field in 63 pressure grid
      real,    intent(in) :: down      !field value under 6.9Pa
      real,    intent(out):: aout(levs)!output in levs pressure grid
!local variable
      real p63(np),z63(np),dz
      integer kref,k,i
!
      DATA p63/6.90775528,  6.57442194, 
     &   6.24108859,  5.90775525,  5.57442191,
     &   5.24108856,  4.90775522,  4.57442188,
     &   4.24108853,  3.90775519,  3.57442185,
     &   3.2410885,   2.90775516,  2.57442182,
     &   2.24108847,  1.90775513,  1.57442179,
     &   1.24108844,  0.9077551,   0.574421757,
     &   0.241088414, -0.0922449296,-0.425578273,
     &   -0.758911616,-1.09224496,  -1.4255783, 
     &   -1.75891165, -2.09224499,  -2.42557833,
     &   -2.75891168, -3.09224502,  -3.42557836,
     &   -3.75891171, -4.09224505,  -4.42557839,
     &   -4.75891174, -5.09224508,  -5.42557842,
     &   -5.75891177, -6.09224511,  -6.42557845,
     &   -6.75891179, -7.09224514,  -7.42557848,
     &   -7.75891182, -8.09224517,  -8.42557851,
     &   -8.75891185, -9.0922452,   -9.42557854,
     &   -9.75891188, -10.0922452,  -10.4255786,
     &   -10.7589119, -11.0922453,  -11.4255786,
     &   -11.7589119, -12.0922453,  -12.4255786,
     &   -12.758912,  -13.0922453,  -13.4255787,
     &   -13.758912/
      do k=1,np
        z63(k)=-p63(k)
      enddo
!
      do k=1,levs
        do i=1,np-1
          if(z(k).ge.z63(i).and.z(k).le.z63(i+1)) then
          dz=(z(k)-z63(i))/(z63(i+1)-z63(i))*(ain(i+1)-ain(i))
           aout(k)=dz+ain(i)
          endif
        enddo
        if(z(k).lt.z63(1))  aout(k)=down
        if(z(k).gt.z63(63)) aout(k)=ain(63)
      
      enddo
      return
      end
!
