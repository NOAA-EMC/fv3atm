module update_ca

use atmosphere_mod,    only: atmosphere_scalar_field_halo
use mersenne_twister,  only: random_setseed,random_gauss,random_stat
use fv_mp_mod, only : mp_reduce_sum,mp_bcst,mp_reduce_min,mp_reduce_max

implicit none

!L. Bengtsson 2017-06
!Evolve the cellular automata in time


contains

subroutine update_cells(kstep,nca,nxc,nyc,nxch,nych,nlon,nlat,CA,ca_plumes,iini,ilives, &
                        nlives,ncells,nfracseed,nseed,nthresh,ca_global, &
                        ca_sgs,nspinup,condition,vertvelhigh,nf,nca_plumes)

implicit none

integer, intent(in) :: kstep,nxc,nyc,nlon,nlat,nxch,nych,nca
integer, intent(in) :: iini(nxc,nyc,nca), ilives(nxc,nyc)
real, intent(in) :: condition(nxc,nyc),vertvelhigh(nxc,nyc)
real, intent(out) :: CA(nlon,nlat)
integer,intent(out) :: ca_plumes(nlon,nlat)
integer, intent(in) :: nlives, ncells, nseed, nspinup, nf
real, intent(in) :: nfracseed, nthresh
logical,intent(in) :: nca_plumes
real, dimension(nlon,nlat) :: frac
integer,allocatable,save :: board(:,:,:), lives(:,:,:)
integer,allocatable :: V(:),L(:)
integer :: inci, incj, i, j, k, iii,sub,spinup,it,halo,k_in,isize,jsize
integer :: haloh, ih, jh,lives_max,kend
real :: thresh,threshc,threshk,wp_max,wp_min,mthresh,kthresh
real, allocatable :: field_in(:,:),board_halo(:,:,:), Wpert_halo(:,:,:),Cpert_halo(:,:,:)
integer, dimension(nxc,nyc) :: neighbours, birth, newlives
integer, dimension(nxc,nyc) :: neg, newcell, oldlives, newval,temp
integer, dimension(ncells,ncells) :: onegrid
real,dimension(nxc,nyc) :: normlives
logical :: ca_global, ca_sgs
real :: Wpert(nxc,nyc),Cpert(nxc,nyc)

!SGS parameters:                                                                                                               
integer(8) :: count, count_rate, count_max, count_trunc
integer(8) :: iscale = 10000000000                                                                               
integer :: count5, count6
type(random_stat) :: rstate
real :: dt, timescale, sigma, skew, kurt, acorr, gamma
real :: E_sq2, E2, g2, B_sq2, B2, sqrtdt,flamx2, tmp1, tmp1a
real, dimension(nxc,nyc) :: NOISE_A, NOISE_B, g2D
real, dimension(nxc*nyc) :: noise1D2, noise1D1
real, allocatable, save :: sgs1(:,:,:),sgs2(:,:,:)
integer, dimension(nxch,nych,1) :: M_halo


!-------------------------------------------------------------------------------------------------
halo=1
isize=nlon+2*halo
jsize=nlat+2*halo
k_in=1

 if (.not. allocated(board))then
 allocate(board(nxc,nyc,nca))
 endif
 if (.not. allocated(lives))then
 allocate(lives(nxc,nyc,nca))
 endif
 if (.not. allocated(sgs1))then
 allocate(sgs1(nxc,nyc,nca))
 endif
 if (.not. allocated(sgs2))then
 allocate(sgs2(nxc,nyc,nca))
 endif
 if (.not. allocated(field_in))then
 allocate(field_in(nxc*nyc,1))
 endif
 if(.not. allocated(board_halo))then                                                                      
 allocate(board_halo(nxch,nych,1))   
 endif
 if(.not. allocated(Wpert_halo))then
 allocate(Wpert_halo(nxch,nych,1))
 endif
 if(.not. allocated(Cpert_halo))then
 allocate(Cpert_halo(nxch,nych,1))
 endif
 
 if(ca_sgs)then  

  if(kstep <= 1)then
 
  do j=1,nyc
   do i=1,nxc
    board(i,j,nf) = 0
    lives(i,j,nf) = 0
   enddo
  enddo
 
  endif

  if(kstep == 2)then !Initiate CA at kstep 2 as physics field is empty at 0 and 1.
 
  do j=1,nyc
   do i=1,nxc
    board(i,j,nf) = iini(i,j,nf)
    lives(i,j,nf) = ilives(i,j)*iini(i,j,nf)
   enddo
  enddo
 
  endif

 !Seed with new CA cells at each nseed step
  if(mod(kstep,nseed) == 0 .and. kstep >= 2)then
 
   do j=1,nyc
    do i=1,nxc
     board(i,j,nf) = iini(i,j,nf)                                                                              
     lives(i,j,nf) = ilives(i,j)*iini(i,j,nf)                                                                             
    enddo
   enddo
 
  endif
 
  if(kstep == 2)then
  spinup=nspinup
  else
  spinup = 1
  endif

 else !ca_global

  if(kstep == 0)then

   do j=1,nyc
    do i=1,nxc
     board(i,j,nf) = iini(i,j,nf)
     lives(i,j,nf) = ilives(i,j)*iini(i,j,nf)
    enddo
   enddo

  endif

 !Seed with new CA cells at each nseed step                                                                                                                                                                                  
  if(mod(kstep,nseed) == 0)then

   do j=1,nyc
    do i=1,nxc
     board(i,j,nf) = iini(i,j,nf)
     lives(i,j,nf) = ilives(i,j)*iini(i,j,nf)
    enddo
   enddo

  endif

  if(kstep == 0)then
  spinup=nspinup
  else
  spinup = 1
  endif
 

 endif !sgs/global

!Step 1 - Solve the stochastic gaussian skewed SGS equation in order to generate
!perturbations to the grid mean model fields.

if (ca_sgs) then
!Compute the SGS and perturb the vertical velocity field:                                                                                                          
!Read these values in from namelist, guided from LES data:                                                                                                                            

dt=900.                                                                                                                                                                   
timescale=21600.0
sigma=0.8
skew=0.8
kurt=2.0
acorr=exp(-dt/timescale)
gamma=-log(acorr)/dt

!calculate coeffients for SGS with auto-correlation and use
!these for predictor-correcting time stepping                                                                                                                                                                  
E_sq2=2.0*gamma/3.0*(kurt-3.0/2.0*skew**2.)/(kurt-skew**2.+2.)
E2=sqrt(E_sq2)
g2D=skew*sigma*(gamma-E_sq2)/(2.0*E2)
if(kstep>=2)then
g2D=g2D+E2*condition(:,:)
endif
B_sq2=2.0*sigma**2*(gamma-E_sq2/2.0-(gamma-E_sq2)**2*skew**2/(8.0*E_sq2))
B2=sqrt(B_sq2)
B2=0.

sqrtdt=sqrt(dt)
flamx2=0.5*E_sq2+gamma

endif

do it=1,spinup

if (ca_sgs) then

 !Random seed for SGS                                                                                                        
 noise1D1 = 0.0
 noise1D2 = 0.0
                                                                               
 call system_clock(count, count_rate, count_max)
 count_trunc = iscale*(count/iscale)
 count5 = count - count_trunc
 count6=count5+9827
 !broadcast to all tasks                                                                                                                                                                                        
 call mp_bcst(count5)
 call mp_bcst(count6)

 call random_setseed(count5,rstate)
 call random_gauss(noise1D1,rstate)

 call random_setseed(count6,rstate)
 call random_gauss(noise1D2,rstate)
 !Put on 2D:                                                                                                                                                                                                    
 do j=1,nyc
  do i=1,nxc
  NOISE_A(i,j)=noise1D1(i+(j-1)*nxc)
  NOISE_B(i,j)=noise1D2(i+(j-1)*nxc)
  enddo
 enddo

tmp1=0.
tmp1a=0.
Wpert=0.
Cpert=0.

if(kstep == 0)then
 do j=1,nyc
  do i=1,nxc
  sgs1(i,j,nf)=sigma*NOISE_A(i,j)
  sgs2(i,j,nf)=sigma*NOISE_B(i,j)
  Cpert(i,j)=sgs1(i,j,nf)
  Wpert(i,j)=sgs2(i,j,nf)
  enddo
 enddo

else

 do j=1,nyc
  do i=1,nxc
      tmp1=sgs1(i,j,nf)-(flamx2*sgs1(i,j,nf)+0.5*E2*g2D(i,j))*dt + (B2*NOISE_A(i,j)*sqrtdt) &
          + (g2D(i,j) + E2 * sgs1(i,j,nf))* NOISE_B(i,j)*sqrtdt
      tmp1a=(tmp1+sgs1(i,j,nf))*0.5
      sgs1(i,j,nf)=sgs1(i,j,nf)-(flamx2*tmp1a+0.5*E2*g2D(i,j))*dt + (B2*NOISE_A(i,j)*sqrtdt) &
          + (g2D(i,j) + E2 * tmp1a  )* NOISE_B(i,j)*sqrtdt

      Cpert(i,j)=condition(i,j)*(1.0 + sgs1(i,j,nf))

  enddo
 enddo

 do j=1,nyc
  do i=1,nxc
      tmp1=sgs2(i,j,nf)-(flamx2*sgs2(i,j,nf)+0.5*E2*g2D(i,j))*dt + (B2*NOISE_A(i,j)*sqrtdt) &
          + (g2D(i,j) + E2 * sgs2(i,j,nf))* NOISE_B(i,j)*sqrtdt
      tmp1a=(tmp1+sgs2(i,j,nf))*0.5

      sgs2(i,j,nf)=sgs2(i,j,nf)-(flamx2*tmp1a+0.5*E2*g2D(i,j))*dt + (B2*NOISE_A(i,j)*sqrtdt) &
          + (g2D(i,j) + E2 * tmp1a  )* NOISE_B(i,j)*sqrtdt

      Wpert(i,j)=vertvelhigh(i,j)*(1.0 + sgs2(i,j,nf))                                                                                                                                                                                               
  enddo
 enddo

endif

endif !ca sgs true

 
!Step 2 - Initialize variables to 0 and extract the halo
 
 neighbours=0
 birth=0
 newlives=0
 neg=0
 newcell=0
 oldlives=0
 newval=0
 frac=0
 board_halo=0
 Wpert_halo=0
 Cpert_halo=0
 field_in=0

!The input to scalar_field_halo needs to be 1D.                                                          
!take the updated board fields and extract the halo
! in order to have updated values in the halo region. 

 if(ca_global)then 
  do j=1,nyc                                                                        
   do i=1,nxc                                                                                           
   field_in(i+(j-1)*nxc,1)=board(i,j,nf)                                                                
   enddo                                                                                                
  enddo                                                                                                  
!Step 3 - Extract the halo                                              
  call atmosphere_scalar_field_halo(board_halo,halo,nxch,nych,k_in,field_in)
 endif
 

 if(ca_sgs)then
  field_in=0
  do j=1,nyc
   do i=1,nxc
   field_in(i+(j-1)*nxc,1)=Wpert(i,j)
   enddo
  enddo

  call atmosphere_scalar_field_halo(Wpert_halo,halo,nxch,nych,k_in,field_in)  

  field_in=0
  do j=1,nyc
   do i=1,nxc
   field_in(i+(j-1)*nxc,1)=Cpert(i,j)
   enddo
  enddo

  call atmosphere_scalar_field_halo(Cpert_halo,halo,nxch,nych,k_in,field_in)
 endif !sgs




!Step 4 - Compute the neighbourhood
if(ca_sgs)then !SGSmethod

 !Count the number of neighbours where perturbed massflux is larger than 
 !a threshold


if(nf==1)then      !Deep convection
 M_halo = 0
 do j=1,nych
  do i=1,nxch
  if(Wpert_halo(i,j,1) < nthresh .and. Cpert_halo(i,j,1) > nthresh)then
  M_halo(i,j,1) = 1
  endif
  enddo
 enddo
 elseif(nf==2)then  !Shallow convection
  M_halo=0
  do j=1,nych
  do i=1,nxch
  if(Wpert_halo(i,j,1) < nthresh .and. Cpert_halo(i,j,1)>nthresh)then
  M_halo(i,j,1) = 1
  endif
  enddo
 enddo
 elseif(nf==3)then  !Turbulence
  M_halo = 0
 do j=1,nych
  do i=1,nxch
  if(Wpert_halo(i,j,1) < nthresh .and. Cpert_halo(i,j,1) > nthresh)then
  M_halo(i,j,1) = 1
  endif
  enddo
 enddo
 elseif(nf==4)then !Radiation
 M_halo = 0
 do j=1,nych
  do i=1,nxch
  if(Cpert_halo(i,j,1)>nthresh)then
  M_halo(i,j,1) = 1
  endif
  enddo
 enddo
 else              !nf=5 Microphysics
 M_halo = 0
 do j=1,nych
  do i=1,nxch
  if(Wpert_halo(i,j,1) < nthresh .and. Cpert_halo(i,j,1)>nthresh)then
  M_halo(i,j,1) = 1
  endif
  enddo
 enddo
endif !nf

 do j=1,nyc
  do i=1,nxc
     ih=i+halo
     jh=j+halo
     neighbours(i,j)=M_halo(ih-1,jh-1,1)+M_halo(ih-1,jh,1)+ &
                        M_halo(ih-1,jh+1,1)+M_halo(ih,jh+1,1)+M_halo(ih+1,jh+1,1)+&
                        M_halo(ih+1,jh,1)+M_halo(ih+1,jh-1,1)+M_halo(ih,jh-1,1)
  enddo
 enddo

 
!CA stand alone method
else !global

  do j=1,nyc
     do i=1,nxc
        ih=i+halo
        jh=j+halo
        neighbours(i,j)=board_halo(ih-1,jh-1,1)+board_halo(ih-1,jh,1)+ &
                        board_halo(ih-1,jh+1,1)+board_halo(ih,jh+1,1)+board_halo(ih+1,jh+1,1)+&
                        board_halo(ih+1,jh,1)+board_halo(ih+1,jh-1,1)+board_halo(ih,jh-1,1)
     enddo
  enddo

endif  !sgs/global
! Step 5 - Check rules; the birth condition differs between SGS and GOL method

if(ca_sgs)then !SGS

 if(nf==1)then
 do j=1,nyc
   do i=1,nxc
     if((Wpert(i,j) < nthresh .and. Cpert_halo(i,j,1) > nthresh) .or. neighbours(i,j)==2 .or. neighbours(i,j)==3)then
     birth(i,j)=1
     endif
   enddo
  enddo
 elseif(nf==2)then
 do j=1,nyc
   do i=1,nxc
     if((Wpert(i,j) < nthresh .and. Cpert_halo(i,j,1) > nthresh) .or. neighbours(i,j)==2 .or. neighbours(i,j)==3)then
     birth(i,j)=1
     endif
   enddo
  enddo
 elseif(nf==3)then
 do j=1,nyc
   do i=1,nxc
     if((Wpert(i,j) < nthresh .and. Cpert_halo(i,j,1) > nthresh) .or. neighbours(i,j)==2 .or. neighbours(i,j)==3)then
     birth(i,j)=1
     endif
   enddo
  enddo
 elseif(nf==4)then
 do j=1,nyc
   do i=1,nxc
     if(Cpert(i,j) > nthresh .or. neighbours(i,j)==2 .or. neighbours(i,j)==3)then
     birth(i,j)=1
     endif
   enddo
  enddo
 else !nf=5
 do j=1,nyc
   do i=1,nxc
     if((Wpert(i,j) < nthresh .and. Cpert_halo(i,j,1) > nthresh) .or. neighbours(i,j)==2 .or. neighbours(i,j)==3)then
     birth(i,j)=1
     endif
   enddo
  enddo
 endif

else !GOL

  do j=1,nyc
   do i=1,nxc
     if(neighbours(i,j)==2 .or. neighbours(i,j)==3)then
     birth(i,j)=1
     endif
   enddo
  enddo

endif

  do j=1,nyc
   do i=1,nxc
     if(neighbours(i,j).ne.4 .or. neighbours(i,j).ne.5)then
     lives(i,j,nf)=lives(i,j,nf) - 1
     endif
   enddo
  enddo

  do j=1,nyc
   do i=1,nxc
   if(lives(i,j,nf)<0)then
     lives(i,j,nf)=0
   endif
   enddo
  enddo

  do j=1,nyc
   do i=1,nxc
    if(birth(i,j)==1 .and. lives(i,j,nf)==0)then
    newcell(i,j)=1
    endif
   enddo
  enddo

 if(ca_sgs)then
  do j=1,nyc
   do i=1,nxc
    lives(i,j,nf)=lives(i,j,nf)+newcell(i,j)*ilives(i,j)
   enddo
  enddo
 else !GOL
  do j=1,nyc
   do i=1,nxc
    lives(i,j,nf)=lives(i,j,nf)+newcell(i,j)*ilives(i,j)!*nlives
   enddo
  enddo
 endif

   do j=1,nyc
   do i=1,nxc
    if(neighbours(i,j)==3 .or. (board(i,j,nf)==1 .and. neighbours(i,j)==2))then
    board(i,j,nf)=1
    else
    board(i,j,nf)=0
    endif
   enddo
  enddo

enddo !spinup

!COARSE-GRAIN BACK TO NWP GRID
 
  inci=ncells
  incj=ncells
  sub=ncells-1
  DO j=1,nlat
     DO i=1,nlon
        frac(i,j)=(SUM(lives(inci-sub:inci,incj-sub:incj,nf)))/(ncells*ncells)
        inci=inci+ncells
     ENDDO
     inci=ncells
     incj=incj+ncells
  ENDDO

lives_max=maxval(ilives)
call mp_reduce_max(lives_max)

  if(ca_sgs)then
   CA(:,:) = (frac(:,:)/lives_max)
  else !global
   CA(:,:) = (frac(:,:)/real(nlives))
  endif


if(nca_plumes) then
!COMPUTE NUMBER OF CLUSTERS (CONVECTIVE PLUMES) IN EACH CA-CELL
!Note, at the moment we only use the count of the plumes found in a grid-cell
!In the future the routine "plumes" can also be used to give the size of 
!each individual plume for better coupling to the convection scheme.

  temp=0
  do j=1,nyc
   do i=1,nxc
     if(lives(i,j,1) > 0)then
      temp(i,j)=1  
     endif
   enddo
  enddo

  kend=ceiling((ncells*ncells)/2.)
  if (.not. allocated(V))then
  allocate(V(kend))
  endif
  if (.not. allocated(L))then
  allocate(L(kend))
  endif
  
  ca_plumes(:,:)=0
  inci=ncells
  incj=ncells
  sub=ncells-1
  DO j=1,nlat
     DO i=1,nlon
        onegrid(1:ncells,1:ncells)=temp(inci-sub:inci,incj-sub:incj)
        call plumes(V,L,onegrid,ncells,ncells,kend)
        do k=1,kend
         if (V(k) == 1)then
         ca_plumes(i,j)=ca_plumes(i,j)+L(k)
         endif
        enddo
        inci=inci+ncells
     ENDDO
     inci=ncells
     incj=incj+ncells
  ENDDO

else

ca_plumes(:,:)=0.

endif ! nca_plumes

end subroutine update_cells

end module update_ca
