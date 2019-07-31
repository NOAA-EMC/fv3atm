subroutine cellular_automata(kstep,Statein,Coupling,Diag,nblks,nlev, &
            nca,ncells,nlives,nfracseed,nseed,nthresh,ca_global,ca_sgs,iseed_ca, &
            ca_smooth,nspinup,blocksize)

use machine
use update_ca,         only: update_cells
use atmosphere_mod,    only: atmosphere_resolution, atmosphere_domain, &
                             atmosphere_scalar_field_halo, atmosphere_control_data
use mersenne_twister,  only: random_setseed,random_gauss,random_stat,random_number
use GFS_typedefs,      only: GFS_Coupling_type, GFS_diag_type, GFS_statein_type
use mpp_domains_mod,   only: domain2D
use block_control_mod, only: block_control_type, define_blocks_packed
use fv_mp_mod,         only : mp_reduce_sum,mp_bcst,mp_reduce_max,is_master


implicit none

!L.Bengtsson, 2017-06

!This program evolves a cellular automaton uniform over the globe given
!the flag ca_global, if instead ca_sgs is .true. it evolves a cellular automata conditioned on 
!perturbed grid-box mean field. The perturbations to the mean field are given by a 
!stochastic gaussian skewed (SGS) distribution.

!If ca_global is .true. it weighs the number of ca (nca) together to produce 1 output pattern
!If instead ca_sgs is given, it produces nca ca:                                                                                                                                                                            
! 1 CA_DEEP = deep convection                                                                                                                                                                                               
! 2 CA_SHAL = shallow convection                                                                                                                                                                                            
! 3 CA_TURB = turbulence                                                                                                                                                                                                    
! 4 CA_RAD = radiation                                                                                                                                                                                                      
! 5 CA_MICRO = microphysics 

!PLEASE NOTE: This is considered to be version 0 of the cellular automata code for FV3GFS, some functionally 
!is missing/limited. 

integer,intent(in) :: kstep,ncells,nca,nlives,nseed,iseed_ca,nspinup
real,intent(in) :: nfracseed,nthresh
logical,intent(in) :: ca_global, ca_sgs, ca_smooth
integer, intent(in) :: nblks,nlev,blocksize
type(GFS_coupling_type),intent(inout) :: Coupling(nblks)
type(GFS_diag_type),intent(inout) :: Diag(nblks)
type(GFS_statein_type),intent(in) :: Statein(nblks)
type(block_control_type)          :: Atm_block
type(random_stat) :: rstate
integer :: nlon, nlat, isize,jsize,nf,k350,k850
integer :: inci, incj, nxc, nyc, nxch, nych
integer :: halo, k_in, i, j, k, iec, jec, isc, jsc
integer :: seed, ierr7,blk, ix, iix, count4,ih,jh
integer :: blocksz,levs,ra,rb,rc
integer(8) :: count, count_rate, count_max, count_trunc
integer(8) :: iscale = 10000000000
integer, allocatable :: iini(:,:,:),ilives(:,:),ca_plumes(:,:)
real(kind=kind_phys), allocatable :: field_out(:,:,:), field_in(:,:),field_smooth(:,:),Detfield(:,:,:)
real(kind=kind_phys), allocatable :: omega(:,:,:),pressure(:,:,:),cloud(:,:),humidity(:,:)
real(kind=kind_phys), allocatable :: vertvelsum(:,:),vertvelmean(:,:),dp(:,:,:),surfp(:,:)
real(kind=kind_phys), allocatable :: CA(:,:),CAstore(:,:),CAavg(:,:),condition(:,:),rho(:,:),cape(:,:)
real(kind=kind_phys), allocatable :: CA_DEEP(:,:),CA_TURB(:,:),CA_RAD(:,:),CA_MICRO(:,:),CA_SHAL(:,:)
real(kind=kind_phys), allocatable :: noise1D(:),vertvelhigh(:,:),noise(:,:,:)
real(kind=kind_phys) :: psum,csum,CAmean,sq_diff,CAstdv,count1,alpha
real(kind=kind_phys) :: Detmax(nca),Detmean(nca),phi,stdev,delt
logical,save         :: block_message=.true.
logical              :: nca_plumes

!nca         :: switch for number of cellular automata to be used.
!ca_global   :: switch for global cellular automata
!ca_sgs      :: switch for cellular automata conditioned on SGS perturbed vertvel. 
!nfracseed   :: switch for number of random cells initially seeded
!nlives      :: switch for maximum number of lives a cell can have
!nspinup     :: switch for number of itterations to spin up the ca
!ncells      :: switch for higher resolution grid e.g ncells=4 
!               gives 4x4 times the FV3 model grid resolution.                
!ca_smooth   :: switch to smooth the cellular automata
!nthresh     :: threshold of perturbed vertical velocity used in case of sgs
!nca_plumes   :: compute number of CA-cells ("plumes") within a NWP gridbox.

k350 = 29
k850 = 13
alpha = 1.5
ra = 201
rb = 2147483648
rc = 4294967296

halo=1
k_in=1

nca_plumes = .false.
!----------------------------------------------------------------------------
! Get information about the compute domain, allocate fields on this
! domain

! WRITE(*,*)'Entering cellular automata calculations'

! Some security checks for namelist combinations:
 if(nca > 5)then
 write(0,*)'Namelist option nca cannot be larger than 5 - exiting'
 stop
 endif

 if(ca_global .and. ca_sgs)then
 write(0,*)'Namelist options ca_global and ca_sgs cannot both be true - exiting'
 stop
 endif

 if(ca_sgs .and. ca_smooth)then
 write(0,*)'Currently ca_smooth does not work with ca_sgs - exiting'
 stop
 endif

 call atmosphere_resolution (nlon, nlat, global=.false.)
 isize=nlon+2*halo
 jsize=nlat+2*halo
 !nlon,nlat is the compute domain - without haloes       
 !mlon,mlat is the cubed-sphere tile size. 

 inci=ncells
 incj=ncells
 
 nxc=nlon*ncells
 nyc=nlat*ncells
 
 nxch=nxc+2*halo
 nych=nyc+2*halo


 !Allocate fields:

 allocate(cloud(nlon,nlat))
 allocate(omega(nlon,nlat,k350))
 allocate(pressure(nlon,nlat,k350))
 allocate(humidity(nlon,nlat))
 allocate(dp(nlon,nlat,k350))
 allocate(rho(nlon,nlat))
 allocate(surfp(nlon,nlat))
 allocate(vertvelmean(nlon,nlat))
 allocate(vertvelsum(nlon,nlat))
 allocate(field_in(nlon*nlat,1))
 allocate(field_out(isize,jsize,1))
 allocate(field_smooth(nlon,nlat))
 allocate(iini(nxc,nyc,nca))
 allocate(ilives(nxc,nyc))
 allocate(vertvelhigh(nxc,nyc))
 allocate(condition(nxc,nyc))
 allocate(cape(nlon,nlat))
 allocate(Detfield(nlon,nlat,nca))
 allocate(CA(nlon,nlat))
 allocate(ca_plumes(nlon,nlat))
 allocate(CAstore(nlon,nlat))
 allocate(CAavg(nlon,nlat))
 allocate(CA_TURB(nlon,nlat))
 allocate(CA_RAD(nlon,nlat))
 allocate(CA_DEEP(nlon,nlat)) 
 allocate(CA_MICRO(nlon,nlat))
 allocate(CA_SHAL(nlon,nlat))
 allocate(noise(nxc,nyc,nca))
 allocate(noise1D(nxc*nyc))
  
 !Initialize:
 Detfield(:,:,:)=0.
 vertvelmean(:,:) =0.
 vertvelsum(:,:)=0.
 cloud(:,:)=0. 
 humidity(:,:)=0.
 condition(:,:)=0.
 cape(:,:)=0.
 vertvelhigh(:,:)=0.
 noise(:,:,:) = 0.0 
 noise1D(:) = 0.0
 iini(:,:,:) = 0
 ilives(:,:) = 0
 CA(:,:) = 0.0
 CAavg(:,:) = 0.0
 ca_plumes(:,:) = 0
 CA_TURB(:,:) = 0.0
 CA_RAD(:,:) = 0.0
 CA_DEEP(:,:) = 0.0
 CA_MICRO(:,:) = 0.0
 CA_SHAL(:,:) = 0.0
 
!Put the blocks of model fields into a 2d array
 levs=nlev
 blocksz=blocksize

 call atmosphere_control_data(isc, iec, jsc, jec, levs)
 call define_blocks_packed('cellular_automata', Atm_block, isc, iec, jsc, jec, levs, &
                              blocksz, block_message)

  isc = Atm_block%isc
  iec = Atm_block%iec
  jsc = Atm_block%jsc
  jec = Atm_block%jec 

 do blk = 1,Atm_block%nblks
  do ix = 1, Atm_block%blksz(blk)
      i = Atm_block%index(blk)%ii(ix) - isc + 1
      j = Atm_block%index(blk)%jj(ix) - jsc + 1
      cape(i,j) = Coupling(blk)%cape(ix) 
      surfp(i,j) = Statein(blk)%pgr(ix)
      humidity(i,j)=Statein(blk)%qgrs(ix,k850,1) !about 850 hpa
      do k = 1,k350 !Lower troposphere: level k350 is about 350hPa 
      omega(i,j,k) = Statein(blk)%vvl(ix,k) ! layer mean vertical velocity in pa/sec
      pressure(i,j,k) = Statein(blk)%prsl(ix,k) ! layer mean pressure in Pa
      enddo
  enddo
 enddo

!Compute layer averaged vertical velocity (Pa/s)
 vertvelsum=0.
 vertvelmean=0.
 do j=1,nlat
  do i =1,nlon
    dp(i,j,1)=(surfp(i,j)-pressure(i,j,1))
    do k=2,k350
     dp(i,j,k)=(pressure(i,j,k-1)-pressure(i,j,k))
    enddo
    count1=0.
    do k=1,k350
     count1=count1+1.
     vertvelsum(i,j)=vertvelsum(i,j)+(omega(i,j,k)*dp(i,j,k))
   enddo
  enddo
 enddo

 do j=1,nlat
  do i=1,nlon
   vertvelmean(i,j)=vertvelsum(i,j)/(surfp(i,j)-pressure(i,j,k350))
  enddo
 enddo

!Generate random number, following stochastic physics code:
do nf=1,nca

  if (iseed_ca == 0) then
    ! generate a random seed from system clock and ens member number
    call system_clock(count, count_rate, count_max)
    ! iseed is elapsed time since unix epoch began (secs)
    ! truncate to 4 byte integer
    count_trunc = iscale*(count/iscale)
    count4 = count - count_trunc + nf*ra 
  else
    ! don't rely on compiler to truncate integer(8) to integer(4) on
    ! overflow, do wrap around explicitly.
    count4 = mod(iseed_ca + rb, rc) - rb + nf*ra
  endif

 !Set seed (to be the same) on all tasks. Save random state.
 call random_setseed(count4,rstate)
 call random_number(noise1D,rstate)
 !Put on 2D:
 do j=1,nyc
  do i=1,nxc
  noise(i,j,nf)=noise1D(i+(j-1)*nxc)
  enddo
 enddo


!Initiate the cellular automaton with random numbers larger than nfracseed
 
   do j = 1,nyc
    do i = 1,nxc
     if (noise(i,j,nf) > nfracseed ) then
      iini(i,j,nf)=1
     else
      iini(i,j,nf)=0
    endif
    enddo
   enddo

 enddo !nf
 
!In case we want to condition the cellular automaton on a large scale field
!we here set the "condition" variable to a different model field depending
!on nf. (this is not used if ca_global = .true.)

CAstore = 0.

do nf=1,nca !update each ca
 
 if(ca_sgs)then

  if(nf==1)then 
  inci=ncells
  incj=ncells
  do j=1,nyc
   do i=1,nxc
     condition(i,j)=cape(inci/ncells,incj/ncells) 
     if(i.eq.inci)then
     inci=inci+ncells
     endif
   enddo
   inci=ncells
   if(j.eq.incj)then
   incj=incj+ncells
   endif
  enddo

   do j = 1,nyc
    do i = 1,nxc
      ilives(i,j)=int(real(nlives)*alpha*noise(i,j,nf))
    enddo
   enddo

 elseif(nf==2)then
  inci=ncells
  incj=ncells
  do j=1,nyc
   do i=1,nxc
     condition(i,j)=humidity(inci/ncells,incj/ncells)
     if(i.eq.inci)then
     inci=inci+ncells
     endif
   enddo
   inci=ncells
   if(j.eq.incj)then
   incj=incj+ncells
   endif
  enddo

   do j = 1,nyc
    do i = 1,nxc
      ilives(i,j)=int(real(nlives)*alpha*noise(i,j,nf))                                                                                                                                               
    enddo
   enddo

 elseif(nf==3)then
  inci=ncells
  incj=ncells
  do j=1,nyc
   do i=1,nxc
     condition(i,j)=cape(inci/ncells,incj/ncells)
     if(i.eq.inci)then
     inci=inci+ncells
     endif
   enddo
   inci=ncells
   if(j.eq.incj)then
   incj=incj+ncells
   endif
  enddo

  do j = 1,nyc
    do i = 1,nxc
      ilives(i,j)=int(real(nlives)*alpha*noise(i,j,nf))
    enddo
   enddo

 elseif(nf==4)then
  inci=ncells
  incj=ncells
  do j=1,nyc
   do i=1,nxc
     condition(i,j)=cape(inci/ncells,incj/ncells)
     if(i.eq.inci)then
     inci=inci+ncells
     endif
   enddo
   inci=ncells
   if(j.eq.incj)then
   incj=incj+ncells
   endif
  enddo

  do j = 1,nyc
    do i = 1,nxc
      ilives(i,j)=int(real(nlives)*alpha*noise(i,j,nf))
    enddo
   enddo

  else
   inci=ncells
   incj=ncells 
  do j=1,nyc
   do i=1,nxc
     condition(i,j)=cape(inci/ncells,incj/ncells)
     if(i.eq.inci)then
     inci=inci+ncells
     endif
   enddo
   inci=ncells
   if(j.eq.incj)then
   incj=incj+ncells
   endif
  enddo

  do j = 1,nyc
    do i = 1,nxc
      ilives(i,j)=int(real(nlives)*alpha*noise(i,j,nf))
    enddo
   enddo


 endif !nf

  
!Vertical velocity has its own variable in order to condition on combination
!of "condition" and vertical velocity.

  inci=ncells
  incj=ncells
  do j=1,nyc
   do i=1,nxc
     vertvelhigh(i,j)=vertvelmean(inci/ncells,incj/ncells)
     if(i.eq.inci)then
     inci=inci+ncells
     endif
   enddo
   inci=ncells
   if(j.eq.incj)then
   incj=incj+ncells
   endif
  enddo


  else !ca_global

   do j = 1,nyc
    do i = 1,nxc
     ilives(i,j)=int(real(nlives)*alpha*noise(i,j,nf))
    enddo
   enddo

 endif !sgs/global

!Calculate neighbours and update the automata                                                                                                                                                            
!If ca-global is used, then nca independent CAs are called and weighted together to create one field; CA                                                                                                                                                                                                                                  

  call update_cells(kstep,nca,nxc,nyc,nxch,nych,nlon,nlat,CA,ca_plumes,iini,ilives, &
                   nlives, ncells, nfracseed, nseed,nthresh, ca_global, &
                   ca_sgs,nspinup, condition, vertvelhigh,nf,nca_plumes)
 
  if(ca_global)then
  CAstore(:,:) = CAstore(:,:) + CA(:,:)
  elseif(ca_sgs)then
   if(nf==1)then
   CA_DEEP(:,:)=CA(:,:)
   elseif(nf==2)then
   CA_TURB(:,:)=CA(:,:)
   elseif(nf==3)then
   CA_SHAL(:,:)=CA(:,:)
   elseif(nf==4)then
   CA_RAD(:,:)=CA(:,:)
   else
   CA_MICRO(:,:)=CA(:,:)
   endif
  else
  write(*,*)'Either sgs or global needs to be selected'
  endif

 enddo !nf (nca)

 if(ca_global)then
 CAavg = CAstore / real(nca)
 endif

!smooth CA field

if (ca_smooth .and. ca_global) then
field_in=0.

!get halo
do j=1,nlat
 do i=1,nlon
 field_in(i+(j-1)*nlon,1)=CAavg(i,j)
 enddo
enddo

field_out=0.

call atmosphere_scalar_field_halo(field_out,halo,isize,jsize,k_in,field_in)

do j=1,nlat
 do i=1,nlon
    ih=i+halo
    jh=j+halo
    field_smooth(i,j)=(4.0*field_out(ih,jh,1)+2.0*field_out(ih-1,jh,1)+ & 
                       2.0*field_out(ih,jh-1,1)+2.0*field_out(ih+1,jh,1)+&
                       2.0*field_out(ih,jh+1,1)+2.0*field_out(ih-1,jh-1,1)+&
                       2.0*field_out(ih-1,jh+1,1)+2.0*field_out(ih+1,jh+1,1)+&
                       2.0*field_out(ih+1,jh-1,1))/20.
 enddo
enddo

do j=1,nlat
 do i=1,nlon
   CAavg(i,j)=field_smooth(i,j)
 enddo
enddo

endif !smooth

!In case of ca_global give data zero mean and unit standard deviation

!if(ca_global == .true.)then

!CAmean=0.
!psum=0.
!csum=0.
!do j=1,nlat
! do i=1,nlon
!  psum=psum+CAavg(i,j)
!  csum=csum+1
! enddo
!enddo

!call mp_reduce_sum(psum)
!call mp_reduce_sum(csum)

!CAmean=psum/csum

!CAstdv=0.
!sq_diff = 0.
!do j=1,nlat
! do i=1,nlon
!  sq_diff = sq_diff + (CAavg(i,j)-CAmean)**2.0
! enddo
!enddo

!call mp_reduce_sum(sq_diff)

!CAstdv = sqrt( sq_diff / (csum-1.0) )

!do j=1,nlat
! do i=1,nlon
! CAavg(i,j)=(CAavg(i,j)-CAmean)/CAstdv
! enddo
!enddo

!endif

!Set the range for the nca individual ca_sgs patterns:
if(ca_sgs)then

Detmax(1)=maxval(CA_DEEP(:,:))
call mp_reduce_max(Detmax(1))

do j=1,nlat
 do i=1,nlon
  if(CA_DEEP(i,j)>0.)then
  CA_DEEP(i,j)=CA_DEEP(i,j)/Detmax(1) !Now the range goes from 0-1                                                                                        
  endif
 enddo
enddo

CAmean=0.
psum=0.
csum=0.
do j=1,nlat
 do i=1,nlon
  if(CA_DEEP(i,j)>0.)then
  psum=psum+(CA_DEEP(i,j))
  csum=csum+1
  endif
 enddo
enddo

call mp_reduce_sum(psum)
call mp_reduce_sum(csum)

CAmean=psum/csum

do j=1,nlat
 do i=1,nlon
 if(CA_DEEP(i,j)>0.)then
 CA_DEEP(i,j)=(CA_DEEP(i,j)-CAmean) !Can we compute the median? 
 endif
 enddo
enddo

!!!


!This is used for coupling with the Chikira-Sugiyama deep 
!cumulus scheme. 
do j=1,nlat
 do i=1,nlon
 if(ca_plumes(i,j)==0)then
 ca_plumes(i,j)=20
  endif
 enddo
enddo

endif

!Put back into blocks 1D array to be passed to physics
!or diagnostics output
  
  do blk = 1, Atm_block%nblks
  do ix = 1,Atm_block%blksz(blk)
     i = Atm_block%index(blk)%ii(ix) - isc + 1
     j = Atm_block%index(blk)%jj(ix) - jsc + 1
     Diag(blk)%ca_out(ix)=CAavg(i,j)
     Diag(blk)%ca_deep(ix)=CA_DEEP(i,j)
     Diag(blk)%ca_turb(ix)=CA_TURB(i,j)
     Diag(blk)%ca_shal(ix)=CA_SHAL(i,j)
     Diag(blk)%ca_rad(ix)=CA_RAD(i,j)
     Diag(blk)%ca_micro(ix)=CA_MICRO(i,j)
     Coupling(blk)%ca_out(ix)=CAavg(i,j)        !Cellular Automata
     Coupling(blk)%ca_deep(ix)=CA_DEEP(i,j)
     Coupling(blk)%ca_turb(ix)=CA_TURB(i,j)
     Coupling(blk)%ca_shal(ix)=CA_SHAL(i,j)
     Coupling(blk)%ca_rad(ix)=CA_RAD(i,j)
     Coupling(blk)%ca_micro(ix)=CA_MICRO(i,j)
  enddo
  enddo

 deallocate(omega)
 deallocate(pressure)
 deallocate(humidity)
 deallocate(dp)
 deallocate(cape)
 deallocate(rho)
 deallocate(surfp)
 deallocate(vertvelmean)
 deallocate(vertvelsum)
 deallocate(field_in)
 deallocate(field_out)
 deallocate(field_smooth)
 deallocate(iini)
 deallocate(ilives)
 deallocate(condition)
 deallocate(Detfield)
 deallocate(CA)
 deallocate(ca_plumes)
 deallocate(CAstore)
 deallocate(CAavg)
 deallocate(CA_TURB)
 deallocate(CA_DEEP)
 deallocate(CA_MICRO)
 deallocate(CA_SHAL)
 deallocate(CA_RAD)
 deallocate(noise)
 deallocate(noise1D)

end subroutine cellular_automata
