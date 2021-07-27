      MODULE IDEA_TRACERS_INPUT
!
!VAY-2015 initial code
!
      use IDEA_IO_UNITS, only : iulog
      use physcons, only : avgd => con_avgd,             
     &                     amo3 => con_amo3,amh2o => con_amw
!
      use idea_composition, only : bz, amo, amo2, amn2,mmr_min,mmr_max
!
       implicit none
       real ::  rbz1000, rkgavgd, rbz
       real ::  rmo, rmo2, rmn2, rmo3, rmh2o
       real ::   mo, mo2, mn2, mo3, mh2o    
       public :: hprofile
       public :: Jprofile
       public :: WAM_GLOBAL_TRACERS
       public :: init_tracer_constants
       contains
!
       subroutine  init_tracer_constants
       implicit none
!
! TO DO: should be moved to "tracer_init"
!
       rbz= 1./bz
       rbz1000= 1000./bz
       rkgavgd =1.e-3/avgd
       mo=amo *rkgavgd
       mo2=amo2 *rkgavgd
       mn2=amn2 *rkgavgd
       mh2o=amh2o *rkgavgd
       mo3=amo3 *rkgavgd

       rmn2  = 1./mn2
       rmo   = 1./mo
       rmo3  = 1./mo3
       rmo2  = 1./mo2
       rmh2o = 1./mh2o
       end subroutine  init_tracer_constants
!
       SUBROUTINE     WAM_GLOBAL_TRACERS(levs, nvmr, vmr_glob)
       use netcdf      
       use idea_composition,  only: mpi_me, mpi_master
!SK    use idea_mpi_def, ONLY:  info, mpi_comm_all
     
       implicit none
!       include 'mpif.h'
!
! Read in WAM_GLOBAL_TRACERS from nc-file
! /scratch3/NCEPDEV/swpc/save/Valery.Yudin/BASE_SVN/BASE_WAM_DATA/WAM_COMP
!
       INTEGER                ::  levs, nvmr
       REAL, intent(out)      ::  vmr_glob(levs, nvmr)
!
       Character(len=128) ::File_glob = 'wam_rm3_globcomp.nc'
!
!       Character(len=128) ::File_glob = '/wam_vmr_globcomp.nc'
! locals
!
        integer ::  ierr
        integer ::  ncid, vid, ierNC
        integer  :: astat, dimidt(2)  
        integer  :: nvars, nlevs
!
! 
       ierNC=NF90_OPEN(trim(File_glob), nf90_nowrite, ncid)   
       if (iernc /=0) then 
           print *,  ncid, 'ncid ', iernc, ' iernc '
           print *, ' VAY-File_glob ', File_glob 
       endif

         iernc=nf90_inq_varid( ncid, 'vmr_wam', vid )
         iernc=nf90_inquire_variable(ncid, vid, dimids=dimidT) 
         iernc = nf90_inquire_dimension(ncid, dimidT(2), len=nvars)
         iernc = nf90_inquire_dimension(ncid, dimidT(1), len=nlevs)

        if (mpi_me.eq.mpi_master) then
         if (nvars.ne.nvmr.or.nlevs.ne.levs) then
           print *,  nvmr, ' nvmr ',  nvars, ' nvars '
           print *,  levs, ' levs ',  nlevs, ' nlevs '
           print *, ' VAY-incorrect dimensions in WAM_GLOBAL_TRACERS'
         endif
        endif
        iernc = nf90_inq_varid( ncid, 'vmr_wam', vid )
        iernc = nf90_get_var( ncid, vid, vmr_glob)
        iernc=nf90_close(ncid) 
        RETURN
!==========================================================   
!         1          5  6  7  8
!"names: h2o o3 n2 o o2 h oh ho2 no n co2 ndens tg mu zkm"
!
!==========================================================
       print *, ' ozone global ', maxval( vmr_glob(1:levs,2)),
     & minval( vmr_glob(1:levs,2)) 
       print *, ' o3p global ', maxval( vmr_glob(1:levs,4)),
     & minval( vmr_glob(1:levs,4)) 
!  
!
       END SUBROUTINE WAM_GLOBAL_TRACERS

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!   
!-------------------------------------------------------------------------
      SUBROUTINE jprofile(levs,f107, J)
! get photo dissociation rate
! vay 2015
!      use idea_solar_input,  only : f107 => wf107_s
      implicit none
   
      integer, parameter :: np=17  !number of pressure levels of orig
      integer, intent(in) :: levs  !number of pressure levels of output 
      real,    intent(in) :: f107  !  as an input due to SWPC-SAIR-WAM drivers
      real,    intent(out):: J(levs)
! local variables
      real JI(np),FHT(np),C(np),J17(np)
      integer k
!
      DATA C/8*0.900,0.680,0.43,0.18,6*-0.066/
      DATA JI/.4e-8,.78e-8,1.5e-8,3.e-8,6.8e-8,.15e-6,.34e-6,.77e-6,    
     &1.07e-6,1.35e-6,1.6e-6,1.81e-6,2.05e-6,2.23e-6,2.36e-6,2.5e-6,    
     &2.57e-6/
      DATA FHT/8*1.2,1.85,2.50,3.150,6*3.8/
! calculate photo dissociation rate (/s) in Tims 17 pressure grid
      do k=1,17                                                   
        J17(k)=JI(k)*((FHT(k)-1.0)*f107/176.+C(k))
      enddo
! interplate to GFS pressure grid
      call z17toz(levs,J17,J,0.)
      return
      end subroutine jprofile
!-------------------------------------------------------------------------
      subroutine z17toz(levs,ain,aout,down)
! interpolate 17 pressure levels (from Tim's grid) to
! idea pressure grid pr(levs)
      use idea_composition, only : pr=> pr_idea
      implicit none
      integer, parameter :: np=17  !number of pressure levels of input
      integer, intent(in) :: levs  !number of pressure levels of output 
      real,    intent(in) :: ain(np)  !input field in 17 pressure grid
      real,    intent(in) :: down     !field value under 5.2285Pa
      real,    intent(out):: aout(levs)!output in levs pressure grid
!local variable
      real p17(np),z17(np),z(levs),dz
      integer kref,k,i
!
      do k=1,np
        p17(k)=5.2285*exp(1.-k)
        z17(k)=-1.*log(p17(k))
      enddo
      do k=1,levs
        z(k)=-1.*log(pr(k)*100.)
      enddo
      do k=1,levs
        kref=0
        do i=1,np-1
          if(z(k).ge.z17(i).and.z(k).le.z17(i+1)) then
            kref=i
            dz=(z(k)-z17(i))/(z17(i+1)-z17(i))
          endif
        enddo
        if(kref.ne.0) then
          aout(k)=dz*ain(kref+1)+(1.-dz)*ain(kref)
        elseif(z(k).lt.z17(1)) then
          aout(k)=down
        elseif(z(k).gt.z17(17)) then
          aout(k)=ain(17)
        endif
      enddo
      return
      end subroutine z17toz
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc! 
!-------------------------------------------------------------------------
      SUBROUTINE hprofile(levs,oh,ho2)
      implicit none
      integer, parameter :: np=65  !number of pressure levels of orig
      integer, intent(in) :: levs  !number of pressure levels of output 
      real,    intent(out):: oh(levs),ho2(levs)
! local variables
      real ohi(np),ho2i(np)               ! units of conc-n 1/cm3 or 1/m3 
      integer k
!
      data ohi/9.5e12,1.3e13,1.6e13,1.8e13,2.0e13,
     &   2.1e13,2.2e13,1.9e13,1.3e13,6.0e12,2.1e12,
     &   1.0e12,3.0e11,1.0e11,4.0e10,1.6e10,7.0e9,
     &   3.2e9,1.2e9,5.0e8,2.0e8,44*0.0/
      data ho2i/7.0e12,9.0e12,1.2e13,1.3e13,1.6e13,
     &   1.7e13,1.7e13,1.3e13,7.0e12,2.5e12,7.0e11,
     &   1.5e11,3.5e10,8.0e9,2.5e9,9.0e8,3.0e8,
     &   1.4e8,6.0e7,2.5e7,1.0e7,44*0.0/
! interplate to GFS pressure grid
      call z65toz(levs,ohi,oh,0.)
      call z65toz(levs,ho2i,ho2,0.)
      return
      end SUBROUTINE hprofile
!
      subroutine z65toz(levs,ain,aout,down)
! interpolate 65 pressure levels (from Tim's grid) to
! idea pressure grid pr(levs)
      use idea_composition, only : pr=> pr_idea
      implicit none
      integer, parameter  :: np=65  !number of pressure levels of input
      integer, intent(in) :: levs  !number of pressure levels of output 
      real,    intent(in) :: ain(np)  !input field in 65 pressure grid
      real,    intent(in) :: down     !field value under 5.2285Pa
      real,    intent(out):: aout(levs)!output in levs pressure grid
!local variable
      real p65(np),z65(np),z(levs),dz
      integer kref,k,i
!
      do k=1,np
        p65(k)=5.2285*exp((1.-k)*.25)
        z65(k)=-log(p65(k))
      enddo
      do k=1,levs
        z(k)=-log(pr(k)*100.)
      enddo
!
      do k=1,levs
        kref=0
        do i=1,np-1
          if(z(k).ge.z65(i).and.z(k).le.z65(i+1)) then
            kref=i
            dz=(z(k)-z65(i))/(z65(i+1)-z65(i))
          endif
        enddo
        if(kref.ne.0) then
          aout(k)=dz*ain(kref+1)+(1.-dz)*ain(kref)
        elseif(z(k).lt.z65(1)) then
          aout(k)=down
        elseif(z(k).gt.z65(65)) then
          aout(k)=ain(np)
        endif
      enddo
      return
      end subroutine z65toz
!
      END MODULE IDEA_TRACERS_INPUT
