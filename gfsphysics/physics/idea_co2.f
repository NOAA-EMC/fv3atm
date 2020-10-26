      subroutine idea_co2(im,ix,levs,nlev,ntrac,cp,adr,adt,        
     &dtdt,cosz,dtdth)
!SK   subroutine idea_co2(im,ix,levs,nlev,ntrac,grav,cp,adr,adt,        
!SK  &dtdt,cosz,dtdth)
!
! Apr 06 2012   Henry Juang, initial implement for nems
! Dec 13 2012   Jun Wang     move init step out of column physics
! Feb 13 2012   Jun Wang     move gravity array gg to idea_compistion module
! Sep 14 2020   Sajal Kar    remove unused grav array
!
      use co2pro_mod, only: co2my
!     use co2c_mod
!     use qnir_mod
      use physcons, only : amo2=>con_amo2, amo3=>con_amo3,                    
     &                     amh2o=>con_amw
      use idea_composition, only : amo, amn2, prlog, k43, mpi_me
!     use module_physics_driver, only: mpi_me
      implicit none
! Argument
      integer, intent(in) :: im  ! number of data points in adt (first dim)
      integer, intent(in) :: ix  ! max data points in adt (first dim)
      integer, intent(in) :: levs   ! number of pressure levels
      integer, intent(in) :: nlev   ! number of pressure levels in calculation
      integer, intent(in) :: ntrac  ! number of tracer
      real, intent(in)    :: adr(ix,levs,ntrac) ! tracer
      real, intent(in)    :: adt(ix,levs)    ! temperature
      real, intent(in)    :: cp(ix,levs)    ! J/kg/k
!SK   real, intent(in)    :: grav(ix,levs)    ! g (m/s2)
      real, intent(in)    :: cosz(im)    !cos solar zenith angle 
!hmhj character*(*), intent(in) ::   dir    ! directory located coef files
      real, intent(out)   :: dtdt(ix,levs)    ! cooling rate k/s
      real, intent(out)   :: dtdth(ix,levs)   ! heating rate k/s
!
      real pmod(levs),q_n2(ix,nlev),ma(ix,nlev)                         
     &,q_o(ix,nlev),q_o2(ix,nlev),hold(levs)
      integer i,k,kk
      integer, parameter :: skprnt = .false.
!
! precalling
      dtdth(:,:)=0.
      dtdt(:,:) =0.
!
!SK2020Sep20
      if (skprnt) then
       if (minval(adt).lt.0.) then
       write(2000+mpi_me,98) minval(adt),maxval(adt),mpi_me
98     format(1x,'min adt=',e15.7,1x,'max adt=',e15.7,' me=',i4)
       endif
       if (minval(adr).lt.0.) then
       write(2000+mpi_me,99) minval(adr),maxval(adr),mpi_me
99     format(1x,'min adr=',e15.7,1x,'max adr=',e15.7,' me=',i4)
       endif
      endif
!
      do i=1,im
        do k=k43,levs
          kk=k-k43+1
          q_n2(i,kk)=1.-adr(i,k,4)-adr(i,k,5)-adr(i,k,1)-adr(i,k,2)
          ma(i,kk)=1./(adr(i,k,4)/amo+adr(i,k,5)/amo2+adr(i,k,1)/amh2o+ 
     &             adr(i,k,2)/amo3+q_n2(i,kk)/amn2)
          q_o(i,kk)=adr(i,k,4)*ma(i,kk)/amo
          q_o2(i,kk)=adr(i,k,5)*ma(i,kk)/amo2
          q_n2(i,kk)=q_n2(i,kk)*ma(i,kk)/amn2
        enddo
      enddo
!     print*,'www2',im,ix,q_o(1:im1,nlev)

      if (skprnt) then
       if (minval(q_n2).lt.0. or. minval(q_o2).lt.0. or.
     &     minval(q_o).lt.0. or. minval(ma).lt.0.) then
       write(1000+mpi_me,100) minval(q_n2),maxval(q_n2),mpi_me
100    format(1x,'min q_n2=',e15.7,1x,'max q_n2=',e15.7,' me=',i4)
       write(1000+mpi_me,101) minval(q_o2),maxval(q_o2),mpi_me
101    format(1x,'min q_o2=',e15.7,1x,'max q_o2=',e15.7,' me=',i4)
       write(1000+mpi_me,102) minval(q_o),maxval(q_o),mpi_me
102    format(1x,'min q_o=',e15.7,1x,'max q_o=',e15.7,' me=',i4)
       write(1000+mpi_me,103) minval(ma),maxval(ma),mpi_me
103    format(1x,'min ma=',e15.7,1x,'max ma=',e15.7,' me=',i4)
       endif
      endif
!     return
!
! CO2 cooling
      call co2cc(ix,im,prlog,adt,levs,prlog(k43),                       
     &           dtdt(1,k43),nlev,ma,q_o,q_o2,q_n2)
!sk2020sep14
!     print*,'sk2020sep14:idea_co2:min(dtdt(1:im,k43:levs))',
!    & minval(dtdt(1:im,k43:levs))
!     print*,'sk2020sep14:idea_co2:max(dtdt(1:im,k43:levs))',
!    & maxval(dtdt(1:im,k43:levs))
!     print*,'sk2020sep14/idea_co2<-co2cc'
!     RETURN
!sk

! J/kg/s to k/s
      do i=1,im
        do k=k43,levs
          dtdt(i,k)=dtdt(i,k)/cp(i,k)
        enddo
!vay-16          dtdt(i,1:k43-1)=0.
      enddo
! CO2 heating
      do i=1,im
        call qnirc(cosz(i),prlog(k43),co2my,hold(k43),nlev)
        do k=k43,levs
!       dtdth(i,k)=hold(k-k43+1)
        dtdth(i,k)=hold(k)
        enddo
!vay-16        dtdth(i,1:k43-1)=0.
      enddo
!sk2020sep14
!     print*,'sk2020sep14:idea_co2:min(dtdth(1:im,k43:levs))',
!    & minval(dtdth(1:im,k43:levs))
!     print*,'sk2020sep14:idea_co2:max(dtdth(1:im,k43:levs))',
!    & maxval(dtdth(1:im,k43:levs))
!     print*,'sk2020sep14/idea_co2<-qnirc'
!     RETURN
!sk
      return
      end
