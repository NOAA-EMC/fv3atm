!-----------------------------------------------------------------------
!
      module module_kinds
!
!-----------------------------------------------------------------------
!
      integer,parameter :: isingle=selected_int_kind(r=9)
      integer,parameter :: idouble=selected_int_kind(r=18)
      integer,parameter :: single=selected_real_kind(p=6,r=37)
      integer,parameter :: double=selected_real_kind(p=13,r=200)

      integer,parameter:: &
       klog=4 &
      ,kint=isingle &
      ,kdin=idouble &
      ,kfpt=double &
      ,kdbl=double
!
!-----------------------------------------------------------------------
!
      end module module_kinds
!
!-----------------------------------------------------------------------
