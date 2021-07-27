      MODULE module_IPE_to_WAM
!
! This module povides the parameters, arrays and control flags
! for the interface of IPE coupling to WAM.

! Weiyu Yang -- August 26, 2017.
!-------------------------------------------------------------

!     lowst_ipe_level is the number of the most low level over global
!     IPE variable fields from the IPE coupling. 
!--------------------------------------------------------------------
      INTEGER :: lowst_ipe_level
      LOGICAL :: ipe_to_wam_coupling

! dx is the vertical extent of the merge layer in log-pressure 
! (scale heights), may be increased for smoother merging.
!-------------------------------------------------------------
      REAL, parameter :: dx=1.0, rdx=1.0/dx

      REAL, DIMENSION(:, :, :), POINTER :: ZMT, MMT, JHR, SHR, O2DR
      
      END MODULE module_IPE_to_WAM
