        SUBROUTINE get_vertical_parameters_for_merge(aa,
     &    im, ix, lowst_ipe_level, levs, plow, phigh, 
     &    xpk_low, xpk_high, prsl)

! Get all necessray vertical related parameters to merge the IPE 
! back coupling arrays into the WAM related arrays.
!---------------------------------------------------------------
        USE module_IPE_to_WAM, ONLY: dx
        IMPLICIT none

        REAL, DIMENSION(ix,lowst_ipe_level:levs)  :: aa
        REAL, DIMENSION(ix, levs)                 :: prsl
        REAL, DIMENSION(ix)    :: xpk_low, xpk_high
        REAL                   :: alowipe
        INTEGER, DIMENSION(ix) :: plow, phigh
        INTEGER                :: im, ix, lowst_ipe_level, levs
        INTEGER                :: i, k, low_ipe_level
        
        do i = 1, im

! Dtermine the number of the lowest level in IPE coupling arrays.
! Pre-set up the default values in the IPE arrays are 1E-50.
!----------------------------------------------------------------
          low_ipe_level = levs - 2
          do k = lowst_ipe_level, LEVS
            alowipe = max(-1E30, aa(i, k))
            IF(alowipe > -1E29) THEN
              low_ipe_level = k
              GO TO 1002
            END IF
          end do
1002      CONTINUE

! Get the merging lower and higher levels.
!-----------------------------------------
          plow(i)    = low_ipe_level-1 ! pressure level above which merging is done
          xpk_low(i) = alog(prsl(i, plow(i))) ! its log-pressure

! log-pressure below which merging is done, not necessarily a model level,
! note also that log-pressures are < 0
!-------------------------------------------------------------------------
          xpk_high(i) = xpk_low(i) - dx

! find model level at or just below xpk_high (a very short loop)
!---------------------------------------------------------------
          do k = low_ipe_level, levs
             if (prsl(i,k) < exp(xpk_high(i))) then
                phigh(i) = k - 1
                GO TO 1003
             endif
          enddo
1003      CONTINUE

        end do

        END SUBROUTINE get_vertical_parameters_for_merge




      SUBROUTINE idea_merge_ipe_to_wam(array_ipe, array_wam,  
     &    im, ix, levs, lev_low, prsl, plow, phigh, 
     &    xpk_low, xpk_high)

! This subroutine is to take care the merging calcualtions
! for all six IPE back coupling variable arrays to merge into
! the corresponding WAM model arrays.

      USE module_IPE_to_WAM, ONLY: rdx
      IMPLICIT NONE

! September 2017 Weiyu Yang, initiAL coding.
! September 2017 Rashid Akmaev, some changes suggested
!-----------------------------------------------------
      REAL, DIMENSION(ix, lev_low:levs) :: array_ipe
      REAL, DIMENSION(ix, levs)         :: array_wam
      REAL, DIMENSION(ix, levs)         :: prsl
      REAL, DIMENSION(ix)               :: xpk_low, xpk_high
      REAL :: xk

      INTEGER, DIMENSION(ix) :: plow, phigh
      INTEGER :: i, k, lev_low, levs, ix, im

      do i=1,im
        do k = plow(i)+1, phigh(i)
          xk  = alog(prsl(i,k))
          array_wam(i, k) = 
     &      (array_ipe(i,k) * (xpk_low(i) - xk) +
     &       array_wam(i,k) * (xk - xpk_high(i))) * rdx
        enddo
        do k = phigh(i) + 1, levs
          array_wam(i, k) = array_ipe(i, k)
        end do
      enddo

      END SUBROUTINE idea_merge_ipe_to_wam
