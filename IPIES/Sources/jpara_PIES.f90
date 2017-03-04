!-----------------------------------------------------------------------
!     Subroutine:    jpara_PIES
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/13/2012
!     Description:   This subroutine calculates the parallel component
!                    of the total current density in the old PIES way.
!-----------------------------------------------------------------------
      SUBROUTINE jpara
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_background
      USE pies_realspace
      USE pies_runtime
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier,u,v,mn
      REAL(rprec) :: a11, a12, a13, a23, a31, a32, a33
      REAL(rprec), ALLOCATABLE :: dmn1(:)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      DO mn = 1, mnmax
         dmn1 = jumnc(mn,:) * ( xm(mn)*rmnc(mn,:)*xn(mn)*rmnc(mn,:) + xm(mn)*zmns(mn,:)*xn(mn)*zmns(mn,:) )
         dmn2 = jvmnc(mn,:) * ( xn(mn)*rmnc(mn,:)*xn(mn)*rmnc(mn,:) + xn(mn)*zmns(mn,:)*xn(mn)*zmns(mn,:) )
         dmn3 = jumnc(mn,:) * ( xm(mn)*rmnc(mn,:)*xm(mn)*rmnc(mn,:) + xm(mn)*zmns(mn,:)*xm(mn)*zmns(mn,:) )
         dmn4 = jvmnc(mn,:) * ( xm(mn)*rmnc(mn,:)*xn(mn)*rmnc(mn,:) + xm(mn)*zmns(mn,:)*xn(mn)*zmns(mn,:) )
         mumnc(mn,:) = ( dmn1 + dmn2 + iota(:) * (dmn3 + dmn4) ) / modb(mn,:)
      END DO

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE jpara
