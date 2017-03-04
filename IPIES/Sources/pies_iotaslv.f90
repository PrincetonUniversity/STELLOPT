!-----------------------------------------------------------------------
!     Subroutine:    pies_iotaslv
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/05/2012
!     Description:   This subroutine calculates iota by minimizing the
!                    function:
!                     f=(theta-iota*phi)^2
!                    along a fieldline.
!-----------------------------------------------------------------------
      SUBROUTINE pies_iotaslv
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_fieldlines
      USE pies_background
      USE pies_realspace
      USE pies_runtime
      USE pies_profile
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL :: fix_test
      INTEGER :: ier, ik, itmax, j
      REAL(rprec) :: reliota, absiota, iotal, iotau, f, fp, pi2
      REAL(rprec) :: dtheta(0:k)
      EXTERNAL :: E04BBF, ftheta
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      pi2 = 8 * ATAN(1._rprec)
      ! We need to fix theta so it doesn't reset at 2pi
      dtheta = 1.0
      DO j = 1, nintw
         dtheta(1:k-1) = thetaln(1:k-1,j)-thetaln(1:k-1,j-1)
         DO ik = 1, k-1
            IF (dtheta(ik) < 0)  thetaln(ik,j:nintw) = thetaln(ik,j:nintw) + pi2
         END DO
      END DO
      
      DO ik = 1, k-1
         iota(ik) = SUM(thetaln(ik,1:nintw)/philn(ik,1:nintw))/nintw
         PRINT *,'ik',iota(ik)
      END DO
      stop
      
      ! Now solve for iota
      iota(0:k) = 0.0
      reliota = 0.0  ! Relative Error
      absiota = 0.0  ! Absolute Error
      PRINT *,''
      DO ik = 1, k-1
         iotal   = 0.0  ! Guess for lower bounds
         iotau   = iotamx  !Guess for upper bounds
         itmax   = 100  ! Maximum iterations
         isurf = ik
         ier     = 0
         CALL E04BBF(ftheta, reliota, absiota, iotal, iotau, itmax, iota(ik), f, fp, ier)
         PRINT *,'ik',ik,iota(ik),f,fp
      END DO
      iota(k) = 2.0*iota(k-1) - iota(k-2)
      iota(0) = 2.0*iota(1)-iota(2)
      WRITE(6,'(6X,f5.3,5X,f5.3)',advance='no') MINVAL(iota),MAXVAL(iota)

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE pies_iotaslv
