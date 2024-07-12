      SUBROUTINE filit(m, n, k, nplus, nmins, mhi, nhi, ier)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER m, n, k, mhi, nhi, ier
      INTEGER, DIMENSION(mhi,nhi,3) :: nplus, nmins
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mp1, np1
C-----------------------------------------------

c-----------------------------------------------------------------------
c     auxilary routine for PROGRAM GIDKES2
c     Author:     W.I. van Rij
c     modified:   HNM   STOP replaced by error PARAMETER IER and RETURN
c-----------------------------------------------------------------------
      ier = 0
      IF (m < 0) THEN
         m = -m
         n = -n
      ENDIF

      mp1 = m + 1
      IF (mp1 > mhi) THEN
         WRITE (6, '(a,i3)') ' error in FILIT:   m+1 exceeds mhi = ',
     1      mhi
         ier = 1
         RETURN
      ENDIF

      IF (n<0 .and. m.eq.0) n = -n

      IF (n >= 0) THEN
         np1 = n + 1
         IF (np1 > nhi) THEN
            WRITE (6, '(a,i3)') ' error in FILIT:   n+1 exceeds nhi = '
     1         , nhi
            ier = 2
            RETURN
         ENDIF
         nplus(mp1,np1,k) = 1
         RETURN

      ELSE

         IF ((-n) > nhi) THEN
            WRITE (6, '(a,i3)') ' error in FILIT:   -n  exceeds nhi = '
     1         , nhi
            ier = 2
            RETURN
         ENDIF
         nmins(mp1,(-n),k) = 1
         RETURN

      ENDIF

      END SUBROUTINE filit
