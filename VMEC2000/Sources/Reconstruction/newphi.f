      SUBROUTINE newphi(phipog)
      USE stel_kinds
      USE vmec_main
      USE realspace
      USE vsvd
      USE xstuff
c-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(nrzt) :: phipog
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: phifold, tnorm
C-----------------------------------------------

!
!     UPDATE PHI SCALE FACTOR (PHIFAC)
!
      phifold = phifac
      phifac = xc(neqs1)

      IF (phifold .eq. zero) STOP 'phifac = 0 in newphi'
      tnorm = phifac/phifold
      phip = tnorm * phip
      phipog = tnorm * phipog

      END SUBROUTINE newphi
