      SUBROUTINE fixaray(pexp, qexp, mexp)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: mexp
      REAL(rprec), INTENT(in) :: pexp, qexp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l, nn0, n, m
      REAL(rprec) :: sgn
C-----------------------------------------------
!
!       This routine stores toroidal and poloidal mode number arrays
!
!       MPOL = NUMBER OF POLOIDAL MODES USED IN CURVE-FIT
!       NPHI = NUMBER OF EQUALLY-SPACED PHI PLANES PER FIELD
!              PERIOD IN WHICH THE CURVE-FITTING IS PERFORMED
!       IMPORTANT: NPHI MUST BE EVEN FOR NON-SYMMETRIC SYSTEMS
!
!       MPNT = NUMBER OF R AND Z MODES APPEARING IN FIT
!       NTHETA = NUMBER OF THETA POINTS IN EACH PHI PLANE
!       N2   = TOTAL NUMBER OF RESIDUALS IN FSQ (PER PHI PLANE)
!            == 2 (m=0 modes) + 2*mrho (rhoc, rhos) + ntheta
!       NFP  = NUMBER OF TOROIDAL FIELD PERIODS (IRRELEVANT)
!
      mrho = mpol-1
      nphi2 = nphi/2
      IF (nphi2 .eq. 0) nphi2 = 1
      n2 = 2*(mrho+1) + ntheta

      l = 0
      ntor = MAX0(1,nphi-1)

      IF (nphi .eq. 2) ntor = 2

      nn0 = 1 - (ntor + 1)/2
      DO n = 1, ntor
         nn(n) = (nn0 + (n - 1))*nfp
      END DO
      DO m = 1, mpol
         mm(m) = m - 1
         DO n = 1, ntor
            IF (mm(m).ne.0 .or. nn(n).ge.0) THEN
               l = l + 1
               m1(l) = mm(m)
               n1(l) = nn(n)
            ENDIF
         END DO
      END DO
      mpnt = l
      dnorm = 2._dp/ntheta

      DO m = 0, mpol
         dm1(m) = m
         IF (m .eq. 0) THEN
            xmpq(m,1) = zero
            xmpq(m,2) = zero
         ELSE
            xmpq(m,1) = dm1(m)**(pexp+qexp)
            xmpq(m,2) = dm1(m)**(pexp)
         ENDIF
      END DO

!
!     FOR SECOND DERIVATIVE-TYPE CONSTRAINT, CHOOSE SGN = -one
!
      sgn = one                            !First derivative-TYPE constraint
!     sgn = -one                      !Second derivative-TYPE constraint
      t1m(1) = one
      DO m = 2, mrho
        t1m(m) = (REAL(m-1, rprec)/m)**mexp
      END DO
      DO m = 1,mrho-2
        t2m(m) = sgn*(REAL(m+1, rprec)/m)**mexp
      END DO
      t2m(mrho-1) = zero
      t2m(mrho) = zero

      END SUBROUTINE fixaray
