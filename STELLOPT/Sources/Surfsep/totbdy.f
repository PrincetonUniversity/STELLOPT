      SUBROUTINE totbdy(ntheta, phiin, r, z,
     1   rmnc, zmns, xm, xn, mnmax)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ntheta, mnmax
      REAL(rprec), DIMENSION(*), INTENT(out) :: r, z
      REAL(rprec), DIMENSION(*), INTENT(in) ::
     1   rmnc, zmns, xm, xn
      REAL(rprec) :: phiin

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nrth, mn, kt, m, iarg
      REAL(rprec) :: pit, xm0, xn0, arg, twopi
C-----------------------------------------------
      twopi = 8*ATAN(1._dp)
      pit = 1._dp/(ntheta - 1)
      nrth = 1*ntheta
      r(:nrth) = 0
      z(:nrth) = 0
         DO mn = 1, mnmax
            m = mn
            IF(rmnc(m).eq.0._dp)cycle
            xm0 = xm(mn)*pit
            xn0 = xn(mn)
            DO kt = 1, ntheta
               arg = xm0*(kt - 1) - xn0*phiin/twopi
               iarg = arg
               arg = twopi*(arg - iarg)
               r(kt) = r(kt) + rmnc(m)*COS(arg)
               z(kt) = z(kt) + zmns(m)*SIN(arg)
            END DO
         END DO
      END SUBROUTINE totbdy
