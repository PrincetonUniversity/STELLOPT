      FUNCTION vnorm (n, v, w)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER n
      REAL(rprec), DIMENSION(n) :: v, w
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: sum0, vnorm
C-----------------------------------------------
clll. optimize
c-----------------------------------------------------------------------
c this function routine computes the weighted root-mean-square norm
c of the vector of length n contained in the array v, with weights
c contained in the array w of length n..
c   vnorm = SQRT( (1/n) * SUM( v(i)*w(i) )**2 )
c-----------------------------------------------------------------------
      sum0 = SUM((v*w)**2)
      vnorm = SQRT(sum0/(n))

      END FUNCTION vnorm
