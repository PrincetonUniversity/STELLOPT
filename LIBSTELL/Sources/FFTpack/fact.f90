      SUBROUTINE fact_g(n, ifax)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER n
      INTEGER, DIMENSION(13) :: ifax
#if !defined(CRAY) || defined(LONESTAR) || defined(MCURIE)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nn, k, l, MAX, inc
!-----------------------------------------------
!     factorization routine that first extracts ALL factors of 4
      IF (n.gt.1) GOTO 10
      ifax(1) = 0
      IF (n.lt.1) ifax(1) = -99
      RETURN
   10 nn=n
      k=1
!     test for factors of 4
   20 IF (MOD(nn,4).ne.0) GOTO 30
      k=k+1
      ifax(k)=4
      nn=nn/4
      IF (nn.eq.1) GOTO 80
      GOTO 20
!     test for extra factor of 2
   30 IF (MOD(nn,2).ne.0) GOTO 40
      k=k+1
      ifax(k)=2
      nn=nn/2
      IF (nn.eq.1) GOTO 80
!     test for factors of 3
   40 IF (MOD(nn,3).ne.0) GOTO 50
      k=k+1
      ifax(k)=3
      nn=nn/3
      IF (nn.eq.1) GOTO 80
      GOTO 40
!     now find remaining factors
   50 l=5
      MAX = SQRT(REAL(nn,rprec))
      inc=2
!     inc alternately takes on values 2 and 4
   60 IF (MOD(nn,l).ne.0) GOTO 70
      k=k+1
      ifax(k)=l
      nn=nn/l
      IF (nn.eq.1) GOTO 80
      GOTO 60
   70 IF (l.gt.max) GOTO 75
      l=l+inc
      inc=6-inc
      GOTO 60
   75 k = k+1
      ifax(k) = nn
   80 ifax(1)=k-1
!     ifax(1) now CONTAINS number of factors
#else
      CALL fact (n, ifax)
#endif
      END SUBROUTINE fact_g
