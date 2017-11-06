      SUBROUTINE getpresprofile
      USE vmec_main
      USE vsvd
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i
      REAL(rprec) :: sigmin
C-----------------------------------------------

!
!       WRITE OVER THOMPSON DATA
!
      lpofr = .false.                   !!these data are at s-half nodes
      lpprof = .false.
      itse = 10
      pthommax = datathom(1)                  !!compute in final version
      sigMIN = 0.03*pthommax

      DO i = 1, itse
         rthom(i) = REAL(i - 1,rprec)/(itse - 1)
         datathom(i) = (pthommax - sigmin)*(1. - rthom(i))**2 + sigmin
         sigma_thom(i) = 0.2*pthommax
      END DO

      END SUBROUTINE getpresprofile
