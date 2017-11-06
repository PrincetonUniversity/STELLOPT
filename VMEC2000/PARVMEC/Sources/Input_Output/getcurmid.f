      SUBROUTINE getcurmid (curmid, izeta, gsqrt, r12)
      USE vmec_input, ONLY: rprec, dp, nzeta
      USE vmec_dim, ONLY: ns, ns1, ntheta2
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec) :: curmid(2*ns)
      REAL(rprec) :: izeta(ns,nzeta,*), gsqrt(ns,nzeta,*),
     1    r12(ns,nzeta,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: midcur(ns)
C-----------------------------------------------
!     THETA = pi, PHI = 0
      midcur(2:ns) = r12(2:ns,1,ntheta2)/gsqrt(2:ns,1,ntheta2)

      curmid(1) = izeta(ns,1,ntheta2)*midcur(ns)
      curmid(2:ns1) = 0.5_dp*izeta(ns1:2:-1,1,ntheta2)*
     1                   (midcur(ns1:2:-1) + midcur(ns:3:-1))

!     THETA = 0, PHI = 0
      midcur(2:ns) = r12(2:ns,1,1)/gsqrt(2:ns,1,1)

      curmid(ns+1:2*ns-1) = 0.5_dp*izeta(2:ns1,1,1)*
     1                   (midcur(2:ns1) + midcur(3:ns))

      curmid(ns) = 0.5_dp*(curmid(ns-1) + curmid(ns+1))
      curmid(2*ns) = 2*curmid(2*ns-1) - curmid(2*ns-2)

      END SUBROUTINE getcurmid
