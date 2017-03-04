      SUBROUTINE theta_f (icoil, th0, th1, th2)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE modular_coils
      USE Vwire
      USE coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, mu, icoil
      REAL(rprec) :: theta, th0, th1, th2, t0, t1, t2
      REAL(rprec) :: rs, rc, sk, ck
!-----------------------------------------------
!
!     Computes u(s) and u', u'', from Eq. 2 in Strickler, et al
!     u(s) = theta/2*pi + SUM(rhoc*cos(m*theta)+rhos*sin(m*theta))
!     where theta=2*pi*s
!
      i = icoil
      theta = th0
      t0 = theta
      t1 = 1
      t2 = 0
      DO mu = 0, nf_rho
         rs = modular(i)%rhos(mu)
         rc = modular(i)%rhoc(mu)
         sk = SIN(mu*theta)
         ck = COS(mu*theta)
         t0 = t0 + rs*sk + rc*ck
         t1 = t1 + mu*(rs*ck - rc*sk)       !1st derivative 
         t2 = t2 - mu**2*(rs*sk + rc*ck)    !2nd derivative (for curvature)
      END DO
      th0 = t0
      th1 = t1
      th2 = t2

      END SUBROUTINE theta_f
