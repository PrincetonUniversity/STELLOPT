      SUBROUTINE getgreen
      USE vsvd
      USE vparams, ONLY: twopi
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: jnm, i
      REAL(rprec), DIMENSION(10) :: ak, bk
      REAL(rprec), DIMENSION(9) :: ae, be
      REAL(rprec), DIMENSION(jngrn) :: ye, yk, SQRT1u
      REAL(rprec):: dqk2, qk2, eta, alg, sum1, SUMa, SUM2,
     1     SUMb, SUM3, SUMc, SUM4, SUMd
C-----------------------------------------------

      data ak/3.0072519903686507E-04_dp, 3.9684709020989806E-03_dp,
     1   1.0795990490591656E-02_dp, 1.0589953620989356E-02_dp,
     2   7.5193867218083799E-03_dp, 8.9266462945564728E-03_dp,
     3   1.4942029142282098E-02_dp, 3.0885173001899746E-02_dp,
     4   9.6573590301742396E-02_dp, 1.3862943611198872e+0_dp/
      data bk/6.6631752464607272E-05_dp, 1.7216147097986537E-03_dp,
     1   9.2811603829686118E-03_dp, 2.0690240005100891E-02_dp,
     2   2.9503729348688723E-02_dp, 3.7335546682286003E-02_dp,
     3   4.8827155048118076E-02_dp, 7.0312495459546653E-02_dp,
     4   1.2499999999764055e-1_dp, 5.0000000000000000e-1_dp/
      data ae/3.2519201550638976E-04_dp, 4.3025377747931137E-03_dp,
     1   1.1785841008733922E-02_dp, 1.1841925995501268E-02_dp,
     2   9.0355277375409049E-03_dp, 1.1716766944657730E-02_dp,
     3   2.1836131405486903E-02_dp, 5.6805223329308374E-02_dp,
     4   4.4314718058336844E-1_dp/
      data be/7.2031696345715643E-05_dp, 1.8645379184063365E-03_dp,
     1   1.0087958494375104E-02_dp, 2.2660309891604169E-02_dp,
     2   3.2811069172721030E-02_dp, 4.2672510126591678E-02_dp,
     3   5.8592707184265347E-02_dp, 9.3749995116366946E-02_dp,
     4   2.4999999999746159E-1_dp/

!
!       Compute "Green's Functions" for Poloidal Flux, 2*pi*R*A-sub-phi,
!       BR, and BZ at point (XT,ZT) due to unit current (mu0*I = 1) at (XS,ZS) ...
!       modified to interpolate on k**2 - 3-34-92 - sph
!
      jnm = jngrn - 1
      odqk2 = (jnm)
      dqk2 = 1.0_dp/odqk2
      DO i = 2, jnm
         qk2 = dqk2*(i - 1)
         qsq(i) = qk2
         eta = 1 - qk2
         alg = log(eta)
         sum1 = ((((ak(1)*eta+ak(2))*eta+ak(3))*eta+ak(4))*eta+ak(5))*
     1      eta + ak(6)
         SUMa = (((sum1*eta + ak(7))*eta+ak(8))*eta+ak(9))*eta + ak(10)
         SUM2 = ((((bk(1)*eta+bk(2))*eta+bk(3))*eta+bk(4))*eta+bk(5))*
     1      eta + bk(6)
         SUMb = (((sum2*eta + bk(7))*eta+bk(8))*eta+bk(9))*eta + bk(10)
         yk(i) = SUMa - alg*sumb
         SUM3 = (((ae(1)*eta+ae(2))*eta+ae(3))*eta+ae(4))*eta
         SUMc = (((((sum3 + ae(5))*eta+ae(6))*eta+ae(7))*eta+ae(8))*eta+
     1      ae(9))*eta
         SUM4 = (((be(1)*eta+be(2))*eta+be(3))*eta+be(4))*eta
         SUMd = (((((sum4 + be(5))*eta+be(6))*eta+be(7))*eta+be(8))*eta+
     1      be(9))*eta
         ye(i) = SUMc - alg*sumd + 1
         yf(i) = ((1 + eta)*yk(i)-2*ye(i))/qk2
      END DO
      ye(1) = 0.25_dp*twopi
      ye(jngrn) = 1
      yk(1) = ye(1)
      yk(jngrn) = 2*yk(jnm) - yk(jngrn-2)
      yf(1) = 0.
      yf(jngrn) = 2*yf(jnm) - yf(jngrn-2)
      qsq(1) = 0
      qsq(jngrn) = 1

      SQRT1u = SQRT(qsq(:jngrn))/twopi
c                                      !Factor of 1/2 from SQRT(4*xs*xt)
      yek(:jngrn) = 0.5_dp*sqrt1u*(ye(:jngrn)-yk(:jngrn))
      yeq(:jngrn) = 0.25_dp*qsq(:jngrn)*sqrt1u*ye(:jngrn)
c                                 !Factor of 2 absorbed by SQRT(4 xt xs)
      yf(:jngrn) = twopi*sqrt1u*yf(:jngrn)
      dyek(:jnm) = (yek(2:jnm+1)-yek(:jnm))*odqk2
      dyeq(:jnm) = (yeq(2:jnm+1)-yeq(:jnm))*odqk2
      dyf(:jnm) = (yf(2:jnm+1)-yf(:jnm))*odqk2

      END SUBROUTINE getgreen
