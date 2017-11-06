      SUBROUTINE grnbfld(xsqr, xs, zs, br, bz, idim, nobs1, nobs2)
      USE vsvd
      USE vparams, ONLY: one
      USE mgrid_mod, ONLY: rbcoil, zbcoil, rbcoilsqr
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: four=4.0_dp, p5=0.5_dp
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER idim, nobs1, nobs2
      REAL(rprec), DIMENSION(idim) :: xsqr, xs, zs, br, bz
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER j, i1, i2
      REAL(rprec) :: xt, zt, xtsqr, xt4, oxt, zrp, zrm, xvv,
     1 fxu, sqrxu, qqp, qqm, delqp, delqm, yeqp, yeqm, brp, brm
C-----------------------------------------------
!
!       COMPUTE BR = (ZT-ZS)/RT/SQRT(4*RT*RS) * F1(k)
!               BZ = 1/SQRT(4*RT*RS)*[RS/RT * F2(k) - F1(k)]
!       WHERE   F1 = k/2pi[ (E(k) - K(k)) + q1(k)*E(k) ]
!               F2 = k/(2pi)  [ q1(k)*E(k) ]
!               q1 = .5*k**2/(1. - k**2)   [Most SINgular piece near k=1]
!               k**2 = 4*RT*RS/[(RT+RS)**2 + (ZT-ZS)**2]
!
      xt = rbcoil(nobs1,nobs2)
      zt = zbcoil(nobs1,nobs2)
      xtsqr = p5/rbcoilsqr(nobs1,nobs2)            !1/2 from symmetrizing
      xt4 = four*xt
      oxt = one/xt

      DO j = 2,idim
        zrp = zt - zs(j)
        zrm = zt + zs(j)
        xvv =(xt + xs(j))**2
        fxu = xs(j)*xt4
        sqrxu = xtsqr/xsqr(j)
        qqp = fxu/(xvv + zrp*zrp)
        qqm = fxu/(xvv + zrm*zrm)
!
!       WHICH INDEX LIES BELOW ?
!
        i1 = INT(qqp*odqk2) + 1
        i2 = INT(qqm*odqk2) + 1
!
!       LINEAR INTERPOLATION
!
        delqp = qqp - qsq(i1)
        delqm = qqm - qsq(i2)
        yeqp = (yeq(i1) + delqp*dyeq(i1))/(one - qqp)
        yeqm = (yeq(i2) + delqm*dyeq(i2))/(one - qqm)
        brp = yek(i1) + delqp*dyek(i1)
        brm = yek(i2) + delqm*dyek(i2)
        br(j) = sqrxu*oxt*(zrp*(brp+yeqp) + zrm*(brm+yeqm))
        bz(j) = sqrxu*((xs(j)*oxt-1.0)*(yeqp + yeqm) - (brp + brm))
      ENDDO

      END SUBROUTINE grnbfld
