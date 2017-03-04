      SUBROUTINE grnflx(xsqr, xs, zs, ansp, idim, nobs)
      USE vsvd
      USE mgrid_mod, ONLY: xobser, zobser, xobsqr
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p5 =0.5_dp, four=4.0_dp
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: idim, nobs
      REAL(rprec), DIMENSION(idim) :: xsqr, xs, zs, ansp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i1, i2, j
      REAL(rprec) :: xt, zt, xtsqr, xt4,
     4   zrp, zrm, xvv, fxu, sqrxu, qqp, qqm
C-----------------------------------------------
!
!       EVALUATES "GREEN'S" FUNCTION FOR POLOIDAL FLUX AT INTERIOR
!       POINTS XS,ZS AND OBSERVATION POINT XT,ZT (ANSP)
!       (RECALL THETA INTEGRATION ONLY FROM ZERO TO PI, SO NEED
!        TO REFLECT ZS TO -ZS, AT LEAST IN UP-DOWN SYMMETRIC CASE)
!
      xt = xobser(nobs)
      zt = zobser(nobs)
      xtsqr = p5*xobsqr(nobs)        !1/2 factor from averaging up,down
      xt4 = four*xt
      DO j = 2,idim
        zrp = zt - zs(j)
        zrm = zt + zs(j)
        xvv =(xt + xs(j))**2
        fxu = xs(j)*xt4
        sqrxu = xsqr(j)*xtsqr
        qqp = fxu/(xvv + zrp*zrp)                !k**2 for zplasma > 0
        qqm = fxu/(xvv + zrm*zrm)                !k**2 for zplasma < 0
!
!       WHICH INDEX LIES BELOW ?
!
        i1 = INT(qqp*odqk2) + 1
        i2 = INT(qqm*odqk2) + 1
!
!       LINEAR INTERPOLATION
!
        ansp(j)  =  sqrxu *( ( yf(i1)+(qqp-qsq(i1))*dyf(i1) )
     >   + ( yf(i2)+(qqm-qsq(i2))*dyf(i2) ) )
      ENDDO

      END SUBROUTINE grnflx
