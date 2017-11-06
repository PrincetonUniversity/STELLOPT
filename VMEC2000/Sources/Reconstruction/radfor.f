      SUBROUTINE radfor(pfac0)
      USE vmec_main
      USE vmec_params, ONLY: signgs
      USE vforces, ONLY : r12=>armn_o, gsqrt=>azmn_o, gor=>clmn_o
      USE realspace
      USE vsvd
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec) pfac0
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec) :: p05 = 0.05, p5 = 0.5_dp, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js
      REAL(rprec), DIMENSION(ns) :: vpres
      REAL(rprec) :: delpres, pedge, t1
C-----------------------------------------------
!
!       COMPUTE VPRES, NEEDED FOR F00 PRESSURE BALANCE
!
      gor(2:nrzt) = gsqrt(2:nrzt) / r12(2:nrzt)
      DO js = 2, ns
         vpres(js) =signgs*SUM(gor(js:nrzt:ns)*wint(js:nrzt:ns))
      END DO

      pedge = c1p5*pres(ns) - p5*pres(ns1)
      pressum0 = DOT_PRODUCT(wint(ns:nrzt:ns)*zu0(ns:nrzt:ns),
     1   r1(ns:nrzt:ns,0)+r1(ns:nrzt:ns,1))
      pressum0 = signgs*pedge*pressum0
      pressum0 = pressum0 + hs*DOT_PRODUCT(vpres(2:ns),pres(2:ns))

      IF (pressum0 .eq. zero) pressum0 = one

      pfac0 = pfac
      IF (iresidue .ge. 3) RETURN                    !!Axis moved by fsqr in residue
!
!       COMPUTE AVERAGE FORCE BALANCE CONSTRAINT FOR FIXING R(0)
!       (INTEGRAL OF RADIAL FORCE BALANCE,M=0,N=0, OVER S)
!

      IF (1.e6_dp*fsq .le. one) THEN
         delpres = 0.
         delpres = -fsqsum0/pressum0
         t1 = ABS(delpres)
         IF (t1 .gt. p05) delpres = p05*delpres/t1   !!Wait til close
         pfac0 = pfac*(one + delpres)
      ENDIF

      END SUBROUTINE radfor
