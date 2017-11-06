      SUBROUTINE splinint(grn, cm, jacob, h, u, u1, w, w1, nk, nots,
     1   ifunc, nmesh)
      USE vparams
      USE vsvd
      USE vspline
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER nots, ifunc, nmesh
      INTEGER, DIMENSION(*) :: nk
      REAL(rprec), DIMENSION(nmesh) :: grn, cm
      REAL(rprec), DIMENSION(nots) :: jacob, h
      REAL(rprec), DIMENSION(nmesh) :: u, u1, w, w1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: j, k, nb, ia, ib, k1, ksp1, nmesh1
      REAL(rprec), DIMENSION(nmesh) :: func
      REAL(rprec), DIMENSION(nots) :: af, bs
C-----------------------------------------------

!
!       COMPUTES af,bs FACTORS IN
!       (ifunc=INTFUN)     Int[ GRN(X) * F(X) ds ]      or
!       (ifunc=INTDER)     Int[ GRN(X) * {d(CmF)/ds} ds ]
!                 =   af(k)*f(k) + bs(k)*f''(k)
!       WHERE f(k),f''(k) ARE SPLINE COEFFICIENTS OF F(X)
!
!       NOTE: FOR ifunc = INTDER, the OHS factor in CmF cancels
!             THE HS factor in the Integral.
!             FOR ifunc = INTFUN, GRN is assumed to be pre-multiplied
!             OUTSIDE this routine by HS factor
!       ALSO, COMPUTES af(k) + (SUM on i)bs(i)*J(i,k) = jacob(k),
!       WHERE J(i,k) = d[g(i)]/d[f(k)], g = f''
!
!       nk(k): Number of smesh-pts in k-th spline interval
!       xknots(k) < smesh <= xknots(k+1),   k = 1,nots-1
!
!       NOTE: The ifunc=INTDER CASE is done by integrating by parts,
!             so that the half-point integration (GRN at half mesh pts)
!             becomes a full-point integration in Cm*F.
!
      nb = ideriv           !Pressure, iota derivatives vanish at origin
      ksp1 = nots - 1
      nmesh1 = nmesh - 1

      IF (ifunc .eq. intder) THEN
!
!       Integrate by parts (in s), func(1) and func(nmesh) are 'surface terms'
!
         func(1) = -cm(1)*grn(2)
         func(2:nmesh1) = cm(2:nmesh1)*(grn(2:nmesh1)-grn(3:nmesh1+1))
         func(nmesh) = cm(nmesh)*grn(nmesh)
      ELSE
         func = grn
      ENDIF

      af(:nots) = zero
      bs(:nots) = zero

      ia = 1
      DO k = 1, ksp1
         IF (nk(k) .ne. 0) THEN
            k1 = k + 1
            ib = ia + nk(k) - 1
            DO j = ia,ib
              af(k)  = af(k)  + func(j)*w(j)
              bs(k)  = bs(k)  + func(j)*u(j)
              af(k1) = af(k1) + func(j)*w1(j)
              bs(k1) = bs(k1) + func(j)*u1(j)
            ENDDO
            ia = ib + 1
         ENDIF
      END DO

      IF (ib .ne. nmesh) STOP 'ib!=nmesh'
      IF (nb .eq. natur) bs(1) = 0.           !Natural boundary conditions
      bs(nots) = 0.
      CALL jacprod (bs, h, nots, nb)         !Returns bs(i)=bs(j)*J(j,i)
      jacob(:nots) = af(:nots) + bs(:nots)

      END SUBROUTINE splinint
