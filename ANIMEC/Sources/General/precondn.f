      SUBROUTINE precondn(lu1, bsq, gsqrt, r12, xs, xu12, xue, xuo,
     1                    xodd, axm, axd, bxm, bxd, 
#ifdef _HBANGLE
     2                    axm2, axd2, bxm2, bxd2,
#endif
     3                    cx, eqfactor, trigmult)
      USE vmec_main
      USE vmec_params, ONLY: signgs
      USE realspace
#ifdef _HBANGLE
      USE angle_constraints, ONLY: cos2u, sin2u
#endif
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(nrzt), INTENT(in) ::
     1  lu1, bsq, gsqrt, r12, xs, xu12, xue, xuo, xodd
      REAL(rprec), DIMENSION(ns+1,2), INTENT(out) ::
     1  axm, axd, bxm, bxd
#ifdef _HBANGLE
      REAL(rprec), DIMENSION(ns+1,2), INTENT(out) ::
     1  axm2, axd2, bxm2, bxd2
#endif
      REAL(rprec), DIMENSION(ns+1), INTENT(out) :: cx
      REAL(rprec), DIMENSION(ns), INTENT(out) :: eqfactor
      REAL(rprec), DIMENSION(nznt), INTENT(in) :: trigmult
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, l, lk
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ax, bx
#ifdef _HBANGLE
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ax2, bx2
#endif
      REAL(rprec) :: temp(ns+1)
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: ptau, ptau2
      REAL(rprec) :: t1, t2, t3, pfactor
C-----------------------------------------------
!
!     COMPUTE PRECONDITIONING MATRIX ELEMENTS FOR R,Z
!     FORCE. NOTE THAT THE NORMALIZATION IS:
!
!     AX(off-diag) ~ <(cosmui cosmu cosnv cosnv) 2(R**2*Xu**2*bsq/gsqrt)>
!     Factor of 2 arising from 1/gsqrt**2 in bsq
!
!     Now, cosmui cosmu ~ mscale(0)**2, cosnv**2 ~ nscale(0)**2
!     Therefore, AX ~ (mscale(0)*nscale(0))**2 2<R**2*Xu**2*bsq/gsqrt>
!                   ~ 2*r0scale**2 <...>
!
      ALLOCATE (ax(ns+1,4), bx(ns+1,4), ptau(nznt), ptau2(nznt))
      ax = 0;  bx = 0;  cx = 0
      temp = 0
      pfactor = -2*r0scale**2       !v8.50
!     pfactor = -4*r0scale**2
#ifdef _HBANGLE
      ALLOCATE (ax2(ns+1,4), bx2(ns+1,4))
      ax2 = 0; bx2 = 0
#endif
      DO 20 js = 2,ns
!
!     COMPUTE DOMINANT (1/DELTA-S)**2 PRECONDITIONING
!     MATRIX ELEMENTS
!
        lk = 0
        DO l = js,nrzt,ns
          lk = lk + 1
          t1 = pfactor*r12(l)*bsq(l)
          ptau2(lk) = r12(l)*t1/gsqrt(l)
          t1 = t1*wint(l)
          temp(js) = temp(js) + t1*trigmult(lk)*xu12(l)
          ptau(lk) = r12(l)*t1/gsqrt(l)
          t1 = xu12(l)*ohs
          t2 = cp25*(xue(l)/shalf(js) + xuo(l))/shalf(js)
          t3 = cp25*(xue(l-1)/shalf(js) + xuo(l-1))/shalf(js)
          ax(js,1) = ax(js,1) + ptau(lk)*t1*t1
          ax(js,2) = ax(js,2) + ptau(lk)*(-t1+t3)*(t1+t2)
          ax(js,3) = ax(js,3) + ptau(lk)*(t1+t2)*(t1+t2)
          ax(js,4) = ax(js,4) + ptau(lk)*(-t1+t3)*(-t1+t3)
#ifdef _HBANGLE
          ptau2(lk) = ptau2(lk)*cos2u(lk)
          ax2(js,1) = ax2(js,1) + ptau2(lk)*t1*t1
          ax2(js,2) = ax2(js,2) + ptau2(lk)*(-t1+t3)*(t1+t2)
          ax2(js,3) = ax2(js,3) + ptau2(lk)*(t1+t2)*(t1+t2)
          ax2(js,4) = ax2(js,4) + ptau2(lk)*(-t1+t3)*(-t1+t3)
#endif
        END DO
!
!       COMPUTE PRECONDITIONING MATRIX ELEMENTS FOR M**2, N**2 TERMS
!
        lk = 0
        DO l = js,nrzt,ns
          lk = lk+1
          t1 = cp5*(xs(l) + cp5*xodd(l)/shalf(js))
          t2 = cp5*(xs(l) + cp5*xodd(l-1)/shalf(js))
          bx(js,1) = bx(js,1) + ptau(lk)*t1*t2
          bx(js,2) = bx(js,2) + ptau(lk)*t1*t1
          bx(js,3) = bx(js,3) + ptau(lk)*t2*t2
          cx(js) = cx(js) + cp25*pfactor*lu1(l)**2*gsqrt(l)*wint(l)
#ifdef _HBANGLE
          bx2(js,1) = bx2(js,1) + ptau2(lk)*t1*t2
          bx2(js,2) = bx2(js,2) + ptau2(lk)*t1*t1
          bx2(js,3) = bx2(js,3) + ptau2(lk)*t2*t2
#endif
        END DO
 20   CONTINUE

      temp(1)=0; temp(2:ns)=temp(2:ns)/vp(2:ns); temp(ns+1)=0
      DO js = 1,ns
        axm(js,1) =-ax(js,1)
        axd(js,1) = ax(js,1) + ax(js+1,1)
        axm(js,2) = ax(js,2) * sm(js) * sp(js-1)
        axd(js,2) = ax(js,3)*sm(js)**2 + ax(js+1,4)*sp(js)**2
        bxm(js,1) = bx(js,1)
        bxm(js,2) = bx(js,1) * sm(js) * sp(js-1)
        bxd(js,1) = bx(js,2) + bx(js+1,3)
        bxd(js,2) = bx(js,2)*sm(js)**2 + bx(js+1,3)*sp(js)**2
#ifdef _HBANGLE
        axm2(js,1) =-ax2(js,1)
        axd2(js,1) = ax2(js,1) + ax2(js+1,1)
        axm2(js,2) = ax2(js,2) * sm(js) * sp(js-1)
        axd2(js,2) = ax2(js,3)*sm(js)**2 + ax2(js+1,4)*sp(js)**2
        bxm2(js,1) = bx2(js,1)
        bxm2(js,2) = bx2(js,1) * sm(js) * sp(js-1)
        bxd2(js,1) = bx2(js,2) + bx2(js+1,3)
        bxd2(js,2) = bx2(js,2)*sm(js)**2 + bx2(js+1,3)*sp(js)**2
#endif
        cx(js)    = cx(js) + cx(js+1)
        temp(js) = signgs*(temp(js) + temp(js+1))
      END DO

      eqfactor(2:ns-1) = axd(2:ns-1,2)*hs*hs/temp(2:ns-1)
      eqfactor(1) = 0;  eqfactor(ns) = 0
      axm(ns+1,:) = 0;  axd(ns+1,:) = 0
      bxm(ns+1,:) = 0;  bxd(ns+1,:) = 0
      DEALLOCATE (ax, bx, ptau, ptau2)
#ifdef _HBANGLE
      axm2(ns+1,:) = 0;  axd2(ns+1,:) = 0
      bxm2(ns+1,:) = 0;  bxd2(ns+1,:) = 0
      DEALLOCATE (ax2, bx2)
#endif
      END SUBROUTINE precondn
