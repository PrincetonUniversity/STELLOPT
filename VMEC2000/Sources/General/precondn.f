#if defined (SKS)
      ! SKS-RANGE: All OUT arrays computed correctly between 
      ! [tlglob, trglob] 
      SUBROUTINE precondn_par(lu1, bsq, gsqrt, r12, xs, xu12, xue, xuo,
     1                    xodd, axm, axd, bxm, bxd, 
     3                    cx, eqfactor, trigmult)
      USE vmec_main
      USE vmec_params, ONLY: signgs
      USE realspace, ONLY: pshalf, pwint
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(nznt,ns), INTENT(in) ::
     1  lu1, bsq, gsqrt, r12, xs, xu12, xue, xuo, xodd
      REAL(rprec), DIMENSION(ns+1,2), INTENT(out) ::
     1  axm, axd, bxm, bxd
      REAL(rprec), DIMENSION(ns+1), INTENT(out) :: cx
      REAL(rprec), DIMENSION(ns), INTENT(out) :: eqfactor
      REAL(rprec), DIMENSION(nznt), INTENT(in) :: trigmult
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, l, lk
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ax, bx
      !REAL(rprec) :: temp(ns+1)
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: temp
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: ptau, ptau2
      REAL(rprec) :: t1, t2, t3, pfactor

      INTEGER :: nsmin, nsmax,  i, j, k, numjs
      INTEGER, DIMENSION(:), ALLOCATABLE :: ldisps, lcounts
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: sbuf

      nsmin=tlglob; nsmax=t1rglob
      ! Correct incoming lu1, bsq, gsqrt,   , xue, xuo, xodd, trigmult 
      ! between [t1lglob, t1rglob]

      ALLOCATE (ax(ns+1,4), bx(ns+1,4), ptau(nznt),
     1                     ptau2(nznt),temp(ns+1))
      ax = 0;  bx = 0;  cx = 0
      temp = 0
      pfactor = -4*r0scale**2        !restored in v8.51

      nsmin=MAX(2,tlglob); nsmax=t1rglob
      DO js = nsmin, nsmax
!
!     COMPUTE DOMINANT (1/DELTA-S)**2 PRECONDITIONING
!     MATRIX ELEMENTS
!
        lk = 0
        DO k = 1, ntheta3
          DO j = 1, nzeta
            lk = lk + 1
            t1 = pfactor*r12(lk,js)*bsq(lk,js)
            ptau2(lk) = r12(lk,js)*t1/gsqrt(lk,js)
            t1 = t1*pwint(lk,js)
            temp(js) = temp(js) + t1*trigmult(lk)*xu12(lk,js)
            ptau(lk) = r12(lk,js)*t1/gsqrt(lk,js)
            t1 = xu12(lk,js)*ohs
            t2=cp25*(xue(lk,js)/pshalf(lk,js)+xuo(lk,js))/pshalf(lk,js)
            t3 = cp25*(xue(lk,js-1)/pshalf(lk,js) + 
     1           xuo(lk,js-1))/pshalf(lk,js)
            ax(js,1) = ax(js,1) + ptau(lk)*t1*t1
            ax(js,2) = ax(js,2) + ptau(lk)*(-t1+t3)*(t1+t2)
            ax(js,3) = ax(js,3) + ptau(lk)*(t1+t2)*(t1+t2)
            ax(js,4) = ax(js,4) + ptau(lk)*(-t1+t3)*(-t1+t3)
          END DO
        END DO
!
!       COMPUTE PRECONDITIONING MATRIX ELEMENTS FOR M**2, N**2 TERMS
!
        lk = 0
        DO k = 1, ntheta3
          DO j = 1, nzeta
            lk = lk+1
            t1 = cp5*(xs(lk,js) + cp5*xodd(lk,js)/pshalf(lk,js))
            t2 = cp5*(xs(lk,js) + cp5*xodd(lk,js-1)/pshalf(lk,js))
            bx(js,1) = bx(js,1) + ptau(lk)*t1*t2
            bx(js,2) = bx(js,2) + ptau(lk)*t1*t1
            bx(js,3) = bx(js,3) + ptau(lk)*t2*t2
            cx(js) = cx(js) + cp25*pfactor*lu1(lk,js)**2
     1               *gsqrt(lk,js)*pwint(lk,js)
          END DO
        END DO

      END DO

      nsmin=MAX(2,tlglob); nsmax=t1rglob
      temp(1)=0; temp(nsmin:nsmax)=temp(nsmin:nsmax)/vp(nsmin:nsmax); 
      temp(ns+1)=0

      nsmin=t1lglob; nsmax=t1rglob
!      nsmin=t1lglob; nsmax=trglob
      DO js = nsmin, nsmax
        axm(js,1) =-ax(js,1)
        axd(js,1) = ax(js,1) + ax(js+1,1)
        axm(js,2) = ax(js,2) * sm(js) * sp(js-1)
        axd(js,2) = ax(js,3)*sm(js)**2 + ax(js+1,4)*sp(js)**2
        bxm(js,1) = bx(js,1)
        bxm(js,2) = bx(js,1) * sm(js) * sp(js-1)
        bxd(js,1) = bx(js,2) + bx(js+1,3)
        bxd(js,2) = bx(js,2)*sm(js)**2 + bx(js+1,3)*sp(js)**2
        cx(js)    = cx(js) + cx(js+1)
        temp(js) = signgs*(temp(js) + temp(js+1))
      END DO

      nsmin=MAX(2,tlglob); nsmax=MIN(ns-1,trglob)
      eqfactor(nsmin:nsmax) = axd(nsmin:nsmax,2)*hs*hs/
     1                        temp(nsmin:nsmax)
      eqfactor(1) = 0;  eqfactor(ns) = 0
      axm(ns+1,:) = 0;  axd(ns+1,:) = 0
      bxm(ns+1,:) = 0;  bxd(ns+1,:) = 0
      DEALLOCATE (ax, bx, ptau, ptau2, temp)

      END SUBROUTINE precondn_par
#endif

      SUBROUTINE precondn(lu1, bsq, gsqrt, r12, xs, xu12, xue, xuo,
     1                    xodd, axm, axd, bxm, bxd, 
     3                    cx, eqfactor, trigmult)
      USE vmec_main
      USE vmec_params, ONLY: signgs
      USE realspace
#if defined (SKS)
      USE parallel_include_module
#endif
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(nrzt), INTENT(in) ::
     1  lu1, bsq, gsqrt, r12, xs, xu12, xue, xuo, xodd
      REAL(rprec), DIMENSION(ns+1,2), INTENT(out) ::
     1  axm, axd, bxm, bxd
      REAL(rprec), DIMENSION(ns+1), INTENT(out) :: cx
      REAL(rprec), DIMENSION(ns), INTENT(out) :: eqfactor
      REAL(rprec), DIMENSION(nznt), INTENT(in) :: trigmult
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, l, lk
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ax, bx
      REAL(rprec) :: temp(ns+1)
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: ptau, ptau2
      REAL(rprec) :: t1, t2, t3, pfactor
#if defined(SKS)
      INTEGER :: nsmin, nsmax,  i, j, k
#endif
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
#if defined(SKS)      
      nsmin=tlglob; nsmax=t1rglob
#endif
      ALLOCATE (ax(ns+1,4), bx(ns+1,4), ptau(nznt), ptau2(nznt))
      ax = 0;  bx = 0;  cx = 0
      temp = 0
!      pfactor = -2*r0scale**2       !v8.50
      pfactor = -4*r0scale**2        !restored in v8.51
#if defined(SKS)      
      nsmin=MAX(2,tlglob); nsmax=t1rglob
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
        cx(js)    = cx(js) + cx(js+1)
        temp(js) = signgs*(temp(js) + temp(js+1))
      END DO

      eqfactor(2:ns-1) = axd(2:ns-1,2)*hs*hs/temp(2:ns-1)
      eqfactor(1) = 0;  eqfactor(ns) = 0
      axm(ns+1,:) = 0;  axd(ns+1,:) = 0
      bxm(ns+1,:) = 0;  bxd(ns+1,:) = 0
      DEALLOCATE (ax, bx, ptau, ptau2)

      END SUBROUTINE precondn
