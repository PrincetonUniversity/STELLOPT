      SUBROUTINE add_tension(amat, wten, hx, tens, tensv, fpoly, 
     1   n, nb, ioff, nmat)
      USE vspline
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER n, nb, ioff, nmat
      REAL(rprec) tens, tensv, fpoly
      REAL(rprec), DIMENSION(nmat,*) :: amat
      REAL(rprec), DIMENSION(nmat) :: wten
      REAL(rprec), DIMENSION(n) :: hx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: k, i, koff
      REAL(rprec), DIMENSION(n) :: wten0, tenshx, work
      REAL(rprec) :: delta_x, tension
C-----------------------------------------------

!
!       THIS SUBROUTINE ADDS THE TENSION TERM WTEN*(DEL(K) - DEL(K-1))
!       TO THE REST OF THE CHI-SQ AMAT, WHERE DEL(K) = (G(K+1)-G(K))/H(K)
!       AND G(K) IS THE SECOND DERIVATIVE AT THE K-TH KNOT
!
!       IOFF ALLOWS FOR ADDING TENSION SEPARATELY TO IOTA (IOFF=0)
!       AND PRESSURE (IOFF=NIOTA) SPLINES
!
!       tens:   spline tension (optionally, at the 1st pt only;see note)
!       tensv:  vbl spline tension for n-1th point (see note)
!       fpoly:  vbl spline tension form factor (note: if tens1<>tens2
!               then tension(i-th point) = tens+(tensv-tens)*(i/n-1))**fpoly)

!
!       BOUNDS CHECKING
!
      IF (n + ioff > nmat) STOP '(n+ioff>nmat)'
      IF (fpoly < 0.) STOP '(fpoly<0)'
      IF (n < 1) STOP '(n < 1)'

!
!       COMPUTE TENSION COEFFICIENTS
!

      delta_x = SUM(hx(:n-1))
      tension = 0.5*(delta_x/(n))**3
      IF (fpoly.eq.0. .or. tens.eq.tensv) THEN
         tenshx(:n-1) = tens
      ELSE
         DO i = 1, n - 1
            tenshx(i) = tens + (tensv - tens)*(REAL(i - 1,rprec)/
     1         (n - 1))**fpoly
         END DO
      ENDIF

      DO i = 1,n-1
        tenshx(i) = tension * tenshx(i) / hx(i)
        work(i) = hx(i)*(wten(i+ioff) + wten(i+ioff+1))
      ENDDO
      DO i = 2,n-1
        wten0(i) = 0.5 * ( work(i) + work(i-1) )/(hx(i) + hx(i-1))
      ENDDO
      wten0(1) = wten0(2)
      wten0(n) = wten(n+ioff)
!
!       COMPUTE, FOR K = 1,N, B(K,L)*JACOBIAN(L,I) = W(K,I),
!       WHERE JACOBIAN = D[G]/D[F] and B is TRIDIAGONAL
!       SEE EQN(27) IN PHYS.PLASMAS 1, p 2277.
!
      DO k = 1, n
         koff = k + ioff
         work(:n-1) = 0
!       SET UP COEFFICIENTS IN [G(K+1)-G(K)]/h(k) - [G(K)-G(K-1)]/h(k-1)
         IF (k .eq. 1) THEN
            work(2) = tenshx(1)*wten0(2)
            work(1) = -work(2)
         ELSE IF (k .eq. n) THEN
            work(n-1) = tenshx(n-1)*wten0(n-1)
         ELSE
            work(k-1) = tenshx(k-1)*wten0(k-1)
            work(k) = -(tenshx(k)+tenshx(k-1))*wten0(k)
            work(k+1) = tenshx(k)*wten0(k+1)
         ENDIF
         IF (nb .eq. natur) work(1) = 0
         work(n) = 0
!
!       COMPUTE work(j) = work(i)*Jacobian(i,j) and add to amat(k,j)
!
         CALL jacprod (work, hx, n, nb)
         amat(koff,1+ioff:n+ioff) = amat(koff,1+ioff:n+ioff) + work(:n)
      END DO

      END SUBROUTINE add_tension
