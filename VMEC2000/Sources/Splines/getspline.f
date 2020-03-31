      SUBROUTINE getspline(amat, splnots, hk, delse, hs, indexs, isort,
     1   ndata0, nots)
      USE vsvd0
      USE vparams, ONLY: zero, epstan
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER ndata0, nots
      REAL(rprec) hs
      INTEGER, DIMENSION(ndata0) :: indexs, isort
      REAL(rprec), DIMENSION(nots,ndata0) :: amat
      REAL(rprec), DIMENSION(nots) :: splnots, hk
      REAL(rprec), DIMENSION(ndata0) :: delse
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, DIMENSION(nots) :: nk
      INTEGER :: i, js, ia, k, k1, ib, j, nb
      REAL(rprec), DIMENSION(ndata0) :: w, w1, u, u1, snodes
      REAL(rprec), DIMENSION(nots,ndata0) :: bmat
C-----------------------------------------------

!
!       ON EXIT, AMAT = AMAT + BMAT, WHERE BMAT IS THE 2nd
!       DERIVATIVE COEFFICIENT MATRIX ARRAY, MULTIPLIED BY JACOBIAN
!       AND AMAT (ON RHS) WAS FUNCTION COEFFICIENT MATRIX ARRAY
!

!
!       SORT KNOT POSITIONS IN ASCENDING ORDER IN S-SPACE
!       USE SQRT(S) KNOT POSITIONS FOR IMPROVED RESOLUTION
!       NOTE: SNODES(I) IS THE VALUE OF SQRT(S) CORRESPONDING TO
!       THE MESH VALUES CORRESPONDING TO DELSE, INDEXS (COMPUTED OUTSIDE
!       THIS PROGRAM)
!

      DO i = 1, ndata0
         js = indexs(i)
         snodes(i) = SQRT(hs*((js - 1) + delse(i)))
!        IF (snodes(i) .le. zero) snodes(i) = epstan
      END DO

!     Avoid roundoff error in SPLININT
      IF( snodes(ndata0) .gt. splnots(nots) )
     1    snodes(ndata0) = splnots(nots)

      CALL sort_data (snodes, isort, ndata0)

!
!       COMPUTE MATRIX COEFFICIENTS RELATING SPLINE AT SPLNOTS
!       TO REAL-SPACE FUNCTION AT SORTED MESH POINTS RMESH
!
      amat(:nots,:ndata0) = zero
      bmat(:nots,:ndata0) = zero

!
!       SETUP SPLINE PARAMETERS AT EACH TIME STEP, SINCE SNODES
!       MAY BE CHANGING DYNAMICALLY IN S-SPACE
!
      CALL setup_int(splnots,snodes,hk,w,w1,u,u1,nk,nots,ndata0)

      ia = 1
      DO k = 1, nots - 1
         IF (nk(k) .gt. 0) THEN
            k1 = k + 1
            ib = ia + nk(k) - 1
            amat(k,ia:ib) = amat(k,ia:ib) + w(ia:ib)
            bmat(k,ia:ib) = bmat(k,ia:ib) + u(ia:ib)
            amat(k1,ia:ib) = amat(k1,ia:ib) + w1(ia:ib)
            bmat(k1,ia:ib) = bmat(k1,ia:ib) + u1(ia:ib)
            ia = ib + 1
         ENDIF
      END DO

      IF (ib .ne. ndata0) STOP 'ib!=ndat'
      nb = ideriv

      DO j = 1, ndata0
         bmat(nots,j) = 0.
         CALL jacprod (bmat(1,j), hk, nots, nb)
         amat(:nots,j) = amat(:nots,j) + bmat(:nots,j)
      END DO

      END SUBROUTINE getspline
