      SUBROUTINE svdsolveuv(m, n, mp, np, wcut, file, a, b, x)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER m, n, mp, np
      REAL(rprec) wcut
      CHARACTER file*(*)
      REAL(rprec), DIMENSION(mp,np) :: a
      REAL(rprec), DIMENSION(n) :: x
      REAL(rprec), DIMENSION(m) :: b
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec) :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: n1, nonzero, nlast, istat, i, j, iunit
      REAL(rprec), ALLOCATABLE :: u(:,:), w(:), v(:,:)
      REAL(rprec) :: wmax, wcuta, udotb
C-----------------------------------------------

c       Solves Matrix equation A X = B for X using SVD method
c       Allows supression of weights, storage and recall of answers
c       Uses svd routines from numerical recipes
c       Prashant Valanju (Sept 1998) pvalanju@mail.utexas.edu

c       Almost same as SvdSolveB except this one stores U, W, V
c       rather than basis vectors.
c       Use this one when you want to change both B and weights for fixed A

c       Inputs:
c       A(M,N) - Matrix with physical dimensions (Mp,Np)
c       B(N)   - R.H.S. of A X = B
c       M      - Number of rows of A to use
c       N      - Number of columns of A to use
c       Mp,NP  - Storage dimensions of A(Mp,Np)
c       Wcut   - Cutoff parameter
c                If Wcut = 0, writes answers into "file", all weights are kept
c                If Wcut < 0, reads previous answers from "file",
c                              and uses -Wcut for weights to keep
c                If 0 < Wcut < 1 : Cuts off weights with w(i)/wmax < wcut
c                      Use this to cut off small weights, does calculation
c                If Wcut = integer > 1, keeps Wcut weights, does calculation
c                      Use this to cut off fixed number of weights
c       Output:
c       X(N)   - Solution of A X = B
c
      n1 = n

c     ALLOCATE local arrays
      ALLOCATE (u(mp,np), w(np), v(np,np), stat=istat)
      IF (istat .eq. 0) THEN                       !have enough memory
         iunit = 41

         x(:n1) = zero                    !zero answer from previous CALL
         wcuta = ABS(wcut)
c...............................
c     If Wcut < 0 use precalculated weights to get answer fast
c        use wcuta = ABS(wcut) for everything
c     If any error happens while reading "file", do full svd calculation
         IF (wcut .lt. 0) THEN         !Read file from previous calculation
            CALL safe_open(iunit, istat, file, 'old', 'formatted')
            IF (istat .ne. 0) GOTO 98
            READ (iunit, *, err=98)
            READ (iunit, *, err=98)
            READ (iunit, *, err=98) m, n, nonzero
            DO i = 1, n     !Read exactly the same way they were written
               READ (iunit, *, err=98) w(i)         !READ i-th weight
               DO j = 1, n
                  READ (iunit, *, err=98) v(j,i)
               END DO
               DO j = 1, m
                  READ (iunit, *, err=98) u(j,i)
               END DO
            END DO
            CLOSE(unit=iunit)

            GOTO 99                !Use U, w, V matrices to calculate X
         ENDIF                     !End of precalculated branch Wcut < 0
   98    CONTINUE
c.......................................
c     If Wcut .ge. 0, do full svd calculation
c     Initialize all to zero to wipe out effects of old call
         DO j = 1, n
            w(j) = zero        !Zero all weights
            u(:m,j) = a(:m,j)  !Because U will be changed by svdcmp
            v(:n,j) = zero
         END DO

         CALL svdcmp (u, m, n, mp, np, w, v)     !Do SVD decomposition

c       Sort weights and matrices with DECREASING weights
c       Permute weight w(i) AND column vectors U(*,i), V(*,i) at the same time
         CALL sortsvd (m, n, mp, np, w, u, v)

c       Find the number of nonzero weights (already dcreasing ordered)
         DO nonzero = n, 1, -1
c                                   !Found first nonzero weight, get out
            IF (w(nonzero) .ne. 0) EXIT
         END DO
         IF (nonzero .le. 0) GOTO 999

c............................................
c        Decide how many weights to keep
   99    CONTINUE
         IF (wcuta .gt. 1) THEN
            nlast = wcuta
         ELSE                                    !cutoff small weights
            wMAX = w(1)                     !weights are already ordered
            nlast = 0
            DO i = 1, nonzero
               IF (w(i) .gt. wmax*wcuta) THEN
                  nlast = nlast + 1              !accept this weight
               ELSE
                  GOTO 96
               END IF
            END DO
         END IF

c       Find solution X and basis vectors (put into columns of V)
c       This is the NR SVBKSB routine
   96    CONTINUE
         IF (nlast.le.0 .or. nlast.gt.nonzero) nlast = nonzero
         DO i = 1, nlast     !Calc i-th coeff (U.b/w) for noNonzero w(i)
c        First calculate i-th coeff (U(i)*b/w(i)) in the SUM in eq 14.3.17 (NR)
c                            !Dot product of i-th column of U with B
c                            !summed over all plasma points
            udotb = SUM(u(:m,i)*b(:m))
c                        !=(Ui.b/wi) in eq 14.3.17, saves many divisions
            udotb = udotb/w(i)
c         Now run DOwn the i-th column of the vector V(j,i)
c                         !j-th entry (row) of vectors in eq 14.3.17
c                         !add it to the j-th entry of X
            x(:n) = x(:n) + udotb*v(:n,i)
c                                  !This did the svdbksb, faster than NR
         END DO                   !and also calculated the basis vectors

c............................................
c       Write ALL weights and U, V into SVD file (IF Wcut=0)
         IF (wcut .eq. 0) THEN
            CALL safe_open(iunit, istat, file, 'unknown', 'formatted')
            WRITE (iunit, *) 'Max w = ', w(1), ', Min w = ', w(nonzero)
            WRITE (iunit, *) 'Ratio = ', w(1)/w(nonzero)
            WRITE (iunit, *) m, n, nonzero
            DO i = 1, n
               WRITE (iunit, *) w(i)                !WRITE i-th weight
               DO j = 1, n
                  WRITE (iunit, *) v(j,i)
               END DO
               DO j = 1, m
                  WRITE (iunit, *) u(j,i)
               END DO
            END DO
            CLOSE(unit=iunit)
         END IF
c................................................
      END IF
  999 CONTINUE

      DEALLOCATE (u, w, v, stat=istat)

      END SUBROUTINE svdsolveuv
