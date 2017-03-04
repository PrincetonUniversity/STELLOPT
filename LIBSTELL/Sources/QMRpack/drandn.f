C**********************************************************************
C
      SUBROUTINE DRANDN (N,DX,SEED)
      USE stel_kinds, ONLY: rprec
      IMPLICIT NONE
C
C     Purpose:
C     Fills the vector DX with random numbers  between 0 and 1.  If the
C     SEED is given, it should be odd and positive.  The generator is a
C     fairly unsophisticated one, from Pearson's  "Numerical methods in
C     engineering and science" book.
C
C     Parameters:
C     N    = the dimension of the vector (input).
C     DX   = the vector to fill with random numbers (output).
C     SEED = the seed for the generator (input).
C
C     Noel M. Nachtigal
C     April 23, 1993
C
C**********************************************************************
C
      INTRINSIC DBLE, ABS, MOD
C
      INTEGER N, SEED
      REAL(rprec) DX(N)
C
C     Local variables.
C
      INTEGER I, J
C
C     Local variables that are saved from one call to the next.
C
      REAL(rprec) DMAX
      INTEGER :: IM=0, IMAX, IS
      SAVE DMAX, IM, IMAX, IS
C      DATA IM/0/
C
C     Initialize the generator data.
C
      IF (IM.EQ.0) THEN
         J  = 0
         IM = 1
         DO 10 I = 1, 31
            J = J + 1
            IF (IM*2.LE.IM) GO TO 20
            IM = IM * 2
 10      CONTINUE
 20      IMAX = (IM-1) * 2 + 1
         DMAX = IMAX
         DO 30 I = 1, MOD(J,3)
            J = J - 1
            IM = IM / 2
 30      CONTINUE
         IM = IM + 5
         IS = ABS(MOD(IM*30107,IMAX))
      END IF
C
C     Check whether we have a new seed.
C
      IF (SEED.GT.0) IS = (SEED / 2) * 2 + 1
C
C     Here goes the rest.
C
      DO 40 I = 1, N
         DX(I) = DBLE(IS) / DMAX
         IS    = ABS(MOD(IM*IS,IMAX))
 40   CONTINUE
C
      RETURN
      END
C
C**********************************************************************
