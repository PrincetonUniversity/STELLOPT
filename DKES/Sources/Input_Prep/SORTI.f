      SUBROUTINE SORTI(IX, NM)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER NM
      INTEGER, DIMENSION(*) :: IX
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ML, MG, IIL, IIG, ML1, MG1, NL, NG, N, II, IX0
C-----------------------------------------------

Ccom*    Sorting of an Integer Array (size)
C----------------------------------------------------------------------
C     Sorting of an INTEGER array. On EXIT, the elements of IX increase
C     with INDEX i.
C
C     IX      INTEGER array
C     NM      no. of elements
C----------------------------------------------------------------------
Cend
C
      IF (NM <= 1) RETURN
      ML = 1
      MG = NM
      IF (NM .ne.2) THEN
C
    1    CONTINUE
         IIL = IX(ML)
         IIG = IIL
         ML1 = ML + 1
         MG1 = MG - 1
         NL = ML
         NG = ML
C
         DO N = ML1, MG
            II = IX(N)
            IF (II < IIL) THEN
               NL = N
               IIL = II
            ENDIF
            IF (II > IIG) THEN
               NG = N
               IIG = II
            ENDIF
         END DO
C
         IF (NL .ne.ML) THEN
            IX0 = IX(NL)
            IX(NL) = IX(ML)
            IX(ML) = IX0
         ENDIF
         IF (NG .ne.MG) THEN
            IF (NG .eq. ML) NG = NL
            IX0 = IX(NG)
            IX(NG) = IX(MG)
            IX(MG) = IX0
         ENDIF
         ML = ML1
         MG = MG1
         IF (MG - ML >= 2) GO TO 1
C
         IF (MG <= ML) RETURN
      ENDIF
      IF (IX(MG) >= IX(ML)) RETURN
      IX0 = IX(MG)
      IX(MG) = IX(ML)
      IX(ML) = IX0

      RETURN

      END SUBROUTINE SORTI
