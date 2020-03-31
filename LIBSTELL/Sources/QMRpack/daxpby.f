**********************************************************************
C
C     Copyright (C) 1992  Roland W. Freund and Noel M. Nachtigal
C     All rights reserved.
C
C     This code is part of a copyrighted package.  For details, see the
C     file `cpyrit.doc' in the top-level directory.
C
C     *****************************************************************
C     ANY USE OF  THIS CODE CONSTITUTES ACCEPTANCE OF  THE TERMS OF THE
C                             COPYRIGHT NOTICE
C     *****************************************************************
C
C**********************************************************************
C
      SUBROUTINE DAXPBY (N,DZ,DA,DX,DB,DY)
      USE stel_kinds, ONLY: rprec, dp
      IMPLICIT NONE
C
C     Purpose:
C     This subroutine computes DZ = DA * DX + DB * DY.  Several special
C     cases are handled separately:
C        DA =  0.0, DB =  0.0 => DZ = 0.0
C        DA =  0.0, DB =  1.0 => DZ = DY  (this is COPY)
C        DA =  0.0, DB = -1.0 => DZ = -DY
C        DA =  0.0, DB =   DB => DZ = DB * DY  (this is SCAL)
C        DA =  1.0, DB =  0.0 => DZ = DX  (this is COPY)
C        DA =  1.0, DB =  1.0 => DZ = DX + DY
C        DA =  1.0, DB = -1.0 => DZ = DX - DY
C        DA =  1.0, DB =   DB => DZ = DX + DB * DY (this is AXPY)
C        DA = -1.0, DB =  0.0 => DZ = -DX
C        DA = -1.0, DB =  1.0 => DZ = -DX + DY
C        DA = -1.0, DB = -1.0 => DZ = -DX - DY
C        DA = -1.0, DB =   DB => DZ = -DX + DB * DY
C        DA =   DA, DB =  0.0 => DZ = DA * DX  (this is SCAL)
C        DA =   DA, DB =  1.0 => DZ = DA * DX + DY  (this is AXPY)
C        DA =   DA, DB = -1.0 => DZ = DA * DX - DY
C        DA =   DA, DB =   DB => DZ = DA * DX + DB * DY
C     DZ may be the same as DX or DY.
C
C     Parameters:
C     N  = the dimension of the vectors (input).
C     DZ = the vector result (output).
C     DA = scalar multiplier for DX (input).
C     DX = one of the vectors (input).
C     DB = scalar multiplier for DY (input).
C     DY = the other vector (input).
C
C     Noel M. Nachtigal
C     March 23, 1993
C
C**********************************************************************
C
      INTEGER N
      REAL(rprec), PARAMETER :: ZERO = 0, ONE = 1
      REAL(rprec) DA, DB, DX(N), DY(N), DZ(N)
C
C     Local variables.
C
      INTEGER I
C
      IF (N.LE.0) RETURN
C
      IF (DA.EQ.ZERO) THEN
         IF (DB.EQ.ZERO) THEN
C           DA = 0.0, DB = 0.0 => DZ = 0.0.
            DO 10 I = 1, N
               DZ(I) = ZERO
 10         CONTINUE
         ELSE IF (DB.EQ.ONE) THEN
C           DA = 0.0, DB = 1.0 => DZ = DY (this is COPY).
            DO 20 I = 1, N
               DZ(I) = DY(I)
 20         CONTINUE
         ELSE IF (DB.EQ.-ONE) THEN
C           DA = 0.0, DB = -1.0 => DZ = -DY.
            DO 30 I = 1, N
               DZ(I) = -DY(I)
 30         CONTINUE
         ELSE
C           DA = 0.0, DB = DB => DZ = DB * DY (this is SCAL).
            DO 40 I = 1, N
               DZ(I) = DB * DY(I)
 40         CONTINUE
         END IF
      ELSE IF (DA.EQ.ONE) THEN
         IF (DB.EQ.ZERO) THEN
C           DA = 1.0, DB = 0.0 => DZ = DX (this is COPY).
            DO 50 I = 1, N
               DZ(I) = DX(I)
 50         CONTINUE
         ELSE IF (DB.EQ.ONE) THEN
C           DA = 1.0, DB = 1.0 => DZ = DX + DY.
            DO 60 I = 1, N
               DZ(I) = DX(I) + DY(I)
 60         CONTINUE
         ELSE IF (DB.EQ.-ONE) THEN
C           DA = 1.0, DB = -1.0 => DZ = DX - DY.
            DO 70 I = 1, N
               DZ(I) = DX(I) - DY(I)
 70         CONTINUE
         ELSE
C           DA = 1.0, DB = DB => DZ = DX + DB * DY (this is AXPY).
            DO 80 I = 1, N
               DZ(I) = DX(I) + DB * DY(I)
 80         CONTINUE
         END IF
      ELSE IF (DA.EQ.-ONE) THEN
         IF (DB.EQ.ZERO) THEN
C           DA = -1.0, DB = 0.0 => DZ = -DX
            DO 90 I = 1, N
               DZ(I) = -DX(I)
 90         CONTINUE
         ELSE IF (DB.EQ.ONE) THEN
C           DA = -1.0, DB = 1.0 => DZ = -DX + DY
            DO 100 I = 1, N
               DZ(I) = -DX(I) + DY(I)
 100        CONTINUE
         ELSE IF (DB.EQ.-ONE) THEN
C           DA = -1.0, DB = -1.0 => DZ = -DX - DY.
            DO 110 I = 1, N
               DZ(I) = -DX(I) - DY(I)
 110        CONTINUE
         ELSE
C           DA = -1.0, DB = DB => DZ = -DX + DB * DY
            DO 120 I = 1, N
               DZ(I) = -DX(I) + DB * DY(I)
 120        CONTINUE
         END IF
      ELSE
         IF (DB.EQ.ZERO) THEN
C           DA = DA, DB = 0.0 => DZ = DA * DX (this is SCAL).
            DO 130 I = 1, N
               DZ(I) = DA * DX(I)
 130        CONTINUE
         ELSE IF (DB.EQ.ONE) THEN
C           DA = DA, DB = 1.0 => DZ = DA * DX + DY (this is AXPY)
            DO 140 I = 1, N
               DZ(I) = DA * DX(I) + DY(I)
 140        CONTINUE
         ELSE IF (DB.EQ.-ONE) THEN
C           DA = DA, DB = -1.0 => DZ = DA * DX - DY.
            DO 150 I = 1, N
               DZ(I) = DA * DX(I) - DY(I)
 150        CONTINUE
         ELSE
C           DA = DA, DB = DB => DZ = DA * DX + DB * DY.
            DO 160 I = 1, N
               DZ(I) = DA * DX(I) + DB * DY(I)
 160        CONTINUE
         END IF
      END IF
C
      RETURN
      END
C
C**********************************************************************
