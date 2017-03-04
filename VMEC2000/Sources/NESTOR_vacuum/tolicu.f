      SUBROUTINE tolicu(torcur)
      USE vparams, ONLY: mu0
      USE vacmod
      USE biotsavart
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), INTENT(in) :: torcur
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, kper, kv
      REAL(rprec) :: current(1)
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: xpts
C-----------------------------------------------
!
!     COMPUTE WIRE SEGMENTS (DO NOT CLOSE LOOP, CLOSURE DONE IN biotsavart ROUTINES)
!
      ALLOCATE (xpts(3,nvp), stat=i)
      IF (i .ne. 0) STOP ' allocation error in tolicu'

      current = torcur/mu0

      i = 0
      DO kper = 1, nvper
         DO kv = 1, nv
            i = i + 1
            xpts(1,i) = raxis_nestor(kv)*(cosper(kper)*cosuv(kv)
     1                -                   sinper(kper)*sinuv(kv))
            xpts(2,i) = raxis_nestor(kv)*(sinper(kper)*cosuv(kv) 
     1                +            cosper(kper)*sinuv(kv))
            xpts(3,i) = zaxis_nestor(kv)
         END DO
      END DO

 
!
!     INITIALIZE COIL-RELATED QUANTITIES
!     
      CALL initialize_biotsavart (current, xpt=xpts)

      END SUBROUTINE tolicu
