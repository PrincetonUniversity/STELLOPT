      SUBROUTINE tolicu (torcur)
      USE vparams, ONLY: mu0
      USE vacmod
      USE biotsavart
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), INTENT(IN) :: torcur
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, kper, kv
      REAL(dp) :: current(1), ttolion, ttolioff
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: xpts
C-----------------------------------------------
!
!     COMPUTE WIRE SEGMENTS (DO NOT CLOSE LOOP, CLOSURE DONE IN biotsavart ROUTINES)
!
      CALL second0(ttolion)

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

      CALL second0(ttolioff)
      s_tolicu_time = s_tolicu_time + (ttolioff - ttolion)

      END SUBROUTINE tolicu
