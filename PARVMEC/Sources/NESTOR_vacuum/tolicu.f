#if defined(SKS)      
      SUBROUTINE tolicu_par(torcur)
      USE vparams, ONLY: mu0
      USE vacmod
      USE biotsavart
      USE parallel_include_module
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
      REAL(rprec) :: skston, skstoff
C-----------------------------------------------
!
!     COMPUTE WIRE SEGMENTS (DO NOT CLOSE LOOP, CLOSURE DONE IN biotsavart ROUTINES)
!
      CALL second0(skston)

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

      CALL second0(skstoff)
      s_tolicu_time = s_tolicu_time + (skstoff - skston)

      END SUBROUTINE tolicu_par
#endif


      SUBROUTINE tolicu(torcur)
      USE vparams, ONLY: mu0
      USE vacmod
      USE biotsavart
#if defined(SKS)
      USE parallel_include_module
#endif
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
#if defined(SKS)
      REAL(rprec) :: skston, skstoff
#endif

C-----------------------------------------------
!
!     COMPUTE WIRE SEGMENTS (DO NOT CLOSE LOOP, CLOSURE DONE IN biotsavart ROUTINES)
!
#if defined(SKS)
      CALL second0(skston)
#endif
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

#if defined(SKS)
      CALL second0(skstoff)
      s_tolicu_time = s_tolicu_time + (skstoff - skston)
#endif

      END SUBROUTINE tolicu
