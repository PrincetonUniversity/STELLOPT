      SUBROUTINE restart(irst)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER irst
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, SAVE :: nresets = 0
C-----------------------------------------------
!
!       This routine either stores an accepted value of the local solution
!       (irst = 1) or reset the PRESENT solution to a previous value (irst = 2)
!
      SELECT CASE (irst)

      CASE DEFAULT
         xstore(:n2) = xvec(:n2)
         RETURN

      CASE (2)
         xdot(:n2) = zero
         xvec(:n2) = xstore(:n2)
         delt = .95*delt
         irst = 1
         nresets = nresets + 1
         IF (nresets .ge. 1E6) THEN
            PRINT *, ' Time step reduced 1E6 times without convergence'
            STOP
         ENDIF
         RETURN
      END SELECT

      END SUBROUTINE restart
