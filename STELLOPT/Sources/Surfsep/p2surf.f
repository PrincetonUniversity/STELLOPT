      SUBROUTINE p2surf(u, v, d, ierr)
      USE stel_kinds
      USE geom
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER ierr
      REAL(dp) u, v, d
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: n = 2
      INTEGER, PARAMETER :: ldfjac = 2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: tries
      REAL(dp), DIMENSION(n) :: xx
      REAL(dp) :: step2, r, step1, x, y, z, phi, factor
      LOGICAL :: check
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      EXTERNAL newt
!-----------------------------------------------
!
      ierr = 0
      CALL random_seed                           ! initialize the random number generator
      factor = 100
      xx(1) = u
      xx(2) = v
!
      check = .FALSE.                            ! presume success to start
      DO tries = 1, 20
!
         CALL newt (xx, n, check)
!
         IF (.not.check) THEN
            EXIT
         ELSE
!            WRITE(*,*) 'check = ', check
            CALL random_number (step1)
            CALL random_number (step2)
            IF (step1 <= 0.5_dp) step1 = -step1
            IF (step2 <= 0.5_dp) step2 = -step2
            xx(1) = (1 + step1)*xx(1)
            xx(2) = (1 + step2)*xx(2)
         ENDIF
      END DO
!
!
!     If successful get point r,z on plasfree at u,v
      IF (.not.check) THEN
         CALL dmnf1 (xx(1), xx(2), r, z)

         u = xx(1)
         v = xx(2)

         phi = alp*v
         x = r*COS(phi)
         y = r*SIN(phi)
         d = SQRT((x - xp)**2 + (y - yp)**2 + (z - zp)**2)
         ierr = 0
      ELSE
         ierr = 1
      ENDIF

      END SUBROUTINE p2surf
