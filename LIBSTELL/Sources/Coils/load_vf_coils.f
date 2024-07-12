      SUBROUTINE load_vf_coils (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vf_coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, k, nsv
      INTEGER :: nvariables
      REAL(rprec) :: xvariables(*)
!-----------------------------------------------

      nsv = 0

      IF (lvfvar) THEN
!        DO i=1, num_vf
         DO i=nvf_fix+1, num_vf
            IF (lvfz) THEN
               nsv = nsv + 1
               zc_vf(i) = xvariables(nsv)
            END IF
            IF (lvfr) THEN
               nsv = nsv + 1
               rc_vf(i) = xvariables(nsv)
            END IF
            IF (nrvf_c .gt. 0) THEN
               DO k=1, nrvf_c
                  nsv = nsv + 1
                  rcfc_vf(i,k) = xvariables(nsv)
                  nsv = nsv + 1
                  rcfs_vf(i,k) = xvariables(nsv)
               END DO
            END IF
         END DO
      END IF

      nvariables = nsv

      END SUBROUTINE load_vf_coils
