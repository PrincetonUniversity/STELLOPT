      SUBROUTINE load_tf_coils (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE tf_coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, nvariables
      REAL(rprec) :: xvariables(*)
!-----------------------------------------------
      nvariables = 1
      i_tfc = xvariables(1)
      tfc_cur(1) = xvariables(1)

      IF (lqos) THEN
         DO i = 1, mtfcoil
            tfc_cur(i) = xvariables(1)/mtfcoil
         END DO
      END IF

      END SUBROUTINE load_tf_coils
