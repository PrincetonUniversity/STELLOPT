      SUBROUTINE load_saddle_coils (nvariables, xvariables)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE saddle_coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, n, nsv, modes
      INTEGER :: nvariables
      REAL(rprec) :: xvariables(*)
!-----------------------------------------------

!     load the unique coil parameters with values from variables in
!     optimization

      nvariables = 0
      nsv = 0

      IF (.not.lsadshape) RETURN

      IF (lspline) THEN

         ! load bspline series representation
         DO i = 1, nsmid
            ! coefficients for v
            DO n = 1, nsad_v - 3
               IF (nvar_vc(i,n).ne.0) THEN
                  nsv = nsv + 1
                  sad_v_c(i,n) = xvariables(nsv)
               END IF
            END DO
            IF (lsplbkp) THEN
               ! breakpoints for v
               DO n = 5, nsad_v
                  nsv = nsv + 1
                  sad_v_s(i,n) = xvariables(nsv)
                  IF (sad_v_s(i,n) .le. sad_v_s(i,n-1)) THEN
                     STOP 'v-spline knots decreasing'
                  END IF
               END DO
            END IF
            IF (nsad_u .gt. 0) THEN
               ! coefficients for u
               DO n = 1, nsad_u - 3
                  IF (nvar_uc(i,n).ne.0) THEN
                     nsv = nsv + 1
                     sad_u_c(i,n) = xvariables(nsv)
                  END IF
               END DO
               IF (lsplbkp) THEN
                  ! breakpoints for u
                  DO n = 5, nsad_u
                     nsv = nsv + 1
                     sad_u_s(i,n) = xvariables(nsv)
                     IF (sad_u_s(i,n) .le. sad_u_s(i,n-1)) THEN
                        STOP 'u-spline knots decreasing'
                     END IF
                  END DO
               END IF
            END IF
         END DO

      ELSE

         ! load fourier series representation
         DO i = 1, nsmid
            modes = 0
            nsv = nsv + 1
            sad_v_c(i,modes) = xvariables(nsv)
            sad_v_s(i,modes) = 0
            DO modes = 1,nsad_v
               nsv = nsv + 1
               sad_v_c(i,modes) = xvariables(nsv)
               nsv = nsv + 1
               sad_v_s(i,modes) = xvariables(nsv)
            END DO

            IF (nsad_u .gt. 0) THEN
               modes = 0
               nsv = nsv + 1
               sad_u_c(i,modes) = xvariables(nsv)
               sad_u_s(i,modes) = 0
               DO modes = 1,nsad_u
                  nsv = nsv + 1
                  sad_u_c(i,modes) = xvariables(nsv)
                  nsv = nsv + 1
                  sad_u_s(i,modes) = xvariables(nsv)
               END DO
            END IF
         END DO

      END IF                 ! IF (lspline)

      nvariables = nsv

      END SUBROUTINE load_saddle_coils
