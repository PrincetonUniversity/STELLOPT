      SUBROUTINE init_saddle_coils (nvariables, xvariables, nfp)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE saddle_coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, n, nsv, nc, modes, nfp
      INTEGER :: nvariables
      REAL(rprec) :: xvariables(*)
!-----------------------------------------------

      nc = nsad_coils_per_period
      nvariables = 0
      nsad_coeffs = 0

      nsad_coils = nc * nfp
      nsmid = nc/2
      nsodd = MOD(nc,2)
      nsad_unique_coils = nsmid

      IF (nsad_coils .le. 0) RETURN

!     initialize the variables to values of unique coil parameters
!     and count the number of variables

      nsv = 0

      IF (lspline) THEN

         ! initialize bspline series representation
         DO i = 1, nsmid
            ! coefficients for v
            DO n = 1, nsad_v - 3
               IF (nvar_vc(i,n).ne.0) THEN
                  nsv = nsv + 1
                  xvariables(nsv) = sad_v_c(i,n)
               END IF
            END DO
            IF (lsplbkp) THEN
               ! breakpoints for v
               DO n = 5, nsad_v
                  nsv = nsv + 1
                  xvariables(nsv) = sad_v_s(i,n)
               END DO
            END IF
            IF (nsad_u .gt. 0) THEN
               ! coefficients for u
               DO n = 1, nsad_u - 3
                  IF (nvar_uc(i,n).ne.0) THEN
                     nsv = nsv + 1
                     xvariables(nsv) = sad_u_c(i,n)
                  END IF
               END DO
               IF (lsplbkp) THEN
                  ! breakpoints for u
                  DO n = 5, nsad_u
                     nsv = nsv + 1
                     xvariables(nsv) = sad_u_s(i,n)
                  END DO
               END IF
            END IF
         END DO

      ELSE

         ! initialize fourier series representation
         DO i = 1, nsmid
            modes = 0
            nsv = nsv + 1
            xvariables(nsv) = sad_v_c(i,modes)
            DO modes = 1,nsad_v
               nsv = nsv + 1
               xvariables(nsv) = sad_v_c(i,modes)
               nsv = nsv + 1
               xvariables(nsv) = sad_v_s(i,modes)
            END DO

            modes = 0
            IF (nsad_u .gt. 0) THEN
               nsv = nsv + 1
               xvariables(nsv) = sad_u_c(i,modes)
               DO modes = 1,nsad_u
                  nsv = nsv + 1
                  xvariables(nsv) = sad_u_c(i,modes)
                  nsv = nsv + 1
                  xvariables(nsv) = sad_u_s(i,modes)
               END DO
            END IF
         END DO

      END IF                 ! IF (lspline)

      nsad_coeffs = nsv

      IF (lsadshape) THEN
         nvariables = nsv
      ELSE
         nvariables = 0
      END IF

      END SUBROUTINE init_saddle_coils
