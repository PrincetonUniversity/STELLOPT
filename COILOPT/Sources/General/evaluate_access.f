      SUBROUTINE evaluate_access
      USE boundary, ONLY: nfp
      USE bcoils_mod
      USE modular_coils
      USE saddle_coils
      USE Vwire
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, j, k, n, nc, ns
      REAL(rprec) :: x_c0, y_c0, z_c0, x_c1, y_c1, z_c1
      REAL(rprec) :: d_cc, ci_min
!-----------------------------------------------

      nc = nmod_coils_per_period*nfp
      ns = nsad_coils_per_period*nfp
      DO n=1,n_access
         ci_min = 1000
         DO i = 1,n_access_pts
            x_c0=x_access(n,i)
            y_c0=y_access(n,i)
            z_c0=z_access(n,i)
!           Distance to ALL modular coils
            IF (lmodular) THEN
               DO j = 1,nc
                  DO k = 1,nwire
                     x_c1=x_mod(k,1,j)
                     y_c1=y_mod(k,1,j)
                     z_c1=z_mod(k,1,j)
                     d_cc = (x_c0-x_c1)**2
     1                    + (y_c0-y_c1)**2
     2                    + (z_c0-z_c1)**2
                     ci_min = MIN(ci_min, d_cc)
                  END DO
               END DO
            END IF
!           Distance to ALL saddle coils
            IF (lsaddle) THEN
               DO j = 1,ns
                  DO k = 1,nwire
                     x_c1=x_sad(k,j,1)
                     y_c1=y_sad(k,j,1)
                     z_c1=z_sad(k,j,1)
                     d_cc = (x_c0-x_c1)**2
     1                    + (y_c0-y_c1)**2
     2                    + (z_c0-z_c1)**2
                     ci_min = MIN(ci_min, d_cc)
                  END DO
               END DO
            END IF
         END DO
         acc_min(n) = SQRT(ci_min)
      END DO

      END SUBROUTINE evaluate_access
