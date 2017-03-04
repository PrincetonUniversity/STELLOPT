      SUBROUTINE mod_coil_distance
      USE boundary, ONLY: nfp
      USE modular_coils
      USE saddle_coils
      USE vf_coils
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
      DO i=1,nmid
        ci_min = HUGE(ci_min)
        DO n = 1,nc
!       Distance to other modular coils
           IF (i .NE. n) THEN
              DO j = 1,nwire
                 x_c0=x_mod(j,1,i)
                 y_c0=y_mod(j,1,i)
                 z_c0=z_mod(j,1,i)
                 DO k = 1,nwire
                    x_c1=x_mod(k,1,n)
                    y_c1=y_mod(k,1,n)
                    z_c1=z_mod(k,1,n)
                    d_cc = (x_c0-x_c1)**2
     1                   + (y_c0-y_c1)**2
     2                   + (z_c0-z_c1)**2
                    ci_min = MIN(ci_min, d_cc)
                 END DO
              END DO
           END IF
        END DO
!       Distance to ALL saddle coils
        IF (lsaddle) THEN
           DO n = 1,ns
              DO j = 1,nwire
                 x_c0=x_mod(j,1,i)
                 y_c0=y_mod(j,1,i)
                 z_c0=z_mod(j,1,i)
                 DO k = 1,nwire
                    x_c1=x_sad(k,n,1)
                    y_c1=y_sad(k,n,1)
                    z_c1=z_sad(k,n,1)
                    d_cc = (x_c0-x_c1)**2
     1                   + (y_c0-y_c1)**2
     2                   + (z_c0-z_c1)**2
                    ci_min = MIN(ci_min, d_cc)
                 END DO
              END DO
           END DO
        END IF
!       Distance to ALL vf coils
        IF (lvf) THEN
           DO n = 1,nvf
              DO j = 1,nwire
                 x_c0=x_mod(j,1,i)
                 y_c0=y_mod(j,1,i)
                 z_c0=z_mod(j,1,i)
                 DO k = 1,nwire
                    x_c1=x_vf(k,1,n)
                    y_c1=y_vf(k,1,n)
                    z_c1=z_vf(k,1,n)
                    d_cc = (x_c0-x_c1)**2
     1                   + (y_c0-y_c1)**2
     2                   + (z_c0-z_c1)**2
                    ci_min = MIN(ci_min, d_cc)
                 END DO
              END DO
           END DO
        END IF
        cc_min(i) = SQRT(ci_min)
      END DO

      END SUBROUTINE mod_coil_distance
