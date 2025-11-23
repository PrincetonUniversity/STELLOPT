      SUBROUTINE spline_bkp_distance
      USE saddle_coils
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, n
      REAL(rprec), PARAMETER :: vs1 = 1, us1 = 1
      REAL(rprec) :: bkp_tst
!-----------------------------------------------

      bkp_min=1000

      DO i = 1, nsmid

         IF (nsad_v .gt. 5) THEN
            DO n = 4, nsad_v - 1
               bkp_tst = sad_v_s(i,n+1) - sad_v_s(i,n)
               bkp_min = MIN(bkp_min, bkp_tst)
            END DO
            bkp_tst = vs1 - sad_v_s(i,nsad_v)
            bkp_min = MIN(bkp_min, bkp_tst)
         END IF

         IF (nsad_u .gt. 5) THEN
            DO n = 4, nsad_u - 1
               bkp_tst = sad_u_s(i,n+1) - sad_u_s(i,n)
               bkp_min = MIN(bkp_min, bkp_tst)
            END DO
            bkp_tst = us1 - sad_u_s(i,nsad_u)
            bkp_min = MIN(bkp_min, bkp_tst)
         END IF

      END DO

      END SUBROUTINE spline_bkp_distance
