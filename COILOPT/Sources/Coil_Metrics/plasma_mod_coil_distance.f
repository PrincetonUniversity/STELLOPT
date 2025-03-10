      SUBROUTINE plasma_mod_coil_distance
      USE boundary
      USE modular_coils
      USE Vwire
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: ic, n, j
      REAL(dp) :: xpl, ypl, zpl, xcl, ycl, zcl, d_cp
!-----------------------------------------------

      p_d_min = HUGE(p_d_min)

      DO n = 1, nedge
         xpl=rb(n)*COS(phib(n))
         ypl=rb(n)*SIN(phib(n))
         zpl=zb(n)
         DO ic=1,nmod_coils
            DO j = 1,nwire
               xcl=x_mod(j,1,ic)
               ycl=y_mod(j,1,ic)
               zcl=z_mod(j,1,ic)
               d_cp=(xcl-xpl)**2+(ycl-ypl)**2+(zcl-zpl)**2
               p_d_min=MIN(p_d_min,d_cp)
            END DO
         END DO
      END DO

      p_d_min = SQRT(p_d_min)

      END SUBROUTINE plasma_mod_coil_distance
