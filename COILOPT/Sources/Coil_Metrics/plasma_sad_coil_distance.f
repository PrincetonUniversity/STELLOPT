      SUBROUTINE plasma_sad_coil_distance
      USE boundary
      USE saddle_coils
      USE Vwire
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, n, j, ntot
      REAL(rprec) :: xpl, ypl, zpl, xcl, ycl, zcl, d_cp
      REAL(rprec) :: xsad(nwire*nsad_coils), 
     1               ysad(nwire*nsad_coils), zsad(nwire*nsad_coils)
!-----------------------------------------------

      p_s_min = HUGE(p_s_min)
      ntot = 0
      DO i = 1, nsad_coils
         DO j = 1,nwire
            ntot = ntot + 1
            xsad(ntot) = x_sad(j,i,1)
            ysad(ntot) = y_sad(j,i,1)
            zsad(ntot) = z_sad(j,i,1)
         END DO
      END DO


      DO n = 1, nedge
         xpl=rb(n)*COS(phib(n))
         ypl=rb(n)*SIN(phib(n))
         zpl=zb(n)
         DO j = 1, ntot
            xcl = xsad(j) 
            ycl = ysad(j)
            zcl = zsad(j)
            d_cp = (xcl-xpl)**2+(ycl-ypl)**2+(zcl-zpl)**2
            p_s_min = MIN(p_s_min,d_cp)
         END DO
      END DO

      p_s_min = SQRT(p_s_min)

      END SUBROUTINE plasma_sad_coil_distance
