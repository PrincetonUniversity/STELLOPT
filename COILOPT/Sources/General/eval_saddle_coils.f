      SUBROUTINE eval_saddle_coils
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      USE boundary, ONLY: nfp
      USE saddle_coils
      USE saddle_surface
      USE Vcoilpts
      USE Vwire
      USE coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, j, k, n, m, nc, np, ierr
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: uc, us
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: vc, vs
      REAL(rprec) :: u,  v,  s, ds, cmax, dlsq, csum, crad
      REAL(rprec), DIMENSION(2,3,5) :: x,  y,  z
      REAL(rprec), DIMENSION(2) :: c
      REAL(rprec) :: r, alf, bet
      LOGICAL :: la, lb
!-----------------------------------------------
      ALLOCATE (uc(0:nsad_u), us(0:nsad_u+4),
     1          vc(0:nsad_v), vs(0:nsad_v+4))

      alf = 1
      bet = 0
      la = .true.
      lb = .false.

      nc = nsad_coils_per_period
      ds = 1.0_dp/nwire
!     loop over unique coils
      DO n=1, nsmid
         DO k=0, nsad_u
            uc(k) = saddle(n)%u_c(k)
            us(k) = saddle(n)%u_s(k)
         END DO
         DO k=0, nsad_v
            vc(k) = saddle(n)%v_c(k)
            vs(k) = saddle(n)%v_s(k)
         END DO
         m = nc + 1 - n
!        loop over number of poloidal points
         cmax = 0
         csum = 0
         DO i=1, nwire+1
            s = (i-1)*ds
            CALL coil_curve( nsad_u, uc, us, nsad_v, vc, vs, la, lb,
     1         lctrlpt, lspline, nfp, numsurf_sad, m_sad, n_sad,
     2         rmn_sad, zmn_sad, nfils, s, alf, bet, deln, delt,
     3         u, v, r, c, x, y, z, ierr )
            cmax = MAX (cmax, c(1))
            csum = csum + c(1)
!           set coil x,y,z points for ALL field periods
            DO j=1, nfp
               np = (j-1)*nc
               u_sad(i,n + np) = u
               v_sad(i,n + np) = v + 1.0*(j-1)
!              stellarator symmetry
               u_sad(nwire + 2 - i,m + np) = 1.0 - u
               v_sad(nwire + 2 - i,m + np) = -v + 1.0*(j-1)
               DO k=1, nfils
                  x_sad(i,n + np,k) = x(1,j,k)
                  y_sad(i,n + np,k) = y(1,j,k)
                  z_sad(i,n + np,k) = z(1,j,k)
!                 stellarator symmetry
                  x_sad(i,m + np,k) = x(2,j,k)
                  y_sad(i,m + np,k) = y(2,j,k)
                  z_sad(i,m + np,k) = z(2,j,k)
               END DO
            END DO
         END DO
         rs_min(n) = 10000.0_dp
         IF(cmax .gt. 0.0_dp) rs_min(n) = 1.0_dp/cmax
         cs_sum(n) = csum*ds
      END DO

!     coil currents
      DO n=1, nsmid
         m = nc + 1 - n
         DO j=1, nfp
            np = (j-1)*nc
            c_sad(n + np) =  saddle(n)%current
            c_sad(m + np) = -saddle(n)%current
         END DO
      END DO

!     coil length
      DO n=1, nsmid
         ymin_sad(n) = 1000
         rmax_sad(n) = 0
         sad_length(n) = 0
         DO i=1, nwire
            dlsq = (x_sad(i+1,n,1) - x_sad(i,n,1))**2
     1           + (y_sad(i+1,n,1) - y_sad(i,n,1))**2
     2           + (z_sad(i+1,n,1) - z_sad(i,n,1))**2
            sad_length(n) = sad_length(n) + SQRT(dlsq)
            IF (ABS(y_sad(i,n,1)) .lt. ymin_sad(n))
     1         ymin_sad(n) = ABS(y_sad(i,n,1))
            crad = SQRT(x_sad(i,n,1)**2+y_sad(i,n,1)**2)
            IF (crad .gt. rmax_sad(n)) rmax_sad(n) = crad
         END DO
      END DO

      DEALLOCATE (uc, us, vc, vs)

      END SUBROUTINE eval_saddle_coils
