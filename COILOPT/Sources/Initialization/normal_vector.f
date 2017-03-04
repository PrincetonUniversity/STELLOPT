      SUBROUTINE normal_vector
      USE stel_constants
      USE boundary
      USE bnorm_mod
      USE mpi_params                                         !mpi stuff
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, ku, kv, mn
      REAL(rprec) :: cosmu, cosnv, sinmu, sinnv, alu, alv,
     1  cosmn1, sinmn1
!-----------------------------------------------
      ALLOCATE (rb(nuv), zb(nuv),
     1  phib(nuv), thetab(nuv), x_p(nuv), y_p(nuv), z_p(nuv),
     2  rb_ph(nuv), rb_th(nuv), zb_ph(nuv), zb_th(nuv),
     3  n_r(nuv), n_phi(nuv), n_z(nuv), d_area(nuv),
     4  theta_d(nuv), phi_d(nuv), bnormal_match(nuv),
     5  b_error(nuv), b_mod(nuv), bmn_error(nuv), bsl_error(nuv),
     6  luv(nuv), stat = i)

      IF( myid .EQ. master ) THEN
         IF (i .ne. 0) STOP 'allocation error in normal_vector'
      END IF

      alu = twopi/nu
      alv = twopi/(nv*nfp)

      z_p = zero
      x_p = zero
      y_p = zero
      d_area = zero

!  as part of calculation the normal vector, the array phib calculated for each
!  of the mnmax boundary points  (nedge  = nuv = nu*nv from vmec or something
!  larger from the assignment statements in the read_wout routine.
!  note that the rb, zb points assume vmec symmetry (m*theta-n*zeta), whereas
!  the winding surface has NESCOIL symmetry (m*theta+n*zeta).

      !IF( myid .EQ. master ) PRINT *, 'nedge = ', nedge
      i = 0
      DO ku = 1, nu
         DO kv = 1, nv
            i = i + 1
            thetab(i) = (ku-1)*alu
            phib(i) = (kv-1)*alv
!  boundary grid coordinates - cylindrical - 
            rb(i)  = 0
            zb(i)  = 0
            rb_th(i)  = 0
            rb_ph(i)  = 0
            zb_th(i)  = 0
            zb_ph(i)  = 0
            DO mn = 1, mnmax
               cosnv = COS(xn_b(mn)*phib(i))
               sinnv = SIN(xn_b(mn)*phib(i))
               cosmu = COS(xm_b(mn)*thetab(i))
               sinmu = SIN(xm_b(mn)*thetab(i))
               cosmn1 = cosmu*cosnv + sinmu*sinnv
               sinmn1 = sinmu*cosnv - cosmu*sinnv
               rb(i) = rb (i) + rmnc_b(mn) * cosmn1
               zb(i) = zb (i) + zmns_b(mn) * sinmn1
!  partial derivatives w.r.t. phi, theta
               rb_th(i) = rb_th(i) - xm_b(mn) * rmnc_b(mn) * sinmn1
               rb_ph(i) = rb_ph(i) + xn_b(mn) * rmnc_b(mn) * sinmn1
               zb_th(i) = zb_th(i) + xm_b(mn) * zmns_b(mn) * cosmn1
               zb_ph(i) = zb_ph(i) - xn_b(mn) * zmns_b(mn) * cosmn1
            END DO
         END DO
      END DO
!  surface normal components
      sum_d_area = 0
      DO i=1,nedge
         n_r(i) = rb(i)*zb_th(i)
         n_phi(i) = rb_th(i)*zb_ph(i) - rb_ph(i)*zb_th(i)
         n_z(i) = -rb(i)*rb_th(i)
      END DO
      DO i=1,nedge
         d_area(i) = SQRT(n_r(i)**2 + n_phi(i)**2 + n_z(i)**2)
         sum_d_area = SUM_d_area + d_area(i)
      END DO
!     PRINT *, sum_d_area
      DO i=1,nedge
         n_r(i) = n_r(i)/d_area(i)
         n_phi(i) = n_phi(i)/d_area(i)
         n_z(i) = n_z(i)/d_area(i)
      END DO
!  boundary grid coordinates - cartesian
      DO i=1,nedge
         x_p(i)=rb(i)*COS(phib(i))
         y_p(i)=rb(i)*SIN(phib(i))
         z_p(i)=zb(i)
      END DO

      END SUBROUTINE normal_vector
