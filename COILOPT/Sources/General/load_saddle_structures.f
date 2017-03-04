      SUBROUTINE load_saddle_structures
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE boundary, ONLY: nfp
      USE saddle_coils
      USE tf_coils
      USE Vcoilpts
      USE Vwire
      USE coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, k, modes, nv, nu, ifail
      REAL(rprec), DIMENSION(4) :: fv, fu
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: vc0, vc, vs
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: uc0, uc, us
      REAL(rprec) :: b1_0_0, b1_1_0, b1_2_0, b2_1_0, b2_2_0, b3_2_0
      REAL(rprec) :: bnm0_0_1, bnm0_1_1, bnm0_2_1
      REAL(rprec) :: bnm1_1_1, bnm1_2_1, bnm2_2_1
      REAL(rprec) :: s0, s1
      REAL(rprec) :: au, bu, cu, av, bv, cv
!-----------------------------------------------
      ALLOCATE (uc0(0:nsad_u), uc(0:nsad_u), us(0:nsad_u+4),
     1          vc0(0:nsad_v), vc(0:nsad_v), vs(0:nsad_v+4))

!     load the coil structures with values from variables in
!     optimization

      DO i = 1, nsmid
         modes = 0
         saddle(i)%v_c(modes) = sad_v_c(i,modes)
         sad_v_s(i,modes) = 0
         saddle(i)%v_s(modes) = 0
         DO modes = 1,nsad_v
            saddle(i)%v_c(modes) = sad_v_c(i,modes)
            saddle(i)%v_s(modes) = sad_v_s(i,modes)
         END DO

         IF (nsad_u .gt. 0) THEN
            modes = 0
            saddle(i)%u_c(modes) = sad_u_c(i,modes)
            sad_u_s(i,modes) = 0
            saddle(i)%u_s(modes) = 0
            DO modes = 1,nsad_u
               saddle(i)%u_c(modes) = sad_u_c(i,modes)
               saddle(i)%u_s(modes) = sad_u_s(i,modes)
            END DO
         END IF
      END DO

!     saddle coil currents for each unique current group
      DO i=1, nsmid
         saddle(i)%current = cursad(nsad_group(i))*csad_scl(i)
      END DO

      IF (lspline) THEN

         IF (lctrlpt) THEN
!     set boundary conditions for control point spline representation

            IF (lsmod) THEN

               DO i = 1, nsmid

                  saddle(i)%u_c(nsad_u) = one + saddle(i)%u_c(1)
                  saddle(i)%u_c(nsad_u - 1) = saddle(i)%u_c(nsad_u)
     1               - saddle(i)%u_c(2) + saddle(i)%u_c(1)
                  saddle(i)%u_c(nsad_u - 2) = saddle(i)%u_c(nsad_u)
     1               - saddle(i)%u_c(3) + saddle(i)%u_c(1)
                  saddle(i)%v_c(nsad_v)     = saddle(i)%v_c(1)
                  saddle(i)%v_c(nsad_v - 1) = 2.0_dp*saddle(i)%v_c(1)
     1               - saddle(i)%v_c(2)
                  saddle(i)%v_c(nsad_v - 2) = 2.0_dp*saddle(i)%v_c(1)
     1               - saddle(i)%v_c(3)

                  sad_u_c(i,nsad_u)     = one + saddle(i)%u_c(1)
                  sad_u_c(i,nsad_u - 1) = saddle(i)%u_c(nsad_u)
     1               - saddle(i)%u_c(2) + saddle(i)%u_c(1)
                  sad_u_c(i,nsad_u - 2) = saddle(i)%u_c(nsad_u)
     1               - saddle(i)%u_c(3) + saddle(i)%u_c(1)
                  sad_v_c(i,nsad_v)     = saddle(i)%v_c(1)
                  sad_v_c(i,nsad_v - 1) = 2.0_dp*saddle(i)%v_c(1)
     1               - saddle(i)%v_c(2)
                  sad_v_c(i,nsad_v - 2) = 2.0_dp*saddle(i)%v_c(1)
     1               - saddle(i)%v_c(3)

               END DO

            ELSE

               DO i = 1, nsmid

                  au = saddle(i)%u_c(1)
                  bu = 0.5_dp*(4.0_dp*saddle(i)%u_c(2)
     1                       - 3.0_dp*saddle(i)%u_c(1)
     2                       - saddle(i)%u_c(3))
                  cu = 0.5_dp*(saddle(i)%u_c(3)
     1                       - 2.0_dp*saddle(i)%u_c(2)
     2                       + saddle(i)%u_c(1))

                  av = saddle(i)%v_c(1)
                  bv = 0.5_dp*(4.0_dp*saddle(i)%v_c(2)
     1                       - 3.0_dp*saddle(i)%v_c(1)
     2                       - saddle(i)%v_c(3))
                  cv = 0.5_dp*(saddle(i)%v_c(3)
     1                       - 2.0_dp*saddle(i)%v_c(2)
     2                       + saddle(i)%v_c(1))

                  saddle(i)%u_c(nsad_u) = au
                  saddle(i)%v_c(nsad_v) = av
                  saddle(i)%u_c(nsad_u - 1) = au - bu + cu
                  saddle(i)%u_c(nsad_u - 2) = au - 2.0_dp*bu + 4.0_dp*cu
                  saddle(i)%v_c(nsad_v - 1) = av - bv + cv
                  saddle(i)%v_c(nsad_v - 2) = av - 2.0_dp*bv + 4.0_dp*cv

                  sad_u_c(i,nsad_u)     = saddle(i)%u_c(nsad_u)
                  sad_v_c(i,nsad_v)     = saddle(i)%v_c(nsad_v)
                  sad_u_c(i,nsad_u - 1) = saddle(i)%u_c(nsad_u - 1)
                  sad_v_c(i,nsad_v - 1) = saddle(i)%v_c(nsad_v - 1)
                  sad_u_c(i,nsad_u - 2) = saddle(i)%u_c(nsad_u - 2)
                  sad_v_c(i,nsad_v - 2) = saddle(i)%v_c(nsad_v - 2)

               END DO

            END IF

         ELSE

!     set boundary conditions for spline representation of v

            nv = nsad_v

            DO i = 1, nsmid
               ! breakpoints
               DO k = 5, nv
                  vs(k) = saddle(i)%v_s(k)
               END DO
               ! boundary conditions on breakpoints
               DO k = 1, 4
                  vs(k) = 0
                  vs(nv + k) = 1
               END DO

               ! boundary conditions on coefficients
               s0 = 0
               s1 = 1
               vc0 = 0

               vc0(1) = 1
               CALL speval (nv, vs, vc0, s0, fv, ifail)
               b1_0_0 = fv(1)
               b1_1_0 = fv(2)
               b1_2_0 = fv(3)
               vc0(1) = 0

               vc0(2) = 1
               CALL speval (nv, vs, vc0, s0, fv, ifail)
               b2_1_0 = fv(2)
               b2_2_0 = fv(3)
               vc0(2) = 0

               vc0(3) = 1
               CALL speval (nv, vs, vc0, s0, fv, ifail)
               b3_2_0 = fv(3)
               vc0(3) = 0

               vc0(nv) = 1
               CALL speval (nv, vs, vc0, s1, fv, ifail)
               bnm0_0_1 = fv(1)
               bnm0_1_1 = fv(2)
               bnm0_2_1 = fv(3)
               vc0(nv) = 0

               vc0(nv-1) = 1
               CALL speval (nv, vs, vc0, s1, fv, ifail)
               bnm1_1_1 = fv(2)
               bnm1_2_1 = fv(3)
               vc0(nv-1) = 0

               vc0(nv-2) = 1
               CALL speval (nv, vs, vc0, s1, fv, ifail)
               bnm2_2_1 = fv(3)
               vc0(nv-2) = 0

               vc = 0
               DO k = 1, 3
                  vc(k) = saddle(i)%v_c(k)
               END DO
               vc(nv) = vc(1)*(b1_0_0/bnm0_0_1)
               vc(nv-1) = (vc(1)*b1_1_0 + vc(2)*b2_1_0
     1                  - vc(nv)*bnm0_1_1) /bnm1_1_1
               vc(nv-2) = (vc(1)*b1_2_0 + vc(2)*b2_2_0 + vc(3)*b3_2_0
     1                  - vc(nv-1)*bnm1_2_1 - vc(nv)*bnm0_2_1)/bnm2_2_1
               DO k = nv-2, nv
                  saddle(i)%v_c(k) = vc(k)
                  sad_v_c(i,k) = vc(k)
               END DO

            END DO     ! i=1,nsmid

!     set boundary conditions for spline representation of u

            nu = nsad_u

            DO i = 1, nsmid
               ! breakpoints
               DO k = 5, nu
                  us(k) = saddle(i)%u_s(k)
               END DO
               ! boundary conditions on breakpoints
               DO k = 1, 4
                  us(k) = 0
                  us(nu + k) = 1
               END DO

               ! boundary conditions on coefficients
               s0 = 0
               s1 = 1
               uc0 = 0

               uc0(1) = 1
               CALL speval (nu, us, uc0, s0, fu, ifail)
               b1_0_0 = fu(1)
               b1_1_0 = fu(2)
               b1_2_0 = fu(3)
               uc0(1) = 0

               uc0(2) = 1
               CALL speval (nu, us, uc0, s0, fu, ifail)
               b2_1_0 = fu(2)
               b2_2_0 = fu(3)
               uc0(2) = 0

               uc0(3) = 1
               CALL speval (nu, us, uc0, s0, fu, ifail)
               b3_2_0 = fu(3)
               uc0(3) = 0

               uc0(nu) = 1
               CALL speval (nu, us, uc0, s1, fu, ifail)
               bnm0_0_1 = fu(1)
               bnm0_1_1 = fu(2)
               bnm0_2_1 = fu(3)
               uc0(nu) = 0

               uc0(nu-1) = 1
               CALL speval (nu, us, uc0, s1, fu, ifail)
               bnm1_1_1 = fu(2)
               bnm1_2_1 = fu(3)
               uc0(nu-1) = 0

               uc0(nu-2) = 1
               CALL speval (nu, us, uc0, s1, fu, ifail)
               bnm2_2_1 = fu(3)
               uc0(nu-2) = 0

               uc = 0
               DO k = 1, 3
                  uc(k) = saddle(i)%u_c(k)
               END DO
               uc(nu) = uc(1)*(b1_0_0/bnm0_0_1)
               uc(nu-1) = (uc(1)*b1_1_0 + uc(2)*b2_1_0
     1                  - uc(nu)*bnm0_1_1)/bnm1_1_1
               uc(nu-2) = (uc(1)*b1_2_0 + uc(2)*b2_2_0 + uc(3)*b3_2_0
     1                  - uc(nu-1)*bnm1_2_1 - uc(nu)*bnm0_2_1)/bnm2_2_1
               DO k = nu-2, nu
                  saddle(i)%u_c(k) = uc(k)
                  sad_u_c(i,k) = uc(k)
               END DO

            END DO     ! i=1,nsmid

         END IF        ! IF (lctrlpt)

      END IF           ! IF (lspline)

      DEALLOCATE (vc0, vc, vs, uc0, uc, us)

      END SUBROUTINE load_saddle_structures
