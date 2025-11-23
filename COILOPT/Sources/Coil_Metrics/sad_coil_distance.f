      SUBROUTINE sad_coil_distance
      USE stel_constants
      USE boundary
      USE saddle_coils
      USE modular_coils
      USE vf_coils
      USE Vwire
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, j, k, l, m, n, nmc, nsc, js, ks, kcount,
     1           kcjl, kckm
      INTEGER, DIMENSION(ncdim) :: nsv, jsv, ksv
      REAL(rprec) :: x_c0, y_c0, z_c0, x_c1, y_c1, z_c1
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: 
     1               x_c, y_c, z_c
      REAL(rprec) :: zpl, xpls, ypls, zpls, cpd
      REAL(rprec) :: xpl(nedge, nfp), ypl(nedge, nfp)
      REAL(rprec) :: dphib, d_cc, ci_min, cin_min, d_cp, d_cxp
!-----------------------------------------------

      dphib = twopi/nfp
      nsc = nsad_coils_per_period*nfp
      nmc = nmod_coils_per_period*nfp

      ALLOCATE (x_c(nwire*nfils, nsc), y_c(nwire*nfils, nsc),
     1          z_c(nwire*nfils, nsc), stat = j)
      IF (j .ne. 0) STOP 'Allocation error in sad_coil_distance'

      DO i = 1, nsc
         kcount = 0
         DO l = 1, nfils
            DO j = 1,nwire
               kcount = kcount+1
               x_c(kcount,i) = x_sad(j,i,1) + 2*(x_sad(j,i,l)
     1                       - x_sad(j,i,1))
               y_c(kcount,i) = y_sad(j,i,1) + 2*(y_sad(j,i,l)
     1                       - y_sad(j,i,1))
               z_c(kcount,i) = z_sad(j,i,1) + 2*(z_sad(j,i,l)
     1                       - z_sad(j,i,1))
            END DO
         END DO
      END DO

!     For each unique coil i ...

      sc_dmin = 0

      BASE_COIL: DO i=1,nsmid
        ci_min = HUGE(ci_min)

        TEST_COIL: DO n = 1,nsc
!          Find the min distance to all other coils n .ne. i
           IF (n .eq. i) CYCLE
           cin_min = HUGE(cin_min)
           kcjl = 0
           COIL_I: DO kcjl = 1, nwire*nfils
              x_c0 = x_c(kcjl,i)
              y_c0 = y_c(kcjl,i)
              z_c0 = z_c(kcjl,i)
!             Find the nearest point k, coil n, to point j on coil i
                 COIL_N: DO kckm = 1, nwire*nfils
                    x_c1 = x_c(kckm,n)
                    y_c1 = y_c(kckm,n)
                    z_c1 = z_c(kckm,n)
                     d_cc = (x_c0-x_c1)**2
     1                    + (y_c0-y_c1)**2
     2                    + (z_c0-z_c1)**2
                    IF (d_cc .lt. cin_min) cin_min = d_cc
                    IF (d_cc .lt. ci_min) THEN
                       ci_min = d_cc
!                      Save indices of points j, k and coil n
                       ns = n
                       js = 1 + MOD(kcjl-1,nwire)        !j
                       ks = 1 + MOD(kckm-1,nwire)        !k
                    END IF
                END DO COIL_N
            END DO COIL_I
            sc_dmin(i, n) = SQRT(cin_min)
         END DO TEST_COIL


!       Distance to all modular coils
        IF (lmodular) THEN
           DO n = 1,nmc
              DO j = 1,nwire
                 x_c0=x_sad(j,i,1)
                 y_c0=y_sad(j,i,1)
                 z_c0=z_sad(j,i,1)
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
           END DO
        END IF
!       Distance to all vf coils
        IF (lvf) THEN
           DO n = 1,nvf
              DO j = 1,nwire
                 x_c0=x_sad(j,i,1)
                 y_c0=y_sad(j,i,1)
                 z_c0=z_sad(j,i,1)
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
        sc_min(i) = SQRT(ci_min)
        nsv(i) = ns
        jsv(i) = js
        ksv(i) = ks
      END DO BASE_COIL

      DEALLOCATE (x_c, y_c, z_c)

      IF (nedge .gt. 0) THEN
!     Find size of min coil-coil vector crossed with min coil-plasma vector
         DO k = 1,nfp
            DO m = 1, nedge
               xpl(m,k) = rb(m)*COS(phib(m) + (k-1)*dphib)
               ypl(m,k) = rb(m)*SIN(phib(m) + (k-1)*dphib)
            END DO
         END DO
         COIL_I1: DO i=1,nsmid
            ns = nsv(i)
            js = jsv(i)
            ks = ksv(i)
            x_c0=x_sad(js,i,1)
            y_c0=y_sad(js,i,1)
            z_c0=z_sad(js,i,1)
            x_c1=x_sad(ks,ns,1)
            y_c1=y_sad(ks,ns,1)
            z_c1=z_sad(ks,ns,1)
            cpd = 1000
!           Min distance to from coil i, point js to plasma
            FP: DO k = 1,nfp
               PLASMA: DO m = 1, nedge
                  zpl=zb(m)
                  d_cp = (x_c0-xpl(m,k))**2
     1                 + (y_c0-ypl(m,k))**2
     2                 + (z_c0-zpl)**2
                  IF (d_cp .lt. cpd) THEN
                     cpd = d_cp
                     xpls = xpl(m,k)
                     ypls = ypl(m,k)
                     zpls = zpl
                  END IF
               END DO PLASMA
            END DO FP
            d_cxp = ((y_c0 - y_c1)*(zpls - z_c0)
     1            -  (ypls - y_c0)*(z_c0 - z_c1))**2
     2            + ((z_c0 - z_c1)*(xpls - x_c0)
     3            -  (zpls - z_c0)*(x_c0 - x_c1))**2
     4            + ((x_c0 - x_c1)*(ypls - y_c0)
     5            -  (xpls - x_c0)*(y_c0 - y_c1))**2
            scxp_min(i) = SQRT(d_cxp/((xpls - x_c0)**2
     1                              + (ypls - y_c0)**2
     2                              + (zpls - z_c0)**2))
         END DO COIL_I1
      END IF

      END SUBROUTINE sad_coil_distance
