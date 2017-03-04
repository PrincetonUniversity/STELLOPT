      SUBROUTINE chisq_3dbd(sigma, ivar, sigma_rms, sigma_max,
     1           TARGET, target_rms, lvv_tgt, ivar_rms, ivar_max,
     2           num, nopt, iflag, lscreen)
      USE stel_kinds
      USE chisq_mod
      USE optim_params, ONLY: ntor_vv, mpol_vv, rbc_vv, zbs_vv,
     1    vv_dist, vv_dist_rms, vv_dist_max, shapeweight,
     2    phi_lim, r_lim, z_lim, nu_vv, nv_vv
      USE optim, ONLY: bigno, nfp_opt, mnmax_opt, iunit_opt_local
      USE boozer_params, ONLY: rmnc_bdy, zmns_bdy, xm_bdy, xn_bdy
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ivar, ivar_rms, ivar_max, nopt
      INTEGER :: num, iflag
      REAL(rprec) :: sigma, sigma_rms, sigma_max, TARGET, target_rms
      LOGICAL, INTENT(in) :: lscreen, lvv_tgt
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i,m,n, ierr, mnmax_vv
      integer, parameter :: ngrid = 64
      real(rprec) :: dist, drms, lim_dist, adj_sigma
      REAL(rprec), DIMENSION(mpol_vv*(2*ntor_vv+1)+ntor_vv+1) ::
     1                    rmnc_vv, zmns_vv, xm_vv, xn_vv
      REAL(rprec), DIMENSION(mnmax_opt) :: xn_pl
      real(rprec), dimension(:,:), allocatable  :: grid_dist
      LOGICAL :: deviance
C-----------------------------------------------
      iflag = 0
      deviance = .false.
      vv_dist = 1.e30
      
      if (abs(sigma) .ge. bigno .and. abs(sigma_rms) .ge.  bigno
     1    .and. abs(sigma_max) .ge. bigno) return
     
      IF(shapeweight .and.
     1   (.not. ABS(sigma_rms) >=  bigno) ) deviance=.true.

      IF (nopt .gt. 0) THEN
         mnmax_vv = mpol_vv*(2*ntor_vv+1)+ntor_vv+1

         allocate( grid_dist(ngrid, ngrid) )

         DO n = 1, mnmax_opt
            xn_pl(n) = NINT(xn_bdy(n)/nfp_opt)       ! convert back to per period
         END DO

         IF ( mnmax_vv > 1) THEN
            i = 0
            DO m=0, mpol_vv
               DO n=-ntor_vv, ntor_vv
                  IF (m.eq.0 .and. n.lt.0) CYCLE
                  i = i+1
                  xm_vv(i) = m
                  xn_vv(i) = n
                  rmnc_vv(i) = rbc_vv(n,m)
                  zmns_vv(i) = zbs_vv(n,m)
               END DO
            END DO

            IF( i .ne. mnmax_vv) THEN
               IF (lscreen) PRINT *,
     1            'mnmax_vv mis-match in chisq_3dbd!:', mnmax_vv, i
               STOP
            END IF

            ierr=0
            IF(.not.deviance) THEN
              CALL surfsep(mnmax_opt, rmnc_bdy, zmns_bdy, xm_bdy, xn_pl,
     1                     mnmax_vv, rmnc_vv, zmns_vv, xm_vv, xn_vv,
     2                     nfp_opt, 64, vv_dist, vv_dist_rms, ierr)
            ELSE
              CALL weighted_surfsep(rmnc_bdy, zmns_bdy, xm_bdy, xn_bdy,
     1                     nfp_opt, mnmax_opt, vv_dist_rms, vv_dist)
   
            END IF
            IF(ierr .ne. 0) THEN
               IF (lscreen)
     1            PRINT *, 'SurfSep returned IERR = ',ierr, '!!!'
               iflag = -21
               RETURN
            END IF

            IF( lscreen) THEN
               WRITE(6,'(a,3es12.3)')  'VV Distances: sep, rms, max:',
     1                   vv_dist, vv_dist_rms, vv_dist_max
            END IF
         END IF

!
!    Calculate distance to a piecewise-linear limiter set
!
         IF( r_lim(1,1) .ne. 0) THEN
            CALL limsep(mnmax_opt, rmnc_bdy, zmns_bdy, xm_bdy, xn_pl,
     1          size(phi_lim), phi_lim, size(r_lim,1), r_lim, z_lim,
     2          nfp_opt, lim_dist, ierr)
     
            IF(ierr .ne. 0) THEN
               PRINT *, 'LimSep returned IERR = ',ierr, '!!!'
               iflag = -21

               DEALLOCATE(grid_dist)
               RETURN
            END IF

            IF( lscreen) THEN
               WRITE(6,'(a,3es12.3)')
     1             'Min. Limiter Separation:', lim_dist
            END IF
            vv_dist = MIN(vv_dist, lim_dist)
         END IF

         IF( sigma .lt. bigno ) THEN
            num = num + 1
            index_array(num) = ivar
            IF( lvv_tgt) THEN
               chisq_target(num) = MAX(target, vv_dist)
            ELSE
               chisq_target(num) = target
            END IF
            chisq_target(num) = target
            chisq_match(num) = vv_dist
         ENDIF

         IF( sigma_rms .lt. bigno .and. .not. deviance) THEN
            adj_sigma = sigma_rms *
     1                  SQRT(FLOAT(count(grid_dist(:,:)<bigno)))

            DO n=1, ngrid
               DO m=1, ngrid
                  num = num + 1
                  index_array(num) = ivar_rms
                  chisq_target(num) = target_rms
                  chisq_descript(num) = descript(ivar_rms)
                  IF( grid_dist(n,m) < bigno) THEN
                     wegt(num) = adj_sigma
                     chisq_match(num) = grid_dist(n,m)
                  ELSE
                     wegt(num) = bigno
                     chisq_match(num) = target_rms
                  END IF
               END DO
            END DO

         ELSE IF( sigma_rms .lt. bigno .and. deviance) THEN
            num = num + 1
            index_array(num) = ivar_rms
            wegt(num) = sigma_rms
            chisq_target(num) = target_rms
            chisq_match(num) = vv_dist_rms
         END IF

         IF( sigma_max .lt. bigno ) THEN
            num = num + 1
            index_array(num) = ivar_max
            wegt(num) = sigma_max
            chisq_target(num) = 0
            chisq_descript(num) = descript(ivar_max)
            chisq_match(num) = vv_dist_max
         END IF

         deallocate(grid_dist)

      ELSE
         IF( sigma .lt. bigno) THEN
            num = num + 1
            IF (nopt .eq. -2) chisq_descript(num) = descript(ivar)
         END IF
         IF( sigma_rms .lt. bigno) THEN
            IF ( .not. deviance) THEN
               IF (nopt .eq. -2) 
     1             chisq_descript(num+1:num+ngrid*ngrid)=
     2                                              descript(ivar_rms)
               num = num + ngrid*ngrid
            ELSE
               num = num + 1
               IF (nopt .eq. -2) chisq_descript(num)=descript(ivar_rms)
            END IF
         END IF
         IF( sigma_max .lt. bigno ) THEN
            num = num + 1
            IF (nopt .eq. -1) chisq_descript(num) = descript(ivar_max)
         END IF
      END IF

      END SUBROUTINE chisq_3dbd
