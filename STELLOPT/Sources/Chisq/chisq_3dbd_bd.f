C----------PPPL Routine writen by M. Zarnstorff
      SUBROUTINE chisq_3dbd_bd(sigma, ivar, sigma_rms, sigma_max,
     1           TARGET, target_rms, ivar_rms, ivar_max,
     2           num, nopt, iflag, lscreen)
      USE stel_kinds
      USE chisq_mod
      USE optim_params, ONLY: ntor_bd, mpol_bd, rbc_bd, zbs_bd,
     1    bd_dist, bd_dist_rms, bd_dist_max, shapeweight,
     2    phi_lim, r_lim, z_lim, nu_bd, nv_bd
      USE optim, ONLY: bigno, nfp_opt, mnmax_opt, iunit_opt_local
      USE boozer_params, ONLY: rmnc_bdy, zmns_bdy, xm_bdy, xn_bdy
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ivar, ivar_rms, ivar_max, nopt
      INTEGER :: num, iflag
      REAL(rprec) :: sigma, sigma_rms, sigma_max, TARGET, target_rms
      LOGICAL, INTENT(in) :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i,m,n, ierr, mnmax_bd
      INTEGER, PARAMETER :: ngrid = 64
      real(rprec) :: dist, drms, lim_dist, adj_sigma
      REAL(rprec), DIMENSION(mpol_bd*(2*ntor_bd+1)+ntor_bd+1) ::
     1                    rmnc_bd, zmns_bd, xm_bd, xn_bd
      REAL(rprec), DIMENSION(mnmax_opt) :: xn_pl
      real(rprec), dimension(:,:), allocatable  :: grid_dist
C-----------------------------------------------
      iflag = 0
      bd_dist = 1.e30
      
      IF (ABS(sigma) .ge. bigno .and. ABS(sigma_rms) .ge.  bigno
     1    .and. ABS(sigma_max) .ge. bigno) RETURN

      IF (nopt .gt. 0) THEN
         mnmax_bd = mpol_bd*(2*ntor_bd+1)+ntor_bd+1

         ALLOCATE( grid_dist(ngrid, ngrid) )
         
         DO n = 1, mnmax_opt
            xn_pl(n) = NINT(xn_bdy(n)/nfp_opt)       ! convert back to per period
         END DO
         
         IF (mnmax_bd > 1) THEN
            i = 0
            DO m=0, mpol_bd
               DO n=-ntor_bd, ntor_bd
                  IF (m.eq.0 .and. n.lt.0) CYCLE
                  i = i+1
                  xm_bd(i) = m
                  xn_bd(i) = n
                  rmnc_bd(i) = rbc_bd(n,m)
                  zmns_bd(i) = zbs_bd(n,m)
               END DO
            END DO

            IF( i .ne. mnmax_bd) THEN
               IF (lscreen) PRINT *,
     1         'mnmax_bd mis-match in chisq_3dbd!:', mnmax_bd, i
               STOP
            END IF

            ierr=0
            CALL surfsep(mnmax_opt, rmnc_bdy, zmns_bdy, xm_bdy,xn_pl,
     1                mnmax_bd, rmnc_bd, zmns_bd, xm_bd, xn_bd,
     2                nfp_opt, ngrid, bd_dist, bd_dist_rms, ierr)
            
            IF(ierr .ne. 0) THEN
               IF (lscreen)
     1            PRINT *, 'SurfSep returned IERR = ',ierr, '!!!'
               iflag = -21
               RETURN
            END IF

            IF (lscreen) THEN
               WRITE(6,'(a,3es12.3)')  'BD Distances: sep, rms, max:',
     1                   bd_dist, bd_dist_rms, bd_dist_max
            END IF
         END IF
         
         IF( sigma .lt. bigno ) THEN
            num = num + 1
            index_array(num) = ivar
            wegt(num) = sigma
            chisq_target(num) = target
            chisq_match(num) = bd_dist
         ENDIF
         
         IF( sigma_rms .lt. bigno) THEN
            adj_sigma = sigma_rms *
     1                  SQRT(FLOAT(count(grid_dist(:,:)<bigno)))

            DO n=1, ngrid
               DO m=1, ngrid
                  num = num + 1
                  index_array(num) = ivar_rms
                  chisq_target(num) = target_rms
                  IF( grid_dist(n,m) < bigno) THEN
                     wegt(num) = adj_sigma
                     chisq_match(num) = grid_dist(n,m)
                  ELSE
                     wegt(num) = bigno
                     chisq_match(num) = target_rms
                  END IF
               END DO
            END DO
         END IF
   
         IF ( sigma_max .lt. bigno ) THEN  !(PPPL)
            num = num + 1
            index_array(num) = ivar_max
            wegt(num) = sigma_max
            chisq_target(num) = 0
            chisq_match(num) = bd_dist_max
         END IF

         DEALLOCATE(grid_dist)

      ELSE
         IF( sigma .lt. bigno) THEN
            num = num + 1
            IF (nopt .eq. -2) chisq_descript(num) = descript(ivar)
         END IF
         IF( sigma_rms .lt. bigno) THEN
            IF (nopt .eq. -2) THEN
               DO i=1, ngrid*ngrid
                  num = num + 1
                  chisq_descript(num) = descript(ivar_rms)
               END DO
            END IF
         END IF
         IF ( sigma_max .lt. bigno ) THEN  !(PPPL)
            num=num + 1
            IF (nopt .eq. -2) chisq_descript(num) = descript(ivar_max)
         END IF         
      END IF

      END SUBROUTINE chisq_3dbd_bd
