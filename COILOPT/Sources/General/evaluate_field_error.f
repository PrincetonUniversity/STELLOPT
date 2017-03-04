      SUBROUTINE evaluate_field_error
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      USE boundary
      USE modular_coils
      USE saddle_coils
      USE tf_coils
      USE vf_coils
      USE bcoils_mod
      USE Vcoilpts
      USE Vwire
      USE biotsavart, coil_temp => single_coil
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, j, k, ks, n, nw, mbw, mtfw
      REAL (dp), DIMENSION(nedge) :: xp, yp, zp, bx, by, bz, 
     1                               sinphi, cosphi
      REAL (dp) :: cur, br, bphi
      REAL (dp) :: b_error_rms, bnorm1
      REAL (dp), DIMENSION(:,:), ALLOCATABLE :: xwire
      REAL (dp), DIMENSION(1:3) :: xpt, bvec
!-----------------------------------------------
      DO n=1, nedge
         xp(n) = x_p(n)
         yp(n) = y_p(n)
         zp(n) = z_p(n)
         sinphi(n) = SIN(phib(n))
         cosphi(n) = COS(phib(n))
      END DO

      bx = 0
      by = 0
      bz = 0

      nw = nwire
      ALLOCATE (xwire(1:3, 1:nw))

      IF (lmodular) THEN
!     compute field due to modulars
         MOD: DO j=1, nmod_coils
            cur = curcon(j)
            DO i=1, nw
               xwire(1,i) = x_mod(i,1,j)
               xwire(2,i) = y_mod(i,1,j)
               xwire(3,i) = z_mod(i,1,j)
            END DO
!     CALL to new Biot-Savart routine (see LIBSTELL MODULE biotsavart)
            CALL bsc_construct(coil_temp, 'floop', '', '',cur,xwire)
            DO n = 1, nedge
               xpt(1) = xp(n);  xpt(2) = yp(n);  xpt(3) = zp(n)
               CALL bsc_b(coil_temp, xpt, bvec)
               bx(n) = bx(n) + bvec(1)
               by(n) = by(n) + bvec(2)
               bz(n) = bz(n) + bvec(3)
            END DO
         END DO MOD
      END IF

      
      IF (lsaddle) THEN
         PRINT *,'in saddle', nsad_coils
!        add field due to saddle coils
         SAD: DO j=1, nsad_coils
               IF (nfils .lt. 4) THEN
!              One - three filament model
                  ks = 1
                  cur = c_sad(j)/nfils
               ELSE
!              Five filament model (no current in central filament 1)
                  ks = 2
                  cur = c_sad(j)/(nfils - 1)
               END IF
               FILS: DO k=ks, nfils
                  DO i=1, nw
                     xwire(1,i) = x_sad(i,j,k)
                     xwire(2,i) = y_sad(i,j,k)
                     xwire(3,i) = z_sad(i,j,k)
                  END DO
                  WRITE(327,*) xwire
                  CALL FLUSH(327)
                  PRINT *,'j,k',j,k, cur, nw, nfils
                  CALL bsc_construct(coil_temp,'floop','','',cur,xwire)
                  PRINT *,'after bsc_construct'
                  DO n = 1, nedge
                     xpt(1) = xp(n);  xpt(2) = yp(n);  xpt(3) = zp(n)
                     CALL bsc_b(coil_temp, xpt, bvec)
                     bx(n) = bx(n) + bvec(1)
                     by(n) = by(n) + bvec(2)
                     bz(n) = bz(n) + bvec(3)
                  END DO
               END DO FILS
         END DO SAD
      END IF

      IF (lvf) THEN
!     add field due to vf coils
         VF: DO j=1, nvf
            cur = cvf(j)
            DO i=1, nw
               xwire(1,i) = x_vf(i,1,j)
               xwire(2,i) = y_vf(i,1,j)
               xwire(3,i) = z_vf(i,1,j)
            END DO
            CALL bsc_construct(coil_temp, 'floop','','',cur,xwire)
            DO n = 1, nedge
               xpt(1) = xp(n);  xpt(2) = yp(n);  xpt(3) = zp(n)
               CALL bsc_b(coil_temp, xpt, bvec)
               bx(n) = bx(n) + bvec(1)
               by(n) = by(n) + bvec(2)
               bz(n) = bz(n) + bvec(3)
            END DO
         END DO VF
      END IF

      DEALLOCATE (xwire)
     
      mtfw = MAX(mtfwire, 2)
      ALLOCATE (xwire(1:3,1:mtfw))
      
      IF (ltfc) THEN
!        add field due to tf coil (1/R)
         TF: DO j=1, mtfcoil
            cur = tfc_cur(j)
            DO i=1, mtfw
               xwire(1,i) = tfc_x(j,i)
               xwire(2,i) = tfc_y(j,i)
               xwire(3,i) = tfc_z(j,i)
            END DO
            CALL bsc_construct(coil_temp, 'floop','','',cur,xwire)
            DO n = 1, nedge
               xpt(1) = xp(n);  xpt(2) = yp(n);  xpt(3) = zp(n)
               CALL bsc_b(coil_temp, xpt, bvec)
               bx(n) = bx(n) + bvec(1)
               by(n) = by(n) + bvec(2)
               bz(n) = bz(n) + bvec(3)
            END DO
         END DO TF
      END IF
      
      DEALLOCATE (xwire)

      IF (lbcoil) THEN
         PRINT *,'got here'
!        add field due to background coils
         BG: DO j=1, mbcoils
            DO n = 1, nedge
               xpt(1) = xp(n);  xpt(2) = yp(n);  xpt(3) = zp(n)
               CALL bsc_b(coil_group(j), xpt, bvec)
               bx(n) = bx(n) + bvec(1)
               by(n) = by(n) + bvec(2)
               bz(n) = bz(n) + bvec(3)
            END DO
!            mbw = MAX(mbwires(j), 2)
!            ALLOCATE (xwire(1:3, mbw))
!            cur = bcoil_cur(j)
!            DO i=1, mbw
!               xwire(1,i) = bcoil_x(j,i)
!               xwire(2,i) = bcoil_y(j,i)
!               xwire(3,i) = bcoil_z(j,i)
!            END DO
!            CALL bsc_construct(coil_temp, 'floop','','',cur,xwire)
!            DO n = 1, nedge
!               xpt(1) = xp(n);  xpt(2) = yp(n);  xpt(3) = zp(n)
!               CALL bsc_b(coil_temp, xpt, bvec)
!               bx(n) = bx(n) + bvec(1)
!               by(n) = by(n) + bvec(2)
!               bz(n) = bz(n) + bvec(3)
!            END DO
!            DEALLOCATE (xwire)
         END DO BG
      END IF

      rbphi_avg = 0
      IF (nedge .NE. nuv) STOP 'nedge != nuv'
      bnorm1 = SQRT(SUM(bnormal_match**2)/nedge)
      IF (bnorm1 .EQ. 0) STOP 'ERROR IN EVALUATE_FIELD_ERROR'

      DO n = 1, nedge
         br =    cosphi(n)*bx(n) + sinphi(n)*by(n)
         bphi = -sinphi(n)*bx(n) + cosphi(n)*by(n)
         b_mod(n) = SQRT (br**2 + bphi**2 + bz(n)**2)
         IF (b_mod(n) .EQ. 0._dp) 
     1      STOP '|B| = 0 in evaluate_field_error!'
         b_error(n) =  br*n_r(n) + bphi*n_phi(n) + bz(n)*n_z(n)
     1              +  bnormal_match(n)
         rbphi_avg = rbphi_avg + d_area(n)*rb(n)*bphi
      END DO

      b_error = b_error/bnorm1
      b_error_rms = SQRT(SUM(b_error(1:nedge)**2)/nedge)
      rbphi_avg = rbphi_avg/sum_d_area

      END SUBROUTINE evaluate_field_error
