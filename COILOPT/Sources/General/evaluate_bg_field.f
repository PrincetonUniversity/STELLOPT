      SUBROUTINE evaluate_bg_field
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
      REAL (rprec) :: cur, bx, by, bz
      REAL (rprec), ALLOCATABLE, DIMENSION(:,:) :: xwire
      REAL (rprec), DIMENSION(1:3) :: bvec, xpt
!-----------------------------------------------
      xpt(1) = 1.4_dp
      xpt(2) = 0
      xpt(3) = 0
      bx = 0
      by = 0
      bz = 0

!     Field due to background coils
      DO n=1, mbcoils
         CALL bsc_b(coil_group(n),xpt,bvec)
         bx = bvec(1) + bx
         by = bvec(2) + by
         bz = bvec(3) + bz
      END DO
!
!!     Field due to background coils
!

!      DO n=1, mbcoils
!         mbw = MAX (mbwires(n), 2)
!         cur = bcoil_cur(n)
!!        PRINT *, 'bcoil_cur = ', bcoil_cur(n)
!         ALLOCATE (xwire(1:3,1:mbw))
!         DO i=1, mbw
!            xwire(1,i) = bcoil_x(n,i)
!            xwire(2,i) = bcoil_y(n,i)
!            xwire(3,i) = bcoil_z(n,i)
!         END DO
!
!!     CALL Biot-Savart routine (see LIBSTELL MODULE bsc)
!         CALL bsc_construct(coil_temp, 'floop', '', '', cur, xwire)
!         CALL bsc_b(coil_temp, xpt, bvec)
!
!         bx = bvec(1) + bx
!         by = bvec(2) + by
!         bz = bvec(3) + bz
!    
!         DEALLOCATE (xwire)
!      END DO

!     Field due to other coils

      nw = nwire
      ALLOCATE (xwire(1:3,1:nw))

      IF (lmodular) THEN
!        compute field due to modulars
         DO j=1, nmod_coils
            cur = curcon(j)
            DO i=1, nw
               xwire(1,i) = x_mod(i,1,j)
               xwire(2,i) = y_mod(i,1,j)
               xwire(3,i) = z_mod(i,1,j)
            END DO
!     CALL to new Biot-Savart routine (see LIBSTELL MODULE biotsavart)
            CALL bsc_construct(coil_temp,'floop', '', '', cur, xwire)
            CALL bsc_b(coil_temp, xpt, bvec)

            bx = bvec(1) + bx
            by = bvec(2) + by
            bz = bvec(3) + bz
         END DO
      END IF

      IF (lsaddle) THEN
!        add field due to saddle coils
         DO j=1, nsad_coils
            IF (nfils .lt. 4) THEN
!           One or three filament model
               ks = 1
               cur = c_sad(j)/nfils
            ELSE
!           Five filament model (no current in central filament 1)
               ks = 2
               cur = c_sad(j)/(nfils - 1)
            END IF
            DO k=ks, nfils
               DO i=1, nw
                  xwire(1,i) = x_sad(i,j,k)
                  xwire(2,i) = y_sad(i,j,k)
                  xwire(3,i) = z_sad(i,j,k)
               END DO
!     CALL to new Biot-Savart routine (see LIBSTELL MODULE biotsavart)
               CALL bsc_construct(coil_temp,'floop','','',cur,xwire)
               CALL bsc_b(coil_temp, xpt, bvec)

               bx = bvec(1) + bx
               by = bvec(2) + by
               bz = bvec(3) + bz
            END DO
         END DO
      END IF

      IF (lvf) THEN
!        add field due to vf coils
         DO j=1, nvf
            cur = cvf(j)
            DO i=1, nw
               xwire(1,i) = x_vf(i,1,j)
               xwire(2,i) = y_vf(i,1,j)
               xwire(3,i) = z_vf(i,1,j)
            END DO
!     CALL to new Biot-Savart routine (see LIBSTELL MODULE biotsavart)
            CALL bsc_construct(coil_temp,'floop','','', cur, xwire)
            CALL bsc_b(coil_temp, xpt, bvec)

            bx = bvec(1) + bx
            by = bvec(2) + by
            bz = bvec(3) + bz
         END DO
      END IF

      DEALLOCATE (xwire)

      mtfw = MAX(mtfwire-1, 2)
      ALLOCATE (xwire(1:3,1:mtfw))

      IF (ltfc) THEN
!        add field due to tf coil (1/R)
         DO j=1, mtfcoil
            cur = tfc_cur(j)
            DO i=1, mtfwire
               xwire(1,i) = tfc_x(j,i)
               xwire(2,i) = tfc_y(j,i)
               xwire(3,i) = tfc_z(j,i)
            END DO
!     CALL to new Biot-Savart routine (see LIBSTELL MODULE biotsavart)
            CALL bsc_construct(coil_temp,'floop','','', cur, xwire)
            CALL bsc_b(coil_temp, xpt, bvec)

            bx = bvec(1) + bx
            by = bvec(2) + by
            bz = bvec(3) + bz
         END DO
      END IF

      DEALLOCATE (xwire)

      PRINT 1010, xpt(1), xpt(2), xpt(3)
      PRINT 1020, bx, by, bz
 1010 FORMAT('          r [m] = ',1p,3e14.4)
 1020 FORMAT('       B(r) [T] = ',1p,3e14.4)

      END SUBROUTINE evaluate_bg_field
