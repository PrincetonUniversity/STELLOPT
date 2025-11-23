      SUBROUTINE eval_vf_coils
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE boundary
      USE vf_coils
      USE Vwire
      USE coils
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, k, n, j1, j2
      REAL(rprec) :: phi, dphi, rvfc, zvfc, argt, rvfmax
      REAL(rprec) :: avf, bvf, pvf
!-----------------------------------------------

      dphi = twopi/nwire
      nvf = 2*num_vf

      DO n = 1, num_vf
         j1 = 2*(n - 1) + 1
         j2 = j1 + 1
         zvfc = zc_vf(n)
         rvfc = rc_vf(n)
         rvfmax = rc_vf(n)
         DO i = 1, nwire + 1
            phi = (i-1)*dphi
            IF (nrvf_c .eq. 1) THEN
            ! Ellipse
               avf = rvfc
               bvf = rvfc*rcfc_vf(n,1)
               pvf = rcfs_vf(n,1)
               x_vf(i,1,j1) =  avf*COS(pvf)*COS(phi)
     1                      +  bvf*SIN(pvf)*SIN(phi)
               y_vf(i,1,j1) = -avf*SIN(pvf)*COS(phi)
     1                      +  bvf*COS(pvf)*SIN(phi)
               z_vf(i,1,j1) =  zvfc
               x_vf(i,1,j2) =  avf*COS(-pvf)*COS(-phi)
     1                      +  bvf*SIN(-pvf)*SIN(-phi)
               y_vf(i,1,j2) = -avf*SIN(-pvf)*COS(-phi)
     1                      +  bvf*COS(-pvf)*SIN(-phi)
               z_vf(i,1,j2) = -zvfc
               IF (avf .gt. rvfmax) rvfmax = avf
               IF (bvf .gt. rvfmax) rvfmax = bvf
            ELSE IF (nrvf_c .gt. 1) THEN
            ! Fourier Series
               DO k = 1, nrvf_c
                  argt = k*nfp*phi
                  rvfc = rvfc + rcfc_vf(n,k)*COS(argt)
     1                        + rcfs_vf(n,k)*SIN(argt)
               END DO
               IF (rvfc .gt. rvfmax) rvfmax = rvfc
               x_vf(i,1,j1) = rvfc*COS(phi)
               y_vf(i,1,j1) = rvfc*SIN(phi)
               z_vf(i,1,j1) = zvfc
               x_vf(i,1,j2) = rvfc*COS(-phi)
               y_vf(i,1,j2) = rvfc*SIN(-phi)
               z_vf(i,1,j2) = -zvfc
            ELSE
            ! Circle
               IF (rvfc .gt. rvfmax) rvfmax = rvfc
               x_vf(i,1,j1) = rvfc*COS(phi)
               y_vf(i,1,j1) = rvfc*SIN(phi)
               z_vf(i,1,j1) = zvfc
               x_vf(i,1,j2) = rvfc*COS(-phi)
               y_vf(i,1,j2) = rvfc*SIN(-phi)
               z_vf(i,1,j2) = -zvfc
            END IF
         END DO
         rvf_max(n) = rvfmax
      END DO

!     vf coil currents
      k = 0
      DO i = 1, nvf, 2
         k = k + 1
         cvf(i)   =  cc_vf(k)
         cvf(i+1) = -cc_vf(k)
      END DO
!     norm of vf current vector
      cvf_ssq = 0
      DO i=1, nvf
         cvf_ssq = cvf_ssq + cc_vf(i)**2
      END DO

      END SUBROUTINE eval_vf_coils
