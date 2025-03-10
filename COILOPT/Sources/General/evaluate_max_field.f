      SUBROUTINE evaluate_max_field (iunit)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      USE bcoils_mod
      USE modular_coils
      USE Vwire
      USE biotsavart, coil_temp => single_coil
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: iunit, i, j, k, n, nw, ncpts, icl
      REAL (rprec) :: cur
      REAL (rprec) :: bxp, byp, bzp, bpmax, bp
      REAL (rprec) :: u0, v0, du, delta
      REAL (rprec) :: tx, ty, tz, crv
      REAL (rprec) :: nx, ny, nz, x, y, z
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE :: xwire
      REAL (rprec), DIMENSION(1:3) :: xpt, bvec
!-----------------------------------------------
      nw = nwire                   
      ncpts = 2*nwire
      delta = 0.06_dp
      du = 1.0_dp/ncpts

      ALLOCATE (xwire(1:3, 1:nw))
! loop over unique modular coils
      
      MOD: DO i=1, nmid
         icl = i
         bpmax = 0
         PTS: DO j=1, ncpts
            u0 = (j-1)*du
            CALL modular_curve (icl, u0, v0, tx, ty, tz, crv)
            CALL modular_surface (u0, v0, x, y, z, nx, ny, nz)

! determine xp, yp, zp

            xpt(1) = x - delta*nx
            xpt(2) = y - delta*ny
            xpt(3) = z - delta*nz

! compute field at xpt

            bxp = 0
            byp = 0
            bzp = 0
            DO n = 1, nmod_coils
               cur = curcon(n)
               DO k = 1, nw
                  xwire(1,k) = x_mod(k,1,n)
                  xwire(2,k) = y_mod(k,1,n)
                  xwire(3,k) = z_mod(k,1,n)
               END DO                            ! over k
!     CALL to new Biot-Savart routine (see LIBSTELL MODULE biotsavart)
               CALL bsc_construct(coil_temp, 'floop','','',cur,xwire)
               CALL bsc_b(coil_temp, xpt, bvec)

               bxp = bxp + bvec(1)
               byp = byp + bvec(2)
               bzp = bzp + bvec(3)
            END DO                               ! over n

            bp = SQRT (bxp**2 + byp**2 + bzp**2)
            bpmax = MAX (bpmax, bp)
            WRITE (iunit, 1010) x, y, z, xpt(1), xpt(2), xpt(3), bp

! add field due to vf''s

         END DO PTS

         WRITE (iunit, 1000)
         b_max(i) = bpmax

      END DO MOD

      DEALLOCATE (xwire)

 1000 FORMAT ("")
 1010 FORMAT (7f11.4)

      END SUBROUTINE evaluate_max_field
