      SUBROUTINE chisq_ellipticity(Target_e, sigma, num, nopt)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: bigno, mnmax => mnmax_opt
      USE vparams, ONLY: zero, one, twopi
      USE boozer_params, ONLY: rmnc_bdy, zmns_bdy, xm_bdy, xn_bdy
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nopt
      INTEGER :: num
      REAL(rprec), INTENT(in) :: target_e, sigma
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: nlen = 100
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, m, mn
      REAL(rprec), DIMENSION(nlen) :: theta, z1
      REAL(rprec) :: ellipticity, rmax0, rmin0, zmax0, zmin0
C-----------------------------------------------

      IF (sigma .ge. bigno) RETURN

      num = num + 1

      IF (nopt .gt. 0) THEN

         DO i = 1,nlen
            theta(i) = (i-1)*(twopi/nlen)
         ENDDO
!
!        Limits waist thickness in symmetry plane
!
         rmax0 = SUM(rmnc_bdy(:mnmax))
         rmin0 = SUM(rmnc_bdy(:mnmax) * (-one)**NINT(xm_bdy(:)))
         z1 = zero
         DO mn = 1, mnmax
            m = NINT(xm_bdy(mn))
            z1 = z1 + zmns_bdy(mn)*SIN(m*theta)
         END DO

         zmax0 = MAXVAL(z1(:))
         zmin0 = MINVAL(z1(:))
         ellipticity = ABS(zmax0 - zmin0)/ABS(rmax0 - rmin0)

         index_array(num) = ivar_ellipticity
         wegt(num) = sigma
         IF (wegt(num) .eq. zero) wegt(num) = one
         chisq_target(num) = target_e
         chisq_match(num) = ellipticity
       
      ELSE

         IF (nopt .eq. -2) chisq_descript(num) = 
     1                     descript(ivar_ellipticity)

      END IF

      END SUBROUTINE chisq_ellipticity
