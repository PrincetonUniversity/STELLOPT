      SUBROUTINE cauchy(rbdy, zbdy, rubdy, zubdy, rlim, zlim, residue,
     1   sep, distmax, ntheta, nlim)
      USE vparams, ONLY: twopi, dp, rprec
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ntheta, nlim
      REAL(rprec), INTENT(in) :: distmax
      REAL(rprec), DIMENSION(ntheta), INTENT(in) ::
     1     rbdy, zbdy, rubdy, zubdy
      REAL(rprec), DIMENSION(nlim) :: rlim, zlim, residue, sep
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero=0, p5=0.5_dp, two=2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: n, i, imin, imax
      REAL(rprec), DIMENSION(ntheta) :: dsq, dsepdu
      REAL(rprec), DIMENSION(ntheta) :: x1u, y1u
      REAL(rprec) :: delu, dmin, delta_d, alpha, gam0
C-----------------------------------------------
c       Check that the points (rlim(i),zlim(i)) are inside boundary surface
c       using Cauchy's theorem in "complex"-plane (for a fixed
c       toroidal plane, nphi=const. It is assumed that rbdy, zbdy are
c       extended around the full interval, 0-2pi, in theta.
c       with rbdy(1) = rbdy(ntheta+1) (i.e., ntheta intervals)
c
c       Because of numerical inaccuracies, instead of testing on
c       res = 0 (outside), res = 1 (inside), we use the test:
c       res >= .5, inside;  res < .5, outside
c**********************************************************************
c
c       LOCAL VARIABLE ARRAYS
c
c       dsq:    Distance squared between limiter point and plasma boundary
c       sep:    Minimum dsq for each limiter point
c    dsepdu:    .5* d(dsq)/d(theta)
c   residue:    Contour integral of 1/(X-rlim)in complex X=(R,Z) plane
c
      delu = twopi/ntheta
      dmin = 1.E-20_dp*distmax

      DO n = 1, nlim
         residue(n) = zero
         x1u = rbdy(:ntheta) - rlim(n)
         y1u = zbdy(:ntheta) - zlim(n)
         dsq(:ntheta) = x1u*x1u + y1u*y1u
         dsepdu(:ntheta) = x1u*rubdy(:ntheta) + y1u*zubdy(:ntheta)
         residue(n) = residue(n) + SUM((x1u*zubdy(:ntheta) -
     1      y1u*rubdy(:ntheta))/(dsq(:ntheta)+dmin))

         residue(n) = residue(n)/ntheta
!
!        Find actual minimum distance from nth limiter point to boundary
!
         sep(n) = distmax
         DO i = 1,ntheta
            IF( dsq(i).le.sep(n) )then
               imin = i
               sep(n) = dsq(i)
            ENDIF
         ENDDO

!        gamu = two*ABS(dsepdu(imin))*delu

         IF (dsepdu(imin) .le. zero) THEN
            imax = 1 + MOD(imin,ntheta)
         ELSE
            imax = imin
            imin = imax - 1
            IF (imin .eq. 0) imin = ntheta
         ENDIF

         delta_d = two*(dsepdu(imax)-dsepdu(imin))
         alpha = delta_d/delu
!        gamu = gamu/delta_d
!        sep(n) = sep(n) - p5*alpha*gamu**2
         IF (alpha .ne. zero)
     1      gam0 = 0.5_dp - (dsq(imax)-dsq(imin))/(alpha*delu**2)
         sep(n) = dsq(imin) - p5*alpha*(gam0*delu)**2
         IF (sep(n) .lt. zero) sep(n) = zero

      END DO

      END SUBROUTINE cauchy
