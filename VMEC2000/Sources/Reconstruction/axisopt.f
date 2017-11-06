      SUBROUTINE axisopt(fsq, r00, iresidue, ivac)
      USE vsvd
      USE vparams, ONLY: zero, one, nthreed
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: iresidue, ivac
      REAL(rprec) :: fsq, r00
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: smax = 0.998_dp
      REAL(rprec), PARAMETER :: smin = 0.985_dp
      CHARACTER(LEN=60), PARAMETER :: optbegin =
     1   'Begin variation of Raxis to MINimize total RMS error'
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: delstep, dedrmax, factor, delerr,
     1   dedr, rstepx1
      REAL(rprec), SAVE :: delstep_old, errmax, errold, rstepx,
     1   raxold, scale
C-----------------------------------------------


      IF (iresidue.lt.1 .or. errsvd*1.e6_dp.lt.one .or.
     1    fsq.gt.fturnon_axis .or. ivac .le.2) RETURN
!
!     MOVE R-AXIS BASED ON dR/dt = (-dEsvd/dR)
!     LIMIT MAXIMUM RSTEPX TO RSTEPX0
!     TRY TO FIND ZERO-CROSSING IN dEsvd/dR (ESTIMATED NUMERICALLY)
!

      IF (iresidue .eq. 1) THEN                    !First time through
         iresidue = 2
         raxold = r00
         errold = errsvd
         errmax = zero
         rstepx = rstepx0
         scale = smax
         IF (iopt_raxis .gt. 0) THEN
            WRITE (*, 115) optbegin
            WRITE (nthreed, 115) optbegin
         ENDIF
      ELSE
         delerr = errsvd - errold                !delta E-svd
         delstep = r00 - raxold                  !delta R-axis
         IF (delerr.ne.zero .and. ABS(delstep).gt.1.e-3_dp*rstepx0) THEN
            dedr = delerr/delstep
            errmax = MAX(errmax,errsvd)
            dedrmax = 2*errmax/rwidth
            rstepx1 = MIN(one,ABS(dedr)/dedrmax)*rsfac*rstepx0
            factor = SIGN(one,(-dedr))        !Move in -dE/dR direction
            rstepx = rstepx1*factor
            scale = smax
            IF (delstep*delstep_old .le. zero) scale = smin
            delstep_old = delstep
            raxold = r00
            errold = errsvd
         ENDIF
      ENDIF
      rsfac = scale*rsfac
c-5/1/96 raxmse = raxmse + rstepx
      raxmse = raxold + rstepx
  115 FORMAT(2x,a)

      END SUBROUTINE axisopt
