!-----------------------------------------------------------------------
!     Function:      fieldlines_find_axis
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/21/2012
!     Description:   This subroutine finds the magnetic axis given an
!                    initial guess.  This can be a bit confusing so here
!                    is the full explanation.  This is achieved by
!                    following a field line for one field period
!                    (follow_single).  The distance between the head and
!                    the tail is then calculated (faxis_nag).  This
!                    distance is then minimized by this routine and
!                    the minimum coordinates returned.  So in the
!                    calling order we have
!                         fieldlines_find_axis
!                                  |
!                                E04CCA
!                                  |
!                               faxis_nag
!                                  |
!                              follow_single
!                                  |
!                              ODE Solvers
!                                  |
!                                fblin's
!
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_find_axis(r0,z0,phi0)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_grid
      USE fieldlines_runtime
!-----------------------------------------------------------------------
!     Input Variables
!          r0         Initial R position
!          phi0       Initial phi position
!          z0         Initial z position
!          phi1       Final phi position
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(inout) :: r0, z0,phi0
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      INTEGER     :: nfp,i
      REAL(rprec) :: r1,z1,rerr,zerr,alpha,phi2
      REAL(rprec) :: q4(4)
      REAL(rprec) :: detq, traq, l1r, l1i, l2r, l2i
      COMPLEX(8)  :: sqrtval
!-----------------------------------------------------------------------
!     External Functions
!          E04CBF               NAG mimizier
!          E04CBK               NAG Dummy function
!-----------------------------------------------------------------------
      EXTERNAL E04FCF, monit_find_axis, faxis_nag, E04FDZ
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      alpha = 0.75
      nfp = (pi2/phimax)
      phi2 = phi0+pi2/nfp
      delta_phi = phimax-phimin
      !phi2 = phi0+pi2
      IF (lverb) WRITE(6,'(A)') '===========AXIS SEARCH=========='
      IF (lverb) WRITE(6,'(A)') '        R0        Z0       Error     '
      IF (lverb) WRITE(6,'(A,F8.5,F8.5)') '     ',r0,z0
      i=0
      DO
         i = i+1
         r1=r0
         z1=z0
         CALL follow_single(r1,phi0,z1,phi2,q4)
         rerr = r0-r1
         zerr = z0-z1
         IF (lverb) WRITE(6,'(A,F8.5,1X,F8.5,1X,ES11.4)') '     ',r0,z0,sqrt(rerr*rerr+zerr*zerr)
         IF (sqrt(rerr*rerr+zerr*zerr) < (follow_tol*10)) EXIT ! Factor of 10 prevents overworking the algorithm.
         q4(1) = q4(1) - 1
         q4(4) = q4(4) - 1
         r0   = r0 + alpha*(q4(4) * rerr - q4(2) * zerr) / (q4(1)*q4(4)-q4(2)*q4(3))
         z0   = z0 + alpha*(-q4(3) * rerr + q4(1) * zerr) / (q4(1)*q4(4)-q4(2)*q4(3))
         IF (i > 200) EXIT
      END DO
      IF (lverb) WRITE(6,'(A,F8.5,1X,F8.5,1X,ES11.4)') '     ',r0,z0,sqrt(rerr*rerr+zerr*zerr)
      ! Now return eigenvalues of q
      ! From http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
      q4(1) = q4(1) + 1
      q4(4) = q4(4) + 1
      detq = (q4(1)*q4(4)-q4(2)*q4(3))
      traq = q4(1)+q4(4)
      sqrtval = sqrt(CMPLX(0.25*traq*traq-detq))
      l1r   = 0.5*traq + REAL(sqrtval)
      l1i   =          + AIMAG(sqrtval)
      l2r   = 0.5*traq - REAL(sqrtval)
      l2i   =          - AIMAG(sqrtval)
      iota0 = 2.*ATAN2(l1i,l1r)
      IF (lverb) THEN
         WRITE(6,'(A,F8.5,SP,F8.5,A)') '      AXIS: EIGN. VAL 1    =  ',l1r,l1i,'i '
         WRITE(6,'(A,F8.5,SP,F8.5,A)') '            EIGN. VAL 2    =  ',l2r,l2i,'i '
         WRITE(6,'(A,F8.5)')           '      IOTA: ',iota0
      END IF
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fieldlines_find_axis
