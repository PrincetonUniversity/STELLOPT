!-----------------------------------------------------------------------
!     Function:      faxis
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/9/2011
!     Description:   This subroutine calculates the distance return map
!                    from the starting point of a field line.  Used in
!                    calculating the location of the magnetic axis.
!-----------------------------------------------------------------------
      REAL FUNCTION faxis(x)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_fieldlines
      USE pies_background
      USE pies_realspace
      USE pies_runtime
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Input Parameters
!          x           Initial guess for magnetic axis
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: x
!-----------------------------------------------------------------------
!     Local Variables
!          ierr        Error flag
!          numint      Number of equations
!          q           Polar Coordiantes (xsi,eta)
!          foltol      Tollerance for field line following
!          phi         Starting Point
!          phiend      Ending Point
!          w           Work Array
!          relab       Error calcualtion type
!-----------------------------------------------------------------------
      INTEGER :: ier,numint
      DOUBLE PRECISION    :: pi2, phi0, foltol
      DOUBLE PRECISION    :: q(2)
      DOUBLE PRECISION, ALLOCATABLE :: w(:)
      CHARACTER*1 :: relab
!-----------------------------------------------------------------------
!     External Functions
!          fblin       RHS of differential equation
!          outm        Output function
!          D02CJW      Dummy routine as we do not require root finder
!          D02CJF      ODE Solver
!-----------------------------------------------------------------------
      EXTERNAL fblin
      EXTERNAL D02CJW,D02CJX,D02CJF
!-----------------------------------------------------------------------
!     Begin Function
!-----------------------------------------------------------------------
      pi2 = 8 * ATAN(1._rprec)
      last_line = 0
      q(1)      = x
      q(2)      = 0.0
      foltol    = follow_tol
      ier       = 1
      phi0      = 0.0
      phiend    = pi2
      numint    = 2
      faxis     = 0.0
      ALLOCATE(w(28+21*numint),STAT=ier)
      relab     = 'M'
      CALL D02CJF(phi0,phiend,numint,q,fblin,foltol,relab,D02CJX,D02CJW,w,ier)
      IF (ier /= 0) CALL handle_err(D02CJF_ERR,'faxis',ier)
      faxis = SIGN(SQRT(q(2)*q(2)+(q(1)-x)*(q(1)-x)),q(2))
      !PRINT *,'x0,faxis',x,faxis
      !PRINT *,'x,y',q(1),q(2)
      DEALLOCATE(w)
      RETURN
!-----------------------------------------------------------------------
!     End Function
!-----------------------------------------------------------------------
      END FUNCTION faxis
