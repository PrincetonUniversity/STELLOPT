      SUBROUTINE spline_it(ndata, xdata, ydata, npts, x, y, i_full)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ndata, npts, i_full
      REAL(rprec), DIMENSION(ndata), INTENT(IN):: xdata, ydata
      REAL(rprec), DIMENSION(npts), INTENT(IN) :: x
      REAL(rprec), DIMENSION(npts), INTENT(OUT):: y
      REAL(rprec), DIMENSION(:), ALLOCATABLE   :: xfull, yfull,
     1                                            dyfull, wk, dy
!-----------------------------------------------
      INTEGER:: iwk, ierr
      LOGICAL:: lspline
!-----------------------------------------------

      ALLOCATE(xfull(ndata+i_full), yfull(ndata+i_full),
     1         dyfull(ndata+i_full), wk(2*(ndata+i_full)),
     2         dy(npts), stat=ierr)
      IF (ierr .ne. 0) STOP 'Allocation error in SPLINE_IT'

      IF (i_full .eq. 1) THEN                                            !! I_FULL= 1, data on HALF mesh
         xfull(1) = -1
         xfull(ndata+1) = 1
         xfull(2:ndata) = 0.5_dp*(xdata(1:ndata-1)+xdata(2:ndata))
         yfull(1) = ydata(1) + (xdata(1) + 1)*
     1    (ydata(2)-ydata(1))/(xdata(2)-xdata(1))
         yfull(2:ndata) = 0.5_dp*(ydata(1:ndata-1)+ydata(2:ndata))
         yfull(ndata+1) = ydata(ndata) + (1 - xdata(ndata))*
     1    (ydata(ndata-1)-ydata(ndata))/(xdata(ndata-1)-xdata(ndata))
      ELSE                                                               !! I_FULL =0, data on FULL mesh;
         xfull(1:ndata) = xdata(1:ndata)
         yfull(1:ndata) = ydata(1:ndata)
      END IF

      lspline = .false.
      wk = 0
      ierr = 0
      iwk = 2*(ndata+i_full)
      dyfull = 0

      CALL pchez(ndata+i_full, xfull, yfull, dyfull, lspline,
     1           wk, iwk, ierr)
      IF(ierr.lt.0) STOP 'LEGENDRE: error in SPLINE_IT'

      CALL pchev(ndata+i_full, xfull, yfull, dyfull,
     1           npts, x, y, dy, ierr)
      IF(ierr.lt.0) STOP 'LEGENDRE: error in EVAL_SPLINE'

      DEALLOCATE(xfull, yfull, dyfull, wk, dy)

      END SUBROUTINE spline_it
