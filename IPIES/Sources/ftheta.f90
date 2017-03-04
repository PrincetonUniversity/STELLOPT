!-----------------------------------------------------------------------
!     Function:      ftheta
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/05/2012
!     Description:   The routine returns the value and derivative of
!                    the function:
!                             f=(theta-x*phi)^2
!                    where x is the rotational transform.  Used in
!                    the calculation of iota (x=iota).
!                     
!-----------------------------------------------------------------------
      SUBROUTINE ftheta(x,f,fprime)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_fieldlines
      USE pies_background
      USE pies_realspace
      USE pies_runtime
!-----------------------------------------------------------------------
!     Input Variables
!          x        Variable
!          f        f(x)
!          fprime   df/dx
!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION :: x, f, fprime
!-----------------------------------------------------------------------
!     Local Variables
!          ierr        Error flag
!-----------------------------------------------------------------------
      INTEGER :: ier, mn
      DOUBLE PRECISION :: df(0:nintw), dfp(0:nintw)
      
      EXTERNAL :: SDOT
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      df(0) = 0.0
      PRINT *,'iotabar',thetaln(isurf,1:10)/philn(isurf,1:10)
      stop
      df(1:nintw) = thetaln(isurf,1:nintw)-x*philn(isurf,1:nintw)/nfp
      dfp = 1.0
      f = SUM(df(0:nintw)*df(0:nintw))
      !f = SDOT(nintw,df,1,df,1)
      dfp(1:nintw) = -2.0*philn(isurf,1:nintw)*dfp(1:nintw)/nfp
      !fprime = SDOT(nintw,df,1,dfp,1)
      fprime = SUM(df(0:nintw)*dfp(0:nintw))
      PRINT *,'     iota=',x,'  f=',f,'  fprime=',fprime
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE ftheta
