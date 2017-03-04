!-----------------------------------------------------------------------
!     Function:      fxmn_lsqmon
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/9/2011
!     Description:   This subroutine is the monitoring funciton for the
!                    iterative xmn calculation
!-----------------------------------------------------------------------
      SUBROUTINE fxmn_lsqmon(m,n,xc,fvecc,fjac,ldfjac,s,igrade,niter,nf,iw,liw,w,lw)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Input Variables
!          m                  Number of residuals
!          n                  Number of variables
!          xc(n)              Variables
!          fvec(m)            Function Values
!          fjac(ldfjac,n)     Jacobian
!          ldfjac             First dimension of fjac
!          s(n)               Singular Values
!          igrade             Grade of the Jacobian Matrix
!          niter              Number of iterations
!          nf                 Number of function evaluations
!          iw(liw)            Workspace
!          liw                Workspace
!          w                  Workspace
!          lw                 Workspace   
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER           igrade, liw, ldfjac, lw, m, n, nf, niter
      INTEGER           iw(liw)
      DOUBLE PRECISION  fjac(ldfjac,n), fvecc(m), s(n), w(lw), xc(n)

      
      PRINT *,'NITER',niter
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fxmn_lsqmon
