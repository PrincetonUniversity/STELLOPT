!-----------------------------------------------------------------------
!     Function:      fxmn
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/9/2011
!     Description:   This subroutine calculates the RHS of iterative
!                    rhomnc calculation
!-----------------------------------------------------------------------
      SUBROUTINE fxmn(n,xc,fc)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_background, ONLY : mnmax,xm,xn,nfp,k,m
      USE pies_magco, ONLY: rhomnc_polar
      USE pies_fieldlines, ONLY : xu_invar,xv_invar, nintw, rholn, k_invar
!-----------------------------------------------------------------------
!     Input Variables
!          n        Number of variables
!          xc       Variables
!          fc       Function Value
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER          :: n
      DOUBLE PRECISION :: fc
      DOUBLE PRECISION :: xc(n)
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      INTEGER          :: i, ik, mn
      REAL(rprec)      :: rho_temp
      REAL(rprec)      :: xu_local(nintw),xv_local(nintw)
      REAL(rprec)      :: rhotemp2(1:nintw)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! First unpack variables
      i = 0
      ! Now transform to rho
      fc = 0.0
      rhotemp2 = 0.0
      DO i = 1, nintw
         DO mn = 1,n
!           commented out to compile other fxmn_lsq
!           rhotemp2(i) = rhotemp2(i) + xc(mn) * &
!                      cos(xm(mn)*xu_invar(i)+ &
!                          xn(mn)*xv_invar(i))
         END DO
      END DO
      !PRINT *,' '
      !PRINT *,'ik',k_invar
      DO i =1, nintw
      !   PRINT *,'rholn',rholn(k_invar,i),'rhotemp2',rhotemp2(i)
         fc = fc + ABS(rholn(k_invar,i) - rhotemp2(i))
      END DO
      !PRINT *,'fc',fc
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fxmn
