!-----------------------------------------------------------------------
!     Function:      fxmn_lsq
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/9/2011
!     Description:   This subroutine calculates the RHS of iterative
!                    rhomnc calculation
!-----------------------------------------------------------------------
      SUBROUTINE fxmn_lsq(iflag,m,n,xc,fvec,iw,liw,w,lw)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_background, ONLY : mnmax,xm,xn,nfp,k
      USE pies_magco, ONLY: rhomnc_polar
      USE pies_fieldlines, ONLY : xu_invar,xv_invar, nintw, rholn, k_invar
!-----------------------------------------------------------------------
!     Input Variables
!          iflag    Communication Flag
!          m        Number of residuals
!          n        Number of variables
!          xc(n)    Variables
!          fvec(m)  Function Values
!          iw(liw)  Workspace
!          liw      Workspace
!          w        Workspace
!          lw       Workspace   
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER          :: iflag, m, n, iw(liw), liw, lw
      DOUBLE PRECISION :: xc(n), fvec(m), w(lw)
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      INTEGER          :: i, j, ik, mn
      REAL(rprec)      :: rho_temp
      REAL(rprec)      :: xu_local(nintw),xv_local(nintw)
      REAL(rprec)      :: rhomnc_temp(1:mnmax,0:k)
      REAL(rprec)      :: rhotemp2(0:k,1:nintw)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! First unpack variables
      i = 0
      DO ik = 0, k-1
         DO mn = 1, mnmax
            i = i + 1
            rhomnc_temp(mn,ik) = xc(i)
         END DO
      END DO
      ! Now transform to rho
      fvec = 0.0
      rhotemp2 = 0.0
      DO ik = 1, k -1
         DO i = 0, nintw
            DO mn = 1,mnmax
               rhotemp2(ik,i) = rhotemp2(ik,i) + rhomnc_temp(mn,ik) * &
                      cos(xm(mn)*xu_invar(i,ik)+ &
                          xn(mn)*xv_invar(i))
            END DO
         END DO
      END DO
      i = 0
      DO ik = 0, k-1
         !PRINT *,' '
         !PRINT *,'ik',ik
         DO j = 0, nintw
            i = i + 1
            fvec(i) = ABS(rholn(ik,j) - rhotemp2(ik,j))
            !PRINT *,ik,j,i,'rholn',rholn(ik,j),'rhotemp2',rhotemp2(ik,j),'fvec(i)',fvec(i)
         END DO
      END DO
      PRINT *,'fvec',SUM(fvec,DIM=1)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fxmn_lsq
