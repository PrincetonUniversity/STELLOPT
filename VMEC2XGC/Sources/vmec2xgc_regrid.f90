!-----------------------------------------------------------------------
!     Subroutine:    vmec2xgc_regrid
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          01/17/2017
!     Description:   This subroutine maps the VMEC wout variables from 
!                    flux to rho space.
!-----------------------------------------------------------------------
      SUBROUTINE vmec2xgc_regrid
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE read_wout_mod
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER     :: ik,mn,mn0,ier
      INTEGER     :: bcs1(2)
      REAL(rprec) :: f0_temp
      REAL(rprec), ALLOCATABLE :: ftemp(:),rho(:),rho_vmec(:),fmn(:)
      TYPE(EZspline1_r8) :: f_spl
      DOUBLE PRECISION, PARAMETER      :: two = 2.0D+00
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ALLOCATE(rho_vmec(ns),rho(ns),ftemp(ns),fmn(ns))
      FORALL(ik = 1:ns) rho_vmec(ik) =  SQRT(REAL(ik-1) / REAL(ns-1))
      FORALL(ik = 1:ns) rho(ik) =  REAL(ik-1) / REAL(ns-1)
      bcs1=(/0,0/)
      CALL EZspline_init(f_spl,ns,bcs1,ier)
      f_spl%isHermite = 1
      f_spl%x1 = rho_vmec
      DO mn = 1, mnmax
         IF (xm(mn) == 0 .and. xn(mn) == 0) mn0 = mn
         ! RMNC
         f0_temp = rmnc(mn,1)
         ftemp(1:ns) = (rmnc(mn,1:ns)-f0_temp)*rho_vmec(1:ns)**(-INT(xm(mn))/two+1)
         ftemp(1) = 0.0
         CALL EZspline_setup(f_spl,ftemp,ier)
         CALL EZspline_interp(f_spl,ns,rho,fmn,ier)
         rmnc(mn,:) = f0_temp+fmn*SQRT(rho)**(INT(xm(mn))-2)
         IF (xm(mn) < 2) rmnc(mn,1) = f0_temp
         ! LMNS
         f0_temp = lmns(mn,1)
         ftemp(1:ns) = (lmns(mn,1:ns)-f0_temp)*rho_vmec(1:ns)**(-INT(xm(mn))/two+1)
         ftemp(1) = 0.0
         CALL EZspline_setup(f_spl,ftemp,ier)
         CALL EZspline_interp(f_spl,ns,rho,fmn,ier)
         lmns(mn,:) = f0_temp+fmn*SQRT(rho)**(INT(xm(mn))-2)
         IF (xm(mn) < 2) lmns(mn,1) = f0_temp
         ! ZMNS
         f0_temp = zmns(mn,1)
         ftemp(1:ns) = (zmns(mn,1:ns)-f0_temp)*rho_vmec(1:ns)**(-INT(xm(mn))/two+1)
         ftemp(1) = 0.0
         CALL EZspline_setup(f_spl,ftemp,ier)
         CALL EZspline_interp(f_spl,ns,rho,fmn,ier)
         zmns(mn,:) = f0_temp+fmn*SQRT(rho)**(INT(xm(mn))-2)
         IF (xm(mn) < 2) zmns(mn,1) = f0_temp
         ! BSUPUMNC
         f0_temp = bsupumnc(mn,1)
         ftemp(1:ns) = (bsupumnc(mn,1:ns)-f0_temp)*rho_vmec(1:ns)**(-INT(xm(mn))/two+1)
         ftemp(1) = 0.0
         CALL EZspline_setup(f_spl,ftemp,ier)
         CALL EZspline_interp(f_spl,ns,rho,fmn,ier)
         bsupumnc(mn,:) = f0_temp+fmn*SQRT(rho)**(INT(xm(mn))-2)
         IF (xm(mn) < 2) bsupumnc(mn,1) = f0_temp
         ! BSUPVMNC
         f0_temp = bsupvmnc(mn,1)
         ftemp(1:ns) = (bsupvmnc(mn,1:ns)-f0_temp)*rho_vmec(1:ns)**(-INT(xm(mn))/two+1)
         ftemp(1) = 0.0
         CALL EZspline_setup(f_spl,ftemp,ier)
         CALL EZspline_interp(f_spl,ns,rho,fmn,ier)
         bsupvmnc(mn,:) = f0_temp+fmn*SQRT(rho)**(INT(xm(mn))-2)
         IF (xm(mn) < 2) bsupvmnc(mn,1) = f0_temp
      END DO
      DEALLOCATE(ftemp,fmn,rho,rho_vmec)
      CALL EZspline_free(f_spl,ier)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE vmec2xgc_regrid

  
