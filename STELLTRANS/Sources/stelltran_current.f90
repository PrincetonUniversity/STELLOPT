!-----------------------------------------------------------------------
!     Subroutine:    stelltran_current
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          08/31/2015
!     Description:   Calculate and compare the current densities
!-----------------------------------------------------------------------
      SUBROUTINE stelltran_current(itime)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stelltran_runtime
      USE stelltran_data
      USE stelltran_vars
      USE stelltran_equilutils
      USE parambs, ONLY: rhoar, dibs, irup
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: itime
      INTEGER :: ik,ier,i1,i2
      REAL(rprec) :: sflx,tval, nval, ival, omega_pe, lambda, nu_ei, E
      REAL(rprec), ALLOCATABLE :: sarr(:), farr(:), ssarr(:), sfarr(:)
      REAL(rprec), PARAMETER :: smooth = 0.25
      
      INTEGER, PARAMETER :: dp = selected_real_kind(15, 307)
      INTEGER, PARAMETER :: n_curfit = 4
      REAL(dp), DIMENSION(n_curfit+1) :: cur_poly
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Save the Old values
      johm_old = johm
      jboot_old = jboot
      jrf_old = jrf
      jbeam_old = jbeam
      
      ! Calculate Ohmic current (and nu_star for bootstrap correction)
      !ALLOCATE(nu_star(prof_length))
      DO ik = 2, prof_length
         sflx = te(ik,1)
         tval = te(ik,2)
         nval = ne(ik,2)
         CALL eval_prof_spline(nrad,rho,iotaf,sflx,ival,ier)
         IF (nval <=0) nval = 1E16
         IF (tval >0) THEN
            omega_pe = SQRT(nval*ec*ec/(me*e0))
            lambda   = SQRT(e0*ec*tval/(nval*ec*ec))
            lambda   = nval*lambda*lambda*lambda
            nu_ei    = omega_pe*log(lambda)/(16*pi2*lambda*pi2) ! extra pi2 to get into Hz
            ! https://en.wikipedia.org/wiki/Collisionality
            E        = SQRT(sflx/aspect) ! SQRT inverse aspect ratio.
            nu_star(ik,itime) = nu_ei*R0*SQRT(me/(ec*tval))/(ival*E*E*E)
            E        = Vloop(itime) / (pi2*R0)
            tval     = nval*ec*ec*E/(me*nu_ei)
         ELSE
            tval = 0
            nu_star(ik,itime) = 0
         END IF
         johm(ik,2) = johm(ik,2) - smooth*(johm(ik,2)-tval)
      END DO
      WHERE (nu_star < 0) nu_star = 0
      nu_star(1,itime) = 2*nu_star(2,itime)-nu_star(3,itime)
      johm_sav(:,itime) = johm(:,2)
      
      ! Bootstrap calculation (with nu* correction)
      CALL stelltran_calcboot
      ! Polynomial fitting
      !cur_poly = polyfit(rhoar,dibs,n_curfit)
      !DO ik = 1, irup
      !   dibs(ik) = polyval(cur_poly,rhoar(ik),n_curfit+1)
      !END DO
      ALLOCATE(sarr(irup+2), farr(irup+2))
      ! Smoothing
      ALLOCATE(sfarr(irup))
      CALL smoothg(dibs,irup,0.75_rprec,sfarr)
      dibs = sfarr
      DEALLOCATE(sfarr)
      sarr(1) = 0; sarr(irup+2) = 1;
      sarr(2:irup+1) = rhoar
      farr(2:irup+1) = dibs
      farr(1) = farr(2) + (sarr(1)-sarr(2))*(farr(3)-farr(2))/(sarr(3)-sarr(2))
      farr(irup+2) = farr(irup+1) + (sarr(irup+2)-sarr(irup+1))*(farr(irup+1)-farr(irup))/(sarr(irup+1)-sarr(irup))
      DO ik = 1, prof_length
         sflx = jboot(ik,1)
         CALL eval_prof_spline(irup+2,sarr,farr,sflx,tval,ier)
         tval = tval/(1.0+nu_star(ik,itime))
         jboot(ik,2) = jboot(ik,2) - smooth*(jboot(ik,2)-tval)
      END DO
      DEALLOCATE(sarr,farr)
      jboot_sav(:,itime) = jboot(:,2)
     
      
      !  ECCD calculation
      CALL stelltran_ecrh(itime)
      jecrh_sav(:,itime) = jrf(:,2)
      
      
      !!!!!!!!!!!!!!!BEAM CODE
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE stelltran_current
