!-----------------------------------------------------------------------
!     Subroutine:    stelltran_updateeq
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          08/22/2015
!     Description:   This routine updates the equilibrium
!-----------------------------------------------------------------------
      SUBROUTINE stelltran_updateeq(itime)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stelltran_runtime
      USE stelltran_data
      USE stelltran_vars
      USE stelltran_equilutils
!-----------------------------------------------------------------------
!     Local Variables
!        ier          Error flag
!        ik           Indexing variable
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: itime
      INTEGER ::  ier, ik
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! First setup internal grid
      IF (itime == 1) THEN
         FORALL (ik = 1:prof_length) ne(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) te(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) ti(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) zeff(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) johm(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) jboot(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) jrf(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) jbeam(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) Qe(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) Xe(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) Qi(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) Xi(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) Vp(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) rho_coord(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) gradrho(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) gradrhosq(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) S_pe(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) S_pi(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) Ge(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) De(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) Gi(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) Di(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) Er(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) pe_rf(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) pi_col(ik,1) = REAL(ik -1)/REAL(prof_length-1)
         FORALL (ik = 1:prof_length) ni(ik,1) = REAL(ik -1)/REAL(prof_length-1)
      END IF
	
      SELECT CASE(TRIM(run_type))
         CASE('analysis')
         CASE('predictive')
            IF (itime == 1) RETURN
            ! Update profiles according to transport equations
            ! Density
            !   d/dt(nV') - d/dr(|grad(r)|^2 V' D dn/dr) = V' (Sgas + Sbeam - Lc-x,recomb)
            ! Energy
            !   d/dt(T) ...
            CALL stelltran_bal_update_asked(itime)
      END SELECT
            ! Put everything on the uniform internal grid
            DO ik = 1, prof_length
               CALL eval_prof_spline(nne,ne_s,ne_f(:,itime),ne(ik,1),ne(ik,2),ier)
            END DO
            DO ik = 1, prof_length
               CALL eval_prof_spline(nte,te_s,te_f(:,itime),te(ik,1),te(ik,2),ier)
            END DO
            DO ik = 1, prof_length
               CALL eval_prof_spline(nti,ti_s,ti_f(:,itime),ti(ik,1),ti(ik,2),ier)
            END DO
            DO ik = 1, prof_length
               CALL eval_prof_spline(nzeff,zeff_s,zeff_f(:,itime),zeff(ik,1),zeff(ik,2),ier)
            END DO 

      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE stelltran_updateeq
