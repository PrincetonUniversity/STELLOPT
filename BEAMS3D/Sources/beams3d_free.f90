!-----------------------------------------------------------------------
!     Module:        beams3d_free
!     Authors:       S. Lazerson (lazerson@pppl.gov), M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          12/15/2014
!     Description:   Deallocate and free all arrays. 
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_free
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE beams3d_runtime
      USE beams3d_grid, ONLY: raxis,phiaxis,zaxis, B_R, B_Z, B_PHI, MODB, &
                                 BR_spl, BZ_spl, BPHI_spl, MODB_spl, &
                                 TE_spl, NE_spl, TI_spl, TE, NE, TI, &
                                 TE_spl_s, NE_spl_s, TI_spl_s, Vp_spl_s,&
                                 S_ARR, S_spl, U_ARR, U_spl, &
                                 POT_ARR, POT_spl, ZEFF_ARR, ZEFF_spl
      USE beams3d_lines, ONLY: R_lines, PHI_lines, Z_lines, vll_lines, &
                               neut_lines, moment_lines, S_lines, U_lines, &
                               PE_lines, PI_lines, shine_through, &
                               ndot_prof, epower_prof, ipower_prof, j_prof,&
                               B_lines
!      USE wall_mod, ONLY: wall_free
      USE EZspline_obj
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier
!-----------------------------------------------------------------------
!     External Functions
!          A00ADF               NAG Detection
!-----------------------------------------------------------------------
!      EXTERNAL A00ADF
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ier = 0
      !IF (lvessel) CALL wall_free(ier)  !Moved to beams3d_follow
      IF (EZspline_allocated(BR_spl))   CALL EZspline_free(BR_spl,ier)
      IF (EZspline_allocated(BZ_spl))   CALL EZspline_free(BZ_spl,ier)
      IF (EZspline_allocated(BPHI_spl)) CALL EZspline_free(BPHI_spl,ier)
      IF (EZspline_allocated(MODB_spl)) CALL EZspline_free(MODB_spl,ier)
      IF (EZspline_allocated(S_spl)) CALL EZspline_free(S_spl,ier)
      IF (EZspline_allocated(U_spl)) CALL EZspline_free(U_spl,ier)
      IF (EZspline_allocated(TE_spl))   CALL EZspline_free(TE_spl,ier)
      IF (EZspline_allocated(NE_spl))   CALL EZspline_free(NE_spl,ier)
      IF (EZspline_allocated(TI_spl))   CALL EZspline_free(TI_spl,ier)
      IF (EZspline_allocated(ZEFF_spl))   CALL EZspline_free(ZEFF_spl,ier)
      IF (EZspline_allocated(POT_spl))   CALL EZspline_free(POT_spl,ier)
      IF (EZspline_allocated(TE_spl_s))   CALL EZspline_free(TE_spl_s,ier)
      IF (EZspline_allocated(NE_spl_s))   CALL EZspline_free(NE_spl_s,ier)
      IF (EZspline_allocated(TI_spl_s))   CALL EZspline_free(TI_spl_S,ier)
      IF (EZspline_allocated(Vp_spl_s))   CALL EZspline_free(Vp_spl_S,ier)
      IF (ALLOCATED(R_lines)) DEALLOCATE(R_lines)
      IF (ALLOCATED(PHI_lines)) DEALLOCATE(PHI_lines)
      IF (ALLOCATED(Z_lines)) DEALLOCATE(Z_lines)
      IF (ALLOCATED(vll_lines)) DEALLOCATE(vll_lines)
      IF (ALLOCATED(neut_lines)) DEALLOCATE(neut_lines)
      IF (ALLOCATED(moment_lines)) DEALLOCATE(moment_lines)
      IF (ALLOCATED(PE_lines)) DEALLOCATE(PE_lines)
      IF (ALLOCATED(PI_lines)) DEALLOCATE(PI_lines)
      IF (ALLOCATED(S_lines)) DEALLOCATE(S_lines)
      IF (ALLOCATED(U_lines)) DEALLOCATE(U_lines)
      IF (ALLOCATED(B_lines)) DEALLOCATE(B_lines)
      IF (ALLOCATED(weight)) DEALLOCATE(weight)
      IF (ALLOCATED(beam)) DEALLOCATE(beam)
      IF (ALLOCATED(raxis)) DEALLOCATE(raxis)
      IF (ALLOCATED(phiaxis)) DEALLOCATE(phiaxis)
      IF (ALLOCATED(zaxis)) DEALLOCATE(zaxis)
      IF (ALLOCATED(B_R)) DEALLOCATE(B_R)
      IF (ALLOCATED(B_PHI)) DEALLOCATE(B_PHI)
      IF (ALLOCATED(B_Z)) DEALLOCATE(B_Z)
      IF (ALLOCATED(MODB)) DEALLOCATE(MODB)
      IF (ALLOCATED(S_ARR)) DEALLOCATE(S_ARR)
      IF (ALLOCATED(U_ARR)) DEALLOCATE(U_ARR)
      IF (ALLOCATED(TE)) DEALLOCATE(TE)
      IF (ALLOCATED(NE)) DEALLOCATE(NE)
      IF (ALLOCATED(TI)) DEALLOCATE(TI)
      IF (ALLOCATED(ZEFF_ARR)) DEALLOCATE(ZEFF_ARR)
      IF (ALLOCATED(POT_ARR)) DEALLOCATE(POT_ARR)
      IF (ALLOCATED(R_start))   DEALLOCATE(R_start)
      IF (ALLOCATED(phi_start)) DEALLOCATE(phi_start)
      IF (ALLOCATED(Z_start))   DEALLOCATE(Z_start)
      IF (ALLOCATED(v_neut))    DEALLOCATE(v_neut)
      IF (ALLOCATED(mass))      DEALLOCATE(mass)
      IF (ALLOCATED(charge))    DEALLOCATE(charge)
      IF (ALLOCATED(mu_start))  DEALLOCATE(mu_start)
      IF (ALLOCATED(Zatom))     DEALLOCATE(Zatom)
      IF (ALLOCATED(t_end))     DEALLOCATE(t_end)
      IF (ALLOCATED(vll_start)) DEALLOCATE(vll_start)
      IF (ALLOCATED(beam))      DEALLOCATE(beam)
      IF (ALLOCATED(weight))    DEALLOCATE(weight)
      IF (ALLOCATED(shine_through))    DEALLOCATE(shine_through)
      IF (ALLOCATED(ndot_prof))    DEALLOCATE(ndot_prof)
      IF (ALLOCATED(epower_prof))    DEALLOCATE(epower_prof)
      IF (ALLOCATED(ipower_prof))    DEALLOCATE(ipower_prof)
      IF (ALLOCATED(j_prof))    DEALLOCATE(j_prof)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_free
