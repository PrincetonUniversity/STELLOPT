!-----------------------------------------------------------------------
!     Module:        beams3d_free
!     Authors:       S. Lazerson (lazerson@pppl.gov), M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          12/15/2014
!     Description:   Deallocate and free all arrays. 
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_free(IN_COMM)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE beams3d_runtime
      USE beams3d_grid
      USE beams3d_lines, ONLY: R_lines, PHI_lines, Z_lines, vll_lines, &
                               neut_lines, moment_lines, S_lines, U_lines, &
                               shine_through, &
                               B_lines, end_state, shine_port, Gfactor, &
                               ndot_prof, epower_prof, ipower_prof, j_prof,&
                               dense_prof, dist5d_prof, dist5d_fida, &
                               win_ndot, win_epower, win_ipower, win_jprof, &
                               win_dense, win_dist5d, win_dist5d_fida
      USE mpi_sharmem
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier,i
      INTEGER, INTENT(INOUT), OPTIONAL :: IN_COMM
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
      DO i = 1, NION
         IF (EZspline_allocated(NI_spl_s(i)))   CALL EZspline_free(NI_spl_s(i),ier)
      END DO
      IF (ALLOCATED(R_lines)) DEALLOCATE(R_lines)
      IF (ALLOCATED(PHI_lines)) DEALLOCATE(PHI_lines)
      IF (ALLOCATED(Z_lines)) DEALLOCATE(Z_lines)
      IF (ALLOCATED(vll_lines)) DEALLOCATE(vll_lines)
      IF (ALLOCATED(neut_lines)) DEALLOCATE(neut_lines)
      IF (ALLOCATED(moment_lines)) DEALLOCATE(moment_lines)
      IF (ALLOCATED(S_lines)) DEALLOCATE(S_lines)
      IF (ALLOCATED(U_lines)) DEALLOCATE(U_lines)
      IF (ALLOCATED(B_lines)) DEALLOCATE(B_lines)
      IF (ALLOCATED(weight)) DEALLOCATE(weight)
      IF (ALLOCATED(beam)) DEALLOCATE(beam)
      IF (ALLOCATED(end_state)) DEALLOCATE(end_state)
      IF (ALLOCATED(Gfactor)) DEALLOCATE(Gfactor)
      !IF (ALLOCATED(energy_fida))  DEALLOCATE(energy_fida)
      !IF (ALLOCATED(pitch_fida))    DEALLOCATE(pitch_fida)
      IF (PRESENT(IN_COMM)) THEN
         IF (ASSOCIATED(req_axis)) CALL mpidealloc(req_axis,win_req_axis)
         IF (ASSOCIATED(zeq_axis)) CALL mpidealloc(zeq_axis,win_zeq_axis)
         IF (ASSOCIATED(raxis))    CALL mpidealloc(raxis,win_raxis)
         IF (ASSOCIATED(phiaxis))  CALL mpidealloc(phiaxis,win_phiaxis)
         IF (ASSOCIATED(zaxis))    CALL mpidealloc(zaxis,win_zaxis)
         IF (ASSOCIATED(raxis_fida))    CALL mpidealloc(raxis_fida,win_raxis_fida)
         IF (ASSOCIATED(phiaxis_fida))  CALL mpidealloc(phiaxis_fida,win_phiaxis_fida)
         IF (ASSOCIATED(zaxis_fida))    CALL mpidealloc(zaxis_fida,win_zaxis_fida)
         IF (ASSOCIATED(hr))       CALL mpidealloc(hr,win_hr)
         IF (ASSOCIATED(hp))       CALL mpidealloc(hp,win_hp)
         IF (ASSOCIATED(hz))       CALL mpidealloc(hz,win_hz)
         IF (ASSOCIATED(hri))      CALL mpidealloc(hri,win_hri)
         IF (ASSOCIATED(hpi))      CALL mpidealloc(hpi,win_hpi)
         IF (ASSOCIATED(hzi))      CALL mpidealloc(hzi,win_hzi)
         IF (ASSOCIATED(B_R))      CALL mpidealloc(B_R,win_B_R)
         IF (ASSOCIATED(B_PHI))    CALL mpidealloc(B_PHI,win_B_PHI)
         IF (ASSOCIATED(B_Z))      CALL mpidealloc(B_Z,win_B_Z)
         IF (ASSOCIATED(MODB))     CALL mpidealloc(MODB,win_MODB)
         IF (ASSOCIATED(S_ARR))    CALL mpidealloc(S_ARR,win_S_ARR)
         IF (ASSOCIATED(U_ARR))    CALL mpidealloc(U_ARR,win_U_ARR)
         IF (ASSOCIATED(X_ARR))    CALL mpidealloc(S_ARR,win_X_ARR)
         IF (ASSOCIATED(Y_ARR))    CALL mpidealloc(U_ARR,win_Y_ARR)
         IF (ASSOCIATED(TE))       CALL mpidealloc(TE,win_TE)
         IF (ASSOCIATED(TI))       CALL mpidealloc(TI,win_TI)
         IF (ASSOCIATED(NE))       CALL mpidealloc(NE,win_NE)
         IF (ASSOCIATED(NI))       CALL mpidealloc(NI,win_NI)
         IF (ASSOCIATED(ZEFF_ARR)) CALL mpidealloc(ZEFF_ARR,win_ZEFF_ARR)
         IF (ASSOCIATED(POT_ARR))  CALL mpidealloc(POT_ARR,win_POT_ARR)
         IF (ASSOCIATED(BR4D))     CALL mpidealloc(BR4D,win_BR4D)
         IF (ASSOCIATED(BPHI4D))   CALL mpidealloc(BPHI4D,win_BPHI4D)
         IF (ASSOCIATED(BZ4D))     CALL mpidealloc(BZ4D,win_BZ4D)
         IF (ASSOCIATED(MODB4D))   CALL mpidealloc(MODB4D,win_MODB4D)
         IF (ASSOCIATED(TE4D))     CALL mpidealloc(TE4D,win_TE4D)
         IF (ASSOCIATED(NE4D))     CALL mpidealloc(NE4D,win_NE4D)
         IF (ASSOCIATED(NI5D))     CALL mpidealloc(NI5D,win_NI5D)
         IF (ASSOCIATED(TI4D))     CALL mpidealloc(TI4D,win_TI4D)
         IF (ASSOCIATED(ZEFF4D))   CALL mpidealloc(ZEFF4D,win_ZEFF4D)
         IF (ASSOCIATED(S4D))      CALL mpidealloc(S4D,win_S4D)
         IF (ASSOCIATED(U4D))      CALL mpidealloc(U4D,win_U4D)
         IF (ASSOCIATED(X4D))      CALL mpidealloc(X4D,win_Y4D)
         IF (ASSOCIATED(Y4D))      CALL mpidealloc(Y4D,win_X4D)
         IF (ASSOCIATED(POT4D))    CALL mpidealloc(POT4D,win_POT4D)
         IF (ASSOCIATED(wall_load))   CALL mpidealloc(wall_load,win_wall_load)
         IF (ASSOCIATED(wall_shine))  CALL mpidealloc(wall_shine,win_wall_shine)
         !IF (ASSOCIATED(ndot_prof))    CALL mpidealloc(ndot_prof,win_ndot)
         !IF (ASSOCIATED(epower_prof))  CALL mpidealloc(epower_prof,win_epower)
         !IF (ASSOCIATED(ipower_prof))  CALL mpidealloc(ipower_prof,win_ipower)
         !IF (ASSOCIATED(j_prof))       CALL mpidealloc(j_prof,win_jprof)
         !IF (ASSOCIATED(dense_prof))   CALL mpidealloc(dense_prof,win_dense)
         IF (ASSOCIATED(ndot_prof))    DEALLOCATE(ndot_prof)
         IF (ASSOCIATED(epower_prof))    DEALLOCATE(epower_prof)
         IF (ASSOCIATED(ipower_prof))    DEALLOCATE(ipower_prof)
         IF (ASSOCIATED(j_prof))    DEALLOCATE(j_prof)
         IF (ASSOCIATED(dense_prof))    DEALLOCATE(dense_prof)
         IF (ASSOCIATED(dist5d_prof)) CALL mpidealloc(dist5d_prof,win_dist5d)
         IF (ASSOCIATED(dist5d_fida)) CALL mpidealloc(dist5d_fida,win_dist5d_fida)
      ELSE
         IF (ASSOCIATED(req_axis)) DEALLOCATE(req_axis)
         IF (ASSOCIATED(zeq_axis)) DEALLOCATE(zeq_axis)
         IF (ASSOCIATED(raxis))    DEALLOCATE(raxis)
         IF (ASSOCIATED(phiaxis))  DEALLOCATE(phiaxis)
         IF (ASSOCIATED(zaxis))    DEALLOCATE(zaxis)
         IF (ASSOCIATED(hr))       DEALLOCATE(hr)
         IF (ASSOCIATED(hp))       DEALLOCATE(hp)
         IF (ASSOCIATED(hz))       DEALLOCATE(hz)
         IF (ASSOCIATED(hri))      DEALLOCATE(hri)
         IF (ASSOCIATED(hpi))      DEALLOCATE(hpi)
         IF (ASSOCIATED(hzi))      DEALLOCATE(hzi)
         IF (ASSOCIATED(B_R))      DEALLOCATE(B_R)
         IF (ASSOCIATED(B_PHI))    DEALLOCATE(B_PHI)
         IF (ASSOCIATED(B_Z))      DEALLOCATE(B_Z)
         IF (ASSOCIATED(MODB))     DEALLOCATE(MODB)
         IF (ASSOCIATED(S_ARR))    DEALLOCATE(S_ARR)
         IF (ASSOCIATED(U_ARR))    DEALLOCATE(U_ARR)
         IF (ASSOCIATED(TE))       DEALLOCATE(TE)
         IF (ASSOCIATED(TI))       DEALLOCATE(TI)
         IF (ASSOCIATED(NE))       DEALLOCATE(NE)
         IF (ASSOCIATED(NI))       DEALLOCATE(NI)
         IF (ASSOCIATED(ZEFF_ARR)) DEALLOCATE(ZEFF_ARR)
         IF (ASSOCIATED(POT_ARR))  DEALLOCATE(POT_ARR)
         IF (ASSOCIATED(BR4D))     DEALLOCATE(BR4D)
         IF (ASSOCIATED(BPHI4D))   DEALLOCATE(BPHI4D)
         IF (ASSOCIATED(BZ4D))     DEALLOCATE(BZ4D)
         IF (ASSOCIATED(MODB4D))   DEALLOCATE(MODB4D)
         IF (ASSOCIATED(TE4D))     DEALLOCATE(TE4D)
         IF (ASSOCIATED(NE4D))     DEALLOCATE(NE4D)
         IF (ASSOCIATED(NI5D))     DEALLOCATE(NI5D)
         IF (ASSOCIATED(TI4D))     DEALLOCATE(TI4D)
         IF (ASSOCIATED(ZEFF4D))   DEALLOCATE(ZEFF4D)
         IF (ASSOCIATED(S4D))      DEALLOCATE(S4D)
         IF (ASSOCIATED(U4D))      DEALLOCATE(U4D)
         IF (ASSOCIATED(POT4D))    DEALLOCATE(POT4D)
         IF (ASSOCIATED(wall_load))    DEALLOCATE(wall_load)
         IF (ASSOCIATED(wall_shine))    DEALLOCATE(wall_shine)
         IF (ASSOCIATED(ndot_prof))    DEALLOCATE(ndot_prof)
         IF (ASSOCIATED(epower_prof))    DEALLOCATE(epower_prof)
         IF (ASSOCIATED(ipower_prof))    DEALLOCATE(ipower_prof)
         IF (ASSOCIATED(j_prof))    DEALLOCATE(j_prof)
         IF (ASSOCIATED(dense_prof))    DEALLOCATE(dense_prof)
         IF (ASSOCIATED(dist5d_prof))   DEALLOCATE(dist5d_prof)
      ENDIF
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
      IF (ALLOCATED(shine_port))    DEALLOCATE(shine_port)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_free
