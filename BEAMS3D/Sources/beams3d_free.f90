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
                            vr_lines, vphi_lines, vz_lines, &
                            shine_through, wall_hit_valid, wall_hit_field_valid, &
                            wall_hit_model, wall_hit_face, &
                            wall_hit_fraction, wall_hit_time, wall_hit_r, &
                            wall_hit_phi, wall_hit_z, wall_hit_vll, &
                            wall_hit_moment, wall_hit_b, wall_hit_s, &
                            wall_hit_u, wall_hit_vr, wall_hit_vphi, &
                            wall_hit_vz, wall_hit_energy, time_lines, &
                            B_lines, end_state, shine_port, Gfactor, &
                            ndot_prof, epower_prof, ipower_prof, j_prof,&
                            dense_prof, dist5d_prof, dist5d_fida, &
                            win_ndot, win_epower, win_ipower, win_jprof, &
                            win_dense, win_dist5d, win_dist5d_fida, &
                            win_end_state
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
   IF (EZspline_allocated(TE_spl_s))   CALL EZspline_free(TE_spl_s,ier)
   IF (EZspline_allocated(NE_spl_s))   CALL EZspline_free(NE_spl_s,ier)
   IF (EZspline_allocated(TI_spl_s))   CALL EZspline_free(TI_spl_S,ier)
   IF (EZspline_allocated(Vp_spl_s))   CALL EZspline_free(Vp_spl_S,ier)
   DO i = 1, NION
      IF (EZspline_allocated(NI_spl_s(i)))   CALL EZspline_free(NI_spl_s(i),ier)
   END DO
   IF (ALLOCATED(R_lines)) DEALLOCATE(R_lines)
   IF (ALLOCATED(wall_hit_valid)) DEALLOCATE(wall_hit_valid)
   IF (ALLOCATED(wall_hit_field_valid)) DEALLOCATE(wall_hit_field_valid)
   IF (ALLOCATED(wall_hit_model)) DEALLOCATE(wall_hit_model)
   IF (ALLOCATED(wall_hit_face)) DEALLOCATE(wall_hit_face)
   IF (ALLOCATED(wall_hit_fraction)) DEALLOCATE(wall_hit_fraction)
   IF (ALLOCATED(wall_hit_time)) DEALLOCATE(wall_hit_time)
   IF (ALLOCATED(wall_hit_r)) DEALLOCATE(wall_hit_r)
   IF (ALLOCATED(wall_hit_phi)) DEALLOCATE(wall_hit_phi)
   IF (ALLOCATED(wall_hit_z)) DEALLOCATE(wall_hit_z)
   IF (ALLOCATED(wall_hit_vll)) DEALLOCATE(wall_hit_vll)
   IF (ALLOCATED(wall_hit_moment)) DEALLOCATE(wall_hit_moment)
   IF (ALLOCATED(wall_hit_b)) DEALLOCATE(wall_hit_b)
   IF (ALLOCATED(wall_hit_s)) DEALLOCATE(wall_hit_s)
   IF (ALLOCATED(wall_hit_u)) DEALLOCATE(wall_hit_u)
   IF (ALLOCATED(wall_hit_vr)) DEALLOCATE(wall_hit_vr)
   IF (ALLOCATED(wall_hit_vphi)) DEALLOCATE(wall_hit_vphi)
   IF (ALLOCATED(wall_hit_vz)) DEALLOCATE(wall_hit_vz)
   IF (ALLOCATED(wall_hit_energy)) DEALLOCATE(wall_hit_energy)
   IF (ALLOCATED(time_lines)) DEALLOCATE(time_lines)
   IF (ALLOCATED(PHI_lines)) DEALLOCATE(PHI_lines)
   IF (ALLOCATED(Z_lines)) DEALLOCATE(Z_lines)
   IF (ALLOCATED(vll_lines)) DEALLOCATE(vll_lines)
   IF (ALLOCATED(neut_lines)) DEALLOCATE(neut_lines)
   IF (ALLOCATED(moment_lines)) DEALLOCATE(moment_lines)
   IF (ALLOCATED(S_lines)) DEALLOCATE(S_lines)
   IF (ALLOCATED(U_lines)) DEALLOCATE(U_lines)
   IF (ALLOCATED(B_lines)) DEALLOCATE(B_lines)
   IF (ALLOCATED(vr_lines)) DEALLOCATE(vr_lines)
   IF (ALLOCATED(vphi_lines)) DEALLOCATE(vphi_lines)
   IF (ALLOCATED(vz_lines)) DEALLOCATE(vz_lines)
   IF (ALLOCATED(Gfactor)) DEALLOCATE(Gfactor)
   IF (ALLOCATED(shine_through))    DEALLOCATE(shine_through)
   IF (ALLOCATED(shine_port))    DEALLOCATE(shine_port)
   IF (PRESENT(IN_COMM)) THEN
      IF (ASSOCIATED(req_axis)) CALL mpidealloc(req_axis,win_req_axis)
      IF (ASSOCIATED(zeq_axis)) CALL mpidealloc(zeq_axis,win_zeq_axis)
      IF (ASSOCIATED(raxis))    CALL mpidealloc(raxis,win_raxis)
      IF (ASSOCIATED(phiaxis))  CALL mpidealloc(phiaxis,win_phiaxis)
      IF (ASSOCIATED(zaxis))    CALL mpidealloc(zaxis,win_zaxis)
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
      IF (ASSOCIATED(RHO_ARR))  CALL mpidealloc(RHO_ARR,win_RHO_ARR)
      IF (ASSOCIATED(U_ARR))    CALL mpidealloc(U_ARR,win_U_ARR)
      IF (ASSOCIATED(XRHO_ARR)) CALL mpidealloc(XRHO_ARR,win_XRHO_ARR)
      IF (ASSOCIATED(YRHO_ARR)) CALL mpidealloc(YRHO_ARR,win_YRHO_ARR)
      IF (ASSOCIATED(TE))       CALL mpidealloc(TE,win_TE)
      IF (ASSOCIATED(TI))       CALL mpidealloc(TI,win_TI)
      IF (ASSOCIATED(NE))       CALL mpidealloc(NE,win_NE)
      IF (ASSOCIATED(NI))       CALL mpidealloc(NI,win_NI)
      IF (ASSOCIATED(ZEFF_ARR)) CALL mpidealloc(ZEFF_ARR,win_ZEFF_ARR)
      IF (ASSOCIATED(POT_ARR))  CALL mpidealloc(POT_ARR,win_POT_ARR)
      IF (ASSOCIATED(OMEG_ARR))  CALL mpidealloc(OMEG_ARR,win_OMEG_ARR)
      IF (ASSOCIATED(BR4D))     CALL mpidealloc(BR4D,win_BR4D)
      IF (ASSOCIATED(BPHI4D))   CALL mpidealloc(BPHI4D,win_BPHI4D)
      IF (ASSOCIATED(BZ4D))     CALL mpidealloc(BZ4D,win_BZ4D)
      IF (ASSOCIATED(MODB4D))   CALL mpidealloc(MODB4D,win_MODB4D)
      IF (ASSOCIATED(TE4D))     CALL mpidealloc(TE4D,win_TE4D)
      IF (ASSOCIATED(NE4D))     CALL mpidealloc(NE4D,win_NE4D)
      IF (ASSOCIATED(NI5D))     CALL mpidealloc(NI5D,win_NI5D)
      IF (ASSOCIATED(TI4D))     CALL mpidealloc(TI4D,win_TI4D)
      IF (ASSOCIATED(ZEFF4D))   CALL mpidealloc(ZEFF4D,win_ZEFF4D)
      IF (ASSOCIATED(RHO4D))    CALL mpidealloc(RHO4D,win_RHO4D)
      IF (ASSOCIATED(U4D))      CALL mpidealloc(U4D,win_U4D)
      IF (ASSOCIATED(XRHO4D))   CALL mpidealloc(XRHO4D,win_XRHO4D)
      IF (ASSOCIATED(YRHO4D))   CALL mpidealloc(YRHO4D,win_YRHO4D)
      IF (ASSOCIATED(POT4D))    CALL mpidealloc(POT4D,win_POT4D)
      IF (ASSOCIATED(OMEG4D))    CALL mpidealloc(OMEG4D,win_OMEG4D)
      IF (ASSOCIATED(wall_load))   CALL mpidealloc(wall_load,win_wall_load)
      IF (ASSOCIATED(wall_shine))  CALL mpidealloc(wall_shine,win_wall_shine)
      IF (ASSOCIATED(dist5d_prof)) CALL mpidealloc(dist5d_prof,win_dist5d)
      IF (ASSOCIATED(dist5d_fida)) CALL mpidealloc(dist5d_fida,win_dist5d_fida)
      IF (ASSOCIATED(ndot_prof))    DEALLOCATE(ndot_prof)
      IF (ASSOCIATED(epower_prof))  DEALLOCATE(epower_prof)
      IF (ASSOCIATED(ipower_prof))  DEALLOCATE(ipower_prof)
      IF (ASSOCIATED(j_prof))       DEALLOCATE(j_prof)
      IF (ASSOCIATED(dense_prof))   DEALLOCATE(dense_prof)
      IF (ASSOCIATED(dist5d_fida)) CALL mpidealloc(dist5d_fida,win_dist5d_fida)
      IF (ASSOCIATED(raxis_fida))    CALL mpidealloc(raxis_fida,win_raxis_fida)
      IF (ASSOCIATED(phiaxis_fida))  CALL mpidealloc(phiaxis_fida,win_phiaxis_fida)
      IF (ASSOCIATED(zaxis_fida))    CALL mpidealloc(zaxis_fida,win_zaxis_fida)
      IF (ASSOCIATED(energy_fida))   CALL mpidealloc(energy_fida,win_energy_fida)
      IF (ASSOCIATED(pitch_fida))    CALL mpidealloc(pitch_fida,win_pitch_fida)    
      IF (ASSOCIATED(NEUTRONS_ARR)) CALL mpidealloc(NEUTRONS_ARR,win_NEUTRONS)
      IF (ASSOCIATED(E_NEUTRONS)) CALL mpidealloc(E_NEUTRONS,win_E_NEUTRONS)
      IF (ASSOCIATED(R_start)) CALL mpidealloc(R_start, win_R_start)
      IF (ASSOCIATED(PHI_start)) CALL mpidealloc(PHI_start, win_PHI_start)
      IF (ASSOCIATED(Z_start)) CALL mpidealloc(Z_start, win_Z_start)
      IF (ASSOCIATED(vr_start)) CALL mpidealloc(vr_start, win_vr_start)
      IF (ASSOCIATED(vphi_start)) CALL mpidealloc(vphi_start, win_vphi_start)
      IF (ASSOCIATED(vz_start)) CALL mpidealloc(vz_start, win_vz_start)
      IF (ASSOCIATED(mass)) CALL mpidealloc(mass, win_mass)
      IF (ASSOCIATED(charge)) CALL mpidealloc(charge, win_charge)
      IF (ASSOCIATED(mu_start)) CALL mpidealloc(mu_start, win_mu_start)
      IF (ASSOCIATED(Zatom)) CALL mpidealloc(Zatom, win_Zatom)
      IF (ASSOCIATED(t_end)) CALL mpidealloc(t_end, win_t_end)
      IF (ASSOCIATED(vll_start)) CALL mpidealloc(vll_start, win_vll_start)
      IF (ASSOCIATED(beam)) CALL mpidealloc(beam, win_beam)
      IF (ASSOCIATED(weight)) CALL mpidealloc(weight, win_weight)
      IF (ASSOCIATED(lgc2fo_start)) CALL mpidealloc(lgc2fo_start, win_lgc2fo_start)
      IF (ASSOCIATED(end_state)) CALL mpidealloc(lgc2fo_start, win_end_state)
   ELSE
      IF (ASSOCIATED(req_axis)) DEALLOCATE(req_axis)
      IF (ASSOCIATED(zeq_axis)) DEALLOCATE(zeq_axis)
      IF (ASSOCIATED(raxis))    DEALLOCATE(raxis)
      IF (ASSOCIATED(phiaxis))  DEALLOCATE(phiaxis)
      IF (ASSOCIATED(zaxis))    DEALLOCATE(zaxis)
      IF (ASSOCIATED(raxis_fida))    DEALLOCATE(raxis_fida)
      IF (ASSOCIATED(phiaxis_fida))  DEALLOCATE(phiaxis_fida)
      IF (ASSOCIATED(zaxis_fida))    DEALLOCATE(zaxis_fida)
      IF (ASSOCIATED(energy_fida))   DEALLOCATE(energy_fida)
      IF (ASSOCIATED(pitch_fida))    DEALLOCATE(pitch_fida)
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
      IF (ASSOCIATED(OMEG_ARR))  DEALLOCATE(OMEG_ARR)
      IF (ASSOCIATED(BR4D))     DEALLOCATE(BR4D)
      IF (ASSOCIATED(BPHI4D))   DEALLOCATE(BPHI4D)
      IF (ASSOCIATED(BZ4D))     DEALLOCATE(BZ4D)
      IF (ASSOCIATED(MODB4D))   DEALLOCATE(MODB4D)
      IF (ASSOCIATED(TE4D))     DEALLOCATE(TE4D)
      IF (ASSOCIATED(NE4D))     DEALLOCATE(NE4D)
      IF (ASSOCIATED(NI5D))     DEALLOCATE(NI5D)
      IF (ASSOCIATED(TI4D))     DEALLOCATE(TI4D)
      IF (ASSOCIATED(ZEFF4D))   DEALLOCATE(ZEFF4D)
      IF (ASSOCIATED(RHO4D))      DEALLOCATE(RHO4D)
      IF (ASSOCIATED(U4D))      DEALLOCATE(U4D)
      IF (ASSOCIATED(POT4D))    DEALLOCATE(POT4D)
      IF (ASSOCIATED(OMEG4D))    DEALLOCATE(OMEG4D)
      IF (ASSOCIATED(wall_load))    DEALLOCATE(wall_load)
      IF (ASSOCIATED(wall_shine))    DEALLOCATE(wall_shine)
      IF (ASSOCIATED(ndot_prof))    DEALLOCATE(ndot_prof)
      IF (ASSOCIATED(epower_prof))    DEALLOCATE(epower_prof)
      IF (ASSOCIATED(ipower_prof))    DEALLOCATE(ipower_prof)
      IF (ASSOCIATED(j_prof))    DEALLOCATE(j_prof)
      IF (ASSOCIATED(dense_prof))    DEALLOCATE(dense_prof)
      IF (ASSOCIATED(dist5d_prof))   DEALLOCATE(dist5d_prof)
      IF (ASSOCIATED(dist5d_fida)) DEALLOCATE(dist5d_fida)
      IF (ASSOCIATED(NEUTRONS_ARR)) DEALLOCATE(NEUTRONS_ARR)     
      IF (ASSOCIATED(E_NEUTRONS)) DEALLOCATE(E_NEUTRONS)     
      IF (ASSOCIATED(R_start))   DEALLOCATE(R_start)
      IF (ASSOCIATED(phi_start)) DEALLOCATE(phi_start)
      IF (ASSOCIATED(Z_start))   DEALLOCATE(Z_start)  
      IF (ASSOCIATED(vR_start))   DEALLOCATE(vR_start)
      IF (ASSOCIATED(vphi_start)) DEALLOCATE(vphi_start)
      IF (ASSOCIATED(vZ_start))   DEALLOCATE(vZ_start)
      IF (ASSOCIATED(mass))      DEALLOCATE(mass)
      IF (ASSOCIATED(charge))    DEALLOCATE(charge)
      IF (ASSOCIATED(mu_start))  DEALLOCATE(mu_start)
      IF (ASSOCIATED(Zatom))     DEALLOCATE(Zatom)
      IF (ASSOCIATED(t_end))     DEALLOCATE(t_end)
      IF (ASSOCIATED(vll_start)) DEALLOCATE(vll_start)
      IF (ASSOCIATED(lgc2fo_start)) DEALLOCATE(lgc2fo_start)
      IF (ASSOCIATED(beam))      DEALLOCATE(beam)
      IF (ASSOCIATED(weight))    DEALLOCATE(weight)
      IF (ASSOCIATED(end_state)) DEALLOCATE(end_state)
   ENDIF
   RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
   END SUBROUTINE beams3d_free
