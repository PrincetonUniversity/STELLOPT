      MODULE timer_mod
      USE stel_kinds
      IMPLICIT NONE
      
      REAL(dp) :: time_total, time_init, time_diag_prec,                &
                  time_block_prec, time_toijsp, time_tomnsp,            &
                  time_update_bfield, time_update_force,                &
                  time_update_pres, time_apply_precon, time_update_upperv,&
                  time_factor_blocks, time_update_state,                &
                  time_current, time_init_state, time_generate_blocks,  &
                  time_funci, conj_grad_time
      REAL     :: time_divb=0, time_divj=0, time_bgradp=0, time_bdotj=0
      REAL(dp) :: gmres_time
      REAL(dp) :: get_force_harmonics_time
      REAL(dp) :: comp_diag_elements_time
      REAL(dp) :: compute_hessian_time 
      REAL(dp) :: diag_add_pert_time=0
      REAL(dp) :: block_add_pert_time=0

      REAL(dp) :: evolve_funct_island_time=0
      REAL(dp) :: evolve_restart_file_time=0
      REAL(dp) :: evolve_add_resistive_E_time=0

#if defined(SKS)
      REAL(dp) :: total_time

      REAL(dp) :: init_data_time
      REAL(dp) :: init_metric_elements_time
      REAL(dp) :: test_fourier_time
      REAL(dp) :: init_quantities_time
      REAL(dp) :: init_timers_time
      REAL(dp) :: read_wout_file_time
      REAL(dp) :: init_evolution_time
      REAL(dp) :: Spline_Fourier_Modes_time
      REAL(dp) :: Add_Ghost_Points_time
      REAL(dp) :: Spline_OneD_Array_time
      REAL(dp) :: LoadRZL_VMEC_time

      REAL(dp) :: converge_diagonal_time
      REAL(dp) :: converge_blocks_time
      
      REAL(dp) :: diag_evolve_time
      REAL(dp) :: block_evolve_time

      REAL(dp) :: construct_hessian_time
      REAL(dp) :: asymmetry_check_time
      REAL(dp) :: block_factorization_time
      REAL(dp) :: hessian_funct_island_time
      REAL(dp) :: cv_current_time
      REAL(dp) :: bhtobf_time
      REAL(dp) :: toijsp_time
      REAL(dp) :: tomnsp_time
      REAL(dp) :: to_full_mesh_time 

      REAL(dp) :: gmres_funct_island_time
      REAL(dp) :: gmres_init_dgmres_time
      REAL(dp) :: gmres_wrap_time
      REAL(dp) :: gmresr_time

      REAL(dp) :: drive_dgmres_time
      REAL(dp) :: ParyAx_time
      REAL(dp) :: dcopy_time 
      REAL(dp) :: apply_precond_time
      REAL(dp) :: dgemv_time
      REAL(dp) :: getnlforce_time
      REAL(dp) :: gmres_wrap_allgather_time
      REAL(dp) :: gmres_wrap_allreduce_time
      REAL(dp) :: matvec_funct_island_time

      REAL(dp) :: sendrecv_time


      REAL(dp)     :: time_total_max, time_total_min 
      REAL(dp)     :: construct_hessian_time_max, construct_hessian_time_min
      REAL(dp)     :: asymmetry_check_time_max, asymmetry_check_time_min
      REAL(dp)     :: block_factorization_time_max, block_factorization_time_min
      REAL(dp)     :: hessian_funct_island_time_max, hessian_funct_island_time_min
      REAL(dp)     :: update_upperv_time_max, update_upperv_time_min
      REAL(dp)     :: init_state_time_max, init_state_time_min
      REAL(dp)     :: update_bfield_time_max, update_bfield_time_min
      REAL(dp)     :: update_pres_time_max, update_pres_time_min
      REAL(dp)     :: update_force_time_max, update_force_time_min
      REAL(dp)     :: residual_funct_island_time_max, residual_funct_island_time_min
      REAL(dp)     :: cv_currents_time_max, cv_currents_time_min
      REAL(dp)     :: get_force_harmonics_time_max, get_force_harmonics_time_min 
      REAL(dp)     :: bhtobf_time_max, bhtobf_time_min  
      REAL(dp)     :: tomnsp_time_max, tomnsp_time_min
      REAL(dp)     :: toijsp_time_max, toijsp_time_min
      REAL(dp)     :: to_full_mesh_time_max, to_full_mesh_time_min
      REAL(dp)     :: total_time_max, total_time_min
      REAL(dp)     :: sendrecv_time_max, sendrecv_time_min

      REAL(dp)     :: gmres_time_max, gmres_time_min
      REAL(dp)     :: gmres_wrap_time_max, gmres_wrap_time_min
      REAL(dp)     :: ParyAx_time_max, ParyAx_time_min

#endif

      CONTAINS

      SUBROUTINE init_timers
!
!       Initialize timer counters
!
      time_diag_prec=0; time_block_prec=0
      time_toijsp=0; time_tomnsp=0
      time_update_bfield=0; time_update_force=0; time_update_pres=0
      time_update_upperv=0; time_update_state=0; time_apply_precon=0
      time_factor_blocks=0; time_init_state=0
      time_generate_blocks=0; time_funci=0;
      time_current=0; gmres_time=0; conj_grad_time=0
      get_force_harmonics_time=0
      comp_diag_elements_time=0;  compute_hessian_time=0

      END SUBROUTINE init_timers

#if defined(SKS)
      SUBROUTINE sks_timers
        total_time=0
        construct_hessian_time=0
        asymmetry_check_time=0
        block_factorization_time=0
        hessian_funct_island_time=0
        cv_current_time=0
        bhtobf_time=0
        toijsp_time=0
        tomnsp_time=0
        to_full_mesh_time=0
        gmres_wrap_time=0
        ParyAx_time=0

        sendrecv_time=0

      END SUBROUTINE sks_timers
#endif

      END MODULE timer_mod
