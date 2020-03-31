      SUBROUTINE DE_preset
      USE stel_kinds
      USE de_mod
      IMPLICIT NONE

      npopsiz = 0
      ngen = 0
      idum = 0
      strategy = 2
      CR_strategy = 0
      f_cross = 0.5_dp
      pcross = 0.3_dp
      parmin = 0.5_dp
      parmax = 2
      ibound = 1
      out_iter = 1
      save_space = .false.

      n_pop = 0
      n_free = 0
      nopt = 0
      IF (ALLOCATED(ui_xc)) DEALLOCATE( ui_xc)

      END SUBROUTINE DE_preset
