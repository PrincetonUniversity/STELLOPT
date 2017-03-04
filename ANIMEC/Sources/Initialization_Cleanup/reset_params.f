      SUBROUTINE reset_params
      USE precon2d, ONLY: ictrl_prec2d
      USE vmec_main, ONLY: iequi, ivac, ftolv, fsqr, fsqz
      USE vsvd, ONLY: pfac, phifac
      USE timer_sub, ONLY: timer
      IMPLICIT NONE

!     2d preconditioner
      ictrl_prec2d = 0

      iequi = 0
      ivac  = -1

      fsqr = 1
      fsqz = 1
      ftolv = fsqr

      pfac   = 1
      phifac = 1
      timer =  0

      END SUBROUTINE reset_params
