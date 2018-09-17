!----------------------
! Module BFGS_params
! Contains parameters and variables used by the BFGS finite difference
! omptimizer
!
! JC Schmitt
! 2018
!
!-----------------------
  MODULE BFGS_PARAMS

! Other modules that used by BFGS
  USE mpi_params
  USE safe_open_mod
  USE stel_kinds
  ! JCS to do: verify this list of flags.
  USE fdjac_mod, ONLY: FDJAC2_MP_QUEUE, FLAG_CLEANUP, &
                       FLAG_CLEANUP_BFGS, FLAG_SINGLETASK, &
                       flip, jac_order, h_order, jac_err, jac_index
  ! flag_singletask = -1, flag_cleanup = -100, flag_cleanup_lev = -101
  ! USE jacfcn_mod, ONLY: fjac_curr, write_jacobian, jac_analytic

  IMPLICIT NONE
! BFGS Variables

!    enable_flip     Enables the 'flipping' of the components of the
!                    vector in the finite difference section
!                    Default: .true.
!    enable_tr_red   Enables the alogrithm to reduce or contract the
!                    'trust region' which is used for finite differences
!                    NOT YET IMPLEMENTED

  LOGICAL :: enable_flip = .true.
  LOGICAL :: enable_tr_red = .false.


  ! Logical flag used to control debugging statements in the BFGS subroutines
  ! .true. turns the debbuging on.
  ! .false. turns the debugging off.
  LOGICAL :: DEBUG_BFGS = .true.
  !LOGICAL :: DEBUG_BFGS = .false.

  ! Exit flags for the Armijo backtracking line search
  INTEGER :: BT_EXIT_FLAG
  INTEGER, PARAMETER :: BT_EXIT_NORMAL = 0
  INTEGER, PARAMETER :: BT_EXIT_STEPTOOSMALL = -1

  END MODULE BFGS_PARAMS

