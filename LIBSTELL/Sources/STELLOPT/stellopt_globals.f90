!-----------------------------------------------------------------------
!     Module:        stellopt_globals
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          05/31/2024
!     Description:   This module contains the STELLOPT global variables
!                    needed by the input namelist.
!-----------------------------------------------------------------------
      MODULE stellopt_globals
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      
!-----------------------------------------------------------------------
!     Module Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      ! Moved from stellopt_runtime
      INTEGER, PARAMETER :: maxwindsurf=32  
      REAL(rprec), PARAMETER :: bigno = 1.0E+10
      LOGICAL :: lcentered_differences, lkeep_mins, lrefit, lcoil_geom, lno_restart, ltriangulate
      INTEGER :: cr_strategy, npopulation, noptimizers, mode, rho_exp
      INTEGER  ::  nfunc_max
      REAL(rprec)  :: ftol, xtol, gtol, epsfcn, factor, refit_param
      CHARACTER(256)           :: opt_type, axis_init_option
      
      CONTAINS

      INTEGER FUNCTION getmaxwindsurf()
      IMPLICIT NONE
      getmaxwindsurf = maxwindsurf
      END FUNCTION getmaxwindsurf

      REAL(rprec) FUNCTION getbigno()
      IMPLICIT NONE
      getbigno = bigno
      END FUNCTION getbigno

      END MODULE stellopt_globals