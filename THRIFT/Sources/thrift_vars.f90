!-----------------------------------------------------------------------
!     Module:        thrift_vars
!     Authors:       L. van Ham
!     Date:          11/XX/2022
!     Description:   This module contains various variables used in the
!                    thrift code.
!-----------------------------------------------------------------------
MODULE thrift_vars
    !-------------------------------------------------------------------
    !     Libraries
    !-------------------------------------------------------------------
    USE stel_kinds, ONLY: rprec
    USE EZspline_obj
    !-------------------------------------------------------------------
    !     Module Variables
    !          lverb         Logical to control screen output
    !-------------------------------------------------------------------
    IMPLICIT NONE

    LOGICAL :: leccd, lnbcd, lohmic, lscreen_subcodes
    LOGICAL, DIMENSION(:), ALLOCATABLE :: lbooz
    INTEGER :: ntimesteps, nrho, npicard, &
               win_thrift_j, win_thrift_jboot, &
               win_thrift_jplasma, win_thrift_jeccd, win_thrift_jnbcd, &
               win_thrift_johmic, win_thrift_rho, win_thrift_t, &
               win_lbooz, win_thrift_jsource
    REAL(rprec) :: tend, jtol, eq_beta, picard_factor
    REAL(rprec), DIMENSION(:), POINTER :: THRIFT_RHO(:), THRIFT_T(:)
    REAL(rprec), DIMENSION(:,:), POINTER :: THRIFT_J, THRIFT_S11, &
                 THRIFT_S12, THRIFT_S22, THRIFT_JBOOT, THRIFT_JPLASMA, &
                 THRIFT_JECCD, THRIFT_JNBCD, THRIFT_JOHMIC, THRIFT_JSOURCE

    TYPE(EZspline1_r8) :: iota_spl


END MODULE thrift_vars
