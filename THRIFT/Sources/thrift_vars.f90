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
    !-------------------------------------------------------------------
    !     Module Variables
    !          leccd            Calc Elec. Cyclo. Current Drive
    !          lnbcd            Calc Neutral Beam Current Drive
    !          lohmic           Calc Ohmic Curren Drive
    !          lscreen_subcodes Subcodes screen output flag
    !          lbooz            List of surfaces to do boozer transform
    !          ntimesteps       Number of simulation timesteps
    !          nrho             Number of radial gridpoints
    !          npicard          Max number of Picard Iterations
    !          win_XX           Shared memory windows
    !          tend             T-end
    !          jtol             Picard iteration tollerance
    !          picard_factor    Picard iteration factor
    !          THRIFT_RHO       Radial Grid sqrt(s)
    !          THRIFT_T         Temporal Grid
    !          THRIFT_J         Total current gradient (dI/ds)
    !          THRIFT_SXX       Susceptance matrix
    !          THRIFT_JXX       Current density sources
    !          THRIFT_JPLASMA   Inductive Plasma respsonse    
    !          THRIFT_JSOURCE   Total source current density (dI/ds)
    !-------------------------------------------------------------------
    IMPLICIT NONE

    LOGICAL :: leccd, lnbcd, lohmic, ldiagno, lscreen_subcodes
    LOGICAL, DIMENSION(:), ALLOCATABLE :: lbooz
    INTEGER :: ntimesteps, nrho, npicard, &
               win_thrift_j, win_thrift_jboot, &
               win_thrift_jplasma, win_thrift_jeccd, win_thrift_jnbcd, &
               win_thrift_johmic, win_thrift_rho, win_thrift_t, &
               win_lbooz, win_thrift_jsource
    REAL(rprec) :: tend, jtol, picard_factor
    REAL(rprec), DIMENSION(:), POINTER :: THRIFT_RHO(:), THRIFT_T(:)
    REAL(rprec), DIMENSION(:,:), POINTER :: THRIFT_J, THRIFT_S11, &
                 THRIFT_S12, THRIFT_S22, THRIFT_JBOOT, THRIFT_JPLASMA, &
                 THRIFT_JECCD, THRIFT_JNBCD, THRIFT_JOHMIC, THRIFT_JSOURCE
    ! For TRAVIS
    INTEGER, PARAMETER :: nsys   = 16
    INTEGER :: nra_ecrh, nphi_ecrh
    INTEGER, DIMENSION(nsys)     :: wmode_ecrh
    REAL(rprec), DIMENSION(nsys) :: freq_ecrh, power_ecrh
    REAL(rprec), DIMENSION(nsys,3)     :: antennaposition_ecrh, &
                 targetposition_ecrh,rbeam_ecrh,rfocus_ecrh

END MODULE thrift_vars
