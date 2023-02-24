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
    !          jtol             Picard iteration tolerance
    !          picard_factor    Picard iteration factor
    !          THRIFT_RHO       Radial Grid sqrt(s)
    !          THRIFT_T         Temporal Grid
    !          THRIFT_J         Total current gradient
    !          THRIFT_SXX       Susceptance matrix
    !          THRIFT_JXX       Current density sources
    !          THRIFT_JPLASMA   Inductive Plasma respsonse    
    !          THRIFT_JSOURCE   Total source current density
    !          THRIFT_UGRID     u = S11*iota+S12
    !          THRIFT_UEDGE     u at rho = 1
    !          THRIFT_I         Total enclosed current
    !          THRIFT_IPLASMA   Total enclosed induced current  
    !          THRIFT_IXXXXX    Total enclosed bootstrap/CD currents         
    !-------------------------------------------------------------------
    IMPLICIT NONE

    LOGICAL :: leccd, lnbcd, lohmic, ldiagno, lscreen_subcodes
    LOGICAL, DIMENSION(:), ALLOCATABLE :: lbooz
    REAL(rprec), DIMENSION(2) :: THRIFT_UEDGE, THRIFT_PHIEDGE
    INTEGER :: ntimesteps, nrho, npicard, &
               win_thrift_j, win_thrift_jboot, &
               win_thrift_jplasma, win_thrift_jeccd, win_thrift_jnbcd, &
               win_thrift_johmic, win_thrift_rho, win_thrift_t, &
               win_lbooz, win_thrift_jsource, win_thrift_ugrid, &
               win_thrift_vp, win_thrift_bav, win_thrift_bsqav, &
               win_thrift_s11, win_thrift_aminor, win_thrift_rmajor
    REAL(rprec) :: tstart, tend, jtol, picard_factor, THRIFT_I, THRIFT_IPLASMA,&
                 THRIFT_IBOOT, THRIFT_IECCD, THRIFT_INBCD, THRIFT_IOHMIC
    REAL(rprec), DIMENSION(:), POINTER :: THRIFT_RHO(:), THRIFT_T(:)
    REAL(rprec), DIMENSION(:,:), POINTER :: THRIFT_J, THRIFT_S11, &
                 THRIFT_S12, THRIFT_S22, THRIFT_JBOOT, THRIFT_JPLASMA, &
                 THRIFT_JECCD, THRIFT_JNBCD, THRIFT_JOHMIC, THRIFT_JSOURCE, &
                 THRIFT_UGRID, THRIFT_VP, THRIFT_BAV, THRIFT_BSQAV, &
                 THRIFT_S11, THRIFT_AMINOR, THRIFT_RMAJOR

                 
    ! For TRAVIS
    INTEGER, PARAMETER :: nsys   = 16
    INTEGER :: nra_ecrh, nphi_ecrh
    INTEGER, DIMENSION(nsys)     :: wmode_ecrh
    REAL(rprec), DIMENSION(nsys) :: freq_ecrh, power_ecrh
    REAL(rprec), DIMENSION(nsys,3)     :: antennaposition_ecrh, &
                 targetposition_ecrh,rbeam_ecrh,rfocus_ecrh

END MODULE thrift_vars
