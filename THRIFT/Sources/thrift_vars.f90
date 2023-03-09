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

    LOGICAL :: leccd, lnbcd, lohmic, ldiagno, lscreen_subcodes, lverbj
    LOGICAL, DIMENSION(:), ALLOCATABLE :: lbooz
    INTEGER :: ntimesteps, nrho, npicard, &
             win_thrift_j, win_thrift_jplasma, win_thrift_jboot, &
             win_thrift_jeccd, win_thrift_jnbcd, win_thrift_johmic, &
             win_thrift_jsource, &
             win_thrift_i, win_thrift_iplasma, win_thrift_iboot, &
             win_thrift_ieccd, win_thrift_inbcd, win_thrift_iohmic, &
             win_thrift_isource, &
             win_thrift_rho, win_thrift_t, win_thrift_rhofull,&
             win_lbooz, win_thrift_ugrid, &
             win_thrift_vp, win_thrift_bav, win_thrift_bsqav, win_thrift_phiedge &
             win_thrift_s11, win_thrift_aminor, win_thrift_rmajor, &
             win_thrift_coeff_a, win_thrift_coeff_b, win_thrift_coeff_c, win_thrift_coeff_d,&
             win_thrift_coeff_bp, win_thrift_coeff_cp, win_thrift_coeff_dp,&
             win_thrift_alpha1, win_thrift_alpha2, win_thrift_alpha3, win_thrift_alpha4,&
             win_thrift_matld, win_thrift_matmd, win_thrift_matud, win_thrift_matrhs
    REAL(rprec) :: tstart, tend, jtol, picard_factor
    REAL(rprec), DIMENSION(:), POINTER :: THRIFT_RHO(:), THRIFT_T(:), THRIFT_RHOFULL(:)
    REAL(rprec), DIMENSION(:,:), POINTER :: THRIFT_J, THRIFT_JPLASMA, THRIFT_JBOOT,&
                 THRIFT_JECCD, THRIFT_JNBCD, THRIFT_JOHMIC, THRIFT_JSOURCE, &
                 THRIFT_I, THRIFT_IPLASMA, THRIFT_IBOOT, THRIFT_IECCD, THRIFT_INBCD, &
                 THRIFT_IOHMIC, THRIFT_ISOURCE, THRIFT_UGRID, THRIFT_VP,&
                 THRIFT_BAV, THRIFT_BSQAV, &
                 THRIFT_AMINOR, THRIFT_RMAJOR, &
                 THRIFT_S11,THRIFT_S12, THRIFT_S22, THRIFT_PHIEDGE, &
                 THRIFT_COEFF_A, THRIFT_COEFF_B, THRIFT_COEFF_C, THRIFT_COEFF_D,&
                 THRIFT_COEFF_BP, THRIFT_COEFF_CP, THRIFT_COEFF_DP,&
                 THRIFT_ALPHA1,THRIFT_ALPHA2,THRIFT_ALPHA3,THRIFT_ALPHA4,&
                 THRIFT_MATLD, THRIFT_MATMD, THRIFT_MATUD, THRIFT_MATRHS

    ! For TRAVIS
    INTEGER, PARAMETER :: nsys   = 16
    INTEGER :: nra_ecrh, nphi_ecrh
    INTEGER, DIMENSION(nsys)     :: wmode_ecrh
    REAL(rprec), DIMENSION(nsys) :: freq_ecrh, power_ecrh
    REAL(rprec), DIMENSION(nsys,3)     :: antennaposition_ecrh, &
                 targetposition_ecrh,rbeam_ecrh,rfocus_ecrh

END MODULE thrift_vars
