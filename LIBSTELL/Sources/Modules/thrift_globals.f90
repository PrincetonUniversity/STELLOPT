!-----------------------------------------------------------------------
!     Module:        thrift_globals
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          10/16/2024
!     Description:   This module contains the THRIFT global variables
!                    needed by the input namelist.
!-----------------------------------------------------------------------
      MODULE thrift_globals
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      
!-----------------------------------------------------------------------
!     Module Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE

      ! Moved from thrift_vars
      LOGICAL ::  lverbj, leccd, lnbcd, lohmic
      INTEGER ::  nrho, ntimesteps, n_eq, npicard, nsj
      REAL(rprec) :: tstart, tend, jtol, picard_factor, boot_factor

      ! Moved from thrift_vars (For ECCD in general)
      INTEGER, PARAMETER :: ntime_ecrh = 200
      REAL(rprec), DIMENSION(ntime_ecrh) :: PECRH_AUX_T, PECRH_AUX_F
      REAL(rprec) :: ecrh_rc, ecrh_w

      ! Moved from thrift_vars (for TRAVIS)
      INTEGER, PARAMETER :: nsys   = 16
      INTEGER :: nra_ecrh, nphi_ecrh
      INTEGER, DIMENSION(nsys)     :: wmode_ecrh
      REAL(rprec), DIMENSION(nsys) :: freq_ecrh, power_ecrh
      REAL(rprec), DIMENSION(nsys,3)     :: antennaposition_ecrh, &
                 targetposition_ecrh, rbeam_ecrh, rfocus_ecrh

      ! Moved from thrift_vars (For DKES)
      INTEGER, PARAMETER :: DKES_NS_MAX = 64
      INTEGER, PARAMETER :: DKES_NSTAR_MAX = 32
      INTEGER :: nruns_dkes
      INTEGER, DIMENSION(:), POINTER :: DKES_rundex
      INTEGER, DIMENSION(DKES_NS_MAX) :: DKES_K
      REAL(rprec), DIMENSION(DKES_NSTAR_MAX) :: DKES_Erstar, DKES_Nustar

      ! Moved from thrift_runtime
      INTEGER :: nparallel_runs, mboz, nboz
      CHARACTER(256) :: bootstrap_type, eccd_type, vessel_ecrh, &
                        mirror_ecrh, targettype_ecrh, antennatype_ecrh, &
                        etapar_type



      CONTAINS

      ! These expose the global variables through ctypes
      INTEGER FUNCTION getmaxtimeecrh()
      IMPLICIT NONE
      getmaxtimeecrh = ntime_ecrh
      END FUNCTION getmaxtimeecrh
      
      INTEGER FUNCTION getmaxsys()
      IMPLICIT NONE
      getmaxsys = nsys
      END FUNCTION getmaxsys
      
      INTEGER FUNCTION getmaxns()
      IMPLICIT NONE
      getmaxns = DKES_NS_MAX
      END FUNCTION getmaxns
      
      INTEGER FUNCTION getmaxnstar()
      IMPLICIT NONE
      getmaxnstar =DKES_NSTAR_MAX
      END FUNCTION getmaxnstar

      END MODULE thrift_globals
