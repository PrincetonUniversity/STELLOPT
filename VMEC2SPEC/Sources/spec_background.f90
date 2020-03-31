!-----------------------------------------------------------------------
!     Module:        spec_background
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/09/2011
!     Description:   This module contains the background coordinate
!                    variables used by SPEC.
!-----------------------------------------------------------------------
      MODULE spec_background
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds
!-----------------------------------------------------------------------
!     Module Variables
!          nfp        Number of Field Periods
!          m          Number of Poloidal Modes
!          n          Number of Toroidal Modes
!          k          Number of Radial Surfaces
!          mnmax      Total number of Fourier Modes
!          xm         Poloidal Mode Number Array
!          xn         Toroidal Mode Number Array
!          iotamn     Minimum value of iota
!          iotamx     Maximum value of iota
!          rho        Radial Points Array [0,1]
!          rmnc       Even R Fourier Modes (cos)
!          rmns       Odd R Fourier Modes (sin)
!          zmnc       Even Z Fourier Modes (cos)
!          zmns       Odd Z Fourier Modes (sin)      
!----------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, PARAMETER       :: nvol_mx = 100
      INTEGER                  :: nfp, m, n, k, mnmax, signgs, nvol
      INTEGER                  :: ni(1:nvol_mx), pl(0:nvol_mx), &
                                  ql(0:nvol_mx), pr(0:nvol_mx), &
                                  qr(0:nvol_mx)
      INTEGER, ALLOCATABLE     :: xm(:), xn(:)
      REAL                     :: iotamn,iotamx,rbtor_pies
      REAL(rprec)              :: tflux(0:nvol_mx), pflux(0:nvol_mx), &
                                  mu(1:nvol_mx), iota(0:nvol_mx)
      REAL(rprec), ALLOCATABLE :: rho(:)
      REAL(rprec), ALLOCATABLE :: rmnc(:,:), rmns(:,:), zmnc(:,:), zmns(:,:)
      
      END MODULE spec_background
