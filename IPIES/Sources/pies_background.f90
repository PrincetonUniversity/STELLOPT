!-----------------------------------------------------------------------
!     Module:        pies_background
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/2/2011
!     Description:   This module contains the background coordinate
!                    variables used by PIES.
!-----------------------------------------------------------------------
      MODULE pies_background
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
!          nu         Number of real space poloidal points
!          nv         Number of real space toroidal points
!          mnmax      Total number of Fourier Modes
!          xm         Poloidal Mode Number Array
!          xn         Toroidal Mode Number Array
!          iotamn     Minimum value of iota
!          iotamx     Maximum value of iota
!          rho        Radial Points Array [0,1]
!          xu         Poloidal Points Number Array [0,1]
!          xv         Toroidal Points Number Array [0,1]
!          rmnc       Even R Fourier Modes (cos)
!          rmns       Odd R Fourier Modes (sin)
!          zmnc       Even Z Fourier Modes (cos)
!          zmns       Odd Z Fourier Modes (sin)
!          bsmnc      Even B^s Fourier Modes
!          bsmns      Odd B^s Fourier Modes
!          bumnc      Even B^u Fourier Modes
!          bumns      Odd B^u Fourier Modes
!          bvmnc      Even B^v Fourier Modes
!          bvmns      Odd B^v Fourier Modes
!          jsmnc      Even J^s Fourier Modes
!          jsmns      Odd J^s Fourier Modes
!          jumnc      Even J^u Fourier Modes
!          jumns      Odd J^u Fourier Modes
!          jvmnc      Even J^v Fourier Modes
!          jvmns      Odd J^v Fourier Modes         
!----------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER                  :: nfp, m, n, k, mnmax
      INTEGER, ALLOCATABLE     :: xm(:), xn(:)
      REAL                     :: iotamn,iotamx
      REAL(rprec), ALLOCATABLE :: rho(:)
      REAL(rprec), ALLOCATABLE :: rmnc(:,:), rmns(:,:), zmnc(:,:), zmns(:,:)
      REAL(rprec), ALLOCATABLE :: bsmnc(:,:), bsmns(:,:)
      REAL(rprec), ALLOCATABLE :: bumnc(:,:), bumns(:,:)
      REAL(rprec), ALLOCATABLE :: bvmnc(:,:), bvmns(:,:)
      REAL(rprec), ALLOCATABLE :: jsmnc(:,:), jsmns(:,:)
      REAL(rprec), ALLOCATABLE :: jumnc(:,:), jumns(:,:)
      REAL(rprec), ALLOCATABLE :: jvmnc(:,:), jvmns(:,:)
      
      END MODULE pies_background
