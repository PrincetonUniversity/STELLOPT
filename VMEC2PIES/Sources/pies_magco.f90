!-----------------------------------------------------------------------
!     Module:        pies_magco
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/6/2011
!     Description:   This module contains the magnetic coordinate arrays.
!-----------------------------------------------------------------------
      MODULE pies_magco
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds
!-----------------------------------------------------------------------
!     Module Variables
!          nu         Number of real space poloidal points
!          nv         Number of real space toroidal points
!          rreal      Real Space R
!          zreal      Real Space Z
!         
!----------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER                  :: nu_magco, nv_magco
      INTEGER                  :: mnmax_polar,nu_polar,nv_polar
      INTEGER, ALLOCATABLE     :: xm_polar(:)
      INTEGER, ALLOCATABLE     :: xn_polar(:)
      REAL(rprec)              :: xv_polar(1)
      REAL(rprec), ALLOCATABLE :: xu_polar(:)
      REAL(rprec), ALLOCATABLE :: rho_polar(:,:,:)
      REAL(rprec), ALLOCATABLE :: rholn_polar(:,:,:)
      REAL(rprec), ALLOCATABLE :: rhomnc_polar(:,:),rhomns_polar(:,:)
      REAL(rprec), ALLOCATABLE :: rmnc_magco(:,:), rmns_magco(:,:)
      REAL(rprec), ALLOCATABLE :: zmnc_magco(:,:), zmns_magco(:,:)
      REAL(rprec), ALLOCATABLE :: rreal_magco(:,:,:), zreal_magco(:,:,:)
      REAL(rprec), ALLOCATABLE :: bsreal_magco(:,:,:)
      REAL(rprec), ALLOCATABLE :: bureal_magco(:,:,:)
      REAL(rprec), ALLOCATABLE :: bvreal_magco(:,:,:)
      REAL(rprec), ALLOCATABLE :: bsmnc_magco(:,:), bsmns_magco(:,:)
      REAL(rprec), ALLOCATABLE :: bumnc_magco(:,:), bumns_magco(:,:)
      REAL(rprec), ALLOCATABLE :: bvmnc_magco(:,:), bvmns_magco(:,:)
      
      END MODULE pies_magco
