!-----------------------------------------------------------------------
!     Module:        torlines_realspace
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/2/2011
!     Description:   This module contains the realspace arrays.
!-----------------------------------------------------------------------
      MODULE torlines_realspace
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds
!-----------------------------------------------------------------------
!     Module Variables
!          nu         Number of real space poloidal points
!          nv         Number of real space toroidal points
!          nuv        nu*nv
!          nuvp       nu*nv*np
!          nv_vmec    Number of toroidal points in vacuum field file
!          xu         Poloidal Points Number Array [0,1]
!          xv         Toroidal Points Number Array [0,1]
!          rreal      Real Space R
!          zreal      Real Space Z
!          bsreal     Real Space B^s
!          bureal     Real Space B^u
!          bvreal     Real Space B^v
!          rs         Real Space dR/ds
!          ru         Real Space dR/du
!          rv         Real Space dR/dv
!          zs         Real Space dZ/ds
!          zu         Real Space dZ/du
!          zv         Real Space dZ/dv
!         
!----------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER                  :: nrho, nu, nv, nuv, nuvp, nv_vmec, ik_sav
      INTEGER                  :: win_xu, win_xv, win_rreal, win_zreal, &
                                  win_bsreal, win_bureal, win_bvreal, &
                                  win_breal, win_brreal, win_bphireal, &
                                  win_bzreal
      REAL(rprec), POINTER :: xu(:), xv(:)
      REAL(rprec), POINTER :: rreal(:,:,:), zreal(:,:,:)
      REAL(rprec), POINTER :: bsreal(:,:,:), bureal(:,:,:), bvreal(:,:,:)
      REAL(rprec), POINTER :: breal(:,:,:)
      REAL(rprec), POINTER :: brreal(:,:,:), bphireal(:,:,:), bzreal(:,:,:)
      REAL(rprec), ALLOCATABLE :: rs(:,:,:),zs(:,:,:)
      REAL(rprec), ALLOCATABLE :: ru(:,:,:),rv(:,:,:)
      REAL(rprec), ALLOCATABLE :: zu(:,:,:),zv(:,:,:)
      REAL(rprec), ALLOCATABLE :: g11(:,:,:),g12(:,:,:),g13(:,:,:)
      REAL(rprec), ALLOCATABLE :: g22(:,:,:),g23(:,:,:),g33(:,:,:)
      REAL(rprec), ALLOCATABLE :: detg(:,:,:)
      REAL(rprec), ALLOCATABLE :: jsreal(:,:,:), jureal(:,:,:)
      REAL(rprec), ALLOCATABLE :: jvreal(:,:,:)
      
      END MODULE torlines_realspace
