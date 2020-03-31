!----------------------------------------------------------------------
!  A module to hold some values shared among MAKEGRID source files.
!----------------------------------------------------------------------
      MODULE makegrid_global
!---------------------------------------------------------------------------
!  Use statements are followed by variables and subroutines used in this file
!  (subroutines are followed by parentheses)
!---------------------------------------------------------------------------
      USE stel_kinds  ! rprec
      USE vsvd0, ONLY: nigroup
      
      IMPLICIT NONE
      
      LOGICAL           :: lscreen
      CHARACTER(LEN=20) :: task

!---------------------------------------------------------------------------
!  Will need some arrays sized for the coil groups
!  Define a parameter that should be large enough.
!---------------------------------------------------------------------------
      INTEGER, PARAMETER :: nextcur_dim = nigroup

!---------------------------------------------------------------------------
!  Variables A and B on the surface of a circular torus
!  (Used to generate input for NIMROD, CTH simulations)
!    TASK: circ_tor_grid
!---------------------------------------------------------------------------

      REAL(rprec) :: rmajor, aminor
      INTEGER :: nphi, ntheta
      REAL(rprec), DIMENSION(nextcur_dim) :: extcur_mgrid

!---------------------------------------------------------------------------
!  Variables for shifts and rotations of coil_groups
!    TASK: mgrid_rs
!---------------------------------------------------------------------------
      REAL(rprec), DIMENSION(nextcur_dim,3) ::  cg_shift_1,                    &
     &   cg_shift_2, cg_rot_xcent
      REAL(rprec), DIMENSION(nextcur_dim) :: cg_rot_theta,                     &
     &   cg_rot_phi, cg_rot_angle
      LOGICAL, DIMENSION(nextcur_dim) :: l_rot_coil_center
      
      END MODULE makegrid_global
