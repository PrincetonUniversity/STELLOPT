!-----------------------------------------------------------------------
!     Module:        beams3d_lines
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/21/2012
!     Description:   This module contains the FIELDLINES field line
!                    variables.
!-----------------------------------------------------------------------
      MODULE beams3d_lines
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      
!-----------------------------------------------------------------------
!     Module Variables
!          myline    Dummy index
!          nparticles    Number of Particles
!          nsteps    Number of integration steps along fieldline
!          R_lines   Radial locations along fieldline [m] (npoinc per field period)
!          Z_lines   Vertical locations along field line [m]
!          PHI_lines Toroidal locations along field line [radians]
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL  ::  ltherm
      LOGICAL, ALLOCATABLE :: lost_lines(:)
      INTEGER  :: ns_prof = 100
      INTEGER  :: nparticles, nsteps, myline, mybeam, mytdex, myend, mystart_save,myend_save
      REAL(rprec) :: xlast,ylast,zlast ! for storing position
      REAL(rprec) :: moment, mycharge, myZ, mymass, myv_neut(3), &
                     B_temp(4), rand_prob, cum_prob, tau, next_t, &
                     dt_out, partvmax
      LOGICAL, ALLOCATABLE     :: neut_lines(:,:)
      REAL(rprec), ALLOCATABLE :: shine_through(:)
      REAL(rprec), ALLOCATABLE :: ndot_prof(:,:),power_prof(:,:),epower_prof(:,:),ipower_prof(:,:),j_prof(:,:)
      REAL(rprec), ALLOCATABLE :: R_lines(:,:),Z_lines(:,:),PHI_lines(:,:),vll_lines(:,:),moment_lines(:,:),&
                                  S_lines(:,:),U_lines(:,:),PI_lines(:,:), PE_lines(:,:),&
                                  B_lines(:,:),j_lines(:,:)

      END MODULE beams3d_lines
