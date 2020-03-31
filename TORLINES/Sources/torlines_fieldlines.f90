!-----------------------------------------------------------------------
!     Module:        pies_fieldlines
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/8/2011
!     Description:   This module contains the variables utilized by the
!                    field line following routines.
!-----------------------------------------------------------------------
      MODULE torlines_fieldlines
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Module Variables
!          last_line    Last surface on which field lines are followed
!          hitsrf       First surface to hit the wall
!          nfollow
!          nintw        Total number of steps in phi for field line following
!          ninteqs      Size of q helper array (number of equations, two for each surface)
!          nu_spline    Number of polodial points in spline
!          nv_spline    Number of toroidal points in spline
!          hitwal       Array indicating if a field line has hit the wall
!          follow_tol   Tollerance for field line following
!          dkmin        Smallest distance to be resolved in Fourier Space
!          nfollow      Number of step in phi
!          magaxis_abs  Percision of magnetic axis solver
!          magaxis_eta  Value of |f(x)| is taken to be a zero
!          r_axis       Location of magnetic axis in phi=0 plane (in background coords)
!          z_axis       Location of magnetic axis in phi=0 plane (in background coords)
!          dphi         Change in phi when following field lines
!          psjmnc       Psi-Jacobian (1/B^v) (cos)
!          psjmns       Psi-Jacobian (1/B^v) (sin)
!          brho         B^s/B^v
!          btheta       B^u/B^v
!          rhobth       rho*B^u/B^v
!          bxsi_array   Array for field line follower
!          beta_array   Array for field line follower
!          rho_spl      Spline for rho
!          theta_spl    Spline for theta
!          beta_spl     Spline for brho
!          bxsi_spl     Spline for rhobth
!          xsiln        XSI of field line (r coordinate)
!          etaln        ETA of field line (z coordinate)
!         
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: BYTE_8_loc = SELECTED_INT_KIND (8)
      INTEGER(KIND=BYTE_8_loc) :: mystart, myend
      INTEGER  :: nlines, nsteps, myline, mylines, myla, mylb, &
                  myldex
      REAL(rprec) :: phimn
      
      LOGICAL, ALLOCATABLE :: goodline(:)
      REAL(rprec), ALLOCATABLE :: R_lines(:,:),Z_lines(:,:),PHI_lines(:,:),B_lines(:,:),U_lines(:,:)
                                                                
      REAL(rprec), ALLOCATABLE :: Rhc_lines(:,:),Zhc_lines(:,:)
      
      
      END MODULE torlines_fieldlines
