!-----------------------------------------------------------------------
!     Module:        pies_fieldlines
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/8/2011
!     Description:   This module contains the variables utilized by the
!                    field line following routines.
!-----------------------------------------------------------------------
      MODULE pies_fieldlines
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
      INTEGER                  :: last_line, hitsrf, nfollow, nintw, &
                                  ninteqs, nu_spline, nv_spline, ik_help
      INTEGER                  :: iuser(3)
      INTEGER, ALLOCATABLE     :: hitwal(:)
      REAL                     :: follow_tol
      REAL                     :: dkmin
      REAL                     :: magaxis_abs
      REAL                     :: magaxis_eta
      REAL                     :: r_axis
      REAL                     :: z_axis
      REAL                     :: dphi
      REAL                     :: phgaus
      REAL                     :: tfrac
      REAL(rprec), ALLOCATABLE :: psjmnc(:,:)
      REAL(rprec), ALLOCATABLE :: psjmns(:,:)
      REAL(rprec), ALLOCATABLE :: bthmnc(:,:)
      REAL(rprec), ALLOCATABLE :: bthmns(:,:)
      REAL(rprec), ALLOCATABLE :: brhomnc(:,:)
      REAL(rprec), ALLOCATABLE :: brhomns(:,:)
      REAL(rprec), ALLOCATABLE :: rhobth(:,:)
      REAL(rprec), ALLOCATABLE :: R_array(:)
      REAL(rprec), ALLOCATABLE :: Z_array(:)
      REAL(rprec), ALLOCATABLE :: bxsi_array(:)
      REAL(rprec), ALLOCATABLE :: beta_array(:)  
      REAL(rprec), ALLOCATABLE :: Rln(:,:)
      REAL(rprec), ALLOCATABLE :: Zln(:,:) 
      REAL(rprec), ALLOCATABLE :: rholn(:,:)
      REAL(rprec), ALLOCATABLE :: thetaln(:,:)
      REAL(rprec), ALLOCATABLE :: philn(:,:)
      REAL(rprec), ALLOCATABLE :: xsiln(:,:)
      REAL(rprec), ALLOCATABLE :: etaln(:,:)
      TYPE(EZspline3_r8)       :: brho_spl
      TYPE(EZspline3_r8)       :: theta_spl
      TYPE(EZspline3_r8)       :: beta_spl
      TYPE(EZspline3_r8)       :: bxsi_spl
      TYPE(EZspline3_r8)       :: R_spl
      TYPE(EZspline3_r8)       :: Z_spl
      TYPE(EZspline3_r8)       :: BS_spl
      TYPE(EZspline3_r8)       :: BU_spl
      TYPE(EZspline3_r8)       :: BV_spl
      
      END MODULE pies_fieldlines
