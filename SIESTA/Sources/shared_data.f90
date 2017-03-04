!>  \brief Module containing data that is shared with other modules
      MODULE shared_data
      USE stel_kinds
      USE stel_constants

      INTEGER, PARAMETER :: ngmres_steps = 100      !<max number of gmres steps (10-100)
      INTEGER, PARAMETER :: ngmres_type = 2         !<=1, no "peek" improvement, =2 "peek" performance improvement
      INTEGER, PARAMETER :: PREDIAG=1               !<diagonal preconditioning
      INTEGER, PARAMETER :: PREBLOCK=2              !<block preconditioning
      INTEGER            :: neqs !<No. elements in xc array
      INTEGER            :: nprecon
      INTEGER            :: nprecon_type
      INTEGER            :: niter

      REAL(rprec),PARAMETER :: fsq_res = 1.E-16_dp  !<threshold force for turning off resistive perturbations
      REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: gc  !<1D Array of Fourier mode MHD force components
      REAL(rprec), DIMENSION(:), ALLOCATABLE         :: xc  !<1D array of Fourier mode displacement components
      REAL(rprec), DIMENSION(:), ALLOCATABLE         :: xc0, gc0

      REAL(rprec)        :: fsq_total, fsq_total1, fsq_lin, etak_tol
      REAL(rprec)        :: wtotal, wtotal0
      REAL(rprec)        :: delta_t
      REAL(rprec)        :: fsqvs, fsqvu, fsqvv
      REAL(rprec)        :: xPosDef
      REAL(rprec)        :: ste(4)                  !<Spectral Truncation RMS Error    
      REAL(rprec)        :: bs0(6), bu0(6)
      REAL(rprec)        :: bsbu_ratio, jsju_ratio

!Logicals for evolving origin and edge
!      LOGICAL, PARAMETER :: l_push_edge=.TRUE.     !<=TRUE, u,v components of forces IS solved at s=1
      LOGICAL, PARAMETER :: l_push_edge=.FALSE.     !<=FALSE, u,v components of forces NOT solved at s=1
      LOGICAL           :: l_push_s, l_push_u, l_push_v !<=FALSE, s,u,v components of forces NOT solved at s=0
!      LOGICAL, PARAMETER :: l_natural=.TRUE.       !<=TRUE use natural bc at s=0 (preferred: maintains symmetry)
      LOGICAL, PARAMETER :: l_natural=.FALSE.       !<=FALSE evolves all forces at origin
      LOGICAL            :: l_linearize
      LOGICAL            :: l_getwmhd
      LOGICAL            :: l_getfsq
      LOGICAL            :: l_ApplyPrecon           !<controls when preconditioner is applied
      LOGICAL            :: l_PrintOriginForces=.FALSE.
      LOGICAL            :: l_init_state
      LOGICAL            :: l_update_state=.FALSE.  !<=TRUE before funct_island called to get ste array
      LOGICAL            :: l_par_state             !<=TRUE if parallel allocation state (in init_state_par)

      END MODULE shared_data
