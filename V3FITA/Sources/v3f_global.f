!     SPH: INTEGER(iprec) -> INTEGER
!*******************************************************************************
!  File v3f_global.f
!  Contains module v3f_global
!  Contains declarations for variables that are globally accessible throughout
!  the v3fit code

!*******************************************************************************
!  MODULE v3f_global
!    (Global Variable declarations, for the V3FIT code)
! SECTION I.    VARIABLE DECLARATIONS

! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE v3f_global

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!  Frequently used mathematical constants, lots of extra precision.
!-------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
     
!-------------------------------------------------------------------------------
!  Use Statements for other structures, ? V3 Utilities ?
!-------------------------------------------------------------------------------
      USE bsc_T
      USE diagnostic_T
      USE geometric_T
      USE signal_T
!     SPH ADDED FOR GLOBAL EQUILIBRIUM STATE VARIABLES
      USE eq_T                   
!      USE eq_interface  ! commented out 2008-08-05
      USE model_T
      USE recon_param_T

!-------------------------------------------------------------------------------
!  Implicit None comes after USE statements, before other declarations
!-------------------------------------------------------------------------------
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Diagnostic Descriptions
!-------------------------------------------------------------------------------

      INTEGER :: na_d_desc_default = 1000
      INTEGER :: na_d_desc = 1000, n_d_desc = 0
      INTEGER :: n_d_desc_mddc = 0
      INTEGER :: n_d_desc_sxrc = 0      ! GJH added 2010-01-21
      INTEGER :: n_d_desc_ipsl = 0
      TYPE(diagnostic_desc), DIMENSION(:), ALLOCATABLE, TARGET ::              &
     &    d_desc

!-------------------------------------------------------------------------------
!  Geometric Descriptions
!-------------------------------------------------------------------------------

      INTEGER :: na_g_desc_default = 1000
      INTEGER :: na_g_desc = 1000, n_g_desc = 0
      INTEGER :: n_g_desc_iso = 0
      TYPE(geometric_desc), DIMENSION(:), ALLOCATABLE, TARGET ::               &
     &    g_desc
      
!-------------------------------------------------------------------------------
!  Signal Descriptions
!-------------------------------------------------------------------------------

      INTEGER :: na_s_desc_default = 1000
      INTEGER :: na_s_desc = 1000, n_s_desc = 0
      INTEGER :: n_s_desc_mddc = 0
      INTEGER :: n_s_desc_sxrc = 0          ! GJH added 2010-01-21
      INTEGER :: n_s_desc_ipsl = 0
      INTEGER :: n_s_desc_geom = 0
      TYPE(signal_desc), DIMENSION(:), ALLOCATABLE, TARGET ::                  &
     &    s_desc
            
!-------------------------------------------------------------------------------
!  Reconstruction Parameters
!-------------------------------------------------------------------------------

      INTEGER :: na_rparam_default = 100
      INTEGER :: na_rparam = 100, n_rparam = 0
      TYPE(recon_param), DIMENSION(:), ALLOCATABLE, TARGET ::                  &
     &    rparam
            
!-------------------------------------------------------------------------------
!  Reconstruction Constraints
!-------------------------------------------------------------------------------

      TYPE(recon_cnstrnts)       ::  rcnstrnts 
      
!-------------------------------------------------------------------------------
!  Observed Signal Data
!-------------------------------------------------------------------------------
!  Arrays for NLI of data and sigma (sdo_data_a and sdo_sigma_a) declared here
!  iw_sdo_verbose   Verbosity level for writing signal_data observe, -1 no write

      INTEGER :: na_sdata_o_default = 1000
      INTEGER :: na_sdata_o = 1000, n_sdata_o = 0
      INTEGER :: iw_sdo_verbose = -1
      REAL(rprec), DIMENSION(1000)    :: sdo_data_a
      REAL(rprec), DIMENSION(1000)    :: sdo_sigma_a
      TYPE(signal_data), DIMENSION(:), ALLOCATABLE, TARGET ::                  &
     &    sd_observe
            
!-------------------------------------------------------------------------------
!  Command Line Parsing (including debug logical)
!-------------------------------------------------------------------------------

      INTEGER :: num_cl_args
      CHARACTER(LEN=80), DIMENSION(10) :: command_line_args
      LOGICAL :: l_debug

!-------------------------------------------------------------------------------
!  File Names
!-------------------------------------------------------------------------------

!  JDH 2009-04-07. Make sure all filenames are initialized.
      CHARACTER(LEN=8)   :: main_nli_filename_default = 'v3fit.in'
      CHARACTER(LEN=300) :: main_nli_filename = ''
      CHARACTER(LEN=300) :: mdsig_list_filename = ''
      CHARACTER(LEN=300) :: intpol_filename = ''
      CHARACTER(LEN=300) :: vmec_nli_filename = ''
      CHARACTER(LEN=300) :: geometric_filename = ''
      CHARACTER(LEN=300) :: runlog_filename = ''
      CHARACTER(LEN=300) :: recout_filename = ''
      CHARACTER(LEN=300) :: sxr_nli_filename=''         !GJH 8/31/09
!      
!-------------------------------------------------------------------------------
!  I/O Unit Numbers
!-------------------------------------------------------------------------------

!  2010-06-07 Changed numbers, added info on VMEC.
!  Use of safe_open should avoid collisions
      INTEGER :: iou_mnli = 30
      INTEGER :: iou_mdsig_list = 31
      INTEGER :: iou_mdsig = 32
      INTEGER :: iou_runlog = 33
      INTEGER :: iou_recout = 34
      INTEGER :: iou_sxr_nli = 35                       !GJH 9/31/09

!  VMEC IOU Unit Numbers
!     8    nfort8                vparams.f
!     9    nthreed0              vparams.f 
!    10    nmac0 (nthreed0+1)    vparams.f
!    11    indata0 (nthreed0+2)  vparams.f
!    12    nwout0 (nthreed0+3)   vparams.f
!    13    jxbout0 (nthreed0+3)  vparams.f
!    18    nfort18               vparams.f
!    24    jxbout text file      jxbforce.f
!    51    nlog0                 vparams.f
!    52    nmercier0             vparams.f

!-------------------------------------------------------------------------------
!  Task information
!-------------------------------------------------------------------------------
      CHARACTER(LEN=20) :: my_task = 'run_vmec'
      LOGICAL :: l_zero_xcdot = .TRUE.
!  l_zero_xcdot   logical, to control zeroing of xcdot using
!     subroutine eq_change_vmi_zero_xcdot
!     After testing on 2008-08-05, changed from default .FALSE.
!     (previous behavior) to default .TRUE.
!     
!-------------------------------------------------------------------------------
!  Model information
!-------------------------------------------------------------------------------
      TYPE (model), SAVE, TARGET     :: model_a
!  SAVE because Lahey Fujitsu on Logjam says
! "'model_a' shall have the SAVE attribute; has type for which component
!  initialization was specified in a MODULE."
!  JDH 11-26-2004. Modified for model_a on 12-14-2004

!-------------------------------------------------------------------------------
!  Work Variables
!-------------------------------------------------------------------------------
!  These variables are to be used as needed by various tasks, so that in the
!  process of code development, one does not need to keep declaring variables.
!  These are in the main NLI, and can be used by any task. 
!  After things have settled down, separate, well-named variables should be
!  used. (Added JDH 06-28-06)

      INTEGER, DIMENSION(100)                 :: i_work
      REAL(rprec), DIMENSION(100)             :: r_work
      CHARACTER(LEN=100), DIMENSION(100)      :: c_work

!-------------------------------------------------------------------------------
! Density NML variables  
!-------------------------------------------------------------------------------
!   These variables are used to specify the density profile for the INITIAL iteration
!   (only) of the model, for use in reconstruction involving Interferometry/Polarimetry 
!   signals.  Note that once the values pass from the NML to the model, these global variables
!   cease to be utilized.  The density profile is assumed to be of form:
!           N = N_max*(1 - s^tau)^kappa + N_ambient       (added 7/24/07 by JMS)

      REAL(rprec) :: density_max_default_g     = 4.0e18
      REAL(rprec) :: density_ambient_default_g = 0.01_rprec
      REAL(rprec) :: density_tau_default_g     = 2.0_rprec
      REAL(rprec) :: density_kappa_default_g   = 1.0_rprec

      REAL(rprec) :: density_max_g
      REAL(rprec) :: density_ambient_g
      REAL(rprec) :: density_tau_g
      REAL(rprec) :: density_kappa_g

!-------------------------------------------------------------------------------
! geometric NLI variables  
!-------------------------------------------------------------------------------

      REAL(rprec) :: r_major_radius            = zero
      REAL(rprec) :: a_minor_radius            = zero
      
      INTEGER, PARAMETER                       :: na_lif = 10
      INTEGER, PARAMETER                       :: na_phi_lif = 10
      INTEGER, PARAMETER, PRIVATE :: n_lif_arz = 25 * na_lif
      INTEGER, PARAMETER, PRIVATE :: n_lif_pd = na_phi_lif * na_lif

      INTEGER                                  :: n_lif = 0
      INTEGER, DIMENSION(na_lif)               :: n_phi_lif = zero
      REAL(rprec), DIMENSION(na_lif,0:4,0:4)   :: lif_arz = zero
      REAL(rprec), DIMENSION(na_lif)           :: lif_rc = zero
      REAL(rprec), DIMENSION(na_lif)           :: lif_zc= zero
      REAL(rprec), DIMENSION(na_lif)           :: lif_sigma = 0.001D00
      REAL(rprec), DIMENSION(na_lif,na_phi_lif):: lif_phi_degree = zero
      LOGICAL, DIMENSION(na_lif)               :: lif_on_edge = .false.
      
      END MODULE v3f_global
!-------------------------------------------------------------------------------
! Changes
!-------------------------------------------------------------------------------
!  
!  JDH 2007-10-06
!    Added recon_cnstrnts named rcnstrnts 
!
!  JDH 2008-08-04
!    Added l_zero_xcdot
!
!  JDH 2009-02-01
!     Added things for geometric
!
!  JDH 2009-04-07
!     Added initialization for file names.
!
!  JDH 2009-04-13
!     Added lif_ geometric NLI variables
! 
!  JDH 2010-06-10
!     iou and filenames for recout and runlog. Added VMEC iou comment
!
!  GJH 2010-01-21 (JDH 2010-06-11)
!     Added sxr_nli_filename for sxr camera input
!     Added iuo_sxr_nli
!     Added n_d_desc_sxrc                   
!     Added n_s_desc_sxrc
!
!  JDH 2010-07-20
!    geometric NLI - eliminated DATA statements - initialized in declarations
!------------------------------------------------------------------------------
