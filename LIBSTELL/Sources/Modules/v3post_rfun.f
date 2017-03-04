      MODULE v3post_rfun
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  BUILD RESPONSE FUNCTION TABLES                **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          Response table variables                                **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          12/08/02..........first created                         **
!**                                                                  **
!**********************************************************************
      USE stel_kinds
!----------------------------------------------------------------------
!--   idrfun            : identification of external coil responses  --
!--   eqtype            : plasma equilibrium identification          --
!--   myid              : user's identification          --
!--   ideqfile          ; identification of plasma equilibrium file  --
!--   eq_file           : equilibrium data file name                 --
!--   n_field_c         : # of coil                                  --
!--   n_field_cg        : # of coil groups                           --
!--   nflx_loop         : # of thin flux loops                       --
!--   nprobe            : # of magnetic probes                       --
!--   n_diagn_c         : # of diagnostic sensors                    --
!----------------------------------------------------------------------
      CHARACTER(len=40)           :: idrfun
      CHARACTER(len=20)           :: eqtype
      CHARACTER(len=20)           :: myid=' '
      CHARACTER(len=50)           :: ideqfile
      CHARACTER(len=70)           :: eq_file
      LOGICAL 			        :: freeb = .true.
      LOGICAL                     :: lsurf = .false.
      INTEGER(iprec)              :: n_field_c
!      INTEGER(iprec)             :: n_field_cg
      INTEGER(iprec)              :: nflx_loop
      INTEGER(iprec)              :: nprobe
!      INTEGER(iprec)             :: n_diagn_c

!  JDH Put some variables here 06.28.2003
      INTEGER(iprec)              :: n_diagn_c
      INTEGER(iprec)              :: n_field_periods


!----------------------------------------------------------------------
!-- response functions due to external coils                         --
!--   rdiag_coilg         : diagnostic coils due to coil groups      --
!----------------------------------------------------------------------
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: rdiag_coilg
!----------------------------------------------------------------------
!-- response functions due to plasma                                 --
!--   name_coils_dot    : name of the 'coils.' file                  --
!--   lstell_sym        : logical - True for stellarator symmetry    --
!--   rmin              : minimum R in m for plasma grid             --
!--   rmax              : maximum R in m for plasma grid             --
!--   zmin              : minimum Z in m for plasma grid             --
!--   zmax              : maximum Z in m for plasma grid             --
!--   ir                : number of grid points in R plasma grid     --
!--   jz                : number of grid points in Z plasma grid     --
!--   kp                : number of phi planes in plasma grid        --
!--   kp_store          : number of phi planes stored in prf grid    --
!--   n_field_periods   : number of field periods                    --
!--   pl_response       : cylindrical components of plasma response  --
!--   rgrid             : R on cylindrical grid in m                 --
!--   zgrid             : Z on cylindrical grid in m                 --
!--   pgrid             : Phi on cylindrical grid in radian          --
!--   delz,delz,delp    : (R,P,Z) grid sizes                         --
!--   sgrid             : radial equilibrium grid                    --
!--   ugrid             : poloidal angular equilibrium grid          --
!--   vgrid             : toroidal angular equilibrium grid          --
!--   dels,delu,delv    : (s,u,v) grid sizes                         --
!--   ju                : number of grid points in u equilibrium grid--
!--   kv                : number of grid points in v equilibrium grid--
!--   rsuv              : R at (s,u,v)                               --
!--   psuv              : Phi at (s,u,v)                             --
!--   zsuv              : Z at (s,u,v)                               --
!--   gsuv              : SQRT(g) at (s,u,v)                         --
!--   currusuv          : J^u at (s,u,v)                             --
!--   currvsuv          : J^v at (s,u,v)                             --
!--   bsubu             : B_u at (s,u,v)                             --
!--   bsubv             : B_v at (s,u,v)                             --
!--   rusuv             : dR/du at (s,u,v)                           --
!--   rvsuv             : dR/dv at (s,u,v)                           --
!--   zusuv             : dZ/du at (s,u,v)                           --
!--   zusuv             : dZ/dv at (s,u,v)                           --
!----------------------------------------------------------------------
      CHARACTER(len=80)           :: name_coils_dot
      CHARACTER(len=120)           :: name_diagnostic_dot
      LOGICAL                     :: lstell_sym
      REAL(rprec)                 :: rmin, rmax, zmin, zmax
      INTEGER(iprec)              :: ir, jz, kp, kp_store
!      INTEGER(iprec)              :: n_field_periods
!
!      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:,:) :: pl_response
!      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:,:) :: pl_respsuv
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: rgrid
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: zgrid
!      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: pgrid
      REAL(rprec)                 :: delr, delz, delp
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: sgrid
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: ugrid
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: vgrid
      REAL(rprec)                 :: dels, delu, delv
      INTEGER(iprec)              :: ju, kv
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: rsuv
!      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: psuv
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: zsuv
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: gsuv
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: currusuv
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: currvsuv
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: bsubu
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: bsubv
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: rusuv
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: zusuv
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: rvsuv
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: zvsuv
!----------------------------------------------------------------------
!--   signals           : signals                                    --
!--     cal             : calculated values                          --
!--     cext            : calculated external values                 --
!----------------------------------------------------------------------
      TYPE signals
        REAL(rprec)       :: cal
        REAL(rprec)       :: cext
      END TYPE signals
!----------------------------------------------------------------------
      TYPE(signals), ALLOCATABLE                   :: signal_diag(:)
      CHARACTER(len=30), ALLOCATABLE, DIMENSION(:) :: shortnames

      END MODULE v3post_rfun
