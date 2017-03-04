!*******************************************************************************
!  File model_T.f
!  Contains module model_T
!  Defines derived-types: model

!*******************************************************************************
!  MODULE model_T
!    (Signal Type Definition, for the V3FIT code)
! SECTION I.    VARIABLE DECLARATIONS
! SECTION II.   DERIVED-TYPE DECLARATIONS
! SECTION III.  INTERFACE BLOCKS
! SECTION IV.   CONSTRUCTION SUBROUTINES
! SECTION V.    DESTRUCTION SUBROUTINES
! SECTION VI.   ASSIGNMENT SUBROUTINES

! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE model_T

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!  Frequently used mathematical constants, lots of extra precision.
!-------------------------------------------------------------------------------
      USE stel_kinds, only : rprec, iprec, cprec
      USE stel_constants, only : pi, twopi, one, zero
     
!-------------------------------------------------------------------------------
!  Use Statements for other structures, V3 Utilities
!-------------------------------------------------------------------------------
      USE bsc_T
      USE v3_utilities
      USE eq_T
      USE density_T

!-------------------------------------------------------------------------------
!  Implicit None comes after USE statements, before other declarations
!-------------------------------------------------------------------------------
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Make type declarations and constants Private, so there are no conflicts.
!-------------------------------------------------------------------------------
      PRIVATE rprec, iprec, cprec, pi, twopi, one, zero

!-------------------------------------------------------------------------------
!  Lengths of Character Variables
!-------------------------------------------------------------------------------
      INTEGER(iprec), PARAMETER, PRIVATE :: type_len=10      
      INTEGER(iprec), PARAMETER, PRIVATE :: sn_len=30      
      INTEGER(iprec), PARAMETER, PRIVATE :: ln_len=80      
      INTEGER(iprec), PARAMETER, PRIVATE :: units_len=30      

!*******************************************************************************
! SECTION II. DERIVED-TYPE DECLARATIONS
! 1)   Model
!         model  
!     A model will contain an equilibrium state, and any other information 
!     that is needed to compute model signals. For example, if the experiment
!     has interferometry and polarimetry diagnsotics, then there will need
!     to be a model density profile.
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Declare type model
!         ---- Equilibrium State ----
!    eqstate         Type eq_state
!
!         ---- Electron Density  -----
!    density          Type e_density
!
!         ----  Other Information ----
!
!-------------------------------------------------------------------------------
      TYPE model
         TYPE (eq_state)          :: eqstate
         TYPE (e_density)         :: density
      END TYPE model

!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for structures
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Generic construct
!-------------------------------------------------------------------------------
 !-------------------------------------------------------------------------------
!  Generic destroy
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Generic write
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Interface block for testing goes here. 
!-------------------------------------------------------------------------------

      CONTAINS
!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a model
!-------------------------------------------------------------------------------
      SUBROUTINE model_construct(this,eqstate, density)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (model), INTENT(inout)              :: this
      TYPE (eq_state), INTENT(in), OPTIONAL    :: eqstate
      TYPE (e_density), INTENT(in), OPTIONAL   :: density

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'model_construct: '

!  Start of executable code

!  Copy the equilibrium state
      IF (PRESENT(eqstate)) THEN
         this % eqstate = eqstate
      ENDIF

!  Copy the density information
      IF (PRESENT(density)) THEN
         this % density = density
      ENDIF
      
      END SUBROUTINE model_construct

      
!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Destroy a model
!-------------------------------------------------------------------------------
      SUBROUTINE model_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (model), INTENT(inout) :: this

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'model_destroy: '

!  Start of executable code

!  As of 12-14-2004, there isn't an eq_destroy for a state.
!      CALL eq_destroy(this % eqstate)

!.........destroy the density info..........................!
      CALL e_density_destroy( this % density)

      END SUBROUTINE model_destroy

!*******************************************************************************
! SECTION VI. ASSIGNMENT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for model
!-------------------------------------------------------------------------------
      SUBROUTINE model_assign(left,right)


      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (model), INTENT (inout) :: left
      TYPE (model), INTENT (in) :: right
      
!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'model_assign: '
         
!  Start of executable code


      left % eqstate = right % eqstate
      left % density = right % density
         
      END SUBROUTINE model_assign

!*******************************************************************************
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!  
!  JDH 12-14-04. 
!     First version, modified signal_T
!  JMS 7-24-07.
!     added density to the model

      END MODULE model_T
