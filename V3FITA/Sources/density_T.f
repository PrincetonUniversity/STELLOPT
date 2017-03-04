!*******************************************************************************
!  File density_T.f
!  Contains module density_T
!  Defines derived-types: e_density

!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  It deals with the interface with the equilibrium code.
!  Right now (6/2004), the VMEC code is the only equilibrium code.


!*******************************************************************************
!  MODULE density_T
! SECTION I.      VARIABLE DECLARATIONS
! SECTION II.     DERIVED TYPE DECLARATIONS
! SECTION III.    INTERFACE BLOCKS
! SECTION IV.     CONSTRUCTION SUBROUTINES
! SECTION V.      DESTRUCTION SUBROUTINES
! SECTION VI.     ASSIGNMENT SUBROUTINES
! SECTION VII.    AUX DEFINITION SUBROUTINES
! SECTION VIII.   OTHER SUBROUTINES
!*******************************************************************************
      MODULE density_T

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!-------------------------------------------------------------------------------
      USE stel_kinds

!-------------------------------------------------------------------------------
!  Frequently used mathematical constants, lots of extra precision.
!-------------------------------------------------------------------------------
      USE stel_constants

!-------------------------------------------------------------------------------
!  Use Statements for  V3 Utilities
!-------------------------------------------------------------------------------
      USE v3_utilities
      
      IMPLICIT NONE

!*******************************************************************************
! SECTION II. DERIVED TYPE (STRUCTURE) DECLARATIONS
!
!   Derived Type to carry information about the parameters that define the density profile
!
!*******************************************************************************
!
!  
!-------------------------------------------------------------------------------
!  Declare type e_density
!    Parameters refer to an electron density profile of form:
!            ne = ne_max * ( 1 - s^tau)^kappa + ne_ambient   
!                  where s = the radial flux coordinate
!
!  ne_max:       maximum e- density (assumed at magnetic axis).  Units = m^(-3)
!  ne_ambient:   ambient e- density in the vicinity surrounding the plasma  [ m^(-3) ]
!  tau:          exponent on the radial flux (dimensionless)
!  kappa:        exponent on the ( 1 - s ^tau) term.   (dimensionless)
!
!-------------------------------------------------------------------------------
      TYPE e_density
!  Variables for a density profile

         REAL(rprec) ::  ne_max
         REAL(rprec) ::  ne_ambient
         REAL(rprec) ::  tau
         REAL(rprec) ::  kappa

      END TYPE e_density


!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************


      CONTAINS
!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct an e_density
!
!-------------------------------------------------------------------------------
      SUBROUTINE e_density_construct(this, ne_max, ne_ambient,                 &
     &                               tau, kappa)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (e_density), INTENT (inout)   :: this
      REAL(rprec), INTENT(in)               :: ne_max
      REAL(rprec), INTENT(in)               :: ne_ambient
      REAL(rprec), INTENT(in)               :: tau
      REAL(rprec), INTENT(in)               :: kappa

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'e_density_construct: '

!  Start of executable code

!      WRITE(*,*) ' now Executing ', sub_name

      this % ne_max     = ne_max
      this % ne_ambient = ne_ambient
      this % tau        = tau
      this % kappa      = kappa
      
      END SUBROUTINE e_density_construct


!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Destroy an e_density
!-------------------------------------------------------------------------------
      SUBROUTINE e_density_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (e_density), INTENT(inout) :: this

!  Start of executable code
!  Get rid of all components

      this % ne_max     = zero
      this % ne_ambient = zero
      this % tau        = zero
      this % kappa      = zero

      END SUBROUTINE e_density_destroy


!*******************************************************************************
! SECTION VI. ASSIGNMENT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for e_density - Default OK as no Pointers
!-------------------------------------------------------------------------------



!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************

           
      END MODULE density_T
