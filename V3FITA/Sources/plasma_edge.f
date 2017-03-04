!     SPH: INTEGER(iprec) -> INTEGER
!============================================================================
      MODULE plasma_edge
!--------------------------------------------------------------------------
!  Module to define data TYPEs that contain information about the outer edge of
!  the plasma for use in the Faraday Rotation module faraday_mod.  It
!     is NOT used in ANY subroutines outside of faraday_mod.
!-------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      USE v3_utilities       ! assert_eq is here
      IMPLICIT NONE



!=============================================================================
      TYPE pv_edge
!-------------------------------------------------------------------------------
! declares variables relating to the plasma/vacuum boundary points along the beam path
!
!   n_edge           :  number of plasma/vacuum boundary points (plus path endpts)
!                        (note that an endpt can be a boundary edge AND and an endpt)
!   etype(n_edge)    :  type of boundary edge ( 1 = vac/plasma, 2 = plasma/vac,
!                      -11 = path endpt INSIDE plasma,  -21 = path endpt OUTSIDE plasma
!      ^ Use as an ALLOCATABLE ARRAY 
!   edist(n_edge)    :  type of boundary edge & beam start point (q0)
!      ^ Use as an ALLOCATABLE ARRAY 
!   xflx(3,n_edge)   :  flux coordinate vector at each boundary point
!      ^ Use as an ALLOCATABLE ARRAY 
!   xcar(3,n_edge)   :  Cartesian vector at each boundary point
!      ^ Use as an ALLOCATABLE ARRAY 
!
!  NOTES:  1) edge arrays include beampath endpoints, even if such points are not actually
!              plasma/vacuum boundaries
!
!          2) The edgepoints are defined to be the last points along the beampath PRIOR
!             to crossing a plasma/vac or vac/plasma boundary.  (ie if the beampath was
!             from left to right, the edgepoints would be just LEFT of the boundary)
!    
!------------------------------------------------------------------------------
      INTEGER                         :: n_edge
      INTEGER, DIMENSION(:), POINTER  :: etype => null()
      REAL(rprec), DIMENSION(:),   POINTER   :: edist => null()
      REAL(rprec), DIMENSION(:,:), POINTER   :: xflx  => null()
      REAL(rprec), DIMENSION(:,:), POINTER   :: xcar  => null()
      END TYPE pv_edge



!----------------------------------------------------------------------------
!  subroutines for use in module plasma_edge
!---------------------------------------------------------------------------
      CONTAINS


!===============================================================================
      SUBROUTINE pv_edge_construct(this, n_edge)
!----------------------------------------------------------------------------
!  allocates memory for a pv_edge variable.  Array values are NOT filled.
!---------------------------------------------------------------------------

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (pv_edge), INTENT (inout)     :: this
      INTEGER, INTENT(in)                :: n_edge

!  Declare local variables
      INTEGER :: nd1, nd2, ier1
      INTEGER, DIMENSION(3)   :: dims
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'pv_edge_construct: '

!  Start of executable code

!      WRITE(*,*) ' Executing pv_edge_construct'

!....if "this" already exists, deallocate the memory and start over....!
!      CALL pv_edge_destroy(this)

      this % n_edge = n_edge

      ALLOCATE(this % etype(n_edge),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc etype')

      ALLOCATE(this % edist(n_edge),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc edist')

      ALLOCATE(this % xflx(3,n_edge),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc xflx')

      ALLOCATE(this % xcar(3,n_edge),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc xcar')

      
      END SUBROUTINE pv_edge_construct



!==============================================================================
      SUBROUTINE pv_edge_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (pv_edge), INTENT(inout) :: this

!  Declare local variables
      INTEGER :: ier1
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'pv_edge_destroy: '

!  Start of executable code
!  Get rid of all components
      
      this % n_edge = 0


      IF ( ASSOCIATED(this % etype) ) THEN
         DEALLOCATE(this % etype,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc etype')
      ENDIF

      IF ( ASSOCIATED(this % edist) ) THEN
         DEALLOCATE(this % edist,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc edist')
      ENDIF

      IF ( ASSOCIATED(this % xflx) ) THEN
         DEALLOCATE(this % xflx,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc xflx')
      ENDIF

      IF ( ASSOCIATED(this % xcar) ) THEN
         DEALLOCATE(this % xcar,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc xcar')
      ENDIF

      END SUBROUTINE pv_edge_destroy

      END MODULE plasma_edge
