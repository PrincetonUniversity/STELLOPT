!*******************************************************************************
!>  @file magnetic.f
!>  @brief Contains module @ref integration_path
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Module is part of the LIBSTELL. This modules contains code to define and
!>  integrate along an arbitray path.
!*******************************************************************************

      MODULE integration_path
      USE stel_kinds
      USE mpi_inc

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) vertex
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  A single point in space defined by an z, y, z coordinate. A vertex is
!>  structured as a singly linked list.
!-------------------------------------------------------------------------------
      TYPE vertex
!>  Position in cartiesian coordinates.
         REAL (rprec), DIMENSION(3) :: position
!>  Reference to the next vertex.
         TYPE (vertex), POINTER     :: next => null()
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for checking the results of the unit tests
!-------------------------------------------------------------------------------
      INTERFACE check
         MODULE PROCEDURE check_log, check_real, check_int
      END INTERFACE

      PUBLIC :: path_construct, path_append_vertex, path_destruct,             &
     &          path_test, path_integrate
      PRIVATE :: check, check_real, check_log, check_int,                      &
     &           integrate, search, test_function

      CONTAINS

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a single @ref vertex.
!>
!>  Allocates memory and initializes a @ref vertex object.
!>
!>  @param[in] position Cartesian position of the vertex object.
!>  @returns A pointer to a constructed @ref vertex object.
!-------------------------------------------------------------------------------
      FUNCTION path_construct(position)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vertex), POINTER                 :: path_construct
      REAL (rprec), DIMENSION(3), INTENT(in) :: position

!  Start of executable code
      ALLOCATE(path_construct)

      path_construct%position = position

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref vertex object.
!>
!>  Deallocates memory and uninitializes a @ref vertex object. This recursively
!>  deconstructed the next vertex until the last in the linked list is found.
!>
!>  @param[inout] this A @ref vertex instance.
!-------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE path_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vertex), POINTER :: this

!  Start of executable code
      IF (ASSOCIATED(this%next)) THEN
         CALL path_destruct(this%next)
         this%next => null()
      END IF

      DEALLOCATE(this)
      this => null()

      END SUBROUTINE

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Append a @ref vertex to a path.
!>
!>  Recursively runs through the next vertex to find the last vertex. Once the
!>  last vertex is found, a new vertex is allocated and appended to the path.
!>  This allow works as a constructer and allocates the first vertex if needed.
!>
!>  @param[inout] this     Vertex to append path to.
!>  @param[in]    position Cartesian position of the vertex object.
!-------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE path_append_vertex(this, position)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (vertex), POINTER                 :: this
      REAL (rprec), DIMENSION(3), INTENT(in) :: position

!  Start of executable code
      IF (ASSOCIATED(this)) THEN
         IF (ASSOCIATED(this%next)) THEN
            CALL path_append_vertex(this%next, position)
         ELSE
            this%next => path_construct(position)
         END IF
      ELSE
         this => path_construct(position)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Integrate along the path.
!>
!>  Recursively runs through the next vertex to find the last vertex. Once the
!>  last vertex is found, integrate alone that path. The integrand is proveded
!>  by means of call back function.
!>
!>  @param[inout] this                 Starting vertex to integrate to the end.
!>  @param        integration_function Function pointer that defines the
!>                                     integrand.
!>  @param[in]    context              Generic object that contains data for the
!>                                     integration function.
!>  @returns The total integrated path to the end.
!-------------------------------------------------------------------------------
      RECURSIVE FUNCTION path_integrate(this, integration_function,            &
     &                                  context) RESULT(total)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                  :: total
      TYPE (vertex), INTENT(inout)  :: this
      CHARACTER (len=1), INTENT(in) :: context(:)
      INTERFACE
         FUNCTION integration_function(context, xcart, dxcart,                 &
     &                                 length, dx)
         USE stel_kinds
         REAL (rprec) :: integration_function
         CHARACTER (len=1), INTENT(in)          :: context(:)
         REAL (rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
         REAL (rprec), INTENT(in)               :: length, dx
         END FUNCTION
      END INTERFACE

!  Start of executable code
      If (ASSOCIATED(this%next)) THEN
         total = path_integrate(this%next, integration_function,               &
     &                          context)
         total = total + integrate(context, this, this%next,                   &
     &                             integration_function)
      ELSE
         total = 0.0
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Search along the path.
!>
!>  Recursively runs through the next vertex to find the last vertex. Once the
!>  last vertex is found, a new vertex is allocated and appended to the path.
!>  The integrand is proveded by means of call back function.
!>
!>  @param[inout] this            Starting vertex to begin search.
!>  @param        search_function Function pointer that defines the search
!>                                criteria.
!>  @param[in]    context         Generic object that contains data for the
!>                                integration function.
!>  @param[out]   found           Signals if the condition was met.
!>  @returns The vertex position along the path where the search condition was
!>           found.
!-------------------------------------------------------------------------------
      RECURSIVE FUNCTION path_search(this, search_function, context,           &
     &                               found) RESULT(xcart)

      REAL (rprec), DIMENSION(3)    :: xcart
      TYPE (vertex), INTENT(inout)  :: this
      CHARACTER (len=1), INTENT(in) :: context(:)
      LOGICAL, INTENT(out)          :: found
      INTERFACE
         FUNCTION search_function(context, xcart1, xcart2)
         USE stel_kinds
         LOGICAL                                :: search_function
         CHARACTER (len=1), INTENT(in)          :: context(:)
         REAL (rprec), DIMENSION(3), INTENT(in) :: xcart1
         REAL (rprec), DIMENSION(3), INTENT(in) :: xcart2
         END FUNCTION
      END INTERFACE

!  Start of executable code
      found = .false.

      IF (ASSOCIATED(this%next)) THEN
         found = search(context, this, this%next, search_function,             &
     &                  xcart)
         IF (.not.found) THEN
            xcart = path_search(this%next, search_function, context,           &
     &                          found)
         END IF
      END IF

      END FUNCTION

!*******************************************************************************
!  PRIVATE
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Line integrate between to points.
!>
!>  This divids the straight line path defined by two vertices and line
!>  integrates the function. The integrand is proveded by means of call back
!>  function.
!>
!>  @param[in] context              Generic object that contains data for the
!>                                  integration function.
!>  @param[in] vertex1              Starting point.
!>  @param[in] vertex2              Ending point.
!>  @param     integration_function Function pointer that defines the
!>                                  integrand.
!>  @returns The path integrated value between the vertex1 and vertex2.
!-------------------------------------------------------------------------------
      FUNCTION integrate(context, vertex1, vertex2,                            &
     &                   integration_function)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                 :: integrate
      CHARACTER(len=1), INTENT(in) :: context(:)
      TYPE (vertex), INTENT(in)    :: vertex1
      TYPE (vertex), INTENT(in)    :: vertex2
      INTERFACE
         FUNCTION integration_function(context, xcart, dxcart,                 &
     &                                 length, dx)
         USE stel_kinds
         REAL (rprec) :: integration_function
         CHARACTER (len=1), INTENT(in)          :: context(:)
         REAL (rprec), DIMENSION(3), INTENT(in) :: xcart, dxcart
         REAL (rprec), INTENT(in)               :: length, dx
         END FUNCTION
      END INTERFACE

!  local variables
      REAL (rprec), DIMENSION(3) :: xcart
      REAL (rprec), DIMENSION(3) :: dxcart
      REAL (rprec)               :: length
      REAL (rprec)               :: dx
      INTEGER                    :: i, nsteps

!  local parameters
      REAL(rprec), PARAMETER     :: dxOptimal = 0.0025

!  Start of executable code
!  Determine the number of integration steps to take by dividing the path length
!  by the step length and rounding to the nearest integer.
      dxcart = vertex2%position - vertex1%position
      length = SQRT(DOT_PRODUCT(dxcart, dxcart))
      nsteps = INT(length/dxOptimal)

!  Choose the actual step size.
      dxcart = dxcart/nsteps
      dx = SQRT(DOT_PRODUCT(dxcart, dxcart))

!  Integrate the length in addition.
      integrate = 0.0
      length = 0.0
      xcart = vertex1%position
      DO i = 1, nsteps
         xcart = xcart + dxcart
         length = length + dx
         integrate = integrate                                                 &
     &             + integration_function(context, xcart, dxcart,              &
     &                                    length, dx)
      ENDDO

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Search line between to points.
!>
!>  This divides the straight line path defined by two vertices and searched for
!>  a condition. The search criteria is proveded by means of call back function.
!>
!>  @param[in] context          Generic object that contains data for the search
!>                              function.
!>  @param[in]  vertex1         Starting point.
!>  @param[in]  vertex2         Ending point.
!>  @param      search_function Function pointer that defines the search
!>                              criteria.
!>  @param[out] xcart           Point where the search criteria was found.
!>  @returns True if the criteria was met between vertex1 and vertex2.
!-------------------------------------------------------------------------------
      FUNCTION search(context, vertex1, vertex2, search_function, xcart)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                                 :: search
      CHARACTER(len=1), INTENT(in)            :: context(:)
      TYPE (vertex), INTENT(in)               :: vertex1
      TYPE (vertex), INTENT(in)               :: vertex2
      REAL (rprec), DIMENSION(3), INTENT(out) :: xcart
      INTERFACE
         FUNCTION search_function(context, xcart1, xcart2)
         USE stel_kinds
         LOGICAL                                :: search_function
         CHARACTER (len=1), INTENT(in)          :: context(:)
         REAL (rprec), DIMENSION(3), INTENT(in) :: xcart1
         REAL (rprec), DIMENSION(3), INTENT(in) :: xcart2
         END FUNCTION
      END INTERFACE

!  local variables
      REAL (rprec), DIMENSION(3) :: dxcart
      REAL (rprec)               :: dx
      INTEGER                    :: i, nsteps

!  local parameters
      REAL(rprec), PARAMETER     :: dxCourse = 0.01
      REAL(rprec), PARAMETER     :: dxFine = 1.0E-20

!  Start of executable code
!  Determine the number of integration steps to take by dividing the path length
!  by the step length and rounding to the nearest integer.
      dxcart = vertex2%position - vertex1%position
      nsteps = INT(SQRT(DOT_PRODUCT(dxcart, dxcart))/dxCourse)

!  Choose the actual step size.
      dxcart = dxcart/nsteps

      search = .false.
      xcart = vertex1%position

!  Linearly search the line until an interval containing the point is detected.
      DO i = 1, nsteps
         search = search_function(context, xcart, xcart + dxcart)
         IF (search) THEN

!  Found an interval. Bisect the interval until the length is machine precision.
            DO WHILE (SQRT(DOT_PRODUCT(dxcart, dxcart)) .gt. dxFine)
               dxcart = dxcart/2.0
               IF (.not.search_function(context, xcart,                        &
     &                                  xcart + dxcart)) THEN
                  xcart = xcart + dxcart
               END IF
            END DO

            xcart = xcart + dxcart/2.0
            RETURN

         END IF

         xcart = xcart + dxcart
      END DO

      END FUNCTION

!*******************************************************************************
!  UNIT TESTS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Path unit test function.
!>
!>  This runs the associated unit tests and returns the result.
!>
!>  @returns True if the tests pass and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION path_test()

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: path_test

!  local variables
      REAL (rprec)                  :: result
      TYPE (vertex), POINTER        :: test_path => null()
      CHARACTER(len=1), ALLOCATABLE :: context(:)

!  Start of executable code
!  Test to make sure the vertices begin in an unallocated state.
      path_test = check(.false., ASSOCIATED(test_path), 1, "ASSOCIATED")
      IF (.not.path_test) RETURN

!  Test to make sure first vertex is created.
      CALL path_append_vertex(test_path,                                       &
     &                        (/ 1.0_rprec, 2.0_rprec, 3.0_rprec /))
      path_test = check(.true., ASSOCIATED(test_path), 1,                      &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(1.0_rprec, test_path%position(1), 2,                   &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(2.0_rprec, test_path%position(2), 3,                   &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(3.0_rprec, test_path%position(3), 4,                   &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN

!  Test to make sure second vertex is appended.
      CALL path_append_vertex(test_path,                                       &
     &                        (/ 4.0_rprec, 5.0_rprec, 6.0_rprec /))
      path_test = check(.true., ASSOCIATED(test_path%next), 5,                 &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(1.0_rprec, test_path%position(1), 6,                   &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(2.0_rprec, test_path%position(2), 7,                   &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(3.0_rprec, test_path%position(3), 8,                   &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(4.0_rprec, test_path%next%position(1), 9,              &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(5.0_rprec, test_path%next%position(2), 10,             &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(6.0_rprec, test_path%next%position(3), 11,             &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN

!  Test to make sure third vertex is appended.
      CALL path_append_vertex(test_path,                                       &
     &                        (/ 7.0_rprec, 8.0_rprec, 9.0_rprec /))
      path_test = check(.true., ASSOCIATED(test_path%next%next),               &
     &                  12, "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(7.0_rprec, test_path%next%next%position(1), 12,        &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(8.0_rprec, test_path%next%next%position(2), 13,        &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN
      path_test = check(9.0_rprec, test_path%next%next%position(3), 14,        &
     &                  "path_append_vertex")
      IF (.not.path_test) RETURN

!  Test to make sure path object is destroyed.
      CALL path_destruct(test_path)
      path_test = check(.false., ASSOCIATED(test_path), 1,                     &
     &                  "path_destruct")
      IF (.not.path_test) RETURN

!  Test path integration.
      CALL path_append_vertex(test_path,                                       &
     &                        (/ 1.0_rprec, 0.0_rprec, 0.0_rprec /))
      CALL path_append_vertex(test_path,                                       &
     &                        (/ 0.0_rprec, 0.0_rprec, 0.0_rprec /))
      result = path_integrate(test_path, test_function, context)
      path_test = check(400.0_rprec, result, 1, "path_integrate")
      CALL path_destruct(test_path)

      END FUNCTION

!*******************************************************************************
!  CHECK FUNCTIONS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Check a logical value.
!>
!>  Checks that the expected value matches the recieved. Otherwise report an
!>  error.
!>
!>  @param[in] expected The known value.
!>  @param[in] received The known test value.
!>  @param[in] testNum  The number of the test.
!>  @param[in] name     The name of the test.
!>  @returns True if the check passes and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION check_log(expected, received, testNum, name)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: check_log
      LOGICAL, INTENT(in)           :: expected
      LOGICAL, INTENT(in)           :: received
      INTEGER, INTENT(in)           :: testNum
      CHARACTER (LEN=*), INTENT(in) :: name

!  Start of executable code
      check_log = expected .eqv. received
      IF (.not.check_log) THEN
         write(*,*) "integration_path.f: ", name, " test", testNum,            &
     &              "failed."
         write(*,*) "Expected", expected, "Received", received
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Check a real value.
!>
!>  Checks that the expected value matches the recieved. Otherwise report an
!>  error.
!>
!>  @param[in] expected The known value.
!>  @param[in] received The known test value.
!>  @param[in] testNum  The number of the test.
!>  @param[in] name     The name of the test.
!>  @returns True if the check passes and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION check_real(expected, received, testNum, name)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: check_real
      REAL(rprec), INTENT(in)       :: expected
      REAL(rprec), INTENT(in)       :: received
      INTEGER, INTENT(in)           :: testNum
      CHARACTER (LEN=*), INTENT(in) :: name

!  Start of executable code
      check_real = expected .eq. received
      IF (.not.check_real) THEN
         write(*,*) "integration_path.f: ", name, " test", testNum,            &
     &              "failed."
         write(*,*) "Expected", expected, "Received", received
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Check an integer value.
!>
!>  Checks that the expected value matches the recieved. Otherwise report an
!>  error.
!>
!>  @param[in] expected The known value.
!>  @param[in] received The known test value.
!>  @param[in] testNum  The number of the test.
!>  @param[in] name     The name of the test.
!>  @returns True if the check passes and false otherwise.
!-------------------------------------------------------------------------------
      FUNCTION check_int(expected, received, testNum, name)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: check_int
      INTEGER, INTENT(in)           :: expected
      INTEGER, INTENT(in)           :: received
      INTEGER, INTENT(in)           :: testNum
      CHARACTER (LEN=*), INTENT(in) :: name

!  Start of executable code
      check_int = expected .eq. received
      IF (.not.check_int) THEN
         write(*,*) "integration_path.f: ", name, " test", testNum,            &
     &              "failed."
         write(*,*) "Expected", expected, "Received", received
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Call back function to test the integration.
!>
!>  This defined a function of f(x) = 1.
!>
!>  @param[in] context Generic object that contains data for the integration
!>                     function.
!>  @param[in] xcart   Current point in the integration.
!>  @param[in] dxcart  Vector direction of the current integration path.
!>  @param[in] length  Length along the current integration.
!>  @param[in] dx      Scaler length of the current integration path.
!>  @returns One
!-------------------------------------------------------------------------------
      FUNCTION test_function(context, xcart, dxcart, length, dx)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec)                           :: test_function
      CHARACTER (len=1), INTENT(in)          :: context(:)
      REAL (rprec), DIMENSION(3), INTENT(in) :: xcart
      REAL (rprec), DIMENSION(3), INTENT(in) :: dxcart
      REAL (rprec), INTENT(in)               :: length
      REAL (rprec), INTENT(in)               :: dx

!  Start of executable code
      test_function = 1.0

      END FUNCTION

      END MODULE
