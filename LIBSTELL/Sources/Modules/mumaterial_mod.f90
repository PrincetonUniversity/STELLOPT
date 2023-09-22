!-----------------------------------------------------------------------
!     Module:        mumaterial_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov), Björn Hamstra
!     Date:          September 2023
!     Description:   This module is designed to help calculate the
!                    magnetic field arrising from ferromagnetic
!                    material
!-----------------------------------------------------------------------
      MODULE mumaterial_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE safe_open_mod
      USE TileNComponents
      USE IterateMagnetSolution
      IMPLICIT NONE

!-----------------------------------------------------------------------
!     Types    
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     Module Variables
!           nvertex:          Number of vertices
!           ntet:             Number of tetrahedrons
!           nstate:           Number of state functions
!           tiles:            Magnetic tiles
!           stateFunction:    Array of state functions
!           maxErr:           Max allowed error for convergence in MagTense iteration
!           maxIter           Max allowed number of MagTense iterations
!           temp:             Temperature of magnetic material [K]
!           vertex:           Vertices [m] (nvertex,3)
!           tet:              Tetrahedron (ntet,4) index into vertex
!           state_dex:        State function for each tetrahedron
!           state_type:       Type of state function (nstate) (1 or 3)
!                                (type 2 not yet supported)
!           constant_mu:      Mu values for constant permeability (nstate)
!           M:                Magnet M Vector for hard magnet (nstate,3)
!
!-----------------------------------------------------------------------
      INTEGER, PRIVATE                    :: nvertex, ntet, nstate
      TYPE(MagTile), PRIVATE, ALLOCATABLE              :: tiles(:)
      TYPE(MagStateFunction), PRIVATE, ALLOCATABLE     :: stateFunction(:)
      DOUBLE PRECISION, PRIVATE           :: maxErr, temp
      INTEGER, PRIVATE                    :: maxIter


      DOUBLE PRECISION, POINTER, PRIVATE :: vertex(:,:)
      INTEGER, POINTER, PRIVATE :: tet(:,:)
      INTEGER, POINTER, PRIVATE :: state_dex(:)
      INTEGER, POINTER, PRIVATE :: state_type(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: constant_mu(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: Mag(:,:)
      INTEGER, PRIVATE            :: win_vertex, win_tet,  &
                                     win_state_dex, win_state_type, &
                                     win_constant_mu, win_m

      CHARACTER(LEN=256), PRIVATE :: machine_string
      CHARACTER(LEN=256), PRIVATE :: date


      INTEGER, PRIVATE                    :: mystart, myend, mydelta
      INTEGER, PRIVATE                    :: shar_rank, shar_size, shar_comm




      
!-----------------------------------------------------------------------
!     Subroutines
!         mumaterial_setd: Overrides default values
!         mumaterial_load: Loads a magnetic material file
!         mumaterial_init: Calculates magnetization of material
!         mumaterial_getb_scalar: Calculates magnetic field at a point in space
!         mumaterial_getb_vector: Calculates magnetic field for multiple points in space
!         mumaterial_output: Output tiles, H-field and points to file
!-----------------------------------------------------------------------
!     Functions
!-----------------------------------------------------------------------
      INTERFACE mumaterial_getb
            MODULE PROCEDURE mumaterial_getb_scalar, mumaterial_getb_vector
      END INTERFACE
      CONTAINS
      




      SUBROUTINE mumaterial_setd(mE, mI, T)
      !-----------------------------------------------------------------------
      ! mumaterial_setd: Overrides default values
      !-----------------------------------------------------------------------
      ! param[in]: mE. New maxErr: max error for MagTense convergence
      ! param[in]: mI. New maxIter: max amount of MagTense iterations
      ! param[in]: T. New temp: temperature of magnetic material in MagTense
      !-----------------------------------------------------------------------
      DOUBLE PRECISION, INTENT(in) :: mE, T
      INTEGER, INTENT(in) :: mI

      WRITE(*,*) "Setting max Error, max Iterations and Temperature to: ", mE, mI, T

      maxErr = mE
      maxIter = mI
      temp = T

      RETURN
      END SUBROUTINE mumaterial_setd





      SUBROUTINE mumaterial_load(filename,istat,comm)
      !-----------------------------------------------------------------------
      ! mumaterial_load: Loads magnetic material file.
      !-----------------------------------------------------------------------
      ! param[in]: filename. The file name to load in
      ! param[in, out]: istat. Integer that shows error if != 0
      ! param[in]: verb: Verbosity. True or false
      ! param[in, out]: comm. MPI communicator, handles shared memory
      !-----------------------------------------------------------------------
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER, INTENT(inout)       :: istat
      INTEGER, INTENT(inout), OPTIONAL :: comm
      LOGICAL :: shared
      INTEGER :: iunit ,ik, i, j
      
      shar_rank = 0; shar_size = 1;
      ! initialize MPI
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_COMM_SPLIT_TYPE(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shar_comm, istat)
         CALL MPI_COMM_RANK( shar_comm, shar_rank, istat )
         CALL MPI_COMM_SIZE( shar_comm, shar_size, istat)
      END IF
#endif
      ! open file, return if fails
      iunit = 327; istat = 0
      CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
      IF (istat/= 0) RETURN
      ! read info
      IF (shar_rank == 0) THEN
         READ(iunit,'(A)') machine_string
         READ(iunit,'(A)') date
         READ(iunit,*) nvertex, ntet, nstate
      END IF
      ! Broadcast info to MPI and allocate vertex and face info
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_Bcast(nvertex,1,MPI_INTEGER,0,shar_comm,istat)
         CALL MPI_Bcast(ntet,1,MPI_INTEGER,0,shar_comm,istat)
         CALL MPI_Bcast(nstate,1,MPI_INTEGER,0,shar_comm,istat)
         CALL mpialloc_2d_dbl(vertex,nvertex,3,shar_rank,0,shar_comm,win_vertex)
         CALL mpialloc_2d_int(tet,ntet,4,shar_rank,0,shar_comm,win_tet)
         CALL mpialloc_1d_int(state_dex,ntet,shar_rank,0,shar_comm,win_state_dex)
         CALL mpialloc_1d_int(state_type,nstate,shar_rank,0,shar_comm,win_state_type)
         CALL mpialloc_1d_dbl(constant_mu,nstate,shar_rank,0,shar_comm,win_constant_mu)
         CALL mpialloc_2d_dbl(Mag,nstate,3,shar_rank,0,shar_comm,win_vertex)
         shared = .true.
      ELSE
#endif
         ! if no MPI, allocate everything on one node
         ALLOCATE(vertex(nvertex,3),tet(ntet,4),state_dex(ntet), &
                  state_type(nstate),constant_mu(nstate),Mag(nstate,3),STAT=istat)
         shared = .false.
#if defined(MPI_OPT)
      END IF
#endif
      WRITE(*,*) "Reading tile data"
      ! read in the mesh
      IF (istat/=0) RETURN
      IF (shar_rank == 0) THEN
         DO ik = 1, nvertex
            READ(iunit,*) vertex(ik,1),vertex(ik,2),vertex(ik,3)
         END DO
         DO ik = 1, ntet
            READ(iunit,*) tet(ik,1),tet(ik,2),tet(ik,3),tet(ik,4),state_dex(ik)
         END DO
         ! This next section will get more complicated in the future for now
         ! we just handle the soft constant and hard magnets
         DO ik = 1, nstate
            READ(iunit,*) state_type(ik)
            IF (state_type(ik) == 1) THEN
               READ(iunit,*) constant_mu(ik)
            ELSEIF (state_type(ik) == 2) THEN
               PRINT *, 'STATE_TYPE == 2 (soft magnet) not yet supported.'
                  ! old function for loading state function, may be applicable, adapted from MagTense Standalone IO
                  ! subroutine loadStateFunctionFortran(stateFunction)
                  !       integer :: i
                  !       type(MagStateFunction), dimension(1), intent(out) :: stateFunction
                  !       OPEN(11, file='stateFunction.dat')
                  !       read(11,*) stateFunction(1)%nT,stateFunction(1)%nH
                        
                  !       allocate( stateFunction(1)%M(stateFunction(1)%nT,stateFunction(1)%nH) )
                  !       allocate( stateFunction(1)%T(stateFunction(1)%nT) )
                  !       allocate( stateFunction(1)%H(stateFunction(1)%nH) )
                        
                  !       stateFunction(1)%T(1) = 300
                  !       DO i=1,stateFunction(1)%nH
                  !       read(11,*) stateFunction(1)%H(i),stateFunction(1)%M(1,i)        
                  !       END DO
                        
                  !       close (11)
                  ! end subroutine loadStateFunctionFortran
            ELSEIF (state_type(ik) == 3) THEN
               READ(iunit,*) constant_mu(ik)
            ELSE
               PRINT *, '!!! UNKNOWN STATE_TYPE == ',state_type(ik)
            END IF
         END DO
      END IF

#if defined(MPI_OPT)
      IF (PRESENT(comm)) CALL MPI_BARRIER(shar_comm,istat)
#endif
      ! close file
      CLOSE(iunit)

      WRITE(*,*) "Finished reading tile data"

      ! set default values
      call MUMATERIAL_SETD(1.0d-5, 100, 300.d0)



      RETURN
      END SUBROUTINE mumaterial_load





      SUBROUTINE mumaterial_info(iunit)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iunit
      WRITE(iunit,*) 'NAME: ',TRIM(machine_string)
      WRITE(iunit,*) 'DATE: ',TRIM(date)
      WRITE(iunit,*) 'NVERTEX: ',nvertex
      WRITE(iunit,*) 'NTET: ',ntet
      WRITE(iunit,*) 'NSTATE: ',nstate
      END SUBROUTINE mumaterial_info
      




      SUBROUTINE mumaterial_init(getBfld, offset, comm)
      !-----------------------------------------------------------------------
      ! mumaterial_init: Calculates magnetization of material
      !-----------------------------------------------------------------------
      ! fcn           : getBfld. Function which returns the vacuum magnetic field
      !                 SUBROUTINE FCN(x,y,z,bx,by,bz)
      ! param[in]: offset. Offset of all tiles from the origin
      ! param[in, out]: comm. MPI communicator, handles shared memory
      !-----------------------------------------------------------------------
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: offset(3)
      INTEGER, INTENT(inout), OPTIONAL :: comm
      EXTERNAL:: getBfld
      INTEGER :: ik
      DOUBLE PRECISION :: mu0, norm, x, y, z, Bx, By, Bz, Bx_n, By_n, Bz_n
      
      mu0 = 16 * atan(1.d0) * 1.d-7

      WRITE(*,*) "Setting up tiles"
      
      allocate(tiles(ntet))

      DO ik = 1, ntet
            ! set vertices
            ! each tetrahedron is a tile
            ! tet(ik,1) contains the relevant vertex index for tile vertex 1
            ! vertex(index,:) returns the relevant coordinates
            tiles(ik)%vert(:,1) = [vertex(tet(ik,1),1)+offset(1), vertex(tet(ik,1),2)+offset(2), vertex(tet(ik,1),3)+offset(3)]
            tiles(ik)%vert(:,2) = [vertex(tet(ik,2),1)+offset(1), vertex(tet(ik,2),2)+offset(2), vertex(tet(ik,2),3)+offset(3)]
            tiles(ik)%vert(:,3) = [vertex(tet(ik,3),1)+offset(1), vertex(tet(ik,3),2)+offset(2), vertex(tet(ik,3),3)+offset(3)]
            tiles(ik)%vert(:,4) = [vertex(tet(ik,4),1)+offset(1), vertex(tet(ik,4),2)+offset(2), vertex(tet(ik,4),3)+offset(3)]

            ! tetrahedron tile
            tiles(ik)%tileType = 5
            tiles(ik)%includeInIteration = 1
            
            ! todo check if necessary
            tiles(ik)%stateFunctionIndex = 0
            tiles(ik)%Mrem = 0
            tiles(ik)%u_ea = 0.0
            tiles(ik)%u_oa1 = 0.0
            tiles(ik)%u_oa2 = 0.0

            ! set applied field [A/m] at centre of tile
            x = (tiles(ik)%vert(1,1) + tiles(ik)%vert(1,2) + tiles(ik)%vert(1,3) + tiles(ik)%vert(1,4))/4.0
            y = (tiles(ik)%vert(2,1) + tiles(ik)%vert(2,2) + tiles(ik)%vert(2,3) + tiles(ik)%vert(2,4))/4.0
            z = (tiles(ik)%vert(3,1) + tiles(ik)%vert(3,2) + tiles(ik)%vert(3,3) + tiles(ik)%vert(3,4))/4.0
            call getBfld(x, y, z, Bx, By, Bz)
            tiles(ik)%Happ = [Bx/mu0, By/mu0, Bz/mu0]

            ! hard magnet or soft magnet with state function or with constant permeability
            tiles(ik)%magnetType = state_type(state_dex(ik))

            IF (tiles(ik)%magnetType == 1) THEN
               tiles(ik)%mu_r_ea = constant_mu(state_dex(ik))
               tiles(ik)%mu_r_oa = constant_mu(state_dex(ik))

            ! todo change this to use actual remanent magnetisation
               ! set remanent magnetization to strength of the field [A/m] and easy axis in the direction of the field
               norm = sqrt(Bx*Bx + By*By + Bz*Bz)
               Bx_n = Bx / norm
               By_n = By / norm
               Bz_n = Bz / norm
               tiles(ik)%Mrem = norm / mu0
               tiles(ik)%u_ea = [Bx_n, By_n, Bz_n]
               
               ! get orthogonal axes
               IF (By/=0 .or. Bz/=0) THEN          ! cross product of u_ea with [1, 0, 0] and cross product of u_ea with cross product
                     tiles(ik)%u_oa1 = [0.d0, Bz_n, -By_n]
                     tiles(ik)%u_oa2 = [-By_n*By_n - Bz_n*Bz_n, Bx_n*By_n, Bx_n*Bz_n]
               ELSE                                ! cross product of u_ea with [0, 1, 0] and cross product of u_ea with cross product
                     tiles(ik)%u_oa1 = [-Bz_n, 0.d0, Bx_n]
                     tiles(ik)%u_oa2 = [Bx_n*By_n, -Bx_n*Bx_n - Bz_n*Bz_n, By_n*Bz_n]
               END IF

            ELSEIF (tiles(ik)%magnetType == 2) THEN
               PRINT *, 'STATE_TYPE == 2 (soft magnet) not yet supported.'
               ! tiles(ik)%stateFunctionIndex = some_array(state_dex(ik))
               ! mu_r should not have to be set here
            ELSEIF (tiles(ik)%magnetType == 3) THEN
               tiles(ik)%mu_r_ea = constant_mu(state_dex(ik))
               tiles(ik)%mu_r_oa = constant_mu(state_dex(ik))
            END IF
      

            


            ! preset magnetization
            !tiles(ik)%M = tiles(ik)%Mrem * tiles(ik)%u_ea

            ! call setupEvaluationPoints(tiles(ik)) ! only for cylindrical tiles
      END DO

      WRITE(*,*) "Running iterations"      
      
      allocate(stateFunction(1)) ! todo temporary until statefunction is properly implemented
      call loadStateFunctionFortran(stateFunction)


      CALL iterateMagnetization(tiles, ntet, stateFunction, size(stateFunction), temp, maxErr, maxIter, 0.d0) ! todo replace size?

      RETURN
      END SUBROUTINE mumaterial_init

      subroutine loadStateFunctionFortran(stateFunction)
            integer :: i
            type(MagStateFunction), dimension(1), intent(out) :: stateFunction
            OPEN(11, file='stateFunction.dat')
            read(11,*) stateFunction(1)%nT,stateFunction(1)%nH
            
            allocate( stateFunction(1)%M(stateFunction(1)%nT,stateFunction(1)%nH) )
            allocate( stateFunction(1)%T(stateFunction(1)%nT) )
            allocate( stateFunction(1)%H(stateFunction(1)%nH) )
            
            stateFunction(1)%T(1) = 300
            DO i=1,stateFunction(1)%nH
            read(11,*) stateFunction(1)%H(i),stateFunction(1)%M(1,i)        
            END DO
            
            close (11)
      end subroutine loadStateFunctionFortran
      


      SUBROUTINE mumaterial_getb_scalar(x, y, z, Bx, By, Bz)
      !-----------------------------------------------------------------------
      ! mumaterial_getb: Calculates magnetic field at a point in space
      !-----------------------------------------------------------------------
      ! param[in]: x. x-coordinate of point where to get the B-field
      ! param[in]: y. y-coordinate of point where to get the B-field
      ! param[in]: z. z-coordinate of point where to get the B-field
      ! param[out]: Bx. x-component of B-field at this point [T]
      ! param[out]: By. y-component of B-field at this point [T]
      ! param[out]: Bz. z-component of B-field at this point [T]
      !-----------------------------------------------------------------------
      DOUBLE PRECISION, INTENT(in) :: x, y, z
      DOUBLE PRECISION, INTENT(out) :: Bx, By, Bz
      DOUBLE PRECISION :: H(1,3), point(1,3), mu0

      mu0 = 16 * atan(1.d0) * 1.d-7

      point(1,:) = [x, y, z]

      call getFieldFromTiles(tiles, H, point, ntet, 1) ! todo likely possible to directly call function related to tetrahedrons

      Bx = H(1,1) * mu0
      By = H(1,2) * mu0
      Bz = H(1,3) * mu0

      RETURN
      END SUBROUTINE mumaterial_getb_scalar





      SUBROUTINE mumaterial_getb_vector(x, y, z, Bx, By, Bz)
      !-----------------------------------------------------------------------
      ! mumaterial_getb_vector: Calculates magnetic field at multiple points in space
      !-----------------------------------------------------------------------
      ! param[in]: x. x-coordinates of points at which to determine the magnetic field
      ! param[in]: y. y-coordinates of points at which to determine the magnetic field
      ! param[in]: z. z-coordinates of points at which to determine the magnetic field
      ! param[out]: Bx. x-value of B-field at required points [T]
      ! param[out]: By. y-value of B-field at required points [T]
      ! param[out]: Bz. z-value of B-field at required points [T]
      !-----------------------------------------------------------------------
      DOUBLE PRECISION, INTENT(in) :: x(:), y(:), z(:)
      DOUBLE PRECISION, INTENT(out), ALLOCATABLE :: Bx(:), By(:), Bz(:)
      INTEGER :: n_points
      DOUBLE PRECISION :: mu0
      DOUBLE PRECISION, ALLOCATABLE :: B(:,:), points(:,:)

      n_points = size(x)

      allocate(points(n_points,3))
      allocate(B(n_points,3))
      allocate(Bx(n_points))
      allocate(By(n_points))
      allocate(Bz(n_points))

      mu0 = 16 * atan(1.d0) * 1.d-7

      points(:,1) = x
      points(:,2) = y
      points(:,3) = z

      call getFieldFromTiles(tiles, B, points, ntet, n_points) ! todo likely possible to directly call function related to tetrahedrons

      B = B * mu0
      Bx = B(:,1)
      By = B(:,2)
      Bz = B(:,3)
      
      RETURN
      END SUBROUTINE mumaterial_getb_vector





      SUBROUTINE mumaterial_output(path, x, y, z)
      !-----------------------------------------------------------------------
      ! mumaterial_output: Outputs tiles, H-field and points to text files
      !-----------------------------------------------------------------------
      ! param[in]: path. Path to store files in
      ! param[in]: x. x-cooridinates of points at which to determine the magnetic field
      ! param[in]: y. y-cooridinates of points at which to determine the magnetic field
      ! param[in]: z. z-cooridinates of points at which to determine the magnetic field
      !-----------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(in) :: path
      DOUBLE PRECISION, INTENT(in) :: x(:), y(:), z(:)
      INTEGER :: n_points
      DOUBLE PRECISION, ALLOCATABLE :: H(:,:), points(:,:)
      INTEGER :: i, j

      WRITE(*,*) "Outputting tiles"

      n_points = size(x)
      allocate(points(n_points,3))
      points(:,1) = x
      points(:,2) = y
      points(:,3) = z

      OPEN(12, file=TRIM(path)//'/tiles.dat')
        WRITE(12, "(I15)") ntet
        DO i = 1, ntet
            WRITE(12,"(F15.7,A,F15.7,A,F15.7,A)") tiles(i)%u_ea(1),',',tiles(i)%u_ea(2),',',tiles(i)%u_ea(3),',   :u_ea_x,u_ea_y,u_ea_z'
            WRITE(12,"(F15.7,A,F15.7,A,F15.7,A)") tiles(i)%Mrem,',',tiles(i)%mu_r_ea,',',tiles(i)%mu_r_oa,',   :Mrem,mu_r_ea,mu_r_oa'
            WRITE(12,"(I1,A,I1,A,I1,A,I1,A)") tiles(i)%tileType,',',tiles(i)%magnetType,',',tiles(i)%stateFunctionIndex,',',tiles(i)%includeInIteration,',   :tileType,magnetType,stFcnIndex,inclIter'
            DO j = 1, 4
                WRITE(12,"(F15.7,A,F15.7,A,F15.7,A)") tiles(i)%vert(1,j), ',', tiles(i)%vert(2,j), ',', tiles(i)%vert(3,j),',   :vert'
            END DO
        END DO
      CLOSE(12)

      WRITE(*,*) "Outputting points"

      OPEN(13, file=TRIM(path)//'/points.dat')
        DO i = 1, n_points
            WRITE(13, "(F15.7,A,F15.7,A,F15.7)") points(i, 1), ',', points(i, 2), ',', points(i, 3)
        END DO
      CLOSE(13)

      WRITE(*,*) "Getting H-field"

      allocate(H(n_points,3))
      call getFieldFromTiles(tiles, H, points, ntet, n_points)

      WRITE(*,*) "Outputting H-field"

      OPEN(14, file=TRIM(path)//'/H.dat')
        DO i = 1, n_points
            WRITE(14, "(E15.7,A,E15.7,A,E15.7)") H(i, 1), ',', H(i, 2), ',', H(i, 3)
        END DO
      CLOSE(14)

      RETURN
      END SUBROUTINE





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Memory Allocation Subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE mpialloc_1d_int(array,n1,subid,mymaster,share_comm,win)
      ! Libraries
      USE MPI
      USE ISO_C_BINDING
      IMPLICIT NONE
      ! Arguments
      INTEGER, POINTER, INTENT(inout) :: array(:)
      INTEGER, INTENT(in) :: n1
      INTEGER, INTENT(in) :: subid
      INTEGER, INTENT(in) :: mymaster
      INTEGER, INTENT(in) :: share_comm
      INTEGER, INTENT(inout) :: win
      ! Variables
      INTEGER :: disp_unit, ier
      INTEGER :: array_shape(1)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: window_size
      TYPE(C_PTR) :: baseptr
      ! Initialization
      ier = 0
      array_shape(1) = n1
      disp_unit = 1
      window_size = 0_MPI_ADDRESS_KIND
      IF (subid == mymaster) window_size = INT(n1,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
      IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
      CALL C_F_POINTER(baseptr, array, array_shape)
      RETURN
      END SUBROUTINE mpialloc_1d_int

      SUBROUTINE mpialloc_1d_dbl(array,n1,subid,mymaster,share_comm,win)
      !-----------------------------------------------------------------------
      ! mpialloc_1d_int: Allocated a 1D double array to shared memory
      ! Taken from LIBSTELL/Sources/Modules/mpi_sharemem.f90
      ! Included here to reduce dependencies
      !-----------------------------------------------------------------------
      ! Libraries
#if defined(MPI_OPT)
      USE mpi
#endif
      USE ISO_C_BINDING
      IMPLICIT NONE
      ! Arguments
      DOUBLE PRECISION, POINTER, INTENT(inout) :: array(:)
      INTEGER, INTENT(in) :: n1
      INTEGER, INTENT(in) :: subid
      INTEGER, INTENT(in) :: mymaster
      INTEGER, INTENT(inout) :: share_comm
      INTEGER, INTENT(inout) :: win
      ! Variables
      INTEGER :: disp_unit, ier
      INTEGER :: array_shape(1)
#if defined(MPI_OPT)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: window_size
#endif
      TYPE(C_PTR) :: baseptr
      ! Initialization
      ier = 0
      array_shape(1) = n1
      disp_unit = 1
#if defined(MPI_OPT)
      window_size = 0_MPI_ADDRESS_KIND
      IF (subid == mymaster) window_size = INT(n1,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
      IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
      CALL C_F_POINTER(baseptr, array, array_shape)
#endif
      RETURN
      END SUBROUTINE mpialloc_1d_dbl

      SUBROUTINE mpialloc_2d_int(array,n1,n2,subid,mymaster,share_comm,win)
      ! Libraries
      USE MPI
      USE ISO_C_BINDING
      IMPLICIT NONE
      ! Arguments
      INTEGER, POINTER, INTENT(inout) :: array(:,:)
      INTEGER, INTENT(in) :: n1
      INTEGER, INTENT(in) :: n2
      INTEGER, INTENT(in) :: subid
      INTEGER, INTENT(in) :: mymaster
      INTEGER, INTENT(in) :: share_comm
      INTEGER, INTENT(inout) :: win
      ! Variables
      INTEGER :: disp_unit, ier
      INTEGER :: array_shape(2)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: window_size
      TYPE(C_PTR) :: baseptr
      ! Initialization
      ier = 0
      array_shape(1) = n1
      array_shape(2) = n2
      disp_unit = 1
      window_size = 0_MPI_ADDRESS_KIND
      IF (subid == mymaster) window_size = INT(n1*n2,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
      IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
      CALL C_F_POINTER(baseptr, array, array_shape)
      RETURN
      END SUBROUTINE mpialloc_2d_int

      SUBROUTINE mpialloc_2d_dbl(array,n1,n2,subid,mymaster,share_comm,win)
      !-----------------------------------------------------------------------
      ! mpialloc_1d_int: Allocated a 2D double array to shared memory
      ! Taken from LIBSTELL/Sources/Modules/mpi_sharemem.f90
      ! Included here to reduce dependencies
      !-----------------------------------------------------------------------
      ! Libraries
#if defined(MPI_OPT)
      USE mpi
#endif
      USE ISO_C_BINDING
      IMPLICIT NONE
      ! Arguments
      DOUBLE PRECISION, POINTER, INTENT(inout) :: array(:,:)
      INTEGER, INTENT(in) :: n1
      INTEGER, INTENT(in) :: n2
      INTEGER, INTENT(in) :: subid
      INTEGER, INTENT(in) :: mymaster
      INTEGER, INTENT(inout) :: share_comm
      INTEGER, INTENT(inout) :: win
      ! Variables
      INTEGER :: disp_unit, ier
      INTEGER :: array_shape(2)
#if defined(MPI_OPT)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: window_size
#endif
      TYPE(C_PTR) :: baseptr
      ! Initialization
      ier = 0
      array_shape(1) = n1
      array_shape(2) = n2
      disp_unit = 1
#if defined(MPI_OPT)
      window_size = 0_MPI_ADDRESS_KIND
      IF (subid == mymaster) window_size = INT(n1*n2,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
      IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
      CALL C_F_POINTER(baseptr, array, array_shape)
#endif
      RETURN
      END SUBROUTINE mpialloc_2d_dbl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Memory Freeing Subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
         SUBROUTINE free_mpi_array1d_dbl(win_local,array_local,isshared)
         IMPLICIT NONE
         LOGICAL, INTENT(in) :: isshared
         INTEGER, INTENT(inout) :: win_local
         DOUBLE PRECISION, POINTER, INTENT(inout) :: array_local(:)
         INTEGER :: istat
         istat=0
#if defined(MPI_OPT)
         IF (isshared) THEN
            CALL MPI_WIN_FENCE(0, win_local,istat)
            CALL MPI_WIN_FREE(win_local,istat)
            IF (ASSOCIATED(array_local)) NULLIFY(array_local)
         ELSE
#endif
            IF (ASSOCIATED(array_local)) DEALLOCATE(array_local)
#if defined(MPI_OPT)
         ENDIF
#endif
         RETURN
         END SUBROUTINE free_mpi_array1d_dbl
   
         SUBROUTINE free_mpi_array2d_dbl(win_local,array_local,isshared)
         IMPLICIT NONE
         LOGICAL, INTENT(in) :: isshared
         INTEGER, INTENT(inout) :: win_local
         DOUBLE PRECISION, POINTER, INTENT(inout) :: array_local(:,:)
         INTEGER :: istat
         istat=0
#if defined(MPI_OPT)
         IF (isshared) THEN
            CALL MPI_WIN_FENCE(0, win_local,istat)
            CALL MPI_WIN_FREE(win_local,istat)
            IF (ASSOCIATED(array_local)) NULLIFY(array_local)
         ELSE
#endif
            IF (ASSOCIATED(array_local)) DEALLOCATE(array_local)
#if defined(MPI_OPT)
         ENDIF
#endif
         RETURN
         END SUBROUTINE free_mpi_array2d_dbl

!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE mumaterial_mod
