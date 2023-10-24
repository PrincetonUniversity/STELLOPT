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
      ! USE TileNComponents
      ! USE IterateMagnetSolution
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
      TYPE stateFunctionType
            DOUBLE PRECISION, PRIVATE, ALLOCATABLE :: H(:), M(:)
      END TYPE stateFunctionType

      TYPE nearestNeighboursType
            INTEGER, PRIVATE, ALLOCATABLE :: idx(:)
      END TYPE nearestNeighboursType

      TYPE NStoreType
            DOUBLE PRECISION, PRIVATE, ALLOCATABLE :: N(:,:,:)
      END TYPE
      
      INTEGER, PRIVATE                    :: nvertex, ntet, nstate, nMH
      TYPE(stateFunctionType), PRIVATE, ALLOCATABLE     :: stateFunction(:)
      TYPE(nearestNeighboursType), PRIVATE, ALLOCATABLE :: nearestNeighbours(:)
      LOGICAL, PRIVATE                    :: neighbours
      DOUBLE PRECISION, PRIVATE           :: maxErr, lambdaStart, lambdaFactor, maxDist
      INTEGER, PRIVATE                    :: maxIter


      DOUBLE PRECISION, POINTER, PRIVATE :: vertex(:,:), tet_cen(:,:)
      INTEGER, POINTER, PRIVATE :: tet(:,:)
      INTEGER, POINTER, PRIVATE :: state_dex(:)
      INTEGER, POINTER, PRIVATE :: state_type(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: constant_mu(:), constant_mu_o(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: M(:,:), Happ(:,:), Mrem(:,:)
      INTEGER, PRIVATE            :: win_vertex, win_tet,  &
                                     win_state_dex, win_state_type, &
                                     win_constant_mu, win_m

      CHARACTER(LEN=256), PRIVATE :: machine_string
      CHARACTER(LEN=256), PRIVATE :: date


      INTEGER, PRIVATE                    :: mystart, myend, mydelta
      INTEGER, PRIVATE                    :: shar_rank, shar_size, shar_comm




      
!-----------------------------------------------------------------------
!     Subroutines
!         mumaterial_setd: Sets default values
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
      




      SUBROUTINE mumaterial_setd(mE, mI, la, laF, mD)
      !-----------------------------------------------------------------------
      ! mumaterial_setd: Sets default values
      !-----------------------------------------------------------------------
      ! param[in]: mE. New maxErr: max error for MagTense convergence
      ! param[in]: mI. New maxIter: max amount of MagTense iterations
      ! param[in]: T. New temp: temperature of magnetic material in MagTense
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: mE, la, laF, mD
      INTEGER, INTENT(in) :: mI

      WRITE(*,'(A18,E10.4,A17,I4,A9,E10.4,A16,E10.4)') "Setting max Error:", mE, " max Iterations: ", mI, " lambda: ", la, " lambda factor: ", laF
      IF (mD .le. 0.0d0) THEN
            WRITE(*,'(A27)') "Nearest neighbours disabled"
      ELSE
            WRITE(*,'(A42,E10.4)') "Nearest neighbours enabled; max distance: ", mD
      END IF

      maxErr = mE
      maxIter = mI
      lambdaStart = la
      lambdaFactor = laF
      maxDist = mD

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
      ! todo: tet_cen, M, Happ, Mrem, constant_mu_o, nMH, stateFunction, vertex, tet, N_store, nearestneighbours, maxDist
         CALL MPI_Bcast(nvertex,1,MPI_INTEGER,0,shar_comm,istat)
         CALL MPI_Bcast(ntet,1,MPI_INTEGER,0,shar_comm,istat)
         CALL MPI_Bcast(nstate,1,MPI_INTEGER,0,shar_comm,istat)
         CALL mpialloc_2d_dbl(vertex,nvertex,3,shar_rank,0,shar_comm,win_vertex)
         CALL mpialloc_2d_int(tet,ntet,4,shar_rank,0,shar_comm,win_tet)
         CALL mpialloc_1d_int(state_dex,ntet,shar_rank,0,shar_comm,win_state_dex)
         CALL mpialloc_1d_int(state_type,nstate,shar_rank,0,shar_comm,win_state_type)
         CALL mpialloc_1d_dbl(constant_mu,nstate,shar_rank,0,shar_comm,win_constant_mu)
      !    CALL mpialloc_2d_dbl(Mag,nstate,3,shar_rank,0,shar_comm,win_vertex)
         shared = .true.
      ELSE
#endif
         ! if no MPI, allocate everything on one node
         ALLOCATE(vertex(3,nvertex),tet(4,ntet),state_dex(ntet), &
                  state_type(nstate),constant_mu(nstate), &
                  tet_cen(3,ntet),M(3,ntet),Happ(3,ntet), &
                  constant_mu_o(nstate),Mrem(3,ntet),stateFunction(nstate), &
                  nearestNeighbours(ntet),STAT=istat)
         shared = .false.
#if defined(MPI_OPT)
      END IF
#endif
      WRITE(*,*) "Reading tile data"
      ! read in the mesh
      IF (istat/=0) RETURN
      IF (shar_rank == 0) THEN
         DO ik = 1, nvertex
            READ(iunit,*) vertex(1,ik),vertex(2,ik),vertex(3,ik)
         END DO

         DO ik = 1, ntet
            READ(iunit,*) tet(1,ik),tet(2,ik),tet(3,ik),tet(4,ik),state_dex(ik)
         END DO

         DO ik = 1, nstate
            READ(iunit,*) state_type(ik)
            IF (state_type(ik) == 1) THEN
               READ(iunit,*) constant_mu(ik), constant_mu_o(ik)
               READ(iunit,*) Mrem(1,ik),Mrem(2,ik),Mrem(3,ik)
            ELSEIF (state_type(ik) == 2) THEN
               READ(iunit,*) nMH
               ALLOCATE(stateFunction(ik)%H(nMH),stateFunction(ik)%M(nMH))
               DO i = 1, nMH
                  READ(iunit,*) stateFunction(ik)%H(i),stateFunction(ik)%M(i)
               END DO
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
      CALL MUMATERIAL_SETD(1.0d-5, 100, 0.7d0, 0.75d0, 0.d0)

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
      INTEGER :: i, j, j_tile, n
      DOUBLE PRECISION :: mu0, norm, x, y, z, Bx, By, Bz, Bx_n, By_n, Bz_n
      DOUBLE PRECISION, ALLOCATABLE :: neighbours_temp(:)
      DOUBLE PRECISION, ALLOCATABLE :: N_store_all(:,:,:,:)
      TYPE(NStoreType), ALLOCATABLE :: N_store(:)
      
      mu0 = 16 * atan(1.d0) * 1.d-7

      WRITE(*,*) "Setting up tiles"
      

      IF (MAXVAL(ABS(offset)) .gt. 0.d0) THEN
            DO i = 1, nvertex
                  vertex(:,i) = vertex(:,i) + offset
            END DO
      END IF

      DO i = 1, ntet
            tet_cen(:,i) = (vertex(:,tet(1,i)) + vertex(:,tet(2,i)) + vertex(:,tet(3,i)) + vertex(:,tet(4,i)))/4.d0

            CALL getBfld(tet_cen(1,i), tet_cen(2,i), tet_cen(3,i), Bx, By, Bz)
            Happ(:,i) = [Bx/mu0, By/mu0, Bz/mu0]
      END DO

      IF (maxDist .gt. 0.d0) THEN ! Double approach since not using neighbours is faster for smaller objects
            WRITE (*,*) "Determining nearest neighbours"
            ALLOCATE(N_store(ntet))
            DO i = 1, ntet
                  ALLOCATE(nearestNeighbours(i)%idx(1))
                  nearestNeighbours(i)%idx(1) = i
            END DO

            DO i = 1, ntet                  
                  IF (MOD(i, 100) .eq. 0) THEN
                        WRITE(*,'(A6,I8,A11,F6.2,A1)') "Tile: ", i, " Complete: ", i*1.d2/ntet, '%'
                  END IF
                  DO j = i+1, ntet
                        IF (i .ne. j .AND. NORM2(tet_cen(:,i)-tet_cen(:,j)) .lt. maxDist) THEN
                              ! ! set j as nearest neighbour for i
                              ! n = SIZE(nearestNeighbours(i)%idx)
                              ! ALLOCATE(neighbours_temp(n))
                              ! neighbours_temp = nearestNeighbours(i)%idx
                              ! DEALLOCATE(nearestNeighbours(i)%idx)
                              ! ALLOCATE(nearestNeighbours(i)%idx(n+1))
                              ! nearestNeighbours(i)%idx(1:n) = neighbours_temp
                              ! nearestNeighbours(i)%idx(n+1) = j
                              ! DEALLOCATE(neighbours_temp)

                              ! ! set i as nearest neighbour for j
                              ! n = SIZE(nearestNeighbours(j)%idx)
                              ! ALLOCATE(neighbours_temp(n))
                              ! neighbours_temp = nearestNeighbours(j)%idx
                              ! DEALLOCATE(nearestNeighbours(j)%idx)
                              ! ALLOCATE(nearestNeighbours(j)%idx(n+1))
                              ! nearestNeighbours(j)%idx(1:n) = neighbours_temp
                              ! nearestNeighbours(j)%idx(n+1) = i
                              ! DEALLOCATE(neighbours_temp)

                              nearestNeighbours(i)%idx = [nearestNeighbours(i)%idx, j]
                              nearestNeighbours(j)%idx = [nearestNeighbours(j)%idx, i]
                        END IF
                  END DO
                  ALLOCATE(N_store(i)%N(3,3,SIZE(nearestNeighbours(i)%idx)))
            END DO

            WRITE(*,*) "Determining N-tensors"
            DO i = 1, ntet
                  IF (MOD(i, 100) .eq. 0) THEN
                        WRITE(*,'(A6,I8,A11,F6.2,A1)') "Tile: ", i, " Complete: ", i*1.d2/ntet, '%'
                  END IF
                  DO j = 1, SIZE(nearestNeighbours(i)%idx)
                        CALL mumaterial_getN(vertex(:,tet(1,nearestNeighbours(i)%idx(j))), vertex(:,tet(2,nearestNeighbours(i)%idx(j))), vertex(:,tet(3,nearestNeighbours(i)%idx(j))), vertex(:,tet(4,nearestNeighbours(i)%idx(j))), tet_cen(:,i), N_store(i)%N(:,:,j))
                  END DO
            END DO
            neighbours = .TRUE.
      ELSE
            WRITE(*,*) "Determining N-tensors"
            ALLOCATE(N_store_all(3,3,ntet,ntet))
            DO i = 1, ntet                  
                  DO j = 1, ntet
                        CALL mumaterial_getN(vertex(:,tet(1,j)), vertex(:,tet(2,j)), vertex(:,tet(3,j)), vertex(:,tet(4,j)), tet_cen(:,i), N_store_all(:,:,j,i))
                  END DO
            END DO
            neighbours = .FALSE.
      END IF



      WRITE(*,*) "Running iterations"      

      IF (neighbours) THEN
            CALL mumaterial_iterate_magnetization(lambdaStart, lambdaFactor, N_store=N_store)
      ELSE
            CALL mumaterial_iterate_magnetization(lambdaStart, lambdaFactor, N_store_all=N_store_all)
      END IF

      RETURN
      END SUBROUTINE mumaterial_init









      
      SUBROUTINE mumaterial_iterate_magnetization(lambdaStart, lambdaFac, N_store_all, N_store)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: lambdaStart, lambdaFac
      INTEGER :: i, j, i_tile, j_tile, count, lambdaCount
      LOGICAL :: run
      DOUBLE PRECISION :: Mnorm(ntet), Mnorm_old(ntet), H(3), N(3,3), chi(ntet)
      DOUBLE PRECISION :: H_old(3), H_new(3), M_new(3,ntet), lambda, lambda_s, error, maxDiff(4), Hnorm, M_tmp_norm
      DOUBLE PRECISION, OPTIONAL :: N_store_all(3,3,ntet,ntet)
      TYPE(NStoreType), OPTIONAL :: N_store(ntet)
      DOUBLE PRECISION :: M_tmp(3), M_tmp_local(3), Mrem_norm, u_ea(3), u_oa_1(3), u_oa_2(3) ! hard magnet

      count = 0
      lambda = lambdaStart
      lambdaCount = 0
      Mnorm_old = 0.d0
      maxDiff = 0.d0
      ! Iterate on magnetization
      DO
            ! Get the field and new magnetization for each tile
            DO i_tile = 1, ntet
                  H = Happ(:,i_tile)

                  ! Get the field from all other tiles
                  IF (neighbours) THEN
                        DO j = 2, SIZE(nearestNeighbours(i_tile)%idx) ! Current tile is not included in H
                                    H = H + MATMUL(N_store(i_tile)%N(:,:,j), M(:,nearestNeighbours(i_tile)%idx(j)))
                        END DO
                        N = N_store(i_tile)%N(:,:,1)  
                  ELSE
                        DO j_tile = 1, ntet
                              IF (i_tile .ne. j_tile) THEN ! Current tile is not included in H
                                    H = H + MATMUL(N_store_all(:,:,j_tile,i_tile), M(:,j_tile))
                              END IF
                        END DO
                        N = N_store_all(:,:,i_tile,i_tile)
                  END IF

                  H_new = H
                  SELECT CASE (state_type(state_dex(i_tile)))
                  CASE (1) ! Hard magnet
                        Mrem_norm = NORM2(Mrem(:,state_dex(i_tile)))
                        u_ea = Mrem(:,state_dex(i_tile))/Mrem_norm ! easy axis assumed parallel to remanent magnetization
                        IF (u_ea(2)/=0 .OR. u_ea(3)/=0) THEN      ! cross product of u_ea with [1, 0, 0] and cross product of u_ea with cross product
                              u_oa_1 = [0.d0, u_ea(3), -u_ea(2)]
                              u_oa_2 = [-u_ea(2)*u_ea(2) - u_ea(3)*u_ea(3), u_ea(1)*u_ea(2), u_ea(1)*u_ea(3)]
                        ELSE                                      ! cross product of u_ea with [0, 1, 0] and cross product of u_ea with cross product
                              u_oa_1 = [-u_ea(3), 0.d0, u_ea(1)]
                              u_oa_2 = [u_ea(1)*u_ea(2), -u_ea(1)*u_ea(1) - u_ea(3)*u_ea(3), u_ea(2)*u_ea(3)]
                        END IF

                        u_oa_1 = u_oa_1/NORM2(u_oa_1)
                        u_oa_2 = u_oa_2/NORM2(u_oa_2)
                        
                        lambda_s = MIN(1/constant_mu(state_dex(i_tile)), 1/constant_mu_o(state_dex(i_tile)), 0.5)
                        DO
                              H_old = H_new

                              ! Determine magnetization taking into account easy axis
                              M_tmp = (Mrem_norm + (constant_mu(state_dex(i_tile)) - 1) * DOT_PRODUCT(H_new, u_ea)) * u_ea &
                                    + (constant_mu_o(state_dex(i_tile)) - 1) * DOT_PRODUCT(H_new, u_oa_1) * u_oa_1 &
                                    + (constant_mu_o(state_dex(i_tile)) - 1) * DOT_PRODUCT(H_new, u_oa_2) * u_oa_2

                              H_new = H + MATMUL(N, M_tmp)
                              H_new = H_old + lambda_s * (H_new - H_old)

                              IF (MAXVAL(ABS((H_new - H_old)/H_old)) .lt. maxErr*lambda_s) THEN
                                    M_new(:,i_tile) = (Mrem_norm + (constant_mu(state_dex(i_tile)) - 1) * DOT_PRODUCT(H_new, u_ea)) * u_ea &
                                                    + (constant_mu_o(state_dex(i_tile)) - 1) * DOT_PRODUCT(H_new, u_oa_1) * u_oa_1 &
                                                    + (constant_mu_o(state_dex(i_tile)) - 1) * DOT_PRODUCT(H_new, u_oa_2) * u_oa_2
                                    EXIT
                              END IF
                        END DO
                  CASE (2) ! Soft magnet using state function
                        DO
                              H_old = H_new
                              Hnorm = NORM2(H_new)
                              IF (Hnorm .ne. 0) THEN
                                    CALL mumaterial_getState(stateFunction(state_dex(i_tile))%H, stateFunction(state_dex(i_tile))%M, Hnorm, M_tmp_norm)
                                    M_tmp = M_tmp_norm * H_new / Hnorm
                              ELSE
                                    M_tmp = 0
                                    M_tmp_norm = 0
                              END IF
                              lambda_s = MIN(Hnorm/M_tmp_norm, 0.5)

                              H_new = H + MATMUL(N, M_tmp)
                              H_new = H_old + lambda_s * (H_new - H_old)

                              IF (MAXVAL(ABS((H_new - H_old)/H_old)) .lt. maxErr*lambda_s) THEN
                                    Hnorm = NORM2(H_new)
                                    IF (Hnorm .ne. 0) THEN
                                          CALL mumaterial_getState(stateFunction(state_dex(i_tile))%H, stateFunction(state_dex(i_tile))%M, Hnorm, M_tmp_norm)
                                          M_new(:,i_tile) = M_tmp_norm * H_new / Hnorm
                                          chi(i_tile) = M_tmp_norm / Hnorm
                                    ELSE
                                          M_new(:,i_tile) = 0
                                          chi(i_tile) = 0
                                    END IF
                                    EXIT
                              END IF
                        END DO
                  CASE (3) ! Soft magnet using constant permeability
                        lambda_s = MIN(1/constant_mu(state_dex(i_tile)), 0.5)
                        DO
                              H_old = H_new
                              H_new = H + (constant_mu(state_dex(i_tile)) - 1) * MATMUL(N, H_new)
                              H_new = H_old + lambda_s * (H_new - H_old)

                              IF (MAXVAL(ABS((H_new - H_old)/H_old)) .lt. maxErr*lambda_s) THEN
                                    M_new(:,i_tile) = (constant_mu(state_dex(i_tile)) - 1) * H_new
                                    EXIT
                              END IF
                        END DO
                  CASE DEFAULT
                        WRITE(*,*) "Unknown magnet type: ", state_type(state_dex(i_tile))
                        STOP
                  END SELECT

                  
            END DO
            M = M + lambda*(M_new - M)
            
            count = count + 1

            error = 0.d0
            DO i = 1, ntet   
                  Mnorm(i) = NORM2(M(:,i))
                  IF (Mnorm_old(i) .ne. 0.d0 .AND. ABS((Mnorm(i) - Mnorm_old(i))/Mnorm_old(i)) .gt. error) THEN
                        error = ABS((Mnorm(i) - Mnorm_old(i))/Mnorm_old(i))
                  END IF
            END DO



            WRITE(*,'(A7,I5,A8,E15.7,A13,E15.7)') 'Count: ', count, ' Error: ', error, ' Max. Error: ', maxErr*lambda

            IF (count .gt. 2 .AND. (error .lt. maxErr*lambda .OR. count .gt. maxIter)) THEN
                  WRITE(*,'(A14,E15.7)') "Average mu_r: ", SUM(chi)/ntet + 1
                  EXIT
            ELSE
                  maxDiff = CSHIFT(maxDiff, -1)
                  maxDiff(1) = error

                  ! Update lambda if there is any increase in the error (maxDiff(1) is the most recent)
                  IF (lambdaCount .gt. 4 .AND. ((maxDiff(2) - maxDiff(1)) .lt. 0.d0 .OR. (maxDiff(3) - maxDiff(2)) .lt. 0.d0 .OR. (maxDiff(4) - maxDiff(3)) .lt. 0.d0)) THEN
                        lambda = lambda * lambdaFac
                        lambdaCount = 1
                  END IF

                  lambdaCount = lambdaCount + 1
                  
                  Mnorm_old = Mnorm
            END IF
      END DO

      WRITE(*,*) "Finished determining magnetization"

      RETURN
      END SUBROUTINE mumaterial_iterate_magnetization





      SUBROUTINE mumaterial_getN(v1, v2, v3, v4, pos, N)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in), DIMENSION(3) :: v1, v2, v3, v4, pos
      DOUBLE PRECISION, INTENT(out) :: N(3,3)

      DOUBLE PRECISION :: N_loc(3,3), v(3,4), v_temp(3), angles(3), P(3,3), Pinv(3,3), D(3), r(3)
      INTEGER :: i, j

      N = 0.d0

      DO i = 1, 4
            ! Shift vertices
            v(:,i) = v1
            v(:,MOD(i,4)+1) = v2
            v(:,MOD(i+1,4)+1) = v3
            v(:,MOD(i+2,4)+1)= v4

            ! todo: ensure vertices are not collinear and v4 is not in plane of v1-3?

            ! Ensure largest angle is for v2
            angles(1) = ACOS(DOT_PRODUCT(v(:,1)-v(:,2),v(:,1)-v(:,3)) / (NORM2(v(:,1)-v(:,2)) * NORM2(v(:,1)-v(:,3))))
            angles(2) = ACOS(DOT_PRODUCT(v(:,2)-v(:,1),v(:,2)-v(:,3)) / (NORM2(v(:,2)-v(:,1)) * NORM2(v(:,2)-v(:,3))))
            angles(3) = ACOS(DOT_PRODUCT(v(:,3)-v(:,2),v(:,3)-v(:,1)) / (NORM2(v(:,3)-v(:,2)) * NORM2(v(:,3)-v(:,1))))

            IF (angles(1) > angles(2) .and. angles(1) > angles(3)) THEN ! v1 and v2 should be interchanged
                  v_temp = v(:,2)
                  v(:,2) = v(:,1)
                  v(:,1) = v_temp
            ELSE IF (angles(3) > angles(1) .and. angles(3) > angles(2)) THEN ! v2 and v3 should be interchanged
                  v_temp = v(:,2)
                  v(:,2) = v(:,3)
                  v(:,3) = v_temp
            END IF

            ! Ensure normal vector is pointing in the right direction
            IF (DOT_PRODUCT(mumaterial_cross(v(:,1) - v(:,3), v(:,2) - v(:,3)), v(:,4) - v(:,2)) .gt. 0) THEN ! normal vector of triangle is pointing towards v4, so v1 and v3 need to be interchanged
                  v_temp = v(:,1)
                  v(:,1) = v(:,3)
                  v(:,3) = v_temp
            END IF

            ! Rotation matrix
            P(:,1) = v(:,1) - v(:,3)
            P(:,1) = P(:,1) / NORM2(P(:,1))

            P(:,3) = mumaterial_cross(P(:,1), v(:,2)-v(:,3))
            P(:,3) = P(:,3) / NORM2(P(:,3))

            P(:,2) = mumaterial_cross(P(:,3), P(:,1))
            P(:,2) = P(:,2) / NORM2(P(:,2))

            ! Inverse rotation matrix, transpose since P is orthogonal
            Pinv = TRANSPOSE(P)

            ! Position of triangle base
            D = DOT_PRODUCT(v(:,3)-v(:,2),v(:,3)-v(:,1)) / (NORM2(v(:,3)-v(:,2)) * NORM2(v(:,3)-v(:,1))) * NORM2(v(:,2) - v(:,3)) * P(:,1) + v(:,3)

            ! Transform evaluation position and vertices to local coordinate frame
            r = MATMUL(Pinv, (pos - D))

            DO j = 1, 3
                  v(:,j) = MATMUL(Pinv, (v(:,j) - D))

                  IF (ABS(r(j)) .lt. 1.0d-20) THEN ! make sure position is not too close to x, y or z = 0
                        r(j) = SIGN(1.0d-20, r(j))
                  END IF
            END DO

            N_loc = 0.d0

            N_loc(1,3) = mumaterial_getNxz(r, v(1,1), v(2,2)) - mumaterial_getNxz(r, v(1,3), v(2,2)) ! todo: exchange indices? also for P
            N_loc(2,3) = mumaterial_getNyz(r, v(1,1), v(2,2)) - mumaterial_getNyz(r, v(1,3), v(2,2))
            N_loc(3,3) = mumaterial_getNzz(r, v(1,1), v(2,2)) - mumaterial_getNzz(r, v(1,3), v(2,2))

            N = N + MATMUL(MATMUL(P, N_loc), Pinv)
      END DO

      RETURN
      END SUBROUTINE mumaterial_getN

      FUNCTION mumaterial_getNxz(r, l, h)
      IMPLICIT NONE
      DOUBLE PRECISION :: mumaterial_getNxz
      DOUBLE PRECISION, INTENT(IN) :: r(3), l, h

            mumaterial_getNxz = -1.d0/(16.d0*ATAN(1.d0)) * (F(r,h,l,h) - F(r,0.d0,l,h) - (G(r,h) - G(r,0.d0)))

            RETURN

      CONTAINS

            FUNCTION F(r, yp, l, h)
            IMPLICIT NONE
            DOUBLE PRECISION :: F
            DOUBLE PRECISION, INTENT(IN) :: r(3), yp, l, h

                  F = h / sqrt(h*h + l*l) * ATANH((l*l - l*r(1) + h*r(2) - h*yp*(1 + l*l/h/h)) / &
                        sqrt((h*h + l*l) * (r(1)*r(1) - 2*r(1)*l + l*l + r(2)*r(2) - 2*(l*l - l*r(1) + h*r(2))*yp/h + &
                        yp*yp*(1 + l*l/h/h) + r(3)*r(3))))

            RETURN
            END FUNCTION F

            FUNCTION G(r, yp)
            IMPLICIT NONE
            DOUBLE PRECISION :: G
            DOUBLE PRECISION, INTENT(IN) :: r(3), yp

                  G = ATANH((r(2) - yp) / sqrt(r(1)*r(1) + r(2)*r(2) - 2*r(2)*yp + yp*yp + r(3)*r(3)))

                  ! todo: check G is finite?

            RETURN
            END FUNCTION G
      END FUNCTION mumaterial_getNxz

      FUNCTION mumaterial_getNyz(r, l, h)
      IMPLICIT NONE
      DOUBLE PRECISION :: mumaterial_getNyz
      DOUBLE PRECISION, INTENT(IN) :: r(3), l, h

            mumaterial_getNyz = -1.d0/(16.d0*ATAN(1.d0)) * (K(r,l,l,h) - K(r,0.d0,l,h) - (Lfunc(r,l) - Lfunc(r,0.d0)))

            RETURN

      CONTAINS

            FUNCTION K(r, xp, l, h)
            IMPLICIT NONE
            DOUBLE PRECISION :: K
            DOUBLE PRECISION, INTENT(IN) :: r(3), xp, l, h

                  K = l / sqrt(h*h + l*l) * ATANH((h*h + l*r(1) - h*r(2) - l*xp*(1 + h*h/l/l)) / &
                        sqrt((h*h + l*l) * (r(2)*r(2) - 2*r(2)*h + h*h + r(1)*r(1) - 2*(h*h + l*r(1) - h*r(2))*xp/l + &
                        xp*xp*(1 + h*h/l/l) + r(3)*r(3))))

            RETURN
            END FUNCTION K

            FUNCTION Lfunc(r, xp)
            IMPLICIT NONE
            DOUBLE PRECISION :: Lfunc
            DOUBLE PRECISION, INTENT(IN) :: r(3), xp

                  Lfunc = ATANH((r(1) - xp) / sqrt(r(1)*r(1) - 2*r(1)*xp + xp*xp + r(2)*r(2) + r(3)*r(3)))

                  ! todo: check Lfunc is finite?

            RETURN
            END FUNCTION Lfunc

      END FUNCTION mumaterial_getNyz

      FUNCTION mumaterial_getNzz(r, l, h)
      IMPLICIT NONE
      DOUBLE PRECISION :: mumaterial_getNzz
      DOUBLE PRECiSION, INTENT(IN) :: r(3), l, h

            mumaterial_getNzz = -1.d0/(16.d0*ATAN(1.d0)) * (P(r,l,l,h) - P(r,0.d0,l,h) - (Q(r,l) - Q(r,0.d0)))

            RETURN
      
      CONTAINS

            FUNCTION P(r, xp, l, h)
            IMPLICIT NONE
            DOUBLE PRECISION :: P
            DOUBLE PRECISION, INTENT(IN) :: r(3), xp, l, h

                  P = ATAN((r(1)*(h - r(2)) - xp*(h*(1 - r(1)/l) - r(2)) - h*(r(1)*r(1) + r(3)*r(3))/l) / &
                        (r(3)*sqrt(r(2)*r(2) - 2*r(2)*h + h*h + r(1)*r(1) + xp*xp*(1 + h*h/l/l) - &
                        2*xp*(h*h + l*r(1) - h*r(2))/l + r(3)*r(3))))

            RETURN
            END FUNCTION P

            FUNCTION Q(r, xp)
            IMPLICIT NONE
            DOUBLE PRECISION :: Q
            DOUBLE PRECISION, INTENT(IN) :: r(3), xp

                  Q = -ATAN((r(1) - xp)*r(2) / (r(3)*sqrt((r(1)*r(1) - 2*r(1)*xp + xp*xp + r(2)*r(2) + r(3)*r(3)))))
            
            RETURN
            END FUNCTION Q

      END FUNCTION mumaterial_getNzz

      FUNCTION mumaterial_cross(a, b)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN), DIMENSION(3) :: a, b
            DOUBLE PRECISION, DIMENSION(3) :: mumaterial_cross

            mumaterial_cross(1) = a(2)*b(3) - a(3)*b(2)
            mumaterial_cross(2) = a(3)*b(1) - a(1)*b(3)
            mumaterial_cross(3) = a(1)*b(2) - a(2)*b(1)

            RETURN
      END FUNCTION mumaterial_cross





      SUBROUTINE mumaterial_getState(fx, fy, x, y)
      !-----------------------------------------------------------------------
      ! mumaterial_getState: Interpolates a function f at x to get a value y using B-splines based on De Boor's algorithm
      !-----------------------------------------------------------------------
      ! param[in]: fx. x-coordinates of function to be interpolated
      ! param[in]: fy. y-coordinates of function to be interpolated
      ! param[in]: x. x-coordinate of evaluation point
      ! param[in]: y. y-coordinate of evaluation point
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: fx(:), fy(:), x
      DOUBLE PRECISION, INTENT(OUT) :: y
      INTEGER :: n, i, k, p, r
      DOUBLE PRECISION :: alpha
      DOUBLE PRECISION, ALLOCATABLE :: t(:), d(:)

      p = 3 ! Degree of the polynomial used
      n = SIZE(fx)

      ALLOCATE(t(2+2*(p-1)))
      ALLOCATE(d(p+1))

      ! Determine left index k
      ! Assume fx in non-decreasing
      IF (x .lt. fx(1)) THEN
            k = 1
      ELSEIF (x .gt. fx(n)) THEN
            k = n
      ELSE
            DO i = 2, n 
                  IF (x .lt. fx(i)) THEN
                        k = i - 1
                        EXIT
                  END IF
            END DO
      END IF

      ! Determine array with x-values, add padding if necessary
      DO i = 1, 2+2*(p-1)
            r = k + i - p
            IF (r .lt. 1) THEN
                  t(i) = fx(1)
            ELSEIF (r .gt. n) THEN
                  t(i) = fx(n)
            ELSE
                  t(i) = fx(r)
            END IF
      END DO

      ! Determine array with coefficients, add padding if necessary
      DO i = 1, p+1
            r = k + i - p
            IF (r .lt. 1) THEN
                  d(i) = fy(1)
            ELSEIF (r .gt. n) THEN
                  d(i) = fy(n)
            ELSE
                  d(i) = fy(r)
            END IF
      END DO

      ! Determine spline coefficients
      DO r = 1, p
            DO i = p, r, -1
                  alpha = (x - t(i)) / (t(i+1+p-r) - t(i))
                  d(i) = (1 - alpha) * d(i) + alpha * d(i+1)
            END DO
      END DO

      ! Set output variable
      y = d(3)

      RETURN
      END SUBROUTINE mumaterial_getState


      


      SUBROUTINE mumaterial_getb_scalar(x, y, z, Bx, By, Bz, getBfld)
      !-----------------------------------------------------------------------
      ! mumaterial_getb: Calculates total magnetic field at a point in space
      !-----------------------------------------------------------------------
      ! param[in]: x. x-coordinate of point where to get the B-field
      ! param[in]: y. y-coordinate of point where to get the B-field
      ! param[in]: z. z-coordinate of point where to get the B-field
      ! param[out]: Bx. x-component of B-field at this point [T]
      ! param[out]: By. y-component of B-field at this point [T]
      ! param[out]: Bz. z-component of B-field at this point [T]
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      EXTERNAL:: getBfld
      DOUBLE PRECISION, INTENT(in) :: x, y, z
      DOUBLE PRECISION, INTENT(out) :: Bx, By, Bz
      DOUBLE PRECISION :: H(3), mu0, N(3,3)
      INTEGER :: i

      mu0 = 16 * atan(1.d0) * 1.d-7
      H = 0.d0

      DO i = 1, ntet
            CALL mumaterial_getN(vertex(:,tet(1,i)), vertex(:,tet(2,i)), vertex(:,tet(3,i)), vertex(:,tet(4,i)), [x, y, z], N)
            H = H + MATMUL(N, M(:,i))
      END DO

      CALL getBfld(x, y, z, Bx, By, Bz)

      Bx = Bx + H(1) * mu0
      By = By + H(2) * mu0
      Bz = Bz + H(3) * mu0

      RETURN
      END SUBROUTINE mumaterial_getb_scalar





      SUBROUTINE mumaterial_getb_vector(x, y, z, Bx, By, Bz, getBfld)
      !-----------------------------------------------------------------------
      ! mumaterial_getb_vector: Calculates total magnetic field at multiple points in space
      !-----------------------------------------------------------------------
      ! param[in]: x. x-coordinates of points at which to determine the magnetic field
      ! param[in]: y. y-coordinates of points at which to determine the magnetic field
      ! param[in]: z. z-coordinates of points at which to determine the magnetic field
      ! param[out]: Bx. x-value of B-field at required points [T]
      ! param[out]: By. y-value of B-field at required points [T]
      ! param[out]: Bz. z-value of B-field at required points [T]
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      EXTERNAL:: getBfld
      DOUBLE PRECISION, INTENT(in) :: x(:), y(:), z(:)
      DOUBLE PRECISION, INTENT(out), ALLOCATABLE :: Bx(:), By(:), Bz(:)
      INTEGER :: n_points, i

      n_points = size(x)

      allocate(Bx(n_points))
      allocate(By(n_points))
      allocate(Bz(n_points))

      DO i = 1, n_points
            CALL mumaterial_getb_scalar(x(i), y(i), z(i), Bx(i), By(i), Bz(i), getBfld)
      END DO
      
      RETURN
      END SUBROUTINE mumaterial_getb_vector



      SUBROUTINE mumaterial_output(path, x, y, z, getBfld)
      !-----------------------------------------------------------------------
      ! mumaterial_output: Outputs tiles, H-field and points to text files
      !-----------------------------------------------------------------------
      ! param[in]: path. Path to store files in
      ! param[in]: x. x-cooridinates of points at which to determine the magnetic field
      ! param[in]: y. y-cooridinates of points at which to determine the magnetic field
      ! param[in]: z. z-cooridinates of points at which to determine the magnetic field
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      EXTERNAL:: getBfld
      CHARACTER(LEN=*), INTENT(in) :: path
      DOUBLE PRECISION, INTENT(in) :: x(:), y(:), z(:)
      INTEGER :: n_points
      DOUBLE PRECISION, ALLOCATABLE :: Bx(:), By(:), Bz(:)
      INTEGER :: i, j

      n_points = size(x)

      WRITE(*,*) "Outputting points"

      OPEN(13, file=TRIM(path)//'/points.dat')
        DO i = 1, n_points
            WRITE(13, "(F15.7,A,F15.7,A,F15.7)") x(i), ',', y(i), ',', z(i)
        END DO
      CLOSE(13)

      WRITE(*,*) "Getting B-field"

      CALL mumaterial_getb(x, y, z, Bx, By, Bz, getBfld)

      WRITE(*,*) "Outputting B-field"

      OPEN(14, file=TRIM(path)//'/B.dat')
        DO i = 1, n_points
            WRITE(14, "(E15.7,A,E15.7,A,E15.7)") Bx(i), ',', By(i), ',', Bz(i)
        END DO
      CLOSE(14)

      RETURN
      END SUBROUTINE





      ! SUBROUTINE mumaterial_field_variation(getBfld, use_demag)
      ! EXTERNAL :: getBfld
      ! LOGICAL, INTENT(in) :: use_demag
      ! INTEGER :: i
      ! DOUBLE PRECISION :: Bnorm, diffNorm1, diffNorm2, diffNorm3, diffNorm4, Bnorm_c, Bxc, Byc, Bzc, theta1, theta2, theta3, theta4, x, y, z
      ! DOUBLE PRECISION :: Htmp(3), mu_0

      ! mu0 = 16 * atan(1.d0) * 1.d-7

      ! DO i = 1, ntet
      !       x = (tiles(ik)%vert(1,1) + tiles(ik)%vert(1,2) + tiles(ik)%vert(1,3) + tiles(ik)%vert(1,4))/4.0
      !       y = (tiles(ik)%vert(2,1) + tiles(ik)%vert(2,2) + tiles(ik)%vert(2,3) + tiles(ik)%vert(2,4))/4.0
      !       z = (tiles(ik)%vert(3,1) + tiles(ik)%vert(3,2) + tiles(ik)%vert(3,3) + tiles(ik)%vert(3,4))/4.0
      !       CALL getBfld(x, y, z, Bxc, Byc, Bzc)
      !       Bnormc = sqrt(Bxc*Bxc + Byc*Byc + Bzc*Bzc)
            
      !       x = tiles(i)%vert(1,1)
      !       y = tiles(i)%vert(2,1)
      !       z = tiles(i)%vert(3,1)
      !       CALL getBfld(x, y, z, Bx, By, Bz)
      !       IF (use_demag) THEN
      !             CALL getFieldFromTetrahedronTile(tiles(i), Htmp, [x, y, z], 1)
      !             Bx = Bx + Htmp(1)*mu_0*tiles(i)%mu_r_ea
      !             By = By + Htmp(2)*mu_0*tiles(i)%mu_r_ea
      !             Bz = Bz + Htmp(3)*mu_0*tiles(i)%mu_r_ea
      !       END IF
      !       Bnorm = sqrt(Bx*Bx + By*By + Bz*Bz)
      !       diffNorm1 = abs(Bnorm-Bnorm_c)
      !       theta1 = 
      ! END DO

      ! RETURN
      ! END SUBROUTINE





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
