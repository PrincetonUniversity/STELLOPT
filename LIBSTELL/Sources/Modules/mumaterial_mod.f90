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
      INTEGER, PRIVATE                    :: nvertex, ntet, nstate
      ! TYPE(MagTile), PRIVATE, ALLOCATABLE              :: tiles(:)
      ! TYPE(MagStateFunction), PRIVATE, ALLOCATABLE     :: stateFunction(:)
      DOUBLE PRECISION, PRIVATE           :: maxErr, lambdaStart, lambdaFactor
      INTEGER, PRIVATE                    :: maxIter


      DOUBLE PRECISION, POINTER, PRIVATE :: vertex(:,:), tet_cen(:,:)
      INTEGER, POINTER, PRIVATE :: tet(:,:)
      INTEGER, POINTER, PRIVATE :: state_dex(:)
      INTEGER, POINTER, PRIVATE :: state_type(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: constant_mu(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: M(:,:), Happ(:,:)
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
      




      SUBROUTINE mumaterial_setd(mE, mI, la, laF)
      !-----------------------------------------------------------------------
      ! mumaterial_setd: Sets default values
      !-----------------------------------------------------------------------
      ! param[in]: mE. New maxErr: max error for MagTense convergence
      ! param[in]: mI. New maxIter: max amount of MagTense iterations
      ! param[in]: T. New temp: temperature of magnetic material in MagTense
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: mE, la, laF
      INTEGER, INTENT(in) :: mI

      WRITE(*,*) "Setting max Error, max Iterations to: ", mE, mI

      maxErr = mE
      maxIter = mI
      lambdaStart = la
      lambdaFactor = la

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
      ! todo: tet_cen, M, Happ
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
         ALLOCATE(vertex(nvertex,3),tet(ntet,4),state_dex(ntet), &
                  state_type(nstate),constant_mu(nstate), &
                  tet_cen(3,ntet),M(3,ntet),Happ(3,ntet),STAT=istat)
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
      call MUMATERIAL_SETD(1.0d-5, 100, 0.7d0, 0.75d0)

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











      SUBROUTINE mumaterial_iterate_magnetization(lambdaStart, lambdaFac)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: lambdaStart, lambdaFac
      INTEGER :: i, i_tile, j_tile, count, lambdaCount
      LOGICAL :: run, first
      DOUBLE PRECISION :: M_old(3), Mnorm(ntet), Mnorm_old(ntet), H(3), N_store(3,3,ntet,ntet)
      DOUBLE PRECISION :: H_old(3), H_new(3), M_new(3), lambda, lambda_s, error, maxDiff(4)

      ! Get initial magnetization
      ! DO i = 1, ntet
      !       SELECT CASE (state_type(state_dex(i)))
      !       CASE (1) ! Hard magnet
      !             WRITE(*,*) "Hard magnet not yet implemented. Is it even necessary?" ! todo
      !       CASE (2) ! Soft magnet using state function
      !             !CALL mumaterial_getM_soft(Happ(:,i), stateFunction(?), M(:,i)) ! todo
      !       CASE (3) 
      !             M(:,i) = Happ(:,i) * (constant_mu(state_dex(i)) - 1)
      !       CASE DEFAULT
      !             WRITE(*,*) "Unkown magnet type" ! todo stop?
      !             STOP
      !       END SELECT
      !       Mnorm(i) = NORM2(M(:,i))
      ! END DO
      count = 0;
      lambda = lambdaStart
      lambdaCount = 0;
      Mnorm_old = 0.d0;
      maxDiff = 0.d0;

      ! Iterate on magnetization
      first = .TRUE.
      DO
            ! Get the field and new magnetization for each tile
            DO i_tile = 1, ntet
                  M_old = M(:,i_tile)
                  H = Happ(:,i_tile)

                  ! Get the field from all other tiles
                  DO j_tile = 1, ntet
                        IF (first) THEN
                              CALL mumaterial_getN(vertex(tet(j_tile,1),:), vertex(tet(j_tile,2),:), vertex(tet(j_tile,3),:), vertex(tet(j_tile,4),:), tet_cen(:,i_tile), N_store(:,:,i_tile,j_tile))
                        END IF

                        IF (i_tile .ne. j_tile) THEN ! Current tile is not included in H
                              H = H + MATMUL(N_store(:,:,i_tile,j_tile),M(:,j_tile))
                        END IF
                  END DO

                  SELECT CASE (state_type(state_dex(i_tile)))
                  CASE (1) ! Hard magnet

                  CASE (2) ! Soft magnet using state function

                  CASE (3) ! Soft magnet using constant permeability
                        lambda_s = MIN(1/constant_mu(state_dex(i_tile)), 0.5)
                        H_new = H
                        DO
                              H_old = H_new
                              H_new = H + (constant_mu(state_dex(i_tile)) - 1) * MATMUL(N_store(:,:,i_tile,i_tile), H_new)
                              H_new = H_old + lambda_s * (H_new - H_old)

                              IF (MAXVAL(ABS((H_new - H_old)/H_old)) .lt. 0.0001*lambda_s) THEN
                                    M_new = MATMUL(N_store(:,:,i_tile,i_tile), H_new)
                                    EXIT
                              END IF
                        END DO
                  CASE DEFAULT
                        WRITE(*,*) "Unknown magnet type: ", state_type(state_dex(i_tile))
                        STOP
                  END SELECT

                  M(:,i_tile) = M_old + lambda*(M_new - M_old)
                  Mnorm(i_tile) = NORM2(M(:,i_tile))
            END DO
            
            
            first = .FALSE.
            count = count + 1

            error = 0.d0
            DO i = 1, ntet   
                  IF (Mnorm_old(i) .ne. 0.d0 .AND. ABS((Mnorm(i) - Mnorm_old(i))/Mnorm_old(i)) .gt. error) THEN
                        error = ABS((Mnorm(i) - Mnorm_old(i))/Mnorm_old(i))
                  END IF
            END DO



            WRITE(*,'(A7,I5,A8,E15.7,A13,E15.7)') 'Count: ', count, ' Error: ', error, ' Max. Error: ', maxErr*lambda

            IF (count .gt. 2 .AND. (error .lt. maxErr*lambda .OR. count .gt. maxIter)) THEN
                  EXIT
            ELSE
                  maxDiff = CSHIFT(maxDiff, 1)
                  maxDiff(1) = error

                  ! Update lambda if there is any increase in the error
                  IF (lambdaCount .gt. 4 .AND. ((maxDiff(2) - maxDiff(1)) .gt. 0.d0 .OR. (maxDiff(3) - maxDiff(2)) .gt. 0.d0 .OR. (maxDiff(4) - maxDiff(3)) .gt. 0.d0)) THEN
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
      

      IF (MAXVAL(ABS(offset)) .gt. 0.d0) THEN
            DO ik = 1, nvertex
                  vertex(ik,:) = vertex(ik,:) + offset
            END DO
      END IF

      DO ik = 1, ntet
            tet_cen(:,ik) = (vertex(tet(ik,1),:) + vertex(tet(ik,2),:) + vertex(tet(ik,3),:) + vertex(tet(ik,4),:))/4.d0
            
            CALL getBfld(tet_cen(1,ik), tet_cen(2,ik), tet_cen(3,ik), Bx, By, Bz)
            Happ(:,ik) = [Bx/mu0, By/mu0, Bz/mu0]

            ! IF (tiles(ik)%magnetType == 1) THEN
            !    tiles(ik)%mu_r_ea = constant_mu(state_dex(ik))
            !    tiles(ik)%mu_r_oa = constant_mu(state_dex(ik))

            ! ! todo change this to use actual remanent magnetisation
            !    ! set remanent magnetization to strength of the field [A/m] and easy axis in the direction of the field
            !    norm = sqrt(Bx*Bx + By*By + Bz*Bz)
            !    Bx_n = Bx / norm
            !    By_n = By / norm
            !    Bz_n = Bz / norm
            !    tiles(ik)%Mrem = norm / mu0
            !    tiles(ik)%u_ea = [Bx_n, By_n, Bz_n]
               
            !    ! get orthogonal axes
            !    IF (By/=0 .or. Bz/=0) THEN          ! cross product of u_ea with [1, 0, 0] and cross product of u_ea with cross product
            !          tiles(ik)%u_oa1 = [0.d0, Bz_n, -By_n]
            !          tiles(ik)%u_oa2 = [-By_n*By_n - Bz_n*Bz_n, Bx_n*By_n, Bx_n*Bz_n]
            !    ELSE                                ! cross product of u_ea with [0, 1, 0] and cross product of u_ea with cross product
            !          tiles(ik)%u_oa1 = [-Bz_n, 0.d0, Bx_n]
            !          tiles(ik)%u_oa2 = [Bx_n*By_n, -Bx_n*Bx_n - Bz_n*Bz_n, By_n*Bz_n]
            !    END IF
      END DO
            
      WRITE(*,*) "Running iterations"      
      
      ! allocate(stateFunction(1)) ! todo temporary until statefunction is properly implemented
      ! call loadStateFunctionFortran(stateFunction)


      ! CALL iterateMagnetization(tiles, ntet, stateFunction, size(stateFunction), temp, maxErr, maxIter, 0.d0) ! todo replace size?

      CALL mumaterial_iterate_magnetization(lambdaStart, lambdaFactor)

      RETURN
      END SUBROUTINE mumaterial_init

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
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x, y, z
      DOUBLE PRECISION, INTENT(out) :: Bx, By, Bz
      DOUBLE PRECISION :: H(3), mu0, N(3,3)
      INTEGER :: i

      mu0 = 16 * atan(1.d0) * 1.d-7
      H = 0.d0

      DO i = 1, ntet
            CALL mumaterial_getN(vertex(tet(i,1),:), vertex(tet(i,2),:), vertex(tet(i,3),:), vertex(tet(i,4),:), [x, y, z], N)
            H = H + MATMUL(N, M(:,i))
      END DO

      Bx = H(1) * mu0
      By = H(2) * mu0
      Bz = H(3) * mu0

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
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x(:), y(:), z(:)
      DOUBLE PRECISION, INTENT(out), ALLOCATABLE :: Bx(:), By(:), Bz(:)
      INTEGER :: n_points, i

      n_points = size(x)

      allocate(Bx(n_points))
      allocate(By(n_points))
      allocate(Bz(n_points))

      DO i = 1, n_points
            CALL mumaterial_getb_scalar(x(i), y(i), z(i), Bx(i), By(i), Bz(i))
      END DO
      
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
      IMPLICIT NONE
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

      CALL mumaterial_getb(x, y, z, Bx, By, Bz)

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
