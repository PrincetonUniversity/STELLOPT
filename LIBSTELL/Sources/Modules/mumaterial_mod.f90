!-----------------------------------------------------------------------
!     Module:        mumaterial_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov), Björn Hamstra
!     Date:          October 2023
!     Description:   This module is designed to help calculate the
!                    magnetic field arrising from ferromagnetic
!                    material
!-----------------------------------------------------------------------
      MODULE mumaterial_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE safe_open_mod
      IMPLICIT NONE

!-----------------------------------------------------------------------
!     Types    
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     Module Variables
!           lverb:            Controls output to screen
!           nvertex:          Number of vertices
!           ntet:             Number of tetrahedrons
!           nstate:           Number of state functions
!           tiles:            Magnetic tiles
!           stateFunction:    Array of state functions
!           maxErr:           Max allowed error for convergence in all iterations
!           maxIter:          Max allowed number of iterations
!           maxNb:            Max number of neighbours, <=0 disables nearest neighbours
!           lambdaStart:      Base value for lambda in iterate
!           lambdaFactor:     Mutiplication factor for lambda
!           vertex:           Vertices [m] (3,nvertex)
!           tet:              Tetrahedron (4,ntet) index into vertex
!           state_dex:        State function for each tetrahedron
!           state_type:       Type of state function (nstate) (1-3)
!           constant_mu:      Mu values for constant permeability (nstate)
!           constant_mu_o:    Mu values for orthogonal axis for hard magnet (nstate)
!           Mrem:             Magnet M Vector for hard magnet (3,nstate)
!           M:                Magnetization for all tetrahedrons (3,ntet)
!           Happ:             Applied magnetic field [A/m] at the centre of all tetrahedrons (3,ntet)
!-----------------------------------------------------------------------
      TYPE stateFunctionType
            DOUBLE PRECISION, PRIVATE, ALLOCATABLE :: H(:), M(:)
      END TYPE stateFunctionType

      LOGICAL, PRIVATE                    :: lverb
      INTEGER, PRIVATE                    :: nvertex, ntet, nstate
      TYPE(stateFunctionType), PRIVATE, ALLOCATABLE     :: stateFunction(:)
      DOUBLE PRECISION, PRIVATE           :: maxErr, lambdaStart, lambdaFactor
      INTEGER, PRIVATE                    :: maxIter, maxNb


      DOUBLE PRECISION, POINTER, PRIVATE :: vertex(:,:), tet_cen(:,:)
      INTEGER, POINTER, PRIVATE :: tet(:,:)
      INTEGER, POINTER, PRIVATE :: state_dex(:)
      INTEGER, POINTER, PRIVATE :: state_type(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: constant_mu(:), constant_mu_o(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: M(:,:), Happ(:,:), Mrem(:,:)
      INTEGER, PRIVATE            :: win_vertex, win_tet, win_tet_cen, &
                                     win_state_dex, win_state_type, &
                                     win_constant_mu, win_m, win_Mrem, &
                                     win_Happ, win_constant_mu_o

      CHARACTER(LEN=256), PRIVATE :: machine_string
      CHARACTER(LEN=256), PRIVATE :: date


      INTEGER, PRIVATE                    :: mystart, myend, mydelta
      INTEGER, PRIVATE                    :: shar_rank, shar_comm




      
!-----------------------------------------------------------------------
!     Subroutines
!         mumaterial_setd: Sets default values
!         mumaterial_load: Loads a magnetic material file
!         mumaterial_init: Calculates magnetization of material
!         mumaterial_getb: Calculates magnetic field in space
!             mumaterial_getb_scalar: Calculates magnetic field at a point in space
!             mumaterial_getb_vector: Calculates magnetic field for multiple points in space
!         mumaterial_output: Output tiles, H-field and points to file
!-----------------------------------------------------------------------
!     Functions
!-----------------------------------------------------------------------
      INTERFACE mumaterial_getb
            MODULE PROCEDURE mumaterial_getb_scalar, mumaterial_getb_vector
      END INTERFACE
      CONTAINS
      

      SUBROUTINE mumaterial_setverb(lverbin)
      !-----------------------------------------------------------------------
      ! mumaterial_setverb: Sets Verbosity
      !-----------------------------------------------------------------------
      ! param[in]: lverbin. Verbosity on
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: lverbin
      lverb = lverbin
      RETURN
      END SUBROUTINE mumaterial_setverb
      

      SUBROUTINE mumaterial_free()
      !-----------------------------------------------------------------------
      ! mumaterial_free: Deallocates memory
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ik
      IF (ASSOCIATED(state_dex))     CALL free_mpi_array1d_int(win_state_dex,state_dex,.TRUE.)
      IF (ASSOCIATED(state_type))    CALL free_mpi_array1d_int(win_state_type,state_type,.TRUE.)
      IF (ASSOCIATED(constant_mu))   CALL free_mpi_array1d_dbl(win_constant_mu,constant_mu,.TRUE.)
      IF (ASSOCIATED(constant_mu_o)) CALL free_mpi_array1d_dbl(win_constant_mu_o,constant_mu_o,.TRUE.)
      IF (ASSOCIATED(tet))           CALL free_mpi_array2d_int(win_tet,tet,.TRUE.)
      IF (ASSOCIATED(vertex))        CALL free_mpi_array2d_dbl(win_vertex,vertex,.TRUE.)
      IF (ASSOCIATED(tet_cen))       CALL free_mpi_array2d_dbl(win_tet_cen,tet_cen,.TRUE.)
      IF (ASSOCIATED(M))             CALL free_mpi_array2d_dbl(win_M,M,.TRUE.)
      ! TODO: Remove once allocated locally (Make sure code works beforehand)
      IF (ASSOCIATED(Mrem))          CALL free_mpi_array2d_dbl(win_Mrem,Mrem,.TRUE.)
      ! TODO: Remove once allocated locally (Make sure code works beforehand)
      IF (ASSOCIATED(Happ))          CALL free_mpi_array2d_dbl(win_Happ,Happ,.TRUE.)
      DO ik = 1, nstate
         IF (ALLOCATED(stateFunction(ik)%H)) DEALLOCATE(stateFunction(ik)%H)
         IF (ALLOCATED(stateFunction(ik)%M)) DEALLOCATE(stateFunction(ik)%M)
      END DO
      IF (ALLOCATED(stateFunction)) DEALLOCATE(stateFunction)
      RETURN
      END SUBROUTINE mumaterial_free


      SUBROUTINE mumaterial_setup(comm, shar_comm, comm_master, istat)
      !-----------------------------------------------------------------------
      ! mumaterial_setup: Sets up communicators for mumaterial
      !-----------------------------------------------------------------------
      ! param[in]: comm. World communicator
      ! param[out]: shar_comm. Shared memory communicator for calculations
      ! param[out]: comm_master. Master communicator handles cross-node stuff
      !-----------------------------------------------------------------------            
      IMPLICIT NONE

      INTEGER, INTENT(inout) :: comm, istat
      INTEGER, INTENT(out) :: shar_comm, comm_master
      INTEGER :: comm_myworld
      INTEGER :: shar_rank
      INTEGER :: color

      CALL MPI_COMM_DUP( comm, comm_myworld, istat )
      CALL MPI_COMM_SPLIT_TYPE( comm_myworld, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shar_comm, istat)
      CALL MPI_COMM_RANK( shar_comm, shar_rank, istat)

      color = MPI_UNDEFINED
      IF (mysharrank.eq.0) color = 0
      CALL MPI_COMM_SPLIT( comm_myworld, color, shar_rank, comm_master )
      RETURN
      END SUBROUTINE mumaterial_setup

      SUBROUTINE mumaterial_setd(mE, mI, la, laF, mNb)
      !-----------------------------------------------------------------------
      ! mumaterial_setd: Sets default values
      !-----------------------------------------------------------------------
      ! param[in]: mE. New maxErr: max error for MagTense convergence
      ! param[in]: mI. New maxIter: max amount of MagTense iterations
      ! param[in]: T. New temp: temperature of magnetic material in MagTense
      ! MPI should not be necessary here, every process calls this subroutine
      !-----------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in) :: mE, la, laF
      INTEGER, INTENT(in) :: mI, mNb

      maxErr = mE
      maxIter = mI
      lambdaStart = la
      lambdaFactor = laF
      maxNb = MIN(mNb,ntet-1)
      IF (mNb .le. 0) maxNb = ntet-1

      RETURN
      END SUBROUTINE mumaterial_setd


      SUBROUTINE mumaterial_load(filename,istat,shar_comm,comm_master)
      !-----------------------------------------------------------------------
      ! mumaterial_load: Loads magnetic material file.
      !-----------------------------------------------------------------------
      ! param[in]: filename. The file name to load in
      ! param[in, out]: istat. Integer that shows error if != 0
      ! param[in, out]: comm. MPI communicator, handles shared memory
      !-----------------------------------------------------------------------
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER, INTENT(inout)       :: istat
      INTEGER, INTENT(inout), OPTIONAL :: shar_comm, comm_master
      INTEGER :: shar_rank, master_rank
      LOGICAL :: lverb_temp, lcomm
      INTEGER :: iunit ,ik, i, j, nMH

      lcomm = (PRESENT(shar_comm).and.PRESENT(comm_master))
      shar_rank = 0
      master_rank = 0

      ! initialize MPI
#if defined(MPI_OPT)
    IF (lcomm) THEN
        CALL MPI_COMM_SPLIT_TYPE(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shar_comm, istat)
        CALL MPI_COMM_RANK( shar_comm, shar_rank, istat )
        IF (shar_rank.eq.0) CALL MPI_COMM_RANK( comm_master, master_rank, istat )
    END IF
#endif

      NULLIFY(vertex, tet, tet_cen, state_dex, state_type, constant_mu, &
              constant_mu_o, Mrem, M, Happ)

      ! open file, return if fails
      iunit = 327; istat = 0
      CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
      IF (istat/= 0) RETURN
      ! master reads info
      IF ((shar_rank.eq.0).and.(master_rank.eq.0)) THEN
         READ(iunit,'(A)') machine_string
         READ(iunit,'(A)') date
         READ(iunit,*) nvertex, ntet, nstate
      END IF

      ! Broadcast info to MPI and allocate vertex and face info
#if defined(MPI_OPT)
      IF (lcomm) THEN
            IF (shar_rank.eq.0) THEN ! master 0 broadcasts to other masters
                  CALL MPI_Bcast(nvertex,1,MPI_INTEGER,0,comm_master,istat)
                  CALL MPI_Bcast(ntet,   1,MPI_INTEGER,0,comm_master,istat)
                  CALL MPI_Bcast(nstate, 1,MPI_INTEGER,0,comm_master,istat)
            END IF
            ! every sharmem master broadcasts to other sharmem processes
            CALL MPI_Bcast(nvertex,1,MPI_INTEGER,0,shar_comm,istat)
            CALL MPI_Bcast(ntet,   1,MPI_INTEGER,0,shar_comm,istat)
            CALL MPI_Bcast(nstate, 1,MPI_INTEGER,0,shar_comm,istat)
            ! allocate on every sharmem island
            CALL mpialloc_2d_dbl(vertex,3,nvertex,    shar_rank,0,shar_comm,win_vertex)
            CALL mpialloc_2d_int(tet,4,ntet,          shar_rank,0,shar_comm,win_tet)
            CALL mpialloc_2d_dbl(tet_cen,3,ntet,      shar_rank,0,shar_comm,win_tet_cen)
            CALL mpialloc_1d_int(state_dex,ntet,      shar_rank,0,shar_comm,win_state_dex)
            CALL mpialloc_1d_int(state_type,nstate,   shar_rank,0,shar_comm,win_state_type)
            CALL mpialloc_1d_dbl(constant_mu,nstate,  shar_rank,0,shar_comm,win_constant_mu)
            CALL mpialloc_1d_dbl(constant_mu_o,nstate,shar_rank,0,shar_comm,win_constant_mu_o)
            ! TODO: Allocate locally
            CALL mpialloc_2d_dbl(Mrem,3,ntet,         shar_rank,0,shar_comm,win_Mrem)       
            CALL mpialloc_2d_dbl(M,3,ntet,            shar_rank,0,shar_comm,win_m)
            ! TODO: Allocate locally
            CALL mpialloc_2d_dbl(Happ,3,ntet,         shar_rank,0,shar_comm,win_Happ)       
            ALLOCATE(stateFunction(nstate))
      ELSE
#endif
         ! if no MPI, allocate everything on one node
         ALLOCATE(vertex(3,nvertex),tet(4,ntet),state_dex(ntet), &
                  state_type(nstate),constant_mu(nstate), &
                  tet_cen(3,ntet),M(3,ntet),Happ(3,ntet), &
                  constant_mu_o(nstate),Mrem(3,ntet),stateFunction(nstate), &
                  STAT=istat)
#if defined(MPI_OPT)
      END IF
#endif
      ! read in the mesh
      IF (istat/=0) RETURN
      
      IF ((shar_rank.EQ.0).AND.(master_rank.EQ.0)) THEN
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
               READ(iunit,*) stateFunction(ik)%H(:)
               READ(iunit,*) stateFunction(ik)%M(:)
            ELSEIF (state_type(ik) == 3) THEN
               READ(iunit,*) constant_mu(ik)
            ELSE
               PRINT *, '!!! UNKNOWN STATE_TYPE == ',state_type(ik)
            END IF
         END DO
      END IF

#if defined(MPI_OPT)
      IF ((lcomm).AND.(shar_rank.EQ.0)) THEN
            CALL MPI_Bcast(vertex,       3*nvertex,MPI_DOUBLE_PRECISION,0,comm_master,istat)
            CALL MPI_Bcast(tet,          4*ntet,   MPI_INTEGER,         0,comm_master,istat)
            CALL MPI_Bcast(state_dex,    ntet,     MPI_INTEGER,         0,comm_master,istat)
            CALL MPI_Bcast(state_type,   nstate,   MPI_INTEGER,         0,comm_master,istat)
            CALL MPI_Bcast(constant_mu,  nstate,   MPI_DOUBLE_PRECISION,0,comm_master,istat)
            CALL MPI_Bcast(constant_mu_o,nstate,   MPI_DOUBLE_PRECISION,0,comm_master,istat)
            ! TODO: Remove once allocated locally (make sure code works beforehand)
            CALL MPI_Bcast(Mrem,         3*ntet,   MPI_DOUBLE_PRECISION,0,comm_master,istat)
      END IF
      IF (lcomm) THEN ! Transfer state functions
            DO ik = 1, nstate
                  IF (shar_rank.EQ.0) THEN ! First from master to managers
                        IF (master_rank.EQ.0) THEN
                            IF (ALLOCATED(stateFunction(ik)%H)) THEN
                                nMH = SIZE(stateFunction(ik)%H)
                            ELSE
                                nMH = -1
                            END IF
                        END IF
                        CALL MPI_Bcast(nMH,1,MPI_INTEGER,0,comm_master,istat)
                        IF (nMH .gt. 0) THEN
                            IF (master_rank .ne. 0) ALLOCATE(stateFunction(ik)%H(nMH),stateFunction(ik)%M(nMH))
                            CALL MPI_Bcast(stateFunction(ik)%H,nMH,MPI_DOUBLE_PRECISION,0,comm_master,istat)
                            CALL MPI_Bcast(stateFunction(ik)%M,nMH,MPI_DOUBLE_PRECISION,0,comm_master,istat)
                        END IF
                  END IF ! Then from managers to processes
                  CALL MPI_Bcast(nMH,1,MPI_INTEGER,0,shar_comm,istat)
                  IF (nMH .gt. 0) THEN
                        IF (shar_rank .ne. 0) ALLOCATE(stateFunction(ik)%H(nMH),stateFunction(ik)%M(nMH))
                        CALL MPI_Bcast(stateFunction(ik)%H,nMH,MPI_DOUBLE_PRECISION,0,shar_comm,istat)
                        CALL MPI_Bcast(stateFunction(ik)%M,nMH,MPI_DOUBLE_PRECISION,0,shar_comm,istat)
                  END IF
            END DO
      END IF
#endif

      ! close file
      CLOSE(iunit)

      ! set default values
      lverb_temp = lverb
      lverb = .FALSE.
      CALL MUMATERIAL_SETD(1.0d-5, 100, 0.7d0, 0.75d0, 0)
      lverb = lverb_temp
      RETURN
      END SUBROUTINE mumaterial_load


      SUBROUTINE mumaterial_info(iunit)
      !-----------------------------------------------------------------------
      ! mumaterial_info: Prints info to iunit
      !-----------------------------------------------------------------------
      ! param[in]: iunit. Unit number to print to
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iunit
      INTEGER :: i,k
      WRITE(iunit,'(A)')         ' -----  Magnetic Material  -----'
      WRITE(iunit,'(3X,A,A)')    'Model Name   : ',TRIM(machine_string)
      WRITE(iunit,'(3X,A,A)')    'Date         : ',TRIM(date)
      WRITE(iunit,'(3X,A,I7)')   'Vertices     : ',nvertex
      WRITE(iunit,'(3X,A,I7)')   'Tetrahedrons : ',ntet
      WRITE(iunit,'(3X,A,I7)')   'State Funcs. : ',nstate
      WRITE(iunit,'(3X,A,I7)')   'Max Iter.    : ',maxIter
      WRITE(iunit,'(3X,A,EN12.3)') 'Max Error    : ',maxErr
      WRITE(iunit,'(3X,A,EN12.3)') 'Lambda       : ',lambdaStart
      WRITE(iunit,'(3X,A,EN12.3)') 'Lambda_fact  : ',lambdaFactor
      IF (maxNb > 0) THEN
         WRITE(iunit,'(3X,A,I7)')   'N-N Tiles    : ',maxNb
      ELSE
         WRITE(iunit,'(3X,A)')      'Full Model Run'
      END IF
      DO i = 1, nstate
         WRITE(iunit,'(6X,A,I3)') 'State Fuction ',i
         IF (state_type(i)==1) THEN
            WRITE(iunit,'(9X,A)') 'Type: Hard Magnet'
            WRITE(iunit,'(9X,A,EN12.3)')    '  Mu   :',constant_mu(i)
            WRITE(iunit,'(9X,A,EN12.3)')    '  Mu_o :',constant_mu_o(i)
            WRITE(iunit,'(9X,A,3(EN12.3))') '  Mrem :',Mrem(:,i)
         ELSEIF (state_type(i)==2) THEN
            k = SIZE(stateFunction(i)%H)
            WRITE(iunit,'(9X,A)')           '  Type : Soft Magnet (H-M)'
            WRITE(iunit,'(9X,A,I3)')        'NKnots :',k
            WRITE(iunit,'(9X,A,2(EN12.3))') '     H :',stateFunction(i)%H(1),stateFunction(i)%H(k)
            WRITE(iunit,'(9X,A,2(EN12.3))') '     M :',stateFunction(i)%M(1),stateFunction(i)%M(k)
         ELSEIF (state_type(i)==3) THEN
            WRITE(iunit,'(9X,A)') 'Type: Soft Magnet (mu constant)'
            WRITE(iunit,'(9X,A,EN12.3)')    '    Mu :',constant_mu(i)
         ELSE
            WRITE(iunit,'(9X,A,I3)') 'Type: UNKNOWN (ERROR) state_type=',state_type(i)
         END IF
      END DO
      END SUBROUTINE mumaterial_info



      SUBROUTINE mumaterial_init_new(getBfld, comm_world, shar_comm, comm_master, offset)
      !-----------------------------------------------------------------------
      ! mumaterial_init: Calculates magnetization of material
      !-----------------------------------------------------------------------
      ! fcn           : getBfld. Function which returns the vacuum magnetic field
      !                 SUBROUTINE FCN(x,y,z,bx,by,bz)
      ! param[in]: offset. Offset of all tiles from the origin
      ! param[in, out]: comms. MPI communicator, handles shared memory
      !-----------------------------------------------------------------------
#if defined(MPI_OPT)
      USE mpi
      USE mpi_params
#endif
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: offset(3)
      INTEGER, INTENT(inout), OPTIONAL :: comm_world, shar_comm, comm_master
      INTEGER :: shar_rank, master_rank
      LOGICAL :: lcomm
      INTEGER :: i, j, k, istat

      INTEGER, DIMENSION(:,:), POINTER :: neighbours
      LOGICAL, ALLOCATABLE :: mask(:)
      DOUBLE PRECISION :: x, y, z, Bx, By, Bz, Bx_n, By_n, Bz_n, mu0
      DOUBLE PRECISION, DIMENSION(:,:,:,:), POINTER :: N_store
      DOUBLE PRECISION, ALLOCATABLE :: dist(:), dx(:,:)

      EXTERNAL:: getBfld
      
      mu0 = 16.0D-7 * ATAN(1.d0)
      shar_rank = 0
      master_rank = 0
      lcomm = ((PRESENT(shar_comm).AND.PRESENT(comm_master)).AND.PRESENT(comm_world))

#if defined(MPI_OPT)
      IF (lcomm) THEN
         CALL MPI_COMM_RANK( shar_comm, shar_rank, istat )
         IF (shar_rank.EQ.0) CALL MPI_COMM_RANK( comm_master, master_rank, istat )
      END IF
#endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Apply offset
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (PRESENT(offset)) THEN
         IF (MAXVAL(ABS(offset)) .gt. 0.d0) THEN
            IF (lverb) WRITE(6,*) "  MUMAT_INIT:  Applying offset"
            mystart = 1; myend = nvertex
#if defined(MPI_OPT)
            IF (lcomm) CALL MPI_CALC_MYRANGE(comm_world, 1, nvertex, mystart, myend)
#endif
            DO i = mystart, myend
               vertex(:,i) = vertex(:,i) + offset
            END DO
            
#if defined(MPI_OPT)
            IF (lcomm) THEN
               CALL mumaterial_sync_array2d_dbl(vertex,3,nvertex,comm_master,shar_comm,mystart,myend,istat)
            END IF
#endif
         END IF
      END IF

#if defined(MPI_OPT)
      IF (lcomm) THEN
         CALL MPI_CALC_MYRANGE(shar_comm, 1, ntet, mystart, myend)
         tet_cen(:,mystart:myend) = 0.0
         ! TODO: Allocate locally [ ALLOCATE(Happ(3,mystart:myend)) ]
         ! since Happ is only required for loop over process assigned tets
         ! (saves time and allocated memory)
         ! BEFORE THAT: make sure code works without this change
         Happ(:,mystart:myend) = 0.0    
         CALL MPI_BARRIER( shar_comm,istat )
      END IF
#endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Calculate the Applied H-Field (do over all nodes)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (lverb) WRITE(6,*) "  MUMAT_INIT:  Calculating Fields at Tet. centers"
      mystart = 1; myend = ntet
#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_CALC_MYRANGE(comm_world, 1, ntet, mystart, myend)
#endif
      DO i = mystart, myend
        tet_cen(:,i) = (vertex(:,tet(1,i)) + vertex(:,tet(2,i)) + vertex(:,tet(3,i)) + vertex(:,tet(4,i)))/4.d0
        CALL getBfld(tet_cen(1,i), tet_cen(2,i), tet_cen(3,i), Bx, By, Bz)
        Happ(:,i) = [Bx/mu0, By/mu0, Bz/mu0]
      END DO

#if defined(MPI_OPT)
      IF (lcomm) THEN
        CALL mumaterial_sync_array2d_dbl(tet_cen,3,ntet,comm_master,shar_comm,mystart,myend,istat)
        ! TODO: Allocate locally, then remove this line
        CALL mumaterial_sync_array2d_dbl(Happ,   3,ntet,comm_master,shar_comm,mystart,myend,istat)
      END IF
#endif

      ! From here on out each thread only works on its own subset
      IF (lcomm) CALL MPI_CALC_MYRANGE(shar_comm, 1, ntet, mystart, myend)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Allocate helpers and Neighbors
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NULLIFY(N_store,neighbours)
      k = maxNb+1
      ALLOCATE(neighbours(maxNb,mystart:myend))
      ALLOCATE(N_store(3,3,maxNb+1,mystart:myend))
      neighbours(:,:) = 0
      N_store(:,:,:,:) = 0.0
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate nearest neighbors
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (lverb) WRITE (6,*) "  MUMAT_INIT:  Determining nearest neighbours"
      ALLOCATE(mask(ntet),dist(ntet),dx(3,ntet))
      DO i = mystart, myend
         ! We need to define helper variables dx(3,ntet)
         dx(1,:) = tet_cen(1,:)-tet_cen(1,i)
         dx(2,:) = tet_cen(2,:)-tet_cen(2,i)
         dx(3,:) = tet_cen(3,:)-tet_cen(3,i)
         dist = NORM2(dx,DIM=1)
         mask = .TRUE.
         mask(i) = .FALSE.
         DO j = 1, maxNb
            k = MINLOC(dist,1,mask)
            neighbours(j,i) = k
            mask(k) = .FALSE.
         END DO
      END DO
      DEALLOCATE(mask,dist,dx)

      IF (lverb) WRITE (6,*) "  MUMAT_INIT:  Calculating N_tensor"
      DO i = mystart, myend
         N_store(:,:,:,i) = 0
         CALL mumaterial_getN(vertex(:,tet(1,i)), vertex(:,tet(2,i)), vertex(:,tet(3,i)), vertex(:,tet(4,i)), tet_cen(:,i), N_store(:,:,maxNb+1,i))
      END DO

      ! Calculate N-tensors of neighbours
      DO i = mystart, myend      
         DO j = 1, maxNb     
            CALL mumaterial_getN(vertex(:,tet(1,neighbours(j,i))), vertex(:,tet(2,neighbours(j,i))), vertex(:,tet(3,neighbours(j,i))), vertex(:,tet(4,neighbours(j,i))), tet_cen(:,i), N_store(:,:,j,i)) 
         END DO               
      END DO

      ! No longer necessary
      CALL free_mpi_array2d_dbl(win_tet_cen,tet_cen,.TRUE.)

      IF (lverb) WRITE (6,*) "  MUMAT_INIT:  Beginning Iterations"
      IF (lcomm) THEN
            CALL mumaterial_iterate_magnetization_new(maxNb, mystart, myend, neighbours, N_store, shar_comm, comm_master)
      ELSE
            CALL mumaterial_iterate_magnetization_new(maxNb, mystart, myend, neighbours, N_store)
      END IF
      IF (lverb) WRITE (6,*) "  MUMAT_INIT:  End Iterations"

      ! DEALLOCATE Helpers
      DEALLOCATE(neighbours)
      DEALLOCATE(N_store)

      RETURN
      END SUBROUTINE mumaterial_init_new


      SUBROUTINE mumaterial_iterate_magnetization_new(N1, iA, iB, neighbours, N_store, shar_comm, comm_master)
      !-----------------------------------------------------------------------
      ! mumaterial_iterate_magnetization: Iterates the magnetic field over all tiles, called by mumaterial_init
      !-----------------------------------------------------------------------
      ! param[in]: N_store. Storage for demagnetization tensors from nearest neighbours
      ! param[in]: neighbours. Indices of nearest neighbours
      ! param[in, out]: comm. MPI communicator, handles shared memory
      !-----------------------------------------------------------------------
#if defined(MPI_OPT)
      USE mpi
      USE mpi_params
#endif
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N1, iA, iB
      INTEGER, OPTIONAL :: neighbours(N1,iA:iB)
      DOUBLE PRECISION, OPTIONAL :: N_store(3,3,N1+1,iA:iB)
      INTEGER, INTENT(inout), OPTIONAL :: shar_comm, comm_master
      INTEGER :: shar_rank

      LOGICAL :: lcomm
      INTEGER :: count, i_tile, j_tile, lambdaCount, istat, iC, iD
      DOUBLE PRECISION, DIMENSION(:), POINTER :: chi, Mnorm, Mnorm_old
      DOUBLE PRECISION, DIMENSION(:,:), POINTER :: M_new

      DOUBLE PRECISION :: H(3), N(3,3)
      DOUBLE PRECISION :: H_old(3), H_new(3), lambda, lambda_s, error, maxDiff(4), Hnorm, M_tmp_norm
      DOUBLE PRECISION :: M_tmp(3), M_tmp_local(3), Mrem_norm, u_ea(3), u_oa_1(3), u_oa_2(3) ! hard magnet

      lcomm = (PRESENT(shar_comm).AND.PRESENT(comm_master))
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         ! Get extent of shared memory area
         CALL MPI_COMM_RANK( shar_comm, shar_rank, istat )
         iC = iA
         iD = iB
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, iC, 1, MPI_INTEGER, MPI_MIN, shar_comm, istat)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, iD, 1, MPI_INTEGER, MPI_MAX, shar_comm, istat)
         PRINT *,shar_rank,iA,iB,iC,iD
      END IF
#endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Allocate Helper Arrays
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE(M_new(3,iA:iB),chi(iA:iB),Mnorm(iA:iB),Mnorm_old(iA:iB))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Defaults
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      count = 0
      lambda = lambdaStart
      lambdaCount = 0
      maxDiff = 0.d0
      chi = 0.0
      Mnorm = 0.0
      Mnorm_old = 1.0E-5 ! Initial value is large
      M(:,iA:iB) = 0.0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Main Iteration Loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO
         M_new = 0.0
         error = 0.d0
         ! Get the field and new magnetization for each tile
         DO i_tile = iA, iB
            H = Happ(:,i_tile)
            ! Get the field from all other tiles
            DO j_tile = 1, maxNb
               H = H + MATMUL(N_store(:,:,j_tile,i_tile), M(:,neighbours(j_tile,i_tile)))
            END DO
            N = N_store(:,:,maxNb+1,i_tile)  

            ! Determine field and magnetization at tile due to all other tiles and itself
            H_new = H
            SELECT CASE (state_type(state_dex(i_tile)))
               CASE (1) ! Hard magnet
                  Mrem_norm = NORM2(Mrem(:,state_dex(i_tile)))
                  u_ea = Mrem(:,state_dex(i_tile))/Mrem_norm ! Easy axis assumed parallel to remanent magnetization
                  IF (u_ea(2)/=0 .OR. u_ea(3)/=0) THEN      ! Cross product of u_ea with [1, 0, 0] and cross product of u_ea with cross product
                     u_oa_1 = [0.d0, u_ea(3), -u_ea(2)]
                     u_oa_2 = [-u_ea(2)*u_ea(2) - u_ea(3)*u_ea(3), u_ea(1)*u_ea(2), u_ea(1)*u_ea(3)]
                  ELSE                                      ! Cross product of u_ea with [0, 1, 0] and cross product of u_ea with cross product
                     u_oa_1 = [-u_ea(3), 0.d0, u_ea(1)]
                     u_oa_2 = [u_ea(1)*u_ea(2), -u_ea(1)*u_ea(1) - u_ea(3)*u_ea(3), u_ea(2)*u_ea(3)]
                  END IF

                  ! Normalize unit vectors
                  u_oa_1 = u_oa_1/NORM2(u_oa_1)
                  u_oa_2 = u_oa_2/NORM2(u_oa_2)
                     
                  lambda_s = MIN(1/constant_mu(state_dex(i_tile)), 1/constant_mu_o(state_dex(i_tile)), 0.5)
                  DO
                     H_old = H_new
                     ! Determine magnetization taking into account easy axis
                     M_tmp = (Mrem_norm + (constant_mu(state_dex(i_tile))   - 1) * DOT_PRODUCT(H_new, u_ea)) * u_ea &
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
                  WRITE(6,*) "Unknown magnet type: ", state_type(state_dex(i_tile))
                  STOP
            END SELECT
            ! Moved from further down
            M(:,i_tile) = M(:,i_tile) + lambda*(M_new(:,i_tile) - M(:,i_tile))
            Mnorm(i_tile) = NORM2(M(:,i_tile))
            error = MAX(ABS((Mnorm(i_tile) - Mnorm_old(i_tile))/Mnorm_old(i_tile)),error)
         END DO

#if defined(MPI_OPT)
         ! Synchronise M
         IF (lcomm) THEN
            CALL mumaterial_sync_array2d_dbl(M,3,ntet,comm_master,shar_comm,iA,iB,istat)
            IF (shar_rank.EQ.0) CALL MPI_ALLREDUCE(MPI_IN_PLACE, error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_master, istat)
            CALL MPI_BARRIER( shar_comm, istat)
         END IF
#endif

         Mnorm_old = Mnorm

         ! Determine change in magnetization
         count = count + 1
         IF (lverb) WRITE(6,'(5X,A7,I5,A8,E15.7,A13,E15.7)') 'Count: ', count, ' Error: ', error, ' Max. Error: ', maxErr*lambda
         CALL FLUSH(6)

         IF (count .gt. 2 .AND. (error .lt. maxErr*lambda .OR. count .gt. maxIter)) EXIT
            
         maxDiff = CSHIFT(maxDiff, -1)
         maxDiff(1) = error

         ! Update lambda if there is any increase in the error (maxDiff(1) is the most recent)
         IF (lambdaCount .gt. 5 .AND. ((maxDiff(2) - maxDiff(1)) .lt. 0.d0 .OR. (maxDiff(3) - maxDiff(2)) .lt. 0.d0 .OR. (maxDiff(4) - maxDiff(3)) .lt. 0.d0)) THEN
               lambda = lambda * lambdaFactor
               lambdaCount = 1
         !ELSE
            !lambda = MIN(lambdaStart,lambda/lambdaFactor)
         END IF
         lambdaCount = lambdaCount + 1
      END DO

      DEALLOCATE(M_new,chi,Mnorm,Mnorm_old)

      RETURN
      END SUBROUTINE mumaterial_iterate_magnetization_new


      SUBROUTINE mumaterial_getN(v1, v2, v3, v4, pos, N)
      !-----------------------------------------------------------------------
      ! mumaterial_getN: Helper function to determine the demagnetization tensor
      !-----------------------------------------------------------------------
      ! param[in]: v1-4: Vertices of the tetrahedron (3)
      ! param[in]: pos: Reference position for which to determine the demagnetization tensor (3)
      ! param[in]: N. Resulting demagnetization tensor (3,3)
      !-----------------------------------------------------------------------
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

            N_loc(1,3) = mumaterial_getNxz(r, v(1,1), v(2,2)) - mumaterial_getNxz(r, v(1,3), v(2,2))
            N_loc(2,3) = mumaterial_getNyz(r, v(1,1), v(2,2)) - mumaterial_getNyz(r, v(1,3), v(2,2))
            N_loc(3,3) = mumaterial_getNzz(r, v(1,1), v(2,2)) - mumaterial_getNzz(r, v(1,3), v(2,2))

            N = N + MATMUL(MATMUL(P, N_loc), Pinv)
      END DO

      RETURN
      END SUBROUTINE mumaterial_getN

      FUNCTION mumaterial_getNxz(r, l, h)
      !-----------------------------------------------------------------------
      ! mumaterial_getNxz: Helper function to determine the x-component of the demagnetization tensor
      !-----------------------------------------------------------------------
      ! param[in]: r: Reference position for which to determine the demagnetization tensor (3)
      ! param[in]: l: Bottom side of the triangle 
      ! param[in]: h. Top side of the triangle
      !-----------------------------------------------------------------------
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
      !-----------------------------------------------------------------------
      ! mumaterial_getNyz: Helper function to determine the y-component of the demagnetization tensor
      !-----------------------------------------------------------------------
      ! param[in]: r: Reference position for which to determine the demagnetization tensor (3)
      ! param[in]: l: Bottom side of the triangle 
      ! param[in]: h. Top side of the triangle
      !-----------------------------------------------------------------------
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
      !-----------------------------------------------------------------------
      ! mumaterial_getNzz: Helper function to determine the z-component of the demagnetization tensor
      !-----------------------------------------------------------------------
      ! param[in]: r: Reference position for which to determine the demagnetization tensor (3)
      ! param[in]: l: Bottom side of the triangle 
      ! param[in]: h. Top side of the triangle
      !-----------------------------------------------------------------------
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

        SUBROUTINE mumaterial_sync_array2d_dbl(array, n1, n2, comm_master, shar_comm, &
            ourstart, ourend, istat)
        IMPLICIT NONE

        INTEGER, INTENT(in) :: n1, n2
        DOUBLE PRECISION, DIMENSION(n1,n2), INTENT(inout) :: array
        INTEGER, INTENT(inout) :: comm_master, shar_comm
        INTEGER, INTENT(in) :: ourstart, ourend
        INTEGER :: shar_rank, istat
        INTEGER :: start, end,

        ! Zero array "above" data to keep
        IF (ourstart.NE.1) THEN
            CALL MPI_CALC_MYRANGE(shar_comm, 1, ourstart-1, start, end)
            array(:,mystart:myend) = 0
        END IF
        ! Zero array "below" data to keep
        IF (ourend.NE.n1) THEN
            CALL MPI_CALC_MYRANGE(shar_comm,ourend+1,n1, start, end)
            array(:,mystart:myend) = 0
        END IF

        ! Make sure the array is properly modified on the shared memory
        CALL MPI_BARRIER( shar_comm, istat ) 
        CALL MPI_COMM_RANK( shar_comm, shar_rank, istat )
        IF (shar_rank.eq.0) THEN
            ! Make sure the other masters are also done with their array
            CALL MPI_BARRIER( comm_master, istat )
            ! Finally, reduce arrays onto all shared memory islands
            CALL MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, istat )
        END IF
        CALL MPI_BARRIER( shar_comm, istat ) 

      END SUBROUTINE mumaterial_sync_array2d_dbl

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
      ! Assume fx is non-decreasing
      IF (x .lt. fx(1)) THEN
            y = fy(1)
            RETURN
      ELSEIF (x .gt. fx(n)) THEN
            y = fy(n)
            RETURN
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
      ! fcn           : getBfld. Function which returns the vacuum magnetic field
      !                 SUBROUTINE FCN(x,y,z,bx,by,bz)
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


      


      SUBROUTINE mumaterial_getbmag_scalar(x, y, z, Bx, By, Bz)
      !-----------------------------------------------------------------------
      ! mumaterial_getbmag: Calculates total magnetic field at a point in space
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

      mu0 = 16.0D-7 * atan(1.d0)
      H = 0.d0

      DO i = 1, ntet
            CALL mumaterial_getN(vertex(:,tet(1,i)), vertex(:,tet(2,i)), vertex(:,tet(3,i)), vertex(:,tet(4,i)), [x, y, z], N)
            H = H + MATMUL(N, M(:,i))
      END DO

      Bx = H(1) * mu0
      By = H(2) * mu0
      Bz = H(3) * mu0

      RETURN
      END SUBROUTINE mumaterial_getbmag_scalar


      SUBROUTINE mumaterial_getb_vector(x, y, z, B, getBfld, comm_world, shar_comm, comm_master)
      !-----------------------------------------------------------------------
      ! mumaterial_getb_vector: Calculates total magnetic field at multiple points in space
      !-----------------------------------------------------------------------
      ! param[in]: x. x-coordinates of points at which to determine the magnetic field
      ! param[in]: y. y-coordinates of points at which to determine the magnetic field
      ! param[in]: z. z-coordinates of points at which to determine the magnetic field
      ! param[out]: B.  B-field at required points [T]
      ! fcn           : getBfld. Function which returns the vacuum magnetic field
      !                 SUBROUTINE FCN(x,y,z,bx,by,bz)
      !-----------------------------------------------------------------------
#if defined(MPI_OPT)
      USE mpi
      USE mpi_params
#endif
      IMPLICIT NONE
      EXTERNAL:: getBfld
      DOUBLE PRECISION, INTENT(in) :: x(:), y(:), z(:)
      DOUBLE PRECISION, INTENT(out), ALLOCATABLE :: B(:,:)
      INTEGER, INTENT(inout), OPTIONAL :: comm_world, shar_comm, comm_master
      LOGICAL :: lcomm
      INTEGER :: i, istat 
      INTEGER :: npoints

      shar_rank = 0;
      lcomm = ((PRESENT(comm_world).AND.PRESENT(shar_comm)).AND.PRESENT(comm_master))

      npoints = size(x)
      mystart = 1; myend = npoints
#if defined(MPI_OPT)
         IF (lcomm) CALL MPI_CALC_MYRANGE(comm_world, 1, npoints, mystart, myend)
#endif

      allocate(B(3,npoints)))
      B = 0

      DO i = mystart, myend
            CALL mumaterial_getb_scalar(x(i), y(i), z(i), B(1,i), B(2,i), B(3,i), getBfld)
      END DO

#if defined(MPI_OPT)
      IF (lcomm) THEN
        CALL mumaterial_sync_array2d_dbl(B,3,npoints,comm_master,shar_comm,mystart,myend,istat)
      END IF
#endif
      
      RETURN
      END SUBROUTINE mumaterial_getb_vector





      SUBROUTINE mumaterial_output(path, x, y, z, getBfld, comm_world, shar_comm, comm_master)
      !-----------------------------------------------------------------------
      ! mumaterial_output: Outputs B-field and points to text files
      !-----------------------------------------------------------------------
      ! param[in]: path. Path to store files in
      ! param[in]: x. x-cooridinates of points at which to determine the magnetic field
      ! param[in]: y. y-cooridinates of points at which to determine the magnetic field
      ! param[in]: z. z-cooridinates of points at which to determine the magnetic field
      ! fcn           : getBfld. Function which returns the vacuum magnetic field
      !                 SUBROUTINE FCN(x,y,z,bx,by,bz)
      !-----------------------------------------------------------------------
#if defined(MPI_OPT)
      USE mpi
#endif      
      IMPLICIT NONE
      EXTERNAL:: getBfld
      CHARACTER(LEN=*), INTENT(in) :: path
      DOUBLE PRECISION, INTENT(in) :: x(:), y(:), z(:)
      INTEGER, INTENT(inout), OPTIONAL :: comm_world, shar_comm, comm_master
      INTEGER :: i, shar_rank, world_rank, master_rank, istat 
      LOGICAL :: lcomm, lismaster
      INTEGER :: npoints
      DOUBLE PRECISION, ALLOCATABLE :: B(:,:)

      shar_rank = 0; master_rank = 0; world_rank = 0; 
      lcomm = ((PRESENT(comm_world).AND.PRESENT(shar_comm)).AND.PRESENT(comm_master))
      lismaster = .FALSE.
      IF (NOT(lcomm)) lismaster = .TRUE.

#if defined(MPI_OPT)
      IF (lcomm) THEN
         CALL MPI_COMM_RANK( shar_comm, shar_rank, istat )
         IF (shar_rank.EQ.0) THEN
            CALL MPI_COMM_RANK( comm_master, master_rank, istat)
            IF (master_rank.EQ.0) lismaster = .TRUE.
         END IF
      END IF
#endif

      IF (lismaster) THEN
            npoints = size(x)
            WRITE(6,*) "Outputting points"
            OPEN(13, file=TRIM(path)//'/points.dat')
            DO i = 1, npoints
                  WRITE(13, "(F15.7,A,F15.7,A,F15.7)") x(i), ',', y(i), ',', z(i)
            END DO
            CLOSE(13)
      END IF

      CALL mumaterial_getb_vector(x, y, z, B, getBfld, comm_world, shar_comm, comm_master)
 

      IF (lismaster) THEN
            WRITE(6,*) "Outputting B-field"
            OPEN(14, file=TRIM(path)//'/B.dat')
            DO i = 1, npoints
                  WRITE(14, "(E15.7,A,E15.7,A,E15.7)") B(1,i), ',', B(2,i), ',', B(3,i)
            END DO
            CLOSE(14)
      END IF

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

      SUBROUTINE mpialloc_4d_int(array,n1,n2,n3,n4,subid,mymaster,share_comm,win)
      ! Libraries
      USE MPI
      USE ISO_C_BINDING
      IMPLICIT NONE
      ! Arguments
      INTEGER, POINTER, INTENT(inout) :: array(:,:,:,:)
      INTEGER, INTENT(in) :: n1
      INTEGER, INTENT(in) :: n2
      INTEGER, INTENT(in) :: n3
      INTEGER, INTENT(in) :: n4
      INTEGER, INTENT(in) :: subid
      INTEGER, INTENT(in) :: mymaster
      INTEGER, INTENT(in) :: share_comm
      INTEGER, INTENT(inout) :: win
      ! Variables
      INTEGER :: disp_unit, ier
      INTEGER :: array_shape(4)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: window_size
      TYPE(C_PTR) :: baseptr
      ! Initialization
      ier = 0
      array_shape(1) = n1
      array_shape(2) = n2
      array_shape(3) = n3
      array_shape(4) = n4
      disp_unit = 1
      window_size = 0_MPI_ADDRESS_KIND
      IF (subid == mymaster) window_size = INT(n1*n2*n3*n4,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
      IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
      CALL C_F_POINTER(baseptr, array, array_shape)
      RETURN
      END SUBROUTINE mpialloc_4d_int

      SUBROUTINE mpialloc_4d_dbl(array,n1,n2,n3,n4,subid,mymaster,share_comm,win)
      ! Libraries
      USE MPI
      USE ISO_C_BINDING
      IMPLICIT NONE
      ! Arguments
      DOUBLE PRECISION, POINTER, INTENT(inout) :: array(:,:,:,:)
      INTEGER, INTENT(in) :: n1
      INTEGER, INTENT(in) :: n2
      INTEGER, INTENT(in) :: n3
      INTEGER, INTENT(in) :: n4
      INTEGER, INTENT(in) :: subid
      INTEGER, INTENT(in) :: mymaster
      INTEGER, INTENT(in) :: share_comm
      INTEGER, INTENT(inout) :: win
      ! Variables
      INTEGER :: disp_unit, ier
      INTEGER :: array_shape(4)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: window_size
      TYPE(C_PTR) :: baseptr
      ! Initialization
      ier = 0
      array_shape(1) = n1
      array_shape(2) = n2
      array_shape(3) = n3
      array_shape(4) = n4
      disp_unit = 1
      window_size = 0_MPI_ADDRESS_KIND
      IF (subid == mymaster) window_size = INT(n1*n2*n3*n4,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
      IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
      CALL C_F_POINTER(baseptr, array, array_shape)
      RETURN
      END SUBROUTINE mpialloc_4d_dbl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Memory Freeing Subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
         SUBROUTINE free_mpi_array1d_int(win_local,array_local,isshared)
         IMPLICIT NONE
         LOGICAL, INTENT(in) :: isshared
         INTEGER, INTENT(inout) :: win_local
         INTEGER, POINTER, INTENT(inout) :: array_local(:)
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
         END SUBROUTINE free_mpi_array1d_int
   
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
   
         SUBROUTINE free_mpi_array2d_int(win_local,array_local,isshared)
         IMPLICIT NONE
         LOGICAL, INTENT(in) :: isshared
         INTEGER, INTENT(inout) :: win_local
         INTEGER, POINTER, INTENT(inout) :: array_local(:,:)
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
         END SUBROUTINE free_mpi_array2d_int
   
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
   
         SUBROUTINE free_mpi_array4d_int(win_local,array_local,isshared)
         IMPLICIT NONE
         LOGICAL, INTENT(in) :: isshared
         INTEGER, INTENT(inout) :: win_local
         INTEGER, POINTER, INTENT(inout) :: array_local(:,:,:,:)
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
         END SUBROUTINE free_mpi_array4d_int
   
         SUBROUTINE free_mpi_array4d_dbl(win_local,array_local,isshared)
         IMPLICIT NONE
         LOGICAL, INTENT(in) :: isshared
         INTEGER, INTENT(inout) :: win_local
         DOUBLE PRECISION, POINTER, INTENT(inout) :: array_local(:,:,:,:)
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
         END SUBROUTINE free_mpi_array4d_dbl

!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE mumaterial_mod
