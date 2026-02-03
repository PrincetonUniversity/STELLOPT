!------------------------------------------------------------------------------
!     Module:        mumaterial_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov), Björn Hamstra,
!                    Lucas van Ham (lucas.van.ham@ipp.mpg.de)
!     Date:          October 2023
!                    Jan-June 2024 [LvH]: Modifications for scaled-up problems,
!                    neighbour formulation, dipole approximation
!     Description:   This module calculates the magnetic response of 
!                    ferromagnetic materials, modelled as a mesh of 
!                    tetrahedrons, to magnetic fields.
!------------------------------------------------------------------------------
      MODULE mumaterial_mod
!------------------------------------------------------------------------------
!     Libraries
!------------------------------------------------------------------------------
      USE safe_open_mod
#if defined(MPI_OPT)
      USE mpi
      USE mpi_params
      USE mpi_sharmem, ONLY: mpialloc, mpidealloc
#endif
      IMPLICIT NONE
!------------------------------------------------------------------------------
!     Types    
!       StateFunctionType: State function for soft magnets. Contains:
!         H: Array of H-field values for material
!         M: Array of corresponding M values        
!------------------------------------------------------------------------------
      TYPE stateFunctionType
            DOUBLE PRECISION, PRIVATE, ALLOCATABLE :: H(:), M(:)
      END TYPE stateFunctionType    

      PROCEDURE(externalFieldFunc), POINTER :: getBfld
      ABSTRACT INTERFACE
        SUBROUTINE externalFieldFunc(x,y,z,Bx,By,Bz)
          DOUBLE PRECISION, INTENT(in)  :: x,y,z
          DOUBLE PRECISION, INTENT(out) :: Bx,By,Bz
        END SUBROUTINE externalFieldFunc
      END INTERFACE
!------------------------------------------------------------------------------
!     Module Variables
!        lverb:      Controls output to screen
!   
!       MPI
!         lcomm:     .TRUE. if code is run with MPI
!         ldosync:   .TRUE. if more than one MPI node is used
!         lismaster: .TRUE. if rank of thread in world-communicator is 0
!
!         comm_shar:   Shared-memory communicator
!         master_comm: Communicator of threads whose rank in comm_shar is 0
!         world_comm:  Communicator with all MPI threads
!         color:       Used to create master_comm    
!         COMM_rank:   Rank of thread in communicator COMM
!         COMM_size:   Number of threads in communicator COMM
!         
!         dom_shar:    Collection of tetrahedrons worked on by comm_shar
!         outmydom: Collection of tetrahedrons NOT worked on by comm_shar
!         ntet_shar:  Size of dom_shar
!         win_OBJ:  MPI shared memory window for OBJ
!
!       Neighbours
!         nbrs:         Array of neighbours for each tetrahedron (:,:)
!         nbrs_count:        Number of neighbours for each tetrahedron (:)
!         nbrs_maxc:     Largest neighbour count in nbrs_count
!         Nb_domidx:  Neighbours indexed by appearance in dom_shar (:,:)
!
!       Mesh
!         ntet:     Number of tetrahedrons in mesh
!         nvertex:  Number of vertices in mesh
!         vertex:   Coordinates for vertices in mesh (3, ntet)
!         tet:      Vertex indices for each tetrahedron (4, ntet)
!         tet_cen:  Coordinates for tetrahedron centers (3, ntet)
!         tet_vol:  Volumes of tetrahedrons (ntet)
!         tet_rad: Equivalent length of edge of tetrahedrons (ntet)
!   
!       Magnetics
!         nstate:           Number of state functions
!         state_dex:        State function for each tetrahedron (ntet)
!         state_type:       Type of state function (nstate) (1-3)
!         constant_mu:      Mu for constant permeability (nstate)
!         constant_mu_o:    Mu for orthogonal axis for hard magnet (nstate)
!         Mrem:             Remanent magnetization for hard magnet (3,nstate)
!         M:                Magnetization for all tetrahedrons (3,ntet)
!         H_app:             Applied H-field at tetrahedron centres (3, ntet)
!         MU0:              Permeability of free space: 4*pi*1E-7 [H/m]
!         N_store:          Demagnetization tensor (3,3,nbrs_maxc,:)
!
!       User settings
!         threshold:           Threshold error for convergence
!         maxIter:         Max allowed number of iterations
!         padFactor:   Affects number of neighbours for each tetrahedron
!         lambdaStart:     Initial value of lambda_n for iterations
!         lambdaFactor:    Multiplication factor for lambda_n
!         lambdaThresh:    Multiply lambda_n if error grows this number of times
!------------------------------------------------------------------------------

      CHARACTER(LEN=256), PRIVATE :: machine_string, file_string
      CHARACTER(LEN=256), PRIVATE :: date

      ! mesh variables
      INTEGER, PRIVATE  ::  ntet, nvertex
      DOUBLE PRECISION, POINTER, PRIVATE :: vertex(:,:), tet_cen(:,:), & 
                                            tet_vol(:), tet_rad(:), &
                                            r_cluster(:,:), d_cluster(:), &
                                            m_cluster(:,:), Q_cluster(:,:)
      INTEGER, POINTER, PRIVATE :: tet(:,:), dom_cluster(:,:)

      ! magnetics variables
      INTEGER, PRIVATE :: nstate
      INTEGER, POINTER, PRIVATE :: state_dex(:), state_type(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: constant_mu(:), constant_mu_o(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: M(:,:), Mrem(:,:)
      DOUBLE PRECISION, DIMENSION(:,:,:,:), POINTER, PRIVATE :: N_store
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: H_app
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: inv_mat_local
      TYPE(stateFunctionType), PRIVATE, ALLOCATABLE :: stateFunction(:)

      ! user settings variables
      DOUBLE PRECISION, PRIVATE :: threshold, padFactor, lambdaStart, lambdaFactor, convCheck
      INTEGER, PRIVATE          :: lambdaThresh, maxIter

      ! neighbour variables
      INTEGER, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: nbrs, Nb_domidx
      INTEGER, DIMENSION(:),   ALLOCATABLE, PRIVATE :: nbrs_count
      INTEGER, PRIVATE                              :: nbrs_maxc

      ! MPI variables
      INTEGER, PRIVATE :: comm_shar,   shar_rank,   shar_size, &
                          comm_master, master_rank, master_size, &
                          comm_world,  world_rank,  world_size, &
                          color
      LOGICAL, PRIVATE :: lcomm, lismaster, ldosync

      ! MPI windows
      INTEGER, PRIVATE :: win_vertex, win_tet, win_tet_cen, &
                          win_tet_vol, win_tet_rad,  &
                          win_state_dex, win_state_type, &
                          win_constant_mu, win_m, win_Mrem, &
                          win_Happ, win_constant_mu_o, & 
                          win_r_cluster, win_m_cluster, &
                          win_d_cluster, win_dom_cluster, &
                          win_Q_cluster
      ! box division variables
      INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: dom_shar, outmydom, dom_proc
      INTEGER, PRIVATE                            :: ntet_shar,odomsize, ntet_proc
      LOGICAL, DIMENSION(:), ALLOCATABLE, PRIVATE :: lisfar
      INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: dom_mid_proc
      INTEGER, PRIVATE                            :: ntet_mid_proc

      ! verbose
      LOGICAL, PRIVATE                    :: lverb

      ! precomputed constants
      DOUBLE PRECISION, PARAMETER, PRIVATE :: PI = 4.0D0*ATAN(1.0D0)
      DOUBLE PRECISION, PARAMETER, PRIVATE :: INVPI = 1.0D0/PI
      DOUBLE PRECISION, PARAMETER, PRIVATE :: INV4PI = 1.0D0/(4.0D0*PI)
      DOUBLE PRECISION, PARAMETER, PRIVATE :: INV8PI = 1.0D0/(8.0D0*PI)
      DOUBLE PRECISION, PARAMETER, PRIVATE :: MU0 = 4.0D-7*PI
      DOUBLE PRECISION, PARAMETER, PRIVATE :: small = 1E-12
!------------------------------------------------------------------------------
!     Subroutines
!       Main flow
!         mumaterial_setup:     Sets up MPI communicators (optional)
!         mumaterial_load:      Loads magnetic material file and sets up MPI stuff
!         mumaterial_setdefs:   Sets default values
!         mumaterial_setverb:   Sets standard verbosity
!         mumaterial_setBfld:   Sets function for external B field
!         mumaterial_info:      Prints information to screen
!         mumaterial_init:      Initializes everything, calls iteration subroutine
!         mumaterial_iterate_M: Main calculation loop
!
!       Helpers
!         mumaterial_gettetvolume:  Calculates volume of a tetrahedron
!         mumaterial_getneighbours: Determines tetrahedron neighbours
!         mumaterial_getN:          Determines demagnetization tensor
!           mumaterial_getNxz: x-component 
!           mumaterial_getNyz: y-component 
!           mumaterial_getNzz: z-component
!         mumaterial_cross:         Cross product of two vectors
!         mumaterial_getState:      Interpolates function 
!
!       MPI 
!         mumaterial_split:            Splits input domain into two subdomains
!         mumaterial_syncM: Syncs magnetization array on shar_mem nodes
!         mumaterial_free:  Frees MPI memory
!       Output
!         mumaterial_output:  Output B-field and points to file
!         mumaterial_getb:    Calculates magnetic field in space
!             mumaterial_getb_scalar:      Single point in space
!               mumaterial_getbmag_scalar: Excludes applied field
!             mumaterial_getb_vector: Multiple points in space
!
!       Debug
!         mumaterial_writedebug: Outputs files for debug
!------------------------------------------------------------------------------
!     Functions
!------------------------------------------------------------------------------
      INTERFACE mumaterial_getb
        MODULE PROCEDURE mumaterial_getb_scalar, mumaterial_getb_vector
      END INTERFACE
      CONTAINS
!------------------------------------------------------------------------------
! mumaterial_setverb: Sets Verbosity
!------------------------------------------------------------------------------
! param[in]: lverbin. Verbosity on
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_setverb(lverbin)

      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: lverbin

      lverb = lverbin
      RETURN

      END SUBROUTINE mumaterial_setverb
 
!------------------------------------------------------------------------------
! mumaterial_free: Deallocates memory and destroys MPI windows
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_free()

      IMPLICIT NONE

      INTEGER :: ik
      IF (ASSOCIATED(state_dex))     CALL mpidealloc(state_dex,win_state_dex)
      IF (ASSOCIATED(state_type))    CALL mpidealloc(state_type,win_state_type)
      IF (ASSOCIATED(constant_mu))   CALL mpidealloc(constant_mu,win_constant_mu)
      IF (ASSOCIATED(constant_mu_o)) CALL mpidealloc(constant_mu_o,win_constant_mu_o)
      IF (ASSOCIATED(tet))           CALL mpidealloc(tet,win_tet)
      IF (ASSOCIATED(vertex))        CALL mpidealloc(vertex,win_vertex)
      IF (ASSOCIATED(tet_cen))       CALL mpidealloc(tet_cen,win_tet_cen)
      IF (ASSOCIATED(tet_vol))       CALL mpidealloc(tet_vol,win_tet_vol)
      IF (ASSOCIATED(tet_rad))       CALL mpidealloc(tet_rad,win_tet_rad)
      IF (ASSOCIATED(M))             CALL mpidealloc(M,win_M)
      ! TODO: Remove once allocated locally (Make sure code works beforehand)
      IF (ASSOCIATED(Mrem))          CALL mpidealloc(Mrem,win_Mrem)
      IF (ASSOCIATED(r_cluster))     CALL mpidealloc(r_cluster,win_r_cluster)
      IF (ASSOCIATED(d_cluster))     CALL mpidealloc(d_cluster,win_d_cluster)
      IF (ASSOCIATED(dom_cluster))   CALL mpidealloc(dom_cluster,win_dom_cluster)
      IF (ASSOCIATED(m_cluster))     CALL mpidealloc(m_cluster,win_m_cluster)
      IF (ASSOCIATED(Q_cluster))     CALL mpidealloc(Q_cluster,win_Q_cluster)

      DO ik = 1, nstate
         IF (ALLOCATED(stateFunction(ik)%H)) DEALLOCATE(stateFunction(ik)%H)
         IF (ALLOCATED(stateFunction(ik)%M)) DEALLOCATE(stateFunction(ik)%M)
      END DO
      IF (ALLOCATED(stateFunction)) DEALLOCATE(stateFunction)

      RETURN
      END SUBROUTINE mumaterial_free

!------------------------------------------------------------------------------
! mumaterial_setup: Sets up communicators for mumaterial
!------------------------------------------------------------------------------
! param[in]:  comm. World communicator from which other comms are born
! param[out]: comm_shar_out. Shared memory communicator for calculations
! param[out]: comm_master_out. Master communicator handles cross-node stuff
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_setup(comm, comm_shar_out, comm_master_out)

      IMPLICIT NONE

      INTEGER, INTENT(inout) :: comm
      INTEGER, INTENT(out) :: comm_shar_out, comm_master_out
      INTEGER :: comm_myworld

      CALL MPI_COMM_DUP( comm, comm_myworld, ierr_mpi )
      CALL MPI_COMM_SPLIT_TYPE( comm_myworld, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, comm_shar_out, ierr_mpi)
      CALL MPI_COMM_RANK( comm_shar_out, shar_rank, ierr_mpi)

      color = MPI_UNDEFINED
      IF (shar_rank.EQ.0) color = 0
      CALL MPI_COMM_SPLIT( comm_myworld, color, shar_rank, comm_master_out, ierr_mpi )

      RETURN

      END SUBROUTINE mumaterial_setup

!------------------------------------------------------------------------------
!       mumaterial_setdefs: Sets default values
!------------------------------------------------------------------------------
! param[in]: mE. threshold: threshold for determining convergence
! param[in]: mI. maxIter: max amount of iterations
! param[in]: la. lambdaStart: initial value of lambda_n
! param[in]: laF. lambdaFactor: multiplicative factor for lambda_n
! param[in]: laT. lambdaThresh: amount of dM>0 before lambda_n is multiplied
! param[in]: padF. padFactor: factor for sphere around tets for neighbours
! param[in]: cc. convCheck: Stop when this percentage of elemnts has converged
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_setdefs(mE, mI, la, laF, laT, padF, cc)
 
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in) :: mE, la, laF, padF, cc
      INTEGER, INTENT(in) :: mI, laT

      threshold = mE
      maxIter = mI
      lambdaStart = la
      lambdaFactor = laF
      lambdaThresh = laT
      padFactor = padF
      convCheck = cc

      RETURN

      END SUBROUTINE mumaterial_setdefs

!------------------------------------------------------------------------------
! mumaterial_load: Loads magnetic material file and sets MPI defaults
!------------------------------------------------------------------------------
! param[in]: filename. The file name to load in
! param[in, out]: istat. Integer that shows  if != 0
! param[in, out]: comm_shar_in. MUMAT shared memory communicator
! param[in, out]: comm_master_in. MUMAT communicator of sharmem masters
! param[in, out]: comm_world_in. MUMAT world communicator
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_load(filename,istat,comm_shar_in,comm_master_in,comm_world_in)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER, INTENT(inout)       :: istat
      INTEGER, INTENT(inout), OPTIONAL :: comm_shar_in, comm_master_in, comm_world_in
      INTEGER :: iunit ,ik, i, j, nMH

      ! Set lcomm
      lcomm = ((PRESENT(comm_shar_in).and.PRESENT(comm_master_in)).AND.PRESENT(comm_world_in))
      IF (lcomm) THEN
        comm_shar   = comm_shar_in
        comm_master = comm_master_in 
        comm_world  = comm_world_in
      END IF

      ! Default parameters for no MPI
      world_rank = 0; shar_rank = 0; master_rank = 0 
      world_size = 1; shar_size = 1; master_size = 1
      lismaster = .TRUE.

      ! Set up MPI parameters properly now
#if defined(MPI_OPT)
      IF (lcomm) THEN
        lismaster = .FALSE.; master_rank = 1
        CALL MPI_COMM_RANK( comm_world, world_rank, ierr_mpi)
        CALL MPI_COMM_SIZE( comm_world, world_size, ierr_mpi)
        CALL MPI_COMM_RANK( comm_shar,  shar_rank,  ierr_mpi)
        CALL MPI_COMM_SIZE( comm_shar,  shar_size,  ierr_mpi)
        IF (shar_rank.eq.0) THEN
          CALL MPI_COMM_RANK( comm_master, master_rank, ierr_mpi )
          CALL MPI_COMM_SIZE( comm_master, master_size, ierr_mpi )
          lismaster = (master_rank.EQ.0)
        END IF
        CALL MPI_Bcast( master_size, 1, MPI_INTEGER, 0, comm_shar, ierr_mpi)
      END IF
#endif

      ! Nullify pointers
      NULLIFY(vertex, tet, tet_cen, tet_vol, tet_rad, state_dex, state_type, &
              constant_mu, constant_mu_o, Mrem, M, N_store, &
              r_cluster, m_cluster, d_cluster, dom_cluster, Q_cluster)

      ! open file, return if fails
      iunit = 327; istat = 0
      file_string = TRIM(filename)
      CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
      IF (istat/= 0) RETURN
      ! master reads info
      IF (lismaster) THEN
         READ(iunit,'(A)') machine_string
         READ(iunit,'(A)') date
         READ(iunit,*) nvertex, ntet, nstate
      END IF

      ! Broadcast info to MPI and allocate vertex and face info
#if defined(MPI_OPT)
      IF (lcomm) THEN
        IF (shar_rank.eq.0) THEN ! world master broadcasts to other masters
          CALL MPI_Bcast(nvertex,1,MPI_INTEGER,0,comm_master,ierr_mpi)
          CALL MPI_Bcast(ntet,   1,MPI_INTEGER,0,comm_master,ierr_mpi)
          CALL MPI_Bcast(nstate, 1,MPI_INTEGER,0,comm_master,ierr_mpi)
        END IF
        ! every sharmem master broadcasts to other sharmem processes
        CALL MPI_Bcast(nvertex,1,MPI_INTEGER,0,comm_shar,ierr_mpi)
        CALL MPI_Bcast(ntet,   1,MPI_INTEGER,0,comm_shar,ierr_mpi)
        CALL MPI_Bcast(nstate, 1,MPI_INTEGER,0,comm_shar,ierr_mpi)
        ! allocate on every sharmem island
        CALL mpialloc(vertex,3,nvertex,    shar_rank,0,comm_shar,win_vertex)
        CALL mpialloc(tet,4,ntet,          shar_rank,0,comm_shar,win_tet)
        CALL mpialloc(tet_cen,3,ntet,      shar_rank,0,comm_shar,win_tet_cen)
        CALL mpialloc(tet_vol,ntet,        shar_rank,0,comm_shar,win_tet_vol)
        CALL mpialloc(tet_rad,ntet,       shar_rank,0,comm_shar,win_tet_rad)
        CALL mpialloc(state_dex,ntet,      shar_rank,0,comm_shar,win_state_dex)
        CALL mpialloc(state_type,nstate,   shar_rank,0,comm_shar,win_state_type)
        CALL mpialloc(constant_mu,nstate,  shar_rank,0,comm_shar,win_constant_mu)
        CALL mpialloc(constant_mu_o,nstate,shar_rank,0,comm_shar,win_constant_mu_o)
        CALL mpialloc(M,            3,ntet,shar_rank,0,comm_shar,win_m)
        CALL mpialloc(Mrem,3,nstate,       shar_rank,0,comm_shar,win_Mrem)  ! TODO: Allocate locally
        CALL mpialloc(r_cluster,3,world_size,  shar_rank,0,comm_shar,win_r_cluster)
        CALL mpialloc(d_cluster,world_size,    shar_rank,0,comm_shar,win_d_cluster)
        CALL mpialloc(m_cluster,3,world_size,  shar_rank,0,comm_shar,win_m_cluster)
        CALL mpialloc(Q_cluster,6,world_size,  shar_rank,0,comm_shar,win_Q_cluster)

        ALLOCATE(stateFunction(nstate))
      ELSE
#endif
         ! if no MPI, allocate everything on one node
         ALLOCATE(vertex(3,nvertex),tet(4,ntet),state_dex(ntet), &
                  state_type(nstate),constant_mu(nstate), &
                  tet_cen(3,ntet),tet_vol(ntet),tet_rad(ntet),M(3,ntet), &
                  constant_mu_o(nstate),Mrem(3,nstate),stateFunction(nstate), &
                  r_cluster(3,1),d_cluster(1),dom_cluster(1,1), &
                  m_cluster(3,1),Q_cluster(6,1), &
                  STAT=istat)
#if defined(MPI_OPT)
      END IF
#endif
      IF (istat/=0) RETURN
      M(:,:) = 0.0

      ! read in the mesh
      IF (lismaster) THEN
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
        CALL MPI_Bcast(vertex,       3*nvertex,MPI_DOUBLE_PRECISION,0,comm_master,ierr_mpi)
        CALL MPI_Bcast(tet,          4*ntet,   MPI_INTEGER,         0,comm_master,ierr_mpi)
        CALL MPI_Bcast(state_dex,    ntet,     MPI_INTEGER,         0,comm_master,ierr_mpi)
        CALL MPI_Bcast(state_type,   nstate,   MPI_INTEGER,         0,comm_master,ierr_mpi)
        CALL MPI_Bcast(constant_mu,  nstate,   MPI_DOUBLE_PRECISION,0,comm_master,ierr_mpi)
        CALL MPI_Bcast(constant_mu_o,nstate,   MPI_DOUBLE_PRECISION,0,comm_master,ierr_mpi)
        CALL MPI_Bcast(Mrem,       3*nstate,   MPI_DOUBLE_PRECISION,0,comm_master,ierr_mpi) ! TODO: Remove once allocated locally (make sure code works beforehand)
      END IF
      IF (lcomm) THEN ! Transfer state functions
        DO ik = 1, nstate
          ! First from master to submasters
          IF (shar_rank.EQ.0) THEN 
            IF (master_rank.EQ.0) THEN
              IF (ALLOCATED(stateFunction(ik)%H)) THEN
                nMH = SIZE(stateFunction(ik)%H)
              ELSE
                nMH = -1
              END IF
            END IF
            CALL MPI_Bcast(nMH,1,MPI_INTEGER,0,comm_master,ierr_mpi)
            IF (nMH .gt. 0) THEN
              IF (master_rank .ne. 0) ALLOCATE(stateFunction(ik)%H(nMH),stateFunction(ik)%M(nMH))
              CALL MPI_Bcast(stateFunction(ik)%H,nMH,MPI_DOUBLE_PRECISION,0,comm_master,ierr_mpi)
              CALL MPI_Bcast(stateFunction(ik)%M,nMH,MPI_DOUBLE_PRECISION,0,comm_master,ierr_mpi)
            END IF
          END IF 
          ! Then from submasters to other threads
          CALL MPI_Bcast(nMH,1,MPI_INTEGER,0,comm_shar,ierr_mpi)
          IF (nMH .gt. 0) THEN
            IF (shar_rank .ne. 0) ALLOCATE(stateFunction(ik)%H(nMH),stateFunction(ik)%M(nMH))
            CALL MPI_Bcast(stateFunction(ik)%H,nMH,MPI_DOUBLE_PRECISION,0,comm_shar,ierr_mpi)
            CALL MPI_Bcast(stateFunction(ik)%M,nMH,MPI_DOUBLE_PRECISION,0,comm_shar,ierr_mpi)
          END IF
        END DO
      END IF
#endif

      ! close file
      CLOSE(iunit)

      ! set default values
      CALL mumaterial_setdefs(1.0d-5, 100, 0.7d0, 0.75d0, 10, 20.d0, 99.d0)

      RETURN

      END SUBROUTINE mumaterial_load

!------------------------------------------------------------------------------
! mumaterial_setBfld: associates external B field function getBfld with func_B
!------------------------------------------------------------------------------
! param[in]: func_B. External B field function.
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_setBfld(func_B)
        procedure(externalFieldFunc) :: func_B
        getBfld => func_B
      END SUBROUTINE mumaterial_setBfld

!------------------------------------------------------------------------------
! mumaterial_info: Prints info to iunit
!------------------------------------------------------------------------------
! param[in]: iunit. Unit number to print to
! param[in]: lnoiter: Whether to do mumat iterations
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_info(iunit, lnoiter)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iunit
      LOGICAL, INTENT(IN) :: lnoiter
      INTEGER :: i,k

          WRITE(iunit,'(A)')           ' -----  MUMAT calculation  -----'
          IF (lnoiter) THEN
            WRITE(iunit,'(3X,A)') '!!! SKIPPING ITERATIONS !!!'
          ELSE
            WRITE(iunit,'(3X,A,F9.3)')  'Pad factor   : ',padFactor
            WRITE(iunit,'(3X,A,I9)')    'Max Iter.    : ',maxIter
            WRITE(iunit,'(3X,A,ES9.3)') 'Max Error    : ',threshold
            WRITE(iunit,'(3X,A,F9.3)')  'Lambda start : ',lambdaStart
            WRITE(iunit,'(3X,A,F9.3)')  'Lambda fact. : ',lambdaFactor
            WRITE(iunit,'(3X,A,I9)')    'Lambda thrsh.: ',lambdaThresh
            WRITE(iunit,'(3X,A,F7.2,A)')'Converged at : ',convCheck,' %'
          END IF
          WRITE(iunit,'(A)')           ' -----  Magnetic structure  ----'
          WRITE(iunit,'(3X,A,A)')      'File: ',TRIM(file_string)
          WRITE(iunit,'(3X,A,A)')      'Model Name   : ',TRIM(machine_string)
          WRITE(iunit,'(3X,A,A)')      'Date         : ',TRIM(date)
          WRITE(iunit,'(3X,A,I9)')     'Vertices     : ',nvertex
          WRITE(iunit,'(3X,A,I9)')     'Tetrahedrons : ',ntet
          WRITE(iunit,'(3X,A,I9)')     'State Funcs. : ',nstate
          DO i = 1, nstate
            WRITE(iunit,'(5X,A,I0)') 'State Function ',i
            IF (state_type(i)==1) THEN
              WRITE(iunit,'(7X,A)') 'Type: Hard Magnet'
              WRITE(iunit,'(7X,A,ES12.3)')    '  Mu   :',constant_mu(i)
              WRITE(iunit,'(7X,A,ES12.3)')    '  Mu_o :',constant_mu_o(i)
              WRITE(iunit,'(7X,A,3(ES12.3))') '  Mrem :',Mrem(:,i)
            ELSEIF (state_type(i)==2) THEN
              k = SIZE(stateFunction(i)%H)
              WRITE(iunit,'(7X,A)')           '  Type : Soft Magnet (H-M)'
              WRITE(iunit,'(7X,A,I3)')        'NKnots :',k
              WRITE(iunit,'(7X,A,2(ES12.3))') '     H :',stateFunction(i)%H(1),stateFunction(i)%H(k)
              WRITE(iunit,'(7X,A,2(ES12.3))') '     M :',stateFunction(i)%M(1),stateFunction(i)%M(k)
            ELSEIF (state_type(i)==3) THEN
              WRITE(iunit,'(7X,A)') 'Type: Soft Magnet (mu constant)'
              WRITE(iunit,'(7X,A,F12.3)')    '    Mu :',constant_mu(i)
            ELSE
              WRITE(iunit,'(7X,A,I3)') 'Type: UNKNOWN (ERROR) state_type=',state_type(i)
            END IF
          END DO
        FLUSH(iunit)

      END SUBROUTINE mumaterial_info

!------------------------------------------------------------------------------
! mumaterial_writedebug: Writes files for debugging
!------------------------------------------------------------------------------
! param[in]: array. Array of values to write out. (n1,n2)
! param[in]: n1, n2. Dimensions of array.
! param[in]: filename. Filename to which program will write out.
! param[in]: displayname (optional). Will write this name to screen.
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_writedebug(array, n1, n2, filename, displayname)

      IMPLICIT NONE

      INTEGER, INTENT(in) :: n1, n2
      DOUBLE PRECISION, INTENT(in) :: array(n1,n2)
      CHARACTER(LEN=*), INTENT(in) :: filename
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: displayname
      INTEGER :: i

      IF (PRESENT(displayname)) WRITE(6,*) '  MUMAT_DEBUG: Outputting ' // TRIM(displayname)
      OPEN(15, file=TRIM(filename))
      IF (n1.EQ.3) THEN
        DO i = 1, n2
          WRITE(15, "(E15.7,A,E15.7,A,E15.7)") array(1,i), ',', array(2,i), ',', array(3,i)
        END DO
      ELSE IF (n1.EQ.1) THEN
        DO i = 1, n2
          WRITE(15, "(E15.7)")                 array(1,i)
        END DO
      END IF
      CLOSE(15)
        
      END SUBROUTINE mumaterial_writedebug

!------------------------------------------------------------------------------
! mumaterial_init: Initial calculations, does MPI, and calls iterations
!------------------------------------------------------------------------------
! param[in]: offset. Offset of all tiles from the origin
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_init_new(offset)

      IMPLICIT NONE
      
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: offset(3)
      INTEGER :: mystart, myend, ourstart, ourend, wr_dex
      LOGICAL :: lwork
      INTEGER :: i, j, k, i_tile, j_tile, stype
      INTEGER :: mstat(MPI_STATUS_SIZE)
      CHARACTER(LEN=6) :: strcount, splitcount

      DOUBLE PRECISION :: Bx, By, Bz
      DOUBLE PRECISION :: cutoff
      DOUBLE PRECISION :: targ
      DOUBLE PRECISION :: MAT(3,3)
      INTEGER :: splits, ydomsize, reci, src, n_proc_targ, a, b
      INTEGER, ALLOCATABLE :: dom_in(:), dom_out_2(:), dom_sizes(:),  &
                              mid_ids(:), cluster(:), temp_dom(:)
      LOGICAL, ALLOCATABLE :: mid_mask(:)
      INTEGER :: n, idx
      INTEGER :: ntet_proc_min, ntet_proc_max, ntet_shar_min, ntet_shar_max
      INTEGER :: nbrs_proc_min, nbrs_proc_max, ntet_mid_proc_min, ntet_mid_proc_max
      INTEGER :: tile1, tile2
      DOUBLE PRECISION, ALLOCATABLE :: d(:)
      DOUBLE PRECISION :: d_cluster_min,  d_cluster_max, d_max, d_worst
      DOUBLE PRECISION :: H_app_norm_min, H_app_norm_max
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Apply offset
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (PRESENT(offset).AND.(NORM2(offset) .GT. 0.d0)) THEN
        IF (lverb) WRITE(6,*) "  MUMAT_INIT:  Applying offset to vertices"
        mystart = 1; myend = nvertex
#if defined(MPI_OPT)
        IF (lcomm) CALL MPI_CALC_MYRANGE(comm_world, 1, nvertex, mystart, myend)
#endif
        DO i = mystart, myend
          vertex(:,i) = vertex(:,i) + offset
        END DO
            
#if defined(MPI_OPT)
        IF (shar_rank.EQ.0) THEN
          CALL MPI_ALLREDUCE( MPI_IN_PLACE, vertex,   3*nvertex, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
        END IF
#endif
      END IF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Calculate tet centers, edges, volumes, then synchronize
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (lverb) WRITE(6,*) "  MUMAT_INIT:  Calculating tetrahedron quantities"; FLUSH(6)
      ! Wipe shared array
      IF (shar_rank.EQ.0) THEN 
        tet_cen = 0  
        tet_vol = 0  
        tet_rad = 0     
      END IF  
      ! Calculate range
      IF (lcomm) THEN 
#if defined(MPI_OPT)
        CALL MPI_BARRIER(comm_world, ierr_mpi)
        CALL MPI_CALC_MYRANGE(comm_world, 1, ntet, mystart, myend) 
#endif
      ELSE
        mystart = 1
        myend = ntet
      END IF

      DO i = mystart, myend
          tet_cen(:,i) = (vertex(:,tet(1,i)) + vertex(:,tet(2,i)) + &
                          vertex(:,tet(3,i)) + vertex(:,tet(4,i))) / 4.d0
          tet_vol(i) = mumaterial_gettetvolume( &
              vertex(:,tet(1,i)),vertex(:,tet(2,i)), vertex(:,tet(3,i)),vertex(:,tet(4,i)))
          tet_rad(i) = SQRT(6.0)/12.d0*(6.d0*SQRT(2.0)*tet_vol(i))**(1.0/3.0)
      END DO

#if defined(MPI_OPT)
      IF (shar_rank.EQ.0) THEN
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tet_cen, 3*ntet, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tet_vol,   ntet, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tet_rad,   ntet, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
      ENDIF
#endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Split domain across MPI nodes
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      color = 0
#if defined(MPI_OPT)
      IF (shar_rank.NE.0) THEN
        color = 1
      ELSE
        Bx = master_size ! master_size needs to be dbl for next line
        color = world_rank*Bx/world_size
      END IF
#endif

      lwork = (color.EQ.0)
      IF (lwork) THEN ! Global master, create first box
        ALLOCATE(dom_shar(ntet))
        DO i = 1, ntet
          dom_shar(i) = i
        END DO
        ntet_shar = SIZE(dom_shar)
      END IF

#if defined(MPI_OPT)   
      IF ((master_size.GT.1).AND.(shar_rank.EQ.0)) THEN         
        splits = NINT(LOG(Bx)/LOG(2.0)) ! log_2(X) = ln(X)/log(2)
        targ = 0.5
        DO
          IF (splits.EQ.0) EXIT ! Reached end

          IF (lwork) THEN 
            ALLOCATE(dom_in(ntet_shar))
            dom_in = dom_shar
            DEALLOCATE(dom_shar)
            CALL mumaterial_split(dom_in,targ,dom_shar,dom_out_2) ! Split box
            DEALLOCATE(dom_in)
            ntet_shar  = SIZE(dom_shar)
            ydomsize = SIZE(dom_out_2)            
            splits = splits-1 

            ! now mail one of new boxes to the appropriate recipient
            reci = color + 2**splits 
            CALL MPI_SEND(ydomsize,      1, MPI_INTEGER, reci, 1234, comm_master, ierr_mpi) 
            CALL MPI_SEND(dom_out_2,ydomsize, MPI_INTEGER, reci, 1235, comm_master, ierr_mpi);  DEALLOCATE(dom_out_2) 
            CALL MPI_SEND(splits,        1, MPI_INTEGER, reci, 1236, comm_master, ierr_mpi)
          ELSE
            CALL MPI_RECV(ntet_shar,     1, MPI_INTEGER, MPI_ANY_SOURCE, 1234, comm_master, mstat, ierr_mpi); ALLOCATE(dom_shar(ntet_shar))
            CALL MPI_RECV(dom_shar, ntet_shar, MPI_INTEGER, MPI_ANY_SOURCE, 1235, comm_master, mstat, ierr_mpi)
            CALL MPI_RECV(splits,      1, MPI_INTEGER, MPI_ANY_SOURCE, 1236, comm_master, mstat, ierr_mpi);  
            lwork = .TRUE. ! Activate node
          END IF
        END DO

        ! Now every master needs to know their range relative to other 
        CALL MPI_SCAN(ntet_shar, ourend, 1, MPI_INTEGER, MPI_SUM, comm_master, ierr_mpi)
        ourstart = ourend-ntet_shar+1
      END IF

      CALL MPI_Bcast(ntet_shar,    1, MPI_INTEGER, 0, comm_shar, ierr_mpi)
      
      IF (shar_rank.NE.master) THEN
        ALLOCATE(dom_shar(ntet_shar))
      END IF
      
      CALL MPI_Bcast(dom_shar,   ntet_shar,  MPI_INTEGER, 0, comm_shar, ierr_mpi)
      CALL MPI_Bcast(ourstart, 1, MPI_INTEGER, 0, comm_shar, ierr_mpi)
      CALL MPI_Bcast(ourend,   1, MPI_INTEGER, 0, comm_shar, ierr_mpi)
      CALL MPI_BARRIER(comm_world, ierr_mpi)
#endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Assign each MPI thread a spatially localized cluster
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (shar_size.EQ.1) THEN
        ntet_proc = SIZE(dom_shar)
        ALLOCATE(dom_proc(ntet_proc))
        dom_proc = dom_shar
      END IF

#if defined(MPI_OPT)   
      IF (lcomm) THEN
        IF (shar_size.GT.1) THEN
          IF (shar_rank.EQ.master) THEN 
            ! Set up work initially
            n_proc_targ = shar_size
            a = n_proc_targ/2
            b = n_proc_targ - a
            targ = DBLE(a)/DBLE(n_proc_targ)
            reci = master + 1
            ! Divide box initially and send to two threads
            CALL mumaterial_split(dom_shar,targ,dom_proc,dom_out_2)
            ntet_proc = SIZE(dom_proc)
            CALL MPI_SEND(SIZE(dom_out_2),          1, MPI_INTEGER, reci, 101, comm_shar, ierr_mpi) 
            CALL MPI_SEND(dom_out_2,  SIZE(dom_out_2), MPI_INTEGER, reci, 102, comm_shar, ierr_mpi)
            CALL MPI_SEND(b,                        1, MPI_INTEGER, reci, 103, comm_shar, ierr_mpi)
            IF (shar_size.GT.2) THEN ! Careful with 2 threads; then need to keep one.
              reci = reci + 1
              CALL MPI_SEND(ntet_proc,       1, MPI_INTEGER, reci, 101, comm_shar, ierr_mpi) 
              CALL MPI_SEND(dom_proc,ntet_proc, MPI_INTEGER, reci, 102, comm_shar, ierr_mpi)
              CALL MPI_SEND(a,               1, MPI_INTEGER, reci, 103, comm_shar, ierr_mpi)
              DEALLOCATE(dom_proc)
            END IF
            DEALLOCATE(dom_out_2)
            ! Now just be the messenger
            IF (shar_size.GT.2) THEN
              DO
                reci = reci + 1
                ! Receive from any source
                CALL MPI_RECV(ntet_proc,        1, MPI_INTEGER, MPI_ANY_SOURCE, 101, comm_shar, mstat, ierr_mpi)
                src = mstat(MPI_SOURCE)
                ALLOCATE(dom_proc(ntet_proc))
                CALL MPI_RECV(dom_proc, ntet_proc, MPI_INTEGER, src, 102, comm_shar, mstat, ierr_mpi)
                CALL MPI_RECV(n_proc_targ,      1, MPI_INTEGER, src, 103, comm_shar, mstat, ierr_mpi)
                IF (reci.EQ.shar_size) THEN
                  EXIT ! That's us!
                ELSE
                  CALL MPI_SEND(ntet_proc,        1, MPI_INTEGER, reci, 101, comm_shar, ierr_mpi) 
                  CALL MPI_SEND(dom_proc, ntet_proc, MPI_INTEGER, reci, 102, comm_shar, ierr_mpi)
                  CALL MPI_SEND(n_proc_targ,      1, MPI_INTEGER, reci, 103, comm_shar, ierr_mpi)
                  DEALLOCATE(dom_proc)
                END IF

              END DO
            END IF
          ELSE
            ! non-master work
            lwork = .FALSE.
            DO
              ! Receive box from master if doesn't yet have
              IF (.NOT.lwork) THEN
                CALL MPI_RECV(ntet_proc,        1, MPI_INTEGER, master, 101, comm_shar, mstat, ierr_mpi) 
                ALLOCATE(dom_proc(ntet_proc))
                CALL MPI_RECV(dom_proc, ntet_proc, MPI_INTEGER, master, 102, comm_shar, mstat, ierr_mpi)       
                CALL MPI_RECV(n_proc_targ,      1, MPI_INTEGER, master, 103, comm_shar, mstat, ierr_mpi)
              END IF
              ! Check if should divide
              IF (n_proc_targ.EQ.1) THEN
                EXIT ! That's us!
              ELSE
                ! Divide box and send 
                a = n_proc_targ/2
                b = n_proc_targ - a
                targ = DBLE(a)/DBLE(n_proc_targ)
                ALLOCATE(temp_dom(ntet_proc))
                temp_dom = dom_proc
                DEALLOCATE(dom_proc)
                CALL mumaterial_split(temp_dom,targ,dom_proc,dom_out_2)
                ntet_proc = SIZE(dom_proc)
                CALL MPI_SEND(SIZE(dom_out_2),          1, MPI_INTEGER, master, 101, comm_shar, ierr_mpi) 
                CALL MPI_SEND(dom_out_2,  SIZE(dom_out_2), MPI_INTEGER, master, 102, comm_shar, ierr_mpi)
                CALL MPI_SEND(b,                        1, MPI_INTEGER, master, 103, comm_shar, ierr_mpi)
                DEALLOCATE(dom_out_2,temp_dom)

                n_proc_targ = a
                lwork = .TRUE.
              END IF
            END DO
          END IF
        END IF
        CALL MPI_BARRIER(comm_shar, ierr_mpi)

      END IF

#endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate and write cluster quantities
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      cutoff = 3.0
      r_cluster = 0.0
      d_cluster = 0.0
      m_cluster = 0.0
      Q_cluster = 0.0
      IF (lcomm) THEN
#if defined(MPI_OPT)
        ALLOCATE(dom_sizes(world_size))
        dom_sizes = 0
        wr_dex = world_rank+1 ! Pesky 0-based indexing
        dom_sizes(wr_dex) = ntet_proc
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, dom_sizes, world_size, MPI_INTEGER, MPI_SUM, comm_world, ierr_mpi)
        CALL mpialloc(dom_cluster,MAXVAL(dom_sizes),world_size,shar_rank,0,comm_shar,win_dom_cluster)
        IF (shar_rank.EQ.master) THEN 
          dom_cluster = 0
        END IF
        CALL MPI_BARRIER(comm_shar, ierr_mpi)
        dom_cluster(1:ntet_proc, wr_dex) = dom_proc(1:ntet_proc)
        IF (shar_rank.EQ.master) THEN
          CALL MPI_ALLREDUCE( MPI_IN_PLACE, dom_cluster, MAXVAL(dom_sizes)*world_size, MPI_INTEGER, MPI_SUM, comm_master, ierr_mpi )
        END IF

        ! Cluster position and diameter
        r_cluster(:,wr_dex) = SUM(tet_cen(:,dom_proc(1:ntet_proc)),DIM=2)/ntet_proc 
        d_cluster(wr_dex) = 2.0 * SQRT( SUM(NORM2(tet_cen(:,dom_proc(1:ntet_proc))-SPREAD(r_cluster(:,wr_dex), DIM=2, NCOPIES=ntet_proc),DIM=1)**2) / ntet_proc)
        CALL MPI_BARRIER(comm_shar, ierr_mpi)
        IF (shar_rank.EQ.master) THEN
          CALL MPI_ALLREDUCE( MPI_IN_PLACE, r_cluster, 3*world_size, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
          CALL MPI_ALLREDUCE( MPI_IN_PLACE, d_cluster,   world_size, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
        END IF

        ! Quadrupole matrix
        CALL mumaterial_calcquad()

        ! Determine which clusters are in the "mid-field" from our cluster
        ALLOCATE(lisfar(world_size),mid_mask(world_size))
        lisfar = NORM2(r_cluster - SPREAD(r_cluster(:,wr_dex),DIM=2, NCOPIES=world_size),DIM=1) > cutoff * d_cluster
        mid_mask = (.NOT.lisfar) .AND. ([(i, i=1, world_size)] .NE. wr_dex)
        ALLOCATE(mid_ids(COUNT(mid_mask)))
        mid_ids = PACK([(i, i=1, world_size)], MASK=mid_mask)
        ntet_mid_proc = SUM(dom_sizes(mid_ids))

        ! Array of which elements are assigned to which cluster
        ALLOCATE(dom_mid_proc(ntet_mid_proc))
        idx = 0
        DO i = 1, SIZE(mid_ids)
          j = mid_ids(i)
          cluster = PACK(dom_cluster(:,j),MASK=dom_cluster(:,j)>0)
          n = SIZE(cluster)
          dom_mid_proc(idx+1:idx+n) = cluster
          idx = idx + n
        END DO
        DEALLOCATE(mid_ids,mid_mask)

        CALL MPI_BARRIER(comm_shar, ierr_mpi)
#endif
      ELSE ! Non-MPI
        ntet_mid_proc = 0
        ALLOCATE(dom_mid_proc(0))
        ALLOCATE(lisfar(1))
        lisfar(1) = .FALSE.
      END IF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Division info print (only for MPI)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(MPI_OPT)
      IF (lcomm) THEN
        CALL MPI_ALLREDUCE(ntet_proc,ntet_proc_min,1,MPI_INTEGER,MPI_MIN,comm_world,ierr_mpi)
        CALL MPI_ALLREDUCE(ntet_proc,ntet_proc_max,1,MPI_INTEGER,MPI_MAX,comm_world,ierr_mpi)
        CALL MPI_ALLREDUCE(ntet_shar,ntet_shar_min,1,MPI_INTEGER,MPI_MIN,comm_world,ierr_mpi)
        CALL MPI_ALLREDUCE(ntet_shar,ntet_shar_max,1,MPI_INTEGER,MPI_MAX,comm_world,ierr_mpi)        
        CALL MPI_ALLREDUCE(ntet_mid_proc,ntet_mid_proc_min,1,MPI_INTEGER,MPI_MAX,comm_world,ierr_mpi)        
        CALL MPI_ALLREDUCE(ntet_mid_proc,ntet_mid_proc_max,1,MPI_INTEGER,MPI_MAX,comm_world,ierr_mpi)        
        IF (lverb) THEN 
          d_cluster_min = MINVAL(d_cluster)
          d_cluster_max = MAXVAL(d_cluster)
          ! Find worst element pair
          idx = MAXLOC(d_cluster,DIM=1)
          n = dom_sizes(idx)
          tile1 = -1
          tile2 = -1
          d_worst = -1.0
          ALLOCATE(d(n))
          DO k = 1, n
            i_tile = dom_cluster(k,idx)
            d = NORM2(tet_cen(:,dom_cluster(1:n,idx))-SPREAD(SOURCE=tet_cen(:,i_tile), DIM=2, NCOPIES=n),DIM=1)
            d_max = MAXVAL(d)
            IF (d_max.GT.d_worst) THEN
              d_worst = d_max
              tile1 = i_tile
              tile2 = dom_cluster(MAXLOC(d,DIM=1),idx)
            END IF
          END DO
          DEALLOCATE(d)
          WRITE(6,*)               ' ------- Domain Division ------'
          WRITE(6,'(3X,A,I7)')        'MPI Nodes    : ',master_size
          WRITE(6,'(3X,A,I0,A,I0,A)') 'Node range   : [',ntet_shar_min,', ',ntet_shar_max,']'
          WRITE(6,'(3X,A,I7)')        'MPI Threads  : ',world_size
          WRITE(6,'(3X,A,I0,A,I0,A)') 'Thread range : [',ntet_proc_min,', ',ntet_proc_max,']'
          WRITE(6,'(3X,A,ES0.3,A,ES0.3,A)') 'Cluster diam.: [',d_cluster_min,', ',d_cluster_max,'] m'
          WRITE(6,'(3X,A,I0,A,I0,A,I0,A,ES0.3,A)') 'Worst pair in rank ', idx, ': [',tile1,', ',tile2,'] (',d_worst,' m)'
          WRITE(6,'(3X,A,I0,A,I0,A)') 'Dipole range : [',ntet_mid_proc_min,', ',ntet_mid_proc_max,']'
          FLUSH(6)         
        END IF
      END IF
      DEALLOCATE(dom_sizes)
#endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Determine nearest neighbors (array includes self)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL mumaterial_getneighbours()  
      nbrs_proc_min = MINVAL(nbrs_count)
      nbrs_proc_max = MAXVAL(nbrs_count)
#if defined(MPI_OPT)
      IF (lcomm) THEN
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,nbrs_proc_min,1,MPI_INTEGER,MPI_MIN,comm_world,ierr_mpi)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,nbrs_proc_max,1,MPI_INTEGER,MPI_MAX,comm_world,ierr_mpi)
      END IF
#endif
      IF (lverb) THEN
        WRITE(6,'(3X,A,I0,A,I0,A)') 'Neighbor range  : [',nbrs_proc_min,', ',nbrs_proc_max,']'
        FLUSH(6)
      END IF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate H_app
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE(H_app(3,ntet_proc))
      H_app(:,:) = 0.0
      DO i = 1, ntet_proc
        i_tile = dom_proc(i)
        CALL getBfld(tet_cen(1,i_tile), tet_cen(2,i_tile), tet_cen(3,i_tile), Bx, By, Bz)
        H_app(:,i) = [Bx/MU0, By/MU0, Bz/MU0]
      END DO
      H_app_norm_min = MINVAL(NORM2(H_app,DIM=1))
      H_app_norm_max = MAXVAL(NORM2(H_app,DIM=1))
#if defined(MPI_OPT)
      IF (lcomm) THEN
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,H_app_norm_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm_world,ierr_mpi)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,H_app_norm_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm_world,ierr_mpi)
      END IF
#endif
      IF (lverb) THEN
        WRITE(6,'(3X,A,ES11.2,A,ES11.2,A)') 'Happ-field range: [',H_app_norm_min,', ',H_app_norm_max,'] A/m'
        FLUSH(6)
      END IF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate N_store
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NULLIFY(N_store)
      ALLOCATE(N_store(3,3,nbrs_maxc,ntet_proc))
      N_store(:,:,:,:) = 0.0
      DO i = 1, ntet_proc
        i_tile = dom_proc(i)
        DO j = 1, nbrs_count(i)
          j_tile = nbrs(j,i)
          CALL mumaterial_getN(vertex(:,tet(1,j_tile)), vertex(:,tet(2,j_tile)), vertex(:,tet(3,j_tile)), vertex(:,tet(4,j_tile)), tet_cen(:,i_tile), N_store(:,:,j,i)) 
        END DO 
      END DO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate inv_mat_local for constant-mu-cases
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE(inv_mat_local(3,3,ntet_proc))
      inv_mat_local = 0.0
      DO i = 1, ntet_proc
        i_tile = dom_proc(i)
        stype = state_type(state_dex(i_tile))
        IF (stype.EQ.3) THEN
            DO j = 1, nbrs_count(i)
              IF (nbrs(j,i).EQ.i_tile) EXIT ! i_tile is always in nbrs
            END DO 
            MAT = -(constant_mu(state_dex(i_tile))- 1)*N_store(:,:,j,i)
            MAT(1,1) = MAT(1,1) + 1.0
            MAT(2,2) = MAT(2,2) + 1.0
            MAT(3,3) = MAT(3,3) + 1.0
            CALL mumaterial_inv33(MAT)
            inv_mat_local(:,:,i) = MAT
        END IF
      END DO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Begin iterations
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_BARRIER(comm_world, ierr_mpi)
#endif
      IF (lverb) WRITE (6,*) ' ------- MUMAT init done -------'
      CALL mumaterial_iterate_M()
      IF (lverb) WRITE (6,*) ' ---- MUMAT iterations done ----'

      ! DEALLOCATE Helpers
      DEALLOCATE(nbrs, nbrs_count)
      DEALLOCATE(N_store,inv_mat_local)
      DEALLOCATE(H_app)
      DEALLOCATE(lisfar)

      RETURN
      END SUBROUTINE mumaterial_init_new

!-----------------------------------------------------------------------
! mumaterial_iterate_M: Iteration loop
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_iterate_M()

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Mnorm
      CHARACTER(LEN=6) :: str
      !------------------------ PRIMARY PICARD LOOP --------------------------!
      INTEGER :: iter_n, i, i_tile, j, j_tile
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: res_M, res_M_prev
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lambda_n
      !----------------------- SECONDARY PICARD LOOP -------------------------!
      INTEGER :: stype
      INTEGER :: iter_2, maxiter_2
      DOUBLE PRECISION :: lambda_k
      DOUBLE PRECISION :: H_i(3), N_self(3,3)
      DOUBLE PRECISION :: H_old(3), H_targ(3), H_new(3), H_norm, res_H(3)
      DOUBLE PRECISION :: M_targ_norm, M_targ(3), M_old(3)
      DOUBLE PRECISION :: res_M_2(3), relM, relH, threshold_2
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: H_prev
      DOUBLE PRECISION :: dH_rel_max, dH(3), dH_norm, H_norm_prev
      LOGICAL :: lgoodsec
      ! Hard magnet:
      DOUBLE PRECISION :: M_rem_norm
      DOUBLE PRECISION :: u_ea(3), u_oa_1(3), u_oa_2(3), mu_ea, mu_oa
      !------------------ BACKGROUND FIELD CALCULATION -----------------------!
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE  :: H_ext
      LOGICAL, DIMENSION(:), ALLOCATABLE             :: is_midfield
      INTEGER, DIMENSION(:), ALLOCATABLE             :: dom_mid_nn
      INTEGER                                        :: n_mid_nn
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: r_vec, r_hat
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: H_dip, mom_nn, H_quad,t
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: q
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: r_norm, mrdotrhat
      DOUBLE PRECISION ::  Bx, By, Bz
      !--------------------------- CONVERGENCE  ------------------------------!
      DOUBLE PRECISION :: converged_proc, converged_global, converged_print
      !--------------------------- DISPLAY ONLY ------------------------------!
      INTEGER ::          rank_bad, i_bad, i_tile_bad
      DOUBLE PRECISION :: pair_in(2), pair_out(2)
      DOUBLE PRECISION :: lambda_bad
      DOUBLE PRECISION :: res_rel, res_rel_loc, res_rel_bad, res_global
      DOUBLE PRECISION :: M_targ_loc, M_targ_bad, H_new_loc, H_bad, M_new_bad
      INTEGER :: mstat(MPI_STATUS_SIZE)
      !-----------------------------------------------------------------------!

      ! Allocate helpers
      ALLOCATE(Mnorm(ntet_proc))
      ALLOCATE(res_M(3,ntet_proc),res_M_prev(3,ntet_proc))
      ALLOCATE(H_prev(3,ntet_proc))
      ALLOCATE(lambda_n(ntet_proc))
      ALLOCATE(H_ext(3,ntet_proc))
      maxiter_2 = 5000
      threshold_2 = 1E-6
      Mnorm = 1.0E-5
      lambda_n = lambdaStart
      res_M = 0.0
      H_prev = 0.0
      H_ext = H_app
      iter_n = 0
      dH_rel_max = 0.1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------ PRIMARY LOOP ---------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO
        iter_n = iter_n + 1        
        converged_proc = 0.0


        ! Reset residues
        res_rel_loc = -1.0 ! Always overwritten since |res|/|M| > 0
        res_M_prev = res_M
        res_M = 0.0
        res_global = 0.0

        !-------------------- START OF ELEMENT LOOP ----------------------!
        DO i = 1, ntet_proc
          ! Get total field without contribution from self
          i_tile = dom_proc(i)  ! Tile index in global array
          H_i = H_ext(:,i) ! Non-neighbors and external sources
          DO j = 1, nbrs_count(i)   ! Get full N.M field from neighbors
            j_tile = nbrs(j,i)
            IF (j_tile.EQ.i_tile) THEN ! Found N_self
              N_self = N_store(:,:,j,i)
              CYCLE
            END IF
            H_i = H_i + MATMUL(N_store(:,:,j,i), M(:,j_tile))
          END DO

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !------------------------- SECONDARY LOOP --------------------------!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
          H_new = H_i
          iter_2 = 0
          stype = state_type(state_dex(i_tile))
          M_targ = 0.0
          lgoodsec = .TRUE.
          ! Select type (hard magnet, soft magnet, linear medium)
          SELECT CASE (stype)
            !-----------------------------------------------------------------!   
            CASE (1) ! Hard magnet
              M_rem_norm = NORM2(Mrem(:,state_dex(i_tile)))
              mu_ea = constant_mu(  state_dex(i_tile))
              mu_oa = constant_mu_o(state_dex(i_tile))

              ! Get easy axis (ea) and off axis (oa); ea assumed parallel to remanent magnetization
              u_ea = Mrem(:,state_dex(i_tile))/M_rem_norm 
              IF (u_ea(2).NE.0 .OR. u_ea(3).NE.0) THEN      ! x-product of u_ea with [1, 0, 0]; x-product of u_ea with u_oa_1
                  u_oa_1 = [0.d0, u_ea(3), -u_ea(2)]
                  u_oa_2 = [-u_ea(2)*u_ea(2) - u_ea(3)*u_ea(3), u_ea(1)*u_ea(2), u_ea(1)*u_ea(3)]
              ELSE                                          ! x-product of u_ea with [0, 1, 0]; x-product of u_ea with u_oa_1
                  u_oa_1 = [-u_ea(3), 0.d0, u_ea(1)]
                  u_oa_2 = [u_ea(1)*u_ea(2), -u_ea(1)*u_ea(1) - u_ea(3)*u_ea(3), u_ea(2)*u_ea(3)]
              END IF
              u_oa_1 = u_oa_1/NORM2(u_oa_1)
              u_oa_2 = u_oa_2/NORM2(u_oa_2)
              
              ! Picard factor
              lambda_k = MIN(1/mu_ea, 1/mu_oa, 0.5)
              DO
                iter_2 = iter_2 + 1
                H_old = H_new
                M_old = M_targ
                ! Determine magnetization taking into account easy axis
                M_targ = (M_rem_norm + (mu_ea-1)*DOT_PRODUCT(H_new,u_ea))*u_ea &
                                     + (mu_oa-1)*DOT_PRODUCT(H_new,u_oa_1)*u_oa_1 &
                                     + (mu_oa-1)*DOT_PRODUCT(H_new,u_oa_2)*u_oa_2
                ! Update H-field
                H_targ = H_i + MATMUL(N_self, M_targ)
                res_H = H_targ - H_old
                H_new = H_old + lambda_k * res_H

                res_M_2 = M_targ - M_old
                relM = NORM2(res_M_2)/MAX(NORM2(M_targ), small)
                relH = NORM2(res_H  )/MAX(NORM2(H_targ), small)

                ! Exit loop if converged or max iter exceeded
                IF (((relH.LE.threshold_2).AND.(relM.LE.threshold_2)) & 
                    .OR.(iter_2.GE.maxiter_2)) THEN
                  ! Cap change in H
                  H_norm_prev = NORM2(H_prev(:,i))
                  dH = H_new-H_prev(:,i)
                  dH_norm = NORM2(dH)
                  IF (H_norm_prev.GE.small) THEN
                    IF (dH_norm/H_norm_prev.GT.dH_rel_max) THEN
                      H_new = H_prev(:,i) + (dH_rel_max*H_norm_prev)*(dH/dH_norm)
                      lgoodsec = .FALSE.
                    END IF
                  END IF
                  H_norm = NORM2(H_new)
                  M_targ = (M_rem_norm + (mu_ea-1)*DOT_PRODUCT(H_new,u_ea))*u_ea &
                                       + (mu_oa-1)*DOT_PRODUCT(H_new,u_oa_1)*u_oa_1 &
                                       + (mu_oa-1)*DOT_PRODUCT(H_new,u_oa_2)*u_oa_2
                  IF (iter_2.GT.maxiter_2)  WRITE(6,*) "  Exceeded maxiter_2 on tile ", i_tile                  
                  EXIT
                END IF
              END DO
              H_prev(:,i) = H_new
            !-----------------------------------------------------------------!   
            CASE (2) ! Soft magnet using state function
              DO
                iter_2 = iter_2 + 1
                H_old = H_new
                M_old = M_targ
                H_norm = NORM2(H_new)
                CALL mumaterial_getState(stateFunction(state_dex(i_tile))%H, stateFunction(state_dex(i_tile))%M, H_norm, M_targ_norm)
                IF (H_norm .GE. small) THEN
                  M_targ = M_targ_norm * H_new / H_norm
                  lambda_k = MIN(H_norm/M_targ_norm, 0.5)
                ELSE
                  M_targ = 0
                  lambda_k = 0.5
                END IF

                ! Update H-field
                H_targ = H_i + MATMUL(N_self, M_targ)
                res_H = H_targ - H_old
                H_new = H_old + lambda_k * res_H
                
                res_M_2 = M_targ - M_old
                relM = NORM2(res_M_2)/MAX(NORM2(M_targ), small)
                relH = NORM2(res_H  )/MAX(NORM2(H_targ), small)
                
                ! Exit loop if converged or max iter exceeded
                IF (((relH.LE.threshold_2).AND.(relM.LE.threshold_2)) & 
                    .OR.(iter_2.GE.maxiter_2)) THEN
                  ! Cap change in H
                  H_norm_prev = NORM2(H_prev(:,i))
                  dH = H_new-H_prev(:,i)
                  dH_norm = NORM2(dH)
                  IF (H_norm_prev.GT.small) THEN
                    IF (dH_norm/H_norm_prev.GT.dH_rel_max.AND.dH_norm.GT.small) THEN
                      H_new = H_prev(:,i) + (dH_rel_max*H_norm_prev)*(dH/dH_norm)
                      lgoodsec = .FALSE.
                    END IF
                  END IF
                  H_norm = NORM2(H_new)
                  ! Recalculate M
                  CALL mumaterial_getState(stateFunction(state_dex(i_tile))%H, stateFunction(state_dex(i_tile))%M, H_norm, M_targ_norm)
                  IF (H_norm .GT. small) THEN
                        M_targ = M_targ_norm * H_new / H_norm
                  ELSE
                        M_targ = 0
                  END IF
                  IF (iter_2.GT.maxiter_2)  WRITE(6,*) "  Exceeded maxiter_2 on tile ", i_tile                  
                  EXIT
                END IF
              END DO
              H_prev(:,i) = H_new
            !-----------------------------------------------------------------!
            CASE (3) ! Constant permeability, solve directly using inverse: 
              ! MAT = (I-(mu_r-1)*N)
              ! H = inv(MAT)*Hext
              mu_ea = constant_mu(state_dex(i_tile))
              M_targ = (mu_ea - 1) * MATMUL(inv_mat_local(:,:,i),H_new)
              H_norm = NORM2(M_targ/(mu_ea - 1))
            !-----------------------------------------------------------------!
            CASE DEFAULT ! Something went wrong.
              WRITE(6,*) "  Unknown magnet type: ", stype
              STOP
          END SELECT
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !----------------------- END SECONDARY LOOP ------------------------!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

          ! Update magnetization
          res_M(:,i) = M_targ - M(:,i_tile) ! target - old
          M(:,i_tile) = M(:,i_tile) + lambda_n(i)*res_M(:,i)

          ! Check convergence progress
          res_rel = NORM2(res_M(:,i))/NORM2(M_targ)
          res_global = res_global + res_rel*tet_vol(i_tile)
	        IF ((res_rel.LT.threshold) .AND. (lgoodsec)) THEN                     
            converged_proc = converged_proc + tet_vol(i_tile)
          END IF

          ! Find worst residual for printing
          IF (res_rel.GT.res_rel_loc) THEN
            res_rel_loc = res_rel
            M_targ_loc = NORM2(M_targ)
            H_new_loc = H_norm
            i_bad = i
          END IF

          ! Dampen evolution during oscillations in res_M
          IF (res_rel.LT.threshold) THEN
            ! Change nothing
          ELSE IF ((DOT_PRODUCT(res_M(:,i),res_M_prev(:,i))<0.0) .AND. (NORM2(res_M(:,i))>NORM2(res_M_prev(:,i)))) THEN 
            lambda_n(i) = lambda_n(i) * lambdaFactor
            lambda_n(i) = MAX(lambda_n(i),0.001)
          ELSE IF (NORM2(res_M(:,i)).LT.NORM2(res_M_prev(:,i))) THEN
            lambda_n(i) = lambda_n(i) * 1.05
            lambda_n(i) = MIN(lambda_n(i), 0.75)
          END IF
        END DO
        !------------------------ END OF ELEMENT LOOP ------------------------!
        !---------------------------------------------------------------------!
        !------------------------- PRINT TO SCREEN ---------------------------!
        ! MPI Communication
        IF (lcomm) THEN
#if defined(MPI_OPT)
          ! Find worst element across all ranks using MPI_ALLREDUCE
          pair_in(1) = res_rel_loc
          pair_in(2) = REAL(world_rank)
          CALL MPI_ALLREDUCE(pair_in,pair_out,1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, comm_world, ierr_mpi)
          res_rel_bad = pair_out(1)
          rank_bad = INT(pair_out(2))
          ! Communicate info to master (only proc that prints)
          IF (world_rank.EQ.rank_bad) THEN ! If this proc has the bad element:
            IF (lismaster) THEN            ! If this proc is also the master, just grab info
              i_tile_bad = dom_proc(i_bad)
              lambda_bad = lambda_n(i_bad) ! I know lambda_n has already been updated at this point.
              M_targ_bad = M_targ_loc
              H_bad = H_new_loc
            ELSE                           ! If this proc is NOT the master, send info to master
              CALL MPI_SEND(dom_proc(i_bad),    1, MPI_INTEGER,          0, 1240, comm_world, ierr_mpi)
              CALL MPI_SEND(lambda_n(i_bad), 1, MPI_DOUBLE_PRECISION, 0, 1241, comm_world, ierr_mpi) 
              CALL MPI_SEND(M_targ_loc,      1, MPI_DOUBLE_PRECISION, 0, 1242, comm_world, ierr_mpi) 
              CALL MPI_SEND(H_new_loc,       1, MPI_DOUBLE_PRECISION, 0, 1243, comm_world, ierr_mpi) 
            END IF
          ELSE IF (lismaster) THEN         ! If this proc does not have the bad element, AND is the master, receive info:
            CALL MPI_RECV(i_tile_bad, 1, MPI_INTEGER,          rank_bad, 1240, comm_world, mstat, ierr_mpi)
            CALL MPI_RECV(lambda_bad, 1, MPI_DOUBLE_PRECISION, rank_bad, 1241, comm_world, mstat, ierr_mpi) 
            CALL MPI_RECV(M_targ_bad, 1, MPI_DOUBLE_PRECISION, rank_bad, 1242, comm_world, mstat, ierr_mpi)
            CALL MPI_RECV(H_bad,      1, MPI_DOUBLE_PRECISION, rank_bad, 1243, comm_world, mstat, ierr_mpi)
          END IF
          CALL MPI_ALLREDUCE(converged_proc, converged_global,  1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_world, ierr_mpi) 
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,         res_global,  1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_world, ierr_mpi) 
          CALL MPI_BARRIER(comm_world, ierr_mpi)
#endif
        ELSE
            i_tile_bad = dom_proc(i_bad)
            lambda_bad = lambda_n(i_bad)
            M_targ_bad = M_targ_loc
            H_bad = H_new_loc            
            converged_global = converged_proc
        END IF
        res_global = res_global/SUM(tet_vol)
        IF (ldosync) CALL mumaterial_syncM()

        IF (lverb) THEN 
          IF (iter_n.EQ.1) THEN
            WRITE(6,*) ''
            WRITE(6,*) '  iter  %good  res(avg)  res(bad) |    tile   H(norm)   M(targ)   M(norm)   '
            WRITE(6,*) '==============================================================================='
          END IF
          M_new_bad = NORM2(M(:,i_tile_bad))
          WRITE(6,'(1X, I6, 2X, F5.1, 1X, ES9.2, 1X, ES9.2, A, &
                    I7, 1X,ES9.2,1X, ES9.2,1X, ES9.2)') & 
                  iter_n, converged_global*100.0/SUM(tet_vol),res_global,res_rel_bad, ' | ',  &
                  i_tile_bad, H_bad, M_targ_bad, M_new_bad
          CALL FLUSH(6)
        END IF

        IF (((converged_global*100.0/SUM(tet_vol) .GE.convCheck).AND.(MOD(iter_n,100).LE.20)).OR.(iter_n.GE.maxIter)) THEN
            IF (lverb) WRITE(6,*) "  MUMAT:  Stopping"
            EXIT
        END IF
        !---------------------------------------------------------------------!
        !--------------------- UPDATE BACKGROUND H_APP -----------------------!
        ALLOCATE(is_midfield(ntet_mid_proc))
        H_ext =  H_app ! Static field using slice
        CALL mumaterial_calcquad()
        DO i = 1, ntet_proc
          i_tile = dom_proc(i)
          !---------------- CONTRIBUTION FROM DISTANT CLUSTERS -----------------!
          ALLOCATE(r_vec(3,world_size),r_norm(world_size),r_hat(3,world_size),&
                   H_dip(3,world_size),mrdotrhat(world_size),H_quad(3,world_size))
          r_vec = SPREAD(tet_cen(:,i_tile),DIM=2,NCOPIES=world_size)-r_cluster
          r_norm = NORM2(r_vec, DIM=1)
          r_hat = r_vec / SPREAD(r_norm, DIM=1, NCOPIES=3)
          ! Dipole contribution: 1/(4*pi*r^3) *  (3*(m.rhat)rhat - m)
          mrdotrhat = SUM(m_cluster*r_hat,DIM=1)
          H_dip = INV4PI*(3.0*SPREAD(mrdotrhat,DIM=1,NCOPIES=3)*r_hat-m_cluster)/SPREAD(r_norm**3,DIM=1,NCOPIES=3)
          ! Quadrupole contribution: 1/(8*pi*r^4) * ((rhat.Q rhat)rhat - 2*Q.rhat)
          ALLOCATE(t(3,world_size),q(world_size)) ! Helpers for vectorization
          t(1,:) = Q_cluster(1,:)*r_hat(1,:) + Q_cluster(4,:)*r_hat(2,:) + Q_cluster(5,:)*r_hat(3,:)
          t(2,:) = Q_cluster(4,:)*r_hat(1,:) + Q_cluster(2,:)*r_hat(2,:) + Q_cluster(6,:)*r_hat(3,:)
          t(3,:) = Q_cluster(5,:)*r_hat(1,:) + Q_cluster(6,:)*r_hat(2,:) + Q_cluster(3,:)*r_hat(3,:)
          q = r_hat(1,:)*t(1,:) + r_hat(2,:)*t(2,:) + r_hat(3,:)*t(3,:)
          H_quad = INV8PI*(5*SPREAD(q, DIM=1, NCOPIES=3)*r_hat-2*t)/SPREAD(r_norm**4,DIM=1,NCOPIES=3)
          DEALLOCATE(t,q)
          ! Corrections
          WHERE (.NOT. SPREAD(lisfar, DIM=1, NCOPIES=3))
            H_dip  = 0.0
            H_quad = 0.0
          END WHERE
          H_ext(:,i) = H_ext(:,i) + SUM(H_dip,DIM=2) + SUM(H_quad,DIM=2)
          DEALLOCATE(r_vec,r_norm,r_hat,H_dip,mrdotrhat,H_quad)
          !---------------- CONTRIBUTION FROM MID-FIELD DIPOLES ----------------!
          ! Get all non-neighbors
          is_midfield = .FALSE.
          DO j = 1, ntet_mid_proc
            j_tile = dom_mid_proc(j)
            is_midfield(j) = .NOT. ANY(nbrs(1:nbrs_count(i), i) .EQ. j_tile)
          END DO
          dom_mid_nn = PACK(dom_mid_proc, MASK=is_midfield)
          n_mid_nn = SIZE(dom_mid_nn)

          IF (n_mid_nn.GT.0) THEN
            ALLOCATE(r_vec(3,n_mid_nn),r_norm(n_mid_nn),r_hat(3,n_mid_nn),&
                     H_dip(3,n_mid_nn),mrdotrhat(n_mid_nn),mom_nn(3,n_mid_nn))      
            r_vec = SPREAD(tet_cen(:,i_tile),DIM=2,NCOPIES=n_mid_nn)-tet_cen(:,dom_mid_nn)
            r_norm = NORM2(r_vec, DIM=1)
            r_hat = r_vec / SPREAD(r_norm, DIM=1, NCOPIES=3)
            mom_nn = M(:,dom_mid_nn)*SPREAD(tet_vol(dom_mid_nn), DIM=1, NCOPIES=3)
            mrdotrhat = SUM(mom_nn*r_hat,DIM=1)
            H_dip = INV4PI*(3.0*SPREAD(mrdotrhat,DIM=1,NCOPIES=3)*r_hat-mom_nn)/SPREAD(r_norm**3,DIM=1,NCOPIES=3)
            H_ext(:,i) = H_ext(:,i) + SUM(H_dip,DIM=2)              
            DEALLOCATE(r_vec,r_norm,r_hat,H_dip,mrdotrhat,mom_nn)
          END IF
          DEALLOCATE(dom_mid_nn)

        END DO  
        DEALLOCATE(is_midfield)
        !---------------------------------------------------------------------!
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------- END PRIMARY LOOP -------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DEALLOCATE(Mnorm,res_M,res_M_prev,H_prev,lambda_n,H_ext)
      RETURN
      END SUBROUTINE mumaterial_iterate_M

      SUBROUTINE mumaterial_inv33(B)
      !-----------------------------------------------------------------------
      ! mumaterial_inv33: Helper function for inverting a 3x3 matrix.
      !-----------------------------------------------------------------------
      ! param[inout]: B: target matrix 
      !-----------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(inout) :: B(3,3)
      DOUBLE PRECISION :: det, INV(3,3)

      DET = B(1,1)*(B(2,2)*B(3,3)-B(2,3)*B(3,2)) + &
            B(1,2)*(B(2,3)*B(3,1)-B(2,1)*B(3,3)) + &
            B(1,3)*(B(2,1)*B(3,2)-B(2,2)*B(3,1))

      INV(1,1) = (B(2,2)*B(3,3)-B(2,3)*B(3,2))/DET
      INV(2,1) = (B(2,3)*B(3,1)-B(2,1)*B(3,3))/DET
      INV(3,1) = (B(2,1)*B(3,2)-B(2,2)*B(3,1))/DET

      INV(1,2) = (B(1,3)*B(3,2)-B(1,2)*B(3,3))/DET
      INV(2,2) = (B(1,1)*B(3,3)-B(1,3)*B(3,1))/DET
      INV(3,2) = (B(1,2)*B(3,1)-B(1,1)*B(3,2))/DET

      INV(1,3) = (B(1,2)*B(2,3)-B(1,3)*B(2,2))/DET
      INV(2,3) = (B(1,3)*B(2,1)-B(1,1)*B(2,3))/DET
      INV(3,3) = (B(1,1)*B(2,2)-B(1,2)*B(2,1))/DET

      B = INV

      END SUBROUTINE mumaterial_inv33

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
        IF (DOT_PRODUCT(mumaterial_cross(v(:,1) - v(:,3), v(:,2) - v(:,3)), v(:,4) - v(:,2)) .gt. 0) THEN 
        ! normal vector of triangle is pointing towards v4, so v1 and v3 need to be interchanged
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
          IF (ABS(r(j)) .lt. 1.0D-6) THEN ! make sure position is not too close to x, y or z = 0
            r(j) = SIGN(1.0D-6, r(j))
          END IF
          IF (ABS(v(1,j)) .lt. 1.0D-6) THEN 
            v(1,j) = SIGN(1.0D-6, v(1,j))
          END IF                  
        END DO

        N_loc = 0.d0
        N_loc(1,3) = mumaterial_getNxz(r, v(1,1), v(2,2)) - mumaterial_getNxz(r, v(1,3), v(2,2))
        N_loc(2,3) = mumaterial_getNyz(r, v(1,1), v(2,2)) - mumaterial_getNyz(r, v(1,3), v(2,2))
        N_loc(3,3) = mumaterial_getNzz(r, v(1,1), v(2,2)) - mumaterial_getNzz(r, v(1,3), v(2,2))
        IF ((ISNAN(N_loc(1,3)).or.ISNAN(N_loc(2,3))).or.ISNAN(N_loc(3,3))) THEN 
              WRITE(6,*) "FOUND A NAN IN N_LOC"
              WRITE(6,*) "POS=",pos(1),pos(2),pos(3)
              WRITE(6,*) "R=",r(1),r(2),r(3)
              WRITE(6,*) "l_1=",v(1,1), "l_2=", v(1,3)
              WRITE(6,*) "h=",v(2,2)
              WRITE(6,*)
        END IF
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
      !-----------------------------------------------------------------------
      ! mumaterial_cross: calculates the cross product a * b.
      !-----------------------------------------------------------------------
      ! param[in]: a. 1x3 vector
      ! param[in]: b: 1x3 vector
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN), DIMENSION(3) :: a, b
      DOUBLE PRECISION, DIMENSION(3) :: mumaterial_cross

      mumaterial_cross(1) = a(2)*b(3) - a(3)*b(2)
      mumaterial_cross(2) = a(3)*b(1) - a(1)*b(3)
      mumaterial_cross(3) = a(1)*b(2) - a(2)*b(1)

      RETURN
      END FUNCTION mumaterial_cross

      FUNCTION mumaterial_gettetvolume(v1,v2,v3,v4)
      !-----------------------------------------------------------------------
      ! mumaterial_gettetvolume: Calculates the volume of an element.
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3), INTENT(in) :: v1, v2, v3, v4
      DOUBLE PRECISION :: mumaterial_gettetvolume

      mumaterial_gettetvolume = ABS(dot_product(v1-v4,mumaterial_cross(v2-v4,v3-v4)))/6.0
      RETURN

      END FUNCTION mumaterial_gettetvolume

      SUBROUTINE mumaterial_getneighbours()
      !-----------------------------------------------------------------------
      ! mumaterial_getneighbours: Finds all UK neighbors of an element.
      !-----------------------------------------------------------------------
      INTEGER :: i, j, k, i_tile
      DOUBLE PRECISION, ALLOCATABLE ::  dist(:)
      LOGICAL, ALLOCATABLE :: mask(:)

      ALLOCATE(nbrs_count(ntet_proc),mask(ntet),dist(ntet))
      nbrs_count = 0
      
      ! Get largest neighbor count for allocation first
      DO i = 1, ntet_proc
        i_tile = dom_proc(i)
        dist = NORM2(tet_cen - SPREAD(SOURCE=tet_cen(:,i_tile), DIM=2, NCOPIES=ntet),DIM=1)
        nbrs_count(i) = COUNT(dist.LE.padFactor*tet_rad)
      END DO

      nbrs_maxc = MAXVAL(nbrs_count)
      ALLOCATE(nbrs(nbrs_maxc,ntet_proc))

      ! Actual neighbour loop
      DO i = 1, ntet_proc
        i_tile = dom_proc(i)
        dist = NORM2(tet_cen - SPREAD(SOURCE=tet_cen(:,i_tile), DIM=2, NCOPIES=ntet),DIM=1)
        mask = dist.LE.padFactor*tet_rad
        j = 0
        DO k = 1, ntet
          IF (mask(k)) THEN
            j = j + 1
            nbrs(j,i) = k
          END IF
        END DO
      END DO
      DEALLOCATE(mask,dist)

      END SUBROUTINE mumaterial_getneighbours

      SUBROUTINE mumaterial_split(boxin,targ,box1,box2)
      !-----------------------------------------------------------------------
      ! mumaterial_split: Binary splits a collection of elements (boxin) into
      ! two spatially localized subdomains (box1, box2) based on the position
      ! of the elements and the target split (targ). box1 has size
      ! nint(targ*size(boxin)), box2 has remainder.
      !-----------------------------------------------------------------------
      ! param[in]: boxin. indices of elements
      ! param[in]: targ. relative size of box1 out compared to boxin.
      ! param[out]: box1. first output subdomain
      ! param[out]: box2. second output subdomain.
      !-----------------------------------------------------------------------
        
      USE qsort ! quicksort

      IMPLICIT NONE

      INTEGER, INTENT(in)               :: boxin(:)
      DOUBLE PRECISION, INTENT(in)      :: targ
      INTEGER, ALLOCATABLE, INTENT(out) :: box1(:), box2(:) 
      INTEGER, DIMENSION(:), ALLOCATABLE :: idx, temp1, temp2
      INTEGER :: i, i_dim, size1, size2, boxsize
      DOUBLE PRECISION :: r_com(3)
      DOUBLE PRECISION :: d1, d2, d, d_best

      ! Check if box contains enough elements; otherwise stop
      boxsize = SIZE(boxin)
      size1 = NINT(targ*boxsize)
      size2 = boxsize-size1
      IF ((size1.LT.1).OR.(size2.LT.1)) THEN
        WRITE(6,"(A,I0,A,F0.3,A,I0,A)") "  MUMAT_SPLIT: RANK ", world_rank, & 
          " CANNOT SPLIT (targ=", targ, ", boxsize=", boxsize, ")"
        WRITE(6,"(A)") "  MUMAT_SPLIT: FORCE STOPPING CALCULATIONS"
#if defined(MPI_OPT)
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr_mpi)
#else
        STOP
#endif
      END IF
      !---------------------------------------
      ALLOCATE(idx(boxsize),temp1(size1),temp2(size2),box1(size1),box2(size2))
      d_best = -1.0 ! Overwritten anyway
      !---------------------------------------
      ! Loop over each direction
      DO i_dim = 1, 3
        DO i = 1, boxsize
          idx(i) = i
        END DO
        CALL quicksort(tet_cen(i_dim,boxin),idx,1,boxsize) ! idx is sorted based on coordinates
        temp1 = boxin(idx(1:size1))
        temp2 = boxin(idx(size1+1:boxsize))
      !---------------------------------------
        ! Evaluate distances from center of masses
        r_com = SUM(tet_cen(:,temp1),DIM=2) / size1
        d1 = SUM( &
              NORM2(tet_cen(:,temp1)-SPREAD(r_com,DIM=2,NCOPIES=size1),DIM=1) &
                ) / size1
        r_com = SUM(tet_cen(:,temp2),DIM=2) / size2
        d2 = SUM( & 
              NORM2(tet_cen(:,temp2)-SPREAD(r_com,DIM=2,NCOPIES=size2),DIM=1) & 
                ) / size2
        d = (d1+d2)/2
      !---------------------------------------
        ! Update if dimension is better
        IF ((d.LT.d_best).OR.(i_dim.EQ.1)) THEN
          box1 = temp1
          box2 = temp2
          d_best = d
        END IF
      END DO
      !---------------------------------------
      DEALLOCATE(temp1, temp2)

      END SUBROUTINE mumaterial_split

      SUBROUTINE mumaterial_calcquad
      !-----------------------------------------------------------------------
      ! mumaterial_calcquad: Recalculates Q tensor for quadrupole contribution
      !-----------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R, mom
      INTEGER :: wr_dex
      DOUBLE PRECISION :: A, B, C

      IF (shar_rank.EQ.master) THEN
        Q_cluster = 0.0
      END IF
      CALL MPI_BARRIER(comm_shar, ierr_mpi)
      
      wr_dex = world_rank+1
      ALLOCATE(R(3,ntet_proc),mom(3,ntet_proc))
      R = tet_cen(:,dom_proc(1:ntet_proc))-SPREAD(r_cluster(:,wr_dex),DIM=2,NCOPIES=ntet_proc)
      mom = M(:,dom_proc(1:ntet_proc))*SPREAD(tet_vol(dom_proc(1:ntet_proc)),DIM=1,NCOPIES=3)

      A = DOT_PRODUCT(R(1,:),mom(1,:))
      B = DOT_PRODUCT(R(2,:),mom(2,:))
      C = DOT_PRODUCT(R(3,:),mom(3,:))

      Q_cluster(1,wr_dex) = 2.0/3.0*(2*A-B-C) ! xx
      Q_cluster(2,wr_dex) = 2.0/3.0*(2*B-A-C) ! yy
      Q_cluster(3,wr_dex) = -(Q_cluster(1,wr_dex)+Q_cluster(2,wr_dex)) ! zz

      Q_cluster(4,wr_dex) = DOT_PRODUCT(R(1,:),mom(2,:))+DOT_PRODUCT(R(2,:),mom(1,:)) ! xy
      Q_cluster(5,wr_dex) = DOT_PRODUCT(R(1,:),mom(3,:))+DOT_PRODUCT(R(3,:),mom(1,:)) ! xz
      Q_cluster(6,wr_dex) = DOT_PRODUCT(R(2,:),mom(3,:))+DOT_PRODUCT(R(3,:),mom(2,:)) ! yz
      DEALLOCATE(R,mom)

      IF (shar_rank.EQ.master) THEN
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, Q_cluster, 6*world_size, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi)
      END IF
      CALL MPI_BARRIER(comm_shar, ierr_mpi)

      END SUBROUTINE mumaterial_calcquad

      SUBROUTINE mumaterial_syncM()
      !-----------------------------------------------------------------------
      ! mumaterial_syncM: Synchronizes M across all MPI nodes and recalculates
      ! magnetic moment of clusters. Called after each iteration.
      !-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: i, i_tile
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: M_local

      ! First recalculate cluster moments
      IF (shar_rank.EQ.master) THEN 
        m_cluster = 0.0 ! master zeroes shared memory window
      END IF
      CALL MPI_BARRIER(comm_shar, ierr_mpi) ! threads wait for zeroing
      m_cluster(:,world_rank+1) = SUM(M(:, dom_proc(1:ntet_proc)) * &
                                    SPREAD(tet_vol(dom_proc(1:ntet_proc)), DIM=1, NCOPIES=3), DIM=2) ! everyone does work
      IF (shar_rank.EQ.master) THEN ! synchronize on all mpi nodes
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, m_cluster, 3*world_size, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
      END IF

      ! Global M array
      CALL MPI_BARRIER( comm_shar, ierr_mpi)
      IF (shar_rank.EQ.master) THEN
        ALLOCATE(M_local(3,ntet))
        M_local = 0.0
        M_local(:, dom_shar(1:ntet_shar)) = M(:, dom_shar(1:ntet_shar))
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, M_local, 3*ntet, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
        M = M_local
        DEALLOCATE(M_local)
      END IF
      CALL MPI_BARRIER( comm_shar, ierr_mpi)

      END SUBROUTINE mumaterial_syncM

      SUBROUTINE mumaterial_getState(fx, fy, xq, yq)
      !-----------------------------------------------------------------------
      ! mumaterial_getState: Interpolates a function f at xq to get a value y 
      ! using B-splines based on De Boor's algorithm
      !-----------------------------------------------------------------------
      ! param[in]:  fx. x-coordinates of function to be interpolated
      ! param[in]:  fy. y-values of function to be interpolated
      ! param[in]:  xq. evaluation point
      ! param[out]: yq. interpolated f(xq)
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: fx(:), fy(:), xq
      DOUBLE PRECISION, INTENT(OUT) :: yq
      INTEGER :: n, i, k, p, r
      DOUBLE PRECISION :: alpha
      DOUBLE PRECISION, ALLOCATABLE :: t(:), d(:)

      p = 3 ! Degree of the polynomial used
      n = SIZE(fx)

      ALLOCATE(t(2+2*(p-1)))
      ALLOCATE(d(p+1))

      ! Determine left index k
      ! Assume fx is non-decreasing
      IF (xq .lt. fx(1)) THEN
        yq = fy(1)
        RETURN
      ELSEIF (xq .gt. fx(n)) THEN
        yq = fy(n)
        RETURN
      ELSE
        DO i = 2, n 
          IF (xq .lt. fx(i)) THEN
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
          alpha = (xq - t(i)) / (t(i+1+p-r) - t(i))
          d(i) = (1 - alpha) * d(i) + alpha * d(i+1)
        END DO
      END DO

      ! Set output variable
      yq = d(3)

      RETURN
      END SUBROUTINE mumaterial_getState
    
      SUBROUTINE mumaterial_getb_scalar(x, y, z, Bx, By, Bz)
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
      DOUBLE PRECISION, INTENT(in) :: x, y, z
      DOUBLE PRECISION :: Bx_mag, By_mag, Bz_mag
      DOUBLE PRECISION, INTENT(out) :: Bx, By, Bz

      CALL mumaterial_getbmag_scalar(x, y, z, Bx_mag, By_mag, Bz_mag)
      CALL getBfld(x, y, z, Bx, By, Bz)

      Bx = Bx + Bx_mag
      By = By + By_mag
      Bz = Bz + Bz_mag

      RETURN
      END SUBROUTINE mumaterial_getb_scalar

      SUBROUTINE mumaterial_getbmag_scalar(x, y, z, Bx, By, Bz)
      !-----------------------------------------------------------------------
      ! mumaterial_getbmag: Calculates magnetic field from magnetizations
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
      DOUBLE PRECISION :: H(3), N(3,3)
      INTEGER :: i

      H = 0.d0
      DO i = 1, ntet
        CALL mumaterial_getN(vertex(:,tet(1,i)), vertex(:,tet(2,i)), vertex(:,tet(3,i)), vertex(:,tet(4,i)), [x, y, z], N)
        H = H + MATMUL(N, M(:,i))
      END DO

      Bx = H(1) * MU0
      By = H(2) * MU0
      Bz = H(3) * MU0

      RETURN
      END SUBROUTINE mumaterial_getbmag_scalar


      SUBROUTINE mumaterial_getb_vector(x, y, z, B)
      !-----------------------------------------------------------------------
      ! mumaterial_getb_vector: Calculates total magnetic field at multiple points in space
      !-----------------------------------------------------------------------
      ! param[in]: x. x-coordinates of points at which to determine the magnetic field
      ! param[in]: y. y-coordinates of points at which to determine the magnetic field
      ! param[in]: z. z-coordinates of points at which to determine the magnetic field
      ! param[out]: B.  B-field at required points [T]
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x(:), y(:), z(:)
      DOUBLE PRECISION, INTENT(out), ALLOCATABLE :: B(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: B_local(:,:)
      INTEGER :: mystart, myend
      INTEGER :: i 
      INTEGER :: npoints

      npoints = size(x)
      mystart = 1; myend = npoints

#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_CALC_MYRANGE(comm_world, 1, npoints, mystart, myend)
#endif

      ALLOCATE(B_local(3,npoints),B(3,npoints))
      B_local = 0; B = 0
      
      DO i = mystart, myend
        CALL mumaterial_getb_scalar(x(i), y(i), z(i), B_local(1,i), B_local(2,i), B_local(3,i))
      END DO
    
#if defined(MPI_OPT)
      IF (lcomm) THEN
        CALL MPI_REDUCE(B_local,B,3*npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm_shar,ierr_mpi)
        IF (shar_rank.EQ.0) CALL MPI_ALLREDUCE( MPI_IN_PLACE,B,3*npoints,MPI_DOUBLE_PRECISION,MPI_SUM,comm_master,ierr_mpi)
      END IF
#endif

      DEALLOCATE(B_local)
      
      RETURN
      END SUBROUTINE mumaterial_getb_vector

      SUBROUTINE mumaterial_readmag(filename)
      !-----------------------------------------------------------------------
      ! mumaterial_readmag: Reads magnetization .dat file
      !-----------------------------------------------------------------------   
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER :: i, istat, iunit

      IF (lismaster) THEN
        WRITE(6,'(A)')           ' -------- MUMAT magfile --------'
        WRITE(6,'(3X,A,A)')     ' File         : ',filename
        ! open file, return if fails
        iunit = 327; istat = 0
        CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
        IF (istat/= 0) THEN
          WRITE(6,*) "ISSUE READING MAG; STOPPING"
          RETURN
        END IF
        DO i = 1, ntet
          READ(iunit, *) M(1,i),M(2,i),M(3,i)
        END DO
        CLOSE(iunit)
      END IF    

#if defined(MPI_OPT)
      IF ((lcomm).AND.(shar_rank.EQ.0)) THEN
        CALL MPI_Bcast(M,3*ntet,MPI_DOUBLE_PRECISION,0,comm_master,ierr_mpi)
      END IF
#endif
      END SUBROUTINE mumaterial_readmag


      SUBROUTINE mumaterial_writemag(str)
      !-----------------------------------------------------------------------
      ! mumaterial_writemag: Outputs magnetization to .dat file
      !-----------------------------------------------------------------------
      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(in) :: str
      CHARACTER(LEN=256) :: filename
      INTEGER :: i

      IF (lismaster) THEN
        filename = './mumat_mag_'//TRIM(str)//'.dat'
        WRITE(6,"(A)") "  MUMAT: Writing magnetization to" // filename
        OPEN(13, file=filename)
        DO i = 1, ntet
          WRITE(13, "(ES15.7,ES15.7,ES15.7)") M(1,i),M(2,i),M(3,i)
        END DO
        CLOSE(13)
      END IF

      END SUBROUTINE mumaterial_writemag

      SUBROUTINE mumaterial_output(x, y, z)
      !-----------------------------------------------------------------------
      ! mumaterial_output: Outputs B-field and points to text files
      !-----------------------------------------------------------------------
      ! param[in]: x. x-cooridinates of points at which to determine the magnetic field
      ! param[in]: y. y-cooridinates of points at which to determine the magnetic field
      ! param[in]: z. z-cooridinates of points at which to determine the magnetic field
      ! param[in]: linclvac. Whether or not vacuum magnetic field should be included.
      !-----------------------------------------------------------------------
    
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x(:), y(:), z(:)
      INTEGER :: i 
      INTEGER :: npoints
      DOUBLE PRECISION, ALLOCATABLE :: B(:,:)

      IF (lismaster) THEN
        npoints = size(x)
        WRITE(6,*) "Outputting points"
        OPEN(13, file='./points.dat')
        DO i = 1, npoints
          WRITE(13, "(F15.7,A,F15.7,A,F15.7)") x(i), ',', y(i), ',', z(i)
        END DO
        CLOSE(13)
      END IF

      CALL mumaterial_getb_vector(x, y, z, B)
 
      IF (lismaster) THEN
        WRITE(6,*) "Outputting B-field"
        OPEN(14, file='./B.dat')
        DO i = 1, npoints
          WRITE(14, "(E15.7,A,E15.7,A,E15.7)") B(1,i), ',', B(2,i), ',', B(3,i)
        END DO
        CLOSE(14)
      END IF

      RETURN
      END SUBROUTINE
!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE mumaterial_mod
