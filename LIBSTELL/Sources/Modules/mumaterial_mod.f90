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
!------------------------------------------------------------------------------
!     Module Variables
!        lverb:      Controls output to screen
!   
!       MPI
!         lcomm:     .TRUE. if code is run with MPI
!         ldosync:   .TRUE. if more than one MPI node is used
!         lismaster: .TRUE. if rank of thread in world-communicator is 0
!         ldebugX:   Debug flags: world-master (m); shar-master (s); thread (t)
!
!         comm_shar:   Shared-memory communicator
!         master_comm: Communicator of threads whose rank in comm_shar is 0
!         world_comm:  Communicator with all MPI threads
!         color:       Used to create master_comm    
!         COMM_rank:   Rank of thread in communicator COMM
!         COMM_size:   Number of threads in communicator COMM
!         
!         mydom:    Collection of tetrahedrons worked on by comm_shar
!         mypdom:   mydom with extra tetrahedrons included (see padFactor)
!         outmydom: Collection of tetrahedrons NOT worked on by comm_shar
!         domsize:  Size of mydom
!         pdomsize: Size of mypdom
!         win_OBJ:  MPI shared memory window for OBJ
!
!       Neighbours
!         Nb:         Array of neighbours for each tetrahedron (:,:)
!         NbC:        Number of neighbours for each tetrahedron (:)
!         maxNbC:     Largest neighbour count in NbC
!         Nb_domidx:  Neighbours indexed by appearance in mydom (:,:)
!         NbC_dom:    Number of neighbours per tet in mydom (:)
!         maxNbC_dom: Largest neighbour count in NbC_dom
!
!       Mesh
!         ntet:     Number of tetrahedrons in mesh
!         nvertex:  Number of vertices in mesh
!         vertex:   Coordinates for vertices in mesh (3, ntet)
!         tet:      Vertex indices for each tetrahedron (4, ntet)
!         tet_cen:  Coordinates for tetrahedron centers (3, ntet)
!         tet_vol:  Volumes of tetrahedrons (ntet)
!         tet_edge: Equivalent length of edge of tetrahedrons (ntet)
!   
!       Magnetics
!         nstate:           Number of state functions
!         state_dex:        State function for each tetrahedron (ntet)
!         state_type:       Type of state function (nstate) (1-3)
!         constant_mu:      Mu for constant permeability (nstate)
!         constant_mu_o:    Mu for orthogonal axis for hard magnet (nstate)
!         Mrem:             Remanent magnetization for hard magnet (3,nstate)
!         M:                Magnetization for all tetrahedrons (3,ntet)
!         Happ:             Applied H-field at tetrahedron centres (3, ntet)
!         mu0:              Permeability of free space: 4*pi*1E-7 [H/m]
!         N_store:          Demagnetization tensor (3,3,maxNbC,:)
!
!       User settings
!         dMmax:           Threshold error for convergence
!         maxIter:         Max allowed number of iterations
!         padFactor:   Affects number of neighbours for each tetrahedron
!         lambdaStart:     Initial value of lambda for iterations
!         lambdaFactor:    Multiplication factor for lambda
!         lambdaThresh:    Multiply lambda if error grows this number of times
!------------------------------------------------------------------------------

      CHARACTER(LEN=256), PRIVATE :: machine_string
      CHARACTER(LEN=256), PRIVATE :: date

      ! mesh variables
      INTEGER, PRIVATE  ::  ntet, nvertex
      DOUBLE PRECISION, POINTER, PRIVATE :: vertex(:,:), tet_cen(:,:), & 
                                            tet_vol(:), tet_edge(:)
      INTEGER, POINTER, PRIVATE :: tet(:,:)

      ! magnetics variables
      INTEGER, POINTER, PRIVATE :: state_dex(:), state_type(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: constant_mu(:), constant_mu_o(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: M(:,:), Happ(:,:), Mrem(:,:)
      DOUBLE PRECISION, DIMENSION(:,:,:,:), POINTER, PRIVATE :: N_store
      DOUBLE PRECISION, PRIVATE :: mu0
      INTEGER, PRIVATE :: nstate
      TYPE(stateFunctionType), PRIVATE, ALLOCATABLE :: stateFunction(:)

      ! user settings variables
      DOUBLE PRECISION, PRIVATE :: dMmax, padFactor, lambdaStart, lambdaFactor
      INTEGER, PRIVATE          :: lambdaThresh, maxIter

      ! neighbour variables
      INTEGER, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: Nb, Nb_domidx
      INTEGER, DIMENSION(:),   ALLOCATABLE, PRIVATE :: NbC, NbC_dom
      INTEGER, PRIVATE                              :: maxNbC, maxNbC_dom

      ! MPI variables
      INTEGER, PRIVATE :: comm_shar,   shar_rank,   shar_size, &
                          comm_master, master_rank, master_size, &
                          comm_world,  world_rank,  world_size, &
                          color, ierr_mpi
      LOGICAL, PRIVATE :: lcomm, lismaster, ldosync

      ! MPI windows
      INTEGER, PRIVATE :: win_vertex, win_tet, win_tet_cen, &
                          win_tet_vol, win_tet_edge,  &
                          win_state_dex, win_state_type, &
                          win_constant_mu, win_m, win_Mrem, &
                          win_Happ, win_constant_mu_o
      ! box division variables
      INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: mydom,   mypdom,  outmydom
      INTEGER, PRIVATE                            :: domsize, pdomsize,odomsize

      ! verbose and debug variables
      LOGICAL, PRIVATE                    :: lverb, ldebugm, ldebugs, ldebugt

!------------------------------------------------------------------------------
!     Subroutines
!       Main flow
!         mumaterial_setup:     Sets up MPI communicators (optional)
!         mumaterial_load:      Loads magnetic material file and sets up MPI stuff
!         mumaterial_setd:      Sets default values
!         mumaterial_setverb:   Sets standard verbosity
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
!         mumaterial_gethdipole:    Calculates dipole field at point from tet
!
!       MPI 
!         mumaterial_split:            Divides domain amongst shar_mem nodes
!         mumaterial_sync_array2d_dbl: Syncs any 2D,DBL array on shar_mem nodes
!         mumaterial_syncM:       Syncs (3,domsize) DBL array on shar_mem nodes
!         mumaterial_free:             Frees MPI memory
!       Output
!         mumaterial_output:  Output B-field and points to file
!         mumaterial_getb:    Calculates magnetic field in space
!             mumaterial_getb_scalar:      Single point in space
!               mumaterial_getbmag_scalar: Excludes applied field
!             mumaterial_getb_vector: Multiple points in space
!
!       Debug
!         mumaterial_debug:      Sets debug verbosity
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
! mumaterial_debug: Enable writing extra information to screen and writing
!                   additional files to folder
!                   [description of input parameters assume debug is enabled]
!------------------------------------------------------------------------------
! param[in]: ldebugmaster: should be .TRUE. if thread has rank 0 in world comm
! param[in]: ldebugsubmaster: should be .TRUE. if thread had rank 0 in sharcomm
! param[in]: ldebugthread: should be .TRUE.
! Set all to .FALSE. to disable debug output
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_debug(ldebugmaster, ldebugsubmaster, ldebugthread)

      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: ldebugmaster, ldebugsubmaster, ldebugthread

      ldebugm = ldebugmaster
      ldebugs = ldebugsubmaster
      ldebugt = ldebugthread
      RETURN

      END SUBROUTINE mumaterial_debug

!------------------------------------------------------------------------------
! mumaterial_free: Deallocates memory and destroys MPI windows
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_free()

      IMPLICIT NONE

      INTEGER :: ik
      IF (ASSOCIATED(state_dex))     CALL free_mpi_array1d_int(win_state_dex,state_dex,.TRUE.)
      IF (ASSOCIATED(state_type))    CALL free_mpi_array1d_int(win_state_type,state_type,.TRUE.)
      IF (ASSOCIATED(constant_mu))   CALL free_mpi_array1d_dbl(win_constant_mu,constant_mu,.TRUE.)
      IF (ASSOCIATED(constant_mu_o)) CALL free_mpi_array1d_dbl(win_constant_mu_o,constant_mu_o,.TRUE.)
      IF (ASSOCIATED(tet))           CALL free_mpi_array2d_int(win_tet,tet,.TRUE.)
      IF (ASSOCIATED(vertex))        CALL free_mpi_array2d_dbl(win_vertex,vertex,.TRUE.)
      IF (ASSOCIATED(tet_cen))       CALL free_mpi_array2d_dbl(win_tet_cen,tet_cen,.TRUE.)
      IF (ASSOCIATED(tet_vol))       CALL free_mpi_array1d_dbl(win_tet_vol,tet_vol,.TRUE.)
      IF (ASSOCIATED(tet_edge))      CALL free_mpi_array1d_dbl(win_tet_edge,tet_edge,.TRUE.)
      IF (ASSOCIATED(M))             CALL free_mpi_array2d_dbl(win_M,M,.TRUE.)
      ! TODO: Remove once allocated locally (Make sure code works beforehand)
      IF (ASSOCIATED(Mrem))          CALL free_mpi_array2d_dbl(win_Mrem,Mrem,.TRUE.)
!      IF (ASSOCIATED(Happ_shar))     CALL free_mpi_array2d_dbl(win_Happ,Happ_shar,.TRUE.)

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

#if defined(MPI_OPT)
      USE mpi
      USE mpi_params
#endif

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
!       mumaterial_setd: Sets default values
!------------------------------------------------------------------------------
! param[in]: mE. dMmax: threshold for determining convergence
! param[in]: mI. maxIter: max amount of iterations
! param[in]: la. lambdaStart: initial value of lambda
! param[in]: laF. lambdaFactor: multiplicative factor for lambda
! param[in]: laT. lambdaThresh: amount of dM>0 before lambda is multiplied
! param[in]: padF. padFactor: factor for sphere around tets for neighbours
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_setd(mE, mI, la, laF, laT, padF)
 
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in) :: mE, la, laF, padF
      INTEGER, INTENT(in) :: mI, laT

      dMmax = mE
      maxIter = mI
      lambdaStart = la
      lambdaFactor = laF
      lambdaThresh = laT
      padFactor = padF

      RETURN

      END SUBROUTINE mumaterial_setd

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

#if defined(MPI_OPT)
      USE mpi
#endif
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
      shar_rank = 0; master_rank = 0; master_size = 1; world_size = 1
      lismaster = .TRUE.; ldosync = .FALSE. 

      ! Set up MPI parameters properly now
#if defined(MPI_OPT)
      IF (lcomm) THEN
        lismaster = .FALSE.; master_rank = 1
        CALL MPI_COMM_RANK( comm_world, world_rank, ierr_mpi)
        CALL MPI_COMM_SIZE( comm_world, world_size, ierr_mpi)
        CALL MPI_COMM_RANK( comm_shar,  shar_rank,  ierr_mpi )
        IF (shar_rank.eq.0) THEN
          CALL MPI_COMM_RANK( comm_master, master_rank, ierr_mpi )
          CALL MPI_COMM_SIZE( comm_master, master_size, ierr_mpi )
          lismaster = (master_rank.EQ.0)
        END IF
        CALL MPI_Bcast( master_size, 1, MPI_INTEGER, 0, comm_shar, ierr_mpi)
        ldosync = (master_size.GE.2) 
      END IF
#endif

      mu0 = 16.0D-7 * ATAN(1.d0)

      ! Nullify pointers
      NULLIFY(vertex, tet, tet_cen, tet_vol, tet_edge, state_dex, state_type, &
              constant_mu, constant_mu_o, Mrem, M, Happ)

      ! open file, return if fails
      iunit = 327; istat = 0
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
        CALL mpialloc_2d_dbl(vertex,3,nvertex,    shar_rank,0,comm_shar,win_vertex)
        CALL mpialloc_2d_int(tet,4,ntet,          shar_rank,0,comm_shar,win_tet)
        CALL mpialloc_2d_dbl(tet_cen,3,ntet,      shar_rank,0,comm_shar,win_tet_cen)
        CALL mpialloc_1d_dbl(tet_vol,ntet,        shar_rank,0,comm_shar,win_tet_vol)
        CALL mpialloc_1d_dbl(tet_edge,ntet,       shar_rank,0,comm_shar,win_tet_edge)
        CALL mpialloc_1d_int(state_dex,ntet,      shar_rank,0,comm_shar,win_state_dex)
        CALL mpialloc_1d_int(state_type,nstate,   shar_rank,0,comm_shar,win_state_type)
        CALL mpialloc_1d_dbl(constant_mu,nstate,  shar_rank,0,comm_shar,win_constant_mu)
        CALL mpialloc_1d_dbl(constant_mu_o,nstate,shar_rank,0,comm_shar,win_constant_mu_o)
        CALL mpialloc_2d_dbl(M,            3,ntet,shar_rank,0,comm_shar,win_m)
        CALL mpialloc_2d_dbl(Mrem,3,ntet,         shar_rank,0,comm_shar,win_Mrem)  ! TODO: Allocate locally
!       CALL mpialloc_2d_dbl(Happ_shar,    3,ntet,shar_rank,0,comm_shar,win_Happ)
        ALLOCATE(stateFunction(nstate))
      ELSE
#endif
         ! if no MPI, allocate everything on one node
         ALLOCATE(vertex(3,nvertex),tet(4,ntet),state_dex(ntet), &
                  state_type(nstate),constant_mu(nstate), &
                  tet_cen(3,ntet),tet_vol(ntet),tet_edge(ntet),M(3,ntet), &
                  constant_mu_o(nstate),Mrem(3,ntet),stateFunction(nstate), &
                  STAT=istat)
#if defined(MPI_OPT)
      END IF
#endif
      IF (istat/=0) RETURN

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
        CALL MPI_Bcast(Mrem,         3*ntet,   MPI_DOUBLE_PRECISION,0,comm_master,ierr_mpi) ! TODO: Remove once allocated locally (make sure code works beforehand)
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
      CALL MUMATERIAL_SETD(1.0d-5, 100, 0.7d0, 0.75d0, 10, 1.d0)

      RETURN

      END SUBROUTINE mumaterial_load

!------------------------------------------------------------------------------
! mumaterial_info: Prints info to iunit
!------------------------------------------------------------------------------
! param[in]: iunit. Unit number to print to
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_info(iunit)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iunit
      INTEGER :: i,k

      WRITE(iunit,'(A)')           ' ---------- MUMAT MPI ----------'
      WRITE(iunit,'(3X,A,I7)')     'MPI Nodes    : ',master_size
      WRITE(iunit,'(3X,A,I7)')     'MPI Threads  : ',world_size
      WRITE(iunit,'(A)')           ' -----  Magnetic Material  -----'
      WRITE(iunit,'(3X,A,A)')      'Model Name   : ',TRIM(machine_string)
      WRITE(iunit,'(3X,A,A)')      'Date         : ',TRIM(date)
      WRITE(iunit,'(3X,A,I7)')     'Vertices     : ',nvertex
      WRITE(iunit,'(3X,A,I7)')     'Tetrahedrons : ',ntet
      WRITE(iunit,'(3X,A,I7)')     'State Funcs. : ',nstate
      WRITE(iunit,'(3X,A,EN12.3)') 'Pad factor   : ',padFactor
      WRITE(iunit,'(3X,A,I7)')     'Max Iter.    : ',maxIter
      WRITE(iunit,'(3X,A,EN12.3)') 'Max Error    : ',dMmax
      WRITE(iunit,'(3X,A,EN12.3)') 'Lambda start : ',lambdaStart
      WRITE(iunit,'(3X,A,EN12.3)') 'Lambda fact. : ',lambdaFactor
      WRITE(iunit,'(3X,A,I7)')     'Lambda thrsh.: ',lambdaThresh
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
! fcn           : getBfld. Function which returns the vacuum magnetic field
!                 SUBROUTINE FCN(x,y,z,bx,by,bz)
! param[in]: offset. Offset of all tiles from the origin
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_init_new(getBfld, offset)

#if defined(MPI_OPT)
      USE mpi
      USE mpi_params
#endif

      IMPLICIT NONE
      
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: offset(3)
      INTEGER :: mystart, myend, ourstart, ourend
      LOGICAL :: lwork
      INTEGER :: i, j, k, i_tile, j_tile
      INTEGER :: mstat(MPI_STATUS_SIZE)
      CHARACTER(LEN=6) :: strcount, splitcount

      DOUBLE PRECISION :: Bx, By, Bz
      DOUBLE PRECISION :: tol, delta, xmin, xmax, ymin, ymax, zmin, zmax, pad
      INTEGER :: splits, dim, ydomsize, reci
      INTEGER, ALLOCATABLE :: domin(:), yourdom(:), tdom(:), idx(:)

      EXTERNAL:: getBfld
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Apply offset
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (PRESENT(offset).AND.(MAXVAL(ABS(offset)) .gt. 0.d0)) THEN
        IF (lverb) WRITE(6,*) "  MUMAT_INIT:  Applying offset to vertices"
        mystart = 1; myend = nvertex
#if defined(MPI_OPT)
        IF (lcomm) CALL MPI_CALC_MYRANGE(comm_world, 1, nvertex, mystart, myend)
#endif
        DO i = mystart, myend
          vertex(:,i) = vertex(:,i) + offset
        END DO
            
#if defined(MPI_OPT)
        IF (ldosync) THEN
          IF (lverb) WRITE(6,*) "  MUMAT_DEBUG:  Synchronising offset vertices"
          CALL mumaterial_sync_array2d_dbl(vertex,3,nvertex,mystart,myend)
        END IF
#endif
        IF (ldebugm) CALL mumaterial_writedebug(vertex,3,nvertex, 'verts.dat','vertices')
      END IF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Calculate tet centers, edges, volumes, then synchronize
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (lverb) WRITE(6,*) "  MUMAT_INIT:  Calculating tetrahedron quantities"; FLUSH(6)
      mystart = 1; myend = ntet   

#if defined(MPI_OPT)
      IF (lcomm) THEN 
        CALL MPI_CALC_MYRANGE(comm_world, 1, ntet, mystart, myend) 
        tet_cen(:,mystart:myend) = 99999.0  ! these values will be overwritten;
        tet_vol(mystart:myend)   = 99999.0  ! big numbers make problems obvious
        tet_edge(mystart:myend)  = 99999.0
      END IF
#endif

      DO i = mystart, myend
        tet_cen(:,i) = (vertex(:,tet(1,i)) + vertex(:,tet(2,i)) + &
                        vertex(:,tet(3,i)) + vertex(:,tet(4,i)))/4.d0
        tet_vol(i) = mumaterial_gettetvolume(vertex(:,tet(1,i)),vertex(:,tet(2,i)), &
                                             vertex(:,tet(3,i)),vertex(:,tet(4,i)))
        tet_edge(i) = SQRT(6.0)/4.d0*(6.d0*SQRT(2.0)*tet_vol(i))**(1.0/3.0) 
      END DO

#if defined(MPI_OPT)
      IF (ldosync) THEN
        IF (lverb) WRITE(6,*) "  MUMAT_INIT:  Synchronising tetrahedron quantities"; FLUSH(6)
        CALL mumaterial_sync_array2d_dbl(tet_cen, 3,ntet,mystart,myend)
        CALL mumaterial_sync_array2d_dbl(tet_vol, 1,ntet,mystart,myend)
        CALL mumaterial_sync_array2d_dbl(tet_edge,1,ntet,mystart,myend)
      END IF
#endif
      IF (ldebugm) THEN
        CALL mumaterial_writedebug(tet_cen, 3, ntet, 'tet_cen.dat')
        CALL mumaterial_writedebug(tet_vol, 1, ntet, 'tet_cen.dat')
        CALL mumaterial_writedebug(tet_edge,1, ntet, 'tet_cen.dat')
      END IF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Domain split
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      color = 0
#if defined(MPI_OPT)
      IF (shar_rank.NE.0) THEN
        color = 1
      ELSE
        Bx = master_size ! master_size needs to be dbl for next line
        color = world_rank*Bx/world_size
        WRITE(strcount, '(I0)') color
      END IF
#endif

      lwork = (color.EQ.0)
      IF (lwork) THEN ! Global master, create first box
        ALLOCATE(mydom(ntet))
        DO i = 1, ntet
          mydom(i) = i
        END DO
        domsize = SIZE(mydom)
        IF (.NOT.(ldosync)) THEN
          ALLOCATE(mypdom(domsize))
          mypdom = mydom
          pdomsize = domsize
        END IF
      END IF

#if defined(MPI_OPT)   
      IF (ldosync.AND.(shar_rank.EQ.0)) THEN                  
        IF (ldebugs) WRITE(6,'(A22,I3,A11,L1)') '  MUMAT_DEBUG: MASTER ', color,': LWORK is ', lwork; FLUSH(6)

        splits = NINT(LOG(Bx)/LOG(2.0)) ! log_2(X) = ln(X)/log(2)
        tol = 0.001
        delta = 1.0

        DO
          IF (splits.EQ.0) THEN
            IF (ldebugs) WRITE(6,'(A22,I3,A13)') '  MUMAT_DEBUG: MASTER ', color,' EXITING LOOP'; FLUSH(6)
            EXIT ! Reached end
          END IF

          IF (lwork) THEN 
            ALLOCATE(domin(domsize))
            domin = mydom
            DEALLOCATE(mydom)

            CALL mumaterial_split(domsize, domin, ntet, tet_cen, tol, delta, mydom, yourdom) ! Split box
            DEALLOCATE(domin)
            domsize  = SIZE(mydom)
            ydomsize = SIZE(yourdom)            
            splits = splits-1 

            IF (ldebugs) THEN
              WRITE(splitcount, '(I0)') splits
              OPEN(color, file='./mydom_' // TRIM(ADJUSTL(strcount)) // '_' // TRIM(ADJUSTL(splitcount)) // '.dat')
              DO i = 1, domsize
                WRITE(color, "(I8)") mydom(i)
              END DO
              CLOSE(color)
              OPEN(color, file='./yourdom_' // TRIM(ADJUSTL(strcount)) // '_' // TRIM(ADJUSTL(splitcount)) // '.dat')
              DO i = 1, ydomsize
                WRITE(color, "(I8)") yourdom(i)
              END DO
              CLOSE(color)
            END IF

            IF (ldebugs) WRITE(6,'(A22,I3,A14,I2)') '  MUMAT_DEBUG: MASTER ', color,' splits left: ', splits; FLUSH(6)

            ! now mail one of new boxes to the appropriate recipient
            reci = color + 2**splits 
            CALL MPI_SEND(ydomsize,      1, MPI_INTEGER, reci, 1234, comm_master, ierr_mpi) 
            CALL MPI_SEND(yourdom,ydomsize, MPI_INTEGER, reci, 1235, comm_master, ierr_mpi);  DEALLOCATE(yourdom) 
            CALL MPI_SEND(splits,        1, MPI_INTEGER, reci, 1236, comm_master, ierr_mpi)
          ELSE
            CALL MPI_RECV(domsize,     1, MPI_INTEGER, MPI_ANY_SOURCE, 1234, comm_master, mstat, ierr_mpi); ALLOCATE(mydom(domsize))
            CALL MPI_RECV(mydom, domsize, MPI_INTEGER, MPI_ANY_SOURCE, 1235, comm_master, mstat, ierr_mpi)
            CALL MPI_RECV(splits,      1, MPI_INTEGER, MPI_ANY_SOURCE, 1236, comm_master, mstat, ierr_mpi);  
            IF (ldebugs) WRITE(6,'(A22,I3,A28,I8,A1)') '  MUMAT_DEBUG: MASTER ', color,' received box [', domsize, ']'; FLUSH(6)
            lwork = .TRUE. ! Activate node
          END IF
        END DO

        IF (ldebugs) THEN
          WRITE(strcount, '(I0)') color
          OPEN(color, file='./mydom_' // TRIM(ADJUSTL(strcount)) // '.dat')
          DO i = 1, domsize
            WRITE(color, "(I8)") mydom(i)
          END DO
          CLOSE(color)
        END IF
        ! Now every master needs to know their range relative to other 
        lwork = (color.EQ.0)
        splits = 0
        DO
          IF (lwork) THEN
            ourstart = splits+1
            ourend   = splits+domsize
            IF (ldebugs) WRITE(6,'(A22,I3,A12,I8,I8,A1)') '  MUMAT_DEBUG: MASTER ', color,' has range [', ourstart, ourend, ']'; FLUSH(6)
            reci = color + 1
            IF (reci.EQ.master_size) EXIT
            CALL MPI_SEND(ourend, 1, MPI_INTEGER, reci, 1234, comm_master, ierr_mpi)             
            EXIT
          ELSE
            CALL MPI_RECV(splits, 1, MPI_INTEGER, MPI_ANY_SOURCE, 1234, comm_master, mstat, ierr_mpi)
            lwork = .TRUE.
          END IF
        END DO
        
        ! Construct outside domain (used for synchronization)
        odomsize = ntet-domsize
        ALLOCATE(outmydom(odomsize))
        i = 0
        DO i_tile = 1, ntet
          IF (.NOT.(ANY(i_tile==mydom))) THEN 
            i = i + 1
            outmydom(i) = i_tile
          END IF
        END DO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Create padded domains
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        xmin = tet_cen(1,mydom(1)); xmax = xmin
        ymin = tet_cen(2,mydom(1)); ymax = ymin
        zmin = tet_cen(3,mydom(1)); zmax = zmin
        pad = 0

        DO i = 2, domsize
          xmin = MIN(tet_cen(1,mydom(i)),xmin); xmax = MAX(tet_cen(1,mydom(i)),xmax)
          ymin = MIN(tet_cen(2,mydom(i)),ymin); ymax = MAX(tet_cen(2,mydom(i)),ymax)
          zmin = MIN(tet_cen(3,mydom(i)),zmin); zmax = MAX(tet_cen(3,mydom(i)),zmax)
          pad = MAX(padFactor*tet_edge(mydom(i)),pad)
        END DO
  
        ALLOCATE(tdom(ntet))
        pdomsize = 0
        DO i = 1, ntet
          IF (((tet_cen(1,i)>=xmin-pad).AND.(tet_cen(1,i)<=xmax+pad)) .AND. &
              ((tet_cen(2,i)>=ymin-pad).AND.(tet_cen(2,i)<=ymax+pad)) .AND. & 
              ((tet_cen(3,i)>=zmin-pad).AND.(tet_cen(3,i)<=zmax+pad))) THEN
                pdomsize = pdomsize + 1
                tdom(pdomsize) = i
          END IF
        END DO
        ALLOCATE(mypdom(pdomsize))
        mypdom = tdom(1:pdomsize)
        DEALLOCATE(tdom)

        IF (ldebugs) THEN
          OPEN(color, file='./mypdom_' // TRIM(ADJUSTL(strcount)) // '.dat')
          DO i = 1, pdomsize
            WRITE(color, "(I8)") mypdom(i)
          END DO
          CLOSE(color)
        END IF
      END IF

      IF (ldebugs) WRITE(6,'(A22,I3,A19,I8,A1)') '  MUMAT_DEBUG: MASTER ', color,' broadcasting box [', domsize, ']'; FLUSH(6)
      CALL MPI_Bcast(domsize,    1, MPI_INTEGER, 0, comm_shar, ierr_mpi)
      CALL MPI_Bcast(pdomsize,   1, MPI_INTEGER, 0, comm_shar, ierr_mpi)
      
      IF (shar_rank.NE.0) THEN
        ALLOCATE(mydom(domsize))
        ALLOCATE(mypdom(pdomsize))
        odomsize = ntet-domsize
        ALLOCATE(outmydom(odomsize))
      END IF
      
      CALL MPI_Bcast(mydom,   domsize,  MPI_INTEGER, 0, comm_shar, ierr_mpi)
      CALL MPI_Bcast(mypdom,  pdomsize, MPI_INTEGER, 0, comm_shar, ierr_mpi)
      CALL MPI_Bcast(outmydom,odomsize, MPI_INTEGER, 0, comm_shar, ierr_mpi)
!      IF (shar_rank.EQ.1) WRITE(6,'(A37,I8,A1)') '  MUMAT_DEBUG: SUBJECT received box [', domsize, ']'; FLUSH(6)
      CALL MPI_Bcast(ourstart, 1, MPI_INTEGER, 0, comm_shar, ierr_mpi)
      CALL MPI_Bcast(ourend,   1, MPI_INTEGER, 0, comm_shar, ierr_mpi)
      CALL MPI_BARRIER(comm_world, ierr_mpi)
#endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Determine nearest Nb (includes self)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (lverb) WRITE (6,*) "  MUMAT_INIT:  Determining nearest Nb"
      IF (lcomm) CALL MPI_CALC_MYRANGE(comm_shar, 1, domsize, mystart, myend)
      CALL mumaterial_getneighbours(mystart, myend)
      IF (ldebugt) WRITE(6,'(3X,A13,I6,A12,I8,I8,A3,I6,I6,A1)') 'MUMAT_DEBUG: ', world_rank, ' NB RANGE: [', MINVAL(NbC), maxNbC, '] [',MINLOC(NbC,1), MAXLOC(NbC,1) ,']'

#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_BARRIER(comm_world, ierr_mpi)
#endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate H_app
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (lverb) WRITE (6,*) "  MUMAT_INIT:  Calculating H_app"
      NULLIFY(Happ)
      ALLOCATE(Happ(3,mystart:myend))
      Happ(:,:) = 0.0
      M(:,:) = 0.0
      DO i = mystart, myend
        i_tile = mydom(i)
        CALL getBfld(tet_cen(1,i_tile), tet_cen(2,i_tile), tet_cen(3,i_tile), Bx, By, Bz)
        Happ(:,i) = [Bx/mu0, By/mu0, Bz/mu0]
      END DO

#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_BARRIER(comm_world, ierr_mpi)
#endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate N_store
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (lverb) WRITE (6,*) "  MUMAT_INIT:  Calculating N_store"

      NULLIFY(N_store)
      ALLOCATE(N_store(3,3,maxNbC,mystart:myend))
      N_store(:,:,:,:) = 0.0
      DO i = mystart, myend
        i_tile = mydom(i)
        DO j = 1, NbC(i)
          j_tile = Nb(j,i)
          CALL mumaterial_getN(vertex(:,tet(1,j_tile)), vertex(:,tet(2,j_tile)), vertex(:,tet(3,j_tile)), vertex(:,tet(4,j_tile)), tet_cen(:,i_tile), N_store(:,:,j,i)) 
        END DO 
      END DO

#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_BARRIER(comm_world, ierr_mpi)
#endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Begin iterations
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (lverb) WRITE (6,*) "  MUMAT_INIT:  Beginning Iterations"
      
      CALL mumaterial_iterate_M(getBfld, mystart, myend)

      IF (lverb) WRITE (6,*) "  MUMAT_INIT:  End Iterations"

      ! DEALLOCATE Helpers
      DEALLOCATE(Nb, NbC)
      ! DEALLOCATE( Nb_domidx, NbC_dom)
      DEALLOCATE(N_store)
      DEALLOCATE(Happ)

      RETURN
      END SUBROUTINE mumaterial_init_new

!-----------------------------------------------------------------------
! mumaterial_iterate_M: Iteration loop
!-----------------------------------------------------------------------
! param[in]: mystart, myend. range of tetrahedrons worked on by thread
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_iterate_M(getBfld, mystart, myend)

#if defined(MPI_OPT)
      USE mpi
      USE mpi_params
#endif

      IMPLICIT NONE

      INTEGER, INTENT(in) :: mystart, myend

      INTEGER :: count, i, i_tile, j, j_tile, k, k_tile
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: chi
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: M_new
      DOUBLE PRECISION :: H(3), N(3,3), Bx, By, Bz
      DOUBLE PRECISION :: H_old(3), H_new(3),  lambda_s,  Hnorm, M_tmp_norm
      DOUBLE PRECISION :: M_tmp(3), M_tmp_local(3), Mrem_norm, u_ea(3), u_oa_1(3), u_oa_2(3) ! hard magnet

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Mnorm, MnormPrev, dM, dMPrev, d2M, d2MPrev
      DOUBLE PRECISION :: lambda, maxdM, maxdMPrev, lambdaIncrF, lambdaMax
      INTEGER          :: lambdaCount, maxd2MC, lastldM, ldMCount
      INTEGER, DIMENSION(:), ALLOCATABLE :: d2MCount
      LOGICAL          :: ldM

      EXTERNAL:: getBfld

      CHARACTER(LEN=6) :: strcount
      CHARACTER(LEN=20) :: filename
    
      ! Allocate helpers
      ALLOCATE(M_new(3,mystart:myend),chi(mystart:myend),Mnorm(mystart:myend),MnormPrev(mystart:myend))
      ALLOCATE(dM(mystart:myend),dMPrev(mystart:myend))
      ALLOCATE(d2M(mystart:myend),d2MPrev(mystart:myend))
      ALLOCATE(d2MCount(mystart:myend))

      count = 0
      lambda = lambdaStart
      lambdaCount = 0
      chi = 0.0
      Mnorm = 1.0E-5
      dM = 0.d0
      d2M = 0.d0
      d2MCount = 0
      lastldM = 0
      ! hard-coded, but could always change later.
!      ldMCount = 5
!      lambdaIncrF = 1.05d0
!      lambdaMax = 0.95d0

      IF (lverb) THEN
        WRITE(6,*) ''
        WRITE(6,*) '  Count     Max Error        Target        Lambda   LC  d2MC'
        WRITE(6,*) '============================================================================'
      END IF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Main Iteration Loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO
        count = count + 1        

        MnormPrev = Mnorm
        M_new = 0.0
        dMPrev = dM
        d2MPrev = d2M
        dM = 0.d0
        d2M = 0.d0
        maxdMPrev = maxdM
        maxdM = 0.d0
        maxd2MC = 0
        ldM = .FALSE.
        
        DO i = mystart, myend ! Get the field and new magnetization for each tile
          i_tile = mydom(i)
          H = Happ(:,i)

          DO j = 1, NbC(i)  ! Full field if neighbour
            j_tile = Nb(j,i)
            IF (j_tile.EQ.i_tile) THEN
              N = N_store(:,:,j,i)
              CYCLE
            END IF
            H = H + MATMUL(N_store(:,:,j,i), M(:,j_tile))
          END DO
          
          ! Determine field and magnetization at tile due to all other tiles and itself
          H_new = H
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

                  IF (MAXVAL(ABS((H_new - H_old)/H_old)) .lt. dMmax*lambda_s) THEN
                    M_new(:,i) = (Mrem_norm + (constant_mu(state_dex(i_tile)) - 1) * DOT_PRODUCT(H_new, u_ea)) * u_ea &
                                                  + (constant_mu_o(state_dex(i_tile)) - 1) * DOT_PRODUCT(H_new, u_oa_1) * u_oa_1 &
                                                  + (constant_mu_o(state_dex(i_tile)) - 1) * DOT_PRODUCT(H_new, u_oa_2) * u_oa_2
                    EXIT
                  END IF
              END DO
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            CASE (2) ! Soft magnet using state function
              DO
                H_old = H_new
                Hnorm = NORM2(H_new)
                IF (Hnorm .ne. 0) THEN
                  CALL mumaterial_getState(stateFunction(state_dex(i_tile))%H, stateFunction(state_dex(i_tile))%M, Hnorm, M_tmp_norm)
                  M_tmp = M_tmp_norm * H_new / Hnorm
                  lambda_s = MIN(Hnorm/M_tmp_norm, 0.5)
                ELSE
                  M_tmp = 0
                  M_tmp_norm = 0
                  lambda_s = 0.5
                END IF
                H_new = H + MATMUL(N, M_tmp)
                H_new = H_old + lambda_s * (H_new - H_old)

                IF (MAXVAL(ABS((H_new - H_old)/H_old)) .lt. dMmax*lambda_s) THEN
                  Hnorm = NORM2(H_new)
                  IF (Hnorm .ne. 0) THEN
                    CALL mumaterial_getState(stateFunction(state_dex(i_tile))%H, stateFunction(state_dex(i_tile))%M, Hnorm, M_tmp_norm)
                    M_new(:,i) = M_tmp_norm * H_new / Hnorm
                    chi(i) = M_tmp_norm / Hnorm
                  ELSE
                    M_new(:,i) = 0
                    chi(i) = 0
                  END IF
                  EXIT
                END IF
              END DO
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            CASE (3) ! Soft magnet using constant permeability
              lambda_s = MIN(1/constant_mu(state_dex(i_tile)), 0.5)
              DO
                  H_old = H_new
                  H_new = H + (constant_mu(state_dex(i_tile)) - 1) * MATMUL(N, H_new)
                  H_new = H_old + lambda_s * (H_new - H_old)
                  IF (MAXVAL(ABS((H_new - H_old)/H_old)) .lt. dMmax*lambda_s) THEN
                    M_new(:,i) = (constant_mu(state_dex(i_tile)) - 1) * H_new
                    EXIT
                  END IF
              END DO
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            CASE DEFAULT
              WRITE(6,*) "Unknown magnet type: ", state_type(state_dex(i_tile))
              STOP
          END SELECT

          M(:,i_tile) = M(:,i_tile) + lambda*(M_new(:,i) - M(:,i_tile))
          Mnorm(i) = NORM2(M(:,i_tile))

          ! "Derivatives" for convergence checks
          dM(i) = ABS((Mnorm(i) - MnormPrev(i))/MnormPrev(i))
          IF (dMPrev(i).GT.1E-12) THEN
            d2M(i) = ABS((dM(i)-dMPrev(i))/dMPrev(i))
            IF ((d2M(i).GT.d2MPrev(i)).AND.(dM(i).GT.d2MPrev(i))) THEN 
              d2MCount(i) = d2MCount(i)+1
            ELSE
              d2MCount(i) = 0
            END IF
          END IF
          maxdM = MAX(dM(i), maxdM) 
          maxd2MC = MAX(d2MCount(i), maxd2MC) 

          ! Dampen evolution when Mnorm has increased AND the rate of increase
          ! has itself increased for a couple of consecutive iters
          IF ((Mnorm(i).GT.MnormPrev(i)) .AND. (dM(i).GT.dMPrev(i)) .AND. (d2MCount(i).GE.3)) THEN 
            d2MCount(i) = 0
            ldM = .TRUE.     
          END IF     
        END DO

#if defined(MPI_OPT)
        ! Synchronise dM
        IF (lcomm) THEN
          CALL MPI_BARRIER(comm_world, ierr_mpi)
          CALL MPI_ALLREDUCE(MPI_IN_PLACE, maxdM, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_world, ierr_mpi)
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,maxd2MC,1, MPI_INTEGER,          MPI_MAX, comm_world, ierr_mpi)
          CALL MPI_ALLREDUCE(MPI_IN_PLACE, ldM,   1, MPI_LOGICAL,          MPI_LOR, comm_world, ierr_mpi) 
        END IF
#endif
        IF (maxdM.LT.maxdMPrev) THEN
          ldM = .FALSE.
        END IF

        IF (lverb) THEN 
          WRITE(6,'(3X,I5,A2,E12.4,A2,E12.4,A2,E12.4,A2,I3,A2,I3)') count, '  ', maxdM, '  ', dMmax*lambda,  '  ', lambda, '  ', lambdaCount, '  ', maxd2MC
          CALL FLUSH(6)
        END IF

        ! Update lambda to prevent oscillations
        IF (count.NE.1) THEN
          IF (ldM) THEN
            lambdaCount = lambdaCount + 1
!            lastldM = 0
            IF (lambdaCount.EQ.lambdaThresh) THEN
              lambda = lambda * lambdaFactor
              lambdaCount = 0
            END IF
          ELSE
            lambdaCount= MAX(lambdaCount-1,0)
!            lastldM = lastldM+1
!            IF ((lastldM.GT.0).AND.(MODULO(lastldM,ldMCount).EQ.0)) THEN
!              lambda = MIN(lambda*lambdaIncrF,lambdaMax)
!            END IF
          END IF
        END IF

        IF (ldebugm) THEN
            WRITE(strcount, '(I0)') count
            CALL mumaterial_writedebug(M,3,ntet,'./M_' // TRIM(ADJUSTL(strcount)) // '.dat')
        END IF

        IF (count.GE.maxIter) THEN
          IF (lverb) WRITE(6,*) "  MUMAT:  Stopping - Exceeded maximum iterations"
          EXIT
        END IF
        IF ( (maxdM.LT.dMmax*lambda) .AND.(count.GT.1)) THEN
          IF (lverb) WRITE(6,*) "  MUMAT:  Stopping - converged"
          EXIT
        END IF
        IF (lambda.LT.1E-4) THEN
          IF (lverb) WRITE(6,*) "  MUMAT:  Stopping - lambda < 1e-4"
          EXIT
        END IF

        ! Synchronize magnetization
        IF (ldosync) CALL mumaterial_syncM(M,ntet,outmydom)
        ! Update H-field from non-Nb
        DO i = mystart, myend
          i_tile = mydom(i)
          CALL getBfld(tet_cen(1,i_tile), tet_cen(2,i_tile), tet_cen(3,i_tile), Bx, By, Bz)
          Happ(:,i) = [Bx/mu0, By/mu0, Bz/mu0]
        
        ! Happ loop logic
        ! ----------
        ! Tetrahedron array e.g.  T=[1 2 ... 12407 12408]
        ! Neighbor array 1  e.g. Nb=[6 48 3874 4838 6792 11240]
        ! 
        ! Want to iterate over all tetrahedrons in T that do not appear in Nb
        ! IF statements for loop over ntet elements is slow so we do this instead
        !
        ! 1. Loop over elements before first neighbour
        !     Nb(1)=6                         => loop over [1 2 3 4 5]
        ! 2. For each pair of neighbours, loop over elements IN BETWEEN 
        !     For j=2, Nb(j-1)=6,  Nb(j)=48   => loop over [7 8 ... 46 47]
        !     For j=3, Nb(j-1)=48, Nb(j)=3874 => loop over [49 ... 3873] etc.
        ! 3. Loop over elements after last neighbour
        !     Nb(end)=11240                   => loop over [11241 ... 12408]

          DO k_tile = 1,Nb(1,i)-1 
              CALL mumaterial_gethdipole(tet_cen(:,k_tile),tet_cen(:,i_tile),M(:,k_tile),tet_vol(k_tile),Happ(:,i))
          END DO
          DO j = 2, NbC(i)                    
            DO k_tile = Nb(j-1,i)+1,Nb(j,i)-1
              CALL mumaterial_gethdipole(tet_cen(:,k_tile),tet_cen(:,i_tile),M(:,k_tile),tet_vol(k_tile),Happ(:,i))
            END DO
          END DO
          DO k_tile = Nb(NbC(i),i)+1,ntet     
              CALL mumaterial_gethdipole(tet_cen(:,k_tile),tet_cen(:,i_tile),M(:,k_tile),tet_vol(k_tile),Happ(:,i))
          END DO
        END DO  


      END DO
      DEALLOCATE(M_new,chi,Mnorm,MnormPrev,dM,dMPrev,d2MCount)

      RETURN
      END SUBROUTINE mumaterial_iterate_M


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

      FUNCTION mumaterial_gettetvolume(v1,v2,v3,v4)
            IMPLICIT NONE
            DOUBLE PRECISION, DIMENSION(3), INTENT(in) :: v1, v2, v3, v4
            DOUBLE PRECISION :: mumaterial_gettetvolume

            mumaterial_gettetvolume = ABS(dot_product(v1-v4,mumaterial_cross(v2-v4,v3-v4)))/6.0
            RETURN

      END FUNCTION mumaterial_gettetvolume


      SUBROUTINE mumaterial_getneighbours(mystart, myend)

      INTEGER, INTENT(in) :: mystart, myend
      INTEGER :: i, j, k, c, i_tile
      DOUBLE PRECISION, ALLOCATABLE ::  dist(:), dx(:,:)
      !DOUBLE PRECISION, ALLOCATABLE :: tet_cen_pdom(:,:)
      LOGICAL, ALLOCATABLE :: mask(:)
      !INTEGER, ALLOCATABLE :: idx(:)
      !INTEGER, ALLOCATABLE :: Nb_temp(:,:), Nb_domidx_temp(:,:)

      !ALLOCATE(Nb_temp(pdomsize, mystart:myend))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Nb_temp = 0
      
      !ALLOCATE(mask(pdomsize), dist(pdomsize), dx(3,pdomsize),tet_cen_pdom(3, pdomsize))
      !DO i = 1, pdomsize
      !  tet_cen_pdom(:,i) = tet_cen(:,mypdom(i))
      !END DO
      ALLOCATE(NbC(mystart:myend),mask(ntet),dist(ntet),dx(3,ntet))
      NbC = 0
      
      ! Get largest neighbor count for allocation first
      DO i = mystart, myend
        i_tile = mydom(i)
        dx(1,:) = tet_cen(1,:) - tet_cen(1,i_tile)
        dx(2,:) = tet_cen(2,:) - tet_cen(2,i_tile)
        dx(3,:) = tet_cen(3,:) - tet_cen(3,i_tile)
        dist = NORM2(dx,DIM=1)
        NbC(i) = COUNT(dist.LE.padFactor*tet_edge)
      END DO

      maxNbC = MAXVAL(NbC)
      ALLOCATE(Nb(maxNbC,mystart:myend))

      ! Actual neighbour loop
      DO i = mystart, myend
        i_tile = mydom(i)
        dx(1,:) = tet_cen(1,:) - tet_cen(1,i_tile)
        dx(2,:) = tet_cen(2,:) - tet_cen(2,i_tile)
        dx(3,:) = tet_cen(3,:) - tet_cen(3,i_tile)
        dist = NORM2(dx,DIM=1)
        mask = dist.LE.padFactor*tet_edge
        j = 0
        DO k = 1, ntet
          IF (mask(k)) THEN
            j = j + 1
            Nb(j,i) = k
          END IF
        END DO
      END DO
      DEALLOCATE(mask,dist,dx)

      ! In case the above needs to be changed
      ! DO i = mystart, myend
      !   i_tile = mydom(i)
      !   dx(1,:) = tet_cen_pdom(1,:) - tet_cen(1,i_tile)
      !   dx(2,:) = tet_cen_pdom(2,:) - tet_cen(2,i_tile)
      !   dx(3,:) = tet_cen_pdom(3,:) - tet_cen(3,i_tile)
      !   dist = NORM2(dx,DIM=1)
      !   mask = .TRUE.
      !   !NbC(i) = COUNT(dist.LE.padFactor*tet_edge(i_tile))
      !   mask = dist.LE.padFactor*tet_edge(i_tile)
      !   NbC(i) = COUNT(mask)
      !   j = 0
      !   DO k = 1, pdomsize
      !     IF (mask(k)) THEN
      !       j = j + 1
      !       Nb_temp(j,i) = mypdom(k)
      !     END IF
      !   END DO
      ! END DO

      !maxNbC = MAXVAL(NbC)

      !DEALLOCATE(tet_cen_pdom)

      !ALLOCATE(Nb(maxNbC,mystart:myend))
      !Nb = Nb_temp(1:maxNbC,mystart:myend)
      !DEALLOCATE(Nb_temp)
      ! No longer necessary
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ALLOCATE(Nb_domidx_temp(maxNbc, mystart:myend), NbC_dom(mystart:myend))
      ! Nb_domidx_temp = 0
      ! NbC_dom = maxNbC
      ! ALLOCATE(idx(maxNbC))
      ! DO i = mystart, myend
      !   c=0
      !   DO j = 1, NbC(i)
      !       k = FINDLOC(mydom, Nb(j,i), DIM=1)
      !       IF (k.NE.0) THEN
      !         c = c + 1
      !         Nb_domidx_temp(c, i) = k
      !       END IF
      !   END DO
      !   NbC_dom(i) = c

      !   WHERE(Nb_domidx_temp(:,i).EQ.0) Nb_domidx_temp(:,i) = ntet+1
      !   idx = 0
      !   CALL SORT(maxNbC,Nb_domidx_temp(:,i),idx)

      ! END DO

      ! DEALLOCATE(idx)
      ! maxNbC_dom = MAXVAL(NbC_dom)

      ! ALLOCATE(Nb_domidx(maxNbC_dom,mystart:myend))
      ! Nb_domidx=Nb_domidx_temp(1:maxNbC_dom, mystart:myend) 
      ! DEALLOCATE(Nb_domidx_temp)

      END SUBROUTINE mumaterial_getneighbours

      SUBROUTINE mumaterial_split(boxsize,boxin,ncoords,coords,tol,delta,box1,box2)
      !-----------------------------------------------------------------------
      ! mumaterial_split: Divides a set of neighboring tetrahedrons into two
      ! approximately equally sized groups of neighboring tetrahedrons
      !-----------------------------------------------------------------------
      ! param[in]: coords. coordinates of tetrahedrons given in boxin
      ! param[in]: dim. dimension over which to split (1:X, 2:Y, 3:Z)
      ! param[in]: tol. allowed deviation from exact 50/50 split
      ! param[in]: delta. initial increment in dim
      ! param[out]: box1. collection of half the input tets.
      ! param[out]: box2. collection of other half of input tets.
      !-----------------------------------------------------------------------
#if defined(MPI_OPT)
      USE mpi
      USE mpi_params
#endif     

      IMPLICIT NONE

      INTEGER, INTENT(in) :: boxsize, ncoords
      DOUBLE PRECISION, DIMENSION(3,ncoords), INTENT(in) :: coords
      INTEGER, INTENT(in)  :: boxin(boxsize)
      DOUBLE PRECISION, INTENT(in) :: tol, delta
      INTEGER, ALLOCATABLE, INTENT(out) :: box1(:), box2(:) 

      INTEGER :: i1, i2, iter, dim, dimc, itermax, ptsinbox1, ptsinbox2, cnochange, idx(boxsize)
      DOUBLE PRECISION, ALLOCATABLE :: xyz(:,:), com(:), dx(:,:)
      DOUBLE PRECISION :: divval, prevdev, currdev, usedelta, pts_dbl, small
      DOUBLE PRECISION :: dist1, dist2, mean, prevmean
      INTEGER, ALLOCATABLE :: tbox1(:), tbox2(:)

      ALLOCATE(xyz(3,boxsize))
      DO i1 = 1, boxsize
        xyz(:,i1) = coords(:,boxin(i1))
      END DO
      itermax = 201
      small = 1E-12
      prevmean = 1E+12
      ALLOCATE(box1(1),box2(1))

      DO dim = 1, 3
        usedelta = delta
        iter = 0
        cnochange = 0
        prevdev = 2 
        divval = (MINVAL(xyz(dim,:))+MAXVAL(xyz(dim,:)))/2 ! initial estimate of dividing value
        
        DO
          iter = iter + 1 
          IF (iter.GE.itermax) EXIT ! exceeding max iteration

          ptsinbox1 = COUNT(xyz(dim,:).LT.divval); pts_dbl = ptsinbox1
          currdev = ABS(pts_dbl/boxsize-0.5) ! deviation from 50/50 split

          IF (ldebugs) WRITE(6,'(I3,A2,E15.7,A2,I9,A2,E15.7,A2,I2)') iter, ', ', divval, ', ', ptsinbox1, ', ', currdev, ', ', cnochange

          IF (ABS(currdev-prevdev).LE.small) THEN 
            cnochange = cnochange + 1
            IF (cnochange.GE.2) EXIT ! no improvement
          ELSE
            cnochange = 0
          END IF
          IF (currdev.LE.tol) EXIT ! within tolerance

          IF (currdev-prevdev>0) usedelta = -0.5*usedelta ! reverse direction, decrease step size
          divval = divval + usedelta
          prevdev = currdev
        END DO

        ptsinbox1 = COUNT(xyz(dim,:).LT.divval)
        ptsinbox2 = boxsize - ptsinbox1
        ALLOCATE(tbox1(ptsinbox1),tbox2(ptsinbox2))

        ! Make boxes
        i1 = 1; i2 = 1
        DO iter = 1, boxsize
          IF (xyz(dim,iter).LT.divval) THEN
            tbox1(i1) = boxin(iter)
            i1 = i1+1
          ELSE
            tbox2(i2) = boxin(iter)
            i2 = i2+1
          END IF
        END DO

        ! Evaluate distances from center of masses
        dist1 = 0; dist2 = 0
        ALLOCATE(com(3), dx(3, ptsinbox1)) 
        com = 0
        DO i1 = 1, ptsinbox1
          com = com + coords(:, tbox1(i1))
        END DO
        com = com/ptsinbox1
        DO i1 = 1, ptsinbox1
          dx(:,i1) = coords(:,tbox1(i1)) - com
        END DO
        dist1 = SUM(NORM2(dx, DIM=1))/ptsinbox1

        DEALLOCATE(dx)
        ALLOCATE(dx(3, ptsinbox2)) 
        com = 0
        DO i2 = 1, ptsinbox2
          com = com + coords(:, tbox2(i2))
        END DO
        com = com/ptsinbox2
        DO i2 = 1, ptsinbox2
          dx(:,i2) = coords(:,tbox2(i2)) - com
        END DO
        dist2 = SUM(NORM2(dx, DIM=1))/ptsinbox2

        mean = (dist1+dist2)/2
        IF (ldebugs) WRITE(6,'(3X,A31,I1,A1,E15.7)') "Mean distance from COM for dim=", dim, ':', mean
        ! Update if dimension is better
        IF (mean.LT.prevmean) THEN
          DEALLOCATE(box1, box2)
          ALLOCATE(box1(ptsinbox1),box2(ptsinbox2))
          box1 = tbox1
          box2 = tbox2
          prevmean = mean
          dimc = dim
        END IF
        DEALLOCATE(tbox1, tbox2, com, dx)
        
      END DO
      IF (ldebugs) WRITE(6,'(3X,A18,I1)') "Dimension chosen: ", dimc
      DEALLOCATE(xyz)

      END SUBROUTINE mumaterial_split

      SUBROUTINE mumaterial_sync_array2d_dbl(array, n1, n2, mystart,myend)

#if defined(MPI_OPT)
        USE mpi
        USE mpi_params
#endif

        IMPLICIT NONE

        INTEGER, INTENT(in) :: n1, n2
        DOUBLE PRECISION, DIMENSION(n1,n2), INTENT(inout) :: array
        INTEGER, INTENT(in) :: mystart,myend
        INTEGER :: ourstart, ourend
        INTEGER :: i

        CALL MPI_REDUCE(mystart, ourstart, 1, MPI_INTEGER, MPI_MIN, 0, comm_shar, ierr_mpi)
        CALL MPI_REDUCE(myend,     ourend, 1, MPI_INTEGER, MPI_MAX, 0, comm_shar, ierr_mpi)
        IF (shar_rank.EQ.0) THEN
            DO i = 1, ourstart-1
                array(:,i) = 0 ! Zero array "above" data to keep
            END DO
            DO i = ourend+1, n2
                array(:,i) = 0 ! Zero array "below" data to keep
            END DO
        END IF
        CALL MPI_BARRIER( comm_shar, ierr_mpi)
        IF (shar_rank.EQ.0) THEN ! Reduce arrays onto all shared memory islands
            CALL MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
        END IF
        CALL MPI_BARRIER( comm_shar, ierr_mpi)

        END SUBROUTINE mumaterial_sync_array2d_dbl

        
        
        SUBROUTINE mumaterial_syncM(array, n, outside)

#if defined(MPI_OPT)
      USE mpi
      USE mpi_params
#endif

      IMPLICIT NONE

      INTEGER, INTENT(in) :: n
      DOUBLE PRECISION, DIMENSION(3,n), INTENT(inout) :: array
      INTEGER, DIMENSION(n-domsize), INTENT(in) :: outside
      INTEGER :: i

      IF (shar_rank.EQ.0) THEN
        DO i = 1, n-domsize
          array(:,outside(i)) = 0
        END DO
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, array, 3*n, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
      END IF
      CALL MPI_BARRIER( comm_shar, ierr_mpi)

      END SUBROUTINE mumaterial_syncM

      SUBROUTINE mumaterial_getState(fx, fy, x, y)
      !-----------------------------------------------------------------------
      ! mumaterial_getState: Interpolates a function f at x to get a value y using B-splines based on De Boor's algorithm
      !-----------------------------------------------------------------------
      ! param[in]: fx. x-coordinates of function to be interpolated
      ! param[in]: fy. y-coordinates of function to be interpolated
      ! param[in]: x. x-coordinate of evaluation point
      ! param[out]: y. y-coordinate of evaluation point
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

      SUBROUTINE mumaterial_gethdipole(pos1, pos2, mag, vol, H)
            
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3), INTENT(in) :: pos1, pos2, mag
      DOUBLE PRECISION, DIMENSION(3), INTENT(inout) :: H
      DOUBLE PRECISION, DIMENSION(3) :: n, mom, r, rhat
      DOUBLE PRECISION :: rnorm, vol, pi

      pi = 4.D0*ATAN(1.D0)

      r = (pos2-pos1)
      rnorm = NORM2(r)
      rhat = r/rnorm

      mom = vol*mag
      H = H + (3*dot_product(mom, rhat)*rhat-mom)/(4*pi*rnorm**3)

      END SUBROUTINE mumaterial_gethdipole
      


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
      DOUBLE PRECISION :: H(3), N(3,3)
      INTEGER :: i

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
      DOUBLE PRECISION :: H(3), N(3,3)
      INTEGER :: i

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


      SUBROUTINE mumaterial_getb_vector(x, y, z, B, getBfld)!, linclvac)
      !-----------------------------------------------------------------------
      ! mumaterial_getb_vector: Calculates total magnetic field at multiple points in space
      !-----------------------------------------------------------------------
      ! param[in]: x. x-coordinates of points at which to determine the magnetic field
      ! param[in]: y. y-coordinates of points at which to determine the magnetic field
      ! param[in]: z. z-coordinates of points at which to determine the magnetic field
      ! param[in]: linclvac. Whether or not vacuum magnetic field should be included.
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
      DOUBLE PRECISION, ALLOCATABLE :: B_local(:,:)
      INTEGER :: mystart, myend
      INTEGER :: i 
      INTEGER :: npoints
!      LOGICAL, OPTIONAL :: linclvac

!      IF (.NOT.(PRESENT(linclvac))) linclvac = .TRUE.


      npoints = size(x)
      mystart = 1; myend = npoints

#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_CALC_MYRANGE(comm_world, 1, npoints, mystart, myend)
#endif

      allocate(B_local(3,npoints),B(3,npoints))
      B_local = 0; B = 0
      
!     IF (linclvac) THEN
      DO i = mystart, myend
        CALL mumaterial_getb_scalar(   x(i), y(i), z(i), B_local(1,i), B_local(2,i), B_local(3,i), getBfld)
      END DO
!     ELSE
!       DO i = mystart, myend
!          CALL mumaterial_getbmag_scalar(x(i), y(i), z(i), B_local(1,i), B_local(2,i), B_local(3,i))
!       END DO
!      END IF
    
#if defined(MPI_OPT)
      IF (lcomm) THEN
        CALL MPI_REDUCE(B_local,B,3*npoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm_shar,ierr_mpi)
        IF (shar_rank.EQ.0) CALL MPI_ALLREDUCE( MPI_IN_PLACE,B,3*npoints,MPI_DOUBLE_PRECISION,MPI_SUM,comm_master,ierr_mpi)
      END IF
#endif

      deallocate(B_local)
      
      RETURN
      END SUBROUTINE mumaterial_getb_vector





      SUBROUTINE mumaterial_output(path, x, y, z, getBfld)!, linclvac)
      !-----------------------------------------------------------------------
      ! mumaterial_output: Outputs B-field and points to text files
      !-----------------------------------------------------------------------
      ! param[in]: path. Path to store files in
      ! param[in]: x. x-cooridinates of points at which to determine the magnetic field
      ! param[in]: y. y-cooridinates of points at which to determine the magnetic field
      ! param[in]: z. z-cooridinates of points at which to determine the magnetic field
      ! param[in]: linclvac. Whether or not vacuum magnetic field should be included.
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
      INTEGER :: i 
      INTEGER :: npoints
      DOUBLE PRECISION, ALLOCATABLE :: B(:,:)
!      LOGICAL, OPTIONAL :: linclvac

!      IF (.NOT.(PRESENT(linclvac))) linclvac = .TRUE.

      IF (lismaster) THEN
        npoints = size(x)
        WRITE(6,*) "Outputting points"
        OPEN(13, file='./points.dat')
        DO i = 1, npoints
          WRITE(13, "(F15.7,A,F15.7,A,F15.7)") x(i), ',', y(i), ',', z(i)
        END DO
        CLOSE(13)
      END IF

      CALL mumaterial_getb_vector(x, y, z, B, getBfld)!, linclvac)
 
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
