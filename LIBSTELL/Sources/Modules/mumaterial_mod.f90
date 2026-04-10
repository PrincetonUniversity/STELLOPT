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
            DOUBLE PRECISION, PRIVATE, ALLOCATABLE :: H(:), M(:), dMdH(:)
      END TYPE stateFunctionType    

      PROCEDURE(externalFieldFunc), POINTER :: getBfld
      ABSTRACT INTERFACE
        SUBROUTINE externalFieldFunc(x,y,z,Bx,By,Bz)
          DOUBLE PRECISION, INTENT(in)  :: x,y,z
          DOUBLE PRECISION, INTENT(out) :: Bx,By,Bz
        END SUBROUTINE externalFieldFunc
      END INTERFACE

    !-------------------------------------------------------------------
    ! Shared memory windows
      DOUBLE PRECISION, POINTER, PRIVATE :: &
        vertex(:,:),  & ! Coordinates of mesh vertices
        tet_cen(:,:), & ! Coordinates of tetrahedron centroids
        tet_vol(:),   & ! Volumes of tetrahedrons
        tet_rad(:)      ! Inradii of tetrahedrons
      DOUBLE PRECISION, POINTER, PRIVATE :: &
        tet_P(:,:,:,:), & ! Rotation matrix of all tetrahedron faces
        tet_D(:,:,:),   & ! Base vectors of all tetrahedron faces
        tet_v(:,:,:,:)! Rotated vertices of all tetrahedron faces
      INTEGER, POINTER, PRIVATE :: &
        tet(:,:)        ! Vertices of tetrahedrons
      DOUBLE PRECISION, POINTER, PRIVATE :: &
        constant_mu(:),   & ! Mu along easy axis for hard magnet (or linear)
        constant_mu_o(:), & ! Mu along off axis
        M_rem(:,:)          ! Remanent magnetization
      INTEGER, POINTER, PRIVATE :: &
        state_dex(:), & ! Index of state type
        state_type(:)   ! Type of state type (see below)
      INTEGER, PRIVATE :: &
        win_vertex, win_tet, win_tet_cen, win_tet_vol, &
        win_tet_P, win_tet_D, win_tet_v, &
        win_tet_rad, win_state_dex, win_state_type,    &
        win_constant_mu, win_constant_mu_o, win_M_rem   ! MPI windows
    !-------------------------------------------------------------------
    ! Mesh info
      DOUBLE PRECISION, PRIVATE :: &
        tet_vol_tot ! Total volume of tetrahedrons
      INTEGER, PRIVATE  :: &
        ntet, & ! Number of tetrahedrons
        nvertex ! Number of vertices
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, PRIVATE :: &
        max_tet_rad, & ! Largest tetrahedron inradius in a cluster
        vol_proc
    !-------------------------------------------------------------------
    ! Mesh division
      INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: &
        dom_sizes,& ! Sizes of each dom_proc across MPI ranks/clusters
        dom_shar, & ! Domain of elements treated by an MPI node
        dom_proc, & ! Local domain of elements treated by an MPI rank
        dom_partial_proc ! Domain of elements necessary for calculations for an MPI rank
      INTEGER, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: &
        dom_cluster ! Array of all dom_procs
      INTEGER, PRIVATE :: &
        ntet_shar, & ! Number of elements in dom_shar
        ntet_proc, & ! Number of elements in dom_proc 
        ntet_partial_proc, & ! Number of elements in dom_partial_proc
        ntet_proc_max! Largest value in dom_sizes
      INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE:: &
        map_g2all, &! Map of global indices to "all" indices 
        map_g2p     ! Map of global indices to "partial" indices
    !-------------------------------------------------------------------
    ! Magnetic material types
      TYPE(stateFunctionType), PRIVATE, ALLOCATABLE :: &
        statefunc(:) ! State function structure. 
      INTEGER, PARAMETER, PRIVATE :: & 
        TYPE_HARD   = 1, &  ! Hard magnet with remanent magnetization
        TYPE_SOFT   = 2, &  ! Soft magnet with state function
        TYPE_LINEAR = 3     ! Linear material with constant mu
      INTEGER, PRIVATE :: &
        nstate, & ! Number of unique state functions
        nlinear, ntet_linear_max, &
        nsoft,   ntet_soft_max, &
        nhard,   ntet_hard_max
      INTEGER, ALLOCATABLE, PRIVATE :: &
        dom_linear(:,:), &
        dom_soft(:,:), &
        dom_hard(:,:), &
        ntet_linear(:), &
        ntet_soft(:), &
        ntet_hard(:), &
        sdex_linear(:), &
        sdex_soft(:), &
        sdex_hard(:)
    !-------------------------------------------------------------------
    ! Magnetizations
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: &
        M_local,   & ! Magnetizations of local elements
        M_global,  & ! Magnetizations of all elements (ordered by element)
        M_partial, & ! Necessary magnetizations for calculations
        M_snapshot,& ! Snapshot of M_partial at certain iteration
        M_all        ! Magnetizations of all elements (ordered by rank)
    !-------------------------------------------------------------------
    ! Fields
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, PRIVATE  :: &
        H_ext, & ! Field from non-neighbors
        H_mid, & ! Field from dipoles
        H_app, & ! Field from currents
        H_prev
    !-------------------------------------------------------------------
    ! Neighbors
      INTEGER, ALLOCATABLE, PRIVATE :: &
        nbrs_proc(:,:),& ! Array of neighbors of local elements
        nbrs_count(:)    ! Number of neighbors for each local element
      DOUBLE PRECISION, ALLOCATABLE, PRIVATE :: &
        N_store(:,:,:,:) ! Demagnetization tensor for neighbors of element
    !-------------------------------------------------------------------
    ! Clusters
      DOUBLE PRECISION, ALLOCATABLE, PRIVATE :: &
        r_cluster(:,:), & ! Position of cluster centroid
        d_cluster(:),   & ! RMS distance of elements to cluster centroid
        m_cluster(:,:), & ! Magnetic moment vector of cluster
        Q_cluster(:,:), & ! Quadrupole tensor of cluster
        tet_cen_proc(:,:)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: &
        R_quad, mom_quad ! Helpers 
      LOGICAL, DIMENSION(:), ALLOCATABLE, PRIVATE :: iscluster_proc
    !-------------------------------------------------------------------
    ! Block solver
      DOUBLE PRECISION, ALLOCATABLE, PRIVATE :: &
        M_block(:), &
        H_block(:), &
        J_block(:,:), &
        N_block(:,:,:,:), &
        R_block(:), &
        M_spline(:,:), &
        dMdH_spline(:), &
        M_norm_proc(:), &
        H_norm_proc(:), &
        dMdH_mat(:,:,:)
    !-------------------------------------------------------------------
    ! Dipoles
      LOGICAL, DIMENSION(:), ALLOCATABLE, PRIVATE :: &
        is_stale      ! Helper mask to determine which dipoles to recompute
      INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: &
        dom_dip, &    ! Helper array to store dipole indices
        stale_dex     ! Indices of stale elements 
      DOUBLE PRECISION, ALLOCATABLE, PRIVATE :: &
        M_dip(:,:), & ! Contiguous helper array of dipole magnetizations
        r_partial(:,:), & ! Contiguous helper array of dipole positions 
        V_partial(:), &   ! Contiguous helper array of dipole volumes
        r_dip(:,:), & ! Contiguous helper array of dipole positions 
        V_dip(:) ! Contiguous helper array of dipole volumes
    !-------------------------------------------------------------------
    ! User settings
      DOUBLE PRECISION, PRIVATE :: &
        threshold   =1.0D-3, & ! Convergence threshold of dW/W_cl
        pF   =20.0D0, & ! Neighbor cutoff
        lambdaStart =0.10D0, & ! Starting lambda value
        lambdaFactor=0.75D0, & ! Multiplicative lambda value
        convCheck=99.9D0      ! Converged at %. No longer used
      INTEGER, PRIVATE :: &
        lambdaThresh=1, & ! No longer used
        maxIter=10000     ! Maximum number of iterations
      DOUBLE PRECISION, PARAMETER, PRIVATE :: &
        cutoff = 3.0d0    ! Cluster cutoff
                          ! Should be made a user variable
    !-------------------------------------------------------------------
    ! MPI variables
      LOGICAL, PRIVATE :: &
        lcomm, &   ! True if running with MPI
        lismaster  ! True if world_rank = 0
      INTEGER, PRIVATE :: &
        comm_world, comm_shar, comm_master, & ! Global, share, and master communicators
        world_size, shar_size, master_size, & ! Global, share, and master sizes
        world_rank, shar_rank, master_rank    ! Global, share, and master ranks
      INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: &
        shar_ranks_arr ! Array of which MPI ranks are in same comm_shar
      INTEGER, PRIVATE :: &
        mstat(MPI_STATUS_SIZE) ! Status for MPI_RECV
    !-------------------------------------------------------------------
    ! Synchronization of magnetizations
      INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: &
        sync_counts, & ! Number of doubles an MPI rank sends in sync
        sync_displs    ! Offset of MPI rank send in sync
      INTEGER, PARAMETER, PRIVATE :: &
        SYNC_READ = 1, & ! Initial broadcast after read
        SYNC_ITER = 2, & ! Standard sync during iterations
        SYNC_DONE = 3    ! Finalize magnetization arrays

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eig_proc
    !-------------------------------------------------------------------
    ! I/O
      LOGICAL, PRIVATE  :: & 
        lverb ! Verbose flag
      CHARACTER(LEN=256), PRIVATE :: &
        file_string,    & ! Mesh file to read in
        machine_string, & ! Mesh name as read in
        date_string       ! Date as read in
      LOGICAL, DIMENSION(:), ALLOCATABLE, PRIVATE :: &
        mask_out ! Helper mask
      INTEGER, DIMENSION(:), ALLOCATABLE, PRIVATE :: &
        nbrs_out ! Elements close to evaluation point (follows neighbor logic)
    !------------------------------------------------------------------------------
    ! Constants
      DOUBLE PRECISION, PARAMETER, PRIVATE :: &
        PI     = 4.0D0*ATAN(1.0d0),& ! pi
        INVPI  = 1.0D0/PI,         & ! 1/pi
        INV4PI = 1.0D0/(4.0D0*PI), & ! 1/(4pi)
        INV8PI = 1.0D0/(8.0D0*PI), & ! 1/(8pi)
        MU0 = 4.0D-7*PI,           & ! Permeability of free space
        small = 1.0D-12,           & ! A small number
        smaller = 1.0D-100,        & ! A smaller number
        I3(3,3) = RESHAPE([1.0d0, 0.0d0, 0.0d0, &
                           0.0d0, 1.0d0, 0.0d0, &
                           0.0d0, 0.0d0, 1.0d0],[3,3]) ! identity matrix
    !------------------------------------------------------------------------------
    ! Functions
      PRIVATE :: CROSS_PRODUCT, TET_VOLUME, INVERT_3x3, OUTER_PRODUCT, &
                 GET_DEMAG, GET_Nxz, GET_Nyz, GET_Nzz, SOLVE_3x3

!------------------------------------------------------------------------------
!     Subroutines
!       Main flow
!         mumaterial_set_comms:     Sets up MPI communicators (optional)
!         mumaterial_load:      Loads magnetic material file and sets up MPI stuff
!         mumaterial_set_user:   Sets default values
!         mumaterial_set_verb:   Sets standard verbosity
!         mumaterial_set_Bfld:   Sets function for external B field
!         mumaterial_info:      Prints information to screen
!         mumaterial_run:      Initializes everything, calls iteration subroutine
!         mumaterial_iterate: Main calculation loop
!
!       Helpers
!         TET_VOLUME:  Calculates volume of a tetrahedron
!         mumaterial_init_neighbors: Determines tetrahedron neighbors
!         GET_DEMAG:          Determines demagnetization tensor
!           GET_Nxz: x-component 
!           GET_Nyz: y-component 
!           GET_Nzz: z-component
!         CROSS_PRODUCT:         Cross product of two vectors
!         mumaterial_getstate_scalar:      Interpolates function 
!
!       MPI 
!         mumaterial_split_sub:            Splits input domain into two subdomains
!         mumaterial_syncmag: Syncs magnetization array on shar_mem nodes
!         mumaterial_free:  Frees MPI memory
!       Output
!         mumaterial_output:  Output B-field and points to file
!         mumaterial_getb:    Calculates magnetic field in space
!             mumaterial_getb_scalar:      Single point in space
!             mumaterial_getb_vector: Multiple points in space
!     Functions
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
! mumaterial_abort: Helper function for code stop. Works with(out) MPI.
!------------------------------------------------------------------------------     
      SUBROUTINE mumaterial_abort()
      IMPLICIT NONE
      IF (lcomm) THEN 
#if defined(MPI_OPT)
          CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr_mpi)
#endif
      ELSE 
        STOP
      ENDIF
      END SUBROUTINE mumaterial_abort

!-----------------------------------------------------------------------
! mumaterial_alloc_init: Allocates M_local, M_partial, etc
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_alloc_init()

      IMPLICIT NONE
      ! All ranks need M_local and M_partial
      ALLOCATE(M_local(3,ntet_proc),M_partial(3,ntet_partial_proc),M_snapshot(3,ntet_partial_proc),M_all(3,ntet))
      M_local = 0.0d0
      M_partial = 0.0d0
      M_snapshot = 0.0d0
      M_all = 0.0d0
      ALLOCATE(M_norm_proc(ntet_proc),H_norm_proc(ntet_proc),dMdH_mat(3,3,ntet_proc))
      ! Solver
      ALLOCATE(M_block(3*ntet_proc),H_block(3*ntet_proc),R_block(3*ntet_proc),J_block(3*ntet_proc,3*ntet_proc),M_spline(3,ntet_proc),dMdH_spline(ntet_proc))
      ! Allocate helpers
      ALLOCATE(R_quad(3,MAXVAL(dom_sizes)),mom_quad(3,MAXVAL(dom_sizes)))
      ALLOCATE(H_ext(3,ntet_proc),H_mid(3,ntet_proc),H_prev(3,ntet_proc))
      ! Intermediate field dipoles
      ALLOCATE(stale_dex(ntet_partial_proc))
      ALLOCATE(is_stale(ntet_partial_proc))
      ALLOCATE(r_dip(3,ntet_partial_proc),V_dip(ntet_partial_proc),M_dip(3,ntet_partial_proc))

      END SUBROUTINE mumaterial_alloc_init

!-----------------------------------------------------------------------
! mumaterial_dealloc_init: Deallocates arrays after iterations
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_dealloc_init()
      IMPLICIT NONE
      INTEGER :: ik

      ! Solver
      DEALLOCATE(M_block,H_block,R_block,J_block,M_spline,dMdH_spline)
      DEALLOCATE(M_norm_proc,H_norm_proc,dMdH_mat)
      ! Fields
      DEALLOCATE(H_app,H_mid,H_ext,H_prev)
      ! Domains
      DEALLOCATE(dom_proc,dom_partial_proc,tet_cen_proc)
      ! Neighbors
      DEALLOCATE(nbrs_proc,nbrs_count,N_store)
      ! Dipoles
      DEALLOCATE(r_partial,V_partial,is_stale,stale_dex)
      DEALLOCATE(r_dip,V_dip,M_dip)
      ! Cluster
      DEALLOCATE(shar_ranks_arr)
      ! Magnetizations (all ranks)
      DEALLOCATE(M_partial,M_local,M_snapshot)
      ! Sync
      DEALLOCATE(sync_counts,sync_displs)
      ! State functions
      DO ik = 1, nstate
        IF (state_type(ik).EQ.TYPE_SOFT) THEN
          DEALLOCATE(statefunc(ik)%H,statefunc(ik)%M,statefunc(ik)%dMdH)
        END IF
      END DO
      DEALLOCATE(statefunc)
      DEALLOCATE(dom_linear,dom_soft,dom_hard)

      RETURN
      END SUBROUTINE mumaterial_dealloc_init

!-----------------------------------------------------------------------
! mumaterial_alloc_out: Allocates helper arrays for output
!                       (Avoids repeated alloc/dealloc)
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_alloc_out()
      IMPLICIT NONE

      ! Neighbors
      ALLOCATE(mask_out(ntet), nbrs_out(ntet))
      ! Dipoles
      ALLOCATE(r_dip(3,ntet),M_dip(3,ntet),V_dip(ntet),dom_dip(ntet))

      RETURN

      END SUBROUTINE mumaterial_alloc_out

!-----------------------------------------------------------------------
! mumaterial_dealloc_out: Dellocates helper arrays after output
!                         (Avoids repeated alloc/dealloc)
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_dealloc_out()
      IMPLICIT NONE

      DEALLOCATE(M_global)
      ! Neighbors
      DEALLOCATE(mask_out, nbrs_out)
      ! Clusters
      DEALLOCATE(iscluster_proc, dom_cluster, dom_sizes, r_cluster, d_cluster)
      DEALLOCATE(map_g2p,m_cluster, Q_cluster)
      DEALLOCATE(R_quad,mom_quad)
      DEALLOCATE(max_tet_rad)
      ! Dipoles
      DEALLOCATE(r_dip,M_dip,V_dip,dom_dip)
      RETURN

      END SUBROUTINE mumaterial_dealloc_out
  

!------------------------------------------------------------------------------
! mumaterial_free: Destroys MPI windows
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_free()
      IMPLICIT NONE

      IF (ASSOCIATED(state_dex))     CALL mpidealloc(state_dex,win_state_dex)
      IF (ASSOCIATED(state_type))    CALL mpidealloc(state_type,win_state_type)
      IF (ASSOCIATED(constant_mu))   CALL mpidealloc(constant_mu,win_constant_mu)
      IF (ASSOCIATED(constant_mu_o)) CALL mpidealloc(constant_mu_o,win_constant_mu_o)
      IF (ASSOCIATED(tet))           CALL mpidealloc(tet,win_tet)
      IF (ASSOCIATED(vertex))        CALL mpidealloc(vertex,win_vertex)
      IF (ASSOCIATED(tet_cen))       CALL mpidealloc(tet_cen,win_tet_cen)
      IF (ASSOCIATED(tet_vol))       CALL mpidealloc(tet_vol,win_tet_vol)
      IF (ASSOCIATED(tet_rad))       CALL mpidealloc(tet_rad,win_tet_rad)
      IF (ASSOCIATED(M_rem))         CALL mpidealloc(M_rem,win_M_rem)
      IF (ASSOCIATED(tet_P))         CALL mpidealloc(tet_P,win_tet_P)
      IF (ASSOCIATED(tet_D))         CALL mpidealloc(tet_D,win_tet_D)
      IF (ASSOCIATED(tet_v))         CALL mpidealloc(tet_v,win_tet_v)

      RETURN
      END SUBROUTINE mumaterial_free

!------------------------------------------------------------------------------
! mumaterial_set_comms: Sets up communicators for mumaterial
!------------------------------------------------------------------------------
! param[in]:  comm. World communicator from which other comms are born
! param[out]: comm_shar_out. Shared memory communicator for calculations
! param[out]: comm_master_out. Master communicator handles cross-node stuff
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_set_comms(comm, comm_shar_out, comm_master_out)

      IMPLICIT NONE

      INTEGER, INTENT(inout) :: comm
      INTEGER, INTENT(out) :: comm_shar_out, comm_master_out
      INTEGER :: comm_myworld, color

      CALL MPI_COMM_DUP( comm, comm_myworld, ierr_mpi )
      CALL MPI_COMM_SPLIT_TYPE( comm_myworld, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, comm_shar_out, ierr_mpi)
      CALL MPI_COMM_RANK( comm_shar_out, shar_rank, ierr_mpi)

      color = MPI_UNDEFINED
      IF (shar_rank.EQ.master) color = 0
      CALL MPI_COMM_SPLIT( comm_myworld, color, shar_rank, comm_master_out, ierr_mpi )

      RETURN

      END SUBROUTINE mumaterial_set_comms

!------------------------------------------------------------------------------
! mumaterial_set_user: Sets default values
!------------------------------------------------------------------------------
! param[in]: mE. threshold: threshold for determining convergence
! param[in]: mI. maxIter: max amount of iterations
! param[in]: la. lambdaStart: initial value of lambda_n
! param[in]: laF. lambdaFactor: multiplicative factor for lambda_n
! param[in]: laT. lambdaThresh: amount of dM>0 before lambda_n is multiplied
! param[in]: padF. pF: factor for sphere around tets for neighbors
! param[in]: cc. convCheck: Stop when this percentage of elemnts has converged
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_set_user(mE, mI, la, laF, laT, padF, cc)
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in) :: mE, la, laF, padF, cc
      INTEGER, INTENT(in) :: mI, laT

      threshold = mE
      maxIter = mI
      lambdaStart = la
      lambdaFactor = laF
      lambdaThresh = laT
      pF = padF
      convCheck = cc

      RETURN

      END SUBROUTINE mumaterial_set_user

!------------------------------------------------------------------------------
! mumaterial_set_verb: Sets Verbosity
!------------------------------------------------------------------------------
! param[in]: lverbin. Verbosity on
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_set_verb(lverbin)
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: lverbin
      lverb = lverbin
      RETURN
      END SUBROUTINE mumaterial_set_verb
  
!------------------------------------------------------------------------------
! mumaterial_set_Bfld: associates external B field function getBfld with func_B
!------------------------------------------------------------------------------
! param[in]: func_B. External B field function.
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_set_Bfld(func_B)
      IMPLICIT NONE
      PROCEDURE(externalFieldFunc) :: func_B
      getBfld => func_B
      RETURN 
      END SUBROUTINE mumaterial_set_Bfld

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
      INTEGER :: iunit ,ik, i, j, nMH, stype, n, ilinear, ihard, isoft

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
      IF (lcomm) THEN
#if defined(MPI_OPT)
        lismaster = .FALSE.; master_rank = 1
        CALL MPI_COMM_RANK( comm_world, world_rank, ierr_mpi)
        CALL MPI_COMM_SIZE( comm_world, world_size, ierr_mpi)
        CALL MPI_COMM_RANK( comm_shar,  shar_rank,  ierr_mpi)
        CALL MPI_COMM_SIZE( comm_shar,  shar_size,  ierr_mpi)
        IF (shar_rank.eq.master) THEN
          CALL MPI_COMM_RANK( comm_master, master_rank, ierr_mpi )
          CALL MPI_COMM_SIZE( comm_master, master_size, ierr_mpi )
          lismaster = (master_rank.EQ.0)
          ALLOCATE(shar_ranks_arr(shar_size))
        ELSE          
          ALLOCATE(shar_ranks_arr(1))
        END IF
        CALL MPI_GATHER(world_rank, 1, MPI_INTEGER, shar_ranks_arr, 1, MPI_INTEGER, master, comm_shar, ierr_mpi)
        CALL MPI_Bcast( master_size, 1, MPI_INTEGER, 0, comm_shar, ierr_mpi)
#endif
      END IF

      ! Nullify pointers
      NULLIFY(vertex, tet, tet_cen, tet_vol, tet_rad, state_dex, state_type, &
              constant_mu, constant_mu_o, M_rem)

      ! open file, return if fails
      iunit = 327; istat = 0
      file_string = TRIM(filename)
      CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
      IF (istat/= 0) RETURN
      ! master reads info
      IF (lismaster) THEN
        READ(iunit,'(A)') machine_string
        READ(iunit,'(A)') date_string
        READ(iunit,*) nvertex, ntet, nstate
      END IF

      ! Broadcast info to MPI and allocate vertex and face info
      IF (lcomm) THEN
#if defined(MPI_OPT)
        CALL MPI_Bcast(nvertex,1,MPI_INTEGER,0,comm_world,ierr_mpi)
        CALL MPI_Bcast(ntet,   1,MPI_INTEGER,0,comm_world,ierr_mpi)
        CALL MPI_Bcast(nstate, 1,MPI_INTEGER,0,comm_world,ierr_mpi)

        CALL mpialloc(vertex,3,nvertex,    shar_rank,0,comm_shar,win_vertex)
        CALL mpialloc(tet,4,ntet,          shar_rank,0,comm_shar,win_tet)
        CALL mpialloc(tet_cen,3,ntet,      shar_rank,0,comm_shar,win_tet_cen)
        CALL mpialloc(tet_vol,ntet,        shar_rank,0,comm_shar,win_tet_vol)
        CALL mpialloc(tet_rad,ntet,        shar_rank,0,comm_shar,win_tet_rad)
        CALL mpialloc(state_dex,ntet,      shar_rank,0,comm_shar,win_state_dex)
        CALL mpialloc(state_type,nstate,   shar_rank,0,comm_shar,win_state_type)
        CALL mpialloc(constant_mu,nstate,  shar_rank,0,comm_shar,win_constant_mu)
        CALL mpialloc(constant_mu_o,nstate,shar_rank,0,comm_shar,win_constant_mu_o)
        CALL mpialloc(M_rem,3,nstate,      shar_rank,0,comm_shar,win_M_rem)
        CALL mpialloc(tet_P,3,3,4,ntet,    shar_rank,0,comm_shar,win_tet_P)
        CALL mpialloc(tet_D,3,3,ntet,      shar_rank,0,comm_shar,win_tet_D)
        CALL mpialloc(tet_v,3,3,4,ntet,    shar_rank,0,comm_shar,win_tet_v)

#endif
      ELSE! if no MPI, allocate everything on one node
          ALLOCATE(vertex(3,nvertex))
          ALLOCATE(tet(4,ntet))
          ALLOCATE(tet_cen(3,ntet))
          ALLOCATE(tet_vol(ntet))
          ALLOCATE(tet_rad(ntet))
          ALLOCATE(state_dex(ntet))
          ALLOCATE(state_type(nstate))
          ALLOCATE(constant_mu(nstate))
          ALLOCATE(constant_mu_o(nstate))
          ALLOCATE(M_rem(3,nstate))
          ALLOCATE(tet_P(3,3,4,ntet))
          ALLOCATE(tet_D(3,4,ntet))
          ALLOCATE(tet_v(3,3,4,ntet))
      END IF
      ALLOCATE(statefunc(nstate))

      IF (lismaster) THEN
        ! Get mesh
        DO ik = 1, nvertex
          READ(iunit,*) vertex(1,ik),vertex(2,ik),vertex(3,ik)
        END DO
        DO ik = 1, ntet
          READ(iunit,*) tet(1,ik),tet(2,ik),tet(3,ik),tet(4,ik),state_dex(ik)
        END DO

        ! Get state functions
        nlinear = 0
        nsoft = 0
        nhard = 0
        DO ik = 1, nstate
          READ(iunit,*) state_type(ik)
          !------------------------------------ 
          ! Hard magnet with remanent magnetization
          IF     (state_type(ik) .EQ. TYPE_HARD) THEN
            READ(iunit,*) constant_mu(ik), constant_mu_o(ik)
            READ(iunit,*) M_rem(1,ik), M_rem(2,ik), M_rem(3,ik)
            nhard = nhard + 1
          !------------------------------------ 
          ! Soft magnet
          ELSEIF (state_type(ik) .EQ. TYPE_SOFT) THEN
            READ(iunit,*) nMH
            ALLOCATE(statefunc(ik)%H(nMH),statefunc(ik)%M(nMH),statefunc(ik)%dMdH(nMH))
            READ(iunit,*) statefunc(ik)%H(:)
            READ(iunit,*) statefunc(ik)%M(:)
            CALL mumaterial_getstate_slopes(statefunc(ik)%H,statefunc(ik)%M,statefunc(ik)%dMdH)
            nsoft = nsoft + 1
          !------------------------------------ 
          ! Linear
          ELSEIF (state_type(ik) .EQ. TYPE_LINEAR) THEN 
            READ(iunit,*) constant_mu(ik)
            nlinear = nlinear + 1
          ELSE
            PRINT *, '!!! UNKNOWN STATE_TYPE == ',state_type(ik)
          END IF
        END DO
        ! Get number of everything
        ALLOCATE(sdex_linear(nlinear), sdex_soft(nsoft), sdex_hard(nhard))
        sdex_linear = 0
        sdex_soft = 0
        sdex_hard = 0

        ilinear = 0
        isoft = 0
        ihard = 0
        DO ik = 1, nstate
          n = COUNT(state_dex==ik)
          stype = state_type(ik)
          SELECT CASE (stype)
            CASE (TYPE_LINEAR)
              ilinear = ilinear + 1
              sdex_linear(ilinear) = n
            CASE (TYPE_SOFT)
              isoft = isoft + 1
              sdex_soft(isoft) = n
            CASE (TYPE_HARD) 
              ihard = ihard  + 1
              sdex_hard(ihard) = n
          END SELECT
        END DO
        ntet_linear_max = MAXVAL(sdex_linear)
        ntet_soft_max = MAXVAL(sdex_soft)
        ntet_hard_max = MAXVAL(sdex_hard)
        DEALLOCATE(sdex_linear,sdex_soft,sdex_hard)
      END IF
      ! Close file
      CLOSE(iunit)

      IF (lcomm) THEN
#if defined(MPI_OPT)
      ! BCast from world master to shar masters
      IF (shar_rank.EQ.master) THEN
        CALL MPI_Bcast(vertex,       3*nvertex,MPI_DOUBLE_PRECISION,0,comm_master,ierr_mpi)
        CALL MPI_Bcast(tet,          4*ntet,   MPI_INTEGER,         0,comm_master,ierr_mpi)
        CALL MPI_Bcast(state_dex,    ntet,     MPI_INTEGER,         0,comm_master,ierr_mpi)
        CALL MPI_Bcast(state_type,   nstate,   MPI_INTEGER,         0,comm_master,ierr_mpi)
        CALL MPI_Bcast(constant_mu,  nstate,   MPI_DOUBLE_PRECISION,0,comm_master,ierr_mpi)
        CALL MPI_Bcast(constant_mu_o,nstate,   MPI_DOUBLE_PRECISION,0,comm_master,ierr_mpi)
        CALL MPI_Bcast(M_rem,      3*nstate,   MPI_DOUBLE_PRECISION,0,comm_master,ierr_mpi) 
      END IF
      CALL MPI_BARRIER(comm_world, ierr_mpi)
      ! Broadcast state functions to non-masters
      DO ik = 1, nstate
        IF (state_type(ik).EQ.TYPE_SOFT) THEN
          IF (lismaster) nMH = SIZE(statefunc(ik)%H)
          CALL MPI_Bcast(nMH,1,MPI_INTEGER,0,comm_world,ierr_mpi)
          IF (world_rank.NE.master) ALLOCATE(statefunc(ik)%H(nMH),statefunc(ik)%M(nMH),statefunc(ik)%dMdH(nMH))
          CALL MPI_Bcast(statefunc(ik)%H,   nMH,MPI_DOUBLE_PRECISION,0,comm_world,ierr_mpi)
          CALL MPI_Bcast(statefunc(ik)%M,   nMH,MPI_DOUBLE_PRECISION,0,comm_world,ierr_mpi)
          CALL MPI_Bcast(statefunc(ik)%dMdH,nMH,MPI_DOUBLE_PRECISION,0,comm_world,ierr_mpi)
        END IF
      END DO
      CALL MPI_BCAST(nlinear, 1, MPI_INTEGER, 0, comm_world, ierr_mpi)
      CALL MPI_BCAST(nsoft,   1, MPI_INTEGER, 0, comm_world, ierr_mpi)
      CALL MPI_BCAST(nhard,   1, MPI_INTEGER, 0, comm_world, ierr_mpi)
      CALL MPI_BCAST(ntet_linear_max, 1, MPI_INTEGER, 0, comm_world, ierr_mpi)
      CALL MPI_BCAST(ntet_soft_max,   1, MPI_INTEGER, 0, comm_world, ierr_mpi)
      CALL MPI_BCAST(ntet_hard_max,   1, MPI_INTEGER, 0, comm_world, ierr_mpi)
#endif
      END IF
      ! set default values
      CALL mumaterial_set_user(1.0d-5, 100, 0.7d0, 0.75d0, 10, 20.d0, 99.d0)
      RETURN

      END SUBROUTINE mumaterial_load

!-----------------------------------------------------------------------
! mumaterial_init_mesh: Offsets tetrahedrons, calculates centroid, 
! volume, inradius, and synchronizes across all ranks.
!-----------------------------------------------------------------------
! param[in]: offset(3)
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_init_mesh(offset)
      IMPLICIT NONE 

      DOUBLE PRECISION, INTENT(in), OPTIONAL :: offset(3)
      INTEGER :: mystart, myend, i
      !-----------------------------------------
      ! Offset vertices and synchronize
      !-----------------------------------------
      IF (PRESENT(offset).AND.(NORM2(offset) .GT. 0.d0)) THEN
        mystart = 1; myend = nvertex
#if defined(MPI_OPT)
        IF (lcomm) CALL MPI_CALC_MYRANGE(comm_world, 1, nvertex, mystart, myend)
#endif
        DO i = mystart, myend
          vertex(:,i) = vertex(:,i) + offset
        END DO
            
#if defined(MPI_OPT)
        CALL MPI_BARRIER(comm_world, ierr_mpi)
        IF (shar_rank.EQ.master) THEN
          CALL MPI_ALLREDUCE( MPI_IN_PLACE, vertex,   3*nvertex, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
        END IF
#endif
      END IF
      !-----------------------------------------
      ! Calculate center, volume, inradius
      !-----------------------------------------
      ! Master wipes shared arrays
      IF (shar_rank.EQ.master) THEN 
        tet_cen = 0; tet_vol = 0; tet_rad = 0     
      END IF  
      ! Calculate range
      mystart = 1; myend = ntet
      IF (lcomm) THEN 
#if defined(MPI_OPT)
        CALL MPI_BARRIER(comm_world, ierr_mpi)
        CALL MPI_CALC_MYRANGE(comm_world, 1, ntet, mystart, myend) 
#endif
      END IF

      DO i = mystart, myend
          tet_cen(:,i) = (vertex(:,tet(1,i)) + vertex(:,tet(2,i)) + &
                          vertex(:,tet(3,i)) + vertex(:,tet(4,i))) / 4.d0
          tet_vol(i) = TET_VOLUME( &
              vertex(:,tet(1,i)),vertex(:,tet(2,i)), vertex(:,tet(3,i)),vertex(:,tet(4,i)))
          tet_rad(i) = SQRT(6.0)/12.d0*(6.d0*SQRT(2.0d0)*tet_vol(i))**(1.0d0/3.0d0)
      END DO
      !-----------------------------------------
      ! Synchronize shared arrays
      !-----------------------------------------
#if defined(MPI_OPT)
      CALL MPI_BARRIER(comm_shar, ierr_mpi)
      IF (shar_rank.EQ.master) THEN
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tet_cen, 3*ntet, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tet_vol,   ntet, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tet_rad,   ntet, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
      ENDIF
      CALL MPI_BARRIER(comm_shar, ierr_mpi)
#endif
      tet_vol_tot = SUM(tet_vol)

      END SUBROUTINE mumaterial_init_mesh
!-----------------------------------------------------------------------
! mumaterial_init_states: Constructs dom_linear, dom_soft, dom_hard, which
! groups mesh elements by their state function. For every dom(i,j) array,
! j corresponds to the state function, and i are elements.
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_init_states()

      IMPLICIT NONE
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: dom_linear_temp, dom_soft_temp, dom_hard_temp
      INTEGER :: itet, j, dex, sdex, stype, free_sdex_lin, free_sdex_soft, free_sdex_hard

      ALLOCATE(dom_linear_temp(ntet_linear_max,nlinear), dom_soft_temp(ntet_soft_max,nsoft), dom_hard_temp(ntet_hard_max,nhard))
      ALLOCATE(sdex_linear(nlinear), sdex_soft(nsoft), sdex_hard(nhard))
      ALLOCATE(ntet_linear(nlinear), ntet_soft(nsoft), ntet_hard(nhard))
      dom_linear_temp = 0;   dom_soft_temp = 0;  dom_hard_temp = 0 ! Element indices
      ntet_linear = 0; ntet_soft = 0; ntet_hard = 0             ! Number of elements per state function
      sdex_linear = 0; sdex_soft = 0; sdex_hard = 0             ! State function index
      free_sdex_lin = 1; free_sdex_soft = 1; free_sdex_hard = 1 ! Next free state function index
      
      DO itet = 1, ntet_proc
        dex = 0
        sdex = state_dex(itet)
        stype = state_type(sdex)
        SELECT CASE (stype)
          ! Linear
          CASE (TYPE_LINEAR)
            DO j = 1, free_sdex_lin-1
              IF (sdex_linear(j) == sdex) THEN
                dex = j
                EXIT
              END IF
            END DO
            IF (dex == 0) THEN
              sdex_linear(free_sdex_lin) = sdex
              dex = free_sdex_lin
              free_sdex_lin = free_sdex_lin + 1
            END IF
            ntet_linear(dex) = ntet_linear(dex) + 1
            dom_linear_temp(ntet_linear(dex),dex) = itet
          ! Soft
          CASE (TYPE_SOFT)
            DO j = 1, free_sdex_soft-1
              IF (sdex_soft(j) == sdex) THEN
                dex = j
                EXIT
              END IF
            END DO
            IF (dex == 0) THEN
              sdex_soft(free_sdex_soft) = sdex
              dex = free_sdex_soft
              free_sdex_soft = free_sdex_soft + 1
            END IF
            ntet_soft(dex) = ntet_soft(dex)+1
            dom_soft_temp(ntet_soft(dex),dex) = itet
          ! Hard
          CASE (TYPE_HARD)
            DO j = 1, free_sdex_hard-1
              IF (sdex_hard(j) == sdex) THEN
                dex = j
                EXIT
              END IF
            END DO
            IF (dex == 0) THEN
              sdex_hard(free_sdex_hard) = sdex
              dex = free_sdex_hard
              free_sdex_hard = free_sdex_hard + 1
            END IF
            ntet_hard(dex) = ntet_hard(dex)+1
            dom_hard_temp(ntet_hard(dex),dex) = itet

        END SELECT
      END DO
      ntet_linear_max  = MAXVAL(ntet_linear)
      ntet_soft_max = MAXVAL(ntet_soft)
      ntet_hard_max = MAXVAL(ntet_hard) 

      ALLOCATE(dom_linear(ntet_linear_max,nlinear), &
               dom_soft(ntet_soft_max,nsoft), &
               dom_hard(ntet_hard_max,nhard))
      dom_linear = dom_linear_temp(1:ntet_linear_max, :)
      dom_soft   = dom_soft_temp(1:ntet_soft_max, :)
      dom_hard   = dom_hard_temp(1:ntet_hard_max, :)
      DEALLOCATE(dom_linear_temp,dom_hard_temp,dom_soft_temp)

      END SUBROUTINE mumaterial_init_states
!-----------------------------------------------------------------------
! mumaterial_init_happ: Calculates H_app for every local element
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_init_happ()

      IMPLICIT NONE

      INTEGER :: i, i_tile
      DOUBLE PRECISION :: H_app_norm_min, H_app_norm_max
      DOUBLE PRECISION :: Bx, By, Bz

      !---------------------------------
      ! Get field from external function
      !---------------------------------
      ALLOCATE(H_app(3,ntet_proc))
      H_app = 0.0d0
      DO i = 1, ntet_proc
        CALL getBfld(tet_cen_proc(1,i), tet_cen_proc(2,i), tet_cen_proc(3,i), Bx, By, Bz)
        H_app(:,i) = [Bx/MU0, By/MU0, Bz/MU0]
      END DO
      !---------------------------------
      ! Display only
      !---------------------------------
      H_app_norm_min = MINVAL(NORM2(H_app,DIM=1))
      H_app_norm_max = MAXVAL(NORM2(H_app,DIM=1))
      IF (lcomm) THEN
#if defined(MPI_OPT)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,H_app_norm_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm_world,ierr_mpi)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,H_app_norm_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm_world,ierr_mpi)
#endif
      END IF

      IF (lverb) THEN
        WRITE(6,'(3X,A,ES11.2,A,/,17X,ES11.2,A)') 'H-field range:',H_app_norm_min,' A/m,',H_app_norm_max,' A/m'
        FLUSH(6)
      END IF
      END SUBROUTINE mumaterial_init_happ
      
!-----------------------------------------------------------------------
! mumaterial_init_demag: Calculates and synchronizes helpers for demagnetization tensor calculation
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_init_demag()

      IMPLICIT NONE
            
      INTEGER :: i, i_tile

#if defined(MPI_OPT)
      IF (shar_rank.EQ.master) THEN
        tet_P = 0.0d0
        tet_D = 0.0d0
        tet_v = 0.0d0
      END IF
      CALL MPI_BARRIER(comm_world, ierr_mpi)
#endif
      DO i = 1, ntet_proc
        i_tile = dom_proc(i)
        CALL GET_DEMAG_HELPERS(vertex(:,tet(1,i)), vertex(:,tet(2,i)), vertex(:,tet(3,i)), vertex(:,tet(4,i)), &
                               tet_P(:,:,:,i_tile), tet_D(:,:,i_tile), tet_v(:,:,:,i_tile))
      END DO

#if defined(MPI_OPT)
      CALL MPI_BARRIER(comm_world, ierr_mpi)
      IF (shar_rank.EQ.master) THEN
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tet_P,    3*3*4*ntet, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tet_D,    3*4*ntet,   MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tet_v,  3*3*4*ntet,   MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
      END IF
      CALL MPI_BARRIER(comm_world, ierr_mpi)
#endif
      
      END SUBROUTINE mumaterial_init_demag

!-----------------------------------------------------------------------
! mumaterial_init_neighbors: Finds all neighbors of an element and calculates N_store
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_init_neighbors()

      INTEGER :: i, j, k, i_tile, j_tile
      INTEGER :: nbrs_proc_min, nbrs_proc_max, nc
      DOUBLE PRECISION ::  x ,y, z, dx, dy, dz, cen_i(3)
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: is_local

      ALLOCATE(is_local(ntet))
      is_local = .FALSE.
      DO i = 1, ntet_proc
        is_local(dom_proc(i)) = .TRUE.
      END DO

      !-----------------------------------------
      ! Get largest neighbor count for allocation first
      ! (ignore local elements)
      !-----------------------------------------
      ALLOCATE(nbrs_count(ntet_proc))
      nbrs_count = 0
      nbrs_proc_max = 0
      DO i = 1, ntet_proc
        nc = 0
        x = tet_cen_proc(1,i)
        y = tet_cen_proc(2,i)
        z = tet_cen_proc(3,i)
        DO j = 1, ntet
          IF (is_local(j)) CYCLE
          dx = x - tet_cen(1,j)
          dy = y - tet_cen(2,j)
          dz = z - tet_cen(3,j)
          IF (dx*dx+dy*dy+dz*dz .LE. (pF*tet_rad(j))**2) nc = nc + 1
        END DO
        nbrs_proc_max = MAX(nbrs_proc_max, nc)
        nbrs_count(i) = nc
      END DO
      ALLOCATE(nbrs_proc(nbrs_proc_max,ntet_proc))
      !-----------------------------------------
      ! Actual neighbor construction
      !-----------------------------------------
      DO i = 1, ntet_proc
        x = tet_cen_proc(1,i)
        y = tet_cen_proc(2,i)
        z = tet_cen_proc(3,i)
        j = 0
        DO k = 1, ntet
          IF (is_local(k)) CYCLE
          dx = x - tet_cen(1,k)
          dy = y - tet_cen(2,k)
          dz = z - tet_cen(3,k)
          IF (dx*dx+dy*dy+dz*dz .LE. (pF*tet_rad(k))**2) THEN
            j = j + 1
            nbrs_proc(j,i) = k
          END IF
        END DO
      END DO
      DEALLOCATE(is_local)

      !-----------------------------------------
      ! Calculate demagnetization tensor
      !-----------------------------------------
      ALLOCATE(N_store(3,3,nbrs_proc_max,ntet_proc),N_block(3,3,ntet_proc,ntet_proc))
      N_store= 0.0d0
      N_block = 0.0d0
      DO i = 1, ntet_proc
        i_tile = dom_proc(i)
        cen_i = tet_cen(:,i_tile)
        DO j = 1, nbrs_count(i)
          j_tile = nbrs_proc(j,i)
          N_store(:,:,j,i) = GET_DEMAG(tet_P(:,:,:,j_tile), &
                                       tet_D(:,:,j_tile), &
                                       tet_v(:,:,:,j_tile), &
                                       cen_i) 
        END DO 
        ! Full cluster
        DO j = 1, ntet_proc
          j_tile = dom_proc(j)
          N_block(:,:,j,i) =  GET_DEMAG(tet_P(:,:,:,j_tile), &
                                        tet_D(:,:,j_tile), &
                                        tet_v(:,:,:,j_tile), &
                                        cen_i) 
        END DO 
      END DO
      !-----------------------------------------
      ! Print to screen only
      !-----------------------------------------
      nbrs_proc_min = MINVAL(nbrs_count)
#if defined(MPI_OPT)
      IF (lcomm) THEN
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,nbrs_proc_min,1,MPI_INTEGER,MPI_MIN,comm_world,ierr_mpi)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,nbrs_proc_max,1,MPI_INTEGER,MPI_MAX,comm_world,ierr_mpi)
      END IF
#endif
      IF (lverb) THEN
        WRITE(6,'(3X,A,I0,A,I0,A)') 'Neighb. range: [',nbrs_proc_min,', ',nbrs_proc_max,']'
        FLUSH(6)
      END IF

      END SUBROUTINE mumaterial_init_neighbors

  !-----------------------------------------------------------------------
  ! mumaterial_init_partial: "Partial index" stuff, where the "partial index"
  ! is the index in dom_partial_proc, a local array of all elements which
  ! are necessary for neighbor or dipole calculations. Excludes all elements
  ! which are treated inside a cluster.
  !-----------------------------------------------------------------------
      SUBROUTINE mumaterial_init_partial
      IMPLICIT NONE

      INTEGER :: i, i_tile, j, j_tile, j_partial, wr_dex, idx, n
      LOGICAL, DIMENSION(:), ALLOCATABLE :: mask, isdipole
      INTEGER, DIMENSION(:), ALLOCATABLE :: nonclusters, cluster, partials_1d, recvcounts, displs
      INTEGER, DIMENSION(:,:), ALLOCATABLE ::  nbrs_temp
      INTEGER :: ntet_mid_proc_min, ntet_mid_proc_max, ntet_dip_min, ntet_dip_max, maxsize, nbrs_maxc
      DOUBLE PRECISION :: N_temp(3,3)

      !-----------------------------------------------------------------------
      ! iscluster_proc(i): TRUE  if cluster i should be treated as a cluster
      !                    FALSE if cluster i should be treated as dipoles
      ! 1. Cluster far enough away from our cluster centroid -> TRUE
      ! 2. Cluster contains any neighbors from any of our local elements -> FALSE
      !-----------------------------------------------------------------------
      ALLOCATE(iscluster_proc(world_size))
      DO i = 1, world_size
        iscluster_proc(i) = NORM2(r_cluster(:,i)-r_cluster(:,world_rank+1))>cutoff*d_cluster(i)
      END DO
      
      ALLOCATE(mask(ntet))
      mask = .FALSE.
      DO i = 1, ntet_proc 
        mask(nbrs_proc(1:nbrs_count(i),i)) = .TRUE. ! mask of any neighbors of elements in A
      END DO
      DO i = 1, world_size
        IF (.NOT.iscluster_proc(i)) CYCLE ! Already not a cluster, skip
        IF (ANY(mask(dom_cluster(1:dom_sizes(i),i)))) iscluster_proc(i) = .FALSE. 
      END DO
      DEALLOCATE(mask)
      !---------------------------------------------------------------------
      ! dom_partial_proc: Domain of all elements which are not treated by
      ! clusters. Used for both neighbor and dipole interactions.
      !---------------------------------------------------------------------
      ALLOCATE(nonclusters(COUNT(.NOT.iscluster_proc)))
      nonclusters = PACK([(i, i=1, world_size)], MASK=.NOT.iscluster_proc)
      ntet_partial_proc = SUM(dom_sizes(nonclusters)) ! 
      
      ALLOCATE(dom_partial_proc(ntet_partial_proc))
      idx = 0
      DO i = 1, SIZE(nonclusters)
        j = nonclusters(i)
        cluster = PACK(dom_cluster(:,j),MASK=dom_cluster(:,j)>0)
        n = SIZE(cluster)
        dom_partial_proc(idx+1:idx+n) = cluster
        idx = idx + n
      END DO
      DEALLOCATE(nonclusters)
      !---------------------------------------------------------------------
      ! map_g2p: Maps global index to index in dom_partial_proc
      !---------------------------------------------------------------------
      ALLOCATE(map_g2p(ntet),map_g2all(ntet))
      map_g2p = -1
      DO i = 1, ntet_partial_proc
        map_g2p(dom_partial_proc(i)) = i 
      END DO
      idx = 0
      DO i = 1, world_size
        DO j = 1, dom_sizes(i)
          map_g2all(dom_cluster(j,i)) = idx + j
        END DO
        idx = idx + dom_sizes(i)
      END DO
      !---------------------------------------------------------------------
      ! nbrs_proc should use partial index, not global index
      !---------------------------------------------------------------------
      nbrs_maxc = MAXVAL(nbrs_count)
      ALLOCATE(nbrs_temp(nbrs_maxc,ntet_proc))
      nbrs_temp = 0
      DO i = 1, ntet_proc
        i_tile = dom_proc(i)
        DO j = 1, nbrs_count(i)
          j_tile = nbrs_proc(j,i) ! Compare i_tile and j_tile; global indices
          j_partial = map_g2p(j_tile)
          IF ((i_tile.EQ.j_tile).AND.(j.NE.1)) THEN ! Found self
            nbrs_temp(j,i) = nbrs_temp(1,i)         ! Move 1st entry to j
            nbrs_temp(1,i) = j_partial              ! Put self in 1st
            ! Shuffle N_store
            N_temp           = N_store(:,:,1,i) ! Self
            N_store(:,:,1,i) = N_store(:,:,j,i)
            N_store(:,:,j,i) = N_temp
          ELSE
            nbrs_temp(j,i) = j_partial
          ENDIF
        END DO
      END DO
      nbrs_proc = nbrs_temp
      DEALLOCATE(nbrs_temp)
      !--------------------------------------------------------------------
      ALLOCATE(sync_counts(world_size),sync_displs(world_size))
      DO i = 1, world_size
        sync_counts(i) = 3*dom_sizes(i)
      END DO
      sync_displs = 0
      DO i = 2, world_size
        sync_displs(i) = sync_displs(i-1)+sync_counts(i-1)
      END DO
      !-----------------------------
      ! Fetch into contiguous arrays
      !-----------------------------  
      ALLOCATE(r_partial(3,ntet_partial_proc),V_partial(ntet_partial_proc))
      DO i = 1, ntet_partial_proc
        r_partial(:,i) = tet_cen(:,dom_partial_proc(i))
        V_partial(i) = tet_vol(dom_partial_proc(i))
      END DO
      !---------------------------------------------------------------------
      ! Display:
      ! ntet_dip_max: maximum # of dipoles for any local element
      ! ntet_dip_min: minimum # of dipoles for any local element
      !---------------------------------------------------------------------
      ntet_dip_max = 0
      ntet_dip_min = ntet_partial_proc
      ALLOCATE(isdipole(ntet_partial_proc))
      DO i = 1, ntet_proc
        IF (nbrs_count(i).EQ.0) THEN
          ntet_dip_max = ntet_partial_proc
          EXIT
        ELSE
          isdipole = .TRUE.
          DO j = 1, nbrs_count(i)
            j_partial = nbrs_proc(j,i)
            isdipole(j_partial) = .FALSE.
          END DO
          n = COUNT(isdipole)
          ntet_dip_max = MAX(ntet_dip_max, n)
          ntet_dip_min = MIN(ntet_dip_min, n)
        END IF
      END DO
      DEALLOCATE(isdipole)
#if defined(MPI_OPT)
        IF (lcomm) THEN
          CALL MPI_ALLREDUCE(ntet_dip_min,ntet_mid_proc_min,1,MPI_INTEGER,MPI_MIN,comm_world,ierr_mpi)        
          CALL MPI_ALLREDUCE(ntet_dip_max,ntet_mid_proc_max,1,MPI_INTEGER,MPI_MAX,comm_world,ierr_mpi) 
        END IF
#endif
        IF (lverb) THEN
          WRITE(6,'(3X,A,I0,A,I0,A)') 'Dipole range : [',ntet_mid_proc_min,', ',ntet_mid_proc_max,']'
          FLUSH(6)
        END IF

      END SUBROUTINE mumaterial_init_partial

  !-----------------------------------------------------------------------
  ! mumaterial_split_world: Splits a domain into disjoint domains dom_shar,
  ! one for each MPI node, and all spatially local and similar sized.
  ! See mumaterial_split_sub.
  !-----------------------------------------------------------------------
      SUBROUTINE mumaterial_split_world

      IMPLICIT NONE

      INTEGER :: color, i, n_targ, size2, r, reci
      LOGICAL :: lactive
      INTEGER, DIMENSION(:), ALLOCATABLE :: dom_in, dom2
      DOUBLE PRECISION, PARAMETER :: targ = 0.5
      INTEGER :: ntet_shar_min, ntet_shar_max
      INTEGER, PARAMETER :: TAG_SIZE = 1234, TAG_DOM = 1235, TAG_SPLT = 1236

      IF (lverb)  WRITE(6,*)  ' ------- Domain Division ------'
      ! Only masters get a color
      color = 0
#if defined(MPI_OPT)
      IF (shar_rank.NE.0) THEN
        color = 1
      ELSE
        color = (world_rank*master_size)/world_size
      END IF
#endif
      !-----------------------------------------
      ! Boot up world_rank = 0 to start; non-MPI
      !-----------------------------------------
      lactive = (color.EQ.0)
      IF (lactive) THEN ! 
        ALLOCATE(dom_shar(ntet))
        DO i = 1, ntet
          dom_shar(i) = i
        END DO
        ntet_shar = SIZE(dom_shar)
        ntet_shar_min = ntet_shar
        ntet_shar_max = ntet_shar
      END IF
      !-----------------------------------------
      ! Recursively boot up other masters
      !-----------------------------------------
#if defined(MPI_OPT)   
      IF ((master_size.GT.1) .AND. (shar_rank.EQ.master)) THEN         
        r = NINT(LOG(DBLE(master_size))/LOG(2.0)) ! log_2(X) = ln(X)/log(2)
        DO WHILE (r.GT.0)
          IF (lactive) THEN 
            ! Free up dom_shar
            ALLOCATE(dom_in(ntet_shar))
            dom_in = dom_shar
            DEALLOCATE(dom_shar)
            ! Split
            n_targ = 2**(r)
            CALL mumaterial_split_sub(dom_in,targ,n_targ,dom_shar,dom2) ! Split box
            DEALLOCATE(dom_in)
            ntet_shar  = SIZE(dom_shar)
            size2 = SIZE(dom2)            
            r = r-1 

            ! now mail one of new boxes to the appropriate recipient
            reci = color + 2**r 
            CALL MPI_SEND(size2,     1, MPI_INTEGER, reci, TAG_SIZE, comm_master, ierr_mpi) 
            CALL MPI_SEND(dom2,  size2, MPI_INTEGER, reci, TAG_DOM,  comm_master, ierr_mpi)
            DEALLOCATE(dom2) 
            CALL MPI_SEND(r,         1, MPI_INTEGER, reci, TAG_SPLT, comm_master, ierr_mpi)
          ELSE
            ! Waiting part
            CALL MPI_RECV(ntet_shar,        1, MPI_INTEGER, MPI_ANY_SOURCE, TAG_SIZE, comm_master, mstat, ierr_mpi)
            ALLOCATE(dom_shar(ntet_shar))
            CALL MPI_RECV(dom_shar, ntet_shar, MPI_INTEGER, MPI_ANY_SOURCE, TAG_DOM,  comm_master, mstat, ierr_mpi)
            CALL MPI_RECV(r,                1, MPI_INTEGER, MPI_ANY_SOURCE, TAG_SPLT, comm_master, mstat, ierr_mpi)
            lactive = .TRUE. ! Activate node
          END IF
        END DO
      END IF
      !-----------------------------------------
      ! Share information with ranks within node
      !-----------------------------------------
      CALL MPI_Bcast(ntet_shar,         1, MPI_INTEGER, 0, comm_shar, ierr_mpi)
      IF (shar_rank.NE.master) ALLOCATE(dom_shar(ntet_shar))      
      CALL MPI_Bcast(dom_shar, ntet_shar,  MPI_INTEGER, 0, comm_shar, ierr_mpi)
      CALL MPI_BARRIER(comm_world, ierr_mpi)
      CALL MPI_ALLREDUCE(ntet_shar,ntet_shar_min,1,MPI_INTEGER,MPI_MIN,comm_world,ierr_mpi)
      CALL MPI_ALLREDUCE(ntet_shar,ntet_shar_max,1,MPI_INTEGER,MPI_MAX,comm_world,ierr_mpi)  
#endif
      !------------------------------------------------------------------------
      ! Done; print to screen
      !------------------------------------------------------------------------
      IF (lverb) THEN
        WRITE(6,'(3X,A,I7)')        'MPI Nodes    : ',master_size
        WRITE(6,'(3X,A,I0,A,I0,A)') 'Node range   : [',ntet_shar_min,', ',ntet_shar_max,']'
        FLUSH(6)
      END IF

      RETURN
      END SUBROUTINE mumaterial_split_world

  !-----------------------------------------------------------------------
  ! mumaterial_split_share: Splits a node's dom_shar into several disjoint
  ! dom_proc's, one for each rank, each spatially localized and similarly sized.
  ! See mumaterial_split_sub.
  !-----------------------------------------------------------------------
      SUBROUTINE mumaterial_split_share

      IMPLICIT NONE

      INTEGER :: n_targ, n_spl_1, n_spl_2, reci, src, i
      DOUBLE PRECISION :: targ
      LOGICAL :: lactive
      INTEGER, ALLOCATABLE :: dom_out_2(:), temp_dom(:)
      INTEGER :: ntet_proc_min
      DOUBLE PRECISION :: min_tet_rad, min_tet_rad_all, max_tet_rad_all
      INTEGER, PARAMETER :: TAG_SIZE = 101, TAG_DOM = 102, TAG_SPLT = 103
      !-----------------------------------------
      ! If only one node (or not MPI)
      !-----------------------------------------
      IF (shar_size.EQ.1) THEN
        ntet_proc = SIZE(dom_shar)
        ALLOCATE(dom_proc(ntet_proc))
        dom_proc = dom_shar
      END IF

      IF (lcomm .AND. (shar_size.GT.1)) THEN
#if defined(MPI_OPT)   
        IF (shar_rank.EQ.master) THEN 
          !-----------------------------------------
          ! Master logic A: Build first 2 domains 
          !-----------------------------------------
          n_targ = shar_size         ! Example: 13
          n_spl_1 = n_targ/2         ! Example: 7
          n_spl_2 = n_targ - n_spl_1 ! Example: 6
          targ = DBLE(n_spl_1)/DBLE(n_targ) ! 7/13 ~ 0.54
          reci = master + 1
          ! Divide box initially and send to two ranks
          CALL mumaterial_split_sub(dom_shar,targ,n_targ,dom_proc,dom_out_2)
          ntet_proc = SIZE(dom_proc)
          CALL MPI_SEND(SIZE(dom_out_2),          1, MPI_INTEGER, reci, TAG_SIZE, comm_shar, ierr_mpi) 
          CALL MPI_SEND(dom_out_2,  SIZE(dom_out_2), MPI_INTEGER, reci, TAG_DOM,  comm_shar, ierr_mpi)
          CALL MPI_SEND(n_spl_2,                 1, MPI_INTEGER,  reci, TAG_SPLT, comm_shar, ierr_mpi)
          IF (shar_size.GT.2) THEN ! Careful with 2 threads; then need to keep one.
            reci = reci + 1
            CALL MPI_SEND(ntet_proc,       1, MPI_INTEGER, reci, TAG_SIZE, comm_shar, ierr_mpi) 
            CALL MPI_SEND(dom_proc,ntet_proc, MPI_INTEGER, reci, TAG_DOM,  comm_shar, ierr_mpi)
            CALL MPI_SEND(n_spl_1,         1, MPI_INTEGER, reci, TAG_SPLT, comm_shar, ierr_mpi)
            DEALLOCATE(dom_proc)
          END IF
          DEALLOCATE(dom_out_2)
          !-----------------------------------------
          ! Master logic B: Receive and send domains
          !-----------------------------------------
          IF (shar_size.GT.2) THEN
            DO
              reci = reci + 1
              ! Receive from any source
              CALL MPI_RECV(ntet_proc,        1, MPI_INTEGER, MPI_ANY_SOURCE, TAG_SIZE, comm_shar, mstat, ierr_mpi)
              src = mstat(MPI_SOURCE)
              ALLOCATE(dom_proc(ntet_proc))
              CALL MPI_RECV(dom_proc, ntet_proc, MPI_INTEGER, src, TAG_DOM,  comm_shar, mstat, ierr_mpi)
              CALL MPI_RECV(n_targ,           1, MPI_INTEGER, src, TAG_SPLT, comm_shar, mstat, ierr_mpi)
              IF (reci.EQ.shar_size) THEN
                EXIT ! That's us!
              ELSE
                CALL MPI_SEND(ntet_proc,        1, MPI_INTEGER, reci, TAG_SIZE, comm_shar, ierr_mpi) 
                CALL MPI_SEND(dom_proc, ntet_proc, MPI_INTEGER, reci, TAG_DOM,  comm_shar, ierr_mpi)
                CALL MPI_SEND(n_targ,           1, MPI_INTEGER, reci, TAG_SPLT, comm_shar, ierr_mpi)
                DEALLOCATE(dom_proc)
              END IF

            END DO
          END IF
          !-----------------------------------------
        ELSE ! Non-masters
          lactive = .FALSE.
          DO
          !-----------------------------------------
          ! Rank logic A: Receive domain from master
          !-----------------------------------------
            IF (.NOT.lactive) THEN
              CALL MPI_RECV(ntet_proc,        1, MPI_INTEGER, master, TAG_SIZE, comm_shar, mstat, ierr_mpi) 
              ALLOCATE(dom_proc(ntet_proc))
              CALL MPI_RECV(dom_proc, ntet_proc, MPI_INTEGER, master, TAG_DOM,  comm_shar, mstat, ierr_mpi)       
              CALL MPI_RECV(n_targ,           1, MPI_INTEGER, master, TAG_SPLT, comm_shar, mstat, ierr_mpi)
            END IF
            ! Check if should divide
            IF (n_targ.EQ.1) EXIT ! That's us!
          !-----------------------------------------
          ! Rank logic B: Split domain then send to master
          !-----------------------------------------
            n_spl_1 = n_targ/2         
            n_spl_2 = n_targ - n_spl_1 
            targ = DBLE(n_spl_1)/DBLE(n_targ) 
            ALLOCATE(temp_dom(ntet_proc))
            temp_dom = dom_proc
            DEALLOCATE(dom_proc)
            CALL mumaterial_split_sub(temp_dom,targ,n_targ,dom_proc,dom_out_2)
            ntet_proc = SIZE(dom_proc)
            CALL MPI_SEND(SIZE(dom_out_2),          1, MPI_INTEGER, master, TAG_SIZE, comm_shar, ierr_mpi) 
            CALL MPI_SEND(dom_out_2,  SIZE(dom_out_2), MPI_INTEGER, master, TAG_DOM,  comm_shar, ierr_mpi)
            CALL MPI_SEND(n_spl_2,                  1, MPI_INTEGER, master, TAG_SPLT, comm_shar, ierr_mpi)
            DEALLOCATE(dom_out_2,temp_dom)

            n_targ = n_spl_1
            lactive = .TRUE.
          END DO
        END IF
      END IF
      CALL MPI_BARRIER(comm_shar, ierr_mpi)
      ALLOCATE(tet_cen_proc(3,ntet_proc),max_tet_rad(world_size),dom_sizes(world_size))
      tet_cen_proc = tet_cen(:, dom_proc(1:ntet_proc))
      dom_sizes = 0; dom_sizes(world_rank+1) = ntet_proc
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, dom_sizes, world_size, MPI_INTEGER, MPI_SUM, comm_world, ierr_mpi )
      ntet_proc_max = MAXVAL(dom_sizes); ntet_proc_min = MINVAL(dom_sizes)
      max_tet_rad = 0.0d0; max_tet_rad(world_rank+1) = MAXVAL(tet_rad(dom_proc)); min_tet_rad = MINVAL(tet_rad(dom_proc))
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,max_tet_rad,world_size,MPI_DOUBLE_PRECISION,MPI_SUM,comm_world,ierr_mpi)
      CALL MPI_ALLREDUCE(min_tet_rad,min_tet_rad_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm_world,ierr_mpi)
      max_tet_rad_all = MAXVAL(max_tet_rad)
      IF (lverb) THEN
        WRITE(6,'(3X,A,I7)')        'MPI Threads  : ',world_size
        WRITE(6,'(3X,A,I0,A,I0,A)') 'Thread range : [',ntet_proc_min,', ',ntet_proc_max,']'
        WRITE(6,'(3X,A,ES9.2,A,/,19X,ES9.2,A)') 'Element size :  ',min_tet_rad_all,' m',max_tet_rad_all,' m'
        FLUSH(6)
      END IF
#endif
      ALLOCATE(vol_proc(ntet_proc))
      DO i = 1, ntet_proc
        vol_proc(i) = tet_vol(dom_proc(i))
      END DO
      RETURN
      END SUBROUTINE mumaterial_split_share

  !-----------------------------------------------------------------------
  ! mumaterial_split_sub: Binary splits a collection of elements (boxin) into
  ! two spatially localized subdomains (box1, box2) based on the position
  ! of the elements and the target split (targ). Box sizes can differ by
  ! a fraction from target to try and minimize spatial extent.
  !-----------------------------------------------------------------------
  ! param[in]: boxin. indices of elements
  ! param[in]: targ. relative size of box1 out compared to boxin.
  ! param[in]: ntarg. number of boxes remaining.
  ! param[out]: box1. first output subdomain
  ! param[out]: box2. second output subdomain.
  !-----------------------------------------------------------------------
    SUBROUTINE mumaterial_split_sub(boxin,targ,ntarg,box1,box2)        
      USE qsort ! quicksort

      IMPLICIT NONE

      INTEGER, INTENT(in)               :: boxin(:)
      DOUBLE PRECISION, INTENT(in)      :: targ
      INTEGER, INTENT(in)               :: ntarg
      INTEGER, ALLOCATABLE, INTENT(out) :: box1(:), box2(:) 

      INTEGER :: boxsize, size1, size2
      INTEGER :: i, j, k
      ! LAPACK
      DOUBLE PRECISION :: cov(3,3), evals(3), evec(3), WORK(9)
      INTEGER          :: INFO, LWORK
      ! Sorting
      DOUBLE PRECISION :: r_com(3)
      INTEGER,          DIMENSION(:), ALLOCATABLE   :: idx
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE   :: proj
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: P
      ! Window slide
      DOUBLE PRECISION :: r, eps
      DOUBLE PRECISION :: pmin_fix(3,2), pmax_fix(3,2), Lmin(3), Lmax(3), Rmin(3), Rmax(3), d2, d2_best
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: pref_min, pref_max, suff_min, suff_max
      INTEGER :: iL_end, iM_beg, iM_end, iR_beg, W, hW,  k_best
      DOUBLE PRECISION, PARAMETER  :: eps_max = 0.1d0 ! 

      ! Range of allowed sizes
      boxsize = SIZE(boxin)
      IF (boxsize.EQ.1) THEN
        CALL mumaterial_abort()
      ELSEIF (boxsize.EQ.2) THEN
        ALLOCATE(box1(1),box2(1))
        box1(1) = boxin(1)
        box2(1) = boxin(2)
        RETURN
      ENDIF
      !---------------------------------------
      ! Covariant matrix for finding vector to split along
      r_com = SUM(tet_cen(:,boxin),DIM=2)/boxsize
      cov = 0.0d0
      DO i = 1, boxsize
          DO j = 1, 3
              DO k = 1, 3
                  cov(j,k) = cov(j,k) + (tet_cen(j,boxin(i)) - r_com(j)) * &
                                        (tet_cen(k,boxin(i)) - r_com(k))
              END DO
          END DO
      END DO
      cov = cov / DBLE(boxsize)
      ! Eigenvectors cov and eigenvalues evals
      LWORK = 9
      CALL DSYEV('V','U',3,cov,3,evals,WORK,LWORK,INFO)
      IF (INFO.NE.0) THEN
        WRITE(6,"(2X,A,I0,A,I0)") "ERROR: Rank ", world_rank, " / DSYEV failed, INFO = ", INFO 
        CALL mumaterial_abort()
      END IF
      ! Project coordinates along "longest" eigenvector, then sort
      ALLOCATE(proj(boxsize),idx(boxsize))
      evec = cov(:,MAXLOC(evals,DIM=1))
      proj = [(DOT_PRODUCT(evec, (tet_cen(:,boxin(i)) - r_com)),i=1,boxsize)]
      idx = [(i, i=1, boxsize)]
      CALL quicksort(proj,idx,1,boxsize)
      DEALLOCATE(proj)
      !---------------------------------------
      ! Try varying the number of elements; 
      ! Determine oriented bounding box for each domain size
      ! and the maximum diameter of each bounding box
      ALLOCATE(P(3,boxsize)) ! Avoid MATMUL
      DO i = 1, boxsize
        P(:,i) = tet_cen(:,boxin(idx(i)))
      END DO
      !--------
      ! Target size
      size1 = NINT(targ*boxsize)
      r = LOG(DBLE(ntarg))/LOG(2.0D0)
      eps = (1+eps_max)**(1/r)-1.0d0
      hW = INT(eps/2.0d0*DBLE(boxsize))
      ! Window slide
      W = 2*hW
      iL_end = size1 - hW
      iM_beg = iL_end + 1
      iM_end = iL_end + W
      iR_beg = iM_end + 1

      ! Fixed windows
      pmin_fix(:, 1) = MINVAL(P(:, 1:iL_end), DIM=2)
      pmax_fix(:, 1) = MAXVAL(P(:, 1:iL_end), DIM=2)
      pmin_fix(:, 2) = MINVAL(P(:, iR_beg:boxsize), DIM=2)
      pmax_fix(:, 2) = MAXVAL(P(:, iR_beg:boxsize), DIM=2)
      
      ! Prefix window
      ALLOCATE(pref_min(3,W), pref_max(3,W))
      pref_min(:,1) = P(:, iM_beg)
      pref_max(:,1) = P(:, iM_beg)
      DO k = 2, W
        pref_min(:,k) = MIN(pref_min(:,k-1), P(:,iM_beg-1+k))
        pref_max(:,k) = MAX(pref_max(:,k-1), P(:,iM_beg-1+k))
      END DO

      ! Suffix window
      ALLOCATE(suff_min(3,W), suff_max(3,W))
      suff_min(:,W) = P(:, iM_end)
      suff_max(:,W) = P(:, iM_end)
      DO k = W-1, 1, -1
        suff_min(:,k) = MIN(suff_min(:,k+1), P(:,iM_beg-1+k))
        suff_max(:,k) = MAX(suff_max(:,k+1), P(:,iM_beg-1+k))
      END DO      

      ! Find best case
      d2_best = HUGE(0.d0)
      k_best = 0
      DO k = 0, W
        ! Left = fixed left ∪ prefix of middle window
        IF (k .GT. 0) THEN
          Lmin = MIN(pmin_fix(:,1), pref_min(:,k))
          Lmax = MAX(pmax_fix(:,1), pref_max(:,k))
        ELSE
          Lmin = pmin_fix(:,1);  Lmax = pmax_fix(:,1)
        END IF
        ! Right = suffix of middle window ∪ fixed right
        IF (k .LT. W) THEN
          Rmin = MIN(suff_min(:,k+1), pmin_fix(:,2))
          Rmax = MAX(suff_max(:,k+1), pmax_fix(:,2))
        ELSE
          Rmin = pmin_fix(:,2);  Rmax = pmax_fix(:,2)
        END IF
        d2 = MAX(SUM((Lmax-Lmin)**2), SUM((Rmax-Rmin)**2))
        IF (d2 .LT. d2_best) THEN
          d2_best = d2; k_best = k
        END IF
      END DO
      DEALLOCATE(P, suff_min, suff_max, pref_min, pref_max)
      !---------------------------------------
      size1 = size1-hW+k_best
      size2 = boxsize - size1
      ! Now allocate the boxes for sending
      ALLOCATE(box1(size1), box2(size2))
      box1 = boxin(idx(1:size1))
      box2 = boxin(idx(size1+1:boxsize))
      DEALLOCATE(idx)

      END SUBROUTINE mumaterial_split_sub

        
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
        WRITE(iunit,'(3X,A,F9.3)')  'Pad factor   : ',pF
        WRITE(iunit,'(3X,A,I9)')    'Max Iter.    : ',maxIter
        WRITE(iunit,'(3X,A,ES9.2)') 'Max Error    : ',threshold
        WRITE(iunit,'(3X,A,F9.3)')  'Lambda start : ',lambdaStart
        WRITE(iunit,'(3X,A,F9.3)')  'Lambda fact. : ',lambdaFactor
        WRITE(iunit,'(3X,A,I9)')    'Lambda thrsh.: ',lambdaThresh
        WRITE(iunit,'(3X,A,F7.2,A)')'Converged at : ',convCheck,' %'
      END IF
      WRITE(iunit,'(A)')           ' -----  Magnetic structure  ----'
      WRITE(iunit,'(3X,A,A)')      'File: ',TRIM(file_string)
      WRITE(iunit,'(3X,A,A)')      'Model Name   : ',TRIM(machine_string)
      WRITE(iunit,'(3X,A,A)')      'Date         : ',TRIM(date_string)
      WRITE(iunit,'(3X,A,I9)')     'Vertices     : ',nvertex
      WRITE(iunit,'(3X,A,I9)')     'Tetrahedrons : ',ntet
      WRITE(iunit,'(3X,A,I9)')     'State Funcs. : ',nstate
      DO i = 1, nstate
        WRITE(iunit,'(5X,A,I0)') 'State Function ',i
        IF (state_type(i)==1) THEN
          WRITE(iunit,'(7X,A)') 'Type: Hard Magnet'
          WRITE(iunit,'(7X,A,ES12.3)')    '  Mu   :',constant_mu(i)
          WRITE(iunit,'(7X,A,ES12.3)')    '  Mu_o :',constant_mu_o(i)
          WRITE(iunit,'(7X,A,3(ES12.3))') '  M_rem :',M_rem(:,i)
        ELSEIF (state_type(i)==2) THEN
          k = SIZE(statefunc(i)%H)
          WRITE(iunit,'(7X,A)')           '  Type : Soft Magnet (H-M)'
          WRITE(iunit,'(7X,A,I3)')        'NKnots :',k
          WRITE(iunit,'(7X,A,2(ES12.3))') '     H :',statefunc(i)%H(1),statefunc(i)%H(k)
          WRITE(iunit,'(7X,A,2(ES12.3))') '     M :',statefunc(i)%M(1),statefunc(i)%M(k)
        ELSEIF (state_type(i)==3) THEN
          WRITE(iunit,'(7X,A)') 'Type: Soft Magnet (mu constant)'
          WRITE(iunit,'(7X,A,F12.3)')    '    mu :',constant_mu(i)
        ELSE
          WRITE(iunit,'(7X,A,I3)') 'Type: UNKNOWN (ERROR) state_type=',state_type(i)
        END IF
      END DO
      FLUSH(iunit)

      END SUBROUTINE mumaterial_info
  
!------------------------------------------------------------------------------
! mumaterial_run: Initial calculations, does MPI, and calls iterations
!------------------------------------------------------------------------------
! param[in]: offset. Offset of all tiles from the origin
! param[in], optional: filename. Magnetization file.
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_run(vert_offset, filename)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in), OPTIONAL :: vert_offset(3)
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: filename
      
      ! Calculate tetrahedron quantities
      CALL mumaterial_init_mesh(vert_offset) 
      ! Split domain across MPI nodes
      CALL mumaterial_split_world()   
      ! Split MPI subdomains across MPI ranks
      CALL mumaterial_split_share() 
      ! Sort elements by material type
      CALL mumaterial_init_states()
      ! Calculate background field for rank-local elements
      CALL mumaterial_init_happ()
      ! Setup for demagnetization tensors
      CALL mumaterial_init_demag()
      ! Get nearest neighbors
      CALL mumaterial_init_neighbors()
      ! Build clusters
      CALL mumaterial_cluster_setup()
      ! Build dom_partial_proc
      CALL mumaterial_init_partial()   
      ! Allocate M_local, M_partial, etc.
      CALL mumaterial_alloc_init() 
      ! Read magnetization file and sync
      IF (PRESENT(filename)) CALL mumaterial_magfile_read(filename) 

#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_BARRIER(comm_world, ierr_mpi)
#endif
      IF (lverb) WRITE (6,*) ' ------- MUMAT init done -------'

      CALL mumaterial_iterate()
      ! Deallocate Helpers
      CALL mumaterial_dealloc_init()
      ! Synchronize M_global for output
      CALL mumaterial_syncmag(SYNC_DONE)

      RETURN
      END SUBROUTINE mumaterial_run



!-----------------------------------------------------------------------
! mumaterial_iterate: Iteration loop
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_iterate()
      IMPLICIT NONE

      ! Picard loop
      INTEGER :: iter, i, j, i_tile
      INTEGER :: stype, sdex
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: res_M, res_M_prev
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lambda_n
      DOUBLE PRECISION :: M_targ(3), chi, M_targ_norm
      LOGICAL :: lfulldipole
      DOUBLE PRECISION, PARAMETER :: lambda_min = 0.1d0
      INTEGER, PARAMETER :: iter_recalc = 1
      ! Convergence, residuals
      DOUBLE PRECISION :: conv_loc, conv_glob
      DOUBLE PRECISION, ALLOCATABLE :: dW(:), W(:), f(:)
      LOGICAL, ALLOCATABLE :: is_conv(:)
      DOUBLE PRECISION :: r_M_loc, r_M_max
      DOUBLE PRECISION :: dW_cl,  W_cl,  r_W_cl 
      DOUBLE PRECISION :: dW_all, W_all, r_W_all
      DOUBLE PRECISION :: M2, c1, c2
      ! Verbose
      INTEGER ::          i_bad
      DOUBLE PRECISION :: pair_in(2), pair_out(2)
      DOUBLE PRECISION :: M_targ_bad, H_norm_bad, M_norm_bad, lambda_bad
      DOUBLE PRECISION :: info_max(6)
      INTEGER, PARAMETER :: TAG_WORST = 4547
      !-----------------------------------------------------------------------!
      ALLOCATE(res_M(3,ntet_proc),res_M_prev(3,ntet_proc))
      ALLOCATE(lambda_n(ntet_proc))
      ALLOCATE(is_conv(ntet_proc))
      ALLOCATE(f(ntet_proc),W(ntet_proc),dW(ntet_proc))
      H_prev = 0.0d0
      lambda_n = lambdaStart
      res_M = 0.0d0
      is_conv = .FALSE.
      f = 0.5d0*MU0*vol_proc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------ PRIMARY LOOP ---------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO iter = 1, maxiter
        !-----------------------------
        ! Update background field
        !-----------------------------   
        lfulldipole = (MOD(iter,iter_recalc).EQ.0)
        CALL mumaterial_applied_field() ! Sets H_ext = H_app
        CALL mumaterial_cluster_field() ! Sets H_ext = H_ext + H_cluster
        CALL mumaterial_dipole_field(lfulldipole) ! Calculates H_mid
        H_ext = H_ext + H_mid
        
        !-----------------------------------
        ! Block solve for local cluster
        !-----------------------------------  
        CALL mumaterial_block_solve() 
        ! Reset residues        
        dW_cl = 0.0d0; W_cl = 0.0d0;
        res_M_prev = res_M
        r_M_max = -1.0d0
        DO i = 1, ntet_proc
          res_M(:,i)   = M_spline(:,i) - M_local(:,i) ! target - old
          M_local(:,i) = M_local(:,i) + lambda_n(i)*res_M(:,i)
        END DO

        ! Magnetization residual
        r_M_loc = 0.0d0
        DO i = 1, ntet_proc
          M2 = DOT_PRODUCT(M_spline(:,i), M_spline(:,i))+small
          r_M_loc = SQRT(DOT_PRODUCT(res_M(:,i),res_M(:,i))/M2)
          is_conv(i) = (r_M_loc.LT.threshold)  
          IF (r_M_loc.GT.r_M_max) THEN
            r_M_max = r_M_loc
            i_bad = i
          END IF
        END DO
        info_max = (/r_M_max, DBLE(dom_proc(i_bad)), NORM2(H_prev(:,i_bad)),  NORM2(M_local(:,i_bad)), NORM2(res_M(:,i_bad)), lambda_n(i_bad)/)
        conv_loc = 0.0
        DO i = 1, ntet_proc
          IF (is_conv(i)) conv_loc = conv_loc + vol_proc(i) ! Track converged volume
        END DO

        ! Energy residual
        dW = f * (H_prev(:,1)*res_M(:,1) + H_prev(:,2)*res_M(:,2) + H_prev(:,3)*res_M(:,3))
        W  = f * (H_prev(:,1)*(H_prev(:,1)+M_local(:,1)) + &
                  H_prev(:,2)*(H_prev(:,2)+M_local(:,2)) + &
                  H_prev(:,3)*(H_prev(:,3)+M_local(:,3)))
        W_cl = SUM(ABS(W))
        dW_cl = SUM(ABS(dW))

        !-----------------------------------
        ! Handle lambda
        !-----------------------------------  
        DO i = 1, ntet_proc
          IF (is_conv(i)) THEN
            ! Change nothing
          ELSE IF ((DOT_PRODUCT(res_M(:,i),res_M_prev(:,i))<0.0) .AND. (NORM2(res_M(:,i))>NORM2(res_M_prev(:,i)))) THEN 
            lambda_n(i) = lambda_n(i)*lambdaFactor
          ELSE IF (NORM2(res_M(:,i)).LT.NORM2(res_M_prev(:,i))) THEN
            lambda_n(i) = lambda_n(i) * 1.05
          END IF
          lambda_n(i) = MAX(MIN(lambda_n(i),1.0d0), lambda_min) ! Clamp
        END DO

        !-----------------------------------
        ! Master prints to screen each iteration
        !----------------------------------- 
        r_W_cl = 0.0d0; r_W_all = 0.0d0
        IF (W_cl.GT.small) r_W_cl = dW_cl/W_cl
        IF (lcomm) THEN
#if defined(MPI_OPT)
          ! Worst element
          pair_in(1) = r_M_max; pair_in(2) = DBLE(world_rank)
          CALL MPI_ALLREDUCE(pair_in,pair_out,1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, comm_world, ierr_mpi)
          CALL MPI_BCAST(info_max, 6, MPI_DOUBLE_PRECISION, INT(pair_out(2)), comm_world, ierr_mpi)
          ! Cluster
          CALL MPI_ALLREDUCE(MPI_IN_PLACE, r_W_cl, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_world, ierr_mpi)
          ! World
          CALL MPI_ALLREDUCE(dW_cl, dW_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_world, ierr_mpi) 
          CALL MPI_ALLREDUCE( W_cl,  W_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_world, ierr_mpi) 
          IF (W_all.GT.small) r_W_all = dW_all/W_all
          ! Convergence
          CALL MPI_ALLREDUCE(conv_loc, conv_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_world, ierr_mpi) 
#endif
        ELSE
          r_W_all = r_W_cl
          conv_glob = conv_loc
        END IF
        IF (lverb) THEN 
          IF (iter.EQ.1) THEN
            WRITE(6,'(/,A)') '  iter %done dW/W_all   dW/W_cl  dM/M_max |    tile     Hnorm     Mnorm    dMnorm'
            WRITE(6,*)       '==============================================================================='
          END IF
          WRITE(6,'(1X,I5,F5.1, 3ES10.3,A,I7,3ES10.3,F7.4)') iter, conv_glob/tet_vol_tot*100.0d0, r_W_all, r_W_cl, info_max(1), ' | ', & 
                      INT(info_max(2)),info_max(3),info_max(4),info_max(5), info_max(6)
          CALL FLUSH(6)
        END IF

        !-----------------------------
        ! Synchronize
        !-----------------------------    
        CALL mumaterial_syncmag(SYNC_ITER)
        !-----------------------------
        ! Stop if converged
        !-----------------------------            
        IF (r_W_cl.LE.threshold) EXIT


        CALL MPI_BARRIER(comm_world, ierr_mpi)
        !---------------------------------------------------------------------!
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------- END PRIMARY LOOP -------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (lverb) THEN
        WRITE(6,*)               ' ------ Iterations done ------'
        IF (iter.GE.maxiter) WRITE(6,"(2X,A,I0,A)") "WARNING: exceeded maximum iterations (", maxiter,")"
        FLUSH(6)
      END IF
      CALL mumaterial_cluster_update(R_quad,mom_quad) ! One last update

      DEALLOCATE(res_M,res_M_prev,lambda_n)
      DEALLOCATE(f,dW,W)

      RETURN
      END SUBROUTINE mumaterial_iterate
!-----------------------------------------------------------------------
! mumaterial_syncmag: Synchronizes M across all MPI nodes
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_syncmag(ctype)

      IMPLICIT NONE

      INTEGER, INTENT(in) :: ctype
      INTEGER :: i, j, n_i, n, offset
      INTEGER :: counts(world_size), displs(world_size)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: M_send

      IF (.NOT.lcomm) RETURN ! Nothing to synchronize
      counts = 0
      displs = 0
#if defined(MPI_OPT)
      SELECT CASE(ctype)
      !-----------------------------------------------------------------------
      ! After reading in, world_master needs to send to other shar_masters 
      ! then send M_local to ranks in comm_shar
      !-----------------------------------------------------------------------
        CASE(SYNC_READ)
          IF (shar_rank.EQ.master) THEN
            ! First send to other masters
            IF (.NOT.lismaster) M_global = 0.0d0 
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, M_global, 3*ntet, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )

            ! Get number of elements to send to comm_shar ranks
            n = 0
            DO j = 1, shar_size
              i = shar_ranks_arr(j)+1
              n = n + dom_sizes(i)
            END DO
            ! Build M_send array
            ALLOCATE(M_send(3,n))
            offset = 0
            DO j = 1, shar_size
              i = shar_ranks_arr(j)+1
              n_i = dom_sizes(i)
              M_send(:,offset+1:offset+n_i) = M_global(:,dom_cluster(1:n_i,i))
              counts(j) = 3*n_i
              displs(j) = 3*offset
              offset = offset+n_i
            END DO
          END IF
          n_i = dom_sizes(world_rank+1)
          ! Now send
          IF (shar_rank.NE.master) ALLOCATE(M_send(1,1))
          CALL MPI_Scatterv(M_send,counts,displs,MPI_DOUBLE_PRECISION,M_local,3*n_i,MPI_DOUBLE_PRECISION,0,comm_shar,ierr_mpi)
          DEALLOCATE(M_send)
      !-----------------------------------------------------------------------
      ! During iterations, do allgetherv on all ranks
      !-----------------------------------------------------------------------
        CASE(SYNC_ITER)
          CALL MPI_ALLGATHERV(M_local, 3*ntet_proc, MPI_DOUBLE_PRECISION, &
                              M_all, sync_counts, sync_displs, MPI_DOUBLE_PRECISION, & 
                              comm_world, ierr_mpi)
          DO i = 1, ntet_partial_proc
            M_partial(:,i) = M_all(:, map_g2all(dom_partial_proc(i)))
          END DO
      !-----------------------------------------------------------------------
      ! After iterations, just restructure.
      !-----------------------------------------------------------------------
        CASE(SYNC_DONE)
          IF (.NOT.ALLOCATED(M_global)) ALLOCATE(M_global(3,ntet))
          DO i = 1, ntet
            M_global(:,i) = M_all(:, map_g2all(i))
          END DO
        END SELECT 
#endif 

      END SUBROUTINE mumaterial_syncmag

! !-----------------------------------------------------------------------
! ! mumaterial_picard_rem: Determines M_targ for given H_back. Limits relative
! ! change in H.
! !-----------------------------------------------------------------------
! ! param[in]: H_back. Input magnetic field
! ! param[inout]: H_prev. Last iteration's magnetic field. overwritten.
! ! param[in]: i. local index
! ! param[out]: M_targ. Target magnetization
! ! param[out]: lgood. True if change in H was not limited.
! !-----------------------------------------------------------------------
!       SUBROUTINE mumaterial_picard_rem(H_back, i, M_targ,  lgood)

!       IMPLICIT NONE

!       DOUBLE PRECISION, INTENT(in)  :: H_back(3)
!       INTEGER, INTENT(in)           :: i ! LOCAL index!
!       DOUBLE PRECISION, INTENT(out) :: M_targ(3)
!       LOGICAL, INTENT(out)          :: lgood

!       INTEGER :: iter, i_tile
!       DOUBLE PRECISION :: lambda_k

!       DOUBLE PRECISION :: M_rem_norm, u_ea(3), u_oa_1(3), u_oa_2(3), mu_ea, mu_oa 
!       DOUBLE PRECISION :: H_next(3), H_new(3), H_targ(3), M_targ_old(3), H_norm
!       DOUBLE PRECISION :: res_M_loc(3), res_rel_M, res_H(3), res_rel_H
!       DOUBLE PRECISION :: H_norm_prev, dH_norm, dH(3)

!       DOUBLE PRECISION, PARAMETER :: r0 = 1.0D-3
!       INTEGER, PARAMETER          :: maxi = 1000
!       DOUBLE PRECISION, PARAMETER :: dH_rel_max = 0.1

!       ! Get easy axis (ea) and off axis (oa); ea assumed parallel to remanent magnetization
!       i_tile = dom_proc(i)
!       M_rem_norm = NORM2(M_rem(:,state_dex(i_tile)))
!       mu_ea = constant_mu(  state_dex(i_tile))
!       mu_oa = constant_mu_o(state_dex(i_tile))
!       u_ea = M_rem(:,state_dex(i_tile))/M_rem_norm 
!       IF (u_ea(2).NE.0 .OR. u_ea(3).NE.0) THEN      ! x-product of u_ea with [1, 0, 0]; x-product of u_ea with u_oa_1
!           u_oa_1 = [0.d0, u_ea(3), -u_ea(2)]
!           u_oa_2 = [-u_ea(2)*u_ea(2) - u_ea(3)*u_ea(3), u_ea(1)*u_ea(2), u_ea(1)*u_ea(3)]
!       ELSE                                          ! x-product of u_ea with [0, 1, 0]; x-product of u_ea with u_oa_1
!           u_oa_1 = [-u_ea(3), 0.d0, u_ea(1)]
!           u_oa_2 = [u_ea(1)*u_ea(2), -u_ea(1)*u_ea(1) - u_ea(3)*u_ea(3), u_ea(2)*u_ea(3)]
!       END IF
!       u_oa_1 = u_oa_1/NORM2(u_oa_1)
!       u_oa_2 = u_oa_2/NORM2(u_oa_2)
      
!       ! Initialize picard
!       lambda_k = MIN(1/mu_ea, 1/mu_oa, 0.5)
!       H_next = H_back
!       H_new = H_next
!       M_targ_old = 0.0d0
!       DO iter = 1, maxi
!         ! Determine magnetization taking into account easy axis
!         M_targ = (M_rem_norm + (mu_ea-1)*DOT_PRODUCT(H_next,u_ea))*u_ea &
!                              + (mu_oa-1)*DOT_PRODUCT(H_next,u_oa_1)*u_oa_1 &
!                              + (mu_oa-1)*DOT_PRODUCT(H_next,u_oa_2)*u_oa_2
!         ! Update H-field
!         H_targ = H_back + MATMUL(N_store(:,:,1,i), M_targ)
!         res_H = H_targ - H_new
!         H_next = H_new + lambda_k * res_H
!         ! Residues
!         res_M_loc = M_targ - M_targ_old
!         res_rel_M = NORM2(res_M_loc)/MAX(NORM2(M_targ), small)
!         res_rel_H = NORM2(res_H    )/MAX(NORM2(H_targ), small)
!         ! Converged?
!         IF ((res_rel_H.LE.r0).AND.(res_rel_M.LE.r0)) EXIT
!         ! Update "old"
!         H_new = H_next
!         M_targ_old = M_targ
!       END DO
!       ! Cap change in H
!       lgood = .TRUE.
!       H_norm_prev = NORM2(H_prev)
!       dH = H_next-H_prev
!       dH_norm = NORM2(dH)
!       IF (H_norm_prev.GE.small) THEN
!         IF (dH_norm/H_norm_prev.GT.dH_rel_max) THEN
!           H_next = H_prev + (dH_rel_max*H_norm_prev)*(dH/dH_norm)
!           lgood = .FALSE.
!         END IF
!       END IF
!       ! Recalculate M
!       M_targ = (M_rem_norm + (mu_ea-1)*DOT_PRODUCT(H_next,u_ea))*u_ea &
!                            + (mu_oa-1)*DOT_PRODUCT(H_next,u_oa_1)*u_oa_1 &
!                            + (mu_oa-1)*DOT_PRODUCT(H_next,u_oa_2)*u_oa_2
!       H_prev = H_next
!       RETURN

!       END SUBROUTINE mumaterial_picard_rem

! !-----------------------------------------------------------------------
! ! mumaterial_picard_soft: Determines M_targ for given H_back. Limits relative
! ! change in H.
! !-----------------------------------------------------------------------
! ! param[in]: H_back. Input magnetic field
! ! param[inout]: H_prev. Last iteration's magnetic field. overwritten.
! ! param[in]: i. local index
! ! param[out]: M_targ. Target magnetization
! ! param[out]: lgood. True if change in H was not limited.
! !-----------------------------------------------------------------------
!       SUBROUTINE mumaterial_picard_soft(H_back, i, M_targ, lgood)
      
!       DOUBLE PRECISION, INTENT(in)  :: H_back(3)
!       INTEGER, INTENT(in)           :: i ! LOCAL index!
!       DOUBLE PRECISION, INTENT(out) :: M_targ(3)
!       LOGICAL, INTENT(out)          :: lgood

!       INTEGER :: iter, i_tile
!       DOUBLE PRECISION :: lambda_k

!       DOUBLE PRECISION :: H_next(3), H_new(3), H_targ(3), H_norm, M_targ_norm
!       DOUBLE PRECISION :: dH(3), dH_rel
!       DOUBLE PRECISION :: H_norm_prev, dH_norm
!       DOUBLE PRECISION :: mu, MAT(3,3)

!       INTEGER, PARAMETER          :: maxi = 1000
!       DOUBLE PRECISION, PARAMETER :: dH_rel_max = 0.1

!       lgood = .TRUE.
!       i_tile = dom_proc(i)
!       ASSOCIATE( sH => statefunc(state_dex(i_tile))%H, &
!                  sM => statefunc(state_dex(i_tile))%M, &
!                  m  => statefunc(state_dex(i_tile))%dMdH, &
!                  N_self => N_store(:,:,1,i))

!       ! Inital guess
!       CALL mumaterial_getstate_scalar_new(sH, sM, m, NORM2(H_back), M_targ_norm, mu) ! mu = chi
!       H_new = MATMUL(INVERT_3x3(I3 - mu*N_self),H_back)
!       !----------------------------------------
!       ! Iterative solver
!       !----------------------------------------
!       DO iter = 1, maxi
!         H_norm = NORM2(H_new)
!         CALL mumaterial_getstate_scalar_new(sH, sM, m, H_norm, M_targ_norm)
!         IF (H_norm .GT. 0.0d0) THEN
!           M_targ = M_targ_norm * H_new / H_norm
!           lambda_k = MIN(H_norm/M_targ_norm, 0.5d0)
!         ELSE
!           M_targ = 0.0d0
!           lambda_k = 0.5d0
!         END IF
!         ! Update H-field
!         H_targ = H_back + MATMUL(N_self, M_targ)
!         dH = H_targ - H_new
!         H_next = H_new + lambda_k * dH
!         ! Residues
!         dH_rel = NORM2(dH)/MAX(NORM2(H_targ), small)
!         ! Converged?
!         IF ((dH_rel.LE.threshold).OR.(NORM2(H_new).LT.small)) EXIT
!         ! Update "old"
!         H_new = H_next
!       END DO
!       !----------------------------------------
!       ! Solver done; cap H
!       !----------------------------------------
!       H_norm_prev = NORM2(H_prev)
!       dH = H_next-H_prev
!       dH_norm = NORM2(dH)
!       IF (H_norm_prev.GE.small) THEN
!         IF (dH_norm/H_norm_prev.GT.dH_rel_max) THEN
!           H_next = H_prev + (dH_rel_max*H_norm_prev)*(dH/dH_norm)
!           lgood = .FALSE.
!         END IF
!       END IF
!       !----------------------------------------
!       ! Recalculate M
!       !----------------------------------------
!       H_norm = NORM2(H_next)
!       CALL mumaterial_getstate_scalar_new(sH, sM, m, H_norm, M_targ_norm)
!       END ASSOCIATE
!       IF (H_norm .GT. small) THEN
!             M_targ = M_targ_norm * H_next / H_norm
!       ELSE
!             M_targ = 0
!       END IF
      
!       IF (iter.GE.maxi)  WRITE(6,'(A,I0,A,ES12.4,A,ES12.4,A)') "[",i_tile,"] H-solve exceed maxi (H=",H_norm, " res=",dH_rel,")"
!       H_prev = H_next
!       RETURN

!       END SUBROUTINE mumaterial_picard_soft

!-----------------------------------------------------------------------
! mumaterial_NR_soft: Newton Rhapson solver for finding self-consistent
! H and M as a function of the total field without self contribution.
!-----------------------------------------------------------------------
! param[in]: H_back. Input magnetic field
! param[inout]: H_prev. Last iteration's magnetic field. overwritten.
! param[in]: i. local index
! param[out]: M_targ. Target magnetization
! param[out]: lgood. True if change in H was not limited.
!-----------------------------------------------------------------------
      ! SUBROUTINE mumaterial_NR_soft(H_back, i, M_targ, lgood)

      ! IMPLICIT NONE

      ! DOUBLE PRECISION, INTENT(in)    :: H_back(3)
      ! INTEGER, INTENT(in)             :: i 
      ! DOUBLE PRECISION, INTENT(out)   :: M_targ(3)
      ! LOGICAL, INTENT(out)            :: lgood

      ! INTEGER :: iter, i_tile
      ! DOUBLE PRECISION :: H_new(3), H_norm
      ! DOUBLE PRECISION :: F(3), r, dMdH_mat(3,3), J(3,3), dN(3), J_reg(3,3)
      ! DOUBLE PRECISION :: alpha,phi,gTd,grad_phi(3)
      ! DOUBLE PRECISION :: deltaH(3), H_norm_prev, dH_norm, dH(3)
      ! DOUBLE PRECISION :: dMdH, M_norm,r0,H_trial(3),F_trial(3)
      ! DOUBLE PRECISION :: phi_trial, M_trial(3)
      ! DOUBLE PRECISION :: phi_best, lambda, MAT(3,3)

      ! DOUBLE PRECISION, PARAMETER :: c = 1.0D-4, rho = 0.50d0, invrho = 2.0d0, alpha_min=1.0D-2, alpha_max = 0.5d0
      ! DOUBLE PRECISION, PARAMETER :: F0 = 1.0D-3
      ! INTEGER, PARAMETER          :: maxi = 1000
      ! DOUBLE PRECISION, PARAMETER :: dH_rel_max = 0.1, lambda_max = 10.0d0, lambda_min = 1.0d-12

      ! i_tile = dom_proc(i)
      ! ASSOCIATE( sH => statefunc(state_dex(i_tile))%H, &
      !            sM => statefunc(state_dex(i_tile))%M, &
      !            m  => statefunc(state_dex(i_tile))%dMdH, &
      !            N  => N_store(:,:,1,i))
      ! H_new = H_back
      ! r0 = threshold
      ! alpha = alpha_max
      ! lambda = 1.0d-3
      ! lgood = .TRUE.
      ! !----------------------------------------
      ! ! Linear step initially
      ! !----------------------------------------
      ! CALL mumaterial_getstate_vector(sH, sM, m, H_new, M_targ, dMdH) ! Get dMdH    
      ! MAT = I3-dMdH*N
      ! H_new = SOLVE_3x3(MAT, H_back, "LINSTEP")

      ! !--------------------------------------
      ! ! Begin iterations
      ! !-------------------------------------- 
      ! DO iter = 1, maxi
      !   !--------------------------------------
      !   ! Compute Jacobian
      !   !-------------------------------------- 
      !   H_norm = NORM2(H_new)
      !   CALL mumaterial_getstate_vector(sH, sM, m, H_new, M_targ, dMdH) ! Get dMdH    
      !   M_norm = NORM2(M_targ)
      !   IF (H_norm.GT.small) THEN
      !     dMdH_mat = (dMdH-M_norm/H_norm)*OUTER_PRODUCT(H_new,H_new)/(H_norm*H_norm)+M_norm/H_norm*I3
      !   ELSE
      !     dMdH_mat = dMdH*I3
      !   END IF
      !   F = H_new - H_back - MATMUL(N,M_targ)
      !   J = (I3 - MATMUL(N,dMdH_mat)) ! Jacobian
      !   !--------------------------------------
      !   ! Calculate with regularized Jacobian
      !   !--------------------------------------
      !   alpha = alpha_max
      !   DO
      !     J_reg = J + MAX(lambda*SQRT(SUM(J**2)),lambda_min)*I3
      !     phi_best = 0.5d0*DOT_PRODUCT(F,F)
      !     !--------------------------------------
      !     ! Newton descent or gradient descent
      !     !--------------------------------------
      !     grad_phi = MATMUL(TRANSPOSE(J_reg),F)
      !     deltaH = -SOLVE_3x3(J_reg, F, "NEWTON")   ! invJ*F
      !     !--------------------------------------
      !     ! Line search
      !     !--------------------------------------
      !     gTd = DOT_PRODUCT(grad_phi, deltaH)
      !     DO  ! Reduce alpha until phi is good
      !       H_trial = H_new+alpha*deltaH
      !       CALL mumaterial_getstate_vector(sH, sM, m, H_trial, M_trial)
      !       F_trial = H_trial - H_back - MATMUL(N, M_trial)
      !       phi = 0.5d0*DOT_PRODUCT(F_trial,F_trial)
      !       IF (phi.LE.phi_best+c*alpha*gTd .OR. alpha.EQ.alpha_min) EXIT
      !       alpha = MAX(rho*alpha,alpha_min)        
      !     END DO

      !     IF (alpha.EQ.alpha_min) THEN
      !       IF (lambda.EQ.lambda_max) THEN
      !         H_new = H_new - 0.01d0*F
      !         EXIT
      !       END IF
      !       lambda = MIN(lambda * 10.d0,lambda_max)
      !       alpha = alpha_max
      !       CYCLE
      !     END IF

      !     IF (alpha.EQ.alpha_max) lambda = lambda*0.5d0
      !     H_new = H_trial
      !     F = F_trial
      !     EXIT
          
      !   END DO

      !   !--------------------------------------
      !   ! Convergence check
      !   !-------------------------------------- 
      !   H_norm = NORM2(H_new)
      !   r = NORM2(F)/max(H_norm,small)
      !   IF (r.LT.r0.OR.NORM2(F).LT.F0) EXIT   
      ! END DO
      ! IF (iter.GE.maxi) WRITE(6,'(A,I0,A,ES10.2,A,ES10.2,A)') "[",i_tile,"] H-solve exceed maxi (H=",H_norm, " res=",r,")"

      ! !----------------------------------------
      ! ! Solver done; cap H
      ! !----------------------------------------
      ! H_norm_prev = NORM2(H_prev)
      ! dH = H_new-H_prev
      ! dH_norm = NORM2(dH)
      ! IF (H_norm_prev.GE.small) THEN
      !   IF (dH_norm/H_norm_prev.GT.dH_rel_max) THEN
      !     H_new = H_prev + (dH_rel_max*H_norm_prev)*(dH/dH_norm)
      !     lgood = .FALSE.
      !   END IF
      ! END IF
      ! !----------------------------------------
      ! ! Recalculate M
      ! !----------------------------------------
      ! CALL mumaterial_getstate_vector(sH, sM, m, H_new, M_targ)
      ! H_prev = H_new
      ! END ASSOCIATE

      ! RETURN

      ! END SUBROUTINE mumaterial_NR_soft




!-----------------------------------------------------------------------
! mumaterial_block_solve: Block NR solver for the H-field of a single cluster.
! Minimizes the residual given by:
! 
!     R = H_ext + N.M - H
!   H_ext: (3,ntet_proc) H-field due to currents (Happ), clusters, dipoles,
!   and neighbors IN OTHER CLUSTERS
!   N: (3,3,ntet_proc,ntet_proc) demagnetization tensor
!   M: (3,ntet_proc) magnetizations
!   H: (3,ntet_proc) total H-field at every element
!
!   
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_block_solve

      IMPLICIT NONE
      INTEGER :: i, i_tile, i1, i2, j, j1, i_nr, j_tile, ik, itet, n
      DOUBLE PRECISION :: dMdH, M_norm, H_norm, res_rel, chi
      DOUBLE PRECISION ::  M_soft(3,ntet_soft_max), dMdH_soft(ntet_soft_max)
      INTEGER :: INFO, IPIV(3*ntet_proc)
      INTEGER, PARAMETER :: maxi_nr = 10
      DOUBLE PRECISION, PARAMETER :: tol_nr2 = 1.0d-8
      INTEGER :: stype, sdex
      LOGICAL :: is_soft(ntet_proc)
      DOUBLE PRECISION :: N11, N12, N13, N21, N22, N23, N31, N32, N33
      DOUBLE PRECISION :: M1, M2, M3, R1, R2, R3, H1, H2, H3
      DOUBLE PRECISION :: dMdH11, dMdH12, dMdH13, dMdH21, dMdH22, dMdH23, dMdH31, dMdH32, dMdH33

      
      ! NR iterations until converged
      DO i_nr = 1, maxi_nr

        !---------------------------
        ! First determine H, M, dMdH
        !---------------------------
        ! Linear materials
        DO ik = 1, nlinear
          sdex = sdex_linear(ik)
          chi = constant_mu(sdex)-1.0d0
          DO i = 1, ntet_linear(ik)
            itet = dom_linear(i,ik)
            M_spline(:,itet) = chi*H_prev(:,itet)
            dMdH_spline(itet) = chi
          END DO
        END DO
        ! Soft materials
        DO ik = 1, nsoft
          sdex = sdex_soft(ik)
          n = ntet_soft(ik)
          CALL mumaterial_getstate_vector_batch(&
            statefunc(sdex)%H, &
            statefunc(sdex)%M, &
            statefunc(sdex)%dMdH, &
            H_prev(:,dom_soft(1:n,ik)), &
            M_soft(:,1:n), n,  &
            dMdH_soft(1:n))

          DO itet = 1, n
            i = dom_soft(itet,ik)
            M_spline(:,i) = M_soft(:,itet)
            dMdH_spline(i) = dMdH_soft(itet)
          END DO
        END DO
        ! Hard materials (TODO)
        DO ik = 1, nhard
        END DO

        !---------------------------
        ! Get norms
        !---------------------------
        DO i = 1, ntet_proc
          H_norm_proc(i) = NORM2(H_prev(:,i))
          M_norm_proc(i) = NORM2(M_spline(:,i))
        END DO
        !---------------------------
        ! Get dMdH_mat
        !---------------------------
        ! Linear
        DO ik = 1, nlinear
          sdex = sdex_linear(ik)
          DO i = 1, ntet_linear(ik)
            itet = dom_linear(i,ik)
            dMdH_mat(:,:,itet) = dMdH_spline(itet)*I3
          END DO
        END DO
        ! Soft
        DO ik = 1, nsoft
          sdex = sdex_soft(ik)
          DO i = 1, ntet_soft(ik)
            itet = dom_soft(i,ik)
            M_norm = M_norm_proc(itet)
            H_norm = H_norm_proc(itet)
            dMdH = dMdH_spline(itet)
            IF (H_norm.GT.small) THEN
              dMdH_mat(:,:,itet) = (dMdH-M_norm/H_norm)*OUTER_PRODUCT(H_prev(:,itet),H_prev(:,itet))/(H_norm*H_norm)+M_norm/H_norm*I3
            ELSE
              dMdH_mat(:,:,itet) = dMdH*I3
            END IF
          END DO
        END DO
        ! Hard (TODO)
        DO ik = 1, nhard
        END DO
        
        !---------------------------
        ! Build system of equations
        !---------------------------
        DO i = 1, ntet_proc
          i1 = 3*(i-1)+1
          H_block(i1)   = H_prev(1,i)
          H_block(i1+1) = H_prev(2,i)
          H_block(i1+2) = H_prev(3,i)
          R1 = H_ext(1,i) - H_prev(1,i)
          R2 = H_ext(2,i) - H_prev(2,i)
          R3 = H_ext(3,i) - H_prev(3,i)
          DO j = 1, ntet_proc
            ! Load stuff into cache first
            N11 = N_block(1,1,j,i); N12 = N_block(1,2,j,i); N13 = N_block(1,3,j,i)
            N21 = N_block(2,1,j,i); N22 = N_block(2,2,j,i); N23 = N_block(2,3,j,i)
            N31 = N_block(3,1,j,i); N32 = N_block(3,2,j,i); N33 = N_block(3,3,j,i)
            M1 = M_spline(1,j); M2 = M_spline(2,j); M3 = M_spline(3,j)
            dMdH11 = dMdH_mat(1,1,j); dMdH12 = dMdH_mat(1,2,j); dMdH13 = dMdH_mat(1,3,j)
            dMdH21 = dMdH_mat(2,1,j); dMdH22 = dMdH_mat(2,2,j); dMdH23 = dMdH_mat(2,3,j)
            dMdH31 = dMdH_mat(3,1,j); dMdH32 = dMdH_mat(3,2,j); dMdH33 = dMdH_mat(3,3,j)
            ! Residual
            R1 = R1 + N11*M1 + N12*M2 + N13*M3
            R2 = R2 + N21*M1 + N22*M2 + N23*M3
            R3 = R3 + N31*M1 + N32*M2 + N33*M3
            ! Jacobian
            j1 = 3*(j-1)+1
            J_block(i1  ,j1  ) =  N11*dMdH11 + N12*dMdH21 + N13*dMdH31
            J_block(i1  ,j1+1) =  N11*dMdH12 + N12*dMdH22 + N13*dMdH32
            J_block(i1  ,j1+2) =  N11*dMdH13 + N12*dMdH23 + N13*dMdH33

            J_block(i1+1,j1  ) =  N21*dMdH11 + N22*dMdH21 + N23*dMdH31
            J_block(i1+1,j1+1) =  N21*dMdH12 + N22*dMdH22 + N23*dMdH32
            J_block(i1+1,j1+2) =  N21*dMdH13 + N22*dMdH23 + N23*dMdH33

            J_block(i1+2,j1  ) =  N31*dMdH11 + N32*dMdH21 + N33*dMdH31
            J_block(i1+2,j1+1) =  N31*dMdH12 + N32*dMdH22 + N33*dMdH32
            J_block(i1+2,j1+2) =  N31*dMdH13 + N32*dMdH23 + N33*dMdH33
          END DO
          R_block(i1)   = -R1
          R_block(i1+1) = -R2
          R_block(i1+2) = -R3
        END DO
        ! Diagonal terms
        DO i = 1, ntet_proc
          i1 = 3*(i-1)+1; i2 = 3*i
          J_block(i1:i2,i1:i2) = J_block(i1:i2,i1:i2) - I3
        END DO

        ! LAPACK solve
        CALL DGESV(3*ntet_proc, 1, J_block, 3*ntet_proc, IPIV, R_block, 3*ntet_proc, INFO)
        IF (INFO /= 0) THEN
          WRITE(6,*) "DGESV failed, INFO=", INFO
          CALL mumaterial_abort()
        END IF

        ! Update H
        res_rel = 0.0d0 
        DO i = 1, ntet_proc
          i1 = 3*(i-1)+1
          R1 = R_block(i1)
          R2 = R_block(i1+1)
          R3 = R_block(i1+2)
          H1 = H_block(i1)
          H2 = H_block(i1+1)
          H3 = H_block(i1+2)
          H_prev(1,i) = H1 + R1
          H_prev(2,i) = H2 + R2
          H_prev(3,i) = H3 + R3
          res_rel = MAX(res_rel, (R1*R1+R2*R2+R3*R3)/MAX(H1*H1+H2*H2+H3*H3,small))
        END DO

        ! Convergence check
        IF (res_rel.LT.tol_nr2) EXIT
      END DO

      RETURN

      END SUBROUTINE mumaterial_block_solve

!-----------------------------------------------------------------------
! mumaterial_cluster_setup: Calculates position, diameter, initial dipole
! and quadrupole moment of clusters. m_cluster and Q_cluster are only
! necessary on shar_rank = 0, r_cluster and d_cluster on all ranks.
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_cluster_setup
      IMPLICIT NONE

      INTEGER :: wr_dex, i, j, tile1, tile2, i_tile, j_tile, tiles(2), rank_worst
      DOUBLE PRECISION :: pair_in(2), pair_out(2), d_ij, d_loc, d_global
      INTEGER, PARAMETER :: TAG_TILES = 201
      
      IF (lcomm) THEN
#if defined(MPI_OPT)
        !-----------------------------------------
        ! dom_cluster: Array of which cluster gets which elements: ~Ntet + number of 0 pads
        !-----------------------------------------    
        ALLOCATE(dom_cluster(ntet_proc_max,world_size),&
                 r_cluster(3,world_size),d_cluster(world_size))
        dom_cluster = 0
        wr_dex = world_rank+1 ! Pesky 0-based indexing
        dom_cluster(1:ntet_proc, wr_dex) = dom_proc(1:ntet_proc)
        ! Cluster position and diameter
        r_cluster = 0.0d0; d_cluster = -0.0d0; 
        r_cluster(:, wr_dex) = SUM(tet_cen_proc,DIM=2)/ntet_proc 
        d_cluster(wr_dex) = SQRT( SUM(NORM2(tet_cen_proc-SPREAD(r_cluster(:,wr_dex), DIM=2, NCOPIES=ntet_proc),DIM=1)**2) / ntet_proc)
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, dom_cluster, ntet_proc_max*world_size, MPI_INTEGER, MPI_SUM, comm_world, ierr_mpi )
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, r_cluster,               3*world_size, MPI_DOUBLE_PRECISION, MPI_SUM, comm_world, ierr_mpi )
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, d_cluster,                 world_size, MPI_DOUBLE_PRECISION, MPI_SUM, comm_world, ierr_mpi )

        !-----------------------------------------
        ! m_cluster, Q_cluster: arrays for dipole and quadrupole tensors
        !-----------------------------------------   
        ALLOCATE(m_cluster(3,world_size),Q_cluster(6,world_size))

        !------------------------------------------------------------
        ! Find worst element pair; only for printing to screen
        tile1 = -1; tile2 = -1; d_loc = -1.0d0
        DO i = 1, ntet_proc-1
          DO j = i+1, ntet_proc
            d_ij = NORM2(tet_cen_proc(:,j)-tet_cen_proc(:,i))
            IF (d_ij.GT.d_loc) THEN
              tile1 = dom_proc(i); tile2 = dom_proc(j); d_loc = d_ij
            END IF
          END DO
        END DO

        ! Send to master
        pair_in = [d_loc, DBLE(world_rank)]
        CALL MPI_ALLREDUCE(pair_in,pair_out,1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,comm_world,ierr_mpi)       
        rank_worst = NINT(pair_out(2)); d_global = pair_out(1)
        tiles = [tile1, tile2]
        IF (lismaster .AND. (world_rank.NE.rank_worst)) THEN
          CALL MPI_RECV(tiles,2,MPI_INTEGER, rank_worst,TAG_TILES, comm_world, mstat, ierr_mpi)
        ELSE
          IF (world_rank.EQ.rank_worst) THEN
            CALL MPI_SEND(tiles,2,MPI_INTEGER, 0, TAG_TILES, comm_world, ierr_mpi)
          END IF 
        END IF

        IF (lverb) THEN
          WRITE(6,'(3X,A,ES9.2,A,/,17X,ES9.2,A)')   'Cluster diam.:',MINVAL(d_cluster),' m,',MAXVAL(d_cluster),' m'
          WRITE(6,'(3X,A,I0,A,I0,A,/,18X,A,ES9.2,A)') 'Worst pair   : [',tiles(1), ',', tiles(2), ']','(', d_global, ' m)'
          FLUSH(6)         
        END IF
        !------------------------------------------------------------
#endif
      END IF
      RETURN

      END SUBROUTINE mumaterial_cluster_setup

!-----------------------------------------------------------------------
! mumaterial_cluster_update: Recalculates Q tensor and moment
!-----------------------------------------------------------------------
! param[inout] R_in, mom_in. Helper 2-D arrays that are allocated elsewhere.
!               passed as inout to avoid repeated de/allocation.
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_cluster_update(R_in,mom_in)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: R_in, mom_in
      INTEGER :: i, j, k
      DOUBLE PRECISION :: A, B, C, R(3,ntet_proc_max), mom(3,ntet_proc_max), Q4, Q5, Q6, m1, m2, m3

      IF (shar_rank.EQ.master) THEN
        m_cluster = 0.0d0
        Q_cluster = 0.0d0
        DO j = 1, shar_size
          i = shar_ranks_arr(j) + 1 ! MPI 0-index
          ASSOCIATE(  sz  => dom_sizes(i), &
                      dom => dom_cluster(1:sz,i))
            DO k = 1, sz
              mom(:,k) = M_all(:,map_g2all(dom(k))) * tet_vol(dom(k))
              R(:,k) = tet_cen(:,dom(k)) - r_cluster(:,i)
            END DO
            
            m1= 0.0d0; m2= 0.0d0; m3= 0.0d0
            A = 0.0d0; B = 0.0d0; C = 0.0d0
            Q4= 0.0d0; Q5= 0.0d0; Q6= 0.0d0
            DO k = 1, sz
              m1 = m1 + mom(1,k)
              m2 = m2 + mom(2,k)
              m3 = m3 + mom(3,k)
              A = A + R(1,k)*mom(1,k)
              B = B + R(2,k)*mom(2,k)
              C = C + R(3,k)*mom(3,k)
              Q4 = Q4 + R(1,k)*mom(2,k) + R(2,k)*mom(1,k)
              Q5 = Q5 + R(1,k)*mom(3,k) + R(3,k)*mom(1,k)
              Q6 = Q6 + R(2,k)*mom(3,k) + R(3,k)*mom(2,k)
            END DO
            ! Dipole tensor terms
            m_cluster(1,i) = m1
            m_cluster(2,i) = m2
            m_cluster(3,i) = m3
            ! Quadrupole tensor terms
            Q_cluster(1,i) = 2.0d0/3.0d0*(2*A-B-C) ! xx
            Q_cluster(2,i) = 2.0d0/3.0d0*(2*B-A-C) ! yy
            Q_cluster(3,i) = 2.0d0/3.0d0*(2*C-A-B) ! zz, traceless
            Q_cluster(4,i) = Q4 ! xy
            Q_cluster(5,i) = Q5 ! xz
            Q_cluster(6,i) = Q6 ! yz
          END ASSOCIATE

        END DO
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, m_cluster, 3*world_size, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, Q_cluster, 6*world_size, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi)
      END IF
      CALL MPI_Bcast(m_cluster,3*world_size,MPI_DOUBLE_PRECISION,0,comm_shar,ierr_mpi)
      CALL MPI_Bcast(Q_cluster,6*world_size,MPI_DOUBLE_PRECISION,0,comm_shar,ierr_mpi)

      END SUBROUTINE mumaterial_cluster_update

      
!-----------------------------------------------------------------------
! mumaterial_applied_field: Overwrites H_ext with H_app.
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_applied_field()
      IMPLICIT NONE

      H_ext = H_app ! Static field using slice

      END SUBROUTINE mumaterial_applied_field

!-----------------------------------------------------------------------
! mumaterial_cluster_field: Adds cluster contributions to H_ext
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_cluster_field()

      IMPLICIT NONE

      INTEGER :: i
      DOUBLE PRECISION :: r_tile(3), hx, hy, hz
      !-----------------------------
      ! Update cluster tensors
      !-----------------------------
      CALL mumaterial_cluster_update(R_quad,mom_quad)
      !-----------------------------
      ! Iterate over elements
      !-----------------------------
      DO i = 1, ntet_proc
        r_tile = tet_cen_proc(:,i)
        hx = H_ext(1,i); hy = H_ext(2,i); hz = H_ext(3,i)
        CALL mumaterial_cluster_field_calc(r_tile(1), r_tile(2), r_tile(3), hx, hy, hz)
        H_ext(1,i) = hx; H_ext(2,i) = hy; H_ext(3,i) = hz
      END DO
      RETURN
      END SUBROUTINE mumaterial_cluster_field

      
!-----------------------------------------------------------------------
! mumaterial_cluster_field_calc: Calculates cluster field at a point
!-----------------------------------------------------------------------
! param[in]: x_in. x-coordinate
! param[in]: y_in. y-coordinate
! param[in]: z_in. z-coordinate
! param[inout]: hx. x-component of H_cluster (additive)
! param[inout]: hy. y-component of H_cluster (additive)
! param[inout]: hz. z-component of H_cluster (additive)
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_cluster_field_calc(x_in, y_in, z_in, hx, hy, hz)
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in) :: x_in, y_in, z_in
      DOUBLE PRECISION, INTENT(inout) :: hx, hy, hz
        
      INTEGER :: i
      DOUBLE PRECISION :: x, y, z, r2, m_x, m_y, m_z, mr, invr5, invr2
      DOUBLE PRECISION :: c1, c2, rq_x, rq_y, rq_z, rqr
      !-----------------------------
      ! Iterate over clusters
      !-----------------------------  
      DO i = 1, world_size
        IF (.NOT.iscluster_proc(i)) CYCLE
        x = x_in - r_cluster(1,i);  y = y_in - r_cluster(2,i); z = z_in - r_cluster(3,i)
        r2 = x*x + y*y + z*z; invr2 = 1.0d0/r2; invr5 = invr2 * invr2 * SQRT(invr2); 
        c1 = INV4PI*invr5; c2 = INV8PI*invr5*invr2 ! 1/(8pi.r^7)
        ! Dipole contribution 
        m_x = M_cluster(1,i); m_y = M_cluster(2,i); m_z = M_cluster(3,i)
        mr = m_x*x + m_y*y + m_z*z
        ! Quadrupole contribution: 1/(8*pi*r^7) * ((r.Q r)r - 2*Q.r*r^2)
        rq_x = Q_cluster(1,i)*x + Q_cluster(4,i)*y + Q_cluster(5,i)*z
        rq_y = Q_cluster(4,i)*x + Q_cluster(2,i)*y + Q_cluster(6,i)*z
        rq_z = Q_cluster(5,i)*x + Q_cluster(6,i)*y + Q_cluster(3,i)*z
        rqr = rq_x*x + rq_y*y +rq_z*z
        hx = hx + c1*(3.0d0*mr*x - r2*m_x) + c2*(5.0d0*rqr*x - 2*rq_x*r2)
        hy = hy + c1*(3.0d0*mr*y - r2*m_y) + c2*(5.0d0*rqr*y - 2*rq_y*r2)
        hz = hz + c1*(3.0d0*mr*z - r2*m_z) + c2*(5.0d0*rqr*z - 2*rq_z*r2)
      END DO
      RETURN
      END SUBROUTINE mumaterial_cluster_field_calc

!-----------------------------------------------------------------------
! mumaterial_dipole_field: Adds dipole contributions to H_ext
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_dipole_field(lfull)

      IMPLICIT NONE
      INTEGER :: i, j, k, i_tile, n_dip
      DOUBLE PRECISION :: r_tile(3), rix, riy, riz, dMx, dMy, dMz, r2
      DOUBLE PRECISION :: hx, hy, hz, M(3)
      LOGICAL, INTENT(in) :: lfull
    
      !-----------------------------
      ! First check for stale elements
      !-----------------------------  
      is_stale = .FALSE.
      IF (.NOT.lfull) THEN
        n_dip = 0
        r2 = 0.25d0*threshold**2 ! (threshold / 2)^2
        DO j = 1, ntet_partial_proc
          dMx = M_partial(1,j) - M_snapshot(1,j) 
          dMy = M_partial(2,j) - M_snapshot(2,j) 
          dMz = M_partial(3,j) - M_snapshot(3,j)
          is_stale(j) = (dMx*dMx+dMy*dMy+dMz*dMz.GT.(r2*DOT_PRODUCT(M_snapshot(:,j),M_snapshot(:,j)))) 
          IF (is_stale(j)) THEN
            n_dip = n_dip + 1
            stale_dex(n_dip) = j
            M_dip(:,n_dip) = M_partial(:,j) - M_snapshot(:,j)
            r_dip(:,n_dip) = r_partial(:,j)
            V_dip(n_dip) = V_partial(j)
          END IF         
        END DO
        IF (n_dip.EQ.0) RETURN
      END IF
      !-----------------------------
      ! Now iterate over ntet_proc
      !-----------------------------  
      DO i = 1, ntet_proc
        i_tile = dom_proc(i)
        r_tile = tet_cen(:,i_tile)
        rix = r_tile(1); riy = r_tile(2); riz = r_tile(3)
        !-----------------------------
        ! Calculate everything; 
        ! If full recompute, realculate H_mid from scratch
        ! If not, calculate change in H_mid and add (linear in dM)
        !-----------------------------
        IF (lfull) THEN
          hx = 0.0d0; hy = 0.0d0; hz = 0.0d0
          CALL mumaterial_dipole_field_calc(rix, riy, riz, ntet_partial_proc, r_partial, V_partial, M_partial, hx, hy, hz)
        ELSE
          hx = H_mid(1,i); hy = H_mid(2,i); hz = H_mid(3,i)
          CALL mumaterial_dipole_field_calc(rix, riy, riz, n_dip, r_dip, V_dip, M_dip, hx, hy, hz)
        END IF

        ! Subtract neighbors
        DO k = 1, nbrs_count(i)
          j = nbrs_proc(k,i)
          IF (.NOT.(lfull.OR.is_stale(j))) CYCLE
          IF (lfull) THEN
            M = -M_partial(:,j) ! Added full contribution before
          ELSE
            M = -(M_partial(:,j)-M_snapshot(:,j)) ! Added correction before
          END IF
          CALL mumaterial_dipole_field_calc_single(rix, riy, riz, r_partial(:,j), V_partial(j), M, hx, hy, hz)
        END DO
        ! Subtract block elements
        DO k = 1, ntet_proc
          j = map_g2p(dom_proc(k))
          IF (.NOT.(lfull.OR.is_stale(j))) CYCLE
          IF (lfull) THEN
            M = -M_partial(:,j) ! Added full contribution before
          ELSE
            M = -(M_partial(:,j)-M_snapshot(:,j)) ! Added correction before
          END IF
          CALL mumaterial_dipole_field_calc_single(rix, riy, riz, r_partial(:,j), V_partial(j), M, hx, hy, hz)
        END DO
        H_mid(1,i) = hx; H_mid(2,i) = hy; H_mid(3,i) = hz
      END DO  
      !-----------------------------
      ! Update snapshots
      !-----------------------------  
      DO j = 1, ntet_partial_proc
        IF (is_stale(j).OR.lfull) M_snapshot(:,j) = M_partial(:,j)
      END DO

      RETURN
      END SUBROUTINE mumaterial_dipole_field

!-----------------------------------------------------------------------
! mumaterial_dipole_field_calc_full: Calculates dipole field at a point due to
! magnetized elements passed in.
!-----------------------------------------------------------------------
! param[in]: x_in. x-coordinate
! param[in]: y_in. y-coordinate
! param[in]: z_in. z-coordinate
! param[in]: n. number of elements
! param[in]: r_arr. coordinates of the elements
! param[in]: V_arr. volumes of the elements
! param[in]: M_arr. magnetizations of the elements
! param[inout]: hx. x-component of H_dipole (additive)
! param[inout]: hy. y-component of H_dipole (additive)
! param[inout]: hz. z-component of H_dipole (additive)
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_dipole_field_calc(x_in, y_in, z_in, n, r_arr, V_arr, M_arr, hx, hy, hz)
      IMPLICIT NONE

      INTEGER, INTENT(in) :: n
      DOUBLE PRECISION, INTENT(in) :: x_in, y_in, z_in, r_arr(3,n), V_arr(n), M_arr(3,n)
      DOUBLE PRECISION, INTENT(inout) :: hx, hy, hz
        
      INTEGER :: i
      DOUBLE PRECISION :: x, y, z, r2, mr, MV_x, MV_y, MV_z, rix, riy, riz, invr5, invr2
      DOUBLE PRECISION :: V, c, dMx, dMy, dMz

      DO i = 1, n
        V = V_arr(i)
        x = x_in - r_arr(1,i);  y = y_in - r_arr(2,i); z = z_in - r_arr(3,i)
        r2 = x*x + y*y + z*z + small; invr2 = 1.0d0/r2; invr5 = invr2 * invr2 * SQRT(invr2); c = INV4PI*invr5
        MV_x = M_arr(1,i)*V;MV_y = M_arr(2,i)*V;MV_z = M_arr(3,i)*V
        mr = MV_x*x+MV_y*y+MV_z*z
        hx = hx + c*(3.0d0*mr*x-r2*MV_x)
        hy = hy + c*(3.0d0*mr*y-r2*MV_y)
        hz = hz + c*(3.0d0*mr*z-r2*MV_z)
      END DO
      RETURN

      END SUBROUTINE mumaterial_dipole_field_calc

!-----------------------------------------------------------------------
! mumaterial_dipole_field_calc_single: Calculates dipole field at a point due to
! a single magnetized element.
!-----------------------------------------------------------------------
! param[in]: x_in. x-coordinate
! param[in]: y_in. y-coordinate
! param[in]: z_in. z-coordinate
! param[in]: r. coordinates of element
! param[in]: V. volume of element
! param[in]: M. magnetization of element
! param[inout]: hx. x-component of H_dipole (additive)
! param[inout]: hy. y-component of H_dipole (additive)
! param[inout]: hz. z-component of H_dipole (additive)
!-----------------------------------------------------------------------
    SUBROUTINE mumaterial_dipole_field_calc_single(x_in, y_in, z_in, r, V, M, hx, hy, hz)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(in) :: x_in, y_in, z_in, r(3), V, M(3)
    DOUBLE PRECISION, INTENT(inout) :: hx, hy, hz
      
    DOUBLE PRECISION :: x, y, z, r2, mr, MV_x, MV_y, MV_z, rix, riy, riz, invr5, invr2
    DOUBLE PRECISION :: c, dMx, dMy, dMz

    x = x_in - r(1);  y = y_in - r(2); z = z_in - r(3)
    r2 = x*x + y*y + z*z + small; invr2 = 1.0d0/r2; invr5 = invr2 * invr2 * SQRT(invr2); c = INV4PI*invr5
    MV_x = M(1)*V; MV_y = M(2)*V; MV_z = M(3)*V
    mr = MV_x*x+MV_y*y+MV_z*z
    hx = hx + c*(3.0d0*mr*x-r2*MV_x)
    hy = hy + c*(3.0d0*mr*y-r2*MV_y)
    hz = hz + c*(3.0d0*mr*z-r2*MV_z)

    RETURN

    END SUBROUTINE mumaterial_dipole_field_calc_single

!-----------------------------------------------------------------------
! mumaterial_getstate_slopes: Gets spline slopes by solving a tridiagonal matrix
!-----------------------------------------------------------------------
! param[in]:  fx. x-coordinates of function to be interpolated
! param[in]:  fy. y-values of function to be interpolated
! param[out]: m. slopes
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_getstate_slopes(fx, fy, m)
        
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: fx, fy
      DOUBLE PRECISION, DIMENSION(:), INTENT(out) :: m
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: h, delta, dl, d, du
      DOUBLE PRECISION :: alpha, beta, s2,  tau
      INTEGER :: n, i, INFO

      n = SIZE(fx)
      ALLOCATE(h(n-1),delta(n-1),dl(n-1),d(n),du(n-1))

      DO i = 1, n-1
        h(i) = fx(i+1) - fx(i)
        delta(i) = (fy(i+1)-fy(i))/h(i)
      END DO

      ! Lower diagonal
      dl(1) = 1.0d0
      DO i = 2, n-1
        dl(i) = h(i)
      END DO

      ! Main diagonal
      d(1) = 2.0d0
      DO i = 2, n-1
        d(i) = 2.0d0 * (h(i-1)+h(i))
      END DO
      d(n) = 2.0d0

      ! Upper diagonal
      DO i = 2, n-1
        du(i) = h(i-1)
      END DO
      du(n-1) = 1.0d0

      ! RHS
      m(1) = 3.0d0*delta(1)
      DO i = 2, n-1
        m(i) = 3.0d0*(h(i)*delta(i-1)+h(i-1)*delta(i))
      END DO
      m(n) = 3.0d0*delta(n-1)

      ! LAPACK
      CALL DGTSV(n, 1, dl, d, du, m, n, INFO)

      
      ! Monotonicity correction (Fritsch-Carlson)
      DO i = 1, n-1
        IF (ABS(delta(i)) .LT. TINY(1.0d0)) THEN
            m(i) = 0.0d0
            m(i+1) = 0.0d0
        ELSE
            alpha = m(i) / delta(i)
            beta  = m(i+1) / delta(i)
            s2 = alpha**2 + beta**2
            IF (s2 .GT. 9.0d0) THEN
                tau = 3.0d0 / SQRT(s2)
                m(i)   = tau * alpha * delta(i)
                m(i+1) = tau * beta  * delta(i)
            END IF
        END IF
      END DO
      DEALLOCATE(h, delta, dl, d, du)

      RETURN

      END SUBROUTINE mumaterial_getstate_slopes

!-----------------------------------------------------------------------
! mumaterial_getstate_scalar_batch: Interpolates a function f at 
! many evaluation points xq using Cubic Hermite polynomials. With derivative.
!-----------------------------------------------------------------------
! param[in]:  fx. x-coordinates of function to be interpolated
! param[in]:  fy. y-values of function to be interpolated
! param[in]:  x. evaluation points
! param[out]: y. interpolated f(x)
! param[out]: dydx. derivative at x
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_getstate_scalar_batch(fx, fy, m, x, y, dydx)
        USE qsort

        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: fx(:), fy(:), m(:), x(:)
        DOUBLE PRECISION, INTENT(OUT) :: y(:)
        DOUBLE PRECISION, INTENT(OUT) :: dydx(:)

        DOUBLE PRECISION :: x_s(SIZE(x)), xt, y_s(SIZE(x)), dydx_s(SIZE(x))
        INTEGER :: i, j, k, n, nq, idx(SIZE(x))

        DOUBLE PRECISION :: t, h, hinv, t2, t3, mk, mk1
        DOUBLE PRECISION :: fxk, fxk1, fyk, fyk1
        DOUBLE PRECISION :: fx1, fxn, fy1, fyn


  
        ! First sort
        nq = SIZE(x)
        DO i = 1, nq
          idx(i) = i
        END DO
        CALL quicksort(x, idx, 1, nq)
        x_s = x(idx(:))

        ! Now loop
        n = SIZE(fx)
        fx1 = fx(1)
        fxn = fx(n)  
        fy1 = fy(1)
        fyn = fy(n)

        k = 1
        DO i = 1, nq
          xt = x_s(i)
          ! Clamp to range
          IF (xt .LE. fx1) THEN
            y_s(i) = fy1
            dydx_s(i) = 0.0d0
            CYCLE
          ELSE IF (xt .GE. fxn) THEN
            y_s(i) = fyn
            dydx_s(i) = 0.0d0
            CYCLE
          END IF
  
          ! Find interval
          DO WHILE (xt .GT. fx(k+1) .AND. k .LT. (n-1))
            k = k + 1
          END DO

          ! Spline x & f(x)
          fxk = fx(k)
          fxk1 = fx(k+1)
          fyk = fy(k)
          fyk1 = fy(k+1)
          ! Spline slopes
          mk = m(k)
          mk1 = m(k+1)
    
          ! Cubic Hermite interpolation
          h = fxk1 - fxk
          hinv = 1.0d0/h
          t = (xt - fxk)*hinv
          t2 = t*t
          t3 = t2*t
    
          y_s(i) = (2.0d0*t3 - 3.0d0*t2 + 1.0d0) * fyk   &
                  + (t3 - 2.0d0*t2 + t)          * h * mk  &
                  + (-2.0d0*t3 + 3.0d0*t2)       * fyk1 &
                  + (t3 - t2)                    * h * mk1
    
          ! Analytical derivative
          dydx_s(i) = (6.0d0*t2 - 6.0d0*t) * fyk * hinv      &
                + (3.0d0*t2 - 4.0d0*t + 1.0d0) * mk      &
                + (-6.0d0*t2 + 6.0d0*t) * fyk1 * hinv    &
                + (3.0d0*t2 - 2.0d0*t) * mk1
        END DO

        ! Scatter back
        DO i = 1, nq
          j = idx(i)
          y(j) = y_s(i)
          dydx(j) = dydx_s(i)
        END DO

        RETURN
  
        END SUBROUTINE mumaterial_getstate_scalar_batch

!-----------------------------------------------------------------------
! mumaterial_getstate_scalar_new: Interpolates a function f at xq to get 
! a value y using Cubic Hermite polynomials.
!-----------------------------------------------------------------------
! param[in]:  fx. x-coordinates of function to be interpolated
! param[in]:  fy. y-values of function to be interpolated
! param[in]:  xq. evaluation point
! param[out]: yq. interpolated f(xq)
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_getstate_scalar_new(fx, fy, m, xq, yq, dydx)

      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: fx(:), fy(:), m(:), xq
      DOUBLE PRECISION, INTENT(OUT) :: yq
      DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: dydx

      INTEGER :: n, k, i
      DOUBLE PRECISION :: t, h, t2, t3, mk, mk1
      LOGICAL :: lderiv

      lderiv = PRESENT(dydx)
      n = SIZE(fx)

      ! Clamp to range
      IF (xq .LE. fx(1)) THEN
        yq = fy(1)
        IF (lderiv) dydx = 0.0d0
        RETURN
      ELSE IF (xq .GE. fx(n)) THEN
        yq = fy(n)
        IF (lderiv) dydx = 0.0d0
        RETURN
      END IF

      ! Find interval
      DO i = 2, n
        IF (xq .LT. fx(i)) THEN
          k = i - 1
          EXIT
        END IF
      END DO

      ! Slopes
      mk = m(k)
      mk1 = m(k+1)

      ! Cubic Hermite interpolation
      h = fx(k+1) - fx(k)
      t = (xq - fx(k)) / h
      t2 = t*t
      t3 = t2*t

      yq = (2.0d0*t3 - 3.0d0*t2 + 1.0d0) * fy(k)   &
          + (t3 - 2.0d0*t2 + t)          * h * mk  &
          + (-2.0d0*t3 + 3.0d0*t2)       * fy(k+1) &
          + (t3 - t2)                    * h * mk1

      ! Analytical derivative
      IF (lderiv) THEN
        dydx = (6.0d0*t2 - 6.0d0*t) * fy(k) / h      &
              + (3.0d0*t2 - 4.0d0*t + 1.0d0) * mk      &
              + (-6.0d0*t2 + 6.0d0*t) * fy(k+1) / h    &
              + (3.0d0*t2 - 2.0d0*t) * mk1
      END IF

      RETURN

      END SUBROUTINE mumaterial_getstate_scalar_new
!-----------------------------------------------------------------------
! mumaterial_getstate_scalar: Interpolates a function f at xq to get a value y 
! using B-splines based on De Boor's algorithm
!-----------------------------------------------------------------------
! param[in]:  fx. x-coordinates of function to be interpolated
! param[in]:  fy. y-values of function to be interpolated
! param[in]:  xq. evaluation point
! param[out]: yq. interpolated f(xq)
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_getstate_scalar(fx, fy, xq, yq, dydx)

      ! Assume polynomial degree p  3
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: fx(:), fy(:), xq
      DOUBLE PRECISION, INTENT(OUT) :: yq
      DOUBLE PRECISION, INTENT(out), OPTIONAL :: dydx
      INTEGER :: n, i, k, p, r
      DOUBLE PRECISION :: alpha
      LOGICAL :: lderiv
      ! size(t) = 2+2*(p-1) = 6; size(d) = p+1 = 4
      DOUBLE PRECISION :: t(6), d(4), dder(4)
      DOUBLE PRECISION, PARAMETER :: eps = TINY(1.0d0)
      
      dder = 0.0d0
      lderiv = PRESENT(dydx)
      p = 3 ! Degree of the polynomial used
      n = SIZE(fx)

      ! Determine left index k
      ! Assume fx is non-decreasing
      IF (xq .lt. fx(1)) THEN
        yq = fy(1)
        IF (lderiv) dydx = 0.0d0
        RETURN
      ELSEIF (xq .gt. fx(n)) THEN
        yq = fy(n)
        IF (lderiv) dydx = 0.0d0
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
          alpha = (xq - t(i)) / (t(i+1+p-r) - t(i) + eps)
          d(i) = (1 - alpha) * d(i) + alpha * d(i+1)
          dder(i) = (d(i+1) - d(i)) / (t(i+1+p-r) - t(i) + eps)
        END DO
      END DO

      ! Set output variable
      yq = d(3)
      IF (lderiv) dydx = dder(3)

      RETURN
      END SUBROUTINE mumaterial_getstate_scalar

!-----------------------------------------------------------------------
! mumaterial_getstate_vector: Vector form of above. Calculates a vector
! magnitude using mumaterial_getstate_vector, then outputs a vector of
! that magnitude aligned with the input vector.
!-----------------------------------------------------------------------
! param[in]:  fx. x-coordinates of function to be interpolated
! param[in]:  fy. y-values of function to be interpolated
! param[in]:  Xq. evaluation vector
! param[out]: Yq. interpolated fy(|Xq|) aligned with Xq
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_getstate_vector(fx, fy, m, Xq, Yq, dydx)
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in)  :: fx(:), fy(:), Xq(3), m(:)
      DOUBLE PRECISION, INTENT(out) :: Yq(3)
      DOUBLE PRECISION, INTENT(out), OPTIONAL :: dydx
      DOUBLE PRECISION :: Xnorm, Ynorm

      Xnorm = NORM2(Xq)
      CALL mumaterial_getstate_scalar_new(fx, fy, m, Xnorm, Ynorm, dydx)
      IF (Xnorm .GT. small) THEN
        Yq = Ynorm * Xq/Xnorm
      ELSE
        Yq = 0.0d0
      END IF
      RETURN

      END SUBROUTINE mumaterial_getstate_vector

!-----------------------------------------------------------------------
! mumaterial_getstate_vector_batch: Vector form of above. Calculates a vector
! magnitude using mumaterial_getstate_vector, then outputs a vector of
! that magnitude aligned with the input vector.
!-----------------------------------------------------------------------
! param[in]:  fx. x-coordinates of function to be interpolated
! param[in]:  fy. y-values of function to be interpolated
! param[in]:  Xq. evaluation vector
! param[out]: Yq. interpolated fy(|Xq|) aligned with Xq
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_getstate_vector_batch(fx, fy, m, x, y, n, dydx)
      IMPLICIT NONE

      INTEGER, INTENT(in) :: n
      DOUBLE PRECISION, INTENT(in)  :: fx(:), fy(:), x(3,n), m(:)
      DOUBLE PRECISION, INTENT(out) :: y(3,n), dydx(n)
      DOUBLE PRECISION :: xnorm(n), ynorm(n), alpha(n)
      INTEGER :: i
      
      DO i = 1, n
        xnorm(i) = SQRT(x(1,i)*x(1,i)+x(2,i)*x(2,i)+x(3,i)*x(3,i))
      END DO

      CALL mumaterial_getstate_scalar_batch(fx, fy, m, xnorm, ynorm, dydx)
      alpha = 0.0d0
      WHERE (xnorm .GT. small)
        alpha = ynorm / xnorm
      END WHERE
      y(1,:) = x(1,:) * alpha
      y(2,:) = x(2,:) * alpha
      y(3,:) = x(3,:) * alpha
     
      RETURN

      END SUBROUTINE mumaterial_getstate_vector_batch

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
      SUBROUTINE mumaterial_getb_scalar(x, y, z, Bx, By, Bz)

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

!-----------------------------------------------------------------------
! mumaterial_getbmag: Calculates magnetic field from magnetizations only
!-----------------------------------------------------------------------
! param[in]: x. x-coordinate of point where to get the B-field
! param[in]: y. y-coordinate of point where to get the B-field
! param[in]: z. z-coordinate of point where to get the B-field
! param[out]: Bx. x-component of B-field at this point [T]
! param[out]: By. y-component of B-field at this point [T]
! param[out]: Bz. z-component of B-field at this point [T]
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_getbmag_scalar(x, y, z, Bx, By, Bz)

      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x, y, z
      DOUBLE PRECISION, INTENT(out) :: Bx, By, Bz
      DOUBLE PRECISION :: H(3), N(3,3), r(3), hx, hy, hz, d

      INTEGER :: i, j, k, n_nbrs, i_tile, n_dip

      !-----------------------------
      ! Set everything up
      !----------------------------- 
      r(1) = x; r(2) = y; r(3) = z

      ! Get clusters
      DO i = 1, world_size
        iscluster_proc(i) = NORM2(r_cluster(:,i)-r)>(cutoff*d_cluster(i) + pF*max_tet_rad(i))
      END DO

      ! Get neighbors
      nbrs_out = 0; mask_out = .FALSE.
      DO i = 1, world_size
        IF (iscluster_proc(i)) CYCLE
        DO j = 1, dom_sizes(i)
          k = dom_cluster(j,i)
          d = (tet_cen(1,k)-r(1))*(tet_cen(1,k)-r(1)) + &
              (tet_cen(2,k)-r(2))*(tet_cen(2,k)-r(2)) + &
              (tet_cen(3,k)-r(3))*(tet_cen(3,k)-r(3))
          mask_out(k) = d < (pF*tet_rad(k))**2
        END DO
      END DO
      n_nbrs = COUNT(mask_out)
      IF (n_nbrs>0) nbrs_out(1:n_nbrs) = PACK([(i, i=1, ntet)], mask_out)

      DO i = 1, world_size
        IF (.NOT.iscluster_proc(i)) CYCLE ! Already not a cluster, skip
        IF (ANY(mask_out(dom_cluster(1:dom_sizes(i),i)))) iscluster_proc(i) = .FALSE.
      END DO

      ! Get all non-clusters
      mask_out = .TRUE. 
      ! Remove clusters
      DO i = 1, world_size
        IF (.NOT.iscluster_proc(i)) CYCLE
        mask_out(dom_cluster(1:dom_sizes(i),i)) = .FALSE.
      END DO
      n_dip = COUNT(mask_out)
      dom_dip(1:n_dip) = PACK([(i, i=1, ntet)], mask_out)
      ! Fetch
      DO i = 1, n_dip
        M_dip(:,i) = M_global(:,dom_dip(i))
        r_dip(:,i) = tet_cen(:,dom_dip(i))
        V_dip(i)   = tet_vol(dom_dip(i))
      END DO

      !-----------------------------
      ! Finally, calculate fields
      !----------------------------- 
      H = 0.d0; hx = 0.0d0; hy = 0.0d0; hz = 0.0d0
      ! Neighbor field
      DO i = 1, n_nbrs
        i_tile = nbrs_out(i)
        N = GET_DEMAG(tet_P(:,:,:,i_tile), tet_D(:,:,i_tile), tet_v(:,:,:,i_tile), [x, y, z])
        H = H + MATMUL(N, M_global(:,i_tile))
        ! Remove dipole contribution
        CALL mumaterial_dipole_field_calc_single(x, y, z, tet_cen(:,i_tile), tet_vol(i_tile), -M_global(:,i_tile), hx, hy, hz)
      END DO
      ! Cluster field
      CALL mumaterial_cluster_field_calc(x, y, z, hx, hy, hz)
      ! Dipole field (includes extra copy of neighbor contributions)
      CALL mumaterial_dipole_field_calc(x, y, z, n_dip, r_dip, V_dip, M_dip, hx, hy, hz)

      Bx = (H(1)+hx) * MU0; By = (H(2)+hy) * MU0; Bz = (H(3)+hz) * MU0

      RETURN
      END SUBROUTINE mumaterial_getbmag_scalar

!-----------------------------------------------------------------------
! mumaterial_getb_vector: Calculates total magnetic field at multiple points in space
!-----------------------------------------------------------------------
! param[in]: x. x-coordinates of points at which to determine the magnetic field
! param[in]: y. y-coordinates of points at which to determine the magnetic field
! param[in]: z. z-coordinates of points at which to determine the magnetic field
! param[out]: B.  B-field at required points [T]
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_getb_vector(x, y, z, B)

      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x(:), y(:), z(:)
      DOUBLE PRECISION, INTENT(out), ALLOCATABLE :: B(:,:)
      INTEGER :: mystart, myend
      INTEGER :: i 
      INTEGER :: npoints

      npoints = size(x)
      mystart = 1; myend = npoints

#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_CALC_MYRANGE(comm_world, 1, npoints, mystart, myend)
#endif

      ALLOCATE(B(3,npoints))
      B = 0
      DO i = mystart, myend
        CALL mumaterial_getb_scalar(x(i), y(i), z(i), B(1,i), B(2,i), B(3,i))
      END DO
    
#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_ALLREDUCE( MPI_IN_PLACE,B,3*npoints,MPI_DOUBLE_PRECISION,MPI_SUM,comm_world,ierr_mpi)
#endif
      
      RETURN
      END SUBROUTINE mumaterial_getb_vector

!-----------------------------------------------------------------------
! mumaterial_magfile_read: Reads magnetization .dat file and synchronizes
!-----------------------------------------------------------------------   
      SUBROUTINE mumaterial_magfile_read(filename)

      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER :: i, istat, iunit

      IF (VERIFY(filename, CHAR(0)//' ') == 0) RETURN

      ! Read file
      IF (lismaster) THEN
        WRITE(6,'(A)')           ' -------- MUMAT magfile --------'
        WRITE(6,'(3X,A,A)')     ' File         : ',filename
        ! open file, return if fails
        iunit = 327; istat = 0
        CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
        IF (istat/= 0) THEN
          WRITE(6,*) "ISSUE READING MAG; STOPPING"
          CALL mumaterial_abort()
        END IF
        DO i = 1, ntet
          READ(iunit, *) M_global(1,i),M_global(2,i),M_global(3,i)
        END DO
        CLOSE(iunit)
      END IF    

      ! Sync
      CALL mumaterial_syncmag(SYNC_READ)
      RETURN

      END SUBROUTINE mumaterial_magfile_read

!-----------------------------------------------------------------------
! mumaterial_magfile_write: Outputs magnetization to .dat file
!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_magfile_write(str)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(in) :: str
      CHARACTER(LEN=256) :: filename
      INTEGER :: i

      IF (lismaster) THEN
        filename = './mumat_mag_'//TRIM(str)//'.dat'
        WRITE(6,"(A)") "  MUMAT: Writing magnetization to" // filename
        OPEN(13, file=filename)
        DO i = 1, ntet
          WRITE(13, "(ES15.7,ES15.7,ES15.7)") M_global(1,i),M_global(2,i),M_global(3,i)
        END DO
        CLOSE(13)
      END IF

      END SUBROUTINE mumaterial_magfile_write

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

      CALL mumaterial_alloc_out()
      CALL mumaterial_getb_vector(x, y, z, B)
      CALL mumaterial_dealloc_out()
 
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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!       PRIVATE HELPER FUNCTIONS      !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------------------------------------------------------------------
! EYE: Identity matrix of size n.
!-----------------------------------------------------------------------
! param[in]: n. size
!-----------------------------------------------------------------------      
      FUNCTION EYE(n) result(aye)

      IMPLICIT NONE
      INTEGER, INTENT(in) :: n
      INTEGER :: aye(n,n), i

      aye = 0
      DO i = 1, n
        aye(i,i) = 1
      END DO

      END FUNCTION EYE

!-----------------------------------------------------------------------
! GET_DEMAG_HELPERS: Helper function to determine the demagnetization tensor
!-----------------------------------------------------------------------
! param[in]: v1-4: Vertices of the tetrahedron
! param[out]: P: Rotation matrix of face
! param[out]: D: Face base 
! param[out]: v_rot: Vertices of each face
!-----------------------------------------------------------------------
      SUBROUTINE GET_DEMAG_HELPERS(v1, v2, v3, v4, P, D, v_rot)

      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in), DIMENSION(3) :: v1, v2, v3, v4
      DOUBLE PRECISION, INTENT(out) :: P(3,3,4), D(3,4), v_rot(3,3,4)
      DOUBLE PRECISION :: v(3,4), v_temp(3), cosalpha(3), d12, d13, d23
      INTEGER :: i, j


      ! Rotate vertices
      DO i = 1, 4
        v(:,i)            = v1
        v(:,MOD(i  ,4)+1) = v2
        v(:,MOD(i+1,4)+1) = v3
        v(:,MOD(i+2,4)+1) = v4

        ! todo: ensure vertices are not colinear and v4 is not in plane of v1-3?

        ! Ensure largest angle is for v2
        d12 = NORM2(v(:,1)-v(:,2))
        d13 = NORM2(v(:,1)-v(:,3))
        d23 = NORM2(v(:,2)-v(:,3))
        cosalpha(1) = MAX(1.0d0, MIN(1.0d0, DOT_PRODUCT(v(:,1)-v(:,2),v(:,1)-v(:,3)) / (d12*d13)))
        cosalpha(2) = MAX(1.0d0, MIN(1.0d0, DOT_PRODUCT(v(:,2)-v(:,1),v(:,2)-v(:,3)) / (d12*d23)))
        cosalpha(3) = MAX(1.0d0, MIN(1.0d0, DOT_PRODUCT(v(:,3)-v(:,2),v(:,3)-v(:,1)) / (d13*d23)))

        IF (cosalpha(1) < cosalpha(2) .and. cosalpha(1) < cosalpha(3)) THEN ! v1 and v2 should be interchanged
          v_temp = v(:,2)
          v(:,2) = v(:,1)
          v(:,1) = v_temp
        ELSE IF (cosalpha(3) < cosalpha(1) .and. cosalpha(3) < cosalpha(2)) THEN ! v2 and v3 should be interchanged
          v_temp = v(:,2)
          v(:,2) = v(:,3)
          v(:,3) = v_temp
        END IF

        ! normal vector of triangle is pointing towards v4, so v1 and v3 need to be interchanged
        IF (DOT_PRODUCT(CROSS_PRODUCT(v(:,1) - v(:,3), v(:,2) - v(:,3)), v(:,4) - v(:,2)) .gt. 0) THEN 
          v_temp = v(:,1)
          v(:,1) = v(:,3)
          v(:,3) = v_temp
        END IF
        v_rot(:,1,i) = v(:,1)
        v_rot(:,2,i) = v(:,2)
        v_rot(:,3,i) = v(:,3)
      END DO

      P = 0.0d0
      D = 0.0d0

      DO i = 1, 4
        ! Rotation matrix
        v(:,1) = v_rot(:,1,i)
        v(:,2) = v_rot(:,2,i)
        v(:,3) = v_rot(:,3,i)

        P(:,1,i) = v(:,1) - v(:,3)
        P(:,1,i) = P(:,1,i) / NORM2(P(:,1,i))

        P(:,3,i) = CROSS_PRODUCT(P(:,1,i), v(:,2)-v(:,3))
        P(:,3,i) = P(:,3,i) / NORM2(P(:,3,i))

        P(:,2,i) = CROSS_PRODUCT(P(:,3,i), P(:,1,i))
        P(:,2,i) = P(:,2,i) / NORM2(P(:,2,i))

        ! Position of triangle base
        d13 = NORM2(v(:,1)-v(:,3))
        d23 = NORM2(v(:,2)-v(:,3))
        D(:,i) = DOT_PRODUCT(v(:,3)-v(:,2),v(:,3)-v(:,1)) / (d23 * d13) * d23 * P(:,1,i) + v(:,3)

      END DO
      RETURN

      END SUBROUTINE GET_DEMAG_HELPERS

!-----------------------------------------------------------------------
! GET_DEMAG: Helper function to determine the demagnetization tensor
!-----------------------------------------------------------------------
! param[in]: pos: Reference position for which to determine the demagnetization tensor (3)
! param[in]: N. Resulting demagnetization tensor (3,3)
!-----------------------------------------------------------------------
      FUNCTION GET_DEMAG(P, D, v_rot, pos) result(N)

      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in), DIMENSION(3) :: P(3,3,4), D(3,4), v_rot(3,3,4), pos(3)
      DOUBLE PRECISION :: N(3,3)

      DOUBLE PRECISION :: N_loc(3,3), v(3,4), v_temp(3), cosalpha(3), Pinv(3,3), r(3), d12, d13, d23
      INTEGER :: i, j

      N = 0.d0

      DO i = 1, 4
        ! Rotated vertices
        v(:,1) = v_rot(:,1,i)
        v(:,2) = v_rot(:,2,i)
        v(:,3) = v_rot(:,3,i)

        Pinv = TRANSPOSE(P(:,:,i))

        ! Transform evaluation position and vertices to local coordinate frame
        r = MATMUL(Pinv, (pos - D(:,i)))

        DO j = 1, 3
          v(:,j) = MATMUL(Pinv, (v(:,j) - D(:,i)))
          IF (ABS(r(j)) .lt. 1.0D-6) THEN ! make sure position is not too close to x, y or z = 0
            r(j) = SIGN(1.0D-6, r(j))
          END IF
          IF (ABS(v(1,j)) .lt. 1.0D-6) THEN 
            v(1,j) = SIGN(1.0D-6, v(1,j))
          END IF                  
        END DO

        N_loc = 0.d0
        N_loc(1,3) = GET_Nxz(r, v(1,1), v(2,2)) - GET_Nxz(r, v(1,3), v(2,2))
        N_loc(2,3) = GET_Nyz(r, v(1,1), v(2,2)) - GET_Nyz(r, v(1,3), v(2,2))
        N_loc(3,3) = GET_Nzz(r, v(1,1), v(2,2)) - GET_Nzz(r, v(1,3), v(2,2))
        IF ((ISNAN(N_loc(1,3)).or.ISNAN(N_loc(2,3))).or.ISNAN(N_loc(3,3))) THEN 
              WRITE(6,*) "FOUND A NAN IN N_LOC"
              WRITE(6,*) "POS=",pos(1),pos(2),pos(3)
              WRITE(6,*) "R=",r(1),r(2),r(3)
              WRITE(6,*) "l_1=",v(1,1), "l_2=", v(1,3)
              WRITE(6,*) "h=",v(2,2)
              WRITE(6,*)
        END IF
        N = N + MATMUL(MATMUL(P(:,:,i), N_loc), Pinv)
      END DO

      RETURN
      END FUNCTION GET_DEMAG

!-----------------------------------------------------------------------
! GET_Nxz: Helper function to determine the x-component of the demagnetization tensor
!-----------------------------------------------------------------------
! param[in]: rp: Reference position for which to determine the demagnetization tensor (3)
! param[in]: lp: Bottom side of the triangle 
! param[in]: hp. Top side of the triangle
!-----------------------------------------------------------------------
      FUNCTION GET_Nxz(rp, lp, hp) result(Nxz)

      IMPLICIT NONE
      DOUBLE PRECISION :: Nxz
      DOUBLE PRECISION, INTENT(IN) :: rp(3), lp, hp     

      Nxz = -INV4PI * (F(rp,hp,lp,hp) - F(rp,0.d0,lp,hp) - (G(rp,hp) - G(rp,0.d0)))

      CONTAINS

        FUNCTION F(r, yp, l, h) result(val)
        IMPLICIT NONE
        DOUBLE PRECISION :: val
        DOUBLE PRECISION, INTENT(IN) :: r(3), yp, l, h
              
        val = h / sqrt(h*h + l*l) * ATANH((l*l - l*r(1) + h*r(2) - h*yp*(1 + l*l/h/h)) / &
                sqrt((h*h + l*l) * (r(1)*r(1) - 2*r(1)*l + l*l + r(2)*r(2) - 2*(l*l - l*r(1) + h*r(2))*yp/h + &
                yp*yp*(1 + l*l/h/h) + r(3)*r(3))))
  
        END FUNCTION F
  
        FUNCTION G(r, yp) result(val)
        IMPLICIT NONE
        DOUBLE PRECISION :: val
        DOUBLE PRECISION, INTENT(IN) :: r(3), yp
              
        val = ATANH((r(2) - yp) / sqrt(r(1)*r(1) + r(2)*r(2) - 2*r(2)*yp + yp*yp + r(3)*r(3)))
              
        END FUNCTION G
      END FUNCTION GET_Nxz

!-----------------------------------------------------------------------
! GET_Nyz: Helper function to determine the y-component of the demagnetization tensor
!-----------------------------------------------------------------------
! param[in]: rp: Reference position for which to determine the demagnetization tensor (3)
! param[in]: lp: Bottom side of the triangle 
! param[in]: hp. Top side of the triangle
!-----------------------------------------------------------------------  
      FUNCTION GET_Nyz(rp, lp, hp) result(Nyz)

      IMPLICIT NONE
      DOUBLE PRECISION :: Nyz
      DOUBLE PRECISION, INTENT(IN) :: rp(3), lp, hp

      Nyz = -INV4PI * (K(rp,lp,lp,hp) - K(rp,0.d0,lp,hp) - (Lfunc(rp,lp) - Lfunc(rp,0.d0)))

      CONTAINS

        FUNCTION K(r, xp, l, h) result(val)
        IMPLICIT NONE
        DOUBLE PRECISION :: val
        DOUBLE PRECISION, INTENT(IN) :: r(3), xp, l, h
  
        val = l / sqrt(h*h + l*l) * ATANH((h*h + l*r(1) - h*r(2) - l*xp*(1 + h*h/l/l)) / &
              sqrt((h*h + l*l) * (r(2)*r(2) - 2*r(2)*h + h*h + r(1)*r(1) - 2*(h*h + l*r(1) - h*r(2))*xp/l + &
              xp*xp*(1 + h*h/l/l) + r(3)*r(3))))
  
        END FUNCTION K
  
        FUNCTION Lfunc(r, xp) result(val)
        IMPLICIT NONE
        DOUBLE PRECISION :: val
        DOUBLE PRECISION, INTENT(IN) :: r(3), xp
        
        val = ATANH((r(1) - xp) / sqrt(r(1)*r(1) - 2*r(1)*xp + xp*xp + r(2)*r(2) + r(3)*r(3)))
              
        END FUNCTION Lfunc

      END FUNCTION GET_Nyz

!-----------------------------------------------------------------------
! GET_Nzz: Helper function to determine the z-component of the demagnetization tensor
!-----------------------------------------------------------------------
! param[in]: rp: Reference position for which to determine the demagnetization tensor (3)
! param[in]: lp: Bottom side of the triangle 
! param[in]: hp. Top side of the triangle
!-----------------------------------------------------------------------
      FUNCTION GET_Nzz(rp, lp, hp) result(Nzz)

      IMPLICIT NONE
      DOUBLE PRECISION :: Nzz
      DOUBLE PRECiSION, INTENT(IN) :: rp(3), lp, hp
            
      Nzz = -INV4PI * (P(rp,lp,lp,hp) - P(rp,0.d0,lp,hp) - (Q(rp,lp) - Q(rp,0.d0)))
      
      CONTAINS

        FUNCTION P(r, xp, l, h) result(val)
        IMPLICIT NONE
        DOUBLE PRECISION :: val
        DOUBLE PRECISION, INTENT(IN) :: r(3), xp, l, h
  
        val = ATAN((r(1)*(h - r(2)) - xp*(h*(1 - r(1)/l) - r(2)) - h*(r(1)*r(1) + r(3)*r(3))/l) / &
              (r(3)*sqrt(r(2)*r(2) - 2*r(2)*h + h*h + r(1)*r(1) + xp*xp*(1 + h*h/l/l) - &
              2*xp*(h*h + l*r(1) - h*r(2))/l + r(3)*r(3))))
  
        END FUNCTION P
  
        FUNCTION Q(r, xp) result(val)
        IMPLICIT NONE
        DOUBLE PRECISION :: val
        DOUBLE PRECISION, INTENT(IN) :: r(3), xp
  
        val = -ATAN((r(1) - xp)*r(2) / (r(3)*sqrt((r(1)*r(1) - 2*r(1)*xp + xp*xp + r(2)*r(2) + r(3)*r(3)))))
        
        END FUNCTION Q
  
      END FUNCTION GET_Nzz

!-----------------------------------------------------------------------
! CROSS_PRODUCT: calculates the cross product c = a * b.
!-----------------------------------------------------------------------
! param[in]: a. 1x3 vector
! param[in]: b: 1x3 vector
!-----------------------------------------------------------------------
      FUNCTION CROSS_PRODUCT(a, b) result(c)

      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN), DIMENSION(3) :: a, b
      DOUBLE PRECISION, DIMENSION(3) :: c

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)

      RETURN
      END FUNCTION CROSS_PRODUCT
!-----------------------------------------------------------------------
! OUTER_PRODUCT: calculates the outer product c = a*b'.
!-----------------------------------------------------------------------
! param[in]: a. 1x3 vector
! param[in]: b: 1x3 vector
!-----------------------------------------------------------------------
      FUNCTION OUTER_PRODUCT(a, b) result(c)

      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN), DIMENSION(3) :: a, b
      DOUBLE PRECISION, DIMENSION(3,3) :: c
      INTEGER :: i, j
      c = 0.0d0
      DO i = 1, 3
        DO j = 1, 3
          c(i,j) = a(i)*b(j)
        END DO
      END DO

      RETURN
      END FUNCTION OUTER_PRODUCT

!-----------------------------------------------------------------------
! TET_VOLUME: Calculates the volume of an element.
!-----------------------------------------------------------------------
      FUNCTION TET_VOLUME(v1,v2,v3,v4) result(V)

      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3), INTENT(in) :: v1, v2, v3, v4
      DOUBLE PRECISION :: V

      V = ABS(DOT_PRODUCT(v1-v4,CROSS_PRODUCT(v2-v4,v3-v4)))/6.0
      RETURN

      END FUNCTION TET_VOLUME

!-----------------------------------------------------------------------
! SOLVE_3x3: Helper function for solving a 3*3 matrix problem with LAPACK
!-----------------------------------------------------------------------
! param[in]: A. 3x3 matrix
! param[in]: B. 1x3 vector (RHS)
! param[out]: C. 1x3 vector satisfying C = A\B
!-----------------------------------------------------------------------
      FUNCTION SOLVE_3x3(A, B, LOC) RESULT(C)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: A(3,3), B(3)
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: LOC
      DOUBLE PRECISION :: A_work(3,3), B_work(3), C(3)
      ! LAPACK
      INTEGER :: INFO, IPIV(3)
      INTEGER, PARAMETER :: N = 3, NRHS = 1, LDA = 3, LDB = 3

      A_work = A
      B_work = B
      CALL DGESV(N, NRHS, A_work, LDA, IPIV, B_work, LDB, INFO)

      IF (INFO.NE.0) THEN
        IF (PRESENT(LOC)) THEN 
          WRITE(6,'(A,I0,A,A,A,I0)') "RANK ", world_rank, ": DGESV FAILED (",LOC ,"), INFO=", INFO
        ELSE
          WRITE(6,'(A,I0,A,I0)')  "RANK ", world_rank, "ERROR: DGESV FAILED, INFO=", INFO
        END IF
        WRITE(6,'(A,I0, 3ES12.4)') "RANK ", world_rank, A(1,1), A(1,2), A(1,3)
        WRITE(6,'(A,I0, 3ES12.4)') "RANK ", world_rank, A(2,1), A(2,2), A(2,3)
        WRITE(6,'(A,I0, 3ES12.4)') "RANK ", world_rank, A(3,1), A(3,2), A(3,3)
        WRITE(6,'(A,I0, 3ES12.4)') "RANK ", world_rank, B(1), B(2), B(3)

        CALL mumaterial_abort()
      END IF
      C = B_work
      RETURN

      END FUNCTION SOLVE_3x3

!-----------------------------------------------------------------------
! INVERT_3x3: Helper function for inverting a 3x3 matrix.
!-----------------------------------------------------------------------
! param[in]: B: target matrix 
!-----------------------------------------------------------------------
      FUNCTION INVERT_3x3(B) result(INV)

        IMPLICIT NONE
  
        DOUBLE PRECISION, INTENT(in) :: B(3,3)
        DOUBLE PRECISION :: det, INV(3,3)
  
        DET = B(1,1)*(B(2,2)*B(3,3)-B(2,3)*B(3,2)) + &
              B(1,2)*(B(2,3)*B(3,1)-B(2,1)*B(3,3)) + &
              B(1,3)*(B(2,1)*B(3,2)-B(2,2)*B(3,1))
  
        INV(1,1) = (B(2,2)*B(3,3)-B(2,3)*B(3,2))/DET
        INV(2,1) =-(B(2,3)*B(3,1)-B(2,1)*B(3,3))/DET
        INV(3,1) = (B(2,1)*B(3,2)-B(2,2)*B(3,1))/DET
  
        INV(1,2) =-(B(1,3)*B(3,2)-B(1,2)*B(3,3))/DET
        INV(2,2) = (B(1,1)*B(3,3)-B(1,3)*B(3,1))/DET
        INV(3,2) =-(B(1,2)*B(3,1)-B(1,1)*B(3,2))/DET
  
        INV(1,3) = (B(1,2)*B(2,3)-B(1,3)*B(2,2))/DET
        INV(2,3) =-(B(1,3)*B(2,1)-B(1,1)*B(2,3))/DET
        INV(3,3) = (B(1,1)*B(2,2)-B(1,2)*B(2,1))/DET
  
        END FUNCTION INVERT_3x3

!-----------------------------------------------------------------------
! COND_3x3: Helper function for getting matrix condition number.
!-----------------------------------------------------------------------
! param[in]: MAT: target matrix 
!-----------------------------------------------------------------------
        FUNCTION COND_3x3(MAT) result(cond)

        IMPLICIT NONE
  
        DOUBLE PRECISION, INTENT(in) :: MAT(3,3)
        DOUBLE PRECISION :: cond, INVMAT(3,3), row_sum, val1,val2
        INTEGER :: i

        val1 = 0.0d0
        DO i = 1, 3
          row_sum = ABS(MAT(i,1))+ABS(MAT(i,2))+ABS(MAT(i,3))
          IF (row_sum.GT.val1) val1= row_sum
        END DO

        INVMAT = INVERT_3x3(MAT)
        val2 = 0.0d0
        DO i = 1, 3
          row_sum = ABS(INVMAT(i,1))+ABS(INVMAT(i,2))+ABS(INVMAT(i,3))
          IF (row_sum.GT.val2) val2 = row_sum
        END DO
        
        cond = val1*val2
  
        END FUNCTION COND_3x3
!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE mumaterial_mod

