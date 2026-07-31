!------------------------------------------------------------------------------
!     Module:        mumaterial_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov), Björn Hamstra,
!                    Lucas van Ham (lucas.van.ham@ipp.mpg.de)
!     Date:          October 2023
!                    Jan-June 2024 [LvH]: Modifications for scaled-up problems,
!                    neighbour formulation, dipole approximation
!                    April 2026 [LvH]: Barnes-Hut
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
      USE MPI_SHARMEM, ONLY: mpialloc, mpidealloc
#endif
      IMPLICIT NONE

      TYPE stateFunctionType
            DOUBLE PRECISION, PRIVATE, ALLOCATABLE :: H(:) !! H-values
            DOUBLE PRECISION, PRIVATE, ALLOCATABLE :: M(:) !! M-values 
            DOUBLE PRECISION, PRIVATE, ALLOCATABLE :: b(:) !! Spline slopes
      END TYPE stateFunctionType

      PROCEDURE(externalFieldFunc), POINTER :: ext_Bfld !! Background field function
      ABSTRACT INTERFACE
        SUBROUTINE externalFieldFunc(x,y,z,Bx,By,Bz)
          DOUBLE PRECISION, INTENT(in)  :: x,y,z
          DOUBLE PRECISION, INTENT(out) :: Bx,By,Bz
        END SUBROUTINE externalFieldFunc
      END INTERFACE

      ! Mesh variables
      INTEGER, PRIVATE  :: ntet    !! Number of elements (tetrahedrons)
      INTEGER, PRIVATE  :: nvertex !! Number of vertices
      INTEGER, POINTER, PRIVATE :: tet(:,:) !! (4,ntet) Vertex indices of an element
      DOUBLE PRECISION, POINTER, PRIVATE :: vertex(:,:) !! (3,nvertex) array of vertices
      DOUBLE PRECISION, POINTER, PRIVATE :: tet_cen(:,:)!! (3,ntet) Geometric center of an element
      DOUBLE PRECISION, POINTER, PRIVATE :: tet_vol(:)  !! (ntet) Geometric volume of an element
      DOUBLE PRECISION, PRIVATE :: tet_vol_tot !! Total geometric volume of all elements
      CHARACTER(LEN=256), PRIVATE :: mesh_name !! Mesh name read from file
      CHARACTER(LEN=256), PRIVATE :: mesh_date !! Mesh date read from file
      DOUBLE PRECISION, PRIVATE ::  padFactor = 10  !! Affects number of neighbours for each tetrahedron

      ! Magnetics variables
      INTEGER, PRIVATE :: nstate !! Number of state functions
      TYPE(stateFunctionType), PRIVATE, ALLOCATABLE :: stateFunction(:) !! (nstate) array of state functions
      INTEGER, POINTER, PRIVATE :: state_dex(:) !! (ntet) index into nstate of an element
      INTEGER, POINTER, PRIVATE :: state_type(:) !! (nstate) array of state types (1-3)
      DOUBLE PRECISION, POINTER, PRIVATE :: constant_mu(:) !! (nstate) Overall relative permeability (stype=3) or along easy axis (stype=1)
      DOUBLE PRECISION, POINTER, PRIVATE :: constant_mu_o(:) !! (nstate) relative permeability along off axis (stype=1)
      DOUBLE PRECISION, POINTER, PRIVATE ::  M(:,:) !! (3,ntet) Magnetization of element
      DOUBLE PRECISION, POINTER, PRIVATE ::  Mrem(:,:) !! (3,nstate) Remanent magnetization of istate (stype=1)
      INTEGER, PRIVATE, PARAMETER :: STATE_HARD   = 1 !! Hard magnet with remanent magnetization
      INTEGER, PRIVATE, PARAMETER :: STATE_SOFT   = 2 !! Soft medium (nonlinear state function)
      INTEGER, PRIVATE, PARAMETER :: STATE_LINEAR = 3 !! Linear medium (constant permeability)

      ! Iterations
      DOUBLE PRECISION, PRIVATE :: eps_max = 1.0d-6 !! Threshold error for convergence
      DOUBLE PRECISION, PRIVATE :: lambdaStart = 0.5d0 !! Initial value of lambda for iterations
      DOUBLE PRECISION, PRIVATE :: lambdaFactor = 0.75d0!! Multiplication factor for lambda
      DOUBLE PRECISION, PRIVATE :: lambdaMax = 0.50d0 !! Maximum value of lambda
      DOUBLE PRECISION, PRIVATE :: lambdaMin = 0.01d0 !! Minimum value of lambda
      DOUBLE PRECISION, PRIVATE :: Pconv_min = 95.0d0 !! Minimum volume convergence percentage (%)
      INTEGER, PRIVATE :: lambdaThresh = 10 !! Multiply lambda if error grows this number of times
      INTEGER, PRIVATE :: maxiter = 10000  !!  Max number of iterations

      ! Persistent (global) tree info variables
      INTEGER, PRIVATE :: nnode !! Number of global nodes (internal or leaves)
      INTEGER, PRIVATE :: nleaf !! Number of global leaves
      INTEGER, PRIVATE :: tree_depth !! Depth of tree (root = 0)
      INTEGER, PRIVATE :: max_leaf_size
      INTEGER, POINTER, PRIVATE :: node_level(:) !! (nnode) Level of node (root=0)
      INTEGER, POINTER, PRIVATE :: node_child(:,:) !! (8,nnode) Node indices of child octants of a node
      INTEGER, POINTER, PRIVATE :: node_parent(:) !! (nnode) Parent node index of a node
      INTEGER, POINTER, PRIVATE :: leaf_to_node(:) !! (nleaf) Map from leaf index to node index
      INTEGER, POINTER, PRIVATE :: leaf_size(:) !! (nleaf) Number of elements in a leaf
      INTEGER, POINTER, PRIVATE :: leaf_tile(:,:) !! (leaf_max_size, nnode) Global indices of elements in a leaf (0 for a node)
      DOUBLE PRECISION, POINTER, PRIVATE ::  node_bounds(:,:) !! (6,nnode) [xmin xmax ymin ymax zmin zmax] of a node
      DOUBLE PRECISION, POINTER, PRIVATE ::  node_cen(:,:) !! (3,nnode) geometric center of a node
      INTEGER, PRIVATE :: leaf_max_size = 4  !! Maximum allowed number of elements in a leaf
      INTEGER, PRIVATE :: tree_max_depth = 20 !! Maximum allowed depth (takes precendence over max leaf size)
      DOUBLE PRECISION, PRIVATE :: tree_theta  = 0.33d0      !! Aggregate acceptance criterion: node_size/distance_to_node < theta
      DOUBLE PRECISION, PRIVATE :: tree_theta2 = 0.33d0**2    !! theta squared

      ! Persistent (local) tree info variables
      ! (local: subset managed by an MPI rank)
      INTEGER, PRIVATE :: ntile_loc !! Number of local elements
      INTEGER, PRIVATE :: nnode_loc !! Number of local nodes (internal or leaves)
      INTEGER, PRIVATE :: nleaf_loc !! Number of local leaves
      INTEGER, PRIVATE :: max_leaf_size_seen !! Largest number of elements seen in a leaf locally
      INTEGER, ALLOCATABLE, PRIVATE :: node_level_loc(:) !! (nnode_loc) Level of node (root=0)
      INTEGER, ALLOCATABLE, PRIVATE :: node_child_loc(:,:)! !! (8,nnode_loc) Global node indices of child octants of a node
      INTEGER, ALLOCATABLE, PRIVATE :: node_local_to_global(:) !! (nnode_loc) Map from local node index to global node index
      INTEGER, ALLOCATABLE, PRIVATE :: leaf_local_to_node_global(:) !! (nleaf_loc) Map from local leaf index to global node index
      INTEGER, ALLOCATABLE, PRIVATE :: leaf_size_loc(:) !! (nleaf_loc) Number of elements in a leaf
      INTEGER, ALLOCATABLE, PRIVATE :: leaf_tile_loc(:,:) !! (max_leaf_size_seen, nleaf_loc) Global indices of elements in a leaf
      INTEGER, ALLOCATABLE, PRIVATE :: leaf_offset_loc(:) !! (nleaf_loc) Local leaf offset; itile_loc(iloc) = leaf_offset_loc(ileaf) + leaf_tile_loc(iloc, ileaf)


      ! Temporary tree building variables
      ! Because nnode is unknown apriori, build variables are initialized with size 10*ntet
      INTEGER, PARAMETER, PRIVATE :: NODE_NOCHILD = -1 !! This entry does not correspond to a child/node
      INTEGER, PARAMETER, PRIVATE :: LEAF_NOHEAD  = -2 !! This leaf has no header element
      INTEGER, PARAMETER, PRIVATE :: TILE_NONEXT  = -3 !! This tile has no next element
      INTEGER, PARAMETER, PRIVATE :: NODE_NOTLEAF = -4 !! This node is not a leaf (internal)
      INTEGER, PARAMETER, PRIVATE :: LEAF_NOTTILE = -5 !! This entry does not correspond to a tile
      INTEGER, PARAMETER, PRIVATE :: NODE_BADPARENT= -6 !! This node has no parent!
      INTEGER, PARAMETER, PRIVATE :: ROOT = 1 !! inode of the root of the tree (1)

      ! Node aggregation variables
      DOUBLE PRECISION, POINTER, PRIVATE :: node_m(:,:) !! (3,nnode) Magnetic dipole moment of a node
      DOUBLE PRECISION, POINTER, PRIVATE :: node_p(:,:)   !! (3,nnode) First weighted position moment
      DOUBLE PRECISION, POINTER, PRIVATE :: node_w(:)     !! (nnode) Scalar weight accumulator
      DOUBLE PRECISION, POINTER, PRIVATE :: node_com(:,:) !! (3,nnode) Center of dipole strength
      DOUBLE PRECISION, POINTER, PRIVATE :: node_q(:,:)!! (6,nnode) quadrupole tensor
      INTEGER, PRIVATE :: win_node_m, win_node_p, win_node_com, win_node_w, win_node_q

      ! Magnetization calculation variables
      DOUBLE PRECISION, ALLOCATABLE, PRIVATE :: Happ_loc(:,:,:) !! (3,leaf_max_size,nleaf_loc) Static background field at element
      DOUBLE PRECISION, ALLOCATABLE, PRIVATE :: Nself_loc(:,:,:) !! (3,3,ntile_loc) Self-demagnetization tensor
      DOUBLE PRECISION, ALLOCATABLE, PRIVATE :: Nloc(:,:,:)   !! (exact interactions,3,3) Demagnetization tensor from other elements
      DOUBLE PRECISION, POINTER, PRIVATE :: tet_P(:,:,:,:) !! (3,3,4,ntet) Rotation matrix of all element faces
      DOUBLE PRECISION, POINTER, PRIVATE :: tet_D(:,:,:)   !! (3,4,ntet)   Base vectors of all element faces
      DOUBLE PRECISION, POINTER, PRIVATE :: tet_v(:,:,:,:) !! (3,3,4,ntet) Permutated vertices in local face coordinates
      INTEGER, ALLOCATABLE, PRIVATE :: nN_per_leaf(:)  !! 


      ! Interaction lists and classifications
      INTEGER, ALLOCATABLE, PRIVATE :: interact_list(:) !! Stores the indices of interacting nodes
      INTEGER, ALLOCATABLE, PRIVATE :: interact_type(:) !! 1 = node, 2 = leaf
      INTEGER, ALLOCATABLE, PRIVATE :: interact_ptr(:) !! (nleaf_loc+1) Pointer for interact_list and interact_type
      INTEGER, PARAMETER, PRIVATE :: INTERACT_LEAF = 1 !! Treat this interaction as a leaf.
      INTEGER, PARAMETER, PRIVATE :: INTERACT_NODE = 0 !! Treat this interaction as an aggregate (internal) node.
      INTEGER, PARAMETER, PRIVATE :: EVAL_IN_SPACE = -1 !! This evaluation is NOT taking place at an element
      DOUBLE PRECISION, PARAMETER, PRIVATE :: WORK_EVAL = 100.0d0
      DOUBLE PRECISION, PARAMETER, PRIVATE :: WORK_NODE = 4.0d0
      DOUBLE PRECISION, PARAMETER, PRIVATE :: WORK_TILE = 1.0d0

      ! Magnetization synchronization
      INTEGER, PRIVATE :: ntile_shar
      LOGICAL, PRIVATE :: ldosync !! True if synchronization logic is necessary (TRUE for MPI)
      INTEGER, ALLOCATABLE, PRIVATE :: sync_M_displs(:)  !! (nnode_shar) Offset of master rank during inter-node communication
      INTEGER, ALLOCATABLE, PRIVATE :: sync_M_rcounts(:) !! (nnode_shar) Counts of master rank during inter-node communication
      INTEGER, ALLOCATABLE, PRIVATE :: sync_M_unpack_map(:) !! (ntet) Maps index in M_sync_buf back to global M
      DOUBLE PRECISION, POINTER, PRIVATE :: sync_M_sendbuf(:) !! Buffer where all ranks in comm_shar populate their changed magnetizations
      DOUBLE PRECISION, POINTER, PRIVATE :: sync_M_recvbuf(:) !! Buffer where all ranks in comm_shar can read the remaining changed magnetizations
      INTEGER, ALLOCATABLE, PRIVATE :: M_offset_shar(:) !! (nnode_loc) Offset of inode_loc in sync_node_sendbuf
      INTEGER, PARAMETER :: M_bufsize = 3 !! 3-D magnetization
      INTEGER, PRIVATE :: win_sync_M_sendbuf, win_sync_M_recvbuf

      ! Sparse M sync
      INTEGER, ALLOCATABLE, PRIVATE :: sparse_M_send_counts(:)  ! master_size
      INTEGER, ALLOCATABLE, PRIVATE :: sparse_M_send_displs(:)  ! master_size
      INTEGER, ALLOCATABLE, PRIVATE :: sparse_M_recv_counts(:)  ! master_size
      INTEGER, ALLOCATABLE, PRIVATE :: sparse_M_recv_displs(:)  ! master_size
      INTEGER, POINTER, PRIVATE :: sparse_M_send_map(:)         ! tiles to send
      INTEGER, POINTER, PRIVATE :: sparse_M_recv_map(:)         ! tiles to receive
      DOUBLE PRECISION, POINTER, PRIVATE :: sparse_M_send_buf(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: sparse_M_recv_buf(:)
      INTEGER, PRIVATE :: sparse_M_pack_start, sparse_M_pack_end
      INTEGER, PRIVATE :: sparse_M_unpack_start, sparse_M_unpack_end
      INTEGER, PRIVATE :: win_sparse_M_send_map, win_sparse_M_recv_map
      INTEGER, PRIVATE :: win_sparse_M_send_buf, win_sparse_M_recv_buf

      ! Node synchronization
      INTEGER, PRIVATE :: nnode_shar
      INTEGER, ALLOCATABLE, PRIVATE :: nnodes_per_level_loc(:) !! (tree_depth) Number of nodes at each level for an MPI rank
      INTEGER, ALLOCATABLE, PRIVATE :: nnodes_per_level_shar(:) !! (tree_depth) Number of nodes at each level in comm_shar
      INTEGER, ALLOCATABLE, PRIVATE :: nodes_level_loc(:,:) !! (max(nnodes_per_level_loc),tree_depth) **LOCAL INDICES** of nodes at each level
      INTEGER, ALLOCATABLE, PRIVATE :: sync_node_displs(:,:) !! (max(nnodes_per_level_loc,tree_depth) Offset of master rank during inter-node communication
      INTEGER, ALLOCATABLE, PRIVATE :: sync_node_rcounts(:,:) !! (max(nnodes_per_level_loc,tree_depth) Counts of master rank during inter-node communication
      INTEGER, ALLOCATABLE, PRIVATE :: sync_node_unpack_map(:,:) !! (max(nnodes_per_level_loc,tree_depth) Maps entry i to global node index
      DOUBLE PRECISION, POINTER, PRIVATE :: sync_node_sendbuf(:) !! Buffer where all ranks in comm_shar populate their changed node values
      DOUBLE PRECISION, POINTER, PRIVATE :: sync_node_recvbuf(:) !! Buffer where all ranks in comm_shar can read the remaining changed node values
      INTEGER, ALLOCATABLE, PRIVATE :: node_offset_shar(:) !! (nnode_loc) Offset of inode_loc in sync_node_sendbuf
      INTEGER, PARAMETER, PRIVATE :: node_bufsize = 15 !! Dipole moment-3, position moments-3, weights-1, center of moment-3, quadrupole - 5
      INTEGER, PRIVATE :: split_level
      INTEGER, PRIVATE :: win_sync_node_sendbuf, win_sync_node_recvbuf

      ! Sparse sync infrastructure
      INTEGER, PRIVATE :: shar_id  !! which compute node this rank belongs to
      INTEGER, ALLOCATABLE, PRIVATE :: sparse_send_counts(:)   ! master_size - how many nodes to send to each
      INTEGER, ALLOCATABLE, PRIVATE :: sparse_send_displs(:)   ! master_size - offsets into send buffer
      INTEGER, POINTER, PRIVATE :: sparse_send_map(:)      ! nodes to send, in order
      INTEGER, ALLOCATABLE, PRIVATE :: sparse_recv_counts(:)   ! master_size - how many nodes to recv from each
      INTEGER, ALLOCATABLE, PRIVATE :: sparse_recv_displs(:)   ! master_size - offsets into recv buffer
      INTEGER, POINTER, PRIVATE :: sparse_recv_map(:)      ! where to unpack received nodes
      DOUBLE PRECISION, POINTER, PRIVATE :: sparse_send_buf(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: sparse_recv_buf(:)
      INTEGER, PRIVATE :: sparse_pack_start, sparse_pack_end
      INTEGER, PRIVATE :: sparse_unpack_start, sparse_unpack_end
      INTEGER, PRIVATE :: win_sparse_send_map, win_sparse_recv_map
      INTEGER, PRIVATE :: win_sparse_send_buf, win_sparse_recv_buf

      ! Grid evaluation
      DOUBLE PRECISION, PRIVATE :: theta_eval
      DOUBLE PRECISION, PRIVATE :: theta2_eval
      INTEGER, PRIVATE :: neval_boxes_loc
      INTEGER, PRIVATE :: max_pts_per_box
      INTEGER, PRIVATE :: box_nx, box_ny, box_nz
      DOUBLE PRECISION, PRIVATE :: box_xmin, box_ymin, box_zmin
      DOUBLE PRECISION, PRIVATE :: eval_box_edge
      INTEGER, ALLOCATABLE, PRIVATE :: box_npts(:)
      INTEGER, ALLOCATABLE, PRIVATE :: box_pts_ptr(:)
      INTEGER, ALLOCATABLE, PRIVATE :: box_pts(:)
      INTEGER, ALLOCATABLE, PRIVATE :: eval_local_boxes(:)  !! local box indices

      ! Timing
      DOUBLE PRECISION, PRIVATE :: t_propagate = 0.0d0
      DOUBLE PRECISION, PRIVATE :: t_sync_M    = 0.0d0
      DOUBLE PRECISION, PRIVATE :: t_sync_nodes= 0.0d0
      DOUBLE PRECISION, PRIVATE :: t_field     = 0.0d0
      DOUBLE PRECISION, PRIVATE :: t_solve     = 0.0d0

      ! MPI variables
      INTEGER, PRIVATE :: comm_shar !! Shared memory communicator
      INTEGER, PRIVATE :: shar_rank !! Rank in comm_shar
      INTEGER, PRIVATE :: shar_size !! Size of comm_shar
      INTEGER, PRIVATE :: comm_master  !! Communicator of all ranks satisfying shar_rank=0
      INTEGER, PRIVATE :: master_rank !! Rank in comm_master
      INTEGER, PRIVATE :: master_size !! Size of comm_master
      INTEGER, PRIVATE :: comm_world !! Global communicator of all MPI ranks
      INTEGER, PRIVATE :: world_rank !! Rank in comm_world
      INTEGER, PRIVATE :: world_size !! Size of comm_world
      LOGICAL, PRIVATE :: lcomm, lismaster

      INTEGER, PRIVATE :: ierr_mpi

      ! MPI windows
      INTEGER, PRIVATE :: win_vertex, win_tet, win_tet_cen, &
                          win_tet_vol, win_tet_rad,  &
                          win_tet_P, win_tet_D, win_tet_v, &
                          win_state_dex, win_state_type, &
                          win_constant_mu, win_m, win_Mrem, &
                          win_Happ, win_constant_mu_o

      INTEGER, PRIVATE :: win_node_level, win_node_child, win_node_parent, &
                          win_node_bounds, win_node_cen, &
                          win_leaf_size, win_leaf_to_node, win_leaf_tile

      ! verbose and debug variables
      LOGICAL, PRIVATE :: lverb

      ! precomputed constants
      DOUBLE PRECISION, PARAMETER, PRIVATE :: PI = 4.0D0*ATAN(1.0D0) !! Pie
      DOUBLE PRECISION, PARAMETER, PRIVATE :: INVPI = 1.0D0/PI       !! Evil pie
      DOUBLE PRECISION, PARAMETER, PRIVATE :: INV4PI = 1.0D0/(4.0D0*PI) !! Evil four pies
      DOUBLE PRECISION, PRIVATE, PARAMETER :: mu0 = 16.0D-7*ATAN(1.0d0) !! Permeability of free space
      DOUBLE PRECISION, PRIVATE, PARAMETER :: invmu0 = 1.0d0/mu0 !! Evil permeability of free space
      DOUBLE PRECISION, PRIVATE, PARAMETER :: zero5(5) = 0.0d0
      DOUBLE PRECISION, PARAMETER, PRIVATE :: I3(3,3) = RESHAPE(&
                                                        [1.0d0, 0.0d0, 0.0d0, &
                                                        0.0d0, 1.0d0, 0.0d0, &
                                                        0.0d0, 0.0d0, 1.0d0], [3,3]) !! 3x3 identity matrix

      PRIVATE :: CROSS_PRODUCT

      CONTAINS
!------------------------------------------------------------------------------
!     Subroutines
!       Main flow
!         mumaterial_set_comms:     Sets up MPI communicators (optional)
!         mumaterial_load:      Loads magnetic material file and sets up MPI stuff
!         mumaterial_set_vars:      Sets default values
!         mumaterial_setverb:   Sets standard verbosity
!         mumaterial_info:      Prints information to screen
!         mumaterial_init:      Initializes everything, calls iteration subroutine
!         mumaterial_iterate_M: Main calculation loop
!
!       Namelist Routines
!         mumaterial_init_nml:  Initializes the namelist variables
!         mumaterial_read_nml:  Reads the namelist from an iunit
!         mumaterial_write_nml: Write the namelist to a file.
!
!       Python Interface Routines
!         mumaterial_load_serial:   Loads magnetic material file
!         mumaterial_get_nvertex:   Returns the nvertex variable
!         mumaterial_get_ntet:      Returns the ntet variable
!         mumaterial_get_nstate:    Returns the nstate variable
!         mumaterial_get_vertex:    Returns the vertex variable
!         mumaterial_get_tet:       Returns the tet variable
!         mumaterial_get_statedex:  Returns the state_dex variable
!         mumaterial_get_statetype: Returns the state_type variable
!
!------------------------------------------------------------------------------

      PURE FUNCTION CROSS_PRODUCT(a, b)
!!-----------------------------------------------------------------------
!! Calculates the cross product of 3-vectors a and b.
!!----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN), DIMENSION(3) :: a, b
      DOUBLE PRECISION, DIMENSION(3) :: CROSS_PRODUCT

      CROSS_PRODUCT(1) = a(2)*b(3) - a(3)*b(2)
      CROSS_PRODUCT(2) = a(3)*b(1) - a(1)*b(3)
      CROSS_PRODUCT(3) = a(1)*b(2) - a(2)*b(1)

      
      END FUNCTION CROSS_PRODUCT
!------------------------------------------------------------------------------
! mumaterial_init_nml: Initializes the mumat_input namelist variables
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_init_nml()
      IMPLICIT NONE

      maxiter = 100
      eps_max   = 1.0D-5
      lambdaStart = 0.7
      lambdaFactor = 0.75
      lambdaThresh = 10
      padFactor = 1.0
      Pconv_min = 99.0
      tree_max_depth = 30
      leaf_max_size = 4
      tree_theta = 0.5d0

      END SUBROUTINE mumaterial_init_nml

!------------------------------------------------------------------------------
! mumaterial_read_nml: Reads Mumaterial namelist from file
!------------------------------------------------------------------------------
! param[in]: filename. File containting mumat_input namelist
! param[inout]: istat. Status Flag
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_read_nml(filename, istat)

      IMPLICIT NONE

      CHARACTER(*), INTENT(in)  :: filename
      INTEGER,      INTENT(out) :: istat

      LOGICAL :: lexist
      INTEGER :: iunit
      CHARACTER(LEN=1000) :: line

      NAMELIST /mumat_input/ maxiter, eps_max, lambdaStart, lambdaFactor, lambdaThresh, padFactor, Pconv_min, tree_max_depth, leaf_max_size, tree_theta

      istat = 0
      iunit = 422
      INQUIRE(FILE=TRIM(filename),EXIST=lexist)
      IF (.not.lexist) STOP "Error: Could not find MUMATERIAL namelist file."
      CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
      IF (istat /= 0) THEN
            WRITE(6,'(A)') 'MUMAT error opening file: ',TRIM(filename)
            CALL FLUSH(6)
            RETURN
      END IF
      READ(iunit,NML=mumat_input,IOSTAT=istat)
      IF (istat /= 0) THEN
         WRITE(6,'(A)') 'ERROR reading namelist MUMAT_INPUT from file: ',TRIM(filename)
         backspace(iunit)
         read(iunit,fmt='(A)') line
         write(6,'(A)') 'Invalid line in namelist: '//TRIM(line)
         CALL FLUSH(6)
         STOP
      END IF

      CLOSE(iunit)

      

      END SUBROUTINE mumaterial_read_nml

!------------------------------------------------------------------------------
! mumaterial_write_nml: Writes Mumaterial namelist to a file
!------------------------------------------------------------------------------
! param[inout]: iunit_out. Unit number to write to.
! param[out]: istat. Status flag.
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_write_nml(iunit_out, istat)

      IMPLICIT NONE

      INTEGER,      INTENT(inout)  :: iunit_out
      INTEGER,      INTENT(out) :: istat

      CHARACTER(LEN=*), PARAMETER :: outint  = "(2X,A,1X,'=',1X,I0)"
      CHARACTER(LEN=*), PARAMETER :: outflt  = "(2X,A,1X,'=',1X,ES22.12E3)"

      WRITE(iunit_out,'(A)') '&MUMAT_INPUT'
      WRITE(iunit_out,outint) 'MAXITER',maxiter
      WRITE(iunit_out,outflt) 'DMMAX',eps_max
      WRITE(iunit_out,outflt) 'LAMBDASTART',lambdaStart
      WRITE(iunit_out,outflt) 'LAMBDAFACTOR',lambdaFactor
      WRITE(iunit_out,outint) 'LAMBDATHRESH',lambdaThresh
      WRITE(iunit_out,outflt) 'PADFACTOR',padFactor
      WRITE(iunit_out,outflt) 'CONVCHECK',Pconv_min
      WRITE(iunit_out,'(A)') '/'
        istat = 0 ! Just zero for now
      

      END SUBROUTINE mumaterial_write_nml

!------------------------------------------------------------------------------
! mumaterial_write_nml_byfile: Writes Mumaterial namelist to a file
!------------------------------------------------------------------------------
! param[in]: filename. File to read from
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_write_nml_byfile(filename)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER :: iunit, istat
      LOGICAL :: lexists

      iunit = 100
      istat = 0
      INQUIRE(FILE=TRIM(filename),exist=lexists)
      IF (lexists) THEN
         OPEN(unit=iunit, file=TRIM(filename), iostat=istat, status="old", position="append")
      ELSE
         OPEN(unit=iunit, file=TRIM(filename), iostat=istat, status="new")
      END IF
      IF (istat .ne. 0) RETURN
      CALL mumaterial_write_nml(iunit,istat)
      CLOSE(iunit)

      

      END SUBROUTINE mumaterial_write_nml_byfile


!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_setverb(lverbin)
!!------------------------------------------------------------------------------
!! Sets verbosity.
!!------------------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: lverbin !! Verbose or not
      lverb = lverbin
      END SUBROUTINE mumaterial_setverb

!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_set_Bfld(func_B)
!!------------------------------------------------------------------------------
!! Associates external B field function ext_Bfld with func_B
!!------------------------------------------------------------------------------
      IMPLICIT NONE
      PROCEDURE(externalFieldFunc) :: func_B !! Pointer to external function
      ext_Bfld => func_B
      END SUBROUTINE mumaterial_set_Bfld

!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_set_vars(pad_factor,max_error, max_iter, min_conv_perc, &
        lambda_factor, lambda_start, lambda_threshold,lambda_max,lambda_min, &
        max_depth, max_leafsize, iter_theta, eval_theta)
!!------------------------------------------------------------------------------
!! Sets any user variables.
!!------------------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in), OPTIONAL :: max_error !! Threshold 
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: lambda_start
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: lambda_factor
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: lambda_max
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: lambda_min
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: pad_factor
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: min_conv_perc
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: iter_theta
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: eval_theta
      INTEGER, INTENT(in), OPTIONAL :: max_iter
      INTEGER, INTENT(in), OPTIONAL :: lambda_threshold
      INTEGER, INTENT(in), OPTIONAL :: max_depth
      INTEGER, INTENT(in), OPTIONAL :: max_leafsize

      IF (PRESENT(pad_factor))        padFactor = pad_factor
      IF (PRESENT(max_error))         eps_max = max_error
      IF (PRESENT(max_iter))          maxiter = max_iter
      IF (PRESENT(min_conv_perc))     Pconv_min = min_conv_perc

      IF (PRESENT(lambda_start))      lambdaStart  = lambda_start
      IF (PRESENT(lambda_factor))     lambdaFactor = lambda_factor
      IF (PRESENT(lambda_threshold))  lambdaThresh = lambda_threshold ! unused
      IF (PRESENT(lambda_max))        lambdaMax = lambda_max
      IF (PRESENT(lambda_min))        lambdaMin = lambda_min


      IF (PRESENT(max_depth))         tree_max_depth = max_depth
      IF (PRESENT(max_leafsize))      leaf_max_size = max_leafsize

      IF (PRESENT(iter_theta)) THEN
                              tree_theta = iter_theta
                              tree_theta2 = iter_theta*iter_theta
      ENDIF 
      IF (PRESENT(eval_theta)) THEN
                              theta_eval = eval_theta
                              theta2_eval = eval_theta*eval_theta
      ENDIF 
      

      END SUBROUTINE mumaterial_set_vars

      SUBROUTINE mumaterial_set_comms(comm, comm_shar_out, comm_master_out)
!!------------------------------------------------------------------------------
!! Set up MUMAT communicators
!!------------------------------------------------------------------------------
      USE mpi_inc
      IMPLICIT NONE

      INTEGER, INTENT(inout) :: comm !! World communicator
      INTEGER, INTENT(out) :: comm_shar_out !! Shared memory communicator
      INTEGER, INTENT(out) :: comm_master_out !! Master communicator
      INTEGER :: comm_myworld, i

      CALL MPI_COMM_DUP( comm, comm_myworld, ierr_mpi )
      CALL MPI_COMM_SPLIT_TYPE( comm_myworld, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, comm_shar_out, ierr_mpi)
      CALL MPI_COMM_RANK( comm_shar_out, shar_rank, ierr_mpi)

      i = MPI_UNDEFINED
      IF (shar_rank.EQ.0) i = 0
      CALL MPI_COMM_SPLIT( comm_myworld, i, shar_rank, comm_master_out, ierr_mpi )

      END SUBROUTINE mumaterial_set_comms

!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_free()
!!------------------------------------------------------------------------------
!! Deallocates memory and destroys MPI windows
!!------------------------------------------------------------------------------
      USE mpi_inc
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
      IF (ASSOCIATED(tet_P))         CALL mpidealloc(tet_P,win_tet_P)
      IF (ASSOCIATED(tet_D))         CALL mpidealloc(tet_D,win_tet_D)
      IF (ASSOCIATED(tet_v))         CALL mpidealloc(tet_v,win_tet_v)
      IF (ASSOCIATED(M))             CALL mpidealloc(M,win_M)
      IF (ASSOCIATED(Mrem))          CALL mpidealloc(Mrem,win_Mrem)

      DO ik = 1, nstate
         IF (ALLOCATED(stateFunction(ik)%H)) DEALLOCATE(stateFunction(ik)%H)
         IF (ALLOCATED(stateFunction(ik)%M)) DEALLOCATE(stateFunction(ik)%M)
         IF (ALLOCATED(stateFunction(ik)%b)) DEALLOCATE(stateFunction(ik)%b)
      END DO
      IF (ALLOCATED(stateFunction)) DEALLOCATE(stateFunction)

      IF (ASSOCIATED(node_level))     CALL mpidealloc(node_level,win_node_level)
      IF (ASSOCIATED(node_child))     CALL mpidealloc(node_child,win_node_child)
      IF (ASSOCIATED(node_parent))    CALL mpidealloc(node_parent,win_node_parent)
      IF (ASSOCIATED(node_bounds))    CALL mpidealloc(node_bounds,win_node_bounds)
      IF (ASSOCIATED(node_cen))       CALL mpidealloc(node_cen,win_node_cen)
      IF (ASSOCIATED(leaf_size))      CALL mpidealloc(leaf_size,win_leaf_size)
      IF (ASSOCIATED(leaf_tile))      CALL mpidealloc(leaf_tile,win_leaf_tile)
      IF (ASSOCIATED(leaf_to_node))   CALL mpidealloc(leaf_to_node,win_leaf_to_node)

      IF (ASSOCIATED(node_m))         CALL mpidealloc(node_m,win_node_m)
      IF (ASSOCIATED(node_p))         CALL mpidealloc(node_p,win_node_p)
      IF (ASSOCIATED(node_com))       CALL mpidealloc(node_com,win_node_com)
      IF (ASSOCIATED(node_w))         CALL mpidealloc(node_w,win_node_w)
      IF (ASSOCIATED(node_q))         CALL mpidealloc(node_q,win_node_q)

      IF (ASSOCIATED(sync_node_sendbuf))  CALL mpidealloc(sync_node_sendbuf,win_sync_node_sendbuf)
      IF (ASSOCIATED(sync_node_recvbuf))  CALL mpidealloc(sync_node_recvbuf,win_sync_node_recvbuf)
      IF (ASSOCIATED(sync_M_sendbuf))     CALL mpidealloc(sync_M_sendbuf,win_sync_M_sendbuf)
      IF (ASSOCIATED(sync_M_recvbuf))     CALL mpidealloc(sync_M_recvbuf,win_sync_M_recvbuf)

      IF (ASSOCIATED(sparse_send_buf))    CALL mpidealloc(sparse_send_buf,win_sparse_send_buf)
      IF (ASSOCIATED(sparse_recv_buf))    CALL mpidealloc(sparse_recv_buf,win_sparse_recv_buf)
      IF (ASSOCIATED(sparse_send_map))    CALL mpidealloc(sparse_send_map,win_sparse_send_map)
      IF (ASSOCIATED(sparse_recv_map))    CALL mpidealloc(sparse_recv_map,win_sparse_recv_map)

      IF (ASSOCIATED(sparse_M_send_map))  CALL mpidealloc(sparse_M_send_map,win_sparse_M_send_map)
      IF (ASSOCIATED(sparse_M_recv_map))  CALL mpidealloc(sparse_M_recv_map,win_sparse_M_recv_map)
      IF (ASSOCIATED(sparse_M_send_buf))  CALL mpidealloc(sparse_M_send_buf,win_sparse_M_send_buf)
      IF (ASSOCIATED(sparse_M_recv_buf))  CALL mpidealloc(sparse_M_recv_buf,win_sparse_M_recv_buf)   
      
      END SUBROUTINE mumaterial_free

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
      USE mpi_inc
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
      shar_rank = 0; shar_size = 1;
      master_rank = 0; master_size = 1;
      world_rank = 0; world_size = 1
      lismaster = .TRUE.; ldosync = .FALSE.

      ! Set up MPI parameters properly now
#if defined(MPI_OPT)
      IF (lcomm) THEN
        lismaster = .FALSE.; master_rank = 1
        CALL MPI_COMM_RANK( comm_world, world_rank, ierr_mpi)
        CALL MPI_COMM_SIZE( comm_world, world_size, ierr_mpi)
        CALL MPI_COMM_RANK( comm_shar,  shar_rank,  ierr_mpi)
        CALL MPI_COMM_SIZE( comm_shar,  shar_size,  ierr_mpi)
        shar_id = world_rank / shar_size
        IF (shar_rank.eq.0) THEN
          CALL MPI_COMM_RANK( comm_master, master_rank, ierr_mpi )
          CALL MPI_COMM_SIZE( comm_master, master_size, ierr_mpi )
          lismaster = (master_rank.EQ.0)
        END IF
        CALL MPI_Bcast( master_size, 1, MPI_INTEGER, 0, comm_shar, ierr_mpi)
        ldosync = (master_size>1)      
      END IF
#endif

      ! Nullify pointers
      NULLIFY(vertex, tet, tet_cen, tet_vol, state_dex, state_type, &
              constant_mu, constant_mu_o, Mrem, M, tet_P, tet_D, tet_V)

      ! open file, return if fails
      iunit = 327; istat = 0
      CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
      IF (istat/= 0) RETURN
      ! master reads info
      IF (lismaster) THEN
         READ(iunit,'(A)') mesh_name
         READ(iunit,'(A)') mesh_date
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
        CALL mpialloc(tet_P,3,3,4,ntet,    shar_rank,0,comm_shar,win_tet_P)
        CALL mpialloc(tet_D,3,4,ntet,      shar_rank,0,comm_shar,win_tet_D)
        CALL mpialloc(tet_v,3,3,4,ntet,    shar_rank,0,comm_shar,win_tet_v)
        CALL mpialloc(state_dex,ntet,      shar_rank,0,comm_shar,win_state_dex)
        CALL mpialloc(state_type,nstate,   shar_rank,0,comm_shar,win_state_type)
        CALL mpialloc(constant_mu,nstate,  shar_rank,0,comm_shar,win_constant_mu)
        CALL mpialloc(constant_mu_o,nstate,shar_rank,0,comm_shar,win_constant_mu_o)
        CALL mpialloc(M,            3,ntet,shar_rank,0,comm_shar,win_m)
        CALL mpialloc(Mrem,3,ntet,         shar_rank,0,comm_shar,win_Mrem)  ! TODO: Allocate locally
        ALLOCATE(stateFunction(nstate))
      ELSE
#endif
         ! if no MPI, allocate everything on one node
         ALLOCATE(vertex(3,nvertex),tet(4,ntet),state_dex(ntet), &
                  state_type(nstate),constant_mu(nstate), &
                  tet_cen(3,ntet),tet_vol(ntet),M(3,ntet), &
                  constant_mu_o(nstate),Mrem(3,ntet),stateFunction(nstate), &
                  STAT=istat)
          ALLOCATE(tet_P(3,3,4,ntet))
          ALLOCATE(tet_D(3,4,ntet))
          ALLOCATE(tet_v(3,3,4,ntet))
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
            IF (state_type(ik) == STATE_HARD) THEN
               READ(iunit,*) constant_mu(ik), constant_mu_o(ik)
               READ(iunit,*) Mrem(1,ik),Mrem(2,ik),Mrem(3,ik)
            ELSEIF (state_type(ik) == STATE_SOFT) THEN
               READ(iunit,*) nMH
               ALLOCATE(stateFunction(ik)%H(nMH), &
                        stateFunction(ik)%M(nMH), &
                        stateFunction(ik)%b(nMH))
               READ(iunit,*) stateFunction(ik)%H(:)
               READ(iunit,*) stateFunction(ik)%M(:)
               CALL get_spline_slopes(stateFunction(ik)%H,&
                                               stateFunction(ik)%M,&
                                               stateFunction(ik)%b)
            ELSEIF (state_type(ik) == STATE_LINEAR) THEN
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
              IF (master_rank .ne. 0) THEN
                ALLOCATE(stateFunction(ik)%H(nMH), &
                         stateFunction(ik)%M(nMH), &
                         stateFunction(ik)%b(nMH))
              END IF
              CALL MPI_Bcast(stateFunction(ik)%H,nMH,MPI_DOUBLE_PRECISION,0,comm_master,ierr_mpi)
              CALL MPI_Bcast(stateFunction(ik)%M,nMH,MPI_DOUBLE_PRECISION,0,comm_master,ierr_mpi)
              CALL MPI_Bcast(stateFunction(ik)%b,nMH,MPI_DOUBLE_PRECISION,0,comm_master,ierr_mpi)
            END IF
          END IF
          ! Then from submasters to other threads
          CALL MPI_Bcast(nMH,1,MPI_INTEGER,0,comm_shar,ierr_mpi)
          IF (nMH .gt. 0) THEN
            IF (shar_rank .ne. 0) THEN
              ALLOCATE(stateFunction(ik)%H(nMH), &
                       stateFunction(ik)%M(nMH), &
                       stateFunction(ik)%b(nMH))
            END IF
            CALL MPI_Bcast(stateFunction(ik)%H,nMH,MPI_DOUBLE_PRECISION,0,comm_shar,ierr_mpi)
            CALL MPI_Bcast(stateFunction(ik)%M,nMH,MPI_DOUBLE_PRECISION,0,comm_shar,ierr_mpi)
            CALL MPI_Bcast(stateFunction(ik)%b,nMH,MPI_DOUBLE_PRECISION,0,comm_shar,ierr_mpi)
          END IF
        END DO
      END IF
#endif

      CLOSE(iunit) ! close file
      
      CONTAINS

      SUBROUTINE get_spline_slopes(fx, fy, b)
!!-----------------------------------------------------------------------
!! Gets spline slopes by solving a tridiagonal matrix
!!-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: fx, fy
      DOUBLE PRECISION, DIMENSION(:), INTENT(out) :: b
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: h, delta, dl, d, du
      DOUBLE PRECISION :: alpha, beta, s2,  tau
      INTEGER :: n, ii, INFO

      n = SIZE(fx)
      ALLOCATE(h(n-1),delta(n-1),dl(n-1),d(n),du(n-1))

      DO ii = 1, n-1
        h(ii) = fx(ii+1) - fx(ii)
        delta(ii) = (fy(ii+1)-fy(ii))/h(ii)
      END DO

      ! Lower diagonal
      dl(1) = 1.0d0
      DO ii = 2, n-1
        dl(ii) = h(ii)
      END DO

      ! Main diagonal
      d(1) = 2.0d0
      DO ii = 2, n-1
        d(ii) = 2.0d0 * (h(ii-1)+h(ii))
      END DO
      d(n) = 2.0d0

      ! Upper diagonal
      DO ii = 2, n-1
        du(ii) = h(ii-1)
      END DO
      du(n-1) = 1.0d0

      ! RHS
      b(1) = 3.0d0*delta(1)
      DO ii = 2, n-1
        b(ii) = 3.0d0*(h(ii)*delta(ii-1)+h(ii-1)*delta(ii))
      END DO
      b(n) = 3.0d0*delta(n-1)

      ! LAPACK
      CALL DGTSV(n, 1, dl, d, du, b, n, INFO)
      ! Force tangents to zero if they point opposite to local slope
      DO ii = 1, n
        IF (ii == 1) THEN
          IF (b(ii) * delta(1) < 0.0d0) b(ii) = 0.0d0
        ELSE IF (ii == n) THEN
          IF (b(ii) * delta(n-1) < 0.0d0) b(ii) = 0.0d0
        ELSE
          IF (b(ii) * delta(ii-1) < 0.0d0 .OR. b(ii) * delta(ii) < 0.0d0) b(ii) = 0.0d0
        END IF
      END DO

      ! Monotonicity correction (Fritsch-Carlson)
      DO ii = 1, n-1
        IF (ABS(delta(ii)) .LT. TINY(1.0d0)) THEN
            b(ii) = 0.0d0
            b(ii+1) = 0.0d0
        ELSE
            alpha = b(ii) / delta(ii)
            beta  = b(ii+1) / delta(ii)
            s2 = alpha**2 + beta**2
            IF (s2 .GT. 9.0d0) THEN
                tau = 3.0d0 / SQRT(s2)
                b(ii)   = tau * alpha * delta(ii)
                b(ii+1) = tau * beta  * delta(ii)
            END IF
        END IF
      END DO
      DEALLOCATE(h, delta, dl, d, du)

      

      END SUBROUTINE get_spline_slopes
        
      END SUBROUTINE mumaterial_load

!------------------------------------------------------------------------------
! mumaterial_load_serial: Loads magnetic material file (no MPI for Python)
!------------------------------------------------------------------------------
! param[in]: filename. The file name to load in
! param[in, out]: istat. Integer that shows  if != 0
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_load_serial(filename,istat)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER, INTENT(inout)       :: istat

      CALL mumaterial_load(filename,istat)

      

      END SUBROUTINE mumaterial_load_serial


! param[in]: iunit. Unit number to print to
! param[in]: lnoiter: Whether to do mumat iterations
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_info(iunit, lnoiter)
!!------------------------------------------------------------------------------
!! Prints mumat info to iunit
!!------------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iunit
      LOGICAL, INTENT(IN) :: lnoiter !! Whether to skip iterations
      INTEGER :: istate,k

      WRITE(iunit,'(A)')           ' ---------- MUMAT MPI ----------'
      WRITE(iunit,'(3X,A,I7)')     '  MPI Nodes    : ',master_size
      WRITE(iunit,'(3X,A,I7)')     '  MPI Threads  : ',world_size
      WRITE(iunit,'(A)')           ' -----  Magnetic Material  -----'
      WRITE(iunit,'(3X,A,A)')      '  Model Name   : ',TRIM(mesh_name)
      WRITE(iunit,'(3X,A,A)')      '  Date         : ',TRIM(mesh_date)
      WRITE(iunit,'(3X,A,I7)')     '  Vertices     : ',nvertex
      WRITE(iunit,'(3X,A,I7)')     '  Tetrahedrons : ',ntet
      WRITE(iunit,'(3X,A,I7)')     '  State Funcs. : ',nstate

      DO istate = 1, nstate
        SELECT CASE (state_type(istate))
          CASE (STATE_HARD) 
            WRITE(iunit,'(6X,I3,A)') istate, '. Hard Magnet'
            WRITE(iunit,'(9X,A,EN12.3)')    '├ mu   :',constant_mu(istate)
            WRITE(iunit,'(9X,A,EN12.3)')    '├ mu_o :',constant_mu_o(istate)
            WRITE(iunit,'(9X,A,3(EN12.3))') '└ Mrem :',Mrem(:,istate)
          CASE (STATE_SOFT)
            k = SIZE(stateFunction(istate)%H)
            WRITE(iunit,'(6X,I3,A)') istate, '. Soft material (M(H))'
            WRITE(iunit,'(9X,A,I3)')        '├ NKnots:',k
            WRITE(iunit,'(9X,A,2(EN12.3))') '├──── H :',stateFunction(istate)%H(1),stateFunction(istate)%H(k)
            WRITE(iunit,'(9X,A,2(EN12.3))') '└──── M :',stateFunction(istate)%M(1),stateFunction(istate)%M(k)
          CASE (STATE_LINEAR) 
            WRITE(iunit,'(6X,I3,A)') istate,'. Linear material'
            WRITE(iunit,'(9X,A,EN12.3)')    '└ mu    :',constant_mu(istate)
          CASE DEFAULT
            WRITE(iunit,'(6X,I3,A,I3)') istate,'. UNKNOWN STATE TYPE:',state_type(istate)
        END SELECT
      END DO

      IF (.NOT.lnoiter) THEN
        WRITE(iunit,'(A)')           ' --------  Iterations  ---------'
        WRITE(iunit,'(3X,A,I12)')     '  Max iter. #  : ',maxiter
        WRITE(iunit,'(3X,A,EN12.3)') '  Max error    :',eps_max
        WRITE(iunit,'(3X,A,F12.6)')   '  Lambda start : ',lambdaStart
        WRITE(iunit,'(3X,A,F12.6)')   '  Lambda min.  : ',lambdaMin
        WRITE(iunit,'(3X,A,F12.6)')   '  Lambda max.  : ',lambdaMax
        WRITE(iunit,'(3X,A,F11.5)')   '  Converged at :  ',Pconv_min
      END IF
      WRITE(iunit,'(A)') ' -------------------------------'
      FLUSH(iunit)

      END SUBROUTINE mumaterial_info

     
!-----------------------------------------------------------------------------
      SUBROUTINE mumaterial_run(Bfld,  offset, lskip, linitM)
!!------------------------------------------------------------------------------
!! Start MUMAT
!!------------------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in), OPTIONAL :: offset(3) !! offset in m
      LOGICAL, INTENT(in), OPTIONAL :: lskip !! True to skip iterations
      LOGICAL, INTENT(in), OPTIONAL :: linitM !! True to initialize M from Happ
      LOGICAL :: lskip_loc, linitM_loc

      EXTERNAL:: Bfld
      lskip_loc = .FALSE.; IF (PRESENT(lskip)) lskip_loc = lskip
      linitM_loc = .FALSE.; IF(PRESENT(linitM)) linitM_loc = linitM

      CALL mumaterial_set_Bfld(Bfld)
      CALL mumaterial_init_mesh(offset)

      ! Build tree
      CALL mumaterial_init_tree()
      CALL mumaterial_init_demag()   ! Get N for nearby leaf-leaf pairs

      ! Synchronization setup
#if defined(MPI_OPT)
      IF (ldosync) THEN
        CALL mumaterial_init_sync_M()
        CALL mumaterial_init_sync_M_sparse()
        CALL mumaterial_init_sync_node()
        CALL mumaterial_init_sync_node_sparse()
      END IF
#endif
      CALL mumaterial_init_Happ(linitM_loc)  ! Calculate static background field at local elements
      CALL mumaterial_propagate_nodes(lverb) ! Calculate moments at every node

      ! Finally, run
      IF (.NOT.lskip_loc) THEN
        IF (lverb) WRITE (6,*) "  MUMAT:  Beginning Iterations"
        CALL mumaterial_iterate_M()
      ELSE
        IF (lverb) WRITE (6,*) "  MUMAT:  Skipping iterations"
      END IF

      END SUBROUTINE mumaterial_run

!-----------------------------------------------------------------------
      SUBROUTINE mumaterial_init_mesh(offset)
!!-----------------------------------------------------------------------
!! Offsets vertices, calculates centroid, volume, inradius, and writes
!! everything to shared memory and synchronizes across all ranks.
!!-----------------------------------------------------------------------
      USE mpi_inc
#if defined (MPI_OPT)
      USE mpi_params, ONLY : MPI_CALC_MYRANGE
#endif
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in), OPTIONAL :: offset(3) !! Offset (m)
      INTEGER :: mystart, myend, i
      !-----------------------------------------
      ! Offset vertices and synchronize
      !-----------------------------------------
      IF (PRESENT(offset)) THEN
        IF(NORM2(offset) .GT. 0.d0) THEN
          IF (lverb) WRITE(6,*) "MUMAT_INIT: Offseting vertices"
          IF (shar_rank==0) THEN
            DO i = 1, nvertex
              vertex(:,i) = vertex(:,i) + offset
            END DO
          END IF
        END IF
#if defined(MPI_OPT)
        CALL MPI_BARRIER(comm_shar, ierr_mpi)
#endif
      END IF
      !-----------------------------------------
      ! Calculate center, volume, inradius
      !-----------------------------------------
      IF (shar_rank.EQ.0) THEN  ! Master wipes shared arrays
        tet_cen = 0; tet_vol = 0
      END IF
      CALL MPI_BARRIER(comm_shar, ierr_mpi) ! Calculate range
      mystart = 1; myend = ntet
      IF (lcomm) THEN
#if defined(MPI_OPT)
        CALL MPI_CALC_MYRANGE(comm_world, 1, ntet, mystart, myend)
#endif
      END IF
      IF (lverb) WRITE(6,*) "  MUMAT_INIT: Calculating tetrahedron quantities"

      DO i = mystart, myend
          tet_cen(:,i) = (vertex(:,tet(1,i)) + vertex(:,tet(2,i)) + vertex(:,tet(3,i)) + vertex(:,tet(4,i))) / 4.d0
          tet_vol(i) = GET_VOL(vertex(:,tet(1,i)),vertex(:,tet(2,i)), vertex(:,tet(3,i)),vertex(:,tet(4,i)))
      END DO
      !-----------------------------------------
      ! Synchronize shared arrays
      !-----------------------------------------
#if defined(MPI_OPT)
      CALL MPI_BARRIER(comm_shar, ierr_mpi)
      IF (shar_rank.EQ.0) THEN
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tet_cen, 3*ntet, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tet_vol,   ntet, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
      ENDIF
      CALL MPI_BARRIER(comm_shar, ierr_mpi)
#endif
      tet_vol_tot = SUM(tet_vol)

      CONTAINS

      PURE FUNCTION GET_VOL(v1,v2,v3,v4)
!!-----------------------------------------------------------
!! Calculates the volume of a tetrahedron with vertices v1-v4.
!!------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3), INTENT(in) :: v1, v2, v3, v4 !! vertices
      DOUBLE PRECISION :: GET_VOL
      GET_VOL = ABS(DOT_PRODUCT(v1-v4,CROSS_PRODUCT(v2-v4,v3-v4)))/6.0d0
      END FUNCTION GET_VOL
      END SUBROUTINE mumaterial_init_mesh
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


      SUBROUTINE mumaterial_init_demag()
!!-----------------------------------------------------------------------------
!! Initializes demagnetization tensor
!!-----------------------------------------------------------------------------
      USE mpi_inc
      IMPLICIT NONE
      INTEGER :: nN, iN, iNself, ileaf_loc, itile_loc, itile_src, tile_src, tile_targ, inode
      INTEGER :: tile_loc_idx
      INTEGER :: csr, ptr1, ptr2
      DOUBLE PRECISION :: N(3,3), pos(3,max_leaf_size_seen)

      CALL init_demag() ! Get helpers first

      !---------------------------------
      ! PASS 1: Determine size of Nloc
      !---------------------------------
      IF (lverb) WRITE (6,*) "  MUMAT_INIT:  Setting up N and Nself"
      nN = 0
      DO ileaf_loc = 1, nleaf_loc ! Loop over local leaves
        ptr1 = interact_ptr(ileaf_loc); ptr2 = interact_ptr(ileaf_loc+1)-1
        DO csr = ptr1, ptr2  ! Loop over interactions
          IF (interact_type(csr) == INTERACT_LEAF) THEN ! Open source leaves
            inode = interact_list(csr)
            nN = nN + leaf_size(inode)*leaf_size_loc(ileaf_loc)
          END IF
        END DO
      END DO
      nN = nN - ntile_loc
      ! Allocate Nloc
      ALLOCATE(Nloc(nN,3,3),Nself_loc(ntile_loc,3,3))
      ALLOCATE(nN_per_leaf(nleaf_loc))
      nN_per_leaf = 0
      !---------------------------------
      ! PASS 2: Fill Nloc
      !---------------------------------
      iN = 0
      iNself = 0
      DO ileaf_loc = 1, nleaf_loc ! Loop over local leaves
        nN = 0
        DO itile_loc = 1, leaf_size_loc(ileaf_loc) ! Populate positions
          pos(:,itile_loc) = tet_cen(:, leaf_tile_loc(itile_loc, ileaf_loc))
        END DO
        ptr1 = interact_ptr(ileaf_loc); ptr2 = interact_ptr(ileaf_loc+1)-1
        DO csr = ptr1, ptr2 ! Loop over interactions
          IF (interact_type(csr) == INTERACT_LEAF) THEN ! Open source leaves
            inode = interact_list(csr)
            DO itile_src = 1, leaf_size(inode) ! Loop over source leaf elements
              tile_src = leaf_tile(itile_src, inode)            ! global index
              DO itile_loc = 1, leaf_size_loc(ileaf_loc)
                tile_targ = leaf_tile_loc(itile_loc, ileaf_loc) ! global index
                N = mumaterial_get_N(tet_P(:,:,:,tile_src), &
                          tet_D(:,:,tile_src), &
                          tet_v(:,:,:,tile_src), &
                          pos(:,itile_loc))
                IF (tile_targ==tile_src) THEN ! Self N goes into different array
                  tile_loc_idx = leaf_offset_loc(ileaf_loc) + itile_loc
                  Nself_loc(tile_loc_idx,:,:) = N
                ELSE ! Increment sequentially
                  iN = iN + 1
                  Nloc(iN,:,:) = N
                  nN = nN+1
                END IF
              END DO
            END DO
          END IF
        END DO
        nN_per_leaf(ileaf_loc) = nN
      END DO
      CONTAINS

      SUBROUTINE init_demag()
!!-----------------------------------------------------------------------
!!  Calculates and synchronizes helpers for demagnetization tensor calculation
!!-----------------------------------------------------------------------
      USE mpi_inc
#if defined (MPI_OPT)
      USE mpi_params, ONLY : MPI_CALC_MYRANGE
#endif
      IMPLICIT NONE

      INTEGER :: i_tile, mystart, myend
      DOUBLE PRECISION :: v(3,4)

      IF (lverb) WRITE (6,*) "  MUMAT_INIT:  Calculating demag helpers"
      mystart = 1; myend = ntet
#if defined(MPI_OPT)
      IF (lcomm) THEN
        CALL MPI_CALC_MYRANGE(comm_world, 1, ntet, mystart, myend)
      END IF
#endif
      IF (shar_rank.EQ.0) THEN ! Zero out
            tet_P = 0.0d0; tet_D = 0.0d0; tet_v = 0.0d0
      END IF
#if defined(MPI_OPT)
      CALL MPI_BARRIER(comm_world, ierr_mpi) ! Wait until master is done
#endif
      ! Now populate
      DO i_tile = mystart, myend
        CALL GET_N_HELPERS( vertex(:,tet(1,i_tile)), &
                            vertex(:,tet(2,i_tile)), &
                            vertex(:,tet(3,i_tile)), &
                            vertex(:,tet(4,i_tile)), &
                            tet_P(:,:,:,i_tile), &
                            tet_D(:,:,i_tile), &
                            tet_v(:,:,:,i_tile))
      END DO
      ! Sync
#if defined(MPI_OPT)
        CALL MPI_BARRIER(comm_world, ierr_mpi)
        IF (shar_rank.EQ.0) THEN
          CALL MPI_ALLREDUCE( MPI_IN_PLACE, tet_P,    3*3*4*ntet, MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
          CALL MPI_ALLREDUCE( MPI_IN_PLACE, tet_D,    3*4*ntet,   MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
          CALL MPI_ALLREDUCE( MPI_IN_PLACE, tet_v,  3*3*4*ntet,   MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr_mpi )
        END IF
        CALL MPI_BARRIER(comm_world, ierr_mpi)
#endif

      
      END SUBROUTINE init_demag

      PURE SUBROUTINE GET_N_HELPERS(v1, v2, v3, v4, P, D, v_loc)
!!-----------------------------------------------------------------------
!! Calculates arrays encoding all the necessary information for calculating
!! the demagnetization tensor for any element and any point in space.
!!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in), DIMENSION(3) :: v1, v2, v3, v4 !! vertices
      DOUBLE PRECISION, INTENT(out) :: P(3,3,4) !! Rotation matrix for each face.
      DOUBLE PRECISION, INTENT(out) :: D(3,4)   !!  Position of triangle base for each face.
      DOUBLE PRECISION, INTENT(out) :: v_loc(3,3,4) !! Vertices of each face in right permutation and local face coordinates

      DOUBLE PRECISION :: v_swap(3), cosalpha(3), d12, d13, d23, angles(3), Pinv(3,3)
      DOUBLE PRECISION :: Ptmp(3,3), Dtmp(3), vtmp(3,4)
      INTEGER :: i_f, i_v
      DOUBLE PRECISION, PARAMETER :: v_min = 1.0d-12

      ! Rotate vertices
      DO i_f = 1, 4
        vtmp(:,i_f)            = v1
        vtmp(:,MOD(i_f  ,4)+1) = v2
        vtmp(:,MOD(i_f+1,4)+1) = v3
        vtmp(:,MOD(i_f+2,4)+1) = v4

        ! todo: ensure vertices are not colinear and v4 is not in plane of v1-3?

        ! Ensure largest angle is for v2
        d12 = NORM2(vtmp(:,1)-vtmp(:,2))
        d13 = NORM2(vtmp(:,1)-vtmp(:,3))
        d23 = NORM2(vtmp(:,2)-vtmp(:,3))
        cosalpha(1) = MAX(-1.0d0, MIN(1.0d0, DOT_PRODUCT(vtmp(:,1)-vtmp(:,2),vtmp(:,1)-vtmp(:,3)) / (d12*d13)))
        cosalpha(2) = MAX(-1.0d0, MIN(1.0d0, DOT_PRODUCT(vtmp(:,2)-vtmp(:,1),vtmp(:,2)-vtmp(:,3)) / (d12*d23)))
        cosalpha(3) = MAX(-1.0d0, MIN(1.0d0, DOT_PRODUCT(vtmp(:,3)-vtmp(:,2),vtmp(:,3)-vtmp(:,1)) / (d13*d23)))

        IF (cosalpha(1) < cosalpha(2) .and. cosalpha(1) < cosalpha(3)) THEN ! v1 and v2 should be interchanged
          v_swap = vtmp(:,2)
          vtmp(:,2) = vtmp(:,1)
          vtmp(:,1) = v_swap
        ELSE IF (cosalpha(3) < cosalpha(1) .and. cosalpha(3) < cosalpha(2)) THEN ! v2 and v3 should be interchanged
          v_swap = vtmp(:,2)
          vtmp(:,2) = vtmp(:,3)
          vtmp(:,3) = v_swap
        END IF
        ! Ensure normal vector is pointing in the right direction
        IF (DOT_PRODUCT(CROSS_PRODUCT(vtmp(:,1) - vtmp(:,3), vtmp(:,2) - vtmp(:,3)), vtmp(:,4) - vtmp(:,2)) .gt. 0) THEN
          ! normal vector of triangle is pointing towards v4, so v1 and v3 need to be interchanged
          v_swap = vtmp(:,1)
          vtmp(:,1) = vtmp(:,3)
          vtmp(:,3) = v_swap
        END IF

        ! Rotation matrix
        Ptmp = 0.0d0
        Ptmp(:,1) = vtmp(:,1) - vtmp(:,3)
        Ptmp(:,1) = Ptmp(:,1) / NORM2(Ptmp(:,1))
        Ptmp(:,3) = CROSS_PRODUCT(Ptmp(:,1), vtmp(:,2)-vtmp(:,3))
        Ptmp(:,3) = Ptmp(:,3) / NORM2(Ptmp(:,3))
        Ptmp(:,2) = CROSS_PRODUCT(Ptmp(:,3), Ptmp(:,1))
        Ptmp(:,2) = Ptmp(:,2) / NORM2(Ptmp(:,2))
        Pinv = TRANSPOSE(Ptmp)

        ! Position of triangle base
        Dtmp = 0.0d0
        d13 = NORM2(vtmp(:,1)-vtmp(:,3))
        d23 = NORM2(vtmp(:,2)-vtmp(:,3))
        Dtmp = DOT_PRODUCT(vtmp(:,3)-vtmp(:,2),vtmp(:,3)-vtmp(:,1)) &
                            /(d23 * d13) * d23 * Ptmp(:,1) + vtmp(:,3)

        ! Vertices in local coordinate frame
        DO i_v = 1, 3
          vtmp(:,i_v) = MATMUL(Pinv, (vtmp(:,i_v) - Dtmp))
          IF (ABS(vtmp(1,i_v)) .lt.v_min) vtmp(1,i_v) = SIGN(v_min, vtmp(1,i_v))
        END DO

        v_loc(:,:,i_f) = vtmp(:,1:3)
        P(:,:,i_f) = Ptmp
        D(:,i_f) = Dtmp

      END DO
      END SUBROUTINE GET_N_HELPERS
        
      END SUBROUTINE mumaterial_init_demag
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      PURE FUNCTION mumaterial_get_N(P, D, v_tile, pos) result(N)
!!-----------------------------------------------------------------------
!! Determines demagnetization tensor for a tile at position using existing helper arrays
!!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: P(3,3,4) !! Face rotation matrix array of tile
      DOUBLE PRECISION, INTENT(in) :: D(3,4) !! Face base array of tile
      DOUBLE PRECISION, INTENT(in) :: v_tile(3,3,4) !! Vertex array of tile
      DOUBLE PRECISION, INTENT(in) :: pos(3) !! Evaluation position
      DOUBLE PRECISION :: N(3,3) !! Result
      DOUBLE PRECISION, PARAMETER :: d_min = 1.0D-6

      DOUBLE PRECISION :: N_loc(3,3), v_loc(3,3), Pinv(3,3), r_in(3), Ptmp(3,3), larr_in(2)
      INTEGER :: i_f, j

      N = 0.d0
      ! Loop over faces
      DO i_f = 1, 4
        ! Evaluation position in local frame
        Ptmp = P(:,:,i_f)
        Pinv = TRANSPOSE(Ptmp)
        r_in = MATMUL(Pinv, (pos - D(:,i_f)))
        v_loc = v_tile(:,:,i_f)

        ! make sure position is not too close to face origin
        DO j = 1, 3
          IF (ABS(r_in(j)) .lt. d_min)  r_in(j) = SIGN(d_min, r_in(j))
        END DO

        ! Difficult math
        larr_in = (/v_loc(1,1), v_loc(1,3)/)
        N_loc = 0.0d0
        N_loc(:,3) =  GET_NLOC(r_in, larr_in, v_loc(2,2))
        N = N + MATMUL(MATMUL(Ptmp, N_loc), Pinv)
      END DO
      CONTAINS

      PURE FUNCTION GET_NLOC(r, larr, h) RESULT(N_temp)
!!-----------------------------------------------------------------------
!! Helper function to determine the demagnetization tensor.
!! Combines earlier functions of get_box_nxz, get_Nyz, get_Nzz to reduce
!! the number of operations.
!!-----------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION :: N_temp(3)
      DOUBLE PRECISION, INTENT(IN) :: r(3), larr(2), h

      DOUBLE PRECISION :: r1, r2, r3, ir3
      DOUBLE PRECISION :: r1_2, r2_2, r3_2
      DOUBLE PRECISION :: rnorm_2, h_2

      DOUBLE PRECISION :: C1, C2, C3, C4, C5, C6
      DOUBLE PRECISION :: l, l_2
      DOUBLE PRECISION :: sqrt1, sqrt3, sqrt6
      DOUBLE PRECISION :: isqrt1, isqrt3, isqrt6, isqrt13, isqrt16
      DOUBLE PRECISION :: div1, div2

      DOUBLE PRECISION :: F1, F2, K1, K2, F12, K12
      DOUBLE PRECISION :: L1, L2, Q1, Q2
      DOUBLE PRECISION :: G1, G2, P1, P2

      DOUBLE PRECISION :: Nlocx, Nlocy, Nlocz, s
      INTEGER :: i

      ! Fixed geometry
      r1 = r(1)
      r2 = r(2)
      r3 = r(3)
      r1_2 = r1 * r1
      r2_2 = r2 * r2
      r3_2 = r3 * r3
      ir3 = 1.0d0/r3
      rnorm_2 = r1_2 + r2_2 + r3_2

      h_2 = h * h

      C6 = rnorm_2 - 2 * r2 * h + h_2
      sqrt6 = sqrt(C6)
      isqrt6 = 1.0d0 / sqrt6
      ! Different l
      Nlocx = 0.0d0
      Nlocy = 0.0d0
      Nlocz = 0.0d0
      DO i = 1, 2
        s = 3.0d0-2.0d0*i
        l = larr(i)
        l_2 = l * l

        C1 = l_2 + h_2
        C2 = l_2 - l * r1 + h * r2
        C3 = rnorm_2 - 2 * r1 * l + l_2
        C4 = h_2 + l * r1 - h * r2
        C5 = h * (r1_2 + r3_2) / l

        sqrt1 = sqrt(C1)
        sqrt3 = sqrt(C3)

        isqrt1 = 1.0d0 / sqrt1
        isqrt3 = 1.0d0 / sqrt3
        isqrt13 = isqrt1*isqrt3
        isqrt16 = isqrt1*isqrt6
        div1 = h * isqrt1
        div2 = l * isqrt1

        F1 = div1 * ATANH((C2 - C1)*isqrt16)
        F2 = div1 * ATANH(C2*isqrt13)

        K1 = div2 * ATANH((C4 - C1)*isqrt13)
        K2 = div2 * ATANH(C4*isqrt16)

        L1 = ATANH((r1 - l) * isqrt3)

        P1 = ATAN((r1 * (h - r2) - (h * (l - r1) - r2 * l) - C5) * ir3 * isqrt3)
        P2 = ATAN((r1 * (h - r2) - C5) * ir3*isqrt6)

        Q1 = -ATAN((r1 - l) * r2  * ir3 * isqrt3)

        Nlocx = Nlocx - INV4PI*s*(F1 - F2)
        Nlocy = Nlocy - INV4PI*s*(K1 - K2 - L1)
        Nlocz = Nlocz - INV4PI*s*(P1 - P2 - Q1)
      END DO

      N_temp = (/Nlocx, Nlocy, Nlocz/)

      

      END FUNCTION GET_NLOC
      
      END FUNCTION mumaterial_get_N
!-----------------------------------------------------------------------
!--------------------------------------------------------------------

      SUBROUTINE mumaterial_init_tree()
!!----------------------------------------------------------------
!! Sets up everything related to just the BH-tree
!!----------------------------------------------------------------
      IMPLICIT NONE

      INTEGER              :: n_subtrees
      INTEGER, ALLOCATABLE :: subtree_map(:)
      DOUBLE PRECISION, ALLOCATABLE :: subtree_work(:)
      INTEGER, ALLOCATABLE :: subtree_to_rank(:)
      INTEGER, ALLOCATABLE :: nsubtrees_per_rank(:)
      INTEGER, ALLOCATABLE :: subtree_order(:)

      IF (lverb) WRITE(6,*) "  MUMAT_INIT: Building the tree"

      ! Build tree
      CALL mumaterial_build_tree()
      CALL BUILD_SUBTREE_MAP()
      ! Assign work
      CALL ASSIGN_SUBTREES()
      CALL ASSIGN_NODES()
      CALL ASSIGN_LEAVES()
      ! Interactions
      CALL SETUP_INTERACTIONS()
      ! Allocate
      CALL ALLOC_NODE_ARRAYS()
      CALL GET_LEVELS()
      ! Diagnostics
      CALL TREE_DIAGNOSTICS()


      CONTAINS
      SUBROUTINE GET_LEVELS()
      IMPLICIT NONE
      INTEGER :: inode_loc, lev
      
      ALLOCATE(nnodes_per_level_loc(0:tree_depth))
      nnodes_per_level_loc = 0
      DO inode_loc = 1, nnode_loc
        lev = node_level_loc(inode_loc)
        nnodes_per_level_loc(lev) = nnodes_per_level_loc(lev) + 1
      END DO
      
      ALLOCATE(nodes_level_loc(MAXVAL(nnodes_per_level_loc), 0:tree_depth))
      nodes_level_loc = NODE_NOCHILD
      nnodes_per_level_loc = 0 
      DO inode_loc = 1, nnode_loc
        lev = node_level_loc(inode_loc)
        nnodes_per_level_loc(lev) = nnodes_per_level_loc(lev) + 1
        nodes_level_loc(nnodes_per_level_loc(lev), lev) = inode_loc
      END DO
      
      END SUBROUTINE GET_LEVELS

      SUBROUTINE BUILD_SUBTREE_MAP()
!!-----------------------------------------------------------------------
!! Determines split level and builds map from node index to subtree index
!!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: inode_loc, isub, inode
      
      split_level = CEILING(LOG(DBLE(world_size)) / LOG(8.0d0))
      n_subtrees  = 8**split_level
      
      ALLOCATE(subtree_map(nnode))
      subtree_map = -1
      isub = 0
      DO inode = 1, nnode
        IF (node_level(inode) == split_level) THEN
          isub = isub + 1
          subtree_map(inode) = isub
        END IF
      END DO
      
      END SUBROUTINE BUILD_SUBTREE_MAP

      SUBROUTINE ASSIGN_SUBTREES()
!!-----------------------------------------------------------------------
!! Counts work per subtree and assigns subtrees to ranks using LPT
!!-----------------------------------------------------------------------
      USE mpi_inc
      IMPLICIT NONE

      INTEGER :: i, irank, isub, isubo, inode, counter
      DOUBLE PRECISION :: work
      DOUBLE PRECISION :: rank_work_loc
      DOUBLE PRECISION, ALLOCATABLE :: rank_work(:)

      ALLOCATE(subtree_work(n_subtrees))
      subtree_work = 0
      
      DO inode = 1, nnode
        IF (node_level(inode) /= split_level) CYCLE
        IF (MOD(inode-1, world_size) /= world_rank) CYCLE
        isub = subtree_map(inode)
        counter = 0
        work = 0.0d0
        CALL COUNT_SUBTREE_WORK(inode, counter, work)
        subtree_work(isub) = work
      END DO
      
#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_ALLREDUCE(MPI_IN_PLACE, subtree_work, n_subtrees, &
                                      MPI_DOUBLE_PRECISION, MPI_SUM, comm_world, ierr_mpi)
#endif
      
      ! Sort subtrees by decreasing work (insertion sort)
      ALLOCATE(subtree_order(n_subtrees))
      DO isub = 1, n_subtrees
        subtree_order(isub) = isub
      END DO
      DO i = 2, n_subtrees
        DO isub = i, 2, -1
          IF (subtree_work(subtree_order(isub)) > subtree_work(subtree_order(isub-1))) THEN
            isubo = subtree_order(isub)
            subtree_order(isub) = subtree_order(isub-1)
            subtree_order(isub-1) = isubo
          ELSE
            EXIT
          END IF
        END DO
      END DO
      
      ! LPT assignment
      ALLOCATE(subtree_to_rank(n_subtrees))
      ALLOCATE(nsubtrees_per_rank(0:world_size-1))
      ALLOCATE(rank_work(0:world_size-1))
      nsubtrees_per_rank = 0
      rank_work = 0.0d0
      
      DO i = 1, n_subtrees
        isubo = subtree_order(i)
        irank = MINLOC(rank_work, 1) - 1
        subtree_to_rank(isubo) = irank
        nsubtrees_per_rank(irank) = nsubtrees_per_rank(irank) + 1
        rank_work(irank) = rank_work(irank) + subtree_work(isubo)
      END DO
      
      DEALLOCATE(subtree_order, rank_work, subtree_work)
      
      END SUBROUTINE ASSIGN_SUBTREES

      RECURSIVE SUBROUTINE COUNT_SUBTREE_WORK(inode, counter, work)
!!-----------------------------------------------------------------------
!! Calculate interactions workload of a subtree with root inode
!!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in) :: inode
      INTEGER, INTENT(inout) :: counter
      DOUBLE PRECISION, INTENT(inout) :: work
      INTEGER :: c, child, cnt
      DOUBLE PRECISION :: wrk
      
      IF (node_child(1,inode) == NODE_NOCHILD) THEN
        cnt = 0 ! Leaf - count its interactions
        wrk = 0.0d0
        CALL mumaterial_interactions_count(node_cen(:,inode), ROOT, cnt, wrk, .FALSE.)
        counter = counter + cnt
        work = work + wrk
      ELSE
        DO c = 1, 8
          CALL COUNT_SUBTREE_WORK(node_child(c,inode), counter, work)
        END DO
      END IF
      
      END SUBROUTINE COUNT_SUBTREE_WORK

      SUBROUTINE ASSIGN_NODES()
!!-----------------------------------------------------------------------
!! Assigns nodes to ranks based on subtree ownership
!!-----------------------------------------------------------------------
      IMPLICIT NONE 

      INTEGER :: iroot, inode, isub, inode_loc
      
      nnode_loc = 0
      DO inode = 1, nnode
        IF (node_level(inode) < split_level) THEN
          IF (world_rank == 0) nnode_loc = nnode_loc + 1
        ELSE
          iroot = SUBTREE_GET_ROOT(inode, split_level)
          isub  = subtree_map(iroot)
          IF (subtree_to_rank(isub) == world_rank) nnode_loc = nnode_loc + 1
        END IF
      END DO
      
      ALLOCATE(node_local_to_global(nnode_loc), &
                node_level_loc(nnode_loc), &
                node_child_loc(8,nnode_loc))
      
      inode_loc = 0
      DO inode = 1, nnode
        IF (node_level(inode) < split_level) THEN
          IF (world_rank == 0) THEN
            inode_loc = inode_loc + 1
            node_local_to_global(inode_loc) = inode
            node_level_loc(inode_loc) = node_level(inode)
            node_child_loc(:,inode_loc) = node_child(:,inode)
          END IF
        ELSE
          iroot = SUBTREE_GET_ROOT(inode, split_level)
          isub  = subtree_map(iroot)
          IF (subtree_to_rank(isub) == world_rank) THEN
            inode_loc = inode_loc + 1
            node_local_to_global(inode_loc) = inode
            node_level_loc(inode_loc) = node_level(inode)
            node_child_loc(:,inode_loc) = node_child(:,inode)
          END IF
        END IF
      END DO
      
      DEALLOCATE(subtree_map, subtree_to_rank, nsubtrees_per_rank)
      
      END SUBROUTINE ASSIGN_NODES

      PURE FUNCTION SUBTREE_GET_ROOT(inode, level_req) RESULT(iroot)
!!-----------------------------------------------------------------------
!! Finds the index of the parent node that is at the requested level
!!-----------------------------------------------------------------------
      INTEGER, INTENT(in) :: inode !! Index we want to know highest parent of
      INTEGER, INTENT(in) :: level_req !! Level of parent requested
      INTEGER :: iroot !! Parent index
      
      iroot = inode
      DO WHILE (node_level(iroot) > level_req)
        iroot = node_parent(iroot) ! Climb up the tree
      END DO
      
      END FUNCTION SUBTREE_GET_ROOT

      SUBROUTINE ASSIGN_LEAVES()
!!-----------------------------------------------------------------------
!! Assign leaves to MPI ranks so that each rank has a roughly equal load 
!!-----------------------------------------------------------------------
      USE mpi_inc
#if defined(MPI_OPT)
      USE mpi_params, ONLY: MPI_CALC_MYRANGE
#endif
      IMPLICIT NONE

      INTEGER :: ileaf_loc, inode, ileaf, irank, counter, i, itmp
      DOUBLE PRECISION :: wval
      DOUBLE PRECISION, ALLOCATABLE ::  work(:)

      INTEGER, ALLOCATABLE :: leaf_order(:), leaf_to_rank(:)
      DOUBLE PRECISION, ALLOCATABLE :: rank_work(:)
      INTEGER :: leaf_start, leaf_end


      INTEGER :: itile
      !-----------------------------
      ! COUNT INTERACTIONS PER LEAF
      !-----------------------------
      ALLOCATE(work(nleaf))
      work = 0

      CALL MPI_CALC_MYRANGE(comm_world, 1, nleaf, leaf_start, leaf_end)
      DO ileaf = leaf_start, leaf_end
        inode = leaf_to_node(ileaf)
        counter = 0
        wval = 0.0d0
        CALL mumaterial_interactions_count(node_cen(:,inode), ROOT, counter, wval, .FALSE.)
        work(ileaf) = wval
      END DO
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, work, nleaf, MPI_DOUBLE_PRECISION, MPI_SUM, comm_world, ierr_mpi)

      !-----------------------------
      ! SORT LEAVES BY DECREASING WORK
      !-----------------------------
      ALLOCATE(leaf_order(nleaf))
      DO ileaf = 1, nleaf
        leaf_order(ileaf) = ileaf
      END DO

      ! Insertion sort
      DO i = 2, nleaf
        DO ileaf = i, 2, -1
          IF (work(leaf_order(ileaf)) > work(leaf_order(ileaf-1))) THEN
            itmp = leaf_order(ileaf)
            leaf_order(ileaf) = leaf_order(ileaf-1)
            leaf_order(ileaf-1) = itmp
          ELSE
            EXIT
          END IF
        END DO
      END DO

      !-----------------------------
      ! LPT ASSIGNMENT
      !-----------------------------
      ALLOCATE(rank_work(0:world_size-1))
      ALLOCATE(leaf_to_rank(nleaf))
      rank_work = 0.0d0

      DO i = 1, nleaf
        ileaf = leaf_order(i)
        irank = MINLOC(rank_work, 1) - 1
        leaf_to_rank(ileaf) = irank
        rank_work(irank) = rank_work(irank) + work(ileaf)
      END DO
      DEALLOCATE(leaf_order, rank_work, work)

      !-----------------------------
      ! BUILD LOCAL LEAF ARRAYS
      !-----------------------------
      nleaf_loc = COUNT(leaf_to_rank == world_rank)

      ALLOCATE(leaf_local_to_node_global(nleaf_loc), &
              leaf_size_loc(nleaf_loc))

      ileaf_loc = 0
      DO ileaf = 1, nleaf
        IF (leaf_to_rank(ileaf) == world_rank) THEN
          ileaf_loc = ileaf_loc + 1
          inode = leaf_to_node(ileaf)
          leaf_local_to_node_global(ileaf_loc) = inode
          leaf_size_loc(ileaf_loc) = leaf_size(inode)
        END IF
      END DO
      DEALLOCATE(leaf_to_rank)

      IF (nleaf_loc > 0) THEN
        max_leaf_size_seen = MAXVAL(leaf_size_loc)
      ELSE
        max_leaf_size_seen = 0
      END IF

      CALL MPI_ALLREDUCE(max_leaf_size_seen, max_leaf_size, 1, MPI_INTEGER, MPI_MAX, comm_world, ierr_mpi)
      ALLOCATE(leaf_tile_loc(max_leaf_size_seen, nleaf_loc))
      DO ileaf_loc = 1, nleaf_loc
        inode = leaf_local_to_node_global(ileaf_loc)
        leaf_tile_loc(:,ileaf_loc) = leaf_tile(1:max_leaf_size_seen, inode)
      END DO
      
      END SUBROUTINE ASSIGN_LEAVES

      SUBROUTINE SETUP_INTERACTIONS()
!!-----------------------------------------------------------------------
!! Builds interaction lists for local leaves
!!-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: ileaf_loc, inode, counter, ninteractions
      INTEGER, ALLOCATABLE :: work_loc(:)
      INTEGER, ALLOCATABLE :: cursor(:)
      DOUBLE PRECISION :: wval
    
      IF (lverb) WRITE(6,*) "  MUMAT_INIT: Filling tree interactions"
      ALLOCATE(work_loc(nleaf_loc))
      ALLOCATE(interact_ptr(nleaf_loc+1), cursor(nleaf_loc))
      
      DO ileaf_loc = 1, nleaf_loc
        inode = leaf_local_to_node_global(ileaf_loc)
        counter = 0
        wval = 0.0d0
        CALL mumaterial_interactions_count(node_cen(:,inode), ROOT, counter, wval, .FALSE.)
        work_loc(ileaf_loc) = counter
      END DO
      
      interact_ptr(1) = 1
      DO ileaf_loc = 1, nleaf_loc
        interact_ptr(ileaf_loc+1) = interact_ptr(ileaf_loc) + work_loc(ileaf_loc)
      END DO
      ninteractions = interact_ptr(nleaf_loc+1) - 1
      DEALLOCATE(work_loc)
      
      ALLOCATE(interact_list(ninteractions), interact_type(ninteractions))
      cursor(:) = interact_ptr(1:nleaf_loc)
      DO ileaf_loc = 1, nleaf_loc
        inode = leaf_local_to_node_global(ileaf_loc)
        CALL mumaterial_interactions_fill(node_cen(:,inode), ROOT, ileaf_loc, cursor,.FALSE.)
      END DO
      DEALLOCATE(cursor)

      END SUBROUTINE SETUP_INTERACTIONS

      SUBROUTINE ALLOC_NODE_ARRAYS()
!!-----------------------------------------------------------------------
!! Allocates shared memory arrays for node quantities
!!-----------------------------------------------------------------------
      INTEGER :: ileaf_loc, offset
      
      IF (lcomm) THEN
#if defined(MPI_OPT)
        CALL mpialloc(node_m,   3, nnode, shar_rank, 0, comm_shar, win_node_m)
        CALL mpialloc(node_p,   3, nnode, shar_rank, 0, comm_shar, win_node_p)
        CALL mpialloc(node_com, 3, nnode, shar_rank, 0, comm_shar, win_node_com)
        CALL mpialloc(node_w,      nnode, shar_rank, 0, comm_shar, win_node_w)
        CALL mpialloc(node_q,   5, nnode, shar_rank, 0, comm_shar, win_node_q)
#endif
      ELSE
        ALLOCATE(node_m(3,nnode), node_p(3,nnode), node_com(3,nnode), &
                  node_w(nnode), node_q(5,nnode))
      END IF
      
      IF (shar_rank==0) THEN
        node_m   = 0.0d0
        node_p   = 0.0d0
        node_com = 0.0d0
        node_w   = 0.0d0
        node_q   = 0.0d0
      END IF

#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_BARRIER(comm_shar, ierr_mpi)
#endif
      
      ntile_loc = SUM(leaf_size_loc(1:nleaf_loc))
      ALLOCATE(leaf_offset_loc(nleaf_loc))
      offset = 0
      DO ileaf_loc = 1, nleaf_loc
        leaf_offset_loc(ileaf_loc) = offset
        offset = offset + leaf_size_loc(ileaf_loc)
      END DO
      END SUBROUTINE ALLOC_NODE_ARRAYS

      SUBROUTINE TREE_DIAGNOSTICS
      USE mpi_inc
      IMPLICIT NONE

      INTEGER :: inode, ileaf, csr
      INTEGER :: maxleaf, sumleaf
      INTEGER :: num_structural_leaves
      INTEGER :: lev, n_lev, l_lev, t_lev
      INTEGER :: n_node_int, n_leaf_int, n_tile_int
      DOUBLE PRECISION :: meanleaf, avg_node_int, avg_leaf_int, avg_tile_int

      ! structural node statistics
      num_structural_leaves = 0
      DO inode = 1, nnode
        IF (node_child(1,inode) == NODE_NOCHILD) &
          num_structural_leaves = num_structural_leaves + 1
      END DO

      ! active leaf statistics
      maxleaf = 0
      sumleaf = 0
      DO ileaf = 1, nleaf
        inode = leaf_to_node(ileaf)
        sumleaf = sumleaf + leaf_size(inode)
        maxleaf = MAX(maxleaf, leaf_size(inode))
      END DO
      meanleaf = MERGE(DBLE(sumleaf)/DBLE(nleaf), 0.0d0, nleaf>0)

      ! interaction statistics
      n_node_int = 0; n_leaf_int = 0; n_tile_int = 0
      DO ileaf = 1, nleaf_loc
        DO csr = interact_ptr(ileaf), interact_ptr(ileaf+1)-1
          IF (interact_type(csr) == INTERACT_NODE) THEN
            n_node_int = n_node_int + 1
          ELSE IF (interact_type(csr) == INTERACT_LEAF) THEN
            n_leaf_int = n_leaf_int + 1
            n_tile_int = n_tile_int + leaf_size(interact_list(csr))
          END IF
        END DO
      END DO
#if defined(MPI_OPT)
      IF (lcomm) THEN
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,n_node_int,1,MPI_INTEGER,MPI_SUM,comm_world,ierr_mpi)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,n_leaf_int,1,MPI_INTEGER,MPI_SUM,comm_world,ierr_mpi)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,n_tile_int,1,MPI_INTEGER,MPI_SUM,comm_world,ierr_mpi)
      END IF
#endif
      avg_node_int = DBLE(n_node_int) / DBLE(nleaf)
      avg_leaf_int = DBLE(n_leaf_int) / DBLE(nleaf)
      avg_tile_int = DBLE(n_tile_int) / DBLE(nleaf)
      IF (lverb) THEN
        WRITE(6,'(A)') "================ TREE SUMMARY ================"
        WRITE(6,'(A,I10)')  "  Split level        : ", split_level
        WRITE(6,'(A,F10.2)')"  Theta              : ", tree_theta
        WRITE(6,'(A,I10)')  "  Total nodes        : ", nnode
        WRITE(6,'(A,I10)')  "  Max depth          : ", tree_depth
        WRITE(6,'(A,I10)')  "  Total leaves       : ", num_structural_leaves
        WRITE(6,'(A,I10)')  "  Non-empty leaves   : ", nleaf
        WRITE(6,'(A,I10)')  "  Total tiles        : ", sumleaf
        WRITE(6,'(A,I10)')  "  Max tiles/leaf     : ", maxleaf
        WRITE(6,'(A,F10.2)')"  Avg tiles/leaf     : ", meanleaf
        WRITE(6,'(A)') "-----------------------------------------------"
        WRITE(6,'(A)') "  Level     Nodes   Leavesact   Tiles"
        WRITE(6,'(A)') "-----------------------------------------------"
        DO lev = 0, tree_depth
          n_lev = 0; l_lev = 0; t_lev = 0
          DO inode = 1, nnode
            IF (node_level(inode) == lev) THEN
              n_lev = n_lev + 1
              IF (node_child(1,inode)==NODE_NOCHILD .AND. leaf_size(inode)>0) THEN
                l_lev = l_lev + 1
                t_lev = t_lev + leaf_size(inode)
              END IF
            END IF
          END DO
          IF (n_lev > 0) WRITE(6,'(2X,I5,3(2X,I8))') lev, n_lev, l_lev, t_lev
        END DO
        WRITE(6,'(A)') "-----------------------------------------------"
        WRITE(6,'(A)') "  Interaction Statistics (per active leaf)"
        WRITE(6,'(A)') "-----------------------------------------------"
        WRITE(6,'(A,F10.2)') "  Avg aggregate nodes : ", avg_node_int
        WRITE(6,'(A,F10.2)') "  Avg exact leaves    : ", avg_leaf_int
        WRITE(6,'(A,F10.2)') "  Avg exact elements  : ", avg_tile_int
        WRITE(6,'(A)') "==============================================="
      END IF
      
      END SUBROUTINE TREE_DIAGNOSTICS
      END SUBROUTINE mumaterial_init_tree
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE mumaterial_build_tree()
!!-----------------------------------------------------------------------
!! One rank constructs the tree and shares tree info with other ranks.
!!-----------------------------------------------------------------------
      USE mpi_inc

      IMPLICIT NONE
      DOUBLE PRECISION :: bounds(6)
      INTEGER :: ii, ileaf, i_node, nnode_max, max_size, i_tile, jtile, ldex
      INTEGER, PARAMETER :: NO_PARENT = 0, ZERO_DEPTH = 0
      INTEGER, ALLOCATABLE :: leafdex(:), ltn(:), ls(:), lt(:,:)

      DOUBLE PRECISION, ALLOCATABLE :: node_bounds_build(:,:) !! (6,10*ntet) [xmin xmax ymin ymax zmin zmax] of a node
      DOUBLE PRECISION, ALLOCATABLE :: node_cen_build(:,:) !! (3,10*ntet) geometric center of a node
      INTEGER, ALLOCATABLE :: leaf_head_build(:) !! (10*ntet) Global element index that is header of leaf linked list.
      INTEGER, ALLOCATABLE :: tile_next_build(:) !! (ntet) Global element index that is next in leaf linked list (-1=end)
      INTEGER, ALLOCATABLE :: node_child_build(:,:) !! (8,10*ntet) Node indices of child octants of a node
      INTEGER, ALLOCATABLE :: leaf_size_build(:) !! (10*ntet) Number of elements in a node (0 if internal node)
      INTEGER, ALLOCATABLE :: node_level_build(:) !! (10*ntet) Level of node (root=0)
      INTEGER, ALLOCATABLE :: node_parent_build(:) !! (10*ntet) Parent of node

      IF (world_rank.EQ.0) THEN
        ! Initialize ; allocate larger-than-necessary arrays first
        nnode = 0
        nnode_max = 10*ntet

        ALLOCATE(node_bounds_build(6,nnode_max),   &
                node_cen_build(3,nnode_max), &
                node_child_build(8,nnode_max), &
                leaf_size_build(nnode_max), &
                node_level_build(nnode_max), &
                node_parent_build(nnode_max))
        ALLOCATE(tile_next_build(ntet), leaf_head_build(nnode_max))

        node_bounds_build   = 0.0d0         ! No bounds
        node_cen_build      = 0.0d0         ! No midpoints
        node_child_build    = NODE_NOCHILD  ! No children
        leaf_size_build     = 0             ! No elements
        leaf_head_build     = LEAF_NOHEAD   ! No header element
        tile_next_build     = TILE_NONEXT   ! Elements not linked
        node_level_build    = NODE_BADPARENT! Bad parent

        ! Get bounds
        DO ii = 1,3
          bounds(2*(ii-1)+1) = MINVAL(tet_cen(ii,:))
          bounds(2*ii)       = MAXVAL(tet_cen(ii,:))
        END DO

        ! Create root
        CALL NODE_CREATE(bounds, NO_PARENT, i_node, ZERO_DEPTH)

        ! Insert all points (builds the tree)
        DO ii = 1, ntet
          CALL NODE_INSERT(i_node, ii, ZERO_DEPTH)
        END DO

        tree_depth = MAXVAL(node_level_build)
        ! Get (non-empty) leaves
        ALLOCATE(leafdex(nnode))
        leafdex = NODE_NOTLEAF
        nleaf = 0
        DO i_node = 1, nnode
          IF (node_child_build(1,i_node)==NODE_NOCHILD.AND.leaf_size_build(i_node)>0) THEN
            nleaf = nleaf+1
            leafdex(i_node) = nleaf
          END IF
        END DO

        ! More useful: leaf_to_node map
        ALLOCATE(ltn(nleaf))
        DO i_node = 1, nnode
          IF (leafdex(i_node) /= NODE_NOTLEAF) THEN
            ltn(leafdex(i_node)) = i_node
          END IF
        END DO
        DEALLOCATE(leafdex)

        ! More useful: leaf_size array of size nleaf
        ALLOCATE(ls(nnode)) ! Leaf size
        ls = 0
        DO ileaf = 1, nleaf
          i_node = ltn(ileaf)
          ls(i_node) = leaf_size_build(i_node)
        END DO
        DEALLOCATE(leaf_size_build)

        ! Flatten linked lists
        max_size = MAXVAL(ls(:))
        ALLOCATE(lt(max_size, nnode))
        lt = LEAF_NOTTILE
        DO ileaf = 1, nleaf
          i_node = ltn(ileaf)
          i_tile = leaf_head_build(i_node)
          DO ii = 1, ls(i_node)
            lt(ii, i_node) = i_tile
            i_tile = tile_next_build(i_tile)
          END DO
        END DO
        DEALLOCATE(leaf_head_build,tile_next_build)
      END IF

      ! Proper allocation; share the tree
      IF (lcomm) THEN
#if defined(MPI_OPT)
        IF (shar_rank.EQ.0) THEN
          CALL MPI_Bcast(nnode,1,MPI_INTEGER,0,comm_master,ierr_mpi)
          CALL MPI_Bcast(nleaf,1,MPI_INTEGER,0,comm_master,ierr_mpi)
          CALL MPI_Bcast(max_size,1,MPI_INTEGER,0,comm_master,ierr_mpi)
          CALL MPI_Bcast(tree_depth,1,MPI_INTEGER,0,comm_master,ierr_mpi)
        END IF
        CALL MPI_Bcast(nnode,1,MPI_INTEGER,0,comm_shar,ierr_mpi)
        CALL MPI_Bcast(nleaf,1,MPI_INTEGER,0,comm_shar,ierr_mpi)
        CALL MPI_Bcast(max_size,1,MPI_INTEGER,0,comm_shar,ierr_mpi)
        CALL MPI_Bcast(tree_depth,1,MPI_INTEGER,0,comm_shar,ierr_mpi)
        ! allocate
        CALL mpialloc(node_bounds,6,nnode,          shar_rank,0,comm_shar,win_node_bounds)
        CALL mpialloc(node_cen,   3,nnode,          shar_rank,0,comm_shar,win_node_cen)
        CALL mpialloc(node_child, 8,nnode,          shar_rank,0,comm_shar,win_node_child)
        CALL mpialloc(node_level,   nnode,          shar_rank,0,comm_shar,win_node_level)
        CALL mpialloc(node_parent,  nnode,          shar_rank,0,comm_shar,win_node_parent)
        CALL mpialloc(leaf_size,    nnode,          shar_rank,0,comm_shar,win_leaf_size)
        CALL mpialloc(leaf_to_node, nleaf,          shar_rank,0,comm_shar,win_leaf_to_node)
        CALL mpialloc(leaf_tile,max_size,    nnode, shar_rank,0,comm_shar,win_leaf_tile)
#endif
      ELSE
        ALLOCATE(node_bounds(6,nnode), node_cen(3,nnode), &
                  node_child(8,nnode),node_level(nnode), node_parent(nnode), &
                  leaf_size(nnode), leaf_to_node(nleaf), leaf_tile(max_size,nnode))
      END IF

      IF (world_rank.EQ.0) THEN
        node_bounds = node_bounds_build(:,1:nnode)
        node_cen = node_cen_build(:,1:nnode)
        node_child = node_child_build(:,1:nnode)
        node_parent = node_parent_build(1:nnode)
        node_level = node_level_build(1:nnode)
        leaf_size = ls
        leaf_to_node = ltn
        leaf_tile = lt
        DEALLOCATE(node_bounds_build, node_cen_build, &
                    node_child_build, node_level_build, node_parent_build, &
                  ls, ltn, lt)
      ENDIF

      IF (lcomm) THEN
#if defined(MPI_OPT)
        IF (shar_rank.EQ.0) THEN ! Share with masters, who write in the window
          CALL MPI_Bcast(node_bounds,6*nnode,MPI_DOUBLE_PRECISION,0,comm_master,ierr_mpi)
          CALL MPI_Bcast(node_cen,   3*nnode,MPI_DOUBLE_PRECISION,0,comm_master,ierr_mpi)
          CALL MPI_Bcast(node_child, 8*nnode,MPI_INTEGER,0,comm_master,ierr_mpi)
          CALL MPI_Bcast(node_parent,  nnode,MPI_INTEGER,0,comm_master,ierr_mpi)
          CALL MPI_Bcast(node_level,   nnode,MPI_INTEGER,0,comm_master,ierr_mpi)
          CALL MPI_Bcast(leaf_size,    nnode,MPI_INTEGER,0,comm_master,ierr_mpi)
          CALL MPI_Bcast(leaf_tile,max_size*nnode,MPI_INTEGER,0,comm_master,ierr_mpi)
          CALL MPI_Bcast(leaf_to_node, nleaf,MPI_INTEGER,0,comm_master,ierr_mpi)
        END IF
        CALL MPI_Barrier(comm_shar, ierr_mpi)
#endif
      END IF
      
      CONTAINS
  
      SUBROUTINE NODE_CREATE(bounds_in, iparent, inode, level)
!!-----------------------------------------------------------------------
!! Creates a new node, appending info to the end of the arrays.
!!-----------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in) :: bounds_in(6) !! [xmin, xmax, ymin, ymax, zmin, zmax]
      INTEGER, INTENT(in)  :: iparent !! Parent index of node
      INTEGER, INTENT(in)  :: level !! Level of new node
      INTEGER, INTENT(out) :: inode !! new index

      nnode = nnode+1
      inode = nnode

      node_bounds_build(:,inode) = bounds_in
      node_cen_build(1,inode) = (bounds_in(1)+bounds_in(2))/2.0d0
      node_cen_build(2,inode) = (bounds_in(3)+bounds_in(4))/2.0d0
      node_cen_build(3,inode) = (bounds_in(5)+bounds_in(6))/2.0d0
      node_child_build(:,inode) = NODE_NOCHILD
      node_level_build(inode) = level
      node_parent_build(inode) = iparent

      END SUBROUTINE NODE_CREATE

      SUBROUTINE NODE_SPLIT(inode, depth)
!!-----------------------------------------------------------------------
!! Splits the given node into octants. Notes on indexing:
!!  c   binary   description
!!  0    000     x <  xmid, y <  ymid, z <  zmid
!!  1    001     x >= xmid, y <  ymid, z <  zmid
!!  2    010     x <  xmid, y >= ymid, z <  zmid
!!  3    011     x >= xmid, y >= ymid, z <  zmid
!!  4    100     x <  xmid, y <  ymid, z >= zmid
!!  5    101     x >= xmid, y <  ymid, z >= zmid
!!  6    110     x <  xmid, y >= ymid, z >= zmid
!!  7    111     x >= xmid, y >= ymid, z >= zmid
!!-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(in) :: inode !! Node to split
      INTEGER, INTENT(in) :: depth !! Depth to create new nodes at
      DOUBLE PRECISION :: bnew(6), b(6), mid(3)
      INTEGER :: c,d,bit,cnode
      INTEGER, PARAMETER :: nchild = 8

      ! Get bounds of inode
      b = node_bounds_build(:,inode)
      mid = node_cen_build(:,inode)

      ! Get bounds of children
      node_child_build(:,inode) = NODE_NOCHILD
      DO c = 0,nchild-1
        bnew = 0.0d0
        DO d = 1, 3
          bit = mod(c/(2**(d-1)),2) ! Bit shift
          IF (bit==0) THEN ! Lower half
            bnew(2*d-1) = b(2*d-1)
            bnew(2*d) = mid(d)
          ELSE             ! Upper half
            bnew(2*d-1) = mid(d)
            bnew(2*d) = b(2*d)
          END IF
        END DO
        CALL NODE_CREATE(bnew, inode, cnode, depth) ! Create node
        node_child_build(c+1,inode) = cnode     ! Add to tree

      END DO

      END SUBROUTINE NODE_SPLIT

      SUBROUTINE NODE_GETCHILD(inode,itile,cnode)
!-----------------------------------------------------------------------
!! Find node child that should accept itile
!!-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(in) :: inode !! Node index
      INTEGER, INTENT(in) :: itile !! Index of element in tet_cen
      INTEGER, INTENT(out) :: cnode ! Child node index
      INTEGER :: c
      DOUBLE PRECISION :: mid(3), pos(3)

      pos = tet_cen(:,itile)
      mid = node_cen_build(:,inode)
      c = 0                           !  000 to start
      IF (pos(1).GE.mid(1)) c = c + 1 ! +001 if x >= xmid
      IF (pos(2).GE.mid(2)) c = c + 2 ! +010 if y >= ymid
      IF (pos(3).GE.mid(3)) c = c + 4 ! +100 if z >= zmid
      cnode = node_child_build(c+1, inode)

      
      END SUBROUTINE NODE_GETCHILD

      RECURSIVE SUBROUTINE NODE_INSERT(inode, itile, depth)
!!-----------------------------------------------------------------------
!! Inserts an element into a node. If the node has no
!! child nodes and new #elements <= max, point is inserted. If it has
!! new # elements > max, node is split into children. Otherwise, point
!! is inserted into the appropriate child node.
!!
!! Elements in a node are stored as a linked list, with new elements
!! prepended to the list.
!!-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(in) :: inode !! Node index
      INTEGER, INTENT(in) :: itile !! Index of element in tet_cen
      INTEGER, INTENT(in) :: depth !! Current depth of node
      INTEGER :: head_tile, cnode, ntmp, i, ktile, jtmp
      INTEGER, ALLOCATABLE :: tmp(:)

      head_tile = leaf_head_build(inode) ! Linked list of old elements

      ! Depth limit
      IF (depth >= tree_max_depth) THEN
        tile_next_build(itile) = head_tile
        leaf_head_build(inode) = itile
        leaf_size_build(inode) = leaf_size_build(inode)+1
        WRITE(6,'(A,I0,9F12.6)') "Max depth reached, inserting in node ", inode, tet_cen(:,itile), node_bounds_build(:,inode)
        RETURN
      END IF

      ! Leaf case (no children)
      IF (node_child_build(1,inode)==NODE_NOCHILD) THEN
        IF (leaf_size_build(inode)<leaf_max_size) THEN ! Not yet full
          tile_next_build(itile) = head_tile
          leaf_head_build(inode) = itile
          leaf_size_build(inode) = leaf_size_build(inode)+1
        ELSE ! Full leaf
          ! Format linked list of indices as array instead
          ntmp = leaf_size_build(inode)+1
          ALLOCATE(tmp(ntmp))
          tmp(1) = itile ! New element
          ktile = leaf_head_build(inode) ! First element in leaf
          DO i = 2, ntmp
            tmp(i) = ktile
            jtmp = tile_next_build(ktile) ! Next element in leaf (temp)
            tile_next_build(ktile) = TILE_NONEXT ! Clear link
            ktile = jtmp            ! Now actually handle next element
          END DO
          ! Remove any elements and split
          leaf_head_build(inode) = LEAF_NOHEAD
          leaf_size_build(inode) = 0
          CALL NODE_SPLIT(inode, depth+1)
          ! Reinsert
          DO i = 1, ntmp
            ktile = tmp(i)
            CALL NODE_GETCHILD(inode, ktile, cnode)
            CALL NODE_INSERT(cnode, ktile, depth+1)
          END DO
          DEALLOCATE(tmp)
        END IF
      ELSE ! Node with children
        CALL NODE_GETCHILD(inode, itile, cnode)
        CALL NODE_INSERT(cnode, itile, depth+1)
      END IF

      END SUBROUTINE NODE_INSERT

      END SUBROUTINE mumaterial_build_tree
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE mumaterial_init_Happ(linitM)
!!-----------------------------------------------------------------------------
!! Initializes static background field at elements using external Bfld function
!!-----------------------------------------------------------------------------
      IMPLICIT NONE

      LOGICAL, INTENT(in) :: linitM !! True to initialize M from Happ
      INTEGER :: ileaf_loc, ntile, itile, leaf(max_leaf_size_seen)

      DOUBLE PRECISION :: pos(3), Bx, By, Bz

      INTEGER :: sdex, tile
      DOUBLE PRECISION :: chi

      IF (lverb) WRITE (6,*) "  MUMAT_INIT:  Calculating H_app"
      ALLOCATE(Happ_loc(3, max_leaf_size_seen, nleaf_loc))
      Happ_loc(:,:,:) = 0.0d0

      DO ileaf_loc = 1, nleaf_loc ! Loop over leaves
        ntile = leaf_size_loc(ileaf_loc)
        leaf = leaf_tile_loc(:,ileaf_loc)
        DO itile = 1, ntile ! Loop over tiles
          pos = tet_cen(:,leaf(itile))
          CALL ext_Bfld(pos(1), pos(2), pos(3), Bx, By, Bz)
          Happ_loc(:,itile,ileaf_loc) = invmu0*(/Bx, By, Bz/)
          ! Initialize M if necessary
          IF (linitM) THEN
            tile = leaf(itile)
            sdex = state_dex(tile)
            CALL mumaterial_get_M(sdex, Happ_loc(:,itile,ileaf_loc), M(:,tile), chi)
          END IF
        END DO
      END DO

#if defined(MPI_OPT)
      IF (ldosync .AND. linitM) CALL mumaterial_sync_M()
#endif 

      END SUBROUTINE mumaterial_init_Happ
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE mumaterial_iterate_M
!!-----------------------------------------------------------------------------
!! Initializes static background field at elements using external Bfld function
!!-----------------------------------------------------------------------------
      USE mpi_inc
      IMPLICIT NONE
      INTEGER :: iter, ileaf_loc, itile_leaf, itile_loc
      INTEGER :: ntiles_leaf, tile
      DOUBLE PRECISION :: r_tiles(3,max_leaf_size_seen), H_tiles(3,max_leaf_size_seen)

      DOUBLE PRECISION :: lambda(ntile_loc), lambda_old, lambda_new
      DOUBLE PRECISION :: M_old(3), M_targ(3), M_new(3), Mnorm, H_targ(3)
      DOUBLE PRECISION :: res(3), res_prev(3,ntile_loc), rM, rM_rel, rMprev_rel(ntile_loc)
      DOUBLE PRECISION :: Vconv, Pconv, Vres, rM_rel_avg
      DOUBLE PRECISION :: rM_bad, pair_in(2), pair_out(2), info_bad(4)
      INTEGER ::  leaf_temp(max_leaf_size_seen) 
      INTEGER :: tile_bad, itile_bad
      INTEGER :: csr_N
      INTEGER :: lev
      DOUBLE PRECISION :: t0, t1

      INTEGER, PARAMETER :: nt = 5
      DOUBLE PRECISION :: t_max(nt), t_min(nt), t_avg(nt), tbuf(nt)

      LOGICAL :: ldone(ntile_loc), lalldone

      lambda = lambdaStart
      rMprev_rel = HUGE(1.0d0); res_prev = HUGE(1.0d0)
      ldone = .FALSE.
      !--------------------------
      ! START OF CONVERGENCE LOOP
      !---------------------------
      DO iter = 1, maxiter
        rM_bad = 0.0d0; 
        Vconv = 0.0d0; Pconv = 0.0d0; Vres = 0.0d0
        tile_bad = 1
        
        !-----------------------------
        ! LOCAL ELEMENT LOOP
        !-----------------------------
        csr_N = 0                ! Reset cursor
        DO ileaf_loc = 1, nleaf_loc ! Loop over local leaves
          ntiles_leaf = leaf_size_loc(ileaf_loc)
          leaf_temp = leaf_tile_loc(:,ileaf_loc)

          ! ! Skip leaf if all done
          ! IF (ALL(ldone(leaf_offset_loc(ileaf_loc)+1 : &
          !       leaf_offset_loc(ileaf_loc)+ntiles_leaf))) THEN
          !   csr_N = csr_N + nN_per_leaf(ileaf_loc) ! Move cursor
          !   DO itile_leaf = 1, ntiles_leaf
          !     tile = leaf_temp(itile_leaf)
          !     itile_loc = leaf_offset_loc(ileaf_loc) + itile_leaf
          !     Vres = Vres + tet_vol(tile)*rMprev_rel(itile_loc)
          !     Vconv = Vconv + tet_vol(tile)
          !   END DO
          !   CYCLE
          ! END IF
          ! Get field at each element
          DO itile_leaf = 1, ntiles_leaf
            tile = leaf_temp(itile_leaf) ! global index
            r_tiles(:,itile_leaf) = tet_cen(:,tile) ! Position
            H_tiles(:,itile_leaf) = Happ_loc(:,itile_leaf,ileaf_loc) ! Reset background field
          END DO
          ! Add contributions from everything
          t0 = MPI_WTIME()
          CALL mumaterial_eval_field_iter(r_tiles, ntiles_leaf, ileaf_loc, H_tiles, csr_N)
          t1 = MPI_WTIME(); t_field = t_field + (t1-t0)

          t0 = t1
          DO itile_leaf = 1, ntiles_leaf
            ! Calculate self-field of each element
            itile_loc = leaf_offset_loc(ileaf_loc) + itile_leaf
            tile = leaf_temp(itile_leaf) ! global index
            CALL mumaterial_solve_MH_NR(H_tiles(:,itile_leaf), Nself_loc(itile_loc,:,:), tile, M_targ, H_targ)

            ! Update magnetization
            M_old = M(:,tile)
            res = M_targ - M_old
            rM = SQRT(res(1)*res(1)+res(2)*res(2)+res(3)*res(3))
            M_new = M_old + lambda(itile_loc)*res
            M(:,tile) = M_new

            ! Residual 
            Mnorm = SQRT(M_targ(1)*M_targ(1)+M_targ(2)*M_targ(2)+M_targ(3)*M_targ(3))
            IF (Mnorm>1.0d-12) THEN
              rM_rel = rM/Mnorm
            ELSE
              rM_rel = 0.0d0
            END IF
            IF (rM_rel > rM_bad) THEN
              rM_bad = rM_rel
              tile_bad = tile
              itile_bad = itile_loc
            END IF
            
            ! Convergence?
            Vres = Vres + tet_vol(tile)*rM_rel
            IF ((rM_rel.LT.eps_max*lambda(itile_loc)).AND.(iter.GT.1)) THEN
              ldone(itile_loc) = .TRUE.
              Vconv = Vconv + tet_vol(tile)
            END IF

            ! Update lambda
            lambda_old = lambda(itile_loc)
            IF (iter>1) THEN
              IF (DOT_PRODUCT(res,res_prev(:,itile_loc))<0.0d0 &
                  .AND.(rM_rel>rMprev_rel(itile_loc))) THEN
                lambda_new = lambda_old*lambdaFactor
              ELSE IF (rM_rel<rMprev_rel(itile_loc)) THEN
                lambda_new = lambda_old*1.05
              ELSE
                lambda_new = lambda_old
              END IF
            ELSE
              lambda_new = lambda_old
            END IF
            lambda(itile_loc) = MAX(MIN(lambda_new, lambdamax),lambdaMin)

            ! Update last residuals
            res_prev(:,itile_loc) = res
            rMprev_rel(itile_loc) = rM_rel
          END DO
          t1 = MPI_WTIME(); t_solve = t_solve + (t1-t0)
        END DO
        !-----------------------------
        ! COMMUNICATION
        !-----------------------------
        info_bad = (/rM_bad, DBLE(tile_bad), lambda(itile_bad), NORM2(M(:,tile_bad))/)
        pair_in(1) = rM_bad; pair_in(2) = DBLE(world_rank)
        IF (lcomm) THEN
#if defined(MPI_OPT)
          ! Worst element
          CALL MPI_ALLREDUCE(pair_in,pair_out,1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, comm_world, ierr_mpi)
          CALL MPI_BCAST(info_bad, 4, MPI_DOUBLE_PRECISION, INT(pair_out(2)), comm_world, ierr_mpi)
          ! Convergence
          pair_in(1) = Vconv; pair_in(2) = Vres
          CALL MPI_ALLREDUCE(MPI_IN_PLACE, pair_in, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_world, ierr_mpi)
#endif
        END IF
        Pconv = pair_in(1)/tet_vol_tot*100.0d0
        rM_rel_avg = pair_in(2)/tet_vol_tot
        IF (lverb) THEN
          IF (iter.EQ.1) THEN
            WRITE(6,'(/,A)') ' Count  %Done  AvgRes      Index    Mnorm       MaxRes      TargRes     Lambda'
            WRITE(6,'(A)')   '================================================================================='
          END IF
          WRITE(6,'(1X,I6,1X,F6.1,1X,ES10.3,1X,I8,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3)') &
            iter, Pconv, rM_rel_avg, INT(info_bad(2)),  info_bad(4), info_bad(1), eps_max*info_bad(3), info_bad(3)
          FLUSH(6)
        END IF

        lalldone = (Pconv>=Pconv_min)

#if defined(MPI_OPT)
        t0 = MPI_WTIME()
        IF (ldosync) CALL mumaterial_sync_M_sparse() ! Synchronize M
        t1 = MPI_WTIME(); t_sync_M = t_sync_M + (t1-t0)
#endif
        CALL mumaterial_propagate_nodes() ! Update node quantities

        IF (lalldone) THEN
          IF (lverb) WRITE(6,*) "  MUMAT:  Stopping"
          EXIT
        END IF

      END DO
#if defined(MPI_OPT)
      IF (ldosync) THEN ! Final M & node sync after iterations finish
        CALL mumaterial_sync_M()
        DO lev = tree_depth, 0, -1
          CALL mumaterial_sync_nodes(lev)
        END DO
      END IF
#endif

      ! Output print timing
      tbuf = (/t_field, t_solve, t_propagate, t_sync_M, t_sync_nodes/)
      IF (lcomm) THEN
#if defined(MPI_OPT)
        CALL MPI_REDUCE(tbuf, t_max, nt, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm_world, ierr_mpi)
        CALL MPI_REDUCE(tbuf, t_min, nt, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm_world, ierr_mpi)
        CALL MPI_REDUCE(tbuf, t_avg, nt, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm_world, ierr_mpi)
        t_avg = t_avg / world_size
#endif
      ELSE
        t_max = tbuf; t_min = tbuf; t_avg = tbuf
      END IF

      IF (lverb) THEN
        WRITE(6,'(/,A)') '  MUMAT Timing Summary:'
        WRITE(6,'(A)') '  Section          Min(s)    Avg(s)    Max(s)'
        WRITE(6,'(A)') '  ============================================'
        WRITE(6,'(A,3F10.2)') '  Field eval:  ', t_min(1), t_avg(1), t_max(1)
        WRITE(6,'(A,3F10.2)') '  MH solve:    ', t_min(2), t_avg(2), t_max(2)
        WRITE(6,'(A,3F10.2)') '  Propagate:   ', t_min(3), t_avg(3), t_max(3)
        WRITE(6,'(A,3F10.2)') '  Sync M:      ', t_min(4), t_avg(4), t_max(4)
        WRITE(6,'(A,3F10.2)') '  Sync nodes:  ', t_min(5), t_avg(5), t_max(5)
      END IF
      

      END SUBROUTINE mumaterial_iterate_M

      SUBROUTINE mumaterial_solve_MH_NR(H_app, N_self, i_tile, M_out, H_out)
!!-------------------------------------------------------------------
!! Uses a Newton-Rhapson approach to solve for the self-consistent M
!! and H of a single element
!!-------------------------------------------------------------------
      USE mpi_inc
      IMPLICIT NONE
      
      DOUBLE PRECISION, INTENT(in)  :: H_app(3) !! Background field (no self-field!)
      DOUBLE PRECISION, INTENT(in)  :: N_self(3,3) !! Self-demag tensor
      INTEGER,          INTENT(in)  :: i_tile !! Index of our tile
      DOUBLE PRECISION, INTENT(out) :: M_out(3) !! Self-consistent M
      DOUBLE PRECISION, INTENT(out) :: H_out(3) !! Self-consistent H
      INTEGER :: iter
      
      INTEGER :: stype, sdex
      DOUBLE PRECISION :: A(3,3), JAC(3,3), dMdH_mat(3,3)
      DOUBLE PRECISION :: F(3), dH(3), H_self(3)
      DOUBLE PRECISION :: chi, Hnorm, Mnorm
      DOUBLE PRECISION :: c1, c2, c3

      DOUBLE PRECISION :: u_ea(3), u_oa_1(3), u_oa_2(3), chi_o, Mrem_norm

      DOUBLE PRECISION, PARAMETER :: eps_NR = 1.0d-8
      INTEGER :: i1, i2
      DOUBLE PRECISION :: alpha

      ! Initial guess of H
      H_out = H_app
        
      ! Material info
      sdex  = state_dex(i_tile)
      stype = state_type(sdex)
      DO iter = 1, maxiter  
        CALL mumaterial_get_M(sdex, H_out, M_out, chi)
        Hnorm = NORM2(H_out)

        SELECT CASE (stype)
          CASE (STATE_LINEAR)
            dMdH_mat = chi * I3
          CASE (STATE_SOFT)
            IF (Hnorm < 1.0d-12) THEN
              dMdH_mat = 0.0d0
            ELSE
              Mnorm = NORM2(M_out)
              c1 = chi - Mnorm/Hnorm  ! anisotropic part
              c2 = Mnorm/Hnorm          ! isotropic part
              c3 = c1/Hnorm**2

              dMdH_mat(1,1) = c3*H_out(1)*H_out(1) + c2
              dMdH_mat(2,1) = c3*H_out(2)*H_out(1)
              dMdH_mat(3,1) = c3*H_out(3)*H_out(1)

              dMdH_mat(1,2) = c3*H_out(1)*H_out(2)
              dMdH_mat(2,2) = c3*H_out(2)*H_out(2) + c2
              dMdH_mat(3,2) = c3*H_out(3)*H_out(2)

              dMdH_mat(1,3) = c3*H_out(1)*H_out(3)
              dMdH_mat(2,3) = c3*H_out(2)*H_out(3)
              dMdH_mat(3,3) = c3*H_out(3)*H_out(3) + c2
            END IF      
          CASE (STATE_HARD)
            ! Get vectors
            Mrem_norm = NORM2(Mrem(:,sdex))
            u_ea = Mrem(:,sdex)/Mrem_norm ! Easy axis assumed parallel to remanent magnetization
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
            ! Get dMdH matrix
            chi    = constant_mu(sdex) - 1.0d0
            chi_o  = constant_mu_o(sdex) - 1.0d0
            DO i1 = 1, 3
              DO i2 = 1, 3
                dMdH_mat(i1,i2) = chi    * u_ea(i1)*u_ea(i2) &
                              + chi_o  * u_oa_1(i1)*u_oa_1(i2) &
                              + chi_o  * u_oa_2(i1)*u_oa_2(i2)
              END DO
            END DO
        END SELECT
      
        ! Jacobian: JAC = I - N_self * dMdH_mat
        JAC(1,1) = 1.0d0 - N_self(1,1)*dMdH_mat(1,1) - N_self(1,2)*dMdH_mat(2,1) - N_self(1,3)*dMdH_mat(3,1)
        JAC(1,2) =       - N_self(1,1)*dMdH_mat(1,2) - N_self(1,2)*dMdH_mat(2,2) - N_self(1,3)*dMdH_mat(3,2)
        JAC(1,3) =       - N_self(1,1)*dMdH_mat(1,3) - N_self(1,2)*dMdH_mat(2,3) - N_self(1,3)*dMdH_mat(3,3)

        JAC(2,1) =       - N_self(2,1)*dMdH_mat(1,1) - N_self(2,2)*dMdH_mat(2,1) - N_self(2,3)*dMdH_mat(3,1)
        JAC(2,2) = 1.0d0 - N_self(2,1)*dMdH_mat(1,2) - N_self(2,2)*dMdH_mat(2,2) - N_self(2,3)*dMdH_mat(3,2)
        JAC(2,3) =       - N_self(2,1)*dMdH_mat(1,3) - N_self(2,2)*dMdH_mat(2,3) - N_self(2,3)*dMdH_mat(3,3)

        JAC(3,1) =       - N_self(3,1)*dMdH_mat(1,1) - N_self(3,2)*dMdH_mat(2,1) - N_self(3,3)*dMdH_mat(3,1)
        JAC(3,2) =       - N_self(3,1)*dMdH_mat(1,2) - N_self(3,2)*dMdH_mat(2,2) - N_self(3,3)*dMdH_mat(3,2)
        JAC(3,3) = 1.0d0 - N_self(3,1)*dMdH_mat(1,3) - N_self(3,2)*dMdH_mat(2,3) - N_self(3,3)*dMdH_mat(3,3)

        ! H_self = N_self * M_out
        H_self(1) = N_self(1,1)*M_out(1) + N_self(1,2)*M_out(2) + N_self(1,3)*M_out(3)
        H_self(2) = N_self(2,1)*M_out(1) + N_self(2,2)*M_out(2) + N_self(2,3)*M_out(3)
        H_self(3) = N_self(3,1)*M_out(1) + N_self(3,2)*M_out(2) + N_self(3,3)*M_out(3)  

        ! Residual
        F = H_out - H_app - H_self

        ! Convergence check
        IF (NORM2(F) < eps_NR * (Hnorm + 1.0d-12)) EXIT
      
        ! Solve J*dH = -F
        CALL SOLVE_3x3(JAC, -F, dH)
        
        alpha = 1.0d0
        DO WHILE (DOT_PRODUCT(H_out + alpha*dH, H_out) < 0.0d0)
          alpha = alpha * 0.5d0
          IF (alpha < 1.0d-10) EXIT
        END DO
        H_out = H_out + alpha * dH
      END DO

      CONTAINS

      PURE SUBROUTINE SOLVE_3x3(MAT, b, x)
!!-------------------------------------------------------------------
!! Solves 3x3 linear system A*x = b using Gaussian elimination
!! with partial pivoting
!!-------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in)  :: MAT(3,3)
      DOUBLE PRECISION, INTENT(in)  :: b(3)
      DOUBLE PRECISION, INTENT(out) :: x(3)

      DOUBLE PRECISION :: AM(3,4)  ! augmented matrix [A|b]
      DOUBLE PRECISION :: tmp, factor
      INTEGER :: i, j, k, p

      ! Build augmented matrix
      AM(1:3,1:3) = MAT
      AM(1,4) = b(1); AM(2,4) = b(2); AM(3,4) = b(3)

      ! Forward elimination with partial pivoting
      DO k = 1, 3
        ! Find pivot row
        p = k
        DO i = k+1, 3
          IF (ABS(AM(i,k)) > ABS(AM(p,k))) p = i
        END DO
        ! Swap rows k and p
        IF (p /= k) THEN
          DO j = k, 4
            tmp = AM(k,j); AM(k,j) = AM(p,j); AM(p,j) = tmp
          END DO
        END IF
        ! Eliminate below pivot
        DO i = k+1, 3
          IF (ABS(AM(k,k)) < 1.0d-300) CYCLE
          factor = AM(i,k) / AM(k,k)
          DO j = k, 4
            AM(i,j) = AM(i,j) - factor*AM(k,j)
          END DO
        END DO
      END DO

      ! Back substitution
      x(3) = AM(3,4) / AM(3,3)
      x(2) = (AM(2,4) - AM(2,3)*x(3)) / AM(2,2)
      x(1) = (AM(1,4) - AM(1,3)*x(3) - AM(1,2)*x(2)) / AM(1,1)

      END SUBROUTINE SOLVE_3x3

      END SUBROUTINE mumaterial_solve_MH_NR
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


      PURE SUBROUTINE mumaterial_get_M(sdex, H, Mloc, chi)
!!-----------------------------------------------------------------
!! Get M as a function of the H-field and material type
!!-----------------------------------------------------------------
      INTEGER, INTENT(in) :: sdex !! State dex
      DOUBLE PRECISION, INTENT(in) :: H(3) !! H-field
      DOUBLE PRECISION, INTENT(out) :: Mloc(3) !! Output M
      DOUBLE PRECISION, INTENT(out) :: chi !! Susceptibility
      INTEGER :: stype
      chi = 0.0d0
      stype = state_type(sdex)

      ! Call relevant function
      SELECT CASE (stype)
        CASE (STATE_LINEAR)
          CALL get_M_linear(H,Mloc,chi)
        CASE (STATE_SOFT)
          CALL get_M_soft(H,Mloc,chi)
        CASE (STATE_HARD)
          CALL get_M_hard(H,Mloc)
      END SELECT

      CONTAINS

      PURE SUBROUTINE get_M_linear(H_in,M_out,chi_out)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: H_in(3)
      DOUBLE PRECISION, INTENT(out) :: M_out(3)
      DOUBLE PRECISION, INTENT(out) :: chi_out

      chi_out = constant_mu(sdex) - 1.0d0
      M_out    = chi_out * H_in
      
      END SUBROUTINE get_M_linear

      PURE SUBROUTINE get_M_hard(H_in, M_out) 
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: H_in(3)
      DOUBLE PRECISION, INTENT(out) :: M_out(3)
      DOUBLE PRECISION :: u_ea(3), u_oa_1(3), u_oa_2(3), Mrem_norm

      Mrem_norm = NORM2(Mrem(:,sdex))
      u_ea = Mrem(:,sdex)/Mrem_norm ! Easy axis assumed parallel to remanent magnetization
      IF (u_ea(2)/=0 .OR. u_ea(3)/=0) THEN     
         ! cross(u_ea ,[1, 0, 0]) and cross (u_ea, that vector)
        u_oa_1 = [0.d0, u_ea(3), -u_ea(2)]
        u_oa_2 = [-u_ea(2)*u_ea(2) - u_ea(3)*u_ea(3), u_ea(1)*u_ea(2), u_ea(1)*u_ea(3)]
      ELSE                                     
         ! Cross(u_ea,[0, 1, 0]) and cross (u_ea, that vector)
        u_oa_1 = [-u_ea(3), 0.d0, u_ea(1)]
        u_oa_2 = [u_ea(1)*u_ea(2), -u_ea(1)*u_ea(1) - u_ea(3)*u_ea(3), u_ea(2)*u_ea(3)]
      END IF

      ! Normalize unit vectors
      u_oa_1 = u_oa_1/NORM2(u_oa_1)
      u_oa_2 = u_oa_2/NORM2(u_oa_2)

      ! Determine magnetization taking into account easy axis
      M_out = (Mrem_norm + (constant_mu(sdex) - 1) * DOT_PRODUCT(H_in, u_ea )) * u_ea &
                        + (constant_mu_o(sdex)- 1) * DOT_PRODUCT(H_in, u_oa_1) * u_oa_1 &
                        + (constant_mu_o(sdex)- 1) * DOT_PRODUCT(H_in, u_oa_2) * u_oa_2

      END SUBROUTINE get_M_hard

      PURE SUBROUTINE get_M_soft(H_in,M_out,chi_out)
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in) :: H_in(3)
      DOUBLE PRECISION, INTENT(out) :: M_out(3)
      DOUBLE PRECISION, INTENT(out) :: chi_out
      DOUBLE PRECISION :: Hn, Mn

      Hn = Norm2(H_in)
      CALL getstate_scalar(stateFunction(sdex)%H, stateFunction(sdex)%M, stateFunction(sdex)%b, Hn, Mn, chi_out)
      IF (Hn>1.0d-12) THEN
        M_out = Mn * H_in / Hn
      ELSE
        M_out = 0.0d0
      END IF

      END SUBROUTINE get_M_soft
      
      PURE SUBROUTINE getstate_scalar(fx, fy, b, x, y, dydx)
!!-------------------------------------------------------------------
!! Evaluates cubic Hermite spline at single point x, returns y and dydx
!!-------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in)  :: fx(:), fy(:), b(:), x
      DOUBLE PRECISION, INTENT(out) :: y, dydx

      INTEGER :: k, n, low, mid, high
      DOUBLE PRECISION :: t, dx, hinv, t2, t3
      DOUBLE PRECISION :: fxk, fxk1, fyk, fyk1, bk, bk1

      n = SIZE(fx)

      ! Clamp to range
      IF (x <= fx(1)) THEN
        y = fy(1); dydx = b(1); RETURN
      END IF
      IF (x >= fx(n)) THEN
        y = fy(n)
        ! To remove discontinuities, taper slope linearly to zero over dx beyond H_sat 
        dx = fx(n) - fx(n-1)  ! use last interval as scale
        IF (x < fx(n) + dx) THEN
          dydx = b(n) * (1.0d0 - (x - fx(n))/dx)
        ELSE
          dydx = 0.0d0
        END IF
        dydx = MAX(dydx,0.0d0)
        RETURN
      END IF

      ! Find interval by binary search
      low  = 1
      high = n
      DO WHILE (high - low > 1)
        mid = (low + high) / 2
        IF (x >= fx(mid)) THEN
          low = mid
        ELSE
          high = mid
        END IF
      END DO
      k = low

      fxk  = fx(k);   fxk1 = fx(k+1)
      fyk  = fy(k);   fyk1 = fy(k+1)
      bk   = b(k);    bk1  = b(k+1)

      dx    = fxk1 - fxk
      hinv = 1.0d0 / dx
      t    = (x - fxk) * hinv
      t2   = t*t
      t3   = t2*t

      y = (2.0d0*t3 - 3.0d0*t2 + 1.0d0)*fyk  &
        + (t3 - 2.0d0*t2 + t)*dx*bk            &
        + (-2.0d0*t3 + 3.0d0*t2)*fyk1         &
        + (t3 - t2)*dx*bk1

      dydx = (6.0d0*t2 - 6.0d0*t)*fyk*hinv   &
          + (3.0d0*t2 - 4.0d0*t + 1.0d0)*bk &
          + (-6.0d0*t2 + 6.0d0*t)*fyk1*hinv &
          + (3.0d0*t2 - 2.0d0*t)*bk1

      END SUBROUTINE getstate_scalar

      END SUBROUTINE mumaterial_get_M
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!                MPI SYNCHRONIZATION ROUTINES                  !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined(MPI_OPT)
      SUBROUTINE mumaterial_init_sync_node()
      USE mpi_inc
      IMPLICIT NONE
      
      INTEGER :: ptr, csr, i, inode_loc, lev
      INTEGER :: max_lev_shar, max_lev_world, nnode_lev, nnode_offset_at_lev_shar
      INTEGER, ALLOCATABLE :: nnodes_per_level_world(:)
      INTEGER, POINTER :: shar_id_buf(:)
      INTEGER :: nscalar_send, nscalar_recv
      INTEGER :: win_sync_temp
      
      IF (lverb) WRITE(6,*) "  MUMAT_INIT:  Initializing node synchronization variables"
      
      !-----------------------------
      ! BUILD LEVEL LOOKUP ARRAYS
      !-----------------------------
      CALL MPI_ALLREDUCE(nnode_loc, nnode_shar, 1, MPI_INT, MPI_SUM, comm_shar, ierr_mpi)
      ALLOCATE(nnodes_per_level_shar(0:tree_depth))
      CALL MPI_ALLREDUCE(nnodes_per_level_loc, nnodes_per_level_shar, tree_depth+1, &
                          MPI_INT, MPI_SUM, comm_shar, ierr_mpi)
      
      !-----------------------------
      ! ALLOCATE SEND/RECV BUFFERS
      !-----------------------------
      IF (shar_rank==0) THEN
        ALLOCATE(nnodes_per_level_world(0:tree_depth))
        CALL MPI_ALLREDUCE(nnodes_per_level_shar, nnodes_per_level_world, tree_depth+1, &
                            MPI_INT, MPI_SUM, comm_master, ierr_mpi)
        max_lev_world = MAXVAL(nnodes_per_level_world)
        nscalar_recv = max_lev_world * node_bufsize
        DEALLOCATE(nnodes_per_level_world)
      END IF
      
      max_lev_shar = MAXVAL(nnodes_per_level_shar)
      nscalar_send = max_lev_shar * node_bufsize
      CALL MPI_BCAST(nscalar_recv, 1, MPI_INTEGER, 0, comm_shar, ierr_mpi)
      CALL mpialloc(sync_node_sendbuf, nscalar_send, shar_rank, 0, comm_shar, win_sync_node_sendbuf)
      CALL mpialloc(sync_node_recvbuf, nscalar_recv, shar_rank, 0, comm_shar, win_sync_node_recvbuf)
      
      !-----------------------------
      ! COMPUTE OFFSETS IN SEND BUFFER (per level)
      !-----------------------------
      CALL mpialloc(shar_id_buf, nnode_shar+1, shar_rank, 0, comm_shar, win_sync_temp)
      ALLOCATE(node_offset_shar(nnode_loc))
      CALL MPI_BCAST(max_lev_world, 1, MPI_INTEGER, 0, comm_shar, ierr_mpi)
      ALLOCATE(sync_node_rcounts(master_size, 0:tree_depth), &
                sync_node_displs(master_size,  0:tree_depth), &
                sync_node_unpack_map(max_lev_world, 0:tree_depth))
      
      DO lev = 0, tree_depth
        nnode_lev = 0
        DO inode_loc = 1, nnode_loc
          IF (node_level_loc(inode_loc) == lev) nnode_lev = nnode_lev + 1
        END DO
        CALL MPI_EXSCAN(nnode_lev, nnode_offset_at_lev_shar, 1, MPI_INTEGER, &
                        MPI_SUM, comm_shar, ierr_mpi)
        IF (shar_rank==0) nnode_offset_at_lev_shar = 0
      
        csr = 0
        DO inode_loc = 1, nnode_loc
          IF (node_level_loc(inode_loc) == lev) THEN
            ptr = nnode_offset_at_lev_shar + csr
            shar_id_buf(ptr+1) = node_local_to_global(inode_loc)
            node_offset_shar(inode_loc) = ptr * node_bufsize
            csr = csr + 1
          END IF
        END DO
      
        CALL MPI_BARRIER(comm_shar, ierr_mpi)
        IF (shar_rank==0) THEN
          CALL MPI_ALLGATHER(nnodes_per_level_shar(lev)*node_bufsize, 1, MPI_INTEGER, &
                              sync_node_rcounts(:,lev), 1, MPI_INTEGER, comm_master, ierr_mpi)
          sync_node_displs(1,lev) = 0
          DO i = 2, master_size
            sync_node_displs(i,lev) = sync_node_displs(i-1,lev) + sync_node_rcounts(i-1,lev)
          END DO
          CALL MPI_ALLGATHERV(shar_id_buf, nnodes_per_level_shar(lev), MPI_INTEGER, &
                              sync_node_unpack_map(:,lev), sync_node_rcounts(:,lev)/node_bufsize, &
                              sync_node_displs(:,lev)/node_bufsize, MPI_INTEGER, &
                              comm_master, ierr_mpi)
        END IF
        CALL MPI_BARRIER(comm_shar, ierr_mpi)
      END DO
      
      CALL MPI_BCAST(sync_node_unpack_map, max_lev_world*(tree_depth+1), MPI_INTEGER, 0, comm_shar, ierr_mpi)
      CALL MPI_BCAST(sync_node_rcounts,      master_size*(tree_depth+1), MPI_INTEGER, 0, comm_shar, ierr_mpi)
      CALL MPI_BCAST(sync_node_displs,       master_size*(tree_depth+1), MPI_INTEGER, 0, comm_shar, ierr_mpi)
      
      CALL mpidealloc(shar_id_buf, win_sync_temp)
      
      END SUBROUTINE mumaterial_init_sync_node

      SUBROUTINE mumaterial_sync_nodes(lev)
!!-----------------------------------------------------------------------------
!! Synchronizes node_m, node_q etc at a given level
!!-----------------------------------------------------------------------------
      USE mpi_inc
      IMPLICIT NONE
      
      INTEGER, INTENT(in) :: lev !! Level to synchronize at
      INTEGER :: inode_loc, i, inode, ptr, nnodes_lev
      
      ! All local ranks fill the send buffer with their nodes at this level
      DO inode_loc = 1, nnode_loc
        IF (node_level_loc(inode_loc) /= lev) CYCLE
        ptr = node_offset_shar(inode_loc)
        inode = node_local_to_global(inode_loc)
        sync_node_sendbuf(ptr+1:ptr+3)   = node_m(:, inode)
        sync_node_sendbuf(ptr+4:ptr+6)   = node_p(:, inode)
        sync_node_sendbuf(ptr+7)         = node_w(inode)
        sync_node_sendbuf(ptr+8:ptr+10)  = node_com(:, inode)
        sync_node_sendbuf(ptr+11:ptr+15) = node_q(:, inode)
      END DO
      
      CALL MPI_BARRIER(comm_shar, ierr_mpi)
      IF (shar_rank == 0) THEN ! Masters communicate
        CALL MPI_ALLGATHERV(sync_node_sendbuf, nnodes_per_level_shar(lev)*node_bufsize, &
                            MPI_DOUBLE_PRECISION, sync_node_recvbuf, &
                            sync_node_rcounts(:,lev), sync_node_displs(:,lev), &
                            MPI_DOUBLE_PRECISION, comm_master, ierr_mpi)
      END IF
      CALL MPI_BARRIER(comm_shar, ierr_mpi)
      
      ! Unpacking
      nnodes_lev = SUM(sync_node_rcounts(:,lev))/node_bufsize
      DO i = 1, nnodes_lev
        inode = sync_node_unpack_map(i,lev)
        ptr = (i-1) * node_bufsize
        node_m(:,inode)   = sync_node_recvbuf(ptr+1:ptr+3)
        node_p(:,inode)   = sync_node_recvbuf(ptr+4:ptr+6)
        node_w(inode)     = sync_node_recvbuf(ptr+7)
        node_com(:,inode) = sync_node_recvbuf(ptr+8:ptr+10)
        node_q(:,inode)   = sync_node_recvbuf(ptr+11:ptr+15)
      END DO
      CALL MPI_BARRIER(comm_shar, ierr_mpi)
      
      END SUBROUTINE mumaterial_sync_nodes
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE mumaterial_init_sync_node_sparse()
!!-----------------------------------------------------------------------------
!! Builds sparse communication pattern for node synchronization.
!! Only communicates nodes that actually appear in interaction lists.
!!-----------------------------------------------------------------------------
      USE mpi_inc
      USE mpi_params, ONLY: MPI_CALC_MYRANGE
      IMPLICIT NONE
      
      INTEGER :: inode, inode_loc, ileaf_loc, csr
      INTEGER :: i, j, ishar
      INTEGER :: n_send_total, n_recv_total
      INTEGER, ALLOCATABLE :: node_needed(:)     ! nnode - flag nodes needed from remote
      INTEGER, ALLOCATABLE :: send_counts(:)     ! master_size - how many to send to each
      INTEGER, ALLOCATABLE :: recv_counts(:)     ! master_size - how many to recv from each
      INTEGER, ALLOCATABLE :: send_counts_buf(:) ! for alltoall
      INTEGER, ALLOCATABLE :: node_owner_shar(:)  ! nnode - which MPI node owns each tree node

      IF (lverb) WRITE(6,*) "  MUMAT_INIT:  Building sparse node sync pattern"

      !-----------------------------
      ! BUILD NODE OWNERSHIP MAP
      !-----------------------------
      ALLOCATE(node_owner_shar(nnode))
      node_owner_shar = -1
      DO inode_loc = 1, nnode_loc
        inode = node_local_to_global(inode_loc)
        node_owner_shar(inode) = shar_id
      END DO

      CALL MPI_ALLREDUCE(MPI_IN_PLACE, node_owner_shar, nnode, &
                          MPI_INTEGER, MPI_MAX, comm_world, ierr_mpi)
      
      !-----------------------------
      ! FIND NECESSARY REMOTE NODES
      !-----------------------------
      ALLOCATE(node_needed(nnode))
      node_needed = 0
      DO ileaf_loc = 1, nleaf_loc
        DO csr = interact_ptr(ileaf_loc), interact_ptr(ileaf_loc+1)-1
          IF (interact_type(csr) == INTERACT_NODE) THEN
            inode = interact_list(csr)
            IF (node_owner_shar(inode) /= shar_id) THEN
              node_needed(inode) = 1
            END IF
          END IF
        END DO
      END DO

      CALL MPI_ALLREDUCE(MPI_IN_PLACE, node_needed, nnode, &
                          MPI_INTEGER, MPI_MAX, comm_world, ierr_mpi)

      !-----------------------------
      ! COUNT NECESSARY NODES PER SHAR DOMAIN
      !-----------------------------
      ALLOCATE(sparse_recv_counts(0:master_size-1))
      ALLOCATE(sparse_send_counts(0:master_size-1))
      sparse_recv_counts = 0
      DO inode = 1, nnode
        IF (node_needed(inode) == 1) THEN
          ishar = node_owner_shar(inode)
          sparse_recv_counts(ishar) = sparse_recv_counts(ishar) + 1
        END IF
      END DO
      
      !-----------------------------
      ! EXCHANGE COUNTS VIA ALLTOALL
      !-----------------------------
      IF (shar_rank == 0) THEN
        CALL MPI_ALLTOALL(sparse_recv_counts, 1, MPI_INTEGER, &
                          sparse_send_counts, 1, MPI_INTEGER, &
                          comm_master, ierr_mpi)
      END IF
      CALL MPI_BCAST(sparse_send_counts, master_size, MPI_INTEGER, 0, comm_shar, ierr_mpi)
      
      !-----------------------------
      ! BUILD RECV MAP (needed remote nodes, grouped by shar domain)
      !-----------------------------
      n_recv_total = SUM(sparse_recv_counts)
      ALLOCATE(sparse_recv_displs(0:master_size-1))
      sparse_recv_displs(0) = 0
      DO i = 1, master_size-1
        sparse_recv_displs(i) = sparse_recv_displs(i-1) + sparse_recv_counts(i-1)
      END DO
      CALL mpialloc(sparse_recv_map, n_recv_total, shar_rank, 0, comm_shar, win_sparse_recv_map)
      
      ! Fill recv map - nodes needed from each shar domain in order
      IF (shar_rank == 0) THEN
        ALLOCATE(send_counts_buf(0:master_size-1))
        send_counts_buf = 0
        DO inode = 1, nnode
          IF (node_needed(inode) == 1) THEN
            ishar = node_owner_shar(inode)
            j = sparse_recv_displs(ishar) + send_counts_buf(ishar) + 1
            sparse_recv_map(j) = inode
            send_counts_buf(ishar) = send_counts_buf(ishar) + 1
          END IF
        END DO
        DEALLOCATE(send_counts_buf)
      END IF
      CALL MPI_BARRIER(comm_shar, ierr_mpi)
      
      !-----------------------------
      ! BUILD SEND MAP
      !-----------------------------
      n_send_total = SUM(sparse_send_counts)
      ALLOCATE(sparse_send_displs(0:master_size-1))
      sparse_send_displs(0) = 0
      DO i = 1, master_size-1
        sparse_send_displs(i) = sparse_send_displs(i-1) + sparse_send_counts(i-1)
      END DO
      
      CALL mpialloc(sparse_send_map, n_send_total, shar_rank, 0, comm_shar, win_sparse_send_map)
      
      IF (shar_rank == 0) THEN
        CALL MPI_ALLTOALLV(sparse_recv_map, sparse_recv_counts, sparse_recv_displs, MPI_INTEGER, &
                            sparse_send_map, sparse_send_counts, sparse_send_displs, MPI_INTEGER, &
                            comm_master, ierr_mpi)
      END IF

      CALL MPI_BARRIER(comm_shar, ierr_mpi)
      
      !-----------------------------
      ! ALLOCATE SEND/RECV BUFFERS
      !-----------------------------
      CALL mpialloc(sparse_send_buf, n_send_total*node_bufsize, &
                            shar_rank, 0, comm_shar, win_sparse_send_buf)
      CALL mpialloc(sparse_recv_buf, n_recv_total*node_bufsize, &
                            shar_rank, 0, comm_shar, win_sparse_recv_buf)
      
      ! Convert counts/displs to scalar units for MPI calls
      sparse_send_counts = sparse_send_counts * node_bufsize
      sparse_send_displs = sparse_send_displs * node_bufsize
      sparse_recv_counts = sparse_recv_counts * node_bufsize
      sparse_recv_displs = sparse_recv_displs * node_bufsize
      
      CALL MPI_CALC_MYRANGE(comm_shar, 1, n_send_total, sparse_pack_start, sparse_pack_end)
      CALL MPI_CALC_MYRANGE(comm_shar, 1, n_recv_total, sparse_unpack_start, sparse_unpack_end)

      DEALLOCATE(node_needed, node_owner_shar)
      
      END SUBROUTINE mumaterial_init_sync_node_sparse

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE mumaterial_sync_nodes_sparse()
!!-----------------------------------------------------------------------------
!! Sparse synchronization - only communicates nodes in interaction lists
!!-----------------------------------------------------------------------------
      USE mpi_inc
      IMPLICIT NONE

      INTEGER :: i, inode, ptr
      
      CALL MPI_BARRIER(comm_shar, ierr_mpi)

      ! All ranks pack their owned nodes into send buffer
      DO i = sparse_pack_start, sparse_pack_end
        inode = sparse_send_map(i)
        ptr = (i-1)*node_bufsize
        sparse_send_buf(ptr+1:ptr+3)   = node_m(:, inode)
        sparse_send_buf(ptr+4:ptr+6)   = node_p(:, inode)
        sparse_send_buf(ptr+7)         = node_w(inode)
        sparse_send_buf(ptr+8:ptr+10)  = node_com(:, inode)
        sparse_send_buf(ptr+11:ptr+15) = node_q(:, inode)
      END DO
      
      CALL MPI_BARRIER(comm_shar, ierr_mpi)
      IF (shar_rank == 0) THEN ! Masters communicate
        CALL MPI_ALLTOALLV(sparse_send_buf, sparse_send_counts, sparse_send_displs, &
                            MPI_DOUBLE_PRECISION, &
                            sparse_recv_buf, sparse_recv_counts, sparse_recv_displs, &
                            MPI_DOUBLE_PRECISION, comm_master, ierr_mpi)
      END IF
      CALL MPI_BARRIER(comm_shar, ierr_mpi)
      
      ! All ranks unpack received nodes
      DO i = sparse_unpack_start, sparse_unpack_end
        inode = sparse_recv_map(i)
        ptr = (i-1)*node_bufsize
        node_m(:,inode)   = sparse_recv_buf(ptr+1:ptr+3)
        node_p(:,inode)   = sparse_recv_buf(ptr+4:ptr+6)
        node_w(inode)     = sparse_recv_buf(ptr+7)
        node_com(:,inode) = sparse_recv_buf(ptr+8:ptr+10)
        node_q(:,inode)   = sparse_recv_buf(ptr+11:ptr+15)
      END DO
      
      CALL MPI_BARRIER(comm_shar, ierr_mpi)
      
      END SUBROUTINE mumaterial_sync_nodes_sparse
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE mumaterial_init_sync_M()
!!-----------------------------------------------------------------------------
!! Initialize arrays for a full synchronization of M (send counts, recv counts,
!! offsets, map from allgather'd indices to global indices)
!!-----------------------------------------------------------------------------
      USE mpi_inc
      IMPLICIT NONE

      INTEGER :: ileaf_loc, csr, itile, i, ptr
      INTEGER :: nscalar_send, nscalar_recv, ntile_offset
      INTEGER, POINTER :: shar_id_buf(:)
      INTEGER :: win_sync_temp

      IF (lverb) WRITE (6,*) "  MUMAT_INIT:  Initializing M synchronization "; FLUSH(6)

      ! Build send and recv buffers for M
      CALL MPI_ALLREDUCE(ntile_loc, ntile_shar, 1, MPI_INT, MPI_SUM, comm_shar, ierr_mpi)
      nscalar_send = 3*ntile_shar
      nscalar_recv = 3*ntet
      CALL mpialloc(sync_M_sendbuf, nscalar_send, shar_rank, 0,comm_shar, win_sync_M_sendbuf)
      CALL mpialloc(sync_M_recvbuf, nscalar_recv, shar_rank, 0,comm_shar, win_sync_M_recvbuf)

      ! Get offset in sendbuf
      CALL mpialloc(shar_id_buf, ntile_shar, shar_rank,      0, comm_shar, win_sync_temp)
      CALL MPI_EXSCAN(ntile_loc, ntile_offset, 1, MPI_INTEGER, &
                      MPI_SUM, comm_shar, ierr_mpi)
      CALL MPI_BARRIER(comm_shar, ierr_mpi)

      IF (shar_rank==0) ntile_offset = 0
      ALLOCATE(M_offset_shar(ntile_loc))

      csr = 0
      DO ileaf_loc = 1, nleaf_loc
        DO itile = 1, leaf_size_loc(ileaf_loc)
          ptr = ntile_offset + csr
          csr = csr + 1
          M_offset_shar(csr) = ptr*M_bufsize
          shar_id_buf(ptr+1) = leaf_tile_loc(itile,ileaf_loc)
        END DO
      END DO

      ! Get master offset & counts at every level
      ALLOCATE(sync_M_rcounts(master_size),&
               sync_M_displs(master_size), &
               sync_M_unpack_map(ntet))
      CALL MPI_BARRIER(comm_shar, ierr_mpi)  

      IF (shar_rank==0) THEN
        CALL MPI_ALLGATHER(nscalar_send, 1, MPI_INTEGER, sync_M_rcounts, &
                1, MPI_INTEGER, comm_master, ierr_mpi)
        sync_M_displs(1) = 0
        DO i = 2, master_size
            sync_M_displs(i) = sync_M_displs(i-1) + sync_M_rcounts(i-1)
        END DO
        CALL MPI_ALLGATHERV(shar_id_buf, ntile_shar, MPI_INTEGER, &
            sync_M_unpack_map, sync_M_rcounts/M_bufsize, &
            sync_M_displs/M_bufsize, MPI_INTEGER, &
            comm_master, ierr_mpi)
      END IF

      CALL MPI_Bcast(sync_M_unpack_map, ntet, MPI_INTEGER, 0, comm_shar, ierr_mpi)
      CALL MPI_Bcast(sync_M_rcounts, master_size, MPI_INTEGER, 0, comm_shar, ierr_mpi)
      CALL MPI_Bcast(sync_M_displs, master_size, MPI_INTEGER, 0, comm_shar, ierr_mpi)
      CALL mpidealloc(shar_id_buf, win_sync_temp)

      
      END SUBROUTINE mumaterial_init_sync_M
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE mumaterial_sync_M()
!!-----------------------------------------------------------------------------
!! Synchronization of updated magnetization; all comm_shar ranks fill a sendbuffer of only
!! their nodes, masters do inter-node communication of sendbuffer, then recvbuffer
!! is unpacked
!!-----------------------------------------------------------------------------
      USE mpi_inc
      IMPLICIT NONE

      INTEGER :: i, tile, itile_loc, ileaf_loc, itile_leaf, ptr

      ! All local ranks fill the send buffer in shared memory
      itile_loc = 0
      DO ileaf_loc = 1, nleaf_loc
        DO itile_leaf = 1, leaf_size_loc(ileaf_loc)
          itile_loc = itile_loc + 1
          ptr = M_offset_shar(itile_loc)
          tile = leaf_tile_loc(itile_leaf, ileaf_loc)
          sync_M_sendbuf(ptr+1:ptr+3) = M(:, tile)
        END DO
      END DO

      CALL MPI_BARRIER(comm_shar, ierr_mpi)
      IF (shar_rank==0) THEN ! Inter-node communication
        CALL MPI_ALLGATHERV(sync_M_sendbuf, ntile_shar*M_bufsize, &
                            MPI_DOUBLE_PRECISION, sync_M_recvbuf, &
                            sync_M_rcounts, sync_M_displs, &
                            MPI_DOUBLE_PRECISION, comm_master, ierr_mpi)
      END IF
      CALL MPI_BARRIER(comm_shar, ierr_mpi)

      ! Unpacking
      DO i = 1, ntet
        tile = sync_M_unpack_map(i)
        M(:, tile) = sync_M_recvbuf(3*(i-1)+1:3*(i-1)+3)
      END DO
      CALL MPI_BARRIER(comm_shar, ierr_mpi)

      
      END SUBROUTINE mumaterial_sync_M
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


      SUBROUTINE mumaterial_init_sync_M_sparse()
!!-----------------------------------------------------------------------------
!! Initialize arrays for a sparse synchronization of M (send counts, recv counts,
!! offsets, map from allgather'd indices to global indices)
!!-----------------------------------------------------------------------------
      USE mpi_inc
      USE mpi_params, ONLY: MPI_CALC_MYRANGE
      IMPLICIT NONE

      INTEGER :: ileaf_loc, inode_loc, inode, itile, jtile, ishar
      INTEGER :: i, j, ptr, csr
      INTEGER :: n_send_total, n_recv_total
      INTEGER, ALLOCATABLE :: tile_owner_shar(:)
      INTEGER, ALLOCATABLE :: tile_needed(:)
      INTEGER, ALLOCATABLE :: send_counts_buf(:)

      IF (lverb) WRITE(6,*) "  MUMAT_INIT:  Building sparse M sync pattern"

      !-----------------------------
      ! BUILD TILE OWNERSHIP MAP
      !-----------------------------
      ALLOCATE(tile_owner_shar(ntet))
      tile_owner_shar = -1
      ptr = 0
      DO ishar = 0, master_size-1
        DO i = 1, sync_M_rcounts(ishar+1)/M_bufsize
          jtile = sync_M_unpack_map(ptr+i)
          tile_owner_shar(jtile) = ishar
        END DO
        ptr = ptr + sync_M_rcounts(ishar+1)/M_bufsize
      END DO

      !-----------------------------
      ! FIND NECESSARY REMOTE TILES
      !-----------------------------
      ALLOCATE(tile_needed(ntet))
      tile_needed = 0
      DO ileaf_loc = 1, nleaf_loc
        DO csr = interact_ptr(ileaf_loc), interact_ptr(ileaf_loc+1)-1
          IF (interact_type(csr) == INTERACT_LEAF) THEN
            inode = interact_list(csr)
            DO itile = 1, leaf_size(inode)
              jtile = leaf_tile(itile, inode)
              IF (tile_owner_shar(jtile) /= shar_id) THEN
                tile_needed(jtile) = 1
              END IF
            END DO
          END IF
        END DO
      END DO

      ! Share within shar domain
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, tile_needed, ntet, &
                        MPI_INTEGER, MPI_MAX, comm_shar, ierr_mpi)

      !-----------------------------
      ! COUNT NECESSARY TILES PER SHAR DOMAIN
      !-----------------------------
      ALLOCATE(sparse_M_recv_counts(0:master_size-1))
      ALLOCATE(sparse_M_send_counts(0:master_size-1))
      sparse_M_recv_counts = 0
      DO jtile = 1, ntet
        IF (tile_needed(jtile) == 1) THEN
          ishar = tile_owner_shar(jtile)
          sparse_M_recv_counts(ishar) = sparse_M_recv_counts(ishar) + 1
        END IF
      END DO

      !-----------------------------
      ! EXCHANGE COUNTS
      !-----------------------------
      IF (shar_rank == 0) THEN
        CALL MPI_ALLTOALL(sparse_M_recv_counts, 1, MPI_INTEGER, &
                          sparse_M_send_counts, 1, MPI_INTEGER, &
                          comm_master, ierr_mpi)
      END IF
      CALL MPI_BCAST(sparse_M_send_counts, master_size, MPI_INTEGER, 0, comm_shar, ierr_mpi)
      !-----------------------------
      ! BUILD RECV MAP
      !-----------------------------
      n_recv_total = SUM(sparse_M_recv_counts)
      ALLOCATE(sparse_M_recv_displs(0:master_size-1))
      sparse_M_recv_displs(0) = 0
      DO i = 1, master_size-1
        sparse_M_recv_displs(i) = sparse_M_recv_displs(i-1) + sparse_M_recv_counts(i-1)
      END DO

      CALL mpialloc(sparse_M_recv_map, n_recv_total, shar_rank, 0, &
                          comm_shar, win_sparse_M_recv_map)

      IF (shar_rank == 0) THEN
        ALLOCATE(send_counts_buf(0:master_size-1))
        send_counts_buf = 0
        DO jtile = 1, ntet
          IF (tile_needed(jtile) == 1) THEN
            ishar = tile_owner_shar(jtile)
            j = sparse_M_recv_displs(ishar) + send_counts_buf(ishar) + 1
            sparse_M_recv_map(j) = jtile
            send_counts_buf(ishar) = send_counts_buf(ishar) + 1
          END IF
        END DO
        DEALLOCATE(send_counts_buf)
      END IF
      CALL MPI_BARRIER(comm_shar, ierr_mpi)

      !-----------------------------
      ! BUILD SEND MAP
      !-----------------------------
      n_send_total = SUM(sparse_M_send_counts)
      ALLOCATE(sparse_M_send_displs(0:master_size-1))
      sparse_M_send_displs(0) = 0
      DO i = 1, master_size-1
        sparse_M_send_displs(i) = sparse_M_send_displs(i-1) + sparse_M_send_counts(i-1)
      END DO

      CALL mpialloc(sparse_M_send_map, n_send_total, shar_rank, 0, &
                          comm_shar, win_sparse_M_send_map)

      IF (shar_rank == 0) THEN
        CALL MPI_ALLTOALLV(sparse_M_recv_map, sparse_M_recv_counts, sparse_M_recv_displs, &
                          MPI_INTEGER, &
                          sparse_M_send_map, sparse_M_send_counts, sparse_M_send_displs, &
                          MPI_INTEGER, comm_master, ierr_mpi)
      END IF
      CALL MPI_BARRIER(comm_shar, ierr_mpi)

      !-----------------------------
      ! ALLOCATE BUFFERS
      !-----------------------------
      CALL mpialloc(sparse_M_send_buf, n_send_total*M_bufsize, &
                          shar_rank, 0, comm_shar, win_sparse_M_send_buf)
      CALL mpialloc(sparse_M_recv_buf, n_recv_total*M_bufsize, &
                          shar_rank, 0, comm_shar, win_sparse_M_recv_buf)

      sparse_M_send_counts = sparse_M_send_counts * M_bufsize
      sparse_M_send_displs = sparse_M_send_displs * M_bufsize
      sparse_M_recv_counts = sparse_M_recv_counts * M_bufsize
      sparse_M_recv_displs = sparse_M_recv_displs * M_bufsize

      CALL MPI_CALC_MYRANGE(comm_shar, 1, n_send_total, sparse_M_pack_start, sparse_M_pack_end)
      CALL MPI_CALC_MYRANGE(comm_shar, 1, n_recv_total, sparse_M_unpack_start, sparse_M_unpack_end)

      DEALLOCATE(tile_needed, tile_owner_shar)

      END SUBROUTINE mumaterial_init_sync_M_sparse
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


      SUBROUTINE mumaterial_sync_M_sparse()
!!-----------------------------------------------------------------------------
!! Sparse synchronization of M; does an ALLTOALL of only the necessary M to avoid
!! a full, expensive synchronization
!!-----------------------------------------------------------------------------
      USE mpi_inc
      IMPLICIT NONE

      INTEGER :: i, jtile, ptr

      CALL MPI_BARRIER(comm_shar, ierr_mpi)

      ! Pack send buffer
      DO i = sparse_M_pack_start, sparse_M_pack_end
        jtile = sparse_M_send_map(i)
        ptr = (i-1)*M_bufsize
        sparse_M_send_buf(ptr+1:ptr+3) = M(:, jtile)
      END DO

      CALL MPI_BARRIER(comm_shar, ierr_mpi)
      IF (shar_rank == 0) THEN ! Masters communicate
        CALL MPI_ALLTOALLV(sparse_M_send_buf, sparse_M_send_counts, sparse_M_send_displs, &
                          MPI_DOUBLE_PRECISION, &
                          sparse_M_recv_buf, sparse_M_recv_counts, sparse_M_recv_displs, &
                          MPI_DOUBLE_PRECISION, comm_master, ierr_mpi)
      END IF
      CALL MPI_BARRIER(comm_shar, ierr_mpi)

      ! Unpacking
      DO i = sparse_M_unpack_start, sparse_M_unpack_end
        jtile = sparse_M_recv_map(i)
        ptr = (i-1)*M_bufsize
        M(:, jtile) = sparse_M_recv_buf(ptr+1:ptr+3)
      END DO

      CALL MPI_BARRIER(comm_shar, ierr_mpi)

      END SUBROUTINE mumaterial_sync_M_sparse

#endif

      SUBROUTINE mumaterial_propagate_nodes(lv)
!!-----------------------------------------------------------------------
!! Propagates node quantities from leaves upwards.
!!-----------------------------------------------------------------------
      USE mpi_inc
      
      LOGICAL, OPTIONAL, INTENT(in) :: lv
      INTEGER :: i, lev, inode, c, child, itile, inode_loc, tile, ntile
      DOUBLE PRECISION :: msum(3), psum(3), wsum, qsum(5)
      DOUBLE PRECISION :: mt(3), mnorm
      DOUBLE PRECISION :: t0, t1
      
      IF (PRESENT(lv).AND.lv) WRITE(6,*) "  MUMAT_INIT:  Propagating node quantities"
#if defined (MPI_OPT)
      t0 = MPI_WTIME()
      IF (lcomm) CALL MPI_BARRIER(comm_shar, ierr_mpi)
#endif
      ! From bottom to top
      DO lev = tree_depth, 0, -1
        DO i = 1, nnodes_per_level_loc(lev)
          inode_loc = nodes_level_loc(i, lev)
          inode = node_local_to_global(inode_loc)
          msum = 0.0d0
          psum = 0.0d0
          wsum = 0.0d0
          qsum = 0.0d0
          !-------------------
          ! LEAF
          !-------------------
          IF (node_child_loc(1,inode_loc)==NODE_NOCHILD) THEN
            ntile = leaf_size(inode)
            DO itile = 1, ntile
              tile = leaf_tile(itile, inode)
              mt = M(:,tile)*tet_vol(tile)
              mnorm = NORM2(mt)
              msum = msum + mt
              psum = psum + tet_cen(:,tile)*mnorm
              wsum = wsum + mnorm
            END DO
          ELSE
          !-------------------
          ! AGGREGATE NODE
          !-------------------
            DO c = 1, 8
              child = node_child_loc(c,inode_loc)
              IF (child == NODE_NOCHILD) CYCLE
              msum = msum + node_m(:, child)
              psum = psum + node_p(:, child)
              wsum = wsum + node_w(child)
            END DO
          END IF
          node_m(:, inode) = msum
          node_p(:, inode) = psum
          node_w(inode) = wsum
          IF (wsum > 0d0) THEN
            node_com(:,inode) = psum / wsum
          ELSE
            node_com(:,inode) = node_cen(:,inode)
          END IF
          !-------------------
          ! QUADRUPOLE TERM
          !-------------------
          IF (node_child_loc(1,inode_loc)==NODE_NOCHILD) THEN
            ntile = leaf_size(inode)
            DO itile = 1, ntile
              tile = leaf_tile(itile, inode)
              CALL accumulate_quadrupole(qsum, zero5, M(:,tile)*tet_vol(tile),tet_cen(:,tile)-node_com(:,inode))
            END DO
          ELSE
            DO c = 1, 8
              child = node_child_loc(c,inode_loc)
              CALL accumulate_quadrupole(qsum, node_q(:,child), node_m(:,child),node_com(:,child)-node_com(:,inode))
            END DO
          END IF
          node_q(:, inode) = qsum
        END DO

#if defined(MPI_OPT)
        t1 = MPI_WTIME(); t_propagate = t_propagate + (t1-t0); t0 = MPI_WTIME()
        IF (ldosync) CALL mumaterial_sync_nodes_sparse()
        t1 = MPI_WTIME(); t_sync_nodes = t_sync_nodes + (t1-t0); t0 = MPI_WTIME()
#endif
      END DO

      CONTAINS
      PURE SUBROUTINE accumulate_quadrupole(qs, qc, mc, r)
!!-----------------------------------------------------------------------
!! Helper subroutine for calculating quadrupole tensor
!!-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(inout) :: qs(5) !! Output Q-tensor entries
      DOUBLE PRECISION, INTENT(in)    :: qc(5) !! Q of the child node (0 if element)
      DOUBLE PRECISION, INTENT(in)    :: mc(3) !! Dipole moment of node/element
      DOUBLE PRECISION, INTENT(in)    :: r(3)  !! Distance vector
      DOUBLE PRECISION :: mr

      mr = mc(1)*r(1) + mc(2)*r(2) + mc(3)*r(3)
      qs(1) = qs(1) + qc(1) + 3.0d0*mc(1)*r(1) - mr
      qs(2) = qs(2) + qc(2) + 3.0d0*mc(2)*r(2) - mr
      qs(3) = qs(3) + qc(3) + 1.5d0*(mc(1)*r(2) + mc(2)*r(1))
      qs(4) = qs(4) + qc(4) + 1.5d0*(mc(1)*r(3) + mc(3)*r(1))
      qs(5) = qs(5) + qc(5) + 1.5d0*(mc(2)*r(3) + mc(3)*r(2))
      

      END SUBROUTINE accumulate_quadrupole
      
      END SUBROUTINE mumaterial_propagate_nodes

!--------------------------------------------------------------------
!--------------------------------------------------------------------

      SUBROUTINE mumaterial_magfile_read(filename)
!!-----------------------------------------------------------------------
!! Reads a magnetization file, containing 3 columns of M of each element
!!-----------------------------------------------------------------------
      USE mpi_inc
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER :: i, istat, iunit

      IF (lismaster) THEN

        iunit = 327; istat = 0
        CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
        IF (istat/= 0) THEN
          WRITE(6,'(/,A)') "  WARNING: Could not read magfile!"
          WRITE(6,'(A,/)') "  WARNING: Defaulting to zero magnetization!"
          RETURN
        END IF
        DO i = 1, ntet
          READ(iunit, *) M(:,i)
        END DO
        CLOSE(iunit)
      END IF

#if defined(MPI_OPT)
      IF ((lcomm).AND.(shar_rank.EQ.0)) THEN  ! Broadcast to masters
        CALL MPI_Bcast(M,3*ntet,MPI_DOUBLE_PRECISION,0,comm_master,ierr_mpi)
      END IF
#endif
      END SUBROUTINE
!--------------------------------------------------------------------
!-----------------------------------------------------------------------

      SUBROUTINE mumaterial_magfile_write(str)
!!-----------------------------------------------------------------------
!! Writes M to a .dat file, containing 3 columns of M of each element
!!-----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: str
      CHARACTER(LEN=256) :: filename
      INTEGER :: i, istat

      IF (lismaster) THEN
        filename = './mumat_mag_'//TRIM(str)//'.dat'
        WRITE(6,"(A)") "  MUMAT: Writing magnetization to " // filename
        OPEN(13, file=filename,iostat=istat)
        IF (istat /= 0) THEN
          WRITE(6,*) "ERROR: Could not open" // filename // "for writing."
          RETURN
        END IF
        DO i = 1, ntet
          WRITE(13, "(3E15.7)") M(:,i)
        END DO
        CLOSE(13)
      END IF

      END SUBROUTINE
!--------------------------------------------------------------------
!--------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! OUTPUT SUBROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE mumaterial_output(path, x, y, z, linclvac)
!!-----------------------------------------------------------------------
!! Outputs the field at the requested grid points.
!!-----------------------------------------------------------------------
      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(in) :: path !! Path to store data 
      DOUBLE PRECISION, INTENT(in) :: x(:), y(:), z(:) !! Coordinates 
      LOGICAL, INTENT(in) :: linclvac !! True to include background field

      INTEGER :: i, npoints, iunit
      DOUBLE PRECISION, ALLOCATABLE :: B(:,:)

      IF (lismaster) THEN
        npoints = size(x)
        iunit=13
        OPEN(iunit, file='./points.dat')
        DO i = 1, npoints
          WRITE(iunit, "(F15.7,A,F15.7,A,F15.7)") x(i), ',', y(i), ',', z(i)
        END DO
        CLOSE(iunit)
      END IF

      CALL mumaterial_getb_vector(x, y, z, B, linclvac)

      IF (lismaster) THEN
        iunit=14
        OPEN(iunit, file='./B.dat')
        DO i = 1, npoints
          WRITE(iunit, *) B(:,i)
        END DO
        CLOSE(iunit)
      END IF

      
      END SUBROUTINE mumaterial_output
              
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        
      SUBROUTINE mumaterial_getb_vector(x, y, z, B, linclvac)
!!-----------------------------------------------------------------------
!! Evaluates the total magnetic field at the requested coordinates
!!-----------------------------------------------------------------------
      USE mpi_inc
#if defined (MPI_OPT)
      USE mpi_params, ONLY : MPI_CALC_MYRANGE
#endif
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in) :: x(:)
      DOUBLE PRECISION, INTENT(in) :: y(:)
      DOUBLE PRECISION, INTENT(in) :: z(:)
      DOUBLE PRECISION, INTENT(out), ALLOCATABLE :: B(:,:)
      LOGICAL, INTENT(in), OPTIONAL :: linclvac

      DOUBLE PRECISION :: r_eval(3), H(3)
      DOUBLE PRECISION :: Bextx, Bexty, Bextz
      INTEGER :: i, j, ipt, ibox_loc, ibox, npoints, npoints_loc, npts_box
      LOGICAL :: ldovac

      DOUBLE PRECISION, ALLOCATABLE :: r_box(:,:), H_box(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: B_tmp(:,:)
      INTEGER, ALLOCATABLE :: pts_loc(:)

#if defined(MPI_OPT)
      INTEGER, ALLOCATABLE :: send_counts(:), recv_displs(:)
      INTEGER, ALLOCATABLE :: pts_send_counts(:), pts_recv_displs(:)
      DOUBLE PRECISION, ALLOCATABLE :: B_recv(:)
      INTEGER, ALLOCATABLE :: pts_recv(:)
      INTEGER :: total_recv
#endif

      ldovac = .FALSE.
      IF (PRESENT(linclvac)) ldovac = linclvac

      npoints = SIZE(x)
      IF (lverb) WRITE(6,'(A,I0,A)') "  MUMAT: Dividing field evaluation work (", &
                                      npoints, " points)"

      !-----------------------------
      ! 1. BUILD BOXES
      !-----------------------------
      CALL mumaterial_init_boxes(x, y, z, npoints)

      !-----------------------------
      ! 2. COUNT LOCAL POINTS
      !-----------------------------
      npoints_loc = 0
      max_pts_per_box = 0
      DO ibox_loc = 1, neval_boxes_loc
        ibox = eval_local_boxes(ibox_loc)
        npts_box = box_pts_ptr(ibox+1) - box_pts_ptr(ibox)
        npoints_loc = npoints_loc + npts_box
        max_pts_per_box = MAX(max_pts_per_box, npts_box)
      END DO

      !-----------------------------
      ! 3. EVALUATE FIELD FOR LOCAL POINTS ONLY
      !-----------------------------
      ALLOCATE(B_tmp(3, npoints_loc), &
               pts_loc(npoints_loc))
      ALLOCATE(r_box(3, max_pts_per_box), &
               H_box(3, max_pts_per_box))
      B_tmp = 0.0d0
      pts_loc = 0

      i = 0
      DO ibox_loc = 1, neval_boxes_loc
        ibox = eval_local_boxes(ibox_loc)
        npts_box = box_pts_ptr(ibox+1) - box_pts_ptr(ibox)
        DO j = 1, npts_box
          ipt = box_pts(box_pts_ptr(ibox)+j-1)
          pts_loc(i+j) = ipt
          r_box(1,j) = x(ipt)
          r_box(2,j) = y(ipt)
          r_box(3,j) = z(ipt)
        END DO
        H_box(:, 1:npts_box) = 0.0d0


        CALL mumaterial_eval_field_output(r_box, npts_box, ibox_loc, H_box)
          
        DO j = 1, npts_box
          ipt = pts_loc(i+j)
          IF (ldovac) THEN
            CALL ext_Bfld(x(ipt), y(ipt), z(ipt), Bextx, Bexty, Bextz)
            B_tmp(1,i+j) = mu0*H_box(1,j) + Bextx
            B_tmp(2,i+j) = mu0*H_box(2,j) + Bexty
            B_tmp(3,i+j) = mu0*H_box(3,j) + Bextz
          ELSE
            B_tmp(:,i+j) = mu0*H_box(:,j)
          END IF
        END DO

        IF (lverb) THEN
          IF (npoints_loc > 0 .AND. &
              MOD(i*100/npoints_loc, 5) == 0) THEN
            WRITE(6,'(A,I3,A)',ADVANCE='NO') CHAR(13)//'  MUMAT: Field evaluation ', &
              (i+npts_box)*100/npoints_loc, '% done'
            CALL FLUSH(6)
          END IF
        END IF
        i = i + npts_box
      END DO

      IF (lverb) WRITE(6,*) "  MUMAT: Field evaluation finished"
      DEALLOCATE(r_box, H_box)
#if defined(MPI_OPT)
      IF (lcomm) THEN
        !-----------------------------
        ! 4. GATHER TO SHAR_RANK==0
        !-----------------------------
        ALLOCATE(send_counts(shar_size))
        ALLOCATE(recv_displs(shar_size))
        ALLOCATE(pts_send_counts(shar_size))
        ALLOCATE(pts_recv_displs(shar_size))

        ! Exchange point counts on ALL ranks
        CALL MPI_GATHER(npoints_loc, 1, MPI_INTEGER, &
                        pts_send_counts, 1, MPI_INTEGER, &
                        0, comm_shar, ierr_mpi)

        ! Broadcast so all ranks can compute displacements
        CALL MPI_BCAST(pts_send_counts, shar_size, MPI_INTEGER, &
                      0, comm_shar, ierr_mpi)

        ! Build displacements on ALL ranks
        pts_recv_displs(1) = 0
        DO i = 2, shar_size
          pts_recv_displs(i) = pts_recv_displs(i-1) + pts_send_counts(i-1)
        END DO
        total_recv = SUM(pts_send_counts)
        send_counts = pts_send_counts * 3
        recv_displs(1) = 0
        DO i = 2, shar_size
          recv_displs(i) = recv_displs(i-1) + send_counts(i-1)
        END DO

        ! Only allocate receive buffers on shar_rank==0
        IF (shar_rank == 0) THEN
          ALLOCATE(B_recv(3*total_recv))
          ALLOCATE(pts_recv(total_recv))
        ELSE
          ALLOCATE(B_recv(1))
          ALLOCATE(pts_recv(1))
        END IF

        ! Gather B_tmp and pts_loc to shar_rank==0
        CALL MPI_GATHERV(B_tmp, npoints_loc*3, MPI_DOUBLE_PRECISION, &
                        B_recv, send_counts, recv_displs, &
                        MPI_DOUBLE_PRECISION, 0, comm_shar, ierr_mpi)

        CALL MPI_GATHERV(pts_loc, npoints_loc, MPI_INTEGER, &
                        pts_recv, pts_send_counts, pts_recv_displs, &
                        MPI_INTEGER, 0, comm_shar, ierr_mpi)

        DEALLOCATE(B_tmp, pts_loc)
        DEALLOCATE(send_counts, recv_displs)
        DEALLOCATE(pts_send_counts, pts_recv_displs)

        !-----------------------------
        ! 5. SCATTER INTO B ON SHAR_RANK==0
        !-----------------------------
        IF (shar_rank == 0) THEN
          ALLOCATE(B(3, npoints))
          B = 0.0d0
          DO i = 1, total_recv
            ipt = pts_recv(i)
            B(:,ipt) = B_recv((i-1)*3+1 : i*3)
          END DO
          DEALLOCATE(B_recv, pts_recv)

          !-----------------------------
          ! 6. ALLREDUCE BETWEEN MASTERS
          !-----------------------------
          CALL MPI_ALLREDUCE(MPI_IN_PLACE, B, 3*npoints, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, &
                            comm_master, ierr_mpi)
        ELSE
          DEALLOCATE(B_recv, pts_recv)
          ALLOCATE(B(3,1))
          B = 0.0d0
        END IF

      ELSE
#endif
        ! Non-MPI case
        ALLOCATE(B(3, npoints))
        B = 0.0d0
        DO i = 1, npoints_loc
          ipt = pts_loc(i)
          B(:,ipt) = B_tmp(:,i)
        END DO
        DEALLOCATE(B_tmp, pts_loc)
#if defined(MPI_OPT)
      END IF
#endif

      ! Cleanup evaluation arrays
      IF (ALLOCATED(box_npts))         DEALLOCATE(box_npts)
      IF (ALLOCATED(box_pts_ptr))      DEALLOCATE(box_pts_ptr)
      IF (ALLOCATED(box_pts))          DEALLOCATE(box_pts)
      IF (ALLOCATED(eval_local_boxes)) DEALLOCATE(eval_local_boxes)
      IF (ALLOCATED(interact_ptr))     DEALLOCATE(interact_ptr)
      IF (ALLOCATED(interact_list))    DEALLOCATE(interact_list)
      IF (ALLOCATED(interact_type))    DEALLOCATE(interact_type)

      
      END SUBROUTINE mumaterial_getb_vector

      SUBROUTINE mumaterial_init_boxes(x, y, z, npoints)
!!-------------------------------------------------------------------
!! Divides evaluation points into boxes, builds interaction lists,
!! and load-balances boxes across MPI ranks
!!-------------------------------------------------------------------
      USE mpi_inc
#if defined (MPI_OPT)
      USE mpi_params, ONLY : MPI_CALC_MYRANGE
#endif
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in) :: x(npoints), y(npoints), z(npoints)
      INTEGER,          INTENT(in) :: npoints

      ! Box grid
      DOUBLE PRECISION :: box_xmax, box_ymax, box_zmax

      ! Per-box data
      DOUBLE PRECISION, ALLOCATABLE :: box_work(:)
      INTEGER, ALLOCATABLE :: box_nint(:)
      DOUBLE PRECISION, ALLOCATABLE :: box_cen(:,:)

      ! MPI
      INTEGER :: irank

      INTEGER, ALLOCATABLE :: cursor(:)
      DOUBLE PRECISION :: cost_targ, cost_loc
      INTEGER :: nboxes
      INTEGER :: ibox, ibox_loc, ipt, val
      DOUBLE PRECISION :: wval
      INTEGER :: ix, iy, iz, i,j
      INTEGER :: neval_interactions
      INTEGER :: mystart, myend, nboxes_occ
 
      INTEGER, ALLOCATABLE :: box_order(:), box_to_rank(:), local_boxes(:)
      DOUBLE PRECISION, ALLOCATABLE :: rank_work(:)
      INTEGER :: itmp
      INTEGER, PARAMETER :: max_boxes = 100000

      DOUBLE PRECISION :: t0

      t0 = MPI_WTIME()
      !-----------------------------
      ! 1. BOUNDING BOX WITHOUT PADDING
      !-----------------------------
      box_xmin = MINVAL(x); box_xmax = MAXVAL(x)
      box_ymin = MINVAL(y); box_ymax = MAXVAL(y)
      box_zmin = MINVAL(z); box_zmax = MAXVAL(z)

      !-----------------------------
      ! 2. ESTIMATE BOX SIZE
      !-----------------------------
      eval_box_edge = (node_bounds(2,ROOT) - node_bounds(1,ROOT)) / 2.0d0**tree_depth
      DO WHILE (.TRUE.)
        box_nx = CEILING((box_xmax - box_xmin) / eval_box_edge)
        box_ny = CEILING((box_ymax - box_ymin) / eval_box_edge)
        box_nz = CEILING((box_zmax - box_zmin) / eval_box_edge)
        IF (box_nx*box_ny*box_nz <= max_boxes) EXIT
        eval_box_edge = eval_box_edge * 2.0d0
      END DO

      ! Now add proper padding
      box_xmin = box_xmin - 0.5d0*eval_box_edge
      box_xmax = box_xmax + 0.5d0*eval_box_edge
      box_ymin = box_ymin - 0.5d0*eval_box_edge
      box_ymax = box_ymax + 0.5d0*eval_box_edge
      box_zmin = box_zmin - 0.5d0*eval_box_edge
      box_zmax = box_zmax + 0.5d0*eval_box_edge

      ! Recompute with padding
      box_nx = CEILING((box_xmax - box_xmin) / eval_box_edge)
      box_ny = CEILING((box_ymax - box_ymin) / eval_box_edge)
      box_nz = CEILING((box_zmax - box_zmin) / eval_box_edge)
      nboxes = box_nx*box_ny*box_nz

      !-----------------------------
      ! 3. ASSIGN POINTS TO BOXES
      !-----------------------------
      ALLOCATE(box_npts(nboxes))
      box_npts = 0

      DO ipt = 1, npoints
        ix = INT((x(ipt)-box_xmin)/eval_box_edge) + 1
        iy = INT((y(ipt)-box_ymin)/eval_box_edge) + 1
        iz = INT((z(ipt)-box_zmin)/eval_box_edge) + 1
        ibox = ix + box_nx*(iy-1) + box_nx*box_ny*(iz-1)
        box_npts(ibox) = box_npts(ibox) + 1
      END DO

      ALLOCATE(box_pts_ptr(nboxes+1))
      box_pts_ptr(1) = 1
      DO ibox = 1, nboxes
        box_pts_ptr(ibox+1) = box_pts_ptr(ibox) + box_npts(ibox)
      END DO

      ALLOCATE(box_pts(npoints))
      box_npts = 0  ! reuse as cursor
      DO ipt = 1, npoints
        ix = INT((x(ipt)-box_xmin)/eval_box_edge) + 1
        iy = INT((y(ipt)-box_ymin)/eval_box_edge) + 1
        iz = INT((z(ipt)-box_zmin)/eval_box_edge) + 1
        ibox = ix + box_nx*(iy-1) + box_nx*box_ny*(iz-1)
        i = box_pts_ptr(ibox) + box_npts(ibox)
        box_pts(i) = ipt
        box_npts(ibox) = box_npts(ibox) + 1
      END DO
      
      !-----------------------------
      ! 4. COMPUTE BOX CENTERS
      !-----------------------------
      ALLOCATE(box_cen(3,nboxes))
      DO iz = 1, box_nz
        DO iy = 1, box_ny
          DO ix = 1, box_nx
            ibox = ix + box_nx*(iy-1) + box_nx*box_ny*(iz-1)
            box_cen(1,ibox) = box_xmin + (ix-0.5d0)*eval_box_edge
            box_cen(2,ibox) = box_ymin + (iy-0.5d0)*eval_box_edge
            box_cen(3,ibox) = box_zmin + (iz-0.5d0)*eval_box_edge
          END DO
        END DO
      END DO

      t0 = MPI_WTIME()
      !-----------------------------
      ! 5. COUNT INTERACTIONS PER BOX
      !-----------------------------
      ALLOCATE(box_work(nboxes),box_nint(nboxes))
      box_work = 0.0d0
      box_nint = 0
      
      ! Get occupied box count for range distribution
      nboxes_occ = COUNT(box_npts > 0)
      mystart = 1; myend = nboxes

#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_CALC_MYRANGE(comm_world, 1, nboxes, mystart, myend)
#endif
      
      DO ibox = mystart, myend
        IF (box_npts(ibox) == 0) CYCLE
        val = 0
        wval = 0.0d0
        CALL mumaterial_interactions_count(box_cen(:,ibox), ROOT, val, wval,.TRUE.)
        box_nint(ibox) = val
        box_work(ibox) = box_npts(ibox)*wval
      END DO
      
#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_ALLREDUCE(MPI_IN_PLACE, box_work, nboxes, &
                                     MPI_DOUBLE_PRECISION, MPI_SUM, comm_world, ierr_mpi)
      IF (lcomm) CALL MPI_ALLREDUCE(MPI_IN_PLACE, box_nint, nboxes, &
                                     MPI_INTEGER, MPI_SUM, comm_world, ierr_mpi)
#endif

      t0 = MPI_WTIME()
      !-----------------------------
      ! 6. LOAD BALANCE BOXES ACROSS MPI RANKS (LPT)
      !-----------------------------

      ALLOCATE(box_to_rank(nboxes))
      box_to_rank = -1  ! -1 = empty box, not assigned

#if defined(MPI_OPT)
      IF (lcomm) THEN
        ! Sort non-empty boxes by decreasing work
        ALLOCATE(box_order(COUNT(box_work > 0)))
        j = 0
        DO ibox = 1, nboxes
          IF (box_work(ibox) > 0) THEN
            j = j + 1
            box_order(j) = ibox
          END IF
        END DO

        ! Insertion sort by decreasing work
        DO i = 2, SIZE(box_order)
          DO ibox = i, 2, -1
            IF (box_work(box_order(ibox)) > box_work(box_order(ibox-1))) THEN
              itmp = box_order(ibox)
              box_order(ibox) = box_order(ibox-1)
              box_order(ibox-1) = itmp
            ELSE
              EXIT
            END IF
          END DO
        END DO

        ! LPT assignment
        ALLOCATE(rank_work(0:world_size-1))
        rank_work = 0.0d0
        DO i = 1, SIZE(box_order)
          ibox = box_order(i)
          irank = MINLOC(rank_work, 1) - 1
          box_to_rank(ibox) = irank
          rank_work(irank) = rank_work(irank) + box_work(ibox)
        END DO
        DEALLOCATE(box_order, rank_work)
      ELSE
        ! Non-MPI: all boxes to rank 0
        DO ibox = 1, nboxes
          IF (box_work(ibox) > 0) box_to_rank(ibox) = 0
        END DO
      END IF
#else
      DO ibox = 1, nboxes
        IF (box_work(ibox) > 0) box_to_rank(ibox) = 0
      END DO
#endif

      ! Load balance
      cost_loc = 0.0d0
      DO ibox = 1, nboxes
        IF (box_to_rank(ibox) == world_rank) cost_loc = cost_loc + box_work(ibox)
      END DO

      t0 = MPI_WTIME()
      !-----------------------------
      ! 7. BUILD INTERACTION LISTS
      !-----------------------------
      neval_boxes_loc = COUNT(box_to_rank == world_rank)

      ! Build local box list
      ALLOCATE(local_boxes(neval_boxes_loc))
      ibox_loc = 0
      DO ibox = 1, nboxes
        IF (box_to_rank(ibox) == world_rank) THEN
          ibox_loc = ibox_loc + 1
          local_boxes(ibox_loc) = ibox
        END IF
      END DO
      DEALLOCATE(box_to_rank)

      IF (ALLOCATED(interact_ptr))  DEALLOCATE(interact_ptr)
      IF (ALLOCATED(interact_list)) DEALLOCATE(interact_list)
      IF (ALLOCATED(interact_type)) DEALLOCATE(interact_type)

      ALLOCATE(interact_ptr(neval_boxes_loc+1))
      interact_ptr(1) = 1
      DO ibox_loc = 1, neval_boxes_loc
        ibox = local_boxes(ibox_loc)
        interact_ptr(ibox_loc+1) = interact_ptr(ibox_loc) + box_nint(ibox)
      END DO

      neval_interactions = interact_ptr(neval_boxes_loc+1) - 1
      ALLOCATE(interact_list(neval_interactions), &
               interact_type(neval_interactions))

      ALLOCATE(cursor(neval_boxes_loc))
      cursor = interact_ptr(1:neval_boxes_loc)
      DO ibox_loc = 1, neval_boxes_loc
        ibox = local_boxes(ibox_loc)
        IF (box_npts(ibox) == 0) CYCLE
        CALL mumaterial_interactions_fill(box_cen(:,ibox), ROOT, ibox_loc, cursor,.TRUE.)
      END DO

      ALLOCATE(eval_local_boxes(neval_boxes_loc))
      eval_local_boxes = local_boxes

      DEALLOCATE(cursor, box_work, box_nint, box_cen, local_boxes)

      END SUBROUTINE mumaterial_init_boxes
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE mumaterial_eval_field_iter(r, neval, ileaf_loc, H_tiles, iN)
!!------------------------------------------------------------------------------
!! Evaluates field at tiles inside leaf using precomputed interaction list
!!------------------------------------------------------------------------------
      
      IMPLICIT NONE
    
      INTEGER,          INTENT(in)    :: neval
      DOUBLE PRECISION, INTENT(in)    :: r(3,max_leaf_size_seen)
      INTEGER,          INTENT(in)    :: ileaf_loc
      DOUBLE PRECISION, INTENT(inout) :: H_tiles(3,max_leaf_size_seen)
      INTEGER,          INTENT(inout) :: iN
    
      INTEGER :: csr, inode, itype
      INTEGER :: tiles_src(max_leaf_size), ntile_src, itile_loc, itile_src, itile, jtile
      INTEGER :: loc_tiles(max_leaf_size_seen)
      DOUBLE PRECISION :: Mtile(3)
      DOUBLE PRECISION :: com(3), mom(3), q(5)

      loc_tiles(1:neval) = leaf_tile_loc(1:neval, ileaf_loc)    
      DO csr = interact_ptr(ileaf_loc), interact_ptr(ileaf_loc+1)-1
        inode = interact_list(csr)
        itype = interact_type(csr)
        !--------------------------
        ! LEAF -> DEMAG TENSOR
        !--------------------------
        IF (itype==INTERACT_LEAF) THEN
          tiles_src = leaf_tile(:, inode)
          ntile_src = leaf_size(inode)
          DO itile_src = 1, ntile_src
            jtile  = tiles_src(itile_src)
            Mtile  = M(:, jtile)
            DO itile_loc = 1, neval
              itile = loc_tiles(itile_loc)
              IF (itile==jtile) CYCLE
              iN = iN + 1
              H_tiles(1,itile_loc) = H_tiles(1,itile_loc) &
                + Nloc(iN,1,1)*Mtile(1) + Nloc(iN,1,2)*Mtile(2) + Nloc(iN,1,3)*Mtile(3)
              H_tiles(2,itile_loc) = H_tiles(2,itile_loc) &
                + Nloc(iN,2,1)*Mtile(1) + Nloc(iN,2,2)*Mtile(2) + Nloc(iN,2,3)*Mtile(3)
              H_tiles(3,itile_loc) = H_tiles(3,itile_loc) &
                + Nloc(iN,3,1)*Mtile(1) + Nloc(iN,3,2)*Mtile(2) + Nloc(iN,3,3)*Mtile(3)
            END DO
          END DO
    
        !---------------------------
        ! AGGREGATE NODE
        !---------------------------
        ELSE IF (itype==INTERACT_NODE) THEN
          mom = node_m(:, inode)
          com = node_com(:, inode)
          q   = node_q(:, inode)
          CALL mumaterial_eval_node(r, neval, com, mom, q, H_tiles)
        END IF
      END DO
      
      END SUBROUTINE mumaterial_eval_field_iter
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE mumaterial_eval_field_output(r, neval, ibox_loc, H)
!!------------------------------------------------------------------------------
!! Evaluates field at points inside box using precomputed interaction list
!!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      DOUBLE PRECISION, INTENT(in)    :: r(3,neval)  !! Evaluation position
      INTEGER,          INTENT(in)    :: ibox_loc   !! Local box index
      INTEGER,          intent(in)    :: neval !! Number of points
      DOUBLE PRECISION, INTENT(inout) :: H(3,neval)       !! Accumulated field
      
      INTEGER :: csr, inode, itype
      INTEGER :: itile, jtile, ntile_src
      DOUBLE PRECISION :: N(3,3), Mtile(3)

      DOUBLE PRECISION :: com(3), mom(3), q(5)
      
      DO csr = interact_ptr(ibox_loc), interact_ptr(ibox_loc+1)-1
        inode = interact_list(csr)
        itype = interact_type(csr)
        !---------------------------
        ! AGGREGATE NODE
        !---------------------------
        IF (itype == INTERACT_NODE) THEN
          com = node_com(:,inode)
          mom = node_m(:,inode)
          q = node_q(:,inode)
          CALL mumaterial_eval_node(r, neval, com, mom, q, H)
        !---------------------------
        ! LEAF
        !---------------------------
        ELSE IF (itype == INTERACT_LEAF) THEN
          CALL mumaterial_eval_leaf(r, neval, inode, H)
        END IF
      
      END DO
      
      END SUBROUTINE mumaterial_eval_field_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! EVALUATION SUBROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PURE LOGICAL FUNCTION mumaterial_accept_node(r, inode, leval)
!!-------------------------------------------------------------------
!! Returns whether a node should be accepted as aggregate, or it should
!! be recursed into.
!!-------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in) :: r(3)  !! Evaluation position
      INTEGER,          INTENT(in) :: inode !! Node index
      LOGICAL,          INTENT(in) :: leval !! Eval point flag

      DOUBLE PRECISION :: b(6), r_node(3), s, d2, theta2

      b = node_bounds(:,inode)
      r_node = node_cen(:,inode)
      s  = MAX(b(2)-b(1), b(4)-b(3), b(6)-b(5))
      d2 = (r(1)-r_node(1))**2 + (r(2)-r_node(2))**2 + (r(3)-r_node(3))**2

      ! Reject if evaluation point is inside node bounds
      IF (r(1)>=b(1).AND.r(1)<=b(2).AND. &
          r(2)>=b(3).AND.r(2)<=b(4).AND. &
          r(3)>=b(5).AND.r(3)<=b(6)) THEN
        mumaterial_accept_node = .FALSE.
        RETURN
      END IF
      IF (leval) THEN ! Different theta depending on case
        theta2 = theta2_eval
      ELSE
        theta2 = tree_theta2
      ENDIF
      mumaterial_accept_node = (s*s < theta2*d2)

      END FUNCTION mumaterial_accept_node
!--------------------------------------------------------------------
!--------------------------------------------------------------------

      RECURSIVE SUBROUTINE mumaterial_interactions_count(r, inode, counter, work, leval)
!!-----------------------------------------------------------------------
!! Calculate number of interactions (aggregate or leaves) at position r,
!! to determine allocation size of interact_*
!!-----------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in) :: r(3) !! Evaluation position
      INTEGER, INTENT(in) :: inode !! Node to check; goes into children if not accepted
      INTEGER, INTENT(inout) :: counter !! Running total of interaction count
      DOUBLE PRECISION, INTENT(inout) :: work !! Running total of workload
      LOGICAL, INTENT(in) :: leval

      DOUBLE PRECISION :: b(6), r_node(3), s, d2
      INTEGER :: c

      !-------------------
      ! 1. NODE IS A LEAF
      !-------------------
      IF (node_child(1,inode)==NODE_NOCHILD) THEN
        IF (leaf_size(inode) > 0) THEN
          counter = counter + 1
          IF (leval) THEN
            work = work + leaf_size(inode)*WORK_EVAL
          ELSE
            work = work + WORK_TILE !! Simple MATMUL
          END IF
        END IF
        RETURN
      END IF
      !-------------------
      ! 2a. NODE IS FAR -> ACCEPT
      !-------------------
      IF (mumaterial_accept_node(r, inode, leval)) THEN
        counter = counter + 1
        work = work + WORK_NODE
        RETURN
      END IF
      !-------------------
      ! 2b. NODE IS TOO CLOSE -> RECURSE
      !-------------------
      DO c = 1, 8
        CALL mumaterial_interactions_count(r, node_child(c,inode),counter,work,leval)
      END DO
      

      END SUBROUTINE mumaterial_interactions_count
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      RECURSIVE SUBROUTINE mumaterial_interactions_fill(r, inode, i, csr, leval)
!!-----------------------------------------------------------------------
!! Fills interact_list and interact_type with node indices and evaluation
!! types (aggregate or leaf).
!!-----------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in) :: r(3) !! Evaluation position
      INTEGER, INTENT(in) :: inode !! Node to check; goes into children if not accepted
      INTEGER, INTENT(in) :: i !! Index in cursor
      INTEGER, INTENT(inout) :: csr(:) !! Cursor for writing into arrays
      LOGICAL, INTENT(in) :: leval

      INTEGER :: c, pos
      DOUBLE PRECISION :: d2,b(6),s,r_node(3)

      !-------------------
      ! 1. NODE IS A LEAF
      !-------------------
      IF (node_child(1,inode)==NODE_NOCHILD) THEN
        IF(leaf_size(inode)>0) THEN ! Not an empty leaf
          pos = csr(i)
          interact_list(pos) = inode
          interact_type(pos) = INTERACT_LEAF
          pos = pos+1
          csr(i) = pos
        END IF
        RETURN
      END IF
      !-------------------
      ! 2a. NODE IS FAR -> ACCEPT
      !-------------------
      IF (mumaterial_accept_node(r, inode, leval)) THEN
        pos = csr(i)
        interact_list(pos) = inode
        interact_type(pos) = INTERACT_NODE
        csr(i) = pos + 1
        RETURN
      END IF
      !-------------------
      ! 2b. NODE IS TOO CLOSE -> RECURSE
      !-------------------
      DO c = 1, 8
        CALL mumaterial_interactions_fill(r, node_child(c,inode), i, csr, leval)
      END DO
      

      END SUBROUTINE mumaterial_interactions_fill
              
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      PURE SUBROUTINE mumaterial_eval_node(r, neval, com, mom, q, H)
      IMPLICIT NONE
      INTEGER,          INTENT(in)    :: neval
      DOUBLE PRECISION, INTENT(in)    :: r(3,neval), com(3), mom(3), q(5)
      DOUBLE PRECISION, INTENT(inout) :: H(3,neval)

      DOUBLE PRECISION :: dr1, dr2, dr3, r2, r2i, r5inv
      DOUBLE PRECISION :: mr, Qr1, Qr2, Qr3, rQr, f1, f2, g
      DOUBLE PRECISION :: negq12
      INTEGER :: i

      negq12 = -(q(1) + q(2))

      DO i = 1, neval
        dr1  = r(1,i) - com(1)
        dr2  = r(2,i) - com(2)
        dr3  = r(3,i) - com(3)
        r2   = dr1*dr1 + dr2*dr2 + dr3*dr3
        r2i  = 1.0d0 / r2
        r5inv = INV4PI * r2i * r2i * SQRT(r2i)

        mr   = mom(1)*dr1 + mom(2)*dr2 + mom(3)*dr3
        Qr1  = q(1)*dr1 + q(3)*dr2 + q(4)*dr3
        Qr2  = q(3)*dr1 + q(2)*dr2 + q(5)*dr3
        Qr3  = q(4)*dr1 + q(5)*dr2 + negq12*dr3
        rQr  = dr1*Qr1 + dr2*Qr2 + dr3*Qr3

        f1 = 3.0d0 * mr  * r5inv
        f2 = 2.5d0 * rQr * r5inv * r2i
        g  = r2 * r5inv

        H(1,i) = H(1,i) + (f1+f2)*dr1 - g*mom(1) - r5inv*Qr1
        H(2,i) = H(2,i) + (f1+f2)*dr2 - g*mom(2) - r5inv*Qr2
        H(3,i) = H(3,i) + (f1+f2)*dr3 - g*mom(3) - r5inv*Qr3
      END DO
      END SUBROUTINE mumaterial_eval_node



      PURE SUBROUTINE mumaterial_eval_leaf(r, npts, inode, H)
      IMPLICIT NONE
    
      INTEGER,          INTENT(in)    :: npts
      INTEGER,          INTENT(in)    :: inode
      DOUBLE PRECISION, INTENT(in)    :: r(3,npts)
      DOUBLE PRECISION, INTENT(inout) :: H(3,npts)
    
      INTEGER :: i, itile, tile
      DOUBLE PRECISION :: Mtile(3), N(3,3)
    
      IF (node_child(1,inode)/=NODE_NOCHILD) RETURN
    
      DO itile = 1, leaf_size(inode)
        tile  = leaf_tile(itile, inode)
        Mtile = M(:, tile)
        DO i = 1, npts
          N = mumaterial_get_N(tet_P(:,:,:,tile), tet_D(:,:,tile), tet_v(:,:,:,tile), r(:,i))
          H(1,i) = H(1,i) + N(1,1)*Mtile(1) + N(1,2)*Mtile(2) + N(1,3)*Mtile(3)
          H(2,i) = H(2,i) + N(2,1)*Mtile(1) + N(2,2)*Mtile(2) + N(2,3)*Mtile(3)
          H(3,i) = H(3,i) + N(3,1)*Mtile(1) + N(3,2)*Mtile(2) + N(3,3)*Mtile(3)
        END DO
      END DO
    
      END SUBROUTINE mumaterial_eval_leaf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! PYTHON SUBROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------
! mumaterial_get_nvertex: Returns nvertex (for python)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      INTEGER FUNCTION mumaterial_get_nvertex()
      IMPLICIT NONE
      mumaterial_get_nvertex = nvertex
      
      END FUNCTION mumaterial_get_nvertex

!------------------------------------------------------------------------------
! mumaterial_get_ntet: Returns ntet (for python)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      INTEGER FUNCTION mumaterial_get_ntet()
      IMPLICIT NONE
      mumaterial_get_ntet = ntet
      
      END FUNCTION mumaterial_get_ntet

!------------------------------------------------------------------------------
! mumaterial_get_nstate: Returns ntet (for python)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      INTEGER FUNCTION mumaterial_get_nstate()
      IMPLICIT NONE
      mumaterial_get_nstate = nstate
      
      END FUNCTION mumaterial_get_nstate

!------------------------------------------------------------------------------
! mumaterial_get_vertex: Returns vertex (for python)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_get_vertex(vertex_out)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3,nvertex), INTENT(INOUT) :: vertex_out
      vertex_out = vertex
      
      END SUBROUTINE mumaterial_get_vertex

!------------------------------------------------------------------------------
! mumaterial_get_tet: Returns tet array (for python)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_get_tet(tet_out)
      IMPLICIT NONE
      INTEGER, DIMENSION(4,ntet), INTENT(INOUT) :: tet_out
      tet_out = tet
      
      END SUBROUTINE mumaterial_get_tet

!------------------------------------------------------------------------------
! mumaterial_get_statedex: Returns state_dex array (for python)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_get_statedex(state_out)
      IMPLICIT NONE
      INTEGER, DIMENSION(ntet), INTENT(INOUT) :: state_out
      state_out = state_dex
      
      END SUBROUTINE mumaterial_get_statedex

!------------------------------------------------------------------------------
! mumaterial_get_statetype: Returns state_type array (for python)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      SUBROUTINE mumaterial_get_statetype(state_out)
      IMPLICIT NONE
      INTEGER, DIMENSION(nstate), INTENT(INOUT) :: state_out
      state_out = state_type
      
      END SUBROUTINE mumaterial_get_statetype

!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE mumaterial_mod
