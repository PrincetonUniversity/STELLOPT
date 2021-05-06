!-----------------------------------------------------------------------
!     Module:        wall_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov),
!                    D.J. Engels (d.j.engels@student.tue.nl)
!     Date:          May 2021
!     Description:   This module handles defining a wall as a set of
!                    triangular facets which can be used to calculate
!                    if and where a particle hits the mesh.
!                    Note that the user must deallocate face and vertex
!                    after calling wall_load_mn.
!-----------------------------------------------------------------------
      MODULE wall_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE safe_open_mod
      IMPLICIT NONE

!-----------------------------------------------------------------------
!     Types
!         Wall_Type  The full wall, information about the vertices and the number of blocks
!         Block      A block is a part of the uniform grid, an area in space with triangles in it      
!-----------------------------------------------------------------------
      TYPE block
         INTEGER :: nfaces  ! number of faces in block
         LOGICAL :: isshared  ! whether or not shared memory active
         ! shared memory pointsers
         INTEGER :: win_d, win_FN, win_A0, win_V0, win_V1, &
                    win_DOT00, win_DOT11, win_DOT01, win_ihit, win_invDenom
         ! bounds of block
         DOUBLE PRECISION :: xmin, xmax, ymin, ymax, zmin, zmax
         ! non-shared memory pointers
         INTEGER,          DIMENSION(:),   POINTER :: ihit_array => null()
         DOUBLE PRECISION, DIMENSION(:),   POINTER :: d => null()
         DOUBLE PRECISION, DIMENSION(:),   POINTER :: DOT00 => null()
         DOUBLE PRECISION, DIMENSION(:),   POINTER :: DOT01 => null()
         DOUBLE PRECISION, DIMENSION(:),   POINTER :: DOT11 => null()
         DOUBLE PRECISION, DIMENSION(:),   POINTER :: invDenom => null()
         DOUBLE PRECISION, DIMENSION(:,:), POINTER :: FN => null()
         DOUBLE PRECISION, DIMENSION(:,:), POINTER :: A0 => null()
         DOUBLE PRECISION, DIMENSION(:,:), POINTER :: V0 => null()
         DOUBLE PRECISION, DIMENSION(:,:), POINTER :: V1 => null()
      END TYPE block

      TYPE wall_type
         ! number of blocks
         ! and the step that has to be done in the list to move one in x/y/z-direction
         integer :: nblocks, xstep, ystep, zstep
         TYPE (block), DIMENSION(:), POINTER :: blocks => null()
      END TYPE wall_type


!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
      LOGICAL            :: lwall_loaded, lwall_acc
      INTEGER            :: nvertex, nface
      INTEGER, POINTER :: face(:,:)
      DOUBLE PRECISION, POINTER   :: vertex(:,:)
      CHARACTER(LEN=256) :: machine_string
      CHARACTER(LEN=256) :: date

      INTEGER, PRIVATE                         :: mystart, myend, mydelta, ik_min
      INTEGER, PRIVATE                         :: shar_rank, shar_size
      INTEGER, PRIVATE                         :: win_vertex, win_face

      DOUBLE PRECISION, PRIVATE, POINTER       :: invDenom(:)
                                                  
      DOUBLE PRECISION, PRIVATE, PARAMETER      :: zero = 0.0D+0
      DOUBLE PRECISION, PRIVATE, PARAMETER      :: one  = 1.0D+0
      DOUBLE PRECISION, PRIVATE, PARAMETER      :: epsilon = 1D-6

      LOGICAL, PRIVATE, PARAMETER               :: lverb = .TRUE.

!------------------ Variables for naive approach
      
      INTEGER, POINTER :: ihit_array(:)
      INTEGER, PRIVATE                         :: win_fn, win_a0, win_v0, win_v1, &
                                                  win_dot00, win_dot01, win_dot11, &
                                                  win_d, win_ihit, win_invDenom

      DOUBLE PRECISION, PRIVATE, POINTER       :: FN(:,:), d(:)
      DOUBLE PRECISION, PRIVATE, POINTER       :: A0(:,:), V0(:,:), V1(:,:), V2(:,:),&
                                                  DOT00(:), DOT01(:), DOT02(:),&
                                                  DOT11(:), DOT12(:)

!------------------ Variables for accelerated approach

      TYPE(wall_type), PRIVATE :: wall
      TYPE(block), PRIVATE:: b
      
!-----------------------------------------------------------------------
!     Subroutines
!         wall_load_txt:   Loads triangular mesh from file
!         wall_load_mn:    Creates wall from harmonics
!         wall_load_seg:   Creates wall from segments
!         wall_dump:       Dumps triangulation data
!         wall_info:       Prints wall info
!         collide:         Calculates collision with wall
!                          Has to implementations, for double and float
!         uncount_wall_hit Reduces hit count for last location by one
!         wall_free:       Frees module memory
!-----------------------------------------------------------------------
!     Functions
!         get_wall_ik      Gets index of last hit
!         get_wall_area    Gets area of certain wall index
!-----------------------------------------------------------------------
      INTERFACE collide
         MODULE PROCEDURE collide_double, collide_float
      END INTERFACE

      INTERFACE free_mpi_array
         MODULE PROCEDURE free_mpi_array1d_int, free_mpi_array1d_flt, free_mpi_array1d_dbl, &
                          free_mpi_array2d_int,                       free_mpi_array2d_dbl
      END INTERFACE

      PRIVATE :: mpialloc_1d_int,mpialloc_1d_dbl,mpialloc_2d_int,mpialloc_2d_dbl
      PRIVATE :: free_mpi_array1d_int, free_mpi_array1d_dbl, free_mpi_array2d_int, free_mpi_array2d_dbl
      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Wall Constructors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE INIT_BLOCK(this,xmin,xmax,ymin,ymax,zmin,zmax,nface,istat,comm,shar_comm)
         TYPE(block) :: this
         DOUBLE PRECISION, INTENT(in) :: xmin, ymin, zmin, xmax, ymax, zmax
         INTEGER, INTENT(in) :: nface
         INTEGER, INTENT(inout) :: istat
         INTEGER, INTENT(inout), OPTIONAL :: comm, shar_comm
         INTEGER :: ik, i, dex1, dex2, dex3

         this%xmin = xmin
         this%xmax = xmax
         this%ymin = ymin
         this%ymax = ymax
         this%zmin = zmin
         this%zmax = zmax
         this%nfaces = nface

         IF (lverb) WRITE(6, *) 'Block intialized', xmin, xmax, ymin, ymax, zmin, zmax, nface

         IF (lverb) WRITE(6, *) 'Allocate block info'

         ! allocate memory for information about the mesh
#if defined(MPI_OPT)
         IF (PRESENT(comm)) CALL MPI_BARRIER(comm,istat)
         IF (istat/=0) RETURN
         IF (PRESENT(comm)) THEN
            CALL mpialloc_2d_dbl(this%A0,nface,3,shar_rank,0,shar_comm,this%win_a0)
            CALL mpialloc_2d_dbl(this%V0,nface,3,shar_rank,0,shar_comm,this%win_v0)
            CALL mpialloc_2d_dbl(this%V1,nface,3,shar_rank,0,shar_comm,this%win_v1)
            CALL mpialloc_2d_dbl(this%FN,nface,3,shar_rank,0,shar_comm,this%win_fn)
            this%isshared = .TRUE.
            mydelta = CEILING(REAL(nface) / REAL(shar_size))
            mystart = 1 + shar_rank*mydelta
            myend   = mystart + mydelta
            IF (myend > nface) myend=nface
         ELSE
#endif
            ALLOCATE(this%A0(nface,3),this%V0(nface,3),this%V1(nface,3),&
                  this%FN(nface,3),STAT=istat)
            this%isshared = .FALSE.
            mystart = 1; myend = nface
#if defined(MPI_OPT)
         END IF
#endif
         IF (lverb) WRITE(6, *) 'Calculate mesh in block info'

         IF (istat/=0) RETURN
         ! Precalculate information about mesh
         ! Calculate the face normal
         ! V  = Vertex1-Vertex0
         ! W  = Vertex2-Vertex0
         ! FN = VxW/|VxW|
         ! d  = -Vertex0.FN (. is dot product) (note weve absorbed the negative)
         DO ik = mystart, myend
            dex1 = face(ik,1)
            dex2 = face(ik,2)
            dex3 = face(ik,3)
            this%A0(ik,:) = vertex(dex1,:)
            this%V0(ik,:) = vertex(dex3,:)-vertex(dex1,:)
            this%V1(ik,:) = vertex(dex2,:)-vertex(dex1,:)
            this%FN(ik,1) = (this%V1(ik,2)*this%V0(ik,3))-(this%V1(ik,3)*this%V0(ik,2))
            this%FN(ik,2) = (this%V1(ik,3)*this%V0(ik,1))-(this%V1(ik,1)*this%V0(ik,3))
            this%FN(ik,3) = (this%V1(ik,1)*this%V0(ik,2))-(this%V1(ik,2)*this%V0(ik,1))
         END DO
#if defined(MPI_OPT)
         IF (PRESENT(comm)) CALL MPI_BARRIER(comm,istat)
#endif
         ! Check for zero area
         IF (ANY(SUM(this%FN*this%FN,DIM=2)==zero)) THEN
            istat=-327
            IF (lverb) WRITE(6, *) 'Found zero area in block'
            RETURN
         END IF

         IF (lverb) WRITE(6, *) 'Allocate pre-calculate in block'
         ! allocate memory for information about mesh triangles 
#if defined(MPI_OPT)
         IF (PRESENT(comm)) THEN
            IF (lverb) WRITE(6, *) shar_rank, shar_comm
            CALL mpialloc_1d_dbl(this%DOT00,nface,shar_rank,0,shar_comm,this%win_dot00)
            CALL mpialloc_1d_dbl(this%DOT01,nface,shar_rank,0,shar_comm,this%win_dot01)
            CALL mpialloc_1d_dbl(this%DOT11,nface,shar_rank,0,shar_comm,this%win_dot11)
            IF (lverb) WRITE(6, *) 'Allocate pre-calculate in block 2'
            CALL mpialloc_1d_dbl(this%invDenom,nface,shar_rank,0,shar_comm,this%win_invDenom)
            IF (lverb) WRITE(6, *) 'Allocate pre-calculate in block 3'
            CALL mpialloc_1d_dbl(this%d,nface,shar_rank,0,shar_comm,this%win_d)
            IF (lverb) WRITE(6, *) 'Allocate pre-calculate in block 4'
            CALL mpialloc_1d_int(this%ihit_array,nface,shar_rank,0,shar_comm,this%win_ihit)
            IF (lverb) WRITE(6, *) 'Allocate pre-calculate in block 5'
            mydelta = CEILING(REAL(nface) / REAL(shar_size))
            mystart = 1 + shar_rank*mydelta
            myend   = mystart + mydelta
            IF (myend > nface) myend=nface
         ELSE
#endif
            ALLOCATE(this%DOT00(nface), this%DOT01(nface), this%DOT11(nface), invDenom(nface), STAT=istat)
            ALLOCATE(this%d(nface),STAT=istat)
            ALLOCATE(this%ihit_array(nface),STAT=istat)
            mystart = 1; myend = nface
#if defined(MPI_OPT)
         END IF
#endif
         ! if no error, calculate information about mesh triangles
         IF (lverb)  WRITE(6, *) 'Pre-calculate in block'
         IF (istat/=0) RETURN
         DO i = mystart, myend
            this%ihit_array(i) = 0
            this%DOT00(i) = this%V0(i,1)*this%V0(i,1) + this%V0(i,2)*this%V0(i,2) + this%V0(i,3)*this%V0(i,3)
            this%DOT01(i) = this%V0(i,1)*this%V1(i,1) + this%V0(i,2)*this%V1(i,2) + this%V0(i,3)*this%V1(i,3)
            this%DOT11(i) = this%V1(i,1)*this%V1(i,1) + this%V1(i,2)*this%V1(i,2) + this%V1(i,3)*this%V1(i,3)
            this%d(i)     = this%FN(i,1)*this%A0(i,1) + this%FN(i,2)*this%A0(i,2) + this%FN(i,3)*this%A0(i,3)
            this%invDenom(i) = one / (this%DOT00(i)*this%DOT11(i) - this%DOT01(i)*this%DOT01(i))
         END DO
         ! Multiply with invDenom to reduce calculations later
         DO i = mystart, myend
            this%DOT00(i) = this%DOT00(i) * this%invDenom(i)
            this%DOT01(i) = this%DOT01(i) * this%invDenom(i)
            this%DOT11(i) = this%DOT11(i) * this%invDenom(i)
         END DO
         ! Clear memory
         CALL free_mpi_array(this%win_invDenom, this%invDenom, this%isshared)
      END SUBROUTINE
      
      SUBROUTINE wall_load_txt(filename,istat,comm)
      !-----------------------------------------------------------------------
      ! wall_load_txt: Loads triangular mesh from file
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
      INTEGER, INTENT(inout), OPTIONAL :: comm
      DOUBLE PRECISION :: xmin, ymin, zmin, xmax, ymax, zmax
      INTEGER :: iunit, ik, i, dex1, dex2, dex3
      INTEGER :: shar_comm
      LOGICAL :: shared, first

      first = .true.
      
      shar_rank = 0; shar_size = 1;
      lwall_loaded = .false.
      ! initialize MPI
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_COMM_SPLIT_TYPE(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shar_comm, istat)
         CALL MPI_COMM_RANK( shar_comm, shar_rank, istat )
         CALL MPI_COMM_SIZE( shar_comm, shar_size, istat)
      END IF
#endif
      ! open file, return if fails
      CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
      IF (istat/=0) RETURN
      ! read info
      READ(iunit,'(A)') machine_string
      READ(iunit,'(A)') date
      READ(iunit,*) nvertex,nface
      ! Check if nvertex nface  = 0 -> accelerated structure
      IF (nvertex == 0 .and. nface == 0) THEN
         WRITE(6, *) 'Acc'
         lwall_acc = .true.
         READ(iunit,*) nvertex,nface
#if defined(MPI_OPT)
         IF (PRESENT(comm)) THEN
            CALL mpialloc_2d_dbl(vertex,nvertex,3,shar_rank,0,shar_comm,win_vertex)
            shared = .true.
         ELSE
#endif
            ! if no MPI, allocate everything on one node
            ALLOCATE(vertex(nvertex,3),STAT=istat)
            shared = .false.
#if defined(MPI_OPT)
         END IF
#endif
         IF (lverb) WRITE(6, *) 'Vertex allocation'
         ! read in the mesh on allocated memory
         IF (istat/=0) RETURN
         IF (shar_rank == 0) THEN
            DO ik = 1, nvertex
               READ(iunit,*) vertex(ik,1),vertex(ik,2),vertex(ik,3)
            END DO
         END IF
         IF (lverb) WRITE(6, *) 'Vertex reading'
         ! Read in uniform grid
         READ(iunit, *) wall%nblocks, wall%xstep, wall%ystep, wall%zstep
         ALLOCATE(wall%blocks(wall%nblocks))

         DO ik=1, wall%nblocks
            READ(iunit, *) xmin,xmax,ymin,ymax,zmin,zmax
            READ(iunit, *) nface
            IF (lverb) WRITE(6, *) 'Block reading: ', ik, nface
#if defined(MPI_OPT)
            IF (PRESENT(comm)) THEN
               CALL mpialloc_2d_int(face,nface,3,shar_rank,0,shar_comm,win_face)
               mydelta = CEILING(REAL(nface) / REAL(shar_size))
               mystart = 1 + shar_rank*mydelta
               myend   = mystart + mydelta
               IF (myend > nface) myend=nface
            ELSE
#endif
               ! if no MPI, allocate everything on one node
               ALLOCATE(face(nface,3),STAT=istat)
               mystart = 1; myend=nface
#if defined(MPI_OPT)
            END IF
#endif
            DO i=1, nface
               READ(iunit,*) face(i,1),face(i,2),face(i,3)
            END DO

            IF (lverb) WRITE(6, *) 'Init block: ', ik

            CALL INIT_BLOCK(wall%blocks(ik),xmin,xmax,ymin,ymax,zmin,zmax,nface,istat,comm,shar_comm)

            IF (lverb) WRITE(6, *) 'Init block done: ', ik

            CALL free_mpi_array(win_face, face, shared)
         END DO
         CALL free_mpi_array(win_vertex, vertex, shared)
         ! close file
         CLOSE(iunit)
         ! sync MPI
#if defined(MPI_OPT)
         IF (PRESENT(comm)) THEN
            CALL MPI_BARRIER(shar_comm, istat)
            CALL MPI_COMM_FREE(shar_comm, istat)
         END IF
#endif
      ELSE
         lwall_acc = .false.
         ! allocate shared memory with of mesh
         ! calculate which part each node has to calculate using mydelta
#if defined(MPI_OPT)
         IF (PRESENT(comm)) THEN
            CALL mpialloc_2d_dbl(vertex,nvertex,3,shar_rank,0,shar_comm,win_vertex)
            CALL mpialloc_2d_int(face,nface,3,shar_rank,0,shar_comm,win_face)
            mydelta = CEILING(REAL(nface) / REAL(shar_size))
            mystart = 1 + shar_rank*mydelta
            myend   = mystart + mydelta
            IF (myend > nface) myend=nface
         ELSE
#endif
            ! if no MPI, allocate everything on one node
            ALLOCATE(vertex(nvertex,3),face(nface,3),STAT=istat)
            mystart = 1; myend=nface
#if defined(MPI_OPT)
         END IF
#endif
         ! read in the mesh on allocated memory
         IF (istat/=0) RETURN
         IF (shar_rank == 0) THEN
            DO ik = 1, nvertex
               READ(iunit,*) vertex(ik,1),vertex(ik,2),vertex(ik,3)
            END DO
            DO ik=1,nface
               READ(iunit,*) face(ik,1),face(ik,2),face(ik,3)
            END DO
         END IF
         ! close file
         CLOSE(iunit)
         ! allocate memory for information about the mesh
#if defined(MPI_OPT)
         IF (PRESENT(comm)) CALL MPI_BARRIER(comm,istat)
         IF (istat/=0) RETURN
         IF (PRESENT(comm)) THEN
            CALL mpialloc_2d_dbl(A0,nface,3,shar_rank,0,shar_comm,win_a0)
            CALL mpialloc_2d_dbl(V0,nface,3,shar_rank,0,shar_comm,win_v0)
            CALL mpialloc_2d_dbl(V1,nface,3,shar_rank,0,shar_comm,win_v1)
            CALL mpialloc_2d_dbl(FN,nface,3,shar_rank,0,shar_comm,win_fn)
            mydelta = CEILING(REAL(nface) / REAL(shar_size))
            mystart = 1 + shar_rank*mydelta
            myend   = mystart + mydelta
            IF (myend > nface) myend=nface
         ELSE
#endif
            ALLOCATE(A0(nface,3),V0(nface,3),V1(nface,3),&
                     FN(nface,3),STAT=istat)
            mystart = 1; myend = nface
#if defined(MPI_OPT)
         END IF
#endif
         IF (istat/=0) RETURN
         ! Precalculate information about mesh
         ! Calculate the face normal
         ! V  = Vertex1-Vertex0
         ! W  = Vertex2-Vertex0
         ! FN = VxW/|VxW|
         ! d  = -Vertex0.FN (. is dot product) (note weve absorbed the negative)
         DO ik = mystart, myend
            dex1 = face(ik,1)
            dex2 = face(ik,2)
            dex3 = face(ik,3)
            A0(ik,:) = vertex(dex1,:)
            V0(ik,:)  = vertex(dex3,:)-vertex(dex1,:)
            V1(ik,:)  = vertex(dex2,:)-vertex(dex1,:)
            FN(ik,1) = (V1(ik,2)*V0(ik,3))-(V1(ik,3)*V0(ik,2))
            FN(ik,2) = (V1(ik,3)*V0(ik,1))-(V1(ik,1)*V0(ik,3))
            FN(ik,3) = (V1(ik,1)*V0(ik,2))-(V1(ik,2)*V0(ik,1))
         END DO
#if defined(MPI_OPT)
         IF (PRESENT(comm)) CALL MPI_BARRIER(comm,istat)
#endif
         ! Check for zero area
         IF (ANY(SUM(FN*FN,DIM=2)==zero)) THEN
            istat=-327
            RETURN
         END IF
         ! allocate memory for information about mesh triangles 
#if defined(MPI_OPT)
         IF (PRESENT(comm)) THEN
            CALL mpialloc_1d_dbl(DOT00,nface,shar_rank,0,shar_comm,win_dot00)
            CALL mpialloc_1d_dbl(DOT01,nface,shar_rank,0,shar_comm,win_dot01)
            CALL mpialloc_1d_dbl(DOT11,nface,shar_rank,0,shar_comm,win_dot11)
            CALL mpialloc_1d_dbl(invDenom,nface,shar_rank,0,shar_comm,win_invDenom)
            CALL mpialloc_1d_dbl(d,nface,shar_rank,0,shar_comm,win_d)
            CALL mpialloc_1d_int(ihit_array,nface,shar_rank,0,shar_comm,win_ihit)
            mydelta = CEILING(REAL(nface) / REAL(shar_size))
            mystart = 1 + shar_rank*mydelta
            myend   = mystart + mydelta
            IF (myend > nface) myend=nface
         ELSE
#endif
            ALLOCATE(DOT00(nface), DOT01(nface),&
                     DOT11(nface), invDenom(nface),&
                     STAT=istat)
            ALLOCATE(d(nface),STAT=istat)
            ALLOCATE(ihit_array(nface),STAT=istat)
            mystart = 1; myend = nface
#if defined(MPI_OPT)
         END IF
#endif
         ! if no error, calculate information about mesh triangles
         IF (istat/=0) RETURN
         DO ik = mystart, myend
            ihit_array(ik) = 0
            DOT00(ik) = V0(ik,1)*V0(ik,1) + V0(ik,2)*V0(ik,2) + V0(ik,3)*V0(ik,3)
            DOT01(ik) = V0(ik,1)*V1(ik,1) + V0(ik,2)*V1(ik,2) + V0(ik,3)*V1(ik,3)
            DOT11(ik) = V1(ik,1)*V1(ik,1) + V1(ik,2)*V1(ik,2) + V1(ik,3)*V1(ik,3)
            d(ik)     = FN(ik,1)*A0(ik,1) + FN(ik,2)*A0(ik,2) + FN(ik,3)*A0(ik,3)
            invDenom(ik) = one / (DOT00(ik)*DOT11(ik) - DOT01(ik)*DOT01(ik))
         END DO
         ! Multiply with invDenom to reduce calculations later
         DO ik = mystart, myend
            DOT00(ik) = DOT00(ik) * invDenom(ik)
            DOT01(ik) = DOT01(ik) * invDenom(ik)
            DOT11(ik) = DOT11(ik) * invDenom(ik)
         END DO
         ! sync MPI
#if defined(MPI_OPT)
         IF (PRESENT(comm)) THEN
            CALL MPI_BARRIER(shar_comm, istat)
            CALL MPI_COMM_FREE(shar_comm, istat)
            CALL MPI_WIN_FENCE(0,win_invdenom,istat)
            CALL MPI_WIN_FREE(win_invdenom,istat)
         ELSE
#endif
            IF (ASSOCIATED(invDenom)) DEALLOCATE(invDenom)
#if defined(MPI_OPT)
         END IF
#endif
      END IF
      ! set wall as loaded and return
      lwall_loaded = .true.
      RETURN
      END SUBROUTINE wall_load_txt

      SUBROUTINE wall_load_mn(Rmn,Zmn,xm,xn,mn,nu,nv,comm)
      !-----------------------------------------------------------------------
      ! wall_load_mn: Creates wall from harmonics
      !-----------------------------------------------------------------------
      ! param[in]: Rmn. Harmonics in R direction
      ! param[in]: Zmn. Harmonics in Z direction
      ! param[in]: xm. 
      ! param[in]: xn. 
      ! param[in]: mn. 
      ! param[in]: nu. 
      ! param[in]: nv. 
      ! param[in, out]: comm. MPI communicator, handles shared memory
      !-----------------------------------------------------------------------
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: Rmn(mn), Zmn(mn), xm(mn), xn(mn)
      INTEGER, INTENT(in) :: mn, nu, nv
      INTEGER, INTENT(inout), OPTIONAL :: comm
      INTEGER :: u, v, i, j, istat, dex1, dex2, dex3, ik, nv2
      INTEGER :: shar_comm
      DOUBLE PRECISION :: pi2, th, zt, pi
      DOUBLE PRECISION, ALLOCATABLE :: r_temp(:,:),z_temp(:,:),x_temp(:,:),y_temp(:,:)

      ! create info usually read from file manually
      machine_string = '          HARMONICS'
      date = '      TODAY'
      pi2 = 8.0D+00 * ATAN(one)
      pi  = 4.0E+00 * ATAN(one)
      nv2 = nv/2
      shar_rank = 0; shar_size = 1;
      ! initialize MPI
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_COMM_SPLIT_TYPE(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shar_comm, istat)
         CALL MPI_COMM_RANK( shar_comm, shar_rank, istat )
         CALL MPI_COMM_SIZE( shar_comm, shar_size, istat)
      END IF
#endif
      ! If shared memory rank not changed above, create temporary locations from harmonic calculations
      IF (shar_rank==0)THEN
         ALLOCATE(r_temp(nu,nv),z_temp(nu,nv),x_temp(nu,nv),y_temp(nu,nv))
         r_temp(:,:) = 0
         z_temp(:,:) = 0
         DO u = 1, nu
            DO v = 1, nv
               DO i = 1, mn
                  th = pi2*DBLE(u-1)/DBLE(nu)
                  zt = pi2*DBLE(v-1)/DBLE(nv)
                  r_temp(u,v) = r_temp(u,v) + Rmn(i)*DCOS(xm(i)*th+xn(i)*zt)
                  z_temp(u,v) = z_temp(u,v) + Zmn(i)*DSIN(xm(i)*th+xn(i)*zt)
               END DO
            END DO
         END DO
         DO v = 1, nv
            zt = pi2*DBLE(v-1)/DBLE(nv)
            x_temp(:,v) = r_temp(:,v) * DCOS(zt)
            y_temp(:,v) = r_temp(:,v) * DSIN(zt)
         END DO
      END IF

      nvertex = nu*nv
      nface   = 2*nu*nv
      istat = 0
      ! allocate shared memory with of mesh
      ! calculate which part each node has to calculate using mydelta
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL mpialloc_2d_dbl(vertex,nvertex,3,shar_rank,0,shar_comm,win_vertex)
         CALL mpialloc_2d_int(face,nface,3,shar_rank,0,shar_comm,win_face)
         mydelta = CEILING(REAL(nface) / REAL(shar_size))
         mystart = 1 + shar_rank*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
#endif
         ! if no MPI, allocate everything on one node
         ALLOCATE(vertex(nvertex,3),face(nface,3),STAT=istat)
         mystart = 1; myend=nface
#if defined(MPI_OPT)
      END IF
#endif
      i = 1  ! Tracks vertex index
      j = 1 ! Tracks face index
      ! Do further calculations to create mesh if shared memory rank is zero
      IF (shar_rank==0)THEN
         DO v = 1, nv-1
            DO u = 1, nu-1
               vertex(i,1) = x_temp(u,v)
               vertex(i,2) = y_temp(u,v)
               vertex(i,3) = z_temp(u,v)
               face(j,1) = i
               face(j,2) = i + nu
               face(j,3) = i + 1
               j = j + 1
               face(j,1) = i + nu
               face(j,2) = i + nu + 1
               face(j,3) = i + 1
               j = j + 1
               i=i+1
            END DO
            u = nu
            vertex(i,1) = x_temp(u,v)
            vertex(i,2) = y_temp(u,v)
            vertex(i,3) = z_temp(u,v)
            face(j,1) = i
            face(j,2) = i+nu
            face(j,3) = i-nu+1
            j = j + 1
            face(j,1) = i + nu
            face(j,2) = i + 1
            face(j,3) = i-nu+1
            j = j + 1
            i=i+1
         END DO
         v = nv
         DO u = 1, nu - 1
            vertex(i,1) = x_temp(u,v)
            vertex(i,2) = y_temp(u,v)
            vertex(i,3) = z_temp(u,v)
            face(j,1) = i
            face(j,2) = i - nu*(nv-1)
            face(j,3) = i + 1
            j = j + 1
            face(j,1) = i - nu*(nv-1)
            face(j,2) = i - nu*(nv-1) + 1
            face(j,3) = i + 1
            j = j + 1
            i=i+1
         END DO
         u = nu
         vertex(i,1) = x_temp(u,v)
         vertex(i,2) = y_temp(u,v)
         vertex(i,3) = z_temp(u,v)
         face(j,1) = i
         face(j,2) = nu
         face(j,3) = i - nu + 1
         j = j + 1
         face(j,1) = nu
         face(j,2) = 1
         face(j,3) = i - nu + 1
         j = j + 1
         i=i+1
         ! remove temperary information
         DEALLOCATE(r_temp,z_temp,x_temp,y_temp)
      END IF
      ! if using MPI, wait here
#if defined(MPI_OPT)
      IF (PRESENT(comm)) CALL MPI_BARRIER(shar_comm,istat)
#endif
      ! allocate memory for information about the mesh
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL mpialloc_2d_dbl(A0,nface,3,shar_rank,0,shar_comm,win_a0)
         CALL mpialloc_2d_dbl(V0,nface,3,shar_rank,0,shar_comm,win_v0)
         CALL mpialloc_2d_dbl(V1,nface,3,shar_rank,0,shar_comm,win_v1)
         CALL mpialloc_2d_dbl(FN,nface,3,shar_rank,0,shar_comm,win_fn)
         mydelta = CEILING(REAL(nface) / REAL(shar_size))
         mystart = 1 + shar_rank*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
#endif
         ALLOCATE(A0(nface,3),V0(nface,3),V1(nface,3),&
                  FN(nface,3),STAT=istat)
         mystart = 1; myend = nface
#if defined(MPI_OPT)
      END IF
#endif
      IF (istat/=0) RETURN
      ! Precalculate information about mesh
      ! Calculate the face normal
      ! W  = Vertex1-Vertex0
      ! W  = Vertex2-Vertex0
      ! FN = VxW/|VxW|
      ! d  = -Vertex0.FN (. is dot product) (not weve dropped the minus in our formulation)
      DO ik = mystart, myend
         dex1 = face(ik,1)
         dex2 = face(ik,2)
         dex3 = face(ik,3)
         A0(ik,:) = vertex(dex1,:)
         V0(ik,:)  = vertex(dex3,:)-vertex(dex1,:)
         V1(ik,:)  = vertex(dex2,:)-vertex(dex1,:)
         FN(ik,1) = (V1(ik,2)*V0(ik,3))-(V1(ik,3)*V0(ik,2))
         FN(ik,2) = (V1(ik,3)*V0(ik,1))-(V1(ik,1)*V0(ik,3))
         FN(ik,3) = (V1(ik,1)*V0(ik,2))-(V1(ik,2)*V0(ik,1))
      END DO
      ! allocate memory for information about mesh triangles 
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_BARRIER(shar_comm, istat)
         CALL mpialloc_1d_dbl(DOT00,nface,shar_rank,0,shar_comm,win_dot00)
         CALL mpialloc_1d_dbl(DOT01,nface,shar_rank,0,shar_comm,win_dot01)
         CALL mpialloc_1d_dbl(DOT11,nface,shar_rank,0,shar_comm,win_dot11)
         CALL mpialloc_1d_dbl(invDenom,nface,shar_rank,0,shar_comm,win_invDenom)
         CALL mpialloc_1d_dbl(d,nface,shar_rank,0,shar_comm,win_d)
         CALL mpialloc_1d_int(ihit_array,nface,shar_rank,0,shar_comm,win_ihit)
         mydelta = CEILING(REAL(nface) / REAL(shar_size))
         mystart = 1 + shar_rank*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
#endif
         ALLOCATE(DOT00(nface), DOT01(nface),&
                  DOT11(nface), invDenom(nface),&
                  STAT=istat)
         ALLOCATE(d(nface),STAT=istat)
         ALLOCATE(ihit_array(nface),STAT=istat)
         mystart = 1; myend = nface
#if defined(MPI_OPT)
      END IF
#endif
      ! if no error, calculate information about mesh triangles
      IF (istat/=0) RETURN
      DO ik = mystart, myend
         ihit_array(ik) = 0
         DOT00(ik) = V0(ik,1)*V0(ik,1) + V0(ik,2)*V0(ik,2) + V0(ik,3)*V0(ik,3)
         DOT01(ik) = V0(ik,1)*V1(ik,1) + V0(ik,2)*V1(ik,2) + V0(ik,3)*V1(ik,3)
         DOT11(ik) = V1(ik,1)*V1(ik,1) + V1(ik,2)*V1(ik,2) + V1(ik,3)*V1(ik,3)
         d(ik)     = FN(ik,1)*A0(ik,1) + FN(ik,2)*A0(ik,2) + FN(ik,3)*A0(ik,3)
         invDenom(ik) = one / (DOT00(ik)*DOT11(ik) - DOT01(ik)*DOT01(ik))
      END DO
      ! Multiply with invDenom to reduce calculations later
      DO ik = mystart, myend
         DOT00(ik) = DOT00(ik) * invDenom(ik)
         DOT01(ik) = DOT01(ik) * invDenom(ik)
         DOT11(ik) = DOT11(ik) * invDenom(ik)
      END DO
      ! sync MPI and clear invDenom
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_BARRIER(shar_comm, istat)
         CALL MPI_COMM_FREE(shar_comm, istat)
         CALL MPI_WIN_FENCE(0,win_invdenom,istat)
         CALL MPI_WIN_FREE(win_invdenom,istat)
      ELSE
#endif
         IF (ASSOCIATED(invDenom)) DEALLOCATE(invDenom)
#if defined(MPI_OPT)
      END IF
#endif
      ! set wall as loaded and return
      lwall_loaded = .true.
      RETURN
      END SUBROUTINE wall_load_mn

      SUBROUTINE wall_load_seg(npts,rseg,zseg,nphi,istat,comm)
      !-----------------------------------------------------------------------
      ! wall_load_seg: Creates wall from segments
      !-----------------------------------------------------------------------
      ! param[in]: npts. Number of points in r and Z direction
      ! param[in]: rseg. Segments in r direction (with npts number of points)
      ! param[in]: zseg. Segmetns in z direction (with npts number of points)
      ! param[in]: nphi. Number of points in phi direction
      ! param[in, out]: istat. Integer that shows error if != 0
      ! param[in, out]: comm. MPI communicator, handles shared memory
      !-----------------------------------------------------------------------
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      INTEGER, INTENT(in) :: npts, nphi
      DOUBLE PRECISION, INTENT(in) :: rseg(npts)
      DOUBLE PRECISION, INTENT(in) :: Zseg(npts)
      INTEGER, INTENT(inout)       :: istat
      INTEGER, INTENT(inout), OPTIONAL :: comm
      INTEGER :: shar_comm
      INTEGER :: nseg, ij, ik, il, im
      DOUBLE PRECISION :: dphi
      INTEGER :: dex1, dex2, dex3
      
      shar_rank = 0; shar_size = 1;
      dphi = 8.0D+00 * ATAN(one)/nphi
      ! initialize MPI
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_COMM_SPLIT_TYPE(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shar_comm, istat)
         CALL MPI_COMM_RANK( shar_comm, shar_rank, istat )
         CALL MPI_COMM_SIZE( shar_comm, shar_size, istat)
      END IF
#endif
      ! create info usually read from file manually
      machine_string = '          SEGMENTS'
      date = '      TODAY'
      nseg = npts-1
      nvertex = nseg * 4 * nphi
      nface   = nseg * 2 * nphi
      ! allocate shared memory with of mesh
      ! calculate which part each node has to calculate using mydelta
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL mpialloc_2d_dbl(vertex,nvertex,3,shar_rank,0,shar_comm,win_vertex)
         CALL mpialloc_2d_int(face,nface,3,shar_rank,0,shar_comm,win_face)
         mydelta = CEILING(REAL(nface) / REAL(shar_size))
         mystart = 1 + shar_rank*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
#endif
         ALLOCATE(vertex(nvertex,3),face(nface,3),STAT=istat)
         mystart = 1; myend=nface
#if defined(MPI_OPT)
      END IF
#endif
      ! if no error, and shared rank is zero, create triangles from input info
      IF (istat/=0) RETURN
      IF (shar_rank == 0) THEN
         il = 1; im = 1
         DO ik = 1, nphi
            DO ij = 1, nseg
               ! create the 4 points
               vertex(il,1) = rseg(ij)*cos(dphi*(ik-1))
               vertex(il,2) = rseg(ij)*sin(dphi*(ik-1))
               vertex(il,3) = zseg(ij)
               vertex(il+1,1) = rseg(ij+1)*cos(dphi*(ik-1))
               vertex(il+1,2) = rseg(ij+1)*sin(dphi*(ik-1))
               vertex(il+1,3) = zseg(ij+1)
               vertex(il+2,1) = rseg(ij)*cos(dphi*ik)
               vertex(il+2,2) = rseg(ij)*sin(dphi*ik)
               vertex(il+2,3) = zseg(ij)
               vertex(il+3,1) = rseg(ij+1)*cos(dphi*ik)
               vertex(il+3,2) = rseg(ij+1)*sin(dphi*ik)
               vertex(il+3,3) = zseg(ij+1)
               ! Create the 2 triangles
               face(im,1) = il
               face(im,2) = il+2
               face(im,3) = il+1
               face(im+1,1) = il+2
               face(im+1,2) = il+3
               face(im+1,3) = il+1
               ! Adjust index
               im = im + 2
               il = il + 4
            END DO 
         END DO
      END IF
      ! if using MPI, wait here
      ! allocate memory for information about the mesh
#if defined(MPI_OPT)
      IF (PRESENT(comm)) CALL MPI_BARRIER(comm,istat)
      IF (istat/=0) RETURN
      IF (PRESENT(comm)) THEN
         CALL mpialloc_2d_dbl(A0,nface,3,shar_rank,0,shar_comm,win_a0)
         CALL mpialloc_2d_dbl(V0,nface,3,shar_rank,0,shar_comm,win_v0)
         CALL mpialloc_2d_dbl(V1,nface,3,shar_rank,0,shar_comm,win_v1)
         CALL mpialloc_2d_dbl(FN,nface,3,shar_rank,0,shar_comm,win_fn)
         mydelta = CEILING(REAL(nface) / REAL(shar_size))
         mystart = 1 + shar_rank*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
#endif
         ALLOCATE(A0(nface,3),V0(nface,3),V1(nface,3),&
                  FN(nface,3),STAT=istat)
         mystart = 1; myend = nface
#if defined(MPI_OPT)
      END IF
#endif
      IF (istat/=0) RETURN
      ! Precalculate information about mesh
      ! Calculate the face normal
      ! V  = Vertex1-Vertex0
      ! W  = Vertex2-Vertex0
      ! FN = VxW/|VxW|
      ! d  = -Vertex0.FN (. is dot product) (note weve absorbed the negative)
      DO ik = mystart, myend
         dex1 = face(ik,1)
         dex2 = face(ik,2)
         dex3 = face(ik,3)
         A0(ik,:) = vertex(dex1,:)
         V0(ik,:)  = vertex(dex3,:)-vertex(dex1,:)
         V1(ik,:)  = vertex(dex2,:)-vertex(dex1,:)
         FN(ik,1) = (V1(ik,2)*V0(ik,3))-(V1(ik,3)*V0(ik,2))
         FN(ik,2) = (V1(ik,3)*V0(ik,1))-(V1(ik,1)*V0(ik,3))
         FN(ik,3) = (V1(ik,1)*V0(ik,2))-(V1(ik,2)*V0(ik,1))
      END DO
#if defined(MPI_OPT)
      IF (PRESENT(comm)) CALL MPI_BARRIER(comm,istat)
#endif
      ! Check for zero area
      IF (ANY(SUM(FN*FN,DIM=2)==zero)) THEN
         istat=-327
         RETURN
      END IF
      ! allocate memory for information about mesh triangles 
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL mpialloc_1d_dbl(DOT00,nface,shar_rank,0,shar_comm,win_dot00)
         CALL mpialloc_1d_dbl(DOT01,nface,shar_rank,0,shar_comm,win_dot01)
         CALL mpialloc_1d_dbl(DOT11,nface,shar_rank,0,shar_comm,win_dot11)
         CALL mpialloc_1d_dbl(invDenom,nface,shar_rank,0,shar_comm,win_invDenom)
         CALL mpialloc_1d_dbl(d,nface,shar_rank,0,shar_comm,win_d)
         CALL mpialloc_1d_int(ihit_array,nface,shar_rank,0,shar_comm,win_ihit)
         mydelta = CEILING(REAL(nface) / REAL(shar_size))
         mystart = 1 + shar_rank*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
#endif
         ALLOCATE(DOT00(nface), DOT01(nface),&
                  DOT11(nface), invDenom(nface),&
                  STAT=istat)
         ALLOCATE(d(nface),STAT=istat)
         ALLOCATE(ihit_array(nface),STAT=istat)
         mystart = 1; myend = nface
#if defined(MPI_OPT)
      END IF
#endif
      ! if no error, calculate information about mesh triangles
      IF (istat/=0) RETURN
      DO ik = mystart, myend
         ihit_array(ik) = 0
         DOT00(ik) = V0(ik,1)*V0(ik,1) + V0(ik,2)*V0(ik,2) + V0(ik,3)*V0(ik,3)
         DOT01(ik) = V0(ik,1)*V1(ik,1) + V0(ik,2)*V1(ik,2) + V0(ik,3)*V1(ik,3)
         DOT11(ik) = V1(ik,1)*V1(ik,1) + V1(ik,2)*V1(ik,2) + V1(ik,3)*V1(ik,3)
         d(ik)     = FN(ik,1)*A0(ik,1) + FN(ik,2)*A0(ik,2) + FN(ik,3)*A0(ik,3)
         invDenom(ik) = one / (DOT00(ik)*DOT11(ik) - DOT01(ik)*DOT01(ik))
      END DO
      ! Multiply with invDenom to reduce calculations later
      DO ik = mystart, myend
         DOT00(ik) = DOT00(ik) * invDenom(ik)
         DOT01(ik) = DOT01(ik) * invDenom(ik)
         DOT11(ik) = DOT11(ik) * invDenom(ik)
      END DO
      ! sync MPI
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_BARRIER(shar_comm, istat)
         CALL MPI_COMM_FREE(shar_comm, istat)
         CALL MPI_WIN_FENCE(0,win_invdenom,istat)
         CALL MPI_WIN_FREE(win_invdenom,istat)
      ELSE
#endif
         IF (ASSOCIATED(invDenom)) DEALLOCATE(invDenom)
#if defined(MPI_OPT)
      END IF
#endif
      ! set wall as loaded and return
      lwall_loaded = .true.
      RETURN
      END SUBROUTINE wall_load_seg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Wall info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE wall_dump(filename,istat)
      !-----------------------------------------------------------------------
      ! wall_dump: Dumps triangulation data
      !-----------------------------------------------------------------------
      ! param[in]: filename. The file name to load in
      ! param[in, out]: istat. Integer that shows error if != 0
      !-----------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER, INTENT(inout)       :: istat
      INTEGER :: iunit, ik
      ! open file
      CALL safe_open(iunit,istat,'wall_dump.'//TRIM(filename),'unknown','formatted')
      ! for every face, output info
      DO ik = 1, nface
         WRITE(iunit,'(13(ES20.10))')  A0(ik,1), A0(ik,2), A0(ik,3), &
                                       FN(ik,1), FN(ik,2), FN(ik,3),&
                                       V0(ik,1), V0(ik,2), V0(ik,3),&
                                       V1(ik,1), V1(ik,2), V1(ik,3), &
                                       d(ik)
      END DO
      WRITE(iunit,*) '#  V0(ik,1), V0(ik,2), V0(ik,3), FN(ik,1), FN(ik,2), FN(ik,3),',&
                        'V(ik,1), V(ik,2), V(ik,3), W(ik,1), W(ik,2), W(ik,3), d(ik)'
      ! close file
      CLOSE(iunit)
      RETURN
      END SUBROUTINE wall_dump

      
      SUBROUTINE wall_info(iunit)
      !-----------------------------------------------------------------------
      ! wall_info: Prints wall info
      !-----------------------------------------------------------------------
      ! param[in]: iunit. Location to print information about wall
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in)    :: iunit
      INTEGER                :: i
      WRITE(iunit,'(A)')         ' -----  Vessel Information  -----'
      WRITE(iunit,'(3X,A,A)')    'Wall Name : ',TRIM(machine_string(10:))
      WRITE(iunit,'(3X,A,A)')    'Date      : ',TRIM(date(6:))
      IF (lwall_acc) THEN
         WRITE(iunit,'(3X,A,I7)')    'Blocks    : ',wall%nblocks
         DO i=1,wall%nblocks
            WRITE(iunit,'(3X,A,I4,A,I7)')    'Block ', i, '  : ', wall%blocks(i)%nfaces
         END DO
         
      ELSE
         WRITE(iunit,'(3X,A,I7)')   'Faces     : ',nface
      END IF      
      RETURN
      END SUBROUTINE wall_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Wall collide
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE collide_float(x0,y0,z0,x1,y1,z1,xw,yw,zw,lhit)
      !-----------------------------------------------------------------------
      ! collide_float: Implementation of collide for floating point values
      !-----------------------------------------------------------------------
      ! param[in]: x0. x-location of first point of the line segment to check
      ! param[in]: y0. y-location of first point of the line segment to check
      ! param[in]: z0. z-location of first point of the line segment to check
      ! param[in]: x1. x-location of second point of the line segment to check
      ! param[in]: y1. y-location of second point of the line segment to check
      ! param[in]: z1. z-location of second point of the line segment to check
      ! param[out]: xw. x-location of hit (if hit has been found)
      ! param[out]: yw. y-location of hit (if hit has been found)
      ! param[out]: zw. z-location of hit (if hit has been found)
      ! param[out]: lhit. Logical that shows if hit has been found or not
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL, INTENT(in) :: x0, y0, z0, x1, y1, z1
      REAL, INTENT(out) :: xw, yw, zw
      LOGICAL, INTENT(out) :: lhit
      DOUBLE PRECISION :: x0d, y0d, z0d, x1d, y1d, z1d
      DOUBLE PRECISION :: xwd, ywd, zwd
      LOGICAL          :: lhit2
      ! function simply converts from floating point to double to help compiler
      xw=zero; yw=zero; zw=zero; lhit=.FALSE.
      x0d=x0; y0d=y0; z0d=z0
      x1d=x1; y1d=y1; z1d=z1
      CALL collide_double(x0d,y0d,z0d,x1d,y1d,z1d,xwd,ywd,zwd,lhit2)
      xw=xwd; yw=ywd; zw=zwd; lhit=lhit2
      RETURN
      END SUBROUTINE collide_float

      SUBROUTINE collide_double(x0,y0,z0,x1,y1,z1,xw,yw,zw,lhit)
      !-----------------------------------------------------------------------
      ! collide_double: Implementation of collide for double precision values
      !-----------------------------------------------------------------------
      ! param[in]: x0. x-location of first point of the line segment to check
      ! param[in]: y0. y-location of first point of the line segment to check
      ! param[in]: z0. z-location of first point of the line segment to check
      ! param[in]: x1. x-location of second point of the line segment to check
      ! param[in]: y1. y-location of second point of the line segment to check
      ! param[in]: z1. z-location of second point of the line segment to check
      ! param[out]: xw. x-location of hit (if hit has been found)
      ! param[out]: yw. y-location of hit (if hit has been found)
      ! param[out]: zw. z-location of hit (if hit has been found)
      ! param[out]: lhit. Logical that shows if hit has been found or not
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x0, y0, z0, x1, y1, z1
      DOUBLE PRECISION, INTENT(out) :: xw, yw, zw
      LOGICAL, INTENT(out) :: lhit
      INTEGER :: ik, k1,k2, b_index
      DOUBLE PRECISION :: drx, dry, drz, V2x, V2y, V2z, DOT02l, DOT12l, tloc, tmin, alphal, betal
      DOUBLE PRECISION :: tDeltaXfwd, tDeltaXbck, tDeltaYfwd, tDeltaYbck, tDeltaZfwd, tDeltaZbck
      DOUBLE PRECISION :: DeltaXfwd, DeltaXbck, DeltaYfwd, DeltaYbck, DeltaZfwd, DeltaZbck
      DOUBLE PRECISION :: tcompx, tcompy, tcompz
      LOGICAL :: xfwd, xbck, yfwd, ybck, zfwd, zbck
      xw=zero; yw=zero; zw=zero; lhit=.FALSE.
      ik_min = zero
      tmin = one + epsilon
      ! Define DR
      drx = x1-x0
      dry = y1-y0
      drz = z1-z0
      IF (lwall_acc) THEN
         xfwd = .false.
         xbck = .false.
         yfwd = .false.
         ybck = .false.
         zfwd = .false.
         zbck = .false.

         tcompx = zero
         tcompy = zero
         tcompz = zero
         !WRITE(6, *) x0, y0, z0, x1, y1, z1
         k1 = 1; k2 = wall%nblocks
         DO b_index = k1,k2
            b = wall%blocks(b_index)
            IF (x0 <= b%xmax .and. x0 >= b%xmin .and. y0 <= b%ymax .and. y0 >= b%ymin .and. z0 <= b%zmax .and. z0 >= b%zmin) THEN
               !WRITE(6, *) 'True for: ', b_index, b%xmax, b%xmin, b%ymax, b%ymin, b%zmax, b%zmin
               EXIT
            END IF
         END DO

         DO WHILE (.true.)
            xw=zero; yw=zero; zw=zero; lhit=.FALSE.
            ik_min = zero
            tmin = one + epsilon
            !WRITE(6, *) b_index, wall%nblocks
            IF (b_index > wall%nblocks .or. b_index < 1) EXIT
            b = wall%blocks(b_index)

            k1 = 1; k2 = b%nfaces
            ! Check every triangle
            ! Based on Badouel's algorithm
            ! Source: https://graphics.stanford.edu/courses/cs348b-98/gg/intersect.html
            DO ik = k1,k2
               ! calculate whether or not this line segment ever hits the plane of the triangle
               alphal = b%FN(ik,1)*drx + b%FN(ik,2)*dry + b%FN(ik,3)*drz
               betal = b%FN(ik,1)*x0 + b%FN(ik,2)*y0 + b%FN(ik,3)*z0
               ! tloc indicated when hit. If hit between r0 and r1, tloc between 0 and 1
               tloc = (b%d(ik)-betal)/alphal
               IF (tloc > one) CYCLE
               IF (tloc <= zero) CYCLE
               ! If the line segment hits the plane of the triangle
               ! calculate if it actually hits on the triangle
               V2x = x0 + tloc*drx - b%A0(ik,1)
               V2y = y0 + tloc*dry - b%A0(ik,2)
               V2z = z0 + tloc*drz - b%A0(ik,3)
               DOT02l = b%V0(ik,1)*V2x + b%V0(ik,2)*V2y + b%V0(ik,3)*V2z
               DOT12l = b%V1(ik,1)*V2x + b%V1(ik,2)*V2y + b%V1(ik,3)*V2z
               alphal = b%DOT11(ik)*DOT02l-b%DOT01(ik)*DOT12l
               betal  = b%DOT00(ik)*DOT12l-b%DOT01(ik)*DOT02l
               ! In that case, these should be false
               IF ((alphal < -epsilon) .or. (betal < -epsilon) .or. (alphal+betal > one + epsilon)) CYCLE
               ! else check if this was the closest hit, and then store
               IF (tloc < tmin) THEN
                  ik_min = ik
                  tmin = tloc
               END IF
            END DO

            DeltaXfwd = (b%xmax - x0)
            DeltaXbck = (b%xmin - x0)
            DeltaYfwd = (b%ymax - y0)
            DeltaYbck = (b%ymin - y0)
            DeltaZfwd = (b%zmax - z0)
            DeltaZbck = (b%zmin - z0)

            tDeltaXfwd = DeltaXfwd / drx
            tDeltaXbck = DeltaXbck / drx
            tDeltaYfwd = DeltaYfwd / dry
            tDeltaYbck = DeltaYbck / dry
            tDeltaZfwd = DeltaZfwd / drz
            tDeltaZbck = DeltaZbck / drz

            !WRITE(6, *) tmin, drx, dry, drz
            !WRITE(6, *) DeltaXfwd, DeltaXbck, DeltaYfwd, DeltaYbck, DeltaZfwd, DeltaZbck
            !WRITE(6, *) tDeltaXfwd, tDeltaXbck, tDeltaYfwd, tDeltaYbck, tDeltaZfwd, tDeltaZbck

            ! check if leaves block before hits. 
            ! Compare with tcomp to make sure that ray always moves in positive times 
            IF (tDeltaXfwd < tmin .and. tDeltaXfwd > tcompx) THEN
               !WRITE(6, *) 'Doing forward x-step'
               b_index = b_index + wall%xstep
               tcompx = tDeltaXfwd + epsilon
            ELSE IF (tDeltaXbck < tmin .and. tDeltaXbck > tcompx) THEN
               !WRITE(6, *) 'Doing backward x-step'
               b_index = b_index - wall%xstep
               tcompx = tDeltaXbck + epsilon
            ELSE IF (tDeltaYfwd < tmin .and. tDeltaYfwd > tcompy) THEN
               !WRITE(6, *) 'Doing forward y-step'
               b_index = b_index + wall%ystep
               tcompy = tDeltaYfwd + epsilon
            ELSE IF (tDeltaYbck < tmin .and. tDeltaYbck > tcompy) THEN
               !WRITE(6, *) 'Doing backward y-step'
               b_index = b_index - wall%ystep
               tcompy = tDeltaYbck + epsilon
            ELSE IF (tDeltaZfwd < tmin .and. tDeltaZfwd > tcompz) THEN
               !WRITE(6, *) 'Doing forward z-step'
               b_index = b_index + wall%zstep
               tcompz = tDeltaZfwd + epsilon
            ELSE IF (tDeltaZbck < tmin .and. tDeltaZbck > tcompz) THEN
               !WRITE(6, *) 'Doing backward z-step'
               b_index = b_index - wall%zstep
               tcompz = tDeltaZbck + epsilon
            ELSE
               ! If no collision found, check numerical margins
               IF (tmin < one) THEN
                  EXIT
               ELSE
                  IF (DeltaXfwd < epsilon .and. DeltaXfwd > -epsilon .and. .not. xbck) THEN
                     !WRITE(6, *) 'Doing forward x-step due to numerical'
                     b_index = b_index + wall%xstep
                     xfwd = .true.
                  ELSE IF (DeltaXbck < epsilon .and. DeltaXbck > -epsilon .and. .not. xfwd) THEN
                     !WRITE(6, *) 'Doing backward x-step due to numerical'
                     b_index = b_index - wall%xstep
                     xbck = .true.
                  ELSE IF (DeltaYfwd < epsilon .and. DeltaYfwd > -epsilon .and. .not. ybck) THEN
                     !WRITE(6, *) 'Doing forward y-step due to numerical'
                     b_index = b_index + wall%ystep
                     yfwd = .true.
                  ELSE IF (DeltaYbck < epsilon .and. DeltaYbck > -epsilon .and. .not. yfwd) THEN
                     !WRITE(6, *) 'Doing backward y-step due to numerical'
                     b_index = b_index - wall%ystep
                     ybck = .true.
                  ELSE IF (DeltaZfwd < epsilon .and. DeltaZfwd > -epsilon .and. .not. zbck) THEN
                     !WRITE(6, *) 'Doing forward z-step due to numerical'
                     b_index = b_index + wall%zstep
                     zfwd = .true.
                  ELSE IF (DeltaZbck < epsilon .and. DeltaZbck > -epsilon .and. .not. zfwd) THEN
                     !WRITE(6, *) 'Doing backward z-step due to numerical'
                     b_index = b_index - wall%zstep
                     zbck = .true.
                  ELSE
                     EXIT
                  END IF
               END IF 
            END IF
         END DO

         ! if any index stored, hit was found, calculate location and increment ihit_array
         IF (ik_min > zero) THEN
            lhit = .TRUE.
            xw   = x0 + tmin*drx
            yw   = y0 + tmin*dry
            zw   = z0 + tmin*drz
            b%ihit_array(ik_min) = b%ihit_array(ik_min) + 1
         END IF


      ELSE
         k1 = 1; k2 = nface
         
         ! Check every triangle
         ! Based on Badouel's algorithm
         ! Source: https://graphics.stanford.edu/courses/cs348b-98/gg/intersect.html
         DO ik = k1,k2
            ! calculate whether or not this line segment ever hits the plane of the triangle
            alphal = FN(ik,1)*drx + FN(ik,2)*dry + FN(ik,3)*drz
            betal = FN(ik,1)*x0 + FN(ik,2)*y0 + FN(ik,3)*z0
            ! tloc indicated when hit. If hit between r0 and r1, tloc between 0 and 1
            tloc = (d(ik)-betal)/alphal
            IF (tloc > one) CYCLE
            IF (tloc <= zero) CYCLE
            ! If the line segment hits the plane of the triangle
            ! calculate if it actually hits on the triangle
            V2x = x0 + tloc*drx - A0(ik,1)
            V2y = y0 + tloc*dry - A0(ik,2)
            V2z = z0 + tloc*drz - A0(ik,3)
            DOT02l = V0(ik,1)*V2x + V0(ik,2)*V2y + V0(ik,3)*V2z
            DOT12l = V1(ik,1)*V2x + V1(ik,2)*V2y + V1(ik,3)*V2z
            alphal = DOT11(ik)*DOT02l-DOT01(ik)*DOT12l
            betal  = DOT00(ik)*DOT12l-DOT01(ik)*DOT02l
            ! In that case, these should be false
            IF ((alphal < -epsilon) .or. (betal < -epsilon) .or. (alphal+betal > one + epsilon)) CYCLE
            ! else check if this was the closest hit, and then store
            IF (tloc < tmin) THEN
               ik_min = ik
               tmin = tloc
            END IF
         END DO
         ! if any index stored, hit was found, calculate location and increment ihit_array
         IF (ik_min > zero) THEN
            lhit = .TRUE.
            xw   = x0 + tmin*drx
            yw   = y0 + tmin*dry
            zw   = z0 + tmin*drz
            ihit_array(ik_min) = ihit_array(ik_min) + 1
         END IF
      END IF
      RETURN
      END SUBROUTINE collide_double

      SUBROUTINE uncount_wall_hit
      !-----------------------------------------------------------------------
      ! uncount_wall_hit: Reduces ihit_array at last found hit location with one
      !-----------------------------------------------------------------------
         IMPLICIT NONE
         if (lwall_acc) THEN
            b%ihit_array(ik_min) = b%ihit_array(ik_min) - 1
         ELSE
            ihit_array(ik_min) = ihit_array(ik_min) - 1
         END IF
      END SUBROUTINE

      INTEGER FUNCTION get_wall_ik()
      !-----------------------------------------------------------------------
      ! get_wall_ik: Gets index of last hit location
      !-----------------------------------------------------------------------
      ! return[integer]: get_wall_ik. Last hit index 
      !-----------------------------------------------------------------------
         IMPLICIT NONE
         get_wall_ik = ik_min
         RETURN
      END FUNCTION

      DOUBLE PRECISION FUNCTION get_wall_area(ik)
      !-----------------------------------------------------------------------
      ! get_wall_ik: Gets index of last hit location
      !-----------------------------------------------------------------------
      ! param[in]: ik. Index of wall location to check
      ! return[double]: get_wall_area. Area of wall location checked 
      !-----------------------------------------------------------------------
         IMPLICIT NONE
         INTEGER, INTENT(in) :: ik
         IF (lwall_acc) THEN
            get_wall_area = 0.5*SQRT(SUM(b%FN(ik,:)*b%FN(ik,:)))
         ELSE
            get_wall_area = 0.5*SQRT(SUM(FN(ik,:)*FN(ik,:)))
         END IF
         RETURN
      END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Wall Destructors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE BLOCK_DESTROY(this)
         IMPLICIT NONE
         TYPE(block), INTENT(inout) :: this
         this % nfaces = 0
         CALL free_mpi_array(this % win_d, this % d, this % isshared)
         CALL free_mpi_array(this % win_FN, this % FN, this % isshared)
         CALL free_mpi_array(this % win_A0, this % A0, this % isshared)
         CALL free_mpi_array(this % win_V0, this % V0, this % isshared)
         CALL free_mpi_array(this % win_V1, this % V1, this % isshared)
         CALL free_mpi_array(this % win_DOT00, this % DOT00, this % isshared)
         CALL free_mpi_array(this % win_DOT01, this % DOT01, this % isshared)
         CALL free_mpi_array(this % win_DOT11, this % DOT11, this % isshared)
         CALL free_mpi_array(this % win_ihit, this % ihit_array, this % isshared)
         this % isshared = .FALSE.
         RETURN
         END SUBROUTINE BLOCK_DESTROY

      SUBROUTINE wall_free(istat,shared_comm)
      !-----------------------------------------------------------------------
      ! wall_free: Removes wall from memory
      !-----------------------------------------------------------------------
      ! param[in, out]: istat. Integer that shows error if != 0
      ! param[in, out]: shared_comm. MPI communicator, handles shared memory
      !-----------------------------------------------------------------------
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: istat
      INTEGER, INTENT(inout), OPTIONAL :: shared_comm
      INTEGER :: i
      ! Check if have to delete acceerated wall or the normal one
      IF (lwall_acc) THEN
         IF (ASSOCIATED(wall%blocks)) THEN
            DO i = 1, wall%nblocks
               CALL BLOCK_DESTROY(wall%blocks(i))
            END DO
            DEALLOCATE(wall%blocks)
            wall%nblocks = -1
         END IF
      ELSE
         IF (PRESENT(shared_comm)) THEN
#if defined(MPI_OPT)
            CALL MPI_WIN_FENCE(0,win_vertex,istat)
            CALL MPI_WIN_FREE(win_vertex,istat)
            CALL MPI_WIN_FENCE(0,win_face,istat)
            CALL MPI_WIN_FREE(win_face,istat)
            CALL MPI_WIN_FENCE(0,win_fn,istat)
            CALL MPI_WIN_FREE(win_fn,istat)
            CALL MPI_WIN_FENCE(0,win_a0,istat)
            CALL MPI_WIN_FREE(win_a0,istat)
            CALL MPI_WIN_FENCE(0,win_v0,istat)
            CALL MPI_WIN_FREE(win_v0,istat)
            CALL MPI_WIN_FENCE(0,win_v1,istat)
            CALL MPI_WIN_FREE(win_v1,istat)
            CALL MPI_WIN_FENCE(0,win_dot00,istat)
            CALL MPI_WIN_FREE(win_dot00,istat)
            CALL MPI_WIN_FENCE(0,win_dot01,istat)
            CALL MPI_WIN_FREE(win_dot01,istat)
            CALL MPI_WIN_FENCE(0,win_dot11,istat)
            CALL MPI_WIN_FREE(win_dot11,istat)
            CALL MPI_WIN_FENCE(0,win_d,istat)
            CALL MPI_WIN_FREE(win_d,istat)
            CALL MPI_WIN_FENCE(0,win_ihit,istat)
            CALL MPI_WIN_FREE(win_ihit,istat)    
            IF (ASSOCIATED(vertex)) NULLIFY(vertex)
            IF (ASSOCIATED(face)) NULLIFY(face)
            IF (ASSOCIATED(FN)) NULLIFY(FN)
            IF (ASSOCIATED(A0)) NULLIFY(A0)
            IF (ASSOCIATED(V0)) NULLIFY(V0)
            IF (ASSOCIATED(V1)) NULLIFY(V1)
            IF (ASSOCIATED(DOT00)) NULLIFY(DOT00)
            IF (ASSOCIATED(DOT01)) NULLIFY(DOT01)
            IF (ASSOCIATED(DOT11)) NULLIFY(DOT11)
            IF (ASSOCIATED(d)) NULLIFY(d)
            IF (ASSOCIATED(ihit_array)) NULLIFY(ihit_array)
         ELSE
#endif
            IF (ASSOCIATED(FN)) DEALLOCATE(FN)
            IF (ASSOCIATED(A0)) DEALLOCATE(A0)
            IF (ASSOCIATED(V0)) DEALLOCATE(V0)
            IF (ASSOCIATED(V1)) DEALLOCATE(V1)
            IF (ASSOCIATED(DOT00)) DEALLOCATE(DOT00)
            IF (ASSOCIATED(DOT01)) DEALLOCATE(DOT01)
            IF (ASSOCIATED(DOT11)) DEALLOCATE(DOT11)
            IF (ASSOCIATED(d)) DEALLOCATE(d)
            IF (ASSOCIATED(vertex)) DEALLOCATE(vertex)
            IF (ASSOCIATED(face)) DEALLOCATE(face)
            IF (ASSOCIATED(ihit_array)) DEALLOCATE(ihit_array)
         END IF
      END IF
      machine_string=''
      date=''
      nface = -1
      nvertex = -1

      lwall_loaded = .false.
      RETURN
      END SUBROUTINE wall_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Memory Allocation Subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE mpialloc_1d_int(array,n1,subid,mymaster,share_comm,win)
      !-----------------------------------------------------------------------
      ! mpialloc_1d_int: Allocated a 1D integer array to shared memory
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
      INTEGER, POINTER, INTENT(inout) :: array(:)
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
      IF (subid == mymaster) window_size = INT(n1,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
      IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
      CALL C_F_POINTER(baseptr, array, array_shape)
#endif
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
      !-----------------------------------------------------------------------
      ! mpialloc_1d_int: Allocated a 2D integer array to shared memory
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
      INTEGER, POINTER, INTENT(inout) :: array(:,:)
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
      IF (subid == mymaster) window_size = INT(n1*n2,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
      IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
      CALL C_F_POINTER(baseptr, array, array_shape)
#endif
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
   
         SUBROUTINE free_mpi_array1d_flt(win_local,array_local,isshared)
         IMPLICIT NONE
         LOGICAL, INTENT(in) :: isshared
         INTEGER, INTENT(inout) :: win_local
         REAL, POINTER, INTENT(inout) :: array_local(:)
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
         END SUBROUTINE free_mpi_array1d_flt
   
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

!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE wall_mod
