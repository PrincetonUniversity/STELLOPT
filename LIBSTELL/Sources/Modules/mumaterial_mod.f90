!-----------------------------------------------------------------------
!     Module:        mumaterial_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          January 2022
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
!         
!-----------------------------------------------------------------------
      INTEGER, PRIVATE            :: nvoxels


      DOUBLE PRECISION, POINTER, PRIVATE :: volume(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: chi(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: mx(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: my(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: mz(:)
      DOUBLE PRECISION, POINTER, PRIVATE :: vertex(:,:)
      DOUBLE PRECISION, POINTER, PRIVATE :: dx(:,:)
      DOUBLE PRECISION, POINTER, PRIVATE :: dy(:,:)
      DOUBLE PRECISION, POINTER, PRIVATE :: dz(:,:)
      DOUBLE PRECISION, POINTER, PRIVATE :: ndotm(:,:)
      INTEGER, PRIVATE            :: win_volume, win_chi,  &
                                     win_dx, win_dy, win_dz, &
                                     win_mx, win_my, win_mz, &
                                     win_vertex,win_ndotm 


      INTEGER, PRIVATE                    :: mystart, myend, mydelta
      INTEGER, PRIVATE                    :: shar_rank, shar_size, shar_comm

      
!-----------------------------------------------------------------------
!     Subroutines
!         mumaterial_load: Loads a magnetic material file
!         mumaterial_init: Calculates mangetic moment of material
!         mumaterial_getb: Calculates magnetic field at a point in space
!-----------------------------------------------------------------------
!     Functions
!-----------------------------------------------------------------------
      CONTAINS
      
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
      CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
      IF (istat/=0) RETURN
      ! read info
      IF (shar_rank == 0) THEN
         READ(iunit,*) nvoxels
      END IF
      ! Broadcast info to MPI and allocate vertex and face info
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_Bcast(nvoxels,1,MPI_INTEGER,0,shar_comm,istat)
         CALL mpialloc_2d_dbl(vertex,nvoxels,3,shar_rank,0,shar_comm,win_vertex)
         CALL mpialloc_1d_dbl(volume,nvoxels,shar_rank,0,shar_comm,win_volume)
         CALL mpialloc_1d_dbl(chi,nvoxels,shar_rank,0,shar_comm,win_chi)
         shared = .true.
      ELSE
#endif
         ! if no MPI, allocate everything on one node
         ALLOCATE(volume(nvoxels),chi(nvoxels),vertex(nvoxels,3),STAT=istat)
         shared = .false.
#if defined(MPI_OPT)
      END IF
#endif

      ! read in the mesh
      IF (istat/=0) RETURN
      IF (shar_rank == 0) THEN
         DO ik = 1, nvoxels
            READ(iunit,*) vertex(ik,1),vertex(ik,2),vertex(ik,3), volume(ik), chi(ik)
         END DO
      END IF

      ! allocate memory for information about the mesh
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL mpialloc_1d_dbl(mx, nvoxels,shar_rank,0,shar_comm,win_mx)
         CALL mpialloc_1d_dbl(my, nvoxels,shar_rank,0,shar_comm,win_my)
         CALL mpialloc_1d_dbl(mz, nvoxels,shar_rank,0,shar_comm,win_mz)
         CALL mpialloc_2d_dbl(dx, nvoxels, nvoxels,shar_rank,0,shar_comm,win_dx)
         CALL mpialloc_2d_dbl(dy, nvoxels, nvoxels,shar_rank,0,shar_comm,win_dy)
         CALL mpialloc_2d_dbl(dz, nvoxels, nvoxels,shar_rank,0,shar_comm,win_dz)
         CALL mpialloc_2d_dbl(ndotm, nvoxels, nvoxels,shar_rank,0,shar_comm,win_ndotm)
         mydelta = CEILING(REAL(nvoxels*nvoxels) / REAL(shar_size))
         mystart = 1 + shar_rank*mydelta
         myend   = mystart + mydelta
         IF (myend > nvoxels) myend=nvoxels
      ELSE
#endif
         ALLOCATE(mx(nvoxels),my(nvoxels),mz(nvoxels),STAT=istat)
         ALLOCATE(dx(nvoxels,nvoxels),dy(nvoxels,nvoxels), &
            dz(nvoxels,nvoxels),ndotm(nvoxels,nvoxels),STAT=istat)
         mystart = 1; myend = nvoxels*nvoxels
#if defined(MPI_OPT)
      END IF
#endif
      IF (istat/=0) RETURN
      ! Precalculate information about mesh
      ! Calculate the face normal
      DO ik = mystart, myend
         i = MOD(ik-1,nvoxels)+1
         j = MOD(ik-1,nvoxels*nvoxels)+1
         j = FLOOR(REAL(j) / REAL(nvoxels))+1
         dx(i,j) = vertex(1,i)-vertex(1,j)
         dy(i,j) = vertex(2,i)-vertex(2,j)
         dz(i,j) = vertex(3,i)-vertex(3,j)
      END DO
#if defined(MPI_OPT)
      IF (PRESENT(comm)) CALL MPI_BARRIER(shar_comm,istat)
#endif
      ! close file
      CLOSE(iunit)

      RETURN
      END SUBROUTINE mumaterial_load
      
      SUBROUTINE mumaterial_init(getBfld, comm)
      !-----------------------------------------------------------------------
      ! mumaterial_load: Loads magnetic material file.
      !-----------------------------------------------------------------------
      ! fcn           : Function which returns the vacuum magnetic field
      !                 SUBROUTINE FCN(x,y,z,bx,by,bz)
      ! param[in, out]: comm. MPI communicator, handles shared memory
      !-----------------------------------------------------------------------
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      INTEGER, INTENT(inout), OPTIONAL :: comm
      EXTERNAL:: getBfld
      LOGICAL :: shared
      INTEGER :: ik, i, j, istat
      INTEGER :: win_mx0, win_my0, win_mz0, win_bx, win_by, win_bz, win_help2d,&
                 win_invd, win_invd3, win_invd4
      DOUBLE PRECISION, POINTER :: mx0(:), my0(:), mz0(:), bx(:), by(:), bz(:)
      DOUBLE PRECISION, POINTER :: help_2d(:,:), invd(:,:), invd3(:,:), invd4(:,:)

      ! Allocate helpers
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_COMM_SPLIT_TYPE(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shar_comm, istat)
         CALL MPI_COMM_RANK( shar_comm, shar_rank, istat )
         CALL MPI_COMM_SIZE( shar_comm, shar_size, istat)
         CALL mpialloc_1d_dbl(mx0, nvoxels,shar_rank,0,shar_comm,win_mx0)
         CALL mpialloc_1d_dbl(my0, nvoxels,shar_rank,0,shar_comm,win_my0)
         CALL mpialloc_1d_dbl(mz0, nvoxels,shar_rank,0,shar_comm,win_mz0)
         CALL mpialloc_1d_dbl(bx, nvoxels,shar_rank,0,shar_comm,win_bx)
         CALL mpialloc_1d_dbl(by, nvoxels,shar_rank,0,shar_comm,win_by)
         CALL mpialloc_1d_dbl(bz, nvoxels,shar_rank,0,shar_comm,win_bz)
         CALL mpialloc_2d_dbl(help_2d, nvoxels, nvoxels,shar_rank,0,shar_comm,win_help2d)
         CALL mpialloc_2d_dbl(invd, nvoxels, nvoxels,shar_rank,0,shar_comm,win_invd)
         CALL mpialloc_2d_dbl(invd3, nvoxels, nvoxels,shar_rank,0,shar_comm,win_invd3)
         CALL mpialloc_2d_dbl(invd4, nvoxels, nvoxels,shar_rank,0,shar_comm,win_invd4)
         shared = .true.
         mydelta = CEILING(REAL(nvoxels) / REAL(shar_size))
         mystart = 1 + shar_rank*mydelta
         myend   = mystart + mydelta
         IF (myend > nvoxels) myend=nvoxels
      ELSE
#endif
         ALLOCATE(mx0(nvoxels),my0(nvoxels),mz0(nvoxels),STAT=istat)
         ALLOCATE(bx(nvoxels),by(nvoxels),bz(nvoxels),STAT=istat)
         ALLOCATE(help_2d(nvoxels,nvoxels),STAT=istat)
         ALLOCATE(invd(nvoxels,nvoxels),invd3(nvoxels,nvoxels),&
                  invd4(nvoxels,nvoxels),STAT=istat)
         shared = .false.
         mystart = 1; myend = nvoxels
#if defined(MPI_OPT)
      END IF
#endif

      ! Setup helpers and IC B-applied
      DO ik = mystart,myend
         ! Calculate 1/r, 1/r^3 and 1/r^4
         invd(ik,:) = 1.0/SQRT(dx(ik,:)*dx(ik,:) + &
                               dy(ik,:)*dy(ik,:) + &
                               dz(ik,:)*dz(ik,:) )
         invd(ik,ik) = 0.0 ! Cross terms are zero
         invd3(ik,:) = invd(ik,:)**3
         invd4(ik,:) = invd3(ik,:)*invd(ik,:)
         ! Get B-applied
         CALL getBfld(vertex(ik,1),vertex(ik,2),vertex(ik,3),&
            mx0(ik), my0(ik), mz0(ik))
         ! Initial guess is that m = B-applied
         mx(ik) = mx0(ik);
         my(ik) = my0(ik);
         mz(ik) = mz0(ik);
      END DO

      ! to make code clearer below
      i = mystart; j = myend
      DO
         ! Calculate ndotm and dx*m
         ! mx term
         help_2d(i:j,:) = SPREAD(mx(i:j),2,nvoxels)
         ndotm(i:j,:) = dx(i:j,:)*help_2d(i:j,:)   
         bx(i:j) = - SUM(help_2d(i:j,:)*invd3(i:j,:),DIM=2)
         ! my term
         help_2d(i:j,:) = SPREAD(my(i:j),2,nvoxels)
         ndotm(i:j,:) = ndotm(i:j,:) + dy(i:j,:)*help_2d(i:j,:)
         by(i:j) = - SUM(help_2d(i:j,:)*invd3(i:j,:),DIM=2)
         ! mz term
         help_2d(i:j,:) = SPREAD(mz(i:j),2,nvoxels)
         ndotm(i:j,:) = ndotm(i:j,:) + dz(i:j,:)*help_2d(i:j,:)
         bz(i:j) = - SUM(help_2d(i:j,:)*invd3(i:j,:),DIM=2)
         ! Apply 1/r to get ndotm
         ndotm(i:j,:) = ndotm(i:j,:) * invd(i:j,:)

         ! 3*n_k*(ndotm) term
         bx(i:j) = bx(i:j) + 3*SUM(dx(i:j,:)*ndotm(i:j,:)*invd4(i:j,:),DIM=2)
         by(i:j) = by(i:j) + 3*SUM(dy(i:j,:)*ndotm(i:j,:)*invd4(i:j,:),DIM=2)
         bz(i:j) = bz(i:j) + 3*SUM(dz(i:j,:)*ndotm(i:j,:)*invd4(i:j,:),DIM=2)

         ! Update Magnetization  (mu0/4*pi)
         mx(i:j) = mx0(i:j) + bx(i:j)*1E-7
         my(i:j) = my0(i:j) + by(i:j)*1E-7
         mz(i:j) = mz0(i:j) + bz(i:j)*1E-7

         ! Exit if tollerance achieved
         IF (ALL(bx(i:j)/mx(i:j)<1.0E-3) .and. &
            ALL(by(i:j)/my(i:j)<1.0E-3) .and. &
            ALL(bz(i:j)/mz(i:j)<1.0E-3)) EXIT

      END DO

#if defined(MPI_OPT)
      IF (PRESENT(comm)) CALL MPI_BARRIER(shar_comm,istat)
#endif

      ! Now make normalization 
      DO ik = mystart, myend
            mx(ik) = mx(ik)*volume(ik)*chi(ik)
            my(ik) = my(ik)*volume(ik)*chi(ik)
            mz(ik) = mz(ik)*volume(ik)*chi(ik)
      END DO

      ! Deallocate helpers
      CALL free_mpi_array1d_dbl(win_mx0,mx0,shared)
      CALL free_mpi_array1d_dbl(win_my0,my0,shared)
      CALL free_mpi_array1d_dbl(win_mz0,mz0,shared)
      CALL free_mpi_array1d_dbl(win_bx,bx,shared)
      CALL free_mpi_array1d_dbl(win_by,by,shared)
      CALL free_mpi_array1d_dbl(win_bz,bz,shared)
      CALL free_mpi_array2d_dbl(win_help2d,help_2d,shared)
      CALL free_mpi_array2d_dbl(win_invd,invd,shared)
      CALL free_mpi_array2d_dbl(win_invd3,invd3,shared)
      CALL free_mpi_array2d_dbl(win_invd4,invd4,shared)


      RETURN
      END SUBROUTINE mumaterial_init
      
      SUBROUTINE mumaterial_getBfield(x,y,z,bx,by,bz)
      !-----------------------------------------------------------------------
      ! mumaterial_load: Loads magnetic material file.
      !-----------------------------------------------------------------------
      ! param[in]: x,y,z Evaluation point [m]
      ! param[out]: bx, by bz Magnetic Field [T]
      !-----------------------------------------------------------------------
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in)  :: x, y, z
      DOUBLE PRECISION, INTENT(out) :: bx, by, bz
      INTEGER :: ik
      DOUBLE PRECISION :: dx, dy, dz, d, invd5, invd3, ddotm

      bx = 0; by =0; bz = 0
      DO ik = 1, nvoxels
            dx = x - vertex(1,ik)
            dy = y - vertex(2,ik)
            dz = z - vertex(3,ik)
            d  = SQRT(dx*dx+dy*dy+dz*dz)
            invd5 = 1.0/(d*d*d*d*d)
            invd3 = d*d*invd5
            ddotm = dx*mx(ik)+dy*my(ik)+dz*mz(ik)*invd5
            bx = bx + 3*dx*ddotm - mx(ik)*invd3
            by = by + 3*dy*ddotm - my(ik)*invd3
            bz = bz + 3*dz*ddotm - mz(ik)*invd3
      END DO
      bx = bx*1E-7
      by = by*1E-7
      bz = bz*1E-7

      RETURN
      END SUBROUTINE mumaterial_getBfield



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
      END MODULE mumaterial_mod