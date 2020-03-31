!-----------------------------------------------------------------------
!     Module:        wall_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/02/2012
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
      
!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL            :: lwall_loaded
      INTEGER            :: nvertex, nface
      INTEGER, POINTER :: face(:,:)
      INTEGER, POINTER :: ihit_array(:)
      DOUBLE PRECISION, POINTER   :: vertex(:,:)
      CHARACTER(LEN=256) :: machine_string
      CHARACTER(LEN=256) :: date


      LOGICAL, PRIVATE, ALLOCATABLE            :: lmask(:)
      INTEGER, PRIVATE                         :: mystart, myend, mydelta, ik_min
      INTEGER, PRIVATE                         :: win_vertex, win_face, &
                                                  win_fn, win_a0, win_v0, win_v1, &
                                                  win_dot00, win_dot01, win_dot11, &
                                                  win_d, win_ihit, win_invDenom, &
                                                  win_xmin, win_ymin, win_zmin, &
                                                  win_xmax, win_ymax, win_zmax
      DOUBLE PRECISION, PRIVATE, POINTER   :: FN(:,:), d(:), t(:), r0(:,:), dr(:,:)
      DOUBLE PRECISION, PRIVATE, POINTER   :: A0(:,:), V0(:,:), V1(:,:), V2(:,:),&
                                                  DOT00(:), DOT01(:), DOT02(:),&
                                                  DOT11(:), DOT12(:),&
                                                  invDenom(:), alpha(:), beta(:), &
                                                  xmin(:),ymin(:),zmin(:), &
                                                  xmax(:),ymax(:),zmax(:)
      DOUBLE PRECISION, PRIVATE, PARAMETER      :: zero = 0.0D+0
      DOUBLE PRECISION, PRIVATE, PARAMETER      :: one  = 1.0D+0
      
!-----------------------------------------------------------------------
!     Subroutines
!         wall_load_txt:   Loads trianular mesh from file
!         wall_load_mn:    Creates wall from harmonics
!         wall_dump:       Dumps triangulation data
!         wall_info:       Prints wall info.
!         wall_collide:    Calculates collision with wall
!         wall_free:       Frees module memory
!-----------------------------------------------------------------------
      INTERFACE collide
         MODULE PROCEDURE collide_double, collide_float
      END INTERFACE

      PRIVATE :: mpialloc_1d_int,mpialloc_1d_dbl,mpialloc_2d_int,mpialloc_2d_dbl
      CONTAINS
      
      SUBROUTINE wall_load_txt(filename,istat,comm)
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER, INTENT(inout)       :: istat
      INTEGER, INTENT(inout), OPTIONAL :: comm
      INTEGER :: iunit, ik, dex1, dex2, dex3, bubble
      INTEGER :: shar_comm, shar_rank, shar_size
      DOUBLE PRECISION :: rt1,rt2,rt3
      DOUBLE PRECISION, DIMENSION(3) :: temp
      DOUBLE PRECISION, ALLOCATABLE :: r_temp(:,:),z_temp(:,:)
      
      shar_rank = 0; shar_size = 1;
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_COMM_SPLIT_TYPE(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shar_comm, istat)
         CALL MPI_COMM_RANK( shar_comm, shar_rank, istat )
         CALL MPI_COMM_SIZE( shar_comm, shar_size, istat)
      END IF
#endif
      CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
      IF (istat/=0) RETURN
      READ(iunit,'(A)') machine_string
      READ(iunit,'(A)') date
      READ(iunit,*) nvertex,nface
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
      IF (istat/=0) RETURN
      IF (shar_rank == 0) THEN
         DO ik = 1, nvertex
            READ(iunit,*) vertex(ik,1),vertex(ik,2),vertex(ik,3)
         END DO
         DO ik=1,nface
            READ(iunit,*) face(ik,1),face(ik,2),face(ik,3)
         END DO
      END IF
      CLOSE(iunit)
#if defined(MPI_OPT)
      IF (PRESENT(comm)) CALL MPI_BARRIER(comm,istat)
      IF (istat/=0) RETURN
      IF (PRESENT(comm)) THEN
         CALL mpialloc_2d_dbl(A0,nface,3,shar_rank,0,shar_comm,win_a0)
         CALL mpialloc_2d_dbl(V0,nface,3,shar_rank,0,shar_comm,win_v0)
         CALL mpialloc_2d_dbl(V1,nface,3,shar_rank,0,shar_comm,win_v1)
         CALL mpialloc_2d_dbl(FN,nface,3,shar_rank,0,shar_comm,win_fn)
         CALL mpialloc_1d_dbl(xmin,nface,shar_rank,0,shar_comm,win_xmin)
         CALL mpialloc_1d_dbl(ymin,nface,shar_rank,0,shar_comm,win_ymin)
         CALL mpialloc_1d_dbl(zmin,nface,shar_rank,0,shar_comm,win_zmin)
         CALL mpialloc_1d_dbl(xmax,nface,shar_rank,0,shar_comm,win_xmax)
         CALL mpialloc_1d_dbl(ymax,nface,shar_rank,0,shar_comm,win_ymax)
         CALL mpialloc_1d_dbl(zmax,nface,shar_rank,0,shar_comm,win_zmax)
         mydelta = CEILING(REAL(nface) / REAL(shar_size))
         mystart = 1 + shar_rank*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
#endif
         ALLOCATE(A0(nface,3),V0(nface,3),V1(nface,3),&
                  FN(nface,3),STAT=istat)
         ALLOCATE(xmin(nface),ymin(nface),zmin(nface),&
            xmax(nface),ymax(nface),zmax(nface),STAT=istat)
         mystart = 1; myend = nface
#if defined(MPI_OPT)
      END IF
#endif
      IF (istat/=0) RETURN
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
         xmin(ik) = min(vertex(dex1,1),vertex(dex2,1),vertex(dex3,1))
         xmax(ik) = max(vertex(dex1,1),vertex(dex2,1),vertex(dex3,1))
         ymin(ik) = min(vertex(dex1,2),vertex(dex2,2),vertex(dex3,2))
         ymax(ik) = max(vertex(dex1,2),vertex(dex2,2),vertex(dex3,2))
         zmin(ik) = min(vertex(dex1,3),vertex(dex2,3),vertex(dex3,3))
         zmax(ik) = max(vertex(dex1,3),vertex(dex2,3),vertex(dex3,3))
      END DO
#if defined(MPI_OPT)
      IF (PRESENT(comm)) CALL MPI_BARRIER(comm,istat)
#endif
      ! Check for zero area
      IF (ANY(SUM(FN*FN,DIM=2)==zero)) THEN
         istat=-327
         RETURN
      END IF
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
      IF (istat/=0) RETURN
      DO ik = mystart, myend
         ihit_array(ik) = 0
         DOT00(ik) = V0(ik,1)*V0(ik,1) + V0(ik,2)*V0(ik,2) + V0(ik,3)*V0(ik,3)
         DOT01(ik) = V0(ik,1)*V1(ik,1) + V0(ik,2)*V1(ik,2) + V0(ik,3)*V1(ik,3)
         DOT11(ik) = V1(ik,1)*V1(ik,1) + V1(ik,2)*V1(ik,2) + V1(ik,3)*V1(ik,3)
         d(ik)     = FN(ik,1)*A0(ik,1) + FN(ik,2)*A0(ik,2) + FN(ik,3)*A0(ik,3)
         invDenom(ik) = one / (DOT00(ik)*DOT11(ik) - DOT01(ik)*DOT01(ik))
      END DO
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_BARRIER(shar_comm, istat)
         CALL MPI_COMM_FREE(shar_comm, istat)
      END IF
#endif
      lwall_loaded = .true.
      RETURN
      END SUBROUTINE wall_load_txt

      SUBROUTINE wall_load_mn(Rmn,Zmn,xm,xn,mn,nu,nv,comm)
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: Rmn(mn), Zmn(mn), xm(mn), xn(mn)
      INTEGER, INTENT(in) :: mn, nu, nv
      INTEGER, INTENT(inout), OPTIONAL :: comm
      INTEGER :: u, v, i, j, istat, dex1, dex2, dex3, ik, bubble, nv2
      INTEGER :: shar_comm, shar_rank, shar_size
      DOUBLE PRECISION :: pi2, th, zt, pi
      DOUBLE PRECISION, DIMENSION(3) :: temp
      DOUBLE PRECISION, ALLOCATABLE :: r_temp(:,:),z_temp(:,:),x_temp(:,:),y_temp(:,:)

      machine_string = '          HARMONICS'
      date = '      TODAY'
      pi2 = 8.0D+00 * ATAN(one)
      pi  = 4.0E+00 * ATAN(one)
      nv2 = nv/2
      shar_rank = 0; shar_size = 1;
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_COMM_SPLIT_TYPE(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shar_comm, istat)
         CALL MPI_COMM_RANK( shar_comm, shar_rank, istat )
         CALL MPI_COMM_SIZE( shar_comm, shar_size, istat)
      END IF
#endif
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
      i = 1  ! Tracks vertex index
      j = 1 ! Tracks face index
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
         DEALLOCATE(r_temp,z_temp,x_temp,y_temp)
      END IF
#if defined(MPI_OPT)
      IF (PRESENT(comm)) CALL MPI_BARRIER(shar_comm,istat)
#endif
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL mpialloc_2d_dbl(A0,nface,3,shar_rank,0,shar_comm,win_a0)
         CALL mpialloc_2d_dbl(V0,nface,3,shar_rank,0,shar_comm,win_v0)
         CALL mpialloc_2d_dbl(V1,nface,3,shar_rank,0,shar_comm,win_v1)
         CALL mpialloc_2d_dbl(FN,nface,3,shar_rank,0,shar_comm,win_fn)
         CALL mpialloc_1d_dbl(xmin,nface,shar_rank,0,shar_comm,win_xmin)
         CALL mpialloc_1d_dbl(ymin,nface,shar_rank,0,shar_comm,win_ymin)
         CALL mpialloc_1d_dbl(zmin,nface,shar_rank,0,shar_comm,win_zmin)
         CALL mpialloc_1d_dbl(xmax,nface,shar_rank,0,shar_comm,win_xmax)
         CALL mpialloc_1d_dbl(ymax,nface,shar_rank,0,shar_comm,win_ymax)
         CALL mpialloc_1d_dbl(zmax,nface,shar_rank,0,shar_comm,win_zmax)
         mydelta = CEILING(REAL(nface) / REAL(shar_size))
         mystart = 1 + shar_rank*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
#endif
         ALLOCATE(A0(nface,3),V0(nface,3),V1(nface,3),&
                  FN(nface,3),STAT=istat)
         ALLOCATE(xmin(nface),ymin(nface),zmin(nface),&
            xmax(nface),ymax(nface),zmax(nface),STAT=istat)
         mystart = 1; myend = nface
#if defined(MPI_OPT)
      END IF
#endif
      IF (istat/=0) RETURN
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
         xmin(ik) = min(vertex(dex1,1),vertex(dex2,1),vertex(dex3,1))
         xmax(ik) = max(vertex(dex1,1),vertex(dex2,1),vertex(dex3,1))
         ymin(ik) = min(vertex(dex1,2),vertex(dex2,2),vertex(dex3,2))
         ymax(ik) = max(vertex(dex1,2),vertex(dex2,2),vertex(dex3,2))
         zmin(ik) = min(vertex(dex1,3),vertex(dex2,3),vertex(dex3,3))
         zmax(ik) = max(vertex(dex1,3),vertex(dex2,3),vertex(dex3,3))
      END DO
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
      IF (istat/=0) RETURN
      DO ik = mystart, myend
         ihit_array(ik) = 0
         DOT00(ik) = V0(ik,1)*V0(ik,1) + V0(ik,2)*V0(ik,2) + V0(ik,3)*V0(ik,3)
         DOT01(ik) = V0(ik,1)*V1(ik,1) + V0(ik,2)*V1(ik,2) + V0(ik,3)*V1(ik,3)
         DOT11(ik) = V1(ik,1)*V1(ik,1) + V1(ik,2)*V1(ik,2) + V1(ik,3)*V1(ik,3)
         d(ik)     = FN(ik,1)*A0(ik,1) + FN(ik,2)*A0(ik,2) + FN(ik,3)*A0(ik,3)
         invDenom(ik) = one / (DOT00(ik)*DOT11(ik) - DOT01(ik)*DOT01(ik))
      END DO
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_BARRIER(shar_comm, istat)
         CALL MPI_COMM_FREE(shar_comm, istat)
      END IF
#endif
      lwall_loaded = .true.
      RETURN
      END SUBROUTINE wall_load_mn

      SUBROUTINE wall_dump(filename,istat)
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER, INTENT(inout)       :: istat
      INTEGER :: iunit, ik
      CALL safe_open(iunit,istat,'wall_dump.'//TRIM(filename),'unknown','formatted')
      DO ik = 1, nface
         WRITE(iunit,'(13(ES20.10))')  A0(ik,1), A0(ik,2), A0(ik,3), &
                                       FN(ik,1), FN(ik,2), FN(ik,3),&
                                       V0(ik,1), V0(ik,2), V0(ik,3),&
                                       V1(ik,1), V1(ik,2), V1(ik,3), &
                                       d(ik)
      END DO
      WRITE(iunit,*) '#  V0(ik,1), V0(ik,2), V0(ik,3), FN(ik,1), FN(ik,2), FN(ik,3),',&
                        'V(ik,1), V(ik,2), V(ik,3), W(ik,1), W(ik,2), W(ik,3), d(ik)'
      CLOSE(iunit)
      RETURN
      END SUBROUTINE wall_dump

      
      SUBROUTINE wall_info(iunit)
      IMPLICIT NONE
      INTEGER, INTENT(in)    :: iunit
      WRITE(iunit,'(A)')         ' -----  Vessel Information  -----'
      WRITE(iunit,'(3X,A,A)')    'Wall Name : ',TRIM(machine_string(10:))
      WRITE(iunit,'(3X,A,A)')    'Date      : ',TRIM(date(6:))
      WRITE(iunit,'(3X,A,I7)')   'Faces     : ',nface
      RETURN
      END SUBROUTINE wall_info

      SUBROUTINE collide_float(x0,y0,z0,x1,y1,z1,xw,yw,zw,lhit)
      IMPLICIT NONE
      REAL, INTENT(in) :: x0, y0, z0, x1, y1, z1
      REAL, INTENT(out) :: xw, yw, zw
      LOGICAL, INTENT(out) :: lhit
      DOUBLE PRECISION :: x0d, y0d, z0d, x1d, y1d, z1d
      DOUBLE PRECISION :: xwd, ywd, zwd
      LOGICAL          :: lhit2
      xw=zero; yw=zero; zw=zero; lhit=.FALSE.
      x0d=x0; y0d=y0; z0d=z0
      x1d=x1; y1d=y1; z1d=z1
      CALL collide_double(x0d,y0d,z0d,x1d,y1d,z1d,xwd,ywd,zwd,lhit2)
      xw=xwd; yw=ywd; zw=zwd; lhit=lhit2
      RETURN
      END SUBROUTINE collide_float

      SUBROUTINE collide_float_vec(x0,y0,z0,x1,y1,z1,xw,yw,zw,lhit)
      IMPLICIT NONE
      REAL, INTENT(in) :: x0, y0, z0, x1, y1, z1
      REAL, INTENT(out) :: xw, yw, zw
      LOGICAL, INTENT(out) :: lhit
      DOUBLE PRECISION :: x0d, y0d, z0d, x1d, y1d, z1d
      DOUBLE PRECISION :: xwd, ywd, zwd
      LOGICAL          :: lhit2
      xw=zero; yw=zero; zw=zero; lhit=.FALSE.
      x0d=x0; y0d=y0; z0d=z0
      x1d=x1; y1d=y1; z1d=z1
      CALL collide_double_vec(x0d,y0d,z0d,x1d,y1d,z1d,xwd,ywd,zwd,lhit2)
      xw=xwd; yw=ywd; zw=zwd; lhit=lhit2
      RETURN
      END SUBROUTINE collide_float_vec

      SUBROUTINE collide_double_vec(x0,y0,z0,x1,y1,z1,xw,yw,zw,lhit)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x0, y0, z0, x1, y1, z1
      DOUBLE PRECISION, INTENT(out) :: xw, yw, zw
      LOGICAL, INTENT(out) :: lhit
      INTEGER :: ik, k1,k2
      DOUBLE PRECISION :: drx, dry, drz, V2x, V2y, V2z, DOT02l, DOT12l, tloc, tmin, alphal, betal
      REAL :: xs,ys,zs,xb,yb,zb
      xw=zero; yw=zero; zw=zero; lhit=.FALSE.
      ik_min = zero
      tmin = 2
      k1 = 1; k2 = nface
      ! Define DR
      drx = x1-x0
      dry = y1-y0
      drz = z1-z0
      ! min/max
      xs = min(x0,x1)
      ys = min(y0,y1)
      zs = min(z0,z1)
      xb = max(x0,x1)
      yb = max(y0,y1)
      zb = max(z0,z1)
      ! Calculate distance along trajectory to hit triangle
      DO ik = k1,k2
         IF (xb<xmin(ik) .or. xs>xmax(ik) .or. &
             yb<ymin(ik) .or. ys>ymax(ik) .or. &
             zb<zmin(ik) .or. zs>zmax(ik)) CYCLE
         alphal = FN(ik,1)*drx + FN(ik,2)*dry + FN(ik,3)*drz
         betal = FN(ik,1)*x0 + FN(ik,2)*y0 + FN(ik,3)*z0
         !IF (alphal < zero) CYCLE  ! we get wrong face
         tloc = (d(ik)-betal)/alphal
         IF (tloc > one) CYCLE
         IF (tloc <= zero) CYCLE
         V2x = x0 + tloc*drx - A0(ik,1)
         V2y = y0 + tloc*dry - A0(ik,2)
         V2z = z0 + tloc*drz - A0(ik,3)
         DOT02l = V0(ik,1)*V2x + V0(ik,2)*V2y + V0(ik,3)*V2z
         DOT12l = V1(ik,1)*V2x + V1(ik,2)*V2y + V1(ik,3)*V2z
         alphal = (DOT11(ik)*DOT02l-DOT01(ik)*DOT12l)*invDenom(ik)
         betal  = (DOT00(ik)*DOT12l-DOT01(ik)*DOT02l)*invDenom(ik)
         IF ((alphal < zero) .or. (betal < zero) .or. (alphal+betal > one)) CYCLE
         IF (tloc < tmin) THEN
            ik_min = ik
            tmin = tloc
         END IF
      END DO
      IF (ik_min > zero) THEN
         lhit = .TRUE.
         xw   = x0 + tmin*drx
         yw   = y0 + tmin*dry
         zw   = z0 + tmin*drz
         ihit_array(ik_min) = ihit_array(ik_min) + 1
      END IF
      RETURN
      END SUBROUTINE collide_double_vec

      SUBROUTINE collide_double(x0,y0,z0,x1,y1,z1,xw,yw,zw,lhit)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x0, y0, z0, x1, y1, z1
      DOUBLE PRECISION, INTENT(out) :: xw, yw, zw
      LOGICAL, INTENT(out) :: lhit
      INTEGER :: ik, k1,k2
      DOUBLE PRECISION :: drx, dry, drz, V2x, V2y, V2z, DOT02l, DOT12l, tloc, tmin, alphal, betal
      xw=zero; yw=zero; zw=zero; lhit=.FALSE.
      ik_min = zero
      tmin = 2
      k1 = 1; k2 = nface
      ! Define DR
      drx = x1-x0
      dry = y1-y0
      drz = z1-z0
      ! Calculate distance along trajectory to hit triangle
      DO ik = k1,k2
         alphal = FN(ik,1)*drx + FN(ik,2)*dry + FN(ik,3)*drz
         betal = FN(ik,1)*x0 + FN(ik,2)*y0 + FN(ik,3)*z0
         !IF (alphal < zero) CYCLE  ! we get wrong face
         tloc = (d(ik)-betal)/alphal
         IF (tloc > one) CYCLE
         IF (tloc <= zero) CYCLE
         V2x = x0 + tloc*drx - A0(ik,1)
         V2y = y0 + tloc*dry - A0(ik,2)
         V2z = z0 + tloc*drz - A0(ik,3)
         DOT02l = V0(ik,1)*V2x + V0(ik,2)*V2y + V0(ik,3)*V2z
         DOT12l = V1(ik,1)*V2x + V1(ik,2)*V2y + V1(ik,3)*V2z
         alphal = (DOT11(ik)*DOT02l-DOT01(ik)*DOT12l)*invDenom(ik)
         betal  = (DOT00(ik)*DOT12l-DOT01(ik)*DOT02l)*invDenom(ik)
         IF ((alphal < zero) .or. (betal < zero) .or. (alphal+betal > one)) CYCLE
         IF (tloc < tmin) THEN
            ik_min = ik
            tmin = tloc
         END IF
      END DO
      IF (ik_min > zero) THEN
         lhit = .TRUE.
         xw   = x0 + tmin*drx
         yw   = y0 + tmin*dry
         zw   = z0 + tmin*drz
         ihit_array(ik_min) = ihit_array(ik_min) + 1
      END IF
      RETURN
      END SUBROUTINE collide_double

      SUBROUTINE uncount_wall_hit
         IMPLICIT NONE
         ihit_array(ik_min) = ihit_array(ik_min) - 1
      END SUBROUTINE

      INTEGER FUNCTION get_wall_ik()
         IMPLICIT NONE
         get_wall_ik = ik_min
         RETURN
      END FUNCTION

      DOUBLE PRECISION FUNCTION get_wall_area(ik)
         IMPLICIT NONE
         INTEGER, INTENT(in) :: ik
         get_wall_area = 0.5*SQRT(SUM(FN(ik,:)*FN(ik,:)))
         RETURN
      END FUNCTION

      SUBROUTINE wall_free(istat,shared_comm)
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: istat
      INTEGER, INTENT(inout), OPTIONAL :: shared_comm
      INTEGER :: ier
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
         CALL MPI_WIN_FENCE(0,win_invdenom,istat)
         CALL MPI_WIN_FREE(win_invdenom,istat)
         CALL MPI_WIN_FENCE(0,win_ihit,istat)
         CALL MPI_WIN_FREE(win_ihit,istat)
         CALL MPI_WIN_FENCE(0,win_xmin,istat)
         CALL MPI_WIN_FREE(win_xmin,istat)
         CALL MPI_WIN_FENCE(0,win_ymin,istat)
         CALL MPI_WIN_FREE(win_ymin,istat)
         CALL MPI_WIN_FENCE(0,win_zmin,istat)
         CALL MPI_WIN_FREE(win_zmin,istat)
         CALL MPI_WIN_FENCE(0,win_xmax,istat)
         CALL MPI_WIN_FREE(win_xmax,istat)
         CALL MPI_WIN_FENCE(0,win_ymax,istat)
         CALL MPI_WIN_FREE(win_ymax,istat)
         CALL MPI_WIN_FENCE(0,win_zmax,istat)
         CALL MPI_WIN_FREE(win_zmax,istat)
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
         IF (ASSOCIATED(xmin)) NULLIFY(xmin)
         IF (ASSOCIATED(ymin)) NULLIFY(ymin)
         IF (ASSOCIATED(zmin)) NULLIFY(zmin)
         IF (ASSOCIATED(xmax)) NULLIFY(xmax)
         IF (ASSOCIATED(ymax)) NULLIFY(ymax)
         IF (ASSOCIATED(zmax)) NULLIFY(zmax)
      ELSE
#endif
         IF (ASSOCIATED(FN)) DEALLOCATE(FN)
         IF (ASSOCIATED(A0)) DEALLOCATE(A0)
         IF (ASSOCIATED(V0)) DEALLOCATE(V0)
         IF (ASSOCIATED(V1)) DEALLOCATE(V1)
         IF (ASSOCIATED(DOT00)) DEALLOCATE(DOT00)
         IF (ASSOCIATED(DOT01)) DEALLOCATE(DOT01)
         IF (ASSOCIATED(DOT11)) DEALLOCATE(DOT11)
         IF (ASSOCIATED(invDenom)) DEALLOCATE(invDenom)
         IF (ASSOCIATED(d)) DEALLOCATE(d)
         IF (ASSOCIATED(vertex)) DEALLOCATE(vertex)
         IF (ASSOCIATED(face)) DEALLOCATE(face)
         IF (ASSOCIATED(ihit_array)) DEALLOCATE(ihit_array)
         IF (ASSOCIATED(xmin)) DEALLOCATE(xmin)
         IF (ASSOCIATED(ymin)) DEALLOCATE(ymin)
         IF (ASSOCIATED(zmin)) DEALLOCATE(zmin)
         IF (ASSOCIATED(xmax)) DEALLOCATE(xmax)
         IF (ASSOCIATED(ymax)) DEALLOCATE(ymax)
         IF (ASSOCIATED(zmax)) DEALLOCATE(zmax)
      END IF
      machine_string=''
      date=''
      nface = -1
      nvertex = -1
      lwall_loaded = .false.
      RETURN
      END SUBROUTINE wall_free

      SUBROUTINE wall_test
      LOGICAL :: lhit
      INTEGER :: i
      DOUBLE PRECISION :: pi2,x0,y0,z0,r0,rho0,x1,y1,z1,xw,yw,zw, rr0,zz0,phi0, rho1, rr1, zz1
      r0 = 10.0   ! R0
      rho0 = 0.975
      rho1 = 1.025 ! rho
      z0 = 0 ! Z0
      pi2 = (8.0 * ATAN(1.0))
      phi0 = 0
      CALL collide(r0+rho0,DBLE(0.0),z0,r0+rho1,DBLE(0.0),z0,xw,yw,zw,lhit)
      WRITE(6,*) xw,yw,zw,lhit
      DO i = 1, 10000
         rr0 = rho0*cos(pi2*(i-1)/10000.)
         zz0 = rho0*sin(pi2*(i-1)/10000.)
         rr1 = rho1*cos(pi2*(i-1)/10000.)
         zz1 = rho1*sin(pi2*(i-1)/10000.)
         x0 = (r0+rr0)*cos(phi0)
         y0 = (r0+rr0)*sin(phi0)
         z0 = z0 + zz1
         x1 = (r0+rr1)*cos(phi0)
         y1 = (r0+rr1)*sin(phi0)
         z1 = z0 + zz1
         CALL collide(x0,y0,z0,x1,y1,z1,xw,yw,zw,lhit)
         !WRITE(327,'(9(1X,E22.12))') x0,y0,z0,x1,y1,z1,xw,yw,zw
         CALL FLUSH(327)
      END DO
      RETURN
      END SUBROUTINE wall_test

      SUBROUTINE mpialloc_1d_int(array,n1,subid,mymaster,share_comm,win)
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

!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE wall_mod
