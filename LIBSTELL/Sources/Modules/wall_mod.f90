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
!DEC$ IF DEFINED (MPI_OPT)
      USE mpi
!DEC$ ENDIF
      
!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER            :: nvertex, nface
      INTEGER, ALLOCATABLE :: face(:,:)
      INTEGER, ALLOCATABLE :: ihit_array(:)
      DOUBLE PRECISION, ALLOCATABLE   :: vertex(:,:)
      CHARACTER(LEN=256) :: machine_string
      CHARACTER(LEN=256) :: date


      LOGICAL, PRIVATE, ALLOCATABLE            :: lmask(:)
      INTEGER, PRIVATE                         :: myid_wall, nproc_wall, mystart, myend, mydelta
      INTEGER, PRIVATE                         :: win_vertex, win_face, win_phi, &
                                                  win_fn, win_a0, win_v0, win_v1, &
                                                  win_dot00, win_dot01, win_dot11, &
                                                  win_d, win_ihit, win_invDenom
      DOUBLE PRECISION, PRIVATE, ALLOCATABLE   :: FN(:,:), d(:), t(:), r0(:,:), dr(:,:)
      DOUBLE PRECISION, PRIVATE, ALLOCATABLE   :: A0(:,:), V0(:,:), V1(:,:), V2(:,:),&
                                                  DOT00(:), DOT01(:), DOT02(:),&
                                                  DOT11(:), DOT12(:),&
                                                  invDenom(:), alpha(:), beta(:), PHI(:)
      DOUBLE PRECISION, PRIVATE, PARAMETER      :: zero = 0.0D+0
      DOUBLE PRECISION, PRIVATE, PARAMETER      :: one  = 1.0D+0
      
!-----------------------------------------------------------------------
!     Subroutines
!         wall_load_txt:   Loads trianular mesh from file
!         wall_load_mn:    Creates wall from harmonics
!         wall_dump:       Dumps triangulation data
!         wall_info:       Prints wall info.
!         wall_collide:    Calculates collision with wall
!         wall_free:       Free's module memory
!-----------------------------------------------------------------------
      INTERFACE collide
         MODULE PROCEDURE collide_double, collide_float
      END INTERFACE
      CONTAINS
      
      SUBROUTINE wall_load_txt(filename,istat,shared_comm)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER, INTENT(inout)       :: istat
      INTEGER, INTENT(inout), OPTIONAL :: shared_comm
      INTEGER :: iunit, ik, dex1, dex2, dex3, bubble
      DOUBLE PRECISION, DIMENSION(3) :: temp
      DOUBLE PRECISION, ALLOCATABLE :: r_temp(:,:),z_temp(:,:)
      
      myid_wall = 0; nproc_wall = 1;
      IF (PRESENT(shared_comm)) THEN
         CALL MPI_COMM_RANK( shared_comm, myid_wall, istat )
         CALL MPI_COMM_SIZE( shared_comm, nproc_wall, istat)
      END IF
      CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
      IF (istat/=0) RETURN
      READ(iunit,'(A)') machine_string
      READ(iunit,'(A)') date
      READ(iunit,*) nvertex,nface
      IF (PRESENT(shared_comm)) THEN
         CALL mpialloc_1d_dbl(PHI,nface,myid_wall,0,shared_comm,win_phi)
         CALL mpialloc_2d_dbl(vertex,nvertex,3,myid_wall,0,shared_comm,win_vertex)
         CALL mpialloc_2d_int(face,nvertex,3,myid_wall,0,shared_comm,win_face)
         mydelta = CEILING(REAL(nface) / REAL(nproc_wall))
         mystart = 1 + (myid_wall-1)*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
         ALLOCATE(PHI(nface),STAT=istat)
         ALLOCATE(vertex(nvertex,3),face(nface,3),STAT=istat)
         mystart = 1; myend=nface
      END IF
      IF (istat/=0) RETURN
      IF (myid_wall == 0) THEN
         DO ik=1,nvertex
            READ(iunit,*) vertex(ik,1),vertex(ik,2),vertex(ik,3)
         END DO
         DO ik=1,nface
            READ(iunit,*) face(ik,1),face(ik,2),face(ik,3)
         END DO
      END IF
      CLOSE(iunit)
      IF (PRESENT(shared_comm)) CALL MPI_WIN_FENCE(0, win_vertex, istat)
      IF (PRESENT(shared_comm)) CALL MPI_WIN_FENCE(0, win_face, istat)
      IF (istat/=0) RETURN
      ! Sort the array by toroidal angle
      DO ik = mystart, myend
         dex1 = face(ik,1)
         PHI(ik) = ATAN2(vertex(dex1,2),vertex(dex1,1))
      END DO
      IF (PRESENT(shared_comm)) CALL MPI_WIN_FENCE(0, win_phi, istat)
      ! Bubble Sort
!      dex3 = nface
!      DO WHILE (dex3 > 1)
!         bubble = 0 !bubble in the greatest element out of order
!         DO ik = 1, (dex3-1)
!            IF (PHI(ik) > PHI(ik+1)) THEN
!               ! Adjust PHI
!               temp(1) = PHI(ik)
!               PHI(ik) = PHI(ik+1)
!               PHI(ik+1) = temp(1)
!               ! Adjust VERTEX
!               dex1 = face(ik,1)
!               dex2 = face(ik+1,1)
!               temp = vertex(dex1,:)
!               vertex(dex1,:) = vertex(dex2,:)
!               vertex(dex2,:) = temp
!               ! Adjust the faces
!               WHERE(face == dex1) face = 0
!               WHERE(face == dex2) face = dex1
!               WHERE(face == 0) face = dex2
!               ! Adjust bubble
!               bubble = ik
!            END IF 
!         END DO
!         dex3 = bubble
!      END DO
!      ALLOCATE(A0(nface,3),V0(nface,3),V1(nface,3),&
!               V2(nface,3),FN(nface,3),STAT=istat)
      IF (PRESENT(shared_comm)) THEN
         CALL mpialloc_2d_dbl(A0,nface,3,myid_wall,0,shared_comm,win_a0)
         CALL mpialloc_2d_dbl(V0,nface,3,myid_wall,0,shared_comm,win_v0)
         CALL mpialloc_2d_dbl(V1,nface,3,myid_wall,0,shared_comm,win_v1)
         CALL mpialloc_2d_dbl(FN,nface,3,myid_wall,0,shared_comm,win_fn)
         mydelta = CEILING(REAL(nface) / REAL(nproc_wall))
         mystart = 1 + (myid_wall-1)*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
         ALLOCATE(A0(nface,3),V0(nface,3),V1(nface,3),&
                  FN(nface,3),STAT=istat)
         mystart = 1; myend = nface
      END IF
      IF (istat/=0) RETURN
      ! Calculate the face normal
      ! V  = Vertex1-Vertex0
      ! W  = Vertex2-Vertex0
      ! FN = VxW/|VxW|
      ! d  = -Vertex0.FN (. is dot product) (note we've absorbed the negative)
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
      IF (PRESENT(shared_comm)) CALL MPI_WIN_FENCE(0, win_a0, istat)
      IF (PRESENT(shared_comm)) CALL MPI_WIN_FENCE(0, win_v0, istat)
      IF (PRESENT(shared_comm)) CALL MPI_WIN_FENCE(0, win_v1, istat)
      IF (PRESENT(shared_comm)) CALL MPI_WIN_FENCE(0, win_fn, istat)
      ! Check for zero area
      IF (ANY(SUM(FN*FN,DIM=2)==zero)) THEN
         istat=-327
         RETURN
      END IF
      !DEALLOCATE(vertex,face)
!      ALLOCATE(DOT00(nface), DOT01(nface), DOT02(nface),&
!               DOT11(nface), DOT12(nface),invDenom(nface),&
!               STAT=istat)
      IF (PRESENT(shared_comm)) THEN
         CALL mpialloc_1d_dbl(DOT00,nface,myid_wall,0,shared_comm,win_dot00)
         CALL mpialloc_1d_dbl(DOT01,nface,myid_wall,0,shared_comm,win_dot01)
         CALL mpialloc_1d_dbl(DOT11,nface,myid_wall,0,shared_comm,win_dot11)
         CALL mpialloc_1d_dbl(invDenom,nface,myid_wall,0,shared_comm,win_invDenom)
         CALL mpialloc_1d_dbl(d,nface,myid_wall,0,shared_comm,win_d)
         CALL mpialloc_1d_int(ihit_array,nface,myid_wall,0,shared_comm,win_ihit)
         mydelta = CEILING(REAL(nface) / REAL(nproc_wall))
         mystart = 1 + (myid_wall-1)*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
         ALLOCATE(DOT00(nface), DOT01(nface),&
                  DOT11(nface), invDenom(nface),&
                  STAT=istat)
         ALLOCATE(d(nface),STAT=istat)
         ALLOCATE(ihit_array(nface),STAT=istat)
         mystart = 1; myend = nface
      END IF
      IF (istat/=0) RETURN
      DO ik = mystart, myend
         ihit_array(ik) = 0
         DOT00(ik) = V0(ik,1)*V0(ik,1) + V0(ik,2)*V0(ik,2) + V0(ik,3)*V0(ik,3)
         DOT01(ik) = V0(ik,1)*V1(ik,1) + V0(ik,2)*V1(ik,2) + V0(ik,3)*V1(ik,3)
         DOT11(ik) = V1(ik,1)*V1(ik,1) + V1(ik,2)*V1(ik,2) + V1(ik,3)*V1(ik,3)
         d(ik)     = FN(ik,1)*A0(ik,1) + FN(ik,2)*A0(ik,2) + FN(ik,3)*A0(ik,3)
      END DO
      IF (PRESENT(shared_comm)) CALL MPI_WIN_FENCE(0, win_dot00, istat)
      IF (PRESENT(shared_comm)) CALL MPI_WIN_FENCE(0, win_dot01, istat)
      IF (PRESENT(shared_comm)) CALL MPI_WIN_FENCE(0, win_dot11, istat)
      IF (PRESENT(shared_comm)) CALL MPI_WIN_FENCE(0, win_d, istat)
      IF (PRESENT(shared_comm)) CALL MPI_WIN_FENCE(0, win_ihit, istat)
      DO ik = mystart, myend
         invDenom(ik) = one / (DOT00(ik)*DOT11(ik) - DOT01(ik)*DOT01(ik))
      END DO
      IF (PRESENT(shared_comm)) CALL MPI_WIN_FENCE(0, win_invDenom, istat)
      RETURN
      END SUBROUTINE wall_load_txt

      SUBROUTINE wall_load_mn(Rmn,Zmn,xm,xn,mn,nu,nv)
      DOUBLE PRECISION, INTENT(in) :: Rmn(mn), Zmn(mn), xm(mn), xn(mn)
      INTEGER, INTENT(in) :: mn, nu, nv
      INTEGER :: u, v, i, j, istat, dex1, dex2, dex3, ik, bubble, nv2
      DOUBLE PRECISION :: pi2, th, zt, pi
      DOUBLE PRECISION, DIMENSION(3) :: temp
      DOUBLE PRECISION, ALLOCATABLE :: r_temp(:,:),z_temp(:,:),x_temp(:,:),y_temp(:,:)

      machine_string = '          HARMONICS'
      date = '      TODAY'
      pi2 = 8.0D+00 * ATAN(one)
      pi  = 4.0E+00 * ATAN(one)
      nv2 = nv/2
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

      nvertex = nu*nv
      nface   = 2*nu*nv
      istat = 0
      ALLOCATE(vertex(nvertex,3),face(nface,3),STAT=istat)
      i = 1  ! Tracks vertex index
      j = 1 ! Tracks face index
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

      ! Sort the array by toroidal angle
      ALLOCATE(PHI(nface),STAT=istat)
      IF (istat/=0) RETURN
      DO ik = 1, nface
         dex1 = face(ik,1)
         PHI(ik) = ATAN2(vertex(dex1,2),vertex(dex1,1))
      END DO
!      dex3 = nface
!      DO WHILE (dex3 > 1)
!         bubble = 0 !bubble in the greatest element out of order
!         DO ik = 1, (dex3-1)
!            IF (PHI(ik) > PHI(ik+1)) THEN
!               ! Adjust PHI
!               temp(1) = PHI(ik)
!               PHI(ik) = PHI(ik+1)
!               PHI(ik+1) = temp(1)
!               ! Adjust VERTEX
!               dex1 = face(ik,1)
!               dex2 = face(ik+1,1)
!               temp = vertex(dex1,:)
!               vertex(dex1,:) = vertex(dex2,:)
!               vertex(dex2,:) = temp
!               ! Adjust the faces
!               WHERE(face == dex1) face = 0
!               WHERE(face == dex2) face = dex1
!               WHERE(face == 0) face = dex2
!               ! Adjust bubble
!               bubble = ik
!            END IF 
!         END DO
!         dex3 = bubble
!      END DO


!      ALLOCATE(A0(nface,3),V0(nface,3),V1(nface,3),&
!               V2(nface,3),FN(nface,3),STAT=istat)
      ALLOCATE(A0(nface,3),V0(nface,3),V1(nface,3),&
               FN(nface,3),STAT=istat)
      ! Calculate the face normal
      ! W  = Vertex1-Vertex0
      ! W  = Vertex2-Vertex0
      ! FN = VxW/|VxW|
      ! d  = -Vertex0.FN (. is dot product) (not we've dropped the minus in our formulation)
      DO ik = 1, nface
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
      !DEALLOCATE(vertex,face)
!      ALLOCATE(DOT00(nface), DOT01(nface), DOT02(nface),&
!               DOT11(nface), DOT12(nface),invDenom(nface),&
!               STAT=istat)
      ALLOCATE(DOT00(nface), DOT01(nface),&
               DOT11(nface), invDenom(nface),&
               STAT=istat)
      DOT00 = SUM(V0*V0,DIM=2)
      DOT01 = SUM(V0*V1,DIM=2)
      DOT11 = SUM(V1*V1,DIM=2)
      invDenom = one/(DOT00*DOT11-DOT01*DOT01)
!      ALLOCATE(d(nface),t(nface),alpha(nface),&
!               beta(nface),lmask(nface),STAT=istat)
      ALLOCATE(d(nface),STAT=istat)
!      ALLOCATE(dr(nface,3),r0(nface,3),STAT=istat)
      d  = SUM(FN*A0,DIM=2)
!      r0 = zero
!      dr = zero
      ALLOCATE(ihit_array(nface),STAT=istat)
      ihit_array = 0
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

      SUBROUTINE collide_double_vec(x0,y0,z0,x1,y1,z1,xw,yw,zw,lhit)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x0, y0, z0, x1, y1, z1
      DOUBLE PRECISION, INTENT(out) :: xw, yw, zw
      LOGICAL, INTENT(out) :: lhit
      INTEGER :: ik, ik_min
      xw=zero; yw=zero; zw=zero; lhit=.FALSE.
      ik_min = -1
      lmask(:) = .TRUE.
!      r0(:,1) = x0
!      r0(:,2) = y0
!      r0(:,3) = z0
!      dr(:,1)=x1-x0
!      dr(:,2)=y1-y0
!      dr(:,3)=z1-z0
!      ! Calculate distance along trajectory to hit triangle
!      alpha = SUM((FN*dr),DIM=2)
!      beta  = SUM((FN*r0),DIM=2)
!      t =  (d - beta)/alpha
!      ! Check length of t
!      WHERE(alpha == 0) t=0           ! Because alpha could be zero
!      WHERE(t > one) lmask=.FALSE.    ! Longer than dr is a miss
!      WHERE(t <= zero) lmask=.FALSE.  ! negative misses triangle
!      IF (ALL(.not.lmask)) RETURN
!      ! Calculate hit point relative to triangle coordinates
!      !DO ik = 1, nface
!      !   V2(ik,:) = r0(ik,:) + (t(ik)*dr(ik,:)) - A0(ik,:)
!      !END DO
!      V2 = r0 - A0
!      V2(:,1) = V2(:,1) + t(:)*dr(:,1)
!      V2(:,2) = V2(:,2) + t(:)*dr(:,2)
!      V2(:,3) = V2(:,3) + t(:)*dr(:,3)
!      ! Calculate parametric coordinates (in plane of triangle)
!      DOT02 = SUM(V0*V2,DIM=2)
!      DOT12 = SUM(V1*V2,DIM=2)
!      alpha = ((DOT11*DOT02)-(DOT01*DOT12))*invDenom
!      beta  = ((DOT00*DOT12)-(DOT01*DOT02))*invDenom
!      ! Remove all negative values (aren't in triangle)
!      WHERE(alpha<zero) lmask = .FALSE.
!      WHERE(beta<zero)  lmask = .FALSE.
!      WHERE((alpha+beta) > 1.0001) lmask = .FALSE. ! If alpha+beta > 1 point is not in triangle
!      IF (ANY(lmask)) ik_min = minloc(t,DIM=1,MASK=lmask)
!      IF (ik_min > 0) THEN
!         lhit=.TRUE.
!         xw = V2(ik_min,1) + A0(ik_min,1)
!         yw = V2(ik_min,2) + A0(ik_min,2)
!         zw = V2(ik_min,3) + A0(ik_min,3)
!         ihit_array(ik_min) = ihit_array(ik_min) + 1
!      END IF
      RETURN
      END SUBROUTINE collide_double_vec

      SUBROUTINE collide_double(x0,y0,z0,x1,y1,z1,xw,yw,zw,lhit)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x0, y0, z0, x1, y1, z1
      DOUBLE PRECISION, INTENT(out) :: xw, yw, zw
      LOGICAL, INTENT(out) :: lhit
      INTEGER :: ik, ik_min, k1,k2
      DOUBLE PRECISION :: drx, dry, drz, V2x, V2y, V2z, DOT02l, DOT12l, tloc, tmin, alphal, betal
      xw=zero; yw=zero; zw=zero; lhit=.FALSE.
      ik_min = zero
      tmin = 2
      !LIMIT DOMAIN
      drx = ATAN2(y0,x0)
      k1 = 1; k2 = nface
      !k1  = MINLOC(ABS(PHI-drx-0.25),DIM=1)
      !k2  = MINLOC(ABS(PHI-drx+0.25),DIM=1)
      IF ((k1 == 1) .or. (k2 == nface)) THEN
         k1 = 1; k2 = nface
      END IF
      IF (k1>k2) THEN
         ik=k1; k1=k2; k2=ik
      END IF
      ! Define DR
      drx = x1-x0
      dry = y1-y0
      drz = z1-z0
      ! Calculate distance along trajectory to hit triangle
      DO ik = k1,k2
         alphal = FN(ik,1)*drx + FN(ik,2)*dry + FN(ik,3)*drz
         betal = FN(ik,1)*x0 + FN(ik,2)*y0 + FN(ik,3)*z0
         IF (alphal < zero) CYCLE
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

      SUBROUTINE wall_free(istat,shared_comm)
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: istat
      INTEGER, INTENT(inout), OPTIONAL :: shared_comm
      INTEGER :: ier
      IF (PRESENT(shared_comm)) THEN
         CALL MPI_WIN_FENCE(0,win_vertex,istat)
         CALL MPI_WIN_FREE(win_vertex,istat)
         CALL MPI_WIN_FENCE(0,win_face,istat)
         CALL MPI_WIN_FREE(win_face,istat)
         CALL MPI_WIN_FENCE(0,win_phi,istat)
         CALL MPI_WIN_FREE(win_phi,istat)
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
      ELSE
         IF (ALLOCATED(FN)) DEALLOCATE(FN)
         IF (ALLOCATED(A0)) DEALLOCATE(A0)
         IF (ALLOCATED(V0)) DEALLOCATE(V0)
         IF (ALLOCATED(V1)) DEALLOCATE(V1)
         IF (ALLOCATED(V2)) DEALLOCATE(V2)
         IF (ALLOCATED(DOT00)) DEALLOCATE(DOT00)
         IF (ALLOCATED(DOT01)) DEALLOCATE(DOT01)
         IF (ALLOCATED(DOT02)) DEALLOCATE(DOT02)
         IF (ALLOCATED(DOT11)) DEALLOCATE(DOT11)
         IF (ALLOCATED(DOT12)) DEALLOCATE(DOT12)
         IF (ALLOCATED(invDenom)) DEALLOCATE(invDenom)
         IF (ALLOCATED(d)) DEALLOCATE(d)
         IF (ALLOCATED(t)) DEALLOCATE(t)
         IF (ALLOCATED(alpha)) DEALLOCATE(alpha)
         IF (ALLOCATED(beta)) DEALLOCATE(beta)
         IF (ALLOCATED(dr)) DEALLOCATE(dr)
         IF (ALLOCATED(r0)) DEALLOCATE(r0)
         IF (ALLOCATED(lmask)) DEALLOCATE(lmask)
         IF (ALLOCATED(vertex)) DEALLOCATE(vertex)
         IF (ALLOCATED(face)) DEALLOCATE(face)
         IF (ALLOCATED(ihit_array)) DEALLOCATE(ihit_array)
         IF (ALLOCATED(PHI)) DEALLOCATE(PHI)
      END IF
      machine_string=''
      date=''
      nface = -1
      nvertex = -1
      RETURN
      END SUBROUTINE wall_free

      SUBROUTINE wall_test
      LOGICAL :: lhit
      INTEGER :: i
      DOUBLE PRECISION :: pi2,x0,y0,z0,r0,rho0,x1,y1,z1,xw,yw,zw, rr0,zz0,phi0
      r0 = 1.58   ! R0
      rho0 = 1 ! rho
      z0 = 0 ! Z0
      pi2 = (8.0 * ATAN(1.0))
      phi0 = pi2 / 3.0 ! phi0
      CALL collide(r0,DBLE(0.0),z0,r0+rho0,DBLE(0.0),z0,xw,yw,zw,lhit)
      WRITE(6,*) xw,yw,zw,lhit
      DO i = 1, 100
         x0 = r0*cos(phi0)
         y0 = r0*sin(phi0)
         rr0 = rho0*cos(pi2*(i-1)/100)
         zz0 = rho0*sin(pi2*(i-1)/100)
         x1 = (r0+rr0)*cos(phi0)
         y1 = (r0+rr0)*sin(phi0)
         z1 = z0 + zz0
         CALL collide(x0,y0,z0,x1,y1,z1,xw,yw,zw,lhit)
         WRITE(327,'(9(1X,E22.12))') x0,y0,z0,x1,y1,z1,xw,yw,zw
         CALL FLUSH(327)
      END DO
      RETURN
      END SUBROUTINE wall_test

!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE wall_mod
