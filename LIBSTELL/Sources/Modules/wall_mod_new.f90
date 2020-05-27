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
      MODULE wall_mod_new
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE safe_open_mod
      IMPLICIT NONE

!-----------------------------------------------------------------------
!     Types
!         wall_group       Defines a bunch of wall elements
!         wall_collection  Defines a collection of wall groups
!-----------------------------------------------------------------------
      TYPE wall_element
         LOGICAL :: isshared
         INTEGER :: ntri
         INTEGER :: win_d, win_t, win_r0, win_FN, win_dr, win_A0, &
                    win_V0, win_V1, win_V2, win_DOT00, win_DOT11, &
                    win_DOT01, win_DOT02, win_DOT12, win_invDenom, &
                    win_imap
         INTEGER,          DIMENSION(:),   POINTER :: imap => null()
         DOUBLE PRECISION, DIMENSION(:),   POINTER :: d => null()
         DOUBLE PRECISION, DIMENSION(:),   POINTER :: DOT00 => null()
         DOUBLE PRECISION, DIMENSION(:),   POINTER :: DOT01 => null()
         DOUBLE PRECISION, DIMENSION(:),   POINTER :: DOT11 => null()
         DOUBLE PRECISION, DIMENSION(:),   POINTER :: invDenom => null()
         DOUBLE PRECISION, DIMENSION(:,:), POINTER :: FN => null()
         DOUBLE PRECISION, DIMENSION(:,:), POINTER :: A0 => null()
         DOUBLE PRECISION, DIMENSION(:,:), POINTER :: V0 => null()
         DOUBLE PRECISION, DIMENSION(:,:), POINTER :: V1 => null()
      END TYPE wall_element

      TYPE wall_collection
         LOGICAL :: isshared
         INTEGER :: nwall
         INTEGER :: win_xmin, win_xmax, win_ymin, win_ymax,&
                    win_zmin, win_zmax, win_sub_walls
         REAL, DIMENSION(:), POINTER :: xmin
         REAL, DIMENSION(:), POINTER :: xmax
         REAL, DIMENSION(:), POINTER :: ymin
         REAL, DIMENSION(:), POINTER :: ymax
         REAL, DIMENSION(:), POINTER :: zmin
         REAL, DIMENSION(:), POINTER :: zmax
         TYPE (wall_element) :: wall 
         TYPE (wall_collection), DIMENSION(:), POINTER :: sub_walls => null()
      END TYPE wall_collection
      
!-----------------------------------------------------------------------
!     Module Variables (PUBLIC)
!         
!-----------------------------------------------------------------------
      LOGICAL            :: lwall_loaded
      INTEGER            :: nvertex, nface, nsubdomains
      INTEGER, POINTER :: face(:,:)
      INTEGER, POINTER :: ihit_array(:)
      DOUBLE PRECISION, POINTER   :: vertex(:,:)
      CHARACTER(LEN=256) :: machine_string
      CHARACTER(LEN=256) :: date
      
!-----------------------------------------------------------------------
!     Module Variables (PRIVATE)
!         
!-----------------------------------------------------------------------

      INTEGER, PRIVATE :: win_vertex, win_face, win_ihit
      INTEGER, PRIVATE :: ik_min
      INTEGER, PRIVATE :: nmin_vtex = 512
      REAL, PRIVATE, PARAMETER                  :: deltal = 0.00
      DOUBLE PRECISION, PRIVATE, PARAMETER      :: zero = 0.0D+0
      DOUBLE PRECISION, PRIVATE, PARAMETER      :: one  = 1.0D+0

      TYPE(wall_collection), PRIVATE :: wall


      
!-----------------------------------------------------------------------
!     Subroutines
!         wall_load_txt:   Loads trianular mesh from file
!         wall_load_mn:    Creates wall from harmonics
!         wall_info:       Prints wall info.
!         wall_collide:    Calculates collision with wall
!         wall_free:       Frees module memory
!-----------------------------------------------------------------------
      INTERFACE collide
         MODULE PROCEDURE collide_double, collide_float
      END INTERFACE

      INTERFACE free_mpi_array
         MODULE PROCEDURE free_mpi_array1d_int, free_mpi_array1d_flt, free_mpi_array1d_dbl, &
                          free_mpi_array2d_int,                       free_mpi_array2d_dbl
      END INTERFACE

      INTERFACE wall_destroy
         MODULE PROCEDURE wall_destroy_element, wall_destroy_collection
      END INTERFACE

      PRIVATE :: collide_double_search
      PRIVATE :: MPI_CALC_MYRANGE
      PRIVATE :: wall_destroy_element,wall_destroy_collection
      PRIVATE :: free_mpi_array1d_int, free_mpi_array1d_dbl, free_mpi_array2d_int, free_mpi_array2d_dbl
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
      INTEGER :: shar_comm, shar_rank, shar_size, ik, iunit
      REAL    :: xmin, xmax, ymin, ymax, zmin, zmax
      
      shar_rank = 0; shar_size = 1;
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_COMM_SPLIT_TYPE(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shar_comm, istat)
         CALL MPI_COMM_RANK( shar_comm, shar_rank, istat )
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
      ELSE
#endif
         ALLOCATE(vertex(nvertex,3),face(nface,3),STAT=istat)
         IF (istat/=0) RETURN
#if defined(MPI_OPT)
      END IF
#endif
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
      IF (PRESENT(comm)) CALL MPI_BARRIER(shar_comm,istat)
#endif
      ! Now do the thing
      !CALL CALC_TREE_SIZE(nface,8,nmin_vtex)
      nsubdomains = 0
      xmin = MINVAL(vertex(:,1))! - deltal
      xmax = MAXVAL(vertex(:,1))! + deltal
      ymin = MINVAL(vertex(:,2))! - deltal
      ymax = MAXVAL(vertex(:,2))! + deltal
      zmin = MINVAL(vertex(:,3))! - deltal
      zmax = MAXVAL(vertex(:,3))! + deltal
      IF (PRESENT(comm)) THEN
         CALL wall_create_new(wall,xmin,xmax,ymin,ymax,zmin,zmax,shar_comm)
         CALL mpialloc_1d_int(ihit_array,nface,shar_rank,0,shar_comm,win_ihit)
      ELSE
         CALL wall_create_new(wall,xmin,xmax,ymin,ymax,zmin,zmax)
         ALLOCATE(ihit_array(nface),STAT=istat)
      END IF
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
      REAL    :: xmax,ymax,zmax,xmin,ymin,zmin
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
      ELSE
#endif
         ALLOCATE(vertex(nvertex,3),face(nface,3),STAT=istat)
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
      ! Now create the wall
      !CALL CALC_TREE_SIZE(nface,64,nmin_vtex)
      nsubdomains = 0
      xmin = MINVAL(vertex(:,1)) - deltal
      xmax = MAXVAL(vertex(:,1)) + deltal
      ymin = MINVAL(vertex(:,2)) - deltal
      ymax = MAXVAL(vertex(:,2)) + deltal
      zmin = MINVAL(vertex(:,3)) - deltal
      zmax = MAXVAL(vertex(:,3)) + deltal
      IF (PRESENT(comm)) THEN
         CALL wall_create_new(wall,xmin,xmax,ymin,ymax,zmin,zmax,shar_comm)
         CALL mpialloc_1d_int(ihit_array,nface,shar_rank,0,shar_comm,win_ihit)
      ELSE
         CALL wall_create_new(wall,xmin,xmax,ymin,ymax,zmin,zmax)
         ALLOCATE(ihit_array(nface),STAT=istat)
      END IF
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_BARRIER(shar_comm, istat)
         CALL MPI_COMM_FREE(shar_comm, istat)
      END IF
#endif
      lwall_loaded = .true.
      RETURN
      END SUBROUTINE wall_load_mn
      
      SUBROUTINE wall_info(iunit)
      IMPLICIT NONE
      INTEGER, INTENT(in)    :: iunit
      WRITE(iunit,'(A)')         ' -----  Vessel Information  -----'
      WRITE(iunit,'(3X,A,A)')    'Wall Name  : ',TRIM(machine_string(10:))
      WRITE(iunit,'(3X,A,A)')    'Date       : ',TRIM(date(6:))
      WRITE(iunit,'(3X,A,I7)')   'Faces      : ',nface
      WRITE(iunit,'(3X,A,I7)')   'Subdomains : ',nsubdomains
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

      SUBROUTINE collide_double(x0,y0,z0,x1,y1,z1,xw,yw,zw,lhit)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x0, y0, z0, x1, y1, z1
      DOUBLE PRECISION, INTENT(out) :: xw, yw, zw
      LOGICAL, INTENT(out) :: lhit
      DOUBLE PRECISION :: tmin
      lhit = .FALSE.
      xw = zero; yw = zero; zw = zero
      CALL collide_double_search(wall,x0,y0,z0,x1,y1,z1,xw,yw,zw,lhit,tmin,ik_min)
      IF (lhit) CALL count_wall_hit
      RETURN
      END SUBROUTINE collide_double

      RECURSIVE SUBROUTINE collide_double_search(this,x0,y0,z0,x1,y1,z1,xw,yw,zw,lhit,tmin,ik_minl)
      TYPE(wall_collection) :: this
      DOUBLE PRECISION, INTENT(in) :: x0, y0, z0, x1, y1, z1
      DOUBLE PRECISION, INTENT(out) :: xw, yw, zw, tmin
      INTEGER, INTENT(out) :: ik_minl
      LOGICAL, INTENT(out) :: lhit
      LOGICAL :: lhit2
      INTEGER :: i, ik_min2
      DOUBLE PRECISION :: xs, ys, zs, xb, yb, zb, drx, dry, drz
      DOUBLE PRECISION :: alphal, betal, tloc, DOT02, DOT12, tmin2, &
                          V2x, V2y, V2z
      ! All INTENT OUT variables need default values
      xw = zero; yw = zero; zw = zero; lhit=.FALSE.; tmin = 2; ik_minl = 0
      IF (ASSOCIATED(this%wall%FN)) THEN ! Look for collision
         drx = x1-x0
         dry = y1-y0
         drz = z1-z0
         DO i = 1, this%wall%ntri
            alphal = SUM(this%wall%FN(i,:)*(/drx,dry,drz/))
            betal  = SUM(this%wall%FN(i,:)*(/x0,y0,z0/))
            tloc   = (this%wall%d(i)-betal)/alphal
!            WRITE(6,*) '*****',i,alphal,betal,tloc
            IF (tloc > one .or. tloc <=0) CYCLE
            V2x    = x0 + tloc*drx - this%wall%A0(i,1)
            V2y    = y0 + tloc*dry - this%wall%A0(i,2)
            V2z    = z0 + tloc*drz - this%wall%A0(i,3)
            DOT02  = this%wall%V0(i,1)*V2x + this%wall%V0(i,2)*V2y + this%wall%V0(i,3)*V2z
            DOT12  = this%wall%V1(i,1)*V2x + this%wall%V1(i,2)*V2y + this%wall%V1(i,3)*V2z
            alphal = (this%wall%DOT11(i)*DOT02 - this%wall%DOT01(i)*DOT12)*this%wall%invDenom(i)
            betal  = (this%wall%DOT00(i)*DOT12 - this%wall%DOT01(i)*DOT02)*this%wall%invDenom(i)
            !WRITE(6,*) '********',i,alphal,betal,alphal+betal,tloc
            IF ((alphal < zero) .or. (betal < zero) .or. (alphal+betal > one)) CYCLE
            IF (tloc < tmin) THEN
               ik_minl = this%wall%imap(i)
               tmin = tloc
            END IF
         END DO
         IF (ik_minl > zero) THEN
            lhit = .TRUE.
            xw = x0 + tmin*drx
            yw = y0 + tmin*dry
            zw = z0 + tmin*drz
         END IF
         RETURN
      END IF
      ! Second do a recrusive search
      IF (ASSOCIATED(this%sub_walls)) THEN 
         lhit2=.FALSE.; ik_min2 = 0; tmin2 = 2
         xs = min(x0,x1)
         ys = min(y0,y1)
         zs = min(z0,z1)
         xb = max(x0,x1)
         yb = max(y0,y1)
         zb = max(z0,z1)
         DO i = 1, this%nwall
            IF (xs > this%xmax(i) .or. xb < this%xmin(i) .or. &
                ys > this%ymax(i) .or. yb < this%ymin(i) .or. &
                zs > this%zmax(i) .or. zb < this%zmin(i)) CYCLE
            !PRINT *,'Seaching ',i,this%xmax(i),this%xmin(i),this%zmax(i),this%zmin(i)
            CALL collide_double_search(this%sub_walls(i),x0,y0,z0,x1,y1,z1,xw,yw,zw,lhit2,tmin2,ik_min2)
            ! IF we have a hit determine if it's closer to x0, y0, z0 than others
            CALL FLUSH(6)
            IF (lhit2) THEN
               IF (tmin2 < tmin) THEN
                 tmin   = tmin2
                 ik_minl = ik_min2
                 lhit   = .TRUE.
               END IF
            END IF
         END DO
      END IF
      RETURN
      END SUBROUTINE collide_double_search

      SUBROUTINE count_wall_hit
         IMPLICIT NONE
         ihit_array(ik_min) = ihit_array(ik_min) + 1
      END SUBROUTINE

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
         INTEGER :: d1, d2, d3
         DOUBLE PRECISION, DIMENSION(3) :: V0, V1, FN
         d1 = face(ik,1)
         d2 = face(ik,2)
         d3 = face(ik,3)
         V0 = vertex(d3,:)-vertex(d1,:)
         V1 = vertex(d2,:)-vertex(d1,:)
         FN = (/V1(2)*V0(3)-V1(3)*V0(2), &
                V1(3)*V0(1)-V1(1)*V0(3), &
                V1(1)*V0(2)-V1(2)*V0(1)/)
         get_wall_area = 0.5*SQRT(SUM(FN*FN))
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
#if defined(MPI_OPT)
      IF (PRESENT(shared_comm)) THEN
         CALL MPI_WIN_FENCE(0,win_vertex,istat)
         CALL MPI_WIN_FREE(win_vertex,istat)
         CALL MPI_WIN_FENCE(0,win_face,istat)
         CALL MPI_WIN_FREE(win_face,istat)
         IF (ASSOCIATED(vertex)) NULLIFY(vertex)
         IF (ASSOCIATED(face)) NULLIFY(face)
      ELSE
#endif
         IF (ASSOCIATED(vertex)) DEALLOCATE(vertex)
         IF (ASSOCIATED(face)) DEALLOCATE(face)
#if defined(MPI_OPT)
      END IF
#endif
      CALL wall_destroy_collection(wall)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Class Constructors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE wall_element_subcopy(this,that,xmin,xmax,ymin,ymax,zmin,zmax,shared_comm)
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      TYPE(wall_element), INTENT(inout) :: this,that
      REAL, INTENT(in) :: xmin, ymin, zmin, xmax, ymax, zmax
      INTEGER, INTENT(inout), OPTIONAL :: shared_comm
      LOGICAL, DIMENSION(:), ALLOCATABLE :: SKIP_MASK
      INTEGER :: i, j, d1, d2, d3, ntri, mystart, myend, shar_rank, istat
      REAL :: ximin,ximax,yimin,yimax,zimin,zimax

      ! Calculate indicies to work over
      IF (PRESENT(shared_comm)) THEN
         CALL MPI_CALC_MYRANGE(shared_comm,1,this%ntri,mystart,myend)
      ELSE
         i=-1
         CALL MPI_CALC_MYRANGE(i,1,this%ntri,mystart,myend)
      END IF

      ! Allocate helper array
      ALLOCATE(SKIP_MASK(mystart:myend))
      SKIP_MASK = .TRUE.

      ! Count points in out domain
      ntri = 0
      DO i = mystart, myend
         ximin = min(this%A0(i,1),this%A0(i,1)+this%V0(i,1),this%A0(i,1)+this%V1(i,1))
         ximax = max(this%A0(i,1),this%A0(i,1)+this%V0(i,1),this%A0(i,1)+this%V1(i,1))
         yimin = min(this%A0(i,2),this%A0(i,2)+this%V0(i,2),this%A0(i,2)+this%V1(i,2))
         yimax = max(this%A0(i,2),this%A0(i,2)+this%V0(i,2),this%A0(i,2)+this%V1(i,2))
         zimin = min(this%A0(i,3),this%A0(i,3)+this%V0(i,3),this%A0(i,3)+this%V1(i,3))
         zimax = max(this%A0(i,3),this%A0(i,3)+this%V0(i,3),this%A0(i,3)+this%V1(i,3))
         IF (ximax < xmin .or. ximin > xmax .or. &
             yimax < ymin .or. yimin > ymax .or. &
             zimax < zmin .or. zimin > zmax) CYCLE
         ntri = ntri + 1
         SKIP_MASK(i) = .FALSE.
      END DO

      ! Set number of triangles in subdomain
      IF (PRESENT(shared_comm)) THEN
#if defined(MPI_OPT)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,ntri,1,MPI_INTEGER,MPI_SUM,shared_comm,istat)
#endif
      END IF
      that%ntri = ntri

      ! If there were no triangles then get out of here
      IF (ntri < 1) THEN
         DEALLOCATE(SKIP_MASK)
         RETURN
      END IF

      ! Allocate the arrays
#if defined(MPI_OPT)
      IF (PRESENT(shared_comm)) THEN
         CALL MPI_COMM_RANK( shared_comm, shar_rank, istat )
         CALL mpialloc_1d_int(that%imap,     ntri, shar_rank, 0, shared_comm, that%win_imap)
         CALL mpialloc_1d_dbl(that%DOT00,    ntri, shar_rank, 0, shared_comm, that%win_DOT00)
         CALL mpialloc_1d_dbl(that%DOT01,    ntri, shar_rank, 0, shared_comm, that%win_DOT01)
         CALL mpialloc_1d_dbl(that%DOT11,    ntri, shar_rank, 0, shared_comm, that%win_DOT11)
         CALL mpialloc_1d_dbl(that%d,        ntri, shar_rank, 0, shared_comm, that%win_d)
         CALL mpialloc_1d_dbl(that%invDenom, ntri, shar_rank, 0, shared_comm, that%win_invDenom)
         CALL mpialloc_2d_dbl(that%A0, ntri, 3, shar_rank, 0, shared_comm, that%win_A0)
         CALL mpialloc_2d_dbl(that%V0, ntri, 3, shar_rank, 0, shared_comm, that%win_V0)
         CALL mpialloc_2d_dbl(that%V1, ntri, 3, shar_rank, 0, shared_comm, that%win_V1)
         CALL mpialloc_2d_dbl(that%FN, ntri, 3, shar_rank, 0, shared_comm, that%win_FN)
         that%isshared = .true.
      ELSE
#endif
         ALLOCATE(that%DOT00(ntri),that%DOT01(ntri),that%DOT11(ntri),that%invDenom(ntri), &
                  that%d(ntri),that%imap(ntri))
         ALLOCATE(that%A0(ntri,3),that%V0(ntri,3),that%V1(ntri,3),that%FN(ntri,3))
         that%isshared = .false.
#if defined(MPI_OPT)
      END IF
#endif

      ! Now copy the important arrays
      ntri = 1
      DO i = mystart, myend
         IF (SKIP_MASK(i)) CYCLE
         that%imap(ntri)     = this%imap(i)
         that%DOT00(ntri)    = this%DOT00(i)
         that%DOT01(ntri)    = this%DOT01(i)
         that%DOT11(ntri)    = this%DOT11(i)
         that%d(ntri)        = this%d(i)
         that%invDenom(ntri) = this%invDenom(i)
         that%A0(ntri,:)     = this%A0(i,:)
         that%V0(ntri,:)     = this%V0(i,:)
         that%V1(ntri,:)     = this%V1(i,:)
         that%FN(ntri,:)     = this%FN(i,:)
         ntri = ntri + 1
      END DO

      ! Deallocate helper
      DEALLOCATE(SKIP_MASK)

      RETURN
      END SUBROUTINE wall_element_subcopy


      SUBROUTINE wall_element_init(this,xmin,xmax,ymin,ymax,zmin,zmax,shared_comm)
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      TYPE(wall_element), INTENT(inout) :: this
      REAL, INTENT(in) :: xmin, ymin, zmin, xmax, ymax, zmax
      INTEGER, INTENT(inout), OPTIONAL :: shared_comm
      LOGICAL, DIMENSION(:), ALLOCATABLE :: SKIP_MASK
      INTEGER :: i, d1, d2, d3, ntri, mystart, myend, shar_rank, istat
      REAL :: ximin,ximax,yimin,yimax,zimin,zimax
      ntri = 0

      ! Calculate indicies to work over
      IF (PRESENT(shared_comm)) THEN
         CALL MPI_CALC_MYRANGE(shared_comm,1,nface,mystart,myend)
      ELSE
         i=-1
         CALL MPI_CALC_MYRANGE(i,1,nface,mystart,myend)
      END IF

      ! Allocate Helper Array
      ALLOCATE(SKIP_MASK(mystart:myend))
      SKIP_MASK = .TRUE.

      ! Calculate subdomain 
      DO i = mystart, myend
         d1 = face(i,1)
         d2 = face(i,2)
         d3 = face(i,3)
         ximin = min(vertex(d1,1),vertex(d2,1),vertex(d3,1))
         ximax = max(vertex(d1,1),vertex(d2,1),vertex(d3,1))
         yimin = min(vertex(d1,2),vertex(d2,2),vertex(d3,2))
         yimax = max(vertex(d1,2),vertex(d2,2),vertex(d3,2))
         zimin = min(vertex(d1,3),vertex(d2,3),vertex(d3,3))
         zimax = max(vertex(d1,3),vertex(d2,3),vertex(d3,3))
         IF (ximax < xmin .or. ximin > xmax .or. &
             yimax < ymin .or. yimin > ymax .or. &
             zimax < zmin .or. zimin > zmax) CYCLE
         SKIP_MASK = .FALSE.
         ntri = ntri + 1
      END DO

      ! Set number of triangles
      IF (PRESENT(shared_comm)) THEN
#if defined(MPI_OPT)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,ntri,1,MPI_INTEGER,MPI_SUM,shared_comm,istat)
#endif
      END IF
      this%ntri = ntri

      ! Allocate data arrays
#if defined(MPI_OPT)
      IF (PRESENT(shared_comm)) THEN
         CALL MPI_COMM_RANK( shared_comm, shar_rank, istat )
         CALL mpialloc_1d_int(this%imap,     ntri, shar_rank, 0, shared_comm, this%win_imap)
         CALL mpialloc_1d_dbl(this%DOT00,    ntri, shar_rank, 0, shared_comm, this%win_DOT00)
         CALL mpialloc_1d_dbl(this%DOT01,    ntri, shar_rank, 0, shared_comm, this%win_DOT01)
         CALL mpialloc_1d_dbl(this%DOT11,    ntri, shar_rank, 0, shared_comm, this%win_DOT11)
         CALL mpialloc_1d_dbl(this%d,        ntri, shar_rank, 0, shared_comm, this%win_d)
         CALL mpialloc_1d_dbl(this%invDenom, ntri, shar_rank, 0, shared_comm, this%win_invDenom)
         CALL mpialloc_2d_dbl(this%A0, ntri, 3, shar_rank, 0, shared_comm, this%win_A0)
         CALL mpialloc_2d_dbl(this%V0, ntri, 3, shar_rank, 0, shared_comm, this%win_V0)
         CALL mpialloc_2d_dbl(this%V1, ntri, 3, shar_rank, 0, shared_comm, this%win_V1)
         CALL mpialloc_2d_dbl(this%FN, ntri, 3, shar_rank, 0, shared_comm, this%win_FN)
         this%isshared = .true.
      ELSE
#endif
         ALLOCATE(this%DOT00(ntri),this%DOT01(ntri),this%DOT11(ntri),this%invDenom(ntri), &
                  this%d(ntri),this%imap(ntri))
         ALLOCATE(this%A0(ntri,3),this%V0(ntri,3),this%V1(ntri,3),this%FN(ntri,3))
         this%isshared = .false.
#if defined(MPI_OPT)
      END IF
#endif

      ! Calculate Basis vectors in subdomain
      ntri = 1
      DO i = mystart, myend
         IF (SKIP_MASK(i)) CYCLE
         d1 = face(i,1)
         d2 = face(i,2)
         d3 = face(i,3)
         this%imap(ntri) = i
         this%A0(ntri,:) = vertex(d1,:)
         this%V0(ntri,:)  = vertex(d3,:)-vertex(d1,:)
         this%V1(ntri,:)  = vertex(d2,:)-vertex(d1,:)
         ntri = ntri + 1
      END DO
      DEALLOCATE(SKIP_MASK)

      ! Calculate working dimensions over subdomain triangles
      IF (PRESENT(shared_comm)) THEN
         CALL MPI_CALC_MYRANGE(shared_comm,1,this%ntri,mystart,myend)
      ELSE
         i=-1
         CALL MPI_CALC_MYRANGE(i,1,this%ntri,mystart,myend)
      END IF

      ! Calculate data arrays
      DO i = mystart, myend
         this%FN(i,1) = this%V1(i,2)*this%V0(i,3)-this%V1(i,3)*this%V0(i,2)
         this%FN(i,2) = this%V1(i,3)*this%V0(i,1)-this%V1(i,1)*this%V0(i,3)
         this%FN(i,3) = this%V1(i,1)*this%V0(i,2)-this%V1(i,2)*this%V0(i,1)
         this%DOT00(i) = this%V0(i,1) * this%V0(i,1) + &
                         this%V0(i,2) * this%V0(i,2) + &
                         this%V0(i,3) * this%V0(i,3)
         this%DOT01(i) = this%V0(i,1) * this%V1(i,1) + &
                         this%V0(i,2) * this%V1(i,2) + &
                         this%V0(i,3) * this%V1(i,3)
         this%DOT11(i) = this%V1(i,1) * this%V1(i,1) + &
                         this%V1(i,2) * this%V1(i,2) + &
                         this%V1(i,3) * this%V1(i,3)
         this%d(i)     = this%FN(i,1) * this%A0(i,1) + &
                         this%FN(i,2) * this%A0(i,2) + &
                         this%FN(i,3) * this%A0(i,3)
         this%invDenom(i) = one / (this%DOT00(i)*this%DOT11(i) - &
                                       this%DOT01(i)*this%DOT01(i))
      END DO
      nsubdomains = nsubdomains+1
      RETURN
      END SUBROUTINE wall_element_init

      RECURSIVE SUBROUTINE wall_create(this,xmin,xmax,ymin,ymax,zmin,zmax,shared_comm)
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      TYPE(wall_collection), INTENT(inout) :: this
      REAL, INTENT(in) :: xmin, ymin, zmin, xmax, ymax, zmax
      INTEGER, INTENT(inout), OPTIONAL :: shared_comm
      INTEGER :: i, nlocal, istat, mystart, myend, shar_rank
      REAL :: xmean, ymean, zmean

      nlocal = COUNT(vertex(:,1) .ge. xmin .and. vertex(:,1) .lt. xmax .and. &
         vertex(:,2) .ge. xmin .and. vertex(:,2) .lt. xmax .and. &
         vertex(:,3) .ge. xmin .and. vertex(:,3) .lt. xmax)

      IF (nlocal < nmin_vtex) THEN ! Create a wall_element
         this%nwall=1
         IF (PRESENT(shared_comm)) THEN
            CALL wall_element_init(this%wall,xmin,xmax,ymin,ymax,zmin,zmax,shared_comm)
         ELSE
            CALL wall_element_init(this%wall,xmin,xmax,ymin,ymax,zmin,zmax)
         END IF
      ELSE ! Divide up the space
         xmean = 0; ymean=0; zmean=0;

         DO i = 1, nvertex
            IF (vertex(i,1)<xmin .or. vertex(i,1)>xmax .or. &
                vertex(i,2)<ymin .or. vertex(i,2)>ymax .or. &
                vertex(i,3)<zmin .or. vertex(i,3)>zmax) CYCLE
            xmean = xmean + vertex(i,1)
            ymean = ymean + vertex(i,2)
            zmean = zmean + vertex(i,3)
         END DO

         xmean = xmean/nlocal
         ymean = ymean/nlocal
         zmean = zmean/nlocal
         ALLOCATE(this%sub_walls(8))
#if defined(MPI_OPT)
         IF (PRESENT(shared_comm)) THEN
            CALL MPI_COMM_RANK( shared_comm, shar_rank, istat )
            CALL mpialloc_1d_flt(this%xmin, 8, shar_rank, 0, shared_comm, this%win_xmin)
            CALL mpialloc_1d_flt(this%ymin, 8, shar_rank, 0, shared_comm, this%win_ymin)
            CALL mpialloc_1d_flt(this%zmin, 8, shar_rank, 0, shared_comm, this%win_zmin)
            CALL mpialloc_1d_flt(this%xmax, 8, shar_rank, 0, shared_comm, this%win_xmax)
            CALL mpialloc_1d_flt(this%ymax, 8, shar_rank, 0, shared_comm, this%win_ymax)
            CALL mpialloc_1d_flt(this%zmax, 8, shar_rank, 0, shared_comm, this%win_zmax)
            this % isshared = .TRUE.
         ELSE
#endif
            shar_rank = 0
            ALLOCATE(this%xmin(8), this%xmax(8), &
                     this%ymin(8), this%ymax(8), &
                     this%zmin(8), this%zmax(8))
#if defined(MPI_OPT)
         END IF
#endif
         IF (shar_rank == 0) THEN
            this%xmin = (/xmin,  xmin,  xmin,  xmin,  xmean, xmean, xmean, xmean/)
            this%xmax = (/xmean, xmean, xmean, xmean, xmax,  xmax,  xmax,  xmax/)
            this%ymin = (/ymin,  ymin,  ymean, ymean, ymin,  ymin,  ymean, ymean/)
            this%ymax = (/ymean, ymean, ymax,  ymax,  ymean, ymean, ymax,  ymax/)
            this%zmin = (/zmin,  zmean, zmin,  zmean, zmin,  zmean, zmin,  zmean/)
            this%zmax = (/zmean, zmax,  zmean, zmax,  zmean, zmax,  zmean, zmax/)
         END IF
         DO i = 1, 8
            IF (PRESENT(shared_comm)) THEN
               CALL wall_create(this%sub_walls(i),this%xmin(i),this%xmax(i), &
                                               this%ymin(i),this%ymax(i), &
                                               this%zmin(i),this%zmax(i),shared_comm)
            ELSE
               CALL wall_create(this%sub_walls(i),this%xmin(i),this%xmax(i), &
                                               this%ymin(i),this%ymax(i), &
                                               this%zmin(i),this%zmax(i))
            END IF
         END DO
         this%nwall = 8
      END IF
      RETURN
      END SUBROUTINE wall_create


      RECURSIVE SUBROUTINE wall_create_new(this,xmin,xmax,ymin,ymax,zmin,zmax,shared_comm)
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      TYPE(wall_collection), INTENT(inout) :: this
      REAL, INTENT(in) :: xmin, ymin, zmin, xmax, ymax, zmax
      INTEGER, INTENT(inout), OPTIONAL :: shared_comm
      INTEGER :: i, j, nlocal, istat, mystart, myend, shar_rank,ier
      REAL :: xmean, ymean, zmean, xtemp, ytemp, ztemp, deltax, deltay, deltaz

      ier = 0
      ! If this is the first level we need to create the wall otherwise it's a subcopy
      IF (.not.ASSOCIATED(this%wall%FN)) THEN
         CALL FLUSH(6)
         IF (PRESENT(shared_comm)) THEN
            CALL wall_element_init(this%wall,xmin,xmax,ymin,ymax,zmin,zmax,shared_comm)
         ELSE
            CALL wall_element_init(this%wall,xmin,xmax,ymin,ymax,zmin,zmax)
         END IF
         nsubdomains = nsubdomains - 1
      END IF

      IF(this%wall%ntri < nmin_vtex) THEN
         nsubdomains = nsubdomains + 1
         RETURN  ! Stop creating subtrees
      END IF

      ! Calculate the mean values Note all triangles are considered
      xmean = 0; ymean=0; zmean=0
      DO i = 1, this%wall%ntri
         xtemp = 3*this%wall%A0(i,1)+this%wall%V0(i,1)+this%wall%V1(i,1)
         ytemp = 3*this%wall%A0(i,2)+this%wall%V0(i,2)+this%wall%V1(i,2)
         ztemp = 3*this%wall%A0(i,3)+this%wall%V0(i,3)+this%wall%V1(i,3)
         xmean = xmean + xtemp
         ymean = ymean + ytemp
         zmean = zmean + ztemp
      END DO
      xmean = xmean / (3*this%wall%ntri)
      ymean = ymean / (3*this%wall%ntri)
      zmean = zmean / (3*this%wall%ntri)

      ! Now Allocate the Subdomains
      ALLOCATE(this%sub_walls(8))
#if defined(MPI_OPT)
      IF (PRESENT(shared_comm)) THEN
         CALL MPI_COMM_RANK( shared_comm, shar_rank, istat )
         CALL mpialloc_1d_flt(this%xmin, 8, shar_rank, 0, shared_comm, this%win_xmin)
         CALL mpialloc_1d_flt(this%ymin, 8, shar_rank, 0, shared_comm, this%win_ymin)
         CALL mpialloc_1d_flt(this%zmin, 8, shar_rank, 0, shared_comm, this%win_zmin)
         CALL mpialloc_1d_flt(this%xmax, 8, shar_rank, 0, shared_comm, this%win_xmax)
         CALL mpialloc_1d_flt(this%ymax, 8, shar_rank, 0, shared_comm, this%win_ymax)
         CALL mpialloc_1d_flt(this%zmax, 8, shar_rank, 0, shared_comm, this%win_zmax)
         this % isshared = .TRUE.
      ELSE
#endif
         shar_rank = 0
         ALLOCATE(this%xmin(8), this%xmax(8), &
                  this%ymin(8), this%ymax(8), &
                  this%zmin(8), this%zmax(8))
#if defined(MPI_OPT)
      END IF
#endif

      ! Define the subdomains
      !deltax = MIN(xmean-xmin,xmax-xmean)*deltal
      !deltay = MIN(ymean-ymin,ymax-ymean)*deltal
      !deltaz = MIN(zmean-zmin,zmax-zmean)*deltal
      IF (shar_rank == 0) THEN
         PRINT *,deltax,deltay,deltaz
         this%xmin = (/xmin,  xmin,  xmin,  xmin,  xmean, xmean, xmean, xmean/) - deltal
         this%xmax = (/xmean, xmean, xmean, xmean, xmax,  xmax,  xmax,  xmax/)  + deltal
         this%ymin = (/ymin,  ymin,  ymean, ymean, ymin,  ymin,  ymean, ymean/) - deltal
         this%ymax = (/ymean, ymean, ymax,  ymax,  ymean, ymean, ymax,  ymax/)  + deltal
         this%zmin = (/zmin,  zmean, zmin,  zmean, zmin,  zmean, zmin,  zmean/) - deltal
         this%zmax = (/zmean, zmax,  zmean, zmax,  zmean, zmax,  zmean, zmax/)  + deltal
      END IF

      ! Now fill subdomains (then recursively call)
      ! must do all 8 to keep arrays linked
      this%nwall = 8
      DO i = 1, this%nwall
         IF (PRESENT(shared_comm)) THEN
            CALL wall_element_subcopy(this%wall,this%sub_walls(i)%wall,this%xmin(i),this%xmax(i), &
                                                                       this%ymin(i),this%ymax(i), &
                                                                       this%zmin(i),this%zmax(i),shared_comm)
            IF (this%sub_walls(i)%wall%ntri == 0) CYCLE 
            CALL wall_create_new(this%sub_walls(i),this%xmin(i),this%xmax(i), &
                                                   this%ymin(i),this%ymax(i), &
                                                   this%zmin(i),this%zmax(i), shared_comm)
         ELSE
            CALL wall_element_subcopy(this%wall,this%sub_walls(i)%wall,this%xmin(i),this%xmax(i), &
                                                                       this%ymin(i),this%ymax(i), &
                                                                       this%zmin(i),this%zmax(i))
            IF (this%sub_walls(i)%wall%ntri == 0) CYCLE 
            CALL wall_create_new(this%sub_walls(i),this%xmin(i),this%xmax(i), &
                                                   this%ymin(i),this%ymax(i), &
                                                   this%zmin(i),this%zmax(i))
         END IF
      END DO
      this%nwall = 8

      ! Now delete the wall at this level
#if defined(MPI_OPT)
      IF (PRESENT(shared_comm)) CALL MPI_BARRIER(shared_comm,ier)
#endif
      CALL wall_destroy_element(this%wall)

      RETURN
      END SUBROUTINE wall_create_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Class Destructors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE wall_destroy_element(this)
      IMPLICIT NONE
      TYPE(wall_element), INTENT(inout) :: this
      this % ntri = 0
      CALL free_mpi_array(this % win_imap, this % imap, this % isshared)
      CALL free_mpi_array(this % win_d, this % d, this % isshared)
      CALL free_mpi_array(this % win_FN, this % FN, this % isshared)
      CALL free_mpi_array(this % win_A0, this % A0, this % isshared)
      CALL free_mpi_array(this % win_V0, this % V0, this % isshared)
      CALL free_mpi_array(this % win_V1, this % V1, this % isshared)
      CALL free_mpi_array(this % win_DOT00, this % DOT00, this % isshared)
      CALL free_mpi_array(this % win_DOT01, this % DOT01, this % isshared)
      CALL free_mpi_array(this % win_DOT11, this % DOT11, this % isshared)
      CALL free_mpi_array(this % win_invDenom, this % invDenom, this % isshared)
      this % isshared = .FALSE.
      RETURN
      END SUBROUTINE wall_destroy_element

      RECURSIVE SUBROUTINE wall_destroy_collection(this)
      IMPLICIT NONE
      TYPE(wall_collection), INTENT(inout) :: this
      INTEGER :: i
      IF (ASSOCIATED(this%wall%FN)) THEN
         CALL wall_destroy(this%wall)
      END IF
      IF (ASSOCIATED(this%sub_walls)) THEN
         DO i = 1, this%nwall
            CALL wall_destroy(this%sub_walls(i))
         END DO
         CALL free_mpi_array(this % win_xmin, this % xmin, this % isshared)
         CALL free_mpi_array(this % win_ymin, this % ymin, this % isshared)
         CALL free_mpi_array(this % win_zmin, this % zmin, this % isshared)
         CALL free_mpi_array(this % win_xmax, this % xmax, this % isshared)
         CALL free_mpi_array(this % win_ymax, this % ymax, this % isshared)
         CALL free_mpi_array(this % win_zmax, this % zmax, this % isshared)
         DEALLOCATE(this%sub_walls)
      END IF
      this % nwall = -1
      this % isshared = .FALSE.
      RETURN
      END SUBROUTINE wall_destroy_collection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Memory Allocation Subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

      SUBROUTINE mpialloc_1d_flt(array,n1,subid,mymaster,share_comm,win)
      ! Libraries
#if defined(MPI_OPT)
      USE mpi
#endif
      USE ISO_C_BINDING
      IMPLICIT NONE
      ! Arguments
      REAL, POINTER, INTENT(inout) :: array(:)
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
      END SUBROUTINE mpialloc_1d_flt

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Utility Subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE MPI_CALC_MYRANGE(comm,n1,n2,mystart,myend)
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: comm
      INTEGER, INTENT(in) :: n1, n2
      INTEGER, INTENT(out) :: mystart, myend
      INTEGER :: delta, local_size, local_rank, istat
      mystart = n1; myend = n2
      IF (COMM < 0) RETURN
#if defined(MPI_OPT)
      CALL MPI_COMM_SIZE( comm, local_size, istat)
      CALL MPI_COMM_RANK( comm, local_rank, istat )
      delta = CEILING(REAL(n2-n1+1)/REAL(local_size))
      mystart = n1 + local_rank*delta
      myend   = mystart + delta
      IF (myend > n2) myend=n2
#endif
      RETURN
      END SUBROUTINE MPI_CALC_MYRANGE

      SUBROUTINE CALC_TREE_SIZE(ntotal,ntree,nout)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ntotal, ntree
      INTEGER, INTENT(out)  :: nout
      nout = 1024
      nout = ntotal/ntree
      nout = 8 ** FLOOR(log(REAL(nout))/log(8.0D+00))
      RETURN
      END SUBROUTINE CALC_TREE_SIZE

!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE wall_mod_new
