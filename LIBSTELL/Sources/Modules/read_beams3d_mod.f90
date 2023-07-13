!-----------------------------------------------------------------------
!     Module:        read_beams3d_mod
!     Authors:       D. Kulla (david.kulla@ipp.mpg.de), S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          07/06/2023
!     Description:   This module stores routines for reading
!                    BEAMS3D data.
!-----------------------------------------------------------------------
MODULE read_beams3d_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
   USE stel_kinds, ONLY: rprec
   USE ez_hdf5
   USE mpi_sharmem
   USE mpi_params
   USE mpi_inc

!-----------------------------------------------------------------------
!     Module Variables
!           NTITLE:
!
!-----------------------------------------------------------------------
   IMPLICIT NONE

   ! HINT Variables
   INTEGER, PARAMETER, PRIVATE :: DTYPE =  SELECTED_REAL_KIND(15)
   INTEGER ::  nr, nz, nphi, nlines, nsteps
   DOUBLE PRECISION, PRIVATE, PARAMETER :: one           = 1.0D0 ! 1.0
   DOUBLE PRECISION, PRIVATE :: pi2, dr, dp, dz
   REAL(DTYPE), DIMENSION(:), POINTER, PRIVATE :: raxis, zaxis, phiaxis, &
      RMAGAXIS, ZMAGAXIS
   REAL(DTYPE), DIMENSION(:), POINTER, PRIVATE :: &
      R_1D, PHI_1D, Z_1D, rminor_1D, X_1D, Y_1D
   REAL(DTYPE), DIMENSION(:,:), POINTER, PRIVATE :: &
      R_lines, PHI_lines, Z_lines, rminor_lines, X_lines, Y_lines
   REAL(DTYPE), DIMENSION(:,:,:), POINTER, PRIVATE :: &
      BR3D, BPHI3D, BZ3D, X3D, Y3D, Rminor3D, U3D, S3D
   INTEGER :: win_raxis, win_phiaxis, win_zaxis, &
      win_BR3D, win_BPHI3D, win_BZ3D, win_U3D, win_S3D,&
      win_rminor_lines, &
      win_X3D, win_Y3D, win_Rminor3D, &
      win_RMAGAXIS, win_ZMAGAXIS, &
      win_R_1D, win_PHI_1D, win_Z_1D, win_rminor_1D, win_X_1D, win_Y_1D

CONTAINS

!-----------------------------------------------------------------------
!     Module Subroutines
!           read_beams3d_mag:  For reading beams3d B-Field Data
!           setup_beams3d_helpers: Setup grid helpers
!           get_beams3d_grid: Returns the grid information
!           get_beams3d_magaxis: Get magnetic axis information
!           get_beams3d_B: Get magnetic field at a point in space
!           get_beams3d_gridB: Get magnetic field at gridpoint
!           read_beams3d_deallocate: Deallocated arrays
!
!
!-----------------------------------------------------------------------

   SUBROUTINE read_beams3d_mag(filename,comm_read,istat)
      IMPLICIT NONE
      CHARACTER(*), INTENT(in) :: filename
      INTEGER, INTENT(inout) :: comm_read
      INTEGER, INTENT(inout) :: istat
      INTEGER :: mylocalid, nlocal
      LOGICAL :: lexist
      istat = 0


      ! MPI Stuff
#if defined(MPI_OPT)
      CALL MPI_BARRIER(comm_read,istat)
      CALL MPI_COMM_RANK(comm_read, mylocalid, istat)
      CALL MPI_COMM_SIZE(comm_read, nlocal, istat)
#endif

      ! Check for the file
      INQUIRE(FILE=TRIM(filename),EXIST=lexist)
      IF (.not.lexist) THEN
         istat=-1
         IF (mylocalid == 0) WRITE(6,*) " ERROR: Could not find file: "//TRIM(filename)
         RETURN
      END IF

      ! Read HDF5 File
#if defined(LHDF5)
      IF (mylocalid == master) THEN
         CALL open_hdf5(TRIM(filename),fid,istat)
         CALL read_scalar_hdf5(fid,'nr',istat,INTVAR=nr)
         CALL read_scalar_hdf5(fid,'nphi',istat,INTVAR=nphi)
         CALL read_scalar_hdf5(fid,'nz',istat,INTVAR=nz)
      END IF
#else
      IF (mylocalid == 0) WRITE(6,*) " ERROR: HDF5 Required for this functionality."
      istat = -5
      RETURN
#endif

      ! Broadcast arrays sizes and allocate arrays.
#if defined(MPI_OPT)
      CALL MPI_BARRIER(comm_read,istat)
      CALL MPI_BCAST(nr,   1, MPI_INTEGER, master, comm_read, istat)
      CALL MPI_BCAST(nphi, 1, MPI_INTEGER, master, comm_read, istat)
      CALL MPI_BCAST(nz,   1, MPI_INTEGER, master, comm_read, istat)
#endif
      CALL mpialloc(raxis,   nr,   mylocalid, master, comm_read, win_raxis)
      CALL mpialloc(phiaxis, nphi, mylocalid, master, comm_read, win_phiaxis)
      CALL mpialloc(zaxis,   nz,   mylocalid, master, comm_read, win_zaxis)
      CALL mpialloc(BR3D,   nr, nphi, nz, mylocalid, master, comm_read, win_BR3D)
      CALL mpialloc(BPHI3D, nr, nphi, nz, mylocalid, master, comm_read, win_BPHI3D)
      CALL mpialloc(BZ3D,   nr, nphi, nz, mylocalid, master, comm_read, win_BZ3D)
      CALL mpialloc(U3D, nr, nphi, nz, mylocalid, master, comm_read, win_U3D)
      CALL mpialloc(S3D,   nr, nphi, nz, mylocalid, master, comm_read, win_S3D)

      ! Read Arrays and close file
#if defined(LHDF5)
      IF (mylocalid == master) THEN
         CALL read_var_hdf5(fid, 'raxis',   nr,           istat, DBLVAR=raxis)
         CALL read_var_hdf5(fid, 'phiaxis', nphi,         istat, DBLVAR=phiaxis)
         CALL read_var_hdf5(fid, 'zaxis',   nz,           istat, DBLVAR=zaxis)
         CALL read_var_hdf5(fid, 'B_R',     nr, nphi, nz, istat, DBLVAR=BR3D)
         CALL read_var_hdf5(fid, 'B_PHI',   nr, nphi, nz, istat, DBLVAR=BPHI3D)
         CALL read_var_hdf5(fid, 'B_Z',     nr, nphi, nz, istat, DBLVAR=BZ3D)
         CALL read_var_hdf5(fid, 'S_ARR',     nr, nphi, nz, istat, DBLVAR=S3D)
         CALL read_var_hdf5(fid, 'U_ARR',   nr, nphi, nz, istat, DBLVAR=U3D)
         CALL close_hdf5(fid,istat)
      END IF
#endif
#if defined(MPI_OPT)
      CALL MPI_BARRIER(comm_read,istat)
#endif
      ! Helpers
      CALL setup_beams3d_helpers

      RETURN

   END SUBROUTINE read_beams3d_mag

   SUBROUTINE setup_beams3d_helpers
      IMPLICIT NONE
      INTEGER :: i, j, k
      INTEGER :: ij(2)
      dp = (phiaxis(nphi)-phiaxis(1))/DBLE(nphi-1)
      dr = (raxis(nr)-raxis(1))/DBLE(nr-1)
      dz = (zaxis(nz)-zaxis(1))/DBLE(nz-1)
      RETURN
   END SUBROUTINE setup_beams3d_helpers

   SUBROUTINE get_beams3d_grid(nr_out,nz_out,nphi_out,rmin_out,rmax_out,zmin_out,zmax_out,phimax_out)
      IMPLICIT NONE
      INTEGER, INTENT(out) :: nr_out, nz_out, nphi_out
      REAL(rprec), INTENT(out) :: rmin_out, rmax_out, zmin_out, zmax_out, phimax_out
      nr_out = nr; nz_out= nz; nphi_out = nphi
      rmin_out = raxis(1); zmin_out = zaxis(1)
      rmax_out = raxis(nr); zmax_out = zaxis(nz)
      phimax_out = phiaxis(nphi)
      RETURN
   END SUBROUTINE get_beams3d_grid

   SUBROUTINE get_beams3d_magaxis(phi_in,raxis_out,zaxis_out)
      IMPLICIT NONE
      REAL(rprec), INTENT(in)  :: phi_in
      REAL(rprec), INTENT(out) :: raxis_out
      REAL(rprec), INTENT(out) :: zaxis_out
      ! Define beam_data or read it from a source

      INTEGER :: phidex
      INTEGER :: position(2)
      REAL(rprec) :: smin
      REAL(rprec), ALLOCATABLE :: S2D(:,:)

      ! Find phidex for the given phi_in value
      phidex = FINDLOC(phiaxis, phi_in,DIM=1)
      ALLOCATE(S2D(nr,nz))
      S2D = S3D(:, phidex, :)
      smin = MINVAL(S2D)
      ! Find minimum value and its indices
      position = FINDLOC(S2D,smin)

      ! Assign output values
      raxis_out = raxis(position(1))
      zaxis_out = zaxis(position(2))
      DEALLOCATE(S2D)

   END SUBROUTINE get_beams3d_magaxis

   SUBROUTINE get_beams3d_B(r,phi,z,br,bp,bz,rminor,theta)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: r, phi, z
      REAL(rprec), INTENT(out) :: br, bp, bz
      REAL(rprec), INTENT(out), OPTIONAL :: rminor, theta
      INTEGER :: im, ip, km, kp, jp, jm
      REAL(rprec) :: di, dk, dj, phim, f1, f2, x, y
      br=0; bp=0; bz=0
      ! Check for bounds
      IF ((r < raxis(1)) .or. (r > raxis(nr)) .or. &
         (z < zaxis(1)) .or. (z > zaxis(nz))) RETURN
      phim = MOD(phi,phiaxis(nphi))
      IF (phim < 0) phim = phim + pi2
      ! Grid helpers
      im = COUNT(raxis < r)
      jm = COUNT(phiaxis < phim)
      km = COUNT(zaxis < z)
      im = MIN( MAX(im, 1), nr-1)
      jm = MIN( MAX(jm, 1), nphi-1)
      km = MIN( MAX(km, 1), nz-1)
      ip = im + 1
      jp = jm + 1
      kp = km + 1

      di = r - raxis(1)   - im*dr
      dj =     phiaxis(1) - jm*dp
      dk = z - zaxis(1)   - km*dz

      f1 = BR3D(im,jm,km)*(one-di)*(one-dj) &
         + BR3D(ip,jm,km)*(di)*(one-dj) &
         + BR3D(im,jp,km)*(one-di)*(dj) &
         + BR3D(ip,jp,km)*di*dj
      f2 = BR3D(im,jm,kp)*(one-di)*(one-dj) &
         + BR3D(ip,jm,kp)*(di)*(one-dj) &
         + BR3D(im,jp,kp)*(one-di)*(dj) &
         + BR3D(ip,jp,kp)*di*dj
      br = f1+ (f2-f1)*dk

      f1 = BPHI3D(im,jm,km)*(one-di)*(one-dj) &
         + BPHI3D(ip,jm,km)*(di)*(one-dj) &
         + BPHI3D(im,jp,km)*(one-di)*(dj) &
         + BPHI3D(ip,jp,km)*di*dj
      f2 = BPHI3D(im,jm,kp)*(one-di)*(one-dj) &
         + BPHI3D(ip,jm,kp)*(di)*(one-dj) &
         + BPHI3D(im,jp,kp)*(one-di)*(dj) &
         + BPHI3D(ip,jp,kp)*di*dj
      bp = f1+ (f2-f1)*dk

      f1 = BZ3D(im,jm,km)*(one-di)*(one-dj) &
         + BZ3D(ip,jm,km)*(di)*(one-dj) &
         + BZ3D(im,jp,km)*(one-di)*(dj) &
         + BZ3D(ip,jp,km)*di*dj
      f2 = BZ3D(im,jm,kp)*(one-di)*(one-dj) &
         + BZ3D(ip,jm,kp)*(di)*(one-dj) &
         + BZ3D(im,jp,kp)*(one-di)*(dj) &
         + BZ3D(ip,jp,kp)*di*dj
      bz = f1+ (f2-f1)*dk
      br = br * bp / r
      bz = bz * bp / r

      IF (PRESENT(rminor)) THEN
         rminor = 0
         f1 = Rminor3D(im,jm,km)*(one-di)*(one-dj) &
            + Rminor3D(ip,jm,km)*(di)*(one-dj) &
            + Rminor3D(im,jp,km)*(one-di)*(dj) &
            + Rminor3D(ip,jp,km)*di*dj
         f2 = Rminor3D(im,jm,kp)*(one-di)*(one-dj) &
            + Rminor3D(ip,jm,kp)*(di)*(one-dj) &
            + Rminor3D(im,jp,kp)*(one-di)*(dj) &
            + Rminor3D(ip,jp,kp)*di*dj
         rminor = f1+ (f2-f1)*dk
      END IF

      IF (PRESENT(theta)) THEN
         theta = 0
         f1 = X3D(im,jm,km)*(one-di)*(one-dj) &
            + X3D(ip,jm,km)*(di)*(one-dj) &
            + X3D(im,jp,km)*(one-di)*(dj) &
            + X3D(ip,jp,km)*di*dj
         f2 = X3D(im,jm,kp)*(one-di)*(one-dj) &
            + X3D(ip,jm,kp)*(di)*(one-dj) &
            + X3D(im,jp,kp)*(one-di)*(dj) &
            + X3D(ip,jp,kp)*di*dj
         x = f1+ (f2-f1)*dk
         f1 = Y3D(im,jm,km)*(one-di)*(one-dj) &
            + Y3D(ip,jm,km)*(di)*(one-dj) &
            + Y3D(im,jp,km)*(one-di)*(dj) &
            + Y3D(ip,jp,km)*di*dj
         f2 = Y3D(im,jm,kp)*(one-di)*(one-dj) &
            + Y3D(ip,jm,kp)*(di)*(one-dj) &
            + Y3D(im,jp,kp)*(one-di)*(dj) &
            + Y3D(ip,jp,kp)*di*dj
         y = f1+ (f2-f1)*dk
         theta = ATAN2(y,x)
      END IF

      RETURN
   END SUBROUTINE get_beams3d_B

   SUBROUTINE get_beams3d_gridB(i,j,k,br,bp,bz,rho,theta)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: i,j,k
      REAL(rprec), INTENT(out) :: br, bp, bz
      REAL(rprec), INTENT(out), OPTIONAL :: rho, theta
      br = 0; bp = 0; bz=0
      IF (i>nr .or. j>nphi .or. k>nz) RETURN
      bp = BPHI3D(i,j,k)
      br = BR3D(i,j,k)
      bz = BZ3D(i,j,k)
      IF (PRESENT(rho)) THEN
         rho = S3D(i,j,k)
      END IF
      IF (PRESENT(theta)) THEN
         theta = U3D(i,j,k)
      END IF
      RETURN
   END SUBROUTINE get_beams3d_gridB

   SUBROUTINE read_beams3d_deallocate
      IMPLICIT NONE
      IF (ASSOCIATED(raxis))        CALL mpidealloc(raxis,        win_raxis)
      IF (ASSOCIATED(phiaxis))      CALL mpidealloc(phiaxis,      win_phiaxis)
      IF (ASSOCIATED(zaxis))        CALL mpidealloc(zaxis,        win_zaxis)
      IF (ASSOCIATED(RMAGAXIS))     CALL mpidealloc(RMAGAXIS,     win_RMAGAXIS)
      IF (ASSOCIATED(ZMAGAXIS))     CALL mpidealloc(ZMAGAXIS,     win_ZMAGAXIS)
      IF (ASSOCIATED(BR3D))         CALL mpidealloc(BR3D,         win_BR3D)
      IF (ASSOCIATED(BPHI3D))       CALL mpidealloc(BPHI3D,       win_BPHI3D)
      IF (ASSOCIATED(BZ3D))         CALL mpidealloc(BZ3D,         win_BZ3D)
      IF (ASSOCIATED(X3D))          CALL mpidealloc(X3D,          win_X3D)
      IF (ASSOCIATED(Y3D))          CALL mpidealloc(Y3D,          win_Y3D)
      IF (ASSOCIATED(S3D))          CALL mpidealloc(S3D,          win_S3D)
      IF (ASSOCIATED(U3D))          CALL mpidealloc(U3D,          win_U3D)
      IF (ASSOCIATED(Rminor3D))     CALL mpidealloc(Rminor3D,     win_Rminor3D)

      ! The 2D arrays are just pointers while the 1D arrays are the actual data
      IF (ASSOCIATED(R_lines))      NULLIFY(R_lines)
      IF (ASSOCIATED(Z_lines))      NULLIFY(Z_lines)
      IF (ASSOCIATED(X_lines))      NULLIFY(X_lines)
      IF (ASSOCIATED(Y_lines))      NULLIFY(Y_lines)
      IF (ASSOCIATED(PHI_lines))    NULLIFY(PHI_lines)
      IF (ASSOCIATED(rminor_lines)) NULLIFY(rminor_lines)
      ! Deallocated globals which we don't anymore
      IF (ASSOCIATED(R_1D))      CALL mpidealloc(R_1D,      win_R_1D)
      IF (ASSOCIATED(Z_1D))      CALL mpidealloc(Z_1D,      win_Z_1D)
      IF (ASSOCIATED(PHI_1D))    CALL mpidealloc(PHI_1D,    win_PHI_1D)
      IF (ASSOCIATED(X_1D))      CALL mpidealloc(X_1D,      win_X_1D)
      IF (ASSOCIATED(Y_1D))      CALL mpidealloc(Y_1D,      win_Y_1D)
      IF (ASSOCIATED(rminor_1D)) CALL mpidealloc(rminor_1D, win_rminor_1D)
      RETURN
   END SUBROUTINE read_beams3d_deallocate

END MODULE read_beams3d_mod
