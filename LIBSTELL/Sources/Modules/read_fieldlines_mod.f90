!-----------------------------------------------------------------------
!     Module:        read_fieldlines_mod
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          05/30/2023
!     Description:   This module stores routines for reading
!                    FIELDLINES data.
!-----------------------------------------------------------------------
MODULE read_fieldlines_mod
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
      BR3D, BPHI3D, BZ3D, X3D, Y3D, Rminor3D
   INTEGER :: win_raxis, win_phiaxis, win_zaxis, &
      win_BR3D, win_BPHI3D, win_BZ3D, &
      win_R_lines, win_PHI_lines, win_Z_lines, &
      win_rminor_lines, win_X_lines, win_Y_lines, &
      win_X3D, win_Y3D, win_Rminor3D, &
      win_RMAGAXIS, win_ZMAGAXIS, &
      win_R_1D, win_PHI_1D, win_Z_1D, win_rminor_1D, win_X_1D, win_Y_1D

CONTAINS

!-----------------------------------------------------------------------
!     Module Subroutines
!           read_fieldlines_mag:  For reading fieldlines B-Field Data
!           setup_fieldlines_rhogrid:  For setting up rho grid
!           setup_fieldlines_helpers: Setup grid helpers
!           get_fieldlines_grid: Returns the grid information
!           get_fieldlines_magaxis: Get magnetic axis information
!           get_fieldlines_B: Get magnetic field at a point in space
!           get_fieldlines_gridB: Get magnetic field at gridpoint
!           read_fieldlines_deallocate: Deallocated arrays
!
!
!-----------------------------------------------------------------------

   SUBROUTINE read_fieldlines_mag(filename,comm_read,istat)
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

      ! Read Arrays and close file
#if defined(LHDF5)
      IF (mylocalid == master) THEN
         CALL read_var_hdf5(fid, 'raxis',   nr,           istat, DBLVAR=raxis)
         CALL read_var_hdf5(fid, 'phiaxis', nphi,         istat, DBLVAR=phiaxis)
         CALL read_var_hdf5(fid, 'zaxis',   nz,           istat, DBLVAR=zaxis)
         CALL read_var_hdf5(fid, 'B_R',     nr, nphi, nz, istat, DBLVAR=BR3D)
         CALL read_var_hdf5(fid, 'B_PHI',   nr, nphi, nz, istat, DBLVAR=BPHI3D)
         CALL read_var_hdf5(fid, 'B_Z',     nr, nphi, nz, istat, DBLVAR=BZ3D)
         CALL close_hdf5(fid,istat)
      END IF
#endif
#if defined(MPI_OPT)
      CALL MPI_BARRIER(comm_read,istat)
#endif
      ! Helpers
      CALL setup_fieldlines_helpers

      RETURN

   END SUBROUTINE read_fieldlines_mag

   SUBROUTINE setup_fieldlines_rhogrid(filename,comm_read,istat)
      IMPLICIT NONE
      CHARACTER(*), INTENT(in) :: filename
      INTEGER, INTENT(inout) :: comm_read
      INTEGER, INTENT(out) :: istat
      INTEGER :: mylocalid, nlocal, s, i, j, k, mystart, myend, npoinc, l, m
      REAL(DTYPE) :: r1,r2,p1,p2,z1,z2
      LOGICAL, DIMENSION(:), ALLOCATABLE :: mask_axis
      REAL(DTYPE), DIMENSION(:), ALLOCATABLE :: R_help, Z_help
      REAL(DTYPE), DIMENSION(:,:), POINTER :: zeta_lines
      REAL(DTYPE), DIMENSION(:), POINTER :: ZETA_1D
      INTEGER, DIMENSION(:,:,:), POINTER :: N3D
      INTEGER :: win_zeta_lines, win_N3D
      INTEGER :: win_zeta_1D
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
         CALL read_scalar_hdf5(fid,'nlines',istat,INTVAR=nlines)
         CALL read_scalar_hdf5(fid,'nsteps',istat,INTVAR=nsteps)
         CALL read_scalar_hdf5(fid,'npoinc',istat,INTVAR=npoinc)
      END IF
#else
      IF (mylocalid == 0) WRITE(6,*) " ERROR: HDF5 Required for this functionality."
      istat = -5
      RETURN
#endif

      ! Broadcast arrays sizes and allocate arrays.
      ! Note that in FIELDLINES _lines are sized (1:nlines,0:nsteps)
      ! But FIELDLINES saves nsteps as nsteps+1
      ! so we ressize them to (1:nlines,1:nsteps)
#if defined(MPI_OPT)
      CALL MPI_BARRIER(comm_read,istat)
      CALL MPI_BCAST(nlines, 1, MPI_INTEGER, master, comm_read, istat)
      CALL MPI_BCAST(nsteps, 1, MPI_INTEGER, master, comm_read, istat)
      CALL MPI_BCAST(npoinc, 1, MPI_INTEGER, master, comm_read, istat)
#endif
      CALL mpialloc(RMAGAXIS,  nphi,          mylocalid, master, comm_read, win_RMAGAXIS)
      CALL mpialloc(ZMAGAXIS,  nphi,          mylocalid, master, comm_read, win_ZMAGAXIS)
      CALL mpialloc(R_1D,      nlines*nsteps, mylocalid, master, comm_read, win_R_1D)
      CALL mpialloc(PHI_1D,    nlines*nsteps, mylocalid, master, comm_read, win_PHI_1D)
      CALL mpialloc(Z_1D,      nlines*nsteps, mylocalid, master, comm_read, win_Z_1D)
      CALL mpialloc(X_1D,      nlines*nsteps, mylocalid, master, comm_read, win_X_1D)
      CALL mpialloc(Y_1D,      nlines*nsteps, mylocalid, master, comm_read, win_Y_1D)
      CALL mpialloc(rminor_1D, nlines*nsteps, mylocalid, master, comm_read, win_rminor_1D)
      CALL mpialloc(ZETA_1D,   nlines*nsteps, mylocalid, master, comm_read, win_ZETA_1D)

      ! Create 2D Pointers
      R_lines(1:nlines,1:nsteps) => R_1D
      Z_lines(1:nlines,1:nsteps) => Z_1D
      X_lines(1:nlines,1:nsteps) => X_1D
      Y_lines(1:nlines,1:nsteps) => Y_1D
      PHI_lines(1:nlines,1:nsteps) => PHI_1D
      rminor_lines(1:nlines,1:nsteps) => rminor_1D
      zeta_lines(1:nlines,1:nsteps) => ZETA_1D

      ! Read Arrays and close file
#if defined(LHDF5)
      IF (mylocalid == master) THEN
         CALL read_var_hdf5(fid, 'R_lines',   nlines*nsteps, istat, DBLVAR=R_1D)
         CALL read_var_hdf5(fid, 'PHI_lines', nlines*nsteps, istat, DBLVAR=PHI_1D)
         CALL read_var_hdf5(fid, 'Z_lines',   nlines*nsteps, istat, DBLVAR=Z_1D)
         CALL close_hdf5(fid,istat)
         ZETA_1D = MODULO(PHI_1D,phiaxis(nphi))
      END IF
#endif
#if defined(MPI_OPT)
      CALL MPI_BARRIER(comm_read,istat)
#endif

      ! Create rminor helper
      IF (mylocalid == master) THEN
         ALLOCATE(R_help(nsteps),Z_help(nsteps))
         DO i = 1, nlines
            ! Should use PHI_lines to figure things out
            IF (PHI_lines(i,1) /= 0.0) THEN
               j = COUNT(ABS(PHI_lines(i,:)) <= phiaxis(nphi))
               l = nsteps - j + 1
               m = nsteps - j + 2
               R_help(1:l)      = R_lines(i,j:nsteps)
               R_help(m:nsteps) = R_lines(i,1:j-1)
               Z_help(1:l)      = Z_lines(i,j:nsteps)
               Z_help(m:nsteps) = Z_lines(i,1:j-1)
            ELSE
               R_help = R_lines(i,:)
               Z_help = Z_lines(i,:)
            END IF
            X_lines(i,:) = R_help - R_lines(1,:)
            Y_lines(i,:) = Z_help - Z_lines(1,:)
         END DO
         rminor_lines = SQRT(X_lines**2 + Y_lines**2)
         DO i = 1, nlines
            rminor_lines(i,:) = SUM(rminor_lines(i,:))/(nsteps)
         END DO
         DEALLOCATE(R_help,Z_help)
      END IF

      ! Check we've read the background grids
      IF (nr <= 0 .OR.  nphi <= 0 .OR. nz <= 0) THEN
         IF (mylocalid==master) WRITE(6,*) "ERROR: NR,NPHI,NZ <= 0 in read_fieldines_lines"
         istat = -6
         RETURN
      END IF
      CALL setup_fieldlines_helpers

      ! Allocate rho grid
      CALL mpialloc(X3D,      nr, nphi, nz, mylocalid, master, comm_read, win_X3D)
      CALL mpialloc(Y3D,      nr, nphi, nz, mylocalid, master, comm_read, win_Y3D)
      CALL mpialloc(Rminor3D, nr, nphi, nz, mylocalid, master, comm_read, win_Rminor3D)
      CALL mpialloc(N3D,      nr, nphi, nz, mylocalid, master, comm_read, win_N3D)

      ! Allocate helper
      ALLOCATE(mask_axis(nsteps))

      ! Default the rho grid to the largest value of Rminor
      IF (mylocalid==master) Rminor3D = MAXVAL(rminor_lines)
#if defined(MPI_OPT)
      CALL MPI_BARRIER(comm_read,istat)
#endif

      ! Create the background helpers
      CALL MPI_CALC_MYRANGE(comm_read,1,nlines*nsteps,mystart,myend)
      DO s = mystart, myend
         i = COUNT((raxis-0.5*dr)   <= R_1D(s))
         j = COUNT((phiaxis-0.5*dp) <= zeta_1D(s))
         k = COUNT((zaxis-0.5*dz)   <= Z_1D(s))
         i = MIN(MAX(i,1),nr)
         j = MIN(MAX(j,1),nphi)
         k = MIN(MAX(k,1),nz)
         Rminor3D(i,j,k) = Rminor3D(i,j,k) + rminor_1D(s)
         X3D(i,j,k) = X3D(i,j,k) + X_1D(s)
         Y3D(i,j,k) = Y3D(i,j,k) + Y_1D(s)
         N3D(i,j,k) = N3D(i,j,k) + 1
      END DO
      CALL MPI_BARRIER(comm_read,istat)
      IF (mylocalid==master) THEN
         WHERE(N3D > 0)
            Rminor3D = Rminor3D / N3D
            X3D = X3D / N3D
            Y3D = Y3D / N3D
         END WHERE
      END IF
      CALL MPI_BARRIER(comm_read,istat)
      IF (ASSOCIATED(N3D))      CALL mpidealloc(N3D,   win_N3D)

      ! Calculate magaxis
      CALL MPI_CALC_MYRANGE(comm_read,1,nphi,mystart,myend)
      DO j = mystart, myend
         p1 = phiaxis(j) - dp * 0.5
         p2 = p1 + dp
         mask_axis = .FALSE.
         WHERE(zeta_lines(1,:) >= p1 .AND. zeta_lines(1,:) < p2) mask_axis = .TRUE.
         RMAGAXIS(j) = SUM(R_lines(1,:), MASK=mask_axis) / COUNT(mask_axis)
         ZMAGAXIS(j) = SUM(Z_lines(1,:), MASK=mask_axis) / COUNT(mask_axis)
      END DO

      ! Deallocated local variables
      DEALLOCATE(mask_axis)

#if defined(MPI_OPT)
      CALL MPI_BARRIER(comm_read,istat)
#endif
      ! Locals
      IF (ASSOCIATED(zeta_lines)) NULLIFY(zeta_lines)
      IF (ASSOCIATED(zeta_1D))   CALL mpidealloc(zeta_1D,   win_zeta_1D)

      ! Globals
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
   END SUBROUTINE setup_fieldlines_rhogrid

   SUBROUTINE setup_fieldlines_helpers
      IMPLICIT NONE
      INTEGER :: i, j, k
      INTEGER :: ij(2)
      dp = (phiaxis(nphi)-phiaxis(1))/DBLE(nphi-1)
      dr = (raxis(nr)-raxis(1))/DBLE(nr-1)
      dz = (zaxis(nz)-zaxis(1))/DBLE(nz-1)
      RETURN
   END SUBROUTINE setup_fieldlines_helpers

   SUBROUTINE get_fieldlines_grid(nr_out,nz_out,nphi_out,rmin_out,rmax_out,zmin_out,zmax_out,phimax_out)
      IMPLICIT NONE
      INTEGER, INTENT(out) :: nr_out, nz_out, nphi_out
      REAL(rprec), INTENT(out) :: rmin_out, rmax_out, zmin_out, zmax_out, phimax_out
      nr_out = nr; nz_out= nz; nphi_out = nphi
      rmin_out = raxis(1); zmin_out = zaxis(1)
      rmax_out = raxis(nr); zmax_out = zaxis(nz)
      phimax_out = phiaxis(nphi)
      RETURN
   END SUBROUTINE get_fieldlines_grid

   SUBROUTINE get_fieldlines_magaxis(phi_in,raxis_out,zaxis_out)
      IMPLICIT NONE
      REAL(rprec), INTENT(in)  :: phi_in
      REAL(rprec), INTENT(out) :: raxis_out
      REAL(rprec), INTENT(out) :: zaxis_out
      INTEGER :: km, kp
      REAL(rprec) :: phim, dk
      raxis_out = -1; zaxis_out = 0
      phim = MOD(phi_in,phiaxis(nphi))
      IF (phim < 0) phim = phim + pi2
      km = FLOOR((phim)/dp)
      kp = CEILING((phim)/dp)
      km = MIN(MAX(km,1),nphi-1)
      kp = MIN(MAX(kp,2),nphi)
      dk = phim  -km*dp
      raxis_out = RMAGAXIS(km) + (RMAGAXIS(kp)-RMAGAXIS(km))*dk
      zaxis_out = ZMAGAXIS(km) + (ZMAGAXIS(kp)-ZMAGAXIS(km))*dk
      RETURN
   END SUBROUTINE get_fieldlines_magaxis

   SUBROUTINE get_fieldlines_B(r,phi,z,br,bp,bz,rminor,theta)
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
   END SUBROUTINE get_fieldlines_B

   SUBROUTINE get_fieldlines_gridB(i,j,k,br,bp,bz,rminor,theta)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: i,j,k
      REAL(rprec), INTENT(out) :: br, bp, bz
      REAL(rprec), INTENT(out), OPTIONAL :: rminor, theta
      br = 0; bp = 0; bz=0
      IF (i>nr .or. j>nphi .or. k>nz) RETURN
      bp = BPHI3D(i,j,k)
      br = BR3D(i,j,k) * bp / raxis(i)
      bz = BZ3D(i,j,k) * bp / raxis(i)
      IF (PRESENT(rminor)) THEN
         rminor = Rminor3D(i,j,k)
      END IF
      IF (PRESENT(theta)) THEN
         theta = ATAN2(Y3D(i,j,k),X3D(i,j,k))
      END IF
      RETURN
   END SUBROUTINE get_fieldlines_gridB

   SUBROUTINE read_fieldlines_deallocate
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
   END SUBROUTINE read_fieldlines_deallocate

END MODULE read_fieldlines_mod
