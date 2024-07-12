!-----------------------------------------------------------------------
!     Module:        read_hint_mod
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          02/24/2021
!     Description:   This module stores routines for reading
!                    HINT2 output data
!-----------------------------------------------------------------------
      MODULE read_hint_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec

!-----------------------------------------------------------------------
!     Module Variables
!           NTITLE:   
!
!-----------------------------------------------------------------------
      IMPLICIT NONE

      ! HINT Variables
      INTEGER, PARAMETER, PRIVATE :: DTYPE =  SELECTED_REAL_KIND(15)
      INTEGER ::  tstep, nr, nz, nphi, nfp
      DOUBLE PRECISION, PRIVATE, PARAMETER :: one           = 1.0D0 ! 1.0
      DOUBLE PRECISION, PRIVATE :: time, rmin, rmax, zmin, zmax, phimax
      DOUBLE PRECISION, PRIVATE :: pi2, dr, dp, dz
      REAL(DTYPE), DIMENSION(:), POINTER, PRIVATE :: raxis, zaxis
      REAL(DTYPE), DIMENSION(:,:,:), POINTER, PRIVATE :: BR3D, BPHI3D, BZ3D
      REAL(DTYPE), DIMENSION(:,:,:), POINTER, PRIVATE :: VR3D, VPHI3D, VZ3D
      REAL(DTYPE), DIMENSION(:,:,:), POINTER, PRIVATE :: PRES

      CONTAINS

!-----------------------------------------------------------------------
!     Module Subroutines
!           read_hint_mag:  For reading the hint mag files.
!
!-----------------------------------------------------------------------

      SUBROUTINE read_hint_mag(filename,istat)
         IMPLICIT NONE
         CHARACTER(*), INTENT(in) :: filename
         INTEGER, INTENT(out) :: istat
         LOGICAL :: lexist
         INTEGER :: iunit, i, j, k
         REAL(DTYPE) :: time_d, rmin_d, rmax_d, zmin_d, zmax_d
         REAL(DTYPE), DIMENSION(:,:,:,:), ALLOCATABLE :: d_hint
         istat = 0

         pi2 = 8.0 * ATAN(1.0)
         ! Check for the file
         INQUIRE(FILE=TRIM(filename),EXIST=lexist)
         IF (.not.lexist) THEN
            istat=-1
            WRITE(6,*) " ERROR: Could not find file: "//TRIM(filename)
            RETURN
         END IF

         ! Open File
         OPEN(UNIT=iunit, STATUS='old', FILE=TRIM(filename), &
          POSITION='rewind', ACTION='read',&
           FORM='unformatted', IOSTAT=istat, CONVERT='big_endian')

         ! Read Grid Information
         READ(iunit) tstep
         READ(iunit) time_d
         READ(iunit) nr, nz, nphi, nfp
         READ(iunit) rmin_d, zmin_d, rmax_d, zmax_d

         ! Correct
         nphi = nphi + 1

         rmin = rmin_d; rmax = rmax_d; zmin = zmin_d; zmax=zmax_d; phimax=8.0 * ATAN(one) / nfp

         ! Read magnetic field
         ALLOCATE(BR3D(nr,nz,nphi),BPHI3D(nr,nz,nphi),BZ3D(nr,nz,nphi))
         READ(iunit) (((BR3D(i,j,k), BPHI3D(i,j,k), BZ3D(i,j,k), i = 1, nr), j = 1, nz), k = 1, nphi-1)

         ! Read velocity field
         ALLOCATE(VR3D(nr,nz,nphi),VPHI3D(nr,nz,nphi),VZ3D(nr,nz,nphi))
         READ(iunit) (((VR3D(i,j,k), VPHI3D(i,j,k), VZ3D(i,j,k), i = 1, nr), j = 1, nz), k = 1, nphi-1)

         ! Read Pressure field
         ALLOCATE(PRES(nr,nz,nphi))
         READ(iunit) (((PRES(i,j,k), i = 1, nr), j = 1, nz), k = 1, nphi-1)
         CLOSE(iunit)
         PRES = PRES / (pi2 *2E-7)

         ! Replicate points
         BR3D(:,:,nphi)   = BR3D(:,:,1)
         BPHI3D(:,:,nphi) = BPHI3D(:,:,1)
         BZ3D(:,:,nphi)   = BZ3D(:,:,1)
         VR3D(:,:,nphi)   = VR3D(:,:,1)
         VPHI3D(:,:,nphi) = VPHI3D(:,:,1)
         VZ3D(:,:,nphi)   = VZ3D(:,:,1)
         PRES(:,:,nphi)   = PRES(:,:,1)

         ! Helpers
         CALL setup_hint_helpers

         RETURN

      END SUBROUTINE read_hint_mag

      SUBROUTINE setup_hint_helpers
         IMPLICIT NONE
         INTEGER :: i, j, k
         INTEGER :: ij(2)
         dp = phimax/DBLE(nphi-1)
         dr = (rmax-rmin)/DBLE(nr-1)
         dz = (zmax-zmin)/DBLE(nz-1)
         IF (ASSOCIATED(PRES)) THEN
            ALLOCATE(raxis(nphi),zaxis(nphi))
            DO k = 1 , nphi
               ij = MAXLOC(PRES(:,:,k))
               raxis(k) = rmin + dr*ij(1)
               zaxis(k) = zmin + dz*ij(2)
            END DO
         END IF
         RETURN
      END SUBROUTINE setup_hint_helpers

      SUBROUTINE get_hint_grid(nr_out,nz_out,nphi_out,rmin_out,rmax_out,zmin_out,zmax_out,phimax_out)
         IMPLICIT NONE
         INTEGER, INTENT(out) :: nr_out, nz_out, nphi_out
         REAL(rprec), INTENT(out) :: rmin_out, rmax_out, zmin_out, zmax_out, phimax_out
         nr_out = nr; nz_out= nz; nphi_out = nphi
         rmin_out = rmin; zmin_out = zmin
         rmax_out = rmax; zmax_out = zmax
         phimax_out = phimax
         RETURN
      END SUBROUTINE get_hint_grid

      SUBROUTINE get_hint_maxp(presmax)
         IMPLICIT NONE
         REAL(rprec), INTENT(out) :: presmax
         presmax = -1
         IF (ASSOCIATED(PRES)) presmax = MAXVAL(MAXVAL(MAXVAL(PRES,DIM=3),DIM=2),DIM=1)
         RETURN
      END SUBROUTINE get_hint_maxp

      SUBROUTINE get_hint_magaxis(phi_in,raxis_out,zaxis_out)
         IMPLICIT NONE
         REAL(rprec), INTENT(in)  :: phi_in
         REAL(rprec), INTENT(out) :: raxis_out
         REAL(rprec), INTENT(out) :: zaxis_out
         INTEGER :: km, kp
         REAL(rprec) :: phim, dk
         raxis_out = -1; zaxis_out = 0
         phim = MOD(phi_in,phimax)
         IF (phim < 0) phim = phim + pi2
         km = FLOOR((phim)/dp) 
         kp = CEILING((phim)/dp)
         km = MIN(MAX(km,1),nphi-1)
         kp = MIN(MAX(kp,2),nphi)
         dk = phim  -km*dp
         raxis_out = raxis(km) + (raxis(kp)-raxis(km))*dk
         zaxis_out = zaxis(km) + (zaxis(kp)-zaxis(km))*dk
         RETURN
      END SUBROUTINE get_hint_magaxis

      SUBROUTINE get_hint_B(r,phi,z,br,bp,bz)
         IMPLICIT NONE
         REAL(rprec), INTENT(in) :: r, phi, z
         REAL(rprec), INTENT(out) :: br, bp, bz
         INTEGER :: im, ip, km, kp, jp, jm
         REAL(rprec) :: di, dk, dj, phim, f1, f2
         br=0; bp=0; bz=0
         ! Check for bounds
         IF ((r < rmin) .or. (r>rmax) .or. &
             (z < zmin) .or. (z> zmax)) RETURN
         phim = MOD(phi,phimax)
         IF (phim < 0) phim = phim + pi2
         ! Grid helpers
         im = FLOOR((r-rmin)/dr) 
         ip = CEILING((r-rmin)/dr) 
         jm = FLOOR((z-zmin)/dz) 
         jp = CEILING((z-zmin)/dz)
         km = FLOOR((phim)/dp) 
         kp = CEILING((phim)/dp)
         im = MIN(MAX(im,1),nr-1); jm = MIN(MAX(jm,1),nz-1) ; km = MIN(MAX(km,1),nphi-1)
         ip = MIN(MAX(ip,2),nr);   jp = MIN(MAX(jp,2),nz) ;   kp = MIN(MAX(kp,2),nphi)

         di = r-rmin-im*dr
         dj = z-zmin-jm*dz
         dk = phim  -km*dp

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
         RETURN
      END SUBROUTINE get_hint_B

      SUBROUTINE get_hint_gridB(i,j,k,br,bp,bz)
         IMPLICIT NONE
         INTEGER, INTENT(in) :: i,j,k
         REAL(rprec), INTENT(out) :: br, bp, bz
         br = 0; bp = 0; bz=0
         IF (i>nr .or. j>nphi .or. k>nz) RETURN
         br = BR3D(i,k,j)
         bp = BPHI3D(i,k,j)
         bz = BZ3D(i,k,j)
         RETURN
      END SUBROUTINE get_hint_gridB

      SUBROUTINE get_hint_press(r,phi,z,press)
         IMPLICIT NONE
         REAL(rprec), INTENT(in) :: r, phi, z
         REAL(rprec), INTENT(out) :: press
         INTEGER :: im, ip, km, kp, jp, jm
         REAL(rprec) :: di, dk, dj, f1, f2, phim
         press=0
         ! Check for bounds
         IF ((r < rmin) .or. (r>rmax) .or. &
             (z < zmin) .or. (z> zmax)) RETURN
         phim = MOD(phi,phimax)
         IF (phim < 0) phim = phim + pi2
         ! Grid helpers
         im = FLOOR((r-rmin)/dr) 
         ip = CEILING((r-rmin)/dr) 
         jm = FLOOR((z-zmin)/dz) 
         jp = CEILING((z-zmin)/dz)
         km = FLOOR((phim)/dp) 
         kp = CEILING((phim)/dp)
         im = MIN(MAX(im,1),nr-1); jm = MIN(MAX(jm,1),nz-1) ; km = MIN(MAX(km,1),nphi-1)
         ip = MIN(MAX(ip,2),nr);   jp = MIN(MAX(jp,2),nz) ;   kp = MIN(MAX(kp,2),nphi)

         di = r-rmin-im*dr
         dj = z-zmin-jm*dz
         dk = phim  -km*dp

         f1 = PRES(im,jm,km)*(one-di)*(one-dj) &
            + PRES(ip,jm,km)*(di)*(one-dj) &
            + PRES(im,jp,km)*(one-di)*(dj) &
            + PRES(ip,jp,km)*di*dj
         f2 = PRES(im,jm,kp)*(one-di)*(one-dj) &
            + PRES(ip,jm,kp)*(di)*(one-dj) &
            + PRES(im,jp,kp)*(one-di)*(dj) &
            + PRES(ip,jp,kp)*di*dj
         press = f1+ (f2-f1)*dk
         RETURN
      END SUBROUTINE get_hint_press

      SUBROUTINE get_hint_gridpress(i,j,k,press)
         IMPLICIT NONE
         INTEGER, INTENT(in) :: i,j,k
         REAL(rprec), INTENT(out) :: press
         press = 0
         IF (i>nr .or. j>nphi .or. k>nz) RETURN
         press = PRES(i,k,j)
         RETURN
      END SUBROUTINE get_hint_gridpress

      SUBROUTINE read_hint_deallocate
         IMPLICIT NONE
         IF (ASSOCIATED(BR3D))   DEALLOCATE(BR3D)
         IF (ASSOCIATED(BPHI3D)) DEALLOCATE(BPHI3D)
         IF (ASSOCIATED(BZ3D))   DEALLOCATE(BZ3D)
         IF (ASSOCIATED(VR3D))   DEALLOCATE(VR3D)
         IF (ASSOCIATED(VPHI3D)) DEALLOCATE(VPHI3D)
         IF (ASSOCIATED(VZ3D))   DEALLOCATE(VZ3D)
         IF (ASSOCIATED(PRES))   DEALLOCATE(PRES)
         IF (ASSOCIATED(raxis))  DEALLOCATE(raxis)
         IF (ASSOCIATED(zaxis))  DEALLOCATE(zaxis)
         RETURN
      END SUBROUTINE read_hint_deallocate

      END MODULE read_hint_mod