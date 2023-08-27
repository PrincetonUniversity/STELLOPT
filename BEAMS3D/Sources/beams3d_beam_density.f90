!-----------------------------------------------------------------------
!     Module:        beams3d_beam_density
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/03/2023
!     Description:   This subroutine computes the neutral beam density
!                    on a grid.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_beam_density
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE beams3d_runtime
      USE beams3d_grid
      USE beams3d_lines
      USE mpi_params ! MPI
      USE mpi_inc
      USE mpi_sharmem
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error Flag
!          iunit        File ID
!          ndist        Number of Vll divisions for dist function
!          ns           Number of flux divisions for current calculation
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL :: lline_in_box
      INTEGER :: ier, l, nl, i, j, k, m, s
      DOUBLE PRECISION :: x0, y0, z0, x1, y1, z1, hr2, hz2, hp2, d, &
                          denbeam, dl, xt, yt, zt, rt, pt, dV, &
                          dx, dy, dz , dgrid
      LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE :: lgrid
#if defined(MPI_OPT)
      INTEGER :: MPI_COMM_LOCAL
      INTEGER :: mystart, mypace
#endif
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      ! Allocate the Grid
      CALL mpialloc(BEAM_DENSITY, nbeams, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_BEAM_DENSITY)
      ALLOCATE(lgrid(nbeams,nr,nphi,nz))

      IF (lverb) WRITE(6,'(A)') '----- BEAM Density Calc. -----'

      dgrid = sqrt(hr(1)**2+hz(1)**2)
      hr2 = 0.5*hr(1)
      hz2 = 0.5*hz(1)
      hp2 = 0.5*hp(1)
      ! Loop over particles
      DO l = mystart_save, myend_save
         x0 = R_lines(0,l)*cos(PHI_lines(0,l))
         y0 = R_lines(0,l)*sin(PHI_lines(0,l))
         z0 = Z_lines(0,l)
         x1 = R_lines(1,l)*cos(PHI_lines(1,l))
         y1 = R_lines(1,l)*sin(PHI_lines(1,l))
         z1 = Z_lines(1,l)
         m  = Beam(l)
         ! Because we're in cylindrical space we just need to brute force this
         ! We have a chord with particles/s (weight) going along it.
         ! Since the velocity is constant then the TOF is just L/v = t
         ! The the total number of particles is just n=W*t = W*L/v
         ! then we jsut need to calc the volume of the voxel to get particles/m^3
         d  = SQRT((x1-x0)**2 + (y1-y0)**2 + (z1-z0)**2)
         denbeam = weight(l)*L/vll_lines(0,l) ! This is the total number of particles
         dl = d/dgrid
         nl = NINT(L/dl)+1
         dl = d/nl
         dx = (x1-x0)/(nl-1)
         dy = (y1-y0)/(nl-1)
         dz = (z1-z0)/(nl-1)
         LGRID = .FALSE.
         DO s = 1, nl
            xt = x0 + (i-1)*dx
            yt = y0 + (i-1)*dy
            zt = z0 + (i-1)*dz
            rt = SQRT(xt*xt+yt*yt)
            IF ( (rt > raxis(nr-1)) .or. (rt < raxis(2)) .or. &
                 (zt > zaxis(nz-1)) .or. (zt < zaxis(2)) ) CYCLE
            pt = ATAN2(yt,xt)
            pt = MODULO(pt, phimax)
            ! Find gridpoint but use half grid
            i = COUNT(raxis-hr2 < rt)
            j = COUNT(phiaxis-hp2 < pt)
            k = COUNT(zaxis-hz2 < zt)
            LGRID(m,i,j,k) = .true.
         END DO
         WHERE(lgrid) BEAM_DENSITY = BEAM_DENSITY + denbeam
      END DO
      DEALLOCATE(lgrid)

      ! Barrier so  calculation is done
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
#endif

      ! Now Divide by the grid volumes
      CALL MPI_CALC_MYRANGE(MPI_COMM_SHARMEM, 1, (nr-1)*(nphi-1)*(nz-1), mystart, myend)
      DO s = mystart,myend
         i = MOD(s-1,nr-1)+1
         j = MOD(s-1,(nr-1)*(nphi-1))
         j = FLOOR(REAL(j) / REAL(nr-1))+1
         k = CEILING(REAL(s) / REAL((nr-1)*(nphi-1)))
         dV = raxis(i)*hr(i)*hp(j)*hz(k)
         BEAM_DENSITY(:,i,j,k) = BEAM_DENSITY(:,i,j,k) / dV
      END DO

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
#endif

      RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!----------------------------------------------------------------------- 
      END SUBROUTINE beams3d_beam_density
