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
                          denbeam, dl, dV
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xl,yl,zl,rl,sl,pl
#if defined(MPI_OPT)
      INTEGER :: MPI_COMM_LOCAL
      INTEGER :: mystart, mypace
#endif
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      ! Allocate the Grid
      CALL mpialloc(BEAM_DENSITY, nbeams, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_BEAM_DENSITY)
      IF (myid_sharmem == 0) BEAM_DENSITY = 0
      nl = nr + nr/2 + 1
      ALLOCATE(xl(nl),yl(nl),zl(nl),rl(nl),pl(nl),sl(nl))
      FORALL(s=1:nl) sl(s) = DBLE(s-0.5)/DBLE(nl) ! half grid


      IF (lverb) WRITE(6,'(A)') '----- BEAM Density Calc. -----'
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
         ! then we just need to calc the volume of the voxel to get particles/m^3
         d  = SQRT((x1-x0)**2 + (y1-y0)**2 + (z1-z0)**2)/DBLE(nl) ! so it's dl
         denbeam = weight(l)*d/vll_lines(0,l) ! This is the total number of particles per step
         xl = x0 + sl*(x1-x0)
         yl = y0 + sl*(y1-y0)
         zl = z0 + sl*(z1-z0)
         rl = sqrt(xl*xl+yl*yl)
         pl = atan2(yl,xl)
         pl = MODULO(pl,phimax)
         DO s = 1, nl
            ! Find gridpoint but use half grid
            i = MIN(MAX(COUNT(raxis-hr2 < rl(s)),1),nr)
            j = MIN(MAX(COUNT(phiaxis-hp2 < pl(s)),1),nphi)
            k = MIN(MAX(COUNT(zaxis-hz2 < zl(s)),1),nz)
            BEAM_DENSITY(m,i,j,k) = BEAM_DENSITY(m,i,j,k) + denbeam
         END DO
      END DO
      DEALLOCATE(xl,yl,zl,rl,pl,sl)

      ! Fix edges which are double counts
      IF (myid_sharmem == 0) THEN
         BEAM_DENSITY(:,1,:,:) = 0
         BEAM_DENSITY(:,nr,:,:) = 0
         BEAM_DENSITY(:,:,:,1) = 0
         BEAM_DENSITY(:,1,:,nz) = 0
      ENDIF


      ! Barrier so  calculation is done
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_SHARMEM,ierr_mpi)
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
      CALL MPI_BARRIER(MPI_COMM_SHARMEM,ierr_mpi)
#endif

      RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!----------------------------------------------------------------------- 
      END SUBROUTINE beams3d_beam_density
