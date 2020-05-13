!-----------------------------------------------------------------------
!     Module:        beams3d_randomize_particles
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          05/07/2020
!     Description:   This subroutine radomizes the particles so that
!                    no processor follows only slow or fast particles.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_randomize_particles
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_lines, ONLY: nparticles
      USE mpi_params
      USE mpi_inc
      USE beams3d_runtime

!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      INTEGER :: i, d1, d2, b_temp
      INTEGER, DIMENSION(:), ALLOCATABLE :: randomindex
      REAL(rprec) :: r_temp, phi_temp, z_temp, vll_temp, w_temp, &
                     za_temp, m_temp, c_temp, mu_temp, t_temp
      REAL(rprec), DIMENSION(3) :: v_temp
      REAL, DIMENSION(:), ALLOCATABLE :: randomnumbers

!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      IF (myworkid == master) THEN
         ALLOCATE(randomindex(nparticles))
         ALLOCATE(randomnumbers(nparticles))
         CALL RANDOM_NUMBER(randomnumbers)
         randomindex = NINT(randomnumbers*nparticles)
         WHERE(randomindex < 1) randomindex = 1
         DEALLOCATE(randomnumbers)
         DO i = 1, nparticles
            d1 = i
            d2 = randomindex(i)
            ! Save swap
            r_temp = R_start(d2)
            phi_temp = phi_start(d2)
            z_temp   = Z_start(d2)
            vll_temp = vll_start(d2)
            w_temp   = weight(d2)
            za_temp  = Zatom(d2)
            m_temp   = mass(d2)
            c_temp   = charge(d2)
            mu_temp  = mu_start(d2)
            t_temp   = t_end(d2)
            v_temp   = v_neut(:,d2)
            b_temp   = Beam(d2)
            !swap
            R_start(d2)   = R_start(d1)
            phi_start(d2) = phi_start(d1)
            Z_start(d2)   = Z_start(d1)
            vll_start(d2) = vll_start(d1)
            weight(d2)    = weight(d1)
            Zatom(d2)     = Zatom(d1)
            mass(d2)      = mass(d1)
            charge(d2)    = charge(d1)
            mu_start(d2)  = mu_start(d1)
            t_end(d2)     = t_end(d1)
            v_neut(:,d2)  = v_neut(:,d1)
            Beam(d2)      = Beam(d1)
            ! Finish
            R_start(d1)   = r_temp
            phi_start(d1) = phi_temp
            Z_start(d1)   = z_temp
            vll_start(d1) = vll_temp
            weight(d1)    = w_temp
            Zatom(d1)     = za_temp
            mass(d1)      = m_temp
            charge(d1)    = c_temp
            mu_start(d1)  = mu_temp
            t_end(d1)     = t_temp
            v_neut(:,d1)  = v_temp
            Beam(d1)      = b_temp
         END DO
         DEALLOCATE(randomindex)
      END IF

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(R_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(phi_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(Z_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(vll_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(weight,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(Zatom,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(mass,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(charge,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(mu_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(t_end,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(v_neut,nparticles*3,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(beam,nparticles,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
#endif

      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_randomize_particles