!-----------------------------------------------------------------------
!     Module:        beams3d_init_restart
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/27/2012
!     Description:   This subroutine loads data from a previous run.
!                    For now we define local variables to control how
!                    this is done.  In the future we will let the user
!                    do this interactively.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_init_restart
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_runtime
      USE beams3d_grid
      USE beams3d_lines
#if defined(LHDF5)
      USE ez_hdf5
#endif
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          npoinc_extract Which save state to extract from file.
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL :: lplasma_old, ldepo_old
      INTEGER :: i, k, ier, npoinc_extract, npoinc_save, state_flag
      INTEGER, DIMENSION(:), ALLOCATABLE :: beam2
      REAL(rprec) :: vpartmax
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: mass2, charge2, Zatom2, &
                                                weight2, t_end2
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      IF (lverb) THEN
         WRITE(6,'(A)')  '----- Reading Restart File -----'
         WRITE(6,'(A)')  '   FILE: '//TRIM(restart_string)
      END IF

      IF (myworkid == master) THEN
         ! Save quantities
         npoinc_save = npoinc
         ! Read the data
         CALL open_hdf5(TRIM(restart_string),fid,ier,LCREATE=.false.)
         IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,TRIM(restart_string),ier)
         CALL read_scalar_hdf5(fid,'nparticles',ier,INTVAR=nparticles)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nparticles',ier)
         CALL read_scalar_hdf5(fid,'npoinc',ier,INTVAR=npoinc)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'npoinc',ier)
         !CALL read_scalar_hdf5(fid,'partvmax',ier,DBLVAR=vpartmax)
         !partvmax = MAX(partvmax,vpartmax)
         !IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'partvmax',ier)
         IF (ALLOCATED(v_neut)) DEALLOCATE(v_neut)
         IF (ALLOCATED(mass)) DEALLOCATE(mass)
         IF (ALLOCATED(charge)) DEALLOCATE(charge)
         IF (ALLOCATED(Zatom)) DEALLOCATE(charge)
         IF (ALLOCATED(beam)) DEALLOCATE(beam)
         IF (ALLOCATED(weight)) DEALLOCATE(weight)
         IF (ALLOCATED(end_state)) DEALLOCATE(end_state)
         IF (ALLOCATED(R_lines)) DEALLOCATE(R_lines)
         IF (ALLOCATED(PHI_lines)) DEALLOCATE(PHI_lines)
         IF (ALLOCATED(Z_lines)) DEALLOCATE(Z_lines)
         IF (ALLOCATED(moment_lines)) DEALLOCATE(moment_lines)
         IF (ALLOCATED(vll_lines)) DEALLOCATE(vll_lines)
         IF (ALLOCATED(neut_lines)) DEALLOCATE(neut_lines)
         ALLOCATE(mass2(nparticles),charge2(nparticles),Zatom2(nparticles),&
            beam2(nparticles), weight2(nparticles), t_end2(nparticles), &
            end_state(nparticles))
         ALLOCATE(R_lines(0:npoinc,nparticles),Z_lines(0:npoinc,nparticles),PHI_lines(0:npoinc,nparticles),&
            vll_lines(0:npoinc,nparticles),neut_lines(0:npoinc,nparticles),moment_lines(0:npoinc,nparticles),&
            S_lines(0:npoinc,nparticles))
         CALL read_var_hdf5(fid,'mass',nparticles,ier,DBLVAR=mass2)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'mass2',ier)
         CALL read_var_hdf5(fid,'charge',nparticles,ier,DBLVAR=charge2)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'charge2',ier)
         CALL read_var_hdf5(fid,'Zatom',nparticles,ier,DBLVAR=Zatom2)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'Zatom2',ier)
         CALL read_var_hdf5(fid,'Weight',nparticles,ier,DBLVAR=weight2)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'weight2',ier)
         CALL read_var_hdf5(fid,'Beam',nparticles,ier,INTVAR=beam2)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'beam2',ier)
         CALL read_var_hdf5(fid,'end_state',nparticles,ier,INTVAR=end_state)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'end_state',ier)
         CALL read_var_hdf5(fid,'t_end',nparticles,ier,DBLVAR=t_end2)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'t_end2',ier)
         CALL read_var_hdf5(fid,'R_lines',npoinc+1,nparticles,ier,DBLVAR=R_lines)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'R_lines',ier)
         CALL read_var_hdf5(fid,'Z_lines',npoinc+1,nparticles,ier,DBLVAR=Z_lines)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'Z_lines',ier)
         CALL read_var_hdf5(fid,'PHI_lines',npoinc+1,nparticles,ier,DBLVAR=PHI_lines)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'PHI_lines',ier)
         CALL read_var_hdf5(fid,'vll_lines',npoinc+1,nparticles,ier,DBLVAR=vll_lines)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'vll_lines',ier)
         CALL read_var_hdf5(fid,'neut_lines',npoinc+1,nparticles,ier,BOOVAR=neut_lines)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'neut_lines',ier)
         CALL read_var_hdf5(fid,'moment_lines',npoinc+1,nparticles,ier,DBLVAR=moment_lines)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'moment_lines',ier)
         CALL read_var_hdf5(fid,'S_lines',npoinc+1,nparticles,ier,DBLVAR=S_lines)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'S_lines',ier)
         CALL close_hdf5(fid,ier)
         IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'beams3d_'//TRIM(restart_string)//'.h5',ier)

         ! Decide what to do
         !   IF depo run then start from initial born particle population
         !   Otherwise start from wall hits
         ldepo_old = .false.
         state_flag = 0
         IF (ANY(end_state==3)) ldepo_old = .true.
         IF (ldepo_old) THEN
            state_flag = 0
            IF (lplasma_only) THEN 
               WHERE(S_lines(1,:) >= 1) end_state = -1
            END IF
         ELSE
            state_flag = 2
         END IF
         k = COUNT(end_state == state_flag)

         ! Allocate the particles
         ALLOCATE(  R_start(k), phi_start(k), Z_start(k), &
                    v_neut(3,k), mass(k), charge(k), &
                    mu_start(k), Zatom(k), t_end(k), vll_start(k), &
                    beam(k), weight(k) )

         ! Now fill the arrays downselecting for non-shinethrough particles
         k = 1
         DO i = 1, nparticles
            IF (end_state(i) /= state_flag) CYCLE
            npoinc_extract = COUNT(R_lines(:,i)>0)
            R_start(k)   = R_lines(npoinc_extract,i)
            Z_start(k)   = Z_lines(npoinc_extract,i)
            phi_start(k) = PHI_lines(npoinc_extract,i)
            vll_start(k) = vll_lines(npoinc_extract,i)
            mu_start(k)  = moment_lines(npoinc_extract,i)
            v_neut(3,k)   = 0.0
            mass(k)      = mass2(i)
            charge(k)   = charge2(i)
            Zatom(k)    = Zatom2(i)
            beam(k)     = beam2(i)
            weight(k)   = weight2(i)
            t_end(k)    = MAXVAL(t_end_in)
            k = k + 1
         END DO
         DEALLOCATE(R_lines, Z_lines, PHI_lines, vll_lines, moment_lines, neut_lines, end_state, S_lines)
         DEALLOCATE(mass2, charge2, Zatom2, beam2, weight2, t_end2)

         ! Restore quantities
         nparticles = k-1
         npoinc = npoinc_save
         nbeams = MAXVAL(beam)
         IF (lverb) THEN
            WRITE(6,'(A,I8)') '   # of Particles: ', nparticles
            WRITE(6,'(A,I6)') '   # of Beams: ', nbeams
         END IF
      END IF

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(nparticles,1,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(nbeams,1,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(partvmax,1,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      IF (myworkid /= master) THEN
         ALLOCATE(  R_start(nparticles), phi_start(nparticles), Z_start(nparticles), &
                    v_neut(3,nparticles), mass(nparticles), charge(nparticles), &
                    mu_start(nparticles), Zatom(nparticles), t_end(nparticles), vll_start(nparticles), &
                    beam(nparticles), weight(nparticles) )
      END IF
      CALL MPI_BCAST(mu_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(t_end,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(mass,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(charge,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(Zatom,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(weight,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(R_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(phi_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(Z_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(vll_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(beam,nparticles,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(v_neut,nparticles*3,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
#endif


!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_init_restart
