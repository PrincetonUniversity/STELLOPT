!-----------------------------------------------------------------------
!     Module:        beams3d_init_beams_bbnbi
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          12/05/2019
!     Description:   This subroutine initializes the beam using the
!                    BBNBI formulation where each beam is modeled using
!                    beamlets.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_init_beams_bbnbi
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_runtime
      USE beams3d_lines, ONLY: nparticles
      USE mpi_params
      USE mpi_inc
      USE hdf5

!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, i, j, k, k1, k2
      INTEGER, DIMENSION(:), ALLOCATABLE :: N_start
      REAL(rprec) :: rtemp
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: Energy, X_start, Y_start
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: X, Y, X_BEAMLET, Y_BEAMLET, Z_BEAMLET, &
                                                  NX_BEAMLET, NY_BEAMLET, NZ_BEAMLET, U, V
      REAL(rprec), PARAMETER   :: E_error = .01 ! 1% energy spread

      ! For HDF5
      INTEGER(HID_T)           :: h5_fid, h5_did, h5_sid
      INTEGER(HSIZE_T), DIMENSION(2)    :: dims, maxdims

!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      CALL init_random_seed

      IF (myworkid == master) THEN
         CALL H5open_f(ier)
         CALL h5fopen_f(TRIM(bbnbi_string), H5F_ACC_RDONLY_F, h5_fid, ier)
         CALL h5dopen_f (h5_fid, 'X_BEAM', h5_did, ier)
         CALL h5dget_space_f(h5_did, h5_sid ,ier)
         CALL h5sget_simple_extent_dims_f(h5_sid, dims, maxdims,ier)
         k1 = dims(1)
         k2 = dims(2)
         ALLOCATE(X_BEAMLET(k1,k2), Y_BEAMLET(k1,k2), Z_BEAMLET(k1,k2),&
                  NX_BEAMLET(k1,k2), NY_BEAMLET(k1,k2), NZ_BEAMLET(k1,k2))
         CALL h5dread_f(h5_did, H5T_NATIVE_DOUBLE, X_BEAMLET, dims, ier)
         CALL h5sclose_f(h5_sid, ier)
         CALL h5dclose_f(h5_did, ier)
         CALL h5dopen_f (h5_fid, 'Y_BEAM', h5_did, ier)
         CALL h5dread_f(h5_did, H5T_NATIVE_DOUBLE, Y_BEAMLET, dims, ier)
         CALL h5dclose_f(h5_did, ier)
         CALL h5dopen_f (h5_fid, 'Z_BEAM', h5_did, ier)
         CALL h5dread_f(h5_did, H5T_NATIVE_DOUBLE, Z_BEAMLET, dims, ier)
         CALL h5dclose_f(h5_did, ier)
         CALL h5dopen_f (h5_fid, 'NX_BEAM', h5_did, ier)
         CALL h5dread_f(h5_did, H5T_NATIVE_DOUBLE, NX_BEAMLET, dims, ier)
         CALL h5dclose_f(h5_did, ier)
         CALL h5dopen_f (h5_fid, 'NY_BEAM', h5_did, ier)
         CALL h5dread_f(h5_did, H5T_NATIVE_DOUBLE, NY_BEAMLET, dims, ier)
         CALL h5dclose_f(h5_did, ier)
         CALL h5dopen_f (h5_fid, 'NZ_BEAM', h5_did, ier)
         CALL h5dread_f(h5_did, H5T_NATIVE_DOUBLE, NZ_BEAMLET, dims, ier)
         CALL h5dclose_f(h5_did, ier)
         CALL h5close_f(ier)

         i = 1
         k2 = 0
         DO
            IF (Dex_beams(i)<1) EXIT
            IF (Dex_beams(i)/=j) THEN
               j = Dex_beams(i)
               k2 = k2+1
            END IF
            i=i+1
         END DO

         IF (lverb) THEN
            WRITE(6, '(A)') '----- INITIALIZING BEAMLET BASED BEAMS -----'
            WRITE(6, '(A,A)')       '   filename: ', TRIM(bbnbi_string)
            WRITE(6, '(A,I4,A,I4)') '   nbeams: ', k2,'/',dims(1)
            WRITE(6, '(A,I4)')      '   nbeamlets: ', dims(2)
            WRITE(6, '(A,I6)')      '   nparticles_start: ', nparticles_start
            CALL FLUSH(6)
         END IF
      END IF

      ! Broadcast and allocate the global variables
      k1 = dims(1)
      k2 = dims(2)
      CALL MPI_BCAST(k1, 1, MPI_INTEGER, master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_BCAST(k2, 1, MPI_INTEGER, master, MPI_COMM_BEAMS, ierr_mpi)
      nparticles = nbeams*nparticles_start
      ALLOCATE(   R_start(nparticles), phi_start(nparticles), Z_start(nparticles), vll_start(nparticles), &
                  v_neut(3,nparticles), mass(nparticles), charge(nparticles), Zatom(nparticles), &
                  mu_start(nparticles), t_end(nparticles), &
                  beam(nparticles), weight(nparticles))
      IF (myworkid==master) THEN
         ALLOCATE(N_start(nparticles_start),X_start(nparticles_start),Y_start(nparticles_start),&
                  Energy(nparticles_start), U(3,nparticles_start), V(3,nparticles_start))
         ! Randomly inidialize beamlets
         CALL RANDOM_NUMBER(X_start)
         N_start = NINT(X_Start*dims(2))
         WHERE(N_start==0) N_start=1
         k1 = 1; k2 = nparticles_start
         weight = 0
         DO i=1,nbeams
            ! Cycle if not using beam
            j = Dex_beams(i)
            IF (lverb) WRITE(6, '(A,I2,A,I4,A,A,I2,A,F7.3,A,A,I2,A,I4)') '            E_BEAM(',i ,'): ',&
                     NINT(E_beams(i)*6.24150636309E15),' [keV]',& !(1.0E-3/ec)
                                     ' P_BEAM(',i ,'): ',(P_beams(i)*1E-6),' [MW]',&
                                     ' DEX_BEAM(',i ,'): ',(Dex_beams(i))

            CALL FLUSH(6)
            ! Basics
            beam(k1:k2)         = i
            mu_start(k1:k2)     = 0
            t_end(k1:k2)        = t_end_in(i)
            mass(k1:k2)         = mass_beams(i)
            charge(k1:k2)       = charge_beams(i)
            Zatom(k1:k2)        = Zatom_beams(i)
            ! Energy distribution
            CALL gauss_rand(nparticles_start, Energy)
            Energy = sqrt( (E_beams(i) + E_error*E_beams(i)*Energy)*(E_beams(i) + E_error*E_beams(i)*Energy) )
            IF (lbeam_simple) Energy = E_beams(i)
            weight(k1:k2)       = P_beams(i)/Energy
            ! Starting Points
            X_start          = X_BEAMLET(j,N_start)
            Y_start          = Y_BEAMLET(j,N_start)
            R_start(k1:k2)   = SQRT(X_start*X_start+Y_start*Y_start)
            PHI_start(k1:k2) = ATAN2(Y_start,X_start)
            Z_start(k1:k2)   = Z_BEAMLET(j,N_start)
            ! Now calculate Divergence (use small angle tan(div)=div here)
            CALL gauss_rand(nparticles_start,X_start)
            CALL gauss_rand(nparticles_start,Y_start)
            U(1,:) = SQRT(X_start*X_start+Y_start*Y_start)
            V(1,:) = ATAN2(Y_start,X_start)
            X_start = Div_beams(i)*U(1,:)*COS(V(1,:))
            Y_start = Div_beams(i)*U(1,:)*SIN(V(1,:))
            ! Calcualte Divergence vectors (assume N are unit vectors)
            U(1,:) = NY_BEAMLET(j,N_start)
            U(2,:) = -NX_BEAMLET(j,N_start)
            U(3,:) = 0.0
            V(1,:) = U(2,:)*NZ_BEAMLET(j,N_start)
            V(2,:) =-U(1,:)*NZ_BEAMLET(j,N_start)
            V(3,:) = U(1,:)*NY_BEAMLET(j,N_start)-U(2,:)*NX_BEAMLET(j,N_start)
            ! Starting Velocity 
            vll_start(k1:k2) = SQRT(2*Energy/mass_beams(i))  ! speed E=0.5*mv^2
            v_neut(1,k1:k2)  = (NX_BEAMLET(j,N_start) + U(1,:)*X_Start + V(1,:)*Y_start)*vll_start(k1:k2)
            v_neut(2,k1:k2)  = (NY_BEAMLET(j,N_start) + U(2,:)*X_Start + V(2,:)*Y_start)*vll_start(k1:k2)
            v_neut(3,k1:k2)  = (NZ_BEAMLET(j,N_start) + U(3,:)*X_Start + V(3,:)*Y_start)*vll_start(k1:k2)
            k1 = k2 + 1
            k2 = k2 + nparticles_start
         END DO
         DEALLOCATE(N_start,X_Start,Y_start,Energy, U, V)
         DEALLOCATE(X_BEAMLET,Y_BEAMLET,Z_BEAMLET,NX_BEAMLET,NY_BEAMLET,NZ_BEAMLET)
         weight = weight/nparticles
      END IF
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
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
!DEC$ ENDIF


!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_init_beams_bbnbi
