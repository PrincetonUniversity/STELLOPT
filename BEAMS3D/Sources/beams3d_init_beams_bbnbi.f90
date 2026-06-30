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
      USE beams3d_lines, ONLY: nparticles, partvmax, is_active, mystart_save, myend_save
      USE beams3d_grid, ONLY: X_BEAMLET, Y_BEAMLET, Z_BEAMLET, &
                           NX_BEAMLET, NY_BEAMLET, NZ_BEAMLET
      USE mpi_params
      USE mpi_inc
      USE hdf5
      USE boxsim_db

!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, i, ibeam, iproc, j, MPI_COMM_LOCAL
      INTEGER :: npart_proc, npart_beam_left, n_from_beam
      INTEGER :: k1, k2, k_rank_start, k_rank_end, n_filled 
      INTEGER, DIMENSION(:), ALLOCATABLE :: N_start
      INTEGER, DIMENSION(:), ALLOCATABLE :: mystart_all, myend_all
      REAL(rprec) :: rtemp
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: Energy, X_start, Y_start
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: X, Y, U, V, v_neut
      REAL(rprec), PARAMETER   :: E_error = .01 ! 1% energy spread
      CHARACTER(LEN=8) :: species
      INTEGER :: part_counts(boxsim_nkinds), charge_int, Z_int, ierr
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

         i = 1; j=0
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
            WRITE(6, '(A,I8)')      '   nparticles_start: ', nparticles_start
            CALL FLUSH(6)
         END IF
         IF (lascot) THEN
            CALL beams3d_write_ascoth5('BBNBI')
            IF (lverb) WRITE(6, '(A,I4)')      '   ASCOT5 File: Updated'
         END IF
      END IF

      ! Broadcast and allocate the global variables
      k1 = dims(1)
      k2 = dims(2)
      CALL MPI_BCAST(k1, 1, MPI_INTEGER, master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_BCAST(k2, 1, MPI_INTEGER, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (myid_sharmem==master) lgc2fo_start = .TRUE.
      ! Tell the master what should be the start and end points
      ALLOCATE(mystart_all(nprocs_beams))
      mystart_all = 0
      mystart_all(myworkid+1) = mystart_save
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, mystart_all, nprocs_beams, MPI_INTEGER, MPI_SUM, MPI_COMM_BEAMS, ierr_mpi)
      IF (myworkid==master) THEN
         npart_proc =  INT(CEILING(DBLE(nbeams*nparticles_start)/DBLE(nprocs_beams))) 
         ALLOCATE(N_start(npart_proc),X_start(npart_proc),Y_start(npart_proc),&
                  Energy(npart_proc), U(3,npart_proc), V(3,npart_proc), &
                  v_neut(3,npart_proc))
         ! Randomly inidialize beamlets
         CALL RANDOM_NUMBER(X_start)
         N_start = NINT(X_Start*dims(2))
         WHERE(N_start==0) N_start=1
         weight = 0

         DO ibeam = 1, nbeams
            IF (lverb) WRITE(6, '(A,I2,A,I4,A,A,I2,A,F7.3,A,A,I2,A,I4)') '            E_BEAM(',ibeam ,'): ',&
            NINT(E_beams(ibeam)*6.24150636309E15),' [keV]',& !(1.0E-3/ec)
                           ' P_BEAM(',ibeam ,'): ',(P_beams(ibeam)*1E-6),' [MW]',&
                           ' DEX_BEAM(',ibeam ,'): ',(Dex_beams(ibeam))
            CALL FLUSH(6)
         END DO
         ! Idea: take particles from beams until empty
         npart_beam_left = nparticles_start
         ibeam = 1 ! beam index
         DO iproc = 1, nprocs_beams
            k_rank_start = mystart_all(iproc)
            k_rank_end = k_rank_start + npart_proc - 1

            n_filled = 0
            DO WHILE (n_filled .LT. npart_proc)
               ! Particles that can be taken
               n_from_beam = MIN(npart_beam_left, npart_proc - n_filled)

               k1 = k_rank_start + n_filled ! if nfilled=0, start from obvious start
               k2 = k_rank_start + n_filled + n_from_beam - 1 ! 
               
               ASSOCIATE( E => Energy(1:n_from_beam), &
                          Nstart => N_start(1:n_from_beam), &
                          Xstart => X_start(1:n_from_beam), &
                          Ystart => Y_start(1:n_from_beam), &
                          U1 => U(1, 1:n_from_beam), U2 => U(2, 1:n_from_beam), U3 => U(3, 1:n_from_beam), &
                          V1 => V(1, 1:n_from_beam), V2 => V(2, 1:n_from_beam), V3 => V(3, 1:n_from_beam), &
                          vneut1 => v_neut(1, 1:n_from_beam), &
                          vneut2 => v_neut(2, 1:n_from_beam), &
                          vneut3 => v_neut(3, 1:n_from_beam))
                          
               ! Basics
               beam(k1:k2)         = ibeam
               mu_start(k1:k2)     = 0
               t_end(k1:k2)        = t_end_in(ibeam)

               IF (lboxsim)  THEN
                  is_active(k1:k2) = .TRUE.    ! set active
                  ! Overwrite data with species if available
                  species = species_beams(ibeam)
                  IF (species/='') THEN
                     boxsim_species(k1:k2) = species_beams(ibeam) ! set species
                     CALL boxsim_parse_species(species, part_counts, charge_int, Z_int, ierr)
                     IF (ierr/=0) THEN
                        WRITE(6,*) 'ERROR: could not parse species string for beam ', ibeam, ': "'//TRIM(species)//'"'
                        STOP
                     END IF
                     mass(k1:k2) = DOT_PRODUCT(part_counts, boxsim_kind_mass)
                     charge(k1:k2) = charge_int* 1.60217662E-19 !e_c
                     Zatom(k1:k2)        = Z_int
                  END IF
               ELSE
                  mass(k1:k2)         = mass_beams(ibeam)
                  charge(k1:k2)       = charge_beams(ibeam)
                  Zatom(k1:k2)        = Zatom_beams(ibeam)
               END IF

               partvmax            = MAX(partvmax,SQRT(2*E_beams(ibeam)/mass_beams(ibeam)))
               ! Energy distribution
               CALL gauss_rand(n_from_beam, E)
               E = sqrt( (E_beams(ibeam) + E_error*E_beams(ibeam)*E)*(E_beams(ibeam) + E_error*E_beams(ibeam)*E) )
               IF (lbeam_simple) E = E_beams(ibeam)
               weight(k1:k2)       = P_beams(ibeam)/E
               ! Starting Points
               j = Dex_beams(ibeam)
               Xstart = X_BEAMLET(j,Nstart)
               Ystart = Y_BEAMLET(j,Nstart)
               R_start(k1:k2)   = SQRT(Xstart*Xstart+Ystart*Ystart)
               PHI_start(k1:k2) = ATAN2(Ystart,Xstart)
               Z_start(k1:k2)   = Z_BEAMLET(j,Nstart)
               ! Now calculate Divergence (use small angle tan(div)=div here)
               CALL gauss_rand(n_from_beam,Xstart)
               CALL gauss_rand(n_from_beam,Ystart)
               U1 = SQRT(Xstart*Xstart+Ystart*Ystart)
               V1 = ATAN2(Ystart,Xstart)
               Xstart = Div_beams(ibeam)*U1*COS(V1)
               Ystart = Div_beams(ibeam)*U1*SIN(V1)
               ! Calcualte Divergence vectors (assume N are unit vectors)
               U1 = NY_BEAMLET(j,Nstart)
               U2 = -NX_BEAMLET(j,Nstart)
               U3 = 0.0
               V1 =U2*NZ_BEAMLET(j,Nstart)
               V2 =-U1*NZ_BEAMLET(j,Nstart)
               V3 =U1*NY_BEAMLET(j,Nstart)-U2*NX_BEAMLET(j,Nstart)
               ! Starting Velocity 
               vll_start(k1:k2) = SQRT(2*E/mass_beams(ibeam))  ! speed E=0.5*mv^2
               vneut1  = (NX_BEAMLET(j,Nstart) +U1*Xstart + V1*Ystart)*vll_start(k1:k2)
               vneut2  = (NY_BEAMLET(j,Nstart) +U2*Xstart + V2*Ystart)*vll_start(k1:k2)
               vneut3  = (NZ_BEAMLET(j,Nstart) +U3*Xstart + V3*Ystart)*vll_start(k1:k2)
               ! To cylindrical coords
               vr_start(k1:k2)   =  vneut1*COS(PHI_start(k1:k2)) + &
                                    vneut2*SIN(PHI_start(k1:k2))
               vphi_start(k1:k2) = -vneut1*SIN(PHI_start(k1:k2)) + &
                                    vneut2*COS(PHI_start(k1:k2))
               vz_start(k1:k2)   =  vneut3

               END ASSOCIATE
               
               npart_beam_left = npart_beam_left - n_from_beam
               IF (npart_beam_left.EQ.0) THEN
                  ibeam = ibeam+1 
                  npart_beam_left = nparticles_start
               END IF
               n_filled = n_filled + n_from_beam
            END DO
         END DO
         DEALLOCATE(N_start,X_Start,Y_start,Energy, U, V, v_neut)
         DEALLOCATE(X_BEAMLET,Y_BEAMLET,Z_BEAMLET,NX_BEAMLET,NY_BEAMLET,NZ_BEAMLET)
         weight = weight/nparticles_start
      END IF
      DEALLOCATE(mystart_all)

#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_LOCAL, ierr_mpi)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(partvmax,1,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      i = MPI_UNDEFINED
      IF (myid_sharmem == master) i = 0
      CALL MPI_COMM_SPLIT( MPI_COMM_BEAMS,i,myworkid,MPI_COMM_LOCAL,ierr_mpi)
      IF (myid_sharmem == master) THEN
         CALL MPI_BCAST(     R_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(   PHI_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(     Z_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(   vll_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(    mu_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(    vr_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(  vphi_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(    vz_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(       t_end, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(        mass, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(      charge, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(       Zatom, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(      weight, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(        beam, nparticles, MPI_INTEGER, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(lgc2fo_start, nparticles, MPI_LOGICAL, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      END IF
      IF (lboxsim) THEN 
         CALL MPI_BCAST(is_active,      nparticles, MPI_LOGICAL,   master, MPI_COMM_BEAMS,ierr_mpi)
         CALL MPI_BCAST(boxsim_species, nparticles, MPI_CHARACTER, master, MPI_COMM_BEAMS,ierr_mpi)
      END IF

#endif


!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_init_beams_bbnbi
