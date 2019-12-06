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
!      REAL(rprec)  :: br
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: Energy, X_start, Y_start
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: X, Y, X_BEAMLET, Y_BEAMLET, Z_BEAMLET, &
                                                  NX_BEAMLET, NY_BEAMLET, NZ_BEAMLET
      REAL(rprec)              :: xbeam(2),ybeam(2),zbeam(2),dxbeam,dybeam,dzbeam,&
                                  dxbeam2,dybeam2,dzbeam2,dlbeam,dxbeam3,dybeam3,dzbeam3
      REAL(rprec)              :: magZ, magX, magV_neut, magV, xx(3), yy(3), zz(3)
      REAL(rprec), PARAMETER   :: X_w7x = .114  !  22.4 duct width
      REAL(rprec), PARAMETER   :: Y_w7x = .2533 ! 50.66 duct height
      REAL(rprec), PARAMETER   :: E_error = .01 ! 1% energy spread

      ! For HDF5
      CHARACTER(256)           :: beamlet_file
      INTEGER(HID_T)           :: h5_fid, h5_did, h5_sid
      INTEGER(HSIZE_T), DIMENSION(2)    :: dims, maxdims

!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      CALL init_random_seed

      ! Open and process the HDF5 file containing the data
      beamlet_file = 'w7x_NI20NI21_beamlet_geo.h5'
      CALL H5open_f(ier)
      CALL h5fopen_f(beamlet_file, H5F_ACC_RDONLY_F, h5_fid, ier)
      CALL h5dopen_f (h5_fid, 'X_BEAM', h5_did, ier)
      CALL h5dget_space_f(h5_did, h5_sid ,ier)
      CALL h5sget_simple_extent_dims_f(h5_sid, dims, maxdims,ier)
      PRINT *,dims, maxdims
      ALLOCATE(X_BEAMLET(dims(1),dims(2)), Y_BEAMLET(dims(1),dims(2)), Z_BEAMLET(dims(1),dims(2)),&
         NX_BEAMLET(dims(1),dims(2)), NY_BEAMLET(dims(1),dims(2)), NZ_BEAMLET(dims(1),dims(2)))
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

      IF (lverb) THEN
         WRITE(6, '(A)') '----- INITIALIZING BEAMLET BASED BEAMS -----'
         WRITE(6, '(A,I4)') '      nbeams: ', nbeams
         WRITE(6, '(A,I6)') '      nparticles_start: ', nparticles_start
         CALL FLUSH(6)
      END IF

      ALLOCATE (X(nparticles_start, nbeams), Y(nparticles_start, nbeams), &
                  & Energy(nparticles_start), weight(nparticles_start, nbeams), STAT=ier )
      IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'X, Y, weight', ier)
      weight = 0

      IF (myworkid == master) THEN
         nparticles = 0
         DO i=1,nbeams
            CALL gauss_rand(nparticles_start, X(:,i))
            CALL gauss_rand(nparticles_start, Y(:,i))
            X(:,i) = X(:,i)/MAXVAL(ABS(X(:,i)),DIM=1)
            Y(:,i) = Y(:,i)/MAXVAL(ABS(Y(:,i)),DIM=1)
         END DO
         ! Renormalize to beam box size
         X = X*X_w7x  ! Width
         Y = Y*Y_w7x  ! Height
         weight = 1
      END IF
      nparticles = nparticles_start*nbeams
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(nparticles,1,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
!DEC$ ENDIF
      ALLOCATE( X_start(nparticles_start), Y_start(nparticles_start), STAT=ier   )
      IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'X,Y_start etc.', ier)
      ALLOCATE(   R_start(nparticles), phi_start(nparticles), Z_start(nparticles), vll_start(nparticles), &
                & v_neut(3,nparticles), mass(nparticles), charge(nparticles), Zatom(nparticles), &
                & mu_start(nparticles), t_end(nparticles), &
                & beam(nparticles), STAT=ier   )
      IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'R,phi,Z _start, etc.', ier)

      ! Handle beam distrbution
      IF (myworkid == master) THEN
         k1 = 1; k2 = nparticles_start
         DO i=1,nbeams
            IF (lverb) WRITE(6, '(A,I2,A,I4,A,A,I2,A,F7.3,A)') '            E_BEAM(',i ,'): ',&
                     NINT(E_beams(i)*6.24150636309E15),' [keV]',& !(1.0E-3/ec)
                                     ' P_BEAM(',i ,'): ',&
                     (P_beams(i)*1E-6),' [MW]'
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
            ! Beamline geometry
            xbeam           = r_beams(i,1:2)*cos(phi_beams(i,1:2))
            ybeam           = r_beams(i,1:2)*sin(phi_beams(i,1:2))
            zbeam           = z_beams(i,1:2)
            dxbeam          = xbeam(2)-xbeam(1)
            dybeam          = ybeam(2)-ybeam(1)
            dzbeam          = zbeam(2)-zbeam(1)
            dlbeam          = SQRT(dxbeam*dxbeam+dybeam*dybeam+dzbeam*dzbeam)
            IF (dlbeam == 0) dlbeam = 1
            dxbeam = dxbeam/dlbeam; dybeam=dybeam/dlbeam; dzbeam=dzbeam/dlbeam
            dxbeam2         = dybeam ! normal direction dbeam x hat(z)
            dybeam2         = -dxbeam
            dzbeam2         = 0
            dlbeam          = SQRT(dxbeam2*dxbeam2+dybeam2*dybeam2) ! because dzbeam2=0
            IF (dlbeam == 0) dlbeam = 1
            dxbeam2 = dxbeam2/dlbeam; dybeam2=dybeam2/dlbeam; ! because dzbeam2=0
            dxbeam3         = dybeam2*dzbeam                  ! because dzbeam2=0
            dybeam3         =                - dxbeam2*dzbeam ! because dzbeam2=0
            dzbeam3         = dxbeam2*dybeam - dybeam2*dxbeam
            dlbeam          = SQRT(dxbeam3*dxbeam3+dybeam3*dybeam3+dzbeam3*dzbeam3)
            IF (dlbeam == 0) dlbeam = 1
            dxbeam3 = dxbeam3/dlbeam; dybeam3=dybeam3/dlbeam; dzbeam3=dzbeam3/dlbeam
            ! Starting Points
            X_start          = xbeam(1) + dxbeam2*X(:,i) + dxbeam3*Y(:,i)
            Y_start          = ybeam(1) + dybeam2*X(:,i) + dybeam3*Y(:,i)
            R_start(k1:k2)   = SQRT(X_start*X_start+Y_start*Y_start)
            PHI_start(k1:k2) = ATAN2(Y_start,X_start)
            Z_start(k1:k2)   = zbeam(1) + dzbeam3*Y(:,i)                   ! because dzbeam2=0
            ! Starting Velocity
            vll_start(k1:k2) = SQRT(2*Energy/mass_beams(i))  ! speed E=0.5*mv^2
            v_neut(1,k1:k2)  = dxbeam*vll_start(k1:k2)
            v_neut(2,k1:k2)  = dybeam*vll_start(k1:k2)
            v_neut(3,k1:k2)  = dzbeam*vll_start(k1:k2)
            k1 = k2 + 1
            k2 = k2 + nparticles_start
         END DO
         WHERE(PHI_start < 0) PHI_start = PHI_start+pi2
      END IF

      DEALLOCATE(X,Y,Energy,X_start,Y_start)
      DEALLOCATE(X_BEAMLET,Y_BEAMLET,Z_BEAMLET,NX_BEAMLET,NY_BEAMLET,NZ_BEAMLET)
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(mu_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(t_end,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(mass,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(charge,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(Zatom,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(weight,nparticles_start*nbeams,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
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
