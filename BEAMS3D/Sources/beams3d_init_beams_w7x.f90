!-----------------------------------------------------------------------
!     Module:        beams3d_init_beams_w7x
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/22/2016
!     Description:   This subroutine initializes the beam; takes the beam
!                    parameters and specifies a particle distribution.
!                    Need to add error handling.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_init_beams_w7x
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_runtime
      USE beams3d_lines, ONLY: nparticles, partvmax
      USE mpi_params
      USE mpi_inc

!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, i, j, k, k1, k2
!      REAL(rprec)  :: br
      REAL(rprec), ALLOCATABLE :: X(:,:), Y(:,:), block(:,:), Energy(:), X_start(:), Y_start(:)
      REAL(rprec)              :: xbeam(2),ybeam(2),zbeam(2),dxbeam,dybeam,dzbeam,&
                                  dxbeam2,dybeam2,dzbeam2,dlbeam,dxbeam3,dybeam3,dzbeam3
      REAL(rprec)              :: magZ, magX, magV_neut, magV, xx(3), yy(3), zz(3)
      REAL(rprec), PARAMETER   :: X_w7x = .114  !  22.4 duct width
      REAL(rprec), PARAMETER   :: Y_w7x = .2533 ! 50.66 duct height
      REAL(rprec), PARAMETER   :: E_error = .01 ! 1% energy spread
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      CALL init_random_seed
      IF (lverb) THEN
         WRITE(6, '(A)') '----- INITIALIZING W7-X BEAMS -----'
         WRITE(6, '(A,I4)') '      nbeams: ', nbeams
         WRITE(6, '(A,I6)') '      nparticles_start: ', nparticles_start
         CALL FLUSH(6)
      END IF

      ALLOCATE (X(nparticles_start, nbeams), Y(nparticles_start, nbeams), &
                  & Energy(nparticles_start), block(nparticles_start, nbeams), STAT=ier )
      IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'X, Y, weight', ier)

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
         block = 1
      END IF
      nparticles = nparticles_start*nbeams
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(nparticles,1,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
#endif
      ALLOCATE( X_start(nparticles_start), Y_start(nparticles_start), STAT=ier   )
      IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'X,Y_start etc.', ier)
      ALLOCATE(   R_start(nparticles), phi_start(nparticles), Z_start(nparticles), vll_start(nparticles), &
                & v_neut(3,nparticles), mass(nparticles), charge(nparticles), Zatom(nparticles), &
                & mu_start(nparticles), t_end(nparticles), &
                & beam(nparticles), weight(nparticles), STAT=ier   )
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
            partvmax            = MAX(partvmax,6*SQRT(2*E_beams(i)/mass_beams(i))/5.0)
            ! Energy distribution
            CALL gauss_rand(nparticles_start, Energy)
            Energy = sqrt( (E_beams(i) + E_error*E_beams(i)*Energy)*(E_beams(i) + E_error*E_beams(i)*Energy) )
            IF (lbeam_simple) Energy = E_beams(i)
            weight(k1:k2)       = P_beams(i)/Energy
            ! Beamline geometry
            xbeam           = r_beams(i,1:2)*cos(phi_beams(i,1:2))
            ybeam           = r_beams(i,1:2)*sin(phi_beams(i,1:2))
            zbeam           = z_beams(i,1:2)
            dxbeam          = xbeam(2)-xbeam(1)
            dybeam          = ybeam(2)-ybeam(1)
            dzbeam          = zbeam(2)-zbeam(1)
            dlbeam          = SQRT(dxbeam*dxbeam+dybeam*dybeam+dzbeam*dzbeam)
            IF (dlbeam == 0.0) dlbeam = 1
            dxbeam = dxbeam/dlbeam; dybeam=dybeam/dlbeam; dzbeam=dzbeam/dlbeam
            dxbeam2         = dybeam ! normal direction dbeam x hat(z)
            dybeam2         = -dxbeam
            dzbeam2         = 0
            dlbeam          = SQRT(dxbeam2*dxbeam2+dybeam2*dybeam2) ! because dzbeam2=0
            IF (dlbeam == 0.0) dlbeam = 1
            dxbeam2 = dxbeam2/dlbeam; dybeam2=dybeam2/dlbeam; ! because dzbeam2=0
            dxbeam3         = dybeam2*dzbeam                  ! because dzbeam2=0
            dybeam3         =                - dxbeam2*dzbeam ! because dzbeam2=0
            dzbeam3         = dxbeam2*dybeam - dybeam2*dxbeam
            dlbeam          = SQRT(dxbeam3*dxbeam3+dybeam3*dybeam3+dzbeam3*dzbeam3)
            IF (dlbeam == 0.0) dlbeam = 1
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
      DEALLOCATE(X,Y,Energy,X_start,Y_start,block)
      weight = weight/nparticles_start

#if defined(MPI_OPT)
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
#endif


!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_init_beams_w7x
