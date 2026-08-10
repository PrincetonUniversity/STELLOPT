!-----------------------------------------------------------------------
!     Module:        beams3d_init_beams
!     Authors:       M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          07/03/2012
!     Description:   This subroutine initializes the beam; takes the beam
!                    parameters and specifies a particle distribution.
!                    Need to add error handling.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_init_beams
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
      INTEGER :: ier, i, j, k, k1, k2, MPI_COMM_LOCAL
      INTEGER, DIMENSION(:), ALLOCATABLE :: N_start
      REAL(rprec) :: rtemp, nx, ny ,nz
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: Energy, X_start, Y_start
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: X, Y, U, V, v_neut
      REAL(rprec), PARAMETER   :: E_error = .01 ! 1% energy spread
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      IF (lverb) THEN
         WRITE(6, '(A)') '----- INITIALIZING BEAMS -----'
         WRITE(6, '(A,I4)') '      nbeams: ', nbeams
         WRITE(6, '(A,I6)') '      nparticles_start: ', nparticles_start
         DO i=1, nbeams
            WRITE(6, '(A,I2,A,I4,A,A,I2,A,F7.3,A)') '            E_BEAM(',i ,'): ',&
                     NINT(E_beams(i)*6.24150636309E15),' [keV]',& !(1.0E-3/ec)
                                     ' P_BEAM(',i ,'): ',&
                     (P_beams(i)*1E-6),' [MW]'
         END DO
         CALL FLUSH(6)
      END IF
      IF (lascot) THEN
         CALL beams3d_write_ascoth5('BBNBI_BASIC')
         IF (lverb) WRITE(6, '(A,I4)')      '   ASCOT5 File: Updated'
      END IF

      IF (myid_sharmem == master) lgc2fo_start(:) = .TRUE.
      IF (myid_sharmem == master) THEN
         ALLOCATE(N_start(nparticles_start),X_start(nparticles_start),Y_start(nparticles_start),&
                  Energy(nparticles_start), U(3,nparticles_start), V(3,nparticles_start), &
                  v_neut(3,nparticles_start))
         k1 = 1; k2 = nparticles_start
         DO i = 1, nbeams
            ! Beam Geometry
            nx                  = R_beams(i,2)*cos(PHI_beams(i,2))-R_beams(i,1)*cos(PHI_beams(i,1))
            ny                  = R_beams(i,2)*sin(PHI_beams(i,2))-R_beams(i,1)*sin(PHI_beams(i,1))
            nz                  = Z_beams(i,2)-Z_beams(i,1)
            rtemp               = 1.0/sqrt(nx*nx+ny*ny+nz*nz)
            nx                  = nx*rtemp
            ny                  = ny*rtemp
            nz                  = nz*rtemp
            ! Basics
            beam(k1:k2)         = i
            mu_start(k1:k2)     = 0
            t_end(k1:k2)        = t_end_in(i)
            mass(k1:k2)         = mass_beams(i)
            charge(k1:k2)       = charge_beams(i)
            Zatom(k1:k2)        = Zatom_beams(i)
            partvmax            = MAX(partvmax,SQRT(2*E_beams(i)/mass_beams(i)))
            ! Energy distribution
            CALL gauss_rand(nparticles_start, Energy)
            Energy = sqrt( (E_beams(i) + E_error*E_beams(i)*Energy)*(E_beams(i) + E_error*E_beams(i)*Energy) )
            IF (lbeam_simple) Energy = E_beams(i)
            weight(k1:k2)       = P_beams(i)/Energy
            R_start(k1:k2)      = R_beams(i,1)
            PHI_start(k1:k2)    = PHI_beams(i,1)
            Z_start(k1:k2)      = Z_beams(i,1)
            ! Now calculate Divergence (use small angle tan(div)=div here)
            CALL gauss_rand(nparticles_start,X_start)
            CALL gauss_rand(nparticles_start,Y_start)
            U(1,:) = SQRT(X_start*X_start+Y_start*Y_start)
            V(1,:) = ATAN2(Y_start,X_start)
            X_start = Div_beams(i)*U(1,:)*COS(V(1,:))
            Y_start = Div_beams(i)*U(1,:)*SIN(V(1,:))
            ! Calcualte Divergence vectors (assume N are unit vectors)
            U(1,:) = ny
            U(2,:) = -nx
            U(3,:) = 0.0
            V(1,:) = U(2,:)*nz
            V(2,:) =-U(1,:)*nz
            V(3,:) = U(1,:)*ny-U(2,:)*nx
            ! Starting Velocity 
            vll_start(k1:k2) = SQRT(2*Energy/mass_beams(i))  ! speed E=0.5*mv^2
            v_neut(1,:)  = (nx + U(1,:)*X_Start + V(1,:)*Y_start)*vll_start(k1:k2)
            v_neut(2,:)  = (ny + U(2,:)*X_Start + V(2,:)*Y_start)*vll_start(k1:k2)
            v_neut(3,:)  = (nz + U(3,:)*X_Start + V(3,:)*Y_start)*vll_start(k1:k2)
            ! To cylindrical coords
            vr_start(k1:k2)   =  v_neut(1,:)*COS(PHI_start(k1:k2)) + &
                                 v_neut(2,:)*SIN(PHI_start(k1:k2))
            vphi_start(k1:k2) = -v_neut(1,:)*SIN(PHI_start(k1:k2)) + &
                                 v_neut(2,:)*COS(PHI_start(k1:k2))
            vz_start(k1:k2)   =  v_neut(3,:)
            k1 = k2 + 1
            k2 = k2 + nparticles_start
         END DO
         DEALLOCATE(N_start,X_Start,Y_start,Energy, U, V, v_neut)
         weight = weight/nparticles_start

      END IF

#if defined(MPI_OPT)
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
#endif


!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_init_beams
