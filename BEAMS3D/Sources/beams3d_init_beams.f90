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
      USE beams3d_lines, ONLY: nparticles
      USE mpi_params
      USE mpi_inc

!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, i, j, k
!      REAL(rprec)  :: br
      REAL(rprec), ALLOCATABLE :: X(:,:), Y(:,:), Energy(:), R_temp(:)
      REAL(rprec)              :: magZ, magX, magV_neut, magV, xx(3), yy(3), zz(3)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      CALL init_random_seed
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

      ALLOCATE (X(nparticles_start, nbeams), Y(nparticles_start, nbeams), &
                  & Energy(nparticles_start), weight(nparticles_start, nbeams), STAT=ier )
      IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'X, Y, weight', ier)
      !ALLOCATE (R_temp(nparticles_start),STAT=ier)
      !IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'R_temp', ier)
      weight = 0

      IF (myworkid == master) THEN
      nparticles = 0
      DO i=1,nbeams
          CALL gauss_rand(nparticles_start, X(:,i))
          CALL gauss_rand(nparticles_start, Y(:,i))
          X(:,i) = Div_beams(i)*Adist_beams(i)*X(:,i)
          Y(:,i) = Div_beams(i)*Adist_beams(i)*Y(:,i)
          ! Limit distribution so all particles make it through
          !R_temp = SQRT(X(:,i)*X(:,i) + Y(:,i)*Y(:,i))
          !X(:,i) = Asize_beams(i)*X(:,i)/R_temp
          !Y(:,i) = Asize_beams(i)*Y(:,i)/R_temp
          DO j = 1,nparticles_start
              IF ((X(j,i)*X(j,i) + Y(j,i)*Y(j,i)) <= Asize_beams(i)*Asize_beams(i)) THEN
                  weight(j,i)    = 1
                  nparticles = nparticles + 1
              END IF
          END DO
      END DO
      END IF
      !DEALLOCATE(R_temp)
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(nparticles,1,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
!DEC$ ENDIF

      ALLOCATE(   R_start(nparticles), phi_start(nparticles), Z_start(nparticles), vll_start(nparticles), &
                & v_neut(3,nparticles), mass(nparticles), charge(nparticles), Zatom(nparticles), &
                & mu_start(nparticles), t_end(nparticles), &
                & beam(nparticles), STAT=ier   )
      IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'R,phi,Z _start, etc.', ier)

      IF (myworkid == master) THEN
      k = 0
!     Make beam distributions
      DO i=1,nbeams
!         Div_beams is total angle; twice horizontal to edge angle.
          CALL gauss_rand(nparticles_start, Energy)
          Energy = sqrt( (E_beams(i) + E_beams(i)/100*Energy)*(E_beams(i) + E_beams(i)/100*Energy) )
          IF (lbeam_simple) Energy = E_beams(i)
          DO j=1,nparticles_start
              IF (weight(j,i) /= 0) THEN
                  k = k + 1
                  beam(k)         = i
                  mu_start(k)     = 0
                  t_end(k)        = t_end_in(i)
                  mass(k)         = mass_beams(i)
                  charge(k)       = charge_beams(i)
                  Zatom(k)        = Zatom_beams(i)

                  magV_neut         = sqrt( 2*Energy(j)/mass(k) )
                  IF (magV_neut == 0) magV_neut = 1

                  magZ              = sqrt( (r_beams(i,2)*cos(phi_beams(i,2))-r_beams(i,1)*cos(phi_beams(i,1))) &
                                             & *(r_beams(i,2)*cos(phi_beams(i,2))-r_beams(i,1)*cos(phi_beams(i,1))) &
                                          & + (r_beams(i,2)*sin(phi_beams(i,2))-r_beams(i,1)*sin(phi_beams(i,1))) &
                                             & *(r_beams(i,2)*sin(phi_beams(i,2))-r_beams(i,1)*sin(phi_beams(i,1))) &
                                          & + (z_beams(i,2)-z_beams(i,1))*(z_beams(i,2)-z_beams(i,1)) )
                  IF (magZ == 0) magZ = 1

                  magX              = sqrt( ((r_beams(i,2)*cos(phi_beams(i,2))-r_beams(i,1)*cos(phi_beams(i,1))) &
                                                &*(z_beams(i,2)-z_beams(i,1)))*&
                                            ((r_beams(i,2)*cos(phi_beams(i,2))-r_beams(i,1)*cos(phi_beams(i,1))) &
                                                &*(z_beams(i,2)-z_beams(i,1))) + &
                                            ((r_beams(i,2)*sin(phi_beams(i,2))-r_beams(i,1)*sin(phi_beams(i,1))) &
                                                &*(z_beams(i,2)-z_beams(i,1)))*&
                                            ((r_beams(i,2)*sin(phi_beams(i,2))-r_beams(i,1)*sin(phi_beams(i,1))) &
                                                &*(z_beams(i,2)-z_beams(i,1))) + &
                                            (1 - (z_beams(i,2)-z_beams(i,1))*(z_beams(i,2)-z_beams(i,1)))*&
                                            (1 - (z_beams(i,2)-z_beams(i,1))*(z_beams(i,2)-z_beams(i,1))) )
                  IF (magX == 0) magX = 1

                  magV              = sqrt( Adist_beams(i)*Adist_beams(i) + X(j,i)*X(j,i) + Y(j,i)*Y(j,i) )
                  IF (magV == 0) magV = 1

                  R_start(k)       = r_beams(i,1)
                  phi_start(k)     = phi_beams(i,1)
                  Z_start(k)       = z_beams(i,1)

                  zz(1)              = (r_beams(i,2)*cos(phi_beams(i,2))-r_beams(i,1)*cos(phi_beams(i,1)))/magZ
                  zz(2)              = (r_beams(i,2)*sin(phi_beams(i,2))-r_beams(i,1)*sin(phi_beams(i,1)))/magZ
                  zz(3)              = (z_beams(i,2)-z_beams(i,1))/magZ
                  xx(1)              = -(r_beams(i,2)*cos(phi_beams(i,2))-r_beams(i,1)*cos(phi_beams(i,1)))*&
                                       &(z_beams(i,2)-z_beams(i,1))/magX
                  xx(2)              = -(r_beams(i,2)*sin(phi_beams(i,2))-r_beams(i,1)*sin(phi_beams(i,1)))*&
                                       &(z_beams(i,2)-z_beams(i,1))/magX
                  xx(3)              = ( 1 - (z_beams(i,2)-z_beams(i,1))*(z_beams(i,2)-z_beams(i,1)) )/magX
                  yy(1)              = zz(2)*xx(3) - zz(3)*xx(2)
                  yy(2)              = zz(3)*xx(1) - zz(1)*xx(3)
                  yy(3)              = zz(1)*xx(2) - zz(2)*xx(1)
                  !v_neut(1,k)      = magV_neut*(Adist_beams(i)*zz(1) + X(j,i)*xx(1) + Y(j,i)*yy(1))/magV
                  !v_neut(2,k)      = magV_neut*(Adist_beams(i)*zz(2) + X(j,i)*xx(2) + Y(j,i)*yy(2))/magV
                  !v_neut(3,k)      = magV_neut*(Adist_beams(i)*zz(3) + X(j,i)*xx(3) + Y(j,i)*yy(3))/magV
                  v_neut(1,k)      = (Adist_beams(i)*zz(1) + X(j,i)*xx(1) + Y(j,i)*yy(1))/magV
                  v_neut(2,k)      = (Adist_beams(i)*zz(2) + X(j,i)*xx(2) + Y(j,i)*yy(2))/magV
                  v_neut(3,k)      = (Adist_beams(i)*zz(3) + X(j,i)*xx(3) + Y(j,i)*yy(3))/magV
                  v_neut(:,k)      = magV_neut*v_neut(:,k)/SQRT(SUM(v_neut(:,k)**2,DIM=1))

                  vll_start(k) = sqrt( v_neut(1,k)*v_neut(1,k) + v_neut(2,k)*v_neut(2,k) + v_neut(3,k)*v_neut(3,k) )
                  !PRINT *, k, weight(j,i), v_neut(:,k)
              END IF
          END DO
      END DO
      END IF
      DEALLOCATE(X,Y,Energy)
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
      END SUBROUTINE beams3d_init_beams
