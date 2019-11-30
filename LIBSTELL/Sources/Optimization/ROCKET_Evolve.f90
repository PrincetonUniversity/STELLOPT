!-----------------------------------------------------------------------
!     Subroutine:    ROCKET_Evolve
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/13/2013
!     Description:   This subroutine is similar to the particle swarm
!                    however all members of the population start at
!                    an initial point and move away from that point
!                    with some initial velocity at the start.
!-----------------------------------------------------------------------
      SUBROUTINE ROCKET_Evolve(fcn, m, n, NP, XCmin, XCmax, x, fvec,&
                            c1, c2, Vscale, ftol, xtol, maxfev)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds
      USE mpi_params
      USE safe_open_mod
      USE mpi_inc
      IMPLICIT NONE
!-----------------------------------------------------------------------
!     Input Variables
!        fcn     Target function one wishes to mimizie
!        m       Dimensionality of the function
!        n       Dimensionality of the function space
!        NP      Population size
!        XCmin   Minimum values of function space
!        XCmax   Maximum values of function space
!        x       Vector specifying mimium location in function space
!        fvec    Vector of function values at minimum
!        c1      Coefficient scaling local acceleration
!        c2      Coefficient scaling global acceleration
!        Vscale  Scaling factor for maximum velocity
!        ftol    Tollerance for global minimum in terms of fnorm
!        xtol    Tollerance for global minimum in terms of x
!        maxfev  Maximum number of function evaluations.
!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: m, n, NP, maxfev
      REAL(rprec), INTENT(in)  :: c1, c2, ftol, xtol, Vscale
      REAL(rprec), INTENT(in)  :: XCmin(n), XCmax(n)
      REAL(rprec), INTENT(inout) :: x(n)
      REAL(rprec), INTENT(out) :: fvec(m)
      EXTERNAL  fcn
      
!-----------------------------------------------------------------------
!     Local Variables
!        lkeep_running  Logical to control continued iteration
!        i,j            Dummy index
!        iproc_min      Index of processor with minimum fcn value
!        dex            Index of X-vector to use
!        generation     Current iteration of particle steps
!        iunit          File unit number
!        update_global  Flag to let que processors to update global min
!        nfeval         Counter for function evaluations
!        ierr           Error flag
!        istat          Status flag
!        exit_flag      Flag to control printing of exit message
!        fnorm          Function normalization SQRT(F.F)
!        fnorm_global   Global Min. function normalization
!        fnorm_personal Local Min. function normalization
!        fnorm_max      Maximum fnorm value encountered
!        x_max          Maximum distance across domain
!        gnorm          Distance to global min
!        pnorm          Distance to local min
!        vnorm          Velocity normalization
!        rand_C1        Random number
!        x_temp         Vector of input variables
!        fnorm_array    Vector of function values at (x_temp)
!        x_global       Vector of input variables (global min)
!        x_personal     Vector of input variables (local min)
!        vel            Velocity Vector
!        temp_fvec      Helper array
!        vmax           Array of maximum velocities
!        x_array        Array of all X's for the population
!----------------------------------------------------------------------
      LOGICAL :: lkeep_running
      INTEGER :: i,j, iproc_min, dex, generation, iunit, update_global, &
                 nfeval, ierr, istat, exit_flag
      REAL(rprec) :: fnorm, fnorm_global, fnorm_personal, gnorm, pnorm, &
                     vnorm, rand_C1, fnorm_max, xmax, l_temp, vmax_norm
      REAL(rprec), ALLOCATABLE :: x_temp(:), fnorm_array(:), &
                                  x_global(:), x_personal(:), &
                                  vel(:), temp_fvec(:), vmax(:), &
                                  dx_personal(:), dx_global(:)
      REAL(rprec), ALLOCATABLE :: x_array(:,:), fval_array(:,:)
      
      REAL(rprec) ::  enorm
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      ! Preform Allocations
      ALLOCATE (x_temp(n), fnorm_array(NP), &
                x_global(n), x_personal(n), &
                vel(n), temp_fvec(m), vmax(n), & 
                dx_personal(n), dx_global(n), stat=ierr)
      ALLOCATE (x_array(NP,n),fval_array(NP,m), stat=ierr)

      ! Calculate VMAX
      xmax = enorm(n,XCmax-XCmin)
      vmax = Vscale*ABS((XCmax - XCmin))
      vmax_norm = enorm(n,vmax)
      
      ! Initialize Starting Positions
      CALL RANDOM_SEED ( )  ! Processor reinitializes the seed
      x_array = 0
      DO i=1, NP
         x_array(i,:) = x(:)
      END DO
      DO j=1,n
         CALL random_number(rand_C1)
         vel(j) = 0.001*rand_C1*(XCmax(j)-XCmin(j))
         IF ((ABS(vel(i)) > vmax(i))) vel(i) = sign(vmax(i),vel(i))
      END DO
      vel = vmax_norm * 0.01 * vel / enorm(n,vel)
      x_temp = x_temp + vel
      
      ! Initialize variables
      lkeep_running = .true.
      dex = myid + 1
      nfeval = 0
      generation = 0
      exit_flag  = 0
      vel = 0.0
      fval_array = 0.0
      fnorm_max      = 0.0
      fnorm_personal = 1.0E30
      fnorm_global   = 1.0E30
      fnorm_array    = 1.0E30
      x_global       = x_array(dex,:)
      x_personal     = x_array(dex,:)
      x_temp         = x_array(dex,:)
      ! Iterate
      DO WHILE (lkeep_running)
         ! Evaluate Equilibria
         nfeval = nfeval + dex
         temp_fvec = 0.0
         istat = dex
         IF (myid == master .and. generation == 0) istat = -1
         CALL fcn(m, n, x_temp, temp_fvec, istat, nfeval)
         fnorm = enorm(m,temp_fvec)
         
         ! Default a couple things on first pass and write to screen
         IF (generation == 0)  THEN
            x_personal = x_temp
            fnorm_personal = fnorm
            IF (myid == master) THEN
               WRITE(6,'(70("="))')
               WRITE(6,*) '      NFUNCT      PROC         CHISQ', &
                          '           CHI_BEST'
               WRITE(6,'(70("="))')
            END IF
            CALL MPI_BARRIER(MPI_COMM_STEL, ierr)
            IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
            WRITE(6,'(2(6X,I5),7X,1ES12.4E3)') nfeval,myid,fnorm**2
         ELSE
            IF (fnorm < fnorm_personal) x_personal = x_temp
            IF (fnorm < fnorm_personal) fnorm_personal = fnorm
            WRITE(6,'(2(6X,I5),3(7X,1ES12.4E3))') &
                             nfeval,myid,fnorm**2,fnorm_personal**2
         END IF
         CALL FLUSH(6)
         
         ! Find the global minimum state
         update_global = 0
         fnorm_array   = 1.0E30
!DEC$ IF DEFINED (MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_STEL, ierr)
         IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
         CALL MPI_ALLGATHER(fnorm, 1, MPI_REAL8, fnorm_array, &
                                1, MPI_REAL8, MPI_COMM_STEL, ierr)
         IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
         iproc_min = MINLOC(fnorm_array, DIM=1)
         IF (myid == (iproc_min-1)) THEN
            IF (fnorm < fnorm_global) THEN
               l_temp = enorm(n,x_temp-x_global)
               IF (ABS(fnorm-fnorm_global)/fnorm_global <= ftol) THEN
                  exit_flag = 1
                  lkeep_running = .false.
               END IF
               IF (l_temp/xmax <= xtol) THEN
                  exit_flag = 2
                  lkeep_running = .false.
               END IF
               IF (generation == 0) THEN ! do this so we don't auto-exit first iteration
                  exit_flag = 0
                  lkeep_running = .true.
               END IF
               update_global = 1
               fnorm_global = fnorm
               x_global = x_temp
               istat = -102
               CALL fcn(m, n, x_temp, temp_fvec, istat, generation)  ! Save the output
               WRITE(6,'(A)') '   '
               WRITE(6,'(A,I5,A,1ES12.4E3)') &
                            '---  New global minimum at iteration', &
                               nfeval,' Chi-Squared:',fnorm_global**2
               WRITE(6,'(A)') '   '
            END IF
         END IF
         CALL MPI_BARRIER(MPI_COMM_STEL, ierr)
         IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
         CALL MPI_BCAST(update_global,1,MPI_INTEGER,iproc_min-1, &
                        MPI_COMM_STEL, ierr)
         IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
         IF (update_global /= 0) THEN
            CALL MPI_BCAST(exit_flag,1,MPI_INTEGER,iproc_min-1, &
                        MPI_COMM_STEL, ierr)
            IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
            CALL MPI_BCAST(x_global,n,MPI_REAL8,iproc_min-1, &
                           MPI_COMM_STEL,ierr)
            IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
            CALL MPI_BCAST(fnorm_global,1,MPI_REAL8,iproc_min-1, &
                           MPI_COMM_STEL,ierr)
            IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
            CALL MPI_BCAST(temp_fvec,m,MPI_REAL8,iproc_min-1, &
                           MPI_COMM_STEL,ierr)
            IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
            CALL MPI_BCAST(lkeep_running,1,MPI_LOGICAL,iproc_min-1, &
                           MPI_COMM_STEL,ierr)
            IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
         END IF
!DEC$ ELSE
!DEC$ ENDIF
         
         ! Output the xvectors
         DO i = 0, numprocs-1
            IF (myid == i) THEN
               CALL safe_open(iunit,istat,'xvec.dat','unknown', &
                          'formatted', ACCESS_IN='APPEND')
               WRITE(iunit,'(2(2X,I5.5))') n,generation
               WRITE(iunit,'(10ES22.12E3)') x_temp(1:n)
               WRITE(iunit,'(ES22.12E3)') fnorm
               CLOSE(iunit)
            END IF
!DEC$ IF DEFINED (MPI_OPT)
            CALL MPI_BARRIER(MPI_COMM_STEL,ierr)
            IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
!DEC$ ENDIF
         END DO
         
         ! Calculate the Velocity
         dx_global = x_temp-x_global
         dx_personal = x_temp - x_personal
         l_temp = enorm(n,dx_global)
         IF (l_temp > 0) vel = vel - c2*dx_global/(l_temp*l_temp*l_temp)
         l_temp = enorm(n,dx_personal)
         IF (l_temp > 0) vel = vel - c1*dx_personal/(l_temp*l_temp*l_temp)
         DO i = 1, n
            IF ((ABS(vel(i)) > vmax(i))) vel(i) = sign(vmax(i),vel(i))
         END DO
         
         ! Kick to new velocity vector
         !IF ((myid == (iproc_min-1)) .and. (update_global/=0)) THEN
         !   DO i = 1, n
         !      CALL random_number(rand_C1)
         !      x_temp(i) = XCmin(i)+rand_C1*(XCmax(i)-XCmin(i))
         !      vel(i)    = 0.0
         !   END DO
         !   x_personal = x_temp
         !   fnorm_personal = 1.0E30
         !END IF
         
         ! Now advance the particle
         x_temp = x_temp + vel
         
         ! Check for bound limits
         ! Set to wall value and make sure we move away from wall (bounce)
         DO i =1, n
            IF (x_temp(i) > XCmax(i)) THEN
               x_temp(i) = XCmax(i)
               vel(i)    = - 0.1 * ABS(vel(i))
            END IF
            IF (x_temp(i) < XCmin(i)) THEN
               x_temp(i) = XCmin(i)
               vel(i)    = 0.1 * ABS(vel(i))
            END IF
         END DO
         !X_array(dex,:) = x_temp

!DEC$ IF DEFINED (MPI_OPT)
         ! Update nfeval
         CALL MPI_BARRIER(MPI_COMM_STEL,ierr)
         IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
         CALL MPI_BCAST(nfeval,1,MPI_INTEGER,numprocs-1, &
                       MPI_COMM_STEL, ierr)
         IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
!DEC$ ENDIF
         generation = generation + 1
         
         ! Test for exit condition
         IF (nfeval >= maxfev) THEN
            lkeep_running = .false.
            exit_flag = 3
         END IF
      END DO
      IF (myid == master) THEN
         WRITE(6,*) ' '
         IF (exit_flag == 1) THEN
            WRITE(6,'(7X,A)') 'Execution Terminated FTOL Achieved'
         ELSE IF (exit_flag == 2) THEN
            WRITE(6,'(7X,A)') 'Execution Terminated XTOL Achieved'
         ELSE IF (exit_flag == 3) THEN
            WRITE(6,'(7X,A)') 'Execution Terminated MAXFEV Achieved'
         ELSE
            WRITE(6,'(7X,A)') 'Execution Terminated for unknown reason'
         END IF
         WRITE(6,*) ' '
      END IF
      x = x_global
      fvec = temp_fvec
      DEALLOCATE (x_temp,fnorm_array,x_global,x_personal,vel,temp_fvec)
      DEALLOCATE (x_array, fval_array)
      
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE ROCKET_Evolve
