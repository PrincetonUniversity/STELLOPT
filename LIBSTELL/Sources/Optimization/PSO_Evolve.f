!-----------------------------------------------------------------------
!     Subroutine:    PSO_Evolve
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/13/2013
!     Description:   This subroutine preformes a particle swarm
!                    optimization of a target function (fcn).  The
!                    starting locations of the population are determined
!                    using a random number generator and are located
!                    between XCmin and XCmax.  One member of the
!                    population will start at the input value of x.
!-----------------------------------------------------------------------
      SUBROUTINE PSO_Evolve(fcn, m, n, NP, XCmin, XCmax, x, fvec,
     1                      c1, c2, Vscale, ftol, xtol, maxfev)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds
      USE mpi_params
      USE safe_open_mod
      USE gade_mod, ONLY: pso_cleanup
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      include 'mpif.h'                                       !mpi stuff
!DEC$ ENDIF
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
      LOGICAL :: lkeep_running, lfirst_pass, lnew_global
      INTEGER :: i,j, iproc_min, dex, generation, iunit,
     1           nfeval, ierr, istat, exit_flag, ntot
      REAL(rprec) :: fnorm, fnorm_global, gnorm, pnorm,
     1               vnorm, rand_C1, fnorm_max, xmax, l_temp
      REAL(rprec), ALLOCATABLE :: x_temp(:), fnorm_array(:), 
     1                            x_global(:), x_personal(:,:),
     2                            vel(:), fvec_temp(:), vmax(:),
     3                            fnorm_personal(:)
      REAL(rprec), ALLOCATABLE :: x_array(:,:), fvec_array(:,:)
      CHARACTER(16) :: temp_string
      
      REAL(rprec) ::  enorm
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      ! Preform Allocations
      ALLOCATE (x_temp(n), fnorm_array(NP), fnorm_personal(NP),
     1          x_global(n),
     2          vel(n),fvec_temp(m), vmax(n), stat=ierr)
      ALLOCATE (x_array(n,NP),fvec_array(m,NP),x_personal(n,NP),
     1          stat=ierr)

      ! Calculate VMAX
      xmax = enorm(n,XCmax-XCmin)
      vmax = Vscale*ABS((XCmax - XCmin))
      
      ! Initialize Starting Positions
      CALL RANDOM_SEED ( )  ! Processor reinitializes the seed
      IF (myid .eq. master) THEN
         x_array=0
         x_array(:,1) = x(:)
         DO i=2,NP
            DO j = 1, n
               CALL random_number(rand_C1)
               x_array(j,i)=XCmin(j)+rand_C1*(XCmax(j)-XCmin(j))
            END DO
         END DO
      END IF

!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr)
      IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
      CALL MPI_BCAST(X_array,NP*n,MPI_REAL8,master,MPI_COMM_STEL,
     1               ierr)
      IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
!DEC$ ENDIF

      lkeep_running   = .true.
      lfirst_pass     = .true.
      x_global        = XCmax
      fnorm_global    = 1E30
      fnorm_personal  = 1E30
      iproc_min       = 1
      generation      = 1
      istat           = -1
      nfeval          = 0
      DO WHILE (lkeep_running)
         ! Evaluate population
         IF (lfirst_pass) THEN ! Execute first
            ! Have proc 0 evaluate
            IF (myid .eq. master) THEN
               x_temp        = x_array(:,1)
               fvec_temp(:)  = 0
               CALL fcn(m, n, x_temp, fvec_temp, istat, generation)
               fvec_array(:,1) = fvec_temp
               fnorm            = SQRT(SUM(fvec_temp*fvec_temp))
               WRITE(6,'(70("="))')
               WRITE(6,*) '   MEMBER      PROC          CHISQ'
               WRITE(6,'(70("="))')
               WRITE(6,'(2X,I6,8X,I3,7X,1ES12.4E2)') 1,0,fnorm**2
               CALL FLUSH(6)
            END IF
!DEC$ IF DEFINED (MPI_OPT)
            CALL MPI_BCAST(istat,1,MPI_INTEGER,0,MPI_COMM_STEL,ierr_mpi)
            CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)
!DEC$ ENDIF
            IF (istat < 0) THEN
               IF (myid == 0) WRITE(6,'(A)') 
     1                      '!!!!! First Eq failed !!!!!'
               DEALLOCATE(x_temp,fvec_temp)
               IF (ALLOCATED(x_array)) DEALLOCATE(x_array)
               IF (ALLOCATED(fvec_array)) DEALLOCATE(fvec_array)
               RETURN
            END IF

            ! Cleanup after first iteration
            istat = -100
            IF (myid .eq. 0) istat = -101
            i=0
            call fcn (m, n, x_temp, fvec_temp, istat, i)

            ntot = NP-1
            i=2; j=NP
         ELSE
            ntot = NP
            i=1; j=NP

            ! Write Generation Number to screeen
            IF (myid == master) THEN
!               WRITE(6,'(A)') '   '
               WRITE(6,'(A,I5)') 
     1                      '---  Starting Generation',generation
!               WRITE(6,'(A)') '   '
            END IF
         END IF
         lfirst_pass = .false.

         ! Evaluate
         CALL eval_x_queued(fcn,m,n,ntot,x_array(1:n,i:j),
     1                  fvec_array(1:m,i:j),nfeval,MPI_COMM_STEL)
         nfeval = nfeval + NP

!DEC$ IF DEFINED (MPI_OPT)
         CALL MPI_BCAST(fvec_array,NP*m,MPI_REAL8,master,MPI_COMM_STEL,
     1               ierr)
         IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
         CALL MPI_BCAST(x_array,NP*n,MPI_REAL8,master,MPI_COMM_STEL,
     1               ierr)
         IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
!DEC$ ENDIF
 
         ! Output the xvector and fevals
         IF (myid == master) THEN
            ! XVEC
            CALL safe_open(iunit,istat,'xvec.dat','unknown',
     1                     'formatted', ACCESS_IN='APPEND')
            DO i = 1, NP
               WRITE(iunit,'(2(2X,I5.5))') n,generation
               fnorm = SQRT(SUM(fvec_array(:,i)*fvec_array(:,i),DIM=1))
               WRITE(iunit,'(10ES22.12E3)') x_array(1:n,i)
               WRITE(iunit,'(ES22.12E3)') fnorm
            END DO
            CLOSE(iunit)
            ! FEVALS
            WRITE(temp_string,'(i6.6)') generation
            iunit = 27; j=0;
            CALL safe_open(iunit,istat,'fevals.'//TRIM(temp_string),
     1                     'new','formatted')
            WRITE(iunit,'(2X,i6,2X,i6)') m,n
            WRITE(iunit,'(1p,4e22.14)') (fvec(i), i=1,m)
            DO i = 1, m
               WRITE(iunit,'(1p,4e22.14)') (fvec_array(i,j), j=1,NP)
            END DO
            CLOSE(iunit)
         END IF

         ! Update values
         fnorm_array = SQRT(SUM(fvec_array*fvec_array,DIM=1))
         lnew_global = .false.
         iproc_min = -1
         DO  i=1, NP
            ! Local
            IF (fnorm_array(i) < fnorm_personal(i)) THEN
               fnorm_personal(i) = fnorm_array(i)
               x_personal(:,i) = x_array(:,i)
            END IF
            ! Global
            IF (fnorm_array(i) < fnorm_global) THEN
               ! Check for exit
               IF ((fnorm_global-fnorm_array(i))/fnorm_global
     1                  < ftol) THEN
                  lkeep_running = .false.
                  exit_flag = 1
               END IF
               IF (MAXVAL(ABS((x_global-x_array(:,i))/x_global))
     1                          < xtol) THEN
                  lkeep_running = .false.
                  exit_flag = 2
               END IF
               ! Update values
               x_global = x_array(:,i)
               fnorm_global = fnorm_array(i)
               lnew_global = .true.
               iproc_min = i
            END IF
         END DO

         ! Print Minimum message
         IF (lnew_global) THEN
            x_temp = x_global
            IF (myid == master) THEN
               istat = myid+1
               CALL fcn(m, n, x_temp, fvec_temp, istat, generation)
               WRITE(6,'(A)') '   '
               WRITE(6,'(A,I5,A,1ES12.4E3)') 
     1                      '---  New global minimum at generation',
     2                   generation,' Chi-Squared:',fnorm_global**2
               WRITE(6,'(A)') '   '
            END IF

            ! Cleanup after first iteration
            istat = -100
            IF (myid .eq. master) istat = pso_cleanup
            call fcn (m, n, x_temp, fvec_temp, istat, generation)
         END IF
         
         ! Move particles
         DO i = 1, NP
            ! Calculate Velocity
            CALL random_number(rand_C1)
            vel = vel + c1*(x_personal(:,i)-x_array(:,i))*rand_C1
     1                + c2*(x_global - x_array(:,i))*rand_C1

            ! Bound Velocity
            WHERE (ABS(vel)>vmax) vel = SIGN(vmax,vel)

            ! Step Particle
            x_array(:,i) = x_array(:,i) + vel

            ! Kick minimum particle
            IF (i==iproc_min) THEN
               DO j = 1, n
                  CALL random_number(rand_C1)
                  x_array(j,i)=XCmin(j)+rand_C1*(XCmax(j)-XCmin(j))
               END DO
               x_personal(:,iproc_min) = x_array(:,iproc_min)
               fnorm_personal(i) = 1E30
            END IF

            ! Bound the particle
            x_temp = x_array(:,i)
            WHERE (x_temp > XCmax) x_temp = XCmax
            WHERE (x_temp < XCmin) x_temp = XCmin
            x_array(:,i) = x_temp
         END DO

         ! Test for max 
         IF (generation == maxfev) THEN
            lkeep_running = .false.
            exit_flag = 3
         END IF

         ! Update Generation Number
         generation = generation + 1
      END DO

      ! PRINT Message
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
      fvec = fvec_array(:,iproc_min)
      DEALLOCATE (x_temp,fnorm_array,fnorm_personal,x_global,vel,
     1            fvec_temp)
      DEALLOCATE (x_array, fvec_array, x_personal)
      
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE PSO_Evolve
