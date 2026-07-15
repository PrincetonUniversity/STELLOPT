!======================================================================
! SUBROUTINE: SA_Evolve
! DESCRIPTION: Parallel Population-based Simulated Annealing optimization
!              designed to interface with the LIBSTELL 'eval_x_queued' 
!              subroutine framework.
!
! GENERATION NOTE: 
! This subroutine was scaffolding-generated and co-authored by Gemini 
! (An AI collaborator by Google), adapting standard Differential Evolution 
! interfaces (DE2_Evolve) to a parallelized Simulated Annealing heuristic.
!
! DATE GENERATED: June 2026
!======================================================================
SUBROUTINE SA_Evolve(fcn, m, n, NP, XCmin, XCmax, x, fvec, &
                     maxfev, F_XC, CR_XC, strategy, CR_strategy, &
                     iWRITE, iRESTART, lrestart)
   USE mpi_params
   USE safe_open_mod
   USE gade_mod, ONLY: gade_cleanup
   USE mpi_inc
   IMPLICIT NONE

   ! --- Arguments ---
   EXTERNAL :: fcn
   INTEGER, INTENT(IN) :: m, n, NP
   INTEGER, INTENT(INOUT) :: maxfev
   REAL(8), INTENT(IN) :: XCmin(n), XCmax(n)
   REAL(8), INTENT(INOUT) :: x(n), fvec(m)
   REAL(8), INTENT(IN) :: F_XC         ! SA Use: Initial Temperature (T0)
   REAL(8), INTENT(IN) :: CR_XC        ! SA Use: Cooling Rate (e.g., 0.95)
   INTEGER, INTENT(IN) :: strategy     ! SA Use: Perturbation strategy (not used)
   INTEGER, INTENT(IN) :: CR_strategy  ! SA Use: Cooling schedule strategy (not used)
   INTEGER, INTENT(IN) :: iWRITE       ! Not used
   INTEGER, INTENT(IN) :: iRESTART     ! For reading/writing restart file.
   LOGICAL, INTENT(IN) :: lrestart

   ! --- Local Variables ---
   INTEGER :: i, j, iter, iunitx, ierr, iflag, ibest
   REAL(8) :: T
   REAL(8), ALLOCATABLE :: x_new(:,:), fvec_new(:,:), x_old(:,:), fvec_old(:,:)
   REAL(8), ALLOCATABLE :: cost_old(:), cost_new(:)
   REAL(8) :: rand_val, delta_cost, step_amp, fvecbest, rand_c1
    
   ! --- Interface for the STELLOPT queued evaluator ---
   ! (Adjust the argument order/types to match your exact LIBSTELL branch)
   INTERFACE
      SUBROUTINE eval_x_queued(fcn, m, n, NP, xvec, fvec, iter, HYPER_COMM)
         IMPLICIT NONE
         EXTERNAL :: fcn
         INTEGER, INTENT(IN) :: m, n, NP
         REAL(8), INTENT(IN) :: xvec(n, NP)
         REAL(8), INTENT(OUT) :: fvec(m, NP)
         INTEGER, INTENT(INOUT) :: iter
         INTEGER, INTENT(INOUT) :: HYPER_COMM
      END SUBROUTINE eval_x_queued
   END INTERFACE

   ! --- Memory Allocation ---
   ALLOCATE(x_new(n, NP))
   ALLOCATE(fvec_new(m, NP))
   ALLOCATE(x_old(n, NP))
   ALLOCATE(fvec_old(m, NP))
   ALLOCATE(cost_old(NP))
   ALLOCATE(cost_new(NP))

   fvecbest = 1.0E20
   iter = 0
   T = F_XC  ! Set initial temperature
   CALL RANDOM_SEED()
      
   ! Initialize vector
   IF (myid == master) THEN
      ! Restart loic
      IF (lrestart) THEN
         WRITE(6,'(A)') 'Restarting from saved file'
         READ(iRESTART,*) i
         IF (i /= n) THEN
            WRITE(6,'(A)') '***SA Restart Error***'
            WRITE(6,'(A)') ' Free parameter (n) mismatch!'
            STOP
         END IF
         DO i = 1,NP
            READ(iRESTART,*) j,cost_old(i),(x_old(j,i),j=1,n)
         END DO
         REWIND(UNIT=iRESTART)
         ibest     = MINLOC(cost_old,DIM = 1)  
         x = x_old(:,ibest)
         fvecbest = cost_old(ibest)
      ELSE
         x_old(:,1) = x(:)
         DO i = 2, NP
            DO j = 1, n 
               CALL random_number(rand_C1)
               x_old(j,i)=XCmin(j)+rand_C1*(XCmax(j)-XCmin(j))
            END DO
         END DO
      END IF


      ! Do the initial run so we know everything works
      iflag = -1
      CALL fcn(m, n, x, fvec_old(:,1), iflag, iter)
      iflag = GADE_CLEANUP
      CALL fcn(m, n, x, fvec_old(:,1), iflag, iter)
      iter = 0

      ! Compute the normal
      cost_old(1) = SUM(fvec_old(:,1)*fvec_old(:,1))
      ibest = 1
      fvecbest = cost_old(ibest)

      ! Write output
      WRITE(6,'(/,A,/)')  ' Beginning Simulated Annealing'
      WRITE(6,'(A,I4)')   ' Number of Processors: ',numprocs
      WRITE(6,'(A,I4)')   '      Population Size: ',NP
      WRITE(6,'(A,F9.4)') '  Initial Temperature: ',F_XC
      WRITE(6,'(A,F9.4)') '         Cooling Rate: ',CR_XC
      WRITE(6,"(70('='),/,2x,A,3x,A,7x,A)") 'Member ID', 'Processor','Chi-Sq'
      WRITE(6, '(2x,i6,8x,i3,7x,1es12.4)') 0, myid, cost_old(1)
   END IF
      
   ! --- MPI Broadcast ---
#if defined(MPI_OPT)
   CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)
   IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
   CALL MPI_BCAST(x_old,NP*n,MPI_REAL8,master,MPI_COMM_STEL, ierr_mpi)
   IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
   CALL MPI_BCAST(fvec_old,NP*m,MPI_REAL8,master,MPI_COMM_STEL, ierr_mpi)
   IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
   CALL MPI_BCAST(cost_old,NP,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
   IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
#endif

   ! --- Initial Evaluation ---
   ! If not restarting, evaluate the starting states
   IF (.NOT. lrestart) THEN
      i=1
      CALL eval_x_queued(fcn, m, n, NP, x_old, fvec_old, i, MPI_COMM_STEL)
      iter = iter + 1
      IF (myid == master) THEN
         cost_old = SUM(fvec_old*fvec_old,DIM=1)
         ibest    = MINLOC(cost_old,DIM = 1)
         x        = x_old(:,ibest)
         CALL safe_open(iunitx,ierr,'xvec.dat','unknown','formatted',ACCESS_IN='APPEND')
         IF (ierr .ne. 0) STOP 'SA_Evolve Error OPEN(xvec.dat)'
         DO i = 1, NP
            WRITE(iunitx,'(2(2X,I5.5))') n,1
            WRITE(iunitx,'(10ES22.12E3)') x_old(1:n,i)
            WRITE(iunitx,'(ES22.12E3)') cost_old(i)
         END DO
         CLOSE(iunitx)
         ! Check for new minimum
         IF (ibest /= 1) THEN
            x        = x_new(:,ibest)
            fvec     = fvec_old(:,ibest)
            fvecbest = cost_old(ibest)
            iflag = 0
            CALL fcn(m, n, x, fvec, iflag, iter)
            iflag = GADE_CLEANUP
            CALL fcn(m, n, x, fvec, iflag, iter)
            WRITE(6,*) ' '
            WRITE(6,*) '  New Minimum at ',ibest,fvecbest
            WRITE(6,*) ' '
         END IF
      END IF
   END IF

   ! --- MPI Broadcast ---
#if defined(MPI_OPT)
   CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)
   IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
   CALL MPI_BCAST(x_old,NP*n,MPI_REAL8,master,MPI_COMM_STEL, ierr_mpi)
   IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
   CALL MPI_BCAST(fvec_old,NP*m,MPI_REAL8,master,MPI_COMM_STEL, ierr_mpi)
   IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
   CALL MPI_BCAST(cost_old,NP,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
   IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
#endif

   ! --- Main Annealing Loop ---
   DO WHILE (iter < maxfev)
        
      ! 0. Progress Logging
      IF (myid == master) THEN
         WRITE(6,*) ' '
         WRITE(6,*) '  Generation ',iter
         WRITE(6,*) '========== Temperature: ',T,'  Minval: ',MINVAL(cost_old)
         WRITE(6,*) ' '
      END IF
        
      ! 1. Generate new candidate states for all Markov chains
      ! The perturbation amplitude decays as the temperature cools
      step_amp = MAX(0.001d0, T / F_XC) 

      IF (myid==master) THEN
         DO i = 1, NP
            DO j = 1, n
               CALL RANDOM_NUMBER(rand_val)

               ! Random step in [-1, 1], scaled by parameter bounds and step_amp
               ! Base max step size is set to 10% of the domain bounds
               x_new(j, i) = x_old(j, i) + (2.0d0 * rand_val - 1.0d0) * &
                             (XCmax(j) - XCmin(j)) * 0.1d0 * step_amp
                   
               ! Enforce rigid boundary constraints
               IF (x_new(j, i) < XCmin(j)) x_new(j, i) = XCmin(j)
               IF (x_new(j, i) > XCmax(j)) x_new(j, i) = XCmax(j)
            END DO
         END DO
      ENDIF

   ! --- MPI Broadcast ---
#if defined(MPI_OPT)
      CALL MPI_BCAST(x_new,NP*n,MPI_REAL8,master,MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
#endif

      ! 2. Batch Evaluation using STELLOPT queued evaluator
      i=1
      CALL eval_x_queued(fcn, m, n, NP, x_new, fvec_new, i, MPI_COMM_STEL)
      iter = iter + 1

      IF (myid == master) THEN

         ! Compute the new cost
         cost_new = SUM(fvec_new*fvec_new,DIM=1)

         ! Save xvec data
         CALL safe_open(iunitx,ierr,'xvec.dat','unknown','formatted',ACCESS_IN='APPEND')
         IF (ierr .ne. 0) STOP 'SA_Evolve Error OPEN(xvec.dat)'
         DO i = 1, NP
            WRITE(iunitx,'(2(2X,I5.5))') n,1
            WRITE(iunitx,'(10ES22.12E3)') x_new(1:n,i)
            WRITE(iunitx,'(ES22.12E3)') cost_new(i)
         END DO
         CLOSE(iunitx)

         ! Check for new minimum and save output
         IF (MINVAL(cost_new) < fvecbest) THEN
            ibest    = MINLOC(cost_new,DIM = 1)
            x        = x_new(:,ibest)
            fvec     = fvec_new(:,ibest)
            fvecbest = cost_new(ibest)
            iflag = 0
            CALL fcn(m, n, x, fvec, iflag, iter)
            iflag = GADE_CLEANUP
            CALL fcn(m, n, x, fvec, iflag, iter)
            WRITE(6,*) ' '
            WRITE(6,*) '  New Minimum at ',ibest,fvecbest
            WRITE(6,*) ' '
         END IF

         ! 3. Metropolis-Hastings Acceptance Criterion
         DO i = 1, NP
            delta_cost = cost_new(i) - cost_old(i)

            IF (delta_cost < 0.0d0) THEN
               ! Downhill move: Always accept improvements
               x_old(:, i) = x_new(:, i)
               fvec_old(:, i) = fvec_new(:, i)
               cost_old(i) = cost_new(i)
            ELSE
               ! Uphill move: Accept worse solutions probabilistically
               CALL RANDOM_NUMBER(rand_val)
                
               ! Gatecheck to prevent floating point underflow in EXP
               IF ((-delta_cost / T) > -50.0d0) THEN
                  IF (EXP(-delta_cost / T) > rand_val) THEN
                     x_old(:, i) = x_new(:, i)
                     fvec_old(:, i) = fvec_new(:, i)
                     cost_old(i) = cost_new(i)
                  END IF
               END IF
            END IF
         END DO
      END IF

   ! --- MPI Broadcast ---
#if defined(MPI_OPT)
      CALL MPI_BCAST(x_old,n*NP,MPI_REAL8,master,MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(fvec_old,m*NP,MPI_REAL8,master,MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(cost_old,NP,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
#endif

      ! 4. Update Temperature
      ! Geometric cooling schedule applied here
      T = T * CR_XC

    END DO

    ! --- Cleanup ---
    DEALLOCATE(x_new, fvec_new, cost_old, cost_new)

   ! --- MPI Broadcast ---
#if defined(MPI_OPT)
      CALL MPI_BCAST(x,n,MPI_REAL8,master,MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(fvec,m,MPI_REAL8,master,MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
#endif

END SUBROUTINE SA_Evolve