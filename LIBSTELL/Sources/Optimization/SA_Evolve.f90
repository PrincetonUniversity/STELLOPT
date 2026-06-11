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
    IMPLICIT NONE

    ! --- Arguments ---
    EXTERNAL :: fcn
    INTEGER, INTENT(IN) :: m, n, NP
    INTEGER, INTENT(INOUT) :: maxfev
    REAL(8), INTENT(IN) :: XCmin(n), XCmax(n)
    REAL(8), INTENT(INOUT) :: x(n, NP), fvec(m, NP)
    REAL(8), INTENT(IN) :: F_XC         ! SA Use: Initial Temperature (T0)
    REAL(8), INTENT(IN) :: CR_XC        ! SA Use: Cooling Rate (e.g., 0.95)
    INTEGER, INTENT(IN) :: strategy     ! SA Use: Perturbation strategy
    INTEGER, INTENT(IN) :: CR_strategy  ! SA Use: Cooling schedule strategy
    INTEGER, INTENT(IN) :: iWRITE
    INTEGER, INTENT(IN) :: iRESTART
    LOGICAL, INTENT(IN) :: lrestart

    ! --- Local Variables ---
    INTEGER :: i, j, fev
    REAL(8) :: T
    REAL(8), ALLOCATABLE :: x_new(:,:), fvec_new(:,:)
    REAL(8), ALLOCATABLE :: cost_old(:), cost_new(:)
    REAL(8) :: rand_val, delta_cost, step_amp
    
    ! --- Interface for the STELLOPT queued evaluator ---
    ! (Adjust the argument order/types to match your exact LIBSTELL branch)
    INTERFACE
        SUBROUTINE eval_x_queued(fcn, m, n, NP, x_new, fvec_new)
            IMPLICIT NONE
            EXTERNAL :: fcn
            INTEGER, INTENT(IN) :: m, n, NP
            REAL(8), INTENT(IN) :: x_new(n, NP)
            REAL(8), INTENT(OUT) :: fvec_new(m, NP)
        END SUBROUTINE eval_x_queued
    END INTERFACE

    ! --- Memory Allocation ---
    ALLOCATE(x_new(n, NP))
    ALLOCATE(fvec_new(m, NP))
    ALLOCATE(cost_old(NP))
    ALLOCATE(cost_new(NP))

    fev = 0
    T = F_XC  ! Set initial temperature

    ! --- Initial Evaluation ---
    ! If not restarting, evaluate the starting states
    IF (.NOT. lrestart) THEN
        CALL eval_x_queued(fcn, m, n, NP, x, fvec)
        fev = fev + NP
    END IF

    ! Compute initial costs. Assuming a least-squares formulation 
    ! standard to STELLOPT (Sum of squared residuals).
    DO i = 1, NP
        cost_old(i) = SUM(fvec(:, i)**2)
    END DO

    ! --- Main Annealing Loop ---
    DO WHILE (fev < maxfev)
        
        ! 1. Generate new candidate states for all Markov chains
        ! The perturbation amplitude decays as the temperature cools
        step_amp = MAX(0.001d0, T / F_XC) 

        DO i = 1, NP
            DO j = 1, n
                CALL RANDOM_NUMBER(rand_val)
                
                ! Random step in [-1, 1], scaled by parameter bounds and step_amp
                ! Base max step size is set to 10% of the domain bounds
                x_new(j, i) = x(j, i) + (2.0d0 * rand_val - 1.0d0) * &
                              (XCmax(j) - XCmin(j)) * 0.1d0 * step_amp
                
                ! Enforce rigid boundary constraints
                IF (x_new(j, i) < XCmin(j)) x_new(j, i) = XCmin(j)
                IF (x_new(j, i) > XCmax(j)) x_new(j, i) = XCmax(j)
            END DO
        END DO

        ! 2. Batch Evaluation using STELLOPT queued evaluator
        CALL eval_x_queued(fcn, m, n, NP, x_new, fvec_new)
        fev = fev + NP

        ! 3. Metropolis-Hastings Acceptance Criterion
        DO i = 1, NP
            cost_new(i) = SUM(fvec_new(:, i)**2)
            delta_cost = cost_new(i) - cost_old(i)

            IF (delta_cost < 0.0d0) THEN
                ! Downhill move: Always accept improvements
                x(:, i) = x_new(:, i)
                fvec(:, i) = fvec_new(:, i)
                cost_old(i) = cost_new(i)
            ELSE
                ! Uphill move: Accept worse solutions probabilistically
                CALL RANDOM_NUMBER(rand_val)
                
                ! Gatecheck to prevent floating point underflow in EXP
                IF ((-delta_cost / T) > -50.0d0) THEN
                    IF (EXP(-delta_cost / T) > rand_val) THEN
                        x(:, i) = x_new(:, i)
                        fvec(:, i) = fvec_new(:, i)
                        cost_old(i) = cost_new(i)
                    END IF
                END IF
            END IF
        END DO

        ! 4. Update Temperature
        ! Geometric cooling schedule applied here
        T = T * CR_XC
        
        ! 5. Progress Logging
        IF (iWRITE > 0 .AND. MOD(fev, NP*10) == 0) THEN
            WRITE(iWRITE, '(A, I8, A, ES12.4, A, ES12.4)') &
            & ' SA FEV: ', fev, ' | Temp: ', T, ' | Best Cost: ', MINVAL(cost_old)
        END IF

    END DO

    ! --- Cleanup ---
    DEALLOCATE(x_new, fvec_new, cost_old, cost_new)

END SUBROUTINE SA_Evolve