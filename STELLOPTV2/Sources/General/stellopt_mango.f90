! To do:
! - optimum namelist should have centered_differences option.
! - Set finite difference step size.
! - Set max function evaluations.

SUBROUTINE stellopt_optimize_mango(used_mango_algorithm)
        ! This function returns true or false corresponding to whether opt_type is a valid MANGO algorithm.
        ! Or, if STELLOPT was not built with MANGO, this function returns false.
 
!DEC$ IF DEFINED (MANGO)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE mango
      USE stellopt_runtime, ONLY: opt_type
      USE stellopt_vars, ONLY: mango_problem_instance
!-----------------------------------------------------------------------
!     Local Variables
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: used_mango_algorithm
      DOUBLE PRECISION :: best_objective_function
      LOGICAL :: ltst
      CHARACTER(256) :: tstr1, tstr2
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      print *,"Hello world from stellopt_optimize_mango"

      IF (.not. mango_does_algorithm_exist(TRIM(opt_type))) THEN
         ! Then opt_type is not a MANGO algorithm
         used_mango_algorithm = .false.
         RETURN
      END IF

      print *,"BBB"
      ! If we make it here, then opt_type must be a MANGO algorithm.
      used_mango_algorithm = .true.

      ! The first few steps must be done by all processes, not only those in MPI_COMM_STEL = mpi_comm_group_leaders.
      ltst  = .false.
      tstr1 = 'mango_init'
      tstr2 = ''
      CALL stellopt_paraexe(tstr1,tstr2,ltst)
      print *,"CCC"

      ! Now come the steps that only need to be done by MPI_COMM_STEL = mpi_comm_group_leaders.
      CALL mango_set_algorithm_from_string(mango_problem_instance, TRIM(opt_type))
      CALL mango_set_output_filename(mango_problem_instance, "mango_out")
      best_objective_function = mango_optimize(mango_problem_instance)

      ! The last few steps must be done by all processes, not only those in MPI_COMM_STEL = mpi_comm_group_leaders.
      ltst  = .false.
      tstr1 = 'mango_finalize'
      tstr2 = ''
      CALL stellopt_paraexe(tstr1,tstr2,ltst)

!DEC$ ELSE
    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: used_mango_algorithm = .false.
!DEC$ ENDIF
  END SUBROUTINE stellopt_optimize_mango


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE stellopt_mango_init()
        ! This subroutine is run via stellopt_paraexe, so it is run by all procs in MPI_COMM_WORLD.

!DEC$ IF DEFINED (MANGO)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE mango
      USE stellopt_runtime, ONLY: nvars, mtargets, vars, targets, sigmas
      USE stellopt_vars, ONLY: mango_problem_instance, best_residual_function
      USE mpi_inc
      USE mpi_params, ONLY: MPI_COMM_MYWORLD, MPI_COMM_STEL
!-----------------------------------------------------------------------
!     Local Variables
!----------------------------------------------------------------------
      IMPLICIT NONE
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      print *,"Hello from stellopt_mango_init"
      print *,"allocated(vars):",allocated(vars)
      print *,"allocated(targets):",allocated(targets)
      print *,"allocated(sigmas):",allocated(sigmas)
      ALLOCATE(best_residual_function(mtargets))
      !print *,"stellopt_mango_init called on proc ",mango_get_mpi_rank_world(mango_problem_instance)
      !if (mango_is_proc0_world(mango_problem_instance)) then
      if (.true.) THEN
         print *,"nvars:",nvars
         print *,"vars:",vars
         print *,"mtargets:",mtargets
         print *,"targets:",targets
         print *,"sigmas:",sigmas
      end if
      CALL mango_problem_create_least_squares(mango_problem_instance, nvars, vars, mtargets, targets, sigmas, best_residual_function, mango_residual_function)
      CALL mango_set_custom_mpi_communicators(mango_problem_instance, MPI_COMM_WORLD, MPI_COMM_STEL, MPI_COMM_MYWORLD)

CONTAINS
  ! We make mango_residual_function contained by stellopt_mango_init so the call to mango_problem_create_least_squares knows the interface to mango_residual_function.

    SUBROUTINE mango_residual_function(N_parameters, x, N_terms, f, failed, problem)
      USE iso_c_binding
      IMPLICIT NONE
      INTEGER(C_int), INTENT(IN) :: N_parameters
      REAL(C_double), INTENT(IN) :: x(N_parameters)
      INTEGER(C_int), INTENT(IN) :: N_terms
      REAL(C_double), INTENT(OUT) :: f(N_terms)
      INTEGER(C_int), INTENT(OUT) :: failed
      TYPE(mango_problem), value, INTENT(IN) :: problem

      INTEGER :: iflag

      print *,"Hello from mango_residual_function on proc",mango_get_mpi_rank_world(problem)

      iflag = mango_get_mpi_rank_group_leaders(problem)
      CALL stellopt_fcn(N_terms, N_parameters, x, f, iflag, mango_get_function_evaluations(problem))
      failed = 0

    END SUBROUTINE mango_residual_function

!DEC$ ENDIF
  END SUBROUTINE stellopt_mango_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE stellopt_mango_finalize()
        ! This subroutine is run via stellopt_paraexe, so it is run by all procs in MPI_COMM_WORLD.

!DEC$ IF DEFINED (MANGO)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE mango
      USE stellopt_vars, ONLY: mango_problem_instance, best_residual_function
!-----------------------------------------------------------------------
!     Local Variables
!----------------------------------------------------------------------
      IMPLICIT NONE
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      print *,"stellopt_mango_finalize called on proc ",mango_get_mpi_rank_world(mango_problem_instance)
      DEALLOCATE(best_residual_function)
      CALL mango_problem_destroy(mango_problem_instance)

!DEC$ ENDIF
    END SUBROUTINE stellopt_mango_finalize
