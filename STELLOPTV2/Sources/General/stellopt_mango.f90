! This file contains procedures for using the MANGO interface to optimization libraries,
! https://github.com/hiddenSymmetries/mango

! MANGO's fortran interface assumes DOUBLE PRECISION type, so DOUBLE PRECISION is used in many places here instead of REAL(rprec).


!DEC$ IF DEFINED (MANGO)
MODULE stellopt_mango_mod
  ! This module defines a few variables used by MANGO that must persist beyond and between calls of the other subroutines in this file.

  USE mango_mod

  IMPLICIT NONE

  TYPE(mango_problem) :: mango_problem_instance
  DOUBLE PRECISION, ALLOCATABLE :: best_residual_function(:), mango_targets(:), mango_sigmas(:)

END MODULE stellopt_mango_mod
!DEC$ ENDIF

SUBROUTINE stellopt_optimize_mango(used_mango_algorithm, N_function_evaluations)
        ! This subroutine sets its argument (used_mango_algorithm) to true or false corresponding to whether opt_type is a valid MANGO algorithm.
        ! Or, if STELLOPT was not built with MANGO, this function returns false.
 
!DEC$ IF DEFINED (MANGO)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE mango_mod
      USE stellopt_mango_mod
      USE stellopt_runtime, ONLY: opt_type, epsfcn, lcentered_differences, vars_min, vars_max
      USE stellopt_vars, ONLY: nfunc_max, mango_bound_constraints
      USE mpi_params, ONLY: myid
!-----------------------------------------------------------------------
!     Local Variables
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: used_mango_algorithm
      INTEGER, INTENT(OUT) :: N_function_evaluations
      DOUBLE PRECISION :: best_objective_function
      LOGICAL :: ltst
      CHARACTER(256) :: tstr1, tstr2
      INTEGER :: m, ncnt, iflag
      DOUBLE PRECISION :: fvec_temp(1)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      IF (.not. mango_does_algorithm_exist(TRIM(opt_type))) THEN
         ! Then opt_type is not a MANGO algorithm
         used_mango_algorithm = .false.
         RETURN
      END IF

      ! If we make it here, then opt_type must be a MANGO algorithm.
      used_mango_algorithm = .true.

      ! The first few steps must be done by all processes, not only those in MPI_COMM_STEL = mpi_comm_group_leaders.
      ltst  = .false.
      tstr1 = 'mango_init'
      tstr2 = ''
      CALL stellopt_paraexe(tstr1,tstr2,ltst)

      ! Now come the steps that only need to be done by MPI_COMM_STEL = mpi_comm_group_leaders.
      CALL mango_set_algorithm_from_string(mango_problem_instance, TRIM(opt_type))
      CALL mango_set_output_filename(mango_problem_instance, "mango_out")
      CALL mango_set_finite_difference_step_size(mango_problem_instance, epsfcn)
      CALL mango_set_centered_differences(mango_problem_instance, lcentered_differences)
      CALL mango_set_max_function_evaluations(mango_problem_instance, nfunc_max)
      CALL mango_set_print_residuals_in_output_file(mango_problem_instance, .false.) ! This line makes the mango_out file much smaller.
      IF (mango_bound_constraints) THEN
         IF (myid==0) THEN
            CALL mango_set_bound_constraints(mango_problem_instance, vars_min, vars_max)
         END IF
      END IF

      best_objective_function = mango_optimize(mango_problem_instance)

      N_function_evaluations = mango_get_function_evaluations(mango_problem_instance)

      ! The last few steps must be done by all processes, not only those in MPI_COMM_STEL = mpi_comm_group_leaders.
      ltst  = .false.
      tstr1 = 'mango_finalize'
      tstr2 = ''
      CALL stellopt_paraexe(tstr1,tstr2,ltst)

!DEC$ ELSE
    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: used_mango_algorithm
    INTEGER, INTENT(OUT) :: N_function_evaluations

    used_mango_algorithm = .false.
    N_function_evaluations = 0
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
      USE mango_mod
      USE stellopt_mango_mod
      USE stellopt_runtime, ONLY: nvars, mtargets, vars
      USE mpi_inc
      USE mpi_params, ONLY: MPI_COMM_MYWORLD, MPI_COMM_STEL
!-----------------------------------------------------------------------
!     Local Variables
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: nvars_allprocs, mtargets_allprocs, ierr
      !EXTERNAL mango_residual_function
      PROCEDURE(vector_function_interface) :: mango_residual_function
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ! In STELLOPT, nvars and mtargets are 0 on procs other than the group leaders.
      ! MANGO may require values for the corresponding terms N_parameters and N_terms that are >0;
      ! the policy on this in MANGO is not yet finalized. So let's broadcast these variables to
      ! make sure they are >0 on all procs.
      nvars_allprocs = nvars
      mtargets_allprocs = mtargets
      CALL MPI_BCAST(nvars_allprocs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(mtargets_allprocs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      ALLOCATE(best_residual_function(mtargets_allprocs))
      ALLOCATE(mango_targets(mtargets_allprocs))
      ALLOCATE(mango_sigmas(mtargets_allprocs))
      ! Shifting and scaling of the terms in the objective function is all handled by stellopt. MANGO's capability to shift and scale the residuals will not be used. Hence we set targets=0 and sigmas=1:
      mango_targets = 0
      mango_sigmas = 1
      CALL mango_problem_create_least_squares(mango_problem_instance, nvars_allprocs, vars, mtargets_allprocs, mango_targets, mango_sigmas, best_residual_function, mango_residual_function)
      CALL mango_mpi_partition_set_custom(mango_problem_instance, MPI_COMM_WORLD, MPI_COMM_STEL, MPI_COMM_MYWORLD)
      CALL mango_mpi_partition_write(mango_problem_instance, "mango_mpi")

!DEC$ ENDIF
  END SUBROUTINE stellopt_mango_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!DEC$ IF DEFINED (MANGO)
    SUBROUTINE mango_residual_function(N_parameters, x, N_terms, f, failed, problem, user_data) BIND(C)
      ! "BIND(C)" is needed above because this subroutine will be called from C++ in MANGO.
      USE iso_c_binding
      USE mango_mod
      USE stellopt_mango_mod
      USE stellopt_runtime, ONLY: bigno, opt_type
      USE mpi_params, ONLY: myworkid, myid
      USE gade_mod, ONLY: GADE_CLEANUP
      USE fdjac_mod, ONLY: flag_cleanup

      IMPLICIT NONE

      INTEGER(C_int), INTENT(IN) :: N_parameters
      REAL(C_double), INTENT(IN) :: x(N_parameters)
      INTEGER(C_int), INTENT(IN) :: N_terms
      REAL(C_double), INTENT(OUT) :: f(N_terms)
      INTEGER(C_int), INTENT(OUT) :: failed
      TYPE(mango_problem), VALUE, INTENT(IN) :: problem
      TYPE(C_ptr), VALUE, INTENT(IN) :: user_data

      INTEGER :: iflag
      INTEGER :: N_function_evaluations

      N_function_evaluations = mango_get_function_evaluations(problem)
      ! The first time this subroutine is called, N_function_evaluations will be 0 (not 1).

      ! There is a slight problem in that N_function_evaluations, which was returned by mango_get_function_evaluations() above, is the number of completed evaluations before now,
      ! whereas the global evaluation number for the present run in mango_out has not yet been assigned: this will only happen once this subroutine completes and control
      ! is returned to mango itself. So the evaluation numbers in stellopt vs mango may not exactly agree (even accounting for the +1 offset in mango_out).     

      ! iflag should be the processor number, except that iflag=-1 has the effect of turning 
      ! on output from VMEC and other codes, as traditionally is done for the first function evaluation.
      !iflag = mango_get_mpi_rank_group_leaders(problem)
      iflag = myid
      if (N_function_evaluations < 1 .and. myid==0) iflag = -1

      CALL stellopt_fcn(N_terms, N_parameters, x, f, iflag, N_function_evaluations)

      ! If a failure of any code occurs during stellopt_fcn, stellopt_fcn sets "fvec(1:m) = 10*SQRT(bigno/m)"
      IF (maxval(f(1:N_terms)) < 10*SQRT(bigno/N_terms)) THEN
         ! If we get here, stellopt_fcn must have succeeded.

         failed = 0

         ! When stellopt_fcn is called with iflag < -2, stellopt_fcn calls stellopt_clean_up (which writes the stellopt.<extension> file) and then stellopt_fcn immediately returns,
         ! without actually evaluating the objective function again. We do this here in order to write the stellopt.<extension> file.
         ! For now, we only record output from worker group 1, since (for now) mango_get_function_evaluations does not return meaningful results on other worker groups.

         ! It is important to do this stellopt_clean_up step only if the previous stellopt_fcn succeeded. Otherwise stellopt_clean_up will
         ! call stellopt_paraexe("parvmec_write") and this will hang.
         iflag = FLAG_CLEANUP ! All procs except master use this value, which has the effect of doing nothing in stellopt_clean_up.
         IF (myid==0) iflag = GADE_CLEANUP
         CALL stellopt_fcn(N_terms, N_parameters, x, f, iflag, N_function_evaluations)

      ELSE
         ! If we get here, stellopt_fcn must have detected a failure.
         failed = 1
      END IF

      IF (myid==0 .and. N_function_evaluations == 0) THEN
         WRITE (6,"(A)") ""
         WRITE (6,"(A)") "Optimization with MANGO (https://github.com/hiddenSymmetries/mango)"
         WRITE (6,"(A,A)") " optimization algorithm: ",TRIM(opt_type)
         WRITE (6,"(A)") ""
         WRITE (6,"(A)") " worker group,  objective function"
         WRITE (6,"(A)") " ---------------------------------"
         FLUSH (6)
      END IF

      WRITE (6,"(I5,es24.15)") myid, SUM(f*f)

    END SUBROUTINE mango_residual_function

!DEC$ ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE stellopt_mango_finalize()
        ! This subroutine is run via stellopt_paraexe, so it is run by all procs in MPI_COMM_WORLD.

!DEC$ IF DEFINED (MANGO)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE mango_mod
      USE stellopt_mango_mod
!-----------------------------------------------------------------------
!     Local Variables
!----------------------------------------------------------------------
      IMPLICIT NONE
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      DEALLOCATE(best_residual_function, mango_targets, mango_sigmas)
      CALL mango_problem_destroy(mango_problem_instance)

!DEC$ ENDIF
    END SUBROUTINE stellopt_mango_finalize
