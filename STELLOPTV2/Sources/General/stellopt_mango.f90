SUBROUTINE stellopt_optimize_mango(used_mango_algorithm)
        ! This subroutine sets its argument (used_mango_algorithm) to true or false corresponding to whether opt_type is a valid MANGO algorithm.
        ! Or, if STELLOPT was not built with MANGO, this function returns false.
 
!DEC$ IF DEFINED (MANGO)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE mango_mod
      USE stellopt_runtime, ONLY: opt_type, targets, sigmas, epsfcn, lcentered_differences, vars_min, vars_max
      USE stellopt_vars, ONLY: mango_problem_instance, nfunc_max, mango_bound_constraints
      USE mpi_params, ONLY: myid
!-----------------------------------------------------------------------
!     Local Variables
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: used_mango_algorithm
      DOUBLE PRECISION :: best_objective_function
      LOGICAL :: ltst
      CHARACTER(256) :: tstr1, tstr2
      INTEGER :: m, ncnt, iflag
      DOUBLE PRECISION :: fvec_temp(1)
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      print *,"Hello world from stellopt_optimize_mango"

      IF (.not. mango_does_algorithm_exist(TRIM(opt_type))) THEN
         ! Then opt_type is not a MANGO algorithm
         used_mango_algorithm = .false.
         RETURN
      END IF

      ! If we make it here, then opt_type must be a MANGO algorithm.
      used_mango_algorithm = .true.

      ! At this point in stellopt's execution, the 'targets' and 'sigmas' arrays have been allocated
      ! but not populated with their values. However, mango_optimize() needs to know the values of 'targets' and 'sigmas'
      ! before it starts evaluating the residual function. So, let's populate 'targets' and 'sigmas' now.
      ! This is done by calling stellopt_load_targets().
      m = -1
      ncnt = 0
      iflag = 0
      CALL stellopt_load_targets(m,fvec_temp,iflag,ncnt)

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
            PRINT *,"Adding bound constraints in mango:"
            print *,"vars_min:",vars_min
            print *,"vars_max:",vars_max
            CALL mango_set_bound_constraints(mango_problem_instance, vars_min, vars_max)
         END IF
      ELSE
         IF (myid==0) PRINT *,"NOT including bound constraints in mango."
      END IF

      best_objective_function = mango_optimize(mango_problem_instance)

      ! The last few steps must be done by all processes, not only those in MPI_COMM_STEL = mpi_comm_group_leaders.
      ltst  = .false.
      tstr1 = 'mango_finalize'
      tstr2 = ''
      CALL stellopt_paraexe(tstr1,tstr2,ltst)

!DEC$ ELSE
    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: used_mango_algorithm

    used_mango_algorithm = .false.
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
      USE stellopt_runtime, ONLY: nvars, mtargets, vars, targets, sigmas
      USE stellopt_vars, ONLY: mango_problem_instance, best_residual_function
      USE mpi_inc
      USE mpi_params, ONLY: MPI_COMM_MYWORLD, MPI_COMM_STEL
!-----------------------------------------------------------------------
!     Local Variables
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: nvars_allprocs, mtargets_allprocs, ierr
      !EXTERNAL mango_residual_function
      procedure(vector_function_interface) :: mango_residual_function
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
      !print *,"stellopt_mango_init called on proc ",mango_get_mpi_rank_world(mango_problem_instance)
      !if (mango_is_proc0_world(mango_problem_instance)) then
      if (.true.) THEN
         print "(a,i4,a,i6)","Hello from stellopt_mango_init. nvars:",nvars," mtargets:",mtargets
         !print *,"vars:",vars
         !print *,"mtargets:",mtargets
         !print *,"targets:",targets
         !print *,"sigmas:",sigmas
      end if
      CALL mango_problem_create_least_squares(mango_problem_instance, nvars_allprocs, vars, mtargets_allprocs, targets, sigmas, best_residual_function, mango_residual_function)
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
      USE stellopt_runtime, ONLY: targets, sigmas
      USE mpi_params, ONLY: myworkid, myid

      IMPLICIT NONE

      INTEGER(C_int), INTENT(IN) :: N_parameters
      REAL(C_double), INTENT(IN) :: x(N_parameters)
      INTEGER(C_int), INTENT(IN) :: N_terms
      REAL(C_double), INTENT(OUT) :: f(N_terms)
      INTEGER(C_int), INTENT(OUT) :: failed
      TYPE(mango_problem), VALUE, INTENT(IN) :: problem
      TYPE(C_ptr), VALUE, INTENT(IN) :: user_data

      INTEGER :: iflag
      INTEGER, SAVE :: N_function_evaluations = -1

      N_function_evaluations = N_function_evaluations + 1
      !print *,"Hello from mango_residual_function on proc",mango_get_mpi_rank_world(problem),", function_evaluations=",mango_get_function_evaluations(problem)
      print "(a,i4,a,i7)","Hello from mango_residual_function on proc",myid,", function_evaluations=",N_function_evaluations

      ! iflag should be the processor number, except that iflag=-1 has the effect of turning 
      ! on output from VMEC and other codes, as traditionally is done for the first function evaluation.
      !iflag = mango_get_mpi_rank_group_leaders(problem)
      iflag = myid
      !if (mango_get_function_evaluations(problem) <= 1 .and. mango_get_proc0_world(problem)) iflag = -1
      if (N_function_evaluations < 1 .and. myid==0) iflag = -1

      !print "(4(a,i3),a,256(es24.15))","Proc",mango_get_mpi_rank_world(problem)," N_terms:",N_terms," N_parameters:",N_parameters," iflag:",iflag," x:",x
      !print *,"size(f):",size(f)," f:",f
      !CALL stellopt_fcn(N_terms, N_parameters, x, f, iflag, mango_get_function_evaluations(problem))
      CALL stellopt_fcn(N_terms, N_parameters, x, f, iflag, N_function_evaluations)

!      print "(a,i3,a,i8,4(a,1(es24.15)))","Proc",mango_get_mpi_rank_world(problem)," iflag:",iflag," f from stellopt_fcn:",f," targets:",targets," sigmas:",sigmas," f for mango:",sigmas*f+targets

      ! Stellopt's convention is that fvec is (values - targets) / sigmas. This can be seen from the line
      ! fvec(1:m) = (vals(1:m)-targets(1:m))/ABS(sigmas(1:m))
      ! at the end of stellopt_load_targets().
      ! However mango's f vector should be the original values, not shifted and scaled by "targets" and "sigmas".
      ! So we need to convert stellopt's f to mango's f here:
      f = sigmas * f + targets

      failed = 0
      if (iflag < 0) failed = 1

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
