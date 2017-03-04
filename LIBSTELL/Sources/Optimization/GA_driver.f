      SUBROUTINE GA_driver(fcn, n_opt, n_var, x, fvec, tol, eps,
     1     num_iter_opt, max_processors, filename, info, lwa, lrestart )
      USE ga_mod
      USE system_mod
      USE safe_open_mod
      USE mpi_params, ONLY: master, myid
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      include 'mpif.h'                                       !mpi stuff
      INTEGER :: ierr
!DEC$ ENDIF
      INTEGER :: n_opt, n_var, info, lwa, num_iter_opt, max_processors
      REAL(rprec), DIMENSION(n_opt), TARGET :: fvec
      REAL(rprec), DIMENSION(n_var) :: x
      REAL(rprec), DIMENSION(n_var) :: partemp
      REAL(rprec) :: tol, eps, chi_sq, tmp                   !ga_evaluate
      EXTERNAL fcn
      CHARACTER(LEN=*) :: filename
      LOGICAL :: lrestart

      INTEGER :: num_iter_max
      INTEGER :: i, iflag, nfev

c ******************************************************************
c  entries for the 'ga' NAMELIST
c
c npopsiz    -  population SIZE
c idum       -  IF < 0, THEN |idum| is used as seed for random-number gen.
c pmutate    -  probability for random jump mutation
c pcross     -  crossover probability
c ielite     -  /=0  make sure best parent is preserved into decendent populations
c icreep     -  creep mutation flag:  ONLY DO creep mutations IF .ne. 0
c pcreep     -  probability for random creep mutation
c iunifrm    -  =0 single point crossover at random chromosome point
c              /=0 uniform crossover
c iniche     - /=0 turn on niching
c nichflg    - array of flags for the free-parameters,
c              each non-zero ENTRY enables niching for that free-parameter
c iskip
c iend
c nchild     - DEFAULT=1; IF =2, THEN each crossover creates 2 children,
c              the second child having the second parents genes
c parmin     - array specifying minimum value for each free-parameter,
c parmax     - array specifying maximum value for each free-parameter,
c ibound     - =1 then interpret parmin and parmax as scale-factors to be
c              multiplied by the initial guess values for each parameter
c nposibl
c nowrite    - =0 then write output during optimization
c microga    - =0 perform random mutations
c             /=0 perform micro-GA
c unique_ind
c itourny
c
c     IMPORTANT: MPI_PARAMS MODULE MUST BE LOADED BY EXTERNAL CALLS
c                TO MPI_COMM_RANK PRIOR TO THIS SUBROUTINE CALL
c ******************************************************************
      info = 0
      num_iter_max = num_iter_opt

      itourny=1
      MAXgen=ngen
      kountmx=maxgen
      nparam=n_var
      num_obj = n_opt

      IF( ibound .eq. 1 ) THEN
         par_max(:n_var) = x(:n_var)*parmax(:n_var)
         par_min(:n_var) = x(:n_var)*parmin(:n_var)

         WHERE (par_max(:n_var) < par_min(:n_var) )
            partemp(:n_var) = par_max(:n_var)
            par_max(:n_var) = par_min(:n_var)
            par_min(:n_var) = partemp(:n_var)
         END WHERE

      ELSE
         par_max = parmax
         par_min = parmin
      END IF

      IF (ALL(nposibl .eq. 0 )) nposibl=15

      nfev=1
      f_obj => fvec

!     IF (myid .eq. master) WRITE(6, nml = ga_de)

      irestrt = 0
      IF( lrestart ) irestrt = 1
c      IF(irestrt .ne. 0) THEN
c        temp="cp ../ga_restart." // TRIM(filename) // " ."
c        CALL system(temp)
c      END IF

c
c      store initial PARAMETER values as a unique individual
       parent = 0
       child = 0
       iparent = 0
       ichild = 0
       IF(unique_ind .gt. 0 ) THEN
         unique_ind=MIN(unique_ind, npopsiz)
         parent(1:nparam,unique_ind) = x(1:nparam)
       END IF


       CALL ga_sp(fcn, n_opt, fvec, chi_sq, filename, nfev, iflag,
     >            max_processors, myid)

!       iflag=1
!       chi_sq = ga_evaluate(fcn, n_opt, fvec, nparam, parent(1,jbest),
!     >           iflag, nfev)

       IF (myid .eq. master) THEN
          WRITE(6,*) "final solution: "
          WRITE(6,*) "best individual : ", jbest
          WRITE(6,*) "x ", (parent(i,jbest),i=1,nparam)
c         WRITE(6,*) "fvec ",(fvec(i),i=1,n_opt)
          WRITE(6,*) "y ", chi_sq
       END IF

       x(1:n_var) = parent(1:n_var,jbest)

       iflag=-100
       CALL fcn(n_opt, npopsiz, parent(1,jbest), fvec, iflag, nfev)

       END SUBROUTINE GA_driver
