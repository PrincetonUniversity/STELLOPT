!-----------------------------------------------------------------------
!     Subroutine:    DE2_Evolve
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          01/17/2014
!     Description:   This subroutine preforms a Differential Evolution
!                    (DE) minimization of a target functional.  This is
!                    a Fortran reimplementation of the differential
!                    evolution (DE) algorithm as implemented by Storn
!                    and Price.  Modifications by Zarnstorff (PPPL) to
!                    implement ALL stragegies and make random
!                    selection process be compatible with c versions
!                    and restart capability.  Multi-processor capability
!                    added by Zarnstorff (PPPL) and Hirshman (ORNL).
!  Refences:
!  Storn, R., and Price, K.V., (1996). Minimizing the REAL FUNCTION of the
!    ICEC''96 contest by differential evolution. IEEE conf. on Evolutionary
!    Computation, 842-844.
!
!=======Choice of strategy=================================================
!  We have tried to come up with a sensible naming-convention: DE/x/y/z
!    DE :  stands for Differential Evolution
!    x  :  a string which denotes the vector to be perturbed
!    y  :  number of difference vectors taken for perturbation of x
!    z  :  crossover method (exp = exponential, bin = binomial)
!
!  There are some simple rules which are worth following:
!  1)  F is usually between 0.5 and 1 (in rare cases > 1)
!  2)  CR is between 0 and 1 with 0., 0.3, 0.7 and 1. being worth to be tried first
!  3)  To start off NP = 10*D is a reasonable choice. Increase NP IF
!      misconvergence happens.
!  4)  If you increase NP, F usually has to be decreased
!  5)  When the DE/best... schemes fail DE/rand... usually works and vice versa
!
!
! Here are Storn''s comments on the different strategies:
!
! (1) DE/best/1/'z'
!     Our oldest strategy but still not bad. However, we have found several
!     optimization problems WHERE misconvergence occurs.
!
! (2) DE/rand/1/'z'
!     This is one of my favourite strategies. It works especially well when the
!     "bestit[]"-schemes EXPerience misconvergence. Try e.g. F=0.7 and CR=0.5
!     as a first guess.
!
! (3) DE/rand-to-best/1/'z'
!     This strategy seems to be one of the best strategies. Try F=0.85 and CR=1.
!     If you get misconvergence try to increase NP. If this does not help you
!     should play around with all three control variables.
!
! (4) DE/best/2/'z' is another powerful strategy worth trying
!
! (5) DE/rand/2/'z' seems to be a robust optimizer for many functions
!
!===========================================================================
!-----------------------------------------------------------------------
      SUBROUTINE DE2_Evolve(fcn, m, n, NP, XCmin, XCmax, x, fvec,
     1                      maxfev,F_XC,CR_XC,strategy,CR_strategy,
     2                      iWRITE,iRESTART,lrestart)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE mpi_params
      USE safe_open_mod
      USE gade_mod, ONLY: gade_cleanup
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      include 'mpif.h'                                       !mpi stuff
!DEC$ ENDIF
!-----------------------------------------------------------------------
!     Input Variables
!        fcn          Target function one wishes to mimizie
!        m            Dimensionality of the function
!        n            Dimensionality of the function space
!        NP           Population size
!        XCmin        Minimum values of function space
!        XCmax        Maximum values of function space
!        x            Vector specifying mimium location in function space
!        fvec         Vector of function values at minimum
!        maxfev       Maximum number of function evaluations
!        F_XC         Mutation scaling factor
!        CR_XC        Crossover factor
!        strategy     Mutation strategy
!        CR_strategy  Cross-Over Strategy
!        iWRITE       Unit number for output
!        iRESTART     Unit number for restart file
!        lrestart     Logical controlling restart
!----------------------------------------------------------------------
      LOGICAL, INTENT(in) :: lrestart
      INTEGER, INTENT(in) :: m, n, NP, maxfev, strategy, CR_strategy,
     1                       iWRITE,iRESTART
      REAL(rprec), INTENT(in)  :: F_XC, CR_XC
      REAL(rprec), INTENT(in)  :: XCmin(n), XCmax(n)
      REAL(rprec), INTENT(inout) :: x(n)
      REAL(rprec), INTENT(inout) :: fvec(m)
      EXTERNAL  fcn
      
!----------------------------------------------------------------------
!     Local Variables
!        i,j            Dummy index
!        i_free         Number of variables in restart array
!        fnorm          Function normalization SQRT(F.F)
!        rand_C1        Random number
!        x_array        Array of all X's for the population
!        x_array        Array of all fvals for the population
!----------------------------------------------------------------------
      LOGICAL, ALLOCATABLE :: lcross(:,:)
      INTEGER :: i,j, iter, iflag, ierr, numsent, iunitx, ibest, i_free
      !INTEGER, DIMENSION(1) :: ibest
      INTEGER, ALLOCATABLE :: a1(:), a2(:), a3(:), a4(:), a5(:)
      REAL(rprec) :: rand_C1, fnorm, fnorm_min
      REAL(rprec), ALLOCATABLE :: fnorm_array(:), fnorm_new(:),
     1                            temp_fvec(:), x_temp(:)
      REAL(rprec), ALLOCATABLE :: x_array(:,:), fval_array(:,:),
     1                            x_new(:,:)
!DEC$ IF DEFINED (MPI_OPT)
      INTEGER :: status(MPI_STATUS_size)                     !mpi stuff
      INTEGER :: sender
!DEC$ ENDIF
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      ! Perform Allocations
      ierr_mpi = 0
      ierr = 0
      iter = 0
      fnorm_min = 1.0E30
      ALLOCATE (fnorm_array(NP),temp_fvec(m),x_temp(n), stat=ierr)
      IF (ierr .ne. 0) STOP 'DE2_Evolve Error ALLOC(1)'
      ALLOCATE (x_array(NP,n),fval_array(NP,m), stat=ierr)
      IF (ierr .ne. 0) STOP 'DE2_Evolve Error ALLOC(2)'
      ALLOCATE (x_new(NP,n), stat=ierr)
      IF (ierr .ne. 0) STOP 'DE2_Evolve Error ALLOC(3)'
      temp_fvec = 0
      
      ! Initialize Starting Positions
      CALL RANDOM_SEED()
      IF (myid == master) THEN
         ALLOCATE (fnorm_new(NP), stat=ierr)
         IF (ierr .ne. 0) STOP 'DE2_Evolve Error ALLOC(4)'
         x_array(1,:) = x(:)
         DO i = 2, NP
            DO j = 1, n
               CALL random_number(rand_C1)
               fnorm=XCmin(j)+rand_C1*(XCmax(j)-XCmin(j))
               x_array(i,j)=fnorm
            END DO
         END DO
         ! Handle restart case
         IF (lrestart) THEN
            WRITE(6,'(A)') 'Restarting from saved file'
            READ(iRESTART,*) i_free
            IF (i_free /= n) THEN
               WRITE(6,'(A)') '***DE Restart Error***'
               WRITE(6,'(A)') ' Free parameter mismatch!'
               STOP
            END IF
            DO i = 1,NP
               READ(iRESTART,*) j,fnorm_array(i),(x_array(i,j),j=1,n)
            END DO
            REWIND(UNIT=iRESTART)
            iter = 0
            ibest     = MINLOC(fnorm_array,DIM = 1)  
            x_temp = x_array(ibest,:)
            CALL fcn(m, n, x_temp, temp_fvec, iflag, iter)
            iflag = GADE_CLEANUP
            CALL fcn(m, n, x_temp, temp_fvec, iflag, iter)
            fnorm = SUM(temp_fvec*temp_fvec)
            fnorm_array(1) = fnorm
            fval_array(1,:) = temp_fvec
            fnorm_min = fnorm
            WRITE (6, 1327) numprocs, NP, strategy, CR_strategy
            WRITE(6, '(2x,i6,8x,i3,7x,1es12.4)') 0, myid, fnorm
         ELSE
            ! Do the initial run so we know everything works
            iflag = -1
            x_temp = x_array(1,:)
            CALL fcn(m, n, x_temp, temp_fvec, iflag, iter)
            iflag = GADE_CLEANUP
            iter = 0
            CALL fcn(m, n, x_temp, temp_fvec, iflag, iter)
            fnorm = SUM(temp_fvec*temp_fvec)
            fnorm_array(1) = fnorm
            fval_array(1,:) = temp_fvec
            fnorm_min = fnorm
            WRITE (6, 1327) numprocs, NP, strategy, CR_strategy
            WRITE(6, '(2x,i6,8x,i3,7x,1es12.4)') 0, myid, fnorm
         END IF
      END IF
      
!DEC$ IF DEFINED (MPI_OPT)
      ! BROADCAST The X_Array
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(x_array,NP*n,MPI_REAL8,master,MPI_COMM_STEL,
     1               ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(temp_fvec,m,MPI_REAL8,master,MPI_COMM_STEL,
     1               ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(fnorm_min,1,MPI_REAL8,master,MPI_COMM_STEL,
     1               ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF
      
      ! Evaluate the Population if necessary
      IF (lrestart) THEN
         IF (myid == master) 
     1        WRITE(6,*) 'Using evaluation from restart file'
         CALL FLUSH(6)
      ELSE
!DEC$ IF DEFINED (MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
         IF (myid == master) THEN
            numsent = 1
            DO j = 1,MIN(numprocs-1,NP-1)
               CALL MPI_SEND(j+1,1,MPI_INTEGER,j,j,MPI_COMM_STEL,
     1                        ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) 
     1                        CALL mpi_stel_abort(ierr_mpi)
               numsent = numsent + 1
            END DO
            DO j = 1, NP-1
               CALL MPI_RECV(temp_fvec,m,MPI_REAL8,MPI_ANY_SOURCE,
     1                       MPI_ANY_TAG,MPI_COMM_STEL,status,
     2                       ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) 
     1                        CALL mpi_stel_abort(ierr_mpi)
               sender = status(MPI_SOURCE)
               i      = status(MPI_TAG)
               fval_array(i,:) = temp_fvec
               fnorm = SUM(temp_fvec*temp_fvec)
               fnorm_array(i) = fnorm
               write (6, '(2x,i6,8x,i3,7x,1es12.4)') i,
     1             sender, fnorm
               CALL FLUSH(6)
               IF (numsent < NP) THEN
                  numsent = numsent + 1
                  CALL MPI_SEND(numsent, 1, MPI_INTEGER,
     1                      sender, numsent, MPI_COMM_STEL, ierr_mpi)
                  IF (ierr_mpi /= MPI_SUCCESS) 
     1                        CALL mpi_stel_abort(ierr_mpi)
               ELSE
                  i = 0
                  CALL MPI_SEND(i, 1, MPI_INTEGER,
     1                      sender, numsent, MPI_COMM_STEL, ierr_mpi)
                  IF (ierr_mpi /= MPI_SUCCESS) 
     1                        CALL mpi_stel_abort(ierr_mpi)
               END IF
            END DO
         ELSE IF (myid <= NP-1) THEN
            DO
               CALL MPI_RECV(j,1,MPI_INTEGER,master,MPI_ANY_TAG,
     1                       MPI_COMM_STEL,status,ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) 
     1                        CALL mpi_stel_abort(ierr_mpi)
               IF (j > 0) THEN
                  x_temp     = x_array(j,:)
                  temp_fvec  = 0
                  iflag = j
                  CALL fcn(m, n, x_temp, temp_fvec, iflag, iter)
                  fnorm = SUM(temp_fvec*temp_fvec)
                  CALL MPI_SEND(temp_fvec,m, MPI_REAL8, master,
     1                       j, MPI_COMM_STEL, ierr_mpi)
                  IF (ierr_mpi /= MPI_SUCCESS) 
     1                        CALL mpi_stel_abort(ierr_mpi)
               ELSE
                  EXIT
               END IF
            END DO
         ELSE
         END IF
         
         CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
         CALL MPI_BCAST(fnorm_array,NP,MPI_REAL8,
     1                  master,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ELSE
         DO i = 1, NP
            x_temp = x_array(i,:)
            iflag = i
            CALL fcn(m, n, x_temp, temp_fvec, iflag, iter)
            fval_array(i,:) = temp_fvec
            fnorm_array(i)  = SUM(temp_fvec*temp_fvec)
         END DO
!DEC$ ENDIF
      END IF
      
      ! Output to the xvec.dat file
      IF (myid == master) THEN
         CALL safe_open(iunitx,ierr,'xvec.dat','unknown',
     1                  'formatted',ACCESS_IN='APPEND')
         IF (ierr .ne. 0) STOP 'DE2_Evolve Error OPEN(xvec.dat)'
         DO i = 1, NP
            WRITE(iunitx,'(2(2X,I5.5))') n,1
            WRITE(iunitx,'(10ES22.12E3)') x_array(i,1:n)
            WRITE(iunitx,'(ES22.12E3)') fnorm_array(i)
         END DO
         CLOSE(iunitx)
      END IF

         
      ! Get the best value
      ibest     = MINLOC(fnorm_array,DIM = 1)   
      fnorm     = fnorm_array(ibest)
      IF (fnorm_min > fnorm) THEN
         fnorm_min = fnorm
         ! Save the output
         IF (myid == master .and. ibest /= 1
     1       .and. .not.lrestart) THEN
            x_temp = x_array(ibest,:)
            iter  = 1
            iflag = 0
            CALL fcn(m, n, x_temp, temp_fvec, iflag, iter)
            iflag = GADE_CLEANUP
            CALL fcn(m, n, x_temp, temp_fvec, iflag, iter)
            WRITE(6,*) ' '
            WRITE(6,*) '  New Minimum at ',ibest,fnorm_min
            WRITE(6,*) ' '
         END IF
      END IF
      
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)  
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi) 
      CALL MPI_BCAST(fnorm_min,1,MPI_REAL8,master,MPI_COMM_STEL,
     1               ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF
      
      
      !Preform Evolutionary Computation
      DO iter = 2, maxfev
         IF (myid == master) THEN
            WRITE(6,*) ' '
            WRITE(6,*) '  Generation ',iter
            WRITE(6,*) ' '
            ! Prep Crossover
            ALLOCATE (lcross(NP,n), stat=ierr)
            IF (ierr .ne. 0) STOP 'DE2_Evolve Error ALLOC(5)'
            IF (CR_XC == 1.0) THEN
               lcross = .TRUE.
            ELSE IF (CR_strategy == 1) THEN
               DO i = 1, NP
                  DO j = 1, n
                     CALL random_number(rand_C1)
                     lcross(i,j) = rand_C1 < CR_XC
                  END DO
               END DO
               DO i = 1, NP
                  IF (.not. ANY(lcross(i,:))) THEN
                     CALL random_number(rand_C1)
                     j = 1 + rand_C1*(n-1)
                     lcross(i,j) = .TRUE.
                  END IF
               END DO
            ELSE
               lcross = .FALSE.
               DO i = 1, NP
                  CALL random_number(rand_C1)
                  j = 1 + rand_C1*(n-1)
                  DO
                     lcross(i,j) = .TRUE.
                     j = MOD(j,n) + 1
                     CALL random_number(rand_C1)
                     IF (rand_C1 >= CR_XC .or. lcross(i,j)) EXIT
                  END DO
               END DO
            END IF
            
            ! Setup Mutation
            ALLOCATE(a1(NP),a2(NP),a3(NP),a4(NP),a5(NP),stat=ierr)
            IF (ierr .ne. 0) STOP 'DE2_Evolve Error ALLOC(6)'
            DO i = 1, NP
               CALL random_number(rand_C1)
               a1(i) = 1 + rand_C1*(NP-1)
            END DO
            a2 = a1
            DO WHILE (ANY(a2==a1))
               DO i =1, NP
                  CALL random_number(rand_C1)
                  IF (a2(i) == a1(i)) a2(i) = 1 + rand_C1*(NP-1)
               END DO
            END DO
            IF( strategy >= 2) THEN
               a3 = a1
               DO WHILE (ANY(a3==a1))
                  DO i =1, NP
                     CALL random_number(rand_C1)
                     IF (a3(i) == a1(i)) a3(i) = 1 + rand_C1*(NP-1)
                  END DO
               END DO
            END IF
            IF( strategy >= 4) THEN
               a4 = a1
               DO WHILE (ANY(a4==a1))
                  DO i =1, NP
                     CALL random_number(rand_C1)
                     IF (a4(i) == a1(i)) a4(i) = 1 + rand_C1*(NP-1)
                  END DO
               END DO
            END IF
            IF( strategy == 5) THEN
               a5 = a1
               DO WHILE (ANY(a4==a1))
                  DO i =1, NP
                     CALL random_number(rand_C1)
                     IF (a5(i) == a1(i)) a5(i) = 1 + rand_C1*(NP-1)
                  END DO
               END DO
            END IF
            
            ! Preload Descendents
            x_new = x_array
            SELECT CASE (strategy)
               CASE(1,3,4)
                  DO i = 1, NP
                     WHERE(lcross(i,:)) x_new(i,:) = x_temp
                  END DO
            END SELECT
            
            ! Select Mutation strategy
            SELECT CASE (strategy)
               CASE(1)
                  WHERE(lcross) x_new = x_new +
     1                         F_XC*(x_array(a1,:)-x_array(a2,:))
               CASE(2)
                  WHERE(lcross) x_new = x_array(a3,:) +
     1                         F_XC*(x_array(a1,:)-x_array(a2,:))
               CASE(3)
                  WHERE(lcross) x_new = x_array +
     1                         F_XC*( x_new - x_array 
     2                         + x_array(a1,:) - x_array(a2,:))
               CASE(4)
                  WHERE(lcross) x_new = x_new +
     1                         F_XC*( x_array(a1,:) - x_array(a2,:) 
     2                         + x_array(a3,:) - x_array(a4,:))
               CASE(5)
                  WHERE(lcross) x_new = x_array(a5,:) +
     1                         F_XC*( x_array(a1,:) - x_array(a2,:) 
     2                         + x_array(a3,:) - x_array(a4,:))
               CASE DEFAULT
                  WHERE(lcross) x_new = x_array(a3,:) +
     1                         F_XC*(x_array(a1,:)-x_array(a2,:))
            END SELECT
      
            ! Confine to domain
            DO i = 1, NP
               x_new(i,:) = MAX(MIN(x_new(i,:),XCmax),XCmin)
            END DO
            
            
            ! DEALLOCATE
            DEALLOCATE(lcross)
            DEALLOCATE(a1,a2,a3,a4,a5)
         END IF
      
!DEC$ IF DEFINED (MPI_OPT)
         ! BROADCAST The X_Array
         CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi) 
         CALL MPI_BCAST(x_new,NP*n,MPI_REAL8,master,MPI_COMM_STEL,
     1               ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF
         
         ! Evaluate Population (new)
!DEC$ IF DEFINED (MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
         IF (myid == master) THEN
            numsent = 0
            DO j = 1,MIN(numprocs-1,NP)
               CALL MPI_SEND(j,1,MPI_INTEGER,j,j,MPI_COMM_STEL,
     1                        ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) 
     1                        CALL mpi_stel_abort(ierr_mpi)
               numsent = numsent + 1
            END DO
            DO j = 1, NP
               CALL MPI_RECV(temp_fvec,m,MPI_REAL8,MPI_ANY_SOURCE,
     1                       MPI_ANY_TAG,MPI_COMM_STEL,status,
     2                       ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) 
     1                        CALL mpi_stel_abort(ierr_mpi)
               sender = status(MPI_SOURCE)
               i      = status(MPI_TAG)
               fval_array(i,:) = temp_fvec
               fnorm = SUM(temp_fvec*temp_fvec)
               fnorm_new(i) = fnorm
               write (6, '(2x,i6,8x,i3,7x,1es12.4)') i,
     1             sender, fnorm
               CALL FLUSH(6)
               IF (numsent < NP) THEN
                  numsent = numsent + 1
                  CALL MPI_SEND(numsent, 1, MPI_INTEGER,
     1                      sender, numsent, MPI_COMM_STEL, ierr_mpi)
                  IF (ierr_mpi /= MPI_SUCCESS) 
     1                        CALL mpi_stel_abort(ierr_mpi)
               ELSE
                  i = 0
                  CALL MPI_SEND(i, 1, MPI_INTEGER,
     1                      sender, numsent, MPI_COMM_STEL, ierr_mpi)
                  IF (ierr_mpi /= MPI_SUCCESS) 
     1                        CALL mpi_stel_abort(ierr_mpi)
               END IF
            END DO
         ELSE IF (myid < NP) THEN
            DO
               CALL MPI_RECV(j,1,MPI_INTEGER,master,MPI_ANY_TAG,
     1                       MPI_COMM_STEL,status,ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) 
     1                        CALL mpi_stel_abort(ierr_mpi)
               IF (j > 0) THEN
                  x_temp     = x_new(j,:)
                  temp_fvec  = 0
                  iflag = j
                  CALL fcn(m, n, x_temp, temp_fvec, iflag, iter)
                  fnorm = SUM(temp_fvec*temp_fvec)
                  CALL MPI_SEND(temp_fvec,m, MPI_REAL8, master,
     1                       j, MPI_COMM_STEL, ierr_mpi)
                  IF (ierr_mpi /= MPI_SUCCESS) 
     1                        CALL mpi_stel_abort(ierr_mpi)
               ELSE
                  EXIT
               END IF
            END DO
         ELSE
         END IF
         
         CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ELSE
         DO i = 1, NP
            x_temp = x_new(i,:)
            iflag = i
            CALL fcn(m, n, x_temp, temp_fvec, iflag, iter)
            fval_array(i,:) = temp_fvec
            fnorm_new(i)  = SUM(temp_fvec*temp_fvec)
         END DO
!DEC$ ENDIF
      
         ! Output to the xvec.dat file
         IF (myid == master) THEN
            CALL safe_open(iunitx,ierr,'xvec.dat','unknown',
     1                     'formatted',ACCESS_IN='APPEND')
            IF (ierr .ne. 0) STOP 'DE2_Evolve Error OPEN(xvec.dat)'
            DO i = 1, NP
               WRITE(iunitx,'(2(2X,I5.5))') n,iter
               WRITE(iunitx,'(10ES22.12E3)') x_new(i,1:n)
               WRITE(iunitx,'(ES22.12E3)') fnorm_new(i)
            END DO
            CLOSE(iunitx)
         END IF

         
         ! Output the best value
         ibest     = MINLOC(fnorm_new,DIM = 1)   
         fnorm     = fnorm_new(ibest)
         IF (fnorm_min > fnorm) THEN
!            fnorm_min = fnorm
            ! Save the output
            IF (myid == master) THEN
               x_temp = x_array(ibest,:)
               iflag = 0
               CALL fcn(m, n, x_temp, temp_fvec, iflag, iter)
               iflag = GADE_CLEANUP
               CALL fcn(m, n, x_temp, temp_fvec, iflag, iter)
               WRITE(6,*) ' '
               WRITE(6,*) '  New Minimum at ',ibest,fnorm_min
               WRITE(6,*) ' '
            END IF
         END IF

         ! Now find the best
         IF (myid == master) THEN
            j = 0
            ibest     = MINLOC(fnorm_new,DIM = 1)
            fnorm     = fnorm_array(ibest)
            DO i = 1, NP
               IF (fnorm_new(i) < fnorm_array(i)) THEN
                  j = j + 1
                  x_array(i,:) = x_new(i,:)
                  fnorm_array(i) = fnorm_new(i)
               END IF
               IF (fnorm_new(i) >= 1.0E12) THEN
                  x_array(i,:) = x_new(ibest,:)
                  fnorm_array(i) = fnorm_new(ibest)
               END IF
            END DO
            WRITE(6,*) ' Selecting,',j,' improvements'
         END IF

!DEC$ IF DEFINED (MPI_OPT)
         ! Now broadcast this stuff
         CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
         CALL MPI_BCAST(x_array,NP*n,MPI_REAL8,
     1                  master,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi) 
         CALL MPI_BCAST(fnorm_array,NP,MPI_REAL8,
     1                  master,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi) 
         CALL MPI_BCAST(ibest,1,MPI_INTEGER,
     1                  master,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi) 
         CALL MPI_BCAST(fnorm,1,MPI_REAL8,
     1                  master,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi) 
         CALL MPI_BCAST(fnorm_min,1,MPI_REAL8,master,MPI_COMM_STEL,
     1               ierr_mpi)
         IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF
         
         ! WRITE RESTART FILE
         IF (myid == master) THEN
            REWIND (unit=irestart)
            WRITE(iRESTART,*) n
            DO i=1,NP
               WRITE(iRESTART,*) i,fnorm_array(i),
     1                           (x_array(i,j),j=1, n)
            END DO
            CALL FLUSH(iRESTART)
            
         END IF
      END DO
      
      x_temp = x_array(ibest,1)
      x = x_temp
      
      
      DEALLOCATE(fnorm_array,temp_fvec)
      DEALLOCATE(x_array)
      IF (myid == master) DEALLOCATE(fnorm_new)
      
      RETURN
 1327 FORMAT (/,' Beginning Differential Evolution II',/,
     1        ' Number of Processors: ',i4,/,
     1        ' Population Size: ',i4,/,
     1        ' Strategy: ',i4,/,
     1        ' Crossover Strategy: ',i4,//,
     2        70('='),/,2x,'Member ID',3x,'Processor',7x,'Chi-Sq')
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE DE2_Evolve
