      SUBROUTINE DE_Evolve(obj, Dim_XC, XCmin, XCmax, VTR, NP, itermax,
     1     F_XC,
     2     CR_XC, strategy, CR_strategy, refresh, iwrite, irestart,
     3     num_proc, lrestart,
     4     best_XC, bestval, nfeval)
!.......................................................................
!
! Differential Evolution for Optimal Control Problems
!
!.......................................................................
!  This Fortran 90 PROGRAM translates from the original MATLAB
!  version of differential evolution (DE). This FORTRAN 90 code
!  has been tested on Compaq Visual Fortran v6.1.
!  Any user new to the DE are encouraged to READ the article of Storn and Price.
!
!  Refences:
!  Storn, R., and Price, K.V., (1996). Minimizing the REAL FUNCTION of the
!    ICEC''96 contest by differential evolution. IEEE conf. on Evolutionary
!    Computation, 842-844.
!
!  This Fortran 90 PROGRAM is written by Dr. Feng-Sheng Wang
!  Department of Chemical Engineering, National Chung Cheng University,
!  Chia-Yi 621, Taiwan, e-mail: chmfsw@ccunix.ccu.edu.tw
!
!  Modified by M. Zarnstorff, PPPL, to implement ALL strategies and make random
!  selection process be compatible with c versions and restart capability.
!  Also, added multi-processor INTERFACE to STELLOPT suite.
!  Modified by S. Hirshman, ORNL to implement MPI INTERFACE to STELLOPT.
!.........................................................................
!                obj : The User provided file for evlauting the objective function.
!                      SUBROUTINE obj(xc,fitness)
!                      WHERE "xc" is the REAL decision PARAMETER vector.(input)
!                            "fitness" is the fitness value.(output)
!             Dim_XC : Dimension of the REAL decision parameters.
!      XCMIN(Dim_XC) : The lower bound of the REAL decision parameters.
!      XCMAX(Dim_XC) : The upper bound of the REAL decision parameters.
!                VTR : The EXPected fitness value to reach.
!                 NP : Population SIZE.
!            itermax : The maximum number of iteration.
!               F_XC : Mutation scaling factor for REAL decision parameters.
!              CR_XC : Crossover factor for REAL decision parameters.
!           strategy : The strategy of the mutation operations is used in HDE.
!        CR_strategy : cross-over strategy 0=exponential, 1=binomial
!            refresh : The intermediate output will be produced after "refresh"
!                      iterations. No intermediate output will be produced IF
!                      "refresh < 1".
!             iWRITE : The unit specfier for writing to an EXTERNAL data file.
!           irestart : unit number for writing (and reading) the restart file
!           num_proc : number of processors to USE in evaluating population
!           lrestart : =T THEN READ in restart file, =F ignore ANY restart files

!    best_XC(Dim_XC) : The best REAL decision parameters.  If non-zero, this
!                      array also CONTAINS an initial guess that is included
!                      in the population
!            bestval : The best objective FUNCTION.
!             nfeval : The number of FUNCTION CALL.
!

      USE stel_kinds
      USE de_mod, ONLY: ui_XC, nopt
      USE gade_mod, ONLY : save_space, gade_cleanup
      USE mpi_params
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      include 'mpif.h'                                       !mpi stuff
!DEC$ ENDIF

      INTEGER, INTENT(in) :: NP, Dim_XC, itermax, strategy,
     1                       CR_strategy, iwrite, irestart, refresh
      REAL(rprec), INTENT(in) :: VTR, CR_XC, F_XC
      REAL(rprec), DIMENSION(Dim_XC), INTENT(in) :: XCmin, XCmax
      REAL(rprec), DIMENSION(Dim_XC), INTENT(inout) :: best_XC
      REAL(rprec), INTENT(out) :: bestval
      INTEGER, INTENT(out) :: nfeval, num_proc
      LOGICAL :: lrestart
      EXTERNAL  obj

      LOGICAL, DIMENSION(NP,Dim_XC) :: lcross
      INTEGER :: ierr, numsent
      INTEGER :: i, j, iter, k, n, r_free, seed_size
      INTEGER, DIMENSION(1) :: ibest
      INTEGER, DIMENSION(NP) :: a1, a2, a3, a4, a5
      INTEGER, ALLOCATABLE :: seed (:)

      REAL(rprec) :: s, val_mean, val_max, temp
      REAL(rprec), DIMENSION(NP) :: val, rand, new_val
      REAL(rprec), DIMENSION(Dim_XC) :: bestmemit_XC
      REAL(rprec), DIMENSION(Dim_XC) :: rand_C1
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: pop_XC, rand_XC
      REAL(rprec), DIMENSION(nopt) :: fvec
!DEC$ IF DEFINED (MPI_OPT)
      REAL(rprec), DIMENSION(Dim_XC) :: temp_var
      INTEGER :: status(MPI_STATUS_size)                     !mpi stuff
      INTEGER :: sender, anstype,iflag
!DEC$ ENDIF

      INTRINSIC MAX, MIN, RANDOM_NUMBER, MOD, ABS, ANY, ALL, MAXLOC,
     1          MINLOC, MAXVAL, MINVAL

      ierr = 0
!DEC$ IF DEFINED (MPI_OPT)
      ierr_mpi = 0
      call MPI_COMM_RANK (MPI_COMM_STEL, myid, ierr_mpi)        !mpi stuff
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)                 !mpi stuff
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF

      ALLOCATE (ui_XC(NP, Dim_XC), pop_XC(NP, Dim_XC),
     1          rand_XC(NP, Dim_XC), stat=ierr)
      IF (ierr /= 0) STOP 'ALLOCATION ERROR IN DE_EVOLVE'

!-----Initialize the random-number seed stuff

      val = 0.0; ui_XC = 0; pop_XC = 0.0; rand_XC = 0.0
      
      CALL RANDOM_SEED ( )  ! Processor reinitializes the seed
      IF (myid .eq. master) THEN
         IF (ierr .ne. 0) STOP 'Allocation error 1 in DE_Evolve'
         CALL RANDOM_SEED (SIZE = seed_size)  ! SIZE of seed array
         ALLOCATE (seed(seed_size), stat=ierr)
         IF (ierr .ne. 0) STOP 'Allocation error 1 in DE_Evolve'
         CALL RANDOM_SEED (GET=seed(1:seed_size)) ! Gets the current seed

!!-----Initialize a population --------------------------------------------!!

         pop_XC=0
         DO i=1,NP
            CALL random_number(rand_C1)
            pop_XC(i,:)=XCmin+rand_C1*(XCmax-XCmin)
         END DO

         val = 0

         IF( ANY(best_XC /= 0._dp)) pop_XC(1,:) = best_XC   ! init with given CASE

!!----- Restart IF requested ----

         IF (lrestart) THEN
            WRITE(6,'(A)') 'Restarting from saved file'

            READ (irestart,*) r_free
            IF (r_free /= Dim_XC) THEN
               WRITE(6,'(A)') '***DE Restart Error***'
               WRITE(6,'(A)')
     1              '***Mismatch on number of free parameters! ***'
               STOP
            END IF

            READ (irestart, *) seed(:seed_size)
            CALL RANDOM_SEED (PUT=seed(1:seed_size)) ! Sets seed from array

            DO i=1,NP
               READ (irestart,*,end=100)
     1              k,val(i),(pop_XC(i,j),j=1, Dim_XC)
            END DO

 100        REWIND (unit=irestart)
            ! We assume a restart file exists in the directory
         ELSE
            j=-1  ! FLAG Singletask
            nfeval=0
            fvec = 0
            CALL obj(nopt, Dim_XC, best_XC, fvec, j, nfeval)
            nfeval=0
            j=gade_cleanup
            CALL obj(nopt, Dim_XC, best_XC, fvec, j, nfeval)
            
            WRITE (6, 1327) num_proc, NP, strategy, CR_strategy
 1327 FORMAT (/,' Beginning Differential Evolution',/,
     1        ' Number of Processors: ',i4,/,
     1        ' Population Size: ',i4,/,
     1        ' Strategy: ',i4,/,
     1        ' Crossover Strategy: ',i4,//,
     2        70('='),/,2x,'Iteration',3x,'Processor',7x,'Chi-Sq')
     
            bestval = SUM(fvec*fvec)
            WRITE(6, '(2x,i6,8x,i3,7x,1es12.4)')
     1                           0, myid, bestval
         END IF
      END IF
      

!DEC$ IF DEFINED (MPI_OPT)
!!ALL processors need same random variables and initial values
      
   
      CALL MPI_BCAST(pop_XC,NP*Dim_XC,MPI_REAL8,master,
     1               MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
     
      CALL MPI_BCAST(seed_size, 1, MPI_INTEGER, master,
     1     MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi) 
      
      IF (.not.ALLOCATED(seed)) ALLOCATE (seed(seed_size), stat=ierr)
      IF (ierr /=0) PRINT *,myid,'ALLOCATION ERROR'
      
      CALL MPI_BCAST(seed, seed_size, MPI_REAL8, master,
     1     MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
      
      CALL MPI_BCAST(val, NP, MPI_REAL8, master,
     1     MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
      
      CALL FLUSH(6); CALL FLUSH(iwrite); CALL FLUSH(irestart)
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi) 
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
      
!DEC$ ENDIF 
      
!!--------------------------------------------------------------------------!!
!!------Evaluate fitness functions and find the best member-----------------!!
      nfeval=1

      ui_XC(:,:) = pop_XC(:,:)    ! prepare for initial evaluation
   

      IF( ANY(val(:NP) == 0._dp) ) THEN
         !CALL DE_Evaluate( num_proc, obj, val, NP, Dim_XC, nfeval)
!DEC$ IF DEFINED (MPI_OPT)
         IF (myid == master) THEN
            DO j = 1,MIN(numprocs-1,np)
               CALL MPI_SEND(j,1,MPI_INTEGER,j,j,MPI_COMM_STEL,
     1                        ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) 
     1                        CALL mpi_stel_abort(ierr_mpi)
               numsent = numsent + 1
            END DO
            DO j = 1, NP
               CALL MPI_RECV(temp,1,MPI_REAL8,MPI_ANY_SOURCE,
     1                       MPI_ANY_TAG,MPI_COMM_STEL,status,
     2                       ierr_mpi)
               sender     = status(MPI_SOURCE)
               anstype    = status(MPI_TAG)
               val(anstype) = temp
               write (6, '(2x,i6,8x,i3,7x,1es12.4)') anstype,
     1             sender, val(anstype)
               CALL FLUSH(6)
               IF (numsent < NP) THEN
                  numsent = numsent + 1
                  CALL MPI_SEND(numsent, 1, MPI_INTEGER,
     1                      sender, numsent, MPI_COMM_STEL, ierr_mpi)
               ELSE
                  anstype = 0
                  CALL MPI_SEND(anstype, 1, MPI_INTEGER,
     1                      sender, numsent, MPI_COMM_STEL, ierr_mpi)
               END IF
            END DO
         ELSE IF (myid <= NP) THEN
            DO
               j=-1
               CALL MPI_RECV(j,1,MPI_INTEGER,master,MPI_ANY_TAG,
     1                       MPI_COMM_STEL,status,ierr_mpi)
               IF (j > 0) THEN
                  temp_var = ui_XC(j,:)
                  fvec     = 0
                  iflag    = j
                  CALL obj(nopt, Dim_XC, temp_var, fvec, j, nfeval)
                  temp = SUM(fvec*fvec)
                  CALL MPI_SEND(temp,1, MPI_REAL8, master,
     1                       j, MPI_COMM_STEL, ierr_mpi)
               ELSE
                  EXIT
               END IF
            END DO
         ELSE
         END IF
         CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi) 
         CALL MPI_BCAST(val, NP, MPI_REAL8, master,
     1                  MPI_COMM_STEL, ierr_mpi)
!DEC$ ENDIF 
      ELSE IF (myid .eq. master) THEN
         WRITE(6,*) 'Using evaluation from restart file'
      END IF

c     DO i=1,NP
c        CALL obj(pop_XC(i,:), val(i))
c        nfeval=nfeval+1
c     END DO

      ibest = MINLOC(val)
      bestval = val(ibest(1))
      bestmemit_XC=pop_XC(ibest(1),:)
      best_XC=bestmemit_XC

      IF (myid .eq. master)
     1    WRITE(6,*)'initial best=',bestval,' at loc ',ibest
     
      ! Write the progress file
      !CALL FSEEK(iwrite,0,2,ierr)
      IF (myid .eq. master) THEN
         WRITE(iwrite,'(A,3(2X,I5.5))') 'ITER ',iter, NP, Dim_XC
         DO i=1,NP
            WRITE(iwrite,*) i,val(i),(pop_XC(i,j),j=1, Dim_XC)
         END DO
      END IF
      
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi) 
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF 

!!--------------------------------------------------------------------------!!
!!------Perform evolutionary computation------------------------------------!!
!!
      DO iter=1, itermax

        IF (myid .ne. master) GOTO 1000
!---- prep for crossover  -----------------------------------------!
        IF( CR_XC == 1.0_dp ) THEN
           lcross = .true.

        ELSE IF( CR_strategy == 1) THEN
           CALL random_number(rand_XC)
           lcross = rand_XC < CR_XC

           DO i=1, NP
              IF( .not. ANY(lcross(i,:))) THEN
                 n = 1 + rand_XC(i,1)*(Dim_XC-1)
                 lcross(i,n) = .true.
              END IF
           END DO

        ELSE   ! IF cr_strategy = 0 or ANYthing ELSE
           lcross = .false.
           CALL random_number(rand)
           DO i=1, NP
              n = 1 + rand(i)*(Dim_XC-1)
              DO
                 lcross(i,n) = .true.
                 n = MOD(n,Dim_XC) + 1
                 CALL random_number(s)
                 IF(s >= CR_XC .or. lcross(i,n) ) EXIT
              END DO
           END DO
        END IF


!!------Setup Mutation parents----------------------------------------------!!
        IF( NP < 2) STOP 'population too small for Diff. Evolution'

        CALL random_number(rand)
        a1 = 1 + rand*(NP-1)

        a2 = a1
        DO WHILE (ANY(a2 == a1) )
           CALL random_number(rand)
           WHERE (a2 == a1) a2 = 1 + rand*(NP-1)
        END DO

        IF( strategy == 2 .or. strategy == 4 .or. strategy == 5) THEN
           IF( NP < 3) STOP 'population too small for DE strategy'
           a3 = a1
           DO WHILE (ANY(a3 == a1 .or. a3 == a2))
              CALL random_number(rand)
              WHERE (a3 == a1 .or. a3 == a2)
     1            a3 = 1 + rand*(NP-1)
           END DO
        END IF

        IF( strategy == 4 .or. strategy == 5) THEN
           IF( NP < 4) STOP 'population too small for DE strategy'
           a4 = a1
           DO WHILE (ANY(a4 == a1 .or. a4 == a2 .or. a4 == a3))
              CALL random_number(rand)
              WHERE (a4 == a1 .or. a4 == a2 .or. a4 == a3)
     1            a4 = 1 + rand*(NP-1)
           END DO
        END IF

        IF( strategy == 5) THEN
           IF( NP < 5) STOP 'population too small for DE strategy'
           a5 = a1
           DO WHILE (ANY(a5 == a1 .or. a5 == a2 .or. a5 == a3
     1             .or. a5 == a4))
              CALL random_number(rand)
              WHERE (a5 == a1 .or. a5 == a2 .or. a5 == a3 .or.
     1             a5 == a4)
     2            a5 = 1 + rand*(NP-1)
           END DO
        END IF


!---- preload new descendents   -----------------------------------------!

        ui_XC = pop_XC       ! start with new pop. = old pop.


!---- for strategies using best previous CASE, preload best CASE into new pop.

        IF( strategy == 1 .or. strategy == 3 .or. strategy == 4) THEN
           DO i=1,NP
              WHERE(lcross(i,:)) ui_XC(i,:) = bestmemit_XC(:)
           END DO
        END IF

!---- SELECT and DO a  mutation strategy-----------------------------------!

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

        SELECT CASE (strategy)

        CASE (1)
           WHERE(lcross)     ! note: in these locations ui_XC = best CASE
     1     ui_XC=ui_XC+F_XC*(pop_XC(a1,:)-pop_XC(a2,:))

        CASE DEFAULT
           WHERE(lcross)
     1     ui_XC=pop_XC(a3,:)+F_XC*(pop_XC(a1,:)-pop_XC(a2,:))

        CASE (3)
           WHERE(lcross)     ! note: in these locations ui_XC = best CASE
     1     ui_XC=pop_XC+F_XC*(ui_XC-pop_XC+pop_XC(a1,:)-pop_XC(a2,:))

        CASE (4)
           WHERE(lcross)     ! note: in these locations ui_XC = best CASE
     1     ui_XC=ui_XC+F_XC*(pop_XC(a1,:)-pop_XC(a2,:)+pop_XC(a3,:)
     2                     -pop_XC(a4,:))

        CASE (5)
           WHERE(lcross)
     1     ui_XC=pop_XC(a5,:)+F_XC*(pop_XC(a1,:)-pop_XC(a2,:)+
     2                              pop_XC(a3,:)-pop_XC(a4,:))

        END SELECT


!------Confine each of feasible individuals in the lower-upper bound-------!!

        DO i=1,NP
           ui_XC(i,:)=MAX(MIN(ui_XC(i,:),XCmax),XCmin)
        END DO

 1000   CONTINUE
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi) 
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
!!ALL processors need same ui_XC values going into DE_Evaluate
      DO i = 1, NP
         IF (myid .eq. master) temp_var = ui_XC(i,:)
         CALL MPI_BCAST(temp_var, Dim_XC, MPI_REAL8, master,
     1        MPI_COMM_STEL, ierr)
         IF (myid .ne. master) ui_XC(i,:) = temp_var
      END DO
!DEC$ ENDIF


!!--------------------------------------------------------------------------!!
!!------Evaluate fitness functions and find the best member-----------------!!

        CALL DE_Evaluate( num_proc, obj, new_val, NP, Dim_XC, nfeval)

        j = 0
        DO i=1,NP
           IF (new_val(i) < val(i)) THEN
              j = j + 1
              pop_XC(i,:)=ui_XC(i,:)
              val(i)=new_val(i)
              IF (val(i) < bestval) THEN
                  bestval=val(i)
                  best_XC=ui_XC(i,:)
              END IF
           END IF
        END DO

        IF (myid .eq. master) WRITE(6,*) 'selecting,',j,' improvements'
!
!        IF (save_space) THEN
!        SPH: NEED THIS CALL SO MASTER HAS best_XC EVALUATION
!
        IF (myid .eq. master) THEN
           j = 1
           nfeval = nfeval + 1
           CALL obj(nopt, Dim_XC, best_XC, fvec, j, nfeval)
           j = gade_cleanup ! Clean up
           CALL obj(nopt, Dim_XC, best_XC, fvec, j, nfeval)
        END IF
!        END IF


!------Write restart file ----------------------------------------

        IF (myid .eq. master) THEN
           WRITE(irestart,*) Dim_XC

           CALL RANDOM_SEED (GET=seed(1:seed_size)) ! Sets seed from array

           WRITE(irestart,*) seed(:seed_size)

           DO i=1,NP
              WRITE(irestart,*) i,val(i),(pop_XC(i,j),j=1, Dim_XC)
           END DO

           REWIND (unit=irestart)
           ! Write the progress file
           WRITE(iwrite,'(A,3(2X,I5.5))') 'ITER ',iter, NP, Dim_XC
           DO i=1,NP
              WRITE(iwrite,*) i,val(i),(pop_XC(i,j),j=1, Dim_XC)
           END DO
           
        END IF

!------Write output summary IF needed------------------------------

        IF( refresh > 0 ) THEN
           IF (MOD(iter,refresh) == 0) THEN
              val_mean = SUM(val(:NP))/NP
              val_max = MAXVAL(val(:NP))

              IF (myid .eq. master) THEN
c          WRITE(unit=iwrite,FMT=203) iter
                 WRITE(6, FMT=203) iter
                 DO i=1,Dim_XC
c                WRITE(unit=iwrite, FMT=202) i, best_XC(i)
                    WRITE(6,FMT=202) i,best_XC(i)
                 END DO
c          WRITE(unit=iwrite, FMT=201) bestval
                 WRITE(6, FMT=201) bestval,val_mean,val_max
              END IF
           END IF
        END IF

        IF ( bestval <= VTR ) THEN
c          WRITE(iwrite, *) ' The best fitness is smaller than VTR'
           IF (myid .eq. master)
     1        WRITE(6, *) 'The best fitness is smaller than VTR'
           EXIT
        END IF

      END DO

      DEALLOCATE (seed, ui_XC, pop_XC, rand_XC)

!!------END the evolutionary computation------------------------------!!
 201  FORMAT(2x, 'bestval=', ES14.7, ' mean=',ES14.7,' MAX=',ES14.7, /)
 202  FORMAT(5x, 'best_XC(', I3, ')=', ES12.5)
 203  FORMAT(2x, 'No. of iteration =', I8)

      END SUBROUTINE DE_Evolve
