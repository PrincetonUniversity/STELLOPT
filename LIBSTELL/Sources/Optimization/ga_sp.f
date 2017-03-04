      SUBROUTINE ga_sp(fcn, nopt, fvec, best, filename, nfev, iflag,
     1                 max_num_processors, myid)
      USE ga_mod
      USE safe_open_mod
      USE mpi_params, ONLY: master
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      include 'mpif.h'                                       !mpi stuff
!DEC$ ENDIF
      EXTERNAL fcn

      INTEGER :: nopt
      REAL(rprec), DIMENSION(nopt) :: fvec

      INTEGER :: kount, npossum, ig2sum, istart, istore
      INTEGER :: ncross, ipick, mate1, mate2, istat
      INTEGER :: i, j, nfev, iflag, max_num_processors, myid
      REAL(rprec), INTENT(in) :: best
      REAL(rprec), SAVE :: fbar, evals
      CHARACTER(LEN=*) :: filename
      CHARACTER(LEN=LEN(filename)+10) :: temp
c
c     CALL input
c
c  Perform necessary initialization and READ the ga.restart file.
      CALL ga_initial(istart,npossum,ig2sum,filename,myid)
c
c  $$$$$ Main generational processing loop. $$$$$
      kount=0
      nfit_eval=nfev
      istore=0
      iunit_ga_out = 24
      IF (myid .eq. master) THEN
        temp = "ga_out." // filename
        CALL safe_open(iunit_ga_out, istat, TRIM(temp),
     1                'unknown', 'formatted')
      END IF


      DO 20 i=istart,maxgen+istart-1
         iflag=-1
         IF (myid .eq. master) THEN
            WRITE (6,1111) i
            WRITE (iunit_ga_out,1111) i
c           WRITE (iunit_ga_out,1050)
c
c  Evaluate the population, assign fitness, establish the best
c  individual, and write output information.
            WRITE(6,*) 'pre ga_evalout', max_num_processors
            WRITE(6,*) fbar,best,nopt,nfev,max_num_processors,iflag
         END IF
         CALL ga_evalout(fbar, best, fcn, nopt, fvec, nfev,
     >        max_num_processors, iflag, myid)
         istore=istore+1
         geni(istore)  = i
         genavg(istore)=fbar
         genmax(istore)=best
         IF (npopsiz.eq.1 .or. iskip.ne.0) THEN
            IF (myid .eq. master) CLOSE(iunit_ga_out)
            CALL ga_restart(i,istart,kount,filename, myid)
            RETURN
         END IF
c
c  niching
         IF (iniche.ne.0) CALL ga_niche(myid)
c
c  selection, crossover and mutation
         ncross=0
         ipick=npopsiz
         DO 45 j=1,npopsiz,nchild
c
c  Perform selection.
            CALL ga_selectn(ipick,j,mate1,mate2)
c
c  Now perform crossover between the randomly selected pair.
            CALL crosovr(ncross,j,mate1,mate2)
 45      CONTINUE

         IF (myid .eq. master) THEN
            WRITE(6,1225) ncross
            WRITE(iunit_ga_out,1225) ncross
         END IF
c
c  Now perform random mutations.  If running micro-GA, skip mutation.
         IF (microga.eq.0) CALL ga_mutate (myid)
c
c  Write child array back into parent array for new generation.  Check
c  to see IF the best parent was replicated.
         CALL ga_newgen(npossum,ig2sum,myid)
c
c  Implement micro-GA if enabled.
         IF (microga.ne.0) CALL ga_micro(i,npossum,ig2sum,myid)
c
c  Write to restart file.
         CALL ga_restart(i,istart,kount,filename,myid)
 20   CONTINUE

c  $$$$$ End of main generational processing loop. $$$$$

      IF (myid .eq. master) THEN
         WRITE(iunit_ga_out,3000)
         DO 100 i=1,maxgen
            evals = npopsiz*geni(i)
            WRITE(iunit_ga_out,3100) geni(i),evals,genavg(i),genmax(i)
 100     CONTINUE
         CLOSE (iunit_ga_out)
      END IF

 1050 FORMAT(1x,'    Binary Code',16x,'Parameter Values and  Fitness')
 1111 FORMAT(//,'#################  Generation',i5,
     1        '  #################')
 1225 FORMAT(/'  Number of Crossovers      =',i5)
 3000 FORMAT(2x//'Summary of Output'/
     +       2x,'Generation   Evaluations   Avg.Fitness   Best Fitness')
 3100 FORMAT(2x,3(e11.4,4x),e12.5)

      END SUBROUTINE ga_sp
