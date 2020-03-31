      SUBROUTINE ga_initial(istart,npossum,ig2sum,filename,myid)
c#######################################################################
c
c  This subroutine sets up the program by generating the g0, g1 and
c  ig2 arrays, and counting the number of chromosomes required for the
c  specified input.  The subroutine also initializes the random number
c  generator, parent and iparent arrays (reads the ga.restart file).
      USE ga_mod
      USE safe_open_mod
      USE mpi_params, ONLY: master, MPI_COMM_STEL
      USE mpi_inc
      IMPLICIT NONE
      INTEGER :: ierr
      INTEGER :: istart, npossum, ig2sum, myid
      INTEGER :: i, j, k, l, itemp, istat
      CHARACTER(LEN=100) :: filename
      CHARACTER(LEN=200) :: temp
      REAL(rprec) :: rand
c
c
      DO i=1,nparam
         g0(i)=par_min(i)
         pardel(i)=par_max(i)-par_min(i)
         itemp=2**nposibl(i)
         g1(i)=pardel(i)/(itemp-1)
      END DO

      DO  i=1,nparam
         ig2(i)=nposibl(i)
      END DO
c
c  Count the total number of chromosomes (bits) required
      nchrome=0
      npossum=0
      ig2sum=0
      DO 9 i=1,nparam
         nchrome=nchrome+ig2(i)
         npossum=npossum+2**nposibl(i)
         ig2sum=ig2sum+(2**ig2(i))
 9    CONTINUE
      IF (nchrome.gt.nchrmax) THEN
         IF (myid .eq. master) THEN
            WRITE(6,1800) nchrome
            WRITE(iunit_ga_out,1800) nchrome
            CLOSE(iunit_ga_out)
         END IF
         STOP
      END IF
c
      IF (npossum.lt.ig2sum .and. microga.ne.0
     1   .and. myid.eq.master) THEN
         WRITE(6,2100)
         WRITE(iunit_ga_out,2100)
      END IF
c
c  Initialize random number generator
      CALL ran3(idum,rand)
c
      IF(irestrt.eq.0) THEN
c  Initialize the random distribution of parameters in the individual
c  parents when irestrt=0.
         istart=1
         DO 10 i=1,npopsiz
            DO 15 j=1,nchrome
               CALL ran3(1,rand)
               iparent(j,i)=1
               IF(rand.lt.0.5d0) iparent(j,i)=0
 15         CONTINUE
 10      CONTINUE
         IF (npossum.lt.ig2sum) CALL ga_possibl(parent,iparent,myid)
c  insert unique individual
         IF (unique_ind .gt. 0) THEN
            DO i=1, nparam
               CALL ga_code(unique_ind, i, parent, iparent)
            END DO
         END IF

      ELSE
c  If irestrt.ne.0, READ from restart file.
         IF (myid .eq. master) THEN
            iunit_ga_restart = 25
            temp = "../ga_restart." // filename
            CALL safe_open(iunit_ga_restart, istat,
     1               TRIM(temp), 'unknown', 'formatted')
            READ (iunit_ga_restart,*) istart,npopsiz
            DO j=1,npopsiz
               READ(iunit_ga_restart,*) k,(iparent(l,j),l=1,nchrome)
            END DO
            CLOSE (iunit_ga_restart)
         END IF
!DEC$ IF DEFINED (MPI_OPT)
         CALL MPI_BCAST(istart, 1, MPI_INTEGER, master,
     1     MPI_COMM_STEL, ierr)
         CALL MPI_BCAST(npopsiz, 1, MPI_INTEGER, master,
     1     MPI_COMM_STEL, ierr)
      DO l = 1, nchrome
         IF (myid .eq. master) fitness = iparent(l,:)
         CALL MPI_BCAST(fitness, indmax, MPI_REAL8, master,
     1        MPI_COMM_STEL, ierr)
         IF (myid .ne. master) iparent(l,:) = fitness
      END DO
!DEC$ ENDIF

      END IF
c
      IF(irestrt.ne.0) CALL ran3(idum-istart,rand)
c
 1800 FORMAT(1x,'ERROR: nchrome > nchrmax.  Set nchrmax = ',i6)
 2000 FORMAT(1x,'ERROR: you have a parameter with a number of '/
     +       1x,'   possibilities > 2**30!  if you really desire this,'/
     +       1x,'   change the do loop 7 statement and recompile.'//
     +       1x,'   you may also need to alter the code to work with'/
     +       1x,'   real numbers rather than integer numbers; fortran'/
     +       1x,'   does not like to compute 2**j when j>30.')
 2100 FORMAT(1x,'WARNING: for some cases, a considerable performance'/
     +       1x,'   reduction has been observed when running a non-'/
     +       1x,'   optimal number of bits with the micro-GA.'/
     +       1x,'   If possible, use values for nposibl of 2**n,'/
     +       1x,'   e.g. 2, 4, 8, 16, 32, 64, etc.  See ReadMe file.')
c
      
      END SUBROUTINE ga_initial
