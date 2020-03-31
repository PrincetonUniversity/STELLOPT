!-------------------------------------------------------------------------------
!>
!!  \author S. P. Hirshman, V. E. Lynch, K.S.Perumalla, and R. Sanchez
!!  \version 1.0
!!  \date 2006-2009
!!  \brief Computes inverse of large block-tridiagonal system using efficient cyclic-reduction
!!  \bug Please report any bugs to hirshmansp@ornl.gov
!<
!-------------------------------------------------------------------------------
      MODULE CYCLIC_RED

!     SERIAL VERSION WRITTEN BY S.P.Hirshman  (10/14/08)
!     PARALLEL VERSION:         V.E.Lynch     (01/01/09)
!     BLOCKS, NON POWERS OF 2 : K.S.Perumalla (03/01/09)
!
!     PERFORMS CYCLIC REDUCTION FACTORIZATION/SOLUTION OF A BLOCK TRI-DIAGONAL 
!     MATRIX WITH LOWER, DIAGONAL, UPPER BLOCKS LBLK, DBLK, UBLK, RESPECTIVELY
!

      USE stel_kinds
      USE stel_constants, ONLY: one, zero
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'                                       !mpi stuff
      INTEGER, PARAMETER :: tag1=1, tag2=4
      INTEGER :: status(MPI_STATUS_SIZE)                     !mpi stuff
      INTEGER :: numtasks, rank                              !mpi stuff

      INTEGER :: mpireqind !<Index of next MPI_Irecv request
      INTEGER :: mpireqcnt !<Number of MPI_Irecv/Isend requests initiated
      INTEGER :: mpireq(6) !<MPI_Irecv/Isend request handle for top or bot nbr
      INTEGER :: mpierr(6) !<MPI_Irecv/Isend error for top or bot nbr
      INTEGER :: mpistats(MPI_STATUS_SIZE,6) !<MPI_Waitall statuses


      INTEGER :: mblk2                                       !<mpi #elements transferred
      INTEGER :: ranknext                                     !<rank of proc for next row
      INTEGER :: rankprev                                     !<rank of proc for prev row
!DEC$ ENDIF
      INTEGER :: comm_on, comm_off, commRate, omp_on, omp_off
      INTEGER :: tot_time_comm, sum_comm, tot_time_comp, tot_time, &
     &           loc_time_comm, max_tot, tot_time_omp
      INTEGER :: ns0, nsn, nblk_rows, nLast, nFirst
      INTEGER :: timeon, timeoff, timeon1, timeoff1, countRate
      REAL(rprec) :: init_memory, added_memory
    
!>
!!  Defined type to store u,l blocks with minimum memory usage
!<
      TYPE ULSTORE
         INTEGER              :: NumElements    !<last dim of UMAT/LMAT
         REAL(rprec), POINTER :: UMAT(:,:,:)    !<upper block
         REAL(rprec), POINTER :: LMAT(:,:,:)    !<lower block
      END TYPE ULSTORE

      TYPE(ULSTORE), ALLOCATABLE :: OddStorage(:)  !<array of ULSTORES, one at each level
      CONTAINS

!-------------------------------------------------------------------------------
!>
!! Is given integer even?
!<
!-------------------------------------------------------------------------------
LOGICAL FUNCTION ISEVEN( num )
  INTEGER, INTENT(IN) :: num !<A block row
  ISEVEN = MOD(num,2) .EQ. 0
END FUNCTION ISEVEN

!-------------------------------------------------------------------------------
!>
!! Is given integer odd?
!<
!-------------------------------------------------------------------------------
LOGICAL FUNCTION ISODD( num )
  INTEGER, INTENT(IN) :: num !<A block row
  ISODD = (.NOT. ISEVEN(num))
END FUNCTION ISODD

!DEC$ IF DEFINED (MPI_OPT)
!-------------------------------------------------------------------------------
!>
!!Determines the global row number of a given local row number at a given level
!<
!-------------------------------------------------------------------------------
FUNCTION LR2GR( locrow, level ) RESULT ( globrow )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: locrow !<local row number of global row, at given level
  INTEGER, INTENT(IN) :: level !<Level at which locrow's position is given
  INTEGER :: globrow !<Returned: Row number in original input (level 1)

  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: i !<Loop variable

  !-----------------------------------------------------------------------------
  ! Invert the rule r=(r+1)/2 backward for all levels
  !-----------------------------------------------------------------------------
  IF (locrow .lt. 1) THEN
     globrow = -1
  ELSE
     globrow = locrow
     LevelLoop: DO i = level-1, 1, -1
        globrow = 2*globrow - 1
     END DO LevelLoop
  END IF

  RETURN
END FUNCTION LR2GR

!-------------------------------------------------------------------------------
!>
!!Determines the local row number of "globrow" global row number when globrow
!!participates at a given level; returns zero if this globrow does not operate
!!at the given level
!<
!-------------------------------------------------------------------------------
FUNCTION GR2LR( globrow, level ) RESULT ( locrow )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: globrow !<Row number in original input (level 1)
  INTEGER, INTENT(IN) :: level !<Level at which globrow's position is needed
  INTEGER :: locrow !<Returned: local row number of global row, at given level

  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: i !<Loop variable

  !-----------------------------------------------------------------------------
  ! Apply the rule r=(r+1)/2 for all odd-numbered rows
  ! Any time r becomes even, give up and return 0
  !-----------------------------------------------------------------------------
  locrow = globrow
  LevelLoop: DO i = 1, level-1, 1
    IF ( ISEVEN(locrow) ) THEN !Even-numbered rows don't go any further in level
      locrow = 0
      EXIT LevelLoop !Stop looping
    END IF
    locrow = (locrow+1) / 2
  END DO LevelLoop

  RETURN
END FUNCTION GR2LR

!-------------------------------------------------------------------------------
!>
!!Determines the rank of the task holding the given global row (at level 1)
!<
!-------------------------------------------------------------------------------
FUNCTION GR2Rank( globrow ) RESULT ( rank )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: globrow !<Row number in original input (level 1)
  INTEGER :: rank !<Returned: Rank holding the given global row

  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: m     !Average integral number of rows per processor
  INTEGER :: spill !Spill over beyond prefect loading of rows to processors
  INTEGER :: P     
  INTEGER :: N     !<Total number of rows in original input (level 1)

  P = numtasks
  N = nblk_rows

  IF ( (globrow .LT. 1) .OR. (globrow .GT. N) ) THEN
    rank = -1
  ELSE
    m = N / P
    spill = MOD( N, P )
    IF ( globrow .LE. (m+1)*spill ) THEN
      rank = (globrow-1)/(m+1)
    ELSE
      rank = (globrow-1 - (m+1)*spill)/m + spill
    ENDIF
  ENDIF

  RETURN
END FUNCTION GR2Rank

!-------------------------------------------------------------------------------
!>
!!Determines the rank of the task holding the given local row (at a given level)
!<
!-------------------------------------------------------------------------------
FUNCTION LR2Rank( locrow, level ) RESULT ( rank )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: locrow !<Row number at a given level
  INTEGER, INTENT(IN) :: level !<Local row number is at this level
  INTEGER :: rank          !<Returned: Rank holding the given global row

  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: globrow !Global row number of given locrow

  globrow = LR2GR( locrow, level )
  rank = GR2Rank( globrow)

  RETURN
END FUNCTION LR2Rank

!-------------------------------------------------------------------------------
!>
!!Determines the first block row handled by this rank at this level
!<
!-------------------------------------------------------------------------------
FUNCTION GetFirstBlockRow( rank, maxrow, level ) RESULT ( blkrow )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: rank  !<the rank (processor ID)
  INTEGER, INTENT(IN) :: maxrow !<Max block row #
  INTEGER, INTENT(IN) :: level !<Local row number is at this level
!>
!!Returned: beginning block row at this level for this rank
!!Returns 0 if no rows at this level for this rank
!<
  INTEGER :: blkrow
  INTEGER :: ns

  blkrow = 0      
  DO ns = 1, maxrow
     IF (LR2Rank(ns, level) .eq. rank) THEN
        blkrow = ns
        EXIT
     END IF
  END DO

END FUNCTION GetFirstBlockRow

!-------------------------------------------------------------------------------
!>
!!Determines the last block row handled by this rank at this level
!<
!-------------------------------------------------------------------------------
FUNCTION GetLastBlockRow( rank, maxrow, level ) RESULT ( blkrow )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: rank  !<the rank (processor ID)
  INTEGER, INTENT(IN) :: maxrow !<Max block row #
  INTEGER, INTENT(IN) :: level !<Local row number is at this level
!>
!!Returned: ending block row at this level for this rank
!!Returns 0 if no rows at this level for this rank
!<
  INTEGER :: blkrow
  INTEGER :: ns

  blkrow = 0      
  DO ns = maxrow, 1, -1
     IF (LR2Rank(ns, level) .eq. rank) THEN
        blkrow = ns
        EXIT
     END IF
  END DO

END FUNCTION GetLastBlockRow

!DEC$ ENDIF
!-------------------------------------------------------------------------------
!>
!! Main parallel block solver routine.
!! Returns sections of solution in brhs and factorized block rows (for multiple
!! calls with the same matrix but different rhs)
!<
!-------------------------------------------------------------------------------
      SUBROUTINE BCYCLIC_SOLVER (lblk, dblk, ublk, ipivot, brhs, mblock, ns)
      INTEGER, INTENT(in)   :: mblock !<linear size of each block
      INTEGER, INTENT(in)   :: ns     !<number of block rows
      INTEGER, TARGET       :: ipivot(:,ns0:) !<row section of pivot array for block factorization
      REAL(rprec), TARGET   :: lblk(:,:,ns0:) !<row section of lower block
      REAL(rprec), TARGET   :: dblk(:,:,ns0:) !<row section of diagonal block
      REAL(rprec), TARGET   :: ublk(:,:,ns0:) !<row section of upper block
      REAL(rprec), TARGET   :: brhs(:,ns0:)   !<row section of rhs
      INTEGER, POINTER                       :: iptr(:,:) !<ptr to pivot array
      INTEGER, ALLOCATABLE                   :: nsLevel(:) !<number of remaining blocks at each level
!>
!!store number of extra blocks at each level, needed for unhatted L,U matrices
!<
      REAL(rprec), POINTER, DIMENSION(:,:,:) :: lptr  !<pointer to subsection of lower block
      REAL(rprec), POINTER, DIMENSION(:,:,:) :: dptr  !<pointer to subsection of diag block
      REAL(rprec), POINTER, DIMENSION(:,:,:) :: uptr  !<pointer to subsection of upper block
      REAL(rprec), POINTER, DIMENSION(:,:)   :: bptr  !<pointer to subsection of rhs/solution
      INTEGER      :: nCycle, nsTemp, ierr, istat,  &
     &                maxLevels, nLevel, ns0_l
      LOGICAL      :: lStored

!DEC$ IF DEFINED (MPI_OPT)
!     NEED TO INITIALLY SYNCHRONIZE PROCESSORS BEFORE SENDRECV STUFF
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!DEC$ ENDIF
      CALL SYSTEM_CLOCK(timeon, countRate)
      timeon1 = timeon
      tot_time_comm = 0
      tot_time_comp = 0
      tot_time_omp  = 0

!     COMPUTE NUMBER OF CYCLIC REDUCTION LEVELS
      nCycle = 1
      maxLevels = 0
      DO WHILE (nCycle .lt. ns)
         nCycle = 2*nCycle
         maxLevels = maxLevels+1
      END DO

      maxLevels = MAX(1, maxLevels)
      lStored = ALLOCATED(OddStorage)
      IF (lStored .AND. (SIZE(OddStorage).NE.maxLevels)) THEN
         CALL ClearStorage
         lStored = .FALSE.
      END IF

      IF (.NOT. lStored) THEN
         added_memory = 0
         init_memory = 3._dp*KIND(dblk)*SIZE(dblk,1)*SIZE(dblk,2)*SIZE(dblk,3)
         ALLOCATE (OddStorage(maxLevels), stat=istat)
         IF (istat .ne. 0) STOP 'Allocation error'
         DO nLevel = 1, maxLevels
            OddStorage(nLevel)%NumElements = 0
         END DO
      END IF

      ALLOCATE (nsLevel(maxLevels), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error'

!
!     PERFORM CYCLIC-REDUCTION FACTORIZATION STEP
!
      nCycle = 1
      nFirst = ns0
      ns0_l= ns0
      nblk_rows = ns            !global variable used in LR2GR, etc...
      nsTemp = nblk_rows        !remaining block rows at n-th level
      nLevel = 1

      lptr => lblk
      dptr => dblk
      uptr => ublk
      iptr => ipivot
      bptr => brhs

!DEC$ IF DEFINED (MPI_OPT)
!      IF (rank .eq. 0) PRINT *,' ns   rank   FACTOR(s)   COMMUN(s)   TIMESTAMP'
!DEC$ ENDIF
      CALL SYSTEM_CLOCK(timeoff1)
      tot_time_comp = tot_time_comp + (timeoff1-timeon1)

      FACTOR_CYCLE: DO

      loc_time_comm = 0
      timeon1 = timeoff1

      nsLevel(nLevel) = nsTemp

!     FACTOR AND STORE BLOCKS AT THIS LEVEL
!     IF MULTIPLE RHSs, ONLY CALL THIS FIRST TIME
      IF (.not.lStored) CALL CYCLE_ONESTEP(nsTemp, nLevel, lptr, dptr, uptr, iptr)

!     STORE INFORMATION TO COMPUTE EFFECTIVE RHS TERMS
!     CALL THIS FOR MULTIPLE RHSs
      CALL CYCLE_RHS(nsTemp, nLevel, dptr, bptr, iptr)

      CALL SYSTEM_CLOCK(timeoff1)
      tot_time_comm = tot_time_comm + loc_time_comm
      tot_time_comp = tot_time_comp + (timeoff1-timeon1)-loc_time_comm
      timeon1 = timeoff1

      IF (nsTemp .eq. 1) THEN
!         PRINT *,'nsTemp == 1, nLevel=',nLevel
         EXIT
      END IF

      nsTemp = (nsTemp+1)/2
      nLevel = nLevel+1
      nCycle = 2*nCycle

!DEC$ IF DEFINED (MPI_OPT)
      nFirst = GetFirstBlockRow(rank,nsn,nLevel)
!     No more rows to process (all processors other than 0 kick out here)
      IF (nFirst .lt. 1) THEN
!         PRINT *,' nFirst: ', nFirst
         EXIT
      END IF
      ns0_l = LR2GR(nFirst, nLevel)
!DEC$ ELSE
      nFirst = nFirst/2+1
!DEC$ ENDIF

      lptr => lblk(:,:,ns0_l:nsn:nCycle)
      dptr => dblk(:,:,ns0_l:nsn:nCycle)
      uptr => ublk(:,:,ns0_l:nsn:nCycle)
      iptr => ipivot(:,ns0_l:nsn:nCycle)
      bptr => brhs(:,ns0_l:nsn:nCycle)

      CALL SYSTEM_CLOCK(timeoff1)
      tot_time_comp = tot_time_comp + (timeoff1-timeon1)
!!DEC$ IF DEFINED (MPI_OPT)
!      WRITE (6,200) 2*ns/nCycle,rank,REAL(timeoff1-timeon1)/countRate, &
!     &     REAL(comm_off-comm_on)/countRate,timeoff1
! 200  FORMAT(2i5,1x,1pe12.4,1pe12.4,i12)
!!DEC$ ENDIF

      IF (nLevel .gt. maxLevels) THEN
!        Processor 0 kicks out here
!         PRINT *,' nLevel: ',nLevel,' nsTemp: ',nsTemp
         EXIT
      END IF

      END DO FACTOR_CYCLE
!
!     COMPUTE FINAL HIGHEST (ODD) LEVEL SOLUTION USING THOMAS ALGORITHM ON RANK=0 PROC
!     ALL OTHER PROCESSORS HAVE ZERO BLOCKS LEFT TO PROCESS AT THIS LEVEL
!     (ACTUALLY, THOMAS REDUCES TO SOLVING A PURELY DIAGONAL Dx = b EQUATION AT THIS LEVEL
!      SINCE THERE IS ONLY 1 BLOCK ROW)
!
      timeon1 = timeoff1
!DEC$ IF DEFINED (MPI_OPT)
      IF (rank .eq. 0) THEN
!DEC$ ENDIF
      IF (SIZE(dptr,3) .ne. 1) STOP 'CYCLIC REDUCTION FAILED!'
      IF (.not.lStored) CALL THOMAS_FACTOR(lptr, dptr, uptr, iptr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     END OF FACTORIZATION. BEGIN BACK-SOLVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     FIRST BACK-SOLVE FOR ROW=1 X's AT HIGHEST LEVEL USING THOMAS SOLVE
!     THIS WILL ALWAYS BE DONE ON rank=0 PROCESSOR
!
      CALL THOMAS_SOLVE(lptr, dptr, uptr, bptr, iptr)

      CALL SYSTEM_CLOCK(timeoff1)
      tot_time_comp = tot_time_comp + (timeoff1-timeon1)
      timeon1 = timeoff1
!      PRINT *,' THOMAS_SOLVE(s) : ',REAL(timeoff1-timeon1)/countRate
!DEC$ IF DEFINED (MPI_OPT)
!      PRINT *,' ns   rank   SOLVES(s)   COMMUN(s)   TIMESTAMP '
      END IF
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      CALL SYSTEM_CLOCK(timeoff1)
      tot_time_comm = tot_time_comm + (timeoff1-timeon1)
!DEC$ ENDIF

!     BACK SOLVE FOR EVEN LEVEL X's
      SOLVE_CYCLE: DO
      timeon1 = timeoff1
      loc_time_comm = 0
      nCycle = nCycle/2
      nLevel = nLevel-1
      IF (nLevel .lt. 1) STOP 'nLevel < 1!'
!DEC$ IF DEFINED (MPI_OPT)
      nFirst = GetFirstBlockRow(rank, nsn, nLevel)
      IF (nFirst .gt. 0) ns0_l = LR2GR(nFirst, nLevel)
!DEC$ ENDIF

      nsTemp = nsLevel(nLevel)
      lptr => lblk(:,:,ns0_l:nsn:nCycle)
      uptr => ublk(:,:,ns0_l:nsn:nCycle)
      bptr => brhs(:,ns0_l:nsn:nCycle)
      CALL CYCLE_SOLVE(nsTemp, nLevel, lptr, uptr, bptr)
      CALL SYSTEM_CLOCK(timeoff1)
      tot_time_comp = tot_time_comp + (timeoff1-timeon1)-loc_time_comm
!DEC$ IF DEFINED (MPI_OPT)
      tot_time_comm = tot_time_comm + loc_time_comm
!      WRITE (6,300) ns/nCycle/2,rank,REAL(timeoff1-timeon1)/countRate, &
!     &     REAL(loc_time_comm)/countRate,timeoff1
! 300  FORMAT(2i5,1x,1pe12.4,1pe12.4,i12)
!DEC$ ENDIF

      IF (nCycle .eq. 1) EXIT

      END DO SOLVE_CYCLE

      DEALLOCATE (nsLevel)

      CALL SYSTEM_CLOCK(timeoff)

      END SUBROUTINE BCYCLIC_SOLVER

!-------------------------------------------------------------------------------
!>
!! Deallocates storage of Odd (unhatted) blocks
!! Should be called after last right-side is processed
!<
!-------------------------------------------------------------------------------
      SUBROUTINE ClearStorage
      IMPLICIT NONE
      INTEGER     :: nLevel, istat

      nLevel = SIZE(OddStorage)
      DO WHILE (nLevel .gt. 0)
         IF (OddStorage(nLevel)%NumElements .gt. 0) THEN
            DEALLOCATE(OddStorage(nLevel) % UMAT, stat=istat)
            DEALLOCATE(OddStorage(nLevel) % LMAT, stat=istat)
         END IF
         nLevel = nLevel-1
      END DO
      
      DEALLOCATE (OddStorage)

      END SUBROUTINE ClearStorage

!-------------------------------------------------------------------------------
!>
!! Performs matrix factorization and storage for one (forward) cyclic-reduction level (nLevel)
!<
!-------------------------------------------------------------------------------
      SUBROUTINE CYCLE_ONESTEP(nblk, nLevel, lblk, dblk, ublk, ipiv)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: nblk, nLevel
      INTEGER, TARGET, INTENT(out) :: ipiv(:,nfirst:)
      REAL(rprec), TARGET, INTENT(inout) :: lblk(:,:,nfirst:),           &
     &                                      dblk(:,:,nfirst:), ublk(:,:,nfirst:)
      INTEGER          :: ns, mblk, ierr, ierr_loc, nmin, kCount, istat
      INTEGER, POINTER :: ipivot(:)
      REAL(rprec), POINTER, DIMENSION(:,:)  :: dmat, umat, lmat, umat0, lmat0
      REAL(rprec), POINTER, DIMENSION(:,:)  :: lnext, unext, lprev, uprev
      REAL(rprec), ALLOCATABLE, TARGET, DIMENSION(:,:) :: lprevr, lnextr, uprevr, unextr
      REAL(rprec), POINTER, DIMENSION(:,:)  :: lprevs, lnexts, uprevs, unexts
      INTEGER NTHREADS, TID, OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM

      mblk   = SIZE(lblk,1)
      nlast  = nblk
!DEC$ IF DEFINED (MPI_OPT)
      mblk2  = mblk*mblk
      nlast  = GetLastBlockRow(rank, nsn, nLevel)
!     EXIT IF THERE ARE NO ROW ELEMENTS FOR THIS PROCESSOR AT THIS LEVEL
      IF (nlast .lt. 1) RETURN

      CALL SYSTEM_CLOCK(comm_on, commRate)

!     PREPARE RECV BUFFERS (FROM PREV/NEXT PROCESSORS) FOR ODD BOUNDARY ROWS
!     ONLY ODD ROWS RECEIVE

      ALLOCATE (lprevr(mblk,mblk), lnextr(mblk,mblk),                          &
                uprevr(mblk,mblk), unextr(mblk,mblk), stat=ierr)
      IF (ierr .ne. 0) STOP 'Allocation of bdy blocks failed'

      CALL INITCOMM(nLevel)

      IF (ISODD(nfirst) .and. rankprev.ge.0) THEN
         mpireqind = mpireqind + 1
         CALL MPI_Irecv(lprevr, mblk2, MPI_REAL8, rankprev, tag1, MPI_COMM_WORLD, &
                     mpireq(mpireqind), mpierr(mpireqind))
         mpireqind = mpireqind + 1
         CALL MPI_Irecv(uprevr, mblk2, MPI_REAL8, rankprev, tag1+1, MPI_COMM_WORLD, &
                     mpireq(mpireqind), mpierr(mpireqind))
      END IF

      IF (ISODD(nlast) .and. ranknext.ge.0) THEN
         mpireqind = mpireqind + 1
         CALL MPI_Irecv(lnextr, mblk2, MPI_REAL8, ranknext, tag2, MPI_COMM_WORLD, &
                     mpireq(mpireqind), mpierr(mpireqind))
         mpireqind = mpireqind + 1
         CALL MPI_Irecv(unextr, mblk2, MPI_REAL8, ranknext, tag2+1, MPI_COMM_WORLD, &
                     mpireq(mpireqind), mpierr(mpireqind))
      END IF

      CALL SYSTEM_CLOCK(comm_off)
      loc_time_comm = loc_time_comm+(comm_off-comm_on)
!DEC$ ENDIF
!
!     STORE ORIGINAL L/U ODD BLOCKS NEEDED FOR CALCULATION OF MULTIPLE SOURCE TERMS
!     (in Eq. III-3b)
     
      nmin = nfirst
      IF (ISEVEN(nmin)) nmin = nmin+1

      kCount = 1+(nlast-nmin)/2
      OddStorage(nLevel)%NumElements = kCount
 
      ALLOCATE (OddStorage(nLevel)%UMAT(mblk,mblk,kCount),              &
                OddStorage(nLevel)%LMAT(mblk,mblk,kCount), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in CYCLE_ONESTEP'

      added_memory = added_memory + 2._dp*kCount*mblk*mblk*             &
                                    KIND(OddStorage(nLevel) % UMAT)

      kCount = 0
      DO ns = nmin, nlast, 2
         kCount = kCount+1
         OddStorage(nLevel) % UMAT(:,:,kCount) = ublk(:,:,ns)
         OddStorage(nLevel) % LMAT(:,:,kCount) = lblk(:,:,ns)
      END DO

      IF (kCount .gt. OddStorage(nLevel)%NumElements)          &
         STOP 'kCount > NumElements IN CYCLE_ONESTEP'
!
!     COMPUTE LU FACTORIZATION OF EVEN DIAGONAL BLOCKS
!     (D_2k inverse in Eq. III-1)
!
      nmin = nfirst
      IF (ISODD(nmin)) nmin = nmin+1
	  ierr = 0

      CALL SYSTEM_CLOCK(omp_on)
!$OMP PARALLEL DO SCHEDULE(RUNTIME) &
!$OMP PRIVATE(ipivot,dmat,umat,lmat,ns,ierr_loc) &
!$OMP SHARED(mblk,nmin,nlast,ipiv,dblk,ublk,lblk,ierr)
      DO ns = nmin, nlast, 2
         ipivot => ipiv(:,ns)
         dmat   => dblk(:,:,ns)
         CALL DGETRF(mblk, mblk, dmat, mblk, ipivot, ierr_loc)
         IF (ierr_loc .ne. 0) ierr = ierr_loc
!
!     COMPUTE NEW L=-D[-1]L, U=-D[-1]U FOR EVEN BLOCKS
!     (Eq III-1, -L and -U hat)
!
         lmat => lblk(:,:,ns)
         CALL DGETRS('n',mblk,mblk,dmat,mblk,ipivot,lmat,mblk,ierr_loc)
         IF (ierr_loc .ne. 0) ierr = ierr_loc
         lmat = -lmat

         umat => ublk(:,:,ns)
         CALL DGETRS('n',mblk,mblk,dmat,mblk,ipivot,umat,mblk,ierr_loc)
         IF (ierr_loc .ne. 0) ierr = ierr_loc
         umat = -umat
      END DO 
!$OMP END PARALLEL DO
      CALL SYSTEM_CLOCK(omp_off)
      tot_time_omp = tot_time_omp + (omp_off-omp_on)

	  IF (ierr .ne. 0) GOTO 200

!DEC$ IF DEFINED (MPI_OPT)
      comm_on = omp_off

!     COMMUNICATE BOUNDARY BLOCKS TO next AND prev PROCESSORS
!     EVEN ROWS SEND DATA

      OUTER_TEST: IF (nlast.gt.0 .and. numtasks.gt.1) THEN

      IF (ISEVEN(nfirst) .and. rankprev.ge.0) THEN
         mpireqind = mpireqind + 1
         lprevs => lblk(:,:,nfirst)
         CALL MPI_Isend(lprevs, mblk2, MPI_REAL8, rankprev, tag2, MPI_COMM_WORLD, &
                        mpireq(mpireqind), mpierr(mpireqind))
         mpireqind = mpireqind + 1
         uprevs => ublk(:,:,nfirst)
         CALL MPI_Isend(uprevs, mblk2, MPI_REAL8, rankprev, tag2+1, MPI_COMM_WORLD, &
                        mpireq(mpireqind), mpierr(mpireqind))
      END IF

      IF (ISEVEN(nlast) .and. ranknext.ge.0) THEN
         mpireqind = mpireqind + 1
         lnexts => lblk(:,:,nlast)
         CALL MPI_Isend(lnexts, mblk2, MPI_REAL8, ranknext, tag1, MPI_COMM_WORLD, &
                        mpireq(mpireqind), mpierr(mpireqind))
         mpireqind = mpireqind + 1
         unexts => ublk(:,:,nlast)
         CALL MPI_Isend(unexts, mblk2, MPI_REAL8, ranknext, tag1+1, MPI_COMM_WORLD, &
                        mpireq(mpireqind), mpierr(mpireqind))
      END IF

      END IF OUTER_TEST

      CALL SYSTEM_CLOCK(comm_off)
      loc_time_comm = loc_time_comm+(comm_off-comm_on)
!DEC$ ENDIF

!     COMPUTE NEW ODD BLOCKS IN TERMS OF EVEN BLOCK FACTORS (Eq. III-1,
!     CALCULATED ABOVE), USING Eq. III-3a 
!

!     COMPUTE ODD HATTED MATRIX ELEMENTS EXCEPT AT BOUNDARIES (Eq. III-3a)
      CALL SYSTEM_CLOCK(omp_on)
      nmin = nfirst
      IF (ISEVEN(nmin)) nmin=nmin+1

!$OMP PARALLEL DO SCHEDULE(RUNTIME) &
!$OMP PRIVATE(dmat,umat,lmat,ns,lnext,unext,lprev,uprev,kCount,umat0,lmat0) &
!$OMP SHARED(nLevel,nmin,nfirst,nlast,nblk,dblk,ublk,lblk,OddStorage)
    DO ns = nmin, nlast, 2
!DEC$ IF DEFINED (MPI_OPT)
       IF (ns.eq.nfirst .or. ns.eq.nlast) CYCLE
!DEC$ ENDIF
       dmat => dblk(:,:,ns); umat => ublk(:,:,ns); lmat => lblk(:,:,ns)
       lnext => lblk(:,:,ns+1); unext => ublk(:,:,ns+1)
       lprev => lblk(:,:,ns-1); uprev => ublk(:,:,ns-1)

       kCount = 1+(ns-nmin)/2
       lmat0 => OddStorage(nLevel)%LMAT(:,:,kCount)
       umat0 => OddStorage(nLevel)%UMAT(:,:,kCount)

       CALL ComputeOddRowHats(nblk, ns, lmat, dmat, umat, lmat0, umat0,  &
                              lnext, unext, lprev, uprev)
    END DO
!$OMP END PARALLEL DO
      CALL SYSTEM_CLOCK(omp_off)
      tot_time_omp = tot_time_omp + (omp_off-omp_on)

!DEC$ IF DEFINED (MPI_OPT)
      CALL WAITFORALL

!     COMPUTE REMAINING (BOUNDARY) ODD HATTED MATRIX ELEMENTS
      CALL SYSTEM_CLOCK(omp_on)
!$OMP PARALLEL DO SCHEDULE(RUNTIME) &
!$OMP PRIVATE(dmat,umat,lmat,ns,lnext,unext,lprev,uprev,kCount,lmat0,umat0) &
!$OMP SHARED(nLevel,nmin,nfirst,nlast,nblk,dblk,ublk,lblk,lnextr,unextr,lprevr,uprevr)
    DO ns = nmin, nlast, 2
       IF (ns.ne.nfirst .and. ns.ne.nlast) CYCLE
       dmat => dblk(:,:,ns); umat => ublk(:,:,ns); lmat => lblk(:,:,ns)
       IF (ns .lt. nblk) THEN
          IF (ns .eq. nlast) THEN
             lnext => lnextr;  unext => unextr
          ELSE
             lnext => lblk(:,:,ns+1); unext => ublk(:,:,ns+1)
          END IF
       END IF
       IF (ns .gt. 1) THEN
          IF (ns .eq. nfirst) THEN
             lprev => lprevr;  uprev => uprevr
          ELSE
             lprev => lblk(:,:,ns-1); uprev => ublk(:,:,ns-1)
          END IF
       END IF

       kCount = 1+(ns-nmin)/2
       lmat0 => OddStorage(nLevel)%LMAT(:,:,kCount)
       umat0 => OddStorage(nLevel)%UMAT(:,:,kCount)

       CALL ComputeOddRowHats(nblk, ns, lmat, dmat, umat, lmat0, umat0,  &
                              lnext, unext, lprev, uprev)
    END DO
!$OMP END PARALLEL DO
      CALL SYSTEM_CLOCK(omp_off)
      tot_time_omp = tot_time_omp + (omp_off-omp_on)

      DEALLOCATE (lprevr, lnextr, uprevr, unextr, stat=ierr)
!DEC$ ENDIF

      RETURN

 200  CONTINUE
!          < 0:  if info = -i, the i-th argument had an illegal value
!          > 0:  if info = i, u(i,i) is exactly zero. the factorization
      WRITE (6, '(x,a,i4)') 'Error factoring matrix in CYCLE_ONESTEP: block = ', ns
      IF (ierr < 0) WRITE (6,'(i4, a)') ierr, 'th argument has illegal value'
      IF (ierr > 0) WRITE (6,'(i4, a)') ierr, 'th diagonal factor exactly zero'
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_ABORT(MPI_COMM_WORLD,ierr)
!DEC$ ELSE
      STOP
!DEC$ ENDIF
      END SUBROUTINE CYCLE_ONESTEP

!-------------------------------------------------------------------------------
!>
!!Compute the matrix-matrix multiplications in the forward solve at the current level
!!The local row (=ns) should be odd at this level. Implementation of Eq. III-3a
!!NOTE: USE dgemm TO AVOID STACK OVERFLOW WHICH CAN OCCUR WITH MATMUL
!<
!-------------------------------------------------------------------------------
      SUBROUTINE ComputeOddRowHats(nblk, ns, lmat, dmat, umat, lmat0, umat0,      &
                                   lnext, unext, lprev, uprev)
      IMPLICIT NONE
      INTEGER, INTENT(in)                         :: ns, nblk
      REAL(rprec), DIMENSION(:,:), INTENT(inout)  :: dmat, umat, lmat
      REAL(rprec), DIMENSION(:,:), INTENT(in)     :: lnext, unext, lprev, uprev,  &
                                                     lmat0, umat0
      INTEGER                                     :: mblk

      mblk = SIZE(dmat,1)

      IF (ns .lt. nblk) THEN
!!       dmat = dmat + MATMUL(umat0, lnext)
         CALL DGEMM('N','N',mblk,mblk,mblk,one,umat0,mblk,lnext,mblk,one,dmat,mblk)

!!       umat = MATMUL(umat0, unext)
         CALL DGEMM('N','N',mblk,mblk,mblk,one,umat0,mblk,unext,mblk,zero,umat,mblk)
      END IF

      IF (ns .gt. 1) THEN
!!       dmat = dmat + MATMUL(lmat0,uprev)
         CALL DGEMM('N','N',mblk,mblk,mblk,one,lmat0,mblk,uprev,mblk,one,dmat,mblk)

!!       lmat = MATMUL(lmat0, lprev)
         CALL DGEMM('N','N',mblk,mblk,mblk,one,lmat0,mblk,lprev,mblk,zero,lmat,mblk)
      END IF  
      
      END SUBROUTINE ComputeOddRowHats

!-------------------------------------------------------------------------------
!>
!! Performs storage of hatted rhs for one forward cyclic-reduction level (nLevel)
!<
!-------------------------------------------------------------------------------
      SUBROUTINE CYCLE_RHS(nblk, nLevel, dblk, brhs, ipiv)
      IMPLICIT NONE
      INTEGER, INTENT(in)                :: nblk, nLevel
      INTEGER, TARGET, INTENT(in)        :: ipiv(:,nfirst:)
      REAL(rprec), TARGET, INTENT(in)    :: dblk(:,:,nfirst:)
      REAL(rprec), TARGET, INTENT(inout) :: brhs(:,nfirst:)
      INTEGER                            :: mblk, ns, ierr, ierr_loc, nmin, kCount
      INTEGER, POINTER                   :: ipivot(:)
      REAL(rprec), POINTER               :: dmat(:,:), umat(:,:), lmat(:,:),     &
                                            xptr(:), bptr(:), bprev(:), bnext(:)
      REAL(rprec), ALLOCATABLE, TARGET, DIMENSION(:) :: bprevr, bnextr

      mblk   = SIZE(dblk,1)
      nLast = nblk
!DEC$ IF DEFINED (MPI_OPT)
      nlast  = GetLastBlockRow(rank, nsn, nLevel)
!     NO ROW ELEMENTS FOR THIS PROCESSOR AT THIS LEVEL
      IF (nlast .lt. 1) RETURN

      CALL SYSTEM_CLOCK(comm_on)

      ALLOCATE (bprevr(mblk), stat=ierr)
      ALLOCATE (bnextr(mblk), stat=ierr)

      CALL INITCOMM(nLevel)

      IF (ISODD(nfirst) .and. rankprev.ge.0) THEN
         mpireqind = mpireqind + 1
         CALL MPI_Irecv(bprevr, mblk, MPI_REAL8, rankprev, tag2, MPI_COMM_WORLD, &
                        mpireq(mpireqind), mpierr(mpireqind))
      END IF

      IF (ISODD(nlast) .and. ranknext.ge.0) THEN
         mpireqind = mpireqind + 1
         CALL MPI_Irecv(bnextr, mblk, MPI_REAL8, ranknext, tag1, MPI_COMM_WORLD, &
                        mpireq(mpireqind), mpierr(mpireqind))
      END IF

      CALL SYSTEM_CLOCK(comm_off)
      loc_time_comm = loc_time_comm+(comm_off-comm_on)
!DEC$ ENDIF
!
!     COMPUTE DBLK[-1]*BRHS FOR EVEN INDICES AND STORE IN BRHS(even) - Eq.III-1
!
      CALL SYSTEM_CLOCK(omp_on)
      nmin = nfirst
      IF (ISODD(nmin)) nmin = nmin+1
	  ierr_loc = 0

!$OMP PARALLEL DO SCHEDULE(RUNTIME) &
!$OMP PRIVATE(ipivot,bptr,dmat,ns,ierr) &
!$OMP SHARED(nmin,nLast,mblk,brhs,dblk,ipiv,ierr_loc)
      DO ns = nmin, nLast, 2
         ipivot => ipiv(:,ns);  dmat => dblk(:,:,ns)
         bptr => brhs(:,ns)
         CALL DGETRS('n', mblk, 1, dmat, mblk, ipivot, bptr, mblk, ierr_loc)
         IF (ierr_loc .ne. 0) ierr = ierr_loc
      END DO
!$OMP END PARALLEL DO
      CALL SYSTEM_CLOCK(omp_off)
      tot_time_omp = tot_time_omp + (omp_off-omp_on)

	  IF (ierr .ne. 0) GOTO 200

!DEC$ IF DEFINED (MPI_OPT)
!     COMMUNICATION LOOP: EVEN ENDPOINT VALUES SENT TO NEXT/PREV PROCESSORS
      comm_on = omp_off
      IF (ISEVEN(nfirst) .and. rankprev.ge.0) THEN
         xptr => brhs(:,nfirst)
         mpireqind = mpireqind + 1
         CALL MPI_Isend(xptr, mblk, MPI_REAL8, rankprev, tag1, MPI_COMM_WORLD, &
                       mpireq(mpireqind), mpierr(mpireqind))
      END IF

      IF (ISEVEN(nlast) .and. ranknext.ge.0) THEN
         xptr => brhs(:,nlast)
         mpireqind = mpireqind + 1
         CALL MPI_Isend(xptr, mblk, MPI_REAL8, ranknext, tag2, MPI_COMM_WORLD, &
                        mpireq(mpireqind), mpierr(mpireqind))
      END IF

      CALL SYSTEM_CLOCK(comm_off)
      loc_time_comm = loc_time_comm + (comm_off-comm_on)
!DEC$ ENDIF

!
!     COMPUTE ODD (HATTED) SOURCES (b-hats, Eq. III-3b) FOR INTERIOR ROWS
!
      CALL SYSTEM_CLOCK(omp_on)
      nmin = nfirst
      IF (ISEVEN(nmin)) nmin = nmin+1
      kCount = 0

!$OMP PARALLEL DO SCHEDULE(RUNTIME) &
!$OMP PRIVATE(bnext,bprev,bptr,lmat,umat,ns,kCount) &
!$OMP SHARED(nblk,mblk,nmin,nfirst,nlast,brhs,nLevel,OddStorage)
      DO ns = nmin, nlast, 2
!DEC$ IF DEFINED (MPI_OPT)
         IF (ns.eq.nfirst .or. ns.eq.nlast) CYCLE
!DEC$ ENDIF
         bptr => brhs(:,ns)
         kCount = 1+(ns-nmin)/2
         umat => OddStorage(nLevel)%UMAT(:,:,kCount)
         bnext => brhs(:,ns+1)
         lmat => OddStorage(nLevel)%LMAT(:,:,kCount)
         bprev => brhs(:,ns-1)
         CALL ComputeMVProds(nblk, ns, lmat, umat, bptr, bnext, bprev)
      END DO 
!$OMP END PARALLEL DO
      CALL SYSTEM_CLOCK(omp_off)
      tot_time_omp = tot_time_omp + (omp_off-omp_on)

!DEC$ IF DEFINED (MPI_OPT)
      CALL WAITFORALL

      CALL SYSTEM_CLOCK(omp_on)
!
!     FINALLY, COMPUTE b-hats AT BOUNDARY POINTS
!$OMP PARALLEL DO SCHEDULE(RUNTIME) &
!$OMP PRIVATE(bnext,bprev,bptr,lmat,umat,ns,kCount) &
!$OMP SHARED(nblk,nmin,nfirst,nlast,brhs,nLevel,bnextr,bprevr,OddStorage)
      DO ns = nmin, nlast, 2
         IF (ns.ne.nfirst .and. ns.ne.nlast) CYCLE
         bptr => brhs(:,ns)
         kCount = 1+(ns-nmin)/2
         IF (ns .lt. nblk) THEN
            umat => OddStorage(nLevel)%UMAT(:,:,kCount)
            IF (ns .eq. nlast) THEN
               bnext => bnextr(:)
            ELSE
               bnext => brhs(:,ns+1)
            END IF
         END IF
         IF (ns .gt. 1) THEN
            lmat => OddStorage(nLevel)%LMAT(:,:,kCount)
            IF (ns .eq. nfirst) THEN
               bprev => bprevr(:)
            ELSE
               bprev => brhs(:,ns-1)
            END IF
         END IF
         CALL ComputeMVProds(nblk, ns, lmat, umat, bptr, bnext, bprev)
      END DO 
!$OMP END PARALLEL DO
      CALL SYSTEM_CLOCK(omp_off)
      tot_time_omp = tot_time_omp + (omp_off-omp_on)
      DEALLOCATE (bprevr, bnextr)
!DEC$ ENDIF

      IF (kCount .gt. OddStorage(nLevel)%NumElements) THEN
         PRINT *,'nLevel: ',nLevel,' kCount: ',kCount,' NumEl: ',  OddStorage(nLevel)%NumElements
         STOP 'kCount > NumElements IN CYCLE_RHS'
      END IF

      RETURN

 200  CONTINUE

      PRINT *, 'Error in CYCLE_RHS: ierr = ', ierr
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_ABORT(MPI_COMM_WORLD,ierr)
!DEC$ ELSE
      STOP
!DEC$ ENDIF

      END SUBROUTINE CYCLE_RHS
      
!-------------------------------------------------------------------------------
!>
!!Compute matrix-vector multiplications in both forward/backward solves at the current level
!!In forward solve,  ns is odd and Eq.III-3b is solved for b-odd (bptr)
!!In backward solve, ns is even and Eq. III-1 (first line) is solved for x-even (bptr)
!<
!-------------------------------------------------------------------------------
      SUBROUTINE ComputeMVProds(nblk, ns, lmat, umat, bptr, bnext, bprev)
      IMPLICIT NONE
      INTEGER, INTENT(in)                         :: ns, nblk
      REAL(rprec), DIMENSION(:,:), INTENT(in)     :: lmat, umat
      REAL(rprec), DIMENSION(:),   INTENT(inout)  :: bptr
      REAL(rprec), DIMENSION(:),   INTENT(in)     :: bnext, bprev
      INTEGER                                     :: mblk
      REAL(rprec)                                 :: sgnone

      mblk = SIZE(bptr)
      sgnone = -one
      IF (ISEVEN(ns)) sgnone = one

      IF (ns .lt. nblk) THEN
!!       bptr = bptr - MATMUL(umat,bptr)
         CALL DGEMV('N',mblk,mblk,sgnone,umat,mblk,bnext,1,one,bptr,1)
      END IF

      IF (ns .gt. 1) THEN
!!       bptr = bptr - MATMUL(lmat,bptr)
         CALL DGEMV('N',mblk,mblk,sgnone,lmat,mblk,bprev,1,one,bptr,1)
      END IF

      END SUBROUTINE ComputeMVProds

!-------------------------------------------------------------------------------
!>
!!     Computes EVEN index solution from the computed (at previous,higher level)
!!     ODD index solutions at THIS level, using Eq.III-1.
!<
!-------------------------------------------------------------------------------
      SUBROUTINE CYCLE_SOLVE(nblk, nLevel, lblk, ublk, brhs)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: nblk, nLevel
      REAL(rprec), TARGET, INTENT(in)                :: lblk(:,:,nfirst:), ublk(:,:,nfirst:)
      REAL(rprec), TARGET, INTENT(inout)             :: brhs(:,nfirst:)
      INTEGER                                        :: mblk, ns, ierr, nmin
      REAL(rprec), POINTER                           :: umat(:,:), lmat(:,:), bprev(:), &
                                                        bnext(:), bptr(:)
      REAL(rprec), ALLOCATABLE, TARGET, DIMENSION(:) :: bprevr, bnextr

      mblk   = SIZE(lblk,1)
      nlast = nblk

!     NOTE AT THIS POINT, THE ODD BRHS VALUES HAVE BEEN REPLACED (AT THE HIGHEST CYCLE)
!     WITH THE SOLUTION VALUES (X), AT SUBSEQUENT (LOWER) CYCLES, THE
!     ODD VALUES ARE REPLACED BY THE EVEN SOLUTIONS AT THE NEXT HIGHEST CYCLE. THE EVEN 
!     BRHS VALUES WERE MULTIPLIED BY D[-1] AND STORED IN CYCLE_RHS
!
!DEC$ IF DEFINED (MPI_OPT)
      nlast  = GetLastBlockRow(rank, nsn, nLevel)
!     NO ROW ELEMENTS FOR THIS PROCESSOR AT THIS LEVEL
      IF (nlast .lt. 1) RETURN

      ALLOCATE (bprevr(mblk), stat=ierr)
      ALLOCATE (bnextr(mblk), stat=ierr)
      CALL DistributeScalar(mblk,nFirst,nLevel,brhs,bprevr,bnextr)
!DEC$ ENDIF
!
!     SOLVE FOR even INDEX VALUES IN TERMS OF (COMPUTED AT THIS POINT) odd INDEX VALUES
!     (Eq. III-1: note here, UMAT, LMAT are -U-hat, -L-hat in paper)

      CALL SYSTEM_CLOCK(omp_on)
      nmin = nfirst
      IF (ISODD(nmin)) nmin=nmin+1
!$OMP PARALLEL DO SCHEDULE(RUNTIME) &
!$OMP PRIVATE(bprev,bnext,bptr,umat,lmat,ns) &
!$OMP SHARED(nmin,nFirst,nLast,nblk,bnextr,bprevr,brhs,ublk,lblk)
      DO ns = nmin, nLast, 2
         bptr => brhs(:,ns)
         lmat => lblk(:,:,ns)
!DEC$ IF DEFINED (MPI_OPT)
         IF (ns .eq. nFirst) THEN
            bprev => bprevr(:)
         ELSE
!DEC$ ENDIF
            bprev => brhs(:,ns-1)
!DEC$ IF DEFINED (MPI_OPT)
         END IF
!DEC$ ENDIF
         IF (ns .lt. nblk) THEN
            umat => ublk(:,:,ns)
!DEC$ IF DEFINED (MPI_OPT)
            IF (ns .eq. nLast) THEN
               bnext => bnextr(:)
            ELSE
!DEC$ ENDIF
               bnext => brhs(:,ns+1)
!DEC$ IF DEFINED (MPI_OPT)
            END IF
!DEC$ ENDIF
         END IF
         CALL ComputeMVProds(nblk, ns, lmat, umat, bptr, bnext, bprev)
      END DO
!$OMP END PARALLEL DO
      CALL SYSTEM_CLOCK(omp_off)
      tot_time_omp = tot_time_omp + (omp_off-omp_on)
!DEC$ IF DEFINED (MPI_OPT)
      DEALLOCATE (bprevr, bnextr)
!DEC$ ENDIF

      RETURN

      END SUBROUTINE CYCLE_SOLVE


!-------------------------------------------------------------------------------
!>
!!     Uses the serial Thomas factorization algorithm to invert a (reduced)
!!     block tridiagonal system (at the final cyclic reduction level)
!<
!-------------------------------------------------------------------------------
      SUBROUTINE THOMAS_FACTOR(lblk, dblk, ublk, ipiv)
      IMPLICIT NONE
      INTEGER, TARGET, INTENT(out) :: ipiv(:,:)
      REAL(rprec), TARGET, INTENT(inout) :: lblk(:,:,:), dblk(:,:,:), ublk(:,:,:)
      INTEGER :: nblk, mblk, ns, ns1, ierr
      INTEGER, POINTER :: ipivot(:)
      REAL(rprec), POINTER    :: dmat(:,:), umat(:,:), lmat(:,:)
      REAL(rprec), ALLOCATABLE :: temp(:,:)
!
!     USE THOMAS ALGORITHM TO COMPUTE TRIDIAGONAL BLOCKS
!

      nblk   = SIZE(lblk,3)
      mblk   = SIZE(lblk,1)
      ipiv = 0

      BLOCKS: DO ns = nblk, 1, -1
!
!     Compute (and save) qblk(ns) = dblk(ns)[-1] * lblk
!
         dmat => dblk(:,:,ns)

         ipivot => ipiv(:,ns)
         CALL DGETRF (mblk, mblk, dmat, mblk, ipivot, ierr)
         IF (ierr .ne. 0) GOTO 200
         IF (ns .eq. 1) EXIT

         lmat => lblk(:,:,ns)

         CALL DGETRS('n', mblk, mblk, dmat, mblk, ipivot,               &
     &                lmat, mblk, ierr)
         IF (ierr .ne. 0) GOTO 305
!
!      Update effective diagonal matrix. Use dgemm: faster AND doesn't overflow stack
!
         ns1 = ns-1 
         umat => ublk(:,:,ns1)
         dmat => dblk(:,:,ns1)
         CALL DGEMM('N','N',mblk,mblk,mblk,-one,umat,mblk,              &
     &              lmat, mblk, one, dmat, mblk)

      END DO BLOCKS

!
!     COMPUTE TRANSPOSES HERE, SINCE REPEATEDLY CALLING MATMUL OPERATION
!     X*At IS FASTER THAN A*X DUE TO UNIT STRIDE
!

      ALLOCATE (temp(mblk,mblk), stat=ierr)
      IF (ierr .ne. 0) STOP 'Allocation error in thomas_factor!'

      DO ns = 1, nblk
         IF (ns .ne. nblk) THEN
            temp = TRANSPOSE(ublk(:,:,ns))
            ublk(:,:,ns) = temp
         END IF
         IF (ns .ne. 1) THEN
            temp = TRANSPOSE(lblk(:,:,ns))
            lblk(:,:,ns) = temp
         END IF
      END DO

      DEALLOCATE (temp)

      RETURN

!  error returns. ------------------------------------------------------

  200 CONTINUE
!          < 0:  if info = -i, the i-th argument had an illegal value
!          > 0:  if info = i, u(i,i) is exactly zero. the factorization
      WRITE (6, '(x,a,i4)') 'Error factoring matrix in THOMAS_FACTOR: block = '   &
     &                        , ns
      IF (ierr < 0)WRITE (6,'(i4, a)') ierr, 'th argument has illegal value'
      IF (ierr > 0)WRITE (6,'(i4, a)') ierr, 'th diagonal factor exactly zero'
      STOP
  301 CONTINUE
!      WRITE (6, '(a,i8)') ' BLK3D:   error in opening file:  ',
!     1   'RECL = ', irecl
  302 CONTINUE
      WRITE (6, '(a)') ' BLK3D:   error in I/O routine WRDISK'
  303 CONTINUE
      WRITE (6, '(a)') ' BLK3D:   error in I/O routine RDDISK'
      ierr = -2
  305 CONTINUE
      WRITE (6, '(2/a,i10,2/)') ' BLK3D-FACTOR:   error detected:   ier =',     &
     &   ierr
      STOP
 
      END SUBROUTINE THOMAS_FACTOR

      
      SUBROUTINE THOMAS_SOLVE(lblk, dblk, ublk, brhs, ipiv)
      IMPLICIT NONE
      INTEGER, TARGET, INTENT(in) :: ipiv(:,:)
      REAL(rprec), TARGET, INTENT(in) :: lblk(:,:,:), dblk(:,:,:), ublk(:,:,:)
      REAL(rprec), TARGET, INTENT(inout)  :: brhs(:,:)
      INTEGER :: nblk, mblk, ns, ns1, ierr
      INTEGER, POINTER :: ipivot(:)
      REAL(rprec), POINTER :: dmat(:,:), umat(:,:), lmat(:,:), bptr(:), x1(:)
!
!     USE THOMAS ALGORITHM TO COMPUTE TRIDIAGONAL BLOCKS
!
      nblk   = SIZE(lblk,3)
      mblk   = SIZE(lblk,1)
      IF (SIZE(brhs,2) .ne. nblk .or. SIZE(brhs,1) .ne. mblk) STOP 'brhs wrong dims'

      BLOCKS: DO ns = nblk, 1, -1

         bptr => brhs(:,ns)
         dmat => dblk(:,:,ns)
         ipivot => ipiv(:,ns);   
         CALL DGETRS('n', mblk, 1, dmat, mblk, ipivot, bptr, mblk, ierr)
         IF (ierr .ne. 0) GOTO 305

         IF (ns .eq. 1) EXIT

!
!        NOTE: IN BLK3D_FACTOR, BP1 AND BM1 WERE TRANSPOSED (AND STORED)
!        TO MAKE FIRST INDEX FASTEST VARYING IN THE FOLLOWING MATMUL OPS
!
         ns1 = ns-1
         umat => ublk(:,:,ns1)
         bptr => brhs(:,ns1)
         x1 => brhs(:,ns)
!         bptr = bptr - MATMUL(x1,umat)  !USE THIS FORM IF TRANSPOSED bp1
!         bptr = bptr - MATMUL(umat,x1)  !UNTRANSPOSED FORM
         CALL DGEMV('T',mblk,mblk,-one,umat,mblk,x1,1,one,bptr,1)

      END DO BLOCKS
!
!  forward (back-substitution) solution sweep for block-rows k = 2 to nblocks
!  now, source contains the solution vector
!
      DO ns = 2, nblk
         ns1 = ns-1
         bptr => brhs(:,ns)
         x1 => brhs(:,ns1)
         lmat => lblk(:,:,ns)
!         bptr = bptr - MATMUL(x1,lmat)  !USE THIS FORM IF TRANSPOSED qblk
!         bptr = bptr - MATMUL(lmat,y1)  !UNTRANSPOSED FORM
         CALL DGEMV('T',mblk,mblk,-one,lmat,mblk,x1,1,one,bptr,1)

      END DO

      RETURN

!  error returns. ------------------------------------------------------
  305 CONTINUE
      WRITE (6, '(2/a,i10,2/)') ' THOMAS_SOLVE:   error detected:   ierr =', ierr
      STOP

      END SUBROUTINE THOMAS_SOLVE

      
!-------------------------------------------------------------------------------
!>
!!     Sends scalar values at row boundaries to previous/next processors
!<
!-------------------------------------------------------------------------------
!DEC$ IF DEFINED (MPI_OPT)
      SUBROUTINE DistributeScalar(mblk,nstart,nLevel,xsol,bprevr, bnextr)
      IMPLICIT NONE
      INTEGER, INTENT(in)              :: nLevel, mblk, nstart
      REAL(rprec), TARGET              :: xsol(:,nstart:)
      REAL(rprec), INTENT(out)         :: bprevr(:), bnextr(:)
      REAL(rprec), POINTER             :: xptr(:)
      
      CALL SYSTEM_CLOCK(comm_on)

      CALL INITCOMM(nLevel)

      IF (rankprev .ge. 0) THEN
         mpireqind = mpireqind + 1
         CALL MPI_Irecv(bprevr, mblk, MPI_REAL8, rankprev, tag2, MPI_COMM_WORLD, &
                     mpireq(mpireqind), mpierr(mpireqind))
         mpireqind = mpireqind + 1
         xptr => xsol(:,nFirst)
         CALL MPI_Isend(xptr, mblk, MPI_REAL8, rankprev, tag1, MPI_COMM_WORLD, &
                        mpireq(mpireqind), mpierr(mpireqind))
      END IF

      IF (ranknext .ge. 0) THEN
         mpireqind = mpireqind + 1
         CALL MPI_Irecv(bnextr, mblk, MPI_REAL8, ranknext, tag1, MPI_COMM_WORLD, &
                     mpireq(mpireqind), mpierr(mpireqind))
         mpireqind = mpireqind + 1
         xptr => xsol(:,nLast)
         CALL MPI_Isend(xptr, mblk, MPI_REAL8, ranknext, tag2, MPI_COMM_WORLD, &
                        mpireq(mpireqind), mpierr(mpireqind))
      END IF

      CALL SYSTEM_CLOCK(comm_off)
      loc_time_comm = loc_time_comm+(comm_off-comm_on)

      CALL WAITFORALL

      END SUBROUTINE DistributeScalar

!-------------------------------------------------------------------------------
!>
!!     Waits for all pending MPI requests to complete
!<
!-------------------------------------------------------------------------------
      SUBROUTINE WAITFORALL
      INTEGER :: mpiwaiterr !<MPI_Waitall/any error
      
      mpireqcnt = mpireqind
      IF (mpireqcnt .gt. 6) STOP 'MPIREQCNT > 6!'
      IF (mpireqcnt .lt. 0) RETURN
         
      CALL SYSTEM_CLOCK(comm_on)
      CALL MPI_Waitall(mpireqcnt, mpireq, mpistats, mpiwaiterr)
      CALL SYSTEM_CLOCK(comm_off)
      loc_time_comm = loc_time_comm + (comm_off-comm_on)

      END SUBROUTINE WAITFORALL

!-------------------------------------------------------------------------------
!>
!!    Initializes MPI request array and finds the prev/next rank(processor)
!<
!-------------------------------------------------------------------------------
      SUBROUTINE INITCOMM(nLevel)
      INTEGER, INTENT(in)   :: nLevel

      DO mpireqind = 1, 6
         mpireq(mpireqind) = MPI_REQUEST_NULL
      END DO

      mpireqind = 0

      rankprev = LR2Rank(nFirst-1,nLevel)
      ranknext = LR2Rank(nLast+1,nLevel)
      
      END SUBROUTINE INITCOMM
!DEC$ ENDIF

      SUBROUTINE CHECK_SOLVER (lblk, dblk, ublk, xsol, brhs, mblock, ns, istep)
      INTEGER, INTENT(in)     :: mblock, ns, istep
      REAL(rprec), TARGET     :: xsol(:,ns0:)
      REAL(rprec), INTENT(in) :: brhs(:,ns0:)
      REAL(rprec), INTENT(in) :: lblk(:,:,ns0:), dblk(:,:,ns0:), ublk(:,:,ns0:)
      REAL(rprec), POINTER    :: xptr(:), xptr2(:,:), x_global(:,:), x_extend(:,:)
      REAL(rprec), TARGET, ALLOCATABLE :: bprevr(:), bnextr(:) 
      INTEGER                 :: msize, nIndex, ierr, istat, k, nLevel, n1, n2, globrow
      REAL(rprec)             :: rms_error, t1
!
!     checks solution Ax = b for each processor chunk
!     writes out full solution vector on fort.66
!
!     ALL ARRAYS (EXCEPT x_global) ARE ON THE BLOCK SUBSECTION FOR THIS PROCESSOR

!     lblk, dblk, ublk:  original lower, diagonal, upper blocks
!     brhs: original right-hand-side
!     xsol: solution, stored in chunks
!     x_global: solution for all blocks, ONLY known on processor 0
!     x_extend: extended xsol to one block row below, above this processor
!
!DEC$ IF DEFINED (MPI_OPT)
      nLevel = 1
      IF (rank .ne. 0) THEN
         nFirst = GetFirstBlockRow(rank, ns, nLevel)
         nLast  = GetLastBlockRow(rank, ns, nLevel)

         xptr2 => xsol(:,nFirst:nLast)
         msize = nLast-nFirst+1
         CALL MPI_Send(xptr2, mblock*msize, MPI_REAL8, 0, rank,    &
                       MPI_COMM_WORLD, ierr)

      ELSE 
         ALLOCATE(x_global(mblock,ns), stat=istat)
         IF (istat .ne. 0) STOP 'Error allocating x_global!'

         k = 0
         nFirst = GetFirstBlockRow(k, ns, nLevel)
         nLast  = GetLastBlockRow(k, ns, nLevel)
         IF (nFirst .ne. 1) PRINT *,' nFirst != 1 on proc 0'

         x_global(:,nFirst:nLast)=xsol(:,nFirst:nLast)
         msize = nLast-nFirst+1
         nIndex = msize+1
         DO k = 1,numtasks-1
             xptr2  => x_global(:,nIndex:)
             CALL MPI_RECV(xptr2, mblock*msize, MPI_REAL8, k,k,    &
                           MPI_COMM_WORLD, mpistats, ierr)
             n1 = GetFirstBlockRow(k, ns, nLevel)
             n2 = GetLastBlockRow(k, ns, nLevel)
             msize = n2-n1+1
             nIndex = nIndex+msize
         END DO
!DEC$ ELSE
       nFirst = 1
       nLast  = ns
       x_global=>xsol
	   x_extend=>xsol
!DEC$ ENDIF
         DO istat=1,ns
            DO k = 1,mblock
               WRITE (33*(1+istep), 100) istat, k, x_global(k,istat)
            END DO
         END DO
!DEC$ IF DEFINED (MPI_OPT)
         DEALLOCATE(x_global)
      END IF

      n1 = MAX(1,nFirst-1)
	  n2 = MIN(ns,nLast+1)
      ALLOCATE(x_extend(mblock,n1:n2), stat=istat)
	  x_extend(:,nFirst:nLast) = xsol(:,nFirst:nLast)
!      ALLOCATE (bprevr(mblock), bnextr(mblock), stat=istat)
!      CALL DistributeScalar(mblock,ns0,nLevel,xsol,bprevr,bnextr)
      CALL DistributeScalar(mblock,ns0,nLevel,xsol,x_extend(:,n1),x_extend(:,n2))
!DEC$ ENDIF

      rms_error = 0

      DO globrow = nFirst,nLast
         DO k = 1,mblock
		    t1 = 0
			IF (globrow .gt. 1) THEN
               t1 = t1 + SUM(lblk(k,:,globrow)*x_extend(:,globrow-1))
			END IF

			t1 = t1 + SUM(dblk(k,:,globrow)*x_extend(:,globrow))

			IF (globrow .lt. ns) THEN
               t1 = t1 + SUM(ublk(k,:,globrow)*x_extend(:,globrow+1))
			END IF

            rms_error = rms_error + (brhs(k,globrow)-t1)**2

         END DO
      END DO


!DEC$ IF .NOT.DEFINED (MPI_OPT)
      PRINT *,' rms_error: ', SQRT(rms_error/(ns*mblock)), &
     &     ' CYCLIC SOLVE(s): ',REAL(timeoff-timeon)/countRate
!DEC$ ELSE
      DEALLOCATE (bprevr, bnextr, stat=istat)
      tot_time = timeoff-timeon
      PRINT 150,'PID: ',rank,' rms_error: ', SQRT(rms_error/(ns*mblock)), &
          ' TOTAL TIME(s): ',REAL(tot_time)/countRate,                    &
          ' COMP (s): ',REAL(tot_time_comp)/countRate,                    &
          ' COMM (s): ',REAL(tot_time_comm)/countRate,                    &
          ' OMP  (s): ',REAL(tot_time_omp)/countRate
      CALL MPI_REDUCE (tot_time_comm,sum_comm,1,MPI_INTEGER,MPI_SUM,0,    &
                       MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE (tot_time,max_tot,1,MPI_2INTEGER,MPI_MAXLOC,        &
                       0,MPI_COMM_WORLD,ierr)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      IF (rank .eq. 0) THEN
	  PRINT 160,' Max  tot  time (s): ',REAL(max_tot)/countRate
	  PRINT 160,' Mean comm time (s): ',REAL(sum_comm)/(numtasks*countRate)
      END IF
!DEC$ ENDIF

 100  FORMAT(2i5,1pe12.4)
 150  FORMAT(a,i5,a,1pe10.3,4(a,1pe9.2))
 160  FORMAT(a,1pe10.3)

      END SUBROUTINE CHECK_SOLVER

      END MODULE CYCLIC_RED


      PROGRAM CYCLICR
      USE stel_kinds
      USE CYCLIC_RED
      IMPLICIT NONE
      INTEGER, PARAMETER        :: NBLOCKS=256
      INTEGER, PARAMETER        :: MBLOCK_SIZE=273
      INTEGER                   :: ns, mblock, istat, k, nmin, irhs
      INTEGER, POINTER          :: ipivot(:,:)
      REAL(rprec), ALLOCATABLE,  TARGET, DIMENSION(:,:,:) ::                     &
                   lblk, dblk, ublk, lblk1, dblk1, ublk1
      REAL(rprec), ALLOCATABLE, TARGET  :: brhs(:,:), gc(:,:)
      REAL(rprec)  :: t1, rms_error
!******************************************
!
!  mpi setup calls: 
!
!DEC$ IF DEFINED (MPI_OPT)
      INTEGER          :: ierr
      CALL MPI_INIT( ierr )
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )       !mpi stuff
      CALL MPI_COMM_size( MPI_COMM_WORLD, numtasks, ierr )   !mpi stuff
!DEC$ ENDIF      

      mblock = MBLOCK_SIZE
      ns = NBLOCKS
!******************************************
!
!I.   DO SERIAL THOMAS FOR TIMING COMPARISON (SKIP IF BLOCKS ARE TOO LARGE)

      IF (MBLOCK_SIZE .le. 300) CALL SERIAL_THOMAS_TEST(mblock, ns)

!II.  DO CYCLIC REDUCTION
!     IMPORTANT: THE BLOCKS AND DATA ARE ASSUMED TO BE DISTRIBUTED CONSECUTIVELY
!                ON NODES AS DESCRIBED IN THE PAPER ....

!*****************************************
      nsn = ns
      ns0 = 1
!DEC$ IF DEFINED (MPI_OPT)
      nblk_rows = ns
      ns0 = GetFirstBlockRow(rank, ns, 1)
      nsn = GetLastBlockRow(rank, ns, 1)
      IF (rank .eq. 0) THEN
         WRITE (6, 100) mblock, ns
 100     FORMAT (2x,'BLOCK SIZE: ',i5,' BLOCK ROWS: ',i5)
      END IF
!DEC$ ENDIF

!     NOTE: IN PRODUCTION RUN, SKIP THIS SINCE ipivot, lblk, dblk, ublk, brhs ARE
!           KNOWN INPUTS
	  ALLOCATE (lblk(mblock,mblock,ns0:nsn), dblk(mblock,mblock,ns0:nsn),   &
                ublk(mblock,mblock,ns0:nsn), brhs(mblock,ns0:nsn),          &
				ipivot(mblock, ns0:nsn), stat=istat)
	  IF (istat .ne. 0) STOP 'Allocation error!'
      CALL INITIALIZE_BLOCKS (lblk, dblk, ublk, brhs, mblock)

!     STORE lblk1, ETC ONLY FOR BACK-SOLVER CHECK
      ALLOCATE (lblk1(mblock,mblock,ns0:nsn), dblk1(mblock,mblock,ns0:nsn),   &
                ublk1(mblock,mblock,ns0:nsn), gc(mblock,ns0:nsn),             &
                stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error!'

      lblk1 = lblk; dblk1 = dblk; ublk1 = ublk
      gc = brhs
!*****************************************

!     SIMULATE MULTIPLE RHS
      MULTIPLE_RHS: DO irhs = 1, 2

!     SOLVE USING CYCLIC REDUCTION
!     SOLUTION OVERWRITE brhs ON EACH PROCESSOR
!DEC$ IF DEFINED (MPI_OPT)
      IF (rank .eq. 0) PRINT '(/,a,i1,/)',' SOLUTION RHS #',irhs
!DEC$ ELSE
      PRINT '(/,a,i1,/)',' SOLUTION RHS #',irhs
!DEC$ ENDIF
      CALL BCYCLIC_SOLVER (lblk, dblk, ublk, ipivot, brhs, mblock, ns)

!     CHECK SOLUTION
      CALL CHECK_SOLVER (lblk1, dblk1, ublk1, brhs, gc, mblock, ns, irhs)

!     IMPORTANT: DO NOT RESET FACTORED BLOCKS
!     RESET brhs (SOLUTION) TO INITIAL VALUE
      brhs(:,ns0:nsn)=gc(:,ns0:nsn)

      END DO MULTIPLE_RHS

!     CLEAN UP ALLOCATED ARRAYS
      CALL ClearStorage
      DEALLOCATE (ipivot, lblk, dblk, ublk, brhs, stat=istat)
      DEALLOCATE (lblk1, ublk1, dblk1, gc, stat=istat)

!     SHUT DOWN MPI 
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      PRINT 200,'PID: ',rank,' INIT MEMORY (Mb): ', init_memory/1E6,  &
                             ' ADDED MEMORY (Mb): ', added_memory/1E6
 200  FORMAT(a,i5,2(a,f10.3))
      CALL MPI_FINALIZE(ierr)
!DEC$ ENDIF

      END PROGRAM CYCLICR

      
	  SUBROUTINE SERIAL_THOMAS_TEST (mblock, ns)
      USE stel_kinds
      USE CYCLIC_RED
      IMPLICIT NONE
      INTEGER, INTENT(in)            :: ns, mblock
      INTEGER                        :: istat, k, nmin, irhs
      INTEGER, TARGET, ALLOCATABLE   :: ipivot(:,:)
      REAL(rprec), ALLOCATABLE, TARGET, DIMENSION(:,:,:) ::                     &
                   lblk, dblk, ublk, lblk1, dblk1, ublk1
      REAL(rprec), ALLOCATABLE, TARGET   :: brhs(:,:), gc(:,:)
	  INTEGER, ALLOCATABLE           :: seed(:)
	  REAL(rprec)                    :: rms_error, t1
!******************************************
!DEC$ IF DEFINED (MPI_OPT)
      IF (rank .ne. 0) RETURN
!DEC$ ENDIF      

      CALL SYSTEM_CLOCK(timeon, countRate)

!******************************************
!      OPEN(UNIT=33, FILE='TEST_CYCLICR2.bin', FORM='UNFORMATTED',         &
!     &     STATUS='OLD', IOSTAT=istat)

!      IF (istat .ne. 0) STOP 'FILE NOT READ'
!      READ (33, iostat=istat) ns, mblock
!      IF (istat .ne. 0) STOP 'Read error!'
!      IF (ns .ne. nblocks .or. mblock .ne. mblock_size) STOP 'Wrong dimensions!'
!      nmin = MIN(ns, 101)

	  ns0 = 1
	  nsn = ns

      ALLOCATE (ipivot(mblock, ns), stat=istat)
      ALLOCATE (lblk(mblock,mblock,ns), dblk(mblock,mblock,ns),   &
                ublk(mblock,mblock,ns), lblk1(mblock,mblock,ns),  &
                ublk1(mblock,mblock,ns), dblk1(mblock,mblock,ns), &
                brhs(mblock,ns), gc(mblock,ns), stat = istat)

!     RANDOM ELEMENTS

      CALL RANDOM_SEED(size=k)
      ALLOCATE(SEED(k), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error for SEED!'

      SEED = 12345
      CALL RANDOM_SEED(PUT=SEED)

      CALL RANDOM_NUMBER(lblk)
      lblk(:,:,ns0:nsn)=-1 + 0.25*lblk(:,:,ns0:nsn)

      CALL RANDOM_NUMBER(dblk)
      dblk(:,:,ns0:nsn)=2 + .5*dblk(:,:,ns0:nsn)

      CALL RANDOM_NUMBER(ublk)
      ublk(:,:,ns0:nsn)=-1 + 0.25*ublk(:,:,ns0:nsn)

      CALL RANDOM_NUMBER(brhs)

      DEALLOCATE (seed)

      CALL SYSTEM_CLOCK(timeoff, countRate)
      PRINT *,' TIME TO CREATE RANDOM DATA: ', REAL(timeoff-timeon)/countRate,' s'

!I.   WITHOUT CYCLIC REDUCTION
      gc = brhs
      lblk1 = lblk; dblk1 = dblk; ublk1 = ublk

      CALL SYSTEM_CLOCK(timeon)
      CALL THOMAS_FACTOR(lblk, dblk, ublk, ipivot)
      CALL THOMAS_SOLVE(lblk, dblk, ublk, brhs, ipivot)
      CALL SYSTEM_CLOCK(timeoff)
      DO istat = 1,ns
      DO k = 1,mblock
      WRITE (33, 100) istat, k, brhs(k,istat)
      END DO
      END DO

 100  FORMAT(2i5,1pe12.4)

!Ia.  BACKSOLVE: CHECK SOLUTION
      rms_error = 0
      DO istat = 1,ns
         DO k = 1,mblock
            IF (istat .eq. 1) THEN
               t1 = SUM(dblk1(k,:,istat)*brhs(:,istat) + ublk1(k,:,istat)*brhs(:,istat+1))
            ELSE IF (istat .eq. ns) THEN
               t1 = SUM(lblk1(k,:,istat)*brhs(:,istat-1) + dblk1(k,:,istat)*brhs(:,istat))
            ELSE
               t1 = SUM(lblk1(k,:,istat)*brhs(:,istat-1) + dblk1(k,:,istat)*brhs(:,istat)  &
     &         + ublk1(k,:,istat)*brhs(:,istat+1))
            END IF
            rms_error = rms_error + (gc(k,istat)-t1)**2
         END DO
      END DO

      PRINT *,' rms_error: ', SQRT(rms_error/(ns*mblock)),    &
     &    ' THOMAS SOLVE(s): ',REAL(timeoff-timeon)/countRate

      DEALLOCATE (ipivot, lblk, dblk, ublk, lblk1, dblk1, ublk1, brhs, stat=istat)
	  DEALLOCATE (gc)

	  END SUBROUTINE SERIAL_THOMAS_TEST

      SUBROUTINE INITIALIZE_BLOCKS(lblk, dblk, ublk, brhs, mblock)
      USE stel_kinds
!DEC$ IF DEFINED (MPI_OPT)
      USE CYCLIC_RED, ONLY: rank, ns0, nsn
!DEC$ ELSE
      USE CYCLIC_RED, ONLY: ns0, nsn
!DEC$ ENDIF
      IMPLICIT NONE
      INTEGER, INTENT(in)      :: mblock
      REAL(rprec), DIMENSION(mblock,mblock,ns0:nsn), INTENT(out) :: lblk, dblk, ublk
      REAL(rprec), INTENT(out) :: brhs(mblock,ns0:nsn)
      INTEGER                  :: istat, k
	  INTEGER, ALLOCATABLE   :: seed(:)
      
      CALL RANDOM_SEED(size=k)
      ALLOCATE(SEED(k), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error for SEED!'

      SEED = 12345
!DEC$ IF DEFINED (MPI_OPT)
      SEED = SEED + 10*rank
!DEC$ ENDIF
      CALL RANDOM_SEED(PUT=SEED)

      CALL RANDOM_NUMBER(lblk)
      lblk=-1 + 0.5*(.5-lblk)

      CALL RANDOM_NUMBER(dblk)
      dblk= 2 + (.5-dblk)

      CALL RANDOM_NUMBER(ublk)
      ublk=-1 + 0.5*(.5-ublk)

      CALL RANDOM_NUMBER(brhs)
	  brhs = 2*(.5-brhs)

      DEALLOCATE (seed)

	  END SUBROUTINE INITIALIZE_BLOCKS

