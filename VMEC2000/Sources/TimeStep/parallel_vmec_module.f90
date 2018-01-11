!------------------------------------------------
MODULE parallel_vmec_module
  !------------------------------------------------

  USE stel_kinds
  USE mpi_inc
  IMPLICIT NONE
  INTEGER :: TOFU
  INTEGER :: grank=0, gnranks=1
  INTEGER :: vrank=0, vnranks=1
  INTEGER :: rank=0, nranks=1
  INTEGER :: last_ns = -1
  INTEGER :: par_ns, mynsnum, par_nzeta, par_ntheta3
  INTEGER :: par_ntmax, par_ntor, par_mpol1, par_nznt, par_nuv, par_nuv3
  INTEGER :: blocksize, ntmaxblocksize
  INTEGER :: tlglob, trglob, t1lglob, t1rglob, t2lglob, t2rglob
  INTEGER :: nuvmin, nuvmax, nuv3min, nuv3max 
  INTEGER :: MPI_ERR
#if defined(MPI_OPT)
  INTEGER :: MPI_STAT(MPI_STATUS_SIZE)
#endif
  INTEGER :: RUNVMEC_COMM_WORLD
  INTEGER :: NS_COMM
  INTEGER :: VAC_COMM
  INTEGER :: TWODCOMM, px, py
  INTEGER :: NS_RESLTN=0
  INTEGER :: num_grids
  INTEGER, ALLOCATABLE, DIMENSION(:) :: grid_procs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: grid_size
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: grid_time
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: f3d_time
  INTEGER, ALLOCATABLE, DIMENSION(:) :: f3d_num

  REAL(dp), ALLOCATABLE, DIMENSION(:) :: vgrid_time

  INTEGER, ALLOCATABLE, DIMENSION(:) :: tlglob_arr, trglob_arr 
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nuv3min_arr,nuv3max_arr 
  INTEGER, ALLOCATABLE, DIMENSION(:) :: trow, brow
  INTEGER, ALLOCATABLE, DIMENSION(:) :: lcol, rcol

  INTEGER, ALLOCATABLE, DIMENSION (:) :: blkrcounts, blkdisp
  INTEGER, ALLOCATABLE, DIMENSION (:) :: ntblkrcounts, ntblkdisp
  INTEGER, ALLOCATABLE, DIMENSION (:) :: nsrcounts, nsdisp
  INTEGER, PRIVATE :: mmax, nmax, tmax

  LOGICAL :: PARVMEC=.FALSE.
  LOGICAL :: LV3FITCALL=.FALSE.
  LOGICAL :: LIFFREEB=.FALSE.
  LOGICAL :: LPRECOND=.FALSE.

  LOGICAL :: lactive=.FALSE.
  LOGICAL :: vlactive=.FALSE.

  CHARACTER*100 :: envvar
  CHARACTER*100 :: envval

  REAL(dp) :: bcyclic_comp_time
  REAL(dp) :: bcyclic_comm_time
  REAL(dp) :: waitall_time
  REAL(dp) :: dgemm_time
  REAL(dp) :: dgetrf_time
  REAL(dp) :: dgetrs_time
  REAL(dp) :: dgemv_time
  REAL(dp) :: ForwardSolveLoop_time
  REAL(dp), DIMENSION(12) :: t 

  REAL(dp) :: allgather_time=0
  REAL(dp) :: allreduce_time=0
  REAL(dp) :: broadcast_time=0
  REAL(dp) :: sendrecv_time=0
  REAL(dp) :: scatter_time=0

!
!     OVERLOADED FUNCTIONS
!
  INTERFACE CopyLastNType
     MODULE PROCEDURE copy4lastntype, copy1lastntype,    &
                      copym4lastntype, copym1lastntype
  END INTERFACE

CONTAINS

  !--------------------------------------------------------------------------
  ! Read in environment variables 
  !--------------------------------------------------------------------------
  SUBROUTINE MyEnvVariables
    
    PARVMEC=.TRUE.
    envvar='PARVMEC'
    CALL GETENV(envvar,envval)
    IF (envval.EQ.'FALSE'.OR.envval.EQ.'false' &
      .OR.envval.EQ.'F'.OR.envval.EQ.'f') THEN
      PARVMEC=.FALSE.
    END IF

    LV3FITCALL=.FALSE.
    envvar='LV3FITCALL'
    CALL GETENV(envvar,envval)
    IF (envval.EQ.'TRUE'.OR.envval.EQ.'true' &
      .OR.envval.EQ.'T'.OR.envval.EQ.'t') THEN
      LV3FITCALL=.TRUE.
    END IF

  END SUBROUTINE  MyEnvVariables
  !--------------------------------------------------------------------------


  !--------------------------------------------------------------------------
  ! Declarations of all timers to be used in the parallel implementation
  !--------------------------------------------------------------------------
  SUBROUTINE InitializeParallel
    
#if defined(MPI_OPT)
    CALL MPI_Init(MPI_ERR)
    CALL MPI_Comm_rank(MPI_COMM_WORLD,grank,MPI_ERR) 
    CALL MPI_Comm_size(MPI_COMM_WORLD,gnranks,MPI_ERR) 
#endif

!    IF (PARVMEC.EQV..FALSE..AND.gnranks.EQ.1) THEN
!      PARVMEC=.FALSE.
!      LV3FITCALL=.FALSE.
!    END IF

  END SUBROUTINE InitializeParallel
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE InitRunVmec(INCOMM, LVALUE)
    INTEGER, INTENT(IN)  :: INCOMM
    LOGICAL, INTENT(IN)  :: LVALUE

    NS_RESLTN = 0
    LIFFREEB = LVALUE
#if defined(MPI_OPT)
    CALL MPI_Comm_dup(INCOMM,RUNVMEC_COMM_WORLD,MPI_ERR)
    CALL MPI_Comm_rank(RUNVMEC_COMM_WORLD,grank,MPI_ERR) 
    CALL MPI_Comm_size(RUNVMEC_COMM_WORLD,gnranks,MPI_ERR) 
#endif
  END SUBROUTINE InitRunVmec
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE InitSurfaceComm(ns, nzeta, ntheta3, ntmax, ntor, mpol1)
    
    INTEGER, INTENT(IN) :: ns, nzeta, ntheta3
    INTEGER, INTENT(IN) :: ntmax, ntor, mpol1
    INTEGER :: i
    LOGICAL :: FIRSTPASS = .TRUE.

      par_ns=ns
      par_nzeta=nzeta
      par_ntheta3=ntheta3
      par_ntmax=ntmax
      par_ntor=ntor
      par_mpol1=mpol1
      par_nznt=par_nzeta*par_ntheta3
      blocksize= (par_mpol1+1)*(par_ntor+1)
      ntmaxblocksize=3*ntmax*blocksize
      mmax=par_mpol1+1; nmax=par_ntor+1; tmax=3*par_ntmax
      NS_RESLTN=NS_RESLTN+1


    IF (LV3FITCALL) THEN
      IF(last_ns.NE.par_ns) THEN
        CALL SetSurfaceCommunicator
        CALL SetSurfacePartitions
        CALL SetSurfacePartitionArrays
        last_ns = par_ns
      END IF

    ELSE

      CALL SetSurfaceCommunicator
      CALL SetSurfacePartitions
      CALL SetSurfacePartitionArrays
      last_ns = par_ns

      IF (lactive) THEN 
        IF (grank.EQ.0) THEN 
          IF (FIRSTPASS) THEN 
            CALL SetOutputFile(rank,nranks,'parvmecinfo')
            WRITE(TOFU,*)"============================================================" 
            WRITE(TOFU,*) 'PARVMEC = ',PARVMEC
            WRITE(TOFU,*) "> available processor count:", gnranks
            WRITE(TOFU,*) '> global rank:', grank
            WRITE(TOFU,*) '> nzeta:      ', par_nzeta
            WRITE(TOFU,*) '> ntheta3:    ', par_ntheta3
            WRITE(TOFU,*) '> ntor:       ', par_ntor
            WRITE(TOFU,*) '> mpol1:      ', par_mpol1
            WRITE(TOFU,*) '> ntmax:      ', par_ntmax
            WRITE(TOFU,*) '> blocksize:  ',(par_mpol1+1)*(par_ntor+1)
            WRITE(TOFU,*)"============================================================" 
            WRITE(TOFU,*)
            CALL FLUSH(TOFU)
            FIRSTPASS = .FALSE.
          END IF
          WRITE(TOFU,*) ">>> grid resolution:   ",par_ns
          WRITE(TOFU,*) ">>> active processors: ",nranks
          WRITE(TOFU,*) ">>> xc/gc size:        ", par_ns*(par_ntor+1)*(par_mpol1+1)*3*ntmax
          WRITE(TOFU,*)"------------------------------------------------------------" 
          WRITE(TOFU,*)
          CALL FLUSH(TOFU)
        END IF
      END IF
    END IF
    50 FORMAT(3(A,I4))
  END SUBROUTINE InitSurfaceComm
  !------------------------------------------------

  !------------------------------------------------
  ! Setting up the communicator for parallel surface
  ! computations
  !------------------------------------------------
  SUBROUTINE SetSurfaceCommunicator
    
    INTEGER :: num_active, color

    num_active = MIN(gnranks,par_ns/2)

    IF(grank.LT.num_active) THEN
      color=1
    ELSE
      color=0
    END IF
#if defined(MPI_OPT)
    MPI_ERR = 0
    CALL MPI_Comm_split(RUNVMEC_COMM_WORLD,color,grank,NS_COMM,MPI_ERR)

    IF (color .eq. 1) THEN
      CALL MPI_Comm_size(NS_COMM,nranks,MPI_ERR)
      IF (nranks.NE.num_active) THEN
        STOP 'num_active != nranks in InitSurfaceCommunicator!'
      END IF
      lactive=.TRUE.
      CALL MPI_Comm_rank(NS_COMM,rank,MPI_ERR)
    ELSE
      nranks = 1
      rank   = 0
      lactive=.FALSE.
    END IF
#endif
  END SUBROUTINE SetSurfaceCommunicator
  !------------------------------------------------

  !------------------------------------------------
  ! Setting up the partitions for parallel surface
  ! computations
  !------------------------------------------------
  SUBROUTINE SetSurfacePartitions
    
    INTEGER :: mypart
#if defined(MPI_OPT)
    IF(par_ns.LT.nranks) THEN
      IF(grank.EQ.0) print *,"NS is less than NRANKS. Aborting!"
      CALL STOPMPI(5645)
    END IF

    mypart=par_ns/nranks
    IF (rank.LT.MOD(par_ns,nranks)) mypart=mypart+1
    IF (MOD(par_ns,nranks).NE.0) THEN
      IF (rank.LT.MOD(par_ns,nranks)) THEN
        tlglob=rank*mypart
      ELSE IF (rank.GE.MOD(par_ns,nranks)) THEN 
        tlglob=MOD(par_ns,nranks)*(mypart+1)+(rank-MOD(par_ns,nranks))*mypart
      END IF
    ELSE
      tlglob=rank*mypart
    END IF

    tlglob=tlglob+1
    trglob=tlglob+mypart-1

    t1lglob=tlglob-1; IF (rank.EQ.0) t1lglob=1
    t1rglob=trglob+1; IF (rank.EQ.nranks-1) t1rglob=par_ns

    IF(mypart.LT.2) THEN
      CALL MPI_Barrier(NS_COMM,MPI_ERR)
      WRITE(TOFU,*) '***********************************************************' 
      WRITE(TOFU,*) '* This version is not yet tested for mynsnum <= 2. Aborting!'
      WRITE(TOFU,*) '***********************************************************' 
      IF (rank.EQ.0) THEN
        WRITE(*,*)
        WRITE(*,*) '***********************************************************' 
        WRITE(*,*) '* This version is not yet tested for mynsnum <= 2. Aborting!'
        WRITE(*,*) '***********************************************************' 
        WRITE(*,*)
      END IF
      CALL MPI_Abort(NS_COMM,MPI_ERR)
    END IF
#endif
  END SUBROUTINE SetSurfacePartitions
  !------------------------------------------------

  !------------------------------------------------
  ! Setting up the partition arrays for parallel surface
  ! computations
  !------------------------------------------------
  SUBROUTINE SetSurfacePartitionArrays
    
    INTEGER, ALLOCATABLE, DIMENSION(:) :: localpart
    INTEGER :: i, smallpart, largepart
    INTEGER, SAVE :: lastns

    smallpart=par_ns/nranks; largepart = smallpart
    IF(MOD(par_ns,nranks).GT.0) largepart = largepart+1

    ALLOCATE (localpart(nranks))
    DO i = 0, nranks-1

      localpart(i+1) = smallpart

      IF(MOD(par_ns,nranks).GT.0) THEN
        IF(i.LT.MOD(par_ns,nranks)) localpart(i+1) = localpart(i+1) + 1
      END IF

    END DO

    IF(ALLOCATED(tlglob_arr)) DEALLOCATE(tlglob_arr)
    IF(ALLOCATED(trglob_arr)) DEALLOCATE(trglob_arr)
    ALLOCATE (tlglob_arr(nranks), trglob_arr(nranks))

    tlglob_arr(1)=1
    DO i = 2, nranks
      tlglob_arr(i) = tlglob_arr(i-1) + localpart(i-1)
    END DO
    DO i = 1, nranks
      trglob_arr(i) = tlglob_arr(i) + localpart(i) - 1
    END DO

    DEALLOCATE (localpart)

    CALL ComputeNSAllGatherParameters(nranks)
    CALL ComputeBlockAllGatherParameters(nranks)
    CALL ComputeNTmaxBlockAllGatherParameters(nranks)

  END SUBROUTINE SetSurfacePartitionArrays
  !------------------------------------------------

  !------------------------------------------------
  ! Compute AllGather vector variant parameters for 
  ! blocksized movements.
  !------------------------------------------------
  SUBROUTINE ComputeNTmaxBlockAllGatherParameters(activeranks)
    
    INTEGER :: activeranks
    INTEGER :: i
    IF(.NOT.ALLOCATED(ntblkrcounts)) ALLOCATE(ntblkrcounts(activeranks))
    IF(.NOT.ALLOCATED(ntblkdisp)) ALLOCATE(ntblkdisp(activeranks))
    DO i=1,activeranks
      ntblkrcounts(i)=(trglob_arr(i)-tlglob_arr(i)+1)*ntmaxblocksize
    END DO
    ntblkdisp(1)=0
    DO i=2,activeranks
      ntblkdisp(i)=ntblkdisp(i-1)+ntblkrcounts(i-1)
    END DO
  END SUBROUTINE ComputeNTmaxBlockAllGatherParameters
  !------------------------------------------------

  !------------------------------------------------
  ! Compute AllGather vector variant parameters for 
  ! blocksized movements.
  !------------------------------------------------
  SUBROUTINE ComputeBlockAllGatherParameters(activeranks)
    
    INTEGER :: activeranks
    INTEGER :: i
    IF(.NOT.ALLOCATED(blkrcounts)) ALLOCATE(blkrcounts(activeranks))
    IF(.NOT.ALLOCATED(blkdisp)) ALLOCATE(blkdisp(activeranks))
    DO i=1,activeranks
      blkrcounts(i)=(trglob_arr(i)-tlglob_arr(i)+1)*blocksize
    END DO
    blkdisp(1)=0
    DO i=2,activeranks
      blkdisp(i)=blkdisp(i-1)+blkrcounts(i-1)
    END DO
  END SUBROUTINE ComputeBlockAllGatherParameters
  !------------------------------------------------

  !------------------------------------------------
  ! Compute AllGather vector variant parameters for 
  ! blocksized movements.
  !------------------------------------------------
  SUBROUTINE ComputeNSAllGatherParameters(activeranks)
    
    INTEGER :: activeranks
    INTEGER :: i
    IF(.NOT.ALLOCATED(nsrcounts)) ALLOCATE(nsrcounts(activeranks))
    IF(.NOT.ALLOCATED(nsdisp)) ALLOCATE(nsdisp(activeranks))
    DO i=1,activeranks
      nsrcounts(i)=(trglob_arr(i)-tlglob_arr(i)+1)
    END DO
    nsdisp(1)=0
    DO i=2,activeranks
      nsdisp(i)=nsdisp(i-1)+nsrcounts(i-1)
    END DO
  END SUBROUTINE ComputeNSAllGatherParameters
  !------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE FinalizeSurfaceComm(INCOMM)
    
    INTEGER, INTENT(INOUT)  :: INCOMM
#if defined(MPI_OPT)
    CALL MPI_Comm_free(INCOMM,MPI_ERR)
    INCOMM=0; lactive = .false.
    IF(ALLOCATED(ntblkrcounts)) DEALLOCATE(ntblkrcounts)
    IF(ALLOCATED(ntblkdisp)) DEALLOCATE(ntblkdisp)
    IF(ALLOCATED(blkrcounts)) DEALLOCATE(blkrcounts)
    IF(ALLOCATED(blkdisp)) DEALLOCATE(blkdisp)
    IF(ALLOCATED(nsrcounts)) DEALLOCATE(nsrcounts)
    IF(ALLOCATED(nsdisp)) DEALLOCATE(nsdisp)
#endif
  END SUBROUTINE FinalizeSurfaceComm
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE FinalizeRunVmec(INCOMM)
    
    INTEGER, INTENT(INOUT)  :: INCOMM
#if defined(MPI_OPT)
    CALL MPI_Comm_free(INCOMM,MPI_ERR)
    IF(LIFFREEB) CALL MPI_Comm_free(VAC_COMM,MPI_ERR)
    INCOMM=0; VAC_COMM = 0
    rank = 0; par_ns = 0; nranks = 1 !SAL
    grank = 0; gnranks = 1; vrank = 0; vnranks = 1; last_ns = -1; !SAL
    NS_RESLTN=0
#endif
  END SUBROUTINE FinalizeRunVmec
  !--------------------------------------------------------------------------

  !------------------------------------------------
  ! Setting up the communicator for parallel vacuum
  ! computations
  ! nuv = nzeta*ntheta1
  ! nuv3 = nzeta*ntheta3
  !------------------------------------------------
  SUBROUTINE SetVacuumCommunicator(nuv, nuv3, mgs)
    
    INTEGER, INTENT(IN) :: nuv, nuv3, mgs
    INTEGER :: num_active, color, mypart

    par_nuv3 = nuv3
    par_nuv = nuv 
#if defined(MPI_OPT)
    num_active = MIN(gnranks,par_nuv3)
!    num_active = MIN(num_active,mgs/2) !SKS: temporary cutoff

    IF(grank.LT.num_active) THEN
      color=1
    ELSE
      color=0
    END IF
    CALL MPI_Comm_split(RUNVMEC_COMM_WORLD,color,grank,VAC_COMM,MPI_ERR)
    IF (color .eq. 1) THEN
      CALL MPI_Comm_rank(VAC_COMM,vrank,MPI_ERR) 
      CALL MPI_Comm_size(VAC_COMM,vnranks,MPI_ERR) 
      CALL SetVacuumPartitions(nuv3,nuv3min, nuv3max)
      CALL Setnuv3PartitionArrays
      vlactive=.TRUE.
    ELSE
      vnranks = 1
      vrank   = 0
      vlactive=.FALSE.
    ENDIF 
#else
    nuv3min = 1;  nuv3max = nuv3
#endif
  END SUBROUTINE SetVacuumCommunicator
  !------------------------------------------------


  !------------------------------------------------
  ! Setting up the partitions for parallel vacuum
  ! computations
  !------------------------------------------------
  SUBROUTINE SetVacuumPartitions(num,left,right)
    
    INTEGER, INTENT(IN) :: num
    INTEGER, INTENT(INOUT) :: left, right
    INTEGER :: mypart

    IF(num.LT.vnranks) THEN
      IF(grank.EQ.0) print *,"NUM is less than VNRANKS. Aborting!"
      CALL STOPMPI(456745)
    END IF

    mypart=num/vnranks
    IF (vrank.LT.MOD(num,vnranks)) mypart=mypart+1
    IF (MOD(num,vnranks).NE.0) THEN
      IF (vrank.LT.MOD(num,vnranks)) THEN
        left=vrank*mypart
      ELSE IF (vrank.GE.MOD(num,vnranks)) THEN 
        left=MOD(num,vnranks)*(mypart+1)+(vrank-MOD(num,vnranks))*mypart
      END IF
    ELSE
      left=vrank*mypart
    END IF

    left=left+1
    right=left+mypart-1

  END SUBROUTINE SetVacuumPartitions
  !------------------------------------------------

  !------------------------------------------------
  ! Setting up the partition arrays for parallel vacuum
  ! computations
  !------------------------------------------------
  SUBROUTINE Setnuv3PartitionArrays
    
    INTEGER, ALLOCATABLE, DIMENSION(:) :: localpart
    INTEGER :: i, smallpart, largepart

    IF (.NOT.ALLOCATED(nuv3min_arr)) THEN
      ALLOCATE(nuv3min_arr(vnranks),nuv3max_arr(vnranks))
    END IF

    smallpart=par_nuv3/vnranks; largepart = smallpart
    IF(MOD(par_nuv3,vnranks).GT.0) largepart = largepart+1

    ALLOCATE (localpart(vnranks))
    DO i = 0, vnranks-1

      localpart(i+1) = smallpart

      IF(MOD(par_nuv3,vnranks).GT.0) THEN
        IF(i.LT.MOD(par_nuv3,vnranks)) localpart(i+1) = localpart(i+1)+1
      END IF

    END DO

    nuv3min_arr(1)=1
    DO i = 2, vnranks
      nuv3min_arr(i) = nuv3min_arr(i-1) + localpart(i-1)
    END DO
    DO i = 1, vnranks
      nuv3max_arr(i) = nuv3min_arr(i) + localpart(i) - 1
    END DO

    DEALLOCATE (localpart)

  END SUBROUTINE Setnuv3PartitionArrays
  !------------------------------------------------

  !------------------------------------------------
  SUBROUTINE FinalizeParallel
    
    INTEGER :: istat

    envvar='LPRECOND'
    CALL GETENV(envvar,envval)
    IF (envval.EQ.'TRUE') THEN
      LPRECOND=.TRUE.
    ELSE IF (envval.EQ.'FALSE') THEN
      LPRECOND=.FALSE.
    END IF

    IF (grank.EQ.0) THEN 
      WRITE(*,*)
      WRITE(*,'(1x,a,i4)') 'NO. OF PROCS:  ',gnranks
      WRITE(*,100)         'PARVMEC     :  ',PARVMEC
      WRITE(*,100)         'LPRECOND    :  ',LPRECOND
      WRITE(*,100)         'LV3FITCALL  :  ',LV3FITCALL
    END IF
 100  FORMAT(1x,a,l4)

#if defined(MPI_OPT)
    CALL MPI_Finalize(istat)
#endif
  END SUBROUTINE FinalizeParallel
  !------------------------------------------------

  !------------------------------------------------
  ! Print debugging output to compare aray values 
  !------------------------------------------------
  SUBROUTINE PrintNSArray (arr, left, right, fileno, ifstop, string)
    
    INTEGER, INTENT(IN) :: fileno, left, right
    LOGICAL, INTENT(IN) :: ifstop
    REAL(dp), DIMENSION (par_ns+1), INTENT(IN) :: arr
    CHARACTER*(*), INTENT(IN) :: string

    INTEGER :: i, j, k, l

#if defined(_WIN32)
    RETURN
#endif

    DO i=left, right
      WRITE(fileno+rank,50) string, i, arr(i)
      CALL FLUSH(fileno+rank)
    END DO
    WRITE(fileno+rank,*) 
    CALL FLUSH(fileno+rank)
    IF(ifstop) STOP 'STOPPED CODE'

    50   FORMAT(A6,1I6,1P,E20.6)

  END SUBROUTINE PrintNSArray
  !------------------------------------------------

  !------------------------------------------------
  SUBROUTINE STOPMPI (code)
    
    INTEGER, INTENT(IN) :: code

#if defined(MPI_OPT)
    CALL MPI_Barrier(MPI_COMM_WORLD,MPI_ERR)
#endif
    WRITE(*,*) 'Stopping program with code:', code
    STOP
  END SUBROUTINE STOPMPI
  !------------------------------------------------

  !------------------------------------------------
  SUBROUTINE Parallel2Serial4X (inarr,outarr)
    
    REAL(dp), INTENT(IN) :: inarr(par_ns*(par_ntor+1)*(par_mpol1+1)*3*(par_ntmax))
    REAL(dp), INTENT(OUT) :: outarr(par_ns*(par_ntor+1)*(par_mpol1+1)*3*(par_ntmax))
    INTEGER :: i, j, k, l, lk, lks, lkp

    lks=0
    DO l=1, 3*par_ntmax
      DO k=0, par_mpol1
        DO j=0, par_ntor
          DO i=1, par_ns
            lks = lks+1
            lkp = j + (par_ntor+1)*k + (par_ntor+1)*(par_mpol1+1)*(i-1) + &
              (par_ntor+1)*(par_mpol1+1)*par_ns*(l-1)+1
            outarr(lks) = inarr(lkp)
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE Parallel2Serial4X
  !------------------------------------------------

  !------------------------------------------------
  SUBROUTINE Serial2Parallel4X (inarr,outarr)
    
    REAL(dp), INTENT(IN) :: inarr(par_ns*(par_ntor+1)*(par_mpol1+1)*3*(par_ntmax))
    REAL(dp), INTENT(OUT) :: outarr(par_ns*(par_ntor+1)*(par_mpol1+1)*3*(par_ntmax))
    INTEGER :: i, j, k, l, lk, lks, lkp

    lks=0
    DO l=1, 3*par_ntmax
      DO k=0, par_mpol1
        DO j=0, par_ntor
          DO i=1, par_ns
            lks = lks+1
            lkp = j + (par_ntor+1)*k + (par_ntor+1)*(par_mpol1+1)*(i-1) + &
              (par_ntor+1)*(par_mpol1+1)*par_ns*(l-1)+1
            outarr(lkp) = inarr(lks)
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE Serial2Parallel4X
  !------------------------------------------------

  !------------------------------------------------
  SUBROUTINE Parallel2Serial2X(inarr, outarr)
    
    REAL(dp), DIMENSION(par_nzeta,par_ntheta3,par_ns), INTENT(IN) :: inarr
    REAL(dp), DIMENSION(par_ns*par_nzeta*par_ntheta3+1), INTENT(OUT) :: outarr
    INTEGER :: i, j, k, l

    l=0
    DO k = 1, par_ntheta3
      DO j = 1, par_nzeta
        DO i = 1, par_ns
          l = l+1
          outarr(l)=inarr(j,k,i)
        END DO
      END DO
    END DO

  END SUBROUTINE Parallel2Serial2X 
  !------------------------------------------------

  !------------------------------------------------
  SUBROUTINE RPrintOutLinearArray(arr, left, right, flag, fileno)
    
    REAL(dp), DIMENSION(:) :: arr
    INTEGER, INTENT(IN) :: fileno, left, right
    LOGICAL, INTENT(IN) :: flag !flag: TRUE for parallel, FALSE for serial
    INTEGER :: i, j, k, l, lk

    REAL(dp), ALLOCATABLE, DIMENSION(:)  :: tmp
    ALLOCATE (tmp(ntmaxblocksize*par_ns))

    CALL tolastntype(arr,tmp)
    lk=0
    DO l=1, 3*par_ntmax
      DO k=0, par_mpol1
        DO j=0, par_ntor
          DO i=1, par_ns
            lk=lk+1
            IF(flag) THEN
              lk = j + nmax*(k + mmax*((i-1) + par_ns*(l-1)))+1
            END IF
            IF (left.LE.i.AND.i.LE.right) THEN
              WRITE(fileno+rank,50) i, j, k, l, tmp(lk)
              CALL FLUSH(fileno+rank)
            END IF
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(tmp)

 50 FORMAT(4I6,1P,E24.14)

  END SUBROUTINE RPrintOutLinearArray
  !------------------------------------------------


  !------------------------------------------------
  SUBROUTINE PrintOutLinearArray(arr, left, right, flag, fileno)
    
    REAL(dp), DIMENSION(ntmaxblocksize*par_ns) :: arr
    INTEGER, INTENT(IN) :: fileno, left, right
    LOGICAL, INTENT(IN) :: flag !flag: TRUE for parallel, FALSE for serial
    INTEGER :: i, j, k, l, lk

    lk=0
    DO l=1, 3*par_ntmax
      DO k=0, par_mpol1
        DO j=0, par_ntor
          DO i=1, par_ns
            lk=lk+1
            IF(flag) THEN
              lk = j + nmax*(k + mmax*((i-1) + par_ns*(l-1)))+1
            END IF
            IF (left.LE.i.AND.i.LE.right) THEN
              WRITE(fileno+rank,50) i, j, k, l, arr(lk)
              CALL FLUSH(fileno+rank)
            END IF
          END DO
        END DO
      END DO
    END DO

    50   FORMAT(4I6,1P,E24.14)

  END SUBROUTINE PrintOutLinearArray
  !------------------------------------------------

  !------------------------------------------------
  ! Prints out [nsmin, nsmax] part of a (nzeta, ntheta,ns) array
  !------------------------------------------------
  SUBROUTINE PrintParallelIJSPArray (arr, left, right, fileno, ifstop, string)
    
    INTEGER, INTENT(IN) :: fileno, left, right
    LOGICAL, INTENT(IN) :: ifstop
    REAL(dp), DIMENSION (par_nzeta,par_ntheta3,par_ns), INTENT(IN) :: arr
    CHARACTER*(*), INTENT(IN) :: string

    INTEGER :: i, j, k, l

    DO i=left, right
      DO j=1, par_nzeta
        DO k=1, par_ntheta3
          WRITE(fileno+rank,50) string, i, j, k, arr(j,k,i)
          CALL FLUSH(fileno+rank)
        END DO
      END DO
    END DO
    WRITE(fileno+rank,*) 
    CALL FLUSH(fileno+rank)
    IF(ifstop) STOP 'STOPPED PARALLEL CODE'
 50 FORMAT(A,3I6,1P,E24.14)

  END SUBROUTINE PrintParallelIJSPArray
  !------------------------------------------------

  !------------------------------------------------
  ! Prints out [nsmin, nsmax] part of a (ns, nzeta, ntheta) array
  !------------------------------------------------
  SUBROUTINE PrintSerialIJSPArray (arr, left, right, fileno, ifstop, string)
    
    INTEGER, INTENT(IN) :: fileno, left, right
    LOGICAL, INTENT(IN) :: ifstop
    REAL(dp), DIMENSION (par_ns,par_nzeta,par_ntheta3), INTENT(IN) :: arr
    CHARACTER*(*), INTENT(IN) :: string

    INTEGER :: i, j, k, l

    DO i=left, right
      DO j=1, par_nzeta
        DO k=1, par_ntheta3
          WRITE(fileno+rank,50) string, i, j, k, arr(i, j, k)
          CALL FLUSH(fileno+rank)
        END DO
      END DO
    END DO
    WRITE(fileno+rank,*) 
    CALL FLUSH(fileno+rank)
    IF(ifstop) STOP 'STOPPED SERIAL CODE'
 50 FORMAT(A,3I6,1P,E24.14)

  END SUBROUTINE PrintSerialIJSPArray
  !------------------------------------------------

  !------------------------------------------------
  ! Prints out [nsmin, nsmax] part of a (ntor,mpol1,ns) array
  !------------------------------------------------
  SUBROUTINE PrintParallelMNSPArray (arr, left, right, fileno, ifstop, string)
    
    INTEGER, INTENT(IN) :: fileno, left, right
    LOGICAL, INTENT(IN) :: ifstop
    REAL(dp), DIMENSION (0:par_ntor,0:par_mpol1,1:par_ns), INTENT(IN) :: arr
    CHARACTER*(*), INTENT(IN) :: string

    INTEGER :: i, j, k, l

    DO i=left, right
      DO j=0, par_mpol1
        DO k=0, par_ntor
          WRITE(fileno+rank,50) string, i, j, k, arr(k,j,i)
          CALL FLUSH(fileno+rank)
        END DO
      END DO
    END DO
    WRITE(fileno+rank,*) 
    CALL FLUSH(fileno+rank)
    IF(ifstop) STOP 'STOPPED PARALLEL CODE'
 50 FORMAT(A,3I6,1P,E20.12)

  END SUBROUTINE PrintParallelMNSPArray
  !------------------------------------------------

  !------------------------------------------------
  ! Prints out [nsmin, nsmax] part of a (ns, ntor, mpol1) array
  !------------------------------------------------
  SUBROUTINE PrintSerialMNSPArray (arr, left, right, fileno, ifstop, string)

    INTEGER, INTENT(IN) :: fileno, left, right
    LOGICAL, INTENT(IN) :: ifstop
    REAL(dp), DIMENSION (1:par_ns,0:par_ntor,0:par_mpol1), INTENT(IN) :: arr
    CHARACTER*(*), INTENT(IN) :: string

    INTEGER :: i, j, k, l

    DO i=left, right
      DO k=0, par_mpol1
        DO j=0, par_ntor
          WRITE(fileno+rank,50) string, i, j, k, arr(i, j, k)
          CALL FLUSH(fileno+rank)
        END DO
      END DO
    END DO
    WRITE(fileno+rank,*) 
    CALL FLUSH(fileno+rank)
    IF(ifstop) STOP 'STOPPED SERIAL CODE'
 50 FORMAT(A,3I6,1P,E20.12)

  END SUBROUTINE PrintSerialMNSPArray
  !------------------------------------------------

  !------------------------------------------------
  SUBROUTINE Gather4XArray(arr)
    !-----------------------------------------------

    !-----------------------------------------------
    REAL(dp), DIMENSION(0:par_ntor,0:par_mpol1,par_ns,3*par_ntmax), INTENT(INOUT) :: arr
    INTEGER :: i, j, k, l, lk
    INTEGER :: blksize, numjs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: counts, disps
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: sendbuf
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: recvbuf
    REAL(dp) :: allgvton, allgvtoff
    !-----------------------------------------------
    IF (nranks.EQ.1 .OR. .NOT.lactive) THEN
      RETURN
    END IF

    CALL second0(allgvton)

    blksize=(par_ntor+1)*(par_mpol1+1)*3*par_ntmax
    numjs=trglob-tlglob+1
    ALLOCATE (sendbuf(0:par_ntor,0:par_mpol1,numjs,1:3*par_ntmax))
    ALLOCATE (recvbuf(par_ns*blksize))
    ALLOCATE(counts(nranks),disps(nranks))

    DO i=1,nranks
      counts(i)=(trglob_arr(i)-tlglob_arr(i)+1)*blksize
    END DO

    disps(1)=0
    DO i=2,nranks
      disps(i)=disps(i-1)+counts(i-1)
    END DO
#if defined(MPI_OPT)
    sendbuf(0:par_ntor,0:par_mpol1,1:numjs,1:3*par_ntmax)=&
      arr(0:par_ntor,0:par_mpol1,tlglob:trglob,1:3*par_ntmax)
    CALL MPI_Allgatherv(sendbuf,numjs*blksize,MPI_REAL8,recvbuf,&
      counts,disps,MPI_REAL8,NS_COMM,MPI_ERR)
    DO i=1, nranks
      numjs=trglob_arr(i)-tlglob_arr(i)+1
      arr(0:par_ntor,0:par_mpol1,tlglob_arr(i):trglob_arr(i),1:3*par_ntmax)&
        =RESHAPE(recvbuf(disps(i)+1:disps(i)+counts(i)),&
        (/par_ntor+1,par_mpol1+1,numjs,3*par_ntmax/))
    END DO
#endif
    DEALLOCATE(sendbuf, recvbuf)
    DEALLOCATE(counts, disps)
    CALL second0(allgvtoff)
    allgather_time = allgather_time + (allgvtoff-allgvton)
  END SUBROUTINE Gather4XArray
  !------------------------------------------------

  !------------------------------------------------
  !
  !------------------------------------------------
  SUBROUTINE GatherReordered4XArray(arr)
    !-----------------------------------------------

    !-----------------------------------------------
    REAL(dp), DIMENSION(0:par_ntor,0:par_mpol1,3*par_ntmax,par_ns), INTENT(INOUT) :: arr
    INTEGER :: i, j, k, l, lk
    INTEGER :: blksize, numjs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: counts, disps
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: sendbuf
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: recvbuf
    REAL(dp) :: allgvton, allgvtoff
    !-----------------------------------------------
    IF (nranks.EQ.1 .OR. .NOT.lactive) THEN
      RETURN
    END IF

    CALL second0(allgvton)

    blksize=(par_ntor+1)*(par_mpol1+1)*3*par_ntmax
    numjs=trglob-tlglob+1
    ALLOCATE (sendbuf(0:par_ntor,0:par_mpol1,1:3*par_ntmax,numjs))
    ALLOCATE (recvbuf(par_ns*blksize))
    ALLOCATE(counts(nranks),disps(nranks))

    DO i=1,nranks
      counts(i)=(trglob_arr(i)-tlglob_arr(i)+1)*blksize
    END DO

    disps(1)=0
    DO i=2,nranks
      disps(i)=disps(i-1)+counts(i-1)
    END DO
#if defined(MPI_OPT)
    sendbuf(0:par_ntor,0:par_mpol1,1:3*par_ntmax,1:numjs)=&
      arr(0:par_ntor,0:par_mpol1,1:3*par_ntmax,tlglob:trglob)
    CALL MPI_Allgatherv(sendbuf,numjs*blksize,MPI_REAL8,recvbuf,&
      counts,disps,MPI_REAL8,NS_COMM,MPI_ERR)
    DO i=1, nranks
      numjs=trglob_arr(i)-tlglob_arr(i)+1
      arr(0:par_ntor,0:par_mpol1,1:3*par_ntmax,tlglob_arr(i):trglob_arr(i))&
        =RESHAPE(recvbuf(disps(i)+1:disps(i)+counts(i)),&
        (/par_ntor+1,par_mpol1+1,3*par_ntmax,numjs/))
    END DO
#endif
    DEALLOCATE(sendbuf, recvbuf)
    DEALLOCATE(counts, disps)
    CALL second0(allgvtoff)
    allgather_time = allgather_time + (allgvtoff-allgvton)
  END SUBROUTINE GatherReordered4XArray
  !------------------------------------------------


  !------------------------------------------------
  SUBROUTINE Gather2XArray(arr)
    
    REAL(dp), DIMENSION(par_nznt,par_ns), INTENT(INOUT) :: arr
    INTEGER :: i, j, k, l, lk
    INTEGER :: par_nsmin, par_nsmax, blksize, numjs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: counts, disps
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: sendbuf
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: recvbuf
    REAL(dp) :: allgvton, allgvtoff
    !-----------------------------------------------
    IF (nranks.EQ.1 .OR. .NOT.lactive) THEN
      RETURN
    END IF

    CALL second0(allgvton)

    blksize=par_nznt
    numjs=trglob-tlglob+1
    ALLOCATE (sendbuf(par_nznt,numjs))
    ALLOCATE (recvbuf(par_ns*blksize))
    ALLOCATE(counts(nranks),disps(nranks))

    DO i=1,nranks
      counts(i)=(trglob_arr(i)-tlglob_arr(i)+1)*blksize
    END DO

    disps(1)=0
    DO i=2,nranks
      disps(i)=disps(i-1)+counts(i-1)
    END DO

#if defined(MPI_OPT)
    sendbuf(1:par_nznt,1:numjs)=arr(1:par_nznt,tlglob:trglob)
    CALL MPI_Allgatherv(sendbuf,numjs*blksize,MPI_REAL8,recvbuf,&
      counts,disps,MPI_REAL8,NS_COMM,MPI_ERR)
    DO i=1, nranks
      numjs=trglob_arr(i)-tlglob_arr(i)+1
      arr(1:par_nznt,tlglob_arr(i):trglob_arr(i))&
        =RESHAPE(recvbuf(disps(i)+1:disps(i)+counts(i)),&
        (/par_nznt,numjs/))
    END DO
#endif
    DEALLOCATE(sendbuf, recvbuf)
    DEALLOCATE(counts, disps)
    CALL second0(allgvtoff)
    allgather_time = allgather_time + (allgvtoff-allgvton)
  END SUBROUTINE Gather2XArray
  !------------------------------------------------

  !------------------------------------------------
  SUBROUTINE Gather1XArray(arr)
    !-----------------------------------------------
    REAL(dp), DIMENSION(par_ns), INTENT(INOUT) :: arr
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: sendbuf, recvbuf

    INTEGER :: i, numjs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: counts, disps
    REAL(dp) :: allgvton, allgvtoff
    !-----------------------------------------------
    IF (nranks.EQ.1 .OR. .NOT.lactive) THEN
      RETURN
    END IF

    CALL second0(allgvton)

    numjs=trglob-tlglob+1
    ALLOCATE (sendbuf(numjs),recvbuf(par_ns))
    ALLOCATE(counts(nranks),disps(nranks))

    DO i=1,nranks
      counts(i)=(trglob_arr(i)-tlglob_arr(i)+1)
    END DO

    disps(1)=0
    DO i=2,nranks
      disps(i)=disps(i-1)+counts(i-1)
    END DO

#if defined(MPI_OPT)
    sendbuf(1:numjs)=arr(tlglob:trglob)
    CALL MPI_Allgatherv(sendbuf,numjs,MPI_REAL8,recvbuf,&
      counts,disps,MPI_REAL8,NS_COMM,MPI_ERR)
#endif
    arr=recvbuf

    DEALLOCATE(sendbuf, recvbuf)
    DEALLOCATE(counts, disps)
    CALL second0(allgvtoff)
    allgather_time = allgather_time + (allgvtoff-allgvton)
  END SUBROUTINE Gather1XArray
  !------------------------------------------------

  !------------------------------------------------
  !------------------------------------------------

  !------------------------------------------------
  ! Setting the output files
  !------------------------------------------------
  SUBROUTINE SetOutputFile (iam, nprocs, prefix)
    
    INTEGER, INTENT(IN) :: iam, nprocs
    CHARACTER (*), INTENT(IN) :: prefix
    INTEGER :: istat
    CHARACTER(100) :: fname, cfname
    CHARACTER(50) :: ciam, cnprocs
    CHARACTER(25) :: cprefix

    WRITE(ciam,*) iam; WRITE(cnprocs,*) nprocs; WRITE(cprefix,*) prefix
    ciam=ADJUSTL(ciam); cnprocs=ADJUSTL(cnprocs); cprefix=ADJUSTL(cprefix)
    TOFU = 4*nprocs+iam+1000

    fname=TRIM(cprefix)//'.txt'
    OPEN(UNIT=TOFU, FILE=fname, STATUS="REPLACE", ACTION="WRITE",&
       FORM="FORMATTED",POSITION="APPEND", IOSTAT=istat)

  END SUBROUTINE SetOutputFile
  !------------------------------------------------

  !------------------------------------------------
      SUBROUTINE tolastns (inarr,outarr)
      
      REAL(dp), INTENT(IN), & 
        DIMENSION(0:par_ntor,0:par_mpol1,par_ns,3*par_ntmax) :: inarr

#if defined(SKS) 
      INTEGER :: i, j, k, l, cnt
      REAL(dp), INTENT(INOUT), & 
     &   DIMENSION(ntmaxblocksize,par_ns) :: outarr

      DO i=t1lglob, t1rglob
        cnt=0
        DO l=1, 3*par_ntmax
          DO k=0, par_mpol1
            DO j=0, par_ntor
              cnt=cnt+1
              outarr(cnt,i)=inarr(j,k,i,l)
            END DO
          END DO
        END DO
      END DO
#else
      INTEGER :: js, t, m, n, lk
      REAL(dp), INTENT(INOUT), & 
        DIMENSION(ntmaxblocksize*par_ns) :: outarr

      DO js=t1lglob, t1rglob
        DO t=1, 3*par_ntmax
          DO m=0, par_mpol1
            DO n=0, par_ntor
              lk = n+nmax*(m+mmax*((t-1)+tmax*(js-1)))+1
              outarr(lk)=inarr(n,m,js,t)
            END DO
          END DO
        END DO
      END DO
#endif
      END SUBROUTINE tolastns
  !------------------------------------------------

  !------------------------------------------------
      SUBROUTINE tolastntype (inarr,outarr)
      REAL(dp), INTENT(INOUT), &
        DIMENSION(0:par_ntor,0:par_mpol1,par_ns,3*par_ntmax) :: outarr

#if defined(SKS) 
      REAL(dp), INTENT(IN), &
     &   DIMENSION(ntmaxblocksize,par_ns) :: inarr
      INTEGER :: i, j, k, l, cnt

      DO i=t1lglob, t1rglob
        cnt=0
        DO l=1, 3*par_ntmax
          DO k=0, par_mpol1
            DO j=0, par_ntor
              cnt=cnt+1
              outarr(j,k,i,l)=inarr(cnt,i)
            END DO
          END DO
        END DO
      END DO
#else
      INTEGER :: js, t, m, n, lk
      REAL(dp), INTENT(IN), & 
     &   DIMENSION(ntmaxblocksize*par_ns) :: inarr
      
      DO js=t1lglob, t1rglob
        DO t=1, 3*par_ntmax
          DO m=0, par_mpol1
            DO n=0, par_ntor
              lk = n+nmax*(m+mmax*((t-1)+tmax*(js-1)))+1
              outarr(n,m,js,t)=inarr(lk)
            END DO
          END DO
        END DO
      END DO
#endif
      END SUBROUTINE tolastntype
  !------------------------------------------------

  !------------------------------------------------
  ! 
  !------------------------------------------------
      SUBROUTINE copylastns(a1, a2)
        REAL(dp), INTENT(IN),  DIMENSION(ntmaxblocksize,par_ns)    :: a1
        REAL(dp), INTENT(INOUT),  DIMENSION(ntmaxblocksize,par_ns) :: a2
        a2(:,t1lglob:t1rglob) = a1(:,t1lglob:t1rglob)
      END SUBROUTINE copylastns
      !------------------------------------------------

      !------------------------------------------------
      SUBROUTINE copy1lastntype(a1, a2)
        REAL(dp), INTENT(IN),    DIMENSION(ntmaxblocksize*par_ns) :: a1
        REAL(dp), INTENT(INOUT), DIMENSION(ntmaxblocksize*par_ns) :: a2
        CALL copy4lastntype(a1, a2)
      END SUBROUTINE copy1lastntype
      !------------------------------------------------

      !------------------------------------------------
      SUBROUTINE copy4lastntype(a1, a2)
        REAL(dp), INTENT(IN),  DIMENSION(0:par_ntor,0:par_mpol1,par_ns,3*par_ntmax) :: a1
        REAL(dp), INTENT(INOUT),  DIMENSION(0:par_ntor,0:par_mpol1,par_ns,3*par_ntmax) :: a2
        a2(:,:,t1lglob:t1rglob,:) = a1(:,:,t1lglob:t1rglob,:)
      END SUBROUTINE copy4lastntype
      !------------------------------------------------

      !------------------------------------------------
      SUBROUTINE copym1lastntype(a1, a2, d1)
        REAL(dp), INTENT(IN),  DIMENSION(ntmaxblocksize*par_ns) :: a1
        REAL(dp), INTENT(INOUT),  DIMENSION(ntmaxblocksize*par_ns) :: a2
        REAL(dp), INTENT(IN) :: d1
        CALL copym4lastntype(a1, a2, d1)
      END SUBROUTINE copym1lastntype
      !------------------------------------------------

      !------------------------------------------------
      SUBROUTINE copym4lastntype(a1, a2, d1)
        REAL(dp), INTENT(IN),  DIMENSION(0:par_ntor,0:par_mpol1,par_ns,3*par_ntmax) :: a1
        REAL(dp), INTENT(INOUT),  DIMENSION(0:par_ntor,0:par_mpol1,par_ns,3*par_ntmax) :: a2
        REAL(dp), INTENT(IN) :: d1
        a2(:,:,t1lglob:t1rglob,:) = d1*a1(:,:,t1lglob:t1rglob,:)
      END SUBROUTINE copym4lastntype
      !------------------------------------------------

      !------------------------------------------------
      SUBROUTINE saxlastns(a, x, vec)
        REAL(dp), INTENT(IN),  DIMENSION(ntmaxblocksize,par_ns) :: a, x
        REAL(dp), INTENT(INOUT),  DIMENSION(ntmaxblocksize,par_ns) :: vec
        vec(:,tlglob:trglob) = a(:,tlglob:trglob)*x(:,tlglob:trglob)
      END SUBROUTINE saxlastns
      !------------------------------------------------

      !------------------------------------------------
      SUBROUTINE saxlastns1(a, x, vec)
        REAL(dp), INTENT(IN),  DIMENSION(ntmaxblocksize,tlglob:trglob) :: a
        REAL(dp), INTENT(IN),  DIMENSION(ntmaxblocksize,par_ns) :: x
        REAL(dp), INTENT(INOUT),  DIMENSION(ntmaxblocksize,tlglob:trglob) :: vec
        vec(:,tlglob:trglob) = a(:,tlglob:trglob)*x(:,tlglob:trglob)
      END SUBROUTINE saxlastns1

      !------------------------------------------------
      SUBROUTINE saxlastntype(a, x, vec)
        REAL(dp), INTENT(IN),  DIMENSION(blocksize,par_ns,3*par_ntmax) :: a, x
        REAL(dp), INTENT(INOUT),  DIMENSION(blocksize,par_ns,3*par_ntmax) :: vec
        vec(:,t1lglob:t1rglob,:) = a(:,t1lglob:t1rglob,:)*x(:,t1lglob:t1rglob,:)
      END SUBROUTINE saxlastntype
      !------------------------------------------------

      !------------------------------------------------
      SUBROUTINE saxpbylastntype(a, x, b, y, vec)
        REAL(dp), INTENT(IN),  DIMENSION(blocksize,par_ns,3*par_ntmax) :: x, y
	REAL(dp), INTENT(INOUT),  DIMENSION(blocksize,par_ns,3*par_ntmax) :: vec
	REAL(dp), INTENT(IN)  :: a, b
        vec(:,t1lglob:t1rglob,:) = a*x(:,t1lglob:t1rglob,:) + b*y(:,t1lglob:t1rglob,:)
      END SUBROUTINE saxpbylastntype
      !------------------------------------------------

      !------------------------------------------------
      SUBROUTINE saxpbylastns(a, x, b, y, vec)
        REAL(dp), INTENT(IN),  DIMENSION(ntmaxblocksize,tlglob:trglob) :: x
        REAL(dp), INTENT(IN),  DIMENSION(ntmaxblocksize,par_ns) :: y
        REAL(dp), INTENT(INOUT),  DIMENSION(ntmaxblocksize,par_ns) :: vec
        REAL(dp), INTENT(IN)  :: a, b
        vec(:,tlglob:trglob) = a*x(:,tlglob:trglob) + b*y(:,tlglob:trglob)
      END SUBROUTINE saxpbylastns
      !------------------------------------------------

      !------------------------------------------------
      SUBROUTINE saxpby1lastns(a, x, b, y, vec)
        REAL(dp), INTENT(IN),  DIMENSION(ntmaxblocksize,par_ns) :: x, y
        REAL(dp), INTENT(INOUT),  DIMENSION(ntmaxblocksize,par_ns) :: vec
        REAL(dp), INTENT(IN)  :: a, b
        vec(:,tlglob:trglob) = a*x(:,tlglob:trglob) + b*y(:,tlglob:trglob)
      END SUBROUTINE saxpby1lastns
      !------------------------------------------------

      !------------------------------------------------
      SUBROUTINE getderivlastns(x1, x2, d1, vec)
        REAL(dp), INTENT(IN),  DIMENSION(ntmaxblocksize,par_ns) :: x1, x2
        REAL(dp), INTENT(INOUT),  DIMENSION(ntmaxblocksize,tlglob:trglob) :: vec
        REAL(dp), INTENT(IN)  :: d1
        IF (d1 .eq. 0) RETURN
        vec = (x1(:,tlglob:trglob) - x2(:,tlglob:trglob))/d1
      END SUBROUTINE getderivlastns
      !------------------------------------------------

      !------------------------------------------------
      SUBROUTINE zerolastntype(a1)
        REAL(dp), INTENT(INOUT),  DIMENSION(blocksize,par_ns,3*par_ntmax) :: a1
        a1(:,t1lglob:t1rglob,:) = 0
      END SUBROUTINE zerolastntype
      !------------------------------------------------

      !------------------------------------------------
      SUBROUTINE CopyParallelLinearSubarray(fromarr,toarr, first, last)
        INTEGER, INTENT(IN) :: first, last
        REAL(dp), INTENT(IN), &
             DIMENSION(blocksize,par_ns,3*par_ntmax) :: fromarr
        REAL(dp), INTENT(INOUT), &
             DIMENSION(blocksize,par_ns,3*par_ntmax) :: toarr
        toarr(:,first:last,:) = fromarr(:,first:last,:)
      END SUBROUTINE CopyParallelLinearSubarray
      !------------------------------------------------

      !------------------------------------------------
      SUBROUTINE PadSides(arr)
        REAL(dp), INTENT(INOUT), &
             DIMENSION(blocksize,par_ns,3*par_ntmax) :: arr
        INTEGER :: left, right, tag1=1
        REAL(dp) :: ton, toff
#if defined(SKS)
        CALL second0(ton)

        left=rank-1;  IF(rank.EQ.0) left=MPI_PROC_NULL
        right=rank+1; IF(rank.EQ.nranks-1) right=MPI_PROC_NULL
        CALL MPI_Sendrecv(arr(:,tlglob,:),ntmaxblocksize,MPI_REAL8,    &
                 left, tag1,arr(:,t1rglob,:),ntmaxblocksize,MPI_REAL8, &
                 right,tag1,NS_COMM, MPI_STAT,MPI_ERR)
        CALL MPI_Sendrecv(arr(:,trglob,:),ntmaxblocksize,MPI_REAL8,    &
                 right,tag1,arr(:,t1lglob,:),ntmaxblocksize,MPI_REAL8, &
                 left, tag1,NS_COMM,MPI_STAT,MPI_ERR)

        CALL second0(toff)
        sendrecv_time = sendrecv_time + (toff - ton)
#endif
      END SUBROUTINE PadSides
      !------------------------------------------------

      !------------------------------------------------
      SUBROUTINE PadSides1X(arr)
        REAL(dp), INTENT(INOUT), DIMENSION(par_ns) :: arr
        INTEGER :: left, right, tag1=1
        REAL(dp) :: ton, toff
#if defined(SKS)
        CALL second0(ton)

        left=rank-1;  IF(rank.EQ.0) left=MPI_PROC_NULL
        right=rank+1; IF(rank.EQ.nranks-1) right=MPI_PROC_NULL

        IF (grank .LT. nranks) THEN
          CALL MPI_Sendrecv(arr(tlglob),1,MPI_REAL8,left,tag1,  &
               arr(t1rglob),1,MPI_REAL8,right,tag1,NS_COMM,     &
               MPI_STAT, MPI_ERR)
        END IF

        CALL second0(toff)
        sendrecv_time = sendrecv_time + (toff - ton)
#endif
      END SUBROUTINE PadSides1X
      !------------------------------------------------

      !------------------------------------------------
      SUBROUTINE CompareEdgeValues(pxc, pxsave)
        REAL(dp), INTENT(IN), DIMENSION(blocksize,par_ns,3*par_ntmax)   &
                       :: pxc, pxsave
    
       IF (rank.EQ.nranks-1 .AND. ANY(pxc(:,par_ns,:) .NE. pxsave(:,par_ns,:)))     &
         PRINT *,' xsave != xc at edge returning from GMRES'

      END SUBROUTINE CompareEdgeValues

      !------------------------------------------------


END MODULE parallel_vmec_module
!------------------------------------------------
