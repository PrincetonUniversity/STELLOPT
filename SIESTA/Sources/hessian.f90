      MODULE hessian
      USE island_params, mpol=>mpol_i, ntor=>ntor_i, ns=>ns_i,           &
          nfp=>nfp_i
      USE timer_mod
      USE DirectAccess
      USE stel_constants, ONLY: pi
      USE descriptor_mod
#if defined(MPI_OPT)
      USE prof_mod
      USE ptrd_mod
#endif
#if defined(SKS)
      USE nscalingtools, ONLY: PARSOLVER, PARFUNCTISL,                   &
        rcounts, disp, startglobrow, endglobrow,                         &
        SKSDBG, TOFU, UPPER, LOWER, DIAG, SAVEDIAG, nrecd,               &
        SYMMETRYCHECK, totRecvs, PACKSIZE, WriteTime, send, receive,     &
        SetBlockTriDataStruct, GetFullSolution
      USE blocktridiagonalsolver, ONLY: Initialize, ForwardSolve,        &
        SetMatrixRHS, BackwardSolve, GetSolutionVector, CheckSymmetry
#else
      USE nscalingtools, ONLY: PARSOLVER, PARFUNCTISL, SKSDBG, TOFU
#endif
      IMPLICIT NONE
#if defined(MPI_OPT)
      INCLUDE 'mpif.h'
#endif
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: jstart(3) = (/1,2,3/)
      INTEGER   :: mblk_size  
      INTEGER   :: iprec_flag = 0
      INTEGER, ALLOCATABLE :: ipiv_blk(:,:)

      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: ublk, dblk, lblk
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: ublk_s, dblk_s, lblk_s
      REAL(rprec), ALLOCATABLE :: mupar_norm(:)
      REAL(rprec)  :: levmarq_param, levmarq_param0, asym_index, mupar, mupar0
      REAL(rprec), PARAMETER :: eps_factor = 1._dp 
      CHARACTER*10 :: string

      LOGICAL :: l_Compute_Hessian
      LOGICAL :: l_Diagonal_Only = .FALSE.
      LOGICAL :: l_backslv = .FALSE. 

      CHARACTER(LEN=3)   :: FlashDrive =""
      CHARACTER(LEN=128) :: ScratchFile
      INTEGER, PARAMETER :: blmin=3, bldia=2, blpls=1
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: DataItem
      INTEGER :: iunit_dacess=10
      LOGICAL :: lswap2disk=.FALSE.
#if defined(SKS)
      INTEGER :: HESSPASS=0
      INTEGER :: myStartBlock, myEndBlock, myNumBlocks 
      INTEGER :: mystart(3), myend(3)
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: SavedDiag
#endif
      PRIVATE :: block_precond             !diag_precond 
      PUBLIC  :: apply_precond
!-----------------------------------------------

      CONTAINS

      SUBROUTINE Compute_Hessian_Blocks (xc, gc, func,                  &
                 l_Hess_sym, ldiagonal)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ns*(1+mpol)*(2*ntor+1)*3),                 &
               INTENT(OUT)    :: xc, gc
      LOGICAL, INTENT(IN)     :: l_Hess_sym
      LOGICAL, INTENT(IN)     :: ldiagonal
      EXTERNAL                :: func
!-----------------------------------------------
      l_Compute_Hessian = .TRUE.
      l_Diagonal_Only   = ldiagonal

      IF (l_Diagonal_Only .AND. iam.EQ.0) THEN
         WRITE (6, '(/,a)') ' Computing diagonal preconditioner!'
         WRITE (33,'(/,a)') ' Computing diagonal preconditioner!'
      ENDIF

      IF(PARSOLVER) THEN
        IF(PARFUNCTISL) THEN
          IF(SKSDBG) WRITE(TOFU,*) 'Executing Compute_Hessian_Blocks_With_No_Col_Redist'; IF(SKSDBG) CALL FLUSH(TOFU)
          CALL Compute_Hessian_Blocks_With_No_Col_Redist (xc, gc, func, l_Hess_sym)
        ELSE
          IF(SKSDBG) WRITE(TOFU,*) 'Executing Compute_Hessian_Blocks_With_Col_Redist'; IF(SKSDBG) CALL FLUSH(TOFU)
          CALL Compute_Hessian_Blocks_With_Col_Redist (xc, gc, func, l_Hess_sym)
        END IF 
      ELSE 
        IF(SKSDBG) WRITE(TOFU,*) 'Executing Compute_Hessian_Blocks_Thomas'; CALL FLUSH(TOFU)
        CALL Compute_Hessian_Blocks_Thomas (xc, gc, func, l_Hess_sym)
      END IF
      l_Diagonal_Only   = .FALSE.
      l_Compute_Hessian = .FALSE.

      END SUBROUTINE Compute_Hessian_Blocks

      SUBROUTINE Compute_Hessian_Blocks_With_No_Col_Redist (xc, gc,     &
                 func, l_Hess_sym)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION((1+mpol)*(2*ntor+1)*3,ns) :: xc, gc
      LOGICAL, INTENT(IN) :: l_Hess_sym
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!----------------------------------------------
      REAL(rprec), PARAMETER    :: zero=0, one=1, p5=0.5_dp
      INTEGER  :: n, m, js, mesh, ntype, istat, iunit, icol
      REAL(rprec) :: eps, ton, toff, bsize, lmdamp
      CHARACTER(LEN=3) :: label
      EXTERNAL    :: func
      INTEGER  :: js1
#if defined(SKS)
      REAL(rprec) :: starttime, endtime, usedtime
      REAL(rprec) :: colstarttime, colendtime, colusedtime
      LOGICAL :: lfirst
      INTEGER :: nsmin, nsmax
      REAL(dp) :: skston, skstoff, temp
!----------------------------------------------

      starttime=0; endtime=0; usedtime=0; 
      colstarttime=0; colendtime=0; colusedtime=0

      !------------------------------------------------------------
      !     COMPUTE (U)pper, (D)iagonal, and (L)ower Block matrices
      !------------------------------------------------------------

      IF (iam .eq. 0 .AND. (.not. l_diagonal_only)) WRITE (6, 50)
      50   FORMAT (/,' Computing block preconditioner')

      CALL second0(ton)

      mblk_size = SIZE(gc,1)
      eps = SQRT(EPSILON(eps))*ABS(eps_factor)

      ALLOCATE (DataItem(mblk_size), stat=istat)
      IF (istat .ne. 0) STOP 'DataItem allocation failed!'
      DataItem = 0
      ALLOCATE (SavedDiag(mblk_size), stat=istat)
      IF (istat .ne. 0) STOP 'SavedDiag allocation failed!'
      SavedDiag = 0

      IF (iprec_flag == 0) THEN
        bsize = 3*ns*KIND(dblk)
        bsize = bsize*REAL(mblk_size*mblk_size,dp) 
        IF (bsize .lt. 1.E6_dp) THEN
          bsize = bsize/1.E3
          label = " Kb"
        ELSE IF (bsize .lt. 1.E9_dp) THEN
          bsize = bsize/1.E6_dp
          label = " Mb"
        ELSE
          bsize = bsize/1.E9_dp
          label = " Gb"
        END IF

        IF (iam .eq. 0) THEN
          DO iunit = 6,33,27
            WRITE (iunit, '(1x,a,i4,a,f6.2,a)') 'Block dim: ', mblk_size, &
            &         '^2  Preconditioner size: ', bsize, TRIM(label)
          END DO
        ENDIF 

        iprec_flag = 1

      END IF

      CALL second0(skston)
      xc = 0
      CALL func
      CALL second0(skstoff)
      hessian_funct_island_time=hessian_funct_island_time+(skstoff-skston)
      nsmin = MAX(1, startglobrow-1);  nsmax = MIN(ns, endglobrow+1)

      IF (HESSPASS .EQ. 0) THEN
        myStartBlock=startglobrow
        myEndBlock=endglobrow
        myNumBlocks=myEndBlock-myStartBlock+1
        CALL Initialize(.FALSE.,ns,mblk_size) 
        DO mesh = 1, 3
           lfirst = .TRUE.
           mystart(mesh) = nsmin
           myend(mesh) = nsmax
           DO js = jstart(mesh), ns, 3
              IF (js .LT. nsmin) CYCLE
              IF (lfirst) mystart(mesh)=js
              lfirst=.FALSE.
              myend(mesh)=MIN(js, nsmax)
              IF (js .GE. nsmax) EXIT
           END DO
        END DO
      ENDIF

      !----------------------------
      ! Column generation begins
      ! REVERSE js ORDER IN XC, GC AND INTERNALLY IN ALL _PAR ROUTINES
      !----------------------------
      CALL second0(colstarttime)

      icol=0
      VAR_TYPE_LG: DO ntype = 1, 3
        VAR_N_LG: DO n = -ntor, ntor
          VAR_M_LG: DO m = 0, mpol

            icol=icol+1

            MESH_3PT_LG: DO mesh = 1, 3

              xc = 0

              DO js = mystart(mesh), myend(mesh), 3
                 xc(icol,js) = eps
              END DO

              INHESSIAN=.TRUE.
              CALL second0(skston)
              CALL func
              CALL second0(skstoff)
              hessian_funct_island_time=hessian_funct_island_time+(skstoff-skston)
              INHESSIAN=.FALSE.

              SKIP3_MESH_LG: DO js = mystart(mesh), myend(mesh), 3

                lmdamp = levmarq_param*gc(icol,js)/eps
                lmdamp = ABS(lmdamp)

                IF (levmarq_param .NE. zero) THEN
                  IF (lmdamp .EQ. zero) THEN
                    lmdamp = levmarq_param
                  ELSE IF (ntype .NE. 3) THEN
                    lmdamp = lmdamp/10
                  END IF
                END IF

                !ublk(js-1)
                js1=js-1
                IF (startglobrow.LE.js1 .AND. js1.LE.endglobrow) THEN
                   DataItem = gc(:,js1)/eps
                   IF (m.EQ.0 .AND. n.LT.0) DataItem = 0
                   IF (l_diagonal_only) THEN
                      temp=DataItem(icol)
                      DataItem=0
                      DataItem(icol)=temp
                   END IF
                   CALL SetBlockTriDataStruct(UPPER,js1,icol,DataItem) 
                END IF

                !dblk(js)

                IF (startglobrow.LE.js .AND. js.LE.endglobrow) THEN
                   DataItem = gc(:,js)/eps
                   IF (l_diagonal_only) THEN
                      temp=DataItem(icol)
                      DataItem=0
                      DataItem(icol)=temp
                   END IF
                   CALL SetBlockTriDataStruct(SAVEDIAG,js,icol,DataItem) 
                   IF (ALL(DataItem .EQ. zero)) THEN
                      DataItem(icol) = -one
                   END IF

                   DataItem(icol) = DataItem(icol) + SIGN(lmdamp,DataItem(icol))
                   CALL SetBlockTriDataStruct(DIAG,js,icol,DataItem) 
                END IF

                !lblk(js+1) 
                js1 = js+1
                IF (startglobrow.LE.js1 .AND. js1.LE.endglobrow) THEN
                   DataItem = gc(:,js1)/eps
                   IF (m.EQ.0 .AND. n.LT.0) DataItem = 0
                   IF (l_diagonal_only) THEN
                      temp=DataItem(icol)
                      DataItem=0
                      DataItem(icol)=temp
                   END IF
                   CALL SetBlockTriDataStruct(LOWER,js1,icol,DataItem) 
                END IF

              END DO SKIP3_MESH_LG
            END DO MESH_3PT_LG
          END DO VAR_M_LG
        END DO VAR_N_LG
      END DO VAR_TYPE_LG

      CALL second0(colendtime)
      construct_hessian_time=construct_hessian_time+(colendtime-colstarttime)

!      STOP 'END HESSIAN LOOP TEST WRITE FORT.3000'
      
      DEALLOCATE (DataItem, stat=istat)
      IF (istat .NE. 0) STOP 'DataItem deallocation error!'
      DEALLOCATE (SavedDiag, stat=istat)
      IF (istat .NE. 0) STOP 'SavedDiag deallocation error!'
      !----------------------------
      ! Column generation ends
      !----------------------------

      !IF(DIAGONALDONE) CALL WriteBlocks (HESSPASS)

      !
      !     CHECK SYMMETRY OF BLOCKS 
      !

      IF (SYMMETRYCHECK) THEN
        CALL second0(skston)
        CALL CheckSymmetry(asym_index)
        CALL second0(skstoff)
        asymmetry_check_time=asymmetry_check_time+(skstoff-skston)
      ENDIF      

      IF (iam .EQ. 0) THEN
        WRITE (6, 100) levmarq_param, muPar, asym_index
        WRITE (33,110) levmarq_param, muPar, asym_index
      ENDIF

 100  FORMAT(' LM parameter: ',1pe9.2,                                  & 
             ' mu||: ',1pe9.2, ' Asym index: ',1pe9.2)
 110  FORMAT (/,' Computing block preconditioner: ',                    &
              ' LM parameter: ',1pe9.2,' mu||: ',1pe9.2,                &
              ' Asym index: ',1pe9.2)

      !
      !     SYMMETRIZE BLOCKS
      !
!      IF (.not. l_Hess_sym .or. lswap2disk) GOTO 1000
!      1000 CONTINUE

      CALL second0(toff)
      IF (l_Diagonal_Only) THEN
         time_diag_prec = time_diag_prec+(toff-ton)
      ELSE
         time_block_prec = time_block_prec+(toff-ton)
      END IF

!      IF (l_Diagonal_Only) RETURN         !ForwardSolve will be called first time but safer this way!

      ton = toff
      skston = ton
      !
      !     FACTORIZE (Invert) HESSIAN
      !
      IF (ALLOCATED(ipiv_blk)) DEALLOCATE(ipiv_blk,stat=istat)
      CALL ForwardSolve
      CALL second0(toff)
      skstoff=toff
      block_factorization_time=block_factorization_time+(skstoff-skston)

      HESSPASS = HESSPASS+1
      time_factor_blocks = time_factor_blocks + (toff-ton)
#endif
      END SUBROUTINE Compute_Hessian_Blocks_With_No_Col_Redist

      SUBROUTINE Compute_Hessian_Blocks_With_Col_Redist (xc, gc, func,   &
                 l_Hess_sym)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION((1+mpol)*(2*ntor+1)*3,ns) :: xc, gc
      LOGICAL, INTENT(IN) :: l_Hess_sym
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!----------------------------------------------
      REAL(rprec), PARAMETER    :: zero=0, one=1, p5=0.5_dp
      INTEGER  :: n, m, js, js1, mesh, ntype, istat, iunit, icol, icolmpi
      REAL(rprec) :: eps, ton, toff, bsize, lmdamp
      CHARACTER(LEN=3) :: label
      EXTERNAL    :: func
#if defined(SKS)
      INTEGER, ALLOCATABLE, DIMENSION (:) :: sendCount, recvCount
      REAL(rprec), ALLOCATABLE, DIMENSION (:) :: recvBuf
      INTEGER :: procID, MPI_ERR
      REAL(rprec) :: starttime, endtime, usedtime
      REAL(rprec) :: colstarttime, colendtime, colusedtime

      LOGICAL :: PROBEFLAG
      REAL(dp) :: skston, skstoff, temp

      INTEGER :: it 
!----------------------------------------------
      starttime=zero
      endtime=zero
      usedtime=zero
!------------------------------------------------------------
!     COMPUTE (U)pper, (D)iagonal, and (L)ower Block matrices
!------------------------------------------------------------
      IF (iam.eq.0 .AND. .NOT.l_Diagonal_Only) WRITE (6, 50)
 50   FORMAT (/,' Computing block preconditioner')

      CALL second0(ton)

      mblk_size = SIZE(gc,1)

      eps = SQRT(EPSILON(eps))*ABS(eps_factor)

      ALLOCATE (DataItem(mblk_size), stat=istat)
      IF (istat .ne. 0) STOP 'DataItem allocation failed!'
      DataItem = 0
      ALLOCATE (SavedDiag(mblk_size), stat=istat)
      IF (istat .ne. 0) STOP 'SavedDiag allocation failed!'
      SavedDiag = zero

      IF (iprec_flag == 0) THEN
         bsize = 3*ns*KIND(dblk)
         bsize = bsize*REAL(mblk_size*mblk_size,dp) 
         IF (bsize .lt. 1.E6_dp) THEN
            bsize = bsize/1.E3
            label = " Kb"
         ELSE IF (bsize .lt. 1.E9_dp) THEN
            bsize = bsize/1.E6_dp
            label = " Mb"
         ELSE
            bsize = bsize/1.E9_dp
            label = " Gb"
         END IF

         IF (iam .eq. 0) THEN
            DO iunit = 6,33,27
            WRITE (iunit, '(1x,a,i4,a,f6.2,a)') 'Block dim: ', mblk_size, &
               '^2  Preconditioner size: ', bsize, TRIM(label)
            END DO
         ENDIF 
         iprec_flag = 1
      END IF

      xc = 0

      icolmpi = 0
      icol=0

      CALL second0(skston)
      CALL func
      CALL second0(skstoff)
      hessian_funct_island_time=hessian_funct_island_time+(skstoff-skston)

      IF (HESSPASS.EQ.0) THEN
        myStartBlock=startglobrow
        myEndBlock=endglobrow
        myNumBlocks=myEndBlock-myStartBlock+1
        CALL Initialize(.FALSE.,ns,mblk_size) 
      ENDIF

      IF(.NOT.ALLOCATED(sendCount)) ALLOCATE(sendCount(nprocs))
      IF(.NOT.ALLOCATED(recvCount)) ALLOCATE(recvCount(nprocs))
      IF(.NOT.ALLOCATED(recvBuf)) ALLOCATE(recvBuf(PACKSIZE))
      sendCount=0
      recvCount=0
      nrecd=0
      totRecvs=0
      PROBEFLAG=.TRUE.

      CALL second0(colstarttime)
      VAR_TYPE: DO ntype = 1, 3
         VAR_N: DO n = -ntor, ntor
            VAR_M: DO m = 0, mpol
                  
               icol=icol+1
               IF(MOD(icol-1,nprocs)==iam) THEN
                  icolmpi = icolmpi+1

                  MESH_3PT: DO mesh = 1, 3

                  xc = 0

                  DO js = jstart(mesh), ns, 3
                     xc(icol,js) = eps
                  END DO

                  INHESSIAN=.TRUE.
                  CALL second0(skston)
                  CALL func
                  CALL second0(skstoff)
                  hessian_funct_island_time=hessian_funct_island_time+(skstoff-skston)
                  INHESSIAN=.FALSE.

                  SKIP3_MESH: DO js = jstart(mesh), ns, 3

                    lmdamp = levmarq_param*gc(icol,js)/eps
                    lmdamp = ABS(lmdamp)

                    IF (levmarq_param .ne. zero) THEN
                      IF (lmdamp .eq. zero) THEN
                        lmdamp = levmarq_param
                      ELSE IF (ntype .ne. 3) THEN
                        lmdamp = lmdamp/10                          !Emulate v|| damping only
                      END IF
                    END IF
 
                     !ublk(js-1)
                     js1 = js-1
                     IF (js1 .gt. 0) THEN
                        DataItem = gc(:,js1)/eps
                        !m=0 constraint for n<0
                        IF (m.eq.0 .and. n.lt.0) DataItem = 0
                        IF (l_diagonal_only) THEN
                           temp=DataItem(icol)
                           DataItem=0
                           DataItem(icol)=temp
                        END IF
                        istat = 1
                        CALL receive(PROBEFLAG)
                        CALL send(DataItem,js1,istat,icol,procID)
                        IF(procID-1.NE.iam) sendCount(procID) = sendCount(procID)+1
                        CALL receive(PROBEFLAG)
                     END IF !js>1

                     !dblk(js)
                     DataItem = gc(:,js)/eps
                     IF (m.eq.0 .and. n.lt.0) DataItem = 0
                     IF (l_diagonal_only) THEN
                        temp=DataItem(icol)
                        DataItem=0
                        DataItem(icol)=temp
                     END IF
                     SavedDiag=DataItem; istat=4
                     CALL receive(PROBEFLAG)
                     CALL send(SavedDiag,js,istat,icol,procID)
                     IF(procID-1.NE.iam) sendCount(procID) = sendCount(procID)+1
                     CALL receive(PROBEFLAG)
!                    Boundary condition at js=1 and ns (CATCH ALL OF THEM HERE)
!                    and m=0,n<0: NEED THIS to avoid (near) singular Hessian
!                    ASSUMES DIAGONALS ARE ALL NEGATIVE
                     IF (ALL(DataItem .eq. zero)) DataItem(icol) = -one
                     DataItem(icol) = DataItem(icol) + SIGN(lmdamp,DataItem(icol))
                     istat = 2; js1 = js
                     CALL receive(PROBEFLAG)
                     CALL send(DataItem,js1,istat,icol,procID)
                     IF(procID-1.NE.iam) sendCount(procID) = sendCount(procID)+1
                     CALL receive(PROBEFLAG)
                   
                     !lblk(js+1)
                     js1 =js+1 
                     IF (js .lt. ns) THEN
                       DataItem = gc(:,js1)/eps
                       !m=0 constraint for n<0
                       IF (m.eq.0 .and. n.lt.0) DataItem=0
                       IF (l_diagonal_only) THEN
                          temp=DataItem(icol)
                          DataItem=0
                          DataItem(icol)=temp
                       END IF
                       istat = 3
                       CALL receive(PROBEFLAG)
                       CALL send(DataItem,js1,istat,icol,procID)
                       IF(procID-1.NE.iam) sendCount(procID) = sendCount(procID)+1
                       CALL receive(PROBEFLAG)
                     END IF !JS < NS

                  END DO SKIP3_MESH
               END DO MESH_3PT
               ENDIF
            END DO VAR_M
         END DO VAR_N
      END DO VAR_TYPE
      CALL second0(colendtime)
      colusedtime=colendtime-colstarttime

      construct_hessian_time=construct_hessian_time+(colendtime-colstarttime)

      CALL second0(skston)
      IF (PARSOLVER) THEN

        PROBEFLAG=.NOT.PROBEFLAG

        CALL MPI_Alltoall(sendCount,1,MPI_INTEGER,recvCount,1,MPI_INTEGER,MPI_COMM_WORLD,MPI_ERR)

        totRecvs=0
        DO js=1,nprocs,1
          totRecvs=totRecvs+recvCount(js)
        END DO

        DO WHILE (nrecd.LT.totRecvs)
          CALL receive(PROBEFLAG)
        END DO
      END IF

      DEALLOCATE (DataItem, stat=istat)
      DEALLOCATE (SavedDiag, stat=istat)
      IF (istat .ne. 0) STOP 'DataItem deallocation error!'

      CALL second0(skstoff)
      construct_hessian_time=construct_hessian_time+(skstoff-skston)
      CALL MPI_Barrier(MPI_COMM_WORLD, MPI_ERR)

      IF (lswap2disk) GO TO 900
!
!     CHECK SYMMETRY OF BLOCKS 
!
      IF (SYMMETRYCHECK) THEN
        IF (PARSOLVER) THEN 
          CALL second0(skston)
          CALL CheckSymmetry(asym_index)
          CALL second0(skstoff)
          asymmetry_check_time=asymmetry_check_time+(skstoff-skston)
        ENDIF
      ENDIF      

 900  CONTINUE

      IF (iam .eq. 0) THEN
         WRITE (6, 100) levmarq_param, muPar, asym_index
         WRITE (33,110) levmarq_param, muPar, asym_index
      ENDIF

 100  FORMAT(' LM parameter: ',1pe9.2,                                  & 
             ' mu||: ',1pe9.2, ' Asym index: ',1pe9.2)
 110  FORMAT (/,' Computing block preconditioner: ',                    &
              ' LM parameter: ',1pe9.2,' mu||: ',1pe9.2,                &
              ' Asym index: ',1pe9.2)
 
!
!     SYMMETRIZE BLOCKS
!
!      IF (.not. l_Hess_sym .or. lswap2disk) GOTO 1000
! 1000 CONTINUE

!
!     FACTORIZE (Invert) HESSIAN
!
      CALL second0(toff)
      IF (l_Diagonal_Only) THEN
         time_diag_prec = time_diag_prec+(toff-ton)
      ELSE
         time_block_prec = time_block_prec+(toff-ton)
      END IF

      ton=toff
      starttime=ton
      skston=ton
      CALL ForwardSolve
      CALL second0(toff)
      endtime=toff
      skstoff=toff
      block_factorization_time=block_factorization_time+(skstoff-skston)
      usedtime=endtime-starttime

      time_factor_blocks = time_factor_blocks + (toff-ton)

      HESSPASS = HESSPASS+1
#endif
      END SUBROUTINE Compute_Hessian_Blocks_With_Col_Redist

      SUBROUTINE Compute_Hessian_Blocks_Thomas (xc, gc, func,           &
                 l_Hess_sym)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION((1+mpol)*(2*ntor+1)*3,ns),                  &
               INTENT(INOUT)  :: xc, gc
      LOGICAL, INTENT(IN) :: l_Hess_sym
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!----------------------------------------------
      REAL(rprec), PARAMETER    :: zero=0, one=1, p5=0.5_dp
      INTEGER  :: n, m, js, js1, mesh, ntype, istat, iunit, icol, icolmpi
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:)  :: gc1, gcswap
      REAL(rprec) :: eps, ton, toff, bsize, lmdamp, skston, skstoff, temp
      CHARACTER(LEN=3) :: label
      EXTERNAL :: func
#if defined(MPI_OPT)
      INTEGER  :: neqtotal, numroc
      REAL(rprec), ALLOCATABLE, DIMENSION (:,:) :: idsendBuf(:,:)
      REAL(dp) :: starttime, endtime, usedtime
      REAL(dp) :: colstarttime, colendtime, colusedtime
      EXTERNAL :: numroc
!------------------------------------------------------------

      starttime=0; endtime=0; usedtime=0
      colstarttime=0; colendtime=0; colusedtime=0
#endif


!------------------------------------------------------------
!     COMPUTE (U)pper, (D)iagonal, and (L)ower Block matrices
!------------------------------------------------------------
      IF (iam .eq. 0) THEN
         IF (.not. l_diagonal_only) WRITE (6, 50)
      END IF
 50   FORMAT (/,' Computing block preconditioner')
      
      CALL second0(ton)

      mblk_size = SIZE(gc,1)

      eps = SQRT(EPSILON(eps))*ABS(eps_factor)

      IF (ALLOCATED(ublk)) DEALLOCATE(ublk, dblk, lblk, stat=istat)
      IF (LSCALAPACK) THEN
#if defined(MPI_OPT)
        CALL blacs_gridinfo(icontxt_1xp,nprow,npcol,myrow,mycol)
        mb = mblk_size
        nb = 1

        rsrc = 0
        csrc = 0
        Locq = numroc( mblk_size, nb, mycol, csrc, npcol )
        Locp = numroc( mblk_size, mb, myrow, rsrc, nprow )
        mblk_size2=max(1,Locq)
        lld = max(1,Locp)
        call descinit(descA_1xp,mblk_size,mblk_size,mb,nb,rsrc,csrc,       &
                      icontxt_1xp,lld,info)
        if (info.ne.0) then
          write(*,*) 'myrow,mycol,nprow,npcol,desc(LLD_),lld ',            &
             myrow,mycol,nprow,npcol,descA_1xp(LLD_),lld
          write(*,*) 'Locp,m,mb ', Locp,mblk_size,mb
        endif

        call assert(info.eq.0,'descinit descA_1xp ',info)
        ineed = max(1,lld*Locq)
        mm = mblk_size*ns
        mb=mm
        nb=nrhs1
        LocqR = numroc( nrhs1, nb, mycol,csrc,npcol)
        LocpR = numroc( mm, mb, myrow,rsrc,nprow)
        lld = max(1,LocpR)
        call descinit(descR_all,mm,nrhs1, mb,nb,rsrc,csrc,icontxt,lld,info)
        call assert(info.eq.0,'test_pdtrd:descinit return info=',info)

        call blacs_gridinfo(icontxt,nprow,npcol,myrow,mycol)

        mb = 50
        nb = mb

        rsrc = 0
        csrc = 0
        Locq = numroc( mblk_size, nb, mycol, csrc, npcol )
        Locp = numroc( mblk_size, mb, myrow, rsrc, nprow )

        lld = MAX(1,Locp)
        call descinit(descA,mblk_size,mblk_size,mb,nb,rsrc,csrc,icontxt,lld,info)
        call assert(info.eq.0,'descinit descA ',info)
        ineed = MAX(1,lld*Locq)
        IF (ALLOCATED(ublkp)) DEALLOCATE(ublkp, dblkp, lblkp, stat=istat)
        IF (.NOT. ALLOCATED(ublkp)) THEN
           ALLOCATE(ublkp(ineed,ns-1),     &
                    dblkp(ineed,ns),       &
                    lblkp(ineed,ns-1),     &
                    stat=istat)
           IF (istat .ne. 0) THEN
             PRINT *,' Not enough memory to store a single block!'
             STOP
           END IF
        END IF

        ALLOCATE(ublk(mblk_size,mblk_size2,ns),    &
                 dblk(mblk_size,mblk_size2,ns),    &
                 lblk(mblk_size,mblk_size2,ns),    &
                 stat=istat)
        IF (istat .ne. 0) THEN
           PRINT *,' Not enough memory to store a single block!'
           STOP
        END IF
#else
        STOP 'LSCALAPACK=T BUT NOT MPI!'
#endif
      ELSE   !.NOT.LSCALAPACK

        ALLOCATE(ublk(mblk_size,mblk_size,ns),        &
                 dblk(mblk_size,mblk_size,ns),        &
                 lblk(mblk_size,mblk_size,ns),        &
                 stat=istat)

        lswap2disk = (istat .ne. 0)

        IF (lswap2disk .and. .not.l_Diagonal_Only) THEN
          IF (iprec_flag == 0) THEN
            PRINT *,'Not enough memory to store Hessian blocks in memory.'
            PRINT *,'Writing blocks to disk file.'
          END IF

          CALL RemoveDAFile

          ScratchFile = "PRCND2A.bin"
          IF (FlashDrive .ne. "") ScratchFile = FlashDrive // ScratchFile
          CALL OpenDAFile(mblk_size, mblk_size**2, 3, ScratchFile,       &
                          iunit_dacess, CreateNew) 
          ScratchFile = "PRCND2B.bin"
          IF (FlashDrive .ne. "") ScratchFile = FlashDrive // ScratchFile

        ELSE
          ScratchFile = ""
        END IF
      
      END IF       !END LSCALAPACK

      IF (ALLOCATED(ublk)) THEN
         ublk = 0; dblk = 0; lblk = 0
      END IF

      ALLOCATE (DataItem(mblk_size), stat=istat)
      IF (istat .ne. 0) STOP 'DataItem allocation failed!'
      DataItem = 0

      IF (iprec_flag == 0) THEN
         bsize = 3*ns*KIND(dblk)
         bsize = bsize*REAL(mblk_size*mblk_size,dp) 
         IF (bsize .lt. 1.E6_dp) THEN
            bsize = bsize/1.E3
            label = " Kb"
         ELSE IF (bsize .lt. 1.E9_dp) THEN
            bsize = bsize/1.E6_dp
            label = " Mb"
         ELSE
            bsize = bsize/1.E9_dp
            label = " Gb"
         END IF

         IF (iam .eq. 0) THEN
            DO iunit = 6,33,27
            WRITE (iunit, '(1x,a,i4,a,f6.2,a)') 'Block dim: ', mblk_size, &
     &         '^2  Preconditioner size: ', bsize, TRIM(label)
            END DO
         ENDIF 
      
         iprec_flag = 1

      END IF

      xc = 0

#if defined(MPI_OPT)
      icolmpi = 0
#endif
      icol=0

      CALL second0(skston)
      CALL func
      CALL second0(skstoff)
#if defined(SKS)
      hessian_funct_island_time=hessian_funct_island_time+(skstoff-skston)
#endif
#if defined(MPI_OPT)
      CALL second0(colstarttime) !Added by SKS for timing comparisons
#endif
      VAR_TYPE: DO ntype = 1, 3
         VAR_N: DO n = -ntor, ntor
            VAR_M: DO m = 0, mpol
                  
               icol=icol+1
#if defined(MPI_OPT)
               IF(MOD(icol-1,nprocs)==iam .OR. .NOT.LSCALAPACK) THEN
                  IF (LSCALAPACK) THEN
                     icolmpi = icolmpi+1
                  ELSE 
                     icolmpi = icol
                  END IF
#else
                  icolmpi = icol
#endif
                  MESH_3PT: DO mesh = 1, 3

                  xc = 0

                  DO js = jstart(mesh), ns, 3
                     xc(icol,js) = eps
                  END DO

                  INHESSIAN=.TRUE.
                  CALL second0(skston)
                  CALL func
                  CALL second0(skstoff)
#if defined(SKS)
                  hessian_funct_island_time=hessian_funct_island_time+(skstoff-skston)
#endif
                  INHESSIAN=.FALSE.
#if defined(MPI_OPT)
#endif

!OFF FOR l_linearize=T                  gc = gc-gc1

!              COMPUTE PRECONDITIONER (HESSIAN) ELEMENTS. LINEARIZED EQUATIONS
!              OF TRI-DIAGONAL (IN S) FORM 
!
!              F(j-1) = a(j-1)x(j-2) + d(j-1)x(j-1) + b(j-1)x(j)
!              F(j)   =                a(j)x(j-1)   + d(j)  x(j)  + b(j)  x(j+1)
!              F(j+1) =                               a(j+1)x(j)  + d(j+1)x(j+1) + b(j+1)x(j+2)
!
!              HESSIAN IS H(k,j) == dF(k)/dx(j); aj == lblk; dj == dblk; bj = ublk

                  SKIP3_MESH: DO js = jstart(mesh), ns, 3

!FORCE RESPONSE AT mp,np,nptype TO m,n,js,ntype VELOCITY PERTURBATION
                     !Levenberg-parameter: lv-param*diagonal(js=4)
                     !diag(1)=0 at this point, can't use it for scaling!
                     lmdamp = levmarq_param*gc(icol,js)/eps

!Avoid possible singularity when evolving both s and u at origin
!!SPH - MUCKS UP CONVERGENCE WHEN lmdamp = 0
!!SPH                     IF (lmdamp.eq.zero .and. js.eq.1 .and. m.eq.1      &
!!SPH     &                    .and. ntype.eq.1) lmdamp = 1.E-6_dp*gc(icol,js)/eps

!NO SCALING          lmdamp = levmarq_param
                     lmdamp = ABS(lmdamp)

                     IF (levmarq_param .ne. zero) THEN
                        IF (lmdamp .eq. zero) THEN
                           lmdamp = levmarq_param
                        ELSE IF (ntype .ne. 3) THEN
                           lmdamp = lmdamp/10                          !Emulate v|| damping only
                        END IF
                     END IF

!                    WRITE (6, 120) lmdamp,js,m,n,ntype
!120 FORMAT('lmdamp=',1p,e10.2,' js: ',i4,' m: ',i4,' n: ',i4,' ntype: ',i4)
 
                     !ublk(js-1)
                     js1 = js-1
                     IF (js1 .GT. 0) THEN
                        DataItem = gc(:,js1)/eps
!                     IF (l_Viscous) DataItem(icol) = DataItem(icol)-p5*lmdamp

                     !m=0 constraint for n<0
                        IF (m.eq.0 .and. n.lt.0) DataItem = 0
                        IF (l_Diagonal_Only) THEN
                           temp=DataItem(icol)
                           DataItem=0
                           DataItem(icol)=temp
                        END IF
                        IF (lswap2disk) THEN
                           CALL WriteDAItem_SEQ(DataItem)
                        ELSE 
                           ublk(:,icolmpi,js1) = DataItem
                        END IF
                     END IF

                     !dblk(js)
                     DataItem = gc(:,js)/eps
                     IF (m.eq.0 .and. n.lt.0) DataItem = 0

!                    Boundary condition at js=1 and ns (CATCH ALL OF THEM HERE)
!                    and m=0,n<0: NEED THIS to avoid (near) singular Hessian
!                    ASSUMES DIAGONALS ARE ALL NEGATIVE

                     IF (ALL(DataItem .eq. zero)) THEN
                        DataItem(icol) = -one
                     END IF

!                    IF (ntype.eq.3 .and. js.eq.ns) lmdamp = levmarq_param0/10

                     DataItem(icol) = DataItem(icol) + SIGN(lmdamp,DataItem(icol))

                     IF (l_Diagonal_Only) THEN
                        temp=DataItem(icol)
                        DataItem=0
                        DataItem(icol)=temp
                     END IF

                     IF (lswap2disk) THEN
                        CALL WriteDAItem_SEQ(DataItem)
                     ELSE 
                        dblk(:,icolmpi,js) = DataItem
                     END IF

                     !lblk(js+1) 
                     js1=js+1
                     IF (js .LT. ns) THEN
                        DataItem = gc(:,js1)/eps
!                     IF (l_Viscous) DataItem(icol) = DataItem(icol)-p5*lmdamp
                     !m=0 constraint for n<0
                        IF (m.eq.0 .and. n.lt.0) DataItem = 0
                        IF (l_Diagonal_Only) THEN
                           temp=DataItem(icol)
                           DataItem=0
                           DataItem(icol)=temp
                        ENDIF
                        IF (lswap2disk) THEN
                           CALL WriteDAItem_SEQ(DataItem)
                        ELSE 
                           lblk(:,icolmpi,js1) = DataItem
                        END IF
                     END IF

                  END DO SKIP3_MESH
               END DO MESH_3PT
#if defined(MPI_OPT)
               ENDIF
#endif
            END DO VAR_M
         END DO VAR_N
      END DO VAR_TYPE
      
      CALL second0(toff)
      
      IF (l_Diagonal_Only) THEN
         time_diag_prec = time_diag_prec+(toff-ton)
      ELSE
         time_block_prec = time_block_prec+(toff-ton)
      END IF

      ton = toff

#if defined(SKS)
      CALL second0(colendtime)
      colusedtime=colendtime-colstarttime
      time_generate_blocks=time_generate_blocks+colusedtime
      construct_hessian_time=construct_hessian_time+(colendtime-colstarttime)
#endif

!DUMPS FULL HESSIAN BLOCKS FOR TESTING (ONLY WORKS FOR NON-MPI)
      IF (.FALSE.) THEN
!       IF (.not. l_Diagonal_Only) THEN
!      IF (iam.eq.0 .and. nprecon .eq. 2) THEN
         OPEN(UNIT=20,FILE='LBLK.txt',STATUS='REPLACE',FORM='FORMATTED')
         WRITE (20, *) mblk_size, ns
         WRITE (20, *) lblk
         CLOSE (20)

         OPEN(UNIT=20,FILE='DBLK.txt',STATUS='REPLACE',FORM='FORMATTED')
         WRITE (20, *) mblk_size, ns
         WRITE (20, *) dblk
         CLOSE (20)

         OPEN(UNIT=20,FILE='UBLK.txt',STATUS='REPLACE',FORM='FORMATTED')
         WRITE (20, *) mblk_size, ns
         WRITE (20, *) ublk
         CLOSE (20)

         OPEN(UNIT=20,FILE='RHS.txt',STATUS='REPLACE',FORM='FORMATTED')
	     WRITE (20, *) gc

         PRINT *, 'HESSIAN BLOCKS WRITTEN OUT'

      END IF

      DEALLOCATE (DataItem, stat=istat)
!      DEALLOCATE (gc1, DataItem, stat=istat)
      IF (istat .NE. 0) STOP 'DataItem deallocation error!'

      IF (lswap2disk) GO TO 900
!
!     CHECK SYMMETRY OF BLOCKS 
!
      CALL second0(skston)
      IF (LSCALAPACK) THEN
         CALL check3d_symmetry(dblk,lblk,ublk,dblkp,lblkp,ublkp,        &
                               mblk_size,mblk_size2,ns,asym_index)
      ELSE
         CALL check3d_symmetry_serial(dblk,lblk,ublk,mpol,ntor,         &
                               mblk_size,ns,asym_index)
      END IF
      CALL second0(skstoff)
#if defined(SKS)
      asymmetry_check_time=asymmetry_check_time+(skstoff-skston)
#endif

 900  CONTINUE

      IF (iam .eq. 0) THEN
         WRITE (6, 100) levmarq_param, muPar, asym_index
         IF (.not. l_diagonal_only)                                      &
         WRITE (33,110) levmarq_param, muPar, asym_index
      ENDIF

 100  FORMAT(' LM parameter: ',1pe9.2,                                 &
             ' mu||: ',1pe9.2, ' Asym index: ',1pe9.2)
 110  FORMAT (/,' Computing block preconditioner: ',                   &
              ' LM parameter: ',1pe9.2,' mu||: ',1pe9.2,               &
              ' Asym index: ',1pe9.2)

!
!     SYMMETRIZE BLOCKS
!
      IF (.not. l_Hess_sym .or. lswap2disk) GOTO 1000
      icol = 0
      ALLOCATE(gc1(mblk_size,ns))
      
      DO m = 0, mpol
         DO n = -ntor, ntor
            DO ntype = 1,3
               icol = icol+1
               gc1(:,1:nsh) = (ublk(icol,:,1:nsh)          &
     &                      +  lblk(:,icol,2:ns))/2
               ublk(icol,:,1:nsh) = gc1(:,1:nsh)
               lblk(:,icol,2:ns)  = gc1(:,1:nsh)
               gc1(:,:) = (dblk(icol,:,:)                  &
     &                 +  dblk(:,icol,:))/2
               dblk(icol,:,:) = gc1(:,:)
               dblk(:,icol,:) = gc1(:,:)
            END DO
         END DO
      END DO

      DEALLOCATE (gc1)
 1000 CONTINUE

!      IF (l_Diagonal_Only) RETURN


!      IF (nprecon .eq. 4) l_backslv = .TRUE.
      IF (l_backslv) THEN
         ALLOCATE(ublk_s(mblk_size,mblk_size,ns),    &
                  dblk_s(mblk_size,mblk_size,ns),    &
                  lblk_s(mblk_size,mblk_size,ns),    &   
                  stat=istat)
         IF (istat .ne. 0) STOP 'Allocation error2 in Compute_Hessian_Blocks'
         ublk_s = ublk;  dblk_s = dblk; lblk_s = lblk
      END IF
!
!     FACTORIZE (Invert) HESSIAN
!
      IF (ALLOCATED(ipiv_blk)) DEALLOCATE(ipiv_blk,stat=istat)
      CALL second0(skston)
      IF (LSCALAPACK) THEN
#if defined(MPI_OPT)
         ineed = numroc(descA(M_),descA(MB_),0,0,nprow) + descA(MB_)
         ALLOCATE( ipiv_blk(ineed, ns), stat=istat )
         IF (istat .ne. 0) STOP 'ipiv_blk allocation error2 in Compute_Hessian_Blocks'
         CALL blk3d_parallel(dblk,lblk,ublk,dblkp,lblkp,ublkp,mblk_size,mblk_size2,ns)
         CALL blacs_gridinfo(icontxt,nprow,npcol,myrow,mycol)
         mm = mblk_size*ns
         LocqR = numroc( nrhs1, nb, mycol,csrc,npcol)
         LocpR = numroc( mm, mb, myrow,rsrc,nprow)
         lld = MAX(1,LocpR)
         CALL descinit(descX,mm,nrhs1, mb,nb,rsrc,csrc,icontxt,lld,info)
         CALL assert(info.eq.0,'test_pdtrd:descinit return info=',info)
         descR(:) = descX

         ineedR = MAX(1,lld*LocqR)
 
         CALL blacs_barrier(icontxt,'All')
         CALL profstart('factor call:123')
         CALL pdtrdf(mblk_size,ns,dblkp,lblkp,ublkp,ipiv_blk,descA)
         CALL blacs_barrier(icontxt,'All')
         CALL profend('factor call:123')
#else
         STOP 'LSCALAPACK=T BUT NOT IN MPI!'
#endif
      ELSE    !.NOT. LSCALAPACK
         ALLOCATE (ipiv_blk(mblk_size,ns), stat=istat)
         IF (istat .ne. 0) STOP 'ipiv_blk allocation error2 in Compute_Hessian_Blocks'
         IF (lswap2disk) THEN
            CALL blk3d_factor_swp(ipiv_blk, mblk_size, ns)
         ELSE
            CALL blk3d_factor(dblk, lblk, ublk, ipiv_blk, mblk_size, ns)
         END IF
      END IF  ! LSCALAPACK
#if defined(SKS)
      CALL second0(skstoff)
      block_factorization_time=block_factorization_time+(skstoff-skston)
#endif

      CALL second0(toff)
      time_factor_blocks = time_factor_blocks + (toff-ton)

      END SUBROUTINE Compute_Hessian_Blocks_Thomas

!!DEC$ DEFINE SVD_ON
!!DEC$ UNDEFINE SVD_ON

      SUBROUTINE blk3d_factor(a, bm1, bp1, ipiv, mblk, nblocks)
      USE stel_constants, ONLY: one, zero
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER  :: bytes_per_rprec2 = 8
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(out) :: ipiv(mblk,nblocks)
      REAL(dp), TARGET, DIMENSION(mblk,mblk,nblocks), INTENT(inout) :: &
     &                       a, bm1, bp1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!      INTEGER :: ibuph, incnow, irecl, incbu, iunit=102, ndisk
      INTEGER :: k, k1, ier
      INTEGER, POINTER :: ipivot(:)
      REAL(dp), POINTER :: amat(:,:), bmat(:,:), cmat(:,:)
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: temp
#if defined(SVD_ON)
      REAL(dp), ALLOCATABLE, DIMENSION(:)   :: w, work
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: v, u
      REAL(dp), PARAMETER :: small = 1.E-300_dp
      INTEGER :: nw
      INTEGER :: lwork
      CHARACTER(LEN=1), PARAMETER :: jobu='A', jobvt='A'
#endif
!-----------------------------------------------
!  modified (June, 2003, ORNL):         S. P. Hirshman
!-----------------------------------------------------------------------
!
!  this subroutine solves for the Q factors of a block-tridiagonal system of equations.
!
!-----------------------------------------------------------------------
!  INPUT
!  mblk                : block size
!  nblocks             : number of blocks
!  a                   : diagonal blocks
!  bp1, bm1            : lower, upper blocks (see equation below)
!
!  OUTPUT
!  ipiv                : pivot elements for kth block
!  a                   : a-1 LU factor blocks
!  bm1                 : q = a-1 * bm1 matrix
!
!  LOCAL VARIABLES
!  iunit               : unit number for block-tridiagonal solution disk file.
!
!  solutions are indexed in m-n fourier-space, legendre-space. the tri-diagonal
!  equation is:
!
!           bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) = source(l)
!
!     GENERAL SOLUTION SCHEME APPLIED TO EACH BLOCK ROW (INDEX L)
!
!     1. Start from row N and solve for x(N) in terms of x(N-1):
!
!        x(N) = -q(N)*x(N-1) + r(N)
!
!        q(N) =  a(N)[-1] * bm1;    r(N) = a(N)[-1] * s(N)
!
!        where a(N)[-1] is the inverse of a(N)
!
!     2. Substitute for lth row to get recursion equation for q(l) and r(l):
!
!        x(l) = -q(l)*x(l-1) + r(l), in general, where:
!
!        q(l) = (a(l) - bp1(l)*q(l+1))[-1] * bm1(l)
!
!        qblk(l) == (a(l) - bp1(l) * q(l+1))[-1] on return
!
!        r(l) = (a(l) - bp1(l)*q(l+1))[-1] * (s(l) - bp1(l)*r(l+1))
!
!     3. At row l = 1, bm1(1) = 0 and get an equation for x(1) corresponding to q(1) = 0:
!
!        x(1) = r(1)
!
!     4. Finally, can back-solve for x(N) in terms of x(N-1) from eqn.(1) above
!
!
!     NUMERICAL IMPLEMENTATION (USING LAPACK ROUTINES)
!
!     IF USING SVD:
!
!     1. CALL dgesvd:   Perform SVD decomposition
!     2. CALL svbksb:   Solve Ax = bi, for i=1,mblk
!
!     OTHERWISE
!
!     1. CALL dgetrf:   Perform LU factorization of diagonal block (A) - faster than sgefa
!     2. CALL dgetrs:   With multiple (mblk) right-hand sides, to do block inversion
!                         operation, Ax = B  (stores result in B; here B is a matrix)
!
!      ndisk = mblk

!  create disk file for doing direct access i/o.

!      incnow = ndisk
!      irecl  = bytes_per_rprec2*incnow
!      incbu = 1 + (ndisk - 1)/incnow
!      ibuph = 0

!      iunit = 10
!      CALL safe_open(iunit, ier, 'NULL', 'scratch', 'unformatted',
!     1     irecl, 'DIRECT')
!      IF (ier .ne. 0) STOP 'Error opening scratch file in blk3d'

!  main loop. load and process (backwards) block-rows nblocks to 1. 

#if defined(SVD_ON)
      lwork = 10*mblk
      ALLOCATE(v(mblk,mblk), w(mblk),   &
     &         u(mblk,mblk), work(lwork), stat=ier)
      IF (ier .ne. 0) STOP 'Allocation error in blk3d_factor!'
#endif
      ipiv = 0
      ALLOCATE (temp(mblk,mblk), stat=ier)
      IF (ier .ne. 0) STOP 'Allocation error in blk3d_factor!'

      BLOCKS: DO k = nblocks, 1, -1
!
!     Compute (and save) qblk(nblocks) = ablk(nblocks)[-1] * bml
!
         amat => a(:,:,k)

#if defined(SVD_ON)
!NOTE: v = vt coming out here
         CALL dgesvd (jobu, jobvt, mblk, mblk, amat, mblk, w, u, mblk,  &
        &             v, mblk, work, lwork, ier)
!Set SVD weights w to 0 for weights below the allowed threshold, so backsolver
!will compute pseudo-inverse
         WRITE (35, '(a,i4,a,1pe12.4)') 'Block: ',k, ' Condition #: ', SQRT(w(1)/w(mblk))
         DO nw = 2, mblk
            IF (w(nw) .le. small*w(1)) w(nw) = 0
         END DO
!
!        STORE svd pseudo-inverse IN AMAT since u,v,w will be deallocated at end
         CALL svdinv2 (amat, u, v, w, mblk)
#else
         ipivot => ipiv(:,k)
         CALL dgetrf (mblk, mblk, amat, mblk, ipivot, ier)
#endif
         IF (ier .ne. 0) GOTO 200
         IF (k .eq. 1) EXIT

         bmat => bm1(:,:,k)

#if defined(SVD_ON)
!
!        Use (pseudo) inverse stored in AMAT = V*1/w*Ut 
!
         temp = bmat
         bmat = MATMUL(amat, temp)
#else
         CALL dgetrs('n', mblk, mblk, amat, mblk, ipivot,               &
     &                bmat, mblk, ier)
         IF (ier .ne. 0) GOTO 305
#endif
!         CALL wrdisk(iunit, ql, ndisk, incnow, ibuph, incbu, ier)
!         IF (ier .ne. 0) GOTO 302

!
!      Update effective diagonal "a" matrix. Use dgemm: faster AND doesn't overflow normal stack
!
         k1 = k-1 
         amat => bp1(:,:,k1)
         cmat => a(:,:,k1)
!         cmat = cmat - MATMUL(amat, bmat)
         CALL dgemm('N','N',mblk,mblk,mblk,-one,amat,mblk,              &
     &              bmat, mblk, one, cmat, mblk)

      END DO BLOCKS

#if defined(SVD_ON)
      DEALLOCATE(v, w, u, work)
#endif
!
!     COMPUTE TRANSPOSES HERE, SINCE REPEATEDLY CALLING MATMUL OPERATION
!     X*At IS FASTER THAN A*X DUE TO UNIT STRIDE
!

      DO k = 1, nblocks
         IF (k .ne. nblocks) THEN
            temp = TRANSPOSE(bp1(:,:,k))
            bp1(:,:,k) = temp
         END IF
         IF (k .ne. 1) THEN
            temp = TRANSPOSE(bm1(:,:,k))
            bm1(:,:,k) = temp
         END IF
      END DO

      DEALLOCATE (temp)

      GOTO 400

!  error returns. ------------------------------------------------------

  200 CONTINUE
!          < 0:  if info = -i, the i-th argument had an illegal value
!          > 0:  if info = i, u(i,i) is exactly zero. the factorization
      WRITE (6, '(2x,a,i4)') 'Error factoring matrix in blk3d: block = '   &
     &                        , k
      IF (ier < 0)WRITE (6,'(i4, a)') ier, 'th argument has illegal value'
      IF (ier > 0)WRITE (6,'(i4, a)') ier, 'th diagonal factor exactly zero'
      STOP
  301 CONTINUE
!      WRITE (6, '(a,i8)') ' BLK3D:   error in opening file:  ',
!     1   'RECL = ', irecl
  302 CONTINUE
      WRITE (6, '(a)') ' BLK3D:   error in I/O routine WRDISK'
  303 CONTINUE
      WRITE (6, '(a)') ' BLK3D:   error in I/O routine RDDISK'
      ier = -2
  305 CONTINUE
      WRITE (6, '(2/a,i10,2/)') ' BLK3D-FACTOR:   error detected:   ier =',     &
     &   ier
      STOP

!  destroy disk file and return. ---------------------------------------
 
  400 CONTINUE

!      CLOSE (iunit)

      END SUBROUTINE blk3d_factor


      SUBROUTINE blk3d_factor_swp(ipiv, mblk, nblocks)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1    
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(out) :: ipiv(mblk,nblocks)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) ::                       &
     &                       amat, bmat, cmat, temp
      INTEGER :: k, k1, ier
      INTEGER, POINTER :: ipivot(:)
!-----------------------------------------------
!  modified (June, 2003, ORNL):         S. P. Hirshman
!  modified (June, 2007, ORNL), added lswap2disk logic
!-----------------------------------------------------------------------
!
!  this subroutine solves for the Q factors of a block-tridiagonal system of equations.
!  see blk3d_factor for more details
!
!-----------------------------------------------------------------------
!  INPUT
!  mblk                : block dimension (elements in a block=mblk X mblk)
!  nblocks             : number of blocks
!
!  OUTPUT
!  ipiv                : pivot elements for kth block
!
!  LOCAL VARIABLES
!  iunit               : unit number for block-tridiagonal solution disk file.
!
!  solutions are indexed in m-n fourier-space, legendre-space. the tri-diagonal
!  equation is:
!
!           bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) = source(l)
!
!
!     NUMERICAL IMPLEMENTATION (USING LAPACK ROUTINES)
!
!     1. CALL dgetrf:   Perform LU factorization of diagonal block (A) - faster than sgefa
!     2. CALL dgetrs:   With multiple (mblk) right-hand sides, to do block inversion
!                       operation, A X = B  (stores result in B; here B is a matrix)
!

!  main loop. load and process (backwards) block-rows nblocks to 1. 

      ipiv = 0

!     CHANGE Direct Access Record length to block size (from individual rows)
      CALL ChangeDAFileParams(mblk**2, mblk**2, 3, ScratchFile, nblocks)

      ALLOCATE(amat(mblk,mblk), bmat(mblk,mblk), cmat(mblk,mblk),       &
     &         temp(mblk,mblk), stat=ier)
      IF (ier .ne. 0) STOP 'Allocation error in blk3d_factor_swp!' 

      CALL ReadDAItem2(temp, nblocks, bldia)

      BLOCKS: DO k = nblocks, 1, -1
!
!     Compute (and save) qblk(k) = ablk(k)[-1] * bml
!
         amat = temp
         ipivot => ipiv(:,k)
         CALL dgetrf (mblk, mblk, amat, mblk, ipivot, ier)
         IF (ier .ne. 0) GOTO 200
!CONFIRM READ-WRITE ALLOWED...OK FOR DA Files!
         CALL WriteDAItem_RA(amat, k, bldia, 1)

         IF (k .eq. 1) EXIT
          
         CALL ReadDAItem2(bmat, k, blmin)
         CALL dgetrs('n', mblk, mblk, amat, mblk, ipivot, bmat, mblk, ier)
         IF (ier .ne. 0) GOTO 305
!
!     COMPUTE TRANSPOSES HERE (and for cmat below), SINCE REPEATEDLY CALLING MATMUL OPERATION
!     X*At IS FASTER THAN A*X DUE TO UNIT STRIDE
!
         temp = TRANSPOSE(bmat)
         CALL WriteDAItem_RA(temp, k, blmin, 1)

!
!      Update effective diagonal "a" matrix. Use dgemm: faster AND doesn't overflow normal stack
!
         k1 = k-1 
         CALL ReadDAItem2(amat, k1, blpls)
         CALL ReadDAItem2(temp, k1, bldia)
!         temp = temp - MATMUL(amat, bmat)
         CALL dgemm('N','N', mblk, mblk, mblk, -one, amat, mblk,        &
     &              bmat, mblk, one, temp, mblk)
         cmat = TRANSPOSE(amat)
         CALL WriteDAItem_RA(cmat, k1, blpls, 1)

      END DO BLOCKS

      GOTO 400

!  error returns. ------------------------------------------------------

  200 CONTINUE
!          < 0:  if info = -i, the i-th argument had an illegal value
!          > 0:  if info = i, u(i,i) is exactly zero. the factorization
      WRITE (6, '(2x,a,i4)') 'Error factoring matrix in blk3d: block = ' &
     &                        , k
      IF (ier < 0)WRITE (6,'(i4, a)') ier, 'th argument has illegal value'
      IF (ier > 0)WRITE (6,'(i4, a)') ier, 'th diagonal factor exactly zero'
      STOP
  305 CONTINUE
      WRITE (6, '(2/a,i10,2/)') ' BLK3D-FACTOR-SWP:   error detected:   ier =', ier
      STOP

  400 CONTINUE

      DEALLOCATE(amat, bmat, cmat, temp, stat=ier)

      CALL CloseDAFile

      END SUBROUTINE blk3d_factor_swp


      SUBROUTINE block_precond(gc)
      USE timer_mod
      IMPLICIT NONE        
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(mblk_size,ns), INTENT(inout) :: gc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER  ::  zero = 0
      INTEGER :: m, n, js, ntype, istat, icol
      REAL(dp) :: t1, error, error_sum, b_sum
      REAL(dp), ALLOCATABLE :: gc_s(:,:)
#if defined(MPI_OPT)
      REAL(dp) :: alpha0,beta0
      INTEGER  :: i,j,k,nn,kk,ia,ja,ix,jx,ic,jc,ib,jb
      INTEGER :: MPI_ERR
#endif
#if defined(SKS)
      REAL(rprec) :: starttime, endtime, usedtime
      INTEGER :: PRECONDPASS
      SAVE PRECONDPASS
      DATA PRECONDPASS /0/
      INTEGER :: globrow, cnt
      !Might need to make it allocatable - just a test right now.
      REAL(dp), DIMENSION (:,:), ALLOCATABLE :: solvec
!-----------------------------------------------
      starttime=0
      endtime=0
#endif
!-----------------------------------------------
!
!     Applies 3D block-preconditioner to forces vector (gc)
!
#if defined(SKS)
      IF (PARSOLVER) THEN

        DO globrow=myStartBlock,myEndBlock
          CALL SetMatrixRHS (globrow,gc(1:mblk_size,globrow))
        END DO
        
        CALL BackwardSolve
       
        ALLOCATE (solvec(mblk_size,myNumBlocks), stat=istat)
        IF (istat .ne. 0) STOP 'Allocation error in block_precond before gather'

        cnt=1
        DO globrow=myStartBlock,myEndBlock
          CALL GetSolutionVector (globrow,solvec(1:mblk_size,cnt))
          cnt = cnt + 1
        END DO

        CALL MPI_Allgatherv(solvec,myNumBlocks*mblk_size,MPI_REAL8,gc, &
             rcounts, disp, MPI_REAL8, MPI_COMM_WORLD, MPI_ERR)

!!        CALL GetFullSolution(temp,gc)
        
        DEALLOCATE (solvec)

      ELSE

        CALL second0(starttime)
#endif
         IF (l_backslv) THEN
            ALLOCATE (gc_s(mblk_size,ns), stat=istat)
            IF (istat .ne. 0) STOP 'Allocation error0 in block_precond'
            gc_s = gc
         END IF

!     Solve Hx = gc, using Hessian (H) LU factors stored in block_... matrices
!     Store solution "x" in gc on output
         IF (LSCALAPACK) THEN
#if defined(MPI_OPT)
         ALLOCATE (tempp(ineedR), stat=istat)
         ia = 1
         ja = 1
         ib = 1
         jb = 1
         beta0 = 0.0d0
         alpha0 = 1.0d0
         mm0=descR(3)
         nrhs0=descR(4)
         call profstart('pdgeadd call:123')
         call pdgeadd( 'N', mm0,nrhs0, alpha0,gc,ia,ja,descR_all,       &
                       beta0,tempp,ib,jb,descR )
         call profend('pdgeadd call:123')
         ir = 1
         jr = 1
         call blacs_barrier(icontxt,'All')
         call profstart('solver call:123')
         call pdtrds(mblk_size,ns,dblkp,lblkp,ublkp,ipiv_blk,descA,     &
                     nrhs1,tempp,ir,jr,descR )
         call blacs_barrier(icontxt,'All')
         call profend('solver call:123')
         gc = 0
         ia = 1
         ja = 1
         ib = 1
         jb = 1
         beta0 = 0.0d0
         alpha0 = 1.0d0
         mm0=descR(3)
         nrhs0=descR(4)
!     descR_all(7)=-1
!     descR_all(8)=-1
         call profstart('pdgeadd2 call:123')
         call pdgeadd( 'N', mm0,nrhs0, alpha0, tempp,ia,ja,descR,       &
                       beta0,gc,ib,jb,descR_all )
         call dgsum2d( icontxt,'A', ' ', mm,1,gc,mm,       -1,-1)
         call profend('pdgeadd2 call:123')

         DEALLOCATE (tempp, stat=istat)
#else
         STOP 'MPI_OPT must be true for LSCALAPACK!'
#endif
         ELSE       !NOT LSCALAPACK

#if defined(SKS)
         CALL second0(endtime)
         usedtime=endtime-starttime
         IF (SKSDBG) WRITE(TOFU,*) 'HessianTag-2 : Time to BackwardSolve:',usedtime,'PRECONDPASS',PRECONDPASS,'in native run'
         IF (SKSDBG) CALL FLUSH(TOFU)
#endif
!        Serial solver
         IF (lswap2disk) THEN
            CALL blk3d_slv_swp(gc, ipiv_blk, mblk_size, ns)
         ELSE
            CALL blk3d_slv(dblk, lblk, ublk, gc, ipiv_blk, mblk_size, ns)
         END IF
         END IF     !END LSCALAPACK

#if defined(SKS)
      END IF !END OF IF(PARSOLVER) CONDITIONAL
      PRECONDPASS = PRECONDPASS + 1
#endif

#if !defined(MPI_OPT)
      IF (l_backslv) THEN
         error_sum = 0;  b_sum = SUM(gc_s*gc_s)
         l_backslv = .false.
         icol = 0

         WRITE (6, *) ' Writing block Hessian check to unit 34'
         WRITE (34, *)
         WRITE (34, *) ' BLK3D FACTORIZATION CHECK: Ax = b ?'
        
         DO n = -ntor, ntor
            WRITE (34, *) ' N = ', n
            DO m = 0, mpol
               WRITE (34, *) ' M = ', m
               DO ntype = 1, 3
                  icol = icol+1
                  WRITE (34, *) ' TYPE = ', ntype
                  WRITE (34, *)                                        &
                  '   js        x             Ax             b     '// &
                  '     Ax - b        RelErr'
                  js = 1
                  t1 = SUM(dblk_s(icol,:,js)*gc(:,js)                  &
                     + ublk_s(icol,:,js)*gc(:,js+1))

                  error = t1 - gc_s(icol,js)
                  error_sum = error_sum + error*error
                  IF (t1 .eq. zero) t1 = EPSILON(t1)
                  IF (ABS(error) .gt. 1.E-10_dp) THEN
                  WRITE (34, 101) js, gc(icol,js), t1,                 &
                                  gc_s(icol,js), error, error/t1
                  ELSE
                  WRITE (34, 100) js, gc(icol,js), t1,                 &
                                  gc_s(icol,js), error, error/t1
                  END IF

               DO js = 2, ns-1
                  t1 = SUM(lblk_s(icol,:,js)*gc(:,js-1)                &
                         + dblk_s(icol,:,js)*gc(:,js)                  &
                         + ublk_s(icol,:,js)*gc(:,js+1))
                  error = t1 - gc_s(icol,js)
                  error_sum = error_sum + error*error
                  IF (t1 .eq. zero) t1 = EPSILON(t1)
                  IF (ABS(error) .gt. 1.E-10_dp) THEN
                  WRITE (34, 101) js, gc(icol,js), t1,                 &
                         gc_s(icol,js), error, error/t1
                  ELSE
                  WRITE (34, 100) js, gc(icol,js), t1,                 &
                         gc_s(icol,js), error, error/t1
                  END IF
               END DO

                 js = ns        
                 t1 = SUM(lblk_s(icol,:,js)*gc(:,js-1)                 &
                        + dblk_s(icol,:,js)*gc(:,js))
                 error = t1 - gc_s(icol,js)
                 error_sum = error_sum + error*error
                 IF (t1 .eq. zero) t1 = EPSILON(t1)
                 IF (ABS(error) .gt. 1.E-10_dp) THEN
                 WRITE (34, 101) js, gc(icol,js), t1,                  &
                        gc_s(icol,js), error, error/t1
                 ELSE
                 WRITE (34, 100) js, gc(icol,js), t1,                  &
                        gc_s(icol,js), error, error/t1
                 END IF
              END DO
            END DO
         END DO

         DEALLOCATE(dblk_s, lblk_s, ublk_s, gc_s, stat=istat)

         PRINT *,' Hessian Error: ', SQRT(error_sum/b_sum)

      END IF
#endif

 100  FORMAT(i6,1p,5e14.4)
 101  FORMAT(i6,1p,5e14.4,'  *')

      END SUBROUTINE block_precond

 
      SUBROUTINE apply_precond(gc)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(mblk_size,ns), INTENT(inout) :: gc
!      LOGICAL, INTENT(in)                              :: ltype
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp) :: ton, toff
!-----------------------------------------------
      CALL second0(ton)
!      IF (ltype) THEN
!         CALL diag_precond(gc)
!      ELSE
         CALL block_precond(gc)
!      END IF 
      CALL second0(toff)
      time_apply_precon = time_apply_precon+(toff-ton)

      END SUBROUTINE apply_precond


      SUBROUTINE blk3d_slv(ablk, qblk, bp1, source, ipiv, mblk, nblocks)
      USE stel_kinds
      USE stel_constants, ONLY: one
!      USE safe_open_mod
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: bytes_per_rprec2 = 8
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(in) :: ipiv(mblk,nblocks)
      REAL(dp), TARGET, DIMENSION(mblk,mblk,nblocks), INTENT(IN) ::     &
     &                           ablk, qblk, bp1
      REAL(dp), DIMENSION(mblk,nblocks), TARGET, INTENT(INOUT) :: source
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, POINTER  :: ipivot(:)
!      INTEGER :: ibuph, incnow, irecl, incbu, iunit=102, ndisk
      INTEGER :: k, ier
      REAL(dp), POINTER :: amat(:,:), x1(:), source_ptr(:)
!-----------------------------------------------
!  modified (June, 2003, ORNL):         S. P. Hirshman
!-----------------------------------------------------------------------
!
!  this subroutine solves a block-tridiagonal system of equations, using 
!  the ABLK, QBLK factors from blk3d_factor,
!
!-----------------------------------------------------------------------
!  INPUT
!  mblk                : block size
!  nblocks             : number of blocks
!  bp1                 : upper blocks (see equation below)
!  ipiv                : pivot elements for kth block
!  ablk                : a-1 blocks
!  qblk                : q = a-1 * bm1
!  source              : input right side
!
!  OUTPUT
!  source              : Solution x of A x = source
! 
!  LOCAL VARIABLES
!  iunit               : unit number for block-tridiagonal solution disk file.
!
!  the tri-diagonal equation is:
!
!           bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) = source(l)
!
!     GENERAL SOLUTION SCHEME APPLIED TO EACH BLOCK ROW (INDEX L)
!
!     1. Start from row N and solve for x(N) in terms of x(N-1):
!
!        x(N) = -q(N)*x(N-1) + r(N)
!
!        q(N) =  a(N)[-1] * bm1;    r(N) = a(N)[-1] * s(N)
!
!        where a(N)[-1] is the inverse of a(N)
!
!     2. Substitute for lth row to get recursion equation fo q(l) and r(l):
!
!        x(l) = -q(l)*x(l-1) + r(l), in general, where:
!
!        q(l) = (a(l) - bp1(l)*q(l+1))[-1] * bm1(l)
!
!        qblk(l) == (a(l) - bp1(l) * q(l+1))[-1] on return
!
!        r(l) = (a(l) - bp1(l)*q(l+1))[-1] * (s(l) - bp1(l)*r(l+1))
!
!     3. At row l = 1, bm1(1) = 0 and get an equation for x(1) corresponding to q(1) = 0:
!
!        x(1) = r(1)
!
!     4. Finally, can back-solve for x(N) in terms of x(N-1) from eqn.(1) above
!
!
!     NUMERICAL IMPLEMENTATION (USING LAPACK ROUTINES)
!
!     1. CALL dgetrs:   With single right hand side (source) to solve A x = b (b a vector)
!                       Faster than dgesl
!      ndisk = mblk

!  create disk file for doing direct access i/o.

!      incnow = ndisk
!      irecl  = bytes_per_rprec2*incnow
!      incbu = 1 + (ndisk - 1)/incnow
!      ibuph = 0

!      iunit = 10
!      CALL safe_open(iunit, ier, 'NULL', 'scratch', 'unformatted',
!     1     irecl, 'DIRECT')
!      IF (ier .ne. 0) STOP 'Error opening scratch file in blk3d'

!  main loop. load and process (backwards) block-rows nblocks to 1. 
!  note: about equal time is spent in calling dgetrs and in performing
!  the two loop sums: on ibm-pc, 2 s (trs) vs 3 s (sums); on linux (logjam),
!  2.4 s (trs) vs 3 s (sums).

!
!     Back-solve for modified sources first
!
      BLOCKS: DO k = nblocks, 1, -1

         source_ptr => source(:,k)
         amat => ablk(:,:,k)
#if defined(SVD_ON)
         !source_sp = MATMUL(amat, source_sp)
#else
         ipivot => ipiv(:,k);   
         CALL dgetrs('n', mblk, 1, amat, mblk, ipivot, source_ptr, mblk, ier)
         IF (ier .ne. 0) GOTO 305
#endif
         IF (k .eq. 1) EXIT
!
!        NOTE: IN BLK3D_FACTOR, BP1 AND BM1 WERE TRANSPOSED (AND STORED)
!        TO MAKE FIRST INDEX FASTEST VARYING IN THE FOLLOWING MATMUL OPS
!
         amat => bp1(:,:,k-1)
         x1 => source(:,k);  source_ptr => source(:,k-1)
!         source_ptr = source_ptr - MATMUL(x1,amat)  !USE THIS FORM IF TRANSPOSED bp1
!         source_ptr = source_ptr - MATMUL(amat,x1)  !UNTRANSPOSED FORM
         CALL dgemv('T',mblk,mblk,-one,amat,mblk,x1,1,one,source_ptr,1)

      END DO BLOCKS
!
!  forward (back-substitution) solution sweep for block-rows k = 2 to nblocks
!  now, source contains the solution vector
!
      DO k = 2, nblocks
!         CALL rddisk (iunit, ql, ndisk, incnow, ibuph, ier)
!         IF (ier .ne. 0) GOTO 303
!         ibuph = ibuph - incbu

         amat => qblk(:,:,k)
         x1 => source(:,k-1);  source_ptr => source(:,k)
!         source_ptr = source_ptr - MATMUL(x1,amat)  !USE THIS FORM IF TRANSPOSED qblk
!         source_ptr = source_ptr - MATMUL(amat,x1)  !UNTRANSPOSED FORM
         CALL dgemv('T',mblk,mblk,-one,amat,mblk,x1,1,one,source_ptr,1)

      END DO

      GOTO 400

!  error returns. ------------------------------------------------------

  301 CONTINUE
!      WRITE (6, '(a,i8)') ' BLK3D:   error in opening file:  ',
!     1   'RECL = ', irecl
  302 CONTINUE
      WRITE (6, '(a)') ' BLK3D:   error in I/O routine WRDISK'
  303 CONTINUE
      WRITE (6, '(a)') ' BLK3D:   error in I/O routine RDDISK'
      ier = -2
  305 CONTINUE
      WRITE (6, '(2/a,i10,2/)') ' BLK3D-SLV:   error detected:   ier =', ier
      STOP

!  destroy disk file and return. ---------------------------------------

  400 CONTINUE

!      CLOSE (iunit)

      END SUBROUTINE blk3d_slv


      SUBROUTINE blk3d_slv_swp(source, ipiv, mblk, nblocks)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_kinds
      USE stel_constants, ONLY: one
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(in) :: ipiv(mblk,nblocks)
      REAL(rprec), DIMENSION(mblk,nblocks), TARGET, INTENT(inout) :: source
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), POINTER :: amat(:,:), x1(:), source_ptr(:)
      INTEGER, POINTER  :: ipivot(:)
      INTEGER :: k, k1, ier
!-----------------------------------------------
!  modified (June, 2003, ORNL):         S. P. Hirshman
!  modified (June, 2007, ORNL), added lswap2disk logic
!-----------------------------------------------------------------------
!
!  this subroutine solves a block-tridiagonal system of equations, using 
!  the ABLK, QBLK factors from blk3d_factor,
!  See blk3d_slv for details
!
!-----------------------------------------------------------------------
!  INPUT
!  mblk                : block size
!  nblocks             : number of blocks
!  ipiv                : pivot elements for kth block
!  source              : input right side
!
!  OUTPUT
!  source              : Solution x of A x = source
! 
!
!  the tri-diagonal equation is:
!
!           bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) = source(l)
!

!  main loop. load and process (backwards) block-rows nblocks to 1. 
!  note: about equal time is spent in calling dgetrs and in performing
!  the two loop sums: on ibm-pc, 2 s (trs) vs 3 s (sums); on linux (logjam),
!  2.4 s (trs) vs 3 s (sums).

      CALL OpenDAFile(mblk**2, mblk**2, 3, ScratchFile, iunit_dacess, OpenExisting) 

      ALLOCATE (amat(mblk,mblk),  stat=ier)
      IF (ier .ne. 0) STOP 'Allocation error in blk3d_slv_swp!' 
!
!     Back-solve for modified sources first
!
      BLOCKS: DO k = nblocks, 1, -1

         source_ptr => source(:,k)
         ipivot => ipiv(:,k)
         CALL ReadDAItem2(amat, k, bldia)
         CALL dgetrs('n', mblk, 1, amat, mblk, ipivot, source_ptr, mblk, ier)

         IF (ier .ne. 0) GOTO 305
         IF (k .eq. 1) EXIT

!
!        NOTE: IN BLK3D_FACTOR, BP1 AND BM1 WERE TRANSPOSED (AND STORED)
!        TO MAKE FIRST INDEX FASTEST VARYING IN THE FOLLOWING MATMUL OPS
!
         k1 = k-1
         CALL ReadDAItem2(amat, k1, blpls) 
         x1 => source(:,k);  source_ptr => source(:,k1)
!         source_ptr = source_ptr - MATMUL(x1,amat)  !USE THIS FORM IF TRANSPOSED bp1
!         source_ptr = source_ptr - MATMUL(amat,x1)  !UNTRANSPOSED FORM
         CALL dgemv('T',mblk,mblk,-one,amat,mblk,x1,1,one,source_ptr,1)

      END DO BLOCKS
!
!  forward (back-substitution) solution sweep for block-rows k = 2 to nblocks
!  now, source contains the solution vector
!
      DO k = 2, nblocks

         CALL ReadDAItem2(amat, k, blmin)
         x1 => source(:,k-1);  source_ptr => source(:,k)
!         source_ptr = source_ptr - MATMUL(x1,amat)  !USE THIS FORM IF TRANSPOSED qblk
!         source_ptr = source_ptr - MATMUL(amat,x1)  !UNTRANSPOSED FORM
         CALL dgemv('T',mblk,mblk,-one,amat,mblk,x1,1,one,source_ptr,1)

      END DO

      WRITE(100,'(1p,6e14.4)') source(:,ns/2)

      GOTO 400

!  error returns. ------------------------------------------------------

  305 CONTINUE
      WRITE (6, '(2/a,i10,2/)') ' BLK3D_SLV_SWP:   error detected:   ier =', ier
      STOP

  400 CONTINUE

      CALL CloseDAFile

      END SUBROUTINE blk3d_slv_swp

      
      SUBROUTINE dealloc_hessian
      INTEGER :: istat

      IF (ALLOCATED(ublk)) DEALLOCATE(ublk, dblk, lblk, stat=istat)
      IF (ALLOCATED(ublkp))DEALLOCATE(dblkp,lblkp,ublkp,stat=istat)
      IF (ALLOCATED(mupar_norm)) DEALLOCATE(mupar_norm, stat=istat)
      IF (ALLOCATED(ipiv_blk)) DEALLOCATE (ipiv_blk)

      END SUBROUTINE dealloc_hessian

#if defined(MPI_OPT)
      SUBROUTINE blk3d_parallel(a,b,c,ap,bp,cp,mblk,mblkc,nblocks)
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER :: ipart,jpart,i,j,k
      INTEGER :: ia,ja,inc1,ib,jb,inc2,nsm
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk, mblkc
        REAL(dp) :: aij2,aij
      REAL(dp), TARGET, DIMENSION(mblk,mblkc,nblocks), INTENT(inout) :: &
     &                       a, b, c
      REAL(dp), TARGET, DIMENSION(Locq,Locp,nblocks), INTENT(inout) :: &
     &                       ap
      REAL(dp), TARGET, DIMENSION(Locq,Locp,nblocks-1), INTENT(inout) :: &
     &                       bp, cp
      REAL(dp), POINTER :: amat(:,:), bmat(:,:), cmat(:,:)
      REAL(dp), POINTER :: amatp(:,:), bmatp(:,:), cmatp(:,:)

      do k=1,nblocks
        call profstart('change context:123')
        call blacs_barrier(descA(CTXT_),'All')
         amat => a(:,:,k)
         amatp => ap(:,:,k)
          call pdgemr2d(mblk,mblk,amat,1,1,descA_1xp,                            &
     &                       amatp,1,1,descA,  icontxt_global)
        if(k<nblocks)then
         bmat => b(:,:,k+1)
         bmatp => bp(:,:,k)
          call pdgemr2d(mblk,mblk,bmat,1,1,descA_1xp,                            &
     &                       bmatp,1,1,descA,  icontxt_global)
         cmat => c(:,:,k)
         cmatp => cp(:,:,k)
          call pdgemr2d(mblk,mblk,cmat,1,1,descA_1xp,                            &
     &                       cmatp,1,1,descA,  icontxt_global)
       endif
        call blacs_barrier(descA(CTXT_),'All')
        call profend('change context:123')
!        do ja=1,mblk
!          call pdelget( 'A',' ',aij, amatp,ja,ja,descA )
!          call pdelget( 'A',' ',aij2, amat,ja,ja,descA_1xp )
!              write(*,*) 'ja,ja,aij ',ja,ja,aij,aij2
!        enddo
       end do
      END SUBROUTINE blk3d_parallel
#endif

      END MODULE hessian
