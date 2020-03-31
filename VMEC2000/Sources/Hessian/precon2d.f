      MODULE precon2d
      USE stel_kinds, ONLY: dp
      USE vmec_dim
      USE vmec_params
      USE vparams, ONLY: nthreed, one, zero
      USE vmec_input, ONLY: ntor, nzeta, lfreeb, lasym
      USE timer_sub
      USE safe_open_mod
      USE directaccess
      USE parallel_include_module

      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PRIVATE, PARAMETER :: sp = dp
      INTEGER, PRIVATE :: ntyptot, m_2d, n_2d, ntype_2d
      INTEGER, PRIVATE, ALLOCATABLE :: ipiv_blk(:,:)
      INTEGER, PRIVATE :: mblk_size
      INTEGER, PRIVATE :: mystart(3), myend(3) 
      LOGICAL, PRIVATE :: FIRSTPASS=.TRUE.
      REAL(sp),PRIVATE, ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: 
     1    block_diag, block_plus, block_mins,
     2    block_dsave, block_msave, block_psave
      REAL(sp),PRIVATE, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) ::
     1    block_diag_sw, block_plus_sw, block_mins_sw
   
      REAL(dp), PRIVATE, DIMENSION(:,:,:,:), ALLOCATABLE :: gc_save
      REAL(dp) :: ctor_prec2d
      INTEGER :: ictrl_prec2d
      LOGICAL :: lHess_exact = .TRUE.,  !FALSE -> FASTER, LESS ACCURATE VACUUM CALCULATION OF HESSIAN
     1           l_backslv = .FALSE.,
     2           l_comp_prec2D = .TRUE., 
     3           l_edge  = .FALSE.,                !=T IF EDGE PERTURBED
     4           edge_mesh(3)

      LOGICAL, PARAMETER, PRIVATE :: lscreen = .FALSE.
      INTEGER, PARAMETER, PRIVATE :: jstart(3) = (/1,2,3/)

      PRIVATE :: swap_forces, reswap_forces 

!
!     Direct-Access (swap to disk) stuff
!      CHARACTER(LEN=3)   :: FlashDrive ="F:\"
      CHARACTER(LEN=3)   :: FlashDrive =""
      CHARACTER(LEN=128) :: ScratchFile=""
      INTEGER, PARAMETER :: blmin=1, bldia=2, blpls=3
      REAL(sp), ALLOCATABLE, DIMENSION(:,:,:) :: DataItem
      INTEGER            :: iunit_dacess=10
      LOGICAL :: lswap2disk = .FALSE.                                 !Set internally if blocks do not fit in memory
      INTEGER, PARAMETER :: LOWER=3,DIAG=2,UPPER=1

C-----------------------------------------------
!
!     SP:        forces single precision for blocks (smaller size)
!     ICTRL_PREC2D: controls initialization and application of 2d block preconditioner
!                   = 0, no preconditioner applied
!                   = 1, apply preconditioner
!                   = 2, initial call of funct3d to set up residue vector, store saved vectors,
!                        and (for .not.lasym case), call LAMBLKS routine
!                   = 3, radial jog vector is being computed to calculate hessian elements
!     L_BACKSLV: if true, test that Hessian is inverted correctly by back-solving
!     LHESS_EXACT : if true, edge value of ctor (in bcovar) is computed as a constant to make
!                   Hessian symmetric. Also, sets ivacskip=0 in call to vacuum in computation of
!                   Hessian. However, the ivacskip=0 option is (very) slow and found not to be necessary
!                   in practice. Set this true primarily for debugging purposes (check Ap ~ -p in MatVec 
!                   routine, for example, in GMRes module)
!

      CONTAINS

      SUBROUTINE swap_forces(gc, temp, mblk, nblocks)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: mblk, nblocks
      REAL(dp), DIMENSION(nblocks,mblk), INTENT(in)  :: gc
      REAL(dp), DIMENSION(mblk,nblocks), INTENT(out) :: temp
C-----------------------------------------------
!
!     reorders forces (gc) array prior to applying 
!     block-tridiagonal pre-conditioner. on exit, temp is the reordered array
!     flip sign so eigenvalue is negative (corresponding to damping)
!
      temp = -TRANSPOSE(gc)
    
      END SUBROUTINE swap_forces

      SUBROUTINE reswap_forces(temp, gc, mblk, nblocks)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: mblk, nblocks
      REAL(dp), DIMENSION(nblocks,mblk), INTENT(inout)  :: gc
      REAL(dp), DIMENSION(mblk,nblocks), INTENT(in) :: temp
C-----------------------------------------------
!
!     Following application of block pre-conditioner, restores original
!     order of forces (gc) array previously ordered by call to "swap_forces"
!
      gc = TRANSPOSE(temp)

      END SUBROUTINE reswap_forces

      SUBROUTINE block_precond(gc)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor,0:mpol1,ntyptot) :: gc
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mblk, m, n, js, ntype, istat
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: temp
      REAL(dp) :: t1, error
C-----------------------------------------------
!
!     Applies 2D block-preconditioner to forces vector (gc)
!
      IF (ntyptot .le. 0) STOP 'ntyptot must be > 0'

      IF (l_backslv) gc_save = gc

      mblk = ntyptot*mnsize
      ALLOCATE (temp(mblk,ns), stat=ntype)
      IF (ntype .ne. 0) STOP 'Allocation error1 in block_precond'

!     Reorder gc(JS,MN) -> temp(MN,JS)
      CALL swap_forces(gc, temp, mblk, ns)

!     Apply preconditioner to temp, using LU factors stored in block_... matrices
      IF (lswap2disk) THEN
         CALL blk3d_slv_swp(temp, ipiv_blk, mblk, ns)
      ELSE
         CALL blk3d_slv(block_diag, block_mins, block_plus, temp, 
     1                  ipiv_blk, mblk, ns)
      END IF

!     Restores original ordering (after preconditioner applied): temp(MN,JS) -> gc(JS,MN)
      CALL reswap_forces(temp, gc, mblk, ns)      

      IF (l_backslv) THEN
         l_backslv = .false.

         WRITE (6, *) ' Writing block Hessian check to unit 34'
         WRITE (34, *)
         WRITE (34, *) ' BLK3D FACTORIZATION CHECK: Ax = b ?'
        
         DO n = 0, ntor
            WRITE (34, *) ' N = ', n
            DO m = 0, mpol1
               WRITE (34, *) ' M = ', m
               DO ntype = 1, ntyptot
                  WRITE (34, *) ' TYPE = ', ntype
                  WRITE (34, *) 
     1            '   js        Ax             b          Ax - b' //
     2            '        RelErr'
                  js = 1
                  t1 = SUM(block_dsave(n,m,ntype,:,:,:,js)*gc(js,:,:,:)
     1               + block_psave(n,m,ntype,:,:,:,js)*gc(js+1,:,:,:))

                  error = t1 + gc_save(js,n,m,ntype)
                  IF (t1 .eq. zero) t1 = EPSILON(t1)
                  WRITE (34, 100) js, t1, -gc_save(js,n,m,ntype),
     2                            error, error/t1

               DO js = 2, ns-1
                  t1 = SUM(
     1                   block_msave(n,m,ntype,:,:,:,js)*gc(js-1,:,:,:)
     2                 + block_dsave(n,m,ntype,:,:,:,js)*gc(js,:,:,:)
     3                 + block_psave(n,m,ntype,:,:,:,js)*gc(js+1,:,:,:))
                  error = t1 + gc_save(js,n,m,ntype)
                  IF (t1 .eq. zero) t1 = EPSILON(t1)
                  WRITE (34, 100) js, t1, -gc_save(js,n,m,ntype),
     2                            error, error/t1
               END DO

               js = ns        
               t1 = SUM(block_msave(n,m,ntype,:,:,:,js)*gc(js-1,:,:,:)
     1            + block_dsave(n,m,ntype,:,:,:,js)*gc(js,:,:,:))
               error = t1 + gc_save(js,n,m,ntype)
               IF (t1 .eq. zero) t1 = EPSILON(t1)
               WRITE (34, 100) js, t1, -gc_save(js,n,m,ntype), 
     2                         error, error/t1
            END DO
         END DO
         END DO

         IF (.not.l_backslv) DEALLOCATE(block_dsave, block_msave, 
     1       block_psave, gc_save, stat=istat)

 100  FORMAT(i6,1p,4e14.4)
      END IF


      DEALLOCATE (temp, stat=istat)

      END SUBROUTINE block_precond

      SUBROUTINE block_precond_par(gc)
      USE blocktridiagonalsolver, ONLY: SetMatrixRHS
      USE blocktridiagonalsolver, ONLY: BackwardSolve 
      USE blocktridiagonalsolver, ONLY: GetSolutionVector
      USE parallel_include_module
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,0:mpol1,ns,ntyptot) :: gc
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mblk, istat, globrow
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: tmp
      REAL(dp) :: ton, toff
      REAL(dp), DIMENSION (:,:), ALLOCATABLE :: solvec
C-----------------------------------------------
      IF (.NOT.lactive) THEN
         RETURN
      END IF
!
!     Applies 2D block-preconditioner to forces vector (gc)
!
      IF (ntyptot .LE. 0) THEN
         STOP 'ntyptot must be > 0'
      END IF

      mblk = ntyptot*mnsize

!     Apply preconditioner to temp, using LU factors stored in block_... matrices

      ALLOCATE (tmp(mblk,ns), stat=istat)
      CALL tolastns(gc, tmp)
      tmp(:,tlglob:trglob) = -tmp(:,tlglob:trglob)
      
      DO globrow=tlglob, trglob
         CALL SetMatrixRHS(globrow,tmp(:,globrow))
      END DO
      DEALLOCATE (tmp, stat=istat)

      CALL second0(ton)
      CALL BackwardSolve
      CALL second0(toff)
      bcyclic_backwardsolve_time=bcyclic_backwardsolve_time+(toff-ton)
       
      ALLOCATE (solvec(mblk,ns), stat=istat)
      IF (istat .NE. 0) THEN
         STOP 'Allocation error in block_precond before gather'
      END IF

      DO globrow=tlglob, trglob
         CALL GetSolutionVector (globrow,solvec(:,globrow))
      END DO

      CALL tolastntype(solvec, gc)

      CALL Gather4XArray(gc)
      
      DEALLOCATE (solvec)
   
      END SUBROUTINE block_precond_par

      SUBROUTINE compute_blocks_par (xc, xcdot, gc)
      USE blocktridiagonalsolver, ONLY: ForwardSolve
      USE parallel_include_module
      USE vmec_main, ONLY: iter2
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp),DIMENSION(0:ntor,0:mpol1,ns,3*ntmax) :: xc, gc, xcdot
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(dp), PARAMETER :: p5 = 0.5_dp
      INTEGER :: m, n, i, ntype, istat, mblk, ibsize, iunit
      REAL(dp) :: time_on, time_off, bsize, tprec2don, tprec2doff
      REAL(dp) :: ton, toff
      CHARACTER(LEN=100):: label
      LOGICAL, PARAMETER :: lscreen = .false.
      INTEGER :: j, k, l
C-----------------------------------------------
!
!     COMPUTES THE JACOBIAN BLOCKS block_mins, block_diag, block_plus
!     USING EITHER A SLOW - BUT RELIABLE - "JOG" TECHNIQUE, OR
!     USING PARTIAL ANALYTIC FORMULAE.
!
!     THE SUBROUTINE lam_blks IS CALLED FROM BCOVAR TO COMPUTE 
!     THE ANALYTIC BLOCK ELEMENTS
!
      CALL second0(tprec2don)
      
      IF (l_backslv .and. sp.ne.dp) THEN
         STOP 'Should set sp = dp!'
      END IF

      ntyptot = SIZE(gc,4)
      IF (ntyptot .NE. 3*ntmax) THEN
         STOP ' NTYPTOT != 3*ntmax'
      END IF
      mblk = ntyptot*mnsize

      bsize = REAL(mblk*mblk, dp)*3*KIND(block_diag)
      IF (bsize .gt. HUGE(mblk)) THEN
         WRITE (6, *) ' bsize: ', bsize, ' exceeds HUGE(int): ',
     &                HUGE(mblk)
!        WRITE (6, *) ' Blocks will be written to disk.'
!        lswap2disk = .TRUE.
      ELSE
         lswap2disk = .FALSE.
      END IF

      bsize = bsize*ns
      IF (bsize .lt. 1.E6_dp) THEN
         ibsize = bsize/1.E1_dp
         label = " Kb"
      ELSE IF (bsize .lt. 1.E9_dp) THEN
         ibsize = bsize/1.E4_dp
         label = " Mb"
      ELSE
         ibsize = bsize/1.E7_dp
         label = " Gb"
      END IF

      DO i = 1,2
         IF (i .eq. 1) THEN
            iunit = 6
         END IF
         IF (i .eq. 2) THEN
            iunit = nthreed
         END IF
         IF (grank.EQ.0) THEN
            WRITE (iunit, '(/,2x,a,i5,a,/,2x,a,i5,a)')
     &         'Initializing 2D block preconditioner at ', iter2,
     &         ' iterations',
     &         'Estimated time to compute Hessian = ',
     &         3*ntyptot*mnsize,' VMEC time steps'
            WRITE (iunit, '(2x,a,i4,a,f12.2,a)') 'Block dim: ', mblk,
     &         '^2  Preconditioner size: ', REAL(ibsize)/100,
     &         TRIM(label)
         END IF
      END DO
!
!     COMPUTE AND STORE BLOCKS (MN X MN) FOR PRECONDITIONER
!
      CALL second0(time_on)

      ALLOCATE (gc_save(0:ntor,0:mpol1,ns,ntyptot), stat=istat)
      IF (istat .NE. 0) THEN
         STOP 'Allocation error: gc_save in compute_blocks'
      END IF

      IF (ALLOCATED(block_diag)) THEN
         DEALLOCATE (block_diag, block_plus, block_mins, stat=istat)
         IF (istat .ne. 0) THEN
            STOP 'Deallocation error in compute blocks'
         END IF
      END IF

!
!     GENERAL (SLOWER BY 2/3 THAN SYMMETRIC VERSION) METHOD: ASSUMES NO SYMMETRIES OF R, Z COEFFICIENTS
!
      CALL sweep3_blocks_par (xc, xcdot, gc)
      IF (lactive) THEN
         CALL compute_col_scaling_par
      END IF

      ictrl_prec2d = 1                 !Signals funct3d (residue) to use block preconditioner

      CALL second0(time_off)
      IF (grank .EQ. 0) THEN
         WRITE (6,1000) time_off - time_on
         WRITE (nthreed,1000) time_off - time_on
      END IF

!
!     FACTORIZE HESSIAN 
!
      CALL second0(time_on)
      IF (ALLOCATED(ipiv_blk)) THEN
         DEALLOCATE(ipiv_blk, stat=ntype)
      END IF
      ALLOCATE (ipiv_blk(mblk,ns), stat=ntype)     
      IF (ntype .ne. 0) THEN
         STOP 'Allocation error2 in block_precond'
      END IF

      CALL second0(ton)
      IF (lactive) THEN
         CALL ForwardSolve
      END IF

      CALL second0(time_off)
      toff = time_off
      bcyclic_forwardsolve_time = bcyclic_forwardsolve_time
     &                          + (toff - ton)

      IF (grank.EQ.0) THEN
         WRITE(6,1001) time_off - time_on
         WRITE(nthreed,1001) time_off - time_on
      END IF

      IF (.NOT.l_backslv) THEN
         DEALLOCATE (gc_save)
      END IF

      CALL second0(tprec2doff)

      timer(tprec2d) = timer(tprec2d) + (tprec2doff - tprec2don)
      compute_blocks_time = compute_blocks_time
     &                    + (tprec2doff - tprec2don)

1000  FORMAT(1x,' Time to compute blocks: ',f10.2,' s')
1001  FORMAT(1x,' Time to factor blocks:  ',f10.2,' s')

      END SUBROUTINE compute_blocks_par

      SUBROUTINE sweep3_blocks_par(xc, xcdot, gc)
      USE vmec_main, ONLY: ncurr, r01, z01, lthreed, chips, delt0r
      USE blocktridiagonalsolver, ONLY: Initialize, SetBlockRowCol
      USE blocktridiagonalsolver, ONLY: WriteBlocks
      USE parallel_vmec_module, ONLY: MPI_STAT
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,0:mpol1,ns,ntyptot) :: xc, xcdot, gc
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, js1, istat, mesh, lamtype, rztype, icol
      INTEGER :: nsmin, nsmax
      INTEGER :: lastrank, left, right
      REAL(dp) :: eps, hj, hj_scale
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: diag_val
      REAL(dp) :: ton, toff
C-----------------------------------------------
!
!     COMPUTE FORCE "RESPONSE" TO PERTURBATION AT EVERY 3rd RADIAL POINT
!     FOR EACH MESH STARTING AT js=1,2,3, RESPECTIVELY
!
      CALL second0(ton)
      ALLOCATE(diag_val(ns), stat=istat)
      diag_val=zero
      IF (grank.EQ.0) THEN
        WRITE (6, *)
     &   " Using non-symmetric sweep to compute Hessian elements"
      END IF

      eps = SQRT(EPSILON(eps))
      eps = eps/10
      rztype = 2*ntmax
      lamtype = rztype + 1

      n_2d = 0
      m_2d = 0

      ALLOCATE(DataItem(0:ntor,0:mpol1,1:3*ntmax), stat=istat)
      IF (istat .ne. 0) THEN
         STOP 'Allocation error in sweep3_blocks'
      END IF

!
!     CALL FUNCT3D FIRST TIME TO STORE INITIAL UN-PRECONDITIONED FORCES
!     THIS WILL CALL LAMBLKS (SO FAR, ONLY IMPLEMENTED FOR lasym = false)
!
      ictrl_prec2d = 2                 !Signals funct3d that preconditioner is being initialized
      CALL funct3d_par(lscreen, istat)
      IF (istat .NE. 0) THEN
         PRINT *,' ier_flag = ', istat,
     1           ' in SWEEP3_BLOCKS_PAR call to funct3d_par'
         STOP
      ENDIF

      nsmin = t1lglob
      nsmax = t1rglob
      xcdot(:,:,nsmin:nsmax,:) = 0
!      PRINT *,'rank: ', rank,' vrank: ', vrank,' nsmin: ',nsmin,
!     1        ' nsmax: ', nsmax

      IF (FIRSTPASS) THEN
         edge_mesh = .FALSE.
         FIRSTPASS = .FALSE.
         mblk_size = (ntor + 1)*(mpol1 + 1)*3*ntmax
         IF (mblk_size .NE. ntmaxblocksize) THEN
            STOP 'wrong mblk_size in precon2d!'
         END IF
         CALL Initialize(.FALSE.,ns,mblk_size)
         myend = nsmax
!Align starting pt in (nsmin,nsmax)
         DO mesh = 1, 3
            icol = MOD(jstart(mesh) - nsmin, 3)
            IF (icol .LT. 0) THEN
               icol = icol + 3
            END IF
            mystart(mesh) = nsmin + icol
            IF (MOD(jstart(mesh) - ns, 3) .EQ. 0) THEN     ! .AND. nsmax.EQ.ns)
               edge_mesh(mesh) = .TRUE.
            END IF
         END DO
      END IF

!     STORE chips in xc 
#if defined(CHI_FORCE)
      IF (ncurr .EQ. 1) THEN
         xc(0,0,nsmin:nsmax,lamtype) = chips(nsmin:nsmax)
      END IF
#endif
      left = rank - 1
      IF (rank .EQ. 0) THEN
         left = MPI_PROC_NULL
      END IF
      right = rank + 1
      IF (rank .EQ. nranks - 1) THEN
         right = MPI_PROC_NULL
      END IF

      CALL PadSides(xc)
      CALL PadSides(gc)
      CALL PadSides(xcdot)

      CALL restart_iter(delt0r)
      gc_save(:,:,nsmin:nsmax,:) = gc(:,:,nsmin:nsmax,:) 
      ictrl_prec2d = 3                 !Signals funct3d that preconditioner is being computed

      CALL MPI_COMM_SIZE(NS_COMM,lastrank,MPI_ERR)
      lastrank = lastrank - 1
!
!     FIRST DO R00 JOG TO LOAD DIAG_VAL ARRAY (DO NOT RELY ON IT BEING THE FIRST JOG)
!
      m_2d=0
      n_2d=0
#ifdef _HBANGLE
      ntype_2d = zsc + ntmax
#else
      ntype_2d = rcc
#endif
!     APPLY JOG
      hj = eps * MAX(ABS(r01(ns)), ABS(z01(ns)))
      IF (nranks .GT. 1) THEN
         CALL MPI_BCAST(hj,1,MPI_REAL8,lastrank,NS_COMM,MPI_ERR)
         CALL MPI_BCAST(edge_mesh,3,MPI_REAL8,lastrank,NS_COMM,MPI_ERR)
      END IF
      DO js = mystart(1), myend(1), 3
         xcdot(n_2d,m_2d,js,ntype_2d) = hj
      END DO

      istat = 0
      CALL funct3d_par (lscreen, istat)

      LACTIVE0: IF (lactive) THEN
         IF (nranks .GT. 1) THEN
            CALL Gather4XArray(gc)
         END IF
         IF (istat .NE. 0) THEN
            STOP 'Error computing Hessian jog!'
         END IF
!     CLEAR JOG AND STORE BLOCKS FOR THIS JOG
         xcdot(:,:,nsmin:nsmax,:) = 0
         DO js = mystart(1), myend(1), 3
            DataItem = (gc(:,:,js,:) - gc_save(:,:,js,:))/hj
            diag_val(js) = DataItem(0,0,ntype_2d)
         END DO

         IF (nranks .GT. 1) THEN
            icol = 0
            IF (trglob_arr(1) .LT. 4) THEN
               icol = 1
            END IF
            CALL MPI_BCAST(diag_val(4), 1, MPI_REAL8, icol, NS_COMM,
     &                     MPI_ERR)
            icol = nranks - 1
            IF (tlglob_arr(nranks) .GT. ns - 3) THEN
               icol = nranks-2
            END IF
            CALL MPI_BCAST(diag_val(ns - 3), 1, MPI_REAL8, icol,
     &                     NS_COMM, MPI_ERR)
         END IF
         IF (diag_val(1) .EQ. zero) THEN
            diag_val(1)  = diag_val(4)
         END IF
         IF (diag_val(ns) .EQ. zero) THEN
            diag_val(ns) = diag_val(ns - 3)
         END IF

         IF (nranks .GT. 1) THEN
            CALL MPI_Sendrecv(diag_val(trglob), 1, MPI_REAL8, right, 1,
     &                        diag_val(t1lglob), 1, MPI_REAL8, left, 1,
     &                        NS_COMM, MPI_STAT, MPI_ERR)
         END IF
         DO js = mystart(2), myend(2), 3
            diag_val(js) = diag_val(js - 1)
         END DO

         hj_scale = MAX(ABS(r01(ns)), ABS(z01(ns)))

         IF (nranks .GT. 1) THEN
            CALL MPI_Sendrecv(diag_val(trglob), 1, MPI_REAL8, right, 1,
     &                        diag_val(t1lglob), 1, MPI_REAL8, left, 1,
     &                        NS_COMM, MPI_STAT, MPI_ERR)
            CALL MPI_BCAST(hj_scale, 1, MPI_REAL8, lastrank, NS_COMM,
     &                   MPI_ERR)
         END IF

         DO js = mystart(3), myend(3), 3
            diag_val(js) = diag_val(js - 1)
         END DO

         IF (ANY(diag_val(tlglob:trglob) .EQ. zero)) THEN
            PRINT *, 'For rank: ', rank, ' some diag_val == 0'
            STOP
         END IF
      END IF LACTIVE0
!
!     PERFORM "JOGS" FOR EACH VARIABLE AT EVERY 3rd RADIAL POINT ACROSS MESH
!     FOR ntyp = (Rcc, Rss, Rsc, Rcs, Zsc, Zcs, Zcc, Zss)
!     AND EVERY n2d (toroidal mode index) and EVERY m2d (poloidal mode index)

      icol=0

      NTYPE2D: DO ntype_2d = 1, ntyptot
         hj = eps
         IF (ntype_2d .LT. lamtype) THEN
            hj = hj*hj_scale
         END IF

         M2D: DO m_2d = 0, mpol1
            
            N2D: DO n_2d = 0, ntor

               icol = icol + 1
            
               MESH_3PT: DO mesh = 1,3
 
!              APPLY JOG TO ACTIVE PROCESSORS
                  IF (lactive) THEN
                     DO js = mystart(mesh), myend(mesh), 3
                        xcdot(n_2d,m_2d,js,ntype_2d) = hj
                     END DO
                     IF (m_2d.GT.0 .AND. mystart(mesh).EQ.1) THEN
                        xcdot(n_2d,m_2d,1,ntype_2d) = 0
                     END IF
                  END IF

                  l_edge = edge_mesh(mesh)
                  CALL funct3d_par (lscreen, istat)
                  IF (istat .NE. 0) STOP 'Error computing Hessian jog!'
              
!
!              COMPUTE PRECONDITIONER (HESSIAN) ELEMENTS. LINEARIZED EQUATIONS
!              OF FORM (FIXED mn FOR SIMPLICITY):
!
!              F(j-1) = a(j-1)x(j-2) + d(j-1)x(j-1) + b(j-1)x(j)
!              F(j)   =                a(j)x(j-1)   + d(j)  x(j) + b(j)  x(j+1)
!              F(j+1) =                               a(j+1)x(j) + d(j+1)x(j+1) + b(j+1)x(j+2)
!
!              HESSIAN IS H(k,j) == dF(k)/dx(j); aj == block_mins; dj == block_diag; bj = block_plus
!
!              THUS, A PERTURBATION (in xc) AT POSITION js PRODUCES THE FOLLOWING RESULTS:
!
!                     d(js)   = dF(js  )/hj(js)
!                     b(js-1) = dF(js-1)/hj(js)
!                     a(js+1) = dF(js+1)/hj(js)
!
!
                  LACTIVE1: IF (lactive) THEN
                     SKIP3_MESH: DO js = mystart(mesh), myend(mesh), 3

!                 CLEAR JOG AND STORE BLOCKS FOR THIS JOG
                        xcdot(n_2d,m_2d,js,ntype_2d) = 0

                  !block_mins(js+1) 
                       js1 = js+1
                        IF (tlglob.LE.js1 .AND. js1.LE.trglob) THEN
                           DataItem =
     &                        (gc(:,:,js1,:)-gc_save(:,:,js1,:))/hj
                           CALL SetBlockRowCol(js1,icol,DataItem,LOWER)
                        END IF

                  !block_diag(js)
                        IF (tlglob.LE.js .AND. js.LE.trglob) THEN
                           DataItem = (gc(:,:,js,:)
     &                              - gc_save(:,:,js,:))/hj

                           IF (rank .EQ. lastrank .AND.
     &                         js   .EQ. ns       .AND.
     &                         .NOT.lfreeb        .AND.
     &                         ANY(DataItem(:,:,
     &                                      1:rztype) .NE. zero)) THEN
                              STOP 'DIAGONAL BLOCK AT EDGE != 0'
                           END IF

!Levenberg-like offset - do NOT apply here if applied in colscaling routine
                           IF (ntype_2d .GE. lamtype) THEN
                              DataItem(n_2d,m_2d,ntype_2d) =
     &                           1.0001_dp*DataItem(n_2d,m_2d,ntype_2d)
                           END IF

                           IF (DataItem(n_2d,m_2d,
     &                                  ntype_2d) .EQ. zero) THEN
                              DataItem(n_2d,m_2d,ntype_2d) =
     &                           diag_val(js)
                           END IF

                           CALL SetBlockRowCol(js,icol,DataItem,DIAG)
                        END IF

                 !block_plus(js-1)
                        js1 = js - 1
                        IF (tlglob .LE. js1 .AND. js1 .LE. trglob) THEN
                           DataItem = (gc(:,:,js1,:) -
     &                                 gc_save(:,:,js1,:))/hj
                           CALL SetBlockRowCol(js1,icol,DataItem,UPPER)
                        END IF
                     END DO SKIP3_MESH
                  END IF LACTIVE1

               END DO MESH_3PT
            END DO N2D
         END DO M2D
      END DO NTYPE2D

      l_edge = .FALSE.

      DEALLOCATE(DataItem, diag_val)
      CALL second0(toff)
      fill_blocks_time=fill_blocks_time + (toff - ton)
 
      END SUBROUTINE sweep3_blocks_par

      SUBROUTINE compute_blocks(xc, xcdot, gc)
      USE vmec_main, ONLY: iter2
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp),DIMENSION(ns,0:ntor,0:mpol1,3*ntmax) :: xc, gc, xcdot
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(dp), PARAMETER :: p5 = 0.5_dp
      INTEGER :: m, n, i, ntype, istat,  
     1           mblk, ibsize, iunit
      REAL(dp) :: time_on, time_off, bsize, tprec2don, tprec2doff
      CHARACTER(LEN=100):: label
      LOGICAL, PARAMETER :: lscreen = .false.
C-----------------------------------------------
!
!     COMPUTES THE JACOBIAN BLOCKS block_mins, block_diag, block_plus
!     USING EITHER A SLOW - BUT RELIABLE - "JOG" TECHNIQUE, OR
!     USING PARTIAL ANALYTIC FORMULAE.
!
!     THE SUBROUTINE lam_blks IS CALLED FROM BCOVAR TO COMPUTE 
!     THE ANALYTIC BLOCK ELEMENTS
!
      CALL second0(tprec2don)
      IF (l_backslv .and. sp.ne.dp) THEN
         STOP 'Should set sp = dp!'
      END IF

      ntyptot = SIZE(gc,4)
      mblk = ntyptot*mnsize

      bsize = REAL(mblk*mblk, dp)*3*ns*KIND(block_diag)
!      IF (bsize .gt. HUGE(mblk)) THEN
!        WRITE (6, *) ' bsize: ', INT(bsize,SELECTED_INT_KIND(12)), 
!     1               ' exceeds HUGE(int): ', HUGE(mblk)
!        WRITE (6, *) ' Blocks will be written to disk.'
!        lswap2disk = .TRUE.
!      ELSE
      lswap2disk = .FALSE.
!      END IF

      IF (bsize .lt. 1.E6_dp) THEN
         ibsize = bsize/1.E1_dp
         label = " Kb"
      ELSE IF (bsize .lt. 1.E9_dp) THEN
         ibsize = bsize/1.E4_dp
         label = " Mb"
      ELSE
         ibsize = bsize/1.E7_dp
         label = " Gb"
      END IF

      WRITE (6, 1000) iter2, 3*ntyptot*mnsize, mblk, REAL(ibsize)/100,
     &                TRIM(label)
      WRITE (nthreed, 1000) iter2, 3*ntyptot*mnsize, mblk,
     &                      REAL(ibsize)/100, TRIM(label)
!
!     COMPUTE AND STORE BLOCKS (MN X MN) FOR PRECONDITIONER
!
      CALL second0(time_on)

      ALLOCATE (gc_save(ns,0:ntor,0:mpol1,ntyptot), stat=istat)
      IF (istat .ne. 0) THEN
         STOP 'Allocation error: gc_save in compute_blocks'
      END IF

      IF (ALLOCATED(block_diag)) THEN
          DEALLOCATE (block_diag, block_plus, block_mins, stat=istat)
          IF (istat .ne. 0) THEN
             STOP 'Deallocation error in compute blocks'
          END IF
      ELSE IF (ALLOCATED(block_diag_sw)) THEN
         DEALLOCATE (block_diag_sw, block_plus_sw, block_mins_sw, 
     &                stat=istat)
         IF (istat .ne. 0) THEN
            STOP 'Deallocation error in compute blocks'
         END IF
      END IF

      ALLOCATE (block_diag(0:ntor,0:mpol1,ntyptot,
     &                     0:ntor,0:mpol1,ntyptot,ns),
     &          block_plus(0:ntor,0:mpol1,ntyptot,
     &                     0:ntor,0:mpol1,ntyptot,ns),
     &          block_mins(0:ntor,0:mpol1,ntyptot,
     &                     0:ntor,0:mpol1,ntyptot,ns),
     &          stat=istat)

      lswap2disk = (istat .NE. 0)

!FOR DEBUGGING, SET THIS TO TRUE
!      lswap2disk = .TRUE.
      
      IF (lswap2disk) THEN
         WRITE (6,'(a,i4,a)') '  Allocation error(1) = ', istat,
     &               ': Not enough memory in compute_blocks'
         WRITE (6,*) ' Writing blocks to disk file'
         ALLOCATE (block_diag_sw(0:ntor,0:mpol1,ntyptot,
     &                           0:ntor,0:mpol1,ntyptot),
     &             block_plus_sw(0:ntor,0:mpol1,ntyptot,
     &                           0:ntor,0:mpol1,ntyptot),
     &             block_mins_sw(0:ntor,0:mpol1,ntyptot,
     &                           0:ntor,0:mpol1,ntyptot),
     &             stat=istat)

         IF (istat .ne. 0) THEN
            WRITE (6,'(a,i4)') ' Allocation error(2) = ', istat
            STOP
         END IF

!        Open DIRECT ACCESS file for writing blocks to disk
!        FIRST, we need to compute one row (in m,n,ntype-space) at a
!        time (not the full block). So we use a record size = mblk
!        with a block size = mblk**2. We will then close this and re-open
!        it with a record size = mblk**2 do deal with full blocks
         ScratchFile = "PRCND2A.bin"
         IF (FlashDrive .ne. "") THEN
            ScratchFile = FlashDrive // ScratchFile
         END IF
         CALL OpenDAFile(mblk, mblk**2, 3, ScratchFile, iunit_dacess, 0) 
         ScratchFile = "PRCND2B.bin"
         IF (FlashDrive .ne. "") THEN
            ScratchFile = FlashDrive // ScratchFile
         END IF
         block_plus_sw = 0
         block_mins_sw = 0
         block_diag_sw = 0
      ELSE
         block_plus = 0
         block_mins = 0
         block_diag = 0
      END IF

!
!     GENERAL (SLOWER BY 2/3 THAN SYMMETRIC VERSION) METHOD: ASSUMES NO SYMMETRIES OF R, Z COEFFICIENTS
!
      CALL sweep3_blocks (xc, xcdot, gc)
      CALL compute_col_scaling
      ictrl_prec2d = 1                 !Signals funct3d (residue) to use block preconditioner

!SPH021014: compute eigenvalues (for small enough matrices)
!      CALL get_eigenvalues(mblk, ns, block_mins, block_diag, block_plus)

      CALL second0(time_off)
      WRITE (6,1001) time_off - time_on
      WRITE (6,1001) nthreed

!     SAVE ORIGINAL (UNFACTORED) BLOCKS FOR CHECKING INVERSE
!     IN L_BACKSLV=TRUE LOOP IN BLOCK_PRECOND
      IF (l_backslv) THEN
         ALLOCATE (block_dsave(0:ntor,0:mpol1,ntyptot,
     &                         0:ntor,0:mpol1,ntyptot,ns),
     &             block_msave(0:ntor,0:mpol1,ntyptot,
     &                         0:ntor,0:mpol1,ntyptot,ns),
     &             block_psave(0:ntor,0:mpol1,ntyptot,
     &                         0:ntor,0:mpol1,ntyptot,ns),
     &             stat = istat)
         IF (istat .ne. 0) THEN
            WRITE (6,*) 'Allocation error in l_backslv block: stat = ',
     &                istat
            l_backslv = .false.
         ELSE
            block_dsave = block_diag;  block_msave = block_mins
            block_psave = block_plus
         END IF
      END IF

!
!     FACTORIZE HESSIAN 
!
      CALL second0(time_on)
      IF (ALLOCATED(ipiv_blk)) THEN
         DEALLOCATE(ipiv_blk, stat=ntype)
      END IF
      ALLOCATE (ipiv_blk(mblk,ns), stat=ntype)     
      IF (ntype .ne. 0) STOP 'Allocation error2 in block_precond'

      IF (lswap2disk) THEN
         CALL blk3d_factor_swp(ipiv_blk, mblk, ns)
      ELSE
         CALL blk3d_factor(block_diag, block_mins, block_plus,  
     1                     ipiv_blk, mblk, ns)
      END IF

      CALL second0(time_off)
      WRITE (6,1002) time_off - time_on
      WRITE (nthreed,1002) time_off - time_on

      IF (.NOT.l_backslv) DEALLOCATE (gc_save)

      CALL second0(tprec2doff)

      timer(tprec2d) = timer(tprec2d) + (tprec2doff - tprec2don)

1000  FORMAT(/,2x,'Initializing 2D block preconditioner at ',i5,
     &       ' iterations',
     &       /,2x,'Estimated time to compute Hessian = ',i5,
     &       ' VMEC time steps',
     &       /,2x,'Block dim: ',i4,'^2  Preconditioner size: ',f12.2,a)
1001  FORMAT(1x,' Time to compute blocks: ',f10.2,' s')
1002  FORMAT(1x,' Time to factor blocks:  ',f10.2,' s')

      END SUBROUTINE compute_blocks

      SUBROUTINE sweep3_blocks(xc, xcdot, gc )
      USE vmec_main, ONLY: ncurr, r01, z01, lthreed, chips, delt0r
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor,0:mpol1,ntyptot) :: xc, xcdot, 
     &                                                  gc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp) :: eps, hj, diag_val(ns)
      INTEGER :: js, istat, mesh, lamtype, rztype

!      INTEGER :: m1, n1, nt1

!-----------------------------------------------
!
!     COMPUTE EVERY 3rd RADIAL POINT FOR EACH MESH STARTING AT js=1,2,3, RESPECTIVELY
!
      WRITE (6, *)
     &   " Using non-symmetric sweep to compute Hessian elements"

      eps = SQRT(EPSILON(eps))
      eps = eps/10
      rztype = 2*ntmax
      lamtype = rztype+1

      xcdot = 0
      n_2d = 0;  m_2d = 0

      ALLOCATE(DataItem(0:ntor,0:mpol1,3*ntmax), stat=istat)
      IF (istat .NE. 0) THEN
         STOP 'Allocation error in sweep3_blocks'
      END IF

      diag_val = 1

!
!     CALL FUNCT3D FIRST TIME TO STORE INITIAL UN-PRECONDITIONED FORCES
!
      ictrl_prec2d = 2                 !Signals funct3d to initialize preconditioner (save state)
      CALL funct3d (lscreen, istat)
      IF (istat .ne. 0) THEN
         PRINT *,' ier_flag = ', istat,
     &           ' in SWEEP3_BLOCKS call to funct3d'
         STOP
      ENDIF

#if defined(CHI_FORCE)
!     STORE chips in xc
      IF (ncurr .EQ. 1) THEN
         xc(1:ns,0,0,lamtype) = chips(1:ns)
      END IF
#endif	
      CALL restart_iter(delt0r)
      ictrl_prec2d = 3                 !Signals funct3d to compute preconditioner elements
      gc_save = gc

!
!     PERFORM R(m=0,n=0) JOG TO LOAD DIAG_VAL ARRAY 
!     (DO NOT RELY ON IT BEING THE FIRST JOG IN LOOP)
!
      m_2d = 0
      n_2d = 0
#ifdef _HBANGLE
      ntype_2d = zsc + ntmax
#else
      ntype_2d = rcc
#endif
      edge_mesh = .FALSE.
      DO mesh = 1, 3
         DO js = jstart(mesh), ns, 3
            IF (js .EQ. ns) THEN
               edge_mesh(mesh) = .TRUE.
            END IF
         END DO
      END DO

!     APPLY JOG
      hj = eps * MAX(ABS(r01(ns)), ABS(z01(ns)))
      DO js = jstart(1), ns, 3
         xcdot(js,n_2d,m_2d,ntype_2d) = hj
      END DO

      CALL funct3d (lscreen, istat)
      IF (istat .NE. 0) THEN
         STOP 'Error computing Hessian jog!'
      END IF
!     CLEAR JOG AND STORE BLOCKS FOR THIS JOG
      xcdot = 0
      DO js = jstart(1), ns, 3
         DataItem = (gc(js,:,:,:) - gc_save(js,:,:,:))/hj
         diag_val(js) = DataItem(0,0,ntype_2d)
      END DO
      IF (diag_val(1) .EQ. zero) THEN
         diag_val(1) = diag_val(4)
      END IF
      IF (diag_val(ns).EQ. zero) THEN
         diag_val(ns) = diag_val(ns-3)
      END IF

      DO js = jstart(2), ns, 3
         diag_val(js) = diag_val(js-1)
      END DO
      DO js = jstart(3), ns, 3
         diag_val(js) = diag_val(js-1)
      END DO

      IF (ANY(diag_val .EQ. zero)) THEN
         STOP 'diag_val == 0'
      END IF
!
!     PERFORM "JOGS" FOR EACH VARIABLE AT EVERY 3rd RADIAL POINT ACROSS MESH
!     FOR ntype_2d = (Rcc, Rss, Rsc, Rcs, Zsc, Zcs, Zcc, Zss)
!     AND EVERY n_2d (toroidal mode index) and EVERY m_2d (poloidal mode index)
!
      NTYPE2D: DO ntype_2d = 1, ntyptot
         IF (ntype_2d .LT. lamtype) THEN
            hj = eps*MAX(ABS(r01(ns)), ABS(z01(ns)))
         ELSE
            hj = eps
         END IF

         M2D: DO m_2d = 0, mpol1
            N2D: DO n_2d = 0, ntor
#ifdef _HBANGLE
               IF (ntype_2d.GT.ntmax .AND. ntype_2d.LE.2*ntmax) THEN
                  IF (m_2d .NE. 0) THEN
                     block_diag(n_2d,m_2d,ntype_2d,
     &                          n_2d,m_2d,ntype_2d,:) = diag_val
                     CYCLE
                  END IF
               END IF
#endif
               MESH_3PT: DO mesh = 1,3
!                 APPLY JOG
                  DO js = jstart(mesh), ns, 3
                     xcdot(js,n_2d,m_2d,ntype_2d) = hj
                  END DO

                  l_edge = edge_mesh(mesh)
                  CALL funct3d(lscreen, istat)
                  IF (istat .NE. 0) THEN
                     STOP 'Error computing Hessian jog!'
                  END IF

!               IF (.NOT.lfreeb)   !NOT NEEDED, gcr, gcz -> 0 in RESIDUE for lfreeb=F 
!     1            gc(ns,:,:,1:rztype) = gc_save(ns,:,:,1:rztype)
!
!              COMPUTE PRECONDITIONER (HESSIAN) ELEMENTS. LINEARIZED EQUATIONS
!              OF FORM (FIXED mn FOR SIMPLICITY):
!
!              F(j-1) = a(j-1)x(j-2) + d(j-1)x(j-1) + b(j-1)x(j)
!              F(j)   =                a(j)x(j-1)   + d(j)  x(j)  + b(j)  x(j+1)
!              F(j+1) =                               a(j+1)x(j)  + d(j+1)x(j+1) + b(j+1)x(j+2)
!
!              HESSIAN IS H(k,j) == dF(k)/dx(j); aj == block_mins; dj == block_diag; bj = block_plus
!
!              THUS, A PERTURBATION (in xc) AT POSITION js PRODUCES THE FOLLOWING RESULTS:
!
!                     d(js)   = dF(js  )/hj(js)
!                     b(js-1) = dF(js-1)/hj(js)
!                     a(js+1) = dF(js+1)/hj(js)
!
!              CLEAR JOG
                  xcdot = 0
!
!              STORE BLOCK ELEMENTS FOR THIS JOG.
!              FOR OFF-DIAGONAL ELEMENTS, NEED TO ADJUST js INDICES +/- 1
!
                  SKIP3_MESH: DO js = jstart(mesh), ns, 3

                  !block_mins(js+1) == a
                     IF (js .lt. ns) THEN
                        DataItem = (gc(js+1,:,:,:) -
     &                              gc_save(js+1,:,:,:))/hj
                     END IF

                     IF (lswap2disk) THEN
                        CALL WriteDAItem_SEQ(DataItem)
                     ELSE IF (js .lt. ns) THEN
                        block_mins(:,:,:,n_2d,m_2d,ntype_2d,js+1) =
     &                     DataItem
                     END IF

                  !block_diag(js) == d
                     DataItem = (gc(js,:,:,:) - gc_save(js,:,:,:))/hj

                     IF (DataItem(n_2d,m_2d,ntype_2d) .EQ. zero) THEN
                        DataItem(n_2d,m_2d,ntype_2d) = diag_val(js)
                     END IF

!Levenberg-like offset - do NOT apply here if applied in colscaling routine
                     IF (ntype_2d .GE. lamtype) THEN
                        DataItem(n_2d,m_2d,ntype_2d) =
     &                     1.0001_dp*DataItem(n_2d,m_2d,ntype_2d)
                     END IF

                     IF (lswap2disk) THEN
                        CALL WriteDAItem_SEQ(DataItem)
                     ELSE
                        block_diag(:,:,:,n_2d,m_2d,ntype_2d,js) =
     &                     DataItem
                     END IF

                 !block_plus(js-1) == b
                     IF (js .GT. 1) THEN
                        DataItem = (gc(js-1,:,:,:) -
     &                              gc_save(js-1,:,:,:))/hj
                     END IF
!no coupling of ALL fixed bdy forces to ANY r,z bdy values
                     IF (lswap2disk) THEN
                        CALL WriteDAItem_SEQ(DataItem)
                     ELSE IF (js .GT. 1) THEN
                        block_plus(:,:,:,n_2d,m_2d,ntype_2d,js-1) =
     &                     DataItem
                     END IF

                  END DO SKIP3_MESH
               END DO MESH_3PT
            END DO N2D
         END DO M2D
      END DO NTYPE2D

      l_edge = .FALSE.
      DEALLOCATE(DataItem)

      END SUBROUTINE sweep3_blocks

      SUBROUTINE blk3d_factor(a, bm1, bp1, ipiv, mblk, nblocks)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(sp), PARAMETER :: zero = 0, one = 1    
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(out) :: ipiv(mblk,nblocks)
      REAL(sp), TARGET, DIMENSION(mblk,mblk,nblocks), INTENT(inout) :: 
     &   a, bm1, bp1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
!      INTEGER :: ibuph, incnow, irecl, incbu, iunit=102, ndisk
      INTEGER :: k, k1, ier
      INTEGER, POINTER :: ipivot(:)
      REAL(sp), POINTER :: amat(:,:), bmat(:,:), cmat(:,:)
      REAL(sp), ALLOCATABLE, DIMENSION(:,:) :: temp
C-----------------------------------------------
c  modified (June, 2003, ORNL):         S. P. Hirshman
c-----------------------------------------------------------------------
c
c  this subroutine solves for the Q factors of a block-tridiagonal system of equations.
c
c-----------------------------------------------------------------------
c  INPUT
c  mblk                : block dimension (elements in a block=mblk X mblk)
c  nblocks             : number of blocks
c  a                   : diagonal blocks
c  bp1, bm1            : lower, upper blocks (see equation below)
c
c  OUTPUT
c  ipiv                : pivot elements for kth block
c  a                   : a-1 LU factor blocks
c  bm1                 : q = a-1 * bm1 matrix
c
c  LOCAL VARIABLES
c  iunit               : unit number for block-tridiagonal solution disk file.
c
c  solutions are indexed in m-n fourier-space, legendre-space. the tri-diagonal
c  equation is:
c
c           bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) = source(l)
c
c     GENERAL SOLUTION SCHEME APPLIED TO EACH BLOCK ROW (INDEX L)
c
c     1. Start from row N and solve for x(N) in terms of x(N-1):
c
c        x(N) = -q(N)*x(N-1) + r(N)
c
c        q(N) =  a(N)[-1] * bm1;    r(N) = a(N)[-1] * s(N)
c
c        where a(N)[-1] is the inverse of a(N)
c
c     2. Substitute for lth row to get recursion equation fo q(l) and r(l):
c
c        x(l) = -q(l)*x(l-1) + r(l), in general, where:
c
c        q(l) = (a(l) - bp1(l)*q(l+1))[-1] * bm1(l)
c
c        qblk(l) == (a(l) - bp1(l) * q(l+1))[-1] on return
c
c        r(l) = (a(l) - bp1(l)*q(l+1))[-1] * (s(l) - bp1(l)*r(l+1))
c
c     3. At row l = 1, bm1(1) = 0 and get an equation for x(1) corresponding to q(1) = 0:
c
c        x(1) = r(1)
c
c     4. Finally, can back-solve for x(N) in terms of x(N-1) from eqn.(1) above
c
c
c     NUMERICAL IMPLEMENTATION (USING LAPACK ROUTINES)
c
c     1. CALL dgetrf:   Perform LU factorization of diagonal block (A) - faster than sgefa
c     2. CALL dgetrs:   With multiple (mblk) right-hand sides, to do block inversion
c                         operation, A X = B  (stores result in B; here B is a matrix)
c

c  main loop. load and process (backwards) block-rows nblocks to 1. 


      BLOCKS: DO k = nblocks, 1, -1
!
!     Compute (and save) qblk(k) = ablk(k)[-1] * bml
!
         amat => a(:,:,k);  ipivot => ipiv(:,k)
         CALL dgetrf(mblk, mblk, amat, mblk, ipivot, ier)
         IF (ier .ne. 0) GOTO 200
         IF (k .eq. 1) EXIT
          
         bmat => bm1(:,:,k)
         CALL dgetrs('n', mblk, mblk, amat, mblk, ipivot, bmat, mblk,
     &               ier)

         IF (ier .ne. 0) GOTO 305

!
!      Update effective diagonal "a" matrix. Use dgemm: faster AND doesn't overflow normal stack
!
         k1 = k-1 
         amat => bp1(:,:,k1)
         cmat => a(:,:,k1)
!         cmat = cmat - MATMUL(amat, bmat)
         CALL dgemm('N','N',mblk,mblk,mblk,-one,amat,mblk,
     1              bmat, mblk, one, cmat, mblk)

      END DO BLOCKS

!
!     COMPUTE TRANSPOSES HERE, SINCE REPEATEDLY CALLING MATMUL OPERATION
!     X*At IS FASTER THAN A*X DUE TO UNIT STRIDE
!
      ALLOCATE (temp(mblk,mblk), stat=k)
      IF (k .ne. 0) STOP 'Allocation error in blk3d_factor!'

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

c  error returns. ------------------------------------------------------

  200 CONTINUE
!          < 0:  if info = -i, the i-th argument had an illegal value
!          > 0:  if info = i, u(i,i) is exactly zero. the factorization
      WRITE (6,1000) k
      IF (ier < 0) THEN
         WRITE (6,1001) ier
      END IF
      IF (ier > 0) THEN
         WRITE (6,1002) ier
      END IF
      STOP
  305 CONTINUE
      WRITE (6, 1003) ier
      STOP


  400 CONTINUE

1000  FORMAT(2x,'Error factoring matrix in blk3d: block = ',i4)
1001  FORMAT(i4,'th argument has illegal value')
1002  FORMAT(i4,'th diagonal factor exactly zero')
1003  FORMAT(2/' BLK3D:   error detected:   ier =',i4,2/)

      END SUBROUTINE blk3d_factor

 
      SUBROUTINE blk3d_factor_swp(ipiv, mblk, nblocks)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(sp), PARAMETER :: zero = 0, one = 1    
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(out) :: ipiv(mblk,nblocks)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(sp), ALLOCATABLE, DIMENSION(:,:) :: 
     &                       amat, bmat, cmat, temp
      INTEGER :: k, k1, ier
      INTEGER, POINTER :: ipivot(:)
C-----------------------------------------------
c  modified (June, 2003, ORNL):         S. P. Hirshman
c  modified (June, 2007, ORNL), added lswap2disk logic
c-----------------------------------------------------------------------
c
c  this subroutine solves for the Q factors of a block-tridiagonal system of equations.
c  see blk3d_factor for more details
c
c-----------------------------------------------------------------------
c  INPUT
c  mblk                : block dimension (elements in a block=mblk X mblk)
c  nblocks             : number of blocks
c
c  OUTPUT
c  ipiv                : pivot elements for kth block
c
c  LOCAL VARIABLES
c  iunit               : unit number for block-tridiagonal solution disk file.
c
c  solutions are indexed in m-n fourier-space, legendre-space. the tri-diagonal
c  equation is:
c
c           bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) = source(l)
c
c
c     NUMERICAL IMPLEMENTATION (USING LAPACK ROUTINES)
c
c     1. CALL dgetrf:   Perform LU factorization of diagonal block (A) - faster than sgefa
c     2. CALL dgetrs:   With multiple (mblk) right-hand sides, to do block inversion
c                         operation, A X = B  (stores result in B; here B is a matrix)
c

c  main loop. load and process (backwards) block-rows nblocks to 1. 

!     CHANGE Direct Access Record length to block size (from individual rows)
      CALL ChangeDAFileParams(mblk**2, mblk**2, 3, ScratchFile, nblocks)

      ALLOCATE(amat(mblk,mblk), bmat(mblk,mblk), cmat(mblk,mblk),
     &         temp(mblk,mblk), stat=ier)
      IF (ier .ne. 0) STOP 'Allocation error in blk3d_factor_swp!' 

      CALL ReadDAItem2(temp, nblocks, bldia)

      BLOCKS: DO k = nblocks, 1, -1
!
!     Compute (and save) qblk(k) = ablk(k)[-1] * bml
!
         amat = temp
         ipivot => ipiv(:,k)
         CALL dgetrf(mblk, mblk, amat, mblk, ipivot, ier)
         IF (ier .ne. 0) GOTO 200
!CONFIRM READ-WRITE ALLOWED...OK FOR DA Files!
         CALL WriteDAItem_RA(amat, k, bldia, 1)

         IF (k .eq. 1) EXIT
          
         CALL ReadDAItem2(bmat, k, blmin)
         CALL dgetrs('n', mblk, mblk, amat, mblk, ipivot, bmat, mblk,
     &               ier)
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
         CALL dgemm('N','N',mblk,mblk,mblk,-one,amat,mblk, bmat, mblk,
     &              one, temp, mblk)
         cmat = TRANSPOSE(amat)
         CALL WriteDAItem_RA(cmat, k1, blpls, 1)

      END DO BLOCKS

      GOTO 400

c  error returns. ------------------------------------------------------

  200 CONTINUE
!          < 0:  if info = -i, the i-th argument had an illegal value
!          > 0:  if info = i, u(i,i) is exactly zero. the factorization
      WRITE (6,1000) k
      IF (ier < 0) THEN
         WRITE (6,1001) ier
      END IF
      IF (ier > 0) THEN
         WRITE (6,1002) ier
      END IF
      STOP
  305 CONTINUE
      WRITE (6, 1003) ier
      STOP

  400 CONTINUE

      DEALLOCATE(amat, bmat, cmat, temp, stat=ier)

      CALL CloseDAFile

1000  FORMAT(2x,'Error factoring matrix in blk3d: block = ',i4)
1001  FORMAT(i4,'th argument has illegal value')
1002  FORMAT(i4,'th diagonal factor exactly zero')
1003  FORMAT(2/' BLK3D:   error detected:   ier =',i4,2/)

      END SUBROUTINE blk3d_factor_swp

      SUBROUTINE blk3d_factor_swp2(ipiv, mblk, nblocks)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(sp), PARAMETER :: zero = 0, one = 1    
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(out) :: ipiv(mblk,nblocks)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(sp), ALLOCATABLE, DIMENSION(:,:) :: 
     1                       amat, bmat, cmat, temp
      INTEGER :: k, k1, ier
      INTEGER, POINTER :: ipivot(:)
C-----------------------------------------------
c  modified (June, 2003, ORNL):         S. P. Hirshman
c  modified (June, 2007, ORNL), added lswap2disk logic
c  modified (July, 2007, ORNL/VA State):Helen Yang - forward sweep for speed
c-----------------------------------------------------------------------
c
c  this subroutine solves for the Q factors of a block-tridiagonal system of equations.
c  see blk3d_factor for more details
c
c-----------------------------------------------------------------------
c  INPUT
c  mblk                : block dimension (elements in a block=mblk X mblk)
c  nblocks             : number of blocks
c
c  OUTPUT
c  ipiv                : pivot elements for kth block
c
c  LOCAL VARIABLES
c  iunit               : unit number for block-tridiagonal solution disk file.
c
c  solutions are indexed in m-n fourier-space, legendre-space. the tri-diagonal
c  equation is:
c
c           bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) = source(l)
c
c
c     NUMERICAL IMPLEMENTATION (USING LAPACK ROUTINES)
c
c     1. CALL dgetrf:   Perform LU factorization of diagonal block (A) - faster than sgefa
c     2. CALL dgetrs:   With multiple (mblk) right-hand sides, to do block inversion
c                         operation, A X = B  (stores result in B; here B is a matrix)
c

c  main loop. load and process (backwards) block-rows nblocks to 1. 

!     CHANGE Direct Access Record length to block size (from individual rows)
      CALL ChangeDAFileParams(mblk**2, mblk**2, 3, ScratchFile, nblocks)

      ALLOCATE(amat(mblk,mblk), bmat(mblk,mblk), cmat(mblk,mblk),
     &         temp(mblk,mblk), stat=ier)
      IF (ier .ne. 0) STOP 'Allocation error in blk3d_factor_swp!' 

      CALL ReadDAItem2(temp, 1, bldia)

      BLOCKS: DO k = 1,nblocks
!
!     Compute (and save) qblk(k) = ablk(k)[-1] * bml
!
         amat = temp
         ipivot => ipiv(:,k)
         CALL dgetrf (mblk, mblk, amat, mblk, ipivot, ier)
         IF (ier .ne. 0) GOTO 200
!CONFIRM READ-WRITE ALLOWED...OK FOR DA Files!
         CALL WriteDAItem_RA(amat, k, bldia, 1)

         IF (k .eq. nblocks) EXIT
          
         CALL ReadDAItem2(bmat, k, blpls)
         CALL dgetrs('n', mblk, mblk, amat, mblk, ipivot, bmat, mblk,
     &               ier)
         IF (ier .ne. 0) GOTO 305
!
!     COMPUTE TRANSPOSES HERE (and for cmat below), SINCE REPEATEDLY CALLING MATMUL OPERATION
!     X*At IS FASTER THAN A*X DUE TO UNIT STRIDE
!
         temp = TRANSPOSE(bmat)
         CALL WriteDAItem_RA(temp, k, blpls, 1)

!
!      Update effective diagonal "a" matrix. Use dgemm: faster AND doesn't overflow normal stack
!
         k1 = k + 1
         CALL ReadDAItem2(amat, k1, blmin)
         CALL ReadDAItem2(temp, k1, bldia)
!         temp = temp - MATMUL(amat, bmat)
         CALL dgemm('N','N',mblk,mblk,mblk,-one,amat,mblk, bmat, mblk,
     &              one, temp, mblk)
         cmat = TRANSPOSE(amat)
         CALL WriteDAItem_RA(cmat, k1, blmin, 1)

      END DO BLOCKS

      GOTO 400

c  error returns. ------------------------------------------------------

  200 CONTINUE
!          < 0:  if info = -i, the i-th argument had an illegal value
!          > 0:  if info = i, u(i,i) is exactly zero. the factorization
      WRITE (6,1000) k
      IF (ier < 0) THEN
         WRITE (6,1001) ier
      END IF
      IF (ier > 0) THEN
         WRITE (6,1002) ier
      END IF
      STOP
  305 CONTINUE
      WRITE (6,1003) ier
      STOP

  400 CONTINUE

      DEALLOCATE(amat, bmat, cmat, temp, stat=ier)

      CALL CloseDAFile

1000  FORMAT(2x,'Error factoring matrix in blk3d: block = ',i4)
1001  FORMAT(i4,'th argument has illegal value')
1002  FORMAT(i4,'th diagonal factor exactly zero')
1003  FORMAT(2/' BLK3D:   error detected:   ier =',i4,2/)

      END SUBROUTINE blk3d_factor_swp2


      SUBROUTINE blk3d_slv(ablk, qblk, bp1, source, 
     1                     ipiv, mblk, nblocks)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE stel_kinds
!      USE safe_open_mod
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(in) :: ipiv(mblk,nblocks)
      REAL(sp), TARGET, DIMENSION(mblk,mblk,nblocks), INTENT(in) :: 
     1                           ablk, qblk, bp1
      REAL(dp), DIMENSION(mblk,nblocks), INTENT(inout) 
     1                  :: source
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, POINTER  :: ipivot(:)
      INTEGER :: k, ier
      REAL(sp), POINTER :: amat(:,:)          !, x1(:), y1(:)
      REAL(sp) :: source_sp(mblk)
C-----------------------------------------------
c  modified (June, 2003, ORNL):         S. P. Hirshman
c-----------------------------------------------------------------------
c
c  this subroutine solves a block-tridiagonal system of equations, using 
c  the ABLK, QBLK factors from blk3d_factor,
c
c-----------------------------------------------------------------------
c  INPUT
c  mblk                : block size
c  nblocks             : number of blocks
c  bp1                 : upper blocks (see equation below)
c  ipiv                : pivot elements for kth block
c  ablk                : a-1 blocks
c  qblk                : q = a-1 * bm1
c  source              : input right side
c
c  OUTPUT
c  source              : Solution x of A x = source
c 
c  LOCAL VARIABLES
c  iunit               : unit number for block-tridiagonal solution disk file.
c
c  solutions are indexed in m-n fourier-space, legendre-space. the tri-diagonal
c  equation is:
c
c           bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) = source(l)
c
c     GENERAL SOLUTION SCHEME APPLIED TO EACH BLOCK ROW (INDEX L)
c
c     1. Start from row N and solve for x(N) in terms of x(N-1):
c
c        x(N) = -q(N)*x(N-1) + r(N)
c
c        q(N) =  a(N)[-1] * bm1;    r(N) = a(N)[-1] * s(N)
c
c        where a(N)[-1] is the inverse of a(N)
c
c     2. Substitute for lth row to get recursion equation fo q(l) and r(l):
c
c        x(l) = -q(l)*x(l-1) + r(l), in general, where:
c
c        q(l) = (a(l) - bp1(l)*q(l+1))[-1] * bm1(l)
c
c        qblk(l) == (a(l) - bp1(l) * q(l+1))[-1] on return
c
c        r(l) = (a(l) - bp1(l)*q(l+1))[-1] * (s(l) - bp1(l)*r(l+1))
c
c     3. At row l = 1, bm1(1) = 0 and get an equation for x(1) corresponding to q(1) = 0:
c
c        x(1) = r(1)
c
c     4. Finally, can back-solve for x(N) in terms of x(N-1) from eqn.(1) above
c
c
c     NUMERICAL IMPLEMENTATION (USING LAPACK ROUTINES)
c
c     1. CALL dgetrs:   With single right hand side (source) to solve A x = b (b a vector)
c                         Faster than dgesl
!      ndisk = mblk*mblk

c  main loop. load and process (backwards) block-rows nblocks to 1. 
!  note: about equal time is spent in calling dgetrs and in performing
!  the two loop sums: on ibm-pc, 2 s (trs) vs 3 s (sums); on linux (logjam),
!  2.4 s (trs) vs 3 s (sums).

!
!     Back-solve for modified sources first
!
      BLOCKS: DO k = nblocks, 1, -1

         source_sp = source(:,k)
         ipivot => ipiv(:,k);   amat => ablk(:,:,k)
         CALL dgetrs('n', mblk, 1, amat, mblk,
     1                    ipivot, source_sp, mblk, ier)
         source(:,k) = source_sp

         IF (ier .ne. 0) GOTO 305
         IF (k .eq. 1) EXIT

!
!        NOTE: IN BLK3D_FACTOR, BP1 AND BM1 WERE TRANSPOSED (AND STORED)
!        TO MAKE FIRST INDEX FASTEST VARYING IN THE FOLLOWING MATMUL OPS
!
         amat => bp1(:,:,k-1)
         source(:,k-1) = source(:,k-1) - MATMUL(source(:,k),amat)  !USE THIS FORM IF TRANSPOSED bp1
!         source(:,k-1) = source(:,k-1) - MATMUL(amat,source(:,k))  !UNTRANSPOSED FORM
!         x1 => source(:,k);  y1 => source(:,k-1)
!         CALL dgemv('T',mblk,mblk,-one,amat,mblk,x1,1,
!     1                  one,y1,1)

      END DO BLOCKS
!
!  forward (back-substitution) solution sweep for block-rows k = 2 to nblocks
!  now, source contains the solution vector
!
      DO k = 2, nblocks
         amat => qblk(:,:,k)
         source(:,k) = source(:,k) - MATMUL(source(:,k-1),amat)  !USE THIS FORM IF TRANSPOSED qblk
!         source(:,k) = source(:,k) - MATMUL(amat,source(:,k-1))  !UNTRANSPOSED FORM
!         x1 => source(:,k-1);  y1 => source(:,k)
!         CALL dgemv('T',mblk,mblk,-one,amat,mblk,x1,1,
!     1                  one,y1,1)

      END DO

      GOTO 400

c  error returns. ------------------------------------------------------

  305 CONTINUE
      WRITE (6, '(2/a,i4,2/)') ' BLK3D:   error detected:   ier =',
     1   ier
      STOP

  400 CONTINUE

      END SUBROUTINE blk3d_slv

 
      SUBROUTINE blk3d_slv_swp(source, ipiv, mblk, nblocks)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE stel_kinds
!      USE safe_open_mod
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(in) :: ipiv(mblk,nblocks)
      REAL(dp), DIMENSION(mblk,nblocks), INTENT(inout) 
     1                  :: source
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(sp), ALLOCATABLE, DIMENSION(:,:) :: amat
      INTEGER, POINTER  :: ipivot(:)
      INTEGER :: k, k1, ier
      REAL(sp) :: source_sp(mblk)
C-----------------------------------------------
c  modified (June, 2003, ORNL):         S. P. Hirshman
c  modified (June, 2007, ORNL), added lswap2disk logic
c-----------------------------------------------------------------------
c
c  this subroutine solves a block-tridiagonal system of equations, using 
c  the ABLK, QBLK factors from blk3d_factor,
c  See blk3d_slv for details
c
c-----------------------------------------------------------------------
c  INPUT
c  mblk                : block size
c  nblocks             : number of blocks
c  ipiv                : pivot elements for kth block
c  source              : input right side
c
c  OUTPUT
c  source              : Solution x of A x = source
c 
c  LOCAL VARIABLES
c  iunit               : unit number for block-tridiagonal solution disk file.
c
c  the tri-diagonal equation is:
c
c           bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) = source(l)
c

c  main loop. load and process (backwards) block-rows nblocks to 1. 
!  note: about equal time is spent in calling dgetrs and in performing
!  the two loop sums: on ibm-pc, 2 s (trs) vs 3 s (sums); on linux (logjam),
!  2.4 s (trs) vs 3 s (sums).

      CALL OpenDAFile(mblk**2, mblk**2, 3, ScratchFile, iunit_dacess, 1) 

      ALLOCATE (amat(mblk,mblk),  stat=ier)
      IF (ier .ne. 0) STOP 'Allocation error in blk3d_slv_swp!' 
!
!     Back-solve for modified sources first
!
      BLOCKS: DO k = nblocks, 1, -1

         source_sp = source(:,k)
         ipivot => ipiv(:,k)
         CALL ReadDAItem2(amat, k, bldia)
         CALL dgetrs('n', mblk, 1, amat, mblk,
     1                    ipivot, source_sp, mblk, ier)
         source(:,k) = source_sp

         IF (ier .ne. 0) GOTO 305
         IF (k .eq. 1) EXIT

!
!        NOTE: IN BLK3D_FACTOR, BP1 AND BM1 WERE TRANSPOSED (AND STORED)
!        TO MAKE FIRST INDEX FASTEST VARYING IN THE FOLLOWING MATMUL OPS
!
         k1 = k-1
         CALL ReadDAItem2(amat, k1, blpls) 
         source(:,k1) = source(:,k1) - MATMUL(source(:,k),amat)  !USE THIS FORM IF TRANSPOSED bp1
!         source(:,k1) = source(:,k1) - MATMUL(amat,source(:,k))  !UNTRANSPOSED FORM
!         x1 => source(:,k);  y1 => source(:,k-1)
!         CALL dgemv('T',mblk,mblk,-one,amat,mblk,x1,1,
!     1                  one,y1,1)

      END DO BLOCKS
!
!  forward (back-substitution) solution sweep for block-rows k = 2 to nblocks
!  now, source contains the solution vector
!
      DO k = 2, nblocks

         CALL ReadDAItem2(amat, k, blmin)
         source(:,k) = source(:,k) - MATMUL(source(:,k-1),amat)  !USE THIS FORM IF TRANSPOSED qblk
!         source(:,k) = source(:,k) - MATMUL(amat,source(:,k-1))  !UNTRANSPOSED FORM
!         x1 => source(:,k-1);  y1 => source(:,k)
!         CALL dgemv('T',mblk,mblk,-one,amat,mblk,x1,1,
!     1                  one,y1,1)

      END DO

      WRITE(100,'(1p,6e14.4)') source(:,ns/2)

      GOTO 400

c  error returns. ------------------------------------------------------

  305 CONTINUE
      WRITE (6, '(2/a,i4,2/)') ' BLK3D:   error detected:   ier =',
     1   ier
      STOP

  400 CONTINUE

      CALL CloseDAFile

      END SUBROUTINE blk3d_slv_swp


      SUBROUTINE blk3d_slv_swp2(source, ipiv, mblk, nblocks)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE stel_kinds
!      USE safe_open_mod
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(in) :: ipiv(mblk,nblocks)
      REAL(dp), DIMENSION(mblk,nblocks), INTENT(inout) 
     1                  :: source
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(sp), ALLOCATABLE, DIMENSION(:,:) :: amat
      INTEGER, POINTER  :: ipivot(:)
      INTEGER :: k, k1, ier
      REAL(sp) :: source_sp(mblk)
C-----------------------------------------------
c  modified (June, 2003, ORNL):         S. P. Hirshman
c  modified (June, 2007, ORNL), added lswap2disk logic
c  modified (July, 2007, ORNL/VA State):Helen Yang - forward sweep for speed
c-----------------------------------------------------------------------
c
c  this subroutine solves a block-tridiagonal system of equations, using 
c  the ABLK, QBLK factors from blk3d_factor,
c  See blk3d_slv for details
c
c-----------------------------------------------------------------------
c  INPUT
c  mblk                : block size
c  nblocks             : number of blocks
c  bp1                 : upper blocks (see equation below)
c  ipiv                : pivot elements for kth block
c  ablk                : a-1 blocks
c  qblk                : q = a-1 * bm1
c  source              : input right side
c
c  OUTPUT
c  source              : Solution x of A x = source
c 
c  LOCAL VARIABLES
c  iunit               : unit number for block-tridiagonal solution disk file.
c
c  the tri-diagonal equation is:
c
c           bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) = source(l)
c

c  main loop. load and process (backwards) block-rows nblocks to 1. 
!  note: about equal time is spent in calling dgetrs and in performing
!  the two loop sums: on ibm-pc, 2 s (trs) vs 3 s (sums); on linux (logjam),
!  2.4 s (trs) vs 3 s (sums).

      CALL OpenDAFile(mblk**2, mblk**2, 3, ScratchFile, iunit_dacess, 1) 

      ALLOCATE (amat(mblk,mblk),  stat=ier)
      IF (ier .ne. 0) STOP 'Allocation error in blk3d_slv_swp!' 
!
!     Back-solve for modified sources first
!
      BLOCKS: DO k = 1, nblocks

         source_sp = source(:,k)
         ipivot => ipiv(:,k)
         CALL ReadDAItem2(amat, k, bldia)
         CALL dgetrs('n', mblk, 1, amat, mblk,
     1                    ipivot, source_sp, mblk, ier)
         source(:,k) = source_sp

         IF (ier .ne. 0) GOTO 305
         IF (k .eq. nblocks) EXIT

!
!        NOTE: IN BLK3D_FACTOR, BP1 AND BM1 WERE TRANSPOSED (AND STORED)
!        TO MAKE FIRST INDEX FASTEST VARYING IN THE FOLLOWING MATMUL OPS
!
         k1 = k+1
         CALL ReadDAItem2(amat, k1, blmin) 
         source(:,k1) = source(:,k1) - MATMUL(source(:,k),amat)  !USE THIS FORM IF TRANSPOSED bp1

      END DO BLOCKS
!
!  backward solution sweep for block-rows k = nblocks-1 to 1
!  now, source contains the solution vector
!
      DO k = nblocks-1, 1, -1

         CALL ReadDAItem2(amat, k, blpls)
         k1 = k+1
         source(:,k) = source(:,k) - MATMUL(source(:,k1),amat)  !USE THIS FORM IF TRANSPOSED qblk

      END DO

      WRITE(100,'(1p,6e14.4)') source(:,ns/2)

      GOTO 400

c  error returns. ------------------------------------------------------

  305 CONTINUE
      WRITE (6, '(2/a,i4,2/)') ' BLK3D:   error detected:   ier =',
     1   ier
      STOP

  400 CONTINUE

      CALL CloseDAFile

      END SUBROUTINE blk3d_slv_swp2


      SUBROUTINE compute_col_scaling_par
      USE xstuff, ONLY: pcol_scale
      USE blocktridiagonalsolver, ONLY: GetColSum, ParallelScaling
      USE parallel_vmec_module, ONLY: ToLastNtype, CopyLastNtype
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nsmin, nsmax
      REAL(dp), ALLOCATABLE  :: tmp(:)
      REAL(dp), ALLOCATABLE  :: colsum(:,:)
      REAL(dp), PARAMETER    :: levmarq_param = 1.E-6_dp
C-----------------------------------------------

!FOR NO COL SCALING - col-scaling not working well yet (8.1.17)
      pcol_scale = 1
      RETURN

!BE SURE TO TURN OFF LEV_MARQ SCALING IN SUBROUTINE sweep3_blocks_par
      nsmin = tlglob;  nsmax = trglob

      ALLOCATE (colsum(mblk_size,nsmin:nsmax))
      CALL GetColSum(colsum)
      CALL VectorCopyPar (colsum, pcol_scale)
      CALL ParallelScaling(levmarq_param,colsum)

      DEALLOCATE(colsum)

!Convert to internal PARVMEC format
      ALLOCATE (tmp(ntmaxblocksize*ns))
      CALL tolastntype(pcol_scale,tmp)
      CALL copylastntype(tmp,pcol_scale)
      DEALLOCATE(tmp)

      END SUBROUTINE compute_col_scaling_par

      SUBROUTINE compute_col_scaling
      USE xstuff, ONLY: col_scale
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------

!FOR NO COL SCALING
      col_scale = 1

      END SUBROUTINE compute_col_scaling

      SUBROUTINE VectorCopyPar (colsum, colscale)
      USE blocktridiagonalsolver, ONLY: rank
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN)  :: colsum(mblk_size,tlglob:trglob)
      REAL(dp), INTENT(OUT) :: colscale(mblk_size,ns)

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER               :: js, M1
      INTEGER               :: MPI_STAT(MPI_STATUS_SIZE)
!-----------------------------------------------
 
      DO js = tlglob, trglob 
        colscale(:,js) = colsum(:,js)
      END DO

      M1 = mblk_size

! Get left boundary elements (tlglob-1)
      IF (rank.LT.nranks-1) THEN
        CALL MPI_Send(colsum(:,trglob),M1,MPI_REAL8,    
     1                rank+1,1,NS_COMM,MPI_ERR)
      END IF
      IF (rank.GT.0) THEN
        CALL MPI_Recv(colscale(:,tlglob-1),M1,      
     1                MPI_REAL8,rank-1,1,NS_COMM,MPI_STAT,MPI_ERR)
      END IF

! Get right boundary elements (trglob+1)
      IF (rank.GT.0) THEN
        CALL MPI_Send(colsum(:,tlglob),M1,MPI_REAL8,
     1                rank-1,1,NS_COMM,MPI_ERR)
      END IF
      IF (rank.LT.nranks-1) THEN
        CALL MPI_Recv(colscale(:,trglob+1),M1,MPI_REAL8,
     1                rank+1,1,NS_COMM,MPI_STAT,MPI_ERR)
      END IF

      END SUBROUTINE VectorCopyPar


      SUBROUTINE free_mem_precon
      INTEGER :: istat

      istat=0
      IF (ALLOCATED(block_diag)) 
     1    DEALLOCATE (block_diag, block_plus, block_mins, stat=istat)
      IF (istat .ne. 0) STOP 'Deallocation error-1 in free_mem_precon'
      
      istat=0
      IF (ALLOCATED(block_diag_sw))
     1   DEALLOCATE (block_diag_sw, block_plus_sw, block_mins_sw, 
     2                stat=istat)
      IF (istat .ne. 0) STOP 'Deallocation error-2 in free_mem_precon'

      istat=0
      IF (ALLOCATED(ipiv_blk)) DEALLOCATE (ipiv_blk, stat=istat)     
      IF (istat .ne. 0) STOP 'Deallocation error-3 in free_mem_precon'

      END SUBROUTINE free_mem_precon

      END MODULE precon2d
