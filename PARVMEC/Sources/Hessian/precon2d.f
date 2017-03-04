      MODULE precon2d
      USE stel_kinds, ONLY: rprec2 => rprec, dp2 => dp
      USE vmec_dim
      USE vmec_params
      USE vparams, ONLY: nthreed, zero
      USE vmec_input, ONLY: ntor, nzeta, lfreeb, lasym 
      USE timer_sub
      USE safe_open_mod
      USE directaccess
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
!     INTEGER, PARAMETER :: sp = KIND(1.0)
      INTEGER, PARAMETER :: sp = rprec2
      INTEGER :: ntyptot
      INTEGER :: ntype_2d, m_2d, n_2d
      INTEGER, ALLOCATABLE :: ipiv_blk(:,:)
      REAL(sp), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: 
     1    block_diag, block_plus, block_mins,
     3    block_dsave, block_msave, block_psave
      REAL(sp), ALLOCATABLE, DIMENSION(:,:,:,:,:,:) ::
     1    block_diag_sw, block_plus_sw, block_mins_sw
   
      REAL(rprec2), DIMENSION(:,:,:,:), ALLOCATABLE :: gc_save
      REAL(rprec2), ALLOCATABLE, DIMENSION(:,:,:) :: 
     1   r1s_save, rus_save, rvs_save, rcons_save,
     2   z1s_save, zus_save, zvs_save, zcons_save,
     3   lus_save, lvs_save,
     1   r1a_save, rua_save, rva_save, rcona_save,
     2   z1a_save, zua_save, zva_save, zcona_save,
     3   lua_save, lva_save
      REAL(rprec2) :: ctor_prec2d
      INTEGER :: ictrl_prec2d
      LOGICAL :: lHess_exact = .false.,
     1           l_backslv = .false.,
     2           l_comp_prec2D = .true.
      PRIVATE :: swap_forces, reswap_forces, gc_save, block_dsave,
     1           block_msave, block_psave

!
!     Direct-Access (swap to disk) stuff
!      CHARACTER(LEN=3)   :: FlashDrive ="F:\"
      CHARACTER(LEN=3)   :: FlashDrive =""
      CHARACTER(LEN=128) :: ScratchFile=""
      INTEGER, PARAMETER :: blmin=1, bldia=2, blpls=3
      REAL(sp), ALLOCATABLE, DIMENSION(:,:,:) :: DataItem
      INTEGER            :: iunit_dacess=10
      LOGICAL :: lswap2disk = .FALSE.                                 !Set internally if blocks do not fit in memory

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
!                   in practice. Set this true primarily for debugging purposes (check Ap ~ I in axb 
!                   routine, for example, in GMRes module)
!

      CONTAINS

      SUBROUTINE swap_forces(gc, temp, mblk, nblocks)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: mblk, nblocks
      REAL(rprec2), DIMENSION(nblocks,mblk), INTENT(in)  :: gc
      REAL(rprec2), DIMENSION(mblk,nblocks), INTENT(out) :: temp
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
      REAL(rprec2), DIMENSION(nblocks,mblk), INTENT(inout)  :: gc
      REAL(rprec2), DIMENSION(mblk,nblocks), INTENT(in) :: temp
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
      REAL(rprec2), DIMENSION(ns,0:ntor,0:mpol1,ntyptot) :: gc
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mblk, m, n, js, ntype, istat
      REAL(rprec2), ALLOCATABLE, DIMENSION(:,:) :: temp
      REAL(rprec2) :: t1, error
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

      SUBROUTINE compute_blocks (xc, xcdot, gc)
      USE realspace, ONLY: sqrts
      USE vmec_main, ONLY: iter2, lthreed, iequi, ncurr
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec2),DIMENSION(ns,0:ntor,0:mpol1,3*ntmax) :: xc, gc, xcdot
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec2), PARAMETER :: p5 = 0.5_rprec2
      INTEGER :: m, n, i, ntype, istat,  
     1           mblk, nmax_jog, ibsize, iunit
      REAL(rprec2) :: time_on, time_off, bsize, tprec2don, tprec2doff
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
      IF (l_backslv .and. sp.ne.rprec2) STOP 'Should set sp = rprec2!'

      ntyptot = SIZE(gc,4)
      mblk = ntyptot*mnsize
      nmax_jog = ntyptot

      bsize = REAL(mblk*mblk, dp2)*3*ns*KIND(block_diag)
      IF (bsize .gt. HUGE(mblk)) THEN
        WRITE (6, *) ' bsize: ', bsize, ' exceeds HUGE(int): ', 
     1                 HUGE(mblk)
        WRITE (6, *) ' Blocks will be written to disk.'
        lswap2disk = .TRUE.
      ELSE
        lswap2disk = .FALSE.
      END IF

      IF (bsize .lt. 1.E6_dp2) THEN
          ibsize = bsize/1.E1_dp2
          label = " Kb"
      ELSE IF (bsize .lt. 1.E9_dp2) THEN
          ibsize = bsize/1.E4_dp2
          label = " Mb"
      ELSE
          ibsize = bsize/1.E7_dp2
          label = " Gb"
      END IF

      DO i = 1,2
         IF (i .eq. 1) iunit = 6
         IF (i .eq. 2) iunit = nthreed
         WRITE (iunit, '(/,2x,a,i5,a,/,2x,a,i5,a)') 
     1         'Initializing 2D block preconditioner at ', iter2,
     2         ' iterations',
     3         'Estimated time to compute Hessian = ',
     4         3*nmax_jog*mnsize,' VMEC time steps'
         WRITE (iunit, '(2x,a,i4,a,f12.2,a)') 'Block dim: ', mblk,
     1         '^2  Preconditioner size: ', REAL(ibsize)/100, 
     2         TRIM(label)
      END DO
!
!     COMPUTE AND STORE BLOCKS (MN X MN) FOR PRECONDITIONER
!
      CALL second0(time_on)

      ALLOCATE (gc_save(ns,0:ntor,0:mpol1,ntyptot), stat=istat)
      IF (istat .ne. 0) 
     1    STOP 'Allocation error: gc_save in compute_blocks'

      IF (ALLOCATED(block_diag)) THEN
          DEALLOCATE (block_diag, block_plus, block_mins, stat=istat)
          IF (istat .ne. 0) STOP 'Deallocation error in compute blocks'
      ELSE IF (ALLOCATED(block_diag_sw)) THEN
         DEALLOCATE (block_diag_sw, block_plus_sw, block_mins_sw, 
     1                stat=istat)
          IF (istat .ne. 0) STOP 'Deallocation error in compute blocks'
      END IF

      ALLOCATE (
     1 block_diag(0:ntor,0:mpol1,ntyptot,0:ntor,0:mpol1,ntyptot,ns),
     2 block_plus(0:ntor,0:mpol1,ntyptot,0:ntor,0:mpol1,ntyptot,ns),
     3 block_mins(0:ntor,0:mpol1,ntyptot,0:ntor,0:mpol1,ntyptot,ns),
     4            stat=istat)

      lswap2disk = (istat .ne. 0)

!FOR DEBUGGING, SET THIS TO TRUE
!      lswap2disk = .TRUE.
      
      IF (lswap2disk) THEN
         WRITE (6,'(a,i4,a)') '  Allocation error(1) = ', istat,
     1               ': Not enough memory in compute_blocks'
         WRITE (6,*) ' Writing blocks to disk file'
         ALLOCATE (
     1 block_diag_sw(0:ntor,0:mpol1,ntyptot,0:ntor,0:mpol1,ntyptot),
     2 block_plus_sw(0:ntor,0:mpol1,ntyptot,0:ntor,0:mpol1,ntyptot),
     3 block_mins_sw(0:ntor,0:mpol1,ntyptot,0:ntor,0:mpol1,ntyptot),
     4               stat=istat)

         nmax_jog = 3*ntmax
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
         IF (FlashDrive .ne. "") ScratchFile = FlashDrive // ScratchFile
         CALL OpenDAFile(mblk, mblk**2, 3, ScratchFile, iunit_dacess, 0) 
         ScratchFile = "PRCND2B.bin"
         IF (FlashDrive .ne. "") ScratchFile = FlashDrive // ScratchFile
         block_plus_sw = 0; block_mins_sw = 0; block_diag_sw = 0
      ELSE
         block_plus = 0; block_mins = 0; block_diag = 0
      END IF


      ALLOCATE(
     1 lus_save(ns*nzeta,ntheta3,0:1), lvs_save(ns*nzeta,ntheta3,0:1), 
     1 rus_save(ns*nzeta,ntheta3,0:1), rvs_save(ns*nzeta,ntheta3,0:1),
     2 r1s_save(ns*nzeta,ntheta3,0:1), rcons_save(ns*nzeta,ntheta3,0:1), 
     3 zus_save(ns*nzeta,ntheta3,0:1), zvs_save(ns*nzeta,ntheta3,0:1),
     4 z1s_save(ns*nzeta,ntheta3,0:1), zcons_save(ns*nzeta,ntheta3,0:1),
     5 stat=istat)
      IF (lasym .and. istat.eq.0) ALLOCATE(
     1 lua_save(ns*nzeta,ntheta3,0:1), lva_save(ns*nzeta,ntheta3,0:1), 
     1 rua_save(ns*nzeta,ntheta3,0:1), rva_save(ns*nzeta,ntheta3,0:1),
     2 r1a_save(ns*nzeta,ntheta3,0:1), rcona_save(ns*nzeta,ntheta3,0:1), 
     3 zua_save(ns*nzeta,ntheta3,0:1), zva_save(ns*nzeta,ntheta3,0:1),
     4 z1a_save(ns*nzeta,ntheta3,0:1), zcona_save(ns*nzeta,ntheta3,0:1), 
     5 stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error(2) in compute_blocks!'

!
!     GENERAL (SLOWER BY 2/3 THAN SYMMETRIC VERSION) METHOD: ASSUMES NO SYMMETRIES OF R, Z COEFFICIENTS
!
      CALL sweep3_blocks (xc, xcdot, gc, nmax_jog)
      ictrl_prec2d = 1                 !Signals funct3d (residue) to use block preconditioner

      DEALLOCATE(lus_save, lvs_save, r1s_save, rus_save, rvs_save,
     1           rcons_save, z1s_save, zus_save, zvs_save, zcons_save,
     2           stat = istat)
      IF (lasym) DEALLOCATE(
     1    lua_save, lva_save, r1a_save, rua_save, rva_save, rcona_save, 
     2    z1a_save, zua_save, zva_save, zcona_save, stat=istat)

      CALL second0(time_off)
      DO m = 1, 2
         IF (m .eq. 1) n = 6
         IF (m .eq. 2) n = nthreed
         WRITE (n,'(1x,a,f10.2,a)')' Time to compute blocks: ', 
     1      time_off - time_on,' s'
      END DO

!     SAVE ORIGINAL (UNFACTORED) BLOCKS FOR CHECKING INVERSE
!     IN L_BACKSLV=TRUE LOOP IN BLOCK_PRECOND
      IF (l_backslv) THEN
         ALLOCATE (
     1   block_dsave(0:ntor,0:mpol1,ntyptot,0:ntor,0:mpol1,ntyptot,ns),
     2   block_msave(0:ntor,0:mpol1,ntyptot,0:ntor,0:mpol1,ntyptot,ns),
     3   block_psave(0:ntor,0:mpol1,ntyptot,0:ntor,0:mpol1,ntyptot,ns),
     4   stat = istat)
         IF (istat .ne. 0) THEN
            WRITE (6,*) 'Allocation error in l_backslv block: stat = ',
     1                istat
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
      IF (ALLOCATED(ipiv_blk)) DEALLOCATE(ipiv_blk, stat=ntype)
      ALLOCATE (ipiv_blk(mblk,ns), stat=ntype)     
      IF (ntype .ne. 0) STOP 'Allocation error2 in block_precond'

      IF (lswap2disk) THEN
         CALL blk3d_factor_swp(ipiv_blk, mblk, ns)
      ELSE
         CALL blk3d_factor(block_diag, block_mins, block_plus,  
     1                     ipiv_blk, mblk, ns)
      END IF

      CALL second0(time_off)
      DO m = 1, 2
          IF (m .eq. 1) n = 6
          IF (m .eq. 2) n = nthreed
          WRITE (n,'(1x,a,f10.2,a)')' Time to factor blocks:  ', 
     1                       time_off - time_on,' s'
      END DO

      IF (.not.l_backslv) DEALLOCATE (gc_save)

      CALL second0(tprec2doff)

      timer(tprec2d) = timer(tprec2d) + (tprec2doff - tprec2don)

      END SUBROUTINE compute_blocks


      SUBROUTINE sweep3_blocks (xc, xcdot, gc, nmax_jog)
      USE vmec_main, ONLY: ncurr, r01, z01, lthreed
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec2), DIMENSION(ns,0:ntor,0:mpol1,ntyptot) :: xc, xcdot, 
     1                                                      gc
      INTEGER :: nmax_jog
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: jstart(3) = (/1,2,3/)
      REAL(rprec2), DIMENSION(:,:,:,:), ALLOCATABLE :: xstore
      REAL(rprec2) :: eps, hj, diag_val(ns), t1
!      INTEGER :: js, i, istat=0, mesh, lamtype, rztype
      INTEGER :: js, istat, mesh, lamtype, rztype
      LOGICAL, PARAMETER :: lscreen = .false.
C-----------------------------------------------
!
!     COMPUTE EVERY 3rd RADIAL POINT FOR EACH MESH STARTING AT js=1,2,3, RESPECTIVELY
!
      WRITE (6, *)
     1   " Using non-symmetric sweep to compute Hessian elements"

      eps = SQRT(EPSILON(eps))
      eps = eps/10
      rztype = 2*ntmax
      lamtype = rztype+1

      xcdot = 0
      n_2d = 0;  m_2d = 0

      ALLOCATE(DataItem(0:ntor,0:mpol1,1:3*ntmax), 
     1         xstore(ns,0:ntor,0:mpol1,ntyptot), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in sweep3_blocks'

      diag_val = 1

!
!     CALL FUNCT3D FIRST TIME TO STORE INITIAL UN-PRECONDITIONED FORCES
!     THIS WILL CALL LAMBLKS (SO FAR, ONLY IMPLEMENTED FOR lasym = false)
!
      ictrl_prec2d = 2                 !Signals funct3d that preconditioner is being initialized
!      STOP 'In precon2d: not parallelized yet'
      CALL funct3d (lscreen, istat)
      IF (istat .ne. 0) THEN
         PRINT *,' ier_flag = ', istat,
     1           ' in SWEEP3_BLOCKS call to funct3d'
         STOP
      ENDIF

      xstore = xc                      !chips is saved in xc(lsc00)
      ictrl_prec2d = 3                 !Signals funct3d that preconditioner is being computed

      gc_save = gc

!
!     PERFORM "JOGS" FOR EACH VARIABLE AT EVERY 3rd RADIAL POINT ACROSS MESH
!     FOR ntyp = (Rcc, Rss, Rsc, Rcs, Zsc, Zcs, Zcc, Zss)
!     AND EVERY n2d (toroidal mode index) and EVERY m2d (poloidal mode index)
!
!     Lambda jogs are implemented analytically for lasym=FALSE
!      

      NTYPE2D: DO ntype_2d = 1, nmax_jog
         IF (ntype_2d .lt. lamtype) THEN
            hj = eps * MAX(ABS(r01(ns)), ABS(z01(ns)))
         ELSE
            hj = eps
         END IF

         M2D: DO m_2d = 0, mpol1
            N2D: DO n_2d = 0, ntor
            MESH_3PT: DO mesh = 1,3
!              APPLY JOG
               DO js = jstart(mesh), ns, 3
                  xc(js,n_2d,m_2d,ntype_2d) = 
     1            xstore(js,n_2d,m_2d,ntype_2d) + hj
                  xcdot(js,n_2d,m_2d,ntype_2d) = hj
               END DO

               istat = 0
               CALL funct3d (lscreen, istat)
               IF (istat .ne. 0) STOP 'Error computing Hessian jog!'

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
!              CLEAR JOG AND STORE BLOCKS FOR THIS JOG
               xc = xstore;  xcdot = 0
!
!              FOR OFF-DIAGONAL ELEMENTS, NEED TO ADJUST js INDICES +/- 1
!
               SKIP3_MESH: DO js = jstart(mesh), ns, 3

                  !block_mins(js+1) 
                  IF (js .lt. ns) THEN
                     DataItem = (gc(js+1,:,:,:)-gc_save(js+1,:,:,:))/hj
                    
!no coupling of ANY r,z fixed bdy forces to ALL xc values
                  IF (.not.lfreeb .and. js.eq.ns-1) 
     1               DataItem(:,:,1:rztype) = 0
                  ELSE IF (lswap2disk) THEN
                     DataItem = 0
                  END IF
                  IF (lswap2disk) THEN
                     !CALL WriteDAItem_RA(DataItem,js+1,blmin,index)
                     CALL WriteDAItem_SEQ(DataItem)
                  ELSE IF (js .lt. ns) THEN
                     block_mins(:,:,:,n_2d,m_2d,ntype_2d,js+1)=DataItem
                  END IF

                  !block_diag(js)
                  IF (js.eq.ns .and. .not.lfreeb) THEN
                     DataItem(:,:,1:rztype) = 0
                     IF (ntype_2d .lt. lamtype) THEN
                        DataItem(:,:,lamtype:) = 0
!     Set diagonal elements of R,Z at bdy (ns) to avoid singular Hessian
                        DataItem(n_2d,m_2d,ntype_2d) = 
     1                  (gc(ns-3,n_2d,m_2d,ntype_2d) - 
     2                   gc_save(ns-3,n_2d,m_2d,ntype_2d))/hj
                     ELSE
                        DataItem(:,:,lamtype:) =  (gc(ns,:,:,lamtype:)
     1                                  - gc_save(ns,:,:,lamtype:))/hj
                     END IF
                  ELSE
                     DataItem = (gc(js,:,:,:) - gc_save(js,:,:,:))/hj
                  END IF

!
!     Finite values for diagonal at JS=1, M>=1 (avoid hessian singularity)
!
                  IF (js.eq.1 .and. m_2d.gt.0) THEN
                     DataItem(n_2d,m_2d,ntype_2d) = 
     1               (gc(js+3,n_2d,m_2d,ntype_2d) - 
     2                gc_save(js+3,n_2d,m_2d,ntype_2d))/hj
                  END IF

!
!     Set diagonals for m=0,n  and m,n=0 sin modes (set diag_val to rcc value)
                  IF (n_2d.eq.0 .and. m_2d.eq.0 .and. ntype_2d.eq.rcc)
     1               diag_val(js) = DataItem(0,0,rcc)
                  IF (m_2d. eq. 0) THEN
                     IF (ntype_2d .eq. zsc+ntmax) 
     1                  DataItem(n_2d,0,ntype_2d) = diag_val(js)
                     IF (ntype_2d .eq. zsc+2*ntmax) THEN 
                        IF (ncurr.eq.0 .or. n_2d.ne.0)             !lsc(0,0) for ncurr=1 stores iotas
     1                     DataItem(n_2d,0,ntype_2d) = 1
                        IF (n_2d .eq. 0 .and. js.eq.1)
     1                     DataItem(n_2d,0,ntype_2d) = 1           !iotas(1), lamcc not evolved
                     END IF
                     IF (lthreed) THEN
                        IF (ntype_2d .eq. rss)
     1                     DataItem(n_2d,0,ntype_2d) = diag_val(js)
                        IF (js.eq.1 .and. jlam(0).gt.1 .and. 
     1                      ntype_2d.eq.zcs+2*ntmax)
     2                     DataItem(n_2d,0,ntype_2d) = 1
                        IF (lasym) THEN
                           IF (ntype_2d .eq. zss+ntmax)
     1                     DataItem(n_2d,0,ntype_2d) = diag_val(js)
                           IF (ntype_2d .eq. zss+2*ntmax)
     1                     DataItem(n_2d,0,ntype_2d) = 1
                        END IF
                     END IF
                     IF (lasym) THEN
                        IF (ntype_2d.eq.rsc)
     1                     DataItem(n_2d,0,ntype_2d) = diag_val(js)
                        IF (n_2d.eq.0 .and. ntype_2d.eq.(zcc+2*ntmax))      !m=0, n=0 lambda not evolved
     1                     DataItem(0,0,ntype_2d) = 1
                        IF (ntype_2d.eq.(zcc+2*ntmax) .and.                 !SPH120808
     1                      jlam(0).ne. 1 .and. js.eq.1)        
     1                     DataItem(n_2d,0,ntype_2d) = 1
                     END IF
                  END IF

                  IF (lthreed .and. n_2d.eq.0) THEN
                     IF (ntype_2d .eq. rss) 
     1                  DataItem(0,m_2d,ntype_2d) = diag_val(js)
                     IF (ntype_2d .eq. zcs+ntmax)
     1                  DataItem(0,m_2d,ntype_2d) = diag_val(js)
                     IF (ntype_2d .eq. zcs+2*ntmax)
     1                  DataItem(0,m_2d,ntype_2d) = 1
                     IF (lasym) THEN
                        IF (ntype_2d .eq. rcs)
     1                     DataItem(0,m_2d,ntype_2d) = diag_val(js)
                        IF (ntype_2d .eq. zss+ntmax)
     1                     DataItem(0,m_2d,ntype_2d) = diag_val(js)
                        IF (ntype_2d .eq. zss+2*ntmax)
     1                     DataItem(0,m_2d,ntype_2d) = 1
                     END IF
                  END IF

!SPH121912
                  IF (DataItem(n_2d,m_2d,ntype_2d) .eq. zero) THEN
                     DataItem(n_2d,m_2d,ntype_2d) = diag_val(js)
                     IF (diag_val(js) .eq. zero) STOP 'diag_val(js) = 0'
                  END IF

!                   DataItem(n_2d,m_2d,ntype_2d) = (1+1.E-4_dp)*
!     1             DataItem(n_2d,m_2d,ntype_2d)

                  IF (lswap2disk) THEN
                     !CALL WriteDAItem_RA(DataItem,js,bldia,index)
                     CALL WriteDAItem_SEQ(DataItem)
                  ELSE
                     block_diag(:,:,:,n_2d,m_2d,ntype_2d,js) = DataItem

!add this for gcz = 0 constraint for m=1, n!=0
               t1 = block_diag(n_2d,m_2d,ntype_2d,n_2d,m_2d,ntype_2d,js)
      IF (t1 .eq. zero .or. (js.eq.1 .and. ntype_2d .ge. lamtype)) THEN 
                        IF (diag_val(js) .eq. zero) diag_val(js) = 1
                   block_diag(n_2d,m_2d,ntype_2d,n_2d,m_2d,ntype_2d,js)=
     1                  diag_val(js)
                     ENDIF
                  END IF

                 !block_plus(js-1)
                 IF (js .gt. 1) THEN
                    DataItem = (gc(js-1,:,:,:)-gc_save(js-1,:,:,:))/hj
!no coupling of ALL fixed bdy forces to ANY r,z bdy values
                    IF (.not.lfreeb .and. js.eq.ns) THEN
                       IF (ntype_2d .le. rztype) DataItem = 0
                    END IF
                 ELSE IF (lswap2disk) THEN
                     DataItem = 0
                 END IF
                 IF (lswap2disk) THEN
                    !CALL WriteDAItem_RA(DataItem,js-1,blpls,index)
                    CALL WriteDAItem_SEQ(DataItem)
                 ELSE IF (js .gt. 1) THEN
                    block_plus(:,:,:,n_2d,m_2d,ntype_2d,js-1)=DataItem
                 END IF

               END DO SKIP3_MESH
               END DO MESH_3PT
            END DO N2D
         END DO M2D
      END DO NTYPE2D

      DEALLOCATE(DataItem, xstore)

      END SUBROUTINE sweep3_blocks


      SUBROUTINE blk3d_factor(a, bm1, bp1, ipiv, mblk, nblocks)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      IMPLICIT NONE
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
     1                       a, bm1, bp1
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
         CALL dgetrf (mblk, mblk, amat, mblk, ipivot, ier)
         IF (ier .ne. 0) GOTO 200
         IF (k .eq. 1) EXIT
          
         bmat => bm1(:,:,k)
         CALL dgetrs('n', mblk, mblk, amat, mblk, ipivot, 
     1                 bmat, mblk, ier)

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
      WRITE (6, '(2x,a,i4)') 'Error factoring matrix in blk3d: block = '
     1                        , k
      IF (ier < 0)WRITE (6,'(i4, a)') 
     1    ier, 'th argument has illegal value'
      IF (ier > 0)WRITE (6,'(i4, a)') 
     1    ier, 'th diagonal factor exactly zero'
      STOP
  305 CONTINUE
      WRITE (6, '(2/a,i4,2/)') ' BLK3D:   error detected:   ier =',
     1   ier
      STOP


  400 CONTINUE

      END SUBROUTINE blk3d_factor

 
      SUBROUTINE blk3d_factor_swp(ipiv, mblk, nblocks)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      IMPLICIT NONE
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
     1         temp(mblk,mblk), stat=ier)
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
         CALL dgetrs('n', mblk, mblk, amat, mblk, ipivot, 
     1                 bmat, mblk, ier)
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
         CALL dgemm('N','N',mblk,mblk,mblk,-one,amat,mblk,
     1              bmat, mblk, one, temp, mblk)
         cmat = TRANSPOSE(amat)
         CALL WriteDAItem_RA(cmat, k1, blpls, 1)

      END DO BLOCKS

      GOTO 400

c  error returns. ------------------------------------------------------

  200 CONTINUE
!          < 0:  if info = -i, the i-th argument had an illegal value
!          > 0:  if info = i, u(i,i) is exactly zero. the factorization
      WRITE (6, '(2x,a,i4)') 'Error factoring matrix in blk3d: block = '
     1                        , k
      IF (ier < 0)WRITE (6,'(i4, a)') 
     1    ier, 'th argument has illegal value'
      IF (ier > 0)WRITE (6,'(i4, a)') 
     1    ier, 'th diagonal factor exactly zero'
      STOP
  305 CONTINUE
      WRITE (6, '(2/a,i4,2/)') ' BLK3D:   error detected:   ier =',
     1   ier
      STOP

  400 CONTINUE

      DEALLOCATE(amat, bmat, cmat, temp, stat=ier)

      CALL CloseDAFile

      END SUBROUTINE blk3d_factor_swp


      SUBROUTINE blk3d_factor_swp2(ipiv, mblk, nblocks)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      IMPLICIT NONE
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
     1         temp(mblk,mblk), stat=ier)
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
         CALL dgetrs('n', mblk, mblk, amat, mblk, ipivot, 
     1                 bmat, mblk, ier)
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
         k1 = k+1 
         CALL ReadDAItem2(amat, k1, blmin)
         CALL ReadDAItem2(temp, k1, bldia)
!         temp = temp - MATMUL(amat, bmat)
         CALL dgemm('N','N',mblk,mblk,mblk,-one,amat,mblk,
     1              bmat, mblk, one, temp, mblk)
         cmat = TRANSPOSE(amat)
         CALL WriteDAItem_RA(cmat, k1, blmin, 1)

      END DO BLOCKS

      GOTO 400

c  error returns. ------------------------------------------------------

  200 CONTINUE
!          < 0:  if info = -i, the i-th argument had an illegal value
!          > 0:  if info = i, u(i,i) is exactly zero. the factorization
      WRITE (6, '(2x,a,i4)') 'Error factoring matrix in blk3d: block = '
     1                        , k
      IF (ier < 0)WRITE (6,'(i4, a)') 
     1    ier, 'th argument has illegal value'
      IF (ier > 0)WRITE (6,'(i4, a)') 
     1    ier, 'th diagonal factor exactly zero'
      STOP
  305 CONTINUE
      WRITE (6, '(2/a,i4,2/)') ' BLK3D:   error detected:   ier =',
     1   ier
      STOP

  400 CONTINUE

      DEALLOCATE(amat, bmat, cmat, temp, stat=ier)

      CALL CloseDAFile

      END SUBROUTINE blk3d_factor_swp2


      SUBROUTINE blk3d_slv(ablk, qblk, bp1, source, 
     1                     ipiv, mblk, nblocks)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE stel_kinds
!      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(in) :: ipiv(mblk,nblocks)
      REAL(sp), TARGET, DIMENSION(mblk,mblk,nblocks), INTENT(in) :: 
     1                           ablk, qblk, bp1
      REAL(rprec2), DIMENSION(mblk,nblocks), INTENT(inout) 
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
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(in) :: ipiv(mblk,nblocks)
      REAL(rprec2), DIMENSION(mblk,nblocks), INTENT(inout) 
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
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(in) :: ipiv(mblk,nblocks)
      REAL(rprec2), DIMENSION(mblk,nblocks), INTENT(inout) 
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
