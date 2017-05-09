#if defined(SKS)
      SUBROUTINE scalfor_par(gcx, axm, bxm, axd, bxd, cx, iflag)
      USE vmec_main
      USE vmec_params
      USE vmec_dim, ONLY: ns 
      USE realspace, ONLY: wint, ru0
      USE parallel_include_module
      USE xstuff, ONLY: pxc, pgc
      IMPLICIT NONE
C-----------------------------------------------
C   Dummy Arguments
C-----------------------------------------------
      INTEGER, INTENT(in) :: iflag
      REAL(rprec), DIMENSION(0:ntor,0:mpol1,ns,ntmax),
     1  INTENT(inout) :: gcx
      REAL(rprec), DIMENSION(ns+1,2), INTENT(in) ::
     1  axm, bxm, axd, bxd
      REAL(rprec), DIMENSION(ns), INTENT(in) :: cx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: ftol_edge = 1.e-9_dp, c1p5 = 1.5_dp,
     1      fac = 0.25_dp, edge_pedestal = 0.05_dp
      INTEGER :: m , mp, n, js, jmax, jmin4(0:mnsize-1)
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: ax, bx, dx
      REAL(rprec) :: mult_fac
      INTEGER :: nsmin, nsmax, i, j, k, l, ier
      INTEGER :: blksize, numjs, left, right
      INTEGER, ALLOCATABLE, DIMENSION(:) :: counts, disps
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:,:) :: send_buf2
      INTEGER :: MPI_STAT(MPI_STATUS_SIZE)
      REAL(rprec) :: allgvton, allgvtoff
      REAL(rprec) :: tridslvton, tridslvtoff
      REAL(rprec) :: scalforton, scalfortoff
      REAL(rprec) :: skston, skstoff
C-----------------------------------------------
      IF (.NOT.lactive) return

      CALL second0(scalforton)
      left=rank-1; IF(rank.EQ.0) left=MPI_PROC_NULL
      right=rank+1; IF(rank.EQ.nranks-1) right=MPI_PROC_NULL

      CALL second0(skston)

      IF(grank.LT.nranks) THEN
        CALL MPI_Sendrecv(axm(tlglob,1),1,MPI_REAL8,left,1,
     1  axm(t1rglob,1),1,MPI_REAL8,right,1,NS_COMM,
     2  MPI_STAT, MPI_ERR)

        CALL MPI_Sendrecv(axm(tlglob,2),1,MPI_REAL8,left,1,
     1  axm(t1rglob,2),1,MPI_REAL8,right,1,NS_COMM,
     2  MPI_STAT, MPI_ERR)

        CALL MPI_Sendrecv(bxm(tlglob,1),1,MPI_REAL8,left,1,
     1  bxm(t1rglob,1),1,MPI_REAL8,right,1,NS_COMM,
     2  MPI_STAT, MPI_ERR)

        CALL MPI_Sendrecv(bxm(tlglob,2),1,MPI_REAL8,left,1,
     1  bxm(t1rglob,2),1,MPI_REAL8,right,1,NS_COMM,
     2  MPI_STAT, MPI_ERR)
      END IF

      CALL second0(skstoff)
      sendrecv_time = sendrecv_time + (skstoff - skston)

      ALLOCATE (ax(0:ntor,0:mpol1,ns), bx(0:ntor,0:mpol1,ns),
     1          dx(0:ntor,0:mpol1,ns))
      ax(:,:,1) = 0; bx(:,:,1) = 0; dx(:,:,1) = 0
      ax(:,:,ns) = 0; bx(:,:,ns) = 0; dx(:,:,ns) = 0

      jmax = ns
      IF (ivac .lt. 1) jmax = ns1

!     FOR SOME 3D PLASMAS, THIS SOMETIME HELPS (CHOOSE mult_fac =1 otherwise)
!     TO AVOID JACOBIAN RESETS BY GIVING A SMOOTH TRANSITION FROM FIXED TO FREE ITERATIONS
!      mult_fac = 1._dp/(1._dp + 10*(fsqr+fsqz))
!      gcx(ns,:,:,:) = mult_fac*gcx(ns,:,:,:)

      DO m = 0, mpol1
        mp = MOD(m,2) + 1
        DO n = 0, ntor
          nsmin=MAX(jmin2(m), tlglob); nsmax=MIN(trglob,jmax)
          DO js = nsmin, nsmax
            ax(n,m,js) = -(axm(js+1,mp) + bxm(js+1,mp)*m**2)
            bx(n,m,js) = -(axm(js,mp) + bxm(js,mp)*m**2)
            dx(n,m,js) = -(axd(js,mp) + bxd(js,mp)*m**2
     1                    + cx(js)*(n*nfp)**2)
          END DO

          IF (m .eq. 1) THEN
            dx(n,m,2) = dx(n,m,2) + bx(n,m,2)
          END IF
        END DO
      END DO

      IF (jmax .ge. ns) THEN
!
!     SMALL EDGE PEDESTAL NEEDED TO IMPROVE CONVERGENCE
!     IN PARTICULAR, NEEDED TO ACCOUNT FOR POTENTIAL ZERO
!     EIGENVALUE DUE TO NEUMANN (GRADIENT) CONDITION AT EDGE
!
         dx(:,0:1,ns)     = (1+edge_pedestal)  *dx(:,0:1,ns)
         dx(:,2:mpol1,ns) = (1+2*edge_pedestal)*dx(:,2:mpol1,ns)
!
!     STABILIZATION ALGORITHM FOR ZC_00(NS)
!     FOR UNSTABLE CASE, HAVE TO FLIP SIGN OF -FAC -> +FAC FOR CONVERGENCE
!     COEFFICIENT OF < Ru (R Pvac)> ~ -fac*(z-zeq) WHERE fac (EIGENVALUE, OR
!     FIELD INDEX) DEPENDS ON THE EQUILIBRIUM MAGNETIC FIELD AND CURRENT,
!     AND zeq IS THE EQUILIBRIUM EDGE VALUE OF Z00
          mult_fac = MIN(fac, fac*hs*15)
          IF (iflag .eq. 1) THEN
!
!     METHOD 1: SUBTRACT (INSTABILITY) Pedge ~ fac*z/hs FROM PRECONDITIONER AT EDGE
!
             dx(0,0,ns) = dx(0,0,ns)*(1-mult_fac)/(1+edge_pedestal)
          END IF

      ENDIF
!
!     ACCELERATE (IMPROVE) CONVERGENCE OF FREE BOUNDARY. THIS WAS ADDED
!     TO DEAL WITH CASES WHICH MIGHT OTHERWISE DIVERGE. BY DECREASING THE
!     FSQ TOLERANCE LEVEL WHERE THIS KICKS IN (FTOL_EDGE), THE USER CAN
!     TURN-OFF THIS FEATURE
!
!     DIAGONALIZE (DX DOMINANT) AND REDUCE FORCE (DX ENHANCED) AT EDGE 
!     TO IMPROVE CONVERGENCE FOR N != 0 TERMS
!

!      ledge = .false.       
!      IF ((fsqr+fsqz) .lt. ftol_edge) ledge = .true.
!      IF ((iter2-iter1).lt.400 .or. ivac.lt.1) ledge = .false.

!      IF (ledge) THEN
!         dx(ns,1:,1:) = 3*dx(ns,1:,1:)
!      END IF

!     FOR DATA MATCHING MODE (0 <= IRESIDUE < 3),
!     MAGNETIC AXIS IS FIXED SO JMIN3(0) => 2 FOR M=0,N=0

      jmin4 = jmin3
      IF (iresidue.GE.0 .and. iresidue.LT.3) jmin4(0) = 2

      IF (nranks.GT.1.AND.grank.LT.nranks) THEN
        IF (THOMAS.AND..NOT.BCYCLIC.AND.SKS_ALLGATHER) THEN
          CALL second0(allgvton)

          numjs=trglob-tlglob+1
          blksize=(ntor+1)*(mpol1+1)

          ALLOCATE (counts(nranks),disps(nranks))

          DO i=1,nranks
            counts(i)=(trglob_arr(i)-tlglob_arr(i)+1)*blksize
          END DO

          disps(1)=0
          DO i=2,nranks
            disps(i)=disps(i-1)+counts(i-1)
          END DO

          CALL MPI_Allgatherv(MPI_IN_PLACE,numjs*blksize,MPI_REAL8,ax,
     1                        counts,disps,MPI_REAL8,NS_COMM,MPI_ERR)
          CALL MPI_Allgatherv(MPI_IN_PLACE,numjs*blksize,MPI_REAL8,bx,
     1                        counts,disps,MPI_REAL8,NS_COMM,MPI_ERR)
          CALL MPI_Allgatherv(MPI_IN_PLACE,numjs*blksize,MPI_REAL8,dx,
     1                        counts,disps,MPI_REAL8,NS_COMM,MPI_ERR)
          DEALLOCATE(counts,disps)

          blksize=(ntor+1)*(mpol1+1)*ntmax

          ALLOCATE (counts(nranks),disps(nranks))

          DO i=1,nranks
            counts(i)=(trglob_arr(i)-tlglob_arr(i)+1)*blksize
          END DO

          disps(1)=0
          DO i=2,nranks
            disps(i)=disps(i-1)+counts(i-1)
          END DO

          CALL MPI_Allgatherv(MPI_IN_PLACE,numjs*blksize,MPI_REAL8,gcx,
     1                        counts,disps,MPI_REAL8,NS_COMM,MPI_ERR)
          DEALLOCATE(counts,disps)

          CALL second0(allgvtoff)
          allgather_time = allgather_time + (allgvtoff-allgvton)

          CALL second0(tridslvton)
          ier = 0
          CALL serial_tridslv_modified(ax,dx,bx,gcx,jmin4,jmax,mnsize-1,
     1                          ns,ntmax,ier)

          CALL MPI_Allreduce(MPI_IN_PLACE,ier,1,MPI_INTEGER,MPI_SUM,
     1                       NS_COMM,MPI_ERR)
          IF (ier /=0) THEN
             lerror_sam = .true.
             DEALLOCATE(ax,bx,dx)
             RETURN
          END IF
          CALL second0(tridslvtoff)
          tridslv_time = tridslv_time + (tridslvtoff-tridslvton)

        ELSE IF(BCYCLIC.AND..NOT.THOMAS) THEN

!      IF (iter2.GE.63) THEN
!        CALL PrintOutLinearArray(pxc,tlglob,trglob,.TRUE.,2000)
!        CALL PrintOutLinearArray(pgc,tlglob,trglob,.TRUE.,3000)
!        DO i=tlglob, trglob
!        DO j=0, ntor
!        DO k=0, mpol1
!          WRITE(40000+rank,"(3F20.8)") ax(j,k,i),bx(j,k,i),dx(j,k,i)
!        END DO
!        END DO
!        END DO
!      ELSE
!        CALL PrintOutLinearArray(pxc,tlglob,trglob,.TRUE.,200)
!        CALL PrintOutLinearArray(pgc,tlglob,trglob,.TRUE.,300)
!        DO i=tlglob, trglob
!        DO j=0, ntor
!        DO k=0, mpol1
!          WRITE(50000+rank,"(3F20.8)") ax(j,k,i),bx(j,k,i),dx(j,k,i)
!        END DO
!        END DO
!        END DO
!      END IF

          IF (grank.LT.nranks) THEN
            CALL second0(tridslvton)
            CALL bst_parallel_tridiag_solver(ax,dx,bx,gcx,jmin4,jmax,
     1              mnsize-1, ns,ntmax)
            CALL second0(tridslvtoff)
            tridslv_time = tridslv_time + (tridslvtoff-tridslvton)

            blksize=(ntor+1)*(mpol1+1)*ntmax
            CALL second0(skston)

            CALL MPI_Sendrecv(gcx(:,:,tlglob,:),blksize,MPI_REAL8,
     1  left,1,gcx(:,:,t1rglob,:),blksize,MPI_REAL8,right,1,NS_COMM,
     2  MPI_STAT, MPI_ERR)

            CALL MPI_Sendrecv(gcx(:,:,trglob,:),blksize,MPI_REAL8,
     1  right,1,gcx(:,:,t1lglob,:),blksize,MPI_REAL8,left,1,NS_COMM,
     2  MPI_STAT, MPI_ERR)
          END IF

!      IF (iter2.GE.63) THEN
!        CALL PrintOutLinearArray(pxc,tlglob,trglob,.TRUE.,8000)
!        CALL PrintOutLinearArray(pgc,tlglob,trglob,.TRUE.,9000)
!        CALL STOPMPI(666)
!      ELSE
!        CALL PrintOutLinearArray(pxc,tlglob,trglob,.TRUE.,800)
!        CALL PrintOutLinearArray(pgc,tlglob,trglob,.TRUE.,900)
!      END IF

          CALL second0(skstoff)
          sendrecv_time = sendrecv_time + (skstoff - skston)

        ELSE

          STOP 'Something wrong with choice of tridiagonal solver'

        END IF

      ELSE ! If only one processor
        CALL second0(tridslvton)
        ier = 0
        CALL serial_tridslv_modified(ax,dx,bx,gcx,jmin4,jmax,mnsize-1,
     1                              ns,ntmax,ier)
        IF (ier /=0) THEN
           lerror_sam = .true.
           DEALLOCATE(ax,bx,dx)
           RETURN
        END IF
        CALL second0(tridslvtoff)
        tridslv_time = tridslv_time + (tridslvtoff-tridslvton)
      END IF

      DEALLOCATE (ax, bx, dx)

      CALL second0(scalfortoff)
      scalfor_time =  scalfor_time + (scalfortoff-scalforton)

      END SUBROUTINE scalfor_par

      SUBROUTINE bst_parallel_tridiag_solver(a, d, b, c, jmin, 
     1           jmax, mnd1, ns, nrhs)
      USE stel_kinds
      USE parallel_include_module
      USE blocktridiagonalsolver_bst, ONLY: SetMatrixRowColL_bst
      USE blocktridiagonalsolver_bst, ONLY: SetMatrixRowColD_bst
      USE blocktridiagonalsolver_bst, ONLY: SetMatrixRowColU_bst
      USE blocktridiagonalsolver_bst, ONLY: ForwardSolve_bst
      USE blocktridiagonalsolver_bst, ONLY: SetMatrixRHS_bst
      USE blocktridiagonalsolver_bst, ONLY: BackwardSolve_bst
      USE blocktridiagonalsolver_bst, ONLY: GetSolutionVector_bst

      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: jmax, mnd1, ns, nrhs
      INTEGER, DIMENSION(0:mnd1), INTENT(in) :: jmin
      REAL(rprec), DIMENSION(0:mnd1,ns) :: a, d, b
      REAL(rprec), DIMENSION(0:mnd1,ns,nrhs), INTENT(inout) :: c
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mn, in, in1, jrhs
      INTEGER :: irow, icol, blklength, i, j
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: tmp
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: tmpv
      REAL(rprec), DIMENSION(0:mnd1) :: psi0
      REAL(rprec) :: t1, t2
C-----------------------------------------------
!     SOLVES B(I)*X(I-1)+D(I)*X(I)+A(I)*X(I+1)=C(I), I=IN,JMAX
!     AND RETURNS ANSWER IN C(I)

      CALL second0(t1)
      IF (jmax .gt. ns) STOP 'jmax>ns in tridslv_par'
      in = MINVAL(jmin)
      DO mn = 0, mnd1
         in1 = jmin(mn)-1
         IF (in1 .ge. in) THEN
            d(mn, in:in1) = 1
            b(mn, in:in1) = 0
            a(mn, in:in1) = 0
            c(mn, in:in1, 1:nrhs) = 0
         END IF
      END DO

      blklength=mnd1+1
      ALLOCATE(tmp(blklength,blklength))
      tmp=zero
      CALL second0(t2)
      init_time = init_time + (t2-t1)
      
      CALL second0(t1)
      DO irow = tlglob, trglob
        
        ! Set up L
        !IF (irow.EQ.ns) THEN
        IF(irow.EQ.ns.AND.jmax.LT.ns) THEN
          b(:,irow) = 0
        END IF
        CALL SetMatrixRowColL_bst(irow,b(:,irow))

        ! Set up D
        IF(irow.EQ.ns.AND.jmax.LT.ns) THEN
          d(:,irow) = 1
        END IF
        CALL SetMatrixRowColD_bst(irow,d(:,irow))

        ! Set up U
        CALL SetMatrixRowColU_bst(irow,a(:,irow))
      END DO
      CALL second0(t2)
      setup_time = setup_time + (t2-t1)

      CALL second0(t1)
      CALL ForwardSolve_bst
      CALL second0(t2)
      forwardsolve_time = forwardsolve_time + (t2-t1)

      ALLOCATE(tmpv(0:mnd1))
      CALL second0(t1)
      DO jrhs = 1, nrhs
        
        ! Set RHS
        DO irow = tlglob, trglob
          tmpv(0:mnd1)=c(:,irow,jrhs)
          IF (irow.EQ.ns.AND.jmax.LT.ns) tmpv(0:mnd1)=0
          CALL SetMatrixRHS_bst(irow,tmpv) 
        END DO

        ! Backward solve
        CALL BackwardSolve_bst

        ! Get solution vector
        DO irow = tlglob, trglob
          CALL GetSolutionVector_bst(irow, tmpv)
          c(:,irow,jrhs)=tmpv(0:mnd1)
        END DO

      END DO
      DEALLOCATE(tmp, tmpv)
      CALL second0(t2)
      backwardsolve_time = backwardsolve_time + (t2-t1)

      END SUBROUTINE bst_parallel_tridiag_solver

      SUBROUTINE serial_tridslv_modified(a, d, b, c, jmin, jmax, mnd1,
     1                                   ns,nrhs,ier)
      USE stel_kinds
      USE vmec_main, ONLY: lfreeb
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: jmax, mnd1, ns, nrhs
      INTEGER, DIMENSION(0:mnd1), INTENT(in) :: jmin
      REAL(rprec), DIMENSION(0:mnd1,ns) :: a, d, b
      REAL(rprec), DIMENSION(0:mnd1,ns,nrhs), INTENT(inout) :: c
      INTEGER, INTENT(inout) :: ier
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mn, in, i0, in1, jrhs
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: alf
      REAL(rprec), DIMENSION(0:mnd1) :: psi0
C-----------------------------------------------
!
!     SOLVES B(I)*X(I-1)+D(I)*X(I)+A(I)*X(I+1)=C(I), I=IN,JMAX
!     AND RETURNS ANSWER IN C(I)
!     ADDED VECTORIZATION ON FOURIER MODE ARGUMENT (01-2000)
!     AND NEW ARGUMENT (NRHS) TO DO MULTIPLE RIGHT SIDES SIMULTANEOUSLY
!
      ier = 0
      IF (jmax .gt. ns) STOP 'jmax>ns in tridslv_par'

      ALLOCATE (alf(0:mnd1,ns), stat=in)
      IF (in .ne. 0) STOP 'Allocation error in tridslv_par'

      in = MINVAL(jmin)
!
!      FILL IN MN BELOW MAX(JMIN) WITH DUMMY VALUES
!      TO ALLOW VECTORIZATION ON MN INDEX
!
      DO mn = 0, mnd1
         in1 = jmin(mn)-1
         IF (in1 .ge. in) THEN
            d(mn, in:in1) = 1
            c(mn, in:in1, 1:nrhs) = 0
            b(mn, in:in1) = 0
            a(mn, in:in1) = 0
         END IF
      END DO

      in1 = in + 1

      psi0(:)= d(:,in)
      IF (ANY(psi0 .eq. zero)) STOP 'psi0 == 0 error in tridslv_par'
      psi0 = one/psi0
      DO jrhs = 1, nrhs
         c(:,in,jrhs) = c(:,in,jrhs)*psi0(:)
      END DO

      DO i0 = in1,jmax
         alf(:,i0-1) = a(:,i0-1)*psi0(:)
         psi0 = d(:,i0) - b(:,i0)*alf(:,i0-1)
!         IF (ANY(ABS(psi0) .le. 1.E-8_dp*ABS(d(:,i0)))) 
!     1       STOP 'psi0/d(i0) < 1.E-8: possible singularity in tridslv'
         IF (ANY(ABS(psi0) .le. 1.E-8_dp*ABS(d(:,i0)))) THEN
            DEALLOCATE(alf)
            ier = -1
            RETURN
         END IF
         psi0  = one/psi0
         DO jrhs = 1, nrhs
            c(:,i0,jrhs) = (c(:,i0,jrhs) - b(:,i0)*c(:,i0-1,jrhs))
     1                   * psi0
         END DO
      END DO

      DO i0 = jmax - 1, in, -1
         DO jrhs = 1,nrhs
            c(:,i0,jrhs) = c(:,i0,jrhs) - alf(:,i0)*c(:,i0+1,jrhs)
         END DO
      END DO

      DEALLOCATE (alf)

      END SUBROUTINE serial_tridslv_modified

#endif      

      SUBROUTINE scalfor(gcx, axm, bxm, axd, bxd, cx, iflag)
      USE vmec_main
      USE vmec_params
      USE vmec_dim, ONLY: ns 
      USE realspace, ONLY: wint, ru0
#if defined(SKS)      
      USE xstuff, ONLY: xc, gc
      USE parallel_include_module
#endif      
      IMPLICIT NONE
C-----------------------------------------------
C   Dummy Arguments
C-----------------------------------------------
      INTEGER, INTENT(in) :: iflag
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax),
     1  INTENT(inout) :: gcx
      REAL(rprec), DIMENSION(ns+1,2), INTENT(in) ::
     1  axm, bxm, axd, bxd
      REAL(rprec), DIMENSION(ns), INTENT(in) :: cx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: ftol_edge = 1.e-9_dp, c1p5 = 1.5_dp,
     1      fac = 0.25_dp, edge_pedestal = 0.05_dp
      INTEGER :: m , mp, n, js, jmax, jmin4(0:mnsize-1)
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: ax, bx, dx
      REAL(rprec) :: mult_fac
!      LOGICAL :: ledge
#if defined(SKS)      
      INTEGER :: nsmin, nsmax, i, j, k, l
      REAL(rprec) :: tridslvton, tridslvtoff
      REAL(rprec) :: scalforton, scalfortoff
#endif
C-----------------------------------------------
#if defined(SKS)      
      CALL second0(scalforton)
#endif
      ALLOCATE (ax(ns,0:ntor,0:mpol1), bx(ns,0:ntor,0:mpol1),
     1          dx(ns,0:ntor,0:mpol1))
      ax(1,:,:) = 0; bx(1,:,:) = 0; dx(1,:,:) = 0
      ax(ns,:,:) = 0; bx(ns,:,:) = 0; dx(ns,:,:) = 0

      jmax = ns
      IF (ivac .lt. 1) jmax = ns1

!     FOR SOME 3D PLASMAS, THIS SOMETIME HELPS (CHOOSE mult_fac =1 otherwise)
!     TO AVOID JACOBIAN RESETS BY GIVING A SMOOTH TRANSITION FROM FIXED TO FREE ITERATIONS
!      mult_fac = 1._dp/(1._dp + 10*(fsqr+fsqz))
!      gcx(ns,:,:,:) = mult_fac*gcx(ns,:,:,:)

      DO m = 0, mpol1
         mp = MOD(m,2) + 1
         DO n = 0, ntor
            DO js = jmin2(m), jmax
               ax(js,n,m) = -(axm(js+1,mp) + bxm(js+1,mp)*m**2)
               bx(js,n,m) = -(axm(js,mp) + bxm(js,mp)*m**2)
               dx(js,n,m) = -(axd(js,mp) + bxd(js,mp)*m**2
     1                    + cx(js)*(n*nfp)**2)
            END DO

            IF (m .eq. 1) THEN
               dx(2,n,m) = dx(2,n,m) + bx(2,n,m)
!OFF 050311        DO js = jmin2(m), jmax
!OFF 050311        ax(js,n,m) = c1p5*ax(js,n,m)
!OFF 050311        bx(js,n,m) = c1p5*bx(js,n,m)
!OFF 050311        dx(js,n,m) = c1p5*dx(js,n,m)
!OFF 050311        END DO
            END IF
         END DO
      END DO

      IF (jmax .ge. ns) THEN
!
!     SMALL EDGE PEDESTAL NEEDED TO IMPROVE CONVERGENCE
!     IN PARTICULAR, NEEDED TO ACCOUNT FOR POTENTIAL ZERO
!     EIGENVALUE DUE TO NEUMANN (GRADIENT) CONDITION AT EDGE
!
         dx(ns,:,0:1)     = (1+edge_pedestal)  *dx(ns,:,0:1)
         dx(ns,:,2:mpol1) = (1+2*edge_pedestal)*dx(ns,:,2:mpol1)
!
!     STABILIZATION ALGORITHM FOR ZC_00(NS)
!     FOR UNSTABLE CASE, HAVE TO FLIP SIGN OF -FAC -> +FAC FOR CONVERGENCE
!     COEFFICIENT OF < Ru (R Pvac)> ~ -fac*(z-zeq) WHERE fac (EIGENVALUE, OR
!     FIELD INDEX) DEPENDS ON THE EQUILIBRIUM MAGNETIC FIELD AND CURRENT,
!     AND zeq IS THE EQUILIBRIUM EDGE VALUE OF Z00
          mult_fac = MIN(fac, fac*hs*15)
          IF (iflag .eq. 1) THEN
!
!     METHOD 1: SUBTRACT (INSTABILITY) Pedge ~ fac*z/hs FROM PRECONDITIONER AT EDGE
!
             dx(ns,0,0) = dx(ns,0,0)*(1-mult_fac)/(1+edge_pedestal)
          END IF
      ENDIF


!
!     ACCELERATE (IMPROVE) CONVERGENCE OF FREE BOUNDARY. THIS WAS ADDED
!     TO DEAL WITH CASES WHICH MIGHT OTHERWISE DIVERGE. BY DECREASING THE
!     FSQ TOLERANCE LEVEL WHERE THIS KICKS IN (FTOL_EDGE), THE USER CAN
!     TURN-OFF THIS FEATURE
!
!     DIAGONALIZE (DX DOMINANT) AND REDUCE FORCE (DX ENHANCED) AT EDGE 
!     TO IMPROVE CONVERGENCE FOR N != 0 TERMS
!

!      ledge = .false.       
!      IF ((fsqr+fsqz) .lt. ftol_edge) ledge = .true.
!      IF ((iter2-iter1).lt.400 .or. ivac.lt.1) ledge = .false.

!      IF (ledge) THEN
!         dx(ns,1:,1:) = 3*dx(ns,1:,1:)
!      END IF

!     FOR DATA MATCHING MODE (0 <= IRESIDUE < 3),
!     MAGNETIC AXIS IS FIXED SO JMIN3(0) => 2 FOR M=0,N=0

      jmin4 = jmin3
      IF (iresidue.GE.0 .and. iresidue.LT.3) jmin4(0) = 2

!     SOLVES BX(I)*X(I-1)+DX(I)*X(I)+AX(I)*X(I+1)=GCX(I), I=JMIN4,JMAX
!     AND RETURNS ANSWER IN GCX(I)

!      IF (iter2.GE.63) THEN
!        CALL PrintOutLinearArray(xc,tlglob,trglob,.FALSE.,2000)
!        CALL PrintOutLinearArray(gc,tlglob,trglob,.FALSE.,3000)
!        DO i=tlglob, trglob
!        DO j=0, ntor
!        DO k=0, mpol1
!          WRITE(40000+rank,"(3F20.8)") ax(i,j,k),bx(i,j,k),dx(i,j,k)
!        END DO
!        END DO
!        END DO
!      ELSE
!        CALL PrintOutLinearArray(xc,tlglob,trglob,.FALSE.,200)
!        CALL PrintOutLinearArray(gc,tlglob,trglob,.FALSE.,300)
!        DO i=tlglob, trglob
!        DO j=0, ntor
!        DO k=0, mpol1
!          WRITE(50000+rank,"(3F20.8)") ax(i,j,k),bx(i,j,k),dx(i,j,k)
!        END DO
!        END DO
!        END DO
!      END IF

#if defined(SKS)      
      CALL second0(tridslvton)
#endif
      CALL serial_tridslv (ax, dx, bx, gcx, jmin4, jmax, mnsize-1, 
     1                     ns, ntmax)
#if defined(SKS)      
      CALL second0(tridslvtoff)
      s_tridslv_time = s_tridslv_time + (tridslvtoff-tridslvton)
#endif
      DEALLOCATE (ax, bx, dx)
#if defined(SKS)      
      CALL second0(scalfortoff)
      s_scalfor_time =  s_scalfor_time + (scalfortoff-scalforton)
#endif
!      CALL STOPMPI(999)

!      IF (iter2.GE.63) THEN
!        CALL PrintOutLinearArray(xc,tlglob,trglob,.FALSE.,8000)
!        CALL PrintOutLinearArray(gc,tlglob,trglob,.FALSE.,9000)
!        CALL STOPMPI(666)
!      ELSE
!        CALL PrintOutLinearArray(xc,tlglob,trglob,.FALSE.,800)
!        CALL PrintOutLinearArray(gc,tlglob,trglob,.FALSE.,900)
!      END IF

      END SUBROUTINE scalfor

      SUBROUTINE serial_tridslv(a, d, b, c, jmin, jmax, mnd1, ns, nrhs)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: jmax, mnd1, ns, nrhs
      INTEGER, DIMENSION(0:mnd1), INTENT(in) :: jmin
      REAL(rprec), DIMENSION(ns,0:mnd1) :: a, d, b
      REAL(rprec), DIMENSION(ns,0:mnd1, nrhs), INTENT(inout) :: c
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mn, in, i0, in1, jrhs
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: alf
      REAL(rprec), DIMENSION(0:mnd1) :: psi0
C-----------------------------------------------
!
!     SOLVES B(I)*X(I-1)+D(I)*X(I)+A(I)*X(I+1)=C(I), I=IN,JMAX
!     AND RETURNS ANSWER IN C(I)
!     ADDED VECTORIZATION ON FOURIER MODE ARGUMENT (01-2000)
!     AND NEW ARGUMENT (NRHS) TO DO MULTIPLE RIGHT SIDES SIMULTANEOUSLY
!
      IF (jmax .gt. ns) STOP 'jmax>ns in tridslv'

      ALLOCATE (alf(ns,0:mnd1), stat = in)
      IF (in .ne. 0) STOP 'Allocation error in tridslv'

      in = MINVAL(jmin)
!
!      FILL IN MN BELOW MAX(JMIN) WITH DUMMY VALUES
!      TO ALLOW VECTORIZATION ON MN INDEX
!
      DO mn = 0, mnd1
         in1 = jmin(mn)-1
         IF (in1 .ge. in) THEN
            d(in:in1, mn) = 1
            c(in:in1, mn, 1:nrhs) = 0
            b(in:in1, mn) = 0
            a(in:in1, mn) = 0
         END IF
      END DO

      in1 = in + 1

      psi0(:)= d(in,:)
      IF (ANY(psi0 .eq. zero)) STOP 'psi0 == 0 error in tridslv'
      psi0 = one/psi0
      DO jrhs = 1, nrhs
         c(in,:,jrhs) = c(in,:,jrhs)*psi0(:)
      END DO

      DO i0 = in1,jmax
         alf(i0-1,:) = a(i0-1,:)*psi0(:)
         psi0 = d(i0,:) - b(i0,:)*alf(i0-1,:)
         IF (ANY(ABS(psi0) .le. 1.E-8_dp*ABS(d(i0,:)))) 
     1       STOP 'psi0/d(i0) < 1.E-8: possible singularity in tridslv'
         psi0  = one/psi0
         DO jrhs = 1, nrhs
            c(i0,:,jrhs) = (c(i0,:,jrhs) - b(i0,:)*c(i0-1,:,jrhs))
     1                   * psi0
         END DO
      END DO

      DO i0 = jmax - 1, in, -1
         DO jrhs = 1,nrhs
            c(i0,:,jrhs) = c(i0,:,jrhs) - alf(i0,:)*c(i0+1,:,jrhs)
         END DO
      END DO

      DEALLOCATE (alf)

      END SUBROUTINE serial_tridslv

