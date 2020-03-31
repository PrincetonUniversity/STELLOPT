      MODULE gmres_mod
      USE vmec_main, ONLY: dp, rprec, neqs, ns, nthreed, 
     1                     one, fsqr, fsqz, fsql
      USE parallel_include_module
      USE parallel_vmec_module, ONLY: CopyLastNtype, SaxpbyLastNtype,
     1                                SaxpbyLastNs, Saxpby1LastNs,
     2                                GetDerivLastNs, SaxLastNtype
      USE precon2d, ONLY: ictrl_prec2d
      IMPLICIT NONE
      INTEGER :: nfcn = 0, lqmrstep = 0
      INTEGER :: ier_flag_res
      LOGICAL :: lqmr, lfirst
      LOGICAL, PARAMETER :: lscreen0 = .FALSE.

!
!     nfcn :  number of calls to function (funct3d)
!     lqmr :  logical, used by external programs to control calling these routines
!
      CONTAINS

      SUBROUTINE matvec_par (ploc, Ap, nloc)
      USE blocktridiagonalsolver, ONLY: ParMatVec
      USE stel_kinds
      USE xstuff, ONLY: pxc, px0=>pxsave, pgc0=>pxcdot, pgc
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN)   :: nloc
      REAL(dp), INTENT(IN)  :: 
     &   ploc(ntmaxblocksize,tlglob:trglob)
      REAL(dp), INTENT(OUT) :: 
     &   Ap(ntmaxblocksize,tlglob:trglob)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER  :: mblk_size, istat
      REAL(dp) :: delta, lmax, gmax, pmax, apmax
C-----------------------------------------------
      mblk_size = neqs/ns
!
!     Computes linearized matrix product A*p = [F(x0+delta*p) - F(x0)]/delta, about point x0
!     Must scale p so delta*max|p| ~ sqrt(epsilon) to get accurate numerical derivative
!     Because the block preconditioner is applied in funct3d, the net result
!     should be Ap ~ -p.

      delta = SQRT(EPSILON(delta))

      LACTIVE0: IF (lactive) THEN
         IF (nloc .NE. (trglob - tlglob + 1)*mblk_size) THEN
            STOP 'nloc wrong in matvec_par'
         END IF
         lmax = SUM(ploc*ploc)
         CALL MPI_ALLREDUCE(lmax, pmax, 1, MPI_REAL8, MPI_SUM, NS_COMM,
     &                      MPI_ERR)
         pmax = SQRT(pmax)
         delta = delta/MAX(delta, pmax)

         CALL SaxpbyLastNs(delta, ploc, one, px0, pxc)

         CALL last_ntype_par
         CALL PadSides(pxc)
      END IF LACTIVE0

      CALL funct3d_par(lscreen0, ier_flag_res)
      
      IF (lactive) THEN
         CALL last_ns_par
         CALL GetDerivLastNs(pgc, pgc0, delta, Ap)
      ENDIF

      IF (ier_flag_res.NE.0 .AND. rank.EQ.0) THEN
         PRINT *,'IN MATVEC_PAR, IER_FLAG = ', ier_flag_res
      END IF

 90   CONTINUE
      nfcn = nfcn + 1

      END SUBROUTINE matvec_par

      SUBROUTINE GetNLForce_par(xcstate, fsq_nl, bnorm)
      USE xstuff, ONLY: pxc, pgc, x0=>pxsave
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN)  :: xcstate(neqs), bnorm
      REAL(dp), INTENT(OUT) :: fsq_nl
!-----------------------------------------------
!undo internal gmres normalization

      LACTIVE0: IF (lactive) THEN
         CALL Saxpby1LastNs(bnorm, xcstate, one, x0, pxc)
         CALL last_ntype_par
         CALL PadSides(pxc)
      END IF LACTIVE0

      CALL funct3d_par(lscreen0, ier_flag_res)
      IF (lactive) CALL last_ns_par

      fsq_nl = fsqr + fsqz + fsql
      nfcn = nfcn + 1

      END SUBROUTINE GetNLForce_par
  !------------------------------------------------
  ! 
  !------------------------------------------------
      SUBROUTINE last_ns_par 
         USE xstuff

         REAL(dp), ALLOCATABLE, DIMENSION(:)  :: tmp
         ALLOCATE (tmp(ntmaxblocksize*ns))

         CALL tolastns(pgc,tmp)
         CALL copylastns(tmp,pgc)

         CALL tolastns(pxcdot,tmp)
         CALL copylastns(tmp,pxcdot)

         CALL tolastns(pxc,tmp)
         CALL copylastns(tmp,pxc)

         CALL tolastns(pxsave,tmp)
         CALL copylastns(tmp,pxsave)

         DEALLOCATE(tmp)

      END SUBROUTINE last_ns_par
  !------------------------------------------------

  !------------------------------------------------
      SUBROUTINE last_ntype_par 
         USE xstuff

         REAL(dp), ALLOCATABLE, DIMENSION(:)  :: tmp
         ALLOCATE (tmp(ntmaxblocksize*ns))

         CALL tolastntype(pgc,tmp)
         CALL copylastntype(tmp,pgc)

         CALL tolastntype(pxcdot,tmp)
         CALL copylastntype(tmp,pxcdot)

         CALL tolastntype(pxc,tmp)
         CALL copylastntype(tmp,pxc)

         CALL tolastntype(pxsave,tmp)
         CALL copylastntype(tmp,pxsave)

         DEALLOCATE(tmp)
      END SUBROUTINE last_ntype_par
  !------------------------------------------------

  !------------------------------------------------
      SUBROUTINE gmres_fun_par (ier_flag, itype, lscreen)
      USE precon2d, ONLY: ictrl_prec2d, block_precond_par
      USE xstuff
      USE vmec_main, ONLY: fsqr, fsqz, fsql, ftolv
      USE gmres_lib, ONLY: gmres_par, gmres_info
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: itype
      INTEGER, INTENT(OUT) :: ier_flag
      LOGICAL, INTENT(IN)  :: lscreen
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      TYPE (gmres_info)    :: gi
      INTEGER              :: n, m, l
      INTEGER              :: icntl(9), info(3)
      REAL(dp)             :: cntl(5), fact, fact_min, fsq_min, fsq2,
     &                        fsqr_min, fsqz_min, fsql_min
      CHARACTER(LEN=*), PARAMETER :: qmr_message = 
     &                              'Beginning GMRES iterations'
      INTEGER, PARAMETER   :: noPrec=0, leftPrec=1, rightPrec=2
!-----------------------------------------------
      LGMRESCALL = .TRUE.
      n = neqs
!
!     CHOOSE TYPE OF SOLVER
!
      IF (itype == 2) THEN
         CALL gmresr_fun (ier_flag)
         RETURN
      ELSE IF (itype == 3) THEN
         CALL qmr_fun
         RETURN
      END IF


      IF (lfirst) THEN
         lfirst = .FALSE.
         IF (grank.EQ.0) THEN
            WRITE (*,'(2x,a,/)') qmr_message
            WRITE (nthreed, '(2x,a,/)') qmr_message
         END IF
      END IF

!******************************************************
!* Initialize the control parameters to default value
!******************************************************

      CALL init_dgmres(icntl,cntl)

!************************
! Tune some parameters
!************************

! Tolerance
      cntl(1) = 1.e-3_dp
!      cntl(1) = 1.e-4_dp
! Write errors to fort.21
!      icntl(1) = 21   !21
! Write warnings to fort.21
      icntl(2) = 0   !21
! Save the convergence history in file fort.20
      icntl(3) = 20
! No preconditioning
      icntl(4) = noPrec
! Left preconditioning (doesn't work well with no col-scaling)     
!      icntl(4) = leftPrec 
! ICGS orthogonalization
      icntl(5) = 3
!      icntl(5) = 0
! Initial guess
      icntl(6) = 0
!      icntl(6) = 1
! Maximum number of iterations at each step (~ns/5)
!      icntl(7) = 15
      icntl(7) = 20
! Stops to peek at progress during rev com loop
      icntl(9) = 1             


!********************************
!* Choose the restart parameter
!********************************
!      write(*,*) 'Restart  <', ldstrt
!      read(*,*) m
!
!     m <= n
!
      m = 20
!
!      Load gmres_info structure
!
      info = 0
      gi%m=m; gi%icntl=icntl; gi%cntl=cntl; gi%info = info
      gi%ftol = ftolv


      ALLOCATE(gi%rcounts(nranks),gi%disp(nranks))
      gi%startglobrow=tlglob
      gi%endglobrow=trglob
      gi%iam=rank
      gi%nprocs=nranks
      gi%rcounts=ntblkrcounts
      gi%disp=ntblkdisp
      gi%mblk_size=ntmaxblocksize

      gi%my_comm = NS_COMM
      gi%my_comm_world = RUNVMEC_COMM_WORLD
      gi%lactive = lactive

      gi%lverbose = lscreen

      l = ictrl_prec2d
      IF (icntl(4) .NE. noPrec) ictrl_prec2d = -1
      CALL funct3d_par(lscreen0, ier_flag_res)
      ictrl_prec2d = l
      nfcn = nfcn+1

!STORE INITIAL POINT AND INITIAL FORCE (AT INIT PT)
!SO DEVIATIONS FROM THEM CAN BE COMPUTED REPEATEDLY
      CALL CopyLastNtype(pxc, pxsave) 
      CALL CopyLastNtype(pgc, pxcdot)

!RHS: RETURN RESULT OF SOLVING LINEARIZED A*x = -gc IN XCDOT 
!     AND IS DISTRIBUTED OVER ALL PROCESSORS
      CALL CopyLastNtype(pgc, pgc, -one)

      CALL last_ns_par

      CALL gmres_par(n, gi, matvec_par, block_precond_par,
     &               getnlforce_par, pxcdot, pgc)

      CALL last_ntype_par

!      ier_flag = gi%info(1)
!      ictrl_prec2d = 1
      ier_flag = 0
      fact = 1;  fact_min = fact
      fsq_min  = gi%ftol
      fsqr_min = fsqr; fsqz_min = fsqz; fsql_min = fsql

!     SIMPLE LINESEARCH SCALING SCAN

 1010 FORMAT(1x,'LINE SEARCH - SCAN ||X|| FOR MIN FSQ_NL',/,
     &       '-------------',/,
     &       1x,'ITER',7x,'FSQ_NL',10x,'||X||',9x,'MAX|X|')
     
      CALL MPI_BCAST(fsq_min, 1, MPI_REAL8, 0, RUNVMEC_COMM_WORLD, 
     &               MPI_ERR)

      DO m = 1, 5
         fact = fact*SQRT(0.5_dp)
         CALL SaxpbyLastNtype(fact, pxcdot, one, pxsave, pxc)
 
         CALL funct3d_par(lscreen0, ier_flag_res)

         fsq2 = fsqr+fsqz+fsql

         IF (fsq2 .LT. fsq_min) THEN
            fsq_min = fsq2
            fact_min = fact
            fsqr_min = fsqr; fsqz_min = fsqz; fsql_min = fsql
            IF (grank .EQ. 0) PRINT 1020, fact, fsq2
         ELSE
            EXIT
         END IF
      END DO

 1020 FORMAT(2x,'GMRES_FUN, TIME_STEP: ',1p,e10.3, ' FSQ_MIN: ',1pe10.3)

      fsqr = fsqr_min; fsqz = fsqz_min; fsql = fsql_min
      IF (ictrl_prec2d .EQ. 1) THEN
         CALL SaxLastNtype(pxcdot,pcol_scale,pxcdot)
      END IF

      CALL SaxpbyLastNtype(fact_min, pxcdot, one, pxsave, pxc)
      CALL COPYLASTNTYPE(pxc, pxsave)

      DEALLOCATE(gi%rcounts,gi%disp)
      LGMRESCALL=.FALSE.

      END SUBROUTINE gmres_fun_par

      SUBROUTINE GetNLForce(xcstate, fsq_nl, bnorm)
      USE xstuff, ONLY: xc, gc, x0=>xsave
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp),INTENT(IN)  :: xcstate(neqs), bnorm
      REAL(dp),INTENT(OUT) :: fsq_nl
!-----------------------------------------------
!undo internal gmres normalization
      xc(1:neqs) = x0(1:neqs) + bnorm*xcstate(1:neqs)

      CALL funct3d(lscreen0, ier_flag_res)
      fsq_nl = fsqr+fsqz+fsql
 
      nfcn = nfcn + 1

      END SUBROUTINE GetNLForce

      SUBROUTINE matvec (p, Ap, ndim)
      USE xstuff, ONLY: xc, x0=>xsave, gc0=>xcdot, gc
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in)                    :: ndim
      REAL(dp), INTENT(in), DIMENSION(ndim)  :: p
      REAL(dp), INTENT(out), DIMENSION(ndim) :: Ap
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER  :: l
      REAL(dp) :: delta, pmax
C-----------------------------------------------
!
!     Computes linearized matrix product A*p = [F(x0+delta*p) - F(x0)]/delta, about point x0
!     Must scale p so delta*max|p| ~ sqrt(epsilon) to get accurate numerical derivative
!
!      Ap = -p
!      GOTO 90

      delta = SQRT(EPSILON(delta))
      pmax = SQRT(SUM(p(1:ndim)**2))
      delta = delta/MAX(delta, pmax)

      xc(1:ndim) = x0(1:ndim) + delta*p(1:ndim)
      CALL funct3d(lscreen0, ier_flag_res)
      Ap = (gc(1:ndim) - gc0(1:ndim))/delta

      IF (ier_flag_res .NE. 0) THEN
         PRINT *,'IN MATVEC, IER_FLAG = ', ier_flag_res
      END IF

 90   CONTINUE
      nfcn = nfcn + 1

      END SUBROUTINE matvec

      SUBROUTINE gmres_fun (ier_flag, itype)
      USE precon2d, ONLY: block_precond
      USE xstuff
      USE vmec_main, ONLY: fsqr, fsqz, fsql, ftolv
      USE gmres_lib, ONLY: gmres_ser, gmres_info
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: itype
      INTEGER, INTENT(OUT) :: ier_flag
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      TYPE (gmres_info)   :: gi
      INTEGER :: n, m, l
      INTEGER :: icntl(9), info(3)
      REAL(dp) :: cntl(5), fact, fact_min, fsq_min, fsq2,
     &            fsqr_min, fsqz_min, fsql_min
      CHARACTER(LEN=*), PARAMETER :: qmr_message = 
     &                              'Beginning GMRES iterations'
      INTEGER, PARAMETER   :: noPrec=0, leftPrec=1, rightPrec=2
!-----------------------------------------------
      n = neqs

!
!     CHOOSE TYPE OF SOLVER
!
      IF (itype == 2) THEN
         CALL gmresr_fun (ier_flag)
         RETURN
      ELSE IF (itype == 3) THEN
         CALL qmr_fun
         RETURN
      END IF


      IF (lfirst) THEN
         lfirst = .FALSE.
         WRITE (*,'(2x,a,/)') qmr_message
         WRITE (nthreed, '(2x,a,/)') qmr_message
      END IF

!******************************************************
!* Initialize the control parameters to default value
!******************************************************

      CALL init_dgmres(icntl,cntl)

!************************
! Tune some parameters
!************************

! Tolerance
      cntl(1) = 1.e-3_dp
!      cntl(1) = 1.e-4_dp
! Write errors to fort.21
!      icntl(1) = 21   !21
! Write warnings to fort.21
      icntl(2) = 0   !21
! Save the convergence history in file fort.20
      icntl(3) = 20
! No preconditioning
      icntl(4) = noPrec
! Left preconditioning (doesn't work well with no col-scaling)     
!      icntl(4) = leftPrec 
! ICGS orthogonalization
      icntl(5) = 3
!      icntl(5) = 0
! Initial guess
      icntl(6) = 0
!      icntl(6) = 1
! Maximum number of iterations at each step (~ns/5)
!      icntl(7) = 15
      icntl(7) = 20
! Stops to peek at progress during rev com loop
      icntl(9) = 1             


!********************************
!* Choose the restart parameter
!********************************
!      write(*,*) 'Restart  <', ldstrt
!      read(*,*) m
!
!     m <= n
!
      m = 20
!
!      Load gmres_info structure
!
      info = 0
      gi%m=m; gi%icntl=icntl; gi%cntl=cntl; gi%info = info
      gi%ftol = ftolv

!Store initial gc
      l = ictrl_prec2d
      IF (icntl(4) .NE. noPrec) ictrl_prec2d = -1
      CALL funct3d(lscreen0, ier_flag_res)
      ictrl_prec2d = l
      nfcn = nfcn+1

!
!STORE INITIAL POINT AND INITIAL FORCE (AT INIT PT)
!SO DEVIATIONS FROM THEM CAN BE COMPUTED REPEATEDLY
      xcdot = gc
      xsave = xc

!RHS: RETURN RESULT OF A*X = -gc IN XCDOT
      gc = -gc

      CALL gmres_ser(n, gi, matvec, block_precond, getnlforce,
     &               xcdot, gc)

100   CONTINUE

!      ier_flag = gi%info(1)
!      ictrl_prec2d = 1
      ier_flag = 0
      fact = 1
      fact_min = fact
      fsq_min  = HUGE(fsq_min)

!     SIMPLE LINESEARCH SCALING SCAN
      DO m = 1, 5

         xc(1:n) = xsave(1:n) + fact*xcdot(1:n)
        
         CALL funct3d(lscreen0, ier_flag_res)

         fsq2 = fsqr+fsqz+fsql

         IF (fsq2 .LT. fsq_min) THEN
            fsq_min = fsq2
            fact_min = fact
            fsqr_min = fsqr; fsqz_min = fsqz; fsql_min = fsql
         ELSE
            EXIT
         END IF
         fact = fact/2._dp

      END DO

!      PRINT 1010, fact_min, fsq_min
! 1010 FORMAT(2x,'GMRES_FUN, TIME_STEP: ',1p,e10.3, ' FSQ_MIN: ',1pe10.3)

      xc(1:n) = xsave(1:n) + fact_min*xcdot(1:n)
      fsqr = fsqr_min; fsqz = fsqz_min; fsql = fsql_min
      xsave = xc

      END SUBROUTINE gmres_fun

      SUBROUTINE gmresr_fun (ier_flag)
      USE xstuff
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ier_flag
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ndim, jtrunc, mgmres, maxits, lwork
      LOGICAL :: oktest
      REAL(dp) :: eps, resid
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: work, delx, brhs
      CHARACTER(len=3), PARAMETER :: stc="rel"
      CHARACTER(LEN=*), PARAMETER :: qmr_message = 
     &                              'Beginning GMRESR iterations'
C-----------------------------------------------
      IF (lfirst) THEN
         lfirst = .false.
         WRITE (*,'(2x,a,/)') qmr_message
         WRITE (nthreed, '(2x,a,/)') qmr_message
      END IF

      oktest = .false.
      ndim   = neqs
      jtrunc = 10
      mgmres = 20
      maxits = 10
      lwork  = ndim*(2*jtrunc + mgmres + 2)
      eps = .3_dp

      ALLOCATE(work(lwork), delx(ndim), brhs(ndim), stat=ier_flag_res)
      IF (ier_flag_res .ne. 0) STOP 'Allocation failed in gmresr'

      brhs = -gc(1:ndim)
      delx = 0

      CALL gmresr(oktest, ndim, jtrunc, mgmres, brhs, delx, work,
     &            eps, stc, maxits, resid, matvec, ier_flag_res)

      xc(1:ndim) = xsave(1:ndim) + delx(1:ndim)

      DEALLOCATE (work, delx, brhs)

      ier_flag = 0

      END SUBROUTINE gmresr_fun

      SUBROUTINE qmr_fun
      USE vmec_dim, ONLY: ns, mpol1, ntor1
      USE vmec_params, ONLY: ntmax
      USE vmec_main, ONLY: lfreeb
      USE xstuff
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ndim, nlen, nlim, ierr, info(4)
      INTEGER :: revcom, colx, colb
      INTEGER :: nty, ntyp, mt, mp, nt, np, jp, js
      REAL(dp), DIMENSION(neqs,9) :: vecs
      REAL(dp) :: tol = 1.E-3_dp
      CHARACTER(LEN=*), PARAMETER :: qmr_message = 
     1                               'Beginning TF-QMR iterations'
      LOGICAL, PARAMETER :: ldump_fort33 = .false.
C-----------------------------------------------
      ndim = SIZE(vecs,1)
      nlen = ndim
      nlim = 10
      xcdot = gc
      xsave = xc

      IF (lfirst) THEN
         lfirst = .false.
         WRITE (*,'(2x,a,/)') qmr_message
         WRITE (nthreed, '(2x,a,/)') qmr_message
         IF (ldump_fort33) THEN

!     CHECK THAT INITIALLY, PRECONDITIONER YIELDS (APPROXIMATE) IDENTITY MATRIX
!     OUTER LOOP: SUM OF XC PERTURBATION
!     INNER LOOP: GETS ALLL THE FORCES IN RESPONSE TO EACH OUTER LOOP PERTURBATION
         ierr = 0
         DO nty = 1, 3*ntmax
            WRITE(33, '(a,i4)') "NTYPE' (XC-pert) = ",nty
         DO mt = 0, mpol1
            WRITE(33, '(a,i4)') "M' = ",mt
         DO nt = 0, ntor1-1
            WRITE(33, '(a,i4)') "N' = ",nt
         DO js = 1, ns
            ierr = ierr+1
            IF (ierr .gt. ndim) EXIT
            IF (MOD(ierr,50).eq.0) PRINT '(2x,a,f8.2,a)', 'Progress: ', 
     1                             REAL(100*ierr)/ndim, ' %'
            IF (js.eq.ns .and. .not.lfreeb) CYCLE
            colx = 1;  colb = 3
            vecs(:,colx) = 0; vecs(ierr,colx) = 1
            CALL matvec(vecs(1,colx), vecs(1,colb), ndim)
            WRITE (33, '(a,i4,2x,a,i5,2x,a,1p,e12.2)') "js' = ", js, 
     1        ' ipert = ',ierr,' Ap[ipert,ipert] = ', vecs(ierr, colb)
            colx = 0
            DO ntyp = 1, 3*ntmax
               DO mp = 0, mpol1
                  DO np = 0, ntor1-1
                     DO jp = 1, ns
                        colx = colx + 1
                        IF (colx .gt. ndim) CYCLE
                        IF (colx.eq.ierr .or. 
     1                     ABS(vecs(colx,colb)).lt.0.05_dp) CYCLE
                        WRITE (33, 123)'ntype = ', ntyp,' m = ',mp,
     1                  ' n = ', np,' js = ', jp,' iforce = ',colx,
     2                  ' Ap[iforce,ipert] = ', vecs(colx,colb)
                     END DO
                  END DO
               END DO
            END DO
         END DO
         END DO
         END DO
         END DO
 
         PRINT '(/,2x,a,/)','Jacobian check in file FORT.33'

         END IF

      END IF
 123  FORMAT(4x,4(a,i4),a,i6,a,1p,e12.3)
!
!     INITIALIZE vecs
!
      vecs(:ndim,2) = -gc(:ndim)
      vecs(:ndim,3) =  gc(:ndim)

      ierr = 100000
      info = 0
      info(1) = ierr

 10   CALL dutfx (ndim,nlen,nlim,vecs,tol,info)
      revcom = info(2)
      colx   = info(3)
      colb   = info(4)
      IF (revcom .eq. 1) THEN
         CALL matvec (vecs(1,colx), vecs(1,colb), ndim)
         GO TO 10
      END IF

      xc(1:ndim) = xsave(1:ndim) + vecs(1:ndim,1)

      END SUBROUTINE qmr_fun

      END MODULE gmres_mod

      SUBROUTINE Truncate(num, iprec0)
      USE stel_kinds, ONLY: dp
      IMPLICIT NONE
!     NEEDED TO RESOLVE CALL IN gmres_par
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)     :: iprec0
      REAL(dp), INTENT(INOUT) :: num
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER*24 :: chnum, tchnum
!-----------------------------------------------
!
!     TRUNCATES double-precision to precision iprec0 digits, keeping exponent range of double
!     WRITE TO INTERNAL FILE TO DO TRUNCATION
!
!UNCOMMENT THESE LINES TO ACTIVATE
!      WRITE (chnum, '(a,i2,a,i2,a)') '(1p,e',iprec0+7,'.',iprec0,')'
!      WRITE (tchnum, chnum) num

!      READ (tchnum, chnum) num

      END SUBROUTINE Truncate
