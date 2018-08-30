      MODULE gmres_mod
      USE vmec_main, ONLY: dp, rprec, neqs, ns, mns, mnmax, nthreed
      IMPLICIT NONE
      INTEGER :: nfcn = 0
      INTEGER :: ier_flag_res
      LOGICAL :: lqmr, lfirst
			LOGICAL, PARAMETER :: lscreen0 = .FALSE.

!
!     nfcn :  number of calls to function (funct3d)
!     lqmr :  logical, used by external programs to control calling these routines
!
      CONTAINS 

      SUBROUTINE matvec (p, Ap, ndim)
      USE stel_kinds
      USE xstuff, ONLY: xc, x0=>xsave, gc0=>xcdot, gc
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in)      :: ndim
      REAL(rprec), INTENT(in), DIMENSION(ndim)  :: p
      REAL(rprec), INTENT(out), DIMENSION(ndim) :: Ap
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      LOGICAL, PARAMETER :: lscreen = .false.
      REAL(rprec) :: delta
C-----------------------------------------------
!
!     Computes linearized matrix product A*p = [F(x0+delta*p) - F(x0)]/delta, about point x0
!     Must scale p so delta*max|p| ~ sqrt(epsilon) to get accurate numerical derivative
!
      delta = SQRT(EPSILON(delta))
      delta = delta/MAX(delta, MAXVAL(ABS(p)))

      xc(1:ndim) = x0(1:ndim) + delta*p(1:ndim)
      CALL funct3d(lscreen, ier_flag_res)
      Ap = (gc(1:ndim) - gc0(1:ndim))/delta

      IF (ier_flag_res .ne. 0) 
     1  PRINT *,' IN 2D PRECONDITIONER MATVEC, IER_FLAG = ', 
     2          ier_flag_res
      nfcn = nfcn + 1

      END SUBROUTINE matvec


      SUBROUTINE gmres_fun (ier_flag, itype)
      USE xstuff
			USE gmres_lib, only: gmres_info, gmres_ser
			USE vmec_main, ONLY: ftolv
			USE precon2d, ONLY: block_precond
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in)  :: itype
      INTEGER, INTENT(out) :: ier_flag
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero=0, one=1
      INTEGER :: n, m
      INTEGER :: icntl(9), info(3)
      REAL(rprec) :: cntl(5)
			TYPE (gmres_info) :: gi
      CHARACTER(LEN=*), PARAMETER :: qmr_message = 
     1                              'Beginning GMRES iterations'
C-----------------------------------------------
      EXTERNAL gmres
C-----------------------------------------------
!
!     STORE INITIAL POINT AND INITIAL RESIDUES (AT INIT PT)
!     SO DEVIATIONS FROM THEM CAN BE COMPUTED REPEATEDLY
!
      xcdot = gc
      xsave = xc
      n     = neqs

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
         lfirst = .false.
         WRITE (*,'(2x,a,/)') qmr_message
         WRITE (nthreed, '(2x,a,/)') qmr_message
      END IF

*******************************************************
** Initialize the control parameters to default value
*******************************************************
*
      CALL init_dgmres(icntl,cntl)
*
*************************
* Tune some parameters
*************************
*
* Tolerance
      cntl(1) = 1.e-3_dp
!      cntl(1) = 1.e-4_dp
! Write errors to fort.21
!      icntl(1) = 21   !21
! Write warnings to fort.21
      icntl(2) = 0   !21
! Save the convergence history in file fort.20
      icntl(3) = 20
! No preconditioning
      icntl(4) = 0
! ICGS orthogonalization
      icntl(5) = 3
!      icntl(5) = 0
! Initial guess
      icntl(6) = 0
!      icntl(6) = 1
! Maximum number of iterations at each step (~ns/5)
!      icntl(7) = 15
      icntl(7) = 20

*********************************
** Choose the restart parameter
*********************************
!      write(*,*) 'Restart  <', ldstrt
!      read(*,*) m
!
!     m <= n
!
      m = 20
      gc = -gc          !RHS
			! Load gmres_info structure
			info = 0
			gi%m = m
			gi%icntl=icntl
			gi%cntl=cntl
			gi%info=info
			gi%ftol=ftolv
      CALL gmres_ser(n,gi,matvec,block_precond,getnlforce,xcdot,gc)

100   CONTINUE

!      ier_flag = info(1)
      ier_flag = 0
      xc(1:n) = xsave(1:n) + xcdot(1:n)

      END SUBROUTINE gmres_fun


      SUBROUTINE gmresr_fun (ier_flag)
      USE xstuff
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ier_flag
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ndim, jtrunc, mgmres, maxits, lwork
      LOGICAL :: oktest
      REAL(rprec) :: eps, resid
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: work, delx, brhs
      CHARACTER(len=3), PARAMETER :: stc="rel"
      CHARACTER(LEN=*), PARAMETER :: qmr_message = 
     1                              'Beginning GMRESR iterations'
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
     1     eps, stc, maxits, resid, matvec, ier_flag_res)

      xc(1:ndim) = xsave(1:ndim) + delx(1:ndim)

      DEALLOCATE (work, delx, brhs)

!     ier_flag = ier_flag_res
      ier_flag = 0

      END SUBROUTINE gmresr_fun


      SUBROUTINE qmr_fun
      USE vmec_dim, ONLY: ns, mpol1, ntor1
      USE vmec_params, ONLY: ntmax
      USE vmec_main, ONLY: lfreeb
      USE xstuff
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ndim, nlen, nlim, ierr, info(4)
      INTEGER :: revcom, colx, colb
      INTEGER :: nty, ntyp, mt, mp, nt, np, jp, js
      REAL(rprec), DIMENSION(neqs,9) :: vecs
      REAL(rprec) :: tol = 1.E-3_dp
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

			SUBROUTINE GetNLForce(xcstate, fsq_nl, bnorm)
      USE xstuff, ONLY: xc, gc, x0=>xsave
			USE vmec_main, ONLY: fsql, fsqr, fsqz
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp),INTENT(IN)  :: xcstate(neqs), bnorm
      REAL(dp),INTENT(OUT) :: fsq_nl
!-----------------------------------------------
!undo internal gmres normalization
      xc(1:neqs) = x0(1:neqs)+bnorm*xcstate(1:neqs)

      CALL funct3d(lscreen0, ier_flag_res)
      fsq_nl = fsqr+fsqz+fsql

      nfcn = nfcn + 1

      END SUBROUTINE GetNLForce

      END MODULE gmres_mod
