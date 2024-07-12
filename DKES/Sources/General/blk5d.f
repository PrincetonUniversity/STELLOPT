      SUBROUTINE blk5d(a, bm1, bp1, cp2, cm2back, cm2, pl2back,
     1   fz1, fz3, srces)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE stel_kinds
      USE Vimatrix
      USE Vnamecl2
      USE dkes_input, ONLY: idisk, lalpha
      USE dkes_realspace, ONLY: mn0, mpnt, mpntsq, diagle, diagl
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: bytes_per_rprec = 8
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(mpnt,mpnt), INTENT(out) ::
     1     a, bp1, bm1, cp2, cm2back, cm2, pl2back
      REAL(rprec), DIMENSION(mpnt,0:lalpha), INTENT(out) :: fz1, fz3
      REAL(rprec), DIMENSION(mpnt,lsource,2,2), INTENT(in) :: srces
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ibuph, incnow, irecl, incbu, iunit
      INTEGER :: kbot, k, ll, mblk2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ipiv
      REAL(rprec) :: fac
      REAL(rprec), DIMENSION(:,:), POINTER :: ql, pl, ql2back
      REAL(rprec), DIMENSION(:,:,:), POINTER :: transf
      LOGICAL :: ldisk
C-----------------------------------------------
c
c  original version:                   W.I. van Rij
c                         in the variational version of DKES code
c
c  modified for CRAY-XMP (Garching):   H. Maassberg       May 89
c                  ---->  FORTRAN SUBROUTINE MC32AD
c                  ---->  simulation of BLAS routines SGEMM/SGEMV
c                         by routine SAXPY
c                  ---->  disk i/o by routines WRDISK / RDDISK
c
c  modified (May, 1999, ORNL):         S. P. Hirshman
c                  Removed MC32AD, replaced with F90 MATMUL routines
c
c  significantly modified (June, 2001, ORNL) S. P. Hirshman
c                  fixed bug associated with particle conservation
c                  rewrote block pentadiagonal solver in readable fashion
c                  and to incorporate l=0 constraints as a "ghost" l=-1 component
c                  rewrote blox routine
c
c-----------------------------------------------------------------------
c
c  this subroutine solves the block-pentadiagonal system of equations.
c  it is called once with "plus" indexing for fz1,fz3 and once with "minus"
c  indexing, corresponding to maximizing, minimizing bounding distributions
c  for use in the variational equations
c
c-----------------------------------------------------------------------
c
c  iunit       : unit number for block-pentadiagonal solution disk file.
c  fz1(mn,l)   : distribution function in response to density gradients (1)
c  fz3(mn,l)   : distribution function in response to parallel electric field (3)
c                Note: the l=0 (real l=-1) component is the associated f(-,+)(l=0)
c                component needed for particle conservation constraint, previously
c                call fo1 and fo3
c
c  srces(mn,l,itype,plus-min)
c              : sources, mn=fourier index;       l=legendre index +1 (that is leg=0 => l=1 here)
c                itype=1,3 (uses 2) (n'',<e.b>);   plus(+) or minus(-) type (max,min)
c
c  distributions are indexed in m-n fourier-space, legendre-space. the penta-diagonal
c  equation is:
c
c  cm2 * f(l-2) + bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) + cp2 * f(l+2) = source(l)
c
c  where:   bp1(l) = bm1(l+1) (transpose);     cp1(l) = cm1(l+2) (transpose)
c
c
c
c     GENERAL SOLUTION SCHEME APPLIED TO EACH BLOCK (L)
c
c     1. CALL dgetrf:   Perform LU factorization of diagonal block (A)
c     2. CALL dgetrs:   With multiple (mpnt) right-hand sides, to do block inversion
c                        operation, A X = B  (store result in X; here B is a matrix)
c     3. CALL dgetrs:   With single right hand side (source) to solve A x = b (b a vector)
c

      mblk2 = 2*mpntsq
      ALLOCATE (transf(mpnt,mpnt,2), ql2back(mpnt,mpnt), ipiv(mpnt),
     1          stat=ier)
      IF (ier .ne. 0) STOP 'Allocation error in blk5'

      pl => transf(:,:,1)
      ql => transf(:,:,2)

      kbot = idisk*lam6 + 1


c  create disk file for doing direct access i/o.

      incnow = mblk2
      irecl  = bytes_per_rprec*incnow
      incbu = 1 + (mblk2 - 1)/incnow
      ibuph = 0

      iunit = 10
      CALL safe_open(iunit, ier, 'NULL', 'scratch', 'unformatted',
     1     irecl, 'DIRECT')
      IF (ier .ne. 0) STOP 'Error opening scratch file in blk5d DKES'

c  load and process first block-row. -----------------------------------

      k = 1
      ll = lalpha

      CALL blox (ll, a, bm1, cm2, fz1(1,ll), fz3(1,ll), srces)
      ql = bm1
      pl = cm2
      cm2back = cm2

!
!     Compute (and save) ql = A-1 ql,  pl = A-1 pl, and source terms A-1 fz1,3 [A-1 == inv(A)]
!
      CALL dgetrf (mpnt, mpnt, a, mpnt, ipiv, ier)
      IF (ier .ne. 0) GOTO 200

      CALL dgetrs ('n', mpnt, mpnt, a, mpnt, ipiv, ql, mpnt, ier)
      CALL dgetrs ('n', mpnt, mpnt, a, mpnt, ipiv, pl, mpnt, ier)

      CALL dgetrs ('n', mpnt, 1, a, mpnt, ipiv, fz1(:,ll), mpnt, ier)
      CALL dgetrs ('n', mpnt, 1, a, mpnt, ipiv, fz3(:,ll), mpnt, ier)
!
!     save pl as pl2back (will use at l+2 iteration)
!
      pl2back = pl
      ql2back = ql

      ldisk = (idisk.eq.0) .or. ((idisk.eq.1) .and. (lalpha.le.7))

      IF (ldisk) CALL wrdisk(iunit, transf, mblk2,
     1          incnow, ibuph, incbu, ier)
      IF (ier .ne. 0) GOTO 302


c  load and process second block-row. ----------------------------------

      k = k+1
      ll = ll-1
      bp1 = TRANSPOSE(bm1)

      CALL blox (ll, a, bm1, cm2, fz1(1,ll), fz3(1,ll), srces)

      a   = a - MATMUL(bp1, ql)
      ql  = bm1 - MATMUL(bp1, pl)
      pl  = cm2

      fz1(:,ll) = fz1(:,ll) - MATMUL(bp1, fz1(:,ll+1))
      fz3(:,ll) = fz3(:,ll) - MATMUL(bp1, fz3(:,ll+1))

      CALL dgetrf (mpnt, mpnt, a, mpnt, ipiv, ier)
      IF (ier .ne. 0) GOTO 200

      CALL dgetrs ('n', mpnt, mpnt, a, mpnt, ipiv, ql, mpnt, ier)
      CALL dgetrs ('n', mpnt, mpnt, a, mpnt, ipiv, pl, mpnt, ier)

      CALL dgetrs ('n', mpnt, 1, a, mpnt, ipiv, fz1(:,ll), mpnt, ier)
      CALL dgetrs ('n', mpnt, 1, a, mpnt, ipiv, fz3(:,ll), mpnt, ier)

      IF (ldisk) CALL wrdisk(iunit, transf, mblk2,
     1          incnow, ibuph, incbu, ier)
      IF (ier .ne. 0) GOTO 302

c  main loop. load and process block-rows 3 to lalpha+1. The last row (k=lalph+1) corresponds
c  to the constraint function F- which has been added to the f vector to satisfy particle
c  conservation

      fac = iswpm - 1

      BLOCKS: DO k = 3, lalpha + 1
         ll = ll - 1

         IF (k .le. lalpha) THEN
            bp1 = TRANSPOSE(bm1)
            cp2 = TRANSPOSE(cm2back)
            cm2back = cm2                                  !!stored 2 blocks back

            CALL blox (ll, a, bm1, cm2, fz1(1,ll), fz3(1,ll), srces)
         ELSE
!
!         particle conservation constraint; note V(l=0)F changes parity (+,-)
!         which is why the -transpose is used here
!
            a = 0
            a(mn0, mn0) = 1
            bp1 =-TRANSPOSE(diagle(:,:,iswpm))
            cp2 =-TRANSPOSE(diagl(:,:,iswpm))                        !f(l=1) part of V(l=0)
            fz1(:,0) = fac*srces(:,1,1,1)
            fz3(:,0) = 0
         END IF

!
!      Update diagonal "a" matrix and source terms and store pl,ql 2 l-steps back
!
         bp1 = bp1 - MATMUL(cp2, ql2back)
         a   = a - MATMUL(bp1, ql)
         a   = a - MATMUL(cp2, pl2back)
         ql2back = ql
         pl2back = pl

         fz1(:,ll) = fz1(:,ll) - MATMUL(bp1,fz1(:,ll+1))
     1                         - MATMUL(cp2,fz1(:,ll+2))
         fz3(:,ll) = fz3(:,ll) - MATMUL(bp1,fz3(:,ll+1))
     1                         - MATMUL(cp2,fz3(:,ll+2))

         IF (k .gt. lam1) THEN
            cm2 = 0
            bm1 = 0
         END IF
!
!        Compute (-,+)V[F(l=0)] contributions from l = 1, l = 0 Legendre moments (of V)
!
         IF (k .eq. lam1) cm2 = diagl(:,:,iswpm)
         IF (k .eq. lalpha) bm1 = diagle(:,:,iswpm)

!
!        Compute a-1; pl = A-1 * pl,  ql = A-1 * ql; sources = A-1 * sources
!
         CALL dgetrf (mpnt, mpnt, a, mpnt, ipiv, ier)
         IF (ier .ne. 0) GOTO 200

         IF (k .le. lalpha) THEN
            ql = bm1 - MATMUL(bp1, pl)
            pl = cm2
            CALL dgetrs('n', mpnt, mpnt, a, mpnt, ipiv, ql, mpnt,ier)
            CALL dgetrs('n', mpnt, mpnt, a, mpnt, ipiv, pl, mpnt,ier)

            CALL dgetrs('n',mpnt,1,a,mpnt, ipiv, fz1(:,ll), mpnt, ier)
            CALL dgetrs('n',mpnt,1,a,mpnt, ipiv, fz3(:,ll), mpnt, ier)

            ldisk = (idisk.eq.0) .or. ((idisk.eq.1) .and. (k.gt.lam6))
            IF (ldisk) CALL wrdisk(iunit, transf, mblk2,
     1          incnow, ibuph, incbu, ier)
            IF (ier .ne. 0) GOTO 302
         ELSE
            CALL dgetrs('n',mpnt,1,a,mpnt, ipiv, fz1(:,ll), mpnt, ier)
            CALL dgetrs('n',mpnt,1,a,mpnt, ipiv, fz3(:,ll), mpnt, ier)
         END IF

      END DO BLOCKS


c  backward solution sweep for block-rows ll = 1 (l=0) to ll = lap1-kbot

      DO k = lalpha, kbot, -1
         ll = lap1 - k
c  read blocks transf => (pl,ql) from disk.

         CALL rddisk (iunit, transf, mblk2, incnow, ibuph, ier)
         IF (ier .ne. 0) THEN
            WRITE (ioout, '(a)')
     1         ' BLK5D:   error in I/O routine RDDISK'
            GOTO 303
         ENDIF
         ibuph = ibuph - incbu


         fz1(:,ll) = fz1(:,ll) - MATMUL(ql,fz1(:,ll-1))
         fz3(:,ll) = fz3(:,ll) - MATMUL(ql,fz3(:,ll-1))

         IF (ll .ge. 2) THEN
            fz1(:,ll) = fz1(:,ll) - MATMUL(pl,fz1(:,ll-2))
            fz3(:,ll) = fz3(:,ll) - MATMUL(pl,fz3(:,ll-2))
         END IF
      END DO

      GOTO 400

c  error returns. ------------------------------------------------------

  200 CONTINUE
      ier = ll-1
      GOTO 400
  301 CONTINUE
      WRITE (ioout, '(a,i8)') ' BLK5D:   error in opening file:  ',
     1   'RECL = ', irecl
  302 CONTINUE
      WRITE (ioout, '(a)') ' BLK5D:   error in I/O routine WRDISK'
  303 CONTINUE
      ier = -2
  305 CONTINUE
      WRITE (ioout, '(2/a,i4,2/)') ' BLK5D:   error detected:   ier =',
     1   ier
      STOP

c  destroy disk file and return. ---------------------------------------

  400 CONTINUE

      CLOSE (iunit)

      DEALLOCATE (transf, ql2back, ipiv)

      END SUBROUTINE blk5d
