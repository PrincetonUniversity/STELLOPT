      SUBROUTINE gmres (n, m, icntl, cntl, yAx, x0, b, info)
      USE stel_kinds, ONLY: rprec
      USE stel_constants, ONLY: one, zero
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: n, m,  icntl(8), info(3)
      REAL(rprec), INTENT(in)    :: b(n)
      REAL(rprec), INTENT(inout) :: x0(n)
      REAL(rprec) :: cntl(5)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: revcom, colx, coly, colz, nbscal
      INTEGER :: irc(5), jcount
      INTEGER, PARAMETER :: matvec=1, precondLeft=2, 
     1                      precondRight=3, dotProd=4
      INTEGER :: nout, lwork
      REAL(rprec)  :: rinfo(2)
      REAL(rprec), TARGET, ALLOCATABLE :: work(:)
      REAL(rprec), POINTER :: sx(:), sy(:), sz(:)
C-----------------------------------------------
      EXTERNAL yAx
C-----------------------------------------------
!
!     EASY-TO-USE WRAPPER FOR GMRES DRIVER CALL
!
!     X0: on input, initial guess if icntl(6) == 1
!         on output, solution of Ax = b
!         NOTE: it is not overwritten UNTIL the end of this routine
!
      lwork = m**2 + m*(n+5) + 5*n + 1
      ALLOCATE (work(lwork), stat=nout)
      IF (nout .ne. 0) STOP 'Allocation error in gmres!'
      work = 0
      IF (icntl(6) .eq. 1) work(1:n) = x0
      work(n+1:2*n) = b(1:n)

*****************************************
** Reverse communication implementation
*****************************************
*
 10   CONTINUE
      CALL drive_dgmres(n,n,m,lwork,work,
     &                  irc,icntl,cntl,info,rinfo)
      revcom = irc(1)
      colx   = irc(2)
      coly   = irc(3)
      colz   = irc(4)
      nbscal = irc(5)
      sx => work(colx:);  sy => work(coly:);  sz => work(colz:)

      IF (revcom.eq.matvec) THEN
* perform the matrix vector product
*        work(colz) <-- A * work(colx)
         CALL yAx (sx, sz, n) 
         GOTO 10
*
      ELSE IF (revcom.eq.precondLeft) THEN
* perform the left preconditioning
*         work(colz) <-- M^{-1} * work(colx)
!         CALL dcopy(n,work(colx),1,work(colz),1)
         CALL dcopy(n,sx,1,sz,1)
         GOTO 10
*
      ELSE IF (revcom.eq.precondRight) THEN
* perform the right preconditioning
!         CALL dcopy(n,work(colx),1,work(colz),1)
         CALL dcopy(n,sx,1,sz,1)
         GOTO 10
*
      ELSE IF (revcom.eq.dotProd) THEN
*      perform the scalar product
*      work(colz) <-- work(colx) work(coly)
*
         CALL dgemv('C',n,nbscal,ONE,sx,n,sy,1,ZERO,sz,1)
         GOTO 10
      ENDIF

*******************************
* dump the solution to a file for debugging
*******************************
!  JDH Commented out below 2008-05-15
!      GOTO 100
!
!      nout = 11
!      open(nout,FILE='sol_dTestgmres',STATUS='unknown')
!      if (icntl(5).eq.0) then
!        write(nout,*) 'Orthogonalization : MGS'
!      elseif (icntl(5).eq.1) then
!        write(nout,*) 'Orthogonalization : IMGS'
!      elseif (icntl(5).eq.2) then
!        write(nout,*) 'Orthogonalization : CGS'
!      elseif (icntl(5).eq.3) then
!        write(nout,*) 'Orthogonalization : ICGS'
!      endif
!      write(nout,*) 'Restart : ', m
!      write(nout,*) 'info(1) = ',info(1),'  info(2) = ',info(2)
!      write(nout,*) 'rinfo(1) = ',rinfo(1),'  rinfo(2) = ',rinfo(2)
!      write(nout,*) 'Optimal workspace = ', info(3)
!      write(nout,*) 'Solution : '
!      do jcount=1,n
!        write(nout,*) work(jcount)
!      enddo
!      write(nout,*) '   '
!*
!100   continue
!*

      x0 = work(1:n)

      DEALLOCATE (work)

      END SUBROUTINE gmres


      SUBROUTINE dgmres(n,m,b,x,H,w,r0,V,yCurrent,xCurrent,rotSin,
     &                  rotCos,irc,icntl,cntl,info,rinfo)
*
*
*  Purpose
*  =======
*  dgmres solves the linear system Ax = b using the
*  Generalized Minimal Residual iterative method
*
* When preconditioning is used we solve :
*     M_1^{-1} A M_2^{-1} y = M_1^{-1} b
*     x = M_2^{-1} y
*
*   Convergence test based on the normwise backward error for
*  the preconditioned system
*
* Written : June 1996
* Authors : Luc Giraud, Serge Gratton, V. Fraysse
*             Parallel Algorithms - CERFACS
*
* Updated : April 1997
* Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
*             Parallel Algorithms - CERFACS
*
* Updated : March 1998
* Purpose : Pb with F90 on DEC ws
*           cure : remove "ZDSCAL" when used to initialize vectors to zero
*
* Updated : May 1998
* Purpose : r0(1) <-- r0'r0 : pb when used with DGEMV for the dot product
*           cure : w(1) <--  r0'r0
*
* Updated : June 1998
* Purpose : Make clear that the warning and error messages come from the
*           dgmres modules.
*
* Updated : February 2001 - L. Giraud
* Purpose : In complex version, initializations to zero performed  in complex
*           arithmetic to avoid implicit conversion by the compiler.
*
* Updated : July 2001 - L. Giraud, J. Langou
* Purpose : Avoid to compute the approximate solution at each step of
*           the Krylov space construction when spA is zero.
*
* Updated : November 2002 - S. Gratton
* Purpose : Use Givens rotations conform to the classical definition.
*           No impact one the convergence history.
*
* Updated : November 2002 - L. Giraud
* Purpose : Properly handle the situation when the convergence is obtained
*           exactly at the "IterMax" iteration
*
* Updated : December 2002 - L. Giraud, J.Langou
* Purpose : Add the capability to avoid explicit residual calculation at restart
*
* Updated : January  2003 - L. Giraud, S. Gratton
* Purpose : Use Givens rotations from BLAS.
*
* Updated : March    2003 - L. Giraud
* Purpose : Set back retlbl to zero, if initial guess is solution
*           or right-hand side is zero
*
*  Arguments
*  =========
*
*  n       (input) INTEGER.
*           On entry, the dimension of the problem.
*           Unchanged on exit.
*
*  m        (input) INTEGER
*           Restart parameter, <= N. This parameter controls the amount
*           of memory required for matrix H (see WORK and H).
*           Unchanged on exit.
*
*  b        (input) real*8/real*8
*           Right hand side of the linear system.
*
*  x        (output) real*8/real*8
*           Computed solution of the linear system.
*
*  H        (workspace)  real*8/real*8
*           Hessenberg matrix built within dgmres
*
*  w        (workspace)  real*8/real*8
*           Vector used as temporary storage
*
*  r0       (workspace)  real*8/real*8
*           Vector used as temporary storage
*
*  V        (workspace)  real*8/real*8
*           Basis computed by the Arnoldi's procedure.
*  
*  yCurrent (workspace) real*8/real*8
*           solution of the current LS
*
*  xCurrent (workspace) real*8/real*8
*           current iterate
*
*  rotSin   (workspace) real*8/real*8
*           Sine of the Givens rotation
*
*  rotCos   (workspace) real*8
*           Cosine of the Givens rotation
*
*  irc      (input/output) INTEGER array. length 3
*             irc(1) : REVCOM   used for reverse communication
*                              (type of external operation)
*             irc(2) : COLX     used for reverse communication
*             irc(3) : COLY     used for reverse communication
*             irc(4) : COLZ     used for reverse communication
*             irc(5) : NBSCAL   used for reverse communication
*
*  icntl    (input) INTEGER array. length 7
*             icntl(1) : stdout for error messages
*             icntl(2) : stdout for warnings
*             icntl(3) : stdout for convergence history
*             icntl(4) : 0 - no preconditioning
*                        1 - left preconditioning
*                        2 - right preconditioning
*                        3 - double side preconditioning
*                        4 - error, default set in Init
*             icntl(5) : 0 - modified Gram-Schmidt
*                        1 - iterative modified Gram-Schmidt
*                        2 - classical Gram-Schmidt
*                        3 - iterative classical Gram-Schmidt
*             icntl(6) : 0 - default initial guess x_0 = 0 (to be set)
*                        1 - user supplied initial guess
*             icntl(7) : maximum number of iterations
*             icntl(8) : 1 - default compute the true residual at each restart
*                        0 - use recurence formula at restart
*
*  cntl     (input) real*8 array, length 5
*             cntl(1) : tolerance for convergence
*             cntl(2) : scaling factor for normwise perturbation on A
*             cntl(3) : scaling factor for normwise perturbation on b
*             cntl(4) : scaling factor for normwise perturbation on the
*                       preconditioned matrix
*             cntl(5) : scaling factor for normwise perturbation on
*                       preconditioned right hand side
*
*  info     (output) INTEGER array, length 2
*             info(1) :  0 - normal exit
*                       -1 - n < 1
*                       -2 - m < 1
*                       -3 - lwork too small
*                       -4 - convergence not achieved after icntl(7) iterations
*                       -5 - precondition type not set by user
*             info(2) : if info(1)=0 - number of iteration to converge
*                       if info(1)=-3 - minimum workspace size necessary
*             info(3) : optimal size for the workspace
*
* rinfo     (output) real*8 array, length 2
*             if info(1)=0 
*               rinfo(1) : backward error for the preconditioned system
*               rinfo(2) : backward error for the unpreconditioned system
*
* Input variables
* ---------------
        integer  n, m, icntl(*)
        real*8 b(*)
        real*8    cntl(*)
*
* Output variables
* ----------------
       integer  info(*)
       real*8    rinfo(*)
*
* Input/Output variables
* ----------------------
       integer  irc(*)
       real*8 x(*), H(m+1,*), w(*), r0(*), V(n,*), yCurrent(*)
       real*8 xCurrent(*), rotSin(*)
       real*8    rotCos(*)
*
* Local variables
* ---------------
       integer  j, jH, iterOut, nOrtho, iterMax, initGuess, iOrthog
       integer  xptr, bptr, wptr, r0ptr, Vptr, Hptr, yptr, xcuptr
       integer  typePrec, leftPrec, rightPrec, dblePrec, noPrec
       integer  iwarn, ihist
       integer  compRsd
       real*8    beta, bn, sA, sb, sPA, sPb, bea, be
       real*8    dloo, dnormw, dnormx, dnormres, trueNormRes
       real*8 dVi, temp, aux
       real*8 auxHjj, auxHjp1j
*
       parameter (noPrec = 0, leftPrec = 1)
       parameter (rightPrec = 2, dblePrec = 3)  
*
       real*8 ZERO, ONE
       parameter (ZERO = 0.0d0, ONE = 1.0d0)
       real*8 DZRO,DONE
       parameter (DZRO = 0.0d0, DONE = 1.0d0)
*
*
* External functions
* ------------------
       real*8    dnrm2
       external dnrm2
*
* Reverse communication variables
* -------------------------------
       integer retlbl
       DATA retlbl /0/
       integer matvec, precondLeft, precondRight, prosca
       parameter(matvec=1, precondLeft=2, precondRight=3, prosca=4)
*
* Saved variables
* ---------------
       save iterOut, jH, beta, bn, dnormres, retlbl, j
       save sA, sb, sPA, sPb, dnormx, trueNormRes, bea, be
       save dloo, nOrtho, compRsd

!     Added by SPH to converge on Arnoldi backward error (bea) rather than be
!     Set to FALSE for original behavior
 	 LOGICAL, PARAMETER :: larnoldi = .TRUE.

*
* Intrinsic function
* ------------------
       intrinsic dabs, dsqrt 
*
*       Executable statements
*
* setup some pointers on the workspace
       xptr     = 1
       bptr     = xptr + n
       r0ptr    = bptr + n
       wptr     = r0ptr + n
       Vptr     = wptr + n
       if (icntl(8).eq.1) then
         Hptr     = Vptr + m*n
       else
         Hptr     = Vptr + (m+1)*n
       endif
       yptr     = Hptr + (m+1)*(m+1)
       xcuptr   = yptr + m
*
       iwarn      = icntl(2)
       ihist      = icntl(3)
       typePrec   = icntl(4)
       iOrthog    = icntl(5)
       initGuess  = icntl(6)
       iterMax    = icntl(7)
*
       if (retlbl.eq.0) then
         compRsd    = icntl(8)
       endif
*
      if (retlbl.ne.0) then
        if (retlbl.eq.5) then
          goto 5
        else if (retlbl.eq.6) then
          goto 6
        else if (retlbl.eq.8) then
          goto 8
        else if (retlbl.eq.11) then
          goto 11
        else if (retlbl.eq.16) then
          goto 16
        else if (retlbl.eq.18) then
          goto 18
        else if (retlbl.eq.21) then
          goto 21
        else if (retlbl.eq.26) then
          goto 26
        else if (retlbl.eq.31) then
          goto 31
        else if (retlbl.eq.32) then
          goto 32
        else if (retlbl.eq.33) then
          goto 33
        else if (retlbl.eq.34) then
          goto 34 
        else if (retlbl.eq.36) then
          goto 36
        else if (retlbl.eq.37) then
          goto 37
        else if (retlbl.eq.38) then
          goto 38
        else if (retlbl.eq.41) then
          goto 41
        else if (retlbl.eq.43) then
          goto 43
        else if (retlbl.eq.46) then
          goto 46
        else if (retlbl.eq.48) then
          goto 48
        else if (retlbl.eq.51) then
          goto 51
        else if (retlbl.eq.52) then
          goto 52
        else if (retlbl.eq.61) then
          goto 61
        else if (retlbl.eq.66) then
          goto 66
        else if (retlbl.eq.68) then
          goto 68
        endif
      endif
*
*
* intialization of various variables
*
      iterOut  = 0
      beta     = DZRO
*
      if (initGuess.eq.0) then
        do j=1,n
          x(j) = ZERO
        enddo
      endif
*
*        bn = dnrm2(n,b,1)
*
      irc(1) = prosca
      irc(2) = bptr
      irc(3) = bptr
      irc(4) = r0ptr
      irc(5) = 1
      retlbl = 5
      return
 5    continue
      bn = dsqrt((r0(1)))
*
      if (bn.eq.DZRO) then
        do j=1,n
          x(j) = ZERO
        enddo  
        if (iwarn.ne.0) then
          write(iwarn,*)
          write(iwarn,*) ' WARNING GMRES : '
          write(iwarn,*) '       Null right hand side'
          write(iwarn,*) '       solution set to zero'
          write(iwarn,*)
        endif
        info(1)  = 0
        info(2)  = 0
        rinfo(1) = DZRO
        rinfo(2) = DZRO
        irc(1)   = 0
        retlbl = 0
        return
      endif
*
* Compute the scaling factor for the backward error on the 
*  unpreconditioned sytem
*
      sA       = cntl(2)
      sb       = cntl(3)
      if ((sA.eq.DZRO).and.(sb.eq.DZRO)) then
        sb = bn
      endif
* Compute the scaling factor for the backward error on the
*  preconditioned sytem
*
       sPA      = cntl(4)
       sPb      = cntl(5)
       if ((sPA.eq.DZRO).and.(sPb.eq.DZRO)) then
         if ((typePrec.eq.noPrec).or.(typePrec.eq.rightPrec)) then
           sPb = bn
         else
           irc(1) = precondLeft
           irc(2) = bptr
           irc(4) = r0ptr
           retlbl = 6
           return
         endif
       endif
 6     continue
       if ((sPA.eq.DZRO).and.(sPb.eq.DZRO)) then
         if ((typePrec.eq.dblePrec).or.(typePrec.eq.leftPrec)) then
*
*           sPb = dnrm2(n,r0,1)
*
           irc(1) = prosca
           irc(2) = r0ptr
           irc(3) = r0ptr
           irc(4) = wptr
           irc(5) = 1
           retlbl = 8
           return
         endif
       endif
 8     continue
       if ((sPA.eq.DZRO).and.(sPb.eq.DZRO)) then
         if ((typePrec.eq.dblePrec).or.(typePrec.eq.leftPrec)) then
           sPb = dsqrt((w(1)))
* 
         endif
       endif
*
*
* Compute the first residual
*           Y = AX : r0 <-- A x
*
* The residual is computed only if the initial guess is not zero
*
       if (initGuess.ne.0) then
         irc(1) = matvec
         irc(2) = xptr
         irc(4) = r0ptr
         retlbl = 11
         return
       endif
 11    continue
       if (initGuess.ne.0) then
         do j=1,n
           r0(j) = b(j)-r0(j)
         enddo
       else
         call dcopy(n,b,1,r0,1)
       endif 
*
* Compute the preconditioned residual if necessary
*      M_1Y = X : w <-- M_1^{-1} r0
*
       if ((typePrec.eq.noPrec).or.(typePrec.eq.rightPrec)) then
         call dcopy(n,r0,1,w,1)
       else
         irc(1) = precondLeft
         irc(2) = r0ptr
         irc(4) = wptr
         retlbl = 16
         return
       endif
 16    continue
*
*
*       beta = dnrm2(n,w,1)
*
*
       irc(1) = prosca
       irc(2) = wptr
       irc(3) = wptr
       irc(4) = r0ptr
       irc(5) = 1
       retlbl = 18
       return
 18    continue
       beta = dsqrt((r0(1)))
*
       if (beta .eq. DZRO) then
*  The residual is exactly zero : x is the exact solution
         info(1) = 0
         info(2) = 0
         rinfo(1) = DZRO
         rinfo(2) = DZRO
         irc(1)   = 0
         retlbl = 0
         if (iwarn.ne.0) then
          write(iwarn,*)
          write(iwarn,*) ' WARNING GMRES : '
          write(iwarn,*) '       Intial residual is zero'
          write(iwarn,*) '       initial guess is solution'
          write(iwarn,*)
         endif
         return
       endif
*
       aux = ONE/beta
       do j=1,n
         V(j,1) = ZERO
       enddo
       call daxpy(n,aux,w,1,V(1,1),1)
*
*       Most outer loop : dgmres iteration
*
*       REPEAT
 7     continue
*
*
       H(1,m+1)=beta
       do j=1,m
         H(j+1,m+1) = ZERO
       enddo
*
*        Construction of the hessenberg matrix WORK and of the orthogonal
*        basis V such that AV=VH 
*
       jH = 1
 10    continue
* Remark : this  do loop has been written with a while do
*          because the
*               " do jH=1,restart "
*         fails with the reverse communication.
*      do  jH=1,restart
*
*
* Compute the preconditioned residual if necessary
*
       if ((typePrec.eq.rightPrec).or.(typePrec.eq.dblePrec)) then  
*
*           Y = M_2^{-1}X : w <-- M_2^{-1} V(1,jH)
*
         irc(1) = precondRight
         irc(2) = vptr + (jH-1)*n
         irc(4) = wptr
         retlbl = 21
         return
       else
         call dcopy(n,V(1,jH),1,w,1)
       endif
 21    continue
*
*           Y = AX : r0 <-- A w
*
       irc(1) = matvec
       irc(2) = wptr
       irc(4) = r0ptr
       retlbl = 26
       return
 26    continue
*
*      MY = X : w <-- M_1^{-1} r0
*
       if ((typePrec.eq.noPrec).or.(typePrec.eq.rightPrec)) then
         call dcopy(n,r0,1,w,1)
       else
         irc(1) = precondLeft
         irc(2) = r0ptr
         irc(4) = wptr
         retlbl = 31
         return
       endif
 31    continue
*
* Orthogonalization using either MGS or IMGS
*  
* initialize the Hessenberg matrix to zero in order to be able to use
*     IMGS as orthogonalization procedure.
       do j=1,jH
         H(j,jH) = ZERO
       enddo
       nOrtho = 0
 19    continue
       nOrtho = nOrtho +1
       dloo   = DZRO
*
       if ((iOrthog.eq.0).or.(iOrthog.eq.1)) then
* MGS
*
*           do j=1,jH
*
         j = 1
*           REPEAT
       endif
 23    continue
       if ((iOrthog.eq.0).or.(iOrthog.eq.1)) then
*
*             dVi     = ddot(n,V(1,j),1,w,1)
*
         irc(1) = prosca
         irc(2) = vptr + (j-1)*n
         irc(3) = wptr
         irc(4) = r0ptr
         irc(5) = 1
         retlbl = 32
         return
       endif
 32    continue
       if ((iOrthog.eq.0).or.(iOrthog.eq.1)) then
         dVi     = r0(1)
         H(j,jH) = H(j,jH) + dVi
         dloo    = dloo + dabs(dVi)**2
         aux = -ONE*dVi
         call daxpy(n,aux,V(1,j),1,w,1)
         j = j + 1
         if (j.le.jH) goto 23
*          enddo_j
       else
* CGS
* produit scalaire groupe
*
*           call dgemv('C',n,jH,ONE,V(1,1),n,w,1,ZERO,r0,1)
*
         irc(1) = prosca
         irc(2) = vptr
         irc(3) = wptr
         irc(4) = r0ptr
         irc(5) = jH
         retlbl = 34
         return
       endif
 34    continue
       if ((iOrthog.eq.2).or.(iOrthog.eq.3)) then
*
         call daxpy(jH,ONE,r0,1,H(1,jH),1)
         call dgemv('N',n,jH,-ONE,V(1,1),n,r0,1,ONE,w,1)
         dloo = dnrm2(jH,r0,1)**2
       endif
*
*         dnormw = dnrm2(n,w,1)
*
       irc(1) = prosca
       irc(2) = wptr
       irc(3) = wptr
       irc(4) = r0ptr
       irc(5) = 1
       retlbl = 33
       return
 33    continue
       dnormw = dsqrt((r0(1)))
*
       if ((iOrthog.eq.1).or.(iOrthog.eq.3)) then
* IMGS / CGS orthogonalisation
         dloo = dsqrt(dloo)
* check the orthogonalization quality
         if ((dnormw.le.dloo).and.(nOrtho.lt.3)) then
           goto 19
         endif
       endif
*
       H(jH+1,jH) = dnormw
       if ((jH.lt.m).or.(icntl(8).eq.0)) then
         aux = ONE/dnormw
         do j=1,n
           V(j,jH+1) = ZERO
         enddo
         call daxpy(n,aux,w,1,V(1,jH+1),1)
       endif
* Apply previous Givens rotations to the new column of H
       do j=1,jH-1
         call drot(1, H(j,jH), 1, H(j+1,jH), 1, rotCos(j), rotSin(j))
       enddo
       auxHjj = H(jH,jH)
       auxHjp1j= H(jH+1,jH)
       call drotg(auxHjj, auxHjp1j, rotCos(jH),rotSin(jH))
* Apply current rotation to the rhs of the least squares problem
       call drot(1, H(jH,m+1), 1, H(jH+1,m+1), 1, rotCos(jH),
     &            rotSin(jH))
*
* zabs(H(jH+1,m+1)) is the residual computed using the least squares
*          solver
* Complete the QR factorisation of the Hessenberg matrix by apply the current
* rotation to the last entry of the collumn
       call drot(1, H(jH,jH), 1, H(jH+1,jH), 1, rotCos(jH), rotSin(jH))
       H(jH+1,jH) = ZERO
*
* Get the Least square residual
*
       dnormres = dabs(H(jH+1,m+1))
       if (sPa.ne.DZRO) then
*
* Compute the solution of the current linear least squares problem
*
         call dcopy(jH,H(1,m+1),1,yCurrent,1)
         call dtrsv('U','N','N',jH,H,m+1,yCurrent,1)
*
* Compute the value of the new iterate 
*
         call dgemv('N',n,jH,ONE,v,n,
     &            yCurrent,1,ZERO,xCurrent,1)
*
         if ((typePrec.eq.rightPrec).or.(typePrec.eq.dblePrec)) then  
*
*         Y = M_2^{-1}X : r0 <-- M_2^{-1} xCurrent
*
           irc(1) = precondRight
           irc(2) = xcuptr
           irc(4) = r0ptr
           retlbl = 36
           return
         else
           call dcopy(n,xCurrent,1,r0,1)
         endif
       endif
 36    continue
*
*
       if (sPa.ne.DZRO) then
* Update the current solution
         call dcopy(n,x,1,xCurrent,1)
         call daxpy(n,ONE,r0,1,xCurrent,1)
*
*         dnormx = dnrm2(n,xCurrent,1)
*
         irc(1) = prosca
         irc(2) = xcuptr
         irc(3) = xcuptr
         irc(4) = r0ptr
         irc(5) = 1
         retlbl = 38
         return
       else
         dnormx    = DONE
       endif
 38    continue
       if (sPa.ne.DZRO) then
         dnormx = dsqrt((r0(1)))
       endif
*
       bea = dnormres/(sPA*dnormx+sPb)
*
* Check the convergence based on the Arnoldi Backward error for the
* preconditioned system
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then  
* 
* The Arnoldi Backward error indicates that dgmres might have converge
* enforce the calculation of the true residual at next restart
         compRsd = 1
*
*  If the update of X has not yet been performed
         if (sPA.eq.DZRO) then
*
* Compute the solution of the current linear least squares problem
*
           call dcopy(jH,H(1,m+1),1,yCurrent,1)
           call dtrsv('U','N','N',jH,H,m+1,yCurrent,1)
*
* Compute the value of the new iterate 
*
           call dgemv('N',n,jH,ONE,v,n,
     &            yCurrent,1,ZERO,xCurrent,1)
*
           if ((typePrec.eq.rightPrec).or.(typePrec.eq.dblePrec)) then
* 
*         Y = M_2^{-1}X : r0 <-- M_2^{-1} xCurrent
*
             irc(1) = precondRight
             irc(2) = xcuptr
             irc(4) = r0ptr 
             retlbl = 37
             return
           else
             call dcopy(n,xCurrent,1,r0,1)
           endif
         endif
       endif
 37    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
         if (sPA.eq.DZRO) then
* Update the current solution
            call dcopy(n,x,1,xCurrent,1)
            call daxpy(n,ONE,r0,1,xCurrent,1)
         endif
*
         call dcopy(n,xCurrent,1,r0,1)
* Compute the true residual, the Arnoldi one may be unaccurate
*
*           Y = AX : w  <-- A r0
*
         irc(1) = matvec
         irc(2) = r0ptr
         irc(4) = wptr
         retlbl = 41
         return
       endif
 41    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
*
         do j=1,n
           w(j) = b(j) - w(j)
         enddo
* Compute the norm of the unpreconditioned residual
*
*        trueNormRes = dnrm2(n,w,1)
*
         irc(1) = prosca
         irc(2) = wptr
         irc(3) = wptr
         irc(4) = r0ptr
         irc(5) = 1
         retlbl = 43
         return
       endif
 43    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
         trueNormRes = dsqrt((r0(1)))
*
         if ((typePrec.eq.leftPrec).or.(typePrec.eq.dblePrec)) then  
*
*      MY = X : r0 <-- M_1^{-1} w 
*
           irc(1) = precondLeft
           irc(2) = wptr
           irc(4) = r0ptr
           retlbl = 46
           return
         else 
           call dcopy(n,w,1,r0,1)
         endif
       endif
 46    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
*
*        dnormres = dnrm2(n,r0,1) 
*
         irc(1) = prosca
         irc(2) = r0ptr
         irc(3) = r0ptr
         irc(4) = wptr
         irc(5) = 1
         retlbl = 48
         return
       endif
 48    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
         dnormres = dsqrt((w(1)))
*
         be = dnormres/(sPA*dnormx+sPb)
* Save the backward error on a file if convergence history requested
         if (ihist.ne.0) then
           write(ihist,'(I5,11x,E9.2,$)') iterOut*m+jH,bea
           write(ihist,'(7x,E9.2)') be
         endif
*
       endif
*
*
* Check again the convergence
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then   
         if (larnoldi .or. 
     &      (be.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then   
* The convergence has been achieved, we restore the solution in x
* and compute the two backward errors.
           call dcopy(n,xCurrent,1,x,1)
*
           if (sA.ne.DZRO) then
*
*            dnormx = dnrm2(n,x,1)
*
             irc(1) = prosca
             irc(2) = xptr
             irc(3) = xptr
             irc(4) = r0ptr
             irc(5) = 1
             retlbl = 51
             return
           endif
         endif
       endif
 51    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
         if (larnoldi .or. 
     &       (be.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
           if (sA.ne.DZRO) then
             dnormx = dsqrt((r0(1)))
*
           else
             dnormx = DONE
           endif
* Return the backward errors
           rinfo(1) = be
           rinfo(2) = trueNormRes/(sA*dnormx+sb)
           if (be.le.cntl(1)) then
             info(1) = 0
             if (ihist.ne.0) then
               write(ihist,*)
               write(ihist,'(A20)') 'Convergence achieved'
             endif
           else if (be.gt.cntl(1)) then
             if (iwarn.ne.0) then
               write(iwarn,*)
               write(iwarn,*) ' WARNING GMRES : '
               write(iwarn,*) '       No convergence after '
               write(iwarn,*) iterOut*m+jH,' iterations '
               write(iwarn,*)
             endif
             if (ihist.ne.0) then
               write(ihist,*)
               write(ihist,*) ' WARNING GMRES :'
               write(ihist,*) '       No convergence after '
               write(ihist,*) iterOut*m+jH,' iterations '
               write(ihist,*)
             endif
             info(1) = -4
           endif
           if (ihist.ne.0) then
             write(ihist,'(A27,$)') 'B.E. on the preconditioned '
             write(ihist,'(A10,E9.2)') 'system:   ', rinfo(1)
             write(ihist,'(A29,$)') 'B.E. on the unpreconditioned '
             write(ihist,'(A8,E9.2)') 'system: ', rinfo(2)
           endif
           info(2) = iterOut*m+jH
           if (ihist.ne.0) then
             write(ihist,'(A10,I2)') 'info(1) = ',info(1)
             write(ihist,'(A32,I5)') 
     &                'Number of iterations (info(2)): ',info(2)  
           endif
           irc(1)  = 0
           retlbl  = 0
           return
         endif
       else
* Save the backward error on a file if convergence history requested
         if (ihist.ne.0) then
           write(ihist,'(I5,11x,E9.2,$)') iterOut*m+jH,bea
           write(ihist,'(9x,A2)') '--'
         endif
*
       endif  
*
       jH = jH + 1
       if (jH.le.m) then
         goto 10
       endif
*
       iterOut = iterOut + 1
*
* we have completed the Krylov space construction, we restart if
* we have not yet exceeded the maximum number of iterations allowed.
*
       if ((sPa.eq.DZRO).and.(bea.gt.cntl(1))) then
*
* Compute the solution of the current linear least squares problem
*
         jH = jH - 1
         call dcopy(jH,H(1,m+1),1,yCurrent,1)
         call dtrsv('U','N','N',jH,H,m+1,yCurrent,1)
*
* Compute the value of the new iterate
*
         call dgemv('N',n,jH,ONE,v,n,
     &            yCurrent,1,ZERO,xCurrent,1)
*
         if ((typePrec.eq.rightPrec).or.(typePrec.eq.dblePrec)) then
*
*         Y = M_2^{-1}X : r0 <-- M_2^{-1} xCurrent
*
           irc(1) = precondRight
           irc(2) = xcuptr
           irc(4) = r0ptr
           retlbl = 52
           return
         else
           call dcopy(n,xCurrent,1,r0,1)
         endif
       endif
 52    continue
       if ((sPa.eq.DZRO).and.(bea.gt.cntl(1))) then
* Update the current solution
         call dcopy(n,x,1,xCurrent,1)
         call daxpy(n,ONE,r0,1,xCurrent,1)
       endif
*
       call dcopy(n,xCurrent,1,x,1)
*
       if (compRsd.eq.1) then
*
* Compute the true residual
*
         call dcopy(n,x,1,w,1)
         irc(1) = matvec
         irc(2) = wptr
         irc(4) = r0ptr
         retlbl = 61
         return
       endif
 61    continue
       if (compRsd.eq.1) then
         do j=1,n
           r0(j) = b(j) - r0(j)
         enddo
*
* Precondition the new residual if necessary
*
         if ((typePrec.eq.leftPrec).or.(typePrec.eq.dblePrec)) then
*
*      MY = X : w <-- M_1^{-1} r0
*
           irc(1) = precondLeft
           irc(2) = r0ptr
           irc(4) = wptr
           retlbl = 66
           return
         else
           call dcopy(n,r0,1,w,1)
         endif
       endif
 66    continue
*
*           beta = dnrm2(n,w,1)
*
       if (compRsd.eq.1) then
         irc(1) = prosca
         irc(2) = wptr
         irc(3) = wptr
         irc(4) = r0ptr
         irc(5) = 1
         retlbl = 68
         return
       endif
 68    continue
       if (compRsd.eq.1) then
         beta = dsqrt((r0(1)))
*
       else
* Use recurrence to approximate the residual at restart
         beta = dabs(H(m+1,m+1))
* Apply the Givens rotation is the reverse order
         do j=m,1,-1
           H(j,m+1)   = ZERO
           call drot(1, H(j,m+1), 1, H(j+1,m+1), 1,
     &               rotCos(j), -rotSin(j))
         enddo
*
* On applique les vecteurs V
*
         call dgemv('N',n,m+1,ONE,v,n,H(1,m+1),1,ZERO,w,1)
*
       endif
       do j=1,n
         V(j,1) = ZERO
       enddo
       aux = ONE/beta
       call daxpy(n,aux,w,1,V(1,1),1)
*
       goto 7
*
       end
