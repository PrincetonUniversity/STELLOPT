C*
C*  Copyright (C) CERFACS 1998
C*
C*  SOFTWARE LICENSE AGREEMENT NOTICE - THIS SOFTWARE IS BEING PROVIDED TO
C*  YOU BY CERFACS UNDER THE FOLLOWING LICENSE. BY DOWN-LOADING, INSTALLING
C*  AND/OR USING THE SOFTWARE YOU AGREE THAT YOU HAVE READ, UNDERSTOOD AND
C*  WILL COMPLY WITH THESE FOLLOWING TERMS AND CONDITIONS.
C*
C*  1 - This software program provided in source code format ("the " Source
C*  Code ") and any associated documentation (the " Documentation ") are
C*  licensed, not sold, to you.
C*j
C*  2 - CERFACS grants you a personal, non-exclusive, non-transferable and
C*  royalty-free right to use, copy or modify the Source Code and
C*  Documentation, provided that you agree to comply with the terms and
C*  restrictions of this agreement. You may modify the Source Code and
C*  Documentation to make source code derivative works, object code
C*  derivative works and/or documentation derivative Works (called "
C*  Derivative Works "). The Source Code, Documentation and Derivative
C*  Works (called " Licensed Software ") may be used by you for personal
C*  and non-commercial use only. " non-commercial use " means uses that are
C*  not or will not result in the sale, lease or rental of the Licensed
C*  Software and/or the use of the Licensed Software in any commercial
C*  product or service. CERFACS reserves all rights not expressly granted
C*  to you. No other licenses are granted or implied.
C*
C*  3 - The Source Code and Documentation are and will remain the sole
C*  property of CERFACS. The Source Code and Documentation are copyrighted
C*  works. You agree to treat any modification or derivative work of the
C*  Licensed Software as if it were part of the Licensed Software itself.
C*  In return for this license, you grant CERFACS a non-exclusive perpetual
C*  paid-up royalty-free license to make, sell, have made, copy, distribute
C*  and make derivative works of any modification or derivative work you
C*  make of the Licensed Software.
C*
C*  4- The licensee shall acknowledge the contribution of the Source Code
C*  (using the reference [1]) in any publication of material dependent upon
C*  upon the use of the Source Code. The licensee shall use reasonable
C*  endeavours to notify the authors of the package of this publication.
C*
C*  [1] V. Frayssé, L. Giraud, S. Gratton, and J. Langou, A set of GMRES 
C*    routines for real and complex arithmetics on high performance
C*    computers, CERFACS Technical Report TR/PA/03/3, public domain software
C*    available on www.cerfacs/algor/Softs, 2003
C*
C*  5- CERFACS has no obligation to support the Licensed Software it is
C*  providing under this license.
C*
C*  THE LICENSED SOFTWARE IS PROVIDED " AS IS " AND CERFACS MAKE NO
C*  REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED. BY WAY OF EXAMPLE,
C*  BUT NOT LIMITATION, CERFACS MAKE NO REPRESENTATIONS OR WARRANTIES OF
C*  MERCHANTIBILY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF
C*  THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE ANY THIRD
C*  PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. CERFACS WILL NOT
C*  BE LIABLE FOR ANY CONSEQUENTIAL, INCIDENTAL, OR SPECIAL DAMAGES, OR ANY
C*  OTHER RELIEF, OR FOR ANY CLAIM BY ANY THIRD PARTY, ARISING FROM YOUR
C*  USE OF THE LICENSED SOFTWARE.
C*
C*  6- For information regarding a commercial license for the Source Code
C*  and Documentation, please contact Mrs Campassens (campasse@cerfacs.fr)
C*
C*  7- This license is effective until terminated. You may terminate this
C*  license at any time by destroying the Licensed Software.
C*
C*    I agree all the terms and conditions of the above license agreement
C*
*
        subroutine drive_dgmres(n,nloc,m,lwork,work,
     &                         irc,icntl,cntl,info,rinfo)
        USE stel_kinds, ONLY: rprec
*
*  Purpose
*  =======
*    drive_dgmres is the driver routine for solving the linear system 
*  Ax = b using the *  Generalized Minimal Residual iterative method 
*  with preconditioning.
*  This solver is implemented with a reverse communication scheme: control
*  is returned to the user for computing the (preconditioned) 
*  matrix-vector product.
*  See the User's Guide for an example of use.
*
*
* Written : June 1996
* Authors : Luc Giraud, Serge Gratton, V. Fraysse
*             Parallel Algorithms - CERFACS
*
* Updated : April 1997
* Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
*             Parallel Algorithms - CERFACS
*
* Updated : June 1998
* Authors : Valerie Fraysse, Luc Giraud, Serge Gratton
*             Parallel Algorithms - CERFACS
* Purpose : Make clear that the warning and error messages come from the
*           dgmres modules.
*
* Updated : December 2002 - L. Giraud, J.Langou
* Purpose : Add the capability to avoid explicit residual calculation at restart
*
*
*  Arguments
*  =========
*
*  n      (input) INTEGER.
*          On entry, the dimension of the problem.
*          Unchanged on exit.
*
*  nloc   (input) INTEGER.
*          On entry, the dimension of the local problem.
*          In a parallel distributed envirionment, this corresponds
*          to the size of the subset of entries of the right hand side
*          and solution allocated to the calling process.
*          Unchanged on exit.
*
*
*  m      (input) INTEGER
*          Restart parameter, <= N. This parameter controls the amount
*          of memory required for matrix H (see WORK and H).
*          Unchanged on exit.
*
*  lwork  (input) INTEGER
*          size of the workspace
*          lwork >= m*m + m*(n+5) + 5*n+1, if icntl(8) = 1
*          lwork >= m*m + m*(n+5) + 6*n+1, if icntl(8) = 0
*
*  work   (workspace) real*8/real*8 array, length lwork
*          work contains the required vector and matrices stored in the 
*          following order :
*            x  (n,1)       : computed solution.
*            b  (n,1)       : right hand side.
*            r0 (n,1)       : vector workspace.
*            w  (n,1)       : vector workspace.
*            V  (n,m)       : Krylov basis.
*            H  (m+1,m+1)   : Hessenberg matrix (full storage).
*            yCurrent (m,1) : solution of the current LS
*            xCurrent (n,1) : current iterate
*            rotSin (m,1)   : Sine of the Givens rotation
*            rotCos (m,1)   : Cosine of the Givens rotation
*
*  irc     (input/output) INTEGER array. length 5
*            irc(1) : REVCOM   used for reverse communication
*                             (type of external operation)
*            irc(2) : COLX     used for reverse communication
*            irc(3) : COLY     used for reverse communication
*            irc(4) : COLZ     used for reverse communication
*            irc(5) : NBSCAL   used for reverse communication
*
*  icntl   (input) INTEGER array. length 7
*            icntl(1) : stdout for error messages
*            icntl(2) : stdout for warnings
*            icntl(3) : stdout for convergence history
*            icntl(4) : 0 - no preconditioning
*                       1 - left preconditioning
*                       2 - right preconditioning
*                       3 - double side preconditioning
*                       4 - error, default set in Init
*            icntl(5) : 0 - modified Gram-Schmidt
*                       1 - iterative modified Gram-Schmidt
*                       2 - classical Gram-Schmidt
*                       3 - iterative classical Gram-Schmidt
*            icntl(6) : 0 - default initial guess x_0 = 0 (to be set)
*                       1 - user supplied initial guess
*            icntl(7) : maximum number of iterations
*            icntl(8) : 0 - use recurence formula at restart
*                       1 - default compute the true residual at each restart
*
*  cntl    (input) real*8 array, length 5
*            cntl(1) : tolerance for convergence
*            cntl(2) : scaling factor for normwise perturbation on A
*            cntl(3) : scaling factor for normwise perturbation on b
*            cntl(4) : scaling factor for normwise perturbation on the
*                      preconditioned matrix
*            cntl(5) : scaling factor for normwise perturbation on 
*                      preconditioned right hand side
*
*  info    (output) INTEGER array, length 3
*            info(1) :  0 - normal exit
*                      -1 - n < 1
*                      -2 - m < 1
*                      -3 - lwork too small
*                      -4 - convergence not achieved after icntl(7) iterations
*                      -5 - precondition type not set by user
*            info(2) : if info(1)=0 - number of iteration to converge
*                      if info(1)=-3 - minimum workspace size necessary
*            info(3) : optimal size for the workspace
*
*  rinfo   (output) real*8 array, length 2
*            if info(1)=0 
*              rinfo(1) : backward error for the preconditioned system
*              rinfo(2) : backward error for the unpreconditioned system
*
* Input variables
* ---------------
       integer n, nloc, lwork, icntl(*)
       real*8   cntl(*)
       real*8   sA, sb, sPA, sPb
* Output variables
* ----------------
       integer  info(*)
       real*8    rinfo(*)
* Input/Output variables
* ----------------------
       integer  m, irc(*)
       real*8 work(*)
* Local variables
* ---------------
       integer xptr, bptr, wptr, r0ptr, Vptr, Hptr
       integer yCurrent,rotSin, rotCos, xCurrent
       integer sizeWrk, newRestart
       integer iwarn, ierr, ihist, compRsd
!SPH       real   rn
       real*8   rn
       real*8 DZRO
       parameter (DZRO = 0.0d0)
*
       integer icheck
       DATA icheck /0/
       save icheck
*
!     SPH: ifix -> int;  float -> real
!       integer   ifix
!       intrinsic ifix
!       real      float
!       real*8    float
!       intrinsic float
*
*       Executable statements :
*
       ierr    = icntl(1)
       iwarn   = icntl(2)
       ihist   = icntl(3)
       compRsd = icntl(8)
*
       if (ierr.lt.0) ierr = 6
*
       if (compRsd.eq.1) then
          sizeWrk  = m*m + m*(nloc+5) + 5*nloc+ 1
       else
          sizeWrk  = m*m + m*(nloc+5) + 6*nloc+ 1
       endif
*
       if (icheck.eq.0) then
* Check the value of the arguments
         if ((n.lt.1).or.(nloc.lt.1)) then
            write(ierr,*)
            write(ierr,*)' ERROR GMRES : '
            write(ierr,*)'     N < 1 '
            write(ierr,*)
            info(1) = -1
            irc(1)  = 0
            return
         endif
         if (m.lt.1) then
            write(ierr,*)
            write(ierr,*)' ERROR GMRES :'
            write(ierr,*)'     M < 1 '
            write(ierr,*)
            info(1) = -2
            irc(1)  = 0
            return
         endif
         if ((icntl(4).ne.0).and.(icntl(4).ne.1).and.
     &     (icntl(4).ne.2).and.(icntl(4).ne.3)) then
            write(ierr,*)
            write(ierr,*)' ERROR GMRES : '
            write(ierr,*)'     Undefined preconditioner '
            write(ierr,*)
            info(1) = -5
            irc(1)  = 0
            return
         endif
*
         if (iwarn.ne.0) then
           write(iwarn,*)
           write(iwarn,*) ' WARNING GMRES : '
           write(iwarn,*) '       For M = ',m,' optimal value '     
           write(iwarn,*) '       for LWORK =  ', sizeWrk
           write(iwarn,*)
         endif
*
         if ((icntl(5).lt.0).or.(icntl(5).gt.3)) then
           icntl(5) = 0
           if (iwarn.ne.0) then
             write(iwarn,*)
             write(iwarn,*) ' WARNING  GMRES : '
             write(iwarn,*) '       Undefined orthogonalisation '  
             write(iwarn,*) '       Default MGS '
             write(iwarn,*)
           endif
         endif
         if ((icntl(6).ne.0).and.(icntl(6).ne.1)) then
           icntl(6) = 0
           if (iwarn.ne.0) then
             write(iwarn,*)
             write(iwarn,*) ' WARNING GMRES : '
             write(iwarn,*) '       Undefined intial guess '
             write(iwarn,*) '       Default x0 = 0 '
             write(iwarn,*)
           endif
         endif
         if (icntl(7).le.0) then
           icntl(7) = n
           if (iwarn.ne.0) then
             write(iwarn,*)
             write(iwarn,*) ' WARNING GMRES :'
             write(iwarn,*) '       Negative max number of iterations'  
             write(iwarn,*) '       Default N '
             write(iwarn,*)
           endif
         endif
         if ((icntl(8).ne.0).and.(icntl(8).ne.1)) then
           icntl(8) = 1
           write(iwarn,*)
           write(iwarn,*) ' WARNING GMRES :'
           write(iwarn,*) '       Undefined strategy for the residual'
           write(iwarn,*) '       at restart'
           write(iwarn,*) '       Default 1 '
           write(iwarn,*)
         endif
* Check if the restart parameter is correct and if the size of the
*  workspace is big enough for the restart.
* If not try to fix correctly the parameters
*
         if ((m .gt. n).or.(lwork.lt.sizeWrk)) then
           if (m .gt. n) then
             m = n
             if (iwarn.ne.0) then
               write(iwarn,*)
               write(iwarn,*) ' WARNING GMRES : '
               write(iwarn,*) '       Parameter M bigger than N'  
               write(iwarn,*) '       New value for M ',m
               write(iwarn,*)
             endif
             if (compRsd.eq.1) then
               sizeWrk = m*m + m*(nloc+5) + 5*nloc+1
             else
               sizeWrk = m*m + m*(nloc+5) + 6*nloc+1
             endif
           endif
           if ((lwork.lt.sizeWrk).and.(n.eq.nloc)) then
* Compute the maximum size of the m according to the memory space
             rn         = real(n,rprec)
* SPH replaced ifix with generic int
             newRestart = int((-5-rn+sqrt((rn+5)**2-4*(5*rn
     &                   +1-real(lwork,rprec))))/2)
             if (compRsd.eq.0) then
                newRestart = newRestart - 1
             endif
             if (newRestart.gt.0) then
               m = newRestart
               if (iwarn.ne.0) then
                 write(iwarn,*)
                 write(iwarn,*)' WARNING GMRES : '
                 write(iwarn,*)'       Workspace too small for M'  
                 write(iwarn,*)'       New value for M ',m
                 write(iwarn,*)
               endif
             else
               write(ierr,*)
               write(ierr,*)' ERROR GMRES : '
               write(ierr,*)'     Not enough space for the problem'
               write(ierr,*)'     the space does not permit any m'
               write(ierr,*)
               info(1) = -3
               irc(1)  = 0
               return
             endif
           endif
           if ((lwork.lt.sizeWrk).and.(n.ne.nloc)) then
              write(ierr,*)
              write(ierr,*)' ERROR GMRES : '
              write(ierr,*)'     Not enough space for the problem'
              write(ierr,*)
              info(1) = -3
              irc(1)  = 0
              return
           endif
         endif
*
         info(3) = sizeWrk
         icheck = 1
*
* save the parameters the the history file
*
         if (ihist.ne.0) then
           write(ihist,'(10x,A39)') 'CONVERGENCE HISTORY FOR GMRES'
           write(ihist,*)
           write(ihist,'(A30,I2)') 'Errors are displayed in unit: ',ierr
           if (iwarn.eq.0) then
             write(ihist,'(A27)') 'Warnings are not displayed:'
           else
             write(ihist,'(A32,I2)') 'Warnings are displayed in unit: ',
     &                               iwarn
           endif 
           write(ihist,'(A13,I7)') 'Matrix size: ',n
           write(ihist,'(A19,I7)') 'Local matrix size: ',nloc
           write(ihist,'(A9,I7)') 'Restart: ',m
           if (icntl(4).eq.0) then
             write(ihist,'(A18)') 'No preconditioning'
           elseif (icntl(4).eq.1) then
             write(ihist,'(A20)') 'Left preconditioning'
           elseif (icntl(4).eq.2) then
             write(ihist,'(A21)') 'Right preconditioning'
           elseif (icntl(4).eq.3) then
             write(ihist,'(A30)') 'Left and right preconditioning'
           endif
           if (icntl(5).eq.0) then
             write(ihist,'(A21)') 'Modified Gram-Schmidt'
           elseif (icntl(5).eq.1) then
             write(ihist,'(A31)') 'Iterative modified Gram-Schmidt'
           elseif (icntl(5).eq.2) then
             write(ihist,'(A22)') 'Classical Gram-Schmidt'
           else
             write(ihist,'(A32)') 'Iterative classical Gram-Schmidt'
           endif
           if (icntl(6).eq.0) then
             write(ihist,'(A29)') 'Default initial guess x_0 = 0'
           else
             write(ihist,'(A27)') 'User supplied initial guess'
           endif
           if (icntl(8).eq.1) then
             write(ihist,'(A33)') 'True residual computed at restart'
           else
             write(ihist,'(A30)') 'Recurrence residual at restart'
           endif
           write(ihist,'(A30,I5)') 'Maximum number of iterations: ',
     &                              icntl(7)
           write(ihist,'(A27,E9.2)') 'Tolerance for convergence: ', 
     &                                cntl(1) 
* 
           write(ihist,'(A53)') 
     &       'Backward error on the unpreconditioned system Ax = b:'
           sA       = cntl(2)
           sb       = cntl(3)
           if ((sA.eq.DZRO).and.(sb.eq.DZRO)) then
             write(ihist,'(A39)') 
     $        '    the residual is normalised by ||b||'
           else
             write(ihist,'(A34,$)') '    the residual is normalised by '
             write(ihist,'(A8,E9.2,$)') '        ', sA
             write(ihist,'(A11,E9.2)') ' * ||x|| + ', sb
           endif
           sPA      = cntl(4)
           sPb      = cntl(5)
           write(ihist,'(A22,$)') 'Backward error on the '
           write(ihist,'(A41)') 
     &            'preconditioned system (P1)A(P2)y = (P1)b:'
           if ((sPA.eq.DZRO).and.(sPb.eq.DZRO)) then
             write(ihist,'(A34,$)') '    the preconditioned residual is'
             write(ihist,'(A24)') ' normalised by ||(P1)b||'
           else
             write(ihist,'(A35,$)') '    the preconditioned residual is'
             write(ihist,'(A15)') ' normalised by ' 
             write(ihist,'(A8,E9.2,A3,$)') '        ', sPA, ' * '
             write(ihist,'(A12,E9.2)') '||(P2)y|| + ', sPb
           endif
*
           write(ihist,'(A31,I7)') 'Optimal size for the workspace:',
     &                              info(3)
           write(ihist,*) 
           write(ihist,'(A32,$)') 'Convergence history: b.e. on the'
           write(ihist,'(A22)') ' preconditioned system'
           write(ihist,'(A11,$)') ' Iteration '
           write(ihist,'(A27)') '  Arnoldi b.e.     True b.e.'
         endif
*
       endif
* setup some pointers on the workspace
       xptr     = 1
       bptr     = xptr + nloc
       r0ptr    = bptr + nloc
       wptr     = r0ptr + nloc
       Vptr     = wptr + nloc
       if (compRsd.eq.1) then
          Hptr     = Vptr + m*nloc
       else
          Hptr     = Vptr + (m+1)*nloc
       endif
       yCurrent = Hptr + (m+1)*(m+1)
       xCurrent = yCurrent + m
       rotSin   = xCurrent + nloc
       rotCos   = rotSin + m
*
       call dgmres(nloc,m,work(bptr),work(xptr),
     &            work(Hptr),work(wptr),work(r0ptr),work(Vptr),
     &            work(yCurrent),work(xCurrent),
     &            work(rotSin),work(rotCos),irc,icntl,cntl,info,rinfo)
*
       if (irc(1).eq.0) then
         icheck = 0
       endif
*
       return
       end
