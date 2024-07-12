      SUBROUTINE gmresr(oktest,n,j,mgmres,b,x,work,eps,stc,
     &                  maxits,resid,matvec,iflag)
C*********************************************************
C GMRESR algorithm to solve linear system Ax = b
C
C  author:
C  m.botchev, utrecht university, december 1996
C  report bugs to botchev@cwi.nl or botchev@excite.com
C
C Copyright (c) 1996 by M.A.Botchev
C Permission to copy all or part of this work is granted,
C provided that the copies are not made or distributed
C for resale, and that the copyright notice and this
C notice are retained.
C
C Details to the algorithm may be found in
C  H.A. van der Vorst, C. Vuik, "GMRESR: a Family of Nested GMRES
C  Methods", Num. Lin. Alg. Appl., vol. 1(4), 369--386 (1994)
C
C parameter list:
C oktest is TRUE if intermediate residuals should be printed
C n      == INTEGER size of problem
C j      == INTEGER truncation parameter (work with j last vectors)
C mgmres == INTEGER dimension of the invoked GMRES
C           if mgmres.eq.0 then we get GCR algorithm, the simplest
C           version of GMRESR 
C b      == real*8 righthand side vector
C x      == real*8 initial guess on input,
C           (approximate) solution on output
C work   == real*8 work space 1 of size n x (2*j+mgmres+2)
C eps    == real*8 tolerance of stopping criterion. 
C           process is stopped as soon as 
C           (1) residual norm has been decreased by factor eps, 
C           i.e.  ||res|| / ||res0|| <= eps   OR
C           (2) maximum number of iterations maxit has been performed
C stc    == CHARACTER*3
C           Determine stopping criterion (||.|| denotes the 2-norm):
C           stc='rel'    -- relative stopping crit.: ||res|| < eps*||res0||
C           stc='abs'    -- absolute stopping crit.: ||res|| < eps
C maxits == INTEGER max. no. outer_iterative_steps/truncation_length on input
C           on output it is the actual number of total iterative steps   
C resid  == real*8 residual measure (depends on stopping criterion)
C           achieved on output 
C iflag  == INTEGER on output 0 - solution found within tolerance
C                             1 - no convergence within maxits
C ----------------------------------------------------------
C subroutines used
C   matvec   == matrix-vector product y <- A*x
C  blas subroutines:
C   dscal
C   daxpy
C  blas functions:
C   ddot
C   dnrm2 
C**********************************************************

      external matvec

C     list of variables: arrays in alphabetical order,
C     then other variables in alphabetical order

      logical oktest
      character*3 stc

      integer i,iflag,its,nits,itsinn,j,k,maxits,mgmres,n

      real*8 b(n),x(n),work(n,0:(2*j+mgmres+2-1)),
     &     alpha,alphai,cknorm,ckres,ddot,dnrm2,eps,epsinn,
     &     res0,resnor,resid

C distribute work space work(n,*) among some virtual variables;
C namely, we think of columns of work as being occupied by 
C c(n,0:j-1), u(n,0:j-1), resid(n), workgmr(n,mgmres+1)
C therefore we define "shifts"

      integer c, u, workres, workgmre
C ----------------------------------------------------------------------------

      if((stc.NE.'rel').and.(stc.NE.'abs'))then
         PRINT *,'Error in VACGMRESR:'
         PRINT *,'PARAMETER STC=',stc,' SHOULD BE rel OR abs.'
         STOP
      endif

C     c occupies j columns 0...j-1:
      c = 0 
C     u occupies j columns j...2*j-1:
      u = j
C     resid occupies 1 column No. 2*j:
      workres = 2*j    
C     workgmre occupies mgmres+1 columns 2*j+1...2*j+mgmres+1:
      workgmre = 2*j+1 
C     so that we can access, say, to the (k)th column of the virtual
C     array c(n,0:j-1) as work(1,c+k),
C     virtual residual vector resid(n) is work(1,workres) and so on ...

C ***Furthermore, we build sequences c_k and u_k, k = 0,...,m-1
C but we store only the last j vectors in the following way:
C Assume j=3, then
C --------------------------------------------------------------
C  k    |  number of column of work | columns of work which are vectors
C       |  in which we store c_k    |  we actually store and use
C  0    |           0               |   c_0             u_0            ...
C  1    |           1               |   c_0  c_1        u_0  u_1       ...
C  2    |           2               |   c_0  c_1  c_2   u_0  u_1  u_2  ...
C  3    |           0               |   c_3  c_1  c_2   u_3  u_1  u_2  ...
C  4    |           1               |   c_3  c_4  c_2   u_3  u_4  u_2  ... 
C  5    |           2               |   c_3  c_4  c_5   u_3  u_4  u_5  ...
C  6    |           0               |   c_6  c_4  c_5   u_6  u_4  u_5  ...
C ...   |           ...             |      ...               ...
C This mapping is performed by function mod(k,j)

C     Reset iteration counter
      nits= 0
      its = 0
      
C     Calculate (initial) residual norm
      call matvec(x,work(1,workres),n)
      alpha = -1
      call daxpy(n,alpha,b,1,work(1,workres),1)
      call dscal(n,alpha,work(1,workres),1)

C     Calculate residual norm and quit if it is zero
      res0 = dnrm2(n,work(1,workres),1)
      resnor = res0
      resid = 0

      if ( res0 .eq. 0.0d0 ) then
         iflag = 0
         maxits = 0
         return  
      end if

      if (stc.eq.'abs') then
         resid=resnor
      else
         resid=resnor/res0
      endif

      if ( resid .le. eps ) then
         iflag = 0
         maxits = 0
         return
      end if 
 
C     Main iterative loop ============================
      k = -1
      do while (.true.)

         if(oktest)write(*,199)its,resid
 199     format('   its =', i4, ' resid =', d20.6)

C        Loop to increment dimension of projection subspace
         k=k+1

C        Number of step (not restart) to be done
         its = its + 1
C        write(*,'(A,i3)') '+++++++++++++++++ its ',its 

C        - - - - - - - - - - - - - - - - - - - - - - - - -
C        This part should deliver 
C        u(1,k) <-- invA * resid
C        where u(1,k) is the k-th column of array u(1:n,0:m) and
C        invA is some reasonable approximation to the inverse of A
C
C        If mgmres=0 then no inner iterations are performed 
C        to get invA, so that invA is just the identity matrix. 
C        In this case algorithm GMRESR is nothing but GCR
C
C        Otherwise for inner iterations we perform ONE restart of GMRES
C        ATTN: for this implementation it is crucial to perform only
C        ONE restart of GMRES
         if (mgmres.eq.0) then
C           u(1,k) := resid  
            call dcopy(n,work(1,workres),1,work(1,u+mod(k,j)),1)
            call matvec(work(1,u+mod(k,j)),work(1,c+mod(k,j)),n)
            nits=nits+1
         else
C           Solve linear system A*u(1,k)=resid by GMRES
C           The stopping criterion for inner iterations is 
C           always absolute but it is automatically adjusted
C           not to be stricter than the stopping criterion for the 
C           outer iterations.  For example, if stop.criterion for
C           the outer iterations is relative than absolute stop.
C           criterion for the inner iterations is (eps*res0)
C           Accuracy for inner iteration:

            if(stc.eq.'abs')then
               epsinn = eps
            else
               epsinn = eps*res0
            endif

C           After envoking gmres0 epsinn and itsinn contain actual achieved
C           accuracy and number of performed iterations respectively

            itsinn=mgmres

            call gmres0(oktest,n,mgmres,
     &           work(1,workres),work(1,u+mod(k,j)),
     &           work(1,c+mod(k,j)),work(1,workgmre),
     &           epsinn,itsinn,matvec)

            nits=nits+itsinn
         end if           
C - - - - - - - - - - - - - - - - - - - - - - - - 
      
C        Inner loop to orthogonalize 
C        c(1,k) with respect to c(1,k-j),...,c(1,k-1)
C        and to update correspondingly 
C        u(1,k) with respect to u(1,k-j),...,u(1,k-1)
C        parameter j is used only here
         do i = max0(0,k-j),k-1
            alphai = ddot(n,work(1,c+mod(i,j)),1,work(1,c+mod(k,j)),1)
            call daxpy(n,-alphai,work(1,c+mod(i,j)),1,
     &           work(1,c+mod(k,j)),1)
            call daxpy(n,-alphai,work(1,u+mod(i,j)),1,
     &           work(1,u+mod(k,j)),1)
         end do

C        Normalize c(1,k) and "normalize" u(1,k)
         cknorm = dnrm2(n,work(1,c+mod(k,j)),1)
         cknorm = 1 / cknorm
         call dscal(n,cknorm,work(1,c+mod(k,j)),1)
         call dscal(n,cknorm,work(1,u+mod(k,j)),1)

C        Update current solution and residual
         ckres = ddot(n,work(1,c+mod(k,j)),1,work(1,workres),1)
         call daxpy(n, ckres,work(1,u+mod(k,j)),1,x,          1)
         call daxpy(n,-ckres,work(1,c+mod(k,j)),1,work(1,workres),1)

C        call show(n,10,x,'GMRESR       ')  

C        Calculate residual norm, check convergence
         resnor = dnrm2(n,work(1,workres),1)

         if (stc.eq.'abs') then
            resid=resnor
         else
            resid=resnor/res0
         endif

         if ( resid .le. eps ) then
            iflag = 0
            maxits = nits
            return
         end if
         if (its .ge. maxits*j) then
            iflag = 1
            maxits = nits
            return
         end if

C        print 11, '            ||res|| = ',resnor 
C 11     format(A,d)
C 13     format(i4,A,d)
      
      end do
C End of infinite iterative loop =================
C End of GMRESR subroutine      
      end 

C=============================================================================
       subroutine gmres0(oktest,n,im,rhs,uu,cc,work0,eps,maxits,matvec)

C This is the modified GMRES routine gmres0 adapted for GMRESR by 
C Mike Botchev, Utrecht University, Dec. 1996
C For detail on how to make GMRES (for GMRESR) cheaper see 
C the above-mentioned paper on GMRESR 
c*************************************************************
C This code was initially written by Youcef Saad
C then revised by Henk A. van der Vorst  
C and Mike Botchev (oct. 1996)
C ************************************************************ 
c gmres algorithm . simple version .  (may 23, 1985)
c parameter list:
c oktest == TRUE for printing intermediate results
c n      == size of problem
c im     == size of krylov subspace:  should not exceed 50 in this
c          version (can be reset in code. looking at comment below)
c rhs    == right hand side
c uu     == initial guess for vector u (see above-mentioned paper on GMRESR)
c           on input, approximate solution on output
c cc     == initial guess for vector c (see above-mentioned paper on GMRESR)
c           on input, approximate solution on output
c work0  == work space of size n x (im+1)
c eps    == tolerance for stopping criterion. process is stopped
c           as soon as ( ||.|| is the euclidean norm):
c           || current residual || <= eps  
c maxits == maximum number of iterations allowed
c           on OUTPUT: actual number of iterations performed
c ----------------------------------------------------------------
c subroutines 
c matvec      == matrix vector multiplication y <- A*x
c
c BLAS:
c dcopy       == y <-- x routine
c ddot        == dot product function
c dnrm2       == euclidean norm function
c daxpy       == y <-- y+ax routine
c dscal       == x <-- ax routine
c dtsrv       == to solve linear system with a triangular matrix
c*************************************************************
c-------------------------------------------------------------
c arnoldi size should not exceed 10 in this version..
c to reset modify maxdim. BUT:             ----------------
c maxdim was set to 10 because of two reasons:
c (1) it is assumed in this implementation that it is cheaper to
c make maxdim vector updates than to make 1 matrix-vector
c multiplication;
c (2) for large maxdim we may lose the orthogonality property
c on which this cheap implementation is based.
c Please keep it in mind changing maxdim
c-------------------------------------------------------------
      integer maxdim,maxd1,md1max
      parameter (maxdim=10, maxd1=maxdim+1, md1max=maxdim*maxd1)
      external matvec

      logical oktest
      integer jjj,jj1
      integer i,i1,im,its,j,k,k1,maxits,n
      real*8 cc,coeff,coef1,dabs,ddot,dnrm2,dsqrt,eps,epsmac,
     &                 gam,rhs(n),ro,uu(n),work0(n,im+1),t     

      real*8 hh(maxd1,maxdim),hh1(maxd1,maxdim),c(maxdim),
     &                 s(maxdim),rs(maxd1),rs1(maxd1)

      data (( hh(jj1,jjj), jj1=1,maxd1), jjj=1,maxdim) / md1max*0.0 / ,
     &      epsmac / 1.d-16 / 
C-----------------------------------------------------------------------------

      if (im .gt. maxdim) then
         im = maxdim
         write (*,'(A,i2)') 'GMRES0: dimension has been reduced to ',im
         write (*,'(A)') ' =&gt; reset MAXDIM if you want it to be more'
         write (*,'(A)') ' BUT read comments near MAXDIM before'
      end if

      its = 0

C     ----------------------------
C     Outer loop starts here.. 
C     BUT here (for GMRESR) only 1 outer loop is allowed
C     Compute initial residual vector 
C     ----------------------------
 10   continue
C        do not calculate initial residual first restart because 
C        initial guess is always zero. 
C        make initial guess zero:
         coeff = 0.0
         call dscal(n,coeff,uu,1)
C        make initial residual right-hand side:
         call dcopy(n,rhs,1,work0,1)

	 ro = dnrm2 (n, work0, 1)
	 if ((ro .eq. 0.0d0).or.(ro .le. eps)) then
            call matvec(uu, cc, n)
            eps = ro
            maxits = its 
            return
         end if

         coeff = 1 / ro
         call dscal(n,coeff,work0,1)

	 if (oktest) write(*, 199) its, ro

c        initialize 1-st term  of rhs of hessenberg system..
	 rs(1) = ro
	 i = 0

 4       continue
            i=i+1
            its = its + 1
            i1 = i + 1
            call  matvec(work0(1,i), work0(1,i1), n)
c           -----------------------------------------
c           modified gram - schmidt...
c           -----------------------------------------
            do j=1, i
               t = ddot(n, work0(1,j),1,work0(1,i1),1)
               hh(j,i) = t
               call daxpy(n, -t, work0(1,j), 1, work0(1,i1), 1)
            end do
            t = dnrm2(n, work0(1,i1), 1)
            hh(i1,i) = t
            if (t .ne. 0.0d0)then
               t = 1 / t
               call dscal(n, t, work0(1,i1), 1)
C              save new column of hh in hh1 to reproduce vector cc later on
               call dcopy(maxd1,hh(1,i),1,hh1(1,i),1)
            endif
c           done with modified gram schmidt and arnoldi step..

c           now  update factorization of hh
            if (i .ne. 1) then
c              perform previous transformations  on i-th column of h
               do k=2,i
                  k1 = k-1
                  t = hh(k1,i)
                  hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
                  hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
               end do
            endif
            gam = dsqrt(hh(i,i)**2 + hh(i1,i)**2)
            if (gam .eq. 0.0d0) gam = epsmac

c           determine next plane rotation
            c(i) = hh(i,i)/gam
            s(i) = hh(i1,i)/gam
            rs(i1) = -s(i)*rs(i)
            rs(i) =  c(i)*rs(i)

c           determine residual norm and test for convergence-
            hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
            ro = dabs(rs(i1))
            if (oktest) write(*, 199) its, ro
	 if ((i .lt. im) .and. (ro .gt. eps))  goto 4
c
c        now compute solution. first solve upper triangular system.
c
C        rs := hh(1:i,1:i) ^-1 * rs

         call dtrsv('U','N','N',i,hh,maxd1,rs,1)   
c        done with back substitution..

c        now form linear combination to get vector uu
	 do j=1, i
            t = rs(j)
	    call daxpy(n, t, work0(1,j), 1, uu,1)
         end do
C        DO NOT restart outer loop EVEN when necessary (that is crucial
C        for this implementation of GMRESR):  NEVER goto 10 !  
C     if (ro .gt. eps .and. its .lt. maxits) goto 10

C     Finally, reproduce vector cc as cc = A*uu = work0*hh1*rs:
C     rs := hh1(1:i1,1:i) * rs
      coeff = 1
      coef1 = 0
      call dgemv('N',i1,i,coeff,hh1,maxd1,rs,1,coef1,rs1,1)

C     now multiply Krylov basis vectors work0 by rs:
C     cc := work0*rs
      call dscal(n,coef1,cc,1)
      do j=1, i1
         t = rs1(j)
	 call daxpy(n, t, work0(1,j), 1, cc,1)
      end do        

 199  format('itsinn =', i4, ' res. norm =', d20.6)

      maxits=its
      eps=ro 
      return
c------------------------------- end of gmres0 ----------------------
      end
