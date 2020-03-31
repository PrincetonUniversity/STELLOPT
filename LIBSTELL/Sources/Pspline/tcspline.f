c  tcspline -- dmc 20 Jan 1999
c
c  set up coefficients for bicubic spline with following BC's:
c  * LHS and RHS handled as in cubspl.for for 1st coordinate
c  * derivatives periodic in second coordinate (use pspline.for)
c
c workspace:
c  if phi bdy cond. is periodic, not-a-knot, df/dphi = 0 everywhere,
c  or d2f/dphi2 = 0 everywhere, then the phi boundary condition is
c  "linear" and a workspace of size at least:
c
c     nwk = 20*inx*inth + 10*max(inx,inth,inph)
c
c  will suffice.
c
c  if the phi bdy cond. involves specification of df/dphi .ne. 0 or
c  d2f/dphi .ne. 0 at any (x,theta) grid point, then, the phi boundary
c  condition is "non-linear", a correction step is needed, and a workspace
c  of size at least:
c
c     nwk = 16*inx*inth*inph
c
c  is required.
c
      subroutine tcspline(x,inx,th,inth,ph,inph,fspl,inf4,inf5,
     >                    ibcxmin,bcxmin,ibcxmax,bcxmax,inb1x,
     >                    ibcthmin,bcthmin,ibcthmax,bcthmax,inb1th,
     >                    ibcphmin,bcphmin,ibcphmax,bcphmax,inb1ph,
     >                    wk,nwk,ilinx,ilinth,ilinph,ier)
c
      real x(inx),th(inth),ph(inph)
      real fspl(4,4,4,inf4,inf5,inph),wk(nwk)
      real bcxmin(inb1x,*),bcxmax(inb1x,*) ! inth x inph defined (if used)
      real bcthmin(inb1th,*),bcthmax(inb1th,*) ! inx x inph defined (if used)
      real bcphmin(inb1ph,*),bcphmax(inb1ph,*) ! inx x inth defined (if used)
c
c  input:
c    x(1...inx) -- abscissae, first dimension of data
c   th(1...inth) -- abscissae, second (periodic) dimension of data
c   ph(1...inph) -- abscissae, third (periodic) dimension of data
c   fspl(1,1,1,1..inx,1..inth,1..inph) -- function values
c   inf4 -- fspl dimensioning, inf4.ge.inx required.
c   inf5 -- fspl dimensioning, inf5.ge.inth required.
c
c  boundary conditions input:
c
c   bc data at xmin, xmax  vs.  theta,phi
c   bc data at thmin, thmax  vs.  x,phi
c   bc data at phmin, phmax  vs.  x,theta
c
c   ibcxmin -- indicator for boundary condition at x(1):
c    bcxmin(...) -- boundary condition data
c     =-1 -- use periodic boundary condition
c     =0 -- use "not a knot", bcxmin(...) ignored
c     =1 -- match slope, specified at x(1),th(ith),ph(iph) by bcxmin(ith,iph)
c     =2 -- match 2nd derivative, specified at x(1),th(ith),ph(iph)
c           by bcxmin(ith,iph
c     =3 -- boundary condition is slope=0 (df/dx=0) at x(1), all th(j)
c     =4 -- boundary condition is d2f/dx2=0 at x(1), all th(j)
c     =5 -- match 1st derivative to 1st divided difference
c     =6 -- match 2nd derivative to 2nd divided difference
c     =7 -- match 3rd derivative to 3rd divided difference
c           (for more detailed definition of BCs 5-7, see the
c           comments of subroutine mkspline)
c   NOTE bcxmin(...) referenced ONLY if ibcxmin=1 or ibcxmin=2
c
c   ibcxmax -- indicator for boundary condition at x(nx):
c    bcxmax(...) -- boundary condition data
c     (interpretation as with ibcxmin, bcxmin)
c     NOTE:  if ibcxmin=-1 then the periodic BC applies on both sides
c            and ibcxmax is ignored.
c   inb1x -- 1st dimension of bcxmin, bcxmax: if ibcxmin or ibcxmax .gt. 0
c            this must be .ge. inth:
c
c   interpretation of ibcthmin,bcthmin,ibcthmax,bcthmax,inb1th
c     is same as with ibcxmin,...
c
c   interpretation of ibcphmin,bcphmin,ibcphmax,bcphmax,inb1ph
c     is same as with ibcxmin,...
c
c   the explicit bdy condition arrays are referenced only if the
c     corresponding "ibc" flag values are set to 1 or 2.
c
c  output:
c   fspl(*,*,*,1..inx,1..inth,1..inph) -- bicubic spline coeffs (4x4)
c   ...fspl(1,1,1,*,*,*) is not replaced.
c
c   ilinx -- =1 on output if x(inx) pts are nearly evenly spaced (tol=1e-3)
c   ilinth-- =1 on output if th(inth) evenly spaced (tol=1e-3)
c   ilinph-- =1 on output if ph(inph) evenly spaced (tol=1e-3)
c
c   ier -- completion code, 0 for normal
c
c  workspace:
c   wk -- must be at least 5*max(inx,inth,inph) large -- or more, see
c         comments, above.
c  nwk -- size of workspace
c
c---------------------------------
c  ** in what follows, f is an abbreviation for fspl **
c
c  compute tricubic spline of 3d function, given values at the
c  grid crossing points, f(1,1,1,i,j,k)=f(x(i),th(j),ph(k)).
c
c  on evaluation:  for point x btw x(i) and x(i+1), dx=x-x(i)
c                       and th btw th(j) and th(j+1), dt=th-th(j),
c                       and ph btw ph(k) and ph(k+1), dp=ph-ph(k),
c
c      spline =
c        f(1,1,1,i,j,k)+dx*f(2,1,1,i,j,k)+dx2*f(3,1,1,i,j,k)+dx3*f(4,1,1,i,j,k)
c  +dt*(f(1,2,1,i,j,k)+dx*f(2,2,1,i,j,k)+dx2*f(3,2,1,i,j,k)+dx3*f(4,2,1,i,j,k))
c +dt2*(f(1,3,1,i,j,k)+dx*f(2,3,1,i,j,k)+dx2*f(3,3,1,i,j,k)+dx3*f(4,3,1,i,j,k))
c +dt3*(f(1,4,1,i,j,k)+dx*f(2,4,1,i,j,k)+dx2*f(3,4,1,i,j,k)+dx3*f(4,4,1,i,j,k))
c        +dp*(
c        f(1,1,2,i,j,k)+dx*f(2,1,2,i,j,k)+dx2*f(3,1,2,i,j,k)+dx3*f(4,1,2,i,j,k)
c  +dt*(f(1,2,2,i,j,k)+dx*f(2,2,2,i,j,k)+dx2*f(3,2,2,i,j,k)+dx3*f(4,2,2,i,j,k))
c +dt2*(f(1,3,2,i,j,k)+dx*f(2,3,2,i,j,k)+dx2*f(3,3,2,i,j,k)+dx3*f(4,3,2,i,j,k))
c +dt3*(f(1,4,2,i,j,k)+dx*f(2,4,2,i,j,k)+dx2*f(3,4,2,i,j,k)+dx3*f(4,4,2,i,j,k)))
c        +dp2*(
c        f(1,1,3,i,j,k)+dx*f(2,1,3,i,j,k)+dx2*f(3,1,3,i,j,k)+dx3*f(4,1,3,i,j,k)
c  +dt*(f(1,2,3,i,j,k)+dx*f(2,2,3,i,j,k)+dx2*f(3,2,3,i,j,k)+dx3*f(4,2,3,i,j,k))
c +dt2*(f(1,3,3,i,j,k)+dx*f(2,3,3,i,j,k)+dx2*f(3,3,3,i,j,k)+dx3*f(4,3,3,i,j,k))
c +dt3*(f(1,4,3,i,j,k)+dx*f(2,4,3,i,j,k)+dx2*f(3,4,3,i,j,k)+dx3*f(4,4,3,i,j,k)))
c        +dp3*(
c        f(1,1,4,i,j,k)+dx*f(2,1,4,i,j,k)+dx2*f(3,1,4,i,j,k)+dx3*f(4,1,4,i,j,k)
c  +dt*(f(1,2,4,i,j,k)+dx*f(2,2,4,i,j,k)+dx2*f(3,2,4,i,j,k)+dx3*f(4,2,4,i,j,k))
c +dt2*(f(1,3,4,i,j,k)+dx*f(2,3,4,i,j,k)+dx2*f(3,3,4,i,j,k)+dx3*f(4,3,4,i,j,k))
c +dt3*(f(1,4,4,i,j,k)+dx*f(2,4,4,i,j,k)+dx2*f(3,4,4,i,j,k)+dx3*f(4,4,4,i,j,k)))
c
c      where dx2=dx**2 and dx3=dx**3.
c      where dt2=dt**2 and dt3=dt**3.
c      where dp2=dp**2 and dp3=dp**3.
c
c---------------------------------
      integer iselect1(10)
      integer iselect2(10)
c
      real z0,z1,ztol
      real zcur(1)
c
      data z0/0.0e0/
      data z1/1.0e0/
      data ztol/1.0e-3/
c
c---------------------------------
c
      ier=0
c
      iflg=0
c
c  check phi bdy condition "linearity"
c
      if(ibcphmin.ne.-1) then
         if((ibcphmin.eq.1).or.(ibcphmin.eq.2)) then
            do ith=1,inth
               do ix=1,inx
                  if(bcphmin(ix,ith).ne.z0) iflg=1
               enddo
            enddo
         endif
         if((ibcphmax.eq.1).or.(ibcphmax.eq.2)) then
            do ith=1,inth
               do ix=1,inx
                  if(bcphmax(ix,ith).ne.z0) iflg=1
               enddo
            enddo
         endif
      endif
c
      itest=10*max(inx,inth,inph) + 20*inx*inth
      if(iflg.eq.1) then
         itest=16*inx*inth*inph
      endif
C
      if(nwk.lt.itest) then
         write(6,9901) nwk,itest
         ier=1
 9901    format(' ?tcspline:  workspace too small.'/
     >      '  user supplied nwk=',i7,'; need at least: ',i7/
     >      '  If no explicit df/dph boundary condition is set,'/
     >      '  nwk = 20*inx*inth + 10*max(inx,inth,inph) can be used.'/
     >      '  If an explicit df/dph or d2f/dph2 boundary condition'/
     >      '  is set, nwk=16*inx*inth*inph is required.')
      endif
      if(inx.lt.2) then
         write(6,'('' ?tcspline:  at least 2 x points required.'')')
         ier=1
      endif
      if(inth.lt.2) then
         write(6,'('' ?tcspline:  need at least 2 theta points.'')')
         ier=1
      endif
      if(inph.lt.2) then
         write(6,'('' ?tcspline:  need at least 2 phi points.'')')
         ier=1
      endif
c
      if((ibcxmin.eq.1).or.(ibcxmax.eq.1).or.(ibcxmin.eq.2).or.
     >   (ibcxmax.eq.2)) then
         if(inb1x.lt.inth) then
            ier=1
            write(6,
     >'('' ?tcspline:  1st dim of bcxmin/max arrays .lt. inth'')')
         endif
      endif
c
      if((ibcthmin.eq.1).or.(ibcthmax.eq.1).or.(ibcthmin.eq.2).or.
     >   (ibcthmax.eq.2)) then
         if(inb1th.lt.inx) then
            ier=1
            write(6,
     >'('' ?tcspline:  1st dim of bcthmin/max arrays .lt. inx'')')
         endif
      endif
c
      if((ibcphmin.eq.1).or.(ibcphmax.eq.1).or.(ibcphmin.eq.2).or.
     >   (ibcphmax.eq.2)) then
         if(inb1ph.lt.inx) then
            ier=1
            write(6,
     >'('' ?tcspline:  1st dim of bphmin/max arrays .lt. inx'')')
         endif
      endif
c
      call ibc_ck(ibcxmin,'tcspline','xmin',-1,7,ier)
      if(ibcxmin.ge.0) call ibc_ck(ibcxmax,'tcspline','xmax',0,7,ier)
c
      call ibc_ck(ibcthmin,'tcspline','thmin',-1,7,ier)
      if(ibcthmin.ge.0) call ibc_ck(ibcthmax,'tcspline','thmax',0,7,ier)
c
      call ibc_ck(ibcphmin,'tcspline','phmin',-1,7,ier)
      if(ibcphmax.ge.0) call ibc_ck(ibcphmax,'tcspline','phmax',0,7,ier)
c
c  check ilinx & x vector
c
      call splinck(x,inx,ilinx,ztol,ierx)
      if(ierx.ne.0) ier=2
c
      if(ier.eq.2) then
         write(6,'('' ?tcspline:  x axis not strict ascending'')')
      endif
c
c  check ilinth & th vector
c
      call splinck(th,inth,ilinth,ztol,ierth)
      if(ierth.ne.0) ier=3
c
      if(ier.eq.3) then
         write(6,'('' ?tcspline:  theta axis not strict ascending'')')
      endif
c
c  check ilinth & th vector
c
      call splinck(ph,inph,ilinph,ztol,ierph)
      if(ierph.ne.0) ier=4
c
      if(ier.eq.4) then
         write(6,'('' ?tcspline:  phi axis not strict ascending'')')
      endif
c
      if(ier.ne.0) return
c
c------------------------------------
c
c  part 1.  compute (x,theta) spline coeffs via an intermediate
c  routine that call bcspline
c
c  workspace addresses
c
      iaspl2=1
      iabcx1=iaspl2+16*inx*inth
      iabcx2=iabcx1+inth
      iabcth1=iabcx2+inth
      iabcth2=iabcth1+inx
      iawk=iabcth2+inx
      inwk=nwk-iawk+1
c
      do iph=1,inph
c
c  copy bc data
c
         do ix=1,inx
            wk(iabcth1+ix-1)=0.0
            wk(iabcth2+ix-1)=0.0
            if((ibcthmin.eq.1).or.(ibcthmin.eq.2)) then
               wk(iabcth1+ix-1)=bcthmin(ix,iph)
            endif
            if((ibcthmin.ne.-1).and.
     >         ((ibcthmax.eq.1).or.(ibcthmax.eq.2))) then
               wk(iabcth2+ix-1)=bcthmax(ix,iph)
            endif
         enddo
         do ith=1,inth
            wk(iabcx1+ith-1)=0.0
            wk(iabcx2+ith-1)=0.0
            if((ibcxmin.eq.1).or.(ibcxmin.eq.2)) then
               wk(iabcx1+ith-1)=bcxmin(ith,iph)
            endif
            if((ibcxmin.ne.-1).and.
     >         ((ibcxmax.eq.1).or.(ibcxmax.eq.2))) then
               wk(iabcx2+ith-1)=bcxmax(ith,iph)
            endif
         enddo
c
c  call 2d spline intermediary routine
c
         call tcsp23(x,inx,th,inth,fspl(1,1,1,1,1,iph),inf4,
     >      ibcxmin,wk(iabcx1),ibcxmax,wk(iabcx2),
     >      ibcthmin,wk(iabcth1),ibcthmax,wk(iabcth2),
     >      wk(iaspl2),wk(iawk),inwk,ilinx,ilinth,ier)
c
         if(ier.ne.0) then
            write(6,*) ' ?tcspline:  error in 2d spline, exiting.'
            return
         endif
c
      enddo
c
c  ok now fspl(*,*,1,*,*,*) have been evaluated and C2 in (x,theta)
c  now need to extend to coeffs in phi direction.
c
      xo2=0.5e0
      xo6=1.0e0/6.0e0
c
c  spline each (x,th) coeff in the phi direction
c
      inpho=4*(inph-1)
      do ith=1,inth-1
         do ix=1,inx-1
c
            do ic1=1,4
               do ic2=1,4
c
c  copy coeff. ordinates in
c
                  do iph=1,inph
                     wk(4*(iph-1)+1)=fspl(ic1,ic2,1,ix,ith,iph)
                  enddo
c
c  use linear BC on this first pass; will correct later if
c  necessary
c
                  wk(2)=0.0
                  wk(3)=0.0
                  wk(inpho+2)=0.0
                  wk(inpho+3)=0.0
c
                  ibcphmina=ibcphmin
                  ibcphmaxa=ibcphmax
                  if(iflg.eq.1) then
                     if((ibcphmin.eq.1).or.(ibcphmin.eq.2)) ibcphmina=0
                     if((ibcphmax.eq.1).or.(ibcphmax.eq.2)) ibcphmaxa=0
                  endif
c
                  call v_spline(ibcphmina,ibcphmaxa,inph,ph,wk,
     >               wk(4*inph+1))
c
c  copy coeffs out
c
                  do iph=1,inph-1
                     fspl(ic1,ic2,2,ix,ith,iph)=wk(4*(iph-1)+2)
                     fspl(ic1,ic2,3,ix,ith,iph)=wk(4*(iph-1)+3)*xo2
                     fspl(ic1,ic2,4,ix,ith,iph)=wk(4*(iph-1)+4)*xo6
                  enddo
c
               enddo                    ! ic2
            enddo                       ! ic1
c
         enddo                          ! ix
      enddo                             ! ith
c
c  if there are "non-linear" BCs requiring correction...
c
c  at each (x(ix),th(ith)) get the d/dph BC's right while preserving C2
c  everywhere...
c
      if(iflg.eq.1) then
c
c  first get BC correction numbers
c
         iabcph1=1
         iabcph2=iabcph1+inx*inth
         iaccoef=iabcph2+inx*inth
         iawk=iaccoef+12*inx*inth*inph
         inwk=nwk-iawk+1
c
         do i=1,10
            iselect1(i)=0
            iselect2(i)=0
         enddo
c
c  note because iflg=1, we know at least one of ibcphmin/max = 1 or 2
c
         iskip1=0
         if(ibcphmin.eq.1) then
            iselect1(4)=1               ! df/dph
         else if(ibcphmin.eq.2) then
            iselect1(7)=1               ! d2f/dph2
         else
            iskip1=1
         endif
c
         iskip2=0
         if(ibcphmax.eq.1) then
            iselect2(4)=1               ! df/dph
         else if(ibcphmax.eq.2) then
            iselect2(7)=1               ! d2f/dph2
         else
            iskip2=1
         endif
c
         ia1=iabcph1-1
         ia2=iabcph2-1
         do ith=1,inth
            do ix=1,inx
               ia1=ia1+1
               ia2=ia2+1
c
               if(iskip1.eq.0) then
                  call tcspeval(x(ix),th(ith),ph(1),iselect1, zcur,
     >               x,inx,th,inth,ph,inph,ilinx,ilinth,ilinph,
     >               fspl,inf4,inf5,ier)
                  if(ier.ne.0) then
                     write(6,*) ' ?? tcspline:  error in tcspeval call'
                     return
                  endif
                  wk(ia1)=bcphmin(ix,ith)-zcur(1) ! correction needed
               else
                  wk(ia1)=z0
               endif
c
               if(iskip2.eq.0) then
                  call tcspeval(x(ix),th(ith),ph(inph),iselect2, zcur,
     >               x,inx,th,inth,ph,inph,ilinx,ilinth,ilinph,
     >               fspl,inf4,inf5,ier)
                  if(ier.ne.0) then
                     write(6,*) ' ?? tcspline:  error in tcspeval call'
                     return
                  endif
                  wk(ia2)=bcphmax(ix,ith)-zcur(1) ! correction needed
               else
                  wk(ia2)=z0
               endif
            enddo
         enddo
c
         call tcspcorr(x,inx,th,inth,ph,inph,fspl,inf4,inf5,
     >      ibcxmin,ibcxmax,ibcthmin,ibcthmax,
     >      ibcphmin,wk(iabcph1),ibcphmax,wk(iabcph2),
     >      wk(iaccoef),wk(iawk),inwk)
c
      endif
c
      return
      end
c-----------------------------------------
      subroutine tcspcorr(x,inx,th,inth,ph,inph,fspl,inf4,inf5,
     >   ibcxmin,ibcxmax,ibcthmin,ibcthmax,
     >   ibcphmin,bcph1,ibcphmax,bcph2,ccorr,wk,nwk)
c
c  intermediary routine for tcspline:
c  do correction needed to get C2 3d spline with phi bdy conditions
c  matched.
c
c  all input unless noted:
c
      real x(inx)                       ! x axis
      real th(inth)                     ! th axis
      real ph(inph)                     ! ph axis
c
      real fspl(4,4,4,inf4,inf5,inph)   ! spline coeffs -- adjusted
c
      integer ibcxmin,ibcxmax           ! x BC flags
      integer ibcthmin,ibcthmax         ! th BC flags
      integer ibcphmin,ibcphmax         ! ph BC flags
c
      real bcph1(inx,inth)              ! ph BC correction array @ phi(1)
      real bcph2(inx,inth)              ! ph BC correction array @ phi(inph)
c
c  workspaces:
c
      real ccorr(3,4,inx,inth,inph)     ! correction coefficients (partial)
c
      real wk(nwk)                      ! workspace
c
c---------------------------
c
      xo2=0.5e0
      xo6=1.0e0/6.0e0
c
c  1.  splines in phi -- fcns are zero everywhere but have non-zero BCs
c
      if(nwk.lt.10*max(inx,inth,inph)) then
         write(6,*) ' ?? programming error in tcspcorr (tcspline)'
         return
      endif
c
      z0=0.0
      iawk2=4*inph+1
c
      inpho=4*(inph-1)
      do ith=1,inth
         do ix=1,inx
c
            do iph=1,inph
               wk(4*(iph-1)+1)=z0
            enddo
c
c  set BC for this 1d spline correction
c
            if(ibcphmin.eq.1) then
               wk(2)=bcph1(ix,ith)
            else if(ibcphmin.eq.2) then
               wk(3)=bcph1(ix,ith)
            endif
c
            if(ibcphmax.eq.1) then
               wk(inpho+2)=bcph2(ix,ith)
            else if(ibcphmax.eq.2) then
               wk(inpho+3)=bcph2(ix,ith)
            endif
c
            call v_spline(ibcphmin,ibcphmax,inph,ph,wk,wk(iawk2))
c
c  copy non-zero coeffs out to ccorr
c
            do iph=1,inph-1
               ccorr(1,1,ix,ith,iph)=wk(4*(iph-1)+2)
               ccorr(2,1,ix,ith,iph)=wk(4*(iph-1)+3)*xo2
               ccorr(3,1,ix,ith,iph)=wk(4*(iph-1)+4)*xo6
            enddo
c
         enddo
      enddo
c
c  2. spline the coeffs in x -- use ibcx flags & zero for derivative
c  bc if necessary
c
      iawk2=4*inx+1
c
      inxo=4*(inx-1)
      do iph=1,inph-1
         do ith=1,inth
c
            do icph=1,3
c
               do ix=1,inx
                  wk(4*(ix-1)+1)=ccorr(icph,1,ix,ith,iph)
               enddo
c
c  zero BC:  correction spline
c
               wk(2)=0.0
               wk(3)=0.0
               wk(inxo+2)=0.0
               wk(inxo+3)=0.0
c
               call v_spline(ibcxmin,ibcxmax,inx,x,wk,wk(iawk2))
c
               do ix=1,inx-1
                  ccorr(icph,2,ix,ith,iph)=wk(4*(ix-1)+2)
                  ccorr(icph,3,ix,ith,iph)=wk(4*(ix-1)+3)*xo2
                  ccorr(icph,4,ix,ith,iph)=wk(4*(ix-1)+4)*xo6
               enddo
c
            enddo
c
         enddo
      enddo
c
c  3.  spline all the ccorr coefs in th -- use ibcth flags & zero for
c      derivative correction BC if necessary
c
c      add the results into fspl
c
      iawk2=4*inth+1
c
      intho=4*(inth-1)
      do iph=1,inph-1
         do ix=1,inx-1
c
            do icx=1,4
               do icph=1,3
c
                  do ith=1,inth
                     wk(4*(ith-1)+1)=ccorr(icph,icx,ix,ith,iph)
                  enddo
c
c  zero BC:  correction spline
c
                  wk(2)=0.0
                  wk(3)=0.0
                  wk(intho+2)=0.0
                  wk(intho+3)=0.0
c
                  call v_spline(ibcthmin,ibcthmax,inth,th,wk,
     >               wk(iawk2))
c
                  do ith=1,inth-1
                     do icth=1,4
                        zfac=1.0
                        if(icth.eq.3) zfac=xo2
                        if(icth.eq.4) zfac=xo6
                        fspl(icx,icth,icph+1,ix,ith,iph)=
     >                     fspl(icx,icth,icph+1,ix,ith,iph)+
     >                     wk(4*(ith-1)+icth)*zfac
                     enddo
                  enddo
c
               enddo
            enddo
c
         enddo
      enddo
c
      return
      end
c-----------------------------------------
      subroutine tcsp23(x,inx,th,inth,fspl,inf4,
     >   ibcxmin,bcxmin,ibcxmax,bcxmax,
     >   ibcthmin,bcthmin,ibcthmax,bcthmax,
     >   fspl2,wk,nwk,ilinx,ilinth,ier)
c
c  intermediary routines
c  call bcspline from tcspline loop
c  to set up 2d splines in each phi plane
c
      real x(inx)                       ! x axis
      real th(inth)                     ! th axis
c
      real fspl(4,4,4,inf4,inth)        ! fspl array at one phi pt.
      real fspl2(4,4,inx,inth)          ! temp fspl array for bcspline
c
      real bcxmin(inth),bcxmax(inth)    ! d/dx BC's @ x(1),x(inx), th(*)
      real bcthmin(inx),bcthmax(inx)    ! d/dth BC's @ th(1),th(inth), x(*)
c
      real wk(nwk)
c
c--------------------
c
c  1.  copy spline data in
c
      do ith=1,inth
         do ix=1,inx
            fspl2(1,1,ix,ith)=fspl(1,1,1,ix,ith)
         enddo
      enddo
c
c  2.  compute the 2d spline
c
      call bcspline(x,inx,th,inth,fspl2,inx,
     >   ibcxmin,bcxmin,ibcxmax,bcxmax,
     >   ibcthmin,bcthmin,ibcthmax,bcthmax,
     >   wk,nwk,ilinx,ilinth,ier)
      if(ier.ne.0) return
c
c  3.  copy spline coeff results out
c
      do ith=1,inth-1
         do ix=1,inx-1
            do j=1,4
               do i=1,4
                  fspl(i,j,1,ix,ith)=fspl2(i,j,ix,ith)
               enddo
            enddo
         enddo
      enddo
c
      return
      end
