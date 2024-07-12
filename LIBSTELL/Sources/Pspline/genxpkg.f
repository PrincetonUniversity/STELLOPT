      subroutine genxpkg(nx,x,xpkg,iper,imsg,itol,ztol,ialg,ier)
c
c  from an x axis assemble a "package":
c
c  MOD DMC Feb 2010: handle very small grids: nx=2, nx=3...
c  interchanged meaning of xpkg(1,4) and xpkg(3,4); 
c  xpkg(3,4) only set if nx.ge.3;
c  xpkg(4,4) only set if nx.ge.4.
c
c  there are corresponding changes in xlookup: simplified code lookup
c  code for the cases nx=2 and nx=3.
c
c     xpkg(j,1) = x(j), j = 1 to nx    ... nx.ge.2
c
c     if(abs(ialg).ne.3) then...
c       for j=1:nx-1
c       xpkg(j,2) = h(j) = x(j+1)-x(j), j=1 to nx-1
c     else
c       for j=1:nx-1
c       xpkg(j,2) = index location, with linear offset, of
c                   xpkg(1,1)+<h>*(j-1) in the original x(1:nx)
c            (piecewise linear indexing function)
c     endif
c     xpkg(nx,2) = <h> = (x(nx)-x(1))/(nx-1)
c
c     xpkg(j,3) = 1/h(j)
c     xpkg(nx,3) = 1/<h>
c
c     xpkg(1,4) = +/- tolerance epsilon for out-of-range warnings/errors
c                 +if message is to be written on out of range errors
c                 -if no message to be written.  In either case, a
c                        warning flag is set.
c
c     xpkg(2,4) = 1.0 if x is *periodic*, else 0.0
c
c  only set if nx.ge.3:
c     xpkg(3,4) = 0.0 if x is *evenly spaced* (to ztol or 5.e-7 rel) else:
c                   = 1.0 => use (1/h) newton's method like search algorithm
c                   = 2.0 => use binary search algorithm
c                   = 3.0 => use piecewise-linear indexing function...
c
c  only set if nx.ge.4:
c     xpkg(4,4) = 0.0 -- do not use past search result as start point
c                        for next search;
c               = 1.0 -- do use past search result as start point for
c                        next search.
c
c  tolerance epsilon means:
c     if xtest is within epsilon*max(abs(x(1)),abs(x(nx))) outside the
c     range [x(1),x(nx)] treat xtest as being at the endpoint x(1) or
c     x(nx) whichever is closer.
c
c input arguments:
c
      integer nx                        ! length of x, .ge.4
      real x(nx)                        ! x axis vector, strict ascending
c
      integer iper                      ! =1 if x is periodic
      integer imsg                      ! =1 for range error messages
      integer itol                      ! =1 to specify tolerance, else default
c
c default tolerance is 5.0e-7
c
      real ztol                         ! range tolerance, if itol=1
c
c lookup algorithm selection for uneven grids:
c
      integer ialg                      ! = +/- 1:  <1/h(j)> "Newton" method
c                                       ! = +/- 2:  binary search
c                                       ! = +/- 3:  indexing function
c
c       to use past search result as init. cond. for next search, give
c       ialg .lt. 0; otherwise give ialg .gt. 0.
c
c output arguments:
c
      real xpkg(nx,4)                   ! xpkg, as described above
      integer ier                       ! completion code, 0=OK
c
c------------------------------------------------
c
      if(nx.lt.2) then
         write(6,*) ' %genxpkg:  nx.ge.2 required!'
         ier=1
         return
      else
         ier=0
      endif
c
      ialgu=ialg
      if(ialgu.eq.0) ialgu=3
      if(iabs(ialgu).gt.3) ialgu=3
c
c  get tolerance for end point range check & even spacing check
c
      if(itol.eq.1) then
         ztolr=abs(ztol)
      else
         ztolr=5.0e-7
      endif
c
      ztola=max(abs(x(1)),abs(x(nx)))*ztolr
c
c  assume even spacing for now...
c
      if(nx.ge.3) then
         xpkg(3,4)=0.0
      endif
c
c  mark if x axis is a periodic coordinate
c
      if(iper.eq.1) then
         xpkg(2,4)=1.0
      else
         xpkg(2,4)=0.0
      endif
c
c  store tolerance parameter
c
      xpkg(1,4)=ztola
c
c  mark if messages are to be written if range lookup errors occur
c
      if(imsg.eq.1) then
         continue                       ! xpkg(1,4) left .gt. 0.0
      else
         xpkg(1,4)=-xpkg(1,4)
      endif
c
c  OK check linearity and spacing
c
      ier=0
c
      xpkg(nx,2)=(x(nx)-x(1))/(nx-1)    ! average spacing
c
      do ix=1,nx
         xpkg(ix,1)=x(ix)
         if((ier.eq.0).and.(ix.lt.nx)) then
            if(x(ix+1).le.x(ix)) then
               ier=1
               write(6,*) ' %genxpkg:  x axis not strict ascending!'
            else
               zh=x(ix+1)-x(ix)
c
c  checking even spacing now...
c
               if(nx.ge.3) then
                  if(abs(zh-xpkg(nx,2)).gt.ztola) xpkg(3,4)=1.0
               endif
c
               xpkg(ix,2)=zh
               xpkg(ix,3)=1.0/zh
c
            endif
         endif
      enddo
c
      if(ier.ne.0) return
c
c  OK, store inverse average spacing
c
      xpkg(nx,3)=1.0/xpkg(nx,2)
c
c  if even spacing is detected, redefine x axis slightly, for
c  improved regularity of behaviour
c
      if(nx.ge.3) then
c
         if(xpkg(3,4).eq.0.0) then
            do ix=1,nx-2
               ixp=ix+1
               xpkg(ixp,1)=xpkg(1,1)+ix*xpkg(nx,2)
            enddo
            if(nx.gt.3) then
               if(ialgu.lt.0) then
                  xpkg(4,4)=1.0 ! check init guess
               else
                  xpkg(4,4)=0.0
               endif
            endif
         endif
c
c  if uneven spacing is detected, must use an algorithm...
c
         if(xpkg(3,4).ne.0.0) then
            xpkg(3,4)=abs(ialgu)
            if(nx.gt.3) then
               if(ialgu.lt.0) then
                  xpkg(4,4)=1.0 ! check init guess
               else
                  xpkg(4,4)=0.0
               endif
            endif
c
            if(abs(ialgu).eq.3) then
c
c  construct a piecewise linear indexing function
c
               xpkg(1,2)=1.0
               xtest=xpkg(1,1)
               itest=1
               do i=2,nx-1
                  xtest=xtest+xpkg(nx,2) ! x1 + (i-1)*<h>
 10               continue
                  if((xpkg(itest,1).le.xtest).and.
     >                 (xtest.le.xpkg(itest+1,1))) then
                     xpkg(i,2)=itest+(xtest-xpkg(itest,1))/
     >                    (xpkg(itest+1,1)-xpkg(itest,1))
                  else
                     itest=itest+1
                     go to 10
                  endif
               enddo
c
c  (implicitly, xpkg(nx,2)=nx*1.0; but we leave <h> in xpkg(nx,2)
c   instead...)
c
            endif
         endif
c
      endif  ! nx.ge.3
c
      return
      end
