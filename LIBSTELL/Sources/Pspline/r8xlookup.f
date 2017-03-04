      subroutine r8xlookup(ivec,xvec,nx,xpkg,imode,iv,dxn,hv,hiv,iwarn)
c
c  vector lookup routine
c
c   given a set of x points xvec(...) and an x grid (xpkg, nx x pts)
c   return the vector of indices iv(...) and displacements dxv(...)
c   within the indexed zones, corresponding to each x point.
c
c   if any of the x points in the vector are out of range the warning
c   flag iwarn is set.
c
c  MOD DMC Feb 2010: changes related to supporting nx=2 and nx=3 small grids.
c  Changes are consistent with Feb 2010 changes in genxpkg:
c    meanings of xpkg(1,4) and xpkg(3,4) interchanged.
c
c    if nx.eq.2:  xpkg(3,4) and xpkg(4,4) never referenced;
c    if nx.eq.3:  xpkg(4,4) never referenced.
c
c----------------------------
c  input:
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER inum,istat,ilin,ialg,iper,imsg,init_guess,iprev,i
      INTEGER init,iprob,isrch,i_sign,inc
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 stat,ztola,period,hav,havi,hloci,hloc,zdelta,xfac
      REAL*8 zindx0,zdindx,zindex
!============
      integer ivec                      ! size of vector
      REAL*8 xvec(ivec)                   ! x points to lookup on xpkg grid
c
      integer nx                        ! size of grid
      REAL*8 xpkg(nx,4)                   ! grid data
c
      integer imode                     ! output control flag
c
c  imode=1:  return indices iv(...) and un-normalized displacements dxn(...)
c            ignore hv and hiv
c                 dxn(j)=xvec(j)-xpkg(iv(j),1)
c
c  imode=2   return indices iv(...) and *normalized* displacements dxn(...)
c            and hv(...) and hiv(...)
c                 dxn(j)=(xvec(j)-xpkg(iv(j),1))*hiv(j)
c
c  output:
      integer iv(ivec)                  ! index into grid for each xvec(j)
      !  note: old values of iv(...) may be used as start point for grid
      !  searches, depending on xpkg controls
 
      REAL*8 dxn(ivec)                    ! displacement w/in zone (see imode)
      REAL*8 hv(ivec)                     ! zone width (if imode=2)
      REAL*8 hiv(ivec)                    ! inverse zone width (if imode=2)
c
      integer iwarn                     ! =0: OK; =n:  n points out of range
c
c----------------------------
c
c  xpkg is a "structure" constructed by subroutine genxpkg, which
c  contains the x grid, xpkg(1:nx,1), spacing and error handling
c  information -- see genxpkg.for (& r8genxpkg.for)
c
      REAL*8, dimension(:), allocatable :: xuse
      logical, dimension(:), allocatable :: iok
      integer, dimension(:), allocatable :: imina,imaxa
c
cdbg      data idbg/0/
c
c----------------------------
      if(nx.lt.2) then
         iwarn=1
         write(6,*) ' ?? xlookup: nx.lt.2, nx=',nx
         go to 1100
      endif
 
      inum=ivec
      allocate(xuse(inum),stat=istat)
      if(istat.ne.0) then
         iwarn=1
         write(6,*) ' ?? xlookup "xuse" vector allocation failure!'
         go to 1000
      endif
c
      if(nx.eq.2) then
         ilin=1
      endif
c
      if(nx.gt.2) then
         if(xpkg(3,4).eq.0.0_r8) then
            ilin=1              ! evenly spaced grid
         else
            ilin=0
            if(xpkg(3,4).gt.2.5_r8) then
               ialg=3
            else if(xpkg(3,4).gt.1.5_r8) then
               ialg=2
            else
               ialg=1
            endif
         endif
      endif
c
      if(xpkg(2,4).ne.0.0_r8) then
         iper=1                         ! periodic grid
      else
         iper=0
      endif
c
      ztola=abs(xpkg(1,4))              ! tolerance for range checking
c
      if(xpkg(1,4).ge.0.0_r8) then
         imsg=1                         ! write message on range error
      else
         imsg=0
      endif
c
      init_guess=0
      if(nx.gt.3) then
         if(xpkg(4,4).gt.0.0_r8) then
            init_guess=1
            iprev=min((nx-1),max(1,iv(ivec)))
         else
            init_guess=0
         endif
      endif
c
      iwarn=0
c
c---------------------
c  range check
c
      if(iper.eq.0) then
c
c  check min/max with tolerance
c
         do i=1,ivec
            if(xvec(i).lt.xpkg(1,1)) then
               xuse(i)=xpkg(1,1)
               if((xpkg(1,1)-xvec(i)).gt.ztola) iwarn=iwarn+1
            else if(xvec(i).gt.xpkg(nx,1)) then
               xuse(i)=xpkg(nx,1)
               if((xvec(i)-xpkg(nx,1)).gt.ztola) iwarn=iwarn+1
            else
               xuse(i)=xvec(i)
            endif
         enddo
c
      else
c
c  normalize to interval
c
         period=xpkg(nx,1)-xpkg(1,1)
         do i=1,ivec
            if((xvec(i).lt.xpkg(1,1)).or.(xvec(i).gt.xpkg(nx,1))) then
               xuse(i)=mod(xvec(i)-xpkg(1,1),period)
               if(xuse(i).lt.0.0_r8) xuse(i)=xuse(i)+period
               xuse(i)=xuse(i)+xpkg(1,1)
               xuse(i)=max(xpkg(1,1),min(xpkg(nx,1),xuse(i)))
            else
               xuse(i)=xvec(i)
            endif
         enddo
c
      endif
c
      if((imsg.eq.1).and.(iwarn.gt.0)) then
         write(6,*) ' %xlookup:  ',iwarn,' points not in range: ',
     >      xpkg(1,1),' to ',xpkg(nx,1)
      endif
c
c---------------------
c  actual lookup -- initially, assume even spacing
c
c   initial index guess:  1 + <1/h>*(x-x0) ... then refine
c
      hav=xpkg(nx,2)
      havi=xpkg(nx,3)
      if(ilin.eq.1) then
c
c  faster lookup OK:  even spacing
c
         if(init_guess.eq.0) then
c
c  even spacing lookup, no initial guess from previous iteration
c
            do i=1,ivec
               iv(i)=1+havi*(xuse(i)-xpkg(1,1))
               iv(i)=max(1,min((nx-1),iv(i)))
               if(imode.eq.1) then
                  dxn(i)=(xuse(i)-xpkg(iv(i),1))
               else
                  dxn(i)=(xuse(i)-xpkg(iv(i),1))*havi
                  hiv(i)=havi
                  hv(i)=hav
               endif
            enddo
c
         else
c
c  even spacing lookup, do use initial guess from previous iteration
c
            do i=1,ivec
               if((xpkg(iprev,1).le.xuse(i)).and.
     >            (xuse(i).le.xpkg(iprev+1,1))) then
                  iv(i)=iprev
               else
                  iv(i)=1+havi*(xuse(i)-xpkg(1,1))
                  iv(i)=max(1,min((nx-1),iv(i)))
                  iprev=iv(i)
               endif
               if(imode.eq.1) then
                  dxn(i)=(xuse(i)-xpkg(iv(i),1))
               else
                  dxn(i)=(xuse(i)-xpkg(iv(i),1))*havi
                  hiv(i)=havi
                  hv(i)=hav
               endif
            enddo
c
         endif
c
         go to 1000
c
      endif
c
 4    continue
c
      init=-1
      allocate(iok(inum),stat=istat)
      if(istat.ne.0) then
         iwarn=1
         write(6,*) ' ?? xlookup "iok" vector allocation failure!'
         go to 1000
      endif
c
      if(ialg.lt.3) then
         allocate(imina(inum),stat=istat)
         if(istat.ne.0) then
            iwarn=1
            write(6,*) ' ?? xlookup "imina" vector allocation failure!'
            go to 1000
         endif
c
         allocate(imaxa(inum),stat=istat)
         if(istat.ne.0) then
            iwarn=1
            write(6,*) ' ?? xlookup "imaxa" vector allocation failure!'
            go to 1000
         endif
      endif
c
 5    continue                          ! re-entry -- hit problem cases...
c
      iprob=0                           ! count "problem" cases...
      init=init+1
c
      if(init_guess.eq.0) then
         go to (100,200,300) ialg
      else
         go to (150,250,350) ialg
      endif
c-----------------------------------------------------------
c  Newton like algorithm:  use local spacing to estimate
c   step to next zone guess; don't use prior guess
c
 100  continue
      do i=1,ivec
cdbg         jdbg=0
c
c  first iteration
c
         if(init.eq.0) then
            iok(i)=.FALSE.
            imina(i)=1
            imaxa(i)=nx-1
c
            iv(i)=1+havi*(xuse(i)-xpkg(1,1))
cdbg            jdbg=jdbg+1
cdbg            if(jdbg.le.idbg) iv(i)=1
            iv(i)=max(imina(i),min(imaxa(i),iv(i)))
c
            if(xuse(i).lt.xpkg(iv(i),1)) then
               isrch=0
               i_sign=-1
               imaxa(i)=max(1,(iv(i)-1))
            else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
               isrch=1
               i_sign=+1
               imina(i)=min((nx-1),(iv(i)+1))
            else
               iok(i)=.TRUE.
            endif
c
         endif
c
c  second iteration
c
         if(.not.iok(i)) then
            hloci=xpkg(iv(i),3)
            hloc=xpkg(iv(i),2)
            zdelta=(xuse(i)-xpkg(iv(i)+isrch,1))
            if(i_sign*zdelta.le.hloc) then
               inc=i_sign
            else
               inc=zdelta*hloci
               inc=inc+i_sign
            endif
c
            iv(i)=iv(i)+inc
cdbg            jdbg=jdbg+1
cdbg            if(jdbg.le.idbg) iv(i)=1
            iv(i)=max(imina(i),min(imaxa(i),iv(i)))
c
            if(xuse(i).lt.xpkg(iv(i),1)) then
               isrch=0
               i_sign=-1
               imaxa(i)=max(1,(iv(i)-1))
            else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
               isrch=1
               i_sign=+1
               imina(i)=min((nx-1),(iv(i)+1))
            else
               iok(i)=.TRUE.
            endif
c
c  third iteration
c
            if(.not.iok(i)) then
               hloci=xpkg(iv(i),3)
               hloc=xpkg(iv(i),2)
               zdelta=(xuse(i)-xpkg(iv(i)+isrch,1))
               if(i_sign*zdelta.le.hloc) then
                  inc=i_sign
               else
                  inc=zdelta*hloci
                  inc=inc+i_sign
               endif
c
               iv(i)=iv(i)+inc
cdbg            jdbg=jdbg+1
cdbg            if(jdbg.le.idbg) iv(i)=1
               iv(i)=max(imina(i),min(imaxa(i),iv(i)))
c
               if(xuse(i).lt.xpkg(iv(i),1)) then
                  isrch=0
                  i_sign=-1
                  imaxa(i)=max(1,(iv(i)-1))
               else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                  isrch=1
                  i_sign=+1
                  imina(i)=min((nx-1),(iv(i)+1))
               else
                  iok(i)=.TRUE.
               endif
c
c  fourth iteration
c
               if(.not.iok(i)) then
                  hloci=xpkg(iv(i),3)
                  hloc=xpkg(iv(i),2)
                  zdelta=(xuse(i)-xpkg(iv(i)+isrch,1))
                  if(i_sign*zdelta.le.hloc) then
                     inc=i_sign
                  else
                     inc=zdelta*hloci
                     inc=inc+i_sign
                  endif
c
                  iv(i)=iv(i)+inc
cdbg            jdbg=jdbg+1
cdbg            if(jdbg.le.idbg) iv(i)=1
                  iv(i)=max(imina(i),min(imaxa(i),iv(i)))
c
                  if(xuse(i).lt.xpkg(iv(i),1)) then
                     isrch=0
                     i_sign=-1
                     imaxa(i)=max(1,(iv(i)-1))
                  else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                     isrch=1
                     i_sign=+1
                     imina(i)=min((nx-1),(iv(i)+1))
                  else
                     iok(i)=.TRUE.
                  endif
c
c  fifth iteration
c
                  if(.not.iok(i)) then
                     hloci=xpkg(iv(i),3)
                     hloc=xpkg(iv(i),2)
                     zdelta=(xuse(i)-xpkg(iv(i)+isrch,1))
                     if(i_sign*zdelta.le.hloc) then
                        inc=i_sign
                     else
                        inc=zdelta*hloci
                        inc=inc+i_sign
                     endif
c
                     iv(i)=iv(i)+inc
cdbg            jdbg=jdbg+1
cdbg            if(jdbg.le.idbg) iv(i)=1
                     iv(i)=max(imina(i),min(imaxa(i),iv(i)))
c
                     if(xuse(i).lt.xpkg(iv(i),1)) then
                        isrch=0
                        i_sign=-1
                        imaxa(i)=max(1,(iv(i)-1))
                     else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                        isrch=1
                        i_sign=+1
                        imina(i)=min((nx-1),(iv(i)+1))
                     else
                        iok(i)=.TRUE.
                     endif
c
                     if(.not.iok(i)) iprob=iprob+1
c
c  end chain of iteration if-then-else blocks
c
                  endif
               endif
            endif
         endif
c
c  end of loop
c
      enddo
c
      go to 500
c-----------------------------------------------------------
c  Newton like algorithm:  use local spacing to estimate
c   step to next zone guess; DO use prior guess
c
 150  continue
c
      do i=1,ivec
c
c  first iteration
c
         if(init.eq.0) then
            if((xpkg(iprev,1).le.xuse(i)).and.
     >         (xuse(i).le.xpkg(iprev+1,1))) then
               iok(i)=.TRUE.
               iv(i)=iprev
            else
               iok(i)=.FALSE.
               imina(i)=1
               imaxa(i)=nx-1
c
               iv(i)=1+havi*(xuse(i)-xpkg(1,1))
               iv(i)=max(imina(i),min(imaxa(i),iv(i)))
c
               if(xuse(i).lt.xpkg(iv(i),1)) then
                  isrch=0
                  i_sign=-1
                  imaxa(i)=max(1,(iv(i)-1))
               else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                  isrch=1
                  i_sign=+1
                  imina(i)=min((nx-1),(iv(i)+1))
               else
                  iok(i)=.TRUE.
               endif
            endif
c
         endif
c
c  second iteration
c
         if(.not.iok(i)) then
            hloci=xpkg(iv(i),3)
            hloc=xpkg(iv(i),2)
            zdelta=(xuse(i)-xpkg(iv(i)+isrch,1))
            if(i_sign*zdelta.le.hloc) then
               inc=i_sign
            else
               inc=zdelta*hloci
               inc=inc+i_sign
            endif
c
            iv(i)=iv(i)+inc
            iv(i)=max(imina(i),min(imaxa(i),iv(i)))
c
            if(xuse(i).lt.xpkg(iv(i),1)) then
               isrch=0
               i_sign=-1
               imaxa(i)=max(1,(iv(i)-1))
            else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
               isrch=1
               i_sign=+1
               imina(i)=min((nx-1),(iv(i)+1))
            else
               iok(i)=.TRUE.
            endif
c
c  third iteration
c
            if(.not.iok(i)) then
               hloci=xpkg(iv(i),3)
               hloc=xpkg(iv(i),2)
               zdelta=(xuse(i)-xpkg(iv(i)+isrch,1))
               if(i_sign*zdelta.le.hloc) then
                  inc=i_sign
               else
                  inc=zdelta*hloci
                  inc=inc+i_sign
               endif
c
               iv(i)=iv(i)+inc
               iv(i)=max(imina(i),min(imaxa(i),iv(i)))
c
               if(xuse(i).lt.xpkg(iv(i),1)) then
                  isrch=0
                  i_sign=-1
                  imaxa(i)=max(1,(iv(i)-1))
               else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                  isrch=1
                  i_sign=+1
                  imina(i)=min((nx-1),(iv(i)+1))
               else
                  iok(i)=.TRUE.
               endif
c
c  fourth iteration
c
               if(.not.iok(i)) then
                  hloci=xpkg(iv(i),3)
                  hloc=xpkg(iv(i),2)
                  zdelta=(xuse(i)-xpkg(iv(i)+isrch,1))
                  if(i_sign*zdelta.le.hloc) then
                     inc=i_sign
                  else
                     inc=zdelta*hloci
                     inc=inc+i_sign
                  endif
c
                  iv(i)=iv(i)+inc
                  iv(i)=max(imina(i),min(imaxa(i),iv(i)))
c
                  if(xuse(i).lt.xpkg(iv(i),1)) then
                     isrch=0
                     i_sign=-1
                     imaxa(i)=max(1,(iv(i)-1))
                  else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                     isrch=1
                     i_sign=+1
                     imina(i)=min((nx-1),(iv(i)+1))
                  else
                     iok(i)=.TRUE.
                  endif
c
c  fifth iteration
c
                  if(.not.iok(i)) then
                     hloci=xpkg(iv(i),3)
                     hloc=xpkg(iv(i),2)
                     zdelta=(xuse(i)-xpkg(iv(i)+isrch,1))
                     if(i_sign*zdelta.le.hloc) then
                        inc=i_sign
                     else
                        inc=zdelta*hloci
                        inc=inc+i_sign
                     endif
c
                     iv(i)=iv(i)+inc
                     iv(i)=max(imina(i),min(imaxa(i),iv(i)))
c
                     if(xuse(i).lt.xpkg(iv(i),1)) then
                        isrch=0
                        i_sign=-1
                        imaxa(i)=max(1,(iv(i)-1))
                     else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                        isrch=1
                        i_sign=+1
                        imina(i)=min((nx-1),(iv(i)+1))
                     else
                        iok(i)=.TRUE.
                     endif
c
                     if(.not.iok(i)) iprob=iprob+1
c
c  end chain of iteration if-then-else blocks
c
                  endif
               endif
            endif
         endif
c
c  end of loop
c
         iprev=iv(i)
      enddo
c
      go to 500
c-----------------------------------------------------------
c  Binary search algorithm
c
 200  continue
      do i=1,ivec
c
c  first iteration
c
         if(init.eq.0) then
            iok(i)=.FALSE.
            imina(i)=1
            imaxa(i)=nx-1
c
            iv(i)=nx/2
c
            if(xuse(i).lt.xpkg(iv(i),1)) then
               imaxa(i)=max(1,(iv(i)-1))
            else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
               imina(i)=min((nx-1),(iv(i)+1))
            else
               iok(i)=.TRUE.
            endif
c
         endif
c
c  second iteration
c
         if(.not.iok(i)) then
c
            iv(i)=(imina(i)+imaxa(i))/2
c
            if(xuse(i).lt.xpkg(iv(i),1)) then
               imaxa(i)=max(1,(iv(i)-1))
            else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
               imina(i)=min((nx-1),(iv(i)+1))
            else
               iok(i)=.TRUE.
            endif
c
c  third iteration
c
            if(.not.iok(i)) then
c
               iv(i)=(imina(i)+imaxa(i))/2
c
               if(xuse(i).lt.xpkg(iv(i),1)) then
                  imaxa(i)=max(1,(iv(i)-1))
               else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                  imina(i)=min((nx-1),(iv(i)+1))
               else
                  iok(i)=.TRUE.
               endif
c
c  fourth iteration
c
               if(.not.iok(i)) then
c
                  iv(i)=(imina(i)+imaxa(i))/2
c
                  if(xuse(i).lt.xpkg(iv(i),1)) then
                     imaxa(i)=max(1,(iv(i)-1))
                  else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                     imina(i)=min((nx-1),(iv(i)+1))
                  else
                     iok(i)=.TRUE.
                  endif
c
c  fifth iteration
c
                  if(.not.iok(i)) then
c
                     iv(i)=(imina(i)+imaxa(i))/2
c
                     if(xuse(i).lt.xpkg(iv(i),1)) then
                        imaxa(i)=max(1,(iv(i)-1))
                     else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                        imina(i)=min((nx-1),(iv(i)+1))
                     else
                        iok(i)=.TRUE.
                     endif
c
                     if(.not.iok(i)) iprob=iprob+1
c
c  end chain of iteration if-then-else blocks
c
                  endif
               endif
            endif
         endif
c
c  end of loop
c
      enddo
c
      go to 500
c-----------------------------------------------------------
c  Binary search algorithm
c
 250  continue
c
      do i=1,ivec
c
c  first iteration
c
         if(init.eq.0) then
            if((xpkg(iprev,1).le.xuse(i)).and.
     >         (xuse(i).le.xpkg(iprev+1,1))) then
               iok(i)=.TRUE.
               iv(i)=iprev
            else
               iok(i)=.FALSE.
               imina(i)=1
               imaxa(i)=nx-1
c
               iv(i)=nx/2
c
               if(xuse(i).lt.xpkg(iv(i),1)) then
                  imaxa(i)=max(1,(iv(i)-1))
               else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                  imina(i)=min((nx-1),(iv(i)+1))
               else
                  iok(i)=.TRUE.
               endif
            endif
c
         endif
c
c  second iteration
c
         if(.not.iok(i)) then
c
            iv(i)=(imina(i)+imaxa(i))/2
c
            if(xuse(i).lt.xpkg(iv(i),1)) then
               imaxa(i)=max(1,(iv(i)-1))
            else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
               imina(i)=min((nx-1),(iv(i)+1))
            else
               iok(i)=.TRUE.
            endif
c
c  third iteration
c
            if(.not.iok(i)) then
c
               iv(i)=(imina(i)+imaxa(i))/2
c
               if(xuse(i).lt.xpkg(iv(i),1)) then
                  imaxa(i)=max(1,(iv(i)-1))
               else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                  imina(i)=min((nx-1),(iv(i)+1))
               else
                  iok(i)=.TRUE.
               endif
c
c  fourth iteration
c
               if(.not.iok(i)) then
c
                  iv(i)=(imina(i)+imaxa(i))/2
c
                  if(xuse(i).lt.xpkg(iv(i),1)) then
                     imaxa(i)=max(1,(iv(i)-1))
                  else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                     imina(i)=min((nx-1),(iv(i)+1))
                  else
                     iok(i)=.TRUE.
                  endif
c
c  fifth iteration
c
                  if(.not.iok(i)) then
c
                     iv(i)=(imina(i)+imaxa(i))/2
c
                     if(xuse(i).lt.xpkg(iv(i),1)) then
                        imaxa(i)=max(1,(iv(i)-1))
                     else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                        imina(i)=min((nx-1),(iv(i)+1))
                     else
                        iok(i)=.TRUE.
                     endif
c
                     if(.not.iok(i)) iprob=iprob+1
c
c  end chain of iteration if-then-else blocks
c
                  endif
               endif
            endif
         endif
c
c  end of loop
c
         iprev=iv(i)
      enddo
c
      go to 500
c-----------------------------------------------------------
c  algorithm:  piecewise linear indexing function lookup & correction
c
 300  continue
      do i=1,ivec
c
c  first iteration
c
         if(init.eq.0) then
            iok(i)=.FALSE.
c
c  piecewise linear indexing function on even spaced grid
c  (not same grid as x axis itself)
c
            iv(i)=1+havi*(xuse(i)-xpkg(1,1))
            iv(i)=max(1,min(nx-1,iv(i)))
            xfac=(xuse(i)-(xpkg(1,1)+(iv(i)-1)*hav))*havi
            zindx0=xpkg(iv(i),2)
            if(iv(i).lt.nx-1) then
               zdindx=xpkg(iv(i)+1,2)-zindx0
            else
               zdindx=nx-zindx0
            endif
            zindex=zindx0+xfac*zdindx
c
            iv(i)=zindex
            iv(i)=max(1,min(nx-1,iv(i)))
c
            if(xuse(i).lt.xpkg(iv(i),1)) then
               i_sign=-1
            else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
               i_sign=+1
            else
               iok(i)=.TRUE.
            endif
c
         endif
c
c  second iteration
c
         if(.not.iok(i)) then
            iv(i)=iv(i)+i_sign
            if(xuse(i).lt.xpkg(iv(i),1)) then
               i_sign=-1
            else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
               i_sign=+1
            else
               iok(i)=.TRUE.
            endif
c
c  third iteration
c
            if(.not.iok(i)) then
               iv(i)=iv(i)+i_sign
               if(xuse(i).lt.xpkg(iv(i),1)) then
                  i_sign=-1
               else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                  i_sign=+1
               else
                  iok(i)=.TRUE.
               endif
c
c  fourth iteration
c
               if(.not.iok(i)) then
                  iv(i)=iv(i)+i_sign
                  if(xuse(i).lt.xpkg(iv(i),1)) then
                     i_sign=-1
                  else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                     i_sign=+1
                  else
                     iok(i)=.TRUE.
                  endif
c
c  fifth iteration
c
                  if(.not.iok(i)) then
                     iv(i)=iv(i)+i_sign
                     if(xuse(i).lt.xpkg(iv(i),1)) then
                        i_sign=-1
                     else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                        i_sign=+1
                     else
                        iok(i)=.TRUE.
                     endif
c
                     if(.not.iok(i)) iprob=iprob+1
c
c  end chain of iteration if-then-else blocks
c
                  endif
               endif
            endif
         endif
c
c  end of loop
c
      enddo
c
      go to 500
c
c-----------------------------------------------------------
c  algorithm:  piecewise linear indexing function lookup & correction
c
 350  continue
      do i=1,ivec
c
c  first iteration
c
         if(init.eq.0) then
            if((xpkg(iprev,1).le.xuse(i)).and.
     >         (xuse(i).le.xpkg(iprev+1,1))) then
               iok(i)=.TRUE.
               iv(i)=iprev
            else
               iok(i)=.FALSE.
c
c  piecewise linear indexing function on even spaced grid
c  (not same grid as x axis itself)
c
               iv(i)=1+havi*(xuse(i)-xpkg(1,1))
               iv(i)=max(1,min(nx-1,iv(i)))
               xfac=(xuse(i)-(xpkg(1,1)+(iv(i)-1)*hav))*havi
               zindx0=xpkg(iv(i),2)
               if(iv(i).lt.nx-1) then
                  zdindx=xpkg(iv(i)+1,2)-zindx0
               else
                  zdindx=nx-zindx0
               endif
               zindex=zindx0+xfac*zdindx
c
               iv(i)=zindex
               iv(i)=max(1,min(nx-1,iv(i)))
c
               if(xuse(i).lt.xpkg(iv(i),1)) then
                  i_sign=-1
               else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                  i_sign=+1
               else
                  iok(i)=.TRUE.
               endif
            endif
c
         endif
c
c  second iteration
c
         if(.not.iok(i)) then
            iv(i)=iv(i)+i_sign
            if(xuse(i).lt.xpkg(iv(i),1)) then
               i_sign=-1
            else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
               i_sign=+1
            else
               iok(i)=.TRUE.
            endif
c
c  third iteration
c
            if(.not.iok(i)) then
               iv(i)=iv(i)+i_sign
               if(xuse(i).lt.xpkg(iv(i),1)) then
                  i_sign=-1
               else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                  i_sign=+1
               else
                  iok(i)=.TRUE.
               endif
c
c  fourth iteration
c
               if(.not.iok(i)) then
                  iv(i)=iv(i)+i_sign
                  if(xuse(i).lt.xpkg(iv(i),1)) then
                     i_sign=-1
                  else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                     i_sign=+1
                  else
                     iok(i)=.TRUE.
                  endif
c
c  fifth iteration
c
                  if(.not.iok(i)) then
                     iv(i)=iv(i)+i_sign
                     if(xuse(i).lt.xpkg(iv(i),1)) then
                        i_sign=-1
                     else if(xuse(i).gt.xpkg(iv(i)+1,1)) then
                        i_sign=+1
                     else
                        iok(i)=.TRUE.
                     endif
c
                     if(.not.iok(i)) iprob=iprob+1
c
c  end chain of iteration if-then-else blocks
c
                  endif
               endif
            endif
         endif
c
c  end of loop
c
         iprev=iv(i)
      enddo
c
      go to 500
c
c--------------------------------------------------------------------
c  any "problems" left? if so, re-enter the loop...
c
 500  continue
      if(iprob.gt.0) go to 5
c
c  OK -- all zones found; complete output stats
c
      if(imode.eq.1) then
         dxn=(xuse-xpkg(iv,1))          ! un-normalized
      else if(ialg.ne.3) then
         dxn=(xuse-xpkg(iv,1))*xpkg(iv,3) ! normalized to 1/h in each zone
         hv=xpkg(iv,2)
         hiv=xpkg(iv,3)
      else
         dxn=(xuse-xpkg(iv,1))*xpkg(iv,3) ! normalized to 1/h in each zone
         hv=xpkg(iv+1,1)-xpkg(iv,1)
         hiv=xpkg(iv,3)
      endif
c
      deallocate(iok)
      if(ialg.lt.3) then
         deallocate(imina)
         deallocate(imaxa)
      endif
c
c  all done -- generalized lookup
c
 1000 continue
      deallocate(xuse)
c
 1100 continue
      return
      end
