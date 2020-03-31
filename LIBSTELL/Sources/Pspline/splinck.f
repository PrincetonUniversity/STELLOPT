      subroutine splinck(x,inx,ilinx,ztol,ier)
C
C  check if a grid is strictly ascending and if it is evenly spaced
C  to w/in ztol
C
      real x(inx)                       ! input -- grid to check
C
      integer ilinx                     ! output -- =1 if evenly spaced =2 O.W.
C
      real ztol                         ! input -- spacing check tolerance
C
      integer ier                       ! output -- =0 if OK
C
C  ier=1 is returned if x(1...inx) is NOT STRICTLY ASCENDING...
C
C-------------------------------
C
      ier=0
      ilinx=1
      if(inx.le.1) return
c
      dxavg=(x(inx)-x(1))/(inx-1)
      zeps=abs(ztol*dxavg)
c
      do ix=2,inx
         zdiffx=(x(ix)-x(ix-1))
         if(zdiffx.le.0.0) ier=2
         zdiff=zdiffx-dxavg
         if(abs(zdiff).gt.zeps) then
            ilinx=2
         endif
      enddo
 10   continue
c
      return
      end
