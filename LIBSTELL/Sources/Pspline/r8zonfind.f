      subroutine r8zonfind(x,nx,zxget,i)
c
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nx,nxm,i1,i2,ij,ii
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 dx
!============
      REAL*8 x(nx),zxget
      integer i
c
c  find index i such that x(i).le.zxget.le.x(i+1)
c
c  x(1...nx) is strict increasing and x(1).le.zxget.le.x(nx)
c  (this is assumed to already have been checked -- no check here!)
c
      nxm=nx-1
      if((i.lt.1).or.(i.gt.nxm)) then
         i1=1
         i2=nx-1
         go to 10
      endif
c
      if(x(i).gt.zxget) then
c  look down
         dx=x(i+1)-x(i)
         if((x(i)-zxget).gt.4*dx) then
            i1=1
            i2=i-1
            go to 10
         else
            i2=i-1
            do ij=i2,1,-1
               if((x(ij).le.zxget).and.(zxget.le.x(ij+1))) then
                  i=ij
                  return
               endif
            enddo
            i=1
            return
         endif
      else if(x(i+1).lt.zxget) then
c  look up
         dx=x(i+1)-x(i)
         if((zxget-x(i+1)).gt.4*dx) then
            i1=i+1
            i2=nxm
            go to 10
         else
            i2=i+1
            do ij=i2,nxm
               if((x(ij).le.zxget).and.(zxget.le.x(ij+1))) then
                  i=ij
                  return
               endif
            enddo
            ij=nxm
            return
         endif
      else
c  already there...
         return
      endif
c
c---------------------------
c  binary search
c
 10   continue
c
      if(i1.eq.i2) then
c found by proc. of elimination
         i=i1
         return
      endif
c
      ii=(i1+i2)/2
c
      if(zxget.lt.x(ii)) then
         i2=ii-1
      else if(zxget.gt.x(ii+1)) then
         i1=ii+1
      else
c found
         i=ii
         return
      endif
c
      go to 10
c
      end
 
