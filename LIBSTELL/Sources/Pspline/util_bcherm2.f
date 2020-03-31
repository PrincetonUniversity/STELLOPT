      subroutine util_bcherm2(fherm,idimx1,idimx2,
     >     jbcx1a,jbcx1b,   jbcx2a,jbcx2b,
     >     zbcx1a,zbcx1b,   zbcx2a,zbcx2b,
     >     x1, x2)

C...  insert BCs as needed for Hermite interpolation setup

      implicit NONE

      integer :: idimx1,idimx2  ! array dimensions
      real :: fherm(0:3,idimx1,idimx2)

      integer :: jbcx1a,jbcx1b  ! x1 BC controls
      real :: zbcx1a(idimx2)    ! x1(1) BC data
      real :: zbcx1b(idimx2)    ! x1(idimx1) BC data

      integer :: jbcx2a,jbcx2b  ! x2 BC controls
      real :: zbcx2a(idimx1)    ! x2(1) BC data
      real :: zbcx2b(idimx1)    ! x2(idimx2) BC data

      real :: x1(idimx1),x2(idimx2)

C-------------------------------------------------------------
      integer :: ix
      real :: zdx
C-------------------------------------------------------------

      if((jbcx1a.eq.1).or.(jbcx1b.eq.1)) then

         if(jbcx1a.eq.1) then
            fherm(1,1,1:idimx2)=zbcx1a(1:idimx2)
         else
            zdx = x1(2)-x1(1)
            do ix=1,idimx2
               fherm(1,1,ix)=(fherm(0,2,ix)-fherm(0,1,ix))/zdx
            enddo
         endif

         if(jbcx1b.eq.1) then
            fherm(1,idimx1,1:idimx2)=zbcx1b(1:idimx2)
         else
            zdx = x1(idimx1)-x1(idimx1-1)
            do ix=1,idimx2
               fherm(1,idimx1,ix)=
     >              (fherm(0,idimx1,ix)-fherm(0,idimx1-1,ix))/zdx
            enddo
         endif

      endif

      if((jbcx2a.eq.1).or.(jbcx2b.eq.1)) then

         if(jbcx2a.eq.1) then
            fherm(2,1:idimx1,1)=zbcx2a(1:idimx1)
         else
            zdx=x2(2)-x2(1)
            do ix=1,idimx1
               fherm(2,ix,1)=(fherm(0,ix,2)-fherm(0,ix,1))/zdx
            enddo
         endif

         if(jbcx2b.eq.1) then
            fherm(2,1:idimx1,idimx2)=zbcx2b(1:idimx1)
         else
            zdx=x2(idimx2)-x2(idimx2-1)
            do ix=1,idimx1
               fherm(2,ix,idimx2)=
     >              (fherm(0,ix,idimx2)-fherm(0,ix,idimx2-1))/zdx
            enddo
         endif
      endif

      end subroutine util_bcherm2
