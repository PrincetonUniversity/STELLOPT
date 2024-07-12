      subroutine util_bcherm1(fherm,idimx1,
     >     jbcxa,jbcxb,
     >     zbcxa,zbcxb,
     >     x1)

C...  insert BCs as needed for Hermite interpolation setup

      implicit NONE

      integer :: idimx1         ! array dimensions
      real :: fherm(0:1,idimx1)

      integer :: jbcxa,jbcxb    ! x BC controls
      real :: zbcxa            ! x(1) BC data
      real :: zbcxb            ! x(idimx1) BC data

      real :: x1(idimx1)

C-------------------------------------------------------------
      integer :: ix
      real :: zdx
C-------------------------------------------------------------

      if((jbcxa.eq.1).or.(jbcxb.eq.1)) then

         if(jbcxa.eq.1) then
            fherm(1,1)=zbcxa
         else
            zdx = x1(2)-x1(1)
            fherm(1,1)=(fherm(0,2)-fherm(0,1))/zdx
         endif

         if(jbcxb.eq.1) then
            fherm(1,idimx1)=zbcxb
         else
            zdx = x1(idimx1)-x1(idimx1-1)
            fherm(1,idimx1)=(fherm(0,idimx1)-fherm(0,idimx1-1))/zdx
         endif

      endif

      end subroutine util_bcherm1
