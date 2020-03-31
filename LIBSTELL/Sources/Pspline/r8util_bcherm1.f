      subroutine r8util_bcherm1(fherm,idimx1,
     >     jbcxa,jbcxb,
     >     zbcxa,zbcxb,
     >     x1)
 
C...  insert BCs as needed for Hermite interpolation setup
 
 
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer :: idimx1         ! array dimensions
      REAL*8 :: fherm(0:1,idimx1)
 
      integer :: jbcxa,jbcxb    ! x BC controls
      REAL*8 :: zbcxa            ! x(1) BC data
      REAL*8 :: zbcxb            ! x(idimx1) BC data
 
      REAL*8 :: x1(idimx1)
 
C-------------------------------------------------------------
      integer :: ix
      REAL*8 :: zdx
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
 
      end subroutine r8util_bcherm1
