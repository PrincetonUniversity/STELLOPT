      subroutine r8util_bcherm3(fherm,idimx1,idimx2,idimx3,
     >     jbcx1a,jbcx1b,   jbcx2a,jbcx2b,   jbcx3a,jbcx3b,
     >     zbcx1a,zbcx1b,   zbcx2a,zbcx2b,   zbcx3a,zbcx3b,
     >     x1, x2, x3)
 
C... insert BCs as needed for Hermite interpolation setup
 
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer :: idimx1,idimx2,idimx3 ! array dimensions
      REAL*8 :: fherm(0:7,idimx1,idimx2,idimx3)
 
      integer :: jbcx1a,jbcx1b  ! x1 BC controls
      REAL*8 :: zbcx1a(idimx2,idimx3) ! x1(1) BC data
      REAL*8 :: zbcx1b(idimx2,idimx3) ! x1(idimx1) BC data
 
      integer :: jbcx2a,jbcx2b  ! x2 BC controls
      REAL*8 :: zbcx2a(idimx1,idimx3) ! x2(1) BC data
      REAL*8 :: zbcx2b(idimx1,idimx3) ! x2(idimx2) BC data
 
      integer :: jbcx3a,jbcx3b  ! x3 BC controls
      REAL*8 :: zbcx3a(idimx2,idimx3) ! x3(1) BC data
      REAL*8 :: zbcx3b(idimx2,idimx3) ! x3(idimx3) BC data
 
      REAL*8 :: x1(idimx1),x2(idimx2),x3(idimx3)
 
C-----------------------------------------------------------------
      integer :: ix,iy
      REAL*8 :: zdx
C-----------------------------------------------------------------
 
      if((jbcx1a.eq.1).or.(jbcx1b.eq.1)) then
 
         if(jbcx1a.eq.1) then
            fherm(1,1,1:idimx2,1:idimx3)=zbcx1a(1:idimx2,1:idimx3)
         else
            zdx = x1(2)-x1(1)
            do iy=1,idimx3
               do ix=1,idimx2
                  fherm(1,1,ix,iy)=
     >                 (fherm(0,2,ix,iy)-fherm(0,1,ix,iy))/zdx
               enddo
            enddo
         endif
 
         if(jbcx1b.eq.1) then
            fherm(1,idimx1,1:idimx2,1:idimx3)=zbcx1b(1:idimx2,1:idimx3)
         else
            zdx = x1(idimx1)-x1(idimx1-1)
            do iy=1,idimx3
               do ix=1,idimx2
                  fherm(1,idimx1,ix,iy)=
     >                 (fherm(0,idimx1,ix,iy)-fherm(0,idimx1-1,ix,iy))/
     >                 zdx
               enddo
            enddo
         endif
 
      endif
 
      if((jbcx2a.eq.1).or.(jbcx2b.eq.1)) then
 
         if(jbcx2a.eq.1) then
            fherm(2,1:idimx1,1,1:idimx3)=zbcx2a(1:idimx1,1:idimx3)
         else
            zdx=x2(2)-x2(1)
            do iy=1,idimx3
               do ix=1,idimx1
                  fherm(2,ix,1,iy)=
     >                 (fherm(0,ix,2,iy)-fherm(0,ix,1,iy))/zdx
               enddo
            enddo
         endif
 
         if(jbcx2b.eq.1) then
            fherm(2,1:idimx1,idimx2,1:idimx3)=zbcx2b(1:idimx1,1:idimx3)
         else
            zdx=x2(idimx2)-x2(idimx2-1)
            do iy=1,idimx3
               do ix=1,idimx1
                  fherm(2,ix,idimx2,iy)=
     >                 (fherm(0,ix,idimx2,iy)-fherm(0,ix,idimx2-1,iy))/
     >                 zdx
               enddo
            enddo
         endif
      endif
 
      if((jbcx3a.eq.1).or.(jbcx3b.eq.1)) then
 
         if(jbcx3a.eq.1) then
            fherm(3,1:idimx1,1:idimx2,1)=zbcx3a(1:idimx1,1:idimx2)
         else
            zdx=x3(2)-x3(1)
            do iy=1,idimx2
               do ix=1,idimx1
                  fherm(3,ix,iy,1)=
     >                 (fherm(0,ix,iy,2)-fherm(0,ix,iy,1))/zdx
               enddo
            enddo
         endif
 
         if(jbcx3b.eq.1) then
            fherm(3,1:idimx1,1:idimx2,idimx3)=zbcx3b(1:idimx1,1:idimx2)
         else
            zdx=x3(idimx3)-x3(idimx3-1)
            do iy=1,idimx2
               do ix=1,idimx1
                  fherm(3,ix,iy,idimx3)=
     >                 (fherm(0,ix,iy,idimx3)-fherm(0,ix,iy,idimx3-1))/
     >                 zdx
               enddo
            enddo
         endif
      endif
 
      end subroutine r8util_bcherm3
