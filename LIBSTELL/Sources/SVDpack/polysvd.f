!CALL polysvd(kpp,strial,ptrial,SIZE(strial),cfnt(1:kpp))
      SUBROUTINE polysvd(nfitin,x,y,m,b)
! USAGE SVD regression of input vectors x, y of length m
! producess nfitin polynoimial coefficients b
      USE stel_kinds
      USE stel_constants
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m, nfitin
      REAL (rprec), DIMENSION(m), INTENT(in) :: x,y
      REAL(rprec) ,DIMENSION(nfitin), INTENT(out) :: b
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL (rprec) :: cutoff=1.e-7_dp
      REAL (rprec), DIMENSION(:,:), ALLOCATABLE
     1          :: amatrix, vv ,uu, wwd
      REAL (rprec), DIMENSION(:) , ALLOCATABLE :: ww
      REAL (rprec), DIMENSION(:), ALLOCATABLE, SAVE :: apar
      INTEGER , DIMENSION(:) , ALLOCATABLE :: pwr
      INTEGER ::  i, j, n=11
      n=nfitin
      ALLOCATE(amatrix(m,n))
      ALLOCATE(uu(m,n))
      ALLOCATE(pwr(n))
      ALLOCATE(apar(n))
      ALLOCATE(wwd(n,n))
      ALLOCATE(ww(n))
      ALLOCATE(vv(n,n))
      DO i=1,n
         pwr(i)=i-1
      ENDDO
! make amatrix
      DO i=1,m
         DO j=1,n
           IF(x(i).eq.0.) THEN
              amatrix(i,j)=zero
           ELSE
              amatrix(i,j)=x(i)**pwr(j)
           ENDIF
         ENDDO 
      ENDDO
      uu=amatrix
      CALL svdcmp(uu,m,n,m,n,ww,vv) ! m rows > n columns
      DO i=1,n
         DO j=1,n
           wwd(i,j)=0  ! reset after matmul
         ENDDO      
      ENDDO      
      DO i=1,n
         wwd(i,i)=1/ww(i)
      ENDDO
      DO i= 1,SIZE(ww)
         IF(ww(i) .lt. cutoff) wwd(i,i)=0
      ENDDO
      apar=0.
      apar=matmul(vv,matmul(wwd,matmul(TRANSPOSE(uu),y)))
      b(1:n)=apar(1:n)
      DEALLOCATE(amatrix)
      DEALLOCATE(uu)
      DEALLOCATE(pwr)
      DEALLOCATE(apar)
      DEALLOCATE(wwd)
      DEALLOCATE(ww)
      DEALLOCATE(vv)
      END SUBROUTINE polysvd

!      PROGRAM driver
! TESTED against IDL with this PARAMETER set. Results agree.
!      USE stel_kinds
!      USE stel_constants
!      INTEGER, PARAMETER :: ns=50, kpp=5
!      REAL(rprec), DIMENSION(ns) :: s, y
!      REAL(rprec), DIMENSION(kpp) :: b
!      DO i=1,ns 
!       s(i)=one*(i-1)/(ns-1)
!      ENDDO
!      y=COS(s*pio2)**2
!      CALL polysvd(kpp,s,y,SIZE(s),b)
!      WRITE(6,109)b
!      WRITE(6,109)s
!      WRITE (6,109)y
!109   FORMAT(/,'    =[',3(x,1pe14.6,1h,),'$')
!      END PROGRAM


