
SUBROUTINE spl2d(nx,ny,hx,hy,mx,my,f,spl)

! Makes a 2-dimensional cubic spline of function f(x,y)
!
! Input:  nx, ny              number of values in x and y
!         hx, hy              step size in x and y (aequidistant)
!         mx, my              spline mode (0: standard, 1: periodic)
!         f(nx,ny)            f(x,y)-values
! Output: spl                 Array with spline parameters

  USE neo_precision

  IMPLICIT NONE

  INTEGER,                             INTENT(in)  :: nx, ny, mx, my
  REAL(kind=dp),                       INTENT(in)  :: hx, hy
  REAL(kind=dp), DIMENSION(nx,ny)    , INTENT(in)  :: f
  REAL(kind=dp), DIMENSION(4,4,nx,ny), INTENT(out) :: spl

  REAL(kind=dp), DIMENSION(:),         ALLOCATABLE :: bi, ci, di, s
  INTEGER                                          :: i, j, k, l

  ALLOCATE ( bi(nx), ci(nx), di(nx), s(nx) )
  DO j = 1,ny
     DO i = 1,nx
        s(i) = f(i,j)
     END DO
     IF (mx .EQ. 0) THEN
        CALL splreg(nx,hx,s,bi,ci,di)
     ELSE
        CALL splper(nx,hx,s,bi,ci,di)
     ENDIF
     DO i = 1,nx
        spl(1,1,i,j) = s(i)
        spl(2,1,i,j) = bi(i)
        spl(3,1,i,j) = ci(i)
        spl(4,1,i,j) = di(i)
     END DO
  END DO
  DEALLOCATE ( bi, ci, di, s )

  ALLOCATE ( bi(ny), ci(ny), di(ny), s(ny) )
  DO k = 1,4
     DO i = 1,nx
        DO j = 1,ny
           s(j) = spl(k,1,i,j)
        END DO
        IF (my .EQ. 0) THEN
           CALL splreg(ny,hy,s,bi,ci,di)
        ELSE
           CALL splper(ny,hy,s,bi,ci,di)
        ENDIF
        DO j=1,ny
           spl(k,2,i,j)=bi(j)
           spl(k,3,i,j)=ci(j)
           spl(k,4,i,j)=di(j)
        END DO
     END DO
  END DO
  DEALLOCATE ( bi, ci, di, s )

  RETURN
END SUBROUTINE spl2d
