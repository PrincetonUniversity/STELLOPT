

!ref: http://www.femto-st.fr/~daniau/fvn.html
      SUBROUTINE fvn_akima(n,x,y,br,co)
       USE precision
          IMPLICIT NONE
          INTEGER, INTENT(in)  :: n
          REAL(rprec), INTENT(in) :: x(n)
          REAL(rprec), INTENT(in) :: y(n)
          REAL(rprec), INTENT(out) :: br(n)
          REAL(rprec), INTENT(out) :: co(4,n)
          
          REAL(rprec), ALLOCATABLE :: var(:),z(:)
          REAL(rprec) :: wi_1,wi
          INTEGER :: i
          REAL(rprec) :: dx,a,b
          
          ! br is just a copy of x
          br(:)=x(:)
      
          ALLOCATE(var(n+3))
          ALLOCATE(z(n))
          ! evaluate the variations
          DO i=1, n-1
              var(i+2)=(y(i+1)-y(i))/(x(i+1)-x(i))
          ENDDO
          var(n+2)=2.d0*var(n+1)-var(n)
          var(n+3)=2.d0*var(n+2)-var(n+1)
          var(2)=2.d0*var(3)-var(4)
          var(1)=2.d0*var(2)-var(3)
        
          DO i = 1, n
          wi_1=ABS(var(i+3)-var(i+2))
          wi=ABS(var(i+1)-var(i))
            IF ((wi_1+wi) .EQ. 0.d0) THEN
              z(i)=(var(i+2)+var(i+1))/2.d0
            ELSE
              z(i)=(wi_1*var(i+1)+wi*var(i+2))/(wi_1+wi)
            ENDIF
          ENDDO
          
          DO i=1, n-1
              dx=x(i+1)-x(i)
              a=(z(i+1)-z(i))*dx ! coeff INTermediaires pour calcul wd
              b=y(i+1)-y(i)-z(i)*dx ! coeff INTermediaires pour calcul wd
              co(1,i)=y(i)
              co(2,i)=z(i)
              co(3,i)=(3.d0*var(i+2)-2.d0*z(i)-z(i+1))/dx   
              co(4,i)=(z(i)+z(i+1)-2.d0*var(i+2))/dx**2  !
          ENDDO
          co(1,n)=y(n)
          co(2,n)=z(n)
          co(3,n)=0.d0
          co(4,n)=0.d0
      
          DEALLOCATE(z)
          DEALLOCATE(var)
      
      END SUBROUTINE

      FUNCTION fvn_spline_eval(x,n,br,co)
! VALUES
       USE precision
          IMPLICIT NONE
          REAL(rprec) fvn_spline_eval
          REAL(rprec), INTENT(in) :: x           
! x must be br(1)<= x <= br(n+1) otherwise value is extrapolated
          INTEGER, INTENT(in) :: n        ! number of INTervals
          REAL(rprec), INTENT(in) :: br(n+1)  ! breakpoints (abcissas)
          REAL(rprec), INTENT(in) :: co(4,n+1)! spline coeeficients
          REAL(rprec) :: fvn_d_spline_eval
          
          INTEGER :: i
          REAL(rprec) :: dx
          
          
          IF (x<=br(1)) THEN
              i=1
          ELSE IF (x>=br(n+1)) THEN
              i=n
          ELSE
            i=1
          DO WHILE(x>=br(i))
              i=i+1
          ENDDO
            i=i-1
          ENDIF
          
          dx=x-br(i)
          fvn_spline_eval=co(1,i)+co(2,i)*dx+co(3,i)*dx**2+co(4,i)*dx**3
      
      END FUNCTION fvn_spline_eval


      FUNCTION fvn_spline_devaldx(x,n,br,co)
! 1ST DERIVATIVE
       USE precision
          IMPLICIT NONE
          REAL(rprec) fvn_spline_devaldx
          REAL(rprec), INTENT(in) :: x           
! x must be br(1)<= x <= br(n+1) otherwise value is extrapolated
          INTEGER, INTENT(in) :: n        ! number of INTervals
          REAL(rprec), INTENT(in) :: br(n+1)  ! breakpoints
          REAL(rprec), INTENT(in) :: co(4,n+1)! spline coeeficients
          REAL(rprec) :: fvn_d_spline_eval
          
          INTEGER :: i
          REAL(rprec) :: dx
          
          
          IF (x<=br(1)) THEN
              i=1
          ELSE IF (x>=br(n+1)) THEN
              i=n
          ELSE
            i=1
          DO WHILE(x>=br(i))
              i=i+1
          ENDDO
            i=i-1
          ENDIF
          
          dx=x-br(i)
          fvn_spline_devaldx=co(2,i)+2*co(3,i)*dx+3*co(4,i)*dx**2
      
      END FUNCTION fvn_spline_devaldx

      FUNCTION fvn_spline_d2evaldx2(x,n,br,co)
! 2ND DERIVATIVE
       USE precision
          IMPLICIT NONE
          REAL(rprec) fvn_spline_d2evaldx2
          REAL(rprec), INTENT(in) :: x           
! x must be br(1)<= x <= br(n+1) otherwise value is extrapolated
          INTEGER, INTENT(in) :: n        ! number of INTervals
          REAL(rprec), INTENT(in) :: br(n+1)  ! breakpoints
          REAL(rprec), INTENT(in) :: co(4,n+1)! spline coeeficients
          
          INTEGER :: i
          REAL(rprec) :: dx
          
          
          IF (x<=br(1)) THEN
              i=1
          ELSE IF (x>=br(n+1)) THEN
              i=n
          ELSE
            i=1
          DO WHILE(x>=br(i))
              i=i+1
          ENDDO
            i=i-1
          ENDIF
          
          dx=x-br(i)
          fvn_spline_d2evaldx2=2*co(3,i)+6*co(4,i)*dx
      
      END FUNCTION fvn_spline_d2evaldx2

