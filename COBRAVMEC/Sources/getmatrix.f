      SUBROUTINE getmatrix(nsurf, ad, asup, asub, h, n, ik, xf,
     1  xh, pf, ph, qf, qh, rf, rh, inc, xmin)
      USE stel_kinds
      USE ballooning_data
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: nsurf, n, ik, inc
      REAL(rprec), INTENT(in):: h, xmin
      REAL(rprec), INTENT(out):: ad(n-2+inc), asup(n-3+inc),
     1  asub(n-2+inc)
      REAL(rprec),INTENT(inout):: xf(n), xh(n-1), pf(n), ph(n-1),
     1  qf(n), qh(n-1), rf(n), rh(n-1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec) :: h2
      INTEGER :: k
!------------------------------------------------------------------------

      h2=h*h
      IF(ik.eq.1)then                                              ! IF first CALL THEN.....
         xf=(/(xmin+(k-1)*h,k=1,n)/)                               ! .... build full mesh and ...
         CALL coeffs(nsurf, n, xf, pf, qf, rf)                     ! ...evaluate coefficients on full mesh
         IF (lfail_balloon) RETURN
      ENDIF

      xh(1:n-1)=xf(2:n)-h/2                                        ! build half mesh and ...
      CALL coeffs(nsurf, n-1, xh, ph, qh, rh)                      ! ...evaluate coefficients on half mesh
      IF (lfail_balloon) RETURN


      IF(tsymm .eq. 0)then

         ad(1)=(-2*ph(1)+qf(1)*h2)/(rf(1)*h2)                      !  B.C.: F(1)=F(-1) , i.e. F'(0)=0
         asup(1)=2*ph(1)/(rf(1)*h2)                                !  SIZE of matrix is THEN NPTS-1

         ad(2:n-1)=(-ph(2:n-1)-ph(1:n-2)+qf(2:n-1)*h2)
     1            /(rf(2:n-1)*h2)                                  ! build A-diagonal
         asup(2:n-2)=ph(2:n-2)/(rf(2:n-2)*h2)                      ! build A-supdiagonal
         asub(1:n-2)=ph(1:n-2)/(rf(2:n-1)*h2)                      ! build A-subdiagonal

      ELSE                                                         !  SIZE of matrix is NPTS-2

         ad(1:n-2)=(-ph(2:n-1)-ph(1:n-2)+qf(2:n-1)*h2)
     1            /(rf(2:n-1)*h2)                                  ! build A-diagonal
         asup(1:n-3)=ph(2:n-2)/(rf(2:n-2)*h2)                      ! build A-supdiagonal
         asub(1:n-3)=ph(2:n-2)/(rf(3:n-1)*h2)                      ! build A-subdiagonal

      ENDIF

      END SUBROUTINE getmatrix
