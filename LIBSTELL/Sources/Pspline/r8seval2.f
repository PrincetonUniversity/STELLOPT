C******************** START FILE SEVAL2.FOR ; GROUP TRKRLIB ************
C
      REAL*8 FUNCTION r8seval2(N, U, X, Y, B, C, D, deriv)
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 zslop
!============
      INTEGER N
      REAL*8  U, X(N), Y(N), B(N), C(N), D(N), deriv
C
CCCCCCCCCCCCCCC
CCCCCCCCCCCCCCC
C  THIS SUBROUTINE EVALUATES THE CUBIC SPLINE FUNCTION
C
C    SEVAL2 = Y(I) + B(I)*(U-X(I)) + C(I)*(U-X(I))**2 + D(I)*(U-X(I))**3
C
C  and the derivative
C
C    deriv = B(i) + 2*C(i)*(U-X(i)) + 3*D(i)*(U-X(i))**2
C
C  WHERE  X(I) .LT. U .LT. X(I+1), USING HORNER'S RULE
C
C  IF  U .LT. X(1) THEN  I = 1  IS USED.
C  IF  U .GE. X(N) THEN  I = N  IS USED.
C
C  INPUT..
C
C    N = THE NUMBER OF DATA POINTS
C    U = THE ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED
C    X,Y = THE ARRAYS OF DATA ABSCISSAS AND ORDINATES
C    B,C,D = ARRAYS OF SPLINE COEFFICIENTS COMPUTED BY SPLINE
C
C  IF  U  IS NOT IN THE SAME INTERVAL AS THE PREVIOUS CALL, THEN A
C  BINARY SEARCH IS PERFORMED TO DETERMINE THE PROPER INTERVAL.
C
      INTEGER I, J, K, ILIN
      REAL*8 DX
      DATA I/1/
      DATA ILIN/0/
C
      IF ( I .GE. N ) I = 1
      IF ( U .LT. X(I) ) GO TO 10
      IF ( U .LE. X(I+1) ) GO TO 30
C
C  BINARY SEARCH
C
   10 I = 1
      J = N+1
   20 K = (I+J)/2
      IF ( U .LT. X(K) ) J = K
      IF ( U .GE. X(K) ) I = K
      IF ( J .GT. I+1 ) GO TO 20
C
C  EVALUATE SPLINE
C
   30 DX = U - X(I)
      IF(ILIN.EQ.0) THEN
        r8seval2 = Y(I) + DX*(B(I) + DX*(C(I) + DX*D(I)))
        deriv = B(i) + DX*(2.0_r8*C(i) + DX*3.0_r8*D(i))
      ELSE
        IF(I.EQ.N) THEN
          ZSLOP=(Y(N)-Y(N-1))/(X(N)-X(N-1))
        ELSE
          ZSLOP=(Y(I+1)-Y(I))/(X(I+1)-X(I))
        ENDIF
        r8seval2=Y(I)+DX*ZSLOP
        deriv=ZSLOP
      ENDIF
C
      RETURN
      END
C******************** END FILE SEVAL2.FOR ; GROUP TRKRLIB **************
