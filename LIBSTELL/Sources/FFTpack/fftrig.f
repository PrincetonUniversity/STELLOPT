      SUBROUTINE fftrig_g (trigs, n, mode)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER n, mode
      REAL(rprec), DIMENSION(*) :: trigs
!DEC$ IF .NOT.DEFINED (CRAY) .OR. DEFINED(LONESTAR) .OR. DEFINED(MCURIE)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1, two = 2, p5 = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: imode, nn, l, i, nh, la
      REAL(rprec) :: pi, del, angle
C-----------------------------------------------
      pi=two*ASIN(one)
      imode=IABS(mode)
      nn=n
      IF (imode.gt.1.and.imode.lt.6) nn=n/2
      del=(pi+pi)/nn
      l=nn+nn
      DO 10 i=1,l,2
      angle=(p5*del)*(i-1)
      trigs(i)=COS(angle)
      trigs(i+1)=SIN(angle)
   10 CONTINUE
      IF (imode.eq.1) RETURN
      IF (imode.eq.8) RETURN
      del=del/2
      nh=(nn+1)/2
      l=nh+nh
      la=nn+nn
      DO 20 i=1,l,2
      angle=(i-1)*(del/2)
      trigs(la+i)=COS(angle)
      trigs(la+i+1)=SIN(angle)
   20 CONTINUE
      IF (imode.le.3) RETURN
      del=del/2
      la=la+nn
      IF (mode.eq.5) GOTO 40
      DO 30 i=2,nn
      angle=(i-1)*del
      trigs(la+i)=two*SIN(angle)
   30 CONTINUE
      RETURN
   40 CONTINUE
      del=del/2
      DO 50 i=2,n
      angle = (i-1)*del
      trigs(la+i)=SIN(angle)
   50 CONTINUE
!DEC$ ELSE
      CALL fftrig (trigs, n, mode)
!DEC$ ENDIF
      END SUBROUTINE fftrig_g
