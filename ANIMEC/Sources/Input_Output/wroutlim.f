      SUBROUTINE wroutlim(nwout)
      USE vmec_main
      USE vsvd
      USE mgrid_mod, ONLY: nlim_max, nlim, limitr, seplim,
     1                     rlim, zlim, reslim
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: nwout
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER l1, l2, ltouch, limpts, n, i
      INTEGER, DIMENSION(nlim_max*nlim) :: l1min,l2min
      REAL(rprec) :: sepmin
      CHARACTER, DIMENSION(nlim_max,nlim) :: closest
C-----------------------------------------------


      IF (nlim .le. 0) RETURN

      DO l1 = 1, nlim
         DO l2 = 1, limitr(l1)
            closest(l2,l1) = ' '
            seplim(l2,l1) = SQRT(ABS(seplim(l2,l1)))
         END DO
      END DO

      sepmin = seplim(1,1)
      DO l1 = 1, nlim
         DO l2 = 1, limitr(l1)
            IF (seplim(l2,l1) .lt. sepmin)
     1      sepmin = seplim(l2,l1)
         END DO
      END DO

!     Do it this way to catch multiple MINs
      ltouch = 0
      DO l1 = 1, nlim
         DO l2 = 1, limitr(l1)
            IF (ABS(seplim(l2,l1)-sepmin) .lt. 1.e-5_dp) THEN
               closest(l2,l1) = '*'
               ltouch = ltouch + 1
               l2min(ltouch) = l2
               l1min(ltouch) = l1
            ENDIF
         END DO
      END DO
      WRITE (nwout, 702) ltouch
      WRITE (nwout, 704) (l1min(i),l2min(i),rlim(l2min(i),l1min(i)),
     1   zlim(l2min(i),l1min(i)),seplim(l2min(i),l1min(i)),i=1,ltouch)
  702 FORMAT(8i10)
  704 FORMAT(2i5,3e20.13)
      IF (lrecon) THEN
      WRITE (nthreed, 705)
  705 FORMAT(/,' PLASMA BOUNDARY-LIMITER PROXIMITY'/
     1   ' POINTS OUTSIDE PLASMA HAVE RESIDUE < 0.5')
      WRITE (nthreed, 707)
  707 FORMAT(/,13x,' nlim    n       R         Z        Residue',
     1   '          MIN |d|        nearest'/13x,' ----',4x,'-',7x,'-',9x
     2   ,'-',8x,'-------',10x,'-------',8x,'-------')

      DO n = 1, nlim
         limpts = limitr(n)
         DO i = 1, limpts
            WRITE (nthreed, 708) n, i, rlim(i,n), zlim(i,n), reslim(i,n)
     1         , seplim(i,n), closest(i,n)
         END DO
      END DO
  708 FORMAT(13x,i5,i5,2f10.3,1p,e15.4,5x,e12.4,8x,a)
      END IF

      END SUBROUTINE wroutlim
