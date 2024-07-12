      SUBROUTINE header
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vnamecl2
      USE Vimatrix, ONLY: ioout
      USE dkes_input
      USE dkes_realspace, ONLY: mvalue, nvalue, mpnt, borbi1
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: m, n, mm, mn, mnmax
      INTEGER :: index1(mpnt)
C-----------------------------------------------
      WRITE (ioout, 10) dashes, chip, psip, btheta, bzeta, vp,
     1   nzperiod, mpolb, ntorb, dashes, blabl(ibbi), blabl(ibbi),
     2   blabl(ibbi), blabl(ibbi)
 10   FORMAT(50x,'VARIATIONAL DKES-II CODE - 06/2001'/,
     1       54x,'UNITS ARE ASSUMED TO BE MKS',//,1x,131a,/,8x,
     1   'CHIP',8x,'PSIP',6x,'BTHETA',7x,'BZETA',10x,"V'",4x,
     2   'NZPERIOD',7x,'MPOLB',7x,'NTORB',/
     3   1p,5e12.4,3i12/1x,131a,/,4x,4('N',4x,'M',14x,a3,9x))

      mnmax = 0
      DO mn = 1, mpnt
         IF (borbi1(mn) .eq. zero) CYCLE
         mnmax = mnmax + 1
         index1(mnmax) = mn
      END DO
      WRITE (ioout, 20) (nvalue(index1(mn)), mvalue(index1(mn)),
     1   borbi1(index1(mn)), mn = 1,mnmax)
 20   FORMAT(1p,4(2i5,5x,e12.4,5x))
      WRITE (ioout, '(1x,131a)') dashes

c  output numerical parameters and Fourier spectrum

      WRITE (ioout, 30) mpol, ntor, lalpha, mpnt, meshtz, idisk,
     1   ipmb, dashes
 30   FORMAT(8x,'MPOL',8x,'NTOR',6x,'LALPHA',2x,'BLOCK SIZE',
     1   6x,'MESHTZ',7x,'IDISK',9x,'IPMB'/7i12/1x,131a/4x,'N')
      DO n = 1, ntor
         mm = mmnn(2,n)
         WRITE (ioout, 40) mmnn(1,n), (mmnn(2+m,n),m=1,mm)
      END DO
 40   FORMAT(i5,'    M =',20i6/12x,20i6/12x,20i6)
      WRITE (ioout, '(1x,131a)') dashes
      END SUBROUTINE header
