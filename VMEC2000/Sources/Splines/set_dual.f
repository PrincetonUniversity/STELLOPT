      SUBROUTINE set_dual(data1, hi, yi, yi2, hp, yp, yp2, wten, alsq,
     1   niota, npres, nots)
      USE vspline
      USE vmec_input, ONLY: tensi, tensi2, tensp, fpolyi
      USE vparams, ONLY: zero
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: niota, npres, nots, info
      REAL(rprec), DIMENSION(nots) :: data1
      REAL(rprec), DIMENSION(niota) :: hi, yi, yi2
      REAL(rprec), DIMENSION(npres) :: hp, yp, yp2
      REAL(rprec), DIMENSION(nots) :: wten
      REAL(rprec), DIMENSION(nots,nots) :: alsq
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER nb, ioff
C-----------------------------------------------
!
!       ADD TENSION TO DIAGONAL BLOCKS
!
      nb = ideriv
      ioff = 0
      CALL add_tension (alsq, wten, hi, tensi, tensi2, fpolyi, niota, nb
     1   , ioff, nots)
      CALL add_tension (alsq, wten, hp, tensp, zero, zero, npres, nb,
     1   niota, nots)

!
!       FREEZE EDGE PRESSURE IF NO PRESSURE SPECIED
!       COMPUTE SOLUTION FOR SPLINES
!
      CALL solver (alsq, data1, nots, 1, info)
      yi(:niota) = data1(:niota)
      yp(:npres) = data1(1+niota:npres+niota)

!
!       COMPUTE SECOND DERIVATIVES
!
      CALL gety2 (yi, yi2, hi, niota, nb)
      CALL gety2 (yp, yp2, hp, npres, nb)

      END SUBROUTINE set_dual
