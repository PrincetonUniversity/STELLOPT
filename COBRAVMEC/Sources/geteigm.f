      SUBROUTINE geteigm(ad, asup, asub, n, eigm, eigf)
      USE stel_kinds
      USE ballooning_data
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: n
      REAL(rprec), INTENT(IN), DIMENSION(n) :: ad
      REAL(rprec), INTENT(IN), DIMENSION(n-1) ::asup, asub
      REAL(rprec), INTENT(OUT) :: eigm
      REAL(rprec), INTENT(OUT),DIMENSION(n) :: eigf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: j
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: w, w2
      REAL(rprec):: MAXm, MINm, norm
!-------------------------------------------------------------------------------------
      ALLOCATE (w(n), w2(4*n), stat = j)
      IF (j .ne. 0) STOP 'Allocation error in COBRA geteigm'

      CALL TVAL(eigm, -kth, asub, ad, asup, n, w)
      CALL TVECT(eigm, eigf, asub, ad, asup, n, w2)

      DEALLOCATE (w, w2)

      MAXm=MAXVAL(eigf)
      MINm=MINVAL(eigf)
      norm=MAX(ABS(maxm),ABS(minm))

      IF(norm.eq.ABS(maxm))norm=maxm
      IF(norm.eq.ABS(minm))norm=minm
      eigf=eigf/norm                                               ! normalize eigenfunction to unity at maximum
      eigm=-eigm                                                   ! convert to Sturm-Liouville eigenvalue

      END SUBROUTINE geteigm
