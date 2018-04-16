      MODULE xstuff
      USE stel_kinds, ONLY: rprec
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1    gc, xcdot, xsave, xstore, scalxc, col_scale
      REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET ::
     1    xc

      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1   pgc, pxcdot, pxsave, pxstore, pscalxc, pcol_scale
      REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET ::
     1    pxc
C-----------------------------------------------
      END MODULE xstuff
