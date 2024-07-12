      MODULE ga_mod
      USE gade_mod
      IMPLICIT NONE

      INTEGER, DIMENSION(nparmax) :: ig2
      INTEGER, DIMENSION(nchrmax,indmax) :: iparent, ichild
      INTEGER, DIMENSION(nchrmax) :: ibest
      REAL(rprec), DIMENSION(nparmax,indmax) :: parent, child
      REAL(rprec), DIMENSION(indmax) :: fitness
      REAL(rprec), DIMENSION(nparmax) :: g0, g1,
     +                     pardel, par_max, par_min
      REAL(rprec), DIMENSION(max_gen) :: geni, genavg, genmax
      INTEGER :: nparam, nchrome, num_obj
      INTEGER :: jbest,irestrt
      INTEGER :: MAXgen,nfit_eval,
     +                 kountmx
      INTEGER :: iunit_ga_restart, iunit_ga_out
      REAL(rprec), DIMENSION(:), POINTER :: f_obj

      END MODULE ga_mod
