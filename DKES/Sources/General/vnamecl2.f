      MODULE Vnamecl2
      USE Vcmain, ONLY: lsource, rprec, dp
      INTEGER :: lap1, lam1, lam3, lam6, iswpm, ier
      REAL(rprec), PARAMETER :: zero = 0, one = 1, two = 2
      CHARACTER*3, DIMENSION(2), PARAMETER :: blabl = (/'  B','1/B'/)

      REAL(rprec) :: p2, root2, rthalf, rt3o2, rt5o2
      REAL(rprec) :: b00, bsqav
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1   cols, al1, al2, al3, al4, bl1,
     2   bl2, bl3, bl4, cl1, cl2, cl3, cl4
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1   cols0, omgl, al01, al02, al03, al04, bl01,
     2   bl02, bl03, bl04, cl01, cl02, cl03, cl04
      REAL(rprec) :: cmul1, efield1, wtov, weov, wcyclo, vthermi, vp,
     1   bpfac, s1cs10, s1cs1,
     2   rsd1p, rsd1m, rsd3p, rsd3m, g11p, g11m, g33p, g33m, g31p,
     3   g33s, g31m, g13p, g13m, crs1p, crs3p, crs1m, crs3m

      END MODULE Vnamecl2
