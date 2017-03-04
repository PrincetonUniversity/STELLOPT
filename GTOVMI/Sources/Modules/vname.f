
      MODULE Vname0
      INTEGER, PARAMETER :: mu = 20
      INTEGER, PARAMETER :: nu = 258
      INTEGER, PARAMETER :: nv = 20
      INTEGER, PARAMETER :: mnd = mu*(nv + 3)
      END MODULE Vname0
      MODULE Vname1
      USE precision
      USE Vname0
      REAL(rprec), PARAMETER :: one = 1, zero = 0
      INTEGER :: mpol, mpol1, mrho, nfp, ntheta, nphi, nphi2, n2, mpnt
      INTEGER :: mm(mu), nn(nv+1), m1(mnd), n1(mnd)
      REAL(rprec), DIMENSION(nv) :: r0n, z0n, raxis, zaxis
      REAL(rprec), DIMENSION(0:mu) :: dm1
      REAL(rprec), DIMENSION(0:mu+2) :: t1m, t2m
      REAL(rprec), DIMENSION(0:mu,2) :: xmpq
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xvec, gvec,
     1   xdot, xstore
      REAL(rprec) :: gnorm, specw, delt, deltf, dnorm, elongate,
     1   twopi, r10, HB_Parameter
      END MODULE Vname1

