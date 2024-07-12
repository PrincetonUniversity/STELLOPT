      MODULE vmec_main
      USE vmec_dim
      USE vmec_input
      USE vmec_persistent
      USE vmec_params, ONLY: ndamp
      USE vparams
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(dp), DIMENSION(:,:), ALLOCATABLE ::
     1    ard, arm, brd, brm, azd, azm, bzd, bzm, bmin, bmax
      REAL(dp), DIMENSION(:), ALLOCATABLE ::
     1    crd, iotaf, phipf, chipf, mass, phi, presf, beta_vol, 
     2    jcuru, jcurv, jdotb, 
#ifdef _ANIMEC
!WAC: ANISOTROPIC ARRAYS
     2    phot, pmap, pppr, papr, tpotb, pd,
#endif
     2    buco, bvco, bdotgradv, equif, specw, tcon, 
     3    psi, yellip, yinden, ytrian, yshift, ygeo, overr, 
     4    sm, sp, iotas, phips, chips, pres, vp, jpar2, jperp2, bdotb, 
     5    blam, clam, dlam, icurv, vpphi, presgrad, 
     6    r01, z01, bdamp, bucof, bvcof, chi
      REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: faclam, faclam0
      REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: pfaclam
      REAL(dp), ALLOCATABLE :: xcl0(:)

      REAL(dp), DIMENSION(0:mpol1d,3) :: xmpq
      REAL(dp), DIMENSION(0:mpol1d) :: faccon
      REAL(dp) :: dcon, currv, aspect, hs, ohs, voli, 
     1   signiota, rc0mse, r00, r0scale, z00, dkappa, fsqsum0,
     2   pressum0, fnorm, fsqr=1, fsqz=1, fsql=1, fnorm1, fnorml, 
     3   fsqr1, fsqz1, fsql1, fsq, fedge, wb, wp, r00b, z00b, fz00_edge
#ifdef _ANIMEC
     4  ,wpar, wper
#endif
      REAL(dp), DIMENSION(nstore_seq) :: fsqt, wdot
      REAL(dp) :: ftolv, otav, alphaR, alphaZ
      REAL(dp), DIMENSION(ndamp) :: otau
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::
     1    rmn_bdy, zmn_bdy
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: bsqsav
      REAL(dp), DIMENSION(:), ALLOCATABLE :: bsubu0, dbsq, rbsq
#ifdef _ANIMEC
      REAL(dp), DIMENSION(:), ALLOCATABLE :: pperp_ns
      REAL(dp) :: medge, phedg
#endif
      REAL(dp) :: rbtor, rbtor0, ctor, delbsq, res0, res1, delt0r   !DO NOT remove res0, delt0r -> V3FIT
      REAL(dp), DIMENSION(ndatafmax) ::
     1  spfa, spfa2, hp, sifa, sifa2, hi
      LOGICAL :: lthreed, lconm1
      INTEGER, DIMENSION(:), ALLOCATABLE :: ireflect
      INTEGER :: multi_ns_grid, iequi, itype_precon, irst,
     1    iter1, iter2, iterc=1, ijacob, itfsq, iresidue, neqs,
     2    irzloff, ivac, ndatap, ndatai

      REAL(dp) :: router
      REAL(dp) :: rinner
C-----------------------------------------------
      END MODULE vmec_main
