      MODULE vsvd
      USE vmec_input
C      USE vmec_input, indata1 => indata, mseprofile1 => mseprofile
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, DIMENSION(:), ALLOCATABLE :: indexr, imid
      INTEGER, DIMENSION(:), ALLOCATABLE :: nk_ia, nk_pa, nk_ib,
     1   nk_pb, indexs2, indexu2, indexs1, indexu1, isortr, isorts
      INTEGER :: nmeasurements, imovephi,
     1   imse2, icurrout, islope, itse2, iphidiam, ipresin, ipresout,
     2   nchistp
      INTEGER, DIMENSION(mstp) :: nchi2
      INTEGER :: nchisaddle, nchitot, nchidia, nchipres, nchimse
      REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET ::
     1      datamse, qmid, shear, presmid, alfa, curmid,
     2      curint, psimid, ageo, volpsi, phimid
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1   current, rm2, vrm2, ovrm2, ochip, presph,
     2   presint, w_ia, w1_ia, u_ia, u1_ia, w_pa, w1_pa, u_pa, u1_pa,
     3   w_ib, w1_ib, u_ib, u1_ib, w_pb, w1_pb, u_pb, u1_pb, rmid,
     4   isplinef, isplineh, psplinef, psplineh, sthom, delse2, delso2,
     5   pcalc, delse1, delso1, starkcal, qmeas, qcalc, fpsical,
     6   stark_weight, rsort, rsort0
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: pm, im
      REAL(rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: pmb, imb
      REAL(rprec), DIMENSION(jngrn) ::
     1   yf, dyf, qsq, yek, yeq, dyek, dyeq
      REAL(rprec) :: odqk2, rstarkmin, rstarkmax, torflux, ppeak,
     1   pfac, phifac, phifsave, rthompeak, pthommax, rthommax,
     2   rthommin, delphid, dlim_min, rlim_min, zlim_min, router,
     3   rinner, apres, aminor, grmse, gphifac, rstepx0,
     4   rsfac, raxmse, rwidth, errsvd
      REAL(rprec), DIMENSION(jchix1) :: chisqerr
      REAL(rprec), DIMENSION(jchix1,mstp) :: chi2
      REAL(rprec) :: total_chi_square_n, total_chisq_n0,
     1   total_chi_square, total_saddle_chi, total_b_chi,
     2   total_pres_chi, total_mse_chi, total_chi_cur, total_chi_dia,
     3   scstark, scthom, flmwgt, bcwgt, tswgt, msewgt
      LOGICAL :: lpprof
      END MODULE vsvd
