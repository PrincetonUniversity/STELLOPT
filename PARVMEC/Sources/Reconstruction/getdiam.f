      SUBROUTINE getdiam(amat_i, amat_p, data_array, kcdiam)
      USE vmec_main
      USE vmec_params, ONLY: signgs
      USE realspace
      USE vforces, ONLY : r12=>armn_o, ru12=>azmn_e
      USE vsvd
      USE vspline
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER kcdiam
      REAL(rprec), DIMENSION(isnodes,*) :: amat_i
      REAL(rprec), DIMENSION(ipnodes,*) :: amat_p
      REAL(rprec), DIMENSION(*) :: data_array
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: ilimit = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, lk, l
      REAL(rprec), DIMENSION(ns) :: gp, gi, gip
      REAL(rprec), DIMENSION(isnodes) :: amat2_i
      REAL(rprec) :: wdiam, z12, tv, ti, t2, sum1
C-----------------------------------------------


      kcdiam = 0
      IF (iphidiam.eq.0 .or. iresidue.lt.ilimit) RETURN
!
!       COMPUTE FIT TO DIAMAGNETIC SIGNAL, USING EQUILIBRIUM RELATION
!       (modified 7/96 by SPH)
!
!       PHI-DIAMAG = 2*pi*INT[ Gp dp/ds + Gi d(<Bu>)/ds ]
!
!       WHERE
!
!       Gp = Ru * Z * <SQRT(g)> /(R * phip)
!       Gi = Ru * Z * iota / R
!

      kcdiam = kcdiam + 1
c-7/96  dNewPhiedge = signgs*twopi*hs*Ssum_1(ns1,phip(2),1)
c-7/96  VacPhiedge  = signgs*bsubvvac*hs*Ssum_1(ns1,vrm2(2),1)
c-7/96  delphid0    = VacPhiedge - dNewPhiedge

      wdiam = one/sigma_delphid
      gp(1) = zero
      gi(1) = zero
      DO js = 2, ns
         gp(js) = zero
         DO lk = 1, nznt
            l = js + ns*(lk - 1)
            z12 = .5_dp*(z1(l,0)+z1(l-1,0)+shalf(l)*(z1(l,1)+z1(l-1,1)))
            gp(js) = gp(js) + ru12(l)*z12/r12(l)*wint(l)
         END DO
      END DO
!
!       NOTE: gip terms comes from linearizing the iota*d/ds[current*iota]
!             terms
!
      DO js = 2, ns
         tv = twopi*vp(js)/phip(js)
         ti = -gp(js)*signgs*wdiam
         gi(js) = ti*iotas(js)
         gp(js) = -gp(js)*tv*wdiam
         gip(js) = ti*(current(js)*iotaf(js)-current(js-1)*iotaf(js-1))
      END DO

      CALL splinint (gi, current, amat_i(1,kcdiam), hstark, u_ib, u1_ib
     1   , w_ib, w1_ib, nk_ib, isnodes, intder, ns)

      CALL splinint (gip(2), current(2), amat2_i, hstark, u_ia, u1_ia,
     1   w_ia, w1_ia, nk_ia, isnodes, intfun, ns1)

      CALL splinint (gp, presint, amat_p(1,kcdiam), hthom, u_pb, u1_pb,
     1   w_pb, w1_pb, nk_pb, ipnodes, intder, ns)

      amat_i(:isnodes,kcdiam) = amat_i(:isnodes,kcdiam) + amat2_i(:
     1   isnodes)

      t2 = mu0*pthommax                        !!*pfac moved to getthom
      amat_p(:ipnodes,kcdiam) = t2*amat_p(:ipnodes,kcdiam)
      sum1 = SUM(iotas(2:ns)*gip(2:ns))
      data_array(kcdiam) = wdiam*phidiam + sum1
      IF (iequi .eq. 0) THEN
!
!       Eliminate p variation until well-converged
!
!@        DO i = 1,ipnodes
!@          data_array(kcdiam) = data_array(kcdiam) -
!@     >    amat_p(i,kcdiam)*ythom(i)
!@          amat_p(i,kcdiam) = 0.
!@        END DO

!
!       FINAL OUTPUT (ALSO USE FOR DEBUGGING)
!
      ELSE
!
!       Integrate by parts
!
         delphid = gp(ns)*presf(ns) + gi(ns)*iotaf(ns)*current(ns) -
     1      gp(2)*presf(1) - gi(2)*current(1)*iotaf(1)
         DO js = 2, ns1
            delphid = delphid - presf(js)*(gp(js+1)-gp(js)) - iotaf(js)*
     1         current(js)*(gi(js+1)-gi(js))
         END DO
         delphid = delphid/wdiam
      ENDIF
!@        DO js = 2,ns
!@        END DO
!@
!@        sumi = SUM(amat_i(:isnodes,kcdiam)*ystark(:isnodes))
!@     >       - SUM(amat2_i(isnodes)*ystark(:isnodes))
!@        sump = bla_dot(ipnodes,amat_p(1,kcdiam),1,ythom,1)
!@
!@        WRITE(*,1212)delphid,(sumi+sump)/wdiam,delphid0
!@ 1212   FORMAT(' DelPhid = ',1pe10.3,' PhiD = ',1pe10.3,
!@     >   ' DelPhid0 = ',1pe10.3)

      END SUBROUTINE getdiam
