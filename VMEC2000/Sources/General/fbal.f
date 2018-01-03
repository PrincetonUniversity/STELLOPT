      MODULE fbal
      USE stel_kinds, ONLY: dp
      REAL(dp), DIMENSION(:), ALLOCATABLE :: rzu_fac, rru_fac,
     1  frcc_fac, fzsc_fac

      CONTAINS

#if defined(SKS)

      SUBROUTINE calc_fbal_par(bsubu, bsubv)
      USE vmec_main, ONLY: buco, bvco, equif, iequi,
     1                     jcurv, jcuru, chipf, vp, pres, 
     2                     phipf, vpphi, presgrad, ohs
      USE vmec_params, ONLY: signgs
      USE vmec_dim, ONLY: ns, nrzt, nznt, ns1
      USE realspace, ONLY: pwint, phip
      USE vmec_input, ONLY: nzeta
      USE vmec_dim, ONLY: ntheta3
      USE parallel_include_module 
      IMPLICIT NONE
!-----------------------------------------------
      REAL(dp), INTENT(in) :: bsubu(nznt,ns),
     1                        bsubv(nznt,ns)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: js, lk
      INTEGER :: nsmin, nsmax
!-----------------------------------------------

      nsmin=t1lglob; nsmax=t1rglob
      DO js = nsmin, nsmax
        buco(js) = SUM(bsubu(:,js)*pwint(:,js))
        bvco(js) = SUM(bsubv(:,js)*pwint(:,js))
      END DO

      CALL Gather1XArray(bvco)
      CALL Gather1XArray(buco)

!     FROM AMPERE'S LAW, JcurX are angle averages of jac*JsupX, so
!                        JcurX = (dV/ds)/twopi**2 <JsupX> where <...> is flux surface average
      !nsmin=MAX(2,t1lglob); nsmax=MIN(t1rglob,ns-1)
      nsmin=MAX(2,tlglob); nsmax=MIN(trglob,ns-1)
      DO js = nsmin, nsmax
         jcurv(js) = (signgs*ohs)*(buco(js+1) - buco(js))
         jcuru(js) =-(signgs*ohs)*(bvco(js+1) - bvco(js))
         vpphi(js) = (vp(js+1) + vp(js))/2
         presgrad(js) = (pres(js+1) - pres(js))*ohs
         equif(js) = (-phipf(js)*jcuru(js) + chipf(js)*jcurv(js))
     1                /vpphi(js) + presgrad(js)
      END DO
      equif(1) = 0
      equif(ns) = 0

      !SKS-RANGE: All LHS's computed correctly in [t1lglob, trglob]

      END SUBROUTINE calc_fbal_par
#endif

      SUBROUTINE calc_fbal(bsubu, bsubv)
      USE vmec_main, ONLY: buco, bvco, equif, 
     1                     jcurv, jcuru, chipf, vp, pres, 
     2                     phipf, vpphi, presgrad, ohs
#ifdef _ANIMEC
     3                    ,pmap, pd, phot, tpotb, zero
#endif
      USE vmec_params, ONLY: signgs
      USE vmec_dim, ONLY: ns, nrzt, nznt, ns1
      USE realspace, ONLY: wint, phip
#ifdef _ANIMEC
     1                    ,pperp, ppar, onembc, sigma_an,
     2                     pp1, pp2, pp3
      USE vforces, gsqrt => azmn_o
#endif
#if defined (SKS)      
      USE parallel_include_module 
#endif      
      IMPLICIT NONE
!-----------------------------------------------
      REAL(dp), INTENT(in) :: bsubu(1:nrzt), bsubv(1:nrzt)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: js, lk
      INTEGER :: nsmin, nsmax
#ifdef _ANIMEC
      REAL(dp) :: t4, t5
#endif
!-----------------------------------------------
      DO js = 2, ns
         buco(js) = SUM(bsubu(js:nrzt:ns)*wint(js:nrzt:ns))
         bvco(js) = SUM(bsubv(js:nrzt:ns)*wint(js:nrzt:ns))
      END DO

!     FROM AMPERE'S LAW, JcurX are angle averages of jac*JsupX, so
!                        JcurX = (dV/ds)/twopi**2 <JsupX> where <...> is flux surface average
      DO js = 2, ns1
         jcurv(js) = (signgs*ohs)*(buco(js+1) - buco(js))
         jcuru(js) =-(signgs*ohs)*(bvco(js+1) - bvco(js))
!FOR RFP vpphi(js) = (vp(js+1)/phip(js+1) + vp(js)/phip(js))/2
         vpphi(js) = (vp(js+1) + vp(js))/2
         presgrad(js) = (pres(js+1) - pres(js))*ohs
         equif(js) = (-phipf(js)*jcuru(js) + chipf(js)*jcurv(js))
     1                /vpphi(js) + presgrad(js)
      END DO
      equif(1) = 0
      equif(ns) = 0

      END SUBROUTINE calc_fbal

#ifdef _ANIMEC
      SUBROUTINE bimax_ppargrad(pp1, pp2, pp3, ppar, onembc, pres, phot,
     1                          tpotb)
      USE vmec_main, ONLY: zero
      USE vmec_dim, ONLY: ns, nrzt, nznt, ns1
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), INTENT(in) :: pp3(nrzt), ppar(nrzt), phot(ns),
     1                        onembc(nrzt) ,pres(ns), tpotb(ns)
      REAL(dp), INTENT(out):: pp1(nrzt), pp2(nrzt)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER     :: js, lk, l
      REAL(dp), ALLOCATABLE :: tmp0(:), tmp2(:)
      REAL(dp) :: eps
C-----------------------------------------------
!********0*********0*********0*********0*********0*********0*********0**
!     Model with Anisotropic Pressure from Bi-Maxwellian Distribution. *
!     EVALUATION OF AMPLITUDES OF partial p_parallel/partial s AT FIXED B.
!     Compute hot particle parallel pressure gradient at fixed B times PP3
!     PP3 can be sqrt(g), sqrt(g)/V' or unity
!*********0*********0*********0*********0*********0*********0*********0*
      eps = EPSILON(eps)
      ALLOCATE (tmp0(nrzt), tmp2(nrzt))

      pp1(1:nrzt:ns) = 0;  pp2(1:nrzt:ns) = 0
      DO js = 2, ns
         tmp0(js:nrzt:ns) = pp3(js:nrzt:ns)*(ppar(js:nrzt:ns)-pres(js))
         tmp2(js:nrzt:ns) = tmp0(js:nrzt:ns)/(pres(js)*phot(js)+eps)
      END DO

      DO l = 1,nrzt-1
         pp2(l) = 0.5_dp*(tmp2(l)+tmp2(l+1))
      END DO
      DO l = 1, nrzt
         tmp2(l) = pp3(l)*(1-onembc(l))
      END DO

      DO js = 2,ns
         DO lk = 1,nznt
            l = js +(lk-1)*ns
            IF (onembc(l) <= zero) THEN
               pp1(l)= tmp0(l)/
     &                (1 - tpotb(js)   * onembc(l)+eps)
            ELSE
               pp1(l)=(2*tpotb(js  )*tmp0(l)*onembc(l)
     &               + pres(js  ) * phot(js  ) * tmp2(l)*
     &             ( 1 - 5*(tpotb(js  )*onembc(l))**1.5_dp))
     &          /  (1 - (tpotb(js  )*onembc(l))**2 + eps)
            END IF
         END DO
      END DO

      pp1 = pp1 * onembc
      DO l=1, nrzt-1
         pp1(l) = 0.5_dp*(pp1(l)+pp1(l+1))
      END DO

      DEALLOCATE (tmp0, tmp2)

      END SUBROUTINE bimax_ppargrad

      SUBROUTINE mirror_crit(taumir, bsq)
      USE stel_kinds, ONLY: rprec, dp
      USE realspace, ONLY: sigma_an, pperp, ppar, onembc
!      USE vforces,   ONLY: bsq => bzmn_o
      USE vmec_main, ONLY: tpotb, pppr, pmap, pres, papr, vp, wb,
     &                     wp, wpar, wper, ns, nznt, nrzt,
     &                     zero, one, nthreed
!
!     WAC (11/30/07): See description of anisotropic pressure below
!
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(nrzt), INTENT(out)            :: taumir
      REAL(dp), DIMENSION(nrzt), INTENT(inout)          :: bsq
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER     :: js   , lk   , l
      REAL(dp) :: whpar, whper, bhotmx, taumin, ppprht, parpht,
     &        sigmin, bhotc, betath, betaht, betato, betapa,
     &        betape, betadi, betaew, es, betprc, betprx
C-----------------------------------------------
!     CALCULATION OF MIRROR STABILITY CRITERION
!********0*********0*********0*********0*********0*********0*********0**
!                                                                      *
!                 Anisotropic Pressure Model                           *
!                  specific to case where:                             *
!          a Bi-Maxwellian distribution is considered (by J. Graves)   *
!          p_parallel(s,B) = pth(1 + phot(s)*H(s,B))                     *
!      H(s,B)=(B/B_crit)/[1-(T_perp/T_par)(1-B/B_crit)] for B>B_crit   *
!      For B<B_crit,                                                   *
!      H(s,B)=H(s,B>B_crit){1-2[(T_perp/T_par)(1-B/B_crit)]^(5/2) /    *
!                                    [1+(T_perp/T_par)(1-B/B_crit)]}   *
!                                                                      *
!********0*********0*********0*********0*********0*********0*********0**
!
!********0*********0*********0*********0*********0*********0*********0**
!   1.  Change BSQ from total pressure to magnetic pressure            *
!********0*********0*********0*********0*********0*********0*********0**
c
      bsq = bsq - pperp
!
!********0*********0*********0*********0*********0*********0*********0**
!   2.  Compute Tau-minimum, Sigma-minimum, Peak Hot Particle Beta.    *
!********0*********0*********0*********0*********0*********0*********0**
c
      whpar = wpar - wp
      whper = wper - wp
      bhotmx = -HUGE(bhotmx)
      taumin =  HUGE(taumin)
      sigmin =  HUGE(sigmin)

      DO 57 js = 2,ns
         DO 55 lk = 1,nznt
            l = (lk-1) * ns + js
            taumir(l) = one + (pperp(l) - pres(js)) / bsq(l)*
     &      (one - tpotb(js)) / (one - tpotb(js) * onembc(l))
!              taut = taumir(l) + pperp(l)/bsq(l) *0.25*tpotb(js)
            if(onembc(l) .gt. zero)
     &        taumir(l) = taumir(l)
     &         + (pperp(l) - pres(js)) / bsq(l) * 0.25_dp*tpotb(js)
     &         * (one-onembc(l))*sqrt(tpotb(js)*onembc(l))
     &         / (one + tpotb(js) * onembc(l))
     &         * (15 - 5*tpotb(js) * onembc(l)
     &         - 7*(tpotb(js)*onembc(l))**2
     &         - 3*(tpotb(js)*onembc(l))**3 )
     &         / ((one + tpotb(js)*onembc(l))**2
     &         - (tpotb(js)*onembc(l))**1.5_dp * (5
     &         - (tpotb(js)*onembc(l))**2))
55       end do
57    end do
      do 62 js = 2,ns
         do 60 lk = 1,nznt
            l = (lk-1) * ns + js
            taumin = MIN(taumir(l),taumin)
            sigmin = MIN(sigma_an(l),sigmin)
            bhotc = (2*pperp(l)+ppar(l)-3*pres(js))/bsq(l)
            bhotmx = MAX(bhotc,bhotmx)
 60     end do
 62   end do
      betprx = HUGE(betprx)
      do 64 js=2,ns
         betprc = 0.76786580_dp - 0.29708540_dp * tpotb(js)
     &      + 0.054249860_dp * tpotb(js)**2
     &      - 0.0054254849_dp * tpotb(js)**3
     &      + 0.00030947525_dp * tpotb(js)**4
     &      - 9.7144781e-6_dp * tpotb(js)**5
     &      + 1.3844199e-7_dp * tpotb(js)**6
     &      - 1.4328673e-11_dp * (one + tpotb(js))*tpotb(js)**7
        betprx=min(betprx,betprc)
 64   end do
      bhotmx = bhotmx/3
      write (nthreed,101) taumin, sigmin, bhotmx
 101  format(" taumin=",1pe15.6," sigmin=",1pe15.6,
     &       "  peak hot beta=",1pe15.6)
      betath = wp/wb
      betaht = (2*whper+whpar)/(3*wb)
      betato = (2*wper+wpar)/(3*wb)
      betapa = whpar/wb
      betape = whper/wb
      betadi = betath + betape
      betaew = betath + 0.5*(betapa+betape)
      write (nthreed,102) betath,betaht,betato
 102  format(" thermal beta =",1pe13.4," beta-hot =",1pe13.4,
     &       " total beta =",1pe13.4)
      write (nthreed,103) betapa,betape,betaew-betath
 103  format(' hot parallel beta =',1pe13.4,
     &  ' hot perpendicular beta =',1pe13.4,' beta-hot(EW) =',1pe13.4)
      write (nthreed,104) betadi,betaew,betprx
 104  format(' diamagnetic beta =',1pe13.4,' beta(equal weighting) ='
     &     , 1pe13.4,' maximum hot perp. beta_c =',1pe10.3)

      IF (ABS(betaht) .GT. 1.E-12_dp) THEN
         write (nthreed,105)
 105  format(/,' js',8x,'s',9x,'th. press',3x,'par. pres',2x,
     &         'perp. pres')
      DO 79 js = 2,ns
         es = REAL(js - 1.5,dp) / (ns-1)
         ppprht = pppr(js) - pres(js)
         parpht = papr(js) - pres(js)
         WRITE (nthreed,106) js, es, pres(js), parpht, ppprht
 79   END DO
      END IF
 106  format (i3,1p,1e15.6,1p,3e12.3)

      END SUBROUTINE mirror_crit
#endif

      END MODULE fbal
