      SUBROUTINE analyt(grpmn, bvec, ivacskip, ndim)
      USE vacmod
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: ivacskip, ndim
      REAL(rprec), INTENT(out) :: grpmn(nuv2*mnpd2)
      REAL(rprec), INTENT(out) :: bvec(mnpd2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: l, n, m
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1   r0p, r1p, r0m, r1m, sqrtc, sqrta, tlp2, tlp1, tlp, tlm2,
     2    tlm1, tlm, adp, adm, cma, ra1p, ra1m, slm, slp, tlpm, slpm,
     3    delt1u, azp1u, azm1u, cma11u, sqad1u, sqad2u
      REAL(rprec) :: fl, fl1, sign1
!-----------------------------------------------
      ALLOCATE (r0p(nuv2), r1p(nuv2), r0m(nuv2), r1m(nuv2),
     1          sqrtc(nuv2), sqrta(nuv2), tlp2(nuv2), tlp1(nuv2),
     2          tlp(nuv2), tlm2(nuv2), tlm1(nuv2), tlm(nuv2), adp(nuv2),
     3          adm(nuv2), cma(nuv2), ra1p(nuv2), ra1m(nuv2), slm(nuv2),
     4          slp(nuv2), tlpm(nuv2), slpm(nuv2), delt1u(nuv2),
     5          azp1u(nuv2), azm1u(nuv2), cma11u(nuv2), sqad1u(nuv2),
     6          sqad2u(nuv2), stat = l)
      IF (l .ne. 0) STOP 'Allocation error in SUBROUTINE analyt'

!
!     ALL EQUATIONS REFER TO THE PAPER BY P. MERKEL (PKM)
!     IN J. COMPUT. PHYSICS 66, p83 (1986)
!
!     IN GOING BETWEEN THE COMPLEX NOTATION OF (PKM) AND OUR REAL FORM,
!     NOTE THAT THE INTEGRALS (APPENDIX, PKM) Imn AND Kmn ARE BOTH REAL.
!     THUS, THE SIN(mu-nv) INTEGRALS OF THE SINGULAR PIECE (ANALYTIC CONTRIBUTION)
!     VANISHES.
!
!     THE REQUIRED SOURCE-TERM INTEGRALS ARE (Eq.2.16-2.17):
!
!     BVECS(m,n) = Int< SIN(mu' - nv') han(u',v') >
!     BVECC(m,n) = Int< COS(mu' - nv') han(u',v') >
!
!     Where Int<...> means integration over u (theta) and v (zeta) and
!     summation over field periods. These can be written in terms of PKM integrals
!     Imn(a,b,c), where a(u,v) = guu (g theta-theta), etc.:
!
!     BVECS(m,n) = ALP * Int<SIN(mu' - nv') * F * Im,-n(a,b,c)>
!     BVECC(m,n) = ALP * Int<COS(mu' - nv') * F * Im,-n(a,b,c)>
!
!     Here, F = - BNORM(u',v') is defined in Eq.(2.13), and ALP = (2*pi/nfp).
!
!     Similarly, the analytic part of the matrix A(m,n;m',n') can be written:
!
!     A(m,n;m',n') = (2*pi/nfp) * Int<SIN(mu' - nv')*SIN(m'u' - n'v')
!                              [Km,-n](a',b',c';A',B',C')>
!
!     On EXIT, GRPMN(ip,m,n) = ALP * SIN(ip,m,n) * K[m,-n](ip)
!
!
!     COMPUTE ALL QUANTITIES INDEPENDENT OF THE MODE INDICES L,M,N
!     NOTE: 2b = guv_b HAS FACTOR OF 2 BUILT IN (see SUBROUTINE SURFACE)
!
!     ADP(M): a +(-)2b + c
!     CMA:    c - a
!     DELTA:  4*(ac - b**2)
!     AZP(M): A +(-)2*B + C
!     CMA1:   C - A
!     R1P(M): Coefficient of l*Tl+(-) in eq (A17)
!     R0P(M): Coefficient of l*T(l-1)+(-) in eq (A17)
!     RA1P(M):Coefficient of Tl+(-) in eq (A17)
!
      adp  = guu_b  + guv_b  + gvv_b
      adm  = guu_b  - guv_b  + gvv_b
      cma  = gvv_b  - guu_b
      sqrtc  = two*SQRT(gvv_b)
      sqrta  = two*SQRT(guu_b)

      IF (ivacskip .eq. 0) THEN
         grpmn = 0
         delt1u  = adp*adm  - cma*cma
         azp1u  = auu  + auv  + avv
         azm1u  = auu  - auv  + avv
         cma11u  = avv  - auu
         r1p  = (azp1u*(delt1u - cma*cma)/adp
     1        -  azm1u*adp + two*cma11u*cma)/delt1u
         r1m  = (azm1u*(delt1u - cma*cma)/adm
     1        -  azp1u*adm + two*cma11u*cma)/delt1u
         r0p  = (-azp1u*adm*cma/adp - azm1u*cma
     1        + two*cma11u*adm)/delt1u
         r0m  = (-azm1u*adp*cma/adm - azp1u*cma
     1        + two*cma11u*adp)/delt1u
         ra1p = azp1u/adp
         ra1m = azm1u/adm
      ENDIF

!
!     INITIALIZE VECTORS
!
      bvec = 0
!
!     INITIALIZE T0+ and T0-
!
!     TLP(M): TL+(-)
!     TLP(M)1:T(L-1)+(-)
!     TLP(M)2:T(L-2)+(-)
!
      sqad1u = SQRT(adp)
      sqad2u = SQRT(adm)
      tlp1 = 0
      tlm1 = 0
      tlp  = one/sqad1u*log((sqad1u*sqrtc + adp + cma)
     1                     /(sqad1u*sqrta - adp + cma))
      tlm  = one/sqad2u*log((sqad2u*sqrtc + adm + cma)
     1                     /(sqad2u*sqrta - adm + cma ))
      tlpm = tlp  + tlm
!
!     BEGIN L-SUM IN EQ (A14) TO COMPUTE Imn (and Kmn) INTEGRALS
!     NOTE THAT IN THE LOOP OVER L BELOW: L == |m - n| + 2L_A14
!     THUS, L BELOW IS THE INDEX OF THE T+- (S+-)
!
      sign1 = 1
      fl1 = 0

      LLOOP: DO l = 0, mf + nf
         fl = fl1
!
!     COMPUTE SL+ and SL- , Eq (A17)
!     SLP(M): SL+(-)
!
         IF (ivacskip .eq. 0) THEN
            slp = (r1p*fl + ra1p)*tlp + r0p*fl*tlp1 - (r1p + r0p)/sqrtc
     1          + sign1*(r0p - r1p)/sqrta
            slm = (r1m*fl + ra1m)*tlm + r0m*fl*tlm1 - (r1m + r0m)/sqrtc
     1          + sign1*(r0m - r1m)/sqrta
            slpm = slp + slm
         ENDIF
!
!     BEGIN MODE NUMBER (m,n) LOOP
!
         DO n = 0, nf
            DO m = 0, mf
               IF (cmns(l,m,n) .eq. zero) CYCLE

               IF (n.eq.0 .or. m.eq.0) THEN
!
!       1. n = 0 and  m >= 0  OR n > 0 and m = 0
!
                 CALL analysum (grpmn, bvec, slpm, tlpm, m, n, l,
     1               ivacskip, ndim)

               ELSE
!
!       2. n>=1  and  m>=1
!
                 CALL analysum2 (grpmn, bvec, slm, tlm, slp, tlp,
     1               m, n, l, ivacskip, ndim)
               ENDIF
            END DO
         END DO
!
!     UPDATE "TL's" (FOR L -> L+1) USING EQ (A15)
!
         tlp2 = tlp1
         tlm2 = tlm1
         tlp1 = tlp
         tlm1 = tlm
         fl1  = fl1 + 1                       !next l
         sign1 = -sign1                       !(-1)**l (next l now)
         tlp = ((sqrtc + sign1*sqrta) - (two*fl1 - one)*cma*tlp1 -
     1      fl*adm*tlp2)/(adp*fl1)
         tlm = ((sqrtc + sign1*sqrta) - (two*fl1 - one)*cma*tlm1 -
     1      fl*adp*tlm2)/(adm*fl1)
         tlpm = tlp + tlm

      END DO LLOOP

      DEALLOCATE (r0p, r1p, r0m, r1m, sqrtc, sqrta, tlp2, tlp1,
     1          tlp, tlm2, tlm1, tlm, adp, adm, cma, ra1p, ra1m, slm,
     2          slp, tlpm, slpm, delt1u, azp1u, azm1u, cma11u, sqad1u,
     3          sqad2u, stat = l)
      END SUBROUTINE analyt
