      SUBROUTINE analyt (grpmn, bvec, ivacskip, ndim)
      USE vacmod
      USE parallel_include_module
      USE timer_sub
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: ivacskip, ndim
      REAL(dp), INTENT(OUT) :: grpmn(mnpd2,nuv3)
      REAL(dp), INTENT(OUT) :: bvec(mnpd,ndim)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l, n, m, i, q, j, k, ll, blksize, mn
      REAL(dp), DIMENSION(:), ALLOCATABLE ::
     1   r0p, r1p, r0m, r1m, sqrtc, sqrta, tlp2, tlp1, tlp, tlm2,
     2    tlm1, tlm, adp, adm, cma, ra1p, ra1m, slm, slp, tlpm, slpm,
     3    delt1u, azp1u, azm1u, cma11u, sqad1u, sqad2u
      REAL(dp) :: fl, fl1, fl2, sign1, tanalon, tanaloff
C-----------------------------------------------
      CALL second0(tanalon)

      ALLOCATE (r0p(nuv3), r1p(nuv3), r0m(nuv3), r1m(nuv3),
     1          sqrtc(nuv3), sqrta(nuv3), tlp2(nuv3), tlp1(nuv3),
     2          tlp(nuv3), tlm2(nuv3), tlm1(nuv3), tlm(nuv3), adp(nuv3),
     3          adm(nuv3), cma(nuv3), ra1p(nuv3), ra1m(nuv3), slm(nuv3),
     4          slp(nuv3), tlpm(nuv3), slpm(nuv3), delt1u(nuv3),
     5          azp1u(nuv3), azm1u(nuv3), cma11u(nuv3), sqad1u(nuv3),
     6          sqad2u(nuv3), stat = l)
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
      DO k = nuv3min, nuv3max
         adp(k) = guu_b(k) + guv_b(k) + gvv_b(k) 
         adm(k) = guu_b(k) - guv_b(k) + gvv_b(k) 
         cma(k) = gvv_b(k) - guu_b(k) 
         sqrtc(k)   = two*SQRT(gvv_b(k))
         sqrta(k)   = two*SQRT(guu_b(k))
      END DO

      IF (ivacskip .EQ. 0) THEN

         grpmn(:,nuv3min:nuv3max) = 0

         DO k = nuv3min, nuv3max
         delt1u(k) = adp(k)*adm(k) - cma(k)*cma(k)
         azp1u(k) = auu(k) + auv(k) + avv(k)
         azm1u(k) = auu(k) - auv(k) + avv(k)
         cma11u(k)= avv(k) - auu(k)
         r1p(k) = (azp1u(k)*(delt1u(k) - cma(k)*cma(k))/adp(k)
     1          -  azm1u(k)*adp(k) + two*cma11u(k)*cma(k))/delt1u(k)
         r1m(k) = (azm1u(k)* (delt1u(k) - cma(k)*cma(k))/adm(k)
     1          -  azp1u(k)*adm(k) + two*cma11u(k)*cma(k))/delt1u(k)
         r0p(k) = (-azp1u(k)*adm(k)*cma(k)/adp(k) - azm1u(k)*cma(k)
     1          +  two*cma11u(k)*adm(k))/delt1u(k)
         r0m(k)  = (-azm1u(k)*adp(k)*cma(k)/adm(k) - azp1u(k)*cma(k)
     1           + two*cma11u(k)*adp(k))/delt1u(k)
         ra1p(k) = azp1u(k)/adp(k)
         ra1m(k) = azm1u(k)/adm(k)
         END DO
      ENDIF

!
!     INITIALIZE VECTORS
!
!      bvec = 0
!
!     INITIALIZE T0+ and T0-
!
!     TLP(M): TL+(-)
!     TLP(M)1:T(L-1)+(-)
!     TLP(M)2:T(L-2)+(-)
!
      DO k = nuv3min,nuv3max
      sqad1u(k) = SQRT(adp(k))
      sqad2u(k) = SQRT(adm(k))
      tlp1(k) = 0
      tlm1(k) = 0
      tlp(k)  = one/sqad1u(k)*log((sqad1u(k)*sqrtc(k) 
     1        + adp(k) + cma(k))/(sqad1u(k)*sqrta(k) 
     2        - adp(k) + cma(k)))
      tlm(k)  = one/sqad2u(k)*log((sqad2u(k)*sqrtc(k) 
     1        + adm(k) + cma(k))/(sqad2u(k)*sqrta(k) 
     2        - adm(k) + cma(k)))
      tlpm(k) = tlp(k) + tlm(k)
      END DO
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
            DO k = nuv3min,nuv3max
            slp(k) = (r1p(k)*fl + ra1p(k))*tlp(k) + r0p(k)*fl*tlp1(k)
     1             - (r1p(k) + r0p(k))/sqrtc(k) 
     2             + sign1*(r0p(k) - r1p(k))/sqrta(k)
            slm(k) = (r1m(k)*fl + ra1m(k))*tlm(k) + r0m(k)*fl*tlm1(k)
     1             - (r1m(k) + r0m(k))/sqrtc(k) 
     2             + sign1*(r0m(k) - r1m(k))/sqrta(k)
            slpm(k) = slp(k) + slm(k)
            END DO
         ENDIF
!
!     BEGIN MODE NUMBER (m,n) LOOP
!
         DO n = 0, nf
            DO m = 0, mf

               IF (l .EQ. 0) THEN
                  mn = m + mf1*(n+nf) + 1
                  bvec(mn,:) = 0
                  mn = m + mf1*(nf-n) + 1
                  bvec(mn,:) = 0
               END IF

               IF (cmns(l,m,n) .eq. zero) CYCLE
               
               IF (n.eq.0 .or. m.eq.0) THEN
!
!       1. n = 0 and  m >= 0  OR n > 0 and m = 0
!
                 CALL analysum (grpmn, bvec, slpm, tlpm, m, n, l,
     1                          ivacskip, ndim)

               ELSE
!
!       2. n>=1  and  m>=1
!
                 CALL analysum2 (grpmn, bvec, slm, tlm, slp, tlp,
     1                           m, n, l, ivacskip, ndim)

               ENDIF
            END DO
         END DO

!
!     UPDATE "TL's" (FOR L -> L+1) USING EQ (A15)
!
         fl1  = fl1 + 1                       !next l
         fl2  = 2*fl1-1
         sign1 = -sign1                       !(-1)**l (next l now)
         DO k = nuv3min, nuv3max
         tlp2(k) = tlp1(k)
         tlm2(k) = tlm1(k)
         tlp1(k) = tlp(k)
         tlm1(k) = tlm(k)
         tlp(k) = ((sqrtc(k) + sign1*sqrta(k)) - fl2*
     1       cma(k)*tlp1(k) - fl*adm(k)*tlp2(k))/(adp(k)*fl1)
         tlm(k) = ((sqrtc(k) + sign1*sqrta(k)) - fl2*
     1       cma(k)*tlm1(k) - fl*adp(k)*tlm2(k))/(adm(k)*fl1)
         tlpm(k) = tlp(k) + tlm(k)
         END DO

      END DO LLOOP

      DEALLOCATE (r0p, r1p, r0m, r1m, sqrtc, sqrta, tlp2, tlp1,
     1          tlp, tlm2, tlm1, tlm, adp, adm, cma, ra1p, ra1m, slm,
     2          slp, tlpm, slpm, delt1u, azp1u, azm1u, cma11u, sqad1u,
     3          sqad2u, stat = l)

      CALL second0(tanaloff)
      timer_vac(tanal) = timer_vac(tanal) + (tanaloff-tanalon)
      analyt_time = timer_vac(tanal)


      END SUBROUTINE analyt

