      SUBROUTINE symforce(ars, brs, crs, azs, bzs, czs, bls, 
     1 cls, rcs, zcs, ara, bra, cra, aza, bza, cza, bla, cla, rca, zca)
      USE vmec_main, p5 => cp5
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1),
     1   INTENT(inout) :: ars, brs, crs, azs, bzs, czs,
     2   bls, cls, rcs, zcs
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(out) ::
     1   ara, bra, cra, aza, bza, cza, bla, cla, rca, zca
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mpar, ir, i, jk, jka
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: ars_0, brs_0, azs_0, 
     1              bzs_0, bls_0, rcs_0, zcs_0, crs_0, czs_0, cls_0
C-----------------------------------------------
      i = ns*nzeta
      ALLOCATE (ars_0(i), brs_0(i), azs_0(i), bzs_0(i), bls_0(i),
     1          rcs_0(i), zcs_0(i), crs_0(i), czs_0(i), cls_0(i),
     2          stat=ir)

!
!       SYMMETRIZE FORCES ON RESTRICTED THETA INTERVAL (0 <= u <= pi)
!       SO COS,SIN INTEGRALS CAN BE PERFORMED. FOR EXAMPLE,
!
!       ARS(v,u) = .5*( ARS(v,u) + ARS(-v,-u) )     ! * COS(mu - nv)
!       ARA(v,u) = .5*( ARS(v,u) - ARS(-v,-u) )     ! * SIN(mu - nv)
!
!
      DO mpar = 0, 1
         DO i = 1, ntheta2
            ir = ntheta1 + 2 - i                 !-theta
            IF (i .eq. 1) ir = 1
            DO jk = 1, ns*nzeta
               jka = ireflect(jk)                !-zeta
               ara(jk,i,mpar) = p5*(ars(jk,i,mpar)-ars(jka,ir,mpar))
               ars_0(jk)      = p5*(ars(jk,i,mpar)+ars(jka,ir,mpar))
               bra(jk,i,mpar) = p5*(brs(jk,i,mpar)+brs(jka,ir,mpar))
               brs_0(jk)      = p5*(brs(jk,i,mpar)-brs(jka,ir,mpar))
               aza(jk,i,mpar) = p5*(azs(jk,i,mpar)+azs(jka,ir,mpar))
               azs_0(jk)      = p5*(azs(jk,i,mpar)-azs(jka,ir,mpar))
               bza(jk,i,mpar) = p5*(bzs(jk,i,mpar)-bzs(jka,ir,mpar))
               bzs_0(jk)      = p5*(bzs(jk,i,mpar)+bzs(jka,ir,mpar))
               bla(jk,i,mpar) = p5*(bls(jk,i,mpar)-bls(jka,ir,mpar))
               bls_0(jk)      = p5*(bls(jk,i,mpar)+bls(jka,ir,mpar))
               rca(jk,i,mpar) = p5*(rcs(jk,i,mpar)-rcs(jka,ir,mpar))
               rcs_0(jk)      = p5*(rcs(jk,i,mpar)+rcs(jka,ir,mpar))
               zca(jk,i,mpar) = p5*(zcs(jk,i,mpar)+zcs(jka,ir,mpar))
               zcs_0(jk)      = p5*(zcs(jk,i,mpar)-zcs(jka,ir,mpar))
            END DO

            ars(:,i,mpar) = ars_0(:)
            brs(:,i,mpar) = brs_0(:)
            azs(:,i,mpar) = azs_0(:)
            bzs(:,i,mpar) = bzs_0(:)
            bls(:,i,mpar) = bls_0(:)
            rcs(:,i,mpar) = rcs_0(:)
            zcs(:,i,mpar) = zcs_0(:)

            IF (lthreed) THEN
               DO jk = 1, ns*nzeta
                  jka = ireflect(jk)
                  cra(jk,i,mpar) = p5*(crs(jk,i,mpar)+crs(jka,ir,mpar))
                  crs_0(jk)      = p5*(crs(jk,i,mpar)-crs(jka,ir,mpar))
                  cza(jk,i,mpar) = p5*(czs(jk,i,mpar)-czs(jka,ir,mpar))
                  czs_0(jk)      = p5*(czs(jk,i,mpar)+czs(jka,ir,mpar))
                  cla(jk,i,mpar) = p5*(cls(jk,i,mpar)-cls(jka,ir,mpar))
                  cls_0(jk)      = p5*(cls(jk,i,mpar)+cls(jka,ir,mpar))
               END DO

               crs(:,i,mpar) = crs_0(:)
               czs(:,i,mpar) = czs_0(:)
               cls(:,i,mpar) = cls_0(:)
            ENDIF

         END DO
      END DO

      DEALLOCATE (ars_0, brs_0, azs_0, bzs_0, bls_0,
     1            rcs_0, zcs_0, crs_0, czs_0, cls_0, stat=ir)


      END SUBROUTINE symforce

      SUBROUTINE symoutput (bsq     , gsqrt , bsubu , bsubv ,bsupu ,  
     1                      bsupv   , bsubs ,
#ifdef _ANIMEC 
     2                      ppar  , pperp ,densit,
     3                      sigma_an, tau_an, pbprim, ppprim,
#endif
     4                      bsqa    , gsqrta, bsubua, bsubva,bsupua,
     5                      bsupva  , bsubsa
#ifdef _ANIMEC
     6                     ,ppara , pperpa,densita,
     7                      sigma_ana,tau_ana,pbprima,ppprima
#endif
     6                      )

      USE vmec_main, p5 => cp5
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3), INTENT(inout) ::
     1   bsq, gsqrt, bsubu, bsubv, bsupu, bsupv, bsubs
#ifdef _ANIMEC
     2  ,ppar, pperp, sigma_an, tau_an, pbprim, ppprim, densit
#endif
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3), INTENT(out)   ::
     1   bsqa,gsqrta,bsubua,bsubva,bsupua,bsupva,bsubsa
#ifdef _ANIMEC
     2  ,ppara ,pperpa, sigma_ana, tau_ana, pbprima, ppprima ,
     3   densita
#endif
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ir, i, jk, jka
      REAL(rprec), DIMENSION(ns*nzeta) :: bsq_0, gsqrt_0, bsubu_0, 
     1    bsubv_0, bsupu_0, bsupv_0, bsubs_0
#ifdef _ANIMEC
     2  , ppar_0, pperp_0,
     2    sigma_an0 , tau_an0 , pbprim_0, ppprim_0, densit_0
#endif
C-----------------------------------------------

!
!       SYMMETRIZE FORCES ON RESTRICTED THETA INTERVAL (0 <= u <= pi)
!       SO COS,SIN INTEGRALS CAN BE PERFORMED. FOR EXAMPLE,
!
!       BSQ-S(v,u) = .5*( BSQ(v,u) + BSQ(-v,-u) )     ! * COS(mu - nv)
!       BSQ-A(v,u) = .5*( BSQ(v,u) - BSQ(-v,-u) )     ! * SIN(mu - nv)
!    
!       FOR BSUBS, THIS IS REVERSED, S-PIECE ~ SIN, A-PIECE ~ COS
!
!
      DO i = 1, ntheta2
         ir = ntheta1 + 2 - i                 !-theta
         IF (i == 1) ir = 1
         DO jk = 1, ns*nzeta
            jka = ireflect(jk)                !-zeta
            bsqa(jk,i)      = p5*(bsq(jk,i)     -bsq(jka,ir))
            bsq_0(jk)       = p5*(bsq(jk,i)     +bsq(jka,ir))
            gsqrta(jk,i)    = p5*(gsqrt(jk,i)   -gsqrt(jka,ir))
            gsqrt_0(jk)     = p5*(gsqrt(jk,i)   +gsqrt(jka,ir))
            bsubua(jk,i)    = p5*(bsubu(jk,i)   -bsubu(jka,ir))
            bsubu_0(jk)     = p5*(bsubu(jk,i)   +bsubu(jka,ir))
            bsubva(jk,i)    = p5*(bsubv(jk,i)   -bsubv(jka,ir))
            bsubv_0(jk)     = p5*(bsubv(jk,i)   +bsubv(jka,ir))
            bsupua(jk,i)    = p5*(bsupu(jk,i)   -bsupu(jka,ir))
            bsupu_0(jk)     = p5*(bsupu(jk,i)   +bsupu(jka,ir))
            bsupva(jk,i)    = p5*(bsupv(jk,i)   -bsupv(jka,ir))
            bsupv_0(jk)     = p5*(bsupv(jk,i)   +bsupv(jka,ir))
#ifdef _ANIMEC
            sigma_ana(jk,i) = p5*(sigma_an(jk,i)-sigma_an(jka,ir))
            sigma_an0(jk)   = p5*(sigma_an(jk,i)+sigma_an(jka,ir))
            tau_ana(jk,i)   = p5*(tau_an(jk,i)  -tau_an(jka,ir))
            tau_an0(jk)     = p5*(tau_an(jk,i)  +tau_an(jka,ir))
            ppara(jk,i)     = p5*(ppar(jk,i)    -ppar(jka,ir))
            ppar_0(jk)      = p5*(ppar(jk,i)    +ppar(jka,ir))	    	    
            pperpa(jk,i)    = p5*(pperp(jk,i)   -pperp(jka,ir))
            pperp_0(jk)     = p5*(pperp(jk,i)   +pperp(jka,ir))
            pbprima(jk,i)   = p5*(pbprim(jk,i)  -pbprim(jka,ir))
            pbprim_0(jk)    = p5*(pbprim(jk,i)  +pbprim(jka,ir))
            ppprima(jk,i)   = p5*(ppprim(jk,i)  -ppprim(jka,ir))
            ppprim_0(jk)    = p5*(ppprim(jk,i)  +ppprim(jka,ir))
            densita(jk,i)   = p5*(densit(jk,i)  -densit(jka,ir))
            densit_0(jk)    = p5*(densit(jk,i)  +densit(jka,ir))
#endif
! Dominant symmetry reversed
            bsubsa(jk,i)    = p5*(bsubs(jk,i)   +bsubs(jka,ir))  
            bsubs_0(jk)     = p5*(bsubs(jk,i)   -bsubs(jka,ir))
         END DO

         bsq(:,i)      = bsq_0(:)
         gsqrt(:,i)    = gsqrt_0(:)
         bsubu(:,i)    = bsubu_0(:)
         bsubv(:,i)    = bsubv_0(:)
         bsupu(:,i)    = bsupu_0(:)
         bsupv(:,i)    = bsupv_0(:)
         bsubs(:,i)    = bsubs_0(:)
#ifdef _ANIMEC
         sigma_an(:,i) = sigma_an0(:)
         tau_an(:,i)   = tau_an0(:)
         ppar(:,i)     = ppar_0(:)
         pperp(:,i)    = pperp_0(:)	 
         pbprim(:,i)   = pbprim_0(:)
         ppprim(:,i)   = ppprim_0(:)	 
         densit(:,i)   = densit_0(:)
#endif

      END DO

      END SUBROUTINE symoutput

!  Put the surface routines in a separate subroutine since these the quantites
!  these work on only exist on free boundary runs.
      SUBROUTINE symoutput_sur(bsubu, bsubv, bsupu, bsupv,                     &
     &                         bsubua, bsubva, bsupua, bsupva)
      USE vmec_main, p5 => cp5

      IMPLICIT NONE

C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(nzeta,ntheta3), INTENT(inout) ::
     1   bsubu, bsubv, bsupu, bsupv
      REAL(rprec), DIMENSION(nzeta,ntheta2), INTENT(out)   ::
     1   bsubua, bsubva, bsupua, bsupva

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ir, i, jk, jka
      REAL(rprec), DIMENSION(nzeta) :: bsubu_0, bsubv_0,
     1                                 bsupu_0, bsupv_0
C-----------------------------------------------

!
!       SYMMETRIZE FORCES ON RESTRICTED THETA INTERVAL (0 <= u <= pi)
!       SO COS,SIN INTEGRALS CAN BE PERFORMED. FOR EXAMPLE,
!
!       BSQ-S(v,u) = .5*( BSQ(v,u) + BSQ(-v,-u) )     ! * COS(mu - nv)
!       BSQ-A(v,u) = .5*( BSQ(v,u) - BSQ(-v,-u) )     ! * SIN(mu - nv)
!
!       FOR BSUBS, THIS IS REVERSED, S-PIECE ~ SIN, A-PIECE ~ COS
!
!

      ir = 1 !-theta
      DO i = 1, ntheta2
         jka = 1 !-zeta
         DO jk = 1, nzeta
            bsubua(jk,i) = p5*(bsubu(jk,i) - bsubu(jka,ir))
            bsubu_0(jk)  = p5*(bsubu(jk,i) + bsubu(jka,ir))
            bsubva(jk,i) = p5*(bsubv(jk,i) - bsubv(jka,ir))
            bsubv_0(jk)  = p5*(bsubv(jk,i) + bsubv(jka,ir))
            bsupua(jk,i) = p5*(bsupu(jk,i) - bsupu(jka,ir))
            bsupu_0(jk)  = p5*(bsupu(jk,i) + bsupu(jka,ir))
            bsupva(jk,i) = p5*(bsupv(jk,i) - bsupv(jka,ir))
            bsupv_0(jk)  = p5*(bsupv(jk,i) + bsupv(jka,ir))
            jka = nzeta - jk + 1

         END DO
         ir = ntheta3 - i + 1

         bsubu(:,i) = bsubu_0(:)
         bsubv(:,i) = bsubv_0(:)
         bsupu(:,i) = bsupu_0(:)
         bsupv(:,i) = bsupv_0(:)

      END DO

      END SUBROUTINE symoutput_sur
