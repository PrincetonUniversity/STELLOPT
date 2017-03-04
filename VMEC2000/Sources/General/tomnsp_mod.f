      MODULE tomnsp_mod

      CONTAINS

      SUBROUTINE tomnsps(frzl_array, armn, brmn, crmn, azmn, 
     1   bzmn, czmn, blmn, clmn, arcon, azcon)
      USE realspace, ONLY: wint, phip
      USE vmec_main, p5 => cp5
      USE vmec_params, ONLY: jlam, jmin2, ntmax, rcc, rss, zsc, zcs,
     1                       nscale
      USE fbal, ONLY: rru_fac, rzu_fac, frcc_fac, fzsc_fac
      USE precon2d, ONLY: ictrl_prec2d
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(out) :: frzl_array
      REAL(rprec), DIMENSION(ns*nzeta*ntheta3,0:1), INTENT(in) ::
     1   armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn, arcon, azcon
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: jmax, m, mparity, i, n, k, l, nsz
      INTEGER :: ioff, joff, mj, ni, nsl, j2, j2l, jl, jll, jmaxl 
      REAL(rprec), DIMENSION(:,:,:), POINTER :: 
     1           frcc, frss, fzcs, fzsc, flcs, flsc
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: work1
      REAL(rprec), DIMENSION(:), ALLOCATABLE   :: tempr, tempz
      REAL(rprec)  :: t1
!-----------------------------------------------
      frcc => frzl_array(:,:,:,rcc)               !!COS(mu) COS(nv)
      fzsc => frzl_array(:,:,:,zsc+ntmax)         !!SIN(mu) COS(nv)
      flsc => frzl_array(:,:,:,zsc+2*ntmax)       !!SIN(mu) COS(nv)
      IF (lthreed) THEN 
         frss => frzl_array(:,:,:,rss)               !!SIN(mu) SIN(nv)
         fzcs => frzl_array(:,:,:,zcs+ntmax)         !!COS(mu) SIN(nv)
         flcs => frzl_array(:,:,:,zcs+2*ntmax)       !!COS(mu) SIN(nv)
      END IF

      nsz = ns*nzeta

      ALLOCATE (work1(nsz,12), tempr(nsz), tempz(nsz),
     1          stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC2000 tomnsps'

      ioff = LBOUND(frcc,2)
      joff = LBOUND(frcc,3)

      frzl_array = 0

      jmax = ns
      IF (ivac .lt. 1) jmax = ns1

!
!     BEGIN FOURIER TRANSFORM
!
!       FRmn = ARmn - d(BRmn)/du + d(CRmn)/dv
!       FZmn = AZmn - d(BZmn)/du + d(CZmn)/dv
!       FLmn =      - d(BLmn)/du + d(CLmn)/dv
!
!       NOTE: sinmumi = -m sin(mu),  sinnvn = -n sin(nv)
!
      DO m = 0, mpol1
         mparity = MOD(m,2)
         mj = m+joff
         j2 = jmin2(m)
         jl = jlam(m)
         work1 = 0
!
!        DO THETA (U) INTEGRATION FIRST ON HALF INTERVAL (0 < U < PI)
!
         l = 0
         DO i = 1, ntheta2
            jll = l+1;  nsl = nsz+l
            l = l+nsz
            tempr(:) = armn(jll:nsl,mparity) 
#ifndef _HBANGLE
     1               + xmpq(m,1)*arcon(jll:nsl,mparity)
#endif
            tempz(:) = azmn(jll:nsl,mparity) 
#ifndef _HBANGLE
     1               + xmpq(m,1)*azcon(jll:nsl,mparity)
#endif
            work1(:,1) = work1(:,1) + tempr(:)*cosmui(i,m) 
     1                              + brmn(jll:nsl,mparity)*sinmumi(i,m)
            work1(:,7) = work1(:,7) + tempz(:)*sinmui(i,m)
     1                              + bzmn(jll:nsl,mparity)*cosmumi(i,m)
            work1(:,11)= work1(:,11)+ blmn(jll:nsl,mparity)*cosmumi(i,m)
 
            IF (.not.lthreed) CYCLE

            work1(:,2) = work1(:,2) - crmn(jll:nsl,mparity)*cosmui(i,m)
            work1(:,3) = work1(:,3) + tempr(:)*sinmui(i,m)
     1                              + brmn(jll:nsl,mparity)*cosmumi(i,m)
            work1(:,4) = work1(:,4) - crmn(jll:nsl,mparity)*sinmui(i,m)
            work1(:,5) = work1(:,5) + tempz(:)*cosmui(i,m)
     1                              + bzmn(jll:nsl,mparity)*sinmumi(i,m)
            work1(:,6) = work1(:,6) - czmn(jll:nsl,mparity)*cosmui(i,m)
            work1(:,8) = work1(:,8) - czmn(jll:nsl,mparity)*sinmui(i,m)

            work1(:,9) = work1(:,9) + blmn(jll:nsl,mparity)*sinmumi(i,m)
            work1(:,10) =work1(:,10)- clmn(jll:nsl,mparity)*cosmui(i,m)
            work1(:,12) =work1(:,12)- clmn(jll:nsl,mparity)*sinmui(i,m)
         END DO
!
!        NEXT, DO ZETA (V) TRANSFORM
!
         DO n = 0, ntor
            ni = n+ioff
            l = 0
            DO k = 1, nzeta
               j2l = j2+l; jmaxl = jmax+l; jll = jl+l; nsl = ns+l
               l = l+ns
               frcc(j2:jmax,ni,mj) = frcc(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,1)*cosnv(k,n)
               fzsc(j2:jmax,ni,mj) = fzsc(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,7)*cosnv(k,n)
               flsc(jl:ns,ni,mj) = flsc(jl:ns,ni,mj)
     1                           + work1(jll:nsl,11)*cosnv(k,n)

               IF (.not.lthreed) CYCLE

               frcc(j2:jmax,ni,mj) = frcc(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,2)*sinnvn(k,n)
               fzsc(j2:jmax,ni,mj) = fzsc(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,8)*sinnvn(k,n)
               frss(j2:jmax,ni,mj) = frss(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,3)*sinnv(k,n) 
     2                             + work1(j2l:jmaxl,4)*cosnvn(k,n)
               fzcs(j2:jmax,ni,mj) = fzcs(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,5)*sinnv(k,n)
     2                             + work1(j2l:jmaxl,6)*cosnvn(k,n)

               flsc(jl:ns,ni,mj) = flsc(jl:ns,ni,mj)
     1                           + work1(jll:nsl,12)*sinnvn(k,n)
               flcs(jl:ns,ni,mj) = flcs(jl:ns,ni,mj)
     1                           + work1(jll:nsl,9)*sinnv(k,n)
     2                           + work1(jll:nsl,10)*cosnvn(k,n)
            END DO
         END DO
      END DO
!
!     COMPUTE IOTA EVOLUTION EQUATION [STORED IN LMNSC(0,0) COMPONENT]
!
      IF (ictrl_prec2d.gt.0 .and. ncurr.eq.1) THEN
         ni = 0+ioff;  mj = 0+joff
         t1 = r0scale
         DO jl = 2, ns
            flsc(jl, ni, mj) = -t1*(buco(jl) - icurv(jl))
         END DO
      END IF

!
!     MAKE R,Z(m=1,n=0) SATISFY AVERAGE FORCE BALANCE EXACTLY
!     NOTE: for m=1, FR ~ Z1*(f0 + f2), FZ ~ R1*(f0 - f2), WHERE
!     f0 is the m=0 component of frho, f2 is m=2 component.
      IF (lforbal) THEN
         mj = 1 + joff
         ni = 0 + ioff
         t1 = nscale(0)*r0scale    !/4    !!v8.52

         DO jl = 2, ns-1
            work1(jl,1) = frcc_fac(jl)*frcc(jl,ni,mj)
     1                  + fzsc_fac(jl)*fzsc(jl,ni,mj)
            frcc(jl,ni,mj) = rzu_fac(jl)*(t1*equif(jl) + work1(jl,1))
            fzsc(jl,ni,mj) = rru_fac(jl)*(t1*equif(jl) - work1(jl,1))
         END DO
      END IF

      DEALLOCATE (work1, tempr, tempz)

      END SUBROUTINE tomnsps

#ifdef _TEST_FOURIER
!     TEST ROUTINE TO MAKE SURE FFTs ARE CONSISTENT
      SUBROUTINE tomnsps_t(rzl_array, rzl_array_in, r1, ru, rv, z1, zu, 
     1                     zv)
      USE realspace, ONLY: wint, phip
      USE vmec_main, p5 => cp5
      USE vmec_params, ONLY: jlam, jmin2, ntmax, rcc, rss, zsc, zcs
      USE precon2d, ONLY: ictrl_prec2d
      USE totzsp_mod, ONLY: convert_sym
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(out) :: rzl_array
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   INTENT(inout) :: rzl_array_in
      REAL(rprec), DIMENSION(ns*nzeta*ntheta3,0:1), INTENT(in) ::
     1   r1, ru, rv, z1, zu, zv
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: jmax, m, mparity, i, n, k, l, nsz
      INTEGER :: ioff, joff, mj, ni, nsl, j2, j2l, jl, jll, jmaxl 
      REAL(rprec), DIMENSION(:,:,:), POINTER :: 
     1           rmncc, rmnss, zmncs, zmnsc
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: work1
!-----------------------------------------------
      rmncc => rzl_array(:,:,:,rcc)               !!COS(mu) COS(nv)
      zmnsc => rzl_array(:,:,:,zsc+ntmax)         !!SIN(mu) COS(nv)
      IF (lthreed) THEN 
         rmnss => rzl_array(:,:,:,rss)               !!SIN(mu) SIN(nv)
         zmncs => rzl_array(:,:,:,zcs+ntmax)         !!COS(mu) SIN(nv)
      END IF

      nsz = ns*nzeta

      ALLOCATE (work1(nsz,12), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC2000 tomnsps_t'

      ioff = LBOUND(rmncc,2)
      joff = LBOUND(rmncc,3)

      rzl_array = 0

      jmax = ns
!
!     BEGIN FOURIER TRANSFORM
!
!     NOTE: sinmumi = -m sin(mu),  sinnvn = -n sin(nv)
!
      DO m = 0, mpol1
         mparity = MOD(m,2)
         mj = m+joff
         j2 = jmin2(m)
         jl = jlam(m)
         work1 = 0
!
!        DO THETA (U) INTEGRATION FIRST ON HALF INTERVAL (0 < U < PI)
!
         l = 0
         DO i = 1, ntheta2
            jll = l+1;  nsl = nsz+l
            l = l+nsz
            work1(:,1) = work1(:,1) + r1(jll:nsl,mparity)*cosmui(i,m)
            work1(:,7) = work1(:,7) + z1(jll:nsl,mparity)*sinmui(i,m)
 
            IF (.not.lthreed) CYCLE

!            work1(:,2) = work1(:,2) - crmn(jll:nsl,mparity)*cosmui(i,m)
            work1(:,3) = work1(:,3) + r1(jll:nsl,mparity)*sinmui(i,m)
!            work1(:,4) = work1(:,4) - crmn(jll:nsl,mparity)*sinmui(i,m)
            work1(:,5) = work1(:,5) + z1(jll:nsl,mparity)*cosmui(i,m)
!            work1(:,6) = work1(:,6) - czmn(jll:nsl,mparity)*cosmui(i,m)
!           work1(:,8) = work1(:,8) - czmn(jll:nsl,mparity)*sinmui(i,m)

!            work1(:,10) =work1(:,10)- clmn(jll:nsl,mparity)*cosmui(i,m)
!            work1(:,12) =work1(:,12)- clmn(jll:nsl,mparity)*sinmui(i,m)
         END DO
!
!        NEXT, DO ZETA (V) TRANSFORM
!
         DO n = 0, ntor
            ni = n+ioff
            l = 0
            DO k = 1, nzeta
               j2l = j2+l; jmaxl = jmax+l; jll = jl+l; nsl = ns+l
               l = l+ns
               rmncc(j2:jmax,ni,mj) = rmncc(j2:jmax,ni,mj)
     1                              + work1(j2l:jmaxl,1)*cosnv(k,n)
               zmnsc(j2:jmax,ni,mj) = zmnsc(j2:jmax,ni,mj)
     1                              + work1(j2l:jmaxl,7)*cosnv(k,n)

               IF (.not.lthreed) CYCLE

!               frcc(j2:jmax,ni,mj) = frcc(j2:jmax,ni,mj)
!     1                             + work1(j2l:jmaxl,2)*sinnvn(k,n)
!               fzsc(j2:jmax,ni,mj) = fzsc(j2:jmax,ni,mj)
!     1                             + work1(j2l:jmaxl,8)*sinnvn(k,n)
               rmnss(j2:jmax,ni,mj) = rmnss(j2:jmax,ni,mj)
     1                              + work1(j2l:jmaxl,3)*sinnv(k,n) 
!     2                             + work1(j2l:jmaxl,4)*cosnvn(k,n)
               zmncs(j2:jmax,ni,mj) = zmncs(j2:jmax,ni,mj)
     1                              + work1(j2l:jmaxl,5)*sinnv(k,n)
!     2                             + work1(j2l:jmaxl,6)*cosnvn(k,n)

!               flsc(jl:ns,ni,mj) = flsc(jl:ns,ni,mj)
!     1                           + work1(jll:nsl,12)*sinnvn(k,n)
            END DO
         END DO
      END DO

      DEALLOCATE (work1)

      IF (lthreed) CALL convert_sym(rzl_array_in(:,:,:,rss), 
     1                              rzl_array_in(:,:,:,ntmax+zcs))

      IF (lasym) RETURN
         
      IF (ANY(ABS(rzl_array(:,:,:,1:2*ntmax) - 
     1            rzl_array_in(:,:,:,1:2*ntmax)) .GT. 1.E-12_dp)) THEN
         STOP 'FFT ERROR!'
      END IF

      END SUBROUTINE tomnsps_t
#endif

      SUBROUTINE tomnspa(frzl_array, armn, brmn, crmn, azmn, bzmn,
     1   czmn, blmn, clmn, arcon, azcon)
      USE vmec_main
      USE vmec_params, ONLY: jlam, jmin2, ntmax, rsc, rcs, zcc, zss
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(inout) :: frzl_array
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(in) ::
     1   armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn, arcon, azcon
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: jmax, m, mparity, i, n, k, l, nsz
      INTEGER :: ioff, joff, mj, ni, nsl, j2, j2l, jl, jll, jmaxl 
      REAL(rprec), DIMENSION(:,:,:), POINTER :: 
     1           frcs, frsc, fzcc, fzss, flcc, flss
      REAL(rprec), DIMENSION(:), ALLOCATABLE   :: temp1, temp3
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work1
!-----------------------------------------------
      frsc => frzl_array(:,:,:,rsc)               !!R-SIN(mu) COS(nv)
      fzcc => frzl_array(:,:,:,zcc+ntmax)         !!Z-COS(mu) COS(nv)
      flcc => frzl_array(:,:,:,zcc+2*ntmax)       !!L-COS(mu) COS(nv)
      IF (lthreed) THEN
         frcs => frzl_array(:,:,:,rcs)            !!R-COS(mu) SIN(nv)
         fzss => frzl_array(:,:,:,zss+ntmax)      !!Z-SIN(mu) SIN(nv)
         flss => frzl_array(:,:,:,zss+2*ntmax)    !!L-SIN(mu) SIN(nv)
      END IF

      nsz = ns*nzeta
      ALLOCATE (work1(nsz,12), temp1(nsz), temp3(nsz), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC tomnspa'

      ioff = LBOUND(frsc,2)
      joff = LBOUND(frsc,3)

      jmax = ns
      IF (ivac .lt. 1) jmax = ns1
!
!     BEGIN FOURIER TRANSFORM
!
      DO m = 0, mpol1
         mparity = MOD(m,2)
         mj = m+joff
         j2 = jmin2(m)
         jl = jlam(m)
         work1 = 0
!
!        DO THETA (U) TRANSFORM FIRST
!
         DO i = 1, ntheta2
            temp1(:) = armn(:,i,mparity)
#ifndef _HBANGLE
     1                + xmpq(m,1)*arcon(:,i,mparity)
#endif
            temp3(:) = azmn(:,i,mparity)
#ifndef _HBANGLE
     1               + xmpq(m,1)*azcon(:,i,mparity)
#endif
            work1(:,3) = work1(:,3) + temp1(:)*sinmui(i,m) 
     1                              + brmn(:,i,mparity)*cosmumi(i,m)
            work1(:,5) = work1(:,5) + temp3(:)*cosmui(i,m)
     1                              + bzmn(:,i,mparity)*sinmumi(i,m)
            work1(:,9) = work1(:,9) + blmn(:,i,mparity)*sinmumi(i,m)

            IF (.not.lthreed) CYCLE

            work1(:,1) = work1(:,1) + temp1(:)*cosmui(i,m)
     1                              + brmn(:,i,mparity)*sinmumi(i,m)
            work1(:,2) = work1(:,2) - crmn(:,i,mparity)*cosmui(i,m)
            work1(:,4) = work1(:,4) - crmn(:,i,mparity)*sinmui(i,m)
            work1(:,6) = work1(:,6) - czmn(:,i,mparity)*cosmui(i,m)
            work1(:,7) = work1(:,7) + temp3(:)*sinmui(i,m)
     1                              + bzmn(:,i,mparity)*cosmumi(i,m)
            work1(:,8) = work1(:,8) - czmn(:,i,mparity)*sinmui(i,m)
            work1(:,10) = work1(:,10) - clmn(:,i,mparity)*cosmui(i,m)
            work1(:,11) = work1(:,11) + blmn(:,i,mparity)*cosmumi(i,m)
            work1(:,12) = work1(:,12) - clmn(:,i,mparity)*sinmui(i,m)
         END DO
!
!        NEXT, DO ZETA (V) TRANSFORM
!
         DO n = 0, ntor
            ni = n+ioff
            DO k = 1, nzeta
               l = ns*(k-1)
               j2l = j2+l; jmaxl = jmax+l; jll = jl+l; nsl = ns+l
               frsc(j2:jmax,ni,mj) = frsc(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,3)*cosnv(k,n)
               fzcc(j2:jmax,ni,mj) = fzcc(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,5)*cosnv(k,n)
               flcc(jl:ns,ni,mj) = flcc(jl:ns,ni,mj)
     1                           + work1(jll:nsl,9)*cosnv(k,n)

               IF (.not.lthreed) CYCLE

               frsc(j2:jmax,ni,mj) = frsc(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,4)*sinnvn(k,n)
               fzcc(j2:jmax,ni,mj) = fzcc(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,6)*sinnvn(k,n)
               frcs(j2:jmax,ni,mj) = frcs(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,1)*sinnv(k,n) 
     2                             + work1(j2l:jmaxl,2)*cosnvn(k,n)
               fzss(j2:jmax,ni,mj) = fzss(j2:jmax,ni,mj)
     1                             + work1(j2l:jmaxl,7)*sinnv(k,n)
     2                             + work1(j2l:jmaxl,8)*cosnvn(k,n)
               flcc(jl:ns,ni,mj) = flcc(jl:ns,ni,mj)
     1                           + work1(jll:nsl,10)*sinnvn(k,n)
               flss(jl:ns,ni,mj) = flss(jl:ns,ni,mj)
     1                           + work1(jll:nsl,11)*sinnv(k,n)
     2                           + work1(jll:nsl,12)*cosnvn(k,n)
            END DO
         END DO
      END DO

      DEALLOCATE (work1, temp1, temp3)

      END SUBROUTINE tomnspa

#ifdef _TEST_FOURIER
      SUBROUTINE tomnspa_t(rzl_array, rzl_array_in, r1, ru, rv, z1, zu, 
     1                     zv)
      USE realspace, ONLY: wint, phip
      USE vmec_main, p5 => cp5
      USE vmec_params, ONLY: jlam, jmin2, ntmax, rsc, rcs, zcc, zss
      USE precon2d, ONLY: ictrl_prec2d
      USE totzsp_mod, ONLY: convert_asym
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(out) :: rzl_array
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   INTENT(inout) :: rzl_array_in
      REAL(rprec), DIMENSION(ns*nzeta*ntheta3,0:1), INTENT(in) ::
     1   r1, ru, rv, z1, zu, zv
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: jmax, m, mparity, i, n, k, l, nsz
      INTEGER :: ioff, joff, mj, ni, nsl, j2, j2l, jl, jll, jmaxl 
      REAL(rprec), DIMENSION(:,:,:), POINTER :: 
     1           rmnsc, rmncs, zmncc, zmnss
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: work1
!-----------------------------------------------
      rmnsc => rzl_array(:,:,:,rsc)               !!COS(mu) COS(nv)
      zmncc => rzl_array(:,:,:,zcc+ntmax)         !!SIN(mu) COS(nv)
      IF (lthreed) THEN 
         rmncs => rzl_array(:,:,:,rcs)               !!SIN(mu) SIN(nv)
         zmnss => rzl_array(:,:,:,zss+ntmax)         !!COS(mu) SIN(nv)
      END IF

      nsz = ns*nzeta

      ALLOCATE (work1(nsz,12), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC2000 tomnsps_t'

      ioff = LBOUND(rmnsc,2)
      joff = LBOUND(rmnsc,3)

      jmax = ns
!
!     BEGIN FOURIER TRANSFORM
!
!     NOTE: sinmumi = -m sin(mu),  sinnvn = -n sin(nv)
!
      DO m = 0, mpol1
         mparity = MOD(m,2)
         mj = m+joff
         j2 = jmin2(m)
         jl = jlam(m)
         work1 = 0
!
!        DO THETA (U) INTEGRATION FIRST ON HALF INTERVAL (0 < U < PI)
!
         l = 0
         DO i = 1, ntheta2
            jll = l+1;  nsl = nsz+l
            l = l+nsz
            work1(:,1) = work1(:,1) + r1(jll:nsl,mparity)*sinmui(i,m)
            work1(:,7) = work1(:,7) + z1(jll:nsl,mparity)*cosmui(i,m)
 
            IF (.not.lthreed) CYCLE

!            work1(:,2) = work1(:,2) - crmn(jll:nsl,mparity)*cosmui(i,m)
            work1(:,3) = work1(:,3) + r1(jll:nsl,mparity)*cosmui(i,m)
!            work1(:,4) = work1(:,4) - crmn(jll:nsl,mparity)*sinmui(i,m)
            work1(:,5) = work1(:,5) + z1(jll:nsl,mparity)*sinmui(i,m)
!            work1(:,6) = work1(:,6) - czmn(jll:nsl,mparity)*cosmui(i,m)
!           work1(:,8) = work1(:,8) - czmn(jll:nsl,mparity)*sinmui(i,m)

!            work1(:,10) =work1(:,10)- clmn(jll:nsl,mparity)*cosmui(i,m)
!            work1(:,12) =work1(:,12)- clmn(jll:nsl,mparity)*sinmui(i,m)
         END DO
!
!        NEXT, DO ZETA (V) TRANSFORM
!
         DO n = 0, ntor
            ni = n+ioff
            l = 0
            DO k = 1, nzeta
               j2l = j2+l; jmaxl = jmax+l; jll = jl+l; nsl = ns+l
               l = l+ns
               rmnsc(j2:jmax,ni,mj) = rmnsc(j2:jmax,ni,mj)
     1                              + work1(j2l:jmaxl,1)*cosnv(k,n)
               zmncc(j2:jmax,ni,mj) = zmncc(j2:jmax,ni,mj)
     1                              + work1(j2l:jmaxl,7)*cosnv(k,n)

               IF (.not.lthreed) CYCLE

!               frcc(j2:jmax,ni,mj) = frcc(j2:jmax,ni,mj)
!     1                             + work1(j2l:jmaxl,2)*sinnvn(k,n)
!               fzsc(j2:jmax,ni,mj) = fzsc(j2:jmax,ni,mj)
!     1                             + work1(j2l:jmaxl,8)*sinnvn(k,n)
               rmncs(j2:jmax,ni,mj) = rmncs(j2:jmax,ni,mj)
     1                              + work1(j2l:jmaxl,3)*sinnv(k,n) 
!     2                             + work1(j2l:jmaxl,4)*cosnvn(k,n)
               zmnss(j2:jmax,ni,mj) = zmnss(j2:jmax,ni,mj)
     1                              + work1(j2l:jmaxl,5)*sinnv(k,n)
!     2                             + work1(j2l:jmaxl,6)*cosnvn(k,n)

!               flsc(jl:ns,ni,mj) = flsc(jl:ns,ni,mj)
!     1                           + work1(jll:nsl,12)*sinnvn(k,n)
            END DO
         END DO
      END DO

!     IF THE SYMMETRIZED MODE USED, NEED EXTRA FACTOR OF 2 
!     IF ntheta3 USED INSTEAD OF ntheta3, DO NOT NEED THIS FACTOR
      rzl_array = 2*rzl_array

      CALL convert_asym(rzl_array_in(:,:,:,rsc), 
     1                  rzl_array_in(:,:,:,zcc+ntmax))

      IF (ANY(ABS(rzl_array(:,:,:,1:2*ntmax) - 
     1            rzl_array_in(:,:,:,1:2*ntmax)) .GT. 1.E-12_dp)) THEN
         STOP 'FFT ERROR 1!'
      END IF
            
      DEALLOCATE (work1)

      END SUBROUTINE tomnspa_t
#endif
      END MODULE tomnsp_mod
