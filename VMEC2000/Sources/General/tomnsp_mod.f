      MODULE tomnsp_mod
      USE timer_sub
      IMPLICIT NONE

      CONTAINS

      SUBROUTINE tomnsps_par(frzl_array, armn, brmn, crmn, azmn, 
     &                       bzmn, czmn, blmn, clmn, arcon, azcon)
      USE realspace, ONLY: wint, phip
      USE vmec_main, p5 => cp5
      USE vmec_params, ONLY: jlam, jmin2, ntmax, rcc, rss, zsc, zcs,
     &                       nscale
      USE fbal, ONLY: rru_fac, rzu_fac, frcc_fac, fzsc_fac
      USE precon2d, ONLY: ictrl_prec2d
      USE parallel_include_module
      USE xstuff
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,0:mpol1,ns,3*ntmax),
     &   TARGET, INTENT(out) :: frzl_array
      REAL(dp), DIMENSION(nzeta,ntheta3,ns,0:1), INTENT(INout) ::
     &   armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn, arcon, azcon
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: m0 = 0, m1 = 1, n0 = 0
      INTEGER :: jmax, m, mparity, i, n, k, l, nsz
      INTEGER :: ioff, joff, mj, ni, nsl, j2, j2l, jl, jll, jmaxl 
      REAL(dp), DIMENSION(:,:,:), POINTER :: 
     &           frcc, frss, fzcs, fzsc, flcs, flsc
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: work1
      REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: tempr, tempz
      REAL(dp)  :: t1
      INTEGER :: j, nsmin, nsmax, ub1, lb1, ub2, lb2, js
!-----------------------------------------------
      CALL second0 (tffton)

      frcc => frzl_array(:,:,:,rcc)               !!COS(mu) COS(nv)
      fzsc => frzl_array(:,:,:,zsc+ntmax)         !!SIN(mu) COS(nv)
      flsc => frzl_array(:,:,:,zsc+2*ntmax)       !!SIN(mu) COS(nv)
      IF (lthreed) THEN 
         frss => frzl_array(:,:,:,rss)               !!SIN(mu) SIN(nv)
         fzcs => frzl_array(:,:,:,zcs+ntmax)         !!COS(mu) SIN(nv)
         flcs => frzl_array(:,:,:,zcs+2*ntmax)       !!COS(mu) SIN(nv)
      END IF

      nsz = ns*nzeta

      nsmin = tlglob
      nsmax = trglob

      ALLOCATE (work1(12,nzeta,nsmin:nsmax), stat=i)
      ALLOCATE (tempr(nzeta,nsmin:nsmax), stat=i)
      ALLOCATE (tempz(nzeta,nsmin:nsmax), stat=i)
      IF (i .ne. 0) THEN
         STOP 'Allocation error in VMEC2000 tomnsps'
      END IF

      ioff = LBOUND(frcc,1)
      joff = LBOUND(frcc,2)

      jmax = ns
      IF (ivac .LT. 1) THEN
         jmax = ns1
      END IF

!
!     BEGIN FOURIER TRANSFORM
!
!     FRmn = ARmn - d(BRmn)/du + d(CRmn)/dv
!     FZmn = AZmn - d(BZmn)/du + d(CZmn)/dv
!     FLmn =      - d(BLmn)/du + d(CLmn)/dv
!
!     NOTE: sinmumi = -m sin(mu),  sinnvn = -n sin(nv)
!
      DO js = nsmin, nsmax
         frzl_array(:,:,js,:) = 0
         DO m = 0, mpol1
            mparity = MOD(m,2)
            work1(:,:,js) = 0

!        DO THETA (U) INTEGRATION FIRST ON HALF INTERVAL (0 < U < PI)
            DO i = 1, ntheta2
               DO k = 1, nzeta
                  tempr(k,js) = armn(k,i,js,mparity)
#ifndef _HBANGLE
     &                        + xmpq(m,1)*arcon(k,i,js,mparity)
#endif
                  tempz(k,js) = azmn(k,i,js,mparity)
#ifndef _HBANGLE
     &                        + xmpq(m,1)*azcon(k,i,js,mparity)
#endif
                  work1(1,k,js) = work1(1,k,js)
     &                          + tempr(k,js)*cosmui(i,m)
     &                          + brmn(k,i,js,mparity)*sinmumi(i,m)
                  work1(7,k,js) = work1(7,k,js)
     &                          + tempz(k,js)*sinmui(i,m)
     &                          + bzmn(k,i,js,mparity)*cosmumi(i,m)
                  work1(11,k,js) = work1(11,k,js)
     &                           + blmn(k,i,js,mparity)*cosmumi(i,m)

                  IF (.NOT.lthreed) CYCLE

                  work1(2,k,js) = work1(2,k,js)
     &                          - crmn(k,i,js,mparity)*cosmui(i,m)
                  work1(3,k,js) = work1(3,k,js)
     &                          + tempr(k,js)*sinmui(i,m)
     &                          + brmn(k,i,js,mparity)*cosmumi(i,m)
                  work1(4,k,js) = work1(4,k,js)
     &                          - crmn(k,i,js,mparity)*sinmui(i,m)
                  work1(5,k,js) = work1(5,k,js)
     &                          + tempz(k,js)*cosmui(i,m)
     &                          + bzmn(k,i,js,mparity)*sinmumi(i,m)
                  work1(6,k,js) = work1(6,k,js)
     &                          - czmn(k,i,js,mparity)*cosmui(i,m)
                  work1(8,k,js) = work1(8,k,js)
     &                          - czmn(k,i,js,mparity)*sinmui(i,m)
                  work1(9,k,js) = work1(9,k,js)
     &                          + blmn(k,i,js,mparity)*sinmumi(i,m)
                  work1(10,k,js) = work1(10,k,js)
     &                           - clmn(k,i,js,mparity)*cosmui(i,m)
                  work1(12,k,js) = work1(12,k,js)
     &                           - clmn(k,i,js,mparity)*sinmui(i,m)
               END DO
            END DO

!
!        NEXT, DO ZETA (V) TRANSFORM
            mj = m + joff
            j2 = jmin2(m)
            jl = jlam(m)

            lb1 = MAX(tlglob,j2)
            ub1 = MIN(trglob,jmax)
            lb2 = MAX(tlglob,jl)
            ub2 = trglob


            DO n = 0, ntor
               ni = n+ioff
               DO k = 1, nzeta

                  IF (lb1 .LE. js .AND. js .LE. ub1) THEN
                     frcc(ni,mj,js) = frcc(ni,mj,js)
     &                              + work1(1,k,js)*cosnv(k,n)

                     fzsc(ni,mj,js) = fzsc(ni,mj,js)
     &                              + work1(7,k,js)*cosnv(k,n)
                  END IF

                  IF (lb2 .LE. js .AND. js .LE. ub2) THEN
                     flsc(ni,mj,js) = flsc(ni,mj,js)
     &                              + work1(11,k,js)*cosnv(k,n)
                  END IF

                  IF (.NOT.lthreed) CYCLE

                  IF (lb1 .LE. js .AND. js .LE. ub1) THEN
                     frcc(ni,mj,js) = frcc(ni,mj,js)
     &                              + work1(2,k,js)*sinnvn(k,n)

                     fzsc(ni,mj,js) = fzsc(ni,mj,js)
     &                              + work1(8,k,js)*sinnvn(k,n)

                     frss(ni,mj,js) = frss(ni,mj,js)
     &                              + work1(3,k,js)*sinnv(k,n)
     &                              + work1(4,k,js)*cosnvn(k,n)

                     fzcs(ni,mj,js) = fzcs(ni,mj,js)
     &                              + work1(5,k,js)*sinnv(k,n)
     &                              + work1(6,k,js)*cosnvn(k,n)
                  END IF

                  IF (lb2 .LE. js .AND. js .LE. ub2) THEN
                     flsc(ni,mj,js) = flsc(ni,mj,js)
     &                              + work1(12,k,js)*sinnvn(k,n)

                     flcs(ni,mj,js) = flcs(ni,mj,js)
     &                              + work1(9,k,js)*sinnv(k,n)
     &                              + work1(10,k,js)*cosnvn(k,n)
                  END IF
               END DO
            END DO
         END DO
      END DO

!
!     COMPUTE IOTA EVOLUTION EQUATION [STORED IN LMNSC(0,0) COMPONENT]
!
!SPH071017
#if defined(CHI_FORCE)
      IF (ictrl_prec2d .NE. 0 .AND. ncurr .EQ. 1) THEN
         ni = n0 + ioff
         mj = m0 + joff
         t1 = r0scale
         nsmin = MAX(2,tlglob)
         nsmax = trglob
         DO js = nsmin, nsmax
            flsc(ni, mj, js) = -t1*(buco(js) - icurv(js))
         END DO
      END IF
#endif
!
!     MAKE R,Z(m=1,n=0) SATISFY AVERAGE FORCE BALANCE EXACTLY
!     NOTE: for m=1, FR ~ Z1*(f0 + f2), FZ ~ R1*(f0 - f2), WHERE
!     f0 is the m=0 component of frho, f2 is m=2 component.
      IF (lforbal) THEN
         ni = m0 + ioff
         mj = m1 + joff
         t1 = nscale(n0)*r0scale !/4    !!v8.52
         nsmin = MAX(2,tlglob)
         nsmax = MIN(trglob,ns-1)
         DO jl = nsmin, nsmax
            DO k = 1, nzeta
               work1(k,1,jl) = frcc_fac(jl)*frcc(ni,mj,jl)
     &                       + fzsc_fac(jl)*fzsc(ni,mj,jl)
               frcc(ni,mj,jl) = rzu_fac(jl)*(t1*equif(jl)
     &                        + work1(k,1,jl))
               fzsc(ni,mj,jl) = rru_fac(jl)*(t1*equif(jl)
     &                        - work1(k,1,jl))
            END DO
         END DO
      END IF

      DEALLOCATE (work1, tempr, tempz)

      CALL second0 (tfftoff)
      tomnsps_time = tomnsps_time + (tfftoff - tffton)
      timer(tffi) = timer(tffi) + (tfftoff - tffton)

      END SUBROUTINE tomnsps_par

      SUBROUTINE tomnspa_par(frzl_array, armn, brmn, crmn, azmn, bzmn,
     &                       czmn, blmn, clmn, arcon, azcon)
      USE vmec_main
      USE vmec_params, ONLY: jlam, jmin2, ntmax, rsc, rcs, zcc, zss
      USE parallel_include_module
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,0:mpol1,ns,3*ntmax),
     &   TARGET, INTENT(inout) :: frzl_array
      REAL(dp), DIMENSION(nzeta,ntheta3,ns,0:1), INTENT(in) ::
     &   armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn, arcon, azcon
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: jmax, m, mparity, i, n, k, l
      INTEGER :: ioff, joff, mj, ni, nsl, j2, j2l, jl, jll, jmaxl 
      REAL(dp), DIMENSION(:,:,:), POINTER :: 
     &           frcs, frsc, fzcc, fzss, flcc, flss
!      REAL(dp), DIMENSION(ns*nzeta) :: temp1, temp3
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: temp1, temp3
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: work1
      INTEGER :: j, nsmin, nsmax, ub1, lb1, ub2, lb2, js
!-----------------------------------------------
      CALL second0(tffton)

      frsc => frzl_array(:,:,:,rsc)               !!R-SIN(mu) COS(nv)
      fzcc => frzl_array(:,:,:,zcc+ntmax)         !!Z-COS(mu) COS(nv)
      flcc => frzl_array(:,:,:,zcc+2*ntmax)       !!L-COS(mu) COS(nv)
      IF (lthreed) THEN
         frcs => frzl_array(:,:,:,rcs)            !!R-COS(mu) SIN(nv)
         fzss => frzl_array(:,:,:,zss+ntmax)      !!Z-SIN(mu) SIN(nv)
         flss => frzl_array(:,:,:,zss+2*ntmax)    !!L-SIN(mu) SIN(nv)
      END IF

      nsmin = tlglob
      nsmax = trglob
      ALLOCATE (work1(12,nzeta,nsmin:nsmax),
     &          temp1(nzeta,nsmin:nsmax),
     &          temp3(nzeta,nsmin:nsmax), stat=i)
      IF (i .NE. 0) THEN
         STOP 'Allocation error in VMEC tomnspa'
      END IF

      ioff = LBOUND(frsc,1)
      joff = LBOUND(frsc,2)

      jmax = ns
      IF (ivac .LT. 1) THEN
         jmax = ns1
      END IF

!
!     BEGIN FOURIER TRANSFORM
!
      DO js = nsmin, nsmax
         DO m = 0, mpol1
            mparity = MOD(m,2)
            mj = m + joff
            j2 = jmin2(m)
            jl = jlam(m)
            work1(:,:,js) = 0
!
!        DO THETA (U) TRANSFORM FIRST
!
            DO i = 1, ntheta2
               DO k = 1, nzeta
                  temp1(k,js) = armn(k,i,js,mparity)
#ifndef _HBANGLE
     &                        + xmpq(m,1)*arcon(k,i,js,mparity)
#endif
                  temp3(k,js) = azmn(k,i,js,mparity)
#ifndef _HBANGLE
     &                        + xmpq(m,1)*azcon(k,i,js,mparity)
#endif
                  work1(3,k,js) = work1(3,k,js)
     &                          + temp1(k,js)*sinmui(i,m)
     &                          + brmn(k,i,js,mparity)*cosmumi(i,m)
                  work1(5,k,js) = work1(5,k,js)
     &                          + temp3(k,js)*cosmui(i,m)
     &                          + bzmn(k,i,js,mparity)*sinmumi(i,m)
                  work1(9,k,js) = work1(9,k,js)
     &                          + blmn(k,i,js,mparity)*sinmumi(i,m)

                  IF (.not.lthreed) CYCLE

                  work1(1,k,js) = work1(1,k,js)
     &                          + temp1(k,js)*cosmui(i,m)
     &                          + brmn(k,i,js,mparity)*sinmumi(i,m)
                  work1(2,k,js) = work1(2,k,js)
     &                          - crmn(k,i,js,mparity)*cosmui(i,m)
                  work1(4,k,js) = work1(4,k,js)
     &                          - crmn(k,i,js,mparity)*sinmui(i,m)
                  work1(6,k,js) = work1(6,k,js)
     &                          - czmn(k,i,js,mparity)*cosmui(i,m)
                  work1(7,k,js) = work1(7,k,js)
     &                          + temp3(k,js)*sinmui(i,m)
     &                          + bzmn(k,i,js,mparity)*cosmumi(i,m)
                  work1(8,k,js) = work1(8,k,js)
     &                          - czmn(k,i,js,mparity)*sinmui(i,m)
                  work1(10,k,js) = work1(10,k,js)
     &                           - clmn(k,i,js,mparity)*cosmui(i,m)
                  work1(11,k,js) = work1(11,k,js)
     &                           + blmn(k,i,js,mparity)*cosmumi(i,m)
                  work1(12,k,js) = work1(12,k,js)
     &                           - clmn(k,i,js,mparity)*sinmui(i,m)
               END DO
            END DO
!
!        NEXT, DO ZETA (V) TRANSFORM
!

            lb1 = MAX(tlglob,j2)
            ub1 = MIN(trglob,jmax)
            lb2 = MAX(tlglob,jl)
            ub2 = trglob

            DO n = 0, ntor
               ni = n + ioff
               DO k = 1, nzeta

                  IF (lb1 .LE. js .AND. js .LE. ub1) THEN
                     frsc(ni,mj,js) = frsc(ni,mj,js)
     &                              + work1(3,k,js)*cosnv(k,n)
                     fzcc(ni,mj,js) = fzcc(ni,mj,js)
     &                              + work1(5,k,js)*cosnv(k,n)
                  END IF

                  IF (lb2 .LE. js .AND. js .LE. ub2) THEN
                     flcc(ni,mj,js) = flcc(ni,mj,js)
     &                              + work1(9,k,js)*cosnv(k,n)
                  END IF

                  IF (.not.lthreed) CYCLE

                  IF (lb1 .LE. js .AND. js .LE. ub1) THEN
                     frsc(ni,mj,js) = frsc(ni,mj,js)
     &                              + work1(4,k,js)*sinnvn(k,n)
                     fzcc(ni,mj,js) = fzcc(ni,mj,js)
     &                              + work1(6,k,js)*sinnvn(k,n)
                     frcs(ni,mj,js) = frcs(ni,mj,js)
     &                              + work1(1,k,js)*sinnv(k,n)
     &                              + work1(2,k,js)*cosnvn(k,n)
                     fzss(ni,mj,js) = fzss(ni,mj,js)
     &                              + work1(7,k,js)*sinnv(k,n)
     &                              + work1(8,k,js)*cosnvn(k,n)
                  END IF

                  IF (lb2 .LE. js .AND. js .LE. ub2) THEN
                     flcc(ni,mj,js) = flcc(ni,mj,js)
     &                              + work1(10,k,js)*sinnvn(k,n)
                     flss(ni,mj,js) = flss(ni,mj,js)
     &                              + work1(11,k,js)*sinnv(k,n)
     &                              + work1(12,k,js)*cosnvn(k,n)
                  END IF
               END DO
            END DO
         END DO
      END DO

!     IF THE SYMMETRIZED MODE USED, NEED EXTRA FACTOR OF 2 
!     IF ntheta3 USED INSTEAD OF ntheta3, DO NOT NEED THIS FACTOR
!      frzl_array(:,:,nsmin:nsmax,:) = 2*frzl_array(:,:,nsmin:nsmax,:)

      DEALLOCATE (work1, temp1, temp3)
      CALL second0(tfftoff)
      tomnspa_time = tomnspa_time + (tfftoff - tffton)
      timer(tffi) = timer(tffi) + (tfftoff - tffton)

      END SUBROUTINE tomnspa_par

      SUBROUTINE tomnsps(frzl_array, armn, brmn, crmn, azmn, 
     1   bzmn, czmn, blmn, clmn, arcon, azcon)
      USE realspace, ONLY: wint, phip
      USE vmec_main, p5 => cp5
      USE vmec_params, ONLY: jlam, jmin2, ntmax, rcc, rss, zsc, zcs,
     1                       nscale
      USE fbal, ONLY: rru_fac, rzu_fac, frcc_fac, fzsc_fac
      USE precon2d, ONLY: ictrl_prec2d
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(out) :: frzl_array
      REAL(dp), DIMENSION(ns*nzeta*ntheta3,0:1), INTENT(in) ::
     1   armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn, arcon, azcon
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: jmax, m, mparity, i, n, k, l, nsz
      INTEGER :: ioff, joff, mj, ni, nsl, j2, j2l, jl, jll, jmaxl 
      REAL(dp), DIMENSION(:,:,:), POINTER :: 
     1           frcc, frss, fzcs, fzsc, flcs, flsc
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: work1
      REAL(dp), DIMENSION(:), ALLOCATABLE   :: tempr, tempz
      REAL(dp)  :: t1
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
      IF (i .ne. 0) THEN
         STOP 'Allocation error in VMEC2000 tomnsps'
      END IF

      ioff = LBOUND(frcc,2)
      joff = LBOUND(frcc,3)

      frzl_array = 0

      jmax = ns
      IF (ivac .lt. 1) THEN
         jmax = ns1
      END IF

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
         work1 = 0
!        DO THETA (U) INTEGRATION FIRST ON HALF INTERVAL (0 < U < PI)
!
         l = 0
         DO i = 1, ntheta2
            jll = l + 1
            nsl = nsz + l
            l = l + nsz
            tempr(:) = armn(jll:nsl,mparity) 
#ifndef _HBANGLE
     &               + xmpq(m,1)*arcon(jll:nsl,mparity)
#endif
            tempz(:) = azmn(jll:nsl,mparity) 
#ifndef _HBANGLE
     &               + xmpq(m,1)*azcon(jll:nsl,mparity)
#endif
            work1(:,1) = work1(:,1) + tempr(:)*cosmui(i,m) 
     &                              + brmn(jll:nsl,mparity)*sinmumi(i,m)
            work1(:,7) = work1(:,7) + tempz(:)*sinmui(i,m)
     &                              + bzmn(jll:nsl,mparity)*cosmumi(i,m)
            work1(:,11)= work1(:,11)+ blmn(jll:nsl,mparity)*cosmumi(i,m)
 
            IF (.not.lthreed) CYCLE

            work1(:,2) = work1(:,2) - crmn(jll:nsl,mparity)*cosmui(i,m)
            work1(:,3) = work1(:,3) + tempr(:)*sinmui(i,m)
     &                              + brmn(jll:nsl,mparity)*cosmumi(i,m)
            work1(:,4) = work1(:,4) - crmn(jll:nsl,mparity)*sinmui(i,m)
            work1(:,5) = work1(:,5) + tempz(:)*cosmui(i,m)
     &                              + bzmn(jll:nsl,mparity)*sinmumi(i,m)
            work1(:,6) = work1(:,6) - czmn(jll:nsl,mparity)*cosmui(i,m)
            work1(:,8) = work1(:,8) - czmn(jll:nsl,mparity)*sinmui(i,m)

            work1(:,9) = work1(:,9) + blmn(jll:nsl,mparity)*sinmumi(i,m)
            work1(:,10) =work1(:,10)- clmn(jll:nsl,mparity)*cosmui(i,m)
            work1(:,12) =work1(:,12)- clmn(jll:nsl,mparity)*sinmui(i,m)
         END DO

!
!        NEXT, DO ZETA (V) TRANSFORM
         mj = m + joff
         j2 = jmin2(m)
         jl = jlam(m)

         DO n = 0, ntor
            ni = n+ioff
            l = 0
            DO k = 1, nzeta
               j2l = j2 + l
               jmaxl = jmax + l
               jll = jl + l
               nsl = ns + l
               l = l + ns
               frcc(j2:jmax,ni,mj) = frcc(j2:jmax,ni,mj)
     &                             + work1(j2l:jmaxl,1)*cosnv(k,n)
               fzsc(j2:jmax,ni,mj) = fzsc(j2:jmax,ni,mj)
     &                             + work1(j2l:jmaxl,7)*cosnv(k,n)
               flsc(jl:ns,ni,mj) = flsc(jl:ns,ni,mj)
     &                           + work1(jll:nsl,11)*cosnv(k,n)

               IF (.not.lthreed) CYCLE

               frcc(j2:jmax,ni,mj) = frcc(j2:jmax,ni,mj)
     &                             + work1(j2l:jmaxl,2)*sinnvn(k,n)
               fzsc(j2:jmax,ni,mj) = fzsc(j2:jmax,ni,mj)
     &                             + work1(j2l:jmaxl,8)*sinnvn(k,n)
               frss(j2:jmax,ni,mj) = frss(j2:jmax,ni,mj)
     &                             + work1(j2l:jmaxl,3)*sinnv(k,n)
     &                             + work1(j2l:jmaxl,4)*cosnvn(k,n)
               fzcs(j2:jmax,ni,mj) = fzcs(j2:jmax,ni,mj)
     &                             + work1(j2l:jmaxl,5)*sinnv(k,n)
     &                             + work1(j2l:jmaxl,6)*cosnvn(k,n)

               flsc(jl:ns,ni,mj) = flsc(jl:ns,ni,mj)
     &                           + work1(jll:nsl,12)*sinnvn(k,n)
               flcs(jl:ns,ni,mj) = flcs(jl:ns,ni,mj)
     &                           + work1(jll:nsl,9)*sinnv(k,n)
     &                           + work1(jll:nsl,10)*cosnvn(k,n)
            END DO
         END DO
      END DO

!
!     COMPUTE IOTA EVOLUTION EQUATION [STORED IN LMNSC(0,0) COMPONENT]
!
#if defined(CHI_FORCE)
      IF (ictrl_prec2d.gt.0 .and. ncurr.eq.1) THEN
         ni = 0 + ioff
         mj = 0 + joff
         t1 = r0scale
         DO jl = 2, ns
            flsc(jl, ni, mj) = -t1*(buco(jl) - icurv(jl))
         END DO
      END IF
#endif
!
!     MAKE R,Z(m=1,n=0) SATISFY AVERAGE FORCE BALANCE EXACTLY
!     NOTE: for m=1, FR ~ Z1*(f0 + f2), FZ ~ R1*(f0 - f2), WHERE
!     f0 is the m=0 component of frho, f2 is m=2 component.
      IF (lforbal) THEN
         mj = 1 + joff
         ni = 0 + ioff
         t1 = nscale(0)*r0scale !/4    !!v8.52
         DO jl = 2, ns - 1
            work1(jl,1) = frcc_fac(jl)*frcc(jl,ni,mj)
     &                  + fzsc_fac(jl)*fzsc(jl,ni,mj)
            frcc(jl,ni,mj) = rzu_fac(jl)*(t1*equif(jl) + work1(jl,1))
            fzsc(jl,ni,mj) = rru_fac(jl)*(t1*equif(jl) - work1(jl,1))
         END DO
      END IF

      DEALLOCATE (work1, tempr, tempz)

      END SUBROUTINE tomnsps

      SUBROUTINE tomnspa(frzl_array, armn, brmn, crmn, azmn, bzmn,
     &   czmn, blmn, clmn, arcon, azcon)
      USE vmec_main
      USE vmec_params, ONLY: jlam, jmin2, ntmax, rsc, rcs, zcc, zss
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     &   TARGET, INTENT(inout) :: frzl_array
      REAL(dp), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(in) ::
     &   armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn, arcon, azcon
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: jmax, m, mparity, i, n, k, l
      INTEGER :: ioff, joff, mj, ni, nsl, j2, j2l, jl, jll, jmaxl 
      REAL(dp), DIMENSION(:,:,:), POINTER :: 
     &           frcs, frsc, fzcc, fzss, flcc, flss
      REAL(dp), DIMENSION(ns*nzeta) :: temp1, temp3
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: work1
!-----------------------------------------------
      frsc => frzl_array(:,:,:,rsc)               !!R-SIN(mu) COS(nv)
      fzcc => frzl_array(:,:,:,zcc+ntmax)         !!Z-COS(mu) COS(nv)
      flcc => frzl_array(:,:,:,zcc+2*ntmax)       !!L-COS(mu) COS(nv)
      IF (lthreed) THEN
         frcs => frzl_array(:,:,:,rcs)            !!R-COS(mu) SIN(nv)
         fzss => frzl_array(:,:,:,zss+ntmax)      !!Z-SIN(mu) SIN(nv)
         flss => frzl_array(:,:,:,zss+2*ntmax)    !!L-SIN(mu) SIN(nv)
      END IF

      ALLOCATE (work1(ns*nzeta,12), stat=i)
      IF (i .ne. 0) THEN
         STOP 'Allocation error in VMEC tomnspa'
      END IF

      ioff = LBOUND(frsc,2)
      joff = LBOUND(frsc,3)

      jmax = ns
      IF (ivac .lt. 1) jmax = ns1
!
!     BEGIN FOURIER TRANSFORM
!
      DO m = 0, mpol1
         mparity = MOD(m,2)
         mj = m + joff
         j2 = jmin2(m)
         jl = jlam(m)
         work1 = 0
!
!        DO THETA (U) TRANSFORM FIRST
!
         DO i = 1, ntheta2
            temp1(:) = armn(:,i,mparity)
#ifndef _HBANGLE
     &               + xmpq(m,1)*arcon(:,i,mparity)
#endif
            temp3(:) = azmn(:,i,mparity)
#ifndef _HBANGLE
     &               + xmpq(m,1)*azcon(:,i,mparity)
#endif
            work1(:,3) = work1(:,3) + temp1(:)*sinmui(i,m) 
     &                              + brmn(:,i,mparity)*cosmumi(i,m)
            work1(:,5) = work1(:,5) + temp3(:)*cosmui(i,m)
     &                              + bzmn(:,i,mparity)*sinmumi(i,m)
            work1(:,9) = work1(:,9) + blmn(:,i,mparity)*sinmumi(i,m)

            IF (.not.lthreed) CYCLE

            work1(:,1) = work1(:,1) + temp1(:)*cosmui(i,m)
     &                              + brmn(:,i,mparity)*sinmumi(i,m)
            work1(:,2) = work1(:,2) - crmn(:,i,mparity)*cosmui(i,m)
            work1(:,4) = work1(:,4) - crmn(:,i,mparity)*sinmui(i,m)
            work1(:,6) = work1(:,6) - czmn(:,i,mparity)*cosmui(i,m)
            work1(:,7) = work1(:,7) + temp3(:)*sinmui(i,m)
     &                              + bzmn(:,i,mparity)*cosmumi(i,m)
            work1(:,8) = work1(:,8) - czmn(:,i,mparity)*sinmui(i,m)
            work1(:,10) = work1(:,10) - clmn(:,i,mparity)*cosmui(i,m)
            work1(:,11) = work1(:,11) + blmn(:,i,mparity)*cosmumi(i,m)
            work1(:,12) = work1(:,12) - clmn(:,i,mparity)*sinmui(i,m)
         END DO
!
!        NEXT, DO ZETA (V) TRANSFORM
!
         DO n = 0, ntor
            ni = n + ioff
            DO k = 1, nzeta
               l = ns*(k - 1)
               j2l = j2 + l
               jmaxl = jmax + l
               jll = jl + l
               nsl = ns + l
               frsc(j2:jmax,ni,mj) = frsc(j2:jmax,ni,mj)
     &                             + work1(j2l:jmaxl,3)*cosnv(k,n)
               fzcc(j2:jmax,ni,mj) = fzcc(j2:jmax,ni,mj)
     &                             + work1(j2l:jmaxl,5)*cosnv(k,n)
               flcc(jl:ns,ni,mj) = flcc(jl:ns,ni,mj)
     &                           + work1(jll:nsl,9)*cosnv(k,n)

               IF (.not.lthreed) CYCLE

               frsc(j2:jmax,ni,mj) = frsc(j2:jmax,ni,mj)
     &                             + work1(j2l:jmaxl,4)*sinnvn(k,n)
               fzcc(j2:jmax,ni,mj) = fzcc(j2:jmax,ni,mj)
     &                             + work1(j2l:jmaxl,6)*sinnvn(k,n)
               frcs(j2:jmax,ni,mj) = frcs(j2:jmax,ni,mj)
     &                             + work1(j2l:jmaxl,1)*sinnv(k,n)
     &                             + work1(j2l:jmaxl,2)*cosnvn(k,n)
               fzss(j2:jmax,ni,mj) = fzss(j2:jmax,ni,mj)
     &                             + work1(j2l:jmaxl,7)*sinnv(k,n)
     &                             + work1(j2l:jmaxl,8)*cosnvn(k,n)
               flcc(jl:ns,ni,mj) = flcc(jl:ns,ni,mj)
     &                           + work1(jll:nsl,10)*sinnvn(k,n)
               flss(jl:ns,ni,mj) = flss(jl:ns,ni,mj)
     &                           + work1(jll:nsl,11)*sinnv(k,n)
     &                           + work1(jll:nsl,12)*cosnvn(k,n)
            END DO
         END DO
      END DO

!     IF THE SYMMETRIZED MODE USED, NEED EXTRA FACTOR OF 2 
!     IF ntheta3 USED INSTEAD OF ntheta3, DO NOT NEED THIS FACTOR
!      frzl_array = 2*frzl_array

      DEALLOCATE (work1)

      END SUBROUTINE tomnspa

      END MODULE tomnsp_mod
