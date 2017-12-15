      MODULE totzsp_mod
      USE vmec_main
      USE timer_sub
      IMPLICIT NONE

      INTEGER, PARAMETER, PRIVATE :: m0=0, m1=1, n0=0
      REAL(dp), ALLOCATABLE, PRIVATE :: work1(:,:,:), work2(:,:)
      REAL(dp), PRIVATE :: cosmux, sinmux

      CONTAINS

#if defined(SKS)
      SUBROUTINE totzsps_par(rzl_array, r11, ru1, rv1, z11, zu1, zv1,
     1                       lu1, lv1, rcn1, zcn1, ier_flag)
      USE vmec_params, ONLY: jmin1, jlam, ntmax, rcc, rss, zsc, zcs,
     1                       r01_bad_value_flag
      USE precon2d, ONLY: ictrl_prec2d
      USE parallel_include_module
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,0:mpol1,ns,3*ntmax),
     1   TARGET, INTENT(INOUT) :: rzl_array
      REAL(dp), DIMENSION(nzeta,ntheta3,ns,0:1),
     1   INTENT(out) :: r11, ru1,
     1   rv1, z11, zu1,  zv1, lu1, lv1, rcn1, zcn1
      INTEGER, INTENT(inout) :: ier_flag
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: n, m, mparity, k, i, l
      INTEGER :: ioff, joff, mj, ni, nsz
      INTEGER :: nsmin, nsmax, js
      REAL(dp), DIMENSION(:,:,:), POINTER ::
     1           rmncc, rmnss, zmncs, zmnsc, lmncs, lmnsc
      REAL(dp) :: tbroadon, tbroadoff
!-----------------------------------------------
      CALL second0(tffton)

      nsmin=t1lglob; nsmax=t1rglob

      rmncc=>rzl_array(:,:,:,rcc)                  !COS(mu) COS(nv)
      zmnsc=>rzl_array(:,:,:,zsc+ntmax)            !SIN(mu) COS(nv)
      lmnsc=>rzl_array(:,:,:,zsc+2*ntmax)          !SIN(mu) COS(nv)
      IF (lthreed) THEN
         rmnss=>rzl_array(:,:,:,rss)               !SIN(mu) SIN(nv)
         zmncs=>rzl_array(:,:,:,zcs+ntmax)         !COS(mu) SIN(nv)
         lmncs=>rzl_array(:,:,:,zcs+2*ntmax)       !COS(mu) SIN(nv)
      END IF
      rzl_array(:,m1,1,:) = rzl_array(:,m1,2,:)

      ioff = LBOUND(rmncc,1)
      joff = LBOUND(rmncc,2)
      IF (lthreed) CALL convert_sym_par (rmnss(:,m1+joff,:), 
     1             zmncs(:,m1+joff,:), nsmin, nsmax)

!
!     ORIGIN EXTRAPOLATION OF M=0 MODES FOR LAMBDA 
!
      IF (lthreed .AND. jlam(m0).GT.1)
     1 lmncs(:,m0+joff,1) = lmncs(:,m0+joff,2)

!
!     EVOLVE CHIPS BY FORCES IN TOMNSPS WHEN NCURR=1, ICTRL_PREC2D != 0
!
!SPH071017
#if defined(CHI_FORCE)
      IF (ncurr .EQ. 1) THEN
         IF (ictrl_prec2d .EQ. 2) THEN
            lmnsc(n0+ioff,m0+joff,nsmin:nsmax) = chips(nsmin:nsmax)
         ELSE IF (ictrl_prec2d .NE. 0) THEN
            chips(nsmin:nsmax) = lmnsc(n0+ioff,m0+joff,nsmin:nsmax)  
         END IF
      END IF
#endif
     
      ALLOCATE (work1(nzeta,12,ns), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC2000 totzsps'

      DO js = nsmin, nsmax
        r11(:,:,js,:) = 0;   ru1(:,:,js,:) = 0;   rv1(:,:,js,:) = 0
        rcn1(:,:,js,:) = 0;  zcn1(:,:,js,:) = 0
        z11(:,:,js,:) = 0;   zu1(:,:,js,:) = 0;   zv1(:,:,js,:) = 0
        lu1(:,:,js,:) = 0;   lv1(:,:,js,:) = 0
        DO m = 0, mpol1
          mparity = MOD(m,2)
          mj = m+joff
          work1(:,:,js) = 0
!
!        INVERSE TRANSFORM IN N-ZETA, FOR FIXED M
!
          DO n = 0, ntor
            ni = n+ioff
            DO k = 1, nzeta
              work1(k,1,js) = work1(k,1,js) 
     1                         + rmncc(ni,mj,js)*cosnv(k,n)
              work1(k,6,js) = work1(k,6,js)
     1                         + zmnsc(ni,mj,js)*cosnv(k,n)
              work1(k,10,js) = work1(k,10,js) 
     1                         + lmnsc(ni,mj,js)*cosnv(k,n)

              IF (.NOT.lthreed) CYCLE

              work1(k,4,js) = work1(k,4,js) 
     1                        + rmnss(ni,mj,js)*cosnvn(k,n)
              work1(k,7,js) = work1(k,7,js) 
     1                        + zmncs(ni,mj,js)*cosnvn(k,n)
              work1(k,11,js) = work1(k,11,js)
     1                        + lmncs(ni,mj,js)*cosnvn(k,n)

              work1(k,2,js) = work1(k,2,js) 
     1                        + rmnss(ni,mj,js)*sinnv(k,n)
              work1(k,5,js) = work1(k,5,js) 
     1                        + zmncs(ni,mj,js)*sinnv(k,n)
              work1(k,9,js) = work1(k,9,js) 
     1                        + lmncs(ni,mj,js)*sinnv(k,n)

              work1(k,3,js) = work1(k,3,js) 
     1                        + rmncc(ni,mj,js)*sinnvn(k,n)
              work1(k,8,js) = work1(k,8,js) 
     1                        + zmnsc(ni,mj,js)*sinnvn(k,n)
              work1(k,12,js) = work1(k,12,js)
     1                        + lmnsc(ni,mj,js)*sinnvn(k,n)
            END DO
          END DO

!
!        INVERSE TRANSFORM IN M-THETA, FOR ALL RADIAL, ZETA VALUES
!
          l = 0
          DO i = 1, ntheta2
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)

            r11(:,i,js,mparity)=r11(:,i,js,mparity)
     1                            + work1(:,1,js)*cosmu(i,m)
            ru1(:,i,js,mparity)=ru1(:,i,js,mparity)
     1                            + work1(:,1,js)*sinmum(i,m)
            rcn1(:,i,js,mparity)=rcn1(:,i,js,mparity)
     1                            + work1(:,1,js)*cosmux

            z11(:,i,js,mparity)=z11(:,i,js,mparity)  
     1                            + work1(:,6,js)*sinmu(i,m)
            zu1(:,i,js,mparity)=zu1(:,i,js,mparity)  
     1                            + work1(:,6,js)*cosmum(i,m)
            zcn1(:,i,js,mparity)=zcn1(:,i,js,mparity) 
     1                            + work1(:,6,js)*sinmux

            lu1(:,i,js,mparity)=lu1(:,i,js,mparity)  
     1                            + work1(:,10,js)*cosmum(i,m)

            IF (.not.lthreed) CYCLE

            r11(:,i,js,mparity)=r11(:,i,js,mparity)  
     1                            + work1(:,2,js)*sinmu(i,m)
            ru1(:,i,js,mparity)=ru1(:,i,js,mparity)  
     1                            + work1(:,2,js)*cosmum(i,m)
            rcn1(:,i,js,mparity)=rcn1(:,i,js,mparity) 
     1                            + work1(:,2,js)*sinmux

            rv1(:,i,js,mparity)=rv1(:,i,js,mparity)  
     1                            + work1(:,3,js)*cosmu(i,m) 
     1                            + work1(:,4,js)*sinmu(i,m)
            z11(:,i,js,mparity)=z11(:,i,js,mparity)  
     1                            + work1(:,5,js)*cosmu(i,m)

            zu1(:,i,js,mparity)=zu1(:,i,js,mparity)  
     1                            + work1(:,5,js)*sinmum(i,m)
            zcn1(:,i,js,mparity)=zcn1(:,i,js,mparity) 
     1                            + work1(:,5,js)*cosmux
            zv1(:,i,js,mparity)=zv1(:,i,js,mparity)  
     1                            + work1(:,7,js)*cosmu(i,m) 
     1                            + work1(:,8,js)*sinmu(i,m)

            lu1(:,i,js,mparity)=lu1(:,i,js,mparity)  
     1                            + work1(:,9,js)*sinmum(i,m)
            lv1(:,i,js,mparity)=lv1(:,i,js,mparity)  
     1                            - (work1(:,11,js)*cosmu(i,m) 
     1                            + work1(:,12,js)*sinmu(i,m))
          END DO
        END DO
      END DO

      DEALLOCATE (work1)

      z01(nsmin:nsmax) = zmnsc(n0+ioff,m1+joff,nsmin:nsmax)
      r01(nsmin:nsmax) = rmncc(n0+ioff,m1+joff,nsmin:nsmax)
      IF (lactive) THEN
        IF (rank.EQ.0 .AND. r01(1).EQ.zero) THEN
          ier_flag = r01_bad_value_flag
        ELSE IF (rank.EQ.0 .AND. r01(1).NE.zero) THEN
           dkappa = z01(1)/r01(1)
        END IF
        CALL second0(tbroadon)
        CALL MPI_Bcast(dkappa,1, MPI_REAL8,0,NS_COMM,MPI_ERR)
        CALL second0(tbroadoff)
        broadcast_time = broadcast_time + (tbroadoff - tbroadon)
      END IF

      CALL second0(tfftoff)
      totzsps_time = totzsps_time + (tfftoff - tffton) 
      timer(tfft) = timer(tfft) + (tfftoff - tffton)

      END SUBROUTINE totzsps_par

      SUBROUTINE totzspa_par(rzl_array, r11, ru1, rv1, z11, zu1, zv1,
     1                       lu1, lv1, rcn1, zcn1)
      USE vmec_params, ONLY: jmin1, jlam, ntmax, rcs, rsc, zcc, zss
      USE precon2d, ONLY: ictrl_prec2d
      USE parallel_include_module
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,0:mpol1,ns,3*ntmax),
     1   TARGET, INTENT(inout) :: rzl_array
      REAL(dp), DIMENSION(nzeta,ntheta3,ns,0:1), INTENT(out) ::
     1   r11, ru1, rv1, z11, zu1, zv1, lu1, lv1, rcn1, zcn1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: m, n, mparity, k, i, l, j1
      INTEGER :: ioff, joff, mj, ni
      INTEGER :: nsmin, nsmax, js
      REAL(dp), DIMENSION(:,:,:), POINTER ::
     1           rmncs, rmnsc, zmncc, zmnss, lmncc, lmnss
C-----------------------------------------------
      CALL second0(tffton)
      nsmin=t1lglob; nsmax=t1rglob

      rmnsc => rzl_array(:,:,:,rsc)               !!SIN(mu) COS(nv)
      zmncc => rzl_array(:,:,:,zcc+ntmax)         !!COS(mu) COS(nv)
      lmncc => rzl_array(:,:,:,zcc+2*ntmax)       !!COS(mu) COS(nv)
      IF (lthreed) THEN
         rmncs => rzl_array(:,:,:,rcs)               !!COS(mu) SIN(nv)
         zmnss => rzl_array(:,:,:,zss+ntmax)         !!SIN(mu) SIN(nv)
         lmnss => rzl_array(:,:,:,zss+2*ntmax)       !!SIN(mu) SIN(nv)
      END IF

!
!     CONVERT FROM INTERNAL XC REPRESENTATION FOR m=1 MODES, R+(at rsc) = .5(rsc + zcc),
!     R-(at zcc) = .5(rsc - zcc), TO REQUIRED rsc, zcc FORMS
!
      ioff = LBOUND(rmnsc,1)
      joff = LBOUND(rmnsc,2)
      CALL convert_asym_par (rmnsc(:,m1+joff,:), 
     1                       zmncc(:,m1+joff,:), nsmin, nsmax)

      z00b = zmncc(ioff,joff,ns)

      ALLOCATE (work1(nzeta,12,ns), stat=i)
      IF (i .NE. 0) STOP 'Allocation error in VMEC totzspa'

!
!     INITIALIZATION BLOCK
!

      IF (jlam(m0) .gt. 1) lmncc(:,m0+joff,1) = lmncc(:,m0+joff,2)

      DO js=nsmin, nsmax
        r11(:,:,js,:) = 0;   ru1(:,:,js,:) = 0;   rv1(:,:,js,:) = 0
        rcn1(:,:,js,:) = 0;  zcn1(:,:,js,:) = 0
        z11(:,:,js,:) = 0;   zu1(:,:,js,:) = 0;   zv1(:,:,js,:) = 0
        lu1(:,:,js,:) = 0;   lv1(:,:,js,:) = 0
        DO m = 0, mpol1
          mparity = MOD(m,2)
          mj = m+joff
          work1(:,:,js) = 0
          j1 = jmin1(m)

          DO n = 0, ntor
            ni = n+ioff
            DO k = 1, nzeta
              work1(k,1,js) = work1(k,1,js)
     1                          + rmnsc(ni,mj,js)*cosnv(k,n)
              work1(k,6,js) = work1(k,6,js) 
     1                          + zmncc(ni,mj,js)*cosnv(k,n)
              work1(k,10,js) = work1(k,10,js)
     1                          + lmncc(ni,mj,js)*cosnv(k,n)

              IF (.NOT.lthreed) CYCLE

              work1(k,2,js) = work1(k,2,js)
     1                          + rmncs(ni,mj,js)*sinnv(k,n)
              work1(k,3,js) = work1(k,3,js) 
     1                          + rmnsc(ni,mj,js)*sinnvn(k,n)
              work1(k,4,js) = work1(k,4,js) 
     1                          + rmncs(ni,mj,js)*cosnvn(k,n)
              work1(k,5,js) = work1(k,5,js) 
     1                          + zmnss(ni,mj,js)*sinnv(k,n)
              work1(k,7,js) = work1(k,7,js)
     1                          + zmnss(ni,mj,js)*cosnvn(k,n)
              work1(k,8,js) = work1(k,8,js) 
     1                          + zmncc(ni,mj,js)*sinnvn(k,n)
              work1(k,9,js) = work1(k,9,js) 
     1                          + lmnss(ni,mj,js)*sinnv(k,n)
              work1(k,11,js) = work1(k,11,js)
     1                          + lmnss(ni,mj,js)*cosnvn(k,n)
              work1(k,12,js) = work1(k,12,js) 
     1                          + lmncc(ni,mj,js)*sinnvn(k,n)
            END DO
          END DO

!
!        INVERSE TRANSFORM IN M-THETA
!
          DO i = 1, ntheta2
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)

            r11(:,i,js,mparity) = r11(:,i,js,mparity) + work1(:,1,js)*
     1            sinmu(i,m)
            ru1(:,i,js,mparity) = ru1(:,i,js,mparity) + work1(:,1,js)*
     1            cosmum(i,m)
            z11(:,i,js,mparity) = z11(:,i,js,mparity) + work1(:,6,js)*
     1            cosmu(i,m)
            zu1(:,i,js,mparity) = zu1(:,i,js,mparity) + work1(:,6,js)*
     1            sinmum(i,m)
            lu1(:,i,js,mparity) = lu1(:,i,js,mparity) + work1(:,10,js)*
     1            sinmum(i,m)
            rcn1(:,i,js,mparity) = rcn1(:,i,js,mparity)+ work1(:,1,js)*
     1            sinmux
            zcn1(:,i,js,mparity) = zcn1(:,i,js,mparity)+ work1(:,6,js)*
     1            cosmux

            IF (.not.lthreed) CYCLE

            r11(:,i,js,mparity) = r11(:,i,js,mparity) + work1(:,2,js)*
     1               cosmu(i,m)
            ru1(:,i,js,mparity) = ru1(:,i,js,mparity) + work1(:,2,js)*
     1               sinmum(i,m)
            z11(:,i,js,mparity) = z11(:,i,js,mparity) + work1(:,5,js)*
     1               sinmu(i,m)
            zu1(:,i,js,mparity) = zu1(:,i,js,mparity) + work1(:,5,js)*
     1               cosmum(i,m)
            lu1(:,i,js,mparity) = lu1(:,i,js,mparity) + work1(:,9,js)*
     1               cosmum(i,m)
            rcn1(:,i,js,mparity) = rcn1(:,i,js,mparity)+ work1(:,2,js)*
     1               cosmux
            zcn1(:,i,js,mparity) = zcn1(:,i,js,mparity)+ work1(:,5,js)*
     1               sinmux
            rv1(:,i,js,mparity) = rv1(:,i,js,mparity) + work1(:,3,js)*
     1               sinmu(i,m) + work1(:,4,js)*cosmu(i,m)
            zv1(:,i,js,mparity) = zv1(:,i,js,mparity) + work1(:,7,js)*
     1               sinmu(i,m) + work1(:,8,js)*cosmu(i,m)
            lv1(:,i,js,mparity) = lv1(:,i,js,mparity)-(work1(:,11,js)*
     1               sinmu(i,m)+work1(:,12,js)*cosmu(i,m))
          END DO
        END DO
      END DO

      DEALLOCATE (work1)

      CALL second0(tfftoff)
      totzspa_time = totzspa_time + (tfftoff - tffton)
      timer(tfft) = timer(tfft) + (tfftoff - tffton)

      END SUBROUTINE totzspa_par

      SUBROUTINE convert_sym_par(rmnss, zmncs, nsmin, nsmax)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: nsmin, nsmax
      REAL(dp), DIMENSION(0:ntor,ns), INTENT(INOUT) :: rmnss, zmncs
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,nsmin:nsmax) :: temp
C-----------------------------------------------
!
!     CONVERT FROM INTERNAL REPRESENTATION TO "PHYSICAL" RMNSS, ZMNCS FOURIER FORM
!     (for lconm1, rss = zmncs)
!
      IF (.NOT.lconm1) RETURN

      temp(:,nsmin:nsmax) = rmnss(:,nsmin:nsmax)
      rmnss(:,nsmin:nsmax) = temp(:,nsmin:nsmax)
     1                     + zmncs(:,nsmin:nsmax)
      zmncs(:,nsmin:nsmax) = temp(:,nsmin:nsmax)
     1                     - zmncs(:,nsmin:nsmax)

      END SUBROUTINE convert_sym_par

      SUBROUTINE convert_asym_par(rmnsc, zmncc, nsmin, nsmax)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: nsmin, nsmax
      REAL(dp), DIMENSION(0:ntor,ns), INTENT(INOUT) :: 
     1                                           rmnsc, zmncc
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,nsmin:nsmax) :: temp
C-----------------------------------------------
!
!     CONVERT FROM INTERNAL REPRESENTATION TO RMNSC, ZMNCC FOURIER FORM
!
      IF (.NOT.lconm1) RETURN

      temp(:,nsmin:nsmax) = rmnsc(:,nsmin:nsmax)
      rmnsc(:,nsmin:nsmax) = temp(:,nsmin:nsmax) 
     1                     + zmncc(:,nsmin:nsmax)
      zmncc(:,nsmin:nsmax) = temp(:,nsmin:nsmax) 
     1                     - zmncc(:,nsmin:nsmax)

      END SUBROUTINE convert_asym_par
#endif

      SUBROUTINE totzsps(rzl_array, r11, ru1, rv1, z11, zu1, zv1,
     1                   lu1, lv1, rcn1, zcn1)
      USE vmec_params, ONLY: jmin1, jlam, ntmax, rcc, rss, zsc, zcs
      USE precon2d, ONLY: ictrl_prec2d
#if defined(SKS)
      USE realspace
      USE vforces
      USE parallel_include_module
#endif
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(INOUT) :: rzl_array
      REAL(dp), DIMENSION(ns*nzeta*ntheta3,0:1),
     1   INTENT(OUT) :: r11, ru1,
     1   rv1, z11, zu1,  zv1, lu1, lv1, rcn1, zcn1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: n, m, mparity, k, i, j1, l, j1l, nsl
      INTEGER :: ioff, joff, mj, ni, nsz
      REAL(dp), DIMENSION(:,:,:), POINTER ::
     1           rmncc, rmnss, zmncs, zmnsc, lmncs, lmnsc
!-----------------------------------------------
!
!     WORK1 Array of inverse transforms in toroidal angle (zeta), for all radial positions
!     NOTE: ORDERING OF LAST INDEX IS DIFFERENT HERE THAN IN PREVIOUS VMEC2000 VERSIONS
!
!     CONVERT FROM INTERNAL XC REPRESENTATION FOR m=1 MODES, R+(stored at rss) = .5(rss + zcs),
!     R-(stored at zcs) = .5(rss - zcs), TO EXTERNAL ("PHYSICAL") rss, zcs FORMS. NEED THIS EVEN 
!     WHEN COMPUTING HESSIAN FOR FREE BOUNDARY (rmnss, zmncs at JS=NS needed in vacuum call)
!
      CALL second0(tffton)

      rmncc => rzl_array(:,:,:,rcc)               !!COS(mu) COS(nv)
      zmnsc => rzl_array(:,:,:,zsc+ntmax)         !!SIN(mu) COS(nv)
      lmnsc => rzl_array(:,:,:,zsc+2*ntmax)       !!SIN(mu) COS(nv)
      IF (lthreed) THEN
         rmnss => rzl_array(:,:,:,rss)               !!SIN(mu) SIN(nv)
         zmncs => rzl_array(:,:,:,zcs+ntmax)         !!COS(mu) SIN(nv)
         lmncs => rzl_array(:,:,:,zcs+2*ntmax)       !!COS(mu) SIN(nv)
      END IF

      ioff = LBOUND(rmncc,2)
      joff = LBOUND(rmncc,3)
#ifndef _HBANGLE
      IF (lthreed) 
     1   CALL convert_sym (rmnss(:,:,m1+joff), zmncs(:,:,m1+joff))
#endif

!v8.50: Norm for preconditioned R,Z forces: scale to boundary value only
!v8.51  Restore hs dependence (1:ns, not just ns)
!     fnorm1 = one/SUM(rzl_array(1:ns,:,m1:,1:2*ntmax)**2)    

!
!     ORIGIN EXTRAPOLATION (JS=1) FOR M=1 MODES
!     NOTE: PREVIOUS VERSIONS OF VMEC USED TWO-POINT EXTRAPOLATION 
!           FOR R,Z. HOWEVER,THIS CAN NOT BE USED TO COMPUTE THE 
!           TRI-DIAG 2D PRECONDITIONER
!
      rzl_array(1,:,m1,:)  = rzl_array(2,:,m1,:)

!
!     ORIGIN EXTRAPOLATION OF M=0 MODES FOR LAMBDA 
!
      IF (lthreed .and. jlam(m0).gt.1) 
     1   lmncs(1,:,m0+joff) = lmncs(2,:,m0+joff)

!
!     EVOLVE CHIPS BY FORCES IN TOMNSPS WHEN NCURR=1, ICTRL_PREC2D != 0
!
#if defined(CHI_FORCE)
      IF (ncurr .EQ. 1) THEN
         IF (ictrl_prec2d .EQ. 2) THEN
            lmnsc(2:ns,n0+ioff,m0+joff) = chips(2:ns)
         ELSE IF (ictrl_prec2d .NE. 0) THEN
            chips(2:ns) = lmnsc(2:ns,n0+ioff,m0+joff)  
         END IF
      END IF
#endif     
      nsz = ns*nzeta
      ALLOCATE (work2(nsz,12), stat=i)
      IF (i .NE. 0) STOP 'Allocation error in VMEC2000 totzsps'

      r11 = 0;  ru1 = 0;  rv1 = 0;  rcn1 = 0
      z11 = 0;  zu1 = 0;  zv1 = 0;  zcn1 = 0
      lu1 = 0;  lv1 = 0

!
!     COMPUTE R, Z, AND LAMBDA IN REAL SPACE
!     NOTE: LU = d(Lam)/du, LV = -d(Lam)/dv
!

      DO m = 0, mpol1
         mparity = MOD(m,2)
         mj = m+joff
         work2 = 0
         j1 = jmin1(m)
!
!        INVERSE TRANSFORM IN N-ZETA, FOR FIXED M
!
         DO n = 0, ntor
            ni = n+ioff
            DO k = 1, nzeta
               l = ns*(k-1)
               j1l = j1+l;  nsl = ns+l
               work2(j1l:nsl,1) = work2(j1l:nsl,1) 
     1                          + rmncc(j1:ns,ni,mj)*cosnv(k,n)
               work2(j1l:nsl,6) = work2(j1l:nsl,6)
     1                          + zmnsc(j1:ns,ni,mj)*cosnv(k,n)
               work2(j1l:nsl,10) = work2(j1l:nsl,10) 
     1                          + lmnsc(j1:ns,ni,mj)*cosnv(k,n)

               IF (.not.lthreed) CYCLE
               
               work2(j1l:nsl,4) = work2(j1l:nsl,4) 
     1                          + rmnss(j1:ns,ni,mj)*cosnvn(k,n)
               work2(j1l:nsl,7) = work2(j1l:nsl,7) 
     1                          + zmncs(j1:ns,ni,mj)*cosnvn(k,n)
               work2(j1l:nsl,11) = work2(j1l:nsl,11)
     1                          + lmncs(j1:ns,ni,mj)*cosnvn(k,n)

               work2(j1l:nsl,2) = work2(j1l:nsl,2) 
     1                          + rmnss(j1:ns,ni,mj)*sinnv(k,n)
               work2(j1l:nsl,5) = work2(j1l:nsl,5) 
     1                          + zmncs(j1:ns,ni,mj)*sinnv(k,n)
               work2(j1l:nsl,9) = work2(j1l:nsl,9) 
     1                          + lmncs(j1:ns,ni,mj)*sinnv(k,n)

               work2(j1l:nsl,3) = work2(j1l:nsl,3) 
     1                          + rmncc(j1:ns,ni,mj)*sinnvn(k,n)
               work2(j1l:nsl,8) = work2(j1l:nsl,8) 
     1                          + zmnsc(j1:ns,ni,mj)*sinnvn(k,n)
               work2(j1l:nsl,12) = work2(j1l:nsl,12)
     1                          + lmnsc(j1:ns,ni,mj)*sinnvn(k,n)
            END DO
         END DO
!
!        INVERSE TRANSFORM IN M-THETA, FOR ALL RADIAL, ZETA VALUES
!
         l = 0
         DO i = 1, ntheta2
            j1l = l+1;  nsl = nsz+l
            l = l + nsz
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)
           
            r11(j1l:nsl,mparity)  = r11(j1l:nsl,mparity)  
     1                            + work2(1:nsz,1)*cosmu(i,m)
            ru1(j1l:nsl,mparity)  = ru1(j1l:nsl,mparity)  
     1                            + work2(1:nsz,1)*sinmum(i,m)
#ifndef _HBANGLE
            rcn1(j1l:nsl,mparity) = rcn1(j1l:nsl,mparity) 
     1                            + work2(1:nsz,1)*cosmux
#endif

            z11(j1l:nsl,mparity)  = z11(j1l:nsl,mparity)  
     1                            + work2(1:nsz,6)*sinmu(i,m)
            zu1(j1l:nsl,mparity)  = zu1(j1l:nsl,mparity)  
     1                            + work2(1:nsz,6)*cosmum(i,m)
#ifndef _HBANGLE
            zcn1(j1l:nsl,mparity) = zcn1(j1l:nsl,mparity) 
     1                            + work2(1:nsz,6)*sinmux
#endif
             
            lu1(j1l:nsl,mparity)  = lu1(j1l:nsl,mparity)  
     1                            + work2(1:nsz,10)*cosmum(i,m)

            IF (.not.lthreed) CYCLE

            r11(j1l:nsl,mparity)  = r11(j1l:nsl,mparity)  
     1                            + work2(1:nsz,2)*sinmu(i,m)
            ru1(j1l:nsl,mparity)  = ru1(j1l:nsl,mparity)  
     1                            + work2(1:nsz,2)*cosmum(i,m)
#ifndef _HBANGLE
            rcn1(j1l:nsl,mparity) = rcn1(j1l:nsl,mparity) 
     1                            + work2(1:nsz,2)*sinmux
#endif
     
            rv1(j1l:nsl,mparity)  = rv1(j1l:nsl,mparity)  
     1                            + work2(1:nsz,3)*cosmu(i,m) 
     1                            + work2(1:nsz,4)*sinmu(i,m)
            z11(j1l:nsl,mparity)  = z11(j1l:nsl,mparity)  
     1                            + work2(1:nsz,5)*cosmu(i,m)

            zu1(j1l:nsl,mparity)  = zu1(j1l:nsl,mparity)  
     1                            + work2(1:nsz,5)*sinmum(i,m)
#ifndef _HBANGLE
            zcn1(j1l:nsl,mparity) = zcn1(j1l:nsl,mparity) 
     1                            + work2(1:nsz,5)*cosmux
#endif
            zv1(j1l:nsl,mparity)  = zv1(j1l:nsl,mparity)  
     1                            + work2(1:nsz,7)*cosmu(i,m) 
     1                            + work2(1:nsz,8)*sinmu(i,m)

            lu1(j1l:nsl,mparity)  = lu1(j1l:nsl,mparity)  
     1                            + work2(1:nsz,9)*sinmum(i,m)
            lv1(j1l:nsl,mparity)  = lv1(j1l:nsl,mparity)  
     1                            - (work2(1:nsz,11)*cosmu(i,m) 
     1                            + work2(1:nsz,12)*sinmu(i,m))
         END DO
      END DO

      DEALLOCATE (work2)

      z01(1:ns) = zmnsc(1:ns,n0+ioff,m1+joff)
      r01(1:ns) = rmncc(1:ns,n0+ioff,m1+joff)
      IF (r01(1) .eq. zero) STOP 'r01(0) = 0 in totzsps_SPH'
      dkappa = z01(1)/r01(1)

      CALL second0(tfftoff)
      timer(tfft) = timer(tfft) + (tfftoff - tffton)
#if defined(SKS)
      s_totzsps_time=timer(tfft)
#endif

      END SUBROUTINE totzsps

      SUBROUTINE convert_sym(rmnss, zmncs)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor), INTENT(INOUT) :: rmnss, zmncs
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor) :: temp
C-----------------------------------------------
!
!     CONVERT FROM INTERNAL REPRESENTATION TO "PHYSICAL" RMNSS, ZMNCS FOURIER FORM
!
      IF (.NOT.lconm1) RETURN

      temp = rmnss
      rmnss = temp + zmncs
      zmncs = temp - zmncs

      END SUBROUTINE convert_sym


      SUBROUTINE totzspa(rzl_array, r11, ru1, rv1, z11, zu1, zv1, lu1,
     1   lv1, rcn1, zcn1)
      USE vmec_params, ONLY: jmin1, jlam, ntmax, rcs, rsc, zcc, zss
      USE precon2d, ONLY: ictrl_prec2d
#if defined(SKS)
      USE parallel_include_module
#endif
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(inout) :: rzl_array
      REAL(dp), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(out) ::
     1   r11, ru1, rv1, z11, zu1, zv1, lu1, lv1, rcn1, zcn1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: m, n, mparity, k, i, l, j1, j1l, nsl
      INTEGER :: ioff, joff, mj, ni
      REAL(dp), DIMENSION(:,:,:), POINTER ::
     1           rmncs, rmnsc, zmncc, zmnss, lmncc, lmnss
C-----------------------------------------------
      CALL second0(tffton)
      rmnsc => rzl_array(:,:,:,rsc)               !!SIN(mu) COS(nv)
      zmncc => rzl_array(:,:,:,zcc+ntmax)         !!COS(mu) COS(nv)
      lmncc => rzl_array(:,:,:,zcc+2*ntmax)       !!COS(mu) COS(nv)
      IF (lthreed) THEN
         rmncs => rzl_array(:,:,:,rcs)               !!COS(mu) SIN(nv)
         zmnss => rzl_array(:,:,:,zss+ntmax)         !!SIN(mu) SIN(nv)
         lmnss => rzl_array(:,:,:,zss+2*ntmax)       !!SIN(mu) SIN(nv)
      END IF

!
!     CONVERT FROM INTERNAL XC REPRESENTATION FOR m=1 MODES, R+(at rsc) = .5(rsc + zcc),
!     R-(at zcc) = .5(rsc - zcc), TO REQUIRED rsc, zcc FORMS
!
      ioff = LBOUND(rmnsc,2)
      joff = LBOUND(rmnsc,3)
      CALL convert_asym (rmnsc(:,:,m1+joff), zmncc(:,:,m1+joff))

      z00b = zmncc(ns,ioff,joff)

      IF (jlam(m0) .GT. 1) lmncc(1,:,m0+joff) = lmncc(2,:,m0+joff)

      IF (.FALSE.) RETURN
!      IF (ictrl_prec2d .eq. 3) RETURN

      ALLOCATE (work2(ns*nzeta,12), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC totzspa'

!
!     INITIALIZATION BLOCK
!
      r11 = 0;  ru1 = 0;  rv1 = 0;  z11 = 0;  zu1 = 0
      zv1 = 0;  lu1 = 0;  lv1 = 0;  rcn1 = 0; zcn1 = 0

      DO m = 0, mpol1
         mparity = MOD(m,2)
         mj = m+joff
         work2 = 0
         j1 = jmin1(m)
         DO n = 0, ntor
            ni = n+ioff
            DO k = 1, nzeta
               l = ns*(k-1)
               j1l = j1+l;  nsl = ns+l
               work2(j1l:nsl,1) = work2(j1l:nsl,1)
     1                          + rmnsc(j1:ns,ni,mj)*cosnv(k,n)
               work2(j1l:nsl,6) = work2(j1l:nsl,6) 
     1                          + zmncc(j1:ns,ni,mj)*cosnv(k,n)
               work2(j1l:nsl,10) = work2(j1l:nsl,10)
     1                          + lmncc(j1:ns,ni,mj)*cosnv(k,n)

               IF (.not.lthreed) CYCLE

               work2(j1l:nsl,2) = work2(j1l:nsl,2)
     1                          + rmncs(j1:ns,ni,mj)*sinnv(k,n)
               work2(j1l:nsl,3) = work2(j1l:nsl,3) 
     1                          + rmnsc(j1:ns,ni,mj)*sinnvn(k,n)
               work2(j1l:nsl,4) = work2(j1l:nsl,4) 
     1                          + rmncs(j1:ns,ni,mj)*cosnvn(k,n)
               work2(j1l:nsl,5) = work2(j1l:nsl,5) 
     1                          + zmnss(j1:ns,ni,mj)*sinnv(k,n)
               work2(j1l:nsl,7) = work2(j1l:nsl,7)
     1                          + zmnss(j1:ns,ni,mj)*cosnvn(k,n)
               work2(j1l:nsl,8) = work2(j1l:nsl,8) 
     1                          + zmncc(j1:ns,ni,mj)*sinnvn(k,n)
               work2(j1l:nsl,9) = work2(j1l:nsl,9) 
     1                          + lmnss(j1:ns,ni,mj)*sinnv(k,n)
               work2(j1l:nsl,11) = work2(j1l:nsl,11)
     1                          + lmnss(j1:ns,ni,mj)*cosnvn(k,n)
               work2(j1l:nsl,12) = work2(j1l:nsl,12) 
     1                          + lmncc(j1:ns,ni,mj)*sinnvn(k,n)
            END DO
         END DO

!
!        INVERSE TRANSFORM IN M-THETA
!
         DO i = 1, ntheta2
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)
            r11(:,i,mparity) = r11(:,i,mparity) + work2(:,1)*
     1            sinmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work2(:,1)*
     1            cosmum(i,m)
            z11(:,i,mparity) = z11(:,i,mparity) + work2(:,6)*
     1            cosmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work2(:,6)*
     1            sinmum(i,m)
            lu1(:,i,mparity) = lu1(:,i,mparity) + work2(:,10)*
     1            sinmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work2(:,1)*
     1            sinmux
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work2(:,6)*
     1            cosmux

            IF (.not.lthreed) CYCLE
               
            r11(:,i,mparity) = r11(:,i,mparity) + work2(:,2)*
     1               cosmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work2(:,2)*
     1               sinmum(i,m)
            z11(:,i,mparity) = z11(:,i,mparity) + work2(:,5)*
     1               sinmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work2(:,5)*
     1               cosmum(i,m)
            lu1(:,i,mparity) = lu1(:,i,mparity) + work2(:,9)*
     1               cosmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work2(:,2)*
     1               cosmux
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work2(:,5)*
     1               sinmux
            rv1(:,i,mparity) = rv1(:,i,mparity) + work2(:,3)*
     1               sinmu(i,m) + work2(:,4)*cosmu(i,m)
            zv1(:,i,mparity) = zv1(:,i,mparity) + work2(:,7)*
     1               sinmu(i,m) + work2(:,8)*cosmu(i,m)
            lv1(:,i,mparity) = lv1(:,i,mparity)-(work2(:,11)*
     1               sinmu(i,m)+work2(:,12)*cosmu(i,m))
         END DO
      END DO

      DEALLOCATE (work2)
 
      CALL second0(tfftoff)
      timer(tfft) = timer(tfft) + (tfftoff - tffton)
#if defined(SKS)
	s_totzspa_time=timer(tfft)
#endif
      END SUBROUTINE totzspa


      SUBROUTINE convert_asym(rmnsc, zmncc)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor), INTENT(INOUT) :: rmnsc, zmncc
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor) :: temp
C-----------------------------------------------
!
!     CONVERT FROM INTERNAL REPRESENTATION TO RMNSC, ZMNCC FOURIER FORM
!
      IF (.NOT.lconm1) RETURN

      temp = rmnsc
      rmnsc = temp + zmncc
      zmncc = temp - zmncc

      END SUBROUTINE convert_asym

      END MODULE totzsp_mod
