      SUBROUTINE weighted_surfsep(rmnc,zmns,xm,xn,nfp,mnmax,rms,unrms)
C-----------------------------------------------
C EAL 9/00
C The boundary is a FUNCTION of theta, the TARGET a FUNCTION of angle
C d=SQRT( (rb-rt)**2 + (zb-zt)**2)
C Dd= d (d)/d(angle)
C search angle ~ theta for a zero crossing of Dd, find nearest theta values
C symmetric about theta0
C find the root Dd=0 via quadratic interpolation
C IF quadratic formula fails USE linear interpolation
C on EXIT angle CONTAINS nu values which make the nu points (Rt,Zt) nearest
C to the nu points (Rb,Zb) evaluated at theta. (For each point)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE stel_kinds
      USE optim, ONLY: rbc_vv, zbs_vv, mpol_vv, ntor_vv,
     1  theta0=>theta0_bw, phi0=>phi0_bw, wtheta=>wtheta_bw,
     2  wphi=>wphi_bw, amplw=>amplw_bw, planes=>planes_bw
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(*), INTENT(in) ::
     1   rmnc, zmns, xm, xn
      REAL(rprec),  INTENT(out) :: rms, unrms
      INTEGER, INTENT(in) :: mnmax, nfp
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      REAL(rprec) , EXTERNAL :: weightfcn
C-----------------------------------------------


C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: mu = 12     !       modes - poloidal
      INTEGER, PARAMETER :: nu = 1024   !       points - poloidal
      INTEGER, PARAMETER :: nv = 16     !       slices - toroidal
      INTEGER, PARAMETER :: nplanes1 = 4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nphi_d, nphi2,
     1   ntheta_d, i, k, j, m, n, mpol_d, mpol_d1,
     2  k1, k2, ijk
      REAL(rprec) :: twopi, SUMwgt, wgtampl(3),
     1   arg, rnow, znow, v, u, xxm, xxn,
     2   cosa, sina, a, b, c, fa, fb, fc, sol,
     3  r1, r2, z1, z2, rp2, zp2, sp, val, rpnow, zpnow
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: dummy1, dummy2
      INTEGER, DIMENSION(:), ALLOCATABLE :: INDEX, mask
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rin, zin,
     1  trin, tzin, tprin,tpzin, theta, angle, wgt
      REAL(rprec) :: rmin, rmax, zmin, zmax
C-----------------------------------------------
      sp(r1,r2,z1,z2,rp2,zp2)=1/SQRT((r1-r2)*(r1-r2)+(z1-z2)*(z1-z2))*
     .(rp2*(r2-r1)+zp2*(z2-z1))
C-----------------------------------------------
C-----------------------------------------------
      twopi = 8*ATAN(1._dp)
      wgtampl=amplw
      ntheta_d=nu;nphi_d=nv
      mpol_d = mu
      mpol_d1 = mpol_d - 1
      nphi2 = 1 + nphi_d/2
      IF(.not.ALLOCATED(rin))
     .ALLOCATE(rin(ntheta_d),zin(ntheta_d))
      IF(.not.ALLOCATED(trin))
     .ALLOCATE(trin(ntheta_d),tzin(ntheta_d))
      IF(.not.ALLOCATED(tprin))
     .ALLOCATE(tprin(ntheta_d),tpzin(ntheta_d))
      IF(.not.ALLOCATED(theta))ALLOCATE(wgt(ntheta_d),
     .theta(ntheta_d),angle(ntheta_d),INDEX(ntheta_d),mask(ntheta_d))
      IF(.not.ALLOCATED(dummy1))ALLOCATE(dummy1(ntheta_d))
      IF(.not.ALLOCATED(dummy2))ALLOCATE(dummy2(ntheta_d))

10001 CONTINUE
CEAL      v=twopi*(nfp-1)*(kz-1)/nplanes1/2
      rms=0
      unrms=0
      SUMwgt=0
      DO ijk=1,nplanes1-1
      v=planes(ijk)/nfp
      trin=0;tzin=0;rin=0;zin=0;tprin=0;tpzin=0
         CALL totbdy(ntheta_d,v,
     .          rin(1:ntheta_d),
     .          zin(1:ntheta_d),
     .          rmnc,zmns,xm,xn,mnmax)
       DO j=1,ntheta_d
        u=twopi*(j-1)/(ntheta_d-1)
        theta(j)=u
        wgt(j)= weightfcn(theta(j), v, theta0(ijk), phi0,
     .   wtheta, wphi, amplw, Nfp)
        rnow=0;znow=0;rpnow=0;zpnow=0
        DO n=-ntor_vv,ntor_vv
         DO m=0,mpol_vv
          xxm=m;xxn=n
          arg=xxm*u-xxn*planes(ijk)
          cosa=COS(arg);sina=SIN(arg)
          rnow=rnow+rbc_vv(n,m)*cosa
          znow=znow+zbs_vv(n,m)*sina
          rpnow=rpnow-rbc_vv(n,m)*xxm*sina
          zpnow=zpnow+zbs_vv(n,m)*xxm*cosa
         ENDDO
        ENDDO
        trin(j)=rnow
        tzin(j)=znow
        tprin(j)=rpnow
        tpzin(j)=zpnow
       ENDDO
 ! grid search
      INDEX=(/(k,k=1,ntheta_d)/)
      i=0
      angle=0
      DO WHILE (i < SIZE(theta))
c        k=NINT(10.*ntheta_d/360.)
c        k=NINT(20.*ntheta_d/360.)      !       seems better for li383
        k=NINT(30.*ntheta_d/360.)       !       still better
        i=i+1
        mask=cshift(INDEX,-k)
        k1=mask(i)
        mask=cshift(INDEX,k)
        k2=mask(i)
        fa=sp(rin(i),trin(k2),zin(i),tzin(k2),tprin(k2),tpzin(k2))
        fb=sp(rin(i),trin(k1),zin(i),tzin(k1),tprin(k1),tpzin(k1))
        IF(fb.eq.0.)angle(i)=theta(k1)
        IF(fb.eq.0.)cycle
        val=fa/fb
        DO WHILE(val > 0 .and. k>0)
         k=k-1
         mask=cshift(INDEX,-k)
         k1=mask(i)
         mask=cshift(INDEX,k)
         k2=mask(i)
         fa=sp(rin(i),trin(k2),zin(i),tzin(k2),tprin(k2),tpzin(k2))
         fb=sp(rin(i),trin(k1),zin(i),tzin(k1),tprin(k1),tpzin(k1))
         val=fa/fb
        ENDDO   !       WHILE(val < 0 .and. k>0)
          IF(val.gt.0)then
           angle(j)=theta(j)
           CYCLE
          ENDIF
        DO WHILE(val < 0 .and. k>0)
         k=k-1
         mask=cshift(INDEX,-k)
         k1=mask(i)
         mask=cshift(INDEX,k)
         k2=mask(i)
         fa=sp(rin(i),trin(k2),zin(i),tzin(k2),tprin(k2),tpzin(k2))
         fb=sp(rin(i),trin(k1),zin(i),tzin(k1),tprin(k1),tpzin(k1))
         IF(fa.eq.0.)angle(i)=theta(k2)
         IF(fb.eq.0.)angle(i)=theta(k1)
         IF(fb.eq.0.)cycle
         IF(fa.eq.0.)cycle
         val=fa/fb
        ENDDO   !       WHILE(val < 0 .and. k>0)
        IF(fb.eq.0.)cycle
        k=k+1
        mask=cshift(INDEX,-k)
        k1=mask(i)
        mask=cshift(INDEX,k)
        k2=mask(i)
!       find zero crossing
        a=theta(k1)
        b=theta(i)
        c=theta(k2)
        IF(b < a) a=a-twopi
        IF(c < b) c=c+twopi
        IF(c < a) c=c+twopi

        fa=sp(rin(i),trin(k1),zin(i),tzin(k1),tprin(k1),tpzin(k1))
        fb=sp(rin(i),trin(i),zin(i),tzin(i),tprin(i),tpzin(i))
        fc=sp(rin(i),trin(k2),zin(i),tzin(k2),tprin(k2),tpzin(k2))
        IF(fa.eq.0.)angle(i)=a
        IF(fb.eq.0.)angle(i)=b
        IF(fc.eq.0.)angle(i)=c
        IF((fa.eq.0.).or.(fb.eq.0.).or.(fc.eq.0.))cycle

10002   IF(fc/fb < 0)then
            sol=-(c*fb-b*fc)/(fc-fb)
        ELSEIF(fa/fb < 0)then
            sol=(a*fb-b*fa)/(fb-fa)
        ELSE
            sol=theta(i)
        ENDIF   !       IF(fc/fb < 0)
10003   angle(i)=sol
      ENDDO     !       WHILE(val < 0 .and. k>0)
      WHERE(angle < 0.)angle=angle+twopi
      WHERE(angle > twopi)angle=angle-twopi

! evaluate trin,tzin on new mesh
       DO j=1,ntheta_d
        u=angle(j)
        rnow=0;znow=0
        DO n=-ntor_vv,ntor_vv
         DO m=0,mpol_vv
          xxm=m;xxn=n
          arg=xxm*u-xxn*planes(ijk)
          cosa=COS(arg);sina=SIN(arg)
          rnow=rnow+rbc_vv(n,m)*cosa
          znow=znow+zbs_vv(n,m)*sina
         ENDDO
        ENDDO
        trin(j)=rnow
        tzin(j)=znow
       ENDDO
        dummy1(1:ntheta_d)=wgt(1:ntheta_d)*SQRT(
     &  (rin(1:ntheta_d)-trin(1:ntheta_d))**2 +
     &  (zin(1:ntheta_d)-tzin(1:ntheta_d))**2
     &)
        dummy2(1:ntheta_d)=SQRT(
     &  (rin(1:ntheta_d)-trin(1:ntheta_d))**2 +
     &  (zin(1:ntheta_d)-tzin(1:ntheta_d))**2
     &)
!      rms=rms+SUM(dummy1(1:ntheta_d))/SUM(wgt)
!       need normalization to get length, but DO not remove planar
!        weighting
      rms=rms+SUM(dummy1(1:ntheta_d))/SUM(wgt)*wgtampl(ijk)/SUM(wgtampl)

      unrms=unrms+SUM(dummy2(1:ntheta_d))/ntheta_d
      SUMwgt=sumwgt+SUM(wgt)
      ENDDO     !       ijk
      ijk=ijk-1
      unrms=unrms/ijk
      rmin = MINVAL(rin(1:ntheta_d))
      rmax = MAXVAL(rin(1:ntheta_d))
      rmin=MIN(rmin,MINVAL(trin(1:ntheta_d)))
      rmax=MAX(rmax,MAXVAL(trin(1:ntheta_d)))
      zmin = MINVAL(zin(1:ntheta_d))
      zmax = MAXVAL(zin(1:ntheta_d))
      zmin=MIN(zmin,MINVAL(tzin(1:ntheta_d)))
      zmax=MAX(zmax,MAXVAL(tzin(1:ntheta_d)))
! DEALLOCATE
       IF(ALLOCATED(dummy1))DEALLOCATE(dummy1)
       IF(ALLOCATED(dummy2))DEALLOCATE(dummy2)
       IF(ALLOCATED(tzin))DEALLOCATE(tzin)
       IF(ALLOCATED(trin))DEALLOCATE(trin)
       IF(ALLOCATED(tzin))DEALLOCATE(tpzin)
       IF(ALLOCATED(trin))DEALLOCATE(tprin)
       IF(ALLOCATED(rin))DEALLOCATE(rin)
       IF(ALLOCATED(zin))DEALLOCATE(zin)
       IF(ALLOCATED(zin))DEALLOCATE(theta)
       IF(ALLOCATED(zin))DEALLOCATE(angle)
       IF(ALLOCATED(zin))DEALLOCATE(INDEX)
      RETURN
 2100 FORMAT(' RBC(',i3,',',i2,') = ',1p,e12.4,3x,
     1       ' ZBS(',i3,',',i2,') = ',e12.4)
      END SUBROUTINE weighted_surfsep
