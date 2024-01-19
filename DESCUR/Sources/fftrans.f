      SUBROUTINE fftrans(r0c, z0c, rhoc, rhos, rbc, zbs, rbs, zbc,
     1   rmnaxis, zmnaxis)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vname0
      USE Vname1
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(nphi) :: r0c, z0c
      REAL(rprec), DIMENSION(0:mrho-1,nphi) :: rhoc, rhos
      REAL(rprec), DIMENSION(0:mpol-1,-nphi2:nphi2) ::
     1  rbc, zbs, rbs, zbc
      REAL(rprec), DIMENSION(0:nphi2) :: rmnaxis, zmnaxis
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, mn, mreal, nreal
      REAL(rprec), DIMENSION(nv) :: intgrate, argi
      REAL(rprec) :: delphi,dn,rmc_p,zms_p,rms_p,zmc_p,
     1  arg,tcosn,tsinn
C-----------------------------------------------
!
!       PERFORM FOURIER TRANSFORM IN phi
!
      delphi = one/nphi
      DO i = 1, nphi
         intgrate(i) = delphi
         argi(i) = twopi*(i - 1)/REAL(nphi*nfp,rprec)
      END DO

      rbc = 0;   zbs = 0;   rbs = 0;   zbc = 0

      DO mn = 1, mpnt
         mreal = m1(mn)
         nreal = n1(mn)/nfp
         dn = REAL(n1(mn))
         DO i = 1, nphi
            CALL getrz(rmc_p,rms_p,zmc_p,zms_p,r0c(i),z0c(i),
     1         rhoc(0,i),rhos(0,i),mreal,mrho)
            arg = dn*argi(i)
            tcosn = COS(arg)
            tsinn = SIN(arg)
            rbc(mreal,nreal) = rbc(mreal,nreal) + intgrate(i)*(tcosn*
     1         rmc_p + tsinn*rms_p)
            zbs(mreal,nreal) = zbs(mreal,nreal) + intgrate(i)*(tcosn*
     1         zms_p - tsinn*zmc_p)
            zbc(mreal,nreal) = zbc(mreal,nreal) + intgrate(i)*(tcosn*
     1         zmc_p + tsinn*zms_p)
            rbs(mreal,nreal) = rbs(mreal,nreal) + intgrate(i)*(tcosn*
     1         rms_p - tsinn*rmc_p)
     
!            zbc(mreal, nreal) = 0
!            rbs(mreal, nreal) = 0
     
         END DO
         IF (mreal.eq.0 .and. nreal.eq.0) THEN
            rmnaxis(0) = DOT_PRODUCT(intgrate(:nphi),raxis(:nphi))
            zmnaxis(0) = zero
         ELSE IF (mreal.eq.0 .and. nreal.gt.0) THEN
            rbc(0,nreal) = 2*rbc(0,nreal)
            rbs(0,nreal) = 2*rbs(0,nreal)
            zbc(0,nreal) = 2*zbc(0,nreal)
            zbs(0,nreal) = 2*zbs(0,nreal)
            rmnaxis(nreal) = zero
            zmnaxis(nreal) = zero
            rmnaxis(nreal) = rmnaxis(nreal) + SUM(2*intgrate(:nphi)*
     1         raxis(:nphi)*COS(dn*argi(:nphi)))
            zmnaxis(nreal) = zmnaxis(nreal) - SUM(2*intgrate(:nphi)*
     1         zaxis(:nphi)*SIN(dn*argi(:nphi)))
         ENDIF
      END DO

      END SUBROUTINE fftrans
