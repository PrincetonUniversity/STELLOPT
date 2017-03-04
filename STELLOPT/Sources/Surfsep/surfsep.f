      SUBROUTINE surfsep(mnmax_pl, r_pl, z_pl, xm_pl, xn_pl, mnmax_vv,
     1    r_vv, z_vv, xm_vv, xn_vv, nper, igrid, distance, drms, ierr)
! *************************************************************************
!     SURFSEP SUBS SUITE - nested surface separation SUBROUTINEs
!     modified from SURFSEP4 by A. Brooks to compute the sign of
!     the distance to determine IF the point on the free surface
!     lies inside the fixed surface Note that this could ONLY be
!     compiled with optimization at -O 1 when it was a standalone
!     routine. This problem does not exist in this version so we
!     compile with -O (= -O 2) In this set of routines 2 refers to
!     the plasma surface while 1 refers to the vacuum vessel or
!     "fixed" surface.
!     R.E. Hatcher - May 2000
! *************************************************************************
      USE stel_kinds
      USE cmnf1
      USE cmnf2
      USE geom
      USE optim_params, ONLY: nu_vv, nv_vv

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER mnmax_pl, mnmax_vv, nper, igrid, ierr
      REAL(rprec) distance, drms
      REAL(rprec), DIMENSION(mnmax_pl) :: r_pl, z_pl, xm_pl, xn_pl
      REAL(rprec), DIMENSION(mnmax_vv) :: r_vv, z_vv, xm_vv, xn_vv
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: iperr,ngood2,i,isurf,ip,it,np,nt,ngood
      REAL(rprec) :: xun, yun, zun, davg, xvec, yvec, zvec, xf, yf,
     1   zf, r1, z1, d, DOtp, du, dv, phi, u1, v1, u2, v2, r2, z2,
     2   dmean, std
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: dmat

      INTEGER :: jp, jt, kp, kt
      REAL(rprec) :: u1_sav, v1_sav, dist_test

!-----------------------------------------------
!
      ierr = 0                                   ! presume success to start
!
!     define some angular constants
      q1 = nper
      q2 = nper
      pi2 = 4*ASIN(1._dp)
      alp = pi2/nper
!
!     set np and nt to igrid
      np = igrid
      nt = igrid
!
!     ALLOCATE an array for distances
      ALLOCATE(dmat(np*nt), stat=ierr)
!
!     compute the total number of modes for each surface
      nmn2 = mnmax_pl
      nmn1 = mnmax_vv
!
!     equate the common block variables with the arguments
!     Surface 1 - Vacuum Vessel
!
      ALLOCATE (xm1(nmn1), xn1(nmn1), rmn1(nmn1), zmn1(nmn1),
     1          xm2(nmn2), xn2(nmn2), rmn2(nmn2), zmn2(nmn2),stat=ierr)
      IF (ierr .ne. 0) RETURN
!
!
      xm1(:nmn1) = xm_vv(:nmn1)
      xn1(:nmn1) = xn_vv(:nmn1)
      rmn1(:nmn1) = r_vv(:nmn1)
      zmn1(:nmn1) = z_vv(:nmn1)
!
!     equate the common block variables with the arguments
!     Surface 2 - Plasma
!
      xm2(:nmn2) = xm_pl(:nmn2)
      xn2(:nmn2) = xn_pl(:nmn2)
      rmn2(:nmn2) = r_pl(:nmn2)
      zmn2(:nmn2) = z_pl(:nmn2)
!
!     WRITE (6,*) 'done copying data', nmn2,nmn1
!
!     foreach each point on plasfree,calculate dist to plas fix
      distance = 1.0e30_dp
      ngood = 0

      du = 1.0_dp/np
      dv = 1.0_dp/nt
      DO ip = 1, np
         DO it = 1, nt
            u2 = (ip - 1)*du
            v2 = (it - 1)*dv
!           ! get point (r2,z2) on plasma surface corrseponding to u2,v2
            CALL dmnf2 (u2, v2, r2, z2)
            phi = alp*v2
            xp = r2*COS(phi)
            yp = r2*SIN(phi)
            zp = z2
!           ! assume that u2,v2 is a good starting estimate for u1,v1
            u1 = u2
            v1 = v2

!           un-vectorized loop to find a better starting u1, v1       LPK-011002

            u1_sav = u1
            v1_sav = v1
            dist_test=1.0e30_dp

!           The following neighborhood SIZE should be made as part of input
!
            DO jt=-nv_vv, nv_vv
            DO jp=-nu_vv, nu_vv
            kt = it + jt
            kp = ip + jp
            IF ( kt .le. 0 ) kt = kt + nt
            IF ( kp .le. 0 ) kp = kp + np
            u1 = (kp -1)*du
            v1 = (kt -1)*dv
            CALL dmnf1( u1, v1, r1, z1 )
            phi = alp*v1
            xf = r1*COS(phi)
            yf = r1*SIN(phi)
            zf = z1
            d  = (xf-xp)**2+(yf-yp)**2+(zf-zp)**2
            IF( d .lt. dist_test ) THEN
            u1_sav = u1
            v1_sav = v1
            dist_test = d
            END IF
            END DO
            END DO
            u1 = u1_sav
            v1 = v1_sav
            d  = dist_test
!
!           END initial search of minimum distance

!           ! find the point of a minimum distance from the plasma
!           ! to the vacuum vessel
            iperr = 0
            CALL p2surf (u1, v1, d, iperr)
            IF (iperr == 1) CYCLE                ! could not find a minimum distance
!           ! check for whether free surface lies inside
!           ! the fixed surface by checking the sign of the
!           ! dot product between the normal and the vector
!           ! from (u2,v2) to (u1,v1)
!           ! get r1 and z1
            CALL dmnf1 (u1, v1, r1, z1)
            phi = alp*v1
            xf = r1*COS(phi)
            yf = r1*SIN(phi)
            zf = z1
!           ! compute vector from surface 2 to surface 1
            xvec = xf - xp
            yvec = yf - yp
            zvec = zf - zp
!           ! compute the normal at (u1,v1) on the fixed surface
            isurf = 1
            CALL surnorm (isurf, u1, v1, xun, yun, zun)
!           ! compute the dot product of the normal and the vector
            DOtp = xvec*xun + yvec*yun + zvec*zun
            IF (dotp < 0._dp) THEN
!              WRITE(*,*) 'ip =',ip,'it = ',it
               d = -d
            ELSE IF (dotp == 0._dp) THEN
               distance = 0
               ierr = 1
               RETURN
            ENDIF
            distance = MIN(distance,d)
            ngood = ngood + 1
            dmat(ngood) = d
         END DO
      END DO
!
c     compute the mean and the standard deviation
      dmean = SUM(dmat(1:ngood)) / ngood
      std = SQRT(SUM((dmat(1:ngood) - dmean)**2)/ngood)
c     davg is the mean of the absolute values
      davg = SUM(ABS(dmat(1:ngood))) / ngood
      ngood2 = 0
      drms = 0
      DO i = 1, ngood
         IF(ABS(dmat(i) - dmean) .le. 2*std) THEN
            ngood2 = ngood2 + 1
            drms = drms + dmat(i)**2
         ENDIF
      ENDDO
      drms = SQRT(drms/ngood2)
!
      DEALLOCATE(dmat, xm1, xn1, rmn1, zmn1, xm2, xn2, rmn2, zmn2)

      END SUBROUTINE surfsep
