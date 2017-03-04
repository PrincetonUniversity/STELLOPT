      SUBROUTINE modbooz(bmod, nfp, icontrol, nsval)
      USE boozer_params
      USE read_boozer_mod, ONLY: bmn_b=>bmnc_b, ixm_b, ixn_b
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nfp, icontrol, nsval
      REAL(rprec), DIMENSION(nunv), INTENT(out) :: bmod
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mn, m, n, i, lk, lz, lt, isgn
      REAL(rprec), PARAMETER :: zero = 0, one = 1
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: cost, sint
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: thgrd, ztgrd
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE, SAVE ::
     1     cosmm, sinmm, cosnn, sinnn
      REAL(rprec) :: twopi, dth, dzt
C-----------------------------------------------
!
!     PERFORMS INVERSE TRANSFORM OF |BMN_B| INTO BOOZER ANGLE SPACE
!     CAN BE EASILY EXPANDED TO INVERSE TRANSFORM OTHER QUANTITIES (R, Z, p)
!
      twopi = 8*ATAN(one)

      IF (.not. ALLOCATED(cosmm)) THEN
         ALLOCATE (cosmm(nunv,0:mboz), sinmm(nunv,0:mboz),
     1             cosnn(nunv,0:nboz), sinnn(nunv,0:nboz),
     2             thgrd(nunv), ztgrd(nunv), stat = i)
      IF (i .ne. 0) STOP 'allocation error in modbooz'

!
!     COMPUTE POLOIDAL (thgrd) AND TOROIDAL (ztgrd) ANGLES
!
!        dth = twopi/REAL(nu_boz,rprec)         !USE THIS FOR FULL 2-pi
         dth = twopi/REAL(2*(nu2_b-1),rprec)    !Half-around in theta
         dzt = twopi/REAL(nv_boz,rprec)
         lk = 0

         DO lt = 1, nu2_b
            DO lz = 1, nv_boz
               lk = lk + 1
               thgrd(lk) = (lt-1)*dth
               ztgrd(lk) = (lz-1)*dzt
            END DO
         END DO

         IF (lk .ne. nunv) STOP 'lk != nunv in modbooz'

         cosmm(:,0) = one
         sinmm(:,0) = zero
         cosmm(:,1) = COS(thgrd)
         sinmm(:,1) = SIN(thgrd)

         cosnn(:,0) = one
         sinnn(:,0) = zero
         cosnn(:,1) = COS(ztgrd)
         sinnn(:,1) = SIN(ztgrd)

         DEALLOCATE (thgrd, ztgrd)

         DO m = 2,mboz
            cosmm(:,m) = cosmm(:,m-1)*cosmm(:,1)
     1                 - sinmm(:,m-1)*sinmm(:,1)
            sinmm(:,m) = sinmm(:,m-1)*cosmm(:,1)
     1                 + cosmm(:,m-1)*sinmm(:,1)
         END DO

         DO n = 2,nboz
            cosnn(:,n) = cosnn(:,n-1)*cosnn(:,1)
     1                 - sinnn(:,n-1)*sinnn(:,1)
            sinnn(:,n) = sinnn(:,n-1)*cosnn(:,1)
     1                 + cosnn(:,n-1)*sinnn(:,1)
         END DO
      END IF

      ALLOCATE (cost(nunv), sint(nunv), stat = i)

      bmod = zero

      DO mn = 1,mnboz
        m = ixm_b(mn)
        n = ABS(ixn_b(mn))/nfp
        IF (m .gt. mboz) STOP 'm > mboz in modbooz'
        IF (n .gt. nboz) STOP 'n > nboz in modbooz'
        isgn = sign(1,ixn_b(mn))
        cost = cosmm(:,m)*cosnn(:,n)
     1       + sinmm(:,m)*sinnn(:,n)*isgn
c       sint = sinmm(:,m)*cosnn(:,n)
c    1       - cosmm(:,m)*sinnn(:,n)*isgn
        bmod = bmod + bmn_b(mn, nsval)*cost
c       r12  = r12  + rmnb(mn)*cost
c       z12  = z12  + zmnb(mn)*sint
c       p12  = p12  + pmnb(mn)*sint
      END DO

      DEALLOCATE (cost, sint)
      IF (icontrol .eq. -1) DEALLOCATE
     1   (cosmm, sinmm, cosnn, sinnn)

      END SUBROUTINE modbooz
