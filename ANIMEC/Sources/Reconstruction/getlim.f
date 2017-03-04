      SUBROUTINE getlim
      USE vmec_main
      USE realspace
      USE vsvd
      USE mgrid_mod, ONLY: nlim, limitr, rlim, zlim, reslim, seplim
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: resup = 3.0_dp
      REAL(rprec), PARAMETER :: resdn = 0.5_dp
      REAL(rprec), PARAMETER :: eps = 0.005_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ntheta_2pi, nphi_plane, i,
     1   nthtot, iexpand, ishrink, ionlim, n,
     2   limpts, nonlim, nexpand, nshrink, ilim0, nlim0
      REAL(rprec), DIMENSION(2*ntheta1) ::
     1   rbdy, zbdy, rubdy, zubdy
      REAL(rprec) :: fshrink, distmax, fexpand
C-----------------------------------------------

c
c     DETERMINES WHEN PLASMA TOUCHES LIMITER
c     USE DOUBLE THE NO. OF THETA POINTS FOR INCREASED RESOLUTION
c
      ntheta_2pi = ntheta1
      nphi_plane = 1                 !Pick a phi plane (phi = 0 for now)

      rbdy(:ntheta3*2-1:2) = r1(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     1   nphi_plane,0) + r1(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     2   nphi_plane,1)
      zbdy(:ntheta3*2-1:2) = z1(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     1   nphi_plane,0) + z1(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     2   nphi_plane,1)
      rubdy(:ntheta3*2-1:2) = ru(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     1   nphi_plane,0) + ru(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     2   nphi_plane,1)
      zubdy(:ntheta3*2-1:2) = zu(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     1   nphi_plane,0) + zu(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     2   nphi_plane,1)


      IF (.not.lasym) THEN
!     FOR NOW, THIS ONLY WORKS FOR NZETA=1 (PHI=0 PLANE)
!     TO EXTEND TO OTHER PHI PLANES, MUST USE IREFLECT(NZETA)
      DO i = 1, ntheta_2pi - ntheta2
         rbdy(2*(ntheta2+i)-1) = rbdy(2*(ntheta1-ntheta2-i)+3)
      END DO
      DO i = 1, ntheta_2pi - ntheta2
         zbdy(2*(ntheta2+i)-1) = -zbdy(2*(ntheta1-ntheta2-i)+3)
      END DO
      DO i = 1, ntheta_2pi - ntheta2
         rubdy(2*(ntheta2+i)-1) = -rubdy(2*(ntheta1-ntheta2-i)+3)
      END DO
      DO i = 1, ntheta_2pi - ntheta2
         zubdy(2*(ntheta2+i)-1) = zubdy(2*(ntheta1-ntheta2-i)+3)
      END DO
      END IF

!
!     FIND EVEN INDEXED POINTS BY INTERPOLATION
!
      nthtot = 2*ntheta_2pi
      rbdy(nthtot) = .5_dp*(rbdy(1)+rbdy(nthtot-1))
      zbdy(nthtot) = .5_dp*(zbdy(1)+zbdy(nthtot-1))
      rubdy(nthtot) = .5_dp*(rubdy(1)+rubdy(nthtot-1))
      zubdy(nthtot) = .5_dp*(zubdy(1)+zubdy(nthtot-1))
      rbdy(2:(ntheta_2pi-1)*2:2) = .5_dp*(rbdy(3:ntheta_2pi*2-1:2) +
     1      rbdy(:ntheta_2pi*2-3:2))
      zbdy(2:(ntheta_2pi-1)*2:2) = .5_dp*(zbdy(3:ntheta_2pi*2-1:2) +
     1      zbdy(:ntheta_2pi*2-3:2))
      rubdy(2:(ntheta_2pi-1)*2:2) = .5_dp*(rubdy(3:ntheta_2pi*2-1:2)+
     1   rubdy(:ntheta_2pi*2-3:2))
      zubdy(2:(ntheta_2pi-1)*2:2) = .5_dp*(zubdy(3:ntheta_2pi*2-1:2)+
     1   zubdy(:ntheta_2pi*2-3:2))

      fshrink = 0.0_dp
      distmax = SUM(rbdy(:nthtot)**2) + SUM(zbdy(:nthtot)**2)
      fexpand = distmax
      iexpand = 0
      ishrink = 0
      ionlim = 0

      DO n = 1, nlim
         limpts = limitr(n)
         CALL cauchy (rbdy, zbdy, rubdy, zubdy, rlim(:,n), zlim(:,n),
     1      reslim(:,n), seplim(:,n), distmax, nthtot, limpts)

           DO i = 1,limpts
c       LIMITER POINT ON PLASMA
            IF( (ABS(reslim(i,n)-resdn).lt.eps) )then
!    .gt.      .or. (ABS(reslim(i,n)).gt.resup) )then
              ionlim = i
              nonlim = n
c       LIMITER POINT OUTSIDE PLASMA
            ELSE IF( reslim(i,n).lt.RESDN )then
              IF( seplim(i,n).le.fexpand )then
                fexpand = seplim(i,n)
                iexpand = i
                nexpand = n
              ENDIF
c       LIMITER POINT INSIDE PLASMA
            ELSE IF( reslim(i,n).ge.RESDN )then
              IF( seplim(i,n).gt.fshrink )then
                fshrink =  seplim(i,n)
                ishrink = i
                nshrink = n
              ENDIF
            ENDIF
          ENDDO
      END DO

c
c       LOGIC: IF THERE IS A LIMITER POINT INSIDE PLASMA, THEN MUST
c       SHRINK CONTOUR. OTHERWISE, IF THERE IS AT LEAST ONE LIMITER
c       POINT ON CONTOUR, AND ALL THE REST OUTSIDE, DO NOT CHANGE ANYTHING.
c       FINALLY, IF ALL LIMITER POINTS ARE OUTSIDE PLASMA, EXPAND PLASMA
c       TO OSCULATE WITH LIMITER
c
      IF (ishrink .gt. 0) THEN
         gphifac = -SQRT(fshrink)
         ilim0 = ishrink
         nlim0 = nshrink
      ELSE IF (ionlim .gt. 0) THEN
         gphifac = 0.
         ilim0 = ionlim
         nlim0 = nonlim
      ELSE IF (iexpand .gt. 0) THEN
         gphifac = SQRT(fexpand)
         ilim0 = iexpand
         nlim0 = nexpand
      ENDIF

      dlim_min = gphifac
      rlim_min = rlim(ilim0,nlim0)
      zlim_min = zlim(ilim0,nlim0)

c     overALL damping in time, rsfac/SQRT(rsfac) = SQRT(rsfac)
      gphifac = gphifac/r01(ns)
      IF (ABS(gphifac) .gt. 0.04*one)
     1   gphifac = 0.04*gphifac/ABS(gphifac)
      gphifac = 0.20*gphifac/SQRT(rsfac)

      END SUBROUTINE getlim
