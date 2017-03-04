      MODULE angle_constraints
      USE vmec_main, ONLY: ns, mpol, ntor, dp, mpol1, lthreed, lasym
      USE vmec_params, ONLY: signgs, ntmax, rcc, rss, zsc, zcs, rsc, rcs, zss, zcc
      USE precon2d, ONLY: ictrl_prec2d
      IMPLICIT NONE
      INTEGER, PARAMETER :: pexp=4, m0=0, m1=1, m2=2, m3=3
      LOGICAL, PARAMETER :: lorigin=.FALSE.
      INTEGER            :: mrho, m, istat
      REAL, PARAMETER    :: p5=0.5_dp, zero=0
      REAL(dp), ALLOCATABLE :: t1m(:), t2m(:), cos_HB(:), sin_HB(:)
      REAL(dp), ALLOCATABLE :: rz_array0(:,:,:), xtempa(:)
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: arhod, arhom, brhod, brhom,   &
         ard2, arm2, azd2, azm2, arhod2, arhom2,                             &
         brd2, brm2, bzd2, bzm2, brhod2, brhom2
      REAL(dp), DIMENSION(:), ALLOCATABLE   :: crhod, sin2u, cos2u, sfact
      REAL(dp) :: sqp5
#ifdef _HBANGLE
      CONTAINS

!     CALLED FROM FIXARAY
      SUBROUTINE init_multipliers
      IMPLICIT NONE
      REAL(dp) :: dnorm, t0
!     NOTE: rho = SUM(m) rhomn*max(m,1)**pexp cos(mu-nv), etc
!           the extra m**pexp factor is to give good scaling

      IF (ALLOCATED(t1m)) RETURN

      sqp5 = SQRT(p5)
      mrho = mpol1      
!HB Angle constraint      
!      mrho = mpol1-1

      ALLOCATE(t1m(0:mpol), t2m(0:mpol), cos_HB(0:mpol1), sin_HB(0:mpol1), stat=istat)
      IF (istat .NE. 0) STOP 'Allocation error in init_multipliers!'

!
!     REF: Hirshman/Breslau paper
!
      t1m(m0) = 0; t2m(m0) = 0
      DO m = 1, mpol1+1
         t0 = MAX(1,m-1)
         t1m(m) = t0/m
         t0 = m+1
         t2m(m) = t0/m
         t1m(m) = t1m(m)**pexp
         t2m(m) = t2m(m)**pexp
      END DO

!      t1m(3) = 0
!SPH-TEST (IMPROVES CONDITION #?)
!      t1m(m2) = t2m(m2)

!HB Constraint  
!      t2m(mpol1) = 0; t2m(mpol) = 0

!      dnorm = MAXVAL(t1m); dnorm = MAX(dnorm,MAXVAL(t2m))
!      t1m = t1m/dnorm; t2m = t2m/dnorm

      cos_HB(m0) = 1;   sin_HB(m0) = 0

      DO m = 1, mrho   
         dnorm = SQRT(t1m(m+1)**2 + t2m(m-1)**2)
         IF (dnorm .EQ. zero) CYCLE
         cos_HB(m) = t1m(m+1)/dnorm
         sin_HB(m) = t2m(m-1)/dnorm
      END DO

      IF (mrho .NE. mpol1) THEN
         cos_HB(mpol1) = 0;  sin_HB(mpol1) = 1
      END IF

!DOESN'T CONVERGE IF THIS GOES BEFORE cos,sin_HB calculations (F_0 too larger???)  
!      t2m(m0) = t1m(m2)

      END SUBROUTINE init_multipliers

      SUBROUTINE free_multipliers
      IF (ALLOCATED(t1m)) DEALLOCATE (t1m, t2m, cos_HB, sin_HB)
      IF (ALLOCATED(rz_array0)) DEALLOCATE (rz_array0)
      IF (ALLOCATED(xtempa)) DEALLOCATE (xtempa)
      IF (ALLOCATED(arhod)) DEALLOCATE(arhod, arhom, brhod, brhom, crhod,      &
                                       ard2, arm2, azd2, azm2, arhod2, arhom2, &
                                       brd2, brm2, bzd2, bzm2, brhod2, brhom2)
      IF (ALLOCATED(cos2u)) DEALLOCATE(cos2u, sin2u)
      IF (ALLOCATED(sfact)) DEALLOCATE(sfact)

      END SUBROUTINE free_multipliers

      SUBROUTINE store_init_array(rzl_array)
      USE vmec_main, ONLY: neqs2, nznt, nzeta, cosmui, sinmui, cosmu, sinmu,   &
                           nzeta, ntheta2, ntheta3, hs
      REAL(dp), DIMENSION(ns*(1+ntor),0:mpol1,3*ntmax), INTENT(inout) :: rzl_array
      INTEGER :: nsp1, l, js
      
      IF (ALLOCATED(rz_array0)) DEALLOCATE (rz_array0)
      ALLOCATE(rz_array0(ns*(1+ntor),0:mpol1,2*ntmax))
      rz_array0 = rzl_array(:,:,1:2*ntmax)
      rzl_array(:,:,1:2*ntmax) = 0

      CALL init_multipliers
      CALL get_rep_mismatch(rz_array0, rzl_array)

      IF (ALLOCATED(xtempa)) DEALLOCATE (xtempa)
      ALLOCATE (xtempa(neqs2))            !Used as temp storage of xc in funct3d

      nsp1=ns+1
      IF (ALLOCATED(arhod)) DEALLOCATE(arhod, arhom, brhod, brhom, crhod,      &
                                       ard2, arm2, azd2, azm2, arhod2, arhom2, &
                                       brd2, brm2, bzd2, bzm2, brhod2, brhom2)
      ALLOCATE(arhod(nsp1,2),arhom(nsp1,2),brhod(nsp1,2),brhom(nsp1,2),   &
               ard2(nsp1,2), arm2(nsp1,2), arhod2(nsp1,2), arhom2(nsp1,2),&
               azd2(nsp1,2), azm2(nsp1,2),                                &
               brd2(nsp1,2), brm2(nsp1,2), brhod2(nsp1,2), brhom2(nsp1,2),&
               bzd2(nsp1,2), bzm2(nsp1,2), crhod(nsp1))
      IF (ALLOCATED(cos2u)) DEALLOCATE(cos2u, sin2u)
      ALLOCATE (cos2u(nznt), sin2u(nznt))
      DO l=1,nzeta
         cos2u(l:nznt:nzeta) = cosmui(:,m2)
         sin2u(l:nznt:nzeta) = sinmui(:,m2)
      END DO

      IF (ALLOCATED(sfact)) DEALLOCATE(sfact)
      ALLOCATE(sfact(ns))
      DO js=1,ns
         sfact(js) = hs*(js-1)
      END DO

      END SUBROUTINE store_init_array

      SUBROUTINE getrz (rz_array)
!      USE vmec_main, ONLY: iter2
      USE xstuff, ONLY: xc, xstore
      IMPLICIT NONE
      REAL(dp), DIMENSION(ns*(1+ntor),0:mpol1,2*ntmax), INTENT(inout) :: rz_array
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: rhocc, rhoss, rhocs, rhosc
      REAL(dp), DIMENSION(:), ALLOCATABLE   :: r0c, r0s, z0c, z0s
      INTEGER :: nsnt, mrho1, istat

!     INPUT:  rho quasi-polar components (Fourier modes) stored in R 
!             R,Z(m=0) centroid components stored in Z(m=0)
!     OUTPUT: Cylindrical component Fourier modes of R, Z
!     REF: HB Paper, Eq 32
!
!     FOR STELLARATOR SYMMETRY, rho = cos(mu-nv), R' = sum(m^pexp Rmn cos(mu-nv)), Z' = sum(m^pexp Zmn sin(mu-nv)
!     R' = (rhocc + rhoss) cosu
!     Z' = (rhocc + rhoss) sinu
!
!     FOR ASYMMETRY, ADD rho = sin(mu-nv) TERMS
!     R' = (rhosc + rhocs) cosu
!     Z' = (rhosc + rhocs) sinu
!

      nsnt = ns*(1+ntor)
      mrho1 = mrho+1
!      IF (mrho .NE. mpol1) STOP 'mrho != mpol1'

!     Enforce asymptotic behavior near axis (HB paper after Eq (34))
      IF (lorigin) THEN
         rz_array(2:nsnt:ns,m2,1:ntmax) = sqp5*rz_array(3:nsnt:ns,m2,1:ntmax)
         IF (mpol1 .GE. m3) THEN
            rz_array(2:nsnt:ns,m1:m3:2,1:ntmax) = p5*rz_array(3:nsnt:ns,m1:m3:2,1:ntmax)
            rz_array(2:nsnt:ns,m3+1:,1:ntmax) = 0
         END IF
      END IF

      ALLOCATE (rhocc(nsnt,0:mpol1+1), rhoss(nsnt,0:mpol1+1),               &
                r0c(nsnt), z0s(nsnt),stat=istat)
      IF (istat .NE. 0) STOP 'Allocation Error #1 in GETRZ'
      rhocc(:,0:mrho) = rz_array(:,0:mrho,rcc)
      r0c = rz_array(:,m0,zsc+ntmax)     !Used for temp storage
      rhocc(:,mrho1) = 0; rhocc(:,mpol1+1) = 0

      IF (lthreed) THEN
         rhoss(:,1:mrho) = rz_array(:,1:mrho,rss)
         rhoss(:,m0) = 0
         rhoss(:,mrho1) = 0; rhoss(:,0) = 0
         z0s = rz_array(:,m0,zcs+ntmax)
      END IF

      MODES: DO m = 0, mpol1
      IF (m .EQ. m0) THEN
         rz_array(:,m0,rcc) = r0c + t2m(m0)*rhocc(:,m1)
         rz_array(:,m0,zsc+ntmax) = 0
         IF (lthreed) THEN
            rz_array(:,m0,zcs+ntmax) = z0s - t2m(m0)*rhoss(:,m1)*signgs
            rz_array(:,m0,rss) = 0
         END IF
      ELSE 
         rz_array(:,m,rcc)       = (t1m(m)*rhocc(:,m-1) + t2m(m)*rhocc(:,m+1))
         rz_array(:,m,zsc+ntmax) =-(t1m(m)*rhocc(:,m-1) - t2m(m)*rhocc(:,m+1))*signgs
         IF (lthreed) THEN
            rz_array(:,m,rss)       = (t1m(m)*rhoss(:,m-1) + t2m(m)*rhoss(:,m+1))
            rz_array(:,m,zcs+ntmax) = (t1m(m)*rhoss(:,m-1) - t2m(m)*rhoss(:,m+1))*signgs
         END IF
      ENDIF
      END DO MODES

      DEALLOCATE (rhocc, rhoss, r0c, z0s)

      IF (.NOT.lasym) GOTO 1002

      ALLOCATE (rhosc(nsnt,0:mrho1), rhocs(nsnt,0:mrho1),               &
                r0s(nsnt), z0c(nsnt), stat=istat)
      IF (istat .NE. 0) STOP 'Allocation Error #2 in GETRZ'

      rhosc(:,0:mrho) = rz_array(:,0:mrho,rsc)
      rhosc(:,mrho1) = 0  
      z0c = rz_array(:,m0,zcc+ntmax)
      IF (lthreed) THEN
         rhocs(:,0:mrho) = rz_array(:,0:mrho,rcs)
         rhocs(:,mrho1) = 0; rhocs(:,0) = 0
         r0s = rz_array(:,m0,zss+ntmax)
      END IF

      MODEA: DO m = 0, mpol1
      IF (m .EQ. 0) THEN
         rz_array(:,m0,zcc+ntmax) = z0c - t2m(m0)*rhosc(:,m1)*signgs  
         rz_array(:,m0,rsc) = 0
         IF (lthreed) THEN
            rz_array(:,m0,rcs) = r0s + t2m(m0)*rhocs(:,m1)
            rz_array(:,m0,zss+ntmax) = 0
         END IF
      ELSE 
         rz_array(:,m,rsc)       = (t1m(m)*rhosc(:,m-1) + t2m(m)*rhosc(:,m+1))
         rz_array(:,m,zcc+ntmax) = (t1m(m)*rhosc(:,m-1) - t2m(m)*rhosc(:,m+1))*signgs
         IF (lthreed) THEN
            rz_array(:,m,rcs)       = (t1m(m)*rhocs(:,m-1) + t2m(m)*rhocs(:,m+1))
            rz_array(:,m,zss+ntmax) =-(t1m(m)*rhocs(:,m-1) - t2m(m)*rhocs(:,m+1))*signgs
         END IF
      ENDIF
      END DO MODEA

      DEALLOCATE (rhosc, rhocs, r0s, z0c)

 1002 CONTINUE

      !Add initial (scaled) boundary values if needed
      IF (ictrl_prec2d .NE. 3) rz_array = rz_array + rz_array0 

      END SUBROUTINE getrz

      SUBROUTINE getfrho (frho_array)
      USE vmec_main, ONLY: iter2, fnorm, hs
      REAL(dp), DIMENSION(ns*(1+ntor),0:mpol1,2*ntmax), INTENT(inout) ::   &
                frho_array
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: frcc, frss, fzcs, fzsc,     &
                                               frsc, frcs, fzcc, fzss
      INTEGER  :: nsnt, mrho1, istat

!     INPUT:  Frho contains cylindrical components of FR, FZ
!     OUTPUT: Frho(1:ntmax) contains Quasi-polar components the MHD forces
!             Frho(ntmax+1:2*ntmax) contains centroid (m=0) forces
      nsnt = ns*(1+ntor)
      mrho1 = mrho+1

      ALLOCATE (frcc(nsnt,0:mrho1), frss(nsnt,0:mrho1),                 &
                fzsc(nsnt,0:mrho1), fzcs(nsnt,0:mrho1), stat=istat)
      IF (istat .NE. 0) STOP 'Allocation error #1 in getfrho!'
      
      frcc(:,0:mrho) = frho_array(:,0:mrho,rcc);   frcc(:,mrho1)=0
      frho_array(:,:,rcc)=0

      fzsc(:,0:mrho) = frho_array(:,0:mrho,zsc+ntmax);  fzsc(:,mrho1) = 0
      frho_array(:,:,zsc+ntmax)=0

      IF (lthreed) THEN
         frss(:,0:mrho) = frho_array(:,0:mrho,rss); frss(:,mrho1) = 0
         frho_array(:,:,rss)=0
         fzcs(:,0:mrho) = frho_array(:,0:mrho,zcs+ntmax); fzcs(:,mrho1) = 0
         frho_array(:,:,zcs+ntmax)=0
      END IF

      MODES: DO m = 0, mrho
      IF (m .EQ. m0) THEN
         frho_array(:,m0,zsc+ntmax) = frcc(:,m0)        !storage
         frho_array(:,m0,rcc)       = cos_HB(m0)*(frcc(:,m1) - signgs*fzsc(:,m1))
         IF (.NOT.lthreed) CYCLE
         frho_array(:,m0,zcs+ntmax)  = fzcs(:,m0)
      ELSE 
         frho_array(:,m,rcc)       =  cos_HB(m)*(frcc(:,m+1) - signgs*fzsc(:,m+1))    &
                                   +  sin_HB(m)*(frcc(:,m-1) + signgs*fzsc(:,m-1))
         IF (.NOT.lthreed) CYCLE
         frho_array(:,m,rss)       =  cos_HB(m)*(frss(:,m+1) + signgs*fzcs(:,m+1))    &
                                   +  sin_HB(m)*(frss(:,m-1) - signgs*fzcs(:,m-1))
      ENDIF
      END DO MODES

      DEALLOCATE (frcc, frss, fzsc, fzcs)
      IF (lthreed) THEN
         IF (ANY(frho_array(:,m0,rss) .NE. zero)) STOP 'FRHO(m0,rss) != 0'
      END IF

      IF (.NOT. lasym) GOTO 1000

      ALLOCATE (frsc(nsnt,0:mrho1), frcs(nsnt,0:mrho1),  &
                fzcc(nsnt,0:mrho1), fzss(nsnt,0:mrho1),  stat=istat)
      IF (istat .NE. 0) STOP 'Allocation error #2 in getfrho!'

      frsc(:,0:mrho) = frho_array(:,0:mrho,rsc);  frsc(:,mrho1) = 0
      frho_array(:,:,rsc) = 0
      fzcc(:,0:mrho) = frho_array(:,0:mrho,zcc+ntmax);  fzcc(:,mrho1) = 0
      frho_array(:,:,zcc+ntmax) = 0
      IF (lthreed) THEN
         frcs(:,0:mrho) = frho_array(:,0:mrho,rcs); frcs(:,mrho1) = 0
         frho_array(:,:,rcs)=0
         fzss(:,0:mrho) = frho_array(:,0:mrho,zss+ntmax); fzss(:,mrho1) = 0 
         frho_array(:,:,zss+ntmax)=0
      END IF

      MODEA: DO m = 0, mrho
      IF (m .EQ. m0) THEN
         frho_array(:,m0,zcc+ntmax) = fzcc(:,m0)
         IF (.NOT.lthreed) CYCLE
         frho_array(:,m0,zss+ntmax) = frcs(:,m0)
         frho_array(:,m0,rcs)       =cos_HB(m0)*(frcs(:,m1) - signgs*fzss(:,m1))
      ELSE 
         frho_array(:,m,rsc)       = cos_HB(m)*(frsc(:,m+1) + signgs*fzcc(:,m+1))    &
                                   + sin_HB(m)*(frsc(:,m-1) - signgs*fzcc(:,m-1))
         IF (.NOT.lthreed) CYCLE
         frho_array(:,m,rcs)       = cos_HB(m)*(frcs(:,m+1) - signgs*fzss(:,m+1))    &
                                   + sin_HB(m)*(frcs(:,m-1) + signgs*fzss(:,m-1))
      ENDIF
      END DO MODEA

      DEALLOCATE (frsc, frcs, fzcc, fzss)

      IF (lasym) THEN
         IF (ANY(frho_array(:,m0,rsc) .NE. zero)) STOP 'FRHO(m=0,rsc) != 0'
      END IF

 1000 CONTINUE

!     Use asymptotic behavior near axis, rather than evolution (SPH010214)
      IF (lorigin) THEN
         frho_array(2:nsnt:ns,m1:,1:ntmax) = 0
      END IF

!ADD UNIQUE POLAR AXIS "CONSTRAINT"
      frho_array(:,m1,1:ntmax) = 0                   !Unique angle


      END SUBROUTINE getfrho


      SUBROUTINE scalfor_rho(gcr, gcz)
      USE vmec_main
      USE vmec_params
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) ::    &
                gcr, gcz
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER :: edge_pedestal=0.05_dp
      INTEGER :: m , mp, n, js, jmax, jmin4(0:mnsize-1)
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: arho, brho, drho,      &
                                                 arho2,brho2,drho2
      REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: acen, bcen, dcen
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: gcen
      REAL(dp) :: tar, taz, tc, tc1
!-----------------------------------------------
      jmax = ns
      IF (ivac .lt. 1) jmax = ns1

      ALLOCATE (acen(ns,0:ntor), bcen(ns,0:ntor), dcen(ns,0:ntor),      &
                gcen(ns,0:ntor,4))
!
!     FIRST, SCALE m=0 R,Z COMPONENTS
!
      gcen=0
      gcen(:,:,1) = gcz(:,:,m0,zsc)        !used for r0n storage
      IF (lasym) gcen(:,:,3) = gcz(:,:,m0,zcc)
      IF (lthreed) THEN
         gcen(:,:,2) = gcz(:,:,m0,zcs)
         IF (lasym) gcen(:,:,4) = gcz(:,:,m0,zss)
      END IF
  
      IF (SUM(gcz*gcz) .NE. SUM(gcen*gcen)) STOP 'ERROR #1 IN SCALFOR_RHO'

      CALL scalaxis(arm, ard, crd, acen, bcen, dcen, jmax, gcen(:,:,1), 0)
      IF (lasym) CALL scalaxis(azm, azd, crd, acen, bcen, dcen, jmax, gcen(:,:,3), 0)
      IF (lthreed) THEN
         CALL scalaxis(azm, azd, crd, acen, bcen, dcen, jmax, gcen(:,:,2), 1)
         IF (lasym) CALL scalaxis(arm, ard, crd, acen, bcen, dcen, jmax, gcen(:,:,4), 1)
      END IF

      gcz = 0
!
!     RESTORE SCALED CENTROID VALUES
!
      gcz(:,:,m0,zsc) = gcen(:,:,1)                !R0n-cc
      IF (lasym) gcz(:,:,m0,zcc) = gcen(:,:,3)     !Z0n-cc
      IF (lthreed) THEN
         gcz(:,:,m0,zcs) = gcen(:,:,2)             !Z0n-cs
         IF (lasym) gcz(:,:,m0,zss) = gcen(:,:,4)  !R0n-cs
      END IF

!
!     NEXT, SCALE QUASI-POLAR COMPONENTS
!
      ALLOCATE (arho(ns,0:ntor,0:mpol1), brho(ns,0:ntor,0:mpol1),       &
                drho(ns,0:ntor,0:mpol1), arho2(ns,0:ntor,0:mpol1),      &
                brho2(ns,0:ntor,0:mpol1),drho2(ns,0:ntor,0:mpol1))

      DO m = 0, mpol1
         mp = MOD(m,2)+1
         DO n = 0, ntor
            DO js = 1, ns
               arho(js,n,m) = -(arhom(js+1,mp) + brhom(js+1,mp)*m**2)
               brho(js,n,m) = -(arhom(js,mp)   + brhom(js,mp)  *m**2)
               drho(js,n,m) = -(arhod(js,mp)   + brhod(js,mp)  *m**2    &
                            +   crhod(js)*(n*nfp)**2)
               arho2(js,n,m)= -p5*(arhom2(js+1,mp)+ brhom2(js+1,mp)*m**2)
               brho2(js,n,m)= -p5*(arhom2(js,mp)  + brhom2(js,mp)  *m**2)
               drho2(js,n,m)= -p5*(arhod2(js,mp)  + brhod2(js,mp)  *m**2)
            END DO
         END DO
      END DO

      CALL avg_rho(arho, arho2)
      CALL avg_rho(brho, brho2)
      CALL avg_rho(drho, drho2)

      IF (jmax .ge. ns) THEN
!
!     SMALL EDGE PEDESTAL NEEDED TO IMPROVE CONVERGENCE
!     IN PARTICULAR, NEEDED TO ACCOUNT FOR POTENTIAL ZERO
!     EIGENVALUE DUE TO NEUMANN (GRADIENT) CONDITION AT EDGE
!
         drho(ns,:,:) = (1+edge_pedestal)*drho(ns,:,:)
!         drho(ns,:,3:)  = (1+2*edge_pedestal)*drho(ns,:,3:)
      END IF

      jmin4 = 2
      CALL tridslv (arho,drho,brho,gcr(:,:,:,rcc),jmin4,jmax,           &
                    mnsize-1,ns,1)
      IF (lthreed) CALL tridslv (arho,drho,brho,gcr(:,:,:,rss),         &
                   jmin4,jmax,mnsize-1,ns,1)
      IF (lasym) THEN
      CALL tridslv (arho,drho,brho,gcr(:,:,:,rsc),jmin4,jmax,           &
                    mnsize-1,ns,1)
      IF (lthreed) CALL tridslv (arho,drho,brho,gcr(:,:,:,rcs),         &
                   jmin4,jmax,mnsize-1,ns,1)
      END IF

      DEALLOCATE (arho, brho, drho, arho2, brho2, drho2, gcen)

      END SUBROUTINE scalfor_rho

      SUBROUTINE scalaxis(axm, axd, cxd, acen, bcen, dcen, jmax, gcen, iflag)
      USE vmec_main
      USE vmec_params
      REAL(dp), PARAMETER :: edge_pedestal=0.15_dp, fac=0.25_dp
!      REAL(dp), PARAMETER :: edge_pedestal=0.05_dp, fac=0.25_dp
      REAL(dp), INTENT(IN), DIMENSION(ns+1,2) ::  axm, axd
      REAL(dp), INTENT(IN), DIMENSION(ns+1)   ::  cxd
      REAL(dp), DIMENSION(ns,0:ntor) ::  acen, bcen, dcen
      REAL(dp), INTENT(INOUT) :: gcen(ns,0:ntor)
      REAL(dp) :: mult_fac
      INTEGER, INTENT(IN)     :: jmax, iflag
      INTEGER :: mp, n, js, mnsize0, jmin4(0:mnsize-1)

      mp=1
      DO n = 0, ntor
         DO js = 1, jmax
            acen(js,n) = -axm(js+1,mp)
            bcen(js,n) = -axm(js,mp)
            dcen(js,n) = -(axd(js,mp)+cxd(js)*(n*nfp)**2)
         END DO
      END DO

      IF (jmax .ge. ns) THEN
!
!     SMALL EDGE PEDESTAL NEEDED TO IMPROVE CONVERGENCE
!     IN PARTICULAR, NEEDED TO ACCOUNT FOR POTENTIAL ZERO
!     EIGENVALUE DUE TO NEUMANN (GRADIENT) CONDITION AT EDGE
!
         dcen(ns,:)     = (1+edge_pedestal)*dcen(ns,:)
!
!     STABILIZATION ALGORITHM FOR ZC_00(NS)
!     FOR UNSTABLE CASE, HAVE TO FLIP SIGN OF -FAC -> +FAC FOR CONVERGENCE
!     COEFFICIENT OF < Ru (R Pvac)> ~ -fac*(z-zeq) WHERE fac (EIGENVALUE, OR
!     FIELD INDEX) DEPENDS ON THE EQUILIBRIUM MAGNETIC FIELD AND CURRENT,
!     AND zeq IS THE EQUILIBRIUM EDGE VALUE OF Z00
          IF (iflag .eq. 1) THEN
!
!     METHOD 1: SUBTRACT (INSTABILITY) Pedge ~ fac*z/hs FROM PRECONDITIONER AT EDGE
!
             mult_fac = MIN(fac, fac*hs*15)
             dcen(ns,0) = dcen(ns,0)*(1-mult_fac)/(1+edge_pedestal)
          END IF

      ENDIF


      mnsize0 = 1+ntor
      jmin4 = jmin3
      IF (iresidue.GE.0 .AND. iresidue.LT.3) jmin4(0)=2

      CALL tridslv(acen,dcen,bcen,gcen,jmin4,jmax,mnsize0-1,ns,1)

      END SUBROUTINE scalaxis

      SUBROUTINE avg_rho(ax, ax2)
      USE vmec_main
      IMPLICIT NONE
      REAL(dp), INTENT(inout), DIMENSION(ns*(1+ntor), 0:mpol1) :: ax, ax2
      REAL(dp), ALLOCATABLE   :: ax1(:,:)
      INTEGER :: m

      ALLOCATE (ax1(ns*(1+ntor),0:mpol1))

      ax1(:,m0) = cos_HB(m0)*t1m(m1)*(ax(:,m1)+ax2(:,m1))
      DO m=1,mpol1-1
         ax1(:,m) = cos_HB(m)*(t1m(m+1)*ax(:,m+1)+t2m(m-1)*ax2(:,m-1))   &
                  + sin_HB(m)*(t2m(m-1)*ax(:,m-1)+t1m(m+1)*ax2(:,m+1))
      END DO
!      ax1(:,m2) = ax1(:,m2)+sin_HB(m2)*t2m(m1)*ax2(:,m2-1)
      ax1(:,mpol1) = sin_HB(mpol1)*t2m(mpol1-1)*ax(:,mpol1-1)            &
                   + cos_HB(mpol1)*t2m(mpol1-1)*ax2(:,mpol1-1)
      ax = ax1

      DEALLOCATE (ax1)

      END SUBROUTINE avg_rho

      SUBROUTINE precondn_rho
      USE vmec_main, ONLY: ard, arm, brd, brm, azd, azm, bzd, bzm, crd
      USE fbal, ONLY: rzu_fac, rru_fac, frcc_fac, fzsc_fac
      REAL(dp), PARAMETER :: one=1
      
      arhod = ard + azd
      arhom = arm + azm
      brhod = brd + bzd
      brhom = brm + bzm
      crhod = crd + crd

      arhod2 = ard2-azd2; brhod2 = (brd2-bzd2)
      arhom2 = arm2-azm2; brhom2 = (brm2-bzm2)

      rzu_fac = (rzu_fac+rru_fac); rru_fac=rzu_fac
      frcc_fac(2:ns-1) = one/rzu_fac(2:ns-1)
      fzsc_fac(2:ns-1) = -frcc_fac(2:ns-1)

      END SUBROUTINE precondn_rho

      SUBROUTINE get_rep_mismatch(rz0_array, rho_array)
      USE vmec_main, ONLY: irst, hs
      IMPLICIT NONE
      REAL(dp), PARAMETER :: p5=0.5_dp, p25=p5*p5
      REAL(dp), DIMENSION(ns*(1+ntor),0:mpol1,2*ntmax), INTENT(in)  :: rz0_array
      REAL(dp), DIMENSION(ns*(1+ntor),0:mpol1,2*ntmax), INTENT(out) :: rho_array
      INTEGER  :: m, nsnt, mrho1, js, n, ntc
      REAL(dp) :: match, delta, t1(ns*(1+ntor)), t2(ns*(1+ntor)),    &
                  temp1(ns*(1+ntor)), temp2(ns*(1+ntor)), es
!
!     COMPUTES mismatch BETWEEN INITIAL R-Z REPRESENTATION AND QPOLAR FORM
!
      nsnt = SIZE(rz0_array,1)
      mrho1 = MAX(mrho+1,mpol1)
      temp1 = 0;  temp2 = 0

!
!     STORE m=0 AXIS data
!
      rho_array(:,m0,zsc+ntmax) = rz0_array(:,m0,rcc)
      IF (lthreed) THEN
         rho_array(:,m0,zcs+ntmax) = rz0_array(:,m0,zcs+ntmax)
      END IF

      DO m = 0, mrho
         IF (m .LT. mpol1-1) THEN
            t1 = (rz0_array(:,m+1,rcc) - signgs*rz0_array(:,m+1,zsc+ntmax))/t1m(m+1)
            IF (m .LE. 1) THEN
               t2 = t1
            ELSE
               t2 = (rz0_array(:,m-1,rcc) + signgs*rz0_array(:,m-1,zsc+ntmax))/t2m(m-1)
            END IF
         ELSE IF (m .GT. 1) THEN
            t2 = (rz0_array(:,m-1,rcc) + signgs*rz0_array(:,m-1,zsc+ntmax))/t2m(m-1)
            t1 = t2
         END IF
         rho_array(:,m,rcc)  = p25*(t1 + t2)
         temp1 = temp1 + (t1 - t2)**2
         temp2 = temp2 + (t1 + t2)**2

         IF (.NOT.lthreed) CYCLE

         IF (m .LT. mpol1-1) THEN
            t1 = (rz0_array(:,m+1,rss) + signgs*rz0_array(:,m+1,zcs+ntmax))/t1m(m+1)
            IF (m .LE. 1) THEN
               t2 = t1
            ELSE
               t2 = (rz0_array(:,m-1,rss) - signgs*rz0_array(:,m-1,zcs+ntmax))/t2m(m-1)
            END IF
         ELSE
            t2 = (rz0_array(:,m-1,rss) - rz0_array(:,m-1,rcs+ntmax))/t2m(m-1)
            t1 = t2
         END IF
         rho_array(:,m,rss)  = p25*(t1 + t2)
         temp1 = temp1 + (t1 - t2)**2
         temp2 = temp2 + (t1 + t2)**2
      END DO

      rho_array(:,m0,zsc+ntmax) = rho_array(:,m0,zsc+ntmax) - t2m(m0)*rho_array(:,m1,rcc)
      IF (lthreed) THEN
         rho_array(:,m0,zcs+ntmax) = rho_array(:,m0,zcs+ntmax) + t2m(m0)*rho_array(:,m1,rss)*signgs
      END IF

!
!     NON-SYMMETRIC CONTRIBUTIONS
!
      IF (.NOT.lasym) GOTO 1000

      rho_array(:,m0,zcc+ntmax) = rz0_array(:,m0,zcc+ntmax)
      IF (lthreed) THEN
         rho_array(:,m0,zss+ntmax) = rz0_array(:,m0,rcs)
      END IF

      DO m = 0, mrho
         IF (m .LT. mpol1-1) THEN
            t1 = (rz0_array(:,m+1,rsc) + signgs*rz0_array(:,m+1,zcc+ntmax))/t1m(m+1)
            IF (m .LE. 1) THEN
               t2 = t1
            ELSE
               t2 = (rz0_array(:,m-1,rsc) - signgs*rz0_array(:,m-1,zcc+ntmax))/t2m(m-1)
            END IF
         ELSE
            t2 = (rz0_array(:,m-1,rsc) - signgs*rz0_array(:,m-1,zcc+ntmax))/t2m(m-1)
            t1 = t2
         END IF
         rho_array(:,m,rsc)  = p25*(t1 + t2)
         temp1 = temp1 + (t1 - t2)**2
         temp2 = temp2 + (t1 + t2)**2
         
         IF (.not.lthreed) CYCLE

         IF (m .LT. mpol1-1) THEN
            t1 = (rz0_array(:,m+1,rcs) - signgs*rz0_array(:,m+1,zss+ntmax))/t1m(m+1)
            IF (m .LE. 1) THEN
               t2 = t1
            ELSE
               t2 = (rz0_array(:,m-1,rcs) + signgs*rz0_array(:,m-1,zss+ntmax))/t2m(m-1)
            END IF
         ELSE
            t2 = (rz0_array(:,m-1,rcs) + signgs*rz0_array(:,m-1,zss+ntmax))/t2m(m-1)
            t1 = t2
         END IF
         rho_array(:,m,rcs)  = p25*(t1 + t2)
         temp1 = temp1 + (t1 - t2)**2
         temp2 = temp2 + (t1 + t2)**2
      END DO

      rho_array(:,m0,zcc+ntmax) = rho_array(:,m0,zcc+ntmax) + t2m(m0)*rho_array(:,m1,rsc)
      IF (lthreed) THEN
         rho_array(:,m0,zss+ntmax) = rho_array(:,m0,zss+ntmax) - t2m(m0)*rho_array(:,m1,rcs)*signgs
      END IF

 1000 CONTINUE

      IF (irst .EQ. 1) GOTO 2000

      DO js = 2, ns-1
         es = (js-1)*hs
         DO n = 0, ntor
            ntc = js+ns*ntor
            DO m = 0, mrho
               IF (MOD(m,2) .EQ. 0) THEN
                  rho_array(ntc, m, 1:ntmax) = SQRT(es)*rho_array(ns*(1+ntor), m, 1:ntmax)
               ELSE
                  rho_array(ntc, m, 1:ntmax) = es*rho_array(ns*(1+ntor), m, 1:ntmax)
               END IF
            END DO
         END DO
      END DO

 2000 CONTINUE

      WRITE (*, '(a)') '  Quasi-polar representation mismatch vs radius'
      DO js=2,ns,ns-2
         delta = SUM(temp1(js:nsnt:ns))
         match = SUM(temp2(js:nsnt:ns))
         IF (match .NE. 0._dp) WRITE (*, '(1x,i4,1p,e10.2)')js,delta/match
      END DO

!      rho_array=0  !Use this if delta at the edge is not good
      rz_array0=0   !Initial jacobian can change sign: check

!      IF (ALLOCATED(rho1)) DEALLOCATE(rho1)
!      ALLOCATE (rho1(nsnt,ntmax))
!      rho1 = rho_array(:,m1,1:ntmax)

      END SUBROUTINE get_rep_mismatch
#endif
      END MODULE angle_constraints
