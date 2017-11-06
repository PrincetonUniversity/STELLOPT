#if defined(SKS)
      SUBROUTINE residue_par (gcr, gcz, gcl)
      USE vmec_main, p5 => cp5
      USE vmec_params, ONLY: rss, zcs, rsc, zcc,                        &
                             meven, modd, ntmax, signgs
      USE realspace, ONLY: phip
      USE vsvd
      USE xstuff
      USE precon2d
      USE parallel_include_module
      USE parallel_vmec_module, ONLY: tlglob_arr, trglob_arr,           &
                                      lactive, SAXLASTNTYPE
      USE blocktridiagonalsolver, ONLY: L_COLSCALE
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,0:mpol1,ns,ntmax), INTENT(INOUT) :: &
        gcr, gcz, gcl
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: n0=0, m0=0, m1=1
      INTEGER, PARAMETER :: n3d=0, nasym=1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nsfix, jedge, delIter
      REAL(dp) :: r1, tnorm, fac, tmp, tmp2(ns), ftotal

      INTEGER  :: i, j, k, l, m, blksize, left, right
      INTEGER, ALLOCATABLE, DIMENSION(:) :: counts, disps
      INTEGER :: MPI_STAT(MPI_STATUS_SIZE)
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: send_buf
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: recv_buf
      REAL(dp) :: tredon, tredoff
!-----------------------------------------------
      CALL second0 (treson)

#ifdef _HBANGLE
!FREE-BDY RFP MAY NEED THIS TO IMPROVE CONVERGENCE (SPH 022514)
      IF (lfreeb .AND. lrfp) THEN
         fac = 0
         IF (ictrl_prec2d .EQ. 0) fac = 1.E-1_dp
         gcz(ns,0,m0,:) = fac*gcz(ns,0,m0,:)
      END IF
#else
!
!     SYMMETRIC PERTURBATIONS (BASED ON POLAR RELATIONS):
!        Rss(n) = Zcs(n), n != 0
!     ASYMMETRIC PERTURBATIONS:
!        Rsc(n) = Zcc(n), ALL n
!
!     INTERNALLY:
!        XC(rss) = .5*(Rss + Zcs), XC(zcs) = .5*(Rss - Zcs) -> 0
!        XC(rsc) = .5*(Rsc + Zcc), XC(zcc) = .5*(Rsc - Zcc) -> 0
!     THIS IMPLIES THE CONSTRAINT
!        3D ONLY : GC(zcs) = 0;  
!        ASYM:     GC(zcc) = 0
!

      IF (lthreed) CALL constrain_m1_par(gcr(:,m1,:,rss), gcz(:,m1,:,zcs))
      IF (lasym)   CALL constrain_m1_par(gcr(:,m1,:,rsc), gcz(:,m1,:,zcc))

      !print *, ">>>>>",lfreeb.AND.lrfp
      IF (lfreeb .AND. lrfp) THEN
         fac = 0
         IF (ictrl_prec2d .EQ. 0) fac = 1.E-1_dp
         gcr(0,m0,ns,:) = fac*gcr(0,m0,ns,:)
         gcz(0,m0,ns,:) = fac*gcz(0,m0,ns,:)
      END IF
#endif

!     PRECONDITIONER MUST BE CALCULATED USING RAW (UNPRECONDITIONED) FORCES
      IF (ictrl_prec2d.GE.2 .OR. ictrl_prec2d.EQ.-1) RETURN

!
!     PUT FORCES INTO PHIFSAVE UNITS USED BY PRECONDITIONERS, FNORM
!
      IF (phifac .eq. zero) THEN
         STOP 'phifac = 0 in residue'
      ELSE
         tnorm = phifsave/phifac           !put all forces into phifac=phifsave units
      END IF

      IF (lrecon) THEN

          CALL second0(tredon)
          tmp = SUM(gcr(n0,m0,tlglob:trglob,1))
          CALL MPI_Allreduce(tmp,r1,1,MPI_REAL8,MPI_SUM,NS_COMM,MPI_ERR)
          CALL second0(tredoff)
          allreduce_time = allreduce_time + (tredoff - tredon)
        fsqsum0 = signgs*hs*r1/r0scale
        nsfix = 1                   !fix origin for reconstruction mode
        gcr(:,:,tlglob:trglob,:) = gcr(:,:,tlglob:trglob,:) * tnorm**2
        gcz(:,:,tlglob:trglob,:) = gcz(:,:,tlglob:trglob,:) * tnorm**2
        gcl(:,:,tlglob:trglob,:) = gcl(:,:,tlglob:trglob,:) * tnorm
        IF (iopt_raxis.GT.0 .AND. iresidue.EQ.2                        &
           .AND. fsq.LT.fopt_axis) iresidue = 3
        IF (iresidue .LT. 3) gcr(n0,m0,nsfix,1) = zero
      ELSE
!
!     ADJUST PHIEDGE
!
         IF (imovephi .GT. 0) CALL movephi1 (gphifac)
      ENDIF
      gc(neqs1) = gphifac

!
!     COMPUTE INVARIANT RESIDUALS
!
      r1 = one/(2*r0scale)**2
      jedge = 0    
      delIter = iter2-iter1

!      IF (l_v3fit) THEN
!         IF (iter2-iter1.lt.50) jedge = 1
!      ELSE
         IF (delIter.lt.50 .and.                                        &
            (fsqr+fsqz).LT.1.E-6_dp) jedge = 1
!      ENDIF

      CALL getfsq_par (gcr, gcz, fsqr, fsqz, r1*fnorm, jedge)

      CALL second0(tredon)
      tmp = SUM(gcl(:,:,tlglob:trglob,:)*gcl(:,:,tlglob:trglob,:))
      CALL MPI_Allreduce(tmp,ftotal,1,MPI_REAL8,MPI_SUM,NS_COMM,MPI_ERR)
      CALL second0(tredoff)
      allreduce_time = allreduce_time + (tredoff - tredon)
      fsql = fnormL*ftotal
      IF(rank .EQ. nranks-1) &
        fedge = r1*fnorm*SUM(gcr(:,:,ns,:)**2 + gcz(:,:,ns,:)**2)
!
!     PERFORM PRECONDITIONING AND COMPUTE RESIDUES
!
#if defined(PTESTING)
      IF (ictrl_prec2d.EQ.1 .AND. lactive) THEN
      CALL Gather4XArray(pgc)
      IF (rank.eq.0) THEN
         CALL PrintOutLinearArray(pgc, 1, ns, .FALSE., 1000+nranks)
         PRINT *,'BEFORE precond_par'
      END IF
      CALL MPI_BARRIER(NS_COMM,MPI_ERR)
      END IF
#endif
      IF (ictrl_prec2d .EQ. 1) THEN
         
         IF (l_colscale .AND. lactive) CALL SAXLASTNTYPE(pgc, pcol_scale, pgc)

         LRESIDUECALL=.TRUE.
         CALL block_precond_par(pgc)
         LRESIDUECALL=.FALSE.

         IF (.NOT.lfreeb .AND. ANY(gcr(:,:,ns,:) .NE. zero))            &
            STOP 'gcr(ns) != 0 for fixed boundary in residue'
         IF (.NOT.lfreeb .AND. ANY(gcz(:,:,ns,:) .NE. zero))            &
            STOP 'gcz(ns) != 0 for fixed boundary in residue'
         IF (ANY(gcl(1:,m0,:,zsc) .NE. zero))                           &
            STOP 'gcl(m=0,n>0,sc) != 0 in residue'
         IF (lthreed .AND. ANY(gcl(n0,:,:,zcs) .NE. zero))              &
            STOP 'gcl(n=0,m,cs) != 0 in residue'

         fsqr1 = SUM(gcr*gcr)
         fsqz1 = SUM(gcz*gcz)
         fsql1 = SUM(gcl*gcl)

      ELSE
!        m = 1 constraint scaling

         IF (lthreed) CALL scale_m1_par(gcr(:,m1,:,rss), gcz(:,m1,:,zcs))
         IF (lasym)   CALL scale_m1_par(gcr(:,m1,:,rsc), gcz(:,m1,:,zcc))

         jedge = 0
         CALL scalfor_par (gcr, arm, brm, ard, brd, crd, jedge)
         jedge = 1
         CALL scalfor_par (gcz, azm, bzm, azd, bzd, crd, jedge)

         CALL getfsq_par (gcr, gcz, fsqr1, fsqz1, fnorm1, m1)

         DO l = tlglob, trglob
            gcl(:,:,l,:) = pfaclam(:,:,l,:)*gcl(:,:,l,:)
            tmp2(l) = SUM(gcl(:,:,l,:)**2)
         END DO
         CALL Gather1XArray(tmp2)
         ftotal = SUM(tmp2(2:ns))
         fsql1 = hs*ftotal

         CALL PadSides(pgc)  

      ENDIF

#if defined(PTESTING)
         IF (ictrl_prec2d.eq.1 .and. lactive) THEN
            IF (ictrl_prec2d .NE. 1) CALL Gather4XArray(pgc)
            IF (rank.eq.0)   &
              CALL PrintOutLinearArray(pgc, 1, ns, .TRUE., 2000+nranks)
            CALL MPI_BARRIER(NS_COMM,MPI_ERR)
            STOP 'AFTER precond_par'
         END IF
#endif
      CALL second0 (tresoff)
      residue_time = residue_time + (tresoff-treson)

      END SUBROUTINE residue_par

      SUBROUTINE constrain_m1_par(gcr, gcz)
      USE vmec_main
      USE parallel_include_module
      USE precon2d, ONLY: ictrl_prec2d
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,ns), INTENT(INOUT) :: gcr, gcz
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(dp), PARAMETER :: FThreshold = 1.E-6_dp
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: temp
!-----------------------------------------------
!
!     COMPUTE INTERNAL gr, gz
!     NOTE: internal gz => 0 for both values of lconm1 (although gz is different)
!     FOR lconm1=T, gcr(internal) = gcr+gcz, gcz(internal) = gcr-gcz->0
!
      ALLOCATE(temp(0:ntor,ns))
      IF (lconm1) THEN
         temp(:,tlglob:trglob) = gcr(:,tlglob:trglob)
         gcr(:,tlglob:trglob) = osqrt2*(gcr (:,tlglob:trglob) + gcz(:,tlglob:trglob))
         gcz(:,tlglob:trglob) = osqrt2*(temp(:,tlglob:trglob) - gcz(:,tlglob:trglob))
      END IF

!v8.50: ADD iter2<2 so reset=<WOUT_FILE> works
      IF (fsqz.LT.FThreshold .OR. iter2.LT.2 .OR. ictrl_prec2d.NE.0)   &
         gcz(:,tlglob:trglob) = 0
 
      DEALLOCATE(temp)
      END SUBROUTINE constrain_m1_par

      SUBROUTINE scale_m1_par(gcr, gcz)
      USE vmec_main
      USE parallel_include_module
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,ns), INTENT(inout) :: gcr, gcz
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: nodd=2
      INTEGER :: n
      REAL(dp) :: fac(ns)
!-----------------------------------------------
      IF (.not.lconm1) RETURN

      fac(tlglob:trglob) = (ard(tlglob:trglob,nodd)+brd(tlglob:trglob,nodd))/ &
            (ard(tlglob:trglob,nodd)+brd(tlglob:trglob,nodd)+ &
            azd(tlglob:trglob,nodd)+bzd(tlglob:trglob,nodd))
      DO n = 0, ntor
         gcr(n,tlglob:trglob) = fac(tlglob:trglob)*gcr(n,tlglob:trglob)
      END DO

      fac(tlglob:trglob) = (azd(tlglob:trglob,nodd)+bzd(tlglob:trglob,nodd))/ &
            (ard(tlglob:trglob,nodd)+brd(tlglob:trglob,nodd)+  &
            azd(tlglob:trglob,nodd)+bzd(tlglob:trglob,nodd))
      DO n = 0, ntor
         gcz(n,tlglob:trglob) = fac(tlglob:trglob)*gcz(n,tlglob:trglob)
      END DO
 
      END SUBROUTINE scale_m1_par
#endif

      SUBROUTINE residue (gcr, gcz, gcl)
      USE vmec_main, p5 => cp5
      USE vmec_params, ONLY: rss, zcs, rsc, zcc,                        &
                             meven, modd, ntmax, signgs
      USE realspace, ONLY: phip
      USE vsvd
      USE xstuff
      USE precon2d
#if defined(SKS)
      USE parallel_include_module
#endif
#ifdef _HBANGLE
      USE angle_constraints, ONLY: scalfor_rho
#endif
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) :: &
        gcr, gcz, gcl
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: n0=0, m0=0, m1=1
      INTEGER, PARAMETER :: n3d=0, nasym=1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nsfix, jedge, delIter
      REAL(dp) :: r1, tnorm, fac
      INTEGER  :: i, j, k, l
!-----------------------------------------------
!
!     IMPOSE M=1 MODE CONSTRAINT TO MAKE THETA ANGLE
!     INVARIANT TO PHI-SHIFTS (AND THETA SHIFTS FOR ASYMMETRIC CASE)
!     (ZCS = RSS, ZSS = RCS ARE THE CORRECT POLAR RELATIONS)
!

#if defined(SKS)      
      CALL second0 (treson)
#endif

#ifdef _HBANGLE
!FREE-BDY RFP MAY NEED THIS TO IMPROVE CONVERGENCE (SPH 022514)
      IF (lfreeb .AND. lrfp) THEN
         fac = 0
         IF (ictrl_prec2d .EQ. 0) fac = 1.E-1_dp
         gcz(ns,0,m0,:) = fac*gcz(ns,0,m0,:)
      END IF
#else
!
!     SYMMETRIC PERTURBATIONS (BASED ON POLAR RELATIONS):
!        RSS(n) = ZCS(n), n != 0
!     ASYMMETRIC PERTURBATIONS:
!        RSC(n) = ZCC(n), ALL n
!
!     INTERNALLY:
!        XC(rss) = .5*(Rss + Zcs), XC(zcs) = .5*(Rss - Zcs) -> 0
!        XC(rsc) = .5*(Rsc + Zcc), XC(zcc) = .5*(Rsc - Zcc) -> 0
!     THIS IMPLIES THE CONSTRAINT
!        3D ONLY : GC(zcs) = 0;  
!        ASYM:     GC(zcc) = 0
!

      IF (lthreed) CALL constrain_m1(gcr(:,:,m1,rss), gcz(:,:,m1,zcs))
      IF (lasym)   CALL constrain_m1(gcr(:,:,m1,rsc), gcz(:,:,m1,zcc))

!FREE-BDY RFP MAY NEED THIS TO IMPROVE CONVERGENCE (SPH 022514)
      IF (lfreeb .AND. lrfp) THEN
!      IF (lfreeb .AND. lrfp .AND. iter2 .gt. 2000) THEN
         fac = 0
         IF (ictrl_prec2d .EQ. 0) fac = 1.E-1_dp
         gcr(ns,0,m0,:) = fac*gcr(ns,0,m0,:)
         gcz(ns,0,m0,:) = fac*gcz(ns,0,m0,:)
      END IF
#endif

!     PRECONDITIONER MUST BE CALCULATED USING RAW (UNPRECONDITIONED) FORCES
      IF (ictrl_prec2d.GE.2 .OR. ictrl_prec2d.EQ.-1) RETURN

!
!     PUT FORCES INTO PHIFSAVE UNITS USED BY PRECONDITIONERS, FNORM
!
      IF (phifac .eq. zero) THEN
         STOP 'phifac = 0 in residue'
      ELSE
         tnorm = phifsave/phifac           !put all forces into phifac=phifsave units
      END IF

      IF (lrecon) THEN
!
!       MOVE R(n=0,m=0) TO SATISFY LIMITER OR AXIS POSITION
!       USE XC(NEQS2) TO STORE THIS PERTURBATION
!       TO SATISFY FORCE BALANCE AT JS=1, ADJUST PFAC IN RADFOR
!       ALSO, SCALE TOROIDAL FLUX AT EDGE TO MATCH MINOR RADIUS

        r1 = SUM(gcr(:ns,n0,m0,1))
        fsqsum0 = signgs*hs*r1/r0scale
        nsfix = 1                   !fix origin for reconstruction mode
        gcr = gcr * tnorm**2
        gcz = gcz * tnorm**2
        gcl = gcl * tnorm
        IF (iopt_raxis.gt.0 .and. iresidue.eq.2                        &
           .and. fsq.lt.fopt_axis) iresidue = 3
        IF (iresidue .lt. 3) gcr(nsfix,n0,m0,1) = zero
      ELSE
!
!     ADJUST PHIEDGE
!
         IF (imovephi .gt. 0) CALL movephi1 (gphifac)
      ENDIF
      gc(neqs1) = gphifac
!
!     COMPUTE INVARIANT RESIDUALS
!
      r1 = one/(2*r0scale)**2
      jedge = 0    
!SPH-JAH013108: MUST INCLUDE EDGE FORCE (INITIALLY) FOR V3FITA TO WORK
!ADD A V3FIT RELATED FLAG? ADD fsq criterion first
      delIter = iter2-iter1

      IF (l_v3fit) THEN
!  Coding for when run by V3FIT. Needed for correct computation
!  of partial derivatives
         IF (iter2-iter1.lt.50) jedge = 1
      ELSE
!  Coding for VMEC2000 run stand-alone
         IF (delIter.lt.50 .and.                                        &
            (fsqr+fsqz).lt.1.E-6_dp) jedge = 1
      ENDIF

      CALL getfsq (gcr, gcz, fsqr, fsqz, r1*fnorm, jedge)

      fsql = fnormL*SUM(gcl*gcl)
      fedge = r1*fnorm*SUM(gcr(ns,:,:,:)**2 + gcz(ns,:,:,:)**2)

!
!     PERFORM PRECONDITIONING AND COMPUTE RESIDUES
!

      IF (ictrl_prec2d .EQ. 1) THEN

         CALL block_precond(gc)

         IF (.not.lfreeb .and. ANY(gcr(ns,:,:,:) .ne. zero))            &
            STOP 'gcr(ns) != 0 for fixed boundary in residue'
         IF (.not.lfreeb .and. ANY(gcz(ns,:,:,:) .ne. zero))            &
            STOP 'gcz(ns) != 0 for fixed boundary in residue'
         IF (ANY(gcl(:,1:,0,zsc) .ne. zero))                            &
            STOP 'gcl(m=0,n>0,sc) != 0 in residue'
         IF (lthreed) THEN
            IF (ANY(gcl(:,n0,:,zcs) .ne. zero))                         &
            STOP 'gcl(n=0,m,cs) != 0 in residue'
         END IF

         fsqr1 = SUM(gcr*gcr)
         fsqz1 = SUM(gcz*gcz)
         fsql1 = SUM(gcl*gcl)

      ELSE
#ifdef _HBANGLE
         CALL scalfor_rho(gcr, gcz)
#else

!        m = 1 constraint scaling
         IF (lthreed) CALL scale_m1(gcr(:,:,m1,rss), gcz(:,:,m1,zcs))
         IF (lasym)   CALL scale_m1(gcr(:,:,m1,rsc), gcz(:,:,m1,zcc))

         jedge = 0
         CALL scalfor (gcr, arm, brm, ard, brd, crd, jedge)
         jedge = 1
         CALL scalfor (gcz, azm, bzm, azd, bzd, crd, jedge)

#endif

!SPH: add fnorm1 ~ 1/R**2, since preconditioned forces gcr,gcz ~ Rmn or Zmn
         CALL getfsq (gcr, gcz, fsqr1, fsqz1, fnorm1, m1)
#ifdef _HBANGLE
!
!     TO IMPROVE CONVERGENCE, REDUCE FORCES INITIALLY IF THEY ARE TOO LARGE
!
         fac = .5_dp
         IF ((iter2-iter1).LT.25 .AND. (fsqr+fsqz).GT.1.E-2_dp)         &
         fac = fac / SQRT(1.E2_dp*(fsqr+fsqz))
         gcr = fac*gcr; gcz = fac*gcz 
#endif

!SPH: THIS IS NOT INVARIANT UNDER PHIP->A*PHIP, AM->A**2*AM IN PROFIL1D
!     (EXTCUR -> A*EXTCUR for FREE BOUNDARY)

         gcl = faclam*gcl
         fsql1 = hs*SUM(gcl*gcl)
!030514      fsql1 = hs*lamscale**2*SUM(gcl*gcl)

      ENDIF

#if defined(SKS)      
      CALL second0 (tresoff)
      s_residue_time = s_residue_time + (tresoff-treson)
#endif

      END SUBROUTINE residue

      SUBROUTINE constrain_m1(gcr, gcz)
      USE vmec_main
      USE precon2d, ONLY: ictrl_prec2d
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor), INTENT(INOUT) :: gcr, gcz
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(dp), PARAMETER :: FThreshold = 1.E-6_dp
      REAL(dp) :: temp(ns,0:ntor)
!-----------------------------------------------
!
!     COMPUTE INTERNAL gr, gz
!     NOTE: internal gz => 0 for both values of lconm1 (although gz is different)
!     FOR lconm1=T, gcr(internal) = gcr+gcz, gcz(internal) = gcr-gcz->0
!
      IF (lconm1) THEN
         temp = gcr
         gcr = osqrt2*(gcr + gcz)
         gcz = osqrt2*(temp - gcz)
      END IF

!v8.50: ADD iter2<2 so reset=<WOUT_FILE> works
      IF (fsqz.LT.FThreshold .OR. iter2.LT.2 .OR. ictrl_prec2d.NE.0)   &
        gcz = 0
 
      END SUBROUTINE constrain_m1

      SUBROUTINE scale_m1(gcr, gcz)
      USE vmec_main
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor), INTENT(inout) :: gcr, gcz
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: nodd=2
      INTEGER :: n
      REAL(dp) :: fac(ns)
!-----------------------------------------------
      IF (.not.lconm1) RETURN

      fac = (ard(:,nodd)+brd(:,nodd))/                                        &
            (ard(:,nodd)+brd(:,nodd)+azd(:,nodd)+bzd(:,nodd))
      DO n = 0, ntor
         gcr(:,n) = fac*gcr(:,n)
      END DO

      fac = (azd(:,nodd)+bzd(:,nodd))/                                        &
            (ard(:,nodd)+brd(:,nodd)+azd(:,nodd)+bzd(:,nodd))
      DO n = 0, ntor
         gcz(:,n) = fac*gcz(:,n)
      END DO
 
      END SUBROUTINE scale_m1

