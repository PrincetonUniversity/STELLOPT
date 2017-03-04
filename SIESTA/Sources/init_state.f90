      SUBROUTINE init_state_par(lcurrent_only)
!     
!     WRITTEN 01-12-07 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes initial values of equilibrium pressure, fields
!              needed to evaluate forces to update those quantities
!              NOT called inside of linearization loops (where these are
!              FIXED quantities)
!
#if defined(SKS)
       USE stel_kinds
       USE quantities
       USE diagnostics_mod
       USE evolution, ONLY: ste
       USE shared_data, ONLY: fsq_total, l_getfsq, l_update_state
       USE fourier, ONLY: toijsp_par, tomnsp_par
       USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
       USE perturbation, ONLY: buv_res
       USE hessian, ONLY: muPar
       IMPLICIT NONE
       INCLUDE 'mpif.h'
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       LOGICAL, INTENT(IN)      :: lcurrent_only
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER :: ifull=0, ihalf=1, m0=0, m1=1
       LOGICAL, PARAMETER :: lpar=.TRUE.
       INTEGER            :: istat, m, nblock
       REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: work1,             &
                    pijch_du, pijch_dv 
       REAL(rprec) :: ton, toff, stes(2), rtes(2)
       INTEGER     :: i, nsmin, nsmax, ns_span, n1, n2        ,js
       LOGICAL     :: ladd_pert, lcurr
!-----------------------------------------------
       istat = 0
       nsmin=MAX(1,startglobrow-1); nsmax=MIN(endglobrow+2,ns)
       ns_span = nsmax-nsmin+1

       CALL Init_Allocate_Arrays (lpar)

! Calculate unperturbed half mesh jac*bsupX
! SPH030410: need to get values at js=1, where UPDATE_BFIELD stored origin values
       DO i=nsmin, nsmax
          jbsupsmnsh_p(:,:,i) = jbsupsmnsh(:,:,i)
       END DO
       CALL toijsp_par(jbsupsmnsh_p, bsupsijsh0, 0, 0, 1, ifull, 1, ns_span)
       DO i=nsmin, nsmax
          jbsupumnch_p(:,:,i) = jbsupumnch(:,:,i)
       END DO
       CALL toijsp_par(jbsupumnch_p, bsupuijch0, 0, 0, 0, ifull, 1, ns_span)
       DO i=nsmin, nsmax
          jbsupvmnch_p(:,:,i) = jbsupvmnch(:,:,i)
       END DO
       CALL toijsp_par(jbsupvmnch_p, bsupvijch0, 0, 0, 0, ifull, 1, ns_span)
       DO i=nsmin, nsmax
          jpmnch_p(:,:,i) = jpmnch(:,:,i)
       END DO

! Convert jac*bsubX to bsupX in real space


       bsupsijsh0 = bsupsijsh0/jacobh_p(:,:,nsmin:nsmax)
       bsupvijch0 = bsupvijch0/jacobh_p(:,:,nsmin:nsmax)

!      On exit, bsupXF valid on [nsmin, nsmax-1]
       CALL bhalftobfull_par (bsupsijsh0, bsupuijch0, bsupvijch0,       &
                              bsupsijsf0, bsupuijcf0, bsupvijcf0,       &
                              nsmin, nsmax)

!Compute and store unperturbed current components on full mesh, KsupX = sqrt(g)*JsupX
       ladd_pert = ALLOCATED(buv_res)
       IF (ladd_pert) THEN
          nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow+1,ns)
       END IF

       lcurr = lcurrent_only .OR. (.NOT.ladd_pert)
       CALL cv_currents_par (bsupsijsh0(:,:,nsmin:nsmax),               &
                             bsupuijch0(:,:,nsmin:nsmax),               &
                             bsupvijch0(:,:,nsmin:nsmax),               &
                             ksupsijsf0(:,:,nsmin:nsmax),               &
                             ksupuijcf0(:,:,nsmin:nsmax),               &
                             ksupvijcf0(:,:,nsmin:nsmax),               &
                             l_getfsq, lcurr )
       IF (ladd_pert) THEN
          nsmin=MAX(1,startglobrow-1); nsmax=MIN(endglobrow+2,ns)
       END IF

       IF (lcurrent_only) RETURN

!Need this for resistive diffusion (|| resistivity) and mupar damping
       IF (ladd_pert .OR. muPar.NE.zero) THEN
          nsmax=MIN(endglobrow+1,ns)
          CALL tolowerf_par(bsupsijsf0,bsupuijcf0,bsupvijcf0,           & ! Calculate K || B
                            bsubsijsf, bsubuijcf, bsubvijcf,            &
                            nsmin, nsmax)  
!     IN UPDATE-BFIELD, KSUB = jacob*JSUB: PARALLEL CURRENT
          IF (ladd_pert) THEN
          ksubsijsf(:,:,nsmin:nsmax) = jacobf_p(:,:,nsmin:nsmax)*bsubsijsf(:,:,nsmin:nsmax)
          ksubuijcf(:,:,nsmin:nsmax) = jacobf_p(:,:,nsmin:nsmax)*bsubuijcf(:,:,nsmin:nsmax)
          ksubvijcf(:,:,nsmin:nsmax) = jacobf_p(:,:,nsmin:nsmax)*bsubvijcf(:,:,nsmin:nsmax)
          END IF
          nsmax=MIN(endglobrow+2,ns)
       END IF

! Calculate pressure in real space (half mesh) - need js=1 point here
       CALL toijsp_par(jpmnch_p, pijch, 0, 0, 0, ifull, 1, ns_span)
       
!      Get real pressure from pijch = jacobh * real-pressure
       pijch = pijch/jacobh_p(:,:,nsmin:nsmax)
 
!      pressure on full mesh (extrap to origin stored in pijch(1)
!      On Exit, pijcf valid on [nsmin:nsmax-1]
       nblock=SIZE(pijch,1)*SIZE(pijch,2)
       CALL to_full_mesh_par(pijch, pijcf, nblock, nsmin, nsmax)           ! Get pressure on full mesh

!      Get true pressure p on half grid needed for dp/du and dp/dv
       ALLOCATE (work1(0:mpol,-ntor:ntor,nsmin:nsmax),                  &
                 pijch_du(ntheta,nzeta,nsmin:nsmax),                    &
                 pijch_dv(ntheta,nzeta,nsmin:nsmax), stat=istat)
       IF (istat .NE. 0) STOP 'Allocation error in INIT_STATE_PAR'

!      work1 now contains filtered (to MPOL, NTOR) fourier harmonics of p on half mesh
       CALL tomnsp_par(pijch, work1(:,:,nsmin:nsmax), 0, ifull, 1, ns_span)

       IF (l_update_state) THEN
          pijch_du = pijch                                                ! unfiltered
          CALL toijsp_par(work1, pijch_dv, 0, 0, 0, ihalf, 1, ns_span)    ! p (fourier filtered)
          n1=MAX(1,startglobrow); n2=MIN(endglobrow,ns)
          IF (n1 .EQ. 1) THEN
             pijch_du(:,:,1) = 0;   pijch_dv(:,:,1) = 0
          END IF
          stes(1) = SUM((pijch_du(:,:,n1:n2)-pijch_dv(:,:,n1:n2))**2)
          stes(2) = SUM((pijch_du(:,:,n1:n2)+pijch_dv(:,:,n1:n2))**2)
          CALL MPI_ALLREDUCE(stes,rtes,2,MPI_REAL8,MPI_SUM,              &
                             MPI_COMM_WORLD,MPI_ERR)
          IF (rtes(2) .NE. zero) ste(1) = SQRT(rtes(1)/rtes(2))
       END IF

!      compute radial derivative at full mesh points
!      On Exit, pijcf_ds valid on [nsmin:nsmax-1]
       CALL GradientFull_par(pijcf_ds, pijch, nblock, nsmin, nsmax)

!SPH:  one-sided (m=1) derivative at origin yields factor of 2 (hs/2)
!!       CALL getorigin(pijch, m1, 0)
       IF (nsmin .EQ. 1) THEN
          pijcf_ds(:,:,1) = 2*pijch(:,:,1)*ohs
          pijch(:,:,1) = 0 
       END IF

!      spectrally filtered angle derivatives
       CALL toijsp_par(work1, pijch_du, 1, 0, 0, ihalf, 1, ns_span)      ! dp/dtheta in real space (half mesh)     
       CALL toijsp_par(work1, pijch_dv, 0, 1, 0, ifull, 1, ns_span)      ! dp/dzeta in real space (half mesh, include js=1 extrap)

       DEALLOCATE(work1)                 
       
!SPH   AVERAGE pijch_du,dv to FULL mesh and set origin value (SIN SYMMETRY)
!      On Exit, pijcf_du,v valid on [nsmin:nsmax-1]
       CALL to_full_mesh_par(pijch_du, pijcf_du, nblock, nsmin, nsmax)
       CALL to_full_mesh_par(pijch_dv, pijcf_dv, nblock, nsmin, nsmax)

       DEALLOCATE(pijch_du, pijch_dv)                 
#endif
       END SUBROUTINE init_state_par

       SUBROUTINE init_state(lcurrent_only)
!     
!     WRITTEN 01-12-07 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes initial values of equilibrium pressure, fields
!              needed to evaluate forces to update those quantities
!              NOT called inside of linearization loops (where these are
!              FIXED quantities)
!
       USE stel_kinds
       USE quantities
       USE diagnostics_mod
       USE evolution, ONLY: ste
       USE shared_data, ONLY: fsq_total, l_getfsq, l_update_state
       USE perturbation, ONLY: buv_res
       USE hessian, ONLY: muPar
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       LOGICAL, INTENT(in)      :: lcurrent_only
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER :: ifull=0, ihalf=1, m0=0, m1=1
       LOGICAL, PARAMETER :: lpar=.FALSE.
       INTEGER            :: istat, m, nblock                                   !, mpol_s, ntor_s
       REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: work1,             &
                    pijch_du, pijch_dv 
       LOGICAL            :: ladd_pert, lcurr
!-----------------------------------------------
       istat = 0

       CALL Init_Allocate_Arrays (lpar)

! Calculate unperturbed half mesh jac*bsupX
! SPH030410: need to get values at js=1, where UPDATE_BFIELD stored origin values

       CALL toijsp(jbsupsmnsh, bsupsijsh0, 0, 0, 1, ifull)
       CALL toijsp(jbsupumnch, bsupuijch0, 0, 0, 0, ifull)
       CALL toijsp(jbsupvmnch, bsupvijch0, 0, 0, 0, ifull)

! Convert jac*bsubX to bsupX in real space

       bsupsijsh0 = bsupsijsh0/jacobh
       bsupvijch0 = bsupvijch0/jacobh

       CALL bhalftobfull (bsupsijsh0, bsupuijch0, bsupvijch0,           &
                          bsupsijsf0, bsupuijcf0, bsupvijcf0)
       IF (ANY(bsupsijsf0(:,:,ns) .NE. zero)) PRINT *,'bsupsijsf0(ns) != 0'


!Compute and store unperturbed current components on full mesh, KsupX = sqrt(g)*JsupX
       ladd_pert = ALLOCATED(buv_res)
       lcurr = lcurrent_only .OR. (.NOT.ladd_pert)
       CALL cv_currents (bsupsijsh0, bsupuijch0, bsupvijch0,            &
                         ksupsijsf0, ksupuijcf0, ksupvijcf0,            &
                         l_getfsq, lcurr)

       IF (lcurrent_only) RETURN

!Need this for resistive diffusion (|| resistivity) and mupar damping
       IF (ladd_pert .OR. muPar.NE.zero) THEN
          CALL tolowerf(bsupsijsf0, bsupuijcf0, bsupvijcf0,             & ! Calculate K || B
                        bsubsijsf,  bsubuijcf,  bsubvijcf)  
!     IN UPDATE-BFIELD, KSUB = jacob*JSUB
          IF (ladd_pert) THEN
             ksubsijsf = jacobf*bsubsijsf
             ksubuijcf = jacobf*bsubuijcf
             ksubvijcf = jacobf*bsubvijcf
          END IF
       END IF

! Calculate pressure in real space (half mesh) - need js=1 point here
       CALL toijsp(jpmnch, pijch, 0, 0, 0, ifull)
              
!      Get real pressure from pijch = jacobh * real-pressure
       pijch = pijch/jacobh

!      pressure on full mesh (extrap to origin stored in pijch(1)
       nblock=SIZE(pijch,1)*SIZE(pijch,2)
       CALL to_full_mesh(pijch, pijcf, nblock)                               ! Get pressure on full mesh
       pijcf(:,:,1) = pijch(:,:,1)   
       pijcf(:,:,ns) = pijch(:,:,ns)                                    ! 1pt "extrap"

!      Get true pressure p on half grid needed for dp/du and dp/dv
       ALLOCATE (work1(0:mpol,-ntor:ntor,ns),                           &
                 pijch_du(ntheta,nzeta,ns),                             &
                 pijch_dv(ntheta,nzeta,ns), stat=istat)
       IF (istat .ne. 0) STOP 'Allocation error in INIT_STATE'

!      work1 now contains filtered (to MPOL, NTOR) fourier harmonics of p on half mesh
       CALL tomnsp(pijch, work1, 0, ifull)

!      spectral filtering of p (compute spectral truncation error diagnostic)
       IF (l_update_state) THEN
          pijch_du = pijch                                                 ! unfiltered
          CALL toijsp(work1, pijch_dv, 0, 0, 0, ihalf)                     ! p (fourier filtered)
          pijch_du(:,:,1) = 0;   pijch_dv(:,:,1) = 0
          ste(1) = SUM((pijch_du+pijch_dv)**2)
          IF (ste(1) .NE. zero) ste(1) = SQRT(SUM((pijch_du-pijch_dv)**2)                     &
                                       /      ste(1))
       END IF

!      compute radial derivative at full mesh points
       CALL GradientFull(pijcf_ds, pijch, nblock)
!!       CALL rad_deriv_sp(pijch, pijcf_ds)
!SPH:  one-sided (m=1) derivative at origin yields factor of 2 (hs/2)
!!       CALL getorigin(pijch, m1, 0)
       pijcf_ds(:,:,1) = 2*pijch(:,:,1)*ohs
       pijch(:,:,1) = 0 

!      spectrally filtered angle derivatives
       CALL toijsp(work1, pijch_du, 1, 0, 0, ihalf)                     ! dp/dtheta in real space (half mesh)     
       CALL toijsp(work1, pijch_dv, 0, 1, 0, ifull)                     ! dp/dzeta in real space (half mesh, include js=1 extrap)
       DEALLOCATE(work1)                 
       
!SPH   AVERAGE pijch_du,dv to FULL mesh and set origin value (SIN SYMMETRY)
       CALL to_full_mesh(pijch_du, pijcf_du, nblock)
       pijcf_du(:,:,ns) = pijch_du(:,:,ns)
       pijcf_du(:,:,1) = pijch_du(:,:,2)        !!This is zero, but mult by v^u is finite for m=1

       CALL to_full_mesh(pijch_dv, pijcf_dv, nblock)
       pijcf_dv(:,:,ns) = pijch_dv(:,:,ns)      ! evaluate edge forces at Ns - 1/2 point
       pijcf_dv(:,:,1) = pijch_dv(:,:,1)

       pijch(:,:,1) = 0 

       DEALLOCATE(pijch_du, pijch_dv)                 

       END SUBROUTINE init_state

       SUBROUTINE Init_Allocate_Arrays(lpar)
       USE quantities
       USE shared_data, ONLY: l_par_state
       USE nscalingtools, ONLY: startglobrow, endglobrow
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       LOGICAL, INTENT(in)     :: lpar
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       LOGICAL                 :: lrealloc, lalloc
       INTEGER                 :: istat, n1, n2, n3, nsmin, nsmax,      &
                                  nsmin1, nsmax1
!-----------------------------------------------
       lalloc = ALLOCATED(jvsupsijcf) 
       lrealloc = .TRUE.
       istat = 0
       nsmin =MAX(1,startglobrow);   nsmax =MIN(endglobrow+1,ns)
       nsmin1=MAX(1,startglobrow-1); nsmax1=MIN(endglobrow+2,ns)

       IF (lpar) THEN
          n1 = ntheta; n2 = nzeta; n3 = nsmax1-nsmin1+1
          IF (.NOT.ALLOCATED(jbsupsmnsh_p)) THEN
             ALLOCATE(jbsupsmnsh_p(0:mpol,-ntor:ntor,nsmin1:nsmax1),    &
                      jbsupumnch_p(0:mpol,-ntor:ntor,nsmin1:nsmax1),    &
                      jbsupvmnch_p(0:mpol,-ntor:ntor,nsmin1:nsmax1),    &
                      jpmnch_p    (0:mpol,-ntor:ntor,nsmin1:nsmax1),    &
                      stat=istat)
             IF (istat.ne.0) STOP 'Allocation failed in ALLOC_QUANTITIES'    
          END IF
          l_par_state=.TRUE.
       ELSE
          n1 = ntheta; n2 = nzeta; n3 = ns
          l_par_state=.FALSE.
       END IF

       IF (lalloc) THEN
          IF (SIZE(bsupsijsh0,3) .EQ. n3) lrealloc=.FALSE.
          IF (lrealloc) THEN
          DEALLOCATE (jvsupsijcf, jvsupuijsf, jvsupvijsf,               &
                      bsupsijsf0, bsupuijcf0, bsupvijcf0,               &
                      bsupsijsh0, bsupuijch0, bsupvijch0,               &
                      bsubsijsf,  bsubuijcf,  bsubvijcf,                &
                      ksupsijsf0, ksupuijcf0, ksupvijcf0,               &
                      ksubsijsf,  ksubuijcf,  ksubvijcf,                &
                      pijcf, pijch, pijcf_ds, pijcf_du, pijcf_dv,       &
                      stat=istat)
          IF (istat .NE. 0) STOP 'Deallocation failed in ALLOC_QUANTITIES'
          END IF
       END IF

       IF (lrealloc) THEN
          n3 = ns
          ALLOCATE (jvsupsijcf(n1,n2,n3), jvsupuijsf(n1,n2,n3),         &  ! Full mesh quantities (real space)
                    jvsupvijsf(n1,n2,n3))                                  ! V_s (cosine), V_u (sine), V_v (sine)
          IF (lpar) THEN
          ALLOCATE (ksupsijsf0(n1,n2,nsmin1:nsmax1),                    &  ! K^s (sine), K^u (cosine), K^v (cosine)
                    ksupuijcf0(n1,n2,nsmin1:nsmax1),                    &
                    ksupvijcf0(n1,n2,nsmin1:nsmax1),                    &
                    bsupsijsf0(n1,n2,nsmin1:nsmax1),                    &  ! B^s (sin), B^u (cos), B^v (cos)
                    bsupuijcf0(n1,n2,nsmin1:nsmax1),                    &
                    bsupvijcf0(n1,n2,nsmin1:nsmax1),                    &
                    bsupsijsh0(n1,n2,nsmin1:nsmax1),                    &  ! B^s (sin), B^u (cos), B^v (cos)
                    bsupuijch0(n1,n2,nsmin1:nsmax1),                    &
                    bsupvijch0(n1,n2,nsmin1:nsmax1),                    &
                    ksubsijsf(n1,n2,nsmin1:nsmax),                      &  ! Full mesh quantities (real space)
                    ksubuijcf(n1,n2,nsmin1:nsmax),                      &  ! K^s (sine), K^u (cosine), K^v (cosine)
                    ksubvijcf(n1,n2,nsmin1:nsmax),                      &
                    bsubsijsf(n1,n2,nsmin1:nsmax),                      &  ! Full mesh quantities (real space)
                    bsubuijcf(n1,n2,nsmin1:nsmax),                      &  ! K^s (sine), K^u (cosine), K^v (cosine)
                    bsubvijcf(n1,n2,nsmin1:nsmax),                      &
                    pijch(n1,n2,nsmin1:nsmax1),                         &
                    pijcf(n1,n2,nsmin1:nsmax1),                         &
                    pijcf_ds(n1,n2,nsmin1:nsmax1),                      &
                    pijcf_du(n1,n2,nsmin1:nsmax1),                      &
                    pijcf_dv(n1,n2,nsmin1:nsmax1), stat=istat)
          ELSE
          ALLOCATE (ksupsijsf0(n1,n2,n3), ksupuijcf0(n1,n2,n3),         &
                    ksupvijcf0(n1,n2,n3),                               &
                    ksubsijsf(n1,n2,n3),                                &
                    ksubuijcf(n1,n2,n3),                                &
                    ksubvijcf(n1,n2,n3),                                &
                    bsubsijsf(n1,n2,n3),                                &
                    bsubuijcf(n1,n2,n3),                                &
                    bsubvijcf(n1,n2,n3),                                &
                    pijcf(n1,n2,n3), pijch(n1,n2,n3),                   &
                    pijcf_ds(n1,n2,n3), pijcf_du(n1,n2,n3),             &
                    pijcf_dv(n1,n2,n3), bsupsijsf0(n1,n2,n3),           &
                    bsupuijcf0(n1,n2,n3),bsupvijcf0(n1,n2,n3),          &
                    bsupsijsh0(n1,n2,n3),                               &
                    bsupuijch0(n1,n2,n3),bsupvijch0(n1,n2,n3),          &
                    stat=istat)
          END IF
          IF (istat .EQ. 0) THEN
             jvsupsijcf = 0; jvsupuijsf = 0; jvsupvijsf = 0                ! Need in add_resistivity loop
          ELSE
             STOP 'Allocation error in Allocate_Arrays'
          END IF
       END IF


       lalloc = ALLOCATED(ksupsmnsf) 
       lrealloc = .TRUE.

       IF (lpar) THEN
          n1 = mpol+1; n2 = 2*ntor+1; n3 = nsmax-nsmin1+1
       ELSE
          n1 = mpol+1; n2 = 2*ntor+1; n3 = ns
       END IF

       IF (lalloc) THEN
          IF (SIZE(ksupsmnsf,3) .EQ. n3) lrealloc=.FALSE.
          IF (lrealloc) THEN
          DEALLOCATE (ksupsmnsf, ksupumncf, ksupvmncf,                  &
                      djpmnch, djbsupsmnsh, djbsupumnch, djbsupvmnch,   &
                      stat=istat)
          IF (istat .NE. 0) STOP 'Deallocation failed in ALLOC_QUANTITIES'
          END IF
       END IF

       IF (lrealloc) THEN
          IF (lpar) THEN
             ALLOCATE (ksupsmnsf(0:mpol,-ntor:ntor,nsmin1:nsmax),       & 
                       ksupumncf(0:mpol,-ntor:ntor,nsmin1:nsmax),       &
                       ksupvmncf(0:mpol,-ntor:ntor,nsmin1:nsmax),       &
                       djpmnch(0:mpol,-ntor:ntor,nsmin:nsmax),          &
                       djbsupsmnsh(0:mpol,-ntor:ntor,nsmin:nsmax),      &
                       djbsupumnch(0:mpol,-ntor:ntor,nsmin:nsmax),      &
                       djbsupvmnch(0:mpol,-ntor:ntor,nsmin:nsmax),      &
                       stat=istat)
          ELSE
             ALLOCATE (ksupsmnsf(0:mpol,-ntor:ntor,n3),                 & 
                       ksupumncf(0:mpol,-ntor:ntor,n3),                 &
                       ksupvmncf(0:mpol,-ntor:ntor,n3),                 & 
                       djpmnch(0:mpol,-ntor:ntor,n3),                   &
                       djbsupsmnsh(0:mpol,-ntor:ntor,n3),               &
                       djbsupumnch(0:mpol,-ntor:ntor,n3),               &
                       djbsupvmnch(0:mpol,-ntor:ntor,n3),               &
                       stat=istat)
          END IF
          IF (istat .EQ. 0) THEN
             djpmnch = 0; djbsupsmnsh = 0; djbsupumnch = 0; djbsupvmnch = 0 
          ELSE
             STOP 'Allocation error in Allocate_Arrays'
          END IF
       END IF


       END SUBROUTINE Init_Allocate_Arrays
