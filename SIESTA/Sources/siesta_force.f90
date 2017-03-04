!>  \brief Module contained subroutines for updating FMHD=JXB - gradp as part of the SIESTA project.
!!   Upadates FMHD using the advanced values of B and p obtained from calls to UPDATE_PRES and UPDATE_BFIELD
!!   Stores values of FSUB?MNH in module quantities
!!   Linearized force is Flin ~ delta_v has ONLY the linear part (not DC part and not non-linear)
!!  \author S. P. Hirshman and R. Sanchez
!!  \date Aug 21, 2006
      MODULE siesta_force

      USE stel_kinds
      USE stel_constants
      USE fourier 
      USE quantities
      USE island_params, ONLY: nu_i, hs_i, dnorm_i
      USE timer_mod
      USE shared_data, ONLY: l_linearize, l_getwmhd, xPosDef
      USE hessian, ONLY: l_Compute_Hessian

      CONTAINS
#if defined(SKS)
!>  \brief Parallel subroutine to update MHD force
      SUBROUTINE update_force_par
!     
!     PURPOSE: UPDATES FORCE using the advanced values of B and p obtained
!              from calls to UPDATE_PRES and UPDATE_BFIELD
!              Updates values of FSUB?MNH and stores in QUANTITIES MODULE
!
!              Linearized force is Flin ~ delta_v has ONLY the linear part (not DC part and not non-linear)
!
       USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
       IMPLICIT NONE
       INCLUDE 'mpif.h'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER  :: istat, l
       INTEGER, PARAMETER :: ifull=0, ihalf=1, m1=1, m0=0, pcos=0, psin=1
       REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) ::                    &
         pmncf_ds, pmncf, pmnch, jpmnch_i,                              &
         bsupsmnsh_i, bsupumnch_i, bsupvmnch_i,                         &
         ksupsijsf, ksupuijcf, ksupvijcf,                               &
         bsupsijsf1,bsupuijcf1,bsupvijcf1,                              &
         work4, work5, work6

       REAL(rprec), POINTER, DIMENSION(:,:,:) :: work1, work2, work3
       REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:), TARGET ::            &
         bsupsijsh, bsupuijch, bsupvijch
       REAL(rprec) :: ton, toff, skston, skstoff, xfact, wp0
       INTEGER     :: nsmin, nsmax, ns_span, nblock
!-----------------------------------------------
       CALL second0(skston)
       nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow+1,ns)
       ns_span = nsmax-nsmin+1
!
!     Allocate space for TOTAL values of BSUPXMN?H_i and pmnch_i
!     (they are actually updated in UPDATE_STATE, not here, so UNPERTURBED
!      values are not changed until UPDATE_STATE is called)
!     increments are calculated in UPDATE_BFIELD and UPDATE_PRES
!     
       ALLOCATE(bsupsmnsh_i(0:mpol,-ntor:ntor,nsmin:nsmax),             &
                bsupumnch_i(0:mpol,-ntor:ntor,nsmin:nsmax),             &
                bsupvmnch_i(0:mpol,-ntor:ntor,nsmin:nsmax),             &  
                jpmnch_i   (0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'
!
!      IF nonlinear force, add perturbation to DC (initial) value
!      IF linear force, compute JUST perturbation here (add gc0 back in funct_island)
!
       xfact = 1
       IF (l_linearize) xfact = 0
       bsupsmnsh_i = xfact*jbsupsmnsh_p(:,:,nsmin:nsmax) + djbsupsmnsh
       bsupumnch_i = xfact*jbsupumnch_p(:,:,nsmin:nsmax) + djbsupumnch
       bsupvmnch_i = xfact*jbsupvmnch_p(:,:,nsmin:nsmax) + djbsupvmnch
       jpmnch_i    = xfact*jpmnch_p(:,:,nsmin:nsmax)     + djpmnch

!SPH:  Convert to real covariant components

       ALLOCATE(bsupsijsh(ntheta,nzeta,nsmin:nsmax),                   &
                bsupuijch(ntheta,nzeta,nsmin:nsmax),                   &
                bsupvijch(ntheta,nzeta,nsmin:nsmax), stat=istat)  
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE_PAR'

!NEED VALUES AT js=1 (STORED VALUES AT ORIGIN IN UPDATE_BFIELD)

       CALL toijsp_par(bsupsmnsh_i, bsupsijsh, 0, 0, psin, ifull, 1, ns_span)             ! Calculate BSUP on half mesh
       CALL toijsp_par(bsupumnch_i, bsupuijch, 0, 0, pcos, ifull, 1, ns_span)
       CALL toijsp_par(bsupvmnch_i, bsupvijch, 0, 0, pcos, ifull, 1, ns_span)

       bsupsijsh = bsupsijsh/jacobh_p(:,:,nsmin:nsmax)
       bsupvijch = bsupvijch/jacobh_p(:,:,nsmin:nsmax)

       DEALLOCATE(bsupsmnsh_i, bsupumnch_i, bsupvmnch_i, stat=istat)
!
!      Calculate contravariant components of J = curl B (Ksup = gsqrt*J)
!      and update magnetic energy if nonlinear force is being computed
!                             
       ALLOCATE(ksupsijsf(ntheta,nzeta,nsmin:nsmax),                    & ! Full mesh quantities (real space)
                ksupuijcf(ntheta,nzeta,nsmin:nsmax),                    & ! K_s (sine), K_u (cosine), K_v (cosine)
                ksupvijcf(ntheta,nzeta,nsmin:nsmax),                    &
                bsupsijsf1(ntheta,nzeta,nsmin:nsmax),                   &
                bsupuijcf1(ntheta,nzeta,nsmin:nsmax),                   &
                bsupvijcf1(ntheta,nzeta,nsmin:nsmax), stat=istat)
       IF (istat.ne.0) STOP 'Allocation failed in UPDATE_FORCE'

!
!      Calculate updated FULL-MESH components of BSUPX in real-space
!
       CALL bhalftobfull_par (bsupsijsh, bsupuijch, bsupvijch,          &
                              bsupsijsf1,bsupuijcf1,bsupvijcf1,         &
                              nsmin, nsmax)

!Origin currents are X1/2 (1/2 factor in grad-p too, in update_pres)
       CALL cv_currents_par(bsupsijsh, bsupuijch, bsupvijch,            &
                            ksupsijsf, ksupuijcf, ksupvijcf,            &
                           (.NOT.l_linearize .OR. l_getwmhd), .FALSE.)

!      alias work arrays for use as temp work space
       work1=>bsupsijsh; work2=>bsupuijch;  work3=>bsupvijch

       nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow,ns)
       ns_span = nsmax-nsmin+1

       IF (l_linearize) THEN
!      dJ X B0 (note here, ksupXij?f is the perturbed current)
          work1(:,:,nsmin:nsmax) =  bsupvijcf0(:,:,nsmin:nsmax)*ksupuijcf(:,:,nsmin:nsmax) - &
            bsupuijcf0(:,:,nsmin:nsmax)*ksupvijcf(:,:,nsmin:nsmax)
          work2(:,:,nsmin:nsmax) =  bsupsijsf0(:,:,nsmin:nsmax)*ksupvijcf(:,:,nsmin:nsmax) - &
            bsupvijcf0(:,:,nsmin:nsmax)*ksupsijsf(:,:,nsmin:nsmax)
          work3(:,:,nsmin:nsmax) =  bsupuijcf0(:,:,nsmin:nsmax)*ksupsijsf(:,:,nsmin:nsmax) - &
            bsupsijsf0(:,:,nsmin:nsmax)*ksupuijcf(:,:,nsmin:nsmax)

!      J0 X dB (Ignored for pos-def preconditioner AND linearized forces)
!      Now, bsupXij?f1 = perturbed B-field
          xfact = 1
!!        IF (l_Compute_Hessian) xfact = xfact-xPosDef
          xfact = xfact-xPosDef

          work1(:,:,nsmin:nsmax) = work1(:,:,nsmin:nsmax) + xfact*                &
                     (bsupvijcf1(:,:,nsmin:nsmax)*ksupuijcf0(:,:,nsmin:nsmax) -   &
                      bsupuijcf1(:,:,nsmin:nsmax)*ksupvijcf0(:,:,nsmin:nsmax))
          work2(:,:,nsmin:nsmax) = work2(:,:,nsmin:nsmax) + xfact*                &
                     (bsupsijsf1(:,:,nsmin:nsmax)*ksupvijcf0(:,:,nsmin:nsmax) -   &
                      bsupvijcf1(:,:,nsmin:nsmax)*ksupsijsf0(:,:,nsmin:nsmax))
          work3(:,:,nsmin:nsmax) = work3(:,:,nsmin:nsmax) + xfact*                &
                     (bsupuijcf1(:,:,nsmin:nsmax)*ksupsijsf0(:,:,nsmin:nsmax) -   &
                      bsupsijsf1(:,:,nsmin:nsmax)*ksupuijcf0(:,:,nsmin:nsmax))

       ELSE

          work1(:,:,nsmin:nsmax) = bsupvijcf1(:,:,nsmin:nsmax)*ksupuijcf(:,:,nsmin:nsmax) - &
             bsupuijcf1(:,:,nsmin:nsmax)*ksupvijcf(:,:,nsmin:nsmax)           ! WORK1 has cosine parity (IPARITY = 0)
          work2(:,:,nsmin:nsmax) = bsupsijsf1(:,:,nsmin:nsmax)*ksupvijcf(:,:,nsmin:nsmax) - &
             bsupvijcf1(:,:,nsmin:nsmax)*ksupsijsf(:,:,nsmin:nsmax)           ! WORK2 has sine parity (IPARITY = 1)
          work3(:,:,nsmin:nsmax) = bsupuijcf1(:,:,nsmin:nsmax)*ksupsijsf(:,:,nsmin:nsmax) - &
             bsupsijsf1(:,:,nsmin:nsmax)*ksupuijcf(:,:,nsmin:nsmax)           ! WORK3 has sine parity (IPARITY = 1)

       END IF

       DEALLOCATE(ksupsijsf, ksupuijcf, ksupvijcf,                      &
                  bsupsijsf1,bsupuijcf1,bsupvijcf1, stat=istat)

       ALLOCATE(work4(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
                work5(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
                work6(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'
       
       CALL tomnsp_par(work1, work4, pcos, ifull, 1, ns_span)             ! Harmonics corresponding to..
       CALL tomnsp_par(work2, work5, psin, ifull, 1, ns_span)             ! Lorentz Force on full grid
       CALL tomnsp_par(work3, work6, psin, ifull, 1, ns_span)

       DEALLOCATE(bsupuijch, bsupvijch, stat=istat)                       ! aliased work2-3 deallocated now

       nsmin=MAX(1,startglobrow); nsmax=MIN(endglobrow+1,ns)
       ns_span = nsmax-nsmin+1
!
!      Calculate radial derivative of pressure gradient 
!
       ALLOCATE(pmncf_ds(0:mpol,-ntor:ntor,nsmin:nsmax),                &
                pmncf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_FORCE'

       CALL toijsp_par(jpmnch_i, work1, 0, 0, pcos, ihalf, 1, ns_span)                    ! Calculate jac*pressure in real space (half mesh)      

       DEALLOCATE (jpmnch_i, stat=istat)
!
!      UPDATE THERMAL (PRESSURE) ENERGY
!
       IF (l_getwmhd) THEN
          IF (INHESSIAN) STOP 'l_getwmhd must be set to FALSE in Hessian'
          nsmax=MIN(endglobrow,ns)
          wp = 0
          DO l = 1, nu_i
             wp = wp + SUM(work1(l,:,MAX(2,nsmin):nsmax))*cosmui(l,m0)
          END DO
          wp0 = 2*signjac*(twopi*pi * hs_i)*wp                            !Factor of 2 if 0-pi only in theta
          nsmax=MIN(endglobrow+1,ns)
          CALL MPI_ALLREDUCE(wp0,wp,1,MPI_REAL8, MPI_SUM,               &
                             MPI_COMM_WORLD,MPI_ERR)
       END IF

!SPH   NOW CONVERT THIS TO REAL PRESSURE
!SPH   Evolved pressure is jacobh * real-pressure
       work1 = work1/jacobh_p(:,:,nsmin:nsmax)
       ALLOCATE(pmnch(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'

       CALL tomnsp_par(work1, pmnch, pcos, ihalf, 1, ns_span)

       DEALLOCATE(bsupsijsh, stat=istat)
       
       nblock=SIZE(pmnch,1)*SIZE(pmnch,2)
       CALL to_full_mesh_par(pmnch, pmncf, nblock, nsmin, nsmax)
       CALL GradientFull_par(pmncf_ds, pmnch, nblock, nsmin, nsmax)

!
!      Calculate harmonics of net force (i.e., F = JxB - grad p) on full mesh
!      THIS IS THE TRUE NONLINEAR JXB-grad(p) FORCE when l_linearize=FALSE
!      AND THE PERTURBED FORCE WHEN l_linearize=TRUE

!SPH041408 : NEED FOR grad-p - apply 1/2 factor in get_force_harmonics to pmncf_ds(1)
!      NEED pmncf(1,m1) so limit at r0=hs_i/2 in dp/du is taken correctly
       IF (nsmin .EQ. 1) THEN
          pmncf(:,:,1) = 0
          pmncf(m0:m1,:,1) = pmnch(m0:m1,:,2)
 !         pmncf_ds(:,:,1) = 0
          pmncf_ds(m1,:,1) = 2*ohs*pmnch(m1,:,2)
       END IF

       nsmax=MIN(endglobrow,ns)
       work4 = work4 - pmncf_ds(:,:,nsmin:nsmax)

       CALL get_force_harmonics_par(pmncf, work4, work5, work6,   &
                                    fsubsmncf, fsubumnsf, fsubvmnsf)

       DEALLOCATE(work4, work5, work6, pmncf, pmncf_ds, pmnch, stat=istat)

       CALL second0(skstoff)
       time_update_force = time_update_force + (skstoff-skston)

      END SUBROUTINE update_force_par
#endif

!>  \brief Serial subroutine to update MHD force
      SUBROUTINE update_force
#if defined(SKS)
       USE descriptor_mod, ONLY: DIAGONALDONE, INHESSIAN
#endif
       IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER  :: istat, l, nblock
       INTEGER, PARAMETER :: ifull=0, ihalf=1, m1=1, m0=0, pcos=0, psin=1
       REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) ::                    &
         pmncf_ds, pmncf, pmnch, jpmnch_i,                              &
         bsupsmnsh_i, bsupumnch_i, bsupvmnch_i,                         &
         ksupsijsf, ksupuijcf, ksupvijcf,                               &
         bsupsijsf1,bsupuijcf1,bsupvijcf1,                              &
         work4, work5, work6

       REAL(rprec), POINTER, DIMENSION(:,:,:) :: work1, work2, work3
       REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:), TARGET ::            &
         bsupsijsh, bsupuijch, bsupvijch
       REAL(rprec) :: ton, toff, xfact, skston, skstoff 
!-----------------------------------------------
       CALL second0(ton)
       skston=ton
!
!     Allocate space for TOTAL values of BSUPXMN?H_i and pmnch_i
!     (they are actually updated in UPDATE_STATE, not here, so UNPERTURBED
!      values are not changed until UPDATE_STATE is called)
!     increments are calculated in UPDATE_BFIELD and UPDATE_PRES
!     
       ALLOCATE(bsupsmnsh_i(0:mpol,-ntor:ntor,ns),                      &
                bsupumnch_i(0:mpol,-ntor:ntor,ns),                      &
                bsupvmnch_i(0:mpol,-ntor:ntor,ns),                      &  
                jpmnch_i   (0:mpol,-ntor:ntor,ns), stat=istat)  
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'
!
!      IF nonlinear force, add perturbation to DC (initial) value
!      IF linear force, compute JUST perturbation here (add gc0 back in funct_island)
!
       xfact = 1
       IF (l_linearize) xfact = 0
       bsupsmnsh_i = xfact*jbsupsmnsh + djbsupsmnsh
       bsupumnch_i = xfact*jbsupumnch + djbsupumnch
       bsupvmnch_i = xfact*jbsupvmnch + djbsupvmnch
       jpmnch_i    = xfact*jpmnch     + djpmnch

!SPH:  Convert to real covariant components

       ALLOCATE(bsupsijsh(ntheta,nzeta,ns),                             &
                bsupuijch(ntheta,nzeta,ns),                             &
                bsupvijch(ntheta,nzeta,ns), stat=istat)  
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'

!NEED VALUES AT js=1 (STORED VALUES AT ORIGIN IN UPDATE_BFIELD)

       CALL toijsp(bsupsmnsh_i, bsupsijsh, 0, 0, psin, ifull)             ! Calculate BSUP on half mesh
       CALL toijsp(bsupumnch_i, bsupuijch, 0, 0, pcos, ifull)
       CALL toijsp(bsupvmnch_i, bsupvijch, 0, 0, pcos, ifull)

       bsupsijsh = bsupsijsh/jacobh
!       bsupuijch = bsupuijch/jacobh
       bsupvijch = bsupvijch/jacobh

       DEALLOCATE(bsupsmnsh_i, bsupumnch_i, bsupvmnch_i)
!
!      Calculate contravariant components of J = curl B (Ksup = gsqrt*J)
!      and update magnetic energy if nonlinear force is being computed
!                             
       ALLOCATE(ksupsijsf(ntheta,nzeta,ns),                             & ! Full mesh quantities (real space)
                ksupuijcf(ntheta,nzeta,ns),                             & ! K_s (sine), K_u (cosine), K_v (cosine)
                ksupvijcf(ntheta,nzeta,ns),                             &
                bsupsijsf1(ntheta,nzeta,ns),                            &
                bsupuijcf1(ntheta,nzeta,ns),                            &
                bsupvijcf1(ntheta,nzeta,ns), stat=istat)
       IF (istat.ne.0) STOP 'Allocation failed in UPDATE_FORCE'

!
!      Calculate updated FULL-MESH components of BSUPX in real-space
!
       CALL bhalftobfull (bsupsijsh,  bsupuijch,  bsupvijch,            &
                          bsupsijsf1, bsupuijcf1, bsupvijcf1)
       IF (ANY(bsupsijsf1(:,:,ns) .NE. zero)) PRINT *,'bsupsijsf1(ns) != 0'

!Origin currents are X1/2 (1/2 factor in grad-p too, in update_pres)
       CALL cv_currents (bsupsijsh, bsupuijch, bsupvijch,               &
                         ksupsijsf, ksupuijcf, ksupvijcf,               &
                         .NOT.l_linearize, .FALSE.)

!      alias work arrays for use as temp work space
       work1=>bsupsijsh; work2=>bsupuijch;  work3=>bsupvijch

       IF (l_linearize) THEN
!      dJ X B0 (note here, ksupXij?f is the perturbed current)
          work1 =  bsupvijcf0*ksupuijcf - bsupuijcf0*ksupvijcf
          work2 =  bsupsijsf0*ksupvijcf - bsupvijcf0*ksupsijsf
          work3 =  bsupuijcf0*ksupsijsf - bsupsijsf0*ksupuijcf

!      J0 X dB (Ignored for pos-def preconditioner AND linearized forces)
!      Now, bsupXij?f1 = perturbed B-field
          xfact = 1
!!        IF (l_Compute_Hessian) xfact = xfact-xPosDef
          xfact = xfact-xPosDef

          work1 = work1 + xfact*                                        &
      &          (bsupvijcf1*ksupuijcf0 - bsupuijcf1*ksupvijcf0)
          work2 = work2 + xfact*                                        &
      &          (bsupsijsf1*ksupvijcf0 - bsupvijcf1*ksupsijsf0)
          work3 = work3 + xfact*                                        &
      &          (bsupuijcf1*ksupsijsf0 - bsupsijsf1*ksupuijcf0)

       ELSE

          work1 = bsupvijcf1*ksupuijcf - bsupuijcf1*ksupvijcf           ! WORK1 has cosine parity (IPARITY = 0)
          work2 = bsupsijsf1*ksupvijcf - bsupvijcf1*ksupsijsf           ! WORK2 has sine parity (IPARITY = 1)
          work3 = bsupuijcf1*ksupsijsf - bsupsijsf1*ksupuijcf           ! WORK3 has sine parity (IPARITY = 1)

       END IF
      
       DEALLOCATE(ksupsijsf, ksupuijcf, ksupvijcf,                      &
                  bsupsijsf1,bsupuijcf1,bsupvijcf1)

       ALLOCATE(work4(0:mpol,-ntor:ntor,ns),                            &
                work5(0:mpol,-ntor:ntor,ns),                            &
                work6(0:mpol,-ntor:ntor,ns), stat = istat)  
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'
       
       CALL tomnsp(work1, work4, pcos, ifull)                             ! Harmonics corresponding to..
       CALL tomnsp(work2, work5, psin, ifull)                             ! Lorentz Force on full grid
       CALL tomnsp(work3, work6, psin, ifull)
       DEALLOCATE(bsupuijch, bsupvijch)                                   ! aliased work2-3 deallocated now

!
!      Calculate radial derivative of pressure gradient 
!
       ALLOCATE(pmncf_ds(0:mpol,-ntor:ntor,ns),                         &
                pmncf(0:mpol,-ntor:ntor,ns), stat=istat)
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'
       CALL toijsp(jpmnch_i, work1, 0, 0, pcos, ihalf)                    ! Calculate jac*pressure in real space (half mesh)      
       DEALLOCATE (jpmnch_i)

       IF (.not.l_linearize) THEN
!
!      UPDATE THERMAL (PRESSURE) ENERGY
!
          wp = 0
          DO l = 1, nu_i
             wp = wp + SUM(work1(l,:,2:))*cosmui(l,m0)
          END DO

          wp = 2*signjac*(twopi*pi * hs_i)*wp                            !Factor of 2 if 0-pi only in theta

       END IF

!SPH   NOW CONVERT THIS TO REAL PRESSURE
!SPH   Evolved pressure is jacobh * real-pressure
       work1 = work1/jacobh
       ALLOCATE(pmnch(0:mpol,-ntor:ntor,ns), stat=istat)  
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_FORCE'
       CALL tomnsp(work1, pmnch, pcos, ihalf)
       DEALLOCATE(bsupsijsh)

       nblock=SIZE(pmnch,1)*SIZE(pmnch,2)
       CALL to_full_mesh(pmnch, pmncf, nblock)
       CALL GradientFull(pmncf_ds, pmnch, nblock)
!
!      Calculate harmonics of net force (i.e., F = JxB - grad p) on full mesh
!      THIS IS THE TRUE NONLINEAR JXB-grad(p) FORCE when l_linearize=FALSE
!      AND THE PERTURBED FORCE WHEN l_linearize=TRUE

!SPH041408 : NEED FOR grad-p - apply 1/2 factor in get_force_harmonics to pmncf_ds(1)
!      NEED pmncf(1,m1) so limit at r0=hs_i/2 in dp/du is taken correctly
       pmncf(:,:,1) = 0
       pmncf(m0:m1,:,1) = pmnch(m0:m1,:,2)
       pmncf_ds(:,:,1) = 0
       pmncf_ds(m1,:,1) = 2*ohs*pmnch(m1,:,2)

       CALL get_force_harmonics(pmncf, pmncf_ds, work4, work5, work6,   &
                                fsubsmncf, fsubumnsf, fsubvmnsf)

       DEALLOCATE(work4, work5, work6, pmncf, pmncf_ds, pmnch)

       CALL second0(toff)
       time_update_force = time_update_force + (toff-ton)

       CALL second0(skstoff)
#if defined(SKS)
       IF(DIAGONALDONE.AND.INHESSIAN) time_update_force=time_update_force+(skstoff-skston)
#endif
      END SUBROUTINE update_force

      END MODULE siesta_force
