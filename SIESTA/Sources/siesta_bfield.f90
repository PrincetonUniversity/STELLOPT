!>  \brief Module contained subroutines for updating magnetic fields as part of the SIESTA project
!!  \author S. P. Hirshman and R. Sanchez
!!  \date Aug 21, 2006
       MODULE siesta_bfield

       USE stel_kinds
       USE stel_constants
       USE fourier 
       USE quantities
       USE island_params, ONLY: nu_i
       USE perturbation, ONLY: lresistive, eta_factor, buv_res
       USE timer_mod
       USE shared_data, ONLY: delta_t, fsq_total, l_natural, niter, fsq_res

       CONTAINS
#if defined(SKS)      
!>  \brief Parallel routine for updating B^s, B^u, B^v for ideal (and resistive) perturbations
       SUBROUTINE update_bfield_par (l_add_res)
!     
!     PURPOSE: ADVANCES MAGNETIC FIELD FROM T TO T+DELTA_T
!              Updates values of BSUPMNF and BSUPIJF in QUANTITIES MODULE
!
       USE nscalingtools, ONLY: startglobrow, endglobrow 
       IMPLICIT NONE
       INCLUDE 'mpif.h'
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       LOGICAL, INTENT(in)     :: l_add_res
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER :: m0=0, m1=1, m2=2, pcos=0, psin=1, ifull=0
       INTEGER :: m, n, l, istat, js, nblock
       LOGICAL :: ladd_pert
       REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE ::                    &
         esubsijsf, esubuijcf, esubvijcf,                               &
         esubsmnsf, esubumncf, esubvmncf, esubumnch, esubvmnch,         &
         esubsmnsh, work1, work2, resistivity
       REAL(rprec) :: ton, toff, rho
       REAL(rprec) :: delt_cfl, eta_prof, r0
       INTEGER :: nsmin, nsmax, ns_span
 !-----------------------------------------------
       CALL second0(ton)

       nsmin=MAX(1,startglobrow-1); nsmax=MIN(ns,endglobrow+1)
       ns_span = nsmax-nsmin+1
       ALLOCATE(esubsijsf(ntheta,nzeta,nsmin:nsmax),                    &
                esubuijcf(ntheta,nzeta,nsmin:nsmax),                    &
                esubvijcf(ntheta,nzeta,nsmin:nsmax), stat=istat)   
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_BFIELD' 
!
!      Calculate ideal covariant components of electric field E-ideal = -v X B on full-space mesh
!  
       esubsijsf = -(jvsupuijsf(:,:,nsmin:nsmax)*bsupvijcf0(:,:,nsmin:nsmax)    &
                   - jvsupvijsf(:,:,nsmin:nsmax)*bsupuijcf0(:,:,nsmin:nsmax))    ! -E_s = V^uB^v - V^vB^u  (SINE): jv -> jac*v
       esubuijcf = -(jvsupvijsf(:,:,nsmin:nsmax)*bsupsijsf0(:,:,nsmin:nsmax)    &
                    -jvsupsijcf(:,:,nsmin:nsmax)*bsupvijcf0(:,:,nsmin:nsmax))    ! -E_u = V^vB^s - V^sB^v  (COSINE)
       esubvijcf = -(jvsupsijcf(:,:,nsmin:nsmax)*bsupuijcf0(:,:,nsmin:nsmax)    &
                   - jvsupuijsf(:,:,nsmin:nsmax)*bsupsijsf0(:,:,nsmin:nsmax))    ! -E_v = V^sB^u - V^uB^s  (COSINE)

!      confirm boundary condition: esubu,v(s=1) = 0 (tangential E vanishes at bdy)
       IF (nsmax .EQ. ns) THEN
         IF (ANY(esubuijcf(:,:,ns) .NE. zero)) THEN
           STOP 'esubuijcf(ns) != 0 in UPDATE_BFIELD_par'
         ELSE IF (ANY(esubvijcf(:,:,ns) .NE. zero)) THEN
           STOP 'esubvijcf(ns) != 0 in UPDATE_BFIELD_par'
         END IF
       END IF

!
!      Note this will LOWER the energy due to eta*|J|||**2 heating
!      esubX(resistive) = eta(JdotB/B^2)*BsubX
!      so the magnetic energy DECREASES due to this term. Note ksubX=jac*JsubX are the
!      covariant components of the current TIMES the jac factor
!
       ladd_pert = ALLOCATED(buv_res)
       IF (lresistive .AND. (l_add_res .or. ladd_pert)) THEN
          delt_cfl = hs_i**2*ABS(eta_factor)
          IF (fsq_total < fsq_res) delt_cfl = delt_cfl*SQRT(fsq_total/fsq_res)

          ALLOCATE(resistivity(ntheta,nzeta,nsmin:nsmax), stat=istat)
          IF (istat .NE. 0) STOP 'ALLOCATION ERROR IN update_bfield'

          DO js = nsmin, nsmax
             rho = hs_i*(js-1)
             eta_prof = rho*rho*(1-rho)
             resistivity(:,:,js) = delt_cfl*eta_prof
          END DO

!Note: need to divide current by 1/jac
          resistivity = resistivity/jacobf_p(:,:,nsmin:nsmax)
          IF (ladd_pert) THEN
             DO js = nsmin, nsmax
                resistivity(:,:,js) = resistivity(:,:,js)*buv_res(:,:,js)
             END DO
          END IF

!Isotropic resistivity, E ~ eta*J
          DO js = nsmin, nsmax
             esubsijsf(:,:,js) = esubsijsf(:,:,js) + resistivity(:,:,js)*ksubsijsf(:,:,js)
             esubuijcf(:,:,js) = esubuijcf(:,:,js) + resistivity(:,:,js)*ksubuijcf(:,:,js)
             esubvijcf(:,:,js) = esubvijcf(:,:,js) + resistivity(:,:,js)*ksubvijcf(:,:,js)
          END DO

          DEALLOCATE (resistivity)
       ENDIF 

!
!   Calculate harmonics of electric field
! 
       ALLOCATE(esubumncf(0:mpol,-ntor:ntor,nsmin:nsmax),               & 
                esubvmncf(0:mpol,-ntor:ntor,nsmin:nsmax),               &
                esubsmnsf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_BFIELD'

       CALL tomnsp_par(esubsijsf, esubsmnsf, psin, ifull, 1, ns_span)
       CALL tomnsp_par(esubuijcf, esubumncf, pcos, ifull, 1, ns_span)
       CALL tomnsp_par(esubvijcf, esubvmncf, pcos, ifull, 1, ns_span)

       DEALLOCATE(esubsijsf, esubuijcf, esubvijcf, stat=istat)
       
!      CANCEL ROUND-OFF ERROR BUT KEEP m=1 COMPONENT OF esubu(1)
       IF (nsmin .EQ. 1) THEN
          esubsmnsf(m0,:,1) = 0;  esubsmnsf(m2:,:,1) = 0
          esubumncf(m0,:,1) = 0;  esubumncf(m2:,:,1) = 0
          esubvmncf(m1:,:,1)= 0;
       END IF

       ALLOCATE(esubsmnsh(0:mpol,-ntor:ntor,nsmin:nsmax),               &
                esubumnch(0:mpol,-ntor:ntor,nsmin:nsmax),               &
                esubvmnch(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_BFIELD'  
       
       !!On exit, esubXh valid on nsmin+1:nsmax grid
       nblock = SIZE(esubsmnsh,1)*SIZE(esubsmnsh,2)
       CALL to_half_mesh_par(esubsmnsf, esubsmnsh, nblock, nsmin, nsmax)
       CALL to_half_mesh_par(esubumncf, esubumnch, nblock, nsmin, nsmax)
       CALL to_half_mesh_par(esubvmncf, esubvmnch, nblock, nsmin, nsmax)

       ALLOCATE(work1(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
                work2(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_BFIELD'

!      Gradients of esubu,v on half grid j-half=2,ns: note esubu,v(1) used here!
       !!On exit, work1,2 valid on nsmin+1:nsmax grid
       CALL GradientHalf_par(work1, esubumncf, nblock, nsmin, nsmax)
       CALL GradientHalf_par(work2, esubvmncf, nblock, nsmin, nsmax)

!      Calculate r.h.s. of the evolution equations of all components of jacobh*BSUP?
!      Compute contravariant components of -curl(E) for right side of Faraday's Equation
!

       r0 = hs_i/2
       nsmin=MAX(1,startglobrow)

!       IF (nsmin .EQ. 1) THEN
!!!!!!!IMPOSE NATURAL BOUNDARY CONDITION AT FIRST HALF-GRID POINT
!          djbsupsmnsh(m1,:,1) = (esubsmnsh(m1,:,2)*r0-esubumnch(m1,:,2))/2
!          esubsmnsh(m1,:,2) = esubsmnsh(m1,:,2) - djbsupsmnsh(m1,:,1)/r0
!          esubumnch(m1,:,2) = esubumnch(m1,:,2) + djbsupsmnsh(m1,:,1)
!       END IF

       DO m = 0, mpol
         DO n = -ntor, ntor
           djbsupsmnsh(m,n,nsmin:nsmax) =      m*esubvmnch(m,n,nsmin:nsmax)               & ! -m and -n*NFP from deriv. COS(m*theta+n*NFP*phi)
                                          -n*nfp*esubumnch(m,n,nsmin:nsmax)
           djbsupumnch(m,n,nsmin:nsmax) = -n*nfp*esubsmnsh(m,n,nsmin:nsmax) + work2(m,n,nsmin:nsmax) ! n*NFP from deriv. SIN(m*theta+n*NFP*phi)             
           djbsupvmnch(m,n,nsmin:nsmax) =      m*esubsmnsh(m,n,nsmin:nsmax) - work1(m,n,nsmin:nsmax) ! m from deriv. SIN(m*theta+n*NFP*phi) 
          ENDDO
       ENDDO

!CONSTRAINT AT ORIGIN [jbsups(m=1) = r0*jbsupu(m=1) <=> esubumnch = r0*esubsmnsh]
!Origin bc for m=1: e_u = r0*e_s (Follows form v X B at origin)
!This suffices to make B^s = r0*B^u
       IF (nsmin .EQ. 1) THEN

          djbsupsmnsh(:,:,1) = 0
          djbsupumnch(:,:,1) = 0
          djbsupvmnch(:,:,1) = 0
      
!      SYMMETRIC FORM
          djbsupumnch(m1,:,1) = djbsupumnch(m1,:,2)
          IF (l_natural) THEN
             djbsupsmnsh(m1,:,1) = r0*djbsupumnch(m1,:,1)
          ELSE
             djbsupsmnsh(m1,:,1) = djbsupsmnsh(m1,:,2)
          END IF

          djbsupvmnch(m0,:,1) = djbsupvmnch(m0,:,2)
 
       END IF
 
       DEALLOCATE(esubumncf, esubvmncf, esubsmnsf)
       DEALLOCATE(work1, work2, esubsmnsh, esubumnch, esubvmnch)    
!
!      Calculate INCREMENT of magnetic field harmonics: use Euler scheme with DT given by VSUBMN advance equations.
!      Results are DBSUP?IJ?F, which is used in the calculation of the force and to advance the B's in UPDATE_STATE
!
!      B^S (SIN) 
       djbsupsmnsh(:,:,nsmin:nsmax) = delta_t*djbsupsmnsh(:,:,nsmin:nsmax)
!
!      B^U (COS) 
       djbsupumnch(:,:,nsmin:nsmax) = delta_t*djbsupumnch(:,:,nsmin:nsmax)
!
!      B^V (COS) 
       djbsupvmnch(:,:,nsmin:nsmax) = delta_t*djbsupvmnch(:,:,nsmin:nsmax)

       CALL second0(toff)
       time_update_bfield = time_update_bfield + (toff-ton)

       END SUBROUTINE update_bfield_par
#endif

!>  \brief Serial routine for updating B^s, B^u, B^v for ideal (and resistive) perturbations
       SUBROUTINE update_bfield (l_add_res)
!     
!     PURPOSE: ADVANCES MAGNETIC FIELD FROM T TO T+DELTA_T
!              Updates values of BSUPMNF and BSUPIJF in QUANTITIES MODULE
!
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       LOGICAL, INTENT(in)     :: l_add_res
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER :: m0=0, m1=1, m2=2, pcos=0, psin=1, ifull=0
       INTEGER :: m, n, l, istat, js, nblock
       LOGICAL :: ladd_pert
       REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE ::                    &
         esubsijsf, esubuijcf, esubvijcf,                               &
         esubsmnsf, esubumncf, esubvmncf, esubumnch, esubvmnch,         &
         esubsmnsh, work1, work2, resistivity
       REAL(rprec) :: ton, toff, rho 
       REAL(rprec) :: delt_cfl, eta_prof, r0
!-----------------------------------------------
       CALL second0(ton)

       ALLOCATE(esubsijsf(ntheta,nzeta,ns),                             &
                esubuijcf(ntheta,nzeta,ns),                             &
                esubvijcf(ntheta,nzeta,ns), stat=istat)   
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_BFIELD' 

!
!      Calculate ideal covariant components of electric field E-ideal = -v X B on full-space mesh
!  
       esubsijsf = -(jvsupuijsf*bsupvijcf0 - jvsupvijsf*bsupuijcf0)        ! -E_s = V^uB^v - V^vB^u  (SINE): jv -> jac*v
       esubuijcf = -(jvsupvijsf*bsupsijsf0 - jvsupsijcf*bsupvijcf0)        ! -E_u = V^vB^s - V^sB^v  (COSINE)
       esubvijcf = -(jvsupsijcf*bsupuijcf0 - jvsupuijsf*bsupsijsf0)        ! -E_v = V^sB^u - V^uB^s  (COSINE)

#if defined(EPAR)
!
!      Add resistive components E-resistive = eparB B: why doesn't it work if I set 0*jeparB?
!
       esubsijsf = esubsijsf + jeparBijcf*bsubsijsf
       esubuijcf = esubuijcf + jeparBijcf*bsubuijcf
       esubvijcf = esubvijcf + jeparBijcf*bsubvijcf
#endif

!      confirm boundary condition: esubu,v(s=1) = 0 (tangential E vanishes at bdy)
       IF (ANY(esubuijcf(:,:,ns) .NE. zero)) THEN
          STOP 'esubuijcf(ns) != 0 in UPDATE_BFIELD'
       ELSE IF (ANY(esubvijcf(:,:,ns) .NE. zero)) THEN
          STOP 'esubvijcf(ns) != 0 in UPDATE_BFIELD'
       END IF
!
!      Note this will LOWER the energy due to eta*|J|||**2 heating
!      esubX(resistive) = eta(JdotB/B^2)*BsubX
!      so the magnetic energy DECREASES due to this term. Note ksubX=jac*JsubX are the
!      covariant components of the current TIMES the jac factor
!
       ladd_pert = ALLOCATED(buv_res)
       IF (lresistive .AND. (l_add_res .OR. ladd_pert)) THEN
          delt_cfl = hs_i**2*ABS(eta_factor)
          IF (fsq_total < fsq_res) delt_cfl = delt_cfl*SQRT(fsq_total/fsq_res)

          ALLOCATE(resistivity(ntheta,nzeta,ns), stat=istat)
          IF (istat .ne. 0) STOP 'ALLOCATION ERROR IN update_bfield'

          DO js = 1, ns
             rho = hs_i*(js-1)
             eta_prof = rho*rho*(1-rho)
             resistivity(:,:,js) = delt_cfl*eta_prof
          END DO

!Note: need to divide current through by 1/jac
          resistivity = resistivity/jacobf
          IF (ladd_pert) resistivity = resistivity*buv_res

!Isotropic resistivity, E ~ eta*J
          esubsijsf = esubsijsf + resistivity*ksubsijsf
          esubuijcf = esubuijcf + resistivity*ksubuijcf
          esubvijcf = esubvijcf + resistivity*ksubvijcf

          DEALLOCATE (resistivity)
       ENDIF 
!
!   Calculate harmonics of electric field
! 
       ALLOCATE(esubumncf(0:mpol,-ntor:ntor,ns),                        & 
                esubvmncf(0:mpol,-ntor:ntor,ns),                        &
                esubsmnsf(0:mpol,-ntor:ntor,ns), stat=istat)  
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_BFIELD'

       CALL tomnsp(esubsijsf, esubsmnsf, psin, ifull)
       CALL tomnsp(esubuijcf, esubumncf, pcos, ifull)
       CALL tomnsp(esubvijcf, esubvmncf, pcos, ifull)
       DEALLOCATE(esubsijsf, esubuijcf, esubvijcf, stat=istat)
       
!      CANCEL ROUND-OFF ERROR BUT KEEP m=1 COMPONENT OF esubu(1)
       esubsmnsf(m0,:,1) = 0;  esubsmnsf(m2:,:,1) = 0
       esubumncf(m0,:,1) = 0;  esubumncf(m2:,:,1) = 0
       esubvmncf(m1:,:,1)= 0

       ALLOCATE(esubsmnsh(0:mpol,-ntor:ntor,ns),                        &
                esubumnch(0:mpol,-ntor:ntor,ns),                        &
                esubvmnch(0:mpol,-ntor:ntor,ns), stat=istat)  
       IF (istat .ne. 0) STOP 'Allocation failed in UPDATE_BFIELD'  

       nblock = SIZE(esubsmnsf,1)*SIZE(esubsmnsf,2)
       CALL to_half_mesh(esubsmnsf, esubsmnsh, nblock)
       CALL to_half_mesh(esubumncf, esubumnch, nblock)
       CALL to_half_mesh(esubvmncf, esubvmnch, nblock)

       ALLOCATE(work1(0:mpol,-ntor:ntor,ns),                            &
                work2(0:mpol,-ntor:ntor,ns), stat=istat)  
       IF (istat .NE. 0) STOP 'Allocation failed in UPDATE_BFIELD'

!      Gradients of esubu,v on half grid j-half=2,ns: note esubu,v(1) used here!
       CALL GradientHalf(work1, esubumncf, nblock)
       CALL GradientHalf(work2, esubvmncf, nblock)
!
!      Calculate r.h.s. of the evolution equations of all components of jacobh*BSUP?
!      Compute contravariant components of -curl(E) for right side of Faraday's Equation
!
       r0 = hs_i/2
!!!!!!!IMPOSE NATURAL BOUNDARY CONDITION AT FIRST HALF-GRID POINT
!       djbsupsmnsh(m1,:,1) = (esubsmnsh(m1,:,2)*r0-esubumnch(m1,:,2))/2
!       esubsmnsh(m1,:,2) = esubsmnsh(m1,:,2) - djbsupsmnsh(m1,:,1)/r0
!       esubumnch(m1,:,2) = esubumnch(m1,:,2) + djbsupsmnsh(m1,:,1)

       DO m = 0, mpol
         DO n = -ntor, ntor
           djbsupsmnsh(m,n,2:) =      m*esubvmnch(m,n,2:)               & ! -m and -n*NFP from deriv. COS(m*theta+n*NFP*phi)
                                 -n*nfp*esubumnch(m,n,2:)
           djbsupumnch(m,n,2:) = -n*nfp*esubsmnsh(m,n,2:) + work2(m,n,2:) ! n*NFP from deriv. SIN(m*theta+n*NFP*phi)             
           djbsupvmnch(m,n,2:) =      m*esubsmnsh(m,n,2:) - work1(m,n,2:) ! m from deriv. SIN(m*theta+n*NFP*phi) 
          ENDDO
       ENDDO

!CONSTRAINT AT ORIGIN [jbsups(m=1) = r0*jbsupu(m=1) <=> esubumnch = r0*esubsmnsh]
!Origin bc for m=1: e_u = r0*e_s (Follows form v X B at origin)
!This suffices to make B^s = r0*B^u
       djbsupsmnsh(:,:,1) = 0
       djbsupumnch(:,:,1) = 0
       djbsupvmnch(:,:,1) = 0

!      SYMMETRIC FORM
       djbsupumnch(m1,:,1) = djbsupumnch(m1,:,2)
       IF (l_natural) THEN
          djbsupsmnsh(m1,:,1) = r0*djbsupumnch(m1,:,1)
       ELSE
          djbsupsmnsh(m1,:,1) = djbsupsmnsh(m1,:,2)
       END IF

       djbsupvmnch(m0,:,1) = djbsupvmnch(m0,:,2)

       DEALLOCATE(esubumncf, esubvmncf, esubsmnsf)
       DEALLOCATE(work1, work2, esubsmnsh, esubumnch, esubvmnch)    
!
!      Calculate INCREMENT of magnetic field harmonics: use Euler scheme with DT given by VSUBMN advance equations.
!      Results are DBSUP?IJ?F, which is used in the calculation of the force and to advance the B's in UPDATE_STATE
!
!      B^S (SIN) 
!
       djbsupsmnsh = delta_t*djbsupsmnsh
!
!      B^U (COS) 
!    
       djbsupumnch = delta_t*djbsupumnch
!
!      B^V (COS) 
!       
       djbsupvmnch = delta_t*djbsupvmnch

       CALL second0(toff)
       time_update_bfield = time_update_bfield + (toff-ton)

       END SUBROUTINE update_bfield

       END MODULE siesta_bfield
