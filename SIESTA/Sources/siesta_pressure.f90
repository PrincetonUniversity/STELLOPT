!>  \brief Module contained subroutines for updating pressure as part of the SIESTA project
!!  \author S. P. Hirshman and R. Sanchez
!!  \date Aug 18, 2006
      MODULE siesta_pressure

       USE stel_kinds
       USE stel_constants
       USE fourier       
       USE quantities
       USE hessian, ONLY: l_Compute_Hessian
       USE shared_data, ONLY: l_linearize, l_push_s, l_push_u, xPosDef
       USE timer_mod, ONLY: time_update_pres

       CONTAINS
#if defined(SKS) 
!>     \brief Parallel subroutine to advance pressure from t to t+delta_t by updating jpmnch (jacobh*presh)
       SUBROUTINE update_pres_par (delta_t)
       USE nscalingtools, ONLY: startglobrow, endglobrow
       IMPLICIT NONE
       INCLUDE 'mpif.h'
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(rprec), INTENT(in)  :: delta_t
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER     :: ifull=0, ihalf=1, pcos=0, psin=1, m0=0, m1=1
       LOGICAL, PARAMETER     :: l_upwind=.FALSE.
       INTEGER                :: m, n, istat, js, nblock
       REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE ::                     &
         work1, work2, work3, work4, work5, work6,                       &
         jvsupuijsh, jvsupvijsh
       REAL(rprec)            :: ton, toff, gamm1
       INTEGER                :: nsmin, nsmax, ns_span
!-----------------------------------------------
       CALL second0(ton)
       nsmin=MAX(1,startglobrow-1); nsmax=MIN(ns,endglobrow+1)
       ns_span=nsmax-nsmin+1
!
!      First, calculate first term (~gamma) on r.h.s. of pressure evolution equation       
! 
       ALLOCATE(work1(ntheta,nzeta,nsmin:nsmax),                        & ! Allocate work space
                work2(ntheta,nzeta,nsmin:nsmax),                                 &   
                work3(ntheta,nzeta,nsmin:nsmax), stat=istat)   
       IF (istat .ne. 0) STOP 'Allocation1 failed in UPDATE_PRES'

!      Div(p*Flux) terms
!      Note that near axis, jvsups ~ s cos(u+nv), jvsupv -> s sin(nv) due to jacobian element
!      In particular, pijcf(1)*jvsupsijcf makes no contribution to forces (integrates to zero,
!      since jvsupsijcf(1) ~ cos(u) + O(s)!)
!      jvsupu may be finite, although dp/du -> 0 at axis
!      Note: WORK1 has COSINE parity, nsh = ns_i-1


!Upwind to suppress grid separation in radial direction
       IF (l_upwind) THEN
          work1(:,:,nsmin:nsmax) = -0.5_dp*signjac*hs_i*                &
             pijcf_ds(:,:,nsmin:nsmax)*ABS(jvsupsijcf(:,:,nsmin:nsmax))
          work2(:,:,nsmin:nsmax) = work1(:,:,nsmin:nsmax) +             &
             pijch(:,:,nsmin:nsmax)*jvsupsijcf(:,:,nsmin:nsmax) 
       ELSE
          work2(:,:,nsmin:nsmax) = pijcf(:,:,nsmin:nsmax)*jvsupsijcf(:,:,nsmin:nsmax)
       END IF

       !!On exit, work1 valid on nsmin+1:nsmax grid
       nblock = SIZE(work1,1)*SIZE(work1,2)
       CALL GradientHalf_par(work1, work2, nblock, nsmin, nsmax)

       ALLOCATE (jvsupuijsh(ntheta,nzeta,nsmin:nsmax),                  &   
                 jvsupvijsh(ntheta,nzeta,nsmin:nsmax), stat=istat)
       IF (istat .ne. 0) STOP 'Allocation1 failed in UPDATE_PRES'

       !!On exit, jsupXh valid on nsmin+1:nsmax grid
       CALL to_half_mesh_par(jvsupuijsf(:,:,nsmin:nsmax), jvsupuijsh, nblock, nsmin, nsmax)
       CALL to_half_mesh_par(jvsupvijsf(:,:,nsmin:nsmax), jvsupvijsh, nblock, nsmin, nsmax)

       ALLOCATE(work4(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
                work5(0:mpol,-ntor:ntor,nsmin:nsmax),                   &
                work6(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
       IF (istat .ne. 0) STOP 'Allocation2 failed in UPDATE_PRES'

       nsmin=nsmin+1   !MAX(1,startglobrow)
       work2(:,:,nsmin:nsmax) = pijch(:,:,nsmin:nsmax)*jvsupuijsh(:,:,nsmin:nsmax)
       work3(:,:,nsmin:nsmax) = pijch(:,:,nsmin:nsmax)*jvsupvijsh(:,:,nsmin:nsmax)
       
!      FFT of Div(p*Flux) terms (multiply below by m,n*nfp)
       CALL tomnsp_par(work2, work5, psin, ihalf, 2, ns_span)
       CALL tomnsp_par(work3, work6, psin, ihalf, 2, ns_span)
       CALL tomnsp_par(work1, work4, pcos, ihalf, 2, ns_span)
              
       DO m = 0, mpol
          DO n = -ntor, ntor
             work4(m,n,nsmin:nsmax) = work4(m,n,nsmin:nsmax) +          &
             m*work5(m,n,nsmin:nsmax) + n*nfp*work6(m,n,nsmin:nsmax)      !m: Pol. deriv of SIN(m*theta+n*NFP*phi)
          ENDDO                                                           !n: Tor. deriv of SIN(m*theta+n*NFP*phi)
       ENDDO

       djpmnch(:,:,nsmin:nsmax) = -gamma*work4(:,:,nsmin:nsmax)           !Gamma = adiabatic factor
       
       DEALLOCATE(work5, work6)
!   
!      Now, calculate next term ~(gamma-1) on r.h.s. of pressure equation
!      the u,v derivatives should be consistent with div(pv) term and combine
!      to cancel if div(v) = 0
       nsmin=MAX(1,startglobrow-1)
       work1(:,:,nsmin:nsmax) =                                         &
         pijcf_ds(:,:,nsmin:nsmax)*jvsupsijcf(:,:,nsmin:nsmax) +        &
         pijcf_du(:,:,nsmin:nsmax)*jvsupuijsf(:,:,nsmin:nsmax) +        &
         pijcf_dv(:,:,nsmin:nsmax)*jvsupvijsf(:,:,nsmin:nsmax)

       !!On exit, work2 valid on nsmin+1:nsmax grid
       CALL to_half_mesh_par(work1, work2, nblock, nsmin, nsmax)

!      MIXED HALF/FULL: DIAGONAL IN ANGULAR DERIVATIVES
!!       work1 = pijcf_ds*jvsupsijcf
!!       CALL to_half_mesh_par(work1,work2,nblock,nsmin,nsmax)
!!       work2 = work2+jvsupuijsh*pijch_du+jvsupvijsh*pijch_dv

       DEALLOCATE(jvsupuijsh, jvsupvijsh)

!Upwinding correction (makes 1st order streaming dominantly diagonal)
!NOTE: since it enters conservatively, it doesn't effect MHD force, only HESSIAN
!Can multiply this by a positive factor to increase numerical diffusion
       IF (l_upwind) THEN
          work1(:,:,nsmin:nsmax) = -0.5_dp*signjac*                     &
             pijcf_ds(:,:,nsmin:nsmax)*ABS(jvsupsijcf(:,:,nsmin:nsmax))
          IF (nsmin .EQ. 1) work1(:,:,1) = 0
          CALL GradientHalf_par(work3, work1, nblock, nsmin, nsmax)
          work2(:,:,nsmin:nsmax) = work2(:,:,nsmin:nsmax) +             &
                              hs_i*work3(:,:,nsmin:nsmax)
       END IF
       
       nsmin=nsmin+1
       CALL tomnsp_par(work2, work4, pcos, ihalf, 2, ns_span)            ! Calculate Fourier transform of WORK2

!      MAKES PRES-PERTURB A POSITIVE DEFINITE CONTRIBUTION TO HESSIAN: 
!      When computing hessian, guarantee positive-definiteness this way
!      (sound wave compressional energy contribution)
       gamm1 = gamma-1
!!       IF (l_Compute_Hessian) gamm1 = gamm1 + xPosDef
       IF (l_linearize) gamm1 = gamm1 + xPosDef

       djpmnch(:,:,nsmin:nsmax) = djpmnch(:,:,nsmin:nsmax) +            &
              gamm1*work4(:,:,nsmin:nsmax)                                ! Harmonics of r.h.s. of Pres. Eq. (full mesh) 

!SPH   Evolve harmonics of (jacobh * pres)mn to conserve pressure
       djpmnch(:,:,nsmin:nsmax) = delta_t*djpmnch(:,:,nsmin:nsmax)    ! calculate increment of pressure harmonics.
!SPH   Suppress radial mesh separation together with Orszag 3/2 integration rule in metrics
       IF (nsmin .EQ. 2) THEN
!          djpmnch(m1:,:,2) = 0             
          djpmnch(:,:,1)  = 0
          djpmnch(m0,:,1) = djpmnch(m0,:,2)                           ! store js=1, m=0 value
       END IF

       DEALLOCATE(work1, work2, work3, work4)

       CALL second0(toff)
       time_update_pres = time_update_pres + (toff-ton)


       END SUBROUTINE update_pres_par
#endif

!>    \brief Serial subroutine to advance pressure from t to t+delta_t by updating jpmnch (jacobh*presh)
       SUBROUTINE update_pres (delta_t)
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       REAL(rprec), INTENT(in)  :: delta_t
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER, PARAMETER     :: ifull=0, ihalf=1, pcos=0, psin=1, m0=0, m1=1
       LOGICAL, PARAMETER     :: l_upwind=.FALSE.
       INTEGER :: m, n, istat, js, nblock
       REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE ::                     &
         work1, work2, work3, work4, work5, work6,                       &
         jvsupuijsh, jvsupvijsh
       REAL(rprec)            :: ton, toff, gamm1
!-----------------------------------------------
       CALL second0(ton)
!
!      First, calculate first term (~gamma) on r.h.s. of pressure evolution equation       
! 
       ALLOCATE(work1(ntheta,nzeta,ns),                                 & ! Allocate work space
                work2(ntheta,nzeta,ns),                                 &   
                work3(ntheta,nzeta,ns), stat=istat)   
       IF (istat .ne. 0) STOP 'Allocation1 failed in UPDATE_PRES'

!      Div(p*Flux) terms
!      Note that near axis, jvsups ~ s cos(u+nv), jvsupv -> s sin(nv) due to jacobian element
!      In particular, pijcf(1)*jvsupsijcf makes no contribution to forces (integrates to zero,
!      since jvsupsijcf(1) ~ cos(u) + O(s)!)
!      jvsupu may be finite, although dp/du -> 0 at axis
!      Note: WORK1 has COSINE parity, nsh = ns_i-1


!Upwind to suppress grid separation in radial direction
       IF (l_upwind) THEN
          work1 = -0.5_dp*signjac*hs_i*pijcf_ds*ABS(jvsupsijcf)
          work2 = pijch*jvsupsijcf + work1
       ELSE
          work2 = pijcf*jvsupsijcf
       END IF

       nblock = SIZE(work1,1)*SIZE(work1,2)
       CALL GradientHalf(work1, work2, nblock)
!!     work1(:,:,2:ns) = ohs*(work2(:,:,2:ns) - work2(:,:,1:nsh))          ! BC: JVSUPS(NS) = 0

       ALLOCATE (jvsupuijsh(ntheta,nzeta,ns),                           &   
                 jvsupvijsh(ntheta,nzeta,ns), stat=istat)
       IF (istat .NE. 0) STOP 'Allocation1 failed in UPDATE_PRES'

       CALL to_half_mesh(jvsupuijsf, jvsupuijsh, nblock)
       CALL to_half_mesh(jvsupvijsf, jvsupvijsh, nblock)
       work2 = pijch*jvsupuijsh
       work3 = pijch*jvsupvijsh
       
       ALLOCATE(work4(0:mpol,-ntor:ntor,ns),                            &
                work5(0:mpol,-ntor:ntor,ns),                            &
                work6(0:mpol,-ntor:ntor,ns), stat=istat)
       IF (istat .NE. 0) STOP 'Allocation2 failed in UPDATE_PRES'


!      FFT of Div(p*Flux) terms (multiply below by m,n*nfp)
       CALL tomnsp(work2, work5, psin, ihalf)
       CALL tomnsp(work3, work6, psin, ihalf)
       CALL tomnsp(work1, work4, pcos, ihalf)
              
       DO m = 0, mpol
          DO n = -ntor, ntor
             work4(m,n,:) = work4(m,n,:)                               &
                          + m*work5(m,n,:) + n*nfp*work6(m,n,:)         !m: Pol. deriv of SIN(m*theta+n*NFP*phi)
          ENDDO                                                         !n: Tor. deriv of SIN(m*theta+n*NFP*phi)
       ENDDO

       djpmnch = -gamma*work4                                           !Gamma = adiabatic factor
       
       DEALLOCATE(work5)
!   
!      Now, calculate next term ~(gamma-1) on r.h.s. of pressure equation
!      the u,v derivatives should be consistent with div(pv) term and combine
!      to cancel if div(v) = 0

       work1 = pijcf_ds*jvsupsijcf+pijcf_du*jvsupuijsf+pijcf_dv*jvsupvijsf
       CALL to_half_mesh(work1, work2, nblock)

!      MIXED HALF/FULL: DIAGONAL IN ANGULAR DERIVATIVES
!!       work1 = pijcf_ds*jvsupsijcf
!!       CALL to_half_mesh(work1,work2)
!!       work2 = work2+jvsupuijsh*pijch_du+jvsupvijsh*pijch_dv

       DEALLOCATE(jvsupuijsh, jvsupvijsh)
       
!Upwinding correction (makes 1st order streaming dominantly diagonal)
!NOTE: since it enters conservatively, it doesn't effect MHD force, only Hessian
!Can multiply this by a positive factor to increase numerical diffusion
       IF (l_upwind) THEN
          work1 = -0.5_dp*signjac*pijcf_ds*ABS(jvsupsijcf)
          work1(:,:,1) = 0;  work1(:,:,ns) = 0
          CALL GradientHalf(work3, work1, nblock)
          work2(:,:,2:ns) = work2(:,:,2:ns) + hs_i*work3(:,:,2:ns)
       END IF

       CALL tomnsp(work2, work4, pcos, ihalf)                           ! Calculate Fourier transform of WORK2


!      MAKES PRES-PERTURB A POSITIVE DEFINITE CONTRIBUTION TO HESSIAN: 
!      When computing hessian, guarantee positive-definiteness this way
!      (sound wave compressional energy contribution)
       gamm1 = gamma-1
!!       IF (l_Compute_Hessian) gamm1 = gamm1 + xPosDef
       IF (l_linearize) gamm1 = gamm1 + xPosDef
  
       djpmnch = djpmnch + gamm1*work4                                  ! Harmonics of r.h.s. of Pres. Eq. (full mesh) 

!SPH   Evolve harmonics of (jacobh * pres)mn to conserve pressure
       djpmnch = delta_t*djpmnch                                        ! calculate increment of pressure harmonics.
       djpmnch(:,:,1)  = 0
       djpmnch(m0,:,1) = djpmnch(m0,:,2)                                ! store js=1, m=0 value

!SPH   Suppress radial mesh separation together with Orszag 3/2 integration rule in metrics
!       djpmnch(m1:,:,2) = 0

       DEALLOCATE(work1, work2, work3, work4, work6)

       CALL second0(toff)
       time_update_pres = time_update_pres + (toff-ton)

       END SUBROUTINE update_pres

       END MODULE siesta_pressure
