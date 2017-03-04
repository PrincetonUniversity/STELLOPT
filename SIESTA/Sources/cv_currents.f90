#if defined(SKS)
      SUBROUTINE cv_currents_par (bsupsijsh, bsupuijch, bsupvijch,      &
                                  ksupsijsf, ksupuijcf, ksupvijcf,      &
                                  lmagen, lcurr)
!     
!     WRITTEN 02-26-07 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes contravariant currents on full-grid from 
!              covariant B-fields on half-grid
!
      USE stel_kinds
      USE stel_constants
      USE fourier 
      USE quantities
      USE island_params, ONLY: nu_i
      USE evolution, ONLY: l_natural, jsju_ratio, ste 
      USE shared_data, ONLY: l_update_state
      USE timer_mod
      USE Hessian, ONLY: mupar_norm
      USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
      IMPLICIT NONE
      INCLUDE 'mpif.h'
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ntheta,nzeta,ns), INTENT(IN)  ::           &  !FOR NOW, ASSUME THESE ARE ARRAYS FOR ALL NS!
           bsupsijsh, bsupuijch, bsupvijch
      REAL(rprec), DIMENSION(ntheta,nzeta,ns), INTENT(OUT) ::           &
           ksupsijsf, ksupuijcf, ksupvijcf
      LOGICAL, INTENT(in)  :: lmagen, lcurr
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: ifull=0, ihalf=1, m0=0, m1=1, m2=2
      INTEGER :: m, n, istat, l, js, nblock
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) ::                     &
         bsubsijsh, bsubuijch, bsubvijch,                               & 
         bsubsmnsh, bsubumnch, bsubvmnch,                               &
         bsubsmnsf, bsubumncf, bsubvmncf,                               &
         work1, work2, bsq_loc
      REAL(rprec) :: beta, r0, temp(-ntor:ntor), skston, skstoff,       &
                     ton, toff, stes(6), rtes(6), wb0
      INTEGER     :: k, nsmin, nsmax, ns_span, n1, n2
!-----------------------------------------------
      CALL second0(ton)
      skston=ton

      nsmin=MAX(1, startglobrow); nsmax=MIN(ns,endglobrow+1)
      n1=MAX(1, startglobrow); n2=MIN(ns,endglobrow)

      IF (lcurr) THEN
         nsmin=MAX(1, startglobrow-1); nsmax=MIN(ns,endglobrow+2)
      END IF

      ns_span = nsmax-nsmin+1
!
! Calculate covariant BsubXh on half mesh 
!       
      ALLOCATE(bsubsijsh(ntheta,nzeta,nsmin:nsmax),                     &
               bsubuijch(ntheta,nzeta,nsmin:nsmax),                     &
               bsubvijch(ntheta,nzeta,nsmin:nsmax), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation1 failed in CV_CURRENTS'
     
      CALL tolowerh_par(bsupsijsh, bsupuijch, bsupvijch,                &
                        bsubsijsh, bsubuijch, bsubvijch,                &
                        nsmin, nsmax)  

      IF (lmagen) THEN
         ALLOCATE (bsq_loc(ntheta,nzeta,nsmin:nsmax), stat=istat)
         bsq_loc = bsupsijsh(:,:,1:ns_span)*bsubsijsh                   & 
                 + bsupuijch(:,:,1:ns_span)*bsubuijch                   &
                 + bsupvijch(:,:,1:ns_span)*bsubvijch
         IF (n1 .EQ. 1) bsq_loc(:,:,1) = 0
         IF (n1 .GT. 1) THEN
            IF(ANY(bsq_loc(:,:,n1:n2) .LE. zero)) STOP 'BSQ <= 0 in cv_currents_par!'
         END IF

         wb = 0
         DO l = 1, nu_i
            wb = wb + SUM(bsq_loc(l,:,n1:n2)*jacobh_p(l,:,n1:n2))*cosmui(l,m0)
         END DO         
         wb0 = 2*signjac*(twopi*pi*hs_i)*wb/2                              !2X : theta 0-pi
         CALL MPI_ALLREDUCE(wb0,wb,1,MPI_REAL8, MPI_SUM,                &
                            MPI_COMM_WORLD,MPI_ERR)

!COMPUTE PARALLEL DAMPING SCALING COEFFICIENT <Bsupv**2/B**2> ~ k^2
         IF (.NOT. ALLOCATED(mupar_norm)) THEN
            ALLOCATE(mupar_norm(nsmin:nsmax), stat=istat)
            IF (istat .NE. 0) STOP 'MUPAR_NORM ALLOCATION FAILED!'
            mupar_norm = 0; k=0
            DO js = nsmin, nsmax
               k = k+1
               IF (js .EQ. 1) CYCLE
               DO l = 1, nu_i
                  mupar_norm(js) = mupar_norm(js) +                     &
                  SUM(bsupvijch(l,:,k)**2/bsq_loc(l,:,js))*cosmui(l,m0)
               END DO
            END DO
            beta = MAX(1.E-5_dp, wp_i/wb_i)
            IF (nsmin .EQ. 1) mupar_norm(1) = mupar_norm(2)
            mupar_norm = -signjac*beta*mupar_norm
         END IF

         DEALLOCATE (bsq_loc)
      END IF

!
! Calculate BsubXmnh Fourier coefficients
!       
      ALLOCATE(bsubsmnsh(0:mpol,-ntor:ntor,nsmin:nsmax),                & ! First, get BSUBMN on half mesh
               bsubumnch(0:mpol,-ntor:ntor,nsmin:nsmax),                &
               bsubvmnch(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)   
      IF (istat .ne. 0) STOP 'Allocation2 failed in CV_CURRENTS'

      CALL tomnsp_par(bsubsijsh, bsubsmnsh, 1, ifull, 1, ns_span)    ! Calculate BSUB?MN?H
      CALL tomnsp_par(bsubuijch, bsubumnch, 0, ifull, 1, ns_span)
      CALL tomnsp_par(bsubvijch, bsubvmnch, 0, ifull, 1, ns_span)
      CALL second0(toff)

!
!     Calculate full mesh, contravariant current components
!     KsupXF = sqrt(g)*JsupXF 
!        
                    
      ALLOCATE(bsubsmnsf(0:mpol,-ntor:ntor,nsmin:nsmax),                &  ! Next, get BSUBMN on full mesh
               bsubumncf(0:mpol,-ntor:ntor,nsmin:nsmax),                &
               bsubvmncf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
      IF (istat .ne. 0) STOP 'Allocation3 failed in CV_CURRENTS'
 
      nblock=SIZE(bsubsmnsh,1)*SIZE(bsubsmnsh,2)
      CALL to_full_mesh_par(bsubsmnsh, bsubsmnsf, nblock, nsmin, nsmax)                       
      CALL to_full_mesh_par(bsubumnch, bsubumncf, nblock, nsmin, nsmax)
      CALL to_full_mesh_par(bsubvmnch, bsubvmncf, nblock, nsmin, nsmax)
      CALL second0(toff)

!     NOTE: bsubumncf is zero, even for m=1, since guu ~ s^2, guv ~ s
      IF (nsmin .EQ. 1) THEN
         bsubsmnsf(:,:,1) = 0
         bsubumncf(:,:,1) = 0
         bsubsmnsf(m1,:,1) = bsubsmnsh(m1,:,1)
         bsubumncf(m1,:,1) = bsubumnch(m1,:,1)

!     VARIATIONAL FORM BELOW: DOES NOT WORK WELL - PHYSICALLY, ERRONEOUS COUPLING
!     TO M > 1 MODES IN BSUPX WITH METRIC ELEMENTS OCCURS

!     NOTE: the contributions from gvu B^u + gvs B^s cancel at the origin,
!     so there is ONLY an m=0 contribution to B_v but need m=1 too for correct
!     limit at js=1 for k^s,u(m1)
         bsubvmncf(:,:,1) = 0
         bsubvmncf(m0:m1,:,1) = bsubvmnch(m0:m1,:,2)
      END IF
!
!     Compute K^X == sqrt(g)J^X contravariant current components in Fourier space
!             
!                   
      ALLOCATE(work1(0:mpol,-ntor:ntor,nsmin:nsmax),                    &
               work2(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)  
      IF (istat .ne. 0) STOP 'Allocation4 failed in CV_CURRENTS'
      CALL second0(toff)

!     Radial gradients of bsubmnu,v on full grid (excluding origin: treat separately)
      CALL GradientFull_par(work1, bsubumnch, nblock, nsmin, nsmax)
      CALL GradientFull_par(work2, bsubvmnch, nblock, nsmin, nsmax)

      CALL second0(toff)

!     Approximate end values...(factor of 2 at origin due to 1/2*hs interval)
!     Forces will be multiplied by 1/2 in get_force_harmonics
      IF (nsmin .EQ. 1) THEN
         r0 = hs_i/2
         work1(m0:m1,:,1) = 2*ohs*bsubumnch(m0:m1,:,2)
         bsubvmncf(m1,:,1) = bsubvmnch(m1,:,2)
         work2(m1,:,1)    = 2*ohs*bsubvmnch(m1,:,2)         !--> 0
         bsubumncf(m1,:,1) = r0*bsubsmnsf(m1,:,1)           !--> 0
      END IF

!     NOTE: ksupsmnsf(ns) is REALLY ksupsmnsh(ns) a half-grid point in from edge
!           to satisfy k^s = r0*k^u at s=r0, for m=1, must have b_u = r0*b_s for m=1
      CALL second0(toff)

      IF (lcurr) THEN
         nsmax=MIN(ns,endglobrow+1)
      ELSE
         nsmax=MIN(ns,endglobrow)
      END IF
      ns_span = nsmax-nsmin+1
      DO n = -ntor, ntor
         DO m = 0, mpol
           ksupsmnsf(m,n,nsmin:nsmax) = -m*bsubvmncf(m,n,nsmin:nsmax) + n*nfp*bsubumncf(m,n,nsmin:nsmax)
           ksupumncf(m,n,nsmin:nsmax) =  n*nfp*bsubsmnsf(m,n,nsmin:nsmax) - work2(m,n,nsmin:nsmax)                     
           ksupvmncf(m,n,nsmin:nsmax) = -m*bsubsmnsf(m,n,nsmin:nsmax)     + work1(m,n,nsmin:nsmax)
         ENDDO
      ENDDO
      CALL second0(toff)

      IF (nsmin .EQ. 1) THEN
!Diagnostic output: K^s-r0*K^u at origin BEFORE bc applied
         temp = (ksupsmnsf(m1,:,1) + r0*ksupumncf(m1,:,1))/2
         IF (l_update_state) THEN
         jsju_ratio = SUM((2*temp)**2)
            IF (jsju_ratio .GT. zero) THEN
               jsju_ratio = SUM((ksupsmnsf(m1,:,2) - r0*ksupumncf(m1,:,2))**2) / jsju_ratio
               jsju_ratio = SQRT(jsju_ratio)
            END IF
         END IF

!Apply exact bc K^s = r0*K^u for m=1 at origin
         IF (l_natural) ksupsmnsf(m1,:,1) = r0*ksupumncf(m1,:,1)
         
         ksupsmnsf(m0,:,1) = 0;  ksupsmnsf(m2:,:,1) = 0
         ksupumncf(m0,:,1) = 0;  ksupumncf(m2:,:,1) = 0
         ksupvmncf(m1:,:,1)= 0
      END IF

      DEALLOCATE (work1, work2, bsubsmnsf, bsubumncf, bsubvmncf)

!
!     Compute spectral truncation error measure
!     (Use ksupsijcf for temporary FILTERED bsubXijsh values)
!
      IF (l_update_state) THEN
         ns_span = n2-n1+1
         CALL toijsp_par(bsubsmnsh(:,:,n1:n2), ksupsijsf, 0, 0, 1, 0, 1, ns_span) 
         CALL toijsp_par(bsubumnch(:,:,n1:n2), ksupuijcf, 0, 0, 0, 0, 1, ns_span)
         CALL toijsp_par(bsubvmnch(:,:,n1:n2), ksupvijcf, 0, 0, 0, 0, 1, ns_span)
         stes(1) = SUM((ksupsijsf(:,:,1:ns_span)-bsubsijsh(:,:,n1:n2))**2)
         stes(2) = SUM((ksupsijsf(:,:,1:ns_span)+bsubsijsh(:,:,n1:n2))**2)         
         stes(3) = SUM((ksupuijcf(:,:,1:ns_span)-bsubuijch(:,:,n1:n2))**2)
         stes(4) = SUM((ksupuijcf(:,:,1:ns_span)+bsubuijch(:,:,n1:n2))**2)         
         stes(5) = SUM((ksupvijcf(:,:,1:ns_span)-bsubvijch(:,:,n1:n2))**2)
         stes(6) = SUM((ksupvijcf(:,:,1:ns_span)+bsubvijch(:,:,n1:n2))**2)         
         CALL MPI_ALLREDUCE(stes,rtes,6,MPI_REAL8,MPI_SUM,              &
                            MPI_COMM_WORLD,MPI_ERR)
         IF (rtes(2) .NE. zero) ste(2) = SQRT(rtes(1)/rtes(2))
         IF (rtes(4) .NE. zero) ste(3) = SQRT(rtes(3)/rtes(4))
         IF (rtes(6) .NE. zero) ste(4) = SQRT(rtes(5)/rtes(6))
         ns_span = nsmax-nsmin+1
      END IF

      DEALLOCATE (bsubsmnsh, bsubumnch, bsubvmnch,                      &
                  bsubsijsh, bsubuijch, bsubvijch) 

      CALL second0(toff)

      CALL toijsp_par(ksupsmnsf(:,:,nsmin:), ksupsijsf, 0, 0, 1, 0, 1, ns_span) 
      CALL toijsp_par(ksupumncf(:,:,nsmin:), ksupuijcf, 0, 0, 0, 0, 1, ns_span)
      CALL toijsp_par(ksupvmncf(:,:,nsmin:), ksupvijcf, 0, 0, 0, 0, 1, ns_span)

      IF (lcurr) THEN
         CALL tolowerf_par(ksupsijsf, ksupuijcf, ksupvijcf,             & 
                           ksubsijsf, ksubuijcf, ksubvijcf, nsmin, nsmax)
!DON'T KNOW WHY THIS WOULD BE NEEDED?
!         IF (nsmax .EQ. ns) THEN
!            ksubsijsf(:,:,ns) = 0
!            ksubuijcf(:,:,ns) = 0
!            ksubvijcf(:,:,ns) = 0
!         END IF
         IF (nsmin .EQ. 1) ksubuijcf(:,:,1) = 0 
      ENDIF


      CALL second0(toff)
      skstoff = toff
      cv_current_time=cv_current_time+(skstoff-skston)

      time_current = time_current+(skstoff-skston)

      END SUBROUTINE cv_currents_par
#endif

      SUBROUTINE cv_currents (bsupsijsh, bsupuijch, bsupvijch,          &
                              ksupsijsf, ksupuijcf, ksupvijcf,          &
                              lmagen, lcurr)
!     
!     WRITTEN 02-26-07 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes contravariant currents on full-grid from 
!              contravariant B-fields on half-grid
!
      USE stel_kinds
      USE stel_constants
      USE fourier 
      USE quantities
      USE island_params, ONLY: nu_i
      USE evolution, ONLY: l_natural, jsju_ratio, ste
      USE shared_data, ONLY: l_update_state
      USE Hessian, ONLY: mupar_norm
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ntheta,nzeta,ns), INTENT(in)  ::           &
           bsupsijsh, bsupuijch, bsupvijch
      REAL(rprec), DIMENSION(ntheta,nzeta,ns), INTENT(out) ::           &
           ksupsijsf, ksupuijcf, ksupvijcf
      LOGICAL, INTENT(in)  :: lmagen, lcurr
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: ifull=0, ihalf=1, m0=0, m1=1, m2=2
      INTEGER :: m, n, istat, l, js, i, j, k, nblock
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) ::                     &
         bsubsijsh, bsubuijch, bsubvijch,                               & 
         bsubsmnsh, bsubumnch, bsubvmnch,                               &
         bsubsmnsf, bsubumncf, bsubvmncf,                               &
         work1, work2
      REAL(rprec) :: beta, r0, temp(-ntor:ntor)
      REAL(rprec) :: skston, skstoff
      INTEGER :: nsmin, nsmax
!-----------------------------------------------
!
! Calculate BsubXh on half mesh 
!       
      CALL second0(skston)

      ALLOCATE(bsubsijsh(ntheta,nzeta,ns),                              &
               bsubuijch(ntheta,nzeta,ns),                              &
               bsubvijch(ntheta,nzeta,ns), stat=istat)
      IF (istat .NE. 0) STOP 'Allocation1 failed in CV_CURRENTS'
     
      CALL tolowerh(bsupsijsh, bsupuijch, bsupvijch,                    &
                    bsubsijsh, bsubuijch, bsubvijch)  

      IF (lmagen) THEN
         bsq = bsupsijsh*bsubsijsh + bsupuijch*bsubuijch                &
             + bsupvijch*bsubvijch
         IF (ANY(bsq(:,:,2) .le. zero)) STOP 'BSQ <= 0 in cv_currents!'

         wb = zero
         DO l = 1, nu_i
            wb = wb + SUM(bsq(l,:,2:ns)*jacobh(l,:,2:ns))*cosmui(l,m0)
         END DO         
         wb = 2*signjac*(twopi*pi*hs_i)*wb/2                              !2X : theta 0-pi

!COMPUTE PARALLEL DAMPING SCALING COEFFICIENT <Bsupv**2/B**2> ~ k^2
         IF (.NOT. ALLOCATED(mupar_norm)) THEN
            ALLOCATE(mupar_norm(ns), stat=istat)
            IF (istat .NE. 0) STOP 'MUPAR_NORM ALLOCATION FAILED!'
            mupar_norm = 0
            DO js = 2, ns
               DO l = 1, nu_i
                  mupar_norm(js) = mupar_norm(js) +                     &
                  SUM(bsupvijch(l,:,js)**2/bsq(l,:,js))*cosmui(l,m0)
               END DO
            END DO
            beta = MAX(1.E-5_dp, wp_i/wb_i)
            mupar_norm(1) = mupar_norm(2)
            mupar_norm = -signjac*beta*mupar_norm
         END IF
      END IF

!
! Calculate BsubXmnh Fourier coefficients
!       
      ALLOCATE(bsubsmnsh(0:mpol,-ntor:ntor,ns),                         & ! First, get BSUBMN on half mesh
               bsubumnch(0:mpol,-ntor:ntor,ns),                         &
               bsubvmnch(0:mpol,-ntor:ntor,ns), stat=istat)   
      IF (istat .ne. 0) STOP 'Allocation2 failed in CV_CURRENTS'

      CALL tomnsp(bsubsijsh, bsubsmnsh, 1, ifull)                         ! Calculate BSUB?MN?H
      CALL tomnsp(bsubuijch, bsubumnch, 0, ifull)
      CALL tomnsp(bsubvijch, bsubvmnch, 0, ifull)
!
!     Calculate full mesh, contravariant current components
!     KsupXF = sqrt(g)*JsupXF 
!        
      ALLOCATE(bsubsmnsf(0:mpol,-ntor:ntor,ns),                         &  ! Next, get BSUBMN on full mesh
               bsubumncf(0:mpol,-ntor:ntor,ns),                         &
               bsubvmncf(0:mpol,-ntor:ntor,ns), stat=istat)  
      IF (istat .ne. 0) STOP 'Allocation3 failed in CV_CURRENTS'
 
      nblock=SIZE(bsubsmnsh,1)*SIZE(bsubsmnsh,2)
      CALL to_full_mesh(bsubsmnsh, bsubsmnsf, nblock)                       
      CALL to_full_mesh(bsubumnch, bsubumncf, nblock)
      CALL to_full_mesh(bsubvmnch, bsubvmncf, nblock)

!     NOTE: bsubumncf is zero, even for m=1, since guu ~ s^2, guv ~ s
      bsubsmnsf(:,:,1) = 0
      bsubumncf(:,:,1) = 0
      bsubsmnsf(m1,:,1) = bsubsmnsh(m1,:,1)
      bsubumncf(m1,:,1) = bsubumnch(m1,:,1)

!     VARIATIONAL FORM BELOW: DOES NOT WORK WELL - PHYSICALLY, ERRONEOUS COUPLING
!     TO M > 1 MODES IN BSUPX WITH METRIC ELEMENTS OCCURS
!      bsubsmnsf(m1,:,1) = bsubsmnsh(m1,:,2)
!      bsubumncf(m1,:,1) = bsubumnch(m1,:,2)

      r0 = hs_i/2
!      bsubumncf(m1,-ntor:-1,1) = r0*bsubsmnsf(m1,-ntor:-1,1)
!      bsubumncf(m1,1:ntor,1)   = r0*bsubsmnsf(m1,1:ntor,1)

!     NOTE: the contributions from gvu B^u + gvs B^s cancel at the origin,
!     so there is ONLY an m=0 contribution to B_v but need m=1 too for correct
!     limit at js=1 for k^s,u(m1)
      bsubvmncf(:,:,1) = 0
      bsubvmncf(m0:m1,:,1) = bsubvmnch(m0:m1,:,2)

      bsubsmnsf(:,:,ns) = bsubsmnsh(:,:,ns)
      bsubumncf(:,:,ns) = bsubumnch(:,:,ns)
      bsubvmncf(:,:,ns) = bsubvmnch(:,:,ns)

!
!     Compute K^X == sqrt(g)J^X contravariant current components in Fourier space
!             
!                   
      ALLOCATE(work1(0:mpol,-ntor:ntor,ns),                             &
               work2(0:mpol,-ntor:ntor,ns), stat=istat)  
      IF (istat .ne. 0) STOP 'Allocation4 failed in CV_CURRENTS'

!     Radial gradients of bsubmnu,v on full grid (excluding origin: treat separately)
      CALL GradientFull(work1, bsubumnch, nblock)
      CALL GradientFull(work2, bsubvmnch, nblock)

!     Approximate end values...(factor of 2 at origin due to 1/2*hs interval)
!     Forces will be multiplied by 1/2 in get_force_harmonics
      work1(m0:m1,:,1) = 2*ohs*bsubumnch(m0:m1,:,2)
      work2(m1,:,1) = 2*ohs*bsubvmnch(m1,:,2)
      bsubvmncf(m1,:,1) = bsubvmnch(m1,:,2)
      bsubumncf(m1,:,1) = r0*bsubsmnsf(m1,:,1)

!     NOTE: ksupsmnsf(ns) is REALLY ksupsmnsh(ns) a half-grid point in from edge
!           to satisfy k^s = r0*k^u at s=r0, for m=1, must have b_u = r0*b_s for m=1
      DO m = 0, mpol
         DO n = -ntor, ntor
           ksupsmnsf(m,n,:) = -m*bsubvmncf(m,n,:)                       &
                            +  n*nfp*bsubumncf(m,n,:)
           ksupumncf(m,n,:) =  n*nfp*bsubsmnsf(m,n,:) - work2(m,n,:)                     
           ksupvmncf(m,n,:) = -m*bsubsmnsf(m,n,:)     + work1(m,n,:)
         ENDDO
      ENDDO

!Diagnostic output: K^s-r0*K^u at origin BEFORE bc applied
      temp = (ksupsmnsf(m1,:,1) + r0*ksupumncf(m1,:,1))/2
      IF (l_update_state) THEN
        jsju_ratio = SUM((2*temp)**2)
        IF (jsju_ratio .GT. zero) THEN
          jsju_ratio = SUM((ksupsmnsf(m1,:,1) - r0*ksupumncf(m1,:,1))**2) / jsju_ratio
          jsju_ratio = SQRT(jsju_ratio)
        END IF
      END IF

!Apply exact bc K^s = r0*K^u for m=1 at origin
      IF (l_natural) ksupsmnsf(m1,:,1) = r0*ksupumncf(m1,:,1)

      ksupsmnsf(m0,:,1) = 0;  ksupsmnsf(m2:,:,1) = 0
      ksupumncf(m0,:,1) = 0;  ksupumncf(m2:,:,1) = 0
      ksupvmncf(m1:,:,1)= 0

      DEALLOCATE (work1, work2)
      DEALLOCATE(bsubsmnsf, bsubumncf, bsubvmncf)
!
!     Compute spectral truncation error measure
!     (Use ksupsijcf for temporary FILTERED bsubXijsh values)
!
      IF (l_update_state) THEN
         CALL toijsp(bsubsmnsh, ksupsijsf, 0, 0, 1, 0) 
         CALL toijsp(bsubumnch, ksupuijcf, 0, 0, 0, 0)
         CALL toijsp(bsubvmnch, ksupvijcf, 0, 0, 0, 0)
         r0 = SUM((ksupsijsf+bsubsijsh)**2)
         IF (r0 .NE. zero) ste(2) = SQRT(SUM((ksupsijsf-bsubsijsh)**2)/r0)
         r0 = SUM((ksupuijcf+bsubuijch)**2)
         IF (r0 .NE. zero) ste(3) = SQRT(SUM((ksupuijcf-bsubuijch)**2)/r0)
         r0 = SUM((ksupvijcf+bsubvijch)**2)         
         IF (r0 .NE. zero) ste(4) = SQRT(SUM((ksupvijcf-bsubvijch)**2)/r0)
      END IF

      DEALLOCATE(bsubsmnsh, bsubumnch, bsubvmnch)
      DEALLOCATE(bsubsijsh, bsubuijch, bsubvijch) 

      CALL toijsp(ksupsmnsf, ksupsijsf, 0, 0, 1, 0) 
      CALL toijsp(ksupumncf, ksupuijcf, 0, 0, 0, 0)
      CALL toijsp(ksupvmncf, ksupvijcf, 0, 0, 0, 0)

!     NEED THESE FOR RESISTIVE CONTRIBUTIONS (DIVIDE OUT jac FACTOR LATER)
!     ksubX = jac*JsubX;  ksupX = jac*JsupX
!
      IF (lcurr) THEN
         CALL tolowerf(ksupsijsf, ksupuijcf, ksupvijcf,                  & 
                       ksubsijsf, ksubuijcf, ksubvijcf)

!        ksubsijsf(:,:,ns) = 0; ksubuijcf(:,:,ns) = 0
         ksubuijcf(:,:,1)  = 0  
!        ksubvijcf(:,:,ns) = 0
      ENDIF

      CALL second0(skstoff)
      time_current = time_current+(skstoff-skston)

      END SUBROUTINE cv_currents
