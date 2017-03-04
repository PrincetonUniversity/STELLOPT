      MODULE quantities
!     
!     WRITTEN 08-18-06 BY R. SANCHEZ AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: This module contains the dynamical variables evolved: 
!        Namely, pressure, velocity and magnetic field.
!
      USE stel_kinds
      USE stel_constants      
      USE island_params, ns=>ns_i, ntheta=>nu_i, nzeta=>nv_i,           &
          mpol=>mpol_i, ntor=>ntor_i, nuv=>nuv_i, mnmax=>mnmax_i,       &
          ohs=>ohs_i, nfp=>nfp_i    !, mpol32=>mpol32_i, ntor32=>ntor32_i
      IMPLICIT NONE

!
!       Convention: The name of the variables may contain:
!           i.   The first one or two letters tells the quantity. So, "p"=pressure, "B" = magnetic field,
!                "jB" = magnetic field*jacobian, "v" = velocity, "jv" = velocity*jacobian 
!           ii.  The next three letters maybe "SUB" or "SUP", to distinguish covariant/contravariant 
!           ii.  The next letter maybe "s", "u" or "v" to distinguish radial/poloidal/toroidal components
!           iii. The next two letters maybe "MN" or "IJ", to distinguish Fourier harmonics or angular values 
!           iv.  The next letter maybe "C" or "S" to distinguish COSINE or SINE parity
!           v.   The last letter maybe "F" or "H" to distinguish FULL or HALF radial mesh
!
!       Radial meshes
!           velocity, force components       :: full mesh
!           p, contravariant components of B :: half mesh
!           jacobh (jacobian to curvilinear coordinates) :: half mesh (stores <jacobh(2)> in jacobh(1))
!     VARIABLE DECLARATIONS  
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       REAL(rprec) :: b_factor, w_factor, p_factor                     ! Scale internal energy to 1
       REAL(rprec) :: signjac
       REAL(rprec) :: wp, wb
       REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:)::                 &  ! Fourier harmonics= BSUP?, VSUB?, P
        fsupsmncf, fsupumnsf, fsupvmnsf,                            &
        jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch,                 & 
        djpmnch, djbsupsmnsh, djbsupumnch, djbsupvmnch,             &  ! Those starting with "d" represent increments.  
        jbsupsmnsh_p, jbsupumnch_p, jbsupvmnch_p, jpmnch_p             ! _p == PARALLEL (ns index reversed)        

!      Pointers are subsections of xc, gc (copied, do NOT use pointers!)
       REAL(rprec), POINTER, DIMENSION(:,:,:) ::                    &  ! Fourier harmonics= JVSUP?mn?F
!      Evolved quantities
        jvsupsmncf, jvsupumnsf, jvsupvmnsf,                         &
!      Computed forces
        fsubsmncf, fsubumnsf, fsubvmnsf                           

      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: jvsupsijcf,     &  ! all these are needed to compute convolutions
        jvsupuijsf, jvsupvijsf, jacobh, jacobf,                     &  ! in pressure and magn. field evolution equations.
        jacobh_p, jacobf_p

      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) ::                 &
        bsupsijsf0, bsupuijcf0, bsupvijcf0,                         &
        bsupsijsh0, bsupuijch0, bsupvijch0,                         &
        bsubsijsf,  bsubuijcf,  bsubvijcf, bsq

      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) ::                 &
        pijch, pijcf, pijcf_ds, pijcf_du, pijcf_dv

      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) ::                 &
        ksupsijsf0, ksupuijcf0, ksupvijcf0,                         &
        ksupsmnsf,  ksupumncf,  ksupvmncf

      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) ::                 &
        ksubsijsf, ksubuijcf, ksubvijcf                                 

      REAL(rprec) :: fbdy(13)                
!
!    Note that: DBSUP?MNF:              calculated in UPDATE_BFIELD;    used in UPDATE_FORCE and UPDATE_STATE
!               DPCMNF:                 calculated in UPDATE_PRES;      used in UPDATE_FORCE and UPDATE_STATE
!               BSUP?MN:                updated by UPDATE_STATE;        used in UPDATE_FORCE   
!               PMN:                    updated in UPDATE_STATE;        used in UPDATE_FORCE     
!               JVSUP?IJ:               updated in UPDATE_UPPERV;       used in UPDATE_BFIELD + UPDATE_PRES 
!               FSUB?MN:                updated in UPDATE_FORCE;        used in UPDATE_VEL
!               VSUB?MN:                updated by UPDATE_VEL;          used in UPDATE_UPPERV
!-----------------------------------------------

      CONTAINS

        SUBROUTINE init_quantities
        USE fourier, ONLY: tomnsp, toijsp
        IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER, PARAMETER :: m0=0
        INTEGER :: istat, js, m
        REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE:: presif, presih
!-----------------------------------------------
        w_factor = 1._dp/ABS(wb_i)
        p_factor = w_factor
        b_factor = SQRT(w_factor)
        gnorm_i  = gnorm_i/w_factor
!
!       Allocates and initializes dependent variables
!
        CALL alloc_quantities
!
!       Initialize full mesh covariant components of jac*Bsup-i (Only MN values are passed to other subroutines).
!       
        CALL init_bcovar
!
!       Initialize half mesh pressure (Only MN Values are passed to other subroutines)
!
        ALLOCATE(presih(ntheta,nzeta,ns),                               &
       &         presif(ntheta,nzeta,ns), stat = istat)
        IF (istat .ne. 0) STOP "Allocate error #1 in init_quantities"        

        DO js = 1, ns 
          presif(:,:,js) = p_factor*presf_i(js)
        ENDDO

        CALL to_half_mesh(presif, presih, nuv)
       
!SPH    Push jacob*(real pressure)
        presih = presih * jacobh

        CALL tomnsp(presih, jpmnch, 0, 1)                                 ! Initialize pressure harmonics (half mesh)
        jpmnch(:,:,1) = 0
        jpmnch(m0,:,1) = jpmnch(m0,:,2)

        DEALLOCATE(presif, presih, stat=istat)
        IF (istat .ne. 0) STOP "Deallocate error #2 in init_quantities"
       
        END SUBROUTINE init_quantities


        SUBROUTINE init_bcovar
        USE metrics, ONLY: lmns_i, sqrtg, iflipj
        USE fourier, ONLY: getorigin, tomnsp, toijsp
        IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER, PARAMETER :: m0=0, m1=1
        INTEGER  :: istat, js, l, m, n
        REAL(rprec) :: sum1, dif1
        REAL(rprec), DIMENSION(:), ALLOCATABLE :: phiph, chiph
!-----------------------------------------------
        jacobh(:,:,:) = RESHAPE(sqrtg, SHAPE(jacobh))                           ! Reshapes sqrtg in local form

!AVOID DIVIDE BY ZERO AND STORE AVERAGED VALUE (OVER THETA) 
!AT FIRST 1/2 GRID PT AT ORIGIN IN JACOBH(1)
        CALL getorigin (jacobh, m0, 0)
        signjac = jacobh(1,1,2)/ABS(jacobh(1,1,2))

        CALL to_full_mesh(jacobh, jacobf, nuv)
!Must be consistent with variational form of pressure evolution equation (DO NOT USE 2pt EXTRAPOLATIONS!)
        jacobf(:,:,1) = jacobh(:,:,1)
        jacobf(:,:,ns)= jacobh(:,:,ns)

        DO l=1,ns
           jacobf_p(:,:,l) = jacobf(:,:,l)
           jacobh_p(:,:,l) = jacobh(:,:,l)
        END DO

        ALLOCATE (vp_f(ns), stat=istat)
        IF (istat .ne. 0) STOP 'Allocation error in init_bcovar'        
        vp_f = 0
        DO js = 1, ns
           DO l = 1, ntheta
!              vp_h(js) = vp_h(js) + SUM(jacobh(l,:,js))*cosmui(l,m0)
              vp_f(js) = vp_f(js) + SUM(jacobf(l,:,js))*cosmui(l,m0)
           END DO                                                         
        END DO

!        vp_h = vp_h/(ntheta*nzeta)
        vp_f = vp_f/(ntheta*nzeta)

        DEALLOCATE (sqrtg)

!SPH 9/27/06 Convert to -jac*bsupu, bsupv components on HALF mesh
        
!Initially, lmns_i is on the FULL mesh, so convert to half mesh
!use jbsupvmnch to store half-grid lambdas
!!        CALL to_half_mesh_mn(lmns_i, jbsupvmnch)
!!        jbsupvmnch(1,m1,:) = 0
        jbsupsmnsh = 0


        ALLOCATE(phiph(ns), chiph(ns))

        phipf_i = signjac*iflipj*phipf_i
        chipf_i = signjac*iflipj*chipf_i
        phiph(2:ns) = (phipf_i(2:ns) + phipf_i(1:nsh))/2
        chiph(2:ns) = (chipf_i(2:ns) + chipf_i(1:nsh))/2
        phiph(1) = 0;  chiph(1) = 0

!REMEMBER: lmns_i is coefficient of sin(m*u+n*nfp*v)! (originally allocated to 1:ns+1)
!          and jbsupu ~ chip - phip*d(lambda)/dv); jbsupv ~ phip + phip*d(lambda)/du)

        DO m = 0,mpol
           DO n = -ntor, ntor
              jbsupumnch(m,n,:) = phiph(:)*(-n*nfp*lmns_i(m,n,1:ns))
              jbsupvmnch(m,n,:) = phiph(:)*( m    *lmns_i(m,n,1:ns))
           END DO
        END DO

        jbsupumnch(0,0,:) = chiph(:)
        jbsupvmnch(0,0,:) = phiph(:)

        jbsupumnch(m1,:,1) = 0                            ! Need to make jbsupsmnsh - r0jbsupumnch = 0 at js=2

        jbsupumnch = b_factor*jbsupumnch 
        jbsupvmnch = b_factor*jbsupvmnch 

!STORE m=0 component of jbsupvmnch at origin in half-grid js=1
        jbsupvmnch(m0,:,1) = jbsupvmnch(m0,:,2) 
        jbsupvmnch(m1:,:,1) = 0

        DEALLOCATE (phiph, chiph, stat=istat)
        IF (istat .ne. 0) STOP "Deallocate error in init_bcovar"
        
!        CALL update_current                             ! Initialize values of KSUP

        END SUBROUTINE init_bcovar


        SUBROUTINE alloc_quantities
!
!       Allocates variables that will be needed through the advancing subroutines. 
!       More maybe added in the future. 
!       All variables allocated here are then deallocated in DEALLOC_QUANTITIES, 
!       which should be called at the end of the run.
! 
        IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER :: istat
!-----------------------------------------------
!       Half mesh quantities (harmonics) evolved: 3 magnetic field components and pressure                
        ALLOCATE(jbsupsmnsh(0:mpol,-ntor:ntor,ns),                      &  ! Half mesh quantities (harmonics)
                 jbsupumnch(0:mpol,-ntor:ntor,ns),                      &  ! B^s (sine), B^u (cosine), B^v (cosine)
                 jbsupvmnch(0:mpol,-ntor:ntor,ns),                      & 
                 jpmnch    (0:mpol,-ntor:ntor,ns), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation failed in ALLOC_QUANTITIES'     ! CONTRAVARIANT MAGNETIC FIELD
        jbsupsmnsh = zero; jbsupumnch = zero; jbsupvmnch = zero        
        jpmnch = zero

        ALLOCATE(bsq(ntheta,nzeta,ns), stat=istat)  
        IF (istat .NE. 0) STOP 'Allocation failed in ALLOC_QUANTITES'

        ALLOCATE(fsupsmncf(0:mpol,-ntor:ntor,ns),                       &  ! Full mesh quantities (harmonics)
                 fsupumnsf(0:mpol,-ntor:ntor,ns),                       &  ! F_s (cosine), F_u (sine), F_v (sine)
                 fsupvmnsf(0:mpol,-ntor:ntor,ns), stat=istat) 
        IF (istat .NE. 0) STOP 'Allocation failed in ALLOC_QUANTITIES'     ! COVARIANT MHD FORCE
                
        ALLOCATE(jacobh(ntheta,nzeta,ns),                               &
                 jacobf(ntheta,nzeta,ns),                               &
                 jacobh_p(ntheta,nzeta,ns),                             &
                 jacobf_p(ntheta,nzeta,ns), stat=istat)
        IF (istat .NE. 0) STOP 'Allocation failed in ALLOC_QUANTITIES'

        END SUBROUTINE alloc_quantities
        
        
        SUBROUTINE dealloc_quantities
        USE metrics, ONLY: lmns_i
!
!       Deallocates all variables allocated in ALLOC_QUANTITIES
!        
        IMPLICIT NONE
        INTEGER :: istat
        
        DEALLOCATE(jpmnch, jbsupsmnsh, jbsupumnch, jbsupvmnch, djpmnch, &
          djbsupsmnsh, djbsupumnch, djbsupvmnch, jacobh, jacobf, vp_f,  &
          jacobh_p, jacobf_p, jvsupsijcf, jvsupuijsf, jvsupvijsf,       &
          ksubsijsf, ksubuijcf, ksubvijcf, pijch, pijcf,                &
          pijcf_ds, fsupsmncf, fsupumnsf, fsupvmnsf, bsupuijcf0,        &
          bsupvijcf0, bsupsijsf0, bsq, pijcf_du, pijcf_dv,              &
          ksupsmnsf,  ksupumncf, ksupvmncf, stat=istat)

        IF(ALLOCATED(jpmnch_p))   DEALLOCATE(jpmnch_p, jbsupsmnsh_p,    &
                                             jbsupumnch_p, jbsupvmnch_p)
        IF(ALLOCATED(bsupsijsf0)) DEALLOCATE(bsupsijsf0, bsupuijcf0, bsupvijcf0)
        IF(ALLOCATED(bsupsijsh0)) DEALLOCATE(bsupsijsh0, bsupuijch0, bsupvijch0)
        IF(ALLOCATED(ksupsijsf0)) DEALLOCATE(ksupsijsf0, ksupuijcf0, ksupvijcf0)
        IF(ALLOCATED(bsubsijsf))  DEALLOCATE(bsubsijsf,  bsubuijcf, bsubvijcf)
        IF(ALLOCATED(lmns_i)) DEALLOCATE(lmns_i)
 
        IF (istat .NE. 0) STOP 'Deallocation failed in Quantities!'
        

        END SUBROUTINE dealloc_quantities
        
               
      END MODULE quantities
          
