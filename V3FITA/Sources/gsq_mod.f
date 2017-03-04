!     SPH: INTEGER(iprec) -> INTEGER
!*******************************************************************************
!  File gsq_mod.f
!  Contains module gsq_mod
!  Defines derived-type: (none at this time)
!    This deals with the function g^2, which V3FITA minimizes to do the
!    reconstruction. It is closely related to a chi-squared.
!    Now (9/2006) the specification of the function is pretty simple - 
!       all the signals are used, and the weights are all 1.
!    NB. At this time (10-2006) all signals are assumed to only contribute
!        one scalar value to the g^2. This may change in the future.
!
!       MODULE VARIABLES
!  gsq_nsig               integer, number of signals
!  gsq_nrp                integer, number of (reconstruction) parameters
!  gsq_nsv                integer, number of singular values. MIN(nsig,nrp)
!  gsq_g2                 real, value of g^2 function
!  gsq_sweight            real array, weights for signals
!  gsq_ssigma             real array, signal sigmas, combined observed and model
!  gsq_evector            real array, normalized error vector (_observe - _model)
!  gsq_ppi                real array, normalizing factors for parameters
!  gsq_jacobian_observe   real 2d array, jacobian of the observed signals wrt
!                         parameters
!  gsq_jacobian_model     real 2d array, jacobian of the model signals wrt
!                         parameters
!  gsq_jacobian           real 2d array, combined jacobian (_model - _observe)
!  gsq_jacobian_norm      real 2d array, combined jacobian, normalized. (A matrix)
!  gsq_jsvd_u             real 2d array, U portion of SVD of A (jacobian_norm)
!  gsq_jsvd_vt            real 2d array, V Transpose portion of SVD of A
!  gsq_jsvd_w             real vector, singular values of SVD of A
!  gsq_ata                real 2d array, A-Transpose . A, where A is the
!                            gsq_jacobian_norm matrix
!  gsq_beta               real vector, - A-Transpose . gsq_evector
!  gsq_delta_a_sd         real vector, steepest descent step
!  gsq_delta_a_svd        real 2d array, list of k-SVD steps
!  gsq_g2exp              real 1d array, list of expected g2 values
!  gsq_dg2exp             real 1d array, list of expected changes in g2 values
!  gsq_dg2exp_lin         real 1d array, list of linear portion of dg2exp
!  gsq_dg2exp_quad        real 1d array, list of quadratic portion of dg2exp
!  gsq_delta_a_len        real 1d array, lengths of delta_a_svd vectors
!  gsq_exp_eff            real 1d array, expected step efficiencies
!  gsq_marg_exp_eff       real 1d array, expected marginal step efficiencies

!       MODULE SUBROUTINES
!    gsq_initialize
!    gsq_evaluate_g2(sdm_a_arg,sdo_a_arg)
!    gsq_evaluate_jac(model_arg,sdm_a_arg,sdo_a_arg,rparam_a_arg)
!    gsq_evaluate_dar(svdc_arg,delta_a_arg,delta_r_arg)
!    gsq_write_jac(identifier,unit,verbose)
!    gsq_write_conf(rpa,identifier,unit,verbose)
!    gsq_write_mr(rpa,identifier,unit,verbose)

!*******************************************************************************
!  MODULE gsq_mod
!    (g^2 function module, for the V3FIT code)
! SECTION I.    VARIABLE DECLARATIONS
! SECTION II.   DERIVED-TYPE DECLARATIONS
! SECTION III.  INTERFACE BLOCKS
! SECTION IV.   INITIALIZATION SUBROUTINES
! SECTION V .   EVALUATION SUBROUTINES
! SECTION VII.  OUTPUT SUBROUTINES

! SECTION XII.  AUXILIARY SUBROUTINES
! SECTION XIII. DEBUGGING SUBROUTINES
! SECTION XV.   DUPLICATE CODING FOR TESTING
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE gsq_mod

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!  Frequently used mathematical constants, lots of extra precision.
!-------------------------------------------------------------------------------
      USE stel_kinds, only : rprec, cprec
      USE stel_constants, only : pi, twopi, one, zero
     
!-------------------------------------------------------------------------------
!  Use Statements for other structures, V3 Utilities
!    Module also USEs:
!       v3f_global               in    _initialize
!       signal_mc                in    _evaluate_jac
!       eq_interface             in    _evaluate_jac
!       v3f_global               in    _evaluate_jac
!       recon_param_model        in    _evaluate_jac
!-------------------------------------------------------------------------------
      USE v3_utilities
      USE signal_T
      USE recon_param_T
      USE model_T
      USE eq_interface

!-------------------------------------------------------------------------------
!  Implicit None comes after USE statements, before other declarations
!-------------------------------------------------------------------------------
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Make type declarations and constants Private, so there are no conflicts.
!-------------------------------------------------------------------------------
      PRIVATE rprec, cprec, pi, twopi, one, zero

!-------------------------------------------------------------------------------
!  Variables for g^2 specification and computation
!-------------------------------------------------------------------------------
!
      INTEGER :: gsq_nsig, gsq_nrp, gsq_nsv
      REAL(rprec) :: gsq_g2

!  _ns size arrays
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: gsq_sweight,                   & 
     &   gsq_ssigma, gsq_evector

!  _np size array
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: gsq_ppi

!  _ns by _np size arrays (jacobians)
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: gsq_jacobian_observe,        &
     &   gsq_jacobian_model, gsq_jacobian, gsq_jacobian_norm

!  Jacobian SVD arrays
!     u is ns x ns
!     vt is np by np
!     w is a vector
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: gsq_jsvd_u,                  &
     &   gsq_jsvd_vt
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: gsq_jsvd_w

!  Other matrices and vectors
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: gsq_ata
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: gsq_beta

!  Steps, lists of things related to steps
!  Use zero row for Steepest Descent step
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: gsq_delta_a_sd
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: gsq_delta_a_svd
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: gsq_g2exp, gsq_dg2exp,         &
     &   gsq_dg2exp_lin, gsq_dg2exp_quad, gsq_delta_a_len, gsq_exp_eff,        &
     &   gsq_marg_exp_eff

!*******************************************************************************
! SECTION II. DERIVED-TYPE DECLARATIONS
!*******************************************************************************

!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************

      CONTAINS
!*******************************************************************************
! SECTION IV. GSQ INITIALIZATION
!*******************************************************************************
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      SUBROUTINE gsq_initialize
!  Subroutine to initialize variables and arrays that will be needed.
!  Size determination comes from variables in v3f_global
      
      USE v3f_global, ONLY: n_s_desc, n_rparam

      IMPLICIT NONE

!  Declare Arguments 

!  Declare local variables
      INTEGER                          :: ier1, ier2
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'gsq_initialize: '

!  Start of executable code

!  Define array lengths
      gsq_nsig = n_s_desc
      gsq_nrp = n_rparam

! If allocated, deallocate.
      IF (ALLOCATED(gsq_sweight)) THEN
         DEALLOCATE(gsq_sweight,gsq_ssigma,gsq_ppi,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'deallocation problem 1')
      ENDIF

      IF (ALLOCATED(gsq_evector)) THEN
         DEALLOCATE(gsq_evector,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'deallocation problem 2')
      ENDIF

      IF (ALLOCATED(gsq_jacobian_observe)) THEN
         DEALLOCATE(gsq_jacobian_observe,gsq_jacobian_model,
     &      gsq_jacobian,gsq_jacobian_norm,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'deallocation problem 3')
      ENDIF

      IF (ALLOCATED(gsq_jsvd_u)) THEN
         DEALLOCATE(gsq_jsvd_u,gsq_jsvd_vt,gsq_jsvd_w,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'deallocation problem 4')
      ENDIF

      IF (ALLOCATED(gsq_ata)) THEN
         DEALLOCATE(gsq_ata,gsq_beta,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'deallocation problem 5')
      ENDIF

! Allocate
!  _ns size arrays
      IF (gsq_nsig .gt. 0) THEN
         ALLOCATE(gsq_sweight(gsq_nsig),gsq_ssigma(gsq_nsig),                  &
     &      gsq_evector(gsq_nsig),STAT=ier1)  
         CALL assert_eq(0,ier1,sub_name // 'allocation problem 1')
       ELSE
          CALL err_fatal(sub_name // 'gsq_nsig .le. 0',int=gsq_nsig)
       ENDIF

!  _np size array
      IF (gsq_nrp .gt. 0) THEN
         ALLOCATE(gsq_ppi(gsq_nrp),STAT=ier1)  
         CALL assert_eq(0,ier1,sub_name // 'allocation problem 2')
      ENDIF

!  Two d arrays - jacobian
!  Note that both u and v are square, while w is a vector
      IF ((gsq_nsig .gt. 0) .and. (gsq_nrp .gt. 0)) THEN
         gsq_nsv = MIN(gsq_nsig,gsq_nrp)
         ALLOCATE(gsq_jacobian_observe(gsq_nsig,gsq_nrp),                      &
     &      gsq_jacobian_model(gsq_nsig,gsq_nrp),                              &
     &      gsq_jacobian(gsq_nsig,gsq_nrp),                                    &
     &      gsq_jacobian_norm(gsq_nsig,gsq_nrp),STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'allocation problem 3')
         ALLOCATE(gsq_jsvd_u(gsq_nsig,gsq_nsig),                               &
     &      gsq_jsvd_vt(gsq_nrp,gsq_nrp),                                      &
     &      gsq_jsvd_w(gsq_nsv),STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'allocation problem 4')
      ENDIF

!  Other arrays
      IF (gsq_nrp .gt. 0) THEN
         ALLOCATE(gsq_ata(gsq_nrp,gsq_nrp),STAT=ier1)
         ALLOCATE(gsq_beta(gsq_nrp),STAT=ier2)
         CALL assert_eq(0,ier1,ier2,sub_name // 'allocation problem 5')
      ENDIF

!  Arrays related to various steps
      IF (gsq_nrp .gt. 0) THEN
         ALLOCATE(gsq_delta_a_sd(gsq_nrp),STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'allocation problem 6')
      ENDIF

!  Note zero row used for Steepest Descent step.
      IF ((gsq_nrp .gt. 0) .and. (gsq_nsv .gt. 0)) THEN
         ALLOCATE(gsq_delta_a_svd(0:gsq_nsv,gsq_nrp),STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'allocation problem 7')
      ENDIF

      IF (gsq_nsv .gt. 0) THEN
         ALLOCATE(gsq_g2exp(0:gsq_nsv), gsq_dg2exp(0:gsq_nsv),                 &
     &      gsq_dg2exp_lin(0:gsq_nsv), gsq_dg2exp_quad(0:gsq_nsv),             &
     &      gsq_delta_a_len(0:gsq_nsv), gsq_exp_eff(0:gsq_nsv),                &
     &      gsq_marg_exp_eff(0:gsq_nsv), STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'allocation problem 8')
      ENDIF

!  Initialize the weights to 1
      gsq_sweight = one
      
      RETURN
      END SUBROUTINE gsq_initialize
      
!*******************************************************************************
! SECTION V. GSQ EVALUATION
!*******************************************************************************
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      SUBROUTINE gsq_evaluate_g2(sdm_a_arg,sdo_a_arg)

      IMPLICIT NONE

!  Declare Arguments
!  sdm_a_arg   Signal Data, Model - Array - ARGument
!  sdo_a_arg   Signal Data, Observe - Array - ARGument

      TYPE (signal_data), DIMENSION(:), INTENT(in) :: sdm_a_arg,               &
     &   sdo_a_arg

!  Declare local variables
      INTEGER                          :: nsm, nso, i
      LOGICAL                                 :: l_sigma_problem
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'gsq_evaluate_g2: '

!  Start of executable code

!  Check Allocation, array lengths
      nsm = SIZE(sdm_a_arg)
      nso = SIZE(sdo_a_arg)
      IF (.not. ALLOCATED(gsq_sweight)) THEN
         CALL gsq_initialize
      ENDIF
      
      CALL assert_eq(gsq_nsig,nsm,nso,sub_name // 'array lengths off')

!  Define sigma to use - sqrt(sigma_observe **2 + sigma_model**2)

      l_sigma_problem = .FALSE.
      DO i = 1,gsq_nsig
         gsq_ssigma(i) = SQRT(sdm_a_arg(i) % sigma(1) ** 2 +                   &
     &      sdo_a_arg(i) % sigma(1) ** 2)
         IF (gsq_ssigma(i) .le. zero) l_sigma_problem = .TRUE.
      END DO
      IF (l_sigma_problem) THEN
         CALL err_warn(sub_name // 'sigma problem')
         DO i = 1,gsq_nsig
            WRITE(*,*) i, sdm_a_arg(i) % sigma(1),                             &
     &         sdo_a_arg(i) % sigma(1)
         END DO
         CALL err_fatal(sub_name // 'sigma problem')
      END IF

!  Compute Function
      DO i = 1,gsq_nsig
         gsq_evector(i) = SQRT(gsq_sweight(i)) * ((sdo_a_arg(i) %              &
     &      data(1) - sdm_a_arg(i) % data(1)) / gsq_ssigma(i))
      END DO
      gsq_g2 = dot_product(gsq_evector,gsq_evector)
      
      RETURN
      END SUBROUTINE gsq_evaluate_g2

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      SUBROUTINE gsq_evaluate_jac(model_arg,sdm_a_arg,sdo_a_arg,               &
     &   rparam_a_arg,rcons)

!  Modules USEd
      USE signal_mc
      USE eq_interface
      USE v3f_global
      USE recon_param_model

      IMPLICIT NONE

!  Declare Arguments
      TYPE (model), INTENT (inout), TARGET :: model_arg
      TYPE (signal_data), DIMENSION(:), INTENT(inout) :: sdm_a_arg,            &
     &   sdo_a_arg
      TYPE (recon_param), DIMENSION(:), INTENT(inout) :: rparam_a_arg
      TYPE (recon_cnstrnts), INTENT(inout) :: rcons

!  Declare local variables
      TYPE(eq_state), POINTER :: state_vec => null()
      TYPE(eq_state) :: eq_state_0
      TYPE (recon_param), DIMENSION(:), ALLOCATABLE :: rp_initial_a
      TYPE (recon_param) :: rp_delta
      TYPE (signal_data), DIMENSION(:), ALLOCATABLE :: sdm_temp
      INTEGER  :: numiter, iter_v, isig, irp
      INTEGER :: ier1, ier2
      INTEGER  :: nsig, nrp, nrow, ncol, minrowcol, l_work_svd,                &
     &   info_svd
      LOGICAL :: lconv
      REAL(rprec) :: fn, delta_rp, g2_initial
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: work_svd
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: evector_initial
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'gsq_evaluate_jac: '

      INTEGER :: nda, ndr, nw, iw, ir, k, icount
      INTEGER :: ier_flag_vmec_local, delta_iter2_local
      CHARACTER(len=80) vmec_flag_name_local
      REAL(rprec) :: fd_factor
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: temp_np, temp_ns,              &
     &   wi_temp
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: temp_jac
      INTEGER, DIMENSION(2) :: dimlens

!  Start of executable code

!  Initialization stage
!-------------------------------------------------------------------------------
!  Get the reconstruction parameters of the State
      state_vec => model_arg % eqstate
      nrp = SIZE(rparam_a_arg)
      ALLOCATE(rp_initial_a(nrp),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc rp_initial_a')
      rp_initial_a = rparam_a_arg
      CALL recon_param_model_get(rp_initial_a,model_arg)
!
!  Iterate the state to equilibrium.
      CALL eq_history_set(i2=0)
      numiter = state_vec % fixp % niter
      CALL eq_step(numiter,state_vec,lconverged = lconv,                       &
     &   iter_es = iter_v)
      CALL assert(lconv,sub_name // ' no initial convergence')
      eq_state_0 = state_vec

!  Compute the model signal data. Assumes that sdm_a_arg has been created.
      nsig = SIZE(sdm_a_arg)
      CALL assert_eq(nsig,n_s_desc,sub_name // 'sdm_a_arg size')
      ALLOCATE(sdm_temp(nsig),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'sdm_temp alloc')
!  (need call for array)
      DO isig = 1,nsig
         CALL signal_mc_model_compute(sdm_a_arg(isig),s_desc(isig),            &
     &      model_arg)
      END DO

!  Evaluate the g2 function. (Will call _initialize if needed) 
      CALL gsq_evaluate_g2(sdm_a_arg,sdo_a_arg)

!  Allocate space for the evector temporary array
      IF (ALLOCATED(evector_initial)) THEN
         DEALLOCATE(evector_initial,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc evector_initial')
      ENDIF
      ALLOCATE(evector_initial(gsq_nsig),STAT=ier1)  
      CALL assert_eq(0,ier1,sub_name // 'alloc evector_initial')

!  Check array sizes
      CALL assert_eq(gsq_nsig,nsig,sub_name // 'array nsig ?')
      CALL assert_eq(gsq_nrp,nrp,sub_name // 'array nrp ?')

!  Define the pphi values. _ppi allocated in _initialize
      DO irp = 1,nrp
         gsq_ppi(irp) = rparam_a_arg(irp) % vrnc
         IF (gsq_ppi(irp) .le. zero) THEN
            gsq_ppi(irp) = one
         ENDIF
      END DO

!  set Jacobian arrays to zero
      gsq_jacobian_observe = zero
      gsq_jacobian_model = zero
      gsq_jacobian = zero
      gsq_jacobian_norm = zero
      gsq_jsvd_u = zero
      gsq_jsvd_vt = zero
      gsq_jsvd_w = zero

!  Jacobian for observe signals
!-------------------------------------------------------------------------------
!  09-28-06. Leave gsq_jacobian_observe = zero. 

!  Finite Difference stage, for model signals
!-------------------------------------------------------------------------------
!  Save initial g2 value, and initial gsq_evector
      g2_initial = gsq_g2
      evector_initial = gsq_evector
! [JDH commented out below 2008-08-08]
!      WRITE(59,*)
!      WRITE(59,*) ' In gsq_evaluate_jac'
!      WRITE(59,*) ' irp   icount   delta_iter2   fd_factor'
      DO irp = 1,nrp   ! loop over the reconstruction parameters
         icount = 0   !  Count the number of attempts with this parameter
100      icount = icount + 1

!  Put in a test to see if the count is too high. If so, neutralize the recon_param
               
!  Put a test to see if this recon_param has been temporarily neutralized
!    If so, set the jacobian column to zero, and cycle.

!  Start from the initial reconstruction parameters, and increment the
!  parameter for the finite difference. Note that the recon_param keeps track 
!  of its own finite difference increment.
         rp_delta = rp_initial_a(irp)
         CALL recon_param_fd_increment(rp_delta,delta_rp)

!  JDH 2010-11-23. Switch Order. Before, call to recon_param_model_put was
!    before eq_change_vmi_cp_xc
!    Should be eq_change_vmi_cp_xc before recon_param_model_put because
!    eq_change_vmi_cp_xc also changes internal state variables.
!  Restore eq_state_0 xc vector to the VMEC Internal State
         CALL eq_change_vmi_cp_xc(eq_state_0)

!  Propagate this change to the VMEC Internal State
         CALL recon_param_model_put(rp_delta,model_arg,rcons,'VMI')

!  Try to equilibrate
         CALL eq_history_set(i2 = irp)
         numiter = state_vec % fixp % niter
         CALL eq_step(numiter,state_vec,lconverged = lconv,                    &
     &      iter_es = iter_v,ier_flag_vmec=ier_flag_vmec_local,                &
     &      vmec_flag_name=vmec_flag_name_local,                               &
     &      delta_iter2=delta_iter2_local)

!  Check for convergence
         IF (.not. lconv) THEN
            CALL eq_history_print
            CALL assert(lconv,sub_name // ' no convergence, irp')
         ENDIF

!  Coding below commented out JDH 2008-08-08
!         IF ((delta_iter2_local .le. 1) .AND. (icount .le. 3)) THEN
!            CALL recon_param_fdf_get(rp_initial_a(irp),fd_factor)
!            fd_factor = 3 * fd_factor
!            WRITE(*,*) '     jac: irp, fd_factor:', irp, fd_factor
!            CALL recon_param_fdf_put(rp_initial_a(irp),fd_factor)
!            GO TO 100
!         END IF
            
!  Compute model signals at displaced parameters
         DO isig = 1,nsig
            CALL signal_mc_model_compute(sdm_temp(isig),s_desc(isig),          &
     &         model_arg)
         END DO

!  Update the jacobian
         DO isig = 1,nsig
            gsq_jacobian_model(isig,irp) = (sdm_temp(isig) % data(1) -         &
     &            sdm_a_arg(isig) % data(1)) / delta_rp
         END DO
         
!  Compute displaced gsq
         CALL gsq_evaluate_g2(sdm_temp,sdo_a_arg)

!  Restore the initial parameter
         CALL recon_param_model_put(rp_initial_a(irp),model_arg,rcons,         &
     &      'VMI')

!  print out to 59 for diagnostic [JDH commented out below 2008-08-08]
!         CALL recon_param_fdf_get(rp_delta,fd_factor)
!         WRITE(59,1000) irp, icount, delta_iter2_local, fd_factor
      END DO ! loop over the reconstruction parameters

!  Restore initial g2 and evector values. (reconstruction parameter argument is unchanged)
      gsq_g2 = g2_initial
      gsq_evector = evector_initial
      
!  Calculation of other jacobian matrices
      
      gsq_jacobian = gsq_jacobian_model -  gsq_jacobian_observe
      DO irp = 1,nrp
         DO isig = 1,nsig
            fn = SQRT(gsq_sweight(isig)) * gsq_ppi(irp) /                      &
     &         gsq_ssigma(isig)
            gsq_jacobian_norm(isig,irp) = gsq_jacobian(isig,irp) * fn
         END DO
      END DO

!-------------------------------------------------------------------------------
!  Singular Value Decomposition of normalized Jacobian
!-------------------------------------------------------------------------------
      nrow = nsig
      ncol = nrp
      minrowcol = min(nrow,ncol)
      CALL assert_eq(gsq_nsv,minrowcol,sub_name // 'array nsv ?')

!  JDH 12-30-2006. Fixed up and uncommented
!  Allocate space for the svd work array
      l_work_svd = 5 * MAX(nrow,ncol)
      IF (ALLOCATED(work_svd)) THEN
         DEALLOCATE(work_svd)
      ENDIF
      ALLOCATE(work_svd(l_work_svd), stat = ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocate work_svd')
!  Allocate space for the temporary jacobian array
      IF (ALLOCATED(temp_jac)) THEN
         DEALLOCATE(temp_jac)
      ENDIF
      dimlens = SHAPE(gsq_jacobian_norm)
      ALLOCATE(temp_jac(dimlens(1),dimlens(2)), stat = ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocate temp_jac')
      temp_jac = gsq_jacobian_norm

      CALL dgesvd('All','All',nrow,ncol,temp_jac,nrow,                         &
     &   gsq_jsvd_w,gsq_jsvd_u,nrow,gsq_jsvd_vt,ncol,work_svd,                 &
     &   l_work_svd,info_svd)
      CALL assert_eq(0,info_svd,sub_name // 'dgesvd problem')

!  JDH 12-30-06. Commented out below. 
!  Use SVDCMP, from LIBSTELL/Sources/SVDpack. (Numerical Recipes + mods)
!      gsq_jsvd_u = gsq_jacobian_norm
!      CALL svdcmp(gsq_jsvd_u,nrow,ncol,nrow,ncol,gsq_jsvd_w,                   &
!     &   gsq_jsvd_vt)
!  I don't think SVDCMP returns sorted singular values. Use SORTSVD from
!  LIBSTELL/Sources/SVDpack
!  JDH 12-30-2006. It appears that sortsvd has a problem when ncol=1
!      CALL sortsvd(nrow,ncol,nrow,ncol,gsq_jsvd_w,gsq_jsvd_u,                  &
!     &   gsq_jsvd_vt)
!  SVDCMP returns V, not Transpose of V, so transpose it.
!      gsq_jsvd_vt = TRANSPOSE(gsq_jsvd_vt)

!  Hessian ( alpha, A-TRANSPOSE . A) and beta vector
      gsq_ata = MATMUL(TRANSPOSE(gsq_jacobian_norm),gsq_jacobian_norm)
      gsq_beta = MATMUL(TRANSPOSE(gsq_jacobian_norm),                          &
     &   gsq_evector)

!-------------------------------------------------------------------------------
!  Various Steps in normalized parameter (a) space
!-------------------------------------------------------------------------------
!  Steepest Descent step, and copy to zero row of svd steps
!      IF (gsq_jsvd_w(1) .gt. zero) THEN
!         gsq_delta_a_sd = gsq_beta / (gsq_jsvd_w(1) ** 2)
!      ENDIF
!  2007-06-28 - Put in minimizing length of SD step
      gsq_delta_a_sd = gsq_beta * DOT_PRODUCT(gsq_beta,gsq_beta) /             &
     &   DOT_PRODUCT(gsq_beta,MATMUL(gsq_ata,gsq_beta))
      gsq_delta_a_svd(0,:) = gsq_delta_a_sd

!  k-SVD steps
!  First, Define a temporary vector for - U-Tranpose dot e
      ALLOCATE(temp_ns(gsq_nsig),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // ' allocate temp_ns')
      temp_ns = MATMUL(TRANSPOSE(gsq_jsvd_u),gsq_evector)

!  Second, Define a temporary vector the inverse singular values
      ALLOCATE(wi_temp(gsq_nsv),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // ' allocate wi')
      DO k = 1,gsq_nsv
         IF (gsq_jsvd_w(k) .gt. zero) THEN
            wi_temp(k) = one / gsq_jsvd_w(k)
         ELSE
            wi_temp(k) = zero
         ENDIF
      END DO

!  Third, Define a temporary vector to hold a gsq_nrp size vector
      ALLOCATE(temp_np(gsq_nrp),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // ' allocate temp_np')

!  Perform pseudo-inversion for successive numbers of singular values retained
      DO k = 1, gsq_nsv
         DO ir = 1,gsq_nrp
            IF (ir .le. k) THEN
               temp_np(ir) = wi_temp(ir) * temp_ns(ir)
            ELSE
               temp_np(ir) = zero
            ENDIF
         END DO
         gsq_delta_a_svd(k,:) = MATMUL(TRANSPOSE(gsq_jsvd_vt),temp_np)
      END DO

!   Calculations related to change in gsq. Use temp_ns for delta_e
      DO k = 0, gsq_nsv
         temp_ns = MATMUL(gsq_jacobian_norm,gsq_delta_a_svd(k,:))
         gsq_dg2exp_lin(k) = - 2 * dot_product(gsq_evector,temp_ns)
         gsq_dg2exp_quad(k) = dot_product(temp_ns,temp_ns)
         gsq_dg2exp(k) = gsq_dg2exp_lin(k) + gsq_dg2exp_quad(k)
         gsq_g2exp(k) = gsq_g2 + gsq_dg2exp(k)
         gsq_delta_a_len(k) = SQRT(DOT_PRODUCT(gsq_delta_a_svd(k,:),           &
     &      gsq_delta_a_svd(k,:)))
         gsq_exp_eff(k) = ABS(gsq_dg2exp_lin(k) + gsq_dg2exp_quad(k))          &
     &      / gsq_delta_a_len(k)
      END DO
!  Although marginal efficiencies are not strictly defined for the
!  Steepest Descent (index 0) and just one singular value cases, for convenience
!  I define them here. JDH 2008-06-04
      gsq_marg_exp_eff(0:1) = gsq_exp_eff(0:1)
      DO k = 2,gsq_nsv
         gsq_marg_exp_eff(k) =                                                 &
     &      (ABS(gsq_dg2exp_lin(k) + gsq_dg2exp_quad(k)) -                     &
     &       ABS(gsq_dg2exp_lin(k-1) + gsq_dg2exp_quad(k-1))) /                & 
     &      (gsq_delta_a_len(k) - gsq_delta_a_len(k-1))
      END DO
      
      RETURN
1000  FORMAT(1x,i3,2x,i2,2x,i5,2x,es10.3)
      END SUBROUTINE gsq_evaluate_jac
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      SUBROUTINE gsq_evaluate_dar(step_cntrl_a_arg,delta_a_arg,                &
     &  delta_r_arg,kuse_arg)
!  Subroutine to evaluate the change in the parameters, using an SVD 
!  pseudo-inverse of the normalized jacobian matrix.
!  Assumes that gsq_evaluate_jac has been called, and all the relevant 
!  gsq_ variables have been appropriately defined.
!  The array of step_cntrl array determines how many singular values to retain,
!  and the step size.

      USE v3f_global, ONLY: iou_runlog

      IMPLICIT NONE

!  Declare Arguments
      REAL(rprec), DIMENSION(6), INTENT(in) :: step_cntrl_a_arg
      REAL(rprec), DIMENSION(:), INTENT(out) :: delta_a_arg
      REAL(rprec), DIMENSION(:), INTENT(out) :: delta_r_arg
      INTEGER, INTENT(out) :: kuse_arg

!  Declare local variables
      INTEGER :: kuse, iw, k
      REAL(rprec) :: w_one, fracstep
      REAL(rprec) :: cut_svd, cut_eff, cut_marg_eff, cut_delta_a,              &
     &  cut_dg2, astep_max
      LOGICAL, DIMENSION(gsq_nsv) :: l_svd, l_eff, l_marg_eff,                 &
     &  l_delta_a, l_dg2, l_all
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'gsq_evaluate_dar: '

!  Start of executable code

!  Algorithm to choose step. Uses lots of cutoffs
      cut_svd = step_cntrl_a_arg(1)
      cut_eff = step_cntrl_a_arg(2)
      cut_marg_eff = step_cntrl_a_arg(3)
      cut_delta_a = step_cntrl_a_arg(4)
      cut_dg2 = step_cntrl_a_arg(5)
      astep_max = step_cntrl_a_arg(6)
!      WRITE(*,*) 'cut_svd  cut_eff  cut_marg_eff  cut_delta_a  cut_dg2'
!      WRITE(*,300) cut_svd, cut_eff, cut_marg_eff, cut_delta_a, cut_dg2
!      WRITE(*,*) 'astep_max = ', astep_max
      WRITE(iou_runlog,*)                                                      & 
     &   'cut_svd  cut_eff  cut_marg_eff  cut_delta_a  cut_dg2'
      WRITE(iou_runlog,300)                                                    & 
     &   cut_svd, cut_eff, cut_marg_eff, cut_delta_a, cut_dg2
      WRITE(iou_runlog,*) 'astep_max = ', astep_max
300   FORMAT(5(2x,es9.2))

      w_one = gsq_jsvd_w(1)
      IF (w_one .le. zero) THEN
         w_one = one
      ENDIF
      l_svd = gsq_jsvd_w / w_one .GE. cut_svd
      l_eff = gsq_exp_eff(1:gsq_nsv) .GE. cut_eff
      l_marg_eff = gsq_marg_exp_eff(1:gsq_nsv) .GE. cut_marg_eff
      l_delta_a = gsq_delta_a_len(1:gsq_nsv) .GE. cut_delta_a
      l_dg2 = ABS(gsq_dg2exp_lin(1:gsq_nsv) +                                  &
     &   gsq_dg2exp_quad(1:gsq_nsv)).GE. cut_dg2
      l_all = l_svd .AND. l_eff .AND. l_marg_eff .AND. l_delta_a .AND.         &
     &   l_dg2

      kuse = 0
!      WRITE(*,*) ' k  l_svd  l_eff l_marg_eff l_delta_a l_dg2 l_all'
      WRITE(iou_runlog,*)                                                      & 
     &   ' k  l_svd  l_eff l_marg_eff l_delta_a l_dg2 l_all'
      DO k = 1,gsq_nsv
!         WRITE(*,320) k, l_svd(k), l_eff(k), l_marg_eff(k),                    &
         WRITE(iou_runlog,320) k, l_svd(k), l_eff(k), l_marg_eff(k),           &
     &      l_delta_a(k), l_dg2(k), l_all(k)
320   FORMAT(i3,4x,6(l1,7x))
         IF (l_all(k))  kuse = k
      END DO
      
      WRITE(*,*) 'Number of singular values used for step is ',kuse
      WRITE(iou_runlog,*)                                                      & 
     &   'Number of singular values used for step is ',kuse
      kuse_arg = kuse

!  Check size of step
      fracstep = one
      IF ((astep_max .GE. 0.) .AND.                                            &
     &    (gsq_delta_a_len(kuse) .GT. astep_max)) THEN
        fracstep = astep_max / gsq_delta_a_len(kuse)
        WRITE(*,*) ' Step Size Limited to ', astep_max
        WRITE(*,*) ' Expected g2 = ', gsq_g2 + gsq_dg2exp(kuse) *              &
     &     (fracstep * (2. - fracstep))
        WRITE(iou_runlog,*) ' Step Size Limited to ', astep_max
        WRITE(iou_runlog,*) ' Expected g2 = ',                                 &
     &     gsq_g2 + gsq_dg2exp(kuse) * (fracstep * (2. - fracstep))
      ENDIF
      
      delta_a_arg = fracstep * gsq_delta_a_svd(kuse,:)

!  Unnormalize parameters
      delta_r_arg = delta_a_arg * gsq_ppi

      RETURN
      END SUBROUTINE gsq_evaluate_dar

*******************************************************************************
! SECTION XVI.  OUTPUT SUBROUTINES
!*******************************************************************************

      SUBROUTINE gsq_write_jac(identifier,unit,verbose)
      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in), OPTIONAL :: identifier
      INTEGER, INTENT(in), OPTIONAL :: unit
      INTEGER, INTENT(in), OPTIONAL :: verbose
!  identifier   character variable, also written out
!  unit         I/O unit number to write to
!  verbose      integer, to specify verbosity level of write

!  Declare local variables and constants
      INTEGER      :: iv_default = 1
      INTEGER      :: iv
      INTEGER      :: iou_default = 6
      INTEGER      :: iou
      CHARACTER (len=60)  :: id
      INTEGER      :: i, icol_per_page, kmax, k, jmin, jmax, j
      REAL(rprec)        :: w_one, temp
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'gsq_write_jac: '
      
!  Do Formats with statement labels. Format Arrays are too clumsy to modify.
200   FORMAT(/'**** Subroutine GSQ_WRITE_JAC  START  ******')
210   FORMAT(' Called with id = ',a)
220   FORMAT(' gsq_nsig, gsq_nrp, gsq_nsv = ',i5,2x,i3,2x,i3)
230   FORMAT(' g squared = ',es12.5)
240   FORMAT(/' Jacobian (normalized) SVD U(i,j), i down')
250   FORMAT(5x,20(2x,i5,2x))
260   FORMAT(i5,20(1x,f8.5))
270   FORMAT(/' Jacobian (normalized) SVD VT(i,j), i down')
273   FORMAT(/' Jacobian (normalized) A(i,j), i down')
276   FORMAT(/' Jacobian (model - observe) J(i,j), i down')
280   FORMAT(/' Signals:')
290   FORMAT(8x,'Weight',8x,'Sigma',9x,'evector',7x,'ev-norm') 
300   FORMAT(i5,20(2x,es12.5))
310   FORMAT(/' Various Steps. gsq_g2 now = ',es12.5,/                         &
     &   4x,'(dg2: lin = -2 * dg2, quad = +1 * dg2)'/)
320   FORMAT(19x,'Relative',38x,'Efficiency  Mrgnl Eff'/                       &
     &   9x,'SV',12x,'SV    g2-expect',5x,'dg2',7x,'d_|a|',                   &
     &   6x,'|dg2|/d_a',3x,'ddg2/dd_a')
!320   FORMAT(9x,'SV',7x,'Relative SV  g2-expect',5x,'dg2',7x,                  &
!     &   'd_|a|')
330   FORMAT('StDe',6x,'-----',6x,'-----',3(2x,es10.3),2(2x,es9.2))
340   FORMAT(i3,2x,es10.3,2x,es9.2,3(2x,es10.3),2(2x,es9.2))
!350   FORMAT(/7x,'Efficiency  Marginal Eff'/                                   &
!     &   7x,'dg2/d_a',5x,'ddelg2/ddel_a')
!360   FORMAT('StDe',2(2x,es10.3))
!365   FORMAT(i4,2x,es10.3,7x,'-----')
!370   FORMAT(i4,20(2x,es10.3))
380   FORMAT(/'**** Subroutine GSQ_WRITE_JAC  END  ******'/)

!  start of executable code
!  Check for arguments present
      IF (PRESENT(identifier)) THEN
         id = identifier
      ELSE
         id = ' '
      END IF

      IF (PRESENT(unit)) THEN
         iou = unit
      ELSE
         iou = iou_default
      END IF

      IF (PRESENT(verbose)) THEN
         iv = verbose
      ELSE
         iv = iv_default
      END IF

!  Identifier, Integers, g2
!  All verbosity levels
      WRITE(iou,200)
      WRITE(iou,210) id
      WRITE(iou,220) gsq_nsig, gsq_nrp, gsq_nsv
      WRITE(iou,230) gsq_g2
      
!  Jacobian SVD matrices
!  Only verbosity >= 3
      IF (iv .GE. 3) THEN

!  U first, nsignal by nsignal. k loops over "pages" of printout
!  Only print out nsv columns. Others are irrelevant
         icol_per_page = 7
         kmax = (gsq_nsv - 1) / icol_per_page + 1
         WRITE(iou,240)
         DO k = 1,kmax
            jmin = (k - 1) * icol_per_page + 1
            jmax = MIN(jmin + icol_per_page - 1,gsq_nsv)
            WRITE(iou,250) (j,j=jmin,jmax) 
            DO i = 1,gsq_nsig
               WRITE(iou,260) i, (gsq_jsvd_u(i,j),j=jmin,jmax)
            END DO
         END DO

!  V-Transpose, nrp by nrp. k loops over "pages" of printout
         icol_per_page = 7
         kmax = (gsq_nrp - 1) / icol_per_page + 1
         WRITE(iou,270)
         DO k = 1,kmax
            jmin = (k - 1) * icol_per_page + 1
            jmax = MIN(jmin + icol_per_page - 1,gsq_nrp)
            WRITE(iou,250) (j,j=jmin,jmax)
            DO i = 1,gsq_nrp
               WRITE(iou,260) i, (gsq_jsvd_vt(i,j),j=jmin,jmax)
            END DO
         END DO
      END IF

!  Normalized Jacobian, nsignal by nrp 
!  Only verbosity >= 3
      IF (iv .GE. 3) THEN
         icol_per_page = 7
         kmax = (gsq_nrp - 1) / icol_per_page + 1
         WRITE(iou,273)
         DO k = 1,kmax
            jmin = (k - 1) * icol_per_page + 1
            jmax = MIN(jmin + icol_per_page - 1,gsq_nrp)
            WRITE(iou,250) (j,j=jmin,jmax) 
            DO i = 1,gsq_nsig
               WRITE(iou,260) i, (gsq_jacobian_norm(i,j),j=jmin,jmax)
            END DO
         END DO
      END IF

!  Jacobian, nsignal by nrp 
!  Only verbosity >= 2
      IF (iv .GE. 2) THEN
         icol_per_page = 7
         kmax = (gsq_nrp - 1) / icol_per_page + 1
         WRITE(iou,276)
         DO k = 1,kmax
            jmin = (k - 1) * icol_per_page + 1
            jmax = MIN(jmin + icol_per_page - 1,gsq_nrp)
            WRITE(iou,250) (j,j=jmin,jmax) 
            DO i = 1,gsq_nsig
               WRITE(iou,260) i, (gsq_jacobian(i,j),j=jmin,jmax)
            END DO
         END DO
      END IF
      
!  Vectors of length gsq_nsig
!  Only verbosity .GE. 1
      IF (iv .GE. 1) THEN
         WRITE(iou,280)
         WRITE(iou,290)
         temp = MAXVAL(ABS(gsq_evector(1:gsq_nsig)))
         DO i = 1,gsq_nsig
            WRITE(iou,300) i, gsq_sweight(i), gsq_ssigma(i),                   &
     &         gsq_evector(i), gsq_evector(i) / temp
         END DO
      END IF

!  Various Steps
!  Only verbosity .GE. 1
      IF (iv .GE. 1) THEN
         WRITE(iou,310) gsq_g2
         WRITE(iou,320)
         WRITE(iou,330) gsq_g2exp(0), -gsq_dg2exp_quad(0),                     &
     &         gsq_delta_a_len(0), gsq_exp_eff(0), gsq_marg_exp_eff(0)
         w_one = gsq_jsvd_w(1)
         IF (w_one .le. zero) THEN
            w_one = one
         ENDIF
         DO i = 1,gsq_nsv
            WRITE(iou,340) i, gsq_jsvd_w(i), gsq_jsvd_w(i) / w_one,            &
     &         gsq_g2exp(i), -gsq_dg2exp_quad(i),                              &
     &         gsq_delta_a_len(i), gsq_exp_eff(i), gsq_marg_exp_eff(i)
         END DO
!         WRITE(iou,350)
!         WRITE(iou,360) gsq_exp_eff(0), gsq_marg_exp_eff(0)
!         DO i = 1,gsq_nsv
!            WRITE(iou,370) i, gsq_exp_eff(i), gsq_marg_exp_eff(i)                                &  
!         END DO
      END IF
      WRITE(iou,380)

      RETURN
      END SUBROUTINE gsq_write_jac
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      SUBROUTINE gsq_write_conf(rpa,identifier,unit,verbose)
      IMPLICIT NONE

!  Subroutine to write out information about the confidence limits on the
!  individual reconstruction parameters
!  Assumes that the subroutine gsq_evaluate_jac has been called.
!  JDH 2007-06-28 First Version
!  JDH 2009-10-06 Significant Modifications

!  Declare Arguments

      TYPE(recon_param), DIMENSION(:), INTENT(in) :: rpa
      CHARACTER (len=*), INTENT(in), OPTIONAL :: identifier
      INTEGER, INTENT(in), OPTIONAL :: unit
      INTEGER, INTENT(in), OPTIONAL :: verbose
!  rpa          type recon_param, array of parameter values
!  identifier   character variable, also written out
!  unit         I/O unit number to write to
!  verbose      integer, to specify verbosity level of write

!  Declare local variables and constants
      INTEGER      :: iv_default = 1
      INTEGER      :: iv
      INTEGER      :: iou_default = 6
      INTEGER      :: iou
      CHARACTER (len=60)  :: id
      INTEGER      :: i, icol_per_page, kmax, k, jmin, jmax, j
      INTEGER      :: mrp, ksvr, lsv, irp, jsv, ier1
      INTEGER      :: ier2, ier3
      REAL(rprec)         :: w_one, temp
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: gsq_1pcl
      
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: gsq_b_1pcl,gsq_b_pace
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: gsq_b, gsq_b_u,              &
     &  gsq_b_vt
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: gsq_b_w

      INTEGER  :: nrow, ncol, minrowcol, l_work_svd,                    &
     &   info_svd
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: work_svd
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: temp_jac
      INTEGER, DIMENSION(2) :: dimlens

      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: sem_a
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: psigma_v

      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'gsq_write_conf: '
      
!  Do Formats with statement labels. Format Arrays are too clumsy to modify.
200   FORMAT(/'**** Subroutine GSQ_WRITE_CONF  START  ******')
210   FORMAT(' Called with id = ',a)
220   FORMAT(' gsq_nsig, gsq_nrp, gsq_nsv = ',i5,2x,i3,2x,i3)
240   FORMAT(/' One-parameter 68% confidence limits.',a/                       &
     &   '    Down - reconstruction parameter'/                                &
     &   '    Across - Number singular values retained')
250   FORMAT(a,t36,20(3x,i5,4x))
260   FORMAT(a,t35,20(2x,es10.3))
270   FORMAT(a,a)
280   FORMAT(i4,2x,a10,2x,i4,20(2x,es10.3))

290   FORMAT(/' Signal Effectiveness Matrix',                                  &
     &   ' partial log(psigma) wrt log(ssigma)'/,                              &
     &   '    Down - reconstruction parameter',                                &
     &   '    Across - signal number')
300   FORMAT(/' Principal Axes of the Delta-g^2=1 ellipsoid',                  &
     &   ' in (unnormalized) parameter space'/,                                &
     &   '    Down - reconstruction parameter',                                &
     &   '    Across - Principal Axis Number (ordered by SV)')
310   FORMAT(a,t36,20(3x,i5,4x))
320   FORMAT(a)
330   FORMAT(i4,2x,a10,2x,i4,20(2x,es10.3))

350   FORMAT(/' Signal Effectiveness Matrix - Each Parameter',/                &
     &   ' irp  %p_type   %index    %value     psigma',                        &
     &   '    SUM    MAX   MAXLOC (Signal index)')
354   FORMAT(i4,2x,a10,2x,i4,4(2x,es10.3),2x,i4)

360   FORMAT(/' Signal Effectiveness Matrix - Each Signal',/                   &
     &   ' isig  SUM    MAX   MAXLOC (Parameter index)') 
364   FORMAT(i4,2x,2(2x,es10.3),2x,i4)

400   FORMAT(/'**** Subroutine GSQ_WRITE_CONF  END  ******'/)

!  start of executable code
!  Check for arguments present
      IF (PRESENT(identifier)) THEN
         id = identifier
      ELSE
         id = ' '
      END IF

      IF (PRESENT(unit)) THEN
         iou = unit
      ELSE
         iou = iou_default
      END IF

      IF (PRESENT(verbose)) THEN
         iv = verbose
      ELSE
         iv = iv_default
      END IF

!  Identifier, Integers, g2
      WRITE(iou,200)
      WRITE(iou,210) id
      WRITE(iou,220) gsq_nsig, gsq_nrp, gsq_nsv
      
!  Form the matrix of one-parameter confidence limits, gsq_1pcl
      ALLOCATE(gsq_1pcl(gsq_nrp,gsq_nsv),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'allocate gsq_1pcl')
      DO mrp = 1,gsq_nrp
         DO ksvr = 1,gsq_nsv
            temp = zero
            DO lsv = 1,ksvr
               temp = temp + (gsq_jsvd_vt(lsv,mrp) /                           &
     &                        gsq_jsvd_w(lsv)) ** 2
            END DO
            gsq_1pcl(mrp,ksvr) = gsq_ppi(mrp) * SQRT(temp)
         END DO
      END DO
!  JDH 2009-10-06. Only want the last column of gsq_1pcl.
!  Delay printout

!  Now generate the partially normalized (B) Jacobian.
!  Weights and Signal sigmas are used, but the parameter pi's are not used
!  Generated from the normalized Jacobian using PI(-1)
!  Allocate space for _b_ arrays
      ALLOCATE(gsq_b_1pcl(gsq_nrp,gsq_nsv),STAT=ier1)
      ALLOCATE(gsq_b_pace(gsq_nrp,gsq_nsv),STAT=ier2)
      ALLOCATE(gsq_b(gsq_nsig,gsq_nrp),STAT=ier3)
      CALL assert_eq(0,ier1,ier2,ier3,sub_name // 'allocation-b')
      ALLOCATE(gsq_b_u(gsq_nsig,gsq_nsig),                                     &
     &      gsq_b_vt(gsq_nrp,gsq_nrp),                                         &
     &      gsq_b_w(gsq_nsv),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'allocation-bsvd')

!  Generate B array
      DO i = 1,gsq_nsig
         DO j = 1,gsq_nrp
            gsq_b(i,j) = gsq_jacobian_norm(i,j) / gsq_ppi(j)
         END DO
      END DO

!  Perform SVD on the B array
      nrow = gsq_nsig
      ncol = gsq_nrp
      minrowcol = min(nrow,ncol)
      CALL assert_eq(gsq_nsv,minrowcol,sub_name // 'array nsv ?')

!  Allocate space for the svd work array
      l_work_svd = 5 * MAX(nrow,ncol)
      IF (ALLOCATED(work_svd)) THEN
         DEALLOCATE(work_svd)
      ENDIF
      ALLOCATE(work_svd(l_work_svd), stat = ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocate work_svd')
!  Allocate space for the temporary jacobian array
      IF (ALLOCATED(temp_jac)) THEN
         DEALLOCATE(temp_jac)
      ENDIF
      dimlens = SHAPE(gsq_b)
      ALLOCATE(temp_jac(dimlens(1),dimlens(2)), stat = ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocate temp_jac')
      temp_jac = gsq_b

      CALL dgesvd('All','All',nrow,ncol,temp_jac,nrow,                         &
     &   gsq_b_w,gsq_b_u,nrow,gsq_b_vt,ncol,work_svd,                          &
     &   l_work_svd,info_svd)
      CALL assert_eq(0,info_svd,sub_name // 'dgesvd problem')

!  Form the matrix of one-parameter confidence limits, gsq_b_1pcl
      DO mrp = 1,gsq_nrp
         DO ksvr = 1,gsq_nsv
            temp = zero
            DO lsv = 1,ksvr
               temp = temp + (gsq_b_vt(lsv,mrp) /                              &
     &                        gsq_b_w(lsv)) ** 2
            END DO
            gsq_b_1pcl(mrp,ksvr) = SQRT(temp)
         END DO
      END DO
!  JDH 2009-10-06. Only want the last column of gsq_b_1pcl.
!  Delay Printout

!  Now call new subroutine, gsq_sem
      ALLOCATE(psigma_v(gsq_nrp),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'allocate psigma_v')
      ALLOCATE(sem_a(gsq_nrp,gsq_nsig),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'allocate sem_a')
      CALL gsq_sem(gsq_ssigma,gsq_jacobian,psigma_v,sem_a)

!  Write Out Posterior Parameter Variances
!  Computed 3 ways - gsq_1pcl(:,nrp), gsq_b_1pcl(:,nrp), psigma_v(:)
      WRITE(iou,270) ' irp  %p_type   %index    %value    psigma_1',           &
     &   '   psigma_2    psigma_3'
      DO i = 1,gsq_nrp
         WRITE(iou,280) i,rpa(i) % p_type,rpa(i) % index,                      &
     &      rpa(i) % value, gsq_1pcl(i,gsq_nrp), gsq_b_1pcl(i,gsq_nrp),        &
     &      psigma_v(i)
      END DO

!  Write Out Signal Effectiveness Matrix
!  Write out sem_a, gsq_nrp rows by gsq_nsig columns
!  k loops over "pages" of printout
      icol_per_page = 6
      kmax = (gsq_nsig - 1) / icol_per_page + 1
      WRITE(iou,290)
      DO k = 1,kmax
         jmin = (k - 1) * icol_per_page + 1
         jmax = MIN(jmin + icol_per_page - 1,gsq_nsig)
         WRITE(iou,310) ' Signal Number',(j,j=jmin,jmax) 
         WRITE(iou,320) ' irp  %p_type   %index    %value '
         DO i = 1,gsq_nrp
            WRITE(iou,330) i,rpa(i) % p_type,rpa(i) % index,                   &
     &         rpa(i) % value, (sem_a(i,j),j=jmin,jmax)
         END DO
      END DO

!  Write out information about SEM, each parameter
      WRITE(iou,350)
      DO i = 1,gsq_nrp
         WRITE(iou,354) i,rpa(i) % p_type,rpa(i) % index,                      &
     &      rpa(i) % value, psigma_v(i), SUM(sem_a(i,:)),                      &
     &      MAXVAL(sem_a(i,:)), MAXLOC(sem_a(i,:))
      END DO

!  Write out information about SEM, each signal
      WRITE(iou,360)
      DO i = 1,gsq_nsig
         WRITE(iou,364) i, SUM(sem_a(:,i)),                                    &
     &      MAXVAL(sem_a(:,i)), MAXLOC(sem_a(:,i))
      END DO      
      
!  Form the matrix of the principal axes of the Delta-g^2=1 ellipsoid, gsq_b_pace
      DO irp = 1,gsq_nrp
         DO jsv = 1,gsq_nsv
            gsq_b_pace(irp,jsv) = gsq_b_vt(jsv,irp) /
     &         MAX(gsq_b_w(jsv),1.e-20_dp)
         END DO
      END DO

!  Write out gsq_pace, gsq_nrp rows by gsq_nsv columns
!  k loops over "pages" of printout
      icol_per_page = 5
      kmax = (gsq_nsv - 1) / icol_per_page + 1
      WRITE(iou,300)
      DO k = 1,kmax
         jmin = (k - 1) * icol_per_page + 1
         jmax = MIN(jmin + icol_per_page - 1,gsq_nsv)
         WRITE(iou,310) ' Principal Axis Number',(j,j=jmin,jmax) 
         WRITE(iou,320) ' irp  %p_type   %index    %value '
         DO i = 1,gsq_nrp
            WRITE(iou,330) i,rpa(i) % p_type,rpa(i) % index,                   &
     &         rpa(i) % value, (gsq_b_pace(i,j),j=jmin,jmax)
         END DO
      END DO

      WRITE(iou,400)

!  Clean Up
      DEALLOCATE(gsq_1pcl,gsq_b_1pcl,gsq_b_pace,gsq_b, gsq_b_u,                &
     &   gsq_b_vt,gsq_b_w,STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'Deallocate')

      RETURN
      END SUBROUTINE gsq_write_conf
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      SUBROUTINE gsq_sem(ssigma_v,jac_a,psigma_v,sem_a) 
!  Subroutine to compute the posterior sigmas and the Signal Effectiveness Matrix
!  Written to be a stand-alone subroutine, even though it is located in the gsq
!  module.
!  It does use variables from the stell_kinds and stell_constants modules
!  The Singular Value Decomposition is done with the LAPACK dgesvd subroutine
!    (Used for the inversion of the inverse posterior covariance matrix)
!  As of 2009-10-06, the signal covariance matrix is assumed to be diagonal.
!  JDH 2009-10-06

!            Arguments
!   ssigma_v      1d array of signal sigma values, nsig long
!   jac_a         2d jacobian, nsig by nrp
!   psigma_v      1d array of posterior sigma values, nrp long
!   sem_a         2d Signal Effectiveness Matrix, partial derivative
!                  of ln(psigma_v(j)) wrt ln(ssigma_v(i))
!                  nrp by nsig

!          Local 2d Arrays.  tr - transpose      pi - pseudoinverse
!   jac_local     nsig x nrp, local copy of the jacobian. Argument to dgesvd.
!   jac_dot_cp    nsig x nrp, Jacobian dotted with posterior covariance
!   csignal       nsig x nsig, signal covariance matrix - diagonal
!   csignalinv    nsig x nsig, Inverse of signal covariance matrix - diagonal
!   cposterior    nrp x nrp, posterior parameter covariance matrix
!   cpinv         nrp x nrp, inverse of posterior parameter covariance matrix
!   svd_u         nsig x nsig, U matrix from SVD of jacobian. Argument to dgesvd
!   svd_vt        nrp x nrp, VT matrix from SVD of jacobian. Argument to dgesvd
!   wpi           nrp x nsig, pseudoinverse of W matrix

!            Local 1d Arrays
!   svd_w         Min(nsig,nrp), Vector of SV of jacobian. Argument to dgesvd
  

      IMPLICIT NONE

!  Declare Arguments
      REAL(rprec), INTENT(in), DIMENSION(:) :: ssigma_v           
      REAL(rprec), INTENT(in), DIMENSION(:,:) :: jac_a
      REAL(rprec), INTENT(inout), DIMENSION(:) :: psigma_v
      REAL(rprec), INTENT(inout), DIMENSION(:,:) :: sem_a

!  Declare local variables
      INTEGER  :: nsig, nsig1, nsig2, nsig3
      INTEGER  :: nrp, nrp1, nrp2, nrp3
      INTEGER  :: i, j, ier1
      INTEGER  :: nrow, ncol, minrowcol, l_work_svd, info_svd

!      REAL(rprec)  :: tiny=1.D-20
      REAL(rprec)  :: tiny=1.D-40 !  JDH 2010-06-11. Kludge for sxrc test

      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: jac_local,                   &
     &   csignal, cposterior, svd_u, svd_vt, cpinv, wpi,                       &
     &   csignalinv, jac_dot_cp
     
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: svd_w

      REAL(rprec), DIMENSION(:), ALLOCATABLE :: work_svd
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'gsq_sem: '

!  Start of executable code
!-------------------------------------------------------------------------------

!  Check Dimensions
!-------------------------------------------------------------------------------
      nsig1 = SIZE(ssigma_v)
      nsig2 = SIZE(jac_a,1)
      nrp1 = SIZE(jac_a,2)
      nrp2 = SIZE(psigma_v)
      nrp3 = SIZE(sem_a,1)
      nsig3 = SIZE(sem_a,2)
      CALL assert_eq(nsig1,nsig2,nsig3,sub_name // 'nsig different')
      CALL assert_eq(nrp1,nrp2,nrp3,sub_name // 'nrp different')
      nsig = nsig1
      nrp = nrp1
      IF (nsig .lt. nrp) THEN
         CALL err_warn(sub_name ,'nsig < nrp. nsig=',int=nsig)
         CALL err_warn(sub_name ,'nsig < nrp. nrp=',int=nrp)
      ENDIF

!  Allocate space for local arrays. Allocation for SVD arrays is near SVD call
!-------------------------------------------------------------------------------
      ALLOCATE(jac_local(nsig,nrp),jac_dot_cp(nsig,nrp),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc jac_local')
      
      ALLOCATE(csignal(nsig,nsig),csignalinv(nsig,nsig),                       &
     &   cposterior(nrp,nrp),cpinv(nrp,nrp),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc csignal')

!  Generate signal covariance matrix
!-------------------------------------------------------------------------------
      csignal = zero
      csignalinv = zero
      DO i = 1,nsig
         csignal(i,i) = ssigma_v(i) ** 2
         csignalinv(i,i) = one / ssigma_v(i) ** 2
      END DO

!  Generate inverse of posterior covariance matrix
!-------------------------------------------------------------------------------
      jac_local = jac_a
      cpinv = MATMUL(TRANSPOSE(jac_a),MATMUL(csignalinv,jac_a))
      
!-------------------------------------------------------------------------------
!  Singular Value Decomposition of cpinv (posterior covariance inverse)
!-------------------------------------------------------------------------------
      nrow = nrp
      ncol = nrp
      minrowcol = min(nrow,ncol)

      l_work_svd = 5 * MAX(nrow,ncol)
      ALLOCATE(work_svd(l_work_svd), stat = ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocate work_svd')
      
      ALLOCATE(svd_u(nrow,nrow),svd_vt(ncol,ncol),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc svd 2d arrays')

      ALLOCATE(svd_w(minrowcol),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc svd w array')

      CALL dgesvd('All','All',nrow,ncol,cpinv,nrow,                            &
     &   svd_w,svd_u,nrow,svd_vt,ncol,work_svd,                                &
     &   l_work_svd,info_svd)
      CALL assert_eq(0,info_svd,sub_name // 'dgesvd problem')

!  Generate cposterior - pseudoinverse of cpinv, and the posterior sigmas
!-------------------------------------------------------------------------------
      ALLOCATE(wpi(ncol,nrow),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc wpi array')
      wpi = zero
      DO i = 1,minrowcol
         IF (svd_w(i) .lt. tiny) THEN
            CALL err_warn(sub_name ,'SV small',svd_w(i),i)
         ENDIF
         wpi(i,i) = one / MAX(svd_w(i),tiny)
      END DO
      cposterior = MATMUL(TRANSPOSE(svd_vt),                                   &
     &   MATMUL(wpi,TRANSPOSE(svd_u)))

      DO j = 1,nrp
         psigma_v(j) = SQRT(cposterior(j,j))
      END DO

!  Generate SEM matrix
!-------------------------------------------------------------------------------
      jac_dot_cp = MATMUL(jac_a,cposterior)
      DO j = 1,nrp
         DO i = 1,nsig
            sem_a(j,i) = (jac_dot_cp(i,j) /                                    &
     &            (psigma_v(j) * ssigma_v(i))) ** 2
         END DO
      END DO
      
      RETURN

      END SUBROUTINE gsq_sem
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE gsq_write_mr(rpa,identifier,unit,verbose)
      IMPLICIT NONE

!  Subroutine to write out the most-redundant column variables of the
!  normalized jacobian.
!  Uses subroutine most_redundant, from module v3_utilities
!  Assumes that the subroutine gsq_evaluate_jac has been called.

!  Declare Arguments
      TYPE(recon_param), DIMENSION(:), INTENT(in) :: rpa
      CHARACTER (len=*), INTENT(in), OPTIONAL :: identifier
      INTEGER, INTENT(in), OPTIONAL :: unit
      INTEGER, INTENT(in), OPTIONAL :: verbose
!  rpa          type recon_param, array of parameter values
!  identifier   character variable, also written out
!  unit         I/O unit number to write to
!  verbose      integer, to specify verbosity level of write

!  Declare local variables and constants
      INTEGER      :: iv_default = 1
      INTEGER      :: iv
      INTEGER      :: iou_default = 6
      INTEGER      :: iou
      CHARACTER (len=60)  :: id
      INTEGER      :: i, ier1, jce
      INTEGER, DIMENSION(2) :: dimlens
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: amat
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: svrat_a
      INTEGER, DIMENSION(:), ALLOCATABLE :: ncol_a, j_col_elim_a

      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'gsq_write_mr: '
     
!  Do Formats with statement labels. Format Arrays are too clumsy to modify.
200   FORMAT(/'**** Subroutine GSQ_WRITE_MR  START  ******')
210   FORMAT(' Called with id = ',a)
215   FORMAT(' gsq_nrp = 1, so return')
220   FORMAT(' gsq_nsig (rows), gsq_nrp (columns), gsq_nsv = ',                &
     &   i5,2x,i3,2x,i3)
230   FORMAT('i',3x,'col-remain',3x,'svprod-1',3x,                             &
     &   'j_col_elim',3x,'rp_type',3x,'rp_index')
240   FORMAT(i3,4x,i3,3x,es12.3,2x,i3,2x,a10,1x,i3)
250   FORMAT(//' Last Two Reconstruction Parameters',                          &
     &   'j_col_elim',3x,'rp_type',3x,'rp_index')
260   FORMAT(i3,3x,a10,3x,i3)
380   FORMAT(/'**** Subroutine GSQ_WRITE_MR  END  ******'/)

!  start of executable code
!  Check for arguments present
      IF (PRESENT(identifier)) THEN
         id = identifier
      ELSE
         id = ' '
      END IF

      IF (PRESENT(unit)) THEN
         iou = unit
      ELSE
         iou = iou_default
      END IF

      IF (PRESENT(verbose)) THEN
         iv = verbose
      ELSE
         iv = iv_default
      END IF

!  Identifier, Integers, g2
!  All verbosity levels
      WRITE(iou,200)
      WRITE(iou,210) id
      
!  Quit if only one reconstruction parameter
      IF (gsq_nrp .eq. 1) THEN
         WRITE(iou,215)
         WRITE(iou,380)
         RETURN
      ENDIF
      
!  Allocate Space for the arguments to most_redundant
      dimlens = SHAPE(gsq_jacobian_norm)
      CALL assert_eq(gsq_nsig,dimlens(1),sub_name // 'dimlens(1) ??')
      CALL assert_eq(gsq_nrp,dimlens(2),sub_name // 'dimlens(2) ??')
      ALLOCATE(amat(gsq_nsig,gsq_nrp), stat = ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocate amat')
      ALLOCATE(svrat_a(gsq_nrp), stat = ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocate svrat_a')
      ALLOCATE(ncol_a(gsq_nrp), stat = ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocate ncol_a')
      ALLOCATE(j_col_elim_a(gsq_nrp), stat = ier1)
      CALL assert_eq(0,ier1,sub_name // 'Allocate j_col_elim_a')

!  Copy normalized jacobian
      amat = gsq_jacobian_norm

!  Write out number of rows and columns
      WRITE(iou,220) gsq_nsig, gsq_nrp, gsq_nsv

!  Write Column headings
      WRITE(iou,230)

!  Call most_redundant
      CALL most_redundant(amat,ncol_a,svrat_a,j_col_elim_a)

!  Write out
      WRITE(iou,240) 1, ncol_a(1), svrat_a(1), j_col_elim_a(1)
      
      WRITE(iou,240) (i, ncol_a(i), svrat_a(i), j_col_elim_a(i),               &
     &   rpa(j_col_elim_a(i)) % p_type,rpa(j_col_elim_a(i)) % index,           &
     &   i=2,gsq_nrp - 1)

!  write out the last two parameters, appropriately
      WRITE(iou,250)
      
      jce = j_col_elim_a(gsq_nrp)
      WRITE(iou,260) jce, rpa(jce) % p_type,rpa(jce) % index

      jce = svrat_a(gsq_nrp)
      WRITE(iou,260) jce, rpa(jce) % p_type,rpa(jce) % index

      WRITE(iou,380)      

!  Clean Up
      DEALLOCATE(amat,svrat_a,ncol_a,j_col_elim_a,STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'Deallocate')

      RETURN
      END SUBROUTINE gsq_write_mr
      
*******************************************************************************
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!  
!  JDH 09-20-2006
!     First version of module
!
!  JDH 09-27-2006 - 10-03-06
!     Revised structure, added first version of jacobian calculation
!
!  JDH 9 October 2006
!    Added gsq_evaluate_dar
!
!  JDH 10 October 2006
!    Clean up. Fix bug in _evector calculation.
!
!  JDH 4 December 2006
!    Moved SVD calculation to _evaluate_jac
!    Added gsq_jac_write, to write out details of jacobian matrix stuff.
!
!  JDH 18 December 2006
!    Improvements to gsq_jac_write
!
!  JDH 29 December 2006
!    Tweaks to gsq_jac_write
!
!  JDH 2007-06-22
!    SUBROUTINE gsq_evaluate_dar - added optional argument kuse_arg
!
!  JDH 2007-06-24
!    Slight modifications to gsq_jac_write - At minimum, linear change in g2 
!    is always -2 times the quadratic change.
!
!  JDH 2007-06-28
!    Corrected SD step length. Changed sign of definition of Jacobian - see notes
!    "V3FIT Minimization" of today.
!
!  JDH 2007-06-29
!     Added gsq_write_conf - print out of confidence limit information.
!
!  JDH 2007-06-30
!     Fixing problems (conceptual and coding) in gsq_write_conf
!
!  JDH 2007-10-06
!    Added rcons argument to gsq_evaluate_jac - reconstruction constraints
!
!  SPH 2008-01-?? (JDH writing this)
!    INTEGER(iprec) -> INTEGER - always want to use default integers.
!    2008-01-20 JDH eliminated unneeded references to iprec.
!
!  JDH 2008-01-19 - Added => null() to pointer declarations
!
!  JDH 2008-01-20
!    Revise gsq_write_jac
!
!  JDH 2008-03-15 (?) and 2008-05-19 JDH
!    Add gsq_write_mr
!
!  JDH 2008-06-04
!     Revise gsq_evaluate_dar, so that there are more cutoffs considered.
!
!  JDH 2008-08-08
!    Cleaned up logic a bit in loop over reconstruction parameters in
!    gsq_evaluate_jac. Commented out some coding having to do with
!    changing the finite-difference step size - the problem that prompted
!    this coding has been fixed.
!
!  JDH 2009-01-15
!    Changed dimension on step_cntrl_a_arg from 5 to 6. Oops.
!
!  JDH 2009-10-06
!    Significant Changes to gsq_write_conf. Now calls gsq_sem.
!    Added gsq_sem - computes the Signal Effectiveness Matrix.
!    
!  JDH 2009-10-08
!    Continue revisions to gsq_sem and gsq_write_conf

      END MODULE gsq_mod
