
!>  \brief Module for controlling iteration ("time") evolution of MHD convergence sequence.
!!  Updates the CONTRAVARIANT Fourier components of the "velocity" field (displacement X delt)
      MODULE evolution
      USE stel_kinds
      USE stel_constants
      USE shared_data
      USE shared_functions
      USE island_params, mpol=>mpol_i, ntor=>ntor_i, ns=>ns_i,          &
          nfp=>nfp_i, mnmax=>mnmax_i, hs=>hs_i
      USE timer_mod
      USE perturbation, ONLY: eta_factor
      USE descriptor_mod, ONLY: iam, DIAGONALDONE, nprocs,              &
                                in_hess_nfunct, out_hess_nfunct
      USE nscalingtools, ONLY: SKSDBG, TOFU, PARSOLVER, PARFUNCTISL,    &
                 startglobrow, endglobrow, leftproc, rightproc, MPI_ERR
      IMPLICIT NONE

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: n_bsupu_full=1
      INTEGER, PARAMETER :: ndamp = 20              !<number smoothing steps for conj grad
      INTEGER, PARAMETER :: ncongrad_steps = 100    !<max number of conj grad steps
      LOGICAL, PARAMETER :: l_Hess_sym=.FALSE.      !<forces symmetrization of Hessian if TRUE
      LOGICAL            :: l_EPAR                  !<=TRUE when applying parallel tearing perturbations 
      LOGICAL            :: l_Gmres                 !<=TRUE performs GMRES iterations 
      LOGICAL            :: l_conjmin               !<=TRUE if conjugate gradient lowers fsq_min

      INTEGER            :: nprint_step=1
      INTEGER            :: nfun_calls
      INTEGER            :: js, moff, noff

      REAL(rprec)        :: tau_damp, delt_cg=1
      REAL(rprec)        :: fsqprev, fsqprev1
      REAL(rprec)        :: fsq_min, fsq_ratio, fsq_ratio1 
      REAL(rprec)        :: otau(ndamp)
      REAL(rprec), DIMENSION(:), ALLOCATABLE         :: xcdot
            
#if defined(PETSC_OPT)
      integer(kind=selected_int_kind(10)) Apc
      integer(kind=selected_int_kind(10)) subApc(1048)
      integer(kind=selected_int_kind(10)) ksp
      integer(kind=selected_int_kind(10)) subksp(1048)
#endif

#if defined(SKS)
      INTEGER            :: diagnum=1, blocknum=1
      INTEGER            :: EVOLVEPASS=0
#endif
!-----------------------------------------------

      CONTAINS

!>  \brief Performs initial convergence of force residuals using diagonal preconditioner only.
!!  Terminates when the force residual stops changing by more than a factor of 2
      SUBROUTINE converge_diagonal (wout_file, ftol)
      USE perturbation, ONLY: add_perturb, ladd_pert, l_output_alliter,  &   !R. Sanchez, added Jul '10: L_OUTPUT_ALLITER Dump info on file at each timestep
                              lprecon_diag
      USE dumping_mod, ONLY: write_output
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*), INTENT(in)  :: wout_file  !<name of VMEC wout file that SIESTA reads
      REAL(rprec), INTENT(in)    :: ftol       !<force residual tolerance
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec),PARAMETER :: fsq_prec = 1.E-10_dp !<threshold force for turning on block preconditioner
      LOGICAL               :: l_iterate            !<call evolve until fsq_ratio stops decreasing
      REAL(dp)              :: fsq_block, t1, t2
!-----------------------------------------------
!
!     Written by S.P.Hirshman 060308
!
      fsq_block = MAX(ftol, fsq_prec)
      l_iterate = .TRUE.
      l_Gmres   = .TRUE.
      nprecon_type = PREDIAG
      IF (.NOT.lprecon_diag) nprecon_type = PREBLOCK
      niter = 1
      nprecon = 0

      DO WHILE (l_iterate)
         in_hess_nfunct=0; out_hess_nfunct=0
#if defined(SKS)
         CALL second0(t1)
#endif         
         CALL evolve
#if defined(SKS)
         CALL second0(t2)
         diag_evolve_time=diag_evolve_time+(t2-t1)
#endif         
         IF (l_output_alliter) CALL write_output (wout_file, niter) 

         l_iterate = (fsq_ratio.LE.0.5_dp .AND. fsq_total.GT.fsq_block)

         IF (ladd_pert .AND. fsq_total.lt.100*fsq_block) THEN
            l_init_state=.TRUE.
            CALL second0(t1)
            CALL add_perturb(xc, getwmhd)
            CALL second0(t2)
            diag_add_pert_time=diag_add_pert_time+(t2-t1)
            ladd_pert = .FALSE.
         END IF
         niter = niter+1
      END DO

      END SUBROUTINE converge_diagonal


!>  \brief Performs convergence of force residuals using block preconditioner.
!!  Terminates when the force residual drops below ftol or the number of
!!  iterations exceeds a specified value (hard coded for now).
!!  Applies an external island perturbation if not previously done.
      SUBROUTINE converge_blocks (wout_file, ftol)
      USE PERTURBATION, ONLY: niter_max, l_output_alliter, write_restart_file, &     ! R. Sanchez, added Jul '10: Dump info on file at each timestep
                              ladd_pert, add_perturb
      USE dumping_mod, ONLY: write_output
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*), INTENT(in)  :: wout_file  !<name of VMEC wout file that SIESTA reads
      REAL(rprec), INTENT(in)    :: ftol       !<force residual tolerance
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nrow
      REAL(rprec),PARAMETER :: fsq_prec = 1.E-10_dp !<threshold force for turning on block preconditioner
      LOGICAL            :: l_iterate          !<controls iteration sequence
      REAL(dp) :: t1, t2
!-----------------------------------------------
!
!     Written by S.P.Hirshman 060308
!
      l_Gmres   = .FALSE.
      nprecon_type = PREBLOCK
      nprecon = 1
      niter_max = niter+niter_max
      l_iterate = (niter.lt.niter_max) .and. (fsq_total1 .gt. ftol) 
      nrow = 0

      DO WHILE (l_iterate)
         in_hess_nfunct=0; out_hess_nfunct=0
#if defined(SKS)
         CALL second0(t1)
#endif         
         CALL evolve
#if defined(SKS)
         CALL second0(t2)
         block_evolve_time=block_evolve_time+(t2-t1)
#endif         

         IF (l_output_alliter) CALL write_output (wout_file, niter) 

         l_iterate = (niter.LT.niter_max) .AND. (fsq_total1.GT.ftol) 

         IF (.NOT.l_conjmin .AND. .NOT.l_gmres) THEN
            IF (IAM .EQ. 0) PRINT *,'CONJ GRADIENT UNABLE TO MAKE PROGRESS'
            l_iterate=.FALSE.
         ELSE IF (nprecon.GT.6 .AND. fsq_ratio.GT.1.E5_dp) THEN
            IF (IAM .EQ. 0) PRINT 100, fsq_ratio
            l_iterate=.FALSE.
         END IF
!STOP IF LACK OF PROGRESS
         IF (nprecon.GT.6 .AND. ABS(1-fsq_ratio1).LT.1.E-2_dp) THEN
            nrow = nrow+1
            IF (nrow .EQ. 2) THEN
               l_iterate=.FALSE.
               IF (IAM .EQ. 0) PRINT 100, fsq_ratio1
            END IF
         ELSE 
            nrow = 0
         END IF
!In case we didn't add it in diag loop
         IF (ladd_pert .AND. fsq_total.LT.100*fsq_prec) THEN
           l_init_state=.TRUE.
           CALL second0(t1)
           CALL add_perturb(xc, getwmhd)
           CALL second0(t2)
           block_add_pert_time=block_add_pert_time+(t2-t1)
           ladd_pert = .FALSE.
         END IF
         nprecon = nprecon+1
         niter = niter+1
      END DO

100   FORMAT (' FSQ RATIO: ',1p,e10.2, ' SIESTA TERMINATING!')


      END SUBROUTINE converge_blocks

!>   \brief  Initializes variables and pointers prior to calling evolve
      SUBROUTINE init_evolution
      USE hessian, ONLY: levmarq_param0, levmarq_param, l_Compute_Hessian
      USE quantities, ONLY: fsubsmncf, fsubumnsf, fsubvmnsf,            &
                            jvsupsmncf,jvsupumnsf,jvsupvmnsf
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: istat, n1
!-----------------------------------------------
!     "natural" boundary condition (preserves hessian symmetry): push ONLY m=1 u-flux at origin
!     and correctly calculates B^s/r0_B^u for m=1 at origin
      IF (l_natural) THEN
         l_push_s = .FALSE.
         l_push_u = .FALSE.   !.TRUE.  
         l_push_v = .FALSE.   !.TRUE.
      ELSE
!     evolve all velocities (m=1 components) at s=0
!     results in J^s = r0 J^u (m=1,s=0) in equilibrium
         l_push_s = .FALSE.   !.TRUE.
         l_push_u = .FALSE.   !.TRUE.
         l_push_v = .FALSE.   !.TRUE.
      END IF

      ste = 0
      l_EPAR = .FALSE.

      nfun_calls = 0
      fsqprev = -1;  fsqprev1 = -1
      fsq_min = 1.E20_dp
      tau_damp = 0
      niter    = 0
      nprecon  = 0
      wtotal0  = -1
      etak_tol = 1.E-01_dp            !initial Newton tolerance parameter
      delta_t = 1
      l_linearize = .FALSE.
      l_getfsq = .TRUE.
      l_getwmhd= .FALSE.              !set in getwmhd function
      levmarq_param = levmarq_param0
      fsq_total = 1
      fsq_total1= 1

      IF (mnmax .ne. (1+mpol)*(2*ntor+1)) STOP 'WRONG mnmax VALUE!'
      n1 = ns*mnmax
      neqs = 3*n1

      ALLOCATE(xc(neqs), stat=istat)
      IF (istat .ne. 0) STOP 'Allocate xc failed!'
      xc = 0
      CALL init_ptrs (xc, jvsupsmncf, jvsupumnsf, jvsupvmnsf)

      ALLOCATE(gc(neqs), stat=istat)
      IF (istat .ne. 0) STOP 'Allocate gc failed!'
      CALL init_ptrs (gc, fsubsmncf, fsubumnsf, fsubvmnsf)

      l_Compute_Hessian = .FALSE.

      END SUBROUTINE init_evolution

!>    \brief Initializes pointers for xc (dependent variables) and gc (force) arrays

      SUBROUTINE init_ptrs (xtarget, ptr1, ptr2, ptr3)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), TARGET, INTENT(inout)     :: xtarget(0:mpol,-ntor:ntor,3,ns) !<input 3D, 3 component array
      REAL(rprec), POINTER, DIMENSION(:,:,:) :: ptr1  !<first (s) array subsection
      REAL(rprec), POINTER, DIMENSION(:,:,:) :: ptr2  !<second (u) array subsection
      REAL(rprec), POINTER, DIMENSION(:,:,:) :: ptr3  !<third (v) array subsection
!-----------------------------------------------
      xtarget = 0
      
      ptr1 => xtarget(:,:,1,:)
      ptr2 => xtarget(:,:,2,:)
      ptr3 => xtarget(:,:,3,:)

      END SUBROUTINE init_ptrs

!>   \brief Evolves one step toward MHD equilibrium
      SUBROUTINE evolve
      USE quantities, ONLY: jvsupsmncf, jvsupumnsf, jvsupvmnsf,         &
                            jbsupsmnsh, jbsupumnch, w_factor, signjac
      USE gmres, ONLY: gmres_fun
      USE hessian
      USE perturbation, ONLY: eta_factor, lresistive, ftol,             &
                              write_restart_file, niter_max, lPosDef
      USE diagnostics_mod, ONLY: bgradp_rms, toroidal_flux0, bgradp,    &
                                 bgradp_par, tflux
      USE siesta_state, ONLY: update_state, update_state_par
#if defined(SKS)
      USE timer_mod, ONLY: gmres_time
#endif
      IMPLICIT NONE
!
!     Evolves the equilibrium equations one time step and updates delta_t
!     
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: ns0=3, nm0=1, nn0=2
      REAL(dp)    :: skston, skstoff
      INTEGER     :: iprint, m1, i
      INTEGER     :: imaxloc(3), ilbound(3)
      REAL(rprec) :: v2tot, f1, r0, fmax, decay_exp
      REAL(rprec) :: vmaxs, vmaxu, vmaxv, t1, t2
      LOGICAL     :: lprint, ldiagonal, lblock
      CHARACTER(18) :: str_resistive
      INTEGER, DIMENSION (2) :: snum_funct, rnum_funct 
      REAL(dp), PARAMETER  :: levm_ped=1.E-10_dp, mu_ped=1.E-6_dp
!-----------------------------------------------
#if defined(SKS)
      EVOLVEPASS=EVOLVEPASS+1
#endif
!
!     Update COVARIANT Fourier components of velocity field
!
!     THIS IS THE TIME-STEP ALGORITHM. IT IS ESSENTIALLY A CONJUGATE
!     GRADIENT METHOD, WITHOUT THE LINE SEARCHES (FLETCHER-REEVES),
!     BASED ON A METHOD GIVEN BY P. GARABEDIAN
!
      l_init_state = .TRUE.
!
!     Compute preconditioner
!
      IF (nprecon .GT. 1) THEN
!        PSEUDO-TRANSIENT CONTINUATION: TURN OFF FOR CONSTANT levmarq_param, muPar option
         f1 = fsq_total*1.E10_dp
         IF (ngmres_type .eq. 1) THEN
            decay_exp = 0.5_dp
         ELSE IF (ngmres_type .eq. 2) THEN
            decay_exp = 0.75_dp  
         END IF
         f1 = f1**decay_exp
         levmarq_param = (levmarq_param0-levm_ped)*MIN(one, f1) + levm_ped
         muPar = (muPar0-mu_ped)*MIN(one,f1) + mu_ped
     
!         IF (levmarq_param0 .EQ. zero) levmarq_param=0            !!Improves convergence a bit, levm_ped=10^(-8)
         IF (muPar0 .EQ. zero) muPar=0

         fmax = fsq_total1*1.E12_dp               !SPH: raise this if convergence problem occurs
         xPosDef = MAX(1.E-3_dp, fmax/(1+fmax))

      ELSE
         xPosDef = 1
         levmarq_param = levmarq_param0
         muPar = muPar0
      END IF

      IF (.NOT.lposDef) xPosDef=0

!
!     EXCEPT FOR THE INITIAL STEP (WHERE ONLY INITIAL UNPRECONDITIONED FORCES
!     ARE NEEDED) COMPUTE FULL (OR DIAGONAL APPROXIMATION) HESSIAN
!
      IF (niter .GT. 1) THEN
         l_linearize = .TRUE.
         l_getfsq = .FALSE.
         l_ApplyPrecon = .FALSE.
         l_getwmhd = .FALSE.
         ldiagonal=(nprecon_type.eq.PREDIAG .AND. niter.eq.2)
         lblock   =(nprecon_type .EQ. PREBLOCK)

         IF (ldiagonal .OR. lblock) THEN

            CALL second0(skston)
            IF (PARFUNCTISL) THEN
               CALL Compute_Hessian_Blocks (xc, gc, funct_island_par, l_Hess_sym, ldiagonal)
            ELSE
               CALL Compute_Hessian_Blocks (xc, gc, funct_island, l_Hess_sym, ldiagonal)
            END IF
            CALL second0(skstoff)
            IF (ldiagonal) comp_diag_elements_time=comp_diag_elements_time+(skstoff-skston)
            IF (lblock)    compute_hessian_time=compute_hessian_time+(skstoff-skston)

         END IF
      END IF

!
!     Reset run-control logicals
!     (may have been reset to false in Compute_Hessian=>funct_island call)

      l_init_state = .TRUE.   
      l_ApplyPrecon = .FALSE.
      l_linearize = .FALSE.
      l_getfsq = .TRUE.
      fsq_lin = -1

      l_Gmres = (fsq_total.LE.fsq_res .OR. nprecon_type.EQ.PREDIAG)

!THIS USUALLY CONVERGES BETTER AND REQUIRES FEW FUNCTION EVALUATIONS
      l_Gmres = .TRUE.

!
!     Choose time-stepping algorithm
!
      IF (niter .EQ. 1) THEN
         l_init_state=.TRUE.
         IF (PARFUNCTISL) THEN
            l_getwmhd=.TRUE.
            CALL second0(skston)
            CALL funct_island_par
            CALL second0(skstoff)
            evolve_funct_island_time=evolve_funct_island_time+(skstoff-skston)
            l_getwmhd=.FALSE.
         ELSE
            CALL second0(skston)
            CALL funct_island
            CALL second0(skstoff)
            evolve_funct_island_time=evolve_funct_island_time+(skstoff-skston)
         END IF
      ELSE IF (nprecon .EQ. 1) THEN
         IF (l_Gmres) THEN
           CALL second0(skston)
           CALL gmres_fun
           CALL second0(skstoff)
           gmres_time=gmres_time+(skstoff-skston)
         END IF
      ELSE 
         etak_tol = MIN(0.1_dp, 1.E8_dp*fsq_total)
         IF (fsq_total <= 0) etak_tol = 0.1_dp
         IF (l_Gmres) THEN
           CALL second0(skston)
           CALL gmres_fun
           CALL second0(skstoff)
           gmres_time=gmres_time+(skstoff-skston)
         END IF
      END IF

!SPH042808  TRY TO IMPROVE FSQ_TOTAL: BEFORE RESISTIVE E-field added!
      lprint = (MOD(niter, nprint_step).EQ.0 .OR. niter.EQ.1)
      lblock = (nprecon_type.EQ.PREBLOCK .AND. niter.GT.1)
      IF (lblock) THEN 
         CALL second0(skston)
         CALL Conjugate_Grad(ftol)
         CALL second0(skstoff)
         conj_grad_time=conj_grad_time+(skstoff-skston)
      END IF

      v2tot = SQRT(hs*SUM(xc*xc))
!      IF (v2tot.lt.1.E-5_dp .and. .not.l_Gmres .and. niter.gt.1) THEN
!      END IF
        
      lprint = MOD(niter, nprint_step).eq.0 .or. niter.eq.1
      l_update_state=.TRUE.
      IF (PARFUNCTISL) THEN
         IF (.NOT.l_par_state) CALL init_state_par(.FALSE.)
         CALL update_state_par(lprint, fsq_total1, zero)
      ELSE
         IF (l_par_state) CALL init_state(.FALSE.)
         CALL update_state(lprint, fsq_total1, zero)
      END IF
      l_update_state=.FALSE.

!      IF (v2tot.lt.1.E-5_dp .and. .not.l_Gmres .and. niter.gt.1) THEN
!      END IF
        
!
!     Compute force component residues  |Fsubi|**2
!
      IF (fsq_total < zero) PRINT *, 'Fsq_Total = ', fsq_total,' < 0!'
      IF (fsqprev .GT. zero) THEN
         fsq_ratio=fsq_total/fsqprev;  fsq_ratio1=fsq_total1/fsqprev1
      ELSE
         fsq_ratio=0;  fsq_ratio1=0;
      END IF
      fsqprev=fsq_total; fsqprev1=fsq_total1

!     IF v-step is too large, recompute Hessian
!      IF (v2tot*delta_t .ge. one) fsq_ratio = 1

!
!     Convert output to REAL (VMEC) units, divide B-fields by b_factor
!     p by p_factor, WMHD by w_factor

      out_hess_nfunct = out_hess_nfunct/nprocs
#if defined(SKS)
      snum_funct(1)=in_hess_nfunct; snum_funct(2)=out_hess_nfunct
      CALL MPI_Allreduce(snum_funct,rnum_funct,2,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,MPI_ERR)
      nfun_calls=nfun_calls+(rnum_funct(1)+rnum_funct(2))/ns
#else 
      nfun_calls=nfun_calls+(in_hess_nfunct+out_hess_nfunct)/ns
#endif

     IF (niter .EQ. 1) THEN
        IF (PARFUNCTISL) THEN
           CALL init_state_par(.TRUE.)
           CALL bgradp_par
        ELSE
           CALL bgradp
        END IF
        CALL tflux
     END IF

     IF (iam.EQ.0 .AND. lprint) THEN

        IF (niter .GT. 1) THEN
          ilbound = LBOUND(jvsupsmncf)-1
          imaxloc = MAXLOC(jvsupsmncf**2) + ilbound
          vmaxs = jvsupsmncf(imaxloc(1),imaxloc(2),imaxloc(3))
          imaxloc = imaxloc - ilbound
          IF (vmaxs .NE. zero)                                          &
          PRINT 45, ' VS-max: ',vmaxs,imaxloc(ns0),imaxloc(nm0)-1,imaxloc(nn0)-ntor-1
          imaxloc = MAXLOC(jvsupumnsf**2) + ilbound
          vmaxu = jvsupumnsf(imaxloc(1),imaxloc(2),imaxloc(3))
          imaxloc = imaxloc - ilbound
          IF(vmaxu .ne. zero)                                           &
          PRINT 45, ' VU-max: ',vmaxu,imaxloc(ns0),imaxloc(nm0)-1,imaxloc(nn0)-ntor-1
          imaxloc = MAXLOC(jvsupvmnsf**2) + ilbound
          vmaxv = jvsupvmnsf(imaxloc(1),imaxloc(2),imaxloc(3))
          imaxloc = imaxloc - ilbound
          IF (vmaxv .ne. zero)                                          &
          PRINT 45, ' VV-max: ',vmaxv,imaxloc(ns0),imaxloc(nm0)-1,imaxloc(nn0)-ntor-1
          PRINT 110, bsbu_ratio, jsju_ratio,                            &
                    (bs0(i), i=2,6), (bu0(i), i=2,6)
        END IF

 45   FORMAT(a,1pe10.2,' AT JS: ',i4,' M: ',i4,' N: ',i4)

        DO iprint = 6, 33, 27
           IF (niter .EQ. 1) THEN
              WRITE (iprint, *) 'INITIAL PHYSICS PARAMETERS'
              WRITE (iprint, 50) wtotal/w_factor, wp_i/wb_i,            &
                                 toroidal_flux0, bgradp_rms
              IF (lresistive) THEN
                 str_resistive = "RESISTIVE RUN "
              ELSE
                 str_resistive = "NON-RESISTIVE RUN "
              END IF
              WRITE (iprint, 55) str_resistive, eta_factor, l_push_s,   &
                                 l_push_u, l_push_v, l_push_edge
              WRITE (iprint, 60) delta_t, delt_cg, etak_tol, l_Hess_sym
              WRITE (iprint, 65) levmarq_param0, mupar0, lposdef
              WRITE (iprint, 70) neqs, ns, mpol, ntor, nu_i, nv_i,      &
                                 ngmres_steps
           END IF
            
           IF (fsq_lin .eq. -1) THEN
              WRITE (iprint, 100)                                       &
                              niter, (wtotal-wtotal0)*1.E6_dp/wtotal0,  &
                              fsq_total1, fsqvs, fsqvu,                 &
                              fsqvv, v2tot*delta_t, nfun_calls
           ELSE
              WRITE (iprint, 102)                                       &
                              niter, (wtotal-wtotal0)*1.E6_dp/wtotal0,  &
                              fsq_total1, fsq_lin, fsqvs, fsqvu,        &
                              fsqvv, v2tot*delta_t, nfun_calls
           END IF
        END DO
     END IF

 50  FORMAT(' WMHD: ', 1pe13.6,' <BETA>: ',1pe11.4,                     &
            ' TFLUX: ',1pe11.4,' B.GRAD-P (rms): ', 1pe11.4,/,21('-'))
 55  FORMAT(1x, a,' ETA_FACTOR: ', 1pe10.2,' L_PUSH_S: ',l1,            &
            ' L_PUSH_U: ',l1,' L_PUSH_V: ',l1,' L_PUSH_EDGE: ',l1)
 60  FORMAT(' DELTA_T: ',1pe10.2,' DELT_CG: ',1pe10.2, ' ETA_K: ',      &
             1pe10.2,' HESSIAN SYM: ', l1)
 65  FORMAT(' LEVMARQ_PARAM: ',1pe10.2,' MU_PAR: ',1pe10.2,             &
            ' LPOSDEF: ',L2)
 70  FORMAT(' NEQS: ', i6,' NS: ',i4,' MPOL: ',i4,' NTOR: ',i4,         &
             ' NTHETA: ', i4,' NZETA: ',i4, ' NGMRES-STEPS: ', i4,//,   &
             ' NITER (W-W0)/W0*1E6    F2(MHD)    F2(LIN)',              &
             '    F2SUBS     F2SUBU     F2SUBV     |V|rms    NCALLS')
 100 FORMAT(1x,i5,1x, 1pe12.4, 2x, 1x,1pe10.3, 2x,9 ('-'),              &
           4(1x,1pe10.3),3x,i5)
 102 FORMAT(1x,i5,1x, 1pe12.4, 2x, 6(1x,1pe10.3),3x,i5)
 110 FORMAT(' |B^s-r0*B^u|/|B^s+r0*B^u| (m=1,r->0)  : ',1pe10.3,/,      &
             ' |J^s-r0*J^u|/|J^s+r0*J^u| (m=1,r->0)  : ',1pe10.3,/,     &
             ' JBSUPSH(JS=2-6,M=1)/R0 : ',1p5e10.3,/,                   &
             ' JBSUPUH(JS=2-6,M=1)    : ',1p5e10.3,/)

     CALL second0(skston)
     IF (fsq_total1 .LT. fsq_min) THEN
        fsq_min = fsq_total1
        CALL write_restart_file
     END IF
     CALL second0(skstoff)
     evolve_restart_file_time=evolve_restart_file_time+(skstoff-skston)

!    Add resistive part of B perturbation
     CALL second0(skston)
     IF (lresistive .AND. nprecon.GE.1 .AND. fsq_total1.GT.fsq_res) THEN
        l_getfsq = .FALSE.
        CALL add_resistive_E
     END IF
     CALL second0(skstoff)
     evolve_add_resistive_E_time=evolve_add_resistive_E_time+(skstoff-skston)

     END SUBROUTINE evolve


!>  \brief Parallel routine to evolve state toward equilibrium, using Conjugate Gradients

      SUBROUTINE Conjugate_Grad (ftol)
      USE stel_kinds
      USE hessian, ONLY: mblk_size, ns
#if defined(SKS)      
      USE nscalingtools, ONLY: PARFUNCTISL, MPI_ERR, startglobrow,         &
                            endglobrow, rcounts, disp
#endif
      IMPLICIT NONE
#if defined(MPI_OPT)
      INCLUDE 'mpif.h'
#endif
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec) :: ftol
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER     :: icount, nloc, myrowstart, myrowend
      REAL(rprec) :: b1, fsq_start, fsq_min0, fsq_save,                 &
                     wmhd_min, wmhd_start
      REAL(rprec), ALLOCATABLE :: xcsave(:), gc1(:)
      LOGICAL     :: lpar=.FALSE.
!-----------------------------------------------
      l_getfsq = .TRUE.
      l_linearize = .FALSE.
      l_ApplyPrecon = .TRUE.
      l_init_state = .TRUE.

#if defined(SKS)      
      lpar=PARFUNCTISL
      IF (lpar) THEN
         l_getwmhd=.TRUE.
         nloc=(endglobrow-startglobrow+1)*mblk_size
         myrowstart=(startglobrow-1)*mblk_size+1
         myrowend=myrowstart+nloc-1
         ALLOCATE(gc1(neqs))
      END IF
#endif

!
!     Linearize perturbations around B0 and p0
!     dB = curl(xc X B0)
!     dp = dp(p0,xc)
!
!     NON-LINEAR forces 
!     J = J0+dJ,  B = B0+dB, p = p0+dp
!
      IF (nprecon .EQ. 0) THEN
         CALL funct_island
         RETURN
      END IF

!TESTING      lpar=.FALSE.

      ALLOCATE(xcdot(neqs), xc0(neqs), xcsave(neqs), stat=icount) 
      IF (icount .NE. 0) STOP 'ALLOCATION FAILED IN CONJ_GRAD'

      xcsave = xc
      xcdot = 0
      xc = 0
      xc0 = 0

      delt_cg = 0.5_dp
      otau = 0.85_dp*delt_cg

!      IF (IAM.EQ.0) PRINT *,'IN CONJ-GRAD'

      DO icount = 1, ncongrad_steps
!
!     LOOP OVER delt_cg TO FIND MINIMUM F(xc)
!
#if defined(SKS)
         IF (lpar) THEN
            CALL funct_island_par
            gc1 = gc
            CALL MPI_Allgatherv(gc1(myrowstart:myrowend),nloc,MPI_REAL8,   &
                     gc,rcounts,disp,MPI_REAL8,MPI_COMM_WORLD,MPI_ERR)
         ELSE
#endif
            CALL funct_island
#if defined(SKS)
         END IF
#endif
         CALL update_taudamp (fsq_save)

         IF (icount .EQ. 1) THEN
            fsq_start  = fsq_total1
            fsq_min0   = fsq_total1
            fsq_save   = fsq_total
            wmhd_start = wtotal
            wmhd_min   = wtotal
         END IF

         IF (fsq_total1 .LT. fsq_min0) THEN
            fsq_min0 = fsq_total1
            xc0 = xc
         END IF
!         IF (wtotal .lt. wmhd_min) THEN
!            wmhd_min = wtotal
!            xc0 = xc
!         END IF
!         PRINT ('(i4,1p,e12.3)'), icount, fsq_total1

         IF (fsq_total1 .GT. 1.5_dp*fsq_min0) THEN
            delt_cg = 0.95_dp*delt_cg
            IF (fsq_total1 .GT. 10*fsq_min0) delt_cg = delt_cg/2
            IF (fsq_total1.GT.100*fsq_min0 .or. delt_cg.LT.0.02_dp) EXIT
            xc = xc0
            xcdot = 0
         END IF

         IF (l_linearize) STOP 'LINEARIZATION TURNED ON!'

         b1 = 1-tau_damp
!
!     THIS IS THE TIME-STEPPING ALGORITHM. IT IS ESSENTIALLY A CONJUGATE
!     GRADIENT METHOD, WITHOUT THE LINE SEARCHES (FLETCHER-REEVES),
!     BASED ON A METHOD BY P. GARABEDIAN
!
!        CONJ GRAD 
         xcdot = xcdot*b1 - delt_cg*gc
         xc    = xc + delt_cg*xcdot
!        STEEPEST DESCENT (b1->0)
!         xc = xc*b1 - delt_cg*gc

      END DO

!      IF (fsq_min0.lt.fsq_start .or. wmhd_min.lt.wmhd_start) THEN
      IF (.NOT.l_Gmres) l_PrintOriginForces = .TRUE.
      l_conjmin = (fsq_min0 .lt. fsq_start)

      IF (l_conjmin) THEN
         xc=xc0
      ELSE
         xc=0
      END IF

      l_init_state=.TRUE.
      l_update_state=.TRUE.
      IF (lpar) THEN
         CALL funct_island_par
         DEALLOCATE(gc1)
      ELSE
         CALL funct_island
      END IF
      l_update_state=.FALSE.

      IF (l_conjmin) THEN
         IF (iam .EQ. 0) THEN
            IF (fsq_min0.NE.fsq_total1)                                &
               PRINT *, 'ERROR in CONJ_GRADIENT CALC OF FSQ_MIN0!'
            WRITE (6, 100) icount, fsq_start, fsq_min0
            WRITE (33,100) icount, fsq_start, fsq_min0
         END IF
      ELSE
         xc = xcsave
      END IF

      l_PrintOriginForces = .FALSE.

      DEALLOCATE (xcdot, xc0, xcsave)
 100  FORMAT(' ITER_CG: ', i3,' FSQ (START CG): ',1p,e12.3,' FSQ (END CG): ', 1pe12.3)

      END SUBROUTINE Conjugate_Grad

!>  \brief Updates damping parameter (tau) used by conjugate gradient routine

      SUBROUTINE update_taudamp (fsq_save)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), INTENT(INOUT) :: fsq_save
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: mintau = 0.15_dp
      REAL(rprec), PARAMETER :: b = 0.998_dp
      REAL(rprec)            :: ratio
!-----------------------------------------------

      IF (fsq_save .GT. zero) THEN
         ratio = ABS(fsq_total/fsq_save)
         ratio = MIN(b, ratio)
         fsq_save = fsq_total
      ELSE
         ratio = b
      END IF

      tau_damp = 1 - ratio

      END SUBROUTINE update_taudamp

      END MODULE evolution

