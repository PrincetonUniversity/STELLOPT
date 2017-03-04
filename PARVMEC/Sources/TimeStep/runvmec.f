      SUBROUTINE runvmec(ictrl_array, input_file0, 
     1                   lscreen, COMM_WORLD, reset_file_name)
      USE vmec_main
      USE vmec_params, ONLY: bad_jacobian_flag, more_iter_flag,
     1                       norm_term_flag, successful_term_flag,
     2                       restart_flag, readin_flag, output_flag,
     3                       timestep_flag, ns_error_flag, 
     4                       reset_jacdt_flag, lamscale
#if defined (SKS)      
      USE realspace
      USE vmec_params, ONLY: ntmax
      USE vacmod, ONLY: nuv, nuv3
#endif      
      USE vsvd
      USE timer_sub
      USE parallel_include_module
      USE parallel_vmec_module, ONLY: MyEnvVariables
      USE parallel_vmec_module, ONLY: InitRunVmec
      USE parallel_vmec_module, ONLY: FinalizeRunVmec
      USE parallel_vmec_module, ONLY: InitSurfaceComm
      USE parallel_vmec_module, ONLY: FinalizeSurfaceComm
      USE parallel_vmec_module, ONLY: SetVacuumCommunicator
      USE parallel_vmec_module, ONLY: Scatterv4XArray
      USE blocktridiagonalsolver_bst, ONLY: Initialize_bst, Finalize_bst
      USE xstuff
      USE mpi_inc
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(inout), TARGET :: ictrl_array(5)
      LOGICAL, INTENT(in) :: lscreen
      CHARACTER(LEN=*), INTENT(in) :: input_file0
      CHARACTER(LEN=*), OPTIONAL :: reset_file_name
      !INTEGER, INTENT(IN) :: COMM_WORLD
      INTEGER :: COMM_WORLD
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, POINTER :: ier_flag
      INTEGER :: ictrl_flag, iseq_count
      INTEGER :: ns_index, ns_min, nsval, ns_old=0, numsteps
      INTEGER :: igrid, index_end, index_dat,
     1           jacob_off, niter_store
      INTEGER, SAVE :: igrid0
      CHARACTER(LEN=120) :: input_file
      LOGICAL :: lreset 
#if defined(SKS)
      REAL(dp) :: rvton, rvtoff
      REAL(dp) :: gridton, gridtoff
      REAL(dp) :: skston, skstoff
      REAL(dp) :: allgvton, allgvtoff
      REAL(dp) :: bcastton, bcasttoff
      INTEGER :: max_grid_size, flag
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: bcastarr
      INTEGER :: blklength, grid_id, i, js, nsmin, nsmax
      CHARACTER(LEN=20) :: fname
#endif
C-----------------------------------------------
!
!     ictrl_flag = ictrl_array(1)
!                  flag that controls calling of various subroutines of vmec code
!                  add together the values beow to utilize several subroutines with one call
!
!            value     flag-name              calls routines to...
!            -----     ---------              ---------------------
!              1       restart_flag           reset internal run-control parameters (for example, if 
!                                             jacobian was bad, to try a smaller time-step)
!              2       readin_flag            read in data from input_file and initialize parameters/arrays
!                                             which do not dependent on radial grid size
!                                             allocate internal grid-dependent arrays used by vmec;
!                                             initialize internal grid-dependent vmec profiles (xc, iota, etc); 
!                                             setup loop for radial multi-grid meshes or, if ns_index = ictrl_array(4)
!                                             is > 0, use radial grid points specified by ns_array[ns_index]
!              4       timestep_flag          iterate vmec either by "niter" time steps or until ftol satisfied,
!                                             whichever comes first. If numsteps (see below) > 0, vmec will return
!                                             to caller after numsteps, rather than niter, steps.
!              8       output_flag            write out output files (wout, jxbout)
!             16       cleanup_flag           cleanup (deallocate arrays) - this terminates present run of the sequence
!                                             This flag will be ignored if the run might be continued. For example, 
!                                             if ier_flag (see below) returns the value more_iter_flag, the cleanup
!                                             code will be skipped even if cleanup_flag is set, so that the run
!                                             could be continued on the next call to runvmec.
!             32       reset_jacdt_flag       Resets ijacobian flag and time step to delt0
!
!                  thus, setting ictrl_flag = 1+2+4+8+16 will perform ALL the tasks thru cleanup_flag
!                  in addition, if ns_index = 0 and numsteps = 0 (see below), vmec will control its own run history
!
!     ier_flag = ictrl_array(2)
!                  specifies vmec error condition (if nonzero)
!     numsteps = ictrl_array(3)
!                  number time steps to evolve the equilibrium. Iterations will stop EITHER if numsteps > 0 and
!                  when the number of vmec iterations exceeds numsteps; OR if the ftol condition is satisfied, 
!                  whichever comes first. The timestep_flag must be set (in ictrl_flag) for this to be in effect.
!                  If numsteps <= 0, then vmec will choose consecutive (and increasing) values from the ns_array 
!                  until ftol is satisfied on each successive multi-grid.
!     ns_index = ictrl_array(4)
!                  if > 0 on entry, specifies index (in ns_array) of the radial grid to be used for the present iteration
!                  phase. If ns_index <= 0, vmec will use the previous value of this index (if the ftol 
!                  condition was not satisfied during the last call to runvmec) or the next value of this index,
!                  and it will iterate through each successive non-zero member of the ns_array until ftol-convergence
!                  occurs on each multigrid.
!                  on exit, contains last value of ns_array index used
!     iseq_count=ictrl_array(5)
!                  specifies a unique sequence label for identifying output files in a sequential vmec run
C-----------------------------------------------
      INTERFACE
         SUBROUTINE initialize_radial(nsval, ns_old, delt0,
     1                               lscreen, reset_file_name)
         USE vmec_main
         IMPLICIT NONE
         INTEGER, INTENT(in) :: nsval
         INTEGER, INTENT(inout) :: ns_old
         CHARACTER(LEN=*), OPTIONAL :: reset_file_name
         LOGICAL, INTENT(in) :: lscreen
         REAL(rprec), INTENT(out) :: delt0
         END SUBROUTINE initialize_radial
      END INTERFACE

#if defined(SKS)
      RUNVMEC_PASS = RUNVMEC_PASS + 1
      CALL second0(rvton)
      CALL MyEnvVariables
      IF(LV3FITCALL) THEN
        CALL Serial2Parallel4X(xc,pxc)
        CALL Serial2Parallel4X(xcdot,pxcdot)
        CALL Serial2Parallel4X(xstore,pxstore)
        CALL InitRunVmec(COMM_WORLD,lfreeb)
        IF (RUNVMEC_PASS.GT.1) THEN
          CALL second0(bcastton)
          CALL MPI_Bcast(pxc,SIZE(pxc),MPI_REAL8,
     1                      0,RUNVMEC_COMM_WORLD,MPI_ERR)
          CALL MPI_Bcast(pxcdot,SIZE(pxcdot),MPI_REAL8,
     1                      0,RUNVMEC_COMM_WORLD,MPI_ERR)
          CALL MPI_Bcast(pxstore,SIZE(pxstore),MPI_REAL8,
     1                      0,RUNVMEC_COMM_WORLD,MPI_ERR)
          CALL MPI_Bcast(iotas,SIZE(iotas),MPI_REAL8,
     1                      0,RUNVMEC_COMM_WORLD,MPI_ERR)
          CALL MPI_Bcast(iotaf,SIZE(iotaf),MPI_REAL8,
     1                      0,RUNVMEC_COMM_WORLD,MPI_ERR)
          CALL MPI_Bcast(phips,SIZE(phips),MPI_REAL8,
     1                      0,RUNVMEC_COMM_WORLD,MPI_ERR)
          CALL MPI_Bcast(phipf,SIZE(phipf),MPI_REAL8,
     1                      0,RUNVMEC_COMM_WORLD,MPI_ERR)
          CALL MPI_Bcast(chips,SIZE(chips),MPI_REAL8,
     1                      0,RUNVMEC_COMM_WORLD,MPI_ERR)
          CALL MPI_Bcast(chipf,SIZE(chipf),MPI_REAL8,
     1                      0,RUNVMEC_COMM_WORLD,MPI_ERR)
          CALL MPI_Bcast(mass,SIZE(mass),MPI_REAL8,
     1                      0,RUNVMEC_COMM_WORLD,MPI_ERR)
          CALL MPI_Bcast(icurv,SIZE(icurv),MPI_REAL8,
     1                      0,RUNVMEC_COMM_WORLD,MPI_ERR)
          CALL MPI_Bcast(lamscale,1,MPI_REAL8,
     1                      0,RUNVMEC_COMM_WORLD,MPI_ERR)
          CALL second0(bcasttoff)
          broadcast_time = broadcast_time + (bcasttoff - bcastton)
        END IF
      END IF
#endif
      ictrl_flag = ictrl_array(1)
      numsteps = ictrl_array(3)
      ier_flag  => ictrl_array(2)
      ns_index   = ictrl_array(4)
      iseq_count = ictrl_array(5)
      CALL second0 (timeon)

!
!     PARSE input_file into path/input.ext
!
      index_dat = INDEX(input_file0,'input.')
      index_end = LEN_TRIM(input_file0)
      IF (index_dat .gt. 0) THEN
         input_file = TRIM(input_file0)
         input_extension  = input_file0(index_dat+6:index_end)
      ELSE
         input_extension = input_file0(1:index_end)
         input_file = 'input.'//TRIM(input_extension)
      END IF

!
!     INITIALIZE PARAMETERS
!
      lreset = (IAND(ictrl_flag, restart_flag) .ne. 0)

      IF (lreset) THEN
#if defined(SKS)
        CALL second0(skston)
#endif
        CALL reset_params
#if defined(SKS)
        CALL second0(skstoff)
        reset_params_time = reset_params_time + (skstoff-skston)
#endif
!       res0 = -1   Done in reset_params
      END IF

      IF (IAND(ictrl_flag, reset_jacdt_flag) .ne. 0) THEN
         ijacob = 0
         delt0r = delt
      END IF

      IF (IAND(ictrl_flag, readin_flag) .ne. 0) THEN
!
!        READ INPUT FILE (INDATA NAMELIST), MGRID_FILE (VACUUM FIELD DATA)
!
#if defined (SKS)
        CALL second0(skston)
#endif
        CALL vsetup (iseq_count)
#if defined(SKS)
        CALL second0(skstoff)
        vsetup_time = vsetup_time + (skstoff-skston)
        CALL second0(skston)
#endif
        CALL readin (input_file, iseq_count, ier_flag, lscreen)
#if defined(SKS)
        CALL second0(skstoff)
        readin_time = readin_time + (skstoff-skston)

        max_grid_size = ns_array(multi_ns_grid)
!        IF(lfreeb) CALL SetVacuumCommunicator(nuv, nuv3,
!     1              max_grid_size)
#endif
        IF (ier_flag .ne. 0) GOTO 1000
!
!        COMPUTE NS-INVARIANT ARRAYS
!
#if defined(SKS)
        CALL second0(skston)
#endif
        CALL fixaray
#if defined(SKS)
        CALL second0(skstoff)
        fixarray_time = fixarray_time + (skstoff-skston)
#endif
      END IF

!      SAL (This should only be called if we actually doing a run)
!      IF(lfreeb) CALL SetVacuumCommunicator(nuv, nuv3,
!     1              max_grid_size)

      IF (lreset) THEN
!
!        COMPUTE INITIAL SOLUTION ON COARSE GRID
!        IF PREVIOUS SEQUENCE DID NOT CONVERGE WELL
!
!        IF (lreseta) THEN    !NOTE: where externally, lreseta = T, set restart_flag bit 
!                                    (ictrl_flag = IOR(ictrl_flag,restart_flag))
         igrid0 = 1
         ns_old = 0
         IF (PRESENT(reset_file_name)) THEN
            IF (LEN_TRIM(reset_file_name) .ne. 0)igrid0 = multi_ns_grid
         END IF
         IF (grank.EQ.0) WRITE (nthreed, 30)
         delt0r = delt
      ENDIF

  30  FORMAT(' FSQR, FSQZ = Normalized Physical Force Residuals',/,
     1   ' fsqr, fsqz = Preconditioned Force Residuals',/,1x,23('-'),/,
     2   ' BEGIN FORCE ITERATIONS',/,1x,23('-'),/)

      IF (ALL(ns_array.eq.0) .and. ns_index.le.0) THEN
         ier_flag = ns_error_flag
         GOTO 1000
      END IF

      imovephi = 0
      jacob_off = 0

      IF (IAND(ictrl_flag, timestep_flag) == 0) GOTO 1000


!      SAL
      IF(lfreeb) CALL SetVacuumCommunicator(nuv, nuv3,
     1              max_grid_size)

  50  CONTINUE
      iequi = 0
      IF (lfreeb .and. jacob_off.eq.1) ivac = 1    !!restart vacuum calculations

      ns_min = 3
#if defined(SKS)

      num_grids = multi_ns_grid
      IF(.NOT.ALLOCATED(grid_procs)) THEN
        ALLOCATE(grid_procs(num_grids))
        ALLOCATE(grid_size(num_grids))
        ALLOCATE(grid_time(num_grids))
        ALLOCATE(f3d_time(num_grids))
        ALLOCATE(f3d_num(num_grids))
        IF (lfreeb) ALLOCATE(vgrid_time(num_grids))
      END IF

      f3d_time=0; f3d_num=0
      blklength=(ntor+1)*(mpol1+1)
      !BEGIN - Main loop that will be parallelized - SKS
      grid_id=1
      old_vacuum_time=0
#endif
      DO igrid = igrid0-jacob_off, multi_ns_grid
#if defined(SKS)
         CALL second0(gridton)
#endif
         IF (igrid .lt. igrid0) THEN
!           TRY TO GET NON-SINGULAR JACOBIAN ON A 3 PT RADIAL MESH
            nsval = 3; ivac = -1
            ftolv = 1.e-4_dp
         ELSE IF (ns_index .gt. 0) THEN
            IF (ns_index .gt. SIZE(ns_array)) THEN
               ier_flag = ns_error_flag
               RETURN
            END IF
            nsval = ns_array(ns_index)
            IF (nsval .le. 0) STOP 'NSVAL <= 0: WRONG INDEX VALUE'
            ftolv = ftol_array(ns_index)
            niter = niter_array(ns_index)
         ELSE
            nsval = ns_array(igrid)
            IF (nsval .lt. ns_min) CYCLE
            ns_min = nsval
            ictrl_array(4) = igrid
            ftolv = ftol_array(igrid)
            niter = niter_array(igrid)
         END IF


#if defined (SKS)
         IF(PARVMEC .AND. NS_RESLTN.GE.1) THEN
           IF (lactive) THEN
             CALL Gather4XArray(pscalxc)
             CALL Gather4XArray(pxc)
           END IF
           CALL FinalizeSurfaceComm(NS_COMM)
           
           CALL second0(skston)

           bcastton = skston
           CALL MPI_Bcast(pscalxc,SIZE(pscalxc),MPI_REAL8,
     1                      0,RUNVMEC_COMM_WORLD,MPI_ERR)
           CALL MPI_Bcast(pxc,SIZE(pxc),MPI_REAL8,
     1                      0,RUNVMEC_COMM_WORLD,MPI_ERR)
           CALL second0(bcasttoff)
           broadcast_time = broadcast_time + (bcasttoff - bcastton)

           IF (.FALSE..AND.LV3FITCALL) THEN
             nsmin=t1lglob; nsmax=t1rglob
             DO js = nsmin, nsmax
               pphip(:,js) = phips(js)
               pchip(:,js) = chips(js)
             END DO

             nsmin=t1lglob; nsmax=MIN(t2rglob,ns+1)
             pres(nsmin:nsmax) = 0
           END IF

           CALL second0(skstoff)
           broadcast_time = broadcast_time + (skstoff -skston)
         END IF

         CALL InitSurfaceComm(nsval, nzeta, ntheta3,
     1        ntmax, ntor, mpol1)
         CALL second0(skstoff)
         init_parallel_time = init_parallel_time + (skstoff-skston)

         grid_size(grid_id)=nsval
         grid_procs(grid_id)=nranks
#endif
         
!  JDH 2012-06-20. V3FIT fix, inserted with change from VMEC 8.48 -> 8.49
!  (Not sure just what in initialize_radial messes up convergence - happens slowly)
!  Logical l_v3fit is declared in vmec_input, available via vmec_main
#if defined(SKS)
         CALL second0(skston)
#endif
         IF (l_v3fit) THEN                            ! V3FIT is running here 
            IF (ns_old .ne. nsval) THEN
               CALL initialize_radial(nsval, ns_old, delt0r,
     1            lscreen, reset_file_name)
            ENDIF
         ELSE                                         ! V3FIT not running here
            IF (ns_old .le. nsval)
     1         CALL initialize_radial(nsval, ns_old, delt0r, 
     2            lscreen, reset_file_name)
         ENDIF
#if defined (SKS)
         CALL second0(skstoff)
         init_radial_time = init_radial_time + (skstoff-skston)
!         IF(grank.GE.nranks) THEN
!           CYCLE
!         END IF
         CALL Initialize_bst(.FALSE.,nsval,blklength)
#endif
!     CONTROL NUMBER OF STEPS
         IF (numsteps > 0) THEN
            niter_store = niter
            niter = numsteps+iter2-1
         END IF

#if defined(SKS)
         CALL second0(skston)
#endif
         CALL eqsolve (ier_flag, lscreen)
#if defined(SKS)
         CALL second0(skstoff)
         eqsolve_time = eqsolve_time + (skstoff-skston)
#endif
         IF (numsteps > 0) THEN
            niter = niter_store
         END IF

         IF (imatch_phiedge .eq. 3) imovephi = 1
         IF (ier_flag.ne.norm_term_flag .and. 
     1       ier_flag.ne.successful_term_flag .and. 
     2       ier_flag.ne.more_iter_flag) EXIT
         IF (numsteps>0 .or. ns_index>0) EXIT

!
! give up if it refuses to converge, M.Drevlak
! it may help to end a vmec run in an optimization environment, if it
! fails to converge in the first iterations of an ns_array sequence
! within the set number of iterations specified by NITER.
! The parameter fgiveup defaults to 30.
!

         IF (lgiveup .and. (fsqr.gt.ftolv*fgiveup .or. 
     1                      fsqz.gt.ftolv*fgiveup .or.
     2                      fsql.gt.ftolv*fgiveup     )) THEN
            print *, "runvmec: giving up due to poor convergence"
            EXIT
         ENDIF

#if defined(SKS)
         CALL Finalize_bst(.FALSE.)
         CALL second0(gridtoff)
         grid_time(grid_id) = gridtoff - gridton
         IF (lfreeb) THEN
           IF (PARVMEC) THEN
             vgrid_time(grid_id) = vacuum_time - old_vacuum_time
             old_vacuum_time = vacuum_time
           ELSE
             vgrid_time(grid_id) = s_vacuum_time - old_vacuum_time
             old_vacuum_time = s_vacuum_time
           END IF
         END IF
         grid_id = grid_id + 1
#endif
      END DO
      !END - Main loop that will be parallelized - SKS

  100 CONTINUE

      IF (ier_flag.eq.bad_jacobian_flag .and. jacob_off.eq.0) THEN
         jacob_off = 1
         GO TO 50
      END IF

      CALL second0 (timeoff)
      timer(tsum) = timer(tsum) + timeoff - timeon
!
!     WRITE OUTPUT TO THREED1, WOUT FILES; FREE MEMORY ALLOCATED GLOBALLY
!
 1000 IF(lmoreiter .AND. ier_flag.EQ.more_iter_flag) THEN  ! J Geiger
        PRINT *, "runvmec: Running some more iterations",
     &           " -> Skipping call to fileout!"
      ELSE IF (ier_flag.NE.more_iter_flag) THEN
#if defined(SKS)
        CALL second0(skston)
#endif
        IF (PARVMEC) THEN 
          CALL fileout_par (iseq_count, ictrl_flag, ier_flag, lscreen)
        ELSE
          CALL fileout (iseq_count, ictrl_flag, ier_flag, lscreen)
        END IF
#if defined(SKS)
        CALL second0(skstoff)
        fileout_time = fileout_time + (skstoff-skston)
#endif      
      ENDIF
#if defined(SKS)
      IF(LV3FITCALL) CALL FinalizeRunVmec(RUNVMEC_COMM_WORLD)
      CALL second0(rvtoff)
      runvmec_time = runvmec_time + (rvtoff - rvton)
#endif
      END SUBROUTINE runvmec
!------------------------------------------------
