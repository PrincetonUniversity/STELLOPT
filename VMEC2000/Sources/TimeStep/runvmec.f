      SUBROUTINE runvmec(ictrl_array, input_file0, 
     1   lscreen, COMM_WORLD, reset_file_name)
      USE vmec_main
      USE vmec_params, ONLY: bad_jacobian_flag, more_iter_flag,
     1                       norm_term_flag, successful_term_flag,
     2                       restart_flag, readin_flag,
     3                       timestep_flag, ns_error_flag, 
     4                       reset_jacdt_flag
      USE vsvd
      USE timer_sub
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(inout), TARGET :: ictrl_array(5)
      LOGICAL, INTENT(in) :: lscreen
      CHARACTER(LEN=*), INTENT(in) :: input_file0
      CHARACTER(LEN=*), OPTIONAL :: reset_file_name
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
!             32       reset_jacdt_flag       Resets ijacobian flag and time step to delt0r
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
     1                                lscreen, reset_file_name)
         USE vmec_main
         IMPLICIT NONE
         INTEGER, INTENT(in) :: nsval
         INTEGER, INTENT(inout) :: ns_old
         CHARACTER(LEN=*), OPTIONAL :: reset_file_name
         LOGICAL, INTENT(in) :: lscreen
         REAL(rprec), INTENT(out) :: delt0
         END SUBROUTINE initialize_radial
      END INTERFACE

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

      IF (lreset) CALL reset_params

      IF (IAND(ictrl_flag, reset_jacdt_flag) .ne. 0) THEN
         ijacob = 0
         delt0r = delt
      END IF

      IF (IAND(ictrl_flag, readin_flag) .ne. 0) THEN
!
!        READ INPUT FILE (INDATA NAMELIST), MGRID_FILE (VACUUM FIELD DATA)
!
         CALL vsetup (iseq_count)
         CALL readin (input_file, iseq_count, ier_flag, lscreen)
         IF (ier_flag .ne. 0) GOTO 1000
!
!        COMPUTE NS-INVARIANT ARRAYS
!
         CALL fixaray

      END IF

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
         WRITE (nthreed, 30)
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

  50  CONTINUE
      iequi = 0
      IF (lfreeb .and. jacob_off.eq.1) ivac = 1    !!restart vacuum calculations

      ns_min = 3

      ITERATIONS: DO igrid = igrid0-jacob_off, multi_ns_grid
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
         
!  JDH 2012-06-20. V3FIT fix, inserted with change from VMEC 8.48 -> 8.49
!    (Not sure just what in initialize_radial messes up convergence - happens slowly)
!  Logical l_v3fit is declared in vmec_input, available via vmec_main
         IF (l_v3fit) THEN                            ! V3FIT is running here 
            IF (ns_old .ne. nsval) THEN
               CALL initialize_radial(nsval, ns_old, delt0r,
     1                                lscreen, reset_file_name)
            ENDIF
         ELSE                                         ! V3FIT not running here
            IF (ns_old .le. nsval)
     1         CALL initialize_radial(nsval, ns_old, delt0r, 
     2                                lscreen, reset_file_name)
         ENDIF

!     CONTROL NUMBER OF STEPS
         IF (numsteps > 0) THEN
            niter_store = niter
            niter = numsteps+iter2-1
         END IF

         CALL eqsolve (ier_flag, lscreen)

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

      END DO ITERATIONS

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
 1000 IF(ier_flag.eq.more_iter_flag) then  ! J Geiger
! 1000 if(lmoreiter .and. ier_flag.eq.more_iter_flag) then  ! J Geiger
        print *, "runvmec: Running some more iterations",
     &           " -> Skipping call to fileout!"
      ELSE
        CALL fileout (iseq_count, ictrl_flag, ier_flag, lscreen)
      ENDIF

      END SUBROUTINE runvmec
