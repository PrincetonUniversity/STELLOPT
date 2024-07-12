      SUBROUTINE initialize_opt (xc, nxc, nvar, numargs)
      USE stel_constants
      USE boundary, ONLY: nfp, nedge, nuv
      USE vmec_input, ONLY: nfp_vmec => nfp
      USE modular_coils
      USE saddle_coils
      USE saddle_surface
      USE bcoils_mod
      USE tf_coils
      USE Vcoilpts
      USE coilsnamin
      USE Vname
      USE Vwire
      USE coils
      USE control_mod
      USE gade_mod
      USE mpi_params                                         !mpi stuff
      USE safe_open_mod
      USE system_mod, ONLY: chdir, system
      USE biotsavart
      IMPLICIT NONE
#if defined(MPI_OPT)
      include 'mpif.h'                                       !mpi stuff
#endif
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: nxc, nvar
      REAL(rprec) :: xc(nxc)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: num_files=6
      INTEGER :: ierr, ierr_vmec, i, k, numargs
      INTEGER :: iunit=24, iunit_vmec
      INTEGER :: nvariables,nextcur
      REAL(rprec) :: xvariables(nxc)
      LOGICAL :: isthere, isfor05there
      CHARACTER :: arg2*100, file_to_move*200
#if defined(WIN32)
      CHARACTER(LEN=*), PARAMETER :: copy = "copy "
      CHARACTER(LEN=*), PARAMETER :: remd = "rmdir /S /Q "
      CHARACTER(LEN=*), PARAMETER :: makd = "mkdir "
      CHARACTER(LEN=*), PARAMETER :: sep  = "\"
#else
      CHARACTER(LEN=*), PARAMETER :: copy = "/bin/cp "
      CHARACTER(LEN=*), PARAMETER :: remd = "/bin/rm -Rf "
      CHARACTER(LEN=*), PARAMETER :: makd = "/bin/mkdir -m 755 "
      CHARACTER(LEN=*), PARAMETER :: sep  = "/"
#endif
      CHARACTER(LEN=400) :: temp, scratch_dir, new_path
      LOGICAL :: isnc
!-----------------------------------------------

!     INITIALIZE CONSTANTS, ARRAYS

      CALL initialize_coilsin

      mmax_bmn = 10
      nmax_bmn = 18
      nsmid = 0
      nsodd = 0
      m_in=20
      n_in=25

      mbwires = 0
      nvar = 0
      nfp  = 0

!     Let's try something better
      isthere = .FALSE.
      INQUIRE(FILE=TRIM(input_data_file),EXIST=isthere)
      IF (isthere) THEN
         CALL safe_open(iunit, ierr, TRIM(input_data_file),'old',
     1     'formatted')
         IF (ierr .ne. 0) STOP 'error OPENING input file'
         CALL read_namelist (iunit, ierr, 'indata')
         IF (ierr .ne. 0) STOP 'error reading NAMELIST indata'
         nfp = nfp_vmec
         CALL read_namelist (iunit, ierr, 'coilsin')
         IF (ierr .ne. 0) STOP 'error reading NAMELIST coilsin'
      ELSE
         CALL safe_open(iunit, ierr, TRIM(for05_file_new),'old', 
     1        'formatted')
         IF (ierr .ne. 0) STOP 'error OPENING for05 file'
         CALL safe_open(iunit_vmec, ierr, TRIM(input_data_file),'old',
     1         'formatted')
         CALL read_namelist (iunit_vmec, ierr, 'indata')
         IF (ierr .ne. 0) STOP 'error reading NAMELIST indata'
         CLOSE(iunit_vmec)
      END IF

!     OPEN INPUT FILE AND LOOK FOR COILSIN NAMELIST IN FOR05_FILE_NEW
!      temp=TRIM(path)//for05_file_new
!      CALL safe_open(iunit, ierr, TRIM(temp),'old', 'formatted')

!     LOOK FOR VALUE OF NFP IN VMEC INPUT_DATA_FILE (IF IT IS THERE)
!      iunit_vmec = iunit+1
!      temp=TRIM(path)//input_data_file
!      CALL safe_open(iunit_vmec, ierr_vmec, TRIM(temp), 
!     1               'old','formatted')

!      IF (ierr_vmec .eq. 0) THEN
!          CALL read_namelist (iunit_vmec, ierr_vmec, 'indata')
!          IF (ierr_vmec .ne. 0) STOP 'error reading NAMELIST indata'
!          nfp = nfp_vmec
!      END IF

!     DETERMINE FROM WHICH FILE (EITHER FOR05_FILE_NEW OR INPUT_DATA_FILE) TO READ COILSIN 
!      IF (ierr .ne. 0) THEN
!         IF (ierr_vmec .eq. 0) THEN
!            iunit = iunit_vmec
!         ELSE
!            STOP '  Both for05_new and vmec input files DO not exist!'
!         END IF
!      ELSE
!         CLOSE (iunit_vmec, iostat=i)
!      END IF

      ! Just write out the namelist
      !temp = 'coilsin'
      !CALL safe_open(iunit, ierr, TRIM(temp),'old', 'formatted')
      !CALL write_coilsin(iunit,ierr)
      !CLOSE(iunit)
      !STOP   
         
!     LOAD XC WITH INITIAL GUESS

!     IF BOUNDS ARE REQUIRED, USE THE UNCONSTRAINED VARIABLE
!     X AND DEFINE Y(X) = X/SQRT(1+X**2). THEN, VARY X AND LET
!     THE TARGET (R, FOR EXAMPLE) BE WRITTEN:

!        R(X) = .5*(RMAX + RMIN) + .5*Y(X)*(RMAX - RMIN)
!      ierr = 0
!      CALL read_coils_namelist(iunit,ierr)
!      !CALL read_namelist (iunit, ierr, 'coilsin')
!      IF (ierr .NE. 0)then
!        IF (myid .EQ. master) THEN                         !mpi stuff
!          PRINT *,' ierr = ',ierr,' in initialize_opt reading coilsin'
!        END IF     !IF( myid .eq. master )                  !mpi stuff
!        STOP
!      ENDIF

!
!     parse wout file name based on extension (could be .nc, .txt file)
!
      wout_file = 'wout'
      CALL parse_extension(wout_file, extension, isnc)

      IF (lgeom_only) THEN

!     READ IN PLASMA BOUNDARY. IF Xcoilopt IS BEING CALLED FROM INSIDE THE OPTIMIZER
!     (LCOIL_GEOM = T mode), THE FINAL BOUNDARY SHAPE IS NOT KNOWN YET AT THIS POINT,
!     SO READING WOUT FILE IS INAPPROPRIATE. 

         lplasma = .FALSE.
         IF (numargs .ge. 2) THEN
            CALL getcarg (2, arg2, numargs)
            IF (arg2(1:2).EQ.'-P' .OR. arg2(1:2).EQ.'-p') lplasma=.true.
         END IF

      ELSE

         bnorm_file = 'bnorm.' // TRIM(extension)
!
!     DO ALL CALCULATIONS IN SCRATCH DIRECTORY
!
!         scratch_dir = "coilopt" // "_" // TRIM(extension)
!         new_path = " " // TRIM(scratch_dir) // TRIM(sep) 
!         IF (myid .eq. master) THEN                                          !START MPI
!            temp = remd // scratch_dir
!            CALL system(temp)
!            temp = makd // scratch_dir
!            CALL system(temp)
!
!
!     COPY INPUT FILES FROM ORIGINAL DIRECTORY
!
!            DO i = 1, num_files
!
!               move_file: SELECT CASE(i) 
!                  CASE (1)  
!                  file_to_move = for05_file
!                  CASE (2)  
!                  file_to_move = for05_file_new
!                  CASE (3)  
!                  file_to_move = input_data_file 
!                  CASE (4)  
!                  file_to_move = wout_file
!                  CASE (5)  
!                  file_to_move = bnorm_file 
!                  CASE (6)  
!                  file_to_move = bcoil_file 
!               END SELECT move_file
!
!            INQUIRE(file=TRIM(path) // file_to_move, exist=isthere)
!            IF (i .eq. 1) isfor05there = isthere
!            IF (isthere) THEN
!               temp = copy // TRIM(path) // TRIM(file_to_move) //
!     1                     TRIM(new_path) // TRIM(file_to_move)
!               CALL system(temp)
!            END IF
!
!            IF (i.eq.2 .and. .not.isthere .and. isfor05there) THEN
!               temp = copy // TRIM(path) // TRIM(for05_file) //
!     1                    TRIM(new_path) // TRIM(for05_file_new)
!               CALL system(temp)
!            END IF
!
!            END DO
!
!            k = chdir(scratch_dir)
!            IF (k .ne. 0) THEN
!               PRINT *,
!     1         'ierr = ',k,': Unable to chdir to working directory!'
!               STOP
!            END IF
!         END IF                                                          !End myid


#if defined(MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi)                         !MPI
         IF (ierr_mpi .ne. 0) STOP 
     1       'MPI_BARRIER error in COILOPT INITIALIZE_OPT'
         IF (myid .ne. master) k = chdir(scratch_dir)                       !MPI
#endif

         lplasma = .TRUE.

      END IF                                                       !End .not.lgeom_only

      
      IF (lplasma) THEN
         CALL read_wout_coilopt                                    !Did not get nfp from input file
         nedge = nuv

!     COMPUTE SURFACE NORMAL AT A FIXED NUMBER (NEDGE) OF MATCHING POINTS
!     (IN THETA=U, ZETA=V SPACE)

         CALL normal_vector
      ELSE
         IF (ierr_vmec .ne. 0) CALL read_wout_coilopt              !Did not get nfp from input file
         nedge = 0
      END IF

!     GENERATE mode NUMBERS FOR B_ERROR SPECTRUM

      IF (nfp .eq. 0) STOP 'NFP = 0 in initialize_opt (coilopt)'
     
      mnmax_bmn = (mmax_bmn+1)*(2*nmax_bmn/nfp+1)
      IF (mnmax_bmn .ne. 0) THEN
        ALLOCATE (xm_bmn(mnmax_bmn), xn_bmn(mnmax_bmn), bmn(mnmax_bmn),
     1    stat = i)
         IF( myid .eq. master ) THEN                         !mpi stuff
            IF (i .ne. 0) STOP 'allocation error in initialize_opt'
         END IF     !IF( myid .eq. master )                  !mpi stuff
        CALL get_modes (mmax_bmn, nmax_bmn, xm_bmn, xn_bmn, mnmax_bmn)
      END IF

      IF (nopt_alg .gt. 0) THEN
         CALL GA_preset
         CALL DE_preset
         CALL read_namelist (iunit, ierr, 'ga_de')
         niter_opt = ngen
         IF (ierr .ne. 0) THEN
            IF( myid .eq. master ) THEN                         !mpi stuff
               PRINT *,' ierr = ',ierr,
     1         ' in initialize_opt reading ga_de'
            END IF     !IF( myid .eq. master )                  !mpi stuff
            STOP
         END IF
      END IF

!     Initialize and allocate coil and wire-dependent arrays

      CALL allocate_wires

!     Initialize and count modular coil variables

      IF (lmodular) THEN
         CALL allocate_modular_coils

         CALL init_modular_coils (nvariables, xvariables, nfp)
         xc(1:nvariables) = xvariables(1:nvariables)
         nvar = nvar + nvariables

         CALL init_modular_currents (nvariables, xvariables)
         xc(nvar+1:nvar+nvariables) = xvariables(1:nvariables)
         nvar = nvar + nvariables
      END IF


!     Initialize and count saddle coil variables (IF lsaddle = true)

      IF (lsaddle) THEN
         CALL allocate_saddle_coils

         CALL init_saddle_coils (nvariables, xvariables, nfp)
         xc(nvar+1:nvar+nvariables) = xvariables(1:nvariables)
         nvar = nvar + nvariables

         CALL init_saddle_currents (nvariables, xvariables)
         xc(nvar+1:nvar+nvariables) = xvariables(1:nvariables)
         nvar = nvar + nvariables
      END IF

!     Initialize and count background coil variables (IF lbcoil = true)

      IF (lbcoil) THEN
         CALL read_bcoils
         IF (lbcoil) THEN
            CALL init_bg_currents (nvariables, xvariables)
            xc(nvar+1:nvar+nvariables) = xvariables(1:nvariables)
            nvar = nvar + nvariables
         END IF
      END IF

!     Initialize and count vf coil variables (IF lvf = true)

      IF (lvf) THEN
         CALL allocate_vf_coils
         IF (lvfvar) THEN
            CALL init_vf_coils (nvariables, xvariables)
            xc(nvar+1:nvar+nvariables) = xvariables(1:nvariables)
            nvar = nvar + nvariables
         END IF

         CALL init_vf_currents (nvariables, xvariables)
         xc(nvar+1:nvar+nvariables) = xvariables(1:nvariables)
         nvar = nvar + nvariables
      END IF

!     Initialize and count 1/R coil variables (IF ltfc = true)

      IF (ltfc) THEN
         CALL init_tf_coils (nvariables, xvariables)
         IF (ltfcv) THEN
            xc(nvar+1:nvar+nvariables) = xvariables(1:nvariables)
            nvar = nvar + nvariables
         END IF
      END IF

!     Initialize modular winding surface variables (IF lsurfv = true)

      IF (lsurfv) THEN
         CALL init_modular_wsurf (nvariables, xvariables)
         xc(nvar+1:nvar+nvariables) = xvariables(1:nvariables)
         nvar = nvar + nvariables
      END IF

!     Initialize saddle winding surface variables (IF lsadsfv = true)

      IF (lsadsfv) THEN
         CALL init_saddle_wsurf (nvariables, xvariables)
         xc(nvar+1:nvar+nvariables) = xvariables(1:nvariables)
         nvar = nvar + nvariables
      END IF

      !PRINT *,NVAR

      IF (nvar .gt. nxc) STOP 'nvar>nxc'

      CLOSE (iunit)

      END SUBROUTINE initialize_opt
