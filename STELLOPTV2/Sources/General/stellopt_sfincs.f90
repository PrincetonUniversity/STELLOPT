!-----------------------------------------------------------------------
!     Subroutine:    stellopt_sfincs
!     Authors:       Matt Landreman, University of Maryland
!     Date:          August 2017
!     Description:   This subroutine obtains the boostrap current by calling the SFINCS code.
!                    The results from this calculation are then used by stellopt_vboot.f90
!                    and/or chisq_bootstrap.f90.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_sfincs(lscreen,iflag)
!DEC$ IF DEFINED (SFINCS)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime, twopi => pi2
      USE stellopt_input_mod
      USE stellopt_vars
      USE equil_utils

      USE sfincs_main, only: sfincs_init, sfincs_prepare, sfincs_run
      USE globalVariables, only: sfincs_inputFilename => inputFilename, sfincs_outputFilename => outputFilename, equilibriumFile, FSABjHat, FSABHat2, dbootstrapdlambda, ms_sensitivity, ns_sensitivity, nmodesadjoint
      USE equil_vals, only: phiedge
!!$      ! BOOTSJ LIBRARIES
!!$      USE bootsj_input
!!$      use parambs, lscreen_bootsj=>lscreen
!!$      use vmec0
!!$      use read_boozer_mod
!!$      use trig
      ! VMEC

      USE, intrinsic :: iso_fortran_env, only : stdout=>output_unit, stderr=>error_unit
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
!DEC$ ENDIF
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!DEC$ IF DEFINED (SFINCS)
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      integer, parameter :: nfax = 13
      INTEGER ::  ier, i, j, ik, ntheta, nzeta, nlis
      INTEGER ::  dex_ion, dex_zeff
      INTEGER :: mystart,myend, chunk, numprocs_local
      INTEGER, ALLOCATABLE :: mnum(:)
      REAL(rprec) :: t1, t2
      integer, dimension(nfax) :: ifaxu, ifaxv
      integer :: ntrigu, ntrigv, i1
      integer :: ntheta_min, nzeta_min, ir, mn, m, n, imn
      integer :: irho, irho1, ierr, iunit, ijbs, ians, ians_plot
      real(rprec), dimension(:), allocatable :: cputimes
      real(rprec) :: time1, timecpu, unit, file, status, err,&
         time2, r, x, al31t, gradbs1, gradbs2,&
         gradbs3, gradbs4,  al31s, aibstot, a
      real(rprec) :: a1, tempe0, tempi0, pres10, pres0
      real(rprec), dimension(:), allocatable :: work
      integer :: ihere = 0
      real(rprec), PARAMETER :: one=1
      CHARACTER(LEN=32) :: temp_str
      integer :: numProcs_world, numProcs_myworld, myRank_world, myRank_myworld, myRank_sfincs, numProcs_sfincs ! MJL
      integer :: Nradii, radius_index_min, radius_index_max, radius_index, color, key
      integer :: MPI_COMM_SFINCS
      integer, parameter :: buffer_length = 1000
      character(len=buffer_length) :: proc_assignments_string, base_directory_string, directory_string, file_line, file_line_lower, file_line_lower2
      character(len=buffer_length) :: working_directory
      integer :: tag, dummy(1), file_status, unit_in, unit_out
      integer :: mpi_status(MPI_STATUS_SIZE)
      LOGICAL :: lfile_check
      LOGICAL :: added_scanType
      REAL(rprec) :: sfincs_ne, sfincs_ni, sfincs_Te, sfincs_Ti, sfincs_d_ne_d_s, sfincs_d_ni_d_s, sfincs_d_Te_d_s, sfincs_d_Ti_d_s, delta_s
      INTEGER :: numProcs_eff, myRank_eff
      REAL(rprec) :: temp1, temp2
			INTEGER :: sfincs_nMaxAdjoint=0
			INTEGER :: sfincs_mMaxAdjoint=0
!!$      REAL(rprec) :: temp1, temp2, ds_fine, scale_factor
!!$      REAL(rprec), DIMENSION(:), ALLOCATABLE :: sfincs_s_with_0
!!$      INTEGER, PARAMETER :: Ns_fine = 1000
!!$      REAL(rprec), DIMENSION(Ns_fine) :: s_fine_full, s_fine_half, J_dot_B_flux_surface_average_fine, d_p_d_s_fine_half
!!$      REAL(rprec), DIMENSION(Ns_fine) :: B_squared_flux_surface_average_fine_half, B_squared_flux_surface_average_fine_full
!!$      REAL(rprec), DIMENSION(Ns_fine) :: integrating_factor_full, integrating_factor_half, d_p_d_s_fine, integrand, isigng_I_F_full, isigng_I_full
!!$      REAL(rprec), DIMENSION(Ns_fine) :: integrating_factor_half_approximate, pressure_fine_half
!!$      REAL(rprec), DIMENSION(21) :: J_dot_B_flux_surface_average_fit, B_squared_flux_surface_average_fit, sfincs_AC_fit
!!$      REAL(rprec), DIMENSION(Ns_fine) :: sfincs_ac_half, sfincs_ac_low_beta_limit, sfincs_ac_fit_results
!!$      CHARACTER(LEN=32) :: B_squared_flux_surface_average_profile_type

!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      real(rprec) , EXTERNAL :: al31
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (lscreen) WRITE(6,'(a)') ' ---------------------  BOOTSTRAP CALCULATION USING SFINCS  ------------------'
      SELECT CASE(TRIM(equil_type))
         CASE('vmec2000','animec','flow','satire','parvmec','paravmec','vboot')

            ! Begin code by MJL
            base_directory_string = 'sfincs_'//TRIM(proc_string)
            !IF (myworkid == master) CALL SYSTEM('mkdir ' //TRIM(base_directory_string)//' 2>/dev/null') ! The 2>/dev/null is used because otherwise there are messages printed to the screen about the directory already existing, which users need not worry about.
            IF (myworkid == master) CALL execute_command_line('mkdir -p ' //TRIM(base_directory_string)) ! The -p blocks a warning message if the directory already exists.
            !print *,"Calling sfincs at the following radii:",sfincs_s
            CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs_world, ierr_mpi)
            CALL MPI_COMM_RANK(MPI_COMM_WORLD, myRank_world, ierr_mpi)
            CALL MPI_COMM_SIZE(MPI_COMM_MYWORLD, numProcs_myworld, ierr_mpi)
            CALL MPI_COMM_RANK(MPI_COMM_MYWORLD, myRank_myworld, ierr_mpi)
            !print "(a,i3,a,i3,a,i3,a,i3,a,i3,a,i3)", "World: rank",myRank_world," of",numProcs_world,". Myworld: rank",myRank_myworld," of",numProcs_myworld,". myworkid=",myworkid,". master=",master
            IF (sfincs_min_procs > numProcs_myworld) THEN
               IF (myworkid==master) THEN
                  PRINT *,"WARNING! The number of procs in MPI_COMM_MYWORLD is smaller than sfincs_min_procs."
                  PRINT *,"# of procs in MPI_COMM_MYWORLD:",numProcs_myWorld
                  PRINT *,"sfincs_min_procs:",sfincs_min_procs
                  PRINT *,"Lowering sfincs_min_procs to",numProcs_myWorld
               END IF
               sfincs_min_procs = numProcs_myWorld
            END IF
            ! Here sfincs_s is initialized to -1
            Nradii = MINLOC(sfincs_s(2:),DIM=1)
            IF (myworkid==master .and. lscreen) PRINT *,"Number of radii for SFINCS:",Nradii
            DO j=1,Nradii
               IF (sfincs_s(j) <= 0) STOP "Error! sfincs_s values must be > 0 (except for 0 values at the end of the array, which are ignored)."
               IF (sfincs_s(j) > 1) STOP "Error! sfincs_s values must be <= 1."
            END DO

            ! Divide up MPI_COMM_MYWORLD over radii:
            ! These next lines did work, but neglect sfincs_min_procs:
!!$            radius_index_min = 1 + (myRank_myWorld*Nradii)/numProcs_myWorld
!!$            if (Nradii >= numProcs_myWorld) then
!!$               radius_index_max = ((myRank_myWorld+1)*Nradii)/numProcs_myWorld
!!$            else
!!$               radius_index_max = radius_index_min
!!$            end if

            numProcs_eff = numProcs_myWorld / sfincs_min_procs
            myRank_eff = myRank_myWorld / sfincs_min_procs
            ! If the # of procs is not an integer multiple of sfincs_min_procs, have the 'extra' procs join the last group:
            if (myRank_eff >= numProcs_eff) myRank_eff = numProcs_eff - 1
            
            radius_index_min = 1 + (myRank_eff*Nradii)/numProcs_eff
            if (Nradii >= numProcs_eff) then
               radius_index_max = ((myRank_eff+1)*Nradii)/numProcs_eff
            else
               ! There are more procs than radii.
               radius_index_max = radius_index_min
            end if
            
            ! Validation
            if (radius_index_max > Nradii) STOP "ERROR! radius_index_max > Nradii"
            if (myRank_myWorld==0 .and. radius_index_min .ne. 1) STOP "Error! proc 0 does not start with radius 1"
            if (Nradii==1 .and. radius_index_min .ne. 1) &
                 STOP "Error! When Nradii=1, all procs should start with radius 1"
            if (Nradii==1 .and. radius_index_max .ne. Nradii) &
                 STOP "Error! When Nradii=1, all procs should end with radius Nradii"
            if (radius_index_max < radius_index_min) STOP "ERROR! radius_index_max < radius_index_min"

            color = radius_index_min ! This choice ensures that all procs handling the same radii get grouped into one MPI communicator.
            key = myworkid
!              CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, color, key, MPI_COMM_SFINCS, ierr_mpi)
            CALL MPI_COMM_SPLIT(MPI_COMM_MYWORLD, color, key, MPI_COMM_SFINCS, ierr_mpi)
            CALL MPI_COMM_SIZE(MPI_COMM_SFINCS, numProcs_sfincs, ierr_mpi)
            CALL MPI_COMM_RANK(MPI_COMM_SFINCS, myRank_sfincs, ierr_mpi)
            WRITE (proc_assignments_string,fmt="(a,i5,a,i5,a,i4,a,i4,a,i5,a,i5,a)"),"Proc ",myRank_myWorld," of",numProcs_myWorld," in MPI_COMM_MYWORLD will handle radius",radius_index_min," to",radius_index_max," and has rank",myRank_sfincs," of",numProcs_sfincs," in MPI_COMM_SFINCS."

            ! Print the processor/radius assignments in a coordinated manner.
            dummy = 0
            tag = 0
            IF (myworkid==master) THEN
               IF (lscreen) WRITE(*,"(a)") TRIM(proc_assignments_string)
               DO i = 1,numProcs_myWorld - 1
                  ! To avoid a disordered flood of messages to the masterProc,
                  ! ping each proc 1 at a time by sending a dummy value:
                  CALL MPI_SEND(dummy,1,MPI_INT,i,tag,MPI_COMM_MYWORLD,ierr_mpi)
                  ! Now receive the message from proc i:
                  CALL MPI_RECV(proc_assignments_string,buffer_length,MPI_CHAR,i,MPI_ANY_TAG,MPI_COMM_MYWORLD,mpi_status,ierr_mpi)
                  IF (lscreen) WRITE(*,"(a)") TRIM(proc_assignments_string)
               END DO
            ELSE
               ! First, wait for the dummy message from proc 0:
               CALL MPI_RECV(dummy,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_MYWORLD,mpi_status,ierr_mpi)
               ! Now send the message to proc 0:
               CALL MPI_SEND(proc_assignments_string,buffer_length,MPI_CHAR,0,tag,MPI_COMM_MYWORLD,ierr_mpi)
            END IF
            
            ! We at least need +1 for a point at s=0. We might also want another +1 if we wish to impose <j dot B>=0 at s=1.
            IF (.not. ALLOCATED(sfincs_J_dot_B_flux_surface_average))   ALLOCATE(sfincs_J_dot_B_flux_surface_average(Nradii+2))
            IF (.not. ALLOCATED(sfincs_B_squared_flux_surface_average)) ALLOCATE(sfincs_B_squared_flux_surface_average(Nradii+2))
            sfincs_J_dot_B_flux_surface_average = 0
            sfincs_B_squared_flux_surface_average = 0
            CALL GET_ENVIRONMENT_VARIABLE('PWD',working_directory)
            DO radius_index = radius_index_min, radius_index_max
               WRITE (temp_str, fmt="(f20.5)") sfincs_s(radius_index)
               directory_string = TRIM(base_directory_string)//'/psiN_'//ADJUSTL(TRIM(temp_str))
               !print *,directory_string
               ier=0 ! This variable needs to be initialized to 0. Otherwise the get_equil_* routines randomly quit.

               ! Evaluate density, temperature, and their gradients at the flux surface of interest:
               ! For now, use finite differencing to get gradients.
               delta_s = 1.0d-4
               CALL get_equil_ne(sfincs_s(radius_index) + delta_s,TRIM(ne_type),temp1,ier)
               CALL get_equil_ne(sfincs_s(radius_index) - delta_s,TRIM(ne_type),temp2,ier)
               sfincs_d_ne_d_s = (temp1-temp2) / (2*delta_s)

               CALL get_equil_Te(sfincs_s(radius_index) + delta_s,TRIM(Te_type),temp1,ier)
               CALL get_equil_Te(sfincs_s(radius_index) - delta_s,TRIM(Te_type),temp2,ier)
               sfincs_d_Te_d_s = (temp1-temp2) / (2*delta_s)

               CALL get_equil_Ti(sfincs_s(radius_index) + delta_s,TRIM(Ti_type),temp1,ier)
               CALL get_equil_Ti(sfincs_s(radius_index) - delta_s,TRIM(Ti_type),temp2,ier)
               sfincs_d_Ti_d_s = (temp1-temp2) / (2*delta_s)

               CALL get_equil_ne(sfincs_s(radius_index),TRIM(ne_type),sfincs_ne,ier)
               CALL get_equil_Te(sfincs_s(radius_index),TRIM(Te_type),sfincs_Te,ier)
               CALL get_equil_Ti(sfincs_s(radius_index),TRIM(Ti_type),sfincs_Ti,ier)
               ! Convert density and temperature to SFINCS normalization.
               sfincs_ne = sfincs_ne / (1.0d+20)
               sfincs_d_ne_d_s = sfincs_d_ne_d_s / (1.0d+20)
               sfincs_ni = sfincs_ne
               sfincs_d_ni_d_s = sfincs_d_ne_d_s
               ! Stellopt temperature normalization is 1 
               sfincs_Te = sfincs_Te / 1000
               sfincs_Ti = sfincs_Ti / 1000
               sfincs_d_Te_d_s = sfincs_d_Te_d_s / 1000
               sfincs_d_Ti_d_s = sfincs_d_Ti_d_s / 1000

               IF (myRank_sfincs == 0) THEN
                  CALL SYSTEM('mkdir -p '//TRIM(directory_string)) ! -p mutes the warning printed if the directory already exists.

                  ! Copy the SFINCS input file into the new directory
                  unit_in = 11
                  unit_out = 12
                  OPEN (unit=unit_in, file='input.'//TRIM(id_string),status='old',action='read',iostat=file_status)
                  IF (file_status /= 0) THEN
                     PRINT *,"Error opening stellopt input file input."//TRIM(id_string)
                     PRINT *,"iostat=",file_status
                     STOP
                  END IF
                  OPEN (unit=unit_out, file=TRIM(directory_string)//'/input.namelist',action='write',iostat=file_status)
                  IF (file_status /= 0) THEN
                     PRINT *,"Error opening new sfincs input.namelist file. iostat=",file_status
                     STOP
                  END IF
                  added_scanType = .false.
                  DO
                     file_line = ''
                     READ(UNIT=unit_in,FMT='(A)',IOSTAT=file_status) file_line
                     IF (file_status .lt. 0) EXIT ! End of file
                     file_line_lower = file_line
                     CALL tolower(file_line_lower)
                     file_line_lower = ADJUSTL(TRIM(file_line_lower))

                     IF (TRIM(file_line_lower)=='&geometryparameters') THEN
                        WRITE(UNIT=unit_out,FMT='(A)') '&geometryParameters'
												WRITE(UNIT=unit_out,FMT='(A)') '  inputRadialCoordinate = 1'
												WRITE(UNIT=unit_out,FMT='(A)') '  inputRadialCoordinateForGradients = 1'
												WRITE(UNIT=unit_out,FMT='(a, es24.14)') '  psiN_wish = ',sfincs_s(radius_index)
												IF (ANY(lsfincs_boozer_bmnc_opt)) THEN
													WRITE(UNIT=unit_out,FMT='(A)') '  geometryScheme = 13'
													WRITE(UNIT=unit_out,FMT='(a, es24.14)') '  iota = ',sfincs_iota(radius_index)
													WRITE(UNIT=unit_out,FMT='(a, es24.14)') '  GHat = ',sfincs_GHat(radius_index)
													WRITE(UNIT=unit_out,FMT='(a, es24.14)') '  IHat = ',sfincs_IHat(radius_index)
													WRITE(UNIT=unit_out,FMT='(a, es24.14)') '  aHat = ',sfincs_aHat(radius_index)
													WRITE(UNIT=unit_out,FMT='(a, es24.14)') '  psiAHat = ',sfincs_psiAHat(radius_index)
													WRITE(UNIT=unit_out,FMT='(A, I3)') '  NPeriods = ',sfincs_nperiods
													! Write bmnc's for this surface
													DO m=0,sfincs_mmax
														DO n=-sfincs_nmax,sfincs_nmax
															IF (sfincs_boozer_bmnc(m,n,radius_index)/=0) THEN
																WRITE(UNIT=unit_out,FMT='(5X,A,I1,A,I1,A,es24.14)') 'BOOZER_BMNC(',m,',',n,') = ',sfincs_boozer_bmnc(m,n,radius_index)
															END IF
														END DO
													END DO
												ELSE
													WRITE(UNIT=unit_out,FMT='(A)') '  geometryScheme = 5'
													WRITE(UNIT=unit_out,FMT='(A)') '  VMECRadialOption = 0'
													WRITE(UNIT=unit_out,FMT='(A)') '  equilibriumFile = "'//TRIM(working_directory)//'/wout_'//TRIM(proc_string)//'.nc"'
												END IF
                        CYCLE
                     END IF

                     IF (TRIM(file_line_lower)=='&speciesparameters') THEN
                        WRITE (UNIT=unit_out,FMT='(A)') '&speciesParameters'
                        WRITE (UNIT=unit_out,FMT='(a, es24.14, es24.14, a)') '  nHats = ',sfincs_ne, sfincs_ni,' ! From stellopt'
                        WRITE (UNIT=unit_out,FMT='(a, es24.14, es24.14, a)') '  THats = ',sfincs_Te, sfincs_Ti,' ! From stellopt'
                        WRITE (UNIT=unit_out,FMT='(a, es24.14, es24.14, a)') '  dnHatdpsiNs = ',sfincs_d_ne_d_s, sfincs_d_ni_d_s,' ! From stellopt'
                        WRITE (UNIT=unit_out,FMT='(a, es24.14, es24.14, a)') '  dTHatdpsiNs = ',sfincs_d_Te_d_s, sfincs_d_Ti_d_s,' ! From stellopt'
                        CYCLE
                     END IF

										 IF (lsfincs_bootstrap_analytic) THEN
												IF (TRIM(file_line_lower)=='&general') THEN
													WRITE (UNIT=unit_out,FMT='(A)') '&general'
													WRITE (UNIT=unit_out,FMT='(A)') '   RHSMode = 4'
												ELSE IF (TRIM(file_line_lower)=='&sensitivityoptions') THEN
													DO m = 0, sfincs_mmax
														DO n = -sfincs_nmax,sfincs_nmax
															IF (lsfincs_boozer_bmnc_opt(m,n,radius_index)) THEN
																IF (m > sfincs_mMaxAdjoint) THEN
																	sfincs_mMaxAdjoint = m
																END IF
																IF (ABS(n) > sfincs_nMaxAdjoint) THEN
																	sfincs_nMaxAdjoint = ABS(n)
																END IF
															END IF
														END DO
													END DO
													WRITE (UNIT=unit_out,FMT='(A)') '&sensitivityOptions'
													WRITE (UNIT=unit_out,FMT='(A)') '   adjointBootstrapOption = .true.'
													WRITE (UNIT=unit_out,FMT='(A)') '		nMinAdjoint = 0'
													WRITE (UNIT=unit_out,FMT='(A,I3)') '		nMaxAdjoint = ',sfincs_nMaxAdjoint
													WRITE (UNIT=unit_out,FMT='(A)') '		mMinAdjoint = 0'
													WRITE (UNIT=unit_out,FMT='(A,I3)') '		mMaxAdjoint = ',sfincs_mMaxAdjoint
												END IF
										 END IF

                     ! Change any parameters as needed for this particular radius.
                     IF (file_line_lower(1:3)=='!ss') THEN
                        file_line_lower2 = ADJUSTL(file_line_lower(4:))
                        IF (file_line_lower2(1:8)=='scantype') CYCLE ! Remove any previous scanType setting.
                     ELSE
                        IF ((.not. file_line_lower(1:1)=='!') .and. (.not.added_scanType)) THEN
                           added_scanType = .true.
                           WRITE(UNIT=unit_out,FMT='(A)') '!ss scanType = 4' ! Set scanType=4 on the first non-comment line.
                        END IF
                     END IF

                     ! Handle variables that should NOT be copied:
                     IF (file_line_lower(1:7)=='rhsmode') CYCLE
										 IF (file_line_lower(1:8)=='&general') CYCLE
										 IF (file_line_lower(1:19)=='&sensitivityoptions') CYCLE
                     IF (file_line_lower(1:14)=='geometryscheme') CYCLE
                     IF (file_line_lower(1:15)=='equilibriumfile') CYCLE
                     IF (file_line_lower(1:16)=='vmecradialoption') CYCLE
                     IF (file_line_lower(1:22)=='inputradialcoordinate ' .or. file_line_lower(1:22)=='inputradialcoordinate=') CYCLE
                     IF (file_line_lower(1:33)=='inputradialcoordinateforgradients') CYCLE
                     IF (file_line_lower(1: 9)=='psin_wish') CYCLE
                     IF (file_line_lower(1:11)=='psihat_wish') CYCLE
                     IF (file_line_lower(1: 7)=='rn_wish') CYCLE
                     IF (file_line_lower(1: 9)=='rhat_wish') CYCLE
                     IF (file_line_lower(1: 5)=='nhats') CYCLE
                     IF (file_line_lower(1: 5)=='thats') CYCLE
                     IF (file_line_lower(1:11)=='dnhatdpsins') CYCLE
                     IF (file_line_lower(1:11)=='dthatdpsins') CYCLE
                     IF (file_line_lower(1:13)=='dnhatdpsihats') CYCLE
                     IF (file_line_lower(1:13)=='dthatdpsihats') CYCLE
                     IF (file_line_lower(1: 9)=='dnhatdrns') CYCLE
                     IF (file_line_lower(1: 9)=='dthatdrns') CYCLE
                     IF (file_line_lower(1:11)=='dnhatdrhats') CYCLE
                     IF (file_line_lower(1:11)=='dthatdrhats') CYCLE
                     IF (file_line_lower(1: 3)=='er ' .or. file_line_lower(1:3)=='er=') CYCLE
                     IF (file_line_lower(1:10)=='dphihatdrn') CYCLE
                     IF (file_line_lower(1:12)=='dphihatdrhat') CYCLE
                     IF (file_line_lower(1:14)=='dphihatdpsihat') CYCLE
										 IF (file_line_lower(1:22)=='adjointbootstrapoption') CYCLE
										 IF (file_line_lower(1:11)=='nminadjoint') CYCLE
										 IF (file_line_lower(1:11)=='mminadjoint') CYCLE
										 IF (file_line_lower(1:11)=='nmaxadjoint') CYCLE
										 IF (file_line_lower(1:11)=='mmaxadjoint') CYCLE

                     ! If we've made it this far, copy the line from the old to new file:
                     WRITE(UNIT=unit_out,FMT='(A)',IOSTAT=file_status) TRIM(file_line)
                     IF (file_status .gt. 0) PRINT *,'Error writing new sfincs input.namelist file. IOSTAT=',file_status
                  END DO
                  CLOSE (unit_in)
                  CLOSE (unit_out)

                  ! For use with sfincsScanPlot, copy one of the sfincs input namelists to the base directory for the radial scan:
                  CALL copy_txtfile(TRIM(directory_string)//'/input.namelist', TRIM(base_directory_string)//'/input.namelist')
               END IF ! If myRank_sfincs==0
               CALL MPI_BARRIER(MPI_COMM_SFINCS,ierr_mpi) ! Don't let slave procs proceed until master has finished writing the input file.

               ! SFINCS prints lots to stdout. Redirect stdout and stderr to files.
               OPEN(unit=stdout, file=TRIM(directory_string)//"/stdout")
               OPEN(unit=stderr, file=TRIM(directory_string)//"/stderr")

               !CALL SYSTEM('cd '//TRIM(directory_string))
               sfincs_inputFilename  = TRIM(directory_string)//"/input.namelist"
               sfincs_outputFilename = TRIM(directory_string)//"/sfincsOutput.h5"
               ! Call SFINCS:
               CALL sfincs_init(MPI_COMM_SFINCS)
               !equilibriumFile = 'wout_'//TRIM(proc_string)//'.nc'
               CALL sfincs_prepare()
               CALL sfincs_run()
               CALL sfincs_finalize()
               !CALL SYSTEM('cd ../..')
               ! Return stdout and stderr to normal:
               CLOSE(stdout)
               CLOSE(stderr)
               OPEN(unit=stdout, file='/dev/stdout')
               OPEN(unit=stderr, file='/dev/stderr')

               IF (myRank_sfincs==0) THEN
                  ! Convert <j dot B> from sfincs normalization to SI units (Tesla Amperes / m^2):
                  sfincs_J_dot_B_flux_surface_average(radius_index+1) = FSABjHat * 437695 * 1e20 * 1.602177e-19 ! vBar * nBar * e
                  ! Above, the +1 is so there is a point with <j dot B>=0 at s=0.
                  sfincs_B_squared_flux_surface_average(radius_index+1) = FSABHat2
									IF (lsfincs_bootstrap_analytic) THEN
										DO m = 0, sfincs_mmax
											DO n = -sfincs_nmax,sfincs_nmax
												DO imn = 1,NModesAdjoint
													IF (ns_sensitivity(imn)==n .and. ms_sensitivity(imn)==m) THEN
													! First dimension of dbootstrapdlambda is NLambdas (=1 for boozer)
														sfincs_dBootstrapdBmnc(m,n,radius_index) = dbootstrapdlambda(1,imn) * 437695 * 1e20 * 1.602177e-19
													END IF
												END DO
											END DO
										END DO
										IF (ALLOCATED(dbootstrapdlambda)) DEALLOCATE(dbootstrapdlambda)
										IF (ALLOCATED(ns_sensitivity)) DEALLOCATE(ns_sensitivity)
										IF (ALLOCATED(ms_sensitivity)) DEALLOCATE(ms_sensitivity)
									END IF
               END IF

            END DO ! Loop over radii

            ! Send results from all procs to master so master can compute the new radial current profile:
            IF (myworkid == master) THEN

               CALL MPI_REDUCE(MPI_IN_PLACE,sfincs_J_dot_B_flux_surface_average,Nradii+1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,sfincs_B_squared_flux_surface_average,Nradii+1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               ! Extrapolate to get <B^2> on axis:
               IF (Nradii<2) THEN
                  sfincs_B_squared_flux_surface_average(1) = sfincs_B_squared_flux_surface_average(2) ! If only 1 radius, you can't linearly extrapolate.
               ELSE
                  sfincs_B_squared_flux_surface_average(1) = sfincs_B_squared_flux_surface_average(2) &
                       + (sfincs_B_squared_flux_surface_average(3)-sfincs_B_squared_flux_surface_average(2)) * (-sfincs_s(1)) / (sfincs_s(2) - sfincs_s(1))
               END IF
               !PRINT *,"Final <j dot B> profile from sfincs (in Tesla Amperes / m^2):",sfincs_J_dot_B_flux_surface_average(1:(Nradii+1))
							 IF (TRIM(bootcalc_type)=='sfincs') THEN
								 ierr=0
								 CALL safe_open(unit_out,ierr,'sfincs_results_before_profile_fitting.'//TRIM(proc_string),'REPLACE','formatted')
								 WRITE (unit_out,"(a)") "s (normalized toroidal flux), <j dot B> (Tesla Amperes / meters^2), <B^2> (Tesla^2)"
								 WRITE (unit_out,*) 0.0d+0,0.0d+0,sfincs_B_squared_flux_surface_average(1)
								 DO radius_index = 1, Nradii
										WRITE (unit_out,*) sfincs_s(radius_index), sfincs_J_dot_B_flux_surface_average(radius_index+1), sfincs_B_squared_flux_surface_average(radius_index+1)
								 END DO
								 CLOSE (unit_out)
							END IF
            ELSE
               CALL MPI_REDUCE(sfincs_J_dot_B_flux_surface_average,sfincs_J_dot_B_flux_surface_average,Nradii+1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(sfincs_B_squared_flux_surface_average,sfincs_B_squared_flux_surface_average,Nradii+1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            END IF

            CALL MPI_COMM_FREE(MPI_COMM_SFINCS,ierr_mpi)
            CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi) ! Weird things might happen if some procs move on while others are still running sfincs.

!!$            ! Here we use the same convention as in VMEC: half-mesh quantities use arrays with the same size as full-mesh quantities,
!!$            ! but the first array element is 0.
!!$            s_fine_full = [( (j-1.0d+0)/(Ns_fine-1), j=1, Ns_fine )]
!!$            ds_fine = s_fine_full(2) - s_fine_full(1)
!!$            s_fine_half(1) = 0
!!$            s_fine_half(2:Ns_fine) = s_fine_full(1:(Ns_fine-1)) + ds_fine/2
!!$
!!$            ! Fit a function to <j dot B>. We will use the same functional type as bootj_type so that
!!$            ! bootj_aux_f can provide a good initial guess for <j dot B>.
!!$            ALLOCATE(sfincs_s_with_0(Nradii+1))
!!$            sfincs_s_with_0 = 0
!!$            sfincs_s_with_0(2:) = sfincs_s(1:Nradii)
!!$            ! For an initial guess at the fit coefficients for <j dot B>, use the previous dI/ds profile:
!!$            ! AC = -d I / d s = -2 pi (d psi / d s) <j dot B> / <B^2> + (small d p / d s term)
!!$            J_dot_B_flux_surface_average_fit = bootj_aux_f(1:21)
!!$            IF (myworkid==master) PRINT *,"fit before scaling:",J_dot_B_flux_surface_average_fit
!!$            scale_factor = -SUM(sfincs_B_squared_flux_surface_average(2:(Nradii+1)))/(Nradii*(-phiedge))
!!$            IF (myworkid==master) PRINT *,"Scale factor:",scale_factor
!!$            CALL scale_profile(bootj_type, J_dot_B_flux_surface_average_fit, scale_factor)
!!$            IF (myworkid==master) PRINT *,"fit after scaling:",J_dot_B_flux_surface_average_fit
!!$            CALL fit_profile(bootj_type, Nradii+1, sfincs_s_with_0, sfincs_J_dot_B_flux_surface_average, 21, &
!!$                 J_dot_B_flux_surface_average_fit)
!!$            IF (myworkid==master) PRINT *,"New fit to <j dot B>:",J_dot_B_flux_surface_average_fit
!!$            DEALLOCATE(sfincs_s_with_0)
!!$
!!$            ! Fit a profile to <B^2>:
!!$            B_squared_flux_surface_average_fit = 0
!!$            B_squared_flux_surface_average_profile_type = 'power_series'
!!$            B_squared_flux_surface_average_fit(1) = SUM(sfincs_B_squared_flux_surface_average(2:(Nradii+1))) / Nradii ! Initial guess for constant term in the polynomial fit.
!!$            B_squared_flux_surface_average_fit(2:3) = B_squared_flux_surface_average_fit(1) / 10 ! Set to a small nonzero value, so 'fit_profile' detects the correct polynomial order.
!!$            IF (myworkid==master) PRINT *,"Guess for fit to <B^2>:",B_squared_flux_surface_average_fit,"ier=",ier
!!$            CALL fit_profile(B_squared_flux_surface_average_profile_type, Nradii, sfincs_s, sfincs_B_squared_flux_surface_average(2:(Nradii+1)), 21, &
!!$                 B_squared_flux_surface_average_fit)
!!$            IF (myworkid==master) PRINT *,"New fit to <B^2>:",B_squared_flux_surface_average_fit,"ier=",ier
!!$            DO radius_index = 1,Ns_fine
!!$               CALL eval_profile(s_fine_half(radius_index), B_squared_flux_surface_average_profile_type, &
!!$                    B_squared_flux_surface_average_fine_half(radius_index), bootj_aux_s, B_squared_flux_surface_average_fit, ier)
!!$               CALL eval_profile(s_fine_full(radius_index), B_squared_flux_surface_average_profile_type, &
!!$                    B_squared_flux_surface_average_fine_full(radius_index), bootj_aux_s, B_squared_flux_surface_average_fit, ier)
!!$            END DO
!!$
!!$            ! Evaluate dp/ds:
!!$            DO radius_index = 2,Ns_fine
!!$               CALL get_equil_p(s_fine_half(radius_index), pressure_fine_half(radius_index), ier, d_p_d_s_fine_half(radius_index))
!!$            END DO
!!$            CALL get_equil_p(s_fine_full(1), temp1, ier)
!!$            d_p_d_s_fine_half(1) = 0
!!$            ! Sanity test: compute the integrating factor in the approximation that <B^2> is constant:
!!$            integrating_factor_half_approximate = exp(mu0 * (pressure_fine_half-temp1) / (SUM(sfincs_B_squared_flux_surface_average(2:(Nradii+1))) / Nradii))
!!$
!!$            ! Compute the integrating factor, eq (19):
!!$            integrand = 0
!!$            DO radius_index = 2,Ns_fine
!!$               integrand(radius_index) = integrand(radius_index-1) + ds_fine * d_p_d_s_fine_half(radius_index) / B_squared_flux_surface_average_fine_half(radius_index)
!!$            END DO
!!$            integrating_factor_full = exp(mu0 * integrand)
!!$            integrating_factor_half(1) = 0
!!$            integrating_factor_half(2:Ns_fine) = 0.5d+0 * (integrating_factor_full(1:(Ns_fine-1)) + integrating_factor_full(2:Ns_fine))
!!$
!!$            ! Form the integrand to eq (20) on the half mesh:
!!$            DO radius_index = 2,Ns_fine
!!$               CALL eval_profile(s_fine_half(radius_index), bootj_type, j_dot_B_flux_surface_average_fine(radius_index), bootj_aux_s, J_dot_B_flux_surface_average_fit, ier) ! <J dot B>
!!$            END DO
!!$            j_dot_B_flux_surface_average_fine(1) = 0
!!$            B_squared_flux_surface_average_fine_half(1) = 0
!!$            ! Here comes the integrand of eq (20). There are 2 factors of isigng that multiply together to give +1:
!!$            !   * The fact that phi = 2 pi psi isigng.
!!$            !   * The fact that AC = isigng * dI/ds
!!$            sfincs_AC_low_beta_limit = phiedge * j_dot_B_flux_surface_average_fine / B_squared_flux_surface_average_fine_half
!!$            integrand = sfincs_AC_low_beta_limit * integrating_factor_half
!!$            integrand(1) = 0
!!$
!!$            ! Integrate to get isigng*I(s)*F(s) on the full mesh:
!!$            isigng_I_F_full = 0
!!$            DO j=2,Ns_fine
!!$               isigng_I_F_full(j) = isigng_I_F_full(j-1) + integrand(j) * ds_fine
!!$            END DO
!!$
!!$            isigng_I_full = isigng_I_F_full / integrating_factor_full
!!$            ! Differentiate isigng * I(s) to get AC(s):
!!$            sfincs_AC_half = 0
!!$            DO j=2,Ns_fine
!!$               sfincs_AC_half(j) = (isigng_I_full(j) - isigng_I_full(j-1)) / ds_fine
!!$            END DO
!!$
!!$            ! Fit the resulting AC profile
!!$            sfincs_AC_fit = bootj_aux_f(1:21)
!!$            CALL fit_profile(bootj_type, Ns_fine-1, s_fine_half(2:Ns_fine), sfincs_AC_half(2:Ns_fine), 21, sfincs_AC_fit)
!!$            ! Evaluate the fit:
!!$            DO radius_index = 2,Ns_fine
!!$               CALL eval_profile(s_fine_half(radius_index), bootj_type, sfincs_AC_fit_results(radius_index), bootj_aux_s, sfincs_AC_fit, ier)
!!$            END DO
!!$
!!$
!!$            IF (myworkid == master) THEN
!!$               ierr=0
!!$               CALL safe_open(unit_out,ierr,'sfincs_results_after_profile_fitting.'//TRIM(proc_string),'REPLACE','formatted')
!!$               WRITE (unit_out,"(a)") &
!!$                    "s (normalized toroidal flux), " // &
!!$                    "Toroidal current derivative AC (Amperes), " // &
!!$                    "Low beta approximation of toroidal current derivative AC (Amperes), " // &
!!$                    "Fit to toroidal current derivative AC (Amperes), " // &
!!$                    "Fit to <j dot B> (Tesla Amperes / meters^2), " // &
!!$                    "Fit to <B^2> (Tesla^2), " // &
!!$                    "d pressure / d s (Pascals), " // &
!!$                    "Integrating factor (eq (19), dimensionless), " // &
!!$                    "Constant-<B^2> approximation for the integrating factor (dimensionless), " // &
!!$                    "Integrand of eq (20)"
!!$               DO radius_index = 2, Ns_fine
!!$                  WRITE (unit_out,*) s_fine_half(radius_index), sfincs_AC_half(radius_index), sfincs_AC_low_beta_limit(radius_index), &
!!$                       sfincs_AC_fit_results(radius_index), j_dot_B_flux_surface_average_fine(radius_index), B_squared_flux_surface_average_fine_half(radius_index), &
!!$                       d_p_d_s_fine_half(radius_index), integrating_factor_half(radius_index), integrating_factor_half_approximate(radius_index), integrand(radius_index)
!!$               END DO
!!$               CLOSE (unit_out)
!!$            END IF

            !DEALLOCATE(sfincs_J_dot_B_flux_surface_average, sfincs_B_squared_flux_surface_average)
            ! Original code by Sam follows.











!!$            ! Get the data into bootsj
!!$            call second0 (time1)
!!$            timecpu = time1
!!$            lscreen_bootsj = lscreen
!!$            ! Open Files for Output
!!$            IF (myworkid == master) THEN
!!$               CALL read_boozer('')
!!$               ians = 12
!!$               CALL safe_open(ians, ierr, 'answers.'//trim(proc_string), 'replace','formatted')
!!$               ians_plot = 14
!!$               CALL safe_open(ians_plot, ierr, 'answers_plot.'//trim(proc_string),'replace', 'formatted')
!!$               ijbs = 15
!!$               CALL safe_open(ijbs, ierr, 'jBbs.'//trim(proc_string),'replace', 'formatted')
!!$            ELSE
!!$               ians = 12
!!$               WRITE(temp_str,'(I3.3)') myworkid
!!$               CALL safe_open(ians, ierr, 'answers_'//TRIM(temp_str)//'.'//trim(proc_string),'replace', 'formatted')
!!$            END IF
!!$!DEC$ IF DEFINED (MPI_OPT)
!!$            CALL MPI_COMM_SIZE( MPI_COMM_MYWORLD, numprocs_local, ierr_mpi )
!!$            CALL MPI_BCAST(mnboz_b,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(mboz_b,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(nboz_b,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(nfp_b,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(ns_b,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(irdim,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(irup,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(periods,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(lasym_b,1,MPI_LOGICAL,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(lasym_bootsj,1,MPI_LOGICAL,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            IF (myworkid /= master) CALL allocate_radial
!!$            CALL MPI_BCAST(flux,irup,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(qsafety,irup,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(aipsi,irup,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(idx,irup,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(gpsi,irup,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(pres1,irup,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(betar,irup,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(phip,irup,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(sign_jacobian,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            IF (.not. ALLOCATED(ixm_b)) ALLOCATE(ixm_b(mnboz_b))
!!$            IF (.not. ALLOCATED(ixn_b)) ALLOCATE(ixn_b(mnboz_b))
!!$            IF (.not. ALLOCATED(bmnc_b)) ALLOCATE(bmnc_b(mnboz_b,ns_b))
!!$            IF (.not. ALLOCATED(bmns_b) .and. lasym_b) ALLOCATE(bmns_b(mnboz_b,ns_b))
!!$            CALL MPI_BCAST(ixm_b,mnboz_b,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(ixn_b,mnboz_b,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            ik = SIZE(bmnc_b)
!!$            CALL MPI_BCAST(bmnc_b,ik,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            IF (lasym_b) CALL MPI_BCAST(bmns_b,ik,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            ! Handle unbroadcast variables
!!$            dibs=0; aibs = 0; dibst = 0; aibst = 0; bsdense = 0; bsdensi = 0; bstempe = 0
!!$            bstempi = 0; bsdenste = 0; bsdensti = 0; bstempte = 0; bstempti = 0
!!$            capr = 0; caps = 0; h2 = 0; ftrapped = 0; fpassing =0; epsttok = 0
!!$            fttok = 0; gbsnorm = 0; aiterm1 = 0; other1 = 0; rhoar = 0; bsnorm = 0
!!$            fptok = 0; amain = 0; bmax1 = 0; thetamax = 0; zetahmax = 0; ajBbs = 0
!!$            d_rho = 0; b2avg = 0 
!!$
!!$            ! Divide up work
!!$            IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
!!$            ALLOCATE(mnum(numprocs_local))
!!$            mnum=0
!!$            i = 1
!!$            DO
!!$               IF (SUM(mnum,DIM=1) == irup) EXIT
!!$               IF (i > numprocs_local) i = 1
!!$               mnum(i) = mnum(i) + 1
!!$               i=i+1
!!$            END DO
!!$            mystart = 1
!!$            DO i = 1, myworkid
!!$               mystart = SUM(mnum(1:i))+1
!!$            END DO
!!$            myend = mystart + mnum(myworkid+1) - 1
!!$            IF (myend < mystart) myend = mystart
!!$            IF (mnum(myworkid+1) == 0) mystart = myend + 1
!!$            DEALLOCATE(mnum)
!!$
!!$!DEC$ ENDIF
!!$            ! Assume the bootsj namelist has been read in
!!$            if(damp_bs .lt. 0.0) then !in this case no damp_bs was read in
!!$              if(damp .gt. 0.0) then
!!$                 damp_bs = sqrt(damp)
!!$              else
!!$                 damp_bs = 0.001_dp      !in this case no damp was read in
!!$              endif
!!$            endif
!!$            ! Dimensionality Check
!!$            if(nboz_b .lt. nbuse) nbuse = nboz_b
!!$            if(mboz_b .lt. mbuse) mbuse = mboz_b
!!$            nzeta_min = 2*nbuse + 1
!!$            ntheta_min = 2*mbuse + 1
!!$            do i = 0, 6
!!$               nzetah = 4*2**i
!!$               if(nzetah .gt. nzeta_min) exit
!!$               nzetah = 2*2**i * 3
!!$               if(nzetah .gt. nzeta_min) exit
!!$            enddo
!!$            do i = 0, 6
!!$               nthetah = 4*2**i
!!$               if(nthetah .gt. ntheta_min) exit
!!$               nthetah = 2*2**i * 3
!!$               if(nthetah .gt. ntheta_min) exit
!!$            enddo
!!$            ! Convert bmn's to amnfit's
!!$            amnfit = 0.0
!!$            lsurf = .false.
!!$            status = tiny(a1)
!!$            do ir = 1, irup
!!$               do mn = 1,mnboz_b
!!$                  m = ixm_b(mn)
!!$                  n = ixn_b(mn)/nfp_b
!!$                  if (m .gt. mboz_b) stop 'boozmn indexing conflict, m'
!!$                  if (abs(n) .gt. nboz_b) stop 'boozmn indexing conflict, n'
!!$                  if (n.lt.0) then
!!$                     m = -m
!!$                     n = -n
!!$                  end if
!!$                  if (m.eq.0 .and. n.eq.0 .and. bmnc_b(mn,ir).gt.status) lsurf(ir) = .true.
!!$                  amnfit(ir,m,n) = bmnc_b(mn,ir+1)     !!2nd half grid == 1st grid pt. here
!!$                  IF (lasym_b) amnfit2(ir,m,n) = bmns_b(mn,ir+1) 
!!$               end do
!!$            end do
!!$            ! Note we don't deallocate the boozer coordinates
!!$            zeff1 = max(1.0_rprec,zeff1)
!!$            ! Setup Mesh
!!$            psimax = maxval(abs(flux))
!!$            if(flux(irup) .lt. 0.0) psimax = -psimax
!!$            do ir = 1, irup
!!$              rhoar(ir) = 0.5_dp*(flux(ir) + flux(ir+1))/psimax
!!$              d_rho(ir) = (flux(ir+1) - flux(ir))/psimax
!!$            end do
!!$            if (iotasign .lt. 0) then
!!$               qsafety(:irup) = iotasign*qsafety(:irup)
!!$            endif
!!$            aiogar(:irup) = aipsi(:irup)/(gpsi(:irup)+1.0e-36_dp)
!!$            call positiv (pres1, irup, 2) !to be sure that arrays are positive
!!$            call positiv (betar, irup, 2)
!!$            ! Now we get TI and TE on VMEC mesh (our way)
!!$            ! Notes: This part of the code wants quantities in 10^20 [m^-3] and
!!$            !        [keV]
!!$            IF (myworkid == master) THEN
!!$               IF (tempres < 0) THEN
!!$                  DO ir = 1, irup
!!$                     CALL get_equil_te(rhoar(ir),TRIM(te_type),tempe1(ir),ier)
!!$                     CALL get_equil_ne(rhoar(ir),TRIM(ne_type),dense(ir),ier)
!!$                     IF (teti <= 0.0) THEN
!!$                        CALL get_equil_ti(rhoar(ir),TRIM(ti_type),tempi1(ir),ier)
!!$                     ELSE
!!$                        tempi1(ir) = tempe1(ir)/teti
!!$                     END IF
!!$                     dense(ir) = dense(ir) + 1.E-36_dp
!!$                  END DO
!!$               ELSE
!!$                  ! Setup some variables
!!$                  tempe0 = ate(0)            !central electron temperature in keV
!!$                  tempi0 = ati(0)                 !central ion temperature in keV
!!$                  pres0 = 1.6022E-19_DP           !Normalization of P=N*ec Note we want eV and m^-3 at this point
!!$                  ! Mimic behavior (if ate/i >0 then use as central values and scale density)
!!$                  if (tempe0.le.0 .or. tempi0.le.0) tempres = abs(tempres)
!!$                  tempres = min(one,tempres)
!!$                  ! Calculate te/ti profiles
!!$                  a = one + one/(zeff1*teti)
!!$                  tempe0 = pres1(1)/(a*pres0*dens0*1E20) ! note dens0 in 1E20 m^-3
!!$                  tempi0 = tempe0/teti
!!$                  a      = tempe0/(pres1(1)**tempres)
!!$                  a1     = tempi0/(pres1(1)**tempres)
!!$                  tempe1 = pres1**tempres
!!$                  tempi1 = tempe1
!!$                  tempe1 = tempe1*a
!!$                  tempi1 = tempi1*a1
!!$                  dense = pres1/(pres0*(tempe1+tempi1/zeff1)+1.E-36_dp)
!!$               END IF
!!$            END IF
!!$!DEC$ IF DEFINED (MPI_OPT)
!!$            CALL MPI_BCAST(tempe1,irdim,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(dense,irdim,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(tempi1,irdim,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$!DEC$ ENDIF
!!$            tempe1 = tempe1/1000.         ! [eV] to [keV]
!!$            tempi1 = tempi1/1000.         ! [eV] to [keV]
!!$            dense  = dense/(1.0E+20)       ! [m^-3] to 10^20 [m^-3]
!!$            tempe0 = tempe1(1)            !central electron temperature in keV
!!$            tempi0 = tempi1(1)                 !central ion temperature in keV
!!$            pres10 = pres1(1)                    !central total pressure in Pa
!!$            pres0 = 1.6022E4_DP
!!$            call positiv (tempe1, irup, 2)
!!$            call positiv (tempi1, irup, 2)
!!$            IF (ALLOCATED(work)) DEALLOCATE(work)
!!$            allocate(work(irdim))
!!$            call smooth1 (dense, 1, irup, work, 0.0)
!!$            call positiv (dense, irup, 2)
!!$            DEALLOCATE(work)
!!$            i1 = irup - 1
!!$            a = tempe1(irup) + tempi1(irup)/zeff1
!!$            a1 = tempe1(i1) + tempi1(i1)/zeff1
!!$            dense(irup) = dense(i1)*a1*betar(irup)/(a*betar(i1)+1.E-36_dp)
!!$            dex_zeff = MINLOC(zeff_aux_s(2:),DIM=1)
!!$            IF (myworkid == master) THEN
!!$               IF (dex_zeff > 4) THEN
!!$                  DO ir = 1, irup
!!$                     CALL get_equil_zeff(rhoar(ir),TRIM(zeff_type),zeff1,ier)
!!$                     densi(ir) = dense(ir)/zeff1
!!$                  END DO
!!$                  CALL get_equil_zeff(rhoar(1),TRIM(zeff_type),zeff1,ier)
!!$               ELSE
!!$                  densi(:irup) = dense(:irup)/zeff1
!!$               END IF
!!$            END IF
!!$!DEC$ IF DEFINED (MPI_OPT)
!!$            CALL MPI_BCAST(zeff1,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            CALL MPI_BCAST(densi,irdim,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$!DEC$ ENDIF
!!$            dens0 = dense(1)                         !central electron density
!!$            ntrigu = 3*nthetah/2 + 1
!!$            ntrigv = 2*nzetah
!!$            IF (ALLOCATED(cputimes)) DEALLOCATE(cputimes)
!!$            IF (ALLOCATED(dmn)) DEALLOCATE(dmn)
!!$            IF (ALLOCATED(fmn)) DEALLOCATE(fmn)
!!$            IF (ALLOCATED(rfmn)) DEALLOCATE(rfmn)
!!$            IF (ALLOCATED(alpha1mn)) DEALLOCATE(alpha1mn)
!!$            IF (ALLOCATED(trigsv)) DEALLOCATE(trigsv)
!!$            IF (ALLOCATED(trigsu)) DEALLOCATE(trigsu)
!!$            allocate (cputimes(irup))
!!$            allocate (dmn(-mbuse:mbuse,0:nbuse), fmn(-mbuse:mbuse,0:nbuse),&
!!$                rfmn(-mbuse:mbuse,0:nbuse),alpha1mn(-mbuse:mbuse,0:nbuse),&
!!$                trigsv(ntrigv),trigsu(ntrigu), stat=irho)
!!$            if (irho .ne. 0) stop 'allocation error in bootsj main'
!!$            cputimes=0; dmn=0; fmn=0; rfmn=0; alpha1mn=0; trigsv=0; trigsu=0
!!$            ! Now we need to calculate indexes
!!$            ! Now loop over surfaces
!!$            ihere = 0
!!$            idx = 0
!!$            DO irho=1,irup
!!$               IF (lbooz(irho+1)) idx(irho) = 1
!!$            END DO
!!$            l_boot_all = .TRUE.
!!$            do irho=1, irup
!!$               if(idx(irho) .eq. 0) l_boot_all = .false.
!!$            enddo
!!$            call fftfax_g (nthetah, ifaxu, trigsu) ! Initialize trigsu
!!$            call cftfax_g (nzetah, ifaxv, trigsv)  ! Initialize trigsv
!!$            IF (lscreen) THEN
!!$               WRITE(6,'(A)')             '==========================================='
!!$               WRITE(6,'(A,A,A)')         '=========  B O O T S J (v.',version_(1:4),')  =========='
!!$               WRITE(6,'(2X,2(A,I3.3))')   'M_CALC = ',mbuse,';   M_BOOZER = ',mboz 
!!$               WRITE(6,'(2X,2(A,I3.3))')   'N_CALC = ',nbuse,';   N_BOOZER = ',nboz
!!$               WRITE(6,'(2X,2(A,I3.3))')   'NTHETA = ',nthetah,';   NZETA    = ',nzetah
!!$               WRITE(6,'(2X,A,F10.4)')   'ZEFF =',zeff1
!!$               WRITE(6,'(2X,A,F10.4,A)') 'TE0  =',tempe0,' [keV]'
!!$               WRITE(6,'(2X,A,F10.4,A)') 'TI0  =',tempi0,' [keV]'
!!$               WRITE(6,'(2X,A,F10.4,A)') 'NE0  =',dense(1),'x10^20 [m^-3]'
!!$               WRITE(6,'(2X,A,F10.4,A)') 'NI0  =',densi(1),'x10^20 [m^-3]'
!!$               IF (l_boot_all) WRITE(6,'(2X,A)') '<FULL CURRENT CALCULATTION>'
!!$               WRITE(6,'(A)') '-------------------------------------------'
!!$               WRITE(6,'(A)') '   dex      rho      Te[keV]     Ti[keV]       Ne           Ni        BETA       J_BOOT     TOK_FRAC'
!!$               CALL FLUSH(6)             
!!$            END IF
!!$            bsnorm =0; capr = 0; caps = 0; ftrapped =0; h2 =0;
!!$            amain =0; aiterm1 = 0; other1 =0;
!!$            cputimes = 0; dibs=0; aibs=0; gbsnorm =0;
!!$            bsdenste = 0; bsdensti = 0; bstempte = 0; bstempti = 0;
!!$            bsdense = 0; bsdensi = 0; bstempe = 0; bstempi = 0;
!!$            dibst = 0; aibst = 0; amain = 0; ajBbs = 0;
!!$            DO irho = mystart, myend
!!$               IF (idx(irho) .eq. 0) CYCLE
!!$!DEC$ IF DEFINED (MPI_OPT)
!!$               time1 = MPI_WTIME()
!!$!DEC$ ELSE
!!$               call second0 (time1)
!!$!DEC$ ENDIF
!!$               !timecpu = time2 - time1
!!$               !cputimes(irho) = timecpu
!!$               irho1 = irho -1
!!$               r     = SQRT(rhoar(irho) + 1.0E-36_dp)
!!$               CALL bongrid(irho,ians,ihere)
!!$               x     = fttok(irho)/(fptok(irho)+1.0E-36_dp)
!!$               al31t = al31(x,zeff1,alphae,alphai)
!!$               CALL grad(gradbs1,gradbs2,gradbs3,gradbs4,irho)
!!$               bsdenste(irho) = gradbs1*al31t          !due to dens gradient
!!$               bsdensti(irho) = gradbs2*al31t          !due to dens gradient
!!$               bstempte(irho) = gradbs3*al31t          !due to temp gradient
!!$               bstempti(irho) = gradbs4*al31t          !due to temp gradient
!!$               dibst(irho) = bsdenste(irho) + bsdensti(irho) + bstempte(irho) + bstempti(irho) !total Jbst
!!$               IF (l_boot_all) THEN
!!$                  IF (irho .eq. 1) THEN
!!$                     aibst(1) = dibst(1)*d_rho(1)
!!$                  ELSE
!!$                     aibst(irho) = aibst(irho1)+dibst(irho)*d_rho(irho)
!!$                  END IF
!!$               END IF
!!$               CALL denmf (trigsu, trigsv, ifaxu, ifaxv, irho)
!!$               CALL caprsh2(irho)
!!$               IF (h2(irho) == 0.0) CYCLE
!!$               CALL woflam (trigsu, trigsv, ifaxu, ifaxv, irho)
!!$               CALL othersums(irho)
!!$               amain(irho) = 0.5_rprec*(1.0_rprec - aiogar(irho)/qsafety(irho)) + 0.5_rprec*(1.0_rprec + aiogar(irho)/qsafety(irho))*h2(irho)
!!$               gbsnorm(irho) = amain(irho) + other1(irho) +  aiterm1(irho)
!!$               x = ftrapped(irho)/(fpassing(irho)+1.E-36_dp)
!!$               al31s = al31(x,zeff1,alphae,alphai)
!!$               call grad (gradbs1, gradbs2, gradbs3, gradbs4, irho)
!!$               bsdense(irho) = gbsnorm(irho)*gradbs1*al31s
!!$               bsdensi(irho) = gbsnorm(irho)*gradbs2*al31s
!!$               bstempe(irho) = gbsnorm(irho)*gradbs3*al31s
!!$               bstempi(irho) = gbsnorm(irho)*gradbs4*al31s
!!$               dibs(irho) = bsdense(irho) + bsdensi(irho) + bstempe(irho) + bstempi(irho)
!!$               ajBbs(irho) = (2.0e6_dp)*dmu0*dibs(irho)*(pres1(irho)/betar(irho))/psimax 
!!$               IF (l_boot_all) THEN
!!$                  IF (irho .eq. 1) THEN
!!$                     aibs(1) = dibs(1)*d_rho(1)
!!$                  ELSE
!!$                     aibs(irho) = aibs(irho1)+dibs(irho)*d_rho(irho)
!!$                  END IF
!!$               END IF
!!$               bsnorm(irho) = dibs(irho)/(dibst(irho)+1.E-36_dp)
!!$!DEC$ IF DEFINED (MPI_OPT)
!!$               time2 = MPI_WTIME()
!!$!DEC$ ELSE
!!$               call second0 (time2)
!!$!DEC$ ENDIF
!!$               timecpu = time2 - time1
!!$               cputimes(irho) = time2 - time1
!!$               IF (lscreen_bootsj) WRITE(6,'(2X,I3,8(2X,E11.4))') irho,rhoar(irho),tempe1(irho),tempi1(irho),dense(irho),densi(irho),betar(irho),ajBbs(irho),bsnorm(irho)
!!$               CALL FLUSH(6) 
!!$            END DO
!!$!DEC$ IF DEFINED (MPI_OPT)
!!$            CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
!!$            IF (myworkid == master) THEN
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,dibs,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,ajBbs,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,bsnorm,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,capr,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,caps,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,ftrapped,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,h2,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,amain,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,aiterm1,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,other1,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,gbsnorm,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,thetamax,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,zetahmax,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,bmax1,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,fttok,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,fptok,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,fpassing,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,cputimes,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,aibs,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,bsdense,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,bsdensi,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,bstempe,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,bstempi,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(MPI_IN_PLACE,aibst,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$            ELSE
!!$               CALL MPI_REDUCE(dibs,dibs,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(ajBbs,ajBbs,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(bsnorm,bsnorm,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(capr,capr,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(caps,caps,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(ftrapped,ftrapped,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(h2,h2,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(amain,amain,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(aiterm1,aiterm1,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(other1,other1,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(gbsnorm,gbsnorm,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(thetamax,thetamax,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(zetahmax,zetahmax,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(bmax1,bmax1,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(fttok,fttok,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(fptok,fptok,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(fpassing,fpassing,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(cputimes,cputimes,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(aibs,aibs,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(bsdense,bsdense,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(bsdensi,bsdensi,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(bstempe,bstempe,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(bstempi,bstempi,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               CALL MPI_REDUCE(aibst,aibst,irup,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
!!$               IF (ALLOCATED(cputimes)) DEALLOCATE(cputimes)
!!$               IF (ALLOCATED(bmnc_b)) DEALLOCATE(bmnc_b)
!!$               IF (ALLOCATED(bmns_b)) DEALLOCATE(bmns_b)
!!$               CLOSE(UNIT=ians,STATUS='DELETE')
!!$               CALL deallocate_all
!!$               RETURN
!!$            END IF
!!$!DEC$ ENDIF
!!$            ! Recompute AIBS
!!$            IF (l_boot_all) THEN
!!$               aibs(1) = dibs(1)*d_rho(1)
!!$               DO irho = 2,irup
!!$                  aibs(irho) = aibs(irho-1) + dibs(irho)*d_rho(irho)
!!$               END DO
!!$            END IF
!!$            IF (lscreen_bootsj .and. l_boot_all) WRITE(6,'(1X,A,E22.14,A)') 'Total bootstrap current: ',aibs(irup),' [MA]'
!!$            l_boot_all = .FALSE. ! To Prevent output from writing it to the screen again.
!!$            IF (lscreen) THEN
!!$               WRITE(6,'(A)')         '==========================================='
!!$            END IF
!!$            CALL output(cputimes,aibstot,ijbs,ians,ians_plot)
!!$            IF (ALLOCATED(cputimes)) DEALLOCATE(cputimes)
!!$            CLOSE(ians)
!!$            CLOSE(ians_plot)
!!$            CLOSE(ijbs)
         CASE('spec')
      END SELECT
      IF (lscreen) WRITE(6,'(a)') ' --------------------  SFINCS BOOTSTRAP CALCULATION DONE  --------------------'
      RETURN
  90  format(5e16.8)
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
!DEC$ ELSE
      STOP "Error! stellopt_sfincs was called, but STELLOPT was compiled without SFINCS."
!DEC$ ENDIF
    END SUBROUTINE stellopt_sfincs
