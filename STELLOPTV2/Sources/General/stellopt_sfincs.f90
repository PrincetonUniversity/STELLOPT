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
      USE mpi_params
      USE mpi_inc

      USE sfincs_main, only: sfincs_init, sfincs_prepare, sfincs_run
      USE globalVariables, only: sfincs_inputFilename => inputFilename, sfincs_outputFilename => outputFilename, equilibriumFile, FSABjHat, FSABHat2
      USE equil_vals, only: phiedge
      USE, intrinsic :: iso_fortran_env, only : stdout=>output_unit, stderr=>error_unit
!-----------------------------------------------------------------------
!     Subroutine Parameters
!----------------------------------------------------------------------
!DEC$ ENDIF
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!DEC$ IF DEFINED (SFINCS)
!-----------------------------------------------------------------------
!     Local Variables
!----------------------------------------------------------------------
      integer, parameter :: nfax = 13
      INTEGER ::  ier, i, j, ik, ntheta, nzeta, nlis
      INTEGER ::  dex_ion, dex_zeff
      INTEGER :: mystart,myend, chunk, numprocs_local
      INTEGER, ALLOCATABLE :: mnum(:)
      REAL(rprec) :: t1, t2
      integer, dimension(nfax) :: ifaxu, ifaxv
      integer :: ntrigu, ntrigv, i1
      integer :: ntheta_min, nzeta_min, ir, mn, m, n
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
      integer :: numProcs_world, numProcs_myworld, myRank_world, myRank_myworld, myRank_sfincs, numProcs_sfincs
      integer :: Nradii, radius_index_min, radius_index_max, radius_index, color, key
      integer :: MPI_COMM_SFINCS
      integer, parameter :: buffer_length = 1000
      character(len=buffer_length) :: proc_assignments_string, base_directory_string, directory_string, file_line, file_line_lower, file_line_lower2
      character(len=buffer_length) :: working_directory
      integer :: tag, file_status, unit_in, unit_out
      integer :: my_mpi_status(MPI_STATUS_SIZE)
      LOGICAL :: lfile_check
      LOGICAL :: added_scanType
      REAL(rprec) :: sfincs_ne, sfincs_ni, sfincs_Te, sfincs_Ti, sfincs_d_ne_d_s, sfincs_d_ni_d_s, sfincs_d_Te_d_s, sfincs_d_Ti_d_s, delta_s
      INTEGER :: numProcs_eff, myRank_eff
      REAL(rprec) :: temp1, temp2, sfincs_dPhiHatdpsiN

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (lscreen) WRITE(6,'(a)') ' ---------------------  BOOTSTRAP CALCULATION USING SFINCS  ------------------'
      SELECT CASE(TRIM(equil_type))
         CASE('vmec2000','animec','flow','satire','parvmec','paravmec','vboot','vmec2000_oneeq')

            base_directory_string = 'sfincs_'//TRIM(proc_string)
            IF (myworkid == master) CALL execute_command_line('mkdir -p ' //TRIM(base_directory_string)) ! The -p blocks a warning message if the directory already exists.
            CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs_world, ierr_mpi)
            CALL MPI_COMM_RANK(MPI_COMM_WORLD, myRank_world, ierr_mpi)
            CALL MPI_COMM_SIZE(MPI_COMM_MYWORLD, numProcs_myworld, ierr_mpi)
            CALL MPI_COMM_RANK(MPI_COMM_MYWORLD, myRank_myworld, ierr_mpi)

            IF (sfincs_min_procs > numProcs_myworld) THEN
               IF (myworkid==master) THEN
                  PRINT *,"WARNING! The number of procs in MPI_COMM_MYWORLD is smaller than sfincs_min_procs."
                  PRINT *,"# of procs in MPI_COMM_MYWORLD:",numProcs_myWorld
                  PRINT *,"sfincs_min_procs:",sfincs_min_procs
                  PRINT *,"Lowering sfincs_min_procs to",numProcs_myWorld
               END IF
               sfincs_min_procs = numProcs_myWorld
            END IF
            Nradii = MINLOC(sfincs_s(2:),DIM=1)
            IF (myworkid==master .and. lscreen) PRINT *,"Number of radii for SFINCS:",Nradii
            DO j=1,Nradii
               IF (sfincs_s(j) <= 0) STOP "Error! sfincs_s values must be > 0 (except for 0 values at the end of the array, which are ignored)."
               IF (sfincs_s(j) > 1) STOP "Error! sfincs_s values must be <= 1."
            END DO

            ! Divide up MPI_COMM_MYWORLD over radii:
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
            CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, color, key, MPI_COMM_SFINCS, ierr_mpi)
            CALL MPI_COMM_SIZE(MPI_COMM_SFINCS, numProcs_sfincs, ierr_mpi)
            CALL MPI_COMM_RANK(MPI_COMM_SFINCS, myRank_sfincs, ierr_mpi)
            WRITE (proc_assignments_string,fmt="(a,i5,a,i5,a,i4,a,i4,a,i5,a,i5,a)"),"Proc ",myRank_myWorld," of",numProcs_myWorld," in MPI_COMM_MYWORLD will handle radius",radius_index_min," to",radius_index_max," and has rank",myRank_sfincs," of",numProcs_sfincs," in MPI_COMM_SFINCS."

            ! Print the processor/radius assignments in a coordinated manner.
            IF (myworkid==master) THEN
               IF (lscreen) WRITE(*,"(a)") TRIM(proc_assignments_string)
               DO i = 1,numProcs_myWorld - 1
                  tag = i
                  CALL MPI_RECV(proc_assignments_string,buffer_length,MPI_CHAR,i,tag,MPI_COMM_MYWORLD,my_mpi_status,ierr_mpi)
                  IF (lscreen) WRITE(*,"(a)") TRIM(proc_assignments_string)
               END DO
            ELSE
               tag = myworkid
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

               ! Handle radial electric field
               CALL tolower(sfincs_Er_option)
               select case (trim(sfincs_Er_option))
               case ("zero")
                  sfincs_dPhiHatdpsiN = 0.001 ! Sfincs sometimes has roundoff error problems when Er is exactly 0, so use a tiny nonzero value instead.
               case ("estimate")
                  sfincs_dPhiHatdpsiN = -((sfincs_Ti / sfincs_ni) * sfincs_d_ni_d_s + 1.34 * sfincs_d_Ti_d_s)
                  ! The number 1.34 comes from eq (74a) in Ho & Kulsrud: 1.34 = (4.63/1.63) - (3/2).
               case default
                  print *,"Error! Unrecognized sfincs_Er_option:",sfincs_Er_option
                  print *,"Allowed settings are 'zero' or 'estimate'."
                  stop
               end select

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
                        WRITE(UNIT=unit_out,FMT='(A)') '  geometryScheme = 5'
                        WRITE(UNIT=unit_out,FMT='(A)') '  VMECRadialOption = 0'
                        WRITE(UNIT=unit_out,FMT='(A)') '  inputRadialCoordinate = 1'
                        WRITE(UNIT=unit_out,FMT='(A)') '  inputRadialCoordinateForGradients = 1'
                        WRITE(UNIT=unit_out,FMT='(a, es24.14)') '  psiN_wish = ',sfincs_s(radius_index)
                        WRITE(UNIT=unit_out,FMT='(A)') '  equilibriumFile = "'//TRIM(working_directory)//'/wout_'//TRIM(proc_string)//'.nc"'
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

                     IF (TRIM(file_line_lower)=='&physicsparameters') THEN
                        WRITE (UNIT=unit_out,FMT='(A)') '&physicsParameters'
                        WRITE (UNIT=unit_out,FMT='(a, es24.14, a)') '  dPhiHatdpsiN = ', sfincs_dPhiHatdpsiN,' ! From stellopt'
                        CYCLE
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
                     IF (file_line_lower(1: 7)=='rhsmode') CYCLE
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
                     IF (file_line_lower(1:12)=='dphihatdpsin') CYCLE
                     IF (file_line_lower(1:14)=='dphihatdpsihat') CYCLE

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

               sfincs_inputFilename  = TRIM(directory_string)//"/input.namelist"
               sfincs_outputFilename = TRIM(directory_string)//"/sfincsOutput.h5"
               ! Call SFINCS:
               CALL sfincs_init(MPI_COMM_SFINCS)
               CALL sfincs_prepare()
               CALL sfincs_run()
               CALL sfincs_finalize()
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
               ierr=0
               CALL safe_open(unit_out,ierr,'sfincs_results_before_profile_fitting.'//TRIM(proc_string),'REPLACE','formatted')
               WRITE (unit_out,"(a)") "s (normalized toroidal flux), <j dot B> (Tesla Amperes / meters^2), <B^2> (Tesla^2)"
               WRITE (unit_out,*) 0.0d+0,0.0d+0,sfincs_B_squared_flux_surface_average(1)
               DO radius_index = 1, Nradii
                  WRITE (unit_out,*) sfincs_s(radius_index), sfincs_J_dot_B_flux_surface_average(radius_index+1), sfincs_B_squared_flux_surface_average(radius_index+1)
               END DO
               CLOSE (unit_out)
            ELSE
               CALL MPI_REDUCE(sfincs_J_dot_B_flux_surface_average,sfincs_J_dot_B_flux_surface_average,Nradii+1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(sfincs_B_squared_flux_surface_average,sfincs_B_squared_flux_surface_average,Nradii+1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            END IF

            CALL MPI_COMM_FREE(MPI_COMM_SFINCS,ierr_mpi)
            CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi) ! Weird things might happen if some procs move on while others are still running sfincs.

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
