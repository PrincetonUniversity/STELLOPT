!-----------------------------------------------------------------------
!     Subroutine:    thrift_travis
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          11/24/2022
!     Description:   This subroutine calculates envokes the TRAVIS
!                    ECRH Raytracing code to calculate the current
!                    drive
!-----------------------------------------------------------------------
      SUBROUTINE thrift_travis(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_input_mod
      USE thrift_vars
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(inout) :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      INTEGER ::  i, nra, nphi
      REAL(8)        :: hgrid,dphi,B_scale,B0_ref,B_dir,phi_ref,phibx
      REAL(8), ALLOCATABLE :: te_prof(:),ne_prof(:),z_prof(:),rho_prof(:)
      CHARACTER(4)   :: antennaCoordType,targetPositionType
      CHARACTER(256) :: B0type,labelType,equiname,interpType
      INTEGER, parameter :: nprof_travis = 64

      ! Old
      integer, parameter :: nfax = 13
      INTEGER ::  i,j,n,narr,ier, iunit_rzuv
      ! FOR TRAVIS
      INTEGER(4)     :: nra,nphi,QOemul,wmode,ECEerror
      INTEGER(4)     :: travisoutTEMP,travisScrOut
      INTEGER(4)     :: maxSteps,stopray,npass,maxHarm,nrho,nu
      INTEGER(8)     :: mconf8=0,mVessel=0,mMirror=0
      REAL(8)        :: maxlength,minStepSize,maxStepSize,odetolerance,umax
      REAL(8)        :: antennaPosition(3),targetPosition(3),rbeam(2),rfocus(2)
      CHARACTER(256) :: hamilt,dieltensor
      INTEGER, ALLOCATABLE :: mnum(:)
      INTEGER       :: mystart,myend, chunk, numprocs_local,s

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  ECRH (TRAVIS) CALCULATION  -------------------------'
      IF (lvmec) THEN
#if defined(TRAVIS)

#if defined(MPI_OPT)
         ! Broadcast some information
         CALL MPI_COMM_RANK(MPI_COMM_MYWORLD, myworkid, ierr_mpi)
         CALL MPI_BCAST(nprof_travis,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(proc_string,256,MPI_CHARACTER,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(Aminor,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(Baxis,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
#endif
 
         ! Setup profiles
         ALLOCATE(te_prof(nprof_travis),ne_prof(nprof_travis),z_prof(nprof_travis),&
            rho_prof(nprof_travis))
         IF (myworkid == master) THEN
            CALL FLUSH(6)
            DO i = 1,nprof_travis
               rho_prof(i) = DBLE(i-1)/DBLE(nprof_travis-1)
               CALL get_prof_ne(rho_prof(i),mytdex,ne_prof(i))
               CALL get_prof_te(rho_prof(i),mytdex,te_prof(i))
               CALL get_prof_zeff(rho_prof(i),mytdex,z_prof(i))
            END DO
            te_prof = te_prof*1.0E-3 ! Must be in keV
         ELSE
            lscreen = .FALSE.
         END IF

#if defined(MPI_OPT)
         ! Broadcast the profiles
         CALL MPI_BCAST(rho_prof,nprof_travis,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(ne_prof,nprof_travis,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(te_prof,nprof_travis,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(z_prof,nprof_travis,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
#endif            
         ! Setup Beams
         antennaCoordType     = antennatype_ecrh(1:4)
         targetPositionType     = targettype_ecrh(1:4)
         nra = nra_ecrh
         nphi = nphi_ecrh
         nbeams = 0
         DO i = 1,nsys
            IF (ANY(antennaposition_ecrh(i,:) .ne. 0)) nbeams = nbeams + 1
         END DO

         ! Divide up work
         ! Note we should fix this in the end so we divide up over system ans well
         CALL MPI_CALC_MYRANGE(MPI_COMM_MYWORLD, 1, nbeams, mystart, myend)

            ! Everyone allocates to make the reduce simpler
            !IF (ALLOCATED(radto_ece)) DEALLOCATE(radto_ece)
            !IF (ALLOCATED(radtx_ece)) DEALLOCATE(radtx_ece)
            !ALLOCATE(radto_ece(nsys,narr))
            !ALLOCATE(radtx_ece(nsys,narr))
            !radto_ece = 0.0; radtx_ece = 0.0

         IF (mystart <= myend) THEN

            ! Handle Output
            CALL Disable_Output_f77
            CALL Disable_printout_f77

            ! Obligatory setup
            CALL closeErrorReport_f77()
            CALL set_caller_name_f77('Stellopt','ECRH')
            CALL set_hedi_f77('ECRH')

            ! Setup mirrors and vessel
            vessel_ecrh = TRIM(vessel_ecrh)
            mirror_ecrh = TRIM(mirror_ecrh)
            CALL set_VesselMirrors_f77(vessel_ecrh,mirror_ecrh)

            ! Init Equilibrium
            hgrid   = Aminor/30
            dphi    = 2
            B0type  = 'noscale'
            phi_ref = 0
            B_scale = 1
            B0_ref  = ABS(Baxis)
            B_dir   = 0
            CALL initEquilibr_f77(hgrid, dphi, B0_ref, phi_ref, B_dir, B_scale, B0type)

            ! Set Plasma size
            CALL set_plasmaSize_f77(rho_prof(nrad), rho_prof(nrad)-rho_prof(nrad-1))

            ! Load the equilibrium
            equiname = 'wout_'//TRIM(proc_string)//'.nc'
            CALL set_EquiFile_f77(equiname)

            ! Set Profiles
            labelType  = 'tor_rho'
            interpType = 'spline'
            CALL set_nTZ_profs_f77(nrad, rho, ne_prof, te_prof, z_prof, labelType, interpType)

            ! Loop over ECRH systems
            n=1
            antennaCoordType   = antennatype_ece(1:4)
            targetPositionType = targettype_ece(1:4)
            nra                = nra_ece
            nphi               = nphi_ece
            DO i = mystart,myend

               ! Set Configuration
               maxSteps     = 5000   ! upper number of steps on the trajectory
               maxlength    = 5      ! upper limit of the trajectory length [m]
               maxStepSize  = 10     ! max step-size for ODE-integrator [wave_length]
               minStepSize  = 1e-5   ! min step-size for ODE-integrator [wave_length]
               ODEtolerance = 1e-5   ! tolerance for ODE-solver
               stopRay      = 1      ! trajectory stops if optical depth > 10 (1/0)
               Npass        = 3      ! number of passes through the plasma to be done
               hamilt       = 'West' ! Hamiltonian for tracing ('West' or 'Tokman')
               dieltensor   = 'cold' ! diel. tensor model for tracing ('cold' or 'weakly_relativistic')
               maxHarm      = 0      ! upper harmonic for sum. in diel.tensor (if 0 then automatic)
               umax         = 7      ! upper velocity for integr. along the res.line in u-space
               nu           = 700    ! number of grid points in u-space
               nrho         = 100    ! internal grid for the flux-surface label
               CALL set_TracingData_f77(maxSteps,  maxlength, &
                             maxStepSize, minStepSize, odetolerance, &
                             stopray, npass, &
                             maxHarm, umax, nu, nrho, hamilt, dieltensor)

               ! Set Antenna Position
               antennaPosition(1:3) = antennaPosition_ecrh(i,1:3)
               targetPosition(1:3)  = targetPosition_ecrh(i,1:3)
               rbeam(1:2)  = rbeam_ece(i,1:2)
               QOemul      = 0
               IF (rbeam_ece(i,3) > 0) QOemul = 1
               Rfocus(1:2) = rfocus_ece(i,1:2)
               phibx       = rfocus_ece(i,3)
               CALL set_Antenna_f77( rbeam, rfocus, phibx, QOemul, nra, nphi, &
                                     antennaPosition, targetPosition, &
                                     antennaCoordType,targetPositionType)
                  
               ! Need this here
               CALL Disable_Output_f77

               ! Run BEAM
               j=1
               wmode = 1 ! X-mode
               CALL run_ECRH_Beam_f77m(n,freq_ecrh(i),wmode,power_ecrh(i))
               wmode = 0 ! O-mode
               CALL run_ECRH_Beam_f77m(n,freq_ecrh(i),wmode,power_ecrh(i))

            END DO

               ! Free the stuff we loaded.
               CALL Free_MagConfig_f77()

            END IF

            ! Deallocate variables
            DEALLOCATE(rho_prof,te_prof,ne_prof,z_prof)

            ! Print to screen
            IF (lscreen) THEN
               WRITE(6,'(5X,A,3X,A,3X,A,3X,A)') 'Beam','Freq [GHz]','Trad (0) [keV]','Trad (X) [keV]'
               DO i = 1, nsys
                  IF (ALL(sigma_ece(i,:) >= bigno)) CYCLE
                  DO j = 1, narr
                     IF (sigma_ece(i,j) >= bigno) CYCLE
                     WRITE(6,'(5X,i3,3x,f8.2,2(5x,f10.3))') i,freq_ece(i,j),radto_ece(i,j),radtx_ece(i,j)
                  END DO
               END DO
            END IF
#endif
         CASE('spec')
      END SELECT
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  ECE (TRAVIS) CALCULATION DONE  ---------------------'
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_travis
