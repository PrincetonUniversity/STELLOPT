!-----------------------------------------------------------------------
!     Subroutine:    stellopt_travis
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/28/2016
!     Description:   This subroutine calculates envokes the TRAVIS
!                    ECE Raytracing code to calculated radiated power.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_travis(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE equil_utils
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
      integer, parameter :: nfax = 13
      INTEGER ::  i,j,n,narr,ier, iunit_rzuv
      ! FOR TRAVIS
      INTEGER(4)     :: nra,nphi,QOemul,wmode,ECEerror
      INTEGER(4)     :: travisoutTEMP,travisScrOut
      INTEGER(4)     :: maxSteps,stopray,npass,maxHarm,nrho,nu
      INTEGER(8)     :: mconf8=0,mVessel=0,mMirror=0
      REAL(8)        :: hgrid,dphi,B_scale,B0_ref,B_dir,phi_ref,phibx
      REAL(8)        :: maxlength,minStepSize,maxStepSize,odetolerance,umax
      REAL(8)        :: antennaPosition(3),targetPosition(3),rbeam(2),rfocus(2)
      REAL(8), ALLOCATABLE :: te_prof(:),ne_prof(:),z_prof(:)
      CHARACTER(4)   :: antennaCoordType,targetPositionType
      CHARACTER(256) :: B0type,labelType,equiname,interpType
      CHARACTER(256) :: hamilt,dieltensor
      INTEGER, ALLOCATABLE :: mnum(:)
      INTEGER       :: mystart,myend, chunk, numprocs_local,s

      INTERFACE
         INTEGER(4) FUNCTION mcload(mconf,name) 
            INTEGER(8)   :: mconf
            CHARACTER(*) :: name
         END FUNCTION mcLoad
         SUBROUTINE mcFree(mconf) 
            INTEGER(8)   :: mconf
         END SUBROUTINE mcFree
         SUBROUTINE VesselFree(mconf) 
            INTEGER(8)   :: mconf
         END SUBROUTINE VesselFree
      END INTERFACE
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (lscreen) WRITE(6,'(a)') ' ----------------------------  ECE (TRAVIS) CALCULATION  -------------------------'
      SELECT CASE(TRIM(equil_type))
         CASE('vmec2000','animec','flow','satire','parvmec','paravmec','vboot','vmec2000_oneeq')
#if defined(TRAVIS)

            ! Count work to be done
            narr = 0
            DO i = 1, nsys
               DO j = 1, nprof
                  IF (sigma_ece(i,j)< bigno) narr = j
               END DO
            END DO

            ! Broadcast some information
#if defined(MPI_OPT)
            CALL MPI_COMM_RANK(MPI_COMM_MYWORLD, myworkid, ierr_mpi)
            CALL MPI_BCAST(nrad,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(Aminor,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(Baxis,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
#endif
 
            ! Setup profiles
            ALLOCATE(te_prof(nrad),ne_prof(nrad),z_prof(nrad))
            IF (myworkid == master) THEN
               CALL FLUSH(6)
               DO i = 1,nrad
                  CALL get_equil_ne(shat(i),TRIM(ne_type),ne_prof(i),ier)
                  CALL get_equil_te(shat(i),TRIM(te_type),te_prof(i),ier)
                  CALL get_equil_zeff(shat(i),TRIM(zeff_type),z_prof(i),ier)
               END DO
               te_prof = te_prof*1.0E-3 ! Must be in keV
            ELSE
               lscreen = .FALSE.
               ALLOCATE(rho(nrad))
            END IF

            ! Broadcast the profiles
#if defined(MPI_OPT)
            CALL MPI_BCAST(rho,nrad,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(ne_prof,nrad,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(te_prof,nrad,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(z_prof,nrad,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
#endif

            ! Divide up work
            ! Note we should fix this in the end so we divide up over system ans well
            CALL MPI_CALC_MYRANGE(MPI_COMM_MYWORLD, 1, narr, mystart, myend)

            ! Everyone allocates to make the reduce simpler
            IF (ALLOCATED(radto_ece)) DEALLOCATE(radto_ece)
            IF (ALLOCATED(radtx_ece)) DEALLOCATE(radtx_ece)
            ALLOCATE(radto_ece(nsys,narr))
            ALLOCATE(radtx_ece(nsys,narr))
            radto_ece = 0.0; radtx_ece = 0.0

            IF (mystart <= myend) THEN

               ! Handle Output
               CALL Disable_Output_f77
               CALL Disable_printout_f77

               ! Obligatory setup
               CALL set_hedi_f77('ece')

               ! Setup mirrors and vessel
               vessel_ece = TRIM(vessel_ece)
               mirror_ece = TRIM(mirror_ece)
               CALL set_VesselMirrors_f77(vessel_ece,mirror_ece)

               ! Set Profiles
               labelType  = 'tor_rho'
               interpType = 'spline'
               CALL set_nTZ_profs_f77(nrad, rho, ne_prof, te_prof, z_prof, labelType, interpType)

               ! Init Equilibrium
               hgrid   = Aminor/30
               dphi    = 2
               B0type  = 'at angle on magn.axis'
               B0type  = 'noscale'
               phi_ref = 0
               B_scale = 1
               B0_ref  = Baxis
               B_dir   = SIGN(1.0D0, Baxis)
               CALL initEquilibr_f77(hgrid, dphi, B0_ref, phi_ref, B_dir, B_scale, B0type)

               ! Load the equilibrium
               equiname = 'wout_'//TRIM(proc_string)//'.nc'
               CALL set_EquiFile_f77(equiname)

               ! Loop over ECE systems
               n=1
               antennaCoordType   = antennatype_ece(1:4)
               targetPositionType = targettype_ece(1:4)
               nra                = nra_ece
               nphi               = nphi_ece
               DO i = 1,nsys

                  ! Don't evaluate if not set or not in [mystart,myend]
                  IF (ALL(sigma_ece(i,:) >= bigno)) CYCLE
                  IF (sigma_ece(i,mystart) >= bigno) CYCLE

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
                  antennaPosition(1:3) = antennaPosition_ece(i,1:3)
                  targetPosition(1:3)  = targetPosition_ece(i,1:3)
                  rbeam(1:2)  = rbeam_ece(i,1:2)
                  QOemul      = 0
                  IF (rbeam_ece(i,3) > 0) QOemul = 1
                  Rfocus(1:2) = rfocus_ece(i,1:2)
                  phibx       = rfocus_ece(i,3)
                  CALL set_Antenna_f77( rbeam, rfocus, phibx, QOemul, nra, nphi, &
                                        antennaPosition, targetPosition, &
                                        antennaCoordType,targetPositionType)

                  ! Cycle over frequencies
                  DO j = mystart,myend
                     IF (sigma_ece(i,j) >= bigno) CYCLE
                     wmode = 0 ! O-mode
                     CALL run_ECE_Beam_f77m(n,freq_ece(i,j), wmode, radto_ece(i,j),ECEerror)
                     wmode = 1 ! X-mode
                     CALL run_ECE_Beam_f77m(n,freq_ece(i,j), wmode, radtx_ece(i,j),ECEerror)
                     PRINT *,myworkid,i,j,ECEerror
                  END DO

               ! Free the vessel here (only loaded if we get to here)
               !CALL VesselFree(mVessel)
               END DO
            END IF

#if defined(MPI_OPT)
            IF (myworkid == master) THEN
               CALL MPI_REDUCE(MPI_IN_PLACE,radto_ece,nsys*narr,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,radtx_ece,nsys*narr,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            ELSE
               CALL MPI_REDUCE(radto_ece,radto_ece,nsys*narr,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(radtx_ece,radtx_ece,nsys*narr,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            END IF
#endif

            ! Deallocate variables
            !CALL MCFREE(mconf8)
            DEALLOCATE(te_prof,ne_prof,z_prof)
            IF (myworkid /= master) DEALLOCATE(radto_ece,radtx_ece,rho)

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
      END SUBROUTINE stellopt_travis
