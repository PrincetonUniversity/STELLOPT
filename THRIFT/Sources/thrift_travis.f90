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
      USE thrift_equil
      USE thrift_profiles_mod
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
      INTEGER ::  i, n
      REAL(8)        :: hgrid,dphi,B_scale,B0_ref,B_dir,phi_ref,phibx,timenow
      REAL(8), ALLOCATABLE :: te_prof(:),ne_prof(:),z_prof(:),rho_prof(:)
      CHARACTER(4)   :: antennaCoordType,targetPositionType
      CHARACTER(256) :: B0type,labelType,equiname,interpType
      INTEGER, parameter :: nprof_travis = 64

      ! For ECCD
      REAL(8) :: rho, dPdV,Pabs,jbb,jcdt,Icdt
      REAL(8), ALLOCATABLE :: Itotal(:)

      ! Old
      !integer, parameter :: nfax = 13
      !INTEGER ::  i,j,n,narr,ier, iunit_rzuv
      ! FOR TRAVIS
      !INTEGER(4)     :: nra,nphi,QOemul,wmode,ECEerror
      !INTEGER(4)     :: travisoutTEMP,travisScrOut
      INTEGER(4)     :: nra, nphi, QOemul, nbeams, wmode
      INTEGER(4)     :: maxSteps,stopray,npass,maxHarm,nrho_travis,nu
      !INTEGER(8)     :: mconf8=0,mVessel=0,mMirror=0
      REAL(8)        :: maxlength,minStepSize,maxStepSize,odetolerance,umax
      REAL(8)        :: antennaPosition(3),targetPosition(3),rbeam(2),rfocus(2)
      CHARACTER(256) :: hamilt,dieltensor
      !INTEGER, ALLOCATABLE :: mnum(:)
      INTEGER       :: mystart,myend, chunk, numprocs_local,s

!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  ECCD (TRAVIS) CALCULATION  -------------------------'
      IF (lvmec) THEN
#if defined(TRAVIS)

#if defined(MPI_OPT)
         ! Broadcast some information
         CALL MPI_COMM_RANK(MPI_COMM_MYWORLD, myworkid, ierr_mpi)
         !CALL MPI_BCAST(nprof_travis,1,MPI_INTEGER,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(proc_string,256,MPI_CHARACTER,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(eq_Aminor,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
         CALL MPI_BCAST(bt0,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
#endif
 
         ! Setup profiles
         ALLOCATE(te_prof(nprof_travis),ne_prof(nprof_travis),z_prof(nprof_travis),&
            rho_prof(nprof_travis))
         IF (myworkid == master) THEN
            CALL FLUSH(6)
            timenow = THRIFT_T(mytimestep)
            DO i = 1,nprof_travis
               rho_prof(i) = DBLE(i-1)/DBLE(nprof_travis-1)
               CALL get_prof_ne(rho_prof(i),timenow,ne_prof(i))
               CALL get_prof_te(rho_prof(i),timenow,te_prof(i))
               CALL get_prof_zeff(rho_prof(i),timenow,z_prof(i))
            END DO
            te_prof = te_prof*1.0E-3 ! Must be in keV
         ELSE
            lscreen = .FALSE.
         END IF

#if defined(MPI_OPT)
         ! Broadcast the profiles
         CALL MPI_BARRIER(MPI_COMM_MYWORLD,ierr_mpi)
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
            hgrid   = eq_Aminor/30
            dphi    = 2
            B0type  = 'noscale'
            phi_ref = 0
            B_scale = 1
            B0_ref  = ABS(bt0)
            B_dir   = 0
            CALL initEquilibr_f77(hgrid, dphi, B0_ref, phi_ref, B_dir, B_scale, B0type)

            ! Set Plasma size
            !=======================================================================
            !SUBROUTINE set_plasmaSize_f77(rho_pl,drho_bnd)
            ! Sets plasma-size, i.e. extra-surface beyond LCMS, rho_pl >= 1,
            ! and width of boundary layer for forced decay of plasma profiles:
            ! if drho_bnd = 0, no changes in plasma profiles.
            !=======================================================================
            CALL set_plasmaSize_f77(1.0, 0)

            ! Load the equilibrium
            equiname = 'wout_'//TRIM(proc_string)//'.nc'
            CALL set_EquiFile_f77(equiname)

            ! Set Profiles
            !=======================================================================
            ! SUBROUTINE set_nTZ_profs_f77(ns,rho,ne,Te,Zeff,labelType,InterpolType)
            ! Sets plasma profiles.
            ! Inputs:
            ! ns        - dimension of the arrays
            ! rho (ns)  - normalized effective radii (ascending ordered, not necessary uniform)
            ! ne  (ns)  - density, 1/m**3
            ! Te  (ns)  - temperature, keV
            ! Zeff(ns)  - effective charge
            ! labelType    = 'tor_rho' or 'pol_rho' i.e. with tor. or pol. flux 
            ! InterpolType = 'lin' or 'quad' or 'spline': kind of interpol. to be applied 
            !=======================================================================
            labelType  = 'tor_rho'
            interpType = 'spline'
            CALL set_nTZ_profs_f77(nprof_travis, rho_prof, ne_prof, te_prof, z_prof, labelType, interpType)

            ! Loop over ECRH systems
            !antennaCoordType   = antennatype_ece(1:4)
            !targetPositionType = targettype_ece(1:4)
            !nra                = nra_ece
            !nphi               = nphi_ece
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
               nrho_travis         = 100    ! internal grid for the flux-surface label
               CALL set_TracingData_f77(maxSteps,  maxlength, &
                             maxStepSize, minStepSize, odetolerance, &
                             stopray, npass, &
                             maxHarm, umax, nu, nrho_travis, hamilt, dieltensor)

               ! Set Antenna Position
               antennaPosition(1:3) = antennaPosition_ecrh(i,1:3)
               targetPosition(1:3)  = targetPosition_ecrh(i,1:3)
               rbeam(1:2)  = rbeam_ecrh(i,1:2)
               QOemul      = 0
               IF (rbeam_ecrh(i,3) > 0) QOemul = 1
               Rfocus(1:2) = rfocus_ecrh(i,1:2)
               phibx       = rfocus_ecrh(i,3)
               CALL set_Antenna_f77( rbeam, rfocus, phibx, QOemul, nra, nphi, &
                                     antennaPosition, targetPosition, &
                                     antennaCoordType,targetPositionType)
                  
               ! Need this here
               CALL Disable_Output_f77

               ! Run BEAM
               CALL run_ECRH_Beam_f77m(i,freq_ecrh(i)*1E-9,wmode_ecrh(i),power_ecrh(i)*1E-6)

            END DO

            ! Free the stuff we loaded.
            ! CALL Free_MagConfig_f77()

         END IF

         ! Allocatel total current [A]
         ALLOCATE(Itotal(nrho))
         Itotal = 0

         DO i = 1, nrho
            rho = THRIFT_RHO(i)
            CALL get_ECRH_deposition_f77(rho, dPdV,Pabs,jbb,jcdt,Itotal(i))
         END DO

#if defined(MPI_OPT)
         IF (myworkid == master) THEN
            CALL MPI_REDUCE(MPI_IN_PLACE,Itotal,nrho,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
         ELSE
            CALL MPI_REDUCE(Itotal,Itotal,nrho,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
         END IF
#endif

            ! Deallocate variables
            DEALLOCATE(rho_prof,te_prof,ne_prof,z_prof)

            ! Print to screen
            IF (lscreen) THEN
               WRITE(6,*) Itotal
!               WRITE(6,'(5X,A,3X,A,3X,A)') 'Beam','Freq [GHz]','I [kA]'
!               DO i = 1, nsys
!                  IF (ALL(sigma_ece(i,:) >= bigno)) CYCLE
!                  DO j = 1, narr
!                     IF (sigma_ece(i,j) >= bigno) CYCLE
!                     WRITE(6,'(5X,i3,3x,f8.2,1(5x,f10.3))') i,freq_ece(i,j),Itotal(i)
!                  END DO
!               END DO
            END IF
#endif
      END IF

      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  ECCD (TRAVIS) CALCULATION DONE  ---------------------'
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_travis