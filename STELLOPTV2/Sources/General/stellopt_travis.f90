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
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'
!DEC$ ENDIF
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      integer, parameter :: nfax = 13
      INTEGER ::  i,j,n,narr,ier, iunit_rzuv
      ! FOR TRAVIS
      INTEGER(4)     :: nra,nphi,QOemul,wmode
      INTEGER(4)     :: travisoutTEMP,travisScrOut
      INTEGER(4)     :: maxSteps,stopray,npass,maxHarm,nrho,nu
      INTEGER(8)     :: mconf8=0,mVessel=0
      REAL(8)        :: hgrid,dphi,B_scale,B0_ref,phi_ref,phibx
      REAL(8)        :: maxlength,minStepSize,maxStepSize,odetolerance,umax
      REAL(8)        :: antennaPosition(3),targetPosition(3),rbeam(2),rfocus(2)
      REAL(8), ALLOCATABLE :: te_prof(:),ne_prof(:),z_prof(:)
      CHARACTER(4)   :: antennaCoordType,targetPositionType
      CHARACTER(256) :: B0type,labelType,equiname
      CHARACTER(256) :: hamilt,dieltensor
      INTEGER, ALLOCATABLE :: mnum(:)
      INTEGER       :: mystart,myend, chunk, numprocs_local

      INTERFACE
         INTEGER(4) FUNCTION mcload(mconf,name) 
            INTEGER(8)   :: mconf
            CHARACTER(*) :: name
         END FUNCTION mcLoad
      END INTERFACE
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (iflag < 0) RETURN
      IF (lscreen) WRITE(6,'(a)') ' ----------------------------  ECE (TRAVIS) CALCULATION  -------------------------'
      SELECT CASE(TRIM(equil_type))
         CASE('vmec2000','animec','flow','satire','parvmec','paravmec','vboot','vmec2000_oneeq')
!DEC$ IF DEFINED (TRAVIS)

            ! Load equilibrium
            equiname = 'wout_'//TRIM(proc_string)//'.nc'
            ier = mcLoad(mconf8,equiname)
            IF (ier /=1) THEN
               iflag = -1
               RETURN
            END IF
 
            ! Setup profiles
            ALLOCATE(te_prof(nrad),ne_prof(nrad),z_prof(nrad))
            IF (myworkid == master) THEN
               CALL FLUSH(6)
               DO i = 1,nrad
                  CALL get_equil_ne(rho(i),TRIM(ne_type),ne_prof(i),ier)
                  CALL get_equil_te(rho(i),TRIM(te_type),te_prof(i),ier)
                  CALL get_equil_zeff(rho(i),TRIM(zeff_type),z_prof(i),ier)
               END DO
               te_prof = te_prof*1.0E-3 ! Must be in keV
            END IF

!DEC$ IF DEFINED (MPI_OPT)
            CALL MPI_BCAST(rho,nrad,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(ne_prof,nrad,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(te_prof,nrad,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(z_prof,nrad,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(Aminor,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
            CALL MPI_BCAST(Baxis,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_MYWORLD,ierr_mpi)
!DEC$ ENDIF


            ! Handle Output
            !travisOutTemp = 0; travisScrOut=0;
            !IF (lscreen) travisScrOut = 0;
            !CALL set_TracingOutput_f77(travisOutTemp,travisScrOut)
            CALL Disable_Output_f77

            ! Handle Vessel and Mirrors
            vessel_ece = TRIM(vessel_ece)
            mirror_ece = TRIM(mirror_ece)
            CALL set_VesselMirrors_f77(vessel_ece,mirror_ece)

            ! Allocate Arrays
            narr = SUM(COUNT(sigma_ece < bigno,DIM=2))
            IF (ALLOCATED(radto_ece)) DEALLOCATE(radto_ece)
            IF (ALLOCATED(radtx_ece)) DEALLOCATE(radtx_ece)
            ALLOCATE(radto_ece(nsys,nprof))
            ALLOCATE(radtx_ece(nsys,nprof))
            radto_ece = 0.0; radtx_ece = 0.0
 
!DEC$ IF DEFINED (MPI_OPT)
            ! Divide up work
            CALL MPI_COMM_SIZE( MPI_COMM_MYWORLD, numprocs_local, ierr_mpi )
            IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
            ALLOCATE(mnum(numprocs_local))
            mnum=0
            i = 1
            DO
               IF (SUM(mnum,DIM=1) == nprof) EXIT  ! Have to use ns_b because of logic
               IF (i > numprocs_local) i = 1
               mnum(i) = mnum(i) + 1
               i=i+1
            END DO
            mystart = 1
            DO i = 1, myworkid
               mystart = SUM(mnum(1:i))+1
            END DO
            myend = mystart + mnum(myworkid+1) - 1
            IF (myend < mystart) myend = mystart
            IF (mnum(myworkid+1) == 0) mystart = myend + 1
            DEALLOCATE(mnum)
!DEC$ ENDIF

            ! Loop over ECE systems
            n=1
            antennaCoordType     = antennatype_ece(1:4)
            targetPositionType     = targettype_ece(1:4)
            nra = nra_ece
            nphi = nphi_ece
            DO i = 1,nsys
               IF (ALL(sigma_ece(i,:) >= bigno)) CYCLE
               ! Set Antenna Position
               antennaPosition(1:3) = antennaPosition_ece(i,1:3)
               targetPosition(1:3) = targetPosition_ece(i,1:3)
               rbeam(1:2) = rbeam_ece(i,1:2)
               QOemul = 0
               IF (rbeam_ece(i,3) > 0) QOemul = 1
               Rfocus(1:2) = rfocus_ece(i,1:2)
               phibx       = rfocus_ece(i,3)
               ! TEST BEAM!!!!!!!!!!!!
               !antennaPosition(1:3) = (/9.0, 0.0, 0.01/)
               !targetPosition(1:3)  = (/0.5,  0.0, 0.05/)
               !antennaCoordType     = 'cyl'
               !targetPositionType   = 'cyl'
               !Rfocus(1:2) = 10   ! beam focal lengths 
               !QOemul = 0    ! whether to emulate quasi-optical definition
               ! TEST BEAM!!!!!!!!!!!
               CALL set_Antenna_f77( rbeam,rfocus,phibx,QOemul, nra, nphi, &
                                     antennaPosition, targetPosition, &
                                     antennaCoordType,targetPositionType)

               ! Init Equilibrium
               hgrid   = Aminor/30; dphi = 2; B0type  = 'at angle on magn.axis'; phi_ref = 0; B_scale = 1; B0_ref = Baxis
               CALL initEquilibr_f77(hgrid, dphi, B0_ref, phi_ref, B_scale, B0type)

               ! Set Configuration
               CALL SET_MAGCONFIG_F77(mConf8,mVessel)

               ! Set Profiles
               labelType = 'tor_rho'
               CALL set_nTZ_profs_f77(nrad, sqrt(rho), ne_prof, te_prof, z_prof, labelType)

               ! Set Configuration
               !maxSteps     = 5000   ! upper number of steps on the trajectory
               !maxlength    = 5      ! upper limit of the trajectory length [m]
               !maxStepSize  = 10     ! max step-size for ODE-integrator [wave_length]
               !minStepSize  = 1e-5   ! min step-size for ODE-integrator [wave_length]
               !ODEtolerance = 1e-5   ! tolerance for ODE-solver
               !stopRay      = 1      ! trajectory stops if optical depth > 10 (1/0)
               !Npass        = 3      ! number of passes through the plasma to be done
               !hamilt       = 'West' ! Hamiltonian for tracing ('West' or 'Tokman')
               !dieltensor   = 'cold' ! diel. tensor model for tracing ('cold' or 'weakly_relativistic')
               !maxHarm      = 0      ! upper harmonic for sum. in diel.tensor (if 0 then automatic)
               !umax         = 7      ! upper velocity for integr. along the res.line in u-space
               !nu           = 700    ! number of grid points in u-space
               !nrho         = 100    ! internal grid for the flux-surface label
               !CALL set_TracingData_f77(maxSteps,  maxlength, maxStepSize, minStepSize, &
               !           odetolerance, stopray, npass, hamilt, dieltensor, &
               !           maxHarm, umax, nu, nrho)

               ! Print HEADER
               IF (lscreen .and. i==1) WRITE(6,'(5X,A,3X,A,3X,A,3X,A)') 'Beam','Freq [GHz]','Trad (0) [keV]','Trad (X) [keV]'

               ! Cycle over frequencies
               DO j = mystart,myend
               !DO j = 1, nprof
                  IF (sigma_ece(i,j) >= bigno) CYCLE
                  wmode = 0 ! O-mode
                  CALL run_ECE_Beam_f77m(n,freq_ece(i,j), wmode, radto_ece(i,j))
                  wmode = 1 ! X-mode
                  CALL run_ECE_Beam_f77m(n,freq_ece(i,j), wmode, radtx_ece(i,j))
                  !IF (lscreen) WRITE(6,'(5X,i3,3x,f8.2,2(5x,f10.3))') n,freq_ece(i,j),radto_ece(n),radtx_ece(n)
               END DO

            END DO

!DEC$ IF DEFINED (MPI_OPT)
            IF (myworkid == master) THEN
               CALL MPI_REDUCE(MPI_IN_PLACE,radto_ece,nsys*nprof,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(MPI_IN_PLACE,radti_ece,nsys*nprof,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            ELSE
               CALL MPI_REDUCE(radto_ece,radto_ece,nsys*nprof,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
               CALL MPI_REDUCE(radti_ece,radti_ece,nsys*nprof,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_MYWORLD,ierr_mpi)
            END IF
!DEC$ ENDIF

            IF (lscreen) THEN
              DO i = 1, nsys
                 DO j = 1, nprof
                    WRITE(6,'(5X,i3,3x,f8.2,2(5x,f10.3))') n,freq_ece(i,j),radto_ece(i,j),radtx_ece(i,j)
                 END DO
              END DO
            END IF

            DEALLOCATE(te_prof,ne_prof,z_prof)
!DEC$ ENDIF
         CASE('spec')
      END SELECT
      IF (lscreen) WRITE(6,'(a)') ' ---------------------------  ECE (TRAVIS) CALCULATION DONE  ---------------------'
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_travis
