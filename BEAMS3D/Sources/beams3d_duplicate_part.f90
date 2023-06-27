!-----------------------------------------------------------------------
!     Module:        beams3d_duplicate_part
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          12/16/2021
!     Description:   This subroutine replicates gyrocenters.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_duplicate_part
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_runtime
      USE beams3d_lines
      USE beams3d_physics_mod, ONLY: beams3d_gc2fo, beams3d_part2gc
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          npoinc_extract Which save state to extract from file.
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL :: lfullorbit_run
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: ltemp
      INTEGER :: nparticles_new, i, j, k, mystart 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: itemp
      DOUBLE PRECISION :: time0
      DOUBLE PRECISION, DIMENSION(6) :: q
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: rtemp
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Assume we have all the IC quantiteis loaded and we want more
!      lfullorbit_run = (ANY(VR_start .ne. 0) .or. &
!                        ANY(VPHI_start .ne. 0) .or. &
!                        ANY(VZ_start .ne. 0))
      nparticles_new = nparticles*duplicate_factor
      IF (lverb) THEN
         WRITE(6,'(A)')  '----- Multiplying Particles -----'
         WRITE(6,'(A,I8)') '   DUPLICATE_FACTOR: ', duplicate_factor
         WRITE(6,'(A,I8)') '   NPARTICLES:       ', nparticles
         WRITE(6,'(A,I8)') '   NPARTICLES_NEW:   ', nparticles_new
!         IF (lfullorbit_run) THEN
!            WRITE(6,'(A)') '   Multiplying Full Orbit Run'
!         ELSE
!            WRITE(6,'(A)') '   Multiplying Gyrocenter Run'
!         END IF
      END IF

      ! Copy to new grids
      IF (myworkid == master) THEN
         ALLOCATE(ltemp(nparticles_new),itemp(nparticles_new),rtemp(nparticles_new,13))

         ! Make Duplicates
         DO i = 1, nparticles
            j = (i-1)*duplicate_factor+1
            k = j+duplicate_factor-1
            rtemp(j:k,1) = R_start(i)
            rtemp(j:k,2) = PHI_start(i)
            rtemp(j:k,3) = Z_start(i)
            rtemp(j:k,4) = VR_start(i)
            rtemp(j:k,5) = VPHI_start(i)
            rtemp(j:k,6) = VZ_start(i)
            rtemp(j:k,7) = mass(i)
            rtemp(j:k,8) = charge(i)
            rtemp(j:k,9) = mu_start(i)
            rtemp(j:k,10) = Zatom(i)
            rtemp(j:k,11) = t_end(i)
            rtemp(j:k,12) = vll_start(i)
            rtemp(j:k,13) = weight(i)/duplicate_factor

            itemp(j:k)    = beam(i)
            ltemp(j:k)    = lgc2fo_start(i)
         END DO

         ! Deallocate
         DEALLOCATE(R_start,  PHI_start,  Z_start, &
                    VR_start, VPHI_start, VZ_start, &
                    mass, charge, mu_start, Zatom, &
                    t_end, vll_start, beam, weight, lgc2fo_start)

         ! Reallocate
         k = nparticles_new
         ALLOCATE(  R_start(k), phi_start(k), Z_start(k), &
                    vr_start(k), vphi_start(k), vz_start(k), &
                    mass(k), charge(k), &
                    mu_start(k), Zatom(k), t_end(k), vll_start(k), &
                    beam(k), weight(k), end_state(k), lgc2fo_start(k))

         ! Load arrays
         R_start    = rtemp(:,1)
         PHI_start  = rtemp(:,2)
         Z_start    = rtemp(:,3)
         VR_start   = rtemp(:,4)
         VPHI_start = rtemp(:,5)
         VZ_start   = rtemp(:,6)
         mass       = rtemp(:,7)
         charge     = rtemp(:,8)
         mu_start   = rtemp(:,9)
         Zatom      = rtemp(:,10)
         t_end      = rtemp(:,11)
         VLL_start  = rtemp(:,12)
         weight     = rtemp(:,13)
         beam       = itemp(:)
         lgc2fo_start = ltemp(:)

         ! Deallocate helpers
         DEALLOCATE(itemp,rtemp)

         ! Redefine number of particles
         nparticles = nparticles_new

      ELSE ! Everyone else deallocate

         DEALLOCATE(R_start,  PHI_start,  Z_start, &
                    VR_start, VPHI_start, VZ_start, &
                    mass, charge, mu_start, Zatom, &
                    t_end, vll_start, beam, weight, lgc2fo_start)
      END IF

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(nparticles,1,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
      IF (myworkid /= master) THEN
         ALLOCATE(  R_start(nparticles), phi_start(nparticles), Z_start(nparticles), &
                    vr_start(nparticles), vphi_start(nparticles), vz_start(nparticles), &
                    mass(nparticles), charge(nparticles), &
                    mu_start(nparticles), Zatom(nparticles), t_end(nparticles), vll_start(nparticles), &
                    beam(nparticles), weight(nparticles), end_state(nparticles), lgc2fo_start(nparticles) )
      END IF
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(mu_start,  nparticles, MPI_REAL8,   master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_BCAST(t_end,     nparticles, MPI_REAL8,   master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_BCAST(mass,      nparticles, MPI_REAL8,   master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_BCAST(charge,    nparticles, MPI_REAL8,   master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_BCAST(Zatom,     nparticles, MPI_REAL8,   master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_BCAST(weight,    nparticles, MPI_REAL8,   master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_BCAST(R_start,   nparticles, MPI_REAL8,   master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_BCAST(phi_start, nparticles, MPI_REAL8,   master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_BCAST(Z_start,   nparticles, MPI_REAL8,   master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_BCAST(vll_start, nparticles, MPI_REAL8,   master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_BCAST(beam,      nparticles, MPI_INTEGER, master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_BCAST(VR_start,     nparticles, MPI_REAL8,   master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_BCAST(VPHI_start,   nparticles, MPI_REAL8,   master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_BCAST(VZ_start,     nparticles, MPI_REAL8,   master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_BCAST(lgc2fo_start, nparticles, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
#endif

      ! Now we randomize the gyrophase for paritcles which are FO
      !    We use ltemp as the masking array so that the
      !    first FO particle retains it's initial condition
      IF (myworkid == master) THEN
         i = 1
         ltemp = .not. ltemp 
         DO WHILE (i < nparticles-1)
            IF (ltemp(i) .and. ltemp(i+1)) THEN
                  ltemp(i) = .FALSE.
                  i = i + duplicate_factor
            ELSE
               i = i + 1
            END IF
         END DO
      ELSE
         ALLOCATE(ltemp(nparticles))
      END IF
      CALL MPI_BCAST(ltemp, nparticles, MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_CALC_MYRANGE(MPI_COMM_BEAMS, 1, nparticles, mystart, myend)
      ! zero out previous data
      DO i = 1, mystart-1
         R_START(i)    = 0
         PHI_START(i)  = 0
         Z_START(i)    = 0
         VR_START(i)   = 0
         VPHI_START(i) = 0
         VZ_START(i)   = 0
      END DO
      DO i = myend+1, nparticles
         R_START(i)    = 0
         PHI_START(i)  = 0
         Z_START(i)    = 0
         VR_START(i)   = 0
         VPHI_START(i) = 0
         VZ_START(i)   = 0
      END DO
      ! Now process particles
      DO i = mystart, myend
         IF (ltemp(i)) THEN
            q(1) = R_start(i)
            q(2) = PHI_start(i)
            q(3) = Z_start(i)
            q(4) = VR_start(i)
            q(5) = VPHI_start(i)
            q(6) = VZ_start(i)
            myline = i
            mymass = mass(i)
            mycharge = charge(i)
            CALL beams3d_part2gc(q)
            q(5) = moment
            CALL beams3d_gc2fo(q)
            R_start(i) = q(1)
            PHI_start(i) = q(2)
            Z_start(i) = q(3)
            VR_start(i) = q(4)
            VPHI_start(i) = q(5)
            VZ_start(i) = q(6)
         END IF
      END DO

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      ! Now broadcast to everyone
      IF (myworkid == master) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE, R_start,    nparticles, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, PHI_start,  nparticles, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, Z_start,    nparticles, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, VR_start,   nparticles, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, VPHI_start, nparticles, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, VZ_start,   nparticles, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
      ELSE
         CALL MPI_REDUCE(R_start,      R_start,    nparticles, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(PHI_start,    PHI_start,  nparticles, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(Z_start,      Z_start,    nparticles, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(VR_start,     VR_start,   nparticles, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(VPHI_start,   VPHI_start, nparticles, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(VZ_start,     VZ_start,   nparticles, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
      END IF
      CALL MPI_BCAST(R_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(PHI_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(Z_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(vr_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(vphi_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(vz_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
#endif

      DEALLOCATE(end_state, ltemp)

      RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_duplicate_part
