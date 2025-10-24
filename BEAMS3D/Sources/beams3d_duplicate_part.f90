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
#if defined(MPI_OPT)
      INTEGER :: numprocs_local, mylocalid, mylocalmaster
      INTEGER :: MPI_COMM_LOCAL
#endif
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      STOP 'This is broken'
      ! Divide up Work
      mylocalid = myworkid
      numprocs_local = 1
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_LOCAL, ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_LOCAL, mylocalid, ierr_mpi )              ! MPI
      CALL MPI_COMM_SIZE( MPI_COMM_LOCAL, numprocs_local, ierr_mpi )          ! MPI
#endif
      mylocalmaster = master
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
      END IF

      ! Deallocate
      CALL mpidealloc(R_start, win_R_start)
      CALL mpidealloc(PHI_start, win_PHI_start)
      CALL mpidealloc(Z_start, win_Z_start)
      CALL mpidealloc(vr_start, win_vr_start)
      CALL mpidealloc(vphi_start, win_vphi_start)
      CALL mpidealloc(vz_start, win_vz_start)
      CALL mpidealloc(mass, win_mass)
      CALL mpidealloc(charge, win_charge)
      CALL mpidealloc(mu_start, win_mu_start)
      CALL mpidealloc(Zatom, win_Zatom)
      CALL mpidealloc(t_end, win_t_end)
      CALL mpidealloc(vll_start, win_vll_start)
      CALL mpidealloc(beam, win_beam)
      CALL mpidealloc(weight, win_weight)
      CALL mpidealloc(lgc2fo_start, win_lgc2fo_start)
      CALL mpidealloc(end_state, win_end_state)

      ! Allocate
      CALL MPI_BCAST(k, 1, MPI_INTEGER, master, MPI_COMM_BEAMS, ierr_mpi)
      CALL mpialloc(R_start,      k, myid_sharmem, 0, MPI_COMM_SHARMEM, win_R_start)
      CALL mpialloc(PHI_start,    k, myid_sharmem, 0, MPI_COMM_SHARMEM, win_PHI_start)
      CALL mpialloc(Z_start,      k, myid_sharmem, 0, MPI_COMM_SHARMEM, win_Z_start)
      CALL mpialloc(vr_start,     k, myid_sharmem, 0, MPI_COMM_SHARMEM, win_vr_start)
      CALL mpialloc(vphi_start,   k, myid_sharmem, 0, MPI_COMM_SHARMEM, win_vphi_start)
      CALL mpialloc(vz_start,     k, myid_sharmem, 0, MPI_COMM_SHARMEM, win_vz_start)
      CALL mpialloc(mass,         k, myid_sharmem, 0, MPI_COMM_SHARMEM, win_mass)
      CALL mpialloc(charge,       k, myid_sharmem, 0, MPI_COMM_SHARMEM, win_charge)
      CALL mpialloc(mu_start,     k, myid_sharmem, 0, MPI_COMM_SHARMEM, win_mu_start)
      CALL mpialloc(Zatom,        k, myid_sharmem, 0, MPI_COMM_SHARMEM, win_Zatom)
      CALL mpialloc(t_end,        k, myid_sharmem, 0, MPI_COMM_SHARMEM, win_t_end)
      CALL mpialloc(vll_start,    k, myid_sharmem, 0, MPI_COMM_SHARMEM, win_vll_start)
      CALL mpialloc(beam,         k, myid_sharmem, 0, MPI_COMM_SHARMEM, win_beam)
      CALL mpialloc(weight,       k, myid_sharmem, 0, MPI_COMM_SHARMEM, win_weight)
      CALL mpialloc(lgc2fo_start, k, myid_sharmem, 0, MPI_COMM_SHARMEM, win_lgc2fo_start)
      CALL mpialloc(end_state,    k, myid_sharmem, 0, MPI_COMM_SHARMEM, win_end_state)

         ! Deallocate
         !DEALLOCATE(R_start,  PHI_start,  Z_start, &
         !           VR_start, VPHI_start, VZ_start, &
         !           mass, charge, mu_start, Zatom, &
         !           t_end, vll_start, beam, weight, lgc2fo_start)

         ! Reallocate
         !k = nparticles_new
         !ALLOCATE(  R_start(k), phi_start(k), Z_start(k), &
         !           vr_start(k), vphi_start(k), vz_start(k), &
         !           mass(k), charge(k), &
         !           mu_start(k), Zatom(k), t_end(k), vll_start(k), &
         !           beam(k), weight(k), end_state(k), lgc2fo_start(k))

      IF (myworkid == master) THEN
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

      !ELSE ! Everyone else deallocate
      !   DEALLOCATE(R_start,  PHI_start,  Z_start, &
      !              VR_start, VPHI_start, VZ_start, &
      !              mass, charge, mu_start, Zatom, &
      !              t_end, vll_start, beam, weight, lgc2fo_start)
      END IF

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(nparticles,1,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
      i = MPI_UNDEFINED
      IF (myid_sharmem == master) i = 0
      CALL MPI_COMM_SPLIT( MPI_COMM_BEAMS,i,myworkid,MPI_COMM_LOCAL,ierr_mpi)
      IF (myid_sharmem == master) THEN
         CALL MPI_BCAST(     R_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(   PHI_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(     Z_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(   vll_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(    mu_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(    vr_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(  vphi_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(    vz_start, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(       t_end, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(        mass, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(      charge, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(       Zatom, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(      weight, nparticles,   MPI_REAL8, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(        beam, nparticles, MPI_INTEGER, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_BCAST(lgc2fo_start, nparticles, MPI_LOGICAL, master, MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      END IF
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
      !DO i = 1, mystart-1
      !   R_START(i)    = 0
      !   PHI_START(i)  = 0
      !   Z_START(i)    = 0
      !   VR_START(i)   = 0
      !   VPHI_START(i) = 0
      !   VZ_START(i)   = 0
      !END DO
      !DO i = myend+1, nparticles
      !   R_START(i)    = 0
      !   PHI_START(i)  = 0
      !   Z_START(i)    = 0
      !   VR_START(i)   = 0
      !   VPHI_START(i) = 0
      !   VZ_START(i)   = 0
      !END DO
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
