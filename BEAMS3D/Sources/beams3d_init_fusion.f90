!-----------------------------------------------------------------------
!     Module:        beams3d_init_fusion
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          09/30/2020
!     Description:   This subroutine initializes the particle startring
!                    points based on nuclear reaction rates.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_init_fusion
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_runtime
      USE beams3d_grid, ONLY: nr, nphi, nz, hr, hp, hz, raxis, zaxis, &
                              phiaxis
      USE beams3d_lines, ONLY: nparticles, partvmax
      USE beams3d_physics_mod, ONLY: beams3d_DTRATE, beams3d_MODB
      USE mpi_sharmem
      USE mpi_params
      USE mpi_inc

!-----------------------------------------------------------------------
!     Local Variables
!          s/i/j/k/k1/k2    Index variables
!          maxrate          Maximum DT rate
!          Ntotal           Total number of particles generated
!          dV               Differential volume of voxel
!          q                Array helper
!          w_temp           Weight helper
!          X_rand           Random number arrays (X,Y,Z,P)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: s,i,j,k,k1,k2
      REAL(rprec) :: maxrateDT, dV, Ntotal, w_total, vpart
      REAL(rprec), DIMENSION(3) :: q
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: X_rand, Y_rand, Z_rand, &
                                                P_rand
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: temp3d

#if defined(MPI_OPT)
      INTEGER :: numprocs_local, mylocalid, mylocalmaster, mystart, myend
      INTEGER :: MPI_COMM_LOCAL
#endif

      ! Shared memory variables
      INTEGER :: win_rateDT, win_B_help
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: B_help
      DOUBLE PRECISION, POINTER, DIMENSION(:,:,:) :: rateDT

      DOUBLE PRECISION, PARAMETER :: e_charge = 1.60217662E-19 !e_c
      DOUBLE PRECISION, PARAMETER :: mHe4     = 6.6464731D-27
      DOUBLE PRECISION, PARAMETER :: mT       = 5.0082671D-27
      DOUBLE PRECISION, PARAMETER :: mHe3     = 5.0082335D-27

!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      ! Divide up Work
      mylocalid = myworkid
      numprocs_local = 1
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_LOCAL, ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_LOCAL, mylocalid, ierr_mpi )              ! MPI
      CALL MPI_COMM_SIZE( MPI_COMM_LOCAL, numprocs_local, ierr_mpi )          ! MPI
#endif
      mylocalmaster = master

      CALL init_random_seed
      IF (lverb) THEN
         WRITE(6, '(A)') '----- INITIALIZING FUSION REACTIONS -----'
         WRITE(6, '(A,I6)') '      nparticles_start: ', nparticles_start
         CALL FLUSH(6)
      END IF
     
      ! We need to first define the reaction rate over the grid
      CALL mpialloc(rateDT, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_rateDT)

      ! Initialize helpers
      IF (mylocalid == mylocalmaster) THEN
         rateDT = 0
      END IF
      Ntotal = 0; maxrateDT = 0
      
      ! Break up the Work
      CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, nr*nphi*nz, mystart, myend)
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
#endif
      DO s = mystart, myend
         i = MOD(s-1,nr)+1
         j = MOD(s-1,nr*nphi)
         j = FLOOR(REAL(j) / REAL(nr))+1
         k = CEILING(REAL(s) / REAL(nr*nphi))
         IF ((i==nr) .or. (j==nphi) .or. (k==nz)) CYCLE
         q = (/raxis(i), phiaxis(j), zaxis(k)/)+0.5*(/hr(i),hp(j),hz(k)/) ! Half grid
         CALL beams3d_DTRATE(q,rateDT(i,j,k))
         maxrateDT = MAX(rateDT(i,j,k),maxrateDT)
         ! We really want the total number of particles in the box (dV=rdrdpdz)
         dV = (raxis(i)+0.5*hr(i))*hr(i)*hp(j)*hz(k)
         rateDT(i,j,k) = rateDT(i,j,k) * dV
         Ntotal  = Ntotal + rateDT(i,j,k)
      END DO
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,Ntotal,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,maxrateDT,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_LOCAL,ierr_mpi)
#endif

      ! Calculate a new nparticles_start
      IF (myworkid == master) THEN
         DO s = 1, nr*nphi*nz
            i = MOD(s-1,nr)+1
            j = MOD(s-1,nr*nphi)
            j = FLOOR(REAL(j) / REAL(nr))+1
            k = CEILING(REAL(s) / REAL(nr*nphi))
            WRITE(327,*) i,j,k,rateDT(i,j,k)
         END DO
         PRINT *,NTOTAL
         rateDT = nparticles_start*rateDT/Ntotal ! Particles per box
         !WHERE(rateDT > 0) rateDT = MAX(rateDT,1.0)
         rateDT = FLOOR(rateDT)
         nparticles_start = SUM(SUM(SUM(rateDT,3),2),1)
         IF (lverb) THEN
            WRITE(6, '(A,ES11.4,A,I8)') '          max D-T rate ', maxrateDT,'; nparticles = ',nparticles_start
         END IF
      END IF

#if defined(MPI_OPT)
      CALL MPI_BCAST(nparticles_start,1,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
#endif


      ! Particle weighting
      w_total = Ntotal/nparticles_start !part/s

      ! From this point forward we act like a NB calc in the initialization sense
      ! Eventually we'll include DD-T and DD-He3 as BEAM(2 and (3))
      nbeams = 1
      nparticles = nbeams*nparticles_start
      ALLOCATE(   R_start(nparticles), phi_start(nparticles), Z_start(nparticles), vll_start(nparticles), &
                  mass(nparticles), charge(nparticles), Zatom(nparticles), &
                  mu_start(nparticles), t_end(nparticles), &
                  beam(nparticles), weight(nparticles), v_neut(3,nparticles))
      v_neut = 0
      E_BEAMS(1) = 3.52E6*e_charge !D+T -> He4
      E_BEAMS(2) = 1.01E6*e_charge !D+D -> T (+ p) !fast Tritium
      E_BEAMS(3) = 3.02E6*e_charge !D+D -> (T) + p !fast proton
      E_BEAMS(4) = 0.82E6*e_charge !D+D -> He3
      ! We need to define the reaction rate 
      IF (myworkid == master) THEN
         ALLOCATE(X_rand(nparticles_start),Y_rand(nparticles_start),Z_rand(nparticles_start),P_rand(nparticles_start))
         ! Randomize on uniform grid
         CALL RANDOM_NUMBER(X_RAND)
         CALL RANDOM_NUMBER(Y_RAND)
         CALL RANDOM_NUMBER(Z_RAND)
         ! Gausian about 0
         CALL gauss_rand(nparticles_start, P_rand)
         k1 = 1
         ! Normalize the rate
         DO s = 1, nr*nphi*nz
            i = MOD(s-1,nr)+1
            j = MOD(s-1,nr*nphi)
            j = FLOOR(REAL(j) / REAL(nr))+1
            k = CEILING(REAL(s) / REAL(nr*nphi))
            IF (rateDT(i,j,k)==0) CYCLE
            vpart = sqrt(2*E_BEAMS(1)/mHe4)
            ! Randomize starting locations
            k2 = k1+rateDT(i,j,k)
            R_start(k1:k2)   =   raxis(i) + X_rand(k1:k2)*hr(i)
            phi_start(k1:k2) = phiaxis(j) + Y_rand(k1:k2)*hp(j)
            Z_start(k1:k2)   =   zaxis(k) + Z_rand(k1:k2)*hz(k)
            vll_start(k1:k2) = vpart*P_rand(k1:k2) ! Need to add scaling factor here
            mu_start(k1:k2)  = E_BEAMS(1) ! Total energy for now
            beam(k1:k2)      = 1
            mass(k1:k2)      = mHe4
            charge(k1:k2)    = 2*e_charge
            Zatom(k1:k2)     = 2
            t_end(k1:k2)     = t_end_in(1)
            weight(k1:k2)    = w_total
            k1 = k2+1
         END DO
         DEALLOCATE(X_rand,Y_rand,Z_rand,P_rand)

         partvmax = SQRT(MAXVAL(2*mu_start/mass)) ! Partvmax from max energy
      END IF

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(partvmax,1,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(mu_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(t_end,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(mass,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(charge,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(Zatom,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(weight,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(R_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(phi_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(Z_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(vll_start,nparticles,MPI_REAL8, master, MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_BCAST(beam,nparticles,MPI_INTEGER, master, MPI_COMM_BEAMS,ierr_mpi)
#endif
      CALL mpidealloc(rateDT,win_rateDT)


      ! Now we calcualte values over particles (share the work)
      CALL mpialloc(B_help, nparticles, myid_sharmem, 0, MPI_COMM_SHARMEM, win_B_help)

      ! Initialize helpers
      IF (mylocalid == mylocalmaster) THEN
         B_help = 0
      END IF
      
      ! Break up the Work
      CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, nparticles, mystart, myend)
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
#endif
      DO s = mystart, myend
         q = (/R_start(s), phi_start(s), Z_start(s)/)
         CALL beams3d_MODB(q,B_help(s))
      END DO
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
#endif

      mu_start = (mu_start - 0.5*mass*vll_start*vll_start) / B_help

      CALL mpidealloc(B_help,win_B_help)

      RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_init_fusion
