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
      USE beams3d_physics_mod, ONLY: beams3d_DTRATE, beams3d_MODB, &
                                     beams3d_SFLX, beams3d_DDTRATE, &
                                     beams3d_DDHe3RATE
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
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: n3d
      REAL(rprec) :: maxrateDT, maxrateDDT, maxrateDDHe, &
                     dV, w_total, vpart, dr, dphi, &
                     dz, X1_rand, Y1_rand, Z1_rand, sval, sfactor
      REAL(rprec), DIMENSION(3) :: q
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: X_rand, Y_rand, Z_rand, P_rand
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: temp3d

#if defined(MPI_OPT)
      INTEGER :: numprocs_local, mylocalid, mylocalmaster, mystart, myend
      INTEGER :: MPI_COMM_LOCAL
#endif

      ! Shared memory variables
      INTEGER :: win_rateDT, win_rateDDT, win_rateDDHe, win_B_help, win_l3d
      LOGICAL,          POINTER, DIMENSION(:,:,:) :: l3d
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: B_help
      DOUBLE PRECISION, POINTER, DIMENSION(:,:,:) :: rateDT, rateDDT, rateDDHe

      REAL(rprec),      PARAMETER :: S_LIM_MAX = 1.0 ! Max S if lplasma_only
      DOUBLE PRECISION, PARAMETER :: e_charge  = 1.60217662E-19 !e_c
      DOUBLE PRECISION, PARAMETER :: mHe4      = 6.6464731D-27
      DOUBLE PRECISION, PARAMETER :: mT        = 5.0082671D-27
      DOUBLE PRECISION, PARAMETER :: mp        = 1.6726219D-27
      DOUBLE PRECISION, PARAMETER :: mHe3      = 5.0082335D-27

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
      CALL mpialloc(l3d,      nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_l3d)
      CALL mpialloc(rateDT,   nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_rateDT)
      CALL mpialloc(rateDDT,  nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_rateDDT)
      CALL mpialloc(rateDDHe, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_rateDDHe)

      ! Initialize helpers
      IF (mylocalid == mylocalmaster) THEN
         rateDT = 0; rateDDT = 0; rateDDHe = 0; l3d = .false.
      END IF
      maxrateDT = 0; maxrateDDT = 0; maxrateDDHe = 0

      ! Setup masking
      sfactor = 2
      IF (lplasma_only) sfactor = S_LIM_MAX
      
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
         ! S on full grid
         IF ((i==nr) .or. (j==nphi) .or. (k==nz)) CYCLE
         q = (/raxis(i), phiaxis(j), zaxis(k)/)+0.5*(/hr(i),hp(j),hz(k)/) ! Half grid
         CALL beams3d_DTRATE(q,rateDT(i,j,k))
         maxrateDT = MAX(rateDT(i,j,k),maxrateDT)
         CALL beams3d_DDTRATE(q,rateDDT(i,j,k))
         maxrateDDT = MAX(rateDDT(i,j,k),maxrateDDT)
         CALL beams3d_DDHe3RATE(q,rateDDHe(i,j,k))
         maxrateDDHe = MAX(rateDDHe(i,j,k),maxrateDDHe)
         CALL beams3d_SFLX(q,sval)
         IF (sval < sfactor) l3d(i,j,k) = .true.
         ! We really want the total number of particles in the box (dV=rdrdpdz)
         dV = (raxis(i)+0.5*hr(i))*hr(i)*hp(j)*hz(k)
         rateDT(i,j,k) = rateDT(i,j,k) * dV
         rateDDT(i,j,k) = rateDDT(i,j,k) * dV
         rateDDHe(i,j,k) = rateDDHe(i,j,k) * dV
      END DO
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,maxrateDT,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,maxrateDDT,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,maxrateDDHe,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_LOCAL,ierr_mpi)
#endif

      ! Output reaction rate info
      IF (myworkid == master) THEN
         IF (lverb) THEN
            WRITE(6, '(A,ES11.4,A)') '      max D + T -> He4   rate: ', maxrateDT,' [part/(m^3 s)]'
            IF (.not. lfusion_alpha) THEN
               WRITE(6, '(A,ES11.4,A)') '      max D + D -> T + p rate: ', maxrateDDT,' [part/(m^3 s)]'
               WRITE(6, '(A,ES11.4,A)') '      max D + D -> He3   rate: ', maxrateDDHe,' [part/(m^3 s)]'
            END IF
         END IF
      END IF

      ! From this point forward we act like a NB calc in the initialization sense
      ! Eventually we'll include DD-T and DD-He3 as BEAM(2 and (3))
      nbeams = 4
      IF (lfusion_alpha) nbeams = 1 ! Do alphas only
      nparticles = nbeams*nparticles_start
      ALLOCATE(   R_start(nparticles), phi_start(nparticles), Z_start(nparticles), vll_start(nparticles), &
                  mass(nparticles), charge(nparticles), Zatom(nparticles), &
                  mu_start(nparticles), t_end(nparticles), &
                  beam(nparticles), weight(nparticles), v_neut(3,nparticles))
      v_neut = 0
      ! BEAM 1 D+T -> He4
      E_BEAMS(1) = 3.52E6*e_charge*fusion_scale
      MASS_BEAMS(1) = mHe4
      CHARGE_BEAMS(1) = 2*e_charge
      IF (.not. lfusion_alpha) THEN
         ! BEAM 2 D+D -> T (+p) fast Tritium
         E_BEAMS(2) = 1.01E6*e_charge*fusion_scale
         MASS_BEAMS(2) = mT
         CHARGE_BEAMS(2) = e_charge
         ! BEAM 3 D+D -> (T) +p fast proton
         E_BEAMS(3) = 3.02E6*e_charge*fusion_scale
         MASS_BEAMS(3) = mp
         CHARGE_BEAMS(3) = e_charge
         ! BEAM 4 D+D -> He3
         E_BEAMS(4) = 0.82E6*e_charge*fusion_scale
         MASS_BEAMS(4) = mHe3
         CHARGE_BEAMS(4) = 2*e_charge
      END IF

      ! Now we need to monte-carlo our way to get particles inside the grid
      dr = raxis(nr)-raxis(1)
      dphi = phiaxis(nphi)-phiaxis(1)
      dz = zaxis(nz)-zaxis(1)
      IF (myworkid == master) THEN
         ! Handle that l3d is on half grid
         l3d(:,nphi,:) = l3d(:,1,:)
         ALLOCATE(n3d(nr,nphi,nz))
         n3d = 0
         DO s = 1,nparticles_start
            DO
               CALL RANDOM_NUMBER(X1_RAND)
               CALL RANDOM_NUMBER(Y1_RAND)
               CALL RANDOM_NUMBER(Z1_RAND)
               q(1) = X1_RAND*dr+raxis(1)
               q(2) = Y1_RAND*dphi+phiaxis(1)
               q(3) = Z1_RAND*dz+zaxis(1)
               i = MIN(MAX(COUNT(raxis < q(1)),2),nr-1)
               j = MIN(MAX(COUNT(phiaxis < q(2)),1),nphi-1)
               k = MIN(MAX(COUNT(zaxis < q(3)),2),nz-1)
               IF (l3d(i,j,k) .and. l3d(i+1, j,   k)   .and. l3d(i-1, j,   k) &
                              .and. l3d(i,   j,   k+1) .and. l3d(i,   j,   k-1) &
                              .and. l3d(i+1, j,   k+1) .and. l3d(i+1, j,   k-1) &
                              .and. l3d(i-1, j,   k+1) .and. l3d(i-1, j,   k-1) &
                              .and. l3d(i+1, j+1, k)   .and. l3d(i-1, j+1, k) &
                              .and. l3d(i  , j+1, k+1) .and. l3d(i,   j+1, k-1)) &
                  EXIT
            END DO
            n3d(i,j,k) = n3d(i,j,k) + 1
         END DO

         ! Now for each voxel we calcualte the starting points
         ALLOCATE(X_rand(nparticles),Y_rand(nparticles),Z_rand(nparticles),P_rand(nparticles))
         ! Randomize on uniform grid
         CALL RANDOM_NUMBER(X_RAND)
         CALL RANDOM_NUMBER(Y_RAND)
         CALL RANDOM_NUMBER(Z_RAND)
         ! Uniform about 0
         CALL RANDOM_NUMBER(P_rand)
         !P_rand = 2*(P_rand-0.5)
         P_rand = cos(pi*P_rand)
         k1 = 1
         ! Do the D-T -> H4 reaction
         vpart = sqrt(2*E_BEAMS(1)/mHe4)
         DO s = 1,nr*nphi*nz
            i = MOD(s-1,nr)+1
            j = MOD(s-1,nr*nphi)
            j = FLOOR(REAL(j) / REAL(nr))+1
            k = CEILING(REAL(s) / REAL(nr*nphi))
            IF (n3d(i,j,k)==0) CYCLE
            k2 = n3d(i,j,k)+k1-1
            R_start(k1:k2)   =   raxis(i) + X_rand(k1:k2)*hr(i)
            phi_start(k1:k2) = phiaxis(j) + Y_rand(k1:k2)*hp(j)
            Z_start(k1:k2)   =   zaxis(k) + Z_rand(k1:k2)*hz(k)
            vll_start(k1:k2) = vpart*P_rand(k1:k2) 
            mu_start(k1:k2)  = E_BEAMS(1) ! Total energy for now
            beam(k1:k2)      = 1
            mass(k1:k2)      = MASS_BEAMS(1)
            charge(k1:k2)    = CHARGE_BEAMS(1)
            Zatom(k1:k2)     = 2
            t_end(k1:k2)     = t_end_in(1)
            weight(k1:k2)    = rateDT(i,j,k)/n3d(i,j,k)
            k1 = k2+1
         END DO
         IF (.not. lfusion_alpha) THEN
            ! Do the D-D -> T reaction
            vpart = sqrt(2*E_BEAMS(2)/mT)
            DO s = 1,nr*nphi*nz
               i = MOD(s-1,nr)+1
               j = MOD(s-1,nr*nphi)
               j = FLOOR(REAL(j) / REAL(nr))+1
               k = CEILING(REAL(s) / REAL(nr*nphi))
               IF (n3d(i,j,k)==0) CYCLE
               k2 = n3d(i,j,k)+k1-1
               R_start(k1:k2)   =   raxis(i) + X_rand(k1:k2)*hr(i)
               phi_start(k1:k2) = phiaxis(j) + Y_rand(k1:k2)*hp(j)
               Z_start(k1:k2)   =   zaxis(k) + Z_rand(k1:k2)*hz(k)
               vll_start(k1:k2) = vpart*P_rand(k1:k2) 
               mu_start(k1:k2)  = E_BEAMS(2) ! Total energy for now
               beam(k1:k2)      = 2
               mass(k1:k2)      = MASS_BEAMS(2)
               charge(k1:k2)    = CHARGE_BEAMS(2)
               Zatom(k1:k2)     = 1
               t_end(k1:k2)     = t_end_in(1)
               weight(k1:k2)    = rateDDT(i,j,k)/n3d(i,j,k)
               k1 = k2+1
            END DO
            ! Do the D-D -> p reaction
            vpart = sqrt(2*E_BEAMS(3)/mp)
            DO s = 1,nr*nphi*nz
               i = MOD(s-1,nr)+1
               j = MOD(s-1,nr*nphi)
               j = FLOOR(REAL(j) / REAL(nr))+1
               k = CEILING(REAL(s) / REAL(nr*nphi))
               IF (n3d(i,j,k)==0) CYCLE
               k2 = n3d(i,j,k)+k1-1
               R_start(k1:k2)   =   raxis(i) + X_rand(k1:k2)*hr(i)
               phi_start(k1:k2) = phiaxis(j) + Y_rand(k1:k2)*hp(j)
               Z_start(k1:k2)   =   zaxis(k) + Z_rand(k1:k2)*hz(k)
               vll_start(k1:k2) = vpart*P_rand(k1:k2) 
               mu_start(k1:k2)  = E_BEAMS(3) ! Total energy for now
               beam(k1:k2)      = 3
               mass(k1:k2)      = MASS_BEAMS(3)
               charge(k1:k2)    = CHARGE_BEAMS(3)
               Zatom(k1:k2)     = 1
               t_end(k1:k2)     = t_end_in(1)
               weight(k1:k2)    = rateDDT(i,j,k)/n3d(i,j,k)
               k1 = k2+1
            END DO
            ! Do the D-D -> He3 reaction
            vpart = sqrt(2*E_BEAMS(4)/mHe3)
            DO s = 1,nr*nphi*nz
               i = MOD(s-1,nr)+1
               j = MOD(s-1,nr*nphi)
               j = FLOOR(REAL(j) / REAL(nr))+1
               k = CEILING(REAL(s) / REAL(nr*nphi))
               IF (n3d(i,j,k)==0) CYCLE
               k2 = n3d(i,j,k)+k1-1
               R_start(k1:k2)   =   raxis(i) + X_rand(k1:k2)*hr(i)
               phi_start(k1:k2) = phiaxis(j) + Y_rand(k1:k2)*hp(j)
               Z_start(k1:k2)   =   zaxis(k) + Z_rand(k1:k2)*hz(k)
               vll_start(k1:k2) = vpart*P_rand(k1:k2) 
               mu_start(k1:k2)  = E_BEAMS(4) ! Total energy for now
               beam(k1:k2)      = 4
               mass(k1:k2)      = MASS_BEAMS(4)
               charge(k1:k2)    = CHARGE_BEAMS(4)
               Zatom(k1:k2)     = 2
               t_end(k1:k2)     = t_end_in(1)
               weight(k1:k2)    = rateDDHe(i,j,k)/n3d(i,j,k)
               k1 = k2+1
            END DO
         END IF
         ! Deallocate and cleanup
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
      CALL mpidealloc(l3d,win_l3d)
      CALL mpidealloc(rateDT,win_rateDT)
      CALL mpidealloc(rateDDT,win_rateDDT)
      CALL mpidealloc(rateDDHe,win_rateDDHe)


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
