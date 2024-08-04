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
                              phiaxis, S4D, U4D, dexionHe3, &
                              NEUTRONS_ARR, win_NEUTRONS, &
                              E_NEUTRONS, win_E_NEUTRONS
      USE beams3d_lines, ONLY: nparticles, partvmax
      USE beams3d_physics_mod, ONLY: beams3d_DTRATE, beams3d_MODB, &
                                     beams3d_SFLX, beams3d_DDTRATE, &
                                     beams3d_DDHe3RATE, beams3d_DHe3RATE
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
      LOGICAL :: lfusion_DHe3
      INTEGER :: s,i,j,k,k1,k2, nr1, nphi1, nz1, l, nfp
      INTEGER :: minik(2)
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: n3d
      REAL(rprec) :: maxrateDT, maxrateDDT, maxrateDDHe, maxrateDHe3, &
                     dV, w_total, vpart, dr, dphi, &
                     dz, X1_rand, Y1_rand, Z1_rand, sval, sfactor, &
                     rtemp, roffset, utemp, rerr, uerr, s_fullorbit
      REAL(rprec), DIMENSION(3) :: q
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: X_rand, Y_rand, Z_rand, &
                                                P_rand, R_axis, Z_axis
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: f2d
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: temp3d

#if defined(MPI_OPT)
      INTEGER :: numprocs_local, mylocalid, mylocalmaster, mystart, myend
      INTEGER :: MPI_COMM_LOCAL
#endif

      ! Shared memory variables
      INTEGER :: win_rateDT, win_rateDDT, win_rateDDHe, win_B_help, &
                 win_l3d, win_rateDHe3
      LOGICAL,          POINTER, DIMENSION(:,:,:) :: l3d
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: B_help
      DOUBLE PRECISION, POINTER, DIMENSION(:,:,:) :: rateDT, rateDDT, rateDDHe, rateDHe3

      REAL(rprec),      PARAMETER :: S_LIM_MAX = 0.90 ! Max S if lplasma_only
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

      ! We Work with half grid
      nr1   = nr - 1
      nphi1 = nphi - 1
      nz1   = nz - 1
      nfp = NINT(pi2/phiaxis(nphi))
     
      ! We need to first define the reaction rate over the grid
      CALL mpialloc(l3d,      nr1, nphi1, nz1, myid_sharmem, 0, MPI_COMM_SHARMEM, win_l3d)
      CALL mpialloc(rateDT,   nr1, nphi1, nz1, myid_sharmem, 0, MPI_COMM_SHARMEM, win_rateDT)
      CALL mpialloc(rateDDT,  nr1, nphi1, nz1, myid_sharmem, 0, MPI_COMM_SHARMEM, win_rateDDT)
      CALL mpialloc(rateDDHe, nr1, nphi1, nz1, myid_sharmem, 0, MPI_COMM_SHARMEM, win_rateDDHe)
      CALL mpialloc(rateDHe3, nr1, nphi1, nz1, myid_sharmem, 0, MPI_COMM_SHARMEM, win_rateDHe3)

      ! Initialize helpers
      IF (mylocalid == mylocalmaster) THEN
         rateDT = 0; rateDDT = 0; rateDDHe = 0; rateDHe3 = 0; l3d = .false.
      END IF
      maxrateDT = 0; maxrateDDT = 0; maxrateDDHe = 0; maxrateDHe3 = 0
      lfusion_DHe3 = (dexionHe3 > 0 .and. lfusion_He3)

      ! Setup masking
      sfactor = 1
      IF (lplasma_only) sfactor = S_LIM_MAX
      
      ! Break up the Work
      CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, nr1*nphi1*nz1, mystart, myend)
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
#endif

      ! Calc DT rate
      IF (lfusion_alpha) THEN
         DO s = mystart, myend
            i = MOD(s-1,nr1)+1
            j = MOD(s-1,nr1*nphi1)
            j = FLOOR(REAL(j) / REAL(nr1))+1
            k = CEILING(REAL(s) / REAL(nr1*nphi1))
            q = (/raxis(i), phiaxis(j), zaxis(k)/)+0.5*(/hr(i),hp(j),hz(k)/) ! Half grid
            CALL beams3d_DTRATE(q,rateDT(i,j,k))
            maxrateDT = MAX(rateDT(i,j,k),maxrateDT)
         END DO
      END IF

      ! Calc DD -> T + p
      IF (lfusion_tritium .or. lfusion_proton) THEN
         DO s = mystart, myend
            i = MOD(s-1,nr1)+1
            j = MOD(s-1,nr1*nphi1)
            j = FLOOR(REAL(j) / REAL(nr1))+1
            k = CEILING(REAL(s) / REAL(nr1*nphi1))
            q = (/raxis(i), phiaxis(j), zaxis(k)/)+0.5*(/hr(i),hp(j),hz(k)/) ! Half grid
            CALL beams3d_DDTRATE(q,rateDDT(i,j,k))
            maxrateDDT = MAX(rateDDT(i,j,k),maxrateDDT)
         END DO
      END IF

      ! Calc DD -> He3
      IF (lfusion_He3) THEN
         DO s = mystart, myend
            i = MOD(s-1,nr1)+1
            j = MOD(s-1,nr1*nphi1)
            j = FLOOR(REAL(j) / REAL(nr1))+1
            k = CEILING(REAL(s) / REAL(nr1*nphi1))
            q = (/raxis(i), phiaxis(j), zaxis(k)/)+0.5*(/hr(i),hp(j),hz(k)/) ! Half grid
            CALL beams3d_DDHe3RATE(q,rateDDHe(i,j,k))
            maxrateDDHe = MAX(rateDDHe(i,j,k),maxrateDDHe)
         END DO
      END IF

      ! Calc DHe3 -> He4 + p
      IF (lfusion_DHe3) THEN ! should probably be it's own thing
         DO s = mystart, myend
            i = MOD(s-1,nr1)+1
            j = MOD(s-1,nr1*nphi1)
            j = FLOOR(REAL(j) / REAL(nr1))+1
            k = CEILING(REAL(s) / REAL(nr1*nphi1))
            q = (/raxis(i), phiaxis(j), zaxis(k)/)+0.5*(/hr(i),hp(j),hz(k)/) ! Half grid
            CALL beams3d_DHe3RATE(q,rateDHe3(i,j,k))
            maxrateDHe3 = MAX(rateDHe3(i,j,k),maxrateDHe3)
         END DO
      END IF

      ! l3d array and volume normalization
      DO s = mystart, myend
         i = MOD(s-1,nr1)+1
         j = MOD(s-1,nr1*nphi1)
         j = FLOOR(REAL(j) / REAL(nr1))+1
         k = CEILING(REAL(s) / REAL(nr1*nphi1))
         q = (/raxis(i), phiaxis(j), zaxis(k)/)+0.5*(/hr(i),hp(j),hz(k)/) ! Half grid
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
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,maxrateDHe3,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_LOCAL,ierr_mpi)
#endif

      ! Output reaction rate info
      IF (myworkid == master) THEN
         IF (lverb) THEN
            IF (lfusion_alpha) THEN
               WRITE(6, '(A,ES11.4,A)') '      max D + T -> He4   rate: ', maxrateDT,' [part/(m^3 s)]'
               WRITE(6, '(A,F7.2,A)')   '          He4(3.52 MeV) Power: ', SUM(rateDT)*(3.52)*e_charge*nfp,' [MW]'
               WRITE(6, '(A,F7.2,A)')   '           n(14.06 MeV) Power: ', SUM(rateDT)*(14.06)*e_charge*nfp,' [MW]'
            END IF
            IF (lfusion_tritium .or. lfusion_proton) THEN
               WRITE(6, '(A,ES11.4,A)') '      max D + D -> T + p rate: ', maxrateDDT,' [part/(m^3 s)]'
               WRITE(6, '(A,F7.2,A)')   '            T(1.01 MeV) Power: ', SUM(rateDDT)*(1.01)*e_charge*nfp,' [MW]'
               WRITE(6, '(A,F7.2,A)')   '            p(3.02 MeV) Power: ', SUM(rateDDT)*(3.02)*e_charge*nfp,' [MW]'
            END IF
            IF (lfusion_He3) THEN
               WRITE(6, '(A,ES11.4,A)') '      max D + D -> He3   rate: ', maxrateDDHe,' [part/(m^3 s)]'
               WRITE(6, '(A,F7.2,A)')   '          He3(0.82 MeV) Power: ', SUM(rateDDHe)*(0.82)*e_charge*nfp,' [MW]'
               WRITE(6, '(A,F7.2,A)')   '            n(2.45 MeV) Power: ', SUM(rateDDHe)*(2.45)*e_charge*nfp,' [MW]'
            END IF
            IF (lfusion_DHe3) THEN
               WRITE(6, '(A,ES11.4,A)') '      max D + He3 -> He4 + p   rate: ', maxrateDHe3,' [part/(m^3 s)]'
               WRITE(6, '(A,F7.2,A)')   '                 He4(3.6 MeV) Power: ', SUM(rateDHe3)*(3.60)*e_charge*nfp,' [MW]'
               WRITE(6, '(A,F7.2,A)')   '                  p(14.7 MeV) Power: ', SUM(rateDHe3)*(14.7)*e_charge*nfp,' [MW]'
            END IF
         END IF
      END IF

      ! From this point forward we act like a NB calc in the initialization sense
      nbeams = 0
      IF (lfusion_alpha) nbeams = nbeams + 1
      IF (lfusion_tritium) nbeams = nbeams + 1
      IF (lfusion_proton) nbeams = nbeams + 1
      IF (lfusion_He3) nbeams = nbeams + 1
      IF (lfusion_DHe3) nbeams = nbeams + 2 
      nparticles = nbeams*nparticles_start
      ALLOCATE(   R_start(nparticles), phi_start(nparticles), Z_start(nparticles), vll_start(nparticles), &
                  mass(nparticles), charge(nparticles), Zatom(nparticles), &
                  mu_start(nparticles), t_end(nparticles), &
                  beam(nparticles), weight(nparticles), &
                  vr_start(nparticles), vphi_start(nparticles), vz_start(nparticles), &
                  lgc2fo_start(nparticles))
      CALL mpialloc(NEUTRONS_ARR, 2, nr1, nphi1, nz1, myid_sharmem, 0, MPI_COMM_SHARMEM, win_NEUTRONS)
      CALL mpialloc(E_NEUTRONS, 2, myid_sharmem, 0, MPI_COMM_SHARMEM, win_E_NEUTRONS)
      IF (myworkid == master) THEN
         NEUTRONS_ARR(1,1:nr1,1:nphi1,1:nz1) = rateDT
         NEUTRONS_ARR(2,1:nr1,1:nphi1,1:nz1) = rateDDHe
         E_NEUTRONS(1) = 14.06E6
         E_NEUTRONS(2) = 2.45E6
      END IF
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
#endif
      CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, nr1*nphi1*nz1, mystart, myend)
      ! Add volume back into neutrons
      DO s = mystart, myend
         i = MOD(s-1,nr1)+1
         j = MOD(s-1,nr1*nphi1)
         j = FLOOR(REAL(j) / REAL(nr1))+1
         k = CEILING(REAL(s) / REAL(nr1*nphi1))
         ! We have the n/s but we need to add the volume back in (dV=rdrdpdz)
         dV = (raxis(i)+0.5*hr(i))*hr(i)*hp(j)*hz(k)
         NEUTRONS_ARR(:,i,j,k) = NEUTRONS_ARR(:,i,j,k)/dV
      END DO

      ! We set this true because we assume all particles generated are gyrocenters. Should a
      ! particle have RHO > RHO_FO then we don't follow the GC and follow_fo takes care of
      ! the GC2FO step.  Otherwise we follow GC like normal.
      lgc2fo_start(:) = .TRUE.
      vr_start=0; vphi_start=0; vz_start=0; i = 0
      IF (lfusion_alpha) THEN
         ! BEAM D+T -> He4
         i = i + 1
         E_BEAMS(i) = 3.52E6*e_charge*fusion_scale
         MASS_BEAMS(i) = mHe4
         CHARGE_BEAMS(i) = 2*e_charge
      END IF
      IF (lfusion_tritium) THEN
         ! BEAM D+D -> T (+p) fast Tritium
         i = i + 1
         E_BEAMS(i) = 1.01E6*e_charge*fusion_scale
         MASS_BEAMS(i) = mT
         CHARGE_BEAMS(i) = e_charge
      END IF
      IF (lfusion_proton) THEN
         ! BEAM D+D -> (T) +p fast proton
         i = i + 1
         E_BEAMS(i) = 3.02E6*e_charge*fusion_scale
         MASS_BEAMS(i) = mp
         CHARGE_BEAMS(i) = e_charge
      END IF
      IF (lfusion_He3) THEN
         ! BEAM D+D -> He3
         i = i + 1
         E_BEAMS(i) = 0.82E6*e_charge*fusion_scale
         MASS_BEAMS(i) = mHe3
         CHARGE_BEAMS(i) = 2*e_charge
      END IF
      IF (lfusion_DHe3) THEN
         ! BEAM D+He3 -> He4
         i = i + 1
         E_BEAMS(i) = 3.6E6*e_charge*fusion_scale
         MASS_BEAMS(i) = mHe4
         CHARGE_BEAMS(i) = 2*e_charge
         ! BEAM D+He3 -> p
         i = i + 1
         E_BEAMS(i) = 14.7E6*e_charge*fusion_scale
         MASS_BEAMS(i) = mp
         CHARGE_BEAMS(i) = e_charge
      END IF

      ! Now we need to monte-carlo our way to get particles inside the grid
      dr = raxis(nr)-raxis(1)
      dphi = phiaxis(nphi)-phiaxis(1)
      dz = zaxis(nz)-zaxis(1)
      IF (myworkid == master) THEN
         ! We'll calc the roffset helper
         roffset=0.0
         ! Handle that l3d is on half grid
         ALLOCATE(n3d(nr1,nphi1,nz1))
         n3d = 0
         ! Do a check to make sure we asked for enough particles
         k1 = COUNT(l3d)
         IF ((nparticles_start < k1) .and. lverb) THEN
            WRITE(6,'(A)') '  Number of markers less than needed '
            WRITE(6,'(A,I7)') '     Suggest: ', k1
            k2 = nparticles_start
         ELSE
            WHERE(l3d) n3d=1
            k2 = nparticles_start-k1
         END IF
         ALLOCATE(f2d(nr1,nz1))
         DO s = 1,k2
         	DO
               ! Random numbers
               CALL RANDOM_NUMBER(X1_rand) ! s
               CALL RANDOM_NUMBER(Y1_rand) ! phi
               CALL RANDOM_NUMBER(Z1_rand) ! u
               ! Calc j dex
               j = MIN(MAX(NINT(Y1_rand*nphi),1),nphi1)
               ! Calc u = [0,2*pi]
               utemp = Z1_rand*pi2
               ! Calc sval =[0,1]
               rtemp = ABS(X1_rand*X1_rand-roffset)
               ! Now we need to find the proper point
               f2d = ((U4D(1,:,j,:)-utemp)**2)*((S4D(1,:,j,:)-rtemp)**2)
               minik = MINLOC(f2d)
               i = MIN(MAX(minik(1),2),nr1); k = MIN(MAX(minik(2),2),nz1)
               IF (l3d(i,j,k)) EXIT
            END DO
            n3d(i,j,k) = n3d(i,j,k) + 1
         END DO
         DEALLOCATE(f2d)

         ! Now for each voxel we calcualte the starting points
         ALLOCATE(X_rand(nparticles),Y_rand(nparticles),Z_rand(nparticles),P_rand(nparticles))
         ! Randomize on uniform grid
         CALL RANDOM_NUMBER(X_RAND)
         CALL RANDOM_NUMBER(Y_RAND)
         CALL RANDOM_NUMBER(Z_RAND)
         ! Uniform about 0
         CALL RANDOM_NUMBER(P_rand)
         P_rand = 2*(P_rand-0.5)
         k1 = 1 ; l = 0
         IF (lfusion_alpha) THEN
            ! Do the D-T -> H4 reaction
            l = l + 1
            vpart = sqrt(2*E_BEAMS(l)/mHe4)
            DO s = 1,nr1*nphi1*nz1
               i = MOD(s-1,nr1)+1
               j = MOD(s-1,nr1*nphi1)
               j = FLOOR(REAL(j) / REAL(nr1))+1
               k = CEILING(REAL(s) / REAL(nr1*nphi1))
               IF (n3d(i,j,k)==0) CYCLE
               k2 = n3d(i,j,k)+k1-1
               R_start(k1:k2)   =   raxis(i) + X_rand(k1:k2)*hr(i)
               phi_start(k1:k2) = phiaxis(j) + Y_rand(k1:k2)*hp(j)
               Z_start(k1:k2)   =   zaxis(k) + Z_rand(k1:k2)*hz(k)
               vll_start(k1:k2) = vpart*P_rand(k1:k2) 
               mu_start(k1:k2)  = E_BEAMS(l) ! Total energy for now
               beam(k1:k2)      = l
               mass(k1:k2)      = MASS_BEAMS(l)
               charge(k1:k2)    = CHARGE_BEAMS(l)
               Zatom(k1:k2)     = 2
               t_end(k1:k2)     = t_end_in(1)
               weight(k1:k2)    = rateDT(i,j,k)/n3d(i,j,k)
               k1 = k2+1
            END DO
         END IF
         IF (lfusion_tritium) THEN
            ! Do the D-D -> T reaction
            l = l + 1
            vpart = sqrt(2*E_BEAMS(l)/mT)
            !E_NEUTRONS(l) = 3.02E6
            DO s = 1,nr1*nphi1*nz1
               i = MOD(s-1,nr1)+1
               j = MOD(s-1,nr1*nphi1)
               j = FLOOR(REAL(j) / REAL(nr1))+1
               k = CEILING(REAL(s) / REAL(nr1*nphi1))
               IF (n3d(i,j,k)==0) CYCLE
               k2 = n3d(i,j,k)+k1-1
               R_start(k1:k2)   =   raxis(i) + X_rand(k1:k2)*hr(i)
               phi_start(k1:k2) = phiaxis(j) + Y_rand(k1:k2)*hp(j)
               Z_start(k1:k2)   =   zaxis(k) + Z_rand(k1:k2)*hz(k)
               vll_start(k1:k2) = vpart*P_rand(k1:k2) 
               mu_start(k1:k2)  = E_BEAMS(l) ! Total energy for now
               beam(k1:k2)      = l
               mass(k1:k2)      = MASS_BEAMS(l)
               charge(k1:k2)    = CHARGE_BEAMS(l)
               Zatom(k1:k2)     = 1
               t_end(k1:k2)     = t_end_in(1)
               weight(k1:k2)    = rateDDT(i,j,k)/n3d(i,j,k)
               k1 = k2+1
            END DO
         END IF
         IF (lfusion_proton) THEN
            ! Do the D-D -> p reaction
            l = l + 1
            vpart = sqrt(2*E_BEAMS(l)/mp)
            !E_NEUTRONS(l) = 2.45E6
            DO s = 1,nr1*nphi1*nz1
               i = MOD(s-1,nr1)+1
               j = MOD(s-1,nr1*nphi1)
               j = FLOOR(REAL(j) / REAL(nr1))+1
               k = CEILING(REAL(s) / REAL(nr1*nphi1))
               IF (n3d(i,j,k)==0) CYCLE
               k2 = n3d(i,j,k)+k1-1
               R_start(k1:k2)   =   raxis(i) + X_rand(k1:k2)*hr(i)
               phi_start(k1:k2) = phiaxis(j) + Y_rand(k1:k2)*hp(j)
               Z_start(k1:k2)   =   zaxis(k) + Z_rand(k1:k2)*hz(k)
               vll_start(k1:k2) = vpart*P_rand(k1:k2) 
               mu_start(k1:k2)  = E_BEAMS(l) ! Total energy for now
               beam(k1:k2)      = l
               mass(k1:k2)      = MASS_BEAMS(l)
               charge(k1:k2)    = CHARGE_BEAMS(l)
               Zatom(k1:k2)     = 1
               t_end(k1:k2)     = t_end_in(1)
               weight(k1:k2)    = rateDDT(i,j,k)/n3d(i,j,k)
               k1 = k2+1
            END DO
         END IF
         IF (lfusion_He3) THEN
            ! Do the D-D -> He3 reaction
            l = l + 1
            vpart = sqrt(2*E_BEAMS(l)/mHe3)
            DO s = 1,nr1*nphi1*nz1
               i = MOD(s-1,nr1)+1
               j = MOD(s-1,nr1*nphi1)
               j = FLOOR(REAL(j) / REAL(nr1))+1
               k = CEILING(REAL(s) / REAL(nr1*nphi1))
               IF (n3d(i,j,k)==0) CYCLE
               k2 = n3d(i,j,k)+k1-1
               R_start(k1:k2)   =   raxis(i) + X_rand(k1:k2)*hr(i)
               phi_start(k1:k2) = phiaxis(j) + Y_rand(k1:k2)*hp(j)
               Z_start(k1:k2)   =   zaxis(k) + Z_rand(k1:k2)*hz(k)
               vll_start(k1:k2) = vpart*P_rand(k1:k2) 
               mu_start(k1:k2)  = E_BEAMS(l) ! Total energy for now
               beam(k1:k2)      = l
               mass(k1:k2)      = MASS_BEAMS(l)
               charge(k1:k2)    = CHARGE_BEAMS(l)
               Zatom(k1:k2)     = 2
               t_end(k1:k2)     = t_end_in(1)
               weight(k1:k2)    = rateDDHe(i,j,k)/n3d(i,j,k)
               k1 = k2+1
            END DO
         END IF
         IF (lfusion_DHe3) THEN
            ! Do the D-He3 -> He4 reaction
            l = l + 1
            vpart = sqrt(2*E_BEAMS(l)/mHe4)
            !E_NEUTRONS(l) = 2.45E6
            DO s = 1,nr1*nphi1*nz1
               i = MOD(s-1,nr1)+1
               j = MOD(s-1,nr1*nphi1)
               j = FLOOR(REAL(j) / REAL(nr1))+1
               k = CEILING(REAL(s) / REAL(nr1*nphi1))
               IF (n3d(i,j,k)==0) CYCLE
               k2 = n3d(i,j,k)+k1-1
               R_start(k1:k2)   =   raxis(i) + X_rand(k1:k2)*hr(i)
               phi_start(k1:k2) = phiaxis(j) + Y_rand(k1:k2)*hp(j)
               Z_start(k1:k2)   =   zaxis(k) + Z_rand(k1:k2)*hz(k)
               vll_start(k1:k2) = vpart*P_rand(k1:k2) 
               mu_start(k1:k2)  = E_BEAMS(l) ! Total energy for now
               beam(k1:k2)      = l
               mass(k1:k2)      = MASS_BEAMS(l)
               charge(k1:k2)    = CHARGE_BEAMS(l)
               Zatom(k1:k2)     = 2
               t_end(k1:k2)     = t_end_in(1)
               weight(k1:k2)    = rateDHe3(i,j,k)/n3d(i,j,k)
               k1 = k2+1
            END DO
            ! Do the D-He3 -> p reaction
            l = l + 1
            vpart = sqrt(2*E_BEAMS(l)/mp)
            !E_NEUTRONS(l) = 2.45E6
            DO s = 1,nr1*nphi1*nz1
               i = MOD(s-1,nr1)+1
               j = MOD(s-1,nr1*nphi1)
               j = FLOOR(REAL(j) / REAL(nr1))+1
               k = CEILING(REAL(s) / REAL(nr1*nphi1))
               IF (n3d(i,j,k)==0) CYCLE
               k2 = n3d(i,j,k)+k1-1
               R_start(k1:k2)   =   raxis(i) + X_rand(k1:k2)*hr(i)
               phi_start(k1:k2) = phiaxis(j) + Y_rand(k1:k2)*hp(j)
               Z_start(k1:k2)   =   zaxis(k) + Z_rand(k1:k2)*hz(k)
               vll_start(k1:k2) = vpart*P_rand(k1:k2) 
               mu_start(k1:k2)  = E_BEAMS(l) ! Total energy for now
               beam(k1:k2)      = l
               mass(k1:k2)      = MASS_BEAMS(l)
               charge(k1:k2)    = CHARGE_BEAMS(l)
               Zatom(k1:k2)     = 2
               t_end(k1:k2)     = t_end_in(1)
               weight(k1:k2)    = rateDHe3(i,j,k)/n3d(i,j,k)
               k1 = k2+1
            END DO
         END IF

         ! Deallocate and cleanup
         DEALLOCATE(X_rand,Y_rand,Z_rand,P_rand,n3d)
         partvmax = SQRT(MAXVAL(2*mu_start/mass)) ! Partvmax from max energy
         weight = weight*NINT(pi2/phiaxis(nphi)) ! NFP factor
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
      CALL MPI_BCAST(lgc2fo_start,nparticles,MPI_LOGICAL, master, MPI_COMM_BEAMS,ierr_mpi)
#endif
      CALL mpidealloc(l3d,win_l3d)
      CALL mpidealloc(rateDT,win_rateDT)
      CALL mpidealloc(rateDDT,win_rateDDT)
      CALL mpidealloc(rateDDHe,win_rateDDHe)
      CALL mpidealloc(rateDHe3,win_rateDHe3)


      ! Now we need to calculate the magnetic moment as mu is perp energy
      CALL mpialloc(B_help, nparticles, myid_sharmem, 0, MPI_COMM_SHARMEM, win_B_help)

      IF (mylocalid == mylocalmaster) B_help = 0
      
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
      
      ! Now calculate mu
      mu_start = (mu_start - 0.5*mass*vll_start*vll_start) / B_help

      CALL mpidealloc(B_help,win_B_help)

      RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_init_fusion
