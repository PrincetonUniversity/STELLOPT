!-----------------------------------------------------------------------
!     Module:        beams3d_diagnostics
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          01/09/2014
!     Description:   This subroutine outputs a diagnostic text file
!                    of the run.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_diagnostics
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_lines
      USE beams3d_grid, ONLY: nr, nphi, nz, B_R, B_PHI, B_Z, raxis, &
                                 zaxis, phiaxis,vp_spl_s
      USE beams3d_runtime, ONLY: id_string, npoinc, t_end, lbeam, lvac,&
                                 lvmec, charge_beams, &
                                 nbeams, beam, e_beams, charge_beams, &
                                 mass_beams, lverb, p_beams, MPI_BARRIER_ERR,&
                                 MPI_BCAST_ERR,nprocs_beams,handle_err, ldepo,&
                                 MPI_REDU_ERR, pi2, weight
      USE safe_open_mod, ONLY: safe_open
      USE EZspline
      USE mpi_params ! MPI
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error Flag
!          iunit        File ID
!          ndist        Number of Vll divisions for dist function
!          ns           Number of flux divisions for current calculation
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, iunit, istat, i, j, k, nhalf, sbeam, ebeam, ninj2
      REAL(rprec) :: maxdist,mindist,v1,v2,dist,ddist,s1,s2, vp_temp, dvll, dvperp, ninj
      LOGICAL, ALLOCATABLE     :: partmask(:), partmask2(:,:), partmask2t(:,:)
      INTEGER, ALLOCATABLE  :: int_mask(:), int_mask2(:,:)
      INTEGER, ALLOCATABLE  :: dist_func(:,:,:)
      REAL, ALLOCATABLE     :: real_mask(:),vllaxis(:),vperpaxis(:), nlost(:)
      REAL, ALLOCATABLE     :: help3d(:,:,:)
      INTEGER, PARAMETER :: ndist = 100
#if defined(MPI_OPT)
      INTEGER :: mystart, mypace
      REAL(rprec), ALLOCATABLE :: buffer_mast(:,:), buffer_slav(:,:)
#endif
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      IF (lverb) WRITE(6,'(A)')  '----- BEAM DIAGNOSTICS -----'


      ! DEALLOCATE stuff we do not need
      IF (ALLOCATED(Z_lines)) DEALLOCATE(Z_lines)
      IF (ALLOCATED(moment_lines)) DEALLOCATE(moment_lines)
      IF (ALLOCATED(U_lines)) DEALLOCATE(U_lines)
      IF (ALLOCATED(B_lines)) DEALLOCATE(B_lines)

      CALL FLUSH(6)

      mystart = mystart_save
      myend = myend_save

      ! Main Allocations
      IF (ALLOCATED(shine_through)) DEALLOCATE(shine_through)
      IF (ALLOCATED(shine_port)) DEALLOCATE(shine_port)
      IF (ALLOCATED(nlost)) DEALLOCATE(nlost)
      ALLOCATE(shine_through(nbeams))
      ALLOCATE(shine_port(nbeams))
      ALLOCATE(nlost(nbeams))
      ALLOCATE(partmask(mystart:myend))
      ALLOCATE(partmask2(0:npoinc,mystart:myend))
      ALLOCATE(partmask2t(0:npoinc,mystart:myend))
      ALLOCATE(int_mask(mystart:myend))
      ALLOCATE(int_mask2(0:npoinc,mystart:myend))
      ALLOCATE(real_mask(mystart:myend))
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= 0) CALL handle_err(MPI_BARRIER_ERR, 'beams3d_follow', ierr_mpi)
#endif
      maxdist = partvmax
      mindist = -partvmax

      ! Setup masking arrays
      FORALL(i=0:npoinc) int_mask2(i,mystart:myend) = beam(mystart:myend)              ! Index of beams
      WHERE(      ( (R_lines(:,mystart:myend)==0.0) .and. (PHI_lines(:,mystart:myend)==-1.0) )&
             .or. (neut_lines(:,mystart:myend)) ) int_mask2(:,mystart:myend) = 0  ! Mask the neutral and lost particles
      int_mask(mystart:myend) = 0
      IF (lbeam) int_mask(mystart:myend) = 3
      int_mask(mystart:myend)  = COUNT(neut_lines(:,mystart:myend),DIM=1)-1                  ! Starting index of every charge particle
      WHERE(int_mask < 0) int_mask = 0
      FORALL(j=mystart:myend) real_mask(j) = S_lines(int_mask(j),j) ! Starting points in s

      ! Do not need R_lines or PHI_lines after this point
      IF (ALLOCATED(R_lines)) DEALLOCATE(R_lines)
      IF (ALLOCATED(PHI_lines)) DEALLOCATE(PHI_lines)
      IF (ALLOCATED(neut_lines)) DEALLOCATE(neut_lines)

      ! Calculate distribution function
      ALLOCATE(dist_func(1:nbeams,1:ndist,0:npoinc))
      dist_func = 0
      dist = maxdist-mindist
      ddist = dist/ndist
      sbeam = MINVAL(beam(mystart:myend), DIM=1)
      ebeam = MAXVAL(beam(mystart:myend), DIM=1)
      DO k = 1, ndist
         v1 = mindist+(k-1)*ddist
         v2 = mindist+(k)*ddist
         partmask2(:,mystart:myend) = ((vll_lines(:,mystart:myend).ge.v1).and.(vll_lines(:,mystart:myend).lt.v2))
         DO i = sbeam, ebeam
            dist_func(i,k,0:npoinc) = COUNT(partmask2(:,mystart:myend).and.(int_mask2(:,mystart:myend)==i),DIM=2)
         END DO
      END DO

      ! Calculate shinethrough and loss
      shine_through = 0
      DO i = 1, nbeams
         nlost(i)          =      SUM(weight(mystart:myend), MASK = (end_state(mystart:myend) == 2 .and. (beam(mystart:myend)==i)))
         shine_through(i)  = 100.*SUM(weight(mystart:myend), MASK = (end_state(mystart:myend) == 3 .and. (beam(mystart:myend)==i)))/SUM(weight,MASK=(beam==i))
         shine_port(i)     = 100.*SUM(weight(mystart:myend), MASK = (end_state(mystart:myend) == 4 .and. (beam(mystart:myend)==i)))/SUM(weight,MASK=(beam==i))
         !shine_through(i) = 100.*COUNT(end_state(mystart:myend) == 3 .and. (beam(mystart:myend)==i),DIM=1)/COUNT(beam==i)
         !shine_port(i) = 100.*COUNT(end_state(mystart:myend) == 4 .and. (beam(mystart:myend)==i),DIM=1)/COUNT(beam==i)
         !nlost(i)         = COUNT(end_state(mystart:myend) == 2 .and. beam(mystart:myend) == i)
      END DO

#if defined(MPI_OPT)
      IF (myworkid == master) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE, dist_func,     nbeams*ndist*(npoinc+1), MPI_INTEGER,          MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, shine_through, nbeams,                  MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, shine_port,    nbeams,                  MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, nlost,         nbeams,                  MPI_REAL,          MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
      ELSE
         CALL MPI_REDUCE(dist_func,     dist_func,     nbeams*ndist*(npoinc+1), MPI_INTEGER,          MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(shine_through, shine_through, nbeams,                  MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(shine_port,    shine_port,    nbeams,                  MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(nlost,         nlost,         nbeams,                  MPI_REAL,          MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
      END IF
#endif

      IF (myworkid == master) THEN
         ! Open the file
         iunit = 10
         CALL safe_open(iunit,istat,'beams3d_diag_'//TRIM(id_string)//'.txt','replace','formatted')
         ! Output number of beams
         WRITE(iunit,'(A,I5)') 'BEAMLINES: ',nbeams
         ! Screen Output
         DO i = 1, nbeams
            ! Output beam information
            WRITE(iunit,'(A)') ' BEAMLINE        ENERGY                 CHARGE                 MASS'
            WRITE(iunit,'((I5,3(2X,E22.12)))') i,E_BEAMS(i),CHARGE_BEAMS(i),MASS_BEAMS(i)

            ! Output beam losses
            ninj  = SUM(weight,MASK=(beam==i))
            ninj2  = COUNT(beam==i)
            WRITE(iunit,'(A)') ' Particles Launched  Particles Lost  Lost(%)  TIME_END'
            WRITE(iunit,'(6X,I10,11X,I5,7x,F5.1,6x,E22.12)') ninj2, NINT(ninj2*nlost(i)/ninj), 100.*nlost(i)/ninj, MAXVAL(t_end)
            WRITE(iunit,'(A)') ' '
            CALL FLUSH(iunit)

            ! Screen Output
            IF (lverb) THEN
               IF (i==1) WRITE(6,'(A)')  ' BEAMLINE     ENERGY [keV]   CHARGE [e]   MASS [Mp]   Particles [#]   Lost [%]  Shinethrough [%]  Port [%]'
               WRITE(6,'(I5,3(9X,I5),8X,I8,3(8X,F5.1))') i,NINT(E_BEAMS(i)*6.24150636309E15),NINT(CHARGE_BEAMS(i)*6.24150636309E18),&
                                         NINT(MASS_BEAMS(i)*5.97863320194E26), ninj2, 100.*nlost(i)/ninj, shine_through(i), shine_port(i)
               IF (i==nbeams) WRITE(6,'(A,EN12.3)')     '# of Fast Ions in Distribution:: ',SUM(dist5d_prof)           
               CALL FLUSH(6)
            END IF
            ! Write Distribution Function
            WRITE(iunit,'(A)') ' Parallel Velocity Distribution'
            WRITE(iunit,'(A,100(1X,E22.12))') 'VLL',(mindist+(j+0.5)*ddist,j=0,ndist-1)
            WRITE(iunit,'(A)') '==========================================='
            DO j = 0, npoinc
               WRITE(iunit,'(100(1X,I12))') dist_func(i,1:ndist,j)
            END DO
            WRITE(iunit,'(A)') '==========================================='
            WRITE(iunit,'(A)') ' '
            CALL FLUSH(iunit)
         END DO
         CLOSE(iunit)
      END IF

      DEALLOCATE(dist_func)
      DEALLOCATE(int_mask2,int_mask)
      DEALLOCATE(partmask2,partmask2t)
      DEALLOCATE(partmask,real_mask)

      ! These diagnostics need Vp to be defined
      IF ((.not.ldepo .or. lrestart_grid).and. myworkid == master) THEN
         ! Allocate the parallel and perpendicular velcoity axis
         nhalf = ns_prof4/2
         ALLOCATE(dense_prof(nbeams,ns_prof1),j_prof(nbeams,ns_prof1))
         ALLOCATE(vllaxis(ns_prof4),vperpaxis(ns_prof5))
         ALLOCATE(help3d(nbeams,ns_prof1,ns_prof4))
         FORALL(k = 1:ns_prof4) vllaxis(k) = partvmax*REAL(k-nhalf-0.5)/REAL(nhalf)
         FORALL(k = 1:ns_prof5) vperpaxis(k) = partvmax*REAL(k-0.5)/REAL(ns_prof5)
         ! The DIST5D distribution is not normalized to anything at this point.
         !    We calculate a radial profile of FI density.
         help3d = SUM(SUM(SUM(dist5d_prof,DIM=6),DIM=4),DIM=3)
         dense_prof = SUM(help3d,DIM=3)
         ! Now calculate J_fast
         j_prof = 0
         DO k = 1, ns_prof4
            j_prof = j_prof + help3d(:,:,k)*vllaxis(k)
         END DO
         DO k = 1, ns_prof1
            j_prof(1:nbeams,k) = j_prof(1:nbeams,k)*charge_beams(1:nbeams) ! [A*m]
         END DO
         DEALLOCATE(help3d)
         !dense_prof = SUM(SUM(SUM(SUM(dist5d_prof,DIM=6),DIM=5),DIM=4),DIM=3)
         ! We not apply the volume element for the radial profiles [m^-3]
         DO k = 1, ns_prof1
            s1 = REAL(k-0.5)/REAL(ns_prof1) ! Rho
            s2 = s1*s1
            CALL EZspline_interp(Vp_spl_s,s2,vp_temp,ier)
            vp_temp = vp_temp*2*s1*(1./REAL(ns_prof1))
            epower_prof(:,k) = epower_prof(:,k)/vp_temp
            ipower_prof(:,k) = ipower_prof(:,k)/vp_temp
            ndot_prof(:,k)   =   ndot_prof(:,k)/vp_temp
            dense_prof(:,k)  =  dense_prof(:,k)/vp_temp
            j_prof(:,k)      =      j_prof(:,k)/vp_temp ! [A/m^2]
         END DO
         ! Normalize to velocity space volume element
         ! dvll = partvmax*2/ns_prof4 ! dVll
         ! dvperp = pi2*partvmax/ns_prof5 ! dVperp
         ! DO k = 1, ns_prof5 ! VPERP
         !    !s2 = REAL(k-0.5)/REAL(ns_prof5) ! Vperp_frac
         !    vp_temp = vperpaxis(k)*dvll*dvperp
         !    dist5d_prof(:,:,:,:,:,k) = dist5d_prof(:,:,:,:,:,k)/vp_temp
         ! END DO
         ! DEALLOCATIONS
         DEALLOCATE(vperpaxis,vllaxis)
         CALL beams3d_distnorm
      END IF

      CALL beams3d_write('DIAG')


!-----------------------------------------------------------------------
!     End Subroutine
!----------------------------------------------------------------------- 
      END SUBROUTINE beams3d_diagnostics
