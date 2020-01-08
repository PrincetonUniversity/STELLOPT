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
      USE beams3d_runtime, ONLY: id_string, npoinc, t_end, lbeam, &
                                 nbeams, beam, e_beams, charge_beams, &
                                 mass_beams, lverb, p_beams, MPI_BARRIER_ERR,&
                                 MPI_BCAST_ERR,nprocs_beams,handle_err, ldepo,&
                                 MPI_REDU_ERR
      USE safe_open_mod, ONLY: safe_open
      USE EZspline
!DEC$ IF DEFINED (MPI_OPT)
      USE mpi_params ! MPI
      USE mpi
!DEC$ ENDIF
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error Flag
!          iunit        File ID
!          ndist        Number of Vll divisions for dist function
!          ns           Number of flux divisions for current calculation
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, iunit, istat, i, j, k, ninj, sbeam, ebeam
      REAL(rprec) :: maxdist,mindist,v1,v2,dist,ddist,s1,s2, vp_temp
      LOGICAL, ALLOCATABLE     :: partmask(:), partmask2(:,:), partmask2t(:,:)
      INTEGER, ALLOCATABLE  :: int_mask(:), int_mask2(:,:), nlost(:)
      INTEGER, ALLOCATABLE  :: dist_func(:,:,:)
      REAL, ALLOCATABLE     :: real_mask(:)
      INTEGER, PARAMETER :: ndist = 100
!DEC$ IF DEFINED (MPI_OPT)
      INTEGER :: status(MPI_STATUS_size) !mpi stuff
      INTEGER :: mystart, mypace, sender
      INTEGER, ALLOCATABLE :: revcounts(:), displs(:)
      REAL(rprec), ALLOCATABLE :: buffer_mast(:,:), buffer_slav(:,:)
!DEC$ ENDIF
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      IF (lverb) WRITE(6,'(A)')  '----- BEAM DIAGNOSTICS -----'


      ! DEALLOCATE stuff we don't need
      IF (ALLOCATED(Z_lines)) DEALLOCATE(Z_lines)
      IF (ALLOCATED(moment_lines)) DEALLOCATE(moment_lines)
      IF (ALLOCATED(U_lines)) DEALLOCATE(U_lines)
      IF (ALLOCATED(B_lines)) DEALLOCATE(B_lines)


      IF (.not. lbeam) THEN
         beam = 1
         nbeams = 1
      END IF

      CALL FLUSH(6)

      mystart = mystart_save
      myend = myend_save

      ! Main Allocations
      IF (ALLOCATED(shine_through)) DEALLOCATE(shine_through)
      IF (ALLOCATED(nlost)) DEALLOCATE(nlost)
      ALLOCATE(shine_through(nbeams))
      ALLOCATE(nlost(nbeams))
      ALLOCATE(partmask(mystart:myend))
      ALLOCATE(partmask2(0:npoinc,mystart:myend))
      ALLOCATE(partmask2t(0:npoinc,mystart:myend))
      ALLOCATE(int_mask(mystart:myend))
      ALLOCATE(int_mask2(0:npoinc,mystart:myend))
      ALLOCATE(real_mask(mystart:myend))
      maxdist=MAXVAl(MAXVAL(vll_lines,DIM=2,MASK=(ABS(vll_lines)<1E8)),DIM=1)
      mindist=MINVAl(MINVAL(vll_lines,DIM=2,MASK=(ABS(vll_lines)<1E8)),DIM=1)
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS, ierr_mpi)
      IF (ierr_mpi /= 0) CALL handle_err(MPI_BARRIER_ERR, 'beams3d_follow', ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,maxdist,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /= 0) CALL handle_err(MPI_REDU_ERR, 'beams3d_diagnostic1', ierr_mpi)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,mindist,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /= 0) CALL handle_err(MPI_REDU_ERR, 'beams3d_diagnostic2', ierr_mpi)
!DEC$ ENDIF
      maxdist = max(ABS(mindist),ABS(maxdist))
      mindist = -maxdist

      ! Setup masking arrays
      FORALL(i=0:npoinc) int_mask2(i,mystart:myend) = beam(mystart:myend)              ! Index of beams
      WHERE(      ( (R_lines(:,mystart:myend)==0) .and. (PHI_lines(:,mystart:myend)==-1) )&
             .or. (neut_lines(:,mystart:myend)) ) int_mask2(:,mystart:myend) = 0  ! Mask the neutral and lost particles
      int_mask(mystart:myend) = 0
      IF (lbeam) int_mask(mystart:myend) = 3
      int_mask(mystart:myend)  = COUNT(neut_lines(:,mystart:myend),DIM=1)-1                  ! Starting index of every charge particle
      WHERE(int_mask < 0) int_mask = 0
      FORALL(j=mystart:myend) real_mask(j) = S_lines(int_mask(j),j) ! Starting points in s

      ! Don't need R_lines or PHI_lines after this point
      IF (ALLOCATED(R_lines)) DEALLOCATE(R_lines)
      IF (ALLOCATED(PHI_lines)) DEALLOCATE(PHI_lines)

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
         shine_through(i) = 100.*COUNT(neut_lines(1,mystart:myend) .and. (beam(mystart:myend)==i),DIM=1)/COUNT(beam==i)
         nlost(i)         = COUNT(lost_lines(mystart:myend).and. beam(mystart:myend) == i)
      END DO

!DEC$ IF DEFINED (MPI_OPT)
      IF (myworkid == master) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE, dist_func,     nbeams*ndist*(npoinc+1), MPI_INTEGER,          MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, shine_through, nbeams,                  MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, nlost,         nbeams,                  MPI_INTEGER,          MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
      ELSE
         CALL MPI_REDUCE(dist_func,     dist_func,     nbeams*ndist*(npoinc+1), MPI_INTEGER,          MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(shine_through, shine_through, nbeams,                  MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         CALL MPI_REDUCE(nlost,         nlost,         nbeams,                  MPI_INTEGER,          MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
      END IF
!DEC$ ENDIF

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
            ninj  = COUNT(beam .eq. i)
            WRITE(iunit,'(A)') ' Particles Launched  Particles Lost  Lost(%)  TIME_END'
            WRITE(iunit,'(6X,I10,11X,I5,7x,F5.1,6x,E22.12)') ninj, nlost(i), 100.*nlost(i)/ninj, MAXVAL(t_end)
            WRITE(iunit,'(A)') ' '
            CALL FLUSH(iunit)

            ! Screen Output
            IF (lverb) THEN
               IF (i==1) WRITE(6,'(A)')  ' BEAMLINE     ENERGY [keV]   CHARGE [e]   MASS [Mp]   Particles [#]   Lost [%]  Shinethrough [%]'
               WRITE(6,'(I5,3(11X,I3),8X,I8,2(8X,F5.1))') i,NINT(E_BEAMS(i)*6.24150636309E15),NINT(CHARGE_BEAMS(i)*6.24150636309E18),&
                                         NINT(MASS_BEAMS(i)*5.97863320194E26), ninj, 100.*nlost(i)/ninj, shine_through(i)
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
      ELSE
         DEALLOCATE(dist_func)
      END IF

      ! Don't need R_lines or PHI_lines after this point
      IF (ALLOCATED(neut_lines)) DEALLOCATE(neut_lines)

      ! BEAM DIAGNOSTICS
      IF (lbeam .and. .not.ldepo) THEN

         ! Birth Profiles and Power deposition and current profile
         IF (ALLOCATED(ndot_prof)) DEALLOCATE(ndot_prof)
         IF (ALLOCATED(epower_prof)) DEALLOCATE(epower_prof)
         IF (ALLOCATED(ipower_prof)) DEALLOCATE(ipower_prof)
         IF (ALLOCATED(j_prof)) DEALLOCATE(j_prof)
         ALLOCATE(ndot_prof(1:nbeams,1:ns_prof))
         ALLOCATE(epower_prof(1:nbeams,1:ns_prof))
         ALLOCATE(ipower_prof(1:nbeams,1:ns_prof))
         ALLOCATE(j_prof(1:nbeams,1:ns_prof))

         !Initialize arrays
         ndot_prof = 0
         epower_prof = 0
         ipower_prof = 0
         j_prof = 0

         ! Loop over beam lines
         DO k = 1, ns_prof
            s1 = REAL(k-1)/REAL(ns_prof)
            s2 = REAL(k)/REAL(ns_prof)
            partmask(mystart:myend)  = ((real_mask(mystart:myend) >= s1) .and. (real_mask(mystart:myend) < s2))
            partmask2(:,mystart:myend) = ((S_lines(:,mystart:myend) >= s1) .and. (S_lines(:,mystart:myend) < s2))
            CALL EZspline_interp(Vp_spl_s,s2,vp_temp,ier)
            !vp_temp = (s1+s2)*vp_temp ! because 2*s*vp is the qauntity we want and s=0.5*(s1+s2) ! do this if we want to use the rho grid.
            DO i = 1, nbeams
               partmask2t(:,mystart:myend)=(partmask2(:,mystart:myend).and.(int_mask2(:,mystart:myend)==i))
               ndot_prof(i,k) = COUNT((partmask(mystart:myend).and.(beam(mystart:myend)==i)))/vp_temp
               epower_prof(i,k) = SUM(SUM(PE_lines(:,mystart:myend),DIM=1,MASK=partmask2t(:,mystart:myend)))/vp_temp
               ipower_prof(i,k) = SUM(SUM(PI_lines(:,mystart:myend),DIM=1,MASK=partmask2t(:,mystart:myend)))/vp_temp
               real_mask(mystart:myend)=SUM(vll_lines(:,mystart:myend),MASK=partmask2t(:,mystart:myend),DIM=1)*(t_end-int_mask*dt_out)/(COUNT(partmask2t(:,mystart:myend),DIM=1)+1)
               j_prof(i,k) = SUM(real_mask)
            END DO
         END DO

!DEC$ IF DEFINED (MPI_OPT)
         IF (myworkid == master) THEN
            CALL MPI_REDUCE(MPI_IN_PLACE, ndot_prof,   nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE, epower_prof, nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE, ipower_prof, nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
            CALL MPI_REDUCE(MPI_IN_PLACE, j_prof,      nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         ELSE
            CALL MPI_REDUCE(ndot_prof,   ndot_prof,   nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
            CALL MPI_REDUCE(epower_prof, epower_prof, nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
            CALL MPI_REDUCE(ipower_prof, ipower_prof, nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
            CALL MPI_REDUCE(j_prof,      j_prof,      nbeams*ns_prof, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_BEAMS, ierr_mpi)
         END IF
!DEC$ ENDIF

         ! Was only needed if no weight specified
         IF (myworkid == master) THEN
         !   DO i = 1, nbeams
         !      ninj  = COUNT(beam .eq. i)
         !      ndot_prof(i,:) = (P_beams(i)/E_beams(i))*REAL(ndot_prof(i,:)) / REAL(ninj)
         !      epower_prof(i,:) = (P_beams(i)/E_beams(i))*epower_prof(i,:)/ REAL(ninj)
         !      ipower_prof(i,:) = (P_beams(i)/E_beams(i))*ipower_prof(i,:)/ REAL(ninj)
         !   END DO
            DEALLOCATE(dist_func)
         END IF

      END IF

      ! Give master a full copy of lost_lines (for STELLOPT interface)
      partmask(mystart:myend) = lost_lines(mystart:myend)
      IF (myworkid == master) THEN
           DEALLOCATE(lost_lines)
           ALLOCATE(lost_lines(1:nparticles),revcounts(nprocs_beams),displs(nprocs_beams))
           lost_lines(mystart:myend) = partmask(mystart:myend)
      ELSE
           ALLOCATE(revcounts(nprocs_beams),displs(nprocs_beams))
      END IF
      CALL MPI_GATHER(mystart,1,MPI_INTEGER,displs,1,MPI_INTEGER,master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_GATHER(myend-mystart+1,1,MPI_INTEGER,revcounts,1,MPI_INTEGER,master, MPI_COMM_BEAMS, ierr_mpi)
      CALL MPI_GATHERV(partmask(mystart:myend),myend-mystart+1,MPI_LOGICAL,&
                       lost_lines,revcounts,displs,MPI_LOGICAL, master, MPI_COMM_BEAMS, ierr_mpi)
      IF (ALLOCATED(revcounts)) DEALLOCATE(revcounts)
      IF (ALLOCATED(displs)) DEALLOCATE(displs)
      IF (myworkid /= master) DEALLOCATE(lost_lines)

      CALL beams3d_write('DIAG')

      DEALLOCATE(int_mask2,int_mask)
      DEALLOCATE(partmask2,partmask2t)
      DEALLOCATE(partmask,real_mask)


!-----------------------------------------------------------------------
!     End Subroutine
!----------------------------------------------------------------------- 
      END SUBROUTINE beams3d_diagnostics
