!-----------------------------------------------------------------------
!     Module:        fieldlines_init_coil
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/01/2012
!     Description:   This subroutine reads the VMEC input file and the
!                    coils file and calcultes the vacuum field on the
!                    grid.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_init_coil
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_runtime
      USE fieldlines_grid, ONLY: raxis,phiaxis,zaxis, nr, nphi, nz, &
                                 rmin, rmax, zmin, zmax, phimin, &
                                 phimax, B_R, B_Z, B_PHI
      USE vmec_input,  ONLY: extcur_in => extcur, read_indata_namelist,&
                             nv_in => nzeta, nfp_in => nfp, nigroup
      USE biotsavart
      USE mpi_params
      USE mpi_inc            
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: mylocalid, mylocalmaster
      INTEGER :: MPI_COMM_LOCAL
      INTEGER :: ier, iunit, s, i, j, mystart, myend, k, ik, ig
      REAL(rprec)  :: br, bphi, bz, current, current_first, &
                      br_temp, bphi_temp, bz_temp
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      ! Divide up Work
      mylocalid = master
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_LOCAL, ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_LOCAL, mylocalid, ierr_mpi )              ! MPI
#endif
      mylocalmaster = master
      
      ! Read the input file for the EXTCUR array, NV, and NFP
      IF (mylocalid == mylocalmaster) THEN
         iunit = 11
         OPEN(UNIT=iunit, FILE='input.' // TRIM(id_string), STATUS='OLD', IOSTAT=ier)
         IF (ier /= 0) CALL handle_err(FILE_OPEN_ERR,id_string,ier)
         CALL read_indata_namelist(iunit,ier)
         IF (ier /= 0) CALL handle_err(VMEC_INPUT_ERR,id_string,ier)
         CLOSE(iunit)
      END IF
#if defined(MPI_OPT)
      CALL MPI_BCAST(extcur_in,nigroup,MPI_REAL, mylocalmaster, MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'fieldlines_init_coil',ierr_mpi)
#endif
      
      ! Read the coils file
      CALL parse_coils_file(TRIM(coil_string))
      nextcur = SIZE(coil_group) !SAL
      DO ik = 1, nextcur
         DO j = 1, coil_group(ik) % ncoil
            current = coil_group(ik) % coils(j) % current
            IF (j .eq. 1) current_first = current
            IF (lraw) THEN
               coil_group(ik) % coils(j) % current = current*extcur_in(ik)
            ELSE
               IF (current_first .ne. zero) coil_group(ik) % coils(j) % current = (current/current_first)*extcur_in(ik)
            END IF
         END DO
      END DO
      
      ! Reset the phi grid limit to match mgrid
      phimin = 0.0
      phimax = pi2 / nfp_bs
      FORALL(i = 1:nphi) phiaxis(i) = (i-1)*(phimax-phimin)/(nphi-1) + phimin
      
      IF (lverb) THEN
         WRITE(6,'(A)')   '----- COILS Information -----'
         WRITE(6,'(A,A)') '   FILE: ',TRIM(coil_string)
         WRITE(6,'(A,I3)')'   Coil Periodicity: ',nfp_bs
         WRITE(6,'(A,I3)')'   Current Systems:  ',nextcur
         IF (lraw) THEN
            WRITE(6,'(A)')   '   Current Type:      RAW'
         ELSE
            WRITE(6,'(A)')   '   Current Type:      SCALED'
         END IF
         DO ik = 1, nextcur
            IF (ABS(extcur_in(ik)).ge. 1.0E9) THEN
               WRITE(6,'(A,I4,A,F8.3,A)')'   Num Coils  = ',coil_group(ik)%ncoil,'  EXTCUR = ',extcur_in(ik)/1.0E9,' [GA]'
            ELSE IF (ABS(extcur_in(ik)).ge. 1.0E6) THEN
               WRITE(6,'(A,I4,A,F8.3,A)')'   Num Coils  = ',coil_group(ik)%ncoil,'  EXTCUR = ',extcur_in(ik)/1.0E6,' [MA]'
            ELSE IF (ABS(extcur_in(ik)).ge. 1.0E3) THEN
               WRITE(6,'(A,I4,A,F8.3,A)')'   Num Coils  = ',coil_group(ik)%ncoil,'  EXTCUR = ',extcur_in(ik)/1.0E3,' [kA]'
            ELSE
               WRITE(6,'(A,I4,A,F8.3,A)')'   Num Coils  = ',coil_group(ik)%ncoil,'  EXTCUR = ',extcur_in(ik),' [A]'
            END IF
         END DO
         WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Vacuum Field Calculation [',0,']%'
         CALL FLUSH(6)
      END IF

      ! Break up the Work
      CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL,1, nr*nphi*nz, mystart, myend)

      IF (lafield_only) THEN
         DO s = mystart, myend
            i = MOD(s-1,nr)+1
            j = MOD(s-1,nr*nphi)
            j = FLOOR(REAL(j) / REAL(nr))+1
            k = CEILING(REAL(s) / REAL(nr*nphi))
            br = 0.; bphi = 0.; bz = 0.
            br_temp   = 0.0;
            bphi_temp = 0.0;
            bz_temp   = 0.0;
            DO ig = 1, nextcur
               IF (extcur_in(ig) == 0) CYCLE
               CALL afield(raxis(i), phiaxis(j), zaxis(k), br, bphi, bz, IG = ig)
               br_temp = br_temp + br
               bphi_temp = bphi_temp + bphi
               bz_temp = bz_temp + bz
            END DO
            B_R(i,j,k)   = br_temp
            B_PHI(i,j,k) = bphi_temp
            B_Z(i,j,k)   = bz_temp
            IF (lverb .and. (MOD(s,nr) == 0)) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
               CALL FLUSH(6)
            END IF
         END DO
      ELSE
         DO s = mystart, myend
            i = MOD(s-1,nr)+1
            j = MOD(s-1,nr*nphi)
            j = FLOOR(REAL(j) / REAL(nr))+1
            k = CEILING(REAL(s) / REAL(nr*nphi))
            br = 0.; bphi = 0.; bz = 0.
            br_temp   = 0.0;
            bphi_temp = 0.0;
            bz_temp   = 0.0;
            DO ig = 1, nextcur
               IF (extcur_in(ig) == 0) CYCLE
               CALL bfield(raxis(i), phiaxis(j), zaxis(k), br, bphi, bz, IG = ig)
               br_temp = br_temp + br
               bphi_temp = bphi_temp + bphi
               bz_temp = bz_temp + bz
            END DO
            B_R(i,j,k)   = br_temp
            B_PHI(i,j,k) = bphi_temp
            B_Z(i,j,k)   = bz_temp
            IF (lverb .and. (MOD(s,nr) == 0)) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
               CALL FLUSH(6)
            END IF
         END DO
      END IF

      ! Clean up the progress bar
      IF (lverb) THEN
         CALL backspace_out(6,36)
         CALL FLUSH(6)
         WRITE(6,'(36X)',ADVANCE='no')
         CALL FLUSH(6)
         CALL backspace_out(6,36)
         CALL FLUSH(6)
      END IF    
      
      ! Free Variables
      CALL cleanup_biotsavart

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_coil1',ierr_mpi)
      CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'fieldlines_init_coil: MPI_COMM_LOCAL',ierr_mpi)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_coil',ierr_mpi)
#endif
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_init_coil
