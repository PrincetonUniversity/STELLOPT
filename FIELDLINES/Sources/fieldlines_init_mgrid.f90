!-----------------------------------------------------------------------
!     Module:        fieldlines_init_mgrid
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/21/2012
!     Description:   This subroutine reads the VMEC input file and the
!                    mgrid file, then composes the new grid.  Note that
!                    must redefine our phi grid in FIELDLINES to have
!                    the exact extent of the mgrid file.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_init_mgrid
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
      USE mgrid_field_mod, pi2_mgrid => pi2
      USE mpi_params                       
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: BYTE_8 = SELECTED_INT_KIND (8)
      INTEGER(KIND=BYTE_8),ALLOCATABLE :: mnum(:), moffsets(:)
      INTEGER :: numprocs_local, mylocalid, mylocalmaster
      INTEGER :: MPI_COMM_LOCAL
      INTEGER(KIND=BYTE_8) :: chunk
      INTEGER :: ier, iunit, s, i, j, mystart, myend, k
      REAL(rprec)  :: br, bphi, bz
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------

      ! Divide up Work
      numprocs_local = 1; mylocalid = master
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_LOCAL, ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_LOCAL, mylocalid, ierr_mpi )              ! MPI
      CALL MPI_COMM_SIZE( MPI_COMM_LOCAL, numprocs_local, ierr_mpi )          ! MPI
#endif
      mylocalmaster = master


      ! Read the input file for the EXTCUR array, NV, and NFP
      IF (.not. ALLOCATED(extcur) .and. lmgrid) THEN
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
         IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'beams3d_init_mgrid',ierr_mpi)
         CALL MPI_BCAST(nv_in,1,MPI_INTEGER, mylocalmaster, MPI_COMM_LOCAL,ierr_mpi)
         IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'beams3d_init_mgrid',ierr_mpi)
         CALL MPI_BCAST(nfp_in,1,MPI_INTEGER, mylocalmaster, MPI_COMM_LOCAL,ierr_mpi)
         IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'beams3d_init_mgrid',ierr_mpi)
#endif
         nextcur = nigroup
         ALLOCATE(extcur(nextcur),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'EXTCUR',ier)
         extcur = 0.0
         extcur(1:nextcur) = extcur_in(1:nextcur)
      END IF
      
      ! Read the mgrid file
      CALL mgrid_load(mgrid_string,extcur,nextcur,nv_in,nfp_in,ier,mylocalid,MPI_COMM_LOCAL)
      IF (lverb) THEN
         CALL mgrid_info(6)
         WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Vacuum Field Calculation [',0,']%'
         CALL FLUSH(6)
      END IF
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init_mgrid',ierr_mpi)
      
      ! Check for grid consistency
      IF ((rmin < rminb) .or. (rmax > rmaxb) .or. &
          (zmin < zminb) .or. (zmax > zmaxb)) THEN
            IF (lverb) THEN
               WRITE(6,'(A)') '!!!!!!!!!!!! WARNING !!!!!!!!!!!!'
               WRITE(6,'(A)') '!!  Desired Grid larger than   !!'
               WRITE(6,'(A)') '!!  mgrid extent.  Lower grid  !!'
               WRITE(6,'(A)') '!!  dimensions or use coils    !!'
               WRITE(6,'(A)') '!!  get vacuum field.          !!'
               WRITE(6,'(A)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            END IF
#if defined(MPI_OPT)
         CALL MPI_FINALIZE(ierr_mpi)
         IF (ierr_mpi /=0) CALL handle_err(MPI_FINE_ERR,'fieldlines_init_mgrid',ierr_mpi)
#endif
         stop
      END IF
      
      ! Reset the phi grid limit to match mgrid
      phimin = 0.0
      phimax = pi2 / nfp_in
      FORALL(i = 1:nphi) phiaxis(i) = (i-1)*(phimax-phimin)/(nphi-1) + phimin

      ! Break up the Work
      CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL,1, nr*nphi*nz, mystart, myend)

      IF (lafield_only) THEN
         STOP '!!!!!!!!!lafield_only requires coils file!!!!!!!!!'
      ELSE
         DO s = mystart, myend
            i = MOD(s-1,nr)+1
            j = MOD(s-1,nr*nphi)
            j = FLOOR(REAL(j) / REAL(nr))+1
            k = CEILING(REAL(s) / REAL(nr*nphi))
            CALL mgrid_bcyl(raxis(i), phiaxis(j), zaxis(k), &
                            br, bphi, bz, ier)
            IF (ier /= 0) stop 'mgrid_bcyl error!'
            B_R(i,j,k)   = br
            B_PHI(i,j,k) = bphi
            B_Z(i,j,k)   = bz
            IF (lverb .and. (MOD(s,nr) == 0)) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
               CALL FLUSH(6)
            END IF
         END DO
      END IF
      
      ! Free Variables (put this here to make room)
      CALL mgrid_free(ier,MPI_COMM_LOCAL)
      
      ! Clean up the progress bar
      IF (lverb) THEN
         CALL backspace_out(6,38)
         WRITE(6,'(38X)',ADVANCE='no')
         CALL backspace_out(6,38)
         CALL FLUSH(6)
      END IF    
      

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
      END SUBROUTINE fieldlines_init_mgrid
