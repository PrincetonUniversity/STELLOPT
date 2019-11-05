!-----------------------------------------------------------------------
!     Subroutine:    stellopt_fieldlines
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/05/2019
!     Description:   This subroutine provides and interface to the
!                    FIELDLINES code.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_fieldlines(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime, ONLY: proc_string
      USE fieldlines_grid
      USE fieldlines_lines, ONLY: R_lines, Z_lines, PHI_lines, &
                                  Rhc_lines, Zhc_lines, B_lines
      USE fieldlines_runtime, id_string_fieldlines => id_string, pi_fieldlines => pi, &
                              pi2_fieldlines => pi2, mu0_fieldlines =>mu0, &
                              lverb_fieldlines => lverb
      USE EZspline_obj
      USE EZspline
      USE mpi_params
      USE mpi_inc
      USE mpi_sharmem
      USE wall_mod, ONLY: wall_free
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!        iflag         Error flag
!----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!----------------------------------------------------------------------
      INTEGER ::  ier, nshar
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      myworkid = master
      pi_fieldlines = 4.0 * ATAN(1.0)
      pi2_fieldlines = 8.0 * ATAN(1.0)
      mu0_fieldlines = 16.0E-7 * ATAN(1.0)
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_MYWORLD, MPI_COMM_FIELDLINES, ierr_mpi)
      CALL MPI_COMM_RANK(MPI_COMM_FIELDLINES, myworkid, ierr_mpi)
      CALL MPI_COMM_SIZE(MPI_COMM_FIELDLINES, nprocs_fieldlines, ierr_mpi)
      CALL MPI_COMM_SPLIT_TYPE(MPI_COMM_FIELDLINES, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, MPI_COMM_SHARMEM, ierr_mpi)
      CALL MPI_COMM_RANK(MPI_COMM_SHARMEM, myid_sharmem, ierr_mpi)
      CALL MPI_COMM_SIZE(MPI_COMM_SHARMEM, nshar, ierr_mpi)
#endif

      IF (iflag < 0) RETURN
      lverb_fieldlines = lscreen
      id_string_fieldlines = TRIM(proc_string)
      lvmec    = .false.
      lpies    = .false.
      lspec    = .false.
      lcoil    = .false.
      lmgrid   = .false.
      lmu      = .false.
      lvessel  = .false.
      lvac     = .false.
      lrestart = .false.
      laxis_i  = .false.
      ladvanced = .false.
      lemc3 = .false.
      lerror_field = .false.
      lplasma_only = .false.
      lbfield_only = .false.
      lafield_only = .false.
      lreverse  = .false.
      lhitonly  = .false.
      lraw   = .false.
      lwall_trans = .false.
      ledge_start = .false.
      lnescoil    = .false.
      lmodb       = .false.
      nruntype = runtype_old
      coil_string   = ''
      mgrid_string  = ''
      vessel_string = ''
      restart_string = ''

      IF (lverb_fieldlines) THEN
         WRITE(6,'(a,f5.2)') 'FIELDLINES Version ',FIELDLINES_VERSION
         WRITE(6,'(A)')      '-----  MPI Parameters  -----'
         WRITE(6,'(A,I8)')  '   Nproc_total:  ', nprocs_fieldlines
         WRITE(6,'(A,3X,I5)')  '   Nproc_shared: ', nshar
         CALL FLUSH(6)
      END IF

      CALL fieldlines_init

      IF (lemc3 .or. lbfield_only .or. lafield_only) nruntype=runtype_norun

      SELECT CASE(nruntype)
         CASE(runtype_old)
            CALL fieldlines_follow
         CASE(runtype_full)
            IF (lverb_fieldlines) THEN
               WRITE(6,'(A)') '===========ROUGH GRID=========='
               WRITE(6,'(A,F8.5,A,F8.5,A)') '   EDGE_INNER[R,Z]   = [',&
                     r_start(1),',',z_start(1),']'
               WRITE(6,'(A,F8.5,A,F8.5,A)') '   EDGE_OUTER[R,Z]   = [',&
                     MAXVAL(r_start,MASK = r_start > 0),',',MAXVAL(z_start,MASK = r_start > 0),']'
            END IF
#if defined(MPI_OPT)
            CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
            IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init',ierr_mpi)
#endif
            CALL fieldlines_follow  ! This call on field grid
#if defined(MPI_OPT)
            CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
            IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init',ierr_mpi)
#endif
            CALL fieldlines_init_subgrid
            CALL fieldlines_follow  ! This call on subgrid grid
            !CALL fieldlines_periodic_orbits  ! This call on subgrid grid
         CASE(runtype_norun)
      END SELECT

      ! Output Data
      CALL fieldlines_write


      ! Clean up
      IF (lvessel) CALL wall_free(ier)
      IF (ASSOCIATED(raxis)) CALL mpidealloc(raxis,win_raxis)
      IF (ASSOCIATED(phiaxis)) CALL mpidealloc(phiaxis,win_phiaxis)
      IF (ASSOCIATED(zaxis)) CALL mpidealloc(zaxis,win_zaxis)
      IF (ASSOCIATED(B_R)) CALL mpidealloc(B_R,win_B_R)
      IF (ASSOCIATED(B_PHI)) CALL mpidealloc(B_PHI,win_B_PHI)
      IF (ASSOCIATED(B_Z)) CALL mpidealloc(B_Z,win_B_Z)
      IF (ASSOCIATED(BR4D)) CALL mpidealloc(BR4D,win_BR4D)
      IF (ASSOCIATED(BZ4D)) CALL mpidealloc(BZ4D,win_BZ4D)
      IF (ASSOCIATED(MU4D)) CALL mpidealloc(MU4D,win_MU4D)
      IF (ASSOCIATED(MODB4D)) CALL mpidealloc(MODB4D,win_MODB4D)
      IF (ALLOCATED(R_lines)) DEALLOCATE(R_lines)
      IF (ALLOCATED(Z_lines)) DEALLOCATE(Z_lines)
      IF (ALLOCATED(PHI_lines)) DEALLOCATE(PHI_lines)
      IF (ALLOCATED(Rhc_lines)) DEALLOCATE(Rhc_lines)
      IF (ALLOCATED(Zhc_lines)) DEALLOCATE(Zhc_lines)
      IF (ALLOCATED(B_lines)) DEALLOCATE(B_lines)
      IF (EZspline_allocated(BR_spl)) CALL EZspline_free(BR_spl,ier)
      IF (EZspline_allocated(BZ_spl)) CALL EZspline_free(BZ_spl,ier)
      IF (EZspline_allocated(MU_spl)) CALL EZspline_free(MU_spl,ier)
      IF (EZspline_allocated(MODB_spl)) CALL EZspline_free(MODB_spl,ier)

      CALL MPI_COMM_FREE(MPI_COMM_SHARMEM,ierr_mpi)
      CALL MPI_COMM_FREE(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (lverb_fieldlines) THEN
         write(6,*)'==========  FIELDLINES Complete  =========='
         write(6,*)'==========================================='
         CALL FLUSH(6)
      END IF
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_fieldlines
