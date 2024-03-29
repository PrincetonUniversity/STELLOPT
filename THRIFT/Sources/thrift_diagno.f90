!-----------------------------------------------------------------------
!     Subroutine:    thrift_diagno
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/07/2023
!     Description:   This subroutine handles calculating magnetic
!                    diagnostics.
!-----------------------------------------------------------------------
      SUBROUTINE thrift_diagno(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USe thrift_runtime
      USE diagno_input_mod, ONLY: bfield_points_file, bprobes_file,&
                                  mirnov_file, seg_rog_file, &
                                  flux_diag_file, BCAST_DIAGNO_INPUT
      USE diagno_runtime, ONLY: lverb_diagno => lverb, &
                                id_string_diagno => id_string, &
                                lcoil_diagno => lcoil, &
                                diagno_coil_string => coil_string,&
                                DIAGNO_VERSION, nprocs_diagno, lvac
      USE biotsavart, ONLY: cleanup_biotsavart
      USE virtual_casing_mod, ONLY: free_virtual_casing, virtual_casing_surf_dump
      USE mpi_params
      USE mpi_inc
      
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
!        iunit       File unit number
!----------------------------------------------------------------------
      INTEGER ::  ier, ik
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      myworkid = master
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_MYWORLD, MPI_COMM_DIAGNO, ierr_mpi)
      CALL MPI_COMM_RANK(MPI_COMM_DIAGNO, myworkid, ierr_mpi)
      CALL MPI_COMM_SIZE(MPI_COMM_DIAGNO, nprocs_diagno, ierr_mpi)
      CALL MPI_COMM_SPLIT_TYPE(MPI_COMM_DIAGNO, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, MPI_COMM_SHARMEM, ierr_mpi)
      CALL MPI_COMM_RANK(MPI_COMM_SHARMEM, myid_sharmem, ierr_mpi)
#endif

      IF (iflag < 0) RETURN
      lverb_diagno = lscreen
      id_string_diagno = TRIM(proc_string)
      lcoil_diagno = .false.
      lvac = .false.
      IF (LEN_TRIM(magdiag_coil) > 1) THEN
         lcoil_diagno = .true.
         diagno_coil_string = TRIM(magdiag_coil)
         INQUIRE(FILE=TRIM(diagno_coil_string),EXIST=lcoil_diagno,IOSTAT=iflag)
         IF (iflag /= 0) RETURN
      END IF
      CALL BCAST_DIAGNO_INPUT(master,MPI_COMM_DIAGNO,ierr_mpi)
      IF (lverb_diagno) write(6,'(A)')'==========================================='
      IF (lverb_diagno) write(6,'(A,F5.2,A)')'=========  D I A G N O  (v.',DIAGNO_VERSION,')  ========='
      IF (lverb_diagno) write(6,'(A)')' - Loading equilibrium surface'
      IF (lverb_diagno) CALL FLUSH(6)
      IF (lvmec) CALL diagno_init_vmec
      IF (.false.) CALL virtual_casing_surf_dump(327)
      ! All codes interface the same way
      IF (lcoil_diagno) CALL diagno_init_coil
      IF (lverb_diagno) write(6,*)' - Calculating diagnostic responses'
      IF (LEN_TRIM(bfield_points_file) > 1) CALL diagno_bfield
      IF (LEN_TRIM(bprobes_file) > 1) CALL diagno_bprobes
      IF (LEN_TRIM(mirnov_file) > 1) CALL diagno_mirnov
      IF (LEN_TRIM(seg_rog_file) > 1) CALL diagno_rogowski_new
      IF (LEN_TRIM(flux_diag_file) > 1) CALL diagno_flux
      ! Clean up
      IF (lcoil_diagno) CALL cleanup_biotsavart
      CALL free_virtual_casing(MPI_COMM_SHARMEM)
      CALL MPI_COMM_FREE(MPI_COMM_SHARMEM,ierr_mpi)
      CALL MPI_COMM_FREE(MPI_COMM_DIAGNO,ierr_mpi)
      IF (lverb_diagno) write(6,*)'============  DIAGNO Complete  ============'
      IF (lverb_diagno) write(6,*)'==========================================='
      IF (lverb_diagno) CALL FLUSH(6)
      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_diagno
