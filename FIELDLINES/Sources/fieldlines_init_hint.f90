!-----------------------------------------------------------------------
!     Subroutine:    fieldlines_init_hint
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          07/27/2021
!     Description:   This subroutine handles setting up the magnetic
!                    field based on HINT equilibira.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_init_hint
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_runtime, extcur_fieldlines => extcur
      USE fieldlines_grid, ONLY: raxis_g => raxis, phiaxis, &
                                 zaxis_g => zaxis, nr, nphi, nz, &
                                 rmin, rmax, zmin, zmax, phimin, &
                                 phimax, vc_adapt_tol, B_R, B_Z, B_PHI,&
                                 BR_spl, BZ_spl
      USE read_hint_mod, ONLY: get_hint_grid, get_hint_B, get_hint_press, &
                               get_hint_maxp, read_hint_deallocate, &
                               get_hint_magaxis, get_hint_gridB
      USE mpi_params                                                    ! MPI
!      USE mpi
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: BYTE_8 = SELECTED_INT_KIND (8)
#if defined(MPI_OPT)
      INTEGER(KIND=BYTE_8),ALLOCATABLE :: mnum(:), moffsets(:)
      INTEGER :: numprocs_local, mylocalid, mylocalmaster
      INTEGER :: MPI_COMM_LOCAL
#endif
      INTEGER :: ier, s, i, j, k, mystart, myend
      REAL(rprec) :: br, bphi, bz
      INTEGER :: nrh,nzh,nph
      REAL(rprec) :: rmin_hint, rmax_hint, zmin_hint, zmax_hint, &
                     pmax_hint, pres_max, betatot
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Handle what we do

      ! Divide up Work
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_LOCAL, ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_LOCAL, mylocalid, ierr_mpi )              ! MPI
      CALL MPI_COMM_SIZE( MPI_COMM_LOCAL, numprocs_local, ierr_mpi )          ! MPI
#endif
      mylocalmaster = master

      ! Get the maximum pressure
      CALL get_hint_maxp(pres_max)

      ! Write info to screen
      IF (lverb) THEN
         betatot = 0
         CALL get_hint_grid(nrh,nzh,nph,rmin_hint,rmax_hint,zmin_hint,zmax_hint,pmax_hint)
         WRITE(6,'(A)')               '----- HINT Information -----'
         WRITE(6,'(A,F9.5,A,F9.5,A,I4)') '   R   = [',rmin_hint,',',rmax_hint,'];  NR:   ',nrh
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   PHI = [',0.0,',',pmax_hint,'];  NPHI: ',nph
         WRITE(6,'(A,F8.5,A,F8.5,A,I4)') '   Z   = [',zmin_hint,',',zmax_hint,'];  NZ:   ',nzh
         WRITE(6,'(A,F7.2,A)')           '   PRES_MAX = ',pres_max*1E-3,' [kPa]'
      END IF
      
      IF (lverb) THEN
         WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Plasma Field Calculation [',0,']%'
         CALL FLUSH(6)
      END IF

      ! Break up the Work
      CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL,1, nr*nphi*nz, mystart, myend)


      IF (lafield_only) THEN
      ELSE
         DO s = mystart, myend
            i = MOD(s-1,nr)+1
            j = MOD(s-1,nr*nphi)
            j = FLOOR(REAL(j) / REAL(nr))+1
            k = CEILING(REAL(s) / REAL(nr*nphi))

            ! Bfield
            CALL get_hint_gridB(i,j,k,br,bphi,bz)
            B_R(i,j,k)   = B_R(i,j,k)   + br
            B_PHI(i,j,k) = B_PHI(i,j,k) + bphi
            B_Z(i,j,k)   = B_Z(i,j,k)   + bz

            IF (lverb .and. (MOD(s,nr) == 0)) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
               CALL FLUSH(6)
            END IF
         END DO
      END IF
      
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
#endif
      
      ! Free variables
      CALL read_hint_deallocate
      
      IF (lverb) THEN
         CALL backspace_out(6,36)
         CALL FLUSH(6)
         WRITE(6,'(36X)',ADVANCE='no')
         CALL FLUSH(6)
         CALL backspace_out(6,36)
         CALL FLUSH(6)
      END IF    
      
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'beams3d_init_vmec',ierr_mpi)
#endif
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_init_hint
