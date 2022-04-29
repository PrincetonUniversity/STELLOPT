!-----------------------------------------------------------------------
!     Module:        fieldlines_init_backflow
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          04/27/2021
!     Description:   This subroutine sets up running the fieldlines in
!                    reverse after wall hits. 
!                    First the code allocates twice the number of
!                    fieldlines it needs, but only uses the first half
!                    We then store those results, overwrite them with
!                    the reverse trace and then sort everything out
!                    after calling follow a second time.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_init_backflow
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_runtime
      USE fieldlines_lines
      USE wall_mod, ONLY: ihit_array, nface
      USE mpi_sharmem
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: MPI_COMM_LOCAL
      INTEGER :: mystart, mynewstart, mynewend, i, &
                 win_R_local, win_Z_local, win_phis_local, &
                 win_phie_local, win_L_local

      REAL(rprec), POINTER :: R_local(:), Z_local(:), phis_local(:), &
                              phie_local(:), L_local(:)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! If not lhitonly get out of here
      IF (.not. lhitonly) THEN
         IF (lverb) WRITE(6,'(A)') '!!!!! Backtrace only with lhitonly!!!!!!!'
         RETURN
      END IF

      ! Divide up the work
      !CALL MPI_CALC_MYRANGE(MPI_COMM_FIELDLINES,1, nlines, mystart, myend)
      mystart = LBOUND(R_lines,1)
      myend   = UBOUND(R_lines,1)

      ! Allocate the local variables
      CALL mpialloc(R_local, nlines*2, myid_sharmem, 0, MPI_COMM_SHARMEM, win_R_local)
      CALL mpialloc(Z_local, nlines*2, myid_sharmem, 0, MPI_COMM_SHARMEM, win_Z_local)
      CALL mpialloc(phis_local, nlines*2, myid_sharmem, 0, MPI_COMM_SHARMEM, win_phis_local)
      CALL mpialloc(phie_local, nlines*2, myid_sharmem, 0, MPI_COMM_SHARMEM, win_phie_local)
      CALL mpialloc(L_local, nlines, myid_sharmem, 0, MPI_COMM_SHARMEM, win_L_local)

      ! Default to zero
      IF (myid_sharmem==master) THEN
         R_local = 0; Z_local=0; phis_local = 0; phie_local=0; L_local = 0
      END IF

      ! Overwrite and start a new run (note _lines will be deallocated)
      L_local(mystart:myend) = L_lines(mystart:myend)
      IF (lhitonly) THEN
         ! These are the reruns of the hit
         R_local(mystart:myend)    =  R_lines(mystart:myend,0)
         Z_local(mystart:myend)    =  Z_lines(mystart:myend,0)
         phis_local(mystart:myend) =  PHI_lines(mystart:myend,0)
         phie_local(mystart:myend) =  PHI_lines(mystart:myend,2)
         ! Now load the back trace
         ! We use the zero point to avoid hitting the wall
         mynewstart = mystart + nlines
         mynewend   = myend   + nlines
         R_local(mynewstart:mynewend)      =  R_lines(mystart:myend,0)
         Z_local(mynewstart:mynewend)      =  R_lines(mystart:myend,0)
         phis_local(mynewstart:mynewend)   =  PHI_lines(mystart:myend,0)
         phie_local(mynewstart:mynewend)   =  PHI_lines(mystart:myend,0)-PHI_end(mystart:myend)
      ELSE
         PRINT *,'lhitonly must be active'
      END IF

      ! Now we exchange data
      i = MPI_UNDEFINED
      IF (myid_sharmem == master) i = 0
      CALL MPI_COMM_SPLIT( MPI_COMM_FIELDLINES,i,myworkid,MPI_COMM_LOCAL,ierr_mpi)
      IF (myid_sharmem == master) THEN
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,R_local,nlines*2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,Z_local,nlines*2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,phis_local,nlines*2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,phie_local,nlines*2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,L_local,nlines,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      END IF
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES, ierr_mpi)

      ! Restore data
      R_start = R_local
      Z_start = Z_local
      PHI_start = phis_local
      PHI_end   = phie_local

      ! now deallocate the helpers
      IF (ASSOCIATED(R_local))    CALL mpidealloc(R_local,    win_R_local)
      IF (ASSOCIATED(Z_local))    CALL mpidealloc(Z_local,    win_Z_local)
      IF (ASSOCIATED(phis_local)) CALL mpidealloc(phis_local, win_phis_local)
      IF (ASSOCIATED(phie_local)) CALL mpidealloc(phie_local, win_phie_local)

      ! Change some things to get it right
      lmu = .FALSE. ! Turns off diffusion
      nlines = nlines*2 ! Twice the number of lines
      ldex_default = 1
      IF (myid_sharmem == master) ihit_array = 0 ! reset it to zero

      ! Now run the reverse trace
      IF (lverb) WRITE(6,'(A)') '===========BackFlow of Wall Hits=========='
      CALL fieldlines_follow

      ! Reset a thing or two
      lmu = .TRUE.

      ! Fix L_lines
      mystart = LBOUND(R_lines,1)
      myend   = UBOUND(R_lines,1)
      i = nlines/2
      IF (myend > i) myend = i
      IF (mystart <= i) L_lines(mystart:myend) = L_local(mystart:myend)
      IF (ASSOCIATED(L_local)) CALL mpidealloc(L_local, win_L_local)


!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_init_backflow