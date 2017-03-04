!-----------------------------------------------------------------------
!     Module:        beams3d_flux
!     Authors:       M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          07/02/2012
!     Description:   This subroutine converts the particle trajectory 
!                     data from R,PHI,Z to Boozer (s,u,v) coordinates.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_flux
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_lines, b_temp_mask => b_temp
      USE beams3d_grid, ONLY: nr, nphi, nz, B_R, B_PHI, B_Z
      USE beams3d_runtime, ONLY: lverb, lflux, nprocs_beams, npoinc, &
                                 MPI_BARRIER_ERR, MPI_RECV_ERR, &
                                 MPI_SEND_ERR, BEAMS3D_TRANSMIT_2DDBL,&
                                 ALLOC_ERR
      USE vmec_utils
      USE read_wout_mod
!DEC$ IF DEFINED (MPI_OPT)
      USE mpi_params ! MPI
!DEC$ ENDIF
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error Flag
!          iunit        File ID
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, iunit, j, k
      REAL(rprec) :: r_temp, phi_temp, z_temp, s_temp, u_temp, v_temp
      REAL(rprec) :: b_temp(3)
!DEC$ IF DEFINED (MPI_OPT)
!      INCLUDE 'mpif.h' ! MPI
      INTEGER :: mystart, i, sender
      !INTEGER :: status(MPI_STATUS_size) !mpi stuff
      !DOUBLE PRECISION, ALLOCATABLE :: buffer_slav(:,:), buffer_mast(:,:)
      INTEGER,ALLOCATABLE :: mnum(:), moffsets(:)
!DEC$ ENDIF
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      IF (.not.lflux) RETURN

!DEC$ IF DEFINED (MPI_OPT)
      mystart = mystart_save
      myend = myend_save

      IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
      IF (ALLOCATED(moffsets)) DEALLOCATE(moffsets)
      ALLOCATE(mnum(nprocs_beams), moffsets(nprocs_beams))
      CALL MPI_ALLGATHER((myend-mystart+1)*(npoinc+1),1,MPI_INTEGER,mnum,1,MPI_INTEGER,MPI_COMM_BEAMS,ierr_mpi)
      CALL MPI_ALLGATHER((mystart-1)*(npoinc+1),1,MPI_INTEGER,moffsets,1,MPI_INTEGER,MPI_COMM_BEAMS,ierr_mpi)
!DEC$ ENDIF

!DEC$ IF DEFINED (OLD_FLUX)
      IF (lverb) WRITE(6,'(A)')  '----- CONVERTING TO FLUX COORDINATES -----'
      IF (ALLOCATED(moment_lines)) DEALLOCATE(moment_lines)
      IF (myworkid == master) THEN
        ALLOCATE(S_lines(0:npoinc, nparticles), U_lines(0:npoinc, nparticles), &
                 V_lines(0:npoinc, nparticles) )
      ELSE
        ALLOCATE(S_lines(0:npoinc, mystart:myend), U_lines(0:npoinc, mystart:myend), &
                 V_lines(0:npoinc, mystart:myend) )
      END IF
      S_lines = 1.5
      U_lines = 0.0
      V_lines = 0.0
      IF (lverb) THEN
         WRITE(6, '(5X,A,I3,A)', ADVANCE = 'no') 'Coordinate Calculation [', 0, ']%'
         CALL FLUSH(6)
      END IF
      IF (mystart <= nparticles) THEN
         DO j=mystart,myend
            DO k=0,npoinc
               r_temp   = R_lines(k,j)
               phi_temp = PHI_lines(k,j)
               z_temp   = Z_lines(k,j)
               ier = 0
               IF ((phi_temp==-1) .and. (r_temp==0) .and. (z_temp==0)) EXIT
               s_temp = 0.1; u_temp=0.5
               CALL GetBcyl(r_temp,phi_temp,z_temp,&
                               b_temp(1),b_temp(2),b_temp(3),&
                               SFLX=s_temp,UFLX=u_temp,INFO=ier)
               IF (ier == 0) THEN
                  v_temp=phi_temp
                  S_lines(k,j) = s_temp
                  U_lines(k,j) = u_temp
!                  V_lines(k,j) = v_temp
               END IF      
            END DO
            IF (lverb) THEN
               CALL backspace_out(6,6)
               WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT(100.*j/myend),']%'
               CALL FLUSH(6)
            END IF
         END DO
      END IF

!     Clean up progress bar
      IF (lverb) THEN
         CALL backspace_out(6, 33)
         WRITE(6, '(33X)', ADVANCE = 'no')
         CALL backspace_out(6, 33)
         CALL FLUSH(6)
      END IF


!DEC$ ENDIF

!DEC$ IF DEFINED (MPI_OPT)
      IF (myworkid==master) THEN
         mystart = 1; myend=nparticles
      END IF
      CALL BEAMS3D_TRANSMIT_2DDBL(0,npoinc,mystart,myend,S_lines(0:npoinc,mystart:myend),&
                                  nprocs_beams,mnum,moffsets,myworkid,master,MPI_COMM_BEAMS,ier)
!      IF (myworkid==master) PRINT *,'S_LINES'
      CALL BEAMS3D_TRANSMIT_2DDBL(0,npoinc,mystart,myend,U_lines(0:npoinc,mystart:myend),&
                                  nprocs_beams,mnum,moffsets,myworkid,master,MPI_COMM_BEAMS,ier)
!      IF (myworkid==master) PRINT *,'U_LINES'
!      CALL BEAMS3D_TRANSMIT_2DDBL(0,npoinc,mystart,myend,V_lines(0:npoinc,mystart:myend),&
!                                  nprocs_beams,mnum,moffsets,myworkid,master,MPI_COMM_BEAMS,ier)
!      IF (myworkid==master) PRINT *,'V_LINES'
!    CALL MPI_BARRIER(MPI_COMM_BEAMS, ierr_mpi)
!    IF (ierr_mpi /= 0) CALL handle_err(MPI_BARRIER_ERR, 'beams3d_follow', ierr_mpi)
!    IF (myworkid == master) THEN
!        ALLOCATE(buffer_mast(0:npoinc, 3), STAT=ier)
!        IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'buffer_mast', ier)
!        DO i = myend + 1, nparticles
!            CALL MPI_RECV(buffer_mast, 3 * (npoinc + 1), MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
!                          MPI_ANY_TAG, MPI_COMM_BEAMS, status, ierr_mpi)
!            IF (ierr_mpi /= 0) CALL handle_err(MPI_RECV_ERR, 'beams3d_follow', ierr_mpi)
!            sender = status(MPI_SOURCE)
!            j = status(MPI_TAG)
!            S_lines(:, j) = buffer_mast(:, 1)
!            U_lines(:, j) = buffer_mast(:, 2)
!            V_lines(:, j) = buffer_mast(:, 3)
!        END DO
!        DEALLOCATE(buffer_mast)
!    ELSE
!        IF (mystart <= nparticles) THEN
!        ALLOCATE(buffer_slav(0:npoinc, 3), STAT=ier)
!        IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'buffer_slav', ier)
!        DO j = mystart, myend
!            buffer_slav(:, 1) = S_lines(:, j)
!            buffer_slav(:, 2) = U_lines(:, j)
!            buffer_slav(:, 3) = V_lines(:, j)
!            CALL MPI_SEND(buffer_slav, 3 * (npoinc + 1), MPI_DOUBLE_PRECISION, master, j, MPI_COMM_BEAMS, ierr_mpi)
!            IF (ierr_mpi /= 0) CALL handle_err(MPI_SEND_ERR, 'beams3d_follow', ierr_mpi)
!        END DO
!        DEALLOCATE(buffer_slav)
!        END IF
!    END IF
    DEALLOCATE(mnum)
    DEALLOCATE(moffsets)
    CALL MPI_BARRIER(MPI_COMM_BEAMS, ierr_mpi)
    IF (ierr_mpi /= 0) CALL handle_err(MPI_BARRIER_ERR, 'beams3d_follow', ierr_mpi)
!DEC$ ENDIF

      IF (myworkid .ne. master) CALL read_wout_deallocate

      CALL beams3d_write('RHO_TRAJECTORY')

      IF (ALLOCATED(U_lines)) DEALLOCATE(U_lines)
!      IF (ALLOCATED(V_lines)) DEALLOCATE(V_lines)
      IF (ALLOCATED(Z_lines)) DEALLOCATE(Z_lines)


!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_flux
