      SUBROUTINE de_mpi(np, fcn, funcval)
      USE de_mod, ONLY: rprec, n_free,nopt,ui_XC,nfev
      USE mpi_params
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      include 'mpif.h'                                       !mpi stuff
!DEC$ ENDIF
      INTEGER, INTENT(in) :: np
      REAL(rprec),INTENT(inout) :: funcval(np)
      EXTERNAL fcn

!DEC$ IF DEFINED (MPI_OPT)
      INTEGER :: status(MPI_STATUS_size)                     !mpi stuff
      INTEGER :: i,  j, iflag, ikey                                !mpi stuff
      INTEGER :: numsent, sender, ierr                       !mpi stuff
      INTEGER :: anstype, column                             !mpi stuff
      REAL(rprec) :: chisq_temp
      REAL(rprec), DIMENSION(n_free) :: x
      REAL(rprec), DIMENSION(nopt) :: fvec

      
      fvec = 0; x=0; ierr = 0; ierr_mpi = 0; iflag = 0
!
!******************************************
!
!  mpi : set barrier so ALL processors get here before starting
!
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)                 !mpi stuff
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
      call MPI_COMM_RANK( MPI_COMM_STEL, myid, ierr_mpi )       !mpi stuff
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
      call MPI_COMM_SIZE( MPI_COMM_STEL, numprocs, ierr_mpi )   !mpi stuff
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)

!******************************************
!
!     ****Master node controls this portion of the code****
!     
      IF (myid .eq. master) THEN
         numsent = 0
         DO j = 1,MIN(numprocs-1,np)
             CALL MPI_SEND(j,1,MPI_INTEGER,j,j,MPI_COMM_STEL,ierr_mpi)
             IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
             numsent = numsent + 1
         END DO

         ! Now Enter the que loop
         DO j = 1, np
            CALL MPI_RECV(chisq_temp,1, MPI_REAL8, MPI_ANY_SOURCE,
     1                    MPI_ANY_TAG,MPI_COMM_STEL,status,ierr_mpi)
            IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
            sender     = status(MPI_SOURCE)
            anstype    = status(MPI_TAG)
            IF (anstype .gt. np) STOP 'ANSTYPE > NP IN de_mpi'

            funcval(anstype) = chisq_temp
            write (6, '(2x,i6,8x,i3,7x,1es12.4)') anstype,
     1             sender, funcval(anstype)
            CALL FLUSH(6)

            ! Send more work
            IF (numsent .lt. np) THEN
               numsent = numsent+1
               CALL MPI_SEND(numsent, 1, MPI_INTEGER,
     1                      sender, numsent, MPI_COMM_STEL, ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS)
     1              CALL mpi_stel_abort(ierr_mpi)
            ELSE
               column = 0
               CALL MPI_SEND(column, 1, MPI_INTEGER,
     1                      sender, numsent, MPI_COMM_STEL, ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS) 
     1             CALL mpi_stel_abort(ierr_mpi)
            ENDIF
         END DO
c     ****Worker portion of the code****
c        Skip this when processor id exceeds work to be done
      ELSE IF (myid .le. np) THEN
         DO
            j=-1
            CALL MPI_RECV(j,1,MPI_INTEGER,master,MPI_ANY_TAG,
     1                 MPI_COMM_STEL,status,ierr_mpi)
            PRINT *,j
            IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
            IF (j .gt. 0) THEN
               x(:) = ui_XC(j,:)
               fvec(:) = 0
               iflag   = j
               !PRINT *,'myid',nopt, n_free, x, fvec, iflag, nfev
               !CALL fcn(nopt, n_free, x, fvec, iflag, nfev)
               chisq_temp = SUM(fvec*fvec)
               chisq_temp = 5*i+chisq_temp*.25
               CALL MPI_SEND(chisq_temp,1, MPI_REAL8, master,
     1                       j, MPI_COMM_STEL, ierr_mpi)
               IF (ierr_mpi /= MPI_SUCCESS)
     1                  CALL mpi_stel_abort(ierr_mpi)
            ELSE 
               EXIT
            END IF
         END DO
      END IF

!
!     Broadcast the funcval array to all processors FROM master
!
      PRINT *,myid,' AT BARRIER'; CALL FLUSH(6)
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr_mpi)                 !mpi stuff
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(funcval, np, MPI_REAL8, master,
     1     MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL mpi_stel_abort(ierr_mpi)
      IF (ierr_mpi .ne. 0) GOTO 100

      RETURN

 100  CONTINUE
      PRINT *,' MPI_BCAST error in de_mpi: IERR=', ierr

      RETURN

 300  CONTINUE
      PRINT *,' IFLAG = ', iflag, ' in de_mpi CALL to fcn'
      STOP
!DEC$ ENDIF
      END SUBROUTINE de_mpi
