!-----------------------------------------------------------------------
!     Subroutine:    XVEC_EVAL
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/05/2014
!     Description:   This subroutine reads the xvector arrays from an
!                    xvec.dat file and evalutes each equilibria, saving
!                    it to individual files.
!-----------------------------------------------------------------------
      SUBROUTINE XVEC_EVAL(fcn,nvar,mvar,filename)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds
      USE mpi_params
      USE gade_mod, ONLY: pso_cleanup
      USE safe_open_mod
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      include 'mpif.h'                                       !mpi stuff
!DEC$ ENDIF
!-----------------------------------------------------------------------
!     Input Variables
!        fcn      Target function one wishes to evaluate.
!        nvar     Size of xvec array
!        filename Filename of file to read for xvec.
!-----------------------------------------------------------------------
      INTEGER :: nvar, mvar
      CHARACTER(*) :: filename
      EXTERNAL  fcn
!-----------------------------------------------------------------------
!     Local Variables
!        i           Dummy index
!        iunit       Unit number for file access
!        istat       Status dummy var
!        iter        Iteration dummy
!        chunk       Chunk of work to do
!        mystart     Starting Position
!        myend       Ending Position
!        iter_array  Iteration number array
!        n           XVEC size dummy
!        x_temp      X-vector temporary array
!        fvec_temp   F_vector temporary array
!        xvec        X-vector Array
!        chisq       Chi-Squared
!        nexec       Number of equilibria to evaluate
!-----------------------------------------------------------------------
      INTEGER :: iunit, istat, iter, n, nexec, chunk, mystart, myend, i
      INTEGER, ALLOCATABLE :: iter_array(:)
      REAL(rprec) :: chisq
      REAL(rprec), ALLOCATABLE :: x_temp(:), fvec_temp(:), xvec(:,:)
!-----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!-----------------------------------------------------------------------

      ! Initialize from file
      ALLOCATE(x_temp(nvar),fvec_temp(mvar),STAT=istat)

      CALL safe_open(iunit,istat,filename,'unknown','formatted')
      nexec=0
      DO
         READ(iunit,FMT=*,IOSTAT=istat) n,iter
         IF (istat .ne. 0) EXIT
         READ(iunit,FMT=*,IOSTAT=istat) x_temp(1:n)
         IF (istat .ne. 0) EXIT
         READ(iunit,FMT=*,IOSTAT=istat) chisq
         IF (istat .ne. 0) EXIT
         nexec=nexec+1
      END DO

      ALLOCATE(iter_array(nexec),STAT=istat)
      ALLOCATE(xvec(nvar,nexec),STAT=istat)
      REWIND(iunit)
      DO i = 1, nexec
         READ(iunit,*) n,iter_array(i)
         READ(iunit,*) xvec(1:n,i)
         READ(iunit,*) chisq
      END DO

      CLOSE(iunit)

      ! Evaluate the equilibria
      chunk = FLOOR(REAL(nexec)/REAL(numprocs))
      mystart = myid*chunk+1
      myend   = mystart + chunk -1
      DO i = mystart, myend
         ! ITERATE
         x_temp(:) = xvec(:,i)
         istat = myid + 1
         iter = iter_array(i)
         CALL fcn(mvar, nvar, x_temp, fvec_temp, istat, iter)
         ! Save
         istat = pso_cleanup
         CALL fcn(mvar, nvar, x_temp, fvec_temp, istat, iter)
      END DO

      ! DEALLOCATE
      DEALLOCATE(x_temp, fvec_temp)
      DEALLOCATE(iter_array, xvec)

      RETURN
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      END SUBROUTINE XVEC_EVAL