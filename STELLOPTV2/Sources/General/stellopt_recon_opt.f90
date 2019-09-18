!-----------------------------------------------------------------------
!     Subroutine:    stellopt_recon_opt
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   This subroutine preforms a reconstruction
!                    type optimization
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_recon_opt(fcn,m,n,x,fvec,ier)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE equil_vals, ONLY: phiedge
      USE safe_open_mod, ONLY: safe_open
      USE mpi_params
      USE mpi_inc
      IMPLICIT NONE

!-----------------------------------------------------------------------
!     Interface of external callback function (from "stellopt_fcn.f90")
!-----------------------------------------------------------------------
      INTERFACE
         SUBROUTINE stellopt_fcn(m, n, x, fvec, iflag, ncnt) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND=C_INT), INTENT(in)      ::  m, n, ncnt
            INTEGER(KIND=C_INT), INTENT(inout)   :: iflag
            REAL(KIND=C_DOUBLE), INTENT(inout)  :: x(n)
            REAL(KIND=C_DOUBLE), INTENT(out) :: fvec(m)
         END SUBROUTINE stellopt_fcn
      END INTERFACE
 
!-----------------------------------------------------------------------
!     Input Variables    NONE
!-----------------------------------------------------------------------
      EXTERNAL fcn
      INTEGER, INTENT(in) :: m
      INTEGER, INTENT(in) :: n
      REAL(rprec), INTENT(inout) :: x(n)
      REAL(rprec), INTENT(inout) :: fvec(m)
      INTEGER, INTENT(inout) :: ier
      
!-----------------------------------------------------------------------
!     PARAMETERS     NONE    
!-----------------------------------------------------------------------
      INTEGER, PARAMETER :: NREFIT_MAX = 2
   
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      INTEGER :: iflag, ncnt, iproc_min, i, nrefit_stop
      REAL(rprec) :: chisq, chisq_old, delta_chisq, phiedge_local, &
                     pscale_local, chisq_new
      REAL(rprec) :: fvec_local(m)
      REAL(rprec) :: chisq_array(numprocs)
   
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF
   
      lrefit = .FALSE.
      chisq_old = bigno
      chisq_new = bigno
      IF (myid == master) THEN
         ! On first iteration call stellopt_fcn
         iflag = -1
         ncnt  = 0
         fvec_local = 0.0
         CALL stellopt_fcn(m,n,x,fvec_local,iflag,ncnt)
         CALL FLUSH(6)
         chisq_new = SQRT(SUM(fvec_local*fvec_local))
         delta_chisq = 1.0
         iflag = -102
         fvec_local = 0.0
         CALL fcn(m,n,x,fvec_local,iflag,ncnt)
         WRITE(6,'(A)') '--- Initial Profile Fitting----' 
         WRITE(6,'(2x,i6,8x,i3,7x,1es12.4)') ncnt, myid, chisq_new**2
         
         ! Now loop over equilibria till we've settled on a refit value
         iflag = 0
         lrefit = .TRUE.
         nrefit_stop = ncnt + NREFIT_MAX
         DO WHILE ((delta_chisq > 0.01) .and. (ncnt < nrefit_stop))
           ncnt = ncnt + 1
           chisq_old = chisq_new
           iflag = 0
           fvec_local = 0.0
           CALL stellopt_fcn(m,n,x,fvec_local,iflag,ncnt)
           chisq_new = SQRT(SUM(fvec_local*fvec_local))
           delta_chisq = ABS(chisq_new-chisq_old)/chisq_new
           WRITE(6,'(8x,8x,3x,7x,1es12.4)') chisq_new**2
         END DO
         WRITE(6,'(2x,i6,8x,i3,7x,1es12.4)') ncnt, myid, chisq_new**2
         
         ! Now output the file
         iflag = -102
         fvec_local = 0.0
         CALL fcn(m,n,x,fvec_local,iflag,ncnt)
         phiedge_local = phiedge
      END IF
      
      ! Make sure everyone has the proper values
!DEC$ IF DEFINED (MPI_OPT)
      ierr_mpi = 0
      CALL MPI_BARRIER(MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(chisq_new,1,MPI_REAL8,master,MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(ncnt,1,MPI_INTEGER,master,MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(ne_aux_s,ndatafmax,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(ne_aux_f,ndatafmax,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(te_aux_s,ndatafmax,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(te_aux_f,ndatafmax,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(ti_aux_s,ndatafmax,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(ti_aux_f,ndatafmax,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF
   
      ! Now do a spread in PHIEDGE
      lrefit = .TRUE.
      IF (myid == master) WRITE(6,'(A)') '--- Spread in PHIEDGE ----' 
      DO i = 1, n
         IF (arr_dex(i,2) == norm_dex) cycle
         IF (var_dex(i) == iphiedge) THEN
            phiedge_local = x(i) * (1.0 + 0.5*REAL(myid)/(numprocs-1))
            x(i) = phiedge_local
         END IF
      END DO
      ncnt = ncnt + myid + 1
      chisq_old = chisq_new
      iflag = myid
      fvec_local = 0.0
      CALL stellopt_fcn(m,n,x,fvec_local,iflag,ncnt)
      chisq_new = SQRT(SUM(fvec_local*fvec_local))
      chisq_array(:) = bigno
      chisq_array(1) = chisq_new
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_ALLGATHER(chisq_new, 1, MPI_REAL8, chisq_array,1, MPI_REAL8, MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF
      iproc_min = MINLOC(chisq_array,DIM=1)
      IF (myid == master) THEN
         DO i = 1, numprocs
           IF (i == iproc_min) THEN
              WRITE(6,'(8x,8x,i3,7x,1es12.4,A)') i, chisq_array(i)**2,'*'
           ELSE
              WRITE(6,'(8x,8x,i3,7x,1es12.4)') i, chisq_array(i)**2
           END IF
         END DO
      END IF
      IF (myid == iproc_min-1) THEN
         iflag = -102
         CALL fcn(m,n,x,fvec_local,iflag,ncnt)
      END IF
!DEC$ IF DEFINED (MPI_OPT)
      IF ((chisq_array(iproc_min) < chisq_old)) THEN
         CALL MPI_BCAST(phiedge_local,1,MPI_REAL8,iproc_min-1,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
         CALL MPI_BCAST(chisq_new,1,MPI_REAL8,iproc_min-1,MPI_COMM_STEL, ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
         CALL MPI_BCAST(ne_aux_s,ndatafmax,MPI_REAL8,iproc_min-1,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
         CALL MPI_BCAST(ne_aux_f,ndatafmax,MPI_REAL8,iproc_min-1,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
         CALL MPI_BCAST(te_aux_s,ndatafmax,MPI_REAL8,iproc_min-1,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
         CALL MPI_BCAST(te_aux_f,ndatafmax,MPI_REAL8,iproc_min-1,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
         CALL MPI_BCAST(ti_aux_s,ndatafmax,MPI_REAL8,iproc_min-1,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
         CALL MPI_BCAST(ti_aux_f,ndatafmax,MPI_REAL8,iproc_min-1,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      END IF
!DEC$ ENDIF
      DO i = 1, n
         IF (arr_dex(i,2) == norm_dex) cycle
         IF (var_dex(i) == iphiedge) THEN
            x(i) = phiedge_local
         END IF
      END DO
      ncnt = ncnt - myid + numprocs
      
      ! Now Refit
      lrefit = .TRUE.
      IF (myid == master) THEN
         WRITE(6,'(A)') '--- Profile Fitting ----' 
         ! Now loop over equilibria till we've settled on a refit value
         iflag = 0
         nrefit_stop = ncnt + NREFIT_MAX
         delta_chisq = 1.0
         DO WHILE ((delta_chisq > 0.01) .and. (ncnt < nrefit_stop))
           ncnt = ncnt + 1
           chisq_old = chisq_new
           iflag = 0
           CALL stellopt_fcn(m,n,x,fvec_local,iflag,ncnt)
           chisq_new = SQRT(SUM(fvec_local*fvec_local))
           delta_chisq = ABS(chisq_new-chisq_old)/chisq_new
           WRITE(6,'(8x,8x,3x,7x,1es12.4)') chisq_new**2
         END DO
         WRITE(6,'(2x,i6,8x,i3,7x,1es12.4)') ncnt, myid, chisq_new**2
         
         ! Now output the file
         iflag = -102
         CALL fcn(m,n,x,fvec_local,iflag,ncnt)
      END IF
      
      ! Make sure everyone has the proper values
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(chisq_new,1,MPI_REAL8,master,MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(ncnt,1,MPI_INTEGER,master,MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(ne_aux_s,ndatafmax,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(ne_aux_f,ndatafmax,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(te_aux_s,ndatafmax,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(te_aux_f,ndatafmax,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(ti_aux_s,ndatafmax,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(ti_aux_f,ndatafmax,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF
   
      ! Now do a spread in PRES_SCALE
      lrefit = .FALSE.
      IF (myid == master) WRITE(6,'(A)') '--- Spread in PRES_SCALE ----' 
      DO i = 1, n
         IF (arr_dex(i,2) == norm_dex) cycle
         IF (var_dex(i) == ipscale) THEN
            pscale_local = x(i)
            x(i) = pscale_local + pscale_local*(-0.5 + REAL(myid)/(numprocs-1))
         END IF
      END DO
      ncnt = ncnt + myid + 1
      chisq_old = chisq_new
      iflag = myid
      fvec_local = 0.0
      CALL stellopt_fcn(m,n,x,fvec_local,iflag,ncnt)
      chisq_new = SQRT(SUM(fvec_local*fvec_local))
      chisq_array(:) = bigno
      chisq_array(1) = chisq_new
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_ALLGATHER(chisq_new, 1, MPI_REAL8, chisq_array,1, MPI_REAL8, MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF
      iproc_min = MINLOC(chisq_array,DIM=1)
      IF (myid == master) THEN
         DO i = 1, numprocs
           IF (i == iproc_min) THEN
              WRITE(6,'(8x,8x,i3,7x,1es12.4,A)') i, chisq_array(i)**2,'*'
           ELSE
              WRITE(6,'(8x,8x,i3,7x,1es12.4)') i, chisq_array(i)**2
           END IF
         END DO
      END IF
      IF (myid == iproc_min-1) THEN
         iflag = -102
         CALL fcn(m,n,x,fvec_local,iflag,ncnt)
      END IF
!DEC$ IF DEFINED (MPI_OPT)
      IF ((chisq_array(iproc_min) < chisq_old)) THEN
         CALL MPI_BCAST(pscale_local,1,MPI_REAL8,iproc_min-1,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
         CALL MPI_BCAST(chisq_new,1,MPI_REAL8,iproc_min-1,MPI_COMM_STEL, ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
         CALL MPI_BCAST(ne_aux_s,ndatafmax,MPI_REAL8,iproc_min-1,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
         CALL MPI_BCAST(ne_aux_f,ndatafmax,MPI_REAL8,iproc_min-1,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
         CALL MPI_BCAST(te_aux_s,ndatafmax,MPI_REAL8,iproc_min-1,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
         CALL MPI_BCAST(te_aux_f,ndatafmax,MPI_REAL8,iproc_min-1,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
         CALL MPI_BCAST(ti_aux_s,ndatafmax,MPI_REAL8,iproc_min-1,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
         CALL MPI_BCAST(ti_aux_f,ndatafmax,MPI_REAL8,iproc_min-1,MPI_COMM_STEL,ierr_mpi)
         IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      END IF
!DEC$ ENDIF
      DO i = 1, n
         IF (arr_dex(i,2) == norm_dex) cycle
         IF (var_dex(i) == ipscale) THEN
            x(i) = pscale_local
         END IF
      END DO
      ncnt = ncnt - myid + numprocs
      
      ! Now Refit
      lrefit = .TRUE.
      IF (myid == master) THEN
         WRITE(6,'(A)') '--- Profile Fitting ----' 
         ! Now loop over equilibria till we've settled on a refit value
         iflag = 0
         nrefit_stop = ncnt + NREFIT_MAX
         delta_chisq = 1.0
         DO WHILE ((delta_chisq > 0.01) .and. (ncnt < nrefit_stop))
           ncnt = ncnt + 1
           chisq_old = chisq_new
           iflag = 0
           fvec_local = 0.0
           CALL stellopt_fcn(m,n,x,fvec_local,iflag,ncnt)
           chisq_new = SQRT(SUM(fvec_local*fvec_local))
           delta_chisq = ABS(chisq_new-chisq_old)/chisq_new
           WRITE(6,'(8x,8x,3x,7x,1es12.4)') chisq_new**2
         END DO
         WRITE(6,'(2x,i6,8x,i3,7x,1es12.4)') ncnt, myid, chisq_new**2
         
         ! Now output the file
         iflag = -102
         CALL fcn(m,n,x,fvec_local,iflag,ncnt)
      END IF
      
      ! Make sure everyone has the proper values
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(chisq_new,1,MPI_REAL8,master,MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(ncnt,1,MPI_INTEGER,master,MPI_COMM_STEL, ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(ne_aux_s,ndatafmax,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(ne_aux_f,ndatafmax,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(te_aux_s,ndatafmax,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(te_aux_f,ndatafmax,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(ti_aux_s,ndatafmax,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(ti_aux_f,ndatafmax,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
      CALL MPI_BCAST(x,n,MPI_REAL8,master,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi .ne. 0) CALL mpi_stel_abort(ierr_mpi)
!DEC$ ENDIF
      
      ier = ncnt
      RETURN
      
      
      
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_recon_opt
