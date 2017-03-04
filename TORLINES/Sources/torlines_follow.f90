!-----------------------------------------------------------------------
!     Subroutine:    torlines_follow
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/2/2011
!     Description:   This subroutine calls the ODE routines to follow
!                    fieldlines.  First we assume we want to take steps 
!                    in the toroidal direction which gauruntee we can
!                    reconstruct our highest n mode.  Then we need to
!                    make enough transits through a field period that
!                    we resolve the poloidal structure.  
!-----------------------------------------------------------------------
      SUBROUTINE torlines_follow
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE torlines_fieldlines
      USE torlines_background
      USE torlines_realspace
      USE torlines_runtime
      USE ez_hdf5
      USE mpi_params 
!-----------------------------------------------------------------------
!     Local Variables
!          ierr        Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h' 
      INTEGER :: status(MPI_STATUS_size), color, key
      INTEGER,ALLOCATABLE :: mnum(:), moffsets(:)
      REAL(rprec), ALLOCATABLE :: buffer_mast(:,:), buffer_slav(:,:)
!DEC$ ENDIF  
      INTEGER :: mypace, i,j, sender
      INTEGER     :: ier, l, neqs_nag, itol, itask, &
                     istate, iopt, lrw, liw, mf, out
      INTEGER, ALLOCATABLE :: iwork(:)
      REAL        :: dist
      REAL(rprec) :: phif_max, phif_min
      REAL(rprec), ALLOCATABLE :: w(:), q(:), RZ_send(:,:)
      DOUBLE PRECISION :: phif_nag, eps_temp, phi_nag, phi1_nag,&
                          tol_nag, rtol
      DOUBLE PRECISION :: atol(2), rwork(52)
      DOUBLE PRECISION :: rkh_work(2,2)
      CHARACTER*1 :: relab
!-----------------------------------------------------------------------
!     External Functions
!          fphif       Error estimate given a field line length
!          C05AJF      Zero of a function using a continuation method
!-----------------------------------------------------------------------
      EXTERNAL fphif,fblin,out_torlines,fblin_lsode,jacobian_lsode
      EXTERNAL C05AJF,D02CJF,D02CJX,D02CJW
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Initializations
      phif_max = MAXVAL(ABS(phi_end),MASK=(R_start >= 0))
      IF (lemc3) THEN
         phif_min = MINVAL(ABS(phi_start),MASK=(R_start >= 0))
         dphi = (phif_max-phif_min)/npoinc
         nsteps = npoinc
      ELSE
         dphi = phmx/npoinc
         nsteps = INT(phif_max/dphi)
      END IF
      tol_nag = follow_tol
      neqs_nag = 2
      relab = "M"
      mf = 10
      ! Break up the Work
      mypace = FLOOR(REAL(nlines) / REAL(numprocs))
      mystart = myid*mypace + 1
      myend = mystart + mypace - 1
      ! This section sets up the work so we can use GATHERV
!DEC$ IF DEFINED (MPI_OPT)
      IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
      IF (ALLOCATED(moffsets)) DEALLOCATE(moffsets)
      ALLOCATE(mnum(numprocs), moffsets(numprocs))
      CALL MPI_ALLGATHER(mypace,1,MPI_INTEGER,mnum,1,MPI_INTEGER,MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_ALLGATHER(mystart,1,MPI_INTEGER,moffsets,1,MPI_INTEGER,MPI_COMM_FIELDLINES,ierr_mpi)
      WHERE (moffsets > nlines) moffsets = nlines
      l = 1
      DO
         IF ((moffsets(numprocs)+mnum(numprocs)-1) == nlines) EXIT
         IF (l == numprocs) l = 1
         mnum(l) = mnum(l) + 1
         moffsets(l+1:numprocs) = moffsets(l+1:numprocs) + 1
         l=l+1
      END DO
      mystart = moffsets(myid+1)
      mypace  = mnum(myid+1)
      myend   = mystart + mypace - 1
!DEC$ ENDIF
      
      ! Deallocations
      IF (ALLOCATED(q)) DEALLOCATE(q)
      IF (ALLOCATED(w)) DEALLOCATE(w)
      IF (ALLOCATED(R_lines)) DEALLOCATE(R_lines)
      IF (ALLOCATED(Z_lines)) DEALLOCATE(Z_lines)
      IF (ALLOCATED(PHI_lines)) DEALLOCATE(PHI_lines)
      IF (ALLOCATED(U_lines)) DEALLOCATE(U_lines)
      IF (ALLOCATED(goodline)) DEALLOCATE(goodline)
      
      ! Allocations
      ALLOCATE(q(neqs_nag),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'Q',ier)
      IF (myid == master) THEN
         ALLOCATE(R_lines(nlines,0:nsteps),Z_lines(nlines,0:nsteps),PHI_lines(nlines,0:nsteps),&
                  B_lines(nlines,0:nsteps),U_lines(nlines,0:nsteps),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R_LINES, Z_LINES PHI_LINES',ier)
         ALLOCATE(goodline(nlines),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'hitwall',ier)
         R_lines=0.0
         Z_lines=0.0
         PHI_lines=-1.0
         goodline(:) = .true.
         R_lines(mystart:myend,0)   = R_start(mystart:myend)
         Z_lines(mystart:myend,0)   = Z_start(mystart:myend)
         PHI_lines(mystart:myend,0) = phi_start(mystart:myend)
         U_lines(mystart:myend,0)   = 0
         B_lines(mystart:myend,0)   = 0
      ELSE
         IF (mystart <= nlines) THEN
            ALLOCATE(R_lines(mystart:myend,0:nsteps),Z_lines(mystart:myend,0:nsteps),PHI_lines(mystart:myend,0:nsteps),&
                     B_lines(mystart:myend,0:nsteps),U_lines(mystart:myend,0:nsteps),STAT=ier)
            IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R_LINES, Z_LINES PHI_LINES',ier)
            ALLOCATE(goodline(mystart:myend),STAT=ier)
            IF (ier /= 0) CALL handle_err(ALLOC_ERR,'hitwall',ier)
            R_lines=0.0
            Z_lines=0.0
            PHI_lines=-1.0
            goodline(:) = .true.
            R_lines(mystart:myend,0)   = R_start(mystart:myend)
            Z_lines(mystart:myend,0)   = Z_start(mystart:myend)
            PHI_lines(mystart:myend,0) = phi_start(mystart:myend)
            U_lines(mystart:myend,0)   = 0
            B_lines(mystart:myend,0)   = 0
         END IF
      END IF
         
      ! Output some stuff
      IF (lverb) THEN
         WRITE(6,'(A)')                '  ----- FOLLOWING FIELD LINES -----'
         WRITE(6,'(A,A)')              '      Method: ',TRIM(int_type)
         WRITE(6,'(A,I6)')             '       Lines: ',nlines
         WRITE(6,'(A,I7,A,E11.4)')     '       Steps: ',nsteps,'   Delta-phi: ',dphi
         SELECT CASE(TRIM(int_type))
            CASE("NAG")
               WRITE(6,'(A,E11.4,A,A1)')     '         Tol: ',follow_tol,'  Type: ',relab
            CASE("LSODE")
               WRITE(6,'(A,E11.4,A,I2)')     '         Tol: ',follow_tol,'  Type: ',mf
         END SELECT
         WRITE(6,'(5X,A,I3,A)',ADVANCE='no') 'Fieldline Calculation [',0,']%'
         CALL FLUSH(6)
      END IF
      
      
      ! Follow Fieldlines
      IF (mystart <= nlines) THEN
         SELECT CASE (TRIM(int_type))
            CASE ("NAG")
!DEC$ IF DEFINED (NAG)
	       ALLOCATE(w(neqs_nag*21+28),STAT=ier)
               IF (ier /= 0) CALL handle_err(ALLOC_ERR,'W',ier)
               DO l = mystart, myend
                  q(1)    = R_start(l)
                  q(2)    = Z_start(l)
                  phi_nag = phi_start(l)
                  phif_nag = phi_end(l)
                  dphi    = ABS(dphi)
                  if (phi_nag > phif_nag) dphi = -dphi
                  ier     = 1
                  myline  = l
                  myldex  = 0
                  CALL D02CJF(phi_nag,phif_nag,neqs_nag,q,fblin,tol_nag,relab,out_torlines,D02CJW,w,ier)
                  IF (ier < 0) CALL handle_err(D02CJF_ERR,'torlines_follow',ier)
               END DO
!DEC$ ELSE
               ier = -1
               CALL handle_err(NAG_ERR,'torlines_follow',ier)
!DEC$ ENDIF 
            CASE ("LSODE","DLSODE")
               ! Set tolerances to match NAG
               itol = 2; rtol = follow_tol; atol(:) = follow_tol
               iopt = 1 ! Optional settings used
               lrw = 20 + 16 * neqs_nag
               liw = 20
               ALLOCATE(w(lrw))
               ALLOCATE(iwork(liw))
               ier = 0
               DO l = mystart, myend
                  ! Setup LSODE parameters
                  w = 0; iwork = 0; itask = 1; istate = 1;
                  iwork(6) = 2000 !Maximum number of steps
                  iwork(7) = 1    !Maximum messages per problem
                  myline  = l
                  myldex  = 0
                  ier     = 0
                  q(1)    = R_start(l)
                  q(2)    = Z_start(l)
                  phi_nag = phi_start(l)
                  phif_nag = phi_start(l)
                  CALL out_torlines(phif_nag,q)
                  dphi    = ABS(dphi)
                  if (phi_nag > phif_nag) dphi = -dphi
                  DO 
                     CALL DLSODE(fblin_lsode,neqs_nag,q,phi_nag,phif_nag,itol,rtol,atol,itask,istate,&
                                iopt,w,lrw,iwork,liw,jacobian_lsode,mf)
                     IF (istate < -1) CALL handle_err(LSODE_ERR,'torlines_follow',istate)
                     IF (istate == -1) istate = 2
                     CALL out_torlines(phif_nag,q)
                     IF (myldex > nsteps) EXIT
                  END DO
               END DO
         END SELECT
      END IF

!DEC$ IF DEFINED (MPI_OPT)
      ! Now add algather statement
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init',ierr_mpi)
      CALL MPI_GATHERV(goodline,mypace,MPI_LOGICAL,&
                       goodline,mnum,moffsets-1,MPI_LOGICAL,master,&
                       MPI_COMM_FIELDLINES,ierr_mpi)
!       ! Adjust indexing to send 2D arrays
!       mnum = mnum*(nsteps+1)
!       mypace = mypace*(nsteps+1)
!       moffsets = (moffsets-1)*(nsteps+1)+1
!       CALL MPI_GATHERV(R_lines,mypace,MPI_DOUBLE_PRECISION,&
!                        R_lines,mnum,moffsets,MPI_DOUBLE_PRECISION,&
!                        MPI_COMM_FIELDLINES,ierr_mpi)
!       CALL MPI_GATHERV(Z_lines,mypace,MPI_DOUBLE_PRECISION,&
!                        Z_lines,mnum,moffsets,MPI_DOUBLE_PRECISION,&
!                        MPI_COMM_FIELDLINES,ierr_mpi)
!       CALL MPI_GATHERV(PHI_lines,mypace,MPI_DOUBLE_PRECISION,&
!                        PHI_lines,mnum,moffsets,MPI_DOUBLE_PRECISION,&
!                        MPI_COMM_FIELDLINES,ierr_mpi)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init',ierr_mpi)
      IF (myid == master) THEN
         ALLOCATE(buffer_mast(5,0:nsteps))
         DO i = myend+1, nlines
            CALL MPI_RECV(buffer_mast,5*(nsteps+1),MPI_DOUBLE_PRECISION,&
                          MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_FIELDLINES,status,ierr_mpi)
            IF (ierr_mpi /=0) CALL handle_err(MPI_RECV_ERR,'fieldlines_init_mgrid',ierr_mpi)
            sender = status(MPI_SOURCE)
            j      = status(MPI_TAG)
            R_lines(j,:)   = buffer_mast(1,:) 
            Z_lines(j,:)   = buffer_mast(2,:)
            PHI_lines(j,:) = buffer_mast(3,:)
            B_lines(j,:) = buffer_mast(4,:)
            U_lines(j,:) = buffer_mast(5,:)
         END DO
         DEALLOCATE(buffer_mast)
      ELSE
         IF (mystart <= nlines) THEN
            ALLOCATE(buffer_slav(5,0:nsteps))
            DO j = mystart, myend
               buffer_slav(1,:) = R_lines(j,:)
               buffer_slav(2,:) = Z_lines(j,:)
               buffer_slav(3,:) = PHI_lines(j,:)
               buffer_slav(4,:) = B_lines(j,:)
               buffer_slav(5,:) = U_lines(j,:)
               CALL MPI_SEND(buffer_slav,5*(nsteps+1),MPI_DOUBLE_PRECISION,master,j,MPI_COMM_FIELDLINES,ierr_mpi)
               IF (ierr_mpi /=0) CALL handle_err(MPI_SEND_ERR,'fieldlines_init_mgrid',ierr_mpi)
            END DO
            DEALLOCATE(buffer_slav)
         END IF
         IF (ALLOCATED(R_lines)) DEALLOCATE(R_lines)
         IF (ALLOCATED(Z_lines)) DEALLOCATE(Z_lines)
         IF (ALLOCATED(PHI_lines)) DEALLOCATE(PHI_lines)
         IF (ALLOCATED(B_lines)) DEALLOCATE(B_lines)
         IF (ALLOCATED(U_lines)) DEALLOCATE(U_lines)
      END IF
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'fieldlines_init',ierr_mpi)
!DEC$ ENDIF

      IF (ALLOCATED(q)) DEALLOCATE(q)
      IF (ALLOCATED(w)) DEALLOCATE(w)
      IF (ALLOCATED(iwork)) DEALLOCATE(iwork)
      
      
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE torlines_follow
