!-----------------------------------------------------------------------
!     Module:        fieldlines_follow
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/21/2012
!     Description:   This subroutine follows the fieldlines through
!                    the grid.  The general ODE which must be solved
!                    can be written:
!                        dX        B_X
!                       ----   = -------
!                       dphi      B_PHI
!                    where X=(R,Z).  This allows fieldlines trajectories
!                    to be parametrized in terms of the toroidal angle.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_follow
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE vessel_mod
      USE fieldlines_runtime
      USE fieldlines_lines
      USE fieldlines_grid, ONLY: phimin, phimax, delta_phi,&
                                 BR_spl, BZ_spl, MU_spl, MODB_spl
      USE wall_mod, ONLY: ihit_array, nface, wall_free
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          status       MPI stats indicator
!          mystart      Starting fieldline for each processor
!          myend        Ending fieldline for each processor
!          mypace       Number of fieldlines each processor should follow
!          i            Dummy index
!          ier          Error flag
!          l            Dummy index
!          neqs_nag     Number of ODE's to solve (not limited to NAG routines)
!          l2           Index helper
!          itol         LSODE tolerance flag
!          itask        LSODE flag (overshoot and interpolate)  
!          istate       LSODE restart flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER,ALLOCATABLE :: mnum(:)
      INTEGER :: mystart,mypace, i
      INTEGER     :: ier, l, neqs_nag, l2, itol, itask, &
                     istate, iopt, lrw, liw, mf
      INTEGER, ALLOCATABLE :: iwork(:)
      INTEGER :: MPI_COMM_LOCAL
      REAL(rprec) :: phif_max
      DOUBLE PRECISION, ALLOCATABLE :: w(:), q(:)
      DOUBLE PRECISION :: phif_nag, phi_nag,&
                          tol_nag, rtol
      DOUBLE PRECISION :: atol(2)
      DOUBLE PRECISION :: rkh_work(2,2)
      CHARACTER*1 :: relab
!-----------------------------------------------------------------------
!     External Functions
!          fblin_nag            RHS of ODE integrator (for NAG)
!          fblin_lsode          RHS of ODE integrator (for LSODE)
!          fblin_rkh68          RHS of ODE integrator (for RKH68)
!          out_fieldlines_nag   Fieldline output (for NAG)
!          D02CJF               NAG ODE Solver
!          D02CJW               NAG Dummy function
!          jacobian_lsode       Jacobian function (for LSODE, not currently utilized)
!-----------------------------------------------------------------------
      EXTERNAL fblin_nag, out_fieldlines_nag, D02CJF, D02CJW
      EXTERNAL fblin_lsode, jacobian_lsode
      EXTERNAL fblin_rkh68
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Initializations (not all points start at the same spot so find the longest path)
      !    later we let the sign of dphi depend on phi_end and phi_start. (note reverse
      !    done in init routine.
      phif_max = MAXVAL(ABS(phi_end-phi_start))
      dphi = phimax/npoinc
      nsteps = CEILING(phif_max/dphi)
      delta_phi = phimax-phimin
      tol_nag = follow_tol
      neqs_nag = 2
      relab = "M"
      !mf = 10
      mf=21

      ! Break up the work
      mypace  = nlines
      mystart = 1
      myend   = nlines
      IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
      ALLOCATE(mnum(nprocs_fieldlines))
      mnum = FLOOR(REAL(nlines)/REAL(nprocs_fieldlines))
      DO
         IF (SUM(mnum) == nlines) EXIT
         DO i = 1, nprocs_fieldlines
            mnum(i) = mnum(i) + 1
            IF (SUM(mnum) == nlines) EXIT
         END DO
         IF (SUM(mnum) == nlines) EXIT
      END DO
      mypace = mnum(myworkid+1)
      IF (myworkid == 0) THEN
         mystart = 1
      ELSE
         mystart = SUM(mnum(1:myworkid))+1
      END IF
      myend   = mystart+mypace-1
      DEALLOCATE(mnum)

      ! Set the helper arrays
      IF (lhitonly) nsteps=2
      
      ! Deallocations
      IF (ALLOCATED(q)) DEALLOCATE(q)
      IF (ALLOCATED(w)) DEALLOCATE(w)
      IF (ALLOCATED(R_lines)) DEALLOCATE(R_lines)
      IF (ALLOCATED(Z_lines)) DEALLOCATE(Z_lines)
      IF (ALLOCATED(PHI_lines)) DEALLOCATE(PHI_lines)
      IF (ALLOCATED(B_lines)) DEALLOCATE(B_lines)
      IF (ALLOCATED(L_lines)) DEALLOCATE(L_lines)
      
      ! Allocations
      ALLOCATE(q(neqs_nag),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'Q',ier)
      ALLOCATE(R_lines(mystart:myend,0:nsteps),Z_lines(mystart:myend,0:nsteps),PHI_lines(mystart:myend,0:nsteps),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R_LINES, Z_LINES PHI_LINES',ier)
      ALLOCATE(B_lines(mystart:myend,0:nsteps),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'B_LINES',ier)
      ALLOCATE(L_lines(mystart:myend),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'L_LINES',ier)
      R_lines=0.0
      Z_lines=0.0
      L_lines=0.0
      PHI_lines=-1.0
      B_lines = -1
      R_lines(mystart:myend,0)   = R_start(mystart:myend)
      Z_lines(mystart:myend,0)   = Z_start(mystart:myend)
      PHI_lines(mystart:myend,0) = phi_start(mystart:myend)

      ! Output some stuff
      IF (lverb) THEN
         WRITE(6,'(A)')                '----- FOLLOWING FIELD LINES -----'
         WRITE(6,'(A,A)')              '      Method: ',TRIM(int_type)
         WRITE(6,'(A,I7)')             '       Lines: ',nlines
         IF (lhitonly) THEN
            WRITE(6,'(A)')                '       HITONLY: only recording strike points'
         ELSE
            WRITE(6,'(A,I7,A,E11.4)')     '       Steps: ',nsteps,'   Delta-phi: ',dphi
         END IF
         SELECT CASE(TRIM(int_type))
            CASE("NAG")
               WRITE(6,'(A,ES11.4,A,A1)')     '         Tol: ',follow_tol,'  Type: ',relab
            CASE("LSODE")
               WRITE(6,'(A,ES11.4,A,I2)')     '         Tol: ',follow_tol,'  Type: ',mf
         END SELECT
         IF (lreverse) WRITE(6,'(A)')                '       REVERSE: following reversed path'
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
                  dphi = SIGN(dphi,phi_end(l)-phi_start(l))  ! This should work 
                  ier     = 1
                  myline  = l
                  myldex  = 0
                  CALL D02CJF(phi_nag,phif_nag,neqs_nag,q,fblin_nag,tol_nag,relab,out_fieldlines_nag,D02CJW,w,ier)
                  IF (ier < 0) CALL handle_err(D02CJF_ERR,'fieldlines_follow',ier)
               END DO
!DEC$ ELSE
               ier = -1
               CALL handle_err(NAG_ERR,'fieldlines_follow',ier)
!DEC$ ENDIF  
            CASE ("RKH68")
               ier = 0
               DO l = mystart, myend
                  q(1)    = R_start(l)
                  q(2)    = Z_start(l)
                  phi_nag = phi_start(l)
                  ier     = 0
                  myline  = l
                  l2 = 0
                  myldex  = 0
                  dphi = SIGN(dphi,phi_end(l)-phi_start(l))  ! This should work 
                  CALL out_fieldlines_nag(phi_nag,q)
                  phi_nag = phi_start(l)             ! Because out_fieldlines_nag increments phi_nag
                  DO 
                     CALL drkhvg(phi_nag,q,neqs_nag,dphi,2,fblin_rkh68,rkh_work,iopt,ier)
                     IF (ier < 0) CALL handle_err(RKH68_ERR,'fieldlines_follow',ier)
                     q(1)=rkh_work(1,2)
                     q(2)=rkh_work(2,2)
                     phi_nag = phi_nag+dphi  ! Because phi_nag wasn't incremented by drkhvg
                     CALL out_fieldlines_nag(phi_nag,q)
                     IF (phi_nag > phi_end(l)) EXIT
                  END DO
               END DO
            CASE ("LSODE","DLSODE")
               ! Set tolerances to match NAG
               itol = 2; rtol = follow_tol; atol(:) = follow_tol
               iopt = 1 ! Optional output
               IF (mf == 10) THEN
                  lrw = 20 + 16 * neqs_nag
                  liw = 20
               ELSE IF (mf <23) THEN
                  lrw = 22 +  9 * neqs_nag +  neqs_nag**2
                  liw = 20 +  neqs_nag
               ELSE
                  STOP 'mf > 23, banded jacobian is not possible for this problem'
               END IF
               ALLOCATE(w(lrw))
               ALLOCATE(iwork(liw))
               ier = 0
               CALL xsetf(ier)
               DO l = mystart, myend
                  ! Setup LSODE parameters
                  w = 0; iwork = 0; itask = 1; istate = 1;
                  myline  = l
                  myldex  = 0
                  ier     = 0
                  l2 = 0
                  q(1)    = R_start(l)
                  q(2)    = Z_start(l)
                  phi_nag = phi_start(l)
                  phif_nag = phi_start(l)
                  dphi = SIGN(dphi,phi_end(l)-phi_start(l))  ! This should work 
                  DO 
                     IF (lmu) istate = 1
                     CALL DLSODE(fblin_lsode,neqs_nag,q,phi_nag,phif_nag,itol,rtol,atol,itask,istate,&
                                iopt,w,lrw,iwork,liw,jacobian_lsode,mf)
                     IF (istate < -2) CALL handle_err(LSODE_ERR,'fieldlines_follow',istate)
                     CALL out_fieldlines_nag(phif_nag,q)
                     iwork(11) = 0; iwork(12)=0; iwork(13)=0;
                     IF ((istate == -1) .or. (istate ==-2) .or. (ABS(phif_nag) > ABS(phi_end(l))) ) EXIT
                  END DO
               END DO
            CASE ("TEST")
               DO l = 0, nsteps
                  R_lines(mystart:myend,l) = REAL(l)
               END DO
               B_lines(mystart:myend,0:nsteps) = REAL(myworkid)
         END SELECT
      END IF
      
      IF (myworkid == master) WRITE(6,*) ' '

      ! Deallocations
      IF (ALLOCATED(q)) DEALLOCATE(q)
      IF (ALLOCATED(w)) DEALLOCATE(w)
      IF (ALLOCATED(iwork)) DEALLOCATE(iwork)

      ! Handle WALL Heat Map
#if defined(MPI_OPT)
    IF (ASSOCIATED(ihit_array)) THEN
      i = MPI_UNDEFINED
      IF (myid_sharmem == master) i = 0
      CALL MPI_COMM_SPLIT( MPI_COMM_FIELDLINES,i,myworkid,MPI_COMM_LOCAL,ierr_mpi)
      IF (myid_sharmem == master) THEN
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,ihit_array,nface,MPI_INTEGER,MPI_SUM,MPI_COMM_LOCAL,ierr_mpi)
         CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      END IF
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES, ierr_mpi)
    END IF
#endif
   
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_follow
