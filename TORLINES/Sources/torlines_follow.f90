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
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ierr        Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
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
      CALL MPI_CALC_MYRANGE(MPI_COMM_FIELDLINES, 1, nlines, mystart, myend)
      
      ! Deallocations
      IF (ALLOCATED(q)) DEALLOCATE(q)
      IF (ALLOCATED(w)) DEALLOCATE(w)
      IF (ALLOCATED(R_lines)) DEALLOCATE(R_lines)
      IF (ALLOCATED(Z_lines)) DEALLOCATE(Z_lines)
      IF (ALLOCATED(PHI_lines)) DEALLOCATE(PHI_lines)
      IF (ALLOCATED(U_lines)) DEALLOCATE(U_lines)
      IF (ALLOCATED(B_lines)) DEALLOCATE(B_lines)
      IF (ALLOCATED(goodline)) DEALLOCATE(goodline)
      
      ! Allocations
      ALLOCATE(q(neqs_nag),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'Q',ier)
      ALLOCATE(R_lines(mystart:myend,0:nsteps), &
               Z_lines(mystart:myend,0:nsteps), &
               PHI_lines(mystart:myend,0:nsteps),&
               B_lines(mystart:myend,0:nsteps), &
               U_lines(mystart:myend,0:nsteps),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R_LINES, Z_LINES PHI_LINES',ier)
      ALLOCATE(goodline(mystart:myend),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'hitwall',ier)

      ! Initializations
      R_lines=0.0
      Z_lines=0.0
      PHI_lines=-1.0
      goodline(:) = .true.
      R_lines(mystart:myend,0)   = R_start(mystart:myend)
      Z_lines(mystart:myend,0)   = Z_start(mystart:myend)
      PHI_lines(mystart:myend,0) = phi_start(mystart:myend)
      U_lines(mystart:myend,0)   = 0
      B_lines(mystart:myend,0)   = 0
         
      ! Output some stuff
      IF (lverb) THEN
         WRITE(6,'(A)')                '  ----- FOLLOWING FIELD LINES -----'
         WRITE(6,'(A,A)')              '      Method: ',TRIM(int_type)
         WRITE(6,'(A,I6)')             '       Lines: ',nlines
         WRITE(6,'(A,I7,A,ES11.4)')     '       Steps: ',nsteps,'   Delta-phi: ',dphi
         SELECT CASE(TRIM(int_type))
            CASE("NAG")
               WRITE(6,'(A,ES11.4,A,A1)')     '         Tol: ',follow_tol,'  Type: ',relab
            CASE("LSODE")
               WRITE(6,'(A,ES11.4,A,I2)')     '         Tol: ',follow_tol,'  Type: ',mf
         END SELECT
         WRITE(6,'(5X,A,I3,A)',ADVANCE='no') 'Fieldline Calculation [',0,']%'
         CALL FLUSH(6)
      END IF
      
      
      ! Follow Fieldlines
      IF (mystart <= nlines) THEN
         SELECT CASE (TRIM(int_type))
            CASE ("NAG")
#if defined(NAG)
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
#else
               ier = -1
               CALL handle_err(NAG_ERR,'torlines_follow',ier)
#endif
            CASE ("LSODE","DLSODE")
               ! Set tolerances to match NAG
               itol = 2; rtol = follow_tol; atol(:) = follow_tol
               iopt = 1 ! Optional settings used
               lrw = 20 + 16 * neqs_nag
               liw = 20
               ALLOCATE(w(lrw))
               ALLOCATE(iwork(liw))
               ier = 0
               CALL xsetf(ier)
               DO l = mystart, myend
                  ! Setup LSODE parameters
                  w = 0; iwork = 0; itask = 1; istate = 1;
                  !iwork(6) = 2000 !Maximum number of steps
                  !iwork(7) = 1    !Maximum messages per problem
                  myline  = l
                  myldex  = 0
                  ier     = 0
                  q(1)    = R_start(l)
                  q(2)    = Z_start(l)
                  phi_nag = phi_start(l)
                  phif_nag = phi_start(l)
                  dphi = SIGN(dphi,phi_end(l)-phi_start(l))  ! This should work 
                  !CALL out_torlines(phif_nag,q)
                  !dphi    = ABS(dphi)
                  !if (phi_nag > phif_nag) dphi = -dphi
                  DO 
                     CALL DLSODE(fblin_lsode,neqs_nag,q,phi_nag,phif_nag,itol,rtol,atol,itask,istate,&
                                iopt,w,lrw,iwork,liw,jacobian_lsode,mf)
                     IF (istate < -2) CALL handle_err(LSODE_ERR,'torlines_follow',istate)
                     !IF (istate == -1) istate = 2
                     CALL out_torlines(phif_nag,q)
                     iwork(11) = 0; iwork(12)=0; iwork(13)=0;
                     IF ((istate == -1) .or. (istate ==-2) .or. (ABS(phif_nag) > ABS(phi_end(l))) ) EXIT
                     !IF (myldex > nsteps) EXIT
                  END DO
               END DO
         END SELECT
      END IF

      ! Deallocations
      IF (ALLOCATED(q)) DEALLOCATE(q)
      IF (ALLOCATED(w)) DEALLOCATE(w)
      IF (ALLOCATED(iwork)) DEALLOCATE(iwork)
      
      
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE torlines_follow
