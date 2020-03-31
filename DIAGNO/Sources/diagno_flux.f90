!-----------------------------------------------------------------------
!     Module:        diagno_flux
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/02/2012
!     Description:   This subroutine calculates the flux through a
!                    simulated flux loop.
!                    The file defining the flux loop must have the
!                    following format
!  NLOOPS
!  NL     IFLFLG   IDIAFLG   TITLE1
!  X1     Y1     Z1   
!  X2     Y2     Z2           
!   .      .      .         
!   .      .      .         
!   .      .      .
!  XNL-1  YNL-1  ZNL-1
!  XNL    YNL    ZNL
!  NL     IFLFLG   IDIAFLG   TITLE2
!  X1     Y1     Z1   
!  X2     Y2     Z2           
!   .      .      .         
!   .      .      .         
!   .      .      .
!  XNL-1  YNL-1  ZNL-1
!  XNL    YNL    ZNL
!
!  WHERE
!  NLOOPS  : Total Number of flux loops in file
!  NL      : Number of elements in a loop
!  IFLFLG  : Flux Loop Flag (set to 1 if a segment should be repeated)
!  IDIAFLG : Diamagnetic Flag (set to 1 if PHIEDGE should be subtracted)
!  TITLE   : Name of Flux Loop
!  X/Y/Z   : Coordinates of points defining loop (should not close)
!
!-----------------------------------------------------------------------
      SUBROUTINE diagno_flux
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE read_wout_mod, ONLY: nfp
      USE diagno_runtime
      USE safe_open_mod
      USE biotsavart
      USE mpi_params
      USE mpi_inc
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: BYTE_8 = SELECTED_INT_KIND (8)
      INTEGER(KIND=BYTE_8),ALLOCATABLE :: mnum(:), moffsets(:)
      INTEGER(KIND=BYTE_8) :: icount, chunk
      INTEGER :: ier, i, j, k, iunit, nfl, nfl2, nsegmx, nseg, ig,&
                 i1, i2, nfl_mut, ncg, nhelp
      INTEGER, ALLOCATABLE :: nseg_ar(:), iflflg(:), idia(:), workdex(:)
      REAL(rprec) :: xp, yp, zp, xp1, yp1, zp1,int_fac
      DOUBLE PRECISION :: dtemp
      DOUBLE PRECISION, ALLOCATABLE :: flux(:)
      DOUBLE PRECISION, ALLOCATABLE :: xfl(:,:), yfl(:,:), zfl(:,:), &
                                  dx(:,:), dy(:,:), dz(:,:),&
                                  flux_mut(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: FLUX_HELP(:), MUT_HELP(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: XP_HELP(:), YP_HELP(:), ZP_HELP(:)
      DOUBLE PRECISION, ALLOCATABLE :: XP1_HELP(:), YP1_HELP(:), ZP1_HELP(:)
      CHARACTER(256) :: format_flux_headers
      CHARACTER(48), DIMENSION(:), ALLOCATABLE :: title
!      TYPE(EZspline1_r8) :: x_spl,y_spl,z_spl
!      INTEGER     ::  bcs0(2) = (/ 0, 0/)
      
!-----------------------------------------------------------------------
!     External Functions
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      INTERFACE
         REAL FUNCTION dflux_bode(x0,y0,z0,x1,y1,z1,cg)
         USE stel_kinds, ONLY: rprec
         INTEGER, intent(in)     :: cg
         REAL(rprec), intent(in) :: x0,y0,z0,x1,y1,z1
         END FUNCTION dflux_bode
         
         REAL FUNCTION dflux_midpoint(x0,y0,z0,x1,y1,z1,cg)
         USE stel_kinds, ONLY: rprec
         INTEGER, intent(in)     :: cg
         real(rprec), intent(in) :: x0,y0,z0,x1,y1,z1
         END FUNCTION dflux_midpoint
      
         REAL FUNCTION dflux_simpson(x0,y0,z0,x1,y1,z1,cg)
         USE stel_kinds, ONLY: rprec
         INTEGER, intent(in)     :: cg
         REAL(rprec), intent(in) :: x0,y0,z0,x1,y1,z1
         END FUNCTION dflux_simpson
      END INTERFACE
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Basic copy of MPI_COMM_DIANGO
      ncg = SIZE(coil_group)
      IF (ncg == 0) ncg = 2 ! Just so we don't allocate a zero size array

      if(lverb) write(6,*)' --Calculating Fluxloop Values'
      int_fac=1./int_step
      !Define format control variables
      format_flux_headers='(5X,A14,2X,I3,4X,I3,6X,I1,9X,I1,6X,F5.1)'
      IF (lverb) THEN
         WRITE(6,'(5X,A,2(9X,A,9X))') 'Loop','nseg','flux'
      END IF
      
      ! Read the diagnostics
      IF (myworkid == master) THEN
         iunit = 36
         CALL safe_open(iunit,ier,TRIM(flux_diag_file),'old','formatted')
         READ(iunit,'(1I6)') nfl
         ALLOCATE(nseg_ar(nfl), title(nfl), flux(nfl),&
                  iflflg(nfl), idia(nfl), STAT=ier)
         ALLOCATE(flux_mut(nfl,ncg), STAT=ier)

         nseg_ar = 0
         DO i = 1, nfl
            READ(iunit,'(3I6,A48)') nseg_ar(i), iflflg(i), idia(i), title(i)
            DO j = 1, nseg_ar(i)
               READ(iunit,*) xp,yp,zp
            END DO
         END DO
         nsegmx=MAXVAL(nseg_ar)+1
         ALLOCATE(xfl(nfl,0:nsegmx),yfl(nfl,0:nsegmx),zfl(nfl,0:nsegmx), STAT=ier)
         xfl = 0; yfl = 0; zfl = 0; flux=0; flux_mut=0;
         REWIND(UNIT=iunit)
         READ(iunit,'(1I6)') nfl
         DO i = 1, nfl
            READ(iunit,'(3I6,A48)') nseg_ar(i), iflflg(i), idia(i), title(i)
            DO j = 1, nseg_ar(i)
               READ(iunit,*) xfl(i,j), yfl(i,j), zfl(i,j)
            END DO
         END DO
         CLOSE(iunit)
      END IF

#if defined(MPI_OPT)
      ! Broadcast Array sizes
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_flux1',ierr_mpi)
      CALL MPI_BCAST(nfl,1,MPI_INTEGER, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_flux1',ierr_mpi)
      CALL MPI_BCAST(nsegmx,1,MPI_INTEGER, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_flux1',ierr_mpi)
    
      ! Allocate Arrays
      IF (myworkid /= master) THEN
         ALLOCATE(nseg_ar(nfl), title(nfl), flux(nfl),&
                  iflflg(nfl), idia(nfl),flux_mut(nfl,ncg), STAT=ier)
         ALLOCATE(xfl(nfl,0:nsegmx),yfl(nfl,0:nsegmx),zfl(nfl,0:nsegmx),STAT=ier)
         xfl = 0; yfl = 0; zfl = 0; flux=0; flux_mut=0;
      END IF

      ! Broadcast Arrays (1D)
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_flux2',ierr_mpi)
      CALL MPI_BCAST(nseg_ar,nfl,MPI_INTEGER, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_flux2',ierr_mpi)
      CALL MPI_BCAST(iflflg,nfl,MPI_INTEGER, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_rogow2',ierr_mpi)
      CALL MPI_BCAST(idia,nfl,MPI_INTEGER, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_rogow2',ierr_mpi)

      ! Broadcast Arrays (2D)
      CALL MPI_BCAST(xfl(1,0),nfl*(nsegmx+1),MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_rogow2',ierr_mpi)
      CALL MPI_BCAST(yfl(1,0),nfl*(nsegmx+1),MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_rogow2',ierr_mpi)
      CALL MPI_BCAST(zfl(1,0),nfl*(nsegmx+1),MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_rogow2',ierr_mpi)
#endif

      
      ! Calculate the endpoints
      ALLOCATE(dx(nfl,0:nsegmx),dy(nfl,0:nsegmx),dz(nfl,0:nsegmx),STAT=ier)
      dx=0; dy=0; dz=0
      DO i = 1, nfl
         nseg = nseg_ar(i)
         IF (iflflg(i) <= 0) THEN
            xfl(i,0) = xfl(i,nseg)
            yfl(i,0) = yfl(i,nseg)
            xfl(i,nseg+1) = xfl(i,1)
            yfl(i,nseg+1) = yfl(i,1)
         ELSE
            xfl(i,0) = xfl(i,nseg)*cos(pi2/nfp)+yfl(i,nseg)*sin(pi2/nfp)
            yfl(i,0) = yfl(i,nseg)*cos(pi2/nfp)-xfl(i,nseg)*sin(pi2/nfp)
            xfl(i,nseg+1) = xfl(i,1)*cos(pi2/nfp) - yfl(i,1)*sin(pi2/nfp)
            yfl(i,nseg+1) = yfl(i,1)*cos(pi2/nfp) + xfl(i,1)*sin(pi2/nfp)
         END IF
         zfl(i,0) = zfl(i,nseg)
         zfl(i,nseg+1) = zfl(i,1)
         dx(i,1:nseg) = xfl(i,2:nseg+1) - xfl(i,1:nseg)
         dy(i,1:nseg) = yfl(i,2:nseg+1) - yfl(i,1:nseg)
         dz(i,1:nseg) = zfl(i,2:nseg+1) - zfl(i,1:nseg)
         IF (iflflg(i) == 0) THEN
            dx(i,0) = dx(i,nseg)
            dy(i,0) = dy(i,nseg)
            dx(i,nseg+1) = dx(i,1)
            dy(i,nseg+1) = dy(i,1)
         ELSE
            dx(i,0) = dx(i,nseg)*cos(pi2/nfp)+dy(i,nseg)*sin(pi2/nfp)
            dy(i,0) = dy(i,nseg)*cos(pi2/nfp)-dx(i,nseg)*sin(pi2/nfp)
            dx(i,nseg+1) = dx(i,1)*cos(pi2/nfp) - dy(i,1)*sin(pi2/nfp)
            dy(i,nseg+1) = dy(i,1)*cos(pi2/nfp) + dx(i,1)*sin(pi2/nfp)
         END IF
         dz(i,0) = dz(i,nseg)
         dz(i,nseg+1) = dz(i,1)
      END DO

!!!!!!!!!!!!!!     NEW PART   !!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Count amount of work to do
      nhelp = 0
      DO i = 1,nfl
         IF (lskip_flux(i)) CYCLE
         nhelp = nhelp + int_step*nseg_ar(i)
      END DO

      ! Allocate the helper
      ALLOCATE(FLUX_HELP(nhelp), MUT_HELP(nhelp,ncg))
      ALLOCATE(XP_HELP(nhelp),  YP_HELP(nhelp),  ZP_HELP(nhelp))
      ALLOCATE(XP1_HELP(nhelp), YP1_HELP(nhelp), ZP1_HELP(nhelp))
      FLUX_HELP = 0; MUT_HELP = 0

      ! Calculate helpers
      i2 = 1
      DO i = 1, nfl
         IF (lskip_flux(i)) CYCLE
         DO j = 1, nseg_ar(i)
            DO k = 1, int_step
               XP_HELP(i2)  = xfl(i,j) + (k-1)*int_fac*dx(i,j)
               YP_HELP(i2)  = yfl(i,j) + (k-1)*int_fac*dy(i,j)
               ZP_HELP(i2)  = zfl(i,j) + (k-1)*int_fac*dz(i,j)
               XP1_HELP(i2) = xfl(i,j) + k*int_fac*dx(i,j)
               YP1_HELP(i2) = yfl(i,j) + k*int_fac*dy(i,j)
               ZP1_HELP(i2) = zfl(i,j) + k*int_fac*dz(i,j)
               i2 = i2 + 1
            END DO
         END DO
      END DO

      ! Divide up work
      chunk = FLOOR(REAL(nhelp) / REAL(nprocs_diagno))
      mystart = myworkid*chunk + 1
      myend = mystart + chunk - 1
#if defined(MPI_OPT)
      IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
      IF (ALLOCATED(moffsets)) DEALLOCATE(moffsets)
      ALLOCATE(mnum(nprocs_diagno), moffsets(nprocs_diagno))
      CALL MPI_ALLGATHER(chunk,1,MPI_INTEGER,mnum,1,MPI_INTEGER,MPI_COMM_DIAGNO,ierr_mpi)
      CALL MPI_ALLGATHER(mystart,1,MPI_INTEGER,moffsets,1,MPI_INTEGER,MPI_COMM_DIAGNO,ierr_mpi)
      i = 1
      DO
         IF ((moffsets(nprocs_diagno)+mnum(nprocs_diagno)-1) == i2) EXIT
         IF (i == nprocs_diagno) i = 1
         mnum(i) = mnum(i) + 1
         moffsets(i+1:nprocs_diagno) = moffsets(i+1:nprocs_diagno) + 1
         i=i+1
      END DO
      mystart = moffsets(myworkid+1)
      chunk  = mnum(myworkid+1)
      myend   = mystart + chunk - 1
      IF (myend > nhelp) myend = nhelp
#endif

      ! Loop over work (need to imlement skipping)
      DO i2 = mystart, myend
         SELECT CASE (int_type)
            CASE('midpoint')
               IF (lcoil) THEN
                  DO ig = 1, ncg
                     IF (.not.luse_extcur(ig)) CYCLE
                     MUT_HELP(i2,ig) = dflux_midpoint(XP_HELP(i2),YP_HELP(i2), ZP_HELP(i2), &
                                                   XP1_HELP(i2),YP1_HELP(i2), ZP1_HELP(i2),ig)
                  END DO
               END IF
               IF (lvac) CYCLE
               ig = 0
               FLUX_HELP(i2) = dflux_midpoint(XP_HELP(i2),YP_HELP(i2), ZP_HELP(i2), &
                                           XP1_HELP(i2),YP1_HELP(i2), ZP1_HELP(i2),ig)
            CASE('simpson')
               IF (lcoil) THEN
                  DO ig = 1, ncg
                     IF (.not.luse_extcur(ig)) CYCLE
                     MUT_HELP(i2,ig) = dflux_simpson(XP_HELP(i2),YP_HELP(i2), ZP_HELP(i2), &
                                                   XP1_HELP(i2),YP1_HELP(i2), ZP1_HELP(i2),ig)
                  END DO
               END IF
               IF (lvac) CYCLE
               ig = 0
               FLUX_HELP(i2) = dflux_simpson(XP_HELP(i2),YP_HELP(i2), ZP_HELP(i2), &
                                           XP1_HELP(i2),YP1_HELP(i2), ZP1_HELP(i2),ig)
            CASE('bode')
               IF (lcoil) THEN
                  DO ig = 1, ncg
                     IF (.not.luse_extcur(ig)) CYCLE
                     MUT_HELP(i2,ig) = dflux_bode(XP_HELP(i2),YP_HELP(i2), ZP_HELP(i2), &
                                                   XP1_HELP(i2),YP1_HELP(i2), ZP1_HELP(i2),ig)
                  END DO
               END IF
               IF (lvac) CYCLE
               ig = 0
               FLUX_HELP(i2) = dflux_bode(XP_HELP(i2),YP_HELP(i2), ZP_HELP(i2), &
                                           XP1_HELP(i2),YP1_HELP(i2), ZP1_HELP(i2),ig)
         END SELECT
      END DO
      DEALLOCATE(XP_HELP, YP_HELP, ZP_HELP, XP1_HELP, YP1_HELP, ZP1_HELP)

      ! Now handle the arrays
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_flux3',ierr_mpi)
      IF (myworkid == master) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE,FLUX_HELP,nhelp,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
      ELSE
         CALL MPI_REDUCE(FLUX_HELP,FLUX_HELP,nhelp,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
         DEALLOCATE(FLUX_HELP)
      END IF
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_flux4',ierr_mpi)
      IF (myworkid == master) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE,MUT_HELP,nhelp*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
      ELSE
         CALL MPI_REDUCE(MUT_HELP,MUT_HELP,nhelp*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
         DEALLOCATE(MUT_HELP)
      END IF
      IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
      IF (ALLOCATED(moffsets)) DEALLOCATE(moffsets)
#endif

      ! Reasemble the arrays
      i = 1
      IF (myworkid == master) THEN
         DO i2 = 1, nfl
            flux(i2) = 0
            nseg = nseg_ar(i2)
            IF (lskip_flux(i2)) CYCLE
            DO j = 1, nseg
               DO k = 1, int_step
                  flux(i2)=flux(i2) + FLUX_HELP(i)
                  DO ig = 1, ncg
                     flux_mut(i2,ig)=flux_mut(i2,ig) + MUT_HELP(i,ig)
                  END DO
                  i = i + 1
               END DO
            END DO
         END DO
         DEALLOCATE(FLUX_HELP, MUT_HELP)
      END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Get the total flux
      flux = flux + SUM(flux_mut, DIM=2)

      
      ! Write files
      IF (myworkid == master) THEN
         iunit = 30
         IF (lmut) THEN
            CALL safe_open(iunit,ier,TRIM(flux_mut_file),'replace','formatted')
            WRITE(iunit,'(2(I4,1X))') nfl,ncg
            DO i = 1, nfl
               DO ig = 1, ncg
                  WRITE(iunit,'(2(I6,1X),1ES22.12E3)') i,ig, flux_mut(i,ig)
               END DO
            END DO
            DO i = 1, nfl
               WRITE(iunit,'(I4,A)') i, title(i)
            END DO
         ELSE
            iunit = 36
            CALL safe_open(iunit,ier,'diagno_flux.'//TRIM(id_string),'replace','formatted')
            WRITE(iunit,'(1X,I4.4)') nfl
            DO i = 1, nfl
               IF (iflflg(i) == 1) flux(i) = flux(i) * nfp_diagno
               IF (idia(i) == 1) flux(i) = flux(i) + phiedge * eq_sgns
               IF (idia(i) < 0) flux(i) = flux(i) - flux(ABS(idia(i))) ! Do this so we can do flux subtraction
               flux(i) = flux(i) * flux_turns(i)
               IF (lskip_flux(i)) flux(i) = 0
               WRITE(iunit,'(1X,1ES22.12E3)') flux(i)
            END DO
            CALL FLUSH(iunit)
         END IF
         DO i = 1, nfl
            IF (lverb) WRITE(6,'(2X,i6,10X,I5,10X,E18.8)') i, nseg_ar(i), flux(i)
            WRITE(iunit,'(I4,A)') i, title(i)
         END DO
         CLOSE(iunit)
      END IF
      
      ! Clean up
      IF (ALLOCATED(flux_mut)) DEALLOCATE(flux_mut)
      DEALLOCATE(nseg_ar, title, flux, iflflg, idia)
      DEALLOCATE(xfl,yfl,zfl,dx,dy,dz)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE diagno_flux
