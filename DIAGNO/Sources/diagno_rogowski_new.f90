!-----------------------------------------------------------------------
!     Module:        diagno_rogowski_new
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/02/2012
!     Description:   This subroutine calculates the response of a
!                    simulated Rogowski loop.
!                    The file defining the Rogowski coil must have the
!                    following format
!  NLOOPS
!  NL     IFLFLG   IDIAFLG   TITLE1
!  X1     Y1     Z1    EFF_AREA1   
!  X2     Y2     Z2    EFF_AREA2           
!   .      .      .         .
!   .      .      .         .
!   .      .      .         .
!  XNL-1  YNL-1  ZNL-1 EFF_AREA NL-1
!  XNL    YNL    ZNL   EFF_AREA NL
!  NL     IFLFLG   IDIAFLG   TITLE2
!  X1     Y1     Z1    EFF_AREA1   
!  X2     Y2     Z2    EFF_AREA2           
!   .      .      .         .
!   .      .      .         .
!   .      .      .         .
!  XNL-1  YNL-1  ZNL-1 EFF_AREA NL-1
!  XNL    YNL    ZNL   EFF_AREA NL
!
!  WHERE
!  NLOOPS  : Total Number of Rogowski segements in file
!  NL      : Number of elements in a segment
!  IFLFLG  : Flux Loop Flag (set to 1 if a segment should be repeated)
!  IDIAFLG : Diamagnetic Flag (set to 1 if PHIEDGE should be subtracted)
!  TITLE   : Name of Flux Loop
!  X/Y/Z   : Coordinates of points defining loop (should not close)
!  EFF_AREA: Coefficient to multiply B.dl by (set to zero to start next
!                   segment)
!
!-----------------------------------------------------------------------
      SUBROUTINE diagno_rogowski_new
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE diagno_runtime
      USE biotsavart
      USE safe_open_mod
      USE virtual_casing_mod
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
      INTEGER(KIND=BYTE_8) :: chunk
      INTEGER :: ier, i, j, k, iunit, nfl, nsegmx, nseg, ig,&
                 i1, i2, nfl_mut, ncg, nhelp
      INTEGER, ALLOCATABLE :: nseg_ar(:), iflflg(:), idia(:), workdex(:)
      REAL(rprec) :: xp, yp, zp, xp1, yp1, zp1, int_fac
      DOUBLE PRECISION :: dtemp
      DOUBLE PRECISION, ALLOCATABLE :: flux(:)
      DOUBLE PRECISION, ALLOCATABLE :: xfl(:,:), yfl(:,:), zfl(:,:), &
                                  dx(:,:), dy(:,:), dz(:,:),&
                                  flux_mut(:,:), eff_area(:,:)
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
         REAL FUNCTION db_bode(x0,y0,z0,x1,y1,z1,cg,ldb)
         USE stel_kinds, ONLY: rprec
         INTEGER, intent(in)     :: cg
         REAL(rprec), intent(in) :: x0,y0,z0,x1,y1,z1
         LOGICAL, intent(in), optional :: ldb
         END FUNCTION db_bode
         
         REAL FUNCTION db_midpoint(x0,y0,z0,x1,y1,z1,cg,ldb)
         USE stel_kinds, ONLY: rprec
         INTEGER, intent(in)     :: cg
         real(rprec), intent(in) :: x0,y0,z0,x1,y1,z1
         LOGICAL, intent(in), optional :: ldb
         END FUNCTION db_midpoint
      
         REAL FUNCTION db_simpson(x0,y0,z0,x1,y1,z1,cg,ldb)
         USE stel_kinds, ONLY: rprec
         INTEGER, intent(in)     :: cg
         REAL(rprec), intent(in) :: x0,y0,z0,x1,y1,z1
         LOGICAL, intent(in), optional :: ldb
         END FUNCTION db_simpson
      END INTERFACE
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Basic copy of MPI_COMM_DIANGO
      ncg = SIZE(coil_group)
      IF (ncg == 0) ncg = 2 ! Just so we don't allocate a zero size array

      if(lverb) write(6,*)' --Calculating Rogowski Coil Values'
      int_fac=REAL(1)/REAL(int_step)
      !Define format control variables
      format_flux_headers='(5X,A14,2X,I3,4X,I3,6X,I1,9X,I1,6X,F5.1)'
      IF (lverb) THEN
         WRITE(6,'(5X,A,2(9X,A,9X))') 'Loop','nseg','Signal'
      END IF
      
      ! Read the diagnostics
      IF (myworkid == master) THEN
         iunit = 36
         CALL safe_open(iunit,ier,TRIM(seg_rog_file),'old','formatted')
         READ(iunit,'(1I6)') nfl
         ALLOCATE(nseg_ar(nfl), title(nfl), flux(nfl),&
                  iflflg(nfl), idia(nfl),flux_mut(nfl,ncg), STAT=ier)
         nseg_ar = 0
         DO i = 1, nfl
            READ(iunit,'(3I6,A48)') nseg_ar(i), iflflg(i), idia(i), title(i)
            DO j = 1, nseg_ar(i)
               READ(iunit,*) xp,yp,zp,xp1
            END DO
         END DO
         nsegmx=MAXVAL(nseg_ar)
         ALLOCATE(xfl(nfl,nsegmx),yfl(nfl,nsegmx),zfl(nfl,nsegmx),&
                  eff_area(nfl,nsegmx),STAT=ier)
         xfl = 0; yfl = 0; zfl = 0; eff_area=0; flux =0; flux_mut=0;
         REWIND(UNIT=iunit)
         READ(iunit,'(1I6)') nfl
         DO i = 1, nfl
            READ(iunit,'(3I6,A48)') nseg_ar(i), iflflg(i), idia(i), title(i)
            DO j = 1, nseg_ar(i)
               READ(iunit,*) xfl(i,j), yfl(i,j), zfl(i,j), eff_area(i,j)
            END DO
         END DO
         CLOSE(iunit)
      END IF

#if defined(MPI_OPT)
      ! Broadcast Array sizes
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_rogow1',ierr_mpi)
      CALL MPI_BCAST(nfl,1,MPI_INTEGER, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_rogow1',ierr_mpi)
      CALL MPI_BCAST(nsegmx,1,MPI_INTEGER, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_rogow1',ierr_mpi)
    
      ! Allocate Arrays
      IF (myworkid /= master) THEN
         ALLOCATE(nseg_ar(nfl), title(nfl), flux(nfl),&
                  iflflg(nfl), idia(nfl),flux_mut(nfl,ncg), STAT=ier)
         ALLOCATE(xfl(nfl,nsegmx),yfl(nfl,nsegmx),zfl(nfl,nsegmx),&
                  eff_area(nfl,nsegmx),STAT=ier)
         xfl = 0; yfl = 0; zfl = 0; eff_area=0; flux=0; flux_mut=0;
      END IF

      ! Broadcast Arrays (1D)
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_rogow2',ierr_mpi)
      CALL MPI_BCAST(nseg_ar, nfl, MPI_INTEGER, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_rogow2', ierr_mpi)
      CALL MPI_BCAST(iflflg,  nfl, MPI_INTEGER, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_rogow2', ierr_mpi)
      CALL MPI_BCAST(idia,    nfl, MPI_INTEGER, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_rogow2', ierr_mpi)

      ! Broadcast Arrays (2D)
      CALL MPI_BCAST(xfl,      nfl*nsegmx,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_rogow2',ierr_mpi)
      CALL MPI_BCAST(yfl,      nfl*nsegmx,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_rogow2',ierr_mpi)
      CALL MPI_BCAST(zfl,      nfl*nsegmx,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_rogow2',ierr_mpi)
      CALL MPI_BCAST(eff_area, nfl*nsegmx,MPI_DOUBLE_PRECISION, master, MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_rogow2',ierr_mpi)
#endif

      ! Calculate the endpoints
      ALLOCATE(dx(nfl,1:nsegmx),dy(nfl,1:nsegmx),dz(nfl,1:nsegmx),STAT=ier)
      dx=0; dy=0; dz=0;
      DO i = 1, nfl
         nseg = nseg_ar(i)
         dx(i,1:nseg-1) = xfl(i,2:nseg)-xfl(i,1:nseg-1)
         dy(i,1:nseg-1) = yfl(i,2:nseg)-yfl(i,1:nseg-1)
         dz(i,1:nseg-1) = zfl(i,2:nseg)-zfl(i,1:nseg-1)
      END DO


!!!!!!!!!!!!!!     NEW PART   !!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Count amount of work to do
      nhelp = 0
      DO i = 1,nfl
         IF (lskip_rogo(i)) CYCLE
         nhelp = nhelp + int_step*(nseg_ar(i)-1)
      END DO

      ! Allocate the helper
      ALLOCATE(FLUX_HELP(nhelp), MUT_HELP(nhelp,ncg))
      ALLOCATE(XP_HELP(nhelp),  YP_HELP(nhelp),  ZP_HELP(nhelp))
      ALLOCATE(XP1_HELP(nhelp), YP1_HELP(nhelp), ZP1_HELP(nhelp))
      FLUX_HELP = 0; MUT_HELP = 0

      ! Calculate helpers
      i2 = 1
      DO i = 1, nfl
         IF (lskip_rogo(i)) CYCLE
         DO j = 1, nseg_ar(i)-1
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
                     MUT_HELP(i2,ig) = db_midpoint(XP_HELP(i2),YP_HELP(i2), ZP_HELP(i2), &
                                                   XP1_HELP(i2),YP1_HELP(i2), ZP1_HELP(i2),ig)
                  END DO
               END IF
               IF (lvac) CYCLE
               ig = 0
               FLUX_HELP(i2) = db_midpoint(XP_HELP(i2),YP_HELP(i2), ZP_HELP(i2), &
                                           XP1_HELP(i2),YP1_HELP(i2), ZP1_HELP(i2),ig)
            CASE('simpson')
               IF (lcoil) THEN
                  DO ig = 1, ncg
                     IF (.not.luse_extcur(ig)) CYCLE
                     MUT_HELP(i2,ig) = db_simpson(XP_HELP(i2),YP_HELP(i2), ZP_HELP(i2), &
                                                   XP1_HELP(i2),YP1_HELP(i2), ZP1_HELP(i2),ig)
                  END DO
               END IF
               IF (lvac) CYCLE
               ig = 0
               FLUX_HELP(i2) = db_simpson(XP_HELP(i2),YP_HELP(i2), ZP_HELP(i2), &
                                           XP1_HELP(i2),YP1_HELP(i2), ZP1_HELP(i2),ig)
            CASE('bode')
               IF (lcoil) THEN
                  DO ig = 1, ncg
                     IF (.not.luse_extcur(ig)) CYCLE
                     MUT_HELP(i2,ig) = db_bode(XP_HELP(i2),YP_HELP(i2), ZP_HELP(i2), &
                                                   XP1_HELP(i2),YP1_HELP(i2), ZP1_HELP(i2),ig)
                  END DO
               END IF
               IF (lvac) CYCLE
               ig = 0
               FLUX_HELP(i2) = db_bode(XP_HELP(i2),YP_HELP(i2), ZP_HELP(i2), &
                                           XP1_HELP(i2),YP1_HELP(i2), ZP1_HELP(i2),ig)
         END SELECT
      END DO
      DEALLOCATE(XP_HELP, YP_HELP, ZP_HELP, XP1_HELP, YP1_HELP, ZP1_HELP)

      ! Now handle the arrays
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_rogo3',ierr_mpi)
      IF (myworkid == master) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE,FLUX_HELP,nhelp,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
      ELSE
         CALL MPI_REDUCE(FLUX_HELP,FLUX_HELP,nhelp,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
         DEALLOCATE(FLUX_HELP)
      END IF
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_rogo4',ierr_mpi)
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
            IF (lskip_rogo(i2)) CYCLE
            DO j = 1, nseg-1
               DO k = 1, int_step
                  flux(i2)=flux(i2) + FLUX_HELP(i)*eff_area(i2,j)
                  DO ig = 1, ncg
                     flux_mut(i2,ig)=flux_mut(i2,ig) + MUT_HELP(i,ig)*eff_area(i2,j)
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
         ier = 0
         IF (lmut) THEN
            CALL safe_open(iunit,ier,TRIM(rog_mut_file),'replace','formatted')
            WRITE(iunit,'(2(I4,1X))') nfl,ncg
            DO i = 1, nfl
               DO ig = 1, ncg
                  WRITE(iunit,'(2(I6,1X),1ES22.12E3)') i,ig, flux_mut(i,ig)
               END DO
            END DO
         ELSE
            CALL safe_open(iunit,ier,'diagno_seg.'//TRIM(id_string),'replace','formatted')
            WRITE(iunit,'(1X,I4.4)') nfl
            CALL FLUSH(iunit)
            DO i = 1, nfl
               flux(i) = flux(i) * segrog_turns(i)
               WRITE(iunit,'(1X,1ES22.12E3)') flux(i)
               CALL FLUSH(iunit)
            END DO
         END IF
         DO i = 1, nfl
            IF (lverb) WRITE(6,'(2X,i6,10X,I5,10X,E18.8)') i, nseg_ar(i), flux(i)
            WRITE(iunit,'(A)') title(i)
            CALL FLUSH(6)
            CALL FLUSH(iunit)
         END DO
         CLOSE(iunit)
      END IF
      
      ! Clean up
      DEALLOCATE(nseg_ar, title, flux, iflflg, idia, flux_mut)
      DEALLOCATE(xfl,yfl,zfl,eff_area)
      DEALLOCATE(dx,dy,dz)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE diagno_rogowski_new
