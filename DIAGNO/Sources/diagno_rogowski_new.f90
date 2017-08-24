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
!      USE EZspline_obj
!      USE EZspline
      USE mpi_params
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: BYTE_8 = SELECTED_INT_KIND (8)
#if defined(MPI_OPT)
      INCLUDE 'mpif.h'
      INTEGER(KIND=BYTE_8),ALLOCATABLE :: mnum(:), moffsets(:)
#endif
      INTEGER(KIND=BYTE_8) :: chunk
      INTEGER :: ier, i, j, k, iunit, nfl, nsegmx, nseg, ig,&
                 i1, i2, nfl_mut, ncg
      INTEGER, ALLOCATABLE :: nseg_ar(:), iflflg(:), idia(:), workdex(:)
      REAL(rprec) :: xp, yp, zp, xp1, yp1, zp1, int_fac
      DOUBLE PRECISION :: dtemp
      DOUBLE PRECISION, ALLOCATABLE :: flux(:)
      DOUBLE PRECISION, ALLOCATABLE :: xfl(:,:), yfl(:,:), zfl(:,:), &
                                  dx(:,:), dy(:,:), dz(:,:),&
                                  flux_mut(:,:), eff_area(:,:)
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
      !CALL MPI_BCAST(title,nfl*48,MPI_CHARACTER, master, MPI_COMM_DIAGNO,ierr_mpi)
      !IF (ierr_mpi /=0) CALL handle_err(MPI_BCAST_ERR,'diagno_rogow2',ierr_mpi)

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

      ! Only work on correct number of flux loops
      i = COUNT(lskip_rogo(1:nfl))
      j = nfl - i
      ALLOCATE(workdex(j))
      k=1
      DO i = 1, nfl
         IF (lskip_rogo(i)) CYCLE
         workdex(k) = i
         k=k+1
      END DO
 
      ! Divide up the work
      i2 = nfl - COUNT(lskip_rogo(1:nfl))
      chunk = FLOOR(REAL(i2) / REAL(nprocs_diagno))
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
#endif

      ! Read the mutual induction file if present
      iunit = 36
      IF (luse_mut) THEN
         CALL safe_open(iunit,ier,TRIM(rog_mut_file),'old','formatted')
         READ(iunit,'(2I4)') nfl_mut, nextcur
         IF (ALLOCATED(flux_mut)) DEALLOCATE(flux_mut)
         ALLOCATE(flux_mut(nfl,nextcur), STAT=ier)
         DO i = 1, nfl
            DO ig = 1, nextcur
               READ(iunit,'(1X,I4,1X,I4,E22.14)') i1,i2,xp
               flux_mut(i1,i2) = xp
            END DO
            IF (i<mystart .or. i>myend) flux_mut(i,:)=0
         END DO
         CLOSE(iunit)
         DO ig = 1, nextcur
            IF (luse_extcur(ig)) THEN
               flux_mut(:,ig) = flux_mut(:,ig) * extcur(ig)
            ELSE
               flux_mut(:,ig) = 0
            END IF
         END DO
         ncg = nextcur+1
      END IF

      flux = 0
      ! Calculate the fields
      DO i2 = mystart, myend
         i = workdex(i2)
         flux(i) = 0
         !IF (lskip_rogo(i)) CYCLE
         nseg = nseg_ar(i)
         DO j = 1, nseg-1
            DO k = 1, int_step
               xp = xfl(i,j) + (k-1)*int_fac*dx(i,j)
               yp = yfl(i,j) + (k-1)*int_fac*dy(i,j)
               zp = zfl(i,j) + (k-1)*int_fac*dz(i,j)
               xp1 = xfl(i,j) + k*int_fac*dx(i,j)
               yp1 = yfl(i,j) + k*int_fac*dy(i,j)
               zp1 = zfl(i,j) + k*int_fac*dz(i,j)
               SELECT CASE (int_type)
                  CASE('midpoint')
                     IF (lcoil) THEN
                        DO ig = 1, ncg
                           IF (.not.luse_extcur(ig)) CYCLE
                           dtemp = db_midpoint(xp,yp,zp,xp1,yp1,zp1,ig)
                           flux_mut(i,ig)=flux_mut(i,ig) + dtemp*eff_area(i,j)
                        END DO
                     END IF
                     IF (lvac) CYCLE
                     ier = 0
                     dtemp = db_midpoint(xp,yp,zp,xp1,yp1,zp1,ier)
                     flux(i)=flux(i) + dtemp*eff_area(i,j)
                  CASE('simpson')
                     IF (lcoil) THEN
                        DO ig = 1, ncg
                           IF (.not.luse_extcur(ig)) CYCLE 
                           dtemp = db_simpson(xp,yp,zp,xp1,yp1,zp1,ig)
                           flux_mut(i,ig)=flux_mut(i,ig) + dtemp*eff_area(i,j)
                        END DO 
                     END IF
                     IF (lvac) CYCLE
                     ier = 0
                     !IF (myworkid == master) WRITE(6,*) myworkid,i,j,k,xp,yp,zp,xp1,yp1,zp1,eff_area(i,j); CALL FLUSH(6)
                     dtemp = db_simpson(xp,yp,zp,xp1,yp1,zp1,ier)
                     flux(i)=flux(i) + dtemp*eff_area(i,j)
                  CASE('bode')
                     IF (lcoil) THEN
                        DO ig = 1, ncg
                           IF (.not.luse_extcur(ig)) CYCLE
                           dtemp = db_bode(xp,yp,zp,xp1,yp1,zp1,ig)
                           flux_mut(i,ig)=flux_mut(i,ig) + dtemp*eff_area(i,j)
                        END DO
                     END IF
                     IF (lvac) CYCLE
                     ier = 0
                     dtemp = db_bode(xp,yp,zp,xp1,yp1,zp1,ier)
                     flux(i)=flux(i) + dtemp*eff_area(i,j)
               END SELECT
            END DO
         END DO
         !WRITE(6,*) i,flux(i); CALL FLUSH(6)
      END DO
      DEALLOCATE(workdex)

      ! Now handle the arrays
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_DIAGNO,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'diagno_rogo3',ierr_mpi)
      IF (lmut) THEN ! Vacuum part
         IF (myworkid == master) THEN
            CALL MPI_REDUCE(MPI_IN_PLACE,flux_mut,nfl*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
         ELSE
            CALL MPI_REDUCE(flux_mut,flux_mut,nfl*ncg,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
         END IF
      ELSE ! Plasma part
         IF (myworkid == master) THEN
            CALL MPI_REDUCE(MPI_IN_PLACE,flux,nfl,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
         ELSE
            CALL MPI_REDUCE(flux,flux,nfl,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_DIAGNO,ierr_mpi)
         END IF
      END IF
      IF (ALLOCATED(mnum)) DEALLOCATE(mnum)
      IF (ALLOCATED(moffsets)) DEALLOCATE(moffsets)
#endif

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
