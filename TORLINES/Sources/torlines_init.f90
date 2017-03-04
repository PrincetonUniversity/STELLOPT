!-----------------------------------------------------------------------
!     Subroutine:    torlines_init
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/2/2011
!     Description:   This subroutine controls initialization of the code.
!-----------------------------------------------------------------------
      SUBROUTINE torlines_init
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE torlines_runtime
      USE torlines_input_mod
      USE torlines_fieldlines
      USE torlines_background, nfp_bg => nfp
      USE torlines_realspace
      USE torlines_runtime
      USE vmec_input, extcur_in => extcur
      USE EZspline_obj
      USE EZspline
      USE ez_hdf5
      USE wall_mod
      USE safe_open_mod
      USE mpi_params
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i,iunit, ierr
      INTEGER :: ier,u,v,ik,form,mn
      INTEGER :: bcs1(2),bcs2(2),bcs3(2)
      REAL(rprec) :: factor
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
   
      ! Handle Reading the files and initializations
      IF (lvmec) THEN
         ! Read the TORLINES namelist
         CALL safe_open(iunit,ierr,'input.' // TRIM(id_string),'old','formatted')
         CALL read_torlines_input(iunit,ierr)
         CLOSE(iunit)
         ! Read the INDATA namelist
         CALL safe_open(iunit,ierr,'input.' // TRIM(id_string),'old','formatted')
         CALL read_indata_namelist(iunit,ierr)
         CLOSE(iunit)
         ! Setup the EXTCUR ARRAY
         IF (lcoil .or. lmgrid) THEN
            IF (.not. ALLOCATED(extcur)) THEN
               DO i = 1, nigroup
                  IF (ABS(extcur_in(i)) > 0) nextcur = i
               END DO
               ALLOCATE(extcur(nextcur+1),STAT=ierr)
               IF (ierr /= 0) CALL handle_err(ALLOC_ERR,'EXTCUR',ierr)
               extcur = 0.0
               extcur(1:nextcur) = extcur_in(1:nextcur)
            END IF
         END IF
      ELSE IF (lspec) THEN
      ELSE IF (lpies) THEN
      END IF
      phimn = MAX(0.0_rprec,phi_start(1))
      
      IF (lverb) THEN
         WRITE(6,*) '-----TORLINES File Parameters-----'
         WRITE(6,'(A,I3)')      '          k: ',k
         WRITE(6,'(A,I3,A,I3)') '         nu: ',nu,'   nv: ',nv
         WRITE(6,'(A,F8.3)')    '  bound_sep: ',bound_separation
         IF (lvac) WRITE(6,'(A)')      '   !!!!! Vacuum Fields Only  !!!!! '
      END IF   
      
      ! Load vessel into memory if provided
      IF (lvessel) THEN
         CALL wall_load_txt(TRIM(vessel_string),ier)
         IF (lverb) CALL wall_info(6)
      END IF
         
      ! Now Initialize
      IF (lvmec) THEN
         CALL torlines_init_wout
      ELSE IF (lspec) THEN
      ELSE IF (lpies) THEN
      END IF
      
      ! Opent the output file so we get the grids
      CALL torlines_write('GRID_INIT')
      
      ! Now finish off the external grid
      IF (bound_separation > 1.0 .or. lvac) CALL torlines_init_external
      
      
      ! Output the fields
      CALL torlines_write('FIELD')
      
      ! Now create the splines
      ! We have BS, BU, anv BV on a background grid but we need to construct
      ! the spline quantities
      bsreal(:,:,1) = bsreal(:,:,nv)
      bureal(:,:,1) = bureal(:,:,nv)
      bvreal(:,:,1) = bvreal(:,:,nv)
      bcs1 = (/0,0/)   ! Not-a-knot
      bcs2 = (/-1,-1/) ! Periodic
      bcs3 = (/-1,-1/) ! Periodic
      ! Now setup the splines
      CALL EZspline_init(R_spl,k,nu,nv,bcs1,bcs2,bcs3,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/torlines_init',ier)
      CALL EZspline_init(Z_spl,k,nu,nv,bcs1,bcs2,bcs3,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/torlines_init',ier)
      CALL EZspline_init(bxsi_spl,k,nu,nv,bcs1,bcs2,bcs3,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/torlines_init',ier)
      CALL EZspline_init(beta_spl,k,nu,nv,bcs1,bcs2,bcs3,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/torlines_init',ier)
      CALL EZspline_init(B_spl,k,nu,nv,bcs1,bcs2,bcs3,ier)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/torlines_init',ier)
      R_spl%x1    = rho
      Z_spl%x1    = rho
      bxsi_spl%x1 = rho
      beta_spl%x1 = rho
      B_spl%x1    = rho
      R_spl%x2    = xu(1:nu)*pi2
      Z_spl%x2    = xu(1:nu)*pi2
      bxsi_spl%x2 = xu(1:nu)*pi2
      beta_spl%x2 = xu(1:nu)*pi2
      B_spl%x2    = xu(1:nu)*pi2
      R_spl%x3    = xv(1:nv)*pi2/nfp_bg
      Z_spl%x3    = xv(1:nv)*pi2/nfp_bg
      bxsi_spl%x3 = xv(1:nv)*pi2/nfp_bg
      beta_spl%x3 = xv(1:nv)*pi2/nfp_bg
      B_spl%x3   = xv(1:nv)*pi2/nfp_bg
      R_spl%isHermite    = 1
      Z_spl%isHermite    = 1
      bxsi_spl%isHermite = 1
      beta_spl%isHermite = 1
      B_spl%isHermite    = 1
      thmx = R_spl%x2(nu)
      phmx = R_spl%x3(nv)
      CALL EZspline_setup(R_spl,rreal,ier,exact_dim=.true.)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/torlines_init/R_spl',ier)
      CALL EZspline_setup(Z_spl,zreal,ier,exact_dim=.true.)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/torlines_init/Z_spl',ier)
      CALL EZspline_setup(bxsi_spl,bsreal/bvreal,ier,exact_dim=.true.)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/torlines_init/bxsi_spl',ier)
      CALL EZspline_setup(beta_spl,bureal/bvreal,ier,exact_dim=.true.)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/torlines_init/beta_spl',ier)
      CALL EZspline_setup(B_spl,breal,ier,exact_dim=.true.)
      IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/torlines_init/beta_spl',ier)
      
      
      DEALLOCATE(rreal,zreal,bsreal,bureal,bvreal,breal)
      ! Do this so we have something
      nlines = 256
      DO ik = 1, nlines
        r_start(ik) = REAL(ik-1)/(nlines-1)
        !r_start(ik) = (REAL(vsurf+1)/REAL(k))*REAL(ik-1)/(nlines-1)
        z_start(ik) = 0.0
        phi_start(ik) = 0.0
        phi_end(ik)   = 1000*phmx
      END DO
         
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE torlines_init
