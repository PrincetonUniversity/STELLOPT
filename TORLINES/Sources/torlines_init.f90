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
      USE vmec_input, extcur_in => extcur
      USE EZspline_obj
      USE EZspline
      USE ez_hdf5
      USE wall_mod
      USE safe_open_mod
      USE mpi_params          
      USE mpi_inc
      USE mpi_sharmem
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
         WRITE(6,'(A,I3)')      '       nrho: ',nrho
         WRITE(6,'(A,I3,A,I3)') '         nu: ',nu,'   nv: ',nv
         WRITE(6,'(A,F8.3)')    '  bound_sep: ',bound_separation
         IF (lvac) WRITE(6,'(A)')      '   !!!!! Vacuum Fields Only  !!!!! '
      END IF   

      ! ALLOCATE grids 1D
      CALL mpialloc(rho, nrho, myid_sharmem, 0, MPI_COMM_SHARMEM, win_rho)
      CALL mpialloc(xu,    nu, myid_sharmem, 0, MPI_COMM_SHARMEM, win_xu)
      CALL mpialloc(xv,    nv, myid_sharmem, 0, MPI_COMM_SHARMEM, win_xv)
      ! ALLOCATE grids 3D
      CALL mpialloc(rreal,     nrho, nu, nv, myid_sharmem, 0, MPI_COMM_SHARMEM, win_rreal)
      CALL mpialloc(zreal,     nrho, nu, nv, myid_sharmem, 0, MPI_COMM_SHARMEM, win_zreal)
      CALL mpialloc(bsreal,    nrho, nu, nv, myid_sharmem, 0, MPI_COMM_SHARMEM, win_bsreal)
      CALL mpialloc(bureal,    nrho, nu, nv, myid_sharmem, 0, MPI_COMM_SHARMEM, win_bureal)
      CALL mpialloc(bvreal,    nrho, nu, nv, myid_sharmem, 0, MPI_COMM_SHARMEM, win_bvreal)
      CALL mpialloc(breal,     nrho, nu, nv, myid_sharmem, 0, MPI_COMM_SHARMEM, win_breal)
      CALL mpialloc(brreal,    nrho, nu, nv, myid_sharmem, 0, MPI_COMM_SHARMEM, win_brreal)
      CALL mpialloc(bzreal,    nrho, nu, nv, myid_sharmem, 0, MPI_COMM_SHARMEM, win_bzreal)
      CALL mpialloc(bphireal,  nrho, nu, nv, myid_sharmem, 0, MPI_COMM_SHARMEM, win_bphireal)
      !ALLOCATE(rreal(nrho,nu,nv), zreal(nrho,nu,nv), breal(nrho,nu,nu), &
      !         bsreal(nrho,nu,nv), bureal(nrho,nu,nv), bvreal(nrho,nu,nv), STAT=ier)
      !IF (ier /= 0) CALL handle_err(ALLOC_ERR,'REAL SPACE',ier)
      
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
      CALL MPI_BARRIER(MPI_COMM_SHARMEM, ierr_mpi)
      CALL torlines_write('FIELD')

      ! Allocated Shared Helpers
      CALL mpialloc(R4D,    8, nrho, nu, nv, myid_sharmem, 0, MPI_COMM_SHARMEM, win_R4D)
      CALL mpialloc(Z4D,    8, nrho, nu, nv, myid_sharmem, 0, MPI_COMM_SHARMEM, win_Z4D)
      CALL mpialloc(BXSI4D, 8, nrho, nu, nv, myid_sharmem, 0, MPI_COMM_SHARMEM, win_BXSI4D)
      CALL mpialloc(BETA4D, 8, nrho, nu, nv, myid_sharmem, 0, MPI_COMM_SHARMEM, win_BETA4D)
      CALL mpialloc(B4D,    8, nrho, nu, nv, myid_sharmem, 0, MPI_COMM_SHARMEM, win_B4D)

      ! Master Creates and Copies Splines
      IF (myid_sharmem == master) THEN
         xu   = xu*pi2
         xv   = xv*pi2/nfp_bg
         bsreal(:,:,nv) = bsreal(:,:,1)
         bureal(:,:,nv) = bureal(:,:,1)
         bvreal(:,:,nv) = bvreal(:,:,1)
         bcs1 = (/0,0/)   ! Not-a-knot
         bcs2 = (/-1,-1/) ! Periodic
         bcs3 = (/-1,-1/) ! Periodic
         ! Now setup the splines
         CALL EZspline_init(R_spl,nrho,nu,nv,bcs1,bcs2,bcs3,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/torlines_init',ier)
         CALL EZspline_init(Z_spl,nrho,nu,nv,bcs1,bcs2,bcs3,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/torlines_init',ier)
         CALL EZspline_init(bxsi_spl,nrho,nu,nv,bcs1,bcs2,bcs3,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/torlines_init',ier)
         CALL EZspline_init(beta_spl,nrho,nu,nv,bcs1,bcs2,bcs3,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/torlines_init',ier)
         CALL EZspline_init(B_spl,nrho,nu,nv,bcs1,bcs2,bcs3,ier)
         IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/torlines_init',ier)
         ! Set the parameters
         R_spl%x1    = rho
         Z_spl%x1    = rho
         bxsi_spl%x1 = rho
         beta_spl%x1 = rho
         B_spl%x1    = rho
         R_spl%x2    = xu
         Z_spl%x2    = xu
         bxsi_spl%x2 = xu
         beta_spl%x2 = xu
         B_spl%x2    = xu
         R_spl%x3    = xv
         Z_spl%x3    = xv
         bxsi_spl%x3 = xv
         beta_spl%x3 = xv
         B_spl%x3    = xv
         R_spl%isHermite    = 1
         Z_spl%isHermite    = 1
         bxsi_spl%isHermite = 1
         beta_spl%isHermite = 1
         B_spl%isHermite    = 1
         ! Create the Splines
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
         ! Copy to shared memory structures
         R4D = R_spl%fspl
         Z4D = Z_spl%fspl
         BETA4D = beta_spl%fspl
         BXSI4D = bxsi_spl%fspl
         B4D = B_spl%fspl
         ! Free the shapred memory
         CALL EZspline_free(R_spl,ier)
         CALL EZspline_free(Z_spl,ier)
         CALL EZspline_free(beta_spl,ier)
         CALL EZspline_free(bxsi_spl,ier)
         CALL EZspline_free(B_spl,ier)
      END IF

      ! These are helpers for range
      thmx = xu(nu)
      phmx = xv(nv)
      eps1 = MAXVAL(rho)*small
      eps2 = MAXVAL(xu)*small
      eps3 = MAXVAL(xv)*small

      ! DEALLOCATE Background grids
      CALL MPI_BARRIER(MPI_COMM_SHARMEM, ierr_mpi)
      CALL mpidealloc(rreal,  win_rreal)
      CALL mpidealloc(zreal,  win_zreal)
      CALL mpidealloc(bsreal, win_bsreal)
      CALL mpidealloc(bureal, win_bureal)
      CALL mpidealloc(bvreal, win_bvreal)
      CALL mpidealloc(breal,  win_breal)
      CALL mpidealloc(brreal,  win_brreal)
      CALL mpidealloc(bzreal,  win_bzreal)
      CALL mpidealloc(bphireal,  win_bphireal)

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
