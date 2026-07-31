!-----------------------------------------------------------------------
!     Module:        beams3d_init_mumat
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          09/26/2012
!     Description:   This subroutine reads a soft iron or permanent
!                    magnet file, calculates the magnetic response
!                    using the MUMAT library, and adds the resultant
!                    magnetic field to our total magnetic field.
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_init_mumat
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_runtime
      USE beams3d_grid, ONLY: raxis,phiaxis,zaxis, nr, nphi, nz, &
                                 rmin, rmax, zmin, zmax, phimin, &
                                 phimax, B_R, B_Z, B_PHI, &
                                 BR4D, BPHI4D, BZ4D, &
                                 win_BR4D, win_BPHI4D, win_BZ4D, &
                                 small, eps1, eps2, eps3
      USE beams3d_physics_mod, ONLY: beams3d_BCART
      USE mumaterial_mod
      USE mpi_params  
      USE mpi_inc      
      USE mpi_sharmem
      USE EZspline
      USE EZspline_obj
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: iunit, i, j, k, s, istat, mystart, myend, ier, &
                 ourstart, ourend, debugt
      INTEGER :: bcs1(2), bcs2(2), bcs3(2)
      REAL(rprec)  :: bx_temp, by_temp, bz_temp, x_temp, &
                      y_temp, z_temp, br_temp, bphi_temp
      REAL(rprec) :: offset(3)
      INTEGER :: numprocs_local, mylocalid, mymasterid
      INTEGER :: MPI_COMM_MUSHARE, MPI_COMM_MUMASTER
      LOGICAL :: lismaster, lissubmaster
      INTEGER :: npoints_beams
      DOUBLE PRECISION, ALLOCATABLE :: x_out(:), y_out(:), z_out(:)
      DOUBLE PRECISION, ALLOCATABLE :: B_beams(:,:)
      
      TYPE(EZspline3_r8) :: BR_spl, BPHI_spl, BZ_spl
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      istat = 0; ier = 0; iunit = 327
      lismaster = .TRUE.; lissubmaster = .TRUE.
      ! Divide up Work
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_MUSHARE, ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_MUSHARE, mylocalid, ierr_mpi )              ! MPI
      CALL MPI_COMM_SIZE( MPI_COMM_MUSHARE, numprocs_local, ierr_mpi )          ! MPI
      lismaster = .FALSE.; lissubmaster = .FALSE.
      i = MPI_UNDEFINED
      IF (mylocalid.EQ.0) THEN 
        i = 0; lissubmaster = .TRUE.
      END IF
      CALL MPI_COMM_SPLIT( MPI_COMM_BEAMS, i, mylocalid, MPI_COMM_MUMASTER, ierr_mpi)

      ! Locate main master
      IF (lissubmaster) THEN
        CALL MPI_COMM_RANK( MPI_COMM_MUMASTER, mymasterid, ierr_mpi)
        lismaster = (mymasterid.EQ.0)
      END IF
#endif

    ! Set mumaterial verbosity
      CALL mumaterial_setverb(lismaster)

      ! Read the mu materials file
      CALL mumaterial_load(TRIM(mumat_string),istat, MPI_COMM_MUSHARE, MPI_COMM_MUMASTER, MPI_COMM_BEAMS)

      ! Set parameters
      CALL mumaterial_set_vars(max_error=mumaterial_tol, max_iter=mumaterial_niter, lambda_start=mumaterial_lambda, &
                           lambda_min=mumaterial_lambdamin, lambda_max=mumaterial_lambdamax, &
                           lambda_factor=mumaterial_lamfactor, min_conv_perc=mumaterial_convcheck, &
                           max_depth=INT(mumaterial_depth), max_leafsize=INT(mumaterial_leaf), &
                           iter_theta=mumaterial_theta_iter,eval_theta=mumaterial_theta_eval) 
      ! Load magnetization file
      IF (lmumat_readmag) CALL mumaterial_magfile_read(TRIM(mumat_magfile))

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_MUSHARE,  ierr_mpi)
#endif
      
      IF (lverb) CALL mumaterial_info(6, lmumat_skipiter)

      ! Create the Splines 
      IF (lverb) WRITE(6,*) "  BEAMS3D: Creating B-splines"
      IF (lissubmaster) THEN
         bcs1=(/ 0, 0/)
         bcs2=(/-1,-1/)
         bcs3=(/ 0, 0/)
         CALL EZspline_init(BR_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init_mumat:BR_spl',ier)
         CALL EZspline_init(BPHI_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init_mumat:BPHI_spl',ier)
         CALL EZspline_init(BZ_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init_mumat:BZ_spl',ier)
         BR_spl%isHermite   = 1
         BR_spl%x1   = raxis
         BR_spl%x2   = phiaxis
         BR_spl%x3   = zaxis
         BPHI_spl%isHermite = 1
         BPHI_spl%x1 = raxis
         BPHI_spl%x2 = phiaxis
         BPHI_spl%x3 = zaxis
         BZ_spl%isHermite   = 1
         BZ_spl%x1   = raxis
         BZ_spl%x2   = phiaxis
         BZ_spl%x3   = zaxis
         CALL EZspline_setup(BR_spl,B_R,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init_mumat:BR_spl',ier)
         CALL EZspline_setup(BPHI_spl,B_PHI,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init_mumat:BPHI_spl',ier)
         CALL EZspline_setup(BZ_spl,B_Z,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init_mumat:BZ_spl',ier)
      END IF
      CALL MPI_BARRIER(MPI_COMM_MUSHARE, ier)
      CALL mpialloc(BR4D,   8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_MUSHARE, win_BR4D)
      CALL mpialloc(BPHI4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_MUSHARE, win_BPHI4D)
      CALL mpialloc(BZ4D,   8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_MUSHARE, win_BZ4D)
      IF (lissubmaster) THEN
         BR4D = BR_SPL%fspl
         BPHI4D = BPHI_SPL%fspl
         BZ4D = BZ_SPL%fspl
         CALL EZspline_free(BR_spl,ier)
         CALL EZspline_free(BPHI_spl,ier)
         CALL EZspline_free(BZ_spl,ier)
      END IF
      eps1 = (rmax-rmin)*small
      eps2 = (phimax-phimin)*small
      eps3 = (zmax-zmin)*small

      ! Initialize the magnetic calculation
      offset = 0.0d0
      CALL MUMATERIAL_RUN(beams3d_BCART, offset, lmumat_skipiter, .NOT.lmumat_readmag)
      ! Output magnetics file
      IF (lmumat_writemagfile) CALL mumaterial_magfile_write(id_string)

      ! Pack coordinates
      npoints_beams = nr*nphi*nz
      ALLOCATE(x_out(npoints_beams), y_out(npoints_beams), z_out(npoints_beams))
      DO s = 1, npoints_beams
      i = MOD(s-1,nr)+1
      j = MOD(s-1,nr*nphi)/nr+1
      k = (s-1)/(nr*nphi)+1
      x_out(s) = raxis(i)*cos(phiaxis(j))
      y_out(s) = raxis(i)*sin(phiaxis(j))
      z_out(s) = zaxis(k)
      END DO

      ! Batch evaluate
      CALL mumaterial_getb_vector(x_out, y_out, z_out, B_beams, linclvac=.TRUE.)

      ! Unpack
      IF (lissubmaster) THEN
        DO s = 1, npoints_beams
          i = MOD(s-1,nr)+1
          j = MOD(s-1,nr*nphi)/nr+1
          k = (s-1)/(nr*nphi)+1
          B_R(i,j,k)   = B_beams(1,s)*cos(phiaxis(j)) + B_beams(2,s)*sin(phiaxis(j))
          B_PHI(i,j,k) = B_beams(2,s)*cos(phiaxis(j)) - B_beams(1,s)*sin(phiaxis(j))
          B_Z(i,j,k)   = B_beams(3,s)
        END DO
      END IF
      DEALLOCATE(x_out, y_out, z_out, B_beams)

      ! Clear communicators
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_MUSHARE, ierr_mpi)
      IF (lissubmaster) THEN
         CALL MPI_COMM_FREE(MPI_COMM_MUSHARE,ierr_mpi)
         CALL MPI_COMM_FREE(MPI_COMM_MUMASTER,ierr_mpi)
      END IF
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
#endif

      ! Free memory
      CALL mpidealloc(BR4D,win_BR4D)
      CALL mpidealloc(BPHI4D,win_BPHI4D)
      CALL mpidealloc(BZ4D,win_BZ4D)
      CALL MUMATERIAL_FREE()

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'beams3d_init_mumat',ierr_mpi)
#endif
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_init_mumat
