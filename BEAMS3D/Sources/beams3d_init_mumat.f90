!-----------------------------------------------------------------------
!     Module:        beams3d_init_mumat
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          09/26/2012
!     Description:   This subroutine reads a soft iron or permanent
!                    magnet file, calculates the magnetic response
!                    using the MAGTENSE library, and adds the resultant
!                    magnetic field to our total magnetic field.
!                    https://www.magtense.org/
!-----------------------------------------------------------------------
#if defined(MAGTENSE)
#endif
      SUBROUTINE beams3d_init_mumat
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_runtime
      USE beams3d_grid, ONLY: raxis,phiaxis,zaxis, nr, nphi, nz, &
                                 rmin, rmax, zmin, zmax, phimin, &
                                 phimax, B_R, B_Z, B_PHI, &
                                 BR_SPL, BPHI_SPL, BZ_SPL, &
                                 BR4D, BPHI4D, BZ4D, &
                                 win_BR4D, win_BPHI4D, win_BZ4D, &
                                 small, eps1, eps2, eps3
      USE beams3d_physics_mod, ONLY: beams3d_BCART
      USE mumaterial_mod, ONLY: mumaterial_load, mumaterial_init_new, &
                                mumaterial_info, mumaterial_getbmag_scalar,&
                                mumaterial_setverb, mumaterial_setd, &
                                mumaterial_free
      USE mpi_params  
      USE mpi_inc      
      USE mpi_sharmem
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: iunit, i, j, k, s, istat, mystart, myend, ier
      INTEGER :: bcs1(2), bcs2(2), bcs3(2)
      REAL(rprec)  :: bx_temp, by_temp, bz_temp, x_temp, &
                      y_temp, z_temp, br_temp, bphi_temp
      REAL(rprec) :: offset(3)
      INTEGER :: numprocs_local, mylocalid, mylocalmaster
      INTEGER :: MPI_COMM_LOCAL
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      istat = 0; ier = 0; iunit = 327

      ! Divide up Work
#if defined(MPI_OPT)
      CALL MPI_COMM_DUP( MPI_COMM_SHARMEM, MPI_COMM_LOCAL, ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_LOCAL, mylocalid, ierr_mpi )              ! MPI
      CALL MPI_COMM_SIZE( MPI_COMM_LOCAL, numprocs_local, ierr_mpi )          ! MPI
#endif
      
      ! Set mumaterial verbosity
      !CALL mumaterial_setverb(lverb)
      CALL mumaterial_setverb(.FALSE.)
      IF (mylocalid == 0) CALL mumaterial_setverb(.TRUE.)

      ! Read the mu materials file
      CALL mumaterial_load(TRIM(mumat_string),istat,MPI_COMM_BEAMS)

      ! Set parameters
      CALL MUMATERIAL_SETD(1.0d-5, 100, 0.7d0, 0.75d0, 100, MPI_COMM_LOCAL)

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL, ierr_mpi)
#endif
      
      IF (lverb) THEN
         CALL mumaterial_info(6)
         WRITE(6,'(A,A)') '   FILE: ',TRIM(mumat_string)
         CALL FLUSH(6)
      END IF

      ! Create the Splines 
      IF (myid_sharmem == master) THEN
         bcs1=(/ 0, 0/)
         bcs2=(/-1,-1/)
         bcs3=(/ 0, 0/)
         CALL EZspline_init(BR_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init_mumag:BR_spl',ier)
         CALL EZspline_init(BPHI_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init_mumag:BPHI_spl',ier)
         CALL EZspline_init(BZ_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ier)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init_mumag:BZ_spl',ier)
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
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init_mumag:BR_spl',ier)
         CALL EZspline_setup(BPHI_spl,B_PHI,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init_mumag:BPHI_spl',ier)
         CALL EZspline_setup(BZ_spl,B_Z,ier,EXACT_DIM=.true.)
         IF (ier /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_init_mumag:BZ_spl',ier)
      END IF
      CALL MPI_BARRIER(MPI_COMM_SHARMEM, ier)
      CALL mpialloc(BR4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_BR4D)
      CALL mpialloc(BPHI4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_BPHI4D)
      CALL mpialloc(BZ4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_BZ4D)
      IF (myid_sharmem == master) THEN
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
      offset = 0.0
      CALL MUMATERIAL_INIT_NEW(beams3d_BCART, MPI_COMM_BEAMS, offset)

      IF (lverb) THEN
         WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Magnetic Field Calculation [',0,']%'
         CALL FLUSH(6)
      END IF

      ! Break up the Work
      CALL MPI_CALC_MYRANGE(MPI_COMM_LOCAL, 1, nr*nphi*nz, mystart, myend)
      
      ! Get the fields
      DO s = mystart, myend
         i = MOD(s-1,nr)+1
         j = MOD(s-1,nr*nphi)
         j = FLOOR(REAL(j) / REAL(nr))+1
         k = CEILING(REAL(s) / REAL(nr*nphi))
         br_temp   = 0
         bphi_temp = 0
         bz_temp   = 0
         x_temp    = raxis(i)*cos(phiaxis(j))
         y_temp    = raxis(i)*sin(phiaxis(j))
         z_temp    = zaxis(k)
         CALL mumaterial_getbmag_scalar(x_temp,y_temp, z_temp, bx_temp, by_temp, bz_temp)
         br_temp = bx_temp*cos(phiaxis(j))+by_temp*sin(phiaxis(j))
         bphi_temp = by_temp*cos(phiaxis(j)) - bx_temp*sin(phiaxis(j))
         B_R(i,j,k) = B_R(i,j,k) + br_temp
         B_PHI(i,j,k) = B_PHI(i,j,k) + bphi_temp
         B_Z(i,j,k) = B_Z(i,j,k) + bz_temp
         IF (lverb .and. (MOD(s,nr) == 0)) THEN
            CALL backspace_out(6,6)
            WRITE(6,'(A,I3,A)',ADVANCE='no') '[',INT((100.*s)/(myend-mystart+1)),']%'
            CALL FLUSH(6)
         END IF
      END DO
      
      ! Clean up the progress bar
      IF (lverb) THEN
         CALL backspace_out(6,36)
         WRITE(6,'(38X)',ADVANCE='no')
         CALL backspace_out(6,36)
         WRITE(6,*)
         CALL FLUSH(6)
      END IF    

      ! Clean up
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
      CALL mpidealloc(BR4D,win_BR4D)
      CALL mpidealloc(BPHI4D,win_BPHI4D)
      CALL mpidealloc(BZ4D,win_BZ4D)
      CALL MUMATERIAL_FREE()



#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_COMM_FREE(MPI_COMM_LOCAL,ierr_mpi)
      CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'beams3d_init_coil',ierr_mpi)
#endif
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_init_mumat
