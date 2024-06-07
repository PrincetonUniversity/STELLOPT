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
                                mumaterial_free, mumaterial_debug
      USE mpi_params  
      USE mpi_inc      
      USE mpi_sharmem
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
      !CALL mumaterial_debug(lismaster,lissubmaster,.TRUE.)
      CALL mumaterial_debug(.FALSE.,.FALSE.,.FALSE.)

      ! Read the mu materials file
      CALL MUMATERIAL_LOAD(TRIM(mumat_string),istat, MPI_COMM_MUSHARE, MPI_COMM_MUMASTER, MPI_COMM_BEAMS)

      ! Set parameters
      CALL MUMATERIAL_SETD(mumaterial_tol, mumaterial_niter, mumaterial_lambda, &
                           mumaterial_lamfactor, mumaterial_lamthresh, & 
                           mumaterial_padfactor,mumaterial_syncinterval) 

      

#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_MUSHARE,  ierr_mpi)
#endif
      
      IF (lverb) THEN
         CALL mumaterial_info(6)
         WRITE(6,'(A,A)') '   FILE: ',TRIM(mumat_string)
         CALL FLUSH(6)
      END IF

      ! Create the Splines 
      IF (lissubmaster) THEN
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
      offset = 0.0
      CALL MUMATERIAL_INIT_NEW(beams3d_BCART, MPI_COMM_BEAMS, MPI_COMM_MUSHARE, MPI_COMM_MUMASTER, offset)

      ! Break up the Work
      IF (lverb) WRITE(6,*) 'Calculating range'
      CALL MPI_CALC_MYRANGE(MPI_COMM_BEAMS, 1, nr*nphi*nz, mystart, myend)

      ! Find largest mystart in local
      IF (lverb) WRITE(6,*) 'Calculating ourstart and ourend'
      CALL MPI_ALLREDUCE(mystart, ourstart, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_MUSHARE, ierr_mpi)
      CALL MPI_ALLREDUCE(myend,     ourend, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_MUSHARE, ierr_mpi)
      IF (lverb) WRITE(6,*) 'Range calculated'

      ! Zero out non-work areas
      IF (lissubmaster) THEN
        IF (lverb) WRITE(6,*) 'Zeroing out non-work areas'
         DO s = 1, ourstart-1
            i = MOD(s-1,nr)+1
            j = MOD(s-1,nr*nphi)
            j = FLOOR(REAL(j) / REAL(nr))+1
            k = CEILING(REAL(s) / REAL(nr*nphi))
            B_R(i,j,k) = 0.0
            B_PHI(i,j,k) = 0.0
            B_Z(i,j,k) = 0.0
         END DO
         DO s = ourend+1, nr*nphi*nz
            i = MOD(s-1,nr)+1
            j = MOD(s-1,nr*nphi)
            j = FLOOR(REAL(j) / REAL(nr))+1
            k = CEILING(REAL(s) / REAL(nr*nphi))
            B_R(i,j,k) = 0.0
            B_PHI(i,j,k) = 0.0
            B_Z(i,j,k) = 0.0
         END DO
      END IF

#if defined(MPI_OPT)
      IF (lverb) WRITE(6,*) 'Waiting at barrier'
      CALL MPI_BARRIER(MPI_COMM_MUSHARE,ierr_mpi)
      IF (lverb) WRITE(6,*) 'Passed barrier'
#endif
      IF (lverb) WRITE(6,*) 'TEST'
      ! Start progress 
      IF (lverb) THEN
         WRITE(6,'(5X,A,I3.3,A)',ADVANCE='no') 'Magnetic Field Calculation [',0,']%'
      END IF
      CALL FLUSH(6)
      
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
         br_temp   = bx_temp*cos(phiaxis(j)) + by_temp*sin(phiaxis(j))
         bphi_temp = by_temp*cos(phiaxis(j)) - bx_temp*sin(phiaxis(j))
         B_R(i,j,k)   = B_R(i,j,k)   + br_temp
         B_PHI(i,j,k) = B_PHI(i,j,k) + bphi_temp
         B_Z(i,j,k)   = B_Z(i,j,k)   + bz_temp
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

      ! Now have master threads share results
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_MUSHARE, ierr_mpi)
      IF (lissubmaster) THEN
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, B_R,   nr*nphi*nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_MUMASTER, ierr_mpi)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, B_PHI, nr*nphi*nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_MUMASTER, ierr_mpi)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, B_Z,   nr*nphi*nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_MUMASTER, ierr_mpi)
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
      IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'beams3d_init_coil',ierr_mpi)
#endif
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE beams3d_init_mumat
