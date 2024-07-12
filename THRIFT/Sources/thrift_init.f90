!-----------------------------------------------------------------------
!     Subroutine:    thrift_init
!     Authors:       L. van Ham
!     Date:          11/XX/2022
!     Description:   This subroutine initialzies the code for performing
!                    a run.
!-----------------------------------------------------------------------
      SUBROUTINE thrift_init
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_input_mod
      USE thrift_vars
      USE thrift_profiles_mod
      USE diagno_input_mod, ONLY: read_diagno_input
      USE safe_open_mod
      USE mpi_params
      USE mpi_inc
      USE mpi_sharmem
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        ltst        logical for supressing screen output
!        tstr1/2     String for calling paraexe
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL        :: ltst
      INTEGER        :: ier, i, iunit
      CHARACTER(256) :: tstr1,tstr2
      REAL(rprec)    :: dt
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------

      ! Read the Input Namelist for THRIFT
      IF (lverb) WRITE(6,'(A)') '----- THRIFT Input Parameters -----'
      CALL init_thrift_input

      ! Read the THRIFT input
      IF (lvmec) THEN
         CALL read_thrift_input('input.' // TRIM(id_string),ier)
      END IF

      ! Read diagno file
      IF (ldiagno) THEN
         CALL read_diagno_input('input.' // TRIM(id_string),ier)
      END IF

      ! Output to screen
      IF (lverb) THEN 
         WRITE(6,'(A)') '   FILE:             input.' // TRIM(id_string)
         WRITE(6,'(A)') '   BOOTSTRAP MODEL:  ' // TRIM(bootstrap_type)
         IF (leccd) WRITE(6,'(A)') '   ECCD MODEL:       ' // TRIM(eccd_type)
         IF (lnbcd) WRITE(6,'(A)') '   NBCD MODEL:       ' // TRIM(nbcd_type)
         WRITE(6,'(A)') ''
         WRITE(6,'(A11,I5)') '   NRHO:   ', nrho
         WRITE(6,'(A11,I5)') '   NS:     ', nsj
         WRITE(6,'(A11,I5)') '   NT:     ', ntimesteps
         WRITE(6,'(A11,F8.4)') '   TSTART: ', tstart
         WRITE(6,'(A11,F8.4)') '   TEND:   ', tend
      END IF

      ! Create the worker pool

      ! Grid allocations
      CALL mpialloc(THRIFT_RHO,       nrho,    myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_rho)
      CALL mpialloc(THRIFT_RHOFULL, nrho+2,    myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_rhofull)
      CALL mpialloc(THRIFT_S,          nsj,    myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_s)
      CALL mpialloc(THRIFT_SNOB,     nsj-2,    myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_snob)
      CALL mpialloc(THRIFT_T,   ntimesteps,    myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_t)
      CALL mpialloc(THRIFT_PHIEDGE, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_phiedge)

      ! Current densities
      CALL mpialloc(THRIFT_J,         nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_j)
      CALL mpialloc(THRIFT_JBOOT,     nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_jboot)
      CALL mpialloc(THRIFT_JPLASMA,   nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_jplasma)
      CALL mpialloc(THRIFT_JECCD,     nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_jeccd)
      CALL mpialloc(THRIFT_JNBCD,     nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_jnbcd)
      CALL mpialloc(THRIFT_JOHMIC,    nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_johmic)
      CALL mpialloc(THRIFT_JSOURCE,   nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_jsource)     
      ! Total currents
      CALL mpialloc(THRIFT_I,       nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_i)
      CALL mpialloc(THRIFT_IBOOT,   nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_iboot)
      CALL mpialloc(THRIFT_IPLASMA, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_iplasma)
      CALL mpialloc(THRIFT_IECCD,   nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_ieccd)
      CALL mpialloc(THRIFT_INBCD,   nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_inbcd)
      CALL mpialloc(THRIFT_IOHMIC,  nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_iohmic)
      CALL mpialloc(THRIFT_ISOURCE, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_isource)
      CALL mpialloc(THRIFT_UGRID,   nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_ugrid)  
      
      ! Profile variables
      CALL mpialloc(THRIFT_ETAPARA,nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_etapara)  
      CALL mpialloc(THRIFT_PPRIME ,nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_pprime)  
      CALL mpialloc(THRIFT_P      ,nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_p)  
    
      ! Magnetic variables
      CALL mpialloc(THRIFT_S11,    nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_s11)     
      CALL mpialloc(THRIFT_S12,    nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_s12)   
      CALL mpialloc(THRIFT_BAV,    nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_bav)     
      CALL mpialloc(THRIFT_BSQAV,  nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_bsqav)    
      CALL mpialloc(THRIFT_IOTA,   nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_iota)    
      CALL mpialloc(THRIFT_AMINOR, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_aminor)     
      CALL mpialloc(THRIFT_RMAJOR, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_rmajor)  
      CALL mpialloc(THRIFT_VP,     nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_vp)      
      CALL mpialloc(THRIFT_BVAV,   nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_bvav)     
      ! ABCD
      CALL mpialloc(THRIFT_COEFF_A, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_coeff_a)
      CALL mpialloc(THRIFT_COEFF_B, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_coeff_b)
      CALL mpialloc(THRIFT_COEFF_C, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_coeff_c)
      CALL mpialloc(THRIFT_COEFF_D, nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_coeff_d)
      CALL mpialloc(THRIFT_COEFF_BP,nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_coeff_bp)
      CALL mpialloc(THRIFT_COEFF_CP,nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_coeff_cp)
      CALL mpialloc(THRIFT_COEFF_DP,nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_coeff_dp)
      
      ! Alphas
      CALL mpialloc(THRIFT_ALPHA1,   nsj-2, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_alpha1)
      CALL mpialloc(THRIFT_ALPHA2,   nsj-2, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_alpha2)
      CALL mpialloc(THRIFT_ALPHA3,   nsj-2, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_alpha3)
      CALL mpialloc(THRIFT_ALPHA4,   nsj-2, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_alpha4)
      
      ! System of equations
      CALL mpialloc(THRIFT_MATLD,    nsj-1, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_matld)
      CALL mpialloc(THRIFT_MATMD,      nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_matmd)
      CALL mpialloc(THRIFT_MATUD,    nsj-1, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_matud)
      CALL mpialloc(THRIFT_MATRHS,     nsj, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_matrhs)       
      
      ! Read the Bootstrap input
      CALL tolower(bootstrap_type)
      SELECT CASE (TRIM(bootstrap_type))
         CASE('bootsj')
            ! Read BOOTSJ NAMELIST
            CALL safe_open(iunit,ier,'input.'//TRIM(id_string),'old','formatted')
            CALL read_namelist (iunit, ier, 'bootin')
            IF (ier < 0 .and. myid == master) THEN
               WRITE(6,*) '!!!!!!!!!!!!ERRROR!!!!!!!!!!!!!!'
               WRITE(6,*) '  BOOTIN Namelist not found     '
               WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               STOP
            END IF
            CLOSE(iunit)
         CASE('sfincs')
      END SELECT

      ! Now setup the profiles
      CALL read_thrift_profh5(TRIM(prof_string))

      ! Define grids
      dt = (tend-tstart)/(ntimesteps-1)
      IF (myid_sharmem == master) THEN
        FORALL(i = 1:nrho) THRIFT_RHO(i) = DBLE(i-0.5)/DBLE(nrho) ! (half) rho grid
        FORALL(i = 1:nsj)  THRIFT_S(i)   = DBLE(i-1)/DBLE(nsj-1)  ! (full)  s  grid
        FORALL(i = 1:ntimesteps) THRIFT_T(i) = tstart + (i-1)*dt  !       time grid
      END IF

      ! Extra variables (used in debugging process)
      THRIFT_RHOFULL(1) = 0.0
      THRIFT_RHOFULL(2:nrho+1) = THRIFT_RHO
      THRIFT_RHOFULL(nrho+2) = 1.0
      THRIFT_SNOB = THRIFT_S(2:nsj-1)

      ! Split off workers
      CALL thrift_init_mpisubgroup
      IF (myworkid .ne. master) THEN
         ltst  = .false.
         tstr1 = ''
         tstr2 = ''
         ier_paraexe = 0
         CALL thrift_paraexe(tstr1,tstr2,ltst)
         RETURN
      END IF
      ! - From this point on only the main thread of each run executes

      ! Initialize the equilbrium code
      IF (lvmec) THEN
         ltst = .false.
         tstr1 = 'parvmec_init'
         tstr2 = id_string
         CALL thrift_paraexe(tstr1,tstr2,ltst)
      END IF

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_init

