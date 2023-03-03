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
         WRITE(6,'(A)') '   ECCD MODEL:       ' // TRIM(eccd_type)
         WRITE(6,'(A)') '   NBCD MODEL:       ' // TRIM(nbcd_type)
         WRITE(6,'(A)') ''
         WRITE(6,'(A11,I5)') '   NRHO:   ', nrho
         WRITE(6,'(A11,I5)') '   NT:     ', ntimesteps
         WRITE(6,'(A11,F8.4)') '   TSTART: ', tstart
         WRITE(6,'(A11,F8.4)') '   TEND:   ', tend
      END IF

      ! Create the worker pool

      ! Allocate the Current Grid
      CALL mpialloc(THRIFT_RHO, nrho,       myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_rho)
      CALL mpialloc(THRIFT_T,   ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_t)
      CALL mpialloc(THRIFT_J,        nrho, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_j)
      CALL mpialloc(THRIFT_JBOOT,    nrho, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_jboot)
      CALL mpialloc(THRIFT_JPLASMA,  nrho, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_jplasma)
      CALL mpialloc(THRIFT_JECCD,    nrho, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_jeccd)
      CALL mpialloc(THRIFT_JNBCD,    nrho, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_jnbcd)
      CALL mpialloc(THRIFT_JOHMIC,   nrho, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_johmic)
      CALL mpialloc(THRIFT_JSOURCE,  nrho, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_jsource)
      CALL mpialloc(THRIFT_I,        nrho, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_i)
      CALL mpialloc(THRIFT_IBOOT,    nrho, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_iboot)
      CALL mpialloc(THRIFT_IPLASMA,  nrho, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_iplasma)
      CALL mpialloc(THRIFT_IECCD,    nrho, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_ieccd)
      CALL mpialloc(THRIFT_INBCD,    nrho, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_inbcd)
      CALL mpialloc(THRIFT_IOHMIC,   nrho, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_iohmic)
      CALL mpialloc(THRIFT_ISOURCE,  nrho, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_isource)
      CALL mpialloc(THRIFT_UGRID,  nrho+2,          2, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_ugrid)           
      CALL mpialloc(THRIFT_VP,     nrho+2,          2, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_vp)     
      CALL mpialloc(THRIFT_BAV,    nrho+2,          2, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_bav)     
      CALL mpialloc(THRIFT_BSQAV,  nrho+2,          2, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_bsqav)     
      CALL mpialloc(THRIFT_S11,    nrho+2,          2, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_s11)     
      CALL mpialloc(THRIFT_AMINOR, nrho+2,          2, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_aminor)     
      CALL mpialloc(THRIFT_RMAJOR, nrho+2,          2, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_rmajor)     
      CALL mpialloc(THRIFT_COEFF_A,nrho+2, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_coeff_a)
      CALL mpialloc(THRIFT_COEFF_B,nrho+2, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_coeff_b)
      CALL mpialloc(THRIFT_COEFF_C,nrho+2, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_coeff_c)
      CALL mpialloc(THRIFT_COEFF_D,nrho+2, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_coeff_d)
      CALL mpialloc(THRIFT_ALPHA1,   nrho, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_alpha1)
      CALL mpialloc(THRIFT_ALPHA2,   nrho, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_alpha2)
      CALL mpialloc(THRIFT_ALPHA3,   nrho, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_alpha3)
      CALL mpialloc(THRIFT_ALPHA4,   nrho, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_alpha4)
      CALL mpialloc(THRIFT_MATLD, nrho+2, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_matld)
      CALL mpialloc(THRIFT_MATMD, nrho+2, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_matmd)
      CALL mpialloc(THRIFT_MATUD, nrho+2, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_matud)
      CALL mpialloc(THRIFT_MATRHS,nrho+2, ntimesteps, myid_sharmem, 0, MPI_COMM_SHARMEM, win_thrift_matrhs)       
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
        FORALL(i = 1:nrho) THRIFT_RHO(i) = DBLE(i-0.5)/DBLE(nrho)
        FORALL(i = 1:ntimesteps) THRIFT_T(i) = tstart + (i-1)*dt
      END IF
      ! NB: If grids are changed to be non-uniform, change definition of 
      ! time/space steps in thrift_jinductive.f90 or things will break
      
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

