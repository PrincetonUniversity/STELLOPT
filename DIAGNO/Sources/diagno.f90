!-----------------------------------------------------------------------
!     Program:       DIAGNO
!     Authors:       S. Lazerson (PPPL)
!     Date:          07/23/2012
!     Description:   The DIAGNO code calculates sythetic magnetic
!                    diagnostics for a given equilibrium using a
!                    virtual casing principle or a volume integral over
!                    the plasma current density.
!     References:    S.A. Lazerson, PPCF, Vol. 55, 025014 (2013)
!                    S.A. Lazerson, PPCF, Vol. 54, 122002 (2012)
!                    H.J. Gardner, Nuclear Fusion, Vol. 30, No. 8 (1990)
!                    V.D. Shfranov and L.E. Zakharov, Nuclear Fusion,
!                         Vol. 12 (1972)
!                    L.E. Zakharov, Nuclear Fusion, Vol. 13 (1973)
!                    M. Drevlak et al., Nuclear Fusion, Vol. 45 (2005)
!                    J.D. Hanson and S.P. Hirshman, Phys. Plas.,
!                         Vol 9, No. 10 (2002)
!                    
!-----------------------------------------------------------------------
      program diagno
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE diagno_runtime
      USE diagno_input_mod
      USE mpi_params
      USE virtual_casing_mod, ONLY: free_virtual_casing
!-----------------------------------------------------------------------
!     Local Variables
!          numargs      Number of input arguments
!          iseq         Index
!          arg_len      Length of input strings
!          arg1         Input file
!          args         Input arguments
!-----------------------------------------------------------------------   
      implicit none
#if defined(MPI_OPT)
      INCLUDE 'mpif.h'
#endif
      logical :: lbench
      integer numargs,i,ier
      integer, parameter :: arg_len =256
      character*(arg_len) :: arg1
      character*(arg_len),allocatable,dimension(:) :: args
!-----------------------------------------------------------------------
!     Begin Program
!-----------------------------------------------------------------------
    myworkid = master
#if defined(MPI_OPT)
    CALL MPI_INIT(ierr_mpi)
    CALL MPI_COMM_DUP( MPI_COMM_WORLD, MPI_COMM_DIAGNO, ierr_mpi)
    CALL MPI_COMM_RANK(MPI_COMM_DIAGNO, myworkid, ierr_mpi)
    CALL MPI_COMM_SIZE(MPI_COMM_DIAGNO, nprocs_diagno, ierr_mpi)
    CALL MPI_COMM_SPLIT_TYPE(MPI_COMM_DIAGNO, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, MPI_COMM_SHARMEM, ierr_mpi)
    CALL MPI_COMM_RANK(MPI_COMM_SHARMEM, myid_sharmem, ierr_mpi)
#endif
    IF (myworkid == master) THEN
       numargs=0
       arg1=''
       lbench  = .false.
       lverb   = .true.
       lvmec   = .false.
       lpies   = .false.
       lspec   = .false.
       lcoil   = .false.
       lvac    = .false.
       lmut    = .false.
       lskip_flux  = .FALSE.
       lskip_rogo  = .FALSE.
      
       ! First Handle the input arguments
       CALL GETCARG(1, arg1, numargs)
       ALLOCATE(args(numargs))
       ! Cycle through Arguments
       i=1
       DO WHILE (i <= numargs)
          call GETCARG(i,args(i),numargs)
          select case (args(i))
             case ("-noverb")  ! No Verbose Output
                lverb=.false.
             case ("-mutual")  ! Calc Mutual Induction
                lmut=.true.
             case ("-bench")  ! Calc Mutual Induction
                lbench=.true.
             case ("-vac")  ! Vacuum Fields only
                lvac=.true.
             case ("-vmec")
                i = i + 1
                lvmec = .true.
                lpies = .false.
                lspec = .false.
                CALL GETCARG(i,id_string,numargs)
             case ("-pies")
                i = i + 1
                lpies = .true.
                lvmec = .false.
                lspec = .false.
                CALL GETCARG(i,id_string,numargs)
             case ("-spec")
                i = i + 1
                lspec = .true.
                lpies = .false.
                lvmec = .false.
                CALL GETCARG(i,id_string,numargs)
             case ("-coil")
                i = i + 1
                lcoil  = .true.
                CALL GETCARG(i,coil_string,numargs)
             case ("-help","-h") ! Output Help message
               write(6,'(A,F5.2,A)')' Magnetic Diagnostic Code (v.',DIAGNO_VERSION,')'
               write(6,*)' Usage: xdiagno <options>'
               write(6,*)'    <options>'
               write(6,*)'     -vmec ext:    VMEC input/wout extension'
               !write(6,*)'     -pies ext:    PIES input extension (must have &INDATA namelist)'
               !write(6,*)'     -spec ext:    SPEC input extension (must have &INDATA namelist)'
               write(6,*)'     -coil file:   Coils. File'
               write(6,*)'     -vac:         Vacuum Field Only'
               write(6,*)'     -mutual:      Mutual Induction Calc'
               write(6,*)'     -noverb:      Supress all screen output'
               write(6,*)'     -help:        Output help message'
               stop
         end select
         i = i + 1
       END DO
       DEALLOCATE(args)
       id_string = TRIM(id_string)
       id_string = ADJUSTL(id_string)
       coil_string = TRIM(coil_string)
       coil_string = ADJUSTL(coil_string)
    ELSE
       lverb = .false. ! Shutup the workers
       lskip_flux  = .FALSE.
       lskip_rogo  = .FALSE.
    END IF

#if defined(MPI_OPT)
    CALL MPI_BCAST(id_string, 256, MPI_CHARACTER, master, MPI_COMM_DIAGNO, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'diagno_main', ierr_mpi)
    CALL MPI_BCAST(coil_string, 256, MPI_CHARACTER, master, MPI_COMM_DIAGNO, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'diagno_main', ierr_mpi)
    CALL MPI_BCAST(lbench, 1, MPI_LOGICAL, master, MPI_COMM_DIAGNO, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'diagno_main', ierr_mpi)
    CALL MPI_BCAST(lvmec, 1, MPI_LOGICAL, master, MPI_COMM_DIAGNO, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'diagno_main', ierr_mpi)
    CALL MPI_BCAST(lpies, 1, MPI_LOGICAL, master, MPI_COMM_DIAGNO, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'diagno_main', ierr_mpi)
    CALL MPI_BCAST(lspec, 1, MPI_LOGICAL, master, MPI_COMM_DIAGNO, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'diagno_main', ierr_mpi)
    CALL MPI_BCAST(lcoil, 1, MPI_LOGICAL, master, MPI_COMM_DIAGNO, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'diagno_main', ierr_mpi)
    CALL MPI_BCAST(lvac, 1, MPI_LOGICAL, master, MPI_COMM_DIAGNO, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'diagno_main', ierr_mpi)
    CALL MPI_BCAST(lmut, 1, MPI_LOGICAL, master, MPI_COMM_DIAGNO, ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR, 'diagno_main', ierr_mpi)
#endif
      
      ! Read the Input file
      if(lverb) write(6,'(A)')'============================================'
      if(lverb) write(6,'(A,F5.2,A)')'=========  D I A G N O  (v.',DIAGNO_VERSION,')  ========='
      if(lverb) write(6,'(A)')' - Reading input file'
      call read_diagno_input(id_string,ier)        !-----  Read diagno.control File -----
      
      IF (lbench) THEN
         CALL diagno_bench
      ELSE
         
         ! Read the equilibrium file
         IF (lvmec .and. .not. lvac) THEN
            if(lverb) write(6,'(A)')' - Reading equilibrium surface'
            CALL diagno_init_vmec
         ELSE IF (lpies .and. .not. lvac) THEN
            if(lverb) write(6,'(A)')' - Reading equilibrium surface'
            !CALL diagno_init_pies
            stop '!!!!!!! PIES Support pending !!!!!!'
         ELSE IF (lspec .and. .not. lvac) THEN
            if(lverb) write(6,'(A)')' - Reading equilibrium surface'
            CALL diagno_init_spec
         END IF
      
         ! Read the Coils file
         IF (lcoil) CALL diagno_init_coil
      
         ! Calculate diagnostic response
         if(lverb) write(6,*)' - Calculating diagnostic responses'
         IF (LEN_TRIM(bfield_points_file) > 1) CALL diagno_bfield
         IF (LEN_TRIM(bprobes_file) > 1) CALL diagno_bprobes
         IF (LEN_TRIM(mirnov_file) > 1) CALL diagno_mirnov
         IF (LEN_TRIM(seg_rog_file) > 1) CALL diagno_rogowski_new
         IF (LEN_TRIM(flux_diag_file) > 1) CALL diagno_flux

      END IF
      
      
      ! Do Cleanup
      IF (lvmec) THEN
         CALL free_virtual_casing(MPI_COMM_SHARMEM)
      ELSE IF (lpies) THEN
      ELSE IF (lspec) THEN
      ELSE IF (lcoil) THEN
      ELSE
      END IF
      if(lverb) write(6,*)'============  DIAGNO Complete  ==========='
      if(lverb) write(6,*)'=========================================='
#if defined(MPI_OPT)
    ierr_mpi=0
    CALL MPI_FINALIZE(ierr_mpi)
    IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_FINE_ERR, 'diagno_main', ierr_mpi)
#endif
!-----------------------------------------------------------------------
!     End Program
!-----------------------------------------------------------------------
      end program diagno
