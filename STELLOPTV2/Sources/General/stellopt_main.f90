!-----------------------------------------------------------------------
!     Program:       STELLOPT (v2.0)
!     Authors:       S. Lazerson
!     Date:          05/24/2012
!     Description:   The STELLOPT v2.0 code is a versatile optimizer
!                    which fits MHD equilibria to target parameters.
!     References:
!-----------------------------------------------------------------------
      PROGRAM STELLOPT
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE mpi_params                                                    ! MPI
!-----------------------------------------------------------------------
!     Local Variables
!          numargs      Number of input arguments
!          i            Index
!          arg_len      Length of input strings
!          arg1         Input file
!          args         Input arguments
!-----------------------------------------------------------------------
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'                                                          ! MPI
!DEC$ ENDIF  
      logical                                      :: ltst
      integer                                      :: numargs,i,ier,color,key
      integer, parameter                           :: arg_len =256
      character*(arg_len)                          :: arg1
      character*(arg_len),allocatable,dimension(:) :: args
      character(256)                               :: tstr1,tstr2
!-----------------------------------------------------------------------
!     Begin Program
!-----------------------------------------------------------------------
      
      myid = master
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_INIT( ierr_mpi )                                         ! MPI
      color = 0
      CALL MPI_COMM_SPLIT( MPI_COMM_WORLD,color,myid,MPI_COMM_STEL,ierr_mpi)
      CALL MPI_COMM_RANK( MPI_COMM_STEL, myid, ierr_mpi )              ! MPI
      CALL MPI_COMM_SIZE( MPI_COMM_STEL, numprocs, ierr_mpi )          ! MPI
      CALL MPI_ERRHANDLER_SET(MPI_COMM_WORLD,MPI_ERRORS_RETURN,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_main',ierr_mpi)
!DEC$ ENDIF
      write(0,*)'stel_main. myid,numprocs=',myid,numprocs  !hm-12/5/18. (7l23h)
      pi = 4.0 * ATAN(1.0)
      pi2 = 8.0 * ATAN(1.0)
      mu0 = 16.0E-7 * ATAN(1.0)
      lneed_output = .false.
      lrestart = .false.
      ltriangulate = .false.
      lno_restart = .false.
      lauto_domain = .false.
      pct_domain = 0.05
      xvec_file = 'xvec.dat'
      INQUIRE(UNIT=6,NAME=screen_str) ! Store STDOUT
      IF (myid == master) THEN
         !OPEN(6,CARRIAGECONTROL='fortran')
         numargs=0
         i=0
         arg1=''
         lverb    = .true.
         id_string     = ''
         
         ! First Handle the input arguments
         CALL GETCARG(1, arg1, numargs)
         ALLOCATE(args(numargs))
         ! Get filename
         CALL GETCARG(1,id_string,numargs)
         id_string=TRIM(id_string)
         id_string=ADJUSTL(id_string)
         id_string=TRIM(id_string)
         ! Cycle through Arguments
         i=2
         DO WHILE (i <= numargs)
            call GETCARG(i,args(i),numargs)
            select case (args(i))
               case ("-restart") ! Use restart file in directory
                  lrestart = .true.
               case ("-noverb")  ! No Verbose Output
                  lverb=.false.
               case ("-log")
                  screen_str = "log."//id_string(7:)
                  CLOSE(UNIT=6)
                  OPEN(UNIT=6,FILE=TRIM(screen_str))
               case ("-autodomain")  ! No Verbose Output
                  lauto_domain = .true.
                  i=i+1
                  call GETCARG(i,args(i),numargs)
                  READ(args(i),*) pct_domain
               case ("-tri")
                  i=i+1
                  call GETCARG(i,args(i),numargs)
                  i=i+1
                  call GETCARG(i,args(i),numargs)
                  CALL stellopt_triangulate(args(i-1),args(i))
                  ltriangulate = .true.
               case ("-xvec_file")
                  i=i+1
                  call GETCARG(i,args(i),numargs)
                  xvec_file = args(i)
               case ("-help","-h") ! Output Help message
                  write(6,*)' STELLOPT Optimizer'
                  write(6,*)' Usage: xstellopt input_file <options>'
                  write(6,*)'    <options>'
                  write(6,*)'     -restart          Restart a run from reset file'
                  write(6,*)'     -autodomain pct   Automatically calculate min-max domain'
                  write(6,*)'     -noverb           Supress all screen output'
                  write(6,*)'     -log              Output screen to log file'
                  write(6,*)'     -tri file1 file2  Triangulation files'
                  write(6,*)'     -xvec_file file   X_VEC filename (OPT_TYPE: EVAL_XVEC)'
                  write(6,*)'     -help:            Output help message'
!DEC$ IF DEFINED (MPI_OPT)
                  CALL MPI_FINALIZE(ierr_mpi)   
                  IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_FINE_ERR,'stellot_main',ierr_mpi)
!DEC$ ENDIF
            end select
            i = i + 1
         END DO
         DEALLOCATE(args)
         WRITE(6,'(a,f5.2)') 'STELLOPT Version ',STELLOPT_VERSION
      ELSE 
         lverb=.false.   ! Shutup the slaves
      END IF
      ! Broadcast variables
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BCAST(lrestart,1,MPI_LOGICAL, master, MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'stellopt_main',ierr_mpi)
      CALL MPI_BCAST(ltriangulate,1,MPI_LOGICAL, master, MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'stellopt_main',ierr_mpi)
      CALL MPI_BCAST(lauto_domain,1,MPI_LOGICAL, master, MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'stellopt_main',ierr_mpi)
      CALL MPI_BCAST(pct_domain,1,MPI_REAL8, master, MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'stellopt_main',ierr_mpi)
      CALL MPI_BCAST(id_string,256,MPI_CHARACTER, master, MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'stellopt_main',ierr_mpi)
      CALL MPI_BARRIER( MPI_COMM_STEL, ierr_mpi )                   ! MPI
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'stellopt_main',ierr_mpi)
!DEC$ ENDIF

      write(0,*)'hm-9/12/18-1.stellopt_main. unit=0.'
      ! Initialize the Calculation
      CALL stellopt_init

!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER( MPI_COMM_STEL, ierr_mpi )                   ! MPI
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'stellopt_main',ierr_mpi)

      CALL FLUSH(6)
      ! Now we create workers
      IF ((noptimizers == -1) .and. lverb)  THEN
         WRITE(6,*) '   ========Parallel Code Execution Info======='
         WRITE(6,*) '   Number of Processors:            ',numprocs
         WRITE(6,*) '   Number of Optimization Threads:  ',noptimizers
         WRITE(6,*) '   Workers per optimizer thread:      ',numprocs/noptimizers
      END IF
      IF (noptimizers == -1) noptimizers = numprocs + 1
      color = MOD(myid,noptimizers)
      key = myid
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, color, key,MPI_COMM_MYWORLD, ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_main',ierr_mpi)
      CALL MPI_COMM_RANK(MPI_COMM_MYWORLD,myworkid,ierr_mpi)

      ! Now we need to define MPI_COMM_STEL
      CALL MPI_COMM_FREE(MPI_COMM_STEL, ierr_mpi)
      IF (myworkid /= master) THEN
         myid = -1 ! they're not part of MPI_COMM_STEL
         color = MPI_UNDEFINED
      ELSE
         color = 0
      END IF
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, color,key,MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_main',ierr_mpi)
      IF (myworkid == master)THEN
         CALL MPI_COMM_RANK( MPI_COMM_STEL, myid, ierr_mpi )              ! MPI
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_main',ierr_mpi)
         CALL MPI_COMM_SIZE( MPI_COMM_STEL, numprocs, ierr_mpi )          ! MPI
         IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_main',ierr_mpi)
      END IF
      IF (lverb)  WRITE(6,*) '   Number of Optimizer Threads:    ',numprocs
!DEC$ ENDIF

      ! Run optimization
      IF (myworkid .ne. master) THEN
         ltst  = .false.
         tstr1 = ''
         tstr2 = ''
         ier_paraexe = 0
         CALL stellopt_paraexe(tstr1,tstr2,ltst)
      ELSE
         CALL stellopt_optimize
         ltst  = .false.
         tstr1 = 'exit'
         tstr2 = ''
         CALL stellopt_paraexe(tstr1,tstr2,ltst)
      END IF

      ! Clean up
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_COMM_FREE(MPI_COMM_STEL, ierr_mpi)
      CALL MPI_FINALIZE(ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_FINE_ERR,'stellot_main',ierr_mpi)
!DEC$ ENDIF
      IF (lverb) WRITE(6,'(A)')'----- STELLOPT DONE -----'
     
!-----------------------------------------------------------------------
!     End Program
!-----------------------------------------------------------------------
      END PROGRAM STELLOPT
