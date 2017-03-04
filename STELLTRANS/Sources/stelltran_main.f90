!-----------------------------------------------------------------------
!     Program:       STELLTRAN (v0.01)
!     Authors:       S. Lazerson
!     Date:          08/21/2012
!     Description:   The is STELLTRAN, a stellarator transport code.
!     References:
!-----------------------------------------------------------------------
      PROGRAM STELLTRAN
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stelltran_runtime
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
      pi = 4.0 * ATAN(1.0)
      pi2 = 8.0 * ATAN(1.0)
      mu0 = 16.0E-7 * ATAN(1.0)
      lrestart = .false.
      IF (myid == master) THEN
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
               case ("-restart") ! restart a run
                  lrestart = .true.
               case ("-noverb")  ! No Verbose Output
                  lverb=.false.
               case ("-help","-h") ! Output Help message
                  write(6,*)' STELTRAN 3D Transport Code'
                  write(6,*)' Usage: xstelltran input_file <options>'
                  write(6,*)'    <options>'
                  write(6,*)'     -restart:      Restart a run'
                  write(6,*)'     -noverb:       Supress all screen output'
                  write(6,*)'     -help:         Output help message'
!DEC$ IF DEFINED (MPI_OPT)
                  CALL MPI_FINALIZE(ierr_mpi)   
                  IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_FINE_ERR,'stellot_main',ierr_mpi)
!DEC$ ENDIF
            end select
            i = i + 1
         END DO
         DEALLOCATE(args)
         WRITE(6,'(a,f5.2)') 'STELLTRAN Version ',STELLTRAN_VERSION
      ELSE 
         lverb=.false.   ! Shutup the slaves
      END IF
      ! Broadcast variables
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BCAST(lrestart,1,MPI_LOGICAL, master, MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'stellot_main',ierr_mpi)
      CALL MPI_BCAST(id_string,256,MPI_CHARACTER, master, MPI_COMM_STEL,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BCAST_ERR,'stellot_main',ierr_mpi)
      CALL MPI_BARRIER( MPI_COMM_STEL, ierr_mpi )                   ! MPI
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'stellot_init',ierr_mpi)
!DEC$ ENDIF
      ! Initialize the Calculation
      CALL stelltran_init

!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER( MPI_COMM_STEL, ierr_mpi )                   ! MPI
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_BARRIER_ERR,'stellot_init',ierr_mpi)

      CALL FLUSH(6)
      ! Now we create workers
      IF ((npopulation == -1) .and. lverb)  THEN
         WRITE(6,*) '   ========Parallel Code Execution Info======='
         WRITE(6,*) '   Number of Processors:            ',numprocs
         WRITE(6,*) '   Number of Optimization Threads:  ',npopulation
         WRITE(6,*) '   Workers per optimizer thread:      ',numprocs/npopulation
      END IF
      npopulation = -1
      IF (npopulation == -1) npopulation = numprocs + 1
      color = MOD(myid,npopulation)
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
         !CALL stelltran_paraexe(tstr1,tstr2,ltst)
      ELSE
         CALL stelltran_cycle
         ltst  = .false.
         tstr1 = 'exit'
         tstr2 = ''
         !CALL stelltran_paraexe(tstr1,tstr2,ltst)
      END IF
      ! Clean up
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_COMM_FREE(MPI_COMM_STEL, ierr_mpi)
      CALL MPI_FINALIZE(ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_FINE_ERR,'stellot_main',ierr_mpi)
!DEC$ ENDIF
      IF (lverb) WRITE(6,'(A)')'----- STELLTRAN DONE -----'
     
!-----------------------------------------------------------------------
!     End Program
!-----------------------------------------------------------------------
      END PROGRAM STELLTRAN
