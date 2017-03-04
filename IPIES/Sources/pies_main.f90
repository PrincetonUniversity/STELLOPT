!-----------------------------------------------------------------------
!     Program:       iPIES
!     Authors:       H. Greenside
!                    A. Reiman
!                    D. Monticello
!                    S. Hudson
!                    S. Lazerson
!     Date:          11/2/2011
!     Description:   PIES solves for ideal MHD equilibria with islands
!                    and stochastic regions.  Magnetic field lines are
!                    followed to create a set of magnetic surfaces.
!                    The perpendicular current is calculated from
!                    MHD force balance jxB=grad(p).  The parallel
!                    current is calculated from the divergence
!                    property of j.  div(j) = F where F=0 for ideal
!                    MHD.  This version of PIES is an attempt to
!                    interface PIES with the VMEC equilibrium code and
!                    makes use of vectorization to improve speed and
!                    performance.
!     References:
!-----------------------------------------------------------------------
      PROGRAM PIES
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE pies_runtime, ONLY: iter, lverb, ldone, lrestart, lwout, &
                                id_string, PIES_VERSION
!-----------------------------------------------------------------------
!     Local Variables
!          numargs      Number of input arguments
!          i            Index
!          arg_len      Length of input strings
!          arg1         Input file
!          args         Input arguments
!-----------------------------------------------------------------------  
      implicit none
      integer                                      :: numargs,i
      integer, parameter                           :: arg_len =256
      character*(arg_len)                          :: arg1
      character*(arg_len),allocatable,dimension(:) :: args
!-----------------------------------------------------------------------
!     Begin Program
!-----------------------------------------------------------------------
      !OPEN(6,CARRIAGECONTROL='fortran')
      WRITE(6,'(a,f4.2)') 'iPIES Version ',PIES_VERSION
      numargs=0
      i=0
      arg1=''
      lverb=.true.
      lrestart=.false.
      lwout=.false.
      ! First Handle the input arguments
      id_string = ''
      CALL GETCARG(1, arg1, numargs)
      ALLOCATE(args(numargs))
      args(1)=arg1
      ! Cycle through Arguments
      do i=1,numargs
         call GETCARG(i,args(i),numargs)
         select case (args(i))
            case ("-noverb")  ! No Verbose Output
                lverb=.false.
            case ("-restart")  ! No Verbose Output
                lrestart=.true.
            case ("-wout")  ! No Verbose Output
                lwout=.true.
            case ("-help","-h") ! Output Help message
               write(6,*)' Usage: xipies input_suffix <options>'
               write(6,*)'    input_file:  Output of VMEC or PIES'
               write(6,*)'    <options>'
               write(6,*)'     -noverb: Supress all screen output'
               write(6,*)'     -help:   Output help message'
               stop
         end select
      enddo
      ! First Argument should always be the input file
      if (numargs .gt. 0) then
         id_string = trim(arg1)
      else
         write(6,*)' Usage: xipies input_file <options>'
         write(6,*)'    input_file:  Output of VMEC or PIES'
         write(6,*)'    <options>'
         write(6,*)'     -noverb: Supress all screen output'
         write(6,*)'     -help:   Output help message'
         stop
      END IF
      DEALLOCATE(args)
      ! Initialize the Calculation
      CALL pies_init
      WRITE(6,'(a)')'---------- EXECUTION ----------'
      WRITE(6,'(a)')'ITERATION     R_AXIS    Z_AXIS   NSTEPS   BAD_SURFS IOTAMIN IOTAMAX'      
      !DO WHILE (.not. ldone)
         iter = iter + 1
         ! Follow Fieldlines and calculate magnetic surfaces
         WRITE(6,'(1X,i5,3X)',advance='no') iter
         CALL FLUSH(6)
         CALL follow
         ! Calculate the perpendicular current density
         !CALL jperp
         ! Calculate the parallel current density
         !CALL jpara
         ! Initialize quantities for B-Field solver
         !CALL solver_init
         ! Solve for magnetic Field
         !CALL b_solve
         ! Blend old and new solution
         !CALL blendb
         ! Output restart data
         !CALL write_pies_hdf5
         ! Test for convergence
         !CALL converged(ldone)
      !END DO
     
!-----------------------------------------------------------------------
!     End Program
!-----------------------------------------------------------------------
      END PROGRAM PIES
