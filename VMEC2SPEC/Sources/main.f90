!-----------------------------------------------------------------------
!     Program:       VMEC2SPEC
!                    S. Lazerson
!     Date:          04/06/2012
!     Description:   This program creates a SPEC input file from VMEC
!                    input or output files.  The code reads an
!                    additional input name list 'SPEC_INDATA' from the
!                    VMEC input file.  This namelist controls the
!                    creation of the SPEC interfaces.
!
!     Execution:     The xvmec2spec executable takes the VMEC runtime
!                    extension as an input parameters, reads the
!                    VMEC input file for the INDATA and SPEC_INDATA
!                    namelists and creates a SPEC input file.
!                    >xvmec2pies test
!
!                    The code may construct the interfaces from a
!                    VMEC output file by specifying the -wout
!                    parameter on the command line.
!                    >xvmec2spec test -wout
!
!                    Note: If using and input file, the magnetic fields
!                          will be defaulted to an analytic guess based
!                          on iota or I'.
!     References:
!-----------------------------------------------------------------------
      PROGRAM VMEC2SPEC
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE spec_runtime, ONLY: lverb, lwout, id_string, m_new,n_new, &
                              VMEC2SPEC_VERSION, lflipped
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
      numargs=0
      i=0
      arg1=''
      lverb=.true.
      id_string = ''
      m_new = 0
      n_new = -1
      lflipped = .false.
      ! First Handle the input arguments
      CALL GETCARG(1, arg1, numargs)
      ALLOCATE(args(numargs))
      ! Cycle through Arguments
      i=1
      DO WHILE (i .le. numargs)
         call GETCARG(i,args(i),numargs)
         select case (args(i))
            case ("-noverb")  ! No Verbose Output
                lverb=.false.
            case ("-wout")
                lwout=.true.
            case ("-flip")
                lflipped=.true.
            case ("-help","-h") ! Output Help message
               write(6,*)' VMEC2SPEC - Version: ',VMEC2SPEC_VERSION
               write(6,*)' Usage: xvmec2spec input_file <options>'
               write(6,*)'    input_file:  VMEC input file'
               write(6,*)'    <options>'
               WRITE(6,*)'     -wout        Use VMEC Equilibria'
               WRITE(6,*)'     -flip        Flip Toroidal Harmonics'
               write(6,*)'     -noverb      Supress all screen output'
               write(6,*)'     -help        Output help message'
               stop
         end select
         i = i + 1
      END DO
      ! First Argument should always be the input file
      if (numargs .gt. 0) then
         id_string = trim(arg1)
      else
         write(6,*)' VMEC2SPEC - Version: ',VMEC2SPEC_VERSION
         write(6,*)' Usage: xvmec2spec input_file <options>'
         write(6,*)'    input_file:  VMEC input suffix'
         write(6,*)'    <options>'
         write(6,*)'     -help:   Output help message'
         stop
      END IF
      DEALLOCATE(args)
      ! Initialize the Calculation
      CALL spec_init
      CALL write_spec_input
!-----------------------------------------------------------------------
!     End Program
!-----------------------------------------------------------------------
      END PROGRAM VMEC2SPEC
