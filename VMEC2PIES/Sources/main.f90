!-----------------------------------------------------------------------
!     Program:       VMEC2PIES
!                    S. Lazerson
!     Date:          01/18/2012
!     Description:   This program creates a PIES input file from VMEC
!                    input or output files.  The code has the ability to
!                    calculate a boundary larger than the VMEC boundary
!                    in order to facilitate free boundary PIES runs.  
!                    A virtual casing principle is utilized to calculate
!                    the plasma response in the extended region.  The
!                    vacuum field is obtained from the mgrid file used
!                    for the VMEC run.  It also creates a coil_data file
!                    from a supplied coils file using the currents
!                    supplied by the mgrid file.
!
!     Execution:     In it's most basic mode the xvmec2pies executable
!                    takes a VMEC wout output file and generates a set
!                    of background coordinates with resoultions similar
!                    to those found in the VMEC equilibria.  A
!                    conversion from VMEC radial (flux) coordinates to
!                    the PIES coordinates is preformed.
!                    >xvmec2pies wout.test
!
!                    The code can also calculate PIES input file from
!                    a VMEC input file.  The resulting equilibrium has
!                    similar resolutions to those in the input file.
!                    >xvmec2pies input.test
!
!                    The user may override the VMEC radial and fourier
!                    resolutions in the VMEC equlibria via the command
!                    line
!                    >xvmec2pies wout.test -k 50 -m 8 -n 6
!                    (These are the selected modes for PIES)
!
!                    The user may also specify that a set number of
!                    external surfaces should be created.  In this case,
!                    k becomes the total number of surfaces.  The VMEC
!                    surfaces are then splined to k-extsurfs radial
!                    surfaces and extsurfs number of radial surfaces are
!                    extrapolated in real space.  Fields in the
!                    extrapolated region are obtained from the mgrid
!                    file and virtual casing.
!                    >xvmec2pies wout.test -k 50 -extsurf 5
!
!                    A coil_data file can also be produced using the
!                    EXTCUR array from the wout or input file:
!                    >xvmec2pies wout.test -c coils.test_machine
!
!                    The user may default the PIES equilibria via the
!                    command line.
!                    >xvmec2pies wout.test -fixed
!                    or
!                    >xvmec2pies wout.test -free
!                    (default is to use the value in the wout or input
!                     file)
!
!                    Note: If using and input file, the magnetic fields
!                          will be defaulted to an analytic guess based
!                          on iota or I'.  If expanding a fixed boundary
!                          VMEC equilibria, the field on the VMEC
!                          boundary will be used for all external
!                          surfaces.
!     References:
!-----------------------------------------------------------------------
      PROGRAM VMEC2PIES
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE pies_runtime, ONLY: extsurfs,lverb, id_string, m_new, n_new, &
                              free_override, VMEC2PIES_VERSION, &
                              coils_file, lmake_coils, lfreeb, mustochf,&
                              lbrho
      USE pies_background, ONLY: k
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
      extsurfs = 0
      m_new = -1
      n_new = -1
      free_override=-1
      id_string = ''
      coils_file = ''
      lmake_coils = .false.
      lbrho = .true.
      k=-1
      mustochf = 0
      ! First Handle the input arguments
      CALL GETCARG(1, arg1, numargs)
      ALLOCATE(args(numargs))
      ! Cycle through Arguments
      i=1
      DO WHILE (i .le. numargs)
         call GETCARG(i,args(i),numargs)
         select case (args(i))
            case ("-nobrho")  ! No Brho
                lbrho = .false.
            case ("-noverb")  ! No Verbose Output
                lverb=.false.
            case ("-m")
                i = i + 1
                call GETCARG(i,args(i),numargs)
                READ(args(i),'(i5)') m_new
            case ("-n")
                i = i + 1
                call GETCARG(i,args(i),numargs)
                READ(args(i),'(i5)') n_new
            case ("-k")
                i = i + 1
                call GETCARG(i,args(i),numargs)
                READ(args(i),'(i5)') k
            case ("-mustochf")
                i = i + 1
                call GETCARG(i,args(i),numargs)
                READ(args(i),'(i5)') mustochf
            case ("-c")
                i = i + 1
                call GETCARG(i,args(i),numargs)
                READ(args(i),'(a)') coils_file
                lmake_coils = .true.
            case ("-fixed")
                free_override=0
            case ("-free")
                free_override=1
            case ("-extsurf")
                i = i + 1
                call GETCARG(i,args(i),numargs)
                READ(args(i),'(i5)') extsurfs
            case ("-help","-h") ! Output Help message
               write(6,*)' VMEC2PIES - Version: ',VMEC2PIES_VERSION
               write(6,*)' Usage: xvmec2pies input_file <options>'
               write(6,*)'    input_file:  VMEC wout or input file'
               write(6,*)'    <options>'
               WRITE(6,*)'     -n ntor      Overrides VMEC Value'
               WRITE(6,*)'     -m mpol      Overrides VMEC Value'
               WRITE(6,*)'     -k surfs     Overrides VMEC value'
               WRITE(6,*)'     -fixed       Overrides VMEC value'
               WRITE(6,*)'     -free        Overrides VMEC value'
               write(6,*)'     -extsurfs i  External Surfaces'
               WRITE(6,*)'     -c coil      Creates coil_data file from coil'
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
         write(6,*)' VMEC2PIES - Version: ',VMEC2PIES_VERSION
         write(6,*)' Usage: xipies input_file <options>'
         write(6,*)'    input_file:  VMEC wout or input file'
         write(6,*)'    <options>'
         write(6,*)'     -help:   Output help message'
         stop
      END IF
      DEALLOCATE(args)
      ! Initialize the Calculation
      CALL pies_init
      CALL write_pies_input
      IF (lmake_coils) THEN
         CALL write_pies_coil
      ELSE IF (lfreeb .and. lverb) THEN
      	 WRITE(6,*)' ==================NOTICE=========================='
      	 WRITE(6,*)' =  Free Boundary PIES requires a coil_data file  ='
      	 WRITE(6,*)' =  Rerun this code with a -c file to create a    ='
      	 WRITE(6,*)' =  coil_data file for this run.  Currents will   ='
      	 WRITE(6,*)' =  be supplied by the input or wout file.        ='
      	 WRITE(6,*)' =================================================='
      END IF
!-----------------------------------------------------------------------
!     End Program
!-----------------------------------------------------------------------
      END PROGRAM VMEC2PIES
