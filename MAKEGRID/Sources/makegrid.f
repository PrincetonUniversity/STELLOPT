!******************************************************************************
!  File makegrid.f
!*******************************************************************************
!**********************************************************************
!**                                                                  **
!**   MAIN PROGRAM:  Make Grid                                       **
!**                                                                  **
!**                                                                  **
!**     PROGRAM DESCRIPTION:                                         **
!**        MAKEGRID reads in a coils-dot file, and generates         **
!**	     an MGRID file.                                            **
!**                                                                  **
!**     Author: Steve Hirshman                                       **
!**             James Hanson                                         **
!**             Jonathan Hebert                                      **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**                                                                  **
!**********************************************************************
!-------------------------------------------------------------------------------
!   DEPENDENCIES
!-------------------------------------------------------------------------------
!
!  This file uses the following modules:
!    stel_kinds
!       located in LIBSTELL/Sources/Modules
!
!    stel_constants
!       located in LIBSTELL/Sources/Modules
!
!
!-------------------------------------------------------------------------------
!   CHANGE HISTORY
!-------------------------------------------------------------------------------
!
!  See Section V, at end of file
!
!-------------------------------------------------------------------------------
!   USAGE
!-------------------------------------------------------------------------------
!
!  Executing 'xgrid -h' will printout a help message
!
!    INPUT FILES
!
!
!   OUTPUT FILES
!      
!-------------------------------------------------------------------------------
!   COMMENTS
!-------------------------------------------------------------------------------
!
!
!*******************************************************************************
!  CODE MAKEGRID
!    
! SECTION I.	Main Program
!
! SECTION II.	Initialization Subroutines
!    interactive_input
!    namelist_input
!
! SECTION III.	TASK SUBROUTINES
!    task_mgrid
!    task_mgrid_rs
!    task_circ_tor_grid
!
! SECTION IV. AUXILIARY SUBROUTINES
!    coil_group_report
!    fft_local
!
! SECTION V.	SUBROUTINES FOR TESTNG
!
! SECTION VI.	COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************

!*******************************************************************************
! SECTION I. MAIN PROGRAM
!*******************************************************************************

      PROGRAM makegrid
!
!     THIS CODE (MAKEGRID) GENERATES B-FIELD COMPONENTS ON R,Z, PHI GRID
!     AND COMPUTES THE POLOIDAL FLUX AT NOBS OBSERVATION POINTS
!
!     NOTE TO USER: EXPERIENCE SHOWS THAT A GRID SPACING OF
!     DEL R = DEL Z <= 1/50 WORKS WELL.  HERE, DEL R = rmax-rmin.
!     GRID SPACING MUCH LARGER THAN THIS MAY ADVERSELY EFFECT CONVERGENCE
!     OF VACUUM CODE.
!
!     BOX DIMENSIONS: rmin <= R <= rmax,  zmin <= Z <= zmax
!                     KP = NO. TOROIDAL PLANES/FIELD PERIOD (MUST EQUAL VMEC VALUE)
!-------------------------------------------------------------------------------
!
!     grid dimensions:
!          ir  = no. radial (r) points in box
!          jz  = no. z points in box
!
!          suggest choosing hr == (rmax-rmin)/(ir-1) equal to
!                           hz == (zmax-zmin)/(jz-1)
!
!
!     THIS CAN BE RUN EITHER FROM THE COMMAND LINE (INTERACTIVELY) OR FROM A COMMAND FILE:
!
!     xgrid < filename.cmd
!
!     WHERE the command file (filename.cmd) has entries:
!     coils file extension
!     stell_sym (T/F)
!     rmin value
!     rmax value
!     zmin value
!     zmax value
!

!
!     SET UP GRID DIMENSIONS. USE EITHER INTERACTIVE (OR DRIVER FILE)
!     OR COMMAND LINE ARGUMENT LIST
!
!-------------------------------------------------------------------------------
!  CHANGE LOG
!-------------------------------------------------------------------------------
!     1-18-2010  JDHe 
!     Added task structure to makegrid (does not change functionality of old
!        code)
!     Added namelist input support (usable with all tasks)
!     Added help message
!     Added task to calculate vector potential on inside of vacuum vessel
!        for use with NIMROD (task CIRC_TOR_GRID)
!     Previous functionality available as called before (no command line
!        inputs prompts for mgrid data, 10 command line inputs runs mgrid with
!        those inputs) as well as through the namelist (task MGRID)
!     Added 2 arguments to command line input to reflect interactive input
!---------------------------------------------------------------------------
!  Use statements are followed by variables and subroutines used in this file
!  (subroutines are followed by parentheses)
!---------------------------------------------------------------------------

      USE write_mgrid, only: mgrid_ext, mgrid_mode, lstell_sym,                & 
     &   rmin, rmax, zmin, zmax, kp, ir, jz

      USE makegrid_global, only: task

      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: numargs, i
      CHARACTER(LEN=100) :: arg1
      CHARACTER(LEN=20) :: task_use
      
      CHARACTER(LEN=57), DIMENSION(23), PARAMETER  :: help_message=
     & (/' Program makegrid command line help                      ',
     &   '   -h    Display this message.                           ',
     &   ' Number of Command Line Arguments:                       ',
     &   '  1 (not -h):                                            ',
     &   '   Run makegrid using Namelist input file supplied in    ',
     &   '   argument                                              ',
     &   '  10:                                                    ',
     &   '   Use command line arguments as input to makegrid in    ',
     &   '   order (only for task MGRID):                          ',
     &   '     mgrid_ext                                           ',
     &   '     mgrid_mode                                          ',
     &   '     lstell_sym                                          ',
     &   '     rmin                                                ',
     &   '     rmax                                                ',
     &   '     zmin                                                ',
     &   '     zmax                                                ',
     &   '     kp                                                  ',
     &   '     ir  (optional,default=121)                          ',
     &   '     jz  (optional,default=121)                          ',
     &   '  0:                                                     ',
     &   '   Get command line input interactively (only for task   ',
     &   '   MGRID)                                                ',
     &   ' END HELP MESSAGE                                        '/)

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      
      CALL getcarg(1, arg1, numargs)

      SELECT CASE(numargs)

      CASE(1)
!        INTERACTIVE (USE REDIRECTED DRIVER FILE ALSO, XGRID < DRIVER)
         IF (arg1 .EQ. '-h' .OR. arg1 .EQ. '--h' .OR. arg1 .EQ. '-H' 
     1      .OR. arg1 .EQ. '--H') THEN
            DO i=1, SIZE(help_message)
               WRITE(6,*) help_message(i)
            END DO
            STOP
         ELSE
            CALL namelist_input_makegrid(arg1)
         END IF

      CASE(0)
!  Only used for task 'MGRID'
         CALL interactive_input() !sets task='MGRID'

      CASE(8:10)
!  Only used for task 'MGRID'
         arg1 = ADJUSTL(arg1)
         mgrid_ext = TRIM(arg1)
         CALL getcarg(2, arg1, numargs)
         IF (arg1(1:1) == 'R' .or. arg1(1:1) == 'r')
     1      mgrid_mode = 'R'
         CALL getcarg(3, arg1, numargs)
         IF (arg1(1:1) == 'Y' .or. arg1(1:1) == 'y' .or.
     1       arg1(1:1) == 'T' .or. arg1(1:1) == 't')
     1       lstell_sym = .true.
         CALL getcarg(4, arg1, numargs)
         READ (arg1, *) rmin
         CALL getcarg(5, arg1, numargs)
         READ (arg1, *) rmax
         CALL getcarg(6, arg1, numargs)
         READ (arg1, *) zmin
         CALL getcarg(7, arg1, numargs)
         READ (arg1, *) zmax
         CALL getcarg(8, arg1, numargs)
         READ (arg1, *) kp
         IF (numargs .GE. 9) THEN
            CALL getcarg(9, arg1, numargs)
            READ (arg1, *) ir
         END IF
         IF (numargs .EQ. 10) THEN
            CALL getcarg(10, arg1, numargs)
            READ (arg1, *) jz
         END IF
         task='MGRID'

      CASE DEFAULT
         STOP 'Unknown number of arguments, type "xgrid -h" for help.'
      
      END SELECT

!-----------------------------------------------------------------------
!  TASKS:
!  'makegrid' can now run two separate tasks:
!    'B_GRID':         find the magnetic field on a grid of points inside the
!                      plasma volume (for VMEC, V3FIT, etc.)
!    'CIRC_TOR_GRID':  find the vector potential on the outer-most surface
!                      (for NIMROD)
!-----------------------------------------------------------------------

!  Convert to lower case, avoid simple misinterpretations.
      task_use = task
      CALL tolower(task_use)
      
      SELECT CASE(TRIM(ADJUSTL(task_use)))
      
      CASE('mgrid')

         CALL task_mgrid()

      CASE('mgrid_rs')
      
         CALL task_mgrid_rs()

      CASE('circ_tor_mgrid')
      
         CALL task_circ_tor_grid()
          
      CASE DEFAULT
      
         WRITE(*,*) 'Unknown task ', TRIM(ADJUSTL(task_use))
      
      END SELECT

      END PROGRAM makegrid

