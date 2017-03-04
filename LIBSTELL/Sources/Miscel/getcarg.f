      SUBROUTINE getcarg(narg, arg, numargs)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in)  :: narg
      INTEGER, INTENT(out) :: numargs
      CHARACTER(LEN=*), INTENT(out) :: arg
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: numchars
C-----------------------------------------------
!DEC$ IF DEFINED (WIN32)
      INTEGER :: nargs
      numargs = nargs() - 1
      CALL getarg(narg, arg, numchars)
!DEC$ ELSEIF DEFINED (LINUX)
      INTEGER iargc
      numargs = iargc()
      CALL getarg(narg, arg)
!DEC$ ELSEIF DEFINED (VMS)
      CALL lib$get_foreign(arg,,numchars)
      numargs = MIN(1,numchars)
!DEC$ ELSEIF DEFINED (CRAY)
      INTEGER :: ier
      INTEGER ipxfargc
      numargs = ipxfargc()
      CALL pxfgetarg(narg, arg, numchars, ier)
!DEC$ ELSEIF DEFINED (MACOSX)
      numargs = command_argument_count()
      CALL get_command_argument(narg,arg)
!DEC$ ELSE
      INTEGER iargc, getarg
      numargs = iargc()
      numchars = getarg(narg, arg)
!DEC$ ENDIF

      END SUBROUTINE getcarg
