      INTEGER FUNCTION vmec_chdir(new_path)
      IMPLICIT NONE
      CHARACTER*(*), INTENT(in) :: new_path
!DEC$ IF DEFINED (CRAY)
      INTEGER :: ilen
      iLEN = 0
      CALL pxfchdir(new_path, ilen, vmec_chdir)
!DEC$ ELSE
      INTEGER, EXTERNAL :: chdir
 
!DEC$ IF DEFINED (SUNOS)
      vmec_chdir = chdir(TRIM(new_path))
!DEC$ ELSE
!      vmec_chdir = chdir(TRIM(new_path) // char(0))
!DEC$ ENDIF
!DEC$ ENDIF
      END FUNCTION vmec_chdir
