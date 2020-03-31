      INTEGER FUNCTION vmec_chdir(new_path)
      IMPLICIT NONE
      CHARACTER*(*), INTENT(in) :: new_path
#if defined(CRAY)
      INTEGER :: ilen
      iLEN = 0
      CALL pxfchdir(new_path, ilen, vmec_chdir)
#else
      INTEGER :: chdir
      vmec_chdir = chdir(TRIM(new_path))
#endif
      END FUNCTION vmec_chdir
