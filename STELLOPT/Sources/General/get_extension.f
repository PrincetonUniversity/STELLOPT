      SUBROUTINE get_extension(in_file)
      USE optim_params, ONLY: seq_ext
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      CHARACTER(LEN=*) :: in_file
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      CHARACTER(LEN=6), PARAMETER :: input_ext = 'input.'
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: index_end, index_dat
C-----------------------------------------------
      index_dat = INDEX(in_file,input_ext)
      IF (index_dat .eq. 0) THEN
         in_file = input_ext // TRIM(in_file)
         index_dat = 1
      END IF
      index_dat = index_dat + LEN(input_ext)
      index_end = LEN_TRIM(in_file)
      seq_ext   = in_file(index_dat:index_end)

      END SUBROUTINE get_extension
