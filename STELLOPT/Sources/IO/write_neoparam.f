      SUBROUTINE write_neoparam (no_fluxs, fluxs_arr, extension,
     1    write_progress, istat)
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: no_fluxs, fluxs_arr(*), istat, write_progress
      CHARACTER(LEN=*) :: extension
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: max_m_mode = 0, max_n_mode = 0     !LET NEO COMPUTE INTERNALLY
      INTEGER, PARAMETER :: write_output_files = 0, WRITE_integrate = 0,
     1     write_diagnostic = 0, calc_cur = 0, WRITE_cur_inte = 0
      INTEGER :: w_neo = 10
      CHARACTER(LEN=200) :: arg1
C-----------------------------------------------

      arg1 = "neo_in." // TRIM(extension)

      CALL safe_open(w_neo, istat, arg1, 'replace', 'formatted')
      IF (istat .ne. 0) RETURN

      WRITE (w_neo,*) "#"
      WRITE (w_neo,*) "#"
      WRITE (w_neo,*) "#"
!DEC$ IF DEFINED(NETCDF)
      WRITE (w_neo,'(a)') "boozmn." // TRIM(extension)
!DEC$ ELSE
      WRITE (w_neo,'(a)') "boozmn_" // TRIM(extension) // ".nc"
!DEC$ ENDIF
      WRITE (w_neo,'(a)') "neo_out." // TRIM(extension)

      WRITE (w_neo,*) no_fluxs

      IF (no_fluxs .le. 0) THEN
         WRITE (w_neo,*) "#"
      ELSE
         WRITE (w_neo,*) fluxs_arr(1:no_fluxs)
      END IF

      WRITE (w_neo,*) 100
      WRITE (w_neo,*) 100
      WRITE (w_neo,*) max_m_mode
      WRITE (w_neo,*) max_n_mode
      WRITE (w_neo,*) 50
      WRITE (w_neo,*) 1
      WRITE (w_neo,*) 0.01
      WRITE (w_neo,*) 100
      WRITE (w_neo,*) 50
      WRITE (w_neo,*) 500
      WRITE (w_neo,*) 3000
      WRITE (w_neo,*) 0
      WRITE (w_neo,*) 1
      WRITE (w_neo,*) 0
      WRITE (w_neo,*) 0
      WRITE (w_neo,*) 2
      WRITE (w_neo,*) write_progress
      WRITE (w_neo,*) write_output_files
      WRITE (w_neo,*) 0
      WRITE (w_neo,*) write_integrate
      WRITE (w_neo,*) write_diagnostic
      WRITE (w_neo,*) "#"
      WRITE (w_neo,*) "#"
      WRITE (w_neo,*) "#"
      WRITE (w_neo,*) calc_cur
      WRITE (w_neo,*) "Bmns_cur.dat"
      WRITE (w_neo,*) 200
      WRITE (w_neo,*) 2
      WRITE (w_neo,*) write_cur_inte

      CLOSE (unit=w_neo)

      END SUBROUTINE write_neoparam
