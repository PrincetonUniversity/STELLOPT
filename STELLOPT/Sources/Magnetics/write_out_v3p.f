      SUBROUTINE write_out_v3p
!**     SUBPROGRAM DESCRIPTION:                                      **
!**        write output netcdf file                                  **
!**                                                                  **
      USE stel_kinds
      USE stel_constants
      USE v3post_rfun
!DEC$ IF DEFINED (MPI_OPT)
      USE read_response
!DEC$ ELSE
      USE read_response_nompi
!DEC$ ENDIF
      USE ezcdf
      IMPLICIT NONE
!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
      INTEGER(iprec) :: istat, ncdf, index_nc, ndiag
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: signal_temp
      LOGICAL :: lnc
      CHARACTER (LEN=100) :: cdffil
      CHARACTER (LEN=*), PARAMETER ::
     1  vn_cal = 'signal_diag_cal',
     2  vn_cext = 'signal_diag_cext',
     3  vn_plasma = 'signal_diag_plasma',
     4  vn_sname = 'signal_diag_sname'
      CHARACTER (LEN=*), PARAMETER, DIMENSION(1) ::
     1             d1dim = (/'num_diagno'/)
      CHARACTER*(*), DIMENSION(2), PARAMETER ::
     1             d2dim = (/'str_len   ','num_diagno'/)
!----------------------------------------------------------------------
      cdffil='v3post' // '_' // TRIM(ideqfile)
      index_nc = INDEX(cdffil,".nc",BACK=.TRUE.)
      lnc = (index_nc .eq. (LEN_TRIM(cdffil)-2))
      IF(.not.lnc) cdffil = TRIM(cdffil) // '.nc'

      IF (.not. ALLOCATED(shortnames)) THEN             ! as in  vactest
        ALLOCATE(shortnames(SIZE(signal_diag%cal)))
        shortnames='UNDEFINED'
      ENDIF

      ndiag = SIZE(signal_diag%cal)
      
      ALLOCATE(signal_temp(ndiag))
      signal_temp = 0

      CALL cdf_open(ncdf, TRIM(cdffil), 'w', istat)
      CALL cdf_define(ncdf, vn_sname, shortnames, dimname=d2dim)
      CALL cdf_define(ncdf, vn_cal, signal_diag%cal, dimname=d1dim)
      CALL cdf_define(ncdf, vn_cext, signal_diag%cext, dimname=d1dim)
      CALL cdf_define(ncdf, vn_plasma, signal_temp, dimname=d1dim)
      CALL cdf_write(ncdf, vn_sname, shortnames)
      signal_temp = signal_diag%cal
      CALL cdf_write(ncdf, vn_cal, signal_temp)
      signal_temp = signal_diag%cext
      CALL cdf_write(ncdf, vn_cext, signal_temp)
      signal_temp = signal_diag%cal - signal_diag%cext
      CALL cdf_write(ncdf, vn_plasma, signal_temp)
      CALL cdf_close(ncdf)

      DEALLOCATE (signal_temp)
      IF (ALLOCATED(shortnames)) DEALLOCATE (shortnames)

      END SUBROUTINE write_out_v3p
