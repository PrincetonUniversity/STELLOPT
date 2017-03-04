!-----------------------------------------------------------------------
!     Subroutine:    write_pies_netcdf
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/3/2011
!     Description:   This subroutine outputs data to a netCDF file.
!-----------------------------------------------------------------------
      SUBROUTINE write_pies_netcdf
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE pies_runtime
      USE pies_background
      USE ezcdf
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(LEN=*), PARAMETER :: &
         nc_version = 'version_', &
         nc_lwout = 'lwout',&
         nc_lasym = 'lasym',&
         nc_lfreeb = 'lfreeb',&
         nc_nextcur = 'nextcur',&
         nc_iter = 'iter'
      CHARACTER(LEN=*), PARAMETER :: &
         ln_version = 'iPIES Version Number',&
         ln_lwout = 'Start from wout file?',&
         ln_lasym = 'Non-stellarator symmetry?',&
         ln_lfreeb = 'Free Boundary Run?',&
         ln_nextcur = 'Number of coil current groups',&
         ln_iter = 'Iteration Number'
      
      
      
      
      INTEGER :: ier, ncid
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------  
      CALL cdf_open(ncid,'test.nc','w',ier)
      IF (ier /= 0) CALL handle_err(NETCDF_OPEN_ERR,'test.nc',ier)
      !-----  DEFINE VARIABLES
      CALL cdf_define(ncid,nc_version,PIES_VERSION)
      CALL cdf_define(ncid,nc_lwout,lwout)
      CALL cdf_define(ncid,nc_lasym,lasym)
      CALL cdf_define(ncid,nc_lfreeb,lfreeb)
      CALL cdf_define(ncid,nc_nextcur,nextcur)
      CALL cdf_define(ncid,nc_iter,iter)
      !-----  DEFINE ATTRIBUTES
      CALL cdf_setatt(ncid,nc_version,ln_version)
      CALL cdf_setatt(ncid,nc_lwout,ln_lwout)
      CALL cdf_setatt(ncid,nc_lasym,ln_lasym)
      CALL cdf_setatt(ncid,nc_lfreeb,ln_lfreeb)
      CALL cdf_setatt(ncid,nc_iter,ln_iter)
      !-----  WRITE VARIABLES
      CALL cdf_write(ncid,nc_version,PIES_VERSION)
      CALL cdf_write(ncid,nc_lwout,lwout)
      CALL cdf_write(ncid,nc_lasym,lasym)
      CALL cdf_write(ncid,nc_lfreeb,lfreeb)
      CALL cdf_write(ncid,nc_iter,iter)
      !-----  CLOSE NETCDF FILE
      CALL cdf_close(ncid)
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE write_pies_netcdf
