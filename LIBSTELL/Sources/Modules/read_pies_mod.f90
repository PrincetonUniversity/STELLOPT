!-----------------------------------------------------------------------
!     Module:        read_pies_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/26/2012
!     Description:   This subroutine reads the PIES netCDF file.  In order
!                    to ease integration into existing routines the
!                    module initializes the varialbes in the
!                    read_wout_mod VMEC output interface module.
!-----------------------------------------------------------------------
      MODULE read_pies_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE read_wout_mod, ONLY: nfp, ns, mpol,ntor, mnmax, rmax_surf, &
                               rmin_surf, zmax_surf, Rmajor, Itor, &
                               rmnc, rmns, zmnc, zmns, lmnc, lmns, &
                               bsupumnc, bsupumns, bsupvmnc, bsupvmns, &
                               iotaf, presf, phi, xm, xn, lasym, &
                               lthreed, input_extension, rprec
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, i, j, k
      REAL(rprec), ALLOCATABLE :: bsupsmns(:,:)
!-----------------------------------------------------------------------
!     Subroutines
!          read_spec_file     Read the SPEC output file.
!-----------------------------------------------------------------------
      
      INTERFACE read_pies_file
          MODULE PROCEDURE read_pies_hdf5
      END INTERFACE
      
      CONTAINS
      
      SUBROUTINE read_pies_hdf5(id_string)
!DEC$ IF DEFINED (NETCDF)
      include 'netcdf.inc'
!DEC$ ENDIF
      CHARACTER(LEN=*), INTENT(in) :: id_string
      INTEGER :: i,m,n,mn,ierr,i_alloc,ncid,status,& 
     &           BMID, BNID, BKID, BLID, &
     &           XMID, XNID, XKID, XLID, &
     &           XID, RELBUPID, NPERID,RMAJID,BXBYID,&
     &           BXMID, BXNID, BXKID, BXLID,&
     &           IOTAID,PSIID,&
     &           mpol,ntor,k,n2, dex
      REAL(rprec) :: temp
      REAL(rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: relbup,x,bxby
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Check for extension or full name
      input_extension = ''
!DEC$ IF DEFINED (NETCDF)
      dex = INDEX(id_string,'nc',BACK=.TRUE.)
      IF (dex < LEN(id_string)) THEN
         input_extension(1:dex-2) = id_string(1:dex-2)
      ELSE 
         input_extension = TRIM(id_string)
      END IF
      ! Default some stuff
      lasym = .FALSE.
      ! Open SPEC file
      status=nf_open(trim(id_string),0,ncid)
	   if (status .ne. NF_NOERR) then
	      write(*,*) '*****Error Opening netCDF file*****'
		   write(*,*) '     Filename=',trim(id_string)
		   write(*,*) '     Status=',NF_STRERROR(status)
		   stop
      endif
      status=nf_inq_dimid(ncid,'dim_x_1',XMID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_dimid(ncid,'dim_x_2',XNID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_dimid(ncid,'dim_x_3',XKID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_dimid(ncid,'dim_x_4',XLID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_dimid(ncid,'dim_relbup_old_1',BMID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_dimid(ncid,'dim_relbup_old_2',BNID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_dimid(ncid,'dim_relbup_old_3',BKID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_dimid(ncid,'dim_relbup_old_4',BLID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_dimid(ncid,'dim_bxby_1',BXMID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_dimid(ncid,'dim_bxby_2',BXNID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_dimid(ncid,'dim_bxby_3',BXKID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_dimid(ncid,'dim_bxby_4',BXLID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
!     Get the Variable ID's
      status=nf_inq_varid(ncid,'nper',NPERID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_varid(ncid,'rmaj',RMAJID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_varid(ncid,'iota',IOTAID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_varid(ncid,'psi',PSIID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_varid(ncid,'x',XID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_varid(ncid,'relbup_old',RELBUPID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_varid(ncid,'bxby',BXBYID)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
!     Get the Dimensions
      status=nf_inq_dimlen(ncid,XMID,mpol)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_dimlen(ncid,XNID,ntor)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_inq_dimlen(ncid,XKID,k)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
!     Allocate the bup and x arrays
      i_alloc=0
      ALLOCATE(iotaf(k),phi(k))
      ALLOCATE(relbup(mpol,ntor,k,3),x(mpol,ntor,k,2),&
     &         bxby(mpol,ntor,k,2))
      x=0; relbup=0;
!     Get the Variables
      status=nf_get_var_int(ncid,NPERID,nfp)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_get_var_double(ncid,RMAJID,Rmajor)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_get_var_double(ncid,IOTAID,iotaf)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_get_var_double(ncid,PSIID,phi)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_get_var_double(ncid,XID,x)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_get_var_double(ncid,RELBUPID,relbup)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
      status=nf_get_var_double(ncid,BXBYID,bxby)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
!     Close the netCDF File
	   status=nf_close(ncid)
	   if (status .ne. NF_NOERR) then
	      write(*,*) '*****Error Closing the netCDF file*****'
		   write(*,*) '     Filename=',trim(id_string)
		   write(*,*) '     Status=',NF_STRERROR(status)
		   stop
      endif
      ! Define some dimensional quantities
      mnmax = (2*ntor+1)*(mpol + 1)
      ns = k ! Note k is the dimension length not the varaible k
      ALLOCATE(xm(mnmax),xn(mnmax),phi(ns),iotaf(ns))
      ALLOCATE(rmnc(mnmax,ns),zmns(mnmax,ns),&
               bsupsmns(mnmax,ns),bsupumnc(mnmax,ns),bsupvmnc(mnmax,ns))
      mn = 1
      DO m = 0, mpol
         DO n = -ntor, ntor
            n2 = n+ntor+1
            xm(mn) = mpol
            xn(mn) = ntor*nfp
            rmnc(mn,:) =  x(m+1,n2,:,1)
            zmns(mn,:) = -x(m+1,n2,:,2)
            IF ((xm(mn) == 0) .and. (xn(mn) == 0)) rmnc(mn,:) = rmnc(mn,:) + Rmajor
            bsupsmns(mn,:) = -relbup(m+1,n2,:,1)
            bsupumnc(mn,:) =  relbup(m+1,n2,:,2)
            bsupvmnc(mn,:) =  relbup(m+1,n2,:,3)
            mn = mn + 1
         END DO
      END DO
      RETURN
!DEC$ ELSE
      STOP 'No netCDF Support, cannot read PIES output'
!DEC$ ENDIF
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE read_pies_hdf5

      SUBROUTINE HANDLE_ERR(STATUS)
!***********************************************************************
!  This subroutine handles netCDF errors.  Lifted straight from:
!        https://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f77/NF_005fSTRERROR.html#NF_005fSTRERROR
!
!  Written by: S. Lazerson (lazerson@pppl.gov)
!  Version:    2.0
!  Date:       03/26/12
!***********************************************************************
!DEC$ IF DEFINED (NETCDF)
      INCLUDE 'netcdf.inc'
!DEC$ ENDIF
      INTEGER STATUS
!DEC$ IF DEFINED (NETCDF)
      IF (STATUS .NE. NF_NOERR) THEN
        PRINT *, NF_STRERROR(STATUS)
        STOP '!!!!!!DIAGNO v2.0!!!!!'
      ENDIF
!DEC$ ENDIF
      END SUBROUTINE HANDLE_ERR
!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE read_pies_mod
