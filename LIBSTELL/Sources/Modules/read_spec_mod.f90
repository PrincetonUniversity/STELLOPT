!-----------------------------------------------------------------------
!     Module:        read_spec_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/26/2012
!     Description:   This subroutine reads the SPEC HDF5 file.  In order
!                    to ease integration into existing routines the
!                    module initializes the varialbes in the
!                    read_wout_mod VMEC output interface module.
!-----------------------------------------------------------------------
      MODULE read_spec_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE read_wout_mod, ONLY: nfp, ns, mpol,ntor, mnmax, rmax_surf, &
                               rmin_surf, zmax_surf, Rmajor, Itor, &
                               betatot, &
                               rmnc, rmns, zmnc, zmns, lmnc, lmns, &
                               bsupumnc, bsupumns, bsupvmnc, bsupvmns, &
                               iotaf, presf, phi, xm, xn, lasym, &
                               lthreed, input_extension, version_
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PRIVATE :: ier, i, j, k
!-----------------------------------------------------------------------
!     Subroutines
!          read_spec_file     Read the SPEC output file.
!-----------------------------------------------------------------------
      
      INTERFACE read_spec_file
          MODULE PROCEDURE read_spec_hdf5
      END INTERFACE
      
      CONTAINS
      
      SUBROUTINE read_spec_hdf5(id_string)
#ifdef LHDF5
      USE ez_hdf5
#endif
      CHARACTER(LEN=*), INTENT(in) :: id_string
      INTEGER :: ier, dex
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Check for extension or full name
      input_extension = ''
#ifdef LHDF5
      dex = INDEX(id_string,'h5',BACK=.TRUE.)
      IF (dex < LEN(id_string) .and. dex /= 0) THEN
         input_extension(1:dex-2) = id_string(1:dex-2)
      ELSE 
         input_extension = TRIM(id_string)
      END IF
      ! Default some stuff
      lasym = .FALSE.
      ! Open SPEC file
      CALL open_hdf5(TRIM(input_extension)//'.h5',fid,ier)
      CALL read_scalar_hdf5(fid,'Nvol',ier,INTVAR=ns)
      CALL read_scalar_hdf5(fid,'mn',ier,INTVAR=mnmax)
      CALL read_scalar_hdf5(fid,'version',ier,DBLVAR=version_)
      ALLOCATE(xm(mnmax),xn(mnmax),phi(ns+1),iotaf(ns+1))
      ALLOCATE(rmnc(mnmax,ns+1),zmns(mnmax,ns+1),bsupumnc(mnmax,1),bsupvmnc(mnmax,1))
      CALL read_var_hdf5(fid,'im',mnmax,ier,DBLVAR=xm)
      CALL read_var_hdf5(fid,'in',mnmax,ier,DBLVAR=xn)
      CALL read_var_hdf5(fid,'tflux',ns+1,ier,DBLVAR=phi)
      CALL read_var_hdf5(fid,'iota',ns+1,ier,DBLVAR=iotaf)
      CALL read_var_hdf5(fid,'Rbc',mnmax,ns,ier,DBLVAR=rmnc)
      CALL read_var_hdf5(fid,'Zbs',mnmax,ns,ier,DBLVAR=zmns)
      !CALL read_var_hdf5(fid,'Bsupumn',mnmax,1,ier,DBLVAR=bsupumnc)
      !CALL read_var_hdf5(fid,'Bsupvmn',mnmax,1,ier,DBLVAR=bsupvmnc)

      ! Handle some other values
      nfp = 0
      IF (ANY(xn > 0)) THEN
         nfp = MINVAL(xn, MASK = xn > 0)
      ELSE
         nfp = 1
      END IF
      mpol    = MAXVAL(xm)+1  ! VMEC convention
      ntor    = MAXVAL(xn) / nfp
      ! Close File
      CALL close_hdf5(fid,ier)
      RETURN
#else
      STOP 'No HDF5 Support, cannot read SPEC output'
#endif
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE read_spec_hdf5
!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE read_spec_mod
