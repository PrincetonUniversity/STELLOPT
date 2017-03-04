!-----------------------------------------------------------------------
!     Subroutine:    write_pies_hdf5
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/3/2011
!     Description:   This subroutine outputs data to a HDF5 file.
!-----------------------------------------------------------------------
      SUBROUTINE write_pies_hdf5
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE pies_runtime
      USE pies_background
      USE pies_fieldlines
      USE pies_profile
      USE ez_hdf5
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER        :: ier
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      !-----  OPEN HDF5 FILE
         PRINT *,''
         PRINT *,'phiend',phiend
         PRINT *,'dphi',dphi
      CALL open_hdf5('test.h5',fid,ier)
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      !-----  WRITE SCALARS
      CALL write_scalar_hdf5(fid,'version',ier,FLTVAR=REAL(PIES_VERSION),ATT='Version Number')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_scalar_hdf5(fid,'iter',ier,INTVAR=iter,ATT='Iteration Number')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_scalar_hdf5(fid,'lasym',ier,BOOVAR=lasym,ATT='Non-Stellarator Symmetry Flag')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_scalar_hdf5(fid,'nextcur',ier,INTVAR=nextcur,ATT='Number of External Currents')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_scalar_hdf5(fid,'nfp',ier,INTVAR=nfp,ATT='Number of Field Periods')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_scalar_hdf5(fid,'m',ier,INTVAR=m,ATT='Poloidal modes')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_scalar_hdf5(fid,'n',ier,INTVAR=n,ATT='Toroidal modes')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_scalar_hdf5(fid,'k',ier,INTVAR=k,ATT='Radial surfaces')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_scalar_hdf5(fid,'nintw',ier,INTVAR=nintw,ATT='Number of field line punctures')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      
      CALL write_scalar_hdf5(fid,'phiend',ier,DBLVAR=phiend,ATT='Fieldline distance')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_scalar_hdf5(fid,'dphi',ier,FLTVAR=dphi,ATT='Fieldline distance interval')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      
      CALL write_scalar_hdf5(fid,'mnmax',ier,INTVAR=mnmax)
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_scalar_hdf5(fid,'curtor',ier,FLTVAR=REAL(curtor),ATT='Total Encolsed Toroidal Current')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_scalar_hdf5(fid,'phiedge',ier,FLTVAR=REAL(torflux_edge),ATT='Edge Toroidal Flux')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      !-----  WRITE VECTORS
      CALL write_vector_hdf5(fid,'xm',mnmax,ier,INTVAR=xm)
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_vector_hdf5(fid,'xn',mnmax,ier,INTVAR=xn)
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_vector_hdf5(fid,'rho',k+1,ier,DBLVAR=rho)
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_vector_hdf5(fid,'press',k+1,ier,DBLVAR=press)
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_vector_hdf5(fid,'torflux',k+1,ier,DBLVAR=torflux)
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_vector_hdf5(fid,'iprime',k+1,ier,DBLVAR=iprime)
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      !-----  WRITE ARRAYS
      CALL write_arr2d_hdf5(fid,'rmnc',mnmax,k+1,ier,DBLVAR=rmnc,ATT='R Fourier Harmonics (cos) [m]')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_arr2d_hdf5(fid,'zmns',mnmax,k+1,ier,DBLVAR=zmns,ATT='Z Fourier Harmonics (sin) [m]')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_arr2d_hdf5(fid,'bsmns',mnmax,k+1,ier,DBLVAR=bsmns,ATT='B^s Fourier Harmonics (sin) [T]')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_arr2d_hdf5(fid,'bumnc',mnmax,k+1,ier,DBLVAR=bumnc,ATT='B^u Fourier Harmonics (cos) [T]')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_arr2d_hdf5(fid,'bvmnc',mnmax,k+1,ier,DBLVAR=bvmnc,ATT='B^v Fourier Harmonics (cos) [T]')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_arr2d_hdf5(fid,'jsmns',mnmax,k+1,ier,DBLVAR=jsmns,ATT='J^s Fourier Harmonics (sin) [A/m^2]')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_arr2d_hdf5(fid,'jumnc',mnmax,k+1,ier,DBLVAR=jumnc,ATT='J^u Fourier Harmonics (cos) [A/m^2]')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      CALL write_arr2d_hdf5(fid,'jvmnc',mnmax,k+1,ier,DBLVAR=jvmnc,ATT='J^v Fourier Harmonics (cos) [A/m^2]')
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      IF (lasym) THEN
         CALL write_arr2d_hdf5(fid,'rmns',mnmax,k+1,ier,FLTVAR=REAL(rmns),ATT='R Fourier Harmonics (sin) [m]')
         IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
         CALL write_arr2d_hdf5(fid,'zmnc',mnmax,k+1,ier,FLTVAR=REAL(zmnc),ATT='Z Fourier Harmonics (cos) [m]')
         IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
         CALL write_arr2d_hdf5(fid,'bsmnc',mnmax,k+1,ier,FLTVAR=REAL(bsmnc),ATT='B^s Fourier Harmonics (cos) [T]')
         IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
         CALL write_arr2d_hdf5(fid,'bumns',mnmax,k+1,ier,FLTVAR=REAL(bumns),ATT='B^u Fourier Harmonics (sin) [T]')
         IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
         CALL write_arr2d_hdf5(fid,'bvmns',mnmax,k+1,ier,FLTVAR=REAL(bvmns),ATT='B^v Fourier Harmonics (sin) [T]')
         IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
         CALL write_arr2d_hdf5(fid,'jsmnc',mnmax,k+1,ier,FLTVAR=REAL(jsmnc),ATT='J^s Fourier Harmonics (cos) [A/m^2]')
         IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
         CALL write_arr2d_hdf5(fid,'jumns',mnmax,k+1,ier,FLTVAR=REAL(jumns),ATT='J^u Fourier Harmonics (sin) [A/m^2]')
         IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
         CALL write_arr2d_hdf5(fid,'jvmns',mnmax,k+1,ier,FLTVAR=REAL(jvmns),ATT='J^v Fourier Harmonics (sin) [A/m^2]')
         IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      END IF
      IF (ALLOCATED(xsiln)) THEN
         CALL write_arr2d_hdf5(fid,'xsiln',k+1,nintw+1,ier,DBLVAR=xsiln,ATT='Field Line location xsi [m]')
         IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
         CALL write_arr2d_hdf5(fid,'etaln',k+1,nintw+1,ier,DBLVAR=etaln,ATT='Field Line location eta [m]')
         IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
      END IF
      

      !-----  CLOSE HDF5 FILE
      CALL close_hdf5(fid,ier)
      IF (ier /= 0) CALL handle_err(HDF5_ERR,'write_pies_hdf5',ier)
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE write_pies_hdf5
