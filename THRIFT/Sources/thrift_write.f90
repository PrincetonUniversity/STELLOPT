!-----------------------------------------------------------------------
!     Subroutine:    thrift_write
!     Authors:       S. Lazerson
!     Date:          12/20/2022
!     Description:   This subroutine writes the HDF5 output for THRIFT.
!-----------------------------------------------------------------------
      SUBROUTINE thrift_write
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE thrift_runtime
      USE thrift_vars
#if defined(LHDF5)
      USE ez_hdf5
#endif
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, nfg
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      IF (myworkid == master) THEN
#if defined(LHDF5)
         CALL open_hdf5('thrift_'//TRIM(id_string)//'.h5',fid,ier,LCREATE=.true.)
         IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,'thrift_'//TRIM(id_string)//'.h5',ier)
         ! Version
         CALL write_scalar_hdf5(fid,'VERSION',ier,DBLVAR=THRIFT_VERSION,ATT='Version Number',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'VERSION',ier)
         ! Logicals
         CALL write_scalar_hdf5(fid,'lvmec',ier,BOOVAR=lvmec,ATT='VMEC input',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lvmec',ier)
         CALL write_scalar_hdf5(fid,'leccd',ier,BOOVAR=leccd,ATT='ECCD Present',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'leccd',ier)
         CALL write_scalar_hdf5(fid,'lnbcd',ier,BOOVAR=lnbcd,ATT='NBCD Present',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lnbcd',ier)
         CALL write_scalar_hdf5(fid,'lohmic',ier,BOOVAR=lohmic,ATT='Ohmic CD Present',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'lohmic',ier)
         ! Integers
         CALL write_scalar_hdf5(fid,'ntimesteps',ier,INTVAR=ntimesteps,ATT='Number of Time Steps',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'ntimesteps',ier)
         CALL write_scalar_hdf5(fid,'nrho',ier,INTVAR=nrho,ATT='Number of Radial Gridpoints',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nrho',ier)
         CALL write_scalar_hdf5(fid,'npicard',ier,INTVAR=npicard,ATT='Maximum number of Picard Iterations',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'npicard',ier)
         ! Floats
         CALL write_scalar_hdf5(fid,'jtol',ier,DBLVAR=jtol,ATT='j convergence factor [%]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'jtol',ier)
         CALL write_scalar_hdf5(fid,'picard_factor',ier,DBLVAR=picard_factor,ATT='Picard Iteration Factor',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'picard_factor',ier)
         ! 1D Floats
         CALL write_var_hdf5(fid,'THRIFT_RHO',nrho,ier,DBLVAR=THRIFT_RHO,ATT='Radial Grid (r/a)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_RHO',ier)
         CALL write_var_hdf5(fid,'THRIFT_RHOFULL',nrho+2,ier,DBLVAR=THRIFT_RHOFULL,ATT='Extended Radial Grid (r/a)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_RHOFULL',ier)
         CALL write_var_hdf5(fid,'THRIFT_T',ntimesteps,ier,DBLVAR=THRIFT_T,ATT='Time Slice Grid [s]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_T',ier)
         ! 2D Floats
         ! Current densities
         nfg = nrho+2
         CALL write_var_hdf5(fid,'THRIFT_J',nrho,ntimesteps,ier,DBLVAR=THRIFT_J,ATT='Total Current Density [A/m^2]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_J',ier)
         CALL write_var_hdf5(fid,'THRIFT_JBOOT',nrho,ntimesteps,ier,DBLVAR=THRIFT_JBOOT,ATT='Total Bootstrap Current Density [A/m^2]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_JBOOT',ier)
         CALL write_var_hdf5(fid,'THRIFT_JECCD',nrho,ntimesteps,ier,DBLVAR=THRIFT_JECCD,ATT='Total ECCD Current Density [A/m^2]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_JECCD',ier)
         CALL write_var_hdf5(fid,'THRIFT_JNBCD',nrho,ntimesteps,ier,DBLVAR=THRIFT_JNBCD,ATT='Total NBCD Current Density [A/m^2]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_JNBCD',ier)
         CALL write_var_hdf5(fid,'THRIFT_JOHMIC',nrho,ntimesteps,ier,DBLVAR=THRIFT_JOHMIC,ATT='Total Ohmic Current Density [A/m^2]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_JOHMIC',ier)
         CALL write_var_hdf5(fid,'THRIFT_JPLASMA',nrho,ntimesteps,ier,DBLVAR=THRIFT_JPLASMA,ATT='Total Plasma Response Current Density [A/m^2]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_JPLASMA',ier)
         CALL write_var_hdf5(fid,'THRIFT_JSOURCE',nrho,ntimesteps,ier,DBLVAR=THRIFT_JSOURCE,ATT='Total Source Current Density [A/m^2]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_JSOURCE',ier)
         ! Currents
         CALL write_var_hdf5(fid,'THRIFT_UGRID',nfg,ntimesteps,ier,DBLVAR=THRIFT_UGRID,ATT='U on the full grid',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_UGRID',ier)
         CALL write_var_hdf5(fid,'THRIFT_I',nfg,ntimesteps,ier,DBLVAR=THRIFT_I,ATT='Total Enclosed Current [A]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_I',ier)
         CALL write_var_hdf5(fid,'THRIFT_IBOOT',nfg,ntimesteps,ier,DBLVAR=THRIFT_IBOOT,ATT='Total Enclosed Bootstrap Current [A]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_IBOOT',ier)
         CALL write_var_hdf5(fid,'THRIFT_IECCD',nfg,ntimesteps,ier,DBLVAR=THRIFT_IECCD,ATT='Total Enclosed ECCD Current [A]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_IECCD',ier)
         CALL write_var_hdf5(fid,'THRIFT_INBCD',nfg,ntimesteps,ier,DBLVAR=THRIFT_INBCD,ATT='Total Enclosed NBCD Current [A]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_INBCD',ier)
         CALL write_var_hdf5(fid,'THRIFT_IOHMIC',nfg,ntimesteps,ier,DBLVAR=THRIFT_IOHMIC,ATT='Total Enclosed Ohmic Current [A]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_IOHMIC',ier)
         CALL write_var_hdf5(fid,'THRIFT_IPLASMA',nfg,ntimesteps,ier,DBLVAR=THRIFT_IPLASMA,ATT='Total Enclosed Plasma Response Current [A]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_IPLASMA',ier)
         CALL write_var_hdf5(fid,'THRIFT_ISOURCE',nfg,ntimesteps,ier,DBLVAR=THRIFT_ISOURCE,ATT='Total Enclosed Source Current [A]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_ISOURCE',ier)
         ! Coefficients
         CALL write_var_hdf5(fid,'THRIFT_COEFF_A',nfg,ntimesteps,ier,DBLVAR=THRIFT_COEFF_A,ATT='Coefficient A',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_COEFF_A',ier)
         CALL write_var_hdf5(fid,'THRIFT_COEFF_B',nfg,ntimesteps,ier,DBLVAR=THRIFT_COEFF_B,ATT='Coefficient B',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_COEFF_B',ier)
         CALL write_var_hdf5(fid,'THRIFT_COEFF_C',nfg,ntimesteps,ier,DBLVAR=THRIFT_COEFF_C,ATT='Coefficient C',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_COEFF_C',ier)
         CALL write_var_hdf5(fid,'THRIFT_COEFF_D',nfg,ntimesteps,ier,DBLVAR=THRIFT_COEFF_D,ATT='Coefficient D',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_COEFF_D',ier)
         
         CALL write_var_hdf5(fid,'THRIFT_ALPHA1',nrho,ntimesteps,ier,DBLVAR=THRIFT_ALPHA1,ATT='Coefficient alpha_1',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_ALPHA1',ier)
         CALL write_var_hdf5(fid,'THRIFT_ALPHA2',nrho,ntimesteps,ier,DBLVAR=THRIFT_ALPHA2,ATT='Coefficient alpha_2',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_ALPHA2',ier)
         CALL write_var_hdf5(fid,'THRIFT_ALPHA3',nrho,ntimesteps,ier,DBLVAR=THRIFT_ALPHA3,ATT='Coefficient alpha_3',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_ALPHA3',ier)
         CALL write_var_hdf5(fid,'THRIFT_ALPHA4',nrho,ntimesteps,ier,DBLVAR=THRIFT_ALPHA4,ATT='Coefficient alpha_4',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_ALPHA4',ier)
         
         CALL write_var_hdf5(fid,'THRIFT_MATLD',nfg,ntimesteps,ier,DBLVAR=THRIFT_MATLD,ATT='Matrix equation lower diagonal',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_MATLD',ier)
         CALL write_var_hdf5(fid,'THRIFT_MATMD',nfg,ntimesteps,ier,DBLVAR=THRIFT_MATMD,ATT='Matrix equation main diagonal',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_MATMD',ier)
         CALL write_var_hdf5(fid,'THRIFT_MATUD',nfg,ntimesteps,ier,DBLVAR=THRIFT_MATUD,ATT='Matrix equation upper diagonal',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_MATUD',ier)
         CALL write_var_hdf5(fid,'THRIFT_MATRHS',nfg,ntimesteps,ier,DBLVAR=THRIFT_MATRHS,ATT='Matrix equation RHS',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_MATRHS',ier)
                  
         CALL close_hdf5(fid,ier)
         IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,'thrift_'//TRIM(id_string)//'.h5',ier)
#else
         WRITE(6,'(A)')  '   FILE: '//'thrift_'//TRIM(id_string)//'.bin'
         CALL safe_open(iunit,ier,'thrift_'//TRIM(id_string)//'.bin','replace','unformatted')
         ! Version
         WRITE(iunit,*) THRIFT_VERSION
         ! Logicals
         WRITE(iunit,*) lvmec, leccd, lnbcd, lohmic
         ! Integers
         WRITE(iunit,*) ntimesteps, nrho, npicard
         ! Floats
         WRITE(iunit,*) jtol, picard_factor
         ! 1D Floats
         WRITE(iunit,*) THRIFT_RHO
         WRITE(iunit,*) THRIFT_T
         ! 2D Floats
         WRITE(iunit,*) THRIFT_J
         WRITE(iunit,*) THRIFT_JBOOT
         WRITE(iunit,*) THRIFT_JECCD
         WRITE(iunit,*) THRIFT_JNBCD
         WRITE(iunit,*) THRIFT_JOHMIC
         WRITE(iunit,*) THRIFT_JPLASMA
         WRITE(iunit,*) THRIFT_JSOURCE

         WRITE(iunit,*) THRIFT_I
         WRITE(iunit,*) THRIFT_IBOOT
         WRITE(iunit,*) THRIFT_IECCD
         WRITE(iunit,*) THRIFT_INBCD
         WRITE(iunit,*) THRIFT_IOHMIC
         WRITE(iunit,*) THRIFT_IPLASMA
         WRITE(iunit,*) THRIFT_ISOURCE
         CLOSE(iunit)
#endif
      END IF

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE thrift_write

