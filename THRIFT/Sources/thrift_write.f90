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
         CALL write_scalar_hdf5(fid,'nssize',ier,INTVAR=nsj,ATT='Number of Radial Gridpoints in s space',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'nssize',ier)
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
         CALL write_var_hdf5(fid,'THRIFT_S',nsj,ier,DBLVAR=THRIFT_S,ATT='s Grid (r/a)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_S',ier)
         CALL write_var_hdf5(fid,'THRIFT_SNOB',nsj-2,ier,DBLVAR=THRIFT_SNOB,ATT='s grid (no boundaries) (r/a)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_SNOB',ier)
         CALL write_var_hdf5(fid,'THRIFT_RHOFULL',nrho+2,ier,DBLVAR=THRIFT_RHOFULL,ATT='Extended Radial Grid (r/a)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_RHOFULL',ier)
         CALL write_var_hdf5(fid,'THRIFT_T',ntimesteps,ier,DBLVAR=THRIFT_T,ATT='Time Slice Grid [s]',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_T',ier)
         CALL write_var_hdf5(fid,'THRIFT_PHIEDGE',ntimesteps,ier,DBLVAR=THRIFT_PHIEDGE,ATT='Toroidal magnetic flux at plasma edge [Wb] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_PHIEDGE',ier)
         ! 2D Floats
         ! Current densities
         nfg = nrho+2
         CALL write_var_hdf5(fid,'THRIFT_J',      nsj,ntimesteps,ier,DBLVAR=THRIFT_J,ATT='Total Current Density [A/m^2] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_J',ier)
         CALL write_var_hdf5(fid,'THRIFT_JBOOT',  nsj,ntimesteps,ier,DBLVAR=THRIFT_JBOOT,ATT='Bootstrap Current Density [A/m^2] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_JBOOT',ier)
         CALL write_var_hdf5(fid,'THRIFT_JECCD',  nsj,ntimesteps,ier,DBLVAR=THRIFT_JECCD,ATT='ECCD Current Density [A/m^2] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_JECCD',ier)
         CALL write_var_hdf5(fid,'THRIFT_JNBCD',  nsj,ntimesteps,ier,DBLVAR=THRIFT_JNBCD,ATT='NBCD Current Density [A/m^2] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_JNBCD',ier)
         CALL write_var_hdf5(fid,'THRIFT_JOHMIC', nsj,ntimesteps,ier,DBLVAR=THRIFT_JOHMIC,ATT='Ohmic Current Density [A/m^2] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_JOHMIC',ier)
         CALL write_var_hdf5(fid,'THRIFT_JPLASMA',nsj,ntimesteps,ier,DBLVAR=THRIFT_JPLASMA,ATT='Plasma Response Current Density [A/m^2] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_JPLASMA',ier)
         CALL write_var_hdf5(fid,'THRIFT_JSOURCE',nsj,ntimesteps,ier,DBLVAR=THRIFT_JSOURCE,ATT='Source Current Density [A/m^2] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_JSOURCE',ier)
         ! Currents
         CALL write_var_hdf5(fid,'THRIFT_UGRID',nsj,ntimesteps,ier,DBLVAR=THRIFT_UGRID,ATT='U (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_UGRID',ier)
         CALL write_var_hdf5(fid,'THRIFT_I',nsj,ntimesteps,ier,DBLVAR=THRIFT_I,ATT='Enclosed Total Current [A] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_I',ier)
         CALL write_var_hdf5(fid,'THRIFT_IBOOT',nsj,ntimesteps,ier,DBLVAR=THRIFT_IBOOT,ATT='Enclosed Bootstrap Current [A] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_IBOOT',ier)
         CALL write_var_hdf5(fid,'THRIFT_IECCD',nsj,ntimesteps,ier,DBLVAR=THRIFT_IECCD,ATT='Enclosed ECCD Current [A] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_IECCD',ier)
         CALL write_var_hdf5(fid,'THRIFT_INBCD',nsj,ntimesteps,ier,DBLVAR=THRIFT_INBCD,ATT='Enclosed NBCD Current [A] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_INBCD',ier)
         CALL write_var_hdf5(fid,'THRIFT_IOHMIC',nsj,ntimesteps,ier,DBLVAR=THRIFT_IOHMIC,ATT='Enclosed Ohmic Current [A] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_IOHMIC',ier)
         CALL write_var_hdf5(fid,'THRIFT_IPLASMA',nsj,ntimesteps,ier,DBLVAR=THRIFT_IPLASMA,ATT='Enclosed Plasma Response Current [A] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_IPLASMA',ier)
         CALL write_var_hdf5(fid,'THRIFT_ISOURCE',nsj,ntimesteps,ier,DBLVAR=THRIFT_ISOURCE,ATT='Enclosed Source Current [A] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_ISOURCE',ier)
         ! Profile variables
         CALL write_var_hdf5(fid,'THRIFT_ETAPARA',nsj,ntimesteps,ier,DBLVAR=THRIFT_ETAPARA,ATT=' Parallel electric resistivity (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_ETAPARA',ier)
         CALL write_var_hdf5(fid,'THRIFT_P',nsj,ntimesteps,ier,DBLVAR=THRIFT_P,ATT=' Pressure (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_P',ier)
         CALL write_var_hdf5(fid,'THRIFT_PPRIME',nsj,ntimesteps,ier,DBLVAR=THRIFT_PPRIME,ATT=' Radial derivative of pressure (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_PPRIME',ier)
         ! Magnetic variables
         CALL write_var_hdf5(fid,'THRIFT_VP',nsj,ntimesteps,ier,DBLVAR=THRIFT_VP,ATT='dV/dPhi (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_VP',ier)
         CALL write_var_hdf5(fid,'THRIFT_S11',nsj,ntimesteps,ier,DBLVAR=THRIFT_S11,ATT='Susceptance matrix element S11 (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_S11',ier)
         CALL write_var_hdf5(fid,'THRIFT_S12',nsj,ntimesteps,ier,DBLVAR=THRIFT_S12,ATT='Susceptance matrix element S12 (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_S12',ier)
         CALL write_var_hdf5(fid,'THRIFT_BAV',nsj,ntimesteps,ier,DBLVAR=THRIFT_BAV,ATT='<B> [T] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_BAV',ier)
         CALL write_var_hdf5(fid,'THRIFT_BVAV',nsj,ntimesteps,ier,DBLVAR=THRIFT_BVAV,ATT='<Bv> [T] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_BVAV',ier)
         CALL write_var_hdf5(fid,'THRIFT_BSQAV',nsj,ntimesteps,ier,DBLVAR=THRIFT_BSQAV,ATT='<B^2> [T^2] (s-space) ',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_BSQAV',ier)
         CALL write_var_hdf5(fid,'THRIFT_IOTA',nsj,ntimesteps,ier,DBLVAR=THRIFT_IOTA,ATT='Rotational transform (s-space) ',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_IOTA',ier)
         CALL write_var_hdf5(fid,'THRIFT_RMAJOR',nsj,ntimesteps,ier,DBLVAR=THRIFT_RMAJOR,ATT='Effective major radius [m] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_RMAJOR',ier)
         CALL write_var_hdf5(fid,'THRIFT_AMINOR',nsj,ntimesteps,ier,DBLVAR=THRIFT_AMINOR,ATT='Effective minor radius  [m] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_AMINOR',ier)
         ! Electric Field
         CALL write_var_hdf5(fid,'THRIFT_EPARB',nsj,ntimesteps,ier,DBLVAR=THRIFT_EPARB,ATT='<E.B> [V.T/m] (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_EPARB',ier)
         ! ABCD
         CALL write_var_hdf5(fid,'THRIFT_COEFF_A',nsj,ntimesteps,ier,DBLVAR=THRIFT_COEFF_A,ATT='Coefficient A (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_COEFF_A',ier)
         CALL write_var_hdf5(fid,'THRIFT_COEFF_B',nsj,ntimesteps,ier,DBLVAR=THRIFT_COEFF_B,ATT='Coefficient B (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_COEFF_B',ier)
         CALL write_var_hdf5(fid,'THRIFT_COEFF_C',nsj,ntimesteps,ier,DBLVAR=THRIFT_COEFF_C,ATT='Coefficient C (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_COEFF_C',ier)
         CALL write_var_hdf5(fid,'THRIFT_COEFF_D',nsj,ntimesteps,ier,DBLVAR=THRIFT_COEFF_D,ATT='Coefficient D (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_COEFF_D',ier)
         CALL write_var_hdf5(fid,'THRIFT_COEFF_BP',nsj,ntimesteps,ier,DBLVAR=THRIFT_COEFF_BP,ATT='Derivative of coefficient B (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_COEFF_BP',ier)
         CALL write_var_hdf5(fid,'THRIFT_COEFF_CP',nsj,ntimesteps,ier,DBLVAR=THRIFT_COEFF_CP,ATT='Derivative of coefficient C (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_COEFF_CP',ier)
         CALL write_var_hdf5(fid,'THRIFT_COEFF_DP',nsj,ntimesteps,ier,DBLVAR=THRIFT_COEFF_DP,ATT='Derivative of coefficient D (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_COEFF_DP',ier)
         ! Alphas
         CALL write_var_hdf5(fid,'THRIFT_ALPHA1',nsj-2,ntimesteps,ier,DBLVAR=THRIFT_ALPHA1,ATT='Coefficient alpha_1 (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_ALPHA1',ier)
         CALL write_var_hdf5(fid,'THRIFT_ALPHA2',nsj-2,ntimesteps,ier,DBLVAR=THRIFT_ALPHA2,ATT='Coefficient alpha_2 (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_ALPHA2',ier)
         CALL write_var_hdf5(fid,'THRIFT_ALPHA3',nsj-2,ntimesteps,ier,DBLVAR=THRIFT_ALPHA3,ATT='Coefficient alpha_3 (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_ALPHA3',ier)
         CALL write_var_hdf5(fid,'THRIFT_ALPHA4',nsj-2,ntimesteps,ier,DBLVAR=THRIFT_ALPHA4,ATT='Coefficient alpha_4 (s-space)',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_ALPHA4',ier)
         ! System of equations
         CALL write_var_hdf5(fid,'THRIFT_MATLD',nsj-1,ntimesteps,ier,DBLVAR=THRIFT_MATLD,ATT='Matrix equation lower diagonal',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_MATLD',ier)
         CALL write_var_hdf5(fid,'THRIFT_MATMD',nsj,ntimesteps,ier,DBLVAR=THRIFT_MATMD,ATT='Matrix equation main diagonal',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_MATMD',ier)
         CALL write_var_hdf5(fid,'THRIFT_MATUD',nsj-1,ntimesteps,ier,DBLVAR=THRIFT_MATUD,ATT='Matrix equation upper diagonal',ATT_NAME='description')
         IF (ier /= 0) CALL handle_err(HDF5_WRITE_ERR,'THRIFT_MATUD',ier)
         CALL write_var_hdf5(fid,'THRIFT_MATRHS',nsj,ntimesteps,ier,DBLVAR=THRIFT_MATRHS,ATT='Matrix equation RHS',ATT_NAME='description')
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

