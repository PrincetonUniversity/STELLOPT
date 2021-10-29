
import os
import numpy as np
        
def read_finalFields(input_file, dim_kx, dim_ky):
    ''' Read the *.final_fields file: [z, z-thet0, ky, kx, real(phi), imag(phi), real(apar), imag(apar), z_eqarc-thet0]'''

    # Find the final fields file corresponding to the input file
    fields_file = input_file.with_suffix('.final_fields')
    if not os.path.isfile(fields_file):
        return False

    # Read the *.final_fields file: [z   z-thet0  ky  kx  real(phi)  imag(phi)  real(apar) imag(apar) z_eqarc-thet0]
    try:    fields_data  = np.loadtxt(fields_file, dtype='float').reshape(-1, dim_kx, dim_ky, 10) # Old code has 10 columns
    except: fields_data  = np.loadtxt(fields_file, dtype='float').reshape(-1, dim_kx, dim_ky, 11) # New code has 11 columns
    
    # Return the data
    return {"z" :         fields_data[:,:,:,0], \
            "zthet0" :    fields_data[:,:,:,1], \
            "ky" :        fields_data[:,:,:,2], \
            "kx" :        fields_data[:,:,:,3], \
            "phi_real" :  fields_data[:,:,:,4], \
            "phi_imag" :  fields_data[:,:,:,5], \
            "apar_real" : fields_data[:,:,:,6], \
            "apar_imag" : fields_data[:,:,:,7], \
            "z_eqarc" :   fields_data[:,:,:,8]} 


def get_finalFields(self):
    ''' Read the .final_fields file and attach phi and z to the simulation object. 
    
    Returns
    -------
    data : dict[z; kx; ky; phi_real; phi_imag; phi2; phi_end - phi_start; nfield_periods]
        Returns a dictionary with the most important information from the omega and netcdf files.
    
    data[z] : 1D array
    data[kx] : 1D array
    data[ky] : 1D array
    data[phi2] : 3D array with axis [z, kx, ky]
    data[phi_real] : 3D array with axis [z, kx, ky]
    data[phi_imag] : 3D array with axis [z, kx, ky]
    data[phi_end - phi_start] : 2D array with axis [kx, ky]
    nfield_periods : 2D array with axis [kx, ky]
    '''

    # Initiate the attributes wtih axis [z, kx, ky]
    self.phi_real   = np.empty((self.dim_z, self.dim_kx, self.dim_ky)); self.phi_real[:,:,:] = np.NaN
    self.phi_imag   = np.empty((self.dim_z, self.dim_kx, self.dim_ky)); self.phi_imag[:,:,:] = np.NaN
    self.phi2       = np.empty((self.dim_z, self.dim_kx, self.dim_ky)); self.phi2[:,:,:] = np.NaN
    self.z_poloidal = np.empty((self.dim_z, self.dim_kx, self.dim_ky)); self.phi2[:,:,:] = np.NaN
    self.zeta       = np.empty((self.dim_z, self.dim_kx, self.dim_ky)); self.phi2[:,:,:] = np.NaN
        
    # Read the length of one poloidal turn around the stellerator
    # 1 poloidal turn in w7x: nfieldperiods(s) = #poloidalturns * #fieldperiods / iota(s)
    z_one_turn = 1 * 5 / self.iota     
    
    # Get the data for the plots
    for input_file in self.input_files:
        
        # Read the number of field periods
        try:    nfield_periods = self.inputParameters['vmec_parameters']['nfield_periods']
        except: nfield_periods = self.inputs[input_file]['vmec_parameters']['nfield_periods']
        
        # Read zeta
        geometry_path = input_file.with_suffix(".geometry")
        geometry_data = np.loadtxt(geometry_path, dtype='float', skiprows=4).reshape(-1, 13)
        
        # Read the data from the .final_fields file
        fields_data = read_finalFields(input_file, self.dim_kxPerFile[input_file], self.dim_kyPerFile[input_file])
        
        # Save the zeta data per (kx,ky)
        for kx in self.vec_kxPerFile[input_file]:
            for ky in self.vec_kyPerFile[input_file]:
                
                # Get the indices of the mode and the dimensions of this file
                i_kx     = self.vec_kxPerFile[input_file].index(kx)
                i_ky     = self.vec_kyPerFile[input_file].index(ky)
                index_kx = self.vec_kx.index(kx)
                index_ky = self.vec_ky.index(ky) 
                dim_z    = self.dim_zPerFile[input_file]
                
                # Save the zeta data per (kx,ky)
                self.zeta[0:dim_z, index_kx,index_ky] = geometry_data[:,2]
            
                # Z in the netcdf file ranges from [-pi, pi], put this z in the number of poloidal turns
                z_normalized    = (self.z_kxky[:, index_kx, index_ky]+np.pi)/(2*np.pi)
                self.z_poloidal = z_normalized*nfield_periods/z_one_turn    

                # Save the potential
                if fields_data:
                    phi_real = fields_data["phi_real"]  
                    phi_imag = fields_data["phi_imag"]  
                    phi_complex = phi_real[:, i_kx, i_ky] + 1j*phi_imag[:, i_kx, i_ky]
                    self.phi_real[0:dim_z, index_kx, index_ky] = phi_real[:, i_kx, i_ky]
                    self.phi_imag[0:dim_z, index_kx, index_ky] = phi_imag[:, i_kx, i_ky]
                    self.phi2[0:dim_z, index_kx, index_ky] = abs(phi_complex*phi_complex.conjugate())
    return 
