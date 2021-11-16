
import numpy as np
import os, configparser

#===========================================
#  READ THE OMEGA FILE
#===========================================
        
def read_omega(input_file, dim_time, dim_kx, dim_ky):
    ''' Read the *.omega file: [time,  ky,  kx,  Re(om),  Im(om),  Re(omavg), Im(omavg)]'''
    
    # Find the omega file corresponding to the input file
    omega_file = input_file.with_suffix(".omega")
    
    # Make sure we have the same time axis as in the netcdf file
    # Sometimes ".omega" is printed nefore ".out.nc" and it has more time values
    omega_data  = np.loadtxt(omega_file, dtype='float').reshape(-1, dim_kx, dim_ky, 7)
    omega_data  = omega_data[0:dim_time, :, :, :]
    
    # Return the data
    return {"omega" : omega_data[:,:,:,3], \
            "gamma" : omega_data[:,:,:,4]} 

#==================================================
#  ATTACH THE OMEGA DATA TO THE SIMULATION OBJECT
#==================================================

def get_omega(self):
    ''' Attach omega and gamma with axis [time, kx, ky] to the simulation object.'''
    
    # Initiate the attributes
    self.omega_kxky = np.empty((self.dim_time, self.dim_kx, self.dim_ky)); self.omega_kxky[:,:,:] = np.NaN
    self.gamma_kxky = np.empty((self.dim_time, self.dim_kx, self.dim_ky)); self.gamma_kxky[:,:,:] = np.NaN
    
    # Read the data for each input file
    for input_file in self.input_files:
        
        # Show the reading progress
        if self.Progress: i = self.input_files.index(input_file); length = len(self.input_files)
        if self.Progress: self.Progress.move(i/length*100,"Reading omega and gamma ("+str(i)+"/"+str(length)+")")
        
        # Read the omega file
        dim_time   = self.dim_timePerFile[input_file]
        dim_kx     = self.dim_kxPerFile[input_file]
        dim_ky     = self.dim_kyPerFile[input_file]
        omega_data = read_omega(input_file, dim_time, dim_kx, dim_ky)
        
        # Save the omega data per (kx,ky)
        for kx in self.vec_kxPerFile[input_file]:
            for ky in self.vec_kyPerFile[input_file]:
                i_kx     = self.vec_kxPerFile[input_file].index(kx)
                i_ky     = self.vec_kyPerFile[input_file].index(ky)
                index_kx = self.vec_kx.index(kx)
                index_ky = self.vec_ky.index(ky) 
                self.omega_kxky[0:dim_time, index_kx, index_ky] = omega_data['omega'][:, i_kx, i_ky]
                self.gamma_kxky[0:dim_time, index_kx, index_ky] = omega_data['gamma'][:, i_kx, i_ky]
                # Multiply omega with sgn(b0) to make sure it has the correct sign 
                # This can only be done correctly if the wout file is present in the folder! 
                self.omega_kxky[:, index_kx, index_ky] = self.omega_kxky[:, index_kx, index_ky]*self.sign_B
    return



#------------------------
def write_linearOmega(self):
    ''' Write the linear data to a configuration file to save calculation time. '''
    
    # The section headers to save the data
    omega_header = 'OMEGA AT LAST TIME VALUE; AVERAGE OF OMEGA OVER LAST 10%; ERROR ON AVERAGE'
    gamma_header = 'GAMMA AT LAST TIME VALUE; AVERAGE OF GAMMA OVER LAST 10%; ERROR ON AVERAGE'
    
    # If the section already existed, remove it
    if omega_header in self.simulation_file: self.simulation_file.remove_section(omega_header)
    if gamma_header in self.simulation_file: self.simulation_file.remove_section(gamma_header)
    
    # Create the empty sections
    self.simulation_file[omega_header] = {}
    self.simulation_file[gamma_header] = {}
    
    # Initiate the attributes: matrixes with dimensions (kx,ky)
    self.omega_last    = np.empty((self.dim_kx, self.dim_ky)); self.omega_last[:,:] = np.NaN
    self.gamma_last    = np.empty((self.dim_kx, self.dim_ky)); self.gamma_last[:,:] = np.NaN
    self.omega_avg     = np.empty((self.dim_kx, self.dim_ky)); self.omega_avg[:,:]  = np.NaN
    self.gamma_avg     = np.empty((self.dim_kx, self.dim_ky)); self.gamma_avg[:,:]  = np.NaN
    self.gamma_fromPhi = np.empty((self.dim_kx, self.dim_ky)); self.gamma_fromPhi[:,:] = np.NaN
    self.omega_min     = np.empty((self.dim_kx, self.dim_ky)); self.omega_avg[:,:]  = np.NaN
    self.omega_max     = np.empty((self.dim_kx, self.dim_ky)); self.omega_avg[:,:]  = np.NaN
    self.gamma_min     = np.empty((self.dim_kx, self.dim_ky)); self.gamma_avg[:,:]  = np.NaN
    self.gamma_max     = np.empty((self.dim_kx, self.dim_ky)); self.gamma_avg[:,:]  = np.NaN

    # Iterate over each (kx,ky) couple and grab omega(t_last); gamma(t_last) and fluxes(t_last)
    # Here t_last is the last time value of the simulation, assuming the simulation has converged at this point.
    for i in range(self.dim_kx):
        for j in range(self.dim_ky):
            
            # Show the reading progress
            if self.Progress: a = i*self.dim_ky + j; length = self.dim_kx*self.dim_ky
            if self.Progress: self.Progress.move(a/length*100,"Calculating linear omega and gamma ("+str(a)+"/"+str(length)+")")
                 
            # If no wave is examined then omega=0 and gamma=0 and we divide by 0
            if not (np.isclose(self.vec_kx[i],0) and np.isclose(self.vec_ky[j],0)): 
                
                # Isolate all the relevant variables to calculate the linear data                                  
                omega  = self.omega_kxky[:,i,j]
                gamma  = self.gamma_kxky[:,i,j]
                
                # Remove the infinite and NaN values 
                gamma  = gamma[np.isfinite(omega)]
                omega  = omega[np.isfinite(omega)]
    
                # Grab omega and gamma at the last time value where omega is finite. 
                omega_last = omega[-1]
                gamma_last = gamma[-1]
                
                # Average over the last 10% since before it converges it osicillates around its convergent point.
                omega_avg = np.average(omega[int(np.size(omega)*90/100):])
                gamma_avg = np.average(gamma[int(np.size(gamma)*90/100):])
                
                # Calculate the error bars over the last 10%:
                omega_min = np.abs(np.nanmin(omega[int(np.size(omega)*90/100):])-omega_avg)
                omega_max = np.abs(np.nanmax(omega[int(np.size(omega)*90/100):])-omega_avg)
                gamma_min = np.abs(np.nanmin(gamma[int(np.size(gamma)*90/100):])-gamma_avg)
                gamma_max = np.abs(np.nanmax(gamma[int(np.size(gamma)*90/100):])-gamma_avg)
    
                # Don't add the data from the removed modes
                if (self.vec_kx[i] not in self.removed_kx) and (self.vec_ky[j] not in self.removed_ky):
                    self.omega_last[i,j]    = omega_last
                    self.gamma_last[i,j]    = gamma_last
                    self.omega_avg[i,j]     = omega_avg
                    self.gamma_avg[i,j]     = gamma_avg
                    self.omega_min[i,j]     = omega_min
                    self.omega_max[i,j]     = omega_max
                    self.gamma_min[i,j]     = gamma_min
                    self.gamma_max[i,j]     = gamma_max
            
            # Save the linear data: [omega_last, omega_avg, omega_min, omega_max]
            kx = self.vec_kx[i];  ky = self.vec_ky[j]; 
            mode_indices = '(' + str(kx) + ", " + str(ky) + ')'
            omega_data = [ str(k) for k in [omega_last, omega_avg, omega_min, omega_max] ]
            gamma_data = [ str(k) for k in [gamma_last, gamma_avg, gamma_min, gamma_max] ]
            self.simulation_file[omega_header][mode_indices] = "[" + ", ".join(omega_data) + "]"
            self.simulation_file[gamma_header][mode_indices] = "[" + ", ".join(gamma_data) + "]"
    
    # Write the simulation file
    self.simulation_file.write(open(self.simulation_path, 'w'))
    return


#------------------------
def get_linearOmega(self):
    ''' Read the linear data from the simulation file. '''
    
    # If the linear data isn't written in the simulation file then write it.
    omega_header = 'OMEGA AT LAST TIME VALUE; AVERAGE OF OMEGA OVER LAST 10%; ERROR ON AVERAGE'
    gamma_header = 'GAMMA AT LAST TIME VALUE; AVERAGE OF GAMMA OVER LAST 10%; ERROR ON AVERAGE'
    if (omega_header not in self.simulation_file) or (gamma_header not in self.simulation_file):
        print("WRITE OMEGA AGAIN BECAUSE THE HEADER WAS MISSING")
        write_linearOmega(self)
        return
    
    # Initiate the attributes: matrixes with dimensions (kx,ky)
    self.omega_last    = np.empty((self.dim_kx, self.dim_ky)); self.omega_last[:,:] = np.NaN
    self.gamma_last    = np.empty((self.dim_kx, self.dim_ky)); self.gamma_last[:,:] = np.NaN
    self.omega_avg     = np.empty((self.dim_kx, self.dim_ky)); self.omega_avg[:,:]  = np.NaN
    self.gamma_avg     = np.empty((self.dim_kx, self.dim_ky)); self.gamma_avg[:,:]  = np.NaN
    self.omega_min     = np.empty((self.dim_kx, self.dim_ky)); self.omega_avg[:,:]  = np.NaN
    self.omega_max     = np.empty((self.dim_kx, self.dim_ky)); self.omega_avg[:,:]  = np.NaN
    self.gamma_min     = np.empty((self.dim_kx, self.dim_ky)); self.gamma_avg[:,:]  = np.NaN
    self.gamma_max     = np.empty((self.dim_kx, self.dim_ky)); self.gamma_avg[:,:]  = np.NaN

    # Iterate over each (kx,ky) couple and grab omega(t_last); gamma(t_last) and fluxes(t_last)
    # Here t_last is the last time value of the simulation, assuming the simulation has converged at this point.
    try:
        for i in range(self.dim_kx):
            for j in range(self.dim_ky):
                
                # Show the reading progress
                if self.Progress: a = i*self.dim_ky + j; length = self.dim_kx*self.dim_ky
                if self.Progress: self.Progress.move(a/length*100,"Reading linear omega and gamma ("+str(a)+"/"+str(length)+")")
                
                # Don't add the data from the removed modes
                if (self.vec_kx[i] not in self.removed_kx) and (self.vec_ky[j] not in self.removed_ky):
                        
                    # Get the mode indices which is the key in the simulation file
                    kx = self.vec_kx[i];  ky = self.vec_ky[j]; 
                    mode_indices = '(' + str(kx) + ", " + str(ky) + ')'
                         
                    # Read the linear data: [omega_last, omega_avg, omega_min, omega_max]
                    omega_data = self.simulation_file[omega_header][mode_indices].split("[")[-1].split("]")[0].split(", ") 
                    gamma_data = self.simulation_file[gamma_header][mode_indices].split("[")[-1].split("]")[0].split(", ") 
            
                    # Make sure we have floats
                    omega_data = [ float(x) for x in omega_data ]
                    gamma_data = [ float(x) for x in gamma_data ]
                    
                    # Now put this data in the (kx,ky) matrixes
                    self.omega_last[i,j] = omega_data[0]
                    self.omega_avg[i,j]  = omega_data[1]
                    self.omega_min[i,j]  = omega_data[2]
                    self.omega_max[i,j]  = omega_data[3]
                    self.gamma_last[i,j] = gamma_data[0]
                    self.gamma_avg[i,j]  = gamma_data[1]
                    self.gamma_min[i,j]  = gamma_data[2]
                    self.gamma_max[i,j]  = gamma_data[3]
        
    except:
        print("WRITE OMEGA AGAIN BECAUSE THE SOME MODES WERE MISSING")
        write_linearOmega(self)
    
    
                 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
