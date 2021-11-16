   
#################################################################
#                        READ FLUXES
#################################################################

import numpy as np
from stellapy.utils.decorators import verbose_wrapper
 
#===========================================
#  READ THE FLUXES FILE
#===========================================

@verbose_wrapper
def read_fluxes(input_files, dim_species, sign_B):
    ''' Read the heat, momentum and particle fluxes calculated by stella.

    Return fluxes[time, q_flux_0, p_flux_0, v_flux_0, ...] 
    The data is read from the ".fluxes" files in <folder>.
    The "wout*" file is required to give the fluxes the correct sign based on sign(B).
    '''

    # Store the fluxes in a dictionary
    fluxes = {'restart_times' : []} 
    
    # Iterate over the folders
    for input_file in input_files:
        
        # Get the fluxes file
        flux_file = input_file.with_suffix(".fluxes")

        # Read the *.fluxes file: [time  pflx_1 pflx_2  vflx_1 vflx_2  qflx_1 qflx_2]
        flux_data = np.loadtxt(flux_file, dtype='float').reshape(-1, 1+3*dim_species) 

        # Store the data {time, q_flux, p_flux, v_flux} in a dictionary
        if input_files.index(input_file)==0:
            fluxes['time'] = flux_data[1:-1,0] 
            fluxes['restart_times'].append(fluxes['time'][0])
            for specie in range(dim_species):
                fluxes['p_flux_'+str(specie)] = abs(flux_data[1:-1,0*dim_species+specie+1]*sign_B)
                fluxes['v_flux_'+str(specie)] = abs(flux_data[1:-1,1*dim_species+specie+1]*sign_B)
                fluxes['q_flux_'+str(specie)] = abs(flux_data[1:-1,2*dim_species+specie+1]*sign_B)

        # If there are more runs in <runX> subfolders, add them to the exising fluxes dictionary
        if input_files.index(input_file)!=0:
            vec_time = flux_data[1:-1,0]
            fluxes['restart_times'].append(vec_time[0])
            fluxes['time'] = np.concatenate((fluxes['time'], vec_time))
            for specie in range(dim_species):
                vec_part = abs(flux_data[1:-1,0*dim_species+specie+1]*sign_B)   
                vec_momt = abs(flux_data[1:-1,1*dim_species+specie+1]*sign_B)   
                vec_heat = abs(flux_data[1:-1,2*dim_species+specie+1]*sign_B)   
                fluxes['p_flux_'+str(specie)] = np.concatenate((fluxes['p_flux_'+str(specie)], vec_part))
                fluxes['v_flux_'+str(specie)] = np.concatenate((fluxes['v_flux_'+str(specie)], vec_momt))
                fluxes['q_flux_'+str(specie)] = np.concatenate((fluxes['q_flux_'+str(specie)], vec_heat)) 

    # Put the time axis in an ascending order
    sorted_indexes = np.argsort(fluxes['time'])
    fluxes['time'] = fluxes['time'][sorted_indexes]
    for specie in range(dim_species):
        fluxes['p_flux_'+str(specie)] = fluxes['p_flux_'+str(specie)][sorted_indexes]
        fluxes['v_flux_'+str(specie)] = fluxes['v_flux_'+str(specie)][sorted_indexes]
        fluxes['q_flux_'+str(specie)] = fluxes['q_flux_'+str(specie)][sorted_indexes]
        
    # Remove the smallest restart time since it is the actual starting time
    fluxes['restart_times'] = sorted(fluxes['restart_times'])[1:]
        
    # Return fluxes[time, q_flux_0, p_flux_0, v_flux_0, q_flux_1, p_flux_1, v_flux_1, ...]
    return fluxes


#==================================================
#  ATTACH THE FLUXES DATA TO THE SIMULATION OBJECT
#==================================================

def get_nonlinearFluxes(self):
     
    # Read the fluxes data
    fluxes = read_fluxes(self.input_files, self.dim_species, self.sign_B)
    
    # Get the dimension of time, since these can be multiple merged netcdf files
    self.dim_timeFluxes = len(fluxes['p_flux_0'])
    
    # Save the restart times to investigate bugs at the restart
    self.restart_times = fluxes['restart_times']
    
    # Initiate the attributes
    self.vec_time = np.empty((self.dim_timeFluxes)); self.vec_time[:] = np.NaN
    self.p_flux = np.empty((self.dim_timeFluxes, self.dim_species)); self.p_flux_kxky[:,:] = np.NaN
    self.q_flux = np.empty((self.dim_timeFluxes, self.dim_species)); self.q_flux_kxky[:,:] = np.NaN
    self.v_flux = np.empty((self.dim_timeFluxes, self.dim_species)); self.v_flux_kxky[:,:] = np.NaN

    # Save the time axis
    self.vec_time = fluxes['time']

    # Save the fluxes data as attributes
    for i in range(self.dim_species):
        self.p_flux[:,i] = fluxes['p_flux_'+str(i)]
        self.q_flux[:,i] = fluxes['q_flux_'+str(i)]
        self.v_flux[:,i] = fluxes['v_flux_'+str(i)]
        
    # Clear the memory
    del fluxes
    return 
    
#-----------------------------
def get_linearFluxes(self): 
    ''' Get the fluxes from the fluxes file for each mode, assuming no restart files. '''
    
    # Initiate the attributes
    self.p_flux_kxky = np.empty((self.dim_time, self.dim_species, self.dim_kx, self.dim_ky)); self.p_flux_kxky[:,:,:,:] = np.NaN
    self.q_flux_kxky = np.empty((self.dim_time, self.dim_species, self.dim_kx, self.dim_ky)); self.q_flux_kxky[:,:,:,:] = np.NaN
    self.v_flux_kxky = np.empty((self.dim_time, self.dim_species, self.dim_kx, self.dim_ky)); self.v_flux_kxky[:,:,:,:] = np.NaN
     
    # The .fluxes file are only printed for ONE mode, therefore
    # they only make sense when we simulate only ONE mode per simulation
    for input_file in self.input_files:
        if (self.dim_kxPerFile[input_file] > 1) or (self.dim_kyPerFile[input_file] > 1):
            return
    
    # If we have one mode per simulation we can read the fluxes
    for input_file in self.input_files:
        
        # Show the reading progress
        if self.Progress: i = self.input_files.index(input_file); length = len(self.input_files)
        if self.Progress: self.Progress.move(i/length*100,"Reading the fluxes files ("+str(i)+"/"+str(length)+")")
        
        # Get the information about the mode
        kx = self.vec_kxPerFile[input_file][0]
        ky = self.vec_kyPerFile[input_file][0]
        index_kx = self.vec_kx.index(kx)
        index_ky = self.vec_ky.index(ky) 
         
        # Get the data required to read the fluxes
        dim_species = self.dim_species
        dim_time    = self.dim_timePerFile[input_file]
         
        # Read the fluxes file
        flux_file   = input_file.with_suffix(".fluxes")
        flux_data   = np.loadtxt(flux_file, dtype='float')[1:,:] 
        flux_data   = flux_data[0:dim_time, :]
        
        # Seperate the fluxes by heat; momentum and particle flux per species
        # First column is the time, then: p_flux1, p_flux2, p_flux3, v_flux1, v_flux2, v_flux3, q_flux1, q_flux2, q_flux3
        for i in range(dim_species):
            self.q_flux_kxky[0:dim_time,i,index_kx,index_ky] = flux_data[:,0*dim_species+i+1]
            self.v_flux_kxky[0:dim_time,i,index_kx,index_ky] = flux_data[:,1*dim_species+i+1]
            self.p_flux_kxky[0:dim_time,i,index_kx,index_ky] = flux_data[:,2*dim_species+i+1]
         
    return 

#-------------------------
def get_scaledLinearFluxes(self):
    
    # Initiate the attributes
    self.p_fluxQL_kxky = np.empty((self.dim_time, self.dim_species, self.dim_kx, self.dim_ky)); self.p_fluxQL_kxky[:,:,:,:] = np.NaN
        
    # The .fluxes file are only printed for ONE mode, therefore
    # they only make sense when we simulate only ONE mode per simulation
    for input_file in self.input_files:
        if (self.dim_kxPerFile[input_file] > 1) or (self.dim_kyPerFile[input_file] > 1):
            return
    
    # Go through the input files
    for input_file in self.input_files:
        
        # Get the information about the mode
        kx = self.vec_kxPerFile[input_file][0]
        ky = self.vec_kyPerFile[input_file][0]
        index_kx = self.vec_kx.index(kx)
        index_ky = self.vec_ky.index(ky) 
         
        # Get the data required to calculate the fluxes
        dim_species = self.dim_species 
        dim_time    = self.dim_timePerFile[input_file]
        omega       = self.omega_kxky[:,index_kx,index_ky]      
        gamma       = self.gamma_kxky[:,index_kx,index_ky]       
        phi2        = self.phi2_kxky[:,index_kx,index_ky]    
        p_flux      = self.p_flux_kxky[:,index_kx,index_ky]  
    
        # Divide density by the density of the ions.
        density = self.density_kxky[:,index_kx,index_ky]/self.density_kxky[0,index_kx,index_ky]
    
        # Calculate the quasi linear impurity ion flux based on literature
        for i in range(dim_species):  
            # 'mik' applies formula (1) of Mikkelen et al. PoP 21 082302 (2014)
            if self.formula_quasiLinearFluxes == 'mik':
                self.p_fluxQL_kxky[0:dim_time,i,index_kx,index_ky] = p_flux*gamma/phi2/density[i]/ky**2.0
            # 'per' applies formula (13)-(14) of Helander and Zocco PPCF 60 084006 (2018)
            elif self.formula_quasiLinearFluxes == 'per':
                try:    self.p_fluxQL_kxky[0:dim_time,i,index_kx,index_ky] = p_flux*(gamma**2.0+omega**2.0)/phi2/density[i]/ky**2.0/gamma
                except: self.p_fluxQL_kxky[0:dim_time,i,index_kx,index_ky] = np.NaN
    return

#-----------------------------
def get_quasiLinearFluxes(self):
    
    # Initiate the attributes: matrixes with dimensions (kx,ky,species)
    dim_species        = self.dim_species
    self.p_flux_last   = np.empty((self.dim_kx, self.dim_ky, dim_species)); self.p_flux_last[:,:,:] = np.NaN
    self.q_flux_last   = np.empty((self.dim_kx, self.dim_ky, dim_species)); self.q_flux_last[:,:,:] = np.NaN
    self.v_flux_last   = np.empty((self.dim_kx, self.dim_ky, dim_species)); self.v_flux_last[:,:,:] = np.NaN
    self.p_fluxQL_last = np.empty((self.dim_kx, self.dim_ky, dim_species)); self.p_fluxQL_last[:,:,:] = np.NaN
    self.p_flux_avg    = np.empty((self.dim_kx, self.dim_ky, dim_species)); self.p_flux_avg[:,:,:] = np.NaN
    self.q_flux_avg    = np.empty((self.dim_kx, self.dim_ky, dim_species)); self.q_flux_avg[:,:,:] = np.NaN
    self.v_flux_avg    = np.empty((self.dim_kx, self.dim_ky, dim_species)); self.v_flux_avg[:,:,:] = np.NaN
    self.p_fluxQL_avg  = np.empty((self.dim_kx, self.dim_ky, dim_species)); self.p_fluxQL_avg[:,:,:] = np.NaN
        
    # Isolate all the relevant variables to calculate the linear data                              
    dim_kx      = self.dim_kx 
    vec_kx      = self.vec_ky 
    dim_ky      = self.dim_ky 
    vec_ky      = self.vec_ky 
    dim_species = self.dim_species 
        
    # Iterate over each (kx,ky) couple and grab omega(t_last); gamma(t_last) and fluxes(t_last)
    # Here t_last is the last time value of the simulation, assuming the simulation has converged at this point.
    for index_kx in range(dim_kx):
        for index_ky in range(dim_ky):
            
            # Get the modes
            kx = vec_kx[index_kx]
            ky = vec_ky[index_ky]
            
            # Get the fluxes
            p_flux   = self.p_flux_kxky[:,:,index_kx,index_ky]
            q_flux   = self.q_flux_kxky[:,:,index_kx,index_ky]
            v_flux   = self.v_flux_kxky[:,:,index_kx,index_ky]
            p_fluxQL = self.p_fluxQL_kxky [:,:,index_kx,index_ky]
            
            # If no wave is examined then omega=0 and gamma=0 and we divide by 0
            if not (np.isclose(kx,0) and np.isclose(ky,0)): 
                
                # Initiate an array for the fluxes of the different species
                p_flux_last = np.empty((dim_species));  self.p_flux_last[:] = np.NaN
                q_flux_last = np.empty((dim_species));  self.q_flux_last[:] = np.NaN   
                v_flux_last = np.empty((dim_species));  self.v_flux_last[:] = np.NaN  
                p_flux_avg  = np.empty((dim_species));  self.p_flux_avg[:] = np.NaN
                q_flux_avg  = np.empty((dim_species));  self.q_flux_avg[:] = np.NaN
                v_flux_avg = np.empty((dim_species));   self.v_flux_avg[:] = np.NaN
                p_fluxQL_last = np.empty((dim_species)); self.p_fluxQL_last[:] = np.NaN
                p_fluxQL_avg = np.empty((dim_species));  self.p_fluxQL_avg[:] = np.NaN
                
                # Check whether all values are NaN
                if not np.all(~np.isfinite(p_flux)):
                
                    # Remove the infinite and NaN values  
                    notnan  = np.isfinite(p_flux[:,0])
                    for i in range(dim_species): 
                        notnan  = notnan & np.isfinite(p_flux[:,i])
                        notnan  = notnan & np.isfinite(q_flux[:,i])
                        notnan  = notnan & np.isfinite(v_flux[:,i])
    
                    # Iterate over the species
                    for i in range(dim_species):
                        # Grab omega, gamma and quasi-Linear fluxes at the last time value where the fluxes are finite. 
                        p_flux_last[i]   = p_flux[:,i][notnan][-1]
                        q_flux_last[i]   = q_flux[:,i][notnan][-1]
                        v_flux_last[i]   = v_flux[:,i][notnan][-1]
                        p_fluxQL_last[i] = p_fluxQL[:,i][notnan][-1]
                        # Average over the last 10% since before it converges it osicillates around its convergent point.
                        p_flux_avg[i]    = p_flux[:,i][notnan][int(np.size(p_flux)*90/100):]
                        q_flux_avg[i]    = q_flux[:,i][notnan][int(np.size(q_flux)*90/100):]
                        v_flux_avg[i]    = v_flux[:,i][notnan][int(np.size(v_flux)*90/100):]
                        p_fluxQL_avg[i]  = p_fluxQL[:,i][notnan][int(np.size(p_fluxQL)*90/100):]

            # And the fluxes in the (kx,ky,species) matrixes
            for i in range(dim_species):
                self.p_flux_last[index_kx,index_ky,i]   = p_flux_last[i]
                self.q_flux_last[index_kx,index_ky,i]   = q_flux_last[i]
                self.v_flux_last[index_kx,index_ky,i]   = v_flux_last[i]
                self.p_fluxQL_last[index_kx,index_ky,i] = p_fluxQL_last[i]
                self.p_flux_avg[index_kx,index_ky,i]    = p_flux_avg[i]
                self.q_flux_avg[index_kx,index_ky,i]    = q_flux_avg[i]
                self.v_flux_avg[index_kx,index_ky,i]    = v_flux_avg[i]
                self.p_fluxQL_avg[index_kx,index_ky,i]  = p_fluxQL_avg[i]
    return

