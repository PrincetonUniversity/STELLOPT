
# External modules
import numpy as np
import os, configparser, pathlib, time

# Configuration and utilities
from stellapy.utils import  get_filesInFolder
from stellapy.utils.decorators import printv
from stellapy.config import turnOnVerboseWrapper_configurationFile, turnOffVerboseWrapper_configurationFile

# Functions to create the simulations
from stellapy.simulations.utils.load_labelsLinesMarkers import load_labelsLinesMarkers
from stellapy.simulations.utils.get_simulationIdentifier import get_simulationIdentifier
from stellapy.simulations.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime
from stellapy.simulations.utils.remove_simulationsWithoutOutncFile import remove_simulationsWithoutOutncFile

# Read data in the netcdf file
from stellapy.data.read_netcdf import get_netcdfDimensions
from stellapy.data.read_netcdf import get_netcdfTimeVsKxKy
from stellapy.data.read_netcdf import get_netcdfZVsKxKy
from stellapy.data.read_netcdf import get_netcdfPhi2VsKxKy
from stellapy.data.read_netcdf import get_netcdfFluxesKxKy
from stellapy.data.read_netcdf import get_netcdfDensity

# Read data from the fluxes file
from stellapy.data.read_fluxes import get_linearFluxes 
from stellapy.data.read_fluxes import get_quasiLinearFluxes
from stellapy.data.read_fluxes import get_scaledLinearFluxes
from stellapy.data.read_fluxes import get_nonlinearFluxes 

# Read data from the omega file
from stellapy.data.read_omega import get_omega 
from stellapy.data.read_omega import get_linearOmega 
from stellapy.data.read_omega import write_linearOmega

# Read data in the input_file, profile_file, vmec_geo file, ...
from stellapy.data.read_wout import get_woutData
from stellapy.data.read_vmecgeo import get_vmecgeoData
from stellapy.data.read_profile import get_profileData
from stellapy.data.read_finalFields import get_finalFields
from stellapy.data.read_referenceUnits import get_referenceUnits 
from stellapy.data.read_inputParameters import read_inputParameters


#################################################################
#                        CREATE SIMULATIONS
#################################################################

def create_simulations(folders=None, input_files=None, ignore_resolution=True, number_variedVariables=1, folderIsExperiment=False):
    """ Will take multiple folders and input files and will group the input files
    together that only differ in (kx,ky) since they are in essence the same simulation. """
     
    # Make sure we have lists of PosixPath objects
    if isinstance(folders, str): folders = [folders]
    if isinstance(input_files, str): input_files = [input_files]
    if isinstance(folders, pathlib.PurePath): folders = [folders]
    if isinstance(input_files, pathlib.PurePath): input_files = [input_files]
    if input_files is not None and not isinstance(input_files[0], pathlib.PurePath):
        input_files = [pathlib.Path(i) for i in input_files]
    if folders is not None and not isinstance(folders[0], pathlib.PurePath):
        folders = [pathlib.Path(f) for f in folders]
           
    # If <input_files>=None read all input_files in <folders>
    if input_files == None:
        input_files = get_filesInFolder(folders, end=".in")
             
    # If <folders> and <input_files> are given they match
    elif folders!=None and input_files != None:
        input_files = [folders[i] / input_files[i] for i in range(len(folders))]
    
    # Only look at input files that have a netcdf output file
    input_files = remove_simulationsWithoutOutncFile(input_files) 
    if input_files == []:
        print("WARNING: THERE WERE NO NETCDF FILES FOUND!")
        
    # Group the simulations together and save them as simulations[id] = [input_files] 
    simulation_ids = {} 
    
    # Save the inputs that we deleted
    inputs = {}
    
    # Go through the PosixPaths of the input files
    for input_file in input_files:
        
        # Get the simulation identifier and the input parameters except the (kx,ky) section
        simulation_id   = get_simulationIdentifier(input_file)
        inputParameters = read_inputParameters(input_file)
        
        # Save the input parameters that we might delete
        inputs[input_file] = {}
        inputs[input_file]['kt_grids_range_parameters'] = inputParameters['kt_grids_range_parameters'].copy()
        inputs[input_file]['vpamu_grids_parameters'] = inputParameters['vpamu_grids_parameters'].copy()
        inputs[input_file]['zgrid_parameters'] = inputParameters['zgrid_parameters'].copy()
        inputs[input_file]['vmec_parameters'] = inputParameters['vmec_parameters'].copy()
        inputs[input_file]['stella_diagnostics_knobs'] = inputParameters['stella_diagnostics_knobs'].copy()
        inputs[input_file]['init_g_knobs'] = inputParameters['init_g_knobs'].copy()
        inputs[input_file]['knobs'] = inputParameters['knobs'].copy()
        
        # The (kx,ky) values are allowed to be different in the "same" simulation
        del inputParameters['kt_grids_range_parameters']
        
        # Remove other parameters that are allowed to differ
        del inputParameters['knobs']['mat_gen']
        del inputParameters['knobs']['nstep']
        del inputParameters['knobs']['delt_option']
        del inputParameters['init_g_knobs']#['ginit_option']
        del inputParameters['stella_diagnostics_knobs']
        if inputParameters['physics_flags']['nonlinear']==True: del inputParameters['knobs']['delt']
         
        # If ignore_resolution = False: the resolution can differ because convergence studies are performed
        # If ignore_resolution = True: the resolution can differ to make sure each mode is converged
        if ignore_resolution:
            del inputParameters['vmec_parameters']['nfield_periods']
            del inputParameters['zgrid_parameters']['nzed']
            del inputParameters['vpamu_grids_parameters']['nvgrid']
            del inputParameters['vpamu_grids_parameters']['nmu']
            if inputParameters['physics_flags']['nonlinear']==False: del inputParameters['knobs']['delt']
        
        # Assume we have a new simulation and it is in a new folder
        newInputFileIsUnique = True
        newFolder = {}
        for simulation in simulation_ids.keys():
            newFolder[simulation]=False
        
        # We can sort by "1 folder = 1 experiment"
        if number_variedVariables<0 or folderIsExperiment==True:
            for simulation in simulation_ids.keys():
                newFolder[simulation]=True
            parent_newInput = input_file.parent
            if ("run" in parent_newInput.name) and parent_newInput.name.split("run")[-1].isdigit():
                parent_newInput = parent_newInput.parent
            for simulation in simulation_ids.keys():
                old_input = simulation_ids[simulation]['input files'][0]
                parent_oldInput = old_input.parent
                if ("run" in parent_oldInput.name) and parent_oldInput.name.split("run")[-1].isdigit():
                    parent_oldInput = parent_oldInput.parent
                if parent_newInput == parent_oldInput:
                    newFolder[simulation] = False
            
        # Add the input file to the simulation if it has the same input parameters.
        # If the input parameters are different, create a new simulation for this input file
        if len(list(simulation_ids.keys()))!=0:
            for simulation in simulation_ids.keys():
                if newFolder[simulation]==False and simulation_ids[simulation]['input parameters'] == inputParameters: 
                    newInputFileIsUnique = False
                    simulation_ids[simulation]['input files'].append(input_file)
        if newInputFileIsUnique==True:
            simulation_ids[simulation_id] = {'input files' : [input_file], 'input parameters' : inputParameters}
     
    # For each simulation create a simulation object
    simulation_objects = []
    for simulation_id in simulation_ids.keys():
        input_files = simulation_ids[simulation_id]['input files']
        inputParameters = simulation_ids[simulation_id]['input parameters']
        simulation_objects.append(Simulation(simulation_id, input_files, inputParameters, inputs))
    return simulation_objects

    # Dont collapse header on the next line
    if True: return
 
#################################################################
#                        CLASS SIMULATION
#################################################################
 
class Simulation:
    ''' Saves all the information from a simulation run with stella.
    
    Notes
    -----
    A simulation is a combination of input_files that are identical but have different modes, 
    or that are identical but are restarted. These files need to be in the same folder!
        
    '''
     
    def __init__(self, simulation_id, input_files, inputParameters, inputs, Progress=None):
         
        # Save the information that defines the simulation
        self.id = simulation_id                      
        self.input_files = input_files               
        self.inputParameters = inputParameters  
        self.inputs = inputs
        self.formula_quasiLinearFluxes = 'per'
        
        # Get the parent folder of these simulations
        if len(self.input_files)>1:  self.parent = os.path.commonpath(self.input_files)
        if len(self.input_files)==1: self.parent = self.input_files[0].parent
        
        # Print information on the GUI (set from the GUI itself)
        self.Progress = Progress
        
        # Now read the simulation file or write it if it doesn't exist
        self.read_simulationFile()   
        
        # If we reran some modes with a different resolution we need to rewrite the lineardata!
        self.rewrite_simulationFile = False   
        
        # Get some default values for the markers and lines
        self.set_labelsLinesMarkers()           
        
        # Keep track of which modes we isolated for the kx or ky scan (see get_modesForAOneDimensionalScan)
        self.scan = None
        self.k_fixed = None
        self.k_range = None
        self.plotted_modes = None
        
        # Dont collapse header on the next line
        if True: return 

#################################################################
#                        METHODS
#################################################################

#============================================
#  SAVE SIMULATION DATA TO SIMULATION.INI
#============================================

    def read_simulationFile(self):
        ''' Check whether <simulation.ini> exists, if not write it. '''
        
        # Update the GUI
        if self.Progress: self.Progress.move(0,"Reading the simulation file.")
        
        # Create a configuration object and read the simulation file
        if len(self.input_files)>1:  
            parent = pathlib.Path("/".join(os.path.commonprefix(self.input_files).split('/')[0:-1]))
        if len(self.input_files)==1: 
            parent = self.input_files[0].parent
        self.simulation_file = configparser.ConfigParser() 
        self.simulation_path = parent / "simulation.ini"
        self.simulation_file.read(self.simulation_path)
        
        # The configuration file didn't exist when it was read
        if ("GENERAL" not in self.simulation_file):
            self.create_defaultSimulationFile()
 
        # If the file did exist make sure the input files and modes are up to date
        elif ("GENERAL" in self.simulation_file):
            input_files = [str(i) for i in self.input_files]
            input_files = "\t\n" + ("\t\n").join(input_files)
            modes_kx = [ str(k) for k in self.vec_kx ]
            modes_kx = "[" + ", \n\t\t".join([", ".join(modes_kx[i:i+5]) for i in range(0,len(modes_kx),5)]) + "]"
            modes_ky = [ str(k) for k in self.vec_ky ]
            modes_ky = "[" + ", \n\t\t".join([", ".join(modes_ky[i:i+5]) for i in range(0,len(modes_ky),5)]) + "]"
            self.simulation_file['GENERAL']['input_files'] = input_files
            self.simulation_file['GENERAL']['kx_modes'] = modes_kx
            self.simulation_file['GENERAL']['ky_modes'] = modes_ky
            self.simulation_file.write(open(self.simulation_path, 'w')) 
            
        # Save the removed kx and ky modes  
        self.removed_kx = self.simulation_file["GENERAL"]["removed_kx_modes"]
        self.removed_kx = self.removed_kx.split("[")[-1].split("]")[0].split(", ")
        self.removed_ky = self.simulation_file["GENERAL"]["removed_ky_modes"]
        self.removed_ky = self.removed_ky.split("[")[-1].split("]")[0].split(", ")
        self.removed_kx = [ float(kx) for kx in self.removed_kx ] if self.removed_kx != [''] else []
        self.removed_ky = [ float(ky) for ky in self.removed_ky ] if self.removed_ky != [''] else []
            
        # If the last modification data is different in the simulation file, 
        # then rewrite the lineardata/stability/convergence data
        if self.simulation_file['GENERAL']['last_modified'] != self.get_lastModification()\
        or self.simulation_file['GENERAL']['simulation_id'] != self.id:
            
            # Save the new latest modification time and the included modes
            self.simulation_file['GENERAL']['last_modified'] = self.get_lastModification()
            
            # Now rewrite the data in the simulation file
            if self.inputParameters['physics_flags']['nonlinear'] == False:
                print("WRITE AGAIN BECAUSE THE TIME CHANGED OR THE WRONG SIMULATION FILE WAS READ")
                self.simulation_file['GENERAL']['simulation_id'] = self.id
                self.write_unstableModes()
                write_linearOmega(self) 
        return
                    
    #------------------------------------------      
    def create_defaultSimulationFile(self):
        ''' Write a default configuration file and save it as "simulations.ini". '''
        
        # Rewrite the input_files and modes to make them look orderly
        input_files = [str(i) for i in self.input_files]
        input_files = "\t\n" + ("\t\n").join(input_files)
        modes_kx = [ str(k) for k in self.vec_kx ]
        modes_kx = "[" + ", \n\t\t".join([", ".join(modes_kx[i:i+5]) for i in range(0,len(modes_kx),5)]) + "]"
        modes_ky = [ str(k) for k in self.vec_ky ]
        modes_ky = "[" + ", \n\t\t".join([", ".join(modes_ky[i:i+5]) for i in range(0,len(modes_ky),5)]) + "]"
        
        # Remember which simulation this is   
        self.simulation_file['GENERAL'] = {
            'simulation_id' : self.id,\
            'last_modified' : self.get_lastModification(),\
            'kx_modes' : modes_kx,\
            'ky_modes' : modes_ky,\
            'removed_kx_modes' : "[]",\
            'removed_ky_modes' : "[]",\
            'input_files': input_files}
        
        # Style the simulation
        self.simulation_file['OVERRIDE COLORS AND FONTS'] = {
            'line_label'   : 'USE DEFAULT',\
            'line_style'   : 'USE DEFAULT',\
            'line_color'   : 'USE DEFAULT',\
            'marker_label' : 'USE DEFAULT',\
            'marker_style' : 'USE DEFAULT',\
            'marker_color' : 'USE DEFAULT'} 
        
        # Write the configuration file
        self.simulation_file.write(open(self.simulation_path, 'w'))
        return
        
    #--------------------------------
    def get_lastModification(self):
        ''' Calculate the date on which the files were last modified.
        # This allows us to rewrite data if a simulation has been rerun or added. '''
        folder = self.input_files[0].parent
        files = [ f for f in os.scandir(folder) \
                 if (not f.name.startswith('.') \
                 and not f.name.endswith('.h5') \
                 and not f.name.endswith('~') \
                 and not f.name=="simulation.ini" )]
        return time.ctime(max(entry.stat().st_mtime for entry in files if not entry.name.startswith('.'))) 

        # Dont collapse header on the next line
        if True: return
        
#============================================
#  CALCULATE MODE INDICES AND PROPERTIES
#============================================

    def get_indicesOfMode(self, scan, k_fixed, k_range):
        ''' Return (i_kx, i_ky) of the mode to access e.g. omega[i_kx, i_ky] = float '''
        
        # Scan over kx at fixed ky
        if scan=="kx":
            i_kx = self.vec_kx.index(k_range)
            i_ky = self.vec_ky.index(k_fixed)
            
        # Scan over ky at fixed kx
        if scan=="ky":
            i_kx = self.vec_kx.index(k_fixed)
            i_ky = self.vec_ky.index(k_range)

        return i_kx, i_ky
    
    #----------------------------------------
    def get_indicesOfModes(self, scan, k_fixed, modes):
        ''' Return a mask (i_kx, i_ky) of the modes to access e.g. omega[i_kxky] = array'''

        # Scan over kx at fixed ky
        if scan=="kx":
            i_kxky = np.outer(np.isin(self.vec_kx, modes), np.isin(self.vec_ky, k_fixed))
            
        # Scan over ky at fixed kx
        if scan=="ky":
            i_kxky = np.outer(np.isin(self.vec_kx, k_fixed), np.isin(self.vec_ky, modes))

        return i_kxky
    
    #--------------------------------
    def write_unstableModes(self):
        ''' Sort the modes by stable/unstable and converged/unconverged. '''
        
        # The section headers to save the data
        mode_header = 'STABLE MODE; CONVERGED MODE; PHI_END - PHI_START, RELATIVE FLUCTUATION OVER THE LAST 10% OF OMEGA AND GAMMA' 
        
        # If the section already existed, remove it
        if mode_header in self.simulation_file: self.simulation_file.remove_section(mode_header) 
        
        # Create the empty section
        self.simulation_file[mode_header] = {} 
            
        # Build a matrix to hold whether the mode is stable/unstable and converged/unconverged
        self.stable_modes    = np.empty((len(self.vec_kx), len(self.vec_ky))); self.stable_modes[:,:] = np.NaN
        self.converged_modes = np.empty((len(self.vec_kx), len(self.vec_ky))); self.converged_modes[:,:] = np.NaN
        
        # Go through the modes
        for i in range(len(self.vec_kx)):
            for j in range(len(self.vec_ky)):
            
                # Update the progress bar of the GUI
                if self.Progress: c = i*len(self.vec_ky) + j; length = len(self.vec_kx)*len(self.vec_ky)
                if self.Progress: self.Progress.move(c/length*100,"Determining the stability ("+str(c)+"/"+str(length)+")")
                
                # An unstable mode will grow in phi2 while a stable mode will decrease in phi2, therefore:
                # if phi_end < phi_start the mode is stable, if phi_end > phi_start, the mode is unstable.
                phi2_noNaN      = self.phi2_kxky[:,i,j][np.isfinite(self.phi2_kxky[:,i,j])]
                first_10percent = int(np.size(phi2_noNaN)*1/10)
                last_10percent  = int(np.size(phi2_noNaN)*9/10)
                phi_start = np.mean(phi2_noNaN[0:first_10percent])
                phi_end   = np.mean(phi2_noNaN[last_10percent:-1]) 
                 
                # Add <True> to <stable_modes> if the mode is stable, otherwise add <False>
                if phi_end - phi_start < 0: self.stable_modes[i,j] = True
                if phi_end - phi_start > 0: self.stable_modes[i,j] = False
                if len(phi2_noNaN) < 15:    self.stable_modes[i,j] = False
                
#                 # The above criteria doesn't always work correctly: stable modes have noise of gamma around zero!
#                 gamma_noNaN     = self.gamma_kxky[:,i,j][np.isfinite(self.gamma_kxky[:,i,j])] 
#                 last_10percent  = int(np.size(gamma_noNaN)*9/10) 
#                 gamma_10percent = gamma_noNaN[last_10percent:-1]
#                 sign_changes    = np.where(np.diff(np.sign(gamma_10percent)))[0]
#                 if len(sign_changes) >= : self.stable_modes[i,j] = True
                
                # Only if the mode is unstable, look if it is converged
                fluctuation_o = np.NaN; fluctuation_g = np.NaN
                if self.stable_modes[i,j] == False:
    
                    # Look at omega and gamma to determine whether the mode has converged
                    omega  = abs(self.omega_kxky)[:,i,j]
                    gamma  = abs(self.gamma_kxky)[:,i,j]
            
                    # Remove the nan and inifinite values
                    omega  = omega[np.isfinite(gamma)]
                    gamma  = gamma[np.isfinite(gamma)]
                    
                    # If the mode has converged it approaches a fixed value in time
                    # Consinder the mode converged if gamma/omega doesn't change more than 2% over the last 10% of time
                    did_yConverge = np.all( np.isclose(omega[int(len(omega)*0.90):], omega[-1], rtol=0.02) )
                    did_yConverge = np.all( np.isclose(gamma[int(len(gamma)*0.90):], gamma[-1], rtol=0.02) ) & did_yConverge
                    
                    # Add <True> to <converged_modes> if the mode has converged, otherwise add <False>
                    if did_yConverge == True:  self.converged_modes[i,j] = True
                    if did_yConverge == False: self.converged_modes[i,j] = False 
                    
                    # Add some information to the command prompt
                    if did_yConverge == False: 
                        kx = self.vec_kx[i]; ky = self.vec_ky[j]
                        fluctuation_o = (max(omega[int(len(omega)*0.90):])-min(omega[int(len(omega)*0.90):]))/omega[-1]*100
                        fluctuation_g = (max(gamma[int(len(gamma)*0.90):])-min(gamma[int(len(gamma)*0.90):]))/gamma[-1]*100
                        printv("\nThe simulation for (kx,kx) = ("+str(kx)+", "+str(ky)+") has not converged yet:")
                        printv("    The relative fluctuation of omega between (0.9*t, t) is:  "+str(fluctuation_o)+"%.")
                        printv("    The relative fluctuation of gamma between (0.9*t, t) is:  "+str(fluctuation_g)+"%.")
        
                # Save the linear data: [omega_last, omega_avg, omega_min, omega_max]
                kx = self.vec_kx[i];  ky = self.vec_ky[j]; 
                mode_indices = '(' + str(kx) + ", " + str(ky) + ')'
                mode_data = [ str(k) for k in [self.stable_modes[i,j], self.converged_modes[i,j], phi_end - phi_start, fluctuation_o, fluctuation_g] ] 
                self.simulation_file[mode_header][mode_indices] = "[" + ", ".join(mode_data) + "]"
    
        # Write the simulation file
        self.simulation_file.write(open(self.simulation_path, 'w'))
        return 
    
    #-----------------------------------
    def get_unstableModes(self):
        ''' Read whether the modes are stable/unstable and converged/unconverged. '''
        
        # Update the progress bar of the GUI
        if self.Progress: self.Progress.move(0,"Reading the stability.")
                 
        # If the stability isn't written in the simulation file then write it.
        mode_header = 'STABLE MODE; CONVERGED MODE; PHI_END-PHI_START; RELATIVE FLUCTUATION OVER THE LAST 10% OF OMEGA AND GAMMA' 
        if mode_header not in self.simulation_file:
            print("WRITE STABILITY AGAIN BECAUSE THE HEADER WAS MISSING")
            self.write_unstableModes()
            return
            
        # Build a matrix to hold whether the mode is stable/unstable and converged/unconverged
        self.stable_modes    = np.empty((len(self.vec_kx), len(self.vec_ky))); self.stable_modes[:,:] = np.NaN
        self.converged_modes = np.empty((len(self.vec_kx), len(self.vec_ky))); self.converged_modes[:,:] = np.NaN
        
        # Go through the modes
        try:
            for kx in self.vec_kx:
                for ky in self.vec_ky:
                    
                    # Get the mode indices  
                    i = self.vec_kx.index(kx)
                    j = self.vec_ky.index(ky)
                    
                    # Get the mode identifier for the simulation.ini file
                    mode_indices = '(' + str(kx) + ", " + str(ky) + ')'
                         
                    # Read the linear data: [omega_last, omega_avg, omega_min, omega_max]
                    mode_data = self.simulation_file[mode_header][mode_indices].split("[")[-1].split("]")[0].split(", ") 
                    mode_data = [ float(x) for x in mode_data ]
                
                    # Now put this data in the (kx,ky) matrixes
                    self.stable_modes[i,j]    = mode_data[0]
                    self.converged_modes[i,j] = mode_data[1] 
                    
        except:
            print("WRITE STABILITY AGAIN BECAUSE THE FILE WAS MISSING MODES")
            self.write_unstableModes()
                
        return 
    
    #-----------------------------------
    def get_modesForAOneDimensionalScan(self, scan, k_fixed, k_range, plotted_modes):
        ''' Reduces the vec_kx and vec_ky to vec_k depending on which axis we want to scan. '''
    
        # Check with the previous parameters to avoid recalculating
        if self.scan != scan or self.k_fixed != k_fixed or self.k_range != k_range:
            self.scan = scan 
            self.k_fixed = k_fixed
            self.k_range = k_range
            self.plotted_modes = None
            
            # Scan over kx at fixed ky
            if scan=="kx":
                self.vec_kAll    = self.vec_kx
                stable_modes     = self.stable_modes[:,self.vec_ky.index(k_fixed)]
                converged_modes  = self.converged_modes[:,self.vec_ky.index(k_fixed)]
            
            # Scan over ky at fixed kx
            if scan=="ky":
                self.vec_kAll    = self.vec_ky
                stable_modes     = self.stable_modes[self.vec_kx.index(k_fixed),:]
                converged_modes  = self.converged_modes[self.vec_kx.index(k_fixed),:]
                
            # Collect the modes that are stable/unstable and converged/unconverged
            self.vec_kStable      = list(np.array(self.vec_kAll)[stable_modes==1])
            self.vec_kUnstable    = list(np.array(self.vec_kAll)[stable_modes==0])
            self.vec_kConverge    = list(np.array(self.vec_kAll)[converged_modes==1])
            self.vec_kNotConverge = list(np.array(self.vec_kAll)[converged_modes==0])
                
            # Only look at the modes within the correct range
            self.vec_kStable      = [ k for k in self.vec_kStable if (k >= k_range[0] and k <= k_range[1])]
            self.vec_kUnstable    = [ k for k in self.vec_kUnstable if (k >= k_range[0] and k <= k_range[1])]
            self.vec_kConverge    = [ k for k in self.vec_kConverge if (k >= k_range[0] and k <= k_range[1])]
            self.vec_kNotConverge = [ k for k in self.vec_kNotConverge if (k >= k_range[0] and k <= k_range[1])]
    
        # For the plots we only plot one of the following groups of modes
        if self.plotted_modes != plotted_modes:
            if plotted_modes=="stable":      self.vec_k = self.vec_kStable
            if plotted_modes=="unstable":    self.vec_k = self.vec_kUnstable
            if plotted_modes=="converged":   self.vec_k = self.vec_kConverge
            if plotted_modes=="unconverged": self.vec_k = self.vec_kNotConverge
            self.plotted_modes = plotted_modes
        
        # Return the total amount of modes and those to plot
        return self.vec_k
        if True: return

#################################################################
#                        READ RAW DATA
#################################################################    

#=======================================
#        DATA IN THE NETCDF FILE
#=======================================   
    @calculate_attributeWhenReadFirstTime
    def vec_kx(self):           get_netcdfDimensions(self);     return self.vec_kx
    @calculate_attributeWhenReadFirstTime
    def vec_ky(self):           get_netcdfDimensions(self);     return self.vec_ky 
    @calculate_attributeWhenReadFirstTime
    def vec_kxPerFile(self):    get_netcdfDimensions(self);     return self.vec_kxPerFile 
    @calculate_attributeWhenReadFirstTime
    def vec_kyPerFile(self):    get_netcdfDimensions(self);     return self.vec_kyPerFile 
    @calculate_attributeWhenReadFirstTime
    def dim_kx(self):           get_netcdfDimensions(self);     return self.dim_kx
    @calculate_attributeWhenReadFirstTime
    def dim_ky(self):           get_netcdfDimensions(self);     return self.dim_ky
    @calculate_attributeWhenReadFirstTime
    def dim_time(self):         get_netcdfDimensions(self);     return self.dim_time
    @calculate_attributeWhenReadFirstTime
    def dim_z(self):            get_netcdfDimensions(self);     return self.dim_z
    @calculate_attributeWhenReadFirstTime 
    def dim_species(self):      get_netcdfDimensions(self);     return self.dim_species   
    @calculate_attributeWhenReadFirstTime
    def dim_timePerFile(self):  get_netcdfDimensions(self);     return self.dim_timePerFile   
    @calculate_attributeWhenReadFirstTime
    def dim_zPerFile(self):     get_netcdfDimensions(self);     return self.dim_zPerFile   
    @calculate_attributeWhenReadFirstTime
    def dim_kxPerFile(self):    get_netcdfDimensions(self);     return self.dim_kxPerFile 
    @calculate_attributeWhenReadFirstTime
    def dim_kyPerFile(self):    get_netcdfDimensions(self);     return self.dim_kyPerFile
    @calculate_attributeWhenReadFirstTime
    def time_kxky(self):        get_netcdfTimeVsKxKy(self);     return self.time_kxky
    @calculate_attributeWhenReadFirstTime
    def phi2_kxky(self):        get_netcdfPhi2VsKxKy(self);     return self.phi2_kxky
    @calculate_attributeWhenReadFirstTime
    def z_kxky(self):           get_netcdfZVsKxKy(self);        return self.z_kxky
    @calculate_attributeWhenReadFirstTime 
    def pflx_kxky(self):        get_netcdfFluxesKxKy(self);     return self.pflx_kxky
    @calculate_attributeWhenReadFirstTime 
    def qflx_kxky(self):        get_netcdfFluxesKxKy(self);     return self.qflx_kxky
    @calculate_attributeWhenReadFirstTime 
    def vflx_kxky(self):        get_netcdfFluxesKxKy(self);     return self.vflx_kxky
    @calculate_attributeWhenReadFirstTime 
    def density_kxky(self):     get_netcdfDensity(self);        return self.density_kxky
     
#=======================================
#        DATA IN THE OMEGA FILE
#=======================================
    @calculate_attributeWhenReadFirstTime 
    def omega_kxky(self):       get_omega(self);                return self.omega_kxky
    @calculate_attributeWhenReadFirstTime 
    def gamma_kxky(self):       get_omega(self);                return self.gamma_kxky

#=======================================
#        DATA IN THE VMECGEO FILE
#=======================================
    @calculate_attributeWhenReadFirstTime 
    def ref_a(self):            get_vmecgeoData(self);          return self.ref_a
    @calculate_attributeWhenReadFirstTime 
    def ref_B(self):            get_vmecgeoData(self);          return self.ref_B
        
#=======================================
#        DATA IN THE WOUT FILE
#=======================================
    @calculate_attributeWhenReadFirstTime 
    def sign_B(self):           get_woutData(self);             return self.sign_B    
    @calculate_attributeWhenReadFirstTime 
    def iota(self):             get_woutData(self);             return self.iota
      
#=======================================
#        PROFILE FILE
#=======================================
    @calculate_attributeWhenReadFirstTime 
    def prof_n(self):           get_profileData(self);          return self.prof_n
    @calculate_attributeWhenReadFirstTime 
    def prof_T(self):           get_profileData(self);          return self.prof_T
    
#=======================================
#        DATA IN THE FLUXES FILE
#=======================================
    @calculate_attributeWhenReadFirstTime 
    def p_flux_kxky(self):      get_linearFluxes(self);         return self.p_flux_kxky
    @calculate_attributeWhenReadFirstTime 
    def q_flux_kxky(self):      get_linearFluxes(self);         return self.q_flux_kxky
    @calculate_attributeWhenReadFirstTime 
    def v_flux_kxky(self):      get_linearFluxes(self);         return self.v_flux_kxky
    @calculate_attributeWhenReadFirstTime 
    def vec_time(self):         get_nonlinearFluxes(self);      return self.vec_time
    @calculate_attributeWhenReadFirstTime 
    def p_flux(self):           get_nonlinearFluxes(self);      return self.p_flux
    @calculate_attributeWhenReadFirstTime 
    def q_flux(self):           get_nonlinearFluxes(self);      return self.q_flux
    @calculate_attributeWhenReadFirstTime 
    def v_flux(self):           get_nonlinearFluxes(self);      return self.v_flux
    @calculate_attributeWhenReadFirstTime 
    def restart_times(self):    get_nonlinearFluxes(self);      return self.restart_times
     
#=======================================
#          FINAL FIELDS
#=======================================   
    @calculate_attributeWhenReadFirstTime 
    def phi_real(self):         get_finalFields(self);      return self.phi_real
    @calculate_attributeWhenReadFirstTime 
    def phi_imag(self):         get_finalFields(self);      return self.phi_imag
    @calculate_attributeWhenReadFirstTime 
    def phi2(self):             get_finalFields(self);      return self.phi2
    @calculate_attributeWhenReadFirstTime 
    def z_poloidal(self):       get_finalFields(self);      return self.z_poloidal
    @calculate_attributeWhenReadFirstTime 
    def zeta(self):             get_finalFields(self);      return self.zeta
    @calculate_attributeWhenReadFirstTime 
    def dont_foldHeader(self):
        if True: return

#################################################################
#                  PRODUCE PROCESSED DATA
#################################################################

#=======================================
#        REFERENCE UNITS
#=======================================
    @calculate_attributeWhenReadFirstTime 
    def referenceUnits(self):  get_referenceUnits(self);       return self.referenceUnits

#=======================================
#     CALCULATE QUASI LINEAR FLUXES 
#=======================================
    @calculate_attributeWhenReadFirstTime 
    def p_fluxQL_kxky(self):   get_scaledLinearFluxes(self);   return self.p_fluxQL_kxky
    
#=======================================
#        CALCULATE LINEAR DATA
#=======================================
    @calculate_attributeWhenReadFirstTime 
    def omega_last(self):      get_linearOmega(self);          return self.omega_last
    @calculate_attributeWhenReadFirstTime 
    def gamma_last(self):      get_linearOmega(self);          return self.gamma_last
    @calculate_attributeWhenReadFirstTime 
    def omega_avg(self):       get_linearOmega(self);          return self.omega_avg
    @calculate_attributeWhenReadFirstTime 
    def gamma_avg(self):       get_linearOmega(self);          return self.gamma_avg
    @calculate_attributeWhenReadFirstTime 
    def p_flux_last(self):     get_quasiLinearFluxes(self);    return self.p_flux_last
    @calculate_attributeWhenReadFirstTime 
    def q_flux_last(self):     get_quasiLinearFluxes(self);    return self.q_flux_last
    @calculate_attributeWhenReadFirstTime 
    def v_flux_last(self):     get_quasiLinearFluxes(self);    return self.v_flux_last
    @calculate_attributeWhenReadFirstTime 
    def p_fluxQL_last(self):   get_quasiLinearFluxes(self);    return self.p_fluxQL_last
    @calculate_attributeWhenReadFirstTime 
    def p_flux_avg(self):      get_quasiLinearFluxes(self);    return self.p_flux_avg
    @calculate_attributeWhenReadFirstTime 
    def q_flux_avg(self):      get_quasiLinearFluxes(self);    return self.q_flux_avg
    @calculate_attributeWhenReadFirstTime 
    def v_flux_avg(self):      get_quasiLinearFluxes(self);    return self.v_flux_avg
    @calculate_attributeWhenReadFirstTime 
    def p_fluxQL_avg(self):    get_quasiLinearFluxes(self);    return self.p_fluxQL_avg
     
#=======================================
#    STABLE AND UNSTABLE MODES
#=======================================
    @calculate_attributeWhenReadFirstTime 
    def stable_modes(self):    self.get_unstableModes();       return self.stable_modes
    @calculate_attributeWhenReadFirstTime 
    def converged_modes(self): self.get_unstableModes();       return self.converged_modes
    @calculate_attributeWhenReadFirstTime 
    def dont_foldHeader2(self):
        if True: return
        
#################################################################
#                  STYLE AND CONFIGURATION
#################################################################

#=======================================
#       STYLE OF THE SIMULATION
#=======================================

    def set_labelsLinesMarkers(self, line_label=None, marker_label=None, line_style=None, line_color=None, marker_style=None, marker_color=None):
        
        # If the styles are given apply them
        if line_label is not None:
            self.line_label = line_label
            self.marker_label = marker_label
            self.line_style = line_style
            self.line_color = line_color
            self.marker_style = marker_style
            self.marker_color = marker_color
        
        # Otherise load some basic styles
        else: 
            self.line_label,\
            self.marker_label,\
            self.line_style,\
            self.line_color,\
            self.marker_style,\
            self.marker_color = load_labelsLinesMarkers([self.id], [self.id])
        
        # If the simulation was ran before, there is a text file, then read the manually set styles
        style = self.simulation_file['OVERRIDE COLORS AND FONTS']
        if style['line_label']   != 'USE DEFAULT':  self.line_label   = style['line_label']
        if style['line_color']   != 'USE DEFAULT':  self.line_color   = style['line_color']
        if style['line_style']   != 'USE DEFAULT':  self.line_style   = style['line_style']
        if style['marker_label'] != 'USE DEFAULT':  self.marker_label = style['marker_label']
        if style['marker_color'] != 'USE DEFAULT':  self.marker_color = style['marker_color']
        if style['marker_style'] != 'USE DEFAULT':  self.marker_style = style['marker_style']

        # Save the style of the simulation as an indicator for manual inputs
        self.simulation_file['COLORS AND FONTS'] = {
            'line_label'   : self.line_label,\
            'line_style'   : self.line_style,\
            'line_color'   : self.line_color,\
            'marker_label' : self.marker_label,\
            'marker_style' : self.marker_style,\
            'marker_color' : self.marker_color} 
        
        # Write the new data to the configuration file
        self.simulation_file.write(open(self.simulation_path, 'w'))

        # Dont collapse header on the next line
        if True: return



    