
import numpy as np
from stellapy.utils.decorators import verbose_wrapper, printv
from stellapy.data.load_defaultInputParameters import load_defaultInputParameters
from stellapy.data.get_gridDivisionsAndSize import get_jtwist, get_dkx, get_shat, get_naky, get_nakx

#=========================================
# READ ALL PARAMETERS FROM THE INPUT FILE
#=========================================

@verbose_wrapper
def read_inputParameters(input_file):
    ''' Read "*.in" file and return all the input parameters.
    
    Parameters
    ----------
    input_file : PosixPath
        Absolute path of the *.in file that needs to be read.
    
    Returns
    -------
    inputParameters = dict[knobs][variable]
        Returns the namelists from the stella code with their values.
    
    
    Knobs and variables
    -------------------
    zgrid_parameters:  
        nzed; nperiod; ntubes; boundary_option; zed_equal_arc; shat_zero; nzgrid"
        
    geo_knobs:
       geo_option; overwrite_bmag; overwrite_gradpar; overwrite_gds2; overwrite_gds21; overwrite_gds22; 
       overwrite_gds23; overwrite_gds24; overwrite_gbdrift; overwrite_cvdrift; overwrite_gbdrift0; geo_file
       
    vmec_parameters
        vmec_filename; alpha0; zeta_center; nfield_periods; torflux; surface_option;
        verbose; zgrid_scalefac; zgrid_refinement_factor 
        
    parameters
        beta; vnew_ref; rhostar; zeff; tite; nine 
        
    vpamu_grids_parameters
        nvgrid; vpa_max; nmu; vperp_max; equally_spaced_mu_grid 
        
    dist_fn_knobs
        adiabatic_option 
        
    time_advance_knobs
        explicit_option; xdriftknob; ydriftknob; wstarknob; flip_flo

    kt_grids_knobs
        grid_option 
        
    kt_grids_box_parameters
        nx; ny; dkx; dky; jtwist; y0; naky; nakx
        
    kt_grids_range_parameters
        nalpha; naky; nakx; aky_min; aky_max; akx_min; akx_max; theta0_min; theta0_max
        
    physics_flags
        full_flux_surface; include_mirror; nonlinear; include_parallel_nonlinearity; include_parallel_streaming

    init_g_knobs
        tstart; scale; ginit_option; width0; refac; imfac; den0; par0; tperp0; den1; upar1; tpar1; tperp1; 
        den2; upar2; tpar2; tperp2; phiinit; zf_init; chop_side; left; even; restart_file; restart_dir; read_many
        
    knobs
        nstep; delt; fapar; fbpar; delt_option; zed_upwind; vpa_upwind; time_upwind; avail_cpu_time; cfl_cushion
        delt_adjust; mat_gen; mat_read; fields_kxkyz; stream_implicit; mirror_implicit; driftkinetic_implicit
        mirror_semi_lagrange; mirror_linear_interp; maxwellian_inside_zed_derivative; stream_matrix_inversion 

    species_knobs
        nspec; species_option 

    species_parameters_1; species_parameters_2; species_parameters_3; ...
        z; mass; dens; temp; tprim; fprim; d2ndr2; d2Tdr2; type 

    stella_diagnostics_knobs
        nwrite; navg; nmovie; nsave; save_for_restart; write_omega; write_phi_vs_time; write_gvmus
        write_gzvs; write_kspectra; write_moments; flux_norm; write_fluxes_kxky 

    millergeo_parameters
        rhoc; rmaj; shift; qinp; shat; kappa; kapprim; tri; triprim; rgeo; betaprim; betadbprim
        d2qdr2; d2psidr2; nzed_local; read_profile_variation; write_profile_variation 

    layouts_knobs
        xyzs_layout; vms_layout 

    neoclassical_input
        include_neoclassical_terms; nradii; drho; neo_option
        
    sfincs_input
        read_sfincs_output_from_file; nproc_sfincs; irad_min; irad_max; calculate_radial_electric_field
        includeXDotTerm; includeElectricFieldTermInXiDot; magneticDriftScheme; includePhi1
        includePhi1InKineticEquation; geometryScheme; VMECRadialOption; equilibriumFile
        coordinateSystem; inputRadialCoordinate; inputRadialCoordinateForGradients; aHat
        psiAHat; Delta; nu_n; dPhiHatdrN; Er_window; nxi; nx; Ntheta; Nzeta
    '''	
    
    # Read the *in file and get the parameters for species_1
    input_data   = open(input_file, 'r')
    input_text   = input_data.read().replace(' ', '')

    # Initiate the dictionary: load the default parameters
    input_parameters = load_defaultInputParameters()

    # Get the number of species
    nspec = read_integerInput(input_text, 'nspec')
    if nspec=="USE DEFAULT": nspec=1
    
    # Add more default species if nspec>2
    for i in range(2, nspec+1):
        knob = "species_parameters_" + str(i)
        keys =  list(input_parameters["species_parameters_1"].keys())
        input_parameters[knob] = {}
        for key in keys:
            input_parameters[knob][key] = input_parameters["species_parameters_1"][key]

    # Overwrite the default values if they have been changed in the <input files>
    knobs = list(input_parameters.keys()); knobs.remove(" ")
    for knob in knobs:
        try: # Only look at the text of the specific knob
            input_text_knob = input_text.split("&"+knob)[1].split("/")[0]
            parameters = list(input_parameters[knob].keys())
            for parameter in parameters:
                default_value = input_parameters[knob][parameter]
                input_parameters[knob][parameter] = read_parameterFromInputFile(input_text_knob, parameter, default_value)
        except: # If the knob is not present just pass
            pass

    # Fill in the species know for the adiabatic electrons: custom knob for the GUI
    if input_parameters["species_knobs"]["nspec"] == 1:
        input_parameters["species_parameters_a"]["nine"] = input_parameters["parameters"]["nine"]
        input_parameters["species_parameters_a"]["tite"] = input_parameters["parameters"]["tite"]
        input_parameters["species_parameters_a"]["dens"] = round(1/input_parameters["parameters"]["nine"], 4)
        input_parameters["species_parameters_a"]["temp"] = round(1/input_parameters["parameters"]["tite"], 4)

    # Calculate some extra variables
    input_parameters["vmec_parameters"]["rho"] = np.sqrt(input_parameters["vmec_parameters"]["torflux"] )

    # Now calculate some values manually
    if input_parameters["sfincs_input"]["irad_min"] == "-nradii/2":
        nradii = input_parameters["neoclassical_input"]["nradii"]
        input_parameters["sfincs_input"]["irad_min"] = -nradii/2 
        input_parameters["sfincs_input"]["irad_max"] = nradii/2  
    else:
        printv("WARNING: <irad_min> was set by the input file but it should be calculated indirectly through <nradii>.")

    if input_parameters["zgrid_parameters"]["nzgrid"] == "nzed/2 + (nperiod-1)*nzed":
        nzed    = input_parameters["zgrid_parameters"]["nzed"]
        nperiod = input_parameters["zgrid_parameters"]["nperiod"]
        input_parameters["zgrid_parameters"]["nzgrid"] = nzed/2 + (nperiod-1)*nzed
    else:
        printv("WARNING: <nzgrid> was set by the input file but it should be calculated indirectly through <nzed> and <nperiod>.")
 
    if input_parameters["vmec_parameters"]["zgrid_scalefac"] == "2.0 if zed_equal_arc else 1.0":
        zed_equal_arc = input_parameters["zgrid_parameters"]["zed_equal_arc"]
        input_parameters["vmec_parameters"]["zgrid_scalefac"]           = 2.0 if zed_equal_arc else 1.0
        input_parameters["vmec_parameters"]["zgrid_refinement_factor"]  = 4 if zed_equal_arc else 1
    else:
        printv("WARNING: <zgrid_scalefac> was set by the input file but it should be calculated indirectly through <zed_equal_arc>.")

    if input_parameters["physics_flags"]["nonlinear"] == True:
        if input_parameters["kt_grids_box_parameters"]["naky"] == "(ny-1)/3 + 1":
            ny = input_parameters["kt_grids_box_parameters"]["ny"]
            nx = input_parameters["kt_grids_box_parameters"]["nx"]
            input_parameters["kt_grids_box_parameters"]["naky"] = get_naky(ny)
            input_parameters["kt_grids_box_parameters"]["nakx"] = get_nakx(nx)
        if input_parameters["kt_grids_box_parameters"]["y0"] == -1.0:
            if input_parameters["physics_flags"]["full_flux_surface"] == False:
                print("WARNING: When simulating a flux tube, y0 needs to be set in the input file.")
            if input_parameters["physics_flags"]["full_flux_surface"] == True:
                input_parameters["kt_grids_box_parameters"]["y0"] = "1./(rhostar*geo_surf%rhotor)"            

        # Calculate some extra variables for nonlinear runs, note that svalue=psitor=torflux??
        nfield      = input_parameters["vmec_parameters"]["nfield_periods"]
        y0          = input_parameters["kt_grids_box_parameters"]["y0"]
        vmec_file   = input_parameters["vmec_parameters"]["vmec_filename"]
        psitor      = input_parameters["vmec_parameters"]["torflux"]
        nx          = input_parameters["kt_grids_box_parameters"]["nx"]
        input_parameters["kt_grids_box_parameters"]["dky"] = round(1./y0,4) # get the grid spacing in ky and then in kx using twist-and-shift BC
        try:
            vmec_file = input_file.parent + "/" + vmec_file
            if input_parameters["kt_grids_box_parameters"]["jtwist"] == -1:
                input_parameters["kt_grids_box_parameters"]["jtwist"] = get_jtwist(vmec_file, psitor, nfield)
            input_parameters["kt_grids_box_parameters"]["dkx"] = get_dkx(vmec_file, psitor, nfield, y0)
            input_parameters["kt_grids_box_parameters"]["shat"] = get_shat(vmec_file, psitor)
        except:
            try:
                vmec_file = str(vmec_file).split(".nc")[0]+".h5"
                if input_parameters["kt_grids_box_parameters"]["jtwist"] == -1:
                    input_parameters["kt_grids_box_parameters"]["jtwist"] = round(get_jtwist(vmec_file, psitor, nfield), 4)
                input_parameters["kt_grids_box_parameters"]["dkx"] = round(get_dkx(vmec_file, psitor, nfield, y0), 4)
                input_parameters["kt_grids_box_parameters"]["shat"] = round(get_shat(vmec_file, psitor), 4)
            except:
                printv('WARNING: vmec file was not found to calculate jtwist.')

    if input_parameters["physics_flags"]["nonlinear"] == False:
        for key in input_parameters["kt_grids_box_parameters"].keys():
            input_parameters["kt_grids_box_parameters"][key] = "Linear simulation"

    input_data.close()

    # Return the input_parameters
    return input_parameters
    if True: return

#============================================
# READ SPECIFIC PARAMETER FROM THE INPUT FILE
#=============================================

def read_parameterFromInputFile(input_text, parameter, default_value):
    ''' Read the variables if they exist and convert them to the correct datatype. ''' 

    if default_value == None:
        input_value = read_floatInput(input_text, parameter)
        return input_value if input_value != "USE DEFAULT" else default_value        

    if isinstance(default_value, bool):
        input_value = read_booleanInput(input_text, parameter)
        return input_value if input_value != "USE DEFAULT" else default_value

    if isinstance(default_value, int):
        input_value = read_integerInput(input_text, parameter)
        return input_value if input_value != "USE DEFAULT" else default_value

    if isinstance(default_value, float):
        input_value = read_floatInput(input_text, parameter)
        return input_value if input_value != "USE DEFAULT" else default_value

    if isinstance(default_value, str):
        input_value = read_stringInput(input_text, parameter)
        return input_value if input_value != "USE DEFAULT" else default_value
    return

#--------------------------------------------
def read_integerInput(input_text, variable):
    ''' Read the [variable] from the "*.in" file, if it doesn't exist, return NaN. '''
    if variable not in input_text:
        return "USE DEFAULT"
    try:    
        return int(input_text.split(variable)[1].split('=')[1].split('\n')[0])
    except: 
        print("WARNING: ", variable, "is not an integer so failed to read input file.")
        return np.NaN
    return

#--------------------------------------------
def read_floatInput(input_text, variable):
    ''' Read the [variable] from the "*.in" file, if it doesn't exist, return NaN. '''
    if variable not in input_text:
        return "USE DEFAULT"
    try:    
        return float(input_text.split(variable)[1].split('=')[1].split('\n')[0])
    except: 
        print("WARNING: ", variable, "is not a float so failed to read input file.")
        return np.NaN
    return

#--------------------------------------------
def read_booleanInput(input_text, variable):
    ''' Read the [variable] from the "*.in" file and convert .false. to False '''
    if variable not in input_text:
        return "USE DEFAULT"
    try:    
        value = input_text.split(variable)[1].split('=')[1].split('\n')[0]
    except: 
        print("WARNING: ", variable, "is not a boolean so failed to read input file.")
        value = np.NaN
    if value == ".false.": value = False
    if value == ".true." : value = True
    return value

#--------------------------------------------
def read_stringInput(input_text, variable):
    ''' Read the [variable] from the "*.in" file and encode it to a np.string_ or b'' string '''
    if variable not in input_text:
        return "USE DEFAULT"
    try:    
        value = str(input_text.split(variable)[1].split('=')[1].split('\n')[0])
    except: 
        print("WARNING: ", variable, "is not a string so failed to read input file.")
        value = str(np.NaN)
    value = value.replace(" ","")
    value = value.replace("'","")
    value = value.replace('"','') 
    return value
