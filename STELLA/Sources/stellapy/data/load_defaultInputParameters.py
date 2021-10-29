

def load_defaultInputParameters():
    ''' Load default input parameters of the stella code.
    
    Returns
    -------
    inputParameters = dict[knobs][variable]
        Returns the namelists from the stella code with their default values.
        
    Knobs and variables
    -------------------
    zgrid_parameters:  
        nzed; nperiod; ntubes; boundary_option; zed_equal_arc; shat_zero; nzgrid"
        
    geo_knobs:
       geo_option; overwrite_bmag; overwrite_gradpar; overwrite_gds2; overwrite_gds21; overwrite_gds22; q_as_x 
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
        tstart; scale; ginit_option; width0; refac; imfac; den0; par0; tperp0; den1; upar1; tpar1; tperp1; kxmin; kxmax;
        den2; upar2; tpar2; tperp2; phiinit; zf_init; chop_side; left; even; restart_file; restart_dir; read_many
        
    knobs
        nstep; delt; fapar; fbpar; delt_option; zed_upwind; vpa_upwind; time_upwind; avail_cpu_time; cfl_cushion
        delt_adjust; mat_gen; mat_read; fields_kxkyz; stream_implicit; mirror_implicit; driftkinetic_implicit
        mirror_semi_lagrange; mirror_linear_interp; maxwellian_inside_zed_derivative; stream_matrix_inversion 

    species_knobs
        nspec; species_option 

    species_parameters_1; species_parameters_2; species_parameters_3; ...
        z; mass; dens; tprim; fprim; d2ndr2; d2Tdr2; type 

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

    # Initiate dictionaries
    inputParameters = {} 

    inputParameters["zgrid_parameters"] = { 
            "nzed"              : 24        ,\
            "nperiod"           : 1         ,\
            "ntubes"            : 1         ,\
            "boundary_option"   : 'default' ,\
            "zed_equal_arc"     : False     ,\
            "shat_zero"         : 1.e-5     ,\
            "nzgrid"            : "nzed/2 + (nperiod-1)*nzed"} # Derived input variable

    inputParameters["geo_knobs"] = {
            "geo_option"        : 'local'   ,\
            "overwrite_bmag"    : False     ,\
            "overwrite_gradpar" : False     ,\
            "overwrite_gds2"    : False     ,\
            "overwrite_gds21"   : False     ,\
            "overwrite_gds22"   : False     ,\
            "overwrite_gds23"   : False     ,\
            "overwrite_gds24"   : False     ,\
            "overwrite_gbdrift" : False     ,\
            "overwrite_cvdrift" : False     ,\
            "overwrite_gbdrift0": False     ,\
            "geo_file"          :'input.geometry',\
            "q_as_x"            :'radial_variation'} # New parameter 17/12/2020; True by default in radial variation runs

    inputParameters["vmec_parameters"] = { #9
            "vmec_filename"     : 'equilibria/wout_w7x_standardConfig.nc',\
            "alpha0"            : 0.0       ,\
            "zeta_center"       : 0.0       ,\
            "nfield_periods"    : -1.0      ,\
            "torflux"           : 0.6354167 ,\
            "surface_option"    : 0         ,\
            "verbose"           : True       ,\
            "zgrid_scalefac"    : "2.0 if zed_equal_arc else 1.0",\
            "zgrid_refinement_factor" : "4 if zed_equal_arc else 1"}

    inputParameters["parameters"] = {
            "beta"              : 0.0       ,\
            "vnew_ref"          : -1.0      ,   # various input options will override this value if it is negative
            "rhostar"           : -1.0      ,   # = m_ref * vt_ref / (e * B_ref * a_ref), with refs in SI
            "zeff"              : 1.0       ,\
            "tite"              : 1.0       ,\
            "nine"              : 1.0       ,\
            "g_exb"             : 0.0       ,   # New parameter 17/12/2020
            "g_exbfac"          : 1.0       ,   # New parameter 17/12/2020      
            "omprimfac"         : 1.0}          # New parameter 17/12/2020

    inputParameters["vpamu_grids_parameters"] = {
           "nvgrid"             : 24        ,\
           "vpa_max"            : 3.0       ,\
           "nmu"                : 12        ,\
           "vperp_max"          : 3.0       ,\
           "equally_spaced_mu_grid" : False  }

    inputParameters["dist_fn_knobs"] = {
           "adiabatic_option"   : 'default'}

    inputParameters["time_advance_knobs"] = {
           "explicit_option"    : 'default' ,\
           "xdriftknob"         : 1.0       ,\
           "ydriftknob"         : 1.0       ,\
           "wstarknob"          : 1.0       ,\
           "flip_flop"          : False     }

    inputParameters["kt_grids_knobs"] = {
           "grid_option" : 'default'}

    inputParameters["kt_grids_box_parameters"] = {
            "nx"                : 1         ,\
            "ny"                : 1         ,\
            "dkx"               : -1.0      ,           # Derived input variable
            "dky"               : -1.0      ,           # Derived input variable
            "shat"              : -1.0      ,           # Derived input variable
            "jtwist"            : -1.0      ,           # Note that jtwist and y0 will possibly be modified later if they are negative
            "jtwistfac"         : 1.0       ,           # New parameter 17/12/2020
            "centered_in_rho"   : True      ,           # New parameter 17/12/2020
            "periodic_variation": False     ,           # New parameter 17/12/2020
            "y0"                : -1.0      ,\
            "naky"              : "(ny-1)/3 + 1",       # Derived input variable: get the number of de-aliased modes in y and x
            "nakx"              : "2*((nx-1)/3) +  1"}  # Derived input variable: get the number of de-aliased modes in y and x

    inputParameters["kt_grids_range_parameters"] = {
            "nalpha"            : 1         ,\
            "naky"              : 1         ,\
            "nakx"              : 1         ,\
            "aky_min"           : 0.0       ,\
            "aky_max"           : 0.0       ,\
            "akx_min"           : 0.0       ,           # set these to be nonsense values so can check later if they've been set
            "akx_max"           : -1.0      ,           # set these to be nonsense values so can check later if they've been set
            "theta0_min"        : 0.0       ,           # set these to be nonsense values so can check later if they've been set
            "theta0_max"        : -1.0      }           # set these to be nonsense values so can check later if they've been set

    inputParameters["physics_flags"] = {
           "full_flux_surface"  :  False    ,\
           "include_mirror"     :  True     ,\
           "nonlinear"          :  False    ,\
           "include_parallel_nonlinearity" : False  ,\
           "radial_variation"   :  False    ,    # New parameter 17/12/2020
           "include_parallel_streaming" : True,  
           "include_pressure_variation" : True,  # New parameter 17/12/2020
           "include_geometric_variation" : True} # New parameter 17/12/2020

    inputParameters["init_g_knobs"] = { #28
            "tstart"            : 0.        ,\
            "scale"             : 1.0       ,\
            "ginit_option"      : "default" ,\
            "width0"            : -3.5      ,\
            "refac"             : 1.        ,\
            "imfac"             : 0.        ,\
            "den0"              : 1.        ,\
            "upar0"             : 0.        ,\
            "tpar0"             : 0.        ,\
            "tperp0"            : 0.        ,\
            "den1"              : 0.        ,\
            "upar1"             : 0.        ,\
            "tpar1"             : 0.        ,\
            "tperp1"            : 0.        ,\
            "den2"              : 0.        ,\
            "upar2"             : 0.        ,\
            "tpar2"             : 0.        ,\
            "tperp2"            : 0.        ,\
            "phiinit"           : 1.0       ,\
            "zf_init"           : 1.0       ,\
            "chop_side"         : False     ,\
            "left"              : True      ,\
            "even"              : True      ,\
            "restart_file"      : "trim(run_name)//'.nc'",\
            "restart_dir"       : "./"      ,\
            "read_many"         : None      ,\
            "kx_min"            : 0.        ,   # New parameter 17/12/2020
            "kx_max"            : 1.e100     }  # New parameter 17/12/2020

    inputParameters["knobs"] = { #22
           "nstep"              : None      ,\
           "delt"               : None      ,\
           "fphi"               : 1.0       ,\
           "fapar"              : 1.0       ,\
           "fbpar"              : -1.0      ,\
           "delt_option"        : 'default' ,\
           "zed_upwind"         : 0.02      ,\
           "vpa_upwind"         : 0.02      ,\
           "time_upwind"        : 0.02      ,\
           "avail_cpu_time"     : 1.e10     ,\
           "cfl_cushion"        : 0.5       ,\
           "delt_adjust"        : 2.0       ,\
           "mat_gen"            : True      ,\
           "mat_read"           : False     ,\
           "fields_kxkyz"       : False     ,\
           "stream_implicit"    : True      ,\
           "mirror_implicit"    : True      ,\
           "driftkinetic_implicit" : False  ,\
           "mirror_semi_lagrange" : True    ,\
           "mirror_linear_interp" : False   ,\
           "maxwellian_inside_zed_derivative" : False  ,\
           "stream_matrix_inversion" : False,\
           "rng_seed"           : -1        , # New parameter 17/12/2020  !negative values use current time as seed
           "ky_solve_radial"    : 0         , # New parameter 17/12/2020
           "ky_solve_real"      : False}      # New parameter 17/12/2020


    inputParameters["species_knobs"] = {
           "nspec"              : 2         ,\
           "species_option"     : 'stella'  ,\
           "read_profile_variation" : False , # New parameter 17/12/2020
           "write_profile_variation" : False} # New parameter 17/12/2020

    inputParameters["species_parameters_1"] = {
           "z"                  : 1.0       ,\
           "mass"               : 1.0       ,\
           "dens"               : 1.0       ,\
           "temp"               : 1.0       ,\
           "tprim"              : -999.9    ,\
           "fprim"              : -999.9    ,\
           "d2ndr2"             : 0.0       ,       # This is (1/n_s)*d^2 n_s / drho^2
           "d2Tdr2"             : 0.0       ,       # This is (1/T_s)*d^2 T_s / drho^2
           "type"               : "default" }

    inputParameters["species_parameters_2"] = inputParameters["species_parameters_1"].copy()
    inputParameters["species_parameters_3"] = inputParameters["species_parameters_1"].copy()
    inputParameters["species_parameters_4"] = inputParameters["species_parameters_1"].copy()
    inputParameters["species_parameters_5"] = inputParameters["species_parameters_1"].copy()
    inputParameters["species_parameters_6"] = inputParameters["species_parameters_1"].copy()
    inputParameters["species_parameters_7"] = inputParameters["species_parameters_1"].copy()
    inputParameters["species_parameters_8"] = inputParameters["species_parameters_1"].copy()
    inputParameters["species_parameters_9"] = inputParameters["species_parameters_1"].copy()
    inputParameters["species_parameters_10"] = inputParameters["species_parameters_1"].copy()

    inputParameters["species_parameters_a"] = { # Custom knob for the adiabtic electrons
           "z"                  : -1.0      ,\
           "mass"               : 5.43867E-4,\
           "dens"               : 1.0       ,\
           "temp"               : 1.0       ,\
           "tite"               : -999.9    ,\
           "nine"               : -999.9    ,\
           "tprim"              : -999.9    ,\
           "fprim"              : -999.9    ,\
           "type"               : "adiabatic" }

    inputParameters["stella_diagnostics_knobs"] = {
           "nwrite"             : 50        ,\
           "navg"               : 50        ,\
           "nmovie"             : 10000     ,\
           "nsave"              : -1        ,\
           "save_for_restart"   : False     ,\
           "write_omega"        : False     ,\
           "write_phi_vs_time"  : False     ,\
           "write_gvmus"        : False     ,\
           "write_gzvs"         : False     ,\
           "write_kspectra"     : False     ,\
           "write_moments"      : False     ,\
           "flux_norm"          : True      ,\
           "write_radial_fluxes": "radial_variation", # New parameter 17/12/2020
           "write_fluxes_kxky"  : False}        # Old parameter: has been removed! should be readded!

    inputParameters["millergeo_parameters"] = { #17
           "rhoc"               : None      ,\
           "rmaj"               : None      ,\
           "shift"              : None      ,\
           "qinp"               : None      ,\
           "shat"               : None      ,\
           "kappa"              : None      ,\
           "kapprim"            : None      ,\
           "tri"                : None      ,\
           "triprim"            : None      ,\
           "rgeo"               : None      ,\
           "betaprim"           : None      ,\
           "betadbprim"         : None      ,\
           "d2qdr2"             : None      ,\
           "d2psidr2"           : None      ,\
           "nzed_local"         : None      ,\
           "read_profile_variation" : None  ,\
           "write_profile_variation": None  }

    inputParameters["layouts_knobs"] = {
           "xyzs_layout"        : 'xyzs'    ,\
           "vms_layout"         : 'vms'     }

    inputParameters["neoclassical_input"] = {
           "include_neoclassical_terms" : False,    # Set to .true. to include neoclassical terms in GK equation
           "nradii"             : 5         ,       # Number of radial points used for radial derivatives of neoclassical quantities
           "drho"               : 0.01      ,       # Spacing in rhoc between points used for radial derivatives
           "neo_option"         : 'sfincs'  }       # Option for obtaining neoclassical distribution function and potential

    inputParameters["sfincs_input"] = { # 26
            "read_sfincs_output_from_file" : False,  # If True will try to read in Phi1Hat and delta_f from pre-saved file named sfincs.output
                                                     # otherwise, run sfincs to compute these quantities on the fly
            "nproc_sfincs"      :  1        ,        # Number of processors to use for sfincs calculation
            "irad_min"          : "-nradii/2",       # Minimum and maximum radial index (irad=0 corresponds to central radius)
            "irad_max"          : "nradii/2",        # Minimum and maximum radial index (irad=0 corresponds to central radius)
            "calculate_radial_electric_field" : True,# If True then scan in radial electric field to find value for which ambipolarity is satisfied, 
                                                     # and then use this value to obtain neoclassical fluxes, distribution function, and potential
            "includeXDotTerm"   : True      ,        # Do not include radial electric field term if set to False
            "includeElectricFieldTermInXiDot" : True,\
            "magneticDriftScheme": 0        ,        # No poloidal or toroidal magnetic drifts
            "includePhi1"       : True      ,        # Combo of next two variables means phi1 will be calculated via quasineutrality
            "includePhi1InKineticEquation" : False ,\
            "geometryScheme"    : 1         ,        # Will be overridden by geometric quantities unless geometryScheme = 5 (vmec equilibrium)
            "VMECRadialOption"  : 0         ,        # Only relevant if geometryScheme = 5 radial option to use for vmec equilibrium 0 corresponds 
                                                     # to using radial interpolation to get desired surface 1 corresponds to using nearest surface on
                                                     # VMEC HALF grid 2 corresponds to using nearest surface on VMEC FULL grid should not change this
                                                     # unless self-consistently change in the vmec input namelist 
            "equilibriumFile"   : 'wout_161s1.nc',   # Path of vmec equilibrium file
            "coordinateSystem"  : 3         ,        # Seems to be a nonsensical option
            "inputRadialCoordinate" : 3     ,        # Option 3: using sqrt of toroidal flux normalized by toroidal flux enclosed by the LCFS
            "inputRadialCoordinateForGradients" :3,  # Option 3: same choice when calculating gradients of density, temperature, and potential
            "aHat"              : 1.0       ,        # Dorresponds to r_LCFS as reference length in sfincs  only used in sfincs when geometryScheme=1
            "psiAHat"           : "geo_surf%psitor_lcfs",   # psitor_LCFS / (B_ref * a_ref^2)
            "Delta"             : -1.0      ,        # Delta is rho* = mref*vt_ref/(e*Bref*aref), with reference quantities given in SI units unless
                                                     # geometryScheme = 5, in which case Bref=1T and aref = 1m (these are hardwired in sfincs) set
                                                     # negative to allow check later to see if any value given in input file
            "nu_n"              : -1.0      ,        # nu_n = nu_ref * aref/vt_ref; nu_ref = 4*sqrt(2*pi)*nref*e**4*loglam/(3*sqrt(mref)*Tref**3/2);
                                                     # (with nref, Tref, and mref in Gaussian units) set negative to check later for input file
            "dPhiHatdrN"        : -9999.9   ,        # Radial derivative of normalized phi
            "Er_window"         : 0.3       ,\
            "nxi"               : 48        ,        # Number of spectral coefficients in pitch angle
            "nx"                : 12        ,        # Number of speeds
            "Ntheta"            : 65        ,        # Number of poloidal angles
            "Nzeta"             : 1         }        # Number of toroidal angles, 1 is appropriate for tokamak

    inputParameters["multibox_parameters"] : { # New namelist 17/12/2020
            "boundary_size"     : 4         ,\
            "krook_size"        : 0         ,\
            "phi_bound"         : 0         ,\
            "phi_pow"           : 2         ,\
            "krook_exponent"    : 0.0       ,\
            "nu_krook_mb"       : 0.0       ,\
            "mb_debug_step"     : 1000      ,\
            "smooth_ZFs"        : False     ,\
            "comm_at_init"      : False     ,\
            "RK_step"           : False     ,\
            "zf_option"         : 'default' ,\
            "krook_option"      : 'default' ,\
            "LR_debug_option"   : 'default'}

    inputParameters["dissipation"] : { # New namelist 17/12/2020
           "include_collisions" : False     ,\
           "include_krook_operator" : False ,\
           "collisions_implicit": True      ,\
           "collision_model"    : "dougherty" ,\
           "momentum_conservation" : True   ,\
           "energy_conservation": True      ,\
           "vpa_operator"       : True      ,\
           "mu_operator"        : True      ,\
           "hyper_dissipation"  : False     ,\
           "remove_zero_projection" : False ,\
           "D_hyper"            : 0.05      ,\
           "nu_krook"           : 0.05      ,\
           "delay_krook"        : 0.02      ,\
           "ikxmax_source"      : 2         ,       # kx:0 and kx:1
           "krook_odd"          : True      ,       # damp only the odd mode that can affect profiles
           "cfac"               : 1}
        
    
    inputParameters[" "] = {
            " "                 : " "       }        # To display empty labels on the GUI

    return inputParameters
