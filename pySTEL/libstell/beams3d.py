def read_beams3d(file):
    import h5py
    import numpy as np
    beams_data = {}
    #read file
    with h5py.File(file,'r') as f:
        # Logicals
        for temp in ['lascot', 'lbeam', 'lbeam_simple', 'lcoil', 'lcollision', 'ldepo', 'lflux', 'lmgrid', 'lpies', 'lspec',\
                     'lvac', 'lvessel', 'lvmec']:
            if temp in f:
                beams_data[temp] = (f[temp][0]).astype(int)
        # Integers
        for temp in ['nbeams', 'nface', 'nparticles', 'nphi', 'npoinc', 'nz', 'nsteps', 'nvertex', 'nz']:
            if temp in f:
                beams_data[temp] = (f[temp][0]).astype(int)
        # Floats
        for temp in ['VERSION','iota0']:
            if temp in f:
                beams_data[temp] = (f[temp][0]).astype(float)
        # Arrays
        for temp in ['B_PHI', 'B_R', 'B_Z', 'B_lines', 'NE', 'PE_lines', 'PHI_lines', 'PI_lines', 'POT_ARR', 'R_lines',\
                    'S_ARR', 'S_lines', 'Shinethrough', 'TE', 'TI', 'U_ARR', 'U_lines', 'ZEFF_ARR', 'Z_lines', 'Zatom', \
                    'charge','mass','moment_lines','neut_lines','phiaxis','raxis','t_end','vll_lines','wall_faces','wall_strikes',\
                    'wall_vertex','zaxis','end_state','wall_load','wall_shine','epower_prof','ipower_prof','j_prof',\
                    'ndot_prof','dense_prof','dist_prof','Weight','Beam']:
            if temp in f:
                beams_data[temp] = np.array(f[temp][:])
    # Make derived arrays
    beams_data['X_lines'] = beams_data['R_lines']*np.cos(beams_data['PHI_lines'])
    beams_data['Y_lines'] = beams_data['R_lines']*np.sin(beams_data['PHI_lines'])
    return beams_data