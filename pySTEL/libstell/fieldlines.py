def read_fieldlines(file):
    import h5py
    import numpy as np
    fieldline_data = {}
    #read file
    with h5py.File(file,'r') as f:
        # Logicals
        for temp in ['ladvanced', 'laxis_i', 'lcoil', 'lmgrid', 'lmu', 'lpies', 'lreverse', 'lspec', 'lvac', 'lvessel', 'lvmec']:
            fieldline_data[temp] = np.int(f[temp][0])
        # Integers
        for temp in ['nlines', 'nphi', 'npoinc', 'nr', 'nsteps', 'nz']:
            fieldline_data[temp] = np.int(f[temp][0])
        # Floats
        for temp in ['VERSION','iota0']:
            fieldline_data[temp] = np.float(f[temp][0])
        # 1D Arrays
        for temp in ['phiaxis', 'raxis', 'zaxis', 'B_lines', 'PHI_lines', 'R_lines', 'Z_lines', 'B_PHI', 'B_R', 'B_Z']:
            fieldline_data[temp] = np.array(f[temp][:])
    return fieldline_data