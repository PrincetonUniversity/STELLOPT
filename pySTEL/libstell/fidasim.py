def read_fidasim(file):
    import h5py
    import numpy as np
    fidasim_data = {}
    #read file
    with h5py.File(file,'r') as f:
        # Integers
        for temp in ['nenergy', 'npitch', 'nr', 'nz', 'nphi','type']:
            if temp in f:
                fidasim_data[temp] = np.int64(f[temp][0])
        # Floats
        for temp in ['time']:
            if temp in f:
                fidasim_data[temp] = np.float64(f[temp][0])
        # Arrays
        for temp in ['r','phi','z','denf','f','energy','pitch']:
            if temp in f:
                fidasim_data[temp] = np.array(f[temp][:])

    return fidasim_data