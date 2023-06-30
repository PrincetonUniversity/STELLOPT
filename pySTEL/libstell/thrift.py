def read_thrift(file):
    import h5py
    import numpy as np
    thrift_data = {}
    #read file
    with h5py.File(file,'r') as f:
        # Logicals
        for temp in ['leccd','lohmic','lvmec','lnbcd']:
            if temp in f:
                thrift_data[temp] = np.int64(f[temp][0])
        # Integers
        for temp in ['npicard','nrho','nssize','ntimesteps']:
            if temp in f:
                thrift_data[temp] = np.int64(f[temp][0])
        # Floats
        for temp in ['VERSION','jtol','picard_factor']:
            if temp in f:
                thrift_data[temp] = np.float64(f[temp][0])
        # Arrays
        for temp in ['THRIFT_ALPHA1','THRIFT_ALPHA2','THRIFT_ALPHA3','THRIFT_ALPHA4','THRIFT_AMINOR',\
        'THRIFT_BAV','THRIFT_BSQAV','THRIFT_BVAV','THRIFT_COEFF_A','THRIFT_COEFF_B','THRIFT_COEFF_BP',\
        'THRIFT_COEFF_C','THRIFT_COEFF_CP','THRIFT_COEFF_D','THRIFT_COEFF_DP','THRIFT_ETAPARA','THRIFT_I',\
        'THRIFT_IBOOT','THRIFT_IECCD','THRIFT_INBCD','THRIFT_IOHMIC','THRIFT_IOTA','THRIFT_IPLASMA','THRIFT_ISOURCE',\
        'THRIFT_J','THRIFT_JBOOT','THRIFT_JECCD','THRIFT_JNBCD','THRIFT_JOHMIC','THRIFT_JPLASMA','THRIFT_JSOURCE',\
        'THRIFT_MATLD','THRIFT_MATMD','THRIFT_MATRHS','THRIFT_MATUD','THRIFT_P','THRIFT_PHIEDGE','THRIFT_PPRIME',\
        'THRIFT_RHO','THRIFT_RHOFULL','THRIFT_RMAJOR','THRIFT_S','THRIFT_S11','THRIFT_S12','THRIFT_SNOB','THRIFT_T','\
        THRIFT_UGRID','THRIFT_VP']:
            if temp in f:
                thrift_data[temp] = np.array(f[temp][:])
    return thrift_data