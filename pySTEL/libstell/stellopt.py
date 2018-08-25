def read_stellopt_namelist(iunit,istat):
    import os, sys
    import numpy.ctypeslib as npct
    import numpy as np
    # Load Libraries
    # Setup Arrays
    stellopt_namelist={}
    # Will delete this just using for testing
    stellopt_namelist['NFUNC_MAX'] = 5000
    stellopt_namelist['OPT_TYPE'] = 'LMDIF'
    stellopt_namelist['EQUIL_TYPE'] = 'VMEC2000'
    stellopt_namelist['FTOL'] = 1.0E-6
    stellopt_namelist['XTOL'] = 1.0E-6
    stellopt_namelist['GTOL'] = 1.0E-30
    stellopt_namelist['EPSFCN'] = 1.0E-4
    stellopt_namelist['MODE'] = 1
    stellopt_namelist['FACTOR'] = 100
    stellopt_namelist['CR_STRATEGY'] = 0
    stellopt_namelist['NPOPULATION'] = -1
    stellopt_namelist['NOPTIMIZERS'] = -1
    # Vars
    for name in ['PHIEDGE','PRES_SCALE','CURTOR']:
        stellopt_namelist['L'+name+'_OPT'] = 0
        stellopt_namelist['D'+name+'_OPT'] = 1.0
        stellopt_namelist[    name+'_MIN'] = -1E30
        stellopt_namelist[    name+'_MAX'] = +1E30
    # Arrays
    for name in ['AM','AC','AI','EXTCUR']:
        stellopt_namelist['L'+name+'_OPT'] = np.ndarray(20)
        stellopt_namelist['D'+name+'_OPT'] = np.ndarray(20)
        stellopt_namelist[    name+'_MIN'] = np.ndarray(20)
        stellopt_namelist[    name+'_MAX'] = np.ndarray(20)
        for i in range(20):
            stellopt_namelist['L'+name+'_OPT'][i] = 0
            stellopt_namelist['D'+name+'_OPT'][i] = 1.0
            stellopt_namelist[    name+'_MIN'][i] = -1E10
            stellopt_namelist[    name+'_MAX'][i] = +1E10
    # Matrices
    arr_size1=2*101+1
    arr_size2=100+1
    stellopt_namelist['LBOUND_OPT'] = np.ndarray((arr_size1,arr_size2))
    stellopt_namelist['DBOUND_OPT'] = np.ndarray((arr_size1,arr_size2))
    stellopt_namelist['BOUND_MIN'] = np.ndarray((arr_size1,arr_size2))
    stellopt_namelist['BOUND_MAX'] = np.ndarray((arr_size1,arr_size2))
    stellopt_namelist['LRHO_OPT'] = np.ndarray((arr_size1,arr_size2))
    stellopt_namelist['DRHO_OPT'] = np.ndarray((arr_size1,arr_size2))
    stellopt_namelist['LDELTAMN_OPT'] = np.ndarray((arr_size1,arr_size1))
    stellopt_namelist['DDELTAMN_OPT'] = np.ndarray((arr_size1,arr_size1))
    stellopt_namelist['DELTA_MIN'] = np.ndarray((arr_size1,arr_size1))
    stellopt_namelist['DELTA_MAX'] = np.ndarray((arr_size1,arr_size1))
    for n in range(arr_size1):
        for m in range(arr_size2):
            stellopt_namelist['LBOUND_OPT'][n,m]=0
            stellopt_namelist['DBOUND_OPT'][n,m]=1.0
            stellopt_namelist['BOUND_MIN'][n,m]=-1.0E10
            stellopt_namelist['BOUND_MAX'][n,m]=+1.0E10
            stellopt_namelist['LRHO_OPT'][n,m]=0
            stellopt_namelist['DRHO_OPT'][n,m]=1.0
    for n in range(arr_size1):
        for m in range(arr_size1):
            stellopt_namelist['LDELTAMN_OPT'][n,m]=0
            stellopt_namelist['DDELTAMN_OPT'][n,m]=1.0
            stellopt_namelist['DELTA_MIN'][n,m]=-1.0E10
            stellopt_namelist['DELTA_MAX'][n,m]=+1.0E10

    return stellopt_namelist;

def write_stellopt_namelist(filename,stellopt_namelist):
    import os, sys
    import ctypes as ct
    import numpy.ctypeslib as npct
    import numpy as np
    # Load Libraries
    try:
        libstell = ct.cdll.LoadLibrary(os.environ["STELLOPT_PATH"]+"/LIBSTELL/Release/libstell.so")
    except KeyError:
        print("Please set environment variable STELLOPT_PATH")
        sys.exit(1)
    # Handle interface
    #read_stellopt_input = getattr(libstell,'__vmec_input_MOD_read_indata_namelist')
    #SUBROUTINE read_stellopt_input (iunit, istat)
    #read_stellopt_input.argtypes = [ct.POINTER(ct.c_int),ct.POINTER(ct.c_int)]
    #read_stellopt_input.restype=None
    #iunit_temp = ct.c_int(iunit)
    #istat_temp = ct.c_int(istat)
    #read_stellopt_input(ct.byref(iunit_temp),ct.byref(istat_temp))
    #istat = istat_temp
    #iunit = iunit_temp
    # Setup Arrays
    stellopt_namelist={}
    # Will delete this just using for testing
    for i,item in enumerate(stellopt_namelist):
        print(item)
    return stellopt_namelist;


