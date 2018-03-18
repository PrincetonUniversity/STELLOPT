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
    return stellopt_namelist;