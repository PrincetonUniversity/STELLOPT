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