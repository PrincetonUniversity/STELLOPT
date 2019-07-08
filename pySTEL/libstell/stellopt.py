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

def read_stellopt(filename):
    import numpy as np
    #    import numpy as np
    file_handle = open(filename,'r')
    stel_data={}
    niter = 0
    for line in file_handle:
        if 'ITER' in line:
            niter=niter+1
    stel_data['ITER'] = np.ndarray((niter,1));
    file_handle.seek(0)
    line = file_handle.readline()
    ttype,wh=line.split()
    stel_data[ttype] = float(wh)

    # Enter Loop
    citer = -1
    while True:
        line = file_handle.readline()
        if line == '':
            break
        ttype,hw = line.split(' ',1)
        if ttype == 'ITER':
            citer = citer+1
            stel_data[ttype][citer] = int(hw)
            continue
        elif ttype == 'VERSION':
            stel_data[ttype][citer] = float(hw)
            continue
        elif ttype == 'TARGETS':
            h,w = hw.split()
            h = int(h)
            w = int(w)
            if h == 1:
                h=w
                w=1
            line = file_handle.readline()
        elif ttype == 'SIGMAS':
            h,w = hw.split()
            h = int(h)
            w = int(w)
            if h == 1:
                h=w
                w=1
            line = file_handle.readline()
        elif ttype == 'VALS':
            h,w = hw.split()
            h = int(h)
            w = int(w)
            if h == 1:
                h=w
                w=1
            line = file_handle.readline()
        else:
            h,w = hw.split()
            h = int(h)
            w = int(w)
            line = file_handle.readline()
        if ttype not in stel_data:
            stel_data[ttype]=np.ndarray((niter,h,w))
        for i in range(h):
            line = file_handle.readline()
            val = np.fromstring(line,sep=' ')
            stel_data[ttype][citer,i,:] = val       
    file_handle.close()
    for item in list(stel_data):
        #print(item)
        if 'VERSION' == item:
            continue
        elif 'ITER' == item:
            continue
        elif item in ['ASPECT','ASPECT_MAX','BETA','CURTOR','KAPPA','PHIEDGE', \
                    'VOLUME','WP','RBTOR','R0','Z0','BETATOR','BETAPOL']:
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
        elif item == 'BALLOON':
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_grate'] = np.squeeze(stel_data[item][:,:,3])
            stel_data[item+'_theta'] = np.squeeze(stel_data[item][:,:,4])
            stel_data[item+'_zeta'] = np.squeeze(stel_data[item][:,:,5])
            stel_data[item+'_k'] = np.squeeze(stel_data[item][:,:,6])
        elif item == 'B_PROBES':
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,4])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,5])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,6])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_X'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_Y'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_Z'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_MODB'] = np.squeeze(stel_data[item][:,:,3])
        elif item in ['FLUXLOOPS','SEGROG']:
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
        elif item == 'EXTCUR':
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_dex'] = np.squeeze(stel_data[item][:,:,3])
        elif item in ['SEPARATRIX','LIMITER']:
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_R'] = np.squeeze(stel_data[item][:,:,3])
            stel_data[item+'_PHI'] = np.squeeze(stel_data[item][:,:,4])
            stel_data[item+'_Z'] = np.squeeze(stel_data[item][:,:,5])
        elif item in ['NE','TI','TE','IOTA','VPHI','PRESS','VACIOTA']:
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,4])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,5])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,6])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_R'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_PHI'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_Z'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_s'] = np.squeeze(stel_data[item][:,:,3])
        elif item in ['NELINE','TELINE','TILINE','FARADAY','SXR','XICS','XICS_BRIGHT','XICS_W3','XICS_V']:
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_R0'] = np.squeeze(stel_data[item][:,:,3])
            stel_data[item+'_PHI0'] = np.squeeze(stel_data[item][:,:,4])
            stel_data[item+'_Z0'] = np.squeeze(stel_data[item][:,:,5])
            stel_data[item+'_R1'] = np.squeeze(stel_data[item][:,:,6])
            stel_data[item+'_PHI1'] = np.squeeze(stel_data[item][:,:,7])
            stel_data[item+'_Z1'] = np.squeeze(stel_data[item][:,:,8])
        elif item == 'MSE':
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,4])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,5])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,8])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_R'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_PHI'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_Z'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_s'] = np.squeeze(stel_data[item][:,:,3])
            stel_data[item+'_ER'] = np.squeeze(stel_data[item][:,:,6])
            stel_data[item+'_EZ'] = np.squeeze(stel_data[item][:,:,7])
        elif item == 'BOOTSTRAP':
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_s'] = np.squeeze(stel_data[item][:,:,3])
            stel_data[item+'_avg_jdotb'] = np.squeeze(stel_data[item][:,:,4])
            stel_data[item+'_beam_jdotb'] = np.squeeze(stel_data[item][:,:,5])
            stel_data[item+'_boot_jdotb'] = np.squeeze(stel_data[item][:,:,6])
            stel_data[item+'_jBbs'] = np.squeeze(stel_data[item][:,:,7])
            stel_data[item+'_facnu'] = np.squeeze(stel_data[item][:,:,8])
            stel_data[item+'_bsnorm'] = np.squeeze(stel_data[item][:,:,9])
        elif item == 'HELICITY':
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_bnorm'] = np.squeeze(stel_data[item][:,:,3])
        elif item == 'HELICITY_FULL':
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_bnorm'] = np.squeeze(stel_data[item][:,:,3])
            stel_data[item+'_k'] = np.squeeze(stel_data[item][:,:,4])
            stel_data[item+'_m'] = np.squeeze(stel_data[item][:,:,5])
            stel_data[item+'_n'] = np.squeeze(stel_data[item][:,:,6])
        elif item == 'TXPORT':
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_s'] = np.squeeze(stel_data[item][:,:,3])
        elif item == 'KINK':
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_wp'] = np.squeeze(stel_data[item][:,:,3])
            stel_data[item+'_wk'] = np.squeeze(stel_data[item][:,:,4])
            stel_data[item+'_omega'] = np.squeeze(stel_data[item][:,:,5])
        elif item == 'COIL_BNORM':
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_U'] = np.squeeze(stel_data[item][:,:,3])
            stel_data[item+'_V'] = np.squeeze(stel_data[item][:,:,4])
            stel_data[item+'_BNEQ'] = np.squeeze(stel_data[item][:,:,5])
            stel_data[item+'_BNF'] = np.squeeze(stel_data[item][:,:,6])
        elif item == 'ORBIT':
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_s'] = np.squeeze(stel_data[item][:,:,3])
        elif item == 'J_STAR':
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_AVGJSTAR'] = np.squeeze(stel_data[item][:,:,3])
            stel_data[item+'_TRAPSJSTAR'] = np.squeeze(stel_data[item][:,:,4])
            stel_data[item+'_UJSTAR'] = np.squeeze(stel_data[item][:,:,5])
            stel_data[item+'_K'] = np.squeeze(stel_data[item][:,:,6])
            stel_data[item+'_IJSTAR'] = np.squeeze(stel_data[item][:,:,7])
        elif item == 'NEO':
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_k'] = np.squeeze(stel_data[item][:,:,3])
        elif item == 'JDOTB':
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_s'] = np.squeeze(stel_data[item][:,:,3])
        elif item == 'JTOR':
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_s'] = np.squeeze(stel_data[item][:,:,3])
        elif item == 'DKES':
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_s'] = np.squeeze(stel_data[item][:,:,3])
            stel_data[item+'_nu'] = np.squeeze(stel_data[item][:,:,4])
            stel_data[item+'_er'] = np.squeeze(stel_data[item][:,:,5])
            stel_data[item+'_L11p'] = np.squeeze(stel_data[item][:,:,6])
            stel_data[item+'_L11m'] = np.squeeze(stel_data[item][:,:,7])
            stel_data[item+'_L33p'] = np.squeeze(stel_data[item][:,:,8])
            stel_data[item+'_L33m'] = np.squeeze(stel_data[item][:,:,9])
            stel_data[item+'_L31p'] = np.squeeze(stel_data[item][:,:,10])
            stel_data[item+'_L31m'] = np.squeeze(stel_data[item][:,:,11])
            stel_data[item+'_scal11'] = np.squeeze(stel_data[item][:,:,12])
            stel_data[item+'_scal33'] = np.squeeze(stel_data[item][:,:,13])
            stel_data[item+'_scal31'] = np.squeeze(stel_data[item][:,:,14])
        elif item == 'ECEREFLECT':
            stel_data[item+'_target'] = np.squeeze(stel_data[item][:,:,0])
            stel_data[item+'_sigma'] = np.squeeze(stel_data[item][:,:,1])
            stel_data[item+'_equil'] = np.squeeze(stel_data[item][:,:,2])
            stel_data[item+'_chisq'] = ((stel_data[item+'_target'] - stel_data[item+'_equil'])/stel_data[item+'_sigma'])**2
            stel_data[item+'_freq'] = np.squeeze(stel_data[item][:,:,3])
            stel_data[item+'_tradx'] = np.squeeze(stel_data[item][:,:,4])
            stel_data[item+'_trado'] = np.squeeze(stel_data[item][:,:,5])
            stel_data[item+'_mix'] = np.squeeze(stel_data[item][:,:,6])

    return stel_data;


