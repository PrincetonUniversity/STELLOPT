# LIBSTELL Module

import os as _os
import _sys, _platform
import ctypes as _ct
import numpy.ctypeslib as npct
import numpy as _np
import matplotlib.pyplot as _plt

os = _platform.system()
is_64bits = _sys.maxsize > 2**32

if os=='Windows':
    libstell_path = _os.environ["STELLOPT_PATH"]+"/LIBSTELL/Release/libstell.dll"

    s1 = "__"
    s2 = "MOD"
    s3 = ""
elif os=='Linux' or os=='Darwin':
    libstell_path = _os.environ["STELLOPT_PATH"]+"/LIBSTELL/Release/libstell.so"

    s1 = ""
    s2 = "mp"
    s3 = "__"
# end if

if is_64bits:   # 64-bit architecture
    ls_handle = _ct.c_longlong
else:
    ls_handle = _ct.c_long
# end if


def LoadLibrary():
    try:
        libstell = _ct.cdll.LoadLibrary(libstell_path)
    except KeyError:
        print("Please set environment variable STELLOPT_PATH")
        _sys.exit(1)
    # end try
    return libstell

def vname(module, var):
    """Returns the DLL path to the given module and variable"""
    return '%s%s_%s_%s%s'%(s1, module, s2, var, s3)


def set_module_var(module,var,val):
    # Load Libraries
    libstell = LoadLibrary()

    if type(val) == bool:
        f = _ct.c_bool
    elif type(val) == int:
        f = _ct.c_int
    elif type(val) == float:
        f = _ct.c_double
    elif type(val) == str:
        n = len(val)
        f = _ct.c_char*n
    elif type(val) == _np.ndarray:
        if type(val[0]) == bool:
            tt = _ct.c_bool
        elif type(val[0]) == _np.int32:
            tt = _ct.c_int
        elif type(val[0]) == _np.float64:
            tt = _ct.c_double
        else:
            print('   Unrecognized type:',type(val[0]))
            return
        n = val.ndim
        f = tt*val.size
        #print(n,val.size)
    else:
        print('   Unrecognized type:',type(val))
        return
    # end if
    # temp=f.in_dll(libstell,s1+''+module+'_'+s2+'_'+var)
    temp=f.in_dll(libstell, vname(module, var))

    if type(val) == _np.ndarray:
        if n==1:
            for i,col in enumerate(val):
                temp[i] = val[i]
    elif type(val) == str:
        temp.value = val.encode('UTF-8')
    else:
        temp.value = val
    return


def read_vmec(file):
    # Load Libraries
    libstell = LoadLibrary()
    # Read File
    read_wout = getattr(libstell,vname('read_wout_mod_', '_readw_and_open'))
    read_wout.argtypes=[_ct.c_char_p, _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_int), _ct.c_long]
    read_wout.restype=None
    ierr = _ct.c_int(0)
    iopen = _ct.c_int(0)
    read_wout(file.encode('UTF-8'), _ct.byref(ierr), _ct.byref(iopen), len(file))

    # Setup Arrays
    vmec_data={}
    # Check
    if not (ierr.value == 0):
        return vmec_data
    # Logical
    varlist=['lasym','lthreed','lwout_opened']
    for temp in varlist:
        vmec_data[temp]=_ct.c_bool.in_dll(libstell, vname('read_wout_mod_', temp)).value
    # end for

    # Integers
    varlist=['ns','nfp','mpol','ntor','mnmax','mnmax_nyq','iasym','ierr_vmec']
    for temp in varlist:
        vmec_data[temp]=_ct.c_int.in_dll(libstell, vname('read_wout_mod_', temp)).value
    # end for

    # Doubles
    varlist=['wb','wp','gamma','pfac','rmax_surf','rmin_surf','zmax_surf',\
             'aspect','betatot','betapol','betator','betaxis','b0','version_',\
             'ionlarmor','volavgb','fsql','fsqr','fsqz','ftolv','aminor','rmajor',\
             'volume','rbtor','rbtor0','itor','machsq']
    for temp in varlist:
        vmec_data[temp]=_ct.c_double.in_dll(libstell, vname('read_wout_mod_', temp)).value
    # end for

    # REAL Arrays (ns)
    varlist = ['iotas','iotaf','presf','phipf','chipf','chi','phi','mass',\
               'pres','beta_vol','phip','buco','bvco','vp','overr','jcuru',\
               'jcurv','specw','jdotb','dmerc','dwell','dcurr','dgeod','equif']
    arr_size = vmec_data['ns']
    ftemp = _ct.POINTER(_ct.c_double)
    for temp in varlist:
        vmec_data[temp] = npct.as_array(ftemp.in_dll(libstell, vname('read_wout_mod_', temp)),(arr_size,1))
    # end for

    # REAL Arrays (mnmax)
    ftemp = _ct.POINTER(_ct.c_double)
    vmec_data['xm'] = npct.as_array(ftemp.in_dll(libstell, vname('read_wout_mod_', '_xm')),(vmec_data['mnmax'],1))
    vmec_data['xn'] = npct.as_array(ftemp.in_dll(libstell, vname('read_wout_mod_', '_xn')),(vmec_data['mnmax'],1))
    vmec_data['xm_nyq'] = npct.as_array(ftemp.in_dll(libstell, vname('read_wout_mod_', '_xm_nyq')),(vmec_data['mnmax_nyq'],1))
    vmec_data['xn_nyq'] = npct.as_array(ftemp.in_dll(libstell, vname('read_wout_mod_', '_xn_nyq')),(vmec_data['mnmax_nyq'],1))

    ## Array values 1D
    ftemp=_ct.POINTER(_ct.c_double)
    ns = vmec_data['ns']
    mnmax = vmec_data['mnmax']
    mnmax_nyq = vmec_data['mnmax_nyq']
    # ns_size = (ns,1)
    # mn_size = (mnmax,1)
    # mnnyq_size = (mnmax_nyq,1)

    ## 2D Arrays
    mn2d_size = (ns, mnmax)
    mn2d_nyq_size = (ns, mnmax_nyq)
    fmn=_ct.POINTER(_ct.c_double)

    # 2D spatial harmonics of size ns, mnmax format
    keys = ['_rmnc', '_zmns', '_lmns']
    if vmec_data['iasym']:
        keys += ['_rmns', '_zmnc', '_lmnc']
    # end if
    for key in keys:
        vmec_data[key.strip('_')] = npct.as_array(fmn.in_dll(libstell, vname('read_wout_mod_', key)), mn2d_size) #ns,mnmax format
    # end for

    # 2D magnetic field, g-factor, current (sub/sup) harmonics of size mn2d_nyq_size
    keys = ['_bmnc', '_gmnc', '_bsupumnc', '_bsupvmnc', '_bsubsmns', '_bsubumnc', '_bsubvmnc', '_currumnc', '_currvmnc']

    if vmec_data['iasym']:
        keys += ['_bmns', '_gmns', '_bsupumns', '_bsupvmns', '_bsubsmnc', '_bsubumns', '_bsubvmns', '_currumns', '_currvmns']
    # end if
    for key in keys:
        vmec_data[key.strip('_')] = npct.as_array(fmn.in_dll(libstell, vname('read_wout_mod_', key)), mn2d_nyq_size) #ns,mnmax format
    # end for

    # Free memory (don't do this as python accesses this memory)
    #read_wout_dealloc = getattr(libstell,s1+'read_wout_mod_'+s2+'_read_wout_deallocate'+s3)
    #read_wout_dealloc()

    # Correct Arrays (mn-nv) to (mn+nv)
    vmec_data['xn'] = -vmec_data['xn']

    # Put on full grid
    vmec_data['buco'] = h2f(vmec_data['buco'],ns)
    vmec_data['bvco'] = h2f(vmec_data['bvco'],ns)
    vmec_data['vp'] = h2f(vmec_data['vp'],ns)
    vmec_data['overr'] = h2f(vmec_data['overr'],ns)
    vmec_data['specw'] = h2f(vmec_data['specw'],ns)
    # Put matrix quantities on full grid
    for key in ['bmnc','gmnc','lmns','bsupumnc','bsupvmnc','bsubsmns','bsubumnc','bsubvmnc']:
        vmec_data[key][0,:] = 1.5 * vmec_data[key][1,:] - 0.5 * vmec_data[key][2,:]
        vmec_data[key][1:ns-2,:] = 0.5 * (vmec_data[key][1:ns-2,:] + vmec_data[key][2:ns-1,:])
        vmec_data[key][ns-1,:] = 2.0 * vmec_data[key][ns-2,:] - vmec_data[key][ns-3,:]
    if vmec_data['iasym']:
        for key in ['bmns','gmns','lmnc','bsupumns','bsupvmns','bsubsmnc','bsubumns','bsubvmns']:
            vmec_data[key][0,:] = 1.5 * vmec_data[key][1,:] - 0.5 * vmec_data[key][2,:]
            vmec_data[key][1:ns-2,:] = 0.5 * (vmec_data[key][1:ns-2,:] + vmec_data[key][2:ns-1,:])
            vmec_data[key][ns-1,:] = 2.0 * vmec_data[key][ns-2,:] - vmec_data[key][ns-3,:]
    # Return the data dictionary

    return vmec_data


def h2f(var,ns):
    temp = _np.zeros((ns,1))
    temp[0] = 1.5 * var[0] - 0.5 * var[1]
    temp[1:ns-1] = 0.5* (var[0:ns-2] + var[1:ns-1])
    temp[ns-1] = 1.5 * var[ns-2] - 0.5 * var[ns-3]
    return temp


def cfunct(theta,zeta,fmnc,xm,xn):
    f=0
    (ns,mn)=fmnc.shape
    lt = len(theta)
    lz = len(zeta)
    mt =_np.matmul(xm,theta.T)
    nz =_np.matmul(xn,zeta.T)
    cosmt = _np.cos(mt)
    sinmt = _np.sin(mt)
    cosnz = _np.cos(nz)
    sinnz = _np.sin(nz)
    f = _np.zeros((ns,lt,lz))
    fmn = _np.ndarray((mn,lt))
    for k in range(ns):
        fmn = _np.broadcast_to(fmnc[k,:],(lt,mn)).T
        fmncosmt =(fmn*cosmt).T
        fmnsinmt =(fmn*sinmt).T
        f[k,:,:] = _np.matmul(fmncosmt, cosnz)-_np.matmul(fmnsinmt, sinnz)
    return f


def sfunct(theta,zeta,fmnc,xm,xn):
    f=0
    (ns,mn)=fmnc.shape
    lt = len(theta)
    lz = len(zeta)
    mt = _np.matmul(xm,theta.T)
    nz = _np.matmul(xn,zeta.T)
    cosmt = _np.cos(mt)
    sinmt = _np.sin(mt)
    cosnz = _np.cos(nz)
    sinnz = _np.sin(nz)
    f =  _np.zeros((ns,lt,lz))
    fmn =  _np.ndarray((mn,lt))
    for k in range(ns):
        fmn = _np.broadcast_to(fmnc[k,:],(lt,mn)).T
        f[k,:,:] = _np.matmul((fmn*sinmt).T,cosnz) + _np.matmul((fmn*cosmt).T,sinnz)
    return f


def torocont(r,z,val,s):
    h = _plt.axes(xlabel='R [m]',ylabel='Z [m]',aspect='equal')
    _plt.pcolormesh(r[:,:,s],z[:,:,s],val[:,:,s],cmap='jet',shading='gouraud',axes=h)
    _plt.show()
    return h


def toroslice(r,zeta,z,s):
    h = _plt.axes(xlabel='R [m]',ylabel='Z [m]',aspect='equal')
    if (s[0] == 0):
        _plt.plot(r[0,0,zeta],z[0,0,zeta],'+',color='black',axes=h)
        _plt.plot( _np.transpose(r[s[1:],:,zeta]), _np.transpose(z[s[1:],:,zeta]),color='black',axes=h)
    else:
        for i in range(0,):
            _plt.plot( _np.transpose(r[s,:,zeta]), _np.transpose(z[s,:,zeta]),color='black',axes=h)
    _plt.show()
    return h


def isotoro(r,z,zeta,svals,*args,**kwargs):
    import mpl_toolkits.mplot3d as mplot3d # analysis:ignore
    import math as math
    import matplotlib.tri as mtri
    #from mayavi import mlab
    nr =  _np.size(svals)
    if (nr == 1):
        s= [svals]
        nr = 1
    else:
        s=svals
    nt = _np.size(r,1)
    nz = _np.size(r,2)
    vertex = _np.zeros((nt*nz,3,nr))
    for k in range(0,nr):
        ivertex = 0
        # ifaces = 0
        for j in range(0,nz):
            for i in range(0,nt):
                vertex[ivertex, 0, k] = r[s[k],i,j]*math.cos(zeta[j])
                vertex[ivertex, 1, k] = r[s[k],i,j]*math.sin(zeta[j])
                vertex[ivertex, 2, k] = z[s[k],i,j]
                ivertex = ivertex + 1
            # end for
        # end for
    # end for
    u =  _np.linspace(0, 1, endpoint=True, num=nt)
    v =  _np.linspace(0, 1, endpoint=True, num=nz)
    u, v =  _np.meshgrid(u, v)
    u, v = u.flatten(), v.flatten()
    tri = mtri.Triangulation(u, v)
    test=len(kwargs)
    fig=kwargs.pop('fig',_plt.figure())
    h=kwargs.pop('axes',fig.add_subplot(111,projection='3d'))
    for k in range(0,nr):
        if (len(args)==0):
            tsurf = h.plot_trisurf(vertex[:,0,k], vertex[:,1,k], vertex[:,2,k], triangles=tri.triangles,color='red',shade='yes',linewidths=0.0)
            #tsurf=mlab.triangular_mesh(vertex[:,0,k], vertex[:,1,k], vertex[:,2,k], tri.triangless)
        else:
            # Matplotlib way (SLOW)
            vals = args[0][s[k],:,:].T.flatten()
            colors =  _np.mean(vals[tri.triangles], axis=1)
            tsurf = h.plot_trisurf(vertex[:,0,k], vertex[:,1,k], vertex[:,2,k], triangles=tri.triangles,cmap='jet',shade='yes',linewidths=0.0)
            tsurf.set_array(colors)
            tsurf.autoscale()
            #MAYAVI Way (need to figure out how to embed)
            #h    = mlab.figure()
            #vals = args[0][s[k],:,:].T.flatten()
            #tsurf=mlab.triangular_mesh(vertex[:,0,k], vertex[:,1,k], vertex[:,2,k], tri.triangles, scalars=vals, colormap='jet',figure=h)
            #print(type(tsurf))
        # end if
    # end for

    if (test==0):
        _plt.show()
    # end if
    return h


def calc_jll(vmec_data, theta, zeta ):
    # CALC_JLL(vmec_data,theta,zeta) Calculates the parallel current density.
    # This funciton takes a VMEC data structure (as read by read_vmec) and
    # theta/zeta arrays as input and outputs the parallel current density.

    # Example usage (Matlab)
    #      theta=0:2*pi/359:2*pi;
    #      zeta=0:2*pi/63:2*pi;
    #      data=read_vmec('wout.test');        % Reads VMEC wout file
    #      jll=calc_jll(vmec_data,theta,zeta); % Calculate the current
    #
    # Example usage (Python)
    #      theta= _np.linspace(0, 2* _np.pi, 360)
    #      zeta= _np.linspace(0, 2* _np.pi, 64)
    #      vmec_data=read_vmec('wout.nc')
    #      jll=calc_jll(vmec_data, theta, zeta)


    # Maintained by: Samuel Lazerson (lazerson@pppl.gov)
    # Version:       1.00

    b =cfunct(theta,zeta,vmec_data['bmnc'],    vmec_data['xm_nyq'],vmec_data['xn_nyq'])
    g =cfunct(theta,zeta,vmec_data['gmnc'],    vmec_data['xm_nyq'],vmec_data['xn_nyq'])
    bu=cfunct(theta,zeta,vmec_data['bsubumnc'],vmec_data['xm_nyq'],vmec_data['xn_nyq'])
    bv=cfunct(theta,zeta,vmec_data['bsubvmnc'],vmec_data['xm_nyq'],vmec_data['xn_nyq'])
    ju=cfunct(theta,zeta,vmec_data['currumnc'],vmec_data['xm_nyq'],vmec_data['xn_nyq'])
    jv=cfunct(theta,zeta,vmec_data['currvmnc'],vmec_data['xm_nyq'],vmec_data['xn_nyq'])

    if (vmec_data['iasym']):
        b += sfunct(theta,zeta,vmec_data['bmns'],    vmec_data['xm_nyq'],vmec_data['xn_nyq'])
        g += sfunct(theta,zeta,vmec_data['gmns'],    vmec_data['xm_nyq'],vmec_data['xn_nyq'])
        bu += sfunct(theta,zeta,vmec_data['bsubumns'],vmec_data['xm_nyq'],vmec_data['xn_nyq'])
        bv += sfunct(theta,zeta,vmec_data['bsubvmns'],vmec_data['xm_nyq'],vmec_data['xn_nyq'])
        ju += sfunct(theta,zeta,vmec_data['currumns'],vmec_data['xm_nyq'],vmec_data['xn_nyq'])
        jv += sfunct(theta,zeta,vmec_data['currvmns'],vmec_data['xm_nyq'],vmec_data['xn_nyq'])
    # end if

    return (bu*ju+bv*jv)/(g*b)



def safe_close(iunit):
    # Load Libraries
    libstell = LoadLibrary()

    # Handle interface
    safe_close_h = getattr(libstell, vname('safe_open_mod_', '_safe_close'))
    safe_close_h.restype=None
    iunit_temp = _ct.c_int(iunit)
    safe_close_h(_ct.byref(iunit_temp))
    return


def safe_open(iunit,istat,filename,filestat,fileform,record_in,access_in,delim_in):
    # Load Libraries
    libstell = LoadLibrary()

    # Handle interface
    safe_open_h = getattr(libstell, vname('safe_open_mod_', '_safe_open'))

    # SUBROUTINE safe_open(int iunit, int istat, char filename, char filestat, char fileform, int record_in, char access_in, char delim_in)
    safe_open_h.argtypes= [ _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_int), _ct.c_char_p, _ct.c_char_p, _ct.c_char_p, \
        _ct.POINTER(_ct.c_int), _ct.c_char_p, _ct.c_char_p, \
        _ct.c_long, _ct.c_long, _ct.c_long, _ct.c_long, _ct.c_long]
    safe_open_h.restype=None
    iunit_temp = _ct.c_int(iunit)
    istat_temp = _ct.c_int(istat)
    record_in_temp = _ct.c_int(record_in)
    # opt1 = _ct.c_bool(True)
    # opt2 = _ct.c_bool(True)
    # opt3 = _ct.c_bool(True)
    safe_open_h(_ct.byref(iunit_temp),_ct.byref(istat_temp), \
        filename.encode('UTF-8'),filestat.encode('UTF-8'),fileform.encode('UTF-8'),\
        _ct.byref(record_in_temp),access_in.encode('UTF-8'),delim_in.encode('UTF-8'), \
        len(filename),len(filestat),len(fileform),len(access_in),len(delim_in))
    istat = istat_temp
    iunit = iunit_temp
    istat = istat_temp
    return istat


def read_indata_namelist(iunit,istat):
    # Load Libraries
    libstell = LoadLibrary()

    # Handle interface
    read_indata_namelist = getattr(libstell, vname('vmec_input_', '_read_indata_namelist'))

    #SUBROUTINE read_indata_namelist (iunit, istat)
    read_indata_namelist.argtypes = [_ct.POINTER(_ct.c_int),_ct.POINTER(_ct.c_int)]
    read_indata_namelist.restype=None
    iunit_temp = _ct.c_int(iunit)
    istat_temp = _ct.c_int(istat)
    read_indata_namelist(_ct.byref(iunit_temp),_ct.byref(istat_temp))
    istat = istat_temp
    iunit = iunit_temp

    # Setup Arrays
    indata_namelist={}

    # Logicals
    varlist=['lpofr','lmac','lfreeb','lrecon','loldout','ledge_dump','lasym','lforbal','lrfp',\
             'lmovie','lmove_axis','lwouttxt','ldiagno','lmoreiter','lfull3d1out','l_v3fit',\
             'lspectrum_dump','loptim','lgiveup','lbsubs','lgiveup']
    for temp in varlist:
        # indata_namelist[temp]=_ct.c_bool.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3).value
        indata_namelist[temp]=_ct.c_bool.in_dll(libstell, vname('vmec_input_', temp)).value
    # end for

    # Integers
    varlist=['nfp','ncurr','nsin','niter','nstep','nvacskip','mpol','ntor','ntheta','nzeta', \
             'mfilter_fbdy','nfilter_fbdy','max_main_iterations','omp_num_threads',\
             'imse','isnodes','itse','ipnodes','iopt_raxis','imatch_phiedge','nflxs']
    for temp in varlist:
        # indata_namelist[temp]=_ct.c_int.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3).value
        indata_namelist[temp]=_ct.c_int.in_dll(libstell, vname('vmec_input_', temp)).value
    # end for

    # Reals
    varlist=['time_slice','curtor','delt','ftol','tcon0','gamma','phiedge','phidiam',\
             'sigma_current','sigma_delphid','tensi','tensp','tensi2','fpolyi','presfac',\
             'mseangle_offset','pres_offset','mseangle_offsetm','spres_ped','bloat',\
             'pres_scale','prec2d_threshold','bcrit','fgiveup']
    for temp in varlist:
        # indata_namelist[temp]=_ct.c_double.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3).value
        indata_namelist[temp]=_ct.c_double.in_dll(libstell, vname('vmec_input_', temp)).value
    # end for

    # Get integers defined elsewhere (Hardcode for now, not sure how to get them)
    #indata_namelist['nbsetsp']=_ct.c_int.in_dll(libstell,s1+'vsvd0_'+s2+'_nbsetsp').value
    # Integer Arrays (100)
    varlist = ['ns_array','niter_array']
    arr_size=100
    ftemp = _ct.c_int*arr_size
    for temp in varlist:
        # indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3),(arr_size,1))
        indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell, vname('vmec_input_', temp)),(arr_size,1))
    # end for

    # Note that we skip some arrays related to recon stuff because we don't need them and we
    # need to figure out how to pull stuff from other modules see the above issue.
    # Real 2D Arrays (ntord=101,mpol1d=100)
    varlist = ['rbs','zbc','rbc','zbs']
    arr_size1=2*101+1
    arr_size2=100+1
    ftemp = _ct.c_double*arr_size1*arr_size2
    for temp in varlist:
        # indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3),(arr_size1,arr_size2))
        indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell, vname('vmec_input_', temp)),(arr_size1,arr_size2))
    # end for

    # REAL Arrays (21)
    varlist = ['am','ai','ac','ah','at']
    arr_size=21
    ftemp = _ct.c_double*arr_size
    for temp in varlist:
        # indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3),(arr_size,1))
        indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell, vname('vmec_input_', temp)),(arr_size,1))
    # end for

    # REAL Arrays (20)
    varlist = ['aphi']
    arr_size=20
    ftemp = _ct.c_double*arr_size
    for temp in varlist:
        # indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3),(arr_size,1))
        indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell, vname('vmec_input_', temp)),(arr_size,1))
    # end for

    # REAL Arrays (ndatafmax=101)
    varlist = ['am_aux_s','am_aux_f','ac_aux_s','ac_aux_f','ai_aux_s','ai_aux_f',\
               'ah_aux_s','ah_aux_f','at_aux_s','at_aux_f']
    arr_size=101
    ftemp = _ct.c_double*arr_size
    for temp in varlist:
        # indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3),(arr_size,1))
        indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell, vname('vmec_input_', temp)),(arr_size,1))
    # end for

    # REAL Arrays (ntord+1=102)
    varlist = ['raxis','zaxis','raxis_cc','raxis_cs','zaxis_cc','zaxis_cs']
    arr_size=102
    ftemp = _ct.c_double*arr_size
    for temp in varlist:
        # indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3),(arr_size,1))
        indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell, vname('vmec_input_', temp)),(arr_size,1))
    # end for

    # REAL Arrays (100)
    varlist = ['ftol_array']
    arr_size=100
    ftemp = _ct.c_double*arr_size
    for temp in varlist:
        # indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3),(arr_size,1))
        indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell, vname('vmec_input_', temp)),(arr_size,1))
    # end for

    # REAL Arrays (nigroup=300)
    varlist = ['extcur']
    arr_size=300
    ftemp = _ct.c_double*arr_size
    for temp in varlist:
        # indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3),(arr_size,1))
        indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell, vname('vmec_input_', temp)),(arr_size,1))
    # end for

    # Charater arrays
    varlist = ['pcurr_type','piota_type','pmass_type','pt_type','ph_type']
    arr_size=20
    ftemp = _ct.c_char*arr_size
    for temp in varlist:
        # indata_namelist[temp]=ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3).value.decode('UTF-8')
        indata_namelist[temp]=ftemp.in_dll(libstell, vname('vmec_input_', temp)).value.decode('UTF-8')
    # end for

    ftemp = _ct.c_char*200
    indata_namelist['mgrid_file']=ftemp.in_dll(libstell, vname('vmec_input_', '_mgrid_file')).value.decode('UTF-8')
    # indata_namelist['mgrid_file']=ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_mgrid_file'+s3).value.decode('UTF-8')

    ftemp = _ct.c_char*200
    indata_namelist['trip3d_file']=ftemp.in_dll(libstell, vname('vmec_input_', '_trip3d_file')).value.decode('UTF-8')
    # indata_namelist['trip3d_file']=ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_trip3d_file'+s3).value.decode('UTF-8')

    ftemp = _ct.c_char*10
    indata_namelist['precon_type']=ftemp.in_dll(libstell, vname('vmec_input_', '_precon_type')).value.decode('UTF-8')
    # indata_namelist['precon_type']=ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_precon_type'+s3).value.decode('UTF-8')

    ftemp = _ct.c_char*120
    # indata_namelist['arg1']=ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_arg1'+s3).value.decode('UTF-8')
    indata_namelist['arg1']=ftemp.in_dll(libstell, vname('vmec_input_', '_arg1')).value.decode('UTF-8')

    ftemp = _ct.c_char*100
    indata_namelist['input_extension']=ftemp.in_dll(libstell, vname('vmec_input_', '_input_extension')).value.decode('UTF-8')
    # indata_namelist['input_extension']=ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_input_extension'+s3).value.decode('UTF-8')

    return indata_namelist


def write_indata_namelist(iunit,istat):
    # Load Libraries
    libstell = LoadLibrary()

    # Handle interface
    write_indata_namelist = getattr(libstell,s1+'vmec_input_'+s2+'_write_indata_namelist'+s3)
    #SUBROUTINE read_indata_namelist (iunit, istat)
    write_indata_namelist.argtypes = [_ct.POINTER(_ct.c_int),_ct.POINTER(_ct.c_int)]
    write_indata_namelist.restype=None
    iunit_temp = _ct.c_int(iunit)
    istat_temp = _ct.c_int(istat)
    write_indata_namelist(_ct.byref(iunit_temp),_ct.byref(istat_temp))
    #istat = istat_temp
    #iunit = iunit_temp
    return

def pcurr(xx):
    # Load Libraries
    libstell = LoadLibrary()

    pcurr_func = getattr(libstell,'pcurr_')
    #SUBROUTINE pcurr (xx)
    pcurr_func.argtypes = [_ct.POINTER(_ct.c_double)]
    pcurr_func.restype=_ct.c_double
    xx_temp = _ct.c_double(xx)
    val = pcurr_func(_ct.byref(xx_temp))
    return val;

def pmass(xx):
    # Load Libraries
    libstell = LoadLibrary()

    pmass_func = getattr(libstell,'pmass_')
    #SUBROUTINE piota (xx)
    pmass_func.argtypes = [_ct.POINTER(_ct.c_double)]
    pmass_func.restype=_ct.c_double
    xx_temp = _ct.c_double(xx)
    val = pmass_func(_ct.byref(xx_temp))
    return val;

def piota(xx):
    # Load Libraries
    libstell = LoadLibrary()

    piota_func = getattr(libstell,'piota_')
    #SUBROUTINE piota (xx)
    piota_func.argtypes = [_ct.POINTER(_ct.c_double)]
    piota_func.restype=_ct.c_double
    xx_temp = _ct.c_double(xx)
    val = piota_func(_ct.byref(xx_temp))
    return val;


# # ========================================================================= #
# # ===============================  read_wout_mod ========================== #
# # ========================================================================= #

# # @[0-9]{4} DATA
# read_wout_mod_internals = [
#     'ac',
#     'ac_aux_f',
#     'ac_aux_s',
#     'ai',
#     'ai_aux_f',
#     'ai_aux_s',
#     'alfa',
#     'am',
#     'am_aux_f',
#     'am_aux_s',
#     'aminor',
#     'anglemse',
#     'aspect',
#     'b0',
#     'bbc',
#     'bcwgt',
#     'bdotb',
#     'bdotgradv',
#     'beta_vol',
#     'betapol',
#     'betator',
#     'betatot',
#     'betaxis',
#     'bmnc',
#     'bmns',
#     'bsubsmnc',
#     'bsubsmns',
#     'bsubumnc',
#     'bsubumnc_sur',
#     'bsubumns',
#     'bsubumns_sur',
#     'bsubvmnc',
#     'bsubvmnc_sur',
#     'bsubvmns',
#     'bsubvmns_sur',
#     'bsupumnc',
#     'bsupumnc_sur',
#     'bsupumns',
#     'bsupumns_sur',
#     'bsupvmnc',
#     'bsupvmnc_sur',
#     'bsupvmns',
#     'bsupvmns_sur',
#     'buco',
#     'bvco',
#     'chi',
#     'chipf',
#     'compute_currents',
#     'curmid',
#     'currumnc',
#     'currumns',
#     'currvmnc',
#     'currvmns',
#     'datastark',
#     'datathom',
#     'dcurr',
#     'delphid',
#     'dgeod',
#     'dmerc',
#     'dshear',
#     'dsiobt',
#     'dwell',
#     'equif',
#     'extcur',
#     'flmwgt',
#     'fsql',
#     'fsqr',
#     'fsqt',
#     'fsqz',
#     'ftolv',
#     'gamma',
#     'gmnc',
#     'gmns',
#     'hotdmnc',
#     'hotdmns',
#     'iasym',
#     'ierr_vmec',
#     'imatch_phiedge',
#     'imse',
#     'input_extension',
#     'ionlarmor',
#     'iotaf',
#     'iotas',
#     'ipnodes',
#     'ireconstruct',
#     'isigng',
#     'isnodes',
#     'itfsq',
#     'itor',
#     'itse',
#     'jcuru',
#     'jcurv',
#     'jdotb',
#     'lasym',
#     'lmnc',
#     'lmns',
#     'loadrzl',
#     'lthreed',
#     'lwout_opened',
#     'machsq',
#     'mass',
#     'mgrid_file',
#     'mnmax',
#     'mnmax_nyq',
#     'mnmaxpot',
#     'mnyq',
#     'mpol',
#     'msewgt',
#     'nfp',
#     'niter',
#     'nnyq',
#     'ns',
#     'nstore_seq',
#     'ntmax',
#     'ntor',
#     'omega',
#     'overr',
#     'pbprmnc',
#     'pbprmns',
#     'pcurr_type',
#     'pfac',
#     'phi',
#     'phidiam',
#     'phip',
#     'phipf',
#     'piota_type',
#     'pknots',
#     'pmap',
#     'pmass_type',
#     'potcos',
#     'potsin',
#     'pparmnc',
#     'pparmns',
#     'ppermnc',
#     'ppermns',
#     'ppprmnc',
#     'ppprmns',
#     'pres',
#     'presf',
#     'presmid',
#     'protmnc',
#     'protmns',
#     'protrsqmnc',
#     'protrsqmns',
#     'prprmnc',
#     'prprmns',
#     'qfact',
#     'qmeas',
#     'qmid',
#     'raxis',
#     'rbtor',
#     'rbtor0',
#     'read_wout_deallocate',
#     'readw_and_open',
#     'readw_only',
#     'rmajor',
#     'rmax_surf',
#     'rmid',
#     'rmin_surf',
#     'rmnc',
#     'rmns',
#     'rstark',
#     'rthom',
#     'rzl_local',
#     'shear',
#     'sigmnc',
#     'sigmns',
#     'sknots',
#     'specw',
#     'taumnc',
#     'taumns',
#     'tosuvspace',
#     'tpotb',
#     'tswgt',
#     'version_',
#     'vmec_type',
#     'volavgb',
#     'volume',
#     'vp',
#     'wb',
#     'wdot',
#     'wp',
#     'write_wout_file',
#     'xm',
#     'xm_nyq',
#     'xmpot',
#     'xn',
#     'xn_nyq',
#     'xnpot',
#     'y2stark',
#     'y2thom',
#     'ystark',
#     'ythom',
#     'zaxis',
#     'zmax_surf',
#     'zmnc',
#     'zmns',
#     ]

# ========================================================================= #
# ============================== Base Class =============================== #
# ========================================================================= #


class __LibstellBase(object):

    def __init__(self):
        self.load_mod()
    # end def


    def load_mod(self):
        # Load Libraries
        libstell = LoadLibrary()
        self.mconf = libstell
    # end def load_mod


    def unload_mod(self):
        pass
    # end def unload_mod


    def vname(self, attname):
        return vname(self.module, attname)


    # # Integer Arrays (100)
    # varlist = ['ns_array','niter_array']
    # arr_size=100
    #
    # # REAL Arrays (nigroup=300)
    # varlist = ['extcur']
    # arr_size=300
    # ftemp = _ct.c_double*arr_size
    #
    # # Real 2D Arrays (ntord=101,mpol1d=100)
    # varlist = ['rbs','zbc','rbc','zbs']
    # arr_size1=2*101+1
    # arr_size2=100+1
    # ftemp = _ct.c_double*arr_size1*arr_size2
    #
    # # Charater arrays
    # varlist = ['pcurr_type','piota_type','pmass_type','pt_type','ph_type']
    # arr_size=20
    # ftemp = _ct.c_char*arr_size
    #

    # ==================================================== #
    # Getter, setter, and updating methods
    # ==================================================== #


    @property
    def module(self):
        return self._module
    @module.setter
    def module(self, value):
        self._module = value
    @module.deleter
    def module(self):
        del self._module


    # ============


# end class

# ========================================================================= #
# ============================== Biot-Savart ============================== #
# ========================================================================= #

# @[0-9]{4} DATA
class BiotSavart(__LibstellBase):
    """
        an interface to the Biot-Savart module of libstell
    """
    _module = 'biotsavart'

    __params = [
        'nfp_bs',               # integer, number of field periods
        'coil_group',           # structure element with the coil / current information
        'single_coil',          # like coil_group but only a single coil
        ]
    __funcs = [
        'afield',                # Calc. A-vector cylindrical coords (rp, phi, zp, ar, ap az, ig), calculate ar, ap, az
        'bfield',                # Calc. B-vector cylindrical coords (rp, phi, zp, br, bp bz, ig), calculate br, bp, bz
        'cleanup_biotsavart',    # deallocate coil information from memory
        'initialize_biotsavart', # read the coil file (calls parse_coils_file) and deallocates. Outputs the "coil_group"
        'parse_coils_file',     # read the coil file (calls read_coils_pass1 / 2), writing to the "coil_group" variable
        'read_coils_pass1',     # scan the coil file to find the numebr of coil groups / maximum number of nodes in any one coil
        'read_coils_pass2',     # loop the through the coil file again and append a fil_loop to the coil_group
        'write_coils_file',     # write the coil file given an extension, coils.(extension), from the "coilgroup" in memory
        ]


        # extcur = (_ct.c_int * len(extcur_in))(*extcur_in)
        # extension = _ct.c_char_p
        # xpt
        # scaled


    def initialize_biotsavart(self, extcur_in=[], extension='', xpt=None, scaled=False):
        """
        Reads in a coils file if extension is given, otherwise it generates
        a model coil based on the nodes in 'xpt' variable. The name of the
        coils file needs to be 'coils.(extension)'

        inputs:
            extcur_in - [A or Amp-turn?] - current through each model coil
            extension - [str] - string extension of the coils file to load "coils.(extension)"
            xpt       - x,y,z array [m] - nodes for the coil to generate size(3, nodes)
            scaled    - logical - scale currents in coils by input current?
        outputs:
            None



            during initialization, a switch determines whether the input
            'extension' or 'xpt' are used.
            - If extension is present, then the coil file is loaded, read,
              and the current in each coil is scaled as necessary
            - If xpt is specified it generates a model coil based on the
              input current

        """
        extcur_in = _np.atleast_1d(extcur_in)
        ncoils = len(extcur_in)

        fname = self.vname('initialize_biotsavart')
        initialize_biotsavart = getattr(self.mconf, fname)
        initialize_biotsavart.restype = None

        initialize_biotsavart.argtypes=[
                _ct.c_double*ncoils,                # extcur_in
                _ct.c_char_p,                       # extension
                _ct.c_double*arr_size1*arr_size2,   # xpt
                _ct.c_bool                          # scaled
                ]


        # # ==========
        # initialize_biotsavart.argtypes = []
        # inputs = []

        # initialize_biotsavart.argtypes.append(_ct.c_double*ncoils)      # extcur_in
        # inputs.append(extcur_in.tolist())

        # if extension != '':
        #     initialize_biotsavart.argtypes.append(_ct.c_char_p)
        #     inputs.append([extension])
        # if xpt is not None:
        #     arr_size1 = 3
        #     arr_size2 = _np.shape(xpt, 2)
        #     initialize_biotsavart.argtypes.append(_ct.c_double*arr_size1*arr_size2)
        #     inputs.append([xpt])
        # if scaled:
        #     initialize_biotsavart.argtypes.append(_ct.c_bool)
        #     inputs.append([scaled])
        # # end if




        ftemp = _ct.POINTER(_ct.c_double)
        extcur = _npct.as_array(ftemp.in_dll(libstell,s1+'read_wout_mod_'+s2+'_xm'+s3),(len(extcur_in),1))


        ierr = _ct.c_int(0)
        iopen = _ct.c_int(0)
        read_wout(file.encode('UTF-8'), _ct.byref(ierr), _ct.byref(iopen), len(file))
        # Setup Arrays
        vmec_data={}
        # Check
        if not (ierr.value == 0):
            return vmec_data
        # Logical
        varlist=['lasym','lthreed','lwout_opened']
        for temp in varlist:
            vmec_data[temp]=_ct.c_bool.in_dll(libstell,s1+'read_wout_mod_'+s2+'_'+temp+s3).value
        #

    # end def


    def __inittypes__(self, lstell, mc_handle=c_longlong):
        """
        This subfunction loads up all of the C-types/functions from the library
        ... note the short description for each function that is available

        setup the return typs and argument types for libstell library functions
        and exposes them to self.
        """
        ### vec3
        # vec3 is a custom data time that plays well with matrix math
        vec3 = npct.ndpointer(dtype=_np.float64, ndim=1, flags='CONTIGUOUS')
        self.vec3 = vec3

        # ====================== File loading functions ==================== #

    #     """
    #     asdf
    #     """
    #     # lstell.MCload.restype = None
    #     # lstell.MCload.argtypes = []

    #     """
    #     asdf
    #     """
        lstell.initialize_biotsavart.restype = None
        lstell.initialize_biotsavart.argtypes=[
                _ct.c_double*ncoils,                # extcur_in
                _ct.c_char_p,                       # extension
                _ct.c_double*arr_size1*arr_size2,   # xpt
                _ct.c_bool                          # scaled
                ]

    # end def

# # end def

# # ========================================================================= #
# # ===============================  ez_hdf5 ================================ #
# # ========================================================================= #

# class ez_hdf5(object):



#     def open_hdf5(self):
#         """
#             SUBROUTINE open_hdf5(filename,file_id,ier, lcreate,comm,info)

#         """
#         pass
#     # end def open_hdf5


#     def close_hdf5(self):
#         """
#             SUBROUTINE close_hdf5(file_id,ier)

#         """
#         pass
#     # end def close_hdf5

#     def fid(self):
#         pass
#     # end def fid

#     def plistid(self):
#         pass
#     # end def plistid


#     def read_arr2d_hdf5(self):
#         """
#           SUBROUTINE read_arr2d_hdf5(file_id,var,n1,n2,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)

#         """
#         pass
#     # end def read_arr2d_hdf5

#     def read_arr3d_hdf5(self):
#         """
#           SUBROUTINE read_arr3d_hdf5(file_id,var,n1,n2,n3,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)

#         """
#         pass
#     # end def read_arr3d_hdf5


#     def read_arr4d_hdf5(self):
#         """
#           SUBROUTINE read_arr4d_hdf5(file_id,var,n1,n2,n3,n4,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)

#         """
#         pass
#     # end def read_arr4d_hdf5


#     def read_arr5d_hdf5(self):
#         """
#           SUBROUTINE read_arr5d_hdf5(file_id,var,n1,n2,n3,n4,n5,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)

#         """
#         pass
#     # end def read_arr5d_hdf5


#     def read_arr6d_hdf5(self):
#         """
#           SUBROUTINE read_arr6d_hdf5(file_id,var,n1,n2,n3,n4,n5,n6,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)

#         """
#         pass
#     # end def read_arr6d_hdf5


#     def read_scalar_hdf5(self):
#         """
#           SUBROUTINE read_scalar_hdf5(file_id,var,ierr,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)

#         """
#         pass
#     # end def read_scalar_hdf5


#     def read_vector_hdf5(self):
#         """
#           SUBROUTINE read_vector_hdf5(file_id,var,n,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)

#         """
#         pass
#     # end def read_vector_hdf5


#     def write_arr2d_hdf5(self):
#         """
#           SUBROUTINE write_arr2d_hdf5(file_id,var,n1,n2,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)

#         """
#         pass
#     # end def write_arr2d_hdf5

#     def write_arr3d_hdf5(self):
#         """
#           SUBROUTINE write_arr3d_hdf5(file_id,var,n1,n2,n3,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)

#         """
#         pass
#     # end def write_arr3d_hdf5


#     def write_arr4d_hdf5(self):
#         """
#           SUBROUTINE write_arr4d_hdf5(file_id,var,n1,n2,n3,n4,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)

#         """
#         pass
#     # end def write_arr4d_hdf5


#     def write_arr5d_hdf5(self):
#         """
#           SUBROUTINE write_arr5d_hdf5(file_id,var,n1,n2,n3,n4,n5,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)

#         """
#         pass
#     # end def write_arr5d_hdf5


#     def write_arr6d_hdf5(self):
#         """
#           SUBROUTINE write_arr6d_hdf5(file_id,var,n1,n2,n3,n4,n5,n6,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)

#         """
#         pass
#     # end def write_arr6d_hdf5


#     def write_att_hdf5(self):
#         """
#             SUBROUTINE write_att_hdf5(loc_id,ATT_NAME,ATT,ierr)

#         """
#         pass
#     # end def write_att_hdf5


#     def write_att_int_hdf5(self):
#         """
#           SUBROUTINE write_att_int_hdf5(loc_id,ATT_NAME,ATT,ierr)

#         """
#         pass
#     # end def write_att_int_hdf5


#     def write_scalar_hdf5(self):
#         """
#           SUBROUTINE write_scalar_hdf5(file_id,var,ierr,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)

#         """
#         pass
#     # end def write_scalar_hdf5


#     def write_vector_hdf5(self):
#         """
#           SUBROUTINE write_vector_hdf5(file_id,var,n,ier,BOOVAR,INTVAR,FLTVAR,DBLVAR,ATT,ATT_NAME)

#         """
#         pass
#     # end def write_vector_hdf5

# def write_var_hdf5(self):
#     pass
# # end def write_var_hdf5























