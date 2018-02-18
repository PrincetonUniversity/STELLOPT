# LIBSTELL Module

def read_vmec(file):
    import os, sys
    import ctypes as ct
    import numpy.ctypeslib as npct
    import numpy as np
    # Load Libraries
    try:
        libstell = ct.cdll.LoadLibrary(os.environ["STELLOPT_PATH"]+"/LIBSTELL/Release/libstell.so")
        qtCreatorPath=os.environ["STELLOPT_PATH"]
    except KeyError:
        print("Please set environment variable STELLOPT_PATH")
        sys.exit(1)
    #libstell = ct.cdll.LoadLibrary("/u/slazerso/src/STELLOPT_GCC/LIBSTELL/Release/libstell.so")
    #libstell = ct.cdll.LoadLibrary("/u/slazerso/bin/libstell.so")
    #libstell = ct.cdll.LoadLibrary("/home/jonathan/bin/libstell.so")
    # Read File
    read_wout = getattr(libstell,'__read_wout_mod_MOD_readw_and_open')
    read_wout.argparse=[ct.c_char_p, ct.c_int, ct.c_int, ct.c_int]
    read_wout.restype=None
    ierr = ct.c_int(0)
    iopen = ct.c_int(0)
    read_wout(file.encode('UTF-8'), ct.byref(ierr), iopen, len(file))
    # Setup Arrays
    vmec_data={}
    vmec_data['ns']=ct.c_int.in_dll(libstell,'__read_wout_mod_MOD_ns').value
    vmec_data['nfp']=ct.c_int.in_dll(libstell,'__read_wout_mod_MOD_nfp').value
    vmec_data['mpol']=ct.c_int.in_dll(libstell,'__read_wout_mod_MOD_mpol').value
    vmec_data['ntor']=ct.c_int.in_dll(libstell,'__read_wout_mod_MOD_ntor').value
    vmec_data['mnmax']=ct.c_int.in_dll(libstell,'__read_wout_mod_MOD_mnmax').value
    vmec_data['mnmax_nyq']=ct.c_int.in_dll(libstell,'__read_wout_mod_MOD_mnmax_nyq').value
    vmec_data['iasym']=ct.c_int.in_dll(libstell,'__read_wout_mod_MOD_iasym').value
    vmec_data['ierr_vmec']=ct.c_int.in_dll(libstell,'__read_wout_mod_MOD_ierr_vmec').value
    vmec_data['wb']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_wb').value
    vmec_data['wp']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_wp').value
    vmec_data['gamma']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_gamma').value
    vmec_data['pfac']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_pfac').value
    vmec_data['rmax_surf']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_rmax_surf').value
    vmec_data['rmin_surf']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_rmin_surf').value
    vmec_data['zmax_surf']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_zmax_surf').value
    vmec_data['aspect']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_aspect').value
    vmec_data['betatot']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_betatot').value
    vmec_data['betapol']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_betapol').value
    vmec_data['betator']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_betator').value
    vmec_data['betaxis']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_betaxis').value
    vmec_data['b0']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_b0').value
    vmec_data['version']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_version_').value
    vmec_data['IonLarmor']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_ionlarmor').value
    vmec_data['VolAvgB']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_volavgb').value
    vmec_data['fsql']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_fsql').value
    vmec_data['fsqr']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_fsqr').value
    vmec_data['fsqz']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_fsqz').value
    vmec_data['ftolv']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_ftolv').value
    vmec_data['Aminor']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_aminor').value
    vmec_data['Rmajor']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_rmajor').value
    vmec_data['Volume']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_volume').value
    vmec_data['RBtor']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_rbtor').value
    vmec_data['RBtor0']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_rbtor0').value
    vmec_data['Itor']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_itor').value
    vmec_data['machsq']=ct.c_double.in_dll(libstell,'__read_wout_mod_MOD_machsq').value
    ## Array values 1D
    ftemp=ct.POINTER(ct.c_double)
    ns = vmec_data['ns']
    mnmax = vmec_data['mnmax']
    mnmax_nyq = vmec_data['mnmax_nyq']
    ns_size = (ns,1)
    mn_size = (mnmax,1)
    mnnyq_size = (mnmax_nyq,1)
    vmec_data['iotas']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_iotas'),ns_size)
    vmec_data['iotaf']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_iotaf'),ns_size)
    vmec_data['presf']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_presf'),ns_size)
    vmec_data['phipf']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_phipf'),ns_size)
    #vmec_data['qfact']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_qfact'),ns_size)
    vmec_data['chipf']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_chipf'),ns_size)
    vmec_data['chi']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_chi'),ns_size)
    vmec_data['phi']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_phi'),ns_size)
    vmec_data['mass']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_mass'),ns_size)
    vmec_data['pres']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_pres'),ns_size)
    vmec_data['betavol']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_beta_vol'),ns_size)
    vmec_data['xm']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_xm'),mn_size)
    vmec_data['xn']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_xn'),mn_size)
    vmec_data['xm_nyq']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_xm_nyq'),mnnyq_size)
    vmec_data['xn_nyq']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_xn_nyq'),mnnyq_size)
    vmec_data['phip']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_phip'),ns_size)
    vmec_data['buco']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_buco'),ns_size)
    vmec_data['bvco']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_bvco'),ns_size)
    vmec_data['vp']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_vp'),ns_size)
    vmec_data['overr']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_overr'),ns_size)
    vmec_data['jcuru']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_jcuru'),ns_size)
    vmec_data['jcurv']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_jcurv'),ns_size)
    vmec_data['specw']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_specw'),ns_size)
    vmec_data['jdotb']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_jdotb'),ns_size)
    vmec_data['Dmerc']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_dmerc'),ns_size)
    vmec_data['Dshear']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_dshear'),ns_size)
    vmec_data['Dwell']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_dwell'),ns_size)
    vmec_data['Dcurr']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_dcurr'),ns_size)
    vmec_data['Dgeod']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_dgeod'),ns_size)
    vmec_data['equif']=npct.as_array(ftemp.in_dll(libstell,'__read_wout_mod_MOD_equif'),ns_size)
    ## 2D Arrays
    mn2d_size = (ns, mnmax)
    mn2d_nyq_size = (ns, mnmax_nyq)
    fmn=ct.POINTER(ct.c_double)
    vmec_data['rmnc']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_rmnc'),mn2d_size) #ns,mnmax format
    vmec_data['zmns']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_zmns'),mn2d_size) #ns,mnmax format
    vmec_data['lmns']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_lmns'),mn2d_size) #ns,mnmax format
    vmec_data['bmnc']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_bmnc'),mn2d_nyq_size) #ns,mnmax format
    vmec_data['gmnc']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_gmnc'),mn2d_nyq_size) #ns,mnmax format
    vmec_data['bsupumnc']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_bsupumnc'),mn2d_nyq_size) #ns,mnmax format
    vmec_data['bsupvmnc']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_bsupvmnc'),mn2d_nyq_size) #ns,mnmax format
    vmec_data['bsubsmns']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_bsubsmns'),mn2d_nyq_size) #ns,mnmax format
    vmec_data['bsubumnc']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_bsubumnc'),mn2d_nyq_size) #ns,mnmax format
    vmec_data['bsubvmnc']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_bsubvmnc'),mn2d_nyq_size) #ns,mnmax format
    vmec_data['currumnc']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_currumnc'),mn2d_nyq_size) #ns,mnmax format
    vmec_data['currvmnc']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_currvmnc'),mn2d_nyq_size) #ns,mnmax format
    if vmec_data['iasym']:
        vmec_data['rmns']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_rmns'),mn2d_size) #ns,mnmax format
        vmec_data['zmnc']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_zmnc'),mn2d_size) #ns,mnmax format
        vmec_data['lmnc']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_lmnc'),mn2d_size) #ns,mnmax format
        vmec_data['bmns']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_bmns'),mn2d_nyq_size) #ns,mnmax format
        vmec_data['gmns']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_gmns'),mn2d_nyq_size) #ns,mnmax format
        vmec_data['bsupumns']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_bsupumns'),mn2d_nyq_size) #ns,mnmax format
        vmec_data['bsupvmns']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_bsupvmns'),mn2d_nyq_size) #ns,mnmax format
        vmec_data['bsubsmnc']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_bsubsmnc'),mn2d_nyq_size) #ns,mnmax format
        vmec_data['bsubumns']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_bsubumns'),mn2d_nyq_size) #ns,mnmax format
        vmec_data['bsubvmns']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_bsubvmns'),mn2d_nyq_size) #ns,mnmax format
        vmec_data['currumns']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_currumns'),mn2d_nyq_size) #ns,mnmax format
        vmec_data['currvmns']=npct.as_array(fmn.in_dll(libstell,'__read_wout_mod_MOD_currvmns'),mn2d_nyq_size) #ns,mnmax format
    # Free memory (don't do this as python accesses this memory)
    #read_wout_dealloc = getattr(libstell,'__read_wout_mod_MOD_read_wout_deallocate')
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
    import numpy as np
    temp = np.zeros((ns,1))
    temp[0] = 1.5 * var[0] - 0.5 * var[1]
    temp[1:ns-1] = 0.5* (var[0:ns-2] + var[1:ns-1])
    temp[ns-1] = 1.5 * var[ns-2] - 0.5 * var[ns-3]
    return temp

def cfunct(theta,zeta,fmnc,xm,xn):
    import numpy as np
    f=0
    (ns,mn)=fmnc.shape
    lt = len(theta)
    lz = len(zeta)
    mt=np.matmul(xm,theta.T)
    nz=np.matmul(xn,zeta.T)
    cosmt=np.cos(mt)
    sinmt=np.sin(mt)
    cosnz=np.cos(nz)
    sinnz=np.sin(nz)
    f = np.zeros((ns,lt,lz))
    
    fmn = np.ndarray((mn,lt))
    for k in range(ns):
        fmn = np.broadcast_to(fmnc[k,:],(lt,mn)).T
        fmncosmt=(fmn*cosmt).T
        fmnsinmt=(fmn*sinmt).T
        f[k,:,:]=np.matmul(fmncosmt, cosnz)-np.matmul(fmnsinmt, sinnz)
    return f
    
def sfunct(theta,zeta,fmnc,xm,xn):
    import numpy as np
    f=0
    (ns,mn)=fmnc.shape
    lt = len(theta)
    lz = len(zeta)
    mt=np.matmul(xm,theta.T)
    nz=np.matmul(xn,zeta.T)
    cosmt=np.cos(mt)
    sinmt=np.sin(mt)
    cosnz=np.cos(nz)
    sinnz=np.sin(nz)
    f = np.zeros((ns,lt,lz))
    fmn = np.ndarray((mn,lt))
    for k in range(ns):
        fmn = np.broadcast_to(fmnc[k,:],(lt,mn)).T
        f[k,:,:]=np.matmul((fmn*sinmt).T,cosnz)+np.matmul((fmn*cosmt).T,sinnz)
    return f

def torocont(r,z,val,s):
    import numpy as np
    import matplotlib.pyplot as pyplot
    h=pyplot.axes(xlabel='R [m]',ylabel='Z [m]',aspect='equal')
    pyplot.pcolormesh(r[:,:,s],z[:,:,s],val[:,:,s],cmap='jet',shading='gouraud',axes=h)
    pyplot.show()
    return h

def toroslice(r,zeta,z,s):
    import numpy as np
    import matplotlib.pyplot as pyplot
    h=pyplot.axes(xlabel='R [m]',ylabel='Z [m]',aspect='equal')
    if (s[0] == 0):
        pyplot.plot(r[0,0,zeta],z[0,0,zeta],'+',color='black',axes=h)
        pyplot.plot(np.transpose(r[s[1:],:,zeta]),np.transpose(z[s[1:],:,zeta]),color='black',axes=h)
    else:
        for i in range(0,):
            pyplot.plot(np.transpose(r[s,:,zeta]),np.transpose(z[s,:,zeta]),color='black',axes=h)
    pyplot.show()
    return h

def isotoro(r,z,zeta,svals,*args,**kwargs):
    import numpy as np
    import matplotlib.pyplot as pyplot
    import mpl_toolkits.mplot3d as mplot3d
    import math as math
    import matplotlib.tri as mtri
    from mayavi import mlab
    nr = np.size(svals)
    if (nr == 1):
        s= [svals]
        nr = 1
    else:
        s=svals
    nt = np.size(r,1)
    nz = np.size(r,2)
    vertex = np.zeros((nt*nz,3,nr))
    for k in range(0,nr):
        ivertex = 0
        ifaces = 0
        for j in range(0,nz):
            for i in range(0,nt):
                vertex[ivertex,0,k]=r[s[k],i,j]*math.cos(zeta[j])
                vertex[ivertex,1,k]=r[s[k],i,j]*math.sin(zeta[j])
                vertex[ivertex,2,k]=z[s[k],i,j]
                ivertex = ivertex + 1
    u = np.linspace(0, 1, endpoint=True, num=nt)
    v = np.linspace(0, 1, endpoint=True, num=nz)
    u, v = np.meshgrid(u, v)
    u, v = u.flatten(), v.flatten()
    tri = mtri.Triangulation(u, v)
    test=len(kwargs)
    fig=kwargs.pop('fig',pyplot.figure())
    h=kwargs.pop('axes',fig.add_subplot(111,projection='3d'))
    for k in range(0,nr):
        if (len(args)==0):
            tsurf=h.plot_trisurf(vertex[:,0,k],vertex[:,1,k],vertex[:,2,k], triangles=tri.triangles,color='red',shade='yes',linewidths=0.0)
            #tsurf=mlab.triangular_mesh(vertex[:,0,k],vertex[:,1,k],vertex[:,2,k], tri.triangless)
        else:
            # Matplotlib way (SLOW)
            vals = args[0][s[k],:,:].T.flatten()
            colors = np.mean(vals[tri.triangles], axis=1)
            tsurf=h.plot_trisurf(vertex[:,0,k],vertex[:,1,k],vertex[:,2,k], triangles=tri.triangles,cmap='jet',shade='yes',linewidths=0.0)
            tsurf.set_array(colors)
            tsurf.autoscale()
            #MAYAVI Way (need to figure out how to embed)
            #vals = args[0][s[k],:,:].T.flatten()
            #tsurf=mlab.triangular_mesh(vertex[:,0,k],vertex[:,1,k],vertex[:,2,k], tri.triangles, scalars=vals, colormap='jet')
    if (test==0):
        pyplot.show()
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
    #      theta=np.linspace(0, 2*np.pi, 360)
    #      zeta=np.linspace(0, 2*np.pi, 64)
    #      vmec_data=read_vmec('wout.nc')
    #      jll=calc_jll(vmec_data, theta, zeta)
    
    
    # Maintained by: Samuel Lazerson (lazerson@pppl.gov)
    # Version:       1.00
    
    b =cfunct(theta,zeta,vmec_data['bmnc'],    vmec_data['xm'],vmec_data['xn'])
    g =cfunct(theta,zeta,vmec_data['gmnc'],    vmec_data['xm'],vmec_data['xn'])
    bu=cfunct(theta,zeta,vmec_data['bsubumnc'],vmec_data['xm'],vmec_data['xn'])
    bv=cfunct(theta,zeta,vmec_data['bsubvmnc'],vmec_data['xm'],vmec_data['xn'])
    ju=cfunct(theta,zeta,vmec_data['currumnc'],vmec_data['xm'],vmec_data['xn'])
    jv=cfunct(theta,zeta,vmec_data['currvmnc'],vmec_data['xm'],vmec_data['xn'])
    
    if (vmec_data['iasym']):
        b =b +sfunct(theta,zeta,vmec_data['bmns'],    vmec_data['xm'],vmec_data['xn'])
        g =g +sfunct(theta,zeta,vmec_data['gmns'],    vmec_data['xm'],vmec_data['xn'])
        bu=bu+sfunct(theta,zeta,vmec_data['bsubumns'],vmec_data['xm'],vmec_data['xn'])
        bv=bv+sfunct(theta,zeta,vmec_data['bsubvmns'],vmec_data['xm'],vmec_data['xn'])
        ju=ju+sfunct(theta,zeta,vmec_data['currumns'],vmec_data['xm'],vmec_data['xn'])
        jv=jv+sfunct(theta,zeta,vmec_data['currvmns'],vmec_data['xm'],vmec_data['xn'])
    
    
    jll = (bu*ju+bv*jv)/(g*b)
    return jll



    
    
    
#def safe_open(file,iunit):
#    import ctypes as ct
#    # Load Libraries
#    libstell = ct.cdll.LoadLibrary("/Users/slazerso/Sims_PPPL/STELLOPT/LIBSTELL/Release/libstell.so")
#    read_input = getattr(libstell,'__safe_open_mod_MOD_safe_open')
#    read_input.argparse=[ct.c_int,ct.c_int,ct.c_char_p,ct.c_int,ct.c_char_p,ct.c_int,ct.c_char_p,ct.c_int]
#    read_input.restype=None
#    ierr = ct.c_int(0)
#    istat = ct.c_int(0)
#    read_input(ct.byref(ierr),ct.byref(istat),file.encode('UTF-8'))


#def read_indata(file):
#    import ctypes as ct
#    import numpy.ctypeslib as npct
#    import numpy as np
#    # Load Libraries
#    libstell = ct.cdll.LoadLibrary("/Users/slazerso/Sims_PPPL/STELLOPT/LIBSTELL/Release/libstell.so")
#    # Read File
#    read_input = getattr(libstell,'__vmec_input_MOD_read_indata_namelist')
#    read_wout.argparse=[ct.c_int,ct.c_int]
#    read_wout.restype=None
#    ierr = ct.c_int(0)
#    iopen = ct.c_int(0)
#    read_wout(file.encode('UTF-8'),ct.byref(ierr),iopen,len(file))
