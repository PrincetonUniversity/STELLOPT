# -*- coding: utf-8 -*-
"""

 This module reads a VMEC wout file in netCDF format and calculates the currents.
 The code is heavily base on Sam Lazerson's libstell.py and the famous
 read_wout_mod.f by Steve Hirshman.

 History:
  2017-06-18/jons created file
  2017-06-30/jons added stellarator-asymmetric current calculation

 corresponding author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""
# ===================================================================== #
# ===================================================================== #

from __future__ import absolute_import, with_statement, absolute_import, \
                       division, print_function, unicode_literals

from netCDF4 import Dataset
import numpy as _np
from utils import calc_curr

# ========================================================================= #
# ========================================================================= #


def read_vmec(filename):
#    from netCDF4 import Dataset
#    import numpy as _np
#    from utils import calc_curr

    rootgrp = Dataset(filename, 'r')
    vmec_data=dict()

    vmec_data['ns']=rootgrp['/ns'][0]
    vmec_data['nfp']=rootgrp['/nfp'][0]
    vmec_data['mpol']=rootgrp['/mpol'][0]
    vmec_data['ntor']=rootgrp['/ntor'][0]
    vmec_data['mnmax']=rootgrp['/mnmax'][0]
    vmec_data['mnmax_nyq']=rootgrp['/mnmax_nyq'][0]
    vmec_data['iasym']=rootgrp['/lasym__logical__'][0]
    vmec_data['ierr_vmec']=rootgrp['/ier_flag'][0]
    vmec_data['wb']=rootgrp['/wb'][0]
    vmec_data['wp']=rootgrp['/wp'][0]
    vmec_data['gamma']=rootgrp['/gamma'][0]
    vmec_data['rmax_surf']=rootgrp['/rmax_surf'][0]
    vmec_data['rmin_surf']=rootgrp['/rmin_surf'][0]
    vmec_data['zmax_surf']=rootgrp['/zmax_surf'][0]
    vmec_data['aspect']=rootgrp['/aspect'][0]
    vmec_data['betatot']=rootgrp['/betatotal'][0]
    vmec_data['betapol']=rootgrp['/betapol'][0]
    vmec_data['betator']=rootgrp['/betator'][0]
    vmec_data['betaxis']=rootgrp['/betaxis'][0]
    vmec_data['b0']=rootgrp['/b0'][0]
    vmec_data['version']=rootgrp['/version_'][0]
    vmec_data['IonLarmor']=rootgrp['/IonLarmor'][0]
    vmec_data['VolAvgB']=rootgrp['/volavgB'][0]
    #vmec_data['fsql']=rootgrp['/fsql'][0]
    #vmec_data['fsqr']=rootgrp['/fsqr'][0]
    #vmec_data['fsqz']=rootgrp['/fsqz'][0]
    #vmec_data['ftolv']=rootgrp['/ftolv'][0]
    vmec_data['Aminor']=rootgrp['/Aminor_p'][0]
    vmec_data['Rmajor']=rootgrp['/Rmajor_p'][0]
    vmec_data['Volume']=rootgrp['/volume_p'][0]
    vmec_data['RBtor']=rootgrp['/rbtor'][0]
    vmec_data['RBtor0']=rootgrp['/rbtor0'][0]
    vmec_data['Itor']=rootgrp['/ctor'][0]
    vmec_data['isigng']=rootgrp['/signgs'][0] # sign of sqrt(g)

    # since machsq and pfac is declared in read_wout_mod.f, it is assigned a value
    # at the fortran-to-python interface in libstell.
    # Therefore, create this variable here as well even if no data is present.
    for key in ['machsq', 'pfac']:
        if key in rootgrp.variables.keys():
            vmec_data[key]=rootgrp['/'+key][0]
        else:
            vmec_data[key]=0.0

    # Array values 1D
    ns = vmec_data['ns']
    mnmax = vmec_data['mnmax']
    mnmax_nyq = vmec_data['mnmax_nyq']
    ns_size = (ns,1)
    mn_size = (mnmax,1)
    mnnyq_size = (mnmax_nyq,1)
    vmec_data['iotas']  =_np.reshape(rootgrp['/iotas'][:], ns_size)
    vmec_data['iotaf']  =_np.reshape(rootgrp['/iotaf'][:], ns_size)
    vmec_data['presf']  =_np.reshape(rootgrp['/presf'][:], ns_size)
    vmec_data['phipf']  =_np.reshape(rootgrp['/phipf'][:], ns_size)
    #vmec_data['qfact'] =_np.reshape(rootgrp['/qfact'][:], ns_size) ## was commented out in libstell.py as well
    #vmec_data['chipf'] =_np.reshape(rootgrp['/chipf'][:], ns_size) ## random data ???
    #vmec_data['chi']   =_np.reshape(rootgrp['/chi'][:], ns_size) ## random data ???
    vmec_data['phi']    =_np.reshape(rootgrp['/phi'][:], ns_size)
    vmec_data['mass']   =_np.reshape(rootgrp['/mass'][:], ns_size)
    vmec_data['pres']   =_np.reshape(rootgrp['/pres'][:], ns_size)
    vmec_data['betavol']=_np.reshape(rootgrp['/beta_vol'][:], ns_size)
    vmec_data['xm']     =_np.reshape(rootgrp['/xm'][:], mn_size)
    vmec_data['xn']     =_np.reshape(rootgrp['/xn'][:], mn_size)
    vmec_data['xm_nyq'] =_np.reshape(rootgrp['/xm_nyq'][:], mnnyq_size)
    vmec_data['xn_nyq'] =_np.reshape(rootgrp['/xn_nyq'][:], mnnyq_size)
    vmec_data['phip']   =_np.reshape(rootgrp['/phips'][:], ns_size)
    vmec_data['buco']   =_np.reshape(rootgrp['/buco'][:], ns_size)
    vmec_data['bvco']   =_np.reshape(rootgrp['/bvco'][:], ns_size)
    vmec_data['vp']     =_np.reshape(rootgrp['/vp'][:], ns_size)
    #vmec_data['overr'] =_np.reshape(rootgrp['/over_r'][:], ns_size) ## random data ???
    vmec_data['jcuru']  =_np.reshape(rootgrp['/jcuru'][:], ns_size)
    vmec_data['jcurv']  =_np.reshape(rootgrp['/jcurv'][:], ns_size)
    vmec_data['specw']  =_np.reshape(rootgrp['/specw'][:], ns_size)
    vmec_data['jdotb']  =_np.reshape(rootgrp['/jdotb'][:], ns_size)
    vmec_data['Dmerc']  =_np.reshape(rootgrp['/DMerc'][:], ns_size)
    vmec_data['Dshear'] =_np.reshape(rootgrp['/DShear'][:], ns_size)
    vmec_data['Dwell']  =_np.reshape(rootgrp['/DWell'][:], ns_size)
    vmec_data['Dcurr']  =_np.reshape(rootgrp['/DCurr'][:], ns_size)
    vmec_data['Dgeod']  =_np.reshape(rootgrp['/DGeod'][:], ns_size)
    vmec_data['equif']  =_np.reshape(rootgrp['/equif'][:], ns_size)

    for key in ['qfact', 'chipf', 'chi', 'overr']:
        if key in rootgrp.variables.keys():
            vmec_data[key]=_np.reshape(rootgrp['/'+key][:], ns_size)
        else:
            vmec_data[key]=_np.zeros(ns_size)

    vmec_data['rmnc']=rootgrp['/rmnc'][:,:]
    vmec_data['zmns']=rootgrp['/zmns'][:,:]

    vmec_data['lmns']=rootgrp['/lmns'][:,:]
    vmec_data['bmnc']=rootgrp['/bmnc'][:,:]
    vmec_data['gmnc']=rootgrp['/gmnc'][:,:]

    vmec_data['bsupumnc']=rootgrp['/bsupumnc'][:,:]
    vmec_data['bsupvmnc']=rootgrp['/bsupvmnc'][:,:]
    vmec_data['bsubsmns']=rootgrp['/bsubsmns'][:,:]
    vmec_data['bsubumnc']=rootgrp['/bsubumnc'][:,:]
    vmec_data['bsubvmnc']=rootgrp['/bsubvmnc'][:,:]

    # not tested yet...
    if (vmec_data['iasym']):
        vmec_data['rmns']=rootgrp['/rmns'][:,:]
        vmec_data['zmnc']=rootgrp['/zmnc'][:,:]

        vmec_data['lmnc']=rootgrp['/lmnc'][:,:]
        vmec_data['bmns']=rootgrp['/bmns'][:,:]
        vmec_data['gmns']=rootgrp['/gmns'][:,:]

        vmec_data['bsupumns']=rootgrp['/bsupumns'][:,:]
        vmec_data['bsupvmns']=rootgrp['/bsupvmns'][:,:]
        vmec_data['bsubsmnc']=rootgrp['/bsubsmnc'][:,:]
        vmec_data['bsubumns']=rootgrp['/bsubumns'][:,:]
        vmec_data['bsubvmns']=rootgrp['/bsubvmns'][:,:]
    else:
        vmec_data['rmns']=None
        vmec_data['zmnc']=None

        vmec_data['lmnc']=None
        vmec_data['bmns']=None
        vmec_data['gmns']=None

        vmec_data['bsupumns']=None
        vmec_data['bsupvmns']=None
        vmec_data['bsubsmnc']=None
        vmec_data['bsubumns']=None
        vmec_data['bsubvmns']=None

    # we are finished reading the file
    rootgrp.close()

    # calc currumnc, currvmnc
    # important: do this _BEFORE_ you transform the half-mesh quantities onto
    # the full mesh !!!!!!!!!
    currents=calc_curr(vmec_data)
    vmec_data['currumnc']=currents['currumnc']
    vmec_data['currvmnc']=currents['currvmnc']
    if (vmec_data['iasym']):
        vmec_data['currumns']=currents['currumns']
        vmec_data['currvmns']=currents['currvmns']
    else:
        vmec_data['currumns']=None
        vmec_data['currvmns']=None

    # new definition => adapt to NESCOIL ???
#    vmec_data['xn'] = -vmec_data['xn']

#    # Put on full grid
#    vmec_data['buco'] = h2f(vmec_data['buco'],ns)
#    vmec_data['bvco'] = h2f(vmec_data['bvco'],ns)
#    vmec_data['vp'] = h2f(vmec_data['vp'],ns)
#    vmec_data['overr'] = h2f(vmec_data['overr'],ns)
#    vmec_data['specw'] = h2f(vmec_data['specw'],ns)

#    # Put matrix quantities on full grid
#    for key in ['bmnc','gmnc','lmns','bsupumnc','bsupvmnc','bsubsmns','bsubumnc','bsubvmnc']:
#        vmec_data[key][0,:] = 1.5 * vmec_data[key][1,:] - 0.5 * vmec_data[key][2,:]
#        vmec_data[key][1:ns-2,:] = 0.5 * (vmec_data[key][1:ns-2,:] + vmec_data[key][2:ns-1,:])
#        vmec_data[key][ns-1,:] = 2.0 * vmec_data[key][ns-2,:] - vmec_data[key][ns-3,:]
#    if vmec_data['iasym']:
#        for key in ['bmns','gmns','lmnc','bsupumns','bsupvmns','bsubsmnc','bsubumns','bsubvmns']:
#            vmec_data[key][0,:] = 1.5 * vmec_data[key][1,:] - 0.5 * vmec_data[key][2,:]
#            vmec_data[key][1:ns-2,:] = 0.5 * (vmec_data[key][1:ns-2,:] + vmec_data[key][2:ns-1,:])
#            vmec_data[key][ns-1,:] = 2.0 * vmec_data[key][ns-2,:] - vmec_data[key][ns-3,:]

    return vmec_data

# ========================================================================= #
# ========================================================================= #


def read_vmec_netcdf_old(filname):
    mu0=4*_np.pi*1e-7

    f = read_netcdf(filname, 'strip', 'flipdim')

    # Now fix named fields so they match the old way of doing things
    f.ierr_vmec = f.ierflag
    if f.ierr_vmec != 0:
        return f
    # end if
    f.input_extension = _np.copy(f.inputextension)
    f.mgrid_file = _np.copy(str(f.mgridfile))
    f.rmax_surf = _np.copy(f.rmaxsurf)
    f.rmin_surf = _np.copy(f.rminsurf)
    f.zmax_surf = _np.copy(f.zmaxsurf)
    f.ireconstruct = _np.copy(f.lreconlogical)
    f.imse = -1
    f.itse = -1
    f.RBtor = _np.copy(f.rbtor)
    f.Rmajor = _np.copy(f.Rmajorp)
    f.Aminor = _np.copy(f.Aminorp)
    f.betatot = _np.copy(f.betatotal)
    f.Volume  = _np.copy(f.volumep)
    f.VolAvgB = _np.copy(f.volavgB)
    if hasattr(f, 'betavol'):
        f.beta_vol = _np.copy(f.betavol.T)
    if hasattr(f, 'specw'):
        f.specw = _np.copy(f.specw.T)
    if not hasattr(f, 'iasym'):
        f.iasym = 0
    # end if
    f.iasym = _np.copy(f.lasymlogical)
    f.freeb = _np.copy(f.lfreeblogical)
    f.lfreeb = _np.copy(f.freeb)
    f.Itor = _np.copy(f.ctor)
    f.Dmerc = _np.copy(f.DMerc)
    f.Dwell = _np.copy(f.DWell)
    f.Dshear = _np.copy(f.DShear)
    f.Dcurr = _np.copy(f.DCurr)
    f.Dgeod = _np.copy(f.DGeod)

    # Cast some values
    f.ntor = _np.float64(f.ntor)
    f.mpol = _np.float64(f.mpol)
    f.nfp = _np.float64(f.nfp)
    f.ns = _np.float64(f.ns)

    # Calculate Currents
    f.currumnc = _np.zeros((f.mnmaxnyq,f.ns), dtype=_np.float64)
    f.currvmnc = _np.zeros((f.mnmaxnyq,f.ns), dtype=_np.float64)
    ohs = f.ns-1
    hs = 1.0/_np.float64(ohs)
    ns = f.ns

    shalf = _np.zeros((ns,), dtype=_np.float64)
    sfull = _np.zeros((ns,), dtype=_np.float64)
    for ii in range(1, ns):  # 2:ns
        shalf[ii] = _np.sqrt(hs*(ii-1.5))
        sfull[ii] = _np.sqrt(hs*(ii-1))
    # end for

    js1 = range(2, ns)
    js = range(1, ns-1)
    for mn in range(f.mnmaxnyq):
#        if f.xmnyq%2 == 1:  # mod(f.xmnyq, 2) == 1
        if _np.mod(f.xmnyq, 2) == 1:  # mod(f.xmnyq, 2) == 1
            t1 = 0.5 * (shalf[js1] * f.bsubsmns[mn, js1] +
                shalf[js]*f.bsubsmns[mn,js]) / sfull[js]
            bu0 = f.bsubumnc[mn,js] / shalf[js]
            bu1 = f.bsubumnc[mn,js1] / shalf[js1]

            t2 = ohs *(bu1-bu0) * sfull[js] + 0.25*(bu0+bu1)/sfull[js]
            bv0 = f.bsubvmnc[mn,js] / shalf[js]
            bv1 = f.bsubvmnc[mn,js1] / shalf[js1]
            t3 = ohs *(bv1-bv0) * sfull[js] + 0.25*(bv0+bv1)/sfull[js]
        else:
            t1 = 0.5*(f.bsubsmns[mn, js1] + f.bsubsmns[mn, js])
            t2 = ohs*(f.bsubumnc[mn, js1] - f.bsubumnc[mn, js])
            t3 = ohs*(f.bsubvmnc[mn, js1] - f.bsubvmnc[mn, js])
        # end if
        f.currumnc[mn, js] = -_np.float64(f.xnnyq[mn]) * t1 - t3
        f.currvmnc[mn, js] = -_np.float64(f.xmnyq[mn]) * t1 + t2
    # end for
#    # OLD way
#    for ii in range(1,f.ns-1): # i=2:f.ns-1
#        f.currumnc[:, ii] = -_np.float64(f.xnnyq).T * f.bsubsmns[:, ii] - (f.ns-1) * (f.bsubvmnc[:, ii+1] - f.bsubvmnc[:, ii] )
#        f.currvmnc[:, ii] = -_np.float64(f.xmnyq).T * f.bsubsmns[:, ii] + (f.ns-1) * (f.bsubumnc[:, ii+1] - f.bsubumnc[:, ii] )
#    # end

    f.currumnc[:, 0] = 0.0
    f.currvmnc[:, 0] = 0.0
    for ii in range(f.mnmaxnyq): # i=1:f.mnmaxnyq
        if ( f.xmnyq[ii]<=1 ):
            f.currumnc[ii,0] = 2*f.currumnc[ii,1] - f.currumnc[ii,2]
            f.currvmnc[ii,0] = 2*f.currvmnc[ii,1] - f.currvmnc[ii,2]
        # end if
    # end for
    f.currumnc[:, -1] = 2*f.currumnc[:, -2] - f.currumnc[:, -3]
    f.currvmnc[:, -1] = 2*f.currvmnc[:, -2] - f.currvmnc[:, -3]
    f.currumnc /= mu0
    f.currvmnc /= mu0
    if f.iasym:
        f.currumns=_np.zeros( (f.mnmaxnyq,f.ns), dtype=_np.float64)
        f.currvmns=_np.zeros( (f.mnmaxnyq,f.ns), dtype=_np.float64)
        for mn in range(f.mnmaxnyq): # i=1:f.mnmaxnyq
#            if (f.xmnyq%2 == 1):  # mod(f.xmnyq, 2) == 1
            if _np.mod(f.xmnyq, 2) == 1:  # mod(f.xmnyq, 2) == 1
                t1 = 0.5 * (shalf[js1] * f.bsubsmnc[mn, js1] +
                    shalf[js]*f.bsubsmnc[mn,js]) / sfull[js]
                bu0 = f.bsubumns[mn,js] / shalf[js]
                bu1 = f.bsubumns[mn,js1] / shalf[js1]

                t2 = ohs *(bu1-bu0) * sfull[js] + 0.25*(bu0+bu1)/sfull[js]
                bv0 = f.bsubvmns[mn,js] / shalf[js]
                bv1 = f.bsubvmns[mn,js1] / shalf[js1]
                t3 = ohs *(bv1-bv0) * sfull[js] + 0.25*(bv0+bv1)/sfull[js]
            else:
                t1 = 0.5*(f.bsubsmnc[mn, js1] + f.bsubsmnc[mn, js])
                t2 = ohs*(f.bsubumns[mn, js1] - f.bsubumns[mn, js])
                t3 = ohs*(f.bsubvmns[mn, js1] - f.bsubvmns[mn, js])
            # end if
            f.currumns[mn, js] = -_np.float64(f.xnnyq[mn]) * t1 - t3
            f.currvmns[mn, js] = -_np.float64(f.xmnyq[mn]) * t1 + t2
        # end for
#        # OLD way
#        for ii in range(1,f.ns-1): # i=2:f.ns-1
#            f.currumns[:, ii] = -_np.float64(f.xnnyq).T * f.bsubsmnc[:, ii] - (f.ns-1) * (f.bsubvmns[:, ii+1] - f.bsubvmns[:, ii] )
#            f.currvmns[:, ii] = -_np.float64(f.xmnyq).T * f.bsubsmnc[:, ii] + (f.ns-1) * (f.bsubumns[:, ii+1] - f.bsubumns[:, ii] )
#        # end

        f.currumns[:, 0] = 0.0
        f.currvmns[:, 0] = 0.0
        for ii in range(f.mnmaxnyq): # i=1:f.mnmaxnyq
            if (f.xmnyq[ii]<=1):
                f.currumns[ii,0] = 2*f.currumns[ii,1] - f.currumns[ii,2]
                f.currvmns[ii,0] = 2*f.currvmns[ii,1] - f.currvmns[ii,2]
            # end if
        # end for
        f.currumns[:, -1] = 2*f.currumns[:, -2] - f.currumns[:, -3]
        f.currvmns[:, -1] = 2*f.currvmns[:, -2] - f.currvmns[:, -3]
        f.currumns /= mu0
        f.currvmns /= mu0
    # end if iasym

    # Remove Renamed Fields
    del f.inputextension # rmfield(self, 'inputextension')
    del f.mgridfile # rmfield(self, 'mgridfile')
    del f.ierflag # rmfield(self, 'ierflag')
    del f.Rmajorp # rmfield(self, 'Rmajorp')
    del f.Aminorp # rmfield(self, 'Aminorp')
    del f.betatotal # rmfield(self, 'betatotal')
    del f.volumep # rmfield(self, 'volumep')
    del f.volavgB # rmfield(self, 'volavgB')
    del f.betavol # rmfield(self, 'betavol')
    del f.lasymlogical # rmfield(self, 'lasymlogical')
    del f.lfreeblogical # rmfield(self, 'lfreeblogical')
    del f.lreconlogical # rmfield(self, 'lreconlogical')
#    f.close()
    return f
# end  def read_netcdf_old()

# ========================================================================= #
# ========================================================================= #









