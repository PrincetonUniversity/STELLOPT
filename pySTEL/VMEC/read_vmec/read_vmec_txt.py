# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 20:42:27 2019

    converted from
    from stellopt\LIBSTELL\Sources\Modules\read_wout_mod.f90
@author: weir

note that numpy uses the opposite indexing to fortran, so we have to transpose
arrays when using numpy
"""
# ===================================================================== #
# ===================================================================== #

#from __future__ import absolute_import, with_statement, absolute_import, \
#                       division, print_function, unicode_literals
import numpy as _np
from utils import calc_curr

# ===================================================================== #
# ===================================================================== #


def read_vmec(filename, delim=' ', verbose=True):
    try:
        fid = open(filename, 'r')

        vmec_data = read_open_wout_text(fid, delim=delim, verbose=verbose)
        fid.close()
    except:
        raise
    finally:
        try:    fid.close()
        except: pass
    # end try
    return vmec_data


def read_in_constant(vmec_data, iunit, keys, dtype, delim=' ', verbose=True):
    if dtype[keys[0]] == float:
        dtype = float
    elif dtype[keys[0]] == int:
        dtype = int
    # end if
#    pos = iunit.tell()

    data = _np.fromfile(iunit, dtype=dtype, count=len(keys), sep=delim)
#    data = _np.fromfile(iunit, dtype=dt, sep=delim)  # breaking somehow
#    for key in keys:
    for ii, key in enumerate(keys):
        try:
            vmec_data[key] = _np.copy(data[ii])
    #        vmec_data[key] = _np.copy(data[key])
        except:
            if verbose:
                print('missing data! for %s'%(key,))
    #        pass
            raise
        # end try
    # end for
#    pos2 = iunit.tell()
#    pass
# end def

def read_in_index(vmec_data, iunit, keys, dtype, index=0, delim=' ', verbose=True):
    if dtype[keys[0]] == float:
        dtype = float
    elif dtype[keys[0]] == int:
        dtype = int
    # end if

#    data = _np.fromfile(iunit, dtype=dt, sep=delim)  # breaking somehow
    data = _np.fromfile(iunit, dtype=dtype, count=len(keys), sep=delim)

#    if len(data)==0:
#        break
#    for key in keys:
    for ii, key in enumerate(keys):
        try:
            vmec_data[key][index] = _np.copy(data[ii])
    #        vmec_data[key][index] = _np.copy(data[key])
        except:
            if verbose:
                print('missing data! for %s'%(key,))
                print(('index:',index))
            pass
    #        raise
        # end try
    # end for
# end def

def read_open_wout_text(iunit, delim=' ', verbose=True):
    vmec_data = {}
#------------------------------------------------
#   D u m m y   A r g u m e n t s
#------------------------------------------------
#      INTEGER :: iunit, ierr
#------------------------------------------------
#   L o c a l   P a r a m e t e r s
#------------------------------------------------
    eps_w = 1.e-4
    mu0 = 4.0*_np.pi*1.0e-7
#------------------------------------------------
#   L o c a l   V a r i a b l e s
#------------------------------------------------
    i, j, k, js, m, n, n1, mn, nparts_in, i_animec, i_flow = \
        tuple(int(0) for ii in range(11))   # 464
    ierr = 0
    istat = _np.zeros((15,), dtype=int)

    vmec_version = ''
#------------------------------------------------
#
#     THIS SUBROUTINE READS THE TEXT FILE WOUT CREATED BY THE VMEC CODE
#     AND STORES THE INFORMATION IN THE read_WOUT MODULE
#
#     CALL read_wout_file - GENERIC INTERFACE - CAN BE CALLED WITH EITHER UNIT NO. OR FILENAME
#
#     RMNC, ZMNS: FULL-GRID
#     LMNS      : HALF-GRID
#
    vmec_data['lcurr'] = False
    vmec_data['ierr_vmec'] = 0
    vmec_data['nextcur'] = 0

    # ===================================================================== #
    # Read the version info off of the first line
    vmec_version = iunit.readline()

    ii = vmec_version.find('=')
    # === ADDED BY SAL
    i_animec = vmec_version.find('_ANIMEC')
    i_flow = vmec_version.find('_FLOW')
    vmec_type = 0
    if (i_animec > -1):
       vmec_type = 1 # ANIMEC
       vmec_version = vmec_version[range(i_animec-1)]
    # end if
    if (i_flow > -1):
       vmec_type = 2 # FLOW
       vmec_version = vmec_version[range(i_flow-1)]
    # end if
    # === END SAL Addition

    if (ii >= 0):
       version_ = vmec_version[ii+1:vmec_version.find('\n')]
    else:
       version_ = -1.0
    # end if
    version_ = float(version_)
    vmec_data['vmec_version'] = version_

    # ===================================================================== #
#    data = iunit.readline().strip().split(delim)
#    data = [dat.strip() for dat in data if len(dat)>0]
    if (version_ <= (5.10 + eps_w)):
        # ==================== version less than 5.10 ==================== #
        keys = ['wb', 'wp', 'gamma', 'pfac', 'nfp', 'ns', 'mpol', 'ntor', 'mnmax',
                'itfsq', 'niter', 'iasym', 'ireconstruct']

#        dt = _np.dtype([(key, _np.float) if ii<4 else (key, _np.int)
#            for ii, key in enumerate(keys)])
        data = _np.fromfile(iunit, dtype=float, count=len(keys), sep=delim)
        for ii, key in enumerate(keys):
            vmec_data[key] = data[ii]
            if ii>=4:
                vmec_data[key] = int(vmec_data[key])
            # end if
        # end for

    else:
        # ==================== version less than 6.54 ==================== #
        keys = ['wb','wp','gamma','pfac','rmax_surf','rmin_surf']

        # ================== version greater than 6.54 =================== #
        if (version_ >= 6.54-eps_w):
            keys.append('zmax_surf')
        # end if

        dt = _np.dtype([(key, _np.float) for key in keys])
        read_in_constant(vmec_data, iunit, keys, dtype=dt, delim=delim, verbose=verbose)

        # ======================================== #

        if (vmec_type == 2): # === SAL
            vmec_data['machsq'] = _np.fromfile(iunit, dtype=_np.float, count=1, sep=delim)
        # end if

        if (version_ <= (8.0+eps_w)):
            # ================== version less than 8.00 =================== #
            keys = ['nfp', 'ns', 'mpol', 'ntor', 'mnmax', 'itfsq', 'niter',
                    'iasym', 'ireconstruct', 'ierr_vmec']
        else:
            # ================= version greater than 8.00 ================= #
            keys = ['nfp', 'ns', 'mpol', 'ntor', 'mnmax', 'mnmax_nyq', 'itfsq', 'niter',
                    'iasym', 'ireconstruct', 'ierr_vmec']
        # end if

        dt = _np.dtype([(key, _np.int) for key in keys])
        read_in_constant(vmec_data, iunit, keys, dtype=dt, delim=delim, verbose=verbose)

        if 'mnmax_nyq' not in vmec_data:
            # mnmax_nyq only in vmec versions greater than 8.00
            vmec_data['mnmax_nyq'] = vmec_data['mnmax']
        # end if
    # end if
    lasym = (vmec_data['iasym'] > 0)

    # ===================================================================== #

    keys = ['imse', 'itse', 'nbsets', 'nobd', 'nextcur']
    # ================= version greater than 6.20 ================= #
    if (version_ > (6.20+eps_w)):
        keys.append('nstore_seq')
    # end if

    dt = _np.dtype([(key, _np.int) for key in keys])
    read_in_constant(vmec_data, iunit, keys, dtype=dt, delim=delim, verbose=verbose)

    # ================= version up to 6.20 ================= #
    if 'nstore_seq' not in vmec_data:
        vmec_data['nstore_seq'] = 100
    # end if

    # ===================================================================== #

    def Continue1000(ierr=int(0)):
        data = iunit.readline().strip().split(delim)
        data = [dat.strip() for dat in data if len(dat)>0]
        if len(data)==0 or len(data[0])==0:
            vmec_data['mgrid_mode'] = 'N'
            return True, ierr
        vmec_data['mgrid_mode'] = data

        if (ierr != 0):   # analysis:ignore
            ierr = 0
            vmec_data['mgrid_mode'] = 'N'
        if (istat[1] != 0):
            vmec_data['ierr_vmec'] = 1
        # end if

        exit_program = False
        for m in range(15):
            if (istat[m] > 0):
                print(' Error No.%i '%(m+1,)+' in READ_WOUT, iostat = %i'%(istat[m],))
                ierr = m
                exit_program = True
                return exit_program, ierr
            # end if
        # end for
        return exit_program, ierr
    # end def Continue1000

#    if ((vmec_data['ierr_vmec'] != norm_term_flag) and \
#        (vmec_data['ierr_vmec'] != more_iter_flag) and \
#        (vmec_data['ierr_vmec'] != jac75_flag)):
#        ierr = -2
#        exit_flag, ierr = Continue1000(ierr)
#        if exit_flag:
#            return vmec_data
#        # end if
#    # end if

    # ===================================================================== #
    # ==================================================================== #
#    if (vmec_data['nextcur'] > vmec_data['nigroup']):
#        istat[15-1] = -1
#    # end if
    # pre-allocate the first part of arrays to read in
    AllocateVMECArrays(vmec_data, vmec_type, lasym, verbose=verbose)

    # ==================================================================== #

    if vmec_data['nbsets'] >= 0:
        vmec_data['nbfld'] = _np.fromfile(iunit, dtype=_np.float, count=vmec_data['nbsets'], sep=delim)
    # end if
    vmec_data['mgrid_file'] = iunit.readline().strip()

    for js in range(vmec_data['ns']):   # js = 1, ns
        for mn in range(vmec_data['mnmax']):  # 1, mnmax
            if js == 0:
#                dt = _np.dtype(('m',_np.int), ('n', _np.int))
                m, n = tuple(_np.fromfile(iunit, dtype=int, count=2, sep=delim))
                vmec_data['xm'][mn] = float(m)
                vmec_data['xn'][mn] = float(n)
            # end if

            if (version_ <= (6.20+eps_w)):
                keys = ['rmnc', 'zmns', 'lmns', 'bmnc', 'gmnc', 'bsubumnc',
                       'bsubvmnc', 'bsubsmns', 'bsupumnc', 'bsupvmnc', 'currvmnc']
                #fmt = 730: '5e20.13'
            elif (version_ <= (8.0+eps_w)):
                keys = ['rmnc', 'zmns', 'lmns', 'bmnc', 'gmnc', 'bsubumnc',
                       'bsubvmnc', 'bsubsmns', 'bsupumnc', 'bsupvmnc', 'currvmnc']
                #fmt = '*'
            else:
                keys = ['rmnc', 'zmns', 'lmns']
                #fmt = '*'
            # end if

            dt = _np.dtype([(key, _np.float) for key in keys])
            read_in_index(vmec_data, iunit, keys, dtype=dt, index=(mn,js), delim=delim, verbose=verbose)

            if lasym:
                if (version_ <= (8.0+eps_w)):
                    keys = ['rmns', 'zmnc', 'lmnc', 'bmns', 'gmns', 'bsubumns',
                            'bsubvmns', 'bsubsmnc', 'bsupumns', 'bsupvmns']
                    #fmt = '*'
                else:
                    keys = ['rmns', 'zmnc', 'lmnc']
                    #fmt = '*'
                # end if
                dt = _np.dtype([(key, _np.float) for key in keys])
                read_in_index(vmec_data, iunit, keys, dtype=dt, index=(mn,js), delim=delim, verbose=verbose)
            # end if
            if (js == 0) and (m ==0):
                n1 = int(abs(n/vmec_data['nfp']))-1
                if n1 <= vmec_data['ntor']:
                    vmec_data['raxis'][n1, 0] = _np.copy(vmec_data['rmnc'][mn, 0])
                    vmec_data['zaxis'][n1, 0] = _np.copy(vmec_data['zmns'][mn, 0])
                    if lasym:
                        vmec_data['raxis'][n1, 1] = _np.copy(vmec_data['rmns'][mn, 0])
                        vmec_data['zaxis'][n1, 1] = _np.copy(vmec_data['zmnc'][mn, 0])
                    # end if
                # end if
            # end if
        # end do

        if (version_ <= (8.0+eps_w)):
            continue
        # end if

        for mn in range(vmec_data['mnmax_nyq']):
            if js == 0:
                m, n = tuple(_np.fromfile(iunit, dtype=int, count=2, sep=delim))
                vmec_data['xm_nyq'][mn] = float(m)
                vmec_data['xn_nyq'][mn] = float(n)
            # end if

            if vmec_type == 1:  # SAL (ELSE statement below is orriginal)
                keys = ['bmnc', 'gmnc', 'bsubumnc', 'bsubvmnc', 'bsubsmns',
                        'bsupumnc', 'bsupvmnc', 'pparmnc', 'ppermnc', 'hotdmnc'
                        'pbprmnc', 'ppprmnc', 'sigmnc', 'taumnc']

                dt = _np.dtype([(key, _np.float) for key in keys])
                read_in_index(vmec_data, iunit, keys, dtype=dt, index=(mn,js), delim=delim, verbose=verbose)

                if lasym:
                    keys = ['bmns', 'gmns', 'bsubumns', 'bsubvmns', 'bsubsmnc',
                            'bsupumns', 'bsupvmns', 'pparmns', 'ppermns', 'hotdmns',
                            'pbprmns', 'ppprmns', 'sigmns', 'taumns']

                    dt = _np.dtype([(key, _np.float) for key in keys])
                    read_in_index(vmec_data, iunit, keys, dtype=dt, index=(mn,js), delim=delim, verbose=verbose)
                # end if
            elif vmec_type == 2:
                keys = ['bmnc', 'gmnc', 'bsubumnc', 'bsubvmnc', 'bsubsmns',
                        'bsupumnc', 'bsupvmnc', 'protmnc', 'protrsqmnc', 'prprmnc']

                dt = _np.dtype([(key, _np.float) for key in keys])
                read_in_index(vmec_data, iunit, keys, dtype=dt, index=(mn,js), delim=delim, verbose=verbose)

                if lasym:
                    keys = ['bmns', 'gmns', 'bsubumns','bsubvmns', 'bsubsmnc',
                            'bsupumns', 'bsupvmns', 'protmns', 'protrsqmns', 'prprmns']

                    dt = _np.dtype([(key, _np.float) for key in keys])
                    read_in_index(vmec_data, iunit, keys, dtype=dt, index=(mn,js), delim=delim, verbose=verbose)
                # end if
            else:
                keys = ['bmnc', 'gmnc', 'bsubumnc', 'bsubvmnc', 'bsubsmns',
                        'bsupumnc', 'bsupvmnc']

                dt = _np.dtype([(key, _np.float) for key in keys])
                read_in_index(vmec_data, iunit, keys, dtype=dt, index=(mn,js), delim=delim, verbose=verbose)

                if lasym:
                    keys = ['bmns', 'gmns', 'bsubumns', 'bsubvmns', 'bsubsmnc',
                            'bsupumns', 'bsupvmns']
                    dt = _np.dtype([(key, _np.float) for key in keys])
                    read_in_index(vmec_data, iunit, keys, dtype=dt, index=(mn,js), delim=delim, verbose=verbose)
                # end if
            # end if
        # end for
    # end for

    # Compute current coefficients on full mesh
    if (version_ > 8.0+eps_w):
        # The shared function expects some arrays transposed (fortran ordering) and a few as well
        keys = ['bsubsmns','bsubumnc','bsubvmnc']
        if vmec_data['iasym']:
            keys = keys + ['bsubsmnc','bsubumns','bsubvmns']
        # end if
        for key in keys:
            vmec_data[key] = vmec_data[key].T

        # Calculate the currents
        vmec_data = calc_curr(vmec_data)

        # Now undo that mess
        keys = keys + ['currumnc', 'currvmnc']
        if vmec_data['iasym']:
            keys = keys + ['currumns', 'currvmns']
        for key in keys:
            vmec_data[key] = vmec_data[key].T
        # end if
    # end if

    vmec_data['mnyq'] = int(_np.max(vmec_data['xm_nyq']))
    vmec_data['nnyq'] = int(_np.max(_np.abs(vmec_data['xn_nyq'])))/vmec_data['nfp']

    # ===================================================================== #
    """
     Read FULL AND HALF-MESH QUANTITIES

     NOTE: In version_ <= 6.00, mass, press were written out in INTERNAL (VMEC) units
     and are therefore multiplied here by 1/mu0 to transform to pascals. Same is true
     for ALL the currents (jcuru, jcurv, jdotb). Also, in version_ = 6.10 and
     above, PHI is the true (physical) toroidal flux (has the sign of jacobian correctly
     built into it)
    """
    keys = ['iotas', 'mass', 'pres', 'phip', 'buco', 'bvco', 'vp', 'overr', 'specw', 'beta_vol']
    for key in keys:
        vmec_data[key][0] = 0
    # end for
    if version_ <= (6.05+eps_w):
        keys = ['iotas', 'mass', 'pres', 'phip', 'buco', 'bvco', 'phi', 'vp', 'overr',
                'jcuru', 'jcurv', 'specw']
        # 730 FORMAT(5e20.13)
        dt = _np.dtype([(key, _np.float) for key in keys])
        for js in range(1, vmec_data['ns']):
            read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
        # end for

        keys = ['aspect', 'betatot', 'betapol', 'betator', 'betaxis', 'b0']
        dt = _np.dtype([(key, _np.float) for key in keys])
        read_in_constant(vmec_data, iunit, keys, dtype=dt, delim=delim, verbose=verbose)
    elif version_ <= 6.20+eps_w:
        keys = ['iotas', 'mass', 'pres', 'beta_vol', 'phip', 'buco', 'bvco',
                'phi', 'vp', 'overr', 'jcuru', 'jcurv', 'specw']
        # 730 FORMAT(5e20.13)
        dt = _np.dtype([(key, _np.float) for key in keys])
        for js in range(1, vmec_data['ns']):
            read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
        # end for

        keys = ['aspect', 'betatot', 'betapol', 'betator', 'betaxis', 'b0']
        dt = _np.dtype([(key, _np.float) for key in keys])
        read_in_constant(vmec_data, iunit, keys, dtype=dt, delim=delim, verbose=verbose)
    elif version_ <= (6.95+eps_w):
        keys = ['iotas', 'mass', 'pres', 'beta_vol', 'phip', 'buco', 'bvco',
                'phi', 'vp', 'overr', 'jcuru', 'jcurv', 'specw']
        # *
        dt = _np.dtype([(key, _np.float) for key in keys])
        for js in range(1, vmec_data['ns']):
            read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
        # end for
        keys = ['aspect', 'betatot', 'betapol', 'betator', 'betaxis', 'b0']
        dt = _np.dtype([(key, _np.float) for key in keys])
        read_in_constant(vmec_data, iunit, keys, dtype=dt, delim=delim, verbose=verbose)
    else:
        keys = ['iotaf', 'presf', 'phipf', 'phi', 'jcuru', 'jcurv']
        # *
        dt = _np.dtype([(key, _np.float) for key in keys])
        for js in range(vmec_data['ns']):
            read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
        # end for

        keys = ['iotas', 'mass', 'pres', 'beta_vol', 'phip', 'buco', 'bvco',
                'vp', 'overr', 'specw']
        # *
        dt = _np.dtype([(key, _np.float) for key in keys])
        for js in range(1, vmec_data['ns']):
            read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
        # end for

        keys = ['aspect', 'betatot', 'betapol', 'betator', 'betaxis', 'b0']
        dt = _np.dtype([(key, _np.float) for key in keys])
        read_in_constant(vmec_data, iunit, keys, dtype=dt, delim=delim, verbose=verbose)
     # end if

    if version_ > (6.10+eps_w):
        keys = ['isign']
        dt = _np.dtype([(key, _np.int) for key in keys])
        read_in_constant(vmec_data, iunit, keys, dtype=dt, delim=delim, verbose=verbose)
#        isigng = vmec_data['isign']

        keys = ['input_extension']
        vmec_data[keys[0]] = iunit.readline().strip()

        keys = ['IonLarmor', 'VolAvgB', 'RBtor0', 'RBtor', 'Itor', 'Aminor',
                'Rmajor', 'Volume']
        dt = _np.dtype([(key, _np.float) for key in keys])
        read_in_constant(vmec_data, iunit, keys, dtype=dt, delim=delim, verbose=verbose)
    # end if

    # ============================================== #
    #               MERCIER CRITERION                #
    # ============================================== #
    if (version_ > (5.10+eps_w)) and (version_ < (6.20-eps_w)):
        keys = ['Dmerc', 'Dshear', 'Dwell', 'Dcurr', 'Dgeod', 'equif']
        # fortran fmt 730
        dt = _np.dtype([(key, _np.float) for key in keys])
        for js in range(1, vmec_data['ns']-1):
            read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
        # end for
    elif (version_ >= (6.20-eps_w)):
        keys = ['Dmerc', 'Dshear', 'Dwell', 'Dcurr', 'Dgeod', 'equif']
        # fortran fmt *
        dt = _np.dtype([(key, _np.float) for key in keys])
        for js in range(1, vmec_data['ns']-1):
            read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
        # end for
    # end if

    if (vmec_data['nextcur'] > 0):
        if version_ <= (6.20+eps_w):
            data = _np.fromfile(iunit, dtype=_np.float, count=vmec_data['nextcur'], sep=delim)
            vmec_data['extcur'] = _np.copy(data)
        else:
            data = _np.fromfile(iunit, dtype=_np.float, count=vmec_data['nextcur'], sep=delim)
            vmec_data['extcur'] = _np.copy(data)
        # end if

        data = iunit.readline()
        vmec_data['lcurr'] = True if data.find('T')>-1 else False
        if (vmec_data['lcurr']):
            # Figure out how many current labels exist per line (max)
            ncur = len(vmec_data['extcur'])
            vmec_data['curlabel'] = []
            nremaining = ncur
            tstcount = 0
            while nremaining>0 and tstcount<200:
                data = iunit.readline().strip().split(delim)
                data = [dat.strip() for dat in data if len(dat)>1]
                data = [dat.replace("\'", "") for dat in data]
                # store the ones you have already grabbed
                for js in range(len(data)):
                    vmec_data['curlabel'].append(data[js])
                # end for
                tstcount += 1
                nremaining -= len(data)
            # end while
            vmec_data['curlabel'] = _np.asarray(vmec_data['curlabel'])
        # end if
    # end if

    if (version_ <= (6.20+eps_w)):
        keys = ['fsqt', 'wdot']
        dt = _np.dtype([(key, _np.float) for key in keys])
        for js in range(vmec_data['nstore_seq']):
            read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
        # end for
#         READ (iunit, 730, iostat=istat(14))(fsqt(i), wdot(i), i=1,nstore_seq)
    else:
        keys = ['fsqt', 'wdot']
        dt = _np.dtype([(key, _np.float) for key in keys])
        for js in range(vmec_data['nstore_seq']):
            read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
        # end for
#         READ (iunit, *, iostat=istat(14))(fsqt(i), wdot(i), i=1,nstore_seq)
    # end if

    if (version_>=6.20-eps_w) and (version_ < 6.50-eps_w) and (istat[14-1]==0):
        keys = ['jdotb', 'bdotgradv', 'bdotb']
        dt = _np.dtype([(key, _np.float) for key in keys])
        for js in range(vmec_data['ns']):
            read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
        # end for
#         READ (iunit, 730, iostat=istat(14), err=1000)                  &
#           (jdotb(js), bdotgradv(js), bdotb(js), js=1,ns)
    elif (version_ >= (6.50-eps_w)):
        keys = ['jdotb', 'bdotgradv', 'bdotb']
        dt = _np.dtype([(key, _np.float) for key in keys])
        for js in range(vmec_data['ns']):
            read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
        # end for
#         READ (iunit, *, iostat=istat(14), err=1000)                    &
#           (jdotb(js), bdotgradv(js), bdotb(js), js=1,ns)
    else:
         istat[14-1] = 0
    # end if

    vmec_data['chipf'] = vmec_data['iotaf']*vmec_data['phipf']
    #
    #     CONVERT FROM INTERNAL UNITS TO PHYSICAL UNITS IF NEEDED
    #
    if (version_ <= (6.05+eps_w)):
         vmec_data['mass'] /= mu0
         vmec_data['pres'] /= mu0
         vmec_data['jcuru'] /= mu0
         vmec_data['jcurv'] /= mu0
         vmec_data['jdotb'] /= mu0
         vmec_data['phi']   *= -1.0
    # end if

    #-----------------------------------------------
    #     DATA AND MSE FITS
    #-----------------------------------------------
    if vmec_data['ireconstruct'] > 0:
        n1 = AllocateDATAArrays(vmec_data, vmec_type, lasym, verbose=verbose)

        if vmec_data['imse'] >= 2 or vmec_data['itse'] > 0:
            keys = ['tswgt', 'msewgt']
            dt = _np.dtype([(key, _np.float) for key in keys])
            read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
#            READ (iunit, *) tswgt, msewgt

            vmec_data['isnodes'] = iunit.fromfile(dtype=_np.int, count=1, sep=delim)
            keys = ['sknots', 'ystark', 'y2stark']
            dt = _np.dtype([(key, _np.float) for key in keys])
            for js in range(vmec_data['isnodes']):
                read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
            # end for
#            READ (iunit, *) isnodes, (sknots(i),ystark(i),y2stark(i),   &
#               i=1,isnodes)

            vmec_data['ipnodes'] = iunit.fromfile(dtype=_np.dtype('ipnodes', _np.int), count=1, sep=delim)
            keys = ['pknots', 'ythom', 'y2thom']
            dt = _np.dtype([(key, _np.float) for key in keys])
            for js in range(vmec_data['ipnodes']):
                read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
            # end for
#            READ (iunit, *) ipnodes, (pknots(i), ythom(i),              &
#               y2thom(i),i=1,ipnodes)

            keys = ['anglemse', 'rmid', 'qmid', 'shear', 'presmid', 'alfa', 'curmid']
            dt = _np.dtype([(key, _np.float) for key in keys])
            for js in range(2*vmec_data['ns'-1]):
                read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
            # end for
#            READ(iunit, *)(anglemse(i),rmid(i),qmid(i),shear(i),        &
#            presmid(i),alfa(i),curmid(i),i=1,2*ns-1)

            keys = ['rstark', 'datastark', 'qmeas']
            dt = _np.dtype([(key, _np.float) for key in keys])
            for js in range(vmec_data['imse']):
                read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
            # end for
#            READ(iunit, *)(rstark(i),datastark(i),qmeas(i),i=1,imse)

            keys = ['rtom', 'datathom']
            dt = _np.dtype([(key, _np.float) for key in keys])
            for js in range(vmec_data['itse']):
                read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
            # end for
#            READ(iunit, *)(rthom(i),datathom(i),i=1,itse)
        # end if

        if vmec_data['nobd']>0:
            keys = ['dsiext', 'plflux', 'dsiobt']
            dt = _np.dtype([(key, _np.float) for key in keys])
            for js in range(vmec_data['nobd']):
                read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
            # end for
#            READ (iunit, *) (dsiext(i),plflux(i),dsiobt(i),i=1,nobd)

            keys = ['flmwgt']
            dt = _np.dtype([(key, _np.float) for key in keys])
            read_in_constant(vmec_data, iunit, keys, dtype=dt, delim=delim, verbose=verbose)
#            READ (iunit, *) flmwgt
         # end if

        vmec_data['nbfldn'] = _np.sum(vmec_data['nbfld'][:vmec_data['nbsets']])
        if (vmec_data['nbfldn'] > 0):
            keys = ['bcoil', 'plbfld', 'bbc']
            dt = _np.dtype([(key, _np.float) for key in keys])
            for n in range(vmec_data['nbsets']):  # 1, nbsets
                for jn in range(vmec_data['nbfldn'][n]):
                    read_in_index(vmec_data, iunit, keys, dtype=dt, index=(n,jn), delim=delim, verbose=verbose)
                # end for
#               READ (iunit, *) (bcoil(i,n),plbfld(i,n),bbc(i,n),        &
#                  i=1,nbfld(n))
            # end for

            keys = ['bcwgt']
            dt = _np.dtype([(key, _np.float) for key in keys])
            read_in_constant(vmec_data, iunit, keys, dtype=dt, delim=delim, verbose=verbose)
            # READ (iunit, *) bcwgt
        # end if

        keys = ['phidiam', 'delphid']
        dt = _np.dtype([(key, _np.float) for key in keys])
        read_in_constant(vmec_data, iunit, keys, dtype=dt, delim=delim, verbose=verbose)
#         READ (iunit, *) phidiam, delphid

        #
        #     READ Limiter & Prout plotting specs
        #
        keys = ['nsets', 'npars_in', 'nlim']
        dt = _np.dtype([(key, _np.int) for key in keys])
        read_in_constant(vmec_data, iunit, keys, dtype=dt, delim=delim, verbose=verbose)
#         READ (iunit, *) nsets, nparts_in, nlim

        vmec_data['nsetsn'] = _np.zeros((vmec_data['nsets'],), dtype=_np.int)
#         ALLOCATE (nsetsn(nsets))
        keys = ['nsetsn']
        dt = _np.dtype([(key, _np.int) for key in keys])
        for js in range(vmec_data['nsets']):
            read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
#         READ (iunit, *) (nsetsn(i),i=1,nsets)
        # end for

        n1 = int(_np.max(vmec_data['nsetsn'][:vmec_data['nsets']]))
        vmec_data['pfcspec'] = _np.zeros((vmec_data['nparts_in'],n1,vmec_data['nsets']), dtype=_np.float)
        vmec_data['limitr'] = _np.zeros((vmec_data['nlim'],1), dtype=_np.float)
        # ALLOCATE (pfcspec(nparts_in,n1,nsets), limitr(nlim))

        keys = ['pfcspec']
        dt = _np.dtype([(key, _np.float64) for key in keys])
        for kk in range(vmec_data['nsets']):
            for jj in range(vmec_data['nsetsn'][kk]):
                for ii in range(vmec_data['nparts_in']):
#                    vmec_data[keys[0]][ii,jj,kk] = _np.fromfile(iunit, dtype=dt, sep=delim)
                    vmec_data[keys[0]][ii,jj,kk] = _np.fromfile(iunit, dtype=float, count=len(keys), sep=delim)
                # end for
            # end for
        # end for
#        READ (iunit, *) (((pfcspec(i,j,k),i=1,nparts_in),              &
#            j=1,nsetsn(k)),k=1,nsets)

        keys = ['limitr']
        dt = _np.dtype([(key, _np.float64) for key in keys])
        for js in range(vmec_data['nlim']):
            read_in_index(vmec_data, iunit, keys, dtype=dt, index=js, delim=delim, verbose=verbose)
        # end for
#         READ (iunit, *) (limitr(i), i=1,nlim)

        m  = _np.max(vmec_data['limitr'][:vmec_data['nlim']])
        vmec_data['rlim'] = _np.zeros((m,vmec_data['nlim']), dtype=_np.float)
        vmec_data['zlim'] = _np.zeros((m,vmec_data['nlim']), dtype=_np.float)
#        ALLOCATE (rlim(m,nlim), zlim(m,nlim))

        keys = ['rlim', 'zlim']
        dt = _np.dtype([(key, _np.float64) for key in keys])
        for jj in range(vmec_data['nlim']):
            for ii in range(vmec_data['limitr'][jj]):
                read_in_index(vmec_data, iunit, keys, dtype=dt, index=(ii,jj), delim=delim, verbose=verbose)
            # end for
        # end for
#         READ (iunit, *) ((rlim(i,j),zlim(i,j),i=1,limitr(j)),          &
#            j=1,nlim)

        keys = ['nrgrid', 'nzgrid']
        dt = _np.dtype([(key, _np.int) for key in keys])
        read_in_constant(vmec_data, iunit, keys, dtype=dt, delim=delim, verbose=verbose)
#         READ (iunit, *) nrgrid, nzgrid

        keys = ['tokid']
        dt = _np.dtype([(key, _np.float) for key in keys])
        read_in_constant(vmec_data, iunit, keys, dtype=dt, delim=delim, verbose=verbose)
#         READ (iunit, *) tokid

        keys = ['rx1', 'rx2', 'zy1', 'zy2', 'condif']
        dt = _np.dtype([(key, _np.float) for key in keys])
        read_in_constant(vmec_data, iunit, keys, dtype=dt, delim=delim, verbose=verbose)
#         READ (iunit, *) rx1, rx2, zy1, zy2, condif

        keys = ['imatch_phiedge']
        dt = _np.dtype([(key, _np.float) for key in keys])
        read_in_constant(vmec_data, iunit, keys, dtype=dt, delim=delim, verbose=verbose)
#         READ (iunit, *) imatch_phiedge
    # end if

    exit_flag, ierr = Continue1000(ierr)

    for key in vmec_data:
        if type(vmec_data[key]) != type('') and len(_np.atleast_1d(vmec_data[key]))>1:
            if len(_np.atleast_1d(vmec_data[key]).shape) == 1:
                vmec_data[key] = _np.atleast_2d(vmec_data[key])  # change (ns,) to (ns,1)
                vmec_data[key] = vmec_data[key].T  # change (1,ns) to (ns,1)
                if verbose:
                    print(('transposing array %s to shape:'%(key,),vmec_data[key].shape))
                # end if
            # end if
            vmec_data[key] = _np.asfortranarray(vmec_data[key])
            if verbose:
                print('Switching array %s to fortran ordering:'%(key,))
            # end if
    return vmec_data

    #  720 FORMAT(8i10)
    #  730 FORMAT(5e20.13)
    #  740 FORMAT(a)
    #  790 FORMAT(i5,/,(1p,3e12.4))
# end def SUBROUTINE read_wout_text

def AllocateVMECArrays(vmec_data, vmec_type, lasym, verbose=True):
    # ========================================================== #
    # ================== Pre-allocation of arrays ============== #
    vmec_data['xm'] = _np.zeros((vmec_data['mnmax'],1), dtype=float)
    vmec_data['xn'] = _np.zeros_like(vmec_data['xm'])
    vmec_data['xm_nyq'] = _np.zeros((vmec_data['mnmax_nyq'],1), dtype=float)
    vmec_data['xn_nyq'] = _np.zeros_like(vmec_data['xm_nyq'])

    keys = ['rmnc','zmns','lmns','bmnc','gmnc','bsubumnc','bsubvmnc','bsubsmns',
            'bsupumnc','bsupvmnc','currvmnc','currumnc']
    for ii, key in enumerate(keys):
        vmec_data[key] = _np.zeros((vmec_data['mnmax_nyq'],vmec_data['ns']), dtype=float)
    # end for

    keys = ['iotas', 'mass', 'pres', 'beta_vol', 'phip', 'buco', 'bvco',
            'phi', 'iotaf', 'presf', 'phipf', 'chipf', 'vp', 'overr',
            'jcuru', 'jcurv', 'specw', 'Dmerc','Dshear', 'Dwell', 'Dcurr',
            'Dgeod', 'equif', 'jdotb', 'bdotb', 'bdotgradv']
    for ii, key in enumerate(keys):
        vmec_data[key] = _np.zeros((vmec_data['ns'],1), dtype=float)
    # end for

    vmec_data['raxis'] = _np.zeros((vmec_data['ntor'],2), dtype=float)
    vmec_data['zaxis'] = _np.zeros_like(vmec_data['raxis'])

    vmec_data['am'], vmec_data['ac'], vmec_data['ai'] \
        = tuple(_np.zeros((20,), dtype=float) for _ in range(3))

    vmec_data['fsqt'] = _np.zeros((vmec_data['nstore_seq'],), dtype=float)
    vmec_data['wdot'] = _np.zeros_like(vmec_data['fsqt'])

#    stat = istat[6-1]

    if vmec_data['nextcur'] > 0:
        vmec_data['extcur'] = _np.zeros((vmec_data['nextcur'],), dtype=float)
        vmec_data['curlabel'] = _np.zeros((vmec_data['nextcur'],), dtype=str)
#        stat = istat[6-1]
    # end if
    if lasym:
        keys = ['rmns', 'zmnc', 'lmnc']
        for key in keys:
            vmec_data[key] = _np.zeros((vmec_data['mnmax'], vmec_data['ns']), dtype=float)
        # end for

        keys = ['bmns', 'gmns', 'bsubumns','bsubvmns','bsubsmnc','bsupumns',
                'bsupvmns','currumns','currvmns']
        for key in keys:
            vmec_data[key] = _np.zeros((vmec_data['mnmax_nyq'], vmec_data['ns']), dtype=float)
        # end for
#        stat=istat[6-1]
    # end if

    if vmec_type == 1:
        keys = ['pparmnc', 'ppermnc', 'hotdmnc', 'pbprmnc', 'ppprmnc', 'sigmnc', 'taumnc']
        for key in keys:
            vmec_data[key] = _np.zeros((vmec_data['mnmax_nyq'],vmec_data['ns']), dtype=float)
        # end for
#        stat=istat[6-1]

        if lasym:
            keys = ['pparmns', 'ppermns', 'hotdmns', 'pbprmns', 'ppprmns', 'sigmns', 'taumns']
            for key in keys:
                vmec_data[key] = _np.zeros((vmec_data['mnmax_nyq'],vmec_data['ns']), dtype=float)
            # end for
#            stat=istat[6-1]
        # end if
    elif vmec_type == 2:
        vmec_data['pmap'] = _np.zeros((vmec_data['ns'],1), dtype=float)
        vmec_data['omega'] = _np.zeros_like(vmec_data['pmap'])
        vmec_data['tpotb'] = _np.zeros_like(vmec_data['pmap'])

        keys = ['protmnc', 'protrsqmnc', 'prprmnc']
        for key in keys:
            vmec_data[key] = _np.zeros((vmec_data['mnmax_nyq'], vmec_data['ns']), dtype=float)
        # end for
#        stat=istat[6-1]

        if lasym:
            keys = ['protmns', 'protrsqmns', 'prprmns']
            for key in keys:
                vmec_data[key] = _np.zeros((vmec_data['mnmax_nyq'],vmec_data['ns']), dtype=float)
            # end for
#            stat=istat[6-1]
        # end if
    # end if

def AllocateDATAArrays(vmec_data, vmec_type, lasym, verbose=True):
    # ========================================================== #
    # ================== Pre-allocation of arrays ============== #

    n1 = int(_np.max(vmec_data['nbfld'][:vmec_data['nbsets']+1]))

#        keys = ['sknots', 'ystark', 'y2stark', 'pknots', 'ythom', 'y2thom']
#        for key in keys:
#            vmec_data[key] = _np.zeros((isnodes,1), dtype=float)
#        # end for
#
#        keys = ['anglemse','rmid', 'qmid', 'shear', 'presmid', 'alfa', 'curmid']
#        for key in keys:
#            vmec_data[key] = _np.zeros((2*vmec_data['ns'],1), dtype=float)
#        # end for
#
#        vmec_data['rstark'] = _np.zeros((imse,1), dtype=float)
#        vmec_data['datastark'] = _np.zeros_like(vmec_data['rstark'])
#
#        vmec_data['rthom'] = _np.zeros((itse,1), dtype=float)
#        vmec_data['datathom'] = _np.zeros_like(vmec_data['rthom'])
#
#        vmec_data['dsiext'] = _np.zeros((nobd,1), dtype=float)
#        vmec_data['plflux'] = _np.zeros_like(vmec_data['dsiext'])
#        vmec_data['dsiobt'] = _np.zeros_like(vmec_data['dsiext'])
#
#        vmec_data['bcoil'] = _np.zeros((n1,nbsets), dtype=float)
#        vmec_data['plbfld'] = _np.zeros_like(vmec_data['bcoil'])
#        vmec_data['bbc'] = _np.zeros_like(vmec_data['bcoil'])

    return n1