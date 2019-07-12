#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 18 14:14:10 2017

@author: jonathan
"""

import numpy as np
from libstell.libstell import read_vmec as libstell_read_vmec
from read_vmec import read_vmec as pyVMEC_read_vmec

# file to use for comparison
# stellarator-symmetric
#filename='/home/IPP-HGW/jons/wout_w7x_ref_60.nc'

# stellarator-asymmetric
filename='/mnt/jons/03_Masterarbeit/aug_asy.example01.01.01/wout_aug_asy.example01.01.01.nc'

# wether to print the positions at which mismatches occur
print_mismatch_positions=True

### do not change below here ######

libstell_data=libstell_read_vmec(filename)

# use python interface to wout from pyVMEC
pyVMEC_data=pyVMEC_read_vmec(filename)

all_ok=True

for key in libstell_data.keys():
    
    #print("checking key "+key)
    
    # ignore these, since they contain random data
    # maybe an error in read_wout_mod.f?
    if key not in ['overr', 'chi', 'chipf']:

        if (key in pyVMEC_data.keys()):
            # check for shape
            ls_shape=np.shape(libstell_data[key])
            pv_shape=np.shape(pyVMEC_data[key])
            if ls_shape != pv_shape:
                print("shape mismatch in "+key+": expected "+np.str(ls_shape)+", got "+np.str(pv_shape))
                all_ok=False
                
            # check for values
            ls_data=libstell_data[key]
            pv_data=pyVMEC_data[key]
            if np.any(ls_data != pv_data):
                print("mismatch in "+key)
                if print_mismatch_positions:
                    print("at positions")
                    print(np.str(np.transpose(np.where(libstell_data[key]!=pyVMEC_data[key]))))
                all_ok=False
                
        else:
            print("key "+key+" not found")
            all_ok=False
            
if all_ok:
    print("all ok :-)")