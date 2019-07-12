#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 16:29:23 2017

@author: jons
"""

import numpy as np
from osa import Client

vmec = Client('http://esb.ipp-hgw.mpg.de:8280/services/vmec_v6?wsdl')

#outpath = "/mnt/jons/03_Masterarbeit/data/all_wout/"
outpath = "Z:\\03_Masterarbeit\\data\\"

# load id list
#failed_ids=np.loadtxt("pyVMEC_failed_ids.txt", dtype='str')
#for vmec_id in failed_ids:
    
#all_ids = vmec.service.listIdentifiers().VmecIds
all_ids=["w7x_ref_163"]
for vmec_id in all_ids:
    if vmec.service.isReady(vmec_id) and vmec.service.wasSuccessful(vmec_id):
        print(vmec_id)
        
        outfile=outpath+"wout_"+vmec_id+".nc"
        
        # read VMEC wout file onto disk
        wout_bytes=vmec.service.getVmecOutputNetcdf(np.str(vmec_id))
        newFileByteArray = bytearray(wout_bytes)
        tempfile=open(outfile, "wb")
        tempfile.write(newFileByteArray)
        tempfile.close()