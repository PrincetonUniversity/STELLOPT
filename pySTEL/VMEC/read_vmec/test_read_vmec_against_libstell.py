# -*- coding: utf-8 -*-

from osa import Client
import numpy as np
from read_vmec import read_vmec
from libstell.libstell import read_vmec as read_vmec_libstell
import sys, os

vmec = Client('http://esb.ipp-hgw.mpg.de:8280/services/vmec_v7?wsdl')

TOL = 1e-12
DIFF=0.0
TOTAL_DIFF=0.0

failed_ids=[]

filename="/home/IPP-HGW/jons/temp_wout.nc"
#filename="/mnt/jons/03_Masterarbeit/data/w7x_ref_84/wout_w7x_ref_84.nc"


all_vmec_runs = vmec.service.listIdentifiers().VmecIds
num_ids = len(all_vmec_runs)
i=1
#all_vmec_runs = ['w7x_ref_84']
for vmec_id in all_vmec_runs:
    sys.stdout.write(np.str(i)+ " / " + np.str(num_ids) + " is "+vmec_id + ": ")
    sys.stdout.flush()
    if vmec.service.isReady(vmec_id) and vmec.service.wasSuccessful(vmec_id):
        DIFF=0.0
        
        # read VMEC wout file onto disk
        wout_bytes=vmec.service.getVmecOutputNetcdf(vmec_id)
        newFileByteArray = bytearray(wout_bytes)
        os.remove(filename)
        tempfile=open(filename, "wb")
        tempfile.write(newFileByteArray)
        #print(tempfile.tell())
        tempfile.close()
        
        # import data from wout
        vmec_data_reduced=read_vmec(filename)
        vmec_data_libstell=read_vmec_libstell(filename)
        

        currumnc_libstell = vmec_data_libstell['currumnc']
        currvmnc_libstell = vmec_data_libstell['currvmnc']
        
        currumnc_reduced = vmec_data_reduced['currumnc']
        currvmnc_reduced = vmec_data_reduced['currvmnc']
        
        if (np.shape(currumnc_reduced) == np.shape(currumnc_libstell)):
            DIFF += np.sum(np.abs(currumnc_reduced-currumnc_libstell))
            DIFF += np.sum(np.abs(currvmnc_reduced-currvmnc_libstell))
                
            # check also for stellarator-asymmetric terms
            if vmec_data_libstell['iasym']:
                
                currumns_libstell = vmec_data_libstell['currumns']
                currvmns_libstell = vmec_data_libstell['currvmns']
                
                currumns_reduced = vmec_data_reduced['currumns']
                currvmns_reduced = vmec_data_reduced['currvmns']
                
                DIFF += np.sum(np.abs(currumns_reduced-currumns_libstell))
                DIFF += np.sum(np.abs(currvmns_reduced-currvmns_libstell))
        
            if (DIFF > TOL):
                sys.stdout.write("\x1b[36;41m missed by "+np.str(DIFF)+" \x1b[0m\n")
                failed_ids.append(vmec_id)
            else:
                sys.stdout.write("\x1b[37;42m ok \x1b[0m\n")
        else:
            sys.stdout.write("\x1b[36;41m shape mismatch: "+np.str(np.shape(currumnc_reduced))+" != "+np.str(np.shape(currumnc_libstell))+" \x1b[0m\n")
            failed_ids.append(vmec_id)
        sys.stdout.flush()
        
        TOTAL_DIFF += DIFF
    else:
        sys.stdout.write("run failed, not testing...\n")
        sys.stdout.flush()
    i=i+1
            
if (TOTAL_DIFF<TOL):
    print("all ok!")
else:
    print("failed:")
    print(failed_ids)
            
            
            
            