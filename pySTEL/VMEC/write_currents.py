# -*- coding: utf-8 -*-

import numpy as np
from netCDF4 import Dataset
from libstell.libstell import read_vmec as libstell_read_vmec

from osa import Client

vmec = Client('http://esb.ipp-hgw.mpg.de:8280/services/vmec_v6?wsdl')

# stellarator-symmetric
#in_filename='/mnt/jons/03_Masterarbeit/wout_w7x_ref_60.nc'
#out_filename='/mnt/jons/03_Masterarbeit/libstell_currents_w7x_ref_60.nc'

# stellarator-asymmetric
#in_filename='/mnt/jons/03_Masterarbeit/aug_asy.example01.01.01/wout_aug_asy.example01.01.01.nc'
#out_filename='/mnt/jons/03_Masterarbeit/aug_asy.example01.01.01/libstell_currents_aug_asy.example01.01.01.nc'

outpath = "/mnt/jons/03_Masterarbeit/data/all_wout/"

# load id list
#failed_ids=np.loadtxt("pyVMEC_failed_ids.txt", dtype='str')
#for vmec_id in failed_ids:
all_ids = vmec.service.listIdentifiers().VmecIds
for vmec_id in all_ids:
    if vmec.service.isReady(vmec_id) and vmec.service.wasSuccessful(vmec_id):
        print(vmec_id)
        
        wout_file=outpath+"wout_"+vmec_id+".nc"
        out_filename = outpath+"ls_currents_"+vmec_id+".nc"
        
        # read data using the libstell imported routine
        #libstell_data=libstell_read_vmec(in_filename)
        libstell_data=libstell_read_vmec(wout_file)
        
        # create the new Dataset to store the currents in
        ds_currents=Dataset(out_filename, 'w')
        
        # create Dimensions of the data to be stored
        ds_currents.createDimension("ns", libstell_data['ns'])
        ds_currents.createDimension("mnmax_nyq", libstell_data['mnmax_nyq'])
        
        # create Variables to store the currents in 
        ds_currents.createVariable("currumnc", np.double, ['ns', 'mnmax_nyq'], fill_value=0.0)
        ds_currents.createVariable("currvmnc", np.double, ['ns', 'mnmax_nyq'], fill_value=0.0)
        
        # write the currents into the Variables
        ds_currents['/currumnc'][:,:]=libstell_data['currumnc']
        ds_currents['/currvmnc'][:,:]=libstell_data['currvmnc']
        
        # stellarator-asymmetric currents
        if (libstell_data['iasym']):
            # create Variables to store the currents in 
            ds_currents.createVariable("currumns", np.double, ['ns', 'mnmax_nyq'], fill_value=0.0)
            ds_currents.createVariable("currvmns", np.double, ['ns', 'mnmax_nyq'], fill_value=0.0)
            
            # write the currents into the Variables
            ds_currents['/currumns'][:,:]=libstell_data['currumns']
            ds_currents['/currvmns'][:,:]=libstell_data['currvmns']
        
        # write the data to the disk and close the file
        ds_currents.close()