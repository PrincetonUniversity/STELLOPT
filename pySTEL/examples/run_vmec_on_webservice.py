# -*- coding: utf-8 -*-
# Jonathan Schilling <jonathan.schilling@ipp.mpg.de>
# 2019-04-30

# run a VMEC calculation on the webservice

import time
from osa import Client
client =  Client("http://esb.ipp-hgw.mpg.de:8280/services/vmec_v8?wsdl")

# input file (&INDATA namelist)
inputfile = "/home/jonathan/Uni/04_PhD/01_analysis/11_ECCD_crashes/input.20180816_022_realCurrents"

# id under which the run should appear on the webservice
vmec_id = "20180816_022_realCurrents"

# check if the run already exists
runExists = client.service.vmecIdentifierExists(vmec_id)
if runExists:
    print("run '"+vmec_id+"' already exists on server")
else:
    # read input file
    with open(inputfile, "r") as f:
        content_input = f.read()
        
        # execute VMEC on webservice
        print("starting VMEC run '"+vmec_id+"'...")
        client.service.execVmecString(content_input, vmec_id)
        
        # wait for run to finish
        while not client.service.isReady(vmec_id):
            time.sleep(1)
        print("VMEC run '"+vmec_id+"' done!")
        