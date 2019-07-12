# -*- coding: utf-8 -*-

from osa import Client
import numpy as np
from read_vmec import read_vmec
import sys, os

# some arbitrary temporary file
temp_filename = "wout_temp.nc"

vmec = Client('http://esb.ipp-hgw.mpg.de:8280/services/vmec_v7?wsdl')

TOL = 1e-12
DIFF=0.0
TOTAL_DIFF=0.0

failed_ids=[]

def vmec_to_ptr(ns, mpol, ntor, vmec_in):
    result = np.zeros([mpol, 2*ntor+1, ns])
    mn=0

    for js in np.arange(ns):
        for m in np.arange(mpol):
            if m==0:
                n0 = ntor
                nn = ntor+1
            else:
                n0 = 0
                nn = 2*ntor+1
            
            for n in np.arange(nn):
                mn = m*nn - ntor + n + n0
                
                result[m][n+n0][js] = vmec_in[js, mn]

    return result	

all_vmec_runs = vmec.service.listIdentifiers().VmecIds
num_ids = len(all_vmec_runs)
i=1
#all_vmec_runs = ['w7x_ref_84']
for vmec_id in all_vmec_runs:
    sys.stdout.write(np.str(i)+ " / " + np.str(num_ids) + " is "+vmec_id + ": ")
    sys.stdout.flush()
    if vmec.service.isReady(vmec_id) and vmec.service.wasSuccessful(vmec_id):
        DIFF=0.0
        
        # read VMEC wout file into memory
        wout_bytes=vmec.service.getVmecOutputNetcdf(vmec_id)
        newFileByteArray = bytearray(wout_bytes)
        
        # write into temporary file
        if (os.path.isfile(temp_filename)):
            os.remove(temp_filename)
        tempfile=open(temp_filename, "wb")
        tempfile.write(newFileByteArray)
        tempfile.close()
        
        # import data from wout
        vmec_data=read_vmec(temp_filename)
        
        ns = vmec_data['ns']
        nfp = vmec_data['nfp']
        mpol = vmec_data['mpol']
        ntor = vmec_data['ntor']
        mnmax = vmec_data['mnmax']
        
        currumnc = vmec_to_ptr(ns, mpol, ntor, vmec_data['currumnc'])
        currvmnc = vmec_to_ptr(ns, mpol, ntor, vmec_data['currvmnc'])
        
        # import current coefficients from webservice method
        CurrUCos = vmec.service.getFourierCoefficients(vmec_id, 'CurrUCos')
        CurrVCos = vmec.service.getFourierCoefficients(vmec_id, 'CurrVCos')
        
        numPol = len(CurrUCos.poloidalModeNumbers)
        numTor = len(CurrUCos.toroidalModeNumbers)
        numRad = CurrUCos.numRadialPoints
        
        # unflatten coeficcients as gotten from webservice
        CurrUCos_PTR = np.ndarray(shape=(numPol,numTor, numRad), dtype='double')
        CurrVCos_PTR = np.ndarray(shape=(numPol,numTor, numRad), dtype='double')
        
        index = 0
        for p in range(0, numPol) :
            for t in range(0, numTor):
                for r in range(0,numRad):
                    CurrUCos_PTR[p][t][r] = CurrUCos.coefficients[index]
                    CurrVCos_PTR[p][t][r] = CurrVCos.coefficients[index]
                    index = index + 1
        
        # check if shape and contents match
        if (np.shape(currumnc) == np.shape(CurrUCos_PTR)):
            DIFF += np.sum(np.abs(currumnc-CurrUCos_PTR))
            DIFF += np.sum(np.abs(currvmnc-CurrVCos_PTR))
                
            # check also for stellarator-asymmetric terms
            if vmec_data['iasym']:
                currumns = vmec_to_ptr(ns, mpol, ntor, vmec_data['currumns'])
                currvmns = vmec_to_ptr(ns, mpol, ntor, vmec_data['currvmns'])
                
                CurrUSin = vmec.service.getFourierCoefficients(vmec_id, 'CurrUSin')
                CurrVSin = vmec.service.getFourierCoefficients(vmec_id, 'CurrVSin')
            
                CurrUSin_PTR = np.ndarray(shape=(numPol,numTor, numRad), dtype='double')
                CurrVSin_PTR = np.ndarray(shape=(numPol,numTor, numRad), dtype='double')
            
                index = 0
                for p in range(0, numPol) :
                    for t in range(0, numTor):
                        for r in range(0,numRad):
                            CurrUSin_PTR[p][t][r] = CurrUSin.coefficients[index]
                            CurrVSin_PTR[p][t][r] = CurrVSin.coefficients[index]
                            index = index + 1
            
                DIFF += np.sum(np.abs(currumns-CurrUSin_PTR))
                DIFF += np.sum(np.abs(currvmns-CurrVSin_PTR))
        
            if (DIFF > TOL):
                sys.stdout.write("\x1b[36;41m missed by "+np.str(DIFF)+" \x1b[0m\n")
                failed_ids.append(vmec_id)
            else:
                sys.stdout.write("\x1b[37;42m ok \x1b[0m\n")
        else:
            sys.stdout.write("\x1b[36;41m shape mismatch: "+np.str(np.shape(currumnc))+" != "+np.str(np.shape(CurrUCos_PTR))+" \x1b[0m\n")
            failed_ids.append(vmec_id)
        sys.stdout.flush()
        
        TOTAL_DIFF += DIFF
    else:
        sys.stdout.write("run failed, not testing...\n")
        sys.stdout.flush()
    i=i+1

# remove leftovers    
os.remove(temp_filename)
            
if (TOTAL_DIFF<TOL):
    sys.stdout.write("\n")
    sys.stdout.flush()
    print("all ok!")
else:
    print("failed:")
    print(failed_ids)
            
            
            
            