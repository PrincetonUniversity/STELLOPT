#!/usr/bin/env python3
import sys, os
#import shelve
import json
from pathlib import Path
from argparse import ArgumentParser
sys.path.insert(0, '../../pySTEL/')
import numpy as np
_BENCH_FILE_ = 'BENCH_DATA.json'
# This part changes
from libstell.vmec import VMEC

# Main routine
if __name__=="__main__":
    # Parse Arguments
    parser = ArgumentParser(description= 
        '''Provides tool for benchmarking files.''')
    parser.add_argument("--file", dest="run_name",
        help="File to process", default = None)
    parser.add_argument("--make", dest="lmake_db", action='store_true',
        help="Update the database file.", default = False)

    #
    args = parser.parse_args()
    run_name = args.run_name

    # Do nothing if no filename
    if not run_name: sys.exit(0)
    lfail = False
    failtol = 5.0
    vmec = VMEC()
    vmec.read_wout(f'wout_{run_name}.nc')

    # Extract values
    data={}
    data['aspect']=vmec.aspect
    data['b0']=vmec.b0
    data['wp']=vmec.wp
    data['betatot']=vmec.betatot
    data['volume']=vmec.volume
    data['presf']=vmec.presf.flatten().tolist()
    data['iotaf']=vmec.iotaf.flatten().tolist()
    data['jcurv']=vmec.jcurv.flatten().tolist()
    data['rmnc']=vmec.rmnc[round(vmec.ns/2),:].flatten().tolist()
    data['zmns']=vmec.zmns[round(vmec.ns/2),:].flatten().tolist()
    data['lmns']=vmec.lmns[round(vmec.ns/2),:].flatten().tolist()
    data['bmnc']=vmec.bmnc[round(vmec.ns/2),:].flatten().tolist()

    # Read or write to the database file.
    if args.lmake_db:
        data_out = {}
        my_file = Path(_BENCH_FILE_)
        if my_file.exists():
            f = open(_BENCH_FILE_,"r")
            data_out = json.load(f)
            f.close()
        f = open(_BENCH_FILE_,"w")
        data_out[run_name] = data
        json.dump(data_out, f, ensure_ascii=False, indent=4)
        f.close()
        print('  ADDED: '+run_name)
        sys.exit(0)
    else:
        f = open(_BENCH_FILE_,"r")
        d = json.load(f)
        varlist = d[run_name]
        f.close()
#    if args.lmake_db:
#        my_file = Path(_BENCH_FILE_)
#        if my_file.exists():
#            shelf = shelve.open(_BENCH_FILE_, flag="a")
#        else:
#            shelf = shelve.open(_BENCH_FILE_, flag="c")
#        shelf[run_name] = data
#        shelf.close()
#        print('  ADDED: '+run_name)
#        sys.exit(0)
#    else:
#        shelf = shelve.open(_BENCH_FILE_, flag='r')
#        varlist = shelf[run_name]

    print(f'VMEC VERSION: {vmec.version_:4.2f}')
    for temp in varlist:
        act = varlist[temp]
        cal = data[temp]
        if np.isscalar(act) == 1:
            perct = 100*abs(act-cal)/act
            print(f'  {temp} {cal:7.6f} {act:7.6f} {round(perct)}')
        else:
            cal[act==0] = 0.0
            div = act
            div[div==0] = 1.0
            print(f'  Quantity: {temp} -- CODE -- REF. -- %')
            for i in range(len(act)):
                perct = 100.0*abs(act[i]-cal[i])/div[i]
                print(f'  {i} {cal[i]:7.6f} {act[i]:7.6f} {round(perct)}')
        if perct > failtol:
            lfail = True
        print('=================')

    # Error Status
    if lfail:
        print('  STATUS: FAIL!!!!!')
        sys.exit(-1)
    else:
        print('  STATUS: PASS')
        sys.exit(0)





