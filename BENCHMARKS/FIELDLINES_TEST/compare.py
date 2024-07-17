#!/usr/bin/env python3
import sys, os
import json
from pathlib import Path
from argparse import ArgumentParser
sys.path.insert(0, '../../pySTEL/')
import numpy as np
_BENCH_FILE_ = 'BENCH_DATA.json'
# This part changes
from libstell.fieldlines import FIELDLINES

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
    fieldlines = FIELDLINES()
    try:
        fieldlines.read_fieldlines(f'fieldlines_{run_name}.h5')
    except:
        print(f'  ERROR: Cannot find run {run_name}')
        sys.exit(-1)
    version_str = f'FIELDLINES VERSION: {fieldlines.VERSION:4.2f}'

    # Extract values
    data={}
    if run_name == 'NCSX_s1':
        data['iota0'] = fieldlines.iota0

        [_,iota,_] = fieldlines.calc_iota()
        data['iota'] = iota.tolist()
    if run_name == 'NCSX_s1_coll':
        data['wall_strikes'] = fieldlines.wall_strikes.flatten().tolist()
    data['b_r'] = fieldlines.B_R[:,0,64].flatten().tolist()
    data['b_phi'] = fieldlines.B_PHI[:,0,64].flatten().tolist()
    data['b_z'] = fieldlines.B_Z[:,0,64].flatten().tolist()

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
        f.close()
        if run_name in d.keys():
            varlist = d[run_name]
        else:
            print(f'  ERROR: Cannot find {run_name} in {_BENCH_FILE_}')
            sys.exit(-1)

    print(version_str)
    print('=================')
    for temp in varlist:
        act = varlist[temp]
        cal = data[temp]
        if np.isscalar(act) == 1:
            perct = 100*abs(act-cal)/act
            print(f'  {temp} {cal:7.6f} {act:7.6f} {round(perct)}')
        elif temp == 'wall_strikes':
            act = np.array(act)
            cal = np.array(cal)
            cal = np.where(act==0,0,cal)
            div = np.where(act==0,1,act)
            perct = 100*sum(abs(act-cal)/div)
            print(f'  {temp} {max(cal):7.6f} {max(act):7.6f} {round(perct)}')
        else:
            act = np.array(act)
            cal = np.array(cal)
            cal = np.where(act==0,0,cal)
            div = np.where(act==0,1,act)
            print(f'  Quantity: {temp} -- CODE -- REF. -- %')
            for i in range(len(act)):
                perct = 100*abs(act[i]-cal[i])/div[i]
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





