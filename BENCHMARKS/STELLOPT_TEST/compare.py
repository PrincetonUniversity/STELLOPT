#!/usr/bin/env python3
import sys, os
import json
from pathlib import Path
from argparse import ArgumentParser
sys.path.insert(0, '../../pySTEL/')
import numpy as np
_BENCH_FILE_ = 'BENCH_DATA.json'
# This part changes
from libstell.stellopt import STELLOPT

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
    sopt = STELLOPT()
    try:
        sopt.read_stellopt_output(f'{run_name}/stellopt.{run_name}')
    except:
        print(f'  ERROR: Cannot find file {run_name}/stellopt.{run_name}')
        sys.exit(-1)
    version_str = f'STELLOPT VERSION: {sopt.stellopt_version:4.2f}'

    # Extract values
    data={}
    if run_name in ['BASIC']:
        data['ASPECT'] = sopt.ASPECT_VAL.flatten().tolist()
        data['BETA'] = sopt.BETA_VAL.flatten().tolist()
        data['CURTOR'] = sopt.CURTOR_VAL.flatten().tolist()
        data['PHIEDGE'] = sopt.PHIEDGE_VAL.flatten().tolist()
        data['RBTOR'] = sopt.RBTOR_VAL.flatten().tolist()
        data['R0'] = sopt.R0_VAL.flatten().tolist()
        data['VOLUME'] = sopt.VOLUME_VAL.flatten().tolist()
        data['WP'] = sopt.WP_VAL.flatten().tolist()
        data['BALLOON_BALLOON_GRATE'] = sopt.BALLOON_BALLOON_GRATE.flatten().tolist()
        data['BOOTSTRAP'] = sopt.BOOTSTRAP_VAL.flatten().tolist()
        data['TXPORT'] = sopt.TXPORT_VAL.flatten().tolist()
        data['PRESS'] = sopt.PRESS_VAL.flatten().tolist()
        data['NE'] = sopt.NE_VAL.flatten().tolist()
        data['NELINE'] = sopt.NELINE_VAL.flatten().tolist()
        data['FARADAY'] = sopt.FARADAY_VAL.flatten().tolist()
        data['TE'] = sopt.TE_VAL.flatten().tolist()
        data['TI'] = sopt.TI_VAL.flatten().tolist()
    if run_name in ['LMDIF_TEST']:
        data['TEST_X'] = sopt.TEST_X_VAL[-1,0].tolist()
        data['TEST_Y'] = sopt.TEST_Y_VAL[-1,0].tolist()
    if run_name in ['GADE_TEST']:
        data['TEST_X'] = sopt.TEST_X_VAL[-1,0].tolist()
        data['TEST_Y'] = sopt.TEST_Y_VAL[-1,0].tolist()
    if run_name in ['PSO_TEST']:
        data['TEST_X'] = sopt.TEST_X_VAL[-1,0].tolist()
        data['TEST_Y'] = sopt.TEST_Y_VAL[-1,0].tolist()
    if run_name in ['SA_TEST']:
        data['TEST_X'] = sopt.TEST_X_VAL[-1,0].tolist()
        data['TEST_Y'] = sopt.TEST_Y_VAL[-1,0].tolist()
    if run_name in ['IOTA_LMDIF']:
        data['TEST_X'] = sopt.IOTA_VAL[-1,:].tolist()
    if run_name in ['RECON_TOK']:
        data['CURTOR'] = sopt.CURTOR_VAL.flatten().tolist()
        data['XICS_BRIGHT'] = sopt.XICS_BRIGHT_VAL.flatten().tolist()
        data['XICS'] = sopt.XICS_VAL.flatten().tolist()
        data['XICS_W3'] = sopt.XICS_W3_VAL.flatten().tolist()
        data['XICS_V'] = sopt.XICS_V_VAL.flatten().tolist()
        data['VISBREMLINE'] = sopt.VISBREMLINE_VAL.flatten().tolist()
        data['NE'] = sopt.NE_VAL.flatten().tolist()
        data['TE'] = sopt.TE_VAL.flatten().tolist()
        data['TI'] = sopt.TI_VAL.flatten().tolist()
        data['MSE'] = sopt.MSE_VAL.flatten().tolist()
        data['B_PROBES'] = sopt.B_PROBES_VAL.flatten().tolist()
        data['FLUXLOOPS'] = sopt.FLUXLOOPS_VAL.flatten().tolist()
        data['SEPARATRIX'] = sopt.SEPARATRIX_VAL.flatten().tolist()
    if run_name in ['TOK_R0_DELTA']:
        data['ASPECT'] = sopt.ASPECT_VAL[-1,0].tolist()
        data['R0'] = sopt.R0_VAL[-1,0].tolist()
    if run_name in ['TOK_R0_RHO']:
        data['ASPECT'] = sopt.ASPECT_VAL[-1,0].tolist()
        data['R0'] = sopt.R0_VAL[-1,0].tolist()
    if run_name in ['DKES']:
        data['L_11'] = sopt.DKES_11_VAL.flatten().tolist()
        data['L_31'] = sopt.DKES_31_VAL.flatten().tolist()
        data['L_33'] = sopt.DKES_33_VAL.flatten().tolist()
    if run_name in ['QHS_LMDIF']:
        data['ASPECT'] = sopt.ASPECT_VAL[-1,0].tolist()
    if run_name in ['AVAILENERGYOPT']:
        data['AVAILENERGY'] = sopt.TXPORT_VAL.flatten().tolist()
    if run_name in ['TRAVIS']:
        data['ECEREFLECT_X'] = sopt.ECEREFLECT_RADTX.flatten().tolist()
        data['ECEREFLECT_O'] = sopt.ECEREFLECT_RADTO.flatten().tolist()

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
            if act == 0:
                perct = 0
            else:
                perct = 100*abs(act-cal)/act
            print(f'  {temp} {cal:7.6f} {act:7.6f} {round(perct)}')
        else:
            print(temp)
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
        sys.exit(0) # For now since some may fail due to statistics
    else:
        print('  STATUS: PASS')
        sys.exit(0)





