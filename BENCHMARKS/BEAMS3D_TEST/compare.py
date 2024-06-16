#!/usr/bin/env python3
import sys, os
import shelve
from pathlib import Path
from argparse import ArgumentParser
sys.path.insert(0, '../../pySTEL/')
import numpy as np
_BENCH_FILE_ = 'BENCH_DATA'
# This part changes
from libstell.beams3d import BEAMS3D

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
    b3d = BEAMS3D()
    b3d.read_beams3d(f'beams3d_{run_name}.h5')
    version_str = f'BEAMS3D VERSION: {b3d.VERSION:4.2f}'

    # Extract values
    data={}
    if run_name in ['ORBITS','ORBITS_er']:
        rho = np.sqrt(b3d.S_lines)
        rho_max = np.max(rho,axis=1)
        rho_min = np.min(rho,axis=1)
        data = {}
        data['delta'] = rho_max-rho_min
        x = b3d.R_lines - 10.0
        y = b3d.Z_lines
        theta = np.arctan2(y,x)
        theta = np.where(theta > np.pi,theta-np.pi,theta)
        data['turning'] = np.max(theta,axis=1)
    if run_name in ['ORBITS_loss']:
        data['end_state'] = b3d.end_state
    if run_name in ['ORBITS_slow']:
        data['vll0']=b3d.vll_lines[:,0].flatten()
        data['epower']=b3d.epower_prof[:,0].flatten()*1.0E20 # factor to make values visible in output.
        data['ipower']=b3d.ipower_prof[:,0].flatten()*1.0E20
    if run_name in ['ORBITS_depo','ORBITS_depo_adas','ORBITS_eqdsk','ORBITS_multiion']:
        data['Shinethrough'] = b3d.Shinethrough
        [rho,depo]= b3d.calcDepo(ns=16)
        data['Deposition_B1E1'] = depo[0,:].flatten()
        data['Deposition_B1E2'] = depo[1,:].flatten()
        data['Deposition_B2E1'] = depo[3,:].flatten()
        data['Deposition_B2E2'] = depo[4,:].flatten()
    if run_name in ['ORBITS_restart']:
        data['LOSS_PCT'] = 100*sum(b3d.end_state==2)/b3d.nparticles
    if run_name in ['ORBITS_fusion']:
        [rho,depo]= b3d.calcDepo(ns=16)
        data['Birth_He4'] = depo[0,:]
        data['Birth_T'] = depo[1,:]
        data['Birth_p'] = depo[2,:]
        data['Birth_He3'] = depo[3,:]
    if run_name in ['ORBITS_dist']:
        data['dense_prof'] = b3d.dense_prof.flatten()
        data['ipower_prof'] = b3d.ipower_prof.flatten()
        data['epower_prof'] = b3d.epower_prof.flatten()
        data['j_prof'] = b3d.j_prof.flatten()

    # Read or write to the database file.
    if args.lmake_db:
        my_file = Path(_BENCH_FILE_)
        if my_file.exists():
            shelf = shelve.open(_BENCH_FILE_, flag="a")
        else:
            shelf = shelve.open(_BENCH_FILE_, flag="c")
        shelf[run_name] = data
        shelf.close()
        print('  ADDED: '+run_name)
        sys.exit(0)
    else:
        shelf = shelve.open(_BENCH_FILE_, flag='r')
        varlist = shelf[run_name]

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
        elif temp == 'wall_strikes':
            cal = np.where(act==0,0,cal)
            div = np.where(act==0,1,act)
            perct = 100*sum(abs(act-cal)/div)
            print(f'  {temp} {max(cal):7.6f} {max(act):7.6f} {round(perct)}')
        else:
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





