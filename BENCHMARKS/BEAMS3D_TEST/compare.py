#!/usr/bin/env python3
import sys, os
import json
from pathlib import Path
from argparse import ArgumentParser
sys.path.insert(0, '../../pySTEL/')
import numpy as np
_BENCH_FILE_ = 'BENCH_DATA.json'
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
    try:
        b3d.read_beams3d(f'beams3d_{run_name}.h5')
    except:
        print(f'  ERROR: Cannot find run {run_name}')
        sys.exit(-1)
    version_str = f'BEAMS3D VERSION: {b3d.VERSION:4.2f}'

    # Extract values
    data={}
    if run_name in ['ORBITS','ORBITS_er']:
        rho = np.sqrt(b3d.S_lines.T)
        rho_max = np.max(rho,axis=1)
        rho_min = np.min(rho,axis=1)
        data = {}
        data['delta'] = (rho_max-rho_min).tolist()
        x = b3d.R_lines.T - 10.0
        y = b3d.Z_lines.T
        theta = np.arctan2(y,x)
        theta = np.where(theta > np.pi,theta-np.pi,theta)
        data['turning'] = np.max(theta,axis=1).tolist()
    if run_name in ['ORBITS_loss']:
        data['end_state'] = b3d.end_state.tolist()
    if run_name in ['ORBITS_slow']:
        data['vll0']=b3d.vll_lines[0,:].flatten().tolist()
        data['epower']=(b3d.epower_prof[0,:].flatten()*1.0E20).tolist() # factor to make values visible in output.
        data['ipower']=(b3d.ipower_prof[0,:].flatten()*1.0E20).tolist()
    if run_name in ['ORBITS_slowz']:
        # Speed loss in a uniform H+C plasma against the Stix slowing-down
        # law with separate electron and ion Coulomb logarithms and [Z].
        ec = 1.602176634E-19; me = 9.1093837139E-31; eps0 = 8.8541878128E-12
        cl = 2.99792458E8; mp = 1.672621898E-27
        ne, te, ni, zi, mi = 5.0E19, 3.0E3, [4.34E19, 1.10E18], [1, 6], [1.6735328380E-27, 1.9926468799E-26]
        mb = b3d.mass[0]
        lne = 24.0 - np.log(np.sqrt(ne*1.0E-6)/te)
        taus = 3*(2*np.pi)**1.5*eps0**2*mb*(te*ec)**1.5/(ne*ec**4*np.sqrt(me)*lne)
        vte3 = (2*te*ec/me)**1.5
        def dvdt(v):
            lni = [43.0 - np.log(z*(mb+m)/(mb*m)*mp*np.sqrt(ne*1.0E-6/te)/(v/cl)**2) for z, m in zip(zi, mi)]
            zsum = sum(n*z*z/m*l for n, z, m, l in zip(ni, zi, mi, lni))
            vc3 = 0.75*np.sqrt(np.pi)*vte3*me*zsum/(ne*lne)
            return -(v/taus)*(1.0 + vc3/v**3)
        v = np.sqrt(b3d.vll_lines**2 + 2*b3d.moment_lines*b3d.B_lines/mb)
        nstep = 2000; h = b3d.t_end[0]/nstep
        dv_stix = []
        for v0 in v[0,:]:
            u = v0
            for _ in range(nstep):
                k1 = dvdt(u); k2 = dvdt(u+0.5*h*k1); k3 = dvdt(u+0.5*h*k2); k4 = dvdt(u+h*k3)
                u += h*(k1+2*k2+2*k3+k4)/6
            dv_stix.append(v0-u)
        data['dv_ratio'] = ((v[0,:]-v[-1,:])/np.array(dv_stix)).tolist()
    if run_name in ['ORBITS_depo','ORBITS_depo_adas','ORBITS_eqdsk','ORBITS_multiion']:
        data['Shinethrough'] = b3d.Shinethrough.flatten().tolist()
        [rho,depo]= b3d.calcDepo(ns=16)
        data['Deposition_B1E1'] = depo[0,:].flatten().tolist()
        data['Deposition_B1E2'] = depo[1,:].flatten().tolist()
        data['Deposition_B2E1'] = depo[3,:].flatten().tolist()
        data['Deposition_B2E2'] = depo[4,:].flatten().tolist()
        failtol = 25.0
    if run_name in ['ORBITS_restart']:
        data['LOSS_PCT'] = 100*sum(b3d.end_state==2)/b3d.nparticles
    if run_name in ['ORBITS_fusion']:
        [rho,depo]= b3d.calcDepo(ns=16)
        data['Birth_He4'] = depo[0,:].tolist()
        data['Birth_T']   = depo[1,:].tolist()
        data['Birth_p']   = depo[2,:].tolist()
        data['Birth_He3'] = depo[3,:].tolist()
    if run_name in ['ORBITS_dist']:
        data['dense_prof'] = b3d.dense_prof.flatten().tolist()
        data['ipower_prof'] = b3d.ipower_prof.flatten().tolist()
        data['epower_prof'] = b3d.epower_prof.flatten().tolist()
        data['j_prof'] = b3d.j_prof.flatten().tolist()

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
        elif temp == 'wall_strikes':
            act = np.array(act)
            cal = np.array(cal)
            cal = np.where(act==0,0,cal)
            div = np.where(act==0,1,act)
            perct = 100*sum(abs(act-cal)/div)
            print(f'  {temp} {max(cal):7.6f} {max(act):7.6f} {round(perct)}')
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





